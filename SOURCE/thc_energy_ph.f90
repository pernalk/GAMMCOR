
module thc_energy_ph

      ! =============================================================================
      ! ph-AC0 correlation energy with THC integrals. The ph analogue of
      ! thc_energy (THC_energy_loop), called at the end of ACPH0_THC.
      !
      ! THEORY (see .skills/phtheory_erpa_and_w1.md Sec. 7 and
      ! .skills/phtheory-E-contrib.md): W^(1) splits into nine contributions, one
      ! per pair of ph pair spaces (vo, va, ao, aa) minus aa-aa. Contribution
      ! <C>-<R> sums over the integral (pq|rs) with (pq) in C and (rs) in R and
      ! needs the single A1 block phA1<R><C> - name and block id read in opposite
      ! directions by construction.
      !
      !   contribution  printed    integral class   A1 block
      !   E_aoao        (21)(21)   (IA|IA)          phA1aoao
      !   E_aaao        (22)(21)   (AA|IA)          phA1aoaa
      !   E_voao        (31)(21)   (IV|IA)          phA1aovo
      !   E_voaa        (31)(22)   (IV|AA)          phA1aavo
      !   E_vovo        (31)(31)   (IV|IV)          none - closed form
      !   E_vaao        (32)(21)   (VA|IA)          phA1aova
      !   E_vaaa        (32)(22)   (VA|AA)          phA1aava
      !   E_vavo        (32)(31)   (VA|IV)          none - closed form
      !   E_vava        (32)(32)   (VA|VA)          none - closed form
      !
      ! CLOSED FORMS. phA1vovo, phA1vova and phA1vava are never built. These
      ! three contributions have Xt = Yt on both sides, so they need only
      ! A_plus^(1) - A_min^(1) = 2 B^(1), and B^(1) is a single integral term
      ! (.skills/phtheory-E-contrib.md, "Three of the nine collapse to
      ! integrals"), evaluated here on the fly from the same integral class:
      !   vo-vo  B_rs,pq = -2 (pq|rs) + (qr|ps)
      !   va-vo  B_rs,pq = n_q [ -2 (pq|rs) + (qr|ps) ]
      !   va-va  B_rs,pq = - sum_{tu in act} G_tuqs (ur|pt)
      !
      ! The printed labels are those of the (Mu)(Nu) table of AC0CAS_FOFO
      ! (ab0fofo.f90), 1 = inactive, 2 = active, 3 = virtual.
      !
      ! NOTATION. o means INACTIVE here, as in phac0_thc and the ph* notes. In
      ! thc_energy the same letter means occupied = inactive + active (RVO, ROO).
      !
      ! INTERMEDIATES (.skills/phcode_plan_energy_terms.md). All metric (c)
      ! factors live in three scaled eigenvector sets, Xp = Xt/d+, Ym = Yt/d-,
      ! Yb = d+ Yt (TvecXp, TvecYm, TvecYb, built by build_vecXt + build_vecs
      ! before the integral loop), so the A1 blocks are used raw. The numerator
      ! of each A1-based contribution, Xp^T A+raw Xp - Ym^T A-raw Ym, is
      ! projected, divided by the denominator and back-transformed into one
      ! loop-ready pair-space matrix G per contribution (TprojXAp, TprojYAm for
      ! the two terms with a vo side, TprojXApX, TprojYAmY for the four without,
      ! then TdS and TG; build_G_<term>), contracted with the integrals by
      ! G_contr.
      !
      ! All nine contributions are computed. Closed forms: E_vovo (vovo_contr),
      ! E_vavo (vavo_contr), E_vava (vava_contr). Via G and G_contr: E_voao,
      ! E_voaa, E_vaao, E_vaaa, E_aoao, E_aaao. For SA-CASSCF the A1 blocks
      ! aova, aoaa and aava arrive here already symmetrized with their
      ! transposed partners (symmetrize_A1_ph in phac0_thc,
      ! .skills/phcode_plan_energy_terms.md Sec. 8).
      ! =============================================================================

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals
      use types
      use lin
      use clock

      use math_constants
      use ppac_types
      use print_utils
      use thc_exact_ints, only: exact_V_axby, exact_V_xay
      use omp_lib

      implicit none

      integer, parameter, private :: tdebug = 0
      integer, parameter, private :: mdebugg = 50
      double precision, parameter, private :: ttoeV = 27.211399d+0

      ! ph block ids. Must match phac0_thc, which cannot be used here (it uses
      ! this module); pp duplicates its ids in thc_energy the same way.
      integer, parameter, private :: phA0vo=1, phA0aa=2, phA0va=3, phA0oa=4
      integer, parameter, private :: phA1aavo=5, phA1aova=6, phA1aovo=7
      integer, parameter, private :: phA1aoao=8, phA1aava=9, phA1aoaa=10
      ! transposed partners; folded into aova / aoaa / aava and freed by
      ! symmetrize_A1_ph before the energy loop - never read here
      integer, parameter, private :: phA1vaao=11, phA1aaao=12, phA1vaaa=13

      ! Xt eigenvectors of the two A0 blocks where Xt /= Yt (ao and aa); for
      ! vo and va Xt = Yt. Rows = pairs, columns = modes, in the layout of the
      ! Yt they are built from (see build_vecXt). Invalid modes are zero columns.
      !   ao(NA, NA, NI0:NI)  per correlated inactive s: rows r - NI, columns mode slot
      !   aa(D_aa, D_aa)   rows pos(r,s), columns mu
      type TvecXt
            double precision, dimension(:,:,:), allocatable :: ao
            double precision, dimension(:,:), allocatable :: aa
      end type TvecXt

      ! The three scaled eigenvector sets that carry all the metric factors of
      ! the energy, so that the A1 blocks are used raw (see build_vecs):
      !   Xp = Xt / d+      Ym = Yt / d-      Yb = d+ Yt       (per pair row)
      ! For a pair (1,2):  d+ = c_1 + c_2,   d- = c_1 - c_2,   with
      !   c_i = 1 (inactive),  c_v = 0 (virtual),  c_a = sign(n_a - 1/2) sqrt(n_a) (active)
      ! Same layout as TvecXt, plus va:
      !   va(NA, NA, NV)   per virtual p: rows q - NI, columns mode slot
      ! vo sides carry no vectors (Xp = 1, Ym = -1, Yb = 1, the 1/sqrt(2) go
      ! into the prefactor).
      type TvecXp
            double precision, dimension(:,:,:), allocatable :: ao, va
            double precision, dimension(:,:), allocatable :: aa
      end type TvecXp

      type TvecYm
            double precision, dimension(:,:,:), allocatable :: ao, va
            double precision, dimension(:,:), allocatable :: aa
      end type TvecYm

      type TvecYb
            double precision, dimension(:,:,:), allocatable :: ao, va
            double precision, dimension(:,:), allocatable :: aa
      end type TvecYb

      ! One-sided projections, for the two contributions with a vo column side
      ! (the vo side has no vectors). Rows = modes of the row space R, columns =
      ! vo pairs:
      !   projXAp = Xp_R^T A+raw          projYAm = Ym_R^T A-raw Ym_vo = - Ym_R^T A-raw
      type TprojXAp
            double precision, dimension(:,:), allocatable :: voao, voaa
      end type TprojXAp

      type TprojYAm
            double precision, dimension(:,:), allocatable :: voao, voaa
      end type TprojYAm

      ! Two-sided projections, for the four A1-based contributions without a vo
      ! side. Rows = modes of R, columns = modes of C:
      !   projXApX = Xp_R^T A+raw Xp_C          projYAmY = Ym_R^T A-raw Ym_C
      type TprojXApX
            double precision, dimension(:,:), allocatable :: vaao, vaaa, aoao, aaao
      end type TprojXApX

      type TprojYAmY
            double precision, dimension(:,:), allocatable :: vaao, vaaa, aoao, aaao
      end type TprojYAmY

      ! dS = (projX - projY) / (omega_k + omega_n), mode space, shape of the block
      type TdS
            double precision, dimension(:,:), allocatable :: voao, voaa, vaao, vaaa, aoao, aaao
      end type TdS

      ! G = f Yb_R dS Yb_C^T, pair space, loop-ready (prefactor f included);
      ! contracted with the integrals by G_contr
      type TG
            double precision, dimension(:,:), allocatable :: voao, voaa, vaao, vaaa, aoao, aaao
      end type TG

contains

      ! =============================================================================
      ! V_axby = (a x|b y) for one fixed (a, b), x and y local to their space.
      ! What it is in terms of the (pq|rs) of each contribution:
      !
      !   (IA|IA)  E_aoao   a=q (I)  x=p (A)   b=s (I)  y=r (A)
      !   (AA|IA)  E_aaao   a,x = p,q (A)      b=s (I)  y=r (A)
      !   (IV|IA)  E_voao   a=q (I)  x=p (V)   b=s (I)  y=r (A)
      !   (IV|AA)  E_voaa   a=q (I)  x=p (V)   b,y = r,s (A)
      !   (IV|IV)  E_vovo   a=q (I)  x=p (V)   b=s (I)  y=r (V)
      !   (VA|IA)  E_vaao   a=p (V)  x=q (A)   b=s (I)  y=r (A)
      !   (VA|AA)  E_vaaa   a=p (V)  x=q (A)   b,y = r,s (A)
      !   (VA|IV)  E_vavo   a=p (V)  x=q (A)   b=s (I)  y=r (V)   - all a at once, V_xay
      !   (VA|VA)  E_vava   a=p (V)  x=q (A)   b=r (V)  y=s (A)
      !
      ! The (a, b) loops run over the full ranges, i.e. over both orientations
      ! of the two pairs where both are in the same space - except (IV|IV) and
      ! (VA|VA), which run b <= a (see vovo_contr, vava_contr), and (VA|IV),
      ! which builds all virtual a in one gemm per b into V_xay(x,a,y) (see
      ! vavo_contr). Inactive a, b run NI0..NI = NFrozen+1..NI (frozen core);
      ! RIA and RIV are built for those orbitals only, with bounds NI0:NI.
      ! =============================================================================
      subroutine THC_energy_loop_ph(THCData, ACB, Flags, pos, AuxData, rdm2, E_corr_out, DeltaH_out)
            use THC_Gammcor

            type(TTHCData), intent(in) :: THCData
            type(TAC0Block_ph), dimension(:), intent(inout) :: ACB
            type(FlagsData), intent(in) :: Flags
            integer, dimension(:,:), intent(in) :: pos
            type(TACppData), intent(in) :: AuxData
            ! spin-summed active 2-RDM, the rdm2 of ACPH0_THC (E_vava)
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            ! E_corr of the reference, and Delta H_JI for J = trdm_st(:,1:NJ)
            ! (size 0 unless MSAC0)
            double precision, intent(out) :: E_corr_out
            double precision, dimension(:), allocatable, intent(out) :: DeltaH_out

            double precision, dimension(:), allocatable :: V_axby
            double precision, dimension(:,:,:), allocatable :: V_xay
            double precision, dimension(:,:,:), allocatable :: RII, RIA, RAA
            double precision, dimension(:,:,:), allocatable :: RIV, RVA
            integer, dimension(6) :: TabChol
            type(TLoop), dimension(6) :: loops
            type(TvecXt) :: vecXt
            type(TvecXp) :: vecXp
            type(TvecYm) :: vecYm
            type(TvecYb) :: vecYb
            ! MS-AC0: Y^(J<-I) per partner J, layout of vecYb, %aa unallocated
            type(TvecYb), dimension(:), allocatable :: vecYJ
            type(TG) :: G
            ! MS-AC0: G with one Yb replaced by YJ, per partner J (voaa stays empty)
            type(TG), dimension(:), allocatable :: GJ
            ! G_contr index vectors: column / row pair of G for slab index x / y
            integer, dimension(:), allocatable :: colidx, rowidx
            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            double precision, parameter :: AcAlpha = One
            integer :: NI, NA, NV, NIA, NBasis, NI0
            integer :: a, b, i, ii, x, y, j
            type (tclock) :: timer, timer0

            ! The nine contributions to W^(1), E_<space of pq><space of rs>.
            ! There is no E_aaaa. Each is to hold the COMPLETE contribution, the
            ! factor 2 of the six off-diagonal ones included, so that the printed
            ! numbers compare one to one with AC0CAS_FOFO and E_corr is their sum.
            double precision :: E_vovo, E_vavo, E_voao, E_voaa, E_vaao
            double precision :: E_vava, E_vaaa, E_aoao, E_aaao
            double precision :: E_corr
            ! MS-AC0: the same contributions with one Yb replaced by YJ, per
            ! partner J (E_vovo and E_voaa have none); DeltaH = Delta H_JI
            double precision, dimension(:), allocatable :: EJ_voao, EJ_vaao, EJ_vaaa
            double precision, dimension(:), allocatable :: EJ_aoao, EJ_aaao, EJ_vavo, EJ_vava
            double precision, dimension(:), allocatable :: DeltaH
            integer :: NJ

            ! Accumulated wall times per integral class, as in THC_int_loop_ph:
            ! t_<class>_int = building V_axby, t_<class>_ker = the kernels using it.
            double precision :: t_IAIA_int, t_IAIA_ker, t_AAIA_int, t_AAIA_ker
            double precision :: t_IVIA_int, t_IVIA_ker, t_IVAA_int, t_IVAA_ker
            double precision :: t_IVIV_int, t_IVIV_ker, t_VAIA_int, t_VAIA_ker
            double precision :: t_VAAA_int, t_VAAA_ker, t_VAIV_int, t_VAIV_ker
            double precision :: t_VAVA_int, t_VAVA_ker

            ! no associate for Zgk/Xga: they are unallocated in the exact-integral
            ! run, where no R block is built
            associate(NChol=>THCData%NChol)

              call print_section('phAC0 energy evaluation')

              NA = AuxData%NA
              NI = AuxData%NI
              NV = AuxData%NV
              NIA = NI + NA
              ! frozen core: inactive members of ph pairs run NI0..NI. Nothing
              ! here builds the reference (that stays in THC_int_loop_ph), so
              ! every inactive loop below is correlation side.
              NI0 = AuxData%NFrozen + 1
              NBasis = AuxData%NBasis

              ! Step 0.1: Calculate Xt for the ao and aa A0 blocks
              ! For ao and aa, Xt ≠ Yt
              ! For vo and va, Xt = Yt
              call clock_start(timer)
              call build_vecXt(vecXt, ACB(phA0oa), ACB(phA0aa), NA, NI0, NI)
              call tmsg('TIME FOR Xt (ao, aa)', timer, tdebug)

              ! Steps 0.2 + 0.3: drop invalid modes, scaled vector sets
              ! Xp, Ym, Yb. Consumes vecXt.
              call clock_start(timer)
              call build_vecs(vecXp, vecYm, vecYb, vecXt, ACB(phA0oa), ACB(phA0va), &
                    ACB(phA0aa), AuxData)
              call tmsg('TIME FOR Xp, Ym, Yb', timer, tdebug)

              ! MS-AC0: transition amplitudes Y^(J<-I), one vector set per partner
              ! J. NJ = 0 unless MSAC0, so the J loops below are then empty.
              call clock_start(timer)
              NJ = AuxData%NJ
              allocate(vecYJ(NJ), GJ(NJ))
              do j = 1, NJ
                    call build_YJ(vecYJ(j), vecXp, vecYm, ACB(phA0oa), ACB(phA0va), &
                          AuxData%trdm(:,:,j), NA, NV, NI0, NI)
              end do
              call tmsg('TIME FOR YJ', timer, tdebug)

              ! Step 2: G%voao from phA1aovo (consumes the A1 block)
              call clock_start(timer)
              call build_G_voao(G, ACB(phA1aovo), ACB(phA0oa), ACB(phA0vo), &
                    vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
              call tmsg('TIME FOR G voao', timer, tdebug)

              ! Step 3: G%voaa from phA1aavo (consumes the A1 block)
              call clock_start(timer)
              call build_G_voaa(G, ACB(phA1aavo), ACB(phA0aa), ACB(phA0vo), &
                    vecXp, vecYm, vecYb)
              call tmsg('TIME FOR G voaa', timer, tdebug)

              ! Step 4: G%vaao from phA1aova (consumes the A1 block)
              call clock_start(timer)
              call build_G_vaao(G, ACB(phA1aova), ACB(phA0oa), ACB(phA0va), &
                    vecXp, vecYm, vecYb, NA, NV, NI0, NI, vecYJ, GJ)
              call tmsg('TIME FOR G vaao', timer, tdebug)

              ! Step 5: G%vaaa from phA1aava (consumes the A1 block)
              call clock_start(timer)
              call build_G_vaaa(G, ACB(phA1aava), ACB(phA0aa), ACB(phA0va), &
                    vecXp, vecYm, vecYb, NA, NV, vecYJ, GJ)
              call tmsg('TIME FOR G vaaa', timer, tdebug)

              ! Step 6: G%aoao from phA1aoao (consumes the A1 block)
              call clock_start(timer)
              call build_G_aoao(G, ACB(phA1aoao), ACB(phA0oa), &
                    vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
              call tmsg('TIME FOR G aoao', timer, tdebug)

              ! Step 7: G%aaao from phA1aoaa (consumes the A1 block)
              call clock_start(timer)
              call build_G_aaao(G, ACB(phA1aoaa), ACB(phA0oa), ACB(phA0aa), &
                    vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
              call tmsg('TIME FOR G aaao', timer, tdebug)

              allocate(colidx(max(NA, NV)), rowidx(max(NA, NV)))

              E_vovo = zero
              E_vavo = zero
              E_voao = zero
              E_voaa = zero
              E_vaao = zero
              E_vava = zero
              E_vaaa = zero
              E_aoao = zero
              E_aaao = zero

              allocate(EJ_voao(NJ), EJ_vaao(NJ), EJ_vaaa(NJ), EJ_aoao(NJ), EJ_aaao(NJ), &
                    EJ_vavo(NJ), EJ_vava(NJ), DeltaH(NJ))
              EJ_voao = zero
              EJ_vaao = zero
              EJ_vaaa = zero
              EJ_aoao = zero
              EJ_aaao = zero
              EJ_vavo = zero
              EJ_vava = zero

              t_IAIA_int = zero
              t_IAIA_ker = zero
              t_AAIA_int = zero
              t_AAIA_ker = zero
              t_IVIA_int = zero
              t_IVIA_ker = zero
              t_IVAA_int = zero
              t_IVAA_ker = zero
              t_IVIV_int = zero
              t_IVIV_ker = zero
              t_VAIA_int = zero
              t_VAIA_ker = zero
              t_VAAA_int = zero
              t_VAAA_ker = zero
              t_VAIV_int = zero
              t_VAIV_ker = zero
              t_VAVA_int = zero
              t_VAVA_ker = zero

              ! One integral class per contribution (table in the module header).
              do i = 1, 6
                    loops(i)%run_main = .false.
                    loops(i)%run_inner = .false.
              end do

              loops(loopIA)%run_main = .true.
              loops(loopIA)%run_inner(2) = .true.    ! (IA|IA)  E_aoao
              loops(loopAA)%run_main = .true.
              loops(loopAA)%run_inner(2) = .true.    ! (AA|IA)  E_aaao
              loops(loopIV)%run_main = .true.
              loops(loopIV)%run_inner(2) = .true.    ! (IV|IA)  E_voao
              loops(loopIV)%run_inner(3) = .true.    ! (IV|AA)  E_voaa
              loops(loopIV)%run_inner(4) = .true.    ! (IV|IV)  E_vovo
              loops(loopVA)%run_main = .true.
              loops(loopVA)%run_inner(2) = .true.    ! (VA|IA)  E_vaao
              loops(loopVA)%run_inner(3) = .true.    ! (VA|AA)  E_vaaa
              loops(loopVA)%run_inner(4) = .true.    ! (VA|IV)  E_vavo
              loops(loopVA)%run_inner(5) = .true.    ! (VA|VA)  E_vava

              ! No class above has two inactive indices in one pair, so RII is
              ! not needed. Slot 6 (RVV) is never used in ph.
              TabChol = (/ 0, 1, 1, 1, 1, 0 /)
              ! exact-integral run: every slab comes from THCData%TwoEl
              if (THCData%exact_ints) TabChol = 0

              allocate(V_axby(max(NA, NI, NV)**2))

              ! (VA|IV) for one inactive b and ALL virtual a at once,
              ! V_xay(x,a,y) = (a x|b y): E_vavo pairs the slabs (p,s) and (r,s)
              ! (see vavo_contr). NA*NV**2 doubles.
              if (loops(loopVA)%run_inner(4)) then
                    allocate(V_xay(NA, NV, NV))
                    call CalcMem3(V_xay, 'V_xay')
              end if

              if (TabChol(1) == 1) then
                    allocate(RII(NChol, NI, NI))
                    call CalcMem3(RII, 'RII')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RII, THCData%Xga(:,1:NI), THCData%Xga(:,1:NI), &
                          THCData%Zgk, NI, NI, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RII', timer, tdebug)
              end if
              if (TabChol(2) == 1) then
                    allocate(RIA(NChol, NA, NI0:NI))
                    call CalcMem3(RIA, 'RIA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,NI0:NI), &
                          THCData%Zgk, NA, NI-NI0+1, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RIA', timer, tdebug)
              end if
              if (TabChol(3) == 1) then
                    allocate(RAA(NChol, NA, NA))
                    call CalcMem3(RAA, 'RAA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RAA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,NI+1:NIA), &
                          THCData%Zgk, NA, NA, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RAA', timer, tdebug)
              end if
              if (TabChol(4) == 1) then
                    allocate(RIV(NChol, NV, NI0:NI))
                    call CalcMem3(RIV, 'RIV')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIV, THCData%Xga(:,NIA+1:NBasis), THCData%Xga(:,NI0:NI), &
                          THCData%Zgk, NV, NI-NI0+1, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RIV', timer, tdebug)
              end if
              if (TabChol(5) == 1) then
                    allocate(RVA(NChol, NA, NV))
                    call CalcMem3(RVA, 'RVA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RVA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,NIA+1:NBasis), &
                          THCData%Zgk, NA, NV, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RVA', timer, tdebug)
              end if

              call clock_start(timer0)

              caseloop: do ii = 1, 6
                    if (loops(ii)%run_main)then
                          select case(ii)
                          case(loopIA)
                                do a = NI0, NI
                                      if (loops(ii)%run_inner(2)) then ! (IA|IA)
                                            do b = NI0, NI
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  end if
                                                  t_IAIA_int = t_IAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_aoao: (pq|rs) = V_axby(x,y), p = NI+x, q = a, r = NI+y, s = b.
                                                  ! Full a and b ranges: G%aoao is symmetric only for optimized orbitals.
                                                  do x = 1, NA
                                                        colidx(x) = pos(NI+x, a)     ! ao pair (p,q)
                                                  end do
                                                  do y = 1, NA
                                                        rowidx(y) = pos(NI+y, b)     ! ao pair (r,s)
                                                  end do
                                                  call G_contr(E_aoao, G%aoao, colidx, rowidx, NA, NA, V_axby)
                                                  do j = 1, NJ
                                                        call G_contr(EJ_aoao(j), GJ(j)%aoao, colidx, rowidx, NA, NA, V_axby)
                                                  end do
                                                  t_IAIA_ker = t_IAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          case(loopAA)
                                do a = NI+1, NIA
                                      if (loops(ii)%run_inner(2)) then ! (AA|IA)
                                            do b = NI0, NI
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  end if
                                                  t_AAIA_int = t_AAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_aaao: (pq|rs) = V_axby(x,y), (p,q) = (a, NI+x), r = NI+y, s = b
                                                  ! aa pair (p,q) once: p = a > q = NI+x; pos = 0 if ThrSelAct dropped it
                                                  do x = 1, NA
                                                        if (NI+x < a) then
                                                              colidx(x) = pos(a, NI+x)
                                                        else
                                                              colidx(x) = 0
                                                        end if
                                                  end do
                                                  do y = 1, NA
                                                        rowidx(y) = pos(NI+y, b)     ! ao pair (r,s)
                                                  end do
                                                  call G_contr(E_aaao, G%aaao, colidx, rowidx, NA, NA, V_axby)
                                                  do j = 1, NJ
                                                        call G_contr(EJ_aaao(j), GJ(j)%aaao, colidx, rowidx, NA, NA, V_axby)
                                                  end do
                                                  t_AAIA_ker = t_AAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          case(loopIV)
                                do a = NI0, NI
                                      if (loops(ii)%run_inner(2)) then ! (IV|IA)
                                            do b = NI0, NI
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                  end if
                                                  t_IVIA_int = t_IVIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_voao: (pq|rs) = V_axby(x,y), p = NIA+x, q = a, r = NI+y, s = b
                                                  do x = 1, NV
                                                        colidx(x) = pos(NIA+x, a)    ! vo pair (p,q)
                                                  end do
                                                  do y = 1, NA
                                                        rowidx(y) = pos(NI+y, b)     ! ao pair (r,s)
                                                  end do
                                                  call G_contr(E_voao, G%voao, colidx, rowidx, NV, NA, V_axby)
                                                  do j = 1, NJ
                                                        call G_contr(EJ_voao(j), GJ(j)%voao, colidx, rowidx, NV, NA, V_axby)
                                                  end do
                                                  t_IVIA_ker = t_IVIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(3)) then ! (IV|AA)
                                            do b = NI+1, NIA
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                  end if
                                                  t_IVAA_int = t_IVAA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_voaa: (pq|rs) = V_axby(x,y), p = NIA+x, q = a, (r,s) = (b, NI+y)
                                                  do x = 1, NV
                                                        colidx(x) = pos(NIA+x, a)    ! vo pair (p,q)
                                                  end do
                                                  ! aa pair (r,s) once: r = b > s = NI+y; pos = 0 if ThrSelAct dropped it
                                                  do y = 1, NA
                                                        if (NI+y < b) then
                                                              rowidx(y) = pos(b, NI+y)
                                                        else
                                                              rowidx(y) = 0
                                                        end if
                                                  end do
                                                  call G_contr(E_voaa, G%voaa, colidx, rowidx, NV, NA, V_axby)
                                                  t_IVAA_ker = t_IVAA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(4)) then ! (IV|IV)
                                            ! b <= a: slab (b,a) gives the same E_vovo as (a,b)
                                            do b = NI0, a
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, NIA+1, NBasis, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                  end if
                                                  t_IVIV_int = t_IVIV_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  call vovo_contr(E_vovo, ACB(phA0vo), pos, a, b, NIA+1, NBasis, V_axby)
                                                  t_IVIV_ker = t_IVIV_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          case(loopVA)
                                ! (VA|IV): all virtual a in one gemm per inactive b, V_xay(x,a,y) = (a x|b y)
                                if (loops(ii)%run_inner(4) .and. NA*NV > 0) then
                                      do b = NI0, NI
                                            call clock_start(timer)
                                            if (THCData%exact_ints) then
                                                  call exact_V_xay(V_xay, THCData%TwoEl, NI+1, NIA, NIA+1, NBasis, NIA+1, NBasis, b, AcAlpha)
                                            else
                                                  call real_aTb_x(V_xay, NA*NV, RVA, NChol, RIV(:,:,b), NChol, NA*NV, NV, NChol, AcAlpha)
                                            end if
                                            t_VAIV_int = t_VAIV_int + clock_readwall(timer)
                                            call clock_start(timer)
                                            call vavo_contr(E_vavo, ACB(phA0va), ACB(phA0vo), vecYb%va, pos, b, AuxData, &
                                                  NA, NV, V_xay, vecYJ, EJ_vavo)
                                            t_VAIV_ker = t_VAIV_ker + clock_readwall(timer)
                                      end do
                                end if
                                do a = NIA+1, NBasis
                                      if (loops(ii)%run_inner(2)) then ! (VA|IA)
                                            do b = NI0, NI
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  end if
                                                  t_VAIA_int = t_VAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_vaao: (pq|rs) = V_axby(x,y), p = a, q = NI+x, r = NI+y, s = b
                                                  do x = 1, NA
                                                        colidx(x) = pos(a, NI+x)     ! va pair (p,q)
                                                  end do
                                                  do y = 1, NA
                                                        rowidx(y) = pos(NI+y, b)     ! ao pair (r,s)
                                                  end do
                                                  call G_contr(E_vaao, G%vaao, colidx, rowidx, NA, NA, V_axby)
                                                  do j = 1, NJ
                                                        call G_contr(EJ_vaao(j), GJ(j)%vaao, colidx, rowidx, NA, NA, V_axby)
                                                  end do
                                                  t_VAIA_ker = t_VAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(3)) then ! (VA|AA)
                                            do b = NI+1, NIA
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                  end if
                                                  t_VAAA_int = t_VAAA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_vaaa: (pq|rs) = V_axby(x,y), p = a, q = NI+x, (r,s) = (b, NI+y)
                                                  do x = 1, NA
                                                        colidx(x) = pos(a, NI+x)     ! va pair (p,q)
                                                  end do
                                                  ! aa pair (r,s) once: r = b > s = NI+y; pos = 0 if ThrSelAct dropped it
                                                  do y = 1, NA
                                                        if (NI+y < b) then
                                                              rowidx(y) = pos(b, NI+y)
                                                        else
                                                              rowidx(y) = 0
                                                        end if
                                                  end do
                                                  call G_contr(E_vaaa, G%vaaa, colidx, rowidx, NA, NA, V_axby)
                                                  do j = 1, NJ
                                                        call G_contr(EJ_vaaa(j), GJ(j)%vaaa, colidx, rowidx, NA, NA, V_axby)
                                                  end do
                                                  t_VAAA_ker = t_VAAA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(5)) then ! (VA|VA)
                                            ! b <= a: slab (b,a) gives the same E_vava as (a,b)
                                            do b = NIA+1, a
                                                  call clock_start(timer)
                                                  if (THCData%exact_ints) then
                                                        call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                  else
                                                        call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RVA(:,:,b-NIA), NChol, NA, NA, NChol, AcAlpha)
                                                  end if
                                                  t_VAVA_int = t_VAVA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  call vava_contr(E_vava, ACB(phA0va), vecXp%va, vecYb%va, rdm2, &
                                                        AuxData%map, a, b, NIA, NA, V_axby, vecYJ, EJ_vava)
                                                  t_VAVA_ker = t_VAVA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          end select
                    end if
              end do caseloop

              call tmsg('ph energy loop', timer0, tdebug)
              call xmsg('t_IAIA_int', t_IAIA_int, mdebugg)
              call xmsg('t_IAIA_ker', t_IAIA_ker, mdebugg)
              call xmsg('t_AAIA_int', t_AAIA_int, mdebugg)
              call xmsg('t_AAIA_ker', t_AAIA_ker, mdebugg)
              call xmsg('t_IVIA_int', t_IVIA_int, mdebugg)
              call xmsg('t_IVIA_ker', t_IVIA_ker, mdebugg)
              call xmsg('t_IVAA_int', t_IVAA_int, mdebugg)
              call xmsg('t_IVAA_ker', t_IVAA_ker, mdebugg)
              call xmsg('t_IVIV_int', t_IVIV_int, mdebugg)
              call xmsg('t_IVIV_ker', t_IVIV_ker, mdebugg)
              call xmsg('t_VAIA_int', t_VAIA_int, mdebugg)
              call xmsg('t_VAIA_ker', t_VAIA_ker, mdebugg)
              call xmsg('t_VAAA_int', t_VAAA_int, mdebugg)
              call xmsg('t_VAAA_ker', t_VAAA_ker, mdebugg)
              call xmsg('t_VAIV_int', t_VAIV_int, mdebugg)
              call xmsg('t_VAIV_ker', t_VAIV_ker, mdebugg)
              call xmsg('t_VAVA_int', t_VAVA_int, mdebugg)
              call xmsg('t_VAVA_ker', t_VAVA_ker, mdebugg)

              E_corr = E_aoao + E_aaao + E_voao + E_voaa + E_vovo &
                    + E_vaao + E_vaaa + E_vavo + E_vava

              ! Same order and labels as the (Mu)(Nu) table of AC0CAS_FOFO.
              call print_section('Contributions to phAC0 from (Mu)(Nu) pairs of blocks')
              write(*,'(2X,A10,2X,A18)') 'Block pair', 'Contribution [Eh]'
              write(*,'(2X,A)') repeat('-', 30)
              call print_block_contribution('(21)(21)', E_aoao)
              call print_block_contribution('(22)(21)', E_aaao)
              call print_block_contribution('(31)(21)', E_voao)
              call print_block_contribution('(31)(22)', E_voaa)
              call print_block_contribution('(31)(31)', E_vovo)
              call print_block_contribution('(32)(21)', E_vaao)
              call print_block_contribution('(32)(22)', E_vaaa)
              call print_block_contribution('(32)(31)', E_vavo)
              call print_block_contribution('(32)(32)', E_vava)
              if (Flags%JobType == JOB_TYPE_MSAC0) then
                    write(*, '(1x, a, 2(1x, I0, ".", I0), F20.10)') 'H_mn', &
                          AuxData%nst(1), AuxData%nst(2), AuxData%nst(1), AuxData%nst(2), &
                          AuxData%ECAS_THC - AuxData%ENuc + E_corr
              end if

              call print_section('Final energies')
              call print_energy('CASSCF energy (one-electron)', AuxData%ECAS_oneelectr, &
                    ttoeV * AuxData%ECAS_oneelectr)
              call print_energy('ECASSCF_THC', AuxData%ECAS_THC, ttoeV * AuxData%ECAS_THC)
              call print_energy('E_phAC0', E_corr, ttoeV * E_corr)
              call print_energy('E_total', AuxData%ECAS_THC + E_corr, &
                    ttoeV * (AuxData%ECAS_THC + E_corr))
              print*, ''

              ! MS-AC0: one row of H^eff, Delta H_JI per partner J, I = reference
              do j = 1, NJ
                    DeltaH(j) = EJ_aoao(j) + EJ_aaao(j) + EJ_voao(j) &
                          + EJ_vaao(j) + EJ_vaaa(j) + EJ_vavo(j) + EJ_vava(j)
                    write(*, '(/, 1x, a, I0, ".", I0, a, I0, ".", I0)') &
                          'MS-AC0: contributions to Delta H_JI, J = ', &
                          AuxData%trdm_st(1, j), AuxData%trdm_st(2, j), &
                          ', I = ', AuxData%nst(1), AuxData%nst(2)
                    write(*,'(2X,A10,2X,A18)') 'Block pair', 'Contribution [Eh]'
                    write(*,'(2X,A)') repeat('-', 30)
                    call print_block_contribution('(21)(21)', EJ_aoao(j))
                    call print_block_contribution('(22)(21)', EJ_aaao(j))
                    call print_block_contribution('(31)(21)', EJ_voao(j))
                    call print_block_contribution('(31)(22)', zero)
                    call print_block_contribution('(31)(31)', zero)
                    call print_block_contribution('(32)(21)', EJ_vaao(j))
                    call print_block_contribution('(32)(22)', EJ_vaaa(j))
                    call print_block_contribution('(32)(31)', EJ_vavo(j))
                    call print_block_contribution('(32)(32)', EJ_vava(j))
                    write(*, '(1x, a, 2(1x, I0, ".", I0), F20.10)') 'H_mn', &
                          AuxData%trdm_st(1, j), AuxData%trdm_st(2, j), &
                          AuxData%nst(1), AuxData%nst(2), DeltaH(j)
              end do
              if (NJ > 0) print*, ''

              E_corr_out = E_corr
              call move_alloc(DeltaH, DeltaH_out)

            end associate

      end subroutine THC_energy_loop_ph


      ! =============================================================================
      ! Rebuild Xt for the two A0 blocks where Xt /= Yt
      ! (.skills/phtheory-E-contrib.md Sec. 2):
      !
      !   ao, per inactive s:  Xt(:,m)  = MiniAMinus(s) MiniAPlus(s)(:,m) / Eig(m)
      !   aa:                  Xt(:,mu) = AMinus APlus(:,mu) / Eig(mu)
      !
      ! After pherpa_symm_block_ph, APlus / MiniAPlus hold Yt and AMinus /
      ! MiniAMinus the symmetrized scaled A_-^(0). A mode with Eig <= 0 gets a
      ! zero column: that is an empty mini-block slot (Eig = 0, outside ActList)
      ! or an unstable mode (Eig = -sqrt|omega^2|, vector unnormalized). The
      ! mini-block matrices are zero outside ActList x ActList, so so is Xt%ao.
      ! =============================================================================
      subroutine build_vecXt(vecXt, OA, AA, NA, NI0, NI)
            type(TvecXt), intent(out) :: vecXt
            type(TAC0Block_ph), intent(in) :: OA, AA
            integer, intent(in) :: NA, NI0, NI
            integer :: s, m, nd

            allocate(vecXt%ao(NA, NA, NI0:NI))
            vecXt%ao = zero
            do s = NI0, NI
                  if (OA%MiniBlocks(s)%NDim == 0) cycle
                  associate(MB=>OA%MiniBlocks(s))
                    call real_ab(vecXt%ao(:,:,s), MB%MiniAMinus, MB%MiniAPlus)
                    do m = 1, NA
                          if (MB%Eig(m) > zero) then
                                vecXt%ao(:, m, s) = vecXt%ao(:, m, s) / MB%Eig(m)
                          else
                                vecXt%ao(:, m, s) = zero
                          end if
                    end do
                  end associate
            end do

            nd = AA%NDim
            allocate(vecXt%aa(nd, nd))
            if (nd == 0) return
            call real_ab(vecXt%aa, AA%AMinus, AA%APlus)
            do m = 1, nd
                  if (AA%Eig(m) > zero) then
                        vecXt%aa(:, m) = vecXt%aa(:, m) / AA%Eig(m)
                  else
                        vecXt%aa(:, m) = zero
                  end if
            end do

      end subroutine build_vecXt


      ! =============================================================================
      ! Steps 0.2 + 0.3 of .skills/phcode_plan_energy_terms.md. Builds the three
      ! scaled vector sets from Yt (APlus / MiniAPlus) and Xt (vecXt):
      !
      !   Xp = Xt / d+        Ym = Yt / d-        Yb = d+ Yt          (per pair row)
      !
      ! With A_+-^c = A_+-raw / (d+- d+-), Xt^T A_+^c Xt = Xp^T A_+raw Xp and
      ! Yt^T A_-^c Yt = Ym^T A_-raw Ym, and the d+_pq d+_rs of the energy goes
      ! into Yb. So the A1 blocks are used raw and no kernel needs a c factor.
      !
      !   space  pair    d+          d-          Xt
      !   ao     (r,s)   c_r + 1     c_r - 1     vecXt%ao
      !   va     (p,q)   c_q         -c_q        = Yt
      !   aa     (r,s)   c_r + c_s   c_r - c_s   vecXt%aa
      !
      ! 0.2, mode validity: a mode is used only if Eig > 0. That excludes the
      ! empty mini-block slots (Eig = 0, outside ActList) and the unstable modes
      ! pherpa_symm_block_ph flags with Eig = -sqrt|omega^2|. Their columns stay
      ! zero in all three sets, so they drop out of every term.
      !
      ! Division only on ActList rows (va, ao) and IndN rows (aa), where d is
      ! bounded away from 0 by ThrQVirt / ThrQInact / ThrSelAct; all other rows
      ! of Yt and Xt are zero. vecXt is consumed: its arrays move into vecXp.
      ! =============================================================================
      subroutine build_vecs(vecXp, vecYm, vecYb, vecXt, OA, VA, AA, AuxData)
            type(TvecXp), intent(out) :: vecXp
            type(TvecYm), intent(out) :: vecYm
            type(TvecYb), intent(out) :: vecYb
            type(TvecXt), intent(inout) :: vecXt
            type(TAC0Block_ph), intent(in) :: OA, VA, AA
            type(TACppData), intent(in) :: AuxData

            double precision, dimension(:), allocatable :: c
            double precision :: dp, dm
            integer :: k, m, jj, i, x, r, s, NI0

            associate(Occ=>AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, nd=>AA%NDim)

              NI0 = AuxData%NFrozen + 1

              allocate(c(NBasis))
              do i = 1, NBasis
                    c(i) = sqrt(Occ(i))
                    if (Occ(i) .lt. frac12) c(i) = -c(i)
              end do

              ! ---- ao: per correlated inactive k = NI0..NI, rows r - NI, r in ActList ----
              call move_alloc(vecXt%ao, vecXp%ao)
              allocate(vecYm%ao(NA, NA, NI0:NI), vecYb%ao(NA, NA, NI0:NI))
              vecYm%ao = zero
              vecYb%ao = zero
              do k = NI0, NI
                    associate(MB=>OA%MiniBlocks(k))
                      do m = 1, NA
                            ! invalid: Xt column already zero (build_vecXt)
                            if (MB%Eig(m) <= zero) cycle
                            do jj = 1, MB%NDim
                                  r = MB%ActList(jj)
                                  x = r - NI
                                  dp = c(r) + one
                                  dm = c(r) - one
                                  vecXp%ao(x, m, k) = vecXp%ao(x, m, k) / dp
                                  vecYm%ao(x, m, k) = MB%MiniAPlus(x, m) / dm
                                  vecYb%ao(x, m, k) = MB%MiniAPlus(x, m) * dp
                            end do
                      end do
                    end associate
              end do

              ! ---- va: per virtual k, rows q - NI with q in ActList, Xt = Yt ----
              allocate(vecXp%va(NA, NA, NV), vecYm%va(NA, NA, NV), vecYb%va(NA, NA, NV))
              vecXp%va = zero
              vecYm%va = zero
              vecYb%va = zero
              do k = 1, NV
                    associate(MB=>VA%MiniBlocks(k))
                      do m = 1, NA
                            if (MB%Eig(m) <= zero) cycle
                            do jj = 1, MB%NDim
                                  r = MB%ActList(jj)
                                  x = r - NI
                                  dp = c(r)
                                  dm = -c(r)
                                  vecXp%va(x, m, k) = MB%MiniAPlus(x, m) / dp
                                  vecYm%va(x, m, k) = MB%MiniAPlus(x, m) / dm
                                  vecYb%va(x, m, k) = MB%MiniAPlus(x, m) * dp
                            end do
                      end do
                    end associate
              end do

              ! ---- aa: rows i = pos(r,s), r > s, all IndN pairs ----
              call move_alloc(vecXt%aa, vecXp%aa)
              allocate(vecYm%aa(nd, nd), vecYb%aa(nd, nd))
              vecYm%aa = zero
              vecYb%aa = zero
              do m = 1, nd
                    if (AA%Eig(m) <= zero) cycle
                    do i = 1, nd
                          r = AA%IndN(1, i)
                          s = AA%IndN(2, i)
                          dp = c(r) + c(s)
                          dm = c(r) - c(s)
                          vecXp%aa(i, m) = vecXp%aa(i, m) / dp
                          vecYm%aa(i, m) = AA%APlus(i, m) / dm
                          vecYb%aa(i, m) = AA%APlus(i, m) * dp
                    end do
              end do

            end associate

      end subroutine build_vecs


      ! =============================================================================
      ! MS-AC0: transition amplitude Y^(J<-I) of the ao and va modes of the
      ! reference I, from gamma = gamma^JI(p,q) = 1/2 <J|E_qp|I> in the active NOs
      ! of I (.skills/phtheory_MS-AC0.md eqs. 21, 22). With the untilded
      ! X = (Xp - Ym)/2, Y = (Xp + Ym)/2:
      !
      !   ao, per correlated inactive k:  YJ = gamma^T Y - gamma X
      !   va, per virtual k:              YJ = gamma^T X        (Y_va = 0: Ym%va = -Xp%va)
      !
      ! Layout of vecYb (rows r - NI, columns mode slot); YJ replaces vecYb on one
      ! side of the off-diagonal energy. No aa: aa modes are excluded for J /= I,
      ! as in INCORE and FOFO. Rows outside ActList are zeroed: gamma spreads
      ! weight onto them, but those pairs are not in the ERPA of I.
      ! =============================================================================
      subroutine build_YJ(vecYJ, vecXp, vecYm, OA, VA, gamma, NA, NV, NI0, NI)
            type(TvecYb), intent(out) :: vecYJ
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TAC0Block_ph), intent(in) :: OA, VA
            double precision, dimension(:,:), intent(in) :: gamma
            integer, intent(in) :: NA, NV, NI0, NI
            double precision, dimension(:,:), allocatable :: Xm, Ymu
            logical, dimension(:), allocatable :: keep
            integer :: k, jj, x

            allocate(Xm(NA, NA), Ymu(NA, NA), keep(NA))

            allocate(vecYJ%ao(NA, NA, NI0:NI))
            do k = NI0, NI
                  Xm = (vecXp%ao(:,:,k) - vecYm%ao(:,:,k)) / two
                  Ymu = (vecXp%ao(:,:,k) + vecYm%ao(:,:,k)) / two
                  call real_aTb(vecYJ%ao(:,:,k), gamma, Ymu)
                  call real_ab(Ymu, gamma, Xm)
                  vecYJ%ao(:,:,k) = vecYJ%ao(:,:,k) - Ymu

                  associate(MB=>OA%MiniBlocks(k))
                    keep = .false.
                    do jj = 1, MB%NDim
                          keep(MB%ActList(jj) - NI) = .true.
                    end do
                  end associate
                  do x = 1, NA
                        if (.not. keep(x)) vecYJ%ao(x, :, k) = zero
                  end do
            end do

            allocate(vecYJ%va(NA, NA, NV))
            do k = 1, NV
                  call real_aTb(vecYJ%va(:,:,k), gamma, vecXp%va(:,:,k))

                  associate(MB=>VA%MiniBlocks(k))
                    keep = .false.
                    do jj = 1, MB%NDim
                          keep(MB%ActList(jj) - NI) = .true.
                    end do
                  end associate
                  do x = 1, NA
                        if (.not. keep(x)) vecYJ%va(x, :, k) = zero
                  end do
            end do

            deallocate(Xm, Ymu, keep)

      end subroutine build_YJ


      ! =============================================================================
      ! Step 2: G%voao for E_voao, contribution vo-ao, (pq) in vo = column side,
      ! (rs) in ao = row side, block phA1aovo (rows ao, columns vo). f = 2: both
      ! orientations, the two 1/sqrt(2) of the vo vectors and the 2 of FOFO
      ! step 3 (.skills/phcode_plan_energy_terms.md Sec. 1).
      !
      !   projXAp(k,pq) =   sum_rs Xp_ao(rs,k) A+raw(rs,pq)
      !   projYAm(k,pq) = - sum_rs Ym_ao(rs,k) A-raw(rs,pq)         (Ym_vo = -1)
      !   dS(k,pq)      = [projXAp - projYAm](k,pq) / (omega_k + omega_pq)
      !   G(rs,pq)      = 2 sum_k Yb_ao(rs,k) dS(k,pq)
      !
      ! The vo side has no vectors, so only rows are transformed, one NA-row
      ! slice per correlated inactive s: rows l0 = (s-NI0)*NA+1 .. l0+NA-1, in
      ! which pair (r,s) sits at r - NI and mode (s,nu) at its slot. Invalid
      ! modes (Eig <= 0) have zero vector columns and get dS = 0.
      !
      ! Storage chain, each stage in place in the previous one, peak = the A1
      ! block plus one NA x D_vo buffer:
      !   A1%APlus  -> projXAp%voao -> dS%voao -> G%voao
      !   A1%AMinus -> projYAm%voao, freed after the combine
      ! The A1 block is empty on return.
      !
      ! MS-AC0, per partner J: only the orientation with (pq) in ao survives
      ! (YJ_vo = 0), f = 1:   GJ(rs,pq) = sum_k YJ_ao(rs,k) dS(k,pq)
      ! =============================================================================
      subroutine build_G_voao(G, A1, OA, VO, vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: OA, VO
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            integer, intent(in) :: NA, NI0, NI
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            type(TG), dimension(:), intent(inout) :: GJ
            type(TprojXAp) :: projXAp
            type(TprojYAm) :: projYAm
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: buf
            double precision, dimension(:), allocatable :: om
            integer :: s, l0, l1, nrow, ncol, k, pq, jp

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXAp%voao)
            call move_alloc(A1%AMinus, projYAm%voao)

            allocate(buf(NA, ncol), om(nrow))

            ! rows -> ao modes, per correlated inactive s
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  om(l0:l1) = OA%MiniBlocks(s)%Eig(:)
                  call real_aTb_x(buf, NA, vecXp%ao(1,1,s), NA, projXAp%voao(l0,1), nrow, NA, ncol, NA)
                  projXAp%voao(l0:l1, :) = buf
                  call real_aTb_x(buf, NA, vecYm%ao(1,1,s), NA, projYAm%voao(l0,1), nrow, NA, ncol, NA, -one)
                  projYAm%voao(l0:l1, :) = buf
            end do

            ! combine, in place in projXAp
            !$omp parallel do default(shared) private(pq, k)
            do pq = 1, ncol
                  do k = 1, nrow
                        if (om(k) > zero) then
                              projXAp%voao(k, pq) = (projXAp%voao(k, pq) - projYAm%voao(k, pq)) &
                                    / (om(k) + VO%Eig(pq))
                        else
                              projXAp%voao(k, pq) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAm%voao)
            call move_alloc(projXAp%voao, dS%voao)

            ! MS-AC0: GJ from dS, before the back-transform below overwrites it
            do jp = 1, size(GJ)
                  allocate(GJ(jp)%voao(nrow, ncol))
                  do s = NI0, NI
                        l0 = (s - NI0) * NA + 1
                        call real_ab_x(GJ(jp)%voao(l0,1), nrow, vecYJ(jp)%ao(1,1,s), NA, &
                              dS%voao(l0,1), nrow, NA, ncol, NA)
                  end do
            end do

            ! back-transform: ao modes -> ao pairs, with f = 2
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_ab_x(buf, NA, vecYb%ao(1,1,s), NA, dS%voao(l0,1), nrow, NA, ncol, NA, two)
                  dS%voao(l0:l1, :) = buf
            end do
            call move_alloc(dS%voao, G%voao)

      end subroutine build_G_voao


      ! =============================================================================
      ! Step 3: G%voaa for E_voaa, contribution vo-aa, (pq) in vo = column side,
      ! (rs) in aa = row side, block phA1aavo (rows aa, columns vo). f = 2, as
      ! for E_voao. The aa twin of build_G_voao:
      !
      !   projXAp(mu,pq) =   sum_rs Xp_aa(rs,mu) A+raw(rs,pq)
      !   projYAm(mu,pq) = - sum_rs Ym_aa(rs,mu) A-raw(rs,pq)       (Ym_vo = -1)
      !   dS(mu,pq)      = [projXAp - projYAm](mu,pq) / (omega_mu + omega_pq)
      !   G(rs,pq)       = 2 sum_mu Yb_aa(rs,mu) dS(mu,pq)
      !
      ! Rows are aa pairs pos(r,s), r > s, and aa modes mu, both D_aa = AA%NDim.
      ! The aa side has no mini-blocks, so each row transform is one dense
      ! D_aa x D_aa gemm. Columns are independent, so it is done in place over
      ! column chunks of width nchunk with a D_aa x nchunk buffer. Invalid modes
      ! (Eig <= 0) have zero vector columns and get dS = 0.
      !
      ! No symmetrization is needed for SA-CASSCF here: an aa pair and a vo pair
      ! cannot share an orbital, so phA1aavo is symmetric with its transpose
      ! for any orbitals (.skills/phcode_plan_energy_terms.md Sec. 8).
      !
      ! Storage chain, as in build_G_voao:
      !   A1%APlus  -> projXAp%voaa -> dS%voaa -> G%voaa
      !   A1%AMinus -> projYAm%voaa, freed after the combine
      ! =============================================================================
      subroutine build_G_voaa(G, A1, AA, VO, vecXp, vecYm, vecYb)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: AA, VO
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            type(TprojXAp) :: projXAp
            type(TprojYAm) :: projYAm
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: buf
            integer, parameter :: nchunk_max = 1024
            integer :: nrow, ncol, nchunk, c0, nc, mu, pq

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXAp%voaa)
            call move_alloc(A1%AMinus, projYAm%voaa)

            nchunk = max(1, min(ncol, nchunk_max))
            allocate(buf(nrow, nchunk))

            ! rows -> aa modes, per column chunk
            if (nrow > 0) then
                  do c0 = 1, ncol, nchunk
                        nc = min(nchunk, ncol - c0 + 1)
                        call real_aTb_x(buf, nrow, vecXp%aa, nrow, projXAp%voaa(1,c0), nrow, nrow, nc, nrow)
                        projXAp%voaa(:, c0:c0+nc-1) = buf(:, 1:nc)
                        call real_aTb_x(buf, nrow, vecYm%aa, nrow, projYAm%voaa(1,c0), nrow, nrow, nc, nrow, -one)
                        projYAm%voaa(:, c0:c0+nc-1) = buf(:, 1:nc)
                  end do
            end if

            ! combine, in place in projXAp
            !$omp parallel do default(shared) private(pq, mu)
            do pq = 1, ncol
                  do mu = 1, nrow
                        if (AA%Eig(mu) > zero) then
                              projXAp%voaa(mu, pq) = (projXAp%voaa(mu, pq) - projYAm%voaa(mu, pq)) &
                                    / (AA%Eig(mu) + VO%Eig(pq))
                        else
                              projXAp%voaa(mu, pq) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAm%voaa)
            call move_alloc(projXAp%voaa, dS%voaa)

            ! back-transform: aa modes -> aa pairs, with f = 2
            if (nrow > 0) then
                  do c0 = 1, ncol, nchunk
                        nc = min(nchunk, ncol - c0 + 1)
                        call real_ab_x(buf, nrow, vecYb%aa, nrow, dS%voaa(1,c0), nrow, nrow, nc, nrow, two)
                        dS%voaa(:, c0:c0+nc-1) = buf(:, 1:nc)
                  end do
            end if
            call move_alloc(dS%voaa, G%voaa)

      end subroutine build_G_voaa


      ! =============================================================================
      ! Step 4: G%vaao for E_vaao, contribution va-ao, (pq) in va = column side,
      ! (rs) in ao = row side, block phA1aova (rows ao, columns va). f = 4: both
      ! orientations and the 2 of FOFO step 3, no vo side.
      !
      !   projXApX(k,n) = sum_{rs,pq} Xp_ao(rs,k) A+raw(rs,pq) Xp_va(pq,n)
      !   projYAmY(k,n) = sum_{rs,pq} Ym_ao(rs,k) A-raw(rs,pq) Ym_va(pq,n)
      !   dS(k,n)       = [projXApX - projYAmY](k,n) / (omega_(s,nu) + omega_(p,mu))
      !   G(rs,pq)      = 4 sum_{k,n} Yb_ao(rs,k) dS(k,n) Yb_va(pq,n)
      !
      ! Both sides are mini-blocked, so both transforms go slice by slice, as
      ! pp's AZvaao:
      !   rows    per correlated inactive s: rows l0 = (s-NI0)*NA+1 .. l0+NA-1,
      !           pair (r,s) at r - NI, mode (s,nu) at its slot;
      !   columns per virtual p (k = p - NIA): columns (k-1)*NA+1 .. k*NA,
      !           pair (p,q) at q - NI, mode (p,mu) at its slot.
      ! The minus-channel sign of the va side is inside Ym_va = -Yt/c. A mode
      ! pair with Eig <= 0 on either side gets dS = 0.
      !
      ! SA-CASSCF: phA1aova is NOT symmetric with its transpose unless the
      ! state-specific orbital gradient vanishes (shared active orbital, r = q,
      ! dB term). symmetrize_A1_ph has already replaced it by its symmetric
      ! part, so the doubling of the two orientations here is exact.
      !
      ! Storage chain, peak = the A1 block plus an NA x D_va and a D_ao x NA
      ! buffer:
      !   A1%APlus  -> projXApX%vaao -> dS%vaao -> G%vaao
      !   A1%AMinus -> projYAmY%vaao, freed after the combine
      !
      ! MS-AC0, per partner J: both orientations survive, f = 2 each:
      !   GJ = 2 [ YJ_ao dS Yb_va^T + Yb_ao dS YJ_va^T ]
      ! The first term needs dS itself (one extra D_ao x D_va buffer); the second
      ! reuses 4 Yb_ao dS from the row back-transform, scaled by 1/2.
      ! =============================================================================
      subroutine build_G_vaao(G, A1, OA, VA, vecXp, vecYm, vecYb, NA, NV, NI0, NI, vecYJ, GJ)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: OA, VA
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            integer, intent(in) :: NA, NV, NI0, NI
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            type(TG), dimension(:), intent(inout) :: GJ
            type(TprojXApX) :: projXApX
            type(TprojYAmY) :: projYAmY
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: bufr, bufc, tmp
            double precision, dimension(:), allocatable :: omr, omc
            integer :: s, k, l0, l1, c0, c1, nrow, ncol, i, j, jp

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXApX%vaao)
            call move_alloc(A1%AMinus, projYAmY%vaao)

            allocate(bufr(NA, ncol), bufc(nrow, NA), omr(nrow), omc(ncol))

            ! mode energies: rows ao (s,nu), columns va (p,mu)
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  omr(l0:l0+NA-1) = OA%MiniBlocks(s)%Eig(:)
            end do
            do k = 1, NV
                  c0 = (k - 1) * NA + 1
                  omc(c0:c0+NA-1) = VA%MiniBlocks(k)%Eig(:)
            end do

            ! rows -> ao modes, per correlated inactive s
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_aTb_x(bufr, NA, vecXp%ao(1,1,s), NA, projXApX%vaao(l0,1), nrow, NA, ncol, NA)
                  projXApX%vaao(l0:l1, :) = bufr
                  call real_aTb_x(bufr, NA, vecYm%ao(1,1,s), NA, projYAmY%vaao(l0,1), nrow, NA, ncol, NA)
                  projYAmY%vaao(l0:l1, :) = bufr
            end do

            ! columns -> va modes, per virtual k (nrow = 0 if all inactive are
            ! frozen; BLAS rejects lda = 0)
            if (nrow > 0) then
                  do k = 1, NV
                        c0 = (k - 1) * NA + 1
                        c1 = c0 + NA - 1
                        call real_ab_x(bufc, nrow, projXApX%vaao(1,c0), nrow, vecXp%va(1,1,k), NA, nrow, NA, NA)
                        projXApX%vaao(:, c0:c1) = bufc
                        call real_ab_x(bufc, nrow, projYAmY%vaao(1,c0), nrow, vecYm%va(1,1,k), NA, nrow, NA, NA)
                        projYAmY%vaao(:, c0:c1) = bufc
                  end do
            end if

            ! combine, in place in projXApX
            !$omp parallel do default(shared) private(i, j)
            do j = 1, ncol
                  do i = 1, nrow
                        if (omr(i) > zero .and. omc(j) > zero) then
                              projXApX%vaao(i, j) = (projXApX%vaao(i, j) - projYAmY%vaao(i, j)) &
                                    / (omr(i) + omc(j))
                        else
                              projXApX%vaao(i, j) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAmY%vaao)
            call move_alloc(projXApX%vaao, dS%vaao)

            ! MS-AC0: GJ = 2 YJ_ao dS Yb_va^T, from dS before it is overwritten
            do jp = 1, size(GJ)
                  allocate(GJ(jp)%vaao(nrow, ncol))
            end do
            if (size(GJ) > 0 .and. nrow > 0) then
                  allocate(tmp(nrow, ncol))
                  do jp = 1, size(GJ)
                        do s = NI0, NI
                              l0 = (s - NI0) * NA + 1
                              call real_ab_x(tmp(l0,1), nrow, vecYJ(jp)%ao(1,1,s), NA, &
                                    dS%vaao(l0,1), nrow, NA, ncol, NA, two)
                        end do
                        do k = 1, NV
                              c0 = (k - 1) * NA + 1
                              call real_abT_x(GJ(jp)%vaao(1,c0), nrow, tmp(1,c0), nrow, &
                                    vecYb%va(1,1,k), NA, nrow, NA, NA)
                        end do
                  end do
                  deallocate(tmp)
            end if

            ! back-transform rows: ao modes -> ao pairs, with f = 4
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_ab_x(bufr, NA, vecYb%ao(1,1,s), NA, dS%vaao(l0,1), nrow, NA, ncol, NA, four)
                  dS%vaao(l0:l1, :) = bufr
            end do

            ! MS-AC0: GJ += 2 Yb_ao dS YJ_va^T = 1/2 (4 Yb_ao dS) YJ_va^T
            if (nrow > 0) then
                  do jp = 1, size(GJ)
                        do k = 1, NV
                              c0 = (k - 1) * NA + 1
                              call real_abT_x(GJ(jp)%vaao(1,c0), nrow, dS%vaao(1,c0), nrow, &
                                    vecYJ(jp)%va(1,1,k), NA, nrow, NA, NA, frac12, one)
                        end do
                  end do
            end if

            ! back-transform columns: va modes -> va pairs, dS Yb_va^T
            if (nrow > 0) then
                  do k = 1, NV
                        c0 = (k - 1) * NA + 1
                        c1 = c0 + NA - 1
                        call real_abT_x(bufc, nrow, dS%vaao(1,c0), nrow, vecYb%va(1,1,k), NA, nrow, NA, NA)
                        dS%vaao(:, c0:c1) = bufc
                  end do
            end if
            call move_alloc(dS%vaao, G%vaao)

      end subroutine build_G_vaao


      ! =============================================================================
      ! Step 5: G%vaaa for E_vaaa, contribution va-aa, (pq) in va = column side,
      ! (rs) in aa = row side, block phA1aava (rows aa, columns va). f = 4, as
      ! for E_vaao.
      !
      !   projXApX(mu,n) = sum_{rs,pq} Xp_aa(rs,mu) A+raw(rs,pq) Xp_va(pq,n)
      !   projYAmY(mu,n) = sum_{rs,pq} Ym_aa(rs,mu) A-raw(rs,pq) Ym_va(pq,n)
      !   dS(mu,n)       = [projXApX - projYAmY](mu,n) / (omega_mu + omega_(p,mu'))
      !   G(rs,pq)       = 4 sum_{mu,n} Yb_aa(rs,mu) dS(mu,n) Yb_va(pq,n)
      !
      ! Rows are aa pairs pos(r,s), r > s, and aa modes, both D_aa; the aa side
      ! is dense. Columns are va pairs, per virtual p (k = p - NIA) the slice
      ! (k-1)*NA+1 .. k*NA, pair (p,q) at q - NI, mode (p,mu) at its slot. The
      ! two sides are transformed together, slice by slice:
      !   A(:,slice) <- Xp_aa^T A(:,slice) Xp_va(:,:,k)
      ! which keeps both buffers at D_aa x NA. The minus-channel sign of the va
      ! side is inside Ym_va = -Yt/c. A mode pair with Eig <= 0 on either side
      ! gets dS = 0.
      !
      ! SA-CASSCF: phA1aava is NOT symmetric with its transpose unless the
      ! state-specific orbital gradient vanishes (shared active orbital, s = q
      ! dA term and r = q dB term). symmetrize_A1_ph has already replaced it by
      ! its symmetric part, so the doubling of the two orientations here is exact.
      !
      ! Storage chain:
      !   A1%APlus  -> projXApX%vaaa -> dS%vaaa -> G%vaaa
      !   A1%AMinus -> projYAmY%vaaa, freed after the combine
      !
      ! MS-AC0, per partner J: only the orientation with (pq) in va survives
      ! (YJ_aa = 0), f = 2:   GJ = 2 Yb_aa dS YJ_va^T
      ! =============================================================================
      subroutine build_G_vaaa(G, A1, AA, VA, vecXp, vecYm, vecYb, NA, NV, vecYJ, GJ)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: AA, VA
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            integer, intent(in) :: NA, NV
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            type(TG), dimension(:), intent(inout) :: GJ
            type(TprojXApX) :: projXApX
            type(TprojYAmY) :: projYAmY
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: buf1, buf2
            double precision, dimension(:), allocatable :: omc
            integer :: k, c0, c1, nrow, ncol, i, j, jp

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXApX%vaaa)
            call move_alloc(A1%AMinus, projYAmY%vaaa)

            allocate(buf1(nrow, NA), buf2(nrow, NA), omc(ncol))

            ! mode energies of the columns, va (p,mu)
            do k = 1, NV
                  c0 = (k - 1) * NA + 1
                  omc(c0:c0+NA-1) = VA%MiniBlocks(k)%Eig(:)
            end do

            ! both sides -> modes, per virtual k (nrow = D_aa = 0 possible:
            ! BLAS rejects lda = 0)
            if (nrow > 0) then
                  do k = 1, NV
                        c0 = (k - 1) * NA + 1
                        c1 = c0 + NA - 1
                        call real_aTb_x(buf1, nrow, vecXp%aa, nrow, projXApX%vaaa(1,c0), nrow, nrow, NA, nrow)
                        call real_ab_x(buf2, nrow, buf1, nrow, vecXp%va(1,1,k), NA, nrow, NA, NA)
                        projXApX%vaaa(:, c0:c1) = buf2
                        call real_aTb_x(buf1, nrow, vecYm%aa, nrow, projYAmY%vaaa(1,c0), nrow, nrow, NA, nrow)
                        call real_ab_x(buf2, nrow, buf1, nrow, vecYm%va(1,1,k), NA, nrow, NA, NA)
                        projYAmY%vaaa(:, c0:c1) = buf2
                  end do
            end if

            ! combine, in place in projXApX
            !$omp parallel do default(shared) private(i, j)
            do j = 1, ncol
                  do i = 1, nrow
                        if (AA%Eig(i) > zero .and. omc(j) > zero) then
                              projXApX%vaaa(i, j) = (projXApX%vaaa(i, j) - projYAmY%vaaa(i, j)) &
                                    / (AA%Eig(i) + omc(j))
                        else
                              projXApX%vaaa(i, j) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAmY%vaaa)
            call move_alloc(projXApX%vaaa, dS%vaaa)

            ! back-transform both sides, per virtual k, with f = 4:
            !   G(:,slice) = 4 Yb_aa dS(:,slice) Yb_va(:,:,k)^T
            ! MS-AC0: GJ(:,slice) = 1/2 (4 Yb_aa dS(:,slice)) YJ_va(:,:,k)^T
            do jp = 1, size(GJ)
                  allocate(GJ(jp)%vaaa(nrow, ncol))
            end do
            if (nrow > 0) then
                  do k = 1, NV
                        c0 = (k - 1) * NA + 1
                        c1 = c0 + NA - 1
                        call real_ab_x(buf1, nrow, vecYb%aa, nrow, dS%vaaa(1,c0), nrow, nrow, NA, nrow, four)
                        do jp = 1, size(GJ)
                              call real_abT_x(GJ(jp)%vaaa(1,c0), nrow, buf1, nrow, &
                                    vecYJ(jp)%va(1,1,k), NA, nrow, NA, NA, frac12)
                        end do
                        call real_abT_x(buf2, nrow, buf1, nrow, vecYb%va(1,1,k), NA, nrow, NA, NA)
                        dS%vaaa(:, c0:c1) = buf2
                  end do
            end if
            call move_alloc(dS%vaaa, G%vaaa)

      end subroutine build_G_vaaa


      ! =============================================================================
      ! Step 6: G%aoao for E_aoao, contribution ao-ao, (pq) and (rs) both in ao,
      ! block phA1aoao (rows ao, columns ao). Diagonal combination: f = 2, the
      ! 2 of FOFO step 3 only.
      !
      !   projXApX(k,n) = sum_{rs,pq} Xp_ao(rs,k) A+raw(rs,pq) Xp_ao(pq,n)
      !   projYAmY(k,n) = sum_{rs,pq} Ym_ao(rs,k) A-raw(rs,pq) Ym_ao(pq,n)
      !   dS(k,n)       = [projXApX - projYAmY](k,n) / (omega_(s,nu) + omega_(q,mu))
      !   G(rs,pq)      = 2 sum_{k,n} Yb_ao(rs,k) dS(k,n) Yb_ao(pq,n)
      !
      ! Rows and columns have the same layout: per correlated inactive s the
      ! slice (s-NI0)*NA+1 .. +NA, pair (r,s) at r - NI, mode (s,nu) at its
      ! slot. Two passes, as in build_G_vaao: rows per s, then columns per q.
      ! A mode pair with Eig <= 0 on either side gets dS = 0.
      !
      ! SA-CASSCF: no correction needed. The energy sums over all (pq) and all
      ! (rs) of this one block, so it only sees the symmetric part of A^(1).
      !
      ! Storage chain:
      !   A1%APlus  -> projXApX%aoao -> dS%aoao -> G%aoao
      !   A1%AMinus -> projYAmY%aoao, freed after the combine
      !
      ! MS-AC0, per partner J: YJ on the (pq) = column side, f = 2:
      !   GJ = 2 Yb_ao dS YJ_ao^T
      ! =============================================================================
      subroutine build_G_aoao(G, A1, OA, vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: OA
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            integer, intent(in) :: NA, NI0, NI
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            type(TG), dimension(:), intent(inout) :: GJ
            type(TprojXApX) :: projXApX
            type(TprojYAmY) :: projYAmY
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: bufr, bufc
            double precision, dimension(:), allocatable :: om
            integer :: s, l0, l1, nrow, ncol, i, j, jp

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXApX%aoao)
            call move_alloc(A1%AMinus, projYAmY%aoao)

            allocate(bufr(NA, ncol), bufc(nrow, NA), om(nrow))

            ! mode energies, ao (s,nu); rows and columns share them
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  om(l0:l0+NA-1) = OA%MiniBlocks(s)%Eig(:)
            end do

            ! rows -> ao modes, per correlated inactive s
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_aTb_x(bufr, NA, vecXp%ao(1,1,s), NA, projXApX%aoao(l0,1), nrow, NA, ncol, NA)
                  projXApX%aoao(l0:l1, :) = bufr
                  call real_aTb_x(bufr, NA, vecYm%ao(1,1,s), NA, projYAmY%aoao(l0,1), nrow, NA, ncol, NA)
                  projYAmY%aoao(l0:l1, :) = bufr
            end do

            ! columns -> ao modes, per correlated inactive q (loop empty if
            ! nrow = 0)
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_ab_x(bufc, nrow, projXApX%aoao(1,l0), nrow, vecXp%ao(1,1,s), NA, nrow, NA, NA)
                  projXApX%aoao(:, l0:l1) = bufc
                  call real_ab_x(bufc, nrow, projYAmY%aoao(1,l0), nrow, vecYm%ao(1,1,s), NA, nrow, NA, NA)
                  projYAmY%aoao(:, l0:l1) = bufc
            end do

            ! combine, in place in projXApX
            !$omp parallel do default(shared) private(i, j)
            do j = 1, ncol
                  do i = 1, nrow
                        if (om(i) > zero .and. om(j) > zero) then
                              projXApX%aoao(i, j) = (projXApX%aoao(i, j) - projYAmY%aoao(i, j)) &
                                    / (om(i) + om(j))
                        else
                              projXApX%aoao(i, j) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAmY%aoao)
            call move_alloc(projXApX%aoao, dS%aoao)

            ! back-transform rows: ao modes -> ao pairs, with f = 2
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_ab_x(bufr, NA, vecYb%ao(1,1,s), NA, dS%aoao(l0,1), nrow, NA, ncol, NA, two)
                  dS%aoao(l0:l1, :) = bufr
            end do

            ! MS-AC0: GJ = (2 Yb_ao dS) YJ_ao^T, before the column pass below
            do jp = 1, size(GJ)
                  allocate(GJ(jp)%aoao(nrow, ncol))
                  do s = NI0, NI
                        l0 = (s - NI0) * NA + 1
                        call real_abT_x(GJ(jp)%aoao(1,l0), nrow, dS%aoao(1,l0), nrow, &
                              vecYJ(jp)%ao(1,1,s), NA, nrow, NA, NA)
                  end do
            end do

            ! back-transform columns: dS Yb_ao^T
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  l1 = l0 + NA - 1
                  call real_abT_x(bufc, nrow, dS%aoao(1,l0), nrow, vecYb%ao(1,1,s), NA, nrow, NA, NA)
                  dS%aoao(:, l0:l1) = bufc
            end do
            call move_alloc(dS%aoao, G%aoao)

      end subroutine build_G_aoao


      ! =============================================================================
      ! Step 7: G%aaao for E_aaao, contribution aa-ao, (pq) in aa = column side,
      ! (rs) in ao = row side, block phA1aoaa (rows ao, columns aa). f = 4: both
      ! orientations and the 2 of FOFO step 3, no vo side.
      !
      !   projXApX(k,mu) = sum_{rs,pq} Xp_ao(rs,k) A+raw(rs,pq) Xp_aa(pq,mu)
      !   projYAmY(k,mu) = sum_{rs,pq} Ym_ao(rs,k) A-raw(rs,pq) Ym_aa(pq,mu)
      !   dS(k,mu)       = [projXApX - projYAmY](k,mu) / (omega_(s,nu) + omega_mu)
      !   G(rs,pq)       = 4 sum_{k,mu} Yb_ao(rs,k) dS(k,mu) Yb_aa(pq,mu)
      !
      ! The mirror of build_G_vaaa: rows are mini-blocked, per correlated
      ! inactive s the slice (s-NI0)*NA+1 .. +NA, pair (r,s) at r - NI, mode
      ! (s,nu) at its slot; columns are aa pairs pos(p,q), p > q, and aa modes,
      ! both D_aa, dense. The two sides are transformed together, row slice by
      ! row slice:
      !   A(slice,:) <- Xp_ao(:,:,s)^T A(slice,:) Xp_aa
      ! with two NA x D_aa buffers. A mode pair with Eig <= 0 on either side
      ! gets dS = 0.
      !
      ! SA-CASSCF: phA1aoaa is NOT symmetric with its transpose unless the
      ! state-specific orbital gradient vanishes (shared active orbital, r = p
      ! dA term and r = q dB term). symmetrize_A1_ph has already replaced it by
      ! its symmetric part, so the doubling of the two orientations here is exact.
      !
      ! Storage chain:
      !   A1%APlus  -> projXApX%aaao -> dS%aaao -> G%aaao
      !   A1%AMinus -> projYAmY%aaao, freed after the combine
      !
      ! MS-AC0, per partner J: only the orientation with (pq) in ao survives
      ! (YJ_aa = 0), f = 2:   GJ = 2 YJ_ao dS Yb_aa^T
      ! =============================================================================
      subroutine build_G_aaao(G, A1, OA, AA, vecXp, vecYm, vecYb, NA, NI0, NI, vecYJ, GJ)
            type(TG), intent(inout) :: G
            type(TAC0Block_ph), intent(inout) :: A1
            type(TAC0Block_ph), intent(in) :: OA, AA
            type(TvecXp), intent(in) :: vecXp
            type(TvecYm), intent(in) :: vecYm
            type(TvecYb), intent(in) :: vecYb
            integer, intent(in) :: NA, NI0, NI
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            type(TG), dimension(:), intent(inout) :: GJ
            type(TprojXApX) :: projXApX
            type(TprojYAmY) :: projYAmY
            type(TdS) :: dS
            double precision, dimension(:,:), allocatable :: buf1, buf2
            double precision, dimension(:), allocatable :: omr
            integer :: s, l0, l1, nrow, ncol, i, j, jp

            nrow = A1%NDim1
            ncol = A1%NDim2

            call move_alloc(A1%APlus, projXApX%aaao)
            call move_alloc(A1%AMinus, projYAmY%aaao)

            allocate(buf1(NA, ncol), buf2(NA, ncol), omr(nrow))

            ! mode energies of the rows, ao (s,nu)
            do s = NI0, NI
                  l0 = (s - NI0) * NA + 1
                  omr(l0:l0+NA-1) = OA%MiniBlocks(s)%Eig(:)
            end do

            ! both sides -> modes, per correlated inactive s (ncol = D_aa = 0
            ! possible: BLAS rejects ldb = 0)
            if (ncol > 0) then
                  do s = NI0, NI
                        l0 = (s - NI0) * NA + 1
                        l1 = l0 + NA - 1
                        call real_aTb_x(buf1, NA, vecXp%ao(1,1,s), NA, projXApX%aaao(l0,1), nrow, NA, ncol, NA)
                        call real_ab_x(buf2, NA, buf1, NA, vecXp%aa, ncol, NA, ncol, ncol)
                        projXApX%aaao(l0:l1, :) = buf2
                        call real_aTb_x(buf1, NA, vecYm%ao(1,1,s), NA, projYAmY%aaao(l0,1), nrow, NA, ncol, NA)
                        call real_ab_x(buf2, NA, buf1, NA, vecYm%aa, ncol, NA, ncol, ncol)
                        projYAmY%aaao(l0:l1, :) = buf2
                  end do
            end if

            ! combine, in place in projXApX
            !$omp parallel do default(shared) private(i, j)
            do j = 1, ncol
                  do i = 1, nrow
                        if (omr(i) > zero .and. AA%Eig(j) > zero) then
                              projXApX%aaao(i, j) = (projXApX%aaao(i, j) - projYAmY%aaao(i, j)) &
                                    / (omr(i) + AA%Eig(j))
                        else
                              projXApX%aaao(i, j) = zero
                        end if
                  end do
            end do
            !$omp end parallel do
            deallocate(projYAmY%aaao)
            call move_alloc(projXApX%aaao, dS%aaao)

            ! MS-AC0: GJ(slice,:) = 2 YJ_ao(:,:,s) dS(slice,:) Yb_aa^T, from dS
            ! before the back-transform below overwrites it
            do jp = 1, size(GJ)
                  allocate(GJ(jp)%aaao(nrow, ncol))
                  if (ncol == 0) cycle
                  do s = NI0, NI
                        l0 = (s - NI0) * NA + 1
                        call real_ab_x(buf1, NA, vecYJ(jp)%ao(1,1,s), NA, dS%aaao(l0,1), nrow, NA, ncol, NA, two)
                        call real_abT_x(GJ(jp)%aaao(l0,1), nrow, buf1, NA, vecYb%aa, ncol, NA, ncol, ncol)
                  end do
            end do

            ! back-transform both sides, per correlated inactive s, with f = 4:
            !   G(slice,:) = 4 Yb_ao(:,:,s) dS(slice,:) Yb_aa^T
            if (ncol > 0) then
                  do s = NI0, NI
                        l0 = (s - NI0) * NA + 1
                        l1 = l0 + NA - 1
                        call real_ab_x(buf1, NA, vecYb%ao(1,1,s), NA, dS%aaao(l0,1), nrow, NA, ncol, NA, four)
                        call real_abT_x(buf2, NA, buf1, NA, vecYb%aa, ncol, NA, ncol, ncol)
                        dS%aaao(l0:l1, :) = buf2
                  end do
            end if
            call move_alloc(dS%aaao, G%aaao)

      end subroutine build_G_aaao


      ! =============================================================================
      ! Step 0.4: the one loop kernel of the six A1-based contributions. For one
      ! integral slab V_axby(x,y) = (a x|b y),
      !
      !   E += sum_{x,y} V_axby(x,y) * G(rowidx(y), colidx(x))
      !
      ! G is the loop-ready pair-space matrix of the contribution (prefactor and
      ! all c factors already inside), rows = row pairs (rs), columns = column
      ! pairs (pq). colidx(x) is the column pair slab index x belongs to,
      ! rowidx(y) the row pair; 0 skips the element (aa pair taken once, or
      ! removed by ThrSelAct).
      ! =============================================================================
      subroutine G_contr(E, G, colidx, rowidx, nx, ny, V_axby)
            double precision, intent(inout) :: E
            double precision, dimension(:,:), intent(in) :: G
            integer, intent(in) :: nx, ny
            integer, dimension(nx), intent(in) :: colidx
            integer, dimension(ny), intent(in) :: rowidx
            double precision, dimension(nx, ny), intent(in) :: V_axby
            double precision :: val
            integer :: x, y

            val = zero
            !$omp parallel do default(shared) private(x, y) reduction(+:val)
            do x = 1, nx
                  if (colidx(x) == 0) cycle
                  do y = 1, ny
                        if (rowidx(y) == 0) cycle
                        val = val + V_axby(x, y) * G(rowidx(y), colidx(x))
                  end do
            end do
            !$omp end parallel do

            E = E + val

      end subroutine G_contr


      ! =============================================================================
      ! E_vovo from one (IV|IV) slab, V_axby(x,y) = (x a|y b) = (pq|rs) with
      ! p = x, q = a, r = y, s = b. The ph twin of vvoo_contr (thc_energy): both
      ! A0 sides are diagonal, so no eigenvectors and no intermediate, and the A1
      ! element is made of integrals from the same slab. Closed form in the code
      ! convention (.skills/phtheory-E-contrib.md Sec. 3.1):
      !
      !   E_vovo = sum_{pq, rs}  (pq|rs) [ (qr|ps) - 2 (pq|rs) ] / (omega_pq + omega_rs)
      !
      ! with (qr|ps) = (y a|x b) = V_axby(y,x) and omega_pq = f_pp - f_qq =
      ! VO%Eig(pos(p,q)). This is the closed-shell MP2 expression for the
      ! inactive -> virtual doubles with the CAS orbital energies.
      !
      ! Slab (b,a) is the transpose of slab (a,b), and after x <-> y it gives
      ! exactly the same sum (no A1 symmetry is assumed - it is all integrals).
      ! So the caller runs b <= a and the a /= b slabs count twice.
      ! The sum itself is vovo_slab.
      ! =============================================================================
      subroutine vovo_contr(E_vovo, VO, pos, a, b, x0, x1, V_axby)
            double precision, intent(inout) :: E_vovo
            type(TAC0Block_ph), intent(in) :: VO
            integer, dimension(:,:), intent(in) :: pos
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1
            double precision, dimension(x0:x1, x0:x1), intent(in) :: V_axby
            double precision, dimension(x0:x1) :: om_a, om_b
            integer :: x

            ! vo has no filtered pairs, pos(x,a) > 0 always
            do x = x0, x1
                  om_a(x) = VO%Eig(pos(x, a))
                  om_b(x) = VO%Eig(pos(x, b))
            end do

            call vovo_slab(E_vovo, om_a, om_b, a, b, x0, x1, V_axby)

      end subroutine vovo_contr


      ! =============================================================================
      ! The sum of vovo_contr for one (IV|IV) slab, with the vo omegas given:
      !
      !   E_vovo += sum_{x,y} V_axby(x,y) [V_axby(y,x) - 2 V_axby(x,y)] / (om_a(x) + om_b(y))
      !
      ! om_a(x) = omega of the vo pair (x,a), om_b(y) = omega of (y,b). Doubled
      ! for a /= b (the caller runs b <= a). Also used by mp2_energy_thc (mp.f90),
      ! which takes omega = fvw - fij directly instead of the vo A0 block.
      ! =============================================================================
      subroutine vovo_slab(E_vovo, om_a, om_b, a, b, x0, x1, V_axby)
            double precision, intent(inout) :: E_vovo
            integer, intent(in) :: x0, x1
            double precision, dimension(x0:x1), intent(in) :: om_a, om_b
            integer, intent(in) :: a, b
            double precision, dimension(x0:x1, x0:x1), intent(in) :: V_axby
            double precision :: val
            integer :: x, y

            val = zero
            !$omp parallel do default(shared) private(x, y) reduction(+:val)
            do y = x0, x1
                  do x = x0, x1
                        val = val + V_axby(x,y) * (V_axby(y,x) - two * V_axby(x,y)) &
                              / (om_a(x) + om_b(y))
                  end do
            end do
            !$omp end parallel do

            if (a /= b) val = two * val
            E_vovo = E_vovo + val

      end subroutine vovo_slab


      ! =============================================================================
      ! E_vavo for one inactive s, from V_xay(x,a,y) = (a x|s y) over all virtual
      ! a, y. Closed form in the code convention (.skills/phtheory-E-contrib.md
      ! Sec. 3.2), p, r virtual, q, b active:
      !
      !   E_va-vo = 2 sum c_q (pq|rs) sum_mu sum_b c_b Yt_(p,mu)(q) Yt_(p,mu)(b)
      !               [ (br|ps) - 2 (pb|rs) ] / (omega_(p,mu) + omega_rs)
      !
      ! and E_vavo is twice that (va-vo plus vo-va). With the c-scaled vectors of
      ! the va mini-block of p, Yb(x,m) = c_x Yt_(p,m)(x) = vecYb%va(x,m,p):
      !
      !   U(m,r) = sum_x Yb(x,m) (p x|s r) = sum_x Yb(x,m) V_xay(x,p,r)  - (pq|rs), (pb|rs)
      !   T(m,r) = sum_x Yb(x,m) (r x|s p) = sum_x Yb(x,m) V_xay(x,r,p)  - (br|ps)
      !   E_vavo += 4 sum_{p,m,r} U(m,r) [ T(m,r) - 2 U(m,r) ] / (omega_(p,m) + omega_rs)
      !
      ! T takes slab a = r at column p, U slab a = p - the reason the caller builds
      ! all virtual a at once for one s. m runs over the packed mode slots
      ! ActList(jj) - NI of the mini-block; the other slots are empty, and the
      ! columns of invalid modes are zero in Yb (build_vecs).
      !
      ! MS-AC0, per partner J: only the va-vo orientation survives (YJ_vo = 0).
      ! YJ replaces Yb in U only; T - 2U is the B^(1) numerator and keeps Yb:
      !   Ubar(m,r) = sum_x YJ_va(x,m) V_xay(x,p,r)
      !   EJ(J)    += 2 sum_{p,m,r} Ubar(m,r) [ T(m,r) - 2 U(m,r) ] / (omega_(p,m) + omega_rs)
      ! =============================================================================
      subroutine vavo_contr(E_vavo, VA, VO, Yb, pos, s, AuxData, NA, NV, V_xay, vecYJ, EJ)
            double precision, intent(inout) :: E_vavo
            type(TAC0Block_ph), intent(in) :: VA, VO
            integer, intent(in) :: s, NA, NV
            double precision, dimension(NA, NA, NV), intent(in) :: Yb
            integer, dimension(:,:), intent(in) :: pos
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(NA, NV, NV), intent(in) :: V_xay
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            double precision, dimension(:), intent(inout) :: EJ
            double precision, dimension(:,:), allocatable :: U, T, Ub
            double precision, dimension(:), allocatable :: om_s, valJ
            double precision :: val
            integer :: p, r, m, jj, jp

            associate(NI=>AuxData%NI, NIA=>AuxData%NIA)

              allocate(om_s(NV), U(NA, NV), T(NA, NV), Ub(NA, NV), valJ(size(EJ)))

              ! vo has no filtered pairs, pos(r,s) > 0 always
              do r = 1, NV
                    om_s(r) = VO%Eig(pos(NIA+r, s))
              end do

              val = zero
              valJ = zero
              do p = 1, NV
                    if (VA%MiniBlocks(p)%NDim == 0) cycle
                    associate(MB=>VA%MiniBlocks(p))
                      call real_aTb_x(U, NA, Yb(1, 1, p), NA, V_xay(1, p, 1), NA*NV, NA, NV, NA)
                      call real_aTb_x(T, NA, Yb(1, 1, p), NA, V_xay(1, 1, p), NA, NA, NV, NA)
                      do r = 1, NV
                            do jj = 1, MB%NDim
                                  m = MB%ActList(jj) - NI
                                  if (MB%Eig(m) <= zero) cycle
                                  val = val + U(m, r) * (T(m, r) - two * U(m, r)) &
                                        / (MB%Eig(m) + om_s(r))
                            end do
                      end do
                      do jp = 1, size(EJ)
                            call real_aTb_x(Ub, NA, vecYJ(jp)%va(1, 1, p), NA, V_xay(1, p, 1), NA*NV, NA, NV, NA)
                            do r = 1, NV
                                  do jj = 1, MB%NDim
                                        m = MB%ActList(jj) - NI
                                        if (MB%Eig(m) <= zero) cycle
                                        valJ(jp) = valJ(jp) + Ub(m, r) * (T(m, r) - two * U(m, r)) &
                                              / (MB%Eig(m) + om_s(r))
                                  end do
                            end do
                      end do
                    end associate
              end do

              E_vavo = E_vavo + four * val
              EJ = EJ + two * valJ

            end associate

      end subroutine vavo_contr


      ! =============================================================================
      ! E_vava from one (VA|VA) slab, V_axby(x,y) = (a x|b y) = (pq|rs) with
      ! p = a, q = x, r = b, s = y. Closed form in the code convention
      ! (.skills/phtheory-E-contrib.md Sec. 3.6), p, r virtual, q, s, b, c active:
      !
      !   E_vava = -4 sum c_q c_s (pq|rs) sum_{mu,nu} Yt_(p,mu)(q) Yt_(r,nu)(s)
      !              / (omega_(p,mu) + omega_(r,nu))
      !              * sum_bc Yt_(r,nu)(b) Yt_(p,mu)(c) / (c_b c_c) * sum_tu G_tucb (ur|pt)
      !
      ! Diagonal combination, no factor 2. With Yb = c Yt and Xp = Yt / c of the
      ! va mini-blocks of p and r (build_vecs) it factorizes per slab:
      !
      !   M1(mu,nu) = sum_qs Yb^(p)(q,mu) V(q,s) Yb^(r)(s,nu)     = Yb^(p)T V Yb^(r)
      !   K(c,b)    = sum_tu G_tucb (ur|pt) = sum_tu rdm2(t,u,c,b) V(t,u)
      !   M2(mu,nu) = sum_cb Xp^(p)(c,mu) K(c,b) Xp^(r)(b,nu)     = Xp^(p)T K Xp^(r)
      !   E_vava   += -4 sum_{mu,nu valid} M1(mu,nu) M2(mu,nu) / (omega_(p,mu) + omega_(r,nu))
      !
      ! (ur|pt) = (pt|ru) = V(t,u). In the convention of D6_vava, G_tucb =
      ! rdm2(u,t,b,c) = rdm2(t,u,c,b) by the real-RDM symmetry
      ! rdm2(1,2,3,4) = rdm2(2,1,4,3), so K is one ddot over NA**2 per (c,b).
      !
      ! Slab (b,a) is the transpose of slab (a,b); with the same symmetry K and
      ! M1, M2 transpose with it and give the same sum. So the caller runs
      ! b <= a and the a /= b slabs count twice.
      !
      ! MS-AC0, per partner J: YJ replaces one Yb in M1 (M2 is the numerator and
      ! is unchanged). Slab (b,a) is no longer the same sum, so for a /= b both
      ! slabs are written out instead of doubling:
      !   a /= b:  M1J = YJ^(p)T V Yb^(r) + Yb^(p)T V YJ^(r)
      !   a  = b:  M1J = YJ^(p)T V Yb^(p)
      !   EJ(J)  += -4 sum M1J(mu,nu) M2(mu,nu) / (omega_(p,mu) + omega_(r,nu))
      ! =============================================================================
      subroutine vava_contr(E_vava, VA, Xp, Yb, rdm2, map, a, b, NIA, NA, V_axby, vecYJ, EJ)
            double precision, intent(inout) :: E_vava
            type(TAC0Block_ph), intent(in) :: VA
            integer, intent(in) :: a, b, NIA, NA
            double precision, dimension(:,:,:), intent(in) :: Xp, Yb
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, dimension(:), intent(in) :: map
            double precision, dimension(NA, NA), intent(in) :: V_axby
            type(TvecYb), dimension(:), intent(in) :: vecYJ
            double precision, dimension(:), intent(inout) :: EJ
            double precision, dimension(NA, NA) :: W, WY, M1, K, M2, M1J
            double precision :: val, valJ
            integer :: p, r, c, d, mu, nu, NI, jp

            p = a - NIA
            r = b - NIA
            NI = NIA - NA

            associate(MBp=>VA%MiniBlocks(p), MBr=>VA%MiniBlocks(r))

              if (MBp%NDim == 0 .or. MBr%NDim == 0) return

              ! M1 = Yb^(p)T V Yb^(r)
              call real_ab(WY, V_axby, Yb(:,:,r))
              call real_aTb(M1, Yb(:,:,p), WY)

              ! K(c,d) = sum_tu rdm2(t,u,c,d) V(t,u)
              !$omp parallel do default(shared) private(c, d) collapse(2)
              do d = 1, NA
                    do c = 1, NA
                          call real_vw_x(K(c, d), rdm2(:,:,map(NI+c),map(NI+d)), V_axby, NA**2)
                    end do
              end do
              !$omp end parallel do

              ! M2 = Xp^(p)T K Xp^(r)
              call real_ab(W, K, Xp(:,:,r))
              call real_aTb(M2, Xp(:,:,p), W)

              ! invalid modes: zero columns in Xp, Yb, but Eig <= 0 in the denominator
              val = zero
              do nu = 1, NA
                    if (MBr%Eig(nu) <= zero) cycle
                    do mu = 1, NA
                          if (MBp%Eig(mu) <= zero) cycle
                          val = val + M1(mu, nu) * M2(mu, nu) / (MBp%Eig(mu) + MBr%Eig(nu))
                    end do
              end do

              do jp = 1, size(EJ)
                    call real_aTb(M1J, vecYJ(jp)%va(:,:,p), WY)
                    if (a /= b) then
                          call real_ab(W, V_axby, vecYJ(jp)%va(:,:,r))
                          call real_aTb(M1, Yb(:,:,p), W)
                          M1J = M1J + M1
                    end if
                    valJ = zero
                    do nu = 1, NA
                          if (MBr%Eig(nu) <= zero) cycle
                          do mu = 1, NA
                                if (MBp%Eig(mu) <= zero) cycle
                                valJ = valJ + M1J(mu, nu) * M2(mu, nu) / (MBp%Eig(mu) + MBr%Eig(nu))
                          end do
                    end do
                    EJ(jp) = EJ(jp) - four * valJ
              end do

            end associate

            if (a /= b) val = two * val
            E_vava = E_vava - four * val

      end subroutine vava_contr

end module thc_energy_ph
