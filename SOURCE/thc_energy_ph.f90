
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
      ! of each A1-based contribution, Xp^T A+raw Xp - Ym^T A-raw Ym, is to be
      ! projected, divided by the denominator and back-transformed into one
      ! loop-ready pair-space matrix G per contribution, contracted with the
      ! integrals by G_contr (not written yet: TprojXAp, TprojYAm, TprojXApX,
      ! TprojYAmY, TdS, TG).
      !
      ! CURRENT SCOPE: Step 0 of the plan (Xt, Xp, Ym, Yb, G_contr), the R??
      ! vectors and the integrals of the nine classes are built; E_vovo
      ! (vovo_contr), E_vavo (vavo_contr) and E_vava (vava_contr) have kernels,
      ! the other six are still zero.
      ! =============================================================================

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals
      use types
      use lin
      use clock

      use math_constants
      use ppac_types
      use print_utils
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

      ! Xt eigenvectors of the two A0 blocks where Xt /= Yt (ao and aa); for
      ! vo and va Xt = Yt. Rows = pairs, columns = modes, in the layout of the
      ! Yt they are built from (see build_vecXt). Invalid modes are zero columns.
      !   ao(NA, NA, NI)   per inactive s: rows r - NI, columns mode slot
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
      ! vavo_contr).
      ! =============================================================================
      subroutine THC_energy_loop_ph(THCData, ACB, Flags, pos, AuxData, rdm2)
            use THC_Gammcor

            type(TTHCData), intent(in) :: THCData
            type(TAC0Block_ph), dimension(:), intent(inout) :: ACB
            type(FlagsData), intent(in) :: Flags
            integer, dimension(:,:), intent(in) :: pos
            type(TACppData), intent(in) :: AuxData
            ! spin-summed active 2-RDM, the rdm2 of ACPH0_THC (E_vava)
            double precision, dimension(:,:,:,:), intent(in) :: rdm2

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
            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            double precision, parameter :: AcAlpha = One
            integer :: NI, NA, NV, NIA, NBasis
            integer :: a, b, i, ii
            type (tclock) :: timer, timer0

            ! The nine contributions to W^(1), E_<space of pq><space of rs>.
            ! There is no E_aaaa. Each is to hold the COMPLETE contribution, the
            ! factor 2 of the six off-diagonal ones included, so that the printed
            ! numbers compare one to one with AC0CAS_FOFO and E_corr is their sum.
            double precision :: E_vovo, E_vavo, E_voao, E_voaa, E_vaao
            double precision :: E_vava, E_vaaa, E_aoao, E_aaao
            double precision :: E_corr

            ! Accumulated wall times per integral class, as in THC_int_loop_ph:
            ! t_<class>_int = building V_axby, t_<class>_ker = the kernels using it.
            double precision :: t_IAIA_int, t_IAIA_ker, t_AAIA_int, t_AAIA_ker
            double precision :: t_IVIA_int, t_IVIA_ker, t_IVAA_int, t_IVAA_ker
            double precision :: t_IVIV_int, t_IVIV_ker, t_VAIA_int, t_VAIA_ker
            double precision :: t_VAAA_int, t_VAAA_ker, t_VAIV_int, t_VAIV_ker
            double precision :: t_VAVA_int, t_VAVA_ker

            associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga, &
                  NChol=>THCData%NChol, NTHC=>THCData%NTHC)

              call print_section('phAC0 energy evaluation')

              NA = AuxData%NA
              NI = AuxData%NI
              NV = AuxData%NV
              NIA = NI + NA
              NBasis = AuxData%NBasis

              ! Step 0.1: Calculate Xt for the ao and aa A0 blocks
              ! For ao and aa, Xt ≠ Yt
              ! For vo and va, Xt = Yt
              call clock_start(timer)
              call build_vecXt(vecXt, ACB(phA0oa), ACB(phA0aa), NA, NI)
              call tmsg('TIME FOR Xt (ao, aa)', timer, tdebug)

              ! Steps 0.2 + 0.3: drop invalid modes, scaled vector sets
              ! Xp, Ym, Yb. Consumes vecXt.
              call clock_start(timer)
              call build_vecs(vecXp, vecYm, vecYb, vecXt, ACB(phA0oa), ACB(phA0va), &
                    ACB(phA0aa), AuxData)
              call tmsg('TIME FOR Xp, Ym, Yb', timer, tdebug)

              E_vovo = zero
              E_vavo = zero
              E_voao = zero
              E_voaa = zero
              E_vaao = zero
              E_vava = zero
              E_vaaa = zero
              E_aoao = zero
              E_aaao = zero

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
                    call thc_gammcor_Rkab_2(RII, Xga(:,1:NI), Xga(:,1:NI), Zgk, NI, NI, NChol, NTHC)
                    call tmsg('TIME FOR RII', timer, tdebug)
              end if
              if (TabChol(2) == 1) then
                    allocate(RIA(NChol, NA, NI))
                    call CalcMem3(RIA, 'RIA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIA, Xga(:,NI+1:NIA), Xga(:,1:NI), Zgk, NA, NI, NChol, NTHC)
                    call tmsg('TIME FOR RIA', timer, tdebug)
              end if
              if (TabChol(3) == 1) then
                    allocate(RAA(NChol, NA, NA))
                    call CalcMem3(RAA, 'RAA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RAA, Xga(:,NI+1:NIA), Xga(:,NI+1:NIA), Zgk, NA, NA, NChol, NTHC)
                    call tmsg('TIME FOR RAA', timer, tdebug)
              end if
              if (TabChol(4) == 1) then
                    allocate(RIV(NChol, NV, NI))
                    call CalcMem3(RIV, 'RIV')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIV, Xga(:,NIA+1:NBasis), Xga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
                    call tmsg('TIME FOR RIV', timer, tdebug)
              end if
              if (TabChol(5) == 1) then
                    allocate(RVA(NChol, NA, NV))
                    call CalcMem3(RVA, 'RVA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RVA, Xga(:,NI+1:NIA), Xga(:,NIA+1:NBasis), Zgk, NA, NV, NChol, NTHC)
                    call tmsg('TIME FOR RVA', timer, tdebug)
              end if

              call clock_start(timer0)

              caseloop: do ii = 1, 6
                    if (loops(ii)%run_main)then
                          select case(ii)
                          case(loopIA)
                                do a = 1, NI
                                      if (loops(ii)%run_inner(2)) then ! (IA|IA)
                                            do b = 1, NI
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  t_IAIA_int = t_IAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_aoao kernel
                                                  t_IAIA_ker = t_IAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          case(loopAA)
                                do a = NI+1, NIA
                                      if (loops(ii)%run_inner(2)) then ! (AA|IA)
                                            do b = 1, NI
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  t_AAIA_int = t_AAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_aaao kernel
                                                  t_AAIA_ker = t_AAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                end do
                          case(loopIV)
                                do a = 1, NI
                                      if (loops(ii)%run_inner(2)) then ! (IV|IA)
                                            do b = 1, NI
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                  t_IVIA_int = t_IVIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_voao kernel
                                                  t_IVIA_ker = t_IVIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(3)) then ! (IV|AA)
                                            do b = NI+1, NIA
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                  t_IVAA_int = t_IVAA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_voaa kernel
                                                  t_IVAA_ker = t_IVAA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(4)) then ! (IV|IV)
                                            ! b <= a: slab (b,a) gives the same E_vovo as (a,b)
                                            do b = 1, a
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
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
                                      do b = 1, NI
                                            call clock_start(timer)
                                            call real_aTb_x(V_xay, NA*NV, RVA, NChol, RIV(:,:,b), NChol, NA*NV, NV, NChol, AcAlpha)
                                            t_VAIV_int = t_VAIV_int + clock_readwall(timer)
                                            call clock_start(timer)
                                            call vavo_contr(E_vavo, ACB(phA0va), ACB(phA0vo), vecYb%va, pos, b, AuxData, NA, NV, V_xay)
                                            t_VAIV_ker = t_VAIV_ker + clock_readwall(timer)
                                      end do
                                end if
                                do a = NIA+1, NBasis
                                      if (loops(ii)%run_inner(2)) then ! (VA|IA)
                                            do b = 1, NI
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                  t_VAIA_int = t_VAIA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_vaao kernel
                                                  t_VAIA_ker = t_VAIA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(3)) then ! (VA|AA)
                                            do b = NI+1, NIA
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                  t_VAAA_int = t_VAAA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  ! E_vaaa kernel
                                                  t_VAAA_ker = t_VAAA_ker + clock_readwall(timer)
                                            end do
                                      end if
                                      if (loops(ii)%run_inner(5)) then ! (VA|VA)
                                            ! b <= a: slab (b,a) gives the same E_vava as (a,b)
                                            do b = NIA+1, a
                                                  call clock_start(timer)
                                                  call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RVA(:,:,b-NIA), NChol, NA, NA, NChol, AcAlpha)
                                                  t_VAVA_int = t_VAVA_int + clock_readwall(timer)
                                                  call clock_start(timer)
                                                  call vava_contr(E_vava, ACB(phA0va), vecXp%va, vecYb%va, rdm2, &
                                                        AuxData%map, a, b, NIA, NA, V_axby)
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

              call print_section('Final energies')
              call print_energy('CASSCF energy (one-electron)', AuxData%ECAS_oneelectr, &
                    ttoeV * AuxData%ECAS_oneelectr)
              call print_energy('ECASSCF_THC', AuxData%ECAS_THC, ttoeV * AuxData%ECAS_THC)
              call print_energy('E_phAC0', E_corr, ttoeV * E_corr)
              call print_energy('E_total', AuxData%ECAS_THC + E_corr, &
                    ttoeV * (AuxData%ECAS_THC + E_corr))
              print*, ''

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
      subroutine build_vecXt(vecXt, OA, AA, NA, NI)
            type(TvecXt), intent(out) :: vecXt
            type(TAC0Block_ph), intent(in) :: OA, AA
            integer, intent(in) :: NA, NI
            integer :: s, m, nd

            allocate(vecXt%ao(NA, NA, NI))
            vecXt%ao = zero
            do s = 1, NI
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
            integer :: k, m, jj, i, x, r, s

            associate(Occ=>AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, nd=>AA%NDim)

              allocate(c(NBasis))
              do i = 1, NBasis
                    c(i) = sqrt(Occ(i))
                    if (Occ(i) .lt. frac12) c(i) = -c(i)
              end do

              ! ---- ao: per inactive k, rows r - NI with r in ActList ----
              call move_alloc(vecXt%ao, vecXp%ao)
              allocate(vecYm%ao(NA, NA, NI), vecYb%ao(NA, NA, NI))
              vecYm%ao = zero
              vecYb%ao = zero
              do k = 1, NI
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
      ! =============================================================================
      subroutine vovo_contr(E_vovo, VO, pos, a, b, x0, x1, V_axby)
            double precision, intent(inout) :: E_vovo
            type(TAC0Block_ph), intent(in) :: VO
            integer, dimension(:,:), intent(in) :: pos
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1
            double precision, dimension(x0:x1, x0:x1), intent(in) :: V_axby
            double precision, dimension(x0:x1) :: om_a, om_b
            double precision :: val
            integer :: x, y

            ! vo has no filtered pairs, pos(x,a) > 0 always
            do x = x0, x1
                  om_a(x) = VO%Eig(pos(x, a))
                  om_b(x) = VO%Eig(pos(x, b))
            end do

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

      end subroutine vovo_contr


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
      ! =============================================================================
      subroutine vavo_contr(E_vavo, VA, VO, Yb, pos, s, AuxData, NA, NV, V_xay)
            double precision, intent(inout) :: E_vavo
            type(TAC0Block_ph), intent(in) :: VA, VO
            integer, intent(in) :: s, NA, NV
            double precision, dimension(NA, NA, NV), intent(in) :: Yb
            integer, dimension(:,:), intent(in) :: pos
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(NA, NV, NV), intent(in) :: V_xay
            double precision, dimension(:,:), allocatable :: U, T
            double precision, dimension(:), allocatable :: om_s
            double precision :: val
            integer :: p, r, m, jj

            associate(NI=>AuxData%NI, NIA=>AuxData%NIA)

              allocate(om_s(NV), U(NA, NV), T(NA, NV))

              ! vo has no filtered pairs, pos(r,s) > 0 always
              do r = 1, NV
                    om_s(r) = VO%Eig(pos(NIA+r, s))
              end do

              val = zero
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
                    end associate
              end do

              E_vavo = E_vavo + four * val

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
      ! =============================================================================
      subroutine vava_contr(E_vava, VA, Xp, Yb, rdm2, map, a, b, NIA, NA, V_axby)
            double precision, intent(inout) :: E_vava
            type(TAC0Block_ph), intent(in) :: VA
            integer, intent(in) :: a, b, NIA, NA
            double precision, dimension(:,:,:), intent(in) :: Xp, Yb
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, dimension(:), intent(in) :: map
            double precision, dimension(NA, NA), intent(in) :: V_axby
            double precision, dimension(NA, NA) :: W, M1, K, M2
            double precision :: val
            integer :: p, r, c, d, mu, nu, NI

            p = a - NIA
            r = b - NIA
            NI = NIA - NA

            associate(MBp=>VA%MiniBlocks(p), MBr=>VA%MiniBlocks(r))

              if (MBp%NDim == 0 .or. MBr%NDim == 0) return

              ! M1 = Yb^(p)T V Yb^(r)
              call real_ab(W, V_axby, Yb(:,:,r))
              call real_aTb(M1, Yb(:,:,p), W)

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

            end associate

            if (a /= b) val = two * val
            E_vava = E_vava - four * val

      end subroutine vava_contr

end module thc_energy_ph
