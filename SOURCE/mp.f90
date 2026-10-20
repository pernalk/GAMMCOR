module mp

      ! =============================================================================
      ! MP2, ROMP2 and their CBS corrections with THC integrals (JobType MP2; cbs 0, 1
      ! = CBS[DFT], 2 = CBS[H]) on the HF/ROHF dump of the PySCF interface.
      ! Frozen core: orbitals 1..NCoreOrb stay in the reference (fij, Fock) and are
      ! left out of the pairs (a, b run NI0..NI).
      ! Code: .secret_sauce/.skills/mp2/code-mp2.md; theory: .secret_sauce/.skills/mp2/theory-mp2.md.
      ! =============================================================================

      use ppac_types
      use types
      use math_constants
      use real_linalg !from gammcor integrals
      use THC_Gammcor
      use sort  !from gammcor integrals
      use clock
      use print_utils

      implicit none

      integer, parameter, private :: tdebug = 0
      double precision, parameter, private :: ttoeV = 27.211399d+0

      ! The six R vector sets of the dressed slab for one spin, all
      ! (NChol or NCholErf, virtual of the spin, occupied of the spin).
      ! 1: mu-weighted virtual orbital, 2: mu-weighted occupied orbital.
      type, private :: TDressedR
            double precision, dimension(:,:,:), allocatable :: FR, FR1, FR2, LR, LR1, LR2
      end type TDressedR

contains

      subroutine mp2_driver(THCData, AuxData, CAONO, Flags)

            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData

            double precision :: E_mp2, E_srmp2, AvMU, ECorrMD
            double precision, dimension(:,:), allocatable :: XMuMat
            double precision :: E_aa, E_bb, E_ab, E_1, E_doubles
            double precision :: D_aa, D_bb, D_ab, D_doubles, E_sromp2
            double precision, dimension(:,:), allocatable :: Xga_a, Xga_b, fov_a, fov_b, U_a, U_b
            double precision, dimension(:), allocatable :: e_a, e_b

            if (.not. allocated(THCData%Xga)) then
                  print*, 'ERROR: MP2 needs TwoMoInt FOFO with Cholesky THC'
                  stop
            end if

            select case (Flags%IDBBSC)
            case (0)
                  if (AuxData%NA > 0) then
                        call romp2_semican(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b)
                        call romp2_energy_thc(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b, &
                              E_aa, E_bb, E_ab, E_1)
                        E_doubles = E_aa + E_bb + E_ab
                        E_mp2 = E_doubles + E_1

                        call print_section('Final energies')
                        call print_energy('E_HF (ROHF)', AuxData%EROHF, ttoeV * AuxData%EROHF)
                        call print_energy('E_aa (doubles, alpha-alpha)', E_aa, ttoeV * E_aa)
                        call print_energy('E_bb (doubles, beta-beta)', E_bb, ttoeV * E_bb)
                        call print_energy('E_ab (doubles, alpha-beta)', E_ab, ttoeV * E_ab)
                        call print_energy('E_doubles', E_doubles, ttoeV * E_doubles)
                        call print_energy('E_1 (singles)', E_1, ttoeV * E_1)
                        call print_energy('E_ROMP2', E_mp2, ttoeV * E_mp2)
                        call print_energy('E_total', AuxData%EROHF + E_mp2, &
                              ttoeV * (AuxData%EROHF + E_mp2))
                        write(*, '(/1X,"HF+ENuc, ROMP2-Corr   , Total",6X,3F15.8)') &
                              AuxData%EROHF, E_mp2, AuxData%EROHF + E_mp2
                        print*, ''
                        return
                  end if

                  call mp2_energy_thc(THCData, AuxData, E_mp2)

                  call print_section('Final energies')
                  call print_energy('E_HF', AuxData%EROHF, ttoeV * AuxData%EROHF)
                  call print_energy('E_MP2', E_mp2, ttoeV * E_mp2)
                  if (AuxData%HasEMP2_pyscf) then
                        call print_energy('E_MP2 (PySCF reference)', AuxData%EMP2_pyscf, &
                              ttoeV * AuxData%EMP2_pyscf)
                  end if
                  call print_energy('E_total', AuxData%EROHF + E_mp2, &
                        ttoeV * (AuxData%EROHF + E_mp2))
                  write(*, '(/1X,"HF+ENuc, MP2-Corr     , Total",6X,3F15.8)') &
                        AuxData%EROHF, E_mp2, AuxData%EROHF + E_mp2
                  print*, ''

            case (1)
                  call loc_mu_thc(THCData, AuxData, CAONO, Flags, XMuMat, AvMU, ECorrMD)

                  if (AuxData%NA > 0) then
                        call romp2_semican(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b)
                        call romp2_energy_thc(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b, &
                              E_aa, E_bb, E_ab, E_1)
                        E_doubles = E_aa + E_bb + E_ab
                        E_mp2 = E_doubles + E_1

                        call print_section('Final energies')
                        call print_energy('E_HF (ROHF)', AuxData%EROHF, ttoeV * AuxData%EROHF)
                        call print_energy('E_aa (doubles, alpha-alpha)', E_aa, ttoeV * E_aa)
                        call print_energy('E_bb (doubles, beta-beta)', E_bb, ttoeV * E_bb)
                        call print_energy('E_ab (doubles, alpha-beta)', E_ab, ttoeV * E_ab)
                        call print_energy('E_doubles', E_doubles, ttoeV * E_doubles)
                        call print_energy('E_1 (singles)', E_1, ttoeV * E_1)
                        call print_energy('E_ROMP2', E_mp2, ttoeV * E_mp2)
                        call print_energy('CBS[DFT] correction', ECorrMD, ttoeV * ECorrMD)
                        call print_energy('E_ROMP2-CBS[DFT]', E_mp2 + ECorrMD, ttoeV * (E_mp2 + ECorrMD))
                        call print_energy('E_total', AuxData%EROHF + E_mp2 + ECorrMD, &
                              ttoeV * (AuxData%EROHF + E_mp2 + ECorrMD))
                        write(*, '(/1X,"HF+ENuc, ROMP2-Corr   , Total",6X,3F15.8)') &
                              AuxData%EROHF, E_mp2, AuxData%EROHF + E_mp2
                        write(*, '(1X,"HF+ENuc, ROMP2-CBS[DFT], Total",5X,3F15.8)') &
                              AuxData%EROHF, E_mp2 + ECorrMD, AuxData%EROHF + E_mp2 + ECorrMD
                        print*, ''
                        return
                  end if

                  call mp2_energy_thc(THCData, AuxData, E_mp2)

                  call print_section('Final energies')
                  call print_energy('E_HF', AuxData%EROHF, ttoeV * AuxData%EROHF)
                  call print_energy('E_MP2', E_mp2, ttoeV * E_mp2)
                  call print_energy('CBS[DFT] correction', ECorrMD, ttoeV * ECorrMD)
                  call print_energy('E_MP2-CBS[DFT]', E_mp2 + ECorrMD, ttoeV * (E_mp2 + ECorrMD))
                  call print_energy('E_total', AuxData%EROHF + E_mp2 + ECorrMD, &
                        ttoeV * (AuxData%EROHF + E_mp2 + ECorrMD))
                  write(*, '(/1X,"HF+ENuc, MP2-Corr     , Total",6X,3F15.8)') &
                        AuxData%EROHF, E_mp2, AuxData%EROHF + E_mp2
                  write(*, '(1X,"HF+ENuc, MP2-CBS[DFT], Total",6X,3F15.8)') &
                        AuxData%EROHF, E_mp2 + ECorrMD, AuxData%EROHF + E_mp2 + ECorrMD
                  print*, ''

            case (2)
                  if (AuxData%NA > 0) then
                        call romp2_semican(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b, U_a, U_b)
                        call sromp2_energy_thc(THCData, AuxData, CAONO, Flags, Xga_a, Xga_b, U_a, U_b, &
                              e_a, e_b, fov_a, fov_b, E_aa, E_bb, E_ab, E_1, D_aa, D_bb, D_ab)
                        E_doubles = E_aa + E_bb + E_ab
                        D_doubles = D_aa + D_bb + D_ab
                        E_mp2 = E_doubles + E_1
                        E_sromp2 = D_doubles + E_1

                        call print_section('Final energies')
                        call print_energy('E_HF (ROHF)', AuxData%EROHF, ttoeV * AuxData%EROHF)
                        call print_energy('E_aa (doubles, FR)', E_aa, ttoeV * E_aa)
                        call print_energy('E_bb (doubles, FR)', E_bb, ttoeV * E_bb)
                        call print_energy('E_ab (doubles, FR)', E_ab, ttoeV * E_ab)
                        call print_energy('E_aa (doubles, CBS[H])', D_aa, ttoeV * D_aa)
                        call print_energy('E_bb (doubles, CBS[H])', D_bb, ttoeV * D_bb)
                        call print_energy('E_ab (doubles, CBS[H])', D_ab, ttoeV * D_ab)
                        call print_energy('E_1 (singles, FR)', E_1, ttoeV * E_1)
                        call print_energy('E_ROMP2', E_mp2, ttoeV * E_mp2)
                        call print_energy('E_ROMP2-CBS[H]', E_sromp2, ttoeV * E_sromp2)
                        call print_energy('CBS[H] correction', E_sromp2 - E_mp2, ttoeV * (E_sromp2 - E_mp2))
                        call print_energy('E_total', AuxData%EROHF + E_sromp2, &
                              ttoeV * (AuxData%EROHF + E_sromp2))
                        write(*, '(/1X,"HF+ENuc, ROMP2-Corr   , Total",6X,3F15.8)') &
                              AuxData%EROHF, E_mp2, AuxData%EROHF + E_mp2
                        write(*, '(1X,"HF+ENuc, ROMP2-CBS[H]-Corr, Total",2X,3F15.8)') &
                              AuxData%EROHF, E_sromp2, AuxData%EROHF + E_sromp2
                        print*, ''
                        return
                  end if

                  call srmp2_energy_thc(THCData, AuxData, CAONO, Flags, E_mp2, E_srmp2)

                  call print_section('Final energies')
                  call print_energy('E_HF', AuxData%EROHF, ttoeV * AuxData%EROHF)
                  call print_energy('E_MP2', E_mp2, ttoeV * E_mp2)
                  call print_energy('E_MP2-CBS[H]', E_srmp2, ttoeV * E_srmp2)
                  call print_energy('CBS[H] correction', E_srmp2 - E_mp2, ttoeV * (E_srmp2 - E_mp2))
                  call print_energy('E_total', AuxData%EROHF + E_srmp2, &
                        ttoeV * (AuxData%EROHF + E_srmp2))
                  write(*, '(/1X,"HF+ENuc, MP2-CBS[H]-Corr, Total",4X,3F15.8)') &
                        AuxData%EROHF, E_srmp2, AuxData%EROHF + E_srmp2
                  print*, ''

            case default
                  print*, 'ERROR: MP2 supports cbs 0, 1 or 2, got ', Flags%IDBBSC
                  stop
            end select

      end subroutine mp2_driver


      ! =============================================================================
      ! Closed-shell MP2 correlation energy = E_vovo of THC_energy_loop_ph with
      ! omega_(v,i) = fvw(v) - fij(i):
      !
      !   E_mp2 = sum_{pq, rs}  (pq|rs) [ (qr|ps) - 2 (pq|rs) ] / (omega_pq + omega_rs)
      !
      ! p, r virtual, q, s inactive NI0..NI. One (IV|IV) slab per (a, b),
      ! V_axby(x,y) = (x a|y b) from RIV, b <= a (vovo_slab doubles a /= b).
      ! =============================================================================
      subroutine mp2_energy_thc(THCData, AuxData, E_mp2)

            type(TTHCData), intent(in) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, intent(out) :: E_mp2

            double precision, dimension(:), allocatable :: V_axby
            double precision, dimension(:,:,:), allocatable :: RIV
            double precision, dimension(:), allocatable :: w_a, w_b
            double precision, parameter :: AcAlpha = One
            integer :: NI, NV, NIA, NBasis, NI0
            integer :: a, b, x
            type (tclock) :: timer

            associate(NChol=>THCData%NChol)

              call print_section('MP2 energy evaluation')

              NI = AuxData%NI
              NV = AuxData%NV
              NIA = AuxData%NIA
              NBasis = AuxData%NBasis
              NI0 = AuxData%NCoreOrb + 1
              call print_info('Frozen core orbitals', NI0 - 1)
              call print_info('Correlated inactive orbitals', NI - NI0 + 1)
              call print_info('Virtual orbitals', NV)

              allocate(RIV(NChol, NV, NI0:NI))
              call CalcMem3(RIV, 'RIV')
              call clock_start(timer)
              call thc_gammcor_Rkab_2(RIV, THCData%Xga(:,NIA+1:NBasis), THCData%Xga(:,NI0:NI), &
                    THCData%Zgk, NV, NI-NI0+1, NChol, THCData%NTHC)
              call tmsg('TIME FOR RIV', timer, tdebug)

              allocate(V_axby(NV**2))
              allocate(w_a(NIA+1:NBasis), w_b(NIA+1:NBasis))

              call clock_start(timer)
              E_mp2 = zero
              do a = NI0, NI
                    do x = NIA+1, NBasis
                          w_a(x) = THCData%fvw(x-NIA) - THCData%fij(a)
                    end do
                    ! b <= a: slab (b,a) gives the same sum as (a,b)
                    do b = NI0, a
                          do x = NIA+1, NBasis
                                w_b(x) = THCData%fvw(x-NIA) - THCData%fij(b)
                          end do
                          call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                          call vovo_slab(E_mp2, w_a, w_b, a, b, NIA+1, NBasis, V_axby)
                    end do
              end do
              call tmsg('TIME FOR MP2 energy loop', timer, tdebug)

            end associate

      end subroutine mp2_energy_thc


      ! =============================================================================
      ! The MP2 sum for one (IV|IV) slab, the same sum as vovo_contr
      ! (thc_energy_ph.f90) but with the omegas given:
      !
      !   E_vovo += sum_{x,y} V_axby(x,y) [V_axby(y,x) - 2 V_axby(x,y)] / (om_a(x) + om_b(y))
      !
      ! om_a(x) = fvw(x) - fij(a), om_b(y) = fvw(y) - fij(b). Doubled for
      ! a /= b, because the caller runs b <= a only.
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


      ! ROMP2 part 1: semicanonical alpha and beta orbitals of the ROHF reference.
      ! The Fock sums start at 1, not NI0: frozen core is not correlated but is in
      ! the Fock operator. See .secret_sauce/.skills/mp2/code-mp2.md Sec. 4.
      subroutine romp2_semican(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b, Uout_a, Uout_b)
            !use mp_ext, only: mp_fock_jk
            use THCFock

            type(TTHCData), intent(in) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), allocatable, intent(out) :: Xga_a, Xga_b
            double precision, dimension(:), allocatable, intent(out) :: e_a, e_b
            double precision, dimension(:,:), allocatable, intent(out) :: fov_a, fov_b
            double precision, dimension(:,:), allocatable, intent(out), optional :: Uout_a, Uout_b

            double precision, dimension(:,:,:), allocatable :: JK, Cpi
            double precision, dimension(:,:), allocatable :: Zgh, f_a, f_b, U_a, U_b
            integer, dimension(2) :: NOcc
            integer :: NI, NIA, NBasis, NI0, i

            call print_section('ROMP2 semicanonical orbitals')

            NI = AuxData%NI
            NIA = AuxData%NIA
            NBasis = AuxData%NBasis
            NI0 = AuxData%NCoreOrb + 1

            ! occupied orbitals per spin as identity columns of the MO basis
            allocate(Cpi(NBasis, NIA, 2))
            Cpi = zero
            do i = 1, NIA
                  Cpi(i, i, 1) = one
            end do
            do i = 1, NI
                  Cpi(i, i, 2) = one
            end do
            NOcc = (/ NIA, NI /)

            allocate(Zgh(THCData%NTHC, THCData%NTHC))
            call real_abT(Zgh, THCData%Zgk, THCData%Zgk)

            ! JK(:,:,s) = J[D_a + D_b] - K[D_s]
            allocate(JK(NBasis, NBasis, 2))
            !call mp_fock_jk(JK, Cpi, Zgh, THCData%Xga, NOcc)
            call thc_Fock_JK(JK, Cpi, Zgh, THCData%Xga, NOcc, .true., .true., 1d0)

            allocate(f_a(NBasis, NBasis), f_b(NBasis, NBasis))
            f_a = AuxData%HNO0 + JK(:,:,1)
            f_b = AuxData%HNO0 + JK(:,:,2)

            allocate(U_a(NBasis, NBasis), U_b(NBasis, NBasis))
            allocate(e_a(NBasis), e_b(NBasis))
            call semican_block(f_a, U_a, e_a, NI0, NIA, NBasis)
            call semican_block(f_b, U_b, e_b, NI0, NI, NBasis)

            allocate(Xga_a(THCData%NTHC, NBasis), Xga_b(THCData%NTHC, NBasis))
            call real_ab(Xga_a, THCData%Xga, U_a)
            call real_ab(Xga_b, THCData%Xga, U_b)

            ! f in the rotated basis, U^T f U; keep the occupied-virtual block
            call rotate_fock(f_a, U_a, NBasis)
            call rotate_fock(f_b, U_b, NBasis)
            allocate(fov_a(NI0:NIA, NIA+1:NBasis), fov_b(NI0:NI, NI+1:NBasis))
            fov_a = f_a(NI0:NIA, NIA+1:NBasis)
            fov_b = f_b(NI0:NI, NI+1:NBasis)

            call print_info('Frozen core orbitals', NI0 - 1)
            call print_info('alpha occ/vir', (/ NIA - NI0 + 1, NBasis - NIA /))
            call print_info('beta  occ/vir', (/ NI - NI0 + 1, NBasis - NI /))
            if (NIA >= NI0) call print_info('e_a highest occupied', e_a(NIA))
            if (NIA < NBasis) call print_info('e_a lowest virtual', e_a(NIA+1))
            if (NI >= NI0) call print_info('e_b highest occupied', e_b(NI))
            if (NI < NBasis) call print_info('e_b lowest virtual', e_b(NI+1))
            call print_info('||fov_a||', norm2(fov_a))
            call print_info('||fov_b||', norm2(fov_b))

            if (present(Uout_a)) call move_alloc(U_a, Uout_a)
            if (present(Uout_b)) call move_alloc(U_b, Uout_b)

      contains

            ! U = 1 and e = diag(f), then diagonalize f on the occupied block
            ! occ0..occ1 and on the virtual block occ1+1..n.
            subroutine semican_block(f, U, e, occ0, occ1, n)
                  double precision, dimension(:,:), intent(in) :: f
                  double precision, dimension(:,:), intent(out) :: U
                  double precision, dimension(:), intent(out) :: e
                  integer, intent(in) :: occ0, occ1, n
                  integer :: k

                  U = zero
                  do k = 1, n
                        U(k, k) = one
                        e(k) = f(k, k)
                  end do
                  call diag_block(f, U, e, occ0, occ1)
                  call diag_block(f, U, e, occ1+1, n)
            end subroutine semican_block

            subroutine diag_block(f, U, e, p0, p1)
                  double precision, dimension(:,:), intent(in) :: f
                  double precision, dimension(:,:), intent(inout) :: U
                  double precision, dimension(:), intent(inout) :: e
                  integer, intent(in) :: p0, p1
                  double precision, dimension(:,:), allocatable :: W
                  double precision, dimension(:), allocatable :: w_eig
                  integer :: m

                  m = p1 - p0 + 1
                  if (m <= 0) return
                  allocate(W(m, m), w_eig(m))
                  ! symmetrize: f is symmetric up to the THC error
                  W = (f(p0:p1, p0:p1) + transpose(f(p0:p1, p0:p1))) / two
                  call symmetric_eigenproblem(w_eig, W, m, .true.)
                  U(p0:p1, p0:p1) = W
                  e(p0:p1) = w_eig
            end subroutine diag_block

            subroutine rotate_fock(f, U, n)
                  double precision, dimension(:,:), intent(inout) :: f
                  double precision, dimension(:,:), intent(in) :: U
                  integer, intent(in) :: n
                  double precision, dimension(:,:), allocatable :: W

                  allocate(W(n, n))
                  call real_ab(W, f, U)
                  call real_aTb(f, U, W)
            end subroutine rotate_fock

      end subroutine romp2_semican


      ! ROMP2 part 2: E_aa, E_bb, E_ab and E_1 in the semicanonical orbitals. Same
      ! spin over j < i without 1/2. Check: NA = 0 reproduces mp2_energy_thc.
      ! See .secret_sauce/.skills/mp2/code-mp2.md Sec. 4, theory-mp2.md Sec. 3.
      subroutine romp2_energy_thc(THCData, AuxData, Xga_a, Xga_b, e_a, e_b, fov_a, fov_b, &
            E_aa, E_bb, E_ab, E_1)

            type(TTHCData), intent(in) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: Xga_a, Xga_b
            double precision, dimension(:), intent(in) :: e_a, e_b
            ! allocatable, so that the global index bounds set by romp2_semican are kept
            double precision, dimension(:,:), allocatable, intent(in) :: fov_a, fov_b
            double precision, intent(out) :: E_aa, E_bb, E_ab, E_1

            double precision, dimension(:,:,:), allocatable :: RA, RB
            double precision, dimension(:,:), allocatable :: V
            double precision, dimension(:), allocatable :: w_i, w_j
            integer :: NI, NA, NV, NIA, NBasis, NI0
            integer :: i, j
            type (tclock) :: timer

            associate(NChol=>THCData%NChol)

              call print_section('ROMP2 energy evaluation')

              NI = AuxData%NI
              NA = AuxData%NA
              NV = AuxData%NV
              NIA = AuxData%NIA
              NBasis = AuxData%NBasis
              NI0 = AuxData%NCoreOrb + 1

              ! RA(K,A,i): alpha, A = NIA+1..NBasis, i = NI0..NIA
              ! RB(K,B,j): beta,  B = NI+1..NBasis,  j = NI0..NI
              allocate(RA(NChol, NV, NI0:NIA), RB(NChol, NA+NV, NI0:NI))
              call CalcMem3(RA, 'RA')
              call CalcMem3(RB, 'RB')
              call clock_start(timer)
              call thc_gammcor_Rkab_2(RA, Xga_a(:,NIA+1:NBasis), Xga_a(:,NI0:NIA), &
                    THCData%Zgk, NV, NIA-NI0+1, NChol, THCData%NTHC)
              call thc_gammcor_Rkab_2(RB, Xga_b(:,NI+1:NBasis), Xga_b(:,NI0:NI), &
                    THCData%Zgk, NA+NV, NI-NI0+1, NChol, THCData%NTHC)
              call tmsg('TIME FOR RA, RB', timer, tdebug)

              call clock_start(timer)

              ! opposite spin: every alpha i with every beta j
              allocate(V(NV, NA+NV), w_i(NV), w_j(NA+NV))
              E_ab = zero
              do i = NI0, NIA
                    w_i = e_a(NIA+1:NBasis) - e_a(i)
                    do j = NI0, NI
                          w_j = e_b(NI+1:NBasis) - e_b(j)
                          call real_aTb_x(V, NV, RA(:,:,i), NChol, RB(:,:,j), NChol, NV, NA+NV, NChol)
                          call romp2_os_slab(E_ab, w_i, w_j, V)
                    end do
              end do
              deallocate(V, w_i, w_j)

              call same_spin(E_aa, RA, e_a, NIA+1)
              call same_spin(E_bb, RB, e_b, NI+1)

              ! singles
              E_1 = zero
              do i = NI0, NIA
                    E_1 = E_1 - sum(fov_a(i,:)**2 / (e_a(NIA+1:NBasis) - e_a(i)))
              end do
              do i = NI0, NI
                    E_1 = E_1 - sum(fov_b(i,:)**2 / (e_b(NI+1:NBasis) - e_b(i)))
              end do

              call tmsg('TIME FOR ROMP2 energy loop', timer, tdebug)

            end associate

      contains

            ! E_ss for one spin. R(K,A,i) holds that spin's vectors, with i
            ! running over its occupied orbitals (the bounds of R, dim 3) and A
            ! over its virtual orbitals v0..NBasis; e are its orbital energies.
            subroutine same_spin(E_ss, R, e, v0)
                  double precision, intent(out) :: E_ss
                  double precision, dimension(:,:,:), allocatable, intent(in) :: R
                  double precision, dimension(:), intent(in) :: e
                  integer, intent(in) :: v0
                  double precision, dimension(:,:), allocatable :: W
                  double precision, dimension(:), allocatable :: w_i, w_j
                  integer :: nvir, nk, p, q

                  nvir = NBasis - v0 + 1
                  nk = size(R, 1)
                  allocate(W(nvir, nvir), w_i(nvir), w_j(nvir))
                  E_ss = zero
                  do p = lbound(R, 3), ubound(R, 3)
                        w_i = e(v0:NBasis) - e(p)
                        do q = lbound(R, 3), p - 1
                              w_j = e(v0:NBasis) - e(q)
                              call real_aTb_x(W, nvir, R(:,:,p), nk, R(:,:,q), nk, nvir, nvir, nk)
                              call romp2_ss_slab(E_ss, w_i, w_j, W)
                        end do
                  end do
            end subroutine same_spin

      end subroutine romp2_energy_thc


      ! =============================================================================
      ! Opposite-spin ROMP2 sum for one alpha occupied i and one beta occupied j:
      !
      !   E += - sum_{A,B} V(A,B)^2 / (w_i(A) + w_j(B))
      !
      ! V(A,B) = (iA|jB), A alpha virtual, B beta virtual,
      ! w_i(A) = e_a(A) - e_a(i), w_j(B) = e_b(B) - e_b(j).
      ! =============================================================================
      subroutine romp2_os_slab(E, w_i, w_j, V)
            double precision, intent(inout) :: E
            double precision, dimension(:), intent(in) :: w_i, w_j
            double precision, dimension(:,:), intent(in) :: V
            double precision :: val
            integer :: x, y

            val = zero
            do y = 1, size(w_j)
                  do x = 1, size(w_i)
                        val = val + V(x,y)**2 / (w_i(x) + w_j(y))
                  end do
            end do

            E = E - val

      end subroutine romp2_os_slab


      ! Same-spin ROMP2 sum for one pair i > j:
      ! E += - sum_{A,B} V(A,B) [V(A,B) - V(B,A)] / (w_i(A) + w_j(B)).
      subroutine romp2_ss_slab(E, w_i, w_j, V)
            double precision, intent(inout) :: E
            double precision, dimension(:), intent(in) :: w_i, w_j
            double precision, dimension(:,:), intent(in) :: V
            double precision :: val
            integer :: x, y

            val = zero
            do y = 1, size(w_j)
                  do x = 1, size(w_i)
                        val = val + V(x,y) * (V(x,y) - V(y,x)) / (w_i(x) + w_j(y))
                  end do
            end do

            E = E - val

      end subroutine romp2_ss_slab


      ! Local mu(r) and XMuMat for an HF/ROHF reference (NAct = 0 case of
      ! LOC_MU_CBS_CHOL, with THC factors). Frozen core left out (rows/cols of
      ! XMuMat zero). CAONO must be the orbitals of Xga. ECorrMD = CBS[DFT].
      ! See .secret_sauce/.skills/mp2/code-mp2.md, theory-mp2.md Sec. 5.
      subroutine loc_mu_thc(THCData, AuxData, CAONO, Flags, XMuMat, AvMU, ECorrMD)
            !use units
            use grid_internal

            type(TTHCData), intent(in) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:,:), allocatable, intent(out) :: XMuMat
            double precision, intent(out) :: AvMU, ECorrMD

            double precision, dimension(:,:,:), allocatable :: PhiGGA
            double precision, dimension(:), allocatable :: WGrid
            double precision, dimension(:), allocatable :: FPsiB, OnTop, XMuLoc
            double precision, dimension(:), allocatable :: RhoGrid, Sigma, Zk
            double precision, dimension(:), allocatable :: RhoA, RhoB
            double precision, dimension(:,:), allocatable :: phib, A, B, t, ta, work
            double precision :: rho_a, rho_b, gx_a, gy_a, gz_a, gx_b, gy_b, gz_b
            double precision :: XMu, XEl, Bet, OnTopC, SPi, Const
            double precision :: XMuNegMin, XMuNegMax, XElNeg
            integer :: NI, NIA, NBasis, NI0, NGrid, nb, r0, r1, n, r, i, q, ir, NMuNeg
            type (tclock) :: timer

            double precision, parameter :: ThrOnTop = 1d-8
            double precision, parameter :: BigE = 1.d2
            ! grid points per batch: A, B are NTHC x nb
            integer, parameter :: MaxBatchSize = 2**24

            associate(Xga=>THCData%Xga, Zgk=>THCData%Zgk, NChol=>THCData%NChol, NTHC=>THCData%NTHC)

              call print_section('CBS Correction with local Mu (THC)')
              call clock_start(timer)

              NI = AuxData%NI
              NIA = AuxData%NIA
              NBasis = AuxData%NBasis
              NI0 = AuxData%NCoreOrb + 1
              call print_info('Frozen core orbitals', NI0 - 1)
              if (NIA > NI) call print_info('Singly occupied orbitals', NIA - NI)

              ! orbitals and gradients on the grid, PhiGGA(r, p, 1:4)
              !call mp_gga_orbgrid(Flags, CAONO, WGrid, PhiGGA, NGrid, NBasis)

              ! Grid in input.inp: SG1, MEDIUM (default), FINE or XFINE
              if (Flags%IGridType > GRID_PARAMS_XFINE) then
                    print*, 'ERROR: CBS with local mu needs an internal grid (Grid SG1/MEDIUM/FINE/XFINE)'
                    stop
              end if

              call internal_gga_no_orbgrid(Flags, Flags%IGridType, Flags%ORBITAL_ORDERING, CAONO, &
                  WGrid, PhiGGA, NGrid, NBasis, NBasis, Flags%IUnits)

              call print_info('The number of Grid Points', NGrid)
              ! WGrid keeps the size of the initial Becke grid; only 1:NGrid is used

              ! f(r): FPsiB(r) = 2 sum_K t_K^a(r) t_K^b(r). t holds t^b (occupied
              ! NI0..NI); t^a = t^b + the part of the singly occupied NI+1..NIA.
              allocate(FPsiB(NGrid))
              nb = max(1, min(NGrid, MaxBatchSize / NTHC))
              allocate(phib(nb, NBasis), A(NTHC, nb), B(NTHC, nb), t(NChol, nb))
              if (NIA > NI) allocate(ta(NChol, nb))
              do r0 = 1, NGrid, nb
                    r1 = min(NGrid, r0 + nb - 1)
                    n = r1 - r0 + 1
                    phib(1:n, :) = PhiGGA(r0:r1, :, 1)
                    call real_abT(A(:,1:n), Xga(:,NI0:NI), phib(1:n,NI0:NI))
                    call real_abT(B(:,1:n), Xga, phib(1:n,:))
                    A(:,1:n) = A(:,1:n) * B(:,1:n)
                    call real_aTb(t(:,1:n), Zgk, A(:,1:n))
                    if (NIA > NI) then
                          call real_abT(A(:,1:n), Xga(:,NI+1:NIA), phib(1:n,NI+1:NIA))
                          A(:,1:n) = A(:,1:n) * B(:,1:n)
                          call real_aTb(ta(:,1:n), Zgk, A(:,1:n))
                          ta(:,1:n) = ta(:,1:n) + t(:,1:n)
                          FPsiB(r0:r1) = two * sum(ta(:,1:n) * t(:,1:n), dim=1)
                    else
                          FPsiB(r0:r1) = two * sum(t(:,1:n)**2, dim=1)
                    end if
              end do
              deallocate(phib, A, B, t)
              if (allocated(ta)) deallocate(ta)
              call print_info('Total FPsiB', norm2(FPsiB))

              ! density, its gradient, on-top density and mu(r)
              allocate(OnTop(NGrid), XMuLoc(NGrid), RhoGrid(NGrid), Sigma(NGrid))
              ! rho_s = sum_{i in occ_s} phi_i^2, g_s = sum_{i in occ_s} phi_i grad phi_i,
              ! so grad rho = 2 (g_a + g_b). Closed shell: rho_a = rho_b.
              ! NA > 0: RhoA, RhoB are kept only for the DEBUG print below.
              if (NIA > NI) allocate(RhoA(NGrid), RhoB(NGrid))
              NMuNeg = 0
              XMuNegMin = zero
              XMuNegMax = -BigE
              XElNeg = zero
              do r = 1, NGrid
                    rho_b = zero
                    gx_b = zero
                    gy_b = zero
                    gz_b = zero
                    do i = NI0, NI
                          rho_b = rho_b + PhiGGA(r,i,1)**2
                          gx_b = gx_b + PhiGGA(r,i,1) * PhiGGA(r,i,2)
                          gy_b = gy_b + PhiGGA(r,i,1) * PhiGGA(r,i,3)
                          gz_b = gz_b + PhiGGA(r,i,1) * PhiGGA(r,i,4)
                    end do
                    rho_a = rho_b
                    gx_a = gx_b
                    gy_a = gy_b
                    gz_a = gz_b
                    do i = NI+1, NIA
                          rho_a = rho_a + PhiGGA(r,i,1)**2
                          gx_a = gx_a + PhiGGA(r,i,1) * PhiGGA(r,i,2)
                          gy_a = gy_a + PhiGGA(r,i,1) * PhiGGA(r,i,3)
                          gz_a = gz_a + PhiGGA(r,i,1) * PhiGGA(r,i,4)
                    end do
                    RhoGrid(r) = rho_a + rho_b
                    Sigma(r) = (two * (gx_a + gx_b))**2 + (two * (gy_a + gy_b))**2 &
                          + (two * (gz_a + gz_b))**2
                    OnTop(r) = two * rho_a * rho_b
                    if (NIA > NI) then
                          RhoA(r) = rho_a
                          RhoB(r) = rho_b
                    end if

                    XMuLoc(r) = BigE
                    if (OnTop(r) > ThrOnTop) XMuLoc(r) = sqrt(PI) / two * FPsiB(r) / OnTop(r)
                    if (XMuLoc(r) < zero) then
                          NMuNeg = NMuNeg + 1
                          XMuNegMin = min(XMuNegMin, XMuLoc(r))
                          XMuNegMax = max(XMuNegMax, XMuLoc(r))
                          XElNeg = XElNeg + RhoGrid(r) * WGrid(r)
                    end if
                    ! AT 18.10.2026 open shell: mu < 0 where rho_b -> 0 blows up exp(-mu)
                    if (XMuLoc(r) < zero .or. XMuLoc(r) > BigE) XMuLoc(r) = BigE
              end do
              call print_info('OnTop', norm2(OnTop))
              call print_info('XMuLoc', norm2(XMuLoc))
              call print_info('Points with mu < 0 (set to BigE)', NMuNeg)
              if (NMuNeg > 0) then
                    call print_info('  their mu, min', XMuNegMin)
                    call print_info('  their mu, max', XMuNegMax)
                    call print_info('  their electrons', XElNeg)
              end if

              ! DEBUG: negative mu(r) makes exp(-mu) blow up XMuMat
              ir = minloc(XMuLoc, dim=1)
              print*, 'DEBUG points with FPsiB < 0, mu < 0:', count(FPsiB < zero), count(XMuLoc < zero)
              print*, 'DEBUG min mu, at point', XMuLoc(ir), ir
              print*, 'DEBUG   FPsiB, OnTop, WGrid:', FPsiB(ir), OnTop(ir), WGrid(ir)
              if (NIA > NI) print*, 'DEBUG   RhoA, RhoB:', RhoA(ir), RhoB(ir)
              ir = maxloc(WGrid(1:NGrid) * exp(-XMuLoc), dim=1)
              print*, 'DEBUG max WGrid*exp(-mu), mu, at point', WGrid(ir) * exp(-XMuLoc(ir)), XMuLoc(ir), ir

              ! XMuMat(p,q) = sum_r phi_p(r) phi_q(r) WGrid(r) exp(-XMuLoc(r))
              allocate(work(NGrid, NBasis))
              do q = 1, NBasis
                    work(:,q) = PhiGGA(:,q,1) * WGrid(1:NGrid) * exp(-XMuLoc)
              end do
              allocate(XMuMat(NBasis, NBasis))
              call real_aTb(XMuMat, PhiGGA(:,:,1), work)
              deallocate(work)
              XMuMat(1:NI0-1, :) = zero
              XMuMat(:, 1:NI0-1) = zero
              call print_info('XMuMAT', norm2(XMuMat))
              print*, 'DEBUG max |XMuMat| at (p,q):', maxval(abs(XMuMat)), maxloc(abs(XMuMat))

              ! CBS[DFT] correction and the average mu
              allocate(Zk(NGrid))
              call PBECor(RhoGrid, Sigma, Zk, NGrid)
              call print_info('RhoGrid', norm2(RhoGrid))
              call print_info('Zk', norm2(Zk))

              ! single-precision pi, as in LOC_MU_CBS_CHOL
              SPi = sqrt(3.141592653589793)
              Const = three / two / SPi / (one - sqrt(two))

              AvMU = zero
              ECorrMD = zero
              XEl = zero
              do r = 1, NGrid
                    XMu = XMuLoc(r)
                    AvMU = AvMU + XMu * RhoGrid(r) * WGrid(r)
                    XEl = XEl + RhoGrid(r) * WGrid(r)

                    OnTopC = zero
                    if (XMu > 1.d-8) OnTopC = OnTop(r) / (one + two / SPi / XMu)
                    if (OnTopC /= zero) then
                          Bet = Const * Zk(r) / OnTopC
                          ECorrMD = ECorrMD + Zk(r) / (one + Bet * XMu**3) * WGrid(r)
                    end if
              end do
              AvMU = AvMU / XEl

              write(*, '(/," CBS[DFT] correction, average Mu",F15.8,F15.3/)') ECorrMD, AvMU
              call tmsg('TIME FOR loc_mu_thc', timer, tdebug)

            end associate

      end subroutine loc_mu_thc


      ! MP2-CBS[H] closed shell: mp2_energy_thc with the slab dressed FR + dSR
      ! (dSR from XMuMat of loc_mu_thc). Needs the erf THC factors (cbs 2).
      ! See .secret_sauce/.skills/mp2/code-mp2.md, theory-mp2.md Sec. 5.
      subroutine srmp2_energy_thc(THCData, AuxData, CAONO, Flags, E_mp2, E_srmp2)

            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            double precision, intent(out) :: E_mp2, E_srmp2

            double precision, dimension(:,:), allocatable :: XMuMat
            double precision, dimension(:,:,:), allocatable :: RVI_FR, RVI_FR1, RVI_FR2
            double precision, dimension(:,:,:), allocatable :: RVI_LR, RVI_LR1, RVI_LR2
            double precision, dimension(:,:), allocatable :: V_FR, V_axby
            double precision, dimension(:), allocatable :: w_a, w_b
            double precision :: AvMU, ECorrMD
            integer :: NI, NV, NIA, NBasis, NI0, NOcc
            integer :: a, b, x
            type (tclock) :: timer

            if (.not. allocated(THCData%XgaErf)) then
                  print*, 'ERROR: MP2-CBS[H] needs the erf THC factors (cbs 2 in input)'
                  stop
            end if

            NI = AuxData%NI
            NV = AuxData%NV
            NIA = AuxData%NIA
            NBasis = AuxData%NBasis
            NI0 = AuxData%NCoreOrb + 1
            NOcc = NI - NI0 + 1

            call loc_mu_thc(THCData, AuxData, CAONO, Flags, XMuMat, AvMU, ECorrMD)
            call real_abT(THCData%TXga, THCData%Xga, XMuMat)
            call real_abT(THCData%TXgaErf, THCData%XgaErf, XMuMat)

            associate(Xga=>THCData%Xga, TXga=>THCData%TXga, Zgk=>THCData%Zgk, &
                  NChol=>THCData%NChol, NTHC=>THCData%NTHC, &
                  XgaErf=>THCData%XgaErf, TXgaErf=>THCData%TXgaErf, ZgkErf=>THCData%ZgkErf, &
                  NCholErf=>THCData%NCholErf, NTHCErf=>THCData%NTHCErf)

              call print_section('MP2-CBS[H] energy evaluation')
              call print_info('Frozen core orbitals', NI0 - 1)
              call print_info('Correlated inactive orbitals', NOcc)
              call print_info('Virtual orbitals', NV)

              allocate(RVI_FR(NChol, NV, NI0:NI), RVI_FR1(NChol, NV, NI0:NI), RVI_FR2(NChol, NV, NI0:NI))
              allocate(RVI_LR(NCholErf, NV, NI0:NI), RVI_LR1(NCholErf, NV, NI0:NI), RVI_LR2(NCholErf, NV, NI0:NI))
              call CalcMem3(RVI_FR, 'RVI_FR')
              call CalcMem3(RVI_LR, 'RVI_LR')

              call clock_start(timer)
              call thc_gammcor_Rkab_2(RVI_FR,  Xga(:,NIA+1:NBasis),  Xga(:,NI0:NI),  Zgk, NV, NOcc, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_FR1, TXga(:,NIA+1:NBasis), Xga(:,NI0:NI),  Zgk, NV, NOcc, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_FR2, Xga(:,NIA+1:NBasis),  TXga(:,NI0:NI), Zgk, NV, NOcc, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_LR,  XgaErf(:,NIA+1:NBasis),  XgaErf(:,NI0:NI),  ZgkErf, &
                    NV, NOcc, NCholErf, NTHCErf)
              call thc_gammcor_Rkab_2(RVI_LR1, TXgaErf(:,NIA+1:NBasis), XgaErf(:,NI0:NI),  ZgkErf, &
                    NV, NOcc, NCholErf, NTHCErf)
              call thc_gammcor_Rkab_2(RVI_LR2, XgaErf(:,NIA+1:NBasis),  TXgaErf(:,NI0:NI), ZgkErf, &
                    NV, NOcc, NCholErf, NTHCErf)
              call tmsg('TIME FOR RVI_FR, RVI_LR', timer, tdebug)

              allocate(V_FR(NV, NV), V_axby(NV, NV))
              allocate(w_a(NIA+1:NBasis), w_b(NIA+1:NBasis))

              call clock_start(timer)
              E_mp2 = zero
              E_srmp2 = zero
              do a = NI0, NI
                    do x = NIA+1, NBasis
                          w_a(x) = THCData%fvw(x-NIA) - THCData%fij(a)
                    end do
                    ! b <= a: slab (b,a) gives the same sum as (a,b)
                    do b = NI0, a
                          do x = NIA+1, NBasis
                                w_b(x) = THCData%fvw(x-NIA) - THCData%fij(b)
                          end do

                          call real_aTb_x(V_FR, NV, RVI_FR(:,:,a), NChol, RVI_FR(:,:,b), NChol, NV, NV, NChol)
                          V_axby = V_FR

                          call real_aTb_x(V_axby, NV, RVI_FR1(:,:,a), NChol, RVI_FR(:,:,b), NChol, NV, NV, NChol, FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_FR2(:,:,a), NChol, RVI_FR(:,:,b), NChol, NV, NV, NChol, FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_FR(:,:,a), NChol, RVI_FR1(:,:,b), NChol, NV, NV, NChol, FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_FR(:,:,a), NChol, RVI_FR2(:,:,b), NChol, NV, NV, NChol, FRAC14, One)

                          call real_aTb_x(V_axby, NV, RVI_LR1(:,:,a), NCholErf, RVI_LR(:,:,b), NCholErf, &
                                NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_LR2(:,:,a), NCholErf, RVI_LR(:,:,b), NCholErf, &
                                NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_LR(:,:,a), NCholErf, RVI_LR1(:,:,b), NCholErf, &
                                NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(V_axby, NV, RVI_LR(:,:,a), NCholErf, RVI_LR2(:,:,b), NCholErf, &
                                NV, NV, NCholErf, -FRAC14, One)

                          call vovo_slab(E_mp2, w_a, w_b, a, b, NIA+1, NBasis, V_FR)
                          call vovo_slab(E_srmp2, w_a, w_b, a, b, NIA+1, NBasis, V_axby)
                    end do
              end do
              call tmsg('TIME FOR MP2-CBS[H] energy loop', timer, tdebug)

            end associate

      end subroutine srmp2_energy_thc


      ! ROMP2-CBS[H]: romp2_energy_thc with the dressed slab, per spin; every THC
      ! factor rotated with U_s. Singles stay FR. Needs cbs 2.
      ! See .secret_sauce/.skills/mp2/code-mp2.md, theory-mp2.md Sec. 5.
      subroutine sromp2_energy_thc(THCData, AuxData, CAONO, Flags, Xga_a, Xga_b, U_a, U_b, &
            e_a, e_b, fov_a, fov_b, E_aa, E_bb, E_ab, E_1, D_aa, D_bb, D_ab)

            type(TTHCData), intent(in) :: THCData
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:,:), intent(in) :: Xga_a, Xga_b, U_a, U_b
            double precision, dimension(:), intent(in) :: e_a, e_b
            double precision, dimension(:,:), allocatable, intent(in) :: fov_a, fov_b
            ! E_*: undressed ROMP2, D_*: dressed (CBS[H]) doubles
            double precision, intent(out) :: E_aa, E_bb, E_ab, E_1, D_aa, D_bb, D_ab

            double precision, dimension(:,:), allocatable :: XMuMat, TX, TXE, Y1, Y2, Y3
            double precision, dimension(:,:), allocatable :: V_FR, V
            double precision, dimension(:), allocatable :: w_i, w_j
            type(TDressedR) :: SA, SB
            double precision :: AvMU, ECorrMD
            integer :: NI, NA, NV, NIA, NBasis, NI0
            integer :: i, j
            type (tclock) :: timer

            if (.not. allocated(THCData%XgaErf)) then
                  print*, 'ERROR: ROMP2-CBS[H] needs the erf THC factors (cbs 2 in input)'
                  stop
            end if

            NI = AuxData%NI
            NA = AuxData%NA
            NV = AuxData%NV
            NIA = AuxData%NIA
            NBasis = AuxData%NBasis
            NI0 = AuxData%NCoreOrb + 1

            call loc_mu_thc(THCData, AuxData, CAONO, Flags, XMuMat, AvMU, ECorrMD)

            associate(Xga=>THCData%Xga, Zgk=>THCData%Zgk, NChol=>THCData%NChol, NTHC=>THCData%NTHC, &
                  XgaErf=>THCData%XgaErf, ZgkErf=>THCData%ZgkErf, &
                  NCholErf=>THCData%NCholErf, NTHCErf=>THCData%NTHCErf)

              call print_section('ROMP2-CBS[H] energy evaluation')
              call print_info('Frozen core orbitals', NI0 - 1)
              call print_info('alpha occ/vir', (/ NIA - NI0 + 1, NBasis - NIA /))
              call print_info('beta  occ/vir', (/ NI - NI0 + 1, NBasis - NI /))

              ! mu-weighted factors in the ROHF basis
              allocate(TX(NTHC, NBasis), TXE(NTHCErf, NBasis))
              call real_abT(TX, Xga, XMuMat)
              call real_abT(TXE, XgaErf, XMuMat)

              ! rotate to each spin and build its six vector sets
              call clock_start(timer)
              allocate(Y1(NTHC, NBasis), Y2(NTHCErf, NBasis), Y3(NTHCErf, NBasis))
              call real_ab(Y1, TX, U_a)
              call real_ab(Y2, XgaErf, U_a)
              call real_ab(Y3, TXE, U_a)
              call dressed_vectors(SA, Xga_a, Y1, Y2, Y3, NIA+1, NI0, NIA)
              call real_ab(Y1, TX, U_b)
              call real_ab(Y2, XgaErf, U_b)
              call real_ab(Y3, TXE, U_b)
              call dressed_vectors(SB, Xga_b, Y1, Y2, Y3, NI+1, NI0, NI)
              deallocate(Y1, Y2, Y3, TX, TXE)
              call CalcMem3(SA%FR, 'RA (x6)')
              call CalcMem3(SB%FR, 'RB (x6)')
              call tmsg('TIME FOR dressed RA, RB', timer, tdebug)

              call clock_start(timer)

              ! opposite spin: every alpha i with every beta j
              allocate(V_FR(NV, NA+NV), V(NV, NA+NV), w_i(NV), w_j(NA+NV))
              E_ab = zero
              D_ab = zero
              do i = NI0, NIA
                    w_i = e_a(NIA+1:NBasis) - e_a(i)
                    do j = NI0, NI
                          w_j = e_b(NI+1:NBasis) - e_b(j)
                          call dressed_slab(V_FR, V, SA, i, SB, j)
                          call romp2_os_slab(E_ab, w_i, w_j, V_FR)
                          call romp2_os_slab(D_ab, w_i, w_j, V)
                    end do
              end do
              deallocate(V_FR, V, w_i, w_j)

              call same_spin(E_aa, D_aa, SA, e_a, NIA+1, NI0, NIA)
              call same_spin(E_bb, D_bb, SB, e_b, NI+1, NI0, NI)

              ! singles, undressed
              E_1 = zero
              do i = NI0, NIA
                    E_1 = E_1 - sum(fov_a(i,:)**2 / (e_a(NIA+1:NBasis) - e_a(i)))
              end do
              do i = NI0, NI
                    E_1 = E_1 - sum(fov_b(i,:)**2 / (e_b(NI+1:NBasis) - e_b(i)))
              end do

              call tmsg('TIME FOR ROMP2-CBS[H] energy loop', timer, tdebug)

            end associate

      contains

            ! The six vector sets of one spin. X: Xga_s, TX: TXga_s, XE: XgaErf_s,
            ! TXE: TXgaErf_s; virtual v0..NBasis, occupied o0..o1.
            subroutine dressed_vectors(S, X, TX, XE, TXE, v0, o0, o1)
                  type(TDressedR), intent(out) :: S
                  double precision, dimension(:,:), intent(in) :: X, TX, XE, TXE
                  integer, intent(in) :: v0, o0, o1
                  integer :: nv, no

                  associate(Zgk=>THCData%Zgk, NChol=>THCData%NChol, NTHC=>THCData%NTHC, &
                        ZgkErf=>THCData%ZgkErf, NCholErf=>THCData%NCholErf, NTHCErf=>THCData%NTHCErf)

                    nv = NBasis - v0 + 1
                    no = o1 - o0 + 1
                    allocate(S%FR(NChol, nv, o0:o1), S%FR1(NChol, nv, o0:o1), S%FR2(NChol, nv, o0:o1))
                    allocate(S%LR(NCholErf, nv, o0:o1), S%LR1(NCholErf, nv, o0:o1), S%LR2(NCholErf, nv, o0:o1))
                    call thc_gammcor_Rkab_2(S%FR,  X(:,v0:NBasis),  X(:,o0:o1),  Zgk, nv, no, NChol, NTHC)
                    call thc_gammcor_Rkab_2(S%FR1, TX(:,v0:NBasis), X(:,o0:o1),  Zgk, nv, no, NChol, NTHC)
                    call thc_gammcor_Rkab_2(S%FR2, X(:,v0:NBasis),  TX(:,o0:o1), Zgk, nv, no, NChol, NTHC)
                    call thc_gammcor_Rkab_2(S%LR,  XE(:,v0:NBasis),  XE(:,o0:o1),  ZgkErf, nv, no, NCholErf, NTHCErf)
                    call thc_gammcor_Rkab_2(S%LR1, TXE(:,v0:NBasis), XE(:,o0:o1),  ZgkErf, nv, no, NCholErf, NTHCErf)
                    call thc_gammcor_Rkab_2(S%LR2, XE(:,v0:NBasis),  TXE(:,o0:o1), ZgkErf, nv, no, NCholErf, NTHCErf)

                  end associate
            end subroutine dressed_vectors

            ! V_FR(A,B) = (iA|jB) and V = its dressed form, i of the spin of P,
            ! j of the spin of Q (the 9 products of srmp2_energy_thc).
            subroutine dressed_slab(V_FR, V, P, i, Q, j)
                  double precision, dimension(:,:), intent(out) :: V_FR, V
                  type(TDressedR), intent(in) :: P, Q
                  integer, intent(in) :: i, j
                  integer :: m, n, k, ke

                  m = size(V, 1)
                  n = size(V, 2)
                  k = THCData%NChol
                  ke = THCData%NCholErf
                  call real_aTb_x(V_FR, m, P%FR(:,:,i), k, Q%FR(:,:,j), k, m, n, k)
                  V = V_FR
                  call real_aTb_x(V, m, P%FR1(:,:,i), k, Q%FR(:,:,j), k, m, n, k, FRAC14, One)
                  call real_aTb_x(V, m, P%FR2(:,:,i), k, Q%FR(:,:,j), k, m, n, k, FRAC14, One)
                  call real_aTb_x(V, m, P%FR(:,:,i), k, Q%FR1(:,:,j), k, m, n, k, FRAC14, One)
                  call real_aTb_x(V, m, P%FR(:,:,i), k, Q%FR2(:,:,j), k, m, n, k, FRAC14, One)
                  call real_aTb_x(V, m, P%LR1(:,:,i), ke, Q%LR(:,:,j), ke, m, n, ke, -FRAC14, One)
                  call real_aTb_x(V, m, P%LR2(:,:,i), ke, Q%LR(:,:,j), ke, m, n, ke, -FRAC14, One)
                  call real_aTb_x(V, m, P%LR(:,:,i), ke, Q%LR1(:,:,j), ke, m, n, ke, -FRAC14, One)
                  call real_aTb_x(V, m, P%LR(:,:,i), ke, Q%LR2(:,:,j), ke, m, n, ke, -FRAC14, One)
            end subroutine dressed_slab

            ! E_ss (undressed) and D_ss (dressed) for one spin, occupied o0..o1,
            ! virtual v0..NBasis, loop j < i as in romp2_energy_thc.
            subroutine same_spin(E_ss, D_ss, S, e, v0, o0, o1)
                  double precision, intent(out) :: E_ss, D_ss
                  type(TDressedR), intent(in) :: S
                  double precision, dimension(:), intent(in) :: e
                  integer, intent(in) :: v0, o0, o1
                  double precision, dimension(:,:), allocatable :: W_FR, W
                  double precision, dimension(:), allocatable :: w_i, w_j
                  integer :: nvir, p, q

                  nvir = NBasis - v0 + 1
                  allocate(W_FR(nvir, nvir), W(nvir, nvir), w_i(nvir), w_j(nvir))
                  E_ss = zero
                  D_ss = zero
                  do p = o0, o1
                        w_i = e(v0:NBasis) - e(p)
                        do q = o0, p - 1
                              w_j = e(v0:NBasis) - e(q)
                              call dressed_slab(W_FR, W, S, p, S, q)
                              call romp2_ss_slab(E_ss, w_i, w_j, W_FR)
                              call romp2_ss_slab(D_ss, w_i, w_j, W)
                        end do
                  end do
            end subroutine same_spin

      end subroutine sromp2_energy_thc

end module mp
