module thc_exact_ints

      ! =============================================================================
      ! Exact-integral debug route of the ph-AC0 THC code (phac0_thc, thc_energy_ph):
      ! Algorithm THC + TwoMoInt INCORE + JobType AC0/MSAC0, PySCF interface.
      ! Plan and tests: implementation_plans/thc_with_exact.md.
      !
      ! The THC code builds every integral slab it uses with one gemm over R
      ! blocks. With THCData%exact_ints set, each of those gemms is replaced by a
      ! fill of the same slab, in the same layout, from the exact PySCF integrals
      ! THCData%TwoEl (packed, TwoEl(NAddr3(p,q,r,s)) = (pq|rs), NO basis):
      !
      !   exact_V_axby  V(x,y)   = (x a|y b)    for real_aTb_x(V, ., Ra(:,:,a), ., Rb(:,:,b), ...)
      !   exact_V_pair  V(x,y)   = (x y|c d)    for real_aTv_x(V, RAA, ., R(:,c,d), ...)
      !   exact_V_xay   V(x,a,y) = (x a|y b)    for V_xay of the (VA|IV) class
      !
      ! Nothing else in the THC code changes, so the run gives the THC code with
      ! exact integrals: any difference from INCORE AC0 is a code error, any
      ! difference from the THC run is the THC error. The R blocks are not built.
      !
      ! The THC-only quantities set by canonicalize in THC_init2 - fij, fvw, the
      ! inactive/virtual orbital rotation and ECAS_THC - come from
      ! canonicalize_exact, which read_PYSCF calls for every INCORE AC0/MSAC0 run
      ! (the PySCF-wrapper INCORE run included).
      ! =============================================================================

      use types
      use ppac_types
      use real_linalg
      use math_constants
      use print_utils

      implicit none

contains

      ! V(x,y) = alpha (x a|y b), x = x0..x1, y = y0..y1, a and b fixed (global
      ! indices) - the slab real_aTb_x builds from Ra(:,:,a) and Rb(:,:,b) when
      ! the first index of Ra runs over x0..x1 and that of Rb over y0..y1.
      subroutine exact_V_axby(V, TwoEl, a, x0, x1, b, y0, y1, alpha)
            integer, intent(in) :: a, x0, x1, b, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(out) :: V
            double precision, dimension(:), intent(in) :: TwoEl
            double precision, intent(in) :: alpha
            integer :: x, y
            integer(8), external :: NAddr3

            do y = y0, y1
                  do x = x0, x1
                        V(x, y) = alpha * TwoEl(NAddr3(x, a, y, b))
                  end do
            end do

      end subroutine exact_V_axby


      ! V(x,y) = alpha (x y|c d), x = x0..x1, y = y0..y1, the pair (c,d) fixed -
      ! the slab real_aTv_x builds from RAA(:,:,:) and one vector R(:,c,d) in the
      ! "together" loops of THC_int_loop_ph.
      subroutine exact_V_pair(V, TwoEl, x0, x1, y0, y1, c, d, alpha)
            integer, intent(in) :: x0, x1, y0, y1, c, d
            double precision, dimension(x0:x1, y0:y1), intent(out) :: V
            double precision, dimension(:), intent(in) :: TwoEl
            double precision, intent(in) :: alpha
            integer :: x, y
            integer(8), external :: NAddr3

            do y = y0, y1
                  do x = x0, x1
                        V(x, y) = alpha * TwoEl(NAddr3(x, y, c, d))
                  end do
            end do

      end subroutine exact_V_pair


      ! V(x,a,y) = alpha (x a|y b), all a = a0..a1 at once, b fixed - the V_xay
      ! that THC_energy_loop_ph builds from RVA and RIV(:,:,b) for E_vavo.
      subroutine exact_V_xay(V, TwoEl, x0, x1, a0, a1, y0, y1, b, alpha)
            integer, intent(in) :: x0, x1, a0, a1, y0, y1, b
            double precision, dimension(x0:x1, a0:a1, y0:y1), intent(out) :: V
            double precision, dimension(:), intent(in) :: TwoEl
            double precision, intent(in) :: alpha
            integer :: x, a, y
            integer(8), external :: NAddr3

            do y = y0, y1
                  do a = a0, a1
                        do x = x0, x1
                              V(x, a, y) = alpha * TwoEl(NAddr3(x, a, y, b))
                        end do
                  end do
            end do

      end subroutine exact_V_xay


      ! =============================================================================
      ! The exact counterpart of canonicalize (interface_pp.f90), on the INCORE
      ! integrals TwoEl (packed, NO basis). Called by read_PYSCF for every INCORE
      ! AC0/MSAC0 run, so the PySCF-wrapper INCORE run (AC0CAS) and the
      ! exact-integral THC run use the same orbitals.
      !
      !   F(p,q) = HNO0(p,q) + sum_{t <= NIA} Occ(t) [ 2 (pq|tt) - (pt|tq) ]
      !
      ! is the Fock matrix of the computed state (the one canonicalize builds with
      ! THC integrals; Occ(t) is enough because the orbitals are NOs). PySCF
      ! canonicalizes with the Fock matrix of its fcisolver density, which for
      ! SA-CASSCF is the state-averaged one, so the inactive and virtual blocks of
      ! F are diagonalized here for every reference:
      !
      !   fij, U_ii  from F(inact, inact)        fvw, U_vv  from F(virt, virt)
      !   U = diag(U_ii, 1_act, U_vv)
      !   HNO0 <- U^T HNO0 U,   TwoEl <- TwoEl rotated by U (TwoNO1),   CAONO <- CAONO U
      !
      ! The RDMs and the MS-AC0 transition 1-RDMs are active only and do not
      ! change. ECAS_THC is the CAS energy from these integrals, by the same
      ! formula as canonicalize. On exit THCData%fij/fvw hold the orbital
      ! energies and THCData%HNO = AuxData%HNO0.
      ! =============================================================================
      subroutine canonicalize_exact(THCData, AuxData, CAONO, TwoEl)
            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(inout) :: CAONO
            double precision, dimension(:), intent(inout) :: TwoEl

            double precision, dimension(:,:), allocatable :: F, U, Ublk, work, workC
            double precision :: ECAS, ECumul, val
            integer :: i, t, u_, v, w, N
            integer(8), external :: NAddr3

            associate(NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  Occ=>AuxData%Occ)

              N = AuxData%NBasis

              call print_section('Exact Fock matrix, canonicalization')

              allocate(F(N, N))
              call exact_fock(F, TwoEl, AuxData%HNO0, Occ, NIA, N)
              call print_info('max |F_ij|, i /= j inactive (before)', max_offdiag(F, 1, NI))
              call print_info('max |F_vw|, v /= w virtual  (before)', max_offdiag(F, NIA+1, N))

              ! CAS energy, as in canonicalize: Fock part + active cumulant
              ECAS = AuxData%ENuc
              do i = 1, NIA
                    ECAS = ECAS + Occ(i) * (AuxData%HNO0(i, i) + F(i, i))
              end do
              ECumul = zero
              do w = 1, NA
                    do v = 1, NA
                          do u_ = 1, NA
                                do t = 1, NA
                                      val = frac12 * AuxData%rdm2_full(t, u_, v, w)
                                      if (t == v .and. u_ == w) val = val - two * Occ(NI+t) * Occ(NI+u_)
                                      if (t == w .and. u_ == v) val = val + Occ(NI+t) * Occ(NI+u_)
                                      ECumul = ECumul + val * TwoEl(NAddr3(NI+t, NI+v, NI+u_, NI+w))
                                end do
                          end do
                    end do
              end do
              AuxData%ECAS_THC = ECAS + ECumul
              call print_energy1('CASSCF energy, exact integrals', AuxData%ECAS_THC)
              call print_energy1('CASSCF energy, PySCF', AuxData%ECAS_read)

              ! U = diag(U_ii, 1, U_vv); eigenvalues are the orbital energies
              allocate(U(N, N))
              U = zero
              do i = 1, N
                    U(i, i) = one
              end do
              if (NI > 0) then
                    allocate(Ublk(NI, NI))
                    Ublk = F(1:NI, 1:NI)
                    call symmetric_eigenproblem(THCData%fij, Ublk, NI, .true.)
                    U(1:NI, 1:NI) = Ublk
                    deallocate(Ublk)
              end if
              if (NV > 0) then
                    allocate(Ublk(NV, NV))
                    Ublk = F(NIA+1:N, NIA+1:N)
                    call symmetric_eigenproblem(THCData%fvw, Ublk, NV, .true.)
                    U(NIA+1:N, NIA+1:N) = Ublk
                    deallocate(Ublk)
              end if

              ! HNO0 <- U^T HNO0 U
              allocate(work(N, N))
              call real_ab(work, AuxData%HNO0, U)
              call real_aTb(AuxData%HNO0, U, work)
              THCData%HNO = AuxData%HNO0

              ! CAONO <- CAONO U (AO x NO)
              if (size(CAONO, 2) == N) then
                    allocate(workC(size(CAONO, 1), N))
                    workC = CAONO
                    call real_ab(CAONO, workC, U)
                    deallocate(workC)
              end if

              ! TwoEl: new p = sum_p' U(p',p) old p'. TwoNO1 takes URe(new,old),
              ! so it gets U^T (as Trans2NO_AO passes CAONO^T)
              work = transpose(U)
              call TwoNO1(TwoEl, work, N, AuxData%NInte2)

              ! check: the rotated data give a diagonal Fock matrix
              call exact_fock(F, TwoEl, AuxData%HNO0, Occ, NIA, N)
              call print_info('max |F_ij|, i /= j inactive (after)', max_offdiag(F, 1, NI))
              call print_info('max |F_vw|, v /= w virtual  (after)', max_offdiag(F, NIA+1, N))

              deallocate(F, U, work)

            end associate

      end subroutine canonicalize_exact


      ! F(p,q) = H0(p,q) + sum_{t <= NIA} Occ(t) [ 2 (pq|tt) - (pt|tq) ], all p, q.
      subroutine exact_fock(F, TwoEl, H0, Occ, NIA, N)
            double precision, dimension(:,:), intent(out) :: F
            double precision, dimension(:), intent(in) :: TwoEl
            double precision, dimension(:,:), intent(in) :: H0
            double precision, dimension(:), intent(in) :: Occ
            integer, intent(in) :: NIA, N
            integer :: t, p, q
            integer(8), external :: NAddr3

            F = H0
            do q = 1, N
                  do p = 1, N
                        do t = 1, NIA
                              F(p, q) = F(p, q) + Occ(t) * (two * TwoEl(NAddr3(p, q, t, t)) &
                                    - TwoEl(NAddr3(p, t, t, q)))
                        end do
                  end do
            end do

      end subroutine exact_fock


      ! Largest |F(p,q)|, p /= q, within the diagonal block i0..i1 (0 if empty).
      function max_offdiag(F, i0, i1) result(m)
            double precision, dimension(:,:), intent(in) :: F
            integer, intent(in) :: i0, i1
            double precision :: m
            integer :: p, q

            m = zero
            do q = i0, i1
                  do p = i0, i1
                        if (p /= q) m = max(m, abs(F(p, q)))
                  end do
            end do

      end function max_offdiag

end module thc_exact_ints
