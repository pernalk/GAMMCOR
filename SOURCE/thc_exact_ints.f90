module thc_exact_ints

      ! =============================================================================
      ! Exact-integral route of the ph-AC0 THC code (Algorithm THC + TwoMoInt INCORE):
      ! each THC slab is filled from the packed exact integrals TwoEl instead of from
      ! R blocks (exact_V_axby, exact_V_pair, exact_V_xay), and canonicalize_exact
      ! replaces canonicalize. Not yet tested.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-exact-ints.md.
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
      ! Exact counterpart of canonicalize, called by read_PYSCF for every INCORE
      ! AC0/MSAC0 run: diagonalizes the inactive and virtual blocks of the Fock
      ! matrix of the computed state (not the SA one PySCF used) and rotates HNO0,
      ! TwoEl and CAONO; sets fij, fvw, ECAS_THC.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-exact-ints.md.
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
