module mp

      ! =============================================================================
      ! MP2 and SR-MP2 (MP2-CBS[H]) with THC integrals, JobType MP2 / SRMP2.
      ! Reference: the HF dump of the PySCF interface (get_rohf_for_gammcor,
      ! load_mp2_aux_data), closed shell only (NA = 0).
      !
      ! THEORY. With NA = 0 the ph-AC0 energy has only one term, (31)(31) =
      ! E_vovo, and that is closed-shell MP2 (thc_energy_ph, vovo_contr). Its vo
      ! ERPA eigenvalues are the orbital energy differences fvw(v) - fij(i)
      ! (PHERPA0_THC_VO), and canonicalize already stores fij, fvw in THCData.
      ! So no ERPA is solved: mp2_energy_thc is the (IV|IV) branch of
      ! THC_energy_loop_ph with these omegas and the same kernel, vovo_slab.
      !
      ! FROZEN CORE. Inactive orbitals 1..NCoreOrb (input) stay in the
      ! reference (fij, Fock) and are left out of the pairs: a, b run NI0..NI,
      ! NI0 = NCoreOrb + 1. AuxData%NFrozen is set only on the ph path.
      !
      ! SRMP2 is still the old code (srmp2_energy_old), to be redone on top of
      ! mp2_energy_thc: .implementation_plans/mp2.md, steps 6-8.
      ! =============================================================================

      use ppac_types
      use types
      use math_constants
      use real_linalg !from gammcor integrals
      use THC_Gammcor
      use sort  !from gammcor integrals
      use clock
      use print_utils
      use thc_energy_ph, only: vovo_slab

      implicit none

      integer, parameter, private :: tdebug = 0
      double precision, parameter, private :: ttoeV = 27.211399d+0

contains

      subroutine mp2_driver(THCData, AuxData, CAONO, Flags)

            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData

            double precision :: E_mp2

            if (.not. allocated(THCData%Xga)) then
                  print*, 'ERROR: MP2/SRMP2 needs TwoMoInt FOFO with Cholesky THC'
                  stop
            end if

            select case (Flags%JobType)
            case (JOB_TYPE_MP2)
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

            case (JOB_TYPE_SRMP2)
                  call srmp2_energy_old(THCData, AuxData, CAONO, Flags)
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
            double precision, dimension(:), allocatable :: om_a, om_b
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
              allocate(om_a(NIA+1:NBasis), om_b(NIA+1:NBasis))

              call clock_start(timer)
              E_mp2 = zero
              do a = NI0, NI
                    do x = NIA+1, NBasis
                          om_a(x) = THCData%fvw(x-NIA) - THCData%fij(a)
                    end do
                    ! b <= a: slab (b,a) gives the same sum as (a,b)
                    do b = NI0, a
                          do x = NIA+1, NBasis
                                om_b(x) = THCData%fvw(x-NIA) - THCData%fij(b)
                          end do
                          call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                          call vovo_slab(E_mp2, om_a, om_b, a, b, NIA+1, NBasis, V_axby)
                    end do
              end do
              call tmsg('TIME FOR MP2 energy loop', timer, tdebug)

            end associate

      end subroutine mp2_energy_thc


      ! =============================================================================
      ! SRMP2: the former mp2_driver, unchanged except that the never-used RII_*
      ! and RAV arrays are gone. XMuMat is not filled (the mu(r) call is
      ! commented out), so TXga/TXgaErf and the result are not meaningful yet.
      ! To be replaced by the dressed slab + vovo_slab, plan step 8.
      ! =============================================================================
      subroutine srmp2_energy_old(THCData, AuxData, CAONO, Flags)

            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData

            double precision, dimension(:,:),  allocatable  :: XMuMat
            double precision, dimension(:,:),  allocatable  :: URe, CAONOT
            double precision :: AvMu, ECorrMD, emp2
            double precision, dimension(:,:,:), allocatable :: RVI_FR, RVI_FR1, RVI_FR2
            double precision, dimension(:,:,:), allocatable :: RVI_LR, RVI_LR1, RVI_LR2
            double precision, dimension(:,:), allocatable :: this_VV, this_FR, this_SR
            double precision :: numerator, denominator
            double precision, dimension(:), allocatable :: occ_frozen
            integer :: a, b, i, j, i0

            associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga,  TXga=>THCData%TXga, &
                  NChol=>THCData%NChol, NTHC=>THCData%NTHC, &
                  XgaErf=>THCData%XgaErf, TXgaErf=>THCData%TXgaErf, &
                  ZgkErf=>THCData%ZgkErf, NCholErf=>THCData%NCholErf, NTHCErf=>THCData%NTHCErf, NBasis=>AuxData%NBasis, &
                  NI=>AuxData%NI, NIA=>AuxData%NIA, NV=>AuxData%NV)

              allocate(occ_frozen(NBasis))
              occ_frozen = AuxData%occ

              allocate(RVI_FR(NChol, NV, NI))
              allocate(RVI_FR1(NChol, NV, NI))
              allocate(RVI_FR2(NChol, NV, NI))
              allocate(RVI_LR(NCholErf, NV, NI))
              allocate(RVI_LR1(NCholErf, NV, NI))
              allocate(RVI_LR2(NCholErf, NV, NI))

              ! For frozen core
              if (AuxData%NCoreOrb > 0)then
                    i0 = 1+ AuxData%NCoreOrb
                    do i = 1, i0-1
                          occ_frozen(i) = zero
                    end do
              else
                    i0 = 1
              end if

              allocate(URe(NBasis, NBasis))
              allocate(CAONOT(NBasis, NBasis))
              URe = zero
              do i = 1, NBasis
                    URe(i, i) = one
              end do
              CAONOT = transpose(CAONO)
              allocate(XMuMat(nbasis, nbasis))

              !Call LOC_MU_CBS_CHOL(Flags,XMuMat,ECorrMD,AvMU,URe, CAONOT,Occ_frozen,NBasis)

              call real_abt(THCData%TXgaErf, THCData%XgaErf, XMuMAT)
              call real_abt(THCData%TXga, THCData%Xga, XMuMAT)

              call thc_gammcor_Rkab_2(RVI_FR, Xga(:,NIA+1:NV), Xga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_FR1, TXga(:,NIA+1:NV), Xga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_FR2, Xga(:,NIA+1:NV), TXga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
              call thc_gammcor_Rkab_2(RVI_LR, XgaErf(:,NIA+1:NV), XgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)
              call thc_gammcor_Rkab_2(RVI_LR1, TXgaErf(:,NIA+1:NV), XgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)
              call thc_gammcor_Rkab_2(RVI_LR2, XgaErf(:,NIA+1:NV), TXgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)

              allocate(this_VV(NV, NV))
              allocate(this_SR(NV, NV))
              allocate(this_FR(NV, NV))

              emp2 = zero
              do i = i0, NI
                    do j = i0, NI

                          call real_aTb_x(this_FR, NV, RVI_FR(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, One)

                          call real_aTb_x(this_SR, NV, RVI_FR1(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, FRAC14)
                          call real_aTb_x(this_SR, NV, RVI_FR2(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, FRAC14, One)
                          call real_aTb_x(this_SR, NV, RVI_FR(:,:,i), NChol, RVI_FR1(:,:,j), NChol, NV, NV, NChol, FRAC14, One)
                          call real_aTb_x(this_SR, NV, RVI_FR(:,:,i), NChol, RVI_FR2(:,:,j), NChol, NV, NV, NChol, FRAC14, One)

                          call real_aTb_x(this_SR, NV, RVI_LR1(:,:,i), NCholErf, RVI_LR(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(this_SR, NV, RVI_LR2(:,:,i), NCholErf, RVI_LR(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(this_SR, NV, RVI_LR(:,:,i), NCholErf, RVI_LR1(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                          call real_aTb_x(this_SR, NV, RVI_LR(:,:,i), NCholErf, RVI_LR2(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)

                          this_VV = this_FR + this_SR

                          do a = 1, NV
                                do b = 1, NV
                                      numerator  =  this_VV(a,b) *(two * this_VV(a,b) - this_VV(b,a))
                                      denominator = THCData%fij(i)+THCData%fij(j)-THCData%fvw(a)-THCData%fvw(b)
                                      emp2  = emp2 + numerator/denominator
                                end do
                          end do
                    end do
              end do
              print*, 'SRMp2', emp2

            end associate
      end subroutine srmp2_energy_old

end module mp
