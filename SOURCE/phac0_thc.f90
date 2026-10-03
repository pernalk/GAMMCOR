module phac0_thc


      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals
      use types
      use lin
      use clock

      use math_constants
      use pp_utils
      use print_utils
      use ppac_types
      use thc_energy
      use tran
      use omp_lib


      implicit none

contains


      subroutine ACPH0_THC(THCData, AuxData, CAONO_IN, Flags)
            use THC_Gammcor

            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO_IN
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData


            if (AuxData%PYSCF==1)then
                  call print_info('External ordering', 'PySCF')
            end if

            associate(NBasis=>AuxData%NBasis, map=>AuxData%map, Occ=>AuxData%Occ, NI=>AuxData%NI, NV=>AuxData%NV)
              ACAlpha = one !because we are simultaneously calculating A1.

              allocate(HNO0(NBasis, NBasis))
              allocate(HNO(NBasis, NBasis))
              allocate(CAONO(NBasis, NBasis))
              CAONO = CAONO_IN
              allocate(HNO0_THC(Nbasis, NBasis))
              allocate(posS(NBasis, NBasis))
              allocate(posT(NBasis, NBasis))
              posS = 0
              posT = 0
              
              call init_AC0Blocks(AuxData, ACBlocks, posS, posT)
              call init_rdm(RdmData, AuxData)
              allocate(Aux3A(NBasis, NBasis))
              allocate(AuxII(NBasis, NBasis))
              allocate(AuxAA(NBasis, NBasis))
              allocate(BuxII(NBasis, NBasis))
              allocate(BuxAA(NBasis, NBasis))

              AuxII = zero
              BuxII = zero
              AuxAA = zero
              BuxAA = zero
              Aux3A = zero

              !call THC_int_loop_ph(THCData, ACBlocks, AuxII, AuxAA, BuxII, BuxAA, Aux3X, Aux3B, &
              !Flags, AuxData, Occ, map, RdmData, posS, posT, AcAlpha)
                    !call update_HNO_AC0(AuxData, HNO, THCData%fij, THCData%fvw, AuxII, BuxII)

              !call PPERPA0_THC_OOVV
              !call eigs(AuxData, ACBlocks(i), i)
              !call PPERPA0_THC_AA
              !call eigs(AuxData, ACBlocks(A0aa), A0aa)
              !call PPERPA0_THC_VA_Eigs(AuxData, ACBlocks(A0va)
              !call PPERPA0_THC_OA_Eigs(AuxData, ACBlocks(A0oa)
              !do i = 9, 13
              ! call PPERPA_THC
              ! call THC_energy_loop
              
            end associate
            
            call msg('ACPH0_THC: not yet implemented')            
            stop

      end subroutine ACPH0_THC


end module phac0_thc
