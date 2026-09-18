module tddft

      use types
      use ppac_types
      use tddft_types
!     use mp_simple

      implicit none

contains

      subroutine dispatch_tddft(THCData, TDDFTData, Flags, TwoEl)
            type(TTHCData), intent(inout) :: THCData
            type(TtddftData), intent(inout) :: TDDFTData
            type(FlagsData), intent(inout) :: Flags
            double precision, dimension(:), allocatable, intent(inout) :: TwoEl

            integer :: i


            select case (Flags%JobType)
            case (JOB_TYPE_MP2)
                  if (Flags%RDMType == RDM_TYPE_RKS) then
                        print*, 'tutaj trzeba wywolac funckje do prostego mp2 z modulu mp simple'
                        print*, ''
                        print*, 'energie orbitalne test'
                        do i = 1, TDDFTData%NBasis
                              print*, i, TDDFTData%eorbs(i)
                        end do
                        
                        stop
                  end if

            case (JOB_TYPE_TDDFT)
                  ! call run_tddft(Args....) albo inna nazwa
                  print*,  'Not yet implemented, exiting'
                  stop
            end select

      end subroutine dispatch_tddft


!      subroutine run_tddft(Args....)
            !
            !
            !
!      end subroutine run_tddft

end module tddft
