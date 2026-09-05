module dispatcher
 
      use types
      use ppac_types
      use ppac_simple
      use mp
      use phac_spinres
      use ppac_dispatch, only: acpp_driver


      implicit none


contains

      subroutine dispatch_acxx(THCData, AuxData, CAONO, Flags, TwoEl, IntsO, IntsD)

            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), allocatable, intent(inout) :: CAONO
            type(FlagsData), intent(inout) :: Flags
            double precision, dimension(:), allocatable, intent(inout) :: TwoEl
            type(TInts), intent(inout) :: IntsO
            type(TInts), intent(inout) :: IntsD

            if (.not.allocated(CAONO)) allocate(CAONO(1,1))
            if (.not.allocated(TwoEl)) allocate(TwoEl(1))
            if (.not.allocated(IntsO%ints2e)) allocate(IntsO%ints2e(1))
            if (.not.allocated(IntsD%ints2e)) allocate(IntsD%ints2e(1))

            select case(Flags%JobType)

            case(JOB_TYPE_PPERPA_RDMDUMP, JOB_TYPE_HHERPA_RDMDUMP)
                  print*, 'acpp_driver_simple'
                  call acpp_driver_simple(THCData, AuxData, CAONO, Flags, TwoEl)
                  
            case(JOB_TYPE_PPERPA, JOB_TYPE_AC0PP, JOB_TYPE_ACPP)

                  if (Flags%Algorithm == ALG_PPSIMPLE)then
                        print*, 'acpp_driver_simple'
                        call acpp_driver_simple(THCData, AuxData, CAONO, Flags, TwoEl)
                  else
                        print*, 'acpp_driver'
                        call acpp_driver(THCData, AuxData, Flags, CAONO)
                  end if

            case(JOB_TYPE_MP2, JOB_TYPE_SRMP2)

                  Flags%IFlCore = 0
                  print*, 'mp2_driver'
                  call mp2_driver(THCData, AuxData, CAONO, Flags)

            case(JOB_TYPE_DUCC)
                  print*, 'acph_driver'
                  call acph_driver(THCData, AuxData, Flags, IntsO, IntsD)

            case default

                  if (Flags%Algorithm == ALG_SPINRES) then
                        print*, 'acph_driver'
                        call acph_driver(THCData, AuxData, Flags, IntsO)
                  end if

            end select

      end subroutine dispatch_acxx

end module dispatcher
