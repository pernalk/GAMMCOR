module dispatcher
 
      use types
      use ppac_types
      use ppac_simple
      use mp
      use phac_spinres, only: acph_spinres_driver
      use ppac_dispatch, only: acpp_driver
      use phac_dispatch, only: acph_driver


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
                  print*, 'ppac_dispatch_simple'
                  call ppac_dispatch_simple(THCData, AuxData, CAONO, Flags, TwoEl)
                  
            case(JOB_TYPE_PPERPA, JOB_TYPE_AC0PP, JOB_TYPE_ACPP)
                  print*, 'Flags%Algorithm', Flags%Algorithm
                  if (Flags%Algorithm == ALG_PPSIMPLE)then
                        print*, 'ppac_dispatch_simple'
                        call ppac_dispatch_simple(THCData, AuxData, CAONO, Flags, TwoEl)
                  else
                        print*, 'acpp_driver'
                        call acpp_driver(THCData, AuxData, Flags, CAONO)
                  end if

            case(JOB_TYPE_MP2)
                  ! Flags%IFlCore = 0 is set in systemdef (fill_Flags)
                  print*, 'mp2_driver'
                  call mp2_driver(THCData, AuxData, CAONO, Flags)

            case(JOB_TYPE_AC0, JOB_TYPE_MSAC0)
                  if (Flags%ITwoEl == TWOMO_FOFO .and. Flags%ICholeskyTHC == 1) then
                        print*, 'acph_driver'
                        call acph_driver(THCData, AuxData, Flags, CAONO)
                  else if (Flags%ITwoEl == TWOMO_INCORE .and. Flags%Algorithm == ALG_THC) then
                        ! the THC code on the exact PySCF integrals: TwoEl (NO basis,
                        ! already canonicalized by read_PYSCF, canonicalize_exact)
                        ! travels in THCData
                        print*, 'acph_driver (exact integrals)'
                        call move_alloc(TwoEl, THCData%TwoEl)
                        THCData%exact_ints = .true.
                        call acph_driver(THCData, AuxData, Flags, CAONO)
                        call move_alloc(THCData%TwoEl, TwoEl)
                  else if (Flags%Algorithm == ALG_SPINRES) then
                        print*, 'acph_spinres_driver'
                        call acph_spinres_driver(THCData, AuxData, Flags, IntsO)
                  end if

            case(JOB_TYPE_DUCC)
                  print*, 'acph_spinres_driver'
                  call acph_spinres_driver(THCData, AuxData, Flags, IntsO, IntsD)

            case default

                  if (Flags%Algorithm == ALG_SPINRES) then
                        print*, 'acph_spinres_driver'
                        call acph_spinres_driver(THCData, AuxData, Flags, IntsO)
                  end if

            end select

      end subroutine dispatch_acxx

end module dispatcher
