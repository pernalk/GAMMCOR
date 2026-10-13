module phac_dispatch

      use iso_fortran_env

      use display
      use types
      use clock

      use math_constants
      use ppac_types
      use phac0_thc
      use print_utils

      implicit none
      private

      public :: acph_driver
      public :: TACppData, FlagsData, TTHCData
contains

      subroutine acph_driver(THCData, AuxData, Flags, CAONO)
            type(TACppData), intent(inout) :: AuxData
            type(FlagsData), intent(in)    :: Flags
            double precision, dimension(:,:), intent(in) :: CAONO
            type(TTHCData), intent(inout) :: THCData
            type (tclock) :: timer
            integer :: ij, ind, i, j, k
            double precision :: E_corr
            double precision, dimension(:), allocatable :: DeltaH



            associate(Occ=>AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux)



              allocate(AuxData%map(NBasis))
              allocate(AuxData%IndX(NBasis**2))
              allocate(AuxData%IndN(2,Nbasis**2))
              allocate(AuxData%IPair(NBasis,NBasis))

              AuxData%IPair = 0


              AuxData%version = VER_PH_RPA_MULTI
              AuxData%general_version = GENVER_RPA_SING
              AuxData%ACType = PHAC

              call print_section('phAC0 orbital-pair selection')
              call print_info('Threshold for quasi-degeneracy', AuxData%ThrSelAct)
              call print_info('Threshold for quasi-virtual orbital', AuxData%ThrQVirt)
              call print_info('Threshold for quasi-inactive orbital', AuxData%ThrQInact)


            ind = 0
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i-1
                        ij = ij + 1

                        if (IAux(i)+IAux(j).ne.0.and.IAux(i)+IAux(j).ne.4) then
                              if ((IAux(i).eq.1).and.(IAux(j).eq.1).and. &
                                    (abs(Occ(i)-Occ(j))/Occ(i).lt.AuxData%ThrSelAct)) then

                                    write(LOUT,'(2X,"Discarding nearly degenerate pair ",2I4)') i, j
                              else
                                    if ((Flags%IFlCore.eq.1) .or. (Flags%IFlCore.eq.0.and.i.gt.NI.and.j.gt.NI)) then
                                          if ((abs(Occ(i) + Occ(j) - Two).gt.AuxData%ThrQInact) .and. &
                                                (abs(Occ(i) + Occ(j)).gt.AuxData%ThrQVirt)) then

                                                ind = ind + 1
                                                AuxData%IndX(ind) = ind
                                                AuxData%IndN(1, ind) = i
                                                AuxData%IndN(2, ind) = j

                                                AuxData%IPair(i,j) = 1
                                                AuxData%IPair(j,i) = 1
                                          end if
                                    end if

                              end if
                        end if

                  end do
            end do

            AuxData%NDim = ind

            AuxData%map = 0



            k = 1
            do i = 1, NI+NA+NV
                  if(AuxData%IndAux(i)==1)then
                        AuxData%map(i) = k
                        k = k+1
                  end if
                  if(AuxData%IndAux(i)==2)then
                        AuxData%map(i) = 0
                  end if
            end do



            call print_info("Original number of pq pairs", ij)
            call print_info("Reduced number of pq pairs", ind)

            write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'Occ(i)', 'Occ(j)'
            write(*, '(A43)') repeat('-', 43)

            do i = 1, min(100, AuxData%NDim)
                  write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                        Occ(AuxData%IndN(1,i)), Occ(AuxData%IndN(2,i))
            end do


            select case(Flags%Jobtype)

            case(JOB_TYPE_AC0, JOB_TYPE_MSAC0)

                  select case(Flags%ITwoEl)

                  case(TWOMO_INCORE)
                        if (Flags%Algorithm == ALG_THC) then
                              ! THCData%TwoEl is set up in dispatch_acxx, the
                              ! canonicalization in read_PYSCF (thc_exact_ints)
                              call print_section('phAC0 calculation')
                              call print_info('Two-electron integrals', 'exact (PySCF), THC code')
                              call clock_start(timer)
                              call ACPH0_THC(THCData, AuxData, CAONO, Flags, E_corr, DeltaH)
                              call tmsg('TIME FOR ACPH0_THC', timer, 0)
                        end if

                  case(3)
                        select case(Flags%ICholeskyTHC)

                        case(1)
                              call print_section('phAC0 calculation')
                              call print_info('Two-electron integrals', 'THC')
                              call clock_start(timer)
                              call ACPH0_THC(THCData, AuxData, CAONO, Flags, E_corr, DeltaH)
                              call tmsg('TIME FOR ACPH0_THC', timer, 0)

                        case default
                              !
                              print*, 'not implemented'
                              stop
                              !
                        end select

                  case default
                        !
                        print*, 'not implemented'
                        stop
                        !
                  end select

            case default
                  !
                  print*, 'not implemented'
                  stop
                  !
            end select
          end associate
      end subroutine acph_driver


end module phac_dispatch
