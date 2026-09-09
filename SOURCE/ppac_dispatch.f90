module ppac_dispatch

      use iso_fortran_env

      use display
      use types
      use clock

      use math_constants
      use ppac_types
      use ppac0
      use ppac_subs
      use ppac0_thc
      use print_utils

      implicit none
      private

      public :: acpp_driver
      public :: TACppData, FlagsData, TTHCData
contains

      subroutine acpp_driver(THCData, AuxData, Flags, CAONO)
            type(TACppData), intent(inout) :: AuxData
            type(FlagsData), intent(in)    :: Flags
            double precision, dimension(:,:), intent(in) :: CAONO
            type(TTHCData), intent(inout) :: THCData
            double precision :: ECorr_AC, ECorr_AC0
            double precision :: ECorr
            double precision :: ETot
            type (tclock) :: timer, timer0
            integer :: ij, ind, indt, i, j, k 

                  
                  
            associate(Occ=>AuxData%Occ, XOne=>AuxData%XOne, ENuc=>AuxData%ENuc, &
                  TwoNO=>AuxData%TwoNO, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis)



              allocate(AuxData%map(NBasis))
              allocate(AuxData%IndMod(NBasis))
              allocate(AuxData%IndX_s(NBasis**2))
              allocate(AuxData%IndN_s(2,Nbasis**2))
              
              allocate(AuxData%IndX_t(NBasis**2))
              allocate(AuxData%IndN_t(2,Nbasis**2))
              
              allocate(AuxData%IndN(2,Nbasis**2))
              
              if (AuxData%Dalton == 1)then
                    call fill_IndAuxMod(AuxData%Occ, AuxData%IndAux, AuxData%IndMod, NBasis, NI, NA, NV)
              end if


              AuxData%version = VER_PP_RPA_MULTI
              AuxData%general_version = GENVER_PP
              AuxData%ACType = PPAC

              call print_section('ppAC0 orbital-pair selection')
              call print_info('ACPP Threshold for orbital pairs', AuxData%ThrPP)
            

            ind = 0
            indt = 0
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1

                        if (AuxData%version==VER_PP_RPA_MULTI)then
                              if ((abs((Occ(i)+Occ(j)-one)) > AuxData%ThrPP))then
                                    ind = ind + 1
                                    AuxData%IndN(1, ind) = i
                                    AuxData%IndN(2, ind) = j

                                    AuxData%IndN_s(1, ind) = i
                                    AuxData%IndN_s(2, ind) = j


                                    AuxData%IndX_s(ind) = ind
                                    if (j.ne.i)then
                                          indt = indt + 1
                                          AuxData%IndX_t(indt) = indt

                                          AuxData%IndN_t(1, indt) = i
                                          AuxData%IndN_t(2, indt) = j
                                     end if

                              end if
                        else !( pp_rpa_sing_determinant)                                                                                                                                               
                              ind = ind + 1
                              AuxData%IndN(1, ind) = i
                              AuxData%IndN(2, ind) = j

                              AuxData%IndN_s(1, ind) = i
                              AuxData%IndN_s(2, ind) = j


                              AuxData%IndX_s(ind) = ind

                              if (j.ne.i)then
                                    indt = indt +	1
                                    AuxData%IndX_t(indt) = indt

                                    AuxData%IndN_t(1, indt) = i
                                    AuxData%IndN_t(2, indt) = j
                              end if

                        end if
                  end do
            end do

            AuxData%NDim = ind
            AuxData%NDim_s = ind
            AuxData%NDim_t = indt

            AuxData%map = zero



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
            call print_info("Reduced singlet pq pairs", ind)
            call print_info("Reduced triplet pq pairs", indt)

            write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'Occ(i)', 'Occ(j)'
            write(*, '(A43)') repeat('-', 43) 

            do i = 1, min(100, AuxData%NDim_s)
                  write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                        Occ(AuxData%IndN(1,i)), Occ(AuxData%IndN(2,i))
            end do


            ECorr = zero
            ECorr_AC0 = zero
            ECorr_AC = zero


            select case(Flags%Jobtype)

            case(JOB_TYPE_AC0PP)

                  select case(Flags%ITwoEl)

                  case(1)
                        call msg("Starting ppAC0 with integrals in RAM")

                        call clock_start(timer0)
                        call ACPP0_fast(AuxData%HType, ECorr_AC0, & 
                              ETot, ENuc,Occ, XOne, TwoNO, AuxData, AuxData%IndN_s, AuxData%IndAux, AuxData%IndMod, &
                              NBasis, NA, NI, NV, NInte1, NInte2, Flags, AuxData%ACType)
                        print*, 'RDSC ACPP0 ECORR', ETot+ECorr_AC0+ENuc
                        print*, 'Czas na fast ac0 version: ', clock_readwall(timer0)

                  case(3)
                        select case(Flags%ICholeskyTHC)

                        case(1)                              
                              call print_section('ppAC0 calculation')
                              call print_info('Two-electron integrals', 'THC')
                              call clock_start(timer)
                              call ACPP0_THC(THCData, AuxData, CAONO, Flags)
                              call tmsg('TIME FOR ACPP0_THC', timer, 0)

                        case default
                              !
                              print*, 'not implemented'
                              stop
                              !
                        end select

                  end select

            case(JOB_TYPE_ACPP)

                  select case(Flags%ITwoEl)

                  case(1)
                        call msg("Starting ppAC with integrals in RAM")

                        select case(AuxData%omegaorders)

                        case(0)
                              call msg("Algorithm version: FULL AC")
                              call clock_start(timer0)

                              call ppAC_full(ETot, ECorr_AC, XOne, TWONO, CAONO, Flags, AuxData)

                              print*, 'RDSC ACPP ECORR', ETot+ECorr_AC+ENuc
                              print*, 'Czas na full version: ', clock_readwall(timer0)
                              stop

                        case default
                              call msg("Algorithm version: Approximated iterative OMEGA AC ")
                              call clock_start(timer0)

                              call ppAC_int_w(AuxData%HType, ECorr_AC, &
                                    ETot, ENuc, Occ, XOne, TwoNO, AuxData, AuxData%IndAux, AuxData%IndMod, &
                                    NBasis, NA, NI, NV, NInte1, NInte2, AuxData%version, Flags, AuxData%ACType)
                              print*, 'RDSC ACPP OMEGA ECORR', ETot+ECorr_AC+ENuc
                              print*, 'Czas na omega version: ', clock_readwall(timer0)                                          
                              stop                                          
                        end select
                  case(3)
                        call msg('NOT YET IMPLEMENTED, exiting')

                        select case(Flags%ICholeskyTHC)
                        case(1)
                              call msg("Starting ppAC with THC integrals THC")

                              select case(AuxData%omegaorders)
                              case(0)
                                    call msg("Algorithm version: FULL AC")
                                    call ACPP_THC(THCData, AuxData, XOne, CAONO, Flags)
                              case default
                                    call msg("Algorithm version: Approximated iterative AC ")
                              end select

                        case(0)
                              call msg("Starting ppAC with FOFO")

                              select case(AuxData%omegaorders)
                              case(0)
                                    call msg("Algorithm version: FULL AC")
                              case default
                                    call msg("Algorithm version: Approximated iterative AC ")
                              end select

                        end select

                  end select
            case(JOB_TYPE_PPERPA)

                  call pp_erpa(XOne, TWONO, CAONO, Flags, AuxData)


            end select
          end associate
      end subroutine acpp_driver


end module ppac_dispatch
