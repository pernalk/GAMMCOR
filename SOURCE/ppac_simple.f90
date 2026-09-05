module ppac_simple

      use ppac_types
      use types
      use math_constants
      use real_linalg 
      use THC_Gammcor
      use sort
      use ppac_simple_subs

      implicit none
contains

      subroutine acpp_driver_simple(THCData, AuxData, CAONO, Flags, TwoEl)


            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData
            double precision, optional, intent(in) :: TwoEl(:)
            integer :: spin_symm
            integer :: nn
            double precision, dimension(:), allocatable :: TwoNOA
            

            associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              !              AuxData%spinsep = .false.
              AuxData%spinsep = .true.
              if (AuxData%spinsep == .true.)then

!                    if ((Flags%Jobtype == JOB_TYPE_AC0PP &
!                          .or. Flags%Jobtype == JOB_TYPE_ACPP)) then 
                    print*, ''
                    print*, 'THIS IS SPIN SEPARATED ADAPTED VERSION'
                    print*, ''
                    call select_pairs(AuxData)
                    
              else
                    print*, ''
                    print*, 'THIS IS RAW SPIN VERSION'
                    print*, ''
                    call select_pairs_nospinsep(AuxData)
                    !call select_pairs_nospinsep2(AuxData)
              end if

              select case(Flags%Jobtype)

              case(JOB_TYPE_PPERPA)
                    select case(Flags%ITwoEl)
                    case(1)                          
                          call msg("Starting PPERPA incore")
                          AuxData%pperpa_print = .true.
                          nn = size(TwoEl, dim=1)                                                                            
                          allocate(TwoNOA(nn))   
                          call pperpa_incore_driver(AuxData, Flags, TwoEl, Flags%ITrpl, TwoNOA)
                          !call pperpa_incore_driver_mix(AuxData, Flags, TwoEl, 0, TwoNOA)

                          !call pperpa_incore_driver(AuxData, TwoEl, 2)
                    end select

              case(JOB_TYPE_PPERPA_RDMDUMP)
                    select case(Flags%ITwoEl)
                    case(1)
                          

                          call msg("----------------------------------------------------")
                          call msg("Starting PPERPA incore")
                          call msg("")
                          call msg("Starting singlets contributions...")
                          nn = size(TwoEl, dim=1)
                          allocate(TwoNOA(nn))
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 0, 0, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions")
                          call msg("Starting triplets contributions ABAB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 1, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions ABAB")
                          call msg("Starting triplets contributions AAAA...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 2, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions AAAA")
                          call msg("Starting triplets contributions BBBB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 3, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions BBBB")

                          call msg("Dumping pp2RDM...")
!                          call pperpa_gen_rdm(AuxData, TwoEl)
!                          print*, ''
                          print*, 'secondo'
                          call dump_ppRDM(AuxData, TwoEl)
                          call msg("DONE")                          
                    end select
              case(JOB_TYPE_AC0PP)
                    select case(Flags%ITwoEl)
                    case(1)
                          call msg("Starting AC0PP incore")
                          AuxData%pperpa_print = .true.

                          if (AuxData%spinsep == .true.)then
                                call ACPP0_incore(AuxData, TwoEl, Flags)
                          else
                                call ACPP0_incore_mix(AuxData, TwoEl, Flags)
                          end if

                    end select

              case(JOB_TYPE_ACPP)
                    select case(Flags%ITwoEl)
                    case(1)
                          call msg("Starting ACPP incore")
                          AuxData%pperpa_print = .true.
                          AuxData%ACType = 0
                          AuxData%HType = H_DYALL
                          call ppac_incore_driver(AuxData, Flags, TwoEl)

                    end select
              ! case(JOB_TYPE_ACHH)
              !       select case(Flags%ITwoEl)
              !       case(1)
              !             call msg("Starting ACHH incore")
              !             AuxData%pperpa_print = .true.
              !             AuxData%ACType = 0
              !             AuxData%HType = H_DYALL
              !             call ppac_incore_driver(AuxData, Flags, TwoEl)

              !       end select

              case(JOB_TYPE_HHERPA_RDMDUMP)
                    select case(Flags%ITwoEl)
                    case(1)
                          

                          ! call msg("----------------------------------------------------")
                          ! call msg("Starting hhERPA incore (RDM dump)")
                          ! call msg("")
                          ! call msg("Starting singlets contributions...")
                          ! nn = size(TwoEl, dim=1)
                          ! allocate(TwoNOA(nn))
                          ! call pperpa_dump_driver(AuxData, Flags, TwoEl, 0, 0, 1, TwoNOA)
                          ! call msg("...")
                          ! call msg("Done with singles contributions")
                          ! call msg("Starting triplets contributions ABAB...")
                          ! call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 1, 1, TwoNOA)
                          ! call msg("...")
                          ! call msg("Done with singles contributions ABAB")

                          ! call msg("Starting triplets contributions AAAA...")
                          ! call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 2, 1, TwoNOA)
                          ! call msg("...")
                          ! call msg("Done with singles contributions AAAA")

                          ! call msg("Starting triplets contributions BBBB...")
                          ! call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 3, 1, TwoNOA)
                          ! call msg("...")
                          ! call msg("Done with singles contributions BBBB")


                          ! call msg("Dumping hh2RDM...")                          
                          ! call dump_hhRDM(AuxData, TwoEl, 1)
                          ! call msg("DONE")

                          !-------------------------------------------------------------------------------------------------------
                          call msg("----------------------------------------------------")
                          call msg("Starting hhERPA incore (RDM dump)")
                          call msg("")
                          call msg("Starting singlets contributions...")
                          nn = size(TwoEl, dim=1)
                          allocate(TwoNOA(nn))
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 0, 0, 0, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions")
                          call msg("Starting triplets contributions ABAB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 1, 0, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions ABAB")

                          call msg("Starting triplets contributions AAAA...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 2, 0, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions AAAA")

                          call msg("Starting triplets contributions BBBB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 3, 0, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions BBBB")


                          call msg("Dumping hh2RDM...")                          
                          call dump_hhRDM(AuxData, TwoEl, 0)
                          call msg("DONE")

                          call msg("----------------------------------------------------")
                          call msg("Starting hhERPA incore (RDM dump)")
                          call msg("")
                          call msg("Starting singlets contributions...")

                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 0, 0, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions")
                          call msg("Starting triplets contributions ABAB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 1, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions ABAB")

                          call msg("Starting triplets contributions AAAA...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 2, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions AAAA")

                          call msg("Starting triplets contributions BBBB...")
                          call pperpa_dump_driver(AuxData, Flags, TwoEl, 1, 3, 1, TwoNOA)
                          call msg("...")
                          call msg("Done with singles contributions BBBB")

                          call msg("Dumping hh2RDM...")                          
                          call dump_hhRDM(AuxData, TwoEl, 1)
                          call msg("DONE")

                          
                    end select
                    
              end select

              stop
            end associate     

      end subroutine acpp_driver_simple

      subroutine select_pairs(AuxData)
            type(TACppData), intent(inout) :: AuxData
            integer :: ind, indt, indt_aa, indt_bb, i, j, ij, k

            associate(Occ=>AuxData%Occ, NBasis=>AuxData%NBasis, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m)


              print*, 'AuxData%ThrPP', AuxData%ThrPP

              allocate(AuxData%map(NBasis))
              allocate(AuxData%IndN(2,Nbasis**2))

              allocate(AuxData%IndX_s(NBasis**2))
              allocate(AuxData%IndN_s(2,Nbasis**2))

              allocate(AuxData%IndX_t(NBasis**2))
              allocate(AuxData%IndN_t(2,Nbasis**2))

              allocate(AuxData%IndX_t_aa(NBasis**2))
              allocate(AuxData%IndN_t_aa(2,Nbasis**2))

              allocate(AuxData%IndX_t_bb(NBasis**2))
              allocate(AuxData%IndN_t_bb(2,Nbasis**2))

              call imsg("NBASIS:", AuxData%NBasis)
              call imsg("NInactive:", AuxData%NI)
              call imsg("NActive:", AuxData%NA)
              call imsg("NVirt:", AuxData%NV)


              ind = 0
              indt = 0
              indt_aa = 0
              indt_bb = 0
              ij = 0

              do i = 1, Nbasis
                    do j = 1, i
                          ij = ij + 1

                          ! write(*,'(A10, 2I3)') 'wybieram pare', i, j
                          ! write(*,'(A7, F12.6, A7, F12.6)') 'occ(i)', occ(i), 'occ(j)', occ(j)
                          ! write(*,'(A25, 2F12.6)')  'Occ(i)+Occ(j)-one', Occ(i)+Occ(j)-one, AuxData%ThrPP

                          if ((abs((Occ(i)+Occ(j)-one)) > AuxData%ThrPP))then
                                ind = ind + 1
                                AuxData%IndN(1, ind) = i
                                AuxData%IndN(2, ind) = j

                                AuxData%IndN_s(1, ind) = i
                                AuxData%IndN_s(2, ind) = j
                                AuxData%IndX_s(ind) = ind
                                !print*, 'take to singlet'
                                if (j.ne.i)then
                                      indt = indt + 1
                                      AuxData%IndX_t(indt) = indt

                                      AuxData%IndN_t(1, indt) = i
                                      AuxData%IndN_t(2, indt) = j

                                      if ((abs(((n_p(i)+n_p(j))-one)) > AuxData%ThrPP))then
                                            indt_aa = indt_aa + 1
                                            AuxData%IndX_t_aa(indt_aa) = indt_aa

                                            AuxData%IndN_t_aa(1, indt_aa) = i
                                            AuxData%IndN_t_aa(2, indt_aa) = j
                                            !write(*,'(A10, F12.6)') 'accepted aaaa', n_p(i)+n_p(j)-one
                                            !print*, 'accepted aaaa', i, j, n_p(i)+n_p(j)-one
                                      end if

                                      if ((abs(((n_m(i)+n_m(j))-one)) > AuxData%ThrPP))then
                                            indt_bb = indt_bb + 1
                                            AuxData%IndX_t_bb(indt_bb) = indt_bb

                                            !write(*,'(A10, F12.6)') 'accepted bbbb', n_m(i)+n_m(j)-one
                                            AuxData%IndN_t_bb(1, indt_bb) = i
                                            AuxData%IndN_t_bb(2, indt_bb) = j
                                      end if
                                      !print*, ''
                                end if

                          end if
                    end do
              end do


              AuxData%NDim = ind
              AuxData%NDim_s = ind
              AuxData%NDim_t = indt
              AuxData%NDim_t_aa = indt_aa
              AuxData%NDim_t_bb = indt_bb

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

              call imsg("Original number of pq pairs", ij)
              call imsg("Reduced number of pq singlet pairs", ind)
              call imsg("Reduced number of pq triplet pairs", indt)
              call imsg("Reduced number of pq_aa triplet pairs", indt_aa)
              call imsg("Reduced number of pq_bb triplet pairs", indt_bb)

              write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'Occ(i)', 'Occ(j)'
              write(*, '(A43)') repeat('-', 43)

              do i = 1, min(100, AuxData%NDim_s)
                    write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                          Occ(AuxData%IndN(1,i)), Occ(AuxData%IndN(2,i))
              end do

            end associate

      end subroutine select_pairs


      subroutine select_pairs_nospinsep(AuxData)
            type(TACppData), intent(inout) :: AuxData
            integer :: ind, indt, indt_aa, indt_bb, i, j, ij, k

            associate(Occ=>AuxData%Occ, NBasis=>AuxData%NBasis, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m)


              print*, 'AuxData%ThrPP', AuxData%ThrPP

              allocate(AuxData%map(NBasis))
              allocate(AuxData%IndN(2,Nbasis**2))
              allocate(AuxData%IndX(NBasis**2))


              allocate(AuxData%IndX_t_aa(NBasis**2))
              allocate(AuxData%IndN_t_aa(2,Nbasis**2))

              allocate(AuxData%IndX_t_bb(NBasis**2))
              allocate(AuxData%IndN_t_bb(2,Nbasis**2))

              call imsg("NBASIS:", AuxData%NBasis)
              call imsg("NInactive:", AuxData%NI)
              call imsg("NActive:", AuxData%NA)
              call imsg("NVirt:", AuxData%NV)


              ind = 0              
              indt_aa = 0
              indt_bb = 0
              ij = 0


                    
              do i = 1, Nbasis
                    do j = 1, NBasis
                          ij = ij + 1

                          ! write(*,'(A10, 2I3)') 'wybieram pare', i, j
                          ! write(*,'(A7, F12.6, A7, F12.6)') 'occ(i)', occ(i), 'occ(j)', occ(j)
                          ! write(*,'(A25, 2F12.6)')  'Occ(i)+Occ(j)-one', Occ(i)+Occ(j)-one, AuxData%ThrPP

                          write(*,'(A10, 2I3, 3F12.8)') 'wybieram pare', i, j, n_p(i),n_m(j), n_p(i)+n_m(j)-one

                          if ((abs(((n_p(i)+n_m(j))-one)) > AuxData%ThrPP))then
                                print*, 'tak'
                                ind = ind + 1
                                AuxData%IndN(1, ind) = i
                                AuxData%IndN(2, ind) = j
                                AuxData%IndX(ind) = ind
                          end if
                          
                          if (i.gt.j)then

                                if ((abs(((n_p(i)+n_p(j))-one)) > AuxData%ThrPP))then
                                      indt_aa = indt_aa + 1
                                      AuxData%IndX_t_aa(indt_aa) = indt_aa
                                      
                                      AuxData%IndN_t_aa(1, indt_aa) = i
                                      AuxData%IndN_t_aa(2, indt_aa) = j
                                      !write(*,'(A10, F12.6)') 'accepted aaaa', n_p(i)+n_p(j)-one
                                end if

                                if ((abs(((n_m(i)+n_m(j))-one)) > AuxData%ThrPP))then
                                      indt_bb = indt_bb + 1
                                      AuxData%IndX_t_bb(indt_bb) = indt_bb
                                      
                                      !write(*,'(A10, F12.6)') 'accepted bbbb', n_m(i)+n_m(j)-one
                                      AuxData%IndN_t_bb(1, indt_bb) = i
                                      AuxData%IndN_t_bb(2, indt_bb) = j
                                end if
                          end if

                    end do
              end do


              AuxData%NDim = ind
              AuxData%NDim_t_aa = indt_aa
              AuxData%NDim_t_bb = indt_bb

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

              call imsg("Original number of pq pairs", ij)
              call imsg("Reduced number of pq ab ", ind)
              call imsg("Reduced number of pq_aa triplet pairs", indt_aa)
              call imsg("Reduced number of pq_bb triplet pairs", indt_bb)

              write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'n_p(i)', 'n_m(j)'
              write(*, '(A43)') repeat('-', 43)

              do i = 1, min(100, AuxData%NDim)
                    write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                          n_p(AuxData%IndN(1,i)), n_m(AuxData%IndN(2,i))
              end do

            end associate

      end subroutine select_pairs_nospinsep


      subroutine select_pairs_nospinsep2(AuxData)
            type(TACppData), intent(inout) :: AuxData
            integer :: ind, indt, indt_aa, indt_bb, i, j, ij, k

            associate(Occ=>AuxData%Occ, NBasis=>AuxData%NBasis, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m)


              print*, 'AuxData%ThrPP', AuxData%ThrPP

              allocate(AuxData%map(NBasis))
              allocate(AuxData%IndN(2,Nbasis**2))
              allocate(AuxData%IndX(NBasis**2))


              allocate(AuxData%IndX_t_aa(NBasis**2))
              allocate(AuxData%IndN_t_aa(2,Nbasis**2))

              allocate(AuxData%IndX_t_bb(NBasis**2))
              allocate(AuxData%IndN_t_bb(2,Nbasis**2))

              call imsg("NBASIS:", AuxData%NBasis)
              call imsg("NInactive:", AuxData%NI)
              call imsg("NActive:", AuxData%NA)
              call imsg("NVirt:", AuxData%NV)


              ind = 0              
              indt_aa = 0
              indt_bb = 0
              ij = 0


                    
              do i = 1, Nbasis
                    do j = 1, NBasis
                          ij = ij + 1

                          ! write(*,'(A10, 2I3)') 'wybieram pare', i, j
                          ! write(*,'(A7, F12.6, A7, F12.6)') 'occ(i)', occ(i), 'occ(j)', occ(j)
                          ! write(*,'(A25, 2F12.6)')  'Occ(i)+Occ(j)-one', Occ(i)+Occ(j)-one, AuxData%ThrPP

                          write(*,'(A10, 2I3, 3F12.8)') 'wybieram pare', i, j, n_p(i),n_m(j), n_p(i)+n_m(i)+n_m(j)+n_p(j)-two

                          if ((abs(((n_p(i)+n_m(i)+n_m(j)+n_p(j))-two)) > AuxData%ThrPP))then
                                print*, 'tak'
                                ind = ind + 1
                                AuxData%IndN(1, ind) = i
                                AuxData%IndN(2, ind) = j
                                AuxData%IndX(ind) = ind
                          end if
                          
                          if (i.gt.j)then

                                if ((abs(((n_p(i)+n_p(j))-one)) > AuxData%ThrPP))then
                                      indt_aa = indt_aa + 1
                                      AuxData%IndX_t_aa(indt_aa) = indt_aa
                                      
                                      AuxData%IndN_t_aa(1, indt_aa) = i
                                      AuxData%IndN_t_aa(2, indt_aa) = j
                                      !write(*,'(A10, F12.6)') 'accepted aaaa', n_p(i)+n_p(j)-one
                                end if

                                if ((abs(((n_m(i)+n_m(j))-one)) > AuxData%ThrPP))then
                                      indt_bb = indt_bb + 1
                                      AuxData%IndX_t_bb(indt_bb) = indt_bb
                                      
                                      !write(*,'(A10, F12.6)') 'accepted bbbb', n_m(i)+n_m(j)-one
                                      AuxData%IndN_t_bb(1, indt_bb) = i
                                      AuxData%IndN_t_bb(2, indt_bb) = j
                                end if
                          end if

                    end do
              end do


              AuxData%NDim = ind
              AuxData%NDim_t_aa = indt_aa
              AuxData%NDim_t_bb = indt_bb

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

              call imsg("Original number of pq pairs", ij)
              call imsg("Reduced number of pq ab ", ind)
              call imsg("Reduced number of pq_aa triplet pairs", indt_aa)
              call imsg("Reduced number of pq_bb triplet pairs", indt_bb)

              write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'n_p(i)', 'n_m(j)'
              write(*, '(A43)') repeat('-', 43)

              do i = 1, min(100, AuxData%NDim)
                    write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                          n_p(AuxData%IndN(1,i)), n_m(AuxData%IndN(2,i))
              end do

            end associate

      end subroutine select_pairs_nospinsep2


      subroutine pperpa_incore_driver(AuxData, Flags, TwoEl, spin_symm, TwoNOA)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(inout), optional :: TwoNOA(:)
            integer, intent(in) :: spin_symm
            double precision, dimension(:), allocatable :: Eig_i
            double precision, dimension(:,:), allocatable :: MxA, MxS
            integer :: Ndim, j, nn
            type (tclock) :: timer, timer0
            double precision :: alpha
!            double precision, dimension(:), allocatable :: TwoNOA
            type(FlagsData), intent(in) :: Flags

!            print*, AuxData%Occ
!            alpha = zero
!            AuxData%alpha = zero
            alpha = one
!            alpha = zero
!            nn = size(TwoEl, dim=1)
!            nn = size(TwoNOA, dim=1)
!            print*, nn, 'nn'
!            allocate(TwoNOA(nn))
!            print*, 'AuxData%Ninte2', AuxData%Ninte2
!            call PPERPA_init_simple(AuxData, alpha, Flags, TwoEl, TwoNOA)
            
            if(spin_symm==0)then
                  NDim = AuxData%Ndim_s

                  allocate(AuxData%vplus_s(NDim))
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))

                  allocate(AuxData%Eigs_s(NDim))
                  allocate(Eig_i(NDim))
                  allocate(AuxData%Eigvec_s(NDim, NDim))
                  call clock_start(timer0)

                  call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNO0, &
                        TwoEl, AuxData%Ndim_s,  AuxData%Ndim_s, AuxData%IndN_s, AuxData%IndN_s, AuxData%IndX_s, spin_symm, 1)
                  !print*, 'Czas na sing pperpa', clock_readwall(timer0)   
                  call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_s, Eig_i, AuxData%Eigvec_s, MxA, MxS, &
                        AuxData%vplus_s, AuxData%IndN_s, AuxData%IndAux, 1)                  

                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA'
                        do j = 1, NDim
!                              if (abs(AuxData%Eigs_s(j)).gt.1.d-5)then
                                    print*,  AuxData%Eigs_s(j), ',', Eig_i(j), ',', AuxData%vplus_s(j)
 !                             end if
                        end do
                  end if
                  deallocate(Eig_i)
            else if (spin_symm==1)then
                  
                  NDim = AuxData%Ndim_t
                  allocate(AuxData%vplus_t(NDim))
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))

                  allocate(AuxData%Eigs_t(NDim))
                  allocate(Eig_i(NDim))
                  allocate(AuxData%Eigvec_t(NDim, NDim))
                  call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNO0, &
                        TwoEl, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, spin_symm, 1)
                  if (AuxData%Ndim_t .gt. 0)then
                        call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t, Eig_i, AuxData%Eigvec_t, MxA, MxS, &
                              AuxData%vplus_t, AuxData%IndN_t, AuxData%IndAux, 1)
                  end if

                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA (Averaged)'
                        do j = 1, NDim
                              !if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
                              !      if (AuxData%vplus_t(j) ==1)then
                                          print*,  AuxData%Eigs_t(j), ',', Eig_i(j), ',',  AuxData%vplus_t(j) 
                               !     end if
                              !end if
                        end do
                  end if

                  print*, 'calling aaaa'
                  call  pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, AuxData%HNO0, &
                        TwoEl, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, 1 )
                  if (AuxData%Ndim_t .gt. 0)then
                        call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t, Eig_i, AuxData%Eigvec_t, MxA, MxS, &
                              AuxData%vplus_t, AuxData%IndN_t, AuxData%IndAux, 1)
                  end if

                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA (AAAA)'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
                                    if (AuxData%vplus_t(j) ==1)then
                                          print*,  AuxData%Eigs_t(j)!, Eig_i(j), AuxData%vplus_t(j)
                                    end if
                              end if
                        end do
                  end if

                  print*, 'calling bbbb'
                  call  pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, AuxData%HNO0, &
                        TwoEl, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, 1 )
                  
                  if (AuxData%Ndim_t .gt. 0)then

                        call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t, Eig_i, AuxData%Eigvec_t, MxA, MxS, &
                              AuxData%vplus_t, AuxData%IndN_t, AuxData%IndAux, 1)
                  end if

                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA (BBBB)'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
                                    if (AuxData%vplus_t(j) ==1)then
                                          print*,  AuxData%Eigs_t(j)!, Eig_i(j), AuxData%vplus_t(j)
                                    end if
                              end if
                        end do
                  end if

                  deallocate(Eig_i)

            end if

            
      end subroutine pperpa_incore_driver


      subroutine pperpa_incore_driver_mix(AuxData, Flags, TwoEl, spin_symm, TwoNOA)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(inout), optional :: TwoNOA(:)
            integer, intent(in) :: spin_symm
            double precision, dimension(:), allocatable :: Eig_i
            double precision, dimension(:,:), allocatable :: MxA, MxS
            integer :: Ndim, j, nn
            type (tclock) :: timer, timer0
            double precision :: alpha
!            double precision, dimension(:), allocatable :: TwoNOA
            type(FlagsData), intent(in) :: Flags

!            print*, AuxData%Occ
!            alpha = zero
!            AuxData%alpha = zero
            alpha = one
!            alpha=zero
!            AuxData%alpha = zero
!            nn = size(TwoEl, dim=1)
!            nn = size(TwoNOA, dim=1)
!            print*, nn, 'nn'
!            allocate(TwoNOA(nn))
!            print*, 'AuxData%Ninte2', AuxData%Ninte2
            call PPERPA_init_simple(AuxData, alpha, Flags, TwoEl, TwoNOA)
            
            if(spin_symm==0)then
                  NDim = AuxData%Ndim_s
                  
                  allocate(AuxData%vplus_s(2*NDim))
                  allocate(MxA(2*NDim, 2*NDim))
                  allocate(MxS(2*NDim, 2*NDim))
                  
                  allocate(AuxData%Eigs_s(2*NDim))
                  allocate(Eig_i(2*NDim))
                  allocate(AuxData%Eigvec_s(2*NDim, 2*NDim))
                  call clock_start(timer0)
                  
                  call  pperpa_incore_opshell_mix(MxA, MxS, AuxData, AuxData%HNO0, &
                        TwoEl, AuxData%Ndim_s,  AuxData%Ndim_s, AuxData%IndN_s, AuxData%IndN_s, AuxData%IndX_s, 1)
                  !print*, 'Czas na sing pperpa', clock_readwall(timer0)   
                  call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_s, Eig_i, AuxData%Eigvec_s, MxA, MxS, &
                        AuxData%vplus_s, AuxData%IndN_s, AuxData%IndAux, 1)                  
                  
                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA'
                        do j = 1, 2*NDim
                              !if (abs(AuxData%Eigs_s(j)).gt.1.d-5)then
                                    print*,  AuxData%Eigs_s(j), ',', Eig_i(j),',',  AuxData%vplus_s(j)
                              !end if
                        end do
                  end if
                  deallocate(Eig_i)
            else if (spin_symm==1)then

                  NDim = AuxData%Ndim_t_aa
                  allocate(AuxData%vplus_t_aa(NDim))
			allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))

                  allocate(AuxData%Eigs_t_aa(NDim))
			allocate(Eig_i(NDim))
                  allocate(AuxData%Eigvec_t_aa(NDim, NDim))
			call  pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, AuxData%HNO0, &
                              TwoEl, AuxData%Ndim_t_aa, AuxData%Ndim_t_aa, AuxData%IndN_t_aa, AuxData%IndN_t_aa, AuxData%IndX_t_aa, 1)
			if (AuxData%Ndim_t_aa .gt. 0)then
                              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_aa, Eig_i, AuxData%Eigvec_t_aa, MxA, MxS, &
                                    AuxData%vplus_t_aa, AuxData%IndN_t_aa, AuxData%IndAux, 1)
			end if

			print*, 'wartosci wlasne PPERPA triplet_AAAA'
                  do j = 1, NDim
                        if (abs(AuxData%Eigs_t_aa(j)).gt.1.d-5)then
                              print*,  AuxData%Eigs_t_aa(j), Eig_i(j)
                        end if
                  end do

			deallocate(Eig_i)


                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA (AAAA)'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
                                    if (AuxData%vplus_t(j) ==1)then
                                          print*,  AuxData%Eigs_t(j)!, Eig_i(j), AuxData%vplus_t(j)
                                    end if
                              end if
                        end do
                  end if
                  deallocate(MxA)
                  deallocate(MxS)
                  deallocate(Eig_i)
                  
                        NDim = AuxData%Ndim_t_bb
                        allocate(AuxData%vplus_t_bb(NDim))
                        allocate(MxA(NDim, NDim))
                        allocate(MxS(NDim, NDim))
                        
                        allocate(AuxData%Eigs_t_bb(NDim))
                        allocate(Eig_i(NDim))
                        allocate(AuxData%Eigvec_t_bb(NDim, NDim))
                        call  pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, AuxData%HNO0, &
                              TwoEl, AuxData%Ndim_t_bb, AuxData%Ndim_t_bb, AuxData%IndN_t_bb, AuxData%IndN_t_bb, AuxData%IndX_t_bb, 1)
                        if (AuxData%Ndim_t_bb .gt. 0)then
                              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_bb, Eig_i, AuxData%Eigvec_t_bb, MxA, MxS, &
                                    AuxData%vplus_t_bb, AuxData%IndN_t_bb, AuxData%IndAux, 1)
                        end if

                        print*, 'wartosci wlasne PPERPA triplet_BBBB'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t_bb(j)).gt.1.d-5)then
                                    print*,  AuxData%Eigs_t_bb(j), Eig_i(j)
                              end if
                        end do
                        
                        deallocate(Eig_i)

            end if

            
      end subroutine pperpa_incore_driver_mix



      subroutine pperpa_dump_driver(AuxData, Flags, TwoEl, spin_symm, ver, z, TwoNOA)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(inout), optional :: TwoNOA(:)
            integer, intent(in) :: spin_symm
            integer, intent(in) ::  ver, z
            double precision, dimension(:), allocatable :: Eig_i
            double precision, dimension(:,:), allocatable :: MxA, MxS
            integer :: Ndim, j, nn, i, k, l
            type (tclock) :: timer, timer0
            double precision :: alpha
            integer, external :: NAddr3            
!            double precision, dimension(:), allocatable :: TwoNOA
            type(FlagsData), intent(in) :: Flags

!            AuxData%pperpa_print = .true.
!            print*, AuxData%Occ
!            alpha = zero
            !            AuxData%alpha = zero
            if (z == 0) then
                  alpha = zero
            else            
                  alpha = one
            end if
            

!            nn = size(TwoEl, dim=1)
!            nn = size(TwoNOA, dim=1)
!            print*, nn, 'nn'
!            allocate(TwoNOA(nn))
            !            print*, 'AuxData%Ninte2', AuxData%Ninte2
            if (z == 0) then
                  print*, 'for z = 0'
                  print*, 'alpha', alpha
                  call pperpa_incore_init(AuxData, alpha, Flags, TwoEl, TwoNOA)        
                  !call PPERPA_init_simple(AuxData, alpha, Flags, TwoEl, TwoNOA)
            else
                  print*, 'for z = 1'
                  print*, 'alpha', alpha
                  AuxData%HNOA  = AuxData%HNO0
                  TWONOA = TwoEl
            end if

            ! print*, 'onel'
            ! do i = 1, AuxData%NBasis
            !       do j = 1, AuxData%NBasis
            !             if (abs(AuxData%HNO0(i,j)).gt.1.d-5)then
            !                   print*, 'hno', i, j, AuxData%HNO0(i,j)
            !             end if
            !       end do
            ! end do

            ! print*, 'twoel'
            ! do i = 1,AuxData%NBasis
            !       do j = 1, AuxData%NBasis
            !             do k = 1, AuxData%NBasis
            !                   do l = 1, AuxData%NBasis
            !                         if (abs(TwoNOA(NAddr3(i, j, k, l)) ).gt.1.d-5)then
            !                               write(*, '(A5, 4I3, F20.10)')'dups', i, j, k, l, TwoNOA(NAddr3(i, j, k, l))
            !                         end if
            !             end do
            !       end do
            ! end do
            ! end do

            if(spin_symm==0)then
                  NDim = AuxData%Ndim_s

                  allocate(AuxData%vplus_s(NDim))
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))

                  allocate(AuxData%Eigs_s(NDim))
                  allocate(Eig_i(NDim))
                  allocate(AuxData%Eigvec_s(NDim, NDim))
                  call clock_start(timer0)

                  call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNOA, &
                        TwoNOA, AuxData%Ndim_s,  AuxData%Ndim_s, AuxData%IndN_s, AuxData%IndN_s, AuxData%IndX_s, spin_symm, 1)
                  !print*, 'Czas na sing pperpa', clock_readwall(timer0)   
                  call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_s, Eig_i, AuxData%Eigvec_s, MxA, MxS, &
                        AuxData%vplus_s, AuxData%IndN_s, AuxData%IndAux, 1)                  

                  !              if (AuxData%pperpa_print) then
                  print*, 'wartosci wlasne PPERPA'
                  do j = 1, NDim
                        if (abs(AuxData%Eigs_s(j)).gt.1.d-5)then
                              print*,  AuxData%Eigs_s(j), Eig_i(j), AuxData%vplus_s(j)
                        end if
                  end do
                  !             end if
                  deallocate(Eig_i)
                  deallocate(MxA)
                  deallocate(MxS)



            else if (spin_symm==1)then

                  if (ver==1) then
                        NDim = AuxData%Ndim_t
                        allocate(AuxData%vplus_t(NDim))
                        allocate(MxA(NDim, NDim))
                        allocate(MxS(NDim, NDim))
                        
                        allocate(AuxData%Eigs_t(NDim))
                        allocate(Eig_i(NDim))
                        allocate(AuxData%Eigvec_t(NDim, NDim))
                        call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNOA, &
                              TwoNOA, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, spin_symm, 1)
                        if (AuxData%Ndim_t .gt. 0)then
                              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t, Eig_i, AuxData%Eigvec_t, MxA, MxS, &
                                    AuxData%vplus_t, AuxData%IndN_t, AuxData%IndAux, 1)
                        end if
                        !                  if (AuxData%pperpa_print) then
                        print*, 'wartosci wlasne PPERPA triplet_ABAB'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
                                    !                                   if (AuxData%vplus_t(j) ==1)then
                                    print*,  AuxData%Eigs_t(j), Eig_i(j), AuxData%vplus_t(j)
                                    !                                  end if
                              end if
                        end do
                        !               end if
                        
                        deallocate(Eig_i)
                        deallocate(MxA)
                        deallocate(MxS)


                  else if (ver==2) then
                        
                        NDim = AuxData%Ndim_t_aa
                        allocate(AuxData%vplus_t_aa(NDim))
                        allocate(MxA(NDim, NDim))
                        allocate(MxS(NDim, NDim))
                        
                        allocate(AuxData%Eigs_t_aa(NDim))
                        allocate(Eig_i(NDim))
                        allocate(AuxData%Eigvec_t_aa(NDim, NDim))
                        call  pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, AuxData%HNOA, &
                              TwoNOA, AuxData%Ndim_t_aa, AuxData%Ndim_t_aa, AuxData%IndN_t_aa, AuxData%IndN_t_aa, AuxData%IndX_t_aa, 1)
                        if (AuxData%Ndim_t_aa .gt. 0)then
                              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_aa, Eig_i, AuxData%Eigvec_t_aa, MxA, MxS, &
                                    AuxData%vplus_t_aa, AuxData%IndN_t_aa, AuxData%IndAux, 1)
                        end if

                        print*, 'wartosci wlasne PPERPA triplet_AAAA'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t_aa(j)).gt.1.d-5)then
                                    print*,  AuxData%Eigs_t_aa(j), Eig_i(j), AuxData%vplus_t_aa(j)
                              end if
                        end do
                        
                        deallocate(Eig_i)
                        deallocate(MxA)
                        deallocate(MxS)


                  else if (ver==3) then

                        NDim = AuxData%Ndim_t_bb
                        allocate(AuxData%vplus_t_bb(NDim))
                        allocate(MxA(NDim, NDim))
                        allocate(MxS(NDim, NDim))
                        
                        allocate(AuxData%Eigs_t_bb(NDim))
                        allocate(Eig_i(NDim))
                        allocate(AuxData%Eigvec_t_bb(NDim, NDim))
                        call  pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, AuxData%HNOA, &
                              TwoNOA, AuxData%Ndim_t_bb, AuxData%Ndim_t_bb, AuxData%IndN_t_bb, AuxData%IndN_t_bb, AuxData%IndX_t_bb, 1)
                        if (AuxData%Ndim_t_bb .gt. 0)then
                              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_bb, Eig_i, AuxData%Eigvec_t_bb, MxA, MxS, &
                                    AuxData%vplus_t_bb, AuxData%IndN_t_bb, AuxData%IndAux, 1)
                        end if

                        print*, 'wartosci wlasne PPERPA triplet_BBBB'
                        do j = 1, NDim
                              if (abs(AuxData%Eigs_t_bb(j)).gt.1.d-5)then
                                    print*,  AuxData%Eigs_t_bb(j), Eig_i(j), AuxData%vplus_t_bb(j)
                              end if
                        end do
                        
                        deallocate(Eig_i)
                        deallocate(MxA)
                        deallocate(MxS)


                  end if
                        

            end if

            
      end subroutine pperpa_dump_driver


      subroutine ppAC_incore_driver(AuxData, Flags, TwoEl)
        type(TACppData), intent(inout) :: AuxData
        double precision, intent(in) :: TwoEl(:)
        type(FlagsData), intent(in) :: Flags
        double precision :: W, ecorr, Wx, AC0_num
        double precision :: W_0, ecorr_min0, ppAC1
        integer :: ngrid, i, nn
        double precision, dimension(:), allocatable :: XGrid, WGrid
        double precision, dimension(:), allocatable :: TwoNOA

        call pperpa_incore_allocate(AuxData)
        nn = size(TwoEl, dim=1)                                                                            
        allocate(TwoNOA(nn))   

        AuxData%alpha = zero!0.000001!zero
        call pperpa_incore_init(AuxData, AuxData%alpha, Flags, TwoEl, TwoNOA)

        call pperpa_incore_iter(AuxData, Flags, TwoNOA)
        call ppac_incore_energy(AuxData, W_0, TwoEl)
        write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X)') &
              'iter', 'ACalpha', 'wgrid(i)', 'W'
        write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X)') &
              0, AuxData%alpha,  0.0, W_0
!        stop
        ngrid = 15
        allocate(xgrid(ngrid))
        allocate(wgrid(ngrid))

        print*, 'ngrid', ngrid

        Call GauLeg(zero, one, xgrid, wgrid, ngrid)

        ecorr = zero
        
        write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X,A15)') &
              'iter', 'ACalpha', 'wgrid(i)', 'W', 'Ecorr_i'
        do i = 1, NGrid
              
              AuxData%alpha = xgrid(i)
              call pperpa_incore_init(AuxData, AuxData%alpha, Flags, TwoEl, TwoNOA)

              call pperpa_incore_iter(AuxData, Flags, TwoNOA)
              call ppac_incore_energy(AuxData, W, TwoEl)
              ecorr = ecorr + (W-W_0) * wgrid(i)
!              print*, 'iter', i, 'alpha', AuxData%alpha
              write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8)') &
                    i, AuxData%alpha,  wgrid(i), W-W_0, ecorr
              if (i==NGrid)then
                    ppAC1 = (W-W_0)/two
              end if
        end do
        
        write(*,'(/,1X,''ECASSCF+ENuc, AC-Corr, ERPA-CASSCF'',6X,3F15.8)'), &
              AuxData%ECas, ecorr, AuxData%ECas + ecorr

        write(*,'(/,1X,''ECASSCF+ENuc, AC1-Corr, ERPA-CASSCF'',6X,3F15.8)'), &
              AuxData%ECas, ppAC1, AuxData%ECas + ppAC1
       
  end subroutine ppAC_incore_driver

  subroutine pperpa_incore_allocate(AuxData)
        
        type(TACppData), intent(inout) :: AuxData
        
        !allocate(Ints%Aints1e_aa(AuxData%NBasis, AuxData%NBasis))        
        !allocate(AuxData%int_alpha(Ints%Ints2e_dim))
        !allocate(AuxData%MxA(AuxData%NDim, AuxData%NDim))
        !allocate(AuxData%MxAt(AuxData%NDim_t, AuxData%NDim_t))
        !allocate(AuxData%MxS(AuxData%NDim, AuxData%NDim))
        !allocate(AuxData%MxSt(AuxData%NDim_t, AuxData%NDim_t))
        allocate(AuxData%Eigs_s(AuxData%NDim))
        allocate(AuxData%Eigs_t(AuxData%NDim_t))
        allocate(AuxData%Eigs_t_aa(AuxData%NDim_t))
        allocate(AuxData%Eigs_t_bb(AuxData%NDim_t))       

        allocate(AuxData%Eigvec_s(AuxData%NDim, AuxData%NDim))
        allocate(AuxData%Eigvec_t(AuxData%NDim_t, AuxData%NDim_t))
        allocate(AuxData%Eigvec_t_aa(AuxData%NDim_t, AuxData%NDim_t))
        allocate(AuxData%Eigvec_t_bb(AuxData%NDim_t, AuxData%NDim_t))

        allocate(AuxData%vplus_s(AuxData%NDim))
        allocate(AuxData%vplus_t(AuxData%NDim_t))
        allocate(AuxData%vplus_t_aa(AuxData%NDim_t))
        allocate(AuxData%vplus_t_bb(AuxData%NDim_t))

  end subroutine pperpa_incore_allocate

  subroutine pperpa_incore_iter(AuxData, Flags, TwoNOA)
        type(TACppData), intent(inout) :: AuxData
        double precision, dimension(:), intent(in) :: TwoNOA
        double precision, dimension(:), allocatable :: Eig_i
        double precision, dimension(:,:), allocatable :: MxA, MxS
        integer :: Ndim, j, nn
        type(FlagsData), intent(in) :: Flags
        
        !
        ! Solve pperpa singlet eigenvalue problem
        !
        NDim = AuxData%Ndim_s
        
        allocate(MxA(NDim, NDim))
        allocate(MxS(NDim, NDim))
        allocate(Eig_i(NDim))

!        print*, 'SINGLET CONTR'
        call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNOA, &
              TwoNOA, AuxData%Ndim_s,  AuxData%Ndim_s, AuxData%IndN_s, AuxData%IndN_s, AuxData%IndX_s, 0, 1)
        
        call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_s, Eig_i, AuxData%Eigvec_s, MxA, MxS, &
              AuxData%vplus_s, AuxData%IndN_s, AuxData%IndAux, 1)
        
        if (AuxData%pperpa_print) then
              ! write(*,'(A40, F15.8)') 'wartosci wlasne Singlet pperpa eigenvalues for alpha', AuxData%alpha
              ! do j = 1, NDim
              !       if (abs(AuxData%Eigs_s(j)).gt.1.d-5)then
              !             if (AuxData%vplus_s(j)==1)then
              !                   print*,  AuxData%Eigs_s(j), ',', Eig_i(j), ',', AuxData%vplus_s(j)
              !             end if
              !       end if
              ! end do
        end if

        deallocate(Eig_i)
        deallocate(MxA)
        deallocate(MxS)

        !
        ! Solve pperpa singlet triplet problem
        !

        NDim = AuxData%Ndim_t
        allocate(MxA(NDim, NDim))
        allocate(MxS(NDim, NDim))
        allocate(Eig_i(NDim))
!        print*, 'TRIPLET CONTR'
        call  pperpa_incore_opshell(MxA, MxS, AuxData, AuxData%HNOA, &
              TwoNOA, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, 1, 1)
        if (AuxData%Ndim_t .gt. 0)then
              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t, Eig_i, AuxData%Eigvec_t, MxA, MxS, &
                    AuxData%vplus_t, AuxData%IndN_t, AuxData%IndAux, 1)
        end if

        ! if (AuxData%pperpa_print) then
        !       write(*,'(A50, F15.8)') 'First 10 Triplet ABAB pperpa eigenvalues for alpha', AuxData%alpha
        !       do j = 1, NDim
        !             if (abs(AuxData%Eigs_t(j)).gt.1.d-5)then
        !                   if (AuxData%vplus_t(j)==AuxData%ACType)then
        !                         print*,  AuxData%Eigs_t(j), ',', Eig_i(j), ',', AuxData%vplus_t(j)
        !                   end if
        !             end if
        !       end do
        ! end if

 !       print*, 'TRIPLET CONTR AAAA'

        call  pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, AuxData%HNOA, &
              TwoNOA, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, 1 )
        if (AuxData%Ndim_t .gt. 0)then
              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_aa, Eig_i, AuxData%Eigvec_t_aa, MxA, MxS, &
                    AuxData%vplus_t_aa, AuxData%IndN_t, AuxData%IndAux, 1)
        end if

        if (AuxData%pperpa_print) then
              ! write(*,'(A50, F15.8)') 'First 10 Triplet AAAA pperpa eigenvalues for alpha', AuxData%alpha
              ! do j = 1, NDim
              !       if (abs(AuxData%Eigs_t_aa(j)).gt.1.d-5)then
              !             print*,  AuxData%Eigs_t_aa(j), ',', Eig_i(j), ',', AuxData%vplus_t_aa(j)
              !       end if
              ! end do
        end if

  !      print*, 'TRIPLET CONTR BBBB'

        call  pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, AuxData%HNOA, &
              TwoNOA, AuxData%Ndim_t, AuxData%Ndim_t, AuxData%IndN_t, AuxData%IndN_t, AuxData%IndX_t, 1 )

        if (AuxData%Ndim_t .gt. 0)then
              call nonsymmetric_eigenproblem_block_simple(AuxData%Eigs_t_bb, Eig_i, AuxData%Eigvec_t_bb, MxA, MxS, &
                    AuxData%vplus_t_bb, AuxData%IndN_t, AuxData%IndAux, 1)
        end if
        if (AuxData%pperpa_print) then
              ! write(*,'(A50, F15.8)') 'First 10 Triplet BBBB pperpa eigenvalues for alpha', AuxData%alpha
              ! do j = 1, NDim
              !       if (abs(AuxData%Eigs_t_bb(j)).gt.1.d-5)then
              !             print*,  AuxData%Eigs_t_bb(j), ',', Eig_i(j), ',', AuxData%vplus_t_bb(j)
              !       end if
              ! end do
        end if

        deallocate(Eig_i)
        deallocate(MxA)
        deallocate(MxS)



  end subroutine pperpa_incore_iter

  subroutine ppac_incore_energy(AuxData, W, TwoEl)
        type(TACppData), intent(inout) :: AuxData
        double precision, intent(out)  :: W
        double precision, dimension(:), intent(in) :: TwoEl

        integer :: p, q, r, s, i, j, k
        integer :: skip_s, skip_t, skip_aa
        double precision :: aux1, numf, contr
        double precision :: W_s, W_t, W_taa
        double precision, dimension(:,:), allocatable :: work, temp_mat
        double precision :: Npqrs
        logical :: cond
        double precision, parameter :: small_e = 1.d-2
        double precision, parameter :: big_e   = 1.d+8
        integer, external :: NAddr3

        associate( &
              NDims => AuxData%NDim, NDimt => AuxData%NDim_t, &
              IAux  => AuxData%IndAux, &
              IndNs => AuxData%IndN_s, IndNt => AuxData%IndN_t, &
              n => AuxData%occ,&          
              eigs_s => AuxData%Eigs_s, eigvec_s => AuxData%Eigvec_s, & 
              eigs_t => AuxData%Eigs_t, eigvec_t => AuxData%Eigvec_t, & 
              eigs_taa => AuxData%Eigs_t_aa,eigvec_taa => AuxData%Eigvec_t_aa, &
              vplus_s =>AuxData%vplus_s, vplus_t =>AuxData%vplus_t, vplus_taa =>AuxData%vplus_t_aa) 

          W    = zero
          W_s  = zero
          W_t  = zero
          W_taa = zero

          ! ----------------------------------------------------------------
          ! singlet contribution  (pairs p>=q, diagonal numf = sqrt(half))
          ! ----------------------------------------------------------------
          allocate(work(NDims, NDims), temp_mat(NDims, NDims))
!          print*, 'actype is', AuxData%ACType
          work = eigvec_s

          do k = 1, NDims
                if (vplus_s(k).ne.AuxData%ACType)then                  
                      work(:, k) = zero
                end if
          end do
          call real_abT(temp_mat, work, work)

          do i = 1, NDims
                p = IndNs(1, i)
                q = IndNs(2, i)


                do j = 1, NDims
                      r = IndNs(1, j)
                      s = IndNs(2, j)
                      numf = one
                      if (p==q)then
                            numf = numf * sqrt(frac12)
                      end if                      
                      if (r==s) then
                            numf = numf * sqrt(frac12)
                      end if
                      Npqrs = (n(p)+n(q)-one) * (n(r)+n(s)-one)

                      call pp_cond(AuxData, p, q, r, s, cond)
                      if (.not. cond) then
                            aux1  = TwoEl(NAddr3(r,p,s,q)) + TwoEl(NAddr3(r,q,s,p))
                            contr = numf * Npqrs * aux1 * temp_mat(i,j)
!                            if (abs(contr).gt.1.d-5)then
                            !       !write(*,'(A10, 4I5, 7F15.10)') 'contro_s', p, q,r,s,contr,Npqrs, aux1, temp_mat(i,j), eigs_s(i), eigs_s(j), W_s
 !                                 write(*,'(A10, 4I5, 6F15.10)') 'contro_s', p, q,r,s,contr,Npqrs, aux1, temp_mat(i,j), numf, W_s
  !                           end if

                            W_s = W_s + contr
                      end if
                end do
          end do
          deallocate(work, temp_mat)

          ! ----------------------------------------------------------------
          ! triplet ab contribution  (pairs p>q, no diagonal, numf=1)
          ! ----------------------------------------------------------------

          allocate(work(NDimt, NDimt), temp_mat(NDimt, NDimt))

          work = eigvec_t
          do k = 1, NDimt
                if (vplus_t(k).ne.AuxData%ACType)then                  
                      work(:, k) = zero
                end if
          end do
          call real_abT(temp_mat, work, work)

          do i = 1, NDimt
                p = IndNt(1, i)
                q = IndNt(2, i)

                do j = 1, NDimt
                      r = IndNt(1, j)
                      s = IndNt(2, j)

                      Npqrs = (n(p)+n(q)-one) * (n(r)+n(s)-one)

                      call pp_cond(AuxData, p, q, r, s, cond)
                      if (.not. cond) then

                            aux1  = TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p))
                            contr = Npqrs * aux1 * temp_mat(i,j)
                            ! if (abs(contr).gt.1.d-5)then
                            !       write(*,'(A10, 4I5, 7F15.10)') 'contro_t', p, q,r,s,contr,Npqrs, aux1, temp_mat(i,j), eigs_t(i), eigs_t(j), W_t
                            ! end if


                            W_t = W_t + contr
                      end if
                end do
          end do

          work = zero
          temp_mat = zero
          deallocate(work)
          deallocate(temp_mat)
          ! ----------------------------------------------------------------
          ! triplet aa contribution  (same pair index as t; factor 2 for bb)
          ! ----------------------------------------------------------------

          allocate(work(NDimt, NDimt), temp_mat(NDimt, NDimt))

          work = eigvec_taa
          do k = 1, NDimt
                if (vplus_taa(k).ne.AuxData%ACType)then                  
                      work(:, k) = zero
                end if
          end do
          call real_abT(temp_mat, work, work)

          do i = 1, NDimt
                p = IndNt(1, i)
                q = IndNt(2, i)

                do j = 1, NDimt
                      r = IndNt(1, j)
                      s = IndNt(2, j)

                      Npqrs = (one-n(p)-n(q)) * (one - n(r)-n(s)) 

                      call pp_cond(AuxData, p, q, r, s, cond)
                      if (.not. cond) then
                            
                            aux1  = TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p))
                            contr = Npqrs * aux1 * temp_mat(i,j)
                            ! if (abs(contr).gt.1.d-5)then
                            !       write(*,'(A10, 4I5, 7F15.10)') 'contro_ta', p, q,r,s,contr,Npqrs, aux1, temp_mat(i,j), eigs_taa(i), eigs_taa(j), W_taa
                            ! end if

                            W_taa = W_taa + contr
                      end if
                end do
          end do

          deallocate(work, temp_mat)

          W = W_s + W_t + two * W_taa

!          write(*,'(A30, 4F20.10)') 'ppac_incore_energy: W_s, W_t, W_aa =', W_s, W_t, W_taa, W


        end associate
  end subroutine ppac_incore_energy

  subroutine pp_cond(AuxData, p, q, r, s, cond)
        type(TACppData), intent(in)           :: AuxData
        integer, intent(in)                   :: p, q, r, s
        logical, intent(out)                  :: cond

        associate(IAux=>AuxData%IndAux)
          if (AuxData%HType == H_DYALL) then

                cond = (IAux(p)==IAux(q) .and. IAux(p)==IAux(r) .and. &
                      IAux(p)==IAux(s) .and. IAux(p)==1)
          else

                cond = (IAux(p)==IAux(q) .and. IAux(p)==IAux(r) .and. IAux(p)==IAux(s))
          end if
        end associate
  end subroutine pp_cond



  subroutine pperpa_incore_init(AuxData, ACAlpha, Flags, TwoEl, TwoNOA)        
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            double precision, dimension(:),      intent(out) :: TwoNOA
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags

            integer :: twoint_dim, r, l
            integer :: i, j, k, kl, ij, t
            double precision :: temp
            integer, external :: NAddr3

            associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux)


              if (.not. allocated(AuxData%HNOA)) allocate(AuxData%HNOA(NBasis, NBasis))

              AuxData%HNOA = ACAlpha * AuxData%HNO0

              do i = 1, NBasis
                    do j = 1, NBasis

                          if (IAux(i)==IAux(j))then
                                AuxData%HNOA(i,j) = AuxData%HNOA(i,j) + (One-ACAlpha) * AuxData%HNO0(i,j)

                                if (IAux(i)==1)then
                                      l = NI
                                else
                                      l = NIA
                                end if
                                temp = zero
                                do t = 1, l
                                      temp = temp + occ(t)* (two*TwoEl(gmap(t,t,i,j))-TwoEl(gmap(t,i,t,j)))
                                end do

                                AuxData%HNOA(i,j) = AuxData%HNOA(i,j) + (One-ACAlpha) * temp

                          end if
                    end do
              end do

            TwoNOA = TwoEL
            ij = 0
            do i = 1, NBasis
                  do j = 1, i
                        ij = ij + 1
                        kl = 0
                        do k = 1, Nbasis
                              do l = 1, k
                                    kl=kl+1
                                     if ((IAux(i)==IAux(j)).and.(IAux(i)==IAux(k)).and.IAux(i)==IAux(l).and.IAux(i)==1)then
                                           TwoNOA(gmap(i, j, k, l)) = TwoEl(gmap(i, j, k, l))
                                    else
                                          TwoNOA(gmap(i, j, k, l)) = ACAlpha * TwoEl(gmap(i, j, k, l))
                                    end if
                              end do
                        end do
                  end do
            end do
          end associate
    end subroutine pperpa_incore_init


    
  ! subroutine erpa_init_incore(AuxData,  Flags, Ints)
  !       type(TACppData), intent(inout) :: AuxData
  !       type(TInts), intent(inout) :: Ints
  !       type(FlagsData), intent(in) :: Flags

  !       integer :: i, j, k,  t, l
  !       integer :: i0, i1
  !       double precision :: temp
  !       integer(I8) :: nb, npair, pq, rs
  !       integer :: p,q,r,s
  !       integer(I8) :: a1,b1,a2,b2, a,b, c,d, idx


  !       associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, &
  !             NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
  !             NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux, &
  !             ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
  !             ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
  !             int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)


  !         Aints1e_aa = alpha * AuxData%HNO0

  !         do i = 1, NBasis
  !               do j = 1, NBasis
  !                     if (IAux(i) == IAux(j)) then
  !                           Aints1e_aa(i, j) = Aints1e_aa(i, j) + (One - alpha) * AuxData%HNO0(i,j)

  !                           if (AuxData%HType == H_DYALL) then
  !                           !      print*, 'dyal'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
  !                                 if (IAux(i) == 1) then
  !                                       i0 = 1
  !                                       i1 = NI
  !                                 else
  !                                       i0 = 1
  !                                       i1 = NIA
  !                                 end if
  !                           else if (AuxData%HType == H_GPF) then
  !                                 if (IAux(i) == 0) then
  !                                       i0 = NI + 1
  !                                       i1 = NIA
  !                                 else if (IAux(i) == 1) then
  !                                       i0 = 1
  !                                       i1 = NI
  !                                 else
  !                                       i0 = 1
  !                                       i1 = NIA
  !                                 end if
  !                           end if

  !                           temp = zero
  !                           do t = i0, i1
  !                                 temp = temp + occ(t) * &
  !                                       ( ints2e_aa(gmap_4fold(t,t,i,j, NBasis)) &
  !                                       + ints2e_ab(gmap_4fold(t,t,i,j, NBasis)) &
  !                                       - ints2e_aa(gmap_4fold(t,i,t,j, NBasis)) )
  !                           end do

  !                           Aints1e_aa(i, j) = Aints1e_aa(i, j) + (One - alpha) * temp
  !                     end if
  !               end do
  !         end do

  !         nb    = int(NBasis, I8)
  !         npair = nb*nb

  !         AuxData%int_alpha = 1

  !         do pq = 1_I8, npair
  !               q = int((pq-1_I8)/nb + 1_I8)
  !               p = int(pq - int(q-1,I8)*nb)

  !               do rs = pq, npair
  !                     s = int((rs-1_I8)/nb + 1_I8)
  !                     r = int(rs - int(s-1,I8)*nb)

  !                     a1 = pq; b1 = rs
  !                     a  = min(a1,b1); b = max(a1,b1)

  !                     a2 = int(q,I8) + (int(p,I8)-1_I8)*nb
  !                     b2 = int(s,I8) + (int(r,I8)-1_I8)*nb
  !                     c  = min(a2,b2); d = max(a2,b2)

  !                     if ( (c < a) .or. (c == a .and. d < b) ) cycle

  !                     idx = gmap_4fold(r,s,p,q,NBasis)

  !                     if (AuxData%HType == H_DYALL) then
  !                           if (.not. (IAux(p) == 1 .and. IAux(q) == 1 .and. &
  !                                 IAux(r) == 1 .and. IAux(s) == 1)) then
  !                                 AuxData%int_alpha(idx) = 0
  !                           end if

  !                     else if (AuxData%HType == H_GPF) then
  !                           if (.not. (IAux(p) == IAux(q) .and. &
  !                                 IAux(q) == IAux(r) .and. &
  !                                 IAux(r) == IAux(s))) then
  !                                 AuxData%int_alpha(idx) = 0
  !                           end if
  !                     end if
  !               end do
  !         end do

  !       end associate
  ! end subroutine erpa_init_incore
end module ppac_simple
