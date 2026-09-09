module ppac0_thc


      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals      
      use types 
      use lin
      use clock

      use math_constants
      use pp_utils
      use print_utils
      use ppac_types
      use thc_auto_p3
      use thc_auto_p0
      use thc_auto_p2
      use thc_energy
      use tran
      use omp_lib

      use interface_pp



      implicit none

      integer, parameter,private  :: A0oo=1, A0vv=2, A0aa=3, A0va=4, A0oa=5
      integer, parameter, private :: A1vvoo=6, A1vvao=7, A1vvaa=8, A1vaoo=9
      integer, parameter, private :: A1vaao=10, A1vaaa=11, A1aaoo=12, A1aaao=13
      integer, parameter, private :: A1vaao2 = 14

      integer, parameter, private :: tdebug = 0
      integer, parameter, private :: mdebug = 0
      integer, parameter, private :: mdebugg = 1
      integer, parameter, private :: mnormal = 50
      integer, parameter, private :: mverbose = 10
      integer, parameter, private :: merror = 100

      integer, parameter, private :: msgthr = 10
      character(len=8), dimension(14), parameter :: BlockName = (/ &
            'A0oo    ', 'A0vv    ', 'A0aa    ', 'A0va    ', 'A0oa    ', &
            'A1vvoo  ', 'A1vvao  ', 'A1vvaa  ', 'A1vaoo  ', 'A1vaao  ', &
            'A1vaaa  ', 'A1aaoo  ', 'A1aaao  ', 'A1vaao2 ' /)

contains


      subroutine ACPP_THC(THCData, AuxData, XOne, CAONO_IN,Flags)
            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use sys_definitions
            use gammcor_integrals
            use Auto2eInterface


            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:), intent(in) :: XOne
            double precision, dimension(:,:), intent(in) :: CAONO_IN
            type(FlagsData), intent(in) :: Flags

            type(TTHCData), intent(inout) :: THCData
            !double precision, allocatable :: Xga(:,:), Zgk(:,:)
            !integer :: NChol, NTHC

            type(TRdmData) :: RdmData
            type(TAC0Block), dimension(1) :: ACBlocks

            integer, dimension(:,:), allocatable :: posS, posT

            integer :: i
            double precision, dimension(:,:), allocatable :: HNO0, HNO, CAONO
            double precision, dimension(:,:), allocatable :: Aux3A, Aux3B, Aux3X
            double precision, dimension(:,:), allocatable :: AuxII, AuxAA, BuxII, BuxAA, work
            double precision, dimension(:,:), allocatable :: HNO0_THC         
            double precision :: ACAlpha
            integer, parameter :: ExternalOrdering = ORBITAL_ORDERING_DALTON                                                                                                                                                         
            !integer, parameter :: ExternalOrdering = ORBITAL_ORDERING_PYSC

            integer q, s

              associate(NBasis=>AuxData%NBasis, map=>AuxData%map, Occ=>AuxData%Occ, NI=>AuxData%NI, NV=>AuxData%NV)
                THCData%ExternalOrdering = ORBITAL_ORDERING_DALTON
                
              ACAlpha = One!0.0356d+0!One !zero

              allocate(HNO0(NBasis, NBasis))
              allocate(HNO(NBasis, NBasis))
              allocate(THCData%fij(NI))
              allocate(THCData%fvw(NV))
              
              allocate(CAONO(NBasis, NBasis))
              allocate(HNO0_THC(NBasis, NBasis))
              call HNO0_init(AcAlpha, XOne, HNO0, NBasis, AuxData%IndAux)
              allocate(posS(NBasis, NBasis))
              allocate(posT(NBasis, NBasis))
              posS = 0
              posT = 0
              call init_ACBlocks(AuxData, ACBlocks, posS, posT)


              call init_rdm(RdmData, AuxData)


              allocate(Aux3A(NBasis, NBasis))
              allocate(Aux3B(NBasis, NBasis))
              allocate(Aux3X(NBasis, NBasis))
              
              allocate(AuxII(NBasis, NBasis))
              allocate(AuxAA(NBasis, NBasis))
              allocate(BuxII(NBasis, NBasis))
              allocate(BuxAA(NBasis, NBasis))

              AuxII = zero
              BuxII = zero
              AuxAA = zero
              BuxAA = zero
              Aux3A = zero
              Aux3B = zero
              Aux3X = zero


              CAONO= CAONO_IN
              call THC_init(Flags, AuxData, THCData, HNO0, HNO0_THC, CAONO, RdmData)

              call THC_int_loop(THCData, ACBlocks, AuxII, AuxAA, BuxII, BuxAA, Aux3X, Aux3B, &
                    Flags, AuxData, AuxData%Occ, map, RdmData, posS, posT, AcAlpha)



              call update_HNO(AuxData, HNO0, HNO, AuxII, AuxAA, BuxII, BuxAA, ACAlpha, 1)
              deallocate(HNO0)
              call PPERPA_THC(AuxData, ACBlocks(1), HNO, posS, posT, AuxII, AuxAA, BuxII, BuxAA,Aux3X, Aux3B, ACAlpha, AuxData%Ndim, AuxData%Ndim, AuxData%IndN, AuxData%IndN)


            end associate


      end subroutine ACPP_THC

      subroutine ACPP0_THC(THCData, AuxData, CAONO_IN, Flags)
            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use sys_definitions
            use gammcor_integrals
            use Auto2eInterface


            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO_IN
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData

            type(TRdmData) :: RdmData
            type(TAC0Block), dimension(14) :: ACBlocks

            integer, dimension(:,:), allocatable :: posS, posT

            integer :: i, j
            double precision, dimension(:,:), allocatable :: HNO0, HNO, CAONO
            double precision, dimension(:,:), allocatable :: HNO0_THC
            double precision, dimension(:,:), allocatable :: Aux3A, Aux3B, Aux3X
            double precision, dimension(:,:), allocatable :: AuxII, AuxAA, BuxII, BuxAA, work            
            double precision :: ACAlpha
            integer :: this, k
            type (tclock) :: timer
            integer :: ExternalOrdering

            !double precision, dimension(:,:), allocatable :: HNO0_PYSCF, CAONO_PYSCF            
            
            
            if (AuxData%PYSCF==1)then
                  !call read_PYSCF(InputData, CAONO_PYSCF)
                  call print_info('External ordering', 'PySCF')
                  !THCData%ExternalOrdering = ORBITAL_ORDERING_PYSCF
            else if (AuxData%ORCA==1)then
                  !THCData%ExternalOrdering = ORBITAL_ORDERING_ORCA
                  call print_info('External ordering', 'ORCA')
            else  if (AuxData%DALTON==1)then
                  THCData%ExternalOrdering = ORBITAL_ORDERING_DALTON
            end if
            

            associate(NBasis=>AuxData%NBasis, map=>AuxData%map, Occ=>AuxData%Occ, NI=>AuxData%NI, NV=>AuxData%NV)

              ACAlpha = one !because we are simultaneously calculating A1. 

              allocate(HNO0(NBasis, NBasis))
              call CalcMem(HNO0, 'HNO0')
              allocate(HNO(NBasis, NBasis))
              call CalcMem(HNO, 'HNO')
              !allocate(THCData%fij(NI))
              
              !allocate(THCData%fvw(NV))
              
              allocate(CAONO(NBasis, NBasis))
              call CalcMem(CAONO, 'CAONO')

              
              CAONO = CAONO_IN
              allocate(HNO0_THC(Nbasis, NBasis))
              !call CalcMem(HNO0_THC, 'HNO_THC')

              allocate(posS(NBasis, NBasis))
              allocate(posT(NBasis, NBasis))
              call CalcMemI(posS, 'posS')
              call CalcMemI(posT, 'posT')
              posS = 0
              posT = 0


              call init_AC0Blocks(AuxData, ACBlocks, posS, posT)


              call init_rdm(RdmData, AuxData)

              allocate(Aux3A(NBasis, NBasis))
              allocate(Aux3B(NBasis, NBasis))
              allocate(Aux3X(NBasis, NBasis))

              allocate(AuxII(NBasis, NBasis))
              allocate(AuxAA(NBasis, NBasis))
              allocate(BuxII(NBasis, NBasis))
              allocate(BuxAA(NBasis, NBasis))

              call CalcMem(Aux3A, 'Auxall6')


              AuxII = zero
              BuxII = zero
              AuxAA = zero
              BuxAA = zero
              Aux3A = zero
              Aux3B = zero
              Aux3X = zero

              call clock_start(timer)

              if (AuxData%PYSCF==1.or.AuxData%ORCA ==1)then
                    !call THC_init(Flags, AuxData, THCData, AuxData%HNO0, HNO0_THC, CAONO, RdmData)
              else
                    call HNO0_init(AcAlpha, AuxData%XOne, HNO0, NBasis, AuxData%IndAux)
                    !call THC_init(Flags, AuxData, THCData, HNO0, HNO0_THC, CAONO, RdmData)
              end if
              
              ! call tmsg('TIME FOR THC INIT', timer, tdebug)
              ! call clock_start(timer)
              

              ! print*, 'hno po transform'
              ! do i = 1, NBasis
              !       do j = 1, NBasis
              !             if (abs(abs(HNO0(i,j)) - abs(HNO0_THC(i, j))).gt.1.d-2)then
              !             !if (abs(HNO0(i,j)).gt.1.d-5.or.abs(HNO0_THC(i, j)).gt.1.d-5)then
              !                   print*, i, j, HNO0(i,j), HNO0_THC(i, j), abs(abs(HNO0(i,j)) - abs(HNO0_THC(i, j)))
              !             end if
              !       end do
              ! end do
              ! stop
              call THC_int_loop(THCData, ACBlocks, AuxII, AuxAA, BuxII, BuxAA, Aux3X, Aux3B, &
                    Flags, AuxData, Occ, map, RdmData, posS, posT, AcAlpha)
              call tmsg('TIME FOR THC LOOPS', timer, tdebug)
              call clock_start(timer)


!              call update_HNO(AuxData, HNO0,  HNO, AuxII, AuxAA, BuxII, BuxAA, ACAlpha, 0)
              !call update_HNO(AuxData, HNO0_THC,  HNO, AuxII, AuxAA, BuxII, BuxAA, ACAlpha, 0)              
              !call update_HNO_AC0(AuxData, HNO0, HNO, fij, fvw, AuxII, BuxII)

              
              call update_HNO_AC0(AuxData, HNO, THCData%fij, THCData%fvw, AuxII, BuxII)
              call tmsg('TIME FOR UPDATE HNO', timer, tdebug)
              call clock_start(timer)

              call print_section('ppAC0 response blocks')

              do i = 1,2
                    call clock_start(timer)
                    call PPERPA0_THC_OOVV(AuxData, ACBlocks(i), HNO, posS, posT, AuxAA, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(i)%NdimS, ACBlocks(i)%IndNS)
                    call tmsg('Hessian '//ACBlocks(i)%name, timer, tdebug)
                    call clock_start(timer)
                    call eigs(AuxData, ACBlocks(i), i)
                    call tmsg('Diagonalization '//ACBlocks(i)%name, timer, tdebug)
              end do

              call clock_start(timer)
              call PPERPA0_THC_AA(AuxData, ACBlocks(A0aa), HNO, posS, posT, AuxAA, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(A0aa)%NdimS, ACBlocks(A0aa)%IndNS)
              call tmsg('Hessian '//ACBlocks(A0aa)%name, timer, tdebug)
              
              call clock_start(timer)
              call eigs(AuxData, ACBlocks(A0aa), A0aa)
              call tmsg('Diagonalization '//ACBlocks(A0aa)%name, timer, tdebug)

              call clock_start(timer)
              call PPERPA0_THC_VA_Eigs(AuxData, ACBlocks(A0va), HNO, AuxAA, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(A0va)%NdimS, ACBlocks(A0va)%IndNS)
              call tmsg('Hessian and diagonalization '//ACBlocks(A0va)%name, timer, tdebug)

              call clock_start(timer)
              call PPERPA0_THC_OA_Eigs(AuxData, ACBlocks(A0oa), HNO, AuxAA, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(A0oa)%NdimS, ACBlocks(A0oa)%IndNS)
              call tmsg('Hessian and diagonalization '//ACBlocks(A0oa)%name, timer, tdebug)


              do i = 9, 13                    
                    call clock_start(timer)
                    !call PPERPA_THC(AuxData, ACBlocks(i), HNO0, posS, posT, AuxII, AuxAA, BuxII, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(i)%Ndim1S, ACBlocks(i)%Ndim2S, ACBlocks(i)%IndN1S, ACBlocks(i)%IndN2S)
                    call PPERPA_THC(AuxData, ACBlocks(i), AuxData%HNO0, posS, posT, AuxII, AuxAA, BuxII, BuxAA,Aux3X, Aux3B, ACAlpha, ACBlocks(i)%Ndim1S, ACBlocks(i)%Ndim2S, ACBlocks(i)%IndN1S, ACBlocks(i)%IndN2S)
                    call tmsg('Hessian '//ACBlocks(i)%name, timer, tdebug)
              end do
              
              call clock_start(timer)
              call THC_energy_loop(THCData, ACBlocks, Flags, posS, posT, AuxData)
              call tmsg('Correlation energy', timer, tdebug)
              !call THC_int_loop(AC0Block, AuxData, CAONO_IN, Flags, RdmData)
              !call 


            end associate


      end subroutine ACPP0_THC


      subroutine init_rdm(RdmData, AuxData)            
            type(TRdmData), intent(out) :: RdmData
            type(TACppData), intent(in) :: AuxData
            integer :: NRDM2Act
            integer :: i, j, k, l
            integer :: iunit

            associate(Occ=>AuxData%Occ, map=>AuxData%map, NI=>AuxData%NI, NIA=>AuxData%NIA, NA=>AuxData%NA, NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux)
              NRDM2Act = NA**2*(NA**2+1)/2



              allocate(RdmData%rdm2_pp_act(NA,NA,NA,NA))
              allocate(RdmData%rdm2_pm_act(NA,NA,NA,NA))

              allocate(RdmData%rdm2_pp_12_act(NA,NA,NA,NA))
              allocate(RdmData%rdm2_pm_12_act(NA,NA,NA,NA))

              allocate(RdmData%rdm2_pp_13_act(NA,NA,NA,NA))
              allocate(RdmData%rdm2_pm_13_act(NA,NA,NA,NA))
              call CalcMem4(RdmData%rdm2_pp_13_act, 'rdm6')


              ! if (AuxData%PYSCF ==1)then
              !       print*, 'reading rdms pyscf'
                    
              !       RdmData%rdm2_pp_act = AuxData%rdm2_pp
              !       RdmData%rdm2_pm_act = AuxData%rdm2_pm

              !       ! RdmData%rdm2_pp_act = zero
              !       ! RdmData%rdm2_pm_act = zero

              !       ! iunit = 20
              !       ! open(unit=iunit, file='rdm2_aaaa.bin', status='old', access='stream', form='unformatted')
              !       ! read(iunit) RdmData%rdm2_pp_act
              !       ! close(iunit)

              !       ! iunit = 21
              !       ! open(unit=iunit, file='rdm2_abab.bin', status='old', access='stream', form='unformatted')
              !       ! read(iunit) RdmData%rdm2_pm_act
              !       ! close(iunit)
              !        do l=NI+1, NIA
              !               do k=NI+1, NIA
              !                     do j=NI+1, NIA
              !                           do i=NI+1, NIA
              !                                 RdmData%rdm2_pp_12_act(map(i), map(k), map(j), map(l)) = RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))
              !                                 RdmData%rdm2_pm_12_act(map(i), map(k), map(j), map(l)) = RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))
                                              
              !                                 RdmData%rdm2_pp_13_act(map(i), map(l), map(k), map(j)) = RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))
              !                                 RdmData%rdm2_pm_13_act(map(i), map(l), map(k), map(j)) = RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))
              !                           end do
              !                     end do
              !               end do
              !         end do
                      
                      

        if (AuxData%ORCA==1 .or. AuxData%PYSCF==1)then
                                                        
              

                      RdmData%rdm2_pp_act = AuxData%rdm2_pp
                      RdmData%rdm2_pm_act = AuxData%rdm2_pm

                      do l=NI+1, NIA
                            do k=NI+1, NIA
                                  do j=NI+1, NIA
                                        do i=NI+1, NIA
                                              RdmData%rdm2_pp_12_act(map(i), map(k), map(j), map(l)) = AuxData%rdm2_pp(map(i), map(j), map(k), map(l))
                                              RdmData%rdm2_pm_12_act(map(i), map(k), map(j), map(l)) = AuxData%rdm2_pm(map(i), map(j), map(k), map(l))

                                              RdmData%rdm2_pp_13_act(map(i), map(l), map(k), map(j)) = AuxData%rdm2_pp(map(i), map(j), map(k), map(l))
                                              RdmData%rdm2_pm_13_act(map(i), map(l), map(k), map(j)) = AuxData%rdm2_pm(map(i), map(j), map(k), map(l))
                                        end do
                                  end do
                            end do
                      end do


                else
                      
                      allocate (RdmData%R00(NRDM2Act))
                      allocate (RdmData%R11(NRDM2Act))

                      
                      RdmData%R00 = Zero
                      RdmData%R11 = Zero


                      print*, 'reading rdms DALTON'
                      call read_2rdm("rdm2.dat", RdmData%R00, NA)
                      call read_2rdm("rdms2.dat", RdmData%R11, NA)
                    
                      associate(R00 => RdmData%R00, R11 =>RdmData%R11)

                      do l=NI+1, NIA
                            do k=NI+1, NIA
                                  do j=NI+1, NIA
                                        do i=NI+1, NIA
                                              RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))  = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 0)
                                              RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))  = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 1)


                                              RdmData%rdm2_pp_12_act(map(i), map(k), map(j), map(l)) = RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))
                                              RdmData%rdm2_pm_12_act(map(i), map(k), map(j), map(l)) = RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))

                                              RdmData%rdm2_pp_13_act(map(i), map(l), map(k), map(j)) = RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))
                                              RdmData%rdm2_pm_13_act(map(i), map(l), map(k), map(j)) = RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))
                                        end do
                                  end do
                            end do
                      end do
                    end associate
              end if

              
              ! print*, 'this is rdmpp'
              ! do l=NI+1, NIA
              !       do k=NI+1, NIA
              !             do j=NI+1, NIA
              !                   do i=NI+1, NIA
              !                         if (AuxData%ORCA==1)then
              !                               if (abs(AuxData%rdm2_pp(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), AuxData%rdm2_pp(map(i), map(j), map(k), map(l))
              !                               end if
              !                         else
              !                               if (abs(RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))
              !                               end if
              !                         end if
              !                   end do
              !             end do
              !       end do
              ! end do
              
              ! print*, 'this is rdmpm'
              ! do l=NI+1, NIA
              !       do k=NI+1, NIA
              !             do j=NI+1, NIA
              !                   do i=NI+1, NIA
              !                         if (AuxData%ORCA==1)then
              !                               if (abs(AuxData%rdm2_pm(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), AuxData%rdm2_pm(map(i), map(j), map(k), map(l))
              !                               end if
              !                         else
                                            
              !                               if (abs(RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l))
              !                               end if
              !                         end if
              !                   end do
              !             end do
              !       end do
              ! end do
              
              ! print*, 'this is sum'
              ! do l=NI+1, NIA
              !       do k=NI+1, NIA
              !             do j=NI+1, NIA
              !                   do i=NI+1, NIA
              !                          if (AuxData%ORCA==1)then
              !                               if (abs(AuxData%rdm2_pm(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), two*(AuxData%rdm2_pp(map(i), map(j), map(k), map(l))+AuxData%rdm2_pm(map(i), map(j), map(k), map(l)))
              !                               end if
              !                         else
              !                               if (abs(RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))).gt.1.d-5)then
              !                                     write(*, '(4I5, 2F20.15)')map(i), map(j), map(k), map(l), two*(RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l))+RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l)))
              !                               end if
              !                         end if
              !                   end do
              !             end do
              !       end do
              ! end do
              

                      print*, ''

                    end associate

      end subroutine init_rdm


      subroutine init_ACBlocks(AuxData, ACBlock, posS, posT)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), dimension(:), intent(inout) :: ACBlock
            integer, dimension(:,:), intent(inout) :: posS, posT
            integer, parameter :: occ = 0, act = 1, virt = 2
            integer :: i, p, q
            integer :: NDimS, NDimT

            NDimS =  AuxData%NDim_s
            NDimT = AuxData%NDim_t


            ACBlock(1)%NdimS = NDimS
            ACBlock(1)%NdimT = NDimT
            ACBlock(1)%order = 0

            allocate(ACBlock(1)%IndNS(2, ACBlock(1)%NdimS))
            allocate(ACBlock(1)%IndNT(2, ACBlock(1)%NdimT))

            allocate(ACBlock(1)%ASing(ACBlock(1)%NdimS, ACBlock(1)%NdimS))
            allocate(ACBlock(1)%ATrip(ACBlock(1)%NdimT, ACBlock(1)%NdimT))
            allocate(ACBlock(1)%ATripA(ACBlock(1)%NdimT, ACBlock(1)%NdimT))
            allocate(ACBlock(1)%vPlusS(ACBlock(1)%NdimS))
            allocate(ACBlock(1)%vPlusT(ACBlock(1)%NdimT))
            allocate(ACBlock(1)%vPlusTA(ACBlock(1)%NdimT))
            allocate(ACBlock(1)%EigS(ACBlock(1)%NdimS))
            allocate(ACBlock(1)%EigT(ACBlock(1)%NdimT))
            allocate(ACBlock(1)%EigTA(ACBlock(1)%NdimT))

            call CalcMemI(ACBlock(1)%IndNS, 'ACBlock(1)%IndNS')
            call CalcMemI(ACBlock(1)%IndNT, 'ACBlock(1)%IndNT')

            call CalcMem(ACBlock(1)%ASing, 'ACBlock(1)%ASing')
            call CalcMem(ACBlock(1)%ATrip, 'ACBlock(1)%ATrip')
            call CalcMem(ACBlock(1)%ATripA, 'ACBlock(1)%ATripA')



            ACBlock(1)%IndNS =AuxData%IndN_s(:, 1:NDimS)
            ACBlock(1)%IndNT = AuxData%IndN_t(:, 1:NDimT)

            ACBlock(1)%ASing = zero
            ACBlock(1)%ATrip = zero
            ACBlock(1)%ATripA = zero

            do i = 1, AuxData%NDim_s
                  p = AuxData%IndN_s(1, i)
                  q = AuxData%IndN_s(2, i)

                  posS(p, q) = i
            end do

            do i = 1,	AuxData%NDim_t
                  p =	AuxData%IndN_t(1, i)
                  q =	AuxData%IndN_t(2, i)
                  posT(p, q) = i
            end do



      end subroutine init_ACBlocks

      subroutine init_AC0Blocks(AuxData, AC0Block, posS, posT)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), dimension(14), intent(inout) :: AC0Block
            integer, dimension(:,:), intent(inout) :: posS, posT
            integer,  parameter :: occ = 0, act = 1, virt = 2
            integer :: i

            ! Blocks with A0 data, eigenvectors and eigenvalues
            do i = 1, 5
                  AC0Block(i)%order = 0
            end do

            ! blocks with A1 data
            do i = 6, 14
                  AC0Block(i)%order = 1
            end do


            posS = zero
            posT = zero
            associate(IndAux=>AuxData%Indaux, IndN=>AuxData%IndN, NDim=>AuxData%NDim)

              call InitA0Blocks(AC0Block(1), AuxData, IndAux, IndN, NDim, posS, posT, occ, occ)
              call InitA0Blocks(AC0Block(2), AuxData, IndAux, IndN, NDim, posS, posT, virt, virt)
              call InitA0Blocks(AC0Block(3), AuxData, IndAux, IndN, NDim, posS, posT, act, act)
              call InitA0Blocks(AC0Block(4), AuxData, IndAux, IndN, NDim, posS, posT, virt, act)
              call InitA0Blocks(AC0Block(5), AuxData, IndAux, IndN, NDim, posS, posT, occ, act)

              do i = 1, 14
                    AC0Block(i)%name = BlockName(i)
              end do

              call print_section('ppAC0 block dimensions')
              write(*,'(2X,A10,2X,A12,2X,A12)') 'Block', 'Singlet', 'Triplet'
              write(*,'(2X,A)') repeat('-', 40)
              do i = 1, 5
                    write(*,'(2X,A10,2X,I12,2X,I12)') trim(AC0Block(i)%name), &
                          AC0Block(i)%NDimS, AC0Block(i)%NDimT
              end do

              ! These blocks are too large and are calculated OTF
              
              !call InitA1Blocks(AC0Block, AuxData, A1vvoo, A0vv, A0oo)
              !call InitA1Blocks(AC0Block, AuxData, A1vvao, A0vv, A0oa)
              !call InitA1Blocks(AC0Block, AuxData, A1vvaa, A0vv, A0aa)
              call InitA1Blocks(AC0Block, AuxData, A1vaoo, A0va, A0oo)
              call InitA1Blocks(AC0Block, AuxData, A1vaao, A0va, A0oa)
              call InitA1Blocks(AC0Block, AuxData, A1vaaa, A0va, A0aa)
              call InitA1Blocks(AC0Block, AuxData, A1aaoo, A0aa, A0oo)
              call InitA1Blocks(AC0Block, AuxData, A1aaao, A0aa, A0oa)
              !call InitA1Blocks(AC0Block, AuxData, A1vaao2, A0va, A0oa)


            end associate

      end subroutine init_AC0Blocks


      subroutine InitA1Blocks(AC0Block, AuxData, A1ind, A0ind1, A0ind2)
            type(TAC0Block), dimension(:), intent(inout) :: AC0Block
            type(TACppData), intent(in) :: AuxData
            integer, intent(in) :: A1ind, A0ind1, A0ind2
            integer :: NDim1S, NDim2S, NDim1T, NDim2T
            integer :: NDim1St, NDim1Tt, NDim2St, Ndim2Tt
            integer, dimension(:,:), allocatable :: indN1, indN2
            integer :: ind, i, j

            associate(NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)


              NDim1St = AC0Block(A0ind1)%NDimS
              NDim1Tt = AC0Block(A0ind1)%NDimT


              NDim2St = AC0Block(A0ind2)%NDimS
              NDim2Tt = AC0Block(A0ind2)%NDimT


              AC0Block(A1ind)%NDim1S = NDim1St
              AC0Block(A1ind)%NDim2S = NDim2St

              AC0Block(A1ind)%NDim1T = NDim1Tt
              AC0Block(A1ind)%NDim2T = NDim2Tt

              allocate(AC0Block(A1ind)%IndN1S(2, NDim1St))
              allocate(AC0Block(A1ind)%IndN2S(2, NDim2St))
              allocate(AC0Block(A1ind)%IndN1T(2, NDim1Tt))
              allocate(AC0Block(A1ind)%IndN2T(2, NDim1Tt))
              call CalcMemI(AC0Block(A1ind)%IndN1S, 'AC0Block(A1ind)%IndN1S')
              call CalcMemI(AC0Block(A1ind)%IndN2S, 'AC0Block(A1ind)%IndN2S')
              call CalcMemI(AC0Block(A1ind)%IndN1T, 'AC0Block(A1ind)%IndN1T')
              call CalcMemI(AC0Block(A1ind)%IndN2T, 'AC0Block(A1ind)%IndN2T')
              

              AC0Block(A1ind)%IndN1S = AC0Block(A0ind1)%IndNS
              AC0Block(A1ind)%IndN1T = AC0Block(A0ind1)%IndNT

              AC0Block(A1ind)%IndN2S = AC0Block(A0ind2)%IndNS
              AC0Block(A1ind)%IndN2T = AC0Block(A0ind2)%IndNT

              AC0Block(A1ind)%type = A1ind

              allocate(AC0Block(A1ind)%ASing(NDim1St, NDim2St))
              allocate(AC0Block(A1ind)%ATrip(NDim1Tt, NDim2Tt))
              allocate(AC0Block(A1ind)%ATripA(NDim1Tt, NDim2Tt))

              call CalcMem(AC0Block(A1ind)%ASing, 'AC0Block(A1ind)%ASing')
              call CalcMem(AC0Block(A1ind)%ATrip, 'AC0Block(A1ind)%ATrip')

              AC0Block(A1ind)%ASing = zero
              AC0Block(A1ind)%ATrip = zero
              AC0Block(A1ind)%ATripA = zero

            end associate

      end subroutine InitA1Blocks


      subroutine InitA0Blocks(A0Block, AuxData, IndAux, IndN, NDim, posS, posT, Ip, Iq)
            type(TAC0Block), intent(inout) :: A0Block
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:), intent(in) :: IndAux
            integer, dimension(:,:), intent(in) :: IndN
            integer, intent(in) :: IP, Iq
            integer, dimension(:,:), intent(inout) :: posS, posT
            integer, intent(in) :: NDim
            integer, parameter :: occ = 0, act = 1, virt = 2

            integer :: inds, indt
            integer :: i, j, p, q, i0, i1, j0, j1, Nbl, Nbl0, Nbl1
            integer :: indas, indat

            associate(NA=>AuxData%NA, NV=>AuxData%NV, NI=>AuxData%NI, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis)

              if (Ip == occ .and. Iq == occ)then
                    A0Block%type = A0oo
                    i0 = 1
                    i1 = NI
                    j0 = 1
                    j1 = NI
              else if (Ip == virt .and. Iq == virt) then
                    A0Block%type = A0vv
                    i0 = NIA+1
                    i1 = NBasis
                    j0 = NIA+1
                    j1 = NBasis
              else if (Ip == act .and. Iq == act) then
                    A0Block%type = A0aa
                    i0 = NI+1
                    i1 = NIA
                    j0 = NI+1
                    j1 = NIA
              else if (Ip == virt .and. Iq == act)then
                    A0Block%type = A0va
                    i0 = NIA + 1
                    i1 = NBasis
                    j0 = NI+1
                    j1 = NIA
                    Nbl = NV
                    Nbl0 = NIA+1
                    Nbl1 = NBasis
                    allocate(A0Block%MiniBlocks(Nbl))
                    do i = 1, Nbl
                          allocate(A0Block%MiniBlocks(i)%ActListS(NA))
                          allocate(A0Block%MiniBlocks(i)%ActListT(NA))
                          allocate(A0Block%MiniBlocks(i)%ActListTA(NA))
                    end do
              else if (Ip == occ .and. Iq == act)then
                    A0Block%type = A0oa
                    ! These are limits only for dim calculation.
                    ! In final version pairs are ordered OA
                    i0 = NI+1
                    i1 = NIA
                    j0 = 1
                    j1 = NI
                    Nbl = NI
                    Nbl0 = 1
                    Nbl1 = NI
                    allocate(A0Block%MiniBlocks(Nbl))
                    do i = 1, Nbl
                          allocate(A0Block%MiniBlocks(i)%ActListS(NA))
                          allocate(A0Block%MiniBlocks(i)%ActListT(NA))
                          allocate(A0Block%MiniBlocks(i)%ActListTA(NA))
                    end do
              end if



              if (A0Block%type .ne. A0oa .and. A0Block%type .ne. A0va)then

!                    print*, 'block not oa not va'
                    inds = 0
                    indt = 0
                    do i = i0, i1
                          do j = j0, min(j1, i)
                                ! print*, AuxData%Occ(i), AuxData%Occ(j)
                                ! print*, AuxData%Occ(i)+AuxData%Occ(j)-one
                                ! print*, AuxData%ThrPP
                                if ((abs((AuxData%Occ(i)+AuxData%Occ(j)-one)) > AuxData%ThrPP))then
                                      inds = inds + 1
                                      if (i.ne.j)then
                                            indt = indt + 1
                                      end if
                                end if
                          end do
                    end do
                    A0Block%NdimS = indS
                    A0Block%NdimT = indT
                    
                    allocate(A0Block%IndNS(2,indS))
                    allocate(A0Block%IndNT(2,indT))
                    
                    A0Block%IndNS = 0
                    A0Block%IndNT = 0

                    ! fill IndN for this Block
                        
                    inds = 0
                    indt = 0
                    
                    do i = i0, i1
                          do j = j0, min(j1, i)
                                if ((abs((AuxData%Occ(i)+AuxData%Occ(j)-one)) > AuxData%ThrPP))then
                                      inds = inds + 1
                                      A0Block%IndNS(1,inds) = i
                                      A0Block%IndNS(2,inds) = j
                                      posS(i, j) = inds
                                      if (i.ne.j)then
                                            indt = indt + 1
                                            A0Block%IndNT(1,indt) = i
                                            A0Block%IndNT(2,indt) = j
                                            posT(i, j) = indt
                                      end if
                                end if
                          end do
                    end do
              else

                    
                    A0Block%NdimS = NBl * NA
                    A0Block%NdimT = NBl * NA
                    allocate(A0Block%IndNS(2,NBl*NA))
                    allocate(A0Block%IndNT(2,NBl*NA))

                    A0Block%IndNS = 0
                    A0Block%IndNT = 0

                    inds = 0
                    indt = 0

                    do i = Nbl0, Nbl1
                          indas = 0
                          indat = 0
                          do j = NI+1, NIA
                                inds = inds + 1

                                if(A0Block%type == A0oa)then
                                      A0Block%IndNS(1,inds) = j
                                      A0Block%IndNS(2,inds) = i
                                      posS(j, i) = inds
                                else
                                      A0Block%IndNS(1,inds) = i
                                      A0Block%IndNS(2,inds) = j
                                      posS(i, j) = inds
                                end if

                                      
                                if ((abs((AuxData%Occ(i)+AuxData%Occ(j)-one)) > AuxData%ThrPP))then
                                      indas = indas + 1
                                      ! zapisz ktore aktywne sa w minibloczku
                                      A0Block%MiniBlocks(i-Nbl0+1)%ActListS(indas) = j
                                      A0Block%MiniBlocks(i-Nbl0+1)%NdimS = indas
                                end if
                                
                                if (i.ne.j)then
                                      
                                      indt = indt + 1
                                      if(A0Block%type == A0oa)then
                                            A0Block%IndNT(1,indt) = j
                                            A0Block%IndNT(2,indt) = i
                                            posT(j, i) = indt
                                      else
                                            A0Block%IndNT(1,indt) = i
                                            A0Block%IndNT(2,indt) = j
                                            posT(i, j) = indt
                                      end if
                                      
                                          if ((abs((AuxData%Occ(i)+AuxData%Occ(j)-one)) > AuxData%ThrPP))then
                                                indat = indat + 1
                                                A0Block%MiniBlocks(i-Nbl0+1)%ActListT(indat) = j
                                                A0Block%MiniBlocks(i-Nbl0+1)%NdimS = indat
                                          end if
                                    end if
                              end do
                        end do
                  end if
                  
            
            if (A0Block%type == A0oa .or. A0Block%type == A0va)then
                  
                  do i = 1, Nbl
                        allocate(A0Block%MiniBlocks(i)%MiniAVS(NA, NA))
                        allocate(A0Block%MiniBlocks(i)%MiniAVT(NA, NA))
                        allocate(A0Block%MiniBlocks(i)%MiniAVTA(NA, NA))
                        A0Block%MiniBlocks(i)%MiniAVS = zero
                        A0Block%MiniBlocks(i)%MiniAVT = zero
                        A0Block%MiniBlocks(i)%MiniAVTA = zero
                        allocate(A0Block%MiniBlocks(i)%EigS(NA))
                        allocate(A0Block%MiniBlocks(i)%EigT(NA))
                        allocate(A0Block%MiniBlocks(i)%EigTA(NA))
                  end do
                  
            else
                  if (A0Block%type == A0oo .or.A0Block%type == A0vv)then
                        !oo and vv blocks are diagonal. To save memory allocate only 1 column.
                        allocate(A0Block%ASing(A0Block%NdimS, 1))
                        allocate(A0Block%ATrip(A0Block%NdimT, 1))
                        allocate(A0Block%ATripA(A0Block%NdimT, 1))
                        A0Block%ASing = zero
                        A0Block%ATrip = zero
                        A0Block%ATripA = zero

                  else
                        allocate(A0Block%ASing(A0Block%NdimS, A0Block%NdimS))
                        allocate(A0Block%ATrip(A0Block%NdimT, A0Block%NdimT))
                        allocate(A0Block%ATripA(A0Block%NdimT, A0Block%NdimT))
                        A0Block%ASing = zero
                        A0Block%ATrip = zero
                        A0Block%ATripA = zero

                  end if

                  allocate(A0Block%EigS(A0Block%NdimS))
                  allocate(A0Block%EigT(A0Block%NdimT))
                  allocate(A0Block%EigTA(A0Block%NdimT))

                  A0Block%EigS = zero
                  A0Block%EigT = zero
                  A0Block%EigTA = zero
                  

                  allocate(A0Block%vPlusS(A0Block%NDimS))
                  allocate(A0Block%vPlusT(A0Block%NDimT))
                  allocate(A0Block%vPlusTA(A0Block%NDimT))
                  A0Block%vPlusS = zero
                  A0Block%vPlusT = zero
                  A0Block%vPlusTA = zero
            end if
            
          end associate
      end subroutine InitA0Blocks

      ! subroutine MapMiniBlocks(AC0Block, AuxData, Nbl)
      !       type(TAC0Block), intent(inout) :: AC0Block
      !       type(TACppData), intent(in) :: AuxData
      !       integer ::  NVused
      !       integer :: i, offset, k
      !       integer :: ThisVirt, ThisSize
            

      !       associate (NV=> AuxData%NV, NIA=>AuxData%NIA, NI=>AuxData%NI)
      !       allocate(AC0Block%MiniMap(2, Nbl))

      !       AC0Block%MiniMap = 0

      !       offset = NIA
      !       NVused = 0
      !       ThisVirt = AC0Block%IndNS(1,1)
      !       ThisSize = 0
      !       AC0Block%MiniBlock(1, 1) = 1

      !       k = 1
      !       do i = 1, AC0Block%NDimS
      !             if (AC0Block%IndNS(1, i) == ThisVirt)then
      !                   ThisSize = ThisSize + 1
      !             else
      !                   AC0Block%MiniBlock(2, k) = ThisSize
      !                   k = k + 1
      !                   ThisSize = 1
      !                   ThisVirt = AC0Block%IndNS(1, i)
      !                   if (k<=NV)then
      !                         AC0Block%MiniBlock(1, k) = i
      !                   end if
      !             end if
      !       end do
      !       AC0Block%MiniBlock(2, k) = ThisSize


      !       ! do k = 1, NV
      !       !       print*, 'idx virt', k+NIA
      !       !       print*, 'start idx virt', AC0Block%MiniBlock(1,k)
      !       !       print*, 'wymiar idx virt', AC0Block%MiniBlock(2,k)
      !       !       print*, ''
      !       ! end do

      !     end associate
      ! end subroutine InitVABlock


      subroutine HNO0_init(ACAlpha, XOne, HNO0, NBasis, IndAux)

            double precision, intent(in) :: ACAlpha
            double precision, dimension(:), intent(in) :: XOne
            double precision, dimension(:,:), intent(inout) :: HNO0
            integer, dimension(:), intent(in) :: IndAux
            integer, intent(in) :: NBasis


            double precision, dimension(:,:), allocatable :: Ure
            double precision, dimension(:), allocatable :: work1, work2
            integer :: i, j

            allocate(Ure(NBasis, NBasis))
            allocate(work1(NBasis**2),work2(NBasis**2))

            URe = Zero
            do i = 1, Nbasis
                  URe(i,i) = One
            end do


            call triang_to_sq(XOne,work1,NBasis)
            call dgemm('N','N',NBasis,NBasis,NBasis,One,URE,NBasis,work1,NBasis,zero,work2,NBasis)
            call dgemm('N','T',NBasis,NBasis,NBasis,One,work2,NBasis,URE,NBasis,zero,HNO0,NBasis)
            call sq_symmetrize(HNO0,NBasis)
            deallocate(ure)
            deallocate(work1)
            deallocate(work2)


      end subroutine HNO0_init


      subroutine THC_init(Flags, AuxData, THCData, HNO0_EXT, HNO0_THC, CAONO_IN, RdmData)

            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use geom_input, only: geom_ReadSystemBasis
            use sys_definitions
            use gammcor_integrals

            
            type(FlagsData), intent(in) :: Flags
            type(TACppData) :: AuxData
            type(TTHCData), intent(inout) :: THCData
            type(TRdmData) :: RdmData
            !integer, intent(out) :: NChol, NTHC
            !double precision, allocatable, intent(out) :: Xga(:,:), Zgk(:,:)
            double precision, dimension(:,:), intent(in) :: HNO0_EXT
            double precision, dimension(:,:), intent(out) :: HNO0_THC
            double precision, dimension(:,:), intent(inout) :: CAONO_IN
            
            double precision, dimension(:,:), allocatable ::  CAONO


            character(:), allocatable :: basis_path, xyz_path, binPath
            logical, parameter :: sort_ang_mom = .true.
            type(TAOBASIS) :: AObasis
            type(TSystem) :: System

            double precision :: CholThr, THCThr
            integer :: CholAccu
            integer :: i, j, NAO, ij, ab
            double precision, allocatable :: Xgp(:,:)
            integer :: NA, NI, NV, NIA
            integer :: units, nbasis
            type (tclock) :: timer

            ! Fock matrix diag
            double precision, dimension(:,:), allocatable :: Fockij, Fockvw
            double precision, dimension(:,:), allocatable :: Cpi_extao, Cpa_extao, Cpv_extao
            double precision, dimension(:,:), allocatable :: work, H0_extao

            !---------------------testing energy
            integer :: p, q, r, s, k, l, a, b, ii, nn
            double precision, allocatable :: Rkab(:,:,:), Rkcd(:,:,:)
            double precision :: this, ETot, etot0, val, val_this


            
            NA = AuxData%NA
            NI = AuxData%NI
            NV = AuxData%NV
            NIA = NI+NA

            units = Flags%IUnits

            basis_path = Flags%BasisSetPath //Flags%BasisSet
            nbasis = AuxData%NBasis
            print*, 'units', units
            
            print*, 'nbasis', nbasis
            print*, 'basis_path', basis_path
            xyz_path = "./input.inp"
            
            call clock_start(timer)
            ! basis/geometry setup goes through geom_input (see geom_input.f90 header);
            ! basis_newAObasis is no longer part of the gammcor-integrals API.
            call geom_ReadSystemBasis(System, AObasis, Flags, sort_ang_mom)

            nao = AObasis%NAOSpher
            print*, 'nao', nao
            call clock_start(timer)

            CholThr = Flags%DCholeskyThr
            THCThr = Flags%DTHCthr
            CholAccu = Flags%ICholeskyAccu

            if (CholThr < zero.or.THCThr <zero)then
                  print*, 'CholAccu', CholAccu
                  call thc_gammcor_XZ(Xgp, THCData%Zgk, AOBasis, System, CholAccu)
            else
                  print*, 'CholeskyThreshold', CholThr
                  print*, 'THCThreshold', THCThr
                  call thc_gammcor_XZ(Xgp, THCData%Zgk, AOBasis, System,CholAccu, CholThr, THCThr)
                  print*, 'CholeskyThreshold', CholThr
                  print*, 'THCThreshold', THCThr
            end if


            allocate(CAONO(nbasis, nbasis))
            allocate(work(nao, nbasis))
            
            call CalcMem(CAONO, 'CAONO')
            call CalcMem(work, 'work')

            if (AuxData%PYSCF==1 .or. AuxData%ORCA==1)then
                  print*, 'here', AuxData%NI, AuxData%NIA, AuxData%NBasis
                  CAONO = CAONO_IN
                  print*, CAONO(1,2), CAONO(2, 1)
            else
                  CAONO = transpose(CAONO_IN)
                  print*, CAONO(1,2), CAONO(2, 1)
            end if

            call canonicalize(CAONO, THCData%fij, THCData%fvw, AuxData, THCData, Xgp, AObasis, System)


            allocate(H0_extao(nao, nao))
            call CalcMem(H0_extao, 'h0_extao')
            call ints1e_gammcor_H0_extao(H0_extao, AObasis, System, THCData%ExternalOrdering)
            print*, 'h01', H0_extao(1,1)

            call real_ab(work, H0_extao, CAONO)
            call real_atb(HNO0_THC, CAONO, work)

 

            !------------------------------------------testing energy--------------------
            THCData%NTHC=size(Xgp,dim=1)
            THCData%NChol=size(THCData%Zgk,dim=2)
            allocate(THCData%Xga(THCData%NTHC,NBasis))
            call CalcMem(THCData%Xga, 'THCData%Xga')

            Call thc_gammcor_Xga(THCData%Xga, Xgp, CAONO,&
                  AOBasis, THCData%ExternalOrdering)

            

            allocate(Rkab(THCData%NChol,nao, nao))
            allocate(Rkcd(THCData%NChol, nao, nao))
            Call thc_gammcor_Rkab_2(Rkab, THCData%Xga, THCData%Xga, THCData%Zgk, NBasis, NBasis,&
                  THCData%NChol, THCData%NTHC)
            Call thc_gammcor_Rkab_2(Rkcd, THCData%Xga, THCData%Xga, THCData%Zgk, NBasis, NBasis,&
                  THCData%NChol, THCData%NTHC)
            
            etot0 = zero
            do i = 1, NBasis                  
                  ETot0 = ETot0 + two* AuxData%Occ(i) * HNO0_EXT(i,i)!AuxData%HNO0(i,i)
                  write(*, '(I5, 2F20.15)') i, HNO0_EXT(i,i), AuxData%Occ(i) 
            end do

            print*, 'etot1 z HNO_ext', ETot0

            ETot0 = zero
            do i = 1, NBasis                  
                  ETot0 = ETot0 + two* AuxData%Occ(i) * HNO0_THC(i,i)
            end do

            print*, 'etot1 z HNO thc', ETot0


            associate(Occ=>AuxData%Occ, IndAux=>AuxData%IndAux, map=>AuxData%map)
              etot = zero
            do p = 1, NI+NA
                  do q = 1, NI+NA
                        do r = 1, NI+NA
                              do s = 1, NI+NA
                                    call real_vw_x(this, Rkab(:, p, r), Rkcd(:, q, s), THCData%NChol)

                                    
                                    val = zero
                                    if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                          val = val+  (RdmData%rdm2_pp_act(map(p), map(q), map(r), map(s)) + RdmData%rdm2_pm_act(map(p), map(q), map(r), map(s)))
                                    else
                                          if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then                                                
                                                val =  val + two * Occ(p) * Occ(q)
                                          end if

                                          if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                val = val +  -Occ(p) * Occ(q)
                                          end if
                                    end if
                                    etot = etot + val * this
                              end do
                        end do
                  end do
            end do
            print*, 'etot0', etot0
            print*, 'etot1', etot
            print*, 'etot w/o enuc', etot+etot0
            print*, 'etot', etot+etot0+AuxData%enuc
            print*, 'ecas', AuxData%Ecas
            !print*, 'etot', etot0, etot+AuxData%enuc,  etotetot, AuxData%Ecas
          end associate
         stop


            
            ! hno0_THC = zero

            ! ! transformigo 1-el H0 THC to NO basis
            ! call real_ab(work, H0_extao, CAONO)
            ! call real_ab(HNO0_THC, transpose(CAONO), work)

            ! print*, 'H0_extao', H0_extao(1,1), H0_extao(2,2)
            ! print*, 'caono', CAONO(1,1), CAONO(1,2), CAONO(1,3)
            ! print*, 'HNO0_tHC', HNO0_THC(1,1), HNO0_THC(2,2)
            ! print*, 'HNO0_ext', HNO0_EXT(1,1), HNO0_EXT(2,2)
!            stop

            ! do i = 1, nbasis
            !       do j = 1, i
            !             if (abs(abs(HNO0_THC(i, j))-abs(HNO0_EXT(i, j))).gt.1.d-3)then
            !                   write(*, '(2I5, 3F20.15)')i, j, HNO0_THC(i, j),HNO0_EXT(i, j)
            !                   print*, 'FAILED H0 test'
            !             !      stop
            !              end if
            !        end do
            !  end do


            

            ! allocate(Fockij(AuxData%NI, AuxData%NI))
            ! allocate(Fockvw(AuxData%NV, AuxData%NV))
            ! allocate(Cpi_extao(AuxData%nbasis, AuxData%NI))
            ! allocate(Cpv_extao(AuxData%nbasis, AuxData%NV))


            ! print*, 'kawalek', 1, AuxData%NI
            ! print*, 'kawalek', AuxData%NI+1, AuxData%NIA
            ! print*, 'kawalek', AuxData%NIA+1, AuxData%NBasis
            ! call thc_gammcor_F(Fockij, Fockvw, CAONO(:, 1:AuxData%NI),&
            !       CAONO(:, AuxData%NI+1:AuxData%NIA), &
            !       CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
            !       AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering)
            
            ! call symmetric_eigenproblem(fij, Fockij, NI, .true.)
            ! call symmetric_eigenproblem(fvw, Fockvw, NV, .true.)


            ! call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
            ! call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)

            ! CAONO(:, 1:AuxData%NI)=Cpi_extao
            ! CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao


            ! print*, 'time for THC XZ', clock_readwall(timer)
            
            ! THCData%NTHC=size(Xgp,dim=1)
            ! NChol=size(Zgk,dim=2)
            ! allocate(Xga(THCData%NTHC,NBasis))

            
            ! call clock_start(timer)
            ! Call thc_gammcor_Xga(Xga, Xgp, CAONO,&
            !       AOBasis, ExternalOrdering)
            ! print*, 'Time for AOMO transformation', clock_readwall(timer)


            ! if (nao.ne.nbasis)then
            !       print*, 'INCOMPATIBLE NAO and NMO, DALTON INTERFACE ASSUMES THEY ARE EQUAL, ERROR, EXITING'
            !       stop
            ! end if
            
            ! Transformin 1-el to diagonalized Fock basis
            ! HNO0_THC = zero
            ! call real_ab(work, H0_extao, CAONO)
            ! call real_atb(HNO0_THC, CAONO, work)

            ! ETot0 = zero
            ! do i = 1, AuxData%NIA
            !       ETot0 = ETot0 + two* AuxData%Occ(i) * HNO0_THC(i,i)!AuxData%HNO0(i,i)
            !       !if (abs(two* AuxData%Occ(i) * HNO0_EXT(i,i)).gt.1.d-5)then
            !             write(*, '(I5, 2F20.15)') i, HNO0_THC(i,i), AuxData%Occ(i) 
            !       !end if
            ! end do
            ! print*, etot0
            ! stop



      end subroutine THC_init

      ! subroutine canonicalize(CAONO, fij, fvw, AuxData, THCData, Xgp, AObasis, System)
      !       use Cholesky_Gammcor
      !       use THC_Gammcor
      !       use OneElectronInts_Gammcor
      !       use basis_sets
      !       use sys_definitions
      !       use gammcor_integrals

      !       double precision, dimension(:,:), intent(inout) :: CAONO
      !       double precision, dimension(:), intent(out) :: fij, fvw            
      !       type(TACppData), intent(in) :: AuxData
      !       type(TTHCData), intent(in) :: THCData
      !       type(TAOBASIS) :: AObasis
      !       type(TSystem) :: System


      !       double precision, dimension(:,:), allocatable :: Cpi_extao, Cpv_extao
      !       double precision, dimension(:,:), allocatable :: Fockij, Fockvw
      !       double precision, dimension(:,:), intent(in) :: Xgp

      !       associate(Zgk=>THCData%Zgk, ExternalOrdering=>THCData%ExternalOrdering)

      !         print*, 'order2', THCData%ExternalOrdering
      !         print*, 'order3', ExternalOrdering
      !       allocate(Fockij(AuxData%NI, AuxData%NI))
      !       allocate(Fockvw(AuxData%NV, AuxData%NV))
      !       allocate(Cpi_extao(AuxData%nbasis, AuxData%NI))
      !       allocate(Cpv_extao(AuxData%nbasis, AuxData%NV))

      !       call CalcMem(Fockij, 'Fockij')
      !       call CalcMem(Fockvw, 'Fockvw')
      !       call CalcMem(Cpi_extao, 'Cpi_extao')
      !       call CalcMem(Cpv_extao, 'Cpv_extao')



      !       call thc_gammcor_F(Fockij, Fockvw, CAONO(:, 1:AuxData%NI),&
      !             CAONO(:, AuxData%NI+1:AuxData%NIA), &
      !             CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
      !             AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering)

      !       call symmetric_eigenproblem(fij, Fockij, AuxData%NI, .true.)
      !       call symmetric_eigenproblem(fvw, Fockvw, AuxData%NV, .true.)


      !       call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
      !       call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)

      !       CAONO(:, 1:AuxData%NI)=Cpi_extao
      !       CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao
      !     end associate
      ! end subroutine canonicalize

      subroutine THC_int_loop(THCData, ACB, AuxII, AuxAA, BuxII, BuxAA, Aux3X, Aux3B, &
            Flags, AuxData, Occ, map, RdmData, posS, posT, AcAlpha)
            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use sys_definitions
            use gammcor_integrals



            type(TTHCData), intent(in) :: THCData
 
            type(TAC0Block), dimension(:), intent(inout) :: ACB
            double precision, dimension(:,:), intent(inout) :: AuxII, AuxAA, BuxII, BuxAA
            double precision, dimension(:,:), intent(inout) :: Aux3X, Aux3B


            double precision, dimension(:), intent(in) :: Occ
            character(:), allocatable :: basis_path, xyz_path, binPath
            type(TAOBASIS) :: AObasis
            type(TSystem) :: System
            type(FlagsData), intent(in) :: Flags
            type(TACppData) :: AuxData
            integer, dimension(:), intent(in) :: map
            type(TRdmData) :: RdmData
            double precision :: Acalpha


            double precision, dimension(:,:), allocatable ::ASingChunk, ATripChunk
            !----THC---                                                                                                                                                                                                                                                                                                                                                                      
            ! integer, parameter :: ExternalOrdering = ORBITAL_ORDERING_DALTON
            ! logical, parameter :: sort_ang_mom = .true.
            ! double precision :: CholThr, THCThr
            ! integer :: CholAccu
            integer :: i, j, NAO
            ! integer :: NChol, NTHC
            ! double precision, allocatable :: Xgp(:,:), Zgk(:,:), Xga(:,:)
            double precision, allocatable :: XXga(:,:)
            double precision, allocatable :: Rkax(:,:), Rkby(:,:)
            integer, dimension(6) :: TabChol


            integer, dimension(:), allocatable :: indaux, ind2
            integer, dimension(:,:), intent(in) :: posS, posT
            integer :: units
            integer, dimension(2, 19) :: integral_dummy
            type (tclock) :: timer, timer0, timer90, timerx
            double precision :: timer1, timer2, timer3, timer4, timer5, timer6, timer7, timer8
            double precision :: timer9, timer10, timer11, timer12, timer13, timer14, timer15, timer16, timer17
            double precision :: timer18, timer19, timer20, timer21, timer22, timer23, timer24, timer25, timer26, timer27
            double precision :: timer28, timer29, timer30, timer31, timer32, timer33, timer34, &
                  timer35, timer36, timer37, timer38, timer39, timer40, timer41, &
                  timer42, timer43, timer44, timer45, timer46, timer47, timer48, &
                  timer49, timer50, timer51, timer52, timer53, timer54, timer55, &
                  timer74, timer75, timer76, timer70, timer71, timer72, timer73, timer100, &                  
                  timer56, timer57, timer58, timer59, timer101, timer102, timer103, timer104, timer203
            double precision :: timer60, timer61, timer62, timer63, timer64, timer65, timer66, &
                  timer67, timer68, timer69
            double precision, dimension(:), allocatable :: V_axby
            integer :: NA, NI, NV, NIA, nbasis
            integer ::  p, q, r, s, pq1, rs1, pq2, rs2
            integer :: a, b, c, d
            integer :: x0, x1, y0, y1
            integer :: c0, c1, d0, d1
            integer :: t, u
            integer :: ii
            double precision :: val, val1, val2
            type(TLoop), dimension(6) :: loops

            integer :: b0II, b1II, b0IA, b1IA, b0AA, b1AA
            integer :: b0IV, b1IV, b0VA, b1VA, b0VV, b1VV
            integer :: BatchDim
            integer :: NBatchI, NBatchA, NBatchV
            integer :: batch
            double precision, dimension(:,:,:), allocatable :: RII, RIA, RAA, RIV, RVA, RVV

            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            double precision :: timer104_wall, timer104_cpu
            double precision, dimension(:,:,:,:), allocatable :: rdm2_sum, rdm2_sum2
            double precision, dimension(:,:,:,:), allocatable :: rdm2_nsum, rdm2_nsum2
            double precision :: start_wall, end_wall
            double precision :: local_timer104_cpu, sttime, endtime
            integer :: pq, x, y
            integer, dimension(:, :), allocatable :: ind_virt, ind_va
            integer :: ndim_virt, k, nva, kk
            integer, allocatable :: ACBlockList(:,:,:)
            integer :: nblock, bi
            integer :: method
            integer, parameter :: AC0 = 0, AC= 1



              
            BatchDim = AuxData%Batchdim
            if (BatchDim==0)then
                  print*, 'Request BatchDim at least 1'
                  stop
            end if

            call ximsg('Batch Dimension for Cholesky vectors',   BatchDim, mnormal)

            associate (rdm2_pp_13_act => RdmData%rdm2_pp_13_act(:,:,:,:), &
                  rdm2_pm_13_act => RdmData%rdm2_pm_13_act(:,:,:,:), &
                  rdm2_pp_12_act => RdmData%rdm2_pp_12_act(:,:,:,:), &
                  rdm2_pm_12_act => RdmData%rdm2_pm_12_act(:,:,:,:), &
                  rdm2_pp_act => RdmData%rdm2_pp_act(:,:,:,:), &
                  rdm2_pm_act => RdmData%rdm2_pm_act(:,:,:,:))
              NA = AuxData%NA
              NI = AuxData%NI
              NV = AuxData%NV
              NIA = NI+NA
              NBasis = AuxData%NBasis

              allocate(ind_virt(2, (NV*(NV+1))/2))
              allocate(ind_va(2, (NV)*NA))
              k = 1
              kk = 1
              do ii = 1, AuxData%NDim_s
                    if (AuxData%IndAux(AuxData%IndN_s(1, ii))==2.and. AuxData%IndAux(AuxData%IndN_s(2,ii))==2)then
                          ind_virt(1, k) =  AuxData%IndN_s(1, ii)
                          ind_virt(2, k) =  AuxData%IndN_s(2, ii)
                          !print*, 'virt',  ind_act(1, k),  ind_act(2, k), k
                          k = k+ 1
                    end if
                    if ((AuxData%IndN_s(1, ii))==2.and. AuxData%IndAux(AuxData%IndN_s(2,ii))==1)then
                          ind_va(1, kk) =  AuxData%IndN_s(1, ii)
                          ind_va(2, kk) =  AuxData%IndN_s(2, ii)
                          kk = kk+ 1
                    end if

              end do
              ndim_virt = k -1
              nva = kk - 1

              associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga, ExternalOrdering=>THCData%ExternalOrdering, NChol=>THCData%NChol, NTHC=>THCData%NTHC)
              allocate(Rkax(THCData%NChol, NBasis))
              allocate(Rkby(NChol, NBasis))
              allocate(XXga(NTHC,NBasis))
              allocate(rdm2_sum(NA, NA, NA, NA))
              allocate(rdm2_sum2(NA, NA, NA, NA))

              allocate(rdm2_nsum(NA, NA, NA, NA))
              allocate(rdm2_nsum2(NA, NA, NA, NA))

              rdm2_sum =  rdm2_pm_13_act + rdm2_pp_13_act
              rdm2_sum2 = rdm2_pm_act + rdm2_pp_act

              rdm2_nsum =  rdm2_pm_act  + rdm2_pp_act
              rdm2_nsum2 = rdm2_pm_13_act + rdm2_pp_13_act
              allocate(V_axby(max(NA, NI, NV)**2))

              if (Flags%Jobtype==JOB_TYPE_AC0PP)then
                    method = 0
                    nblock = 6
                    allocate(ACBlockList(6,6,6))
                    call parse_config(loops, integral_dummy, TabChol, ACBlockList, method)
              else if (Flags%Jobtype==JOB_TYPE_ACPP)then
                    method = 1
                    nblock = 1
                    allocate(ACBlockList(6,6, 1))
                    call parse_config(loops, integral_dummy, TabChol, ACBlockList, method)
              end if
              call clock_start(timer0)

              
              timer1 = zero
              timer2 = zero
              timer3 = zero
              timer4 = zero
              timer5 = zero
              timer6 = zero
              timer7 = zero
              timer8 = zero
              timer9 = zero
              timer10 = zero
              timer11 = zero
              timer12 = zero
              timer13 = zero
              timer14 = zero
              timer15 = zero
              timer16 = zero
              timer17 = zero
              timer18 = zero
              timer19 = zero
              timer20 = zero
              timer21 = zero
              timer22 = zero
              timer23 = zero
              timer24 = zero
              timer25 = zero
              timer26 = zero
              timer27 = zero
              timer28 = zero
              timer29 = zero
              timer30 = zero
              timer31 = zero
              timer32 = zero
              timer33 = zero
              timer34 = zero
              timer35 = zero
              timer36 = zero
              timer37 = zero
              timer38 = zero
              timer39 = zero
              timer40 = zero
              timer41 = zero
              timer42 = zero
              timer43 = zero
              timer44 = zero
              timer45 = zero
              timer46 = zero
              timer47 = zero
              timer48 = zero
              timer49 = zero
              timer50 = zero
              timer51 = zero
              timer52 = zero
              timer53 = zero
              timer54 = zero
              timer55 = zero
              timer56 = zero
              timer57 = zero
              timer58 = zero
              timer59 = zero
              timer60 = zero
              timer61 = zero
              timer62 = zero
              timer63 = zero
              timer64 = zero
              timer65 = zero
              timer66 = zero
              timer67 = zero
              timer68 = zero
              timer69 = zero
              timer70 = zero
              timer71 = zero
              timer72 = zero
              timer73 = zero
              timer74 = zero
              timer75 = zero
              timer76 = zero
              timer100 = zero
              timer101 = zero
              timer102 = zero
              timer103 = zero
              timer203 = zero
              timer104 = zero



              NBatchI = NI / BatchDim
              NBatchA = NA / BatchDim
              NBatchV = NV / BatchDim
         
              nblock = size(ACBlockList, dim=3)
         
              if (modulo(NI, BatchDim)>0)NBatchI = NBatchI + 1
              if (modulo(NA, BatchDim)>0)NBatchA = NBatchA + 1
              if (modulo(NV, BatchDim)>0)NBatchV = NBatchV + 1

              call print_info('NBatchI', NBatchI)
              call print_info('NBatchA', NBatchA)
              call print_info('NBatchV', NBatchV)

              call print_section('THC integral transformation')
         


              if (TabChol(1) == 1) then
                    allocate(RII(Nchol, NI, min(NI, BatchDim)))
                    call CalcMem3(RII, 'RII')
              end if
              if (TabChol(2) == 1) then
                    allocate(RIA(Nchol, NA, min(NI, BatchDim)))
                    call CalcMem3(RIA, 'RIA')
              end if
              if (TabChol(3) == 1)then
                    allocate(RAA(Nchol, NA, min(NA, BatchDim)))
                    call CalcMem3(RAA, 'RAA')
              end if
         
              if (TabChol(4) == 1)then
                    allocate(RIV(Nchol, NV, min(NI, BatchDim)))
                    call CalcMem3(RIV, 'RIV')
              end if
              if (TabChol(5) == 1) then
                    allocate(RVA(Nchol, NA, min(NV, BatchDim)))
                    call CalcMem3(RVA, 'RVA')
              end if
              if (TabChol(6) == 1) then
                    allocate(RVV(Nchol, NV, min(NV, BatchDim)))
                    call CalcMem3(RVV, 'RVV')
              end if
         

              call clock_start(timer)              
              if (TabChol(1) == 1) call thc_gammcor_Rkab_2(RII, Xga(:,1:NI), Xga(:,1:NI), Zgk, NI, NI, NChol, NTHC)
              call tmsg('TIME FOR RII', timer, tdebug)
              call clock_start(timer)              
              if (TabChol(2) == 1) call thc_gammcor_Rkab_2(RIA, Xga(:,NI+1:NIA), Xga(:,1:NI), Zgk, NA, NI, NChol, NTHC)
              call tmsg('TIME FOR RIA', timer, tdebug)
              call clock_start(timer)              

              if (TabChol(3) == 1) call thc_gammcor_Rkab_2(RAA, Xga(:,NI+1:NIA), Xga(:,NI+1:NIA), Zgk, NA, NA, NChol, NTHC)
              call tmsg('TIME FOR RAA', timer, tdebug)
              call clock_start(timer)              
              
              timer102 = timer102 + clock_readwall(timer)

              batchloop: do batch = 1, max(NBatchI, NBatchA, NBatchV)

                    b0IV = 1 + (batch-1) * BatchDim
                    b1IV = min(b0IV+BatchDim-1, NI)

                    ! b0AV = NI + 1 + (batch-1) * BatchDim
                    ! b1AA = min(b0AV+BatchDim-1, NIA)

                    b0VA = NIA + 1 + (batch-1) * BatchDim
                    b1VA = min(b0VA+BatchDim-1, NBasis)

                    b0VV = NIA + 1 + (batch-1) * BatchDim
                    b1VV = min(b0VV+BatchDim-1, NBasis)

                    call clock_start(timer)              
                    if (b0IV <=NI)then
                          if (TabChol(4)==1)then
                                do b = b0IV, b1IV
                                      call thc_gammcor_Rkab_Batch_a_Fixed_b(RIV(:, 1:NV, b-b0IV+1), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, b), Zgk, NV, NChol, NTHC)       
                                end do
                          end if
                    end if
                    if (TabChol(4)==1) call tmsg('TIME FOR RIV', timer, tdebug)
                    call clock_start(timer)
                    
                    if (b0VA <=NBasis)then
                          if (TabChol(5)==1)then

                                do b = b0VA, b1VA
                                      call thc_gammcor_Rkab_Batch_a_Fixed_b(RVA(:, 1:NA, b-b0VA+1), XXga(:, 1:NA), Xga(:, NI+1:NIA), Xga(:, b), Zgk, NA, NChol, NTHC)	
                                end do
                          end if
                    end if
                    if (TabChol(5)==1) call tmsg('TIME FOR RVA', timer, tdebug)
                    call clock_start(timer)

                    if (b0VV <=NBasis)then
                          if (TabChol(6)==1)then
                                do b = b0VV, b1VV
                                      call thc_gammcor_Rkab_Batch_a_Fixed_b(RVV(:, 1:NV, b-b0VV+1), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, b), Zgk, NV, NChol, NTHC)	
                                end do
                          end if
                    end if
                    if (TabChol(6)==1) call tmsg('TIME FOR RVV', timer, tdebug)
                    call clock_start(timer)

                    caseloop: do ii = 1, 6
                          if (loops(ii)%run_main)then
                                select case(ii)
                                case(loopII) 
                                      loopII: do a = 1, NI
                                            if (loops(ii)%run_inner(1)) then ! (II|II)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NI, RII(:,:,a), NChol, RII(:,:,b), NChol, NI, NI, NChol, AcAlpha)
                                                        timer3 = timer3 + clock_readwall(timer)

                                                        call clock_start(timer)


                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_oooo_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_oooo_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3a1_II_oooo_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a2_II_oooo_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a3_II_oooo_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a4_II_oooo_pr_qrsp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3b1_II_oooo_ps_qsrp(ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b2_II_oooo_qr_prsq(ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b3_II_oooo_qs_psrq(ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b4_II_oooo_pr_qrsp(ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3c1_II_oooo_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c3_II_oooo_pr_qspr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c4_II_oooo_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c5_II_oooo_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, 1, NI, 1, NI, posS, posT, Occ, V_axby)

                                                        end if

                                                        call p2_IIII_qstt_qstt(BuxII, a, b, 1, NI, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IIII_qtst_qtst(AuxII, a, b, 1, NI, 1, NI, NI, NIA, NBasis, Occ, V_axby)

                                                        timer4 = timer4 + clock_readwall(timer)
                                                  end do
                                            end if
                                      end do loopII
                                case(loopIA) ! loop2
                                      loopIA: do a = 1, NI
                                            if (loops(ii)%run_inner(1)) then ! (IA|II)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                        timer7 = timer7 + clock_readwall(timer)

                                                        call clock_start(timer)
                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_aooo_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call bare_int_prqs_ooao_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aooo_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_ooao_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3a1_II_aooo_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a1_AI_ooao_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a2_II_ooao_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a2_IA_aooo_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a3_II_aooo_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a3_II_ooao_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a4_AI_ooao_pr_qrsp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a4_IA_aooo_pr_qrsp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3b1_II_aooo_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b1_AI_ooao_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b2_II_ooao_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b2_IA_aooo_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b3_II_aooo_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b3_II_ooao_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b4_AI_ooao_pr_qrsp(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b4_IA_aooo_pr_qrsp(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              
                                                              call p3c1_II_aooo_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c1_II_ooao_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c3_AI_ooao_pr_qspr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c3_IA_aooo_pr_qspr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c4_II_aooo_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c4_AI_ooao_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c5_II_ooao_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3c5_IA_aooo_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                        end if

                                                        call p2_IAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IAII_qstt_sqtt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IAII_qtst_tsqt(AuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IAII_qtst_tqst(AuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)

                                                        timer8 = timer8 + clock_readwall(timer)

                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(2)) then! (IA|IA)

                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                        timer10 = timer10 +  clock_readwall(timer)
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1aaoo
                                                        end if

                                                        call bare_int_prqs_aaoo_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call bare_int_psqr_aaoo_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer12 = timer12 +  clock_readwall(timer)
                                                        call clock_start(timer)

                                                        call p3a1_IA_aaoo_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3a2_IA_aaoo_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3a3_IA_aaoo_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3a4_IA_aaoo_pr_qrsp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        
                                                        call p3b1_IA_aaoo_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b2_IA_aaoo_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b3_IA_aaoo_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b4_IA_aaoo_pr_qrsp(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        call p3c1_IA_aaoo_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3c3_IA_aaoo_pr_qspr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3c4_IA_aaoo_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3c5_IA_aaoo_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer11 = timer11 +  clock_readwall(timer)

                                                        
                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_ooaa_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aoao_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_ooaa_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a1_AI_ooaa_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a2_AI_ooaa_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a3_II_aoao_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a3_AI_ooaa_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a4_AI_ooaa_pr_qrsp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3b1_AI_ooaa_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b2_AI_ooaa_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b3_II_aoao_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b3_AI_ooaa_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b4_AI_ooaa_pr_qrsp(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              
                                                              call p3c1_AI_ooaa_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c3_AI_ooaa_pr_qspr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c4_AI_aoao_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c4_AI_ooaa_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c5_AI_ooaa_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c5_IA_aoao_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              !---------------------P3a-4-P3b2---------------------------------------------------------------
                                                              !xxx - ao ao
                                                              !call P3a4_P3b2(a, b, NI+1, NIA, NI+1, NIA,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .false.)
              call P3ab42_1(a, b, NI+1, NIA, NI+1, NIA, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)

                                                        end if

                                                        call clock_start(timer)
                                                        call p2_IAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IAIA_qtst_tqts(AuxII, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        timer13 = timer13+  clock_readwall(timer)

                                                  end do
                                            end if
                                      end do loopIA


                                case(loopAA)! loop3
                                      loopAA: do a = NI+1, NIA
                                            if (loops(ii)%run_inner(1)) then ! (AA|II)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                        timer14 = timer14 +  clock_readwall(timer)

                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_aoao_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3a1_AI_aoao_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3a2_IA_aoao_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3b1_AI_aoao_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                              call p3b2_IA_aoao_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)

                                                              call p3c1_II_aoao_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, 1, NI, posS, posT, Occ, V_axby)
                                                        end if

                                                        call p2_AAII_qstt_ttqs(BuxAA, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)

                                                        timer14 = timer14 +  clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(2)) then! (AA|IA) for ddot-1, for ddot-2
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                        timer16 = timer16 +  clock_readwall(timer)

                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1aaao
                                                        end if

                                                        call bare_int_prqs_aaao_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call bare_int_psqr_aaao_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer18 = timer18 +  clock_readwall(timer)
                                                        call clock_start(timer)
                                                        call p3a2_IA_aaao_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3a3_IA_aaao_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        call p3b2_IA_aaao_qr_prsq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b3_IA_aaao_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        call p3c1_IA_aaao_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3c5_IA_aaao_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        timer17 = timer17 +  clock_readwall(timer)
                                                        call clock_start(timer)
                                                        !---------------------P3a-1,2-P3b3,4---------------------------------------------------------------                                                      
                                                        c0 = NI+1
                                                        c1 = min(NIA, a)
                                                        d0 = NI+1
                                                        d1 = NIA

                                                        !xxx A1 aaao do pq1 i rs1                                                      
                                                        call P3ab13(a, b, c0, c1, d0, d1,ACB(bi)%ASing,ACB(bi)%ATrip,ACB(bi)%ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        timer19 = timer19 +  clock_readwall(timer)

                                                        !---------------------P3a-4-P3b-2---------------------------------------------------------------
                                                        call clock_start(timer)
                                                        c0 = max(NI+1, a)
                                                        c1 = NIA
                                                        d0 = NI+1
                                                        d1 = NIA
                                                        !xxx A1 aaao                                                      
                                                        !call P3a4_P3b2(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .true.)
                                                        call P3ab42_1(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)

                                                        timer20 = timer20 +  clock_readwall(timer)
                                                        
                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_aoaa_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aoaa_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3a1_AI_aoaa_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a3_AI_aoaa_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3b1_AI_aoaa_ps_qsrp(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b3_AI_aoaa_qs_psrq(ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3c1_AI_aoaa_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c4_AI_aoaa_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)


                                                              !xxxn ic
                                                              c0 = NI+1
                                                              c1 = min(NIA, a)
                                                              d0 = NI+1
                                                              d1 = NIA
                                                              call P3ab24(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)

                                                        end if

                                                        !call P3a1_P3a2_P3b3_P3b_4(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2)

                                                        call clock_start(timer)
                                                        ! call  P45IA(Aux3A, b, a, map, rdm2_sum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        ! call  P45IA(Aux3B, b, a, map, rdm2_sum2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)

                                                        call  P45IA(Aux3X, b, a, map, rdm2_nsum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)

                                                        timer21 = timer21 +  clock_readwall(timer)

                                                        call clock_start(timer)
                                                        call p2_AAIA_qstt_ttqs(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAIA_qstt_ttsq(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAIA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)

                                                        timer22 = timer22 +  clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(3)) then ! (AA|AA) for ddot-1, for ddot-2
                                                  do b = NI+1, NIA
                                                        call clock_start(timer)
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, One)
                                                        timer23 = timer23 +  clock_readwall(timer)

                                                        

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A0aa
                                                        end if

                                                        call clock_start(timer)

                                                        call bare_int_prqs_aaaa_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'prqs', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        call bare_int_psqr_aaaa_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'psqr', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        timer24 = timer24 +  clock_readwall(timer)



                                                        !---------------------P3a-1,2-P3b-3,4---------------------------------------------------------------
                                                        ! do jakich blokow jesli a0
                                                        call clock_start(timer)
                                                        c0 = NI+1
                                                        c1 = min(NIA, a)
                                                        d0 = max(NI+1, b)
                                                        d1 = NIA
                                                        !xxx A0 aaaa
                                                        call P3ab13(a, b, c0, c1, d0, d1, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'p3ab1234_1', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        call P3ab24(a, b, c0, c1, d0, d1, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'p3ab1234_2', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        !call P3a1_P3a2_P3b3_P3b_4(a, b, c0, c1, d0, d1,ACB(bi)%ASing,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2)
                                                        timer25 = timer25 +  clock_readwall(timer)


                                                        !---------------------P3a-3-P3b-1---------------------------------------------------------------
                                                        call clock_start(timer)
                                                        c0 = NI+1
                                                        c1 = min(NIA, a)
                                                        d0 = NI+1
                                                        d1 = min(NIA, b)
                                                        !xxxx A0 aaaa
                                                        !call P3a3_P3b1(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .false.)
                                                        call P3ab31_1(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'p3ab31_1', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        timer26 = timer26 +  clock_readwall(timer)

                                                        !---------------------P3a-4-P3b-2---------------------------------------------------------------
                                                        call clock_start(timer)
                                                        c0 = max(NI+1, a)
                                                        c1 = NIA
                                                        d0 = max(NI+1, b)
                                                        d1 = NIA
                                                        !xxx A0 aaaa
                                                        !call P3a4_P3b2(a, b, c0, c1, d0, d1,ACB(bi)%ASing,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .false.)
                                                        call P3ab42_1(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        ! if (bi==A0aa)then
                                                        !       write(*, '(A15, F20.15)')'p3ab42_1', ACB(bi)%ASing(5,5)
                                                        ! end if

                                                        timer27 = timer27 +  clock_readwall(timer)

                                                        call clock_start(timer)

                                                        ! call P45(Aux3A, a, b, map, rdm2_sum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        ! call P45(Aux3B, a, b, map, rdm2_sum2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)

                                                        call P45(Aux3X, a, b, map, rdm2_nsum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        timer28 = timer28 +  clock_readwall(timer)
                                                        
                                                        call clock_start(timer)
                                                        call p2_AAAA_qstt_qstt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAAA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        timer29 = timer29 +  clock_readwall(timer)

                                                  end do
                                                  call clock_start(timer)
                                                  do r = 1, NA !(AAAA)
                                                        call clock_start(timer)
                                                        call real_aTv_x(V_axby, RAA(:,:,:), NChol, RAA(:,r, a-NI), NChol, NA**2, One, Zero)
                                                        timer30 = timer30 +  clock_readwall(timer)
                                                        
                                                        if (method == AC)then
                                                              bi = 1
                                                        else
                                                              bi = A0aa
                                                        end if
                                                        !xxx nic
                                                        ! A0
                                                        call clock_start(timer)
                                                        call P3c1_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                         ! if (bi==A0aa)then
                                                         !              write(*, '(A15, F20.15)')'P3c1', ACB(bi)%ASing(5,5)
                                                         !        end if
                                                        call P3c2_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        ! if (bi==A0aa)then
                                                        !        write(*, '(A15, F20.15)')'P3c2', ACB(bi)%ASing(5,5)
                                                        !  end if
                                                        call P3c3_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                         ! if (bi==A0aa)then
                                                         !       write(*, '(A15, F20.15)')'P3c3', ACB(bi)%ASing(5,5)
                                                         ! end if

                                                        call P3c4_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                         !if (bi==A0aa)then
                                                         !      write(*, '(A15, F20.15)')'P3c4', ACB(bi)%ASing(5,5)
                                                         !end if

                                                        timer31 = timer31 +  clock_readwall(timer)
                                                  end do

                                            end if
                                      end do loopAA

                                      ! loop for P3c2
                                      do a = 1, NI
                                            do r = 1, NI !(IIAA)
                                                  call clock_start(timer)
                                                  call real_aTv_x(V_axby, RAA(:,:,:), NChol, RII(:,r, a), NChol, NA**2, ACAlpha, Zero)
                                                  timer32 = timer32 +  clock_readwall(timer)

                                                  if (method == AC)then
                                                        bi = 1
                                                        !xxx nic
                                                        call P3c2_gamma(a, r, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                  end if
                                            end do

                                            do r = 1, NA !(IAAA)
                                                  call clock_start(timer)
                                                  call real_aTv_x(V_axby, RAA(:,:,:), NChol, RIA(:,r, a), NChol, NA**2, ACAlpha, Zero)
                                                  timer33 = timer33 +  clock_readwall(timer)
                                                  if (method == AC)then
                                                        bi = 1
                                                        !xxx nic                                                      
                                                        call P3c2_gamma(r+NI, a, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        !xxx nic                                                                                                                                                    
                                                        call P3c4_gamma(r+NI, a, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                  else
                                                        bi = A1aaao
                                                  end if

                                                  call clock_start(timer)
                                                  !xxx A1aaao
                                                  call P3c2_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                  timer34 = timer34 +  clock_readwall(timer)

                                                  !xxx A1 aaao
                                                  call P3c3_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                  timer35 = timer35 +  clock_readwall(timer)
                                            end do
                                      end do
                                case(loopIV)! loop 4
                                      loopIV: do a = 1, NI
                                            call clock_start(timer)
                                            if (a<b0IV.or.a>b1IV)then
                                                  call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NV), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, a), Zgk, NV, NChol, NTHC)
                                            end if
                                            timer36 = timer36 +  clock_readwall(timer)
                                            
                                            if (loops(ii)%run_inner(1)) then ! (IV|II)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        if (a<b0IV.or.a>b1IV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RII(:,:,b), NChol, NV, NI, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RII(:,:,b), NChol, NV, NI, NChol, AcAlpha)
                                                        end if
                                                        timer37 = timer37 +  clock_readwall(timer)	

                                                        call clock_start(timer)

                                                        call p2_IVII_qstt_qstt(BuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IVII_qstt_sqtt(BuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IVII_qtst_tsqt(AuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IVII_qtst_tqst(AuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        timer38 = timer38 +  clock_readwall(timer)	
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(2)) then !(IV|IA)
                                                  do b = 1, NI

                                                        call clock_start(timer)
                                                        if (a<b0IV.or.a>b1IV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                        end if
                                                        timer39 = timer39 +  clock_readwall(timer)	
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1vaoo
                                                        end if
                                                        call bare_int_prqs_vaoo_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call bare_int_psqr_vaoo_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer40 = timer40 +  clock_readwall(timer)
                                                        call clock_start(timer)
                                                        call p3a1_IA_vaoo_ps_qsrp(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3a3_IA_vaoo_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        call p3b1_IA_vaoo_ps_qsrp(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b3_IA_vaoo_qs_psrq(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        call p3c1_IA_vaoo_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3c4_IA_vaoo_ps_qrps(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer41 = timer41 +  clock_readwall(timer)	

                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_oova_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_oova_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3a2_AI_oova_qr_prsq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3a3_AI_oova_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3b2_AI_oova_qr_prsq(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b3_AI_oova_qs_psrq(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3c1_AI_oova_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3c5_AI_oova_qr_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        end if

                                                        call clock_start(timer)
                                                        call p2_IVIA_qtst_tstq(AuxII, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IVIA_qtst_tqts(AuxII, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)

                                                        timer42 = timer42 +  clock_readwall(timer)	
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(3))then ! (IV|AA) for ddot-2
                                                  do b = NI+1, NIA
                                                        call clock_start(timer)
                                                        if (a<b0IV.or.a>b1IV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                        end if
                                                        timer43 = timer43 +  clock_readwall(timer)

                                                        call clock_start(timer)


                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1vaao
                                                        end if
                                                        call bare_int_psqr_vaao_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer44 = timer44 +  clock_readwall(timer)
                                                        !if (bi==A1vaao)then
                                                              !print*, 'qa1', ACB(bi)%ASing(38,2)
                                                              !print*, 'qa1', ACB(bi)%ATripA(17,1)
                                                        !end if
                                                        
                                                        call clock_start(timer)
                                                        call p3a3_IA_vaao_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        call p3b3_IA_vaao_qs_psrq(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer45 = timer45 +  clock_readwall(timer)
                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa2', ACB(bi)%ASing(38,2)
                                                              !print*, 'qa2', ACB(bi)%ATripA(17,1)
                                                         !end if


                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_psqr_aova_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                              call p3a3_AI_aova_qs_psrq(ACB(bi)%ASing, ACB(bi)%ATrip, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call p3b3_AI_aova_qs_psrq(ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)

                                                        end if

                                                        call clock_start(timer)
                                                        call p2_IVAA_qstt_qstt(BuxAA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_IVAA_qstt_sqtt(BuxAA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)

                                                        timer46 = timer46 +  clock_readwall(timer)	
                                                  end do

                                                  do r = 1, NV !(IV|AA)
                                                        call clock_start(timer)
                                                        if (a<b0IV.or.a>b1IV)then
                                                              call  real_aTv_x(V_axby, RAA(:,:,:), NChol, Rkax(:, r),NChol, NA**2, ACAlpha, Zero)
                                                        else
                                                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RIV(:,r, a), NChol, NA**2, ACAlpha, Zero)
                                                        end if
                                                        timer47 = timer47 +  clock_readwall(timer)	
                                                        if (method == AC)then
                                                              bi = 1
                                                              !xxx nic                                                                                                                                         
                                                              call P3c4_gamma(r+NIA, a, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        else
                                                              bi = A1vaao
                                                        end if
                                                        !xxx A1vaao
                                                        call clock_start(timer)
                                                        call P3c3_gamma(a, r+NIA, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        timer48 = timer48 +  clock_readwall(timer)	
                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa3', ACB(bi)%ASing(38,2)
                                                         !     print*, 'qa3', ACB(bi)%ATripA(17,1)
                                                         !end if

                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(4))then  ! (IV|IV)

                                                  do b = b0IV, b1IV

                                                        call clock_start(timer)
                                                        if (a<b0IV.or.a>b1IV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                        end if
                                                        timer49 = timer49 +  clock_readwall(timer)	

                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                              ! else
                                                              !       bi = A1vvoo ! calculated OTF
                                                              
                                                              call bare_int_prqs_vvoo_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NIA+1, NBasis, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_vvoo_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NIA+1, NBasis, posS, posT, Occ, V_axby)

                                                              
                                                              ! if (method == AC) then
                                                              !       bi = 1
                                                              call bare_int_prqs_oovv_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NIA+1, NBasis, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_oovv_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NIA+1, NBasis, posS, posT, Occ, V_axby)

                                                        end if
                                                        call clock_start(timer)
                                                        call p2_IVIV_qtst_tqts(AuxII, a, b, NIA+1, NBasis, NIA+1, NBasis, NI, NIA, NBasis, Occ, V_axby)
                                                        timer50 = timer50 +  clock_readwall(timer)	                                                    
                                                  end do
                                            end if
                                      end do loopIV

                                case(loopVA) ! loop 5
                                      loopVA: do a = NIA+1, NBasis
                                            if (a<b0VA.or.a>b1VA)then
                                                  !call clock_start(timer)
                                                  call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NA), XXga(:, 1:NA), Xga(:, NI+1:NIA), Xga(:, a), Zgk, NA, NChol, NTHC)
                                                  !timer35 = timer35 + clock_readwall(timer)
                                            end if

                                            if (loops(ii)%run_inner(1)) then ! (VA|II)

                                                  do b = 1, NI

                                                        call clock_start(timer)
                                                        if (a<b0VA.or.a>b1VA)then
                                                              call real_aTb_x(V_axby, NA, Rkax, NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                        end if
                                                        timer51 = timer51 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        call p2_VAII_qstt_sqtt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_VAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        timer52 = timer52 + clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(2)) then !(VA|IA)

                                                  do b = 1, NI
                                                        
                                                        call clock_start(timer)

                                                        if (a<b0VA.or.a>b1VA)then
                                                              call real_aTb_x(V_axby, NA, Rkax, NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                        end if
                                                        timer53 = timer53 + clock_readwall(timer)

                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1vaao
                                                        end if

                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa40', ACB(bi)%ASing(38,2)
                                                         !     print*, 'qa40', ACB(bi)%ATripA(17,1)
                                                        !end if

                                                        call clock_start(timer)
                                                        call bare_int_prqs_vaao_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer54 = timer54 + clock_readwall(timer)
                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa4', ACB(bi)%ASing(38,2)
                                                         !     print*, 'qa4', ACB(bi)%ATripA(17,1)
                                                        !end if

                                                        call clock_start(timer)
                                                        call p3c1_IA_vaao_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer57 = timer57 + clock_readwall(timer)
                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa5', ACB(bi)%ASing(38,2)
                                                         !     print*, 'qa5', ACB(bi)%ATripA(17,1)
                                                         !end if

                                                        !---------------------P3a-1,2--P3b3,4---------------------------------------------------------------
                                                        c0 = NI+1
                                                        c1 = NIA
                                                        d0 = NI+1
                                                        d1 = NIA
                                                        !xxx A1vaao do pq1 rs1
                                                        call clock_start(timer)
                                                        call P3ab13(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  &
                                                              map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        timer55 = timer55 + clock_readwall(timer)
                                                        !if (bi==A1vaao)then
                                                              !       print*, 'qa6', ACB(bi)%ASing(38,2)
                                                        !      print*, 'qa6', ACB(bi)%ATripA(17,1)
                                                        ! end if

                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_aova_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              
                                                              call p3c1_AI_aova_qs_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              !xxx nic
                                                              c0 = NI+1
                                                              c1 = NIA
                                                              d0 = NI+1
                                                              d1 = NIA
                                                              call P3ab24(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                              call P3a1_P3a2_P3b3_P3b_4(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        end if

                                                        call clock_start(timer)
                                                        
                                                        call p2_VAIA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_VAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        timer56 = timer56 + clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(3)) then ! (VA|AA) for ddot-1, for ddot-2

                                                  do b = NI+1, NIA
                                                        call clock_start(timer)

                                                        if (a<b0VA.or.a>b1VA)then
                                                              call real_aTb_x(V_axby, NA, Rkax, NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                        end if
                                                        timer101 = timer101 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        else
                                                              bi = A1vaaa
                                                        end if
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa1', ACB(bi)%ATripA(3,1)
                                                        ! end if

                                                        call clock_start(timer)
                                                        call bare_int_prqs_vaaa_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer58 = timer58 + clock_readwall(timer)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa2X', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa2', ACB(bi)%ATripA(3,1)
                                                        ! end if
                                                        call clock_start(timer)
                                                        call bare_int_psqr_vaaa_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        timer59 = timer59 + clock_readwall(timer)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa3z', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa3', ACB(bi)%ATripA(3,1)
                                                        ! end if



                                                        !---------------------P3a-1,2-P3b-3,4---------------------------------------------------------------                                                         
                                                        c0 = NI+1
                                                        c1 = NIA
                                                        d0 = max(NI+1, b)
                                                        d1 = NIA
                                                        !xxx A1vaaa pq1 rs1
                                                        call clock_start(timer)
                                                        call P3ab13(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        timer60 = timer60 + clock_readwall(timer)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa4x', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa4', ACB(bi)%ATripA(3,1)
                                                        ! end if



                                                        !---------------------P3a-3-P3b-1---------------------------------------------------------------                                                      
                                                        c0 = NI+1
                                                        c1 = NIA
                                                        d0 = NI+1
                                                        d1 = min(NIA, b)                                                      
                                                        !xxx A1 vaaa
                                                        call clock_start(timer)
                                                        call P3ab31_2(a, b, c0, c1, d0, d1,ACB(bi)%ASing , ACB(bi)%ATrip, ACB(bi)%ATripA,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        timer61 = timer61 + clock_readwall(timer)
                                                        !call P3a3_P3b1(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .true.)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa5x', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa5', ACB(bi)%ATripA(3,1)
                                                        ! end if


                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_aava_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aava_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)


                                                              c0 = NI+1
                                                              c1 = NIA
                                                              d0 = max(NI+1, b)
                                                              d1 = NIA                                                            
                                                              call P3ab24(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                              !call P3a1_P3a2_P3b3_P3b_4(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2)
                                                              c0 = NI+1
                                                              c1 = NIA
                                                              d0 = NI+1
                                                              d1 = min(NIA, b)                                                            
                                                              call P3ab31_1(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        end if


                                                        call clock_start(timer)
                                                        ! call P45(Aux3A, a, b, map, rdm2_sum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        ! call P45(Aux3B, a, b, map, rdm2_sum2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)

                                                        call P45(Aux3X, a, b, map, rdm2_nsum, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        timer62 = timer62 + clock_readwall(timer)
                                                        call clock_start(timer)
                                                        call p2_VAAA_qstt_sqtt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_VAAA_qstt_qstt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_VAAA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_VAAA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        timer63 = timer63 + clock_readwall(timer)
                                                  end do
                                                  call clock_start(timer)
                                                  do r = 1, NA
                                                        call clock_start(timer)
                                                        if (a<b0VA.or.a>b1VA)then
                                                              call  real_aTv_x(V_axby, RAA(:,:,:), NChol, Rkax(:, r),NChol, NA**2, ACAlpha, Zero)
                                                        else
                                                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RVA(:,r, a-NIA), NChol, NA**2, ACAlpha, Zero)
                                                        end if
                                                        timer64 = timer64 + clock_readwall(timer)
                                                        if (method == AC)then
                                                              bi = 1
                                                              !xxx nic                                                                                                                                         
                                                              call P3c1_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                              ! xxx nic                                                                                                                                            
                                                              call P3c4_gamma(a, r+NI, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        else
                                                              bi = A1vaaa
                                                        end if

                                                        !xxx A1vaaa
                                                        call clock_start(timer)
                                                        call P3c1_gamma(r+NI, a, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        timer65 = timer65 + clock_readwall(timer)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa6x', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa6', ACB(bi)%ATripA(3,1)
                                                        ! end if

                                                        !xxx A1vaaa
                                                        call P3c3_gamma(r+NI, a, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        timer66 = timer66 + clock_readwall(timer)
                                                        ! if (bi==A1vaaa)then
                                                        !       print*, 'qa7x', ACB(bi)%ATrip(3,1)
                                                        !       print*, 'qa7', ACB(bi)%ATripA(3,1)
                                                              
                                                        ! end if

                                                  end do

                                                  timer103 = timer103 + clock_readwall(timer)

                                            end if
                                            if (loops(ii)%run_inner(4)) then ! (VA|IV)

                                                  do b = b0IV, b1IV
                                                        
                                                        call clock_start(timer)

                                                        if (a<b0VA.or.a>b1VA)then
                                                              call real_aTb_x(V_axby, NA, Rkax, NChol, RIV(:,:,b), NChol, NA, NV, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RIV(:,:,b), NChol, NA, NV, NChol, AcAlpha)
                                                        end if

                                                        timer67 = timer67 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                        ! else
                                                        !       bi = A1vvao calculated OTF
                                                        !end if

                                                              call bare_int_prqs_vvao_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NIA+1, NBasis, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_vvao_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NIA+1, NBasis, posS, posT, Occ, V_axby)


                                                        ! if (method == AC) then
                                                              !       bi = 1
                                                              call bare_int_prqs_aovv_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NIA+1, NBasis, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aovv_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NIA+1, NBasis, posS, posT, Occ, V_axby)

                                                        end if

                                                        timer68 = timer68 + clock_readwall(timer)

                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(5)) then ! (VA|VA) for ddot-1

                                                  do b = b0VA, b1VA
                                                        
                                                        call clock_start(timer)

                                                        if (a<b0VA.or.a>b1VA)then
                                                              call real_aTb_x(V_axby, NA, Rkax, NChol, RVA(:,:,b-NIA), NChol, NA, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RVA(:,:,b-NIA), NChol, NA, NA, NChol, AcAlpha)
                                                        end if

                                                        timer69 = timer69 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                              ! else
                                                              !       bi = A1vvaa Calculated OTF
                                                              !end if
                                                              call clock_start(timer)
                                                              call bare_int_prqs_vvaa_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_vvaa_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              
                                                              ! if (method == AC) then
                                                              !       bi = 1
                                                              call bare_int_psqr_vava_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_prqs_aavv_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_aavv_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NI+1, NIA, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              !---------------------P3a-3-P3b-1---------------------------------------------------------------
                                                              c0 = NI+1
                                                              c1 = NIA
                                                              d0 = NI+1
                                                              d1 = NIA
                                                              !xxx nic
                                                              !call P3a3_P3b1(a, b, c0, c1, d0, d1,ACB(bi)%ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, .false.)
                                          call P3ab31_1(a, b, c0, c1, d0, d1,ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                        end if
                                                        timer102 = timer102 + clock_readwall(timer)

                                                        call clock_start(timer)
                                                        call p2_VAVA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        timer70 = timer70 + clock_readwall(timer)
                                                  end do
                                            end if
                                      end do loopVA

                                case(loopVV) !loop 6
                                      loopVV: do a = NIA+1, NBasis
                                            call clock_start(timer)
                                            if (a<b0VV.or.a>b1VV)then
                                                  call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NV), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, a), Zgk, NV, NChol, NTHC)
                                                  !timer51 = timer51 + clock_readwall(timer)
                                            end if
                                            if (loops(ii)%run_inner(1)) then ! (VV|II)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RII(:,:,b), NChol, NV, NI, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RII(:,:,b), NChol, NV, NI, NChol, AcAlpha)
                                                        end if
                                                        timer71 = timer71 + clock_readwall(timer)

                                                        call clock_start(timer)
                                                        call p2_VVII_qstt_qstt(BuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                        timer72 = timer72 + clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(2)) then !(VV|IA)
                                                  do b = 1, NI
                                                        call clock_start(timer)
                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                        end if
                                                        timer100 = timer100 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        ! if (method == AC )then
                                                        !       bi = 1      
                                                        !       call p3c1_AI_VVIA_qs_prqs_18(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, posS, posT, Occ, V_axby)
                                                        !       call p3c1_IA_VVIA_qs_prqs_19(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, posS, posT, Occ, V_axby)
                                                        ! end if
                                                        !FALSE
                                                        timer57 = timer57 + clock_readwall(timer)

                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(3)) then !(VV|AA)   for ddot-2
                                                  do b = NI+1, NIA
                                                        call clock_start(timer)

                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                        end if
                                                        timer73 = timer73 + clock_readwall(timer)

                                                        call clock_start(timer)
                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_vava_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        end if
                                                        call clock_start(timer)
                                                        call p2_VVAA_qstt_qstt(BuxAA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)

                                                        timer74 = timer74 + clock_readwall(timer)
                                                        
                                                  end do


                                                  if (method == AC )then
                                                        bi = 1
                                                        do r = 1, NV
                                                              call clock_start(timer)
                                                              if (a<b0VV.or.a>b1VV)then
                                                                    ! Vw (1:NA, 1:NA) <= (pr|tu) |tu) part or p and r
                                                                    call  real_aTv_x(V_axby, RAA(:,:,:), NChol, Rkax(:, r),NChol, NA**2, ACAlpha, Zero)
                                                              else
                                                                    call real_aTv_x(V_axby, RAA(:,:,:), NChol, RVV(:,r, a-NIA), NChol, NA**2, ACAlpha, Zero)
                                                              end if
                                                              timer75 = timer75 + clock_readwall(timer)
                                                              call P3c1_gamma(a, r+NIA, ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA,map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act,  NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        end do
                                                  end if


                                            end if
                                            if (loops(ii)%run_inner(4)) then ! (VV|IV)

                                                  do b = b0IV, b1IV
                                                        call clock_start(timer)
                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                        end if
                                                        timer62 = timer62 + clock_readwall(timer)
                                                        call clock_start(timer)
                                                        
                                                        ! FALSEEEE
                                                        timer63 = timer63 + clock_readwall(timer)

                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(5)) then ! (VV|VA)

                                                  do b =b0VA, b1VA

                                                        call clock_start(timer)
                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RVA(:,:,b-NIA), NChol, NV, NA, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RVA(:,:,b-NIA), NChol, NV, NA, NChol, AcAlpha)
                                                        end if
                                                        timer68 = timer68 + clock_readwall(timer)
                                                        call clock_start(timer)
                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_vavv_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_prqs_vvva_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_vvva_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                              call bare_int_psqr_vavv_psqr(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NI+1, NIA, posS, posT, Occ, V_axby)
                                                        end if

                                                        timer69 = timer69 + clock_readwall(timer)
                                                  end do
                                            end if
                                            if (loops(ii)%run_inner(6)) then ! (VV|VV)

                                                  do b = b0VV, min(b1VV, a)
                                                        timer64 = timer64 + clock_readwall(timer)
                                                        call clock_start(timer)
                                                        if (a<b0VV.or.a>b1VV)then
                                                              call real_aTb_x(V_axby, NV, Rkax, NChol, RVV(:,:,b-NIA), NChol, NV, NV, NChol, AcAlpha)
                                                        else
                                                              call real_aTb_x(V_axby, NV, RVV(:,:,a-NIA), NChol, RVV(:,:,b-NIA), NChol, NV, NV, NChol, AcAlpha)
                                                        end if
                                                        timer65 = timer65 + clock_readwall(timer)
                                                        call clock_start(timer)

                                                        if (method == AC) then
                                                              bi = 1
                                                              call bare_int_prqs_vvvv_prqs(ACB(bi)%ASing, ACB(bi)%ATrip, ACB(bi)%ATripA, a, b, NIA+1, NBasis, NIA+1, NBasis, ind_virt, ndim_virt, PosS, PosT, Occ, V_axby)
                                                        end if
                                                        timer66 = timer66 + clock_readwall(timer)
                                                  end do
                                            end if
                                      end do loopVV
                                end select
                          end if
                    end do caseloop
              end do batchloop
            end  associate
#ifdef DEBUG
call xmsg('timer3', timer3, mdebugg)
call xmsg('timer3', timer3, mdebugg)
call xmsg('timer4', timer4, mdebugg)
call xmsg('timer7', timer7, mdebugg)
call xmsg('timer8', timer8, mdebugg)
call xmsg('timer10', timer10, mdebugg)
call xmsg('timer12', timer12, mdebugg)
call xmsg('timer11', timer11, mdebugg)
call xmsg('timer13', timer13, mdebugg)
call xmsg('timer104', timer104, mdebugg)
call xmsg('timer13', timer13, mdebugg)

call xmsg('timer16', timer16, mdebugg)
call xmsg('timer18', timer18, mdebugg)
call xmsg('timer17', timer17, mdebugg)
call xmsg('timer19', timer19, mdebugg)
call xmsg('timer20', timer20, mdebugg)
call xmsg('timer21', timer21, mdebugg)
call xmsg('timer22', timer22, mdebugg)

call xmsg('timer23', timer23, mdebugg)
call xmsg('timer24', timer24, mdebugg)
call xmsg('timer25', timer25, mdebugg)
call xmsg('timer26', timer26, mdebugg)
call xmsg('timer27', timer27, mdebugg)
call xmsg('timer28', timer28, mdebugg)
call xmsg('timer29', timer29, mdebugg)
call xmsg('timer30', timer30, mdebugg)
call xmsg('timer31', timer31, mdebugg)

call xmsg('timer32', timer32, mdebugg)
call xmsg('timer33', timer33, mdebugg)
call xmsg('timer34', timer34, mdebugg)
call xmsg('timer35', timer35, mdebugg)

call xmsg('timer37', timer37, mdebugg)
call xmsg('timer38', timer38, mdebugg)

call xmsg('timer39', timer39, mdebugg)
call xmsg('timer40', timer40, mdebugg)
call xmsg('timer41', timer41, mdebugg)
call xmsg('timer42', timer42, mdebugg)
call xmsg('timer43', timer43, mdebugg)

call xmsg('timer43', timer43, mdebugg)
call xmsg('timer44', timer44, mdebugg)
call xmsg('timer45', timer45, mdebugg)
call xmsg('timer46', timer46, mdebugg)
call xmsg('timer47', timer47, mdebugg)
call xmsg('timer48', timer48, mdebugg)

call xmsg('timer49', timer49, mdebugg)
call xmsg('timer50', timer50, mdebugg)

call xmsg('timer51', timer51, mdebugg)
call xmsg('timer52', timer52, mdebugg)

call xmsg('timer53', timer53, mdebugg)
call xmsg('timer54', timer54, mdebugg)
call xmsg('timer57', timer57, mdebugg)
call xmsg('timer55', timer55, mdebugg)
call xmsg('timer56', timer56, mdebugg)

call xmsg('timer101', timer101, mdebugg)
call xmsg('timer58', timer58, mdebugg)
call xmsg('timer59', timer59, mdebugg)
call xmsg('timer60', timer60, mdebugg)
call xmsg('timer61', timer61, mdebugg)
call xmsg('timer62', timer62, mdebugg)
call xmsg('timer63', timer63, mdebugg)
call xmsg('timer64', timer64, mdebugg)
call xmsg('timer65', timer65, mdebugg)
call xmsg('timer66', timer66, mdebugg)

call xmsg('timer67', timer67, mdebugg)
call xmsg('timer68', timer68, mdebugg)

call xmsg('timer69', timer69, mdebugg)
call xmsg('timer102', timer102, mdebugg)
call xmsg('timer70', timer70, mdebugg)

call xmsg('timer71', timer71, mdebugg)
call xmsg('timer72', timer72, mdebugg)
call xmsg('timer100', timer100, mdebugg)
call xmsg('timer73', timer73, mdebugg)
call xmsg('timer74', timer74, mdebugg)
call xmsg('timer75', timer75, mdebugg)
call xmsg('timer76', timer76, mdebugg)
#endif

              ! print*, 't56', timer56
              ! print*, 't57', timer57
              ! print*, 't58', timer58
              ! print*, 't59', timer59
              ! print*, 't60', timer60
              ! print*, 't61', timer61
              ! print*, 't62', timer62
              ! print*, 't63', timer63
              ! print*, 't64', timer64
              ! print*, 't65', timer65
              ! print*, 't66', timer66
              ! print*, 't67', timer67
              ! print*, 't68', timer68
              ! print*, 't69', timer69

              ! print*, ''

              ! print*, 't101', timer101
              ! print*, 't102', timer102
              ! print*, 't103', timer103
              ! print*, 't104', timer104
              ! print*, 'timeall', clock_readwall(timer0)



            end associate

      end subroutine THC_int_loop

      subroutine P45IA(Aux3A, s, u, map, rdm2_sum, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: Aux3A
            integer, intent(in) :: s, u
            integer, intent(in) ::NI, NA
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_sum
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision :: val, temp
            integer :: q, t, v

            do q = NI+1, NI+NA
                  call real_vw_x(val,  rdm2_sum(:,:,map(u),map(q)), V_axby, NA**2)
                  Aux3A(q, s)  = Aux3A(q, s) + val
            end do

      end subroutine P45IA


      subroutine P45(Aux3A, q, u, map, rdm2_sum, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: Aux3A
            integer, intent(in) :: q, u
            integer, intent(in) ::NI, NA
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_sum
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision :: val, val0
            integer :: s, t, v
            
            do s = NI+1, NI+NA
                  call real_vw_x(val,  rdm2_sum(:,:,map(s),map(u)), V_axby, NA**2)
                  Aux3A(s, q)  = Aux3A(s, q) + val
            end do

      end subroutine P45


      subroutine P3c1_gamma(p, r, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NI, NA)
            integer, intent(in) :: p, r, NI, NA
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision ::val1, val2

            integer :: s, q, t, u, pq, rs, pqT, rsT
            integer :: s1, q1

            s1 = min(r, NI+NA)
            q1 = min(p, NI+NA)

            !$omp parallel do default(shared) private(s, q, pq, rs, pqT, rsT, val1, val2)&
            !$omp collapse(2)
            do s = NI+1, s1
                  do q = NI+1, q1
                        pq = posS(p, q)
                        rs = posS(r, s)
                        if ((pq>0.and.rs>0))then
                              call real_vw_x(val1,  rdm2_pm_12_act(:,:,map(s),map(q)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pp_12_act(:,:,map(s),map(q)), V_axby, NA**2)
                              ASing(rs, pq) = ASing(rs, pq) - (val1+val2)
                              pqT = posT(p, q)
                              rsT = posT(r, s)
                              if ((pqT>0.and.rsT>0))then
                                    ATrip(rsT, pqT) = ATrip(rsT, pqT) - (val1+val2)
                                    ATripA(rsT, pqT) = ATripA(rsT, pqT) - (val1+val2)
                              end if
                        end if
                  end do
            end do
            !$omp end parallel do
      end subroutine P3c1_Gamma


      subroutine P3c2_gamma(p, r, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NI, NA)
            integer, intent(in) :: p, r, NI, NA
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision ::val1, val2

            integer :: s, q, t, u, pq, rs, pqT, rsT
            integer :: s0, q0
            
            s0 =  max(r, NI+1)
            q0 = max(p, NI+1)

            !$omp parallel do default(shared) private(s, q, pq, rs, pqT, rsT, val1, val2)&
            !$omp collapse(2)
            do s = s0, NI+NA
                  do q = q0, NI+NA
                        pq = posS(q, p)
                        rs = posS(s, r)
                        if ((pq>0.and.rs>0))then
                              call real_vw_x(val1,  rdm2_pm_12_act(:,:,map(s),map(q)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pp_12_act(:,:,map(s),map(q)), V_axby, NA**2)
                              ASing(rs, pq) = ASing(rs, pq) - (val1+val2)
                              pqT = posT(q, p)
                              rsT = posT(s, r)
                              if ((pqT>0.and.rsT>0))then
                                    ATrip(rsT, pqT) = ATrip(rsT, pqT) - (val1+val2)
                                    ATripA(rsT, pqT) = ATripA(rsT, pqT) - (val1+val2)
                              end if
                        end if
                  end do
            end do
            !$omp end parallel do
      end subroutine P3c2_Gamma

      subroutine P3c3_gamma(q, r, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NI, NA)
            integer, intent(in) :: q, r, NI, NA
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby

            double precision ::val1, val2

            integer :: s, p, t, u, pq, rs, pqT, rsT
            integer :: s1, p0


            s1 = min(NI+NA, r)
            p0 = max(q, NI+1)

            !$omp parallel do default(shared) private(s, p, pq, rs, pqT, rsT, val1, val2)&
            !$omp collapse(2)
            do s = NI+1, s1
                  do p = p0, NI+NA
                        pq = posS(p, q)
                        rs = posS(r, s)
                        if ((pq>0.and.rs>0))then
                              call real_vw_x(val1,  rdm2_pm_12_act(:,:,map(s),map(p)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pp_12_act(:,:,map(s),map(p)), V_axby, NA**2)
                              ASing(rs, pq) = ASing(rs, pq) - (val1+val2)
                              pqT = posT(p, q)
                              rsT = posT(r, s)
                              if ((pqT>0.and.rsT>0))then
                                    ATrip(rsT, pqT) = ATrip(rsT, pqT) + (val1+val2)
                                    ATripA(rsT, pqT) = ATripA(rsT, pqT) + (val1+val2)
                              end if
                        end if
                  end do
            end do
            !$omp end parallel do
      end subroutine P3c3_gamma


      subroutine P3c4_gamma(p, s, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NI, NA)
            integer, intent(in) :: p, s, NI, NA
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby

            double precision ::val1, val2

            integer :: q, r, t, u, pq, rs, pqT, rsT
            integer :: q1, r0

            q1 = min(NI+NA, p)
            r0 = max(s, NI+1)

            !$omp parallel do default(shared) private(q, r, pq, rs, pqT, rsT, val1, val2)&
            !$omp collapse(2)           
            do q = NI+1, q1
                  do r = r0, NI+NA
                        pq = posS(p, q)
                        rs = posS(r, s)
                        if ((pq>0.and.rs>0))then
                              call real_vw_x(val1,  rdm2_pm_12_act(:,:,map(r),map(q)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pp_12_act(:,:,map(r),map(q)), V_axby, NA**2)
                              ASing(rs, pq) = ASing(rs, pq) - (val1+val2)
                              pqT = posT(p, q)
                              rsT = posT(r, s)

                              if(p>q.and.r>s)then
                                    if ((pqT>0.and.rsT>0))then
                                          ATrip(rsT, pqT) = ATrip(rsT, pqT) + (val1+val2)
                                          ATripA(rsT, pqT) = ATripA(rsT, pqT) + (val1+val2)
                                    end if
                              end if
                        end if
                  end do
            end do
            !$omp end parallel do
      end subroutine P3c4_gamma



      subroutine P3a1_P3a2_P3b3_P3b_4(a, b, c0, c1, d0, d1, ASing,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act
             integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision ::val1, val2

            integer :: c, d, pq1, pq2, rs1, rs2

            do c = c0, c1
                  do d = d0, d1
                        pq1 = posS(d, b)
                        rs1 = posS(a, c)
                        pq2 = posS(a, c)
                        rs2 = posS(d, b)

                        if ((pq1>0.and.rs1>0).or.(pq2>0.and.rs2>0))then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2)
                        end if
                        if (pq1>0.and.rs1>0)then
                              ASing(rs1, pq1) = ASing(rs1, pq1) + (val1+val2)
                        end if
                        if (pq2>0.and.rs2>0)then
                              ASing(rs2, pq2) = ASing(rs2, pq2) + (val1+val2)
                        end if

                  end do
            end do
      end subroutine P3a1_P3a2_P3b3_P3b_4


      subroutine P3ab13(a, b, c0, c1, d0, d1, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision ::val1, val2, val2A

            integer :: c, d, pq, rs, pqT, rsT
            integer :: t, u

            !$omp parallel do default(shared) private(c, d, pq, rs, val1, val2,val2A, pqT, rsT ) &
            !$omp  collapse(2)         
            do c = c0, c1
                  do d = d0, d1
                        pq = posS(d, b)
                        rs = posS(a, c)

                        if (pq>0.and.rs>0)then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2) !czarne P3a1
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2) !zielone P3b3
                              call real_vw_x(val2A,  rdm2_pp_12_act(:,:,map(c),map(d)), V_axby, NA**2) !zielone P3b3 AAAA
                        end if
                        if (pq>0.and.rs>0)then
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (d>b.and.a>c)then
                                    pqT = posT(d, b)
                                    rsT = posT(a, c)
                                    ATrip(rsT, pqT) = ATrip(rsT, pqT) + (val1-val2)
                                    ATripA(rsT, pqT) = ATripA(rsT, pqT)  -val2A
                              end if
                        end if
                  end do
            end do
             !$omp end parallel do
      end subroutine P3ab13

      subroutine P3ab24(a, b, c0, c1, d0, d1, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision ::val1, val2, val2A

            integer :: c, d, pq, rs

            do c = c0, c1
                  do d = d0, d1

                        pq = posS(a, c)
                        rs = posS(d, b)

                        if (pq>0.and.rs>0)then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2) !czarne P3a2
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2) !zielone P3b4
                              call real_vw_x(val2A,  rdm2_pp_12_act(:,:,map(c),map(d)), V_axby, NA**2) !zielone P3b4AAAA
                        end if
                        if (pq>0.and.rs>0)then
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (a>c.and.d>b)then
                                    pq = posT(a, c)
                                    rs = posT(d, b)
                                    ATrip(rs, pq) = ATrip(rs, pq) + (val1-val2)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val2A
                              end if

                        end if
                  end do
            end do
      end subroutine P3ab24


      subroutine P3a4_P3b2(a, b, c0, c1, d0, d1, ASing,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, ns)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act
            double precision, dimension(:,:), intent(in) :: V_axby
            double precision, intent(out) ::val1, val2
            logical, intent(in) :: ns

            integer :: c, d, pq1, pq2, rs1, rs2

            do c = c0, c1
                  do d = d0, d1
                        pq1 = posS(d,b)
                        rs1 = posS(c,a)
                        if (pq1>0.and.rs1>0)then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2)
                        end if
                        if (pq1>0.and.rs1>0)then
                              ASing(rs1, pq1) = ASing(rs1, pq1) + (val1+val2)
                              if (ns == .true.) then
                                    ASing(pq1, rs1) = ASing(pq1, rs1) + (val1+val2)
                              end if
                        end if

                  end do
            end do

      end subroutine P3a4_P3b2

      subroutine P3ab42_1(a, b, c0, c1, d0, d1, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            
            double precision :: val1, val2, val2A
            integer :: c, d, pq, rs

            do c = c0, c1
                  do d = d0, d1
                        pq = posS(d,b)
                        rs = posS(c,a)
                        if (pq>0.and.rs>0)then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2) !zielone !P3a4
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2) !czarne  !P3b2
                              call real_vw_x(val2A,  rdm2_pp_12_act(:,:,map(c),map(d)), V_axby, NA**2) !P3b2AAAA
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (rs==5.and.pq==5)then
                                    print*, 'p3a4', val1
                                    print*, 'p3b2', val2
                              end if
                              
                              if (d>b.and.c>a)then
                                    pq = posT(d, b)
                                    rs = posT(c, a)
                                    ATrip(rs, pq) = ATrip(rs, pq) - (val1-val2)
                                    ATripA(rs, pq) = ATripA(rs, pq) + val2A
                              end if

                        end if
                  end do
            end do

      end subroutine P3ab42_1

      subroutine P3ab42_2(a, b, c0, c1, d0, d1, ASing, ATrip, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision :: val1, val2


            integer :: c, d, pq, rs

            do c = c0, c1
                  do d = d0, d1
                        pq = posS(c,a)
                        rs = posS(d,b)
                        if (pq>0.and.rs>0)then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2)!3a4
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2)!3b2
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (rs==5.and.pq==5)then
                                    print*, 'xp3a4', val1
                                    print*, 'xp3b2', val2
                              end if


                              if (d>b.and.c>a)then
                                    pq = posT(c, a)
                                    rs = posT(d, b)
                                    ATrip(rs, pq) = ATrip(rs, pq) - (val1-val2)
                              end if

                        end if
                  end do
            end do

      end subroutine P3ab42_2



      subroutine P3a3_P3b1(a, b, c0, c1, d0, d1, ASing,   map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, V_axby, NA, val1, val2, ns)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act
            double precision, dimension(:,:), intent(in) :: V_axby
            double precision, intent(out) ::val1, val2
            logical, intent(in) :: ns

            integer :: c, d, pq1, pq2, rs1, rs2


            do c = c0, c1
                  do d = d0, d1
                        pq1 = posS(a, c)
                        rs1 = posS(b, d)
                        if (pq1>0.and.rs1>0) then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2)
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2)
                        end if
                        if (pq1>0.and.rs1>0)then
                              ! P3a-3 P3b-1                                 
                              ASing(rs1, pq1) = ASing(rs1, pq1) + (val1+val2)
                              ! P3a-3 P3b-1                                                                                                                                                                      
                              if (ns ==.true.)then
                                    ASing(pq1, rs1) = ASing(pq1, rs1) + (val2-val1)
                              end if
                        end if

                  end do
            end do
      end subroutine P3a3_P3b1

      subroutine P3ab31_1(a, b, c0, c1, d0, d1, ASing,ATrip,ATripA,  map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision :: val1, val2, val2A

            integer :: c, d, pq, rs


            do c = c0, c1
                  do d = d0, d1
                        pq = posS(a, c)
                        rs = posS(b, d)
                        if (pq>0.and.rs>0) then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2) !zielone P3a3 
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2) !czarne P3b1
                              call real_vw_x(val2A,  rdm2_pp_12_act(:,:,map(c),map(d)), V_axby, NA**2) !czarne P3b1
                              ! P3a-3 P3b-1
                              if (rs==5.and.pq==5)then
                                    print*, 'p3a3', val1
                                    print*, 'p3b1', val2
                              end if
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (b>d.and.a>c)then
                                    pq = posT(a, c)
                                    rs = posT(b, d)
                                    ATrip(rs, pq) = ATrip(rs, pq) - (val1-val2)
                                    ATripA(rs, pq) = ATripA(rs, pq) + val2A
                              end if

                        end if
                  end do
            end do
      end subroutine P3ab31_1


      subroutine P3ab31_2(a, b, c0, c1, d0, d1, ASing, ATrip, ATripA, map, posS, posT, rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act, x0, x1, y0, y1, V_axby, NA)
            integer, intent(in) :: a, b, NA
            integer, intent(in) :: c0, c1, d0, d1
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, dimension(:), intent(in) ::  map
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_pm_13_act, rdm2_pm_12_act, rdm2_pp_12_act
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1,y0:y1), intent(in) :: V_axby
            double precision :: val1, val2, val2A

            integer :: c, d, pq, rs


            do c = c0, c1
                  do d = d0, d1
                        pq = posS(b, d)
                        rs = posS(a, c)
                        if (pq>0.and.rs>0) then
                              call real_vw_x(val1,  rdm2_pm_13_act(:,:,map(d),map(c)), V_axby, NA**2) !zielone P3a3
                              call real_vw_x(val2,  rdm2_pm_12_act(:,:,map(c),map(d)), V_axby, NA**2) !czarne P3b1_2
                              call real_vw_x(val2A,  rdm2_pp_12_act(:,:,map(c),map(d)), V_axby, NA**2) !czarne P3b1_2AAAA
                              ASing(rs, pq) = ASing(rs, pq) + (val1+val2)

                              if (b>d.and.a>c)then
                                    pq = posT(b, d)
                                    rs = posT(a, c)
                                    ATrip(rs, pq) = ATrip(rs, pq) - (val1-val2)
                                    ATripA(rs, pq) = ATripA(rs, pq) + val2A
                              end if

                        end if
                  end do
            end do
      end subroutine P3ab31_2


      
      subroutine parse_config(loops, integral_list, TabChol, ACBlockList, method)
            type(TLoop), dimension(6), intent(out) :: loops              
            integer, dimension(2, 19), intent(in) :: integral_list
            integer, dimension(6), intent(out) :: TabChol
            integer, intent(in) :: method
            integer :: i, j, main_loop, inner_loop, k
            integer, dimension(:,:,:), intent(out) :: ACBlockList




            ! method = -1 for future
            ! method = 0 AC0 (A0 + A1)
            ! method = 1 All

            TabChol = 0


            do i = 1, 6
                  do j = 1, i
                        loops(i)%run_main = .false.
                        loops(i)%run_inner(j) = .false.
                  end do
            end do

            if (method == 1) then


                  do i = 1, 6
                        do j = 1, i
                              loops(i)%run_main = .true.
                              loops(i)%run_inner(j) = .true.
                              ACBlockList(i, j, :) = [1]
                        end do
                  end do
                  ! loops(6)%run_inner(4) = .false. !only integrals vviv not possible in pp (pq cannot be vi)
                  ! loops(6)%run_inner(2) = .false.

                  TabChol = 1 ! Calulate all Cholesky vectors

                  

            else  if (method == 0)then

                  do i = 1, 5
                        do j = 1, i
                              loops(i)%run_main = .true.
                              loops(i)%run_inner(j) = .true.
                        end do
                  end do

                  loops(6)%run_inner(2) = .false.
                  loops(6)%run_inner(4) = .false.
                  loops(6)%run_inner(5) = .false.
                  loops(6)%run_inner(6) = .false.

                  loops(5)%run_inner(4) = .false.
                  loops(5)%run_inner(5) = .false.

                  TabChol = 1 ! Calculate all Cholesky vectors
                  TabChol(6) = 0
                  
                  ACBlockList = 0

                  ACBlockList(1, 1, :) = [1, 3, 13, 0, 0, 0]
                  ACBlockList(1, 1, :) = [1, 3, 13, 0, 0, 0]

                  ACBlockList(2, 1, :) = [12, 13, 0, 0, 0, 0]
                  ACBlockList(2, 2, :) = [1, 3, 4, 5, 12, 13]

                  ACBlockList(3, 1, :) = [1, 3, 4, 5, 12, 13]
                  ACBlockList(3, 2, :) = [12, 13, 0, 0, 0, 0]
                  ACBlockList(3, 3, :) = [13, 0, 0, 0, 0, 0]

                  ACBlockList(4, 1, :) = [9, 10, 11, 0, 0, 0]
                  ACBlockList(4, 2, :) = [9, 10, 11, 0, 0, 0]
                  ACBlockList(4, 3, :) = [9, 10, 11, 0, 0, 0]
                  ACBlockList(4, 4, :) = [2, 4, 6, 0, 0, 0]

                  ACBlockList(5, 1, :) = [9, 10, 11, 0, 0, 0]
                  ACBlockList(5, 2, :) = [9, 10, 11, 0, 0, 0]
                  ACBlockList(5, 3, :) = [7, 0, 0, 0, 0, 0]
                  ACBlockList(5, 4, :) = [10, 11, 0, 0, 0, 0]
                  ACBlockList(5, 5, :) = [2, 4, 8, 0, 0, 0]

                  ACBlockList(6, 1, :) = [2, 4, 0, 0, 0, 0]
                  ACBlockList(6, 4, :) = [2, 4, 0, 0, 0, 0]



            else if (method <0) then
                  do i = 1, 19
                        main_loop = integral_list(1, i)
                        inner_loop = integral_list(2, i)

                        if (main_loop == -1) exit

                        loops(main_loop)%run_main = .true.
                        loops(main_loop)%run_inner(inner_loop) = .true.
                  end do

                  do i = 1, 6
                        if (loops(i)%run_main)then
                              select case(i)
                              case(1)
                                    TabChol(1) = 1
                              case(2)
                                    TabChol(2) = 1
                                    if (loops(i)%run_inner(1)) TabChol(1) = 1
                              case(3)
                                    TabChol(3) = 1
                                    if (loops(i)%run_inner(1)) TabChol(1) = 1
                                    if (loops(i)%run_inner(2)) TabChol(2) = 1
                              case(4)
                                    TabChol(4) = 1
                                    if (loops(i)%run_inner(1)) TabChol(1) = 1
                                    if (loops(i)%run_inner(2)) TabChol(2) = 1
                                    if (loops(i)%run_inner(3)) TabChol(3) = 1
                              case(5)
                                    TabChol(5) = 1
                                    if (loops(i)%run_inner(1)) TabChol(1) = 1
                                    if (loops(i)%run_inner(2)) TabChol(2) = 1
                                    if (loops(i)%run_inner(3)) TabChol(3) = 1
                                    if (loops(i)%run_inner(4)) TabChol(4) = 1
                              case(6)
                                    TabChol(6) = 1
                                    if (loops(i)%run_inner(1)) TabChol(1) = 1
                                    if (loops(i)%run_inner(2)) TabChol(2) = 1
                                    if (loops(i)%run_inner(3)) TabChol(3) = 1
                                    if (loops(i)%run_inner(4)) TabChol(4) = 1
                                    if (loops(i)%run_inner(5)) TabChol(5) = 1
                              end select
                        end if
                  end do
            end if

      end subroutine parse_config


      subroutine update_HNO_AC0(AuxData, HNO, wij, wvw, AuxII, BuxII)
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:, :), intent(out) :: HNO
            double precision, dimension(:), intent(in) :: wij, wvw
            double precision, dimension(:,:), intent(in) :: AuxII, BuxII
            integer :: q, s, i
            double precision :: etot

            associate(NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux, NI=>AuxData%NI, NIA=>AuxData%NIA, NV=>AuxData%NV, HNO0=>AuxData%HNO0)

              
              HNO = zero

              do q = 1, NI
                    HNO(q, q) = wij(q)
              end do

              do q = 1, NV
                    HNO(q+NIA, q+NIA) = wvw(q)
              end do
                    
              do q = NI+1, NIA
                    do s = NI+1, NIA

                          HNO(q, s) = HNO0(q, s)

                          HNO(q,s) = HNO(q, s) + (Two * BuxII(q, s)  - AuxII(q,s))
                          
                    end do
              end do




            etot = 0
            do i = 1, NBasis
                  etot=etot + two* AuxData%Occ(i) * HNO0(i, i)
            !      write(*, '(I5, 3F40.20)')i, hno0(i,i), AuxData%Occ(i), AuxData%Occ(i)*two*hno0(i,i)
            end do

            AuxData%ECAS_oneelectr = etot
          end associate
            
        end subroutine update_HNO_AC0

        
      subroutine update_HNO(AuxData, HNO0, HNO, AuxII, AuxAA, BuxII, BuxAA, ACAlpha, method)
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:, :), intent(inout) ::HNO0, HNO
            double precision, dimension(:,:), intent(in) :: AuxII, AuxAA, BuxII, BuxAA
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: method
            double precision :: auxval
            integer :: q, s, i, j


            call xmsg('ACAlpha',   ACAlpha, mverbose)
            call ximsg('method', method, mverbose)


            associate(NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux)


              if (method == 0)then

                    HNO = zero

                    do q = 1, NBasis
                          do s = 1, Nbasis

                                if (IndAux(q)==IndAux(s))then
                                      HNO(q, s) = HNO0(q, s)

                                      ! if (abs(HNO0(q,s)).gt.1.d-5)then
                                      !       print*, 'pluszek', q, s, hno0(q,s)
                                      ! end if

                                      if(IndAux(q)==0.or.IndAux(q)==2)then
                                            HNO(q,s) = HNO(q, s) + ( Two * BuxAA(q,s)- AuxAA(q,s))
                                      end if
                                      HNO(q,s) = HNO(q, s) + (Two * BuxII(q, s)  - AuxII(q,s))
                                      ! if (abs(HNO0(q,s)).gt.1.d-5)then
                                      !       print*, 'pluszek2', q, s, hno(q,s), (Two * BuxII(q, s)),   - AuxII(q,s)
                                      ! end if

                                end if
                          end do
                    end do

                    ! do q = 1, NBasis
                    !       do s = 1, Nbasis
                    !             if (abs(hno(q,s)).gt.1.d-5)then
                    !                   print*, 'haa', q, s, hno(q,s)
                    !             end if
                    !       end do
                    ! end do


              else
                    do q = 1, NBasis
                          do s = 1, Nbasis
                                if (IndAux(q).ne.IndAux(s))then
                                      HNO(q, s) = ACAlpha * HNO0(q, s)
                                else
                                      HNO(q,s) = HNO0(q,s)
                                end if
                          end do
                    end do

                    AuxVal = (One/ACAlpha -One)

                    do q = 1, NBasis
                          do s = 1, Nbasis
                                if (IndAux(q)==IndAux(s))then                                      
                                      if(IndAux(q)==0.or.IndAux(q)==2)then
                                            HNO(q,s) = HNO(q,s) + AuxVal * ( Two * BuxAA(q,s)- AuxAA(q,s))
                                      end if
                                      HNO(q,s) = HNO(q,s) + AuxVal * (Two * BuxII(q, s)  - AuxII(q,s))
                                end if
                          end do
                    end do
              end if



              !               HNO = HNO0
              !               print*, 'hnohno', hno(2,2)
              !                 AuxVal = (One/ACAlpha -One)
              !                 do j=1,AuxData%NBasis
              !                       do i=1,AuxData%NBasis
              !                             if(AuxData%IndAux(i)/=AuxData%IndAux(j)) HNO(i,j) = ACAlpha*HNO(i,j)
              !                       enddo
              !                 enddo
              !                  print*, 'hnohno2',	hno(2,2)
              ! !          end if
              !                     !AuxVal = (One/ACAlpha -One)
              !               associate(NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux)

              !                 do q = 1, NBasis
              !                       do s = 1, Nbasis
              !                             if (IndAux(q)==IndAux(s))then
              !                                   if(IndAux(q)==0.or.IndAux(q)==2)then
              !                                         HNO(q,s) = HNO(q,s) + AuxVal * ( Two * BuxAA(q,s)- AuxAA(q,s))
              !                                   end if
              !                                   HNO(q,s) = HNO(q,s) + AuxVal * (Two * BuxII(q, s)  - AuxII(q,s))
              !                             end if
              !                       end do
              !                 end do
              !                  print*, 'hnohno3',	hno(2,2)

              ! do q = 1, NBasis
              !       do s = 1, Nbasis
              !             if (abs(hno(q,s)).gt.1.d-5)then
              !                   print*, 'haa-zero', q, s, hno(q,s)
              !             end if
              !       end do
              ! end do

              ! do q = 1, NBasis
              !       do s = 1, Nbasis
              !             if (abs(hno0(q,s)).gt.1.d-5)then
              !                   print*, 'haa-one', q, s, hno0(q,s)
              !             end if
              !       end do
              ! end do



            end associate

      end subroutine update_HNO


      subroutine PPERPA0_THC_OOVV(AuxData, ACBlock, HNO, posS, posT, AuxAA, BuxAA, Aux3X, Aux3B, ACAlpha, NDim1, IndN1)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            double precision, dimension(:, :), intent(inout) ::HNO
            double precision, dimension(:,:), intent(inout) :: AuxAA,  BuxAA,Aux3X, Aux3B
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: NDim1
            integer, dimension(:,:), intent(in) ::IndN1

            integer :: p, q
            integer :: pq, pqT
            double precision :: val, vals, valt, val1, val0
            double precision :: NumF, NumH
            integer :: lll


            lll = 1

            j_colloops: do pq = 1, NDim1

                  p = IndN1(1, pq)
                  q = IndN1(2, pq)


                  val = zero
                  vals = zero

                  val = val + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(q)) * HNO(q,q)
                  val = val + (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(p)) * HNO(p,p)

                  vals = val
                  valt = val
                  val = zero                  

                  if (p==q)then

                        val = val - (one- AuxData%Occ(p)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(p)) * HNO(p,p)
                  end if
                  val = two * val

                  vals = vals -val
                  valt = valt + val

                  NumF = one
                  NumF = one/(1-AuxData%Occ(p)-AuxData%Occ(q))

                  NumH = one
                  if (p==q)then
                        NumH = frac12
                  end if

                  ACBlock%ASing(pq, 1) = NumF * NumH * ACBlock%ASing(pq, 1) + NumF * NumH*vals

                  if (p>q)then
                        pqT = posT(p, q)
                        ACBlock%ATrip(pqT, 1) = NumF * ACBlock%ATrip(pqT, 1) + NumF  * valt
                        ACBlock%ATripA(pqT, 1) = NumF * ACBlock%ATripA(pqT, 1) + NumF  * valt
                  end if

                  ! if (abs(ACBlock%ASing(pq, 1)).gt.1.d-5)then
                  !       write(*,'(2I5, A5, 4I5,2F30.16)') pq, q, '  |  ', p, q, p, q, ACBlock%ASing(pq, 1)
                  ! end if

            end do j_colloops


      end subroutine PPERPA0_THC_OOVV



      subroutine PPERPA0_THC_OA_Eigs(AuxData, ACBlock, HNO, AuxAA, BuxAA, Aux3X, Aux3B, ACAlpha, NDim1,  IndN1)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            double precision, dimension(:, :), intent(inout) ::HNO
            double precision, dimension(:,:), intent(inout) :: AuxAA,  BuxAA,Aux3X, Aux3B
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: NDim1
            integer, dimension(:,:), intent(in) ::IndN1

            integer :: p, q, r, s
            integer :: pq, rs
            double precision :: val, vals, valt
            double precision :: NumF, NumH
            integer :: lll, k, i, j, irel, jrel
            integer :: istart, idim, offset
            double precision, dimension(:,:), allocatable :: work, swork
            double precision, dimension(:,:), allocatable :: temp_work, temp_swork
            double precision, dimension(:), allocatable :: eigswork
            integer, dimension(:), allocatable :: vpluswork
            integer :: nd


            associate(NI=>AuxData%NI, NIA=>AuxData%NIA, NA=>AuxData%NA, NV=>AuxData%NV)


              allocate(work(NA, NA))
              allocate(swork(NA, NA))
              allocate(vpluswork(NA))
              allocate(eigswork(NA))
              swork = zero


              do q = 1, NI
                    offset = NA * (q-1)
                    s = q
                    swork = zero
                    do j = 1, ACBlock%MiniBlocks(q)%NdimS
                          p = ACBlock%MiniBlocks(q)%ActListS(j)
                          do i = 1, ACBlock%MiniBlocks(q)%NdimS
                                r = ACBlock%MiniBlocks(q)%ActListS(i)

                                
                                val = zero
                                vals = zero

                                if (p==r) then
                                      val = val + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s)
                                end if

                                val = val+ (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(r)) * HNO(p,r)

                                if (AuxData%IndAux(p)==1.and. AuxData%IndAux(r)==1)then

                                      val = val + two *  BuxAA(p,r)
                                      val = val  - AuxAA(p,r)

                                      ! val = val - frac12 * Aux3A(p, r)
                                      ! val = val - frac12 * Aux3B(r, p)

                                      val = val - frac12 * Aux3X(p, r)
                                      val = val - frac12 * Aux3X(r, p)

                                      if (AuxData%IndAux(q)==0)then
                                            val = val + AuxData%Occ(q)*( AuxAA(p,r)  - two * BuxAA(p, r))
                                      end if
                                end if

                                NumF = one
                                NumF = one/(1-AuxData%Occ(r)-AuxData%Occ(s))
                                if (i==j)then
                                      swork(i, j) = NumF
                                else
                                      swork(i,j)=zero
                                end if

                                work(i, j) = numF * val
                                
                                ! if (abs(work(i,j)).gt.1.d-5)then
                                !       write(*,'(2I5, A5, 4I5,F30.16)') i, j, '  |  ', r, s, p, q, work(i,j)
                                ! end if
                          end do
                    end do

                    nd = ACBlock%MiniBlocks(q)%NdimS

                    call ximsg('wymiar tego bloczku to', nd, mdebug)

                    allocate(temp_work(nd, nd))
                    allocate(temp_swork(nd, nd))

                    temp_work = work(1:nd, 1:nd)
                    temp_swork = Swork(1:nd, 1:nd)
                    call NonSymEigBlockTHC(eigswork(1:nd), temp_work, temp_swork, vpluswork(1:nd), AuxData%IndAux, nd)
                    work(1:nd, 1:nd) = temp_work
                    Swork(1:nd, 1:nd) = temp_swork
                    deallocate(temp_work)
                    deallocate(temp_swork)
                    call ximsg('wartosci wlasne tego bloczku', nd, mdebug)



                    do j = 1, ACBlock%MiniBlocks(q)%NdimS
                          p = ACBlock%MiniBlocks(q)%ActListS(j)-NI
                          ACBlock%MiniBlocks(q)%EigS(p) = eigswork(j)
                          call xmsg('', eigswork(j), mdebug)

                         do i = 1, ACBlock%MiniBlocks(q)%NdimS
                               r = ACBlock%MiniBlocks(q)%ActListS(i)-NI
                               ACBlock%MiniBlocks(q)%MiniAVs(r,p) = work(i,j)
                               
                         end do
                   end do

             end do
             
             ! print*, 'wartosci wlasne bloku A0oaoa'
             ! do i = 1, NI
             !       print*, 'dla bloku', i
             !       do j = 1, NA
             !             print*, ACBlock%MiniBlocks(i)%EigS(j)
             !       end do
                   
             ! end do

             deallocate(work)
             deallocate(swork)
             deallocate(vpluswork)
             deallocate(eigswork)

          end associate

    end subroutine PPERPA0_THC_OA_EIGS


    subroutine PPERPA0_THC_VA_Eigs(AuxData, ACBlock, HNO, AuxAA, BuxAA, Aux3X, Aux3B, ACAlpha, NDim1,  IndN1)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            double precision, dimension(:, :), intent(inout) ::HNO
            double precision, dimension(:,:), intent(inout) :: AuxAA,  BuxAA,Aux3X, Aux3B
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: NDim1
            integer, dimension(:,:), intent(in) ::IndN1

            integer :: p, q, r, s
            integer :: pq, rs
            double precision :: val, vals, valt
            double precision :: NumF, NumH
            integer :: lll, k, i, j, irel, jrel
            integer :: istart, idim, offset
            double precision, dimension(:,:), allocatable :: work, swork
            double precision, dimension(:,:), allocatable :: temp_work, temp_swork
            double precision, dimension(:), allocatable :: eigswork
            integer, dimension(:), allocatable :: vpluswork
            integer :: nd


            associate(NI=>AuxData%NI, NIA=>AuxData%NIA, NA=>AuxData%NA, NV=>AuxData%NV)


              allocate(work(NA, NA))
              allocate(swork(NA, NA))
              allocate(vpluswork(NA))
              allocate(eigswork(NA))
              swork = zero


              do i = 1, NV
                    offset = NA * (i-1)
                    p = NIA + i
                    r = p
                    swork = zero
                    do j = 1, ACBlock%MiniBlocks(i)%NdimS
                          q = ACBlock%MiniBlocks(i)%ActListS(j)
                          do k = 1, ACBlock%MiniBlocks(i)%NdimS
                                s = ACBlock%MiniBlocks(i)%ActListS(k)
                                
                                val = zero
                                vals = zero

                                val = val + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s)
                                val = val + two * (BuxAA(q, s))
                                val = val  - AuxAA(q, s)

                                ! val = val - frac12 * Aux3A(q, s)
                                ! val = val - frac12 * Aux3B(s, q)

                                val = val - frac12 * Aux3X(q, s)
                                val = val - frac12 * Aux3X(s, q)

                                if (q==s)then
                                      val = val+ (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(r)) * HNO(p,r)
                                end if

                                vals = val
                                valt = val

                                NumF = one
                                NumF = one/(1-AuxData%Occ(r)-AuxData%Occ(s))
                                if (k==j)then
                                      swork(k, k) = NumF
                                else
                                      swork(k,j) = zero
                                end if

                                work(k, j) = numF * vals
    
                               !  if (abs(work(k,j)).gt.1.d-5)then
                               !        write(*,'(2I5, A5, 4I5,3F30.16)') i, j, '  |  ', r, s, p, q, numF, vals, work(k,j)
                               ! end if
                                

                          end do
                    end do
                    
                    nd = ACBlock%MiniBlocks(i)%NdimS
                    call ximsg('wymiar tego bloczku to', nd, mdebug)

                    allocate(temp_work(nd, nd))
                    allocate(temp_swork(nd, nd))

                    temp_work = work(1:nd, 1:nd)
                    temp_swork = Swork(1:nd, 1:nd)
                    
                    call NonSymEigBlockTHC(eigswork(1:nd), temp_work(1:nd, 1:nd), temp_swork(1:nd, 1:nd), vpluswork(1:nd), AuxData%IndAux, nd)

                    work(1:nd, 1:nd) = temp_work
                    Swork(1:nd, 1:nd) = temp_swork
                    deallocate(temp_work)
                    deallocate(temp_swork)
                    call ximsg('wartosci wlasne tego bloczku', nd, mdebug)
                    do j = 1, ACBlock%MiniBlocks(i)%NdimS
                         q = ACBlock%MiniBlocks(i)%ActListS(j)-NI
                         ACBlock%MiniBlocks(i)%EigS(q) = eigswork(j)
                         call xmsg('', eigswork(j), mdebug)
                         do k = 1, ACBlock%MiniBlocks(i)%NdimS
                               s = ACBlock%MiniBlocks(i)%ActListS(k)-NI
                               ACBlock%MiniBlocks(i)%MiniAVs(s,q) = work(k,j)
                         end do
                   end do

             end do
             
             ! call ximsg('wartosci wlasne blocku A0vava', nd, mdebug)
!              print*, 'wart wlasne ac0vava'
!              do i = 1, NV

!                    do j = 1, NA
!                          print*, ACBlock%MiniBlocks(i)%EigS(j)
! !                         call xmsg('', ACBlock%MiniBlocks(i)%EigS(j), mdebug)
!                    end do

!              end do

             deallocate(work)
             deallocate(swork)
             deallocate(vpluswork)
             deallocate(eigswork)

          end associate

    end subroutine PPERPA0_THC_VA_EIGS




    subroutine PPERPA0_THC_AA(AuxData, ACBlock, HNO, posS, posT, AuxAA, BuxAA, Aux3X, Aux3B, ACAlpha, NDim1,  IndN1)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            double precision, dimension(:, :), intent(inout) ::HNO
            double precision, dimension(:,:), intent(in) :: AuxAA,  BuxAA,Aux3X, Aux3B
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: NDim1
            integer, dimension(:,:), intent(in) ::IndN1

            integer :: p, q, r, s
            integer :: pq, rs, pqT, rsT
            double precision :: val, vals, valt, val1, val0
            double precision :: NumF, NumH
            double precision :: u1, u2, u3, u4, u5

            integer :: lll


            lll = 1

            j_colloops: do pq = 1, NDim1
                  i_rowloops: do rs = 1, NDim1

                        r = IndN1(1, rs)
                        s = IndN1(2, rs)

                        p = IndN1(1, pq)
                        q = IndN1(2, pq)


                        val = zero
                        vals = zero

                        if (p==r) then

                              val = val + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s)

                              
                              val = val + two * ( BuxAA(q, s)) !P2a-1                             
                              val = val  - AuxAA(q, s)   !P2a-1

                              u1 =  two * ( BuxAA(q, s)) !P2a-1                             
                              u2 =   - AuxAA(q, s)   !P2a-1

                              
                              ! val = val - frac12 * Aux3A(q, s)
                              ! val = val - frac12 * Aux3B(s, q)

                              val = val - frac12 * Aux3X(q, s)
                              val = val - frac12 * Aux3X(s, q)
                              
                              u3 =  - frac12 * Aux3X(q, s)
                              u4 =  - frac12 * Aux3X(s, q)

                              if (rs==5.and.pq==5)then
                                    print*, 'kupsko', (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s)
                                    write(*, '( 4F20.10)') u1, u2, u3, u4
                                    
                              end if

                        end if
                        
                        
                        if (q==s)then
                              val = val+ (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(r)) * HNO(p,r)


                              val = val + two * (BuxAA(p,r))
                              val = val  - AuxAA(p,r)
                              
                              u1 =  two * (BuxAA(p,r))
                              u2 =   - AuxAA(p,r)

                              !val = val - frac12 * Aux3A(p, r)
                              !val = val - frac12 * Aux3B(r, p)

                              val = val - frac12 * Aux3X(p, r)
                              val = val - frac12 * Aux3X(r, p)

                              
                              u3 =  - frac12 * Aux3X(p, r)
                              u4 =  - frac12 * Aux3X(r, p)
                              if (rs==5.and.pq==5)then
                                    print*, 'kupsko2', (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(r)) * HNO(p,r)
                                    write(*, '( 4F20.10)') u1, u2, u3, u4
                                    print*, HNO(p,r)
                              end if

                        end if

                        vals = val
                        valt = val
                        val = zero

                        if (p==s)then

                              val = val - (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(r)) * HNO(q,r)


                              val = val - two     * (BuxAA(q,r ))
                              val = val + AuxAA(q,r)

                              u1 =  - two     * (BuxAA(q,r ))
                              u2 =  + AuxAA(q,r)

                              !val = val + frac12 * Aux3B(r, q)
                              !val = val + frac12 * Aux3A(q, r)

                              val = val + frac12 * Aux3X(r, q)
                              val = val + frac12 * Aux3X(q, r)

                              u3 =  frac12 * Aux3X(r, q)
                              u4 =  frac12 * Aux3X(q, r)

                              if (rs==5.and.pq==5)then
                                    print*, 'kupsko3', (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(r)) * HNO(q,r)
                                    write(*, '( 4F20.10)') u1, u2, u3, u4
                              end if

                        end if



                        if (q==r)then
                              val = val - (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(s)) * HNO(p,s) 

                              val = val - two     * BuxAA(p,s) 
                              val = val  + AuxAA(p,s)

                              u1 =  - two     * BuxAA(p,s) 
                              u2 =   + AuxAA(p,s)

                              ! val = val + frac12 * Aux3B(s, p)
                              ! val = val + frac12 * Aux3A(p, s)

                              val = val + frac12 * Aux3X(s, p)
                              val = val + frac12 * Aux3X(p, s)

                              u3 =  frac12 * Aux3X(s, p)
                              u4 =  frac12 * Aux3X(p, s)

                              if (rs==5.and.pq==5)then
                                    print*, 'kupsko4', (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(s)) * HNO(p,s)
                                    write(*, '( 4F20.10)') u1, u2, u3, u4
                              end if



                        end if


                        vals = vals -val
                        valt = valt + val

                        NumF = one
                        NumF = one/(one-AuxData%Occ(r)-AuxData%Occ(s))


                        NumH = one
                        if (p==q.and.r==s)then
                              NumH = frac12
                        else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
                              NumH = sqrt(frac12)
                        end if
                        
                        ! if (abs(ACBlock%ASing(rs, pq)).gt.1.d-5)then
                        !       write(*,'(A10, 2I5, A5, 4I5,4F30.16)') 'luzluzprzed', rs, pq, '  |  ', p, q, r, s, ACBlock%ASing(rs, pq), AuxData%Occ(r), AuxData%Occ(s), ACBlock%ASing(rs, pq)+vals
                        !       print*, ''
                        ! end if

                        
                        ACBlock%ASing(rs, pq) = NumF * NumH * ACBlock%ASing(rs, pq) + NumF * NumH*vals

                        ! if (abs(ACBlock%ASing(rs, pq)).gt.1.d-5)then
                        !       write(*,'(A10, 2I5, A5, 4I5,5F30.16)') 'luzluz', rs, pq, '  |  ', p, q, r, s, ACBlock%ASing(rs, pq), NumF * NumH*vals, NumF, NumH, vals
                        !       print*, ''
                        ! end if



                        if (p>q .and. r> s)then
                              pqT = posT(p, q)
                              rsT = posT(r, s)
               !                if (abs(ACBlock%ATrip(rsT, pqT)).gt.1.d-5)then
               ! write(*,'(2I5, A5, 4I5, 5F30.16)') rsT, pqT, '  |  ', p, q, r, s, ACBlock%ATrip(rsT, pqT), valt, numF, NumF * ACBlock%ATrip(rsT, pqT) , NumF  * valt
!                              end if

                              ACBlock%ATrip(rsT, pqT) = NumF * ACBlock%ATrip(rsT, pqT) + NumF  * valt
                              ACBlock%ATripA(rsT, pqT) = NumF * ACBlock%ATripA(rsT, pqT) + NumF  * valt
                              ! if (abs(ACBlock%ATrip(rsT, pqT)).gt.1.d-5)then
                              !       write(*,'(2I5, A5, 4I5, 5F30.16)') rsT, pqT, '  |  ', p, q, r, s, ACBlock%ATrip(rsT, pqT), valt, vals
                              ! end if


                        end if


 !                       if (abs(ACBlock%ASing(rs, pq)).gt.1.d-5)then

!            write(*,'(A10, 2I5, A5, 4I5,6F15.6)') 'sing', rs, pq, '  |  ', p, q, r, s, ACBlock%ASing(rs, pq), NumF, NumH, vals, AuxData%Occ(r),AuxData%Occ(s)
  !                      end if


                  end do i_rowloops
            end do j_colloops
            !                !$omp end parallel do

      end subroutine PPERPA0_THC_AA



      subroutine PPERPA_THC(AuxData, ACBlock, HNO, posS, posT, AuxII, AuxAA, BuxII, BuxAA,Aux3X, Aux3B, ACAlpha, NDim1, NDim2, IndN1, IndN2)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            double precision, dimension(:, :), intent(in) ::HNO
            double precision, dimension(:,:), intent(in) :: AuxII, AuxAA, BuxII, BuxAA,Aux3X, Aux3B
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, intent(in) :: ACAlpha
            integer, intent(in) :: NDim1, NDim2
            integer, dimension(:,:), intent(in) ::IndN1, IndN2

            integer :: p, q, r, s
            integer :: pq, rs, pqT, rsT
            double precision :: val, vals, valt, val1, val0
            double precision :: valA
            double precision :: NumF, NumH
            integer :: lll

            lll = 1


            j_colloops: do pq = 1, NDim2
                  i_rowloops: do rs = 1, NDim1

                        r = IndN1(1, rs)
                        s = IndN1(2, rs)

                        p = IndN2(1, pq)
                        q = IndN2(2, pq)
                        val = zero
                        vals = zero
                        valA = zero


                        if (p==r) then
                              !------------------------------------------------------------------

                              
                              val = val + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s) !P1
                              !valA = valA + (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(s)) * HNO(q,s) !P1

                              val = val + two * (BuxII(q, s) + BuxAA(q, s)) !P2a-1
                              val = val - AuxII(q, s) - AuxAA(q, s)   !P2a-1


                              if (AuxData%IndAux(q)==0) val = val  - (BuxII(q,s) + BuxAA(q, s)) + frac12 * (AuxII(q,s) + AuxAA(q, s))
                              if (AuxData%IndAux(q)==1) val = val  - AuxData%Occ(q) * (BuxII(q, s) - frac12 * AuxII(q,s))


                              if (AuxData%IndAux(s)==0) val = val  - (BuxII(q,s) + BuxAA(q, s)) + frac12 * (AuxII(q,s) + AuxAA(q, s))
                              if (AuxData%IndAux(s)==1) val = val  - AuxData%Occ(s) * (BuxII(q, s) - frac12 * AuxII(q,s))


                              if (AuxData%IndAux(p)==0)then
                                    val = val + AuxData%Occ(p)*(AuxII(q,s) + AuxAA(q,s) - two * BuxII(q, s) -two * BuxAA(q, s))  ! p3B-1, P3c-3,4
                              else
                                    val = val + AuxData%Occ(p)*(AuxII(q,s)- two * BuxII(q, s))
                              end if

                              ! if (AuxData%IndAux(q)==1) val = val - frac12 * Aux3A(q, s)
                              ! if (AuxData%IndAux(s)==1) val = val - frac12 * Aux3B(s, q)

                              if (AuxData%IndAux(q)==1) val = val - frac12 * Aux3X(q, s)
                              if (AuxData%IndAux(s)==1) val = val - frac12 * Aux3X(s, q)

                              
                        end if


                        if (q==s)then


                              val = val+ (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(r)) * HNO(p,r) !p1-2
                              
                              
                              val = val + two * (BuxII(p,r) + BuxAA(p,r)) !p2A-2
                              

                              val = val - AuxII(p,r) - AuxAA(p,r) !p2A-2
                              
                              if (AuxData%IndAux(p)==0) val = val  - (BuxII(p,r) + BuxAA(p,r)) + frac12 * (AuxII(p,r) + AuxAA(p,r))
                        

                              if (AuxData%IndAux(p)==1) val = val  - AuxData%Occ(p) * (BuxII(p,r) - frac12 * AuxII(p,r))
                        

                              if (AuxData%IndAux(r)==0) val = val  - (BuxII(p,r) + BuxAA(p,r))   + frac12 * (AuxII(p,r) + AuxAA(p,r))
                        

                              if (AuxData%IndAux(r)==1) val = val  - AuxData%Occ(r) * (BuxII(p,r) - frac12 * AuxII(p,r))
                        
                              if (AuxData%IndAux(q)==0)then
                                    val = val + AuxData%Occ(q)*(AuxII(p,r) + AuxAA(p,r) - two* BuxII(p, r) - two * BuxAA(p, r)) ! p3B-2 p3C-1,2
                              else
                                    val = val + AuxData%Occ(q)*(AuxII(p,r)- two *BuxII(p, r))
                              end if
                              
                              ! if (AuxData%IndAux(p)==1) val = val - frac12 * Aux3A(p, r)
                              ! if (AuxData%IndAux(r)==1) val = val - frac12 * Aux3B(r, p)

                              if (AuxData%IndAux(p)==1) val = val - frac12 * Aux3X(p, r)
                              if (AuxData%IndAux(r)==1) val = val - frac12 * Aux3X(r, p)

                        end if



                        vals = val
                        valt = val
                        val = zero

                        if (p==s)then

                              val = val - (one- AuxData%Occ(p)-frac12 * AuxData%Occ(q) - frac12*AuxData%Occ(r)) * HNO(q,r) !P1-3


                              val = val - two     * (BuxII(q,r) + BuxAA(q,r )) !p2a-3

                              
                              val = val + (AuxII(q,r) + AuxAA(q,r))!P2a3

                              ! if (AuxData%IndAux(r)==1) val = val + frac12 * Aux3B(r, q)
                              ! if (AuxData%IndAux(q)==1) val = val + frac12 * Aux3A(q, r)

                              if (AuxData%IndAux(r)==1) val = val + frac12 * Aux3X(r, q)
                              if (AuxData%IndAux(q)==1) val = val + frac12 * Aux3X(q, r)

                              if (AuxData%IndAux(q)==0) val = val  + (BuxII(q,r) + BuxAA(q,r)) - frac12 * (AuxII(q,r) + AuxAA(q,r))

                              if (AuxData%IndAux(q)==1) val = val  + AuxData%Occ(q) * frac12*BuxII(q,r) 
                              if (AuxData%IndAux(q)==1) val = val  + AuxData%Occ(q) * (frac12 * BuxII(q,r) - frac12 * AuxII(q,r))
                              
                              if (AuxData%IndAux(r)==0) val = val  + (BuxII(q,r) + BuxAA(q,r))   - frac12 * (AuxII(q,r) + AuxAA(q,r))

                              if (AuxData%IndAux(r)==1) val = val  + AuxData%Occ(r) * frac12* BuxII(q,r)

                              if (AuxData%IndAux(r)==1) val = val  + AuxData%Occ(r) * (frac12 * BuxII(q,r) - frac12 * AuxII(q,r))
                              

                              if (AuxData%IndAux(p)==0)then
                                    val = val - AuxData%Occ(p)*(AuxII(q,r) + AuxAA(q,r) - two * BuxII(q, r) - two * BuxAA(q, r)) !p3B-4 p3C-7,8
                              else                                        
                                    val = val - AuxData%Occ(p)*(AuxII(q,r) - two * BuxII(q,r))
                              end if


                        end if


                        if (q==r)then
                              val = val - two     * (BuxII(p,s) + BuxAA(p,s)) !P2a-4
                              val = val + (AuxII(p,s) + AuxAA(p,s))!P2a-4
                              val = val - (one- AuxData%Occ(q)-frac12 * AuxData%Occ(p) - frac12*AuxData%Occ(s)) * HNO(p,s) !P1-4

                              !if (AuxData%IndAux(s)==1) val = val + frac12 * Aux3B(s, p)
                              !if (AuxData%IndAux(p)==1) val = val + frac12 * Aux3A(p, s)

                              if (AuxData%IndAux(s)==1) val = val + frac12 * Aux3X(s, p)
                              if (AuxData%IndAux(p)==1) val = val + frac12 * Aux3X(p, s)

                              if (AuxData%IndAux(p)==0) val = val  + (BuxII(p,s) + BuxAA(p,s)) - frac12 * (AuxII(p,s) + AuxAA(p,s))
                              if (AuxData%IndAux(p)==1) val = val  + AuxData%Occ(p) * (BuxII(p,s) - frac12 * AuxII(p,s))

                              if (AuxData%IndAux(s)==0) val = val  + (BuxII(p,s) + BuxAA(p,s))   - frac12 * (AuxII(p,s) + AuxAA(p,s))
                              if (AuxData%IndAux(s)==1) val = val  + AuxData%Occ(s) * (BuxII(p,s) - frac12 * AuxII(p,s))

                              
                              if (AuxData%IndAux(q)==0)then
                                    val = val - AuxData%Occ(q)*(AuxII(p,s) + AuxAA(p,s) - two * BuxII(p, s) - two * BuxAA(p, s)) !p3B-4 p3C-7,8
                              else                                        
                                    val = val - AuxData%Occ(q)*(AuxII(p,s) - two * BuxII(p, s))
                              end if

                        end if

                        vals = vals -val
                        valt = valt + val

                        NumF = one
                        !          NumF = one/(one-AuxData%Occ(r)-AuxData%Occ(s))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

                        NumH = one
                        if (p==q.and.r==s)then
                              NumH = frac12
                        else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
                              NumH = sqrt(frac12)
                        end if
                        

                              ACBlock%ASing(rs, pq) = NumF * NumH * ACBlock%ASing(rs, pq) + NumF * NumH*vals


                        if (p>q .and. r> s)then
                              pqT = posT(p, q)
                              rsT = posT(r, s)
                              ! if (rsT==17.and.pqT==1)then
                              !       write(*,'(A10, 2I5, A5, 4I5,2F30.16)') 'przedT', rsT, pqT, '  |  ', p, q, r, s, ACBlock%ATripA(rsT, pqT), valt
                              ! end if


                              ACBlock%ATrip(rsT, pqT) = NumF * ACBlock%ATrip(rsT, pqT) + NumF  * valt
                              ACBlock%ATripA(rsT, pqT) = NumF * ACBlock%ATripA(rsT, pqT) + NumF  * valt

                              ! if (abs(ACBlock%ATrip(rsT, pqT)).gt.1.d-5)then
                              !       write(*,'(A10, 2I5, A5, 4I5,2F30.16)') 'zulzul0', rsT, pqT, '  |  ', p, q, r, s, ACBlock%ATrip(rsT, pqT)
                              !       print*, ''
                              ! end if
                              
                              ! if (abs(ACBlock%ATripA(rsT, pqT)).gt.1.d-5)then
                              !       write(*,'(A10, 2I5, A5, 4I5,2F30.16)') 'zulzulT', rsT, pqT, '  |  ', p, q, r, s, ACBlock%ATripA(rsT, pqT)
                              ! end if
                              
                        end if


                        ! if (abs(ACBlock%ASing(rs, pq)).gt.1.d-4)then
                        !       write(*,'(2I5, A5, 4I5,2F30.16)') rs, pq, '  |  ', p, q, r, s, ACBlock%ASing(rs, pq)
                        ! end if
                        

                  end do i_rowloops
            end do j_colloops
            !                !$omp end parallel do                                                    


      end subroutine PPERPA_THC


      subroutine eigs(AuxData, ACBlock, x)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block), intent(inout) :: ACBlock
            integer, intent(in) :: x
            double precision, dimension(:, :), allocatable :: Swork
            integer :: p, q, pq, i, j, ii
            integer, parameter :: A0oo = 1, A0vv = 2, A0aa = 3


            if (x == A0oo.or.x == A0vv)then
                  if (x == A0oo)then
                        ACBlock%vplusS = 0
                        ACBlock%vplusT = 0
                        ACBlock%vplusTA = 0
                  else
                        ACBlock%vplusS = 1
                        ACBlock%vplusT = 1
                        ACBlock%vplusTA = 1
                  end if

                  ! print*, 'Wartosci wlasne bloku', ACBlock%name, ACBlock%NDimS

                  ! do i = 1, ACBlock%NDimS
                  !       if (abs(ACBlock%ASing(i, 1)).gt.1.d-5)then
                  !             print*, ACBlock%ASing(i, 1)
                  !       end if
                  ! end do

                  ! print*, 'Wartosci wlasne bloku', ACBlock%name, ACBlock%NDimT
                  
                  ! do i = 1, ACBlock%NDimT
                  !       if (abs(ACBlock%ATrip(i, 1)).gt.1.d-5)then
                  !             print*, ACBlock%ATrip(i, 1)
                  !       end if
                  ! end do

                  ! print*, 'Wartosci wlasne bloku T', ACBlock%name, ACBlock%NDimT
                                    
                  ! do i = 1, ACBlock%NDimT
                  !       if (abs(ACBlock%ATripA(i, 1)).gt.1.d-5)then
                  !             print*, ACBlock%ATripA(i, 1)
                  !       end if
                  ! end do

            else
                  
                  allocate(Swork(ACBlock%NDimS, ACBlock%NDimS))

                  swork = zero
                  do i = 1, ACBlock%NDimS
                        p = ACBlock%IndNS(1, i)
                        q = ACBlock%IndNS(2, i)
                        Swork(i, i) = One / (One-AuxData%Occ(p)-AuxData%Occ(q))
                  end do


                  if (ACBlock%NDimS>0)then
                        call NonSymEigBlockTHC(ACBlock%EigS, ACBlock%ASing, Swork, ACBlock%vplusS, AuxData%IndAux, ACBlock%NDimS)
                  end if
                        
                  deallocate(Swork)


                  allocate(Swork(ACBlock%NDimT, ACBlock%NDimT))
                                    
                  swork = zero
                  do i = 1, ACBlock%NDimT
                        p = ACBlock%IndNT(1, i)
                        q = ACBlock%IndNT(2, i)
                        Swork(i, i) = One / (One-AuxData%Occ(p)-AuxData%Occ(q))
                  end do

                  if (ACBlock%NDimT>0)then
                        call NonSymEigBlockTHC(ACBlock%EigT, ACBlock%ATrip, Swork, ACBlock%vplusT, AuxData%IndAux, ACBlock%NDimT)
                  end if

                  swork = zero
                  do i = 1, ACBlock%NDimT
                        p = ACBlock%IndNT(1, i)
                        q = ACBlock%IndNT(2, i)
                        Swork(i, i) = One / (One-AuxData%Occ(p)-AuxData%Occ(q))
                  end do
                  if (ACBlock%NDimT>0)then                  
                        call NonSymEigBlockTHC(ACBlock%EigTA, ACBlock%ATripA, Swork, ACBlock%vplusTA, AuxData%IndAux, ACBlock%NDimT)
                  end if

                  if (ACBlock%type == A0aa) then
                        call divide_AABlocks(ACBlock, 0)
                        call divide_AABlocks(ACBlock, 1)
                        call divide_AABlocks(ACBlock, 2)
                  end if
                  
            end if

      end subroutine eigs



      subroutine divide_AABlocks(ACB, st)
            
            type(TAC0Block), intent(inout) :: ACB
            integer, intent(in) :: st

            integer :: i, j, Npl, Nmn, imn, ipl
            integer, parameter :: sing = 0, trip = 1, tripA = 2
            integer, dimension(:), allocatable :: plus_idx, minus_idx
            

            if (st == sing)then

                  Npl = 0
                  Nmn = 0
                  do i = 1, ACB%NDimS
                        if (ACB%vplusS(i)==0) Nmn = Nmn + 1
                        if (ACB%vplusS(i)==1) Npl = Npl + 1	
                  end do


                  allocate(ACB%MiniBlocks(2))
                  ACB%MiniBlocks(1)%NDimS = Nmn
                  ACB%MiniBlocks(2)%NDimS = Npl
                  
                  allocate(ACB%MiniBlocks(1)%MiniAVS(ACB%NDimS, Nmn))
                  allocate(ACB%MiniBlocks(2)%MiniAVS(ACB%NDimS, Npl))

                  allocate(ACB%MiniBlocks(1)%EigS(Nmn))
                  allocate(ACB%MiniBlocks(2)%EigS(Npl))

                  Npl	= 0
                  Nmn	= 0
                  do i = 1, ACB%NDimS

                        if (ACB%vplusS(i)==0) then
                              Nmn = Nmn + 1
                              ACB%MiniBlocks(1)%MiniAVS(:, Nmn) = ACB%ASing(:,i)
                              ACB%MiniBlocks(1)%EigS(Nmn) = ACB%EigS(i)
                        else if (ACB%vplusS(i)==1) then
                              Npl = Npl + 1
                              ACB%MiniBlocks(2)%MiniAVS(:, Npl) = ACB%ASing(:,i)
                              ACB%MiniBlocks(2)%EigS(Npl) = ACB%EigS(i)
                        end if
                  end do

                  deallocate(ACB%ASing)
                  deallocate(ACB%EigS)

            else if (st ==trip)then

                  Npl = 0
                  Nmn = 0
                  do i = 1, ACB%NDimt
                        if (ACB%vplusT(i)==0) Nmn = Nmn + 1
                        if (ACB%vplusT(i)==1) Npl = Npl + 1   
                  end do

                  ACB%MiniBlocks(1)%NDimT = Nmn
                  ACB%MiniBlocks(2)%NDimT = Npl

                  allocate(ACB%MiniBlocks(1)%MiniAVT(ACB%NDimt, Nmn))
                  allocate(ACB%MiniBlocks(2)%MiniAVT(ACB%NDimt, Npl))

                  allocate(ACB%MiniBlocks(1)%EigT(Nmn))
                  allocate(ACB%MiniBlocks(2)%EigT(Npl))

                  Npl   = 0
                  Nmn   = 0
                  do i = 1, ACB%NDimT

                        if (ACB%vplusT(i)==0) then
                              Nmn = Nmn + 1
                              ACB%MiniBlocks(1)%MiniAVT(:, Nmn) = ACB%ATrip(:,i)
                              ACB%MiniBlocks(1)%EigT(Nmn) = ACB%EigT(i)
                        else if (ACB%vplusT(i)==1) then
                              Npl = Npl + 1
                              ACB%MiniBlocks(2)%MiniAVT(:, Npl) = ACB%ATrip(:,i)
                              ACB%MiniBlocks(2)%EigT(Npl) = ACB%EigT(i)
                        end if
                  end do

                  deallocate(ACB%ATrip)
                  deallocate(ACB%EigT)


            else if (st ==tripA)then

                  Npl = 0
                  Nmn = 0
                  do i = 1, ACB%NDimt
                        if (ACB%vplusTA(i)==0) Nmn = Nmn + 1
                        if (ACB%vplusTA(i)==1) Npl = Npl + 1   
                  end do
                  ACB%MiniBlocks(1)%NDimTA = Nmn
                  ACB%MiniBlocks(2)%NDimTA = Npl


                  allocate(ACB%MiniBlocks(1)%MiniAVTA(ACB%NDimt, Nmn))
                  allocate(ACB%MiniBlocks(2)%MiniAVTA(ACB%NDimt, Npl))

                  allocate(ACB%MiniBlocks(1)%EigTA(Nmn))
                  allocate(ACB%MiniBlocks(2)%EigTA(Npl))

                  Npl   = 0
                  Nmn   = 0
                  do i = 1, ACB%NDimT

                        if (ACB%vplusTA(i)==0) then
                              Nmn = Nmn + 1
                              ACB%MiniBlocks(1)%MiniAVTA(:, Nmn) = ACB%ATripA(:,i)
                              ACB%MiniBlocks(1)%EigTA(Nmn) = ACB%EigTA(i)
                        else if (ACB%vplusTA(i)==1) then
                              Npl = Npl + 1
                              ACB%MiniBlocks(2)%MiniAVTA(:, Npl) = ACB%ATripA(:,i)
                              ACB%MiniBlocks(2)%EigTA(Npl) = ACB%EigTA(i)
                        end if
                  end do

                  deallocate(ACB%ATripA)
                  deallocate(ACB%EigTA)

            else
                  print*, 'this spin is not implemented for AA block'
            end if
            
      end subroutine divide_AABlocks


      subroutine NonSymEigBlockTHC(wr, A, S, v_plus, IndAux, n)

            double precision, dimension(:), intent(out), contiguous      :: wr
            double precision, dimension(:, :), intent(inout), contiguous ::  S
            double precision, dimension(:),allocatable       :: wi
            double precision, dimension(:, :), allocatable   :: vr, vl
            double precision, dimension(:, :), intent(inout), contiguous :: A
            integer, dimension(:), intent(out), contiguous :: v_plus            
            integer, dimension(:),intent(in) :: IndAux
            integer, intent(in) :: n
            double precision :: shift, max_hh, min_pp                
            double precision, dimension(1) :: work0
            double precision, dimension(:), allocatable :: work
            integer :: lwork, info
            integer :: i, j, k
            double precision :: dd, dd_comp
            type (tclock) :: timer
            double precision, dimension(:), allocatable :: wr_plus, wr_minus
            integer, dimension(:), allocatable :: dy_plus, dy_minus
            integer :: count, countj_plus, countj_minus, ss
            double precision, dimension(:), allocatable :: S_diag


            external :: dgeev

            !n = size(A, dim=1)

            allocate(vr(n,n))
            allocate(S_diag(n))

            allocate(dy_plus(n))
            allocate(dy_minus(n))
            allocate(wr_plus(n))
            allocate(wi(n))
            allocate(wr_minus(n))

            allocate(vl(1,1))
            vl = zero
            vr = zero



            S_diag = zero
            do i = 1, n
                  S(i, i) = One / S(i, i)
                  S_diag(i) = S(i, i)
            end do
            


            call clock_start(timer)
            lwork = -1
            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work0, lwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))

            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work, lwork, info)
            call ximsg('info',info, mdebug)

            if (info /= 0) then
                  call mmsg("Nonsymmetric matrix eigendecompositino failed with info=", merror)
                  error stop
            end if
            


            !call tmsg('TIME FOR DIAG', timer, tdebug)

            call clock_start(timer)
            v_plus = 0

            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        if (wi(i).gt.0.d+0)then
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)
                              call ddot_norm(vr(:, i+1),S,  vr(:, i+1), n, dd_comp)

                              if (msgthr == 0)then
                                    print*, '>0 dd', i, i+1, dd
                                    print*, '>0 dd_comp', dd_comp
                              end if
                        else
                              call ddot_norm(vr(:, i-1),S,  vr(:, i-1), n, dd)
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd_comp)
                              if (msgthr == 0)then                                                                  
                                    print*, '<0 dd', i-1, i, dd
                                    print*, '<0 dd_comp',   dd_comp
                              end if

                        end if
                        dd = dd + dd_comp
                        v_plus(i)=2
                  else
                        call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)

                        if (dd.gt.zero)then
                              v_plus(i) = 1
                        end if

                  end if
                  if (msgthr == 0)then                                                                  
                        if(abs(wi(i)).gt.1.d-8)then
                              print*, 'i-norm', i, dd
                        end if
                  end if
            end do


            call clock_start(timer)

            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'COMPLEX EIGENVALUES', i, wr(i), wi(i)
                        !stop                        
                  end if
            end do

            wr_plus = wr
            wr_minus = wr
            do i = 1, n
                  dy_plus(i)=i
            end do
            call dsort(wr_plus, dy_plus, n)

            shift = zero
            do i = 1, n
                  if (v_plus(dy_plus(i))== 0 )then
                        max_hh = wr_plus(i)
                  end if
                  if (v_plus(dy_plus(i))== 1 )then
                        min_pp = wr_plus(i)
                        if (i == 1) then
                              shift = min_pp/two
                        else
                              shift = (abs(min_pp) + abs(max_hh))/two

                              call xmsg('max_hh', max_hh, mdebug)
                        end if
                        call xmsg('min_pp' , wr_plus(i), mdebug)
                        call xmsg('shift', shift, mdebug)
                        exit
                  else

                  end if
            end do

            
            wr_plus = zero
            wr_minus = zero
            dy_plus = 0
            dy_minus = 0
            countj_plus=1
            countj_minus=1

            do i = 1, n
                  if (v_plus(i) == 1) then
                        wr_plus(countj_plus) = wr(i)
                        dy_plus(countj_plus) = i
                        countj_plus = countj_plus + 1
                  else if (v_plus(i) == 0) then
                        wr_minus(countj_minus) = wr(i)
                        dy_minus(countj_minus) = i
                        countj_minus = countj_minus + 1

                  end if
            end do
            countj_plus = countj_plus-1
            countj_minus = countj_minus-1

            call dsort(wr_plus(1:countj_plus), dy_plus(1:countj_plus), countj_plus)
            call dsort(wr_minus(1:countj_minus), dy_minus(1:countj_minus), countj_minus)


            call orthogonalize_degen(n, countj_plus, wr_plus, vr, S_diag, dy_plus, 1)
            call orthogonalize_degen(n, countj_minus, wr_minus, vr, S_diag, dy_minus, 0)


            A = vr

      end subroutine NonSymEigBlockTHC


end module ppac0_thc
