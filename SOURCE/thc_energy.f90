
module thc_energy

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals                                                                                                                                                                                                
      use types
      use lin
      use clock

      use math_constants
      use ppac_types
      use print_utils
      use omp_lib

      implicit none

      integer, parameter, private :: tdebug = 0
      integer, parameter, private :: mdebug = 0
      integer, parameter, private :: mdebugg = 1
      integer, parameter, private :: mnormal = 50
      integer, parameter, private :: mverbose = 10
      integer, parameter, private :: merror = 100

      integer, parameter, private :: msgthr = 0
      double precision, parameter, private :: ttoeV = 27.211399d+0



contains

      subroutine THC_energy_loop(THCData, ACB, Flags, posS, posT, AuxData)
            
            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use sys_definitions
            use gammcor_integrals

            type(TTHCData), intent(in) :: THCData
            type(TAC0Block), dimension(:), intent(inout) :: ACB
            type(FlagsData), intent(in) :: Flags
            integer, dimension(:,:), intent(in) :: posS, posT
            type(TACppData), intent(in) :: AuxData
            integer :: Batchdim
            integer, dimension(6) :: TabChol


            type (tclock) :: timer, timer0, timer90, timerx
            double precision :: timer10, timer11, timer12, timer13, timer14, timer15, timer16, timer17
            double precision :: timer18, timer19, timer33, timer34, timer35, timer36, timer21, timer37
            double precision :: timer38, timer39, timer14a, timer18a, timer13x, timer13y, timer13z

            
            double precision, dimension(:), allocatable :: WorkAO, WorkAA, V_axby, V_bxay, workAOT
            double precision, dimension(:,:), allocatable :: AZaaoo, AZvaoo, AZvaaa, AZvaao
            double precision, dimension(:,:), allocatable :: AZminiaa, AZminiao, AZminiao2, AZmini
            double precision, dimension(:,:), allocatable :: AZminiaaao, AZaaao


            double precision, dimension(:,:), allocatable :: AZaaooT, AZvaooT, AZvaaaT, AZvaaoT
            double precision, dimension(:,:), allocatable :: AZaaooTA
            double precision, dimension(:,:), allocatable :: AZminiaaT, AZminiaoT, AZminiao2T, AZminiT
            double precision, dimension(:,:), allocatable :: AZminiaaaoT, AZaaaoT
            double precision, dimension(:,:), allocatable :: AZminiaaaoTA, AZaaaoTA
            double precision, dimension(:,:), allocatable :: AZvaooTA, AZvaaaTA, AZvaaoTA
            double precision, dimension(:,:), allocatable :: AZminiaaTA, AZminiaoTA, AZminiao2TA, AZminiTA

            integer :: NA, NI, NV, NIA, nbasis
            integer ::  p, q, r, s, pq1, rs1, pq2, rs2, rs
            integer :: a, b, c, d
            integer :: x0, x1, y0, y1
            integer :: c0, c1, d0, d1
            integer :: t, u
            integer :: ii
            double precision :: val, val1, val2
            type(TLoop), dimension(6) :: loops


            integer :: NBatchI, NBatchA, NBatchV, NBatchIA
            integer :: batch
            integer :: b0VI, b1VI, b0VO, b1VO

            
            integer :: b0II, b1II, b0IA, b1IA, b0AA, b1AA
            integer :: b0IV, b1IV, b0VA, b1VA, b0VV, b1VV
            double precision, allocatable :: Rkax(:,:), Rkby(:,:)
            double precision, allocatable :: XXga(:,:)
            double precision, dimension(:,:,:), allocatable :: RII, RIA, RAI, RAA, RIV, RVA, RVV
            double precision, dimension(:,:,:), allocatable :: RVI, RVO, ROO
            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            integer :: pq, pqT, rsT, x, y
            integer :: ndim_virt, k, nva, kk
            integer :: nblock, bi
            integer :: method
            double precision :: E_corr
            integer, parameter :: AC0 = 0, AC= 1
            double precision, parameter :: ACAlpha = One
            double precision :: EVVoo, EVVao, EVVaa, EVAoo, EVAao, EVAaa, EAAoo, EAAao, EVAao2
            double precision :: EVVooT, EVVaoT, EVVaaT, EVAooT, EVAaoT, EVAaaT, EAAooT, EAAaoT, EVAao2T
            double precision :: EAAooTA, EAAaoTA
            double precision :: EVVooTA, EVVaoTA, EVVaaTA, EVAooTA, EVAaaTA, EVAaoTA, EVAao2TA
            integer, parameter  :: A0oo=1, A0vv=2, A0aa=3, A0va=4, A0ao=5
            integer, parameter :: A1vvoo=6, A1vvao=7, A1vvaa=8, A1vaoo=9
            integer, parameter :: A1vaao=10, A1vaaa=11, A1aaoo=12, A1aaao=13
            integer, parameter :: A1vaao2 = 14
            character(len=8), dimension(14), parameter :: BlockName = (/ &
                  'A0oo    ', 'A0vv    ', 'A0aa    ', 'A0va    ', 'A0ao    ', &
                  'A1vvoo  ', 'A1vvao  ', 'A1vvaa  ', 'A1vaoo  ', 'A1vaao  ', &
                  'A1vaaa  ', 'A1aaoo  ', 'A1aaao  ', 'A1vaao2 ' /)
            integer :: i, j, l0, l1, n, l, km, q1, q2
            integer :: ds1, ds2, ds3, ds4

            double precision :: this_mem
            double precision :: togb

            togb = 8.d+0/(1024.0)**3

            associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga, ExternalOrdering=>THCData%ExternalOrdering, NChol=>THCData%NChol, NTHC=>THCData%NTHC)

            BatchDim = AuxData%Batchdim
            call print_section('ppAC0 energy evaluation')
            EVVoo = zero
            EVVao = zero
            EVVaa = zero
            EVAoo = zero
            EVAao = zero
            EVAao2 = zero
            EVAaa = zero
            EAAoo = zero
            EAAao = zero

            EVVooT = zero
            EVVaoT = zero
            EVVaaT = zero
            EVAooT = zero
            EVAaoT = zero
            EVAao2T = zero
            EVAaaT = zero
            EAAooT = zero
            EAAaoT = zero
            EAAooTA = zero
            EAAaoTA = zero
            EVVooTA = zero
            EVVaoTA = zero
            EVVaaTA = zero
            EVAooTA = zero
            EVAaaTA = zero
            EVAaoTA = zero
            EVAao2TA = zero

            
            nbasis = AuxData%NBasis
            NI = AuxData%NI
            NIA = AuxData%NIA
            NV = AuxData%NV
            NA = AuxData%NA

            TabChol = 0
            TabChol(4) = 1
            TabChol(5) = 1
            TabChol(3) = 1
            TabChol(2) = 1


            do i = 1, 6
                  do j = 1, i
                        loops(i)%run_main = .false.
                        loops(i)%run_inner(j) = .false.
                  end do
            end do

            loops(4)%run_main = .true.
            ! loops(5)%run_main = .true.
            loops(4)%run_inner(4) = .true.
            loops(4)%run_inner(2) = .true.
            ! loops(5)%run_inner(2) = .true.

            loops(5)%run_main = .true.
            loops(5)%run_inner(4) = .true.
            loops(5)%run_main = .true.
            loops(5)%run_inner(3) = .true.

            loops(2)%run_main = .true.
            loops(2)%run_inner(2) = .true.
            

            allocate(Rkax(NChol, NBasis))
            allocate(Rkby(NChol, NBasis))
            allocate(XXga(NTHC,NBasis))

#ifdef DEBUG
            this_mem = (NChol*NBasis)*togb
            call print_memory('rkax', this_mem)

            this_mem = (NChol*NBasis)*togb
            call print_memory('rkby', this_mem)

            this_mem = (NTHC*NBasis)*togb
            call print_memory('xxga', this_mem)
#endif


            allocate(V_axby(max(NA, NI, NV, NIA)*max(NA, NI, NV, NIA)))
            allocate(V_bxay(max(NA, NI, NV, NIA)*max(NA, NI, NV, NIA)))

#ifdef DEBUG
            this_mem = (NV*NV*2)*togb
            call print_memory('V_axby', this_mem)

            print*, 'naninv', na, ni, nv
            print*, max(na, ni, nv)
            print*, max(NA, NI, NV)*max(NA, NI, NV)
#endif


            allocate(WorkAO(ACB(A0ao)%NDimS))
            allocate(WorkAOT(ACB(A0ao)%NDimS))
            allocate(WorkAA(ACB(A0aa)%NDimS))
            
#ifdef DEBUG
            this_mem = (ACB(A0ao)%NDimS)*togb
            call print_memory('workao', this_mem)
#endif



            allocate(AZaaoo(ACB(A0aa)%MiniBlocks(2)%NdimS, ACB(A0oo)%NdimS))
            call CalcMem(AZaaoo, 'AZaaoo')
            
            allocate(AZvaoo(ACB(A0oo)%NdimS, ACB(A0va)%NdimS))
            allocate(AZvaaa(ACB(A0aa)%MiniBlocks(1)%NdimS, ACB(A0va)%NdimS))
            allocate(AZvaao(ACB(A0ao)%NdimS, ACB(A0va)%NdimS))
            allocate(AZmini(ACB(A0oo)%NdimS, NA))
            allocate(AZminiaa(ACB(A0aa)%MiniBlocks(1)%NdimS, NA))
            allocate(AZminiao(NA, ACB(A0va)%NdimS))
            allocate(AZminiao2(ACB(A0ao)%NdimS, NA))           
            allocate(AZminiaaao(ACB(A0aa)%MiniBlocks(2)%NdimS, NA))
            allocate(AZaaao(ACB(A0aa)%MiniBlocks(2)%NdimS, ACB(A0ao)%NdimS))


            allocate(AZaaooT(ACB(A0aa)%MiniBlocks(2)%NdimT, ACB(A0oo)%NdimT))
            allocate(AZaaooTA(ACB(A0aa)%MiniBlocks(2)%NdimTA, ACB(A0oo)%NdimT))
            call CalcMem(AZaaooTA, 'AZaaooTA')
            
            allocate(AZvaooT(ACB(A0oo)%NdimT, ACB(A0va)%NdimS))
            allocate(AZvaaaT(ACB(A0aa)%MiniBlocks(1)%NdimT, ACB(A0va)%NdimS))
            allocate(AZvaaoT(ACB(A0ao)%NdimS, ACB(A0va)%NdimS))
            allocate(AZminiT(ACB(A0oo)%NdimT, NA))
            allocate(AZminiaaT(ACB(A0aa)%MiniBlocks(1)%NdimT, NA))
            allocate(AZminiaoT(NA, ACB(A0va)%NdimS))
            allocate(AZminiao2T(ACB(A0ao)%NdimS, NA))           
            allocate(AZminiaaaoT(ACB(A0aa)%MiniBlocks(2)%NdimT, NA))
            allocate(AZaaaoT(ACB(A0aa)%MiniBlocks(2)%NdimT, ACB(A0ao)%NdimS))
            allocate(AZminiaaaoTA(ACB(A0aa)%MiniBlocks(2)%NdimTA, NA))

            allocate(AZaaaoTA(ACB(A0aa)%MiniBlocks(2)%NdimTA, ACB(A0ao)%NdimS))
            allocate(AZvaooTA(ACB(A0oo)%NdimT, ACB(A0va)%NdimS))
            allocate(AZvaaaTA(ACB(A0aa)%MiniBlocks(1)%NdimTA, ACB(A0va)%NdimS))
            allocate(AZvaaoTA(ACB(A0ao)%NdimS, ACB(A0va)%NdimS))

            allocate(AZminiTA(ACB(A0oo)%NdimT, NA))
            allocate(AZminiaaTA(ACB(A0aa)%MiniBlocks(1)%NdimTA, NA))
            allocate(AZminiaoTA(NA, ACB(A0va)%NdimS))
            allocate(AZminiao2TA(ACB(A0ao)%NdimS, NA))



            call CalcMem(AZminiaaaoTA, 'AZminiaaaoTA')
            call CalcMem(AZaaaoTA, 'AZaaaoTA')
            call CalcMem(AZvaooTA, 'AZvaooTA')
            call CalcMem(AZvaaaTA, 'AZvaaaTA')
            call CalcMem(AZvaaoTA, 'AZvaaoTA')
            call CalcMem(AZminiTA, 'AZminiTA')
            call CalcMem(AZminiaaTA, 'AZminiaaTA')
            call CalcMem(AZminiaoTA, 'AZminiaoTA')
            call CalcMem(AZminiao2TA, 'AZminiao2TA')

            call CalcMem(AZvaoo, 'AZvaoo')
            call CalcMem(AZvaaa, 'AZvaaa')
            call CalcMem(AZvaao, 'AZvaao')
            call CalcMem(AZmini, 'AZmini')
            call CalcMem(AZminiaa, 'AZminiaa')
            call CalcMem(AZminiao, 'AZminiao')
            call CalcMem(AZminiao2, 'AZminiao2')
            call CalcMem(AZminiaaao, 'AZminiaaao')
            call CalcMem(AZaaao, 'AZaaao')

            call CalcMem(AZaaooT, 'AZaaooT')
            call CalcMem(AZvaooT, 'AZvaooT')
            call CalcMem(AZvaaaT, 'AZvaaaT')
            call CalcMem(AZvaaoT, 'AZvaaoT')
            call CalcMem(AZminiT, 'AZminiT')
            call CalcMem(AZminiaaT, 'AZminiaaT')
            call CalcMem(AZminiaoT, 'AZminiaoT')
            call CalcMem(AZminiao2T, 'AZminiao2T')
            call CalcMem(AZminiaaaoT, 'AZminiaaaoT')
            call CalcMem(AZaaaoT, 'AZaaaoT')

            
            V_axby = zero
            V_bxay = zero

            WorkAO = zero
            WorkAOT = zero
            WorkAA = zero

            AZaaoo = zero
            AZaaooT = zero
            AZaaooTA = zero
            AZvaooTA = zero
            AZvaaaTA = zero
            AZvaaoTA = zero
            AZminiTA = zero
            AZminiaaTA = zero
            AZminiaoTA = zero
            AZminiao2TA = zero
            
            AZvaoo = zero
            AZvaaa = zero

            AZvaao = zero

            AZmini = zero
            AZminiaa = zero

            AZminiao = zero
            AZminiao2 = zero

            AZvaooT = zero
            AZvaaaT = zero

            AZvaaoT = zero

            AZminiT = zero
            AZminiaaT = zero

            AZminiaoT = zero
            AZminiao2T = zero


            
            
            call clock_start(timer)
            !            if (ACB(A0aa)%MiniBlocks(2)%NdimS.gt.0.and. ACB(A0oo)%NdimS.gt.0)then
            if (size_not_zero(AZaaoo, ACB(A0aa)%MiniBlocks(2)%MiniAVS, ACB(A1aaoo)%ASing))then
                  call real_atb(AZaaoo, ACB(A0aa)%MiniBlocks(2)%MiniAVS, ACB(A1aaoo)%ASing)                  
            end if
            call tmsg('TIME FOR AZ AAOO ', timer, tdebug)

#ifdef DEBUG
            print*, 'size(AZaaooT, dim=1)', size(AZaaooT, dim=1)
            print*, 'size(AZaaooT, dim=2)', size(AZaaooT, dim=2)
            
            print*, 'size(ACB(A1aaoo)%ATrip, dim=1)', size(ACB(A1aaoo)%ATrip, dim=1)
            print*, 'size(ACB(A1aaoo)%ATrip, dim=2)', size(ACB(A1aaoo)%ATrip, dim=2)

            print*, 'size(ACB(A0aa)%MiniBlocks(2)%MiniAVT, dim=1)', size(ACB(A0aa)%MiniBlocks(2)%MiniAVT, dim=1)
            print*, 'size(ACB(A0aa)%MiniBlocks(2)%MiniAVT, dim=2)', size(ACB(A0aa)%MiniBlocks(2)%MiniAVT, dim=2)

            
            
            print*, 'size(AZaaooTA, dim=1)', size(AZaaooTA, dim=1)
            print*, 'size(AZaaooTA, dim=2)', size(AZaaooTA, dim=2)

            print*, 'size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=1)', size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=1)
            print*, 'size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=2)', size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=2)

            print*, 'size(ACB(A1aaoo)%ATripA, dim=1)', size(ACB(A1aaoo)%ATripA, dim=1)
            print*, 'size(ACB(A1aaoo)%ATripA, dim=2)', size(ACB(A1aaoo)%ATripA, dim=2)
            print*, 'mnoże C= ^AT.B'
            print*, 'C to jest AZaaooTA o wymiarach', size(AZaaooTA, dim=1), size(AZaaooTA, dim=2)
            print*, 'A to jest ACB(A0aa)%MiniBlocks(2)%MiniAVTA o wymiarach', size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=1), size(ACB(A0aa)%MiniBlocks(2)%MiniAVTA, dim=2)
            print*, 'B to jest ACB(A1aaoo)%ATripA o wymiarach, ', size(ACB(A1aaoo)%ATripA, dim=1), size(ACB(A1aaoo)%ATripA, dim=2)
#endif

            
            call clock_start(timer)
            !            if (ACB(A0aa)%MiniBlocks(2)%NdimT.gt.0.and.ACB(A0oo)%NdimT.gt.0)then
            if (size_not_zero(AZaaooT, ACB(A0aa)%MiniBlocks(2)%MiniAVT, ACB(A1aaoo)%ATrip))then
                  call real_atb(AZaaooT, ACB(A0aa)%MiniBlocks(2)%MiniAVT, ACB(A1aaoo)%ATrip)
            end if
            !            if (ACB(A0aa)%MiniBlocks(2)%NdimTA.gt.0.and.ACB(A0oo)%NdimT.gt.0)then
            if (size_not_zero(AZaaooTA, ACB(A0aa)%MiniBlocks(2)%MiniAVTA, ACB(A1aaoo)%ATripA))then
                  call real_atb(AZaaooTA, ACB(A0aa)%MiniBlocks(2)%MiniAVTA, ACB(A1aaoo)%ATripA)
            end if
            call tmsg('TIME FOR AZ AAOOT ', timer, tdebug)

            !VAOO  intermediates-------------------------------------------------------------------------------------------------
            
            call clock_start(timer)            
            AZvaoo = transpose(ACB(A1vaoo)%ASing) !@@@ od razu konstruowac transponowana


!            if (ACB(A0oo)%NdimS.gt.0.and.NA.gt.0)then
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA
                  ! print*, size(AZmini, dim=1)
                  ! print*, size(AZmini, dim=2)
                  ! print*, size(AZvaoo, dim=1)
                  ! print*, size(AZvaoo, dim=2)
                  ! print*, size(ACB(A0va)%MiniBlocks(p)%MiniAVS, dim=1)
                  ! print*, size(ACB(A0va)%MiniBlocks(p)%MiniAVS, dim=1)
                  ! print*, l0, l1

                  if (size_not_zero(AZmini, AZvaoo(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZmini, AZvaoo(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if
                  
                  do n = 1, NA
                        do k = 1, ACB(A0oo)%NdimS
                              AZmini(k,n) = AZmini(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ASing(k, 1))
                              if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.   ACB(A0oo)%ASing(k, 1)> zero)then 
                                    write(*,'(2A12, 3F40.10)') 'Idenom', 'vaoomini',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ASing(k, 1)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ASing(k, 1)
                              end if
                              if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ASing(k, 1)).lt.1.d-6)then
                                    write(*, '(2A12, 3F40.10)')'Sdenom', 'vaoomini', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ASing(k, 1)),ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ASing(k, 1)
                              end if
                              
                        end do
                  end do
                  AZvaoo(:, l0:l1)  = AZmini
            end do
!            end if

            call tmsg('TIME FOR AZ VAOO ', timer, tdebug)
            call clock_start(timer)

            
            call clock_start(timer)            
            AZvaooT = transpose(ACB(A1vaoo)%ATrip) !@@@ od razu konstruowac transponowana
            AZvaooTA = transpose(ACB(A1vaoo)%ATripA)

#ifdef DEBUG
            print*, size(AZminiT, dim=1), 'size(AZminiT, dim=1)'
            print*, size(AZminiT, dim=2), 'size(AZminiT, dim=2)'
            print*, size(AZvaooT, dim=1), 'size(AZvaooT, dim=1)'
            print*, size(AZvaooT, dim=2), 'size(AZvaooT, dim=2)'
#endif

            ds1 = size(AZminiT, dim=1)
            ds2 = size(AZminiT, dim=2)
            ds3 = size(AZvaooT, dim=1)
            ds4 =  size(AZvaooT, dim=2)

            ! print*, size(ACB(A0va)%MiniBlocks(p)%MiniAVS, dim=1), 'z'
            ! print*, size(ACB(A0va)%MiniBlocks(p)%MiniAVS, dim=2), 'zz'

!            if (ds1.ne.0.and.ds2.ne.0.and.ds3.ne.0.and.ds4.ne.0)then
                  do p = 1, NV
                        l0 = 1 +  (p-1)*NA
                        l1 = NA + (p-1)*NA
                        if (size_not_zero(AZminiT, AZvaooT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiT, AZvaooT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if
                  
                        do n = 1, NA
                              do k = 1, ACB(A0oo)%NdimT
                                    AZminiT(k,n) = AZminiT(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATrip(k, 1))
                                    if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and. ACB(A0oo)%ATrip(k, 1)  > zero)then 
                                          write(*,'(2A12, 3F40.10)') 'Idenom', 'vaoominit',ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATrip(k, 1), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ATrip(k, 1)
                                    end if
                                    if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATrip(k, 1)).lt.1.d-6)then
                                          write(*, '(2A12, 3F40.10)')'Sdenom', 'vaoominit', ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATrip(k, 1), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ATrip(k, 1)
                                    end if
                              end do
                        end do
                        AZvaooT(:, l0:l1)  = AZminiT
                  end do
!            end if
!            if (ds1.ne.0.and.ds2.ne.0.and.ds3.ne.0.and.ds4.ne.0)then
                  do p = 1, NV
                        l0 = 1 +  (p-1)*NA
                        l1 = NA + (p-1)*NA
                        if (size_not_zero(AZminiTA, AZvaooTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiTA, AZvaooTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if
                  
                        do n = 1, NA
                              do k = 1, ACB(A0oo)%NdimT
                                    AZminiTA(k,n) = AZminiTA(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATripA(k, 1))
                                    if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and. ACB(A0oo)%ATripA(k, 1)  > zero)then 
                                          write(*,'(2A12, 3F40.10)') 'Idenom', 'vaoominit',ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATripA(k, 1), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ATripA(k, 1)
                                    end if
                                    if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATripA(k, 1)).lt.1.d-6)then
                                          write(*, '(2A12, 3F40.10)')'Sdenom', 'vaoominit', ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0oo)%ATripA(k, 1), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0oo)%ATripA(k, 1)
                                    end if
                              end do
                        end do
                        AZvaooTA(:, l0:l1)  = AZminiTA
                  end do
!            end if
            call tmsg('TIME FOR AZ VAOOT ', timer, tdebug)
            call clock_start(timer)

            !VAAA  intermediates-------------------------------------------------------------------------------------------------

            !            if(ACB(A0aa)%MiniBlocks(1)%NdimS.gt.0.and.ACB(A0va)%NdimS.gt.0)then
            if (size_not_zero(AZvaaa, ACB(A0aa)%MiniBlocks(1)%MiniAVS,ACB(A1vaaa)%ASing))then
                  call real_atbt(AZvaaa, ACB(A0aa)%MiniBlocks(1)%MiniAVS,ACB(A1vaaa)%ASing)
            end if

            call tmsg('TIME FOR AZ VAAA ', timer, tdebug)            
            call clock_start(timer)

!            if (ACB(A0aa)%MiniBlocks(1)%NdimS.gt.0.and.NA.gt.0)then
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA
                  if (size_not_zero(AZminiaa, AZvaaa(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiaa, AZvaaa(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if

                  do n = 1, NA
                        
                        do k = 1, ACB(A0aa)%MiniBlocks(1)%NdimS
                              AZminiaa(k,n) = AZminiaa(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k))
                               if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.   ACB(A0aa)%MiniBlocks(1)%EigS(k)> zero)then 
                                    write(*,'(2A12, 3F40.10)') 'Idenom', 'vaaamini',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigS(k)
                                    print*, (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigS(k)
                              end if
                              if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k)).lt.1.d-6)then
                                    write(*, '(2A12, 4F20.5)')'Sdenom', 'vaaamini', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigS(k), AZminiaa(k,n)
                              !       print*, 'sdenom',n,p, (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigS(k), AZminiaa(k,n)
                              end if
                        end do
                  end do
                  AZvaaa(:, l0:l1)  = AZminiaa
            end do
 !     end if

            call tmsg('TIME FOR AZ VAAA part2', timer, tdebug)            
            call clock_start(timer)

!            if (ACB(A0aa)%MiniBlocks(1)%NdimT.gt.0.and.ACB(A0va)%NdimS.gt.0)then
            !                  if (size(ACB(A0aa)%MiniBlocks(1)%MiniAVT,dim=1).gt.0)then
            if (size_not_zero(AZvaaaT, ACB(A0aa)%MiniBlocks(1)%MiniAVT,ACB(A1vaaa)%ATrip))then
                  call real_atbt(AZvaaaT, ACB(A0aa)%MiniBlocks(1)%MiniAVT,ACB(A1vaaa)%ATrip)
            end if
            if (size_not_zero(AZvaaaTA, ACB(A0aa)%MiniBlocks(1)%MiniAVTA,ACB(A1vaaa)%ATripA))then
                  call real_atbt(AZvaaaTA, ACB(A0aa)%MiniBlocks(1)%MiniAVTA,ACB(A1vaaa)%ATripA)
                  end if
!            end if

            call tmsg('TIME FOR AZ VAAAT ', timer, tdebug)            
            call clock_start(timer)

!            if (ACB(A0aa)%MiniBlocks(1)%NdimT.gt.0.and.NA.gt.0)then
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  if (size_not_zero(AZminiaaT, AZvaaaT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiaaT, AZvaaaT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if

                  do n = 1, NA                        
                        do k = 1, ACB(A0aa)%MiniBlocks(1)%NdimT
                              AZminiaaT(k,n) = AZminiaaT(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigT(k))
                                if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.  ACB(A0aa)%MiniBlocks(1)%EigT(k) > zero)then 
                      write(*,'(2A12, 3F40.10)') 'Idenom', 'vaaaminit',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigT(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigT(k)
                end if
                if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigT(k)).lt.1.d-6)then
                      write(*, '(2A12, 3F40.10)')'Sdenom', 'vaaaminit', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigT(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigT(k)
                end if
                        end do
                  end do
                  AZvaaaT(:, l0:l1)  = AZminiaaT
            end do
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  if (size_not_zero(AZminiaaTA, AZvaaaTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiaaTA, AZvaaaTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if

                  do n = 1, NA                        
                        do k = 1, ACB(A0aa)%MiniBlocks(1)%NdimT
                              AZminiaaTA(k,n) = AZminiaaTA(k,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigTA(k))
                                if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.  ACB(A0aa)%MiniBlocks(1)%EigTA(k) > zero)then 
                      write(*,'(2A12, 3F40.10)') 'Idenom', 'vaaaminit',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigTA(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigTA(k)
                end if
                if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigTA(k)).lt.1.d-6)then
                      write(*, '(2A12, 3F40.10)')'Sdenom', 'vaaaminit', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0aa)%MiniBlocks(1)%EigTA(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0aa)%MiniBlocks(1)%EigTA(k)
                end if
                        end do
                  end do
                  AZvaaaTA(:, l0:l1)  = AZminiaaTA
            end do
!      end if
            call tmsg('TIME FOR AZ VAAAT part2', timer, tdebug)            
            call clock_start(timer)


            !VAAO  intermediates ----------------------------------------------------------------------------------------------------------------------

!            if (ACB(A0va)%NdimS.gt.0.and.NA.gt.0)then
            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA

                  if (size_not_zero(AZminiao, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ASing(:, l0:l1)))then
                        call real_atbt(AZminiao, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ASing(:, l0:l1))
                  end if
                  AZvaao(l0:l1, :) = AZminiao
            end do
 !     end if
            call tmsg('TIME FOR AZ VAAO part1', timer, tdebug)            
            call clock_start(timer)

!            if (ACB(A0va)%NdimS.gt.0.and.NA.gt.0)then
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  if (size_not_zero(AZminiao2, AZvaao(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiao2, AZvaao(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if

                  do n = 1, NA
                        do r = 1, NI
                              do k = 1, NA
                                    km = (r-1)*NA + k
                                    
                                    AZminiao2(km,n) = AZminiao2(km,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k))
                                     if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.   ACB(A0ao)%MiniBlocks(r)%EigS(k)> zero)then 
                                          write(*,'(2A12, 3F40.10)') 'Idenom', 'vaaomini',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0ao)%MiniBlocks(r)%EigS(k)
                                    end if
                                    if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)).lt.1.d-6)then
                                          write(*, '(2A12, 3F40.10)')'Sdenom', 'vaaomini', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n), ACB(A0ao)%MiniBlocks(r)%EigS(k)
                                    end if
                              end do
                        end do
                  end do
                  AZvaao(:, l0:l1)  = AZminiao2
            end do
!      end if

            call tmsg('TIME FOR AZ VAA0 part2', timer, tdebug)            
            call clock_start(timer)


!            if (ACB(A0va)%NdimS.gt.0.and.NA.gt.0)then
            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA
                  if (size_not_zero(AZminiaoT, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ATrip(:, l0:l1)))then
                        call real_atbt(AZminiaoT, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ATrip(:, l0:l1))
                  end if
                  AZvaaoT(l0:l1, :) = AZminiaoT
            end do
            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA
                  if (size_not_zero(AZminiaoTA, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ATripA(:, l0:l1)))then
                        call real_atbt(AZminiaoTA, ACB(A0ao)%MiniBlocks(r)%MiniAVS,  ACB(A1vaao)%ATripA(:, l0:l1))
                  end if
                  AZvaaoTA(l0:l1, :) = AZminiaoTA
            end do
!      end if
            call tmsg('TIME FOR AZ VAAOT part1', timer, tdebug)            
            call clock_start(timer)

!            if (ACB(A0va)%NdimS.gt.0.and.NA.gt.0)then   
            do p = 1, NV
                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  if (size_not_zero(AZminiao2T, AZvaaoT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiao2T, AZvaaoT(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)                        
                  end if
                  if (size_not_zero(AZminiao2TA, AZvaaoTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS))then
                        call real_ab(AZminiao2TA, AZvaaoTA(:, l0:l1), ACB(A0va)%MiniBlocks(p)%MiniAVS)
                  end if
                    do n = 1, NA
                        do r = 1, NI
                              do k = 1, NA
                                    km = (r-1)*NA + k
                                    
                                    AZminiao2T(km,n) = AZminiao2T(km,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k))
                                    AZminiao2TA(km,n) = AZminiao2TA(km,n) / (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k))
                !                        if ( ACB(A0va)%MiniBlocks(p)%EigS(n)< zero .and.    ACB(A0ao)%MiniBlocks(r)%EigS(k)> zero)then 
                !       write(*,'(2A12, 3F40.10)') 'Idenom', 'vaaominit',(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n),  ACB(A0ao)%MiniBlocks(r)%EigS(k)
                ! end if
                ! if (abs(ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)).lt.1.d-6)then
                !       write(*, '(2A12, 3F40.10)')'Sdenom', 'vaaominit', (ACB(A0va)%MiniBlocks(p)%EigS(n) - ACB(A0ao)%MiniBlocks(r)%EigS(k)), ACB(A0va)%MiniBlocks(p)%EigS(n),  ACB(A0ao)%MiniBlocks(r)%EigS(k)
                ! end if
                              end do
                                            end do
                                      end do
                  AZvaaoT(:, l0:l1)  = AZminiao2T
                  AZvaaoTA(:, l0:l1)  = AZminiao2TA
            end do
  !    end if



            call tmsg('TIME FOR AZ VAA0T part2', timer, tdebug)            
            call clock_start(timer)



            ! VAAO inetermediates----------------------------------------------------------------------------------------------------------------

            if (size_not_zero(AZaaao, ACB(A0aa)%MiniBlocks(2)%MiniAVS, ACB(A1aaao)%ASing))then
                  call real_atb(AZaaao, ACB(A0aa)%MiniBlocks(2)%MiniAVS, ACB(A1aaao)%ASing)
            end if

            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA
                  if (size_not_zero(AZminiaaao, AZaaao(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS))then
                        call real_ab(AZminiaaao, AZaaao(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS)
                  end if

                  do k = 1, NA
                        do n = 1, ACB(A0aa)%MiniBlocks(2)%NdimS
                              AZminiaaao(n,k) = AZminiaaao(n,k) / (ACB(A0aa)%MiniBlocks(2)%EigS(n)-ACB(A0ao)%MiniBlocks(r)%EigS(k))
                        end do
                  end do

                  AZaaao(:, l0:l1) = AZminiaaao
            end do

            if (size_not_zero(AZaaaoT, ACB(A0aa)%MiniBlocks(2)%MiniAVT, ACB(A1aaao)%ATrip))then
                  call real_atb(AZaaaoT, ACB(A0aa)%MiniBlocks(2)%MiniAVT, ACB(A1aaao)%ATrip)
            end if

            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA
                  if (size_not_zero(AZminiaaaoT, AZaaaoT(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS))then
                        call real_ab(AZminiaaaoT, AZaaaoT(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS)
                  end if

                  do k = 1, NA
                        do n = 1, ACB(A0aa)%MiniBlocks(2)%NdimT
                              AZminiaaaoT(n,k) = AZminiaaaoT(n,k) / (ACB(A0aa)%MiniBlocks(2)%EigT(n)-ACB(A0ao)%MiniBlocks(r)%EigS(k))
                        end do
                  end do

                  AZaaaoT(:, l0:l1) = AZminiaaaoT
            end do

            if (size_not_zero(AZaaaoTA, ACB(A0aa)%MiniBlocks(2)%MiniAVTA, ACB(A1aaao)%ATripA))then
                  call real_atb(AZaaaoTA, ACB(A0aa)%MiniBlocks(2)%MiniAVTA, ACB(A1aaao)%ATripA)
            end if

            do r = 1, NI                  
                  l0 = 1 +  (r-1)*NA
                  l1 = NA + (r-1)*NA

                  if (size_not_zero(AZminiaaaoTA, AZaaaoTA(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS))then
                        call real_ab(AZminiaaaoTA, AZaaaoTA(:, l0:l1), ACB(A0ao)%MiniBlocks(r)%MiniAVS)
                  end if

                  do k = 1, NA
                        do n = 1, ACB(A0aa)%MiniBlocks(2)%NdimT
                              AZminiaaaoTA(n,k) = AZminiaaaoTA(n,k) / (ACB(A0aa)%MiniBlocks(2)%EigTA(n)-ACB(A0ao)%MiniBlocks(r)%EigS(k))
                        end do
                  end do

                  AZaaaoTA(:, l0:l1) = AZminiaaaoTA
            end do
!      end if
      !stop stop stop plusz
      
            timer10 = zero
            timer11 = zero
            timer12 = zero
            timer13 = zero
            timer13x = zero
            timer13y = zero
            timer13z = zero
            timer14 = zero
            timer14a = zero
            timer15 = zero
            timer16 = zero
            timer17 = zero
            timer18 = zero
            timer18a = zero
            timer33 = zero
            timer34 = zero
            timer35 = zero
            timer36 = zero
            timer37 = zero
            timer38 = zero
            timer39 = zero
            timer21 = zero





            NBatchI = NI / BatchDim
            NBatchA = NA / BatchDim
            NBatchV = NV / BatchDim
            NBatchIA = NIA / BatchDim


            if (modulo(NI, BatchDim)>0)NBatchI = NBatchI + 1
            if (modulo(NA, BatchDim)>0)NBatchA = NBatchA + 1
            if (modulo(NV, BatchDim)>0)NBatchV = NBatchV + 1
            if (modulo(NIA, BatchDim)>0)NBatchIA = NBatchIA + 1



            if (TabChol(2) == 1) then
                  allocate(RIA(Nchol, NA, min(NI, BatchDim)))
                  call CalcMem3(RIA, 'RIA')
            end if
            if (TabChol(2) == 1) then
                  allocate(RAI(Nchol, NI, min(NA, BatchDim)))
                  call CalcMem3(RAI, 'RAI')
            end if
            if (TabChol(3) == 1)then
                  allocate(RAA(Nchol, NA, min(NA, BatchDim)))
                  call CalcMem3(RAA, 'RAA')
            end if
            if (TabChol(3) == 1)then
                  allocate(ROO(Nchol, NIA, NIA))
                  call CalcMem3(ROO, 'ROO')
            end if
            
            if (TabChol(4) == 1)then
                  allocate(RIV(Nchol, NV, min(NI, BatchDim)))
                  call CalcMem3(RIV, 'RIV')
            end if

            if (TabChol(4) == 1)then
                  allocate(RVI(Nchol, NI, min(NV, BatchDim)))
                  call CalcMem3(RVI, 'RVI')
            end if
            
            if (TabChol(5) == 1)then
                  allocate(RVA(Nchol, NA, min(NV, BatchDim)))
                  call CalcMem3(RVA, 'RVA')
            end if

            if (TabChol(5) == 1)then
                  allocate(RVO(Nchol, NIA, min(NV, BatchDim)))
                  call CalcMem3(RVO, 'RVO')
            end if


            call clock_start(timer)
            if (TabChol(2) == 1) call thc_gammcor_Rkab_2(RIA, Xga(:,NI+1:NIA), Xga(:,1:NI), Zgk, NA, NI, NChol, NTHC)
            call tmsg('TIME FOR RIA', timer, tdebug)            
            call clock_start(timer)


            call clock_start(timer)
            if (TabChol(2) == 1) call thc_gammcor_Rkab_2(RAI, Xga(:,1:NI), Xga(:,NI+1:NIA), Zgk, NI, NA, NChol, NTHC)
            call tmsg('TIME FOR RAI', timer, tdebug)            
            call clock_start(timer)

            
            if (TabChol(3) == 1) call thc_gammcor_Rkab_2(RAA, Xga(:,NI+1:NIA), Xga(:,NI+1:NIA), Zgk, NA, NA, NChol, NTHC)
            call tmsg('TIME FOR RAA ', timer, tdebug)            
            call clock_start(timer)

            if (TabChol(3) == 1) call thc_gammcor_Rkab_2(ROO, Xga(:,1:NIA), Xga(:,1:NIA), Zgk, NIA, NIA, NChol, NTHC)

            call tmsg('TIME FOR ROO ', timer, tdebug)            
            call clock_start(timer)

            call print_info('NBatchI', NBatchI)
            call print_info('NBatchA', NBatchA)
            call print_info('NBatchV', NBatchV)

            batchloop: do batch = 1, max(NBatchI, NBatchA, NBatchV)
 
                  b0IV = 1 + (batch-1) * BatchDim
                  b1IV = min(b0IV+BatchDim-1, NI)

                  b0VI = NIA + 1 + (batch-1) * BatchDim
                  b1VI = min(b0VI+BatchDim-1, NBasis)

                  ! b0AV = NI + 1 + (batch-1) * BatchDim
                  ! b1AA = min(b0AV+BatchDim-1, NIA)

                  b0VA = NIA + 1 + (batch-1) * BatchDim
                  b1VA = min(b0VA+BatchDim-1, NBasis)
                  
                  b0VV = NIA + 1 + (batch-1) * BatchDim
                  b1VV = min(b0VV+BatchDim-1, NBasis)


                  b0VO = NIA + 1 + (batch-1) * BatchDim
                  b1VO = min(b0VO+BatchDim-1, NBasis)

                  call ximsg('b0IV', b0IV, mdebug)
                  call ximsg('b1IV', b1IV, mdebug)

                  if (b0IV <=NI)then
                        if (TabChol(4)==1)then
                              call clock_start(timer)
                              do b = b0IV, b1IV
                                    call thc_gammcor_Rkab_Batch_a_Fixed_b(RIV(:, 1:NV, b-b0IV+1), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, b), Zgk, NV, NChol, NTHC)
                              end do
                              call tmsg('TIME FOR RIV ', timer, tdebug)            
                        end if
                  end if

                  call ximsg('b0VI', b0VI, mdebug)
                  call ximsg('b1VI', b1VI, mdebug)

                  
                              
                  
                  if (b0VI <= NBasis)then
                        if (TabChol(4)==1)then
                              call clock_start(timer)
                              do b = b0VI, b1VI
                                    call thc_gammcor_Rkab_Batch_a_Fixed_b(RVI(:, 1:NI, b-b0VI+1), XXga(:, 1:NI), Xga(:, 1:NI), Xga(:, b), Zgk, NI, NChol, NTHC)                                    
                              end do
                              call tmsg('TIME FOR RVI ', timer, tdebug)            
                        end if
                  end if


                  if (b0VA <=NBasis)then
                        if (TabChol(5)==1)then
                              call clock_start(timer)
                              do b = b0VA, b1VA
                                    call thc_gammcor_Rkab_Batch_a_Fixed_b(RVA(:, 1:NA, b-b0VA+1), XXga(:, 1:NA), Xga(:, NI+1:NIA), Xga(:, b), Zgk, NA, NChol, NTHC)	
                              end do
                              call tmsg('TIME FOR RVA ', timer, tdebug)            
                        end if
                  end if

                  if (b0VO <=NBasis)then
                        if (TabChol(5)==1)then
                              call clock_start(timer)
                              do b = b0VO, b1VO
                                    call thc_gammcor_Rkab_Batch_a_Fixed_b(RVO(:, 1:NIA, b-b0VO+1), XXga(:, 1:NIA), Xga(:, 1:NIA), Xga(:, b), Zgk, NIA, NChol, NTHC)
                              end do
                              call tmsg('TIME FOR RVO ', timer, tdebug)            
                        end if
                  end if

                  

                  caseloop: do ii = 1, 6
                        mainloopscase: if (loops(ii)%run_main)then
                              select case(ii)
                                    
                              case(loopIA) ! loop2

                                    tloopIA: do a = 1, NI                                          
                                          if (loops(ii)%run_inner(2)) then! (IA|IA)                                                
                                                do b = 1, a

                                                      call clock_start(timer)
                                                      call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                      timer10 = timer10 +  clock_readwall(timer)

                                                      call clock_start(timer)
                                                      
                                                      rs = posS(a, b)
                                                      rsT = posT(a, b)
                                                      if (rs > 0)then
                                                            call aaoo_contr(EAAoo, AZaaoo, ACB(A0aa), ACB(A0oo),  rs, a, b, NI+1, NIA, NI+1, NIA, AuxData%Occ, V_axby, WorkAA)
                                                      end if
                                                      if (rsT >0)then
                                                            call aaoo_contrT(EAAooT, AZaaooT, ACB(A0aa), ACB(A0oo),  rsT, a, b, NI+1, NIA, NI+1, NIA, AuxData%Occ, V_axby, WorkAA)
                                                            call aaoo_contrTA(EAAooTA, AZaaooTA, ACB(A0aa), ACB(A0oo),  rsT, a, b, NI+1, NIA, NI+1, NIA, AuxData%Occ, V_axby, WorkAA)
                                                      end if
                                                      timer11 = timer11+  clock_readwall(timer)
                                                      
                                                end do
                                          end if

                                          do b = NI+1, NIA
                                                call clock_start(timer)
                                                call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)

                                                timer10 = timer10 +  clock_readwall(timer)

                                                call clock_start(timer)
                                                      
                                                rs = posS(a, b)
                                                call aaao_contr(EAAao, AZaaao, ACB(A0aa), ACB(A0ao),  rs, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, WorkAA)
                                                rsT = posT(a, b)
                                                call aaao_contrT(EAAaoT, AZaaaoT, ACB(A0aa), ACB(A0ao),  rsT, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, WorkAA)
                                                call aaao_contrTA(EAAaoTA, AZaaaoTA, ACB(A0aa), ACB(A0ao),  rsT, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, WorkAA)
                                                
                                                
                                          end do

                                    end do tloopIA
                                    

                              case(loopIV)! loop 4
                                    tloopIV: do a = 1, NI
                                          call clock_start(timer)
                                          if (a<b0IV.or.a>b1IV)then
                                                call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NV), XXga(:, 1:NV), Xga(:, NIA+1:NBasis), Xga(:, a), Zgk, NV, NChol, NTHC)
                                          end if
                                          timer21 = timer21 +  clock_readwall(timer)


                                          if (loops(ii)%run_inner(2))then ! (IV|IA)
                                                do b = 1, NI
                                                      call clock_start(timer)
                                                      if (a<b0IV.or.a>b1IV)then
                                                            call real_aTb_x(V_axby, NV, Rkax, NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                            call real_aTb_x(V_bxay, NA, RIA(:,:,a), NChol, Rkax, NChol, NA, Nv, NChol, AcAlpha)
                                                      else
                                                            call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha) ! (aV|bA)
                                                            call real_aTb_x(V_bxay, NA, RIA(:,:,a), NChol, RIV(:,:,b), NChol, NA, Nv, NChol, AcAlpha) ! (aA|bV)
                                                      end if                                                      
                                                      timer12 = timer12 +  clock_readwall(timer)
                                                      call clock_start(timer)
                                                      
                                                      rs = posS(a, b)
                                                      if (rs > 0)then
                                                            call vaoo_contr(EVAoo, AZvaoo, ACB(A0va), ACB(A0oo),  rs, a, b, 1, NV, 1, NA, AuxData, V_axby, V_bxay, workAA)
                                                      end if
                                                      timer18a = timer18a+  clock_readwall(timer)
                                                      call clock_start(timer)
                                                      rsT = posT(a, b)
                                                      if (rsT > 0)then
                                                            call vaoo_contrT(EVAooT, AZvaooT, ACB(A0va), ACB(A0oo),  rsT, a, b, 1, NV, 1, NA, AuxData, V_axby, V_bxay, workAA)
                                                            call vaoo_contrTA(EVAooTA, AZvaooTA, ACB(A0va), ACB(A0oo),  rsT, a, b, 1, NV, 1, NA, AuxData, V_axby, V_bxay, workAA)
                                                      end if

                                                      timer18 = timer18+  clock_readwall(timer)
                                                      
                                                end do
                                          end if
                                          
                                          if (loops(ii)%run_inner(4))then  ! (IV|IV)
                                                do b = b0IV, min(b1IV, a)
                                                      
                                                      call clock_start(timer)
                                                      if (a<b0IV.or.a>b1IV)then
                                                            call real_aTb_x(V_axby, NV, Rkax, NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                      else
                                                            call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIV(:,:,b), NChol, NV, NV, NChol, AcAlpha)
                                                            
                                                      end if
                                                      timer33 = timer33 +  clock_readwall(timer)	

                                                      call clock_start(timer)

                                                      pq = posS(a, b)
                                                      if (pq > 0)then
                                                            pqT = posT(a, b)
                                                            call vvoo_contr(EVVoo, EVVooT, EVVooTA, ACB(A0vv), ACB(A0oo),  pq, pqT, a, b, NIA+1, NBasis, NIA+1, NBasis,  V_axby)


                                                      end if
                                                      

                                                      timer34 = timer34 +  clock_readwall(timer)	                                                    
                                                end do
                                          end if
                                    end do tloopIV

                              case(loopVA) ! loop 5
                                    call clock_start(timer)
                                    do a = NIA+1, NBasis
                                          do b = b0VO, min(b1VO, a)
                                                call real_aTb_x(V_axby, NIA, RVO(:,:,a-NIA), NChol, RVO(:,:,b-NIA), NChol, NIA, NIA, NChol, AcAlpha)
                                          end do
                                    end do
                                    timer13x = timer13x + clock_readwall(timer)


                                    call clock_start(timer)
                                    !$omp parallel do default(shared) private(a, b, V_axby)
                                    do a = NIA+1, NBasis
                                          do b = b0VO, min(b1VO, a)
                                                call real_aTb_x(V_axby, NIA, RVO(:,:,a-NIA), NChol, RVO(:,:,b-NIA), NChol, NIA, NIA, NChol, AcAlpha)
                                          end do
                                    end do
                                    !$omp end parallel do 
                                    timer13y = timer13y + clock_readwall(timer)


                                    call clock_start(timer)
                                    !$omp parallel do default(shared) private(a, b, V_axby)&
                                    !$omp  collapse(2)                                    
                                    do a = NIA+1, NBasis
                                          do b = b0VO, b1VO
                                                call real_aTb_x(V_axby, NIA, RVO(:,:,a-NIA), NChol, RVO(:,:,b-NIA), NChol, NIA, NIA, NChol, AcAlpha)
                                          end do
                                    end do
                                    !$omp end parallel do 
                                    timer13z = timer13z + clock_readwall(timer)




                                                
                                    tloopVO: do a = NIA+1, NBasis
                                          if (a<b0VO.or.a>b1VO)then
                                                call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NIA), XXga(:, 1:NIA), Xga(:, 1:NIA), Xga(:, a), Zgk, NIA, NChol, NTHC)
                                          end if
                                          
                                          do b = NI+1, NIA
                                                call clock_start(timer)
                                                if (a<b0VO.or.a>b1VO)then
                                                      call real_aTb_x(V_axby, NIA, Rkax, NChol, RAA(:,:,b-NI), NChol, NIA, NA, NChol, AcAlpha)

                                                else
                                                      call real_aTb_x(V_axby, NIA, RVO(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NIA, NA, NChol, AcAlpha)
                                                      call real_aTb_x(V_bxay, NIA, RVO(:,:,a-NIA), NChol, RAI(:,:,b-NI), NChol, NIA, NI, NChol, AcAlpha)
                                                end if
                                                timer35 = timer35 + clock_readwall(timer)

                                                call clock_start(timer)
                                                pq = posS(a, b)
                                                call vaao2_contr(EVAao2, AZvaao, ACB(A0va), ACB(A0ao), pq, a, b, 1, NIA, NI+1, NIA, AuxData, V_axby, workAO)
                                                timer36 = timer36 + clock_readwall(timer)

                                                call clock_start(timer)
                                                pqT = posT(a, b)
                                                call vaao2_contrT(EVAao2T, AZvaaoT, ACB(A0va), ACB(A0ao), pqT, a, b, 1, NIA, NI+1, NIA, AuxData, V_axby, workAO)
                                                call vaao2_contrTA(EVAao2TA, AZvaaoTA, ACB(A0va), ACB(A0ao), pqT, a, b, 1, NIA, NI+1, NIA, AuxData, V_axby, workAO)
                                                timer37 = timer37 + clock_readwall(timer)

                                                call clock_start(timer)
                                                call vaao_contr(EVAao, AZvaao, ACB(A0va), ACB(A0ao), pq, a, b, 1, NIA, 1, NI, AuxData, V_bxay, workAO)
                                                timer38 = timer38 + clock_readwall(timer)

                                                call clock_start(timer)
                                                call vaao_contrT(EVAaoT, AZvaaoT, ACB(A0va), ACB(A0ao), pqT, a, b, 1, NIA, 1, NI, AuxData, V_bxay, workAO)
                                                call vaao_contrTA(EVAaoTA, AZvaaoTA, ACB(A0va), ACB(A0ao), pqT, a, b, 1, NIA, 1, NI, AuxData, V_bxay, workAO)
                                                timer39 = timer39 + clock_readwall(timer)
                                                
                                          end do

                                          
                                          do b = b0VO, min(b1VO, a) !(VO|VO)
                                                call clock_start(timer)
                                                if (a<b0VO.or.a>b1VO)then
                                                      call real_aTb_x(V_axby, NIA, Rkax, NChol, RVO(:,:,b-NIA), NChol, NIA, NIA, NChol, AcAlpha)

                                                else
                                                      call real_aTb_x(V_axby, NIA, RVO(:,:,a-NIA), NChol, RVO(:,:,b-NIA), NChol, NIA, NIA, NChol, AcAlpha)
                                                end if

                                                timer13 = timer13 + clock_readwall(timer)

                                                call clock_start(timer)
                                                pq = posS(a, b)
                                                if (pq > 0)then
                                                      call vvao_contr(EVVao, EVVaoT, ACB(A0vv), ACB(A0ao), pq, pqT, a, b, 1, NIA, 1, NIA, AuxData, V_axby, workAO, workAOT)
                                                end if
                                                timer14 = timer14 + clock_readwall(timer)

                                                call clock_start(timer)
                                                pqT = posT(a, b)
                                                if (pqT > 0)then
                                                      call vvao_contrT(EVVaoT, ACB(A0vv), ACB(A0ao), pq, pqT, a, b, 1, NIA, 1, NIA, AuxData, V_axby)
                                                      call vvao_contrTA(EVVaoTA, ACB(A0vv), ACB(A0ao), pq, pqT, a, b, 1, NIA, 1, NIA, AuxData, V_axby)
                                                end if
                                                timer14a = timer14a + clock_readwall(timer)

                                                
                                                call clock_start(timer)
                                                pq = posS(a, b)
                                                pqT = posT(a, b)
                                                if (pq > 0)then
                                                      call vvaa_contr(EVVaa, EVVaaT, EVVaaTA, ACB(A0vv), ACB(A0aa), pq, pqT, a, b, 1, NIA, 1, NIA, AuxData%Occ, V_axby, workAA)
                                                      !call vvaa_contrTA(EVVaaTA, ACB(A0vv), ACB(A0aa), pqT, a, b, 1, NIA, 1, NIA, AuxData%Occ, V_axby, workAA)
                                                end if
                                                timer15 = timer15 + clock_readwall(timer)
                                          end do

                                    end do tloopVO

                                    tloopVA: do a = NIA+1, NBasis
                                          if (a<b0VA.or.a>b1VA)then
                                                call clock_start(timer)
                                                call thc_gammcor_Rkab_Batch_a_Fixed_b(Rkax(:, 1:NA), XXga(:, 1:NA), Xga(:, NI+1:NIA), Xga(:, a), Zgk, NA, NChol, NTHC)
                                                timer37 = timer37 + clock_readwall(timer)
                                          end if
                                          do b = NI+1, NIA !(VA|AA)
                                                call clock_start(timer)
                                                if (a<b0VA.or.a>b1VA)then
                                                      call real_aTb_x(V_axby, NA, Rkax, NChol, RAA(:,:,b-NA), NChol, NA, NA, NChol, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                timer16 = timer16 + clock_readwall(timer)

                                                call clock_start(timer)
                                                pq = posS(a, b)
                                                call vaaa_contr(EVAaa, AZvaaa, ACB(A0va), ACB(A0aa), pq, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, workAA)

                                                pqT = posT(a, b)
                                                call vaaa_contrT(EVAaaT, AZvaaaT, ACB(A0va), ACB(A0aa), pqT, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, workAA)
                                                call vaaa_contrTA(EVAaaTA, AZvaaaTA, ACB(A0va), ACB(A0aa), pqT, a, b, NI+1, NIA, NI+1, NIA, AuxData, V_axby, workAA)

                                                timer17 = timer17 + clock_readwall(timer)
                                          end do
                                    end do tloopVA

                              end select
                        end if mainloopscase
                  end do caseloop
            end do batchloop



#ifdef DEBUG
      call xmsg('timer10', timer10, mdebugg)
      call xmsg('timer11', timer11, mdebugg)
      call xmsg('timer12', timer12, mdebugg)
      call xmsg('timer18', timer18, mdebugg)
      call xmsg('timer18a', timer18, mdebugg)

      call xmsg('timer33', timer33, mdebugg)
      call xmsg('timer34', timer34, mdebugg)
      call xmsg('timer35', timer35, mdebugg)
      call xmsg('timer36', timer36, mdebugg)
      call xmsg('timer37', timer37, mdebugg)
      call xmsg('timer38', timer38, mdebugg)
      call xmsg('timer39', timer39, mdebugg)

      
      call xmsg('timer13', timer13, mdebugg)
      call xmsg('timer13x', timer13, mdebugg)
      call xmsg('timer13y', timer13, mdebugg)
      call xmsg('timer13z', timer13, mdebugg)
      
      call xmsg('timer14', timer14, mdebugg)
      call xmsg('timer14a', timer14a, mdebugg)
      call xmsg('timer15', timer15, mdebugg)
      call xmsg('timer16', timer16, mdebugg)
      call xmsg('timer17', timer17, mdebugg)
#endif

            
      ! print*, ''
      ! print*, EVVoo
      ! call dmsg("Energy contribution Evvoo", EVVoo)

      ! print*, ''
      ! print*, EVVao
      ! call dmsg("Energy contribution Evvao", EVVao)

      ! print*, ''
      ! print*, EVVaa
      ! call dmsg("Energy contribution Evvaa", EVVaa)

      ! print*, ''
      ! print*, EAAoo
      ! call dmsg("Energy contribution Eaaoo", EAAoo)

      ! print*, ''
      ! print*, EVAoo
      ! call dmsg("Energy contribution Evaoo", EVAoo)

      ! print*, ''
      ! print*, EVAaa
      ! call dmsg("Energy contribution Evaaa", EVAaa)


      ! print*, ''
      ! print*, EVAao
      ! call dmsg("Energy contribution Evaao", EVAao)

      ! print*, ''
      ! print*, EVAao2
      ! call dmsg("Energy contribution Evaao2", EVAao2)

      ! print*, ''
      ! print*, EAAao
      ! call dmsg("Energy contribution Eaaao", EAAao)



      ! print*, ''
      ! print*, EVVooT
      ! call dmsg("Energy contribution EvvooT", EVVooT)

      ! print*, ''
      ! print*, EVVaoT
      ! call dmsg("Energy contribution EvvaoT", EVVaoT)

      ! print*, ''
      ! print*, EVVaaT
      ! call dmsg("Energy contribution EvvaaT", EVVaaT)

      ! print*, ''
      ! print*, EAAooT
      ! call dmsg("Energy contribution EaaooT", EAAooT)

      ! print*, ''
      ! print*, EVAooT
      ! call dmsg("Energy contribution EvaooT", EVAooT)

      ! print*, ''
      ! print*, EVAaaT
      ! call dmsg("Energy contribution EvaaaT", EVAaaT)


      ! print*, ''
      ! print*, EVAaoT
      ! call dmsg("Energy contribution EvaaoT", EVAaoT)

      ! print*, ''
      ! print*, EVAao2T
      ! call dmsg("Energy contribution Evaao2T", EVAao2T)

      ! print*, ''
      ! print*, EAAaoT
      ! call dmsg("Energy contribution EaaaoT", EAAaoT)



      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(aaoo) = ,", EAAoo+EAAooT, ',', ttoeV * (EAAoo+EAAooT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(aaao) = ,", EAAao+EAAaoT, ',', ttoeV * (EAAao+EAAaoT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vaoo) = ,", EVAoo+EVAooT, ',', ttoeV * (EVAoo+EVAooT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vaao2) = ,", EVAao2+EVAao2T, ',', ttoeV * (EVAao2+EVAao2T)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vvoo) = ,", EVVoo+EVVooT, ',', ttoeV * (EVVoo+EVVooT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vaao) = ,", EVAao+EVAaoT, ',', ttoeV * (EVAao+EVAaoT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vaaa) = ,", EVAaa+EVAaaT, ',', ttoeV * (EVAaa+EVAaaT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vvao) = ,", EVVao+EVVaoT, ',', ttoeV * (EVVao+EVVaoT)
      ! write(*, '(A18, F20.15, A3, F20.15)')"E_contr(vvaa) = ,", EVVaa+EVVaaT, ',', ttoeV * (EVVaa+EVVaaT)

      E_Corr = EAAoo+EAAooT + EAAao+EAAaoT+EVAoo+EVAooT+EVAao+EVAaoT+EVAao2+EVAao2T + &
            EVVoo+EVVooT+EVAaa+EVAaaT+EVVao+EVVaoT+ EVVaa+EVVaaT + two * (EAAooTA + EAAaoTA + EVVooTA + EVVaoTA + EVVaaTA + EVAooTA + EVAaaTA + EVAaoTA + EVAao2TA)

!      write(*, '(A18, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3, F20.15, A3)')&
!            'E_contr_all', ',', EAAoo+EAAooT, ',',EAAao+EAAaoT, ',',EVAoo+EVAooT, ',',EVAao2+EVAao2T, ',',EVVoo+EVVooT, ',',EVAao+EVAaoT, ',',EVAaa+EVAaaT, ',',EVVao+EVVaoT, ',',EVVaa+EVVaaT
      !write(*, '(A15, 2F25.15)')"E_ppAC0 = ", E_corr, ttoeV * E_corr

      call print_section('ppAC0 energy components')
      write(*,'(2X,A12,4(1X,A14))') 'Component', 'Total [Eh]', 's_ABAB [Eh]', &
            't_ABAB [Eh]', '2*t_AAAA [Eh]'
      write(*,'(2X,A)') repeat('-', 72)
      call print_energy_component('E_vvoo', EVVoo+EVVooT+two*EVVooTA, EVVoo, EVVooT, two*EVVooTA)
      call print_energy_component('E_vvao', EVVao+EVVaoT+two*EVVaoTA, EVVao, EVVaoT, two*EVVaoTA)
      call print_energy_component('E_vvaa', EVVaa+EVVaaT+two*EVVaaTA, EVVaa, EVVaaT, two*EVVaaTA)
      call print_energy_component('E_vaoo', EVAoo+EVAooT+two*EVAooTA, EVAoo, EVAooT, two*EVAooTA)
      call print_energy_component('E_vaao', EVAao+EVAaoT+two*EVAaoTA, EVAao, EVAaoT, two*EVAaoTA)
      call print_energy_component('E_vaaa', EVAaa+EVAaaT+two*EVAaaTA, EVAaa, EVAaaT, two*EVAaaTA)
      call print_energy_component('E_aaoo', EAAoo+EAAooT+two*EAAooTA, EAAoo, EAAooT, two*EAAooTA)
      call print_energy_component('E_aaao', EAAao+EAAaoT+two*EAAaoTA, EAAao, EAAaoT, two*EAAaoTA)
      call print_energy_component('E_vaao2', EVAao2+EVAao2T+two*EVAao2TA, EVAao2, EVAao2T, two*EVAao2TA)

      call print_section('Contributions to ppAC0 from (Mu)(Nu) pairs of blocks')
      write(*,'(2X,A10,2X,A18)') 'Block pair', 'Contribution [Eh]'
      write(*,'(2X,A)') repeat('-', 30)
      call print_block_contribution('(22)(11)', EAAoo + EAAooT + two * EAAooTA)
      call print_block_contribution('(22)(21)', EAAao + EAAaoT + two * EAAaoTA)
      call print_block_contribution('(32)(11)', EVAoo + EVAooT + two * EVAooTA)
      call print_block_contribution('(32)(21)*', EVAao2 + EVAao2T + two * EVAao2TA)
      call print_block_contribution('(33)(11)', EVVoo + EVVooT + two * EVVooTA)
      call print_block_contribution('(32)(21)', EVAao  + EVAaoT  + two * EVAaoTA)
      call print_block_contribution('(32)(22)', EVAaa + EVAaaT + two * EVAaaTA)
      call print_block_contribution('(33)(21)', EVVao + EVVaoT + two * EVVaoTA)
      call print_block_contribution('(33)(22)', EVVaa + EVVaaT + two * EVVaaTA)
      print*, ''
      write(*,'(2X,A)') '*  term used as the replacement in ff'

      call print_section('Final energies')
      call print_energy('CASSCF energy (one-electron)', AuxData%ECAS_oneelectr, &
            ttoeV * AuxData%ECAS_oneelectr)
      call print_energy('ECASSCF_THC', AuxData%ECAS_THC, ttoeV * AuxData%ECAS_THC)
      call print_energy('E_ppAC0', E_corr, ttoeV * E_corr)
      call print_energy('E_total', AuxData%ECAS_THC + E_corr, &
            ttoeV * (AuxData%ECAS_THC + E_corr))
      print*, ''

    end associate
end subroutine THC_energy_loop

function size_not_zero(c, a, b)
      logical :: size_not_zero
      double precision, dimension(:,:), intent(in) :: c, a, b
      integer :: k, l, m, n, u, v

      size_not_zero = .true.
      
      k = size(c,dim=1)
      l = size(c,dim=2)

      m = size(a,dim=1)
      n = size(a,dim=2)

      u = size(b,dim=1)
      v = size(b,dim=2)

      if (k==0.or.l==0.or.m==0.or.n==0.or.u==0.or.v==0)then
            size_not_zero = .false.
      end if
      
      
end function size_not_zero


    subroutine vvoo_contr(EVVoo, EVVooT, EVVooTA, VV, OO,  pq, pqT, p, q, x0, x1, y0, y1, V_axby)

          double precision, intent(inout) :: EVVoo, EVVooT, EVVooTA
          type(TAC0Block), intent(in) :: VV, OO
          integer, intent(in) :: p, q, pq, pqT
          integer, intent(in) :: x0, x1, y0, y1
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          integer :: r, s, rs

          double precision :: val, valT, NumH, NumF

          NumH = One
          NumF = One

          if (p == q) NumH = sqrt(frac12)

          do rs = 1, VV%NdimS
                
                r = VV%IndNS(1, rs)
                s = VV%IndNS(2, rs)
                
                if (r==s) then
                      NumF = sqrt(frac12)
                else
                      NumF = One
                end if

                val =  NumF * NumH * (V_axby(r,s)+V_axby(s,r))

                !write(*,'(A10, 4I5, 9F15.8)') 'sniez', r, s, p, q, EVVoo , NumF , NumH , V_axby(r,s),V_axby(s,r), val**2 , V_axby(r,s), VV%ASing(rs, 1) , OO%ASing(pq, 1)
                EVVoo = EVVoo - val**2 / (VV%ASing(rs, 1) - OO%ASing(pq, 1))
                !write(*,'(A10, 4I5, 5F15.8)') 'sniez', r, s, p, q, EVVoo , val**2 , V_axby(r,s), VV%ASing(rs, 1) , OO%ASing(pq, 1)

          end do


          if (pqT>0)then
                do rs = 1, VV%NdimT

                      r = VV%IndNT(1, rs)
                      s = VV%IndNT(2, rs)

                      valT =  (V_axby(r,s)-V_axby(s,r))

                      EVVooT = EVVooT -  valT**2 / (VV%ATrip(rs, 1) - OO%ATrip(pqT, 1))
                      EVVooTA = EVVooTA -  valT**2 / (VV%ATripA(rs, 1) - OO%ATripA(pqT, 1))

                end do
          end if
    
    end subroutine vvoo_contr

    subroutine vvao_contr(EVVao, EVVaoT, VV, AO,  pq, pqT, p, q, x0, x1, y0, y1, AuxData, V_axby, workAO, workAOT)

          double precision, intent(inout) :: EVVao, EVVaoT
          type(TAC0Block), intent(in) :: VV, AO
          integer, intent(in) :: p, q, pq, pqT
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          double precision, dimension(:), intent(inout) :: workAO, workAOT
          integer :: t,u, km

          double precision :: val, valT, NumH, NumF



          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA)



            NumH = One
            NumF = One
            if (p == q) NumH = sqrt(frac12)
            !$omp parallel do default(shared) private(t, km, u, val) &
            !$omp reduction(+: EVVao) collapse(2)
            do t = 1, NI
                  do km = 1, NA
                        val = zero
                        do u =NI+1, NIA
                              val = val + (One-Occ(t)-Occ(u))*NumH * (V_axby(t,u)+V_axby(u,t))*AO%MiniBlocks(t)%MiniAVS(u-NI, km)
                        end do
                        EVVao = EVVao - val**2 / (VV%ASing(pq, 1) - AO%MiniBlocks(t)%EigS(km))
                  end do
            end do
            !$omp end parallel do      
            

          end associate
    end subroutine vvao_contr

    subroutine vvao_contrT(EVVaoT, VV, AO,  pq, pqT, p, q, x0, x1, y0, y1, AuxData, V_axby)

          double precision, intent(inout) :: EVVaoT
          type(TAC0Block), intent(in) :: VV, AO
          integer, intent(in) :: p, q, pq, pqT
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          integer :: t,u, km

          double precision :: val

          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA)

            !$omp parallel do default(shared) private(t, km, u, val) &
            !$omp reduction(+: EVVaoT) collapse(2)
            do t = 1, NI
                  do km = 1, NA
                        val = zero
                        do u =NI+1, NIA
                              val = val + (One-Occ(t)-Occ(u))*(V_axby(t,u)-V_axby(u,t))*AO%MiniBlocks(t)%MiniAVS(u-NI, km)
                        end do
                        EVVaoT = EVVaoT - val**2 / (VV%ATrip(pqT, 1) - AO%MiniBlocks(t)%EigS(km))
                  end do
            end do
            !$omp end parallel do      
          end associate
    end subroutine vvao_contrT


        subroutine vvao_contrTA(EVVaoTA, VV, AO,  pq, pqT, p, q, x0, x1, y0, y1, AuxData, V_axby)

          double precision, intent(inout) :: EVVaoTA
          type(TAC0Block), intent(in) :: VV, AO
          integer, intent(in) :: p, q, pq, pqT
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          integer :: t,u, km

          double precision :: val

          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA)

            !$omp parallel do default(shared) private(t, km, u, val) &
            !$omp reduction(+: EVVaoTA) collapse(2)
            do t = 1, NI
                  do km = 1, NA
                        val = zero
                        do u =NI+1, NIA
                              val = val + (One-Occ(t)-Occ(u))*(V_axby(t,u)-V_axby(u,t))*AO%MiniBlocks(t)%MiniAVS(u-NI, km)
                        end do
                        EVVaoTA = EVVaoTA - val**2 / (VV%ATripA(pqT, 1) - AO%MiniBlocks(t)%EigS(km))
                  end do
            end do
            !$omp end parallel do      
          end associate
    end subroutine vvao_contrTA

    
    subroutine vvaa_contr(EVVaa, EVVaaT, EVVaaTA, VV, AA,  pq, pqT, p, q, x0, x1, y0, y1, Occ, V_axby, WorkAA)

          double precision, intent(inout) :: EVVaa, EVVaaT, EVVaaTA
          type(TAC0Block), intent(in) :: VV, AA
          integer, intent(in) :: p, q, pq, pqT
          integer, intent(in) :: x0, x1, y0, y1
          double precision, dimension(:), intent(in) :: Occ
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          double precision, dimension(:), intent(inout) :: workAA
          integer :: t,u, km

          double precision :: val, valA, NumH, NumF

          NumH = One
          NumF = One
          if (p == q) NumH = sqrt(frac12)

          workAA = zero
          
          do km = 1, AA%NDimS
                t = AA%IndNS(1, km)
                u = AA%IndNS(2, km)
                if (t==u) then
                      NumF = sqrt(frac12)
                else
                      NumF = One
                end if
                WorkAA(km) = (One-Occ(t)-Occ(u))*NumH * NUmF * (V_axby(t,u)+V_axby(u,t))
          end do


          do km = 1, AA%MiniBlocks(1)%NDimS
                call real_vw_x(val,  WorkAA(1:AA%NDimS), AA%MiniBlocks(1)%MiniAVS(:, km), AA%NDimS)
                EVVaa = EVVaa - val**2 / (VV%ASing(pq, 1) - AA%MiniBlocks(1)%EigS(km))
          end do

          if (pqT>0)then
                do km = 1, AA%NDimT
                      t = AA%IndNT(1, km)
                      u = AA%IndNT(2, km)
                      
                      WorkAA(km) = (One-Occ(t)-Occ(u))*(V_axby(t,u)-V_axby(u,t))
                end do
                
                
                do km = 1, AA%MiniBlocks(1)%NDimT
                      call real_vw_x(val,  WorkAA(1:AA%NDimT), AA%MiniBlocks(1)%MiniAVT(:, km), AA%NDimT)
                      EVVaaT = EVVaaT - val**2 / (VV%ATrip(pqT, 1) - AA%MiniBlocks(1)%EigT(km))

                      call real_vw_x(valA,  WorkAA(1:AA%NDimT), AA%MiniBlocks(1)%MiniAVTA(:, km), AA%NDimT)
                      EVVaaTA = EVVaaTA - valA**2 / (VV%ATripA(pqT, 1) - AA%MiniBlocks(1)%EigTA(km))

                end do
          end if

          
    end subroutine vvaa_contr

    subroutine aaoo_contr(EAAoo, AZaaoo, AA, OO,  rs, r, s, x0, x1, y0, y1, Occ, V_axby, WorkAA)

          double precision, intent(inout) :: EAAoo
          double precision, dimension(:,:), intent(in) :: AZaaoo
          type(TAC0Block), intent(in) :: AA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          double precision, dimension(:), intent(in) :: Occ
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          double precision, dimension(:), intent(inout) :: workAA
          integer :: t,u, np, pq, p, q, i

          double precision :: val, NumH, NumF

          NumH = One
          NumF = One
          if (r == s) NumH = sqrt(frac12)

          workAA = zero

          do pq = 1, AA%NDimS
                p = AA%IndNS(1, pq)
                q = AA%IndNS(2, pq)
                if (p==q) then
                      NumF = sqrt(frac12)
                else
                      NumF = One
                end if
                WorkAA(pq) = -(One-Occ(p)-Occ(q))*NumH * NUmF * (V_axby(p,q)+V_axby(q,p))
                
          end do

          do np = 1, AA%MiniBlocks(2)%NDimS

                call real_vw_x(val,  WorkAA(1:AA%NDimS), AA%MiniBlocks(2)%MiniAVS(:, np), AA%NDimS)
                EAAoo = EAAoo - val*AZaaoo(np, rs) / (AA%MiniBlocks(2)%EigS(np) - OO%ASing(rs, 1))
          end do

    end subroutine aaoo_contr


    subroutine aaoo_contrT(EAAoo, AZaaoo, AA, OO,  rs, r, s, x0, x1, y0, y1, Occ, V_axby, WorkAA)

          double precision, intent(inout) :: EAAoo
          double precision, dimension(:,:), intent(in) :: AZaaoo
          type(TAC0Block), intent(in) :: AA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          double precision, dimension(:), intent(in) :: Occ
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          double precision, dimension(:), intent(inout) :: workAA
          integer :: t,u, np, pq, p, q, i

          double precision :: val, NumH, NumF

          NumH = One
          NumF = One

          workAA = zero

          do pq = 1, AA%NDimT
                p = AA%IndNT(1, pq)
                q = AA%IndNT(2, pq)

                WorkAA(pq) = -(One-Occ(p)-Occ(q))* (V_axby(p,q)-V_axby(q,p))
                
          end do

          do np = 1, AA%MiniBlocks(2)%NDimT
                call real_vw_x(val,  WorkAA(1:AA%NDimT), AA%MiniBlocks(2)%MiniAVT(:, np), AA%NDimT)
                EAAoo = EAAoo -  val*AZaaoo(np, rs) / (AA%MiniBlocks(2)%EigT(np) - OO%ATrip(rs, 1))
          end do

    end subroutine aaoo_contrT

        subroutine aaoo_contrTA(EAAooTA, AZaaoo, AA, OO,  rs, r, s, x0, x1, y0, y1, Occ, V_axby, WorkAA)

          double precision, intent(inout) :: EAAooTA
          double precision, dimension(:,:), intent(in) :: AZaaoo
          type(TAC0Block), intent(in) :: AA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          double precision, dimension(:), intent(in) :: Occ
          double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
          double precision, dimension(:), intent(inout) :: workAA
          integer :: t,u, np, pq, p, q, i

          double precision :: val, NumH, NumF

          NumH = One
          NumF = One

          workAA = zero

          do pq = 1, AA%NDimT
                p = AA%IndNT(1, pq)
                q = AA%IndNT(2, pq)

                WorkAA(pq) = -(One-Occ(p)-Occ(q))* (V_axby(p,q)-V_axby(q,p))
                
          end do

          do np = 1, AA%MiniBlocks(2)%NDimT
                call real_vw_x(val,  WorkAA(1:AA%NDimT), AA%MiniBlocks(2)%MiniAVTA(:, np), AA%NDimT)
                EAAooTA = EAAooTA -  val*AZaaoo(np, rs) / (AA%MiniBlocks(2)%EigTA(np) - OO%ATripA(rs, 1))
          end do

    end subroutine aaoo_contrTA

    subroutine vaoo_contr(EVAoo, AZ, VA, OO,  rs, r, s, x0, x1, y0, y1, AuxData, Vpq, Vqp, work)

          double precision, intent(inout) :: EVAoo
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq
          double precision, dimension(y0:y1,x0:x1), intent(in) :: Vqp
          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, pq, p, q, i, l0, l1

          double precision :: val, NumH, NumF

          integer :: up, qp
          double precision :: temp, calk
          

          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)
            
            NumH = One
            NumF = One
            if (r == s) NumH = sqrt(frac12)

            work  =zero


            do p = 1, NV

                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  do q=1, NA
                        work(q) = (One-Occ(r)-Occ(s))*(One-Occ(p+NIA)-Occ(q+NI))*NumH * (Vpq(p,q)+Vqp(q,p))
                  end do
                  
                  do q = 1, NA
                        do np = 1, NA
                              pq = np + (p-1)*NA  
                              
                              EVAoo = EVAoo - work(q) * VA%MiniBlocks(p)%MiniAVS(q,np)* AZ(rs, pq)
                              
                        end do
                  end do
            end do
          end associate

    end subroutine vaoo_contr

    subroutine vaoo_contrT(EVAoo,  AZ, VA, OO,  rs, r, s, x0, x1, y0, y1, AuxData, Vpq, Vqp, work)

          double precision, intent(inout) :: EVAoo
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq
          double precision, dimension(y0:y1,x0:x1), intent(in) :: Vqp
          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, pq, p, q, i, l0, l1

          double precision :: val, NumH, NumF

          integer :: up, qp
          double precision :: temp, calk
          

          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)
            
            work  =zero
            do p = 1, NV

                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  do q=1, NA
                        work(q) = (One-Occ(r)-Occ(s))*(One-Occ(p+NIA)-Occ(q+NI))*(Vpq(p,q)-Vqp(q,p))
                  end do
                  
                  do q = 1, NA
                        do np = 1, NA
                              pq = np + (p-1)*NA  
                              
                              EVAoo = EVAoo -  work(q) * VA%MiniBlocks(p)%MiniAVS(q,np)* AZ(rs, pq)
                              
                        end do
                  end do
            end do
          end associate

    end subroutine vaoo_contrT


        subroutine vaoo_contrTA(EVAooA,  AZ, VA, OO,  rs, r, s, x0, x1, y0, y1, AuxData, Vpq, Vqp, work)

          double precision, intent(inout) :: EVAooA
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, OO
          integer, intent(in) :: r, s, rs
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq
          double precision, dimension(y0:y1,x0:x1), intent(in) :: Vqp
          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, pq, p, q, i, l0, l1

          double precision :: val, NumH, NumF

          integer :: up, qp
          double precision :: temp, calk
          

          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)
            
            work  =zero
            do p = 1, NV

                  l0 = 1 +  (p-1)*NA
                  l1 = NA + (p-1)*NA

                  do q=1, NA
                        work(q) = (One-Occ(r)-Occ(s))*(One-Occ(p+NIA)-Occ(q+NI))*(Vpq(p,q)-Vqp(q,p))
                  end do
                  
                  do q = 1, NA
                        do np = 1, NA
                              pq = np + (p-1)*NA  
                              
                              EVAooA = EVAooA -  work(q) * VA%MiniBlocks(p)%MiniAVS(q,np)* AZ(rs, pq)
                              
                        end do
                  end do
            end do
          end associate

    end subroutine vaoo_contrTA


    subroutine vaaa_contr(EVAaa, AZ, VA, AA,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAaa
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AA
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

                        
            NumH = One
            NumF = One

            work  =zero

            do rs = 1, AA%NDimS
                  r = AA%IndNS(1,rs)
                  s = AA%IndNS(2,rs)
                  if (r == s) then
                        NumH = sqrt(frac12)
                  else
                        NumH = One
                  end if

                  work(rs) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))*NumH * (Vpq(r,s)+Vpq(s,r))
            end do
                  
            do k = 1, AA%MiniBlocks(1)%NdimS
                  call real_vw_x(val,  work(1:AA%NDimS), AA%MiniBlocks(1)%MiniAVS(:, k), AA%NDimS)

                  do np = 1, NA
                        np2 = ((p-NIA)-1)*NA+np


                        EVAaa = EVAaa -AZ(k, np2)*VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np)*val

                  end do
            end do

          end associate

    end subroutine vaaa_contr


    subroutine vaaa_contrT(EVAaa, AZ, VA, AA,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAaa
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AA
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

                        

            work  =zero

            do rs = 1, AA%NDimT
                  r = AA%IndNT(1,rs)
                  s = AA%IndNT(2,rs)

                  work(rs) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))* (Vpq(r,s)-Vpq(s,r))
            end do
                  
            do k = 1, AA%MiniBlocks(1)%NdimT
                  ! print*, 'miniavta T', k
                  ! do rs = 1, AA%NDimT
                  !       write(*, '(I5, F15.8)')rs, AA%MiniBlocks(1)%MiniAVT(rs, k)
                  ! end do

                  call real_vw_x(val,  work(1:AA%NDimT), AA%MiniBlocks(1)%MiniAVT(:, k), AA%NDimT)
                  
                  do np = 1, NA
                        np2 = ((p-NIA)-1)*NA+np
                        ! if (abs(val).gt.1.d-4)then
                        !       write(*, '(A10,5I5,  5F15.8)') 'plusz', k, np2, p, q, np,  AZ(k, np2) , VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) , val,  AZ(k, np2) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) * val, EVAaa
                        ! end if
                        EVAaa = EVAaa - AZ(k, np2) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) * val

                  end do
            end do

          end associate

    end subroutine vaaa_contrT

        subroutine vaaa_contrTA(EVAaaA, AZ, VA, AA,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAaaA
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AA
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            work  =zero

            do rs = 1, AA%NDimT
                  r = AA%IndNT(1,rs)
                  s = AA%IndNT(2,rs)

                  work(rs) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))* (Vpq(r,s)-Vpq(s,r))
            end do
                  
            do k = 1, AA%MiniBlocks(1)%NdimTA
                  ! print*, 'miniavta TA', k
                  ! do rs = 1, AA%NDimT
                  !       write(*, '(I5, F15.8)')rs, AA%MiniBlocks(1)%MiniAVTA(rs, k)
                  ! end do
                  call real_vw_x(val,  work(1:AA%NDimT), AA%MiniBlocks(1)%MiniAVTA(:, k), AA%NDimT)

                  do np = 1, NA
                        np2 = ((p-NIA)-1)*NA+np
                        ! if (abs(val).gt.1.d-4)then
                        !       write(*, '(A10,5I5,  5F15.8)') 'plusz',  k, np2, p, q, np,  AZ(k, np2) , VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) , val,  AZ(k, np2) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) * val, EVAaaA
                        ! end if
                        EVAaaA = EVAaaA - AZ(k, np2) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, np) * val
                        
                        
                  end do
            end do

          end associate

    end subroutine vaaa_contrTA



    subroutine vaao2_contr(EVAao, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            ! NumH = One
            ! NumF = One

           !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
           !$omp reduction(+: EVAao) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp = temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(s,r))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do
                        do n = 1, NA

                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n
                              
                              val = AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp

                              EVAao = EVAao  - val
                        end do
                  end do
            end do
            !$omp end parallel do      


          end associate

    end subroutine vaao2_contr

    subroutine vaao2_contrT(EVAao, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

           !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
           !$omp reduction(+: EVAao) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp   =  temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(s,r))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do
                        do n = 1, NA
                              
                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n
                              
                              val = -AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp
                              
                              EVAao = EVAao  -  val
                              
                        end do
                  end do
            end do
            !$omp end parallel do      

          end associate

    end subroutine vaao2_contrT


    subroutine vaao2_contrTA(EVAaoA, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)

          double precision, intent(inout) :: EVAaoA
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

           !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
           !$omp reduction(+: EVAaoA) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp   =  temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(s,r))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do
                        do n = 1, NA
                              
                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n
                              
                              val = -AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp
                              
                              EVAaoA = EVAaoA  -  val
                              
                        end do
                  end do
            end do
            !$omp end parallel do      
          end associate
    end subroutine vaao2_contrTA

    subroutine vaao_contr(EVAao, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)
              
          double precision, intent(inout) :: EVAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

           !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
           !$omp reduction(+: EVAao) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp = temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(r,s))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do

                        do n = 1, NA
                              
                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n

                              val = AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp

                              EVAao = EVAao  - val
                                    
                        end do
                  end do
            end do
            !$omp end parallel do      

          end associate

    end subroutine vaao_contr


    subroutine vaao_contrT(EVAao, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)
              
          double precision, intent(inout) :: EVAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
            !$omp reduction(+: EVAao) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp = temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(r,s))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do
                        do n = 1, NA                                    
                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n

                              val = AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp
                                    
                              EVAao = EVAao  -  val
                                    
                        end do
                  end do
            end do
            !$omp end parallel do      


          end associate

    end subroutine vaao_contrT


        subroutine vaao_contrTA(EVAaoA, AZ, VA, AO,  pq, p, q, x0, x1, y0, y1, AuxData, Vpq, work)
              
          double precision, intent(inout) :: EVAaoA
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: VA, AO
          integer, intent(in) :: pq, p, q
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vpq

          double precision, dimension(:), intent(inout) :: work
          integer :: t,u, np, rs, r, s, i, l0, l1, np2
          double precision :: val, NumH, NumF
          integer :: k, n, km
          double precision :: temp, calk


          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            !$omp parallel do default(shared) private(s, k, temp, r, n, km, np, val) &
            !$omp reduction(+: EVAaoA) collapse(2)
            do s = 1, NI
                  do k = 1, NA
                        temp = zero
                        do r = NI+1, NIA
                              temp = temp + (One-Occ(p)-Occ(q))*(One-Occ(r)-Occ(s))*(Vpq(r,s))*AO%MiniBlocks(s)%MiniAVS(r-NI, k)
                        end do
                        do n = 1, NA                                    
                              km = (s-1)*NA+k
                              np = ((p-NIA)-1)*NA+n

                              val = AZ(km, np) * VA%MiniBlocks(p-NIA)%MiniAVS(q-NI, n) * temp
                                    
                              EVAaoA = EVAaoA  -  val
                                    
                        end do
                  end do
            end do
            !$omp end parallel do      


          end associate

    end subroutine vaao_contrTA


    subroutine aaao_contr(EAAao, AZ, AA, AO,  rs, s, r, x0, x1, y0, y1, AuxData, Vsr, work)
              
          double precision, intent(inout) :: EAAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: AA, AO
          integer, intent(in) :: rs, r, s
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vsr

          double precision, dimension(:), intent(inout) :: work
          integer :: pq, n, km, km2, p, q
          double precision :: val, NumH, NumF



          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

                        
            NumH = One
            NumF = One

            work  = zero

            
            do pq = 1, AA%NDimS
                  p = AA%IndNS(1,pq)
                  q = AA%IndNS(2,pq)
                  if (p == q) then
                        NumH = sqrt(frac12)
                  else
                        NumH = One
                  end if

                  work(pq) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))*NumH * (Vsr(p,q)+Vsr(q,p))
            end do
                  
            do n = 1, AA%MiniBlocks(2)%NdimS
                  call real_vw_x(val,  work(1:AA%NDimS), AA%MiniBlocks(2)%MiniAVS(:, n), AA%NDimS)

                  do km = 1, NA
                        km2 = ((s)-1)*NA+km

                        EAAao = EAAao - AZ(n, km2) * Ao%MiniBlocks(s)%MiniAVS(r-NI, km) * val

                  end do
            end do



          end associate

    end subroutine aaao_contr



    subroutine aaao_contrT(EAAao, AZ, AA, AO,  rs, s, r, x0, x1, y0, y1, AuxData, Vsr, work)
              
          double precision, intent(inout) :: EAAao
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: AA, AO
          integer, intent(in) :: rs, r, s
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vsr

          double precision, dimension(:), intent(inout) :: work
          integer :: pq, n, km, km2, p, q
          double precision :: val, NumH, NumF



          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            work  = zero
            
            do pq = 1, AA%NDimT
                  p = AA%IndNT(1,pq)
                  q = AA%IndNT(2,pq)

                  work(pq) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))*(Vsr(p,q)-Vsr(q,p))
            end do
                  
            do n = 1, AA%MiniBlocks(2)%NdimT
                  call real_vw_x(val,  work(1:AA%NDimT), AA%MiniBlocks(2)%MiniAVT(:, n), AA%NDimT)

                  do km = 1, NA
                        km2 = ((s)-1)*NA+km

                        EAAao = EAAao + AZ(n, km2) * AO%MiniBlocks(s)%MiniAVS(r-NI, km) * val

                  end do
            end do



          end associate

    end subroutine aaao_contrT

    
    subroutine aaao_contrTA(EAAaoA, AZ, AA, AO,  rs, s, r, x0, x1, y0, y1, AuxData, Vsr, work)
              
          double precision, intent(inout) :: EAAaoA
          double precision, dimension(:,:), intent(in) :: AZ
          type(TAC0Block), intent(in) :: AA, AO
          integer, intent(in) :: rs, r, s
          integer, intent(in) :: x0, x1, y0, y1
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(x0:x1,y0:y1), intent(in) :: Vsr

          double precision, dimension(:), intent(inout) :: work
          integer :: pq, n, km, km2, p, q
          double precision :: val, NumH, NumF



          associate (Occ=> AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)

            work  = zero
            
            do pq = 1, AA%NDimT
                  p = AA%IndNT(1,pq)
                  q = AA%IndNT(2,pq)

                  work(pq) = (One-Occ(r)-Occ(s))*(One-Occ(p)-Occ(q))*(Vsr(p,q)-Vsr(q,p))
            end do
                  
            do n = 1, AA%MiniBlocks(2)%NdimT
                  call real_vw_x(val,  work(1:AA%NDimT), AA%MiniBlocks(2)%MiniAVTA(:, n), AA%NDimT)

                  do km = 1, NA
                        km2 = ((s)-1)*NA+km

                        EAAaoA = EAAaoA + AZ(n, km2) * AO%MiniBlocks(s)%MiniAVS(r-NI, km) * val

                  end do
            end do

          end associate

    end subroutine aaao_contrTA



end module thc_energy
