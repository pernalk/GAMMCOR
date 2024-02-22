module sapt_resp

use types
use timing
use tran
use abmat
use abfofo
use ab0fofo
use sapt_utils
use sapt_inter

implicit none

contains

! subroutine sapt_2trdm(Flags,A,B,NBasis)
! implicit none

! type(FlagsData)    :: Flags
! type(SystemBlock)  :: A,B
! integer,intent(in) :: NBasis
! integer :: pq,ip,iq,ir,is,it,exc,i,j
! integer :: dimOA,dimOB
! integer :: nstatesA,nstatesB
! double precision,allocatable :: RDM2A(:,:,:,:),RDM2B(:,:,:,:),TRDM1(:,:),TRDM1TEST(:,:),RDM1TEST(:,:)

! write(lout,'(/1x,a)') 'entering 2trdm...'
! dimOA  = A%num0+A%num1
! dimOB  = B%num0+B%num1
! nstatesA = A%NStates
! nstatesB = B%NStates
! allocate(RDM2A(dimOA,dimOA,dimOA,dimOA),&
!          RDM2B(dimOB,dimOB,dimOB,dimOB))

! RDM2A(:,:,:,:) = A%rdm24(A%iref1,:,:,:,:)
! RDM2B(:,:,:,:) = B%rdm24(B%iref1,:,:,:,:)
! !print*,"RDM2A"
! !print*,RDM2A


! !allocate(A%trdm24(nstatesA,NBasis,NBasis,NBasis,NBasis))
! !allocate(B%trdm24(nstatesB,NBasis,NBasis,NBasis,NBasis))

! allocate(A%trdm24(4,NBasis,NBasis,NBasis,NBasis))
! allocate(B%trdm24(4,NBasis,NBasis,NBasis,NBasis))

! allocate(TRDM1(A%NDimX,A%NDimX))
! allocate(TRDM1TEST(A%NDimX,A%NDIMX))
! TRDM1=0d0


! do exc=1,4 
!    A%trdm24(exc,:,:,:,:)=0d0
! enddo


! print*, 'A EigX',norm2(A%EigX)
! print*, 'A EigY',norm2(A%EigY)
! print*, 'B EigX',norm2(B%EigX)
! print*, 'B EigY',norm2(B%EigY)
! print*,'NDimX',A%NDimX
! !print*,'A%EigX'
! do pq=1,A%NDimX
!    if (A%EigX(2*A%NDimX+pq)/= 0d0) then
! !   print *,A%EigX(2*A%NDimX+pq)
!    endif
! enddo
! ! 1 OK OK 
! do exc=1,4
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       do ir = 1,dimOA
!          do is = 1,dimOA
!             do it =1,dimOA
!             A%trdm24(exc,ip,is,ir,it) =  A%trdm24(exc,ip,is,ir,it)+A%EigX((exc-1)*A%NDimX+pq)*RDM2A(iq,is,ir,it)
! !            print*, A%EigX((exc-1)*A%NDimX+pq)
!             enddo
!          enddo
!       enddo
!    enddo
!  ! 2  OK 
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       do ir = 1,dimOA
!          do is = 1,dimOA
!             do it =1,dimOA
!             A%trdm24(exc,ir,is,ip,it) =  A%trdm24(exc,ir,is,ip,it)+A%EigX((exc-1)*A%NDimX+pq)*RDM2A(ir,is,iq,it)
!             enddo
!          enddo
!       enddo
!    enddo
! !3   OK
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       do ir = 1,dimOA
!          do is = 1,dimOA
!             do it =1,dimOA
!                if (ip <= dimOA) then
!                   A%trdm24(exc,ir,iq,is,it) =  A%trdm24(exc,ir,iq,is,it)-A%EigX((exc-1)*A%NDimX+pq)*RDM2A(ir,ip,is,it)
!             endif
!             enddo
!          enddo
!       enddo
!    enddo
! !4   OK
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       do ir = 1,dimOA
!          do is = 1,dimOA
!             do it =1,dimOA
!                if(ip<=dimOA) then
!                   A%trdm24(exc,ir,it,is,iq) =  A%trdm24(exc,ir,it,is,iq)-A%EigX((exc-1)*A%NDimX+pq)*RDM2A(ir,it,is,ip)
!                endif   
!             enddo
!          enddo
!       enddo
!    enddo



! !5 OK
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
!          do ir = 1,dimOA
!             do is = 1,dimOA
!                do it =1,dimOA
!                   if(ip <= dimOA) then
!                      A%trdm24(exc,iq,is,ir,it) =  A%trdm24(exc,iq,is,ir,it)+A%EigY((exc-1)*A%NDimX+pq)*RDM2A(ip,is,ir,it)
!                   endif   
!                enddo
!             enddo
!          enddo
!       enddo


! !6  OK    
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
!          do ir = 1,dimOA
!             do is = 1,dimOA
!                do it =1,dimOA
!                      if(ip<= dimOA) then
!                         A%trdm24(exc,ir,is,iq,it) =  A%trdm24(exc,ir,is,iq,it)+A%EigY((exc-1)*A%NDimX+pq)*RDM2A(ir,is,ip,it)
!                      endif     
!                enddo
!             enddo
!          enddo
!       enddo
!  !7   OK  
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
!          do ir = 1,dimOA
!             do is = 1,dimOA
!                do it =1,dimOA
!                A%trdm24(exc,ir,ip,is,it) =  A%trdm24(exc,ir,ip,is,it)-A%EigY((exc-1)*A%NDimX+pq)*RDM2A(ir,iq,is,it)
!                enddo
!             enddo
!          enddo
!       enddo
! !8      OK
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
!          do ir = 1,dimOA
!             do is = 1,dimOA
!                do it =1,dimOA
!                A%trdm24(exc,ir,it,is,ip) =  A%trdm24(exc,ir,it,is,ip)-A%EigY((exc-1)*A%NDimX+pq)*RDM2A(ir,it,is,iq)
!                enddo
!             enddo
!          enddo
!       enddo

! enddo
! print*,"2TRDM"
! do ip=1,10
!    do iq = 1,10
!       do ir = 1,10
!          do is =1,10
!             if(abs(A%trdm24(1,ip,iq,ir,is))>1d-2) then
!                print*,ip,iq,ir,is,A%trdm24(1,ip,iq,ir,is)
!             endif
!          enddo
!       enddo
!    enddo
! enddo
! print*,'Firs few eigenvalues',A%Eig(1),A%Eig(2),A%Eig(3)
! print*,'Dla pewności',A%rdm1(1,:)

! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!    print*,"ip,iq",ip,",",iq
!          TRDM1(ip,iq)=TRDM1(ip,iq)-(A%rdm1(1,ip)-A%rdm1(1,iq))*A%EigX(0*A%NDimX+pq)
!          TRDM1(iq,ip)=TRDM1(iq,ip)-(A%rdm1(1,iq)-A%rdm1(1,ip))*A%EigY(0*A%NDimX+pq)
! enddo
! print*,"1TRDM z X i Y"
! do j=1,10
!    write(lout,'(*(f12.6))') (TRDM1(i,j),i=1,10)
! enddo
! print*, "1TRDM z X i Y norma : ", norm2(TRDM1)
! do i=1,A%NDimX
!    do j=1,A%NDimX
!    if(abs(TRDM1(i,j))>1d-2) then
!       print *,"TRDM1 duze dla :",i, j, "wynosi", TRDM1(i,j)
!    endif
!    enddo
! enddo   
! TRDM1TEST=0d0
! do ip=1,NBasis
!    do ir=1,NBAsis
!       do iq=1,NBasis
!          TRDM1TEST(ip,ir)= TRDM1TEST(ip,ir)+A%trdm24(1,ip,ir,iq,iq)
!       enddo
!    enddo
! enddo

! print*,"1TRDM odtworzny z TRDM"
! do j=1,10
!    write(lout,'(*(f12.6))') (TRDM1TEST(i,j),i=1,10)
! enddo

! print*,"Norma 2TRDM"
! print*,norm2(A%trdm24(1,:,:,:,:))
! print*,"A czy działa dla zwykłego 2RDM?"
! RDM1TEST=0d0
! allocate(RDM1TEST(dimOA,dimOA))
! do ip=1,dimOA
!    do ir=1,dimOA
!       do iq=1,dimOA
!          RDM1TEST(ip,ir)= RDM1TEST(ip,ir)+A%rdm24(1,ip,ir,iq,iq)
!       enddo
!    enddo
! enddo
! print*,"1RDM odtworzny z 2RDM"
! do j=1,10
!    write(lout,'(*(f12.6))') (RDM1TEST(i,j),i=1,10)
! enddo
! print*,"AOCC",A%Occ



! block
!    double precision :: DipXao(NBasis**2),DipYao(NBasis**2),DipZao(NBasis**2)
!    double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
!    double precision:: TSDipXYZ(3)
!    character(:),allocatable     :: mname,dipfile,mofile
   
!       mname  = 'A'
!       dipfile= "DIP_A"
!       mofile = 'MOLPRO_A.MOPUN'


!    call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

! !print*, 'DipX-AO:',norm2(DipXao)
! !print*, 'DipY-AO:',norm2(DipYao)
! !print*, 'DipZ-AO:',norm2(DipZao)

!    call unpack_sym_molpro(DipXao,dipfile,NBasis)
!    call unpack_sym_molpro(DipYao,dipfile,NBasis)
!    call unpack_sym_molpro(DipZao,dipfile,NBasis)

!    call tran2MO(DipXao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipX,NBasis)
!    call tran2MO(DipYao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipY,NBasis)
!    call tran2MO(DipZao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipZ,NBasis)
!    TSDipXYZ = 0d0
!    do j=1,NBasis
!       do i=1,NBasis
!          TSDipXYZ(1) = TSDipXYZ(1) - 2d0*TRDM1(j,i)*DipX(i,j)
!          TSDipXYZ(2) = TSDipXYZ(2) - 2d0*TRDM1(j,i)*DipY(i,j)
!          TSDipXYZ(3) = TSDipXYZ(3) - 2d0*TRDM1(j,i)*DipZ(i,j)
!       enddo
!    enddo
!    print *,'X,Y,Z moments from X and Y matrix'
!    print *,TSDipXYZ
! end block

! deallocate(TRDM1TEST)
! end subroutine sapt_2trdm

subroutine calc_resp_casgvb(Mon,MO,Flags,NBas,EChck)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags
double precision :: MO(:)
integer :: NBas
logical :: EChck
integer :: NSq,NInte1,NInte2
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:)
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:), &
                                 CMAT(:),EMAT(:),EMATM(:), &
                                 DMAT(:),DMATK(:), &
                                 EigVecR(:), Eig(:), &
                                 ABPlusT(:), ABMinT(:)
double precision, allocatable :: Eig0(:),Eig1(:),EigY0(:),EigY1(:)
double precision,allocatable  :: workSq(:,:)

integer :: dimOcc
integer :: i,j,k,l,ii,jj,ij,ij1,ione,itwo
integer :: ip,iq,pq
integer :: iunit
integer :: DimEx
integer :: INegExcit
!test
integer :: iter

double precision :: ACAlpha
double precision :: ECASSCF,ETot,ECorr

character(8) :: label
character(:),allocatable :: onefile,twofile,propfile,rdmfile
character(:),allocatable :: twojfile,twokfile
character(:),allocatable :: propfile0,propfile1
character(:),allocatable :: y01file,xy0file
character(:),allocatable :: abpm0file
character(:),allocatable :: abfile,testfile

double precision,external  :: ddot
double precision,parameter :: One = 1d0, Half = 0.5d0
double precision,parameter :: SmallE=0d0,BigE=1.D20
!! test Cmat
!double precision :: EGOne,NGOcc
!double precision,allocatable :: Pmat(:,:)
integer :: nblk
type(EblockData) :: A0blockIV
type(EblockData),allocatable :: A0block(:)
! test Calpha
integer :: NGOcc
double precision :: OmI
double precision,allocatable :: COMTilde(:)
double precision,allocatable :: ABPlus0(:,:),ABMin0(:,:)

! test iterative
iter = 1

! set filenames
if(Mon%Monomer==1) then
   onefile    = 'ONEEL_A'
   twofile    = 'TWOMOAA'
   twojfile   = 'FFOOAA'
   twokfile   = 'FOFOAA'
   propfile   = 'PROP_A'
   propfile0  = 'PROP_A0'
   propfile1  = 'PROP_A1'
   y01file    = 'Y01_A'
   xy0file    = 'XY0_A'
   abpm0file  = 'A0BLK_A'
   rdmfile    = 'rdm2_A.dat'
   abfile     = 'ABMAT_A'
   testfile   = 'A0MAT_A'
elseif(Mon%Monomer==2) then
   onefile    = 'ONEEL_B'
   twofile    = 'TWOMOBB'
   twojfile   = 'FFOOBB'
   twokfile   = 'FOFOBB'
   propfile   = 'PROP_B'
   propfile0  = 'PROP_B0'
   propfile1  = 'PROP_B1'
   y01file    = 'Y01_B'
   xy0file    = 'XY0_B'
   abpm0file  = 'A0BLK_B'
   rdmfile    = 'rdm2_B.dat'
   abfile     = 'ABMAT_B'
   testfile   = 'A0MAT_B'
endif

! set dimensions
NInte1 = NBas*(NBas+1)/2
NInte2 = 1
if(Mon%TwoMoInt==TWOMO_INCORE) NInte2 = NInte1*(NInte1+1)/2

allocate(work1(NBas**2),work2(NBas**2),XOne(NInte1),URe(NBas,NBas))
!if(Mon%TwoMoInt==TWOMO_INCORE) allocate(TwoMO(NInte2))
allocate(TwoMO(NInte2))

URe = 0d0
do i=1,NBas
   URe(i,i) = 1d0
enddo

! read 1-el
call get_1el_h_mo(XOne,MO,NBas,onefile)

! INCORE: load 2-el integrals
if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBas,NInte2)

! response is obtained either for GVB or CASSCF functions

ACAlpha=One
! GVB
if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

  allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2), &
           EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

! ACAlpha=sqrt(2d0)/2.d0
!   allocate(ABPlusT(Mon%NDim**2),ABMinT(Mon%NDim**2))
!   call ACABMAT0(ABPlusT,ABMinT,URe,Mon%Occ,XOne,TwoMO, &
!   call ACABMAT0(ABPlus,ABMin,URe,Mon%Occ,XOne,TwoMO, &
!                 NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,ACAlpha,1)

  select case(Mon%TwoMoInt)
  case(TWOMO_FOFO)
     call ACABMAT0_FOFO(ABPlus,ABMin,URe,Mon%Occ,XOne,&
                   Mon%IndN,Mon%IndX,Mon%IGem,Mon%CICoef,&
                   Mon%NAct,Mon%NELE,NBas,Mon%NDim,Mon%NDimX,NInte1,Mon%NGem,&
                   twofile,twojfile,twokfile,0,ACAlpha,1)

  case(TWOMO_FFFF,TWOMO_INCORE)

     call ACABMAT0_mithap(ABPlus,ABMin,URe,Mon%Occ,XOne,&
                   Mon%IndN,Mon%IndX,Mon%IGem,Mon%CICoef,&
                   NBas,Mon%NDim,Mon%NDimX,NInte1,Mon%NGem,&
                   twofile,Flags%ISAPT,ACAlpha,1)

  end select

  ! reduce dim
  !call reduce_dim('AB',ABPlus,ABMin,Mon)

  ! TESTY
  !call reduce_dim('AB',ABPlusT,ABMinT,Mon)
  ! write(*,*) 'TEST_PLUS',norm2(ABPlus),norm2(ABPlusT(1:Mon%NDimX**2))
  ! write(*,*) 'TEST_MIN',norm2(ABMin),norm2(ABMinT(1:Mon%NDimX**2))
  ! write(*,*) 'ABPLUS: ',norm2(ABPlus(1:Mon%NDimX**2)-ABPlusT(1:Mon%NDimX**2))
  ! write(*,*) 'ABMIN: ',norm2(ABMin(1:Mon%NDimX**2)-ABMinT(1:Mon%NDimX**2))

  ! do j=1,Mon%NDimX**2
  !   write(*,*) j,ABPlus(j),ABPlusT(j)
  !enddo

  !deallocate(ABMinT,ABPlusT)

  EigVecR = 0
  Eig = 0

  !call ERPASYMM(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
  ! for exch-disp
  allocate(Mon%EigY(Mon%NDimX**2),Mon%EigX(Mon%NDimX**2),&
           Mon%Eig(Mon%NDimX))
  call ERPASYMMXY(Mon%EigY,Mon%EigX,Mon%Eig,ABPlus,ABMin,&
                  Mon%Occ,Mon%IndN,Mon%NDimX,NBas)

  Eig = Mon%Eig
  do j=1,Mon%NDimX
     do i=1,Mon%NDimX
        ip = Mon%IndN(1,i)
        iq = Mon%IndN(2,i)
        EigVecR((j-1)*Mon%NDimX+i)=(Mon%CICoef(ip)-Mon%CICoef(iq))&
                                  *(Mon%EigY((j-1)*Mon%NDimX+i)-Mon%EigX((j-1)*Mon%NDimX+i))
     enddo
  enddo

  if(Mon%TwoMoInt==TWOMO_INCORE) then
     if(EChck) then
        write(LOUT,'(/,1x,a)') 'ERPA-GVB ENERGY CHECK REQUESTED:'
        call EneERPA(ETot,ECorr,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,&
             Mon%Occ,XOne,Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
     endif
  elseif(Mon%TwoMoInt==TWOMO_FFFF) then
      call EneGVB_FFFF(ETot,URe,Mon%Occ,Mon%CICoef,XOne, &
                   Mon%IGem,Mon%IndN,NBas,NInte1,twofile,Mon%NDimX,Mon%NGem)
!     here: need some intermediate steps to get EigVY2,AMAT,BMAT!
!     see: SndOrder subroutine
!      call ECorrAC0GVB_mithap(ECorr0,ECorr,AMAT,BMAT,ABPLUS,
!                             EigVY2,Mon%Occ,Mon%CICoef,Eig,
!                             IndP,Mon%IndN,IndX,Mon%IGem,
!                             twofile,Mon%NDim,Mon%NDimX,NGem,NBas)
  endif

  ! UNCOUPLED
  !
  ! HERE WE ACTUALLY NEED TO PREPARE PROCEDURES THAT 
  ! EXPLOIT THE DIAGONAL STRUCTURE OF ABMAT IN GVB!
  !
  !allocate(EigY0(Mon%NDimX**2),Eig0(Mon%NDimX),Eig1(Mon%NDimX))

  !EigY0 = 0
  !Eig0 = 0
  !Eig1 = 0

  !call Y01GVB(TwoMO,Mon%Occ,URe,XOne, &
  !     EigY0,Eig0,Eig1, &
  !     Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2)

  !! dump uncoupled response
  !call writeresp(EigY0,Eig0,propfile0)
  !if(Flags%IFlag0==0) then
  !   call writeEval(Eig1,propfile1)
  !endif

  !deallocate(Eig1,Eig0,EigY0)

  ! CAS-SCF
  elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then

  allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

  ECASSCF = 0
!
! test semicoupled
!  ACAlpha = 0.000000001
! if(Mon%Monomer==1) then
!   ACAlpha=0.010
! else
!   ACAlpha=0.0000001
! endif
! KP 31.01.2021 instability test
!   if(Mon%Monomer==1) then
!   IF(COMMAND_ARGUMENT_COUNT().Ne.0) THEN
!   CALL GET_COMMAND_ARGUMENT(1,label)
!   READ(label,*)ACAlpha
!   ENDIF
!   print*,'*********************'
!   print*,'*** MONOMER ALPHA ***',Mon%Monomer,ACAlpha
!   endif

  !ACAlpha=sqrt(2d0)/2d0
  ACAlpha=1d-12
  Print*, 'UNCOUPLED,ACAlpha',ACAlpha

  !ACAlpha=0.953089922969332
  !print*, 'ACAlpha',ACAlpha
  select case(Mon%TwoMoInt)
  case(TWOMO_FOFO)
     call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                 Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                 NInte1,twojfile,twokfile,Flags%ICholesky,ACAlpha,.false.)
  case(TWOMO_FFFF)

     call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                 Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                  NInte1,twofile,ACAlpha,.false.)
  case(TWOMO_INCORE)

     call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
                 Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)

     call EKT(URe,Mon%Occ,XOne,TwoMO,NBas,NInte1,NInte2)

  end select
  !write(LOUT,'(/,1x,a,f16.8,a,1x,f16.8)') 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)
  print*, 'ABPlus',norm2(ABPlus)
  print*, 'ABMin',norm2(ABMin)

  EigVecR = 0
  Eig = 0
  ! HERE WORTH CHECKING SOME MORE FOR EXCITED STATE SELECTION...
  if(Mon%InSt(1,1) > 0) then
     write(lout,'(/1x,a,i2,a,i1)') 'Calculating reponse for state:',Mon%InSt(1,1),'.',Mon%InSt(2,1)
  else
     write(lout,'(/1x,a)') 'Calculating reponse for the ground state'
  endif

  ! MH 1 Dec 2020: use of ERPAVECTRANS is temporarily disabled
  !                due to poor quality of deexcitation energies from ERPA (in general)
  !if(Mon%InSt(1,1)<0.or.Mon%InSt(1,1)==1) then

  if(iter==1) then
     ! dump matrices for iterative C-ERPA
     open(newunit=iunit,file=abfile,form='unformatted')
     write(iunit) ABPlus
     write(iunit) ABMin
     close(iunit)
  endif

!! this is only for testing e2disp_CAlphaTilde_full
!if(TWOMO_FOFO) then
!   if(iter==1) then
!      allocate(ABPlus0(Mon%NDimX,Mon%NDImX),&
!               ABMin0(Mon%NDimX,Mon%NDimX))
!      call AB_CAS_FOFO(ABPlus0,ABMin0,ECASSCF,URe,Mon%Occ,XOne, &
!                  Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
!                  NInte1,twojfile,twokfile,Flags%ICholesky,1d-9,.false.)
!      ! dump matrices for iterative C-ERPA
!      open(newunit=iunit,file=testfile,form='unformatted')
!      write(iunit) ABPlus0
!      write(iunit) ABMin0
!      close(iunit)
!      deallocate(ABMin0,ABPlus0)
!   endif
!endif

  write(lout,'(1x,a/)') 'Solving symmetric eigenvalue equation...'

  allocate(Mon%EigY(Mon%NDimX**2),Mon%EigX(Mon%NDimX**2),&
           Mon%Eig(Mon%NDimX))
  call ERPASYMMXY(Mon%EigY,Mon%EigX,Mon%Eig,ABPlus,ABMin,&
                  Mon%Occ,Mon%IndN,Mon%NDimX,NBas)

  !print*, 'EigY:',norm2(Mon%EigY),norm2(Mon%EigX)

  Eig = Mon%Eig
  do j=1,Mon%NDimX
     do i=1,Mon%NDimX
        ip = Mon%IndN(1,i)
        iq = Mon%IndN(2,i)
        EigVecR((j-1)*Mon%NDimX+i)=(Mon%CICoef(ip)-Mon%CICoef(iq))&
                                  *(Mon%EigY((j-1)*Mon%NDimX+i)-Mon%EigX((j-1)*Mon%NDimX+i))
     enddo
  enddo
  !print*, 'EigVecR',norm2(EigVecR)
  !print*, 'Eig    ',norm2(Mon%Eig)

  ! test for e2ind_hf
  !allocate(Mon%PP(Mon%NDimX**2))
  !call AB_CAS_FOFO(ABPlus,Mon%PP,ECASSCF,URe,Mon%Occ,XOne, &
  !              Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
  !              NInte1,twojfile,twokfile,Flags%ICholesky,1d0,.true.)

  !print*, 'Mon%PP',norm2(Mon%PP)

  !do i=1,Mon%NDimX
  !   print*, i, Eig(i)
  !enddo

  !! TEST COUPLED ANDREAS
  ! print*, 'CPLD-ANDREAS:'
  ! Mon%EigX = EigVecR
  ! Mon%EigY = 0


!!    else  ! MH 1 Dec 2020: disable ERPAVECTRANS
!
!      allocate(Mon%EigY(Mon%NDimX**2),Mon%EigX(Mon%NDimX**2),&
!               Mon%Eig(Mon%NDimX))
!       Mon%EigX = 0
!       Mon%EigY = 0
!   !   call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
!
!       write(lout,'(1x,a)') 'Solving full (non-symmetric) eigenvalue equation...'
!       call ERPAVECTRANS(Mon%EigY,Mon%EigX,Mon%Eig,ABPlus,ABMin,&
!                         Mon%Occ,Mon%IndN,Mon%NDimX,NBas)
!
!      Eig = Mon%Eig
!      do j=1,Mon%NDimX
!         do i=1,Mon%NDimX
!            ip = Mon%IndN(1,i)
!            iq = Mon%IndN(2,i)
!            EigVecR((j-1)*Mon%NDimX+i)=(Mon%CICoef(ip)-Mon%CICoef(iq))&
!                                      *(Mon%EigY((j-1)*Mon%NDimX+i)-Mon%EigX((j-1)*Mon%NDimX+i))
!         enddo
!      enddo
!!
!!   endif

  !print*, 'STOP AFTER RESPONSE FOR CHECKS!'
  !Stop

  !print*, 'Entering ERPAVEC...'
  !call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
  !

 !print*, 'EigVecR',norm2(EigVecR)
 !do i=1,size(Eig)
 !   print*, i,Eig(i)
 !enddo

  if(EChck) then
     ECorr=0
     select case(Mon%TwoMoInt)
     case(TWOMO_FOFO)
        call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Mon%Occ, &
                             Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
                             Mon%NDimX,NBas,twokfile,Flags%ICholesky)
     case(TWOMO_FFFF)
        call ACEneERPA_FFFF(ECorr,EigVecR,Eig,Mon%Occ, &
                             Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
                             Mon%NDimX,NBas,twofile)
     case(TWOMO_INCORE)
        call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
                       Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
     end select
     ECorr=Ecorr*0.5d0

     Mon%ECASSCF = ECASSCF+Mon%PotNuc
     write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
          ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr
  else
     write(LOUT,'(1x,a,5x,f15.8)') "CASSCF Energy           ", ECASSCF+Mon%PotNuc
  endif

  !! snippet for testing Cmat
  !EGOne = 0
  !NGOcc = 0
  !ECorr = 0

  if(Flags%ICholesky==1) then
     open(newunit=iunit,file='cholvecs',form='unformatted')
     write(iunit) Mon%NChol
     write(iunit) Mon%FF
     close(iunit)

     ! sub-snippet for testing Pmat
     ! this will not work now, shit -- Pmat!
     ! uncomment this for e2disp_Cmat_Chol / e2disp_Cmat_Chol_proj
     !allocate(Mon%Pmat(Mon%NDimX,Mon%NDimX))
     !call Project_DChol(Mon%PMat,Mon%IndN,NBas,Mon%NDimX)

     !call CIter_FOFO(ECorr,ACAlpha,XOne,URe,Mon%Occ,EGOne,NGOcc,&
     !                Mon%IGem,Mon%NAct,Mon%INAct,Mon%NELE,NBas,NInte1, &
     !                Mon%NDim,Mon%NGem,Mon%IndAux,Mon%IndN,Mon%IndX,Mon%NDimX,&
     !                twojfile,twokfile)
     !deallocate(Pmat)
  endif

  ! UNCOUPLED

  allocate(Mon%IndNT(2,Mon%NDim))
  Mon%IndNT=0
  do i=1,Mon%NDim
     Mon%IndNT(1,i) = Mon%IndN(1,i)
     Mon%IndNT(2,i) = Mon%IndN(2,i)
  enddo

  select case(Mon%TwoMoInt)
  case(TWOMO_FOFO)
     call Y01CAS_FOFO(Mon%Occ,URe,XOne,ABPlus,ABMin,ECASSCF, &
            propfile0,propfile1, &
            y01file,xy0file,     &
            Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
            NBas,Mon%NDimX,NInte1,Mon%NoSt,twofile,twojfile,twokfile,&
            Flags%IFlag0,Flags%ICholesky)

  if(Flags%ICholesky==1) then
     nblk = 1 + NBas - Mon%NAct
     allocate(A0Block(nblk))
     ! maybe just include blocks in Mon%...?
     call AC0BLOCK(Mon%Occ,URe,XOne, &
          Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
          NBas,Mon%NDimX,NInte1,twojfile,twokfile, &
          A0BlockIV,A0Block,nblk,abpm0file,1)
  endif
  case(TWOMO_FFFF)
     call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
            propfile0,propfile1, &
            y01file, &
            Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
            NBas,Mon%NDim,NInte1,Mon%NoSt,twofile,Flags%IFlag0)
  case(TWOMO_INCORE)

     allocate(EigY0(Mon%NDimX**2),EigY1(Mon%NDimX**2),&
              Eig0(Mon%NDimX),Eig1(Mon%NDimX))

      EigY0 = 0
      EigY1 = 0
      Eig0  = 0
      Eig1  = 0

      !! silly AC0 energy test
      ! Call AC0CAS(ECorr,ETot,TwoMO,Mon%Occ,URe,XOne,&
      !             ABPLUS,ABMIN,EigY0,Eig0, &
      !             Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2)
      !       write(*,*) 'ECORR:',ETot,ECorr
      call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
             EigY0,EigY1,Eig0,Eig1, &
             Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0,xy0file)

     !print*, 'UNCOUPLED--PODMIANKA!'
     !Mon%EigX=0
     !Mon%EigY=0
     !open(newunit=iunit,file=testfile,form='UNFORMATTED',&
     !     access='SEQUENTIAL',status='OLD')
     !read(iunit) Mon%EigX
     !read(iunit) Mon%EigY
     !close(iunit)

     !Eig = Eig0
     !Mon%Eig = Eig0
     !do j=1,Mon%NDimX
     !   do i=1,Mon%NDimX
     !      ip = Mon%IndN(1,i)
     !      iq = Mon%IndN(2,i)
     !      EigVecR((j-1)*Mon%NDimX+i)=(Mon%CICoef(ip)-Mon%CICoef(iq))&
     !                                *(Mon%EigY((j-1)*Mon%NDimX+i)-Mon%EigX((j-1)*Mon%NDimX+i))
     !   enddo
     !enddo

     !print*, 'test-podmianka',norm2(Mon%EigY),norm2(Mon%EigX)
     !print*, 'test-podmiank2',norm2(EigVecR)
     !print*, Mon%EigY(1),Mon%EigX(1)
     !do j=1,10
     !do i=3,3
     !   print*, i,j,Mon%EigY((j-1)*Mon%NDimX+i)
     !enddo
     !enddo

      ! dump uncoupled response
      if(Mon%TwoMoInt==TWOMO_INCORE) then
         call writeresp(EigY0,Eig0,propfile0)
         if(Flags%IFlag0==0) then
            call writeresp(EigY1,Eig1,propfile1)
         endif
      endif

      deallocate(Eig1,Eig0,EigY1,EigY0)

  end select

  !! test CAlpha
!   NGOcc = 0
!   OmI = 1.08185673347417
!   call C_AlphaExpand(COMTilde,OmI,XOne,URe,Mon%Occ,NGOcc,&
!                     Mon%IGem,Mon%NAct,Mon%INAct,NBas,NInte1,&
!                     Mon%NDim,Mon%NGem,&
!                     Mon%IndAux,Mon%IndN,Mon%IndX,Mon%NDimX,&
!                     twojfile,twokfile,xy0file,abpm0file)
!   print*, 'COMTilde',norm2(COMTilde)
!   deallocate(COMTilde)

   ! WARNING: MH, September 2021
   ! in previous versions of the code all PINO cases (ISERPA==2)
   ! i.e., TD-GVB, CAS-LinearResponse and PINO-FCI for 2-electron systems
   ! were included here: they were now moved to calc_resp_pino  

endif

! dump response
 call writeresp(EigVecR,Eig,propfile)

 deallocate(work1,work2,XOne,URe)
 !if(Mon%TwoMoInt==1) deallocate(TwoMO)
 deallocate(TwoMO)
 deallocate(ABPlus,ABMin,EigVecR,Eig)

end subroutine calc_resp_casgvb

subroutine calc_ab_cas(Mon,MO,Flags,NBas)
!
!for Cholesky !
! should work with e2disp_CAlpha
!
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags
double precision :: MO(:)
integer :: NBas
integer :: NSq,NInte1,NInte2
double precision, allocatable :: XOne(:),TwoMO(:)
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:)
double precision, allocatable :: work1(:),work2(:)

integer :: dimOcc
integer :: i,j,k,l,ii,jj,ij,ij1,ione,itwo
integer :: ip,iq,pq
integer :: iunit

double precision :: ACAlpha
double precision :: ECASSCF

character(8) :: label
character(:),allocatable :: onefile,twofile
character(:),allocatable :: twojfile,twokfile
character(:),allocatable :: abfile
character(:),allocatable :: abpm0file,xy0file

double precision,parameter :: One = 1d0, Half = 0.5d0

integer :: nblk
type(EblockData)             :: A0blockIV
type(EblockData),allocatable :: A0block(:)

! set filenames
if(Mon%Monomer==1) then
   onefile    = 'ONEEL_A'
   twofile    = 'TWOMOAA'
   twojfile   = 'FFOOAA'
   twokfile   = 'FOFOAA'
   abfile     = 'ABMAT_A'
   abpm0file  = 'A0BLK_A'
   xy0file    = 'XY0_A'
elseif(Mon%Monomer==2) then
   onefile    = 'ONEEL_B'
   twofile    = 'TWOMOBB'
   twojfile   = 'FFOOBB'
   twokfile   = 'FOFOBB'
   abfile     = 'ABMAT_B'
   abpm0file  = 'A0BLK_B'
   xy0file    = 'XY0_B'
endif

! set dimensions
NInte1 = NBas*(NBas+1)/2
NInte2 = 1
if(Mon%TwoMoInt==TWOMO_INCORE) NInte2 = NInte1*(NInte1+1)/2

allocate(work1(NBas**2),work2(NBas**2),XOne(NInte1),URe(NBas,NBas))
allocate(TwoMO(NInte2))

URe = 0d0
do i=1,NBas
   URe(i,i) = 1d0
enddo

! read 1-el
call get_1el_h_mo(XOne,MO,NBas,onefile)

! INCORE: load 2-el integrals
if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBas,NInte2)

! get AB matrices for CAS function

ACAlpha=One
allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2))

ECASSCF = 0

select case(Mon%TwoMoInt)
case(TWOMO_FOFO)
   call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
               NInte1,twojfile,twokfile,Flags%ICholesky,ACAlpha,.false.)
case(TWOMO_FFFF)

   call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                NInte1,twofile,ACAlpha,.false.)
case(TWOMO_INCORE)

   call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
               Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)
end select

print*, 'ABPlus',norm2(ABPlus)
print*, 'ABMin',norm2(ABMin)

! dump matrices for iterative C-ERPA
open(newunit=iunit,file=abfile,form='unformatted')
write(iunit) ABPlus
write(iunit) ABMin
close(iunit)

! UNCOUPLED
select case(Mon%TwoMoInt)
case(TWOMO_FOFO)
   ! are prop-files necessary?
   call Y01CAS_FOFO(Mon%Occ,URe,XOne,ABPlus,ABMin,ECASSCF, &
          'DUMMY','DUMMY',     &
          'DUMMY',xy0file,     &
          Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
          NBas,Mon%NDimX,NInte1,Mon%NoSt,twofile,twojfile,twokfile,&
          Flags%IFlag0,Flags%ICholesky)

   nblk = 1 + NBas - Mon%NAct
   allocate(A0Block(nblk))
   call AC0BLOCK(Mon%Occ,URe,XOne, &
        Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
        NBas,Mon%NDimX,NInte1,twojfile,twokfile, &
        A0BlockIV,A0Block,nblk,abpm0file,1)
end select

deallocate(TwoMO)
deallocate(ABMin,ABPlus)

end subroutine calc_ab_cas

subroutine calc_resp_extrapolate(Mon,Flags,NBasis)
! calculate response for extrapolated SAPT formulas
! i.e., Cubic=.true.
! see Supporting Information in
! 10.1021/acs.jctc.1c00344 

implicit none

type(SystemBlock)            :: Mon
type(FlagsData)              :: Flags
integer,intent(in)           :: NBasis

type(FileNames)              :: FNam
integer                      :: NInte1,NInte2
double precision             :: ACAlpha0,ACAlpha1,ACAlpha2
double precision,allocatable :: XOne(:),TwoNO(:)
double precision,allocatable :: work(:)

NInte1 = NBasis*(NBasis+1)/2
NInte2 = 1
if(Mon%TwoMoInt==TWOMO_INCORE) NInte2 = NInte1*(NInte1+1)/2

call set_filenames(Mon,FNam)

allocate(XOne(NInte1),TwoNO(NInte2))
allocate(work(NBasis*NBasis))

! read 1-el hamiltonian
call get_one_mat('H',work,Mon%Monomer,NBasis)
call sq_to_triang(work,XOne,NBasis)
call tran_matTr(XOne,Mon%CMO,Mon%CMO,NBasis,.true.)

! INCORE: load 2-el integrals
if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoNO,NBasis,NInte2)

print*, 'Monomer',Mon%Monomer
call calc_cas_resp(Mon%ACAlpha0,XOne,TwoNO,Mon,FNam%propfile0,NInte1,NInte2,NBasis,Flags%ICholesky)
call calc_cas_resp(Mon%ACAlpha1,XOne,TwoNO,Mon,FNam%propfile1,NInte1,NInte2,NBasis,Flags%ICholesky)
call calc_cas_resp(Mon%ACAlpha2,XOne,TwoNO,Mon,FNam%propfile2,NInte1,NInte2,NBasis,Flags%ICholesky)

deallocate(work)
deallocate(TwoNO,XOne)

end subroutine calc_resp_extrapolate

subroutine calc_cas_resp(ACAlpha,XOne,TwoMO,Mon,propfile,NInte1,NInte2,NBas,ICholesky)
!
! Two steps: 1) calculate A+B, A-B for alpha=ACAlpha
!            2) solve symmetric ERPA eigenproblem
!
implicit none

type(SystemBlock)            :: Mon
type(FileNames)              :: FNam
integer, intent(in)          :: NInte1,NInte2,NBas
integer, intent(in)          :: ICholesky
double precision, intent(in) :: ACAlpha
double precision, intent(in) :: XOne(NInte1),TwoMO(NInte2)
character(*)                 :: propfile

integer                      :: i
integer                      :: iunit
double precision             :: ECASSCF

double precision, allocatable :: URe(:,:)
double precision, allocatable :: ABPlus(:), ABMin(:)
double precision, allocatable :: EigX(:),EigY(:),Eig(:)

call set_filenames(Mon,FNam)

allocate(URe(NBas,NBas))

URe = 0d0
do i=1,NBas
   URe(i,i) = 1d0
enddo

allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2))
allocate(EigX(Mon%NDimX**2),EigY(Mon%NDimX**2),Eig(Mon%NDimX))

ECASSCF = 0

select case(Mon%TwoMoInt)
case(TWOMO_FOFO)

   call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
               NInte1,FNam%twojfile,FNam%twokfile,ICholesky,ACAlpha,.false.)
case(TWOMO_FFFF)

   call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                NInte1,FNam%twofile,ACAlpha,.false.)
case(TWOMO_INCORE)

   call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
               Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)

   call EKT(URe,Mon%Occ,XOne,TwoMO,NBas,NInte1,NInte2)

end select

if(Mon%InSt(1,1) > 0) then
   write(lout,'(/1x,a,i2,a,i1)') &
   'Calculating reponse for state:',Mon%InSt(1,1),'.',Mon%InSt(2,1)
else
   write(lout,'(1x,a)') 'Calculating reponse for the ground state'
endif

write(lout,'(1x,a/)') 'Solving symmetric eigenvalue equation...'

!if(iter==1) then
!   ! dump matrices for iterative C-ERPA
!   open(newunit=iunit,file=abfile,form='unformatted')
!   write(iunit) ABPlus
!   write(iunit) ABMin
!   close(iunit)
!endif

call ERPASYMMXY(EigY,EigX,Eig,ABPlus,ABMin,&
                Mon%Occ,Mon%IndN,Mon%NDimX,NBas)

! dump response vecs
call writerespXY(EigX,EigY,Eig,propfile)

deallocate(Eig,EigY,EigX)
deallocate(ABMin,ABPlus)
deallocate(URe)

end subroutine calc_cas_resp

subroutine calc_resp_pino(M,MO,Flags,NBas)
! obtain FCI response based on PINO functional
! also: solve full TD-APSG equations
! using APSG_NEST and PINOVECREDXY
implicit none

type(SystemBlock)  :: M
type(FlagsData)    :: Flags
integer,intent(in) :: NBas
double precision   :: MO(NBas*NBas)

integer          :: iPINO
integer          :: i,j,ii,ip,iq,pq
integer          :: NInte1,NInte2,DimEx
integer          :: INegExcit
double precision :: URe(NBas,NBas)
character(:),allocatable      :: onefile,twofile,rdmfile,propfile
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:)
double precision, allocatable :: ABPlus(:),ABMin(:),DMAT(:),DMATK(:),&
                                 CMAT(:),EMAT(:),EMATM(:),EigVecR(:),Eig(:)
! checks
 if(M%TwoMoInt/=TWOMO_INCORE) then
   write(LOUT,'(/,1x,a)') 'ERROR! PINO possible only with TwoMoInt==INCORE!'
   stop
 endif

! set PINO variant
 if(Flags%ICASSCF==1.and.Flags%ISHF==1.and.M%NELE==1.and.Flags%SaptLevel/=10) then
    ! 2-electron FCI
    iPINO = 0
 elseif(Flags%ICASSCF==1.and.M%NELE==1.and.Flags%SaptLevel==10) then
    ! 2-electron e2dispCAS
    iPINO = 1
 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Flags%SaptLevel/=10) then
    ! CAS/LinearResponse
    iPINO = 2
 elseif(Flags%ICASSCF==0.and.Flags%SaptLevel/=10) then
    ! TEST IS MISSING!
    ! GVB
    iPINO = 3
 endif

! set filenames
 if(M%Monomer==1) then
    onefile    = 'ONEEL_A'
    twofile    = 'TWOMOAA'
    rdmfile    = 'rdm2_A.dat'
    propfile   = 'PROP_A'
 elseif(M%Monomer==2) then
    onefile    = 'ONEEL_B'
    twofile    = 'TWOMOBB'
    rdmfile    = 'rdm2_B.dat'
    propfile   = 'PROP_B'
 endif

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

 URe = 0d0
 do i=1,NBas
    URe(i,i) = 1d0
 enddo

 allocate(XOne(NInte1),TwoMO(NInte2))

 ! get integrals
 call get_1el_h_mo(XOne,MO,NBas,onefile)
 call LoadSaptTwoNO(M%Monomer,TwoMO,NBas,NInte2)

 ! Calculate FCI for 2-el systems
 !if(Flags%ISHF==1.and.M%NELE==1) then
 if(iPINO==0.or.iPINO==1) then
    ! FCI
    call calc_fci(NBas,NInte1,NInte2,XOne,URe,TwoMO,M,iPINO)
!   loading 2RDM-incore is not necessary
!   call read2rdm(M,NBas)

    if(iPINO==0) call init_pino(NBas,M,Flags%ICASSCF)

 elseif(iPINO==2) then

    ! CAS
    ! loading 2RDM-incore is not necessary
    !call read2rdm(M,NBas)
    call system('cp '//rdmfile// ' rdm2.dat')

 endif

 M%NDimN=0
 do i=1,NBas
    if(M%Occ(i).gt.0d0) M%NDimN=M%NDimN+1
 enddo
 write(LOUT,'(1x,a,i4)') 'The NDimN dimension: ',M%NDimN

 ! prepare dimensions
 DimEx = M%NDimX+M%NDimN
 M%DimEx = DimEx

 !if((iPINO==2).or.(iPINO==0)) then
 if(iPINO==2) then
    ! CAS/LR -- full linear response for CAS wfn
    allocate(ABPlus(M%NDim**2),ABMin(M%NDim**2), &
             CMAT(M%NDim**2),EMAT(NBas**2),EMATM(NBas**2),&
             DMAT(M%NDim*NBas),DMATK(M%NDim*NBas),&
             EigVecR(2*(M%NDimX+M%NDimN)*2*(M%NDimX+M%NDimN)),&
             Eig(2*(M%NDimX+M%NDimN)))

    CMAT=0
    call APSG_NEST(ABPlus,ABMin,CMAT,EMAT,EMATM,DMAT,DMATK,&
                   URe,M%Occ,XOne,TwoMO,&
                   NBas,M%NDim,NInte1,NInte2,M%NGem,Flags%ISERPA)

    ! begin DE_CAS experiment
    ! this is an experimental piece of code: 
    ! the idea was to correct DMAT and EMAT in ERPA for CAS(n,n)
    ! in a similar way that we include c_p in full TD-GVB 
    ! however, preliminary tests have not given good results
    ! and the idea was abandoned
    !
    !Mon%NDimN = 0
    !  call DE_CAS(DMAT,DMATK,EMAT,EMATM,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
    !       NBas,Mon%NDim,NInte1,NInte2,ACAlpha)

    !  print*, 'DMAT',norm2(DMAT),norm2(DMATK)
    !  print*, 'EMAT',norm2(EMAT),norm2(EMATM)
    !  ij=0
    !  ij1=0
    !  do j = 1,Mon%NDimN
    !     do i = 1,Mon%NDimX
    !        ij  = ij + 1 !(j-1)*Mon%NDimX + i
    !        ij1 = (j-1)*Mon%NDim + Mon%IndXh(i)
    !        DMAT(ij) = DMAT(ij1)
    !        DMATK(ij) = DMATK(ij1)
    !     enddo
    !  enddo
    !  ij=0
    !  ij1=0
    !  do j=1,Mon%NDimN
    !     do i=1,Mon%NDimN
    !        ij  = (j-1)*Mon%NDimN + i
    !        ij1 = (j-1)*Mon%NBasis + i
    !        EMAT(ij) = EMAT(ij1)
    !        EMATM(ij) = EMATM(ij1)
    !     enddo
    !  enddo
    !print*, 'DMAT-R',norm2(DMAT(1:Mon%NDimX*Mon%NDimN)),&
    !     norm2(DMATK(1:Mon%NDimX*Mon%NDimN))
    !print*, 'EMAT-R',norm2(EMAT(1:Mon%NDimN**2)),&
    !     norm2(EMATM(1:Mon%NDimN**2))
    !print*, 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)
    ! end of DE_CAS experiment

    !reduce dimensions
    call reduce_dim('AB',ABPlus,ABMin,M)
    call reduce_dim('D',DMAT,DMATK,M)
    call reduce_dim('E',EMAT,EMATM,M)

    deallocate(CMAT)

 ! GVB/TD-APSG
 elseif(iPINO==3) then

    allocate(ABPlus(M%NDim**2),ABMin(M%NDim**2),  &
             CMAT(M%NDim**2),EMAT(NBas**2),EMATM(NBas**2), &
             DMAT(M%NDim*NBas),DMATK(M%NDim*NBas), &
             EigVecR(2*(M%NDimX+M%NDimN)*2*(M%NDimX+M%NDimN)),&
             Eig(2*(M%NDimX+M%NDimN)))

    CMAT = 0
    call APSG_NEST(ABPlus,ABMin,CMAT,EMAT,EMATM,DMAT,DMATK,&
         URe,M%Occ,XOne,TwoMO,&
         NBas,M%NDim,NInte1,NInte2,M%NGem,Flags%ISERPA)

    !reduce dimensions
    call reduce_dim('AB',ABPlus,ABMin,M)
    call reduce_dim('D',DMAT,DMATK,M)
    call reduce_dim('E',EMAT,EMATM,M)

    deallocate(CMAT)
 endif

 if(iPINO/=0) then
    ! we need eigenprolem (PINOVECREDXY)
    ! for CAS/LR and TD-GVB

    EigVecR = 0
    Eig = 0

    allocate(M%EigY((M%NDimX+M%NDimN)**2),M%EigX((M%NDimX+M%NDimN)**2),&
             M%Eig(M%NDimX+M%NDimN))

    !print*, 'dimernsions',M%NDimX,M%NDimN
    !print*,'ABPlus',norm2(ABPLus)
    !print*,'ABMin ',norm2(ABMin)
    !print*,'DMAT ',norm2(DMAT)
    !print*,'EMAT ',norm2(EMAT)

   !do i=1,NBas
   !   print*, i,M%Occ(i),M%CICoef(i)
   !enddo

    call PINOVECREDXY(M%EigY,M%EigX,M%Eig,INegExcit,&
                      ABPlus,ABMin,DMAT,DMATK,EMAT,EMATM,M%IndN,&
                      NBas,M%NDimX,M%NDimN)

    print*, 'PINOVECREDXY:'
    print*,'EigY',norm2(M%EigY(1:(M%NDimX+M%NDimN)**2))
    print*,'EigX',norm2(M%EigX(1:(M%NDimX+M%NDimN)**2))

    Eig = M%Eig

    allocate(M%IndNx(2,DimEx))

    do pq=1,DimEx
       if(pq<=M%NDimX) then
          M%IndNx(1,pq) = M%IndN(1,pq)
          M%IndNx(2,pq) = M%IndN(2,pq)
       elseif(pq>M%NDimX) then
          M%IndNx(1,pq) = pq - M%NDimX
          M%IndNx(2,pq) = pq - M%NDimX
       endif
    enddo

    ! prepare EigVecR for E2disp
    EigVecR = 0
    do j=1,DimEx
       if(j<=M%NDimX) then

          ip = M%IndN(1,j)
          iq = M%IndN(2,j)

          do i=1,DimEx
             EigVecR((i-1)*DimEx+j) = (M%CICoef(ip)-M%CICoef(iq))* &
                                      (M%EigY((i-1)*DimEx+j)-M%EigX((i-1)*DimEx+j))
          enddo

       elseif(j>M%NDimX) then

          ii = j - M%NDimX
          do i=1,DimEx
             EigVecR((i-1)*DimEx+j) = M%EigY((i-1)*DimEx+j)/M%CICoef(ii)
          enddo

       endif
    enddo

   ! this has not been used in a long time...
   !ETot = 0
   !call EnePINO(ETot,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne, &
   !     Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NDimN)

    ! dump response
    call writeresp(EigVecR,Eig,propfile)

    deallocate(EMAT,EMATM,DMAT,DMATK)
    deallocate(ABMin,ABPlus,EigVecR,Eig)

 endif

 deallocate(TwoMO,XOne)

end subroutine calc_resp_pino

subroutine calc_resp_unc(Mon,MO,Flags,NBas)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

integer,intent(in)            :: NBas
double precision,intent(in)   :: MO(:)

integer          :: i,ione
integer          :: NSq,NInte1,NInte2
double precision :: ETot
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:)
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:),      &
                                 EigY(:), EigY1(:), Eig(:), Eig1(:)
character(8)                  :: label
character(:),allocatable      :: twojfile,twokfile
character(:),allocatable      :: onefile,twofile,propfile0,propfile1,rdmfile
character(:),allocatable      :: y01file,xy0file

! set filenames
if(Mon%Monomer==1) then
   onefile   = 'ONEEL_A'
   twofile   = 'TWOMOAA'
   twojfile  = 'FFOOAA'
   twokfile  = 'FOFOAA'
   propfile0 = 'PROP_A0'
   propfile1 = 'PROP_A1'
   rdmfile   = 'rdm2_A.dat'
   y01file   = 'Y01_A'
   xy0file   = 'XY0_A'
elseif(Mon%Monomer==2) then
   onefile   = 'ONEEL_B'
   twofile   = 'TWOMOBB'
   twojfile  = 'FFOOBB'
   twokfile  = 'FOFOBB'
   propfile0 = 'PROP_B0'
   propfile1 = 'PROP_B1'
   rdmfile   = 'rdm2_B.dat'
   y01file   = 'Y01_B'
   xy0file   = 'XY0_B'
endif

! set dimensions
NSq = NBas**2
NInte1 = NBas*(NBas+1)/2
NInte2 = NInte1*(NInte1+1)/2

allocate(work1(NSq),work2(NSq),XOne(NInte1),URe(NBas,NBas))
if(Mon%TwoMoInt==TWOMO_INCORE) then
   allocate(TwoMO(NInte2))
endif

URe = 0d0
do i=1,NBas
   URe(i,i) = 1d0
enddo

! read 1-el
open(newunit=ione,file=onefile,access='sequential',&
     form='unformatted',status='old')

read(ione)
read(ione)
read(ione) label, work1
if(label=='ONEHAMIL') then
   call tran_oneint(work1,MO,MO,work2,NBas)
   call sq_to_triang(work1,XOne,NBas)
else
   write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
   stop
endif

if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBas,NInte2)

! GVB
if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

   allocate(EigY(Mon%NDimX**2), &
        Eig(Mon%NDimX),Eig1(Mon%NDimX))

   EigY  = 0
   Eig   = 0
   Eig1  = 0

   call Y01GVB(TwoMO,Mon%Occ,URe,XOne, &
        EigY,Eig,Eig1, &
        Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2)

   ! dump uncoupled response
   call writeResp(EigY,Eig,propfile0)
   if(Flags%IFlag0==0) then
      call writeEval(Eig1,propfile1)
   endif

   deallocate(Eig1,Eig,EigY)

! CAS
elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then

  allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2))

  select case(Mon%TwoMoInt)
  case(TWOMO_FOFO)
     call Y01CAS_FOFO(Mon%Occ,URe,XOne,ABPlus,ABMin,ETot, &
            propfile0,propfile1, &
            y01file,xy0file,     &
            Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
            NBas,Mon%NDim,NInte1,Mon%NoSt,twofile,twojfile,twokfile, &
            Flags%IFlag0,Flags%ICholesky)
  case(TWOMO_FFFF)
     call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
            propfile0,propfile1, &
            y01file,&
            Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
            NBas,Mon%NDim,NInte1,Mon%NoSt,twofile,Flags%IFlag0)
  case(TWOMO_INCORE)
     allocate(EigY(Mon%NDimX**2),EigY1(Mon%NDimX**2), &
              Eig(Mon%NDimX),Eig1(Mon%NDimX))

     EigY  = 0
     EigY1 = 0
     Eig   = 0
     Eig1  = 0

     call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
          EigY,EigY1,Eig,Eig1, &
          Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0,xy0file)

     ! dump uncoupled response
     call writeResp(EigY,Eig,propfile0)
     if(Flags%IFlag0==0) then
        call writeResp(EigY1,Eig1,propfile1)
     endif

     deallocate(Eig1,Eig,EigY1,EigY)
  end select

  deallocate(ABMin,ABPlus)

  Mon%ECASSCF = ETot+Mon%PotNuc
  write(LOUT,'(1x,a,5x,f15.8)') "CASSCF Energy           ",Mon%ECASSCF

endif


close(ione)
if(Mon%TwoMoInt==1) deallocate(TwoMO)
deallocate(URe,XOne,work1,work2)

end subroutine calc_resp_unc

subroutine calc_resp_dft(Mon,MO,Flags,NBas)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

double precision :: MO(:)
integer :: NBas
integer :: NSq,NInte1,NInte2,NGrid,NDimKer
integer :: NSymNO(NBas),MultpC(15,15)
double precision, allocatable :: work1(:),work2(:)
double precision, allocatable :: XOne(:), &
                                 TwoMO(:), TwoElErf(:), &
                                 WGrid(:),XKer(:),OrbGrid(:),&
                                 OrbXGrid(:),OrbYGrid(:),OrbZGrid(:),&
                                 SRKer(:),SRKerW(:)
double precision, allocatable :: ABPlus(:),ABMin(:),URe(:,:),VSR(:), &
                                 EigY0(:),EigY1(:),Eig0(:),Eig1(:), &
                                 EigVecR(:), Eig(:)
integer          :: i,j,ip,iq,ii,ione
double precision :: ACAlpha,Omega,EnSR,EnHSR,ECorr,ECASSCF,XVSR
double precision :: Tcpu,Twall
character(8)     :: label
character(:),allocatable :: onefile,aoerfile,twofile,twoerffile,&
                            twojfile,twokfile,twojerf,twokerf,  &
                            propfile,propfile0,propfile1,xy0file,rdmfile
double precision,parameter :: One = 1d0, Half = 0.5d0
logical :: doRSH

! temporary RSH solution
doRSH = .false.
if(Flags%IFunSR==1.or.Flags%IFunSR==2) doRSH = .true.

! set filenames
if(Mon%Monomer==1) then
   onefile    = 'ONEEL_A'
   twofile    = 'TWOMOAA'
   twojfile   = 'FFOOAA'
   twokfile   = 'FOFOAA'
   twoerffile = 'MO2ERFAA'
   twojerf    = 'FFOOERFAA'
   twokerf    = 'FOFOERFAA'
   propfile   = 'PROP_A'
   propfile0  = 'PROP_A0'
   propfile1  = 'PROP_A1'
   xy0file    = 'XY0_A'
   rdmfile    = 'rdm2_A.dat'
   aoerfile   = 'AOERFSORT'
elseif(Mon%Monomer==2) then
   onefile    = 'ONEEL_B'
   twofile    = 'TWOMOBB'
   twojfile   = 'FFOOBB'
   twokfile   = 'FOFOBB'
   twoerffile = 'MO2ERFBB'
   twojerf    = 'FFOOERFBB'
   twokerf    = 'FOFOERFBB'
   propfile   = 'PROP_B'
   propfile0  = 'PROP_B0'
   propfile1  = 'PROP_B1'
   xy0file    = 'XY0_B'
   rdmfile    = 'rdm2_B.dat'
   if(Mon%SameOm) then
      aoerfile = 'AOERFSORT'
   else
      aoerfile = 'AOERFSORTB'
   endif
endif

! set dimensions
NSq = NBas**2
NInte1 = NBas*(NBas+1)/2
NInte2 = NInte1*(NInte1+1)/2

call create_symmats(Mon,MO,NBas)

allocate(work1(NSq),work2(NSq),XOne(NInte1),URe(NBas,NBas),VSR(NInte1))
if(Mon%TwoMoInt==1) then
   allocate(TwoMO(NInte2))
   allocate(TwoElErf(NInte2))
   !if(doRSH) allocate(TwoElErf(NInte2))
endif

URe = 0d0
do i=1,NBas
   URe(i,i) = 1d0
enddo

! read 1-el
! add Coulomb (RSH: sr Coulomb)
call triang_to_sq(Mon%VCoul,work2,NBas)
open(newunit=ione,file=onefile,access='sequential',&
     form='unformatted',status='old')

read(ione)
read(ione)
read(ione) label, work1
if(label=='ONEHAMIL') then
   work1 = work1 + work2
   call tran_oneint(work1,MO,MO,work2,NBas)
   call sq_to_triang(work1,XOne,NBas)
else
   write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
   stop
endif

! load grid
call molprogrid0(NGrid,NBas)
write(LOUT,'()')
write(LOUT,'(1x,a,i8)') "The number of Grid Points =",NGrid

allocate(WGrid(NGrid),SRKer(NGrid),SRKerW(NGrid), &
         OrbGrid(NGrid*NBas), &
         OrbXGrid(NGrid*NBas),OrbYGrid(NGrid*NBas),OrbZGrid(NGrid*NBas))

! load orbgrid and gradients, and wgrid
call transp_mat1dim(MO,work1,NBas)
call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid, &
                               WGrid,work1,NGrid,NBas)

! set/load Omega - range separation parameter
write(LOUT,'(1x,a,f15.8)') "The range-separation parameter =",Mon%Omega

! transform 2-el integrals
select case(Mon%TwoMoInt)
case(TWOMO_INCORE,TWOMO_FFFF)
    ! full
    call tran4_full(NBas,MO,MO,twofile,'AOTWOSORT')
case(TWOMO_FOFO)
   ! transform J and K
   call tran4_gen(NBas,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            NBas,MO,&
            NBas,MO,&
            twojfile,'AOTWOSORT')
   call tran4_gen(NBas,&
            NBas,MO,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            NBas,MO,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            twokfile,'AOTWOSORT')
end select
if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBas,NInte2)

! transform LR integrals
if(doRSH) then
  select case(Mon%TwoMoInt)
  case(TWOMO_INCORE)
     TwoElErf(1:NInte2) = 0
    print*, 'Check: calc_resp_dft:',Mon%SameOm,aoerfile
     call tran4_full(NBas,MO,MO,twoerffile,aoerfile)
  case(TWOMO_FFFF)
     call tran4_full(NBas,MO,MO,twoerffile,aoerfile)
  case(TWOMO_FOFO)
     call tran4_gen(NBas,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             NBas,MO,&
             NBas,MO,&
             twojerf,aoerfile)
     call tran4_gen(NBas,&
             NBas,MO,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             NBas,MO,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             twokerf,aoerfile)
  end select
  if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer+4,TwoElErf,NBas,NInte2)
endif

NSymNO(1:NBas) = 1
call EPotSR(EnSR,EnHSR,VSR,Mon%Occ,URe,MO,.true.,&
           OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,&
!            NSymNO,TwoMO,TwoElErf,&
           NSymNO,Mon%VCoul,&
           Mon%Omega,Flags%IFunSR,&
           NGrid,NInte1,NInte2,NBas)

! MODIFY XOne in PostCAS calculations
if(Mon%PostCAS) then
   write(LOUT,*) 'TEST POSTCAS!'
   XOne = XOne + VSR
endif
write(LOUT,'(1x,a,f15.8)') "SR Energy: ",EnSR

XVSR = 0
do i=1,NBas
   ii = (i*(i+1))/2
   XVSR = XVSR + 2d0*Mon%Occ(i)*VSR(ii)
enddo

ACAlpha=One
! UNCOUPLED-TEST
!ACAlpha=1d-9

allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
         EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

! read 2-RDMs
! HERE-2?
 call system('cp '//rdmfile// ' rdm2.dat')

ECASSCF = 0
! if(doRSH) then

select case(Mon%TwoMoInt)
case(TWOMO_INCORE)
   call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoElErf,Mon%IPair,&
               Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)
   call EKT(URe,Mon%Occ,XOne,TwoElErf,NBas,NInte1,NInte2)
case(TWOMO_FFFF)
   call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
               NInte1,twoerffile,ACAlpha,.false.)
case(TWOMO_FOFO)
   call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
               Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
               NInte1,twojerf,twokerf,Flags%ICholesky,ACAlpha,.false.)
!else
! HERE:: ADD SEPARATE PROCEDURE FOR Kohn-Sham!
!endif
end select
write(LOUT,'(/,1x,a,f16.8,a,1x,f16.8)') 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)

! ADD CONTRIBUTIONS FROM THE (sr)ALDA KERNEL TO AB MATRICES
MultpC(1,1)=1
call GetKerNPT(SRKer,Mon%Occ,URe,OrbGrid,WGrid,NSymNO,MultpC, &
               NBas,NGrid)
call clock('START',Tcpu,Twall)
select case(Mon%TwoMoInt)
case(TWOMO_INCORE)
   call ModABMin(Mon%Occ,SRKer,WGrid,OrbGrid,TwoMO,TwoElErf,ABMin,&
                 Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NInte2,NBas)
case(TWOMO_FFFF)
   call ModABMin_mithap(Mon%Occ,SRKer,WGrid,OrbGrid,ABMin,&
                        Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NBas,&
                        twofile,twoerffile)
case(TWOMO_FOFO)
   call ModABMin_FOFO(Mon%Occ,SRKer,WGrid,OrbGrid,ABMin,&
                      Mon%MultpC,Mon%NSymNO,&
                      Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NBas,&
                      Mon%num0,Mon%num1, &
                      twokfile,twokerf,.false.)
   print*, 'ABMin-MY',norm2(ABMin)
end select
call clock('Mod ABMin',Tcpu,Twall)

!write(LOUT,'(1x,a)') "*** sr-kernel added. ***"
! test true energy
write(LOUT,'(/,1x,a,f15.8)') "Total lrCASSCF+ENuc+srDF Energy", ECASSCF-XVSR+EnSR+Mon%PotNuc

EigVecR = 0
Eig = 0
if(Mon%NoSt==1) then
!  call ERPASYMM1(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)

   allocate(Mon%EigY(Mon%NDimX**2),Mon%EigX(Mon%NDimX**2),&
            Mon%Eig(Mon%NDimX))
   call ERPASYMMXY(Mon%EigY,Mon%EigX,Mon%Eig,ABPlus,ABMin,&
                   Mon%Occ,Mon%IndN,Mon%NDimX,NBas)

   !print*, 'EigY:',norm2(Mon%EigY),norm2(Mon%EigX)

   Eig = Mon%Eig
    do j=1,Mon%NDimX
       do i=1,Mon%NDimX
          ip = Mon%IndN(1,i)
          iq = Mon%IndN(2,i)
          EigVecR((j-1)*Mon%NDimX+i)=(Mon%CICoef(ip)-Mon%CICoef(iq))&
                                    *(Mon%EigY((j-1)*Mon%NDimX+i)-Mon%EigX((j-1)*Mon%NDimX+i))
       enddo
    enddo

 !! TEST COUPLED ANDREAS
 ! print*, 'CPLD-ANDREAS:'
 ! Mon%EigX = EigVecR
 ! Mon%EigY = 0

else
   call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
endif

write(LOUT,'(/," *** LR-CAS-SR-DFT Excitation Energies *** ",/)')
do i=1,10
   write(LOUT,'(i4,4x,e16.6)') i,Eig(i)
enddo

!print*, 'Check! Eig,EigVecR',norm2(Eig),norm2(EigVecR)

 !ECorr=0
 !!call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
 !!     Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
 !  call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Mon%Occ, &
 !                      Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
 !                      Mon%NDimX,NBas,twokfile,Flags%ICholesky)
 !ECorr=Ecorr*0.5d0
 !write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
 !     ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr

! dump response
call writeresp(EigVecR,Eig,propfile)

! uncoupled
!allocate(EigY0(Mon%NDimX**2),Eig0(Mon%NDimX))
allocate(EigY1(Mon%NDimX**2),Eig1(Mon%NDimX))

!EigY0 = 0
!Eig0  = 0
EigY1 = 0
Eig1  = 0

!call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
!       EigY0,EigY1,Eig0,Eig1, &
!       Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
!       NBas,Mon%NDim,NInte1,twoerffile,Flags%IFlag0)

! call Y01CAS(TwoElErf,Mon%Occ,URe,XOne,ABPlus,ABMin, &
!      EigY0,EigY1,Eig0,Eig1, &
!      Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)
! do i=1,NGrid
!    SRKerW(i) = SRKer(i)*WGrid(i)
! enddo
! call Y01CASLR(TwoElErf,Mon%Occ,URe,XOne,ABPlus,ABMin, &
!        EigY0,EigY1,Eig0,Eig1, &
!        Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0, &
!        TwoMO,OrbGrid,SRKerW,NSymNO,MultpC,NGrid)

select case(Mon%TwoMoInt)
case(TWOMO_INCORE)
  write(lout,'(1x,a)') 'Error! Uncoupled not working for INCORE!'
  stop
case(TWOMO_FOFO)

  call Y01CASLR_FOFO(Mon%Occ,URe,XOne,ABPlus,ABMin,&
                 Mon%MultpC,Mon%NSymNO,SRKer,WGrid,OrbGrid,&
                 propfile0,propfile1,xy0file,&
                 Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,&
                 NGrid,Mon%NDimX,NBas,Mon%NDimX,NInte1,Mon%NoSt,&
                 twokfile,Twojerf,twokerf,Flags%ICholesky,1,1)

end select

! dump uncoupled response
!call writeresp(EigY0,Eig0,propfile0)
if(Flags%IFlag0==0) then
   call writeresp(EigY1,Eig1,propfile1)
endif

close(ione)
!deallocate(Eig0,EigY0)
deallocate(Eig1,EigY1)

deallocate(Eig,EigVecR,ABMin,ABPlus)
deallocate(SRKer,SRKerW,OrbZGrid,OrbYGrid,OrbXGrid,OrbGrid,WGrid)
if(Mon%TwoMoInt==1) then
   deallocate(TwoMO)
   ! if(doRSH) deallocate(TwoElErf)
   deallocate(TwoElErf)
endif
deallocate(VSR,URe,XOne,work1,work2)

end subroutine calc_resp_dft

subroutine sapt_ERPA_TRDMs(Flags,Mon,NBasis)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags
integer,intent(in) :: NBasis

type(EBlockData)             :: SBlockIV
type(EBlockData),allocatable :: SBlock(:)
character(:),allocatable :: xy0file
integer :: iblk, nblk,nAA,NU
integer :: i,j, ii ,ipos, ip, iq,ir,is,it, pq
integer :: dimO,exc
double precision, allocatable :: EigX0(:,:),EigY0(:,:)
double precision :: TRDM1(NBasis,NBasis),TRDM1CAS(NBasis,NBasis),XCAS(NBAsis,NBasis),YCAS(NBasis,NBASIs),C(NBasis)
double precision,allocatable :: RDM2(:,:,:,:),montrdm24(:,:,:,:),TRDM1TEST(:,:),TRDM1fromTilde(:,:),RDM1TEST(:,:)
integer :: iref,iref2,iexcited,NumX, IStERPA
double precision :: SSgn,SumNU
double precision, allocatable :: Eig(:),EigY(:,:),EigX(:,:),iaddr(:)
double precision:: OvMax, OvXMax, OvYMax, SumERY, SumERX, SumCAY, SumCAX,SOvY,SOvX, EigOv
double precision:: YER,YCA,XER,XCA,valX,valY
double precision :: CICoef(NBasis)
double precision :: TESTdipMom,TRUEdipMom,N
double precision :: SumTrdmTrdm,TRDM1CASnorm,TRDM1norm
! CHANGE IN FUTURE
CICoef = mon%CICoef
allocate(mon%trdm24(NBasis,NBasis,NBasis,NBasis))


if(Mon%Monomer == 1) then
   xy0file = "XY0_A"
elseif(Mon%Monomer == 2) then
   xy0file = "XY0_B"
endif

dimO  = mon%num0+mon%num1
call read_SBlock(SBlock,SBlockIV,nblk,xy0file)
!
! First, we want to make XCas a YCas vectors
!
! We transform TRDMs to NO basis of IREF1 state

! nAA   = 0
! do i=1,NDimX
!    ip = IndN(1,i)
!    iq = IndN(2,i)
!    if(IGem(ip)==2.and.IGem(iq)==2) then
!       nAA = nAA + 1
!    endif
! enddo

iref   =  mon%IREF1 !1
iref2  =  mon%IREF2 !4
iexcited = ((iref2-2)*(iref2-1))/2+iref

call tran2MO(mon%trdm1(iexcited,:,:),mon%CMONO(iref,:,:),mon%CMONO(iref,:,:), &
TRDM1CAS(:,:),NBasis)
!
! now we make X_Cas and Y_Cas from TRDM1CAS
!

SumNU=0d0
XCAS=0d0
YCAS=0d0
! do ip=1,NBasis
!    do iq=1,ip-1
!       if(mon%Occ(ip)-mon%Occ(iq).ne.0d0) then
!          YCAS(ip,iq)=TRDM1CAS(iq,ip)/(mon%Occ(ip)-mon%Occ(iq))
!          XCAS(ip,iq)=-TRDM1CAS(ip,iq)/(mon%Occ(ip)-mon%Occ(iq))
!       ! if(mon%Occ(ip)-mon%OCC(iq).ne.0d0) then
!       !     YCAS(ip,iq)=TRDM1CAS(iq,ip)/(mon%Occ(ip)-mon%OCC(iq))
!       !     XCAS(ip,iq)=-TRDM1CAS(ip,iq)/(mon%Occ(ip)-mon%OCC(iq))
!       endif
!    enddo
! enddo
block
! assemble X_SA-CAS / Y_SA-CAS like in Kasia's code
! works with skipped XTilde-->X transformation in dump_EBlock
double precision :: CICoef(NBasis)
CICoef = mon%CICoef
do ip=1,NBasis
  do iq=1,ip-1
     If(CICoef(ip)+CICoef(iq).Ne.0d0) then
       YCAS(ip,iq) = (trdm1cas(ip,iq)+trdm1cas(iq,ip))/(CICoef(ip)+CICoef(iq))
     endif
     If(CICoef(ip)-CICoef(iq).Ne.0d0) then
       XCAS(ip,iq) = (trdm1cas(iq,ip)-trdm1cas(ip,iq))/(CICoef(ip)-CICoef(iq))
     endif
  enddo
enddo
end block



do ip=1,NBasis
   do iq=1,ip-1
      SumNU=SumNU+YCAS(ip,iq)*XCAS(ip,iq)
   enddo
enddo

SSgn=1d0
write(6,'(1x,a,i2)') 'Monomer',mon%monomer
write(6,'(X,"SumNu Y*X before normalization",E15.6)') SumNU
If(SumNU.Lt.0d0) SSgn=-1d0
If(Abs(SumNu).Gt.1.D-8) Then
SumNU=2d0/Sqrt(Abs(SumNU))
Else
SumNU=0d0
EndIf
do ip=1,NBasis
   do iq=1,ip-1
      YCAS(ip,iq)=YCAS(ip,iq)*SumNU
      XCAS(ip,iq)=SSgn*XCAS(ip,iq)*SumNU
      If(Abs(YCAS(ip,iq))+Abs(XCAS(ip,iq)).Gt.1.d-6) Then 
         Write(6,'(X,"NORMALIZED Y_SA-CAS X_SA-CAS",2I3,2E15.6)') ip,iq,YCAS(ip,iq),XCAS(ip,iq)
      EndIf
   EndDo
EndDo
!





!
!
!tu skończyliśmy
!
!
!
nAA=SBlock(1)%n
allocate(Eig(nAA),EigY(nAA,nAA),EigX(nAA,nAA),iaddr(nAA))
allocate(TRDM1TEST(NBasis,NBasis))
Eig = 0
EigY = 0
! nAA to liczba aktywnych
   associate(B => Sblock(1)) 
      do i=1,B%n
         print*,"i,B%vec(i),B%matY(1:B%n,i)"
         print*,i,B%vec(i),B%matY(1:B%n,i)

! ZMIANA 11.12.2023
         !   If(B%vec(i).lt.0) Then
         !       Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)') B%vec(i)
         !       B%vec(i)=0
         !       B%matY(1:B%n,i)=0
         !       B%matX(1:B%n,i)=0
         !   EndIf
! KONIEC ZMIANY
           iaddr(i)=B%pos(i)
          ! EigY(i,B%l1:B%l2) = B%matY(i,1:B%n)
          ! EigX(i,B%l1:B%l2) = B%matX(i,1:B%n)
          !!! UPDATE 21.04 !!!
           ip=mon%IndN(1,iaddr(i))
           iq=mon%IndN(2,iaddr(i))
           EigY(i,B%l1:B%l2) = (B%matY(i,1:B%n) - B%matX(i,1:B%n))*(CICoef(ip)-CICoef(iq))
           EigX(i,B%l1:B%l2) = (B%matY(i,1:B%n) + B%matX(i,1:B%n))*(CICoef(ip)+CICoef(iq))
     enddo
     Eig(B%l1:B%l2) = B%vec(1:B%n)
   end associate

!   if(NoSt>1) then
!     Call SortEigXY(1,Eig,EigY,EigX,nAA)
!   endif

   Write(6,'(X,"First 4 Active ERPA Eigenvalues and Eigenvecs")')
   do i=1,4
   !do i=1,nAA
      Write(6,'(X,I4,E15.6)') i,eig(i)
      do iq=1,nAA
         if(abs(eigy(iq,i))+abs(eigx(iq,i)).gt.1.d-7) Write(6,'(X,"Y_ERPA, X_ERPA",2I3,2E15.6)') &
         mon%IndN(1,iaddr(iq)),mon%IndN(2,iaddr(iq)),eigy(iq,i),eigx(iq,i)
      enddo
   enddo
      Write(6,*)
      OvMax=0
      OvXMax=0
      OvYMax=0

      Do NU=1,nAA

      SumERY=0
      SumERX=0
      SumCAY=0
      SumCAX=0

      Do I=1,nAA
      ip=mon%IndN(1,iaddr(i))
      iq=mon%IndN(2,iaddr(i))
      YER=eigy(i,NU)
      YCA=YCAS(ip,iq)
      XER=eigx(i,NU)
      XCA=XCAS(ip,iq)

      SumERX=SumERX+XER*XER
      SumCAX=SumCAX+XCA*XCA
      SumERY=SumERY+YER*YER
      SumCAY=SumCAY+YCA*YCA
      EndDo

      If(Abs(SumERX).Gt.1.D-12) SumERX=1.D0/Sqrt(Abs(SumERX))
      If(Abs(SumCAX).Gt.1.D-12) SumCAX=1.D0/Sqrt(Abs(SumCAX))
      If(Abs(SumERY).Gt.1.D-12) SumERY=1.D0/Sqrt(Abs(SumERY))
      If(Abs(SumCAY).Gt.1.D-12) SumCAY=1.D0/Sqrt(Abs(SumCAY))

      SOvY=0d0
      SOvX=0d0

      Do I=1,nAA
      ip=mon%IndN(1,iaddr(i))
      iq=mon%IndN(2,iaddr(i))
      YER=eigy(i,NU)
      YCA=YCAS(ip,iq)
      XER=eigx(i,NU)
      XCA=XCAS(ip,iq)
      SOvY=SOvY+YER*YCA*SumERY*SumCAY
      SOvX=SOvX+XER*XCA*SumERX*SumCAX
      EndDo

      SOvY=Abs(SOvY)
      SOvX=Abs(SOvX)
      Write(6,'(X,"SA-CAS-Overlap for ERPA vector",I2,3F15.8)') NU,SOvY,SOvX,(SOvY+SovX)/2.D0

      If((SOvY+SovX)/2.D0.Gt.OvMax) Then
      OvMax=(SOvY+SovX)/2.D0
      OvXMax=SOvX
      OvYMax=SOvY
      NUMx=NU
      EndIf

      EndDo

      If((OvMax.Ge.0.5).Or. &
       ( (OvMax.Lt.0.5).And.(OvXMax.Gt.0.5.Or.OvYMax.Gt.0.5) )) Then

         IStERPA=NUMx
         EigOv=Eig(IStERPA)

         Write(6,'(/,X,"ERPA vector best overlapping with SA-CAS is vector no", &
         I2," of excit energy:",F15.8)') IStERPA,Eig(IStERPA)
         NU=IStERPA

         Do I=1,nAA
         ip=mon%IndN(1,iaddr(i))
         iq=mon%IndN(2,iaddr(i))
         If(Abs(eigy(i,NU))+Abs(eigx(i,NU)).Gt.1.d-7) &
         Write(6,'(X,"Y_ERPA, X_ERPA        ",2I3,2E15.6)')ip,iq,eigy(i,NU),eigx(i,NU)
         EndDo

         ! Write(6,'(X,"best-matching SA-CAS Excitation from ", &
         !       I4," to ",I4,23X,F15.8)') NoSt,IH0St,EExcit(IH0St)

         Do I=1,nAA
         ip=mon%IndN(1,iaddr(i))
         iq=mon%IndN(2,iaddr(i))
         If(Abs(YCAS(ip,iq))+Abs(XCAS(ip,iq)).Gt.1.d-7) &
         Write(6,'(X,"Y_SA-CAS, X_SA-CAS    ",2I3,2E15.6)')ip,iq, &
         YCAS(ip,iq),XCAS(ip,iq)
         EndDo

      Else

         IStERPA=0
         Write(6,'(/,X, "No ERPA vector overlaps with SA-CAS State No",I2)')mon%Iref2

      EndIf
deallocate(Eig,EigY,EigX,iaddr)









! If(IA.Ne.IB.And.IPair(IA,IB).Eq.1) Then
! c      If(IA.Ne.IB) Then
! If(C(IA)+C(IB).Ne.Zero) YCAS(IKet,IAB)=
! $ (trdm(IA,IB)+trdm(IB,IA))/(C(IA)+C(IB))
! If(C(IA)-C(IB).Ne.Zero) XCAS(IKet,IAB)=
! $ (trdm(IB,IA)-trdm(IA,IB))/(C(IA)-C(IB))
! EndIf
! C
! If(IPr.Eq.1.And.Abs(trdm(IA,IB))+Abs(trdm(IB,IA)).gt.1.d-6) Then
! Write(6,'(X,"TRDM_ab 1-TRDM_ba",2I3,2E15.6)')IA,IB,
! $ trdm(IA,IB),trdm(IB,IA)
! c      Write(6,'(X,"YCAS, XCAS    ",6X,2E15.6)')
! c     $ YCAS(IKet,IAB),XCAS(IKet,IAB)
! EndIf
! C
! SumNU=SumNU+YCAS(IKet,IAB)*XCAS(IKet,IAB)
! C
! EndDo
! EndDo
! C
! SSgn=One
! If(IPr.Eq.1) 
! $ Write(6,'(X,"SumNu Y*X before normalization",E15.6)')SumNU
! If(SumNU.Lt.Zero) SSgn=-One
! If(Abs(SumNu).Gt.1.D-8) Then
! SumNU=One/Sqrt(Two*Abs(SumNU))
! Else
! SumNU=Zero
! EndIf
! C
! IAB=0
! Do IA=1,NBasis
! Do IB=1,IA
! IAB=IAB+1
! YCAS(IKet,IAB)=YCAS(IKet,IAB)
! $ *SumNU
! XCAS(IKet,IAB)=SSgn*XCAS(IKet,IAB)
! $ *SumNU
! If(IPr.Eq.1.And.Abs(YCAS(IKet,IAB))+Abs(XCAS(IKet,IAB)).Gt.1.d-6)
! $  Then 
!  Write(6,'(X,"NORMALIZED Y_SA-CAS X_SA-CAS",2I3,2E15.6)')
! $ IA,IB,YCAS(IKet,IAB),XCAS(IKet,IAB)
! EndIf
! EndDo
! EndDo





allocate(RDM2(dimO,dimO,dimO,dimO))
RDM2(:,:,:,:) = mon%rdm24(1,:,:,:,:)
allocate(montrdm24(NBasis,NBasis,NBasis,NBasis))
!allocate(mon%trdm24(4,NBasis,NBasis,NBasis,NBasis))
montrdm24(:,:,:,:)=0d0
! exc = NU !!! tak trzeba będzie tu wrocic
TRDM1 = 0d0
allocate(EigY0(mon%NDimX,mon%NDimX),EigX0(mon%NDimX,mon%NDimX))
!block
!! assemble X_SA-CAS / Y_SA-CAS like in Kasia's code
!! works with skipped XTilde-->X transformation in dump_EBlock
!
!
!do ip=1,NBasis
!   do iq=1,ip-1
!      If(CICoef(ip)+CICoef(iq).Ne.0d0) then
!        YCAS(ip,iq) = (trdm1cas(ip,iq)+trdm1cas(iq,ip))/(CICoef(ip)+CICoef(iq))
!      endif
!      If(CICoef(ip)-CICoef(iq).Ne.0d0) then
!        XCAS(ip,iq) = (trdm1cas(iq,ip)-trdm1cas(ip,iq))/(CICoef(ip)-CICoef(iq))
!      endif
!   enddo
!enddo
!end block
allocate(RDM1TEST(diMO,dimO))
RDM1TEST = 0d0
do ip=1,dimO
   do iq=1,dimO
      do ir=1,dimO
         RDM1TEST(ip,iq)=RDM1TEST(ip,iq)+RDM2(ip,iq,ir,ir)
      enddo
   enddo
enddo
N=0d0
do ip=1,dimO
   N=N+RDM1TEST(ip,ip)
enddo
print *,"'Normalizacja",N,2*mon%NELE-1,mon%ZNucl-1
! Zmiana 11.30.2023
RDM1TEST(:,:)=RDM1TEST(:,:)/(mon%ZNucl-1)
!RDM1TEST(:,:)=RDM1TEST(:,:)/(2*mon%NELE)
!KONIEC ZMIANY
!RDM1TEST(:,:)=RDM1TEST(:,:)/(2*mon%NELE-1)
!USUNIĘTY FRAGMENT
print*,"1RDM z 2RDM dla monomeru", Mon%Monomer
do i=1,dimO
   write(lout,'(*(f12.6))') (RDM1TEST(i,j),j=1,dimO)
enddo
print*,"1rdm",mon%rdm1(1,:)

   NU=IStERPA
   ! ZMIANA 11.12.2023
   ! NIE NA STAŁE
 !  NU = 1
   ! KONIEC ZMIANY

   print*,"We take vector number: ", NU
   EigY0 = 0d0
   EigX0 = 0d0
   ! unpack (1)
   do iblk=1,nblk
      associate(B => Sblock(iblk))
         if(B%l1 /= 1) cycle 
         do i=1,B%n
         ipos = B%pos(i)
   !        EigY0(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
   !        EigX0(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
         ip = mon%IndN(1,ipos)
         iq = mon%IndN(2,ipos)
         print *,ip,iq
               TRDM1(ip,iq)=TRDM1(ip,iq)-(mon%rdm1(1,ip)-mon%rdm1(1,iq))*B%matX(i,NU)!EigX0(pq,1)
               TRDM1(iq,ip)=TRDM1(iq,ip)-(mon%rdm1(1,iq)-mon%rdm1(1,ip))*B%matY(i,NU)
               ! WORKS ONLY WITHOUT XTilde-->X CHANGE
               ! TRDM1(ip,iq) = TRDM1(ip,iq)+0.5d0*(B%matY(i,exc) * (CICoef(ip)+CICoef(iq)) - B%matX(i,exc) * (CICoef(ip)-CICoef(iq)))
               ! TRDM1(iq,ip) = TRDM1(iq,ip)+0.5d0*(B%matY(i,exc) * (CICoef(ip)+CICoef(iq)) + B%matX(i,exc) * (CICoef(ip)-CICoef(iq)))
            do ir = 1,dimO
               do is = 1,dimO
                   do it =1,dimO
                  montrdm24(ip,is,ir,it) =  montrdm24(ip,is,ir,it)+B%matX(i,NU)*RDM2(iq,is,ir,it)
                  montrdm24(ir,is,ip,it) =  montrdm24(ir,is,ip,it)+B%matX(i,NU)*RDM2(ir,is,iq,it)
                  montrdm24(ir,ip,is,it) =  montrdm24(ir,ip,is,it)-B%matY(i,NU)*RDM2(ir,iq,is,it)
                  montrdm24(ir,it,is,ip) =  montrdm24(ir,it,is,ip)-B%matY(i,NU)*RDM2(ir,it,is,iq)
                  if (ip <= dimO) then
                     montrdm24(ir,iq,is,it) =  montrdm24(ir,iq,is,it)-B%matX(i,NU)*RDM2(ir,ip,is,it)
                     montrdm24(ir,it,is,iq) =  montrdm24(ir,it,is,iq)-B%matX(i,NU)*RDM2(ir,it,is,ip)
                     montrdm24(iq,is,ir,it) =  montrdm24(iq,is,ir,it)+B%matY(i,NU)*RDM2(ip,is,ir,it)
                     montrdm24(ir,is,iq,it) =  montrdm24(ir,is,iq,it)+B%matY(i,NU)*RDM2(ir,is,ip,it)
                  endif

                  enddo
               enddo
            enddo
         enddo
      end associate
   enddo
!unpack (2, IV part)
! associate(B => SblockIV)
!   do i=1,B%n
!      ii = B%l1+i-1
!      ipos = B%pos(i)
!      EigY0(ipos,ii) = 1d0/sqrt(2d0)
!      EigX0(ipos,ii) = 1d0/sqrt(2d0)
!   enddo
!end associate

! do pq=1,mon%NDimX
!    ip = mon%IndN(1,pq)
!    iq = mon%IndN(2,pq)
!          TRDM1(ip,iq)=TRDM1(ip,iq)-(mon%rdm1(1,ip)-mon%rdm1(1,iq))*EigX0(pq,1)
!          TRDM1(iq,ip)=TRDM1(iq,ip)-(mon%rdm1(1,iq)-mon%rdm1(1,ip))*EigY0(pq,1)
! enddo

print*,"1TRDM z X i Y z procedury sapt_ERPA_TRDMs dla monomeru", Mon%Monomer
do i=1,10
   write(lout,'(*(f12.6))') (TRDM1(i,j),j=1,10)
enddo

! do i=1,NBasis
!    do j=1,NBasis
!    if(abs(TRDM1(i,j))>1d-2) then
!       print *,"TRDM1 duze dla :",i, j, "wynosi", TRDM1(i,j)
!    endif
!    enddo
! enddo   


! print*,"2TRDM"
! do ip=1,10
!    do iq = 1,10
!       do ir = 1,10
!          do is =1,10
!             if(abs(montrdm24(exc,ip,iq,ir,is))>1d-2) then
!                print*,ip,iq,ir,is,montrdm24(exc,ip,iq,ir,is)
!             endif
!          enddo
!       enddo
!    enddo
! enddo
! print*,"2TRDM norm", norm2(montrdm24(exc,:,:,:,:))


TRDM1TEST=0d0
do ip=1,NBasis
   do ir=1,NBasis
      do iq=1,NBasis
         TRDM1TEST(ip,ir)= TRDM1TEST(ip,ir)+montrdm24(ip,ir,iq,iq)
      enddo
   enddo
enddo
!  11.12.2023
TRDM1TEST(:,:)=TRDM1TEST(:,:)/(mon%ZNucl-1)
! KONIEC ZMIANY
!TRDM1TEST(:,:)=TRDM1TEST(:,:)/(2*mon%NELE-1)
!ZMIENIONY KOD
print *,"2*mon%NELE"
print *,2*mon%NELE
print*,"1TRDM odtworzny z TRDM"
do i=1,10
   write(lout,'(*(f12.6))') (TRDM1TEST(i,j),j=1,10)
enddo

print*,"1TRDM CAS"
do i=1,10
   write(lout,'(*(f12.6))') (TRDM1CAS(i,j),j=1,10)
enddo

SumTrdmTrdm=0d0
TRDM1CASnorm=0d0
TRDM1norm=0d0

do i=1,NBasis
   do j=1,NBasis
      SumTrdmTrdm=SumTrdmTrdm+TRDM1CAS(i,j)*TRDM1TEST(i,j)
      TRDM1CASnorm=TRDM1CASnorm+TRDM1CAS(i,j)*TRDM1CAS(i,j)
      TRDM1norm=TRDM1norm+TRDM1TEST(i,j)*TRDM1TEST(i,j)
   enddo
enddo   
TRDM1CASnorm=sqrt(TRDM1CASnorm)
TRDM1norm=sqrt(TRDM1norm)
SumTrdmTrdm=SumTrdmTrdm/(TRDM1CASnorm*TRDM1norm)
print*,"Nakładanie 1 TRDM CAS i 1TRDM z 2TRDM"
write(lout,'(*(f12.6))') (SumTrdmTrdm)

if(sign(1d0,SumTrdmTrdm) .NE. sign(1d0,1d0)) then
   print *,"CHANGING SIGN OF 2TRDM"
   montrdm24(:,:,:,:)= -1d0*montrdm24(:,:,:,:)
   TRDM1TEST(:,:) = -1d0*TRDM1TEST(:,:)
  endif  


  print*,'BLad wzgledny'
  write(lout,'(*(f12.6))') (norm2(Trdm1CAS-TRDM1TEST)/norm2(TRDM1CAS))
  




print*,"Norma 2TRDM"
print*,norm2(montrdm24(:,:,:,:))


block
   double precision :: DipXao(NBasis**2),DipYao(NBasis**2),DipZao(NBasis**2)
   double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
   double precision:: TSDipXYZ(3)
   character(:),allocatable     :: mname,dipfile,mofile
   
   if (Mon%Monomer == 1) then
      mname  = 'A'
      dipfile= "DIP_A"
      mofile = 'MOLPRO_A.MOPUN'
   elseif (Mon%Monomer == 2) then
      mname  = 'B'
      dipfile= "DIP_B"
      mofile = 'MOLPRO_B.MOPUN' 
   endif

   call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

!print*, 'DipX-AO:',norm2(DipXao)
!print*, 'DipY-AO:',norm2(DipYao)
!print*, 'DipZ-AO:',norm2(DipZao)

   call unpack_sym_molpro(DipXao,dipfile,NBasis)
   call unpack_sym_molpro(DipYao,dipfile,NBasis)
   call unpack_sym_molpro(DipZao,dipfile,NBasis)

   call tran2MO(DipXao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipX,NBasis)
   call tran2MO(DipYao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipY,NBasis)
   call tran2MO(DipZao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipZ,NBasis)
   TSDipXYZ = 0d0
   do j=1,NBasis
      do i=1,NBasis
         TSDipXYZ(1) = TSDipXYZ(1) - 2d0*TRDM1TEST(j,i)*DipX(i,j)
         TSDipXYZ(2) = TSDipXYZ(2) - 2d0*TRDM1TEST(j,i)*DipY(i,j)
         TSDipXYZ(3) = TSDipXYZ(3) - 2d0*TRDM1TEST(j,i)*DipZ(i,j)
      enddo
   enddo
   print *,'X,Y,Z moments from X and Y matrix'
   print *,TSDipXYZ
   TESTdipMom=TSDipXYZ(maxloc(abs(TSDipXYZ(1:3)),1))
   print *,'Max dipole moment calculated using 1TRDM made from 2TRDM', TESTdipMom
end block





block
   double precision :: DipXao(NBasis**2),DipYao(NBasis**2),DipZao(NBasis**2)
   double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
   double precision :: montrdm(NBasis,NBasis)
   double precision:: TSDipXYZ(3)
   character(:),allocatable     :: mname,dipfile,mofile
   
   if (Mon%Monomer == 1) then
      mname  = 'A'
      dipfile= "DIP_A"
      mofile = 'MOLPRO_A.MOPUN'
   elseif (Mon%Monomer == 2) then
      mname  = 'B'
      dipfile= "DIP_B"
      mofile = 'MOLPRO_B.MOPUN' 
   endif

   call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

!print*, 'DipX-AO:',norm2(DipXao)
!print*, 'DipY-AO:',norm2(DipYao)
!print*, 'DipZ-AO:',norm2(DipZao)
   call tran2MO(mon%trdm1(iexcited,:,:),mon%CMONO(iref,:,:),mon%CMONO(iref,:,:), &
                        montrdm(:,:),NBasis)


   call unpack_sym_molpro(DipXao,dipfile,NBasis)
   call unpack_sym_molpro(DipYao,dipfile,NBasis)
   call unpack_sym_molpro(DipZao,dipfile,NBasis)

   call tran2MO(DipXao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipX,NBasis)
   call tran2MO(DipYao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipY,NBasis)
   call tran2MO(DipZao,mon%CAONO(mon%iref1,:,:),mon%CAONO(mon%iref1,:,:),DipZ,NBasis)
   TSDipXYZ = 0d0
   do j=1,NBasis
      do i=1,NBasis
         TSDipXYZ(1) = TSDipXYZ(1) - 2d0*montrdm(j,i)*DipX(i,j)
         TSDipXYZ(2) = TSDipXYZ(2) - 2d0*montrdm(j,i)*DipY(i,j)
         TSDipXYZ(3) = TSDipXYZ(3) - 2d0*montrdm(j,i)*DipZ(i,j)
      enddo
   enddo
   print *,'X,Y,Z moments from true trdm'
   print *,TSDipXYZ
   TRUEdipMom=TSDipXYZ(maxloc(abs(TSDipXYZ(1:3)),1))
   print*,'True max dip moment',TRUEdipMom
end block



! ZMIANA 18.12.2023
! if(sign(1d0,TESTdipMom) .NE. sign(1d0,TRUEdipMom)) then
!  print *,"CHANGING SIGN OF 2TRDM"
!  montrdm24(:,:,:,:)= -1d0*montrdm24(:,:,:,:)
!  TRDM1TEST(:,:) = -1d0*TRDM1TEST(:,:)
! endif   
! KONIEC ZMIANY



mon%trdm24(:,:,:,:)=montrdm24(:,:,:,:)


!!added 16.06.2023 FOR TEST PURPOUSE
! block
! double precision :: TRDM1MO(NBasis,NBasis)
! call tran2MOt(TRDM1TEST,mon%CMONO(iref,:,:),mon%CMONO(iref,:,:),TRDM1MO,NBasis)!
!
!mon%trdm1(iexcited,:,:) = TRDM1MO(:,:)
!
! endblock
!! end 16.06.2023





end subroutine sapt_ERPA_TRDMs





subroutine init_pino(NBas,Mon,ICASSCF)
implicit none

integer,intent(in) :: NBas, ICASSCF
type(SystemBlock) :: Mon
integer :: NInte1
integer :: i,j,ij,ind

 NInte1 = NBas*(NBas+1)/2

 ! set IGem
 Mon%NGem=1
 do j=1,NBas
    Mon%IGem(j)=1
    if(Mon%Occ(1)==0) Mon%IGem(j) = 2
    if(Mon%Occ(1)==0) Mon%NGem = 2
 enddo

 call SaptInter(NBas,Mon,ICASSCF)

! look-up table
 do i=1,Mon%NELE
   Mon%IndAux(i)=0
 enddo
 do i=1+Mon%NELE,NBas
   Mon%IndAux(i)=2
 enddo

 Mon%NAct = 1
 if(Mon%NAct.ne.0) then
   Mon%icnt = 0
   do i=1,NBas
      if(Mon%Occ(i).gt.0d0) then
         Mon%IndAux(i) = 1
         Mon%icnt = Mon%icnt + 1
      endif
   enddo
 endif

 ! set generalized "occupied" = num0 + num1
 ! and "virtual" = num1 + num2 indices
 Mon%num0 = 0
 do i=1,nbas
    if(Mon%IndAux(i)/=0) exit
    Mon%num0 = Mon%num0 + 1
 enddo
 Mon%num2 = 0
 do i=nbas,1,-1
    if(Mon%IndAux(i)/=2) exit
       Mon%num2 = Mon%num2 + 1
 enddo
 Mon%num1 = nbas - Mon%num0 - Mon%num2

! do i=1,NBas
!    if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
! enddo

 ij=0
 ind=0
 do i=1,NBas
    do j=1,i-1

    ij=ij+1

    if(Mon%IndAux(i)+Mon%IndAux(j).Ne.0.and.Mon%IndAux(i)+Mon%IndAux(j).Ne.4) then

       if((Mon%IGem(i).ne.Mon%IGem(j)).And.(Mon%IndAux(i).Eq.1).And.(Mon%IndAux(j).Eq.1) &
         .and.(Abs(Mon%Occ(i)-Mon%Occ(j))/Mon%Occ(i).Lt.1.D-2) ) Then

          write(*,*)"Discarding nearly degenerate pair",i,j

       else

       ind=ind+1
       Mon%IndX(ind)=ind!ij
       Mon%IndN(1,ind)=i
       Mon%IndN(2,ind)=j

       endif
    endif

    enddo
 enddo

 Mon%NDimX = ind
! write(LOUT,*) Mon%NDim,Mon%NDimX,Mon%NDimN

! for OptTwoP procedures
 allocate(Mon%IndNT(2,NInte1))
 Mon%IndNT = 0
 ij = 0
 do j=1,NBas
    do i=1,j
       ij = ij + 1
       Mon%IndNT(1,ij) = j
       Mon%IndNT(2,ij) = i
    enddo
 enddo

 write(LOUT,*) 'init PINO!!!'

end subroutine init_pino

subroutine calc_fci(NBas,NInte1,NInte2,XOne,URe,TwoMO,Mon,iPINO)
!
! Purpose: obtain exact 2RDM for 2-electron systems
!
implicit none

integer,intent(in) :: NBas, NInte1, NInte2, iPINO
double precision   :: URe(NBas,NBas),XOne(NInte1),TwoMO(NInte2)
type(SystemBlock)  :: Mon
double precision   :: ETot
double precision, allocatable :: Ctmp(:,:),Dtmp(:,:), &
                                 NSymMO(:), TwoMOt(:)
integer :: EigNum
! testy?
integer :: INU,IAB,IA,IB
double precision :: Factor,PNorm

 write(LOUT,'(/,1x,a)') 'ENTERING FCI FOR 2-EL SYSTEMS...'

 allocate(Ctmp(NBas,NBas),Dtmp(NBas,NBas),NSymMO(NBas),TwoMOt(NInte2))
 ! new for OptTwoP
 allocate(Mon%AP(NInte1,NInte1),Mon%PP(NInte1))

 Mon%AP = 0
 Mon%PP = 0

 NSymMO = 1
 ETot=0

 EigNum = Mon%EigFCI
 !print*, 'EigFCI-A',Mon%EigFCI
 !if(Mon%Monomer==1) then
 !   EigNum=1
 !elseif(Mon%Monomer==2) then
 !   EigNum=2
 !endif

! print*, 'IGEM:',Mon%IGem
! call OptTwo1(ETot,Mon%PotNuc,URe,Mon%Occ,XOne,TwoMO,NSymMO, &
!              Mon%CICoef,NBas,NInte1,NInte2,EigNum)

 call OptTwoP(ETot,Mon%PotNuc,URe,Mon%Occ,Mon%AP,Mon%PP,XOne,TwoMO,NSymMO, &
              Mon%CICoef,NBas,NInte1,NInte2,EigNum)

! Do INU=1,NInte1
! PNorm=0d0
! IAB=0
! Do IA=1,NBas
! Do IB=1,IA
!    IAB=IAB+1
!    Factor=2d0
!    If(IA.Eq.IB) Factor=1d0
!       PNorm=PNorm+Factor*Mon%AP(INU,IAB)**2
! EndDo
! EndDo
! Write(*,*)'INU, ExcitEn, PNorm',INU,Mon%PP(INU),PNorm
! EndDo
!

 print*, 'NInte1',NInte1
 print*, 'AP:', norm2(Mon%AP)
 print*, 'PP:', norm2(Mon%PP)

 ! Transform orbitals from MO to NO
 call dgemm('N','T',NBas,NBas,NBas,1d0,Mon%CMO,NBas,URe,NBas,0d0,Ctmp,NBas)

 Dtmp=0
 ! call get_den(NBas,Ctmp,2d0*Mon%Occ,Dtmp)
 ! print*, 'Trace test',trace(Dtmp,NBas)
 ! call print_mo(Ctmp,NBas,'MONOMER X')

 ! Save NO orbitals
 Mon%CMO = Ctmp

 ! MO->NO transform 2-el ints
 TwoMOt(1:NInte2) = 0
 call TwoNO(TwoMOt,URe,TwoMO,NBas,NInte2)

 ! Save transformed ints
 TwoMO = TwoMOt

 deallocate(TwoMOt,NSymMO,Ctmp,Dtmp)

end subroutine calc_fci

subroutine set_filenames(Mon,FNam)
implicit none

type(SystemBlock) :: Mon
type(FileNames)   :: FNam

 if(Mon%Monomer==1) then
    FNam%sirifile   = 'SIRIFC_A'
    FNam%siriusfile = 'SIRIUS_A.RST'
    FNam%coefile    = 'coeff_A.dat'
    FNam%occfile    = 'occupations_A.dat'
    FNam%onefile    = 'ONEEL_A'
    FNam%twofile    = 'TWOMOAA'
    FNam%twojfile   = 'FFOOAA'
    FNam%twokfile   = 'FOFOAA'
    FNam%abfile     = 'ABMAT_A'
    FNam%propfile   = 'PROP_A'
    FNam%propfile0  = 'PROP_A0'
    FNam%propfile1  = 'PROP_A1'
    FNam%propfile2  = 'PROP_A2'
    FNam%y01file    = 'Y01_A'
    FNam%xy0file    = 'XY0_A'
    FNam%rdmfile    = 'rdm2_A.dat'
    FNam%testfile   = 'A0MAT_A'
 elseif(Mon%Monomer==2) then
    FNam%sirifile   = 'SIRIFC_B'
    FNam%siriusfile = 'SIRIUS_B.RST'
    FNam%coefile    = 'coeff_B.dat'
    FNam%occfile    = 'occupations_B.dat'
    FNam%onefile    = 'ONEEL_B'
    FNam%twofile    = 'TWOMOBB'
    FNam%twojfile   = 'FFOOBB'
    FNam%twokfile   = 'FOFOBB'
    FNam%abfile     = 'ABMAT_B'
    FNam%propfile   = 'PROP_B'
    FNam%propfile0  = 'PROP_B0'
    FNam%propfile1  = 'PROP_B1'
    FNam%propfile2  = 'PROP_B2'
    FNam%y01file    = 'Y01_B'
    FNam%xy0file    = 'XY0_B'
    FNam%rdmfile    = 'rdm2_B.dat'
    FNam%testfile   = 'A0MAT_B'
 endif

end subroutine set_filenames

subroutine get_1el_h_mo(XOne,MO,NBas,onefile)
! read 1-el Hamiltonian and transform with MO coefs
implicit none

integer,intent(in)          :: NBas
character(*),intent(in)     :: onefile
double precision,intent(in) :: MO(NBas*NBas)
double precision,intent(in) :: XOne(NBas*(NBas+1)/2)
double precision,allocatable :: work1(:),work2(:)

integer      :: iunit
character(8) :: label

 allocate(work1(NBas**2),work2(NBas**2))
 open(newunit=iunit,file=onefile,access='sequential',&
      form='unformatted',status='old')

 read(iunit)
 read(iunit)
 read(iunit) label, work1
 if(label=='ONEHAMIL') then
    call tran_oneint(work1,MO,MO,work2,NBas)
    call sq_to_triang(work1,XOne,NBas)
 else
    write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
    stop
 endif

 close(iunit)
 deallocate(work2,work1)

end subroutine get_1el_h_mo

subroutine create_symmats(Mon,MO,NBasis)
implicit none
! create MultpC and NSymNO 
! WARNING!!! THIS PROCEDURE STILL REQUIRES VERIFICATION
! WITH A SAPT LR-PBE JOB!!!

integer,intent(in) :: NBasis
double precision,intent(in) :: MO(NBasis,NBasis)
type(SystemBlock) :: Mon

integer :: i,j,iorb,istart

allocate(Mon%NSymNO(NBasis),Mon%MultpC(15,15))

!print*, 'NSym',Mon%NSym
Mon%MultpC = 0
if(Mon%NSym==1) then
   Mon%MultpC(1,1)=1
else
   do i=1,Mon%NSym
      do j=1,i
         Mon%MultpC(i,j) = ieor(i-1,j-1)+1
         Mon%MultpC(j,i) = Mon%MultpC(i,j)
      enddo
   enddo
endif


Mon%NSymNO(1:NBasis) = 0
istart = 0
do i=1,Mon%NSym
   do j=istart+1,istart+Mon%NumOSym(i)
      do iorb=1,NBasis
         if(abs(MO(j,iorb)).gt.1.d-1) Mon%NSymNO(iorb) = i
      enddo
   enddo
   istart=istart+Mon%NumOSym(i)
enddo

! check
do i=1,Mon%NSym
   j = 0
   do iorb=1,NBasis
      if(Mon%NSymNO(iorb)==i) j = j + 1
   enddo
   if(j/=Mon%NumOSym(i)) then
      write(LOUT,'(1x,a)') 'ERROR IN create_symmats! SYMMETRY OF NO CANNOT BE ESTABLISHED!'
      stop
   endif
enddo
!
end subroutine create_symmats

subroutine reduce_dim(var,matP,matM,Mon)
implicit none

character(*) :: var
double precision :: matP(:),matM(:)
type(SystemBlock) :: Mon
integer :: i,j,ij,ij1

select case(var)
case('AB','ab')
   ! matP = ABPlus
   ! matM = ABMin
   ij=0
   ij1=0
   do j = 1,Mon%NDimX
      do i = 1,Mon%NDimX
          ij  = (j-1)*Mon%NDimX + i
          ij1 = (Mon%IndX(j)-1)*Mon%NDim + Mon%IndX(i)
          matP(ij) = matP(ij1)
          matM(ij) = matM(ij1)
      enddo
   enddo

case('D','d')
   ! matP = DMAT
   ! matM = DMATK
   ij=0
   ij1=0
   do j = 1,Mon%NDimN
      do i = 1,Mon%NDimX
         ij  = ij + 1 !(j-1)*Mon%NDimX + i
         ij1 = (j-1)*Mon%NDim + Mon%IndX(i)
         matP(ij) = matP(ij1)
         matM(ij) = matM(ij1)
      enddo
   enddo

case('E','e')
   ! matP = EMAT
   ! matM = EMATM
   ij=0
   ij1=0
   do j=1,Mon%NDimN
      do i=1,Mon%NDimN
         ij  = (j-1)*Mon%NDimN + i
         ij1 = (j-1)*Mon%NBasis + i
         matP(ij) = matP(ij1)
         matM(ij) = matM(ij1)
      enddo
   enddo

end select

end subroutine reduce_dim

subroutine writeResp(EVecZ,EVal,mon)
implicit none

character(*) :: mon
double precision :: EVecZ(:), EVal(:)
integer :: iunit

 open(newunit=iunit,file=mon,form='unformatted')
 write(iunit) EVecZ
 write(iunit) EVal
 close(iunit)

end subroutine writeResp

subroutine writeRespXY(EVecX,EvecY,EVal,mon)
implicit none

character(*)     :: mon
double precision :: EVecX(:),EVecY(:),EVal(:)

integer :: iunit

open(newunit=iunit,file=mon,form='unformatted')
write(iunit) EVecX
write(iunit) EVecY
write(iunit) EVal
close(iunit)

end subroutine writeRespXY

subroutine writeEval(EVal,mon)
implicit none

character(*) :: mon
double precision :: EVal(:)
integer :: iunit

 open(newunit=iunit,file=mon,form='unformatted')
 write(iunit) EVal
 close(iunit)

end subroutine writeEval

end module sapt_resp
