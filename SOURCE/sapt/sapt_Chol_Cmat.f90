module sapt_Chol_Cmat

use tran, only : ABPM_HALFTRAN_GEN_L
use sapt_utils
use timing

implicit none

contains

subroutine e2disp_Cmat(Flags,A,B,SAPT)
!
! THIS IS FOR TESTING ONLY:
! calculate 2nd order dispersion energy
! use C(omega) obained brute-force (i.e.,
! by inversion: C(om)=(APLUS.AMIN + om^2)^-1.APLUS
! in coupled and uncoupled approximations
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
double precision :: fact,val
double precision :: Omega,Pi,e2d
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: work(:,:),ints(:)

double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

! check Cholesky
if(Flags%ICholesky==1) then
   write(LOUT,'(1x,a)') 'Cholesky Cmat not ready yet! Aborting...'
   return
endif

! set dimensions
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

Pi = 4.0d0*atan(1.0)

! get ABPM = ABPLUS.ABMIN
allocate(ABPMA(A%NDimX,A%NDimX),ABPMB(B%NDimX,B%NDimX),&
         ABPLUSA(A%NDimX,A%NDimX),ABPLUSB(B%NDimX,B%NDimX))
allocate(work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)

deallocate(work)
allocate(work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)

deallocate(work)

! frequency integration
NFreq = 18
allocate(XFreq(NFreq),WFreq(NFreq))

call FreqGrid(XFreq,WFreq,NFreq)

allocate(CB(B%NDimX,B%NDimX))
allocate(ints(NBas**2),work(A%NDimX,B%NDimX))

e2d = 0
do ifreq=1,NFreq

   Omega = XFreq(ifreq)

   allocate(CA(A%NDimX,A%NDimX))

   call get_Cmat(CA,A%CICoef,A%IndN,ABPMA,ABPLUSA,Omega,A%NDimX,NBas)
   call get_Cmat(CB,B%CICoef,B%IndN,ABPMB,ABPLUSB,Omega,B%NDimX,NBas)

   open(newunit=iunit,file='TWOMOAB',status='OLD',&
        access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

   ints = 0
   work = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)

         do i=1,A%NDimX
            work(i,irs) = work(i,irs) + fact*CA(ipq,i)
         enddo

      enddo
   enddo

   deallocate(CA)
   allocate(CA(B%NDimX,B%NDimX))
   CA = 0
   ints = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)
         do j=1,B%NDimX
            CA(irs,j) = CA(irs,j) + fact*work(ipq,j)
         enddo

      enddo

   enddo

   val = 0
   do j=1,B%NDimX
      do i=1,B%NDimX
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo
   e2d = e2d + WFreq(ifreq)*val

   close(iunit)
   deallocate(CA)

enddo ! NFreq

SAPT%e2disp = -8d0/Pi*e2d
e2d = -8d0/Pi*e2d*1d3

call print_en('E2disp',e2d,.true.)

deallocate(ints,work)
deallocate(WFreq,XFreq)
deallocate(CB)
deallocate(ABPMB,ABPMA,ABPLUSB,ABPLUSA)

end subroutine e2disp_Cmat

subroutine e2disp_Cmat_Chol(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES ADAM'S DIAG/PROJECT
! WHICH ASSUME FULL A0(NDIMX,NDIMX) MATRICES
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats
use class_IterAlgorithm
use class_IterAlgorithmDIIS
use class_LambdaCalculator
use class_LambdaCalculatorDiag
use class_LambdaCalculatorProjector

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d
logical          :: project

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: AB0A(:,:),AB0B(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: work(:,:),work1(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()
class(IterAlgorithmDIIS), allocatable, target    :: iterAlgo
!class(LambdaCalculatorDiag), allocatable, target :: LambdaCalcA,LambdaCalcB
class(LambdaCalculatorProjector), allocatable, target :: LambdaCalcA,LambdaCalcB

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! use projection
project=.true.
!project=.false.

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! get Pmat
print*, 'to-do: adapt Pmat procedure for SAPT...'

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
!print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
!print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq
write(lout,'(/1x,a)',advance='no') 'E2disp(Cmat) with projection: '
write(lout,'(1x,a,i3)') 'Adams diag version'

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NDimX))
allocate(AB0A(A%NDimX,A%NDimX),AB0B(B%NDimX,B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

!iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20)
allocate(iterAlgo,SOURCE = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20))

if(project) then
   !LambdaCalcA = LambdaCalculatorProjector(A%NDimX,ABPMA,A%Pmat)
   !LambdaCalcB = LambdaCalculatorProjector(B%NDimX,ABPMB,B%Pmat)
   allocate(LambdaCalcA, SOURCE = LambdaCalculatorProjector(A%NDimX,ABPMA,A%Pmat))
   allocate(LambdaCalcB, SOURCE = LambdaCalculatorProjector(B%NDimX,ABPMB,B%Pmat))
else
 !  LambdaCalcA = LambdaCalculatorDiag(A%NDimX,ABPMA)
 !  LambdaCalcB = LambdaCalculatorDiag(B%NDimX,ABPMB)
endif

iStatsA%maxIterationsLimit = iterAlgo%maxIterations
iStatsB%maxIterationsLimit = iterAlgo%maxIterations

call LambdaCalcA%calculateInitialA()
call LambdaCalcB%calculateInitialA()

do i=1,A%NDimX
   val = ABPMA((i-1)*A%NDimX+i)
   ABPMA((i-1)*A%NDimX+i) = ABPMA((i-1)*A%NDimX+i) - val
enddo
do i=1,B%NDimX
   val = ABPMB((i-1)*B%NDimX+i)
   ABPMB((i-1)*B%NDimX+i) = ABPMB((i-1)*B%NDimX+i) - val
enddo
if(project) then
   allocate(work1(A%NDimX*A%NDimX))
   call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,A%PMat,A%NDimX,ABPMA,A%NDimX,0d0,work1,A%NDimX)
   ABPMA = ABPMA - work1
   deallocate(work1)

   allocate(work1(B%NDimX*B%NDimX))
   call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,B%PMat,B%NDimX,ABPMB,B%NDimX,0d0,work1,B%NDimX)
   ABPMB = ABPMB - work1
   deallocate(work1)
endif

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call LambdaCalcA%calculateLambda(LambdaA, OmI)
   call LambdaCalcB%calculateLambda(LambdaB, OmI)

   if(ifreq==NFreq) call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,ABPTildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   call iterAlgo%iterate(CTildeA, A%NDimX, NCholesky, LambdaA, ABPTildeA, ABPMA, iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,ABPTildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   call iterAlgo%iterate(CTildeB, B%NDimX, NCholesky, LambdaB, ABPTildeB, ABPMB, iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

e2d = -8d0/Pi*e2d*1d3
!print*, 'E2disp(Cmat)',e2d
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(ABPMB,ABPMA)
deallocate(AB0B,AB0A)
deallocate(ABPTildeB,ABPTildeA)
deallocate(LambdaB,LambdaA,XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol

subroutine e2disp_Cmat_Chol_diag(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES MY DIAG PROCEDURES
! WHICH KEEP AND MULTIPLY ONLY NDIMX ELEMENTS
! PMAT PROJECTION NOT AVAILABLE!
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(SaptDIIS)    :: SAPT_DIIS

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: DiagA(:),DiagB(:)
double precision,allocatable :: LamDiaA(:),LamDiaB(:)
double precision,allocatable :: work(:,:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
! making use of diagonal Lambda matrices

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(DiagA(A%NDimX),DiagB(B%NDimX))
allocate(LamDiaA(A%NDimX),LamDiaB(B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

iStatsA%maxIterationsLimit = SAPT_DIIS%maxIter
iStatsB%maxIterationsLimit = SAPT_DIIS%maxIter

call calculateInitialA_diag(DiagA,ABPMA,A%NDimX)
call calculateInitialA_diag(DiagB,ABPMB,B%NDimX)

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_diag(LamDiaA,DiagA,OmI,A%NDimX)
   call calculateLambda_diag(LamDiaB,DiagB,OmI,B%NDimX)

   if(ifreq==NFreq) call MultpDiagMat(LamDiaA,ABPTildeA,0d0,CTildeA,A%NDimX,NCholesky)
   call Cmat_diag_iterDIIS(CTildeA,A%NDimX,NCholesky,LamDiaA,ABPTildeA,ABPMA,iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call MultpDiagMat(LamDiaB,ABPTildeB,0d0,CTildeB,B%NDimX,NCholesky)
   call Cmat_diag_iterDIIS(CTildeB,B%NDimX,NCholesky,LamDiaB,ABPTildeB,ABPMB,iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

SAPT%e2disp = -8d0/Pi*e2d
e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)
deallocate(DiagB,DiagA)
deallocate(LamDiaB,LamDiaA)
deallocate(XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol_diag

subroutine e2disp_Cmat_Chol_block(Flags,A,B,SAPT)
!
! THIS PROCEDURE SHOULD BE IMPROVED:
! A0 MATRICES ARE KEPT IN DIAGONAL BLOCKS
! AND BLOCK MULTIPLICATION IS USED
! (THE ABPM_HALFTRAN_LR PROCEDURE
! IN sapt_utils.f90);
! IDEALLY, THE BLOCK/DIAGONAL MULTIPLICATIONS
! COMMON FOR THE ENTIRE GAMMCOR SHOULD BE USED
! AND DIAG/BLOCK VERSION OF THE ALGORITHM SHOULD
! BE CHOSEN AUTOMATICALLY, E.G., BASED ON NACT
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! using A0 blocks (not diagonal)
! with Cholesky vectors

use class_IterStats

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n,nblkA,nblkB
integer :: maxIter
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: work(:,:)
double precision,allocatable :: WorkA(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)
type(EBlockData)             :: LambdaIVA,LambdaIVB
type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX,A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX,B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
NFreq = 12
write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! get A2 = ABPM - ABPM0
call calculateInitialA_blk(ABPMA,A0BlkA,A0BlkIVA,nblkA,A%NDimX,'A0BLK_A')
call calculateInitialA_blk(ABPMB,A0BlkB,A0BlkIVB,nblkB,B%NDimX,'A0BLK_B')

maxIter = 20
iStatsA%maxIterationsLimit = maxIter
iStatsB%maxIterationsLimit = maxIter

allocate(LambdaA(nblkA),LambdaB(nblkB))

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_blk(LambdaA,LambdaIVA,OmI**2,A%NDimX,nblkA,A0BlkA,A0BlkIVA)
   call calculateLambda_blk(LambdaB,LambdaIVB,OmI**2,B%NDimX,nblkB,A0BlkB,A0BlkIVB)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeA,CTildeA,0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')
   call Cmat_blk_iterDIIS(CTildeA,A%NDimX,NCholesky,nblkA,LambdaA,LambdaIVA,ABPTildeA,ABPMA,iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeB,CTildeB,0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')
   call Cmat_blk_iterDIIS(CTildeB,B%NDimX,NCholesky,nblkB,LambdaB,LambdaIVB,ABPTildeB,ABPMB,iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

   call release_ac0block(LambdaA,LambdaIVA,nblkA)
   call release_ac0block(LambdaB,LambdaIVB,nblkB)

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

SAPT%e2disp  = -8d0/Pi*e2d

e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(CB,CA)
deallocate(CTildeB,CTildeA)
deallocate(WFreq,XFreq)
!deallocate(DCholB,DCholA)
deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)

end subroutine e2disp_Cmat_Chol_block

subroutine e2disp_Cmat_Chol_proj(Flags,A,B,SAPT)
!
! THIS PROCEDURE USES MY PROJECT PROCEDURES
! BUT THERE IS PERHAPS NO SENSE OF DOING THAT?
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! with Cholesky vectors

use class_IterStats
use class_IterAlgorithm
use class_IterAlgorithmDIIS

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(SaptDIIS)    :: SAPT_DIIS

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ij,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d
logical :: diag, blk

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: PMatA(:,:),PmatB(:,:)
double precision,allocatable :: ABPMA(:),ABPMB(:),&
                                AB0A(:,:),AB0B(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:,:),CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: work(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()
class(IterAlgorithmDIIS), allocatable, target    :: iterAlgo

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! diagonal or block
diag = .true.
blk  = .false.
!diag = .false.
!blk  = .true.

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPMA(A%NDimX*A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,A%DChol,NCholesky,0d0,ABPTildeA,A%NDimX)
print*, 'ABPTildeA',norm2(ABPTildeA)

deallocate(ABPLUSA,work)

! monomer B
allocate(ABPMB(B%NDimX*B%NDimX),ABPTildeB(B%NDimX,NCholesky))
allocate(ABPLUSB(B%NDimX,B%NDimX),work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,B%DChol,NCholesky,0d0,ABPTildeB,B%NDimX)
print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
! making use of diagonal Lambda matrices

allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
write(lout,'(/1x,a)',advance='no') 'E2disp(Cmat) with projection: '
if(diag) write(lout,'(1x,a,i3)') 'diagonal version'
if(blk)  write(lout,'(1x,a,i3)') 'block version'

write(lout,'(/1x,a,i3)') 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(AB0A(A%NDimX,A%NDimX),AB0B(B%NDimX,B%NDimX))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

!iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20)
allocate(iterAlgo, SOURCE = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20))

iStatsA%maxIterationsLimit = iterAlgo%maxIterations
iStatsB%maxIterationsLimit = iterAlgo%maxIterations

if(diag) then
  call calculateInitialA_diagP(AB0A,ABPMA,A%PMat,A%NDimX)
  call calculateInitialA_diagP(AB0B,ABPMB,B%PMat,B%NDimX)
elseif(blk) then
  call calculateInitialA_blkP(ABPMA,AB0A,A%PMat,A%NDimX,'A0BLK_A')
  call calculateInitialA_blkP(ABPMB,AB0B,B%PMat,B%NDimX,'A0BLK_B')
endif

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda_Pmat(LambdaA,AB0A,OmI,A%NDimX)
   call calculateLambda_Pmat(LambdaB,AB0B,OmI,B%NDimX)

   if(ifreq==NFreq) call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,ABPTildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   call iterAlgo%iterate(CTildeA, A%NDimX, NCholesky, LambdaA, ABPTildeA, ABPMA, iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,ABPTildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   call iterAlgo%iterate(CTildeB, B%NDimX, NCholesky, LambdaB, ABPTildeB, ABPMB, iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

write(lout,'(/1x,a)') 'C(omega): Monomer A'
call iStatsA%print()
write(lout,'(1x,a)') 'C(omega): Monomer B'
call iStatsB%print()

e2d = -8d0/Pi*e2d*1d3
call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(AB0B,AB0A)
deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)
deallocate(LambdaB,LambdaA)
deallocate(XFreq,WFreq)
deallocate(CB,CA)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol_proj

subroutine get_Cmat(Cmat,CICoef,IndN,ABPM,ABPLUS,Omega,NDimX,NBas)
implicit none

integer,intent(in) :: NDimX,NBas,IndN(2,NBas)
double precision,intent(in)  :: CICoef(NBas),Omega
double precision,intent(in)  :: ABPM(NDimX,NDimX),ABPLUS(NDimX,NDimX)
double precision,intent(out) :: Cmat(NDimX,NDimX)

integer          :: info,lwork
integer          :: i,ipq,irs,ip,iq,ir,is
double precision :: fact_rs,fact_pq
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:,:),work1(:)

allocate(work(NDimX,NDimX),ipiv(NDimX))
work = ABPM
do i=1,NDimX
   work(i,i) = work(i,i) + Omega**2
enddo

! work=(ABPM+omega^2)-1.ABPLUS
!...

! this is slower
!call dgetrf(NDimX,NDimX,work,NDimX,ipiv,info)
!print*, 'info1',info
!allocate(work1(1))
!lwork = -1
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!lwork = int(work1(1))
!print*, 'lwork:',lwork
!deallocate(work1)
!allocate(work1(lwork))
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!print*, 'info2',info
!!
!! Cmat=work.ABPLUS
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,work,NDimX,ABPLUS,NDimX,0d0,Cmat,NDimX)
!
! deallocate(work1)
!
Cmat = ABPlus
call dgesv(NDimX,NDimX,work,NDimX,ipiv,Cmat,NDimX,info)
!print*, 'info',info
!print*, 'Cmat',norm2(Cmat)

!
do irs=1,NDimX
   ir = IndN(1,irs)
   is = IndN(2,irs)
   fact_rs = CICoef(ir)+CICoef(is)
   do ipq=1,NDimX
      ip = IndN(1,ipq)
      iq = IndN(2,ipq)
      fact_pq = CICoef(ip)+CICoef(iq)
      Cmat(ipq,irs) = fact_rs*fact_pq*Cmat(ipq,irs)
   enddo
enddo

deallocate(ipiv,work)

end subroutine get_Cmat


end module sapt_Chol_Cmat
