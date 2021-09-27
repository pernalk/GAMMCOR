module sapt_CD_pol
use types
use tran, only : ABPM_HALFTRAN_GEN_L
use sapt_utils

implicit none

contains

subroutine e1elst_Chol(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i,j,ii,jj
integer :: NBas,NCholesky
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Vabb(:,:),Vbaa(:,:)
double precision :: ea,eb,eab,elst
double precision,external  :: ddot

! set dimensions
 NBas = A%NBasis
 NCholesky = SAPT%NCholesky

 allocate(Va(NBas,NBas),Vb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas))

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

! sum_p n_p v^B_pp
 ea = 0
 do i=1,A%num0+A%num1
    ea = ea + A%Occ(i)*Vbaa(i,i)
 enddo
 ea = 2d0*ea
 !print*, 'ea',ea

! sum_q n_q v^A_qq
 eb = 0
 do j=1,B%num0+B%num1
    eb = eb + B%Occ(j)*Vabb(j,j)
 enddo
 eb = 2d0*eb
 !print*, 'eb',eb

! sum_pq n_p n_q v_{pq}^{pq}
 eab = 0
 do j=1,B%num0+B%num1
    jj = (j-1)*NBas+j
    do i=1,A%num0+A%num1
       ii = (i-1)*NBas+i
       eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%FF(:,ii),1,B%FF(:,jj),1)
       !eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%FO(:,ii),1,B%FO(:,jj),1)
    enddo
 enddo
 eab = 4d0*eab
 !print*, 'eab', eab

 elst = ea + eb + eab + SAPT%Vnn

 call print_en('V_nn',SAPT%Vnn,.false.)
 call print_en('Eelst',elst*1000,.false.)
 SAPT%elst = elst

 deallocate(Vb,Va,Vbaa,Vabb)

end subroutine e1elst_Chol

subroutine e2disp_Chol(Flags,A,B,SAPT)
! calculate 2nd order dispersion energy
! in coupled and uncoupled approximations
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: NCholesky
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
logical,allocatable          :: condOmA(:),condOmB(:)
double precision,allocatable :: OmA(:), OmB(:), &
                                OmA0(:),OmB0(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),&
                                tmp01(:,:),tmp02(:,:)
double precision,allocatable :: work(:)
double precision :: e2d,fact,tmp
double precision :: e2du,dea,deb
double precision :: inv_omega
! for Be ERPA:
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! print thresholds
 if(SAPT%IPrint>1) then
    write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
    write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
    write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
 endif

! set dimensions
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA  = dimOA*dimVA
 nOVB  = dimOB*dimVB

 NCholesky = SAPT%NCholesky

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),  &
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX),  &
          OmA0(A%NDimX),OmB0(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 ! uncoupled - works for CAS only
 if(Flags%ICASSCF==1) then
    allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))

    call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
    call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')
 endif

allocate(work(B%NDimX))

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
        tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX))

! coupled
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative omega A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative omega B!',i,OmB(i)
enddo

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then

 tmp1=0
 tmp01=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

       associate(Y => Y01BlockA(pq))
          tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
       end associate

    enddo
 enddo
 ! coupled
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp2,A%NDimX)

 ! uncoupled
 tmp02=0
 do rs=1,B%NDimX
    associate(Y => Y01BlockB(rs))
      call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,tmp02(:,Y%l1:Y%l2),A%NDimX)
    end associate
 enddo

elseif(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

 tmp1 = 0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j) + &
                    EVecB(rs+(j-1)*B%NDimX)*tmp1(i,rs)
       enddo
    enddo
 enddo

endif ! end GVB select

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then
   ! uncoupled
    e2du = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then


          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          e2du = e2du + tmp02(i,j)**2*inv_omega

          endif
       enddo
    enddo
    SAPT%e2disp_unc = -16d0*e2du

    e2du = -16d0*e2du*1000d0

    call writeampl(tmp02,'PROP_AB0')

endif

 allocate(condOmA(A%NDimX),condOmB(B%NDimX))
 condOmA = (abs(OmA).gt.SmallE.and.abs(OmA).lt.BigE)
 condOmB = (abs(OmB).gt.SmallE.and.abs(OmB).lt.BigE)

 e2d = 0d0
 do j=1,B%NDimX
    if(condOmB(j)) then
       do i=1,A%NDimX
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
!             .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

             if(condOmA(i)) then
                e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))
             endif
       enddo
    endif
 enddo
 SAPT%e2disp  = -16d0*e2d

 e2d  = -16d0*e2d*1000d0

 call print_en('E2disp',e2d,.true.)
 call print_en('E2disp(unc)',e2du,.false.)

 ! write amplitude to a file
 call writeampl(tmp2,'PROP_AB')

 !! calucate semicoupled and dexcitations
 !if(SAPT%SemiCoupled) call e2disp_semi(Flags,A,B,SAPT)

 !! calculate extrapolated E2disp
 !if(A%Cubic.or.B%Cubic) call e2disp_cpld(Flags,A,B,SAPT)

 !! calculate Wterms (deexcitations)
 !if(SAPT%Wexcit) call e2inddisp_dexc(Flags,A,B,SAPT)

 deallocate(work)

 if(Flags%ICASSCF==1) then
    ! deallocate Y01Block
    do i=1,A%NDimX
       associate(Y => Y01BlockA(i))
         deallocate(Y%vec0)
       end associate
    enddo
    do i=1,B%NDimX
       associate(Y => Y01BlockB(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockB,Y01BlockA)
 endif

 deallocate(condOmB,condOmA)
 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,OmA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp_Chol

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
! THIS PROCEDURE SHOULD BE IMPROVED:
! IT IS BASED ON ADAM'S IMPLEMENTATION
! OF \LAMBDA AND ITERATIVE ALGORITHMS
! WHICH ASSUME FULL A0(NDIMX,NDIMX) MATRICES
! SO THAT NEITHER DIAGONAL NOR BLOCK MULTIPLICATION
! ARE USED
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

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: n
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:),&
                                A0A(:,:),A0B(:,:),&
                                ABPTildeA(:,:),ABPTildeB(:,:)
double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: CTildeA(:),CTildeB(:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: work(:,:)
double precision,allocatable :: DCholAT(:,:)

!class(IterDIIS), allocatable, target       :: iterDIISAlgorithm
type(IterStats) :: iStatsA = IterStats()
type(IterStats) :: iStatsB = IterStats()
class(IterAlgorithmDIIS), allocatable, target    :: iterAlgo
class(LambdaCalculatorDiag), allocatable, target :: LambdaCalc
!class(LambdaCalculatorProjector), allocatable, target :: LambdaCalc

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

! get Dmat
allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))

DCholA = 0
do j=1,A%NDimX
   ip = A%IndN(1,j)
   iq = A%IndN(2,j)
   ipq = iq + (ip-1)*NBas
   Cpq = A%CICoef(ip) + A%CICoef(iq)
   do i=1,NCholesky
      DCholA(i,j) = Cpq*A%FF(i,ipq)
   enddo
enddo
DCholB = 0
do j=1,B%NDimX
   ir = B%IndN(1,j)
   is = B%IndN(2,j)
   irs = is + (ir-1)*NBas
   Crs = B%CICoef(ir) + B%CICoef(is)
   do i=1,NCholesky
      DCholB(i,j) = Crs*B%FF(i,irs)
   enddo
enddo

! monomer A
allocate(ABPMA(A%NDimX,A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,DCholA,NCholesky,0.0,ABPTildeA,A%NDimX)

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
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,DCholB,NCholesky,0.0,ABPTildeB,B%NDimX)
!print*, 'ABPTildeB',norm2(ABPTildeB)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner

allocate(CTildeA(A%NDimX*NCholesky),CTildeB(B%NDimX*NCholesky))
allocate(A0A(A%NDimX,A%NDimX),A0B(B%NDimX,B%NDimX))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

NFreq = 12
print*, 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NDimX))

call FreqGrid(XFreq,WFreq,NFreq)

A0A = 0
A0B = 0

iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=6, maxIterations=20)
LambdaCalc = LambdaCalculatorDiag()

iStatsA%maxIterationsLimit = iterAlgo%maxIterations
iStatsB%maxIterationsLimit = iterAlgo%maxIterations

call LambdaCalc%calculateInitialA(A0A, ABPMA, A%NDimX)
call LambdaCalc%calculateInitialA(A0B, ABPMB, B%NDimX)

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call LambdaCalc%calculateLambda(LambdaA, OmI, A%NDimX, A0A)
   call LambdaCalc%calculateLambda(LambdaB, OmI, B%NDimX, A0B)

   if(ifreq==NFreq) call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,ABPTildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   call iterAlgo%iterate(CTildeA, A%NDimX, NCholesky, LambdaA, ABPTildeA, ABPMA, iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,ABPTildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   call iterAlgo%iterate(CTildeB, B%NDimX, NCholesky, LambdaB, ABPTildeB, ABPMB, iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,DCholB,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

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
deallocate(ABPTildeB,ABPTildeA)
deallocate(LambdaB,LambdaA,XFreq,WFreq)
deallocate(CB,CA,A0B,A0A)
deallocate(CTildeB,CTildeA)

end subroutine e2disp_Cmat_Chol

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
double precision,allocatable :: DCholA(:,:),DCholB(:,:)
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

! get Dmat
allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))

DCholA = 0
do j=1,A%NDimX
   ip = A%IndN(1,j)
   iq = A%IndN(2,j)
   ipq = iq + (ip-1)*NBas
   Cpq = A%CICoef(ip) + A%CICoef(iq)
   do i=1,NCholesky
      DCholA(i,j) = Cpq*A%FF(i,ipq)
   enddo
enddo
DCholB = 0
do j=1,B%NDimX
   ir = B%IndN(1,j)
   is = B%IndN(2,j)
   irs = is + (ir-1)*NBas
   Crs = B%CICoef(ir) + B%CICoef(is)
   do i=1,NCholesky
      DCholB(i,j) = Crs*B%FF(i,irs)
   enddo
enddo

! monomer A
allocate(ABPMA(A%NDimX,A%NDimX),ABPTildeA(A%NDimX,NCholesky))
allocate(ABPLUSA(A%NDimX,A%NDimX),work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUSA,A%NDimX,DCholA,NCholesky,0d0,ABPTildeA,A%NDimX)

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
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUSB,B%NDimX,DCholB,NCholesky,0d0,ABPTildeB,B%NDimX)

deallocate(ABPLUSB,work)

! get CTilde in an iterative manner
NFreq = 12
print*, 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! get A2 = ABPM - ABPM0
call calculateInitialA(ABPMA,A0BlkA,A0BlkIVA,nblkA,A%NDimX,'A0BLK_A')
call calculateInitialA(ABPMB,A0BlkB,A0BlkIVB,nblkB,B%NDimX,'A0BLK_B')

maxIter = 20
iStatsA%maxIterationsLimit = maxIter
iStatsB%maxIterationsLimit = maxIter

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

   call calculateLambda(LambdaA,LambdaIVA,OmI**2,A%NDimX,nblkA,A0BlkA,A0BlkIVA)
   call calculateLambda(LambdaB,LambdaIVB,OmI**2,B%NDimX,nblkB,A0BlkB,A0BlkIVB)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeA,CTildeA,0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')
   call Cmat_iterDIIS(CTildeA,A%NDimX,NCholesky,nblkA,LambdaA,LambdaIVA,ABPTildeA,ABPMA,iStatsA)
   call iStatsA%setFreq(OmI)

   if(ifreq==NFreq) call ABPM_HALFTRAN_GEN_L(ABPTildeB,CTildeB,0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')
   call Cmat_iterDIIS(CTildeB,B%NDimX,NCholesky,nblkB,LambdaB,LambdaIVB,ABPTildeB,ABPMB,iStatsB)
   call iStatsB%setFreq(OmI)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,DCholB,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

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
deallocate(DCholB,DCholA)
deallocate(ABPMB,ABPMA)
deallocate(ABPTildeB,ABPTildeA)

end subroutine e2disp_Cmat_Chol_block

subroutine e2disp_CAlphaTilde_block(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in coupled approximation using
! C(omega) obained in an iterative fashion
! using A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: nblkA,nblkB
integer :: N,Max_Cn
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: XFactorial,XN1,XN2
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUS1A(:,:),ABPLUS1B(:,:),&
                                ABMIN1A(:,:),ABMIN1B(:,:),  &
                                A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:),&
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                C1TildeA(:,:),C1TildeB(:,:), &
                                C2TildeA(:,:),C2TildeB(:,:), &
                                CTildeA(:,:), CTildeB(:,:)
double precision,allocatable :: WorkA(:,:),WorkB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

type(EBlockData)             :: LambdaIVA,LambdaIVB
type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)

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

! get Dmat
allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))

DCholA = 0
do j=1,A%NDimX
   ip = A%IndN(1,j)
   iq = A%IndN(2,j)
   ipq = iq + (ip-1)*NBas
   Cpq = A%CICoef(ip) + A%CICoef(iq)
   do i=1,NCholesky
      DCholA(i,j) = Cpq*A%FF(i,ipq)
   enddo
enddo
DCholB = 0
do j=1,B%NDimX
   ir = B%IndN(1,j)
   is = B%IndN(2,j)
   irs = is + (ir-1)*NBas
   Crs = B%CICoef(ir) + B%CICoef(is)
   do i=1,NCholesky
      DCholB(i,j) = Crs*B%FF(i,irs)
   enddo
enddo

! monomer A
allocate(ABPLUS1A(A%NDimX,A%NDimX),ABMIN1A(A%NDimX,A%NDimX))
allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1A
read(iunit) ABMIN1A

close(iunit)

! also: we may use a ridiculous alternative and generate full matrices...

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkA,A0BlkIVA,A%IndN,A%CICoef,nblkA,NBas,A%NDimX,'XY0_A')

! AB1 = AB1 - A0
call subtr_blk_right(ABPLUS1A,A0BlkA,A0BlkIVA,.false.,nblkA,A%NDimX)
call subtr_blk_right(ABMIN1A, A0BlkA,A0BlkIVA,.true., nblkA,A%NDimX)

print*, 'ABPLUS1A',norm2(ABPLUS1A)
print*, 'ABMIN1A ',norm2(ABMIN1A)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1A, A1A,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1A,A1A,1.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'X')
print*, 'A1A',norm2(A1A)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A2A,A%NDimX)
deallocate(ABMIN1A)
print*, 'A2A',norm2(A2A)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(DCholA),ABP0TildeA,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde',norm2(ABP0TildeA)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeA(A%NDimX,NCholesky))
Call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS1A,A%NDimX,DCholA,NCholesky,0.0d0,ABP1TildeA,A%NDimX)
print*, 'APLUS1Tilde',norm2(ABP1TildeA)

deallocate(ABPLUS1A)
call release_ac0block(A0BlkA,A0BlkIVA,nblkA)
deallocate(A0BlkA)

! monomer B
allocate(ABPLUS1B(B%NDimX,B%NDimX),ABMIN1B(B%NDimX,B%NDimX))
allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))

!allocate(ABPMB(B%NDimX,B%NDimX),ABPTildeB(B%NDimX,NCholesky))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1B
read(iunit) ABMIN1B

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkB,A0BlkIVB,B%IndN,B%CICoef,nblkB,NBas,B%NDimX,'XY0_B')

! AB1 = AB1 - A0
call subtr_blk_right(ABPLUS1B,A0BlkB,A0BlkIVB,.false.,nblkB,B%NDimX)
call subtr_blk_right(ABMIN1B, A0BlkB,A0BlkIVB,.true., nblkB,B%NDimX)

print*, 'ABPLUS1B',norm2(ABPLUS1B)
print*, 'ABMIN1B ',norm2(ABMIN1B)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1B, A1B,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1B,A1B,1.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'X')
print*, 'A1B',norm2(A1B)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A2B,B%NDimX)
deallocate(ABMIN1B)
print*, 'A2B',norm2(A2B)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeB(B%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(DCholB),ABP0TildeB,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeB(B%NDimX,NCholesky))
Call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS1B,B%NDimX,DCholB,NCholesky,0.0d0,ABP1TildeB,B%NDimX)
print*, 'APLUS1Tilde-b',norm2(ABP1TildeB)

deallocate(ABPLUS1B)
call release_ac0block(A0BlkB,A0BlkIVB,nblkB)
deallocate(A0BlkB)

! get CAlphaTilde in an iterative manner

NFreq = 12
Max_Cn = 5
print*, 'SAPT%NFreq =', NFreq
print*, 'SAPT%MaxCn =', Max_Cn

allocate(XFreq(NFreq),WFreq(NFreq))

allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(C1TildeA(A%NDimX,NCholesky),C1TildeB(B%NDimX,NCholesky))
allocate(C2TildeA(A%NDimX,NCholesky),C2TildeB(B%NDimX,NCholesky))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(WorkA(A%NDimX,NCholesky),WorkB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)
   !OmI = 1.08185673347417

!  Calc: LAMBDA=(A0+Om^2)^-1
   call calculateLambda(LambdaA,LambdaIVA,OmI**2,A%NDimX,nblkA,A0BlkA,A0BlkIVA)
   call calculateLambda(LambdaB,LambdaIVB,OmI**2,B%NDimX,nblkB,A0BlkB,A0BlkIVB)

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
   call ABPM_HALFTRAN_GEN_L(ABP0TildeA,C0TildeA,0.0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')
   call ABPM_HALFTRAN_GEN_L(ABP0TildeB,C0TildeB,0.0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')
   C0TildeA = 0.5d0*C0TildeA
   C0TildeB = 0.5d0*C0TildeB

   print*, 'C0Tilde-A',OmI,norm2(C0TildeA)
   print*, 'C0Tilde-B',OmI,norm2(C0TildeB)

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,A1A,A%NDimX,C0TildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   CTildeA = 0.5d0*ABP1TildeA - CTildeA
   call ABPM_HALFTRAN_GEN_L(CTildeA,C1TildeA,0.0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')

   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,A1B,B%NDimX,C0TildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   CTildeB = 0.5d0*ABP1TildeB - CTildeB
   call ABPM_HALFTRAN_GEN_L(CTildeB,C1TildeB,0.0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')

   print*, 'C1Tilde-A',OmI,norm2(C1TildeA)
   print*, 'C1Tilde-B',OmI,norm2(C1TildeB)

   ! test uncoupled
   CTildeA = 0
   CTildeB = 0

   CTildeA = C0TildeA
   CTildeB = C0TildeB

   !! test semicoupled
   CTildeA = CTildeA + C1TildeA
   CTildeB = CTildeB + C1TildeB

   XFactorial = 1
   do N=2,Max_Cn

       XFactorial = XFactorial*N
       XN1 = -N
       XN2 = -N*(N-1)

       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN2,A2A,A%NDimX,C0TildeA,A%NDimX,0.0d0,WorkA,A%NDimX)
       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN1,A1A,A%NDimX,C1TildeA,A%NDimX,1.0d0,WorkA,A%NDimX)
       call ABPM_HALFTRAN_GEN_L(WorkA,C2TildeA,0.0d0,LambdaA,LambdaIVA,nblkA,A%NDimX,NCholesky,'X')
       !call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.0d0,LAMBDA,NDimX,WORK0,NDimX,0.0d0,C2TildeA,NDimX)

       CTildeA  = CTildeA + C2TildeA / XFactorial
       C0TildeA = C1TildeA
       C1TildeA = C2TildeA

       print*, 'N,COMTildeA',N,norm2(CTildeA)

       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN2,A2B,B%NDimX,C0TildeB,B%NDimX,0.0d0,WorkB,B%NDimX)
       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN1,A1B,B%NDimX,C1TildeB,B%NDimX,1.0d0,WorkB,B%NDimX)
       call ABPM_HALFTRAN_GEN_L(WorkB,C2TildeB,0.0d0,LambdaB,LambdaIVB,nblkB,B%NDimX,NCholesky,'X')
       !call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.0d0,LAMBDA,NDimX,WORK0,NDimX,0.0d0,C2Tilde,NDimX)

       CTildeB = CTildeB + C2TildeB / XFactorial
       C0TildeB = C1TildeB
       C1TildeB = C2TildeB
       
       print*, 'N,COMTildeB',N,norm2(CTildeB)

   enddo

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,DCholB,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)
!
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

!SAPT%e2disp  = -8d0/Pi*e2d

e2d = -32d0/Pi*e2d*1d3
print*, 'e2d = ',e2d
!call print_en('E2disp(Cmat)',e2d,.false.)

deallocate(WFreq,XFreq)
deallocate(A2A,A1A)
deallocate(ABP1TildeA,ABP0TildeA)
deallocate(CTildeB,CTildeA)
deallocate(C2TildeB,C2TildeA)
deallocate(C1TildeB,C1TildeA)
deallocate(C0TildeB,C0TildeA)
deallocate(WorkB,WorkA)
deallocate(CB,CA)
deallocate(DCholB,DCholA)

end subroutine e2disp_CAlphaTilde_block

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

end module sapt_CD_pol
