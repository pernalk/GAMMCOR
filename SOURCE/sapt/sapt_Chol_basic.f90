module sapt_Chol_basic

use tran, only : ABPM_HALFTRAN_GEN_L
use sapt_utils
use timing

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
   jj = (j-1)*(B%num0+B%num1)+j
   do i=1,A%num0+A%num1
      ii = (i-1)*(A%num0+A%num1)+i
      eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%OO(:,ii),1,B%OO(:,jj),1)
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
subroutine e2disp_Chol_cpld_batch(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in coupled approximation
! Requires (NCholesky,MaxBatchSize) allocations
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas,NCholesky
integer :: i,j
double precision :: e2d,fact,tmp
!
integer :: iunitA,iunitB
integer :: iloopA,iloopB,nloopA,nloopB
integer :: offA,offB
integer :: BatchSizeA,BatchSizeB
integer :: MaxBatchSizeA,MaxBatchSizeB
!
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: tmpA(:,:),tmpB(:,:)
double precision,allocatable :: tmpAB(:,:)
logical,allocatable          :: condOmA(:),condOmB(:)
!
double precision :: Tcpu,Twall
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

call clock('START',Tcpu,Twall)

NCholesky = SAPT%NCholesky

allocate(OmA(A%NDimX),OmB(B%NDimX))

call Open_RespBatch(A%NDimX,MaxBatchSizeA,OmA,iunitA,'PROB_A')
call Open_RespBatch(B%NDimX,MaxBatchSizeB,OmB,iunitB,'PROB_B')
print*, 'Use MaxBatchSizeA = ', MaxBatchSizeA
print*, 'Use MaxBatchSizeB = ', MaxBatchSizeB

nloopA = (A%NDimX - 1) / MaxBatchSizeA + 1
nloopB = (B%NDimX - 1) / MaxBatchSizeB + 1
print*, 'nloopA = ', nloopA
print*, 'nloopB = ', nloopB

allocate(EVecA(A%NDimX,MaxBatchSizeA),EVecB(B%NDimX,MaxBatchSizeB))
allocate(tmpA(NCholesky,MaxBatchSizeA),tmpB(NCholesky,MaxBatchSizeB))
allocate(tmpAB(MaxBatchSizeA,MaxBatchSizeB))

allocate(condOmA(A%NDimX),condOmB(B%NDimX))
condOmA = (abs(OmA).gt.SmallE.and.abs(OmA).lt.BigE)
condOmB = (abs(OmB).gt.SmallE.and.abs(OmB).lt.BigE)

do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative omega A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative omega B!',i,OmB(i)
enddo

e2d = 0d0

offB = 0
do iloopB=1,nloopB

   BatchSizeB = min(MaxBatchSizeB,B%NDimX-offB)
   call Get_RespBatch(B%NDimX,BatchSizeB,EvecB,iunitB)
   call dgemm('N','N',NCholesky,BatchSizeB,B%NDimX,1d0,B%DChol,NCholesky,EvecB,B%NDimX,0d0,tmpB,NCholesky)

   offA = 0
   do iloopA=1,nloopA

      BatchSizeA = min(MaxBatchSizeA,A%NDimX-offA)
      call Get_RespBatch(A%NDimX,BatchSizeA,EvecA,iunitA)
      call dgemm('N','N',NCholesky,BatchSizeA,A%NDimX,1d0,A%DChol,NCholesky,EvecA,A%NDimX,0d0,tmpA,NCholesky)

      call dgemm('T','N',BatchSizeA,BatchSizeB,NCholesky,1d0,tmpA,NCholesky,tmpB,NCholesky,0d0,tmpAB,MaxBatchSizeA)

      do j=1,BatchSizeB
         if(condOmB(offB+j)) then
            do i=1,BatchSizeA
               if(condOmA(offA+i)) then
                  e2d = e2d + tmpAB(i,j)**2/(OmA(offA+i)+OmB(offB+j))
               endif
            enddo
         endif
      enddo

      offA = offA + BatchSizeA

   enddo
   call Rewind_RespBatch(iunitA)

   offB = offB + BatchSizeB

enddo

call Close_RespBatch(iunitA)
call Close_RespBatch(iunitB)

SAPT%e2disp  = -16d0*e2d
e2d  = -16d0*e2d*1000d0

call print_en('E2disp(batch)',e2d,.true.)

call clock('E2dispCholBatch',Tcpu,Twall)

deallocate(condOmB,condOmA)
deallocate(tmpB,tmpA,tmpAB)
deallocate(OmB,OmA)
deallocate(EvecB,EVecA)

end subroutine e2disp_Chol_cpld_batch
subroutine test_Chol_ints(Flags,A,B,SAPT)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: iunit
integer :: i,j,ij,ic,id,cd,irec
integer :: nA,nB,nC,nD,nAB,nCD
integer :: NCholesky,NBasis
double precision :: diff,diffOne,tot
double precision,allocatable :: work(:),workAO(:)

NBasis    = A%NBasis
NCholesky = SAPT%NCholesky

nA = NBasis
nB = NBasis
nC = NBasis
nD = NBasis

nAB = nA*nB
nCD = nC*nD

allocate(work(nAB),workAO(nAB))

open(newunit=iunit,file='FFFFAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*NBasis**2)

tot = 0d0
irec = 0
do id=1,nD
   do ic=1,nC
      cd = ic+(id-1)*NBasis
      irec = irec + 1
      call dgemv('T',NCholesky,nCD,1d0,A%FF,NCholesky,B%FF(1:NCholesky,cd),1,0d0,work,1)
      read(iunit,rec=irec) workAO(1:nAB)

      ij = 0
      do j=1,nB
      do i=1,nA
         ij = ij + 1
         diffOne = abs(work(ij)) - abs(workAO(ij))
         if(abs(diffOne).gt.1d-7) print*, 'i,j',i,j,diffOne
      enddo
      enddo

      diff = norm2(work) - norm2(workAO)
      tot  = tot + abs(diff)
      write(lout,*) 'cd = ',cd,abs(diff)

   enddo
enddo
print*, 'Total : ', tot
deallocate(work,workAO)
close(iunit)

end subroutine test_Chol_ints

subroutine e2disp_Chol_cpld(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in coupled approximation
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas,NCholesky
integer :: i,j
double precision :: e2d
!
double precision,allocatable :: tmpA(:,:),tmpB(:,:)
double precision,allocatable :: tmpAB(:,:)
double precision,allocatable :: OmA(:), OmB(:)
logical,allocatable          :: condOmA(:),condOmB(:)
!
double precision :: Tcpu,Twall
!
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

call clock('START',Tcpu,Twall)

if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

! print thresholds for discarding spurious omega values
if(SAPT%IPrint>1) then
   write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
   write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
   write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
endif

NCholesky = SAPT%NCholesky

! code below uses 1 (NDimX,NDimX) allocation
! intermediate A
allocate(tmpA(NCholesky,A%NDimX),tmpB(NCholesky,B%NDimX))
allocate(tmpAB(A%NDimX,A%NDimX),OmA(A%NDimX))

call readEvecZ(tmpAB,A%NDimX,'PROP_A')
call readEvalZ(OmA,A%NDimX,'PROP_A')

! I(k,mu) = R(k,pq).Z(pq,mu)
call dgemm('N','N',NCholesky,A%NDimX,A%NDimX,1d0,A%DChol,NCholesky,tmpAB,A%NDimX,0d0,tmpA,NCholesky)
deallocate(tmpAB)

! intermediate B
allocate(tmpAB(B%NDimX,B%NDimX),OmB(B%NDimX))
call readEvecZ(tmpAB,B%NDimX,'PROP_B')
call readEvalZ(OmB,B%NDimX,'PROP_B')

! I(k,mu) = R(k,pq).Z(pq,mu)
call dgemm('N','N',NCholesky,B%NDimX,B%NDimX,1d0,B%DChol,NCholesky,tmpAB,B%NDimX,0d0,tmpB,NCholesky)
deallocate(tmpAB)

! final intermediate
allocate(tmpAB(A%NDimX,B%NDimX))
call dgemm('T','N',A%NDimX,B%NDimX,NCholesky,1d0,tmpA,NCholesky,tmpB,NCholesky,0d0,tmpAB,A%NDimX)

allocate(condOmA(A%NDimX),condOmB(B%NDimX))
condOmA = (abs(OmA).gt.SmallE.and.abs(OmA).lt.BigE)
condOmB = (abs(OmB).gt.SmallE.and.abs(OmB).lt.BigE)

e2d = 0d0
do j=1,B%NDimX
   if(condOmB(j)) then
      do i=1,A%NDimX
         if(condOmA(i)) then
            e2d = e2d + tmpAB(i,j)**2/(OmA(i)+OmB(j))
         endif
      enddo
   endif
enddo
e2d  = -16d0*e2d*1000d0
call print_en('E2disp(full)',e2d,.true.)

! write amplitude to a file
call writeampl(tmpAB,'PROP_AB')

deallocate(tmpAB,tmpB,tmpA)
deallocate(OmB,OmA)
deallocate(condOmB,condOmA)

end subroutine e2disp_Chol_cpld

subroutine e2disp_Chol_unc(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in uncoupled approximation
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: NCholesky
integer :: i,j,ik,pq,rs
logical,allocatable          :: condOmA(:),condOmB(:)
double precision,allocatable :: OmA0(:), OmB0(:)
double precision,allocatable :: tmpA(:,:),tmpB(:,:)
double precision,allocatable :: tmpAB(:,:)
double precision :: e2du

double precision :: Tcpu,Twall
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

call clock('START',Tcpu,Twall)

! print thresholds for discarding spurious omega values
if(SAPT%IPrint>1) then
   write(LOUT,'(/,1x,a)') 'Thresholds in E2disp(unc):'
   write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
   write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
endif

NCholesky = SAPT%NCholesky

allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
allocate(OmA0(A%NDimX),OmB0(B%NDimX))

allocate(tmpA(NCholesky,A%NDimX))
call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')

! I(k,mu) = R(k,pq).Z(pq,mu)
tmpA = 0d0
do pq=1,A%NDimX
   associate(Y => Y01BlockA(pq))
      do ik=1,NCholesky
         tmpA(ik,Y%l1:Y%l2) = tmpA(ik,Y%l1:Y%l2) + A%DChol(ik,pq)*Y%vec0(1:Y%n)
      enddo
   end associate
enddo

allocate(tmpB(NCholesky,B%NDimX))
call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')
tmpB = 0d0
do rs=1,B%NDimX
   associate(Y => Y01BlockB(rs))
      do ik=1,NCholesky
         tmpB(ik,Y%l1:Y%l2) = tmpB(ik,Y%l1:Y%l2) + B%DChol(ik,rs)*Y%vec0(1:Y%n)
      enddo
   end associate
enddo

allocate(tmpAB(A%NDimX,B%NDimX))
call dgemm('T','N',A%NDimX,B%NDimX,NCholesky,1d0,tmpA,NCholesky,tmpB,NCholesky,0d0,tmpAB,A%NDimX)

allocate(condOmA(A%NDimX),condOmB(B%NDimX))
condOmA = (abs(OmA0).gt.SmallE.and.abs(OmA0).lt.BigE)
condOmB = (abs(OmB0).gt.SmallE.and.abs(OmB0).lt.BigE)

e2du = 0d0
do j=1,B%NDimX
   if(condOmB(j)) then
      do i=1,A%NDimX
         if(condOmA(i)) then
            e2du = e2du + tmpAB(i,j)**2/(OmA0(i)+OmB0(j))
         endif
      enddo
   endif
enddo
SAPT%e2disp_unc = -16d0*e2du

e2du = -16d0*e2du*1000d0

call print_en('E2disp(unc)',e2du,.true.)

deallocate(OmB0,OmA0)
deallocate(tmpAB,tmpB,tmpA)

end subroutine e2disp_Chol_unc


end module sapt_Chol_basic
