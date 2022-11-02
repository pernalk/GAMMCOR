subroutine chol_ints_fofo(nA,nB,MatAB,nC,nD,MatCD,NCholesky,NBas,fname)
!
! assumes that MatAB(CD) are NChol,FF type
! constructs (FF|OO) or (FO|FO) integrals
!
implicit none

integer,intent(in) :: nA,nB,nC,nD
integer,intent(in) :: NBas,NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatAB(NCholesky,NBas**2), &
                               MatCD(NCholesky,NBas**2)

integer :: iunit
integer :: nAB,nCD,cd
integer :: ic,id,irec
double precision,allocatable :: work(:)

nAB = nA*nB
nCD = nC*nD

allocate(work(nAB))

print*, 'Assemble ',fname,' from Cholesky Vectors'

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nAB)

! (FF|OO)
irec = 0
do id=1,nD
   do ic=1,nC
      irec = irec + 1
      cd = ic+(id-1)*NBas
      call dgemv('T',NCholesky,nAB,1d0,MatAB,NCholesky,MatCD(1:NCholesky,cd),1,0d0,work,1)
      write(iunit,rec=irec) work(1:nAB)

   enddo
enddo

deallocate(work)
close(iunit)

end subroutine chol_ints_fofo

subroutine chol_fofo_batch(nOA,MatFA,nOB,MatFB,NCholesky,NBasis,fname)
!
! (FO|FO) | (FA|FB)
!
implicit none

integer,intent(in) :: nOA,nOB
integer,intent(in) :: NBasis,NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatFA(NCholesky,NBasis**2), &
                               MatFB(NCholesky,NBasis**2)

integer :: iunit
integer :: nAB,nCD,cd
integer :: i,j,k,l,ic,id,irec
integer :: iloop,nloop,off
integer :: iBatch
integer :: BatchSize,MaxBatchSize = 33
double precision,allocatable :: work(:,:),ints(:,:)

nAB = NBasis*nOA
nCD = NBasis*nOB

print*, 'Assemble ',fname,' from Cholesky Vectors in Batches'

allocate(work(nAB,MaxBatchSize),ints(NBasis,NBasis))

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nAB)

nloop = (nCD - 1) / MaxBatchSize + 1

off = 0
k   = 0
l   = 1
do iloop=1,nloop

   ! batch size for each iloop; last one is smaller
   BatchSize = min(MaxBatchSize,nCD-off)

   ! assemble (FO|BatchSize) batch from CholVecs
   call dgemm('T','N',nAB,BatchSize,NCholesky,1d0,MatFA,NCholesky, &
              MatFB(:,off+1:BatchSize),NCholesky,0d0,work,nAB)

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NBasis) then
         k = 1
         l = l + 1
      endif

      do j=1,nOB
         do i=1,NBasis
            ints(i,j) = work((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(l>nOB) cycle
      ints(:,nOB+1:NBasis) = 0

      irec = (l-1)*NBasis + k
      write(iunit,rec=irec) ints(1:NBasis,1:nOB)
   enddo

   off = off + MaxBatchSize

enddo

close(iunit)

deallocate(work,ints)

end subroutine chol_fofo_batch

subroutine chol_ffoo_batch(nOA,MatFA,nOB,MatFB,NCholesky,NBasis,fname)
!
! (FF|OO) | (FF|AB)
!
implicit none

integer,intent(in) :: nOA,nOB
integer,intent(in) :: NBasis,NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatFA(NCholesky,NBasis**2), &
                               MatFB(NCholesky,NBasis**2)

integer :: iunit
integer :: nAB,nFF,cd
integer :: i,j,k,l,kk,ll,kl,irec
integer :: iloop,nloop,off
integer :: iBatch
integer :: BatchSize,MaxBatchSize = 33
double precision,allocatable :: work1(:,:),work2(:,:),ints(:,:)

nFF = NBasis*NBasis
nAB = nOA*nOB

nloop = (nAB - 1) / MaxBatchSize + 1

print*, 'Assemble ',fname,' from Cholesky Vectors in Batches'

allocate(work1(NBasis**2,MaxBatchSize),work2(NCholesky,MaxBatchSize))
allocate(ints(NBasis,NBasis))

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nFF)

off = 0
k   = 0
l   = 1
do iloop=1,nloop
   ! batch size
   BatchSize = min(MaxBatchSize,nAB-off)

   kk = k
   ll = l
   do iBatch=1,BatchSize
      kk = kk + 1
      if(kk>nOB) then
         kk = 1
         ll = ll + 1
      endif
      kl = (ll - 1)*NBasis + kk
      work2(:,iBatch) = MatFB(:,kl)
   enddo

   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,1d0,MatFA,NCholesky, &
              work2,NCholesky,0d0,work1,NBasis**2)

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>nOB) then
         k = 1
         l = l + 1
      endif

      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(k>nOB.or.l>nOB) cycle
      irec = (l - 1)*nOB + k
      write(iunit,rec=irec) ints(1:NBasis,1:NBasis)

   enddo

   off = off + MaxBatchSize

enddo

close(iunit)
deallocate(work2,work1,ints)

end subroutine chol_ffoo_batch

subroutine chol_triang_fofo(nA,nB,MatAB,nC,nD,MatCD,NCholesky,NInte1,NBas,fname)
!
! MatAB and MatCD are (NChol,FFtriang) type
! constructs (FF|OO) or (FO|FO) integrals
! it is the Cholesky counterpart of read4_gen in tran
! this was added to handle ints from Libor
!
implicit none

integer,intent(in) :: nA,nB,nC,nD
integer,intent(in) :: NInte1,NBas,NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatAB(1:NCholesky,1:NInte1), &
                               MatCD(1:NCholesky,1:NInte1)

integer :: iunit
integer :: nAB,nCD,cd
integer :: ic,id,irec
integer :: ia,ib,iab,tab
double precision,allocatable :: work(:),workT(:)
double precision,allocatable :: ints(:,:)

nAB = nA*nB
nCD = nC*nD

allocate(workT(NInte1),work(nAB),ints(NBas,NBas))

print*, 'Assemble ',fname,' from Cholesky Vectors'

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nAB)

! (FF|OO)
irec = 0
do id=1,nD
   do ic=1,nC
      irec = irec + 1
      cd = max(ic,id)*(max(ic,id)-1)/2 + min(ic,id)
      call dgemv('T',NCholesky,NInte1,1d0,MatAB(1:NCholesky,:),NCholesky,MatCD(1:NCholesky,cd),1,0d0,workT,1)

      ! unpack triangle to square
      iab  = 0
      do ib=1,nB
         do ia=1,nA
            iab = iab + 1
            tab = max(ia,ib)*(max(ia,ib)-1)/2 + min(ia,ib)
            work(iab) = workT(tab)
         enddo
      enddo

      write(iunit,rec=irec) work(1:nAB)

   enddo
enddo

end subroutine chol_triang_fofo

subroutine chol_ints_oooo(nA,nB,MatAB,nC,nD,MatCD,NCholesky,fname)
!
! assumes that MatAB(CD) are NChol,OO type
! constructs (OO|OO) integrals
!
implicit none

integer,intent(in)          :: nA,nB,nC,nD
integer,intent(in)          :: NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatAB(NCholesky,nA*nB), &
                               MatCD(NCholesky,nC*nD)

integer :: iunit
integer :: nAB,nCD,cd
integer :: ic,id,irec
double precision,allocatable :: work(:)

nAB = nA*nB
nCD = nC*nD

allocate(work(nAB))

print*, 'Assemble ',fname,' from Cholesky Vectors'

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nAB)
! (OO|OO)
irec = 0
do id=1,nD
   do ic=1,nC
      irec = irec + 1
      cd = ic+(id-1)*nC
      call dgemv('T',NCholesky,nAB,1d0,MatAB,NCholesky,MatCD(1:NCholesky,cd),1,0d0,work,1)
      write(iunit,rec=irec) work(1:nAB)
   enddo
enddo

deallocate(work)
close(iunit)

end subroutine chol_ints_oooo

subroutine FockGen_Chol(Fock,Vecs,OneRdm,XOne,NInte1,NCholesky,NBasis)
!
!     GENERALIZED FOCK MATRIX
!     COMPOSED FROM CHOLESKY VECS
!
use tran
implicit none

integer,intent(in) :: NInte1,NCholesky,NBasis
double precision,intent(in)  :: Vecs(1:NCholesky,1:NInte1),OneRdm(NInte1),XOne(NInte1)
double precision,intent(out) :: Fock(Ninte1)

integer :: iunit,kl,k,l
double precision,allocatable :: OneRdmSq(:,:),FockSq(:,:),ints(:,:),work1(:)

print*, 'Assemble Fock from CholeskyVecs...'

allocate(OneRdmSq(NBasis,NBasis),FockSq(NBasis,NBasis),ints(NBasis,NBasis),work1(NInte1))

call triang_to_sq2(OneRdm,OneRdmSq,NBasis)
call triang_to_sq2(XOne,FockSq,NBasis)

kl = 0
ints=0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      call dgemv('T',NCholesky,NInte1,1d0,Vecs(1:NCholesky,:),NCholesky,Vecs(1:NCholesky,kl),1,0d0,work1,1)
      call triang_to_sq2(work1,ints,NBasis)

      if(k==l) then
        call daxpy(NBasis**2,2.d0*OneRdmSq(k,l),ints,1,FockSq,1)
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,k),1,1.d0,FockSq(:,l),1)
      else
        call daxpy(NBasis**2,2.d0*(OneRdmSq(k,l)+OneRdmSq(l,k)),ints,1,FockSq,1)
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,k),1,1.d0,FockSq(:,l),1)
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,l),1,1.d0,FockSq(:,k),1)
      endif

   enddo
enddo

call sq_to_triang2(FockSq,Fock,NBasis)

deallocate(work1,ints,FockSq,OneRdmSq)

end subroutine FockGen_Chol

subroutine FockGen_CholR(Fock,Vecs,OneRdm,XOne,NInte1,NCholesky,NBasis)
!
!     GENERALIZED FOCK MATRIX
!     COMPOSED FROM CHOLESKY VECS
!
use tran
implicit none

!type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: NInte1,NCholesky,NBasis
double precision,intent(in)  :: Vecs(NCholesky,NInte1),OneRdm(NInte1),XOne(NInte1)
double precision,intent(out) :: Fock(Ninte1)

integer :: k
double precision :: val
double precision,allocatable :: OneRdmSq(:,:),FockSq(:,:),ints(:),work(:)
double precision,external :: ddot

!NCholesky  = CholeskyVecs%NCholesky

allocate(OneRdmSq(NBasis,NBasis),FockSq(NBasis,NBasis),ints(NBasis**2),work(NBasis**2))

call triang_to_sq2(OneRdm,OneRdmSq,NBasis)
call triang_to_sq2(XOne,FockSq,NBasis)

ints = 0
do k=1,NCholesky
    call triang_to_sq(Vecs(k,1:NInte1),ints,NBasis)
    val = ddot(NBasis**2,ints,1,OneRdmSq,1)
    call daxpy(NBasis**2,2d0*val,ints,1,FockSq,1)
    call dgemm('N','N',NBasis,NBasis,NBasis,1d0,ints,NBasis,OneRdmSq,NBasis,0d0,work,NBasis)
    call dgemm('N','N',NBasis,NBasis,NBasis,-1d0,work,NBasis,ints,NBasis,1d0,FockSq,NBasis)
enddo

call sq_to_triang2(FockSq,Fock,NBasis)

deallocate(work,ints,FockSq,OneRdmSq)

end subroutine FockGen_CholR

subroutine make_K_Chol(Vecs,NCholesky,NBasis,X,K)
! prepare K_sp = X_rq (pq|rs)
! this is for test purpose only (v. slow)

use tran
implicit none

integer,intent(in) :: NCholesky,NBasis
double precision,intent(in)    :: Vecs(NCholesky,NBasis),X(NBasis,NBasis)
double precision,intent(inout) :: K(NBasis,NBasis)

integer :: iunit, ntr
integer :: ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)

 print*, 'Assemble Kmat from CholeskyVecs...'

 K = 0
 ntr = NBasis*(NBasis+1)/2

 allocate(work1(NBasis*NBasis),work2(NBasis*NBasis))

 irs = 0
 do is=1,NBasis
    do ir=1,is
    irs = irs + 1
    call dgemv('T',NCholesky,ntr,1d0,Vecs(1:NCholesky,:),NCholesky,Vecs(1:NCholesky,irs),1,0d0,work1,1)

    call triang_to_sq(work1,work2,NBasis)

    if(ir==is) then

       call dgemv('N',NBasis,NBasis,1.d0,work2,NBasis,X(:,is),1,1.d0,K(:,ir),1)

    else

       call dgemv('N',NBasis,NBasis,1.d0,work2,NBasis,X(:,is),1,1.d0,K(:,ir),1)
       call dgemv('N',NBasis,NBasis,1.d0,work2,NBasis,X(:,ir),1,1.d0,K(:,is),1)

    endif

    enddo
 enddo

 deallocate(work2,work1)

end subroutine make_K_Chol

subroutine make_K_CholR(Vecs,NCholesky,NBasis,X,K)
! prepare K_sp = X_rq (pq|rs)

use tran
implicit none

integer,intent(in) :: NCholesky,NBasis
double precision,intent(in)    :: Vecs(1:NCholesky,NBasis*(NBasis+1)/2), X(NBasis,NBasis)
double precision,intent(inout) :: K(NBasis,NBasis)

integer :: iunit, ntr
integer :: ik,ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work(:),workC(:)

 print*, 'Assemble Kmat from CholeskyVecs...'

 K = 0
 ntr = NBasis*(NBasis+1)/2

 allocate(work(NBasis*NBasis),workC(NBasis*NBasis))

 workC = 0
 do ik=1,NCholesky
    call triang_to_sq(Vecs(ik,1:ntr),work,NBasis)
    call dgemm('N','N',NBasis,NBasis,NBasis,1d0,work,NBasis,X,NBasis,0d0,workC,NBasis)
    call dgemm('N','N',NBasis,NBasis,NBasis,1d0,workC,NBasis,work,NBasis,1d0,K,NBasis)
 enddo

 deallocate(workC,work)

end subroutine make_K_CholR

subroutine make_J2_Chol(Vecs,XA,XB,JA,JB,NCholesky,NBasis)
! prepate Coulomb matrices for monomers A and B
! J_rs = X_pq (pq|rs)
! this is for test purpose only (v. slow)

use tran
implicit none

integer,intent(in) :: NCholesky,NBasis
double precision,intent(in)  :: Vecs(NCholesky,NBasis*(NBasis+1)/2), XA(*), XB(*)
double precision,intent(out) :: JA(NBasis,NBasis), JB(NBasis,NBasis)

integer :: iunit, ntr
integer :: ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
double precision,external :: ddot

 print*, 'Assemble Jmat from CholeskyVecs...'

 JA = 0
 JB = 0
 ntr = NBasis*(NBasis+1)/2

 allocate(work1(NBasis*NBasis),work2(NBasis*NBasis))

 irs=0
 do is=1,NBasis
    do ir=1,is
    irs = irs + 1
    call dgemv('T',NCholesky,ntr,1d0,Vecs(1:NCholesky,:),NCholesky,Vecs(1:NCholesky,irs),1,0d0,work1,1)

    call triang_to_sq(work1,work2,NBasis)

    tmp = ddot(NBasis**2,work2,1,XA,1)
    JA(ir,is) = tmp
    JA(is,ir) = tmp

    tmp = ddot(NBasis**2,work2,1,XB,1)
    JB(ir,is) = tmp
    JB(is,ir) = tmp

    enddo
 enddo

 deallocate(work1,work2)

end subroutine make_J2_Chol

subroutine make_J2_CholR(Vecs,XA,XB,JA,JB,NCholesky,NBasis)
! prepate Coulomb matrices for monomers A and B
! J_rs = X_pq (pq|rs)
! works in DCBS

use tran
implicit none

integer,intent(in) :: NCholesky,NBasis
double precision,intent(in)  :: Vecs(NCholesky,NBasis*(NBasis+1)/2), XA(*), XB(*)
double precision,intent(out) :: JA(NBasis,NBasis), JB(NBasis,NBasis)

integer :: ntr
integer :: k
double precision :: valA,valB
double precision,allocatable :: work(:)
double precision,external :: ddot

 print*, 'Assemble Jmat from CholeskyVecs...'

 JA = 0
 JB = 0

 ntr = NBasis*(NBasis+1)/2

 allocate(work(NBasis**2))

 do k=1,NCholesky
    call triang_to_sq(Vecs(k,1:ntr),work,NBasis)
    valA = ddot(NBasis**2,work,1,XA,1)
    valB = ddot(NBasis**2,work,1,XB,1)
    call daxpy(NBasis**2,valA,work,1,JA,1)
    call daxpy(NBasis**2,valB,work,1,JB,1)
 enddo

 deallocate(work)

end subroutine make_J2_CholR

