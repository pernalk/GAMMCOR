module tran
!!! CAREFUL :: unsym procedures yet to be tested!!!
use types,only: LOUT
implicit none


contains

subroutine tran4_unsym(NBas,nA,CA,nB,CB,nC,CC,nD,CD,TNO)
! CAREFUL: CA,..,CD should be passed
! in AOMO form!!!!
implicit none

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
double precision,intent(out) :: TNO(:,:)
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit
integer :: r,s,rs,d

! 4-index transformation to NOs
! integrals stored in TNO 
 
 write(6,'()') 
 write(6,'(1x,a)') 'Transforming integrals for AB dimer'

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(nA*nB,nC))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 TNO = 0
 do s=1,NBas
    work3 = 0

    do r=1,NBas
       rs = min(r,s) + max(r,s)*(max(r,s)-1)/2 ! idx_tiang
       read(iunit,rec=rs) work1(1:NBas*(NBas+1)/2)
       call triang_to_sq(work1,work2,NBas)
       ! work1=CA^T.work2
       call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
       ! work2=work1.CB
       call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
       ! ... 
       call dger(nA*nB,nC,1d0,work2,1,CC(r:nC*NBas),NBas,work3,nA*nB)
   
    enddo
    
    !call dger(nA*nB*nC,nD,1d0,work3,1,CD(s),NBas,TNO,nA*nB*nC)
    do d=1,nD
       TNO(:,(d-1)*nC+1:d*nC) = TNO(:,(d-1)*nC+1:d*nC) + work3*CD(s+(d-1)*NBas)
    enddo

 enddo

 deallocate(work1,work2,work3)
 close(iunit)

end subroutine tran4_unsym

subroutine tran4_unsym2(NBas,nA,CA,nB,CB,nC,CC,nD,CD,TNO)
implicit none

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
double precision,intent(out) :: TNO(:,:)
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit
integer :: r,s,rs

! 4-index transformation to NOs
! integrals stored in TNO 
 
 write(6,'()') 
 write(6,'(1x,a)') 'Transforming integrals for AB dimer'

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(nA*nB,nC))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 TNO = 0
 do s=1,NBas
    work3 = 0

    do r=1,NBas
       rs = min(r,s) + max(r,s)*(max(r,s)-1)/2 ! idx_tiang
       read(iunit,rec=rs) work1(1:NBas*(NBas+1)/2)
       call triang_to_sq(work1,work2,NBas)
       ! work1=CA.work2
       call dgemm('N','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
       ! work2=work1.CB^T
       call dgemm('N','T',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
       ! ... 
       call dger(nA*nB,nC,1d0,work2,1,CC((r-1)*NBas+1),1,work3,nA*nB)
   
    enddo
    
    call dger(nA*nB*nC,nD,1d0,work3,1,CD((s-1)*NBas+1),1,TNO,nA*nB*nC)
!    do d=1,nD
!       TNO(:,(d-1)*nC+1:d*nC) = TNO(:,(d-1)*nC+1:d*nC) + work3*CD(s+(d-1)*NBas)
!    enddo

 enddo

 deallocate(work1,work2,work3)
 close(iunit)

end subroutine tran4_unsym2

subroutine tran3Q_full(NBas,CX,Q,fname)
! 3-index transformation out of core
implicit none

integer,intent(in) :: NBas
!integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CX(*),Q(*)
character(*) :: fname
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: ntr,nsq,nloop
integer,parameter :: cbuf=512
integer :: i,rs,ab

 write(6,'()') 
! write(6,'(1x,a)') 'TRAN3'
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
 write(6,'(1x,a)') 'Transforming 3-ints for '//fname

 ntr = NBas*(NBas+1)/2
 nsq = NBas**2

 ! set no. of triangles in buffer
 nloop = (ntr-1)/cbuf+1

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,ntr))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (pr|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

       read(iunit,rec=rs) work1(1:ntr)
       call triang_to_sq(work1,work2,NBas)
       ! work1=CA^T.work2
       ! work2=work1.CB
       call dgemm('T','N',NBas,NBas,NBas,1d0,CX,NBas,work2,NBas,0d0,work1,NBas) 
       call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,CX,NBas,0d0,work2,NBas)
       call sq_to_triang(work2,work1,NBas)
       ! transpose
       work3(rs-(i-1)*cbuf,1:ntr) = work1(1:ntr)

    enddo

    do ab=1,ntr
       write(iunit2,rec=(i-1)*ntr+ab) work3(1:cbuf,ab)
    enddo

 enddo

 close(iunit)


! |qa) pr-triang, qa-square
 open(newunit=iunit3,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nsq)

 do ab=1,ntr

    do i=1,nloop
       ! get all (ab| for given |rs)
       read(iunit2,rec=(i-1)*ntr+ab) work1((i-1)*cbuf+1:min(i*cbuf,ntr))
    enddo

    call triang_to_sq(work1,work2,NBas)
    ! work1=CX^T.work2
    ! work2=work1.Q 
    call dgemm('T','N',NBas,NBas,NBas,1d0,CX,NBas,work2,NBas,0d0,work1,NBas)
    call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,Q,NBas,0d0,work2,NBas)

    write(iunit3,rec=ab) work2(1:nsq)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran3Q_full

subroutine tran4_full(NBas,CA,CB,fname)
! 4-index transformation out of core
! dumps all integrals on disk in the (triang,triang) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

integer,intent(in) :: NBas
!integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*)
character(*) :: fname
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: ntr,nloop
integer,parameter :: cbuf=512
integer :: i,rs,ab

 write(6,'()') 
! write(6,'(1x,a)') 'TRAN4_SYM_OUT_OF_CORE'
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
 write(6,'(1x,a)') 'Transforming integrals for '//fname

 ntr = NBas*(NBas+1)/2

 ! set no. of triangles in buffer
 nloop = (ntr-1)/cbuf+1

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,ntr))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (ab|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

       read(iunit,rec=rs) work1(1:ntr)
       call triang_to_sq(work1,work2,NBas)
       ! work1=CA^T.work2
       ! work2=work1.CB
       call dgemm('T','N',NBas,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,NBas) 
       call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,CA,NBas,0d0,work2,NBas)
       call sq_to_triang(work2,work1,NBas)
       ! transpose
       work3(rs-(i-1)*cbuf,1:ntr) = work1(1:ntr)

    enddo

    do ab=1,ntr
       write(iunit2,rec=(i-1)*ntr+ab) work3(1:cbuf,ab)
    enddo

 enddo

 close(iunit)

! |cd)
 open(newunit=iunit3,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 do ab=1,ntr

    do i=1,nloop
       ! get all (ab| for given |rs)
       read(iunit2,rec=(i-1)*ntr+ab) work1((i-1)*cbuf+1:min(i*cbuf,ntr))
    enddo

    call triang_to_sq(work1,work2,NBas)
    ! work1=CC^T.work2
    ! work2=work1.CD
    call dgemm('T','N',NBas,NBas,NBas,1d0,CB,NBas,work2,NBas,0d0,work1,NBas)
    call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,CB,NBas,0d0,work2,NBas)
    call sq_to_triang(work2,work1,NBas)
    write(iunit3,rec=ab) work1(1:ntr)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran4_full

subroutine tran4_gen(NBas,nA,CA,nB,CB,nC,CC,nD,CD,fname)
! 4-index transformation out of core
! dumps all integrals on disk in the (triang,triang) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
character(*) :: fname
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: ntr,nAB,nCD,nloop
integer,parameter :: cbuf=512
integer :: i,rs,ab

 write(6,'()') 
! write(6,'(1x,a)') 'TRAN4_SYM_OUT_OF_CORE'
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
 write(6,'(1x,a)') 'Transforming integrals for '//fname

 ntr = NBas*(NBas+1)/2
 nAB = nA*nB
 nCD = nC*nD

 ! set no. of triangles in buffer
 nloop = (ntr-1)/cbuf+1

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,nAB))


 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (ab|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

       read(iunit,rec=rs) work1(1:ntr)
       call triang_to_sq(work1,work2,NBas)
       ! work1=CA^T.work2
       ! work2=work1.CB
       call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
       call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
       ! transpose
       work3(rs-(i-1)*cbuf,1:nAB) = work2(1:nAB)

    enddo

    do ab=1,nAB
       write(iunit2,rec=(i-1)*nAB+ab) work3(1:cbuf,ab)
    enddo

 enddo

 close(iunit)

! |cd)
 open(newunit=iunit3,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nCD)

 do ab=1,nAB

    do i=1,nloop
       ! get all (ab| for given |rs)
       read(iunit2,rec=(i-1)*nAB+ab) work1((i-1)*cbuf+1:min(i*cbuf,ntr))
    enddo

    call triang_to_sq(work1,work2,NBas)
    ! work1=CC^T.work2
    ! work2=work1.CD
    call dgemm('T','N',nC,NBas,NBas,1d0,CC,NBas,work2,NBas,0d0,work1,nC)
    call dgemm('N','N',nC,nD,NBas,1d0,work1,nC,CD,NBas,0d0,work2,nC)
    write(iunit3,rec=ab) work2(1:nCD)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran4_gen

subroutine make_J1(NBas,X,J)
implicit none
integer :: NBas
double precision :: X(*), J(NBas,NBas)
integer :: iunit, ntr
integer :: ir,is,irs
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
double precision,external :: ddot
! works in DCBS

 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 irs=0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    read(iunit,rec=irs) work1(1:ntr)    
    call triang_to_sq(work1,work2,NBas)

    tmp = ddot(NBas**2,work2,1,X,1)
    J(ir,is) = tmp
    J(is,ir) = tmp

    enddo
 enddo

 deallocate(work1,work2)
 close(iunit)

end subroutine make_J1

subroutine make_J2(NBas,XA,XB,JA,JB)
implicit none
integer :: NBas
double precision :: XA(*), XB(*), JA(NBas,NBas), JB(NBas,NBas)
integer :: iunit, ntr
integer :: ir,is,irs
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
double precision,external :: ddot
! works in DCBS

 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 irs=0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    read(iunit,rec=irs) work1(1:ntr)    
    call triang_to_sq(work1,work2,NBas)

    tmp = ddot(NBas**2,work2,1,XA,1)
    JA(ir,is) = tmp
    JA(is,ir) = tmp

    tmp = ddot(NBas**2,work2,1,XB,1)
    JB(ir,is) = tmp
    JB(is,ir) = tmp

    enddo
 enddo


 deallocate(work1,work2)
 close(iunit)

end subroutine make_J2

subroutine make_K(NBas,X,K)
implicit none
integer,intent(in) :: NBas
double precision,intent(in) :: X(NBas,NBas)
double precision,intent(inout) :: K(NBas,NBas)
integer :: iunit, ntr
integer :: ir,is,irs
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
! works in DCBS

 K = 0
 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*ntr)

 irs = 0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    read(iunit,rec=irs) work1(1:ntr)
    call triang_to_sq(work1,work2,NBas)

    if(ir==is) then

       call dgemv('N',NBas,NBas,1.d0,work2,NBas,X(:,is),1,1.d0,K(:,ir),1)

    else

       call dgemv('N',NBas,NBas,1.d0,work2,NBas,X(:,is),1,1.d0,K(:,ir),1)
       call dgemv('N',NBas,NBas,1.d0,work2,NBas,X(:,ir),1,1.d0,K(:,is),1)

    endif

    enddo
 enddo

 deallocate(work2,work1)
 close(iunit)

end subroutine make_K

subroutine tran_oneint(ints,MO1,MO2,work,nbas)
implicit none

integer,intent(in) :: nbas
double precision,intent(inout),contiguous :: ints(:)
double precision,intent(in),contiguous :: MO1(:),MO2(:)
double precision,contiguous :: work(:)

 if(nbas>0) then
    call dgemm('T','N',nbas,nbas,nbas,1d0,MO1,nbas,ints,nbas,0d0,work,nbas)
    call dgemm('N','N',nbas,nbas,nbas,1d0,work,nbas,MO2,nbas,0d0,ints,nbas)
 endif

end subroutine tran_oneint

subroutine tran2MO(MatIn,Ca,Cb,MatOut,nbas)
implicit none

integer,intent(in) :: nbas
double precision,intent(in) :: MatIn(nbas,nbas)
double precision,intent(in) :: Ca(nbas,nbas),Cb(nbas,nbas)
double precision,intent(out) :: MatOut(nbas,nbas)
double precision,allocatable :: tmp(:,:)
integer :: i,j

 allocate(tmp(nbas,nbas))
 ! tmp=Ca^T.MatIn
 ! MatOut=tmp.Cb
 tmp = 0
 call dgemm('T','N',nbas,nbas,nbas,1d0,Ca,nbas,MatIn,nbas,0d0,tmp,nbas)
 call dgemm('N','N',nbas,nbas,nbas,1d0,tmp,nbas,Cb,nbas,0d0,MatOut,nbas)

 deallocate(tmp)

end subroutine tran2MO

subroutine triang_to_sq(matTr,matSq,NBas)
implicit none

double precision :: matTr(:), matSq(:)
integer :: NBas
integer :: t,p,q,pq,qp

 t = 0
 do q=1,NBas
    do p=1,q
       t = t + 1
       pq = p + (q-1)*NBas
       qp = q + (p-1)*NBas
       matSq(pq) = matTr(t)
       matSq(qp) = matTr(t)
    enddo
 enddo

end subroutine triang_to_sq

subroutine triang_to_sq2(matTr,matSq,NBas)
implicit none

double precision :: matTr(:), matSq(:,:)
integer :: NBas
integer :: t,p,q

 t = 0
 do q=1,NBas
    do p=1,q
       t = t + 1
       matSq(p,q) = matTr(t)
       matSq(q,p) = matTr(t)
    enddo
 enddo

end subroutine triang_to_sq2

subroutine sq_to_triang(matSq,matTr,NBas)
implicit none

double precision :: matSq(:), matTr(:)
integer :: NBas
integer :: t,p,q,pq

 t = 0
 do q=1,NBas
    do p=1,q
       t = t + 1
       pq = p + (q-1)*NBas
       matTr(t) = matSq(pq)
    enddo
 enddo

end subroutine sq_to_triang

subroutine sq_to_triang2(matSq,matTr,NBas)
implicit none

double precision :: matSq(:,:), matTr(:)
integer :: NBas
integer :: t,p,q

 t = 0
 do q=1,NBas
    do p=1,q
       t = t + 1
       matTr(t) = matSq(p,q)
    enddo
 enddo

end subroutine sq_to_triang2


subroutine sq_symmetrize(mat,NBas)
implicit none

double precision :: mat(NBas,NBas)
integer :: NBas
integer :: i,j
double precision :: val

do j=2,NBas
   do i=1,j-1
      val = (mat(i,j)+mat(j,i))*0.5d0
      mat(i,j) = val
      mat(j,i) = val
   enddo
enddo

end subroutine sq_symmetrize

subroutine AB_CAS_mithap(ABPLUS,ABMIN,ETot,URe,Occ,XOne,IPair, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDim,NInte1,IntFileName,ACAlpha)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
! 
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! THE FOLLOWING SYMMETRY IS ASSUMED
! RDM2(ij,kl) = RDM2(kl,ij)
! ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
! SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntFileName
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IPair(NBasis,NBasis),IndN(2,NDim),IndX(NDim),IGem(NBasis)
double precision,intent(in)  :: ACAlpha

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: NOccup,NRDM2Act
integer :: Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
ETot   = 0

! set dimensions
NOccup = NAct + INActive
Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)


call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = ACAlpha
            endif
         enddo
      enddo
   enddo
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxI  = 0
AuxIO = 0
WMAT  = 0

open(newunit=iunit,file=trim(IntFileName),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do ll=1,NBasis
   do kk=1,ll
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
      do klround=1,merge(1,2,kk==ll)
         select case(klround)
         case(1)
            k = kk
            l = ll
         case(2)
            k = ll
            l = kk
         end select

! COMPUTE THE ENERGY FOR CHECKING
         if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))

! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
         ! Coulomb
         if(IGem(k)==IGem(l)) then

            val = 0
            do i=1,NBasis
               if(IGem(i)/=IGem(k)) val = val + Occ(i)*ints(i,i)
            enddo
            val = 2*(1-ACAlpha)*val
            HNO(k,l) = HNO(k,l) + val

         ! exchange
         else

            val = (1-ACAlpha)*Occ(k)
            do i=1,NBasis
               if(IGem(i)==IGem(l)) HNO(i,l) = HNO(i,l) - val*ints(i,k)
            enddo

         endif

! AUXILIARY MATRIX AuxI AND AuxIO
         ! Coulomb
         val = 0
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),1,1)
         do i=1,INActive
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxIO(k,l) = AuxIO(k,l) + 2*val
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),2,2)
         do i=INActive+1,NOccup
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxI(k,l) = AuxI(k,l) + 2*val

         ! exchange
         if(k<=INActive) then
            do i=1,NBasis
               AuxIO(i,l) = AuxIO(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
            enddo
         endif
         if(k<=NOccup) then
            do i=1,NBasis
               AuxI(i,l) = AuxI(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
            enddo
         endif

! AUXILIARY MATRIX WMAT
         if(l<=NOccup) then
            do ir=1,NOccup
               val = 0
               val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                    sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                    sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                    sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
               val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                    sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
               WMAT(k,ir) = WMAT(k,ir) + val
            enddo
         endif

! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
         ! T1+T2+T3+T4+T5+T6: AuxInd = 1
         ! Coulomb
         if(k>INActive.and.l<=NOccup) then
            ir = k
            is = l
            irs = pos(ir,is)
            if(irs>0) then

               do iq=1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = 2*AuxVal*val*ints(ip,iq)
                  
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = 2*AuxVal*val*ints(iq,ip)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo

            endif
         endif

         ! exchange
         if(l<=NOccup) then
            if(k<=NOccup) then
               iq = k
               is = l

               do ir=INActive+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then

                           AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                           val = 0
                           if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                           if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                           if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                           if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                           val = - AuxVal*val*ints(ip,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  endif
               enddo
               
            endif
            if(k>INActive) then
               ip = k
               is = l

               do ir=INActive+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then                  
                     do iq=1,NOccup
                        ipq = pos(ip,iq)
                        if(ipq>0) then

                           AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                           val = 0
                           if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                           if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                           if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                           if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                           val = - AuxVal*val*ints(iq,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  endif
               enddo
            endif
         endif

         ! T1 
         if((k<=NOccup).and.(l<=NOccup)) then
            iq = k
            is = l

            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               endif
            enddo
            
         endif
         
         ! T1 p->q
         if((k>INActive).and.(l<=NOccup)) then
            ip = k
            is = l

            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        
                        val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
                        
                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               endif
            enddo
             
         endif

         ! T3
         if((k>INActive).and.(l>INActive)) then
            ip = k
            ir = l

            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
            
         endif

         if((k<=NOccup).and.(l>INActive)) then
            iq = k
            ir = l

            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                        
                     endif
                  enddo
               endif
            enddo
            
         endif
         
         ! T2+T5: AuxInd = 2
         if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
            is = k
            iu = l
            
            ! T2
            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     do ip=INActive+1,NOccup
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                        endif
                     enddo
                  enddo
               endif
            enddo
            ! p->q
            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo

            ! T5
            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
            ! p->q
            do ir=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     do ip=INActive+1,NOccup
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                        endif
                     enddo
                  enddo
               endif
            enddo
         
      endif
            
         ! T4+T6: AuxInd = 2
         if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
            ir = k
            iu = l

            ! T4
            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(INActive+1:NOccup,ip))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
            ! p->q
            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     do ip=INActive+1,NOccup
                        ipq = pos(ip,iq)
                        if(ipq>0) then
                           val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
            
            ! T6
            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     do ip=INActive+1,NOccup
                        ipq = pos(ip,iq)
                        if(ipq>0) then                        
                           val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
                           
                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo
            ! p->q
            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then  
                           val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                                sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(INActive+1:NOccup,ip))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo

         endif
         
      enddo
   enddo
enddo

close(iunit)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = IndX(ICol)
   if(ir/=is) then
      
      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = IndX(IRow)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
      enddo

   endif
enddo

print*, "AB-my",norm2(ABPLUS),norm2(ABMIN)

call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'DUPAmy', norm2(work1(1:NBasis*(NBasis+1)/2))
HNO=transpose(HNO)
call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'DUPCmy', norm2(work1(1:NBasis*(NBasis+1)/2))

write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot

call sq_to_triang2(AuxI,work1,NBasis)
write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
AuxI=transpose(AuxI)
call sq_to_triang2(AuxI,work1,NBasis)
write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

call sq_to_triang2(AuxIO,work1,NBasis)
write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
AuxIO=transpose(AuxIO)
call sq_to_triang2(AuxIO,work1,NBasis)
write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_mithap

end module


