module tran
!!! CAREFUL :: unsym procedures yet to be tested!!!
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

subroutine tran4_sym(NBas,nA,CA,nB,CB,nC,CC,nD,CD,fname)
! 4-index transformation out of core
! dumps all integrals on disk in the (triang,triang) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

integer,intent(in) :: NBas
integer,intent(in) :: nA
integer,intent(in) :: nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*)
double precision,intent(in) :: CB(*), CC(*), CD(*)
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

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,NBas*(NBas+1)/2))

 open(newunit=iunit,file='AOTWOSORT',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 ntr = NBas*(NBas+1)/2

 ! set no. of triangles in buffer
 nloop = (ntr-1)/cbuf+1

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
       call sq_to_triang(work2,work1,NBas)
       work3(rs,:) = work1(1:ntr)

    enddo

    ! transpose 
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
    call dgemm('T','N',nC,NBas,NBas,1d0,CC,NBas,work2,NBas,0d0,work1,nC)
    call dgemm('N','N',nC,nD,NBas,1d0,work1,nC,CD,NBas,0d0,work2,nC)
    call sq_to_triang(work2,work1,NBas)
    write(iunit3,rec=ab) work1(1:ntr)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran4_sym

subroutine tran_oneint(ints,MO1,MO2,NSym,GFunc,work,n)
implicit none

double precision,intent(inout),contiguous :: ints(:)
double precision,intent(in),contiguous :: MO1(:),MO2(:)
integer :: NSym, GFunc(:)
double precision,contiguous :: work(:)
integer :: irep,n,n2
integer :: offset

offset = 0
do irep=1,NSym
   n = GFunc(irep)
   n2 = n**2
   if(n>0) then
      call dgemm('T','N',n,n,n,&
           1d0,MO1(offset+1:offset+n2),n,ints(offset+1:offset+n2),n,&
           0d0,work,n)
      call dgemm('N','N',n,n,n,&
           1d0,work,n,MO2(offset+1:offset+n2),n,&
           0d0,ints(offset+1:offset+n2),n)
   endif
   offset = offset + n2
enddo

end subroutine tran_oneint

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

end module


