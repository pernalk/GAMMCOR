module tran
!!! CAREFUL :: unsym procedures yet to be tested!!!

!use types,only: EblockData
use blocktypes
use print_units
use sorter

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

subroutine tran3MO_Q(NBas,nocc,CX,Q,fname)
! 3-index transformation out of core
implicit none

integer,intent(in) :: NBas,nocc
!integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CX(*),Q(*)
character(*) :: fname
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: ntr,nsq,noccsq,nloop
integer,parameter :: cbuf=512
integer :: i,j,ij,rs,ab

 write(6,'()') 
! write(6,'(1x,a)') 'TRAN3'
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
 write(6,'(1x,a)') 'Transforming 3-ints for '//fname

 ntr = NBas*(NBas+1)/2
 nsq = NBas**2
 noccsq = nocc**2

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
     !access='DIRECT',form='UNFORMATTED',recl=8*nsq)
     access='DIRECT',form='UNFORMATTED',recl=8*noccsq)

 do ab=1,ntr

    do i=1,nloop
       ! get all (ab| for given |rs)
       read(iunit2,rec=(i-1)*ntr+ab) work1((i-1)*cbuf+1:min(i*cbuf,ntr))
    enddo

    call triang_to_sq(work1,work2,NBas)
    ! work1=CX^T.work2
    ! work2=work1.Q 
    !call dgemm('T','N',NBas,NBas,NBas,1d0,CX,NBas,work2,NBas,0d0,work1,NBas)
    !call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,Q,NBas,0d0,work2,NBas)

    !ij = 0
    !do j=1,nocc
    !   do i=1,nocc
    !      ij = ij + 1
    !      work1(ij) = work2((j-1)*NBas+i)
    !   enddo
    !enddo

    call dgemm('T','N',nocc,NBas,NBas,1d0,CX,NBas,work2,NBas,0d0,work1,nocc)
    call dgemm('N','N',nocc,nocc,NBas,1d0,work1,nocc,Q,NBas,0d0,work2,nocc)

    !write(iunit3,rec=ab) work2(1:nsq)
    write(iunit3,rec=ab) work2(1:noccsq)

 enddo

! write(LOUT,'(1x,a)') 'ACHTUNG-BABY! NEW TRAN3Q USED!'

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran3MO_Q

subroutine tran4_full(NBas,CA,CB,fname,srtfile)
! 4-index transformation out of core
! dumps all integrals on disk in the (triang,triang) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NBas
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*)

character(*)      :: fname, srtfile
logical           :: empty
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer           :: ntr,nloop
integer           :: i,rs,ab
integer           :: iunit,iunit2,iunit3
integer,parameter :: cbuf=512

 write(6,'()') 
! write(6,'(1x,a)') 'TRAN4_SYM_OUT_OF_CORE'
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
 write(6,'(1x,a)') 'Transforming integrals for '//fname

 ntr = NBas*(NBas+1)/2

 ! set no. of triangles in buffer
 nloop = (ntr-1)/cbuf+1

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,ntr))

 !open(newunit=iunit,file='AOTWOSORT',status='OLD',&
 !open(newunit=iunit,file=trim(srtfile),status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)
 call reader%open(trim(srtfile))

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (ab|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

       !read(iunit,rec=rs) work1(1:ntr)
       call reader%getTR(rs,work1(1:ntr),empty)
       if(empty) then
          work3(rs-(i-1)*cbuf,1:ntr) = 0
       else
          call triang_to_sq(work1,work2,NBas)
          ! work1=CA^T.work2
          ! work2=work1.CB
          call dgemm('T','N',NBas,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,NBas)
          call dgemm('N','N',NBas,NBas,NBas,1d0,work1,NBas,CA,NBas,0d0,work2,NBas)
          call sq_to_triang(work2,work1,NBas)
          ! transpose
          work3(rs-(i-1)*cbuf,1:ntr) = work1(1:ntr)
       endif

    enddo

    do ab=1,ntr
       write(iunit2,rec=(i-1)*ntr+ab) work3(1:cbuf,ab)
    enddo

 enddo

 !close(iunit)
 call reader%close

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

subroutine tran4_gen(NBas,nA,CA,nB,CB,nC,CC,nD,CD,fname,srtfile)
! 4-index transformation out of core
! dumps all integrals on disk in the (square,square) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
character(*) :: fname,srtfile
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: ntr,nAB,nCD,nloop
integer,parameter :: cbuf=512
integer :: i,rs,ab
logical :: empty
!  test
integer :: l,k,kl

! write(6,'()') 
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

 !!open(newunit=iunit,file='AOTWOSORT',status='OLD',&
 !open(newunit=iunit,file=trim(srtfile),status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 call reader%open(trim(srtfile))

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (ab|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

       work1 = 0
       !read(iunit,rec=rs) work1(1:ntr)
       call reader%getTR(rs,work1,empty)
       if(empty) then

          work3(rs-(i-1)*cbuf,1:nAB) = 0

       else

          call triang_to_sq(work1,work2,NBas)
          ! work1=CA^T.work2
          ! work2=work1.CB
          call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
          call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
          ! transpose
          work3(rs-(i-1)*cbuf,1:nAB) = work2(1:nAB)
          !work1 = 0
          !kl = 0 
          !do l=1,nB
          !do k=1,nA
          !   kl = kl + 1
          !   work1(kl) = work2((k-1)*NBas+l)
          !enddo
          !enddo
          !work3(rs-(i-1)*cbuf,1:nAB) = work1(1:nAB)

       endif

    enddo

    do ab=1,nAB
       write(iunit2,rec=(i-1)*nAB+ab) work3(1:cbuf,ab)
    enddo

 enddo

! close(iunit)
 call reader%close

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
    !work1 = 0
    !kl = 0 
    !do l=1,nD
    !do k=1,nC
    !   kl = kl + 1
    !   work1(kl) = work2((k-1)*NBas+l)
    !enddo
    !enddo
    !write(iunit3,rec=ab) work1(1:nCD)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine tran4_gen

subroutine new_tran4_gen(NBas,nA,CA,nB,CB,nC,CC,nD,CD, &
                         fname,reader,iunit_srt,thread_id)
! 4-index transformation out of core
! dumps all integrals on disk in the (square,square) form
! CAREFUL: C have to be in AOMO form!
!!! CAREFUL: write to simpler form
!!! ie. (NBas,n,C)
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
character(*) :: fname
!the corresponding file is open only once, before the call
integer, intent(in) :: iunit_srt
integer, intent(in), optional:: thread_id
! tmp file
character(len=80) :: tmpfile, thread_str
logical           :: empty
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer           :: i,rs,ab
integer           :: iunit,iunit2,iunit3
integer           :: ntr,nAB,nCD,nloop
integer,parameter :: cbuf=512
!  test
integer :: l,k,kl

! write(6,'()') 
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

 !!!open(newunit=iunit,file='AOTWOSORT',status='OLD',&
 !open(newunit=iunit,file=trim(srtfile),status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 tmpfile = 'TMPMO'
 iunit2 = 75  !no reason for this particular value
 !print *, "DEBUG: ", thread_id
 !print *, "DEBUG: unit", iunit2
 if (present(thread_id)) then
    write(thread_str, "(I0)") thread_id
    tmpfile = trim(tmpfile)//trim(adjustl(thread_str))
    iunit2 = iunit2 + thread_id
 end if
 !print *, "DEBUG: ", thread_id
 !print *, "DEBUG: before open ", tmpfile, "at", iunit2

 ! half-transformed file
 open(iunit2,file=tmpfile,status='REPLACE',&
      access='DIRECT',form='UNFORMATTED',recl=8*cbuf)
 !print *, "DEBUG: writing ", tmpfile, "at", iunit2

! (ab|
 do i=1,nloop
    ! loop over cbuf
    do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)
      !$omp critical
       !read(iunit_srt,rec=rs) work1(1:ntr)
       call reader%getTR(rs,work1,empty)
      !$omp end critical
       if(empty) then
          work3(rs-(i-1)*cbuf,1:nAB) = 0
       else
          call triang_to_sq(work1,work2,NBas)
          ! work1=CA^T.work2
          ! work2=work1.CB
          call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA)
          call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
          ! transpose
          work3(rs-(i-1)*cbuf,1:nAB) = work2(1:nAB)
          !work1 = 0
          !kl = 0
          !do l=1,nB
          !do k=1,nA
          !   kl = kl + 1
          !   work1(kl) = work2((k-1)*NBas+l)
          !enddo
          !enddo
          !work3(rs-(i-1)*cbuf,1:nAB) = work1(1:nAB)
       endif

    enddo

    do ab=1,nAB
       write(iunit2,rec=(i-1)*nAB+ab) work3(1:cbuf,ab)
    enddo

 enddo


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
    !work1 = 0
    !kl = 0
    !do l=1,nD
    !do k=1,nC
    !   kl = kl + 1
    !   work1(kl) = work2((k-1)*NBas+l)
    !enddo
    !enddo
    !write(iunit3,rec=ab) work1(1:nCD)

 enddo

 deallocate(work1,work2,work3)
 close(iunit3)
 close(iunit2,status='DELETE')

end subroutine new_tran4_gen

subroutine read4_gen(NBas,nA,nB,nC,nD,fname,srtfile)
!
! reads 4-index ints from sorted file
! dumps all integrals on disk in the (square,square) form
! CAREFUL: C have to be in AOMO form!
!
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
character(*) :: fname,srtfile
double precision, allocatable :: work1(:), work2(:)
integer :: iunit,iunit2
integer :: ntr,nAB,nCD
integer :: ip,iq,ir,is,pq,rs
logical :: empty

 write(6,'()')
 write(6,'(1x,a)') 'Reading integrals for '//fname

 ntr = NBas*(NBas+1)/2
 nAB = nA*nB
 nCD = nC*nD

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 !open(newunit=iunit,file=trim(srtfile),status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 call reader%open(trim(srtfile))

 open(newunit=iunit2,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nCD)

 do iq=1,nB
    do ip=1,nA
       pq=max(ip,iq)
       pq=min(ip,iq)+pq*(pq-1)/2

       !read(iunit,rec=pq) work1(1:ntr)
       call reader%getTR(pq,work1,empty)
       if(empty) then
          work1(1:nCD) = 0
       else
          call triang_to_sq(work1,work2,NBas)

          rs = 0
          do is=1,nD
             do ir=1,nC
                rs = rs + 1
                work1(rs) = work2(ir+(is-1)*NBas)
             enddo
          enddo
       endif

       pq=ip+(iq-1)*nA
       write(iunit2,rec=pq) work1(1:nCD)

    enddo
 enddo
 
 close(iunit2)
 !close(iunit)
 call reader%close

 deallocate(work1,work2)

end subroutine read4_gen

subroutine ABPM_TRAN(AMAT,AOUT,EBlock,EBlockIV,nblk,NDimX,isPl)
implicit none

integer,intent(in) :: nblk,NDimX
logical,intent(in) :: isPl
double precision,intent(in) :: AMAT(NDimX,NDimX)
double precision,intent(inout) :: AOUT(NDimX,NDimX)

type(EBlockData),intent(in) :: EBlock(nblk),EBlockIV

integer :: i,j,ii,jj,ipos,jpos,iblk,jblk
double precision,allocatable :: ABP(:,:),ABM(:,:)
double precision :: fac

fac = 1.d0/sqrt(2.d0)

AOUT=0

do jblk=1,nblk
   associate( jB => Eblock(jblk) )
   do iblk=1,nblk
      associate( iB => Eblock(iblk))

        allocate(ABP(iB%n,jB%n),ABM(iB%n,jB%n))
        do j=1,jB%n
           jpos = jB%pos(j)
           do i=1,iB%n
              ipos = iB%pos(i)
              ABP(i,j) = AMAT(ipos,jpos)
           enddo
        enddo
        if(isPl) then
           call dgemm('N','N',iB%n,jB%n,jB%n,1d0,ABP,iB%n,jB%matX,jB%n,0d0,ABM,iB%n)
           call dgemm('T','N',iB%n,jB%n,iB%n,1d0,iB%matX,iB%n,ABM,iB%n,0d0,ABP,iB%n)
        else
           call dgemm('N','N',iB%n,jB%n,jB%n,1d0,ABP,iB%n,jB%matY,jB%n,0d0,ABM,iB%n)
           call dgemm('T','N',iB%n,jB%n,iB%n,1d0,iB%matY,iB%n,ABM,iB%n,0d0,ABP,iB%n)
        endif

        AOUT(iB%l1:iB%l2,jB%l1:jB%l2) = ABP
        deallocate(ABM,ABP)

      end associate
   enddo
   end associate
enddo

associate(B => EblockIV)

  if(B%n>0) then
     do iblk=1,nblk
        associate(iB => Eblock(iblk))

          allocate(ABP(iB%n,B%n),ABM(iB%n,B%n))
          do j=1,B%n
             jpos = B%pos(j)
             do i=1,iB%n
                ipos = iB%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo
          if(isPl) then
             call dgemm('T','N',iB%n,B%n,iB%n,fac,iB%matX,iB%n,ABP,iB%n,0d0,ABM,iB%n) 
          else
             call dgemm('T','N',iB%n,B%n,iB%n,fac,iB%matY,iB%n,ABP,iB%n,0d0,ABM,iB%n) 
          endif

          AOUT(iB%l1:iB%l2,B%l1:B%l2) = ABM
          deallocate(ABM,ABP)

        end associate
     enddo

     do jblk=1,nblk
        associate(jB => Eblock(jblk))

          allocate(ABP(B%n,jB%n),ABM(B%n,jB%n))
          do j=1,jB%n
             jpos = jB%pos(j)
             do i=1,B%n
                ipos = B%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo
          if(isPl) then
             call dgemm('N','N',B%n,jB%n,jB%n,fac,ABP,B%n,jB%matX,jB%n,0d0,ABM,B%n)
          else
             call dgemm('N','N',B%n,jB%n,jB%n,fac,ABP,B%n,jB%matY,jB%n,0d0,ABM,B%n)
          endif

          AOUT(B%l1:B%l2,jB%l1:jB%l2) = ABM
          deallocate(ABM,ABP)

        end associate
     enddo

     do j=1,B%n
        jj = B%l1+j-1
        jpos = B%pos(j)
        do i=1,B%n
           ii = B%l1+i-1
           ipos = B%pos(i)
           AOUT(ii,jj) = AMAT(ipos,jpos)*0.5d0
        enddo
     enddo
  endif

end associate

end subroutine ABPM_TRAN

subroutine ABPM_DGEMV_GEN(AMAT,AOUT,EBlock,EBlockIV,&
                         nblk,NDimX,xyvar)
implicit none

integer,intent(in)      :: nblk,NDimX
character(*),intent(in) :: xyvar
double precision,intent(in)    :: AMAT(NDimX)
double precision,intent(inout) :: AOUT(NDimX)

type(EBlockData),intent(in)    :: EBlock(nblk),EBlockIV

integer :: i,j,ii,jj,ipos,jpos,iblk,jblk
double precision :: fact
double precision,allocatable :: ABP(:),ABM(:)

fact = 1d0/sqrt(2d0)

do iblk=1,nblk
   associate( A => Eblock(iblk) )

   allocate(ABP(A%n),ABM(A%n))
   do i=1,A%n
      ipos = A%pos(i)
      ABP(i) = AMAT(ipos)
   enddo

   select case(xyvar)
   case('X','x')
      call dgemv('T',A%n,A%n,1d0,A%matX,A%n,ABP,1,0d0,ABM,1)
   case('Y','y')
      call dgemv('T',A%n,A%n,1d0,A%matY,A%n,ABP,1,0d0,ABM,1)
   case('YX','yx')
      call dgemv('T',A%n,A%n,1d0,A%matY-A%matX,A%n,ABP,1,0d0,ABM,1)
   end select

   AOUT(A%l1:A%l2) = AOUT(A%l1:A%l2) + ABM

   deallocate(ABM,ABP)

   end associate
enddo

associate(A => EblockIV)

  if(A%n>0) then

     select case(xyvar)
     case('X','x')
        do i=1,A%n
           ii = A%l1+i-1
           ipos = A%pos(i)
           AOUT(ii) = AOUT(ii) + AMAT(ipos)*A%matX(i,1)
        enddo
    case('Y','y')
        do i=1,A%n
           ii = A%l1+i-1
           ipos = A%pos(i)
           AOUT(ii) = AOUT(ii) + AMAT(ipos)*A%matY(i,1)
        enddo
    case('YX','yx')
        do i=1,A%n
           ii = A%l1+i-1
           ipos = A%pos(i)
           AOUT(ii) = AOUT(ii) - AMAT(ipos)*fact
        enddo

    end select
 
  endif
end associate

end subroutine ABPM_DGEMV_GEN

subroutine ABPM_TRAN_GEN(AMAT,AOUT,EBlockA,EBlockAIV,EBlockB,EBlockBIV,&
                         nblkA,nblkB,ANDimX,BNDimX,xyvar)
implicit none

integer,intent(in)      :: nblkA,nblkB,ANDimX,BNDimX
character(*),intent(in) :: xyvar
double precision,intent(in)    :: AMAT(ANDimX,BNDimX)
double precision,intent(inout) :: AOUT(ANDimX,BNDimX)

type(EBlockData),intent(in)    :: EBlockA(nblkA),EBlockAIV, &
                                  EBlockB(nblkB),EBlockBIV

integer :: i,j,ii,jj,ipos,jpos,iblk,jblk
double precision,allocatable :: ABP(:,:),ABM(:,:)

do jblk=1,nblkB
   associate( B => EblockB(jblk) )
   do iblk=1,nblkA
      associate( A => EblockA(iblk))

        allocate(ABP(A%n,B%n),ABM(A%n,B%n))
        do j=1,B%n
           jpos = B%pos(j)
           do i=1,A%n
              ipos = A%pos(i)
              ABP(i,j) = AMAT(ipos,jpos)
           enddo
        enddo
        select case(xyvar)
        case('XX','xx')
           call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matX,B%n,0d0,ABM,A%n)
           call dgemm('T','N',A%n,B%n,A%n,1d0,A%matX,A%n,ABM,A%n,0d0,ABP,A%n)
        case('YY','yy')
           call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matY,B%n,0d0,ABM,A%n)
           call dgemm('T','N',A%n,B%n,A%n,1d0,A%matY,A%n,ABM,A%n,0d0,ABP,A%n)
        case('XY','xy')
           call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matY,B%n,0d0,ABM,A%n)
           call dgemm('T','N',A%n,B%n,A%n,1d0,A%matX,A%n,ABM,A%n,0d0,ABP,A%n)
        case('YX','yx')
           call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matX,B%n,0d0,ABM,A%n)
           call dgemm('T','N',A%n,B%n,A%n,1d0,A%matY,A%n,ABM,A%n,0d0,ABP,A%n)
        end select

        AOUT(A%l1:A%l2,B%l1:B%l2) = AOUT(A%l1:A%l2,B%l1:B%l2) + ABP
        deallocate(ABM,ABP)

      end associate
   enddo
   end associate
enddo

associate(B => EblockBIV)

  if(B%n>0) then
     do iblk=1,nblkA
        associate(A => EblockA(iblk))

          allocate(ABP(A%n,B%n),ABM(A%n,B%n))
          do j=1,B%n
             jpos = B%pos(j)
             do i=1,A%n
                ipos = A%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo

          select case(xyvar)
          case('XX','xx')
             call dgemm('T','N',A%n,B%n,A%n,1d0,A%matX,A%n,ABP,A%n,0d0,ABM,A%n)
             do j=1,B%n
                jj = B%l1+j-1
                AOUT(A%l1:A%l2,jj) = AOUT(A%l1:A%l2,jj) + ABM(1:A%n,j)*B%matX(j,1)
             enddo

          case('YY','yy')
             call dgemm('T','N',A%n,B%n,A%n,1d0,A%matY,A%n,ABP,A%n,0d0,ABM,A%n)
             do j=1,B%n
                jj = B%l1+j-1
                AOUT(A%l1:A%l2,jj) = AOUT(A%l1:A%l2,jj) + ABM(1:A%n,j)*B%matY(j,1)
             enddo

          case('XY','xy')
             call dgemm('T','N',A%n,B%n,A%n,1d0,A%matX,A%n,ABP,A%n,0d0,ABM,A%n)
             do j=1,B%n
                jj = B%l1+j-1
                AOUT(A%l1:A%l2,jj) = AOUT(A%l1:A%l2,jj) + ABM(1:A%n,j)*B%matY(j,1)
             enddo

          case('YX','yx')
             call dgemm('T','N',A%n,B%n,A%n,1d0,A%matY,A%n,ABP,A%n,0d0,ABM,A%n)
             do j=1,B%n
                jj = B%l1+j-1
                AOUT(A%l1:A%l2,jj) = AOUT(A%l1:A%l2,jj) + ABM(1:A%n,j)*B%matX(j,1)
             enddo
          end select

          deallocate(ABM,ABP)

        end associate
     enddo
  endif

end associate
associate(A => EblockAIV)

  if(A%n>0) then
     do jblk=1,nblkB
        associate(B => EblockB(jblk))

          allocate(ABP(A%n,B%n),ABM(A%n,B%n))
          do j=1,B%n
             jpos = B%pos(j)
             do i=1,A%n
                ipos = A%pos(i)
                ABP(i,j) = AMAT(ipos,jpos)
             enddo
          enddo

          select case(xyvar)
          case('XX','xx')
             call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matX,B%n,0d0,ABM,A%n)
             do i=1,A%n
                ii = A%l1+i-1
                AOUT(ii,B%l1:B%l2) = AOUT(ii,B%l1:B%l2) + ABM(i,1:B%n)*A%matX(i,1)
             enddo

          case('YY','yy')
             call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matY,B%n,0d0,ABM,A%n)
             do i=1,A%n
                ii = A%l1+i-1
                AOUT(ii,B%l1:B%l2) = AOUT(ii,B%l1:B%l2) + ABM(i,1:B%n)*A%matY(i,1)
             enddo

          case('XY','xy')
             call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matX,B%n,0d0,ABM,A%n)
             do i=1,A%n
                ii = A%l1+i-1
                AOUT(ii,B%l1:B%l2) = AOUT(ii,B%l1:B%l2) + ABM(i,1:B%n)*A%matY(i,1)
             enddo

          case('YX','yx')
             call dgemm('N','N',A%n,B%n,B%n,1d0,ABP,A%n,B%matY,B%n,0d0,ABM,A%n)
             do i=1,A%n
                ii = A%l1+i-1
                AOUT(ii,B%l1:B%l2) = AOUT(ii,B%l1:B%l2) + ABM(i,1:B%n)*A%matX(i,1)
             enddo

          end select

          deallocate(ABM,ABP)

        end associate
     enddo
  endif

end associate

associate(A => EblockAIV,&
          B => EblockBIV )

  if(A%n>0.and.B%n>0) then
     select case(xyvar)
     case('XX','xx')
        do j=1,B%n
           jj = B%l1+j-1
           jpos = B%pos(j)
           do i=1,A%n
              ii = A%l1+i-1
              ipos = A%pos(i)
              AOUT(ii,jj) = AOUT(ii,jj) + AMAT(ipos,jpos)*A%matX(i,1)*B%matX(j,1)
           enddo
        enddo

     case('YY','yy')
        do j=1,B%n
           jj = B%l1+j-1
           jpos = B%pos(j)
           do i=1,A%n
              ii = A%l1+i-1
              ipos = A%pos(i)
              AOUT(ii,jj) = AOUT(ii,jj) + AMAT(ipos,jpos)*A%matY(i,1)*B%matY(j,1)
           enddo
        enddo

     case('XY','xy')
        do j=1,B%n
           jj = B%l1+j-1
           jpos = B%pos(j)
           do i=1,A%n
              ii = A%l1+i-1
              ipos = A%pos(i)
              AOUT(ii,jj) = AOUT(ii,jj) + AMAT(ipos,jpos)*A%matX(i,1)*B%matY(j,1)
           enddo
        enddo

     case('YX','yx')
        do j=1,B%n
           jj = B%l1+j-1
           jpos = B%pos(j)
           do i=1,A%n
              ii = A%l1+i-1
              ipos = A%pos(i)
              AOUT(ii,jj) = AOUT(ii,jj) + AMAT(ipos,jpos)*A%matY(i,1)*B%matX(j,1)
           enddo
        enddo

     end select
  endif

end associate

end subroutine ABPM_TRAN_GEN

subroutine ABPM_HALFTRAN_GEN_L(AMAT,AOUT,fact,EBlock,EBlockIV,nblk,DimL,DimR,xyvar)
! AOUT = EBlock.AMAT
implicit none

integer,intent(in) :: nblk,DimL,DimR
double precision,intent(in) :: fact
character(*),intent(in) :: xyvar
double precision,intent(in)    :: AMAT(DimL,DimR)
double precision,intent(inout) :: AOUT(DimL,DimR)

type(EBlockData),intent(in) :: EBlock(nblk),EBlockIV

integer :: i,ii,ipos,iblk
integer :: maxEBlockN
double precision,allocatable :: ABP(:,:),ABM(:,:)

AOUT = fact*AOUT

maxEBlockN = maxval(EBlock%n)
if(maxEBlockN > 1) allocate(ABP(maxEBlockN,DimR),ABM(maxEBlockN,DimR))

select case(xyvar)
case('Y','y')

do iblk=1,nblk
   associate(iB => Eblock(iblk))

     if(iB%n==1) then

        ipos = iB%pos(1)
        AOUT(ipos,:) = iB%matY(1,1)*AMAT(ipos,:)

     else

        do i=1,iB%n
           ipos = iB%pos(i)
           ABP(i,:) = AMAT(ipos,:)
        enddo

        call dgemm('N','N',iB%n,DimR,iB%n,1d0,iB%matY,iB%n,ABP,maxEBlockN,0d0,ABM,maxEBlockN)

        do i=1,iB%n
           ipos = iB%pos(i)
           AOUT(ipos,:) = ABM(i,:)
        enddo

     endif

   end associate
enddo

case('X','x')

do iblk=1,nblk
   associate(iB => Eblock(iblk))

     if(iB%n==1) then

        ipos = iB%pos(1)
        AOUT(ipos,:) = iB%matX(1,1)*AMAT(ipos,:)

     else

        do i=1,iB%n
           ipos = iB%pos(i)
           ABP(i,:) = AMAT(ipos,:)
        enddo

        call dgemm('N','N',iB%n,DimR,iB%n,1d0,iB%matX,iB%n,ABP,maxEBlockN,0d0,ABM,maxEBlockN)

        do i=1,iB%n
           ipos = iB%pos(i)
           AOUT(ipos,:) = ABM(i,:)
        enddo

     endif

   end associate
enddo

case default
  write(lout,'(a)') 'Error in ABPM_HALFTRAN_LR!'
  stop
end select

if(maxEBlockN > 1) deallocate(ABM,ABP)

associate(B => EblockIV)

  do i=1,B%n
     !ipos = B%pos(i)
     !AOUT(ipos,:) = B%vec(i)*AMAT(ipos,:)
     ! MH 19.10.2022 change below is a bypass of some weird bug
     !               which is to be fixed!
     AOUT(B%pos(i),:) = B%vec(i)*AMAT(B%pos(i),:)
  enddo

end associate

end subroutine ABPM_HALFTRAN_GEN_L

subroutine ABPM_HALFTRAN_GEN_R(AMAT,AOUT,fact,EBlock,EBlockIV,nblk,DimL,DimR,xyvar)
implicit none

integer,intent(in) :: nblk,DimL,DimR
double precision,intent(in) :: fact
character(*),intent(in)     :: xyvar
double precision,intent(in)    :: AMAT(DimL,DimR)
double precision,intent(inout) :: AOUT(DimL,DimR)

type(EBlockData),intent(in) :: EBlock(nblk),EBlockIV

integer :: i,ii,ipos,iblk
integer :: maxEBlockN
double precision,allocatable :: ABP(:,:),ABM(:,:)

AOUT = fact*AOUT

maxEBlockN = maxval(EBlock%n)
if(maxEBlockN > 1) allocate(ABP(DimL,maxEBlockN),ABM(DimL,maxEBlockN))

select case(xyvar)
case('Y','y')

   do iblk=1,nblk
      associate(iB => Eblock(iblk))

        if(iB%n==1) then

           ipos = iB%pos(1)
           AOUT(:,ipos) = AOUT(:,ipos) + iB%matY(1,1)*AMAT(:,ipos)

        else
   
           do i=1,iB%n
              ipos = iB%pos(i)
              ABP(:,i) = AMAT(:,ipos)
           enddo
   
           call dgemm('N','N',DimL,iB%n,iB%n,1d0,ABP,DimL,iB%matY,iB%n,0d0,ABM,DimL)

           do i=1,iB%n
              ipos = iB%pos(i)
              AOUT(:,ipos) = AOUT(:,ipos) + ABM(:,i)
           enddo

        endif

      end associate
   enddo

case('X','x')

   do iblk=1,nblk
      associate(iB => Eblock(iblk))

        if(iB%n==1) then

           ipos = iB%pos(1)
           AOUT(:,ipos) = AOUT(:,ipos) + iB%matX(1,1)*AMAT(:,ipos)

        else

           do i=1,iB%n
              ipos = iB%pos(i)
              ABP(:,i) = AMAT(:,ipos)
           enddo
   
           call dgemm('N','N',DimL,iB%n,iB%n,1d0,ABP,DimL,iB%matX,iB%n,0d0,ABM,DimL)

           do i=1,iB%n
              ipos = iB%pos(i)
              AOUT(:,ipos) = AOUT(:,ipos) + ABM(:,i)
           enddo

        endif

      end associate
   enddo

case default
  write(lout,'(a)') 'Error in ABPM_HALFTRAN_LR!'
  stop
end select

if(maxEBlockN > 1) deallocate(ABM,ABP)

associate(B => EblockIV)

  do i=1,B%n
     !ipos = B%pos(i)
     !AOUT(:,ipos) = AOUT(:,ipos) + B%vec(i)*AMAT(:,ipos)
     ! MH 19.10.2022 change below is a bypass of some weird bug
     !               which is to be fixed!  
     AOUT(:,B%pos(i)) = AOUT(:,B%pos(i)) + B%vec(i)*AMAT(:,B%pos(i))
  enddo

end associate

end subroutine ABPM_HALFTRAN_GEN_R

subroutine MultpDiagMat(Diag,matA,fact,matC,dimA,dimB)
! multiply matC = diagonal*matrix
implicit none

integer,intent(in) :: dimA,dimB
double precision,intent(in)    :: fact
double precision,intent(in)    :: Diag(dimA),matA(dimA,dimB)
double precision,intent(inout) :: matC(dimA,dimB)

integer :: i,j

matC = fact * matC
do j=1,dimB
   do i=1,dimA
      matC(i,j) = matC(i,j) + Diag(i)*matA(i,j)
   enddo
enddo

end subroutine MultpDiagMat

subroutine make_J1(NBas,X,J,intfile)
implicit none

integer :: NBas
double precision :: X(*), J(NBas,NBas)
character(*) :: intfile

type(AOReaderData) :: reader
integer :: iunit, ntr
integer :: ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
double precision,external :: ddot
! works in DCBS

 J = 0
 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 !open(newunit=iunit,file=trim(intfile),status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 !     
 call reader%open(trim(intfile))

 irs=0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    !read(iunit,rec=irs) work1(1:ntr)
    call reader%getTR(irs,work1,empty)
    if(empty) cycle
    call triang_to_sq(work1,work2,NBas)

    tmp = ddot(NBas**2,work2,1,X,1)
    J(ir,is) = tmp
    J(is,ir) = tmp

    enddo
 enddo

 deallocate(work1,work2)
 !close(iunit)
 call reader%close

end subroutine make_J1

subroutine make_J2(NBas,XA,XB,JA,JB)
implicit none

type(AOReaderData) :: reader

integer :: NBas
double precision :: XA(*), XB(*), JA(NBas,NBas), JB(NBas,NBas)
integer :: iunit, ntr
integer :: ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
double precision,external :: ddot
! works in DCBS

 JA = 0
 JB = 0
 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 !open(newunit=iunit,file='AOTWOSORT',status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 !
 call reader%open('AOTWOSORT')

 irs=0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    !read(iunit,rec=irs) work1(1:ntr)
    call reader%getTR(irs,work1,empty)
    if(empty) cycle
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
 !close(iunit)
 call reader%close

end subroutine make_J2

subroutine make_K(NBas,X,K)
implicit none

type(AOReaderData) :: reader
integer,intent(in) :: NBas
double precision,intent(in) :: X(NBas,NBas)
double precision,intent(inout) :: K(NBas,NBas)
integer :: iunit, ntr
integer :: ir,is,irs
logical :: empty
double precision :: tmp
double precision,allocatable :: work1(:),work2(:)
! works in DCBS

 K = 0
 ntr = NBas*(NBas+1)/2

 allocate(work1(NBas*NBas),work2(NBas*NBas))

 !open(newunit=iunit,file='AOTWOSORT',status='OLD',&
 !     access='DIRECT',form='UNFORMATTED',recl=8*ntr)
 !
 call reader%open('AOTWOSORT')

 irs = 0
 do is=1,NBas
    do ir=1,is
    irs = irs + 1
    !read(iunit,rec=irs) work1(1:ntr)
    call reader%getTR(irs,work1,empty)
    if(empty) cycle
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
 !close(iunit)
 call reader%close

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

subroutine tran_matTr(ints,MO1,MO2,nbas,tran)
implicit none

integer,intent(in) :: nbas
logical,intent(in) :: tran
double precision,intent(inout) :: ints(nbas*(nbas+1)/2)
double precision,intent(in) :: MO1(nbas,nbas),MO2(nbas,nbas)
double precision,allocatable :: work1(:),work2(:)

allocate(work1(nbas*nbas),work2(nbas*nbas))
 call triang_to_sq(ints,work1,nbas)

 if(nbas>0) then
    if(tran) then
       call dgemm('T','N',nbas,nbas,nbas,1d0,MO1,nbas,work1,nbas,0d0,work2,nbas)
       call dgemm('N','N',nbas,nbas,nbas,1d0,work2,nbas,MO2,nbas,0d0,work1,nbas)
    else
       call dgemm('N','N',nbas,nbas,nbas,1d0,MO1,nbas,work1,nbas,0d0,work2,nbas)
       call dgemm('N','T',nbas,nbas,nbas,1d0,work2,nbas,MO2,nbas,0d0,work1,nbas)
    endif    
 endif

 call sq_to_triang(work1,ints,nbas)
deallocate(work2,work1)

end subroutine tran_matTr

subroutine tran_AO2MO(Ain,CA,CB,Aout,NAO,NMO)
implicit none

integer,intent(in) :: NAO, NMO
double precision,intent(in)  :: Ain(:),CA(:),CB(:)
double precision,intent(out) :: Aout(:)

double precision             :: work(NAO*NMO)

if(NAO>0.and.NMO>0) then
  call dgemm('T','N',NMO,NAO,NAO,1d0,CA,NAO,Ain,NAO,0d0,work,NMO)
  call dgemm('N','N',NMO,NMO,NAO,1d0,work,NMO,CB,NAO,0d0,Aout,NMO)
endif

end subroutine tran_AO2MO

subroutine tran_AO2MO2(Ain,CA,CB,Aout,NAO,NMO)
implicit none

integer,intent(in) :: NAO, NMO
double precision,intent(in)  :: Ain(:,:),CA(:,:),CB(:,:)
double precision,intent(out) :: Aout(:,:)

double precision :: work(NMO,NAO)

if(NAO>0.and.NMO>0) then
  call dgemm('T','N',NMO,NAO,NAO,1d0,CA,NAO,Ain,NAO,0d0,work,NMO)
  call dgemm('N','N',NMO,NMO,NAO,1d0,work,NMO,CB,NAO,0d0,Aout,NMO)
endif

end subroutine tran_AO2MO2

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

subroutine transp_mat1dim(matIn,matOut,NBas)
implicit none

integer,intent(in) :: NBas
double precision,intent(in) :: matIn(NBas,NBas)
double precision,intent(out) :: matOut(NBas,NBas)
integer :: i,j

 do i=1,NBas
    do j=1,NBas
       matOut(j,i) = matIn(i,j)
    enddo
 enddo

end subroutine transp_mat1dim

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
       matSq(pq) = matTr(t)
    enddo
 enddo

 t = 0
 do q=1,NBas
    do p=1,q-1
       t = t + 1
       qp = q + (p-1)*NBas
       matSq(qp) = matTr(t)
    enddo
    t = t + 1
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

end module

