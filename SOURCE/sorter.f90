module sorter
use types

implicit none

type InfoData
integer :: nelms
integer :: irec
end type InfoData

type BatchData
integer :: nelms
integer,allocatable :: ival(:)
double precision, allocatable :: rval(:)
integer :: ninfo
type(InfoData), allocatable :: Info(:)
end type BatchData

type SortData
integer :: iunit
integer :: irec
integer :: nelms_max = 512
integer :: nBatch
type(BatchData),allocatable :: Batch(:)
integer,allocatable :: ival(:)
double precision, allocatable :: rval(:)
end type SortData

contains 

subroutine test2el(NBas)
! sorts 2-el integrals
! according to rs index
! and dumps them to file
implicit none

integer :: NBas
integer :: iunit,iunit2
integer iunit77,iunit88, iunit99
integer :: maxrep, naos(8), lbuf, nibuf, nbits
type(SortData) :: srt
integer :: nints, INDX
integer,allocatable :: idx_buf(:)
double precision,allocatable :: val_buf(:), mat(:)
integer :: idx_p, idx_q, idx_r, idx_s, pq, rs
logical :: swap_pqrs
integer :: i

 ! start SORT1
 ! open sorter file 
 call open_Sorter(srt,'TMPSORT',NBas*(NBas+1)/2,NBas*(NBas+1)/2)

 ! newunit works with Fortran 2008
 open(newunit=iunit,file='AOTWOINT',status='OLD',&
      access='SEQUENTIAL',form='UNFORMATTED')

 ! read info
 call readlabel(iunit,'BASINFO ')
 read(iunit) maxrep, naos, lbuf, nibuf, nbits 

 write(6,'(1x, a)') 'Dalton two-el. file initialized'
 write(6,'(1x,a,i3,a,i3,a,i3)') 'Buffer size: ', lbuf, &
       & ', integers per index packet: ', nibuf,       &
       & ', bits: ', nbits


 allocate(val_buf(lbuf))
 allocate(idx_buf(lbuf*nibuf))

 call readlabel(iunit,'BASTWOEL')

 select case(nibuf)
 case(1)

    do
       read(iunit) val_buf, idx_buf, nints
       if(nints<0) exit
       do i=1,nints
          INDX = idx_buf(i)
          idx_p = ibits(INDX,0,8)
          idx_q = ibits(INDX,8,8)
          idx_r = ibits(INDX,16,8)
          idx_s = ibits(INDX,24,8)

          ! pq: position in Batch
          ! rs: Batch number
          pq = idx_p + idx_q*(idx_q-1)/2
          rs = idx_r + idx_s*(idx_s-1)/2 
          swap_pqrs = (pq<rs)
          call add_to_Sorter(srt,rs,pq,val_buf(i))
          if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(i))
       enddo
    enddo
   
 case(2)

    do
       read(iunit) val_buf, idx_buf, nints
       if(nints<0) exit
       do i=1,nints
          INDX = idx_buf(i)
          idx_r = ibits(INDX,0,16)
          idx_s = ibits(INDX,16,16)
          INDX = idx_buf(i+lbuf)
          idx_p = ibits(INDX,0,16)
          idx_q = ibits(INDX,16,16)

          pq = idx_p*(idx_p-1)/2+idx_q 
          rs = idx_r*(idx_r-1)/2+idx_s 
          swap_pqrs = (pq<rs)
          call add_to_Sorter(srt,rs,pq,val_buf(i))
          if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(i))

       enddo
    enddo

 end select

 deallocate(val_buf,idx_buf)
! close AOTWOINT
 close(unit=iunit)
 call dump_Sorter(srt)
! end SORT1

! start SORT2
 open(newunit=iunit2,file='AOTWOSORT',status='REPLACE',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 allocate(mat(NBas*(NBas+1)/2))
 do rs=1,NBas*(NBas+1)/2
    call get_from_Sorter(srt,rs,mat)
    ! dump sorted integrals 
    write(iunit2,rec=rs) mat
 enddo

 call free_Sorter(srt)

 deallocate(mat)
 close(unit=iunit2)
! end SORT2

end subroutine test2el

!subroutine readlabel(iunit,text)
!! sets file pointer 
!! to first data after text
!implicit none
!
!integer :: iunit
!integer :: ios
!character(8) :: text, label(4)
!
!rewind(iunit)
!do 
!
!  read(iunit,iostat=ios) label
!  if(ios<0) then
!     write(6,*) 'ERROR!!! Empty section in AOTWOINT!'
!     stop
!  endif
!  if(label(1)=='********') then
!     if(label(4)==text) exit
!  endif
!
!enddo
!
!end subroutine readlabel

subroutine open_Sorter(srt,name,nBatch,BatchSize)
implicit none
type(SortData) :: srt
character(*) :: name
integer :: nBatch,BatchSize
integer :: InfoSize
integer :: i

 ! recl depends on the number of real(8) numbers
 ! that we store in it
 open(newunit=srt%iunit,file=trim(name),status='REPLACE',&
      access='DIRECT',form='UNFORMATTED',recl=8*srt%nelms_max)
 srt%irec = 0

 srt%nBatch = nBatch
 allocate(srt%Batch(srt%nBatch))

 InfoSize = (BatchSize-1)/srt%nelms_max + 1

 ! init Batch
 do i=1,srt%nBatch
    associate(Bucket => srt%Batch(i))
      Bucket%nelms = 0
      allocate(Bucket%ival(srt%nelms_max))
      allocate(Bucket%rval(srt%nelms_max))
      Bucket%ninfo = 0
      allocate(Bucket%Info(InfoSize))
    end associate
 enddo

end subroutine open_Sorter

subroutine add_to_Sorter(srt,ibatch,ipos,VAL)
implicit none

type(SortData) :: srt
integer :: ibatch, ipos
double precision :: VAL

associate(Bucket => srt%Batch(ibatch))

  if(Bucket%nelms==srt%nelms_max) call dump_Batch(srt%iunit,srt%irec,Bucket)
  Bucket%nelms = Bucket%nelms + 1
  Bucket%ival(Bucket%nelms) = ipos
  Bucket%rval(Bucket%nelms) = VAL

end associate

end subroutine add_to_Sorter

subroutine dump_Batch(iunit,irec,Batch)
implicit none 

integer :: iunit, irec
type(BatchData) :: Batch

if(Batch%nelms>0) then

   Batch%ninfo = Batch%ninfo + 1
   Batch%Info(Batch%ninfo)%nelms = Batch%nelms
   Batch%Info(Batch%ninfo)%irec = irec
   
   irec = irec + 1
   write(iunit,rec=irec) Batch%rval(1:Batch%nelms)
   irec = irec + 1
   write(iunit,rec=irec) Batch%ival(1:Batch%nelms) 
  
   Batch%nelms = 0

endif

end subroutine dump_Batch


subroutine dump_Sorter(srt) 
implicit none

type(SortData) :: srt
integer :: i

do i=1,srt%nBatch 

   associate(Bucket => srt%Batch(i))
     call dump_Batch(srt%iunit,srt%irec,Bucket)
     deallocate(Bucket%ival,Bucket%rval)
   end associate

enddo

! init SORT2
allocate(srt%ival(srt%nelms_max))
allocate(srt%rval(srt%nelms_max))

end subroutine dump_Sorter 

subroutine get_from_Sorter(srt,rs,mat)
implicit none

type(SortData) :: srt
integer :: rs
double precision :: mat(:)
integer :: i,j

mat = 0
associate(Bucket => srt%Batch(rs))

  do i=1,Bucket%ninfo
     associate( Info => Bucket%Info(i))
     read(srt%iunit,rec=Info%irec+1) srt%rval(1:Info%nelms)
     read(srt%iunit,rec=Info%irec+2) srt%ival(1:Info%nelms)
     do j=1,Info%nelms
        mat(srt%ival(j)) = srt%rval(j)
     enddo
     end associate
  enddo

end associate

end subroutine get_from_Sorter

subroutine free_Sorter(srt)
implicit none

type(SortData) :: srt
integer :: i

 ! release Info
 do i=1,srt%nBatch
    deallocate(srt%Batch(i)%Info)
 enddo

 deallocate(srt%ival,srt%rval)
 deallocate(srt%Batch)
 ! close SORT 
 close(srt%iunit)

end subroutine free_Sorter

!subroutine tran4_unsym(NBas,nA,CA,nB,CB,nC,CC,nD,CD,TNO)
!implicit none
!
!integer,intent(in) :: NBas
!integer,intent(in) :: nA,nB,nC,nD
!! CA(NBas*nA)
!double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
!double precision,intent(out) :: TNO(:,:)
!double precision, allocatable :: work1(:), work2(:), work3(:,:)
!integer :: iunit
!integer :: r,s,rs,d
!
!! 4-index transformation to NOs
!! integrals stored in TNO 
! 
! write(6,'()') 
! write(6,*) NBas,nA,nB,nC,nD,'NBas...'
! write(6,*) CA(2), CB(2), CC(2), CD(2)
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
!
! allocate(work1(NBas*NBas),work2(NBas*NBas))
! allocate(work3(nA*nB,nC))
!
! open(newunit=iunit,file='AOTWOSORT',status='OLD',&
!      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)
!
! TNO = 0
! do s=1,NBas
!    work3 = 0
!
!    do r=1,NBas
!       rs = min(r,s) + max(r,s)*(max(r,s)-1)/2 ! idx_tiang
!       read(iunit,rec=rs) work1(1:NBas*(NBas+1)/2)
!       call triang_to_sq(work1,work2,NBas)
!       ! work1=CA^T.work2
!       call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
!       ! work2=work1.CB
!       call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
!       ! ... 
!       call dger(nA*nB,nC,1d0,work2,1,CC(r:nC*NBas),NBas,work3,nA*nB)
!   
!    enddo
!    
!    !call dger(nA*nB*nC,nD,1d0,work3,1,CD(s),NBas,TNO,nA*nB*nC)
!    do d=1,nD
!       TNO(:,(d-1)*nC+1:d*nC) = TNO(:,(d-1)*nC+1:d*nC) + work3*CD(s+(d-1)*NBas)
!    enddo
!
! enddo
!
! deallocate(work1,work2,work3)
! close(iunit)
!
!end subroutine tran4_unsym
!
!subroutine tran4_unsym2(NBas,nA,CA,nB,CB,nC,CC,nD,CD,TNO)
!implicit none
!
!integer,intent(in) :: NBas
!integer,intent(in) :: nA,nB,nC,nD
!! CA(NBas*nA)
!double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
!double precision,intent(out) :: TNO(:,:)
!double precision, allocatable :: work1(:), work2(:), work3(:,:)
!integer :: iunit
!integer :: r,s,rs,d
!
!! 4-index transformation to NOs
!! integrals stored in TNO 
! 
! write(6,'()') 
! write(6,'(1x,a)') 'Transforming integrals for AB dimer'
!
! allocate(work1(NBas*NBas),work2(NBas*NBas))
! allocate(work3(nA*nB,nC))
!
! open(newunit=iunit,file='AOTWOSORT',status='OLD',&
!      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)
!
! TNO = 0
! do s=1,NBas
!    work3 = 0
!
!    do r=1,NBas
!       rs = min(r,s) + max(r,s)*(max(r,s)-1)/2 ! idx_tiang
!       read(iunit,rec=rs) work1(1:NBas*(NBas+1)/2)
!       call triang_to_sq(work1,work2,NBas)
!       ! work1=CA.work2
!       call dgemm('N','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA) 
!       ! work2=work1.CB^T
!       call dgemm('N','T',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
!       ! ... 
!       call dger(nA*nB,nC,1d0,work2,1,CC((r-1)*NBas+1),1,work3,nA*nB)
!   
!    enddo
!    
!    call dger(nA*nB*nC,nD,1d0,work3,1,CD((s-1)*NBas+1),1,TNO,nA*nB*nC)
!!    do d=1,nD
!!       TNO(:,(d-1)*nC+1:d*nC) = TNO(:,(d-1)*nC+1:d*nC) + work3*CD(s+(d-1)*NBas)
!!    enddo
!
! enddo
!
! deallocate(work1,work2,work3)
! close(iunit)
!
!end subroutine tran4_unsym2
!
!
!subroutine triang_to_sq(matTr,matSq,NBas)
!implicit none
!
!double precision :: matTr(:), matSq(:)
!integer :: NBas
!integer :: t,p,q,pq,qp
!
! t = 0
! do q=1,NBas
!    do p=1,q
!       t = t + 1
!       pq = p + (q-1)*NBas
!       qp = q + (p-1)*NBas
!       matSq(pq) = matTr(t)
!       matSq(qp) = matTr(t)
!    enddo
! enddo
!
!end subroutine triang_to_sq

end module



