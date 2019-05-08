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

subroutine readtwoint(NBas,aosource,intfile,sortfile,outinfo)
! read DALTON 2-el integrals,
! sort them according to rs index
! and dump to a file
implicit none

integer :: NBas
integer :: aosource
character(*) :: intfile,sortfile
integer(8),optional :: outinfo
integer :: iunit,iunit2
integer :: maxrep, naos(8), lbuf, nibuf, nbits, lenint4
integer :: nsk, nt(8), ntoff(8)
type(SortData) :: srt
integer :: nints, INDX
integer,allocatable :: idx_buf(:)
double precision :: val
double precision,allocatable :: val_buf(:), mat(:)
integer :: idx_p, idx_q, idx_r, idx_s, pq, rs, idx_end
integer :: sym_p, sym_q, sym_r, sym_s, sym_rs
logical :: swap_pqrs
integer(8) :: i

 ! start SORT1
 ! open sorter file 
 call open_Sorter(srt,'TMPSORT',NBas*(NBas+1)/2,NBas*(NBas+1)/2)


 if(aosource==1) then 
   ! DALTON
   open(newunit=iunit,file=intfile,status='OLD',&
        access='SEQUENTIAL',form='UNFORMATTED')
  
   ! read info
   call readlabel(iunit,'BASINFO ')
   read(iunit) maxrep, naos, lbuf, nibuf, nbits, lenint4 
  
   write(6,'()')
   write(6,'(1x, a)') 'Dalton two-el. file initialized'
   write(6,'(1x,4(a,i5))') 'Buffer size: ', lbuf, &
         & ', integers per index packet: ', nibuf,       &
         & ', bits: ', nbits, &
         & ', record lenght in INT*4: ', lenint4
  
   allocate(val_buf(lbuf))
   allocate(idx_buf(nibuf*lbuf))
  
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
            !write(*,*) idx_p,idx_q,idx_r,idx_s 
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
            INDX = idx_buf(2*i-1)
            idx_r = ibits(INDX,0,16)
            idx_s = ibits(INDX,16,16)
            INDX = idx_buf(2*i)
            idx_p = ibits(INDX,0,16)
            idx_q = ibits(INDX,16,16)
  
            pq = idx_p + idx_q*(idx_q-1)/2
            rs = idx_r + idx_s*(idx_s-1)/2 
            swap_pqrs = (pq<rs)
            call add_to_Sorter(srt,rs,pq,val_buf(i))
            if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(i))
  
         enddo
      enddo
  
   end select
  
   deallocate(idx_buf)
   deallocate(val_buf)
  ! close AOTWOINT
   close(unit=iunit)

 elseif(aosource==2) then
 ! MOLPRO

   open(newunit=iunit,file=intfile,status='OLD',&
        access='SEQUENTIAL',form='UNFORMATTED')
  
   ! read info
   read(iunit) nsk
   read(iunit) nt 

   if(sum(nt(1:nsk))/=NBas) then
     write(LOUT,*) 'ERROR while checking AOTWOINT.mol!!!'
     stop
   endif
  
   ! offset
   ntoff = 0
   do i=2,8
      ntoff(i) = ntoff(i-1)+nt(i-1)
   enddo

   lbuf = 3*maxval(nt(1:nsk))**2
   allocate(val_buf(lbuf))

   do sym_p=1,nsk
      if(nt(sym_p)==0) cycle
      read(iunit) 

      do  
        read(iunit) idx_s,idx_r
        if(idx_s<=0.and.idx_r<=0) exit
        nints = idx_r + idx_s*(idx_s-1)/2 
        read(iunit) val_buf(1:nints)

        idx_r = idx_r + ntoff(sym_p)
        idx_s = idx_s + ntoff(sym_p)
        rs = idx_r + idx_s*(idx_s-1)/2

        INDX = 0 
        do idx_q=ntoff(sym_p)+1,idx_s 
           if(idx_q<idx_s) then 
             idx_end = idx_q
           else
             idx_end = idx_r
           endif
           do idx_p=ntoff(sym_p)+1,idx_end

              INDX = INDX + 1
              if(val_buf(INDX)/=0) then
                pq = idx_p + idx_q*(idx_q-1)/2
                swap_pqrs = (pq<rs)
                call add_to_Sorter(srt,rs,pq,val_buf(INDX))
                if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(INDX))
              endif

           enddo
        enddo

      enddo

   enddo

   do sym_r=2,nsk
     do sym_p=1,sym_r-1
       if(nt(sym_p)==0) cycle
       if(nt(sym_r)==0) cycle
       read(iunit)
       
       nints=(nt(sym_p)*(nt(sym_p)+1))/2
       nints=3*nints
       do
         read(iunit) idx_s,idx_r
         if(idx_s<=0.and.idx_r<=0) exit
         read(iunit) val_buf(1:nints)

         idx_r = idx_r + ntoff(sym_r)
         idx_s = idx_s + ntoff(sym_r)

         INDX=0
         do idx_q=ntoff(sym_p)+1,ntoff(sym_p)+nt(sym_p)
           do idx_p=ntoff(sym_p)+1,idx_q
 
             INDX=INDX+1
             if(val_buf(INDX)/=0) then
               pq = idx_p + idx_q*(idx_q-1)/2
               rs = idx_r + idx_s*(idx_s-1)/2
               call add_to_Sorter(srt,rs,pq,val_buf(INDX))
               call add_to_Sorter(srt,pq,rs,val_buf(INDX))
             endif

             INDX=INDX+1
             if(val_buf(INDX)/=0) then
               pq = idx_p + idx_r*(idx_r-1)/2
               rs = idx_q + idx_s*(idx_s-1)/2
               swap_pqrs = (pq<rs)
               call add_to_Sorter(srt,rs,pq,val_buf(INDX))
               if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(INDX))
             endif

             INDX=INDX+1
             if(val_buf(INDX)/=0) then
               pq = idx_q + idx_r*(idx_r-1)/2
               rs = idx_p + idx_s*(idx_s-1)/2
               swap_pqrs = (pq<rs)
               call add_to_Sorter(srt,rs,pq,val_buf(INDX))
               if(swap_pqrs) call add_to_Sorter(srt,pq,rs,val_buf(INDX))
             endif

           enddo
         enddo

       enddo

     enddo
   enddo

   do sym_s=4,nsk
     do sym_r=3,sym_s-1
       sym_rs = ieor(sym_r-1,sym_s-1) + 1
       do sym_q=2,sym_r-1
         sym_p = ieor(sym_q-1,sym_rs-1) + 1
         if(sym_p>=sym_q) cycle
         if(nt(sym_p)==0) cycle
         if(nt(sym_q)==0) cycle
         if(nt(sym_r)==0) cycle
         if(nt(sym_s)==0) cycle
         read(iunit)
         
         nints=nt(sym_p)*nt(sym_q)
         nints=3*nints
         do
           read(iunit) idx_s,idx_r
           if(idx_s<=0.and.idx_r<=0) exit
           read(iunit) val_buf(1:nints)

           idx_r = idx_r + ntoff(sym_r)
           idx_s = idx_s + ntoff(sym_s)

           INDX=0
           do idx_q=ntoff(sym_q)+1,ntoff(sym_q)+nt(sym_q)
             do idx_p=ntoff(sym_p)+1,ntoff(sym_p)+nt(sym_p)
 
               INDX=INDX+1
               if(val_buf(INDX)/=0) then
                 pq = idx_p + idx_q*(idx_q-1)/2
                 rs = idx_r + idx_s*(idx_s-1)/2
                 call add_to_Sorter(srt,rs,pq,val_buf(INDX))
                 call add_to_Sorter(srt,pq,rs,val_buf(INDX))
               endif

               INDX=INDX+1
               if(val_buf(INDX)/=0) then
                 pq = idx_p + idx_r*(idx_r-1)/2
                 rs = idx_q + idx_s*(idx_s-1)/2
                 call add_to_Sorter(srt,rs,pq,val_buf(INDX))
                 call add_to_Sorter(srt,pq,rs,val_buf(INDX))
               endif

               INDX=INDX+1
               if(val_buf(INDX)/=0) then
                 pq = idx_q + idx_r*(idx_r-1)/2
                 rs = idx_p + idx_s*(idx_s-1)/2
                 call add_to_Sorter(srt,rs,pq,val_buf(INDX))
                 call add_to_Sorter(srt,pq,rs,val_buf(INDX))
               endif

             enddo
           enddo

         enddo

       enddo
     enddo
   enddo


   deallocate(val_buf)

   close(iunit)

   !print*, 'Go Werner, go!'

 elseif(aosource==3) then
 ! Eugene
   open(newunit=iunit,file=intfile,status='OLD',&
        access='STREAM',form='UNFORMATTED')
  
      do
         read(iunit) val,idx_p,idx_q,idx_r,idx_s
         if(idx_r+idx_s==0) exit
  
         ! pq: position in Batch
         ! rs: Batch number
         pq = min(idx_p,idx_q) + max(idx_p,idx_q)*(max(idx_p,idx_q)-1)/2
         rs = min(idx_r,idx_s) + max(idx_r,idx_s)*(max(idx_r,idx_s)-1)/2
         call add_to_Sorter(srt,rs,pq,val)
         call add_to_Sorter(srt,pq,rs,val)
      enddo

   ! find position in file
   inquire(unit=iunit,pos=i)
   if(present(outinfo)) outinfo = i-24

   close(iunit)

 endif

 call dump_Sorter(srt)
! end SORT1

! start SORT2
 open(newunit=iunit2,file=trim(sortfile),status='REPLACE',&
      access='DIRECT',form='UNFORMATTED',recl=8*NBas*(NBas+1)/2)

 allocate(mat(NBas*(NBas+1)/2))
 do rs=1,NBas*(NBas+1)/2
    call get_from_Sorter(srt,rs,mat)
    ! dump sorted integrals
    ! test print all integrals:
    !do pq=1, NBas*(NBas+1)/2
    !   write(LOUT,'(1x,2i6,es20.7)') pq,rs,mat(pq)
    !enddo 
    write(iunit2,rec=rs) mat
 enddo
 deallocate(mat)

 call free_Sorter(srt)

 close(unit=iunit2)
! end SORT2

end subroutine readtwoint

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
integer :: InfoSize,IRValSize
integer :: i

 ! recl depends on the number of real(8) numbers
 ! that we store in it
 IRValSize = ( storage_size(srt%ival) + storage_size(srt%rval) ) / 8

 open(newunit=srt%iunit,file=trim(name),status='REPLACE',&
      access='DIRECT',form='UNFORMATTED',recl=IRValSize*srt%nelms_max)
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

   irec = irec + 1
   Batch%ninfo = Batch%ninfo + 1
   Batch%Info(Batch%ninfo)%nelms = Batch%nelms
   Batch%Info(Batch%ninfo)%irec = irec
   
   write(iunit,rec=irec) Batch%rval(1:Batch%nelms),Batch%ival(1:Batch%nelms)
  
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
     read(srt%iunit,rec=Info%irec) srt%rval(1:Info%nelms),srt%ival(1:Info%nelms)
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
 ! close TMPSORT 
 close(srt%iunit,status='DELETE')

end subroutine free_Sorter

end module


