module sorter
use types
implicit none

private
public readtwoint
public AOReaderData

type DumpData
integer(8) :: n,pos
end type DumpData

type,extends(DumpData) ::  DumpData_queue
type(DumpData_queue),pointer :: next
end type DumpData_queue

type AOBucketData
integer(8) :: n,ntot
integer :: ndmp
type(DumpData_queue),pointer :: first,last  ! pass1
type(DumpData),allocatable :: dmp(:)        !       pass2
integer,allocatable :: pqrs(:,:)            ! pass1        -> bytes
double precision,allocatable :: val(:)      ! pass1        -> bytes
contains
procedure :: dump => dump_AOBucket
end type AOBucketData

type AOSorterData
integer :: unit
integer :: nbas,nelm,ncol,nbck,maxdmp
integer(8) :: seg_size,bck_size
type(AOBucketData),allocatable :: bck(:)
double precision,allocatable :: rs_diag(:)  ! pass1 pass2  -> avail
integer,allocatable :: rs_nelm(:)           ! pass1 pass2  -> avail
integer,allocatable :: seg_pq(:)            !       pass2  -> bytes
double precision,allocatable :: seg_val(:)  !       pass2  -> bytes
integer,allocatable :: pqrs(:,:)            !       pass2  -> ???
double precision,allocatable :: val(:)      !       pass2  -> ???
integer(8),allocatable :: rs_pos(:)         !       pass2  -> avail
contains
procedure :: open => open_AOSorter
procedure :: add => add_AOSorter
procedure :: close => close_AOSorter
end type AOSorterData

type AOReaderData
integer :: unit
integer :: nbas,nelm
integer(8) :: pos_diag,pos_nelm,pos_pos
integer,allocatable :: rs_nelm(:)
integer(8),allocatable :: rs_pos(:)
integer,allocatable :: pq(:)
double precision,allocatable :: val(:)
contains
procedure :: open => open_AOReader
procedure :: get => get_AOReader
procedure :: close => close_AOReader
end type AOReaderData

contains

subroutine readtwoint(nbas,aosource,intfile,sortfile,maxmem,outinfo)
implicit none
integer,intent(in) :: nbas
integer,intent(in) :: aosource
character(*),intent(in) :: intfile,sortfile
integer(8),intent(in) :: maxmem
integer(8),intent(out),optional :: outinfo
type(AOSorterData) :: srt
integer :: iunit
integer :: maxrep, naos(8), lbuf, nibuf, nbits, lenint4
integer :: nsk, nt(8), ntoff(8)
integer :: nints, INDX
integer,allocatable :: idx_buf(:)
double precision :: val
double precision,allocatable :: val_buf(:)!, mat(:)
integer :: idx_p, idx_q, idx_r, idx_s, pq, rs, idx_end
integer :: sym_p, sym_q, sym_r, sym_s, sym_rs
integer :: i

call srt%open('TMPSORT',nbas,maxmem)

select case(aosource)

case(1) ! DALTON

   open(newunit=iunit,file=intfile,status='OLD',&
        access='SEQUENTIAL',form='UNFORMATTED')

   ! read info
   call readlabel(iunit,'BASINFO ')
   read(iunit) maxrep, naos, lbuf, nibuf, nbits, lenint4 

!   write(LOUT,'()')
!   write(LOUT,'(1x, a)') 'Dalton two-el. file initialized'
!   write(LOUT,'(1x,4(a,i5))') &
!        'Buffer size: ', lbuf, &
!        ', integers per index packet: ', nibuf, &
!        ', bits: ', nbits, &
!        ', record lenght in INT*4: ', lenint4

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
            !write(*,*) idx_p,idx_q,idx_r,idx_s 

            pq = idx_p + idx_q*(idx_q-1)/2
            rs = idx_r + idx_s*(idx_s-1)/2 
            call srt%add(pq,rs,val_buf(i))
            call srt%add(rs,pq,val_buf(i))

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
            call srt%add(pq,rs,val_buf(i))
            call srt%add(rs,pq,val_buf(i))

         enddo
      enddo

   end select

   deallocate(idx_buf)
   deallocate(val_buf)

   close(unit=iunit)

case(2) ! MOLPRO

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
                  call srt%add(pq,rs,val_buf(INDX))
                  call srt%add(rs,pq,val_buf(INDX))
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
                     call srt%add(pq,rs,val_buf(INDX))
                     call srt%add(rs,pq,val_buf(INDX))
                  endif

                  INDX=INDX+1
                  if(val_buf(INDX)/=0) then
                     pq = idx_p + idx_r*(idx_r-1)/2
                     rs = idx_q + idx_s*(idx_s-1)/2
                     call srt%add(pq,rs,val_buf(INDX))
                     call srt%add(rs,pq,val_buf(INDX))
                  endif

                  INDX=INDX+1
                  if(val_buf(INDX)/=0) then
                     pq = idx_q + idx_r*(idx_r-1)/2
                     rs = idx_p + idx_s*(idx_s-1)/2
                     call srt%add(pq,rs,val_buf(INDX))
                     call srt%add(rs,pq,val_buf(INDX))
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
                        call srt%add(pq,rs,val_buf(INDX))
                        call srt%add(rs,pq,val_buf(INDX))
                     endif

                     INDX=INDX+1
                     if(val_buf(INDX)/=0) then
                        pq = idx_p + idx_r*(idx_r-1)/2
                        rs = idx_q + idx_s*(idx_s-1)/2
                        call srt%add(pq,rs,val_buf(INDX))
                        call srt%add(rs,pq,val_buf(INDX))
                     endif

                     INDX=INDX+1
                     if(val_buf(INDX)/=0) then
                        pq = idx_q + idx_r*(idx_r-1)/2
                        rs = idx_p + idx_s*(idx_s-1)/2
                        call srt%add(pq,rs,val_buf(INDX))
                        call srt%add(rs,pq,val_buf(INDX))
                     endif

                  enddo
               enddo

            enddo

         enddo
      enddo
   enddo

   deallocate(val_buf)

   close(iunit)

case(3) ! Eugene

   open(newunit=iunit,file=intfile,status='OLD',&
        access='STREAM',form='UNFORMATTED')

   do
      read(iunit) val,idx_p,idx_q,idx_r,idx_s
      if(idx_r+idx_s==0) exit

      pq = min(idx_p,idx_q) + max(idx_p,idx_q)*(max(idx_p,idx_q)-1)/2
      rs = min(idx_r,idx_s) + max(idx_r,idx_s)*(max(idx_r,idx_s)-1)/2
      call srt%add(pq,rs,val)
      call srt%add(rs,pq,val)

   enddo

   ! find position in file
   if(present(outinfo)) then
      inquire(unit=iunit,pos=outinfo)
      outinfo = outinfo-24
   endif

   close(iunit)

case(4) ! Libor

   open(newunit=iunit,file=intfile,status='OLD',&
        access='STREAM',form='UNFORMATTED')
   read(iunit)idx_p, idx_q, idx_r

   pq=0
   do idx_p=1,nbas
      do idx_q=1,idx_p
         pq=pq+1

         rs=0
         do idx_r=1,nbas
            do idx_s=1,idx_r
               rs=rs+1

               if(pq.Ge.rs) then
                  read(iunit) val
                  if(val/=0) then
                     call srt%add(pq,rs,val)
                     call srt%add(rs,pq,val)
                  endif
               endif

            enddo
         enddo

      enddo
   enddo

   ! find position in file
   if(present(outinfo)) then
      inquire(unit=iunit,pos=outinfo)
      outinfo = outinfo-24
   endif

   close(iunit)

case default

   write(LOUT,'(a)') 'Unknown AOSOURCE in readtwoint'
   stop

end select

call srt%close(sortfile)

end subroutine readtwoint

subroutine open_AOSorter(srt,tmpname,nbas,maxmem)
implicit none
class(AOSorterData) :: srt
character(*),intent(in) :: tmpname
integer,intent(in) :: nbas
integer(8),intent(in) :: maxmem
integer(8) :: availmem1,availmem2,pass1_bytes,pass2_bytes,col_bytes,bck_bytes
integer :: ibck

srt%nbas = nbas
srt%nelm = (srt%nbas*(srt%nbas+1))/2

availmem1 = maxmem - (8 + 4    )*srt%nelm
availmem2 = maxmem - (8 + 4 + 8)*srt%nelm

pass1_bytes = 2*4 + 8
pass2_bytes =   4 + 8

col_bytes = pass2_bytes*srt%nelm
srt%ncol = min( int(availmem2/col_bytes) , srt%nelm )

srt%nbck = (srt%nelm-1)/srt%ncol+1

srt%seg_size  = srt%ncol
srt%seg_size  = srt%nelm*srt%seg_size

bck_bytes = pass1_bytes*srt%nbck
srt%bck_size = min( availmem1/bck_bytes , srt%seg_size )

srt%maxdmp = int((srt%seg_size-1)/srt%bck_size+1)

write(LOUT,'()')
write(LOUT,'(1x,a,t36,i18)') 'Number of basis functions:',        srt%nbas
write(LOUT,'(1x,a,t36,i18)') 'Number of independent pairs:',      srt%nelm
write(LOUT,'(1x,a,t36,i18)') 'Number of columns in segment:',     srt%ncol
write(LOUT,'(1x,a,t36,i18)') 'Number of segments:',               srt%nbck
write(LOUT,'(1x,a,t36,i18)') 'Size of segment:',                  srt%seg_size
write(LOUT,'(1x,a,t36,i18)') 'Size of bucket:',                   srt%bck_size
write(LOUT,'(1x,a,t36,i18)') 'Maximal number of writes/segment:', srt%maxdmp

allocate(srt%rs_diag(srt%nelm))
allocate(srt%rs_nelm(srt%nelm))
srt%rs_diag = 0
srt%rs_nelm = 0

allocate(srt%bck(srt%nbck))
do ibck=1,srt%nbck
   associate(bck => srt%bck(ibck))
     bck%n    = 0
     bck%ntot = 0
     bck%ndmp = 0
     nullify(bck%first,bck%last)
     allocate(bck%pqrs(2,srt%bck_size))
     allocate(bck%val(srt%bck_size))
   end associate
enddo

open(newunit=srt%unit,file=tmpname,access='stream',status='replace')

end subroutine open_AOSorter

subroutine add_AOSorter(srt,pq,rs,val)
implicit none
class(AOSorterData) :: srt
integer,intent(in) :: pq,rs
double precision,intent(in) :: val

if(pq==rs) then

   srt%rs_diag(rs) = val

else

   srt%rs_nelm(rs) = srt%rs_nelm(rs) + 1

   associate(bck => srt%bck((rs-1)/srt%ncol+1))

     if(bck%n==srt%bck_size) call bck%dump(srt%unit)
     bck%n = bck%n + 1
     bck%pqrs(1,bck%n) = pq
     bck%pqrs(2,bck%n) = rs
     bck%val(bck%n) = val

   end associate

endif

end subroutine add_AOSorter

subroutine close_AOSorter(srt,srtname)
implicit none
class(AOSorterData) :: srt
character(*),intent(in) :: srtname
integer :: unit
integer :: ibck,idmp
integer :: rs,rs_str,rs_end
integer(8) :: max_seg_size,max_bck_size,i
integer(8),allocatable :: seg_str(:)
integer(8) :: pos_diag,pos_nelm,pos_pos

do ibck=srt%nbck,1,-1
   associate(bck => srt%bck(ibck))
     if(bck%n>0) call bck%dump(srt%unit)
     deallocate(bck%val)
     deallocate(bck%pqrs)
   end associate
enddo

max_seg_size = 0
max_bck_size = 0

do ibck=1,srt%nbck
   associate(bck => srt%bck(ibck))
     allocate(bck%dmp(bck%ndmp))
     do idmp=1,bck%ndmp
        bck%dmp(idmp) = bck%first%DumpData
        bck%last => bck%first%next
        deallocate(bck%first)
        bck%first => bck%last
     enddo
     max_seg_size = max(max_seg_size,bck%ntot)
     max_bck_size = max(max_bck_size,maxval(bck%dmp%n))
   end associate
enddo
max_seg_size = max_seg_size + srt%ncol

allocate(srt%seg_pq(max_seg_size))
allocate(srt%seg_val(max_seg_size))
allocate(srt%pqrs(2,max_bck_size))
allocate(srt%val(max_bck_size))
allocate(srt%rs_pos(srt%nelm))

srt%rs_pos = 0
pos_diag = 0
pos_nelm = 0
pos_pos  = 0

open(newunit=unit,file=srtname,access='stream',status='replace')
write(unit,pos=1) srt%nbas,srt%nelm,pos_diag,pos_nelm,pos_pos
inquire(unit,pos=pos_diag)
write(unit) srt%rs_diag
inquire(unit,pos=pos_nelm)
write(unit) srt%rs_nelm
inquire(unit,pos=pos_pos)
write(unit) srt%rs_pos

do ibck=1,srt%nbck
   rs_str = (ibck-1)*srt%ncol + 1
   rs_end = min(ibck*srt%ncol,srt%nelm)
   allocate(seg_str(rs_str:rs_end))
   seg_str(rs_str) = 1
   do rs=rs_str+1,rs_end
      seg_str(rs) = seg_str(rs-1) + srt%rs_nelm(rs-1) + 1
   enddo
   srt%rs_nelm(rs_str:rs_end) = 0

   associate(bck => srt%bck(ibck))
     do idmp=1,bck%ndmp
        associate(dmp => bck%dmp(idmp))
          read(srt%unit,pos=dmp%pos) srt%pqrs(:,1:dmp%n),srt%val(1:dmp%n)
          do i=1,dmp%n
             rs = srt%pqrs(2,i)
             call insert_seg(srt%pqrs(1,i),srt%val(i),&
                  srt%rs_nelm(rs),&
                  srt%seg_pq(seg_str(rs):),srt%seg_val(seg_str(rs):))
          enddo
        end associate
     enddo
   end associate
   do rs=rs_str,rs_end
      if(abs(srt%rs_diag(rs))>0) then
         call insert_seg(rs,srt%rs_diag(rs),&
              srt%rs_nelm(rs),&
              srt%seg_pq(seg_str(rs):),srt%seg_val(seg_str(rs):))
      endif
   enddo

   do rs=rs_str,rs_end
      if(srt%rs_nelm(rs)>0) then
         inquire(unit,pos=srt%rs_pos(rs))
         write(unit) &
              srt%seg_pq(seg_str(rs):seg_str(rs)+srt%rs_nelm(rs)-1),&
              srt%seg_val(seg_str(rs):seg_str(rs)+srt%rs_nelm(rs)-1)
      endif
   enddo

   deallocate(seg_str)
enddo

write(unit,pos=1) srt%nbas,srt%nelm,pos_diag,pos_nelm,pos_pos
write(unit,pos=pos_diag) srt%rs_diag
write(unit,pos=pos_nelm) srt%rs_nelm
write(unit,pos=pos_pos) srt%rs_pos
close(unit,status='keep')

deallocate(srt%rs_pos)
deallocate(srt%val)
deallocate(srt%pqrs)
deallocate(srt%seg_val)
deallocate(srt%seg_pq)

do ibck=1,srt%nbck
   associate(bck => srt%bck(ibck))
     deallocate(bck%dmp)
   end associate
enddo
deallocate(srt%bck)

deallocate(srt%rs_nelm)
deallocate(srt%rs_diag)

close(srt%unit,status='delete')

end subroutine close_AOSorter

subroutine insert_seg(pq,val,nelm,seg_pq,seg_val)
implicit none
integer,intent(in) :: pq
double precision,intent(in) :: val
integer,intent(inout) :: nelm
integer,intent(inout) :: seg_pq(:)
double precision,intent(inout) :: seg_val(:)
integer :: i,j

!simple

nelm = nelm + 1

seg_pq(nelm)  = pq
seg_val(nelm) = val

!-------------------------------------------------------------------------------

!full sort

!!$do i=nelm,1,-1
!!$   if(seg_pq(i)<=pq) then
!!$      if(seg_pq(i)==pq) return
!!$      exit
!!$   endif
!!$enddo
!!$j = i + 1
!!$
!!$do i=nelm,j,-1
!!$   seg_pq(i+1) = seg_pq(i)
!!$enddo
!!$seg_pq(j) = pq
!!$do i=nelm,j,-1
!!$   seg_val(i+1) = seg_val(i)
!!$enddo
!!$seg_val(j) = val
!!$
!!$nelm = nelm + 1

!-------------------------------------------------------------------------------

end subroutine insert_seg

subroutine dump_AOBucket(bck,unit)
implicit none
class(AOBucketData) :: bck
integer,intent(in) :: unit

if(bck%ndmp==0) then
   allocate(bck%first)
   bck%last => bck%first
else
   allocate(bck%last%next)
   bck%last => bck%last%next
endif
nullify(bck%last%next)

associate(dmp => bck%last)
  dmp%n = bck%n
  inquire(unit,pos=dmp%pos)
  write(unit) bck%pqrs(:,1:dmp%n),bck%val(1:dmp%n)
end associate

bck%ntot = bck%ntot + bck%n
bck%n = 0

bck%ndmp = bck%ndmp + 1

end subroutine dump_AOBucket

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine open_AOReader(rdr,srtname)
implicit none
class(AOReaderData) :: rdr
character(*),intent(in) :: srtname
integer :: max_pq

open(newunit=rdr%unit,file=srtname,access='stream',status='old')
read(rdr%unit) rdr%nbas,rdr%nelm,rdr%pos_diag,rdr%pos_nelm,rdr%pos_pos

allocate(rdr%rs_nelm(rdr%nelm))
allocate(rdr%rs_pos(rdr%nelm))

read(rdr%unit,pos=rdr%pos_nelm) rdr%rs_nelm
read(rdr%unit,pos=rdr%pos_pos)  rdr%rs_pos

max_pq = maxval(rdr%rs_nelm)
allocate(rdr%pq(max_pq))
allocate(rdr%val(max_pq))

end subroutine open_AOReader

subroutine get_AOReader(rdr,rs,ints)
implicit none
class(AOReaderData) :: rdr
integer,intent(in) :: rs
double precision,intent(out) :: ints(:)
integer :: nelm,i

nelm = rdr%rs_nelm(rs)
ints(1:rdr%nelm) = 0

if(nelm==0) return

read(rdr%unit,pos=rdr%rs_pos(rs)) rdr%pq(1:nelm),rdr%val(1:nelm)
do i=1,nelm
   ints(rdr%pq(i)) = rdr%val(i)
enddo

end subroutine get_AOReader

subroutine close_AOReader(rdr)
implicit none
class(AOReaderData) :: rdr

deallocate(rdr%val)
deallocate(rdr%pq)

deallocate(rdr%rs_pos)
deallocate(rdr%rs_nelm)

close(rdr%unit)

end subroutine close_AOReader

end module sorter
