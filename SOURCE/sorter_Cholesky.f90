module sorter_Cholesky
!
! this module is used by Cholesky
! which avoids creation of AOTWOSORT
!
use types
implicit none

type AOReaderChol
integer :: unit
integer :: aosource
integer :: nbas,nelm,nibuf
integer,allocatable :: idx_buf(:)
double precision,allocatable :: val_buf(:)
contains
procedure :: open => open_AOReaderChol
procedure :: get_TR => get_TR_AOReaderChol
procedure :: close => close_AOReaderChol
end type

contains

subroutine add_TR(pq,rs,val,y,M)
implicit none
integer,intent(in) :: pq,rs
integer,intent(in) :: y(:)
double precision,intent(in)  :: val
double precision,intent(out) :: M(:,:)

if(y(rs) > 0) M(pq, y(rs)) = val
if(y(pq) > 0) M(rs, y(pq)) = val

end subroutine add_TR

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine open_AOReaderChol(rdr,nbas,aosource,file_path,rs_diag)
! reads the diagonal
implicit none
class(AOReaderChol) :: rdr
integer,intent(in)  :: aosource,nbas
character(*), intent(in) :: file_path
double precision,intent(out) :: rs_diag(:)

double precision :: val
integer :: nsk, nt(8), ntoff(8)
integer :: maxrep, naos(8), lbuf, nibuf, nbits, lenint4
integer :: nints, INDX
integer :: idx_p, idx_q, idx_r, idx_s, pq, rs, idx_end
integer :: sym_p, sym_q, sym_r, sym_s, sym_rs
integer :: i

rdr%aosource = aosource
rdr%nbas = nbas
rdr%nelm = nbas*(nbas+1)/2
rs_diag = 0

select case(rdr%aosource)
case(1) ! DALTON

   open(newunit=rdr%unit,file=file_path,status='OLD',&
        access='SEQUENTIAL',form='UNFORMATTED')

   ! read info
   call readlabel(rdr%unit,'BASINFO ')
   read(rdr%unit) maxrep, naos, lbuf, nibuf, nbits, lenint4

   rdr%nibuf = nibuf

   allocate(rdr%val_buf(lbuf))
   allocate(rdr%idx_buf(nibuf*lbuf))

   call readlabel(rdr%unit,'BASTWOEL')

      select case(rdr%nibuf)
      case(1)

         do
            read(rdr%unit) rdr%val_buf, rdr%idx_buf, nints
            if(nints<0) exit
            do i=1,nints

               INDX = rdr%idx_buf(i)
               idx_p = ibits(INDX,0,8)
               idx_q = ibits(INDX,8,8)
               idx_r = ibits(INDX,16,8)
               idx_s = ibits(INDX,24,8)

               pq = idx_p + idx_q*(idx_q-1)/2
               rs = idx_r + idx_s*(idx_s-1)/2
               if (pq == rs) rs_diag(rs) = rdr%val_buf(i)

            enddo
         enddo

      case(2)

         do
            read(rdr%unit) rdr%val_buf, rdr%idx_buf, nints
            if(nints<0) exit
            do i=1,nints

               INDX = rdr%idx_buf(2*i-1)
               idx_r = ibits(INDX,0,16)
               idx_s = ibits(INDX,16,16)
               INDX = rdr%idx_buf(2*i)
               idx_p = ibits(INDX,0,16)
               idx_q = ibits(INDX,16,16)

               pq = idx_p + idx_q*(idx_q-1)/2
               rs = idx_r + idx_s*(idx_s-1)/2
               if (pq == rs) rs_diag(rs) = rdr%val_buf(i)

            enddo
         enddo

      end select

case(2) ! MOLPRO

   open(newunit=rdr%unit,file=file_path,status='OLD',&
        access='SEQUENTIAL',form='UNFORMATTED')

   ! read info
   read(rdr%unit) nsk
   read(rdr%unit) nt

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
   allocate(rdr%val_buf(lbuf))

      do sym_p=1,nsk
         if(nt(sym_p)==0) cycle
         read(rdr%unit)

         do
            read(rdr%unit) idx_s,idx_r
            if(idx_s<=0.and.idx_r<=0) exit
            nints = idx_r + idx_s*(idx_s-1)/2
            read(rdr%unit) rdr%val_buf(1:nints)

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
                  if(rdr%val_buf(INDX)/=0) then
                     pq = idx_p + idx_q*(idx_q-1)/2
                     if (pq == rs) rs_diag(rs) = rdr%val_buf(INDX)
                  endif

               enddo
            enddo

         enddo

      enddo

      do sym_r=2,nsk
         do sym_p=1,sym_r-1
            if(nt(sym_p)==0) cycle
            if(nt(sym_r)==0) cycle
            read(rdr%unit)

            nints=(nt(sym_p)*(nt(sym_p)+1))/2
            nints=3*nints
            do
               read(rdr%unit) idx_s,idx_r
               if(idx_s<=0.and.idx_r<=0) exit
               read(rdr%unit) rdr%val_buf(1:nints)

               if(idx_s/=idx_r) cycle

               idx_r = idx_r + ntoff(sym_r)
               idx_s = idx_s + ntoff(sym_r)

               INDX=0
               do idx_q=ntoff(sym_p)+1,ntoff(sym_p)+nt(sym_p)
                  do idx_p=ntoff(sym_p)+1,idx_q

                     INDX=INDX+1

                     INDX=INDX+1
                     if(rdr%val_buf(INDX)/=0) then
                        pq = idx_p + idx_r*(idx_r-1)/2
                        rs = idx_q + idx_s*(idx_s-1)/2
                        if (pq == rs) rs_diag(rs) = rdr%val_buf(INDX)
                     endif

                     INDX=INDX+1

                  enddo
               enddo

            enddo

         enddo
      enddo

case(4) ! Libor

      open(newunit=rdr%unit,file=file_path,status='OLD',&
           access='STREAM',form='UNFORMATTED')
      read(rdr%unit,pos=1) idx_p, idx_q, idx_r

      pq=0
      do idx_p=1,nbas
         do idx_q=1,idx_p
            pq=pq+1

            rs=0
            do idx_r=1,nbas
               do idx_s=1,idx_r
                  rs=rs+1

                  if(pq.ge.rs) then
                     read(rdr%unit) val
                     if(pq ==rs) rs_diag(rs) = val
                  endif

               enddo
            enddo

         enddo
      enddo

case default

   write(LOUT,'(a)') 'Unknown AOSOURCE in open_AOReaderChol'
   stop

end select

end subroutine open_AOReaderChol

subroutine get_TR_AOReaderChol(rdr,y,M)
! reader with adder_TR
implicit none
class(AOReaderChol) :: rdr
integer,intent(in)  :: y(:)
double precision,intent(out) :: M(:,:) 

double precision :: val
integer :: lbuf
integer :: nsk, nt(8), ntoff(8)
integer :: nints, INDX
integer :: idx_p, idx_q, idx_r, idx_s, pq, rs, idx_end
integer :: sym_p, sym_q, sym_r, sym_s, sym_rs
integer :: i

M = 0

select case(rdr%aosource)
case(1) ! DALTON

   ! read info
   call readlabel(rdr%unit,'BASTWOEL')

      select case(rdr%nibuf)
      case(1)

         do
            read(rdr%unit) rdr%val_buf, rdr%idx_buf, nints
            if(nints<0) exit
            do i=1,nints

               INDX = rdr%idx_buf(i)
               idx_p = ibits(INDX,0,8)
               idx_q = ibits(INDX,8,8)
               idx_r = ibits(INDX,16,8)
               idx_s = ibits(INDX,24,8)

               pq = idx_p + idx_q*(idx_q-1)/2
               rs = idx_r + idx_s*(idx_s-1)/2
               call add_TR(pq,rs,rdr%val_buf(i),y,M)

            enddo
         enddo

      case(2)

         do
            read(rdr%unit) rdr%val_buf, rdr%idx_buf, nints
            if(nints<0) exit
            do i=1,nints

               INDX = rdr%idx_buf(2*i-1)
               idx_r = ibits(INDX,0,16)
               idx_s = ibits(INDX,16,16)
               INDX = rdr%idx_buf(2*i)
               idx_p = ibits(INDX,0,16)
               idx_q = ibits(INDX,16,16)

               pq = idx_p + idx_q*(idx_q-1)/2
               rs = idx_r + idx_s*(idx_s-1)/2
               call add_TR(pq,rs,rdr%val_buf(i),y,M)

            enddo
         enddo

      end select

case(2) ! MOLPRO

   ! read info
   rewind(rdr%unit)
   read(rdr%unit) nsk
   read(rdr%unit) nt

   ! offset
   ntoff = 0
   do i=2,8
      ntoff(i) = ntoff(i-1)+nt(i-1)
   enddo

   lbuf = 3*maxval(nt(1:nsk))**2

      do sym_p=1,nsk
         if(nt(sym_p)==0) cycle
         read(rdr%unit)

         do
            read(rdr%unit) idx_s,idx_r
            if(idx_s<=0.and.idx_r<=0) exit
            nints = idx_r + idx_s*(idx_s-1)/2
            read(rdr%unit) rdr%val_buf(1:nints)

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
                  if(rdr%val_buf(INDX)/=0) then
                     pq = idx_p + idx_q*(idx_q-1)/2
                     call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                  endif

               enddo
            enddo

         enddo

      enddo

      do sym_r=2,nsk
         do sym_p=1,sym_r-1
            if(nt(sym_p)==0) cycle
            if(nt(sym_r)==0) cycle
            read(rdr%unit)

            nints=(nt(sym_p)*(nt(sym_p)+1))/2
            nints=3*nints
            do
               read(rdr%unit) idx_s,idx_r
               if(idx_s<=0.and.idx_r<=0) exit
               read(rdr%unit) rdr%val_buf(1:nints)

               idx_r = idx_r + ntoff(sym_r)
               idx_s = idx_s + ntoff(sym_r)

               INDX=0
               do idx_q=ntoff(sym_p)+1,ntoff(sym_p)+nt(sym_p)
                  do idx_p=ntoff(sym_p)+1,idx_q

                     INDX=INDX+1
                     if(rdr%val_buf(INDX)/=0) then
                        pq = idx_p + idx_q*(idx_q-1)/2
                        rs = idx_r + idx_s*(idx_s-1)/2
                        call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                     endif

                     INDX=INDX+1
                     if(rdr%val_buf(INDX)/=0) then
                        pq = idx_p + idx_r*(idx_r-1)/2
                        rs = idx_q + idx_s*(idx_s-1)/2
                        call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                     endif

                     INDX=INDX+1
                     if(rdr%val_buf(INDX)/=0) then
                        pq = idx_q + idx_r*(idx_r-1)/2
                        rs = idx_p + idx_s*(idx_s-1)/2
                        call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
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
               read(rdr%unit)

               nints=nt(sym_p)*nt(sym_q)
               nints=3*nints
               do
                  read(rdr%unit) idx_s,idx_r
                  if(idx_s<=0.and.idx_r<=0) exit
                  read(rdr%unit) rdr%val_buf(1:nints)

                  idx_r = idx_r + ntoff(sym_r)
                  idx_s = idx_s + ntoff(sym_s)

                  INDX=0
                  do idx_q=ntoff(sym_q)+1,ntoff(sym_q)+nt(sym_q)
                     do idx_p=ntoff(sym_p)+1,ntoff(sym_p)+nt(sym_p)

                        INDX=INDX+1
                        if(rdr%val_buf(INDX)/=0) then
                           pq = idx_p + idx_q*(idx_q-1)/2
                           rs = idx_r + idx_s*(idx_s-1)/2
                           call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                        endif

                        INDX=INDX+1
                        if(rdr%val_buf(INDX)/=0) then
                           pq = idx_p + idx_r*(idx_r-1)/2
                           rs = idx_q + idx_s*(idx_s-1)/2
                           call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                        endif

                        INDX=INDX+1
                        if(rdr%val_buf(INDX)/=0) then
                           pq = idx_q + idx_r*(idx_r-1)/2
                           rs = idx_p + idx_s*(idx_s-1)/2
                           call add_TR(pq,rs,rdr%val_buf(INDX),y,M)
                        endif

                     enddo
                  enddo

               enddo

            enddo
         enddo
      enddo

case(4) ! Libor

      read(rdr%unit,pos=1) idx_p, idx_q, idx_r

      pq=0
      do idx_p=1,rdr%nbas
         do idx_q=1,idx_p
            pq=pq+1

            rs=0
            do idx_r=1,rdr%nbas
               do idx_s=1,idx_r
                  rs=rs+1

                  if (pq .ge. rs) then
                     read(rdr%unit) val
                     call add_TR(pq,rs,val,y,M)
                  endif

               enddo
            enddo

         enddo
      enddo

case default

   write(LOUT,'(a)') 'Unknown AOSOURCE in get_TR_AOReaderChol'
   stop

end select

end subroutine get_TR_AOReaderChol

subroutine close_AOReaderChol(rdr)
implicit none
class(AOReaderChol) :: rdr

select case(rdr%aosource)
case(1) ! DALTON

   deallocate(rdr%idx_buf)
   deallocate(rdr%val_buf)

case(2) ! MOLPRO

   deallocate(rdr%val_buf)

end select

close(unit=rdr%unit)

end subroutine close_AOReaderChol

end module

