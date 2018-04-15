module sapt_ener
use types
!use sapt_main

implicit none

contains

subroutine e2disp(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
double precision,allocatable :: OmA(:),OmB(:)
!double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: e2d,tmp
double precision :: e2du,dea,deb

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 print*,NBas,NInte1,NInte2,nOVA,nOVB

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

! TEST RESPONSE
! tmp = 0
!do j=1,A%NDimX
!!   do i=1,A%NDimX
!!!!     write(*,*) 'evecA', EVecA(i,i)
!!     tmp = tmp + EVecA(i,j)**2
!!    enddo 
!     tmp = tmp + OmA(j)**2
!     write(*,*) OmA(j)
!enddo
! print*, 'test-resp:',tmp
! tmp = 0
! do j=1,B%NDimX
!!  do i=1,B%NDimX
!!    !write(*,*) 'evecB', EValB(i)
!!     tmp = tmp + EVecB(i,j)**2
!!   enddo 
!     tmp = tmp + OmB(j)**2
!     write(*,*) OmB(j)
! enddo
! print*, 'test-resp:',tmp

! print*, 'iaddr'
! do is=1,NBas
!    do ir=1,is
!       do iq=1,NBas
!          do ip=1,iq
!             write(*,*) ip,iq,ir,is,TwoMO(iaddr(ip,iq,ir,is))
!           enddo
!       enddo
!    enddo
! enddo

! uncoupled
! works with tran4_full
!allocate(work(NInte1))
!open(newunit=iunit,file='TWOMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NInte1)
!
!
!e2du=0d0
!do pq=1,A%NDimX
!   ip = A%IndN(1,pq)
!   iq = A%IndN(2,pq)
!   dea = A%OrbE(ip)-A%OrbE(iq)
!   !print*, iq,ip,iq+ip*(ip-1)/2,NInte1
!   read(iunit,rec=iq+ip*(ip-1)/2) work(1:NInte1)
!   do rs=1,B%NDimX
!      ir = B%IndN(1,rs)
!      is = B%IndN(2,rs)
!      deb = B%OrbE(ir)-B%OrbE(is) 
!      !print*, is,ir,is+ir*(ir-1)/2,NInte1
!      e2du = e2du + work(is+ir*(ir-1)/2)**2/(dea+deb)
!   enddo
!enddo
! write(LOUT,'()') 
! write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0 

! uncoupled
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

e2du=0d0
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)
   dea = A%OrbE(ip)-A%OrbE(iq)
!   print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
   do rs=1,B%NDimX
      ir = B%IndN(1,rs)
      is = B%IndN(2,rs)
      deb = B%OrbE(ir)-B%OrbE(is) 
!      print*, is,ir,is+(ir-B%num0-1)*dimOB,nOVB
      e2du = e2du + work(is+(ir-B%num0-1)*dimOB)**2/(dea+deb)
   enddo
enddo
 write(LOUT,'()') 
 write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0 


allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))

! coupled - 1
! works with tran4_full
! tmp1=0
! do i=1,A%NDimX
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       read(iunit,rec=iq+ip*(ip-1)/2) work(1:NInte1)
!       do rs=1,B%NDimX
!          ir = B%IndN(1,rs)
!          is = B%IndN(2,rs)
!          tmp1(i,rs) = tmp1(i,rs) + & 
!                       (A%CICoef(iq)+A%CICoef(ip)) * &
!                       (B%CICoef(is)+B%CICoef(ir)) * &
!                       EVecA((i-1)*A%NDimX+pq)*work(is+ir*(ir-1)/2)
!       enddo
!    enddo
! enddo
! tmp2=0
! do j=1,B%NDimX
!    do i=1,A%NDimX
!       do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!       tmp2(i,j) = tmp2(i,j) + &
!                    EVecB((j-1)*B%NDimX+rs)*tmp1(i,rs)
!       enddo
!    enddo   
! enddo
!
! e2d = 0d0
! do i=1,A%NDimX
!    do j=1,B%NDimX
!       e2d = e2d + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
!    enddo
! enddo
! e2d = -16d0*e2d*1000d0
! write(LOUT,'(1x,a,5x,f16.8)') 'E2disp = ',e2d

!! coupled -2
! e2d = 0d0
! do i=1,A%NDimX
!    do j=1,B%NDimX
!
!       tmp=0d0
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
!          read(iunit,rec=iq+ip*(ip-1)/2) work(1:NInte1)
!          do rs=1,B%NDimX
!             ir = B%IndN(1,rs)
!             is = B%IndN(2,rs)
!             tmp = tmp + &
!                   (A%CICoef(ip)+A%CICoef(iq)) * &
!                   (B%CICoef(ir)+B%CICoef(is)) * &
!                   EVecA((i-1)*A%NDimX+pq)*work(is+ir*(ir-1)/2)*&
!                   EVecB((j-1)*B%NDimX+rs)
!          enddo
!       enddo
!
!       e2d = e2d  + tmp**2d0/(OmA(i)+OmB(j))
!    enddo
! enddo
! e2d = -16d0*e2d
! 
! print*, 'e2disp: ',e2d*1000

 close(iunit)
 deallocate(work)

 deallocate(tmp1,tmp2)
 deallocate(EVecA,EVecB,OmA,OmB)

end subroutine e2disp

subroutine readresp(EVec,EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVec(NDim,NDim), EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVec
 read(iunit) EVal

 close(iunit)

end subroutine readresp



end module sapt_ener
