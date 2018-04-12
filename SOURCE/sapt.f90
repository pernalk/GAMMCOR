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
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
double precision,allocatable :: OmA(:),OmB(:)
!double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: TwoMO(:)
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

! read EigValA_B
! allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX),&
!          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX),&
!          TwoMO(NInte2))

 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX),&
          TwoMO(NInte2))

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

! read 2-el integrals
  call LoadSaptTwoEl(3,TwoMO,NBas,NInte2)

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

! HERE!!!!
! uncoupled

e2du=0d0
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)
   dea = A%OrbE(ip)-A%OrbE(iq) 
   do rs=1,B%NDimX
      ir = B%IndN(1,rs)
      is = B%IndN(2,rs)
      deb = B%OrbE(ir)-B%OrbE(is) 
      e2du = e2du + & 
            TwoMO(iaddr(ip,iq,ir,is))**2d0/(dea+deb)
       enddo
       !write(*,*) 'tmp1',tmp1(i,pq)
    enddo
 write(*,*)'e2du: ', -4d0*e2du 

! 2-el 
!tmp=0d0
!do pq=1,A%NDimX
!   ip = A%IndN(1,pq)
!   iq = A%IndN(2,pq)
!   do rs=1,B%NDimX
!      ir = B%IndN(1,rs)
!      is = B%IndN(2,rs)
!      tmp = tmp + & 
!            TwoMO(iaddr(ip,iq,ir,is))**2d0
!       enddo
!       !write(*,*) 'tmp1',tmp1(i,pq)
!    enddo
! write(*,*)'TwoEl:', tmp
! tmp=0d0 

! do i=1,A%NDimX
!    write(*,*) A%IndN(1,i)
!    write(*,*) A%IndN(2,i)
! enddo

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))

!do j=1,A%NDimX
!   ir = A%IndN(1,j)
!   is = A%IndN(2,j)
!   write(*,*) is,ir, sqrt(A%Occ(ir)), A%CICoef(ir)
!enddo

!do i=1,A%NDimX
!   ip = A%IndN(1,i)
!   iq = A%IndN(2,i)
!   write(*,*) A%Occ(iq),A%Occ(ip)
!enddo

! coupled - 1
! tmp1=0
! do i=1,A%NDimX
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       do rs=1,B%NDimX
!          ir = B%IndN(1,rs)
!          is = B%IndN(2,rs)
!          tmp1(i,rs) = tmp1(i,rs) + & 
!                       (A%CICoef(iq)+A%CICoef(ip)) * &
!                       (B%CICoef(is)+B%CICoef(ir)) * &
!                      ! EVecA(i,pq)*TwoMO(iaddr(ip,iq,ir,is))
!                       EVecA(pq,i)*TwoMO(iaddr(ip,iq,ir,is))
!       enddo
!       !write(*,*) 'tmp1',tmp1(i,pq)
!    enddo
! enddo
! tmp2=0
! do j=1,B%NDimX
!    do i=1,A%NDimX
!       do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!       tmp2(i,j) = tmp2(i,j) + &
!                   EVecB(j,rs)*tmp1(i,rs)
!       enddo
!    enddo   
! enddo

! e2d = 0d0
! do i=1,A%NDimX
!    do j=1,B%NDimX
!       e2d = e2d + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
!    enddo
! enddo
! e2d = -16d0*e2d
!!
! print*, 'e2disp: ',e2d

! coupled -2 
! e2d = 0d0
! tmp=0d0
! do i=1,A%NDimX
!    do j=1,B%NDimX
!
!       tmp = 0d0
!       do pq=1,A%NDimX
!          ip = A%IndN(1,pq)
!          iq = A%IndN(2,pq)
! 
!          do rs=1,B%NDimX
!             ir = B%IndN(1,rs)
!             is = B%IndN(2,rs)
!             tmp = tmp - &
!                   (A%CICoef(ip)+A%CICoef(iq)) * &
!                   (B%CICoef(ir)+B%CICoef(is)) * &
!                   EVecA(pq,i)*TwoMO(iaddr(ip,iq,ir,is))*&
!                   EVecB(j,rs)
!          enddo
!       enddo
!
!       e2d = e2d  + tmp**2d0/(OmA(i)+OmB(j))
!    enddo
! enddo
! e2d = -16d0*e2d
! 
! print*, 'e2disp: ',e2d*1000

! coupled-3
 e2d = 0d0
 tmp=0d0
 do i=1,A%NDimX
    do j=1,B%NDimX
!
       do pq=1,A%NDimX
          ip = A%IndN(1,pq)
          iq = A%IndN(2,pq)
 
          do rs=1,B%NDimX
             ir = B%IndN(1,rs)
             is = B%IndN(2,rs)
             tmp = tmp + &
                   (A%CICoef(ip)+A%CICoef(iq)) * &
                   (B%CICoef(ir)+B%CICoef(is)) * &
                   EVecA((pq-1)*A%NDimX+i)*TwoMO(iaddr(ip,iq,ir,is))*&
                   EVecB((j-1)*B%NDimX+rs)
          enddo
       enddo

       e2d = e2d  + tmp**2d0/(OmA(i)+OmB(j))
       tmp=0
    enddo
 enddo
 e2d = -16d0*e2d
! 
 print*, 'e2disp: ',e2d*1000



 deallocate(EVecA,EVecB,OmA,OmB)
 deallocate(TwoMO,tmp1,tmp2)

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
