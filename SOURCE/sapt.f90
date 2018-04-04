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
!double precision,allocatable :: EValA(:),EValB(:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: TwoMO(:)
double precision :: e2d,tmp

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
 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX),&
          TwoMO(NInte2))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')
! TEST RESPONSE
 tmp = 0
!do j=1,A%NDimX
!   do i=1,A%NDimX
!!     write(*,*) 'evecA', EVecA(i,i)
!!     tmp = tmp + EVecA(i,j)**2
!    enddo 
!enddo
!do j=1,B%NDimX
 !  do i=1,B%NDimX
     !write(*,*) 'evecB', EValB(i)
!     tmp = tmp + EVecB(i,j)**2
  !  enddo 
!enddo

 print*, tmp

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

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))

!do j=1,B%NDimX
!   write(*,*) B%IndN(1,j), B%IndN(2,j)
!enddo
 tmp1=0
 do i=1,A%NDimX
    do pq=1,A%NDimX
       ip = A%IndN(1,pq)
       iq = A%IndN(2,pq)
       do rs=1,B%NDimX
          ir = B%IndN(1,rs)
          is = B%IndN(2,rs)
          tmp1(i,rs) = tmp1(i,rs) - & 
                       sqrt(A%Occ(iq)-A%Occ(ip)) * &
                       sqrt(B%Occ(is)-B%Occ(ir)) * &
                       !EVecA(i,pq)*TwoMO(iaddr(ip,iq,ir,is))
                       EVecA(pq,i)*TwoMO(iaddr(ip,iq,ir,is))
       enddo
!       write(*,*) 'tmp1',tmp1(i,rs)
    enddo
 enddo
 tmp2=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j)  + &
                   EVecB(j,rs)*tmp1(i,rs)
       enddo
    enddo   
 enddo

 e2d = 0d0
 do i=1,A%NDimX
    do j=1,B%NDimX
       e2d = e2d + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
    enddo
 enddo
 e2d = -16d0*e2d

 print*, 'e2disp: ',e2d

 e2d = 0d0
 do i=1,A%NDimX
    do j=1,B%NDimX

       tmp = 0d0
       do pq=1,A%NDimX
          ip = A%IndN(1,pq)
          iq = A%IndN(2,pq)
 
          do rs=1,B%NDimX
             ir = B%IndN(1,rs)
             is = B%IndN(2,rs)
             tmp = tmp - &
                   sqrt(A%Occ(iq)-A%Occ(ip))* &
                   sqrt(B%Occ(is)-B%Occ(ir))* & 
                   EVecA(pq,i)*TwoMO(iaddr(ip,iq,ir,is))*&
                   EVecB(j,rs)
          enddo
       enddo

       e2d = e2d  + tmp**2d0/(OmA(i)+OmB(j))
    enddo
 enddo
 e2d = -16d0*e2d

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
