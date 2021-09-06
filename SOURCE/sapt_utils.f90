module sapt_utils
use types
use tran
use diis

contains

subroutine read_SBlock(SBlock,SBlockIV,nblk,xy0file)
implicit none

type(EBlockData)             :: SBlockIV
type(EBlockData),allocatable :: SBlock(:)
integer,intent(out)          :: nblk
character(*),intent(in)      :: xy0file

integer                      :: iunit,i,ival,iblk

! read SBLOCK
open(newunit=iunit,file=xy0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
!print*, 'ival,nblk',ival,nblk
allocate(SBlock(nblk))
do iblk=1,nblk
   associate(A => SBlock(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n),A%matY(A%n,A%n),A%vec(A%n))
     read(iunit) A%pos,A%matX,A%matY,A%vec
   end associate
enddo
associate(A => SBlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  allocate(A%matX(A%n,1),A%matY(A%n,1))
  read(iunit) A%pos,A%vec
  read(iunit) A%matX,A%matY
end associate
close(iunit)

end subroutine read_SBlock

subroutine read_ABPM0Block(A0Block,A0BlockIV,nblk,abpm0file)
implicit none

type(EBlockData)             :: A0BlockIV
type(EBlockData),allocatable :: A0Block(:)
integer,intent(out)          :: nblk
character(*),intent(in)      :: abpm0file

integer                      :: iunit,i,ival,iblk

! read ABPM0 = ABPLUS0.ABMIN0
open(newunit=iunit,file=abpm0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
allocate(A0Block(nblk))
do iblk=1,nblk
   associate(A => A0Block(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n))
     read(iunit) A%pos,A%matX
   end associate
enddo
associate(A => A0BlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  read(iunit) A%pos,A%vec
end associate
close(iunit)

end subroutine read_ABPM0Block

subroutine create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig,CICoef,IndN,nblk,NBasis,NDimX)
implicit none

type(EBlockData)   :: Sblock(nblk),SblockIV
type(Y01BlockData) :: Y01Block(NDimX)

integer,intent(in)           :: nblk,NBasis,NDimX
integer,intent(in)           :: IndN(2,NDimX)
double precision,intent(in)  :: CICoef(NBasis)
double precision,intent(out) :: Eig(NDimX)

integer          :: i,ii,ipos,iblk,ip,iq
double precision :: fact

do iblk=1,nblk
   associate(A => SBlock(iblk))

     do i=1,A%n
        ipos = A%pos(i)

        ip   = IndN(1,ipos)
        iq   = IndN(2,ipos)
        fact = CICoef(ip)-CICoef(iq)
        associate(Y => Y01Block(ipos))

          Y%n  = A%n
          Y%l1 = A%l1
          Y%l2 = A%l2
          allocate(Y%vec0(Y%n))
          Y%vec0(1:Y%n) = fact*(A%matY(i,1:Y%n)-A%matX(i,1:Y%n))

        end associate
     enddo
     Eig(A%l1:A%l2) = A%vec(1:A%n)
   end associate
enddo
associate(A => SBlockIV)

  do i=1,A%n
     ii = A%l1+i-1
     ipos = A%pos(i)

     associate(Y => Y01Block(ipos))

       Y%n = 1
       Y%l1 = ii
       Y%l2 = ii
       allocate(Y%vec0(1))
       Y%vec0(1) = 1d0/sqrt(2d0)

     end associate

     Eig(ii) = A%vec(i)
  enddo

end associate

end subroutine create_Y0tilde

subroutine convert_XY0_to_Y01(Mon,Y01Block,Eig0,NBasis,xy0file)
implicit none

type(SystemBlock)              :: Mon
type(EblockData)               :: SBlockIV
type(EblockData),allocatable   :: SBlock(:)
type(Y01BlockData)             :: Y01Block(Mon%NDimX)

integer,intent(in)           :: NBasis
character(*),intent(in)      :: xy0file
double precision,intent(out) :: Eig0(Mon%NDimX)

integer :: iunit
integer :: i,iblk,nblk,ival,ipos

call read_SBlock(SBlock,SBlockIV,nblk,xy0file)

call create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig0,Mon%CICoef,Mon%IndN,nblk,NBasis,Mon%NDimX)

! deallocate SBLOCK
 do iblk=1,nblk
    associate(A => SBlock(iblk))
      deallocate(A%matY,A%matX,A%vec)
      deallocate(A%pos)
    end associate
 enddo
 ! deallocate IV part
 associate(A => SblockIV)
   if(A%n>0) then
      deallocate(A%vec)
      deallocate(A%pos)
   endif
 end associate

end subroutine convert_XY0_to_Y01

subroutine unpack_XY0_full(VecX0,VecY0,Eig0,C,IndN,NDimX,NBasis,xy0file)
implicit none

integer,intent(in)           :: NDimX,NBasis
integer,intent(in)           :: IndN(2,NDimX)
double precision,intent(in)  :: C(NBasis)
character(*),intent(in)      :: xy0file
double precision,intent(out) :: Eig0(NDimX), &
                                VecX0(NDimX,NDimX),VecY0(NDimX,NDimX)

integer                      :: iunit
integer                      :: i,ii,ip,iq,iblk,nblk,ipos
double precision             :: valX,valY
type(EblockData),allocatable :: SBlock(:)
type(EblockData)             :: SBlockIV
double precision,parameter   :: fac = 1d0/sqrt(2d0)

! read SBLOCK
open(newunit=iunit,file=xy0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
allocate(SBlock(nblk))
!print*, 'ival,nblkA',ival,nblkA
do iblk=1,nblk
   associate(A => SBlock(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n),A%matY(A%n,A%n),A%vec(A%n))
     read(iunit) A%pos,A%matX,A%matY,A%vec
   end associate
enddo
associate(A => SBlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  read(iunit) A%pos,A%vec
end associate
close(iunit)

VecX0 = 0
VecY0 = 0
Eig0  = 0

do iblk=1,nblk
   associate(B => Sblock(iblk))

     do i=1,B%n
        ipos = B%pos(i)
        VecY0(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
        VecX0(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
     enddo
     Eig0(B%l1:B%l2) = B%vec(1:B%n)

   end associate
enddo

!unpack (2, IV part)
associate(B => SblockIV)
  if(B%n>0) then
     do i=1,B%n
        ii = B%l1+i-1
        ipos = B%pos(i)
        ip = IndN(1,ipos)
        iq = IndN(2,ipos)

        valX = 1d0/(C(ip)+C(iq))
        valY = 1d0/(C(ip)-C(iq))

        VecX0(ipos,ii) = 0.5d0*fac*(valX-valY)
        VecY0(ipos,ii) = 0.5d0*fac*(valX+valY)
        Eig0(ii) = B%vec(i)
     enddo
  endif

end associate

! deallocate blocks
do iblk=1,nblk
   associate(B => Sblock(iblk))
     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)
   end associate
enddo
! deallocate IV blocks
associate(B => SblockIV)
  deallocate(B%vec)
  deallocate(B%pos)
end associate

end subroutine unpack_XY0_full

subroutine convert_XY_to_Z(EVecZ,CICoef,IndN,NDimX,NBasis,xyfile)
implicit none

integer,intent(in)           :: NDimX,NBasis
integer,intent(in)           :: IndN(2,NDimX)
character(*),intent(in)      :: xyfile
double precision,intent(in)  :: CICoef(NBasis)
double precision,intent(out) :: EVecZ(NDimX,NDimX)

integer                      :: i,j,ip,iq
double precision             :: fact
double precision,allocatable :: EVecX(:,:),EVecY(:,:)

allocate(EVecX(NDimX,NDimX),EVecY(NDimX,NDimX))

call readEVecXY(EVecX,EVecY,NDimX,xyfile)

do j=1,NDimX
   do i=1,NDimX
      ip = IndN(1,i)
      iq = IndN(2,i)
      fact = CICoef(ip) - CICoef(iq)
      EVecZ(i,j) = fact*(EVecY(i,j)-EVecX(i,j))
   enddo
enddo

deallocate(EVecY,EVecX)

end subroutine convert_XY_to_Z

subroutine writeampl(sij,ampfile)
implicit none

character(*) :: ampfile
double precision :: sij(:,:)
integer :: iunit

 open(newunit=iunit,file=ampfile,form='unformatted')
 write(iunit) sij
 close(iunit)

end subroutine writeampl

subroutine solve_cphf(M,WPot,e2indxy,Flags,NBas)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: M
type(DIISData) :: DIISBlock
type(Y01BlockData),allocatable :: Y01BlockM(:)

integer,intent(in) :: NBas
double precision,intent(in)  :: WPot(NBas,NBas)
double precision,intent(out) :: e2indxy

double precision,allocatable :: WxYY(:,:)
double precision,allocatable :: wVecxYY(:)
double precision,allocatable :: amps(:),vecR(:),delta(:)
double precision,allocatable :: OmM0(:)

integer      :: iter,i,pq,ip,iq
logical      :: conv=.FALSE.
character(8) :: nameunc
double precision :: error
integer,parameter          :: MaxIt = 50
double precision,parameter :: ThrDIIS = 1.d-8

 if(M%Monomer==1) nameunc='XY0_A'
 if(M%Monomer==2) nameunc='XY0_B'

 call init_DIIS(DIISBlock,M%NDimX,M%NDimX,Flags%DIISN)

 allocate(WxYY(NBas,NBas))
 allocate(wVecxYY(M%NDimX))
 call tran2MO(WPot,M%CMO,M%CMO,WxYY,NBas)
 wVecxYY = 0
 do pq=1,M%NDimX
    ip = M%IndN(1,pq)
    iq = M%IndN(2,pq)
    wVecxYY(pq) = WxYY(ip,iq)
 enddo

 ! read EigVals
 allocate(OmM0(M%NDimX))
 allocate(Y01BlockM(M%NDimX))
 call convert_XY0_to_Y01(M,Y01BlockM,OmM0,NBas,trim(nameunc))

 allocate(vecR(M%NDimX),amps(M%NDimX),delta(M%NDimX))

 ! zero iter
 call amplitudes_T1(OmM0,wVecxYY,amps,M%NDimX)

 write(LOUT,'(1x,a,5x,a)') 'ITER', 'ERROR'
 iter = 0
 do

 vecR = wVecxYY
 call dgemv('N',M%NDimX,M%NDimX,1.d0,M%PP,M%NDimX,amps,1,1d0,vecR,1)

 error = norm2(vecR)
 conv=error.le.ThrDIIS

 write(LOUT,'(1x,i3,f11.6)') iter,error

 if(conv) then
    exit
 elseif(.not.conv.and.iter.le.MaxIt) then
    iter = iter + 1
 elseif(.not.conv.and.iter.gt.MaxIt) then
   write(*,*) 'Error!!! E2ind DIIS not converged!'
   exit
 endif

 call amplitudes_T1(OmM0,vecR,delta,M%NDimX)
 amps = amps + delta
 if(iter>Flags%DIISOn) call use_DIIS(DIISBlock,amps,vecR)

 enddo

 !E2ind(X--Y)
 e2indxy = 0
 do pq=1,M%NDimX
    ip = M%IndN(1,pq)
    iq = M%IndN(2,pq)
    e2indxy = e2indxy + amps(pq)*WxYY(ip,iq)
 enddo
 e2indxy = 2d0*e2indxy

 call free_DIIS(DIISBlock)

 ! deallocate Y01Block
 do i=1,M%NDimX
    associate(Y => Y01BlockM(i))
      deallocate(Y%vec0)
    end associate
 enddo

 deallocate(delta,amps,vecR)
 deallocate(Y01BlockM)
 deallocate(OmM0)
 deallocate(wVecxYY,WxYY)

end subroutine solve_cphf

subroutine amplitudes_T1(deps,ints,res,NDimX)
implicit none

integer,intent(in) :: NDimX
double precision,intent(in)  :: deps(NDimX),ints(NDimX)
double precision,intent(out) :: res(NDimX)

integer :: pq,ip,iq

res = 0
do pq=1,NDimX
   res(pq) = res(pq) - ints(pq) / deps(pq)
enddo

end subroutine amplitudes_T1

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

subroutine readEval(EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVal

 close(iunit)

end subroutine readEval

subroutine readEvalZ(EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit)
 read(iunit) EVal

 close(iunit)

end subroutine readEvalZ

subroutine readEvalXY(EVal,NDim,fname)
implicit none

integer,intent(in)           :: NDim
character(*),intent(in)      :: fname
double precision,intent(out) :: EVal(NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit)
 read(iunit)
 read(iunit) EVal

 close(iunit)

end subroutine readEvalXY

subroutine readEvecZ(EVecZ,NDim,fname)
implicit none

integer          :: NDim
character(*)     :: fname
double precision :: EVecZ(NDim,NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVecZ

 close(iunit)

end subroutine readEvecZ

subroutine readEvecXY(EVecX,EVecY,NDim,fname)
implicit none

integer,intent(in)           :: NDim
character(*),intent(in)      :: fname
double precision,intent(out) :: EVecX(NDim,NDim),EVecY(NDim,NDim)

integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVecX
 read(iunit) EVecY

 close(iunit)

end subroutine readEvecXY

subroutine unpack_Eig(SBlock,SBlockIV,nblk,Eig,NDimX)
implicit none

type(EBlockData)             :: SBlock(nblk),SBlockIV
integer,intent(in)           :: nblk,NDimX
double precision,intent(out) :: Eig(NDimX)

integer :: i,ii,iblk

do iblk=1,nblk
   associate(B => SBlock(iblk))
     Eig(B%l1:B%l2) = B%vec(1:B%n)
   end associate
enddo
associate(B => SBlockIV)
  do i=1,B%n
     ii = B%l1+i-1
     Eig(ii) = B%vec(i)
  enddo
end associate

end subroutine unpack_Eig

subroutine fill_Fmat(Fmat,Occ,NBas,variant)
implicit none

integer,intent(in)           :: NBas,variant
double precision,intent(in)  :: Occ(NBas)
double precision,intent(out) :: Fmat(NBas,NBas)

integer :: i,j

Fmat = 0
select case(variant)
case(0)
   ! HF
   do j=1,NBas
      do i=1,NBas
         Fmat(i,j) = Occ(i)*Occ(j)
      enddo
   enddo
case(1)
   ! BB functional
   do j=1,NBas
      do i=1,NBas
         Fmat(i,j) = sqrt(Occ(i)*Occ(j))
      enddo
   enddo
case(3)
   ! Power functional
   write(LOUT,*) 'POWER FUNCITONAL NOT READY YET!'
   stop

end select

end subroutine fill_Fmat

subroutine calc_resp(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec(Mon%NDimX,Mon%NDimX), &
                               EVal(Mon%NDimX)
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha(Mon%NDimX,Mon%NDimX)
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8

 Alpha = 0
 do t=1,Mon%NDimX

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if
    if(EVal(t).gt.SmallE.and.EVal(t).lt.BigE) then

    frac = 8d0*EVal(t)/(EVal(t)**2d0+freq)

    do pq=1,Mon%NDimX
       ip = Mon%IndN(1,pq)
       iq = Mon%IndN(2,pq)
       do rs=1,pq
          ir = Mon%IndN(1,rs)
          is = Mon%IndN(2,rs)

          Alpha(pq,rs) = Alpha(pq,rs) + &
                       (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                       (Mon%CICoef(is)+Mon%CICoef(ir)) * &
                       EVec(pq,t)*EVec(rs,t)*frac
          Alpha(rs,pq) = Alpha(pq,rs)

       enddo
    enddo
    endif
 enddo

end subroutine calc_resp

subroutine calc_resp_apsg(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
!double precision,intent(in) :: EVec(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN)), &
!                               EVal(2*(Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: EVec((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN)), &
                               EVal((Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
!double precision,intent(out) :: Alpha(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN))
double precision,intent(out) :: Alpha((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN))
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
integer :: NDimEx
double precision,parameter :: SmallE = 1.d-6

 NDimEx = (Mon%NDimX + Mon%NDimN)

 Alpha = 0
 ! p>q, r>s
 do t=1,NDimEx

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if

    if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
       frac = 8d0*EVal(t)/(EVal(t)**2+freq)

       do pq=1,Mon%NDimX
          ip = Mon%IndN(1,pq)
          iq = Mon%IndN(2,pq)
          do rs=1,pq
             ir = Mon%IndN(1,rs)
             is = Mon%IndN(2,rs)

             Alpha(pq,rs) = Alpha(pq,rs) + &
                          (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                          (Mon%CICoef(is)+Mon%CICoef(ir)) * &
                          EVec(pq,t)*EVec(rs,t)*frac
             Alpha(rs,pq) = Alpha(pq,rs)

          enddo
       enddo
    endif
 enddo

! p>q,r=s
!
  do t=1,NDimEx
 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if

     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)

        do pq=1,Mon%NDimX
           ip = Mon%IndN(1,pq)
           iq = Mon%IndN(2,pq)
           do ir=1,Mon%NDimN

!              Alpha(pq,ir) = Alpha(pq,ir) + &
!                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
!                          Mon%CICoef(ir) * &
!                           EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,pq) = Alpha(pq,ir)

              Alpha(pq,Mon%NDimX+ir) = Alpha(pq,Mon%NDimX+ir) + &
                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                           Mon%CICoef(ir) * &
                           !EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
                           EVec(pq,t)*EVec(Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,pq) = Alpha(pq,Mon%NDimX+ir)

           enddo
        enddo

     endif
  enddo

 print*, 'aaa?',2*Mon%NDimX+Mon%NDimN,NDimEx

 ! p=q,r=s
  do t=1,NDimEx

 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
 !      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if

     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)

        do ip=1,Mon%NDimN
           do ir=1,ip !Mon%NDimN

!              Alpha(ip,ir) = Alpha(ip,ir) + &
!                           Mon%CICoef(ip) * &
!                           Mon%CICoef(ir) * &
!                           EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,ip) = Alpha(ip,ir)

              Alpha(Mon%NDimX+ip,Mon%NDimX+ir) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir) + &
                           Mon%CICoef(ip) * &
                           Mon%CICoef(ir) * &
                           !EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
                           EVec(Mon%NDimX+ip,t)*EVec(Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,Mon%NDimX+ip) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir)


           enddo
        enddo
     endif
  enddo

print*, 'resp,Alpha',norm2(Alpha)

end subroutine calc_resp_apsg

subroutine calc_resp_apsg2(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec((Mon%NDimX+Mon%NDimN)**2), &
                               EVal((Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN))
double precision :: frac,fact
integer,allocatable :: MIndEx(:,:)
double precision,allocatable :: VecEx(:)
integer :: pq,rs,t,ip,iq,ir,is
integer :: i
integer :: NDimEx
integer :: coef
double precision,parameter :: SmallE = 1.d-6, BigE = 1.d20

 NDimEx = (Mon%NDimX + Mon%NDimN)
 coef = 1
 ! with PINOVEC
 !coef = 2

 allocate(VecEx(NDimEx**2),MIndEx(2,NDimEx))

 MIndEx = Mon%IndNx

 VecEx = 0
 ! prepare extended vec
 do pq=1,NDimEx
    if(pq<=Mon%NDimX) then

        ip = Mon%IndN(1,pq)
        iq = Mon%IndN(2,pq)

        fact = Mon%CICoef(iq) + Mon%CICoef(ip)
        do i=1,coef*NDimEx
           VecEx((i-1)*coef*NDimEx+pq) = fact * &
                                          EVec((i-1)*coef*NDimEx+pq)
        enddo

    elseif(pq>Mon%NDimX) then

        ir = pq - Mon%NDimX
        fact = Mon%CICoef(ir)
        do i=1,coef*NDimEx
           VecEx((i-1)*coef*NDimEx+pq) = fact * &
                                         EVec((i-1)*coef*NDimEx+pq)
        enddo

    endif
 enddo

 Alpha = 0d0
 do t=1,NDimEx

    if(EVal(t).gt.SmallE.and.EVal(t).lt.BigE) then

       frac = 8d0*EVal(t)/(EVal(t)**2d0+freq)

       do pq=1,NDimEx
          ip = MIndEx(1,pq)
          iq = MIndEx(2,pq)
          do rs=1,NDimEx
          !do rs=1,pq
             ir = MIndEx(1,rs)
             is = MIndEx(2,rs)

             Alpha(pq,rs) = Alpha(pq,rs) + &
                          !EVec(pq,t)*EVec(rs,t)*frac
                          !VecEx((pq-1)*coef*NDimEx+t)*VecEx((rs-1)*coef*NDimEx+t)*frac
                          frac*VecEx((t-1)*coef*NDimEx+pq)*VecEx((t-1)*coef*NDimEx+rs)
             !Alpha(rs,pq) = Alpha(pq,rs)

          enddo
       enddo
    endif
 enddo

print*, 'resp2,Alpha',norm2(Alpha)

deallocate(MIndEx)
deallocate(VecEx)

end subroutine calc_resp_apsg2

subroutine calculateInitialA(A2,A0Blk,A0BlkIV,nblk,NDimX,abpm0file)
! a) read ABPM0 blocks from file
! b) calculate A2 = ABPM - ABPM0
implicit none

integer,intent(in)  :: NDimX
integer,intent(out) :: nblk
character(*)        :: abpm0file
double precision    :: A2(NDimX,NDimX)
type(EblockData)    :: A0BlkIV
type(EblockData),allocatable :: A0Blk(:)

integer :: i,j,ipos,jpos,iblk,ii
! test
double precision,allocatable :: A0(:,:)

call read_ABPM0Block(A0Blk,A0BlkIV,nblk,abpm0file)

!allocate(A0(NDimX,NDimX))
!A0 = 0
do iblk=1,nblk
   associate(B => A0Blk(iblk))

     ! A2 = ABPM - A0PM
     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           A2(ipos,jpos) = A2(ipos,jpos) - B%matX(i,j)
           !A0(ipos,jpos) = B%matX(i,j)
        enddo
     enddo

   end associate
enddo
!IV block
associate(B => A0BlkIV)
  do i=1,B%n
     ii = B%pos(i)
     A2(ii,ii) = A2(ii,ii) - B%vec(i)
     !A0(ii,ii) = B%vec(i)
  enddo
end associate

!print*, 'test: A0',norm2(A0)
!deallocate(A0)

end subroutine calculateInitialA

subroutine calculateLambda(Lambda,LambdaIV,omega,NDimX,nblk,A0Blk,A0BlkIV)
! calculate Lambda by inversion of blocks ABPM0 blocks
! Lambda = (A0+omega^2)-1
implicit none

integer,intent(in)  :: NDimX,nblk
double precision,intent(in) :: omega
type(EblockData)    :: A0Blk(nblk),A0BlkIV
type(EblockData),allocatable :: Lambda(:)
type(EblockData)             :: LambdaIV

integer :: i,ii,j,ipos,jpos,iblk
integer :: info
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:)
! for testing full Lambda
double precision,allocatable :: A0Inv(:,:)

!allocate(A0Inv(NDimX,NDimX))
!A0Inv=0

Lambda = A0Blk
do iblk=1,nblk
   associate(B => Lambda(iblk), A => A0Blk(iblk))

     allocate(ipiv(B%n),work(B%n))

     do i=1,B%n
        B%matX(i,i) = B%matX(i,i) + omega
     enddo

     call dgetrf(B%n,B%n,B%matX,B%n,ipiv,info)
     call dgetri(B%n,B%matX,B%n,ipiv,work,B%n,info)

     !! test Lambda
     !do j=1,B%n
     !   jpos=B%pos(j)
     !   do i=1,B%n
     !      ipos=B%pos(i)
     !      A0Inv(ipos,jpos) = B%matX(i,j)
     !   enddo
     !enddo

     deallocate(work,ipiv)

   end associate
enddo

! IV block
associate(B => A0blkIV, A => LambdaIV)
  A%n = B%n
  A%l1 = B%l1
  A%l2 = B%l2
  allocate(A%pos(A%n),A%vec(A%n))
  A%pos = B%pos
  do i=1,B%n
     ii = B%pos(i)
     A%vec(i) = 1d0 / (B%vec(i) + omega)
     !A0Inv(ii,ii) = 1d0 / (B%vec(i) + omega)
  enddo
end associate

!print*, 'test2-2:',norm2(A0Inv)
!deallocate(A0Inv)

end subroutine calculateLambda

subroutine Cmat_iterDIIS(CTilde,NDimX,NCholesky,nblk,Lambda,LambdaIV,ABPTilde,A2,iStats)
use diis
use class_IterStats
implicit none

type(SaptDIIS) :: SAPT_DIIS
type(DIISData) :: DIISBlock

type(EblockData)                :: Lambda(nblk),LambdaIV
class(IterStats), intent(inout) :: iStats
integer,intent(in) :: nblk,NDimX,NCholesky
double precision, intent(in)  :: A2(NDimX,NDimX), ABPTilde(NDimX,NCholesky)
double precision, intent(out) :: CTilde(NDimX,NCholesky)
! for testing full Lambda (passed a work)
!double precision, intent(in)  :: work(NDimX,NDimX)

integer :: N
double precision :: A2CTilde(NDimX,NCholesky), CTilde_prev(NDimX,NCholesky)
double precision :: norm
double precision :: Thresh

!print*, 'SAPTDiis: test'
!print*, 'DIISN',SAPT_DIIS%DIISN
!print*, 'DIISON',SAPT_DIIS%DIISOn
!print*, 'maxIter',SAPT_DIIS%maxIter
!print*, 'Thresh',SAPT_DIIS%Threshold

call init_DIIS(DIISBlock,NDimX*NCholesky,NDimX*NCholesky,SAPT_DIIS%DIISN)

CTilde_prev = 0.0d0
N = 1

do
   call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,CTilde,NDimX,0.0d0,A2CTilde,NDimX)
   A2CTilde = ABPTilde - A2CTilde

   ! with full Lambda(=work)
   !call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,work,NDimX,A2CTilde,NDimX,0.0d0,CTilde,NDimX)
   call ABPM_HALFTRAN_LR(A2CTilde,CTilde,Lambda,LambdaIV,nblk,NDimX,NCholesky,0)

   if(N > SAPT_DIIS%DIISOn ) then
       call use_DIIS(DIISBlock, CTilde, CTilde - CTilde_prev)
   endif

   norm = norm2(CTilde - CTilde_prev)
   !print '(a, i3, a, e)', "norm (", N, ") = ", norm

   if ((norm < SAPT_DIIS%Thresh) .or. (N .ge. SAPT_DIIS%maxIter .and. SAPT_DIIS%maxIter .ge. 0)) then
       call iStats%addN(N)
       exit
   endif

   CTilde_prev = CTilde
   N = N + 1

enddo

call free_DIIS(DIISBlock)

end subroutine Cmat_iterDIIS

subroutine print_en(string,val,nline)
implicit none

character(*),intent(in)     :: string
double precision,intent(in) :: val
logical,intent(in)          :: nline

if(nline) then
   write(LOUT,'(/1x,a,t19,a,f16.8)') string, "=", val
else
   write(LOUT,'(1x,a,t19,a,f16.8)')  string, "=", val
endif

end subroutine print_en

end module sapt_utils
