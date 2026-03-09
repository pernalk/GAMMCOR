module sapt_ind
use types
use tran
use sapt_utils
use read_external

implicit none

contains

subroutine e2ind_icerpa(Flags,A,B,SAPT)
use timing
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas
double precision :: e2ba,e2ab,e2ic
double precision :: Tcpu,Twall

 call clock('START',Tcpu,Twall)

 NBas = A%NBasis

 call solve_cphf(A,B%WPot,e2ba,Flags,NBas)
 call solve_cphf(B,A%WPot,e2ab,Flags,NBas)

 e2ic = (e2ab + e2ba)

 write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0
 write(LOUT,'(1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000d0
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0
 SAPT%e2ind = e2ic

 call clock('E2ind',Tcpu,Twall)

end subroutine e2ind_icerpa

subroutine e2ind_resp(Flags,A,B,SAPT)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
! test
integer             :: info
integer,allocatable :: ipiv(:)
double precision,allocatable :: wtest(:),car(:),Work(:)
!
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: e2ba,e2ab,e2iu,e2ic
double precision :: tmp

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! read EigValA_B
 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX), &
          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(AlphaA(A%NDimX,A%NDimX),AlphaB(B%NDimX,B%NDimX), &
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 call calc_resp(EVecA,OmA,AlphaA,0d0,A)
 call calc_resp(EVecB,OmB,AlphaB,0d0,B)

 e2ba=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    do rs=1,A%NDimX
       ir = A%IndN(1,rs)
       is = A%IndN(2,rs)

       e2ba = e2ba + &
            WbAA(ip,iq)*AlphaA(pq,rs)*WbAA(ir,is)

    enddo
 enddo
 e2ba = -0.5d0*e2ba
 !write(LOUT,'(/,1x,a,f16.8)') 'Ind(B-->A)   = ', e2ba*1000d0
 call print_en('Ind(B-->A)',e2ba*1000,.true.)

 e2ab=0
 do pq=1,B%NDimX
    ip = B%IndN(1,pq)
    iq = B%IndN(2,pq)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       e2ab = e2ab + &
            WaBB(ip,iq)*AlphaB(pq,rs)*WaBB(ir,is)

    enddo
 enddo
 e2ab = -0.5d0*e2ab
 !write(LOUT,'(1x,a,f16.8)') 'Ind(A-->B)   = ', e2ab*1000d0
 call print_en('Ind(A-->B)',e2ab*1000,.false.)

 e2ic = (e2ab + e2ba)
 !write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0
 call print_en('E2ind',e2ic*1000,.false.)

 SAPT%e2ind = e2ic

 deallocate(WaBB,WbAA,AlphaB,AlphaA)
 deallocate(OmB,EVecB,Oma,EVecA)

end subroutine e2ind_resp

subroutine e2ind_cpld(Flags,A,B,SAPT)
! calculate only coupled E2ind
! used also for  cubic E2ind
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: tmpA(:),tmpB(:)

integer          :: NAO,NBas,info
integer          :: i,j,pq,ip,iq,rs,ir,is
double precision :: fact
double precision :: e2ba,e2ab,e2iu,e2ic
double precision :: e2ab_unc,e2ba_unc,e2ic_unc
character(:),allocatable   :: propA,propB
double precision,parameter :: BigE   = 1.D10
double precision,parameter :: SmallE = 1.D-6

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif
 NAO = SAPT%NAO

 ! E2IND(B<--A)

 allocate(EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(tmpB(B%NDimX),WaBB(NBas,NBas))

 if(B%Cubic) then
   if(B%ACAlpha==B%ACAlpha0) then
      propB = 'PROP_B0'
   elseif(B%ACAlpha==B%ACAlpha1) then
      propB = 'PROP_B1'
   elseif(B%ACAlpha==B%ACAlpha2) then
      propB = 'PROP_B2'
   endif
   call convert_XY_to_Z(EVecB,B%CICoef,B%IndN,B%NDimX,NBas,propB)
   call readEvalXY(OmB,B%NDimX,propB)
 else
   !call readresp(EVecB,OmB,B%NDimX,'PROP_B')
   call readEvecZ(EVecB,B%NDimX,'PROP_B')
   call readEvalZ(OmB,B%NDimX,'PROP_B')
 endif

 call tran_AO2MO2(A%WPot,B%CMO,B%CMO,WaBB,NAO,NBas)

 tmpB = 0d0
 do rs=1,B%NDimX
    ir = B%IndN(1,rs)
    is = B%IndN(2,rs)

    fact = (B%CICoef(is)+B%CICoef(ir))*WaBB(ir,is)

    do j=1,B%NDimX
       tmpB(j) = tmpB(j) + fact*EVecB(rs,j)
    enddo

 enddo

 e2ab = 0d0
 do j=1,B%NDimX
    if(OmB(j).lt.BigE.and.OmB(j).gt.SmallE) then
       e2ab = e2ab + tmpB(j)**2 / OmB(j)
    endif
 enddo

 e2ab = -4.0d0*e2ab
 write(LOUT,'(/1x,a,f16.8)') 'Ind(B<--A)     = ', e2ab*1000d0

 deallocate(tmpB)
 deallocate(EVecB,OmB)

 ! E2IND(A<--B)

 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX))
 allocate(tmpA(A%NDimX),WbAA(NBas,NBas))

 if(A%Cubic) then
   if(A%ACAlpha==A%ACAlpha0) then
      propA = 'PROP_A0'
   elseif(A%ACAlpha==A%ACAlpha1) then
      propA = 'PROP_A1'
   elseif(A%ACAlpha==A%ACAlpha2) then
      propA = 'PROP_A2'
   endif
   call convert_XY_to_Z(EVecA,A%CICoef,A%IndN,A%NDimX,NBas,propA)
   call readEvalXY(OmA,A%NDimX,propA)
 else
   !call readresp(EVecA,OmA,A%NDimX,'PROP_A')
   call readEvecZ(EVecA,A%NDimX,'PROP_A')
   call readEvalZ(OmA,A%NDimX,'PROP_A')
 endif

 call tran_AO2MO2(B%WPot,A%CMO,A%CMO,WbAA,NAO,NBas)

 tmpA = 0d0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)

    fact = (A%CICoef(iq)+A%CICoef(ip))*WbAA(ip,iq)

    do i=1,A%NDimX
       tmpA(i) = tmpA(i) + fact*EVecA(pq,i)
    enddo

 enddo

 e2ba = 0d0
 do i=1,A%NDimX
    if(OmA(i).lt.BigE.and.OmA(i).gt.SmallE) then
       e2ba = e2ba + tmpA(i)**2 / OmA(i)
    endif
 enddo

 e2ba = -4.0d0*e2ba
 write(LOUT,'(1x,a,f16.8)') 'Ind(A<--B)     = ', e2ba*1000d0

 e2ic = (e2ab + e2ba)

 if(A%Cubic.or.B%Cubic) then

   if(A%ACAlpha==A%ACAlpha0.or.B%ACAlpha==B%ACAlpha0) SAPT%e2ind_a0 = e2ic
   if(A%ACAlpha==A%ACAlpha1.or.B%ACAlpha==B%ACAlpha1) SAPT%e2ind_a1 = e2ic
   if(A%ACAlpha==A%ACAlpha2.or.B%ACAlpha==B%ACAlpha2) SAPT%e2ind_a2 = e2ic

   !write(LOUT,'(1x,a,f16.8)') 'A: Alpha     = ',A%ACAlpha
   !write(LOUT,'(1x,a,f16.8)') 'B: Alpha     = ',B%ACAlpha
   write(LOUT,'(1x,a,f16.8)') 'E2ind(Alpha)   = ',e2ic*1000d0

 else

    SAPT%e2ind = e2ic
    write(LOUT,'(1x,a,f16.8)') 'E2ind          = ', e2ic*1000d0

 endif

 deallocate(tmpA)
 deallocate(WaBB,WbAA)
 deallocate(EVecA,OmA)

end subroutine e2ind_cpld

subroutine e2ind_unc(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

double precision,allocatable :: OmA0(:),OmB0(:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: tmpA(:),tmpB(:)
integer :: NBas,info
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: fact
double precision :: e2ab_unc,e2ba_unc,e2ic_unc
double precision,parameter :: BigE   = 1.D10
double precision,parameter :: SmallE = 1.D-6

 ! uncoupled for CAS only (for now)
 if(Flags%ICASSCF/=1) then
   write(LOUT,'(1x,a)') 'E2ind(unc) implented only for CAS!'
   return
 endif
 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

 ! E2IND(B<--A)
 allocate(Y01BlockB(B%NDimX))
 allocate(tmpB(B%NDimX),OmB0(B%NDimX),WaBB(NBas,NBas))

 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
 call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')

 tmpB = 0d0
 do rs=1,B%NDimX
    ir = B%IndN(1,rs)
    is = B%IndN(2,rs)

    fact = (B%CICoef(is)+B%CICoef(ir))*WaBB(ir,is)

    associate(Y => Y01BlockB(rs))
       tmpB(Y%l1:Y%l2) = tmpB(Y%l1:Y%l2) + fact * Y%vec0(1:Y%n)
    end associate

 enddo

 e2ab_unc = 0d0
 do j=1,B%NDimX
    if(OmB0(j).lt.BigE.and.OmB0(j).gt.SmallE) then
       e2ab_unc = e2ab_unc + tmpB(j)**2 / OmB0(j)
    endif
 enddo

 e2ab_unc = -4.0d0*e2ab_unc
 write(LOUT,'(/1x,a,f16.8)')  'Ind(B<--A,unc) = ', e2ab_unc*1000d0

 do i=1,B%NDimX
    associate(Y => Y01BlockB(i))
      deallocate(Y%vec0)
    end associate
 enddo
 deallocate(Y01BlockB)
 deallocate(OmB0,tmpB)

 ! E2IND(A<--B)

 allocate(Y01BlockA(A%NDimX))
 allocate(tmpA(A%NDimX),OmA0(A%NDimX),WbAA(NBas,NBas))

 call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 tmpA = 0d0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)

    fact = (A%CICoef(iq)+A%CICoef(ip))*WbAA(ip,iq)

    associate(Y => Y01BlockA(pq))
       tmpA(Y%l1:Y%l2) = tmpA(Y%l1:Y%l2) + fact * Y%vec0(1:Y%n)
    end associate

 enddo

 e2ba_unc = 0d0
 do i=1,A%NDimX
    if(OmA0(i).lt.BigE.and.OmA0(i).gt.SmallE) then
       e2ba_unc = e2ba_unc + tmpA(i)**2 / OmA0(i)
    endif
 enddo

 e2ba_unc = -4.0d0*e2ba_unc
 e2ic_unc = (e2ab_unc + e2ba_unc)

 write(LOUT,'(1x,a,f16.8)') 'Ind(A<--B,unc) = ', e2ba_unc*1000d0

 do i=1,A%NDimX
    associate(Y => Y01BlockA(i))
      deallocate(Y%vec0)
    end associate
 enddo
 deallocate(Y01BlockA)
 deallocate(OmA0,tmpA)

 write(LOUT,'(1x,a,f16.8)') 'E2ind(unc)     = ', e2ic_unc*1000d0
 SAPT%e2ind_unc = e2ic_unc

end subroutine e2ind_unc

subroutine e2ind(Flags,A,B,SAPT)

use timing
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

double precision,allocatable :: OmA(:),OmB(:),&
                                OmA0(:),OmB0(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: tmpA(:),tmpB(:)
double precision,allocatable :: tmp01(:)
integer :: NBas,info
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: fact
double precision :: e2ba,e2ab,e2iu,e2ic
double precision :: e2ab_unc,e2ba_unc,e2ic_unc
double precision :: Tcpu,Twall
double precision,parameter :: BigE   = 1.D10
double precision,parameter :: SmallE = 1.D-6

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

 call clock('START',Tcpu,Twall)

 ! E2IND(B<--A)

 allocate(EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(tmpB(B%NDimX),WaBB(NBas,NBas))

 call readresp(EVecB,OmB,B%NDimX,'PROP_B')
 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)

 ! uncoupled - for CAS only
 if(Flags%ICASSCF==1) then
    allocate(Y01BlockB(B%NDimX))
    allocate(tmp01(B%NDimX),OmB0(B%NDimX))
    call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')
 endif

 tmpB = 0d0
 do rs=1,B%NDimX
    ir = B%IndN(1,rs)
    is = B%IndN(2,rs)

    fact = (B%CICoef(is)+B%CICoef(ir))*WaBB(ir,is)

    do j=1,B%NDimX
       tmpB(j) = tmpB(j) + fact*EVecB(rs,j)
    enddo

 enddo

 e2ab = 0d0
 do j=1,B%NDimX
    if(OmB(j).lt.BigE.and.OmB(j).gt.SmallE) then
       e2ab = e2ab + tmpB(j)**2 / OmB(j)
    endif
 enddo

 e2ab = -4.0d0*e2ab
 !write(LOUT,'(/1x,a,f16.8)') 'Ind(B<--A)     = ', e2ab*1000d0
 call print_en('Ind(B<--A)',e2ab*1000,.true.)

 deallocate(tmpB)
 deallocate(EVecB,OmB)

 if(Flags%ICASSCF==1) then

    tmp01 = 0
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (B%CICoef(is)+B%CICoef(ir))*WaBB(ir,is)

       associate(Y => Y01BlockB(rs))
          tmp01(Y%l1:Y%l2) = tmp01(Y%l1:Y%l2) + fact * Y%vec0(1:Y%n)
       end associate

    enddo

    e2ab_unc = 0d0
    do j=1,B%NDimX
       if(OmB0(j).lt.BigE.and.OmB0(j).gt.SmallE) then
          e2ab_unc = e2ab_unc + tmp01(j)**2 / OmB0(j)
       endif
    enddo

    e2ab_unc = -4.0d0*e2ab_unc

    do i=1,B%NDimX
       associate(Y => Y01BlockB(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockB)
    deallocate(OmB0,tmp01)

 endif

 ! E2IND(A<--B)

 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX))
 allocate(tmpA(A%NDimX),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 ! uncoupled - for CAS only
 if(Flags%ICASSCF==1) then
    allocate(Y01BlockA(A%NDimX))
    allocate(tmp01(A%NDimX),OmA0(A%NDimX))
    call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
 endif

 tmpA = 0d0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)

    fact = (A%CICoef(iq)+A%CICoef(ip))*WbAA(ip,iq)

    do i=1,A%NDimX
       tmpA(i) = tmpA(i) + fact*EVecA(pq,i)
    enddo

 enddo

 e2ba = 0d0
 do i=1,A%NDimX
    if(OmA(i).lt.BigE.and.OmA(i).gt.SmallE) then
       e2ba = e2ba + tmpA(i)**2 / OmA(i)
    endif
 enddo

 e2ba = -4.0d0*e2ba
 !write(LOUT,'(1x,a,f16.8)') 'Ind(A<--B)     = ', e2ba*1000d0
 call print_en('Ind(A<--B)',e2ba*1000,.false.)

 e2ic = (e2ab + e2ba)
 SAPT%e2ind = e2ic
 !write(LOUT,'(1x,a,f16.8)') 'E2ind          = ', e2ic*1000d0
 call print_en('E2ind',e2ic*1000,.false.)

 if(Flags%ICASSCF==1) then

    tmp01 = 0d0
    do pq=1,A%NDimX
       ip = A%IndN(1,pq)
       iq = A%IndN(2,pq)

       fact = (A%CICoef(iq)+A%CICoef(ip))*WbAA(ip,iq)

       associate(Y => Y01BlockA(pq))
          tmp01(Y%l1:Y%l2) = tmp01(Y%l1:Y%l2) + fact * Y%vec0(1:Y%n)
       end associate

    enddo

    e2ba_unc = 0d0
    do i=1,A%NDimX
       if(OmA0(i).lt.BigE.and.OmA0(i).gt.SmallE) then
          e2ba_unc = e2ba_unc + tmp01(i)**2 / OmA0(i)
       endif
    enddo

    e2ba_unc = -4.0d0*e2ba_unc
    e2ic_unc = (e2ab_unc + e2ba_unc)
    SAPT%e2ind_unc = e2ic_unc

    !write(LOUT,'(/1x,a,f16.8)') 'Ind(B<--A,unc) = ', e2ab_unc*1000d0
    !write(LOUT,'(1x,a,f16.8)')  'Ind(A<--B,unc) = ', e2ba_unc*1000d0
    !write(LOUT,'(1x,a,f16.8)')  'E2ind(unc)     = ', e2ic_unc*1000d0

    call print_en('Ind(B<--A,unc)',e2ab_unc*1000d0,.true.)
    call print_en('Ind(A<--B,unc)',e2ba_unc*1000d0,.false.)
    call print_en('E2ind(unc)',e2ic_unc*1000d0,.false.)

    do i=1,A%NDimX
       associate(Y => Y01BlockA(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockA)
    deallocate(OmA0,tmp01)
 endif

 ! calculate deexcitations
 !if(SAPT%Wexcit) call e2ind_dexc(Flags,A,B,SAPT)

 deallocate(tmpA)
 deallocate(WaBB,WbAA)
 deallocate(EVecA,OmA)

 call clock('E2ind ',Tcpu,Twall)

end subroutine e2ind

subroutine e2ind_apsg(Flags,A,B,SAPT)

implicit none
type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas,ADimEx,BDimEx
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: termsBA(3), termsAB(3)
integer :: coef,coef2
double precision :: e2ba,e2ab
double precision :: e2iu,e2ic
double precision :: e2tmp, tmp

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

 coef  = 1
 coef2 = 1

 ! with PINOVEC
 ! coef  = 2
 ! coef2 = 4

 ADimEx = A%NDimX + A%NDimN
 BDimEx = B%NDimX + B%NDimN

! read EigValA_B
 allocate(EVecA(coef*ADimEx*coef*ADimEx),OmA(coef*ADimEx),&
          EVecB(coef*BDimEx*coef*BDimEx),OmB(coef*BDimEx))
 allocate(AlphaA(ADimEx,ADimEx),AlphaB(BDimEx,BDimEx), &
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,coef*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,coef*BDimEx,'PROP_B')

 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 call calc_resp_apsg2(EVecA,OmA,AlphaA,0d0,A)
 call calc_resp_apsg2(EVecB,OmB,AlphaB,0d0,B)

 ! test
 termsBA(1)=0
 do pq=1,ADimEx
    ip = A%IndNx(1,pq)
    iq = A%IndNx(2,pq)
    do rs=1,ADimEx
       ir = A%IndNx(1,rs)
       is = A%IndNx(2,rs)

       termsBA(1) = termsBA(1) + &
            WbAA(ip,iq)*AlphaA(pq,rs)*WbAA(ir,is)

    enddo
 enddo
 termsBA(1) = -0.5d0*termsBA(1)

 termsAB(1)=0
 do pq=1,BDimEx
    ip = B%IndNx(1,pq)
    iq = B%IndNx(2,pq)
    do rs=1,BDimEx
       ir = B%IndNx(1,rs)
       is = B%IndNx(2,rs)

       termsAB(1) = termsAB(1) + &
            WaBB(ip,iq)*AlphaB(pq,rs)*WaBB(ir,is)


    enddo
 enddo
 termsAB(1) = -0.5d0*termsAB(1)
! print*, 'testAB',termsAB(1)*1d3
 e2ic=0
 e2ic=(termsBA(1)+termsAB(1))
 write(LOUT,'(/1x,a,f16.8)') 'E2ind      = ', e2ic*1000d0

 SAPT%e2ind = e2ic

 deallocate(EVecA,OmA,EVecB,OmB)
 deallocate(AlphaA,AlphaB,WaBB,WbAA)

end subroutine e2ind_apsg

subroutine e2ind_dexc(Flags,A,B,SAPT)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: tmpA(:),tmpB(:)
double precision,allocatable :: Wij(:),work(:)

integer :: i,j,pq,ip,iq,rs,ir,is
integer :: NBas,info
integer :: nStates,iStA,nStSum,offset
integer :: NumOSym(8),NumStSym(16),IStSy(16), &
           NSym,NStSym
integer,allocatable :: NSymAO(:)

double precision :: fact
double precision :: e2ba,e2ab,e2iu,e2ic
double precision,parameter :: BigE   = 1.D10
double precision,parameter :: SmallE = 1.D-6

NBas = A%NBasis

write(lout,'(/1x,a)') 'W_ij corrections for E2ind in excited states: (A*<-B)'

! safety check
! in the future this subroutine should work
! also for the (A-B*) case
if(B%Wexcit) then
   write(lout,'(/1x,a)') '(B*-A) not avail. Set (A*-B) for W_ij corrections!'
   return
endif

! establish the type of correction: W_0j, W_1j, ...
! get info about monomer A
allocate(NSymAO(NBas))
call sym_inf_molpro('2RDMA',NumOSym,NSym,NumStSym,IStSy,NStSym,NSymAO,NBas)
deallocate(NSymAO)

! number of states summed over irreps
nStates = sum(NumStSym)

! reference state in A
iStA = A%InSt(1,1)

! number of state above the reference
! accessible from SA-CAS calculation
!nStSum = nStates - istA + 1
nStSum = 2

write(lout,'(1x,a,i2)') 'The number of available states is: ', nStates
write(lout,'(1x,a,i2,a,i1/)') 'The reference state is:',A%InSt(1,1),'.',A%InSt(2,1)

! E2IND(A<--B)

allocate(Wij(nStSum))
allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX))
allocate(tmpA(A%NDimX),WbAA(NBas,NBas))

call readresp(EVecA,OmA,A%NDimX,'PROP_A')
call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

! skip negative excitations
offset = 0
do i=1,A%NDimX

   if(OmA(i)<0d0) then
      offset = offset + 1
   elseif(OmA(i)>0d0) then
      exit
   endif

enddo
!write(lout,'(1x,a,i3)') 'offset = ', offset
!write(lout,'(1x,a,f6.4)')  'OmA(1): ', OmA(1+offset)
!write(lout,'(1x,a,f6.4)')  'OmA(2): ', OmA(2+offset)

tmpA = 0d0
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)

   fact = (A%CICoef(iq)+A%CICoef(ip))

   do i=1,A%NDimX
   !do i=1+offset,nStSum+offset
      tmpA(i) = tmpA(i) + fact*EVecA(pq,i)*WbAA(ip,iq)
   enddo

enddo

!print*, 'tmpA',norm2(tmpA)

Wij = 0d0
do i=1+offset,nStSum+offset
   if(OmA(i).lt.BigE.and.OmA(i).gt.SmallE) then
      !print*, 'i,OmA',i,OmA(i),tmpA(i)**2
      Wij(i-offset) = Wij(i-offset) + tmpA(i)**2 / OmA(i)
      !print*, 'Wij',Wij(i-offset)
   endif
enddo
Wij = 4.0d0*Wij

! print results
j = iStA
do i=1,nStSum
   write(lout,'(1x,a,2i1,a,f12.6)') 'Wind_',iStA-1,j,' =',Wij(i)*1000d0
   j = j + 1
enddo

deallocate(Wij)
deallocate(WbAA,tmpA)
deallocate(EVecA,OmA)

end subroutine e2ind_dexc

subroutine e2ind_pino(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: i,j,ij
integer :: ip,iq,ir,is,pq,rs
integer :: coef,coef2,ADimEx,BDimEx
integer,allocatable :: AuxT(:,:)
double precision,allocatable :: tmpA(:),tmpB(:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision :: fact,fpq,frs
double precision :: e2ab,e2ba,e2ic,tmp
! test
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
double precision,allocatable :: AVecEx(:),BVecEx(:)

double precision,parameter :: SmallE = 1.d-20
!double precision,parameter :: SmallE = 1.d-6
double precision,parameter :: BigE = 1.d10

 write(LOUT,'(/1x,a,e12.4)') 'SmallE in E2Ind PINO:',SmallE

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! set dimensions
 ADimEx = A%NDimX + A%NDimN
 BDimEx = B%NDimX + B%NDimN
 NInte1 = NBas*(NBas+1)/2
 dimOA = A%num0+A%num1
 dimFA = NBas
 dimOB = B%num0+B%num1
 dimFB = NBas
 nOFA = dimOA*dimFA
 nOFB = dimOB*dimFB
 coef  = 1
 coef2 = 1

 allocate(AuxT(2,NInte1))
 AuxT = 0
 ij = 0
 do j=1,NBas
    do i=1,j
       ij = ij + 1
       AuxT(1,ij) = j
       AuxT(2,ij) = i
    enddo
 enddo

 do i=1,coef*ADimEx
    if(A%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer A: Negative EVal",i4,f16.8)') i,A%PP(i)

    endif
 enddo
 do i=1,coef*BDimEx
    if(B%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,B%PP(i)
    endif
 enddo

 ! E2IND(B<--A)

 allocate(tmpB(coef*BDimEx),WaBB(NBas,NBas))
 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)

 tmpB=0
 do rs=1,BDimEx
    ir = AuxT(1,rs)
    is = AuxT(2,rs)

    frs = 1d0
    if(ir==is) frs=0.5d0

    fact = frs*(B%CICoef(is)+B%CICoef(ir))*WaBB(ir,is)

       do i=1,coef*BDimEx
          tmpB(i) = tmpB(i) + fact*B%AP(i,rs)
       enddo
 enddo

 e2ab = 0d0
 do i=1,BDimEx
    if(B%PP(i).lt.BigE.and.B%PP(i).gt.SmallE) then
       e2ab = e2ab + tmpB(i)**2 / B%PP(i)
    endif
 enddo

 e2ab = -4.0d0*e2ab
 write(LOUT,'(1x,a,f16.8)') 'Ind(B<--A)     = ', e2ab*1000d0

 deallocate(WaBB,tmpB)

 ! E2IND(A<--B)

 allocate(tmpA(coef*ADimEx),WbAA(NBas,NBas))
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

 tmpA=0
 do pq=1,ADimEx
    ip = AuxT(1,pq)
    iq = AuxT(2,pq)

    fpq = 1d0
    if(ip==iq) fpq=0.5d0

    fact = fpq*(A%CICoef(iq)+A%CICoef(ip))*WbAA(ip,iq)

       do i=1,coef*ADimEx

         ! if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then

          tmpA(i) = tmpA(i) + &
                    fact*A%AP(i,pq)

         ! endif
       enddo
 enddo

 e2ba = 0d0
 do i=1,ADimEx
    if(A%PP(i).lt.BigE.and.A%PP(i).gt.SmallE) then
       e2ba = e2ba + tmpA(i)**2 / A%PP(i)
    endif
 enddo

 e2ba = -4.0d0*e2ba
 write(LOUT,'(1x,a,f16.8)') 'Ind(A<--B)     = ', e2ba*1000d0

 deallocate(WbAA,tmpA)

 e2ic = (e2ab + e2ba)
 SAPT%e2ind = e2ic
 write(LOUT,'(1x,a,f16.8)') 'E2ind          = ', e2ic*1000d0

end subroutine e2ind_pino

subroutine e2ind_o(Flags,A,B,SAPT)
!
! calculate uncoupled and coupled e2ind
! c.f. Eq (17) in https://doi.org/10.1063/1.4758455
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: Nbasis
double precision :: e2iBAa,e2iBAb,e2iABa,e2iABb
double precision :: e2ab_unc,e2ba_unc,e2ind_unc
double precision :: e2ab,e2ba,e2ind
double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)

NBasis = A%NBasis

allocate(Waa(NBasis,NBasis),Wab(NBasis,NBasis))
allocate(Wba(NBasis,NBasis),Wbb(NBasis,NBasis))

call tran2MO(A%WPot,B%UMO(:,:,1),B%UMO(:,:,1),Waa,NBasis)
call tran2MO(A%WPot,B%UMO(:,:,2),B%UMO(:,:,2),Wab,NBasis)

call tran2MO(B%WPot,A%UMO(:,:,1),A%UMO(:,:,1),Wba,NBasis)
call tran2MO(B%WPot,A%UMO(:,:,2),A%UMO(:,:,2),Wbb,NBasis)

! uncoupled
e2iBAa = e2ind_unc_o(Wba,A%UOrbE(:,1),A%NOa,A%NVa,NBasis)
e2iBAb = e2ind_unc_o(Wbb,A%UOrbE(:,2),A%NOb,A%NVb,NBasis)
e2ba_unc = e2iBAa + e2iBAb

e2iABa = e2ind_unc_o(Waa,B%UOrbE(:,1),B%NOa,B%NVa,NBasis)
e2iABb = e2ind_unc_o(Wab,B%UOrbE(:,2),B%NOb,B%NVb,NBasis)
e2ab_unc = e2iABa + e2iABb

e2ind_unc = e2ba_unc + e2ab_unc

call print_en('Ind(B<--A,unc)',e2ab_unc*1000d0,.true.)
call print_en('Ind(A<--B,unc)',e2ba_unc*1000d0,.false.)
call print_en('E2ind(unc)',e2ind_unc*1000d0,.false.)

SAPT%e2ind_unc = e2ind_unc

deallocate(Wbb,Wba)
deallocate(Wab,Waa)

contains

function e2ind_unc_o(Wmat,OrbE,no,nv,n) result(res)
implicit none

integer :: no,nv,n
double precision :: Wmat(n,n), OrbE(n)
double precision :: delta_e, res

integer :: ip, iq

res = 0d0
do iq=1,no
   do ip=1,nv
      delta_e = OrbE(iq) - OrbE(no+ip)
      res = res + Wmat(iq,no+ip)*Wmat(no+ip,iq)/delta_e
   enddo
enddo

end function e2ind_unc_o

end subroutine e2ind_o

end module sapt_ind
