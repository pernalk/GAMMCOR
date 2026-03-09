module sapt_Chol_pol

use tran, only : ABPM_HALFTRAN_GEN_L
use sapt_utils
use timing
use sapt_Chol_basic
use sapt_Chol_Cmat

implicit none

contains

subroutine e2disp_Chol(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! in coupled and uncoupled approximations
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: NCholesky
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
logical,allocatable          :: condOmA(:),condOmB(:)
double precision,allocatable :: OmA(:), OmB(:), &
                                OmA0(:),OmB0(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),&
                                tmp01(:,:),tmp02(:,:)
double precision,allocatable :: work(:)
double precision :: e2d,fact,tmp
double precision :: e2du,dea,deb
double precision :: inv_omega
! for Be ERPA:
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! print thresholds
 if(SAPT%IPrint>1) then
    write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
    write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
    write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
 endif

! set dimensions
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA  = dimOA*dimVA
 nOVB  = dimOB*dimVB

 NCholesky = SAPT%NCholesky

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),  &
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX),  &
          OmA0(A%NDimX),OmB0(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 ! uncoupled - works for CAS only
 if(Flags%ICASSCF==1) then
    allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))

    call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
    call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')
 endif

allocate(work(B%NDimX))

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
        tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX))

! coupled
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative omega A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative omega B!',i,OmB(i)
enddo

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then

 tmp1=0
 tmp01=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

       associate(Y => Y01BlockA(pq))
          tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
       end associate

    enddo
 enddo
 ! coupled
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp2,A%NDimX)

 ! uncoupled
 tmp02=0
 do rs=1,B%NDimX
    associate(Y => Y01BlockB(rs))
      call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,tmp02(:,Y%l1:Y%l2),A%NDimX)
    end associate
 enddo

elseif(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

 tmp1 = 0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    call dgemv('T',NCholesky,B%NDimX,1d0,B%OV,NCholesky,A%OV(:,pq),1,0d0,work,1)

    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(rs)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + &
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j) + &
                    EVecB(rs+(j-1)*B%NDimX)*tmp1(i,rs)
       enddo
    enddo
 enddo

endif ! end GVB select

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then
   ! uncoupled
    e2du = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then


          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          e2du = e2du + tmp02(i,j)**2*inv_omega

          endif
       enddo
    enddo
    SAPT%e2disp_unc = -16d0*e2du

    e2du = -16d0*e2du*1000d0

    call writeampl(tmp02,'PROP_AB0')

endif

 allocate(condOmA(A%NDimX),condOmB(B%NDimX))
 condOmA = (abs(OmA).gt.SmallE.and.abs(OmA).lt.BigE)
 condOmB = (abs(OmB).gt.SmallE.and.abs(OmB).lt.BigE)

 e2d = 0d0
 do j=1,B%NDimX
    if(condOmB(j)) then
       do i=1,A%NDimX
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
!             .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

             if(condOmA(i)) then
                e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))
             endif
       enddo
    endif
 enddo
 SAPT%e2disp  = -16d0*e2d

 e2d  = -16d0*e2d*1000d0

 call print_en('E2disp',e2d,.true.)
 call print_en('E2disp(unc)',e2du,.false.)

 ! write amplitude to a file
 call writeampl(tmp2,'PROP_AB')

 !! calucate semicoupled and dexcitations
 !if(SAPT%SemiCoupled) call e2disp_semi(Flags,A,B,SAPT)

 !! calculate extrapolated E2disp
 !if(A%Cubic.or.B%Cubic) call e2disp_cpld(Flags,A,B,SAPT)

 !! calculate Wterms (deexcitations)
 !if(SAPT%Wexcit) call e2inddisp_dexc(Flags,A,B,SAPT)

 deallocate(work)

 if(Flags%ICASSCF==1) then
    ! deallocate Y01Block
    do i=1,A%NDimX
       associate(Y => Y01BlockA(i))
         deallocate(Y%vec0)
       end associate
    enddo
    do i=1,B%NDimX
       associate(Y => Y01BlockB(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockB,Y01BlockA)
 endif

 deallocate(condOmB,condOmA)
 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,OmA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp_Chol

subroutine e2ind_CAlphaTilde_block(Flags,A,B,SAPT)
!
! calculate 2nd order induction energy
! using expansion of C(0) in alpha around alpha=0, up to Max_Cn order
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NCholesky,NBas
integer :: i,j,ip,iq,ipq,ir,is,irs
integer :: nblkA,nblkB
double precision :: e2ba, e2ab
double precision,allocatable :: WaBB(:,:),WbAA(:,:),   &
                                WaChBB(:),WbChAA(:)
double precision,allocatable :: A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:), &
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                CTildeA(:,:),CTildeB(:,:),   &
                                CA(:,:),CB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

!! set dimensions
!NCholesky = SAPT%NCholesky
!
!! check MCBS/DCBS
!if(A%NBasis.ne.B%NBasis) then
!   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
!   stop
!else
!   NBas = A%NBasis
!endif
!
!allocate(WaBB(NBas,NBas),WbAA(NBas,NBas))
!
!call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
!call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)
!
!allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))
!allocate(ABP0TildeA(A%NDimX,NCholesky), &
!         ABP1TildeA(A%NDimX,NCholesky), &
!         DCholA(A%NDimX,NCholesky))
!
!call prepare_resp_Cmat(A,A1A,A2A,ABP0TildeA,ABP1TIldeA,DCholA,A%NDimX,NCholesky,NBas)
!
!allocate(WbChAA(NCholesky))
!do i=1,NCholesky
!   do ipq=1,A%NDimX
!      ip = A%IndN(1,ipq)
!      iq = A%IndN(2,ipq)
!   
!      WbChAA(i) = WbChAA(i) + WbAA(ip,iq)*DCholA(i,ipq)
!
!   enddo
!enddo
!
!allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))
!allocate(ABP0TildeB(B%NDimX,NCholesky), &
!         ABP1TildeB(B%NDimX,NCholesky), &
!         DCholB(B%NDimX,NCholesky))
!
!call prepare_resp_Cmat(B,A1B,A2B,ABP0TildeB,ABP1TIldeB,DCholB,B%NDimX,NCholesky,NBas)
!
!allocate(CTildeA(A%NDimX,NCholesky),C0TildeA(A%NDimX,NCholesky))
!allocate(CA(A%NDimX,A%NDimX))
!call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
!
!call C_AlphaExpand(CTildeA,C0TildeA,0d0,SAPT%Max_Cn, &
!                   A1A,A2A,ABP0TildeA,ABP1TildeA,    &
!                   A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
!
!print*, 'CTildeA',norm2(CTildeA)
!!call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
!call dgemm('N','N',A%NDimX,A%NDimX,NCholesky,1d0,CTildeA,A%NDimX,DCholA,NCholesky,0d0,CA,A%NDimX)
!
!e2ba = 0
!do ipq=1,A%NDimX
!   ip = A%IndN(1,ipq)
!   iq = A%IndN(2,ipq)
!   do irs=1,A%NDimX
!      ir = A%IndN(1,irs)
!      is = A%IndN(2,irs)
!      e2ba = e2ba + WbAA(ip,iq)*CA(ipq,irs)*WbAA(ir,is)
!   enddo
!enddo
!e2ba = -0.5d0*e2ba
!print*, 'e2ba',e2ba*1000
!
!!call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')
!!call C_AlphaExpand(CTildeB,C0TildeB,OmI,SAPT%Max_Cn,A1B,A2B,ABP0TildeB,ABP1TildeB, &
!!                   A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)
!
!deallocate(DCholB,DCholA)
!deallocate(WbAA,WaBB)
!deallocate(CA)
!deallocate(C0TildeA,CTildeA)
!deallocate(A2B,A2A,A1B,A1A)
!deallocate(ABP0TildeB,ABP0TildeA,ABP1TildeB,ABP1TIldeA)

end subroutine e2ind_CAlphaTilde_block

subroutine e2disp_CAlphaTilde_block(Flags,A,B,SAPT)
!
! calculate 2nd order dispersion energy
! using expansion of C(w) in alpha around alpha=0, up to Max_Cn order
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: nblkA,nblkB
integer :: N,Max_Cn
double precision :: fact,val,val2
double precision :: Cpq,Crs
double precision :: XFactorial,XN1,XN2
double precision :: ACAlpha,OmI,Pi,e2du,e2d
logical :: both

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUS1A(:,:),ABPLUS1B(:,:),&
                                ABMIN1A(:,:),ABMIN1B(:,:),  &
                                A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:),&
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                C1TildeA(:,:),C1TildeB(:,:), &
                                C2TildeA(:,:),C2TildeB(:,:), &
                                CTildeA(:,:), CTildeB(:,:)
double precision,allocatable :: WorkA(:,:),WorkB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

type(EBlockData)             :: LambdaIVA,LambdaIVB
type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)
! test
double precision :: ErrMax
double precision :: Tcpu,Twall

print*, ''
print*, 'Experimental E2disp(CAlpha) procedure...'

! both = coupled + uncoupled
both = SAPT%iCpld

! timing
call clock('START',Tcpu,Twall)

ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPLUS1A(A%NDimX,A%NDimX),ABMIN1A(A%NDimX,A%NDimX))
allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1A
read(iunit) ABMIN1A

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkA,A0BlkIVA,A%IndN,A%CICoef,nblkA,NBas,A%NDimX,'XY0_A')

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1A,A0BlkA,A0BlkIVA,-1d0,.false.,nblkA,A%NDimX)
call add_blk_right(ABMIN1A, A0BlkA,A0BlkIVA,-1d0,.true., nblkA,A%NDimX)

print*, 'ABPLUS1A',norm2(ABPLUS1A)
print*, 'ABMIN1A ',norm2(ABMIN1A)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1A, A1A,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1A,A1A,1.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,A%NDimX,'X')
print*, 'A1A',norm2(A1A)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A2A,A%NDimX)
deallocate(ABMIN1A)
print*, 'A2A',norm2(A2A)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(A%DChol),ABP0TildeA,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde',norm2(ABP0TildeA)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeA(A%NDimX,NCholesky))
Call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS1A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP1TildeA,A%NDimX)
print*, 'APLUS1Tilde',norm2(ABP1TildeA)

deallocate(ABPLUS1A)
call release_ac0block(A0BlkA,A0BlkIVA,nblkA)
deallocate(A0BlkA)

! monomer B
allocate(ABPLUS1B(B%NDimX,B%NDimX),ABMIN1B(B%NDimX,B%NDimX))
allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))

!allocate(ABPMB(B%NDimX,B%NDimX),ABPTildeB(B%NDimX,NCholesky))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1B
read(iunit) ABMIN1B

close(iunit)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkB,A0BlkIVB,B%IndN,B%CICoef,nblkB,NBas,B%NDimX,'XY0_B')

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1B,A0BlkB,A0BlkIVB,-1d0,.false.,nblkB,B%NDimX)
call add_blk_right(ABMIN1B, A0BlkB,A0BlkIVB,-1d0,.true., nblkB,B%NDimX)

print*, 'ABPLUS1B',norm2(ABPLUS1B)
print*, 'ABMIN1B ',norm2(ABMIN1B)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call ABPM_HALFTRAN_GEN_L(ABMIN1B, A1B,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1B,A1B,1.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,B%NDimX,'X')
print*, 'A1B',norm2(A1B)

!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A2B,B%NDimX)
deallocate(ABMIN1B)
print*, 'A2B',norm2(A2B)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeB(B%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(B%DChol),ABP0TildeB,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeB(B%NDimX,NCholesky))
Call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS1B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP1TildeB,B%NDimX)
print*, 'APLUS1Tilde-b',norm2(ABP1TildeB)

deallocate(ABPLUS1B)
call release_ac0block(A0BlkB,A0BlkIVB,nblkB)
deallocate(A0BlkB)

! get CAlphaTilde in an iterative manner

NFreq = 12
Max_Cn = SAPT%Max_Cn
print*, ''
print*, 'CAlphaTilde: '
print*, 'SAPT%NFreq =', NFreq
print*, 'SAPT%MaxCn =', Max_Cn

allocate(XFreq(NFreq),WFreq(NFreq))

allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(C1TildeA(A%NDimX,NCholesky),C1TildeB(B%NDimX,NCholesky))
allocate(C2TildeA(A%NDimX,NCholesky),C2TildeB(B%NDimX,NCholesky))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(WorkA(A%NDimX,NCholesky),WorkB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2d  = 0
e2du = 0
ErrMax = 0d0
do ifreq=1,NFreq

   OmI = XFreq(ifreq)
   print*,'OmI', OmI

   !if(both) then

      ! coupled
      print*, 'ErrMax', ErrMax
      call C_AlphaExpand(CTildeA,C0TildeA,OmI,Max_Cn,A1A,A2A,ABP0TildeA,ABP1TildeA, &
                         A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX,ErrMax,ifreq)
      call C_AlphaExpand(CTildeB,C0TildeB,OmI,Max_Cn,A1B,A2B,ABP0TildeB,ABP1TildeB, &
                         A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX,ErrMax,ifreq)

      call dgemm('T','T',NCholesky,NCholesky,A%NDimX,1d0,CTildeA,A%NDimX,A%DChol,NCholesky,0d0,CA,NCholesky)
      call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

      val = 0
      do j=1,NCholesky
         do i=1,NCholesky
            val = val + CA(i,j)*CB(i,j)
         enddo
      enddo

      e2d = e2d + WFreq(ifreq)*val
   !endif

   !! uncoupled

   !call C_AlphaExpand_unc(C0TildeA,OmI,A1A,A2A,ABP0TildeA,ABP1TildeA, &
   !                   A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
   !call C_AlphaExpand_unc(C0TildeB,OmI,A1B,A2B,ABP0TildeB,ABP1TildeB, &
   !                   A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)

   !call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,DCholA,NCholesky,C0TildeA,A%NDimX,0d0,CA,NCholesky)
   !call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,DCholB,NCholesky,C0TildeB,B%NDimX,0d0,CB,NCholesky)

   !val2 = 0
   !do j=1,NCholesky
   !   do i=1,NCholesky
   !      val2 = val2+ CA(j,i)*CB(i,j)
   !   enddo
   !enddo
   !e2du = e2du + WFreq(ifreq)*val2

enddo

!SAPT%e2disp_unc = -32d0/Pi*e2du
!e2du = -32d0/Pi*e2du*1d3
!call print_en('E2disp(Alph,unc)',e2du,.false.)
!if(both) then
   SAPT%e2disp  = -32d0/Pi*e2d
   e2d = -32d0/Pi*e2d*1d3
   !print*, 'e2d = ',e2d
   call print_en('E2disp(CAlpha)',e2d,.false.)
!endif

call clock('E2disp(CAlpha)',Tcpu,Twall)

deallocate(WFreq,XFreq)
deallocate(A2A,A1A)
deallocate(ABP1TildeA,ABP0TildeA)
deallocate(CTildeB,CTildeA)
deallocate(C2TildeB,C2TildeA)
deallocate(C1TildeB,C1TildeA)
deallocate(C0TildeB,C0TildeA)
deallocate(WorkB,WorkA)
deallocate(CB,CA)
!deallocate(B%DCholB,A%DCholA)

end subroutine e2disp_CAlphaTilde_block

subroutine e2disp_CAlphaTilde_full(Flags,A,B,SAPT)
!
! this is CAlphaTilde procedure which uses
! NDimX*NDimX Lambda (no block structure)
! THIS IS HERE FOR TESTING ONLY
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBas
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
integer :: nblkA,nblkB
integer :: N,Max_Cn
double precision :: fact,val
double precision :: Cpq,Crs
double precision :: XFactorial,XN1,XN2
double precision :: ACAlpha,OmI,Pi,e2d

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUS0A(:,:),ABPLUS0B(:,:),&
                                ABMIN0A(:,:),ABMIN0B(:,:),  &
                                ABPLUS1A(:,:),ABPLUS1B(:,:),&
                                ABMIN1A(:,:),ABMIN1B(:,:),  &
                                A1A(:,:),A1B(:,:),&
                                A2A(:,:),A2B(:,:),&
                                ABP0TildeA(:,:),ABP0TildeB(:,:),&
                                ABP1TildeA(:,:),ABP1TildeB(:,:)
!double precision,allocatable :: DCholA(:,:),DCholB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:), &
                                C1TildeA(:,:),C1TildeB(:,:), &
                                C2TildeA(:,:),C2TildeB(:,:), &
                                CTildeA(:,:), CTildeB(:,:)
double precision,allocatable :: LambdaA(:,:),LambdaB(:,:)
double precision,allocatable :: WorkA(:,:),WorkB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

!type(EBlockData)             :: LambdaIVA,LambdaIVB
!type(EBlockData),allocatable :: LambdaA(:),LambdaB(:)


ACAlpha = 1.0d0
Pi = 4.0d0*atan(1.0)

! set dimensions
NCholesky = SAPT%NCholesky

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

!! get Dmat
!allocate(DCholA(NCholesky,A%NDimX),DCholB(NCholesky,B%NDimX))
!
!DCholA = 0
!do j=1,A%NDimX
!   ip = A%IndN(1,j)
!   iq = A%IndN(2,j)
!   ipq = iq + (ip-1)*NBas
!   Cpq = A%CICoef(ip) + A%CICoef(iq)
!   do i=1,NCholesky
!      DCholA(i,j) = Cpq*A%FF(i,ipq)
!   enddo
!enddo
!DCholB = 0
!do j=1,B%NDimX
!   ir = B%IndN(1,j)
!   is = B%IndN(2,j)
!   irs = is + (ir-1)*NBas
!   Crs = B%CICoef(ir) + B%CICoef(is)
!   do i=1,NCholesky
!      DCholB(i,j) = Crs*B%FF(i,irs)
!   enddo
!enddo

! monomer A
allocate(ABPLUS0A(A%NDimX,A%NDimX),ABMIN0A(A%NDimX,A%NDimX))
allocate(ABPLUS1A(A%NDimX,A%NDimX),ABMIN1A(A%NDimX,A%NDimX))
allocate(A1A(A%NDimX,A%NDimX),A2A(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1A
read(iunit) ABMIN1A

close(iunit)

open(newunit=iunit,file='A0MAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS0A
read(iunit) ABMIN0A

close(iunit)

! AB1 = AB1 - A0
ABPLUS1A = ABPLUS1A - ABPLUS0A
ABMIN1A  = ABMIN1A  - ABMIN0A

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS0A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A1A,A%NDimX)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN0A,A%NDimX,1d0,A1A,A%NDimX)
print*, 'A1A',norm2(A1A)
!
!Calc: A2 = ABP1*ABM1
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUS1A,A%NDimX,ABMIN1A,A%NDimX,0.0d0,A2A,A%NDimX)
deallocate(ABMIN1A)
print*, 'A2A',norm2(A2A)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS0A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP0TildeA,A%NDimX)
print*, 'APLUS0Tilde',norm2(ABP0TildeA)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeA(A%NDimX,NCholesky))
call dgemm('N','T',A%NDimX,NCholesky,A%NDimX,1d0,ABPLUS1A,A%NDimX,A%DChol,NCholesky,0.0d0,ABP1TildeA,A%NDimX)
print*, 'APLUS1Tilde',norm2(ABP1TildeA)

deallocate(ABPLUS1A)
!
! monomer B
allocate(ABPLUS0B(B%NDimX,B%NDimX),ABMIN0B(B%NDimX,B%NDimX))
allocate(ABPLUS1B(B%NDimX,B%NDimX),ABMIN1B(B%NDimX,B%NDimX))
allocate(A1B(B%NDimX,B%NDimX),A2B(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS1B
read(iunit) ABMIN1B

close(iunit)

open(newunit=iunit,file='A0MAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUS0B
read(iunit) ABMIN0B

close(iunit)

! AB1 = AB1 - A0
ABPLUS1B = ABPLUS1B - ABPLUS0B
ABMIN1B  = ABMIN1B  - ABMIN0B

print*, 'ABPLUS1B',norm2(ABPLUS1B)
print*, 'ABMIN1B ',norm2(ABMIN1B)

!Calc: A1 = ABP0*ABM1+ABP1*ABM0
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS0B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A1B,B%NDimX)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN0B,B%NDimX,1d0,A1B,B%NDimX)
print*, 'A1B',norm2(A1B)
!
!Calc: A2 = ABP1*ABM1
Call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUS1B,B%NDimX,ABMIN1B,B%NDimX,0.0d0,A2B,B%NDimX)
deallocate(ABMIN1B)
print*, 'A2B',norm2(A2B)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeB(B%NDimX,NCholesky))
call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS0B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP0TildeB,B%NDimX)
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(ABP1TildeB(B%NDimX,NCholesky))
Call dgemm('N','T',B%NDimX,NCholesky,B%NDimX,1d0,ABPLUS1B,B%NDimX,B%DChol,NCholesky,0.0d0,ABP1TildeB,B%NDimX)
print*, 'APLUS1Tilde-b',norm2(ABP1TildeB)

deallocate(ABPLUS1B)
!
! get CAlphaTilde in an iterative manner

NFreq = 12
Max_Cn = SAPT%Max_Cn
print*, 'SAPT%NFreq =', NFreq
print*, 'SAPT%MaxCn =', Max_Cn

allocate(XFreq(NFreq),WFreq(NFreq))
!
allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(C1TildeA(A%NDimX,NCholesky),C1TildeB(B%NDimX,NCholesky))
allocate(C2TildeA(A%NDimX,NCholesky),C2TildeB(B%NDimX,NCholesky))
allocate(CTildeA(A%NDimX,NCholesky),CTildeB(B%NDimX,NCholesky))
allocate(LambdaA(A%NDimX,A%NDimX),LambdaB(B%NDimX,B%NdimX))
allocate(WorkA(A%NDimX,NCholesky),WorkB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2d = 0
do ifreq=NFreq,1,-1

   OmI = XFreq(ifreq)

!  Calc: LAMBDA=(A0+Om^2)^-1
   Call Inv_AC0Blk(OmI**2,LambdaA,A0BlkA,A0BlkIVA,nblkA,A%NDimX)
   Call Inv_AC0Blk(OmI**2,LambdaB,A0BlkB,A0BlkIVB,nblkB,B%NDimX)

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,0.5d0,LambdaA,A%NDimX,ABP0TildeA,A%NDimX,0.0d0,C0TildeA,A%NDimX)
   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,0.5d0,LambdaB,B%NDimX,ABP0TildeB,B%NDimX,0.0d0,C0TildeB,B%NDimX)

   !print*, 'C0Tilde-A',OmI,norm2(C0TildeA)
   !print*, 'C0Tilde-B',OmI,norm2(C0TildeB)

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,A1A,A%NDimX,C0TildeA,A%NDimX,0.0d0,CTildeA,A%NDimX)
   CTildeA = 0.5d0*ABP1TildeA - CTildeA
   call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.d0,LambdaA,A%NDimX,CTildeA,A%NDimX,0.0d0,C1TildeA,A%NDimX)

   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,A1B,B%NDimX,C0TildeB,B%NDimX,0.0d0,CTildeB,B%NDimX)
   CTildeB = 0.5d0*ABP1TildeB - CTildeB
   call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.d0,LambdaB,B%NDimX,CTildeB,B%NDimX,0.0d0,C1TildeB,B%NDimX)

   !print*, 'C1Tilde-A',OmI,norm2(C1TildeA)
   !print*, 'C1Tilde-B',OmI,norm2(C1TildeB)

   ! test uncoupled
   CTildeA = 0
   CTildeB = 0

   CTildeA = C0TildeA
   CTildeB = C0TildeB

   ! test semicoupled
   CTildeA = CTildeA + C1TildeA
   CTildeB = CTildeB + C1TildeB

   XFactorial = 1
   do N=2,Max_Cn

       XFactorial = XFactorial*N
       XN1 = -N
       XN2 = -N*(N-1)

       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN2,A2A,A%NDimX,C0TildeA,A%NDimX,0.0d0,WorkA,A%NDimX)
       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,XN1,A1A,A%NDimX,C1TildeA,A%NDimX,1.0d0,WorkA,A%NDimX)
       call dgemm('N','N',A%NDimX,NCholesky,A%NDimX,1.0d0,LambdaA,A%NDimX,WorkA,A%NDimX,0.0d0,C2TildeA,A%NDimX)

       CTildeA  = CTildeA + C2TildeA / XFactorial
       C0TildeA = C1TildeA
       C1TildeA = C2TildeA

       !print*, 'N,COMTildeA',N,norm2(CTildeA)

       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN2,A2B,B%NDimX,C0TildeB,B%NDimX,0.0d0,WorkB,B%NDimX)
       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,XN1,A1B,B%NDimX,C1TildeB,B%NDimX,1.0d0,WorkB,B%NDimX)
       call dgemm('N','N',B%NDimX,NCholesky,B%NDimX,1.0d0,LambdaB,B%NDimX,WorkB,B%NDimX,0.0d0,C2TildeB,B%NDimX)

       CTildeB  = CTildeB + C2TildeB / XFactorial
       C0TildeB = C1TildeB
       C1TildeB = C2TildeB

       !print*, 'N,COMTildeB',N,norm2(CTildeB)

   enddo

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,CTildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,CTildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo

   e2d = e2d + WFreq(ifreq)*val

enddo

e2d = -32d0/Pi*e2d*1d3
!print*, '???e2d = ',e2d
call print_en('E2disp(CAlpha)-TEST',e2d,.false.)

deallocate(WFreq,XFreq)
deallocate(A2B,A1B,A2A,A1A)
deallocate(ABP1TildeA,ABP0TildeA,ABP1TildeB,ABP0TildeB)
deallocate(LambdaB,LambdaA)
deallocate(CTildeB,CTildeA)
deallocate(C2TildeB,C2TildeA)
deallocate(C1TildeB,C1TildeA)
deallocate(C0TildeB,C0TildeA)
deallocate(WorkB,WorkA)
deallocate(CB,CA)
!deallocate(DCholB,DCholA)

end subroutine e2disp_CAlphaTilde_full

subroutine e2disp_CAlphaTilde_unc(Flags,A,B,SAPT)
!
! calculate 2nd order uncoupled dispersion energy
! using expansion of C(w) in alpha around alpha=0, up to Max_Cn order
! use A0 blocks (not diagonal)
! with Cholesky vectors
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBasis,NCholesky
integer :: nblkA,nblkB
integer :: i,j
integer :: ifreq,NFreq

double precision :: Pi
double precision :: OmI,val
double precision :: e2du

double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: ABP0TildeA(:,:),ABP0TildeB(:,:)
double precision,allocatable :: C0TildeA(:,:),C0TildeB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)

type(EBlockData)             :: A0BlkIVA,A0BlkIVB
type(EBlockData),allocatable :: A0BlkA(:),A0BlkB(:)

Pi = 4.0d0*atan(1.0)

! set dimensions
NBasis = A%NBasis
NCholesky = SAPT%NCholesky

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkA,A0BlkIVA,A%IndN,A%CICoef,nblkA,NBasis,A%NDimX,'XY0_A')

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(ABP0TildeA(A%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(A%DChol),ABP0TildeA,0.0d0,A0BlkA,A0BlkIVA,nblkA,A%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-a',norm2(ABP0TildeA)

call release_ac0block(A0BlkA,A0BlkIVA,nblkA)
deallocate(A0BlkA)

! get A0PLUS, A0MIN in blocks matY and matX
call Sblock_to_ABMAT(A0BlkB,A0BlkIVB,B%IndN,B%CICoef,nblkB,NBasis,B%NDimX,'XY0_B')

allocate(ABP0TildeB(B%NDimX,NCholesky))
call ABPM_HALFTRAN_GEN_L(transpose(B%DChol),ABP0TildeB,0.0d0,A0BlkB,A0BlkIVB,nblkB,B%NDimX,NCholesky,'Y')
print*, 'APLUS0Tilde-b',norm2(ABP0TildeB)

call release_ac0block(A0BlkB,A0BlkIVB,nblkB)
deallocate(A0BlkB)

! get CAlphaTilde in 0th order
NFreq = 12
print*, ''
print*, 'E2disp, unc from CAlphaTilde(0)'
print*, 'SAPT%NFreq =', NFreq

allocate(XFreq(NFreq),WFreq(NFreq))
allocate(C0TildeA(A%NDimX,NCholesky),C0TildeB(B%NDimX,NCholesky))
allocate(CA(NCholesky,NCholesky),CB(NCholesky,NCholesky))

call FreqGrid(XFreq,WFreq,NFreq)

! read ABPLUS0.ABMIN0 blocks
call read_ABPM0Block(A0BlkA,A0BlkIVA,nblkA,'A0BLK_A')
call read_ABPM0Block(A0BlkB,A0BlkIVB,nblkB,'A0BLK_B')

e2du = 0
do ifreq=1,NFreq

   OmI = XFreq(ifreq)

   call C_AlphaExpand_unc(C0TildeA,OmI,ABP0TildeA, &
                      A0BlkA,A0BlkIVA,nblkA,NCholesky,A%NDimX)
   call C_AlphaExpand_unc(C0TildeB,OmI,ABP0TildeB, &
                      A0BlkB,A0BlkIVB,nblkB,NCholesky,B%NDimX)

   call dgemm('N','N',NCholesky,NCholesky,A%NDimX,1d0,A%DChol,NCholesky,C0TildeA,A%NDimX,0d0,CA,NCholesky)
   call dgemm('N','N',NCholesky,NCholesky,B%NDimX,1d0,B%DChol,NCholesky,C0TildeB,B%NDimX,0d0,CB,NCholesky)

   val = 0
   do j=1,NCholesky
      do i=1,NCholesky
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo
   e2du = e2du + WFreq(ifreq)*val
   write(6, '(1x,"OmI= ",F12.6,", E2d,unc=",F12.6)') OmI,-32d0/Pi*e2du*1d3

enddo

SAPT%e2disp_unc = -32d0/Pi*e2du
e2du = -32d0/Pi*e2du*1d3

call print_en('E2disp,unc,C(0)',e2du,.true.)

deallocate(CB,CA)
deallocate(C0TildeB,C0TildeA)
deallocate(ABP0TildeB,ABP0TildeA)
deallocate(WFreq,XFreq)

end subroutine e2disp_CAlphaTilde_unc


end module sapt_Chol_pol
