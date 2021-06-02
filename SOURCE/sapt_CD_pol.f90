module sapt_CD_pol
use types
use sapt_utils

implicit none

contains

subroutine e1elst_Chol(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i,j,ii,jj
integer :: NBas,NCholesky
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Vabb(:,:),Vbaa(:,:)
double precision :: ea,eb,eab,elst
double precision,external  :: ddot

! set dimensions
 NBas = A%NBasis 
 NCholesky = SAPT%NCholesky

 allocate(Va(NBas,NBas),Vb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas))

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

! sum_p n_p v^B_pp 
 ea = 0
 do i=1,A%num0+A%num1
    ea = ea + A%Occ(i)*Vbaa(i,i) 
 enddo
 ea = 2d0*ea
 !print*, 'ea',ea

! sum_q n_q v^A_qq 
 eb = 0
 do j=1,B%num0+B%num1
    eb = eb + B%Occ(j)*Vabb(j,j) 
 enddo
 eb = 2d0*eb
 !print*, 'eb',eb

! sum_pq n_p n_q v_{pq}^{pq}
 eab = 0
 do j=1,B%num0+B%num1
    jj = (j-1)*NBas+j
    do i=1,A%num0+A%num1
       ii = (i-1)*NBas+i
       eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%FF(:,ii),1,B%FF(:,jj),1)
       !eab = eab + A%Occ(i)*B%Occ(j)*ddot(NCholesky,A%FO(:,ii),1,B%FO(:,jj),1)
    enddo
 enddo
 eab = 4d0*eab
 !print*, 'eab', eab

 elst = ea + eb + eab + SAPT%Vnn 

 call print_en('V_nn',SAPT%Vnn,.false.)
 call print_en('Eelst',elst*1000,.false.)
 SAPT%elst = elst

 deallocate(Vb,Va,Vbaa,Vabb)

end subroutine e1elst_Chol

subroutine e2disp_Chol(Flags,A,B,SAPT)
! calculate 2nd order dispersion energy
! in coupled and uncoupled approximations
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

subroutine e2disp_Cmat(Flags,A,B,SAPT)
! calculate 2nd order dispersion energy
! use C(omega)
! in coupled and uncoupled approximations
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: ifreq,NFreq,NCholesky
integer :: i,j,ipq,irs
integer :: ip,iq,ir,is
integer :: iunit,info
double precision :: fact,val
double precision :: Omega,Pi,e2d
double precision,allocatable :: ABPMA(:,:),ABPMB(:,:),&
                                ABPLUSA(:,:),ABPLUSB(:,:)
double precision,allocatable :: CA(:,:),CB(:,:)
double precision,allocatable :: XFreq(:),WFreq(:)
double precision,allocatable :: work(:,:),ints(:)

double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

! check MCBS/DCBS
if(A%NBasis.ne.B%NBasis) then
   write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
   stop
else
   NBas = A%NBasis
endif

! set dimensions
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

Pi = 4.0d0*atan(1.0)

! get ABPM = ABPLUS.ABMIN
allocate(ABPMA(A%NDimX,A%NDimX),ABPMB(B%NDimX,B%NDimX),&
         ABPLUSA(A%NDimX,A%NDimX),ABPLUSB(B%NDimX,B%NDimX))
allocate(work(A%NDimX,A%NDimX))

open(newunit=iunit,file='ABMAT_A',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSA
read(iunit) work

close(iunit)
call dgemm('N','N',A%NDimX,A%NDimX,A%NDimX,1d0,ABPLUSA,A%NDimX,work,A%NDimX,0d0,ABPMA,A%NDimX)

deallocate(work)
allocate(work(B%NDimX,B%NDimX))

open(newunit=iunit,file='ABMAT_B',status='OLD',&
     access='SEQUENTIAL',form='UNFORMATTED')

read(iunit) ABPLUSB
read(iunit) work

close(iunit)
call dgemm('N','N',B%NDimX,B%NDimX,B%NDimX,1d0,ABPLUSB,B%NDimX,work,B%NDimX,0d0,ABPMB,B%NDimX)

deallocate(work)

!print*, 'ABPLUS-A',norm2(ABPLUSA(1:A%NDimX,1:A%NDimX))
!print*, 'ABPLUS-B',norm2(ABPLUSB)
!print*, 'ABPMA   ',norm2(ABPMA)
!print*, 'ABPMB   ',norm2(ABPMB)

! frequency integration
NFreq = 18
allocate(XFreq(NFreq),WFreq(NFreq))

call FreqGrid(XFreq,WFreq,NFreq)

allocate(CB(B%NDimX,B%NDimX))
allocate(ints(NBas**2),work(A%NDimX,B%NDimX))

e2d = 0
do ifreq=1,NFreq

   Omega = XFreq(ifreq)

   allocate(CA(A%NDimX,A%NDimX))

   !print*, ifreq,WFreq(ifreq),Omega

   if(Flags%ICholesky==1) then
     print*, 'Cholesky Cmat not ready yet! Aborting...'
     !call get_Cmat(CA,A%CICoef,A%IndN,ABPMA,ABPLUSA,Omega,A%NDimX,NBas)
     return
   else
      call get_Cmat(CA,A%CICoef,A%IndN,ABPMA,ABPLUSA,Omega,A%NDimX,NBas)
      call get_Cmat(CB,B%CICoef,B%IndN,ABPMB,ABPLUSB,Omega,B%NDimX,NBas)
   endif

   open(newunit=iunit,file='TWOMOAB',status='OLD',&
        access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

   ints = 0
   work = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)

         do i=1,A%NDimX
            work(i,irs) = work(i,irs) + fact*CA(ipq,i)
         enddo

      enddo
   enddo

   deallocate(CA)
   allocate(CA(B%NDimX,B%NDimX))
   CA = 0
   ints = 0
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) ints(1:nOVB)

      do irs=1,B%NDimX
         ir = B%IndN(1,irs)
         is = B%IndN(2,irs)

         fact = ints(is+(ir-B%num0-1)*dimOB)
         do j=1,B%NDimX
            CA(irs,j) = CA(irs,j) + fact*work(ipq,j)
         enddo

      enddo

   enddo

   val = 0
   do j=1,B%NDimX
      do i=1,B%NDimX
         val = val + CA(j,i)*CB(i,j)
      enddo
   enddo
   e2d = e2d + WFreq(ifreq)*val

   close(iunit)
   deallocate(CA)

enddo ! NFreq

SAPT%e2disp = -8d0/Pi*e2d
e2d = -8d0/Pi*e2d*1d3

call print_en('E2disp',e2d,.true.)

deallocate(ints,work)
deallocate(WFreq,XFreq)
deallocate(CB)
deallocate(ABPMB,ABPMA,ABPLUSB,ABPLUSA)

end subroutine e2disp_Cmat

subroutine get_Cmat(Cmat,CICoef,IndN,ABPM,ABPLUS,Omega,NDimX,NBas)
implicit none

integer,intent(in) :: NDimX,NBas,IndN(2,NBas)
double precision,intent(in)  :: CICoef(NBas),Omega
double precision,intent(in)  :: ABPM(NDimX,NDimX),ABPLUS(NDimX,NDimX)
double precision,intent(out) :: Cmat(NDimX,NDimX)

integer          :: info,lwork
integer          :: i,ipq,irs,ip,iq,ir,is
double precision :: fact_rs,fact_pq
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:,:),work1(:)

allocate(work(NDimX,NDimX),ipiv(NDimX))
work = ABPM
do i=1,NDimX
   work(i,i) = work(i,i) + Omega**2
enddo

! work=(ABPM+omega^2)-1.ABPLUS
!...

! this is slower
!call dgetrf(NDimX,NDimX,work,NDimX,ipiv,info)
!print*, 'info1',info
!allocate(work1(1))
!lwork = -1
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!lwork = int(work1(1))
!print*, 'lwork:',lwork
!deallocate(work1)
!allocate(work1(lwork))
!call dgetri(NDimX,work,NDimX,ipiv,work1,lwork,info)
!print*, 'info2',info
!!
!! Cmat=work.ABPLUS
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,work,NDimX,ABPLUS,NDimX,0d0,Cmat,NDimX)
!
! deallocate(work1)
!
Cmat = ABPlus
call dgesv(NDimX,NDimX,work,NDimX,ipiv,Cmat,NDimX,info)
!print*, 'info',info
!print*, 'Cmat',norm2(Cmat)

!
do irs=1,NDimX
   ir = IndN(1,irs)
   is = IndN(2,irs)
   fact_rs = CICoef(ir)+CICoef(is)
   do ipq=1,NDimX
      ip = IndN(1,ipq)
      iq = IndN(2,ipq)
      fact_pq = CICoef(ip)+CICoef(iq)
      Cmat(ipq,irs) = fact_rs*fact_pq*Cmat(ipq,irs)
   enddo
enddo

deallocate(ipiv,work)

end subroutine get_Cmat

end module sapt_CD_pol