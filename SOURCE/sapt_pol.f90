module sapt_pol
use types
use tran
use sapt_utils

implicit none

contains

subroutine e1elst(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j
integer :: NBas
double precision,allocatable :: PA(:,:),PB(:,:) 
double precision,allocatable :: Va(:,:),Vb(:,:),Ja(:,:) 
double precision,allocatable :: work(:,:)
double precision :: tmp,ea,eb,elst
double precision,parameter :: Half=0.5d0
double precision,external  :: trace

! set dimensions
 NBas = A%NBasis 

 allocate(PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),Ja(NBas,NBas))
 allocate(work(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,2d0,PA)
 call get_den(NBas,B%CMO,B%Occ,2d0,PB)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call make_J1(NBas,PA,Ja,'AOTWOSORT')

! Tr[Pa.Va + Pb.Vb + Pb.Ja]
 work=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Vb,NBas,0d0,work,NBas)
 ea = trace(work,NBas)
! print*, ea
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Va,NBas,0d0,work,NBas)
 eb = trace(work,NBas) 
! print*, eb
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Ja,NBas,0d0,work,NBas)
 ea = ea + trace(work,NBas)
! print*, trace(work,NBas) 
 elst = ea + eb + SAPT%Vnn 

 write(LOUT,'(1x,a,f16.8)') 'V_nn        = ', SAPT%Vnn
 write(LOUT,'(1x,a,f16.8,/)') 'Eelst       = ', elst*1000d0
 SAPT%elst = elst

!! test wabb
! work = 0
! Ja = 0
! call make_J1(NBas,PB,Ja)
!
! call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Ja,NBas,0d0,work,NBas)
! print*, 'Pa.Jb', trace(work,NBas)
!
! call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,B%WPot,NBas,0d0,work,NBas)
! print*, 'Pa.Wb', trace(work,NBas)
!

 deallocate(work)
 deallocate(Ja,Vb,Va,PB,PA) 

end subroutine e1elst

subroutine e2ind_hf_icphf(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas
double precision :: e2ba,e2ab,e2ic

 NBas = A%NBasis 

 call solve_cphf(A,B%WPot,e2ba,Flags,NBas)
 call solve_cphf(B,A%WPot,e2ab,Flags,NBas)

 e2ic = (e2ab + e2ba)

 write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0 
 write(LOUT,'(1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000d0 
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 
 SAPT%e2ind = e2ic

end subroutine e2ind_hf_icphf

subroutine e2ind_resp(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
! test
integer :: info
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

 !! TEST ICPKS
 !write(LOUT,*) 'A%PP?', norm2(A%PP)
 !write(LOUT,*) 'B%PP?', norm2(B%PP)

 !allocate(wtest(A%NDimX))
 !allocate(Work(NBas**2),ipiv(A%NDimX))
 !wtest = 0 
 !do pq=1,A%NDimX
 !   ip = A%IndN(1,pq)
 !   iq = A%IndN(2,pq)
 !   wtest(pq) = WbAA(ip,iq)
 !enddo
 !call dsysv('U',A%NDimX,1,A%PP,A%NDimX,ipiv,wtest,A%NDimX,Work,NBas,info)
 !e2ba = 0
 !do pq=1,A%NDimX
 !   ip = A%IndN(1,pq)
 !   iq = A%IndN(2,pq)
 !   e2ba = e2ba + wtest(pq)*WbAA(ip,iq)
 !enddo
 !print*, 'e2ba ',2d0*e2ba*1000
 !e2ba = 0
 !
 !do rs=1,B%NDimX
 !   ir = B%IndN(1,rs)
 !   is = B%IndN(2,rs)
 !   wtest(rs) = WaBB(ir,is)
 !enddo
 !call dsysv('U',B%NDimX,1,B%PP,B%NDimX,ipiv,wtest,B%NDimX,Work,NBas,info)
 !e2ab = 0
 !do rs=1,B%NDimX
 !   ir = B%IndN(1,rs)
 !   is = B%IndN(2,rs)
 !   e2ab = e2ab + wtest(rs)*WaBB(ir,is)
 !enddo
 !print*, 'e2ab ',2d0*e2ab*1000
 !e2ab = 0
 !
 !deallocate(wtest,Work)
 !deallocate(A%PP,B%PP)
 !!
 !! end TEST-ICPKS

 !tmp=0
 !do i=1,A%NDimX
 !do j=1,A%NDimX
 !   tmp = tmp + AlphaA(i,j)**2
 !enddo
 !enddo
 !print*, 'RMA: ',tmp
 !tmp=0
 !do i=1,B%NDimX
 !do j=1,B%NDimX
 !   tmp = tmp + AlphaB(i,j)**2
 !enddo
 !enddo
 !print*, 'RMB: ',tmp

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
 write(LOUT,'(/,1x,a,f16.8)') 'Ind(B-->A)   = ', e2ba*1000d0 

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
 write(LOUT,'(1x,a,f16.8)') 'Ind(A-->B)   = ', e2ab*1000d0 

 e2ic = (e2ab + e2ba)
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 
 SAPT%e2ind = e2ic

 deallocate(WaBB,WbAA,AlphaB,AlphaA)
 deallocate(OmB,EVecB,Oma,EVecA)

end subroutine e2ind_resp

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
double precision,parameter :: BigE   = 1.D10
double precision,parameter :: SmallE = 1.D-6

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

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
 write(LOUT,'(/1x,a,f16.8)') 'Ind(B<--A)     = ', e2ab*1000d0

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
 write(LOUT,'(1x,a,f16.8)') 'Ind(A<--B)     = ', e2ba*1000d0

 e2ic = (e2ab + e2ba)
 SAPT%e2ind = e2ic
 write(LOUT,'(1x,a,f16.8)') 'E2ind          = ', e2ic*1000d0

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

    write(LOUT,'(/1x,a,f16.8)') 'Ind(B<--A,unc) = ', e2ab_unc*1000d0
    write(LOUT,'(1x,a,f16.8)')  'Ind(A<--B,unc) = ', e2ba_unc*1000d0
    write(LOUT,'(1x,a,f16.8)')  'E2ind(unc)     = ', e2ic_unc*1000d0

    do i=1,A%NDimX
       associate(Y => Y01BlockA(i))
         deallocate(Y%vec0)
       end associate
    enddo
    deallocate(Y01BlockA)
    deallocate(OmA0,tmp01)
 endif

 ! calculate deexcitations
 !if(SAPT%DispExc) call e2ind_dexc(Flags,A,B,SAPT)

 deallocate(tmpA)
 deallocate(WaBB,WbAA)
 deallocate(EVecA,OmA)

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

 !call calc_resp_apsg(EVecA,OmA,AlphaA,0d0,A)
 !call calc_resp_apsg(EVecB,OmB,AlphaB,0d0,B)

! termsBA=0
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!    do rs=1,A%NDimX
!       ir = A%IndN(1,rs)
!       is = A%IndN(2,rs)
!
!       termsBA(1) = termsBA(1) + & 
!            WbAA(ip,iq)*AlphaA(pq,rs)*WbAA(ir,is)
!
!    enddo
! enddo
! termsBA(1) = -0.5d0*termsBA(1)
!! write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A) p>q r>s   = ', termsBA(1)*1000 
!
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!    do ir=1,A%NDimN
!
!      termsBA(2) = termsBA(2) + & 
!            WbAA(ip,iq)*AlphaA(pq,A%NDimX+ir)*WbAA(ir,ir)
!
!    enddo
! enddo
! termsBA(2) = -1d0*termsBA(2)
!
! do ip=1,A%NDimN
!    do ir=1,A%NDimN
!
!       termsBA(3) = termsBA(3) + & 
!            WbAA(ip,ip)*AlphaA(A%NDimX+ip,A%NDimX+ir)*WbAA(ir,ir)
!
!    enddo
! enddo
! termsBA(3) = -0.5d0*termsBA(3)
!
! e2ba = sum(termsBA)
! write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000 
!
!!print*, 'terms BA',termsBA(1)*1d3
! !print*,''
! print*,'ADimEx,BDimEx',ADimEx,BDimEx
! print*,'ANDimX,BNDimX',A%NDimX,B%NDimX
! print*,'NDimN,BNDimN',A%NDimN,B%NDimN
! !!end test

! termsAB=0
! do pq=1,B%NDimX
!    ip = B%IndN(1,pq)
!    iq = B%IndN(2,pq)
!    do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!
!       termsAB(1) = termsAB(1) + & 
!            WaBB(ip,iq)*AlphaB(pq,rs)*WaBB(ir,is)
!
!
!    enddo
! enddo
! termsAB(1) = -0.5d0*termsAB(1)
!! write(LOUT,'(1x,a,f16.8)') 'Ind(A--B) p>q r>s  = ', termsAB(1)*1000d0 
!
! do pq=1,B%NDimX
!    ip = B%IndN(1,pq)
!    iq = B%IndN(2,pq)
!    do ir=1,B%NDimN
!
!       termsAB(2) = termsAB(2) + & 
!            WaBB(ip,iq)*AlphaB(pq,B%NDimX+ir)*WaBB(ir,ir)
!
!    enddo
! enddo
! termsAB(2) = -1.0d0*termsAB(2)
!
! do ip=1,B%NDimN
!    do ir=1,B%NDimN
!
!       termsAB(3) = termsAB(3) + & 
!            WaBB(ip,ip)*AlphaB(B%NDimX+ip,B%NDimX+ir)*WaBB(ir,ir)
!
!    enddo
! enddo
! termsAB(3) = -0.5d0*termsAB(3)
!
! e2ab = sum(termsAB)
! write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0 
!
! e2tmp = (termsBA(1)+termsAB(1))
! write(LOUT,'(1x,a,f16.8)') 'p>q r>s     = ', e2tmp*1000d0 
!
! e2ic = (e2ab + e2ba)
! write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 

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
 write(LOUT,'(/1x,a,f16.8)') 'E2ind*      = ', e2ic*1000d0 

 SAPT%e2ind = e2ic

 deallocate(EVecA,OmA,EVecB,OmB)
 deallocate(AlphaA,AlphaB,WaBB,WbAA)

end subroutine e2ind_apsg

subroutine e2disp_unc(Flags,A,B,SAPT)
! calculate uncoupled and semi-coupled e2disp
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(EBlockData)               :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable   :: SBlockA(:),SBlockB(:)
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB,nblkA,nblkB
integer :: iunit,ival
integer :: i,j,ii,pq,rs,ip,iq,ir,is
integer :: iblk,ipos,kc,ik,ic
double precision,allocatable :: OmA0(:),OmB0(:),&
                                OmA1(:),OmB1(:)
double precision,allocatable :: EVecA1(:),EVecB1(:)
double precision,allocatable :: tmp01(:,:),y0y0(:,:),&
                                y1y0h(:,:),y1y0(:,:),&
                                y0y1(:,:)
double precision,allocatable :: work(:)
double precision           :: fact,dea,deb
double precision           :: e2du,e2d
double precision           :: e2ds,e2sp,e2ds1,e2ds2
double precision           :: e2dw12,e2dw12_sp,e2dsApp
double precision           :: inv_omega, inv_om12,tmp
double precision           :: Alpha, Beta
double precision,parameter :: BigE   = 1.D8 
double precision,parameter :: SmallE = 1.D-3

! Parameter(SmallE=1.D-3,BigE=1.D8)

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 dimOA  = A%num0+A%num1
 dimVA  = A%num1+A%num2
 dimOB  = B%num0+B%num1
 dimVB  = B%num1+B%num2
 nOVA   = dimOA*dimVA
 nOVB   = dimOB*dimVB
 nblkA  = 1+NBas-A%NAct
 nblkB  = 1+NBas-B%NAct

! read EigVals
 allocate(OmA0(A%NDimX),OmB0(B%NDimX))
 ! semi-coupled
 if(Flags%IFlag0==0) then
    allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX),&
             EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))  
 endif           

 ! semi-coupled
 if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
    ! CAS
    call readresp(EVecA1,OmA1,A%NDimX,'PROP_A1')
    call readresp(EVecB1,OmB1,B%NDimX,'PROP_B1')
 elseif(Flags%IFlag0==0.and.Flags%ICASSCF==0) then
    ! GVB
    call readEval(OmA1,A%NDimX,'PROP_A1')
    call readEval(OmB1,B%NDimX,'PROP_B1')
 endif

 allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
 call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
 call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')

 Alpha = 1.000d0
 Beta  = 1.000d0 

! tran4_gen
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

if(Flags%ISHF==1.and.SAPT%HFCheck) then

   write(LOUT,'(/,1x,a)') 'HARTREE-FOCK E2Disp REQUESTED'
   
   allocate(tmp01(A%NDimX,B%NDimX))
   
    tmp01=0
    do pq=1,A%NDimX
       ip = A%IndN(1,pq)
       iq = A%IndN(2,pq)
       ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
       read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
       do rs=1,B%NDimX
          ir = B%IndN(1,rs)
          is = B%IndN(2,rs)
   
          fact = &
                 work(is+(ir-B%num0-1)*dimOB)
   
          do i=1,A%NDimX
   
             tmp01(i,rs) = tmp01(i,rs) + &
                          fact * &
                          EVecA1(pq+(i-1)*A%NDimX)
          enddo
   
       enddo
    enddo
   
      e2du  = 0d0
      e2ds2 = 0d0
      e2ds1 = 0d0
      tmp   = 0d0
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
   
          if(OmA0(pq).gt.SmallE.and.OmB0(rs).gt.SmallE&
           .and.OmA0(pq).lt.BigE.and.OmB0(rs).lt.BigE) then
   
            tmp=0
            do kc=1,B%NDimX
               ik = B%IndN(1,kc)
               ic = B%IndN(2,kc)
               tmp = tmp + work(ic+(ik-B%num0-1)*dimOB)*EVecB1(kc+(rs-1)*B%NDimX)
            enddo
   
            e2du  = e2du  + work(is+(ir-B%num0-1)*dimOB)**2/(dea+deb)
            e2ds2 = e2ds2 + work(is+(ir-B%num0-1)*dimOB)**2*(Alpha*OmA1(pq)+Beta*OmB1(rs))/(dea+deb)**2
            e2ds1 = e2ds1 + (Beta*tmp+Alpha*tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(dea+deb) 
   
          endif
    
         enddo
      enddo
       e2sp = 4d0*(e2ds2-e2du)*1000 
       e2ds = (-4d0*e2du-16/sqrt(2d0)*e2ds1+4d0*e2ds2)*1000
   
       write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0
       write(LOUT,'(1x,a,f16.8)')'E2disp(sp)  = ', e2sp
       write(LOUT,'(1x,a,f16.8)')'E2disp(sc)  = ', e2ds
   
   deallocate(tmp01)
! end hf_check
endif

allocate(tmp01(A%NDimX,B%NDimX),y0y0(A%NDimX,B%NDimX))

if(Flags%IFlag0==0) then
 allocate(y1y0h(A%NDimX,B%NDimX),y1y0(A%NDimX,B%NDimX),&
          y0y1(A%NDimX,B%NDimX))
endif

! uncoupled E20disp

tmp01 = 0
if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
   y1y0h = 0
   y1y0  = 0
   ! unc & semi-cpld
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
  
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
    
         do i=1,A%NDimX
  
            !tmp01(i,rs) = tmp01(i,rs) + & 
            !             fact * &
            !             EVecA0(pq+(i-1)*A%NDimX)
 
 
            y1y0h(i,rs) = y1y0h(i,rs) + & 
                         fact * &
                         Alpha * &
                         EVecA1(pq+(i-1)*A%NDimX)
                            
         enddo

         associate(Y => Y01BlockA(pq))
            tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
         end associate

      enddo
   enddo

else 
! unc
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
  
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
    
         !do i=1,A%NDimX
  
         !   tmp01(i,rs) = tmp01(i,rs) + & 
         !                fact * &
         !                EVecA0(pq+(i-1)*A%NDimX)
         !                   
         !enddo
 
         associate(Y => Y01BlockA(pq))
            tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
         end associate
 
      enddo
   enddo

endif

y0y0 = 0
if(Flags%IFlag0==1.or.(Flags%IFlag0==0.and.Flags%ICASSCF==0)) then

   !do j=1,B%NDimX
   !   do i=1,A%NDimX
   !      do rs=1,B%NDimX
   !      ir = B%IndN(1,rs)
   !      is = B%IndN(2,rs)

   !      y0y0(i,j) = y0y0(i,j) + &
   !                   EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

   !      enddo
   !   enddo  
   !enddo

   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs)) 
        call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,y0y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo

elseif(Flags%IFlag0==0) then

   ! unc
   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs)) 
        call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,y0y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo

   y1y0 = 0
   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs)) 
        call dger(A%NDimX,Y%n,1d0,y1y0h(:,rs),1,Y%vec0,1,y1y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo
   call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB1,B%NDimX,0d0,y0y1,A%NDimX)

   !do j=1,B%NDimX
   !   do i=1,A%NDimX
   !      do rs=1,B%NDimX
   !      ir = B%IndN(1,rs)
   !      is = B%IndN(2,rs)

   !      !y0y0(i,j) = y0y0(i,j) + &
   !      !             EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

   !      !y1y0(i,j) = y1y0(i,j) + &
   !      !            EVecB0(rs+(j-1)*B%NDimX)*y1y0h(i,rs)

   !      !y0y1(i,j) = y0y1(i,j) + &
   !      !            Beta * &
   !      !            EVecB1(rs+(j-1)*B%NDimX)*tmp01(i,rs)

   !      enddo
   !   enddo  
   !enddo

endif

 e2du      = 0d0
 e2sp      = 0d0
 e2ds1     = 0d0
 e2ds2     = 0d0
 e2dw12    = 0d0
 e2dw12_sp = 0d0
 e2dsApp   = 0d0
 
 if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
 ! CAS
    do j=1,B%NDimX
       do i=1,A%NDimX
          ! remove this if later!!!!
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then
        !  if((OmA0(i)+OmA1(i)).gt.SmallE.and.(OmB0(j)+OmB1(j)).gt.SmallE &
        !   .and.(OmA1(i)+OmA0(i)).lt.BigE.and.(OmB1(j)+OmB0(j)).lt.BigE) then
 
          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          inv_om12  = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du  = e2du + y0y0(i,j)**2*inv_omega
          e2ds2 = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(y0y0(i,j)*inv_omega)**2
          e2ds1 = e2ds1 + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_omega

          e2dw12    = e2dw12 + y0y0(i,j)**2*inv_om12
          e2dw12_sp = e2dw12_sp + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_om12
  
          endif
       enddo
    enddo

 elseif(Flags%IFlag0==0.and.Flags%ICASSCF==0) then
 ! GVB
    do j=1,B%NDimX
       do i=1,A%NDimX
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then
        !  if((OmA0(i)+OmA1(i)).gt.SmallE.and.(OmB0(j)+OmB1(j)).gt.SmallE &
        !   .and.(OmA1(i)+OmA0(i)).lt.BigE.and.(OmB1(j)+OmB0(j)).lt.BigE) then
 
          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          inv_om12  = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du   = e2du + y0y0(i,j)**2*inv_omega
          e2ds2  = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(y0y0(i,j)*inv_omega)**2
          e2dw12 = e2dw12 + y0y0(i,j)**2*inv_om12
  
          endif
       enddo
    enddo

 elseif(Flags%IFlag0==1) then

    do j=1,B%NDimX
       do i=1,A%NDimX
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then

          inv_omega = 1d0/(OmA0(i)+OmB0(j))
   
          e2du = e2du + y0y0(i,j)**2*inv_omega
   
          endif
       enddo
    enddo

 endif

 call writeampl(y0y0,'PROP_AB0')

 SAPT%e2disp_unc = -16d0*e2du
 SAPT%e2disp_sp = -16*e2du+16*e2ds2
 SAPT%e2disp_sc = -16*e2du+16*e2ds2-32*e2ds1

 e2du = -16d0*e2du*1000d0
 e2sp = e2du + 16*e2ds2*1000 
 e2ds = e2du + (16*e2ds2-32*e2ds1)*1000

 e2dw12 = -16*e2dw12*1000
 e2dsApp = e2dw12 -32*e2dw12_sp*1000

 write(LOUT,'(/1x,a,f16.8)')'E2disp(unc) = ', e2du
 if(Flags%IFlag0==0) then
    write(LOUT,'(1x,a,f16.8)')'E2disp(sp) =  ',  e2sp
    write(LOUT,'(1x,a,f16.8)')'E2disp(sc) =  ',  e2ds
    write(LOUT,'(/1x,a,f16.8)')'E2disp(w12) = ', e2dw12
    write(LOUT,'(1x,a,f16.8)')'E2disp(sc2) = ',  e2dsApp
 endif

 close(iunit)
 deallocate(work)

 ! add keyword!
 !if(SAPT%DispExc)  call e2disp_unc_dexc(Flags,A,B,SAPT)

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

 if(Flags%IFlag0==0) then

    deallocate(y0y1,y1y0,y1y0h)
    deallocate(OmB1,EVecB1,OmA1,EVecA1)
 endif

 deallocate(y0y0,tmp01)
 deallocate(OmB0,OmA0)

end subroutine e2disp_unc

subroutine e2dispCAS(Flags,A,B,SAPT,NBasis)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(EBlockData)               :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable   :: SBlockA(:),SBlockB(:)
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)
integer,intent(in)           :: NBasis

integer          :: iunit
integer          :: i,j,pq,rs
integer          :: ip,iq,ir,is
integer          :: dimOA,dimVA, &
                    dimOB,dimVB,nOVA,nOVB
double precision :: fact,e2d
integer,allocatable          :: IGemA(:),IGemB(:)
double precision,allocatable :: OmA(:),OmB(:), &
                                EVecA(:),EVecB(:)
double precision,allocatable :: work(:)
double precision,allocatable :: tmp(:,:),tmp1(:,:)
! for Be ERPA:
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE   = 1.D8
double precision,parameter :: SmallE = 1.D-8

! print thresholds
if(SAPT%IPrint>1) then
   write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
   write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
   write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
endif

! set dimensions
dimOA = A%num0+A%num1
dimVA = A%num1+A%num2
dimOB = B%num0+B%num1
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

! fix IGem for A
allocate(IGemA(NBasis),IGemB(NBasis))
do i=1,A%INAct
   IGemA(i) = 1
enddo
do i=A%INAct+1,dimOA
   IGemA(i) = 2
enddo
do i=dimOA+1,NBasis
   IGemA(i) = 3
enddo
! fix IGem for B
do i=1,B%INAct
   IGemB(i) = 1
enddo
do i=B%INAct+1,dimOB
   IGemB(i) = 2
enddo
do i=dimOB+1,NBasis
   IGemB(i) = 3
enddo

allocate(OmA(A%NDimX),OmB(B%NDimX))
allocate(tmp(A%NDimX,B%NDimX),tmp1(A%NDimX,B%NDimX))

!new
allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
call convert_XY0_to_Y01(A,Y01BlockA,OmA,NBasis,'XY0_A')
call convert_XY0_to_Y01(B,Y01BlockB,OmB,NBasis,'XY0_B')

! check for negative eigenvalues
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative B!',i,OmB(i)
enddo

allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

tmp1 = 0
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
   do rs=1,B%NDimX
      ir = B%IndN(1,rs)
      is = B%IndN(2,rs)

      ! dispersion included in supermolecular CAS is computed
      if(IGemA(ip)==2.and.IGemA(iq)==2.and. &
         IGemB(ir)==2.and.IGemB(is)==2) then
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
      else
         fact = 0.d0
      endif

       associate(Y => Y01BlockA(pq))
         tmp1(Y%l1:Y%l2,rs) = tmp1(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
       end associate

   enddo
enddo

close(iunit)
deallocate(work)

! second multiplication
tmp = 0
do rs=1,B%NDimX
   associate(Y => Y01BlockB(rs))
     call dger(A%NDimX,Y%n,1d0,tmp1(:,rs),1,Y%vec0,1,tmp(:,Y%l1:Y%l2),A%NDimX)
   end associate
enddo
tmp1 = tmp

e2d = 0d0
do j=1,B%NDimX
   do i=1,A%NDimX
      if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
         .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

         e2d = e2d + tmp1(i,j)**2/(OmA(i)+OmB(j))

      endif
   enddo
enddo
SAPT%e2dispinCAS = -16d0*e2d
e2d  = -16d0*e2d*1000d0

! for exch-disp
call writeampl(tmp1,'PROP_AB0')

write(LOUT,'(/1x,a,f16.8)') 'E2disp(CAS) = ',e2d
print*, 'E2DISP(CAS) = ',e2d

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

deallocate(IGemB,IGemA)
deallocate(tmp1,tmp)
deallocate(OmB,OmA)

end subroutine e2dispCAS

subroutine e2disp_cpld(Flags,A,B,SAPT)
! CALCULATE ONLY COUPLED E2DISP
! does not work for GVB yet
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
integer :: NBas
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
double precision             :: fact,e2d
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: work(:)

!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE = 1.D8 
double precision,parameter :: SmallE = 1.D-3

 ! NOWE
 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

!check
 if(Flags%ICASSCF==0) then
   write(LOUT,*) 'ERROR! E2DISP_CPLD WORKS WITH CAS!'
   stop
 endif

! print thresholds
 if(SAPT%IPrint>1) then 
    write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
 endif

! set dimensions
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 allocate(EVecA(A%NDimX*A%NDimX))
 print*, 'alloc EVecA'
 allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))
 print*, 'alloc tmp1,tmp2'

 call readEvec(EVecA,A%NDimX,'PROP_A') 

 allocate(work(nOVB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

 tmp1 = 0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(is+(ir-B%num0-1)*dimOB)

       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + & 
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)
       enddo
      
    enddo
 enddo

 close(iunit)
 
 deallocate(EVecA)

 allocate(EVecB(B%NDimX*B%NDimX))
 call readEvec(EVecB,B%NDimX,'PROP_B') 

 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp2,A%NDimX)

 deallocate(EVecB,tmp1)

 allocate(OmA(A%NDimX),OmB(B%NDimX))
 call readEval2(OmA,A%NDimX,'PROP_A')
 call readEval2(OmB,B%NDimX,'PROP_B')

 do i=1,A%NDimX
    if(OmA(i)<0d0) write(LOUT,*) 'Negative A!',i,OmA(i)
 enddo
 do i=1,B%NDimX
    if(OmB(i)<0d0) write(LOUT,*) 'Negative B!',i,OmB(i)
 enddo

 e2d = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
       if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
          .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

          e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))

       endif
    enddo
 enddo
 SAPT%e2disp  = -16d0*e2d

 e2d  = -16d0*e2d*1000d0
 write(LOUT,'(/1x,a,f16.8)') 'E2disp      = ',e2d

 deallocate(OmB,OmA)
 deallocate(tmp2)
 deallocate(work)

end subroutine e2disp_cpld

subroutine e2disp(Flags,A,B,SAPT)
! calculate 2nd order dispersion energy
! in coupled and uncoupled approximations
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: kc,ik,ic
double precision,allocatable :: OmA(:), OmB(:), &
                                OmA0(:),OmB0(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),&
                                tmp01(:,:),tmp02(:,:)
double precision,allocatable :: work(:)
double precision :: e2d,fact,tmp
double precision :: e2du,dea,deb
double precision :: inv_omega
double precision :: Alpha, Beta
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
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
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

 Alpha = 1.000d0
 Beta  = 1.000d0 

! tran4_gen
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

! Hartree-Fock check
!if(Flags%ISHF==1.and.SAPT%HFCheck) then
!
!   allocate(tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX),&
!            tmp03(A%NDimX,B%NDimX))
!    tmp01=0
!    tmp03=0
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       !print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
!      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
!       do rs=1,B%NDimX
!          ir = B%IndN(1,rs)
!          is = B%IndN(2,rs)
!   
!          fact = &
!                 work(is+(ir-B%num0-1)*dimOB)
!   
!          do i=1,A%NDimX
!   
!             tmp01(i,rs) = tmp01(i,rs) + &
!                          fact * &
!                          EVecA1(pq+(i-1)*A%NDimX)
!   
!   !          tmp03(i,rs) = tmp03(i,rs) + &
!   !                       !fact * sqrt(2d0) * &
!   !                       fact * sqrt(2d0) * &
!   !                       EVecA0(pq+(i-1)*A%NDimX)
!   
!          enddo
!   
!       enddo
!    enddo
!   
!      e2du=0d0
!      e2ds2=0d0
!      e2ds1=0d0
!      tmp=0d0
!      do pq=1,A%NDimX
!         ip = A%IndN(1,pq)
!         iq = A%IndN(2,pq)
!         dea = A%OrbE(ip)-A%OrbE(iq)
!      !   print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
!         read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
!         do rs=1,B%NDimX
!            ir = B%IndN(1,rs)
!            is = B%IndN(2,rs)
!            deb = B%OrbE(ir)-B%OrbE(is) 
!      !      print*, is,ir,is+(ir-B%num0-1)*dimOB,nOVB
!   
!          if(OmA0(pq).gt.SmallE.and.OmB0(rs).gt.SmallE&
!           .and.OmA0(pq).lt.BigE.and.OmB0(rs).lt.BigE) then
!   
!            tmp=0
!            do kc=1,B%NDimX
!               ik = B%IndN(1,kc)
!               ic = B%IndN(2,kc)
!               tmp = tmp + work(ic+(ik-B%num0-1)*dimOB)*EVecB1(kc+(rs-1)*B%NDimX)
!            enddo
!   
!            e2du = e2du + work(is+(ir-B%num0-1)*dimOB)**2/(dea+deb)
!            e2ds2 = e2ds2 + work(is+(ir-B%num0-1)*dimOB)**2*(Alpha*OmA1(pq)+Beta*OmB1(rs))/(dea+deb)**2
!           ! e2ds1 = e2ds1 + (tmp02(pq,rs)+tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(OmA0(pq)+OmB0(rs))
!            e2ds1 = e2ds1 + (Beta*tmp+Alpha*tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(dea+deb) !(OmA0(pq)+OmB0(rs))
!   
!          endif
!    
!         enddo
!      enddo
!       write(LOUT,'()') 
!       e2ds = (-4d0*e2du-16/sqrt(2d0)*e2ds1+4d0*e2ds2)*1000
!       print*, 'DS1', -16/sqrt(2d0)*e2ds1*1000
!       write(LOUT,'(1x,a,f16.8)')'SPOLE:      = ', -4d0*e2du*1000d0 + 4d0*e2ds2*1000d0
!       write(LOUT,'(1x,a,f16.8)')'E2disp(sc)  = ', e2ds
!       write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0
!       write(LOUT,'()') 
!       e2ds  = 0
!       e2ds1 = 0
!       e2ds2 = 0 
!   
!   deallocate(tmp01,tmp02,tmp03)
!   ! end Hartree-Fock check
!endif

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
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(is+(ir-B%num0-1)*dimOB)

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
 ! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB0,B%NDimX,0d0,tmp02,A%NDimX)
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
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(is+(ir-B%num0-1)*dimOB)
  
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

 e2d = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
       if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
          .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

          e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))

       endif
    enddo
 enddo
 SAPT%e2disp  = -16d0*e2d

 e2d  = -16d0*e2d*1000d0

 write(LOUT,'(/1x,a,f16.8)')'E2disp      = ', e2d
 write(LOUT,'(/1x,a,f16.8)')'E2disp(unc) = ', e2du
! write(LOUT,*) 'E2DISP',e2d

 ! write amplitude to a file
 call writeampl(tmp2,'PROP_AB')

 ! calucate semicoupled and dexcitations
 if(SAPT%SemiCoupled) call e2disp_semi(Flags,A,B,SAPT)
 ! here: make it a keyword
 !if(SAPT%DispExc)     call e2disp_dexc(Flags,A,B,SAPT)

 close(iunit)
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
 
 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,OmA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp

subroutine e2disp_semi(Flags,A,B,SAPT)
type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: iunit,NBas
integer :: i,j,ip,iq,ipq,ir,is,irs
integer :: dimOA,dimVA,nOVA
integer :: dimOB,dimVB,nOVB

double precision :: fact,inv_omega
double precision :: e2ds1,e2ds2
double precision :: e2du,e2dsp,e2dsc

double precision,allocatable :: EVecA1(:),EVecB1(:)
double precision,allocatable :: OmA0(:), OmB0(:)
double precision,allocatable :: OmA1(:), OmB1(:)
double precision,allocatable :: tmp01(:,:),tmp02(:,:),&
                                sc10a(:,:),sc10b(:,:),sc01b(:,:)
double precision,allocatable :: work(:)

double precision,parameter :: BigE = 1.D8 
double precision,parameter :: SmallE = 1.D-3

! set dimensions
NBas = A%NBasis

dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

! read uncoupled
allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
allocate(OmA0(A%NDimX),OmB0(B%NDimX))

call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')

! read semicoupled
allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX))

call readresp(EVecA1,OmA1,A%NDimX,'PROP_A1')
!call readresp(EVecB1,OmB1,B%NDimX,'PROP_B1')

allocate(tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX))
allocate(sc10a(A%NDimX,B%NDimX),&
         sc01b(A%NDimX,B%NDimX),work(noVB))

open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

tmp01 = 0
tmp02 = 0
sc10a = 0

do ipq=1,A%NDimX
   ip = A%IndN(1,ipq)
   iq = A%IndN(2,ipq)
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
   do irs=1,B%NDimX
      ir = B%IndN(1,irs)
      is = B%IndN(2,irs)

      fact = (A%CICoef(iq)+A%CICoef(ip)) * &
             (B%CICoef(is)+B%CICoef(ir)) * &
              work(is+(ir-B%num0-1)*dimOB)

      do i=1,A%NDimX
         sc10a(i,irs) = sc10a(i,irs) + &
                       fact * &
                       EVecA1(ipq+(i-1)*A%NDimX)
      enddo

       associate(Y => Y01BlockA(ipq))
          tmp01(Y%l1:Y%l2,irs) = tmp01(Y%l1:Y%l2,irs) + fact * Y%vec0(1:Y%n)
       end associate

   enddo
enddo

! uncoupled
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB0,B%NDimX,0d0,tmp02,A%NDimX)
do irs=1,B%NDimX
   associate(Y => Y01BlockB(irs)) 
     call dger(A%NDimX,Y%n,1d0,tmp01(:,irs),1,Y%vec0,1,tmp02(:,Y%l1:Y%l2),A%NDimX)
   end associate
enddo

deallocate(EVecA1)

! semicoupled
allocate(EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))

call readEVec(EVecB1,B%NDimX,'PROP_B1')
call readEVal2(OmB1,B%NDimX,'PROP_B1')

call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB1,B%NDimX,0d0,sc01b,A%NDimX)

deallocate(tmp01)
allocate(sc10b(A%NDimX,B%NDimX))

sc10b = 0
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,sc10a,A%NDimX,EVecB0,B%NDimX,0d0,sc10b,A%NDimX)
do irs=1,B%NDimX
   associate(Y => Y01BlockB(irs)) 
     call dger(A%NDimX,Y%n,1d0,sc10a(:,irs),1,Y%vec0,1,sc10b(:,Y%l1:Y%l2),A%NDimX)
   end associate
enddo

! uncoupled and semicoupled
e2du  = 0d0
e2ds1 = 0d0
e2ds2 = 0d0
do j=1,B%NDimX
   do i=1,A%NDimX

      if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
         .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then

         inv_omega = 1d0/(OmA0(i)+OmB0(j))

         e2du = e2du + tmp02(i,j)**2*inv_omega
         e2ds2 = e2ds2 + (OmA1(i)+OmB1(j))*(tmp02(i,j)*inv_omega)**2
         e2ds1 = e2ds1 + tmp02(i,j)*(sc10b(i,j)+sc01b(i,j))*inv_omega

      endif
   enddo
enddo

SAPT%e2disp_sp  = -16*e2du+16*e2ds2
SAPT%e2disp_sc  = -16*e2du+16*e2ds2-32*e2ds1

e2du  = -16d0*e2du*1000d0
e2dsp = e2du + 16*e2ds2*1000 
e2dsc = e2du + (16*e2ds2-32*e2ds1)*1000

write(LOUT,'(1x,a,f16.8)')'E2disp(sp) =  ', e2dsp
write(LOUT,'(1x,a,f16.8)')'E2disp(sc) =  ', e2dsc

close(iunit)
deallocate(work)

deallocate(tmp02)
deallocate(sc01b,sc10b,sc10a)

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
deallocate(OmB0,OmA0)
deallocate(OmB1,EVecB1)
deallocate(OmA1)

end subroutine e2disp_semi

subroutine e2ind_dexc(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

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

! establish the type of correction: W_0j, W_1j, ...
! get info about monomer A
allocate(NSymAO(NBas))
call sym_inf('2RDMA',NumOSym,NSym,NumStSym,IStSy,NStSym,NSymAO,NBas)
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
write(lout,'(1x,a,i3)') 'offset = ', offset
write(lout,'(1x,a,f6.4)')  'OmA(1): ', OmA(1+offset)
write(lout,'(1x,a,f6.4)')  'OmA(2): ', OmA(2+offset)

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

print*, 'tmpA',norm2(tmpA)

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

subroutine e2disp_dexc(Flags,A,B,SAPT)
type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: iunit,NBas
integer :: nStates,iStA,nStSum,offset
integer :: i,j,ip,iq,ipq,ir,is,irs
integer :: dimOA,dimVA,nOVA
integer :: dimOB,dimVB,nOVB
integer :: NumOSym(8),NumStSym(16),IStSy(16), &
           NSym,NStSym
integer,allocatable :: NSymAO(:)

double precision :: fact,inv_omega,tmp

double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: tmpA(:),WbAA(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: Wij(:),work(:)

double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-4

write(lout,'(/1x,a)') 'W_ij corrections for E2ind and E2disp in excited states: (A*-B)'

! set dimensions
NBas = A%NBasis

dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

! establish the type of correction: W_ij , i<j
! get info about monomer A
allocate(NSymAO(NBas))
call sym_inf('2RDMA',NumOSym,NSym,NumStSym,IStSy,NStSym,NSymAO,NBas)
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

allocate(Wij(nStSum))

! read EigValA_B
allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),  &
         EVecB(B%NDimX*B%NDimX),OmB(B%NDimX))

call readresp(EVecA,OmA,A%NDimX,'PROP_A')
call readresp(EVecB,OmB,B%NDimX,'PROP_B')

! skip negative excitations
offset = 0
do i=1,A%NDimX

   if(OmA(i)<0d0) then 
      offset = offset + 1
   elseif(OmA(i)>0d0) then
      exit
   endif
   
enddo

write(lout,'(1x,a,3x,i3)') 'offset: ', offset
write(lout,'(1x,a,3x,i3)') 'nStSum: ', nStSum

write(lout,'(1x,a,f6.4)')  'OmA(1): ', OmA(1+offset)
write(lout,'(1x,a,f6.4)')  'OmA(2): ', OmA(2+offset)
write(lout,'(1x,a,f6.4/)') 'OmA(3): ', OmA(3+offset)

! induction part
! E2IND(A<--B)
allocate(tmpA(A%NDimX),WbAA(NBas,NBas))

call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)

tmpA = 0d0
!do ipq=1,A%NDimX
!   ip = A%IndN(1,ipq)
!   iq = A%IndN(2,ipq)
!
!   fact = (A%CICoef(iq)+A%CICoef(ip))
!
!   !do i=1,A%NDimX
!   do i=1+offset,nStSum+offset
!      tmpA(i) = tmpA(i) + fact*EVecA(ipq+(i-1)*A%NDimX)*WbAA(ip,iq)
!   enddo
!
!enddo
do i=1+offset,nStSum+offset
!do i=1,A%NDimX
   do ipq=1,A%NDimX
      ip = A%IndN(1,ipq)
      iq = A%IndN(2,ipq)
      fact = (A%CICoef(iq)+A%CICoef(ip))
   
      tmpA(i) = tmpA(i) + fact*EVecA(ipq+(i-1)*A%NDimX)*WbAA(ip,iq)

   enddo 
enddo

Wij = 0d0
tmp=0
do i=1+offset,nStSum+offset
!do i=1,A%NDimX
   if(OmA(i).lt.BigE.and.OmA(i).gt.SmallE) then
      Wij(i-offset) = Wij(i-offset) + tmpA(i)**2 / OmA(i)
      !tmp = tmp + tmpA(i)**2 / OmA(i)
   endif
enddo
Wij = 4.0d0*Wij

! print results
j = iStA
do i=1,nStSum
   write(lout,'(1x,a,2i1,a,f12.6)') 'Wind_',iStA-1,j,'  =',Wij(i)*1000d0
   write(lout,*) Wij(i)
   j = j + 1
enddo

deallocate(WbAA,tmpA)
write(lout,*) ''

! dispersion part
allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))

allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

tmp1 = 0
do ipq=1,A%NDimX
   ip = A%IndN(1,ipq)
   iq = A%IndN(2,ipq)
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
   do irs=1,B%NDimX
      ir = B%IndN(1,irs)
      is = B%IndN(2,irs)

      fact = (A%CICoef(iq)+A%CICoef(ip)) * &
             (B%CICoef(is)+B%CICoef(ir)) * &
              work(is+(ir-B%num0-1)*dimOB)

      do i=1+offset,nStSum+offset
         tmp1(i,irs) = tmp1(i,irs) + &
                      fact * &
                      EVecA(ipq+(i-1)*A%NDimX)
      enddo

   enddo
enddo
call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp2,A%NDimX)

Wij = 0d0
do j=1,B%NDimX
   do i=1+offset,nStSum+offset
      if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
         .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

         inv_omega = 1d0 / (-OmA(i) + OmB(j))

         Wij(i-offset) = Wij(i-offset) + tmp2(i,j)**2*inv_omega

      endif
   enddo
enddo
Wij = -16d0*Wij

! print results
j = iStA
do i=1,nStSum
   write(lout,'(1x,a,2i1,a,f12.6)') 'Wdisp_',iStA-1,j,' =',Wij(i)*1000
   j = j + 1
enddo

close(iunit)

deallocate(work)
deallocate(tmp2,tmp1)
deallocate(Wij)

deallocate(EVecB,OmB)
deallocate(EVecA,OmA)

end subroutine e2disp_dexc

subroutine e2disp_unc_dexc(Flags,A,B,SAPT)
type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

type(EBlockData)               :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable   :: SBlockA(:),SBlockB(:)
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: iunit,NBas
integer :: nStates,iStA,nStSum,offset
integer :: i,j,ip,iq,ipq,ir,is,irs
integer :: dimOA,dimVA,nOVA
integer :: dimOB,dimVB,nOVB
integer :: NumOSym(8),NumStSym(16),IStSy(16), &
           NSym,NStSym
integer,allocatable :: NSymAO(:)

double precision :: fact,inv_omega
double precision :: e2du

double precision,allocatable :: EVecA1(:),EVecB1(:)
double precision,allocatable :: OmA0(:), OmB0(:)
double precision,allocatable :: OmA1(:), OmB1(:)
double precision,allocatable :: tmp01(:,:),y0y0(:,:),&
                                y1y0h(:,:),y1y0(:,:),&
                                y0y1(:,:)
double precision,allocatable :: Wij_unc(:),work(:)

double precision,parameter :: BigE = 1.D8
double precision,parameter :: SmallE = 1.D-4

write(lout,'(/1x,a)') 'W_ij corrections for E2disp(unc,sc) in excited states: (A*-B)'

! set dimensions
NBas = A%NBasis

dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

! establish the type of correction: W_0j, W_1j, ...
! get info about monomer A
allocate(NSymAO(NBas))
call sym_inf('2RDMA',NumOSym,NSym,NumStSym,IStSy,NStSym,NSymAO,NBas)
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

 ! read EigVals
 allocate(OmA0(A%NDimX),OmB0(B%NDimX))
 allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX), &
          EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))

 ! semi-coupled
 if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
    ! CAS
    call readresp(EVecA1,OmA1,A%NDimX,'PROP_A1')
    call readresp(EVecB1,OmB1,B%NDimX,'PROP_B1')

 elseif(Flags%IFlag0==0.and.Flags%ICASSCF==0) then
    ! GVB
    write(lout,*) 'e2disp_dexc not available in GVB!' 
    return 

 endif

 allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
 call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
 call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')

 ! skip negative excitations
 offset = 0
 do i=1,A%NDimX
 
    if(OmA0(i)<0d0) then 
       offset = offset + 1
    elseif(OmA0(i)>0d0) then
       exit
    endif
    
 enddo
 
 print*, 'offset(unc) = ',offset
 write(lout,'(1x,a,f6.4)')  'OmA0(1): ', OmA0(1+offset)
 write(lout,'(1x,a,f6.4)')  'OmA0(2): ', OmA0(2+offset)
 write(lout,'(1x,a,f6.4/)') 'OmA0(3): ', OmA0(3+offset)

 allocate(work(nOVB))
 allocate(y0y0(A%NDimX,B%NDimX),tmp01(A%NdimX,B%NdimX))
 allocate(y1y0h(A%NDimX,B%NDimX),y1y0(A%NDimX,B%NDimX),&
          y0y1(A%NDimX,B%NDimX))

 open(newunit=iunit,file='TWOMOAB',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

 tmp01 = 0
 y1y0h = 0
 y1y0  = 0
 do ipq=1,A%NDimX
    ip = A%IndN(1,ipq)
    iq = A%IndN(2,ipq)
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do irs=1,B%NDimX
       ir = B%IndN(1,irs)
       is = B%IndN(2,irs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
              work(is+(ir-B%num0-1)*dimOB)
 
       do i=1,A%NDimX

       !   !tmp01(i,irs) = tmp01(i,irs) + &
       !   !             fact * &
       !   !             EVecA0(ipq+(i-1)*A%NDimX)

          y1y0h(i,irs) = y1y0h(i,irs) + &
                       fact * &
                       EVecA1(ipq+(i-1)*A%NDimX)

       enddo

       ! change here...
       associate(Y => Y01BlockA(ipq))
          tmp01(Y%l1:Y%l2,irs) = tmp01(Y%l1:Y%l2,irs) + fact * Y%vec0(1:Y%n)
       end associate

   enddo
 enddo

 y0y0 = 0
 ! unc
 do irs=1,B%NDimX
    associate(Y => Y01BlockB(irs)) 
      call dger(A%NDimX,Y%n,1d0,tmp01(:,irs),1,Y%vec0,1,y0y0(:,Y%l1:Y%l2),A%NDimX)
    end associate
 enddo
 
 y1y0 = 0
 do irs=1,B%NDimX
    associate(Y => Y01BlockB(irs)) 
      call dger(A%NDimX,Y%n,1d0,y1y0h(:,irs),1,Y%vec0,1,y1y0(:,Y%l1:Y%l2),A%NDimX)
    end associate
 enddo
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB1,B%NDimX,0d0,y0y1,A%NDimX)

 !do j=1,B%NDimX
 !   do i=1,A%NDimX
 !      do rs=1,B%NDimX
 !      ir = B%IndN(1,rs)
 !      is = B%IndN(2,rs)

 !      !y0y0(i,j) = y0y0(i,j) + &
 !      !             EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

 !      !y1y0(i,j) = y1y0(i,j) + &
 !      !            EVecB0(rs+(j-1)*B%NDimX)*y1y0h(i,rs)

 !      !y0y1(i,j) = y0y1(i,j) + &
 !      !            Beta * &
 !      !            EVecB1(rs+(j-1)*B%NDimX)*tmp01(i,rs)

 !      enddo
 !   enddo  
 !enddo

 allocate(Wij_unc(nStSum))

 ! CAS
 Wij_unc = 0d0
 do j=1,B%NDimX
    do i=1+offset,nStSum+offset
    !do i=1,A%NDimX

       if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
        .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then

       inv_omega = 1d0 / (-OmA0(i)+OmB0(j))

       Wij_unc(i-offset)  = Wij_unc(i-offset) + y0y0(i,j)**2*inv_omega
       !e2ds1 = e2ds1 + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_omega
       !e2ds1 = e2ds1 + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_omega

       endif
    enddo
 enddo
 Wij_unc = -16d0*Wij_unc

 ! print results
 j = iStA
 do i=1,nStSum
    write(lout,'(1x,a,2i1,a,f12.6)') 'W(unc)_',iStA-1,j,' =',Wij_unc(i)*1000
    j = j + 1
 enddo

 close(iunit)
 deallocate(work)

 deallocate(Wij_unc)

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

 deallocate(tmp01)
 deallocate(y1y0h,y1y0,y0y1)
 deallocate(OmB1,EVecB1,OmA1,EVecA1)
 deallocate(OmB0,OmA0)

end subroutine e2disp_unc_dexc

subroutine e2disp_pino(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ipq,irs
integer :: i1,i2,j1,j2,ii,jj,ij
integer :: coef,coef2,ADimEx,BDimEx
integer,allocatable :: AuxT(:,:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact,fpq,frs
double precision :: e2d,tmp
! testy
double precision :: val
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
double precision,allocatable :: AVecEx(:),BVecEx(:)

!double precision,parameter :: SmallE = 1.d-20
double precision,parameter :: SmallE = 1.d-6
!double precision,parameter :: SmallE = 1.d-1

 write(LOUT,'(/1x,a,e12.4)') 'SmallE in E2Disp PINO:',SmallE

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
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimFA = NBas
 dimOB = B%num0+B%num1
 dimFB = NBas 
 nOFA = dimOA*dimFA
 nOFB = dimOB*dimFB
 coef  = 1 
 coef2 = 1

 !if(allocated(A%PP)) then
 !   print*, size(A%PP),ADimEX
 !else
 !   print*, 'A%PP not allocated?!'
 !endif
 !if(allocated(A%PP)) then
 !   print*, size(B%PP),ADimEx
 !endif

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
      !!! test: Zero elements A
      ! do pq=1,ADimEx
      !     A%AP(i,pq) = 0d0 
      ! enddo

    endif
 enddo
 do i=1,coef*BDimEx
    if(B%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,B%PP(i)
      !!! test: Zero elements B
      ! do rs=1,BDimEx
      !     B%AP(i,rs) = 0d0 
      ! enddo
    endif
 enddo

 !! hm? zero He excitations
 !B%PP(4) = 0 
 !print*, 'He-exc energies'
 !do i=1,BDimEx
 !   write(LOUT,'(1x,i4,f16.8)') i,B%PP(i)
 !enddo

 !print*, 'DIMENSIONS ::',ij,NINte1,ADimEx

 !! test NORM
 !do i=1,ADimEx
 !val = 0
 !   do rs=1,ADimEx
 !      ir = AuxT(1,rs)
 !      is = AuxT(2,rs)

 !      frs = 2d0
 !      if(ir==is) frs=1d0
 !      val = val + frs*A%AP(i,rs)**2
 !   enddo
 !print*, 'val:',i,val
 !enddo
 !val = 0

!! ! test: remove diagonal elements...
! print*, 'Diagonal Elements in Mon%AP removed!'
! do i=1,ADimEx
!    do pq=1,ADimEx
!       ip = AuxT(1,pq)
!       iq = AuxT(2,pq)
!
!      if(ip==iq) then 
!         A%AP(i,pq) = 0d0
!      endif
!
!    enddo
! enddo
! do i=1,BDimEx
!    do rs=1,BDimEx
!       ir = AuxT(1,rs)
!       is = AuxT(2,rs)
!
!      if(ir==is) then 
!         B%AP(i,rs) = 0d0
!      endif
!
!   enddo
!enddo

 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)

 allocate(tmp1(coef*ADimEx,coef*BDimEx),tmp2(coef*ADimEx,coef*BDimEx))

 tmp1=0
 do pq=1,ADimEx
    ip = AuxT(1,pq)
    iq = AuxT(2,pq)

    fpq = 1d0
    if(ip==iq) fpq=0.5d0

    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,BDimEx
       ir = AuxT(1,rs)
       is = AuxT(2,rs)

       frs = 1d0
       if(ir==is) frs=0.5d0

       fact = fpq*frs &
             * (A%CICoef(iq)*B%CICoef(is)+A%CICoef(ip)*B%CICoef(ir) &
             + A%CICoef(iq)*B%CICoef(ir)+A%CICoef(ip)*B%CICoef(is)) & 
             * work(is+(ir-1)*dimOB)

       do i=1,coef*ADimEx

         ! if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*A%AP(i,pq)

         ! endif
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,coef*BDimEx

   ! if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
       do i=1,coef*ADimEx
   !       if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
             do rs=1,BDimEx

                tmp2(i,j) = tmp2(i,j) + &
                     B%AP(j,rs)*tmp1(i,rs)

             enddo
   !       endif
       enddo
   ! endif
 enddo

 e2d = 0d0
 do i=1,coef*ADimEx

    ! A: remove negative Evals
    !if(A%PP(i).gt.SmallE.and.A%PP(i).lt.1d20) then
    if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
       do j=1,coef*BDimEx

          ! B: remove negative Evals
          !if(B%PP(j).gt.SmallE.and.B%PP(j).lt.1d20) then
          if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
             e2d = e2d + tmp2(i,j)**2/(A%PP(i)+B%PP(j))
          endif
       enddo
    endif
 enddo

 !print*, 'TEST: ',-16d0*e2d*1000d0
 SAPT%e2disp = -16d0*e2d
 write(LOUT,'(/1x,a,f16.8)') 'E2disp      = ',-16*(e2d)*1000

 call writeampl(tmp2,'PROP_AB')

 deallocate(AuxT)
 deallocate(tmp2,tmp1)
 deallocate(work)

end subroutine e2disp_pino

subroutine e2disp_apsg(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ipq,irs
integer :: i1,i2,j1,j2,ii,jj,ij
integer :: coef,coef2,ADimEx,BDimEx
integer,allocatable :: AuxT(:,:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact
double precision :: e2d1,e2d2,e2d3,e2d4,e2d,tmp
double precision :: e2du,dea,deb
! testy
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
double precision,allocatable :: AVecEx(:),BVecEx(:)

double precision,parameter :: SmallE = 1.d-6
!double precision,parameter :: SmallE = 1.d-2
!double precision,parameter :: SmallE = 1.d-1

 write(LOUT,'(1x,a,e12.4)') 'SmallE in E2Disp PINO:',SmallE

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
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimFA = NBas
 dimOB = B%num0+B%num1
 dimFB = NBas 
 nOFA = dimOA*dimFA
 nOFB = dimOB*dimFB
 coef  = 1 
 coef2 = 1
 ! Use these coefficients for PINOVEC
 ! coef  = 2
 !coef2 = 4

 print*, A%num0,A%num1,A%num2
 print*, B%num0,B%num1,B%num2

 !print*, A%IGem
 !print*, ''
 !print*, B%IGem

! read EigValA_B
 allocate(EVecA(coef2*ADimEx**2),OmA(coef*ADimEx),&
          EVecB(coef2*BDimEx**2),OmB(coef*BDimEx))

 call readresp(EVecA,OmA,coef*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,coef*BDimEx,'PROP_B')

! do i=1,coef*BDimEx
!    print*, 'OmB',i,OmB(i)
! enddo

! print*, norm2(EVecA),norm2(EVecB)
! print*, 'OmA'
! do i=1,size(OmA)
!   ! if(OmA(i).gt.0d0)  write(*,*) i,OmA(i)
!    write(*,*) i,",",OmA(i)
! enddo
!
! print*, 'OmB'
! do i=1,size(OmB)
!    !  if(OmB(i).gt.0d0)  write(*,*) i,OmB(i)
!    write(*,*) i,",",OmB(i)
! enddo

 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)

 allocate(tmp1(coef*ADimEx,coef*BDimEx),tmp2(coef*ADimEx,coef*BDimEx))

! big test of a single loop
 allocate(AIndEx(2,ADimEx),BIndEx(2,BDimEx))
 allocate(AVecEx(coef2*ADimEx**2),BVecEx(coef2*BDimEx**2))

 AIndEx=A%IndNx
 BIndEx=B%IndNx

 ! prepare vectors
 AVecEx = 0
 BVecEx = 0
 do pq=1,ADimEx
    if(pq<=A%NDimX) then

        ip = A%IndN(1,pq)
        iq = A%IndN(2,pq)

        fact = A%CICoef(iq)+A%CICoef(ip)

        !if(fact<1d-8) then
        !   print*, ip,iq,A%CICoef(iq),A%CICoef(ip)
        !endif

        !fact = 1d0
        do i=1,coef*ADimEx
           AVecEx((i-1)*coef*ADimEx+pq) = fact * &
                                          EVecA((i-1)*coef*ADimEx+pq)
        enddo
 
    elseif(pq>A%NDimX) then

        ir = pq - A%NDimX
        !fact = 1d0
        fact = A%CICoef(ir)
        do i=1,coef*ADimEx
           AVecEx((i-1)*coef*ADimEx+pq) = fact * &
                                          EVecA((i-1)*coef*ADimEx+pq)
           !! for PINOVEC
           !AVecEx((i-1)*coef*ADimEx+pq) = fact * &
           !                               EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ir)
        enddo

    endif
 enddo

 do rs=1,BDimEx
     if(rs<=B%NDimX) then

        ir = B%IndN(1,rs)
        is = B%IndN(2,rs)

        fact = B%CICoef(ir)+B%CICoef(is)
        !fact = 1d0
        do i=1,coef*BDimEx
           BVecEx((i-1)*coef*BDimEx+rs) = fact * &
                                          EVecB((i-1)*coef*BDimEx+rs)
        enddo
 
    elseif(rs>B%NDimX) then

        ip = rs - B%NDimX
        fact = B%CICoef(ip)
        !fact = 1d0
        do i=1,coef*BDimEx
           BVecEx((i-1)*coef*BDimEx+rs) = fact * &
                                          EVecB((i-1)*coef*BDimEx+rs)
           !! for PINOVEC
           !BVecEx((i-1)*coef*BDimEx+rs) = fact * &
           !          EVecB((i-1)*coef*BDimEx+coef*B%NDimX+ip)
        enddo

    endif
 enddo

 do i=1,coef*ADimEx
    if(OmA(i)<0d0) then 
      ! !do pq=A%NDimX+1,ADimEx
      ! do pq=1,ADimEx
      !     AVecEx((pq-1)*coef*ADimEx+i) = 0d0 
      ! enddo
       write(LOUT,'(1x,"Monomer A: Negative EVal",i4,f16.8)') i,OmA(i)
    endif
 enddo
 do i=1,coef*BDimEx
    if(OmB(i)<0d0) then

       !!do rs=B%NDimX+1,BDimEx
       !do rs=1,BDimEx
       !    BVecEx((i-1)*coef*BDimEx+rs) = 0d0 
       !enddo
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,OmB(i)
    endif
 enddo


! good old stuff
 tmp1=0
 do pq=1,ADimEx
    ip = AIndEx(1,pq)
    iq = AIndEx(2,pq)
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,BDimEx
       ir = BIndEx(1,rs)
       is = BIndEx(2,rs)

       fact = work(is+(ir-1)*dimOB)

       do i=1,coef*ADimEx
          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*AVecEx((i-1)*coef*ADimEx+pq)

          endif
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,coef*BDimEx
    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
       do i=1,coef*ADimEx
          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
             do rs=1,BDimEx
                ir = BIndEx(1,rs)
                is = BIndEx(2,rs)

                tmp2(i,j) = tmp2(i,j) + &
                     BVecEx((j-1)*coef*BDimEx+rs)*tmp1(i,rs)

             enddo
          endif
       enddo
    endif
 enddo

 e2d1 = 0d0
 do i=1,coef*ADimEx

    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
       do j=1,coef*BDimEx

          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
             e2d1 = e2d1 + tmp2(i,j)**2/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo

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

!print*, 'ij:',ij,NINte1
!
!! try P matrices
! tmp1=0
! do pq=1,ADimEx
!    ip = AuxT(1,pq)
!    iq = AuxT(2,pq)
!
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,BDimEx
!       ir = AuxT(1,rs)
!       is = AuxT(2,rs)
!
!       fact = (A%CICoef(iq)*B%CICoef(is)+A%CICoef(ip)*B%CICoef(ir) &
!               + A%CICoef(iq)*B%CICoef(ir)+A%CICoef(ip)*B%CICoef(is))*work(is+(ir-1)*dimOB)
!       !fact = work(is+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!
!          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + & 
!                  !fact*AVecEx((i-1)*coef*ADimEx+pq)
!                  fact*A%AP(i,pq)
!
!          endif
!       enddo
!
!    enddo
! enddo
!
! tmp2=0
! do j=1,coef*BDimEx
!
!    if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!             do rs=1,BDimEx
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     !BVecEx((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!                     B%AP(j,rs)*tmp1(i,rs)
!
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
! e2d1 = 0d0
! do i=1,coef*ADimEx
!
!    if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!       do j=1,coef*BDimEx
!
!          if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
!             e2d1 = e2d1 + tmp2(i,j)**2/(A%PP(i)+B%PP(j))
!          endif
!       enddo
!    endif
! enddo
!
 print*, ''
 print*, 'TEST: ',-16d0*e2d1*1000d0
 SAPT%e2disp = -16d0*e2d1
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000

 call writeampl(tmp2,'PROP_AB')

 e2d1=0

 deallocate(BVecEx,AVecEx)
 deallocate(BIndEx,AIndEx)

 deallocate(AuxT)

!! old code in parts
!! PART 1: p>q,r>s
! tmp1=0
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
!
!       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
!              (B%CICoef(is)+B%CICoef(ir)) * work(is+(ir-1)*dimOB)
! !      print*, work(is+(ir-1)*dimOB)
!       !print*, fact
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + & 
!                  fact*EVecA((i-1)*coef*ADimEx+pq)
!
!          endif
!       enddo
!    enddo
! enddo
!
! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do rs=1,B%NDimX
!                ir = B%IndN(1,rs)
!                is = B%IndN(2,rs)
!                
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!             
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!!write(*,*) 'tmp1,tmp2:',norm2(tmp1),norm2(tmp2)
!
!! e2d1 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d1 = e2d1 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!             ! print*, OmA(i),OmB(j)
!!          endif
!!       enddo
!!    endif
!! enddo
!!
!!print*, ''
!!print*, 'PART1: ',-16d0*e2d1*1000d0
!!
!! !print*, 'test?',dimOB,B%NDimN
!
!! PART 2: p>q,r=s
! tmp1=0
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do ir=1,B%NDimN
!!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
!
!       fact = B%CICoef(ir)*(A%CICoef(iq)+A%CICoef(ip)) &
!            * work(ir+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,ir) = tmp1(i,ir) + &
!                  fact*EVecA((i-1)*coef*ADimEx+pq)
!
!          endif
!       enddo
!    enddo
! enddo
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do ir=1,B%NDimN
!          
!                tmp2(i,j) = tmp2(i,j) + &
!!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
!                     tmp1(i,ir)
!                
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! e2d2 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d2 = e2d2 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!          endif
!!       enddo
!!    endif
!! enddo
!! write(LOUT,*) 'PART2: ',-16d0*e2d2*1000d0
!
!! PART 3: p=q,r>s
! tmp1=0
! do ip=1,A%NDimN
!    !   ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
!    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!       !      !print*, is,ir,is+(ir-1)*dimOB,nOFB
!       fact = A%CICoef(ip) * &
!            (B%CICoef(is)+B%CICoef(ir)) * &
!            work(is+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + &
!                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)
!
!          endif
!       enddo
!    enddo
! enddo
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do rs=1,B%NDimX
!                ir = B%IndN(1,rs)
!                is = B%IndN(2,rs)
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! e2d3 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d3 = e2d3 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!          endif
!!       enddo
!!    endif
!! enddo
!! write(LOUT,*) 'PART3: ', -16d0*e2d3*1000d0
!
!! PART 4: p=q,r=s
! tmp1=0
! do ip=1,A%NDimN
!    ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
!    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
!    do ir=1,B%NDimN
!       !print*, is,ir,is+(ir-1)*dimOB,nOFB
!       fact = A%CICoef(ip)*B%CICoef(ir)* &
!              work(ir+(ir-1)*dimOB)
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,ir) = tmp1(i,ir) + &
!                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)
!
!          endif
!       enddo
!    enddo
! enddo
! 
!! print*, 'tmp1', norm2(tmp1(1:2*ADimEx,1:B%NDimN))
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do ir=1,B%NDimN
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
!                     tmp1(i,ir)
!                
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! print*, 'tmp2: ',norm2(tmp2)
!
! e2d4 = 0d0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             e2d4 = e2d4 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!            ! print*, tmp2(i,j)**2,OmA(i),OmB(j)
!          endif
!       enddo
!    endif
! enddo
! print*, 'PART4: ',-16*e2d4*1000

! SAPT%e2disp = -16d0*e2d
! write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000

 close(iunit)
 deallocate(work)
 deallocate(tmp2,tmp1)
 deallocate(OmB,EVecB,OmA,EVecA)

end subroutine e2disp_apsg



subroutine c6_dummy(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: kc, ik, ic
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: Adm(:,:,:),Bdm(:,:,:) 
double precision,allocatable :: intA(:,:),intB(:,:)
double precision :: fact_x,fact_y,fact_z,fact 
double precision :: c6val 
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis
 endif

! print thresholds
 if(SAPT%IPrint>1) then 
    write(LOUT,'(/,1x,a)') 'Thresholds in C6:'
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
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

! transform dipole momenets
 allocate(Adm(3,NBas,NBas),Bdm(3,NBas,NBas))

 call tran2MO(A%dipm(1,:,:),A%CMO,A%CMO,Adm(1,:,:),NBas)
 call tran2MO(A%dipm(2,:,:),A%CMO,A%CMO,Adm(2,:,:),NBas)
 call tran2MO(A%dipm(3,:,:),A%CMO,A%CMO,Adm(3,:,:),NBas)

 call tran2MO(B%dipm(1,:,:),B%CMO,B%CMO,Bdm(1,:,:),NBas)
 call tran2MO(B%dipm(2,:,:),B%CMO,B%CMO,Bdm(2,:,:),NBas)
 call tran2MO(B%dipm(3,:,:),B%CMO,B%CMO,Bdm(3,:,:),NBas)

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 allocate(intA(3,A%NDimX),intB(3,B%NDimX))

 !call print_sqmat(A%dipm(1,:,:),NBas) 
 !call print_sqmat(A%dipm(2,:,:),NBas) 
 !call print_sqmat(A%dipm(3,:,:),NBas) 
 print*, norm2(A%dipm)
 print*, norm2(B%dipm)

 intA=0d0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)

    fact_x = -1d0*(A%CICoef(iq)+A%CICoef(ip))*Adm(1,ip,iq)
    fact_y = -1d0*(A%CICoef(iq)+A%CICoef(ip))*Adm(2,ip,iq)
    fact_z = -1d0*(A%CICoef(iq)+A%CICoef(ip))*Adm(3,ip,iq)

    do i=1,A%NDimX

        intA(1,i) = intA(1,i) + &
                   fact_x * EVecA(pq+(i-1)*A%NDimX)
        intA(2,i) = intA(2,i) + &
                   fact_y * EVecA(pq+(i-1)*A%NDimX)
        intA(3,i) = intA(3,i) + &
                   fact_z * EVecA(pq+(i-1)*A%NDimX)

   enddo
 enddo

 print*, 'intA',norm2(intA(1,1:A%NDimX))

 intB=0d0
 do rs=1,B%NDimX
    ir = B%IndN(1,rs)
    is = B%IndN(2,rs)

    fact_x = -1d0*(B%CICoef(is)+B%CICoef(ir))*Bdm(1,ir,is)
    fact_y = -1d0*(B%CICoef(is)+B%CICoef(ir))*Bdm(2,ir,is)
    fact_z = -1d0*(B%CICoef(is)+B%CICoef(ir))*Bdm(3,ir,is)

    do j=1,B%NDimX

        intB(1,j) = intB(1,j) + &
                   fact_x * EVecB(rs+(j-1)*B%NDimX)
        intB(2,j) = intB(2,j) + &
                   fact_y * EVecB(rs+(j-1)*B%NDimX)
        intB(3,j) = intB(3,j) + &
                   fact_z * EVecB(rs+(j-1)*B%NDimX)

   enddo
 enddo

 print*, 'intB',norm2(intB(1,1:B%NDimX))

 c6val = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
       if(OmA(i).gt.SmallE.and.OmB(j).gt.SmallE&
          .and.OmA(i).lt.BigE.and.OmB(j).lt.BigE) then
       !if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
       !   .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

          fact = intA(1,i)*intB(1,j)+intA(2,i)*intB(2,j)-2d0*intA(3,i)*intB(3,j) 
          c6val = c6val + fact**2/(OmA(i)+OmB(j))

       endif
    enddo
 enddo
 write(LOUT,*) 16d0*c6val

deallocate(Bdm,Adm)
deallocate(intB,intA)
deallocate(OmB,EVecB,OmA,EVecA)

end subroutine c6_dummy

end module sapt_pol

