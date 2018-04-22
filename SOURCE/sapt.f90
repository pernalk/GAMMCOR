module sapt_ener
use types
use tran

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
!double precision,allocatable :: denA(:),Asq(:)

! set dimensions
 NBas = A%NBasis 

 allocate(PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),Ja(NBas,NBas))
 allocate(work(NBas,NBas))

 call get_den(NBas,A%CMO,2d0*A%Occ,PA)
 call get_den(NBas,B%CMO,2d0*B%Occ,PB)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call make_J1(NBas,PA,Ja)

! Tr[Pa.Va + Pb.Vb + Pb.Ja]
 work=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Vb,NBas,0d0,work,NBas)
 ea = trace(work,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Va,NBas,0d0,work,NBas)
 eb = trace(work,NBas) 
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Ja,NBas,0d0,work,NBas)
 ea = ea + trace(work,NBas)
 elst = ea + eb + SAPT%Vnn 

 write(LOUT,'(1x,a,f16.8)') 'V_nn        = ', SAPT%Vnn
 write(LOUT,'(1x,a,f16.8)') 'Eelst       = ', elst*1000d0 
 SAPT%elst = elst

 deallocate(work)
 deallocate(Ja,Vb,Va,PB,PA) 

end subroutine e1elst

subroutine e2ind(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
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
 write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000d0 

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
 write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0 

 e2ic = (e2ab + e2ba)
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 

 deallocate(WaBB,WbAA,AlphaB,AlphaA)
 deallocate(OmB,EVecB,Oma,EVecA)

end subroutine e2ind

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
! tran4_full
!allocate(work(NInte1))
!open(newunit=iunit,file='TWOMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NInte1)

! tran4_gen
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

if(Flags%ISHF==1) then
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
endif

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX))

! coupled - 0
 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + & 
                     (A%CICoef(iq)+A%CICoef(ip)) * &
                     (B%CICoef(is)+B%CICoef(ir)) * &
                     EVecA(pq+(i-1)*A%NDimX)* &
                     work(is+(ir-B%num0-1)*dimOB)
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


!! coupled - 1
!! works with tran4_full
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

 e2d = 0d0
 do i=1,A%NDimX
    do j=1,B%NDimX
       e2d = e2d + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
    enddo
 enddo
 SAPT%e2disp = -16d0*e2d
 e2d = -16d0*e2d*1000d0
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',e2d

!! coupled - TEST!
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
integer :: ii,jj
integer :: ADimEx,BDimEx
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact
double precision :: e2d1,e2d2,e2d3,e2d4,e2d,tmp
double precision :: e2du,dea,deb
double precision,parameter :: SmallE = 1.d-6 

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

! read EigValA_B
 allocate(EVecA(2*ADimEx*2*ADimEx),OmA(2*ADimEx),&
          EVecB(2*BDimEx*2*BDimEx),OmB(2*BDimEx))

 call readresp(EVecA,OmA,2*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,2*BDimEx,'PROP_B')

 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)

 allocate(tmp1(2*ADimEx,2*BDimEx),tmp2(2*ADimEx,2*BDimEx))

! PART 1: p>q,r>s
 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                        (A%CICoef(iq)+A%CICoef(ip)) * &
                        (B%CICoef(is)+B%CICoef(ir)) * &
                        EVecA((i-1)*2*ADimEx+pq)* &
                        work(is+(ir-1)*dimOB)

          endif
       enddo
    enddo
 enddo

 tmp2=0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do rs=1,B%NDimX
                ir = B%IndN(1,rs)
                is = B%IndN(2,rs)

                 tmp2(i,j) = tmp2(i,j) + &
                             EVecB((j-1)*2*BDimEx+rs)*tmp1(i,rs)

             enddo
          endif
       enddo   
    endif
 enddo

 e2d1 = 0d0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,2*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d1 = e2d1 + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 print*, -16d0*e2d1*1000d0

! PART 2: p>q,r=s
 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do ir=1,B%NDimN
!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB

       fact = sqrt(2d0)*work(ir+(ir-1)*dimOB)

       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,ir) = tmp1(i,ir) + &
                        (A%CICoef(iq)+A%CICoef(ip)) * &
                        fact*B%CICoef(ir) * &
                        EVecA((i-1)*2*ADimEx+pq)

          endif
       enddo
    enddo
 enddo
!
 tmp2=0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do ir=1,B%NDimN

                 tmp2(i,j) = tmp2(i,j) + &
                             EVecB((j-1)*2*BDimEx+2*B%NDimX+ir) * &
                             tmp1(i,ir)

             enddo
          endif
       enddo   
    endif
 enddo

 e2d2 = 0d0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,2*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d2 = e2d2 + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 print*, -16d0*e2d2*1000d0

! PART 3: p=q,r>s
 tmp1=0
 do ip=1,A%NDimN
!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
      ! fact = ... before the loop!!!
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                        sqrt(2d0)*A%CICoef(ip) * &
                        (B%CICoef(is)+B%CICoef(ir)) * &
                        EVecA((i-1)*2*ADimEx+2*A%NDimX+ip)* &
                        work(is+(ir-1)*dimOB)

          endif
       enddo
    enddo
 enddo
!
 tmp2=0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do rs=1,B%NDimX
                ir = B%IndN(1,rs)
                is = B%IndN(2,rs)

                 tmp2(i,j) = tmp2(i,j) + &
                             EVecB((j-1)*2*BDimEx+rs)*tmp1(i,rs)

             enddo
          endif
       enddo   
    endif
 enddo

 e2d3 = 0d0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,2*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d3 = e2d3 + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 print*, -16d0*e2d3*1000d0

! PART 4: p=q,r=s
 tmp1=0
 do ip=1,A%NDimN
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
    do ir=1,B%NDimN
       !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,ir) = tmp1(i,ir) + &
                        2d0*A%CICoef(ip) * &
                        B%CICoef(ir) * &
                        EVecA((i-1)*2*ADimEx+2*A%NDimX+ip)* &
                        work(ir+(ir-1)*dimOB)

          endif
       enddo
    enddo
 enddo
 
 print*, 'tmp1', norm2(tmp1(1:2*ADimEx,1:B%NDimN))
 
!
 tmp2=0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do ir=1,B%NDimN

                 tmp2(i,j) = tmp2(i,j) + &
                             EVecB((j-1)*2*BDimEx+2*B%NDimX+ir) * &
                             tmp1(i,ir)

             enddo
          endif
       enddo   
    endif
 enddo

 print*, 'tmp2: ',norm2(tmp2)

 e2d4 = 0d0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
       if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             e2d4 = e2d4 + tmp2(i,j)**2/(OmA(i)+OmB(j))
             !print*, tmp2(i,j)**2,OmA(i),OmB(j)
       endif
       enddo
    endif
 enddo
 print*, 'e2d4',-16*e2d4*1000



! SAPT%e2disp = -16d0*e2d
 print*, -16*(e2d2+e2d3+e2d4)*1000
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000


 close(iunit)
 deallocate(work)
 deallocate(tmp2,tmp1)
 deallocate(OmB,EVecB,OmA,EVecA)

end subroutine e2disp_apsg

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

subroutine calc_resp(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec(Mon%NDimX,Mon%NDimX), &
                               EVal(Mon%NDimX)
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha(Mon%NDimX,Mon%NDimX)
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is

 Alpha = 0
 do t=1,Mon%NDimX

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if

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

 enddo

end subroutine calc_resp

end module sapt_ener
