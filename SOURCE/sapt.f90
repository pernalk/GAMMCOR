module sapt_ener
use types
use tran
!use sapt_main

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

! read EigValA_B
 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 print*, 'E2ind-TEST'

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
 do i=1,A%NDimX
    do pq=1,A%NDimX
       ip = A%IndN(1,pq)
       iq = A%IndN(2,pq)
      ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
       read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
       do rs=1,B%NDimX
          ir = B%IndN(1,rs)
          is = B%IndN(2,rs)
          tmp1(i,rs) = tmp1(i,rs) + & 
                       (A%CICoef(iq)+A%CICoef(ip)) * &
                       (B%CICoef(is)+B%CICoef(ir)) * &
                        EVecA((i-1)*A%NDimX+pq)* &
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
                    EVecB((j-1)*B%NDimX+rs)*tmp1(i,rs)
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
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ii,jj
integer :: ADimEx,BDimEx
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: e2d,tmp
double precision :: e2du,dea,deb
double precision,parameter :: SmallE = 0d0!1.d-9

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
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

! read EigValA_B
 allocate(EVecA(2*ADimEx*2*ADimEx),OmA(2*ADimEx),&
          EVecB(2*BDimEx*2*BDimEx),OmB(2*BDimEx))

 call readresp(EVecA,OmA,2*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,2*BDimEx,'PROP_B')

 ! tran4_gen
 allocate(work(nOVB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

 allocate(tmp1(2*ADimEx,2*BDimEx),tmp2(2*ADimEx,2*BDimEx))

 tmp1=0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
    do pq=1,A%NDimX+A%NDimN
       if(pq.le.A%NDimX) then
          ip = A%IndN(1,pq)
          iq = A%IndN(2,pq)
       else
          ip = pq - A%NDimX
          iq = ip
       endif
       if(ip.ge.iq) then
         ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
          read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
          do rs=1,B%NDimX+B%NDimN
             if(rs.le.B%NDimX) then 
                ir = B%IndN(1,rs)
                is = B%IndN(2,rs)
             else
                ir = rs - B%NDimX
                is = ir
             endif
             if(ir.gt.is.and.ip.gt.iq) then   
                !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOVB

                   tmp1(i,rs) = tmp1(i,rs) + & 
                              (A%CICoef(iq)+A%CICoef(ip)) * &
                              (B%CICoef(is)+B%CICoef(ir)) * &
                              EVecA((i-1)*2*ADimEx+pq)* &
                              work(is+(ir-B%num0-1)*dimOB)
            
             elseif(ir.eq.is.and.ip.eq.iq) then
                   
                   tmp1(i,rs+B%NDimX) = tmp1(i,rs+B%NDimX) + &
                          2d0*(A%CICoef(iq)) * &
                              (B%CICoef(is)) * &
                              EVecA((i-1)*2*ADimEx+A%NDimX+pq)* &
                              work(is+(ir-B%num0-1)*dimOB)
!
!              elseif(ir.eq.is.and.ip.gt.iq) then
!
!                   tmp1(i,rs) = tmp1(i,rs) + & 
!                              (A%CICoef(iq)+B%CICoef(ip)) * &
!                              (B%CICoef(is)) * &
!                              EVecA((i-1)*2*ADimEx+pq)* &
!                              work(is+(ir-B%num0-1)*dimOB)
!
!              elseif(ir.gt.is.and.ip.eq.iq) then
!
!                   tmp1(i,rs) = tmp1(i,rs) + & 
!                              (A%CICoef(iq)) * &
!                              (B%CICoef(is)+B%CICoef(ir)) * &
!                              EVecA((i-1)*2*ADimEx+A%NDimX+pq)* &
!                              work(is+(ir-B%num0-1)*dimOB)

             endif
          enddo
       endif
    enddo
    endif
 enddo

 tmp2=0
 do j=1,2*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,2*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do rs=1,B%NDimX+B%NDimN
                if(rs.le.B%NDimX) then 
                   ir = B%IndN(1,rs)
                   is = B%IndN(2,rs)
                else
                   ir = rs - B%NDimX
                   is = ir
                endif

             if(ir.gt.is) then       
                 tmp2(i,j) = tmp2(i,j) + &
                             EVecB((j-1)*2*BDimEx+rs)*tmp1(i,rs)
!            elseif(ir.eq.is) then
!                  tmp2(i,j) = tmp2(i,j) + &
!                            EVecB((j-1)*2*BDimEx+B%NDimX+rs)*tmp1(i,rs+B%NDimX)
             endif

             enddo
          endif
       enddo   
    endif
 enddo

 e2d = 0d0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,2*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d = e2d + tmp2(i,j)**2d0/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 !SAPT%e2disp = -16d0*e2d
 e2d = -16d0*e2d*1000d0
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',e2d


! full-check
 e2d = 0d0
 do i=1,2*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
    do j=1,2*BDimEx
       if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then

       tmp=0d0
       do pq=1,A%NDimX+A%NDimN
          if(pq.le.A%NDimX) then
             ip = A%IndN(1,pq)
             iq = A%IndN(2,pq)
          else
             ip = pq - A%NDimX
             iq = ip
          endif
          if(ip.ge.iq) then          

             read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)

             do rs=1,B%NDimX+B%NDimN
                if(rs.le.B%NDimX) then
                   ir = B%IndN(1,rs)
                   is = B%IndN(2,rs)
                else
                   ir = rs - B%NDimX
                   is = ir
                endif

                if(ir.gt.is.and.ip.gt.iq) then
                   tmp = tmp + &
                      (A%CICoef(ip)+A%CICoef(iq)) * &
                      (B%CICoef(ir)+B%CICoef(is)) * &
                      EVecA((i-1)*2*ADimEx+pq) * &
                      EVecB((j-1)*2*BDimEx+rs) * & 
                      work(is+(ir-B%num0-1)*dimOB)

                elseif(ir.eq.is.and.ip.eq.iq) then
                       !print*, ir,is,rs,ip,iq,pq
                  ! tmp = tmp + &
                  !    2d0*(A%CICoef(ip)) * &
                  !    (B%CICoef(ir)) * &
                  !    EVecA((i-1)*2*ADimEx+pq) * &
                  !    EVecB((j-1)*2*BDimEx+rs) * & 
                  !    work(is+(ir-B%num0-1)*dimOB)

                endif

             enddo

          endif
       enddo

       e2d = e2d  + tmp**2d0/(OmA(i)+OmB(j))
       endif
    enddo
    endif
 enddo
 e2d = -16d0*e2d
 print*, 'e2d:',e2d

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



end module sapt_ener
