module exd_pino
use types
use tran
use timing
use exmisc
use sapt_utils

implicit none

contains

subroutine e2exdisp_pino(Flags,A,B,SAPT)

implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer :: NBas,NInte1
integer :: ADimEx,BDimEx
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: GdimOA,GdimOB
integer :: iunit
integer :: i,j,k,l,kl,ij,ii,jj,pq,rs
integer :: ip,iq,ir,is,ipq,irs
integer :: kc,ik,ic
logical :: approx
logical :: ipropab 
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
integer,allocatable :: posA(:,:),posB(:,:)
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:),&
                                work(:),workSq(:,:),&
                                sij(:,:),&
                                tpqrs(:,:)
double precision :: nelA,nelB
double precision :: fpq,frs
double precision :: fact,val,termZ,termY,termX
! 
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision,allocatable :: ints(:,:) 
double precision :: Tcpu,Twall
double precision,external  :: trace,FRDM2GVB
! test for Be
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 

! set dimensions
 NBas = A%NBasis 
 ADimEx = A%NDimX + A%NDimN
 BDimEx = B%NDimX + B%NDimN
 !dimOA = A%INAct+A%NAct
 GdimOA = A%num0+A%num1
 dimOA = GdimOA
 dimVA = A%num1+A%num2
! dimOB = B%INAct+B%NAct
 GdimOB = B%num0+B%num1
 dimOB = GdimOB
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          Sba(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas))
 allocate(tmp1(NBas,NBas),tmp2(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 !print*, A%Occ
 !print*, 'AAAA' 
 !print*, B%Occ

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

 allocate(RDM2Aval(GdimOA,GdimOA,GdimOA,GdimOA),&
          RDM2Bval(GdimOB,GdimOB,GdimOB,GdimOB))
 ! get RDM2 for FCI!!!

 !print*, 'dimOA,GdimOA',dimOA,GdimOA
 !print*, 'dimOB,GdimOB',dimOB,GdimOB

 do l=1,GdimOA
    do k=1,GdimOA 
       do j=1,GdimOA
          do i=1,GdimOA
             RDM2Aval(i,j,k,l) = FRDM2GVB(i,k,j,l,A%Occ,NBas)
          enddo
       enddo
    enddo
 enddo
 do l=1,GdimOB
    do k=1,GdimOB 
       do j=1,GdimOB
          do i=1,GdimOB
             RDM2Bval(i,j,k,l) = FRDM2GVB(i,k,j,l,B%Occ,NBas)
          enddo
       enddo
    enddo
 enddo

 !print*, 'RDM2A',norm2(RDM2Aval)
 !print*, 'RDM2B',norm2(RDM2Bval)

 allocate(posA(NBas,NBas),posB(NBas,NBas))
 ! old
 !posA = 0
 !do i=1,ADimEx
 !   if(i<=A%NDimX) then
 !      posA(A%IndNx(1,i),A%IndNx(2,i)) = A%IndX(i)
 !   elseif(i>A%NDimX) then
 !      posA(A%IndNx(1,i),A%IndNx(2,i)) = i 
 !   endif
 !enddo
 !posB = 0
 !do i=1,BDimEx
 !   if(i<=B%NDimX) then
 !      posB(B%IndNx(1,i),B%IndNx(2,i)) = B%IndX(i)
 !   elseif(i>B%NDimX) then
 !      posB(B%IndNx(1,i),B%IndNx(2,i)) = i
 !   endif
 !enddo

 ! new
 posA = 0
 do i=1,ADimEx
    posA(A%IndNT(1,i),A%IndNT(2,i)) = i
 enddo
 posB = 0
 do i=1,BDimEx
    posB(B%IndNT(1,i),B%IndNT(2,i)) = i
 enddo

 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 termZ = 2d0*SAPT%e2disp*termZ
! write(LOUT,'(1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3

 deallocate(tmp2,tmp1)
 
 allocate(tmp1(ADimEx,BDimEx),tmp2(ADimEx,BDimEx),&
          tmp3(ADimEx,BDimEx))
 allocate(sij(ADimEx,BDimEx))
!
! ..... TERM X!
 ! term A1
 ! transform J and K
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFOOABAB','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOABBA','AOTWOSORT')
 ! term A3
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAABB','AOTWOSORT')
! term A2
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          NBas,B%CMO,&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFFFABBB','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,A%CMO,&
          NBas,A%CMO,&
          NBas,B%CMO,&
          NBas,A%CMO,&
          'FFFFBAAA','AOTWOSORT')

 allocate(work(NBas**2),ints(NBas,NBas))

 sij=0
 tmp1=0

 call inter_A3_PINO(A,B,tmp1,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                    posB,posA,dimOB,dimOA,NBas)

 ! new
 ! TERMS P1-P2
 !!(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       ir = k
       iq = l

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                do j=1,dimOB
                   do i=1,NBas
                      ints(i,j) = work((j-1)*NBas+i)
                   enddo
                enddo

                fact = ints(ip,is)

                sij(ipq,irs) = sij(ipq,irs) + fact
 
             endif
          enddo
       enddo
 
    enddo
 enddo
 close(iunit)

 tmp2=0
 ! A1:P12
 call A1_Mod_1el_PINO(sij,tmp2,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'P12')

! A1: P34
! (FF|OO):(AB|AB)
 open(newunit=iunit,file='FFOOABAB',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 sij=0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       iq = k
       is = l

       do j=1,NBas
          do i=1,NBas
             ints(i,j) = work((j-1)*NBas+i)
          enddo
       enddo

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                fact = ints(ip,ir)

                sij(ipq,irs) = sij(ipq,irs) + fact
 
             endif
            enddo
         enddo
 
    enddo
 enddo
 close(iunit)

 ! A1:P34
 call A1_Mod_1el_PINO(sij,tmp2,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'P34')

 tmp1 = tmp1 + 2d0*tmp2

 ! A2(A)
 call inter_A2_PINO(A%DimEx,B%DimEx,tmp1, &
                Sab,nelA,Vabb,nelB,Vbab,'FFFFABBB',&
                A%CICoef,B%CICoef,A%IndNT,B%IndNT,posA,posB, &
                dimOA,dimOB,A%DimEx,B%DimEx,NBas,.false.)
 ! A2(B)
 call inter_A2_PINO(A%DimEx,B%DimEx,tmp1, &
                Sba,nelB,Vbaa,nelA,Vaba,'FFFFBAAA',&
                B%CICoef,A%CICoef,B%IndNT,A%IndNT,posB,posA, &
                dimOB,dimOA,B%DimEx,A%DimEx,NBas,.true.)

 call dgemm('N','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%AP,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
 call dgemm('N','T',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%AP,B%DimEx,0d0,tmp3,A%DimEx)

!!! term X

 sij = 0
 inquire(file='PROP_AB',EXIST=ipropab)
 if(ipropab) then
    ! read s_ij
    open(newunit=iunit,file='PROP_AB',form='UNFORMATTED',&
       access='SEQUENTIAL',status='OLD')
    read(iunit) sij
    close(iunit)

 else
 
    ! make s_ij
    !call make_sij_Y_PINO(sij,tmp1,A,B,nOVB,NBas)

 endif

 !print*, 'sij',norm2(sij)
 !print*, 'tmp3',norm2(tmp3)

 ! new
 termX = 0d0
 do j=1,B%DimEx
    do i=1,A%DimEx
        
       if(abs(A%PP(i)).gt.SmallE.and.abs(B%PP(j)).gt.SmallE&
          .and.abs(A%PP(i)).lt.BigE.and.abs(B%PP(j)).lt.BigE) then

          termX = termX + tmp3(i,j)*sij(i,j)/(A%PP(i)+B%PP(j))
          !termX = termX + tmp3(i,j)/(A%PP(i)+B%PP(j))

       endif
    enddo
 enddo
 termX = -4d0*termX

! write(LOUT,*) 'termX ', termX
! write(LOUT,'(/1x,a,f16.8)') 'term X      = ',  termX*1.0d3

! .........

! new
! TERM Y
 call make_tij_Y_PINO(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)

 ! term Y
 termY = 0d0
 do j=1,BDimEx
    do i=1,ADimEx
        
       if(abs(A%PP(i)).gt.SmallE.and.abs(B%PP(j)).gt.SmallE&
          .and.abs(A%PP(i)).lt.BigE.and.abs(B%PP(j)).lt.BigE) then

          termY = termY + sij(i,j)*tmp3(i,j)/(A%PP(i)+B%PP(j))

       endif
    enddo
 enddo


!! old
!! TERM Y
! call make_tij_Y_apsg(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
!
! ! term Y
! termY = 0d0
! do j=1,BDimEx
!    do i=1,ADimEx
!        
!       if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
!          .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then
!
!          termY = termY + sij(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))
!
!       endif
!    enddo
! enddo

!print*, 'Sij-norm',norm2(sij)
!print*, 'tij-norm',norm2(tmp3)
 termY = -8d0*(SAPT%elst-SAPT%Vnn)*termY

 if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3
 if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)') 'term Y      = ',  termY*1.0d3
 if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)') 'term X      = ',  termX*1.0d3

 write(LOUT,'(/1x,a,f16.8)') 'E2exch-disp = ', (termX+termY+termZ)*1.0d3
 SAPT%e2exdisp=(termX+termY+termZ)

 deallocate(work,ints)
 deallocate(sij)
 deallocate(RDM2Aval,RDM2Bval)
 deallocate(posB,posA)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vabb,Vb,Va,PB,PA,Sab,S)

end subroutine e2exdisp_pino

subroutine make_tij_Y_PINO(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas),Sba(NBas,NBas)
double precision,intent(inout) :: tmp1(A%DimEx,B%DimEx),&
                                  tmp2(A%DimEx,B%DimEx),&
                                  tmp3(A%DimEx,B%DimEx)

integer :: i,j,ii,jj,ir,is,irs,ip,iq,ipq
double precision :: fpq,frs
double precision :: fact,fact1,fact2

tmp1=0
do j=1,B%DimEx

   ir = B%IndNT(1,j)
   is = B%IndNT(2,j)
   frs = 1d0
   if(ir==is) frs = 0.5d0

   do i=1,A%DimEx

      ip = A%IndNT(1,i) 
      iq = A%IndNT(2,i)
      fpq = 1d0
      if(ip==iq) fpq = 0.5d0

      fact1 = (A%CICoef(iq)*B%CICoef(is)+A%CICoef(ip)*B%CICoef(ir))
      fact2 = (A%CICoef(ip)*B%CICoef(is)+A%CICoef(iq)*B%CICoef(ir))

      tmp1(i,j) = tmp1(i,j) + fpq*frs &
                * (Sab(iq,ir)*Sba(is,ip)*fact1 &
                + Sab(ip,ir)*Sba(is,iq)*fact2)

   enddo
enddo

!tmp3=0
!do j=1,B%DimEx 
!   do i=1,A%DimEx 
!      do jj=1,B%DimEX
!         do ii=1,A%DimEX
!            tmp3(ii,jj) = tmp3(ii,jj) + A%AP(ii,i)*tmp1(i,j)*B%AP(jj,j)
!         enddo
!      enddo
!   enddo  
!enddo  

call dgemm('N','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%AP,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
call dgemm('N','T',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%AP,B%DimEx,0d0,tmp3,A%DimEx)

end subroutine make_tij_Y_PINO

subroutine make_tij_Y_apsg(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas),Sba(NBas,NBas)
double precision,intent(inout) :: tmp1(A%DimEx,B%DimEx),&
                                  tmp2(A%DimEx,B%DimEx),&
                                  tmp3(A%DimEx,B%DimEx)

integer :: i,j,ir,is,irs,ip,iq,ipq
double precision :: fact,fpq,frs

tmp1=0
do j=1,B%DimEx

   ir = B%IndNx(1,j)
   is = B%IndNx(2,j)
   irs = posB(ir,is)
 
   frs = 1d0
   if(j<=B%NDimX) frs = B%Occ(ir)-B%Occ(is)

   do i=1,A%DimEx

      ip = A%IndNx(1,i)
      iq = A%IndNx(2,i)
      ipq = posA(ip,iq)

      fpq = 1d0
      if(i<=A%NDimX) fpq = A%Occ(ip)-A%Occ(iq)

      fact = fpq*frs 

      tmp1(ipq,irs) = tmp1(ipq,irs) + Sab(iq,ir)*Sba(is,ip)*fact

   enddo
enddo
!X_A.I.X_B
call dgemm('T','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%EigX,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
call dgemm('N','N',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%EigX,B%DimEx,0d0,tmp3,A%DimEx)
!print*, 'Xa.I.Xb:',norm2(tmp3)
!Y_A.I.Y_B
call dgemm('T','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%EigY,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
call dgemm('N','N',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%EigY,B%DimEx,1d0,tmp3,A%DimEx)
!print*, 'Ya.I.Yb:,test1',norm2(tmp3)

tmp1=0
do j=1,B%DimEx

   ir = B%IndNx(1,j)
   is = B%IndNx(2,j)
   irs = posB(ir,is)

   frs = 1d0
   if(j<=B%NDimX) frs = B%Occ(ir)-B%Occ(is)

   do i=1,A%DimEx

      ip = A%IndNx(1,i)
      iq = A%IndNx(2,i)
      ipq = posA(ip,iq)

      fpq = 1d0
      if(i<=A%NDimX) fpq = A%Occ(ip)-A%Occ(iq)

      fact = frs*fpq

      tmp1(ipq,irs) = tmp1(ipq,irs) - Sab(ip,ir)*Sba(is,iq)*fact

   enddo
enddo
!X_A.I.Y_B
call dgemm('T','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%EigX,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
call dgemm('N','N',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%EigY,B%DimEx,1d0,tmp3,A%DimEx)
!Y_A.I.X_B
call dgemm('T','N',A%DimEx,B%DimEx,A%DimEx,1d0,A%EigY,A%DimEx,tmp1,A%DimEx,0d0,tmp2,A%DimEx)
call dgemm('N','N',A%DimEx,B%DimEx,B%DimEx,1d0,tmp2,A%DimEx,B%EigX,B%DimEx,1d0,tmp3,A%DimEx)
!print*, 'test1:',norm2(tmp3)

end subroutine make_tij_Y_apsg

subroutine A1_Mod_1el_apsg(tmp,Vaab,Vbab,Sab,posA,posB,A,B,NBas,inter)
implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: NBas
character(*),intent(in) :: inter
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas),&
                               Vaab(NBas,NBas),Vbab(NBas,NBas)
double precision,intent(inout) :: tmp(A%DimEx,B%DimEx)

integer :: i,j,ip,iq,ipq,ir,is,irs
double precision :: fact,fpq,frs,ielA,ielB

 ielA = 1d0/(2d0*A%XELE)
 ielB = 1d0/(2d0*B%XELE)

select case(inter)
case('XX','xx')

   do j=1,B%DimEx
   
      ir = B%IndNx(1,j)
      is = B%IndNx(2,j)
      irs = posB(ir,is) 

      frs = B%Occ(ir)-B%Occ(is)
      if(j>B%NDimX) frs = 1d0

      do i=1,A%DimEx
   
         ip = A%IndNx(1,i)
         iq = A%IndNx(2,i)
         ipq = posA(ip,iq)
  
         fpq = A%Occ(ip)-A%Occ(iq)
         if(i>A%NDimX) fpq = 1d0 

         !fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
         !            (B%Occ(ir)-B%Occ(is))
         fact = -2d0*fpq*frs

         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(ip,is)*Sab(iq,ir) + fact*ielB*Vbab(iq,ir)*Sab(ip,is)
   
      enddo
   enddo

case('YY','yy')

   do j=1,B%DimEx
   
      ir = B%IndNx(1,j)
      is = B%IndNx(2,j)
      irs = posB(ir,is) 
   
      frs = B%Occ(ir)-B%Occ(is)
      if(j>B%NDimX) frs = 1d0 

      do i=1,A%DimEx
   
         ip = A%IndNx(1,i)
         iq = A%IndNx(2,i)
         ipq = posA(ip,iq)
   
         fpq = A%Occ(ip)-A%Occ(iq)
         if(i>A%NDimX) fpq = 1d0

         fact = -2d0*fpq*frs 
         !fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
         !            (B%Occ(ir)-B%Occ(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(iq,ir)*Sab(ip,is) + fact*ielB*Sab(iq,ir)*Vbab(ip,is)
   
      enddo
   enddo

case('XY','xy')

   do j=1,B%DimEx
   
      ir = B%IndNx(1,j)
      is = B%IndNx(2,j)
      irs = posB(ir,is) 

      frs = B%Occ(ir)-B%Occ(is)
      if(j>B%NDimX) frs = 1d0

      do i=1,A%DimEx
   
         ip = A%IndNx(1,i)
         iq = A%IndNx(2,i)
         ipq = posA(ip,iq)

         fpq = A%Occ(ip)-A%Occ(iq)
         if(i>A%NDimX) fpq = 1d0

         fact = 2d0*fpq*frs
         !
         !fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
         !           (B%Occ(ir)-B%Occ(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(ip,ir)*Sab(iq,is) + fact*ielB*Sab(ip,ir)*Vbab(iq,is)
   
      enddo
   enddo

case('YX','yx')

   do j=1,B%DimEx
   
      ir = B%IndNx(1,j)
      is = B%IndNx(2,j)
      irs = posB(ir,is) 

      frs = B%Occ(ir)-B%Occ(is)
      if(j>B%NDimX) frs = 1d0

      do i=1,A%DimEx
   
         ip = A%IndNx(1,i)
         iq = A%IndNx(2,i)
         ipq = posA(ip,iq)

         fpq = A%Occ(ip)-A%Occ(iq)
         if(i>A%NDimX) fpq = 1d0

         fact = 2d0*fpq*frs
         !
         !fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
         !            (B%Occ(ir)-B%Occ(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(iq,is)*Sab(ip,ir) + fact*ielB*Sab(iq,is)*Vbab(ip,ir)
   
      enddo
   enddo

end select

end subroutine A1_Mod_1el_apsg

subroutine A1_Mod_1el_PINO(twoel,tmp,Vaab,Vbab,Sab,posA,posB,A,B,NBas,inter)
implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: NBas
character(*),intent(in) :: inter
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas),&
                               Vaab(NBas,NBas),Vbab(NBas,NBas)
double precision,intent(in)    :: twoel(A%DimEx,B%DimEx)
double precision,intent(inout) :: tmp(A%DimEx,B%DimEx)

integer :: i,j,ip,iq,ipq,ir,is,irs
double precision :: fact,fpq,frs,ielA,ielB

 ielA = 1d0/(2d0*A%XELE)
 ielB = 1d0/(2d0*B%XELE)

select case(inter)
case('P12','p12')

   do j=1,B%DimEx
 
      ir = B%IndNT(1,j)
      is = B%IndNT(2,j)
      irs = posB(ir,is) 

      frs = 1d0
      if(ir==is) frs = 0.5d0

      do i=1,A%DimEx

         ip = A%IndNT(1,i)
         iq = A%IndNT(2,i)
         ipq = posA(ip,iq)
     
         fpq = 1d0
         if(ip==iq) fpq = 0.5d0
 
         tmp(ipq,irs) = tmp(ipq,irs) - fpq*frs*A%CICoef(iq)*B%CICoef(is) & 
                      * (twoel(ipq,irs)  & 
                      + ielA*Vaab(ip,is)*Sab(iq,ir) + ielB*Vbab(iq,ir)*Sab(ip,is)) &
                      - fpq*frs*A%CICoef(ip)*B%CICoef(ir) &
                      * (twoel(ipq,irs)  &
                      + ielA*Vaab(iq,ir)*Sab(ip,is) + ielB*Sab(iq,ir)*Vbab(ip,is))
 
      enddo
   enddo

case('P34','p34')

   do j=1,B%DimEx
 
      ir = B%IndNT(1,j)
      is = B%IndNT(2,j)
      irs = posB(ir,is) 

      frs = 1d0
      if(ir==is) frs = 0.5d0

      do i=1,A%DimEx
 
         ip = A%IndNT(1,i)
         iq = A%IndNT(2,i)
         ipq = posA(ip,iq)

         fpq = 1d0
         if(ip==iq) fpq = 0.5d0

         tmp(ipq,irs) = tmp(ipq,irs) &
                      - fpq*frs*A%CICoef(ip)*B%CICoef(is) &
                      * (twoel(ipq,irs) &
                      + ielA*Vaab(iq,is)*Sab(ip,ir) + ielB*Sab(iq,is)*Vbab(ip,ir)) &
                      - fpq*frs*A%CICoef(iq)*B%CICoef(ir) &
                      * (twoel(ipq,irs) &
                      + ielA*Vaab(ip,ir)*Sab(iq,is) + ielB*Sab(ip,ir)*Vbab(iq,is))
 
      enddo
   enddo

end select

end subroutine A1_Mod_1el_PINO

subroutine inter_A3_PINO(A,B,intP,Sab,AXELE,Vabb,BXELE,Vbaa,IntKFile,posB,posA,dimOB,dimOA,NBas)
implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dimOA,dimOB,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
character(*) :: IntKFile

double precision,intent(inout) :: intP(A%DimEx,B%DimEx)
double precision,intent(in) :: AXELE,BXELE
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iab,ipq
double precision :: fpq,fab,val
double precision :: Sxx(NBas,NBas),Sba(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

Sba = transpose(Sab)

allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ic = l 
      ia = k

      do ib=1,dimOB

         iab = posB(ia,ib)

         if(iab/=0) then 

            fab = 1d0
            if(ia==ib) fab = 0.5d0

            do i=1,A%DimEx
               
               ip = A%IndNT(1,i)
               iq = A%IndNT(2,i)
               ipq = posA(ip,iq)
        
               fpq = 1d0
               if(ip==iq) fpq = 0.5d0

               ! P1 
               val = 0
               do ir=1,dimOA
                  val = val + A%CICoef(ir)*Sab(ir,ib)*ints(ip,ir)
               enddo

               intP(ipq,iab) = intP(ipq,iab) + fpq*fab*val*B%CICoef(ic)*Sab(iq,ic)

               ! P2 
               val = 0
               do ir=1,dimOA
                  val = val + A%CICoef(ir)*Sab(ir,ib)*ints(iq,ir)
               enddo

               intP(ipq,iab) = intP(ipq,iab) + fpq*fab*val*B%CICoef(ic)*Sab(ip,ic)

            enddo
         endif

      enddo

     if(k<=dimOB) then

        ic = l
        ib = k

        do ia=1,NBas 

           iab = posB(ia,ib)

           if(iab/=0) then

              fab = 1d0
              if(ia==ib) fab = 0.5d0

              do i=1,A%DimEx
                 
                 ip = A%IndNT(1,i)
                 iq = A%IndNT(2,i)
                 ipq = posA(ip,iq)

                 fpq = 1d0
                 if(ip==iq) fpq = 0.5d0

                 ! P4
                 val = 0
                 do ir=1,dimOA
                    val = val + A%CICoef(ir)*ints(ip,ir)*Sab(ir,ia)
                 enddo

                 intP(ipq,iab) = intP(ipq,iab) + fpq*fab*B%CICoef(ic)*val*Sab(iq,ic)

                 ! P3
                 val = 0
                 do ir=1,dimOA
                    val = val + A%CICoef(ir)*ints(iq,ir)*Sab(ir,ia)
                 enddo

                 intP(ipq,iab) = intP(ipq,iab) + fpq*fab*B%CICoef(ic)*val*Sab(ip,ic)

              enddo

           endif
        enddo

     endif 

   enddo
enddo

close(iunit)

deallocate(ints,work)

intP = -2d0*intP

end subroutine inter_A3_PINO

subroutine inter_A2_PINO(dim1,dim2,intP,Smat,AXELE,Vabb,BXELE,Vbab,IntFile,&
                ACICoef,BCICoef,AIndNT,BIndNT,posA,posB, &
                dimOA,dimOB,ADimEx,BDimEx,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimEx,BDimEx,NBas
integer,intent(in) :: AIndNT(2,ADimEx),BIndNT(2,BDimEx),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntFile
double precision,intent(inout) :: intP(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,ACICoef(NBas),BCICoef(NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,is,it,iu,itu,ipq
double precision :: fpq,ftu,fact,val

double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

!(FF|FF):(AB|BB) or (BA|AA)
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

! one loop over integrals
ints = 0
kl = 0
do l=1,NBas
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*NBas)

      !call ints_modify(NBas,dimOB,ints,NBas,work,&
      !                 Smat(l,k)/AXELE,Vabb,Vbab(l,k)/BXELE,Sxx)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Smat,Sxx(k,l)/BXELE,Vbab)


      !do j=1,NBas
      !   do i=1,NBas
      !      ints(i,j) = work((j-1)*NBas+i)
      !   enddo
      !enddo

     ir = l
     ip = k 

     do iq=1,NBas

        ipq = posB(ip,iq)

        if(ipq/=0) then

           fpq = 1d0 
           if(ip==iq) fpq = 0.5d0 

           do i=1,ADimEx

              it = AIndNT(1,i)
              iu = AIndNT(2,i)
              itu = posA(it,iu)

              ftu = 1d0
              if(it==iu) ftu = 0.5d0
             
              fact = -2d0*ftu*fpq

              if(trans) then
                 ! P1
                 intP(ipq,itu) = intP(ipq,itu) + fact*ACICoef(iu)*BCICoef(ir)*ints(iu,iq)*Smat(it,ir)
                 ! P2
                 intP(ipq,itu) = intP(ipq,itu) + fact*ACICoef(it)*BCICoef(ir)*ints(it,iq)*Smat(iu,ir)
              else
                 ! P1
                 intP(itu,ipq) = intP(itu,ipq) + fact*ACICoef(iu)*BCICoef(ir)*ints(iu,iq)*Smat(it,ir)
                 ! P2
                 intP(itu,ipq) = intP(itu,ipq) + fact*ACICoef(it)*BCICoef(ir)*ints(it,iq)*Smat(iu,ir)
              endif

           enddo

        endif
     enddo

     ir = l 
     iq = k

     do ip=1,NBas

        ipq = posB(ip,iq)

        if(ipq/=0) then

           fpq = 1d0  
           if(ip==iq) fpq = 0.5d0

           do i=1,ADimEx

              it = AIndNT(1,i)
              iu = AIndNT(2,i)
              itu = posA(it,iu)

              ftu = 1d0
              if(it==iu) ftu = 0.5d0
             
              fact = -2d0*ftu*fpq
 
              if(trans) then
                 ! P3
                 intP(ipq,itu) = intP(ipq,itu) + fact*ACICoef(iu)*BCICoef(ir)*ints(iu,ip)*Smat(it,ir)
                 ! P4
                 intP(ipq,itu) = intP(ipq,itu) + fact*ACICoef(it)*BCICoef(ir)*ints(it,ip)*Smat(iu,ir)
              else
                 ! P3
                 intP(itu,ipq) = intP(itu,ipq) + fact*ACICoef(iu)*BCICoef(ir)*ints(iu,ip)*Smat(it,ir)
                 ! P4
                 intP(itu,ipq) = intP(itu,ipq) + fact*ACICoef(it)*BCICoef(ir)*ints(it,ip)*Smat(iu,ir)
              endif

           enddo

        endif
     enddo

   enddo
enddo

close(iunit)

deallocate(ints,work)

end subroutine inter_A2_PINO

end module exd_pino
