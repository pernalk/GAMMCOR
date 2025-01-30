#define SAPT_DSRS_DEBUG 11

module sapt_dRS
use types
use tran
use sapt_utils
use read_external
use sapt_inter

implicit none

contains

subroutine elst_dRS(A,B,SAPT)
!
! calculate electrostatic energy
! in degenerate RS (in NO representation)
!
! requires 1-RDMs and 1-TRDMs
!
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO,NBasis,dimOA,dimOB
integer :: iunit
integer :: iref,iref2,irefB,iref2B
integer :: iexcited,iexcited2
integer :: ip,ir,iq,is,j
double precision :: elab, elab2,elab3,tr
double precision :: elabB, elab2B,elab3B
double precision,allocatable :: work(:,:)
double precision,allocatable :: Btrdm(:,:)
double precision,allocatable :: Atrdm(:,:)

integer :: i
double precision,allocatable :: Vb(:,:)
double precision,allocatable :: Vbaa(:,:)
double precision,allocatable :: Vbaa2(:,:)
double precision,allocatable :: Va(:,:)
double precision,allocatable :: Vabb(:,:)
double precision,allocatable :: Vabb2(:,:)
double precision,allocatable :: workB(:,:)
double precision :: ea1,ea2,eas1,eb1,eb2,elstb1,elstb2
double precision :: sumA,sumB
double precision :: elst1,elst2,elstSAPT
double precision :: elst1exp,elst2exp,elstEXPsym1,elstEXPsym2,ElstEXPsymTrAng1,ELstEXPsymTrAng2
double precision :: A1B2A1B2,A2B1A2B1,A1B2A2B1
double precision :: pi
double precision :: gamma,gammadeg

pi = 4d0*atan(1d0)

! set dimensions
print*, 'A%IREF1',A%IREF1
print*, 'A%IREF2',A%IREF2
print*, 'B%IREF1',B%IREF1
print*, 'B%IREF2',B%IREF2
!print*, 'RDM2',A%rdm24
iref   =  A%IREF1 !1
iref2  =  A%IREF2 !4
irefB  =  B%IREF1
iref2B  = B%IREF2
iexcited = ((iref2-2)*(iref2-1))/2+iref
iexcited2 = ((iref2B-2)*(iref2B-1))/2+irefB

!print*, 'iexcited  = ',iexcited
!print*, 'iexcited2 = ',iexcited2

NAO    = SAPT%NAO
NBasis = A%NBasis
dimOA  = A%num0+A%num1
dimOB  = B%num0+B%num1

allocate(Vb(NAO,NAO),Va(NAO,NAO))
allocate(Vbaa(NBasis,NBasis),Vbaa2(NBasis,NBasis))
allocate(Vabb(NBasis,NBasis),Vabb2(NBasis,NBasis))

call get_one_mat('V',Vb,B%Monomer,NAO)
call get_one_mat('V',Va,A%Monomer,NAO)

call tran_AO2MO2(Vb,A%CAONO(iref,:,:), A%CAONO(iref,:,:), Vbaa, NAO,NBasis)
call tran_AO2MO2(Vb,A%CAONO(iref2,:,:),A%CAONO(iref2,:,:),Vbaa2,NAO,NBasis)

call tran_AO2MO2(Va,B%CAONO(irefB,:,:), B%CAONO(irefB,:,:), Vabb, NAO,NBasis)
call tran_AO2MO2(Va,B%CAONO(iref2B,:,:),B%CAONO(iref2B,:,:),Vabb2,NAO,NBasis)

! sum_p n_p v^B_pp
ea1 = 0d0
ea2 = 0d0
   do i=1,A%num0+A%num1
      ea1 = ea1 + A%rdm1(iref,i)*Vbaa(i,i)
      ea2 = ea2 + A%rdm1(iref2,i)*Vbaa2(i,i)
   enddo
   ea1 = 2d0*ea1
   ea2 = 2d0*ea2
eas1 = 0d0
if(irefB .eq. iref2B) then

   allocate(Atrdm(NBasis,NBasis))
   call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
                        Atrdm,NBasis)
   !!print*, 'Atrdm = ',norm2(Atrdm)

   do i=1,A%num0+A%num1
      do j=1,A%num0+A%num1
       eas1 = eas1 + Atrdm(i,j)*vbaa(i,j)
      enddo
   enddo
   eas1 = 2d0*eas1
endif

print*, 'ea1',ea1
print*, 'ea2',ea2
print*, 'eas1',eas1
eb1 = 0d0
eb2 = 0d0
do i=1,B%num0+B%num1
  eb1 = eb1 + B%rdm1(irefB,i)*Vabb(i,i)
  eb2 = eb2 + B%rdm1(iref2B,i)*Vabb2(i,i)
enddo
eb1 = 2d0*eb1
eb2 = 2d0*eb2
print*, 'eb1',eb1
print*, 'eb2',eb2

allocate(work(dimOA,dimOA))
! n_p * n_q * v_pq^pq
open(newunit=iunit,file='OOOOAABB',status='old',access='direct',&
     form='unformatted',recl=8*dimOA**2)
elab = 0
do ip=1,dimOB
   ! get all (AA| integrals for a given |BB) record
   read(iunit,rec=ip+(ip-1)*dimOB) work(1:dimOA,1:dimOA)
   do iq=1,dimOA
      elab = elab + A%rdm1(iref,iq)*B%rdm1(iref2B,ip)*work(iq,iq)
   enddo
enddo
close(iunit)

work=0d0
! n_p * n_q * v_pq^pq
open(newunit=iunit,file='OOOOAABB2',status='old',access='direct',&
     form='unformatted',recl=8*dimOA**2)
elabB = 0
do ip=1,dimOB
   ! get all (AA| integrals for a given |BB) record
   read(iunit,rec=ip+(ip-1)*dimOB) work(1:dimOA,1:dimOA)
   do iq=1,dimOA
      elabB = elabB + A%rdm1(iref2,iq)*B%rdm1(irefB,ip)*work(iq,iq)
   enddo
enddo
close(iunit)

! test
print*, 'elab = ',4d0*elab
print*, 'elabB = ',4d0*elabB

allocate(Btrdm(NBasis,NBasis))
if(irefB .ne. iref2B) then
   allocate(Atrdm(NBasis,NBasis))
   call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
                        Atrdm(:,:),NBasis)
   call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(iref2B,:,:),B%CMONO(iref2B,:,:), &
                        Btrdm(:,:),NBasis)
endif
if(irefB .eq. iref2B) then
Btrdm(:,:) = 0d0
   do i=1,NBasis
      Btrdm(i,i)=B%rdm1(irefB,i)
   enddo
endif

#if SAPT_DSRS_DEBUG > 10   
    print*, 'DTA in NO(iref) '
    do i=1,dimOA
       write(6,'(*(f13.8))') (Atrdm(i,j),j=1,dimOA)
    enddo
    print*, 'DTB in NO(iref2B) '
    do i=1,dimOB
       write(6,'(*(f13.8))') (Btrdm(i,j),j=1,dimOB)
    enddo
#endif

!call tran2MO(A%trdm(iref,:,:),CMONO(RefState,:,:),CMONO(RefState,:,:), &
!                        Atrdm(iref,:,:),NBasis)

open(newunit=iunit,file='OOOOAABB',status='old',access='direct',&
     form='unformatted',recl=8*dimOA**2)
elab2 = 0
elab3 = 0

do ip = 1,dimOB
   do ir = 1,dimOB
      ! get all (AA| integrals for a given |BB) record
      read(iunit,rec=ir+(ip-1)*dimOB) work(1:dimOA,1:dimOA)
      do iq = 1,dimOA
         do is = 1,dimOA
            elab3 = elab3 + Atrdm(iq,is)*Btrdm(ir,ip)*work(iq,is)
            if (ip == ir .and. iq == is) then
               elab2 = elab2 + A%rdm1(iref,iq)*B%rdm1(iref2B,ip)*work(iq,iq)
            endif
         enddo
      enddo
   enddo
enddo
close(iunit)

! test2
elab2 = 4d0*elab2
elab3 = 4d0*elab3
print*, 'elab2 = ',elab2
print *, 'elab3 = ',elab3

open(newunit=iunit,file='OOOOAABB2',status='old',access='direct',&
     form='unformatted',recl=8*dimOA**2)
elab2B = 0

do ip = 1,dimOB
   do ir = 1,dimOB
      ! get all (AA| integrals for a given |BB) record
      read(iunit,rec=ir+(ip-1)*dimOB) work(1:dimOA,1:dimOA)
      do iq = 1,dimOA
         do is = 1,dimOA
            if (ip == ir .and. iq == is) then
               elab2B = elab2B + A%rdm1(iref2,iq)*B%rdm1(irefB,ip)*work(iq,iq)
            endif
         enddo
      enddo
   enddo
enddo
close(iunit)

! test2
elab2B = 4d0*elab2B
print*, 'elab2B = ',elab2B

call print_en('V_nucA_elB*',eb2,.false.)
call print_en('V_nucB_elA', ea1,.false.)
call print_en('V_elA_elB*', elab2,.false.)
call print_en('V_nn',SAPT%Vnn,.false.)

call print_en('V_nucA_elB', eb1,.true.)
call print_en('V_nucB_elA*',ea2,.false.)
call print_en('V_elA*_elB', elab2B,.false.)
call print_en('V_nn',SAPT%Vnn,.false.)

A1B2A1B2=SAPT%Vnn+ea1+eb2+elab2
A2B1A2B1=SAPT%Vnn+eb1+ea2+elab2B
A1B2A2B1=elab3+eas1

write(LOUT,'(/1x,a,f16.8)') '<AB*|V|AB*>      = ', A1B2A1B2*1000d0
write(LOUT,'(1x,a,f16.8)')  '<A*B|V|A*B>      = ', A2B1A2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|V|A*B>      = ', A1B2A2B1*1000d0

gamma=0.5d0*atan2(2*A1B2A2B1,A1B2A1B2-A2B1A2B1)
gammadeg = gamma * 180d0 / pi
write(lout,'(/1x,a,f16.8)') 'gamma      = ', gamma
write(lout,'(1x,"Mixing angle (elst), rad", f12.8)') gamma
write(lout,'(1x,"Mixing angle (elst), deg", f9.3)')  gammadeg

print*,'4*gamma = ',4*gamma
print*,'COS(gamma)^2',COS(gamma)**2

!TRACK MINUS/PLUS SIGN  (NOT SURE ABOUT IT)
elst1exp=COS(gamma)**2*A1B2A1B2+SIN(gamma)**2*A2B1A2B1+2d0*SIN(gamma)*COS(gamma)*A1B2A2B1
elst2exp=SIN(gamma)**2*A1B2A1B2+COS(gamma)**2*A2B1A2B1-2d0*SIN(gamma)*COS(gamma)*A1B2A2B1
elstb1 = eb1+eb2+elab2B+elab3+SAPT%Vnn
elstb2 = eb1+eb2+elab2B-elab3+SAPT%Vnn

print *, 'EdRS(1,exp)= ' , elst1exp*1d3
print *, 'EdRS(2,exp)= ' , elst2exp*1d3
print*, "elstb1",elstb1*1d3
print*, "elstb2",elstb2*1d3

!print *, 'EdRS(1)= ' , (ea1+ea2+elab2+ABS(elab3)+SAPT%Vnn)*1d3
!print *, 'EdRS(2)= ' , (ea1+ea2+elab2-ABS(elab3)+SAPT%Vnn)*1d3

elst1 = ea1+ea2+elab2-elab3+SAPT%Vnn
elst2 = ea1+ea2+elab2+elab3+SAPT%Vnn
elstSAPT = ea1+ea2+elab2+SAPT%Vnn
!Na potrzeby testowe
SAPT%elst=elstSAPT
print *, 'EdRS(1)= ' , elst1*1d3
print *, 'EdRS(2)= ' , elst2*1d3
print *, 'elstSAPT= ', elstSAPT*1d3
print *, 'Vnn', SAPT%Vnn*1d3
!SAPT%elst  = elst1
SAPT%elst1 = elst1exp
SAPT%elst2 = elst2exp
SAPT%gamma = gamma
SAPT%A1B2VA1B2 = A1B2A1B2
SAPT%A2B1VA2B1 = A2B1A2B1
SAPT%A1B2VA2B1 = A1B2A2B1

deallocate(work)
deallocate(Btrdm,Atrdm)

end subroutine elst_dRS

subroutine e1exch_dSRS(A,B,SAPT)
!
! calculate exchange energy
! in degenerate SRS (in NO representation)
!
! describe IRef1 , IRef2 convention!
! (electronic and spatial degeneracy...)
!
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
integer :: NAO,NBasis
integer :: dimOA,dimOB
integer :: iunit
integer :: iref,iref2,irefB,iref2B
integer :: iexcited,iexcited2
integer :: ip,ir,iq,is,iu,it,iv,iw,i,j
double precision, allocatable :: trdm2A(:,:,:,:)
double precision, allocatable :: trdm2B(:,:,:,:)
double precision :: tvk(3),TNa(2),TNb(2),TNaNb,tvktest(3)
double precision, allocatable :: Atrdm(:,:)
double precision, allocatable :: Btrdm(:,:),Btrdm2(:,:)
double precision, allocatable :: Atrdmtest(:,:),Btrdmtest(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Sab2(:,:),Sab3(:,:)
double precision,allocatable :: Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:),Vbab(:,:)
double precision,allocatable :: Vabb2(:,:),Vbab2(:,:),Vaab2(:,:)
double precision, allocatable :: work(:,:)
double precision, allocatable :: work2A(:,:,:,:), work2B(:,:,:,:)
double precision, allocatable :: intA(:,:,:,:), intB(:,:,:,:)
double precision, allocatable :: ints(:,:)
double precision, allocatable :: tmpAB(:,:,:,:)
double precision, allocatable :: ints2(:)
double precision, allocatable :: workTEST(:,:)
double precision :: val, P1, exch1tot1, exch1tot2, VP,test
double precision :: VP2
double precision :: tvk2(3),TNa2(2),TNb2(2),TNaNb2,tvktest2(3)
double precision :: A1B2PA1B2,A2B1PA2B1,A1B2PA2B1
double precision :: A1B2VA1B2,A2B1VA2B1,A1B2VA2B1
double precision :: A1B2VPA1B2,A2B1VPA2B1,A1B2VPA2B1,A2B1VPA1B2
double precision :: P(2,2),V(2,2),VPnb(2,2)
double precision :: pi
double precision :: gamma
double precision :: gamma_exch1,gamma_exch2
double precision :: gamma_exch1deg,gamma_exch2deg
double precision :: tNaNbTEST

pi = 4d0*atan(1d0)

write(lout,'(/1x,a)') "Running e1exch_dSRS..."

print *,"exch w DRS ",SAPT%exchs21
print *,"SAPT%gamma",SAPT%gamma

! print*,"2TRDMA"
! do ip=1,10
!    do iq = 1,10
!       do ir = 1,10
!          do is =1,10
!             if(abs(A%trdm24(1,ip,iq,ir,is))>1d-2) then
!                print*,ip,iq,ir,is,A%trdm24(1,ip,iq,ir,is)
!             endif
!          enddo
!       enddo
!    enddo
! enddo


! print*,"2TRDMB"
! do ip=1,10
!    do iq = 1,10
!       do ir = 1,10
!          do is =1,10
!             if(abs(A%trdm24(1,ip,iq,ir,is))>1d-2) then
!                print*,ip,iq,ir,is,B%trdm24(1,ip,iq,ir,is)
!             endif
!          enddo
!       enddo
!    enddo
! enddo


! print*,"Czy 2TRDMB+2TRDMA=0"
! do ip=1,10
!    do iq = 1,10
!       do ir = 1,10
!          do is =1,10
!             if(abs(A%trdm24(1,ip,iq,ir,is)+B%trdm24(1,ip,iq,ir,is))>1d-2) then
!                print*,ip,iq,ir,is,B%trdm24(1,ip,iq,ir,is)+A%trdm24(1,ip,iq,ir,is)
!             endif
!          enddo
!       enddo
!    enddo
! enddo

!CZYTAMY KOD NA NOWO

iref   =  A%IREF1 !1
iref2  =  A%IREF2 !4
irefB  =  B%IREF1
iref2B  = B%IREF2
iexcited = ((iref2-2)*(iref2-1))/2+iref
iexcited2 = ((iref2B-2)*(iref2B-1))/2+irefB
tvk = 0d0

NAO    = SAPT%NAO
NBasis = A%NBasis
dimOA  = A%num0+A%num1
dimOB  = B%num0+B%num1

allocate(Atrdm(NBasis,NBasis))
allocate(Btrdm(NBasis,NBasis),Btrdm2(NBasis,NBasis))

allocate(S(NAO,NAO),Va(NAO,NAO),Vb(NAO,NAO))
allocate(Sab(NBasis,NBasis),Sab2(NBasis,NBasis),Sab3(NBasis,NBasis))
allocate(Vabb(NBasis,NBasis),Vbaa(NBasis,NBasis),&
         Vaab(NBasis,NBasis),Vbba(NBasis,NBasis),&
         Vbab(NBasis,NBasis),Vbab2(NBasis,NBasis),Vabb2(NBasis,NBasis),Vaab2(NBasis,NBasis))

! get V_ne in atomic orbs
call get_one_mat('V',Va,A%Monomer,NAO)
call get_one_mat('V',Vb,B%Monomer,NAO)

!call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
! not sure all are needed
call tran_AO2MO2(Va,B%CAONO(iref2B,:,:),B%CAONO(iref2B,:,:),Vabb,NAO,NBasis)
call tran_AO2MO2(Va,B%CAONO(irefB,:,:), B%CAONO(irefB,:,:), Vabb2,NAO,NBasis)
call tran_AO2MO2(Vb,A%CAONO(iref,:,:),  A%CAONO(iref,:,:),  Vbaa,NAO,NBasis)
call tran_AO2MO2(Va,A%CAONO(iref,:,:),  B%CAONO(iref2B,:,:),Vaab,NAO,NBasis)
call tran_AO2MO2(Va,A%CAONO(iref,:,:),  B%CAONO(irefB,:,:), Vaab2,NAO,NBasis)
call tran_AO2MO2(Vb,B%CAONO(iref2B,:,:),A%CAONO(iref,:,:),  Vbba,NAO,NBasis)
call tran_AO2MO2(Vb,A%CAONO(iref,:,:),  B%CAONO(iref2B,:,:),Vbab,NAO,NBasis)
call tran_AO2MO2(Vb,A%CAONO(iref,:,:),  B%CAONO(irefB,:,:), Vbab2,NAO,NBasis)

! get overlap S matrix in AO and transform to NOs
call get_one_mat('S',S,A%Monomer,NAO)

!call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
call tran_AO2MO2(S,A%CAONO(iref,:,:), B%CAONO(iref2B,:,:),Sab,NAO,NBasis)
call tran_AO2MO2(S,A%CAONO(iref2,:,:),B%CAONO(irefB,:,:),Sab3,NAO,NBasis)

call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
             Atrdm(:,:),NBasis)
if(irefB .ne. iref2B) then
   call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(iref2B,:,:),B%CMONO(iref2B,:,:), &
             Btrdm(:,:),NBasis)
endif

if(irefB .eq. iref2B) then
   Btrdm = 0d0
   do i=1,NBasis
      Btrdm(i,i) = B%rdm1(irefB,i)
   enddo
   Btrdm2 = Btrdm
endif

allocate(trdm2A(dimOA,dimOA,dimOA,dimOA))
allocate(trdm2B(dimOB,dimOB,dimOB,dimOB))
! trdm2A(:,:,:,:)=A%trdm24(iref2,:,:,:,:)! TO TYLKO DO TESTÓW TRZEBA ZMIENIĆ PÓŹNIEJ
! trdm2B(:,:,:,:)=B%trdm24(iref2B,:,:,:,:)

trdm2A(:,:,:,:)=A%trdm24(:,:,:,:)
if(irefB .ne. iref2B) then
   trdm2B(:,:,:,:)=B%trdm24(:,:,:,:)
endif
if(irefB .eq. iref2B) then
! ZMIANA 11.30.23
   trdm2B(:,:,:,:)=0d0
   trdm2B(1:dimOB,1:dimOB,1:dimOB,1:dimOB) = B%rdm24(irefB,:,:,:,:)
! KONIEC ZMIANY
!   trdm2B(:,:,:,:)=B%rdm24(irefB,:,:,:,:)
!USUNIETO W RAMACH ZMIANY
endif

! "2" denotes basis of B ground-state NOs
call tran_AO2MO2(S,A%CAONO(iref,:,:),B%CAONO(irefB,:,:),Sab2,NAO,NBasis)
deallocate(Vb,Va,S)
if(irefB .ne. iref2B) then
   call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(irefB,:,:),B%CMONO(irefB,:,:), &
             Btrdm2(:,:),NBasis)
endif

#if SAPT_DSRS_DEBUG > 10
print*, 'Btrdm2 MO',norm2(B%trdm1(iexcited2,:,:))
do j=1,dimOB
   write(6,'(*(f12.6))') (B%trdm1(iexcited2,i,j),i=1,dimOB)
enddo
print*, 'Btrdm2 NO',norm2(Btrdm2)
do j=1,dimOB
   write(6,'(*(f12.6))') (Btrdm2(i,j),i=1,dimOB)
enddo
#endif

A1B2PA1B2 = 0d0
do j=1,dimOB
   do i=1,dimOA
      !   nnS2 = nnS2 + A%Occ(i)*B%Occ(j)*Sab(i,j)**2
      A1B2PA1B2 = A1B2PA1B2 + A%rdm1(iref,i)*B%rdm1(iref2B,j)*Sab(i,j)**2
   enddo
enddo

A1B2PA1B2=-2d0*A1B2PA1B2

A2B1PA2B1 = 0d0
do j=1,dimOB
   do i=1,dimOA
         !   nnS2 = nnS2 + A%Occ(i)*B%Occ(j)*Sab(i,j)**2
         A2B1PA2B1 = A2B1PA2B1 + A%rdm1(iref2,i)*B%rdm1(irefB,j)*Sab3(i,j)**2
   enddo
enddo

A2B1PA2B1=-2d0*A2B1PA2B1
   
! !Trik do testowania
! Atrdm=0d0
! Btrdm=0d0
! Btrdm2=0d0
! trdm2A=0d0
! trdm2B=0d0
! do ip=1,dimOA
!    do iq=1,dimOB
!       if(ip == iq) then
!       Atrdm(ip,iq) = A%rdm1(1,ip)
!       Btrdm2(ip,iq) = B%rdm1(1,ip)
!       endif
!    enddo
! enddo
! trdm2A(1:dimOA,1:dimOA,1:dimOA,1:dimOA)=A%rdm24(1,:,:,:,:)
! trdm2B(1:dimOA,1:dimOA,1:dimOA,1:dimOA)=B%rdm24(1,:,:,:,:)                            
! !!!!! KONIEC ZMIAN DO TESTOWANIA




!tvk(1)=- 4sum(trdm(A->A*)_{pq}trdm(B*->B)_{rs}VB_{ps}S_{qr})
!P_1=-2sum(trdm(A->A*)_{pq}trdm(B*->B)_{rs}S_{ps}S_{qr})
P1   = 0d0
tvk  = 0d0
tvk2 = 0d0
do iq=1,dimOA
   do ip=1,dimOA
      do ir=1,dimOB
         do is=1,dimOB
            !     tvk(1) = tvk(1) + A%Occ(ip)*B%Occ(iq)*Vaab(ip,iq)*Sab(ip,iq)
            P1 = P1 + Atrdm(ip,iq)*Btrdm(ir,is)*Sab(ip,ir)*Sab(iq,is) !poprawione
            tvk(2) = tvk(2) + Atrdm(ip,iq)*Btrdm(ir,is)*Vbab(ip,ir)*Sab(iq,is) !poprawione
            tvk(1) = tvk(1) + Atrdm(ip,iq)*Btrdm(ir,is)*Vaab(iq,is)*Sab(ip,ir)
            tvk2(2) = tvk2(2) + Atrdm(iq,ip)*Btrdm(is,ir)*Vbab(ip,ir)*Sab(iq,is) 
            tvk2(1) = tvk2(1) + Atrdm(iq,ip)*Btrdm(is,ir)*Vaab(iq,is)*Sab(ip,ir)
            tvktest(2) = tvktest(2) + Atrdm(ip,iq)*Btrdm2(ir,is)*Vbab2(ip,ir)*Sab2(iq,is)
            tvktest(1) = tvktest(1) + Atrdm(ip,iq)*Btrdm2(ir,is)*Vaab2(iq,is)*Sab2(ip,ir)
         enddo
      enddo
   enddo
enddo
tvk(1) = -2d0*tvk(1)
tvk(2) = -2d0*tvk(2)
tvk2(1) = -2d0*tvk2(1)
tvk2(2) = -2d0*tvk2(2)
tvktest(1) = -2d0*tvktest(1)
tvktest(2) = -2d0*tvktest(2)
P1= -2d0*P1
!print*, 'tvk(1) z dSRS',tvk(1)*1000

! print*,"A occ", A%occ(:)
! print*, "A rdm1",A%rdm1(1,:)

!WORK IN PROGRESS
allocate(work(dimOA,dimOB))
! tvk(2)=- 2sum(trdm(A->A*)_{pq}trdm(B->B*)_{rs}(pr|qs)
open(newunit=iunit,file='OOOOABBA',status='old',access='direct',&
     form='unformatted',recl=8*dimOA*dimOB)
do is=1,dimOB
   do iq=1,dimOA
   ! get all (AB| integrals for a given |BA) record
      read(iunit,rec=is+(iq-1)*dimOB) work(1:dimOA,1:dimOB)
      do ir = 1,dimOB!
         do ip = 1,dimOA
            tvk(3) = tvk(3) + Atrdm(ip,iq)*Btrdm(ir,is)*work(ip,ir) !poprawione
         enddo
      enddo
   enddo
enddo
close(iunit)
tvk(3) = -2d0*tvk(3)
tvk2(3)= tvk(3)
!print*, 'tvk(2) z dSRS',tvk(2)*1000
! +
!allocate(work2(dimOA,dimOB))
!! n_p * n_q * v_pq^pq
!open(newunit=iunit,file='OOOOABBA',status='old',access='direct',&
!     form='unformatted',recl=8*dimOA**2)
!do ip=1,dimOB
!   do iq=1,dimOA
!       ! get all (AB| integrals for a given |BA) record
!      read(iunit,rec=ip+(iq-1)*dimOB) work2(1:dimOA,1:dimOB)
!     ! tvk(3) = tvk(3) + A%rdm1(iref,iq)*B%rdm1(iref2B,ip)*work2(iq,ip)
!      tvk(3) = tvk(3) + Aocc(iq)*Bocc(ip)*work2(iq,ip)
!   enddo
!enddo
!close(iunit)
!

! allocate(trdm2A(dimOA,dimOA,dimOA,dimOA))
! allocate(trdm2B(dimOB,dimOB,dimOB,dimOB))
! allocate(work2A(dimOA,dimOA,dimOA,dimOA))
! allocate(work2B(dimOB,dimOB,dimOB,dimOB))

! print*, "norm before transform",norm2(A%trdm24(iref,:,:,:,:))
! print*, "A%trdm(1,1,1,2,2)", A%trdm24(1,1,1,2,2) 
! do ip=1,dimOA
!    do iq=1,dimOA
!       call tran2MO(A%trdm24(iref,ip,iq,:,:),A%CMONO(iref,1:dimOA,1:dimOA),A%CMONO(iref,1:dimOA,1:dimOA), &
!        work2A(ip,iq,:,:),dimOA)
!    enddo
! enddo

! do ir=1,dimOA
!    do is=1,dimOA
!       call tran2MO(work2A(:,:,ir,is),A%CMONO(iref,1:dimOA,1:dimOA),A%CMONO(iref,1:dimOA,1:dimOA), &
!        trdm2A(:,:,ir,is),dimOA)
!    enddo
! enddo

! print*, "norm after transform", norm2(trdm2A)
! print*, " trdm2A(1,1,2,2)", trdm2A(1,1,2,2)

! switch to ground-state NOs for correct 2-TRDMs

 !
 !calculating TRDM2*TRDM1*V_{ee}*S
 !
 !
! allocate(intA(dimOA,dimOA,dimOA,dimOB),&
!          intB(dimOB,dimOB,dimOA,dimOB))


! print *,"Trdm1A vs TRdm1B in the same basis"
! do ip=1,3
!    do iq=1,3
!       print *,ip,iq,Atrdm(ip,iq),Btrdm2(ip,iq)
!    enddo
! enddo
! print*,"1TRDM A"
! do i=1,15
!    write(lout,'(*(f12.6))') (Atrdm(i,j),j=1,15)
! enddo
! print*,"1TRDM B"
! do i=1,15
!    write(lout,'(*(f12.6))') (Btrdm2(i,j),j=1,15)
! enddo
! call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%trdm24(iref2,:,:,:,:),dimOA**3,Sab2,NBasis,0d0,intA,dimOA**3)
! !call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2Bval,dimOB**3,Sab,NBas,0d0,intB,dimOB**3)

! ! careful! intB(B,B,A,B)
! do is=1,dimOB
!    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%trdm24(iref2B,:,:,:,is),dimOB**2,Sab2,NBasis,0d0,intB(:,:,:,is),dimOB**2)
! enddo
! !!!! WORK IN PROGRESS
!allocate(trdm2A(dimOA,dimOA,dimOA,dimOA))
!allocate(trdm2B(dimOB,dimOB,dimOB,dimOB))
! block
!    double precision :: DipXao(NBasis**2),DipYao(NBasis**2),DipZao(NBasis**2)
!    double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
!    double precision:: TSDipXYZ(3)
!    character(:),allocatable     :: mname,dipfile,mofile
   
!       mname  = 'A'
!       dipfile= "DIP_A"
!       mofile = 'MOLPRO_A.MOPUN'
   
!    call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

! !print*, 'DipX-AO:',norm2(DipXao)
! !print*, 'DipY-AO:',norm2(DipYao)
! !print*, 'DipZ-AO:',norm2(DipZao)

!    call unpack_sym_molpro(DipXao,dipfile,NBasis)
!    call unpack_sym_molpro(DipYao,dipfile,NBasis)
!    call unpack_sym_molpro(DipZao,dipfile,NBasis)

!    call tran2MO(DipXao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipX,NBasis)
!    call tran2MO(DipYao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipY,NBasis)
!    call tran2MO(DipZao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipZ,NBasis)
!    TSDipXYZ = 0d0
!    do j=1,NBasis
!       do i=1,NBasis
!          TSDipXYZ(1) = TSDipXYZ(1) - 2d0*Atrdm(j,i)*DipX(i,j)
!          TSDipXYZ(2) = TSDipXYZ(2) - 2d0*Atrdm(j,i)*DipY(i,j)
!          TSDipXYZ(3) = TSDipXYZ(3) - 2d0*Atrdm(j,i)*DipZ(i,j)
!       enddo
!    enddo
!    print *,'X,Y,Z moments from read Atrdm'
!    print *,TSDipXYZ

!    mname  = 'B'
!    dipfile= "DIP_B"
!    mofile = 'MOLPRO_B.MOPUN'

! call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

! !print*, 'DipX-AO:',norm2(DipXao)
! !print*, 'DipY-AO:',norm2(DipYao)
! !print*, 'DipZ-AO:',norm2(DipZao)

! call unpack_sym_molpro(DipXao,dipfile,NBasis)
! call unpack_sym_molpro(DipYao,dipfile,NBasis)
! call unpack_sym_molpro(DipZao,dipfile,NBasis)

! call tran2MO(DipXao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipX,NBasis)
! call tran2MO(DipYao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipY,NBasis)
! call tran2MO(DipZao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipZ,NBasis)
! TSDipXYZ = 0d0
! do j=1,NBasis
!    do i=1,NBasis
!       TSDipXYZ(1) = TSDipXYZ(1) - 2d0*Btrdm2(j,i)*DipX(i,j)
!       TSDipXYZ(2) = TSDipXYZ(2) - 2d0*Btrdm2(j,i)*DipY(i,j)
!       TSDipXYZ(3) = TSDipXYZ(3) - 2d0*Btrdm2(j,i)*DipZ(i,j)
!    enddo
! enddo
! print *,'X,Y,Z moments from read Btrdm'
! print *,TSDipXYZ

! end block

allocate(ints(dimOA,dimOA))
allocate(Atrdmtest(NBasis,NBasis),Btrdmtest(NBasis,NBasis))
Atrdmtest = 0d0
Btrdmtest = 0d0
do ip=1,NBasis
   do ir=1,NBAsis
      do iq=1,NBasis
         Atrdmtest(ip,ir)=  Atrdmtest(ip,ir)+Trdm2A(ip,ir,iq,iq)
      enddo
   enddo
enddo

do ip=1,NBasis
   do ir=1,NBAsis
      do iq=1,NBasis
         Btrdmtest(ip,ir)=Btrdmtest(ip,ir)+Trdm2B(ip,ir,iq,iq)
      enddo
   enddo
enddo

! print*,"1TRDM A Test"
! do i=1,15
!    write(lout,'(*(f12.6))') (Atrdmtest(i,j),j=1,15)
! enddo
! print*,"1TRDM B test"
! do i=1,15
!    write(lout,'(*(f12.6))') (Btrdmtest(i,j),j=1,15)
! enddo



!!! END OF WORK IN PROGRESS

!!call delfile('OOOOAABB') !! To nie wiem jak działa
! print*,"1 RDM A"
! do i=1,15
!    write(lout,'(*(f12.6))') (Atrdm(i,j),j=1,15)
! enddo
! print*,"1 RDM B"
! do i=1,15
!    write(lout,'(*(f12.6))') (Btrdm(i,j),j=1,15)
! enddo
! print*,"B%CAONO"
! do i=1,15
!    write(lout,'(*(f12.6))') (B%CAONO(1,i,j),j=1,15)
! enddo
! print*,"B%CMO"
! do i=1,15
!    write(lout,'(*(f12.6))') (B%CMO(i,j),j=1,15)
! enddo
! print*,"A%CMO"
! do i=1,15
!    write(lout,'(*(f12.6))') (A%CMO(i,j),j=1,15)
! enddo
! print*,"B%CMONO"
! do i=1,15
!    write(lout,'(*(f12.6))') (B%CMONO(1,i,j),j=1,15)
! enddo
! print*,"A%CMONO"
! do i=1,15
!    write(lout,'(*(f12.6))') (A%CMONO(1,i,j),j=1,15)
! enddo
! print*,"Sab2"
! do i=1,15
!    write(lout,'(*(f12.6))') (Sab2(i,j),j=1,15)
! enddo



! print*,"rdm24(1,2,2,2,2) dla A", trdm2A(2,2,2,2)
! print*,"rdm24(1,2,2,2,2) dla B", trdm2B(2,2,2,2)

call tran4_gen(NAO,&
     B%num0+B%num1,B%CAONO(irefB,1:NAO,1:(B%num0+B%num1)),&
     B%num0+B%num1,B%CAONO(irefB,1:NAO,1:(B%num0+B%num1)),&
     A%num0+A%num1,A%CAONO(iref, 1:NAO,1:(A%num0+A%num1)),&
     A%num0+A%num1,A%CAONO(iref, 1:NAO,1:(A%num0+A%num1)),&
     'OOOOAABB','AOTWOSORT')

TNa  = 0d0
TNa2 = 0d0
test = 0d0
do iu=1,dimOB
   do it=1,dimOB
      do iq=1,dimOA
         do ir=1,dimOA
            do ip=1,dimOA
               do is=1,dimOA
               !   tNa(1) = tNa(1) + B%Occ(it)*intA(ip,ir,iq,it)*Sab(iq,it)*Vbaa(ip,ir)
                  tNa(1) = tNa(1) + Btrdm(it,iu)*Trdm2A(ip,iq,ir,is)*Sab(ip,it)*Sab(iq,iu)*Vbaa(ir,is)
                  test = test + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*Sab2(ip,it)*Sab2(iq,iu)*Vbaa(ir,is)! z tego wynika, że raczej ze zmiany bazy sa spoko
               enddo   
            enddo
         enddo
      enddo
   enddo
enddo
tNa(1)=-2d0*tNa(1)!poprawione
test = -2d0 * test
tNa2(1)=tNa(1)
print *,"tNa1 liczone w innej bazie", test
print*,"1trdmB norm", norm2(Btrdm2)
print*,"1trdmA norm", norm2(Atrdm)
print*,"2trdmB norm", norm2(trdm2B)
print*,"2trdmA norm", norm2(trdm2A)

open(newunit=iunit,file='OOOOAAAB',status='OLD', access='DIRECT', &
  form='unformatted' ,recl=8*dimOA*dimOA)
! one loop over integrals
! get all (AA| integrals for a give |AB) record
ints = 0d0
do iu=1,dimOB
 do is=1,dimOA
    read(iunit,rec=is+(iu-1)*dimOA) ints(1:dimOA,1:dimOA)
  !     print* , "ints(1,1)", ints(1,1)
   do iq =1,dimOA
   do ir=1,dimOA
      do ip=1,dimOA
         do it=1,dimOB
         !tNa(2) = tNa(2) + B%Occ(it)*intA(ip,ir,iq,it)*ints(ip+(ir-1)*NBas)
         tNa(2) = tNa(2) + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
         !tNa(2) = tNa(2) + Btrdmtest(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
!            print *,tvk(3),Btrdm2(it,iu),intA(ip,ir,iq,iu),ints(ip,ir)
         enddo
      enddo
   enddo
  enddo
enddo
enddo
tNa(2) = -2d0*tNa(2) !poprawione
print*, 'tNa(2)',tNa(2)*1000
close(iunit)


open(newunit=iunit,file='OOOOAAAB',status='OLD', access='DIRECT', &
  form='unformatted' ,recl=8*dimOA*dimOA)
! one loop over integrals
! get all (AA| integrals for a give |AB) record
ints = 0d0
do it=1,dimOB
 do ir=1,dimOA
    read(iunit,rec=ir+(it-1)*dimOA) ints(1:dimOA,1:dimOA)
  !     print* , "ints(1,1)", ints(1,1)
   do iq =1,dimOA
   do is=1,dimOA
      do ip=1,dimOA
         do iu=1,dimOB
         !tNa(2) = tNa(2) + B%Occ(it)*intA(ip,ir,iq,it)*ints(ip+(ir-1)*NBas)
         tNa2(2) = tNa2(2) + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(is,iu)
         !tNa(2) = tNa(2) + Btrdmtest(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
!            print *,tvk(3),Btrdm2(it,iu),intA(ip,ir,iq,iu),ints(ip,ir)
         enddo
      enddo
   enddo
  enddo
enddo
enddo
tNa2(2) = -2d0*tNa2(2) !poprawione
print*, 'tNa2(2)',tNa2(2)*1000
close(iunit)

! tNb
tNb = 0d0
tNb2 = 0d0
do ip=1,dimOA
   do iq=1,dimOA
      do ir=1,dimOB
         do is=1,dimOB
            do it=1,dimOB
               do iu=1,dimOB
                  !tNb(1) = tNb(1) + A%Occ(it)*intB(ip,ir,it,iq)*Sab(it,iq)*Vabb(ip,ir)
                  !tNb(1) = tNb(1) + AOcc(it)*intB(ip,ir,it,iq)*Sab(it,iq)*Vabb(ip,ir)
                  tNb(1) = tNb(1)  + Atrdm(ip,iq)*Trdm2B(ir,is,it,iu)*Vabb2(it,iu)*Sab2(ip,ir)*Sab2(iq,is)
                  !tNa(1) = tNa(1) + Btrdm(it,iu)*Trdm2A(ip,iq,ir,is)*Sab(ip,it)*Sab(iq,iu)*Vbaa(ir,is)
               enddo
            enddo   
         enddo
      enddo
   enddo
enddo
tNb(1) = -2d0*tNb(1)
tNb2(1) = tNB(1)

! aded to dSRS
open(newunit=iunit,file='OOOOBBBA',status='OLD', &
     access='DIRECT',form='unformatted',recl=8*dimOB*dimOB)
!get all (BB| integrals for a given |BA) record
! one loop over integrals
ints = 0d0
do ip=1,dimOA
   do it=1,dimOB
      read(iunit,rec=it+(ip-1)*dimOB) ints(1:dimOB,1:dimOB)
   !  do iu=1,dimOB   
   !    read(iunit,rec=iu+(ip-1)*dimOB) ints(1:dimOB,1:dimOB)

 !     print* , "ints(1,1)", ints(1,1)
      do iq=1,dimOA
         do ir=1,dimOB
            do is =1,dimOB
              do iu=1,dimOB
               ! do it=1,dimOB
            !     tNb(2) = tNb(2) + A%Occ(it)*intB(ip,ir,it,iq)*ints(ip+(ir-1)*NBas)
                  tNb(2) = tNb(2) + Atrdm(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,iu)
                  !tNb(2) = tNb(2) + Atrdmtest(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,iu)
                  !TEST!
                  !tNa(2) = tNa(2) + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
                  ! tNb(2) = tNb(2) + Atrdm(iq,ip)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,it)
               enddo   
            enddo   
         enddo
      enddo

   enddo
enddo
tNb(2) = -2d0*tNb(2)
close(iunit)


open(newunit=iunit,file='OOOOBBBA',status='OLD', &
     access='DIRECT',form='unformatted',recl=8*dimOB*dimOB)
!get all (BB| integrals for a given |BA) record
! one loop over integrals
ints = 0d0
do iq=1,dimOA
   do iu=1,dimOB
      read(iunit,rec=iu+(iq-1)*dimOB) ints(1:dimOB,1:dimOB)
   !  do iu=1,dimOB   
   !    read(iunit,rec=iu+(ip-1)*dimOB) ints(1:dimOB,1:dimOB)

 !     print* , "ints(1,1)", ints(1,1)
      do ip = 1,dimOA
         do ir = 1,dimOB
            do is = 1,dimOB
              do it = 1,dimOB
               ! do it=1,dimOB
            !     tNb(2) = tNb(2) + A%Occ(it)*intB(ip,ir,it,iq)*ints(ip+(ir-1)*NBas)
                  tNb2(2) = tNb2(2) + Atrdm(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(ip,it)
                  !tNb(2) = tNb(2) + Atrdmtest(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,iu)
                  !TEST!
                  !tNa(2) = tNa(2) + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
                  ! tNb(2) = tNb(2) + Atrdm(iq,ip)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,it)
               enddo   
            enddo   
         enddo
      enddo

   enddo
enddo
TNb2(2) = -2d0*TNb2(2)
close(iunit)


print*, 'tNb-1',tNb(1)*1000
print*, 'tNb-2',tNb(2)*1000

call delfile('OOOOAABB')

call tran4_gen(NAO,&
     B%num0+B%num1,B%CAONO(irefB,1:NAO,1:(B%num0+B%num1)),&
     B%num0+B%num1,B%CAONO(irefB,1:NAO,1:(B%num0+B%num1)),&
     A%num0+A%num1,A%CAONO(iref, 1:NAO,1:(A%num0+A%num1)),&
     A%num0+A%num1,A%CAONO(iref, 1:NAO,1:(A%num0+A%num1)),&
     'OOOOAABB','AOTWOSORT')

! TERAZ PRACUJEMY TU

allocate(ints2(dimOA*dimOA))
allocate(workTEST(dimOA,dimOA))
allocate(intA(dimOA,dimOA,dimOA,dimOB),&
         intB(dimOB,dimOB,dimOA,dimOB))
intA=0d0
intB=0d0
! load reference 2-RDMs
! RDM2Aval = A%RDM2val
! RDM2Bval = B%RDM2val
! load reference 2-RDMs dSRS
call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,trdm2A(1:dimOA,1:dimOA,1:dimOA,1:dimOA),dimOA**3,Sab2,NBasis,0d0,intA,dimOA**3)

do is=1,dimOB
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,trdm2B(1:dimOB,1:dimOB,1:dimOB,is),dimOB**2,Sab2,NBasis,0d0,intB(:,:,:,is),dimOB**2)
enddo

open(newunit=iunit,file='OOOOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))
   
call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)
   val  = 0
   workTEST = 0

   ints2 = 0
   do ir=1,dimOB
      do ip=1,dimOB
        read(iunit,rec=ip+(ir-1)*dimOB) workTEST(1:dimOA,1:dimOA)
        val = val + sum(workTEST(1:dimOA,1:dimOA)*tmpAB(1:dimOA,1:dimOA,ip,ir))
      enddo
   enddo
   close(iunit)
   tNaNbTEST = -2d0*val
   print*, 'tNaNb from dgemm procedure',tNaNbTEST*1000
   TNaNb=tNaNbTEST
   TNaNb2=TNaNb
deallocate(intA,intB,tmpAB)

! !
! TNaNb = 0d0
! open(newunit=iunit,file='OOOOAABB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)
! do it=1,dimOB
!    do iu= 1,dimOB
!       read(iunit,rec=it+(iu-1)*dimOB) ints(1:dimOA,1:dimOA)
!       do iw=1,dimOB
!          do iv= 1,dimOB
!             do ip=1,dimOA
!               do iq=1,dimOA
!                  do ir = 1,dimOA
!                      do is = 1,dimOA
!                         TNaNb = TNaNB+Trdm2A(ip,iq,ir,is)*Trdm2B(it,iu,iw,iv)*ints(ip,iq)*Sab2(ir,iw)*Sab2(is,iv)
!                     enddo
!                  enddo
!               enddo     
!             enddo
!          enddo
!       enddo   
!    enddo
! enddo
! TNaNb=-2d0*TNaNb
! TNaNb2=TNaNb
! print *,"TNaNb from loop calculations", TNaNb*1000
!allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))

! call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)

!val  = 0d0
!do ir=1,dimOA
!   do ip=1,dimOA
    ! read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOB,1:dimOB)
     !val = val + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))
     !intA(ip,ir,iq,it)*intB(it,iu,is,iq)*ints(it,is)
!   enddo
!enddo

! ints = 0d0
! do ir=1,dimOB
!    do ip=1,dimOB
!      read(iunit,rec=ip+(ir-1)*dimOB) ints(1:dimOA,1:dimOA)
!      val = val + sum(ints(1:dimOA,1:dimOA)*tmpAB(1:dimOA,1:dimOA,ip,ir))
!    enddo
! enddo
! close(iunit)
! tvk(4) = -2d0*val!poprawione

VP  = sum(tvk)  + sum(TNa)  + sum(TNb)  + TNaNb  + P1*SAPT%Vnn
VP2 = sum(tvk2) + sum(TNa2) + sum(TNb2) + TNaNb2 + P1*SAPT%Vnn
exch1tot1 = SAPT%exchs21 - VP + P1*SAPT%elst1
exch1tot2 = SAPT%exchs22 + VP - P1*SAPT%elst2

A1B2PA2B1 = P1
write(LOUT,'(/1x,a,f16.8)') '<AB*|V|AB*>      = ',SAPT%A1B2VA1B2*1000d0
write(LOUT,'(1x,a,f16.8)')  '<A*B|V|A*B>      = ',SAPT%A2B1VA2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|V|A*B>      = ',SAPT%A1B2VA2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|P|AB*>      = ', A1B2PA1B2*1000d0
write(LOUT,'(1x,a,f16.8)')  '<A*B|P|A*B>      = ', A2B1PA2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|P|A*B>      = ', A1B2PA2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|VP|AB*>     = ',SAPT%A1B2VPA1B2*1000d0
write(LOUT,'(1x,a,f16.8)')  '<A*B|VP|A*B>     = ',SAPT%A2B1VPA2B1*1000d0
write(LOUT,'(1x,a,f16.8)')  '<AB*|VP|A*B>     = ',VP*1000d0
write(LOUT,'(1x,a,f16.8)')  '<A*B|VP|AB*>     = ',VP2*1000d0

if (SAPT%IPrint > 10) then
   write(lout,'(/1x,a)') 'Print all <VP> ingregients:'
   write(lout,'(1x,a,f16.10)') 'SAPT%exchs21', SAPT%exchs21
   write(lout,'(1x,a,f16.10)') 'SAPT%elst1  ', SAPT%elst1
   write(lout,'(1x,a,f16.10)') 'P1', P1
   write(lout,'(1x,a,f16.10)') 'SAPT%elst1-SAPT%Vnn', SAPT%elst1-SAPT%Vnn
   write(lout,'(1x,a,f16.10)') 'tvk(1)', tvk(1)
   write(lout,'(1x,a,f16.10)') 'tvk(2)', tvk(2)
   write(lout,'(1x,a,f16.10)') 'tvktest(1)', tvktest(1)
   write(lout,'(1x,a,f16.10)') 'tvktest(2)', tvktest(2)
   write(lout,'(1x,a,f16.10)') 'tvk(3) ', tvk(3)
   write(lout,'(1x,a,f16.10)') 'tvk2(1)',tvk2(1)
   write(lout,'(1x,a,f16.10)') 'tvk2(2)',tvk2(2)
   write(lout,'(1x,a,f16.10)') 'tvk2(3)',tvk2(3)
   write(lout,'(1x,a,f16.10)') 'sum tvk', sum(tvk)
   write(lout,'(1x,a,f16.10)') 'tNa(1) ',TNa(1)
   write(lout,'(1x,a,f16.10)') 'tNa(2) ',TNa(2)
   write(lout,'(1x,a,f16.10)') 'tNa2(1)',TNa2(1)
   write(lout,'(1x,a,f16.10)') 'tNa2(2)',tNa2(2)
   write(lout,'(1x,a,f16.10)') 'tNb(1) ',tNb(1)
   write(lout,'(1x,a,f16.10)') 'tNb(2) ',tNb(2)
   write(lout,'(1x,a,f16.10)') 'tNb2(1)',TNb2(1)
   write(lout,'(1x,a,f16.10)') 'tNb2(2)',TNb2(2)
   write(lout,'(1x,a,f16.10)') 'TNaNb  ',TNaNb
   write(lout,'(1x,a,f16.10)') 'TNaNb2 ',TNaNb2
   write(lout,'(1x,a,f16.10)') 'elst dSRS 1 = ', SAPT%elst1*1000
   write(lout,'(1x,a,f16.10)') 'elst dSRS 2 = ', SAPT%elst2*1000
   write(lout,'(1x,a,f16.10)') 'exchange elst dSRS 1 = ', exch1tot1*1000
   write(lout,'(1x,a,f16.10)') 'exchange elst dSRS 2 = ', exch1tot2*1000
endif

write(lout,'(/1x,a)') 'Now we calculate matrix elements of V in Psi basis:'
gamma = SAPT%gamma
V(1,1)=COS(gamma)**2*SAPT%A1B2VA1B2+SIN(gamma)**2*SAPT%A2B1VA2B1+2*SIN(gamma)*COS(gamma)*SAPT%A1B2VA2B1
V(2,2)=SIN(gamma)**2*SAPT%A1B2VA1B2+COS(gamma)**2*SAPT%A2B1VA2B1-2*SIN(gamma)*COS(gamma)*SAPT%A1B2VA2B1
V(1,2)=-SIN(gamma)*COS(gamma)*SAPT%A1B2VA1B2+COS(gamma)*SIN(gamma)*SAPT%A2B1VA2B1+(COS(gamma)**2-SIN(gamma)**2)*SAPT%A1B2VA2B1
V(2,1)=-SIN(gamma)*COS(gamma)*SAPT%A1B2VA1B2+COS(gamma)*SIN(gamma)*SAPT%A2B1VA2B1+(COS(gamma)**2-SIN(gamma)**2)*SAPT%A1B2VA2B1
do i=1,2
   write(lout,'(*(f12.6))') (V(i,j),j=1,2)
enddo

write(lout, '(/1x,a)') 'Calculate matrix elements of P in Psi basis:'
P(1,1)=COS(gamma)**2*A1B2PA1B2+SIN(gamma)**2*A2B1PA2B1+2*SIN(gamma)*COS(gamma)*A1B2PA2B1
P(2,2)=SIN(gamma)**2*A1B2PA1B2+COS(gamma)**2*A2B1PA2B1-2*SIN(gamma)*COS(gamma)*A1B2PA2B1
P(1,2)=-SIN(gamma)*COS(gamma)*A1B2PA1B2+COS(gamma)*SIN(gamma)*A2B1PA2B1+(COS(gamma)**2-SIN(gamma)**2)*A1B2PA2B1
P(2,1)=-SIN(gamma)*COS(gamma)*A1B2PA1B2+COS(gamma)*SIN(gamma)*A2B1PA2B1+(COS(gamma)**2-SIN(gamma)**2)*A1B2PA2B1
do i=1,2
   write(lout,'(*(f12.6))') (P(i,j),j=1,2)
enddo

write(lout,'(/1x,a)') 'Calculate matrix elements of VP in Psi basis:'
VPnb(1,1)=COS(gamma)**2*SAPT%A1B2VPA1B2+SIN(gamma)**2*SAPT%A2B1VPA2B1+SIN(gamma)*COS(gamma)*(VP+VP2)
VPnb(2,2)=SIN(gamma)**2*SAPT%A1B2VPA1B2+COS(gamma)**2*SAPT%A2B1VPA2B1-SIN(gamma)*COS(gamma)*(VP+VP2)
VPnb(1,2)=-SIN(gamma)*COS(gamma)*SAPT%A1B2VPA1B2+COS(gamma)*SIN(gamma)*SAPT%A2B1VPA2B1+COS(gamma)**2*VP-SIN(gamma)**2*VP2
VPnb(2,1)=-SIN(gamma)*COS(gamma)*SAPT%A1B2VPA1B2+COS(gamma)*SIN(gamma)*SAPT%A2B1VPA2B1+COS(gamma)**2*VP2-SIN(gamma)**2*VP

do i=1,2
   write(lout,'(*(f12.6))') (VPnb(i,j),j=1,2)
enddo

write(lout, '(/1x,a)') "Electrostatic energy from matrix element(in millihartree)"
print *,"First ", V(1,1)*1000d0
print *,"Second", V(2,2)*1000d0

write(lout, '(/1x,a)') "Exch energy from matrix element(in millihartree)"
print *,"First ", (VPnb(1,1)-P(1,1)*V(1,1))*1000d0
print *,"Second", (VPnb(2,2)-P(2,2)*V(2,2))*1000d0
write(lout, '(/1x,a)') "Interaction energy in S2 approx."
print *,"Eint1", (V(1,1) + VPnb(1,1) - P(1,1)*V(1,1))*1000d0
print *,"Eint2", (V(2,2) + VPnb(2,2) - P(2,2)*V(2,2))*1000d0
write(lout, '(/1x,a)') "Coefficients from diagonalizing V:"
print *,"For first state:"
print *,"C11: ",COS(gamma)
print *,"C12: ",SIN(gamma)
print *,"For second state:"
print *,"C21: ", -SIN(gamma)
print *,"C12: ", COS(gamma)

gamma_exch1 = -atan2(VPnb(2,1)-V(1,1)*P(2,1),V(2,2)-V(1,1))
gamma_exch1deg = 180 * gamma_exch1 / pi

print *,"gamma_exch1",gamma_exch1

write(lout,'(/1x,a)') "Coefficients from exchange diagonalization:"
print *,"For first state:"
print *,"C11: ",COS(gamma+gamma_exch1)
print *,"C12: ",SIN(gamma+gamma_exch1)

gamma_exch2=atan2(VPnb(1,2)-V(2,2)*P(1,2),V(1,1)-V(2,2))
gamma_exch2deg = 180 * gamma_exch2 / pi

print *,"gamma_exch2",gamma_exch2

write(lout,'(/1x,"Delta-1 Mixing angle (exch),rad", f12.6)') gamma_exch1
write(lout,'( 1x,"Delta-1 Mixing angle (exch),deg", f9.3)')  gamma_exch1deg
write(lout,'(/1x,"Delta-2 Mixing angle (exch),rad", f12.6)') gamma_exch2
write(lout, '(1x,"Delta-2 Mixing angle (exch),deg", f9.3,/)')  gamma_exch2deg

print *,"For second state:"
print *,"C21: ", -sin(gamma+gamma_exch2)
print *,"C12: ", cos(gamma+gamma_exch2)

if (SAPT%IPrint >= 10) then
   write(lout, '(/1x,a)') "How good is trygonometric approximation?:"
   print *,"We can check is the vec proposed above an eigenvector or not:"
   print *,(V(1,1)+VPnb(1,1))*cos(gamma_exch1)+(V(1,2)+VPnb(1,2))*sin(gamma_exch1)
   print *,(V(2,1)+VPnb(2,1))*cos(gamma_exch1)+(V(2,2)+VPnb(2,2))*sin(gamma_exch1)
   print *,"VERSUS"
   print *,(1+P(1,1))*cos(gamma_exch1)*((VPnb(1,1)-P(1,1)*V(1,1))+V(1,1))+P(1,2)*sin(gamma_exch1)*((VPnb(1,1)-P(1,1)*V(1,1))+V(1,1))
   print *,P(1,2)*(VPnb(1,1)-P(1,1)*V(1,1)+V(1,1))*cos(gamma_exch1)+(1+P(2,2))*sin(gamma_exch1)*((VPnb(1,1)-P(1,1)*V(1,1))+V(1,1))
   print *,"and"
   print *,-(V(1,1)+VPnb(1,1))*sin(gamma_exch2)+(V(1,2)+VPnb(1,2))*cos(gamma_exch2)
   print *,-(V(2,1)+VPnb(2,1))*sin(gamma_exch2)+(V(2,2)+VPnb(2,2))*cos(gamma_exch2)
   print *,"VERSUS"
   print *,-(1+P(1,1))*sin(gamma_exch2)*((VPnb(2,2)-P(2,2)*V(2,2))+V(2,2))+P(1,2)*cos(gamma_exch2)*((VPnb(2,2)-P(2,2)*V(2,2))+V(2,2))
   print *,-P(1,2)*sin(gamma_exch2)*((VPnb(2,2)-P(2,2)*V(2,2))+V(2,2))+(1+P(2,2))*cos(gamma_exch2)*((VPnb(2,2)-P(2,2)*V(2,2))+V(2,2))
endif

block
   integer :: info, i
   double precision:: eigenvalues(2) ! Eigenvalues
   double precision :: Seigenvectors(2,2) ! Eigenvectors
   double precision :: eigenvectors(2,2)
   double precision :: workX(10)
   double precision :: EigI(2)
   double precision :: eigenvectorsL(2,2) 
   double precision :: Sinverse(2,2)
   double precision :: S22(2,2)
   double precision :: S22h(2,2)
   double precision :: Vh(2,2)
   double precision :: Sdet
   double precision :: Check(2,2)
   double precision :: VSinverse(2,2)

   double precision :: eigenvectors2(2,2)
   double precision :: workX2(16)
   double precision :: eigenvectorsL2(2,2) 
   double precision :: alphar(2), alphai(2), beta(2)
   double precision:: eigenvalues2(2)


   ! FOR TESTING PURPOSE angle 90
   ! V(1,1) = 0.009059
   ! V(2,2) = 0.010340
   ! V(1,2) = 0.001488
   ! V(2,1) = 0.002243
   ! P(1,1) = 0d0
   ! P(2,2) = 0d0
   ! P(1,2) = 0.004653
   ! P(2,1) = 0.004653
   ! WITHOUT EXCHANGE
   ! V(1,1) = -0.030221
   ! V(2,2) =  0.001070
   ! V(1,2) =  0.000759
   ! V(2,1) =  0.000759
   ! P(1,1) = 0d0
   ! P(2,2) = 0d0
   ! P(1,2) = 0d0
   ! P(2,1) = 0d0
   ! 150
   ! V(1,1) =  0.008447
   ! V(2,2) =  0.010390
   ! V(1,2) =  0.000761
   ! V(2,1) =  0.001153
   ! P(1,1) = 0d0
   ! P(2,2) = 0d0
   ! P(1,2) = 0.002438
   ! P(2,1) = 0.002438
   ! Helg
   ! V(1,1) = -2.861627
   ! V(2,2) = -2.169652
   ! V(1,2) = -0.825295
   ! V(2,1) = -0.825295
   ! P(1,1) = 0d0
   ! P(2,2) = 0d0
   ! P(1,2) = 0.280609
   ! P(2,1) = 0.280609


   print *,"Diagonalizing full eigenvalue problem"
   S22h(1,1) = 1 + A2B1PA2B1
   S22h(2,2) = 1 + A1B2PA1B2
   S22h(2,1) = A1B2PA2B1
   S22h(1,2) = A1B2PA2B1
   Vh(1,1)   = SAPT%A2B1VA2B1 + SAPT%A2B1VPA2B1
   Vh(2,2)   = SAPT%A1B2VA1B2 + SAPT%A1B2VPA1B2
   Vh(1,2)   = SAPT%A1B2VA2B1 + VP2
   Vh(2,1)   = SAPT%A1B2VA2B1 + VP


   call dggev('N', 'V', 2, Vh, 2, S22h, 2, alphar, alphai, beta ,eigenvectorsL2,2, eigenvectors2,2, workX2, 16, info)
   write(*, *) "Is alphai=0?"
   do i = 1, 2
   write(*, '(F10.7)') alphai(i)
   end do

   write(*, *) "Eigenvalues2 from dgevv(in milihartree)"
   do i = 1, 2
   eigenvalues2(i)= alphar(i)/beta(i)
   write(*, '(F10.7)') eigenvalues2(i)
   end do

   write(*, *) "Eigenvectors2 from dgevv:"
   do i = 1,2
      write(*, '(2F8.4)') eigenvectors2(i,1)/SQRT(eigenvectors2(1,1)**2+eigenvectors2(2,1)**2), eigenvectors2(i,2)/SQRT(eigenvectors2(1,2)**2+eigenvectors2(2,2)**2)
   end do

   ! S22(1,1)=1+P(1,1)
   ! S22(2,2)=1+P(2,2)
   ! S22(1,2)=P(1,2)
   ! S22(2,1)=P(2,1)
   ! Sdet = (1+P(1,1))*(1+P(2,2))-P(1,2)*P(1,2)
   ! Sinverse(1,1) =(1+P(2,2))/Sdet
   ! Sinverse(1,2) =-P(1,2)/Sdet
   ! Sinverse(2,1) =-P(1,2)/Sdet
   ! Sinverse(2,2) =(1+P(1,1))/Sdet
   ! call dgemm('N','N', 2 , 2 , 2 , 1d0, S22, 2, Sinverse, 2, 0,Check ,2)
   ! write(*, *) "Check:"
   ! do i = 1,2
   ! write(*, '(2F8.4)') Check(i,1), Check(i,2)
   ! end do
   ! call dgemm('N','N', 2 , 2 , 2 , 1d0, Sinverse, 2, S22, 2, 0,Check ,2)
   ! do i = 1,2
   !    write(*, '(2F8.4)') Check(i,1), Check(i,2)
   ! end do
   ! call dgemm('N','N', 2 , 2 , 2 , 1d0, V, 2, Sinverse, 2, 0,VSinverse ,2)
   ! call DGEEV('N', 'V', 2, VSinverse , 2, eigenvalues, EigI ,eigenvectorsL,2, Seigenvectors,2, workX, 10, info)
   ! !CALL SGEEV | DGEEV (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
   ! call dgemm('N','N', 2 , 2 , 2 , 1d0, Sinverse , 2, Seigenvectors, 2, 0, eigenvectors ,2)
   ! ! Print the eigenvalues
   ! write(*, *) "Eigenvalues:"
   ! do i = 1, 2
   ! write(*, '(F10.7)') eigenvalues(i)
   ! end do

   ! write(*, *) "Eigenvectors:"
   ! do i = 1,2
   !    write(*, '(2F8.4)') eigenvectors(i,1)/SQRT(eigenvectors(1,1)**2+eigenvectors(2,1)**2), eigenvectors(i,2)/SQRT(eigenvectors(1,2)**2+eigenvectors(2,2)**2)
   ! end do

end block

deallocate(Sab,Sab2,Sab3)
deallocate(Vabb,Vbaa,Vaab,Vbba,Vbab)
deallocate(Vbab2,Vabb2,Vaab2)

deallocate(work)
deallocate(ints)
deallocate(Trdm2B,Trdm2A)
deallocate(Btrdmtest,Atrdmtest)

end subroutine e1exch_dSRS

end module sapt_dRS
