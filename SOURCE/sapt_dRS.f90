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
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBasis,dimOA,dimOB
integer :: iunit
integer :: iref,iref2,irefB,iref2B
integer :: iexcited,iexcited2
integer :: ip,ir,iq,is,j
double precision :: elab, elab2,elab3,tr
double precision,allocatable :: work(:,:)
double precision,allocatable :: Btrdm(:,:)
double precision,allocatable :: Atrdm(:,:)

integer :: i
double precision,allocatable :: Vb(:,:)
double precision,allocatable :: Vbaa(:,:)
double precision,allocatable :: Vbaa2(:,:)
double precision :: ea1,ea2
double precision :: elst1,elst2,elstSAPT

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

NBasis = A%NBasis
dimOA  = A%num0+A%num1
dimOB  = B%num0+B%num1
print*, '1-TRDM MO ='
      do i=1,10
         write(lout,'(*(f12.6))') (A%trdm1(1,i,j),j=1,10)
      enddo
print*,'Atrdm1(1,2)',A%trdm1(1,1,2)
print*,'Atrdm1(1,3)',A%trdm1(1,1,3)

print*,'trace'
tr=0d0
do i=1,NBasis
   tr=tr+A%trdm1(iexcited,i,i)
enddo
print*,tr
! we will need :
!       a) A%rdm1(istate,NBasis)
!       b) A%trdm1(nnstate,NBasis,NBasis)
!       c) v_{pr}^{qs} integrals, pq \in A, rs \in B
!       d) SAPT%Vnn
!
! let us get the integral files here,
!  we will move them outside to sapt_ab_ints later

! first, try to recover regular Elst for A-B (ground-state)
!
! get v_pr^qs

allocate(Vb(NBasis,NBasis),Vbaa(NBasis,NBasis))
allocate(Vbaa2(NBasis,NBasis))
call get_one_mat('V',Vb,B%Monomer,NBasis)

call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
call tran2MO(Vb,A%CAONO(iref2,:,:),A%CAONO(iref2,:,:),Vbaa2,NBasis)
! sum_p n_p v^B_pp
ea1 = 0
ea2 = 0
do i=1,A%num0+A%num1
   ea1 = ea1 + A%rdm1(iref,i)*Vbaa(i,i)
   ea2 = ea2 + A%rdm1(iref2,i)*Vbaa2(i,i)
enddo
ea1 = 2d0*ea1
ea2 = 2d0*ea2
print*, 'ea1',ea1
print*, 'ea2',ea2

!call tran4_gen(NBasis,&
!               B%num0+B%num1,B%CMO,&
!               B%num0+B%num1,B%CMO,&
!               A%num0+A%num1,A%CMO,&
!               A%num0+A%num1,A%CMO,&
!               'OOOOAABB','AOTWOSORT')

!call tran4_gen(NBasis,&
!     B%num0+B%num1,B%CAONO(iref2B,:,:),&
!     B%num0+B%num1,B%CAONO(iref2B,:,:),&
!     A%num0+A%num1,A%CAONO(iref,:,:),&
!     A%num0+A%num1,A%CAONO(iref,:,:),&
!     'OOOOAABB','AOTWOSORT')

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

! test
print*, 'elab = ',4d0*elab
 
allocate(Atrdm(NBasis,NBasis))
allocate(Btrdm(NBasis,NBasis))
! Does it make sense???
!call tran2MO(A%trdm(iref,:,:),CMONO(RefState,:,:),CMONO(RefState,:,:), &
!                        Atrdm(iref,:,:),NBasis)
call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
                        Atrdm(:,:),NBasis)
call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(iref2B,:,:),B%CMONO(iref2B,:,:), &
                             Btrdm(:,:),NBasis)


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
print *, 'elab3 =',elab3
!print *, 'EdRS(1)= ' , (ea1+ea2+elab2+ABS(elab3)+SAPT%Vnn)*1d3
!print *, 'EdRS(2)= ' , (ea1+ea2+elab2-ABS(elab3)+SAPT%Vnn)*1d3

elst1 = ea1+ea2+elab2+elab3+SAPT%Vnn
elst2 = ea1+ea2+elab2-elab3+SAPT%Vnn
elstSAPT = ea1+ea2+elab2+SAPT%Vnn
!Na potrzeby testotwe
SAPT%elst=elstSAPT
print *, 'EdRS(1)= ' , elst1*1d3
print *, 'EdRS(2)= ' , elst2*1d3
print *, 'elstSAPT= ', elstSAPT*1d3
print *, 'Vnn', SAPT%Vnn
!SAPT%elst  = elst1
SAPT%elst1 = elst1
SAPT%elst2 = elst2

deallocate(work)

end subroutine elst_dRS

subroutine e1exch_dSRS(A,B,SAPT)
!
! calculate exchange energy
! in degenerate SRS (in NO representation)
!
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
integer :: NBasis,dimOA,dimOB
integer :: iunit
integer :: iref,iref2,irefB,iref2B
integer :: iexcited,iexcited2
integer :: ip,ir,iq,is,iu,it,iv,iw,i,j
double precision, allocatable :: trdm2A(:,:,:,:)
double precision, allocatable :: trdm2B(:,:,:,:)
double precision :: tvk(3),TNa(2),TNb(2),TNaNb
double precision, allocatable :: Atrdm(:,:)
double precision, allocatable :: Btrdm(:,:),Btrdm2(:,:)
double precision, allocatable :: Atrdmtest(:,:),Btrdmtest(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Sab2(:,:),Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:),Vbab(:,:),Vabb2(:,:)
double precision, allocatable :: work(:,:)
double precision, allocatable :: work2A(:,:,:,:), work2B(:,:,:,:)
! double precision, allocatable :: intA(:,:,:,:), intB(:,:,:,:)
double precision, allocatable :: ints(:,:)
double precision, allocatable :: tmpAB(:,:,:,:)
double precision :: val, P1, exch1tot1, exch1tot2, VP,test
print *, "running exch dSRS"
print *,"exch w DRS ",SAPT%exchs21


print*,"2TRDMA"
do ip=1,10
   do iq = 1,10
      do ir = 1,10
         do is =1,10
            if(abs(A%trdm24(1,ip,iq,ir,is))>1d-2) then
               print*,ip,iq,ir,is,A%trdm24(1,ip,iq,ir,is)
            endif
         enddo
      enddo
   enddo
enddo


print*,"2TRDMB"
do ip=1,10
   do iq = 1,10
      do ir = 1,10
         do is =1,10
            if(abs(A%trdm24(1,ip,iq,ir,is))>1d-2) then
               print*,ip,iq,ir,is,B%trdm24(1,ip,iq,ir,is)
            endif
         enddo
      enddo
   enddo
enddo


print*,"Czy 2TRDMB+2TRDMA=0"
do ip=1,10
   do iq = 1,10
      do ir = 1,10
         do is =1,10
            if(abs(A%trdm24(1,ip,iq,ir,is)+B%trdm24(1,ip,iq,ir,is))>1d-2) then
               print*,ip,iq,ir,is,B%trdm24(1,ip,iq,ir,is)+A%trdm24(1,ip,iq,ir,is)
            endif
         enddo
      enddo
   enddo
enddo

!CZYTAMY KOD NA NOWO

iref   =  A%IREF1 !1
iref2  =  A%IREF2 !4
irefB  =  B%IREF1
iref2B  = B%IREF2
iexcited = ((iref2-2)*(iref2-1))/2+iref
iexcited2 = ((iref2B-2)*(iref2B-1))/2+irefB
tvk = 0d0

NBasis = A%NBasis
dimOA  = A%num0+A%num1
dimOB  = B%num0+B%num1

allocate(Atrdm(NBasis,NBasis))
allocate(Btrdm(NBasis,NBasis),Btrdm2(NBasis,NBasis))

allocate(S(NBasis,NBasis),Sab(NBasis,NBasis),Sab2(NBasis,NBasis))
allocate(Va(NBasis,NBasis),Vb(NBasis,NBasis),&
         Vabb(NBasis,NBasis),Vbaa(NBasis,NBasis),&
         Vaab(NBasis,NBasis),Vbba(NBasis,NBasis),Vbab(NBasis,NBasis),Vabb2(NBasis,NBasis))


! get V_ne in atomic orbs
call get_one_mat('V',Va,A%Monomer,NBasis)
call get_one_mat('V',Vb,B%Monomer,NBasis)



!call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
! not sure all are needed
call tran2MO(Va,B%CAONO(iref2B,:,:),B%CAONO(iref2B,:,:),Vabb,NBasis)
call tran2MO(Va,B%CAONO(irefB,:,:),B%CAONO(irefB,:,:),Vabb2,NBasis)
call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
call tran2MO(Va,A%CAONO(iref,:,:),B%CAONO(iref2B,:,:),Vaab,NBasis)
call tran2MO(Vb,B%CAONO(iref2B,:,:),A%CAONO(iref,:,:),Vbba,NBasis)
call tran2MO(Vb,A%CAONO(iref,:,:),B%CAONO(iref2B,:,:),Vbab,NBasis)
print*,"norm Va before and after transformation",norm2(Va), norm2(Vabb)


! get overlap S matrix in AO and transform to NOs
call get_one_mat('S',S,A%Monomer,NBasis)

!call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
call tran2MO(S,A%CAONO(iref,:,:),B%CAONO(iref2B,:,:),Sab,NBasis)


call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
                        Atrdm(:,:),NBasis)
call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(iref2B,:,:),B%CMONO(iref2B,:,:), &
                             Btrdm(:,:),NBasis)
print*,"norm S before and after transformation",norm2(S), norm2(Sab)
print*,"norm Atrdm before and after transformation",norm2(A%trdm1(iexcited,:,:)), norm2(Atrdm)
!print*, "Sab",Sab

!tvk(1)=- 4sum(trdm(A->A*)_{pq}trdm(B*->B)_{rs}VB_{ps}S_{qr})
!P_1=-2sum(trdm(A->A*)_{pq}trdm(B*->B)_{rs}S_{ps}S_{qr})
P1=0d0
tvk = 0d0
do iq=1,dimOA
   do ip=1,dimOA
      do ir=1,dimOB
         do is=1,dimOB
            !     tvk(1) = tvk(1) + A%Occ(ip)*B%Occ(iq)*Vaab(ip,iq)*Sab(ip,iq)
            P1 = P1 +Atrdm(ip,iq)*Btrdm(ir,is)*Sab(ip,ir)*Sab(iq,is) !poprawione
            tvk(1) = tvk(1) + Atrdm(ip,iq)*Btrdm(is,ir)*Vbab(ip,is)*Sab(iq,ir) !poprawione
            tvk(2) = tvk(2) + Atrdm(ip,iq)*Btrdm(is,ir)*Vaab(ip,is)*Sab(iq,ir)
         Enddo
      enddo
   Enddo
enddo
tvk(1) = -2d0*tvk(1)
tvk(2) = -2d0*tvk(2)
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

 
call tran2MO(S,A%CAONO(iref,:,:),B%CAONO(irefB,:,:),Sab2,NBasis)

call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(irefB,:,:),B%CMONO(irefB,:,:), &
                             Btrdm2(:,:),NBasis)

print *,"Trdm1A vs TRdm1B in the same basis"
do ip=1,3
   do iq=1,3
      print *,ip,iq,Atrdm(ip,iq),Btrdm2(ip,iq)
   enddo
enddo
print*,"1TRDM A"
do i=1,15
   write(lout,'(*(f12.6))') (Atrdm(i,j),j=1,15)
enddo
print*,"1TRDM B"
do i=1,15
   write(lout,'(*(f12.6))') (Btrdm2(i,j),j=1,15)
enddo
! call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%trdm24(iref2,:,:,:,:),dimOA**3,Sab2,NBasis,0d0,intA,dimOA**3)
! !call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2Bval,dimOB**3,Sab,NBas,0d0,intB,dimOB**3)

! ! careful! intB(B,B,A,B)
! do is=1,dimOB
!    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%trdm24(iref2B,:,:,:,is),dimOB**2,Sab2,NBasis,0d0,intB(:,:,:,is),dimOB**2)
! enddo
! !!!! WORK IN PROGRESS
!allocate(trdm2A(dimOA,dimOA,dimOA,dimOA))
!allocate(trdm2B(dimOB,dimOB,dimOB,dimOB))
block
   double precision :: DipXao(NBasis**2),DipYao(NBasis**2),DipZao(NBasis**2)
   double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
   double precision:: TSDipXYZ(3)
   character(:),allocatable     :: mname,dipfile,mofile
   
      mname  = 'A'
      dipfile= "DIP_A"
      mofile = 'MOLPRO_A.MOPUN'
   
   call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

!print*, 'DipX-AO:',norm2(DipXao)
!print*, 'DipY-AO:',norm2(DipYao)
!print*, 'DipZ-AO:',norm2(DipZao)

   call unpack_sym_molpro(DipXao,dipfile,NBasis)
   call unpack_sym_molpro(DipYao,dipfile,NBasis)
   call unpack_sym_molpro(DipZao,dipfile,NBasis)

   call tran2MO(DipXao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipX,NBasis)
   call tran2MO(DipYao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipY,NBasis)
   call tran2MO(DipZao,A%CAONO(A%iref1,:,:),A%CAONO(A%iref1,:,:),DipZ,NBasis)
   TSDipXYZ = 0d0
   do j=1,NBasis
      do i=1,NBasis
         TSDipXYZ(1) = TSDipXYZ(1) - 2d0*Atrdm(j,i)*DipX(i,j)
         TSDipXYZ(2) = TSDipXYZ(2) - 2d0*Atrdm(j,i)*DipY(i,j)
         TSDipXYZ(3) = TSDipXYZ(3) - 2d0*Atrdm(j,i)*DipZ(i,j)
      enddo
   enddo
   print *,'X,Y,Z moments from read Atrdm'
   print *,TSDipXYZ

   mname  = 'B'
   dipfile= "DIP_B"
   mofile = 'MOLPRO_B.MOPUN'

call read_dip_sym_molpro(DipXao,DipYao,DipZao,dipfile,NBasis)

!print*, 'DipX-AO:',norm2(DipXao)
!print*, 'DipY-AO:',norm2(DipYao)
!print*, 'DipZ-AO:',norm2(DipZao)

call unpack_sym_molpro(DipXao,dipfile,NBasis)
call unpack_sym_molpro(DipYao,dipfile,NBasis)
call unpack_sym_molpro(DipZao,dipfile,NBasis)

call tran2MO(DipXao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipX,NBasis)
call tran2MO(DipYao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipY,NBasis)
call tran2MO(DipZao,B%CAONO(B%iref1,:,:),B%CAONO(B%iref1,:,:),DipZ,NBasis)
TSDipXYZ = 0d0
do j=1,NBasis
   do i=1,NBasis
      TSDipXYZ(1) = TSDipXYZ(1) - 2d0*Btrdm2(j,i)*DipX(i,j)
      TSDipXYZ(2) = TSDipXYZ(2) - 2d0*Btrdm2(j,i)*DipY(i,j)
      TSDipXYZ(3) = TSDipXYZ(3) - 2d0*Btrdm2(j,i)*DipZ(i,j)
   enddo
enddo
print *,'X,Y,Z moments from read Btrdm'
print *,TSDipXYZ

end block

















allocate(trdm2A(NBasis,NBasis,NBasis,NBasis))
allocate(trdm2B(NBasis,NBasis,NBasis,NBasis))
! trdm2A(:,:,:,:)=A%trdm24(iref2,:,:,:,:)! TO TYLKO DO TESTÓW TRZEBA ZMIENIĆ PÓŹNIEJ
! trdm2B(:,:,:,:)=B%trdm24(iref2B,:,:,:,:)

trdm2A(:,:,:,:)=A%trdm24(1,:,:,:,:)
trdm2B(:,:,:,:)=B%trdm24(1,:,:,:,:)

print *,"dimOA",dimOA
allocate(ints(dimOA,dimOA))
allocate(Atrdmtest(NBasis,NBasis),Btrdmtest(NBasis,NBasis))
Atrdmtest = 0d0
Btrdmtest = 0d0
do ip=1,NBasis
   do ir=1,NBAsis
      do iq=1,NBasis
         Atrdmtest(ip,ir)=  Atrdmtest(ip,ir)+trdm2A(ip,ir,iq,iq)
      enddo
   enddo
enddo

do ip=1,NBasis
   do ir=1,NBAsis
      do iq=1,NBasis
         Btrdmtest(ip,ir)=Btrdmtest(ip,ir)+trdm2B(ip,ir,iq,iq)
      enddo
   enddo
enddo

print*,"1TRDM A Test"
do i=1,15
   write(lout,'(*(f12.6))') (Atrdmtest(i,j),j=1,15)
enddo
print*,"1TRDM B test"
do i=1,15
   write(lout,'(*(f12.6))') (Btrdmtest(i,j),j=1,15)
enddo



!!! END OF WORK IN PROGRESS

!!call delfile('OOOOAABB') !! To nie wiem jak działa

! Trik do testowania
Atrdm=0d0
Btrdm=0d0
Btrdm2=0d0
trdm2A=0d0
trdm2B=0d0
do ip=1,dimOA
   do iq=1,dimOB
      if(ip == iq) then
      Atrdm(ip,iq) = A%rdm1(1,ip)
      Btrdm2(ip,iq) = B%rdm1(1,ip)
      Btrdm(ip,iq) = B%rdm1(1,ip)
      endif
   enddo
enddo
trdm2A(1:dimOA,1:dimOA,1:dimOA,1:dimOA)=A%rdm24(1,:,:,:,:)
trdm2B(1:dimOA,1:dimOA,1:dimOA,1:dimOA)=B%rdm24(1,:,:,:,:)

print*,"1 RDM A"
do i=1,15
   write(lout,'(*(f12.6))') (Atrdm(i,j),j=1,15)
enddo
print*,"1 RDM B"
do i=1,15
   write(lout,'(*(f12.6))') (Btrdm(i,j),j=1,15)
enddo
print*,"rdm24(1,2,2,2,2) dla A", trdm2A(2,2,2,2)
print*,"rdm24(1,2,2,2,2) dla B", trdm2B(2,2,2,2)

call tran4_gen(NBasis,&
     B%num0+B%num1,B%CAONO(irefB,1:NBasis,1:(B%num0+B%num1)),&
     B%num0+B%num1,B%CAONO(irefB,1:NBasis,1:(B%num0+B%num1)),&
     A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
     A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
     'OOOOAABB','AOTWOSORT')









test = 0d0
TNa = 0d0
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
print *,"tNa1 licznoe w innej bazie", test
print*,"1trdmB norm",norm2(Btrdm2)
print*,"1trdmA norm", norm2(Atrdm)
print*,"2trdmB norm",norm2(trdm2B)
print*,"2trdmA norm", norm2(trdm2A)







open(newunit=iunit,file='OOOOAAAB',status='OLD', access='DIRECT', &
  form='unformatted' ,recl=8*dimOA*dimOA)
! one loop over integrals
! get all (AA| integrals for a give |AB) record
ints = 0d0
do iu=1,dimOB
 do is=1,dimOA
    read(iunit,rec=is+(iu-1)*dimOA) ints(1:dimOA,1:dimOA)
       print* , "ints(1,1)", ints(1,1)
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


! tNb
tNb = 0
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

! aded to dSRS
open(newunit=iunit,file='OOOOBBBA',status='OLD', &
     access='DIRECT',form='unformatted',recl=8*dimOB*dimOB)
!get all (BB| integrals for a given |BA) record
! one loop over integrals
ints = 0d0
do ip=1,dimOA
   do it=1,dimOB
      read(iunit,rec=it+(ip-1)*dimOB) ints(1:dimOB,1:dimOB)
      print* , "ints(1,1)", ints(1,1)
      do iq=1,dimOA
         do ir=1,dimOB
            do is =1,dimOB
               do iu=1,dimOB
            !     tNb(2) = tNb(2) + A%Occ(it)*intB(ip,ir,it,iq)*ints(ip+(ir-1)*NBas)
                  !tNb(2) = tNb(2) + Atrdm(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,iu)! <- to dobrze ale cos nie dziala
                  tNb(2) = tNb(2) + Atrdmtest(ip,iq)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,iu)
                  !TEST!
                  !tNa(2) = tNa(2) + Btrdm2(it,iu)*Trdm2A(ip,iq,ir,is)*ints(ip,iq)*Sab2(ir,it)
                  !tNb(2) = tNb(2) + Atrdm(iq,ip)*Trdm2B(ir,is,it,iu)*ints(ir,is)*Sab2(iq,it)
               enddo   
            enddo   
         enddo
      enddo

   enddo
enddo
tNb(2) = -2d0*tNb(2)
print*, 'tNb-1',tNb(1)*1000
print*, 'tNb-2',tNb(2)*1000

close(iunit)




!
open(newunit=iunit,file='OOOOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))

! call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)

val  = 0d0
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


!!! TO jeszcze nie gotowe



! VP = sum(tvk)
! exch1tot1 = SAPT%exchs21 - VP + P1*(SAPT%elst1)
! exch1tot2 = SAPT%exchs22 - VP - P1*(SAPT%elst2)


! print*, 'tvk(4)',tvk(4)*1000
print*, ' to sum up'
print*, 'tvk(1)', tvk(1)
print*, 'tvk(2)', tvk(2)
print*, 'tvk(3)', tvk(3)
print*, 'sum tvk', sum(tvk)
print*, 'tNa(1)',TNa(1)
print*, 'tNa(2)',TNa(2)
print*, 'tNb(1)',tNb(1)
print*, 'tNb(2)',tNb(2)
print*, 'P1', P1
print*, 'exhange elst dSRS 1', exch1tot1*1000
print*, 'exhange elst dSRS 2', exch1tot2*1000
print*, 'elst dSRS 1', SAPT%elst1*1000
print*, 'elst dSRS 2', SAPT%elst2*1000




end subroutine e1exch_dSRS



end module sapt_dRS
