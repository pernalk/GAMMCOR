module sapt_dRS
use types
use tran
use sapt_utils
use read_external

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
print*, '1-TRDM in NO'
      do i=1,10
         write(lout,'(*(f12.6))') (Atrdm(i,j),j=1,10)
      enddo
      print*,'Atrdm1(1,2)',Atrdm(1,2)
      print*,'Atrdm1(1,3)',Atrdm(1,3)
      




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
integer :: ip,ir,iq,is,iu,it,iv,iw
double precision, allocatable :: trdm2A(:,:,:,:)
double precision, allocatable :: trdm2B(:,:,:,:)
double precision :: tvk(5)
double precision, allocatable :: Atrdm(:,:)
double precision, allocatable :: Btrdm(:,:),Btrdm2(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Sab2(:,:),Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:)
double precision, allocatable :: work(:,:)
double precision, allocatable :: work2A(:,:,:,:), work2B(:,:,:,:)
double precision, allocatable :: intA(:,:,:,:), intB(:,:,:,:)
double precision, allocatable :: ints(:,:)
double precision, allocatable :: tmpAB(:,:,:,:)
double precision :: val, P1, exch1tot1, exch1tot2, VP 
print *, "running exch dSRS"
print *,"exch w DRS ",SAPT%exchs21


print*,"2TRDM"
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
         Vaab(NBasis,NBasis),Vbba(NBasis,NBasis))


! get V_ne in atomic orbs
call get_one_mat('V',Va,A%Monomer,NBasis)
call get_one_mat('V',Vb,B%Monomer,NBasis)



!call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
! not sure all are needed
call tran2MO(Va,B%CAONO(iref2B,:,:),B%CAONO(iref2B,:,:),Vabb,NBasis)
call tran2MO(Vb,A%CAONO(iref,:,:),A%CAONO(iref,:,:),Vbaa,NBasis)
call tran2MO(Va,A%CAONO(iref,:,:),B%CAONO(iref2B,:,:),Vaab,NBasis)
call tran2MO(Vb,B%CAONO(iref2B,:,:),A%CAONO(iref,:,:),Vbba,NBasis)
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
            P1 = P1 +Atrdm(ip,iq)*Btrdm(is,ir)*Sab(iq,ir)*Sab(ip,is) !poprawione
 !!           tvk(1) = tvk(1) + Atrdm(ip,iq)*Btrdm(is,ir)*Vbab(ip,is)*Sab(iq,ir) !poprawione
         Enddo
      enddo
   Enddo
enddo
tvk(1) = -4d0*tvk(1)
P1= -2d0*P1
!print*, 'tvk(1) z dSRS',tvk(1)*1000




!WORK IN PROGRESS
allocate(work(dimOA,dimOB))
! tvk(2)=- 4sum(trdm(A->A*)_{pq}trdm(B*->B)_{rs}(ps|qr)
open(newunit=iunit,file='OOOOABBA',status='old',access='direct',&
     form='unformatted',recl=8*dimOA*dimOB)
do ip=1,dimOB!ir=1,dimOB
   do iq=1,dimOA!iq=1,dimOA
   ! get all (AB| integrals for a given |BA) record
      read(iunit,rec=ip+(iq-1)*dimOB) work(1:dimOA,1:dimOB)
      do ir = 1,dimOB!is=dimOB
         do is = 1,dimOA!ip=dimOA
            tvk(2) = tvk(2) + Atrdm(is,iq)*Btrdm(ir,ip)*work(is,ir) !poprawione
         enddo
      enddo
   enddo
enddo
close(iunit)
tvk(2) = -4d0*tvk(2)
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
 call tran4_gen(NBasis,&
               A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
               B%num0+B%num1,B%CAONO(irefB,1:NBasis,1:(B%num0+B%num1)),&
               A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
               A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
               'OOOOAAAB','AOTWOSORT')
 !
 !calculating TRDM2*TRDM1*V_{ee}*S
 !
 !
allocate(intA(dimOA,dimOA,dimOA,dimOB),&
         intB(dimOB,dimOB,dimOA,dimOB))

 
call tran2MO(S,A%CAONO(iref,:,:),B%CAONO(irefB,:,:),Sab2,NBasis)

call tran2MO(B%trdm1(iexcited2,:,:),B%CMONO(irefB,:,:),B%CMONO(irefB,:,:), &
                             Btrdm2(:,:),NBasis)
 

call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%trdm24(iref2,:,:,:,:),dimOA**3,Sab2,NBasis,0d0,intA,dimOA**3)
!call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2Bval,dimOB**3,Sab,NBas,0d0,intB,dimOB**3)

! careful! intB(B,B,A,B)
do is=1,dimOB
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%trdm24(iref2B,:,:,:,is),dimOB**2,Sab2,NBasis,0d0,intB(:,:,:,is),dimOB**2)
enddo
!!!! WORK IN PROGRESS
print *,"dimOA",dimOA
allocate(ints(dimOA,dimOA))

  
open(newunit=iunit,file='OOOOAAAB',status='OLD', access='DIRECT', &
     form='unformatted' ,recl=8*dimOA*dimOB)
tvk(3)=0d0
! one loop over integrals
ints = 0d0
do it=1,dimOB
    do iq=1,dimOA
       read(iunit,rec=iq+(it-1)*dimOA) ints(1:dimOA,1:dimOA)
!       print* , "ints", ints
      do iu =1,dimOB
      do ir=1,dimOA
         do ip=1,dimOA
            !tNa(2) = tNa(2) + B%Occ(it)*intA(ip,ir,iq,it)*ints(ip+(ir-1)*NBas)
            tvk(3) = tvk(3) + Btrdm2(it,iu)*intA(ip,ir,iq,iu)*ints(ip,ir)
!            print *,tvk(3),Btrdm2(it,iu),intA(ip,ir,iq,iu),ints(ip,ir)
         enddo
      enddo
   enddo
   enddo
enddo
tvk(3) = -4d0*tvk(3) !poprawione
print*, 'tvk(3)',tvk(3)*1000

close(iunit)

!!! END OF WORK IN PROGRESS

!!call delfile('OOOOAABB') !! To nie wiem jak działa
call tran4_gen(NBasis,&
     B%num0+B%num1,B%CAONO(irefB,1:NBasis,1:(B%num0+B%num1)),&
     B%num0+B%num1,B%CAONO(irefB,1:NBasis,1:(B%num0+B%num1)),&
     A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
     A%num0+A%num1,A%CAONO(iref,1:NBasis,1:(A%num0+A%num1)),&
     'OOOOAABB','AOTWOSORT')


!
open(newunit=iunit,file='OOOOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))

call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)

val  = 0d0
!do ir=1,dimOA
!   do ip=1,dimOA
    ! read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOB,1:dimOB)
     !val = val + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))
     !intA(ip,ir,iq,it)*intB(it,iu,is,iq)*ints(it,is)
!   enddo
!enddo

ints = 0d0
do ir=1,dimOB
   do ip=1,dimOB
     read(iunit,rec=ip+(ir-1)*dimOB) ints(1:dimOA,1:dimOA)
     val = val + sum(ints(1:dimOA,1:dimOA)*tmpAB(1:dimOA,1:dimOA,ip,ir))
   enddo
enddo
close(iunit)
tvk(4) = -2d0*val!poprawione
do iu=1,dimOB
   do it=1,dimOB
   do iq=1,dimOA
      do ir=1,dimOA
         do ip=1,dimOA
            !   tNa(1) = tNa(1) + B%Occ(it)*intA(ip,ir,iq,it)*Sab(iq,it)*Vbaa(ip,ir)
            tvk(5) = tvk(5) + Btrdm2(it,iu)*intA(ip,ir,iq,it)*Sab(iq,iu)*Vbaa(ip,ir)
         enddo
      enddo
   enddo
enddo
enddo
tvk(5)= -4d0*tvk(5)!poprawione







VP = sum(tvk)
exch1tot1 = SAPT%exchs21 - VP + P1*(SAPT%elst1)
exch1tot2 = SAPT%exchs22 - VP - P1*(SAPT%elst2)


print*, 'tvk(4)',tvk(4)*1000
print*, ' to sum up'
print*, 'tvk(1)', tvk(1)
print*, 'tvk(2)', tvk(2)
print*, 'tvk(3)', tvk(3)
print*, 'tvk(4)', tvk(4)
print*, 'tvk(5)', tvk(5)
print*, 'sum tvk', sum(tvk)
print*, 'P1', P1
print*, 'exhange elst dSRS 1', exch1tot1*1000
print*, 'exhange elst dSRS 2', exch1tot2*1000
print*, 'elst dSRS 1', SAPT%elst1*1000
print*, 'elst dSRS 2', SAPT%elst2*1000




end subroutine e1exch_dSRS



end module sapt_dRS
