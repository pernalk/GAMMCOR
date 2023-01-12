module sapt_dRS
use types
use tran
use sapt_utils
!use read_external

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
integer :: iref,iref2
integer :: iexcited
integer :: ip,ir,iq,is
double precision :: elab, elab2,elab3
double precision,allocatable :: work(:,:)
double precision,allocatable :: Btrdm(:,:)
double precision,allocatable :: Atrdm(:,:)

integer :: i
double precision,allocatable :: Vb(:,:)
double precision,allocatable :: Vbaa(:,:)
double precision,allocatable :: Vbaa2(:,:)
double precision :: ea1,ea2

! set dimensions
iref   =  1
iref2  =  2
iexcited = 1
NBasis = A%NBasis
dimOA  = A%num0+A%num1
dimOB  = B%num0+B%num1

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

call tran4_gen(NBasis,&
     B%num0+B%num1,B%CAONO(iref2,:,:),&
     B%num0+B%num1,B%CAONO(iref2,:,:),&
     A%num0+A%num1,A%CAONO(iref,:,:),&
     A%num0+A%num1,A%CAONO(iref,:,:),&
     'OOOOAABB','AOTWOSORT')

allocate(work(dimOA,dimOA))
! n_p * n_q * v_pq^pq
open(newunit=iunit,file='OOOOAABB',status='old',access='direct',&
     form='unformatted',recl=8*dimOA**2)
elab = 0
do ip=1,dimOB
   ! get all (AA| integrals for a given |BB) record
   read(iunit,rec=ip+(ip-1)*dimOB) work(1:dimOA,1:dimOA)
   do iq=1,dimOA
      elab = elab + A%rdm1(iref,iq)*B%rdm1(iref2,ip)*work(iq,iq)
   enddo
enddo
close(iunit)

! test
print*, 'elab = ',4d0*elab

allocate(Atrdm(NBasis,NBasis))
allocate(Btrdm(NBasis,NBasis))
!call tran2MO(A%trdm(iref,:,:),CMONO(RefState,:,:),CMONO(RefState,:,:), &
!                        Atrdm(iref,:,:),NBasis)
call tran2MO(A%trdm1(iexcited,:,:),A%CMONO(iref,:,:),A%CMONO(iref,:,:), &
                        Atrdm(:,:),NBasis)
call tran2MO(B%trdm1(iexcited,:,:),B%CMONO(iref2,:,:),B%CMONO(iref2,:,:), &
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
               elab2 = elab2 + A%rdm1(iref,iq)*B%rdm1(iref2,ip)*work(iq,iq)
            endif
         enddo
      enddo
   enddo
enddo
close(iunit)

! test2
elab2 = 4d0*elab2
print*, 'elab2 = ',elab2
print *, 'elab3 =',elab3
print *, 'EdRS(1)+= ' , (ea1+ea2+elab2+elab3+SAPT%Vnn)*1d3
print *, 'EdRS(1)-= ' , (ea1+ea2+elab2-elab3+SAPT%Vnn)*1d3

deallocate(work)

end subroutine elst_dRS

end module sapt_dRS
