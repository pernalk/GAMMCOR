module sapt_dRS
use types
use tran
!use sapt_utils
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
integer :: iref
integer :: ip,ir,iq,is
double precision :: elab
double precision,allocatable :: work(:,:)

! set dimensions
iref   = 1
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
call tran4_gen(NBasis,&
               B%num0+B%num1,B%CMO,&
               B%num0+B%num1,B%CMO,&
               A%num0+A%num1,A%CMO,&
               A%num0+A%num1,A%CMO,&
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
      elab = elab + A%rdm1(iref,iq)*B%rdm1(iref,ip)*work(iq,iq)
   enddo
enddo
close(iunit)

! test
print*, 'elab = ',4d0*elab

deallocate(work)
print*, 'nothing yet in elst_dRS; quitting...'

end subroutine elst_dRS

end module sapt_dRS
