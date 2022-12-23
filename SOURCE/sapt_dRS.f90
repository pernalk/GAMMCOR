module sapt_dRS
use types
!use tran
!use sapt_utils
!use read_external

implicit none

contains

subroutine elst_dRS(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

print*, 'nothing yet in elst_dRS; quitting...'

end subroutine elst_dRS

end module sapt_dRS
