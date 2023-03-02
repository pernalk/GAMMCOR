module rdmcorr

use print_units
use types
use abfofo

contains

subroutine sapt_rdm_corr(Mon,Flags,NAO,NBasis)
implicit none

type(FlagsData)     :: Flags
type(SystemBlock)   :: Mon
integer, intent(in) :: NAO, NBasis

integer :: i,j
integer :: HlpDim

double precision             :: EVal(NBasis)
double precision             :: CNONO(NBasis,NBasis)
double precision             :: AuxMat(NBasis,NBasis)
double precision,allocatable :: work(:)

HlpDim = max(NBasis**2,3*NBasis)

! get correlation contribution to 1-rdm
allocate(work(HlpDim))
allocate(Mon%rdm1c(NBasis,NBasis))

! call ...

Mon%rdm1c = 0d0
do i=1,NBasis
   Mon%rdm1c(i,i) = Mon%Occ(i) 
enddo

CNONO = 0d0
CNONO = Mon%rdm1c

! diagonalize the new 1rdm 
call Diag8(CNONO,NBasis,NBasis,Eval,work)

do i=1,NBasis
   if (EVal(i) < 0d0) then
      print*, 'Negative occupation = ', i, EVal(i)
   endif
enddo
call SortOcc(EVal,CNONO,NBasis)

! canonicalize
print*, 'WARNING! Add canonicalization for 2nd order!'

! save new orbital coefficients 
! M%CMO = C(AO,NO)
! C(AO,NO).C(NOMP2,NO)^T
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,Mon%CMO,NBasis,CNONO,NBasis,0d0,AuxMat,NBasis)

print*, 'CAONO:before'
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (Mon%CMO(i,j),i=1,NBasis)
enddo
print*, 'CAONO:after'
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (AuxMat(i,j),i=1,NBasis)
enddo

! save new occupation numbers
Mon%Occ = EVal

! save new C(AO,NO) coefficients
Mon%CMO = AuxMat

end subroutine sapt_rdm_corr

end module rdmcorr

