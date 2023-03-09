#define RDMCORR_DEBUG -1

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

integer :: i,j,ij
integer :: iunit
integer :: NInte1,HlpDim

double precision             :: val
double precision             :: EVal(NBasis),Eps(NBasis,NBasis)
double precision             :: CNONO(NBasis,NBasis)
double precision             :: URe(NBasis,NBasis),AuxMat(NBasis,NBasis)
double precision,allocatable :: XOne(:),work(:)

character(8) :: label
character(:),allocatable     :: onefile,twojfile,twokfile

HlpDim = max(NBasis**2,3*NBasis)
NInte1 = NBasis*(NBasis+1)/2

! save old NDimX
Mon%NDimX0 = Mon%NDimX

if(Mon%Monomer==1) then
   onefile  = 'ONEEL_A'
   twojfile = 'FFOOAA'
   twokfile = 'FOFOAA'
elseif(Mon%Monomer==2) then
   onefile  = 'ONEEL_B'
   twojfile = 'FFOOBB'
   twokfile = 'FOFOBB'
endif

! prepare RDM2
if(Flags%ICASSCF==1) then
   if(Mon%Monomer==1) call system('cp rdm2_A.dat rdm2.dat')
   if(Mon%Monomer==2) call system('cp rdm2_B.dat rdm2.dat')
endif

! get H0 matrix
allocate(XOne(NInte1))

open(newunit=iunit,file=onefile,access='sequential',&
     form='unformatted',status='old')

read(iunit)
read(iunit)
read(iunit) label, AuxMat
if(label=='ONEHAMIL') then
   call dgemm('T','N',NBasis,NBasis,NBasis,1d0,Mon%CMO,NBasis,AuxMat,NBasis,0d0,Eps,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,Eps,NBasis,Mon%CMO,NBasis,0d0,AuxMat,NBasis)
   ij = 0
   do j=1,NBasis
      do i=1,j
         ij = ij + 1
         XOne(ij) = AuxMat(i,j)
      enddo
   enddo
else
   write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
   stop
endif

close(iunit)

print*, 'XONE : ',norm2(XOne)

! get correlation contribution to 1-rdm
allocate(work(HlpDim))
allocate(Mon%rdm1c(NBasis,NBasis))

URe = 0d0
do i=1,NBasis
   URe(i,i) = 1d0
enddo
AuxMat = transpose(Mon%CMO)

! return MP2-like rdm1c = CAS + corr
call MP2RDM_FOFO(0d0,Eps,Mon%Occ,URe,AuxMat,XOne,&
                 Mon%IndN,Mon%IndX,Mon%IndAux,Mon%IGem,  &
                 Mon%NAct,Mon%INAct,Mon%NDimX,Mon%NDimX,NBasis,NInte1, &
                 twojfile,twokfile,Mon%ThrVirt,Mon%NVZero,Mon%IPrint,  &
                 Mon%rdm1c)

print*, 'RDMcorr ',norm2(Mon%rdm1c)

! test
!Mon%rdm1c = 0d0
!do i=1,NBasis
!   Mon%rdm1c(i,i) = Mon%Occ(i)
!enddo

CNONO = 0d0
CNONO = Mon%rdm1c

! diagonalize the new 1rdm 
call Diag8(CNONO,NBasis,NBasis,Eval,work)

do i=1,NBasis
   if (EVal(i) < 0d0) then
      print*, 'Negative occupation = ', i, EVal(i)
   endif
enddo

#if RDMCORR_DEBUG > 5
   val = 0d0
   if (Mon%IPrint > 10) write(LOUT,'(2x,"MP2",3x,"Unsorted Occupancy")')
   do i=NBasis,1,-1
      write(LOUT,'(X,I3,E16.6,I6)') i,EVal(i)
      val = val + EVal(i)
   enddo
   write(LOUT,'(/,1x,"Sum of MP2 Occupancies: ",F5.2,/)') val
#endif

call SortOcc(EVal,CNONO,NBasis)

val = 0d0
write(LOUT,'(2x,"MP2",3x,"Sorted Occupancy")')
do i=1,NBasis
   write(LOUT,'(X,I3,E16.6,I6)') i,EVal(i)
   val = val + EVal(i)
enddo
write(LOUT,'(/,1x,"Sum of MP2 Occupancies: ",F5.2,/)') val

! canonicalize
print*, 'WARNING! Add canonicalization for 2nd order!'

! save new orbital coefficients 
! M%CMO = C(AO,NO)
! C(AO,NO).C(NOMP2,NO)^T
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,Mon%CMO,NBasis,CNONO,NBasis,0d0,AuxMat,NBasis)

#if RDMCORR_DEBUG > 5
   print*, 'CAONO:before',norm2(Mon%CMO)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (Mon%CMO(i,j),i=1,NBasis)
   enddo
   print*, 'CAONO:after',norm2(AuxMat)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (AuxMat(i,j),i=1,NBasis)
   enddo
#endif

! save new occupation numbers
print*, 'RDMCORR: Occ-1', norm2(mon%Occ)
print*, 'RDMCORR: EVal ', norm2(EVal)
Mon%Occ = EVal
print*, 'RDMCORR: Occ-2', norm2(mon%Occ)

! save new C(AO,NO) coefficients
Mon%CMO = AuxMat

deallocate(XOne)

end subroutine sapt_rdm_corr

end module rdmcorr

