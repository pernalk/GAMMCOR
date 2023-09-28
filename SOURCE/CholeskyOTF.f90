module choleskyOTF_interface

use print_units
use tran
use basis_sets
use sys_definitions
use Auto2eInterface
use CholeskyOTF, only : TCholeskyVecsOTF
use Cholesky_driver, only : chol_Rkpq_OTF, chol_Rkab_OTF, chol_F, chol_H0_mo

contains

subroutine CholeskyOTF_ao_vecs(CholeskyVecsOTF, &
                          AOBasis,System,XYZPath,BasisSetPath, &
                          SortAngularMomenta, Accuracy)
!
!           Generate Cholesky vectors in AO basis on-the-fly
!
!use arithmetic
!use auto2e
!use Cholesky, only: chol_CoulombMatrix, TCholeskyVecs, &
!                    chol_Rkab_ExternalBinary, chol_MOTransf_TwoStep
!use CholeskyOTF, only: chol_CoulombMatrix_OTF, &
!                       TCholeskyVecsOTF, chol_MOTransf_TwoStep_OTF
!use Cholesky_driver
!use basis_sets
!use sys_definitions
!use chol_definitions

            implicit none

            type(TCholeskyVecsOTF), intent(out) :: CholeskyVecsOTF
            type(TAOBasis), intent(out)         :: AOBasis
            type(TSystem), intent(out)          :: System
            character(*), intent(in)            :: XYZPath
            character(*), intent(in)            :: BasisSetPath
            logical, intent(in)                 :: SortAngularMomenta
            integer, intent(in)                 :: Accuracy

            logical, parameter :: SpherAO = .true.

            ! Initialize the two-electron intergrals library
            !
            call auto2e_init()
            !
            ! Read the XYZ coordinates and atom types
            !
            call sys_Read_XYZ(System, XYZPath)
            !
            call sys_Init(System,SYS_TOTAL)
            !
            ! Read the basis set parameters from an EMSL text file
            ! (GAMESS-US format, no need for any edits, just download it straight from the website)
            !
            call basis_NewAOBasis(AOBasis, System, &
                            BasisSetPath, SpherAO, SortAngularMomenta)
            !
            ! Compute Cholesky vectors in AO basis
            !
            call chol_Rkpq_OTF(CholeskyVecsOTF, AOBasis, Accuracy)

end subroutine CholeskyOTF_ao_vecs

subroutine CholeskyOTF_Fock_MO_v1(F_mo,CholeskyVecsOTF, &
                           AOBasis,System,Monomer, & 
                           Cmat,CSAO,H0in,GammaF,  &
                           MemType,MemVal,NInte1,NBasis, &
                           J_mo,K_mo)
!
!     Generate Fock matrix (NBasis,NBasis) in MO basis
!     from Cholesky OTF vectors
!     optional :: compute J and K matrices in MO basis
!
implicit none

type(TCholeskyVecsOTF), intent(in) :: CholeskyVecsOTF
type(TAOBasis), intent(in)         :: AOBasis
type(TSystem), intent(inout)       :: System

integer,intent(in)          :: Monomer
integer,intent(in)          :: NInte1,NBasis
integer,intent(in)          :: MemType,MemVal
double precision,intent(in) :: H0in(NInte1),GammaF(NInte1)
double precision,intent(in) :: Cmat(NBasis,NBasis),CSAO(NBasis,NBasis)
                               
double precision,intent(out) :: F_mo(NBasis,NBasis)
double precision,optional,intent(out) :: J_mo(NBasis,NBasis)
double precision,optional,intent(out) :: K_mo(NBasis,NBasis)

integer          :: i, j
integer          :: MemMOTransfMB
double precision :: val, h0norm
double precision :: H0tr(NInte1)
double precision :: D_mo(NBasis,NBasis), H0_mo(NBasis,NBasis)
double precision :: work(NBasis,NBasis)
double precision, parameter :: ThreshH0 = 1d-6
integer,external :: IndSym

write(6,'(/,1x,a,i2)') "Construct Fock matrix, monomer", Monomer

! Read whether to put ghost functions
if(Monomer==1) then
   call sys_Init(System,SYS_MONO_A)
elseif(Monomer==2) then
   call sys_Init(System,SYS_MONO_B)
else
   call sys_Init(System,SYS_TOTAL)
endif

! transform H0 (in SAO) to MO
H0tr = H0in
call tran_matTr(H0tr,CSAO,CSAO,NBasis,.true.)
call triang_to_sq2(H0tr,H0_mo,NBasis)
h0norm = norm2(H0_mo)

!write(6,*) 'H0   MO (external)'
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (H0_mo(i,j),i=1,NBasis)
!enddo
!write(LOUT,'()')

! prepare density matrix
H0_mo = 0d0
D_mo  = 0d0
do j=1,NBasis
 do i=1,NBasis
    D_mo(I,J) = 2d0 * GammaF(IndSym(I,J))
 enddo
enddo

!      write(6,*) 'DMAT MO',norm2(D_mo)
!      do j=1,NBasis
!         write(LOUT,'(*(f13.8))') (D_mo(i,j),i=1,NBasis)
!      enddo
!      write(LOUT,'()')

!        write(6,*) 'CAOMO  '
!        do j=1,NBasis
!           write(LOUT,'(*(f13.8))') (Cmat(i,j),i=1,NBasis)
!        enddo
!        write(LOUT,'()')

! set memory for Fock transformation
if(MemType == 2) then       !MB
   MemMOTransfMB = MemVal
elseif(MemType == 3) then   !GB
   MemMOTransfMB = MemVal * 1024_8
endif

if(present(J_mo) .and. present(K_mo)) then
! generate also J_mo and K_mo
   call chol_F(F_mo,H0_mo,D_mo,Cmat,0.5d0,CholeskyVecsOTF, &
           AOBasis,System,ORBITAL_ORDERING_MOLPRO,MemMOTransfMB, &
           J_mo,K_mo)
else
   call chol_F(F_mo,H0_mo,D_mo,Cmat,0.5d0,CholeskyVecsOTF, &
               AOBasis,System,ORBITAL_ORDERING_MOLPRO,MemMOTransfMB)
endif

!print*, 'FockF w bazie MO -- new'
!print*, 'Fock = ',norm2(F_mo)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (F_mo(i,j),i=1,NBasis)
!enddo

!print*, 'J mat w bazie MO -- new'
!print*, 'J mat= ',norm2(J_mo)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (J_mo(i,j),i=1,NBasis)
!enddo
!
!print*, 'K mat w bazie MO -- new'
!print*, 'K mat= ',norm2(K_mo)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (K_mo(i,j),i=1,NBasis)
!enddo

!call CholeskyOTF_H0_test(H0_mo,,Nbasis)

end subroutine CholeskyOTF_Fock_MO_v1

subroutine CholeskyOTF_Fock_MO_v2(F_mo,CholeskyVecsOTF, &
                             AOBasis,System,Monomer,Source,&
                             Cmat,CSAO,H0in,GammaF,&
                             MemType,MemVal,NInte1,NBasis,&
                             IH0Test,J_mo,K_mo)
!
!     Generate Fock matrix (NBasis,NBasis) in MO basis
!     from Cholesky OTF vectors
!     optional :: compute J and K matrices in MO basis
!
!     CAREFUL! Involves one 3-ind AO --> MO transformation
!              with FF indices
!

implicit none

type(TCholeskyVecsOTF), intent(in) :: CholeskyVecsOTF
type(TAOBasis), intent(in)         :: AOBasis
type(TSystem), intent(inout)       :: System

integer,intent(in)          :: Monomer
character(6),intent(in)     :: Source
integer,intent(in)          :: NInte1,NBasis
integer,intent(in)          :: MemType,MemVal
integer,intent(in)          :: IH0Test
double precision,intent(in) :: H0in(NInte1),GammaF(NInte1)
double precision,intent(in) :: Cmat(NBasis,NBasis),CSAO(NBasis,NBasis)

double precision,intent(out)          :: F_mo(NBasis,NBasis)
double precision,optional,intent(out) :: J_mo(NBasis,NBasis)
double precision,optional,intent(out) :: K_mo(NBasis,NBasis)

integer          :: i, j
integer          :: NCholesky
integer          :: ORBITAL_ORDERING
integer          :: MemMOTransfMB
double precision :: val
double precision :: H0tr(NInte1)
double precision :: D_mo(NBasis,NBasis), H0_mo(NBasis,NBasis)
double precision :: H0_int(NBasis,NBasis)
double precision :: Jtmp(NBasis,NBasis), Ktmp(NBasis,NBasis)
double precision :: work(NBasis,NBasis)
double precision, allocatable :: ints(:), matFFMO(:,:)
double precision, parameter :: ThreshH0 = 1.0d-6

integer,external :: IndSym
double precision,external :: ddot

write(6,'(/,1x,a,i2)') "Construct Fock ver 2, monomer", Monomer
if(present(J_mo)) print*, 'J_mo present?'

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs

! set orbital ordering
if(trim(Source)=='MOLPRO') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
elseif(trim(Source)=='ORCA  ') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_ORCA
elseif(trim(Source)=='DALTON') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
endif

! Read whether to put ghost functions
if(Monomer==1) then
   call sys_Init(System,SYS_MONO_A)
elseif(Monomer==2) then
   call sys_Init(System,SYS_MONO_B)
else
   call sys_Init(System,SYS_TOTAL)
endif

! transform H0 (in SAO) to MO
H0tr = H0in
if(trim(Source)=='MOLPRO') then
   call tran_matTr(H0tr,CSAO,CSAO,NBasis,.true.)
   call triang_to_sq2(H0tr,H0_mo,NBasis)
elseif(trim(Source)=='ORCA  ') then
   call triang_to_sq2(H0in,H0_mo,NBasis)
endif

! prepare density matrix
D_mo  = 0d0
do J=1,NBasis
   do I=1,NBasis
      D_mo(I,J) = 2d0 * GammaF(IndSym(I,J))
   enddo
enddo

! set memory for Fock transformation
if(MemType == 2) then       !MB
   MemMOTransfMB = MemVal
elseif(MemType == 3) then   !GB
   MemMOTransfMB = MemVal * 1024_8
endif

! obtain H0 and check if they match
call chol_H0_mo(H0_int,Cmat,AOBasis,System,ORBITAL_ORDERING)

if(IH0Test==1) then
   call CholeskyOTF_H0_test(H0_int,H0_mo,NBasis)
elseif(IH0Test==0) then
   write(6,'(/,1x,"Skipping H0 Test...")')
endif

!transform Cholesky vecs to MO
allocate(MatFFMO(NCholesky,NBasis**2),ints(NBasis**2))

call chol_Rkab_OTF(MatFFMO,Cmat,1,NBasis,Cmat,1,NBasis, &
                   MemMOTransfMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

! construct J and K in MO
ints = 0d0
Jtmp = 0d0
Ktmp = 0d0
do i=1,NCholesky
   ints(:) = MatFFMO(i,:)
   val = ddot(NBasis**2,ints,1,D_mo,1)
   call daxpy(NBasis**2,2d0*val,ints,1,Jtmp,1)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,ints,NBasis, &
              D_mo,NBasis,0d0,work,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,-1d0,work,NBasis, &
              ints,NBasis,1d0,Ktmp,NBasis)
enddo

F_mo = Jtmp + Ktmp
F_mo = H0_mo + 0.5d0*F_mo

! return J and K in MO upon request
if(present(J_mo).and.present(K_mo)) then
   J_mo = 0.5d0*Jtmp
   K_mo = 0.5d0*Ktmp
endif

!print*, 'Fock in MO basis'
!print*, 'Fock = ',norm2(F_mo)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (F_mo(i,j),i=1,NBasis)
!enddo
!
!print*, 'J_mo = ',norm2(Jtmp)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (Jtmp(i,j),i=1,NBasis)
!enddo
!
!print*, 'K_mo = ',norm2(Ktmp)
!do j=1,NBasis
!   write(LOUT,'(*(f13.8))') (Ktmp(i,j),i=1,NBasis)
!enddo

deallocate(MatFFMO,ints)

end subroutine CholeskyOTF_Fock_MO_v2

subroutine CholeskyOTF_H0_test(H0_int,H0_mo,NBasis)
!
! check if H0 from Orca/Molpro match the ones
! generated by gammcor-cholesky module:
! if not, basis sets or geometry are incompatible!
!
implicit none

integer,intent(in)          :: NBasis
double precision,intent(in) :: H0_int(NBasis,NBasis)
double precision,intent(in) :: H0_mo(NBasis,NBasis)

integer :: i,j
double precision            :: val1,val2
double precision, parameter :: ThreshH0 = 1d-6

val1 = norm2(H0_int)
val2 = norm2(H0_mo)

write(6,'(1x,a,f12.6,/)') "H0 (internal) vs. H0 (external) = ",abs(val1)-abs(val2)

if(abs(val1)-abs(val2).gt.ThreshH0) then
  print*, 'Difference in H0 norms = '
  print*, abs(val1)-abs(val2)
  print*, 'Threshold = ',ThreshH0
  print*, 'Check for errors in geometry / basis set?'

  write(6,*) 'H0/MO (internal)'
  do j=1,NBasis
     write(LOUT,'(*(f13.8))') (H0_int(i,j),i=1,NBasis)
  enddo

  print*, 'H0/MO (external) '
  do j=1,NBasis
     write(LOUT,'(*(f13.8))') (H0_mo(i,j),i=1,NBasis)
  enddo

  stop
endif

end subroutine CholeskyOTF_H0_test

end module choleskyOTF_interface
