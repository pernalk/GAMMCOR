#define CHOLOTF_DEBUG -1

module CholeskyOTF_interface

use print_units
use tran
use basis_sets
use sys_definitions
use Auto2eInterface
use Cholesky_Gammcor, only : TCholeskyVecsOTF, &
                             chol_gammcor_Rkpq, chol_gammcor_Rkab
use THC_Gammcor, only : thc_gammcor_XZ, thc_gammcor_Xga, thc_gammcor_Rkab_2
use OneElectronInts_Gammcor, only : ints1e_gammcor_H0_mo, ints1e_gammcor_H0_extao
use OneElectronInts, only : ints1e_S

contains

subroutine CholeskyOTF_ao_vecs(CholeskyVecsOTF, &
                          AOBasis,System,Units,XYZPath,BasisSetPath, &
                          SortAngularMomenta, Accuracy, Omega)
!
!           Generate Cholesky vectors in AO basis on-the-fly
!
!use arithmetic
!use auto2e
!use Cholesky, only: chol_CoulombMatrix, TCholeskyVecs, &
!                    chol_Rkab_ExternalBinary, chol_MOTransf_TwoStep
!use CholeskyOTF, only: chol_CoulombMatrix_OTF, &
!use basis_sets
!use sys_definitions
!use chol_definitions

            implicit none

            type(TCholeskyVecsOTF), intent(out)    :: CholeskyVecsOTF
            type(TAOBasis), intent(out)            :: AOBasis
            type(TSystem), intent(out)             :: System
            character(*), intent(in)               :: XYZPath
            character(*), intent(in)               :: BasisSetPath
            logical, intent(in)                    :: SortAngularMomenta
            integer, intent(in)                    :: Units
            integer, intent(in)                    :: Accuracy
            double precision, optional, intent(in) :: Omega

            logical, parameter :: SpherAO = .true.

            ! Initialize the two-electron intergrals library
            !
            call auto2e_init()
            !
            ! Read the XYZ coordinates and atom types
            !
            call sys_Read_XYZ(System, XYZPath,Units)
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
            if (present(Omega)) then
               print*, 'calling chol_gammcor_Rkpq with omega = ',Omega
               call chol_gammcor_Rkpq(CholeskyVecsOTF, AOBasis, Accuracy, Omega)
            else
               call chol_gammcor_Rkpq(CholeskyVecsOTF, AOBasis, Accuracy)
            endif

end subroutine CholeskyOTF_ao_vecs

subroutine THC_ao_vecs(Xgp, Zgk, &
                       AOBasis,System,Units,XYZPath,BasisSetPath, &
                       SortAngularMomenta, Accuracy)

            implicit none

            double precision, dimension(:, :), allocatable, intent(out) :: Xgp
            double precision, dimension(:, :), allocatable, intent(out) :: Zgk
            type(TAOBasis), intent(out)    :: AOBasis
            type(TSystem), intent(out)     :: System
            character(*), intent(in)       :: XYZPath
            character(*), intent(in)       :: BasisSetPath
            logical, intent(in)            :: SortAngularMomenta
            integer, intent(in)            :: Units
            integer, intent(in)            :: Accuracy

            logical, parameter :: SpherAO = .true.

            ! Initialize the two-electron intergrals library
            !
            call auto2e_init()
            !
            ! Read the XYZ coordinates and atom types
            !
            call sys_Read_XYZ(System, XYZPath,Units)
            !
            call sys_Init(System,SYS_TOTAL)
            !
            ! Read the basis set parameters from an EMSL text file
            ! (GAMESS-US format, no need for any edits, just download it straight from the website)
            !
            call basis_NewAOBasis(AOBasis, System, &
                            BasisSetPath, SpherAO, SortAngularMomenta)
            !
            ! Compute THC vectors in AO basis
            !
            call thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, Accuracy)

end subroutine THC_ao_vecs

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
   stop "chol_F not available in gammcor-integrals!"
!! generate also J_mo and K_mo
!   call chol_F(F_mo,H0_mo,D_mo,Cmat,0.5d0,CholeskyVecsOTF, &
!           AOBasis,System,ORBITAL_ORDERING_MOLPRO,MemMOTransfMB, &
!           J_mo,K_mo)
!else
!   call chol_F(F_mo,H0_mo,D_mo,Cmat,0.5d0,CholeskyVecsOTF, &
!               AOBasis,System,ORBITAL_ORDERING_MOLPRO,MemMOTransfMB)
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
                             Xgp,Zgk,NGridTHC,NCholeskyTHC,&
                             MemType,MemVal,NInte1,NBasis,&
                             IH0Test,J_mo,K_mo)
!
!     Generate Fock matrix (NBasis,NBasis) in MO basis
!     from Cholesky OTF vectors in AO basis
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
integer,intent(in)          :: NGridTHC,NCholeskyTHC
integer,intent(in)          :: MemType,MemVal
integer,intent(in)          :: IH0Test
double precision,intent(in) :: GammaF(NInte1),H0in(NInte1)
double precision,intent(in) :: Cmat(NBasis,NBasis),CSAO(NBasis,NBasis)
double precision,intent(in) :: Xgp(NGridTHC,NBasis),Zgk(NGridTHC,NCholeskyTHC)

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
double precision, allocatable :: Xga(:,:)
double precision, parameter :: ThreshH0 = 1.0d-6

integer,external :: IndSym
double precision,external :: ddot

write(6,'(/,1x,a,i2)') "Construct Fock ver 2, monomer", Monomer
!if(present(J_mo)) print*, 'J_mo present?'

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs

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

!! test 1-el MO
!block
!   double precision :: trace
!   trace=0d0
!   ! set inact
!   do j=1,12
!   do i=1,12
!      D_mo(i,j)=0d0
!   enddo
!   enddo
!
!   do j=1,NBasis
!      val=0d0
!      do i=1,NBasis
!         val = val + D_mo(i,j)*H0_mo(j,i)
!      enddo
!      trace = trace + val
!   enddo
!   print*, 'One-electron active enenrgy (MO)', trace
!end block

! set memory for Fock transformation
if(MemType == 2) then       !MB
   MemMOTransfMB = MemVal
elseif(MemType == 3) then   !GB
   MemMOTransfMB = MemVal * 1024_8
endif

!block
!double precision :: H0ao(Nbasis,NBasis)
!
!call ints1e_gammcor_H0_extao(H0ao,AOBasis,System,Orbital_Ordering)
!print*, 'H0 ao symmetry'
!do i=1,NBasis
!   write(LOUT,'(*(f13.8))') (H0ao(i,j),j=1,NBasis)
!enddo
!
!end block

! obtain H0 and check if they match
call ints1e_gammcor_H0_mo(H0_int,Cmat,AOBasis,System,ORBITAL_ORDERING)

! mh : test for DMRG-in-DFT emebedding
!block
!   integer :: col_map(NBasis)
!   real(8) :: sgn_map(NBasis)
!   real(8) :: work(NBasis,NBasis)
!   ! reorder MOs based on H0_int/H0_mo
!   call reorder_by_first_row(H0_int,H0_mo,"A",NBasis,work,col_map,sgn_map,1d-4)
!   print*, 'reorder C(AO,MO*) and adjust phase...'
!   do i=1,NBasis
!      print*, i,col_map(i),sgn_map(i)
!      work(:,col_map(i)) = sgn_map(i)*Cmat(:,i)
!   enddo
!   Cmat = work
!end block

if(IH0Test==1) then
#if CHOLOTF_DEBUG > 4
   print*, 'Cmat = ',norm2(Cmat)
   do i=1,NBasis
      write(LOUT,'(*(f13.8))') (Cmat(i,j),j=1,NBasis)
   enddo
   print*,'H0_mo  ',norm2(H0_mo)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (H0_mo(i,j),i=1,NBasis)
   enddo
   print*,'H0_int',norm2(H0_int)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (H0_int(i,j),i=1,NBasis)
   enddo
#endif
   call CholeskyOTF_H0_test(H0_int,H0_mo,NBasis)
elseif(IH0Test==0) then
   write(6,'(1x,"Skipping H0 Test: use internal H0")')

   call sq_to_triang2(H0_int,H0in,NBasis)
   H0_mo = H0_int
elseif(IH0Test==2) then
   write(6,'(1x,"Skipping H0 Test: use external H0")')
endif


if(NGridTHC .gt. 1) then
   !transform Cholesky or THC vecs to MO
   print*, 'NGridTHC     = ',NGridTHC
   print*, 'NCholeskyTHC = ',NCholeskyTHC

   Print*, 'THC in FockOTF: MO transform Xgp to Xga'
   NCholesky = NCholeskyTHC
   allocate(MatFFMO(NCholesky,NBasis**2),Xga(NGridTHC,NBasis))
   Call thc_gammcor_Xga(Xga,Xgp,Cmat,AOBasis,ORBITAL_ORDERING)
   Print*, 'THC in FockOTF: assemble Cholesky vecs'
   Call thc_gammcor_Rkab_2(MatFFMO,Xga,Xga,Zgk,NBasis,NBasis,NCholeskyTHC,NGridTHC)

else

   allocate(MatFFMO(NCholesky,NBasis**2))
   call chol_gammcor_Rkab(MatFFMO,Cmat,1,NBasis,Cmat,1,NBasis, &
                      MemMOTransfMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING)

endif

allocate(ints(NBasis**2))
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
elseif(present(J_mo)) then
   ! return only J in MO
   J_mo = 0.5d0*Jtmp
endif

! test

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

subroutine CholeskyOTF_H0_test0(AOBasis,System,Monomer,Source,Cmat,H0in,NINte1,NBasis)

implicit none

type(TAOBasis), intent(in)   :: AOBasis
type(TSystem), intent(inout) :: System

integer,intent(in)      :: NInte1,NBasis
integer,intent(in)      :: Monomer
character(6),intent(in) :: Source

double precision,intent(in) :: H0in(NInte1)
double precision,intent(in) :: Cmat(NBasis,NBasis)

integer :: i,j
integer          :: ORBITAL_ORDERING
double precision :: H0_int(NBasis,NBasis),H0_mo(NBasis,NBasis)

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


if(trim(Source)=='DALTON') then
  call triang_to_sq2(H0in,H0_mo,NBasis)
else 
  stop "tran H0 in CholeskyOTF"
endif

#if CHOLOTF_DEBUG > 4
print*, 'H0_mo = ',norm2(H0_mo)
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (H0_mo(i,j),i=1,NBasis)
enddo
#endif

print*, 'ORBITAL_ORDERING =',ORBITAL_ORDERING
print*, 'Cmat =',norm2(Cmat)

call ints1e_gammcor_H0_mo(H0_int,Cmat,AOBasis,System,ORBITAL_ORDERING)
#if CHOLOTF_DEBUG > 4
!print*, 'Cmat = ',norm2(Cmat)
!do i=1,NBasis
!   write(LOUT,'(*(f13.8))') (Cmat(i,j),j=1,NBasis)
!enddo
print*, 'H0_int = ',norm2(H0_int)
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (H0_int(i,j),i=1,NBasis)
enddo
print*, 'int-mo = '
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (H0_int(i,j)-H0_mo(i,j),i=1,NBasis)
enddo

#endif

call CholeskyOTF_H0_test(H0_int,H0_mo,NBasis)

end subroutine CholeskyOTF_H0_test0

subroutine CholeskyOTF_Jmat_MO(J_mo,CholeskyVecsOTF, &
                             AOBasis,System,Monomer,Source,&
                             Cmat,CSAO,H0in,GammaF,&
                             MemType,MemVal,NInte1,NBasis,&
                             IH0Test)
!
!     Generate Coulomb matrix (NBasis,NBasis) in MO basis
!     from Cholesky OTF vectors
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

double precision,intent(out)          :: J_mo(NBasis,NBasis)

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

write(6,'(/,1x,a,i2)') "Construct Jmat from OTF Cholesky vectors"

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs

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
elseif(trim(Source)=='ORCA  ' .or.  &
       trim(Source)=='DALTON' ) then
   call triang_to_sq2(H0in,H0_mo,NBasis)
endif

! prepare density matrix
D_mo  = 0d0
do J=1,NBasis
   do I=1,NBasis
      D_mo(I,J) = 2d0 * GammaF(IndSym(I,J))
   enddo
enddo
#if CHOLOTF_DEBUG > 4
print*, 'Dmat in MO:'
do j=1,NBasis
   write(LOUT,'(*(f13.8))') (D_mo(i,j),i=1,NBasis)
enddo
#endif

! set memory for Fock transformation
if(MemType == 2) then       !MB
   MemMOTransfMB = MemVal
elseif(MemType == 3) then   !GB
   MemMOTransfMB = MemVal * 1024_8
endif

! obtain H0 and check if they match
call ints1e_gammcor_H0_mo(H0_int,Cmat,AOBasis,System,ORBITAL_ORDERING)

if(IH0Test==1) then
   call CholeskyOTF_H0_test(H0_int,H0_mo,NBasis)
elseif(IH0Test==0) then
   write(6,'(1x,"Skipping H0 Test: use internal H0")')
   call sq_to_triang2(H0_int,H0in,NBasis)
   H0_mo = H0_int
elseif(IH0Test==2) then
   write(6,'(1x,"Skipping H0 Test: use external H0")')
endif

!transform Cholesky vecs to MO
allocate(MatFFMO(NCholesky,NBasis**2),ints(NBasis**2))

call chol_gammcor_Rkab(MatFFMO,Cmat,1,NBasis,Cmat,1,NBasis, &
                   MemMOTransfMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

! construct J in MO
ints = 0d0
Jtmp = 0d0
do i=1,NCholesky
   ints(:) = MatFFMO(i,:)
   val = ddot(NBasis**2,ints,1,D_mo,1)
   call daxpy(NBasis**2,2d0*val,ints,1,Jtmp,1)
enddo

J_mo = 0.5d0*Jtmp

deallocate(MatFFMO,ints)

end subroutine CholeskyOTF_Jmat_MO

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
integer :: ierr
double precision            :: val1,val2
double precision, parameter :: ThreshH0 = 1d-6

val1 = norm2(H0_int)
val2 = norm2(H0_mo)

write(6,'(1x,a,f12.6,/)') "H0 (internal) vs. H0 (external) = ",abs(val1)-abs(val2)
write(6,*) "H0 (internal) vs. H0 (external) = ",abs(val1)-abs(val2)

if(abs(val1)-abs(val2).gt.ThreshH0) then
  print*, 'Difference in H0 norms = '
  print*, abs(val1)-abs(val2)
  print*, 'Threshold = ',ThreshH0
  print*, 'Check for errors in geometry / basis set?'

  write(6,*) 'H0/MO (internal)', norm2(H0_int)
  do i=1,NBasis
     write(LOUT,'(*(f13.8))') (H0_int(i,j),j=1,NBasis)
  enddo

  print*, 'H0/MO (external) ', norm2(H0_mo)
  do i=1,NBasis
     write(LOUT,'(*(f13.8))') (H0_mo(i,j),j=1,NBasis)
  enddo

  stop
endif

! norms may be ok but...
ierr = 0
do j=1,NBasis
   do i=1,NBasis
      val1 = abs(H0_int(i,j)) - abs(H0_mo(i,j))
      if (val1 .gt. 1d-4) then
         ierr = ierr + 1 
        ! print*, i,j, H0_int(i,j),H0_mo(i,j)
      endif
   enddo
enddo
if (ierr .gt. 0) then
  write(6,'(1x,"H0 ext and int differ in ", I5, " places!")') ierr
  stop
endif

end subroutine CholeskyOTF_H0_test

subroutine DipMomOTF_ao(Dx_extao,Dy_extao,Dz_extao,BasisSetPath,XYZPath,Units,Source)
!
! calculate dipole moments
!
use Multipoles

character(*), intent(in) :: BasisSetPath
character(*), intent(in) :: XYZPath
character(6),intent(in)  :: Source
integer, intent(in)      :: Units
real(F64), dimension(:, :), allocatable :: Dx_extao, Dy_extao, Dz_extao

type(TSystem)  :: System
type(TAOBasis) :: AOBasis
real(F64), dimension(3) :: Rc
real(F64), dimension(:, :), allocatable :: Dx, Dy, Dz
!
integer :: NAO
integer            :: ORBITAL_ORDERING
logical, parameter :: SortAngularMomenta = .false.
logical, parameter :: SpherAO = .true.

! set orbital ordering
if(trim(Source)=='MOLPRO') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
elseif(trim(Source)=='ORCA  ') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_ORCA
elseif(trim(Source)=='DALTON') then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
endif
print*, 'orbital_ordering',orbital_ordering

call auto2e_init()
!
! Initialize the Boys function interpolation table
! (used for Coulomb integrals evaluation).
!
call boys_init(4 * AUTO2E_MAXL)
call sys_Read_XYZ(System, XYZPath, Units)
call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
NAO = AOBasis%NAOSpher

!
! Compute the charge center of the nuclei
! Dipole moments will be computed with respect to Rc
!
call sys_ChargeCenter(Rc, System)
!

! Calculate x, y, and z electronic dipole moment matrices
! in the spherical AO gaussian Basis. The ordering of orbitals
! follows the gammcor-integrals convention.
!
allocate(Dx(NAO, NAO))
allocate(Dy(NAO, NAO))
allocate(Dz(NAO, NAO))
Rc(1)=0
Rc(2)=0
Rc(3)=0
print*, 'use 0 0 0 ...'
call multi_ElectronicDipole(Dx, Dy, Dz, Rc, AOBasis)
!
! Convert Dx, Dy, Dz matrices
! to the AO basis with Molpro/Dalton's ordering
!
allocate(Dx_extao(NAO, NAO))
allocate(Dy_extao(NAO, NAO))
allocate(Dz_extao(NAO, NAO))
call auto2e_interface_AngFuncTransf(Dx_extao, Dx, .false., .true., AOBasis, ORBITAL_ORDERING)
call auto2e_interface_AngFuncTransf(Dy_extao, Dy, .false., .true., AOBasis, ORBITAL_ORDERING)
call auto2e_interface_AngFuncTransf(Dz_extao, Dz, .false., .true., AOBasis, ORBITAL_ORDERING)

#if CHOLOTF_DEBUG > 4
   call msg("---------------- dipole matrices ---------------------")
   print*, 'DIPX AO',norm2(Dx_extao)
   call msg("Dx")
   call geprn(Dx_extao)
   print*, 'DIPY AO',norm2(Dy_extao)
   call msg("Dy")
   call geprn(Dy_extao)
   print*, 'DIPZ AO',norm2(Dz_extao)
   call msg("Dz")
   call geprn(Dz_extao)
#endif

call boys_free()

end subroutine DipMomOTF_ao

subroutine CompDipMomOTF(AOBasis,System,CAONO,Occ,DipX,DipY,DipZ,NAO,NMO)
!
! calculate dipole moments
! Input  : DipX, DipY, DipZ in AO
!
type(TSystem)  :: System
type(TAOBasis) :: AOBasis
!
integer,intent(in)  :: NAO,NMO
real(F64),intent(in) :: CAONO(NAO,NMO),Occ(NMO)
real(F64),dimension(NAO,NAO),intent(in)  :: DipX,DipY,DipZ

integer :: l
real(F64) :: Qc
real(F64) :: DM_X, DM_Y, DM_Z, DXYZ
real(F64) :: NUC_DMX, NUC_DMY, NUC_DMZ
real(F64), dimension(3) :: Rc
real(F64), dimension(NMO,NMO) :: Dx, Dy, Dz
real(F64), dimension(NMO,NAO) :: AUXM

associate( &
 NAtoms   => System%NAtoms, &
 ZNumbers => System%ZNumbers, &
 AtomCoords => AOBasis%AtomCoords &
 )
NUC_DMX=0; NUC_DMY=0; NUC_DMZ=0
do l = 1,NAtoms
      Qc = real(ZNumbers(l), F64)
      !write(6,'(i3,4f12.6)')l,Rc(1),Rc(2),Rc(3),Qc
      NUC_DMX = NUC_DMX + Qc*AtomCoords(1,l)
      NUC_DMY = NUC_DMY + Qc*AtomCoords(2,l)
      NUC_DMZ = NUC_DMZ + Qc*AtomCoords(3,l)
enddo
end associate

Call dgemm('T','N',NAO,NMO,NMO,1d0,CAONO,NMO,DipX,NMO,0d0,AUXM,NMO)
Call dgemm('N','N',NMO,NMO,NMO,1d0,AUXM,NMO,CAONO,NMO,0d0,Dx,NMO)
!
Call dgemm('T','N',NAO,NMO,NMO,1d0,CAONO,NMO,DipY,NMO,0d0,AUXM,NMO)
Call dgemm('N','N',NMO,NMO,NMO,1d0,AUXM,NMO,CAONO,NMO,0d0,Dy,NMO)
!
Call dgemm('T','N',NAO,NMO,NMO,1d0,CAONO,NMO,DipZ,NAO,0d0,AUXM,NMO)
Call dgemm('N','N',NMO,NMO,NAO,1d0,AUXM,NMO,CAONO,NAO,0d0,Dz,NMO)

DM_X=0d0; DM_Y=0d0; DM_Z=0d0
Do I=1,NMO
   DM_X = DM_X - 2d0*Occ(i)*Dx(i,i)
   DM_Y = DM_Y - 2d0*Occ(i)*Dy(i,i)
   DM_Z = DM_Z - 2d0*Occ(i)*Dz(i,i)
EndDo

Write(6,'(/1X,"Nuclear Dipole Moment   ",3f12.8)')  NUC_DMX,NUC_DMY,NUC_DMZ
Write(6,'(1X,"Electronic Dipole Moment",3f12.8)')   DM_X,DM_Y,DM_Z
Write(6,'(1X,"Total Dipole Moment     ",3f12.8,/)') NUC_DMX+DM_X,NUC_DMY+DM_Y,NUC_DMZ+DM_Z

DXYZ=SQRT((NUC_DMX+DM_X)**2+(NUC_DMY+DM_Y)**2+(NUC_DMZ+DM_Z)**2)

Write(6,'(1X,A,2f12.8,/)') '|dipole moment| a.u./D', DXYZ, DXYZ/0.393456

end subroutine CompDipMomOTF

subroutine reorder_by_first_row(A, B, matrix_to_reorder, N, reordered, column_map, sign_map, tolerance)
! 
! Purpose: reorder either A or B so that columns match
!          with tolerance |A(i,j)|-B(i,j)| < tolerance
!          and match phases change (change sign of the whole column)
!
! Co-Generated by ChatGPT (06.02.25)
!
    implicit none
    integer, intent(in) :: N
    real(8), intent(in) :: A(N, N), B(N, N)
    character(len=*), intent(in) :: matrix_to_reorder
    real(8), intent(out) :: reordered(N, N)
    integer, intent(out) :: column_map(N) ! Mapping of column indices
    real(8), intent(out) :: sign_map(N)   ! Mapping of phases
    real(8), intent(in)  :: tolerance     ! Allowed numerical difference

    integer :: j, k, found
    integer :: i
    real(8) :: ref_matrix(N, N), target_matrix(N, N)
    logical :: match

    ! Determine which matrix needs reordering
    if (trim(matrix_to_reorder) == "A") then
        ref_matrix = B
        target_matrix = A
    else if (trim(matrix_to_reorder) == "B") then
        ref_matrix = A
        target_matrix = B
    else
        print *, "Error: Invalid matrix selection. Choose 'A' or 'B'."
        stop
    end if

#if CHOLOTF_DEBUG > 10
   print*, 'before reorder: i ref target'
   do i=1,N
      write(6,'(1x,i3,2f12.6)') i, ref_matrix(1,i), target_matrix(1,i)
   enddo
#endif

    ! Initialize column mapping
    column_map(:) = -1
    sign_map(:)   = 1d0

    ! Compare only the first row to determine column mapping
    do j = 1, N
        found = 0
        do k = 1, N
            if (column_map(k) == -1) then  ! Check only unused columns
                if (abs(abs(ref_matrix(1, j)) - abs(target_matrix(1, k))) <= tolerance) then
                    column_map(k) = j
                    if ( abs(ref_matrix(1, j) - target_matrix(1, k)) > tolerance ) then
                       sign_map(k) = -1d0
                    endif
                    found = 1
                    exit
                end if
            end if
        end do

        if (found == 0) then
            print *, "Warning: No matching column found for column ", j
        end if
    end do

    ! Apply column permutation with sign change
    do j = 1, N
        if (column_map(j) /= -1) then
            !reordered(:, j) = target_matrix(:, column_map(j))
            reordered(:, column_map(j)) =  sign_map(j) * target_matrix(:, j)
        else
            print *, "Error: Could not reorder column ", j
            print*, 'ref(1,j) =', ref_matrix(1,j)
            stop
        end if

    end do

    !print*, 'col_map, sign_map'
    !do i=1,n
    !   print*, i, column_map(i), sign_map(i)
    !enddo

#if CHOLOTF_DEBUG > 10
   print*, 'after reorder: i ref target'
   do i=1,N
      write(6,'(1x,i3,2f12.6)') i, ref_matrix(1,i), reordered(1,i)
   enddo
#endif

end subroutine reorder_by_first_row

subroutine test_CtSC(C,AOBasis,NAO,NBasis)
!
! check if C^T(AO,MO).S.C(AO,MO) = 1
! working : works with Orca now!
!
implicit none

integer,intent(in) :: NAO,NBasis
double precision,intent(in) :: C(NAO,NBasis)
type(TAOBasis), intent(in)  :: AOBasis

integer :: i,j
double precision :: val
double precision :: S_ao(NAO,NAO)
double precision :: S_extao(NAO,NAO)
double precision :: work(NBasis,NBasis)

call ints1e_S(S_ao, AOBasis)
call auto2e_interface_ApplyOrcaPhases_Matrix(S_ao, AOBasis, .true.)
call auto2e_interface_AngFuncTransf(S_extao, S_ao, .false., .true., AOBasis, ORBITAL_ORDERING_ORCA)
call dgemm('T','N',NBasis,NAO,NAO,1d0,C,NBasis, &
           S_extao,NAO,0d0,work,NBasis)
call dgemm('N','N',NBasis,NBasis,NAO,1d0,work,NBasis, &
           C,NBasis,0d0,S_ao,NBasis)

do j=1,NBasis
   do i=j+1,NBasis
      val = abs(S_ao(i,j))
      if (val .gt. 1d-6) then
         print*, i,j,val
      endif
   enddo
enddo

#if CHOLOTF_DEBUG > 5
   print*, 'UAONO : C^T.S.C norm =', norm2(S_AO)
   do i=1,NBasis
      write(LOUT,'(*(f13.8))') (S_AO(i,j),j=1,NBasis)
   enddo
#endif

end subroutine test_CtSC

end module choleskyOTF_interface
