module grid_internal

use gammcor_integrals

contains

!subroutine internal_lda_ao_grid(BasisSetPath,Xg,Yg,Zg,Wg,NPoints)
!!
!! returns : a) grid weights, Wg(Npoints)
!!           b) grid points, Xg, Yg, Zg
!!
!implicit none
!
!character(:),allocatable :: BasisSetPath
!
!integer,intent(out) :: NPoints
!double precision, dimension(:), allocatable :: Xg, Yg, Zg
!double precision, dimension(:), allocatable :: Wg
!
!character(:),allocatable :: XYZPath
!logical :: SortAngularMomenta
!
!type(TAOBasis)  :: AOBasis
!type(TSystem)   :: System
!
!logical, parameter :: SpherAO = .true.
!integer, parameter :: GridType = BECKE_PARAMS_MEDIUM
!
!XYZPath = "./input.inp"
!SortAngularMomenta = .true.
!
!call auto2e_init()
!call sys_Read_XYZ(System, XYZPath)
!call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
!!if (AOBasis%SpherAO) then
!!      NAO = AOBasis%NAOSpher
!!else
!!      NAO = AOBasis%NAOCart
!!end if
!
!! Molecular grid
!call becke_MolecularGrid(Xg, Yg, Zg, Wg, NPoints, GridType, System, AOBasis)
!
!end subroutine internal_lda_ao_grid

subroutine internal_lda_ao_orbgrid(GridType,BasisSetPath,AOBasis,Wg,Phi,NPoints,NAO,Units)
!
! returns : a) grid weights, Wg(Npoints)
!           b) AO orbitals on a grid, Phi(NPoints,NAO)
!
implicit none

type(TAOBasis)  :: AOBasis
type(TSystem)   :: System

integer,intent(in)  :: GridType
integer,intent(in)  :: NAO
character(*) :: BasisSetPath
!character(:),allocatable :: BasisSetPath
integer,intent(out) :: NPoints
integer,optional :: Units

double precision, dimension(:,:), allocatable :: Phi
double precision, dimension(:), allocatable :: Wg

integer :: Units0, NAOt
double precision, dimension(:), allocatable :: Xg, Yg, Zg
logical :: SortAngularMomenta
character(:),allocatable :: XYZPath

logical, parameter :: SpherAO = .true.
!integer, parameter :: GridType = BECKE_PARAMS_MEDIUM

XYZPath = "./input.inp"
SortAngularMomenta = .true.

! set units
if (present(Units)) then
   Units0 =  Units
else
   Units0 = SYS_UNITS_ANGSTROM
endif

call auto2e_init()
call sys_Read_XYZ(System, XYZPath, Units0)
call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
if (AOBasis%SpherAO) then
      NAOt = AOBasis%NAOSpher
else
      NAOt = AOBasis%NAOCart
end if
if(NAOt /= NAO) then
  print*, 'NAO =',NAO, 'NAOlib',NAOt
  stop "sth wrong with NAO in internal_orbgrid!"
endif

! Molecular grid
call becke_MolecularGrid(Xg, Yg, Zg, Wg, NPoints, GridType, System, AOBasis)

!call internal_lda_ao_grid(BasisSetPath,Xg,Yg,Zg,Wg,NPoints)

! Atomic orbitals on the grid
allocate(Phi(NPoints, NAO))
call gridfunc_Orbitals(Phi, Xg, Yg, Zg, NPoints, NAO, AOBasis)

end subroutine internal_lda_ao_orbgrid

subroutine internal_gga_ao_orbgrid(GridType,BasisSetPath,AOBasis,Wg,Phi,NPoints,NAO,Units)
!
! returns : a) grid weights, Wg(Npoints)
!           b) Atomic orbitals and gradient components on the grid
!              Phi(k, p, i)
!              k grid point
!              p index of AO
!              i = 1 (value), 2 (d/dx), 3 (d/dy), 4 (d/dz)
!
implicit none

type(TAOBasis)  :: AOBasis
type(TSystem)   :: System

integer,intent(in)  :: GridType
integer,intent(in)  :: NAO
character(*)        :: BasisSetPath
!character(:),allocatable :: BasisSetPath
integer,intent(out) :: NPoints

double precision, dimension(:,:,:), allocatable :: Phi
double precision, dimension(:), allocatable :: Wg

integer,optional :: Units

integer :: NAOt, Units0
double precision, dimension(:), allocatable :: Xg, Yg, Zg
character(:),allocatable :: XYZPath
logical :: SortAngularMomenta

logical, parameter :: SpherAO = .true.
!integer, parameter :: GridType = BECKE_PARAMS_SG1
!integer, parameter :: GridType = BECKE_PARAMS_MEDIUM
!integer, parameter :: GridType = BECKE_PARAMS_FINE
!integer, parameter :: GridType = BECKE_PARAMS_XFINE

XYZPath = "./input.inp"
SortAngularMomenta = .true.

! set units
if (present(Units)) then
   Units0 =  Units
else
   Units0 = SYS_UNITS_ANGSTROM
endif

call auto2e_init()
call sys_Read_XYZ(System, XYZPath, Units0)
call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
if (AOBasis%SpherAO) then
      NAOt = AOBasis%NAOSpher
else
      NAOt = AOBasis%NAOCart
end if
if(NAOt /= NAO) then
  print*, 'NAO =',NAO, 'NAOlib',NAOt
  stop "sth wrong with NAO in internal_orbgrid!"
endif

! Molecular grid
call becke_MolecularGrid(Xg, Yg, Zg, Wg, NPoints, GridType, System, AOBasis)

! Atomic orbitals on the grid
allocate(Phi(NPoints, NAO, 4))
call gridfunc_Orbitals_Grad(Phi, Xg, Yg, Zg, AOBasis)
!
! future work: loop over grid points, e.g. 1000 grid points to get many loops?
!call gridfunc_Orbitals_Grad(Phi(1:NBatch,NAO,4), Xg(x0:x1), Yg(x0:x1), Zg(x0:x1), AOBasis)

end subroutine internal_gga_ao_orbgrid

subroutine internal_lda_no_orbgrid(GridType,BasisSetPath,ExternalOrdering,CAONO,Wg,Phi,NPoints,NAO,NBasis,Units)
!
! Compute orbital on grid and perform AO-->NO transformation
!
! returns : a) no of grid points, NPoints
!           b) grid weights, Wg(Npoints)
!           c) NO orbitals on a grid, Phi(NPoints,NO)
!
integer,intent(out) :: NPoints
double precision, dimension(:,:), allocatable, intent(out) :: Phi
double precision, dimension(:), allocatable, intent(out) :: Wg

type(TAOBasis) :: AOBasis

integer,intent(in) :: GridType
integer,intent(in) :: NAO,NBasis
integer,intent(in) :: ExternalOrdering

double precision,intent(in) :: CAONO(NAO,NBasis)
character(*)        :: BasisSetPath
!character(:),allocatable    :: BasisSetPath

integer, optional :: Units

integer :: Units0

! set units
if (present(Units)) then
   Units0 =  Units
else
   Units0 = SYS_UNITS_ANGSTROM
endif

call internal_lda_ao_orbgrid(GridType,BasisSetPath,AOBasis,Wg,Phi,NPoints,NAO,Units0)
call internal_tran_lda_orbgrid(CAONO,Phi,AOBasis,ExternalOrdering,NPoints,NAO,NBasis)

end subroutine internal_lda_no_orbgrid

subroutine internal_gga_no_orbgrid(GridType,BasisSetPath,ExternalOrdering,CAONO, &
                                   Wg,Phi,NPoints,NAO,NBasis,Units)
!
! Compute orbital on grid and perform AO-->NO transformation
!
! returns : a) no of grid points, NPoints
!           b) grid weights, Wg(Npoints)
!           c) NO orbitals on a grid, Phi(NPoints,NO)
!
integer,intent(out) :: NPoints
double precision, dimension(:,:,:), allocatable, intent(out) :: Phi
double precision, dimension(:), allocatable, intent(out) :: Wg

type(TAOBasis) :: AOBasis

integer,intent(in) :: GridType
integer,intent(in) :: NAO,NBasis
integer,intent(in) :: ExternalOrdering
double precision,intent(in) :: CAONO(NAO,NBasis)
character(*) :: BasisSetPath
!character(:),allocatable    :: BasisSetPath

integer, optional :: Units

integer :: Units0

! set units
if (present(Units)) then
   Units0 =  Units
else
   Units0 = SYS_UNITS_ANGSTROM
endif

call internal_gga_ao_orbgrid(GridType,BasisSetPath,AOBasis,Wg,Phi,NPoints,NAO,Units0)
call internal_tran_gga_orbgrid(CAONO,Phi(:,:,1),Phi(:,:,2),Phi(:,:,3),Phi(:,:,4), &
                               AOBasis,ExternalOrdering,NPoints,NAO,NBasis)

!! additional tests
!block
!integer :: k
!double precision :: RhoIntegral,DivRhoIntegral
!double precision,allocatable :: Rho(:,:)
!
!allocate(Rho(NPoints, 4))
!call gridfunc_GGA_Variables(Rho, Phi, CAONO, OccNumbers)
!
!DivRhoIntegral = 0d0
!RhoIntegral = 0d0
!do k = 1, NPoints
!      RhoIntegral = RhoIntegral + Wg(k) * Rho(k, 1)
!      DivRhoIntegral = DivRhoIntegral + Wg(k) * (Rho(k, 2) + Rho(k, 3) + Rho(k, 4))
!end do
!
!print*, 'RhoIntegral    = ',RhoIntegral
!print*, 'DivRhoIntegral = ',DivRhoIntegral
!
!end block

end subroutine internal_gga_no_orbgrid

subroutine internal_tran_lda_orbgrid(CAONO,Phi,AOBasis,ExternalOrdering,NPoints,NAO,NBasis)
!
! input : Phi = (NGrid,AO); CAONO(AO,NO)
! output: Phi(NGrid,NO) = Phi.CAONO
!
! careful :: it is assumed that NAO=NBasis
!
implicit none

type(TAOBasis)     :: AOBasis
integer,intent(in) :: NPoints, NBasis, NAO
integer,intent(in) :: ExternalOrdering
double precision,intent(in)    :: CAONO(NAO,NBasis)
double precision,intent(inout) :: Phi(NPoints,NBasis)

double precision :: C_ao(NAO,NAO)
double precision :: work(NPoints,NAO)

if(NBasis/=NAO) stop "NAO.ne.NBasis in internal_tran_orbgrid!"

! AOs from external program -> AOs in the Auto2e format
call auto2e_interface_C(C_ao, CAONO, AOBasis, ExternalOrdering)

!print*, 'CAONO',norm2(CAONO)
!print*, 'C_ao',norm2(C_ao)

work = Phi
call dgemm('N','N',NPoints,NBasis,NAO,1d0,work,NPoints,C_ao,NAO,0d0,Phi,NPoints)
!print*, 'Phi',norm2(Phi)

end subroutine internal_tran_lda_orbgrid

subroutine internal_tran_gga_orbgrid(CAONO,Phi,PhiX,PhiY,PhiZ,AOBasis,ExternalOrdering,NPoints,NAO,NBasis)
!
! input : Phi,PhiX=d/dx Phi; CAONO(AO,NO)
! output: Phi(NGrid,NO) = Phi.CAONO
!
! careful :: it is assumed that NAO=NBasis
!
implicit none

type(TAOBasis)     :: AOBasis
integer,intent(in) :: NPoints, NBasis, NAO
integer,intent(in) :: ExternalOrdering
double precision,intent(in)    :: CAONO(NAO,NBasis)
double precision,intent(inout) :: Phi(NPoints,NBasis)
double precision,intent(inout) :: PhiX(NPoints,NBasis),PhiY(NPoints,NBasis),PhiZ(NPoints,NBasis)

double precision :: C_ao(NAO,NAO)
double precision :: work(NPoints,NAO)

if(NBasis/=NAO) stop "NAO.ne.NBasis in internal_tran_orbgrid!"

! AOs from external program -> AOs in the Auto2e format
call auto2e_interface_C(C_ao, CAONO, AOBasis, ExternalOrdering)

work = Phi
call dgemm('N','N',NPoints,NBasis,NAO,1d0,work,NPoints,C_ao,NAO,0d0,Phi,NPoints)
work = PhiX
call dgemm('N','N',NPoints,NBasis,NAO,1d0,work,NPoints,C_ao,NAO,0d0,PhiX,NPoints)
work = PhiY
call dgemm('N','N',NPoints,NBasis,NAO,1d0,work,NPoints,C_ao,NAO,0d0,PhiY,NPoints)
work = PhiZ
call dgemm('N','N',NPoints,NBasis,NAO,1d0,work,NPoints,C_ao,NAO,0d0,PhiZ,NPoints)

!print*, 'Phi  ',norm2(Phi)
!print*, 'PhiX ',norm2(PhiX)
!print*, 'PhiY ',norm2(PhiY)
!print*, 'PhiZ ',norm2(PhiZ)

end subroutine internal_tran_gga_orbgrid

end module grid_internal
