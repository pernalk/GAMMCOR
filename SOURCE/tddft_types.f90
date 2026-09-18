module tddft_types

      use types
      implicit none


      type TtddftData
            integer :: NBasis, NI, NV, NEL
            integer :: NCoreOrb = 0
            integer :: NInte1, NInte2
            double precision :: ENuc
            double precision :: EKS
            double precision :: EMP2 = 0.d0
            double precision :: EMP2_pyscf = 0.d0
            double precision :: EDH_pyscf = 0.d0
            double precision :: XMP2 = 0.d0
            double precision :: XHF = 0.d0
            double precision :: Omega = 0.d0
            double precision :: Alpha = 0.d0
            character(len=256) :: XC = ''
            character(len=16) :: XCType = ''
            double precision, dimension(:), allocatable :: eorbs
            double precision, dimension(:), allocatable :: Occ
            double precision, dimension(:, :), allocatable :: CAOMO
            double precision, dimension(:, :), allocatable :: HAO
            double precision, dimension(:, :), allocatable :: HMO
            double precision, dimension(:, :), allocatable :: FockMO
            double precision, dimension(:, :), allocatable :: rdm1_p
            double precision, dimension(:, :), allocatable :: rdm1_m
            double precision, dimension(:, :), allocatable :: rdm1_full
            integer :: NGrid = 0
            integer :: NGridComp = 0
            double precision, dimension(:, :), allocatable :: RGrid
            double precision, dimension(:), allocatable :: WGrid
            double precision, dimension(:, :), allocatable :: OrbGrid
            double precision, dimension(:, :), allocatable :: OrbXGrid
            double precision, dimension(:, :), allocatable :: OrbYGrid
            double precision, dimension(:, :), allocatable :: OrbZGrid
            double precision, dimension(:, :), allocatable :: RhoGrid
            double precision, dimension(:, :), allocatable :: VxcGrid
            double precision, dimension(:, :, :), allocatable :: FxcS
            double precision, dimension(:, :, :), allocatable :: FxcT
            double precision, dimension(:, :), allocatable :: VxcAO
            integer :: NExcS = 0
            integer :: NExcT = 0
            double precision, dimension(:), allocatable :: ExcS
            double precision, dimension(:), allocatable :: ExcT
            double precision, dimension(:), allocatable :: OscS
            double precision, dimension(:, :, :), allocatable :: XS
            double precision, dimension(:, :, :), allocatable :: YS
            double precision, dimension(:, :, :), allocatable :: XT
            double precision, dimension(:, :, :), allocatable :: YT
            double precision, dimension(:, :), allocatable :: AS
            double precision, dimension(:, :), allocatable :: BS
            double precision, dimension(:, :), allocatable :: AT
            double precision, dimension(:, :), allocatable :: BT

      end type TtddftData

end module tddft_types
