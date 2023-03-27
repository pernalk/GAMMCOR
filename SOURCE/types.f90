module types
! written by M. Hapka, M. Modrzejewski
!            K. Pernal

use print_units

integer, parameter :: INTER_TYPE_DAL  = 1
integer, parameter :: INTER_TYPE_MOL  = 2
integer, parameter :: INTER_TYPE_OWN  = 3
integer, parameter :: INTER_TYPE_ORCA = 4
integer, parameter :: INTER_TYPE_TREX = 5

integer, parameter :: TYPE_NO_SYM = 1
integer, parameter :: TYPE_SYM = 0

integer, parameter :: JOB_TYPE_AC    = 1
integer, parameter :: JOB_TYPE_AC0   = 2
integer, parameter :: JOB_TYPE_AC1   = 3
integer, parameter :: JOB_TYPE_ERPA  = 3
integer, parameter :: JOB_TYPE_EERPA = 4
integer, parameter :: JOB_TYPE_SAPT  = 5
integer, parameter :: JOB_TYPE_PDFT  = 6
integer, parameter :: JOB_TYPE_CASPIDFT    = 7
integer, parameter :: JOB_TYPE_CASPIDFTOPT = 8
integer, parameter :: JOB_TYPE_EERPA_OLD   = 9
integer, parameter :: JOB_TYPE_AC0D        = 10
integer, parameter :: JOB_TYPE_AC0DNOSYMM  = 11
integer, parameter :: JOB_TYPE_NLOCCORR    = 12
integer, parameter :: JOB_TYPE_AC0DP       = 13
integer, parameter :: JOB_TYPE_ACFREQ      = 14
integer, parameter :: JOB_TYPE_ACFREQNTH   = 15
integer, parameter :: JOB_TYPE_AC1FREQNTH  = 16
integer, parameter :: JOB_TYPE_RESPONSE    = 17

integer, parameter :: SAPTLEVEL0 = 0
integer, parameter :: SAPTLEVEL1 = 1
integer, parameter :: SAPTLEVEL2 = 2

integer, parameter :: FLAG_CORE = 1
integer, parameter :: FLAG_NOBASIS = 0
integer, parameter :: FLAG_REDVIRT = 0
integer, parameter :: FLAG_RDMCORR = 0
integer, parameter :: FLAG_RDM2TYP = 0
logical, parameter :: FLAG_RESTART = .FALSE.
logical, parameter :: FLAG_TRIPLET = .FALSE.
integer, parameter :: FLAG_PRINT_LEVEL = 0
integer, parameter :: FLAG_DEBUG_FL12  = 1

integer, parameter :: RDM_TYPE_GVB  = 1
integer, parameter :: RDM_TYPE_APSG = 2
integer, parameter :: RDM_TYPE_CAS  = 3
integer, parameter :: RDM_TYPE_DMRG = 4
integer, parameter :: RDM_TYPE_HF   = 5
integer, parameter :: RDM_TYPE_CI   = 6

integer, parameter :: TWOMO_INCORE = 1
integer, parameter :: TWOMO_FFFF   = 2
integer, parameter :: TWOMO_FOFO   = 3

integer, parameter :: FLAG_CHOLESKY = 0
integer, parameter :: CHOL_ACCU_DEFAULT   = 1
integer, parameter :: CHOL_ACCU_TIGHT     = 2
integer, parameter :: CHOL_ACCU_LUDICROUS = 3

integer, parameter :: DF_NONE       = 0
integer, parameter :: DF_SRLDA      = 1
integer, parameter :: DF_SRPBE      = 2
integer, parameter :: DF_PBE        = 3

integer, parameter :: MONOMER_A = 1
integer, parameter :: MONOMER_B = 2

integer, parameter :: RESP_ERPA = 1
integer, parameter :: RESP_APSG = 2
integer, parameter :: RESP_DFT  = 3

logical, parameter :: FLAG_POSTCAS  = .FALSE.

integer,parameter :: maxcen = 500

character(*),parameter :: PossibleInterface(5) = &
[character(8) :: &
'DALTON', 'MOLPRO', 'OWN', 'ORCA','TREXIO']

character(*),parameter :: PossibleJobType(17) = &
[character(9) :: &
'AC', 'AC0', 'ERPA', 'EERPA', 'SAPT', 'PDFT', 'CASPiDFT','CASPiDFTOpt','EERPA-1', & 
'AC0D', 'AC0DNOSYMM', 'NLOCCORR', 'AC0DP', 'ACFREQ','ACFREQNTH','AC1FREQNTH','RESPONSE']

character(*),parameter :: PossibleRDMType(6) = &
[character(8) :: &
'GVB', 'APSG', 'CASSCF', 'DMRG', 'HF', 'CI']

character(*),parameter :: PossibleDFAType(3) = &
[character(8) :: &
'srLDA', 'srPBE', 'PBE']

character(*),parameter :: PossibleCholAccu(3) = &
[character(8) :: &
'DEFAULT', 'TIGHT', 'LUDICROUS']

character(*),parameter :: PossibleMonomers(2) = &
[character(8) :: 'A', 'B']

character(:), allocatable :: InputPath
!InputPath = "./input.inp"

type CalculationBlock
      integer :: InterfaceType = INTER_TYPE_DAL
      integer :: NBasis    = FLAG_NOBASIS
      integer :: JobType   = JOB_TYPE_AC
      integer :: RDMType ! = RDM_TYPE_GVB
      integer :: RDMSource = INTER_TYPE_DAL
      integer :: Response  = RESP_ERPA
      integer :: DFApp     = DF_NONE
      integer :: Kernel    = 1
      integer :: TwoMoInt  = TWOMO_INCORE
      integer :: Cholesky  = FLAG_CHOLESKY 
      integer :: CholeskyAccu = CHOL_ACCU_DEFAULT
      integer :: Core      = FLAG_CORE
      integer :: SymType   = TYPE_NO_SYM
      integer :: SaptLevel = SAPTLEVEL2
      integer :: vdWCoef   = 0
      integer :: RedVirt   = FLAG_REDVIRT
      integer :: RdmCorr   = FLAG_RDMCORR
      integer :: Rdm2Type  = FLAG_RDM2TYP
      integer :: MemVal = 2, MemType = 3 ! default: use 2 GB for 3-ind_tran (Cholesky)
      logical :: Restart    = FLAG_RESTART
      logical :: Triplet    = FLAG_TRIPLET
      logical :: PostCAS    = FLAG_POSTCAS
      integer :: IPrint     = 0
      double precision :: RPAThresh  = 1.0D-6
      double precision :: ThreshVirt = 1.0D-6
      integer :: imon = 1
      character(:), allocatable :: JobTitle
      character(:), allocatable :: IntegralsFilePath
      integer :: Max_Cn = 3
      double precision :: FreqOm = 0.d0
      logical :: CAlpha = .false.
end type CalculationBlock

type SystemBlock
      integer :: NoSt     = 1
      integer :: NStates  = 1
      integer :: ISpinMs2 = 0
      integer :: EigFCI = 1
      integer :: Charge = 0
      integer :: ZNucl  = 0
      integer :: NBasis = 0
      integer :: NChol  = 0
      integer :: Monomer = MONOMER_A
      integer :: NELE
      double precision :: XELE
      double precision :: PotNuc
      double precision :: SumOcc  = 0d0 
      double precision :: ACAlpha = 1d0
      double precision :: Omega   = 1d0
      double precision :: PerVirt = 0d0
      double precision :: ECASSCF = 0d0
      integer :: NSym
      integer :: NSymBas(8),NSymOrb(8)
      integer :: NOrb, NGem
      integer :: NActOrb = 1
      integer :: NAct, INAct
      integer :: ISwitchAct = 0
      integer :: NActS(8), INActS(8)
      integer :: NDim, NDimX
      integer :: NDimN, DimEx
      integer :: NDimX0  = 0
      integer :: NOccup0 = 0
      integer :: NCen = 0
      integer :: UCen = 0
      integer :: NMonBas(8) = 0
      integer :: IPrint = 0
      integer :: IWarn = 0
      integer :: icnt
      integer :: num0,num1,num2
      integer :: NVZero = 0
      integer :: TwoMoInt = TWOMO_INCORE
      logical :: DeclareTwoMo     = .false.
      logical :: DeclareSt        = .false.
      logical :: DeclareTrSt      = .false.
      logical :: DeclareSpin      = .false.
      logical :: DeclareThrSelAct = .false.
      logical :: DeclareThrQVirt  = .false.
      logical :: DeclareThrQInact = .false.
      logical :: ISHF    = .false.
      logical :: Cubic   = .false.
      logical :: Wexcit  = .false.
      logical :: doRSH   = .false., SameOm = .true.
      logical :: PostCAS = .false.
      logical :: NActFromRDM  = .true.
      logical :: reduceV      = .false.
      logical :: Cholesky2RDM = .false.
      ! for RDMCORR SAPT
      integer :: RDModel = -1
      ! for cubic SAPT
      double precision :: ACAlpha0  = 1.d-10
      double precision :: ACAlpha1  = 0.01d0
      double precision :: ACAlpha2  = 0.45d0

      ! ThrGemAct for active geminals selection
      double precision :: ThrGemAct    = 0.992d0
      ! ThrAct for active orbitals selection (CIPSI)
      double precision :: ThrAct = 1.d-9
      ! ThrSelAct selects active orbital pairs
      double precision :: ThrSelAct = 1.d-8
      ! ThrVirt for reduction of virtual orbs
      double precision :: ThrVirt   = 1.d-6
      ! ThrQVirt for quasi-virtual orbs
      double precision :: ThrQVirt  = 1.d-7
      ! ThrQInact for quasi-inactive orbs (0.9998...)
      double precision :: ThrQInact = 1.d-4

      integer,allocatable :: InSt(:,:),InTrSt(:,:)
      integer,allocatable :: IGem(:), IndAux(:)
      integer,allocatable :: IndX(:), IndN(:,:), IPair(:,:)
      integer,allocatable :: MultpC(:,:),NSymNO(:)
      integer,allocatable :: IndXh(:)
      integer,allocatable :: NumOSym(:),IndInt(:)
      integer,allocatable :: IndNx(:,:)
      ! TEST ONLY
      integer,allocatable :: IGem0(:)
      double precision,allocatable :: Occ0(:)

      integer,allocatable :: IndNT(:,:)
      integer,allocatable :: Ind2(:)
      double precision,allocatable :: Occ(:), CICoef(:)
      double precision,allocatable :: OrbE(:)
      double precision,allocatable :: TwoMO(:)
      double precision,allocatable :: CMO(:,:)
      double precision,allocatable :: CNONO(:,:)
      double precision,allocatable :: OV(:,:),OO(:,:), &
                                      FO(:,:),FF(:,:), &
                                      FFAB(:,:),FFBA(:,:), &
                                      OOAB(:,:),OOBA(:,:)
      double precision,allocatable :: DChol(:,:)
      double precision,allocatable :: Pmat(:,:)
      double precision,allocatable :: WPot(:,:),Kmat(:,:)
      double precision,allocatable :: VCoul(:)
      double precision,allocatable :: RDM1c(:,:)
      double precision,allocatable :: RDM2(:)
      double precision,allocatable :: RDM2val(:,:,:,:)
      double precision,allocatable :: Fmat(:,:) 
      double precision,allocatable :: dipm(:,:,:)
      double precision,allocatable :: Eig(:),EigX(:),EigY(:)
      double precision,allocatable :: AP(:,:),PP(:)
      double precision  :: charg(maxcen),xyz(maxcen,3)

      integer :: Max_Cn = 3
      double precision :: FreqOm = 0.d0

      character(:), allocatable :: TrexFile

end type SystemBlock

type FileNames

     character(:),allocatable :: rdmfile
     character(:),allocatable :: onefile
     character(:),allocatable :: sirifile,siriusfile
     character(:),allocatable :: occfile,coefile
     character(:),allocatable :: testfile
     character(:),allocatable :: twofile,twojfile,twokfile
     character(:),allocatable :: y01file,xy0file
     character(:),allocatable :: abfile
     character(:),allocatable :: propfile
     character(:),allocatable :: propfile0,propfile1,propfile2

end type FileNames

type FlagsData
! default setting: ERPA-GVB
     ! mainp.f
     integer :: IDALTON = 1
     integer :: iORCA   = 0
     integer :: iTREXIO = 0
     integer :: IRes    = 0
     integer :: IAO     = 0
     integer :: INO     = 0
     integer :: NoSym   = 1
     integer :: NoSt    = 1
     integer :: IGVB    = 1
     integer :: ITwoEl    = TWOMO_INCORE 
     integer :: IRedVirt  = FLAG_REDVIRT
     integer :: IRDMCorr  = FLAG_RDMCORR
     integer :: IRDM2Typ  = FLAG_RDM2TYP
     integer :: ICholesky = FLAG_CHOLESKY
     integer :: ICholeskyAccu = CHOL_ACCU_DEFAULT
     integer :: IFun      = 13
     integer :: IFunSR    = 0 
     integer :: IFunSRKer = 0
     double precision :: Alpha = 0
     integer :: IModG   = 1
     integer :: NGOcc   = 0
     integer :: ILoc    = 1
     integer :: IFreeze = 0
     integer :: IAPSG   = 1
     integer :: ISERPA  = 0
     integer :: ITrpl   = 0
     integer :: ISAPT   = 0
     integer :: SaptLevel = 0
     integer :: ISHF      = 0
     character(:), allocatable :: JobTitle
     integer :: JobType = 0
     integer :: MemVal  = 2
     integer :: MemType = 3
     ! initia.f
     integer :: IA        = 1
     integer :: ICASSCF   = 0
     integer :: IDMRG     = 0
     integer :: ICI       = 0
     ! interpa.f  
     integer :: IFlAC     = 0
     integer :: IFlSnd    = 0
     integer :: IFlAC0D   = 0
     integer :: IFlAC0DP  = 0
     integer :: IFlACFREQ = 0
     integer :: IFlACFREQNTH = 0
     integer :: IFlAC1FREQNTH = 0
     integer :: IFlRESPONSE = 0
     integer :: ISymmAC0D = 1
     integer :: IFlCore   = 1
     integer :: IFlFrag1  = 0
     integer :: IFl12     = 1
     ! sapt_main.f90
     integer :: IFlag0 = 0
     ! sapt_utils.f90
     integer :: DIISN  = 6
     integer :: DIISOn = 2

end type FlagsData

type InputData

     type(CalculationBlock) :: CalcParams
     type(SystemBlock),allocatable :: SystemInput(:)
     integer :: iflag = 0
     type(FlagsData) :: Flags 

end type InputData

type SaptData

  type(SystemBlock) :: monA,monB
     double precision  :: Vnn,elst,exchs2,e2ind,e2disp
     double precision  :: elALL(3)
     double precision  :: e2disp_sc,e2disp_sp
     double precision  :: e2ind_unc,e2disp_unc
     double precision  :: e2dispR_unc,e2dispR
     double precision  :: e2exdisp_unc,e2exdisp
     double precision  :: e2exdispR_unc,e2exdispR
     double precision  :: e2exind,e2exind_unc
     double precision  :: e2ind_a0,e2ind_a1,e2ind_a2
     double precision  :: e2disp_a0,e2disp_a1,e2disp_a2
     double precision  :: e2exd_a0,e2exd_a1,e2exd_a2
     double precision  :: e2exi_a0,e2exi_a1,e2exi_a2
     double precision  :: esapt2,esapt0
     double precision  :: e2dispinCAS
     double precision  :: exch_part(5)
     double precision,allocatable :: Wind(:),Wdisp(:)
     ! test Pmat
     double precision,allocatable :: CholVecs(:,:)
     integer :: InterfaceType = INTER_TYPE_DAL
     integer :: SaptLevel = SAPTLEVEL2
     integer :: NAO
     integer :: NCholesky
     integer :: Max_Cn = 4
     integer :: ic6 = 0
     integer :: iPINO=-1
     integer :: IPrint = 1000
     logical :: iCpld  = .true.
     logical :: Cubic  = .false.
     logical :: CAlpha = .false.
     ! MH : add keyword!
     logical :: SemiCoupled = .true.
     logical :: Wexcit  = .false.
     logical :: noE2exi = .true.
     logical :: EnChck  = .true., HFCheck=.true.
     logical :: doRSH   = .false., SameOm = .true.
     logical :: reduceV = .false.
     character(:),allocatable :: TrexFile

end type SaptData

type SaptDIIS
     integer :: DIISN  = 6
     integer :: DIISOn = 2
     integer :: maxIter = 20
     double precision :: Thresh=1d-3
end type

contains 

function iaddr(IAddr1,IAddr2,IAddr3,IAddr4) result(NAddr3)
! POINTER FOR TWO-ELECTRON INTEGRALS
 implicit none
! parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
 integer    :: IAddr1, IAddr2, IAddr3, IAddr4
 integer(8) :: IAddr12, IAddr34
 integer(8) :: NAddr3

! CHANGE THE ORDER IF NECESSARY

 IAddr12 = Max(IAddr1,IAddr2)*(Max(IAddr1,IAddr2)-1)/2 + &
          Min(IAddr2,IAddr1)
 IAddr34 = Max(IAddr3,IAddr4)*(Max(IAddr3,IAddr4)-1)/2 + &
          Min(IAddr3,IAddr4)

! GET THE POSITION OF THE ELEMEMT (12|34)

 NAddr3 = Max(IAddr12,IAddr34)*(Max(IAddr12,IAddr34)-1)/2 + &
          Min(IAddr12,IAddr34)

end function iaddr

end module types
