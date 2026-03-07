module globals
! Module version of commons.inc
! Declares all common block variables with explicit types.
! Common blocks are retained for backward compatibility with
! code that still uses 'Include commons.inc'.
!
! New code should use 'use globals' instead of Include.
!
implicit none
save

! --- MISCFL ---
double precision :: XELE
integer :: IFun, IFunSR, IFunSRKer, NELE, IRes, INO, IPrint
integer :: IInteg, IAO, IGVB, IFreeze, ILoc
integer :: NInActOrb, NActOrb, NVirtOrb, NCoreOrb
integer :: NStronglyOccOrb, NElecBEmb
integer :: ISAPT, IWarn, ITwoEl, ITrpl
integer :: ICholesky, ICholeskyBIN, ICholeskyOTF, ICholeskyTHC
integer :: ICholeskyAccu, Monomer, MemVal, MemType
integer :: IH0test, IUnits, InternalGrid, IGridType

common /MISCFL/ XELE, IFun, IFunSR, IFunSRKer, NELE, IRes, INO, IPrint, &
     IInteg, IAO, IGVB, IFreeze, ILoc, NInActOrb, NActOrb, NVirtOrb, NCoreOrb, &
     NStronglyOccOrb, NElecBEmb, ISAPT, IWarn, ITwoEl, ITrpl, &
     ICholesky, ICholeskyBIN, ICholeskyOTF, ICholeskyTHC, &
     ICholeskyAccu, Monomer, MemVal, MemType, &
     IH0test, IUnits, InternalGrid, IGridType

! --- INTF ---
integer :: IDALTON, IMOLPRO, ITREXIO, IPYSCF
integer :: IOrbOrder, IJobType

common /INTF/ IDALTON, IMOLPRO, ITREXIO, IPYSCF, &
     IOrbOrder, IJobType

! --- CAS ---
integer :: ICASSCF, NAcCAS, NInAcCAS, NStates
integer :: InSt(2,100), InTrSt(2,1), ISpinMs2

common /CAS/ ICASSCF, NAcCAS, NInAcCAS, NStates, &
     InSt, InTrSt, ISpinMs2

! --- BBC3 ---
integer :: IType(5000)

common /BBC3/ IType

! --- DFTSR ---
double precision :: Alpha, Cmix
integer :: IFunSR2

common /DFTSR/ Alpha, Cmix, IFunSR2

! --- Sym ---
character(60) :: FMultTab
integer :: NoSym, MxSym

common /Sym/ FMultTab, NoSym, MxSym

! --- APSG ---
integer :: IGem(5000)
double precision :: CICoef(5000)

common /APSG/ IGem, CICoef

! --- DALTON ---
integer :: NISHT_G, NASHT_G, ISAPSG

common /DALTON/ NISHT_G, NASHT_G, ISAPSG

! --- EMBEDD ---
integer :: IFrag, NGemSave
integer :: IAuxGem(5000), IGemSave(5000), IOrbA(5000)
integer :: NGemB, NActive
integer :: IConnect(5000,5000)

common /EMBEDD/ IFrag, NGemSave, IAuxGem, &
     IGemSave, IOrbA, NGemB, NActive, IConnect

! --- AC ---
integer :: IFlAC, IFlSnd, IFlAC0D, IFlAC0DP, IFlCore, IFlFrag1, IFl12
integer :: NoSt, ISymmAC0D, IFlACFREQ, IFlACFREQNTH, IFlAC1FREQNTH, Max_Cn
integer :: IFlRESPONSE, IFlCorrMD, IRedVirt, IOrbRelax, IOrbIncl
double precision :: FreqOm(10)
integer :: IFlFCorr, IDBBSC, IHNO1, IVEMB, NFreqOm

common /AC/ IFlAC, IFlSnd, IFlAC0D, IFlAC0DP, IFlCore, IFlFrag1, IFl12, &
     NoSt, ISymmAC0D, IFlACFREQ, IFlACFREQNTH, IFlAC1FREQNTH, Max_Cn, &
     IFlRESPONSE, IFlCorrMD, IRedVirt, IOrbRelax, IOrbIncl, FreqOm, &
     IFlFCorr, IDBBSC, IHNO1, IVEMB, NFreqOm

! --- THRGEM ---
double precision :: ThrGemAct

common /THRGEM/ ThrGemAct

! --- THR ---
double precision :: ThrSelAct, ThrAct, ThrQVirt, ThrQInact

common /THR/ ThrSelAct, ThrAct, ThrQVirt, ThrQInact

end module globals
