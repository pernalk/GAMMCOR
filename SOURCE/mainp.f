C     A GENERAL CODE FOR CARRYING OUT LR-DMFT/SR-DFT CALCULATIONS
C
C                 PROGRAM INTERFACED WITH OUTPUT FILES FROM MOLPRO DEVELOP VERSION
C                 ONE- AND TWO-ELECTRON INTEGRALS (IN HF MO) ARE OBTAINED
C                 FROM ATMOL AND THEY ARE USED AS A GUESS WITH U=1
C                 GAUSSIAN BASIS SET USED
C
C
C     K.PERNAL 2018
C
      Program PRDMFT
C
      use types
      use inputfill
      use systemdef
      use sapt_main
      use timing
      use read_external
      use git_info
      use build_info
      use omp_lib
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
C
      Real*8, Dimension(:), Allocatable :: Occ
      Real*8, Dimension(:), Allocatable :: URe
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: XKin
      Real*8, Dimension(:), Allocatable :: XNuc
      Real*8, Dimension(:), Allocatable :: TwoEl
      Real*8, Dimension(:), Allocatable :: TwoElErf
      Real*8, Dimension(:), Allocatable :: DipX
      Real*8, Dimension(:), Allocatable :: DipY
      Real*8, Dimension(:), Allocatable :: DipZ
      Real*8, Dimension(:), Allocatable :: EpsHF 
      Real*8, Dimension(:), Allocatable :: UMOAO
      Integer, Dimension(:), Allocatable :: NSymMO
C
      type(InputData) :: Input
      type(FlagsData) :: Flags
      type(SystemBlock) :: System
      type(SaptData) :: Sapt
C
      Include 'commons.inc'
C
C     COMMON BLOCK USED IN OPTCPMFT
C
      Common/CPMFT/ MFrac,MOcc,NFrac
C
C     *************************************************************************
C
      Call git_print_info()
      Call build_print_info()

      Call omp_set_max_active_levels(2)

      Call read_Input(Input)
      Call check_Calc(Input%CalcParams)
      Call fill_Flags(Input,Flags)
      Call create_System(Input,Flags,System,Sapt)
C
      Call free_Input(Input)
C
C     FILL COMMONS AND CONSTANTS
      XELE    = System%XELE
      NELE    = System%NELE
      Charge  = System%Charge
      NBasis  = System%NBasis
      Title   = Flags%JobTitle
      ITwoEl  = Flags%ITwoEl
      ICholesky = Flags%ICholesky
      ICholeskyAccu = Flags%ICholeskyAccu
      MemVal  = Flags%MemVal
      MemType = Flags%MemType
      IWarn   = 0
      Max_Cn  = System%Max_Cn
      ITrpl   = Flags%ITrpl
      FreqOm  = System%FreqOm
      IRedVirt=Flags%IRedVirt
c      ThrVirt=System%ThrVirt
C
C     *************************************************************************
C
C     ICASSCF = 1 - CORRELATION ENE FOR CAS WILL BE COMPUTED
C             = 0 - CORRELATION ENE FOR GVB WILL BE COMPUTED
C
      ICASSCF=Flags%ICASSCF

C     IF IDALTON=1 - READ INTEGRALS AND C COEF, IGem FROM A FILE GENERATED BY DALTON
C
      IDALTON=Flags%IDALTON
      IDMRG=Flags%IDMRG
      ITREXIO=Flags%ITREXIO
C
C     IF IRes=1 - RESTART THE CALCULATIONS FROM A RESTART FILE
C
      IRes=Flags%IRes
C
C     SET IAO TO ONE IF INTEGRALS IN AO ARE OBTAINED FROM MOLPRO
C
      IAO=Flags%IAO
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
      IFlAC=Flags%IFlAC
C
C     IFlACFREQ = 1 - adiabatic connection formula calculation with frequency integration
C                 0 - AC without frequency integration
      IFlACFREQ=Flags%IFlACFREQ
C
C     IFlACFREQNTH = 1 - AC integrand expanded wrt alpha, frequency integration
C                    0 - AC without frequency integration
      IFlACFREQNTH=Flags%IFlACFREQNTH
C
C     IFlAC1FREQNTH = 1 - AC1 via expansion wrt alpha, frequency integration
      IFlAC1FREQNTH=Flags%IFlAC1FREQNTH
C
C
C     IFlSnd  = 1 - run AC0 (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0
      IFlSnd=Flags%IFlSnd
C
C     IFlAC0D = 1 - run AC0D ()
C             = 0 - do not run AC0D
      IFlAC0D=Flags%IFlAC0D
C
C     ISymmAC0D = 1 : AC0D corrections will be computed based on symmetry (dafault)
C                 0 : AC0D corrections will be computed based on overlap with SA-CAS
C                     states (run molpro with symmetry,nosym)
C
      ISymmAC0D=Flags%ISymmAC0D
C
C     IFlAC0DP  = 1 : AC0D prime corrections will be computed based on symmetry (dafault)
C                 0 : AC0D without prime 
C
      IFlAC0DP=Flags%IFlAC0DP
C
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
      IFlCore=Flags%IFlCore
C
      IFlFrag1=Flags%IFlFrag1
      IFl12=Flags%IFl12
C
C     NoSym = 1 - DO NOT IMPOSE SYMMETRY (RUN molpro WITH 'nosym')
C           = 0 - IMPOSE SYMMETRY (RUN molpro WITHOUT 'nosym')
C
      NoSym=Flags%NoSym
C
C     IFlRESPONSE = 1 : compute polarizability tensor for a given frequency FreqOm
C
      IFlRESPONSE=Flags%IFlRESPONSE
C
C     *************************************************************************
C
C     SELECT A LONG-RANGE DMFT FUNCTIONAL
C
C     IFun   = 0 - NO DMFT FUNCTIONAL (PURE DFT CALCULATIONS, EQUIVALENT TO USING mu=0)
C
C     IFun   = 1 - USE THE KUTZELNIGG FUNCTIONAL:
C                  Eee = sqrt(ni nj) <ij|ji> + sqrt(na nb) <ab|ba>
C                      - 2 sqrt(ni na) <ia|ai>
C
C     IFun   = 2 - USE THE BUIJSE-BAERENDS FUNCTIONAL:
C                  Eee = np nq <pq|pq> - sqrt(np nq) <pq|qp>
C
C     IFun   =20 - BB WITH A "+" SIGN FOR THE BONDING-ANTIBONDING PAIR
C
C     IFun   = 3 - USE THE GOEDECKER-UMRIGAR FUNCTIONAL:
C                  Eee = np nq <pq|pq> - sqrt(np nq) <pq|qp>
C                      + (np-np**2) <pp|pp>
C
C     IFun   = 4 - USE THE BBC1 FUNCTIONAL
C
C     IFun   = 5 - USE THE BBC2 FUNCTIONAL
C
C     IFun   = 6 - USE THE BBC3 FUNCTIONAL
C
C     IFun   = 7 - CORRECTED HARTREE-FOCK FUNCTIONAL
C                  Eee = np nq <pq|pq>
C                      - [(np nq) + sqrt(np(1-np)nq(1-nq)) ]<pq|qp>
C     IFun   = 8 - BB AND HF HYBRID FUNCTIONAL
C                  (1-Cmix)*Eee_BB + Cmix*Eee_HF
C     IFun   = 10 - LR-BB+SR-HF+SR-PBE_CORR
C                  Eee_LRBB + Eee_SRPBE
C     IFun   = 11 - CPMFT method of Scuseria et al.
C                   (define the active space: NActive, MActive
C                    - the number of active elecrons and orbitals)
C
C     IFun   = 12 - Hartree-Fock
C
C     IFun   = 13 - APSG
C
C     IGVB   = 1  - APSG WITH ONLY TWO ORBITALS PER GEMINAL
C
      IGVB=Flags%IGVB
C
      IFun=Flags%IFun
C
C     *************************************************************************
C
C     SELECT A SHORT-RANGE DFT FUNCTIONAL
C
C     IFunSR = 0 - do not include a short-range functional
C     IFunSR = 1 - SR-LSDA, Paziani et al.
C     IFunSR = 2 - SR-PBE, Goll et al. PCCP 7, (2005) 3917
C     IFunSR = 3 - Gagliardi-Truhlar with PBE
C     IFunSR = 4 - LR-CAS+SR-PBE+LR-AC/AC0 (with full range CAS RDM's)
C
C     OTHERS
C
C     IFunSR = 5 - CASPiDFT (CASPIDFT procedure)
C     IFunSR = 6 - CASPiDFT (CASPIDFTOPT, integrals not loaded)
C     IFunSR = 7 - VV10 (nonlocal correlation functional, integrals not loaded)
C
      IFunSR=Flags%IFunSR
      IFunSRKer=Flags%IFunSRKer
C
      If(IFunSRKer.Eq.2)Stop'Fatal Error: Kernel for PBE not available!'
      Write(*,'(/,1x,"IFunSR=",I1)')IFunSR
      Write(*,'(/,1x,"IFunSRKer=",I1,/)')IFunSRKer
C
C     *************************************************************************
C
C     CHECK IF SAPT RUN
      ISAPT=Flags%ISAPT
C
      If(ISAPT.Eq.1) Call sapt_driver(Flags,Sapt)
C
C     NBasis READ FROM SIRIUS.RST
      If(IDALTON.Eq.1) Call basinfo(NBasis,'SIRIUS.RST','DALTON')
CC
C     *************************************************************************
C
C     SELECT ELECTRONIC STATE
C     (maybe move selection to LdInteg and ReadDAL?)
C
      NStates = System%NStates
      ISpinMs2= System%ISpinMs2
      InSt(1:2,1:NStates) = System%InSt
      InTrSt(1:2,1:1) = System%InTrSt
C
      If(InSt(2,1).Gt.0) Then
      NoSt = System%InSt(1,1)
C
C      Print*,'VALUE DECLARED IN INPUT: ',NoSt
C
      ElseIf(IDALTON.Eq.0.and.IDMRG.Eq.0) Then
      Call read_NoSt_molpro(NoSt,'2RDM')
      ElseIf(IDALTON.Eq.1) Then
      NoSt = 1
      Write(6,'(/,1x,a)') 'WARNING! ASSUMING RMDs CORRESPOND TO
     $                     GROUND STATE FROM DALTON OUTPUT (NoSt=1)!'
C
      EndIf
C
      Write(6,'(1x,"NoSt=",I2)')NoSt
C
C     SET THRESHOLD FOR ACTIVE ORBITALS IN CAS
      ThrSelAct = System%ThrSelAct
C     SET THRESHOLD FOR ACTIVE ORBITALS IN GVB
      ThrGemAct = System%ThrGemAct
C     SET THRESHOLD FOR ACTIVE ORBITALS IN CIPSI
      ThrAct    = System%ThrAct
C     SET THRESHOLD FOR QUASI-VIRTUAL ORBITALS IN CAS
      ThrQVirt  = System%ThrQVirt
C     SET THRESHOLD FOR QUASI-INACTIVE ORBITALS IN CAS
      ThrQInact = System%ThrQInact
C
C*************************************************************************
C     READ THE INPUT AND PRINT THE INPUT DATA
C
C     OLD INPUT-READ
C      Call RWInput(Title,ZNucl,Charge,NBasis)
C
C     CALCULATE THE DIMENSIONS
      If(IDALTON.Eq.0) then
C        Call CheckNBa(NBasis,Title)
        Call basinfo(NBasis,'AOONEINT.mol','MOLPRO')
      endif
C
      Call DimSym(NBasis,NInte1,NInte2,MxHVec,MaxXV)
C
      If(IFunSR.Ne.0) Then
C      Call GetNGrid(NGrid,Title)
      Call molprogrid0(NGrid,NBasis)
      Else
      NGrid=1
      EndIf
C
C     GET THE VALUE OF THE SEPARATION PARAMETER OM
C
      Alpha = System%Omega
c      If(IFunSR.Ne.0.And.IFunSR.Ne.3.And.IFunSR.Ne.5) Then
      If(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
C      Call GetAlpha(Title)
      Call readalphamolpro(Alpha)
      Else
      Alpha=0.D0
      EndIf
C
C     ALLOCATE THE MATRICES
C
C     FOR TESTS SWITCHIG IT OFF...
      If(ITwoEl.Eq.3) NInte2=1
      If(ITwoEl.Eq.2) NInte2=1
C
      Allocate  (Occ(NBasis))
      Allocate  (URe(NBasis*NBasis))
      Allocate  (XKin(NInte1))
      Allocate  (XNuc(NInte1))
      Allocate  (TwoEl(NInte2))
      Allocate  (UMOAO(NBasis*NBasis))
C
      Occ(1:NBasis)=           0.D0
      URe(1:NBasis*NBasis)=    0.D0
      XKin(1:NInte1)=          0.D0
      XNuc(1:NInte1)=          0.D0
      TwoEl(1:NInte2)=         0.D0
      UMOAO(1:NBasis*NBasis)=  0.D0
C
C     LOAD THE GRID AND THE VALUES OF THE MO ORBITALS ON THE GRID
C
      write(LOUT,'()')
      write(LOUT,'(1x,a)') 'STARTING CALCULATIONS '
      write(LOUT,'(8a10)') ('**********',i=1,8)
C
      Call clock('START',Tcpu,Twall)
C
C     LOAD THE INTEGRALS
C
      If(IDALTON.Eq.0) Then
C
      IFFSR=0
      If(IFunSR.Eq.7) Then
C     temporarily set IFunSR to 6 to avoid loading integrals and their transformation in LdInteg
      IFunSR=6
      IFFSR=7
      EndIf
C
      Call LdInteg(Title,XKin,XNuc,ENuc,Occ,URe,TwoEl,UMOAO,NInte1,
     $ NBasis,NInte2,NGem)
C
C     create cas_ss.molden file with NOs (useful after SA-CAS calculations to inspect 
C     the character of NOs in a state requested in input.inp)
C
      Call MoldenCAS(Occ,UMOAO,NBasis)
C
C     set back IFunSR to IFFSR
      If(IFFSR.Ne.0) IFunSR=IFFSR
C
      Else
C
      Call ReadDAL(XKin,XNuc,ENuc,Occ,URe,TwoEl,UMOAO,
     $ NInte1,NBasis,NInte2,NGem,Flags)
C
      EndIf
C
      If(IFunSR.Eq.5) Then
      Call CASPIDFT(ENuc,URe,UMOAO,Occ,XKin,TwoEl,
     $ NBasis,NInte1,NInte2)
      ElseIf(IFunSR.Eq.6) Then
      Call CASPIDFTOPT(URe,UMOAO,Occ,NBasis)
      ElseIf(IFunSR.Eq.7) Then
      Call VV10(URe,UMOAO,Occ,NBasis)
      Else
      Call DMSCF(Title,URe,Occ,XKin,XNuc,ENuc,UMOAO,
     $ TwoEl,NBasis,NInte1,NInte2,NGem)
      EndIf
C
      If(IWarn.Gt.0) Then
      Write(6,'(/,1x,a,i2,1x,a)') 'CHECK OUTPUT FOR',IWarn,'WARNINGS!'
C      Write(6,'(8a10)') ('**********',i=1,9)
      EndIf
C
      If(ITwoEl.Eq.2)  Call delfile('TWOMO')
      Call delfile('AOTWOSORT')
C
      Call free_System(System)
      Call clock(PossibleJobType(Flags%JobType),Tcpu,Twall)
C      Stop
      End
C
