*Deck DMSCF 
      Subroutine DMSCF 
     $ (ETot,Title,URe,Occ,Eps,XKin,XNuc,ENuc,UMOAO,
     $ DipX,DipY,DipZ,TwoEl,TwoElErf,NBasis,NInte1,NInte2,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ IAPSG,ISERPA,QMAX,IModG,NGOcc,Small,NGem,
     $ Flags)
C
C     !!! XKin CONTAINS BOTH KINETIC AND EL-N CONTRIBUTIONS !!!
C     !!! XNuc IS EMPTY 
C
C     CALLS AN OPTIMIZER FOR THE ORBITALS AND OCCUPANCIES
C
      use types 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FMultTab
C
      Include 'commons.inc'
      Common/CPMFT/ MFrac,MOcc,NFrac
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Parameter (MXIT=50, ETol=1.D-5, DampG=0.2D0)
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),Eps(NBasis),XKin(NInte1),
     $ XNuc(NInte1),DipX(*),DipY(*),DipZ(*),TwoEl(NInte2),TwoElErf(*),
     $ OrbGrid(NBasis,NGrid),
     $ WGrid(NGrid),UMOAO(Nbasis,NBasis),NSymMO(NBasis)
C
      Dimension OrbXGrid(*),OrbYGrid(*),OrbZGrid(*)
C
C     LOCAL ARRAYS
C
      Dimension VSR(NInte1),Gamma(NInte1),
     $ XOne(NInte1),UNOAO(Nbasis,NBasis),GammaOld(NInte1),Work(NBasis),
     $ NSymNO(NBasis),OccSav(NBasis),UReSav(NBasis,NBasis),
     $ IGemSav(NBasis),Work2(NBasis,NBasis)
C
      Real*8, Allocatable :: TNO(:),TNOLR(:)
      type(FlagsData) :: Flags

C
C     FLAGS:  IRes   = 0 ... run the program from the beginning
C                      1 ... restart the program 
C
C     IF IAO=1 ADD XNuc TO XKin AND EMPTY XNuc
C
      If(IAO.Eq.1) Then
      Do I=1,NInte1
      XKin(I)=XKin(I)+XNuc(I)
      XNuc(I)=Zero
      EndDo
      EndIf
C
C     RESTART CALCULATIONS 
C
      If (IRes.Eq.1) Then
C
      If(IFun.Eq.13) Then
C
      Call ReWrAPSG(0,Occ,URe,UReSav,Title,NGem,NBasis)
C
C     check phases of the mo's and change phases of no's if needed
      Do IMO=1,NBasis
      IPhase=1
C
      Do IA=1,NBasis
      If(Abs(UReSav(IMO,IA)).Gt.Half) Then
      If(UReSav(IMO,IA)*UMOAO(IMO,IA).Lt.Zero) IPhase=-1
      EndIf
      EndDo
C
      Do IA=1,NBasis
      URe(IA,IMO)=IPhase*URe(IA,IMO)
      EndDo
C
      EndDo
C
      If(IGVB.Eq.0) Then
      Write(6,'('' ***  RESTARTED APSG CALCULATION  *** '',
     $      25X,''*'')')
      Else
      Write(6,'('' ***  RESTARTED GVB  CALCULATION  *** '',
     $      25X,''*'')')
      EndIf
C   
      Else
C
      Call ReWr(0,Occ,URe,Title,NBasis) 
C
      If(IPrint.Ge.1) Then
C
      Write(6,'('' ***  RESTARTED CALCULATION  *** '',
     $      30X,''*'')')
      Write(6,'('' * DENSITY MATRIX FUNCTIONAL CALCULATION'',
     $      23X,''*'')') 
      Write(6,
     $ '('' * Functional : '',I3 , 44X,''*'')') IFun
      Write(6, '('' **********************************************'',
     $            ''*****************'')') 
C
      EndIf
C
C     CHECK THE NORM
C
      Sum=Zero
      Do 5 I=1,NBasis
    5 Sum=Sum+Occ(I)
      If(Abs(Sum-XELE).Gt.1.D-10) Call NormN(Occ,NBasis,IFAIL)
      If(IFAIL.Eq.1) 
     $ Stop'Occ retrieved from a restart file cannot be normalized!'
      Sum=Zero
      Do 6 I=1,NBasis
    6 Sum=Sum+Occ(I)
C
C     EndIf of IFun.Eq.13)
      EndIf
C
C     EndIf of IRes.Eq.1
      EndIf
C
      If(IDALTON.Eq.1) Then
      Do I=1,NInte1
      XOne(I)=XKin(I)
      EndDo
      EnSR=Zero
      Do I=1,NInte1
      VSR(I)=Zero
      EndDo
      XVSR=Zero
      GoTo 333
      EndIf
C
      If(IRes.Eq.0) Then
C
C     GENERATE A GUESS FOR THE APSG 
C
      Do I=1,NBasis
      IGem(I)=1
      EndDo
C
      If(IFun.Eq.13) Then
C
      If(INO.Eq.0) Then
C
C     SORT ORBITALS ACC. TO THE HF ORBITAL ENERGIES
C
      Do I=1,NBasis
      NSymNO(I)=NSymMO(I)     
      EndDo
      If(ILoc.Eq.0)
     $ Call HF(URe,Eps,Occ,XKin,XNuc,TwoEl,NSymNO,NInte1,NInte2,NBasis,
     $ NELE,IPrint,1)
C
      NGem=NELE
      ICount=0
      Do I=1,NBasis
C
      CICoef(I)=Sqrt(Occ(I))
C
      ICount=ICount+1
      IGem(I)=ICount
      If(ICount.Eq.NGem) ICount=0
C
      EndDo
C
C     IF IGVB=1 ASSIGN ORBITALS TO GEMINALS BY HAND
C
c      If(IGVB.Eq.1) Then
c      NGem=NGem+1
c      Do I=3,NBasis
c      IGem(I)=NGem
c      EndDo
c      EndIf
C
C     DISTRIBUTE ORBITALS AMONG NGOcc+1,...,NGem GEMINALS
C
      If(NGOcc.Ne.0) Then
C
C     IF NGOcc=NELE ASSIGN VIRTUAL ORBITALS TO THE LAST GEMINAL
C
C KP 24.11.14
C     IF NGOcc=NELE -> ASSIGN ALL VIRTUALS TO A FICTITIOUS GEMINAL NELE+1
C
      If(NGOcc.Eq.NELE) NGem=NELE+1
C
      ICount=NGOcc
      Do I=NGOcc+1,NBasis
      ICount=ICount+1
      IGem(I)=ICount
      If(ICount.Eq.NGem) ICount=NGOcc
      EndDo
C
c      Call Partition(UMOAO,NGOcc,NGem,NBasis)
C
      Write (*,'(/,X,
     $ "The number of orbitals with occupancy fixed to 1 is ",I1,/)')
     $ NGOcc 
C
      EndIf
C
      GoTo 91
C
C     else of INO.Eq.0
      Else
C
      Do I=1,NBasis
      NSymNO(I)=0
      Do J=1,NBasis
      If(Abs(URe(I,J)).Gt.0.1) Then
      ISym=NSymMO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym)
     $ Write(6,'("Symm of NO cannot be established",I3)')I
      EndIf
      EndIf
      EndDo
      Write(6,'(I4,7X,I4,2X,E16.6)')I,NSymNO(I),Occ(I)
      EndDo
C
      NGem=NELE
      ICount=0
      Do I=1,NBasis
      CICoef(I)=Sqrt(Occ(I))
      ICount=ICount+1
      IGem(I)=ICount
      If(ICount.Eq.NGem) ICount=0
      EndDo
C
c     endif INO=0
      EndIf 
C
   91 Write(6,'(2X,''NO OF GEMINALS: '',I3)') NGem
      Write(6,'(2X,"Orb",2X,"Occupancy  HF Orb Energies",
     $ X,"Gem",3X,"Sym")')
      Do I=1,NBasis
      Write(6,'(X,I3,E12.2,E14.4,2I6)')I,Occ(I),Eps(I),IGem(I),NSymNO(I)
      EndDo
C
c      If (IFreeze.Eq.1) Then
C
      Call MultpM(UReSav,URe,UMOAO,NBasis)
      NLine=NBasis/10
      If(NLine*10-NBasis.Ne.0)NLine=NLine+1
      Do I=1,NBasis
      Write(*,'(I3)') I
C
      Do LL=0,NLine-1

      NN=NBasis-10*LL
      If(NN.Le.10) Then
      Write(*,98) (UReSav(I,K),K=10*LL+1,NBasis)
      Else
      Write(*,98) (UReSav(I,K),K=10*LL+1,10*(LL+1))
      EndIf
C
      EndDo
      Write(*,*)
      EndDo
C
c      EndIf 
C
c     endif Ifun=13
      EndIf
C
c     endif of IRes.Eq.0
      EndIf
CCCC
C     BEGINNING OF THE ITERATION OF DENSITY/DENSITY MATRIX SCHEME 
CCCC
C
      EOld=Zero
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
      Do II=1,MXIT
C
C     IF IFunSR=1 AND IFun=13 SAVE THE Occ AND URe 
C
      If(IFunSR.Ne.0.And.IFun.Eq.13) Then
      Call CpyM(UReSav,URe,NBasis)
      Call CpyV(OccSav,Occ,NBasis)
      Do I=1,NBasis
      IGemSav(I)=IGem(I)
      EndDo
      EndIf
C
C     AVERAGE THE DENSITY MATRIX
C
      If(II.Gt.1) Then
C
      Do L=1,NInte1
      GammaOld(L)=Gamma(L)
      EndDo
C 
      Call VecTr(Gamma,Occ,URe,NBasis)
C
      Do L=1,NInte1
      Gamma(L)=DampG*Gamma(L)+(One-DampG)*GammaOld(L)
      EndDo
C
C     RESTORE SYMMETRY
C
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
      If(NSymMO(K).Ne.NSymMO(L)) Gamma(KL)=Zero
      EndDo
      EndDo  
C
      Call CpySym(URe,Gamma,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,Work)
C     sometimes very small occ's come out negative so zero them
      Do I=1,NBasis
      Occ(I)=Abs(Occ(I))
      EndDo
C    
      Call SortOcc(Occ,URe,NBasis)
C
      If(MFrac.Ne.NBasis) Then
      Do I=1,NBasis
      If(I.Le.MOcc) Occ(I)=One
      If(I.Gt.MOcc+MFrac) Occ(I)=Zero
      EndDo
      Call NormN(Occ,NBasis,IFAIL)
      EndIf
C
c     endif of If(II.Gt.1) 
      EndIf   
C
C     COMPUTE THE SR POTENTIAL (H+XC) IN THE MO REP FOR THE DENSITY U*Occ*U+
C
      If(IFunSR.Eq.0) Then
C
      EnSR=Zero
      Do I=1,NInte1
      VSR(I)=Zero
      EndDo
      XVSR=Zero
C
      Else 
C
      Call EPotSR(EnSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,NSymMO,TwoEl,TwoElErf,NGrid,NInte1,NInte2,NBasis)
C
C     COMPUTE THE EXPECTATION VALUE OF THE VSR POTENTIAL FOR Gamma
C
      XVSR=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Fac=Four
      If(I.Eq.J) Fac=Two
      XVSR=XVSR+Fac*Gamma(IJ)*VSR(IJ)
      EndDo
      EndDo
C
c     endif of If(IFunSR.Eq.0)
      EndIf
C
C     BEGIN THE ITERATION WITH THE SAVED Occ AND URe
C
      If(IFunSR.Ne.0.And.IFun.Eq.13) Then
      Call CpyM(URe,UReSav,NBasis)
      Call CpyV(Occ,OccSav,NBasis)
      Do I=1,NBasis
      IGem(I)=IGemSav(I)
      EndDo
      EndIf
C
C     ADD THE SR HXC POTENTIAL TO THE ONE-ELECTRON PART 
C
      Do I=1,NInte1
      XOne(I)=XKin(I)+XNuc(I)+VSR(I) 
      EndDo
C
      If(IFunSR.Ne.0) Then
C
      If(IFun.Eq.1) Then
C
c      Call OptTwo(ETot,URe,Occ,XOne,TwoElErf,NSymMO,
c     $ NBasis,NInte1,NInte2,1)
      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,
     $ TwoElErf,UMOAO,NSymMO,
     $ Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG)
C
      ElseIf(IFun.Eq.2.Or.IFun.Eq.8) Then
C
      Call OptGam2(ETot,Occ,URe,XOne,XNuc,TwoElErf,TwoEl,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ NInte1,NInte2,NBasis)
C
      Else
C
      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,TwoElErf,
     $ UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG) 
C
c     endif of If(IFun.Eq.1)
      EndIf
C
C     FOR IFunSR=0 
C
      Else
C
      If(IFun.Eq.1) Then
C
      Call OptTwo(ETot,URe,Occ,XOne,TwoEl,NSymMO,
     $ NBasis,NInte1,NInte2,1)
c      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,TwoEl,
c     $ UMOAO,NSymMO,
c     $ Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG)
C
      ElseIf(IFun.Eq.2) Then
C
      Call OptGam2(ETot,Occ,URe,XOne,XNuc,TwoEl,TwoEl,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ NInte1,NInte2,NBasis)
C
      Else
C
  333 Continue
      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,TwoEl,UMOAO,
     $ NSymMO,Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG)
C
c     endif of If(IFun.Eq.1)
      EndIf
C
c     endif of If(IFunSR.Ne.0)
      EndIf
C
      Call ReWr(1,Occ,URe,Title,NBasis)
C
      ENew=ETot+EnSR-XVSR
      If(IFunSR.Ne.0) Then
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
      Write(6,'(X,''DFT/DMFT ITER'',I3,2X,''ENERGY'',F14.8,2X,
     $ ''ENE DIFF '',E10.3)') II,Enew,ENew-EOld
      Write(6,'(X,''****************************************'',
     $            ''***************************************'',/)')
      EndIf
C
      Err=Abs((EOld-ENew))
C
      If(Err.Lt.ETol.Or.IFunSR.Eq.0) Then
C
C     PRINT THE TOTAL ENERGY (TOGETHER WITH THE NUCLEAR INTERACTION)
C
      If(IDALTON.Ne.1) 
     $ Write
     $ (6,'(/,10X,''CONVERGENCE ATTAINED WITH THE TOTAL ENERGY '',
     $ F14.8)')ENew+ENuc
C
      Call ReWrAPSG(1,Occ,URe,UMOAO,Title,NGem,NBasis)
C
C     PRINT MO's IN AO's REPRESENTATION 
C
c      Write(6,'(/,X,"MOLECULAR ORBITALS IN AO BASIS SET")')
C
   98 Format(X,10F10.6)
      NLine=NBasis/10
      If(NLine*10-NBasis.Ne.0)NLine=NLine+1
      Do I=1,NBasis
C
      Do LL=0,NLine-1

      NN=NBasis-10*LL
      If(NN.Le.10) Then
c      Write(*,98) (UMOAO(I,K),K=10*LL+1,NBasis)
      Else
c      Write(*,98) (UMOAO(I,K),K=10*LL+1,10*(LL+1))
      EndIf
C
      EndDo
c      Write(*,*)
      EndDo
C
C     PRINT FINAL NO's IN AO's REPRESENTATION 
C
      Write(6,'(/,X,"FINAL NATURAL ORBITALS IN AO BASIS SET")')
C
      Call MultpM(UReSav,URe,UMOAO,NBasis)
      NLine=NBasis/10
      If(NLine*10-NBasis.Ne.0)NLine=NLine+1
      Do I=1,NBasis
      Write(*,'(I3)') I      
C
      Do LL=0,NLine-1

      NN=NBasis-10*LL
      If(NN.Le.10) Then
      Write(*,98) (UReSav(I,K),K=10*LL+1,NBasis)
      Else
      Write(*,98) (UReSav(I,K),K=10*LL+1,10*(LL+1))
      EndIf
C
      EndDo
      Write(*,*)
      EndDo
C
c      Call PrintAPSG(Title,UMOAO,UReSav,URe,Occ,
c     $           XOne,TwoEl,NGem,NBasis,NInte1,NInte2)
C
C     PREPARE MOLDEN FILES TO VISUALIZE GEMINAL DENSITIES
C
      If(IDALTON.Eq.0) Call MoldenPrep(Title,Occ,UReSav,NGem,NBasis)
C
c     move on to compute excitations
      GoTo 444
C
c     endif of If(Err.Lt.ETol.Or.IFunSR.Eq.0)
      EndIf
C
      EOld=ENew
C
      EndDo
C
      Write(6,'(/,10X,''NO CONVERGENCE IN DFT/DMFT !!!'')')
C
  444 Continue
C
C     DYNAMIC PROPERTIES 
C
C     SORT OCCUPATION NUMBERS AND ORBITALS      
C
C DO NOT SORT ORBITALS IF FRAGMENT ASSIGNMENT IS READ FROM A FILE
c      Call SortOcc(Occ,URe,NBasis)
c      Write(6,'(/,X,''  ORBITAL OCCUPANCIES AFTER SORTING'')')
c
c the orbitals must be sorted if NActOrb is set different from zero since 
c the last "occupied" orbitals and the first "virtual" orbitals are taken as active 
      Write(6,'(/,X,'' UNSORTED ORBITAL OCCUPANCIES'')')

      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"CICoef",7X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,2E16.6,I6)') I,Occ(I),CICoef(I),IGem(I)
      EndDo
C
      If(IFunSR.Eq.0) Then
C
      EnSR=Zero
      Do I=1,NInte1
      VSR(I)=Zero
      EndDo
C
      Else
C
      Call EPotSR(EnSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,NSymMO,TwoEl,TwoElErf,NGrid,NInte1,NInte2,NBasis)
C
      EndIf
C
      Do I=1,NInte1
      XOne(I)=XKin(I)+XNuc(I)+VSR(I)
      EndDo
C
      NDim=NBasis*(NBasis-1)/2
      NDimKer=NBasis*(1+NBasis)*(2+NBasis)*(3+NBasis)/24
C
      If(IDALTON.Ne.1) Then
C
      Allocate (TNO(NInte2))
      Allocate (TNOLR(NInte2))
      TNOLR(1:NInte2)=Zero
      Call TwoNO(TNO,URe,TwoEl,NBasis,NInte2)
C
      If(IFunSR.Ne.0.And.IFun.Ne.0)
     $ Call TwoNO(TNOLR,URe,TwoElErf,NBasis,NInte2)
C
      IFunSR=IFunSRKer
C
      If(IFunSR.Ne.0) 
     $ Call ERPA(TNO,TNOLR,URe,Occ,XOne,
     $  OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,NDimKer,
     $  NGem,IAPSG,ISERPA,QMAX,Small)
C
      Deallocate (TNO)
      Deallocate (TNOLR)
C
c     If(IDALTON.Ne.1) Then
      EndIf
C
      If(IDALTON.Eq.1) Then
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
      IFlAC=Flags%IFlAC
C
C     IFlSnd  = 1 - run AC0 (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0
      IFlSnd=Flags%IFlSnd
C
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
      IFlCore=Flags%IFlCore
C
      IFlFrag1=Flags%IFlFrag
      IFl12=Flags%IFl12
C
      If(ICASSCF.Eq.1) Then
      Call ACCAS(ETot,ENuc,TwoEl,URe,UReSav,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NGem)
      Return
      EndIf
C
C     TwoEl - two-electron integrals in NO 
      ISERPA=0
      Call INTERPA(ETot,ENuc,TwoEl,URe,UReSav,Occ,XOne,Title,
     $  OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,
     $  NGem,IAPSG,ISERPA,QMAX,NGOcc,Small)
      EndIf
C
      Return
      End

*Deck Partition
      Subroutine Partition(ULoc,NGOcc,NGem,NBasis)
C     
      Implicit Real*8 (A-H,O-Z)
C     
C     DIVIDE A SET OF OCCUPIED LOCALIZED MO ORBITALS INTO ACTIVE (FRAGMENT A) AND INACTIVE ONES
C     VIRTUALS GO TO THE ACTIVE A, ASSIGN ORBITALS TO GEMINALS
C
      Include 'commons.inc'
C     
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C    
C     UMOAO - A TRANSFORMATION MATRIX FROM AO TO LOCALIZED MO
C 
      Dimension ULoc(NBasis,NBasis)
C    
C     LOCAL ARRAYS
C
C     NAt(I) - AN INDEX OF AN ATOMS ON WHICH THE ITH ORBITALS IS CENTERED
C     NBasisAt(I) - A NUMBER OF ATOMIC ORBITALS CENTERED ON THE ITH ATOM
C     NPair(I,1/2) - INDICES OF THE FIRST AND SECOND ATOMS ON WHICH SPECIES I IS LOCALIZED   
C     IOrbA(I) = 1 - ITH ORBITAL IS ACTIVE
C     IOrbA(I) = 0 - ITH ORBTIAL IS INACTIVE 
C     NoAt - THE NUMBER OF ATOMS IN A MOLECULE
C
      Dimension NAt(NBasis),NBasisAt(1000),NPair(1000,2)
C    
C     BEGINNING OF THE INPUT
C c herer!!! - change NBasisAt if basis set changes
      NoAt=3
      NBasisAt(1)=15
      NBasisAt(2)=5
      NBasisAt(3)=5
C
C     DEFINE AN ACTIVE SYSTEM AS A SET OF ACTIVE BONDS AND ATOMS
C     MAKE SURE THAT NPair(1,I) > NPair(2,I)
C
C     H2O : two bonds are active
C
      NoSpecies=2
C
C     OH1 BOND
      NPair(1,1)=2
      NPair(1,2)=1
C     OH2 bond
      NPair(2,1)=3
      NPair(2,2)=1
C
C     END OF INPUT
C
      Do I=1,NELE
      IOrbA(I)=0
      EndDo
C
C     ALL VIRTUALS WILL BE ACTIVE
C
      Do I=NELE+1,NBasis
      IOrbA(I)=1
      EndDo
C
C     CHECK IF THE NUMBER OF AO's SUM UP TO NBasis
C
      NB=0
      Do I=1,NoAt
      NB=NB+NBasisAt(I)
      EndDo
      If(NB.Ne.NBasis)
     $ Stop 'Fatal Error: No of atomic orbitals different from NBasis!'
C
C     SET A NUMBER OF AN ATOM ON WHICH AN AO IS CENTERED
C
      ICount=1
      Do I=1,NoAt
C
      Do J=1,NBasisAt(I)
      NAt(ICount)=I
      ICount=ICount+1
      EndDo
C
      EndDo
C
C     ASSIGN OCC ORBITALS TO AN ACTIVE OR INACTIVE SPACE
C
      NGemA=0
      Do I=1,NELE
C
      ISpec1=0
      ISpec2=0
C
      Do K=1,NoAt
      AtContr=Zero
      Do J=1,NBasis
      If(NAt(J).Eq.K) AtContr=AtContr+ULoc(I,J)**2
      EndDo
C
      If(AtContr.Gt.0.1D0) Then
      If(ISpec1.Eq.0) Then
      ISpec1=K
      ElseIf(ISpec2.Eq.0) Then
      ISpec2=ISpec1
      ISpec1=K
      Else
      Stop 'Problem with selecting active orbitals!'
      EndIf
      EndIf
C
      EndDo
C
      Do J=1,NoSpecies
      If(NPair(J,1).Eq.ISpec1.And.NPair(J,2).Eq.ISpec2) Then
      IOrbA(I)=1
      NGemA=NGemA+1
      EndIf
      EndDo
C
      If(IOrbA(I).Eq.1) Write(6,*) "ORBITAL NO",I," CENTERED ON ATOMS ",
     $ ISpec1,ISpec2,"IS ACTIVE"
      If(IOrbA(I).Eq.0) Write(6,*) "ORBITAL NO",I," CENTERED ON ATOMS ",
     $ ISpec1,ISpec2,"IS INACTIVE"
C
      EndDo
C
C     ASSIGN INACTIVE ORBITALS TO THE FIRST NGem-NGemA GEMINALS
C
      NGem=NELE
      NGemB=NGem-NGemA
      NGOcc=NGemB
      NActive=NBasis-NGemB
      Write(6,*) 'The number of active orbitals is', NActive
C
      ICount=0
      Do I=1,NBasis
      If(IOrbA(I).Eq.0) Then
      ICount=ICount+1
      IGem(I)=ICount
      EndIf
      EndDo
C
C     ASSIGN ACTIVE ORBITALS TO THE LAST NGemB+1 .... NGem GEMINALS
C
      ICount=NGemB
      Do I=1,NBasis
      If(IOrbA(I).Eq.1) Then
      ICount=ICount+1
      IGem(I)=ICount
      EndIf
      If(ICount.Eq.NGem) ICount=NGemB
      EndDo
C
      Return
      End

*Deck MoldenPrep
      Subroutine MoldenPrep(Title,Occ,UAONO,NGem,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FMultTab
      Include 'commons.inc'
C 
      Character*60 FName,FName2
      Character*100 Line,Aux1
      Character*10 Str
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
      Dimension Occ(NBasis),UAONO(NBasis,NBasis)
C
C     OPEN title.molden file
C
      Do I=1,60
      FName(I:I)=' '
      FName2(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      FName2(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+7)='.molden'
      FName2(K:K+11)='_all.molden'
C      
      Open(10,File=FName)
      Open(20,File=Fname2)
C
    2 Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:4).Eq."[MO]") Then
      GoTo 20
      EndIf
      GoTo 2
C
   20 Continue
C      
      Do I=1,NBasis
C
      Do J=1,3
      Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      EndDo
C
      Read(10,'(A100)')Aux1
      Write(20,'(" Occup=",F12.6)') Two*Occ(I)
C
      Do J=1,NBasis
      Read(10,'(A100)')Aux1
      Write(20,'(I3,F18.14)')J,UAONO(I,J)
      EndDo
C
      EndDo 
C
      Close(10)
      Close(20)
C
C     NOW WRITE GEMINALS TO SEPARATE MOLDEN FILES
C
      Do IG=1,NGem
C
      Do I=1,60
      FName2(I:I)=' '
      EndDo  

      Write(Str,'(I10)')IG
      If(IG.Lt.10) FName2=Trim(Title)//'_g'//Str(10:10)//'.molden'
      If(IG.Ge.10) FName2=Trim(Title)//'_g'//Str(9:10)//'.molden'
C
      Open(10,File=FName)
      Open(20,File=Fname2)
C
    3 Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:4).Eq."[MO]") Then
      GoTo 30
      EndIf
      GoTo 3
C
   30 Continue
C
      Do I=1,NBasis
C
      Do J=1,3
      Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      EndDo
C
      Read(10,'(A100)')Aux1
      If(IGem(I).Eq.IG) Then
      Write(20,'(" Occup=",F12.6)') Two*Occ(I)
      Else
      Write(20,'(" Occup=",F12.6)') Zero
      EndIf
C
      Do J=1,NBasis
      Read(10,'(A100)')Aux1
      Write(20,'(I3,F18.14)')J,UAONO(I,J)
      EndDo
C
      EndDo
C
      Close(20)
      Close(10)
C
      EndDo
C
      Return
      End

*Deck PrintAPSG
      Subroutine PrintAPSG(Title,UMOAO,UNOAO,URe,Occ,
     $           XOne,TwoEl,NGem,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName,FMultTab,Aux1
      Include 'commons.inc'
C 
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
      Dimension UMOAO(NBasis,NBasis),UNOAO(NBasis,NBasis),
     $          URe(NBasis,NBasis),Occ(NBasis),
     $           XOne(NInte1),TwoEl(NInte2)
C
C     LOCAL
C
      Dimension Work(NBasis,NBasis),Work1(NBasis,NBasis),
     $ Gamma(NInte1),Work2(NBasis,NBasis)
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+13)='_APSG_all.dat'
C      
      Open(10,File=FName)
C
C     PRINT THE OVERLAP MATRIX S
C
      Open(20,File='s.dat')

      Read(20,'(A10)')Aux1
      Read(20,'(A10)')Aux1
C
      Do I=1,NBasis
      Do J=1,NBasis/5
      Read(20,222) (Work(I,(J-1)*5+K),K=1,5)
  222 Format(5(F15.8,1X))
      EndDo
      MxK=NBasis-(NBasis/5)*5
      IStart=(NBasis/5)*5
      If(MxK.Ne.0) Read(20,222) (Work(I,IStart+K),K=1,MxK)
      EndDo

      Close(20)
C
      Write(10,'(/,X,"*AO OVERLAP MATRIX S")')
      Do I=1,NBasis
      Write(10,*) (Work(I,J),J=1,NBasis)
      EndDo 
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO MO
C
      Write(10,'(/,X,"*MOLECULAR ORBITALS IN AO BASIS SET")')
      Write(10,'(X,"*MOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (UMOAO(I,J),J=1,NBasis)
      EndDo
C
C     TEST0     UMOAO*S*UMOAO^T = 1
c      Do I=1,NBasis
c      Do J=1,NBasis
c      Sum=0.0D0
c      Do K=1,NBasis
c      Do L=1,NBasis
c      Sum=Sum+UMOAO(I,K)*Work(K,L)*UMOAO(J,L)
c      EndDo
c      EndDo
c      Write(*,*)I,J,Sum
c      EndDo
c      EndDo
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO NO
C
      Write(10,'(/,X,"*NATURAL ORBITALS IN AO BASIS SET")')
      Write(10,'(X,"*NOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (UNOAO(I,J),J=1,NBasis)
      EndDo
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO NO
C     
      Write(10,'(/,X,"*NATURAL ORBITALS IN MO BASIS SET")')
      Write(10,'(X,"*NOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (URe(I,J),J=1,NBasis)
      EndDo
C 
C     TEST1 UNOAO=URe*UMOAO
C
c      Call MultpM(Work,URe,UMOAO,NBasis)
c      Call DiffM(Work,UNOAO,Work1,NBasis)
c      write(*,*)work1
C
C    PRINT Occ, CICOef's, AND ASSIGMENT TO GEMINALS 
C
      Write(10,'(/,X,"*NATURAL OCCUPATION NUMBERS")')
      Write(10,'(X,
     $ "*IN THE RANGE [0,2]. NORMALIZED TO THE NO OF ELECTRONS")')
      Write(10,*)(2.0D0*Occ(I),I=1,NBasis)
C
      Write(10,'(/,X,"*EXPANSION COEFFICIENTS FOR GEMINALS")')
      Write(10,
     $ '(X,"*(c_i)^2=n_i/2. NORMALIZED TO 1 FOR EACH GEMINAL")')
      Write(10,*)(CICoef(I),I=1,NBasis)
C
      Write(10,'(/,X,"*ASSIGNMENT OF ORBITALS TO GEMINALS")')
      Write(10,
     $ '(X,"*INDICES OF GEMINALS TO WHICH ORBITALS BELONG TO")')
      Write(10,*)(IGem(I),I=1,NBasis)
C
C     PRINT 1RDM IN MO
C
      Call VecTr(Gamma,2.0D0*Occ,URe,NBasis)
      Call CpySym(Work2,Gamma,NBasis)
      Write(10,'(/,X,"*1RDM IN MO")')
      Write(10,'(X,"*1RDM NORMALIZED TO THE NUMBER OF ELECTRONS")')
      Do I=1,NBasis
      Write(10,*) (Work2(I,J),J=1,NBasis)
      EndDo
C
C     PRINT 1RDM IN AO
C
      Do I=1,NBasis  
      Do J=1,I
      Work(I,J)=UMOAO(J,I)
      Work(J,I)=UMOAO(I,J)
      EndDo
      EndDo
      Call MatTr(Gamma,Work,NBasis)
      Call CpySym(Work,Gamma,NBasis)
      Write(10,'(/,X,"*1RDM IN AO")')
      Do I=1,NBasis
      Write(10,*) (Work(I,J),J=1,NBasis)
      EndDo
C
C     TEST2 1RDM_AO=UNOAO^T*2*Occ*UNOAO
C
c      Call VecTr(Gamma,2*Occ,UNOAO,NBasis) 
c      IJ=0
c      Do I=1,NBasis
c      Do J=1,I
c      IJ=IJ+1
c      Write(*,*)I,J,Gamma(IJ)-Work(I,J)
c      EndDo
c      EndDo
C
C     COMPUTE AND PRINT 2-RDM IN NO
C
C     RDM2 IS IN THE NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
      Do I=1,60
      FName(I:I)=' '
      EndDo
      K=0
    6 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 6
      EndIf
      FName(K:K+14)='_APSG_2RDM.dat'
      Open(30,File=FName)
      Write(30,'(X,"2RDM IN THE NO REPRESENTATION")')
      Write(30,'(X,"P Q R S 2RDM_pqrs")')  
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,NBasis
      IPQ=IPQ+1
      IRS=0
      Do IR=1,NBasis
      Do IS=1,NBasis
      IRS=IRS+1
C
      RDM2=0.D0
      If(IP.Eq.IQ.And.IR.Eq.IS.And.IGem(IP).Eq.IGem(IR))
     $ RDM2=RDM2+CICoef(IP)*CICoef(IR)
      If(IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2+2.0D0*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      Write(30,'(4I3,E14.6)')IP,IQ,IR,IS,RDM2
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Close(10)
      Close(30)
C
      Return
      End
