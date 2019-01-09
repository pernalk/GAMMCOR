*Deck ACCAS
      Subroutine ACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NGem,System)
C
      use types
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
C     DENSITY MATRIX ELEMENTS
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ IndX(NBasis*(NBasis-1)/2),IndN(2,NBasis*(NBasis-1)/2),
     $ IndAux(NBasis),IPair(NBasis,NBasis)
C
      type(SystemBlock) :: System
C
C
C     CONSTRUCT LOOK-UP TABLES
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      ICount=0
C
      Do I=1,NBasis
C
      If(Occ(I).Lt.One.And.Occ(I).Ne.Zero) Then
      IndAux(I)=1
      Write(6,'(X," Active Orbital: ",I4,E14.4)') I, Occ(I)
      ICount=ICount+1
      EndIf
      EndDo
C
      Write(6,'(/,X," In ACCAS: Active Orbitals ",I4,/)')ICount
C
      IPair(1:NBasis,1:NBasis)=0
C
      Write(LOUT,'(2x,a,e15.5)') 'Threshold for active orbitals: ',
     $ System%ThrSelAct
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
C
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
      If((IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.System%ThrSelAct)
     $ ) Then
C
      Write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals  
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
C
      Ind=Ind+1
      IndX(Ind)=Ind
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      IPair(I,J)=1
      IPair(J,I)=1
C
      EndIf
C
      EndIf
C
c     If(IndAux(I)+IndAux(J).Ne.0 ...
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
      Write(6,'(2X,"Number of pairs reduced to:",I6)')Ind
      Write(6,'(2X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,3I5,2E14.4)')I,Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
      If(IFunSR.Eq.0) Then 
C
      Call RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      ElseIf(IFunSR.Eq.3) Then
C
      Call RunDFOnTop(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      Else 
C
      Call RunACCASLR(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      EndIf
C
      Return
      End

* Deck RunACCAS 
      Subroutine RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      use abmat
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem), EGOne(NGem)
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C     IFlSnd  = 1 - run AC0-CAS (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0-CAS
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1)
     $ Write(6,'(/," *** ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1) 
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1) 
     $ included in ERPA correlation ***",/)')
      EndIf
C
C     CALL AC If IFlAC=1 OR IFlSnd=1
C
      If(IFlAC.Eq.1.Or.IFlSnd.Eq.1) Then
      NGOcc=0
      Call ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      Return
      EndIf
C
C     CALCULATE THE A+B AND A-B MATRICES
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      ACAlpha=One
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
      If(ITwoEl.Eq.3) Then
      Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,
     $ NInte1,'FFOO','FOFO',ACAlpha,.false.)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
      EndIf
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,X,"***************************** ")')
      Write(6,'(  X,"*** ERPA-CAS CALCULATIONS *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(NoSt.Eq.1) Then
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Else
      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      EndIf
      Write(6,'(/,
     $ " *** ERPA-CAS Excitation Energies (a.u., eV) *** ")')
C     the purpose of sorting is only to print a few highest (sorted) eigenvectors
      Call SortEig(1,Eig,ABPLUS,EigVecR,NDimX)
      Do I=1,10
      Write(6,'(I4,4X,2E16.6)') I,Eig(I),27.211*Eig(I)
      EndDo
C
      Write(6,'(/," *** Computing ERPA energy *** ",/)')

      If(ITwoEl.Eq.3) Then
      Call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Occ,
     $ IGem,IndN,IndX,NAcCAS+NInAcCAS,
     $ NDimX,NBasis,'FOFO')
C
      ElseIf(ITwoEl.Eq.1) Then
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      EndIf
      ECorr=Ecorr*Half
      Write
     $ (6,'(/,1X,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
C
      Return
      End

* Deck RunACCASLR
      Subroutine RunACCASLR(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      use types
      use sorter
      use abmat
      use timing
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Character*60 FName
      Real*8, Dimension(:), Allocatable :: TwoEl2
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: SRKer
      Real*8, Dimension(:), Allocatable :: Work
      Dimension NSymNO(NBasis),VSR(NInte1),MultpC(15,15),NumOSym(15),
     $ IndInt(NBasis)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis),
     $ VCoul(NInte1),Den(NBasis,NBasis),Work2(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem),EGOne(NGem),
     $ UAux(NBasis,NBasis)
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C     IFlSnd  = 1 - run AC0-CAS (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0-CAS
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1)
     $ Write(6,'(/," *** LR ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1) 
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1) 
     $ included in ERPA correlation ***",/)')
      EndIf
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,1X,"The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,1X,"The number of Grid Points = ",I8)')
     $ NGrid
C
CC     HAP-TEST
C      UAux=transpose(UNOAO)
C      NOccup=NAcCAS+NInAcCAS
C      Call tran4_gen(NBasis,
C     $        NOccup,UAux(1:NBasis,1:(NOccup)),
C     $        NOccup,UAux(1:NBasis,1:(NOccup)),
C     $        NBasis,UAux,
C     $        NBasis,UAux,
C     $        'FFOOERF','AOERFSORT')
C      Call tran4_gen(NBasis,
C     $        NBasis,UAux,
C     $        NOccup,UAux(1:NBasis,1:(NOccup)),
C     $        NBasis,UAux,
C     $        NOccup,UAux(1:NBasis,1:(NOccup)),
C     $        'FOFOERF','AOERFSORT')

      Allocate  (TwoEl2(NInte2))
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NBasis*NGrid))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (Work(NGrid))
C
c      Call EKT(URe,Occ,XOne,TwoNO,NBasis,NInte1,NInte2)
C
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
C     set/load symmetries of NO's
C
C      Open(10,File='indices_int.dat')
C      Read(10,*) MxSym
C      IStart=0
C      Do I=1,MxSym
C      Read(10,*) X
C      NumOSym(I)=X
C      EndDo
C      Close(10)
C     HAP
      Call create_ind('2RDM',NumOSym,IndInt,NSym,NBasis)
C
      NSymNO(1:NBasis)=0
      IStart=0
      Do I=1,MxSym
      Do J=IStart+1,IStart+NumOSym(I)
C
      Do IOrb=1,NBasis
      If(Abs(UNOAO(IOrb,J)).Gt.1.D-1) NSymNO(IOrb)=I
      EndDo
C
      EndDo
      IStart=IStart+NumOSym(I)
      EndDo
C
C     checking
      Do I=1,MxSym
      II=0
      Do IOrb=1,NBasis
      If(NSymNO(IOrb).Eq.I) II=II+1
      EndDo
      If(II.Ne.NumOSym(I)) 
     $ Stop 'Error in RunACCASLR. Symmetry of NO cannot be established!'
      EndDo
C
      If(MxSym.Eq.1) Then
      MultpC(1,1)=1
      Else
      Do I=1,MxSym
      Do J=1,I
      MultpC(I,J)=IEOR(I-1,J-1)+1
      MultpC(J,I)=MultpC(I,J)
      EndDo
      EndDo
      EndIf
C
C     if IFunSR.Gt.0 (but different from 3) TwoEl includes lr-integrals with erf/r 
C     load full-range two-electron integrals and transform to NO
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+10)='.reg.integ'
C      Call Int2_AO(TwoEl2,NumOSym,MultpC,FName,NInte1,NInte2,NBasis)
      Call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT')
      Call LoadSaptTwoEl(3,TwoEl2,NBasis,NInte2)
C
      Write(6,'(" Transforming two-electron erf integrals ...")')
      Call TwoNO1(TwoEl2,UNOAO,NBasis,NInte2)
C
      If(ITwoEl.Eq.3) Then
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
C     TRANSFORM J AND K
      UAux=transpose(UNOAO)
      Call tran4_gen(NBasis,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        NBasis,UAux,
     $        'FFOO','AOTWOSORT')
      Call tran4_gen(NBasis,
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        'FOFO','AOTWOSORT')
C
      EndIf
C     
C     compute sr potential (vsr=xc+hartree)
C     as a byproduct a sr energy (ensr=sr-xc+sr-hartree) is obtained
C
      IFunSave=IFunSR
      If(IFunSR.Eq.4) IFunSR=2
C
      Den=0
      Work2=0
      Do I=1,NBasis
      Call dger(NBasis,NBasis,1d0*Occ(I),
     $          UNOAO(I,:),1,UNOAO(I,:),1,Den,NBasis)
      EndDo
C 
      IJ = 0
      Do J=1,NBasis      
      Do I=1,J
      IJ = IJ + 1
      Work2(IJ) = Den(I,J)
      EndDo
      EndDo
C
      Write(6,'(/,1X,"RANGE PARAMETER ",F8.3)')Alpha
C
      Call PotCoul_mithap(VCoul,Work2,.true.,'AOERFSORT',NBasis)
      Print*, 'VCoul',norm2(VCoul)
      Call EPotSR(EnSR,EnHSR,VSR,Occ,URe,UNOAO,.false.,
     $        OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,
C     $        NSymNO,VCoul,Alpha,IFunSR,
     $        NSymNO,VCoul,TwoEl2,TwoNO,Alpha,IFunSR,
     $        NGrid,NInte1,NInte2,NBasis)
C
C      Call EPotSR(EnSR,EnHSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,
C     $ OrbZGrid,WGrid,NSymNO,TwoEl2,TwoNO,NGrid,NInte1,NInte2,NBasis)
      Write(6,'(" SR xcH Energy",F15.8)')EnSR
C
C     CALCULATE THE SR_XC_PBE ENERGY WITH "TRANSLATED" ALPHA AND BETA DENSITIES
C     [as in Gagliardi J. Chem. Phys. 146, 034101 (2017)]
C
      Call SR_PBE_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(/," SR_xc_PBE with translated densities",F15.8,/)')
     $ EXCTOP
C
      IFunSR=IFunSave
C
      XVSR=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      XVSR=XVSR+Two*Occ(I)*VSR(II)
      EndDo
C
      If(IFunSR.Eq.4) Then
      Print*, 'VSR',norm2(VSR)
      Write(6,'(X,/,
     $"*** ADDING VSR_HXC TO A ONE-ELECTRON HAMILTONIAN*** ",/)')
      Do I=1,NInte1
      XOne(I)=XOne(I)+VSR(I)
      EndDo
      EndIf
C
      Allocate (SRKer(NGrid))
      SRKer(1:NGrid)=Zero
C
C      IFunSRKer = 1
C
      If(IFunSRKer.Eq.1) Then
C
      Write(6,'(/," *** Generating a sr-kernel on a grid ***")')
C
      Call GetSRKer(SRKer,Occ,URe,OrbGrid,WGrid,NBasis,NGrid)
      Do I=1,NGrid
      Work(I)=WGrid(I)*SRKer(I)
      EndDo 
C
      EndIf
C
      If(IFlSnd.Eq.1) Then
C
C     ***** LR-AC0 CALCULATION **********************
C
      Write(6,'(  X,"*** LR-AC0-CAS CALCULATION *** ")')
C
      If(ITwoEl.Eq.3) Then
      Call Y01CASLR_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,
     $ MultpC,NSymNO,
     $ SRKer,WGrid,OrbGrid,
     $ 'PROP0','PROP1',
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,
     $ NGrid,NDimX,NBasis,NDimX,NInte1,
     $ 'EMPTY','FFOOERF','FOFOERF',0,ECASSCF,ECorr)
CC
C      Call Y01CAS_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,
C     $ 'PROP0','PROP1',
C     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,
C     $ NBasis,NDimX,NInte1,'EMPTY','FFOOERF','FOFOERF',0,
C     $ ETot,ECorr)
CC
C      Call AC0CAS(ECorr,ETot,TwoNO,Occ,URe,XOne,
C     $ ABPLUS,ABMIN,EigVecR,Eig,
C     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
CC
      ElseIf(ITWoEl.Eq.1) Then
      Call AC0CASLR(ECorr,ECASSCF,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
     $ TwoEl2,OrbGrid,Work,NSymNO,MultpC,NGrid)
C
      EndIf
C
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC0-Corr, Total'',6X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC0-Corr,Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr
C
      GoTo 777
C
      EndIf 
C
C     ****** LR-AC CALCULATION *******************************************************
C
      If(IFlAC.Eq.1.And.IFlSnd.Eq.0) Then
C
      Write(6,'(  X,"*** LR-AC-CAS CALCULATION *** ")')
C
      If(IFunSRKer.Eq.1) Then 
C
C     GENERATE A SR KERNEL AND DUMP IT
C
      If(ITwoEl.Eq.3) Then
      Call ModABMin_FOFO(Occ,SRKer,WGrid,OrbGrid,ABMIN,
     $                   MultpC,NSymNO,
     $                   IndN,IndX,NDimX,NGrid,NBasis,
     $                   NAcCAS+NInAcCAS,'FOFO','FOFOERF',
     $                   'srdump')
C 
      ElseIf(ITwoEl.Eq.1) Then
      Open(20,File="srdump",Form='UNFORMATTED')
C
      Do I=1,NBasis
      CICoef(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CICoef(I)=-CICoef(I)
      EndDo
C
      Do IRow=1,NDimX
C
      IA=IndN(1,IRow)
      IB=IndN(2,IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NDimX
C
      IC=IndN(1,ICol)
      ID=IndN(2,ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      XKer1234=Zero
C
      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
      ISym=MultpC(I1I2S,I3I4S)
      If(ISym.Eq.1) Then
      Do I=1,NGrid
      XKer1234=XKer1234+Work(I)*
C     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
C     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
     $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
     $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)
      EndDo
      EndIf
C
      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
C
      Write(20) Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      EndDo
      EndDo      
C
      Close(20)
C
C     if TwoEl
      EndIf
C
      EndIf
C
      NGOcc=0
      Call ACECORR(ECASSCF,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX) 
      ECorr=EGOne(1)
C 
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC-Corr, Total'',6X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC-Corr, Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr

C
      GoTo 777
C
      EndIf
C
C     ****** LR-AC1 (aka ERPA1) CALCULATION ******************************************
C
C     compute AB matrices with the lr integrals and a modified potential
C     TwoNO are lr-integrals!  xone already includes vsr (if properly saved in molpro)
C
      ACAlpha=One
C
      If(ITwoEl.Eq.3) Then
      Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,
     $ NInte1,'FFOOERF','FOFOERF',ACAlpha,.false.)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
      EndIf
C
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
C
      Do I=1,NBasis
      CICoef(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CICoef(I)=-CICoef(I)
      EndDo
C
      If(IFunSRKer.Eq.1) Then
C
C     ADD A SR-ALDA KERNEL
C
      Call clock('START',Tcpu,Twall)
C
C      WITHOUT BATCHES:
C
C      Do IRow=1,NDimX
CC
C      IA=IndN(1,IRow)
C      IB=IndN(2,IRow)
C      CA=CICoef(IA)
C      CB=CICoef(IB)
CC
C      Do ICol=1,NDimX
CC
C      IC=IndN(1,ICol)
C      ID=IndN(2,ICol)
C      CC=CICoef(IC)
C      CD=CICoef(ID)
CC
C      XKer1234=Zero
CC
C      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
C      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
C      ISym=MultpC(I1I2S,I3I4S)
CC
C      If(ISym.Eq.1) Then
C      Do I=1,NGrid
C      XKer1234=XKer1234+Work(I)* 
CC     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
CC     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
C     $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
C     $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)
C      EndDo
C      EndIf
CC
C      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
CC
C      ABMIN((ICol-1)*NDimX+IRow)=ABMIN((ICol-1)*NDimX+IRow)
C     $ +Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
CC
C      EndDo
C      EndDo
C
C
      If(ITwoEl.Eq.3) Then
      Call ModABMin_FOFO(Occ,SRKer,WGrid,OrbGrid,ABMIN,
     $                   MultpC,NSymNO,
     $                   IndN,IndX,NDimX,NGrid,NBasis,
     $                   NAcCAS+NInAcCAS,'FOFO','FOFOERF')
C      Print*,'ABMIN-MY',norm2(ABMIN)
C     
      ElseIf(ITwoEl.Eq.1) Then
      Call ModABMinSym(Occ,SRKer,WGrid,OrbGrid,TwoEl2,TwoNO,ABMIN,
     $          MultpC,NSymNO,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
C      Print*,'ABMIN-KA',norm2(ABMIN)
C
      EndIf
C
      Write(6,'(1X," *** sr-kernel added ***")')
      Call clock('sr-kernel',Tcpu,Twall)
C
      EndIf
C
      Write(6,'(/,1X,"*** LR-ERPA-CAS CALCULATION *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(NoSt.Eq.1) Then
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Else
      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      EndIf
      Write(6,'(/,
     $ " *** LR-CAS-SR-DFT Excitation Energies (a.u., eV) *** ")')
C     the purpose of sorting is only to print a few highest (sorted) eigenvectors
      Call SortEig(1,Eig,ABPLUS,EigVecR,NDimX)
C
      Do IS=1,Min(5,MxSym)
      NExcit=0
C
      Do I=1,NDimX
c      If(NExcit.Lt.5) Then
      If(NExcit.Lt.20) Then
C
C     find symmetry
      XMax=Zero
      JMax=0
      Do J=1,NDimX
      If(Abs(EigVecR((I-1)*NDimX+J)).Gt.XMax) Then
      XMax=Abs(EigVecR((I-1)*NDimX+J))
      JMax=J
      EndIf
      EndDo
      IP=IndN(1,JMax)
      IQ=IndN(2,JMax)
      ISym=MultpC(NSymNO(IP),NSymNO(IQ))
C
      If(ISym.Eq.IS) Then
      NExcit=NExcit+1
      Write(6,'("Excit",I2,I2,4X,2E16.6)')NExcit,ISym,
     $ Eig(I),27.211*Eig(I)
      EndIf
C
      EndIf
      EndDo
C
      EndDo
C
      Write(6,'(/," *** Computing LR-ERPA energy *** ",/)')

      If(ITwoEl.Eq.3) Then
      Call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Occ,
     $ IGem,IndN,IndX,NAcCAS+NInAcCAS,
     $ NDimX,NBasis,'FOFOERF')
C
      ElseIf(ITwoEl.Eq.1) Then
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      EndIf
C
      ECorr=Ecorr*Half
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC1-Corr, Total'',6X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC1-Corr,Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr
C
  777 Continue
C
      DeAllocate  (TwoEl2)
      DeAllocate  (WGrid)
      DeAllocate  (OrbGrid)
      DeAllocate  (OrbXGrid)
      DeAllocate  (OrbYGrid)
      DeAllocate  (OrbZGrid)
      DeAllocate  (Work)
      DeAllocate (SRKer)
C
      Return
      End

*Deck GetSRKer
      Subroutine GetSRKer(SRKer,Occ,URe,OrbGrid,WGrid,NBasis,NGrid)
C
C     RETURNS a SR-ALDA KERNEL ON THE GRID 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension SRKer(NGrid),Occ(NBasis),URe(NBasis,NBasis),
C     $ OrbGrid(NBasis,NGrid),WGrid(NGrid)
     $ OrbGrid(NGrid,NBasis),WGrid(NGrid)
C
C     LOCAL ARRAYS
C
      Dimension RhoVec(NGrid)
C
      Do I=1,NGrid
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      RhoVec(I)=Rho
      EndDo
C
      Call RhoKernel(RhoVec,SRKer,2,Alpha,NGrid)
C
      Return
      End

*Deck AC0CASLR
      Subroutine AC0CASLR(ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigY,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
     $ TwoEl2,OrbGrid,SRKerW,NSymNO,MultpC,NGrid)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C   
      use timing
      use abmat
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Integer,Parameter :: Maxlen = 128
C
      Dimension
     $ URe(NBasis,NBasis),XOne(NInte1),Occ(NBasis),TwoNO(NInte2),
     $ IndAux(NBasis),
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ Eig(NDimX),EigY(NDimX*NDimX),IndX(NDim),IndN(2,NDim),
     $ TwoEl2(NInte2),OrbGrid(NBasis*NGrid),SRKerW(NGrid),
     $ NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Integer :: Offset,Batchlen
      Real*8, Allocatable :: RDM2Act(:)
      Double Precision, Allocatable :: Work(:),Batch(:,:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ EigX(NDimX*NDimX),
     $ IEigAddY(2,NDimX),IEigAddInd(2,NDimX),IndBlock(2,NDimX),
     $ XMAux(NDimX*NDimX)
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
C     AUXILIARY STUFF LATER NEEDED TO GET A+ AND A- MATRICES FOR ALPHA=0
C
C     ONE-ELECTRON MATRIX IN A NO REPRESENTATION
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*HNO(II)
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETot=ETot+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)')ETot
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA=0 HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=Zero
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRIX AuxI AND AuxIO
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      AuxIO(IPQ)=Zero
      Do IT=1,NOccup
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.1) Then
       AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndIf
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRIX WMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
      WMAT(I,J)=Zero
      EndDo
      EndDo
C
      Do IP=1,NBasis
      Do IR=1,NOccup
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.1)
     $ WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))
     $ *FRDM2(IW,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IT,IU,IP,IW))
     $ *FRDM2(IW,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-ACTIVE BLOCK
C
      Write(6,'(" *** ACTIVE-ACTIVE BLOCK ***")')
C
      NFree1=1
      NFree2=1
      NoEig=0
C
      NDimB=0
      Do IQQ=1,NAct
      Do IPP=IQQ+1,NAct
      IP=Ind1(IPP)
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IQQ=1,NAct
      Do IPP=IQQ+1,NAct
      IP=Ind1(IPP)
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      Do ISS=1,NAct
      Do IRR=ISS+1,NAct
      IR=Ind1(IRR)
      IS=Ind1(ISS)
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C    
C     ADD A SR KERNEL
C
      If(IFunSRKer.Eq.1) Then
      I1I2S=MultpC(NSymNO(IP),NSymNO(IQ))
      I3I4S=MultpC(NSymNO(IR),NSymNO(IS))
      ISym=MultpC(I1I2S,I3I4S)
      XKer1234=Zero
      If(ISym.Eq.1) Then
      Do I=1,NGrid
      XKer1234=XKer1234+SRKerW(I)*
C     $ OrbGrid(IP+(I-1)*NBasis)*OrbGrid(IQ+(I-1)*NBasis)*
C     $ OrbGrid(IR+(I-1)*NBasis)*OrbGrid(IS+(I-1)*NBasis)
     $ OrbGrid(I+(IP-1)*NGrid)*OrbGrid(I+(IQ-1)*NGrid)*
     $ OrbGrid(I+(IR-1)*NGrid)*OrbGrid(I+(IS-1)*NGrid)
      EndDo
      EndIf  
      TwoSR=TwoEl2(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IQ,IR,IS))
      ABM=ABM+Four*(C(IP)+C(IQ))*(C(IR)+C(IS))*(XKer1234+TwoSR)
      EndIf
C
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABMIN((IRow-1)*NDimB+ICol)=ABM

C
      EndIf
C
      EndIf
      EndDo
      EndDo
C
      EndIf
      EndDo
      EndDo
C
      If(NDimB.Ne.0) Then
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS
C
      Write(6,'(" *** ACTIVE-INACTIVE BLOCKS ***")')
C
      Do IQ=1,INActive
C
      NDimB=0
      Do IPP=1,NAct
      IP=Ind1(IPP)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IPP=1,NAct
      IP=Ind1(IPP)
C
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      IS=IQ
      Do IRR=1,NAct
      IR=Ind1(IRR)
C
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABMIN((IRow-1)*NDimB+ICol)=ABM
C
      EndIf
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      If(NDimB.Ne.0) Then
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     Do IP
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-ACTIVE BLOCKS
C
      Write(6,'(" *** VIRTUAL-ACTIVE BLOCKS ***")')
C
      Do IP=NOccup+1,NBasis
C
      NDimB=0
      Do IQQ=1,NAct
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IQQ=1,NAct
      IQ=Ind1(IQQ)
C
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      IR=IP
      Do ISS=1,NAct
      IS=Ind1(ISS)
C
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABMIN((IRow-1)*NDimB+ICol)=ABM
C
      EndIf
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      If(NDimB.Ne.0) Then
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $    NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     Do IP
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-INACTIVE BLOCKS
C
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
C
      NDimB=0
C
      If(IPair(IP,IQ).Eq.1) Then
C
      NDimB=1
      IndBlock(1,NFree1)=IP
      IndBlock(2,NFree1)=IQ
C
      IEigAddY(1,NFree1)=NFree2
      IEigAddY(2,NFree1)=IEigAddY(1,NFree1)
      IEigAddInd(1,NFree1)=NFree1
      IEigAddInd(2,NFree1)=NFree1
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IP,IQ,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      Eig(NFree1)=ABP
      EigY(NFree2)=One/Sqrt(Two)
      EigX(NFree2)=One/Sqrt(Two)
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
      EndIf
C
      EndDo
      EndDo
C
      Write(6,'(/," *** DONE WITH 0TH-ORDER IN AC0-CASSCF ***")')
C
      Print*, 'Eig,Y,X',norm2(Eig),norm2(EigY),norm2(EigX)
C     DONE 0TH-ORDER CALCULATIONS
C
      Write(6,'(/,
     $" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
      Print*, 'AB1-KA',norm2(ABPLUS),norm2(ABMIN)
C
C      Call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,
C     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,
C     $ NInte1,'FFOOERF','FOFOERF',1d0,.true.)
CC
C      Print*, 'AB1-MY',norm2(ABPLUS),norm2(ABMIN)
C     ADD A SR KERNEL
C
      If(IFunSRKer.Eq.1) Then
C
      Write(6,'(/," *** ADDING THE SR KERNEL ***" )')
C
      Do IRow=1,NoEig
C
      IA=IndBlock(1,IRow)
      IB=IndBlock(2,IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NoEig
C
      IC=IndBlock(1,ICol)
      ID=IndBlock(2,ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      XKer1234=Zero
C
      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
      ISym=MultpC(I1I2S,I3I4S)
C
      If((ISym.Eq.1).And.(IGFact(NAddr3(IA,IB,IC,ID)).Eq.0)) Then
      Do I=1,NGrid
      XKer1234=XKer1234+SRKerW(I)*
C     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
C     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
     $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
     $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)

      EndDo
      EndIf
C
      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
C
      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
     $ +Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      EndDo
      EndDo
C
CC     ADD KERNEL IN BATCHES:
C      call clock('START',Tcpu,Twall)
C      Allocate(Work(Maxlen),Batch(Maxlen,NBasis))
CC
C      Do Offset=0,NGrid,Maxlen
C      Batchlen = min(NGrid-Offset,Maxlen)
C      If(Batchlen==0) exit
CC    
C      Work(1:Batchlen) = SRKerW(Offset+1:Offset+Batchlen)
C      Call FILL_BATCH(OrbGrid,Batch,Batchlen,Offset,NGrid,NBasis)
CC
C      Do IRow=1,NoEig
CC
C      IA=IndBlock(1,IRow)
C      IB=IndBlock(2,IRow)
C      CA=CICoef(IA)
C      CB=CICoef(IB)
CC
C      Do ICol=1,NoEig
CC
C      IC=IndBlock(1,ICol)
C      ID=IndBlock(2,ICol)
C      CC=CICoef(IC)
C      CD=CICoef(ID)
CC
C      XKer1234=Zero
CC
C      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
C      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
C      ISym=MultpC(I1I2S,I3I4S)
CC
C      If((ISym.Eq.1).And.(IGFact(NAddr3(IA,IB,IC,ID)).Eq.0)) Then
C      Do I=1,Batchlen
C      XKer1234=XKer1234+Work(I)*
C     $ Batch(I,IA)*Batch(I,IB)*Batch(I,IC)*Batch(I,ID)
C      endDo
C      EndIf
CC
C      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
C     $ +Four*(CA+CB)*(CD+CC)*XKer1234
CC
C      EndDo
C      EndDo
C      EndDo
CC
CC     ADD INTEGRALS      
CC
C      Do IRow=1,NoEig
CC
C      IA=IndBlock(1,IRow)
C      IB=IndBlock(2,IRow)
C      CA=CICoef(IA)
C      CB=CICoef(IB)
CC
C      Do ICol=1,NoEig
CC
C      IC=IndBlock(1,ICol)
C      ID=IndBlock(2,ICol)
C      CC=CICoef(IC)
C      CD=CICoef(ID)
CC
C      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
CC
C      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
C     $ +Four*(CA+CB)*(CD+CC)*TwoSR
CC
C      EndDo
C      EndDo
CC 
C      Deallocate(Batch,Work)
CC
C      Write(6,'("*** sr-kernel added ***")')
C      call clock('sr-kernel AB(1)',Tcpu,Twall)
C
      EndIf
C
      Deallocate(RDM2Act)
C
      Print*, 'ABM-KA',norm2(ABMIN)
C      
      Write(6,'(/," *** DONE WITH COMPUTING AB(1) MATRICES ***")')
C
C     1ST-ORDER PART
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      XMAux(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      XMAux(NU+(MU-1)*NoEig)=XMAux(NU+(MU-1)*NoEig)
     $ +ABPLUS(NU+(I-1)*NoEig)*EigX(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      ABPLUS(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
      ABPLUS(NU+(MU-1)*NoEig)=ABPLUS(NU+(MU-1)*NoEig)
     $ +EigX(IStart+II)*XMAux(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      XMAux(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      XMAux(NU+(MU-1)*NoEig)=XMAux(NU+(MU-1)*NoEig)
     $ +ABMIN(NU+(I-1)*NoEig)*EigY(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      ABMIN(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
      ABMIN(NU+(MU-1)*NoEig)=ABMIN(NU+(MU-1)*NoEig)
     $ +EigY(IStart+II)*XMAux(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      XMAux(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C KP 03/04/2018
C
      If(Eig(MU).Ne.Zero) Then
C
      Do NU=1,NoEig
C
C KP 03/04/2018
C
      If(Eig(NU).Ne.Zero) Then
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
C
      XMAux(MU+(I-1)*NoEig)=XMAux(MU+(I-1)*NoEig)+Two*
     $ (ABPLUS(MU+(NU-1)*NoEig)-ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)+Eig(NU))*EigY(IStart+II)
C
      II=II+1
      EndDo
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      ABPLUS(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
C
      Do J=1,NoEig
      ABPLUS(I+(J-1)*NoEig)=ABPLUS(I+(J-1)*NoEig)+XMAux(MU+(J-1)*NoEig)
     $ *EigY(IStart+II)
      EndDo
C
      II=II+1
      EndDo
      EndDo
C
C     FINALLY THE ENERGY CORRECTION
C
      EAll=Zero
      EIntra=Zero
C
      Do I=1,NoEig
C
      IP=IndBlock(1,I)
      IR=IndBlock(2,I)
C
      Do J=1,NoEig
C
      IQ=IndBlock(1,J)
      IS=IndBlock(2,J)
C
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      SumY=ABPLUS(I+(J-1)*NoEig)
      Aux=(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      EAll=EAll+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))
     $ EIntra=EIntra+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
C     endinf of If(IP.Gt.IR.And.IQ.Gt.IS)
      EndIf
C
      EndDo
      EndDo
C
      Print*, 'EAll,EIntra',EAll,EIntra
      ECorr=EAll-EIntra
C
      Return
      End

*Deck fil_Batch
      Subroutine FILL_BATCH(OrbGrid,Batch,Batchlen,Offset,NGrid,NBasis)
C
C     FILLS BATCHES FOR GRID CALCULATIONS
C
      Implicit None
C
      Integer :: Batchlen,Offset,NGrid,NBasis
      Double Precision :: Batch(Batchlen,NBasis),OrbGrid(NGrid,NBasis)

      Batch(1:Batchlen,1:NBasis) =
     $ OrbGrid(Offset+1:Offset+Batchlen,1:NBasis)

      End Subroutine FILL_BATCH

*Deck GGA_ONTOP
      Subroutine GGA_ONTOP(EXCTOP,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A GGA_XC ENERGY IF THE DENSITY AND SPIN-DENSITY ARE COMPUTED AS:
C     RHO_A/B = 1/2 ( RHO +/- SQRT(RHO^2 - 2 PI )  )
C     WHERE PI IS THE ON-TOP PAIR DENSITY, see Gagliardi JCP 146, 034101 (2017)
C
C     XCFUN IS USED !!!
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
      Dimension Zk(NGrid),RhoA(NGrid),RhoB(NGrid),
     & SigmaAA(NGrid),SigmaAB(NGrid),SigmaBB(NGrid)
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid,Occ,URe,OrbGrid,NGrid,NBasis)
C
      If(RhoGrid.Gt.1.D-12) Then
C
      OnTop=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop=OnTop
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      Rho=RhoGrid
      R=Two*OnTop/Rho**2
      XFactor=Zero
      If(R.Lt.One) XFactor=SQRT(One-R)
C 
      RhoA(I)=Rho/Two*(One+XFactor)
      RhoB(I)=Rho/Two*(One-XFactor)
C
      RhoXa=RhoX/Two*(One+XFactor)
      RhoXb=RhoX/Two*(One-XFactor)
      RhoYa=RhoY/Two*(One+XFactor)
      RhoYb=RhoY/Two*(One-XFactor)
      RhoZa=RhoZ/Two*(One+XFactor)
      RhoZb=RhoZ/Two*(One-XFactor)
C
      SigmaAA(I)=RhoXa*RhoXa+RhoYa*RhoYa+RhoZa*RhoZa
      SigmaAB(I)=RhoXa*RhoXb+RhoYa*RhoYb+RhoZa*RhoZb
      SigmaBB(I)=RhoXb*RhoXb+RhoYb*RhoYb+RhoZb*RhoZb
C
      Else
      RhoA(I)=Zero
      RhoB(I)=Zero
      SigmaAA(I)=Zero
      SigmaAB(I)=Zero
      SigmaBB(I)=Zero
C
      EndIf
C
      EndDo
C
      Call dfun_PBE_AB(RhoA,RhoB,SigmaAA,SigmaAB,SigmaBB,Zk,NGrid)
C
      EXCTOP=Zero
      Do I=1,NGrid
      EXCTOP=EXCTOP+Zk(I)*WGrid(I)
      EndDo  
C
      Return
      End

*Deck SR_PBE_ONTOP
      Subroutine SR_PBE_ONTOP(EXCTOP,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A SR-PBE XC ENERGY IF THE DENSITY AND SPIN-DENSITY ARE COMPUTED AS:
C     RHO_A/B = 1/2 ( RHO +/- SQRT(RHO^2 - 2 PI )  )
C     WHERE PI IS THE ON-TOP PAIR DENSITY, see Gagliardi JCP 146, 034101 (2017)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C     $ WGrid(NGrid),OrbGrid(NBasis,NGrid),OrbXGrid(NBasis,NGrid),
C     $ OrbYGrid(NBasis,NGrid),OrbZGrid(NBasis,NGrid)
C
! input
      Dimension Zk(NGrid),RhoGrid(NGrid),RhoO(NGrid),
     & Sigma(NGrid),SigmaCO(NGrid),SigmaOO(NGrid)
      logical fderiv,open
! output
      integer igrad
      character*(30) name
      double precision vrhoc(ngrid),vrhoo(ngrid)
      double precision vsigmacc(ngrid),vsigmaco(ngrid),vsigmaoo(ngrid)
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
C
      If(RhoGrid(I).Gt.1.D-12) Then
C
      OnTop=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop=OnTop
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C     $ *OrbGrid(IP,I)*OrbGrid(IQ,I)*OrbGrid(IR,I)*OrbGrid(IS,I)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      Rho=RhoGrid(I)
      R=Two*OnTop/Rho**2
      XFactor=Zero
      If(R.Lt.One) XFactor=SQRT(One-R)
C
      Rhoa=Rho/Two*(One+XFactor)
      Rhob=Rho/Two*(One-XFactor)
      RhoGrid(I)=Rhoa+Rhob
      RhoO(I)=Rhoa-Rhob
C
      RhoXa=RhoX/Two*(One+XFactor)
      RhoXb=RhoX/Two*(One-XFactor)
      RhoYa=RhoY/Two*(One+XFactor)
      RhoYb=RhoY/Two*(One-XFactor)
      RhoZa=RhoZ/Two*(One+XFactor)
      RhoZb=RhoZ/Two*(One-XFactor)
C
      RhoXC=RhoXa+RhoXb
      RhoYC=RhoYa+RhoYb
      RhoZC=RhoZa+RhoZb

      RhoXO=RhoXa-RhoXb
      RhoYO=RhoYa-RhoYb
      RhoZO=RhoZa-RhoZb
C
      Sigma(I)=  RhoXC*RhoXC+RhoYC*RhoYC+RhoZC*RhoZC
      SigmaCO(I)=RhoXC*RhoXO+RhoYC*RhoYO+RhoZC*RhoZO
      SigmaOO(I)=RhoXO*RhoXO+RhoYO*RhoYO+RhoZO*RhoZO
C
      Else
      RhoGrid(I)=Zero
      RhoO(I)=Zero
      Sigma(I)=Zero
      SigmaCO(I)=Zero
      SigmaOO(I)=Zero
C
      EndIf
C
      EndDo
C
      EXCTOP=Zero
      FDeriv=.False.
      Open=.True.
C 
      Call dftfun_exerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,RhoO,
     >                   Sigma,SigmaCO,SigmaOO,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)

      Do I=1,NGrid
      EXCTOP=EXCTOP+Zk(I)*WGrid(I)
      EndDo
C
      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,RhoO,
     >                   Sigma,SigmaCO,SigmaOO,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
      EnC=Zero
      Do I=1,NGrid
      EnC=EnC+Zk(I)*WGrid(I)
      EndDo
C
      EXCTOP=EXCTOP+EnC
C
      Return
      End

*Deck RunDFOnTop
      Subroutine RunDFOnTop(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C 
C     ETot is calculated from MC-PDFT
C     with PBE xc functional
C     see Eq.(6) in Manni, et al. JCTC 10, 3669-3680 (2014)
C     doi: 10.1021/ct500483t 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Character*60 FName
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RDM2Act
      Dimension NSymNO(NBasis),VSR(NInte1),MultpC(15,15),NumOSym(15)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis),Ind1(NBasis),Ind2(NBasis)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem), EGOne(NGem)
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
C     COMPUTE THE ENERGY FOR CHECKING
C
      EOne=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      EOne=EOne+Two*Occ(I)*XOne(II)
      EndDo
      Write(6,'(/,1X,''One-Electron Energy'',6X,F15.8)')EOne
C
      ETot=EOne
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETot=ETot+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(1X,''CASSCF Energy (w/o ENuc)'',X,F15.8)')ETot
      Write(6,'(1X,''Total CASSCF Energy '',5X,F15.8)')ETot+ENuc
C      
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/," The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NBasis*NGrid))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
C
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
      Call GGA_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(/," PBE_xc from xcfun with tr densities",F15.8,/)')
     $ EXCTOP
C
c herer!!! 
c      Alpha=1.D-12
c      Call SR_PBE_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
c     $ OrbZGrid,WGrid,NGrid,NBasis)
c      Write(6,'(/," PBE_xc with translated densities",F15.8,/)')
c     $ EXCTOP
C
      EnH=Zero
      Do I=1,NOccup
      Do J=1,NOccup
      EnH=EnH+Two*Occ(I)*Occ(J)*TwoNO(NAddr3(I,I,J,J))
      EndDo
      EndDo
C
      ETot=EOne+EnH+EXCTOP
      Write(6,'(1X,''E_Hartree                     '',X,F15.8)') EnH 
      Write(6,'(1X,''EOne+EH+E_xc[OnTop] (w/o ENuc)'',X,F15.8)') ETot
      Write(6,'(1X,''EOne+EH+E_xc[OnTop] +ENuc'',5X,F15.8)')ETot+ENuc
C
      DeAllocate (RDM2Act)
      DeAllocate  (WGrid)
      DeAllocate  (OrbGrid)
      DeAllocate  (OrbXGrid)
      DeAllocate  (OrbYGrid)
      DeAllocate  (OrbZGrid)
C
      Return
      End
