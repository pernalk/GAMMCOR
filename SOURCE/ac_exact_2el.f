*Deck ACPINO
      Subroutine ACPINO(ENuc,TwoEl,Occ,XOne,
     $ NBasis,NInte1,NInte2,NDim,NGem,NoEig)
C
C     IT WORKS ONLY FOR TWO-ELECTRON SYSTEMS!!! 
C
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
C     DENSITY MATRIX ELEMENTS
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ Occ(NBasis),TwoEl(NInte2),XOne(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ UReAlph(NBasis,NBasis),OccAlph(NBasis),CIAlph(NBasis),
     $ TwoAlph(NInte2),H1Alph(NInte1),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ CMAT(NDim*NDim),EMAT(NBasis*NBasis),DMAT(NDim*NBasis),
     $ EMATM(NBasis*NBasis),DMATK(NDim*NBasis),
     $ IndX(NDim),IndN(2,NDim),
     $ EigVecR(2*(NDim+NBasis)*2*(NDim+NBasis)),
     $ Eig(2*(NDim+NBasis)),IndAux(NBasis),
     $ XGrid(100), WGrid(100),CISave(NBasis),IGemS(NBasis),
     $ TrGamm(NInte1,NInte1),EExcit(NInte1),GammAl(NInte1)
C
      goto 777
C
      Open(10,File="rdm2.dat",Status='Old')
      Read(10,'(4I4,F19.12)')I,J,K,L,X
      If(.NOT.(I.Eq.1.And.J.Eq.1.And.K.Eq.1.And.L.Eq.1)) Then
      Write(*,*) 
     $ "Fatal Error: Element 1,1,1,1 of 2-RDM not available!"
      Write(*,*)'Skipping stability analysis'
      GoTo 777
      EndIf
C
      CICoef(1)=SQRT(X/2.D0)
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
      X=X/2.D0
      If(I.Eq.1.And.K.Eq.1.And.J.Eq.L) CICoef(J)=X/CICoef(1)
      GoTo 10
   40 Continue
      Close(10)
      Write(6,'(/,3X,"I      SQRT(Occ)       CICoef (from rdm2.dat)")')
      SumC=Zero
      SumO=Zero
      Do I=1,NBasis
      Write(6,'(I4,2X,2E16.6)')I,SQRT(Occ(I)),CICoef(I)
      SumC=SumC+CICoef(I)**2
      SumO=SumO+Occ(I)
      EndDo
      Write(6,'(2X,"SumOcc SumC2",2E16.6)')SumO,SumC
C
C ************* TEST THE STABILITY OF THE CAS SOLUTION **************
C
      Write(6,'(//,"********************* 
     $ STABILITY TEST ************************")')
C     CONSTRUCT A LOOK-UP TABLE
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
      ICount=0
      Do I=1,NBasis
      If(Occ(I).Gt.Zero) Then
      IndAux(I)=1
      ICount=ICount+1
      EndIf
      EndDo
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
      If((IGem(I).Ne.IGem(J)).And.(IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-20) ) Then
      Write(*,*)"Discarding nearly degenerate pair",I,J
      Else
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
      EndIf
      EndDo
      EndDo
      NDimX=Ind
      NDimN=NBasis
C
      Do I=1,NBasis
      Do J=1,NBasis
      If(I.Eq.J) UReAlph(I,J)=One
      If(I.Ne.J) UReAlph(I,J)=Zero
      EndDo
      EndDo
C
      Call APSG_NEST(ABPLUS,ABMIN,CMAT,EMAT,EMATM,DMAT,DMATK,
     $ UReAlph,Occ,XOne,TwoEl,
     $ NBasis,NDim,NInte1,NInte2,NGem,2)
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1)
      EndDo
      EndDo
      Do J=1,NDimN
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(J-1)*NDim+IndX(I)
      DMAT(IJ)=DMAT(IJ1)
      DMATK(IJ)=DMATK(IJ1)
      EndDo
      EndDo
      Do J=1,NDimN
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(J-1)*NBasis+I
      EMAT(IJ)=EMAT(IJ1)
      EMATM(IJ)=EMATM(IJ1)
      EndDo
      EndDo
C
      Call PINOVEC(EigVecR,Eig,ABPLUS,ABMIN,DMAT,DMATK,EMAT,
     $ EMATM,Occ,NBasis,NDimX,NDimN)
C
      Write(6,'(" *********************
     $  END OF STABILITY TEST *****************",//)')
C
C ********************************************************************************
  777 Continue
C
      Do I=1,NBasis
      CISave(I)=CICoef(I)
      IGemS(I)=IGem(I)
      EndDo 
      NGemS=NGem
C
C     CONSTRUCT A LOOK-UP TABLE
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      ICount=0
      Do I=1,NBasis
c      If(Occ(I).Gt.Zero) Then
      IndAux(I)=1
      ICount=ICount+1
c      EndIf
      EndDo
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
C
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate active degenerate orbitals if from different geminals
      If((IGem(I).Ne.IGem(J)).And.(IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-20) ) Then
C
      Write(*,*)"Discarding nearly degenerate pair",I,J
C
      Else
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C      
      EndIf
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
      NDimN=NBasis
c herer!!!
c      NDimN=0
C
C     GENERATE ABSCISSAS AND WEIGHTS FOR GAUSSIAN-LEGENDRE QUADRATURE
C
      NGrid=30
C
      Call GauLeg(Zero,One,XGrid,WGrid,NGrid)
C
c      GoTo 444
c for tests
c      ACAlpha=One
c      Do I=1,NInte2
c      TwoAlph(I)=TwoEl(I)
c      EndDo
c      NoEig=1
c      Call OptTwoATrip(TrGamm,GammAl,EExcit,ETot,ENuc,Occ,XOne,TwoAlph,
c     $ H1Alph,UReAlph,OccAlph,
c     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
c       write(*,*)eexcit
c       stop
C
C     ******************************************
C     THIS SHOULD WORK FOR TRIPLETS AND SINGLETS
C     ***************************************** 
C
      ECorr=Zero
      Delta=Zero
      Do N=1,NGrid
C
      ACAlpha=XGrid(N)
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      Call OptTwoATrip(TrGamm,GammAl,EExcit,ETot,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
       Call ACEneST(ECorrA,DeltaA,TrGamm,GammAl,EExcit,
     $ TwoEl,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NDimN)
      Write(6,'(/,X,''ACAlpha:'',F11.8,3X,"W_ALPHA:",F12.8,
     $ 3X,"Delta_ALPHA:",F12.8)')ACAlpha,ECorrA,DeltaA
C
      ECorr=ECorr+WGrid(N)*ECorrA
      Delta=Delta+WGrid(N)*DeltaA
C
      EndDo
C
C     INTEGRATION WRT ALPHA FINNISHED
C     CALL OptTwoAlpha with ACAlpha=0.0 to find the reference energy
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      ACAlpha=Zero
      Write(6,'(/,X,''********* FIND THE REFERENCE ENERGY (ACAlpha ='',
     $ F5.2,") *********")') ACAlpha
C
      Call OptTwoATrip(TrGamm,GammAl,EExcit,ETot,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
      Write
     $ (6,'(/,X,''ERef+ENuc, W, Delta, Total'',4F15.8)')
     $ ETot+ENuc,ECorr,Delta,ETot+ENuc+ECorr+Delta 
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      ACAlpha=One
      Write(6,'(/,X,''********* FIND THE EXACT ENERGY (ACAlpha ='',
     $ F5.2,") *********")') ACAlpha
C
      Call OptTwoATrip(TrGamm,GammAl,EExcit,ETot,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
c herer!!! (uncomment the line below after tests)
      Return
C
  444 Continue
C
C     *************************************
C     WHAT FOLLOWS WORKS ONLY FOR SINGLETS
C     *************************************
C
      ECorr=Zero
      Delta=Zero
      Do N=1,NGrid
C
      ACAlpha=XGrid(N)
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      Call OptTwoAlpha(ETotAlph,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha) 
      NGem=1
      Do II=1,NBasis
      CICoef(II)=CIAlph(II)
      IGem(II)=1
      If(OccAlph(II).Eq.Zero) IGem(II)=2
      If(OccAlph(II).Eq.Zero) NGem=2
      EndDo
C
      Call APSG_NEST(ABPLUS,ABMIN,CMAT,EMAT,EMATM,DMAT,DMATK,
     $ UReAlph,OccAlph,H1Alph,TwoAlph,
     $ NBasis,NDim,NInte1,NInte2,NGem,2) 
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1)
      EndDo
      EndDo
      Do J=1,NDimN
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(J-1)*NDim+IndX(I)
      DMAT(IJ)=DMAT(IJ1)
      DMATK(IJ)=DMATK(IJ1)
      EndDo
      EndDo
      Do J=1,NDimN
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(J-1)*NBasis+I
      EMAT(IJ)=EMAT(IJ1)
      EMATM(IJ)=EMATM(IJ1)
      EndDo
      EndDo
C
      Call PINOVEC(EigVecR,Eig,ABPLUS,ABMIN,DMAT,DMATK,EMAT,
     $ EMATM,Occ,NBasis,NDimX,NDimN)       
C
      Do I=1,NBasis
      CICoef(I)=CISave(I)
      IGem(I)=IGemS(I)
      EndDo
      NGem=NGemS
C
      Call ACEnePINO(ECorrA,DeltaA,
     $ EigVecR,Eig,TwoEl,Occ,XOne,
     $ CIAlph,UReAlph,IndN,NBasis,NInte1,NInte2,NDimX,NDimN,1)
C
      Write(6,'(/,X,''ACAlpha:'',F11.8,3X,"W_ALPHA:",F12.8,
     $ 3X,"Delta_ALPHA:",F12.8)')ACAlpha,ECorrA,DeltaA
C
      ECorr=ECorr+WGrid(N)*ECorrA
      Delta=Delta+WGrid(N)*DeltaA
C
      EndDo
C
C     INTEGRATION WRT ALPHA FINNISHED
C     CALL OptTwoAlpha with ACAlpha=0.0 to find the reference energy
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      ACAlpha=Zero
      Write(6,'(/,X,''********* FIND THE REFERENCE ENERGY (ACAlpha ='',
     $ F5.2,") *********")') ACAlpha
      Call OptTwoAlpha(ETot,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
      Write
     $ (6,'(/,X,''ERef+ENuc, W, Delta, Total'',4F15.8)')
     $ ETot+ENuc,ECorr,Delta,ETot+ENuc+ECorr+Delta
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      ACAlpha=One
      Write(6,'(/,X,''********* FIND THE EXACT ENERGY (ACAlpha ='',
     $ F5.2,") *********")') ACAlpha
      Call OptTwoAlpha(ETot,ENuc,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
      Return 
      End

*Deck ACRDM
      Subroutine ACRDM(ETot,ENuc,TwoEl,Occ,XOne,
     $ UNOAO,IndN,IndX,IndAux,NDimX,NBasis,NInte1,NInte2,NDim,NGem)
C
C     AC ENERGY USING ERPA TRANSITION DENSITY MATRIX ELEMENTS
C     AND ALPHA-DEPENDENT 1,2-RDMs
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ Occ(NBasis),TwoEl(NInte2),XOne(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ UReAlph(NBasis,NBasis),OccAlph(NBasis),CIAlph(NBasis),
     $ TwoAlph(NInte2),H1Alph(NInte1),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ CMAT(NDim*NDim),EMAT(NBasis*NBasis),DMAT(NDim*NBasis),
     $ EMATM(NBasis*NBasis),DMATK(NDim*NBasis),
     $ IndX(NDim),IndN(2,NDim),
     $ EigVecR(2*(NDim+NBasis)*2*(NDim+NBasis)),
     $ Eig(2*(NDim+NBasis)),IndAux(NBasis),
     $ XGrid(100), WGrid(100),CISave(NBasis),IGemCAS(NBasis),
     $ TrGamm(NInte1,NInte1),EExcit(NInte1),GammAl(NInte1),
     $ IPair(NBasis,NBasis)
C
      Character*10     :: IntFile 
C
c      UReAlph(1:NBasis,1:NBasis)=Zero
c      Do I=1,NBasis
c      UReAlph(I,I)=One
c      EndDo
c      Call MP2RDM(TwoEl,EMAT,Occ,UReAlph,UNOAO,XOne,
c     $           IndN,IndX,IndAux,NDimX,
c     $           NBasis,NDim,NInte1,NInte2,NVirt,
c     $           IntFile,ThrVirt,.true.)
c      stop
C
      Do I=1,NBasis
      CISave(I)=CICoef(I)
      IGemCAS(I)=IGem(I)
      EndDo
      NGemCAS=NGem
C
      NDimN=0
C
C      NGrid=30
C      Call GauLeg(Zero,One,XGrid,WGrid,NGrid)
      Open(10,File="alpha.txt",Status='Old')
      IX=0
   10 Read(10,*,End=40)XGrid(IX+1)
      IX=IX+1
      GoTo 10
   40 Continue
      Close(10)      
      Open(10,File="weights.txt",Status='Old')
      IW=0
   12 Read(10,*,End=42)WGrid(IW+1)
      IW=IW+1
      GoTo 12
   42 Continue
      Close(10) 
      If(IX.Eq.IW) Then
      NGrid=IX
      Else
      Stop 'Fatal Error: Inconsistent NGrid for XGrid and WGrid'
      EndIf
C
      ECorr=Zero
      Delta=Zero
C
      Do N=1,NGrid
C
      ACAlpha=XGrid(N)
C
      Do I=1,NInte2
      TwoAlph(I)=TwoEl(I)
      EndDo
      Call RDMAlpha(IndX,IndN,IndAux,NDimX,NDim,Occ,XOne,TwoAlph,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ IPair,IGemCAS,
     $ NBasis,NInte1,NInte2,ACAlpha,N)
C
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,UReAlph,OccAlph,H1Alph,TwoAlph,
     $ IPair,IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,1.D0)
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1)
      EndDo
      EndDo
C
c herer!!!
c      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
      Do I=1,NBasis
      CICoef(I)=CISave(I)
      IGem(I)=IGemCAS(I)
      EndDo
      NGem=NGemCAS
C
      Call ACEnePINO(ECorrA,DeltaA,
     $ EigVecR,Eig,TwoEl,Occ,XOne,
     $ CIAlph,UReAlph,IndN,NBasis,NInte1,NInte2,NDimX,NDimN,2)
C
      Write(6,'(/,X,''ACAlpha:'',F11.8,3X,"W_ALPHA:",F12.8,
     $ 3X,"Delta_ALPHA:",F12.8)')ACAlpha,ECorrA,DeltaA
C
      ECorr=ECorr+WGrid(N)*ECorrA
      Delta=Delta+WGrid(N)*DeltaA
C
      EndDo
C
      Write
     $ (6,'(/,X,''ERef+ENuc, W, Delta, Total'',4F15.8)')
     $ ETot+ENuc,ECorr,Delta,ETot+ENuc+ECorr+Delta
C
      Return
      End

*Deck RDMAlpha
      Subroutine RDMAlpha(IndX,IndN,IndAux,NDimX,NDim,Occ,XOne,TwoEl,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ IPair,IGemCAS,
     $ NBasis,NInte1,NInte2,ACAlpha,IGrid)
C
C     For a given alpha: construct/read 1,2-RDM(Alpha)
C     find a transformation matix from CAS MO to NO(Alpha)
C     construct Alpha-dependent Hamiltonian elements and transform 
C     to NO(Alpha)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension UReAlph(Nbasis,NBasis),Occ(NBasis),
     $ OccAlph(NBasis),CIAlph(NBasis),XOne(NInte1),TwoEl(NInte2),
     $ H1Alph(NInte1),
     $ IndX(NDim),IndN(2,NDim),IndAux(NBasis),IGemCAS(NBasis),
     $ IPair(NBasis,NBasis)
C
C     read 1-,2-RDMs in CAS orbital representation
C
      Call ReadRDMs(OccAlph,UReAlph,NBasis,NGem,NInte1,IGrid) 
      Do I=1,NBasis
      CIAlph(I)=SQRT(OccAlph(I))
      If(OccAlph(I).Lt.Half) CIAlph(I)=-CIAlph(I)
      EndDo
C
c herer!!! 
      GoTo 222
C
C     CONSTRUCT A LOOK-UP TABLE
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      ICount=0
      Do I=1,NBasis
      If(OccAlph(I).Lt.One.And.OccAlph(I).Ne.Zero) Then
      IndAux(I)=1
c      Write(6,'(X," Active Orbital: ",I4,E14.4)') I, OccAlph(I)
      ICount=ICount+1
      EndIf
      EndDo
C
      IPair(1:NBasis,1:NBasis)=0
C 
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
     $ .And.(Abs(OccAlph(I)-OccAlph(J))/OccAlph(I).Lt.ThrSelAct)
     $ ) Then
C
      Write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.OccAlph(I).Ne.One.And.OccAlph(J).Ne.One)) Then
C
      If(Abs(OccAlph(i)+OccAlph(j)-Two).Gt.1.D-10.And.
     $   Abs(OccAlph(i)+OccAlph(j)).Gt.ThrQVirt) Then
      Ind=Ind+1
      IndX(Ind)=Ind
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      IPair(I,J)=1
      IPair(J,I)=1
      EndIf
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
c herer!!!
  222 Continue
      Do Ind=1,NDimX
      I=IndN(1,Ind)
      J=IndN(2,Ind)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
      Write(6,'(/,X,"NDimX = ",I4)')NDimX
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      Do I=1,NInte1
      H1Alph(I)=XOne(I)
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGemCAS(I).Ne.IGemCAS(J)) Then
C
      H1Alph(IJ)=ACAlpha*H1Alph(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGemCAS(IT).Ne.IGemCAS(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      H1Alph(IJ)=H1Alph(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C     AND TRANSFORM TO NO(Alpha) REPRESENTATION
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
      If(.Not.(
     $IGemCAS(I).Eq.IGemCAS(J).And.IGemCAS(J).Eq.IGemCAS(K).And.
     $ IGemCAS(K).Eq.IGemCAS(L)))
     $ TwoEl(NAdd)=ACAlpha*TwoEl(NAdd)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call TwoNO1(TwoEl,UReAlph,NBasis,NInte2) 
C
      Return
      End 

*Deck ReadRDMs
      Subroutine ReadRDMs(Occ,URe,NBasis,NGem,NInte1,IGrid)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),UAux(NBasis*NBasis),
     $ UCASNO(Nbasis,NBasis)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2(:)
      Dimension Gamma(NInte1),Work(NBasis),PC(NBasis),
     $ AUXM(NBasis,NBasis),AUXM1(NBasis,NBasis)
      Character*15 FName,StrNum
C
      UAux(1:NBasis*NBasis)=Zero
      URe(1:NBasis,1:NBasis)=Zero
      UCASNO(1:NBasis,1:NBasis)=Zero
      Occ(1:NBasis)=Zero
      PC(1:NBasis)=Zero
      Gamma(1:NInte1)=Zero
C
C     FOR A GIVEN ALPHA READ IN 1-RDM AND DIAGONALIZE IT
C
      If(IGrid-1.Lt.10) L=1
      If(IGrid-1.Ge.10) L=2
      Write(StrNum,*)IGrid-1
C
      FName(1:7+L)="G1_"//StrNum(13-L:12)//".bin"
      write(*,*)'G1 file name ',FName(1:7+L)
      Open(10,File=FName(1:7+L),form='unformatted', access='stream',
c       Open(10,File='G1_0.bin',form='unformatted', access='stream',
     $ Status='Old')
      Read(10)i,j,k
      ICount=0
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Read(10,End=61) X
      Ind=I*(I-1)/2+J
      Gamma(Ind)=X/Two 
      ICount=ICount+1
      EndDo
      EndDo
   61 Close(10)
      Call CpySym(AUXM,Gamma,NBasis)
      Call Diag8(AUXM,NBasis,NBasis,PC,Work)
      Call SortOcc(PC,AUXM,NBasis)
C
      Sum=Zero
      NAc=0
      Do I=1,NBasis
      Sum=Sum+PC(I)
      PC(I)=Abs(PC(I))
      If(PC(I).Gt.Zero) NAc=NAc+1
      EndDo
C
      NInAc=NELE-Sum+1.D-1
      Do I=1,NInAc+NAc
      If(I.Le.NInAc) Then
      Occ(I)=One
      Else
      Occ(I)=PC(I-NInAc)
      EndIf
      EndDo
C
      If(NInAc.Eq.0) Then
      NGem=2
      IGem(1:NInAc+NAc)=1
      IGem(NInAc+NAc+1:NBasis)=2
      Else
      NGem=3
      IGem(1:NInAc)=1
      IGem(NInAc+1:NInAc+NAc)=2
      IGem(NInAc+NAc+1:NBasis)=3
      EndIf
C
      NAcCAS=NAc
      NInAcCAS=NInAc
C
      Write(6,'(2X,"No of DMRG inactive and active orbitals: ",2I4)')
     $ NInAcCAS,NAcCAS
C
      Write(6,'(2X,"DMRG",3X,"Occupancy",4X,"Gem")')
      Sum=Zero
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,I6)') I,Occ(I),IGem(I)
      Sum=Sum+Occ(I)
      EndDo
      Write(6,'(2X,"Sum of Occupancies: ",F5.2)') Sum
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
C
C     COPY AUXM TO URe AND OFF SET BY NInAc
      Do I=1,NBasis
      IIAct=I-NInAc
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      JJAct=J-NInAc
      If(IIAct.Gt.0.And.IIAct.Le.NAc.And.JJAct.Gt.0.And.JJAct.Le.NAc)
     $ Then
      URe(I,J)=AUXM(IIAct,JJAct)
      EndIf
      EndDo
      EndDo
C
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      Allocate (RDM2(NRDM2))
      RDM2(1:NRDM2)=Zero
C
      FName(1:7+L)="G2_"//StrNum(13-L:12)//".bin"
      write(*,*)'G2 file name ',FName(1:7+L)
      Open(10,File=FName(1:7+L),form='unformatted',access='stream',
c      Open(10,File='G2_0.bin',form='unformatted',access='stream',
     $ Status='Old')
      Read(10)I,J,K
C
      Do K1=1,NBasis
      Do J1=1,NBasis
      Do L1=1,NBasis
      Do I1=1,NBasis
C
      Read(10) X
C
      I=I1+NInAc
      J=J1+NInAc
      K=K1+NInAc
      L=L1+NInAc
      RDM2(NAddrRDM(L,K,I,J,NBasis))=X
      RDM2(NAddrRDM(K,L,J,I,NBasis))=X
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Close(10)
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAux((J-1)*NBasis+I)=AUXM(I,J)
      EndDo
      EndDo
C
      Call TrRDM2(RDM2,UAux,NBasis,NRDM2)
C     SAVE THE ACTIVE PART IN rdm2.dat
      Open(10,File='rdm2.dat')
      Do I=NInAc+1,NInAc+NAc
      IIAct=I-NInAc
      Do J=NInAc+1,NInAc+NAc
      JJAct=J-NInAc
      IJ=(I-1)*NBasis+J
      Do K=NInAc+1,NInAc+NAc
      KKAct=K-NInAc
      Do L=NInAc+1,NInAc+NAc
      LLAct=L-NInAc
      KL=(K-1)*NBasis+L
      If(IJ.Ge.KL) Write(10,'(4I4,F19.12)')
     $ KKAct,IIAct,LLAct,JJAct,Two*RDM2(NAddrRDM(I,J,K,L,NBasis))
      EndDo
      EndDo
      EndDo
      EndDo
      Close(10)
      Deallocate(RDM2)
C
C
C     CONSTRUCT A TRANSFORMATION MATRIX FROM NO-CAS to THE NEW ALPHA-ORBITALS
C
C     TRANSFORMATION MATRIX FROM CAS TO NO_CAS
C
      open(10,file='ure_casno.dat')
      read(10,*)UCASNO
      close(10)
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAux((J-1)*NBasis+I)=URe(I,J)
      EndDo
      EndDo 
C
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=Zero
      Do K=1,NBasis
      URe(I,J)=URe(I,J)+UAux((K-1)*NBasis+I)*UCASNO(J,K)
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck OptTwoAlpha
      Subroutine OptTwoAlpha(ETotAlph,ENuc,Occ,XOne,TwoEl,
     $ H1Alph,UReAlph,OccAlph,CIAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SYSTEMS
C     THE WHOLE DENISTY MATRIX IS FOUND IN ONE STEP
C
C     ON EXIT TwoEl INCLUDES TWO-ELECTRON INTEGRALS WITH ALPHA 
C     TRANSFORMED BY UReAlph
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension UReAlph(Nbasis,NBasis),Occ(NBasis),
     $ OccAlph(NBasis),CIAlph(NBasis),XOne(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $ H1Alph(NInte1),HNO(NInte1)
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      Do I=1,NInte1
      H1Alph(I)=XOne(I)
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      H1Alph(IJ)=ACAlpha*H1Alph(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      H1Alph(IJ)=H1Alph(IJ)+Aux
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
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ TwoEl(NAdd)=ACAlpha*TwoEl(NAdd)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))
     $                +TwoEl(NAddr3(IA,ID,IC,IB))
C
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETotAlph=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(X,''State no'',I3,'' Total Energy'',F15.8)') NoEig,
     $ ETotAlph+ENuc
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(UReAlph,PWork,NBasis)
      Call Diag8(UReAlph,NBasis,NBasis,OccAlph,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      CIAlph(I)=OccAlph(I)
      OccAlph(I)=OccAlph(I)**2
      Sum=Sum+OccAlph(I)
      EndDo
      Do I=1,NBasis
      OccAlph(I)=OccAlph(I)/Sum
      CIAlph(I)=CIAlph(I)/Sqrt(Sum)
      EndDo
C
      Call SortAlph(OccAlph,CIAlph,UReAlph,NBasis)
C
      Write(6,'('' ORBITAL OCCUPANCIES AND COEFFICIENTS'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,2E16.6)') I,OccAlph(I),CIAlph(I)
C
      Call TwoNO1(TwoEl,UReAlph,NBasis,NInte2) 
C
      Return
C
C     COMPUTE ONE-ELECTRON AND TWO-ELECTRON COMPONENTS OF THE ENERGY
C     (FOR CHECKING MAINLY) 
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HNO(IJ)=HNO(IJ)+UReAlph(I,IA)*UReAlph(J,IB)*H1Alph(IAB)
      EndDo
      EndDo
      EndDo
      EndDo
C
      EOne=Zero
      ETwo=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      EOne=EOne+Two*OccAlph(I)*HNO(II)
      Do J=1,NBasis
      ETwo=ETwo+CIAlph(I)*CIAlph(J)*TwoEl(NAddr3(I,J,I,J))
      EndDo
      EndDo 
      Write(6,'(X,"EOne, ETwo, ETot",3F15.8)') EOne,ETwo,EOne+ETwo+ENuc
C
      Return
      End

*Deck OptTwoATrip
      Subroutine OptTwoATrip(TrGamm,GammAl,
     $ EExcit,ETotAlph,ENuc,Occ,XOne,TwoEl,
     $ H1Alph,UReAlph,OccAlph,
     $ NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SINGLET OR TRIPLET SYSTEMS 
C     ISpin=0 or ISpin=1
C
C     ON EXIT TRANSITION 1-RDM's ARE RETURNED:
C     Gamma_pq+Gamm_qp for p.ne.q 
C     and Gamma_pp
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension UReAlph(Nbasis,NBasis),Occ(NBasis),
     $ OccAlph(NBasis),CIAlph(NBasis),XOne(NInte1),TwoEl(NInte2),
     $ TrGamm(NInte1,NInte1),GammAl(NInte1),EExcit(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $ H1Alph(NInte1),HNO(NInte1),PFull(NBasis,NBasis,NInte1),
     $ ISkip(NInte1)
C
C change to ISpin=0 for singlets or ISpin=1 for triplets
c herer!!!
c      ISpin=0
      ISpin=1
C
      SpinFac=One
      If(ISpin.Eq.1) SpinFac=-One
C
      ISkip(1:NInte1)=0
C
C     CORE ENERGY (IF ANY)
C
      NCore=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.1.And.Occ(I).Eq.One) NCore=NCore+1
      EndDo
      ECore=Zero
      Do I=1,NCore
      II=(I*(I+1))/2
      ECore=ECore+Two*Occ(I)*XOne(II)
      EndDo
      Do IP=1,NCore
      Do IQ=1,NCore
      ECore=ECore+
     $ Two*TwoEl(NAddr3(IP,IP,IQ,IQ))-TwoEl(NAddr3(IP,IQ,IP,IQ))
      EndDo
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      Do I=1,NInte1
      H1Alph(I)=XOne(I)
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
C     ADD INTERACTION WITH CORE
C
      Do IT=1,NCore
      H1Alph(IJ)=H1Alph(IJ)+
     $ Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J))
      EndDo      
C
      If(IGem(I).Ne.IGem(J)) Then
C
      H1Alph(IJ)=ACAlpha*H1Alph(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=NCore+1,NBasis
      If(IGem(IT).Ne.IGem(I)) 
     $  Aux=Aux+(One-ACAlpha)*Occ(IT)*
     $ (Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J)))
      Enddo
C
      H1Alph(IJ)=H1Alph(IJ)+Aux
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
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ TwoEl(NAdd)=ACAlpha*TwoEl(NAdd)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB.And.ISpin.Eq.1) FAB=Zero
      If(IA.Eq.IB.And.ISpin.Eq.0) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID.And.ISpin.Eq.1) FCD=Zero
      If(IC.Eq.ID.And.ISpin.Eq.0) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))
     $               +SpinFac*TwoEl(NAddr3(IA,ID,IC,IB))
C
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+SpinFac*H1Alph(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+SpinFac*H1Alph(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
C     If core is present, put to zero the corresponding elements 
C     of AP to remove core from diagonalization
C
      If(IA.Le.NCore.Or.IB.Le.NCore.Or.IC.Le.NCore.Or.ID.Le.NCore)
     $ AP(IAB,ICD)=Zero
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      write(*,*)'beginning of diag8',ninte1
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      write(*,*)'1','end of diag8' 
      ETotAlph=PP(NoEig)
C
      Do I=1,NInte1
      EExcit(I)=PP(I)-PP(NoEig)
C
C     because the core-core block (if any) is set to 0, 
C     almost-exact-zero-eigenvalues appear, set them exactly to 0
C     to exclude the pertinent eigenvectors from the correlation energy
C
c herer!!!
c      If(Abs(PP(I)).Lt.1.D-10) EExcit(I)=Zero
C     do not account for zero excitation corresponding to noeig
c      If(Abs(PP(I)-PP(NoEig)).Lt.1.D-10) EExcit(I)=Zero
      If(I.Eq.NoEig) EExcit(I)=Zero 
c herer!!! this is for test only
c      if(EExcit(I).lt.zero)EExcit(I)=0
      EndDo 

      Do I=1,12
      write(*,*)i,'ene',pp(i),pp(i)+ecore+enuc
      enddo
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(X,''ACAlpha:'',F11.8,3X,"ETwoEl:",F12.8)')
     $ACAlpha,ETotAlph
      ETotAlph=ETotAlph+ECore
      If(ECore.Ne.Zero) Write(6,'(" Core Energy: ",F15.8)')ECore
      Write(6,'(X,''State no'',I3,'' Total Energy'',F15.8)') NoEig,
     $ ETotAlph+ENuc
C
      Do I=1,NInte1
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      FAB=One
      If(IA.Ne.IB) FAB=SQRT(Half)
      PFull(IA,IB,I)=FAB*AP(I,IAB)
      If(IA.Ne.IB) PFull(IB,IA,I)=FAB*SpinFac*AP(I,IAB)
C
C in addition to proper triplet state, some other eigenvectors (junk) are found 
C this must be removed by checking if the diag element (remember that for triplet P_pq=-P_qp => P_pp=0)
C is not zero (in practice if gt 1.-5)
      If(IA.Eq.IB.And.ISpin.Eq.1.And.Abs(PFull(IB,IA,I)).Gt.1D-5) Then 
      ISkip(I)=1
      EndIf
      EndDo
      EndDo
C
      Sum=Zero
      Do IA=1,NBasis
      Do IB=1,NBasis
      Sum=Sum+PFull(IA,IB,I)**2
      EndDo
      EndDo
      If(Abs(Sum-One).Gt.1.D-10)
     $  Stop 'Error: wrong normalization of P matrix'
C
      EndDo
C
      If(ISpin.Eq.1) Then
C
      Do I=1,NInte1
      If(ISkip(I).Eq.1) Then
      EExcit(I)=Zero
      Do IA=1,NBasis
      Do IB=1,IA
      PFull(IA,IB,I)=Zero
      PFull(IB,IA,I)=Zero
      EndDo
      EndDo
      EndIf
      EndDo
C
      EndIf
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      UReAlph(I,J)=Zero
      Do K=1,NBasis
      UReAlph(I,J)=UReAlph(I,J)+PFull(I,K,NoEig)*PFull(J,K,NoEig)
      EndDo
      UReAlph(J,I)=UReAlph(I,J)
      GammAl(IJ)=UReAlph(I,J)
      EndDo
      EndDo
C
      Call Diag8(UReAlph,NBasis,NBasis,OccAlph,PWork)
C
      Write(6,'('' EXACT ORBITAL OCCUPANCIES FOR THE NoEig STATE'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,E16.6)') I,OccAlph(I)
C
      Do INU=1,NInte1
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      TrGamm(IAB,INU)=Zero
C
      Do I=1,NBasis
      If(IA.Ne.IB) TrGamm(IAB,INU)=TrGamm(IAB,INU)
     $ +PFull(IA,I,NoEig)*PFull(IB,I,INU)
     $ +PFull(IB,I,NoEig)*PFull(IA,I,INU)
      If(IA.Eq.IB) TrGamm(IAB,INU)=TrGamm(IAB,INU)
     $ +PFull(IA,I,NoEig)*PFull(IA,I,INU)
      EndDo
C
      EndDo
      EndDo
C
      EndDo
C
      Return
      End

*Deck SortAlph
      Subroutine SortAlph(Occ,CI,URe,NBasis) 
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Occ(NBasis),CI(NBasis),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension UReOld(NBasis,NBasis),Ind(1000)
C
      If(NBasis.Gt.1000) Stop 'SortOcc: NBasis > 1000'
C
C     SORT THE OCCUPATION NUMBERS IN A DESCENDING ORDER
C
      Do I=1,NBasis
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,NBasis
C
      OccMx=Occ(IStart)
      IndMx=IStart
C
      Do J=IStart,NBasis
      If(Occ(J).Gt.OccMx) Then
      OccMx=Occ(J)
      IndMx=J
      EndIf
      EndDo
C
      Hlp=Occ(IStart)
      IndHlp=Ind(IStart)

      Occ(IStart)=OccMx
      Ind(IStart)=Ind(IndMx)

      Occ(IndMx)=Hlp
      Ind(IndMx)=IndHlp
C
      IStart=IStart+1
C
      EndDo
C
C     SWAP THE ORBITALS
C
      Do I=1,NBasis
      Do J=1,NBasis
      UReOld(J,I)=URe(J,I)
      EndDo
      EndDo
C
      Do J=1,NBasis
      Do I=1,NBasis
      URe(I,J)=UReOld(Ind(I),J)
      EndDo
      EndDo
C
      Do I=1,NBasis
      UreOld(1,I)=CI(I)
      EndDo
      Do I=1,NBasis
      CI(I)=UreOld(1,Ind(I))
      EndDo
C
      Return
      End

*Deck ACEnePINO
      Subroutine ACEnePINO(ECorr,Delta,
     $ EigVecR,Eig,TwoNO,Occ,
     $ HNO,CIAlph,UReAlph,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NDimN,IF12)
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(SmallE=1.D-4,BigE=1.D20)
C
C     ONLY EXCITATIONS > SmallE AND < BigE ARE INCLUDED 
C     
      Include 'commons.inc'
C     
      Dimension EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Eig(2*(NDimX+NDimN)),
     $ Occ(NBasis),TwoNO(NInte2),
     $ CIAlph(NBasis),HNO(NInte1),
     $ UReAlph(NBasis,NBasis),
     $ IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension OccAlph(NBasis),GammaAlph(NInte1),
     $ EigVecTr(2*(NDimX+NDimN)*2*(NDimX+NDimN))
C
      If(IF12.Eq.1) NI=2*(NDimX+NDimN)
      If(IF12.Eq.2) NI=NDimX+NDimN
      EigVecTr(1:NI*NI)=Zero 
C
      Do I=1,NBasis
      OccAlph(I)=CIAlph(I)**2
      EndDo
C
      Call VecTr(GammaAlph,OccAlph,UReAlph,NBasis) 
C
C     TRANSFORM EIGENVECTORS TO THE NATURAL ORBITALS OF THE REFERENCE WF
C
      Do NU=1,NI
      If(Eig(NU).Gt.Zero.And.Eig(NU).Le.SmallE) 
     $ Write(*,*)'Small omega=',Eig(NU),
     $ ' The corresponding eigenvector excluded from Ecorr'
c      If(Eig(NU).Gt.SmallE.And.Eig(NU).Lt.BigE) Then
      If(Abs(Eig(NU)).Gt.Zero) Then
C
      Do I=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
      IP=IndN(1,I)
      IQ=IndN(2,I)
      Else
      IP=I-NDimX
      IQ=IP
      EndIf
C
      GG=Zero
C
      Do J=1,NDimX+NDimN
C
      If(J.Le.NDimX) Then
      IR=IndN(1,J)
      IS=IndN(2,J)
      Else
      IR=J-NDimX
      IS=IR
      EndIf
C
      If(IP.Ne.IQ) Then
C
      If(IR.Ne.IS) GG=GG+
     $ (UReAlph(IR,IP)*UReAlph(IS,IQ)+UReAlph(IR,IQ)*UReAlph(IS,IP))
     $ *(CIAlph(IR)+CIAlph(IS))*EigVecR((NU-1)*NI+J)
      If(IR.Eq.IS) GG=GG+Two*UReAlph(IR,IP)*UReAlph(IR,IQ)*CIAlph(IR)
     $ *EigVecR((NU-1)*NI+NDimX+J)
C
      EndIf
C
      If(IP.Eq.IQ) Then
C
      If(IR.Ne.IS) GG=GG+UReAlph(IR,IP)*UReAlph(IS,IP)
     $ *(CIAlph(IR)+CIAlph(IS))*EigVecR((NU-1)*NI+J)
      If(IR.Eq.IS) GG=GG+UReAlph(IR,IP)**2*CIAlph(IR)
     $ *EigVecR((NU-1)*NI+NDimX+J)
C
      EndIf
C
C     Do J
      EndDo
C
      If(IP.Ne.IQ) EigVecTr((NU-1)*NI+I)=GG
      If(IP.Eq.IQ) EigVecTr((NU-1)*NI+NDimX+I)=GG
C 
C     Do I
      EndDo
C 
C     end of the loop w.r.t NU
C
      EndIf
C
      EndDo
C
      ECorr=Zero
      Do I=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
      IP=IndN(1,I)
      IR=IndN(2,I)
      Else
      IP=I-NDimX
      IR=IP
      EndIf
C
      Do J=1,NDimX+NDimN
C
      If(J.Le.NDimX) Then
      IQ=IndN(1,J)
      IS=IndN(2,J)
      Else
      IQ=J-NDimX
      IS=IQ
      EndIf
C
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))) Then
C
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      SumY=Zero
      Do K=1,NI
c      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE) Then
      SumY=SumY+EigVecTr((K-1)*NI+I)*EigVecTr((K-1)*NI+J)
c      EndIf
      EndDo
C
      Aux=Two*SumY
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux=Aux+
     $ (-Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP)) )
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
      If(IP.Gt.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NI
c      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE) Then
      SumY=SumY+EigVecTr((K-1)*NI+I)*EigVecTr((K-1)*NI+NDimX+J)
c      EndIf
      EndDo
C
      ECorr=ECorr+Four*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      EndIf
C
      If(IP.Eq.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NI
c      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE) Then
      SumY=SumY+EigVecTr((K-1)*NI+NDimX+I)*EigVecTr((K-1)*NI+NDimX+J)
c      EndIf
      EndDo
C
      ECorr=ECorr+Two*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))
C
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
C     COMPUTE DELTA
C
      Delta=Zero
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))) Then
C
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      Aux=Two*GammaAlph(IPR)*GammaAlph(IQS)
      IQR=(Max(IR,IQ)*(Max(IR,IQ)-1))/2+Min(IR,IQ)
      If(IP.Eq.IS) Aux=Aux-GammaAlph(IQR)
      If(IP.Eq.IS.And.IQ.Eq.IR) Aux=Aux+Occ(IQ)
      If(IP.Eq.IR.And.IQ.Eq.IS) Aux=Aux-Two*Occ(IP)*Occ(IQ)
C
      Aux=Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      Delta=Delta+Aux
C
      EndIf 
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      If((IGem(IP).Eq.IGem(IQ)).And.(IGem(IP).Ne.IGem(IR))) Then
C
      IPQ=(Max(IQ,IP)*(Max(IQ,IP)-1))/2+Min(IQ,IP)
      Aux=-GammaAlph(IPQ)
      If(IP.Eq.IQ) Aux=Aux+Occ(IP)
      Aux=Two*Occ(IR)*Aux*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IR))-TwoNO(NAddr3(IP,IR,IQ,IR)))
      Delta=Delta+Aux
C
      EndIf
      EndDo
      EndDo
      EndDo
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      If(IGem(IP).Ne.IGem(IQ)) Then
C
      IPQ=(Max(IQ,IP)*(Max(IQ,IP)-1))/2+Min(IQ,IP)
      Delta=Delta+Two*HNO(IPQ)*GammaAlph(IPQ) 
C
      EndIf
      EndDo
      EndDo
C
      Return
      End 

*Deck ACEneST
      Subroutine ACEneST(ECorr,Delta,TrGamm,GammaAlph,EExcit,
     $ TwoNO,Occ,HNO,IndN,NBasis,NInte1,NInte2,NDimX,NDimN)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension 
     $ Occ(NBasis),TwoNO(NInte2),
     $ HNO(NInte1),
     $ IndN(2,NDimX),
     $ TrGamm(NInte1,NInte1),EExcit(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension GammaAlph(NInte1),H1(NInte1),IOCore(NBasis)
C
      NI=2*(NDimX+NDimN)
C
C     ADD INTERACTON WITH CORE TO HNO
C
      NCore=Zero
      Do I=1,NBasis
      IOCore(I)=0
      If(IGem(I).Eq.1.And.Occ(I).Eq.One) Then
      NCore=NCore+1
      IOCore(I)=1
      EndIf
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      H1(IJ)=HNO(IJ)
      Do IT=1,NCore
      H1(IJ)=H1(IJ)+
     $ Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J))
      EndDo 
      EndDo
      EndDo
C
      ECorr=Zero
      Do I=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
      IP=IndN(1,I)
      IR=IndN(2,I)
      Else
      IP=I-NDimX
      IR=IP
      EndIf
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Do J=1,NDimX+NDimN
C
      If(J.Le.NDimX) Then
      IQ=IndN(1,J)
      IS=IndN(2,J)
      Else
      IQ=J-NDimX
      IS=IQ
      EndIf
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
C
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))
     $ .And.(IOCore(IP)+IOCore(IR)+IOCore(IQ)+IOCore(IS).Eq.0) ) Then
C
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      SumY=Zero
      Do K=1,NInte1
c herer!!!
      If(Abs(EExcit(K)).Gt.1.D-8) Then
c      If(EExcit(K).Lt.Zero) Then
      SumY=SumY+TrGamm(IPR,K)*TrGamm(IQS,K)
      EndIf
      EndDo
C
      Aux=Two*SumY
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux=Aux+
     $ (-Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP)) )
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
      If(IP.Gt.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NInte1
c herer!!!
      If(Abs(EExcit(K)).Gt.1.D-8) Then
c      If(EExcit(K).Lt.Zero) Then
      SumY=SumY+TrGamm(IPR,K)*TrGamm(IQS,K)
      EndIf
      EndDo
C
      ECorr=ECorr+Four*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      EndIf
C
      If(IP.Eq.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NInte1
c herer!!!
      If(Abs(EExcit(K)).Gt.1.D-8) Then 
c      If(EExcit(K).Lt.Zero) Then 
      SumY=SumY+TrGamm(IPR,K)*TrGamm(IQS,K)
      EndIf
      EndDo
C
      ECorr=ECorr+Two*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))
C
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
C     COMPUTE DELTA
C
      Delta=Zero
C
      Do IP=NCore+1,NBasis
      Do IQ=NCore+1,NBasis
      Do IR=NCore+1,NBasis
      Do IS=NCore+1,NBasis
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))) Then
C
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      Aux=Two*GammaAlph(IPR)*GammaAlph(IQS)
      IQR=(Max(IR,IQ)*(Max(IR,IQ)-1))/2+Min(IR,IQ)
      If(IP.Eq.IS) Aux=Aux-GammaAlph(IQR)
      If(IP.Eq.IS.And.IQ.Eq.IR) Aux=Aux+Occ(IQ)
      If(IP.Eq.IR.And.IQ.Eq.IS) Aux=Aux-Two*Occ(IP)*Occ(IQ)
C
      Aux=Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      Delta=Delta+Aux
C
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do IP=NCore+1,NBasis
      Do IQ=NCore+1,NBasis
      Do IR=NCore+1,NBasis
      If((IGem(IP).Eq.IGem(IQ)).And.(IGem(IP).Ne.IGem(IR))) Then
C
      IPQ=(Max(IQ,IP)*(Max(IQ,IP)-1))/2+Min(IQ,IP)
      Aux=-GammaAlph(IPQ)
      If(IP.Eq.IQ) Aux=Aux+Occ(IP)
      Aux=Two*Occ(IR)*Aux*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IR))-TwoNO(NAddr3(IP,IR,IQ,IR)))
      Delta=Delta+Aux
C
      EndIf
      EndDo
      EndDo
      EndDo
C
      Do IP=NCore+1,NBasis
      Do IQ=NCore+1,NBasis
      If(IGem(IP).Ne.IGem(IQ)) Then
C
      IPQ=(Max(IQ,IP)*(Max(IQ,IP)-1))/2+Min(IQ,IP)
      Delta=Delta+Two*H1(IPQ)*GammaAlph(IPQ)
C
      EndIf
      EndDo
      EndDo
C
      Return
      End

*Deck IntegAlpha 
      Subroutine IntegAlpha(Occ,TwoAlph,H1Alph,
     $ XOne,TwoEl,NBasis,NInte1,NInte2,ACAlpha)
C
C     TRANSFORM INTEGRALS TO BE SUITABLE FOR AC
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Dimension Occ(NBasis),XOne(NInte1),
     $ TwoEl(NInte2),TwoAlph(NInte2),H1Alph(NInte1)
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      Do I=1,NInte1
      H1Alph(I)=XOne(I)
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      H1Alph(IJ)=ACAlpha*H1Alph(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      H1Alph(IJ)=H1Alph(IJ)+Aux
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
      TwoAlph(NAdd)=TwoEl(NAdd)
C
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ TwoAlph(NAdd)=ACAlpha*TwoAlph(NAdd)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Return
      End


