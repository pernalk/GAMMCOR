*Deck ACECORR
      Subroutine ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDim,NGOcc,NGem,
     $ IndN,IndX,NDimX)
C
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING AC CORRELATION ENERGY FORMULA
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),
     $ XOne(NInte1),IndAux(NBasis),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ EigVecR(NDim*NDim),
     $ Eig(NDim),EGOne(NGem),
     $ UNOAO(NBasis,NBasis),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension XGrid(100), WGrid(100)
C
      If(IFlSnd.Eq.1) Then
C
      If(ICASSCF.Ne.1) Then
       Call SndOrder(ECorr1,ECorr01,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,Eig,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux)
      Write
     $ (6,'(/,2X,''EGVB+ENuc, 0th+1st-order ECorr, AC0-GVB '',
     $ 4X,3F15.8,/)') ETot+ENuc,ECorr1,ETot+ENuc+ECorr1
      EndIf
C
      If(ICASSCF.Eq.1) Then
C 

      Call AC0CAS(ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
      Write
     $ (6,'(/,X,''ECASSCF+ENuc, AC0-Corr, AC0-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      EndIf
C
      If(IFlFrag1.Eq.1) Then
C
      Write(6,'(/,2X,''*** Embedding-AC0-GVB Calculation ***'',/)')
C
      NFrag=NGem-1
      Call FragEcorr(ETot,ENuc,ECorr2,EGOne,EigVecR,Eig,ABPLUS,ABMIN,
     $  UNOAO,Occ,TwoNO,URe,XOne,IndAux,NBasis,NInte1,
     $           NInte2,NDim,NGem,NGOcc,IFl12,NFrag)
C
      Write(6,'(2/,2X,''*** SUMMARY ***'')')
C
       Write
     $ (6,'(/,2X,''EGVB+ENuc, 0th+1st-order ECorr, AC0-GVB '',
     $ 4X,3F15.8,/)') ETot+ENuc,ECorr1,ETot+ENuc+ECorr1
C
      Write
     $ (6,'(/,2X,''EGVB+ENuc, 0th+1st-order ECorr, Emb-AC0-GVB '',
     $ 4X,3F15.8,/)') ETot+ENuc,ECorr2,ETot+ENuc+ECorr2
C
C     If(IFlFrag1.Eq.0)
      EndIf
C
      Return
C
C     If(IFlSnd.Eq.1)
      EndIf
C
C     GENERATE ABSCISSAS AND WEIGHTS FOR GAUSSIAN-LEGENDRE QUADRATURE
C
      NGrid=5
C
      Call GauLeg(Zero,One,XGrid,WGrid,NGrid)
C 
      ECorr=Zero
      Do I=1,NGrid
C
      ACAlpha=XGrid(I)
C  
      Call ACEInteg(ECorrA,TwoNO,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
      Write(*,*)'ACAlpha ',ACAlpha,' W_ALPHA ',ECorrA
C
      ECorr=ECorr+WGrid(I)*ECorrA
C
      EndDo
C
      If(ICASSCF.Eq.1) Then
C
      ETot=EGOne(1)
      Write
     $ (6,'(/,2X,''ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      Else 
C
      Write
     $ (6,'(/,2X,''EGVB+ENuc, Corr, AC-ERPA-GVB '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      EndIf
C
      stop 
C
C     SCANNING THROUGH ACAlpha
C
      ACAlpha=1.0
C
      Points=21
      Do I=1,Points
C
      ACAlpha=(I-1)/(Points-1.D0)
C
      Call ACEInteg(ECorr,TwoNO,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
      If(ACAlpha.Eq.Zero) Then
      ECorr0=ECorr
      Write(6,'(1X,'' ECorr0 '',4X,F15.8)') ECorr0
      Write(6,'(2X,''Alpha, W_ALPHA '',4X,2F15.8)')ACAlpha,ECorr
C
      Else       
C
      Write(6,'(2X,''Alpha, W_ALPHA '',4X,2F15.8)')ACAlpha,ECorr
C
      EndIf
C
      EndDo
C    
C     SCAN CLOSE TO 1.D0  
C
      ACAlpha0=0.985
C
      Points=30
      Do I=1,Points
C
      ACAlpha=ACAlpha0+I*(1.D0-ACAlpha0)/Points
C
      Call ACEInteg(ECorr,TwoNO,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
      Write(6,'(2X,''Alpha, W_ALPHA '',4X,2F15.8)')ACAlpha,ECorr
C
      EndDo
C
      Return
      End

*Deck AC0CAS
      Subroutine AC0CAS(ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigY,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),XOne(NInte1),Occ(NBasis),TwoNO(NInte2),
     $ IndAux(NBasis),
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ Eig(NDimX),EigY(NDimX*NDimX),IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
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
      If(NDimB.Ne.0)
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)      
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS
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
      If(NDimB.Ne.0)
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
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
      If(NDimB.Ne.0) 
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
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
C     DONE 0TH-ORDER CALCULATIONS  
C
      Write(6,'(/,
     $" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
C
      Deallocate(RDM2Act)
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
      Do NU=1,NoEig
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
      EndDo
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
c herer ???
c works best for tme
c     if(.not. ( (ip.eq.nele.and.ir.le.NInAcCAS).or.
c     $ (iq.eq.nele.and.is.le.NInAcCAS) ) ) then
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
c herer ???
c      endif
c
C     endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
      EndIf
C
      EndDo
      EndDo
C
      ECorr=EAll-EIntra
C
      Return
      End

*Deck Y01CAS
      Subroutine Y01CAS(TwoNO,Occ,URe,XOne,ABPLUS,ABMIN,
     $ EigY,EigY1,Eig,Eig1,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,IFlag0)
C
C     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
C     IN THE 1ST-ORDER APPROXIMATION
C
C     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
C              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),XOne(NInte1),Occ(NBasis),TwoNO(NInte2),
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ Eig(NDimX),Eig1(NDimX),
     $ EigY(NDimX*NDimX),EigY1(NDimX*NDimX),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Dimension C(NBasis),HNO(NInte1),
C     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ EigX(NDimX*NDimX),
     $ IEigAddY(2,NDimX),IEigAddInd(2,NDimX),IndBlock(2,NDimX),
c NEW 11/07/2018
     $ IMatch(NDimX)
c     NEW 25/07/2018
      Real*8 AuxCoeff(3,3,3,3)
      
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
C     PREPARE AuxCoeff
      Do L=1,3
      Do K=1,3
      Do J=1,3
      Do i=1,3
C            
      If((I==J).And.(J==K).And.(K==L)) Then
      AuxCoeff(I,J,K,L) = 1
      Else
      AuxCoeff(I,J,K,L) = 0
      EndIf
C      
      EndDo
      EndDo
      EndDo
      EndDo
C     READ 2RDM
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
C
Cmh      Write(*,*) 'HNO-Ka',norm2(HNO)
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
C     DO NOT USE IGFact!!
c$$$  C
c$$$      IGFact(NAdd)=1
c$$$      If(.Not.(
c$$$     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
c$$$     $ IGFact(NAdd)=0
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
      If(AuxCoeff(IGem(IT),IGem(IT),IGem(IP),IGem(IQ)).Eq.1) Then   
C      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.1) Then
       AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndIf
      EndDo
      EndDo
      EndDo
C
C      Write(*,*) 'AuxI-Ka',norm2(AuxI)
C      Write(*,*) 'AuxIO-Ka',norm2(AuxIO)
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
      If(AuxCoeff(IGem(IT),IGem(IW),IGem(IP),IGem(IU)).Eq.1)
C      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.1)
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
C      Write(*,*) 'WMAT-Ka',norm2(WMAT)
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
C      Write(*,*) 'NDimB-Act',NDimB
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
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,AuxCoeff,    
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
      EndDo
      EndDo
C
      EndIf
      EndDo
      EndDo
C
C      If(NDimB.Ne.0)
C     $ Write(*,*) 'AB0-Ka',norm2(ABPLUS(1:NDimB**2)),
C     $     norm2(ABMIN(1:NDimB**2))
C
C      write(*,*) 'AB+',
C     $     norm2(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])
C     $     -transpose(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])))
C      write(*,*) 'AB-',norm2(ABMIN(1:nAA,1:nAA)-transpose(ABMIN(1:nAA,1:nAA)))
C      
      If(NDimB.Ne.0)
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)      
C
Cmh      Write(*,*) 'EigY,X,val',norm2(EigY(NFree2:NDimB**2)),
Cmh     $ norm2(EigX(NFree2:NDimB**2))
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
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,AuxCoeff,
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
C      If(NDimB.Ne.0)
C     $     Write(*,*) 'ABia-Ka',IQ,norm2(ABPLUS(1:NDimB**2)),
C     $     norm2(ABMIN(1:NDimB**2))
C     
C      write(*,*) 'AB+',
C     $     norm2(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])
C     $     -transpose(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])))
      
C
      If(NDimB.Ne.0)
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
C
C      Write(*,*) 'EigY,X,val',norm2(EigY(NFree2:NFree2+NDimB**2)),
C     $ norm2(EigX(NFree2:NFree2+NDimB**2)),
C     $ norm2(Eig(NFree1:NFree1+NDimB))
C     
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     Do IQ
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
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,AuxCoeff,         
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
c      If(NDimB.Ne.0)
C     $     Write(*,*) 'ABav-Ka',IP,norm2(ABPLUS(1:NDimB**2)),
C     $     norm2(ABMIN(1:NDimB**2))
CC     
C      write(*,*) 'AB+',
C     $     norm2(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])
C     $     -transpose(reshape(ABPLUS(1:NDimB**2),shape=[NDimB,NDimB])))
C 
C     
      If(NDimB.Ne.0)
     $Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)

C      Write(*,*) 'EigY,X,val',norm2(EigY(NFree2:NFree2+NDimB**2)),
C     $ norm2(EigX(NFree2:NFree2+NDimB**2)),
C     $ norm2(Eig(NFree1:NFree1+NDimB))

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
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IP,IQ,Occ,HNO,AuxCoeff, 
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $     NInte1,NInte2,NBasis)
C
C      Write(*,*) 'ABiv-Ka',IP,IQ,ABP
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
      Write(6,'(" *** DONE WITH 0TH-ORDER IN Y01CAS ***")')
C
C     DONE 0TH-ORDER CALCULATIONS  
C
C NEW 11/07/2018
C     Check if NoEig=NDimX - they should be equal!
      If(NoEig.Ne.NDimX) Stop 'Fatal error in Y01CAS: NoEig.Ne.NDimX!'
C
      Do I=1,NDimX
      IP=IndN(1,I)
      IQ=IndN(2,I)
      Do J=1,NDimX
      If(IP.Eq.IndBlock(1,J).And.IQ.Eq.IndBlock(2,J))
     $ IMatch(I)=J
      EndDo
      EndDo
C
      If(IFlag0.Eq.1) Then
C     
      Do I=1,NFree2
      EigY1(I)=EigY(I)
      EndDo
      Do I=1,NoEig
      Do J=1,NoEig
      ABPLUS(I+(J-1)*NoEig)=Zero
      If(I.Eq.J) ABPLUS(I+(J-1)*NoEig)=One
      EndDo
      EndDo 
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      EigY(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      EigY(NU+(MU-1)*NoEig)=EigY(NU+(MU-1)*NoEig)
     $ +ABPLUS(NU+(I-1)*NoEig)*EigY1(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
C NEW 11/07/2017
C
C     RESORT Y0 ACCORDING TO IndN
C
      Call CpyM(ABPLUS,EigY,NDimX)
C
      Do MU=1,NDimX
      Do I=1,NDimX
      EigY((MU-1)*NoEig+I)=ABPLUS((MU-1)*NoEig+IMatch(I))
      EndDo
      EndDo
C
C     RETURN IF ONLY 0TH-ORDER Y REQUESTED
C       write(*,*) 'EigY(r)-Ka',norm2(EigY)     
C       write(*,*) 'Eig(r)-Ka ',norm2(Eig)
C           
      Return
      EndIf
C
      Write(6,'(" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $     RDM2Act,NRDM2Act,
     $     AuxCoeff, 
     $     C,Ind1,Ind2,
     $     IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
C
      Deallocate(RDM2Act)
C
      Write(6,'(" *** DONE WITH COMPUTING AB(1) MATRICES ***")')
C
Cmh      Write(*,*) 'EigY-Ka1',norm2(EigY)
Cmh      Write(*,*) 'EigX-Ka1',norm2(EigX)
C
C     1ST-ORDER PART
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      EigY1(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      EigY1(NU+(MU-1)*NoEig)=EigY1(NU+(MU-1)*NoEig)
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
     $ +EigX(IStart+II)*EigY1(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo

C
Cmh       Write(*,*) 'YABpY-Ka1',norm2(ABPLUS)
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      EigY1(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      EigY1(NU+(MU-1)*NoEig)=EigY1(NU+(MU-1)*NoEig)
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
     $ +EigY(IStart+II)*EigY1(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo
Cmh       Write(*,*) 'YABmY-Ka1',norm2(ABMIN)
C     
      EigY1(1:NoEig*NoEig)=Zero
      Do MU=1,NoEig
C
      If(Eig(MU).Ne.Zero) Then
C
      Do NU=1,NoEig
C
      If(Eig(NU).Ne.Zero) Then
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
C
      Aux1=(ABPLUS(MU+(NU-1)*NoEig)-ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)+Eig(NU))
      Aux2=Zero
c herer!!!
C      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-5)) Aux2=
      If(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12) Aux2=
Cmh      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12)) Aux2=
     $  (ABPLUS(MU+(NU-1)*NoEig)+ABMIN(MU+(NU-1)*NoEig))/
     $     (Eig(MU)-Eig(NU))
C
      EigY1(I+(MU-1)*NoEig)=EigY1(I+(MU-1)*NoEig)+
     $(Aux1+Aux2)*EigY(IStart+II)
C
      II=II+1
      EndDo
C
      EndIf
      EndDo
      EndIf
      EndDo
C
C     1ST-ORDER EIG
C
      Do NU=1,NoEig
      Eig1(NU)=ABPLUS(NU+(NU-1)*NoEig)+ABMIN(NU+(NU-1)*NoEig)
      EndDo
C
      Do I=1,NFree2
      EigX(I)=EigY(I)
      EndDo
      Do I=1,NoEig
      Do J=1,NoEig
      ABPLUS(I+(J-1)*NoEig)=Zero
      If(I.Eq.J) ABPLUS(I+(J-1)*NoEig)=One
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      EigY(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      EigY(NU+(MU-1)*NoEig)=EigY(NU+(MU-1)*NoEig)
     $ +ABPLUS(NU+(I-1)*NoEig)*EigX(IStart+II)
      II=II+1
      EndDo
C
C     Y^(1) ALREADY IN EigY1
C
      EndDo
      EndDo
C NEW 11/07/2018
C     
C     SORT Y0 AND Y1 ACCORDING TO IndN
C
      Call CpyM(ABPLUS,EigY,NDimX) 
      Call CpyM(ABMIN,EigY1,NDimX)
C
      Do MU=1,NDimX
C
      Do I=1,NDimX
      EigY((MU-1)*NoEig+I)=ABPLUS((MU-1)*NoEig+IMatch(I))
      EigY1((MU-1)*NoEig+I)=ABMIN((MU-1)*NoEig+IMatch(I)) 
      EndDo
C
      EndDo
C
      Write(*,*) 'EigY1-Ka1',norm2(EigY1)     
Cmh      Write(*,*) 'EigY-Ka ',norm2(EigY)
C     
      Return
      End



*Deck AB1_CAS
      Subroutine AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $     RDM2Act,NRDM2Act,
     $     AuxCoeff,
     $     C,Ind1,Ind2,
     $     IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
C
C     COMPUTE THE ALPHA-DEPENDENT PARTS OF A+B AND A-B MATRICES
C
C     RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
C     THE FOLLOWING SYMMETRY IS ASSUMED
C     RDM2(ij,kl) = RDM2(kl,ij)
C     ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
C     SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NoEig*NoEig),ABMIN(NoEig*NoEig),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ RDM2Act(NRDM2Act),
     $ AuxCoeff(3,3,3,3),
C     $ IGFact(NInte2),
     $ C(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ IndBlock(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),WMAT(NBasis,NBasis),
     $  AuxI(NInte1),AuxIO(NInte1)
C
      Do I=1,NoEig*NoEig
      ABPLUS(I)=Zero
      ABMIN(I)=Zero
      EndDo
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
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Eq.IGem(J)) Then
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      HNO(IJ)=-Aux
C
      EndIf
C
      EndDo
      EndDo
C
C      Write(*,*) 'HNO-Ka1',norm2(HNO)
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
C     
C     AUXILIARY MATRIX AuxI  
C   
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      AuxIO(IPQ)=Zero
      Do IT=1,NOccup
      If(AuxCoeff(IGem(IT),IGem(IT),IGem(IP),IGem(IQ)).Eq.0) Then
C     If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) Then
      AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndIf
      EndDo
      EndDo
      EndDo
C
C      Write(*,*) 'AuxI-Ka1',norm2(AuxI)
C      Write(*,*) 'AuxIO-Ka1',norm2(AuxIO)
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
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
C
      If(AuxCoeff(IGem(IT),IGem(IW),IGem(IP),IGem(IU)).Eq.0) Then
C     If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) Then
C
      Do IR=1,NOccup
      WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))
     $ *FRDM2(IW,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IT,IU,IP,IW))
     $ *FRDM2(IW,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
      EndDo
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C      Write(*,*) 'WMAT-Ka1',norm2(WMAT)
C
      icount=0
      Call CPU_TIME(START_TIME)
      Do IRow=1,NoEig
C
      IR=IndBlock(1,IRow)
      IS=IndBlock(2,IRow)
C
C      TURN OFF SYMMETRIZE !
Cmh      Do ICol=1,IRow
      Do ICol=1,NoEig     
C
      IPP=IndBlock(1,ICol)
      IQQ=IndBlock(2,ICol)
C
Cmh      If( .NOT.(IGem(IR).Eq.IGem(IS).And.IGem(IR).Eq.IGem(IPP)
Cmh     $ .And.IGem(IR).Eq.IGem(IQQ)) ) Then
CmhC
Cmh      If( (Occ(IR)*Occ(IS).Eq.Zero.And.Occ(IPP)*Occ(IQQ).Eq.Zero
Cmh     $ .And.Abs(TwoNO(NAddr3(IR,IS,IPP,IQQ))).Lt.1.D-7)
Cmh     $.Or.
Cmh     $((Occ(IR).Eq.One.Or.Occ(IS).Eq.One)
Cmh     $ .And.
Cmh     $ (Occ(IPP).Eq.One.Or.Occ(IQQ).Eq.One)
Cmh     $ .And.Abs(TwoNO(NAddr3(IR,IS,IPP,IQQ))).Lt.1.D-7)) Then
CmhC
Cmh      icount=icount+1
C
Cmh      Else  
C
      Do IP=IQQ,IPP,IPP-IQQ
      Do IQ=IQQ,IPP,IPP-IQQ
      If(IP.Ne.IQ) Then
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=Zero
C     If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=One
      If(AuxCoeff(IGem(IP),IGem(IQ),IGem(IR),IGem(IS)).Eq.0)
     $ AuxTwoPQRS=One
      AuxPQRS=AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
C
      If(AuxTwoPQRS.Eq.One)
     $ Arspq=Arspq+Occ(IP)*Occ(IR)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IS),IGem(IQ),IGem(IT),IGem(IU)).Eq.0)
C     If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0)
     $ Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*
     $  FRDM2(IP,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*
     $ FRDM2(IP,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
C
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxIO(IQS)
C
      EndIf
      EndIf
C
C     T3+T4 
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
C
      If(AuxTwoPQRS.Eq.One)
     $ Arspq=Arspq+Occ(IQ)*Occ(IS)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IU),IGem(IT),IGem(IP),IGem(IR)).Eq.0)
C      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) 
     $ Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*
     $ FRDM2(IS,IT,IQ,IU,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*
     $ FRDM2(IS,IT,IU,IQ,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
C
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxIO(IPR)
C
      EndIf
C     
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IP),IGem(IT),IGem(IS),IGem(IU)).Eq.0)
C      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0)
     $ Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*
     $ FRDM2(IT,IU,IQ,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T6      
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IT),IGem(IQ),IGem(IU),IGem(IR)).Eq.0)
C      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0)
     $ Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*
     $ FRDM2(IS,IP,IU,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      If(IR.Gt.IS.And.IP.Gt.IQ) Then
      ABPLUS(IRow+(ICol-1)*NoEig)=ABPLUS(IRow+(ICol-1)*NoEig)+Arspq 
      ABMIN(IRow+(ICol-1)*NoEig)=ABMIN(IRow+(ICol-1)*NoEig)+Arspq
      EndIf
C    
      If(IR.Gt.IS.And.IQ.Gt.IP) Then
      ABPLUS(IRow+(ICol-1)*NoEig)=ABPLUS(IRow+(ICol-1)*NoEig)+Arspq
      ABMIN(IRow+(ICol-1)*NoEig)=ABMIN(IRow+(ICol-1)*NoEig)-Arspq
      EndIf
C
c     If(IP.Ne.IQ)
      EndIf
C
C     end of IP,IQ LOOPS
      EndDo
      EndDo
C
 1000 Continue
Cmh      EndIf
Cmhc     If IGem ....
Cmh      EndIf
C
      EndDo
      EndDo
      Call CPU_TIME(END_TIME)
      Write(6,'(X,"TIME SPENT ON CONSTRUCTING AB(1)"
     $ ,F10.2)')END_TIME-START_TIME
      write(*,*)'icount',icount
C
C     DIVIDE BY C'c AND SYMMETRIZE
C
      Do I=1,NoEig
      IP=IndBlock(1,I)
      IQ=IndBlock(2,I)
C     TURN OFF SYM! for testing
      Do J=1,NoEig
Cmh      Do J=1,I
      IR=IndBlock(1,J)
      IS=IndBlock(2,J)
C
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ ABPLUS(I+(J-1)*NoEig)=ABPLUS(I+(J-1)*NoEig)
     $/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If((C(IP)-C(IQ))*(C(IR)-C(IS)).Ne.Zero)
     $ ABMIN(I+(J-1)*NoEig)=ABMIN(I+(J-1)*NoEig)
     $/(C(IP)-C(IQ))/(C(IR)-C(IS))
C
C    TURN OFF SYM! for testing
Cmh      ABPLUS(J+(I-1)*NoEig)=ABPLUS(I+(J-1)*NoEig)
Cmh      ABMIN(J+(I-1)*NoEig)=ABMIN(I+(J-1)*NoEig)
C
      EndDo
      EndDo
C
      Write(*,*) 'ABPlus-Ka1',norm2(ABPLUS),'ABMin-Ka1',norm2(ABMIN)
C
      Return
      End

*Deck AB0ELEMENT
      Subroutine AB0ELEMENT(ABPL,ABMIN,IR,IS,IPP,IQQ,Occ,HNO,AuxCoeff, 
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,NInte1,
     $ NInte2,NBasis)
C
C     FOR A GIVEN SET OF INDICES IR,IS,IPP,IQQ RETURNS 
C     VALUES OF ABPLUS AND ABMIN MATRICES FOR ALPHA=0
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ Occ(NBasis),HNO(NInte1),
C    $ IGFact(NInte2),
     $ AuxCoeff(3,3,3,3), 
     $ TwoNO(NInte2),AuxI(NInte1),AuxIO(NInte1),
     $ WMAT(NBasis,NBasis),RDM2Act(NRDM2Act),
     $ C(NBasis),Ind1(NBasis),Ind2(NBasis)
C
      ABPL=Zero
      ABMIN=Zero
C
      Do IP=IQQ,IPP,IPP-IQQ
      Do IQ=IQQ,IPP,IPP-IQQ
      If(IP.Ne.IQ) Then
C
      If(IP.Gt.IQ) IPQ=(IP**2-3*IP+2)/2+IQ
      If(IQ.Gt.IP) IQP=(IQ**2-3*IQ+2)/2+IP
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=AuxCoeff(IGem(IP),IGem(IQ),IGem(IR),IGem(IS))
c     AuxTwoPQRS=One
c      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=Zero
C
      AuxPQRS=Zero
      If(AuxTwoPQRS.Eq.One)
     $ AuxPQRS=Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
C
      If(AuxTwoPQRS.Eq.One) Arspq=Arspq+Occ(IP)*Occ(IR)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
C
      Else
C
      If(IGem(IS).Eq.IGem(IQ)) Then
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IS),IGem(IQ),IGem(IT),IGem(IU)).Eq.1)
C     If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.1) 
     $ Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*
     $  FRDM2(IP,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*
     $  FRDM2(IP,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndIf
C
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxIO(IQS)
C
      EndIf
      EndIf
C
C     T3+T4 
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
C
      If(AuxTwoPQRS.Eq.One) Arspq=Arspq+Occ(IQ)*Occ(IS)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
C
      Else
C
      If(IGem(IP).Eq.IGem(IR)) Then
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IU),IGem(IT),IGem(IP),IGem(IR)).Eq.1)     
C      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.1)
     $ Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*
     $ FRDM2(IS,IT,IQ,IU,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*
     $ FRDM2(IS,IT,IU,IQ,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndIf
C
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxIO(IPR)
C
      EndIf
C     
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
C
      Else
C
      If(IGem(IP).Eq.IGem(IS)) Then
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IP),IGem(IT),IGem(IS),IGem(IU)).Eq.1)
C     If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.1) 
     $ Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*
     $ FRDM2(IT,IU,IQ,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndIf
C
      EndIf
      EndIf
C
C     T6      
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
C
      Else
C
      If(IGem(IR).Eq.IGem(IQ)) Then
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      If(AuxCoeff(IGem(IT),IGem(IQ),IGem(IU),IGem(IR)).Eq.1)
C      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.1) 
     $ Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*
     $ FRDM2(IS,IP,IU,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndIf
C
      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      ABPL=ABPL+Arspq
      If(IP.Gt.IQ) Then
      ABMIN=ABMIN+Arspq
      Else
      ABMIN=ABMIN-Arspq
      EndIf
C    
c     If(IP.Ne.IQ)
      EndIf
C     end of IP,IQ LOOPS
      EndDo
      EndDo
C
      IP=IPP
      IQ=IQQ
C
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ ABPL=ABPL/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If((C(IP)-C(IQ))*(C(IR)-C(IS)).Ne.Zero)
     $ ABMIN=ABMIN/(C(IP)-C(IQ))/(C(IR)-C(IS))
C
      Return
      End

*Deck SndOrder
      Subroutine SndOrder(ECorr,ECorr0,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,Eig,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension
     $ URe(NBasis,NBasis),XOne(NInte1),Occ(NBasis),TwoNO(NInte2),
     $ IndAux(NBasis),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ Eig(NDim),EGOne(NGem)
C
C     LOCAL ARRAYS
C
      Dimension IndX(NDim),IndN(2,NDim),C(NBasis),
     $ EigVY2(NBasis*(NBasis-1)),IndP(NBasis,NBasis),
     $ AMAT(NDim*NDim),BMAT(NDim*NDim)
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C
C     0-TH ORDER EIGENVECTORS 
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
C
      If(IGem(IP).Eq.IGem(IQ)) Then
C
      EigVY2(IPQ)=Half
C  
      Else
C
      If(Occ(IP)-Occ(IQ).Ne.Zero) Then 
      If(Occ(IQ).Gt.Occ(IP))EigVY2(IPQ)=Half*(C(IQ)-C(IP))/(C(IQ)+C(IP))
      If(Occ(IP).Gt.Occ(IQ))EigVY2(IPQ)=Half*(C(IP)-C(IQ))/(C(IQ)+C(IP))
      Else
      EigVY2(IPQ)=Zero 
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
C     FIND THE 0TH-ORDER EPSILONS
C
      ACAlpha=Zero
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
      Do I=1,NDim*NDim
      AMAT(I)=ABPLUS(I)
      BMAT(I)=ABMIN(I)
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
      IRSIRS=(IRS-1)*NDim+IRS
      Eig(IRS)=Zero
C
      If(Occ(IR).Ne.Occ(IS).And.IGem(IR).Ne.IGem(IS)) Then
      AA=Half*((C(IR)+C(IS))**2*ABPLUS(IRSIRS)
     $       + (C(IR)-C(IS))**2*ABMIN(IRSIRS))
      Eig(IRS)=ABS(AA/(Occ(IR)-Occ(IS)))
      Else
      Eig(IRS)=ABS(ABPLUS(IRSIRS))
      EndIf
C
      EndDo
      EndDo
C
      ACAlpha=One
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
C     AMAT AND BMAT WILL INCLUDE 1ST-ORDER A+ AND A- MATRICES, RESPECTIVELY
C
      Do I=1,NDim*NDim
      AMAT(I)=ABPLUS(I)-AMAT(I)
      BMAT(I)=ABMIN(I)-BMAT(I)
      EndDo 
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
C
      IJ=IJ+1
      IndP(I,J)=0
      IndP(J,I)=0
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate active degenerate orbitals if from different geminals
      If((IGem(I).Ne.IGem(J)).And.(IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-2) ) Then
      Write(*,*)"Discarding nearly degenerate pair",I,J
      Else
C    
C     If IFlCore=0 do not include core (inactive) orbitals  
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
C
      IndP(I,J)=1
      IndP(J,I)=1
      EndIf
C
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
      ECorr=Zero
      ECorr0=Zero
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
      IPQ=IndX(I)
C
      If(Occ(IP).Gt.Occ(IQ)) Then
      Aux=Occ(IQ)*(Occ(IP)-One)
      Else
      Aux=Occ(IP)*(Occ(IQ)-One)
      EndIf
C
C     0th-order correlation only if (IP,IQ) pair is allowed
      If(IGem(IP).Ne.IGem(IQ)) ECorr0=ECorr0+
     $ Two*Aux*TwoNO(NAddr3(IP,IQ,IP,IQ))*IndP(IP,IQ)
C
      Do J=1,NDimX
C
      IR=IndN(1,J)
      IS=IndN(2,J)
      IRS=IndX(J)
C
      If(IP.Gt.IQ.And.IR.Gt.IS) Then
C
      EPSJI=Eig(IRS)+Eig(IPQ)
      If(EPSJI.Ne.Zero) Then
      EPSJI=One/EPSJI 
      Else
      EPSJI=Zero
      EndIf
C
      IJ=(IRS-1)*NDim+IPQ 
      Aux=(Half*AMAT(IJ)-Two*EigVY2(IPQ)*EigVY2(IRS)*BMAT(IJ))
     $    *EPSJI 
C
C     Save Aux - it may be needed in embedding calculations
      ABPLUS((J-1)*NDimX+I)=(C(IP)+C(IQ))*(C(IR)+C(IS))*Aux
C
C     only if indices not from the same geminals and (IP,IQ) and 
C     (IR,IS) pairs are accepted
      If(.Not.(IGem(IP).Eq.IGem(IQ).And.IGem(IR).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IR))) Then
      ECorr=ECorr+(C(IP)+C(IQ))*(C(IR)+C(IS))*Aux
     $     *TwoNO(NAddr3(IP,IQ,IR,IS))*IndP(IP,IQ)*IndP(IR,IS)
      EndIf
C
      EndIf
      EndDo
      EndDo
C 
      ECorr=ECorr+ECorr0
C
      Write
     $ (6,'(/,2X,''0-ALPHA-ORDER CORRELATION '',2X,F15.8)') ECorr0
C
      Return
      End

*Deck ACEInteg
      Subroutine ACEInteg(ECorr,TwoNO,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension
     $ URe(NBasis,NBasis),Occ(NBasis),TwoNO(NInte2),
     $ XOne(NInte1),UNOAO(NBasis,NBasis),
     $ IndAux(NBasis),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ EigVecR(NDim*NDim),
     $ Eig(NDim),EGOne(NGem),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension IPair(NBasis,NBasis)
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II) 
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
C     CALCULATE THE A+B AND A-B MATRICES
C   
      If(ICASSCF.Eq.0) Then
C
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
      ElseIf(ICASSCF.Eq.1) Then
C
c      Call RDMFT_AB(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
c     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
c      Call Gamma2_AB(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
c     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
c
      EGOne(1)=ECASSCF
C
      EndIf
C
      If(IFlFrag1.Eq.1) Then
C
      ETot=Zero
      ENuc=Zero
C
      NFrag=NGem-1
C
c      Call OneTwoBody(ETot,ENuc,ECorr,EGOne,EigVecR,Eig,ABPLUS,ABMIN,
c     $ Occ,TwoNO,IndAux,NBasis,NInte1,NInte2,NDim,NGem,NGOcc,IFl12)
C
      Call FragEcorr(ETot,ENuc,ECorr,EGOne,EigVecR,Eig,ABPLUS,ABMIN,
     $  UNOAO,Occ,TwoNO,URe,XOne,IndAux,NBasis,NInte1,
     $           NInte2,NDim,NGem,NGOcc,IFl12,NFrag)
C
C     multiply by factor 2 (in OneTwoBody everything has been divided by '2')
      ECorr=Two*ECorr      
C
      Return
C
      EndIf
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
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      Return
      End

*Deck ACEneERPA
      Subroutine ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Parameter(SmallE=1.D-3,BigE=1.D8)
C
C     ONLY EXCITATIONS SMALLER THAN BigE AND GREATER THAN SmallE ARE INCLUDED 
C     
      Include 'commons.inc'
C     
      Dimension EigVecR(NDimX*NDimX),
     $ Eig(NDimX),URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),Skipped(NDimX)
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C
      ECorr=Zero
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IR=IndN(2,I)
C
      Do J=1,NDimX
C
      IQ=IndN(1,J)
      IS=IndN(2,J)
C
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ)).And.IP.Gt.IR.And.IQ.Gt.IS) Then
C
      ISkippedEig=0
      SumY=Zero
      Do K=1,NDimX
      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE) Then
      SumY=SumY+EigVecR((K-1)*NDimX+I)*EigVecR((K-1)*NDimX+J)
      Else
      ISkippedEig=ISkippedEig+1
      Skipped(ISkippedEig)=Eig(K)
      EndIf
      EndDo
C
      Aux=Two*(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      If(IR.Eq.IS.And.IP.Eq.IQ) 
     $ Aux=Aux
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
C     endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
      EndIf
C
      EndDo
      EndDo
C
      If(ISkippedEig.Ne.0) Then
      Write(6,'(/,1X,"The number of discarded eigenvalues is",I4)')
     $  ISkippedEig
      Do II=1,ISkippedEig
      Write(6,*)'Skipped',II,Skipped(II)
      EndDo
      EndIf
C
      Return
      End

*Deck TEST_COMM_HAPKA
      Subroutine COMMTST(NBasis)
C     TESTS COMMONS PASSED BY SAPT
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C     
      Include 'commons.inc'
C
      write(*,*) CICoef(1:NBasis)
      write(*,*) IGem(1:NBasis)
      End

*Deck ACABMAT0
      Subroutine ACABMAT0(AMAT,BMAT,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,IFlag)
C
C     IFlag = 1 - AMAT AND BMAT WILL CONTAIN (A+B)/C+/C+ AND (A-B)/C-/C-, RESPECTIVELY
C             0 - AMAT AND BMAT WILL CONTAIN A ANB B MATRICES, RESPECTIVELY
C
C     ACAlpha - Alpha-connection parameter in AC
C     HNO AND TwoMO are modified to correspond to an alpha-Hamiltonian
C
C     STRAIGHTFORWARD IMPLEMENTATION OF THE AMAT AND BMAT DEFINITIONS (NO SPECIAL CASES CONSIDERED) 
C     
C     COMPUTE THE A+B AND A-B MATRICES IN ERPA WITH APSG APPROXIMATION
C
C     NESTED COMMUTATOR 
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C     
      Include 'commons.inc'
C     
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),C(NBasis),
     $ AuxXC(NBasis,NInte1),AuxH(NBasis,NInte1),IGFact(NInte2)
C
      If(ISAPT.Eq.1) Then
      Write(6,'(1x,a)') 'Computing response'
      Else
      Write(6,'(/,X,"***** COMPUTING AMAT, BMAT IN ACABMAT0 *****",/)')
      EndIf
C
C     TEST URe
C      Do I=1,NBasis
C      Do J=1,NBasis
C      write(6,*) URe(I,I)
C      EndDo
C      EndDo
C      write(*,*) XOne(1:NBasis)
C      Do I=1,NBasis
C      Write(6,*) XOne(I*(I-1)/2+I)
C      EndDo
C
      Write(*,*) 'NGem',NGem
C
      Do I=1,NDim
      Do J=1,NDim
      AMAT(I,J)=Zero
      BMAT(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
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
C      
C     CONSTRUCT ONE-ELECTRON PART OF AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C

      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IT,IT,I,J))-TwoMO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     TESTUMY
      Write(*,*) 'HNO-Ka',norm2(HNO)
C      
C     CONSTRUCT TWO-ELECTRON PART OF AC ALPHA-HAMILTONIAN      
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
C     AUXILIARY MATRICES 
C
      Do I=1,NBasis
C
      ISP=0
      Do IS=1,NBasis
      Do IP=1,IS
      ISP=ISP+1
C
      AuxXC(I,ISP)=Zero
      AuxH(I,ISP)=Zero
C
      Do IT=1,NBasis
C
      AuxTwo=One
      If(IGFact(NAddr3(IT,IS,IT,IP)).Eq.0) AuxTwo=ACAlpha
C
      If(IGem(IT).Eq.IGem(I)) Then
      AuxXC(I,ISP)=AuxXC(I,ISP)+C(IT)*AuxTwo*TwoMO(NAddr3(IT,IS,IT,IP))
      Else
      AuxH(I,ISP)=AuxH(I,ISP)+Occ(IT)*AuxTwo*
     $ (Two*TwoMO(NAddr3(IT,IT,IS,IP))-TwoMO(NAddr3(IT,IS,IT,IP)))
      EndIf
C
      EndDo
C
      EndDo
      EndDo
C
      EndDo
C
C     COMPUTE THE MATRICES A+B, A-B
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR
      If(IR.Gt.IS) IRS=IRS+1
C
      IPQ=0
      Do IP=1,NBasis
C
      IPS=(Max(IS,IP)*(Max(IS,IP)-1))/2+Min(IS,IP)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Do IQ=1,NBasis
      IQR=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)
      IQS=(Max(IQ,IS)*(Max(IQ,IS)-1))/2+Min(IQ,IS)

      If(IP.Gt.IQ) IPQ=IPQ+1
      If(IQ.Gt.IP) IQP=(IQ**2-3*IQ+2)/2+IP
C
      XMAT=Zero
C
      If(IQ.Eq.IR) XMAT=XMAT+(Occ(IR)-Occ(IS))*HNO(IPS)
      If(IS.Eq.IP) XMAT=XMAT-(Occ(IR)-Occ(IS))*HNO(IQR)
C
      AuxTwo=One
      If(IGFact(NAddr3(IQ,IR,IS,IP)).Eq.0) AuxTwo=ACAlpha
C
C     INTERGEMINAL CONTRIBUTIONS
C
      IGemPR=1
      If(IGem(IP).Eq.IGem(IR)) IGemPR=0
      IGemQS=1
      If(IGem(IQ).Eq.IGem(IS)) IGemQS=0
      IGemPS=1
      If(IGem(IP).Eq.IGem(IS)) IGemPS=0
      IGemQR=1
      If(IGem(IQ).Eq.IGem(IR)) IGemQR=0
C
Cmh      XMAT=XMAT+ (-Occ(IR)*Occ(IP)*IGemPR
Cmh     $            +Occ(IS)*Occ(IP)*IGemPS
Cmh     $            +Occ(IR)*Occ(IQ)*IGemQR
Cmh     $            -Occ(IS)*Occ(IQ)*IGemQS)*AuxTwo
Cmh     $ *( Two*TwoMO(NAddr3(IQ,IP,IS,IR))-TwoMO(NAddr3(IQ,IR,IS,IP)) )
C
Cmh      If(IQ.Eq.IR) XMAT=XMAT+Occ(IQ)*AuxH(IQ,IPS)-Occ(IS)*AuxH(IS,IPS)
Cmh      If(IS.Eq.IP) XMAT=XMAT+Occ(IP)*AuxH(IP,IQR)-Occ(IR)*AuxH(IR,IQR)
C
C     INTRAGEMINAL PART
C
      IGemQR=1
      If(IGem(IQ).Ne.IGem(IR)) IGemQR=0
      IGemPS=1
      If(IGem(IP).Ne.IGem(IS)) IGemPS=0
Cmh      XMAT=XMAT+(C(IQ)*C(IR)*IGemQR+C(IP)*C(IS)*IGemPS)*AuxTwo
Cmh     $ *( TwoMO(NAddr3(IP,IS,IQ,IR))+TwoMO(NAddr3(IP,IR,IQ,IS)) )
C
Cmh      If(IS.Eq.IP) XMAT=XMAT-C(IR)*AuxXC(IR,IQR)
Cmh      If(IR.Eq.IQ) XMAT=XMAT-C(IS)*AuxXC(IS,IPS)
Cmh      If(IR.Eq.IP) XMAT=XMAT-C(IP)*AuxXC(IP,IQS)
Cmh      If(IS.Eq.IQ) XMAT=XMAT-C(IQ)*AuxXC(IQ,IPR)
C
      If (IR.Gt.IS.And.IP.Gt.IQ) BMAT(IRS,IPQ)=XMAT
      If (IR.Gt.IS.And.IQ.Gt.IP) AMAT(IRS,IQP)=XMAT
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C
      Print*, "AB-Ka",norm2(AMAT),norm2(BMAT)      
C
Cmh         
      Return
C     
      If(IFlag.Eq.0) Return
C
      IRS=0
      Do IR=1,NBasis
      IGR=IGem(IR)
C
      Do IS=1,IR-1
      IRS=IRS+1
      IGS=IGem(IS)
C
      IPQ=0
      Do IP=1,NBasis
      IGP=IGem(IP)
C
      Do IQ=1,IP-1
      IPQ=IPQ+1
      IGQ=IGem(IQ)
C
      If(Abs(Abs(C(IP))-Abs(C(IQ))).Le.Delta*Abs(C(IQ)).Or.
     $   Abs(Abs(C(IR))-Abs(C(IS))).Le.Delta*Abs(C(IR))) Then
      AMAT(IRS,IPQ)=Zero
      BMAT(IRS,IPQ)=Zero
      Else
      SaveA=AMAT(IRS,IPQ)
      SaveB=BMAT(IRS,IPQ)
C
      AMAT(IRS,IPQ)=(SaveA+SaveB)/(C(IR)+C(IS))/(C(IP)+C(IQ))
      BMAT(IRS,IPQ)=(SaveA-SaveB)/(C(IR)-C(IS))/(C(IP)-C(IQ))
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     TEST FOR SAPT
C      Do I=1,NDim
CC      Write(6,*) AMAT(1,I)
C      Write(6,*) AMAT(1,I),BMAT(1,I)
C      EndDo
C
      Return
      End

*Deck ACABMAT
      Subroutine ACABMAT(AMAT,BMAT,CMAT,EMAT,EMATM,
     $ DMAT,DMATM,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ISERPA,ACAlpha)
C
C     ACAlpha - Alpha-connection parameter in AC
C     HNO AND TwoMO are modified to correspond to an alpha-Hamiltonian
C     
C     COMPUTE THE A+B AND A-B MATRICES IN ERPA WITH APSG APPROXIMATION
C
C     NESTED COMMUTATOR 
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C     
      Include 'commons.inc'
C     
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ URe(NBasis,NBasis),CMAT(NDim,NDim),EMAT(NBasis,NBasis),
     $ EMATM(NBasis,NBasis), 
     $ DMAT(NDim,NBasis),DMATM(NDim,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),C(NBasis),
     $ AuxXC(NBasis,NInte1),AuxH(NBasis,NInte1),
     $ XMu(NBasis),ANDMAT(NDim,NBasis),
     $ ADNMAT(NBasis,NDim),ADDMAT(NBasis,NBasis),
     $ Aux1(NBasis,NDim),Aux2(NBasis,NBasis)
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
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
C     CONSTRUCT ONE-ELECTRON PART OF AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I)) 
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IT,IT,I,J))-TwoMO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF AC ALPHA-HAMILTONIAN      
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
     $ TwoMO(NAdd)=ACAlpha*TwoMO(NAdd)
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRICES 
C
      Do I=1,NBasis
C
      ISP=0
      Do IS=1,NBasis
      Do IP=1,IS
      ISP=ISP+1
C
      AuxXC(I,ISP)=Zero
      AuxH(I,ISP)=Zero
C
      Do IT=1,NBasis
C
      If(IGem(IT).Eq.IGem(I)) Then
      AuxXC(I,ISP)=AuxXC(I,ISP)+C(IT)*TwoMO(NAddr3(IT,IS,IT,IP))
      Else
      AuxH(I,ISP)=AuxH(I,ISP)+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IT,IT,IS,IP))-TwoMO(NAddr3(IT,IS,IT,IP)))
      EndIf
C
      EndDo
C
      EndDo
      EndDo
C
      EndDo
C
C     XMu FOR EACH GEMINAL
C
      Do I=1,NBasis
      XMu(I)=Zero
C
      Do IP=1,NBasis
C
      If(IGem(IP).Eq.IGem(I)) Then
C
      IPP=(IP*(IP+1))/2
      XMu(I)=XMu(I)+Two*Occ(IP)*HNO(IPP)
C
      Do IQ=1,NBasis
C  
      If(IGem(IQ).Eq.IGem(I)) Then
      XMu(I)=XMu(I)+C(IP)*C(IQ)*TwoMO(NAddr3(IP,IQ,IP,IQ))
      Else
      XMu(I)=XMu(I)+Two*Occ(IP)*Occ(IQ)*
     $ (Two*TwoMO(NAddr3(IP,IP,IQ,IQ))-TwoMO(NAddr3(IP,IQ,IP,IQ)))
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      EndDo
C
C     COMPUTE THE MATRICES A+B, A-B, AND C
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR
      If(IR.Gt.IS) IRS=IRS+1
C
      IPQ=0
      Do IP=1,NBasis
C
      IPS=(Max(IS,IP)*(Max(IS,IP)-1))/2+Min(IS,IP)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Do IQ=1,NBasis
      IQR=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)
      IQS=(Max(IQ,IS)*(Max(IQ,IS)-1))/2+Min(IQ,IS)

      If(IP.Gt.IQ) IPQ=IPQ+1
      If(IQ.Gt.IP) IQP=(IQ**2-3*IQ+2)/2+IP
C
      XMAT=Zero
C
      If(IQ.Eq.IR) XMAT=XMAT+(Occ(IR)-Occ(IS))*HNO(IPS)
      If(IS.Eq.IP) XMAT=XMAT-(Occ(IR)-Occ(IS))*HNO(IQR) 
C
C     INTERGEMINAL CONTRIBUTIONS
C
      IGemPR=1
      If(IGem(IP).Eq.IGem(IR)) IGemPR=0
      IGemQS=1
      If(IGem(IQ).Eq.IGem(IS)) IGemQS=0
      IGemPS=1
      If(IGem(IP).Eq.IGem(IS)) IGemPS=0
      IGemQR=1
      If(IGem(IQ).Eq.IGem(IR)) IGemQR=0
C
      XMAT=XMAT+ (-Occ(IR)*Occ(IP)*IGemPR
     $            +Occ(IS)*Occ(IP)*IGemPS
     $            +Occ(IR)*Occ(IQ)*IGemQR
     $            -Occ(IS)*Occ(IQ)*IGemQS)*
     $  ( Two*TwoMO(NAddr3(IQ,IP,IS,IR))-TwoMO(NAddr3(IQ,IR,IS,IP)) )
C
      If(IQ.Eq.IR) XMAT=XMAT+Occ(IQ)*AuxH(IQ,IPS)-Occ(IS)*AuxH(IS,IPS)
      If(IS.Eq.IP) XMAT=XMAT+Occ(IP)*AuxH(IP,IQR)-Occ(IR)*AuxH(IR,IQR)
C
C     INTRAGEMINAL PART
C
      IGemQR=1
      If(IGem(IQ).Ne.IGem(IR)) IGemQR=0
      IGemPS=1
      If(IGem(IP).Ne.IGem(IS)) IGemPS=0
      XMAT=XMAT+(C(IQ)*C(IR)*IGemQR+C(IP)*C(IS)*IGemPS)
     $ *( TwoMO(NAddr3(IP,IS,IQ,IR))+TwoMO(NAddr3(IP,IR,IQ,IS)) )
C
      If(IS.Eq.IP) XMAT=XMAT-C(IR)*AuxXC(IR,IQR)
      If(IR.Eq.IQ) XMAT=XMAT-C(IS)*AuxXC(IS,IPS)
      If(IR.Eq.IP) XMAT=XMAT-C(IP)*AuxXC(IP,IQS)
      If(IS.Eq.IQ) XMAT=XMAT-C(IQ)*AuxXC(IQ,IPR)
C
      If (IR.Gt.IS.And.IP.Gt.IQ) BMAT(IRS,IPQ)=XMAT
      If (IR.Gt.IS.And.IQ.Gt.IP) AMAT(IRS,IQP)=XMAT
      If (IR.Gt.IS.And.IQ.Eq.IP) ANDMAT(IRS,IP)=XMAT
      If (IR.Eq.IS.And.IP.Gt.IQ) ADNMAT(IR,IPQ)=XMAT
      If (IR.Eq.IS.And.IQ.Eq.IP) ADDMAT(IR,IP)=XMAT
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      IGR=IGem(IR)
C
      Do IS=1,IR-1
      IRS=IRS+1 
      IGS=IGem(IS)
C
      IPQ=0
      Do IP=1,NBasis
      IGP=IGem(IP)    
C
      Do IQ=1,IP-1
      IPQ=IPQ+1
      IGQ=IGem(IQ)
C
      If(IGR.Eq.IGS) Then
C
      If(IGP.Eq.IGR) Then
C
      If(IGQ.Eq.IGR) Then
C     Case: IGR=IGS=IGP=IGQ 
      Call Get_APL(IR,IS,IP,IQ,APL,HNO,XMu(IR),TwoMO,AuxH, 
     $ NBasis,NInte1,NInte2)
      AMAT(IRS,IPQ)=APL
      BMAT(IRS,IPQ)=APL
C
      Else
C     Case IGR=IGS=IGP.Ne.IGQ
c      If(Abs(Abs(C(IP))-Abs(C(IQ))).Lt.Delta*Abs(C(IQ))) 
c     $ Write(*,*)'Err 1'
      Call Get_C(IP,IQ,IR,IS,CC,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      AMAT(IRS,IPQ)=CC/(C(IP)+C(IQ))
      Call Get_D(IQ,IP,IR,IS,DD,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      BMAT(IRS,IPQ)=-DD/(C(IP)-C(IQ))
C
C     endif (IGQ.Eq.IGR)
      EndIf
C
      ElseIf(IGQ.Eq.IGR) Then
C     Case IGR=IGS=IGQ.Ne.IGP
c      If(Abs(Abs(C(IP))-Abs(C(IQ))).Lt.Delta*Abs(C(IQ))) 
c     $ Write(*,*) 'Err 1'
      Call Get_C(IQ,IP,IR,IS,CC,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      AMAT(IRS,IPQ)=CC/(C(IP)+C(IQ))
      Call Get_D(IP,IQ,IR,IS,DD,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      BMAT(IRS,IPQ)=DD/(C(IP)-C(IQ))
C 
      Else
C     Case IGR=IGS.Ne.IGP, IGR=IGS.Ne.IGQ
      AMAT(IRS,IPQ)=(C(IR)-C(IS))*(C(IP)-C(IQ))
     $ *(TwoMO(NAddr3(IQ,IR,IS,IP))-TwoMO(NAddr3(IQ,IS,IP,IR)))
      BMAT(IRS,IPQ)=(C(IR)+C(IS))*(C(IP)+C(IQ))
     $ *(Four*TwoMO(NAddr3(IQ,IP,IS,IR))
     $ -TwoMO(NAddr3(IQ,IR,IS,IP))-TwoMO(NAddr3(IQ,IS,IP,IR)))
C
C     endif (IGP.Eq.IGR)
      EndIf
C
C     Else to (IGR.Eq.IGS)
      Else
C
      If(IGP.Eq.IGQ) Then
C
C     Case IGP=IGQ=IGR.Ne.IGS
      If(IGR.Eq.IGP.And.IGS.Ne.IGP) Then
c      If(Abs(Abs(C(IR))-Abs(C(IS))).Lt.Delta*Abs(C(IR))) 
c     $ Write(*,*) 'Err 2'
      Call Get_C(IR,IS,IP,IQ,CC,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      AMAT(IRS,IPQ)=CC/(C(IR)+C(IS))
      Call Get_D(IS,IR,IP,IQ,DD,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      BMAT(IRS,IPQ)=-DD/(C(IR)-C(IS))
      EndIf
C
C     Case IGP=IGQ=IGS.Ne.IGR
      If(IGS.Eq.IGP.And.IGR.Ne.IGP) Then    
c      If(Abs(Abs(C(IR))-Abs(C(IS))).Lt.Delta*Abs(C(IR)))  
c     $ Write(*,*)'Err 2'
      Call Get_C(IS,IR,IP,IQ,CC,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      AMAT(IRS,IPQ)=CC/(C(IR)+C(IS))
      Call Get_D(IR,IS,IP,IQ,DD,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
      BMAT(IRS,IPQ)=DD/(C(IR)-C(IS))
      EndIf

C     Case IGP=IGQ.Ne.IGR and IGP=IGQ.Ne.IGS
      If(IGR.Ne.IGP.And.IGS.Ne.IGP) Then
      AMAT(IRS,IPQ)=(C(IR)-C(IS))*(C(IP)-C(IQ))
     $ *(TwoMO(NAddr3(IQ,IR,IS,IP))-TwoMO(NAddr3(IQ,IS,IP,IR)))
      BMAT(IRS,IPQ)=(C(IR)+C(IS))*(C(IP)+C(IQ))
     $ *(Four*TwoMO(NAddr3(IQ,IP,IS,IR))
     $ -TwoMO(NAddr3(IQ,IR,IS,IP))-TwoMO(NAddr3(IQ,IS,IP,IR)))
      EndIf
C
C     else to (IGP.Eq.IGQ)
      Else
C     Case IGR.Ne.IGS and IGP.Ne.IGQ 
      If(Abs(Abs(C(IP))-Abs(C(IQ))).Lt.Delta*Abs(C(IQ)).Or. 
     $   Abs(Abs(C(IR))-Abs(C(IS))).Lt.Delta*Abs(C(IR))) Then
      AMAT(IRS,IPQ)=Zero
      BMAT(IRS,IPQ)=Zero
      Else
      SaveA=AMAT(IRS,IPQ)
      SaveB=BMAT(IRS,IPQ)
      AMAT(IRS,IPQ)=(SaveA+SaveB)/(C(IR)+C(IS))/(C(IP)+C(IQ))
      BMAT(IRS,IPQ)=(SaveA-SaveB)/(C(IR)-C(IS))/(C(IP)-C(IQ))   
      EndIf
C
C     endif (IGP.Eq.IGQ)
      EndIf
C
C     endif (IGR.Eq.IGS)
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      If(ISERPA.Eq.1.Or.ISERPA.Eq.2) Then
C
      Do IP=1,NBasis
C
      Do IR=1,IP
C
C     COMPUTE EMAT AND EMATM
C
      EMAT(IP,IR)=Zero
      EMATM(IP,IR)=Zero
c      If(C(IP)*C(IR).Ne.Zero)EMAT(IP,IR)=ADDMAT(IP,IR)/Four/C(IP)/C(IR)
      If(IGem(IP).Eq.IGem(IR)) Then
      EMAT(IP,IR)=TwoMO(NAddr3(IP,IR,IP,IR))
      EMATM(IP,IR)=EMAT(IP,IR)
      Else
      EMATM(IP,IR)=EMATM(IP,IR)+Four*C(IP)*C(IQ)*
     $ (Two*TwoMO(NAddr3(IP,IP,IR,IR))-TwoMO(NAddr3(IP,IR,IP,IR)))
      EndIf
C
      If(IP.Eq.IR) Then
C
      IPP=(IP*(IP+1))/2
      EMAT(IP,IR)=EMAT(IP,IR)+Two*(HNO(IPP)+AuxH(IP,IPP))-XMu(IP)
      EMATM(IP,IR)=EMATM(IP,IR)+Two*(HNO(IPP)+AuxH(IP,IPP))-XMu(IP)
C
      EndIf
C
      EMAT(IR,IP)=EMAT(IP,IR)
      EMATM(IR,IP)=EMATM(IP,IR)
C
C     End of IR LOOP
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
      IRRSS=(Max(IS,IR)*(Max(IS,IR)-1))/2+Min(IS,IR)
C
C     COMPUTE DMAT
C
      DMAT(IRS,IP)=Zero
c      If((C(IR)+C(IS))*C(IP).Ne.Zero) 
c     $ DMAT(IRS,IP)=ANDMAT(IRS,IP)/(C(IR)+C(IS))/C(IP)/Two 
C
C     CASE: IGR=IGS=IGP      
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IP).Eq.IGem(IS)) Then
C
      DMAT(IRS,IP)=TwoMO(NAddr3(IP,IR,IP,IS))
      If(IP.Eq.IR.Or.IP.Eq.IS) DMAT(IRS,IP)=DMAT(IRS,IP)+HNO(IRRSS)
     $ +AuxH(IR,IRRSS) 
C
      EndIf
C
C     CASE IGR.Ne.IGS AND IGR=IGP
C
      If(IGem(IR).Ne.IGem(IS).And.IGem(IP).Eq.IGem(IR)) Then
C
      If(C(IR).Ne.Zero) DMAT(IRS,IP)=
     $ TwoMO(NAddr3(IP,IR,IP,IS))/(One+C(IS)/C(IR))
      If(IP.Eq.IR.And.C(IR)+C(IS).Ne.Zero) DMAT(IRS,IP)=
     $ DMAT(IRS,IP)-AuxXC(IP,IRRSS)/(C(IR)+C(IS))
C
      EndIf
C
C     CASE IGR.Ne.IGS AND IGS=IGP
C
      If(IGem(IR).Ne.IGem(IS).And.IGem(IP).Eq.IGem(IS)) Then
C
      If(C(IS).Ne.Zero) DMAT(IRS,IP)=
     $ TwoMO(NAddr3(IP,IR,IP,IS))/(One+C(IR)/C(IS))
      If(IP.Eq.IS.And.C(IR)+C(IS).Ne.Zero) DMAT(IRS,IP)=
     $ DMAT(IRS,IP)-AuxXC(IP,IRRSS)/(C(IR)+C(IS))
C     
      EndIf
C
C     COMPUTE DMATM
C
      DMATM(IRS,IP)=Zero
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IP).Eq.IGem(IS)) Then
C
      DMATM(IRS,IP)=TwoMO(NAddr3(IP,IR,IP,IS))
      If(IP.Eq.IR.Or.IP.Eq.IS) DMATM(IRS,IP)=DMATM(IRS,IP)+HNO(IRRSS)
     $ +AuxH(IR,IRRSS)
C
      Else
C     
      If(IR.Eq.IP) Then
      DMATM(IRS,IP)=DMATM(IRS,IP)+Two*C(IP)*(HNO(IRRSS)+AuxH(IP,IRRSS))
     $ +AuxXC(IP,IRRSS)
      ElseIf(IS.Eq.IP) Then
      DMATM(IRS,IP)=DMATM(IRS,IP)-Two*C(IP)*(HNO(IRRSS)+AuxH(IP,IRRSS))
     $ -AuxXC(IP,IRRSS)
      EndIf
C
      If(IGem(IP).Eq.IGem(IR)) Then
      DMATM(IRS,IP)=DMATM(IRS,IP)+C(IR)*TwoMO(NAddr3(IP,IR,IP,IS))
      Else 
      DMATM(IRS,IP)=DMATM(IRS,IP)+Two*C(IP)*C(IR)**2*
     $ (Two*TwoMO(NAddr3(IP,IP,IR,IS))-TwoMO(NAddr3(IP,IR,IP,IS)))
      EndIf
C
      If(IGem(IP).Eq.IGem(IS)) Then
      DMATM(IRS,IP)=DMATM(IRS,IP)-C(IS)*TwoMO(NAddr3(IP,IR,IP,IS))
      Else
      DMATM(IRS,IP)=DMATM(IRS,IP)-Two*C(IP)*C(IS)**2*
     $ (Two*TwoMO(NAddr3(IP,IP,IR,IS))-TwoMO(NAddr3(IP,IR,IP,IS)))
      EndIf
C
      If(C(IR)-C(IS).Ne.Zero) Then
      DMATM(IRS,IP)=DMATM(IRS,IP)/(C(IR)-C(IS))
      Else
      DMATM(IRS,IP)=Zero
      EndIf 
C
C     EnidIf of IGem(IP).Eq.IGem(IR).And.IGem(IP).Eq.IGem(IS)
      EndIf
C
C     End of IR,IS LOOPS
      EndDo
      EndDo
C
C     End of IP LOOP
      EndDo
C
C     Else of ISERPA.Eq.1.Or.ISERPA.Eq.2
      Else
C
C     COMPUTE THE CMAT MATRIX DEFINED AS 
C     CMAX = ANDMAT (ADDMAT)^-1 ADNMAT
C
C     INVERT THE ADDMAT MATRIX
C
      Call CpyM(Aux2,ADDMAT,NBasis)
      Call Diag8(Aux2,NBasis,NBasis,XMu,Aux1)
C
      Do I=1,NBasis
      If(XMu(I).Gt.1.E-10) Then
      XMu(I)=One/XMu(I)
      Else
      XMu(I)=Zero
      EndIf
      EndDo
C
      Do I=1,NBasis
      Do J=1,I
      ADDMAT(I,J)=Zero
      Do K=1,NBasis
      ADDMAT(I,J)=ADDMAT(I,J)+Aux2(K,I)*XMu(K)*Aux2(K,J)
      EndDo
      ADDMAT(J,I)=ADDMAT(I,J)
      EndDo
      EndDo
C
      Call MultpMN(Aux1,ADDMAT,ADNMAT,NBasis,NBasis,NBasis,NDim)
      Call MultpMN(CMAT,ANDMAT,Aux1,NDim,NBasis,NBasis,NDim)
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
C
      If(Abs(C(IP)+C(IQ)).Lt.Delta*Abs(C(IQ)).Or.
     $   Abs(C(IR)+C(IS)).Lt.Delta*Abs(C(IR))) Then
      CMAT(IRS,IPQ)=Zero
      Else
      CMAT(IRS,IPQ)=CMAT(IRS,IPQ)/(C(IR)+C(IS))/(C(IP)+C(IQ))
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      EndIf
C
      Return
      End

