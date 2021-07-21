*Deck ACECORR
      Subroutine ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDim,NGOcc,NGem,
     $ IndN,IndX,NDimX)
C
      use abmat
      use abfofo
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
      Integer Points
C     HAP
      Double precision,Allocatable :: WorkVec(:),WorkEig(:),MYAP(:) 
C
      If(IFlSnd.Eq.1) Then
C
      If(ICASSCF.Ne.1) Then
       Call SndOrder(ECorr1,ECorr01,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,Eig,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux)
C
      If(ITwoEl.Eq.2) Then
      Call EneGVB_FFFF(ETot,URe,Occ,CICoef,XOne,
     $                 IGem,IndN,NBasis,NInte1,'TWOMO',NDimX,NGem)
      ElseIf(ITwoEl.Eq.3) Then
      Call EneGVB_FOFO(NActive,NELE,ETot,URe,Occ,CICoef,XOne,
     $                 IGem,IndN,NBasis,NInte1,'FOFO',NDimX,NGem)
      EndIf
C
      Write
     $ (6,'(/,2X,''EGVB+ENuc, 0th+1st-order ECorr, AC0-GVB '',
     $ 4X,3F15.8,/)') ETot+ENuc,ECorr1,ETot+ENuc+ECorr1
      EndIf
C
      If(ICASSCF.Eq.1) Then
C
      If(ITwoEl.eq.1) Then
C      
      Call AC0CAS(ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
C 
      ElseIf(ITwoEl.eq.3) Then
C
      Call AC0CAS_FOFO(ECorr,ETot,Occ,URe,XOne,ABPLUS,ABMIN,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDim,NInte1,
     $ NoSt,'FFOO','FOFO')
C
C     now Y01CAS_FOFO is used in SAPT only
C      Call Y01CAS_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,
C     $ 'PROP0','PROP1',
C     $ 'Y01FILE',
C     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,
C     $ NBasis,NDim,NInte1,NoSt,'EMPTY','FFOO',
C     $ 'FOFO',0,ETot,ECorr)
C
C     ITwoEl
      EndIf
C      
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
C This is not going to work
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
      Call DelInts(ITwoEl)
C
      Return
C
C     If(IFlSnd.Eq.1)
      EndIf
C
C     GENERATE ABSCISSAS AND WEIGHTS FOR GAUSSIAN-LEGENDRE QUADRATURE
C
c      NGrid=30
      NGrid=5
C
      Call GauLeg(Zero,One,XGrid,WGrid,NGrid)
C 
c herer!!!
c      h=5.d-4
c      ACAlpha=0*h
c      Call ACEInteg(ECorr0,TwoNO,URe,Occ,XOne,UNOAO,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ EGOne,NGOcc,
c     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
c     $ IndN,IndX,NDimX) 
c      
c      ACAlpha=h
c      Call ACEInteg(ECorr1,TwoNO,URe,Occ,XOne,UNOAO,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ EGOne,NGOcc,
c     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
c     $ IndN,IndX,NDimX) 
c 
c      ACAlpha=2.*h
c      Call ACEInteg(ECorr2,TwoNO,URe,Occ,XOne,UNOAO,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ EGOne,NGOcc,
c     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
c     $ IndN,IndX,NDimX)
c
c      ACAlpha=3.*h
c      Call ACEInteg(ECorr3,TwoNO,URe,Occ,XOne,UNOAO,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ EGOne,NGOcc,
c     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
c     $ IndN,IndX,NDimX)
c
c      ECorr=ECorr0+0.5*(ECorr1-ECorr0)/h
c     $ +1./6.d0*(ECorr2+ECorr0-2.*ECorr1)/h**2
c     $ +1./24.*(-ECorr0+3.*ECorr1-3.*ECorr2+ECorr3)/h**3
c      write(*,*)'ecorr',ECorr
c
c      stop

      ECorr=Zero
      Do I=1,NGrid
C   
      ACAlpha=XGrid(I)
C
      If(ITwoEl.eq.1) Then

      Call ACEInteg(ECorrA,TwoNO,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,
     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
      ElseIf(ITwoEl.eq.3) Then

      If(ICASSCF.Eq.1) Then

      Call ACEInteg_FOFO(ECorrA,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,CICoef,
     $ NBasis,NInte1,NDimX,NGem,IndAux,ACAlpha,
     $ IGem,NAcCAS,NInAcCAS,NELE,IndN,IndX,NDimX,
     $ NoSt,ICASSCF,IFlFrag1,IFunSR,IFunSRKer)

      ElseIf(ICASSCF.Ne.1) Then

      Call ACEInteg_FOFO(ECorrA,URe,Occ,XOne,UNOAO,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ EGOne,NGOcc,CICoef,
     $ NBasis,NInte1,NDim,NGem,IndAux,ACAlpha,
     $ IGem,NActive,NInAcCAS,NELE,IndN,IndX,NDimX,
     $ NoSt,ICASSCF,IFlFrag1,IFunSR,IFunSRKer)

      EndIf
      EndIf
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
      If(IFunSR.Eq.0) Then
      Write
     $ (6,'(/,2X,''ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
      Else
      EGOne(1)=ECorr
      EndIf
C
      Else 
C
      If(ITwoEl.Eq.3) Then
      Call EneGVB_FOFO(NActive,NELE,ETot,URe,Occ,CICoef,XOne,
     $                 IGem,IndN,NBasis,NInte1,'FOFO',NDimX,NGem)
      EndIf
C
      Write
     $ (6,'(/,2X,''EGVB+ENuc, Corr, AC-ERPA-GVB '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      EndIf
C
      Call DelInts(ITwoEl)
C
      return
      GoTo 888
c      stop 
C
C     SCANNING THROUGH ACAlpha
C
      ACAlpha=1.0
C
      Points=210
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
      stop
C    
C     SCAN CLOSE TO X
C
  888 Continue

      ACAlpha=1.D0
C
c      Call ACEInteg(ECorr,TwoNO,URe,Occ,XOne,UNOAO,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ EGOne,NGOcc,
c     $ Title,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
c     $ IndN,IndX,NDimX)
cC
c      Write(6,'(2X,''Alpha, W_ALPHA '',4X,2F15.8)')ACAlpha,ECorr
c      stop

c      ACAlpha0=0.985
      ACAlpha0=0.8
C
      Points=50
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
      stop
      Return
      End

*Deck Y01CAS
      Subroutine Y01CAS(TwoNO,Occ,URe,XOne,ABPLUS,ABMIN,
     $ EigY,EigY1,Eig,Eig1,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,IFlag0,
     $ filename)
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
C UNCOUPLED 01/02/2021
      Character(*),Optional :: filename
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
c NEW 11/07/2018
     $ IMatch(NDimX)
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
   10 Read(10,*,End=40)I,J,K,L,X
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
      Write(6,'(" *** DONE WITH 0TH-ORDER IN Y01CAS ***")')
C
C     DONE 0TH-ORDER CALCULATIONS  
C
C NEW 11/07/2018
C     Check if NoEig=NDimX - they should be equal!
      If(NoEig.Ne.NDimX) Stop 'Fatal error in Y01CAS: NoEig.Ne.NDimX!'
C
C     UNCOUPLED
      Call ConvertXYtilde(EigY,EigX,Eig,CICoef,IEigAddY,IEigAddInd,
     $                    IndN,IndBlock,NFree2,NBasis,NDimX,NDim,
     $                    filename,.true.)
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
C
      Return
      EndIf
C
      Write(6,'(" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
C
      Deallocate(RDM2Act)
C
      Write(6,'(" *** DONE WITH COMPUTING AB(1) MATRICES ***")')
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
C
      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12)) Aux2=
     $ (ABPLUS(MU+(NU-1)*NoEig)+ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)-Eig(NU))
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
      Return
      End

      Subroutine ConvertXYtilde(EigY,EigX,Eig,C,IEigAddY,IEigAddInd,
     $ IndN,IndBlock,NFree2,NBasis,NDimX,NDim,filename,ifdump)
C
      use abfofo,only : reduce_to_XY0CAS
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Integer :: NFree2,NBasis,NDimX,NDim
      Integer :: IndN(2,NDim),IndBlock(2,NDimX),
     $           IEigAddY(2,NDimX),IEigAddInd(2,NDimX)
      Double precision :: EigY(NFree2),EigX(NFree2),
     $                    Eig(NDimX),C(NBasis)
      Logical      :: ifdump
      Character(*) :: filename
      Double precision,allocatable :: Work(:),WorkY(:),WorkX(:)
C
      Integer :: IMatch(NDimX)
C
      NAct=NAcCAS
      INActive=NInAcCAS
C
      NoEig=NDimX
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
      Allocate(Work(NDimX**2),WorkY(NDimX**2),WorkX(NDimX**2))
C
      Do I=1,NoEig
      Do J=1,NoEig
      Work(I+(J-1)*NoEig)=0d0
      If(I.Eq.J) Work(I+(J-1)*NoEig)=1d0
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      WorkY(NU+(MU-1)*NoEig)=0d0
      WorkX(NU+(MU-1)*NoEig)=0d0
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      WorkY(NU+(MU-1)*NoEig)=WorkY(NU+(MU-1)*NoEig)
     $ +Work(NU+(I-1)*NoEig)*EigY(IStart+II)
      WorkX(NU+(MU-1)*NoEig)=WorkX(NU+(MU-1)*NoEig)
     $ +Work(NU+(I-1)*NoEig)*EigX(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C  
c     RESORT ACCORDING TO IndN
C
      Call CpyM(Work,WorkY,NDimX)
      Do MU=1,NDimX
      Do I=1,NDimX
      WorkY((MU-1)*NoEig+I)=Work((MU-1)*NoEig+IMatch(I))
      EndDo
      EndDo
C
      Call CpyM(Work,WorkX,NDimX)
      Do MU=1,NDimX
      Do I=1,NDimX
      WorkX((MU-1)*NoEig+I)=Work((MU-1)*NoEig+IMatch(I))
      EndDo
      EndDo
C
      Do K=1,NDimX
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
C
      X=WorkX((K-1)*NDimX+I)/(C(IP)+C(IQ))
      Y=WorkY((K-1)*NDimX+I)/(C(IP)-C(IQ))
C
      WorkX((K-1)*NDimX+I)=0.5d0*(X-Y)
      WorkY((K-1)*NDimX+I)=0.5d0*(X+Y)
C
      EndDo
      EndDo
C
c      Print*, 'TEST-X:',norm2(WorkX)
c      Print*, 'TEST-Y:',norm2(WorkY)
C
      If(ifdump) Then
C      
C      Open(newunit=iunit,file=filename,form='unformatted')
C      write(iunit) WorkX
C      write(iunit) WorkY
C      Close(iunit)
C
      Call reduce_to_XY0CAS(WorkX,WorkY,Eig,C,
     $  IndN,NAct,INActive,NDimX,NBasis,filename)
      Else
C     RETURN X(0) and Y(0)
C
      EigX = WorkX
      EigY = WorkY
C
      EndIf
C
      Deallocate(WorkX,WorkY,Work)
C
      End Subroutine

*Deck Y01CASLR
      Subroutine Y01CASLR(TwoNO,Occ,URe,XOne,ABPLUS,ABMIN,
     $ EigY,EigY1,Eig,Eig1,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,IFlag0
C WARNING!!!!! STH WRONG WITH THIS PROCEDURE! CHANGE ADDING sr Kernel!
C   
C DFT START 
C
C     it is assumed that erf/r integrals are in TwoNO, 
C     while TwoEl2 stores 1/r integrals
C     SRKerW(I) = WGrid(I)*SRKer(I), SRKer: sr-alda kernel 
C
     $ ,TwoEl2,OrbGrid,SRKerW,NSymNO,MultpC,NGrid)
C DFT END
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
C DFT START
     $ ,TwoEl2(NInte2),OrbGrid(NBasis*NGrid),SRKerW(NGrid),
     $ NSymNO(NBasis),MultpC(15,15)
C DFT END
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
c NEW 11/07/2018
     $ IMatch(NDimX)
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
   10 Read(10,*,End=40)I,J,K,L,X
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
C DFT START
C
C     ADD A SR KERNEL
C
      I1I2S=MultpC(NSymNO(IP),NSymNO(IQ))
      I3I4S=MultpC(NSymNO(IR),NSymNO(IS))
      ISym=MultpC(I1I2S,I3I4S)
      XKer1234=Zero
      If(ISym.Eq.1) Then
      Do I=1,NGrid
      XKer1234=XKer1234+SRKerW(I)*
     $ OrbGrid(IP+(I-1)*NBasis)*OrbGrid(IQ+(I-1)*NBasis)*
     $ OrbGrid(IR+(I-1)*NBasis)*OrbGrid(IS+(I-1)*NBasis)
      EndDo
      EndIf
      TwoSR=TwoEl2(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IQ,IR,IS))
      ABM=ABM+Four*(C(IP)+C(IQ))*(C(IR)+C(IS))*(XKer1234+TwoSR)
C
C DFT END
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
C
      Return
      EndIf
C
      Write(6,'(" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
C
C DFT START
C
C     ADD A SR KERNEL
C
      Write(6,'(/," ***  adding a sr-kernel ***" )')
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
     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
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
      Write(6,'("*** sr-kernel added ***")')
C
C DFT END
C
      Deallocate(RDM2Act)
C
      Write(6,'(" *** DONE WITH COMPUTING AB(1) MATRICES ***")')
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
      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12)) Aux2=
c      If(MU.Ne.NU) Aux2=
     $  (ABPLUS(MU+(NU-1)*NoEig)+ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)-Eig(NU))
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
      Return
      End

*Deck AC0CAS
      Subroutine AC0CAS(ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ ABPLUS,ABMIN,EigY,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
C      use sapt_ener
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
     $ XMAux(NDimX*NDimX),work1(NBasis,NBasis)
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
   10 Read(10,*,End=40)I,J,K,L,X
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
C      Print*, 'ACT-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
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
C      Print*, 'AI-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
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
C      Print*, 'AV-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN(1:NDimB**2))
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
C     Do IP
      EndDo
C      
C     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-INACTIVE BLOCKS

C KP 15.05.2019
      if(idalton.eq.0) then
      open(10,file='fock.dat')
      work1=0
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
      read(10,*)iip,iiq,xx
      work1(iip,iiq)=xx
      EndDo
      EndDo
      close(10)
      endif
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
C
      If(IDALTON.Eq.0) Then
      If(ABS(ABP-work1(IP,IQ)).Gt.1.d-7)
     $Write(*,*)'ABP inconsistent with eps_a-eps_i for',IP,IQ
      EndIf
C
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
      Print*, 'NoEig,NDimX',NoEig,NDimX
C      Print*, 'Eig,Y,X',norm2(Eig(1:NoEig)),
C     $ norm2(EigY(1:NoEig**2)),norm2(EigX(1:NoEig**2))
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
      Print*, 'AB1-Ka',norm2(ABPLUS),norm2(ABMIN)
C ----------------------------------------------------------------    
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
      print*, 'AB-KA',norm2(ABPLUS),norm2(ABMIN)
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
C      Print*, 'FIRST',norm2(XMAux)
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
      Print*, 'ABPLUS-Ka',norm2(ABPLUS(1:NoEig*NoEig))
C
C     FINALLY THE ENERGY CORRECTION
C     TESTY Z sapt.f90 -- remove later! 
C      Call check_loop(ABPLUS,Occ,IndN,IndBlock,
C     $ NAct,INActive,NDimX,NDimX,NBasis)
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
      ECorr=EAll-EIntra
      Print*, 'EAll,EIntra',EAll,EIntra
C
C     MP2 energy (only inactive-virtual enter)
C     to compare with molpro use {mp2;core,0}
c      EMP2=Zero
c      Do IP=NOccup+1,NBasis
c      Do IR=1,INActive
c      Do IQ=NOccup+1,NBasis
c      Do IS=1,INActive
c      EMP2=EMP2-(Two*TwoNO(NAddr3(IP,IR,IQ,IS))
c     $ -TwoNO(NAddr3(IP,IS,IQ,IR)))*TwoNO(NAddr3(IP,IR,IQ,IS))/
c     $ (work1(IP,IR)+work1(IQ,IS))
cC
c      EndDo
c      EndDo
c      EndDo
c      EndDo
c      Write(6,'(/,1x,a,13x,f16.8)') 'EMP2   Energy', EMP2 
C
C ----------------------------------------------------------------    
      Return
      End

*Deck AC0CDEXCIT
      Subroutine AC0CDEXCIT(IH0St,ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ UNOAO,IndN,IndX,NBasis,NDimX,NDimD,NInte1,NInte2) 
C
C     AC0 AND DEEXCITATION CORRECTIONS BASED ON OVERLAP OF ERPA AND SA-CAS TRDM's
C
C      use sapt_ener
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
     $ IndX(NDimX),IndN(2,NDimX),
     $ UNOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ ABPLUS(NDimD*NDimD),ABMIN(NDimD*NDimD),
     $ Eig(NDimD),EigY(NDimD*NDimD),
     $ EigX(NDimD*NDimD),
     $ IEigAddY(2,NDimD),IEigAddInd(2,NDimD),IndBlock(2,NDimD),
     $ XMAux(NDimD*NDimD),work1(NBasis,NBasis),
     $ TrGamm(NInte1,NInte1),EExcit(NInte1),
     $ XCAS(NBasis,NInte1),YCAS(NBasis,NInte1),
     $ GammaS(100,NInte1)
     $ ,IMatch(NDimX),EigY1(NDimD*NDimD),
c EigX1(NDimD*NDimD),
     $ Eig1(NDimD),
     $ DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
C
      IStERPA=0
C
      If(IH0St.Lt.NoSt) Then
      Write(6,'(/,X,"IH0St=",I4,"  NoSt=",I4)')IH0St,NoSt
      Stop 'Stop in AC0CDEXCIT: IH0St is smaller 
     $ than NoSt. AC0 correction would not be reliable '
      EndIf
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
      Call ReadDip(DipX,DipY,DipZ,UNOAO,NBasis)
C
      NoStMx=0
      Write(6,'(X,"**** SA-CAS FROM MOLPRO ****",/)')
      IPr=0
      If(IH0St.Eq.NoSt) IPr=1
      Call RDM_SACAS(GammaS,XCAS,YCAS,EExcit,C,UNOAO,IPair,
     $ DipX,DipY,DipZ,NoSt,NoStMx,NInte1,NBasis,IPr)
      Write(6,'(X,"The number of states in SA-CAS: ",I4,/)')NoStMx
      If(IH0St.Gt.NoStMx) Stop'Stop in AC0CDEXCIT: IH0St is greater
     $ than the number of states in SA-CAS!'
      Do NU=1,NoSt-1
      Write(6,'(X,"SA-CAS Deexcitation from ",I4," to",I4,2E15.6)') 
     $ NoSt,NU,EExcit(NU)
      EndDo
      Do NU=NoSt+1,NoStMx
      Write(6,'(X,"SA-CAS Excitation   from ",I4," to",I4,2E15.6)') 
     $ NoSt,NU,EExcit(NU)
      EndDo
C
      If(IH0St.Ne.NoSt.And.EExcit(IH0St).Eq.Zero) Then
      Write(6,'(/,X,"Transition energy of interest is zero, quitting")')
      ECorr=Zero
      Return
      EndIf
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
   10 Read(10,*,End=40)I,J,K,L,X
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
      Write(6,'(/,X,"*** H0 constructed for SA-CAS state no",I4,
     $ " ***")') IH0St
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
      Do IT=1,NOccup
      Do IR=1,NOccup
      ITR=(Max(IT,IR)*(Max(IT,IR)-1))/2+Min(IT,IR)
      If(IGem(IT).Ne.IGem(I).And.IGem(IR).Ne.IGem(I))
     $ Aux=Aux+GammaS(IH0St,ITR)*
     $ (Two*TwoNO(NAddr3(IT,IR,I,J))-TwoNO(NAddr3(IT,I,IR,J)))
      EndDo
      EndDo
 
c      Do IT=1,NBasis
c      If(IGem(IT).Ne.IGem(I))
c     $ Aux=Aux+Occ(IT)*
c     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
c      EndDo
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
      Write(6,'(/," *** ACTIVE-ACTIVE BLOCK ***")')
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
c      If(IPair(IP,IQ).Eq.1.And.IPair(IR,IS).Eq.1) Then
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABMIN((IRow-1)*NDimB+ICol)=ABM
      Else
      ABPLUS((ICol-1)*NDimB+IRow)=Zero
      ABPLUS((IRow-1)*NDimB+ICol)=Zero
      ABMIN((ICol-1)*NDimB+IRow)=Zero
      ABMIN((IRow-1)*NDimB+ICol)=Zero
c      EndIf
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
C
c      Print*, 'ACT-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
C
      If(NoSt.Eq.1) Then
C
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
C     
c      Write(6,'(X,"Active Eigenvalues")')
c      Do NU=1,NDimB
c      Write(6,'(X,I4,E15.6)') NU,Eig(NFree1-1+NU)
c      EndDo
cC
c      Do NU=1,NDimB
cC
c      Write(6,'(X,I4,E15.6)') NU,Eig(NFree1-1+NU)
c      Do I=1,NDimB
c      if(abs(EigY(NFree2-1+(NU-1)*NDimB+I)).gt.1.d-7)
c     $ write(*,*)IndBlock(1,NFree1-1+i),
c     $ IndBlock(2,NFree1-1+i),
cc     $ (cicoef(IndBlock(1,NFree1-1+i))+cicoef(IndBlock(2,NFree1-1+i)))*
c     $ EigY(NFree2-1+(NU-1)*NDimB+I),EigX(NFree2-1+(NU-1)*NDimB+I)
cC
c      EndDo
cC
c      EndDo
C
C     If(NoSt.Eq.1) Then
      Else
C
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
C
      Call SortEigXY(1,Eig(NFree1),EigY(NFree2),EigX(NFree2),NDimB) 
C
C 15.07.2020 KP
C SET TO ZERO PSEUDO-DEEXCITATIONS AND LOAD THE CORRECT ONES (IN THE ACTIVE SPACE)
C
c      Call TRDM_SACAS(XCAS,YCAS,NoState,EExcit,C,IPair,NInte1,NBasis)
c      If(NoSt.Ne.NoState) Stop 'Fatal error 0 in TRDM_SACAS'
C
C
C     USE TRDM FROM SA-CAS IN AC0
C
      goto 555
      Write(6,'(/)')
      Do NU0=1,NoState-1
C
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Eq.Zero) Then
C
      Eig(NFree1-1+NU)=EExcit(NU0)
      Write(6,'(X,"Including negative excitation",E15.6)')EExcit(NU0)
C
      SumNu=Zero
      Do I=1,NDimB
C
      IA=IndBlock(1,NFree1-1+I)
      IB=IndBlock(2,NFree1-1+I)
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      EigY(NFree2-1+(NU-1)*NDimB+I)=YCAS(NU0,IAB)
      EigX(NFree2-1+(NU-1)*NDimB+I)=XCAS(NU0,IAB)
      If(Abs(YCAS(NU0,IAB))+Abs(XCAS(NU0,IAB)).Gt.1.D-6) 
     $ Write(6,'(X,"EigY EigX",2I3,2E15.6)')IA,IB,
     $ YCAS(NU0,IAB),XCAS(NU0,IAB)
      SumNu=SumNu+YCAS(NU0,IAB)*XCAS(NU0,IAB)
C 
      EndDo
C
      Write(6,'(X,"Norm Y*X for included excit  ",E15.6)')SumNu
C
      GoTo 987
C 
      EndIf
C
C     Do NU
      EndDo
C
  987 Continue
C     Do NU0
      EndDo
      Write(6,'(/)')
  555 continue
C
C FINE DELLA NOVITA
C
C     If(NDimB.Ne.0) Then
      EndIf
C     If(NoSt.Eq.1) Then 
      EndIf
C
      Write(6,'(X,"Active ERPA Eigenvalues and Eigenvecs")')
C
      IPositive=0
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Gt.Zero) Then
      IPositive=IPositive+1
      If(IPositive.Eq.IH0St-NoSt) IStERPA=NU
      EndIf
C
      Write(6,'(X,I4,E15.6)') NU,Eig(NFree1-1+NU)
      Do I=1,NDimB
      If(Abs(EigY(NFree2-1+(NU-1)*NDimB+I))
     $  +Abs(EigX(NFree2-1+(NU-1)*NDimB+I)).Gt.1.d-7)
     $ Write(6,'(X,"Y_ERPA, X_ERPA",2I3,2E15.6)')IndBlock(1,NFree1-1+I),
     $ IndBlock(2,NFree1-1+I),
     $ EigY(NFree2-1+(NU-1)*NDimB+I),EigX(NFree2-1+(NU-1)*NDimB+I)
      EndDo
C
      EndDo
C
C     set negative eigs to zero only if full AC0 correction is computed for a given state 
C
      If(IH0St.Eq.NoSt) Then
C
      Do NU=1,NDimB
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
      EndDo
C
      EndIf
C
      If(IH0St.Ne.NoSt) Then
C
c      If(IStERPA.Ne.Zero) Then
C
c      Write
c     $ (6,'(/,X,''ERPA Excitation to state '',I2,
c     $ " corresponding to ERPA state ",I2," reads",F15.8)')
c     $ IH0St,IStERPA,Eig(NFree1-1+IStERPA)
c      NU=IStERPA
cC
c      Do I=1,NDimB
c      IA=IndBlock(1,NFree1-1+I)
c      IB=IndBlock(2,NFree1-1+I)
c      If(Abs(EigY(NFree2-1+(NU-1)*NDimB+I))
c     $  +Abs(EigX(NFree2-1+(NU-1)*NDimB+I)).Gt.1.d-7)
c     $ Write(6,'(X,"YCAS, XCAS    ",2I3,2E15.6)')IA,IB,
c     $ EigY(NFree2-1+(NU-1)*NDimB+I),EigX(NFree2-1+(NU-1)*NDimB+I) 
c      EndDo
C      
c      EndIf
C
C     Find X,Y of the highest overlap with YCAS(IH0St),XCAS(IH0St)
C
      Write(6,*)
      OvMax=Zero
      OvXMax=Zero
      OvYMax=Zero
      Do NU=1,NDimB
C
c      SOv=Zero
C
      SumERY=Zero
      SumERX=Zero
      SumCAY=Zero
      SumCAX=Zero
      Do I=1,NDimB
      IA=IndBlock(1,NFree1-1+I)
      IB=IndBlock(2,NFree1-1+I)
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
C
      YER=EigY(NFree2-1+(NU-1)*NDimB+I)
      YCA=YCAS(IH0St,IAB)
      XER=EigX(NFree2-1+(NU-1)*NDimB+I)
      XCA=XCAS(IH0St,IAB)
C
      SumERX=SumERX+XER*XER
      SumCAX=SumCAX+XCA*XCA
      SumERY=SumERY+YER*YER
      SumCAY=SumCAY+YCA*YCA
      EndDo
      If(Abs(SumERX).Gt.1.D-12) SumERX=One/Sqrt(Abs(SumERX))
      If(Abs(SumCAX).Gt.1.D-12) SumCAX=One/Sqrt(Abs(SumCAX))
      If(Abs(SumERY).Gt.1.D-12) SumERY=One/Sqrt(Abs(SumERY))
      If(Abs(SumCAY).Gt.1.D-12) SumCAY=One/Sqrt(Abs(SumCAY))
C
      SOvY=Zero
      SOvX=Zero
C
      Do I=1,NDimB
      IA=IndBlock(1,NFree1-1+I)
      IB=IndBlock(2,NFree1-1+I)
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      YER=EigY(NFree2-1+(NU-1)*NDimB+I)
      YCA=YCAS(IH0St,IAB)
      XER=EigX(NFree2-1+(NU-1)*NDimB+I)
      XCA=XCAS(IH0St,IAB)
c      SOv=SOv+ Abs(YER*YCA)+Abs(XER*XCA)
      SOvY=SOvY+YER*YCA*SumERY*SumCAY
      SOvX=SOvX+XER*XCA*SumERX*SumCAX 
      EndDo
C
      SOvY=Abs(SOvY)
      SOvX=Abs(SOvX)
      Write(6,'(X,"SA-CAS-Overlap for ERPA vector",I2,3F15.8)')
     $ NU,SOvY,SOvX,(SOvY+SovX)/Two
C
      If((SOvY+SovX)/Two.Gt.OvMax) Then
      OvMax=(SOvY+SovX)/Two
      OvXMax=SOvX
      OvYMax=SOvY
      NUMx=NU
      EndIf
C
      EndDo  
C
      If((OvMax.Ge.0.5).Or.
     $ ( (OvMax.Lt.0.5).And.(OvXMax.Gt.0.5.Or.OvYMax.Gt.0.5) )) Then
C
      IStERPA=NUMx
C
      Write
     $ (6,'(/,X,"ERPA vector best overlapping with SA-CAS is vector no",
     $ I2," of excit energy:",F15.8)')
     $ IStERPA,Eig(NFree1-1+IStERPA)
      NU=IStERPA
C
      Do I=1,NDimB
      IA=IndBlock(1,NFree1-1+I)
      IB=IndBlock(2,NFree1-1+I)
      If(Abs(EigY(NFree2-1+(NU-1)*NDimB+I))
     $  +Abs(EigX(NFree2-1+(NU-1)*NDimB+I)).Gt.1.d-7)
     $ Write(6,'(X,"Y_ERPA, X_ERPA        ",2I3,2E15.6)')IA,IB,
     $ EigY(NFree2-1+(NU-1)*NDimB+I),EigX(NFree2-1+(NU-1)*NDimB+I)
      EndDo 
C
      Write(6,'(X,"best-matching SA-CAS Excitation from ",
     $ I4," to ",I4,23X,F15.8)') NoSt,IH0St,EExcit(IH0St)
      Do I=1,NDimB
      IA=IndBlock(1,NFree1-1+I)
      IB=IndBlock(2,NFree1-1+I)
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      If(Abs(YCAS(IH0St,IAB))+Abs(XCAS(IH0St,IAB)).Gt.1.d-7)
     $ Write(6,'(X,"Y_SA-CAS, X_SA-CAS    ",2I3,2E15.6)')IA,IB,
     $ YCAS(IH0St,IAB),XCAS(IH0St,IAB)
      EndDo
C
      Else
C
      IStERPA=0
      Write(6,'(/,X,
     $ "No ERPA vector overlaps with SA-CAS State No",I2)')IH0St
C     If(OvMax.Ge.0.7) Then
      EndIf
C
C     If(IH0St.Ne.NoSt) Then
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS
C
      Write(6,'(/," *** ACTIVE-INACTIVE BLOCKS ***")')
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
c      Print*, 'AI-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
C
      EndDo
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
C      Print*, 'AV-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN(1:NDimB**2))
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
C
      EndDo
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
      Print*, 'NoEig,NDimD',NoEig,NDimD
C
C     DONE 0TH-ORDER CALCULATIONS
C
      Write(6,'(/,
     $" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimD,NBasis,NInte1,NInte2)
C
      Deallocate(RDM2Act)
C
      Write(6,'(/," *** DONE WITH COMPUTING AB(1) MATRICES ***")')
C
C ----------------------------------------------------------------
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
      If(IH0St.Ne.NoSt.And.IStERPA.Ne.0) 
     $ Write(6,'(/,X,
     $ "Deexcitation correction is computed for ERPA vector no ",
     $ I2," Eig=",F15.8)')  IStERPA,Eig(IStERPA)
      If(IH0St.Ne.NoSt.And.IStERPA.Eq.0) 
     $ Write(6,'(/," ERPA vector for deexcitation correction could not 
     $ be determined. The correction will be set to 0.")')
C
      Do MU=1,NoEig
      If(Eig(MU).Ne.Zero) Then
C
      Do NU=1,NoEig
      If(Eig(NU).Ne.Zero) Then
C
      If(IH0St.Eq.NoSt.Or.(NU.Eq.IStERPA.Or.MU.Eq.IStERPA)) Then 
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
C     If(IH0St.Ne.NoSt...
      EndIf
C
C     If(Eig(NU).Ne.Zero)
      EndIf
      EndDo
C
      EndIf
      EndDo
C
C     COMPUTE 1st-ORDER Y AND X PART 1
C
      EigY1(1:NoEig*NoEig)=Zero
c      EigX1(1:NoEig*NoEig)=Zero
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
C
      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12)) Aux2=
     $ (ABPLUS(MU+(NU-1)*NoEig)+ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)-Eig(NU))
C
      EigY1(I+(MU-1)*NoEig)=EigY1(I+(MU-1)*NoEig)+
     $(Aux1+Aux2)*EigY(IStart+II)
C
c      EigX1(I+(MU-1)*NoEig)=EigX1(I+(MU-1)*NoEig)+
c     $(Aux2-Aux1)*EigX(IStart+II)
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
C     CONTINUATION OF AC0
C
      ABPLUS(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C
      IStart=IEigAddY(1,MU)
C
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
      ECorr=EAll-EIntra
C
C     TRANSITION DIPOLE MOMENTS ARE COMPUTED ONLY FOR A GROUND STATE (NoSt=1)
C
      If(NoSt.Eq.1) Then
C
C     Y(0) AND X(0) IN COLUMNS
C
      Do I=1,NoEig
      Do J=1,NoEig
      ABPLUS(I+(J-1)*NoEig)=Zero
      If(I.Eq.J) ABPLUS(I+(J-1)*NoEig)=One
      EndDo
      EndDo
C
      Do I=1,NFree2
      ABMIN(I)=EigY(I)
c      XMAux(I)=EigX(I)
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      EigY(NU+(MU-1)*NoEig)=Zero
c      EigX(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      EigY(NU+(MU-1)*NoEig)=EigY(NU+(MU-1)*NoEig)
     $ +ABPLUS(NU+(I-1)*NoEig)*ABMIN(IStart+II)
c      EigX(NU+(MU-1)*NoEig)=EigX(NU+(MU-1)*NoEig)
c     $ +ABPLUS(NU+(I-1)*NoEig)*XMAux(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
C     SORT Y0,X0 AND Y1,X1 ACCORDING TO IndN
C
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
      Call CpyM(ABPLUS,EigY,NDimX)
      Call CpyM(ABMIN,EigY1,NDimX)
C
      Do MU=1,NDimX
      Do I=1,NDimX
      EigY((MU-1)*NoEig+I)=ABPLUS((MU-1)*NoEig+IMatch(I))
      EigY1((MU-1)*NoEig+I)=ABMIN((MU-1)*NoEig+IMatch(I))
      EndDo
      EndDo
C
c      Call CpyM(ABPLUS,EigX,NDimX)
c      Call CpyM(ABMIN,EigX1,NDimX)
C
c      Do MU=1,NDimX
c      Do I=1,NDimX
c      EigX((MU-1)*NoEig+I)=ABPLUS((MU-1)*NoEig+IMatch(I))
c      EigX1((MU-1)*NoEig+I)=ABMIN((MU-1)*NoEig+IMatch(I))
c      EndDo
c      EndDo
C
      If(IStERPA.Ne.0) Then
C ???
      Write(6,'(/,X,"Y(0) corresponding to SA-CAS best-matching
     $ vector no",I2)')IStERPA
      Call TrDipMoms(IStERPA,TDIP2,EigY,
     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,I2,"->",I2,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')
     $ NoSt,IH0St,TDIP2
C
      Write(6,'(X,"Y(0)+Y(1) corresponding to SA-CAS best-matching
     $ vector no",I2)')IStERPA
C
      Do I=1,NdimX
      XMAux((IStERPA-1)*NoEig+I)=EigY((IStERPA-1)*NoEig+I)
     $ +EigY1((IStERPA-1)*NoEig+I)
      EndDo
      Call TrDipMoms(IStERPA,TDIP2,XMAux,
     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,I2,"->",I2,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')
     $ NoSt,IH0St,TDIP2
C
      Else
C
      If(IH0St.Ne.NoSt) Then 
      Write(6,'(X,I2,"->",I2,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')
     $ NoSt,IH0St,Zero
      Write(6,'(X,I2,"->",I2,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')
     $ NoSt,IH0St,Zero 
      EndIf

c      IStERPA=1
C
c      Write(6,'(/,X,"Y(0) corresponding to vector no",I2)')
c     $ IStERPA
c      Call TrDipMoms(IStERPA,TDIP2,EigY,
c     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
c      Write(6,'(X,"Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')TDIP2
c      Write(6,'(X,"Y(0)+Y(1) corresponding to vector no",I2)')
c     $ IStERPA
cC
c      Do I=1,NdimX
c      XMAux((IStERPA-1)*NoEig+I)=EigY((IStERPA-1)*NoEig+I)
c     $ +EigY1((IStERPA-1)*NoEig+I)
c      EndDo 
c      Call TrDipMoms(IStERPA,TDIP2,XMAux,
c     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
c      Write(6,'(X,"Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')TDIP2
C
C     If( IStERPA.Ne.0) Then
      EndIf
C
c print out for testing       
C
c      Do NU=1,5
cC
cc      Write(6,'(X,"Eig(0), Eig(1), Eig(01) ",I3,3E15.6)') 
cc     $ NU,Eig(NU),Eig1(NU),Eig(NU)+Eig1(NU)
cC
c      Sum0=Zero
c      Sum01=Zero
c      Do I=1,NDimX
cC
c      Y0=EigY((NU-1)*NoEig+I)
c      Y1=EigY1((NU-1)*NoEig+I)
c      X0=EigX((NU-1)*NoEig+I)
c      X1=EigX1((NU-1)*NoEig+I)
cc      If(Abs(Y0)+Abs(Y1).Gt.1.D-7) 
cc     $  Write(6,'(X,"Y(0),Y(1),X(0),X(1)",2I3,4E15.6)')
cc     $ IndN(1,I),IndN(2,I),Y0,Y1,X0,X1
cC
c      Sum0=Sum0+X0*Y0
c      Sum01=Sum01+X0*Y1+X1*Y0
cC
c      EndDo
cC
c      if(abs(sum01).gt.1.d-10) Write(6,'(X,"NU Sum01 ",I2,E15.6,/)')
c     $ nu,Sum01
cc      Write(6,'(X,"Sum0, Sum01 ",2E15.6,/)')
cc     $ Sum0,Sum01
cC
c      EndDo
cC
C     If(NoSt.Eq.1) Then
      EndIf
C
      Return
      End

*Deck AB1_CAS
      Subroutine AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
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
     $ RDM2Act(NRDM2Act),IGFact(NInte2),C(NBasis),
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
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) Then
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
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
C
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) Then
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
      icount=0
      Call CPU_TIME(START_TIME)
C
C     HAP-TEST!
C      IFunSR=4
C
      Do IRow=1,NoEig
C
      IR=IndBlock(1,IRow)
      IS=IndBlock(2,IRow)
C
      If(IFunSR.Eq.4) Then
      IColEnd=NoEig
      Else
      IColEnd=IRow
      EndIf
C
      Do ICol=1,IColEnd
C
      IPP=IndBlock(1,ICol)
      IQQ=IndBlock(2,ICol)
C
      If( .NOT.(IGem(IR).Eq.IGem(IS).And.IGem(IR).Eq.IGem(IPP)
     $ .And.IGem(IR).Eq.IGem(IQQ)) ) Then
CC
      If( (Occ(IR)*Occ(IS).Eq.Zero.And.Occ(IPP)*Occ(IQQ).Eq.Zero
     $ .And.Abs(TwoNO(NAddr3(IR,IS,IPP,IQQ))).Lt.1.D-25)
     $.Or.
     $((Occ(IR).Eq.One.Or.Occ(IS).Eq.One)
     $ .And.
     $ (Occ(IPP).Eq.One.Or.Occ(IQQ).Eq.One)
     $ .And.Abs(TwoNO(NAddr3(IR,IS,IPP,IQQ))).Lt.1.D-25)) Then
C
      icount=icount+1
C
      Else  

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
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=One
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
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0)
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
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) 
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
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0)
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
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0)
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
C     end icount
      EndIf
c     If IGem ....
      EndIf
C
      EndDo
      EndDo
      Call CPU_TIME(END_TIME)
      Write(6,'(X,"TIME SPENT ON CONSTRUCTING AB(1)"
     $ ,F10.2)')END_TIME-START_TIME
      write(*,*)'icount',icount
C
      If(IFunSR.Eq.4) Then
C
C     POSTCAS: DIVIDE BY C'c AND SYMMETRIZE
C
      Do I=1,NoEig
      IP=IndBlock(1,I)
      IQ=IndBlock(2,I)
      Do J=1,NoEig
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
      EndDo
      EndDo
C
      Do I=1,NoEig
      Do J=I+1,NoEig
      ABPLUS((J-1)*NoEig+I)=
     $ Half*(ABPLUS((J-1)*NoEig+I)+ABPLUS((I-1)*NoEig+J))
      ABPLUS((I-1)*NoEig+J)=ABPLUS((J-1)*NoEig+I)
      ABMIN((J-1)*NoEig+I)=
     $ Half*(ABMIN((J-1)*NoEig+I)+ABMIN((I-1)*NoEig+J))
      ABMIN((I-1)*NoEig+J)=ABMIN((J-1)*NoEig+I)
      EndDo
      EndDo
C
      Else
C     DIVIDE BY C'c AND COPY TRIANGLE
C
      Do I=1,NoEig
      IP=IndBlock(1,I)
      IQ=IndBlock(2,I)
      Do J=1,I
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
      ABPLUS(J+(I-1)*NoEig)=ABPLUS(I+(J-1)*NoEig)
      ABMIN(J+(I-1)*NoEig)=ABMIN(I+(J-1)*NoEig)
C
      EndDo
      EndDo
C
C     end IFunSR
      EndIf
C
      Return
      End

*Deck AB0ELEMENT
      Subroutine AB0ELEMENT(ABPL,ABMIN,IR,IS,IPP,IQQ,Occ,HNO,IGFact,
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
     $ Occ(NBasis),HNO(NInte1),IGFact(NInte2),
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
      AuxTwoPQRS=One
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=Zero
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
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.1) 
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
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.1)
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
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.1) 
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
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.1) 
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

*Deck Y01GVB
      Subroutine Y01GVB(TwoNO,Occ,URe,XOne,EigY,Eig,Eig1,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2)
C
C     A ROUTINE FOR COMPUTING 0TH-ORDER Y VECTORS 
C     AND 0TH and 1ST-ORDER EIGENVALUES OF ERPA
C     FOR GVB REFERENCE 
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
     $ Eig(NDimX),Eig1(NDimX),
     $ EigY(NDimX*NDimX),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ INN(NBasis,NBasis),AMAT(NDim*NDim),BMAT(NDim*NDim)
C
      EigY(1:NDimX*NDimX)=Zero
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C 
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
      INN(IP,IQ)=IPQ
      INN(IQ,IP)=IPQ
      EndDo
      EndDo 
C
C     0-TH ORDER EIGENVECTORS AND EIGENVALUES
C
      ACAlpha=Zero
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
      Do NU=1,NDimX
C
      NUI=(NU-1)*NDimX+NU
C
      IP=IndN(1,NU)
      IQ=IndN(2,NU)
      IPQ=INN(IP,IQ)
      IPQPQ=(IPQ-1)*NDim+IPQ
C
C     0TH-ORDER EIGENVALUES
C
      Eig(NU)=Zero
C
      If(Occ(IP).Ne.Occ(IQ).And.IGem(IP).Ne.IGem(IQ)) Then
      AA=Half*((C(IP)+C(IQ))**2*ABPLUS(IPQPQ)
     $       + (C(IP)-C(IQ))**2*ABMIN(IPQPQ))
      Eig(NU)=ABS(AA/(Occ(IP)-Occ(IQ)))
      ElseIf(IGem(IP).Eq.IGem(IQ)) Then
      Eig(NU)=ABS(ABPLUS(IPQPQ))
      EndIf
C
C     0TH-ORDER Y EIGENVECTORS
C
      If(IGem(IP).Eq.IGem(IQ)) Then
C
      EigY(NUI)=SQRT(Half)
C
      Else
C
      If(Occ(IP)-Occ(IQ).Ne.Zero) Then
      If(Occ(IQ).Gt.Occ(IP))EigY(NUI)=
     $ SQRT(Half)*(C(IQ)-C(IP))/SQRT(Occ(IQ)-Occ(IP))
      If(Occ(IP).Gt.Occ(IQ))EigY(NUI)=
     $ SQRT(Half)*(C(IP)-C(IQ))/SQRT(Occ(IP)-Occ(IQ))
      Else
      EigY(NUI)=Zero
      EndIf
C
      EndIf
C
      EndDo
C
      Do I=1,NDim*NDim
      AMAT(I)=ABPLUS(I)
      BMAT(I)=ABMIN(I)
      EndDo
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
      Do NU=1,NDimX
C
      NUI=(NU-1)*NDimX+NU
C
      IP=IndN(1,NU)
      IQ=IndN(2,NU)
      IPQ=INN(IP,IQ)
      IPQPQ=(IPQ-1)*NDim+IPQ
C
      Eig1(NU)=Zero
C
      If(Occ(IP).Ne.Occ(IQ).And.IGem(IP).Ne.IGem(IQ)) Then
C
      Y=EigY(NUI)
      If(Occ(IQ).Gt.Occ(IP)) X=Y/(C(IQ)-C(IP))*(C(IQ)+C(IP))
      If(Occ(IP).Gt.Occ(IQ)) X=Y/(C(IP)-C(IQ))*(C(IQ)+C(IP))
      Eig1(NU)=X*AMAT(IPQPQ)*X+Y*BMAT(IPQPQ)*Y
C
      ElseIf(IGem(IP).Eq.IGem(IQ)) Then
C
      Eig1(NU)=AMAT(IPQPQ)*Half+BMAT(IPQPQ)*Half
C
      EndIf
C
      EndDo
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
      use abmat
      use abfofo
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
     $ ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ Eig(NDim),EGOne(NGem)
C
C     LOCAL ARRAYS
C
      Dimension IndX(NDim),IndN(2,NDim),C(NBasis),
     $ EigVY2(NBasis*(NBasis-1)),IndP(NBasis,NBasis),
     $ AMAT(NDim,NDim),BMAT(NDim,NDim)
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
C
      If(ITwoEl.Eq.1) Then
C
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
      ElseIf(ITwoEl.Eq.2) Then
C
      Call LookUp_mithap(Occ,IndAux,IndP,IndN,IndX,NDimX,NDim,NBasis)
      Call ACABMAT0_mithap(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO',0,ACAlpha,1)
C
      ElseIf(ITwoEl.Eq.3) Then
C
      Call LookUp_mithap(Occ,IndAux,IndP,IndN,IndX,NDimX,NDim,NBasis)
C
      Call ACABMAT0_FOFO(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NActive,NELE,NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO','FFOO','FOFO',0,ACAlpha,1)
C
      EndIf
C
      Do J=1,NDim
      Do I=1,NDim
      AMAT(I,J)=ABPLUS(I,J)
      BMAT(I,J)=ABMIN(I,J)
      EndDo
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C      IRSIRS=(IRS-1)*NDim+IRS
      Eig(IRS)=Zero
C
      If(Occ(IR).Ne.Occ(IS).And.IGem(IR).Ne.IGem(IS)) Then
      AA=Half*((C(IR)+C(IS))**2*ABPLUS(IRS,IRS)
     $       + (C(IR)-C(IS))**2*ABMIN(IRS,IRS))
      Eig(IRS)=ABS(AA/(Occ(IR)-Occ(IS)))
      Else
      Eig(IRS)=ABS(ABPLUS(IRS,IRS))
      EndIf
C
      EndDo
      EndDo
C
      ACAlpha=One
      If(ITwoEl.Eq.1) Then
      Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
C
      ElseIf(ITwoEl.Eq.2) Then
C
      Call ACABMAT0_mithap(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO',0,ACAlpha,1)
C
      ElseIf(ITwoEl.Eq.3) Then
C
C      TEST!
C      ACAlpha=sqrt(2d0)/2.d0
C      Print*, 'ACAlpha',ACAlpha
C
      Call ACABMAT0_FOFO(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NActive,NELE,NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO','FFOO','FOFO',0,ACAlpha,1)
C 
      EndIf
C
C     AMAT AND BMAT WILL INCLUDE 1ST-ORDER A+ AND A- MATRICES, RESPECTIVELY
C
      Do J=1,NDim
      Do I=1,NDim
      AMAT(I,J)=ABPLUS(I,J)-AMAT(I,J)
      BMAT(I,J)=ABMIN(I,J)-BMAT(I,J)
      EndDo
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
C     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-2) ) Then
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct) ) Then
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
      If(ITwoEl.Eq.1) Then
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
C      IJ=(IRS-1)*NDim+IPQ 
C      Aux=(Half*AMAT(IJ)-Two*EigVY2(IPQ)*EigVY2(IRS)*BMAT(IJ))
C     $    *EPSJI 
      Aux=(Half*AMAT(IPQ,IRS)-Two*EigVY2(IPQ)*EigVY2(IRS)*BMAT(IPQ,IRS))
     $    *EPSJI
C
C     Save Aux - it may be needed in embedding calculations
C      ABPLUS((J-1)*NDimX+I)=(C(IP)+C(IQ))*(C(IR)+C(IS))*Aux
      ABPLUS(I,J)=(C(IP)+C(IQ))*(C(IR)+C(IS))*Aux
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
      ElseIf(ITwoEl.Eq.2) Then

      Call ECorrAC0GVB_mithap(ECorr0,ECorr,AMAT,BMAT,ABPLUS,
     $                        EigVY2,Occ,C,Eig,
     $                        IndP,IndN,IndX,IGem,
     $                        'TWOMO',NDim,NDimX,NGem,NBasis)
C
      ElseIf(ITwoEl.Eq.3) Then

      Call ECorrAC0GVB_FOFO(ECorr0,ECorr,AMAT,BMAT,ABPLUS,
     $                      EigVY2,Occ,C,Eig,
     $                      IndP,IndN,IndX,IGem,
     $                      'FOFO',NActive,NELE,NDim,NDimX,NGem,NBasis)
C
C      Write(6,'(1x,a)') "SORRY!"
C      Stop
      EndIf

      ECorr=ECorr+ECorr0
C
      Write
     $ (6,'(/,1X,''0-ALPHA-ORDER CORRELATION '',2X,F15.8)') ECorr0
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
C     ADD A SR KERNEL
C
      If(IFunSR.Ne.0.And.IFunSRKer.Eq.1) Then
C
      Open(20,File="srdump",Form='UNFORMATTED')
C
      Do IRow=1,NDimX
C
      IA=IndN(1,IRow)
      IB=IndN(2,IRow)
C
      Do ICol=1,NDimX
C
      IC=IndN(1,ICol)
      ID=IndN(2,ICol)
C
      Read(20) XKer
      If(.Not.(
     $IGem(IA).Eq.IGem(IB).And.IGem(IB).Eq.IGem(IC).
     $ And.IGem(IC).Eq.IGem(ID))) XKer=XKer*ACAlpha
C
      ABMIN((ICol-1)*NDimX+IRow)=ABMIN((ICol-1)*NDimX+IRow)
     $ + XKer
C
      EndDo
      EndDo
C
      Close(20)
C
      EndIf
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(NoSt.Eq.1) Then
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Else
      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      EndIf
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
      ISkippedEig=0
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
      Write(6,'(/,X,"***** COMPUTING AMAT, BMAT IN ACABMAT0 *****",/)')
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
      XMAT=XMAT+ (-Occ(IR)*Occ(IP)*IGemPR
     $            +Occ(IS)*Occ(IP)*IGemPS
     $            +Occ(IR)*Occ(IQ)*IGemQR
     $            -Occ(IS)*Occ(IQ)*IGemQS)*AuxTwo
     $ *( Two*TwoMO(NAddr3(IQ,IP,IS,IR))-TwoMO(NAddr3(IQ,IR,IS,IP)) )
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
      XMAT=XMAT+(C(IQ)*C(IR)*IGemQR+C(IP)*C(IS)*IGemPS)*AuxTwo
     $ *( TwoMO(NAddr3(IP,IS,IQ,IR))+TwoMO(NAddr3(IP,IR,IQ,IS)) )
C
      If(IS.Eq.IP) XMAT=XMAT-C(IR)*AuxXC(IR,IQR)
      If(IR.Eq.IQ) XMAT=XMAT-C(IS)*AuxXC(IS,IPS)
      If(IR.Eq.IP) XMAT=XMAT-C(IP)*AuxXC(IP,IQS)
      If(IS.Eq.IQ) XMAT=XMAT-C(IQ)*AuxXC(IQ,IPR)
C
      If (IR.Gt.IS.And.IP.Gt.IQ) BMAT(IRS,IPQ)=XMAT
      If (IR.Gt.IS.And.IQ.Gt.IP) AMAT(IRS,IQP)=XMAT
C
      EndDo
      EndDo
      EndDo
      EndDo
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

      Subroutine LookUp_mithap(Occ,IndAux,
     $                         IndP,IndN,IndX,NDimX,NDim,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Integer :: IndAux(NBasis)
      Integer :: IndX(NDim),IndN(2,NDim)
      Dimension :: IndP(NBasis,NBasis)
      Dimension :: Occ(NBasis)
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C     
      Include 'commons.inc'
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
C     do not correlate active degenerate orbitals if from different
C     geminals
      If((IGem(I).Ne.IGem(J)).And.(IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-2) ) Then
C      Write(*,*)"Discarding nearly degenerate pair",I,J
      Else
C    
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
      end Subroutine LookUp_mithap

      Subroutine ECorrAC0GVB_mithap(ECorr0,ECorr,AMAT,BMAT,ABPLUS,
     $                 EigVY2,Occ,C,Eig,IndP,IndN,IndX,IGem,
     $                 IntFile,NDim,NDimX,NGem,NBasis)
C
      use tran
C
      Implicit None

      Integer :: NDim,NDimX,NGem,NBasis
      Integer :: IndX(NDim),IndN(2,NDim),IndP(NBasis,NBasis),
     $           IGem(NBasis)
      Character(*) :: IntFile
      Double Precision :: ECorr0,ECorr
      Double Precision :: Occ(NBasis),C(NBasis),Eig(NDim),
     $                    EigVY2(NBasis*(NBasis-1)),
     $                    AMAT(NDimX,NDimX),BMAT(NDimX,NDimX),
     $                    ABPLUS(NDimX,NDimX)
C
      Integer :: iunit
      Integer :: ip,iq,ir,is,ipq,irs
      Integer :: i,j,k,l,kl
      Integer :: pos(NBasis,NBasis)
      Logical :: AuxCoeff(NGem,NGem,NGem,NGem)
      Double Precision :: Cpq,Crs,Aux,EPSJI
      Double Precision,Allocatable :: Work(:),ints(:,:)

      pos = 0
      Do I=1,NDimX
      pos(IndN(1,I),IndN(2,I)) = IndX(I)
      EndDo
C
      AuxCoeff = .true.
      Do l=1,NGem
      Do k=1,NGem
      Do j=1,NGem
      Do i=1,NGem
      If((i==j).and.(j==k).and.(k==l)) Then
         AuxCoeff(i,j,k,l) = .false.
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C
      Allocate(work(NBasis**2),ints(NBasis,NBasis))
C     FULL INTS
      Open(newunit=iunit,file=trim(IntFile),status='OLD',
     $     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
      kl = 0
      do l=1,NBasis
         do k=1,l
            kl = kl + 1
            if(pos(l,k)/=0) then
              irs = pos(l,k)
              ir = l
              is = k
              read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
              call triang_to_sq2(work,ints,NBasis)

              Crs = (C(l)+C(k))*IndP(l,k)

              If(Occ(ir).Gt.Occ(is)) Then
              Aux=Occ(is)*(Occ(ir)-1d0)
              Else
              Aux=Occ(ir)*(Occ(is)-1d0)
              EndIf
  
C             0th-order correlation only if (IP,IQ) pair is allowed
              If(IGem(ir).Ne.IGem(is)) Then
                 ECorr0 = ECorr0 + 2d0*Aux*ints(ir,is)*IndP(ir,is)
              EndIf
 
              do j=1,NBasis
                 do i=1,j
                    if(pos(j,i)/=0) then
                      ipq = pos(j,i)
                      ip = j
                      iq = i
                      Cpq = (C(j)+C(i))*IndP(j,i)

                      EPSJI=Eig(IRS)+Eig(IPQ)
                      If(EPSJI/=0d0) Then
                         EPSJI=1d0/EPSJI 
                      Else
                         EPSJI=0d0
                      EndIf

                      Aux = (0.5d0*AMAT(irs,ipq) 
     $                    -  2d0*EigVY2(irs)*EigVY2(ipq)*BMAT(irs,ipq))
     $                    *  EPSJI

                     if(AuxCoeff(IGem(ip),IGem(iq),
     $                           IGem(ir),IGem(is))) then

                         ECorr = ECorr + Cpq*Crs*Aux*ints(ip,iq)
                      EndIf

                    endif
                 enddo
              enddo 

            endif
         enddo
       enddo

      Close(iunit)

      Deallocate(ints,work)
C
      End Subroutine ECorrAC0GVB_mithap

      Subroutine ECorrAC0GVB_FOFO(ECorr0,ECorr,AMAT,BMAT,ABPLUS,
     $                 EigVY2,Occ,C,Eig,IndP,IndN,IndX,IGem,
     $                 IntKFile,NAct,NElHlf,NDim,NDimX,NGem,NBasis)
C
      use tran
C
      Implicit None

      Integer :: NAct,NElHlf,NDim,NDimX,NGem,NBasis
      Integer :: IndX(NDim),IndN(2,NDim),IndP(NBasis,NBasis),
     $           IGem(NBasis)
      Character(*) :: IntKFile
      Double Precision :: ECorr0,ECorr
      Double Precision :: Occ(NBasis),C(NBasis),Eig(NDim),
     $                    EigVY2(NBasis*(NBasis-1)),
     $                    AMAT(NDimX,NDimX),BMAT(NDimX,NDimX),
     $                    ABPLUS(NDimX,NDimX)
C
      Integer :: iunit
      Integer :: INActive,NOccup
      Integer :: ip,iq,ir,is,ipq,irs
      Integer :: i,j,k,l,kl
      Integer :: pos(NBasis,NBasis)
      Logical :: AuxCoeff(NGem,NGem,NGem,NGem)
      Double Precision :: Cpq,Crs,Aux,EPSJI
      Double Precision,Allocatable :: Work(:),ints(:,:)

      INActive = NElHlf - NAct
      NOccup = 2*NAct + INActive

      pos = 0
      Do I=1,NDimX
      pos(IndN(1,I),IndN(2,I)) = IndX(I)
      EndDo
C
      AuxCoeff = .true.
      Do l=1,NGem
      Do k=1,NGem
      Do j=1,NGem
      Do i=1,NGem
      If((i==j).and.(j==k).and.(k==l)) Then
         AuxCoeff(i,j,k,l) = .false.
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C
      Allocate(work(NBasis**2),ints(NBasis,NBasis))
C     FOFO INTS
      Open(newunit=iunit,file=trim(IntKFile),status='OLD',
     $      access='DIRECT',recl=8*NBasis*NOccup)
      
      kl = 0
      do k=1,NOccup
         do l=1,NBasis
            kl = kl + 1
            if(pos(l,k)/=0) then
              irs = pos(l,k)
              ir = l
              is = k
              read(iunit,rec=kl) work(1:NBasis*NOccup)
              do j=1,NOccup
                 do i=1,NBasis
                    ints(i,j) = work((j-1)*NBasis+i)
                 enddo
              enddo
              ints(:,NOccup+1:NBasis) = 0

              Crs = (C(l)+C(k))*IndP(l,k)

              If(Occ(ir).Gt.Occ(is)) Then
              Aux=Occ(is)*(Occ(ir)-1d0)
              Else
              Aux=Occ(ir)*(Occ(is)-1d0)
              EndIf
  
C             0th-order correlation only if (IP,IQ) pair is allowed
              If(IGem(ir).Ne.IGem(is)) Then
                 ECorr0 = ECorr0 + 2d0*Aux*ints(ir,is)*IndP(ir,is)
              EndIf

              do j=1,NBasis
                 do i=1,j
                    if(pos(j,i)/=0) then
                      ipq = pos(j,i)
                      ip = j
                      iq = i
                      Cpq = (C(j)+C(i))*IndP(j,i)

                      EPSJI=Eig(IRS)+Eig(IPQ)
                      If(EPSJI/=0d0) Then
                         EPSJI=1d0/EPSJI 
                      Else
                         EPSJI=0d0
                      EndIf

                      Aux = (0.5d0*AMAT(irs,ipq) 
     $                    -  2d0*EigVY2(irs)*EigVY2(ipq)*BMAT(irs,ipq))
     $                    *  EPSJI

                     if(AuxCoeff(IGem(ip),IGem(iq),
     $                           IGem(ir),IGem(is))) then

                         ECorr = ECorr + Cpq*Crs*Aux*ints(ip,iq)
                      EndIf

                    endif
                 enddo
              enddo 

            endif
         enddo
      enddo
C
      Close(iunit)
C
      End Subroutine ECorrAC0GVB_FOFO

*Deck MP2RDM
      Subroutine MP2RDM(TwoNO,Eps,Occ,URe,UNOAO,XOne,
     $ IndN,IndX,IndAux,NDimX,NBasis,NDim,NInte1,NInte2,
     $ NVirt,IntFile,ThrVirtIN,iTwoNOout)
C
      use tran,only : tran4_full
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
c      Parameter(ThrVirt=3.D-5)
C
      Dimension
     $ URe(NBasis,NBasis),XOne(NInte1),UNOAO(NBasis,NBasis),
     $ Occ(NBasis),TwoNO(NInte2),
     $ IndAux(NBasis),IndX(NDim),IndN(2,NDim),
     $ Eps(NBasis,NBasis)
      Integer :: NVirt,NInte1,NInte2
      Double Precision,intent(in) :: ThrVirtIN
      Logical          :: iTwoNOout
      Character(*)     :: IntFile
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ Gamma(NInte1),
     $ PC(NBasis),AUXM(NBasis,NBasis),Work(NBasis),Fock(NBasis*NBasis),
     $ AUX2(NBasis*NBasis),work1(NBasis,NBasis),epsi(nbasis)
      CHARACTER(100) :: num1char
C
C      IF(COMMAND_ARGUMENT_COUNT().Eq.0)THEN
C      WRITE(*,*)'ERROR, COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
C      STOP
C      ENDIF
C
C      CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
C      READ(num1char,*)xnum1
C      ThrVirt=xnum1
C
      ThrVirt = ThrVirtIN
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
   10 Read(10,*,End=40)I,J,K,L,X
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
       Write(6,'(/,X,''One-Electron CASSCF Energy'',3X,F15.8)')ETot
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
C     VIRTUAL-INACTIVE
C
      Do IP=NOccup+1,NBasis
      Do IQ=1,NOccup
C
      If(IPair(IP,IQ).Eq.1) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IP,IQ,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      Eps(IP,IQ)=ABP
C
      EndIf
C     
      EndDo
      EndDo
C
      Gamma(1:NInte1)=Zero
C
      IJ=0
      Do I=1,NOccup
      Do J=1,I
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
C
      Gamma(IJ)=Zero
      If(I.Eq.J) Gamma(IJ)=Occ(I)
C
      Do K=1,NOccup
      Do IA=NOccup+1,NBasis
      Do IB=NOccup+1,NBasis
C
      Gamma(IJ)=Gamma(IJ)
     $ - ( Two*TwoNO(NAddr3(I,IA,K,IB))-TwoNO(NAddr3(I,IB,K,IA)) )
     $ *TwoNO(NAddr3(J,IA,K,IB))/
     $ (Eps(IA,I)+Eps(IB,K))/(Eps(IA,J)+Eps(IB,K))
C
      EndDo
      EndDo
      EndDo 
C
      EndDo
      EndDo
C
C
      Do IA=NOccup+1,NBasis
      Do IB=IA,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      Gamma(IAB)=Zero
C
      Do IC=NOccup+1,NBasis
      Do I=1,NOccup
      Do J=1,NOccup
C
      Gamma(IAB)=Gamma(IAB)
     $ + ( Two*TwoNO(NAddr3(I,IA,J,IC))-TwoNO(NAddr3(I,IC,J,IA)) )
     $ *TwoNO(NAddr3(I,IB,J,IC))/
     $ (Eps(IA,I)+Eps(IC,J))/(Eps(IB,I)+Eps(IC,J))
C
      EndDo
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Call CpySym(AUXM,Gamma,NBasis)
      Call Diag8(AUXM,NBasis,NBasis,PC,Work)
C
      Write(6,'(2X,"MP2",3X,"Unsorted Occupancy")')
      Sum=Zero
      Do I=NBasis,1,-1
      Write(6,'(X,I3,E16.6,I6)') I,PC(I)
      Sum=Sum+PC(I)
      EndDo
      Write(6,'(2X,"Sum of MP2 Occupancies: ",F5.2)') Sum
C
C     SET Occ of the robitals belonging to NOccup to 1 before sorting 
C     to make sure that after sorting all NOccup orbitals come first 
C
      Do I=1,NBasis
C
      IOccup=0
      IVirt=0
      Do J=1,NBasis
      If(AUXM(I,J).Ne.Zero.And.J.Le.NOccup) IOccup=1
      If(AUXM(I,J).Ne.Zero.And.J.Gt.NOccup) IVirt=1
      EndDo
C
      If(IOccup.Eq.1.And.IVirt.Eq.0) PC(I)=One
      If(IOccup*IVirt.Eq.1.Or.IOccup+IVirt.Eq.0) Stop
     $ 'MP2 1RDM is messed up. Quitting.'
C
      EndDo
C
      Call SortP(PC,AUXM,NBasis)
C
C     TRANSFORMATION MATRIX TO NEW VIRTUAL ORBITALS
C
      Do I=1,NBasis
      Do K=1,NOccup
      AUXM(K,I)=Zero
      If(I.Eq.K) AUXM(K,I)=One
      EndDo 
      EndDo
C
      NVZero=0
      Do I=NOccup+1,NBasis
      If(PC(I).Le.ThrVirt) Then
      Write(6,'(2X,"Virtual MP2 orbital no",I3," of the occup",
     $ E14.6," removed")') I,PC(I)
      Do K=1,NBasis
      AUXM(I,K)=Zero
      EndDo
      NVZero=NVZero+1
      EndIf
      EndDo 
      Write(6,'(/,2X,"Threshold for virt occup number :",E16.6)')ThrVirt
      Write(6,'(/,2X,"Total number of removed orbitals: ",I3)') NVZero
C
      NVirtOld = NBasis - NOccup
      NVirt    = NBasis - NOccup - NVZero
      Val = (1d0 - dble(NVirt)/dble(NVirtOld) ) * 100d0
      Write(6,'(2x,"NVirt/NVirtOld",18x,": ",i3," / ",i3)') 
     $ NVirt, NVirtOld
      Write(6,'(2x,"Percent of removed orbitals     : ",f7.2,/)') Val

C     If the new are orbitals are to be used in AC0 they must be 
C     cannonicalized
C
c new line
      NVirt=NBasis-NOccup
C
      Aux2=0d0
C
      Do I=1,NVirt
      Do J=1,NVirt
      II=I+NOccup
      JJ=J+NOccup
      IJ=(Max(II,JJ)*(Max(II,JJ)-1))/2+Min(II,JJ)
      FIJ=XOne(IJ)
      Do K=1,NOccup 
      FIJ=FIJ+Occ(K)*
     $ (Two*TwoNO(NAddr3(II,JJ,K,K))-TwoNO(NAddr3(II,K,JJ,K)))
      EndDo
      IJF=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      AuxI(IJF)=FIJ 
      AUX2((J-1)*NVirt+I)=AUXM(II,JJ)
      EndDo
      EndDo
C
      Call MatTr(AuxI,AUX2,NVirt)
C
      Call CpySym(Fock,AuxI,NVirt)
      Call Diag8(Fock,NVirt,NVirt,PC,Work)
C
c new line
      Call SortF(PC,Fock,NVirt) 
C     TEST
C      Print*, 'Fock,PC-Ka',norm2(Fock),norm2(PC(1:NVirt))
C
C
      Do I=1,NVirt
      Do J=1,NVirt
      II=I+NOccup
      JJ=J+NOccup
      URe(II,JJ)=Fock((J-1)*NVirt+I)
      EndDo
      EndDo
C
C     Set elements corresponding to the removed orbitals to zero
C
c new line
      NVirt=NBasis-NOccup-NVZero
C 
      Do I=NOccup+NVirt+1,NBasis
      Do J=1,NBasis
      URe(I,J)=Zero
      URe(J,I)=Zero
      EndDo
      EndDo
C
      Call MultpM(Eps,URe,AUXM,NBasis)
C
C      Print*, 'Eps-Ka',norm2(Eps)
C
C     TRANSFORM INTEGRALS
C     this step should be made more efficient - transform only virtuals!
      Write(6,'(/,2X,"Integral transformation in progress ... ",/)') 
      Call MatTr(XOne,Eps,NBasis)
C      Print*, 'XOne-Ka',norm2(XOne)
C
      If(iTwoNOout) Then
C     TRANSFORM 2-el INTEGRALS INCORE
      Call TwoNO1(TwoNO,Eps,NBasis,NInte2)
C
      Else
C     DUMP TO DISC
C     get AO --> NOMP2(canon) tran mat
      work1=transpose(UNOAO)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work1,NBasis,
     $ Eps,NBasis,0d0,AUXM,NBasis)
       Call tran4_full(NBasis,AUXM,AUXM,IntFile,'AOTWOSORT')
      EndIf
C
C     INTEGRALS ARE TRANSFORMED SO URe IS SET TO A UNIT MATRIX 
C
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=Zero
      If(I.Eq.J) URe(I,J)=One
      EndDo
      EndDo
C
C     AFTER "TRUNCATING" OF THE VIRTUAL SPACE, FIND A NEW SET OF ACCEPTED PAIRS
C
      Call AcceptPair(IndN,IndX,NDimX,IndAux,Occ,NOccup,NVirt,NBasis)
C
      Return
      End

*Deck AcceptPair
      Subroutine AcceptPair(IndN,IndX,NDimX,IndAux,Occ,
     $ NOccup,NVirt,NBasis)
C
      use types 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc' 
C
      Dimension IndX(NBasis*(NBasis-1)/2),IndN(2,NBasis*(NBasis-1)/2),
     $ IndAux(NBasis),Occ(NBasis)
C
      type(SystemBlock) :: System
C
      NOK=NOccup+NVirt
C
C      ThrSelAct=System%ThrSelAct
      Write(LOUT,'(1x,a,2e15.5)') 'Threshold for quasi-degeneracy ',
     $ ThrSelAct
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
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct)
     $ ) Then
C
      Write(6,'(1X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals  
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
C
C     EXCLUDE EXCITATIONS TO REMOVED VIRTUAL ORBITALS
C
      If(I.Le.NOK.And.J.Le.NOK) Then
C
      Ind=Ind+1
      IndX(Ind)=Ind
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
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
      Write(6,'(/1X,"AFTER TRUNCATING VIRTUAL ORBITAL SPACE: ")')
C
      NDimX=Ind
      Write(6,'(1X,"Number of pairs reduced to:",I6)')Ind
C
      If(IPrint.gt.2) Then
      Write(6,'(1X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(1X,3I5,2E14.4)')I,Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
      EndIf
C
      Return
      End


*Deck SortP
      Subroutine SortP(Occ,URe,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Occ(NBasis),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension UReOld(NBasis,NBasis),Ind(1000),IndOcc(1000)
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
C
      Return
      End

*Deck SortF
      Subroutine SortF(Pcc,URe,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Pcc(NBasis),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension UReOld(NBasis,NBasis),Ind(1000),IndPcc(1000)
C
C     SORT THE Pcc NUMBERS IN A DESCENDING ORDER OF THEIR ABS VALUES
C
      Do I=1,NBasis
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,NBasis
C
      PccMx=Pcc(IStart)
      IndMx=IStart
C
      Do J=IStart,NBasis
      If(Abs(Pcc(J)).Gt.Abs(PccMx)) Then
      PccMx=Pcc(J)
      IndMx=J
      EndIf
      EndDo
C
      Hlp=Pcc(IStart)
      IndHlp=Ind(IStart)

      Pcc(IStart)=PccMx
      Ind(IStart)=Ind(IndMx)

      Pcc(IndMx)=Hlp
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
      Return
      End

*Deck RDM_SACAS
      Subroutine RDM_SACAS(GammaS,XCAS,YCAS,EExcit,C,UNOAO,IPair,
     $ DipX,DipY,DipZ,NoState,NoStMx,NInte1,NBasis,IPr)
C
C     READS 1RDMs FOR STATES FROM 1 TO NoStMx AND TRANSFROM THEM 
C     TO THE REPRESENTATION OF NO's OF THE NoState's STATE
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0)
C
      Dimension URe(NBasis,NBasis),GammaS(100,NInte1),EExcit(NInte1)
      Dimension XCAS(NBasis,NInte1),YCAS(NBasis,NInte1),C(NBasis),
     $ IPair(NBasis,NBasis),UNOAO(NBasis,NBasis)
      Double precision,dimension(NBasis,NBasis) :: dipx,dipy,dipz
C
C     LOCAL ARRAYS
C
      Double precision :: GammaAB(NInte1),PC(NBasis)
      Double precision :: rdm(NBasis,NBasis),trdm(NBasis,NBasis),
     $                    AUXM(NBasis,NBasis),
     $                    WorkV(NBasis),AUXM1(NBasis*NBasis)
      Character*32 Str
C
      If(NoStMx.Gt.10) Stop 'Fatal error in TRDM_SACAS InSt(1,1)>10!'
C
C     grep ' !MCSCF STATE' filename.out > sacas_ene.dat
      Open(10,File="sacas_ene.dat",Status='Old')
C
      NoStMx=0
c      Do I=1,NoStMx
   20 Read(10,'(A32,F22.12)',End=40) Str,EExcit(NoStMx+1)
      If(IPr.Eq.1) 
     $ Write(6,'(X,"SA-CAS Energy for state no ",I3,4X,F12.7)')
     $ NoStMx+1,EExcit(NoStMx+1)
      Read(10,*)
c      EndDo
      NoStMx=NoStMx+1
      GoTo 20
   40 Continue
      Close(10)
C
      ERef=EExcit(NoState)
      Do I=1,NoStMx
      EExcit(I)=EExcit(I)-ERef
c      Write(6,'(X,"SA-CAS I->J De-Excit Energy ",2I3,E15.6)')
c     $ NoState,I,EExcit(I)
      EndDo
C
      NAc=NAcCAS
      NInAc=NInAcCAS
      INActive=NInAc
C
C     Prepare MO->NO_NoState
C
      Call read_1rdm_molpro(GammaAB,NoState,InSt(2,1),
     $ '2RDM',IWarn,NBasis)
C
      Call CpySym(AUXM,GammaAB,NBasis)
C
      Do I=1,NAc
      Do J=1,NAc
      AUXM1((J-1)*NAc+I)=AUXM(I,J)
      EndDo   
      EndDo
      Call Diag8(AUXM1,NAc,NAc,PC,WorkV)
C
      Call SortP(PC,AUXM1,NAc)
C
C     FULL TRANSFORMATION
C
      URe=0
      Do I=1,NBasis
      IIAct=I-NInAc
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      JJAct=J-NInAc
      If(IIAct.Gt.0.And.IIAct.Le.NAc.And.JJAct.Gt.0.And.JJAct.Le.NAc)
     $ URe(I,J)=AUXM1(IIAct+(JJAct-1)*NAc) 
      EndDo
      EndDo
C
C     1-RDMs FROM ALL STATES IN THE REPRESENTATION OF NOs OF THE NoState's STATE
C
      Do IS=1,NoStMx
C
      call read_1rdm_molpro(GammaAB,IS,InSt(2,1),'2RDM',IWarn,NBasis)
      Call CpySym(AUXM,GammaAB,NBasis)
C
      rdm=0
C
      Do I=1,NInAc
      rdm(I,I)=One
      Enddo
C
      Do J=1,NAc
      Do I=1,NAc
      rdm(INActive+I,INActive+J)=AUXM(I,J)
      Enddo
      Enddo
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,
     $           rdm,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           URe,NBasis,0d0,rdm,NBasis)
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      GammaS(IS,IAB)=rdm(IA,IB)
      EndDo
      EndDo
C
      EndDo
C
      IBra=NoState
C     IBra is the index of the state of interest.
C     Stop if IBra is greater than NBasis (dimension of XCAS, YCAS would be exceeded)
      If(IBra.Gt.NBasis) Stop 'Fatal error in TRDM_SACAS'
      If(IPr.Eq.1) Write(6,'(/," Read TRDMs FROM STATE NO",I3,/)') IBra
C
      Do IKet=1,NoStMx
      If(IKet.Ne.IBra) Then
C
      Call read_1trdm_molpro(AUXM,IBra,IKet,'2RDM',NBasis)
C
      trdm=0
C
      Do J=1,NAc
      Do I=1,NAc
      trdm(INActive+I,INActive+J) = AUXM(I,J)
      Enddo
      Enddo
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,
     $           trdm,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           URe,NBasis,0d0,trdm,NBasis)
C
      SumNU=Zero
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      YCAS(IKet,IAB)=Zero
      XCAS(IKet,IAB)=Zero
C
      If(IA.Ne.IB.And.IPair(IA,IB).Eq.1) Then
c      If(IA.Ne.IB) Then
      If(C(IA)+C(IB).Ne.Zero) YCAS(IKet,IAB)=
     $ (trdm(IA,IB)+trdm(IB,IA))/(C(IA)+C(IB))
      If(C(IA)-C(IB).Ne.Zero) XCAS(IKet,IAB)=
     $ (trdm(IB,IA)-trdm(IA,IB))/(C(IA)-C(IB))
      EndIf
C
      If(IPr.Eq.1.And.Abs(trdm(IA,IB))+Abs(trdm(IB,IA)).gt.1.d-6) Then
      Write(6,'(X,"TRDM_ab 1-TRDM_ba",2I3,2E15.6)')IA,IB,
     $ trdm(IA,IB),trdm(IB,IA)
c      Write(6,'(X,"YCAS, XCAS    ",6X,2E15.6)')
c     $ YCAS(IKet,IAB),XCAS(IKet,IAB)
      EndIf
C
      SumNU=SumNU+YCAS(IKet,IAB)*XCAS(IKet,IAB)
C
      EndDo
      EndDo
C
      SSgn=One
      If(IPr.Eq.1) 
     $ Write(6,'(X,"SumNu Y*X before normalization",E15.6)')SumNU
      If(SumNU.Lt.Zero) SSgn=-One
      If(Abs(SumNu).Gt.1.D-8) Then
      SumNU=One/Sqrt(Two*Abs(SumNU))
      Else
      SumNU=Zero
      EndIf
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      YCAS(IKet,IAB)=YCAS(IKet,IAB)
     $ *SumNU
      XCAS(IKet,IAB)=SSgn*XCAS(IKet,IAB)
     $ *SumNU
      If(IPr.Eq.1.And.Abs(YCAS(IKet,IAB))+Abs(XCAS(IKet,IAB)).Gt.1.d-6)
     $  Then 
       Write(6,'(X,"NORMALIZED Y_SA-CAS X_SA-CAS",2I3,2E15.6)')
     $ IA,IB,YCAS(IKet,IAB),XCAS(IKet,IAB)
      EndIf
      EndDo
      EndDo
C
C     TRANSITION DIPOLE MOMENTS
C
      TSDipZ=0
      TSDipY=0
      TSDipX=0
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      TSDipZ = TSDipZ + Two*trdm(J,I)*dipz(I,J)
      TSDipY = TSDipY + Two*trdm(J,I)*dipy(I,J)
      TSDipX = TSDipX + Two*trdm(J,I)*dipx(I,J)
      Enddo
      Enddo
      If(IPr.Eq.1) Write(6,
     $ '(X,I2,"->",I2,
     $ " Transition DipMoms <X> <Y> <Z> <X>^2+<Y>^2+<Z>^2 ",
     $ 4F15.8,/)')
     $ IBra,IKet,
     $ TSDipX,TSDipY,TSDipZ, TSDipX**2+TSDipY**2+TSDipZ**2
C
C     If(IKet.Ne.IBra)
      EndIf
C     IKet
      EndDo
C
      Return
      End

*Deck TrDipMoms
      Subroutine TrDipMoms(NU,TDIP2,EigY,C,IndN,DipX,DipY,DipZ,
     $ NDimX,NBasis)
C
C     FOR A GIVEN [X_tilde,Y_tilde] VECTORS TRANSITION DIPOLE MOMENTS ARE COMPUTED 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0)
C
      Dimension EigY(NDimX*NDimX),IndN(2,NDimX),C(NBasis),
     $ DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
C
      TSDipZ=0
      TSDipY=0
      TSDipX=0
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
C
      TRDM=(C(IP)+C(IQ))*EigY(I+(NU-1)*NDimX)
C
      TSDipZ = TSDipZ + Two*TRDM*DipZ(IP,IQ)
      TSDipY = TSDipY + Two*TRDM*DipY(IP,IQ)
      TSDipX = TSDipX + Two*TRDM*DipX(IP,IQ)
C
      Enddo
C
      Write(6,
     $ '(X,"Transition DipMoms <X> <Y> <Z> ",
     $ 3F15.8)') TSDipX,TSDipY,TSDipZ
C
      TDIP2=TSDipX**2+TSDipY**2+TSDipZ**2     
C 
      Return
      End

*Deck TRDM_SACAS
      Subroutine TRDM_SACAS(XCAS,YCAS,NoState,EExcit,C,IPair,
     $ NInte1,NBasis)
C
C     READS TRDM FOR DEEXCITATION FROM molpro SA-CAS (nosymmetry) 
C     AND TRANSFORMS THEM TO X,Y TILDED VECTORS NORMALIZED TO 1/2
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0)
C
      Dimension URe(NBasis,NBasis),C(NBasis),EExcit(NInte1),
     $ IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Double precision :: GammaAB(NInte1),PC(NBasis)
      Double precision :: trdm(NBasis,NBasis),AUXM(NBasis,NBasis),
     $                    WorkV(NBasis),AUXM1(NBasis*NBasis)
      Dimension XCAS(NBasis,NInte1),YCAS(NBasis,NInte1)
      Character*32 Str
C
C     grep ' !MCSCF STATE' filename.out > sacas_ene.dat
      Write(6,'()')
      Open(10,File="sacas_ene.dat",Status='Old')
      NoState=InSt(1,1)
      Do I=1,NoState
      Read(10,'(A32,F22.12)') Str,EExcit(I)
      Write(6,'(X,"SA-CAS Energy for state no ",I3,4X,F12.7)')
     $ I,EExcit(I)
      Read(10,*) 
      EndDo  
      Close(10)
      Do I=1,NoState-1
      EExcit(I)=EExcit(I)-EExcit(NoState)
      Write(6,'(X,"SA-CAS I->J De-Excit Energy ",2I3,E15.6)') 
     $ InSt(1,1),I,EExcit(I)
      EndDo
C
      NAc=NAcCAS
      NInAc=NInAcCAS
      INActive=NInAc
C
C     Prepare MO->NO
C
      Call read_1rdm_molpro(GammaAB,InSt(1,1),InSt(2,1),
     $ '2RDM',IWarn,NBasis)
C
      Call CpySym(AUXM,GammaAB,NBasis)
      Do I=1,NAc
      Do J=1,NAc
      AUXM1((J-1)*NAc+I)=AUXM(I,J)
      EndDo
      EndDo
      Call Diag8(AUXM1,NAc,NAc,PC,WorkV)
      Call SortP(PC,AUXM1,NAc)

c      Call Diag8(AUXM,NBasis,NBasis,PC,WorkV)
c      Call SortP(PC,AUXM,NBasis)
C
C
C     FULL TRANSFORMATION
      URe=0
      Do I=1,NBasis
      IIAct=I-NInAc
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      JJAct=J-NInAc
      If(IIAct.Gt.0.And.IIAct.Le.NAc.And.JJAct.Gt.0.And.JJAct.Le.NAc)
c     $ URe(I,J)=AUXM(IIAct,JJAct)
     $ URe(I,J)=AUXM1(IIAct+(JJAct-1)*NAc)
      EndDo
      EndDo
C
      IBra=InSt(1,1)
C     IBra is the index of the excited state of interest. 
C     Stop if IBra is greater than NBasis (dimension of XCAS, YCAS would be exceeded)
      If(IBra.Gt.NBasis) Stop 'Fatal error in TRDM_SACAS'
      Write(6,'(/," Read TRDMs TO STATES LOWER THAN STATE NO",I3,/)')
     $ IBra
C
      Do IKet=1,IBra-1
C
      Call read_1trdm_molpro(AUXM,IBra,IKet,
     $ '2RDM',NBasis)
C
      trdm=0
C
      Do J=1,NAc
      Do I=1,NAc
      trdm(INActive+I,INActive+J) = AUXM(I,J)
      Enddo
      Enddo
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,
     $           trdm,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           URe,NBasis,0d0,trdm,NBasis)
C
      SumNU=Zero
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      YCAS(IKet,IAB)=Zero
      XCAS(IKet,IAB)=Zero
C
      If(IA.Ne.IB.And.IPair(IA,IB).Eq.1) Then
c      If(IA.Ne.IB) Then
      If(C(IA)+C(IB).Ne.Zero) YCAS(IKet,IAB)=
     $ (trdm(IA,IB)+trdm(IB,IA))/(C(IA)+C(IB))
      If(C(IA)-C(IB).Ne.Zero) XCAS(IKet,IAB)=
     $ (trdm(IB,IA)-trdm(IA,IB))/(C(IA)-C(IB))
      EndIf
C
      If(Abs(trdm(IA,IB))+Abs(trdm(IB,IA)).gt.1.d-6) Then
      Write(6,'(X,"1-TRDM_(ab,ba)",2I3,2E15.6)')IA,IB,
     $ trdm(IA,IB),trdm(IB,IA)
      Write(6,'(X,"YCAS, XCAS    ",6X,2E15.6)')
     $ YCAS(IKet,IAB),XCAS(IKet,IAB)
      EndIf
C
      SumNU=SumNU+YCAS(IKet,IAB)*XCAS(IKet,IAB)
C
      EndDo
      EndDo
C
      SSgn=One
      Write(6,'(X,"SumNu Y*X",E15.6,/)')SumNU
      If(SumNU.Lt.Zero) SSgn=-One
      SumNU=One/Sqrt(Two*Abs(SumNU))
      Do I=1,NInte1
      YCAS(IKet,I)=YCAS(IKet,I)
c    *SumNU
      XCAS(IKet,I)=SSgn*XCAS(IKet,I)
C     *SumNU
      EndDo
C
C     IKet
      EndDo
C
      Return
      End

*Deck DEEXCIT
      Subroutine DEEXCIT(TrGamm,EExcit,Occ,UNOAO,XOne,TwoEl,ENuc,
     $ NGem,NInte1,NInte2,NBasis)
C
C     THIS PROCEDURE HAS BEEN USED MAINLY FOR DIFFERENT TESTS
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)

c
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ UNOAO(NBasis,NBasis),XOne(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Double precision,dimension(NBasis,NBasis) :: dipx,dipy,dipz
      Double precision :: GammaAB(NInte1),PC(NBasis)
      Double precision :: trdm(NBasis,NBasis),AUXM(NBasis,NBasis),
     $                    WorkV(NBasis)
      Dimension TrGamm(NInte1,NInte1)
C
c for testing
c     local
      Dimension GammAl(NInte1),EExcit(NInte1),
     $ H1Alph(NInte1),OccAlph(NBasis),UReAlph(Nbasis,NBasis)
C remove after tests
      Dimension TwoElA(NInte2)
C
      Do I=1,NInte2
      TwoElA(I)=TwoEl(I)
      EndDo
C test
C
      NoEig=InSt(1,1)
C     do not call anything else after OptTwoATrip: TwoEl is modified
c      Call OptTwoATrip(TrGamm,GammAl,
c     $ EExcit,ETotAlph,ENuc,Occ,XOne,TwoElA,
c     $ H1Alph,UReAlph,OccAlph,
c     $ NBasis,NInte1,NInte2,NoEig,0.D0)
c      Do I=1,15
c      write(*,*)i,'excit ene',eexcit(i)
c      EndDo
c      Return
      NDim=NBasis*(NBasis-1)/2
      Call ACPINO(ENuc,TwoEl,Occ,XOne,
     $ NBasis,NInte1,NInte2,NDim,NGem,NoEig)
      stop
C
C     Prepare MO->NO
C
      Call read_1rdm_molpro(GammaAB,InSt(1,1),InSt(2,1),
     $ '2RDM',IWarn,NBasis)
C
      Call CpySym(AUXM,GammaAB,NBasis)
      Call Diag8(AUXM,NBasis,NBasis,PC,WorkV)
      Call SortP(PC,AUXM,NBasis)
C
      Sum=Zero
      NAc=0
      Do I=1,NBasis
      Sum=Sum+PC(I)
      If(PC(I).Gt.Zero) NAc=NAc+1
      EndDo
C
      NInAc=XELE-Sum+1.D-2
      INActive=NInAcCAS
C
C     FULL TRANSFORMATION
      URe=0
      Do I=1,NBasis
      IIAct=I-NInAc
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      JJAct=J-NInAc
      If(IIAct.Gt.0.And.IIAct.Le.NAc.And.JJAct.Gt.0.And.JJAct.Le.NAc)
     $ URe(I,J)=AUXM(IIAct,JJAct)
      EndDo
      EndDo
C
      IBra=InSt(1,1)
      write(*,*)'BRA',ibra
C
      Do IKet=1,IBra-1
      write(*,*)'KET',iket
C
      Call read_1trdm_molpro(AUXM,IBra,IKet,
     $ '2RDM',NBasis)
C
      trdm=0
C
      Do J=1,NAc
      Do I=1,NAc
      trdm(INActive+I,INActive+J) = AUXM(I,J)
      Enddo
      Enddo
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,
     $           trdm,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           URe,NBasis,0d0,trdm,NBasis)
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      If(IA.Ne.IB) Then
      TrGamm(IAB,IKet)=trdm(IA,IB)+trdm(IB,IA)
      Else
      If(IA.Eq.IB) TrGamm(IAB,IKet)=trdm(IA,IB)
      EndIf
C
      if(abs(TrGamm(IAB,iket)).gt.1.d-8)
     $ write(*,*)iket,ia,ib,TrGamm(IAB,IKet)
C
      EndDo
      EndDo
C
C     ***************************************************************************
C     This part is only for checking transition dipole moments with molpro output
C     READ THE DIPOLE MOMENT
c      Call read_dip_molpro(dipx,dipy,dipz,'DIP',NBasis)
cC
c      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
c     $           dipz,NBasis,0d0,AUXM,NBasis)
c      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
c     $           UNOAO,NBasis,0d0,dipz,NBasis)
cC
c      TSDipZ=0
c      IJ=0
c      Do I=1,NBasis
c      Do J=1,NBasis 
c      TSDipZ = TSDipZ + Two*trdm(J,I)*dipz(I,J)
c      Enddo
c      Enddo
c      Write(6,
c     $ '(X,I2,I2,"  Transition State DMZ: ",F15.8)')IBra,IKet,TSDipZ
cC
C     ***************************************************************************
C     do iket
      EndDo
C
      Return
      End

*Deck SortEigXY
      Subroutine SortEigXY(IFlag,Eig,EigY,EigX,N)
C
C     SORT Eig IN AN ASCENDING ORDER OF THE MODULUS OF EIG
C     AND CHANGE THE ORDER OF THE IMAGINARY Eig AND EigVec
C
C     IFlag = 1 - SORT ACC. TO EIG; EIGENVECTORS STORED IN COLUMNS OF EigVec
C           = 0 - SORT ACC. TO EIG; EIGENVECTORS STORED IN ROWS OF EigVec
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Eig(N),EigY(N*N),EigX(N*N)
C
C     LOCAL ARRAYS
C
      Dimension EigVecOldY(N,N),EigVecOldX(N,N),Save(N),Ind(N)
C
      Do I=1,N
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,N
C
c      EMin=Abs(Eig(IStart))
c      If(IFlag.Eq.0) 
      EMin=Eig(IStart)
      IndMin=IStart
C
c      Do J=IStart,N
c      If(Abs(Eig(J)).Lt.EMin) Then
c      EMin=Abs(Eig(J))
c      IndMin=J
c      EndIf
c      EndDo
C
C
c      If(IFlag.Eq.0) Then
      Do J=IStart,N
      If(Eig(J).Lt.EMin) Then
      EMin=Eig(J)
      IndMin=J
      EndIf
      EndDo
c      EndIf
C
      Hlp=Eig(IStart)
      IndHlp=Ind(IStart)

      Eig(IStart)=Eig(IndMin)
      Ind(IStart)=Ind(IndMin)

      Eig(IndMin)=Hlp
      Ind(IndMin)=IndHlp
C
      IStart=IStart+1
C
      EndDo
C
C     SWAP THE EIGENVECTORS
C
      Do I=1,N
      Do J=1,N
      IJ=(J-1)*N+I
      EigVecOldY(I,J)=EigY(IJ)
      EigVecOldX(I,J)=EigX(IJ)
      EndDo
      EndDo
C
      Do I=1,N
C
      If(IFlag.Eq.1) Then
      Do J=1,N
      IJ=(J-1)*N+I
      EigY(IJ)=EigVecOldY(I,Ind(J))
      EigX(IJ)=EigVecOldX(I,Ind(J))
      EndDo
      Else
      Do J=1,N
      IJ=(J-1)*N+I
      EigY(IJ)=EigVecOldY(Ind(I),J)
      EigX(IJ)=EigVecOldX(Ind(I),J)
      EndDo
      EndIf
C
      EndDo
C
      Return
      End

*Deck ReadDip
      Subroutine ReadDip(DipX,DipY,DipZ,UNOAO,NBasis)
C
C     Read dipole moment matrices and transform them to NO
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0)
C
      Dimension DipX(NBasis,NBasis),DipY(NBasis,NBasis),
     $ DipZ(NBasis,NBasis),UNOAO(NBasis,NBasis),AUXM(NBasis,NBasis)
C
      Call read_dip_molpro(DipX,DipY,DipZ,'DIP',NBasis)
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
     $           dipz,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           UNOAO,NBasis,0d0,dipz,NBasis)
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
     $           dipy,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           UNOAO,NBasis,0d0,dipy,NBasis)
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
     $           dipx,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           UNOAO,NBasis,0d0,dipx,NBasis)
C
      Return
      End

*Deck AC0DSYMM
      Subroutine AC0DSYMM(ICAS,NoStMx,ICORR,EExcit,IStCAS,NSym,NSymNO,
     $ MultpC,ECorrSym,
     $ ETot,TwoNO,Occ,URe,XOne,
     $ UNOAO,IndN,IndX,NBasis,NDimX,NInte1,NInte2) 
C
C     AC0 AND DEEXCITATION CORRECTIONS BASED ON SYMMETRY
C
C      use sapt_ener
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
     $ IndX(NDimX),IndN(2,NDimX),
     $ UNOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Integer, Allocatable :: Ind(:)
      Real*8, Allocatable :: AuxY(:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ Eig(NDimX),EigY(NDimX*NDimX),
     $ EigX(NDimX*NDimX),
     $ IEigAddY(2,NDimX),IEigAddInd(2,NDimX),IndBlock(2,NDimX),
     $ XMAux(NDimX*NDimX),work1(NBasis,NBasis),
     $ TrGamm(NInte1,NInte1),
     $ GammaS(100,NInte1)
     $ ,IMatch(NDimX),EigY1(NDimX*NDimX),
     $ Eig1(NDimX),
     $ DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis),
     $ ECorrSym(100),EExcit(100),NSymNO(NBasis),MultpC(8,8),
     $ IStCAS(2,100),ICORR(100),IStERPA(2,100)
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
      If(IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) 
     $ Call ReadDip(DipX,DipY,DipZ,UNOAO,NBasis)
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
   10 Read(10,*,End=40)I,J,K,L,X
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
      Write(6,'(/," *** ACTIVE-ACTIVE BLOCK ***")')
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
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABMIN((IRow-1)*NDimB+ICol)=ABM
      Else
      ABPLUS((ICol-1)*NDimB+IRow)=Zero
      ABPLUS((IRow-1)*NDimB+ICol)=Zero
      ABMIN((ICol-1)*NDimB+IRow)=Zero
      ABMIN((IRow-1)*NDimB+ICol)=Zero
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
C
c      Print*, 'ACT-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
C
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
C
      EndIf
C
      EndIf
C
      Write(6,'(X,"Active ERPA Eigenvalues")')
C
C     set negative eigs (dexcit correction will be added so there is no need to count it twice)
C
      Do NU=1,NDimB
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
      EndDo
C
C     find symmetry of excitation
C
      Do NU=1,NDimB
C
      IStERPA(2,NU)=0
C
      Do I=1,NDimB
C
      If(Abs(EigY(NFree2-1+(NU-1)*NDimB+I))+
     $ Abs(EigY(NFree2-1+(NU-1)*NDimB+I)).Gt.0.1) Then
C
      IP=IndBlock(1,NFree1-1+I)
      IQ=IndBlock(2,NFree1-1+I)
      ISym=MultpC(NSymNO(IP),NSymNO(IQ)) 
C
      If(IStERPA(2,NU).Eq.0) Then
      IStERPA(2,NU)=ISym
      Else
      If(IStERPA(2,NU).Ne.ISym)
     $ Write(6,'("In AC0DSYMM: Symm of act-act excit ",I3,
     $ " cannot be established")')NU
      EndIf
C
C     If(Abs(EigY...
      EndIf
C
      EndDo      
C
      EndDo
C
C     find the order within each irrep     
C
      Allocate(Ind(1:NDimB))
      Do I=1,NDimB
      Ind(I)=I
      EndDo

      Do ISym=1,NSym
C
      Eig1(1:NDimB)=Eig(1:NDimB)
      IStart=1
      ICount=1
C
      Do NU=1,NDimB
      If(IStERPA(2,NU).Eq.ISym) Then
C
      EigMin=Eig1(NU)
      IndMin=IStart
C
      Do MU=IStart,NDimB
      If(IStERPA(2,MU).Eq.ISym.And.Eig1(MU).Lt.EigMin) Then
      EigMin=Eig1(MU)
      IndMin=MU
      EndIf 
      EndDo
C
      Hlp=Eig1(IStart)
      IndHlp=Ind(IStart)

      Eig1(IStart)=EigMin
      Ind(IStart)=Ind(IndMin)

      IStERPA(1,Ind(IndMin))=ICount
      Eig1(IndMin)=Hlp
      Ind(IndMin)=IndHlp
C
      ICount=ICount+1
C
      EndIf
C
      IStart=IStart+1
C     Do NU=1
      EndDo
C     Do ISym=1
      EndDo
C
C     Shift labels in each irrep depending on the number of states in SA of the energy 
C     lower or equal than that of ICAS
C
      Do I=1,NoStMx
      If(ICORR(I).Eq.0.Or.I.Eq.ICAS) Then
C
      ISym=IStCAS(2,I)
      Do NU=1,NDimB
      If(IStERPA(2,NU).Eq.ISym) IStERPA(1,NU)=IStERPA(1,NU)+1
      EndDo 
C
      EndIf 
      EndDo
C
      Do NU=1,NDimB
      Write(6,'(X,I2,2X,I2,".",I1,E15.6)') NU,IStERPA(1,NU),
     $ IStERPA(2,NU),Eig(NFree1-1+NU)
      EndDo
c hererxxx
c      Do I=1,NoStMx
c      Write(6,'(X,"Excit Energy ",I1,".",I1,2F22.12)')
c     $ IStCAS(1,I),IStCAS(2,I),EExcit(I),EExcit(I)-EExcit(ICAS)
c      EndDo
c      stop
C
      NoEigAct=NDimB
C
      Deallocate(Ind)
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS
C
      Write(6,'(/," *** ACTIVE-INACTIVE BLOCKS ***")')
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
c      Print*, 'AI-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN)
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
C
      EndDo
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
C      Print*, 'AV-KA',norm2(ABPLUS(1:NDimB**2)),norm2(ABMIN(1:NDimB**2))
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      Do NU=1,NDimB
C
      If(Eig(NFree1-1+NU).Lt.Zero) Then
      Write(6,'(X,"Setting to Zero Negative Eigs",E15.6)')
     $ Eig(NFree1-1+NU)
      Eig(NFree1-1+NU)=Zero
      Do I=1,NDimB
      EigY(NFree2-1+(NU-1)*NDimB+I)=Zero
      EigX(NFree2-1+(NU-1)*NDimB+I)=Zero
      EndDo
      EndIf
C
      EndDo
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
      Print*, 'NoEig,NDimX',NoEig,NDimX
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
C ----------------------------------------------------------------
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
      Do IDCORR=1,NoStMx
      If (ICORR(IDCORR).Eq.1) Then
C
C     COMPUTE THE AC0 or AC0D CORRELATION ENERGY 
C
      XMAux(1:NoEig*NoEig)=Zero
C
      IERPA=0
      Do NU=1,NoEigAct
      If(IStERPA(1,NU).Eq.IStCAS(1,IDCORR).And.
     $   IStERPA(2,NU).Eq.IStCAS(2,IDCORR)) IERPA=NU
      EndDo       
C
      IAC0=0
      If(IStCAS(1,ICAS).Eq.IStCAS(1,IDCORR).And.
     $ IStCAS(2,ICAS).Eq.IStCAS(2,IDCORR)) IAC0=1
C
      If(IAC0.Eq.0.And.IERPA.Ne.0)
     $ Write(6,'(X,
     $ "Deexcitation correction is computed for ERPA vector no",
     $ I2," Sym=",I1,".",I1," Eig=",F15.8)')  IERPA,
     $ IStERPA(1,IERPA),IStERPA(2,IERPA),Eig(IERPA)
      If(IAC0.Eq.0.And.IERPA.Eq.0)
     $ Write(6,'(/," ERPA vector for deexcitation correction could not
     $ be determined. The correction will be set to 0.")')
      If(IAC0.Eq.1) Write(6,'(X,
     $ "AC0 correction is computed for SA-CAS state",I2," Sym=",
     $ I1,".",I1)')  ICAS, IStCAS(1,ICAS),IStCAS(2,ICAS)
C
      Do MU=1,NoEig
      If(Eig(MU).Ne.Zero) Then
C
      Do NU=1,NoEig
      If(Eig(NU).Ne.Zero) Then
C
      If(IAC0.Eq.1.Or.MU.Eq.IERPA.Or.NU.Eq.IERPA) Then 
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
C     If(Eig(NU).Ne.Zero)
      EndIf
      EndDo
C
C     If(Eig(MU).Ne.Zero)
      EndIf
      EndDo
C
C     COMPUTE 1st-ORDER Y AND X PART 1 ONLY FOR THE 1.1 STATE
C
      If(IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) Then
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
C
      If((MU.Ne.NU).And.(Abs(Eig(MU)-Eig(NU)).Gt.1.D-12)) Aux2=
     $ (ABPLUS(MU+(NU-1)*NoEig)+ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)-Eig(NU))
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
C     If(IStCAS(1,ICAS).Eq.1 ....
      EndIf
C
C     CONTINUATION OF AC0
C
      EigX(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C
      IStart=IEigAddY(1,MU)
C
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
C
      Do J=1,NoEig
      EigX(I+(J-1)*NoEig)=EigX(I+(J-1)*NoEig)+XMAux(MU+(J-1)*NoEig)
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
      SumY=EigX(I+(J-1)*NoEig)
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
      ECorr=EAll-EIntra
      ECorrSym(IDCORR)=ECorr
C
C     TRANSITION DIPOLE MOMENTS ARE COMPUTED ONLY FOR THE 1.1 STATE
C
      If(IStCAS(1,ICAS).Eq.1.And.IStCAS(2,ICAS).Eq.1) Then
C
      Allocate (AuxY(NoEig*NoEig)) 
C
C     Y(0) AND X(0) IN COLUMNS
C
      Do I=1,NoEig
      Do J=1,NoEig
      EigX(I+(J-1)*NoEig)=Zero
      If(I.Eq.J) EigX(I+(J-1)*NoEig)=One
      EndDo
      EndDo
C
      Do I=1,NFree2
      XMAux(I)=EigY(I)
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
c      EigY(NU+(MU-1)*NoEig)=Zero
      AuxY(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      AuxY(NU+(MU-1)*NoEig)=AuxY(NU+(MU-1)*NoEig) 
     $ +EigX(NU+(I-1)*NoEig)*XMAux(IStart+II)
c      EigY(NU+(MU-1)*NoEig)=EigY(NU+(MU-1)*NoEig)
c     $ +EigX(NU+(I-1)*NoEig)*XMAux(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
C     SORT Y0,X0 AND Y1,X1 ACCORDING TO IndN
C
C     Check if NoEig=NDimX - they should be equal!
      If(NoEig.Ne.NDimX) Stop 'Fatal error in AC0DSYMM: NoEig.Ne.NDimX!'
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
      Call CpyM(EigX,AuxY,NDimX) 
c      Call CpyM(EigX,EigY,NDimX)
      Call CpyM(XMAux,EigY1,NDimX)
C
      Do MU=1,NDimX
      Do I=1,NDimX
      AuxY((MU-1)*NoEig+I)=EigX((MU-1)*NoEig+IMatch(I))
c      EigY((MU-1)*NoEig+I)=EigX((MU-1)*NoEig+IMatch(I))
      EigY1((MU-1)*NoEig+I)=XMAux((MU-1)*NoEig+IMatch(I))
      EndDo
      EndDo
C
      If(IERPA.Ne.0) Then
C
      Write(6,'(/,X,"Y(0) corresponding to vector no",I2,
     $ " Sym=",I1,".",I1)')IERPA,IStERPA(1,IERPA),IStERPA(2,IERPA)
      Call TrDipMoms(IERPA,TDIP2,AuxY,
c      Call TrDipMoms(IERPA,TDIP2,EigY,
     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,"1.1->",I1,".",
     $ I1,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA),
     $ IStERPA(2,IERPA),TDIP2
C
      Write(6,'(/,X,"Y(0)+Y(1) corresponding to vector no",I2,
     $ " Sym=",I1,".",I1)')IERPA,IStERPA(1,IERPA),IStERPA(2,IERPA)
C
      Do I=1,NdimX
      XMAux((IERPA-1)*NoEig+I)=AuxY((IERPA-1)*NoEig+I)
c      XMAux((IERPA-1)*NoEig+I)=EigY((IERPA-1)*NoEig+I)
     $ +EigY1((IERPA-1)*NoEig+I)
      EndDo
      Call TrDipMoms(IERPA,TDIP2,XMAux,
     $ C,IndN,DipX,DipY,DipZ,NDimX,NBasis)
      Write(6,'(X,"1.1->",I1,".",
     $ I1,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA),
     $ IStERPA(2,IERPA),TDIP2
C
      Else
C
      If(IAC0.Ne.1) Then 
      Write(6,'(X,"1.1->",I1,".",
     $ I1,"  Y(0) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA),
     $ IStERPA(2,IERPA),Zero
      Write(6,'(X,"1.1->",I1,".",
     $ I1,"  Y(0)+Y(1) <X>^2+<Y>^2+<Z>^2",F15.8,/)')IStERPA(1,IERPA),
     $ IStERPA(2,IERPA),Zero
      EndIf
C
C     If( IERPA.Ne.0) Then
      EndIf
C
      Deallocate(AuxY)
C
C     If(IStCAS(1,ICAS).Eq.1 ....
      EndIf
C
C     If (ICORR(IDCORR).Eq.1)...
      EndIf
C     Do IDCORR ....
      EndDo
C
C end of AC0DSYMM
      Return
      End

*Deck ACECORR
      Subroutine DelInts(ITwoEl)
C
C     DELETE MO INTEGRALS 
C
      Implicit Real*8 (A-H,O-Z)

      If(ITwoel.Eq.2) Then
        Open(newunit=iunit,file='TWOMO',status='OLD')
        Close(iunit,status='DELETE')
      ElseIf(ITwoel.Eq.3) Then
        Open(newunit=iunit,file='FFOO',status='OLD')
        Close(iunit,status='DELETE')
        Open(newunit=iunit,file='FOFO',status='OLD')
        Close(iunit,status='DELETE')
      EndIf

      Return
      End

