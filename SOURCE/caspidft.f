*Deck CASPIDFT
      Subroutine CASPIDFT(ENuc,URe,UNOAO,Occ,XOne,TwoNO,
     $ NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Dimension(:), Allocatable :: OrbGrid(:,:)
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Dimension(:), Allocatable :: Zk 
      Real*8, Allocatable :: RDM2Act(:)
      Real*8, Allocatable :: OnTop(:)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),XOne(NInte1),TwoNO(NInte2)
C
      Write(6,'(/,X," *** CASPIDFT *** ")')
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,X," The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NBasis,NGrid))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (Sigma(NGrid))
      Allocate  (Zk(NGrid))
      Allocate  (OnTop(NGrid))
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
      ETot=ETot+Two*Occ(I)*XOne(II)
      EndDo
      Write(6,'(/,1X,''One-Electron Energy'',6X,F15.8)')ETot 
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
      Write(6,'(1X,''CASSCF Energy (w/o ENuc)'',X,F15.8)')ETot
      Write(6,'(1X,''Total CASSCF Energy '',5X,F15.8)')ETot+ENuc
C      
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)   
C
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop(I)=OnTop(I)
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(IP,I)*OrbGrid(IQ,I)*OrbGrid(IR,I)*OrbGrid(IS,I)
      EndDo
      EndDo
      EndDo
      EndDo 
C
      EndDo
C
      Call LYP(RhoGrid,Sigma,Zk,NGrid)
C
      A=0.35D0
      B=A-One
C
      EDYN=Zero
      ELYP=Zero
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C 
      XX=Two*OnTop(I)/RhoGrid(I)**2
      PX=A*XX/(One+B*XX)
C
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
      ELYP=ELYP+Zk(I)*WGrid(I)
C
      EndIf
C
      EndDo
C
      Write(6,'(/,1X,''LYP Correlation'',6X,F15.8)')ELYP
      Write(6,'(1X,''CASPIDFT Correlation'',X,F15.8)')EDYN
      stop 
C
      Return
      End
