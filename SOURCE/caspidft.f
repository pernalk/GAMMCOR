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
c      Real*8, Dimension(:), Allocatable :: OrbGrid(:,:)
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Dimension(:), Allocatable :: Zk
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:)
      Real*8, Allocatable :: OnTop(:)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ ConCorr(200)
C
      Write(6,'(/,X,"***************** CASPIDFT ***************** ")')
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NGrid,NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (Sigma(NGrid))
      Allocate  (Zk(NGrid))
      Allocate  (RR(3,NGrid))
      Allocate  (OnTop(NGrid))
C
      Call molprogrid1(RR,NGrid)
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
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo 
C
      EndDo
C
      Call LYP(RhoGrid,Sigma,Zk,NGrid)
c      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
C
      ConCorr(1:200)=Zero
      A=0.35D0
      B=A-One
      C=2.6D0
      G=1.5D0
      D=(C-One)/(One-G)**2
C
      EDYN=Zero
      ELYP=Zero
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C 
      XX=Two*OnTop(I)/RhoGrid(I)**2
c      If(Abs(RR(1,I)).Lt.1.D-8.And.Abs(RR(2,I)).Lt.1.D-8)
      If(Abs(RR(1,I)).eq.zero.And.Abs(RR(2,I)).eq.zero)
     $ Write(*,'(7E15.6)')RR(1,I),RR(2,I),RR(3,I),RhoGrid(I),
     $ XX,Zk(I),WGrid(I)
c      If(xx.gt.1.01.and.abs(Zk(I)*WGrid(I)).gt.1.d-5)
c      If(xx.gt.0.98.and.abs(xx-1.d0).gt.1.d-4.and.
c     $ abs(Zk(I)*WGrid(I)).gt.1.d-5)
c     $ Write(*,'(I6,X,6E15.6)')I,RR(1,I),RR(2,I),RR(3,I),RhoGrid(I),
c     $ XX,Zk(I)*WGrid(I)
c     $ Write(*,'(I4,X,5E15.6)')I,RR(1,I),RR(2,I),RR(3,I),RhoGrid(I),XX
C
c herer!!!
      If(XX.Le.One) Then      
      PX=A*XX/(One+B*XX)
      Else
      PX=C*XX**0.25-D*(XX-G)**2
c      PX=One
      EndIf
c herer!!!
      if(xx.lt.0.75.and.xx.gt.0.60) then
      px=0.8
      endif
C
      IXX=Int(XX*100)
      If(IXX.Eq.0) IXX=1
      ConCorr(IXX)=ConCorr(IXX)+Zk(I)*WGrid(I)
C
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
      ELYP=ELYP+Zk(I)*WGrid(I)
c      if(abs(Zk(I)*WGrid(I)).gt.1.d-4) 
c       write(*,'(i4,x,6e15.4)')i,RhoGrid(I),xx,PX*Zk(I)*WGrid(I),
c     $ Zk(I)*WGrid(I),EDYN,ELYP
c      if(XX.lt.0.9.And.XX.Gt.0.5.and.abs(Zk(I)*WGrid(I)).gt.1.d-4)
c     $ write(*,'(i4,x,6e15.4)')i,RhoGrid(I),xx,PX*Zk(I)*WGrid(I),
c     $ Zk(I)*WGrid(I),EDYN,ELYP
C
      EndIf
C
      EndDo
C
      Sum=Zero
      Sum1=Zero
      Sum07=Zero
      Do I=1,120
      Write(*,'(i4,x,e17.6)'),I,ConCorr(I)
      Sum=Sum+ConCorr(I)
      If(I.Gt.100) Sum1=Sum1+ConCorr(I)
      If(I.Gt.59.And.I.Lt.91) Sum07=Sum07+ConCorr(I)
      EndDo
      Write(*,'("Sum:",x,e17.6)')Sum
      Write(*,'("Contribution from X>1 region: ",x,e17.6)')Sum1
      Write(*,'("Contribution from 0.90.0>X>0.60 region: ",x,e17.6)')
     $ Sum07
C
      Write(6,'(/,1X,''LYP Correlation'',7X,F15.8)')ELYP
      Write(6,'(1X,''CASPIDFT Correlation'',2X,F15.8)')EDYN
      Write(6,'(1X,''Total CAS+ LYP Energy'',X,F15.8)')ELYP+ETot+ENuc
      Write(6,'(1X,''Total CASPIDFT Energy'',X,F15.8)')EDYN+ETot+ENuc
C
      Stop 
C
      Return
      End
