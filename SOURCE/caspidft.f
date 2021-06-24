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
      Real*8, Dimension(:), Allocatable :: Zk,Zk1
      Real*8, Dimension(:), Allocatable :: rhoo,sigmaco,sigmaoo,vrhoc,
     $ vrhoo,vsigmacc,vsigmaco,vsigmaoo
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:)
      Real*8, Allocatable :: OnTop(:)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ ConCorr(200)
     $ ,EpsC(NBasis,NBasis),EpsCI(NBasis),EpsDiag(NBasis)
      logical fderiv,open
c      double precision rhoo(ngrid)
c      double precision sigmaco(ngrid),sigmaoo(ngrid)
      integer igrad
      character*(30) name
c      double precision vrhoc(ngrid),vrhoo(ngrid)
c      double precision vsigmacc(ngrid),vsigmaco(ngrid),vsigmaoo(ngrid)
C
c      Call CORRELON(URe,UNOAO,Occ,NBasis)
c      Call DELCORREL(URe,UNOAO,Occ,NBasis)
c      stop
C
      FDeriv=.True.
      Open=.False.
      Alpha=0.4
      CMix=0.5
c      Write(*,'(/,X,"ALPHA CMix",X,2F5.2)')Alpha,CMix
C
      Write(6,'(/,X,"***************** CASPIDFT ***************** ")')
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (rhoo(NGrid))
      Allocate  (sigmaco(NGrid))
      Allocate  (sigmaoo(NGrid))
      Allocate  (vrhoc(NGrid))
      Allocate  (vrhoo(NGrid))
      Allocate  (vsigmacc(NGrid)) 
      Allocate  (vsigmaco(NGrid)) 
      Allocate  (vsigmaoo(NGrid))
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NGrid,NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (Sigma(NGrid))
      Allocate  (Zk(NGrid))
      Allocate  (Zk1(NGrid))
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
      vrhoc(I)=Zero
      vsigmacc(i)=Zero 
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
C     OVERLAP OF ORBITALS
C
      Do I=1,NOccup
      If(Occ(I).Lt.0.8.And.Occ(I).Gt.0.2) Then
      Do J=1,I-1
      If(Occ(J).Lt.0.8.And.Occ(J).Gt.0.2) Then
C
      SIJ=Zero
      Do IG=1,NGrid
      SIJ=SIJ+Abs(OrbGrid(IG,I)*OrbGrid(IG,J))*WGrid(IG)
      EndDo 
C
      Write(6,'(/,X,"Abs Overlap S",2I4,2F8.5,E12.5)')
     $ I,J,Occ(I),Occ(J),SIJ
C
      EndIf
      EndDo
C
      EndIf
      EndDo
C
      Call LYP(RhoGrid,Sigma,Zk,NGrid)
c      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
c      Call GGA_SPIN(Zk1,URe,Occ,
c     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,rhoo,
C     >                   Sigma,sigmaco,sigmaoo,
C     >                   Zk1,vrhoc,vrhoo,
C     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
C
      IVer=1
      Call MCORRECTION(ELSM,Occ,TwoNO,URe,
     $                 OrbGrid,WGrid,
     $                 OrbXGrid,OrbYGrid,OrbZGrid,
     $                 NGrid,NBasis,NInte1,NInte2,IVer)
C
      ConCorr(1:200)=Zero
      A=0.2D0
      B=A-One
      C=2.6D0
      G=1.5D0
      D=(C-One)/(One-G)**2
C
      Write(6,'(/,1X,''Values of A, C, G parameters in CASPIDFT'',
     $ 2X,3F10.3)')A,C,G
C
      EDYN=Zero
      ELYP=Zero
      ESR=Zero
C
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C 
      XX=Two*OnTop(I)/RhoGrid(I)**2
C
c      If(Abs(RR(1,I)).eq.zero.And.Abs(RR(2,I)).eq.zero)
c     $ Write(*,'(5E15.6)')RR(1,I),RR(2,I),RR(3,I),RhoGrid(I),XX
c     $ Zk(I),WGrid(I)
c      If(Abs(RR(1,I)).eq.zero.And.Abs(RR(2,I)-0.8).Lt.0.02)
c     $ Write(*,'(5E15.6)')RR(1,I),RR(2,I),RR(3,I),RhoGrid(I),XX
C
      If(XX.Le.One) Then      
      PX=A*XX/(One+B*XX)
      Else
      PX=C*XX**0.25-D*(XX-G)**2
      EndIf
C
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
c      EDYN=EDYN+PX*Zk1(I)*WGrid(I)
      ELYP=ELYP+Zk(I)*WGrid(I)
c      ESR=ESR+Zk1(I)*WGrid(I)
C
      EndIf
C
      EndDo
C     TEST TRDMs
c      Call TEST_TRDMs(Occ,UNOAO,NInte1,NBasis)
C
      Write(6,'(/,1X,''LYP Correlation'',7X,F15.8)')ELYP
      Write(6,'(1X,''CASPIDFT Correlation'',2X,F15.8)')EDYN
C      Write(6,'(1X,''SR  Correlation'',7X,F15.8)')ESR
      Write(6,'(1X,''Total CAS+ LYP Energy'',X,F15.8)')ELYP+ETot+ENuc
      Write(6,'(1X,''Total CASPIDFT Energy'',X,F15.8)')EDYN+ETot+ENuc
c      Write(6,'(1X,''Total CASPI(M)DFT Energy'',X,F15.8)')
c     $ EDYN+ETot+ENuc+ELSM
C
      Stop 
C
      Return
      End

*Deck CASPIDFTOPT
      Subroutine CASPIDFTOPT(URe,UNOAO,Occ,NBasis)
C
C     CASPiDFT correlation energy is returned (and others, like LYP)
C     Two-electron integrals are not needed
C     If no parameters (order: A B G C D) are provided in the command line, the parameters used are:
C     ground state (state 1.1): 0.2 -0.8 1.5 2.6 6.4
C     excited state (state neq 1.1): 0.19185565 -0.81488859  1.50106868  2.56648139  6.31514676
C     see: Hapka et al., J. Phys. Chem. Lett. 2020, 11, 5883−5889
C
C     WARNING! The P(X) function is discontinuous at X=1. As a result numerical instability 
C     (related to the X>1.D0 condition) may occur. The final PiDFT energy value is reliable ONLY up to the 4th decimal place. 
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Dimension(:), Allocatable :: Zk,Zk1
      Real*8, Dimension(:), Allocatable :: rhoo,sigmaco,sigmaoo,vrhoc,
     $ vrhoo,vsigmacc,vsigmaco,vsigmaoo
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:)
      Real*8, Allocatable :: OnTop(:)
      Real*8, Dimension(:), Allocatable ::RhoAct,OnTopAct,SigmaAct,
     $ ZkAct,ZkPBE
      Dimension OccAct(NBasis)
C
      CHARACTER(100) :: num1char,num2char,num3char,num4char,num5char
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ ConCorr(200)
     $ ,EpsC(NBasis,NBasis),EpsCI(NBasis),EpsDiag(NBasis)
      logical fderiv,open
      integer igrad
      character*(30) name

c      Call DELCORREL(URe,UNOAO,Occ,NBasis)
c      stop 
C
C     IFlagRead = 0 - generate grid data, do not write to a file
C     if parameteres in caspidft are optimized it is useful to use flags different from 0 to avoid generating data on grid more than once
C     IFlagRead = 1 - read grid data (density, ontop, and lyp energy density) from a file
C     IFlagRead = 2 - write grid data to a file
C
      IFlagRead=0
C
      Write(6,'(/,X,"***************** CASPIDFT ******************* ")')
      Write(6,'(X,"STATE : ",I1,".", I1,/)')inst(1,1),inst(2,1)
C
      IF(COMMAND_ARGUMENT_COUNT().Eq.0) THEN
C
      If(inst(1,1).Eq.1.And.inst(2,1).Eq.1) Then
      Write(6,'(/,X,"Ground-state optimal parameteres 
     $ will be used in PiDFT")')
      A1= 0.2D0
      B1=-0.8D0
      G1= 1.5D0
      C1= 2.6D0 
      D1= 6.4D0
      Else
      Write(6,'(/,X,"Excited-state optimal parameteres
     $ will be used in PiDFT")')
      A1= 0.19185565
      B1=-0.81488859  
      G1= 1.50106868  
      C1= 2.56648139  
      D1= 6.31514676
      EndIf
C
      ELSE
      Write(6,'(/,X,"Parameters read from a command line 
     $ will be used in PiDFT")')
C
      CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
      CALL GET_COMMAND_ARGUMENT(2,num2char)
      CALL GET_COMMAND_ARGUMENT(3,num3char)
      CALL GET_COMMAND_ARGUMENT(4,num4char)
      CALL GET_COMMAND_ARGUMENT(5,num5char)

      READ(num1char,*)xnum1                    !then, convert them to REALs
      READ(num2char,*)xnum2
      READ(num3char,*)xnum3
      READ(num4char,*)xnum4
      READ(num5char,*)xnum5
C
      A1=XNum1
      B1=XNum2
      G1=XNum3
      C1=XNum4 
      D1=XNum5
C
      ENDIF 
C
      Call molprogrid0(NGrid,NBasis)
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
      Allocate  (RhoAct(NGrid))
      Allocate  (OnTopAct(NGrid))
      Allocate  (SigmaAct(NGrid)) 
      Allocate  (ZkAct(NGrid))
      Allocate  (ZkPBE(NGrid))
C
c      Call molprogrid1(RR,NGrid)
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
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
      If(IFlagRead.Eq.1) Then
      Open(10,File="data_on_grid.dat",form='unformatted')
      Read(10) (OnTop(I),I=1,NGrid)
      Read(10) (RhoGrid(I),I=1,NGrid)
      Read(10) (Zk(I),I=1,NGrid)
      Close(10)
      GoTo 300
      EndIf
C
C     APPROXIMATE PROJECTION ON CAS(4,4) [used only in definition of ionicity index]
C
      Write(6,'(/,X,
     $ "Occupancies projected on CAS(4,4) for ionicity index")')
      Do I=1,NBasis
      OccAct(I)=Occ(I)
c      If(I.Le.NInAcCAS) OccAct(I)=Zero
C ONLY HOMO-LUMO
c      If(I.Gt.NELE+1) OccAct(I)=Zero
c      If(I.Gt.NInAcCAS.And.I.Lt.NELE) OccAct(I)=One
C HOMO-1 TO LUMO+1
      If(I.Gt.NELE+2) OccAct(I)=Zero
      If(I.Gt.NInAcCAS.And.I.Lt.NELE-1) OccAct(I)=One 
      If(OccAct(I).Ne.Zero) 
     $ Write(6,'(X,I4,3X,E12.5)')I,OccAct(I) 
      EndDo 
C
      Do I=1,NGrid
      Call DenGrid(I,RhoAct(I),OccAct,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,OccAct,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,OccAct,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,OccAct,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      SigmaAct(I)=RhoX**2+RhoY**2+RhoZ**2
      OnTopAct(I)=Zero
      EndDo
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
c      Gam2=FRDM2(IP,IQ,IR,IS,RDM2Act,OccAct,Ind2,NAct,NBasis)
      Gam2=FRDM2R(IP,IQ,IR,IS,RDM2Act,OccAct,Ind2,NAct,NBasis)
      If(Abs(Gam2).Gt.1.D-8) Then
      Do I=1,NGrid
      OnTopAct(I)=OnTopAct(I)
     $ +Two*Gam2*OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
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
C
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
C
      Gam2=FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      If(Abs(Gam2).Gt.1.D-8) Then
C
      Do I=1,NGrid
      OnTop(I)=OnTop(I)
     $ +Two*Gam2*OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     OVERLAP OF ORBITALS
C
c      Do I=1,NOccup
c      If(Occ(I).Lt.0.8.And.Occ(I).Gt.0.2) Then
c      Do J=1,I-1
c      If(Occ(J).Lt.0.8.And.Occ(J).Gt.0.2) Then
cC
c      SIJ=Zero
c      Do IG=1,NGrid
c      SIJ=SIJ+Abs(OrbGrid(IG,I)*OrbGrid(IG,J))*WGrid(IG)
c      EndDo
cC
c      Write(6,'(/,X,"Abs Overlap S",2I4,2F8.5,E12.5)')
c     $ I,J,Occ(I),Occ(J),SIJ
cC
c      EndIf
c      EndDo
cC
c      EndIf
c      EndDo
C
      Call LYP(RhoGrid,Sigma,Zk,NGrid)
C 
      Call LYP(RhoAct,SigmaAct,ZkAct,NGrid)
      Call PBECor(RhoAct,SigmaAct,ZkPBE,NGrid)
C
      If(IFlagRead.Eq.2) Then
      Open(10,File="data_on_grid.dat",form='unformatted')
      Write(10) (OnTop(I),I=1,NGrid)
      Write(10) (RhoGrid(I),I=1,NGrid)
      Write(10) (Zk(I),I=1,NGrid)
      Close(10)
      EndIf
  300 Continue
C
      A=A1
      B=B1
      C=C1
      G=G1
      D=D1
C
      Write(6,'(/,1X,''Values of A, B, C, G, D parameters in PiDFT:'',
     $ 5F12.6)')A,B,C,G,D
C
      EDYN=Zero
      ELYP=Zero
      EIonLYP=Zero
      EIon=Zero
C     
      EDYNAct=Zero       
      EIonAct=Zero
      ELYPAct=Zero
      EIonLYPAct=Zero
      EPBEAct=Zero
      EIonPBEAct=Zero
C
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C
      XX=Two*OnTop(I)/RhoGrid(I)**2
c      If(Abs(RR(1,I)).Lt.1.D-10.And.Abs(RR(2,I)).Lt.1.D-10)
c     $ Write(*,'(F10.4,E13.4,F10.4)') RR(3,I),RhoGrid(I),XX
C
      If(XX.Le.One) Then
C
      PX=A*XX/(One+B*XX)
      If(RhoGrid(I).Lt.1.d-10)PX=Zero
C
      Else
C
      PX=C*XX**0.25-D*(XX-G)**2
c herer!!! condition P(2)=1.48 is the maximum allowed value
c      If(XX.Gt.2.0.Or.PX.Gt.1.48) PX=1.48
c      If(Abs(Zk(I)).Lt.1.D-7.And.XX.Gt.1.5) PX=One
C     
      EndIf
C
      PXA=Zero
      If(RhoAct(I).Gt.1.d-10) Then
      XXA=Two*OnTopAct(I)/RhoAct(I)**2
      If(XXA.Le.One) Then
      PXA=A*XXA/(One+B*XXA)
      Else
      PXA=C*XXA**0.25-D*(XXA-G)**2
      EndIf 
      EndIf
C
      If(XX.Gt.One) EIonLYP=EIonLYP+Zk(I)*WGrid(I)
      If(XX.Gt.One) EIon=EIon+PX*Zk(I)*WGrid(I)
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
      ELYP=ELYP+Zk(I)*WGrid(I)
C 
      EDYNAct=EDYNAct+PXA*ZkAct(I)*WGrid(I)
      If(XXA.Gt.One) EIonAct=EIonAct+PXA*ZkAct(I)*WGrid(I)
      ELYPAct=ELYPAct+ZkAct(I)*WGrid(I)
      If(XXA.Gt.One) EIonLYPAct=EIonLYPAct+ZkAct(I)*WGrid(I)
      EPBEAct=EPBEAct+ZkPBE(I)*WGrid(I)
      If(XXA.Gt.One) EIonPBEAct=EIonPBEAct+ZkPBE(I)*WGrid(I)
C
      EndIf
C
      EndDo
      Write(6,'(/,1X,''LYP Correlation'',12X,F15.8)')ELYP
      Write(6,'(1X,''PiDFT Correlation'',10X,F15.8)')EDYN
C    
      Write(6,'(/,36X,"PiDFT           LYP")') 
      Write(6,'(1X,''Corr Energy, X>1'',11X,2F15.8)') 
     $ EIon,EIonLYP
      Write(6,'(1X,''Tot Ionicity Index'', 5X,F12.1,"%",2X,F12.1,"%")')
     $ EIon/EDYN*100.,EIonLYP/ELYP*100
C
      Write(6,'(/,1X,''Corr Energy CAS(4,4)_P'',5X,2F15.8)')EDYNAct,
     $ ELYPAct
      Write(6,'(1X,''Corr Energy CAS(4,4)_P, X>1'',2F15.8)') EIonAct,
     $ EIonLYPAct
      Write(6,
     $ '(1X,''Ionicity Index CAS(4,4)_P'', F10.1,"%",4X,F10.1,"%")') 
     $ EIonAct/EDYNAct*100.,EIonLYPAct/ELYPAct*100
C
C PBE PRINT
C
c      Write(6,'(/,1X,''PBE Corr CAS(4,4)_P'',8X,F15.8)')EPBEAct
c      Write(6,'(1X,''PBE Corr CAS(4,4)_P, X>1'',3X,F15.8)')EIonPBEAct
c      Write(6,
c     $ '(1X,''Ionicity Index PBE CAS(4,4)_P'',F11.1,"%")')
c     $ EIonPBEAct/EPBEAct*100
C
C
      Return
      End

*Deck GGA_SPIN
      Subroutine GGA_SPIN(Zk,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS CORRELATION ENERGY DENSITY USING SPIN FUNCTIONAL
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
      IOnTop=1
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
      If(IOnTop.Eq.1) Then
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
C     if ontop.eq.1 
      Else
C full spin polarization
      RhoA(I)=RhoGrid
      RhoB(I)=Zero   
      RhoXa=RhoX
      RhoXb=Zero
      RhoYa=RhoY
      RhoYb=Zero
      RhoZa=RhoZ
      RhoZb=Zero

      Fac=1.999
      RhoA(I)=RhoGrid/2.d0*Fac
      RhoB(I)=RhoGrid/2.D0*(Two-Fac)
      RhoXa=RhoX/2.d0*Fac
      RhoXb=RhoX/2.d0*(Two-Fac)
      RhoYa=RhoY/2.d0*Fac
      RhoYb=RhoY/2.d0*(Two-Fac)
      RhoZa=RhoZ/2.d0*Fac
      RhoZb=RhoZ/2.d0*(Two-Fac)

      SigmaAA(I)=RhoXa*RhoXa+RhoYa*RhoYa+RhoZa*RhoZa
      SigmaAB(I)=RhoXa*RhoXb+RhoYa*RhoYb+RhoZa*RhoZb
      SigmaBB(I)=RhoXb*RhoXb+RhoYb*RhoYb+RhoZb*RhoZb
C
      EndIf
C
      Else
C
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
      Call LYP_SPIN(RhoA,RhoB,SigmaAA,SigmaAB,SigmaBB,Zk,NGrid)
C
      Return
      End

*Deck MCORRECTION
      Subroutine MCORRECTION(ELSM,Occ,TwoNO,URe,OrbGrid,WGrid,
     $                       OrbXGrid,OrbYGrid,OrbZGrid,
     $                       NGrid,NBasis,NInte1,NInte2,IVer)
C
C     "MAP" CAS 2-RDM ON GVB 2-RDM BY COUPLING ORBITALS INTO GEMINALS
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),TwoNO(NInte2),URe(NBasis,NBasis)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $          OrbXGrid(NGrid,NBasis),
     $          OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),RDM2(NBasis**2*(NBasis**2+1)/2),Ind1(NBasis)
      Double Precision, Allocatable :: OccGH(:),OrbGem(:,:),
     $                                 Zk(:),Sigma(:),
     $                                 RhoGridG(:),RhoGridH(:)
C
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      RDM2(1:NRDM2)=Zero
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct 
C
C     scaling of 2nd variant
      ACoef = 2.5d0
C
      Do I=1,NAct
      Ind1(I)=INActive+I
      EndDo
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
      RDM2(NAddrRDM(J,L,I,K,NBasis))=Half*X
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      Hlp=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  Hlp=Hlp+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ Hlp=Hlp-Occ(IP)*Occ(IQ)
      If(Hlp.Ne.Zero) RDM2(IAdd)=Hlp
      EndDo
      EndDo
      EndDo
      EndDo
C
      If(NELE-INActive.Ne.NAct-(NELE-INActive)) Then
      Write(6,*) 'Fatal Error: mapping of CAS
     $ on GVB only defined for CAS(m,m)'
c      Stop
      EndIf
C
      NGemSave=NGem
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
      NGem=INActive+NAct/2
      Do I=1,NBasis
      If(Occ(I).Ne.Zero) IGem(I)=0
      If(Occ(I).Eq.Zero) IGem(I)=NGem+1
      EndDo
      Do I=1,INActive
      IGem(I)=I
      EndDo
C
C     COUPLE ORBITALS
C
      Write(6,'(/,X,"Mappinng of CAS(n,n) 2-RDM on GVB-like 2-RDM")')
C
      Allocate(OrbGem(NELE-INActive,2))
C
      II=1
c herer!!!
c      INActive=4 

      Do IP=INActive+1,NELE
      IGem(IP)=IP
C
      XDevMax=Zero
c herer!!!
      Do IQ=NELE+1,NOccup
c herer!!! do not go through all weakly occupied orbitals but only through the strongly-correlated ones
c      Do IQ=NELE+1,NELE+NELE-INActive
C
      XPQPQ=Abs(RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis))-
     $ 2.0D0*Occ(IP)*Occ(IQ))
      XPQQP=Abs(RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis))-
     $ (-Occ(IP)*Occ(IQ)))
      XDev=Half*(XPQPQ+XPQQP)/Occ(IQ)
      If(XDev.Gt.XDevMax) Then
      XDevMax=XDev
      IQMax=IQ
      EndIf
C
      CP=SQRT(Occ(IP))
      If(Occ(IP).Lt.Half) CP=-CP
      CQ=SQRT(Occ(IQ))
      If(Occ(IQ).Lt.Half) CQ=-CQ
      Write(*,*)IP,IQ
      Write(*,*)'PPQQ',RDM2(NAddrRDM(IP,IP,IQ,IQ,NBasis)),CP*CQ
      Write(*,*)'PQPQ',RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis)),
     $ 2.0D0*Occ(IP)*Occ(IQ)
      Write(*,*)'PQQP',RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis)),
     $ -Occ(IP)*Occ(IQ)
C
      EndDo
C
      If(IGem(IQMax).Eq.0) Then
      IGem(IQMax)=IP
      Write(6,'(X,"**** Orbital",I2," coupled with ",I2 )')IP,IQMax
      OrbGem(II,1)=IP
      OrbGem(II,2)=IQMax
      II = II + 1
      Else
      Write(6,*)
     $ "Warning: more than 2 orbitals assigned to a geminal no",IP
      EndIf
C
      EndDo
C
C
      Do I=1,NGem
C
      Write(6,'(/,X,"Geminal no",I2," includes")')I
      Sum=Zero
C
      II=0
      Do J=1,NBasis
      If(IGem(J).Eq.I) Then
      Sum=Sum+Occ(J)
      Write(6,'(X,"Orbital No: ",I4)')J
      EndIf
      EndDo
      Write(6,'(X,"Norm: ",F12.6)')Sum
C
      EndDo
C
C      If(IVer.Eq.0) Then
C
C     COMPUTE THE M CORRECTION AS GIVEN IN Eq.(26), van Meer et al. JCP 148, 104102 (2018)
C
      ELSM=Zero
      Do I=1,NOccup
      Do J=1,NOccup
      If(IGem(I).Ne.IGem(J)) Then
      ELSM=ELSM+FM(Occ(I),Occ(J))*TwoNO(NAddr3(I,J,I,J))
      EndIf
      EndDo
      EndDo
C
      Write(6,'(/,1X,''M Correlation Correction '',5X,F15.8)')ELSM
C
C      ElseIf(IVer.Eq.1) Then
C
C     MODIFIED M CORRECTION
c
      Allocate(Zk(NGrid),RhoGridH(NGrid),RhoGridG(NGrid))
      Allocate(OccGH(NBasis),Sigma(NGrid))
C
      ELSM=Zero
      Do IG=1,NELE-INActive
      Do IH=1,IG-1
C
      IQ=OrbGem(IG,1)
      IP=OrbGem(IG,2)
      IS=OrbGem(IH,1)
      IR=OrbGem(IH,2)
C 
      XX=max(Occ(IQ),Occ(IP))*(One-max(Occ(IQ),Occ(IP)))
      YY=max(Occ(IS),Occ(IR))*(One-max(Occ(IS),Occ(IR)))
      Coef=SQRT(XX*YY)
C
      OccGH=0
      OccGH(IP)=Occ(IP)
      OccGH(IQ)=Occ(IQ)
      OccGH(IR)=Occ(IR)
      OccGH(IS)=Occ(IS)
C 24.04.2020 
C MODIFICATION: IF WEAKLY OCC ORBITALS ARE DEGENERATE INCLUDE THEM ALL IN THE DENSITY!
      Do I=NELE+1,NBasis
C
      If(Occ(I).Gt.1.D-8.And.Abs(Occ(I)-Occ(IP))/Occ(IP).Le.1.D-2.
     $ And.I.Ne.IP)Then
C
C MAKE SURE THAT I IS NOT ASSIGNED TO A GEMINAL   
C
      IOK=1
      Do IGG=1,NELE-INActive
      If(OrbGem(IGG,1).Eq.I.Or.OrbGem(IGG,2).Eq.I) IOK=0
      EndDo
C
      If(IOK.Eq.1) Then
      OccGH(I)=Occ(I)
      Write(*,*)'Include additional orbital',I,'in density for Mmdc'
      EndIf
C
      EndIf
C
      If(Occ(I).Gt.1.D-8.And.Abs(Occ(I)-Occ(IR))/Occ(IR).Le.1.D-2.
     $ And.I.Ne.IR)Then
C
      IOK=1
      Do IGG=1,NELE-INActive
      If(OrbGem(IGG,1).Eq.I.Or.OrbGem(IGG,2).Eq.I) IOK=0
      EndDo
C
      If(IOK.Eq.1) Then
      OccGH(I)=Occ(I)
      Write(*,*)'Include additional orbital',I,'in density for Mmdc'
      EndIf
C
      EndIf 

      EndDo
C
      Do I=1,NGrid
      Call DenGrid(I,RhoGridH(I),OccGH,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,OccGH,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,OccGH,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,OccGH,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
      EndDo
C
      Call LYP(RhoGridH,Sigma,Zk,NGrid)
C
      VAL=0
      Do I=1,NGrid
      VAL = VAL + Zk(I)*WGrid(I)
      EndDo
C   
      ELSM = ELSM + Coef*VAL
C
      EndDo
      EndDo
C
      ELSM=ACoef*ELSM
C
      Deallocate(OrbGem,OccGH)
      Deallocate(Sigma)
      Deallocate(Zk,RhoGridH,RhoGridG)
C
C      EndIf
C
      Write(6,'(/,1X,''MDC Correlation Correction '',3X,F15.8)')ELSM
C
      NGem=NGemSave
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
      Return
      End

*Deck FM
      Real*8 Function FM(X,Y)
C
C     RETURNS A VALUE OF THE FM FUNCTION DEFINED IN Eqs.(27),(28) van Meer et al. JCP 148, 104102 (2018)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Gamma=1500.D0
C
      XX=X*(One-X)
      YY=Y*(One-Y)
C
      PX=(One+16.D0/Gamma)*Gamma*XX**2/(One+Gamma*XX**2)
      PY=(One+16.D0/Gamma)*Gamma*YY**2/(One+Gamma*YY**2)
C
      FM=-PX*PY*SQRT(XX*YY)
C
      Return
      End

*Deck TEST_TRDMS
      Subroutine TEST_TRDMS(Occ,UNOAO,NInte1,NBasis) 
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
     $ UNOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Double precision,dimension(NBasis,NBasis) :: dipx,dipy,dipz
      Double precision :: GammaAB(NInte1),PC(NBasis)
      Double precision :: trdm(NBasis,NBasis),AUXM(NBasis,NBasis),
     $                    WorkV(NBasis)
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
      Call read_1trdm_molpro(AUXM,InSt(1,1),InTrSt(1,1),
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
C     READ DIPOLE MOMENT
      Call read_dip_molpro(dipx,dipy,dipz,'DIP',NBasis)
C
      Call dgemm('N','N',NBasis,NBasis,NBasis,1d0,UNOAO,NBasis,
     $           dipz,NBasis,0d0,AUXM,NBasis)
      Call dgemm('N','T',NBasis,NBasis,NBasis,1d0,AUXM,NBasis,
     $           UNOAO,NBasis,0d0,dipz,NBasis)
C
      TSDipZ=0
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      TSDipZ = TSDipZ + Two*trdm(J,I)*dipz(I,J)
      Enddo
      Enddo
      Write(6,
     $ '(X,"Transition State DMZ: ",F15.8)')TSDipZ
C 
C     REFERENCE STATE
C
      GSDipZ=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,NBasis
      IJ=IJ+1
      If(I.Eq.J) GSDipZ=GSDipZ+Two*Occ(I)*dipz(I,J)
      Enddo
      Enddo
C
      Write(6,
     $ '(X,"Reference  State DMZ: ",F15.8)')GSDipZ 
C
      End

*Deck FRDM2R
      Real*8 Function FRDM2R(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C
C     FOR A GIVEN SET OF INDICES AND THE KNOWN PART OF ACTIVE RDM2
C     RETURNS THE ELEMENT OF RDM2_PQRS FOR CAS
C
C     It does the same as FRDM2 but here it is used to test different options
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension RDM2Act(NAct**2*(NAct**2+1)/2),Occ(NBasis),Ind2(NBasis)
C
      RDM2=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  RDM2=RDM2+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      If(Occ(IP).Eq.One.Or.Occ(IQ).Eq.One.Or.Occ(IR).Eq.One.
     $ Or.Occ(IS).Eq.One) Then
      FRDM2R=RDM2
      Return
      EndIf
C
C     ACTIVE PART
C
      If(Ind2(IP)*Ind2(IQ)*Ind2(IR)*Ind2(IS).Ne.Zero) Then
      RDM2=RDM2Act(NAddrRDM(Ind2(IP),Ind2(IQ),Ind2(IR),Ind2(IS),NAct))
      EndIf
C
      FRDM2R=RDM2
C
      Return
      End

*Deck CORRELON
      Subroutine CORRELON(URe,UNOAO,Occ,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Dimension(:), Allocatable :: Zk,Zk1
      Real*8, Dimension(:), Allocatable :: rhoo,sigmaco,sigmaoo,vrhoc,
     $ vrhoo,vsigmacc,vsigmaco,vsigmaoo
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:),RX(:),RY(:),RZ(:)
      Real*8, Allocatable :: OnTop(:), CCov(:),CIon(:), XX(:)
C new
      Real*8, Dimension(:), Allocatable ::RhoAct,OnTopAct,SigmaAct,
     $ ZkAct,ZkPBE
      Dimension OccAct(NBasis)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ ConCorr(200)
      integer igrad
      character*(30) name
C
      Call molprogrid0(NGrid,NBasis)
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
      Allocate  (RX(NGrid))
      Allocate  (RY(NGrid))
      Allocate  (RZ(NGrid))
      Allocate  (OnTop(NGrid))
      Allocate  (CCov(NGrid))
      Allocate  (CIon(NGrid))
      Allocate  (XX(NGrid))
C new
      Allocate  (RhoAct(NGrid))
      Allocate  (OnTopAct(NGrid))
      Allocate  (SigmaAct(NGrid))
c      Allocate  (ZkAct(NGrid))
c      Allocate  (ZkPBE(NGrid))
C
      Call molprogrid1(RR,NGrid)
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
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
C new
      Do I=1,NBasis
      OccAct(I)=Occ(I)
      If(I.Gt.NELE+2) OccAct(I)=Zero
      If(I.Gt.NInAcCAS.And.I.Lt.NELE-1) OccAct(I)=One
c      If(occact(i).ne.zero) write(*,*)i,occact(i)
      EndDo
C
      Do I=1,NGrid
      Call DenGrid(I,RhoAct(I),OccAct,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,OccAct,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,OccAct,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,OccAct,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      SigmaAct(I)=RhoX**2+RhoY**2+RhoZ**2
      OnTopAct(I)=Zero
      EndDo
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      Gam2=FRDM2R(IP,IQ,IR,IS,RDM2Act,OccAct,Ind2,NAct,NBasis)
      If(Abs(Gam2).Gt.1.D-8) Then
      Do I=1,NGrid
      OnTopAct(I)=OnTopAct(I)
     $ +Two*Gam2*OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
c end of new
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
C
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
C
      Gam2=FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      If(Abs(Gam2).Gt.1.D-8) Then
C
      Do I=1,NGrid
      OnTop(I)=OnTop(I)
     $ +Two*Gam2*OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do I=1,NGrid
C
c      AP=0.01
      AP=0.01
      XX(I)=Zero
      CCov(I)=Zero
      CIon(I)=Zero
      XPade=RhoGrid(I)/(AP+RhoGrid(I))
      If(RhoGrid(I).Gt.1.D-7) Then
C
      XX(I)=Two*OnTop(I)/RhoGrid(I)**2
c      write(*,'(7F8.4)')
c     $ rr(1,i),rr(2,i),rr(3,i),rhogrid(i),xpade,xx,xx*xpade
C
      If(XX(I).Le.One) Then
C
      CIon(I)=Zero
      CCov(I)=SQRT(Half*(One-XX(I))*XPade)
C
      Else
C
      CCov(I)=Zero
      CIon(I)=SQRT(Half*(XX(I)-One)*XPade)
C
      EndIf
      EndIf
C
      EndDo
C
      XNormC=Zero
C
      Do I=1,NGrid
C
      RX(I)=RR(1,I)
      RY(I)=RR(2,I)
      RZ(I)=RR(3,I)
C
      XNormC=XNormC+(CCov(I)**2+CIon(I)**2)
C
      EndDo
      XNormC=SQRT(XNormC)
C
      Call SortVecs(RZ,RX,RY,RhoGrid,XX,CCov,CIon,NGrid)
C
      Write(*,'("   z(bohr)      Rho",10X,"X       cs_cov   cs_ionic")')
      Do I=1,NGrid
      CCov(I)=CCov(I)/XNormC
      CIon(I)=CIon(I)/XNormC
      If(Abs(RX(I)).Lt.1.D-10.And.Abs(RY(I)).Lt.1.D-10)
     $ Write(*,'(F10.4,E13.4,3F10.4)')
     $ RZ(I),RhoGrid(I),XX(I),CCov(I),CIon(I) 
      EndDo
C
      Return
      End 

*Deck DELCORREL
      Subroutine DELCORREL(URe,UNOAO,Occ,NBasis)
C
C     DELTA-CORRELONS
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Allocatable :: RDM2Act(:),OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:),RX(:),RY(:),RZ(:)
      Real*8, Allocatable :: OnTop(:), CRe(:),CIm(:), DELXX(:)
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis)
C
      Call molprogrid0(NGrid,NBasis)
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NGrid,NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (RR(3,NGrid))
      Allocate  (RX(NGrid))
      Allocate  (RY(NGrid))
      Allocate  (RZ(NGrid))
      Allocate  (OnTop(NGrid))
      Allocate  (CRe(NGrid))
      Allocate  (CIm(NGrid))
      Allocate  (DELXX(NGrid))
C
      Call molprogrid1(RR,NGrid) 
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
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
      Do I=1,NGrid
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      OnTop(I)=Zero
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
C
      Gam2=FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      If(Abs(Gam2).Gt.1.D-8) Then
C
      Do I=1,NGrid
      OnTop(I)=OnTop(I)
     $ +Two*Gam2*OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,"STATE : ",I1,".", I1)')inst(1,1),inst(2,1)
C
      AP=0.01
C
      If(inst(1,1).Eq.1.And.inst(2,1).Eq.1) Then
C
C     if ground state: compute X(r)*Pade and save it
C
      Write(6,'(/,X,"**** COMPUTE GROUND STATE X(r) **** ")')
C
      Open(10,File="x_gs.dat")
C
      Do I=1,NGrid
C 
      XX=Zero
      XPade=RhoGrid(I)/(AP+RhoGrid(I))
      If(RhoGrid(I).Gt.1.D-7) XX=Two*OnTop(I)/RhoGrid(I)**2*XPade
      Write(10,*)I,XX,XPade
C
      EndDo
C
      Close(10)
C
C     If(inst(1,1).Eq.1.And.inst(2,1).Eq.1) Then
      Else
C
C     if excited state: compute delta-correlon wavefunction
C
      Write(6,'(/,X,"**** COMPUTE X_exct-X_gs **** ")')
      Open(10,File="x_gs.dat",Status='Old')
C       
      Do I=1,NGrid
C
      Read(10,*)II,XXGS,XPade
C
      XXES=Zero
c      XPade=RhoGrid(I)/(AP+RhoGrid(I))
      If(RhoGrid(I).Gt.1.D-7) XXES=Two*OnTop(I)/RhoGrid(I)**2*XPade
      DELXX(I)=XXES-XXGS
C
      EndDo
C
      Close(10)
C
      XNorm=Zero
      Do I=1,NGrid
      XNorm=XNorm+WGrid(I)*Abs(DELXX(I))
      EndDo
C
      Do I=1,NGrid
C
      RX(I)=RR(1,I)
      RY(I)=RR(2,I)
      RZ(I)=RR(3,I)
C 
      CRe(I)=Zero
      CIm(I)=Zero
C
      If(DELXX(I).Ge.Zero) Then
      CRe(I)=Sqrt(DELXX(I)/XNorm)
      Else
      CIm(I)=Sqrt(-DELXX(I)/XNorm)
      EndIf
C
      EndDo 
C
      Call SortVecs(RZ,RX,RY,RhoGrid,DELXX,CRe,CIm,NGrid)
C
      Write(*,'("   z(bohr)     Rho",10X,"DELX     Psi_Re   Psi_Im")')
      Do I=1,NGrid
c      if(RhoGrid(I).gt.1.d-2.and. CRe(I)+CIm(I).gt.0.05) 
c     $ Write(*,'(F10.4,E13.4,3F10.4)')
c     $ RZ(I),RX(I),RY(I),CRe(I),CIm(I)

      If(Abs(RX(I)).Lt.1.D-10.And.Abs(RY(I)).Lt.1.D-10)
c      If(Abs(RX(I)).Lt.0.5.And.Abs(RX(I)).Gt.0.10)
c.And.
c     $ Abs(RY(I)).Lt.1.D-10)
     $ Write(*,'(F10.4,E13.4,3F10.4)')
     $ RZ(I),RhoGrid(I),DELXX(I),CRe(I),CIm(I)

      EndDo
C
C     If(inst(1,1).Eq.1.And.inst(2,1).Eq.1) Then
      EndIf
C
      Return
      End

*Deck SortVecs
      Subroutine SortVecs(XX,YY1,YY2,YY3,YY4,YY5,YY6,N)
C
C     SORT XX IN AN ASCENDING ORDER 
C     AND CHANGE THE ORDER OF YY's 
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension XX(N),YY1(N),YY2(N),YY3(N),YY4(N),YY5(N),YY6(N)
C
C     LOCAL ARRAYS
C
      Dimension Save1(N),Save2(N),Save3(N),Save4(N),Save5(N),Save6(N),
     $ Ind(N)
C
      Do I=1,N
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,N
C
      EMin=XX(IStart)
      IndMin=IStart
C
      Do J=IStart,N
      If(XX(J).Lt.EMin) Then
      EMin=XX(J)
      IndMin=J
      EndIf
      EndDo
C
      Hlp=XX(IStart)
      IndHlp=Ind(IStart)

      XX(IStart)=XX(IndMin)
      Ind(IStart)=Ind(IndMin)

      XX(IndMin)=Hlp
      Ind(IndMin)=IndHlp
C
      IStart=IStart+1
C
      EndDo
C
C     YY
C
      Do I=1,N
      Save1(I)=YY1(I)
      Save2(I)=YY2(I)
      Save3(I)=YY3(I)
      Save4(I)=YY4(I)
      Save5(I)=YY5(I)
      Save6(I)=YY6(I)
      EndDo
      Do I=1,N
      YY1(I)=Save1(Ind(I))
      YY2(I)=Save2(Ind(I))
      YY3(I)=Save3(Ind(I))
      YY4(I)=Save4(Ind(I))
      YY5(I)=Save5(Ind(I))
      YY6(I)=Save6(Ind(I))
      EndDo
C
      Return
      End
