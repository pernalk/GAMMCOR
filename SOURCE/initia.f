*Deck ReadDAL
      Subroutine ReadDAL(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem)
C
C     READ HAO, 2-EL INTEGRALS IN NO, C_COEFFICIENTS, IGEM FROM A DALTON_GENERATED FILE
C     READ UMOAO FROM DALTON.MOPUN
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 XKin(NInte1),XNuc(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),UMOAO(NBasis,NBasis)
C
      Character*60 Line
      Character*30 Line1
C
      Include 'commons.inc'
C
C     ICASSCF=1 - read occ from occupations.dat
      ICASSCF=1
C
C     IDMRG - integrals and RDM's read from external dmrg files
      IDMRG=0
      If(IDMRG.Eq.1.And.ICASSCF.Eq.0) Stop 'Set ICASSCF TO 1'
C
      If(IDMRG.Eq.1) Then
      Call ReadDMRG(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem)
      Return
      EndIf
C
C     READ UMOAO FROM DALTON.MOPUN
C
      Open(10,File='DALTON.MOPUN',Form='Formatted',Status='Old')
      Read(10,'(A60)') Line
      IEND=0
      Do J=1,NBasis
      IST=IEND+1
      IEND=IEND+NBasis
      Read(10,'(4F18.14)') (UMOAO(J,I),I=1,NBasis)
      EndDo
      Close(10)
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=0.D0
      If(I.Eq.J) URe(I,J)=1.0D0
      EndDo
      EndDo
C
C     READ 1-EL NO HAMILTONIAN FROM AOONEINT
C
      Call read1el(XKin,UMOAO,NBasis,NInte1)
C
C     GET 2-EL NO INTEGRALS AND CICoef
C 
      Call read2el(TwoEl,UMOAO,NBasis,NInte2)
C
      If(ICASSCF.Eq.0) Then
C
      CICoef(1:NBasis)=0.D0
      Open(10,File='coeff.dat',Form='Formatted',Status='Old')
      Read(10,*) NActive
      INActive=NELE-NActive
      Do I=1,INActive
      CICoef(I)=1.D0
      IGem(I)=I
      EndDo
      Read(10,*) (CICoef(I+INActive),I=1,2*NActive)
      Do I=INActive+1,NELE
      IGem(I)=I
      IGem(NELE+I-INActive)=I
      EndDo
      NGem=NELE+1
      Do I=1,NBasis
      If(CICoef(I).Eq.0.D0) IGem(I)=NGem
      Occ(I)=CICoef(I)**2
      EndDo
      Close(10)
C
      ElseIf(ICASSCF.Eq.1) Then
C
      Occ(1:NBasis)=0.D0
C
      Open(10,File='occupations.dat',Form='Formatted',Status='Old') 
C
      Read(10,*) NInAc,NAc
      NInAc=NInAc/2
      Read(10,*) (Occ(I),I=1,NInAc+NAc)
      Sum=0.D0
      Do I=1,NInAc+NAc
      Occ(I)=Occ(I)/2.D0
      Sum=Sum+Occ(I)
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
C      Write(6,'(2X,"No of CAS inactive and active orbitals: ",4X,4I)')
C     $ NInAcCAS,NAcCAS 
      Write(6,'(2x,a,4x,2i3)')
     $ "No of CAS inactive and active orbitals",NInAcCAS,NAcCAS
      Write(6,'(2X,"CASSCF",3X,"Occupancy",4X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,I6)') I,Occ(I),IGem(I)
      EndDo
      Write(6,'(2X,"Sum of Occupancies: ",E16.6)') Sum
      If(Abs(Sum-NELE).Gt.1.D-8) 
     $ Stop "Fatal Error: Occupancies do not sum up to NELE"
C
      Close(10)
C
c     If(ICASSCF.Eq.0)
      EndIf
C
      Open(10,File='enuc.dat',Form='Formatted',Status='Old')
      Read(10,'(A31,F20.12)')Line,ENuc
      Write(6,'(/,"  Nuclear repulsion:",F20.12)')ENuc
      Close(10)
C
      Return
      End

*Deck ReadDMRG
      Subroutine ReadDMRG(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem)
C
C     READ Integrals and 1-RDM - needed for AC-DMRG CALCULATION
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.D0,One=1.D0,Two=2.D0)
C
      Real*8 XKin(NInte1),XNuc(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),UMOAO(NBasis*NBasis)
C
      Character*60 FName,Aux1
C
      Include 'commons.inc'
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2(:),RDMAB2(:)
      Dimension Gamma(NInte1),Work(NBasis),PC(NBasis),
     $ AUXM(NBasis,NBasis)
c      Real*8, Allocatable :: TwoAux(:)
C
      UMOAO(1:NBasis*NBasis)=Zero
      URe(1:NBasis,1:NBasis)=Zero
      Occ(1:NBasis)=Zero
      PC(1:NBasis)=Zero
      Gamma(1:NInte1)=Zero
C
C     READ IN 1-RDM AND DIAGONALIZE IT
C     
      Open(10,File='rdmdump.dat')
   20 Read(10,*,End=30) X,I1,I2,I3,I4
      If((I1+I2.Ne.0).And.(I3+I4.Eq.0)) Then
      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
      Gamma(Ind)=X
      EndIf
      GoTo 20
  30  Close(10)
      Call CpySym(AUXM,Gamma,NBasis)
      Call Diag8(AUXM,NBasis,NBasis,PC,Work)
      Call SortOcc(PC,AUXM,NBasis)
C
      Sum=Zero
      NAc=0
      Do I=1,NBasis
      Sum=Sum+PC(I)
      If(PC(I).Gt.Zero) NAc=NAc+1
      EndDo
C
      NInAc=NELE-Sum+1.D-11
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
C     READ XOne (SAVE IN XKin), TRANSFORM TO NO'S 
      Write(6,'(" Reading in one-ele integrals ...",/)')
      K=1
      FName(K:K+11)='intcoul.dat'
      Call GetENuc(ENuc,FName,NBasis)
      Call Int1(XKin,XNuc,NInte1,FName,Nbasis)
C
      Call MatTr(XKin,URe,NBasis)
C
C     READ 2-EL INTEGRALS AND TRANSFORM TO NO's
      TwoEl(1:NInte2)=Zero
      Write(6,'(" Reading in two-ele integrals ...",/)')
      Call Int2(TwoEl,FName,NInte2,NBasis)
C
      Write(6,'(" Transforming two-electron integrals ...",/)')
      Call TwoNO1(TwoEl,URe,NBasis,NInte2)
C herer!!! ???
c      Write(6,'(" Skipping reading 2-RDM from rdmdump.dat.",/)')
c      Write(6,'(" 2-RDM will be read from rdm2.dat file ...",/)')
c      GoTo 888
C
C     READ ACTIVE 2-RDM AND TRANSFORM TO NO'S  
C
      Write(6,'(" Reading in 2-RDM ...")')
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      Allocate (RDM2(NRDM2))
      RDM2(1:NRDM2)=Zero
      Allocate (RDMAB2(NRDM2))
      RDMAB2(1:NRDM2)=Zero
C
      IAAAA=1
      IABBA=0
      Open(10,File='rdmdump.dat')
   22 Read(10,*,End=33) X,I1,I2,I3,I4
C
      I=I1+NInAc
      J=I2+NInAc
      K=I3+NInAc
      L=I4+NInAc
C
      If(IABBA.Eq.1.And.I1+I2+I3+I4.Eq.0) GoTo 33
      If(IAAAA.Eq.1.And.I1+I2+I3+I4.Eq.0) Then
      IAAAA=0
      GoTo 22
      EndIf
      If(IAAAA.Eq.0.And.IABBA.Eq.0.And.I1+I2+I3+I4.Eq.0) Then
      IABBA=1
      GoTo 22
      EndIf 
C
      If(IAAAA.Eq.1) Then
      RDM2(NAddrRDM(L,K,I,J,NBasis))=X
      RDM2(NAddrRDM(K,L,I,J,NBasis))=-X 
      RDM2(NAddrRDM(L,K,J,I,NBasis))=-X
      RDM2(NAddrRDM(K,L,J,I,NBasis))=X
      EndIf
C
      If(IABBA.Eq.1) Then
      RDMAB2(NAddrRDM(L,K,I,J,NBasis))=X
      RDMAB2(NAddrRDM(K,L,J,I,NBasis))=X
      EndIf
C
      GoTo 22
   33 Close(10)
C
      Do I=1,NRDM2 
      RDM2(I)=RDM2(I)+RDMAB2(I)
      EndDo
C 
      Call TrRDM2(RDM2,URe,NBasis,NRDM2)
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
C     COMPUTE THE ENERGY FOR CHECKING
C
C     THE INACTIVE PART INActive
      ETot=Zero
      Do I=1,INActive
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*xkin(II)
      EndDo
C
      Do IP=1,INActive
      Do IQ=1,INActive
      Do IR=1,INActive
      Do IS=1,INActive
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      ETot=ETot+RDM2(IAdd)*Twoel(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
      Write(*,*)'Inactive (including core) + ENuclear',ETot+ENuc
C
      ETot=Zero
      eact=Zero
      Do I=1,NOccup
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*XKin(II)
C
      If(Occ(I).Ne.1.D0) Then
      sum=zero
      Do J=1,INActive
      sum=sum+
     $ 2.D0*Twoel(NAddr3(I,I,J,J))-Twoel(NAddr3(I,J,I,J))
      EndDo
      eact=eact+Two*Occ(I)*(XKin(II)+sum)
      EndIf
C
      EndDo
      Write(*,*)'Active One-Electron Energy',eact
C
      etot2=zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      ETot=ETot+RDM2(IAdd)*Twoel(NAddr3(IP,IR,IQ,IS))
      if(occ(ip).ne.1.d0.and.occ(ir).ne.1.d0.and.occ(iq).ne.1.d0.and.
     $ occ(is).ne.1.d0) then
      etot2=etot2+RDM2(IAdd)*Twoel(NAddr3(IP,IR,IQ,IS))
      endif
      EndDo
      EndDo
      EndDo
      EndDo
      Write(*,*)'Active Two-Electron Energy',etot2
      write(*,*)'Total DMRG Energy',etot+enuc
c *****************************************

C
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
      Deallocate(RDMAB2)

c herer!!! ???
c      stop
  888 Continue
C
C     INTEGRALS ARE TRANSFORMED SO URe IS SET AS A UNIT MATRIX 
C
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=Zero
      If(I.Eq.J) URe(I,J)=One
      EndDo
      EndDo
C
      Return
      End

*Deck LdInteg
      Subroutine LdInteg(Title,XKin,XNuc,ENuc,Occ,URe,DipX,DipY,DipZ,
     $ TwoEl,TwoElErf,UAOMO,NSymMO,NInte1,NBasis,NInte2)
C
C     READ/WRITE THE ONE- AND TWO-ELECTRON INTEGRALS 
C     INTERFACED WITH MOLPRO (INTEGRALS ARE READ FROM FCIDUMP FILES)
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 XKin(*),XNuc(*),DipX(*),DipY(*),DipZ(*),TwoEl(*),
     $ TwoElErf(*)
      Dimension NSymMO(NBasis)
C
      Character*60 Title,FName,FMultTab
C
      Include 'commons.inc'
C
C     LOCAL ARRAYS
C
      Dimension MultpC(15,15),NumOSym(15),UAOMO(NBasis,NBasis)
C
C     DETERMINE A NUMBER OF ORBITALS OF EACH SYMMETRY
C
      Call DetSym(MultpC,NSymMO,NumOSym,NBasis)
C
C     READ THE TRANSFORMATION MATRIX FROM AO TO MO
C
      Call GetUAOMO(UAOMO,NSymMO,NumOSym,NBasis)
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=1
C
C     IF INTEGRALS ARE IN MO REPRESENTATION THEN
C
      If(IAO.Eq.0) Then
C
      FName(K:K+11)='intcoul.dat'
      Call GetENuc(ENuc,FName,NBasis)
C
      FName(K:K+11)='intcoul.dat'
      Call Int1(XKin,XNuc,NInte1,FName,Nbasis)
C
      FName(K:K+11)='intcoul.dat'
      Call Int2(TwoEl,FName,NInte2,NBasis) 
      FName(K:K+11)='interf.dat'
C
C     LOAD ERF INTEGRALS ONLY IF REQUIRED
C
      If(IFun.Ne.0.And.IFunSR.Ne.0) 
     $ Call Int2(TwoElErf,FName,NInte2,NBasis)   
C
C     Else of (IAO.Eq.0)
      Else
C
      Call GetENuc_AO(ENuc,Title)
C
C     DETERMINE A NUMBER OF ORBITALS OF EACH SYMMETRY
C
      Call DetSym(MultpC,NSymMO,NumOSym,NBasis)
C
C     READ THE TRANSFORMATION MATRIX FROM AO TO MO
C
      Call GetUAOMO(UAOMO,NSymMO,NumOSym,NBasis)
C
      FName(K:K+8)='ekin.dat'
      Call Int1_AO(XKin,NInte1,FName,UAOMO,NSymMO,NumOSym,NBasis)
C
      FName(K:K+8)='epot.dat'
      Call Int1_AO(XNuc,NInte1,FName,UAOMO,NSymMO,NumOSym,NBasis)
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
C
      FName(K:K+10)='.reg.integ'
      Call Int2_AO(TwoEl,UAOMO,NSymMO,NumOSym,MultpC,
     $ FName,NInte1,NInte2,NBasis)
C
      FName(K:K+10)='.erf.integ'
      Call Int2_AO(TwoElErf,UAOMO,NSymMO,NumOSym,MultpC,
     $ FName,NInte1,NInte2,NBasis)
C
      EndIf
C
      Return
      End

*Deck DimSym
      Subroutine DimSym(NBasis,NInte1,NInte2,MxHVec,MaxXV)
C 
      Implicit Real*8 (A-H,O-Z)
C
C     CALCULATES THE ARRAY DIMENSIONS
C     ASSIGNS THE SYMMETRIES TO BASIS FUNCTIONS
C
C     NBasis - maximal index of the integrals (the number of basis functions)
C     NInte1 - dimension of the one-electron integrals array
C     NInte2 - dimension of the two-electron integrals array
C     MxHVec - dimension of the total hessian matrix
C     MaxXV  - the number of the independent rotations (X parameters) 
C
      NInte1=NBasis*(NBasis+1)/2 
C     
      NInte2=NInte1*(NInte1+1)/2
C
      MaxXV=NBasis*(NBasis-1)/2
C
      MxHVec=MaxXV+NBasis
C
      Return
      End

*Deck RWInput
      Subroutine RWInput(Title,ZNucl,Charge,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FMultTab
C
      Include 'commons.inc'
C
C     READ THE JOB TITLE
C
      Open(10,File='qudmft.in',Status='Old')
C
      Read(10,'(A60)') Title
C
C     READ THE KEYWORDS
C
      Call RInput(ZNucl,Charge,IPrint,NBasis,NELE,XELE)
C
      If(IPrint.Ge.0) Then
C
C     PRINT THE INPUT DATA
C
      Write(6, '('' **********************************************'',
     $            ''*****************'')')
C
      Write(6,'(1X,''* '',A60,''*'')') Title
      Write(6,'('' *'',61X,''*'')')
C
      Write(6,
     $ '('' * NUCLEAR CHARGE = '',F9.6,''   CHARGE = '',F9.6,13X,
     $ ''*'')') ZNucl,Charge
C
      Write(6,'('' *'',61X,''*'')')
      Write(6,'('' * NO. OF BASIS SET FUNCTIONS: '',I5,27X,''*'')')
     $ NBasis
C
      Write(6,'('' *'',61X,''*'')')
C
      If(IRes.Ne.1) Then
C
      Write(6,'('' * DENSITY MATRIX FUNCTIONAL CALCULATION'',
     $      23X,''*'')')
      Write(6,
     $ '('' * Functional : '',I3 , 44X,''*'')') IFun
C
      Write(6, '('' **********************************************'',
     $            ''*****************'')')
C
      EndIf
C
      EndIf
C
      Return
      End 

*Deck RInput
      Subroutine RInput(ZNucl,Charge,IPrint,NBasis,NELE,XELE)
C 
      Implicit Real*8 (A-H,O-Z)
C
C     READS THE INPUT KEYWORDS 
C
      Character*80 Line1
      Character*160 Line
      Character Digit(10)
C
      Data Digit/'0','1','2','3','4','5','6','7','8','9'/
C
      Parameter(Zero=0.0D0,Half=0.5D0,Two=2.0D0,Ten=10.0D0)
C
C     READ THE INPUT FILE
C
      Do 20 NLines=0,1
      Read(10,'(A80)',End=15) Line1
C
      Do 20 I=1,80
   20 Line(I+NLines*80:I+NLines*80)=Line1(I:I) 
C
   15 Close(10,Status='Keep')
C
C     *****ZNucl*****
C
      ZNucl=Zero
      K=0
  250 K=K+1
      If (K.Eq.160) Stop 'FATAL ERROR: ZNucl ENTRY MISSING!'
      If (Line(K:K+4).Ne.'znucl') GoTo 250
C
      IFlag=0
      NDig1=0
      NDig2=0
      ICount=0
  245 ICount=ICount+1
      If (Line(K+5+ICount:K+5+ICount).Eq.'.') Then
      IFlag=1
      NDig1=ICount-1
      GoTo 245
      ElseIf (Line(K+5+ICount:K+5+ICount).Eq.' ') Then
      NDig2=ICount-2-NDig1
      Else
      GoTo 245
      EndIf
      If (NDig1.Eq.0) NDig1=ICount-1
      If(IFlag.Eq.0) NDig2=Zero
C
      Do 260 J=1,NDig1
      I=-1
  270 I=I+1
      If (I.Gt.10) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY ZNucl IN THE INPUT FILE!'
      If (Line(K+5+J:K+5+J).Ne.Digit(I+1)) GoTo 270    
  260 ZNucl=ZNucl+Float(I)*Ten**(NDig1-J)
C
      Do 280 J=1,NDig2
      I=-1
  290 I=I+1
      If (I.Gt.10) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY ZNucl IN THE INPUT FILE!'
      If (Line(K+6+NDig1+J:K+NDig1+6+J).Ne.Digit(I+1)) GoTo 290    
  280 ZNucl=ZNucl+Float(I)/Ten**J
C
C     *****Charge*****
C
  295 Charge=Zero
      K=0
  300 K=K+1
      If (K.Eq.160) Stop 'FATAL ERROR: Charge ENTRY MISSING!'
      If (Line(K:K+5).Ne.'charge') GoTo 300
C
      ISg=0
      If (Line(K+7:K+7).Eq.'+') ISg=1
      If (Line(K+7:K+7).Eq.'-') ISg=-1
      K=K+Abs(ISg)
      If (ISg.Eq.0) ISg=1
C
      NDig1=0
      NDig2=0
      ICount=0
  310 ICount=ICount+1
C
      If (Line(K+6+ICount:K+6+ICount).Eq.'.') Then
      NDig1=ICount-1
      GoTo 310
      ElseIf (Line(K+6+ICount:K+6+ICount).Eq.' ') Then
      NDig2=ICount-2-NDig1
      Else
      GoTo 310
      EndIf
      If (NDig1.Eq.0) NDig1=ICount-1
C
      Do 320 J=1,NDig1
      I=-1
  330 I=I+1
      If (I.Gt.10) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY Charge IN THE INPUT FILE!'
      If (Line(K+6+J:K+6+J).Ne.Digit(I+1)) GoTo 330
  320 Charge=Charge+Float(I)*Ten**(NDig1-J)
C
      Do 340 J=1,NDig2
      I=-1
  350 I=I+1
      If (I.Gt.10) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY Charge IN THE INPUT FILE!'
      If (Line(K+7+NDig1+J:K+NDig1+7+J).Ne.Digit(I+1)) GoTo 350    
  340 Charge=Charge+Float(I)/Ten**J
C
      Charge=Float(ISg)*Charge
C
C     *****NBasis*****
C
  495 NBasis=0
      K=0
  400 K=K+1
      If (K.Eq.160) Stop 'FATAL ERROR: NBasis ENTRY MISSING!'
      If (Line(K:K+5).Ne.'nbasis') GoTo 400
C
      NDig1=0
      NDig2=0
      ICount=0
  410 ICount=ICount+1
C
      If (Line(K+6+ICount:K+6+ICount).Eq.' ') Then
      NDig2=ICount-2-NDig1
      Else
      GoTo 410
      EndIf
      If (NDig1.Eq.0) NDig1=ICount-1
C
      Do 420 J=1,NDig1
      I=-1
  430 I=I+1
      If (I.Gt.10) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY NBasis IN THE INPUT FILE!'
      If (Line(K+6+J:K+6+J).Ne.Digit(I+1)) GoTo 430
  420 NBasis=NBasis+Float(I)*Ten**(NDig1-J)
C
C     *****IPrint*****
C
      K=0
  360 K=K+1
      If (K.Eq.160) Stop 'FATAL ERROR: IPrint ENTRY MISSING!'
      If (Line(K:K+5).Ne.'iprint') GoTo 360
C
      IPrint=-1
      If (Line(K+7:K+7).Eq.'0') IPrint=0
      If (Line(K+7:K+7).Eq.'1') IPrint=1
      If (Line(K+7:K+7).Eq.'2') IPrint=2
      If (Line(K+7:K+7).Eq.'3') IPrint=3
      If (Line(K+7:K+7).Eq.'4') IPrint=4
      If (IPrint.Eq.-1) 
     $ Stop 'FATAL ERROR: INCORRECT ENTRY IPrint IN THE INPUT FILE!' 
C
      XELE=Half*(ZNucl-Charge)
      NELE=Int(ZNucl-Charge+Half)/2
C
      Return
      End

*Deck Int1
      Subroutine Int1(XKin,XNuc,NInte1,FName,NBasis)
C
C     READS ONE-ELECTRON INTEGRALS FROM THE MOLPRO FCIDUMP FILE
C     AND STORE THEM IN XKin (XNuc IS EMPTY!!!)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.D0)
C
      Character*60 FName,Aux1
C
      Dimension XKin(NInte1),XNuc(NInte1)
C
      Do I=1,NInte1
      XKin(I)=Zero
      XNuc(I)=Zero
      EndDo
C
C     READ ONE-ELECTRON INTEGRALS IN THE MO REPR      
C
      Open(10,File=FName)
C
    2 Read(10,'(A10)')Aux1
c      If(Aux1(1:5).Eq." &END") GoTo 20
      If(Aux1(1:6).Eq."  ISYM".Or.Aux1(1:5).Eq." ISYM") Then
      Read(10,'(A10)')Aux1
      GoTo 20
      EndIf
      GoTo 2
C 
   20 Continue
      Read(10,*,End=30) X,I1,I2,I3,I4

      If((I1+I2.Ne.0).And.(I3+I4.Eq.0)) Then
      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
      XKin(Ind)=X
      EndIf 
      GoTo 20   
C
  30  Close(10)      
C
      Return
      End

*Deck Int2
      Subroutine Int2(TwoEl,FName,NInte2,NBasis)
C
C     READS TWO-ELECTRON INTEGRALS FROM THE MOLPRO FCIDUMP FILE 
C     (CHEMICAL NOTATION IS USED!)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FName,Aux1
      Dimension TwoEl(NInte2)
C
      Open(10,File=FName)
C
    2 Read(10,'(A10)')Aux1
c      If(Aux1(1:5).Eq." &END") GoTo 20
      If(Aux1(1:6).Eq."  ISYM".Or.Aux1(1:5).Eq." ISYM") Then
      Read(10,'(A10)')Aux1
      GoTo 20
      EndIf
      GoTo 2
C 
   20 Continue
      Read(10,*,End=30) X,I1,I2,I3,I4 
      If(I3+I4.Ne.0) TwoEl(NAddr3(I1,I2,I3,I4))=X
      GoTo 20
C
  30  Close(10)
C 
      Return
      End  

*Deck GetENuc_AO
      Subroutine GetENuc_AO(ENuc,Title)
C
      Implicit Real*8 (A-H,O-Z)
      Character*60 Title,FName
      Character*100 Line,Aux,Aux1,Aux2
C
C     GET THE NUCLEAR REPULSION ENERGY FROM THE MOLPRO OUTPUT
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
      FName(K:K+4)='.out'
C
      Line="grep 'NUCLEAR REPULSION ENERGY' "//FName(1:K+4)//" >tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Read(10,*)Aux,Aux1,Aux2,ENuc
      Close(10)
C
      Return
      End

*Deck Int1_AO
      Subroutine Int1_AO(XOne,NInte1,FName,UAOMO,NSymMO,NumOSym,NBasis)
C
C     READS ONE-ELECTRON INTEGRALS IN AO FROM THE FName FILE
C     AND STORE THEM IN XOne
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.D0)
      Include 'commons.inc'
C
      Character*60 FName,Aux1
C
      Dimension XOne(NInte1),UAOMO(NBasis,NBasis),NSymMO(NBasis),
     $ NumOSym(15)
C
C
      Do I=1,NInte1
      XOne(I)=Zero
      EndDo 
C
      Open(10,File=FName)
C
      Read(10,'(A10)')Aux1
      Read(10,'(A10)')Aux1
C
      Ind=0  
      Do ISym=1,MxSym
C
      Do I=1,NumOSym(ISym)
      Ind1=Ind+I
      Ind2=Ind+K
      Index=(Max(Ind1,Ind2)*(Max(Ind1,Ind2)-1))/2+Min(Ind1,Ind2)
      Read(10,*) (XOne(IndSym(Ind+I,Ind+K)),K=1,NumOSym(ISym))
      EndDo
C
      Ind=Ind+NumOSym(ISym)
C   
      EndDo
      Close(10)
C
C     TRANSFORM XOne FROM AO TO MO
C
      Call MatTr(XOne,UAOMO,NBasis)
C
      Return
      End

*Deck Int2_AO
      Subroutine Int2_AO(TwoEl,UAOMO,NSymMO,NumOSym,MultpC,
     $ FName,NInte1,NInte2,NBasis)
C
C     READS TWO-ELECTRON INTEGRALS IN AO
C
      Implicit Real*8 (A-H,O-Z)
      Include 'commons.inc'
C
      Parameter (Zero=0.D0)
C
      Character*60 FName
C
      Dimension TwoEl(NInte2),UAOMO(NBasis,NBasis),NSymMO(NBasis),
     $ NumOSym(15),MultpC(15,15),Record(NInte1),NTB(8)
C
      Real*8, Allocatable :: TwoAO(:)
C
      Write(6,'(" Loading two-electron integrals in AO ...",/)')
C
      Allocate (TwoAO(NInte2))
C
      TwoAO(1:NInte2)=         0.D0
C
      Open(10,File=FName,Form='Unformatted')
      Do III=1,8
      Read(10)NTB(III)
      EndDo
C
      Do NSymPQ=1,MxSym
C
      Do NSymP=1,MxSym
C
      NSymQ=MultpC(NSymPQ,NSymP)
      If(NSymQ.Gt.NSymP) Cycle 
C
      Do IP=1,NumOSym(NSymP)
C
      If(NSymP.Eq.NSymQ) Then
      NEndQ=IP
      Else
      NEndQ=NumOSym(NSymQ)
      EndIf
C
      Do IQ=1,NEndQ
C
C     READ A RECORD
C 
      Read(10) NRecL
      Read(10) (Record(I),I=1,NRecL)
      ICounter=0
C
C     LOOP OVER IR, IS
C
      NSymRS=NSymPQ     
      Do NSymR=1,MxSym
C
      NSymS=MultpC(NSymRS,NSymR)
      If(NSymS.Gt.NSymR) Cycle
C
      If(NSymR.Eq.NSymS) Then
C
      Do IR=1,NumOSym(NSymR)
      Do IS=1,IR
      ICounter=ICounter+1
      IPP=IndOrb(IP,NSymP,NumOSym)
      IQQ=IndOrb(IQ,NSymQ,NumOSym)
      IRR=IndOrb(IR,NSymR,NumOSym)
      ISS=IndOrb(IS,NSymS,NumOSym)
      NAdd=NAddr3(IPP,IQQ,IRR,ISS)
      TwoAO(NAdd)=Record(ICounter)
      EndDo
      EndDo
C
      Else
C
      NBlock=NTB(NSymR)*NTB(NSymS) 
      ICountRS=0
C
      Do IS=1,NumOSym(NSymS)
      Do IR=1,NumOSym(NSymR)
      ICounter=ICounter+1
      ICountRS=ICountRS+1
      IPP=IndOrb(IP,NSymP,NumOSym)
      IQQ=IndOrb(IQ,NSymQ,NumOSym)
      IRR=IndOrb(IR,NSymR,NumOSym)
      ISS=IndOrb(IS,NSymS,NumOSym)
      NAdd=NAddr3(IPP,IQQ,IRR,ISS)
      TwoAO(NAdd)=Record(ICounter)
      EndDo
      EndDo
C
C     SHIFT THE COUNTER TO SKIP POSSIBLE ZEROS
C
      If(NTB(NSymR).Gt.NumOSym(NSymR)) ICounter=ICounter+NBlock-ICountRS
C
      EndIf
C
c     enddo NSymR
      EndDo
      If(NRecL-ICounter.Gt.10) 
     $ Stop 'Fatal error in initia.f. Check MxSym'
c     enddo IQ
      EndDo
c     enddo IP
      EndDo     
c     enddo symp
      EndDo
c     enddo sympq
      EndDo
C
      Close(10)
C
C     TRANSFORM THE INTEGRALS TO MO
C
      Write(6,'(" Transforming two-electron integrals to MO ...",/)')
      Call TwoNO(TwoEl,UAOMO,TwoAO,NBasis,NInte2)
C
      Deallocate (TwoAO)
C
      Return
C
      End

*Deck GetENuc 
      Subroutine GetENuc(ENuc,FName,NBasis)
C
C     READS THE NUCLEAR ENERGY FROM THE MOLPRO FCIDUMP FILE
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.D0)
C
      Character*60 FName,Aux1
C
      Open(10,File=FName)
    2 Read(10,'(A10)')Aux1
c      If(Aux1(1:5).Eq." &END") GoTo 20
      If(Aux1(1:6).Eq."  ISYM".Or.Aux1(1:5).Eq." ISYM") Then
      Read(10,'(A10)')Aux1
      GoTo 20
      EndIf
      GoTo 2
C 
   20 Continue
      Read(10,*,End=30) X,I1,I2,I3,I4 
      If(I1+I2+I3+I4.Eq.0) ENuc=X
      GoTo 20
C
  30  Close(10)
C
      Return
      End

*Deck DetSym
      Subroutine DetSym(MultpC,NSymMO,NumOSym,NBasis)
C
C     DETRMINE THE NUMBER OF ORBITALS IN EACH SYMMETRY 
C     AND READ THE MULTIPLICATION TABLE
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter (Zero=0.D0)
C
      Dimension MultpC(15,15),NSymMO(NBasis),NumOSym(15)
C
C     READ THE DIRECT MULTIPLICATION TABLE
C
      Open(10,File=FMultTab)
      Do I=1,MxSym
      Read(10,*)(MultpC(I,J),J=1,I)
      Do J=1,I
      MultpC(J,I)=MultpC(I,J)
      EndDo
      EndDo
      Close(10)
C
      Do I=1,MxSym
C
      NumOSym(I)=0
C
      Do J=1,NBasis
      If(NSymMO(J).Eq.I) NumOSym(I)=NumOSym(I)+1
      EndDo
C
      EndDo
C
      Return 
      End

*Deck GetUAOMO
      Subroutine GetUAOMO(UAOMO,NSymMO,NumOSym,NBasis)
C
C     LOAD A TRANSFORMATION MATRIX AO->MO FROM orb_hf.dat
C
      Implicit Real*8 (A-H,O-Z)
      Parameter (Zero=0.D0)
C
      Character*100 Line
      Character*5 Aux
C
      Include 'commons.inc'
C
      Dimension UAOMO(NBasis,NBasis),NSymMO(NBasis),NumOSym(15)
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAOMO(I,J)=Zero
      EndDo
      EndDo
C
      Open(123,File="orb_hf.dat")
C
      Ind=0
      Do ISym=1,MxSym
C
      Do I=1,NumOSym(ISym)
      Read(123,*) (UAOMO(Ind+I,Ind+K),K=1,NumOSym(ISym))
      EndDo
C
      Ind=Ind+NumOSym(ISym)
C
      EndDo
      Close(123)
C
C     IF ILoc=1 READ UAOMO FROM molden FILE
C
      If(ILoc.Eq.1) Then
C
      Write(6,'(/,X,"LOCALIZED MO ORBITALS USED AS GUESS")')
C
      If(MxSym.Ne.1) 
     $ Stop 'Fatal Error. Localized orbitals not available'
C
      Open(123,File="local.molden")
C
   20 Read(123,'(A4)',End=10) Line(1:4)
C
      If(Line(1:4).Eq."[MO]") Then
C
      Do I=1,NBasis
      Read(123,'(A5,F9.1)')Aux(1:5),XX
      IOrb=XX
      Read(123,*)
      Read(123,*)
      Read(123,*)
      Do J=1,NBasis
      Read(123,*)X,Y
      UAOMO(IOrb,J)=Y
      EndDo
      EndDo
C
      EndIf
C
      GoTo 20
C
   10 Continue
      Close(123)
C
      EndIf
C
      Return
      End

*Deck IndSym
      Function IndSym(I,J)
      IndSym=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      Return
      End

*Deck IndOrb
      Function IndOrb(I,NSymI,NumOSym)
C
C     RETURNS THE INDEX OF AN ORBITALS WHICH IS THE ITH ORBITALS OF THE SYMMETRY NSymI
C 
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Dimension NumOSym(15)
C
      ISum=0
      Do K=1,NSymI-1
      ISum=ISum+NumOSym(K)
      EndDo
      IndOrb=I+ISum
C
      Return
      End


****** PIOTR KOWALSKI, 01/2018 ***********************
      subroutine read1el(XKin,UMOAO,NBasis,NInte1)
C
C     Reads 1-el integrals in AO and ttransform to NO
C     Returns XKin in NO
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension XKin(NInte1),UMOAO(NBasis,NBasis)
C
      parameter (mxbuf = 10000)  ! KP
      double precision  dbuf(mxbuf)
      integer ibuf(mxbuf*2)
      integer iunit77,iunit88, iunit99, ndim, norb, nbas, nfone
      logical iprtvc
      integer  maxrep, naos, lbuf, nibuf, nbits
      common /daltwoel/  maxrep, naos(8), lbuf, nibuf, nbits
      integer nxx
C
      XKin(1:NInte1)=         0.D0
C      
      iunit77=77
      iunit88=88
      iunit99=99

      OPEN(UNIT=iunit77,FILE='AOONEINT',STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      call initoneel(iunit77)

10    continue
      call readoneel(iunit77, dbuf, ibuf, lbuf, nibuf, nxx)

C       write(*,*) nxx
C       write(*,*) "-----------"
C       write(*,*) dbuf
C       write(*,*) ibuf

      do i=1,nxx
        XKin(ibuf(i)) = dbuf(i)
      end do

      if (nxx .ge. 0) go to 10
C
C     TRANSFORM THE INTEGRALS
C
      Write(6,'(" Transforming One-electron integrals ...",/)')
        write(55,*) XKin
      Call MatTr(XKin,UMOAO,NBasis)
        write(56,*) XKin

      return
      end

      subroutine readoneel(iunit, dbuf, ibuf, lbuf, nibuf, nxx)
      integer iunit, ibuf, lbuf, nibuf, nxx
      double precision dbuf
      dimension dbuf(lbuf), ibuf(lbuf, nibuf)

      read (iunit) dbuf, ibuf, nxx

      end

      subroutine initoneel(iunit)
      integer iunit
      character*8 ONEHAMIL
      data ONEHAMIL  /'ONEHAMIL'/

      integer  maxrep, naos, lbuf, nibuf, nbits
      common /daltwoel/  maxrep, naos(8), lbuf, nibuf, nbits
      logical findlab

      integer i

      nibuf=1
      lbuf=600

      if (.not. findlab(ONEHAMIL ,iunit)) then
        write(*,*) 'Error finding label ', ONEHAMIL
        stop
C      else
C        write(*,*) 'label ', ONEHAMIL, " was founded" 
C        stop
      end if
      end

****** SUBROUTINES FROM RAFAL PODESZWA, 03/2017 ******
*
*Deck read2el 
      subroutine read2el(TwoEl,UMOAO,NBasis,NInte2)
C
C     Reads 2-el integrals in AO and ttransform to NO
C     Returns TwoEl in NO
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension TwoEl(NInte2),UMOAO(NBasis,NBasis)
C
      parameter (mxbuf = 10000)  ! KP
      double precision  dbuf(mxbuf)
      integer ibuf(mxbuf*2)
      integer iunit77,iunit88, iunit99, ndim, norb, nbas, nfone
      logical iprtvc
      integer  maxrep, naos, lbuf, nibuf, nbits
      common /daltwoel/  maxrep, naos(8), lbuf, nibuf, nbits
      integer nxx
C
      TwoEl(1:NInte2)=         0.D0
C      
      iunit77=77
      iunit88=88
      iunit99=99
      
      OPEN(UNIT=iunit77,FILE='AOTWOINT',STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     
      call inittwoel(iunit77)
      if (nibuf*lbuf .gt. mxbuf) then
        write(*,*) 'Dalton buffer greater than mxbuf. Stop'
        stop
      end if
      
10    continue
      call readtwoel(iunit77, dbuf, ibuf, lbuf, nibuf, nxx)
      
c      write (*,*) nxx, ' two electron integrals read'
      do i=1,nxx
        call unpckdlt(ibuf(i), ibuf(i+lbuf), nibuf, nbits,ip,iq,ir,is)
c        if (dabs(dbuf(i)).gt.1e-8)
c     *          write(*,'(4I4,F12.8)') ip, iq, ir, is, dbuf(i)
      TwoEl(NAddr3(ip,iq,ir,is))=dbuf(i)
      end do
      
      if (nxx .ge. 0) go to 10
C
C     TRANSFORM THE INTEGRALS
C
      Write(6,'(" Transforming two-electron integrals ...",/)')
      Call TwoNO1(TwoEl,UMOAO,NBasis,NInte2)
C
      return  
      end
      
      
      subroutine unpckdlt(ibuf, ibuf2, nibuf, nbits, ip, iq, ir, is)
      integer ibuf, ibuf2, nibuf, nbits, ip, iq, ir, is
      
      
      IF (NIBUF .EQ. 1) THEN
      IF (NBITS .EQ. 8) THEN
            LABEL = ibuf
            ip = IAND(ISHFT(ibuf,-24),255)
            iq = IAND(ISHFT(ibuf,-16),255)
            ir = IAND(ISHFT(ibuf, -8),255)
            is = IAND(       ibuf,    255)
c#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
c      ELSE IF (NBITS .EQ. 16) THEN
c            LABEL = ibuf
c            ip = AND(SHIFTR(ibuf,48),65535)
c            iq = AND(SHIFTR(ibuf,32),65535)
c            ir = AND(SHIFTR(ibuf,16),65535)
c            is = AND(       ibuf,    65535)
c#endif
      ELSE
         write(*,*) 'UnPack Dalton error'
         stop
      END IF
      ELSE
            ip = IAND(ISHFT(ibuf,-16),65535)
            iq = IAND(       ibuf    ,65535)
            ir = IAND(ISHFT(ibuf2,-16),65535)
            is = IAND(       ibuf2    ,65535)
      END IF

      end
     
      subroutine readtwoel(iunit, dbuf, ibuf, lbuf, nibuf, nxx)
      integer iunit, ibuf, lbuf, nibuf, nxx
      double precision dbuf
      dimension dbuf(lbuf), ibuf(lbuf, nibuf)
      
      read (iunit) dbuf, ibuf, nxx
      
      end
       
      subroutine inittwoel(iunit)
      integer iunit
      character*8 BASINFO, BASTWOEL
      data BASINFO /'BASINFO'/, BASTWOEL /'BASTWOEL'/
      
      integer  maxrep, naos, lbuf, nibuf, nbits
      common /daltwoel/  maxrep, naos(8), lbuf, nibuf, nbits
      
      logical findlab
      
      integer i
      
      
      if (.not. findlab(BASINFO,iunit)) then
        write(*,*) 'Error finding label', BASINFO
        stop
      end if
      read (iunit) maxrep,(naos(i),i=1,8),lbuf,nibuf,nbits
      write(*,*) 'Dalton two-el. file initialized'
      write(*,*) 'Buffer size: ',lbuf, ', integers per index packet: '
     &,nibuf, ', bits: ', nbits
      if (.not. findlab(BASTWOEL,iunit)) then
        write(*,*) 'Error finding label', BASTWOEL
        stop
      end if
      end
    
      
      
      logical function findlab(label, labelunit)
      character*8 label, stars, b(4)
      parameter (stars = '********')
      
      rewind(labelunit)
10    read (labelunit, END=100, ERR=50) b
      if (b(1).ne.stars) go to 10
      if (b(4).ne.label) go to 10
      findlab=.true.
      return
50    continue
100   continue
      findlab=.false.
      return
      end

*Deck TrRDM2
      Subroutine TrRDM2(RDM2,URe,NBasis,NRDM2)
C
C     TRANSFORM RDM2 WITH URe
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0)
C
      Dimension URe(NBasis,NBasis),RDM2(NRDM2)
C
C     LOCAL ARRAYS
C 
      Dimension Aux1(NBasis,NBasis,NBasis,NBasis),
     $  Aux2(NBasis,NBasis,NBasis,NBasis)
C
      Write(6,'(/,X,"FCI RDM2 TRANSFORMATION TO NO IN PROCESS...")')
C
C     COPY RDM2 TO Aux1
C
      Do I=1,NBasis
      Do J=1,NBasis
      Do K=1,NBasis
      Do L=1,NBasis
      Aux1(I,J,K,L)=RDM2(NAddrRDM(I,J,K,L,NBasis))
      EndDo
      EndDo
      EndDo
      EndDo
C
C     FIRST INDEX
C
      Do IA=1,NBasis
      Do J=1,NBasis
      Do K=1,NBasis
      Do L=1,NBasis
C
      Aux2(IA,J,K,L)=Zero
      Do I=1,NBasis
      Aux2(IA,J,K,L)=Aux2(IA,J,K,L)+URe(IA,I)*Aux1(I,J,K,L)
      EndDo
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     SECOND INDEX
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      Do K=1,NBasis
      Do L=1,NBasis
C
      Aux1(IA,IB,K,L)=Zero
      Do J=1,NBasis
      Aux1(IA,IB,K,L)=Aux1(IA,IB,K,L)+URe(IB,J)*Aux2(IA,J,K,L)
      EndDo
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     THIRD INDEX
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      Do IC=1,NBasis
      Do L=1,NBasis
C
      Aux2(IA,IB,IC,L)=Zero
      Do K=1,NBasis
      Aux2(IA,IB,IC,L)=Aux2(IA,IB,IC,L)+URe(IC,K)*Aux1(IA,IB,K,L)
      EndDo
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     FOURTH INDEX
C     
      Do IA=1,NBasis
      Do IB=1,NBasis
      Do IC=1,NBasis
      Do ID=1,NBasis
C
      RDM2(NAddrRDM(IA,IB,IC,ID,NBasis))=Zero
      Do L=1,NBasis
      RDM2(NAddrRDM(IA,IB,IC,ID,NBasis))=
     $ RDM2(NAddrRDM(IA,IB,IC,ID,NBasis))+URe(ID,L)*Aux2(IA,IB,IC,L)
      EndDo
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,"DONE WITH FCI RDM2 TRANSFORMATION")')
C
      Return
      End
