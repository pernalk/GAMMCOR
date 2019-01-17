*Deck ReadDAL
      Subroutine ReadDAL(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem,Flags)
C
C     READ HAO, 2-EL INTEGRALS IN NO, C_COEFFICIENTS, IGEM FROM A DALTON_GENERATED FILE
C     READ UMOAO FROM DALTON.MOPUN
C
      use types   
      use sorter
      use tran   
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 XKin(NInte1),XNuc(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),UMOAO(NBasis,NBasis)
      double precision, allocatable :: TMPMO(:,:)
C
      Character*60 Line
      Character*30 Line1
C
      type(FlagsData) :: Flags
C
      Include 'commons.inc'
C
C     ICASSCF=1 - read occ from occupations.dat
      ICASSCF=Flags%ICASSCF
C
C     IDMRG - integrals and RDM's read from external dmrg files
      IDMRG=Flags%IDMRG
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
C      Call readtwoint(NBasis,'AOTWOINT','AOTWOSORT')
      Call read2el(TwoEl,UMOAO,NBasis,NInte2)
!      Call altread2el(TwoEl,UMOAO,NBasis,NInte2)

C
C     TESTING TRAN4
C      allocate(TMPMO(NBasis**2,NBasis**2))
C
C      call tran4_unsym2(NBasis,NBasis,UMOAO,NBasis,UMOAO,
C     $     NBasis,UMOAO,NBasis,UMOAO,TMPMO)
C
C       block
C       integer :: ip,iq,ir,is,irs,ipq 
C       irs=0
C       do is=1,NBasis
C       do ir=1,NBasis
C       irs=irs+1
C       ipq=0
C       do iq=1,NBasis
C       do ip=1,NBasis
C       ipq=ipq+1
CC
C       write(6,*) ip,iq,ir,is,TwoEl(NAddr3(ip,iq,ir,is))
CC       write(6,*) TwoEl(NAddr3(ip,iq,ir,is)),TMPMO(ipq,irs)
CCC       write(6,*) TMPMO(ipq,irs), TMPMO(irs,ipq)
C       enddo
C       enddo
C       enddo
C       enddo
C       end block
C      call tran4_sym(NBasis,NBasis,transpose(UMOAO),
C     $ NBasis,transpose(UMOAO),
C     $ NBasis,transpose(UMOAO),NBasis,transpose(UMOAO),
C     $ 'TWOMOAA')
CC
C      block
C      integer :: ip,iq,ir,is,irs,ipq
C      integer :: iunit
C      double precision :: work1(NBasis*NBasis)
C      double precision :: work2(NBasis*NBasis)
CC
C      open(newunit=iunit,file='TWOMOAA',status='OLD',
C     $ access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
CC
CC      CHECK-1
CC      irs=0
CC      do is=1,NBasis
CC      do ir=1,NBasis
CC      !irs=irs+1
CC      irs = min(ir,is) + max(ir,is)*(max(ir,is)-1)/2
CC      read(iunit,rec=irs) work1(1:NBasis*(NBasis+1)/2)
CC      call triang_to_sq(work1,work2,NBasis)
CC      ipq=0
CC      do iq=1,NBasis
CC      do ip=1,NBasis
CC      ipq=ipq+1
CC      write(6,*) TwoEl(NAddr3(ip,iq,ir,is)), work2(ipq)
CC      enddo
CC      enddo
CC      enddo
CC      enddo
CC     CHECK-2
C      irs=0
C      do is=1,NBasis
C      do ir=1,is
C      irs=irs+1
C      read(iunit,rec=irs) work1(1:NBasis*(NBasis+1)/2)
C      ipq=0
C      do iq=1,NBasis
C      do ip=1,iq
C      ipq = ipq+1
C      write(6,*) TwoEl(NAddr3(ip,iq,ir,is)), work1(ipq)   
C!      write(6,*) ip,iq,ir,is,TwoEl(NAddr3(ip,iq,ir,is))
C      enddo
C      enddo
C      enddo
C      enddo
C
C      close(iunit)
C      end block
C
C      deallocate(TMPMO)
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
      Subroutine LdInteg(Title,XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UAOMO,NInte1,NBasis,NInte2,NGem)
C     $ TwoEl,UAOMO,NInte1,NBasis,NInte2,NGem,NoSt)
C
C     READ/WRITE THE ONE- AND TWO-ELECTRON INTEGRALS 
C     INTERFACED WITH MOLPRO (INTEGRALS ARE READ FROM FCIDUMP FILES)
C
      use sorter
      use tran
      use abmat 
C
      Implicit Real*8 (A-H,O-Z)
      Parameter (Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Real*8 XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $ UAOMO(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),
     $ UAux(NBasis,NBasis),
     $ Tmp(NInte1) 
C      
      Real*8, Allocatable :: RDM2Act(:)
      Real*8, Allocatable :: HlpRDM2(:)
      Dimension Gamma(NInte1),Work(NBasis),PC(NBasis),
     $ AUXM(NBasis,NBasis),AUXM1(NBasis*NBasis),Ind2(NBasis),
     $ IndInt(NBasis),NumOSym(15),MultpC(15,15),Fock(NBasis*NBasis),
     $ GammaF(NInte1),FockF(NInte1),GammaAB(NInte1),
     $ work1(NBasis,NBasis) 
c herer!!! delete after tests
c     $ ,UMOAOInv(NBasis,NBasis),TwoElAO(NInte2)
C
      Character*60 FName,Aux1,Title
C
      Include 'commons.inc'
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
      XKin(1:NInte1)=Zero
      TwoEl(1:NInte2)=Zero
C
      K=1
C
      If(IAO.Eq.0) Then      
C
C     IF INTEGRALS ARE IN MO REPRESENTATION THEN
C
      FName(K:K+11)='intcoul.dat'
      Call GetENuc(ENuc,FName,NBasis)
C
      Open(10,File='indices_int.dat')
      Read(10,*) NN
C
      If(NN.Eq.1) Then      
C
      FName(K:K+11)='intcoul.dat'
      Call Int1(XKin,XNuc,NInte1,FName,Nbasis)
C
      FName(K:K+11)='intcoul.dat'
      Call Int2(TwoEl,FName,NInte2,NBasis) 
C
      Else
C
      Do I=1,NBasis
      Read(10,*) I1,I2
      IndInt(I1)=I2
      EndDo
C
      Open(20,File=FName)
    2 Read(20,'(A10)')Aux1
      If(Aux1(1:6).Eq."  ISYM".Or.Aux1(1:5).Eq." ISYM") Then
      Read(20,'(A10)')Aux1
      GoTo 20
      EndIf
      GoTo 2
C
   20 Continue
      Read(20,*,End=30) X,I1,I2,I3,I4
c
      If((I1+I2.Ne.0).And.(I3+I4.Eq.0)) Then
      IA=IndInt(I1)
      IB=IndInt(I2)
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)        
      XKin(IAB)=X
      EndIf
      If(I3+I4.Ne.0) Then
      IA1=IndInt(I1)        
      IA2=IndInt(I2)
      IA3=IndInt(I3)
      IA4=IndInt(I4)
      TwoEl(NAddr3(IA1,IA2,IA3,IA4))=X
      EndIf        
C
      GoTo 20
C
   30 Close(20)
C      
      EndIf 
      Close(10)     
C
C     If(IAO.Eq.0)      
      EndIf      
C
      If(IAO.Eq.1) Then
C
C     HAP 
C      Call GetENuc_AO(ENuc,Title)
      Call GetEnuc_AOBin(ENuc,'AOONEINT.mol') 
C
C     READ INFORMATION ABOUT SYMMETRY
C              
C      Open(10,File='indices_int.dat')
C      Read(10,*) MxSym
C      Do I=1,MxSym
C      Read(10,*) X
C      NumOSym(I)=X
C      EndDo
CC
C      If(MxSym.Gt.1) Then
C      Do I=1,NBasis
C      Read(10,*) I1,I2
C      IndInt(I1)=I2
C      EndDo     
C      Else
C      Do I=1,NBasis
C      IndInt(I)=I
C      EndDo      
C      EndIf
CC      
C      Close(10)
CC
CC     
C      Print*, 'TEST0!'
C      Do I=1,NBasis
C      Print*, I,IndInt(I)
C      EndDo
 
C     HAP
      Call create_ind('2RDM',NumOSym,IndInt,MxSym,NBasis)

C     LOAD ONE-ELE INTEGS IN AO
      FName(K:K+8)='xone.dat'
C      Call Int1_AO(XKin,NInte1,FName,NumOSym,Nbasis)
C
C     HAP
      Call readoneint_molpro(XKin,'AOONEINT.mol','ONEHAMIL',
     $     .true.,NInte1)
C      Print*, ' '
C      Print*, 'NInte1',NInte1
C      Do I=1,NInte1
C      write(*,*) XKin(I),Tmp(I)
C      EndDo

C     LOAD TWO-ELE INTEGS IN AO
C      Do I=1,60
C      FName(I:I)=' '
C      EndDo
C      K=0
C    5 K=K+1
C      If (Title(K:K).Ne.' ') Then
C      FName(K:K)=Title(K:K)
C      GoTo 5
C      EndIf
C      FName(K:K+10)='.reg.integ'
CC     
C      If(MxSym.Eq.1) Then
C      MultpC(1,1)=1
C      Else      
C      Open(10,File="multip_table.txt")
C      Do I=1,MxSym
C      Read(10,*)(MultpC(I,J),J=1,I)
C      Do J=1,I
C      MultpC(J,I)=MultpC(I,J)
C      EndDo
C      EndDo
C      Close(10)
C      EndIf
CC
C      Call Int2_AO(TwoEl,NumOSym,MultpC,FName,NInte1,NInte2,NBasis)
C
C     HAP
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3) Then
      Call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT')
      If(ITwoEl.Eq.1) Call LoadSaptTwoEl(3,TwoEl,NBasis,NInte2)
      Else
      Call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT')
      If(ITwoEl.Eq.1) Call LoadSaptTwoEl(4,TwoEl,NBasis,NInte2)
      EndIf
C
C     LOAD AO TO CAS_MO ORBITAL TRANSFORMATION MATRIX FROM uaomo.dat
C      
C      Call GetUAOMO(UAOMO,NumOSym,NBasis)
C     HAP 
      Call read_mo_molpro(UAOMO,'MOLPRO.MOPUN','CASORB  ',NBasis)
C      
C     If(IAO.Eq.1)      
      EndIf              
C
      URe(1:NBasis,1:NBasis)=Zero
      Occ(1:NBasis)=Zero
      PC(1:NBasis)=Zero
      Gamma(1:NInte1)=Zero
C
C     READ RDMs: OLD
      Write(6,'(/," Reading in 1-RDM ...")')
C      Call read_1rdm_molpro(Gamma,InSt(1,1),InSt(2,1),
C     $ '2RDM',IWarn,NBasis)
C     READ RDMs: NEW
      Wght=One/Float(NStates)
      Do I=1,NStates
      GammaAB(1:NInte1)=Zero
      Call read_1rdm_molpro(GammaAB,InSt(1,I),InSt(2,I),
     $ '2RDM',IWarn,NBasis)
      Do K=1,NInte1
      Gamma(K)=Gamma(K)+Wght*GammaAB(K)
      EndDo
      EndDo
C   
C      print*,'Gamma from file:', norm2(Gamma)
C
C      Open(10,File='rdmdump.dat')
C      Read(10,*) NStates
C      IStart=0
C   25 Read(10,*,End=35) X,I1,I2,I3,I4
C      If(X.Eq.NoSt.And.I1+I2+I3+I4.Eq.-4) IStart=1
C      If(X.Ne.NoSt.And.I1+I2+I3+I4.Eq.-4) IStart=0
C      If((I1+I2.Ne.0).And.(I3+I4.Eq.0).And.IStart.Eq.1) Then
C      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
C      Gamma(Ind)=X
C      EndIf
C      GoTo 25
C  35  Close(10)
C
      Call CpySym(AUXM,Gamma,NBasis)
C
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
C     KP
      Call read_nact_molpro(nact,'2RDM')
      If(NAc.Ne.nact) Then
      Write(6,'(1x,"WARNING! The number of partially occ orbitals
     $ different from nact read from molpro. Some active orbitals
     $ must be unoccupied.",/)')
      NAc=nact
      ISwitch=1
      IWarn=IWarn+1
      EndIf
C     KP
C
      NInAc=XELE-Sum+1.D-2
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
      Write(6,'(2X,"No of inactive and active orbitals: ",2I4)')
     $ NInAcCAS,NAcCAS
C
      Write(6,'(2X,"MCSCF",3X,"Occupancy",4X,"Gem")')
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
     $ URe(I,J)=AUXM(IIAct,JJAct)
      EndDo
      EndDo
C
C      print*, norm2(URe)
C      Call print_sqmat(URe,NBasis)
C
      GammaF(1:NInte1)=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      If(I.Le.NInAc.And.J.Le.NInAc.And.I.Eq.J) GammaF(IJ)=One
      If((I.Gt.NInAc.And.I.Le.NInAc+NAc).And.
     $ (J.Gt.NInAc.And.J.Le.NInAc+NAc))
     $ GammaF(IJ)=Gamma(IndSym(I-NInAc,J-NInAc))
      EndDo
      EndDo

C
C      ISkip=1
C      If(ISkip.Ne.1) Then
C     FIND CANONICAL INACTIVE AND VIRTUAL ORBITALS 
C      
      If(IAO.Eq.0) Then
C
      Call FockGen(FockF,GammaF,XKin,TwoEl,NInte1,NBasis,NInte2)
C
      Else
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAux(IndInt(I),J)=UAOMO(J,I)
      EndDo
      EndDo
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      GammaAB(IAB)=Zero
      Do I=1,NBasis
      Do J=1,NBasis
      GammaAB(IAB)=GammaAB(IAB)
     $ +UAux(I,IA)*UAux(J,IB)*GammaF(IndSym(I,J))
      EndDo
      EndDo
      EndDo
      EndDo
C
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3) Then
      Call FockGen_mithap(FockF,GammaAB,XKin,NInte1,NBasis,'AOTWOSORT')
      Else
      Call FockGen_mithap(FockF,GammaAB,XKin,NInte1,NBasis,'AOERFSORT')
      EndIf
C
C     TESTY:
C      print*,'AUXM:', norm2(AUXM)
C      print*,'Gamma:', norm2(Gamma)
C      print*,'GammaF:', norm2(GammaF)
C      print*,'Fock:', norm2(FockF)
C
C      Call FockGen(FockF,GammaAB,XKin,TwoEl,NInte1,NBasis,NInte2)
      Call MatTr(FockF,UAux,NBasis)
C      
      EndIf
C
C     INACTIVE
      If(NInAc.Ne.Zero) Then
C      
      Do I=1,NInAc
      Do J=1,NInAc
      IJ=IndSym(I,J)
      Fock((J-1)*NInAc+I)=FockF(IJ)
      EndDo
      EndDo
      Call Diag8(Fock,NInAc,NInAc,PC,Work)
      Do I=1,NInAc
      Do J=1,NInAc
      URe(I,J)=Fock((J-1)*NInAc+I)
      EndDo
      EndDo
C
      EndIf
C
C     VIRTUAL
C
      NVirt=NBasis-NInAc-NAc
      If(NVirt.Ne.Zero) Then

      Do I=1,NVirt
      Do J=1,NVirt
      IJ=IndSym(I+NInAc+NAc,J+NInAc+NAc)
      Fock((J-1)*NVirt+I)=FockF(IJ)
      EndDo
      EndDo      
      Call Diag8(Fock,NVirt,NVirt,PC,Work)
C      Print*, PC(1:5)
      Do I=1,NVirt
      Do J=1,NVirt
      II=I+NInAc+NAc
      JJ=J+NInAc+NAc
      URe(II,JJ)=Fock((J-1)*NVirt+I)
      EndDo
      EndDo
C
      EndIf
C
C     END OF CANONICALIZING
C     ISkip
C      EndIf
C
C     TRANSFORM INTEGRALS TO NO
C
      If(IAO.Eq.0) Then 
C      
      Call MatTr(XKin,URe,NBasis)
      Write(6,'(/," Transforming two-electron integrals ...",/)')
      Call TwoNO1(TwoEl,URe,NBasis,NInte2)
C
C     If(IAO.Eq.0)      
      Else
C
c herer!!!
c start test
c      Do I=1,NBasis
c      Do J=1,NBasis
cc      UAux(IndInt(I),J)=UAOMO(J,I)
c      UAux(I,J)=UAOMO(J,I)
c      EndDo
c      EndDo
cC
c      Call CpyM(UMOAOInv,UAux,NBasis)
c      tol = 1.0d-7
c      call minvr(UMOAOInv,tol,det,ier,nbasis)
c      if (ier.ne.0) then
c        Write(6,'(/,'' ERROR : transformation from ao to '',
c     $       ''mo basis is singular'')')
c        Stop
c      endif
c      
cc      Call TwoNO1(TwoEl,UAux,NBasis,NInte2)
c      Open(10,File="FCIDUMP")
cC
c   22 Read(10,'(A10)')Aux1
c      If(Aux1(1:6).Eq."  ISYM".Or.Aux1(1:5).Eq." ISYM") Then
c      Read(10,'(A10)')Aux1
c      GoTo 21
c      EndIf
c      GoTo 22
c   21 Continue
c      Read(10,*,End=31) X,I1,I2,I3,I4
c      If(I3+I4.Ne.0) then
c      TwoElAO(NAddr3(I1,I2,I3,I4))=X 
cc      diff=abs(TwoEl(NAddr3(I1,I2,I3,I4))-X)
cc      if(diff.gt.1.d-5) 
cc     $ write(*,*) i1,i2,i3,i4,TwoEl(NAddr3(I1,I2,I3,I4)),X
c      endif
c      GoTo 21
cC
c   31 Close(10)
c    
c      Call TwoNO1(TwoElAO,UMOAOInv,NBasis,NInte2)
cC
c      ij=0
c      do i=1,nbasis
c      do j=1,i
c      ij=ij+1
c      kl=0
c      do k=1,nbasis
c      do l=1,k
c      kl=kl+1
c      if(ij.ge.kl) then
c      diff=abs(TwoEl(NAddr3(i,j,k,l))-TwoElAO(NAddr3(i,j,k,l)))
c      if(diff.gt.1.d-5)
c     $ write(*,*) i,j,k,l,TwoEl(NAddr3(I,j,k,l)),
c     $ TwoElAO(NAddr3(i,j,k,l))
c      endif
c      enddo
c      enddo
c      enddo
c      enddo 
c      stop
cc end test
C      
C
      Call MultpM(UAOMO,URe,UAux,NBasis)
      Call MatTr(XKin,UAOMO,NBasis)
C  
C      Call print_sqmat(UAOMO,NBasis)
C
C
C     ITwoEl
      If(ITwoEl.Eq.1) Then
      Write(6,'(/," Transforming two-electron integrals ...",/)')
      Call TwoNO1(TwoEl,UAOMO,NBasis,NInte2)
C     
      ElseIf(ITwoEl.eq.3) Then
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
      If(ISwitch.Eq.1) Num0=NInAC
      If(ISwitch.Eq.1) Num1=NAc
C     TRANSFORM J AND K
      UAux=transpose(UAOMO)
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3) Then
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
C     TEST MITHAP
C      call tran4_full(NBasis,UAux,UAux,'TWOMO','AOTWOSORT')
C
      Else
      Call tran4_gen(NBasis,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        NBasis,UAux,
     $        'FFOOERF','AOERFSORT')
      Call tran4_gen(NBasis,
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        'FOFOERF','AOERFSORT')
C     TEST MITHAP
C      call tran4_full(NBasis,UAux,UAux,'MO2ERF','AOERFSORT')

      EndIF
C
      EndIf
C              
      EndIf  
C
C     READ ACTIVE 2-RDM AND TRANSFORM TO NO'S  
C
      Write(6,'(" Reading in 2-RDM ...")')
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act=NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act),HlpRDM2(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
      HlpRDM2(1:NRDM2Act)=Zero
C
C     READ RDMs: OLD
C      Call read_2rdm_molpro(RDM2Act,InSt(1,1),InSt(2,1),
C     $ '2RDM',IWarn,NAct)
C
C     READ RDMs: NEW
      Wght=One/Float(NStates)
      Do I=1,NStates
      GammaAB(1:NInte1)=Zero
      Call read_2rdm_molpro(HlpRDM2,InSt(1,I),InSt(2,I),
     $ '2RDM',IWarn,NAct)
      Do K=1,NRDM2Act
      RDM2Act(K)=RDM2Act(K)+Wght*HlpRDM2(K)
      EndDo
      EndDo

C
      Deallocate(HlpRDM2)
C
C      Open(10,File='rdmdump.dat')
C      Read(10,*) NStates
C      IStart=0
C   55 Read(10,*,End=65) X,I1,I2,I3,I4
C      If(X.Eq.NoSt.And.I1+I2+I3+I4.Eq.-4) IStart=1
C      If(X.Ne.NoSt.And.I1+I2+I3+I4.Eq.-4) IStart=0
CC
C      If(I3+I4.Gt.0.And.IStart.Eq.1) Then
C      I=I1
C      J=I2
C      K=I3
C      L=I4        
C      RDM2Act(NAddrRDM(J,L,I,K,NAct))=X*Half
C      RDM2Act(NAddrRDM(L,J,K,I,NAct))=X*Half
C      EndIf
CC
C      GoTo 55
C   65 Close(10)
C
C
      Do I=1,NAct
      Do J=1,NAct
      AUXM1((J-1)*NAct+I)=AUXM(I,J)
      EndDo
      EndDo
      Call TrRDM2(RDM2Act,AUXM1,NAct,NRDM2Act)
C
C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero
      Do I=1,NOccup
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*XKin(II)
      EndDo
      Write(6,'(/,1X,''One-electron Energy'',5X,F15.8)')ETot
C
      EOne=ETot
C
C     ITwoEl
      If(ITwoEl.Eq.1) Then 
      ETwo=Zero 
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      Hlp=FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoEl(NAddr3(IP,IR,IQ,IS))
C      ETot=ETot+Hlp
      ETwo=ETwo+Hlp
      EndDo
      EndDo
      EndDo
      EndDo
C     
      ElseIf(ITwoEl.eq.3) Then
      Call TwoEneChck(ETwo,RDM2Act,Occ,INActive,NAct,NBasis)
      EndIf
C
      Write(6,'(/,1X,''Two-electron Energy'',5X,F15.8)')ETwo
      Write(6,'(/,1X,''MCSCF Molpro Energy'',5X,F15.8)')EOne+ETwo+ENuc
C
C     SAVE THE ACTIVE PART IN rdm2.dat
C
      Open(10,File='rdm2.dat')
      Do I=1,NAct
      Do J=1,NAct
      IJ=(I-1)*NAct+J
      Do K=1,NAct
      Do L=1,NAct
      KL=(K-1)*NAct+L
      If(IJ.Ge.KL) Write(10,'(4I4,F19.12)')
     $ K,I,L,J,Two*RDM2Act(NAddrRDM(I,J,K,L,NAct))
      EndDo
      EndDo
      EndDo
      EndDo
C      
      Close(10)
      Deallocate(RDM2Act)
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
      write(*,*) 'RW-vals:', XELE,NELE

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

*Deck GetENuc_AOBin
      Subroutine GetENuc_AOBin(ENuc,infile)
C
      Implicit Real*8 (A-H,O-Z)
      Character(*) infile
      Logical ex

      inquire(file=trim(infile),EXIST=ex)
     
      if(ex) then
         open(newunit=iunit,file=trim(infile),status='OLD', 
     $        access='SEQUENTIAL',form='UNFORMATTED')
      
         read(iunit)
         read(iunit)
         read(iunit) ENuc
 
         close(iunit)

      else
 
        write(LOUT,'(1x,a)') 'WARNING: '// infile //' NOT FOUND!'
        write(LOUT,'(1x,a)') 'CANNOT READ ENuc!'
        stop

      endif     

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
      Subroutine Int1_AO(XOne,NInte1,FName,NumOSym,NBasis)
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
      Dimension XOne(NInte1),NumOSym(15)
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
      Read(10,*) (XOne(IndSym(Ind+I,Ind+K)),K=1,NumOSym(ISym))
      EndDo
C
      Ind=Ind+NumOSym(ISym)
C   
      EndDo
      Close(10)
C
      Return
      End

*Deck Int2_AO
      Subroutine Int2_AO(TwoEl,NumOSym,MultpC,
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
      Dimension TwoEl(NInte2),
     $ NumOSym(15),MultpC(15,15),Record(NInte1),NTB(8)
C
      Write(6,'(" Loading two-electron integrals in AO ...",/)')
C
      TwoEl(1:NInte2)=         0.D0
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
      TwoEl(NAdd)=Record(ICounter)
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
      TwoEl(NAdd)=Record(ICounter)
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
      Subroutine GetUAOMO(UAOMO,NumOSym,NBasis)
C
C     LOAD A TRANSFORMATION MATRIX AO->MO FROM uaomo.dat
C
      Implicit Real*8 (A-H,O-Z)
      Parameter (Zero=0.D0)
C
      Character*100 Aux
C
      Include 'commons.inc'
C
      Dimension UAOMO(NBasis,NBasis),NumOSym(15)
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAOMO(I,J)=Zero
      EndDo
      EndDo
C
      Open(10,File="uaomo.dat")
      Read(10,'(A10)')Aux
      Read(10,'(A10)')Aux
C
      Ind=0
      Do ISym=1,MxSym
C
      Do I=1,NumOSym(ISym)
      Read(10,*) (UAOMO(Ind+I,Ind+K),K=1,NumOSym(ISym))
      EndDo
C
      Ind=Ind+NumOSym(ISym)
C
      EndDo
      Close(10)
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

****** ALTERNATIVE READ HERE
      subroutine altread2el(TwoEl,UMOAO,NBasis,NInte2)
C
C     Reads 2-el integrals in AO and ttransform to NO
C     Returns TwoEl in NO
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension TwoEl(NInte2),UMOAO(NBasis,NBasis)
C
      Integer(8),external :: NAddr3
C
      integer :: iunit,iunit2
      integer iunit77,iunit88, iunit99
      integer :: maxrep, naos(8), lbuf, nibuf, nbits
      integer :: nints, INDX
      integer,allocatable :: idx_buf(:)
      double precision,allocatable :: val_buf(:), mat(:)
      integer :: idx_p, idx_q, idx_r, idx_s, pq, rs
      logical :: swap_pqrs
      integer :: i

      ! newunit works with Fortran 2008
      open(newunit=iunit,file='AOTWOINT',status='OLD',
     $     access='SEQUENTIAL',form='UNFORMATTED')
    
      ! read info
      call readlabel2(iunit,'BASINFO ')
      read(iunit) maxrep, naos, lbuf, nibuf, nbits 
    
      write(6,'()')
      write(6,'(1x, a)') 'Dalton two-el. file initialized'
      write(6,*) lbuf,nibuf,nbits
C      write(6,'(1x,a,i3,a,i3,a,i3)') 'Buffer size: ', lbuf, &
C           & ', integers per index packet: ', nibuf,       &
C           & ', bits: ', nbits
    
      allocate(val_buf(lbuf))
      allocate(idx_buf(lbuf*nibuf))
    
      call readlabel2(iunit,'BASTWOEL')
    
      select case(nibuf)
      case(1)
    
        do
           read(iunit) val_buf, idx_buf, nints
           if(nints<0) exit
           do i=1,nints
              INDX = idx_buf(i)
              idx_p = ibits(INDX,0,8)
              idx_q = ibits(INDX,8,8)
              idx_r = ibits(INDX,16,8)
              idx_s = ibits(INDX,24,8)
    
              ! pq: position in Batch
              ! rs: Batch number
              pq = idx_p + idx_q*(idx_q-1)/2
              rs = idx_r + idx_s*(idx_s-1)/2 
              !write(*,*) idx_p,idx_q,idx_r,idx_s 
              swap_pqrs = (pq<rs)
           TwoEl(NAddr3(idx_p,idx_q,idx_r,idx_s))=val_buf(i)
           enddo
        enddo
       
      case(2)
    
        do
           read(iunit) val_buf, idx_buf, nints
           if(nints<0) exit
           do i=1,nints
              INDX = idx_buf(i)
              idx_r = ibits(INDX,0,16)
              idx_s = ibits(INDX,16,16)
              INDX = idx_buf(i+lbuf)
              idx_p = ibits(INDX,0,16)
              idx_q = ibits(INDX,16,16)
    
              pq = idx_p + idx_q*(idx_q-1)/2
              rs = idx_r + idx_s*(idx_s-1)/2 
              swap_pqrs = (pq<rs)
           TwoEl(NAddr3(idx_p,idx_q,idx_r,idx_s))=val_buf(i)
           enddo
        enddo
    
      end select
    
      deallocate(val_buf,idx_buf)
      close(unit=iunit)

      Write(6,'(" Transforming two-electron integrals ...",/)')
      Call TwoNO1(TwoEl,UMOAO,NBasis,NInte2)
C
      end subroutine altread2el


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
      Integer(8),external :: NAddr3
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
C       call unpckdlt(ibuf(i), ibuf(i+lbuf), nibuf, nbits,ip,iq,ir,is)
       call unpckdlt(ibuf((i-1)*nibuf + 1), ibuf(i*nibuf),
     $ nibuf, nbits,ip,iq,ir,is)
C        if (dabs(dbuf(i)).gt.1e-8)
c     *          write(*,'(4I4,F12.8)') ip, iq, ir, is, dbuf(i)
      TwoEl(NAddr3(ip,iq,ir,is))=dbuf(i)
      end do
      
      if (nxx .ge. 0) go to 10

C     TEST read2el vs. readtwoint
C      block
C
C      integer :: iunit,ip,iq,ir,is
C      double precision :: TMP1, TMP2
C      double precision :: mat(NBasis*(NBasis+1)/2)
C      Integer(8),external :: NAddr3
CC
C      print*, '2EL-TEST / NAddr3'
CC
C      open(newunit=iunit,file='AOTWOSORT',access='DIRECT',
C     $ recl=8*NBasis*(NBasis+1)/2)
CC
C      do is=1,NBasis
C      do ir=1,NBasis
C      read(iunit,rec=min(ir,is)+max(ir,is)*(max(ir,is)-1)/2) mat
C      do iq=1,NBasis
C      do ip=1,NBasis
CC  
CC      write(*,*) ip,iq,ir,is,TwoEl(NAddr3(ip,iq,ir,is)) 
C      TMP1 = TwoEl(NAddr3(ip,iq,ir,is)) 
C      TMP2 = mat(min(ip,iq)+max(ip,iq)*(max(ip,iq)-1)/2)
CC
C      if(TMP1.ne.TMP2) write(6,*) TMP1,TMP2
CC    
C      enddo
C      enddo
C      enddo
C      enddo
C      close(iunit)
CCC
CC      print*, 'Ended 2-el comparison'
C      end block
CC
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

*Deck SaptInter
      Subroutine SaptInter(NBasis,Mon,ICAS)
C     
C     FEEDS COMMONS.INC WITH SAPT VALUES
      Use types
C     
      Implicit Real*8 (A-H,O-Z)
C      
      type(SystemBlock) :: Mon
C
      Include 'commons.inc'

      Do I=1,NBasis
      CICoef(I) = Mon%CICoef(I)
      IGem(I) = Mon%IGem(I)
      EndDo
     
      If(ICAS.Eq.1) Then
      ICASSCF = ICAS
      NELE = Mon%NELE
      NAcCAS  = Mon%NAct
      NInAcCAS= Mon%INAct
      EndIf

C      write(*,*) 'SINTER,NDimX', Mon%NDimX
C      write(*,*) CICoef(1:NBasis)
C      write(*,*) IGem(1:NBasis)
      
      End 

*Deck LoadSaptTwoEl
      Subroutine LoadSaptTwoEl(Mon,TwoNO,NBasis,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C   
      Integer :: Mon, NInte2, NBasis
      Integer :: IRS, IS, IR, IPQ, IQ, IP
      Integer :: iunit
      Integer(8),external :: NAddr3
      Dimension :: TwoNO(NInte2), Work1(NBasis**2)    
      Character*9 :: fname

      If(Mon.Eq.1) Then
      fname='TWOMOAA '
      ElseIf(Mon.Eq.2) Then
      fname='TWOMOBB '
      ElseIf(Mon.Eq.3) Then
      fname='AOTWOSORT'
      ElseIf(Mon.Eq.4) Then
      fname='AOERFSORT'
      ElseIf(Mon.Eq.5) Then
      fname='MO2ERFAA'
      ElseIf(Mon.Eq.6) Then
      fname='MO2ERFBB'
      EndIf
 
      Work1 = 0d0
      TwoNO = 0d0 
C      write(*,*) trim(fname) 
      open(newunit=iunit,file=trim(fname),status='OLD',
     $ access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

      IRS=0
      Do IS=1,NBasis
      Do IR=1,IS
      IRS=IRS+1
      read(iunit,rec=IRS) Work1(1:NBasis*(NBasis+1)/2)
      IPQ=0
      Do IQ=1,NBasis
      Do IP=1,IQ
      IPQ = IPQ + 1
      TwoNO(NAddr3(IP,IQ,IR,IS)) = Work1(IPQ)   
C      Write(6,*) IP,IQ,IR,IS, Work1(IPQ)
C      Write(6,*) IP,IQ,IR,IS, TwoNO(NAddr3(IP,IQ,IR,IS)) 
      EndDo
      EndDo
      EndDo
      EndDo

C      If(Mon.Eq.3) Then
C      Write(6,*) 'TEST!!'
C      Do IS=1,NBasis
C      Do IR=1,IS
C      Do IQ=1,NBasis
C      Do IP=1,IQ
C
C      Write(6,*) ip,iq,ir,is,TwoNO(NAddr3(IP,IQ,IR,IS))  
C
C      EndDo
C      EndDo
C      EndDo
C      EndDo
C      EndIf

      close(iunit)
      End

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

*Deck FockGen      
      Subroutine FockGen(Fock,Gamma,XOne,TwoEl,NInte1,NBasis,NInte2)
C
C     GENERALIZED FOCK MATRIX
C
      Implicit Real*8 (A-H,O-Z)
      Parameter (Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Real*8 Fock(NInte1),TwoEl(NInte2),Gamma(NInte1),XOne(NInte1)
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      FIJ=XOne(IJ)
      Do K=1,NBasis
      Do L=1,NBasis
      KL=IndSym(K,L)
      FIJ=FIJ+Gamma(KL)*
     $ (Two*TwoEl(NAddr3(I,J,K,L))-TwoEl(NAddr3(I,L,J,K)))
      EndDo
      EndDo
      Fock(IJ)=FIJ
C
      EndDo
      EndDo
C      
      Return
      End

      subroutine readlabel2(iunit,text)
      ! sets file pointer 
      ! to first data after text
      implicit none
      
      integer :: iunit
      integer :: ios
      character(8) :: text, label(4)
      
      rewind(iunit)
      do 
      
        read(iunit,iostat=ios) label
        if(ios<0) then
           write(6,*) 'ERROR!!! Empty section in AOTWOINT!'
           stop
        endif
        if(label(1)=='********') then
           if(label(4)==text) exit
        endif
      
      enddo
      
      end subroutine readlabel2


      subroutine prepare_nums(Occ,Num0,Num1,NBasis)
      Implicit Real*8 (A-H,O-Z)
C       
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension Occ(NBasis)
C
C     LOCAL ARRAYS
C
      Dimension IndAux(NBasis)
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      Do I=1,NBasis
C
      If(Occ(I).Lt.One.And.Occ(I).Ne.Zero) Then
      IndAux(I)=1
      EndIf
      EndDo
C
      Num0=0
      Do I=1,NBasis
      If(IndAux(I).ne.0) Exit
      Num0=Num0+1
      EndDo
      Num2=0
      Do I=NBasis,1,-1
      If(IndAux(I).ne.2) Exit
      Num2=Num2+1
      EndDo
      Num1=NBasis-Num0-Num2
C
      end subroutine prepare_nums

      subroutine TwoEneChck(ETwo,RDM2Act,Occ,INActive,NAct,NBasis)
      Implicit Real*8 (A-H,O-Z)
C       
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Integer INActive,NAct,NBasis
      Double Precision ETwo
      Dimension Occ(NBasis),RDM2Act(NAct**2*(NAct**2+1)/2)
C
C     LOCAL ARRAYS
C
      Dimension Ind(NBasis)
      Double Precision, Allocatable :: RDM2val(:,:,:,:),
     $                                 work(:),ints(:,:)
      Character(:),Allocatable :: IntJFile
C  
C     SET FILES
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3) Then
      IntJFile='FFOO'
      Else
      IntJFile='FFOOERF'
      EndIf
C
C     SET DIMENSIONS
      NOccup=NAct+INActive
      Ind=0
      Do I=1,NAct
      Ind(INActive+I)=I
      EndDo
C
      Allocate(work(NBasis**2),ints(NBasis,NBasis))
      Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
C
      Do L=1,NOccup
      Do K=1,NOccup
      Do J=1,NOccup
      Do I=1,NOccup
      RDM2val(I,J,K,L) = FRDM2(I,K,J,L,RDM2Act,Occ,Ind,NAct,NBasis)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Open(newunit=iunit,file=IntJFile,status='OLD',
     $     access='DIRECT',recl=8*NBasis**2)
C
      ETwo=0
C     COULOMB LOOP (FF|OO)
      kl=0
      Do ll=1,NOccup
      Do kk=1,NOccup
      kl=kl+1
      read(iunit,rec=kl) work(1:NBasis**2)
      Do j=1,NBasis
      Do i=1,NBasis
      ints(i,j) = work((j-1)*NBasis+i)
      EndDo
      EndDo
C 
      k = kk
      l = ll
C
      If(k>NOccup.or.l>NOccup) Cycle
C
      ETwo = ETwo + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
C
      EndDo
      EndDo
C
      Close(iunit)
C
      Deallocate(ints,work)
      Deallocate(RDM2val)
C
      end subroutine TwoEneChck

      subroutine TwoEHartree(EnH,RDM2Act,Occ,INActive,NAct,NBasis)
      Implicit Real*8 (A-H,O-Z)
C       
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Integer INActive,NAct,NBasis
      Double Precision EnH
      Dimension Occ(NBasis),RDM2Act(NAct**2*(NAct**2+1)/2)
C
C     LOCAL ARRAYS
C
      Dimension Ind(NBasis)
      Double Precision, Allocatable :: RDM2val(:,:,:,:),
     $                                 work(:),ints(:,:)
      Character(:),Allocatable :: IntJFile
C
C     SET FILES
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3) Then
      IntJFile='FFOO'
      Else
      IntJFile='FFOOERF'
      EndIf
C
C     SET DIMENSIONS
      NOccup=NAct+INActive
      Ind=0
      Do I=1,NAct
      Ind(INActive+I)=I
      EndDo
C
      Allocate(work(NBasis**2),ints(NBasis,NBasis))
      Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
C
      Do L=1,NOccup
      Do K=1,NOccup
      Do J=1,NOccup
      Do I=1,NOccup
      RDM2val(I,J,K,L) = FRDM2(I,K,J,L,RDM2Act,Occ,Ind,NAct,NBasis)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Open(newunit=iunit,file=IntJFile,status='OLD',
     $     access='DIRECT',recl=8*NBasis**2)

C     GET E_HARTREE     
      EnH=0
C
      kl=0
      Do ll=1,NOccup
      Do kk=1,NOccup
      kl=kl+1
      If(ll==kk) Then
      read(iunit,rec=kl) work(1:NBasis**2)
      Do j=1,NBasis
      Do i=1,NBasis
      ints(i,j) = work((j-1)*NBasis+i)
      EndDo
      EndDo
C 
      k = kk
      l = ll
C
      If(k>NOccup.or.l>NOccup) Cycle
C
      Do i=1,NOccup
      EnH = EnH + Occ(k)*Occ(i)*ints(i,i)
      EndDo
C
      EndIf
C
      EndDo
      EndDo
C
      EnH = 2d0*EnH
C
      Close(iunit)
C
      end subroutine TwoEHartree

