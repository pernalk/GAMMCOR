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
      use sapt_main
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 XKin(NInte1),XNuc(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),UMOAO(NBasis,NBasis),
     $ UAux(NBasis,NBasis)
      integer :: ione,NBas(8),NSymBas(8),NSymOrb(8),nrhf(8),ioprhf
      integer(8) :: MemSrtSize
      logical :: exione,ex
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
      iORCA=Flags%iORCA
      If(IDMRG.Eq.1.And.ICASSCF.Eq.0) Stop 'Set ICASSCF TO 1'
C
      If(IDMRG.Eq.1) Then
      Call ReadDMRG(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem,iORCA)
      Return
      EndIf
C
      If(ICASSCF.Eq.1) Then
C
C     READ UMOAO FROM SIRIUS.RST
C
      inquire(file='SIRIUS.RST',EXIST=ex)
      if(ex) then
      open(newunit=iunit,file='SIRIUS.RST',status='OLD',
     $    access='SEQUENTIAL',form='UNFORMATTED')
      call readlabel(iunit,'BASINFO ')
      read (iunit) NSym,NSymBas,NSymOrb,nrhf,ioprhf
      close(iunit)
      else
      stop 'WARNING: SIRIUS.RST NOT FOUND!'
      endif
C
      Call read_mo_dalton(UAux,NBasis,NSym,NSymBas,NSymOrb,
     $            'SIRIUS.RST','DALTON.MOPUN')
C
      Do I=1,NBasis
      Do J=1,NBasis
      UMOAO(J,I)=UAux(I,J)
      EndDo
      EndDo
C
      Call SortOrbDal(UMOAO,Occ,NInAc,NAc,NSym,NSymOrb,NBasis)
      Do I=1,NInAc+NAc
      Occ(I)=Occ(I)/2.D0
      EndDo
C
      Call read1elsym(XKin,UMOAO,NSym,NSymBas,NBasis,NInte1)
C
      Else
C
C     READ UMOAO FROM DALTON.MOPUN
C
      Open(10,File='DALTON.MOPUN',Form='Formatted',Status='Old')
      Read(10,'(A60)') Line
      Do J=1,NBasis
      Read(10,'(4F18.14)') (UMOAO(J,I),I=1,NBasis)
      EndDo
      Close(10)
      Call read1el(XKin,UMOAO,NBasis,NInte1)
C
      EndIf
C
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=0.D0
      If(I.Eq.J) URe(I,J)=1.0D0
      EndDo
      EndDo
C
C     GET 2-EL NO INTEGRALS AND CICoef
C
      If(ITwoEl.Eq.1) Then
      Call read2el(TwoEl,UMOAO,NBasis,NInte2)
      Else

      MemSrtSize=MemVal*1024_8**MemType
      Call readtwoint(NBasis,1,'AOTWOINT','AOTWOSORT',MemSrtSize)
      EndIf
C
      If(ICASSCF.Eq.0) Then
C
C     read geminal coefficients
      Call ReadCGemDal(CICoef,NELE,INActive,NActive,NBasis)
C
C     set IGem and Occ
      Do I=1,INActive
      IGem(I)=I
      EndDo
      Do I=INActive+1,NELE
      IGem(I)=I
      IGem(NELE+I-INActive)=I
      EndDo
      NGem=NELE+1
      Do I=1,NBasis
      If(CICoef(I).Eq.0.D0) IGem(I)=NGem
      Occ(I)=CICoef(I)**2
      EndDo
C
      ElseIf(ICASSCF.Eq.1) Then
C
      Sum=0.D0
      Do I=1,NInAc+NAc
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

C     OUT-OF-CORE INTEGRAL TRANSFORMATIONS
      If(ITwoEl.Ne.1) Then
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
      If(ISwitch.Eq.1) Num0=NInAC
      If(ISwitch.Eq.1) Num1=NAc
      EndIf
C
      If(ITwoEl.Eq.2) Then
C     FULL TRANSFORMATION FOR "mithap" VARIANT
      UAux=transpose(UMOAO)
      Call tran4_full(NBasis,UAux,UAux,'TWOMO','AOTWOSORT')
C
      ElseIf(ITwoEl.Eq.3) Then
C     TRANSFORM J AND K
      UAux=transpose(UMOAO)
c     print*, 'Num0-1',Num0,Num1
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
      If(ITwoEl.Gt.1) Then
C     DELETE SORTED AOTWOINTS FOFO
      Write(6,'(/A)') 'SORTED 2-el AO ints deleted after tran4_gen!'
      Call delfile('AOTWOSORT')
      EndIf
C
      Inquire(File='AOONEINT',EXIST=exione)
      If(exione) Then
C
      Open(newunit=ione,File='AOONEINT',access='SEQUENTIAL',
     $     Form='UNFORMATTED',Status='OLD')
      Read(ione)
      Read(ione) NSym,NBas(1:NSym),ENuc
      Write(6,'(/,"  Nuclear repulsion:",F20.12)')ENuc
      Close(ione)
C
      Else
C
      Open(10,File='enuc.dat',Form='Formatted',Status='Old')
      Read(10,'(A31,F20.12)')Line,ENuc
      Write(6,'(/,"  Nuclear repulsion:",F20.12)')ENuc
      Close(10)
      EndIf
C
      Return
      End

C
C     ReadDMRG: VERSION THAT WORKS WITH EUGENE'S INTS (IEugene=1)
C               IT WORKS WITH ORCA OUTPUTS (IEugene=0)
C
*Deck ReadDMRG
      Subroutine ReadDMRG(XKin,XNuc,ENuc,Occ,URe,
     $ TwoEl,UMOAO,NInte1,NBasis,NInte2,NGem,iORCA)
C
      use types
      use sorter
c     use Cholesky_old  ! create AOTWOSORT file
      use Cholesky
      use tran
      use abmat
C
C     READ Integrals and 1-RDM - needed for AC-DMRG CALCULATION
C
      Implicit Real*8 (A-H,O-Z)
C
!      Parameter (Zero=0.D0,One=1.D0,Two=2.D0)
C
      Real*8 XKin(NInte1),XNuc(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),UMOAO(NBasis*NBasis)
      Type(TCholeskyVecs) :: CholeskyVecs
C
      Character*60 FName,Aux1
C
      Include 'commons.inc'
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2(:),RDMAB2(:)
      Real*8, Allocatable :: MatFF(:,:)
      Dimension Gamma(NInte1),Work(NBasis),PC(NBasis),
     $ AUXM(NBasis,NBasis),AUXM1(NBasis,NBasis),
     $ Fock(NBasis*NBasis),
     $ UAux(NBasis,NBasis),
     $ FockF(NInte1),GammaAB(NInte1),Eps(NBasis,NBasis)
      Integer(8) MemSrtSize,IOutInfo
      Dimension FockF2(NInte1),WorkSq(NBasis,NBasis)

      If(iORCA==1) then
      LiborNew=1
      IEugene=0
      ElseIf(iORCA==0) then
      IEugene=1
      EndIf
C
      UMOAO(1:NBasis*NBasis)=Zero
      URe(1:NBasis,1:NBasis)=Zero
      Occ(1:NBasis)=Zero
      PC(1:NBasis)=Zero
      Gamma(1:NInte1)=Zero
C
C     READ IN 1-RDM AND DIAGONALIZE IT
C
      If(IEugene.Eq.0) Then
C
      If(LiborNew.Eq.1) Then
C
      Open(10,File='G1.bin',form='unformatted', access='stream',
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
C
      Else
C
      Open(10,File='rdmdump.dat',Status='Old')
   20 Read(10,*,End=30) X,I1,I2,I3,I4
      If((I1+I2.Ne.0).And.(I3+I4.Eq.0)) Then
      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
      Gamma(Ind)=X
      EndIf
      GoTo 20
  30  Close(10)
C
C     If(LiborNew.Eq.1) Then
      EndIf
C
      Else
C
      Open(10,File='rdmdump.dat',form='unformatted',access='stream',
     $ Status='Old')
   27 Read(10,End=37) X,I1,I2,I3,I4
      If((I1+I2.Ne.0).And.(I3+I4.Eq.0)) Then
      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
      Gamma(Ind)=X/Two
      EndIf
      GoTo 27
   37 Close(10)
C
      EndIf
C
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
      TwoEl(1:NInte2)=Zero
C
      If(IEugene.Eq.0) Then
C
C IBin=1 - integrals in binary files
C IBin=0 - integrals in text files
C
      IBin=1
c      IBin=0
C
      If(IBin.Eq.0) Then
C
      Write(6,'(" Reading in one-ele integrals ...",/)')
      FName(1:12)='intcoul.dat'
      Call GetENuc(ENuc,FName,NBasis)
      Call Int1(XKin,XNuc,NInte1,FName,Nbasis)
      Write(6,'(" Reading in two-ele integrals ...",/)')
      Call Int2(TwoEl,FName,NInte2,NBasis)
C
      ElseIf(IBin.Eq.1) Then
C
      Open(10,File='NUC_REP.bin',form='unformatted',access='stream',
     $ Status='Old')
      Read(10) X
      ENuc=X
      Close(10)
C
      Open(10,File='FACT.bin',form='unformatted',access='stream',
     $ Status='Old')
      If(LiborNew.Eq.1) Read(10)I,J,K
      ICount=0
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Read(10,End=31) X
      Ind=I*(I-1)/2+J
      XKin(Ind)=X
      ICount=ICount+1
      EndDo
      EndDo
   31 Close(10)
      Write(6,'(" The number of 1-el integrals read vs. expected",
     $ 2I10)') ICount,NInte1
C
      If(ITwoEl.Eq.1) Then
C
      Open(10,File='DPQRS.bin',form='unformatted',access='stream',
     $ Status='Old')
      If(LiborNew.Eq.1) Read(10)I,J,K
      ICount=0
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
      If(IJ.Ge.KL) Then
      Read(10,End=35) X
      TwoEl(NAddr3(I,J,K,L))=X
      ICount=ICount+1
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
   35 Close(10)
      Write(6,'(" The number of 2-el integrals read vs. expected",
     $ 2I15)') ICount,NInte2
C
      ElseIf(ITwoEl.Gt.1) Then
C
      If(ICholesky==0) Then
      MemSrtSize=MemVal*1024_8**MemType
      Call readtwoint(NBasis,4,'DPQRS.bin','AOTWOSORT',
     $                MemSrtSize,IOutInfo)
C
      Else If(ICholesky==1) Then
c     If(ICholesky==1) Then
c     print*, 'use Cholesky-old'
!     Call chol_CoulombMatrix(CholeskyVecs,'AOTWOSORT',ICholeskyAccu)
c     print*, 'use Cholesky-new'
      Call chol_CoulombMatrix(CholeskyVecs,NBasis,'DPQRS.bin',4,
     &                        ICholeskyAccu)
      NCholesky=CholeskyVecs%NCholesky
      EndIf
C
      EndIf
C
c     If(IBin.Eq.0)
      EndIf
C
C else of IEugene.Eq.0
      Else
C
      If(ITwoEl.Eq.1) Then
C
      Open(10,File='intcoul.dat',form='unformatted',access='stream',
     $ Status='Old')
   40 Read(10,End=39) X,I1,I2,I3,I4
C
      If(I3+I4.Ne.0) Then
      TwoEl(NAddr3(I1,I2,I3,I4))=X
      Else
      If(I1+I2.Ne.0) Then
      Ind=(Max(I1,I2)*(Max(I1,I2)-1))/2+Min(I1,I2)
      XKin(Ind)=X
      EndIf
      If(I1+I2.Eq.0) ENuc=X
      EndIf
C
      GoTo 40
   39 Close(10)
C
      ElseIf(ITwoEl.Gt.1) Then
C
      MemSrtSize=MemVal*1024_8**MemType
      Call readtwoint(NBasis,3,'intcoul.dat','AOTWOSORT',
     $                MemSrtSize,IOutInfo)
C      Call CheckSaptTwoEl(3,TwoEl,NBasis,NInte2)
C      Call LoadSaptTwoEl(3,TwoEl,NBasis,NInte2)
      Call readoneint_eugene(XKin,ENuc,'intcoul.dat',NInte1,IOutInfo)
C
      EndIf
C
      EndIf
C
C     FIND CANONICAL INACTIVE ORBITALS
C
      Call CpySym(AUXM1,Gamma,NBasis)
C
      If(ITWoEl.Eq.3) Then
C
      UAux = 0
      Do I=1,NInAc
      UAux(I,I) = Occ(I)
      EndDo
      Do J=1,NAc
      Do I=1,NAc
      UAux(NInAc+I,NInAc+J) = AUXM1(I,J)
      EndDo
      EndDo
      Call sq_to_triang2(UAux,GammaAB,NBasis)
C
      If(ICholesky==0) Then
         Call FockGen_mithap(FockF,GammaAB,XKin,NInte1,NBasis,
     &                       'AOTWOSORT')
      ElseIf(ICholesky==1) Then
         Call FockGen_CholR(FockF,CholeskyVecs%R(1:NCholesky,1:NInte1),
     &                      GammaAB,XKin,NInte1,NCholesky,NBasis)
      EndIf
C
C      WorkSq = 0
C      call triang_to_sq2(FockF,WorkSq,NBasis)
C      Err = 0
C      do J=1,NInAc
C      do I=1,j
C      if(i.ne.j) Err = Err + WorkSq(i,j)**2
C      enddo
C      enddo
C      do J=NInAc+NAc+1,NBasis
C      do I=NInAc+NAc+1,J
C      if(i.ne.j) Err = Err + WorkSq(i,j)**2
C      enddo
C      enddo
C      Print*, 'Err-2',Sqrt(Err)
C
C     INACTIVE
      If(NInAc.Ne.Zero) Then
C
      Do I=1,NInAc
      Do J=1,NInAc
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
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
C      Print*, 'INAct-MY', norm2(URe)
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
C
C      Print*, PC(1:5)
C      Print*, 'VIRT-MY',norm2(Fock)
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
      ElseIf(ITWoEl.Eq.1) Then
C
C     INACTIVE
C
      If(NInAc.Ne.0) Then
C
      Do I=1,NInAc
      Do J=1,NInAc
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      FIJ=XKin(IJ)
      Do K=1,NInAc
      FIJ=FIJ+Occ(K)*(Two*TwoEl(NAddr3(I,J,K,K))-TwoEl(NAddr3(I,K,J,K)))
      EndDo
      Do K=1+NInAc,NInAc+NAc
      Do L=1+NInAc,NInAc+NAc
      FIJ=FIJ+AUXM1(K-NInAc,L-NInAc)*
     $ (Two*TwoEl(NAddr3(I,J,K,L))-TwoEl(NAddr3(I,L,J,K)))
      EndDo
      EndDo
      Fock((J-1)*NInAc+I)=FIJ
      EndDo
      EndDo
C
      Call Diag8(Fock,NInAc,NInAc,PC,Work)
C
      Do I=1,NInAc
      Do J=1,NInAc
      URe(I,J)=Fock((J-1)*NInAc+I)
      EndDo
      EndDo
C
      Do I=1,NInAc
      Do J=1,NBasis
      Eps(J,I)=PC(I)
      EndDo
      EndDo
C
C      Print*, 'INAct-KA',norm2(URe)
      EndIf
C
C     VIRTUAL
C
      NVirt=NBasis-NInAc-NAc
C
      If(NVirt.Ne.0) Then
C
      Do I=1,NVirt
      Do J=1,NVirt
      II=I+NInAc+NAc
      JJ=J+NInAc+NAc
      IJ=(Max(II,JJ)*(Max(II,JJ)-1))/2+Min(II,JJ)
      FIJ=XKin(IJ)
      Do K=1,NInAc
      FIJ=FIJ+Occ(K)*
     $ (Two*TwoEl(NAddr3(II,JJ,K,K))-TwoEl(NAddr3(II,K,JJ,K)))
      EndDo
      Do K=1+NInAc,NInAc+NAc
      Do L=1+NInAc,NInAc+NAc
      FIJ=FIJ+AUXM1(K-NInAc,L-NInAc)*
     $ (Two*TwoEl(NAddr3(II,JJ,K,L))-TwoEl(NAddr3(II,L,JJ,K)))
      EndDo
      EndDo
      Fock((J-1)*NVirt+I)=FIJ
      EndDo
      EndDo
C
      Call Diag8(Fock,NVirt,NVirt,PC,Work)
C
      Open(10,file='fock.dat')
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
      Eps(IP,IQ)=PC(IP-NOccup)-Eps(IP,IQ)
      Write(10,*)IP,IQ,Eps(IP,IQ)
      EndDo
      EndDo
      Close(10)
C
C      Print*, 'VIRT-Ka',norm2(Fock)
C
      Do I=1,NVirt
      Do J=1,NVirt
      II=I+NInAc+NAc
      JJ=J+NInAc+NAc
      URe(II,JJ)=Fock((J-1)*NVirt+I)
      EndDo
      EndDo
C
      EndIf
C     end of ITwoEl==1
      EndIf
C
C     END OF CANONICALIZING
C
C     CHECK IF URe IS A UNIT MATRIX
C
      IUNIT=1
      Err=Zero
      Do I=1,NBasis
      Do J=1,NBasis
      If(I.Eq.J) Err=Err+Abs(One-Abs(URe(I,J)))
      If(I.Ne.J) Err=Err+Abs(URe(I,J))
      EndDo
      EndDo
C     Print*, 'Err',Err
      If(Err.Gt.1.D-5) IUNIT=0
C     If(Err.Gt.1.D-4) IUNIT=0 ! this should work with Cholesky/Ludicrous
C
      If(IUNIT.Eq.1) Then
      Write(6,'(/,X,"URe is a unit matrix up to ",E16.6)') ERR
      Write(6,'(X,"do not transform integrals")')
      Do I=1,NBasis
      Do J=1,NBasis
      If(I.Eq.J) URe(I,J)=One
      If(I.Ne.J) URe(I,J)=Zero
      EndDo
      EndDo
      EndIf
C
      If(IUNIT.Eq.0) Then
C
      Write(6,'(" Transforming two-electron integrals ...",/)')
C
      Call MatTr(XKin,URe,NBasis)
C
      If(ITwoEl.Eq.1) Then
      Call TwoNO1(TwoEl,URe,NBasis,NInte2)
C
      ElseIf(ITwoEl.Eq.3) Then
C
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
      If(ISwitch.Eq.1) Num0=NInAC
      If(ISwitch.Eq.1) Num1=NAc
C     TRANSFORM J AND K
      UAux=transpose(URe)
      If(ICholesky==0) Then
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
      ElseIf(ICholesky==1) Then
      Allocate(MatFF(NCholesky,NBasis**2))
C       print*, 'use chol_MOTransf-old,IUNIT',IUNIT
C       Call chol_MOTransf(MatFF,CholeskyVecs,
C     $              UAux,1,NBasis,
C     $              UAux,1,NBasis)
C
      MemMOTransf = 1200 ! this is fixed for now in MB
                         ! will be changed in next revision
      Call chol_MOTransf_TwoStep(MatFF,CholeskyVecs,
     $              UAux,1,NBasis,
     $              UAux,1,NBasis,
     $              MemMOTransf)
C
      Call chol_ints_fofo(NBasis,Num0+Num1,MatFF,
     $                    NBasis,Num0+Num1,MatFF,
     $                    NCholesky,NBasis,'FOFO')
      Call chol_ints_fofo(NBasis,NBasis,MatFF,
     $                    Num0+Num1,Num0+Num1,MatFF,
     $                    NCholesky,NBasis,'FFOO')

      open(newunit=iunt,file='cholvecs',form='unformatted')
      write(iunt) NCholesky
      write(iunt) MatFF(1:NCholesky,1:NBasis**2)
      close(iunt)

      Deallocate(MatFF)
      EndIf
CC
      EndIf
C
      ElseIf(IUNIT.Eq.1.And.ITwoEl.Eq.3) Then
C
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
      If(ISwitch.Eq.1) Num0=NInAC
      If(ISwitch.Eq.1) Num1=NAc
C     READ J AND K AND DUMP TO DISC
      UAux = 0d0
      Do I=1,NBasis
      UAux(I,I) = 1d0
      EndDo
!
      If(ICholesky==1) then
c      print*, 'use chol_triang_fofo,IUNIT',IUNIT
         call chol_triang_fofo(NBasis,NBasis,
     $                  CholeskyVecs%R(1:NCholesky,1:NInte1),
     $                  Num0+Num1,Num0+Num1,
     $                  CholeskyVecs%R(1:NCholesky,1:NInte1),
     $                  NCholesky,NInte1,NBasis,'FFOO')
         call chol_triang_fofo(NBasis,Num0+Num1,
     $                  CholeskyVecs%R(1:NCholesky,1:NInte1),
     $                  NBasis,Num0+Num1,
     $                  CholeskyVecs%R(1:NCholesky,1:NInte1),
     $                  NCholesky,NInte1,NBasis,'FOFO')

      Allocate(MatFF(NCholesky,NBasis**2))
      do i=1,NCholesky
         call triang_to_sq(CholeskyVecs%R(i,1:NInte1),MatFF(i,:),NBasis)
      enddo
      open(newunit=iunt,file='cholvecs',form='unformatted')
      write(iunt) NCholesky
      write(iunt) MatFF(1:NCholesky,1:NBasis**2)
      close(iunt)
      Deallocate(MatFF)

      Else
         Call read4_gen(NBasis,
     $           Num0+Num1,Num0+Num1,NBasis,NBasis,
     $           'FFOO','AOTWOSORT')
         Call read4_gen(NBasis,
     $           NBasis,Num0+Num1,NBasis,Num0+Num1,
     $           'FOFO','AOTWOSORT')
C
      EndIf ! ICholesky
C
      EndIf
C
C     CHECK IF INACT AND VIRT ORBITALS ARE CANONICAL
C
      If(IBin.Ge.0) Then
      NBSave=NBasis
      NBasis=NAc
      NInAc=0
      EndIf
C
c      Write(6,'(" Skipping reading 2-RDM from rdmdump.dat.",/)')
c      Write(6,'(" 2-RDM will be read from rdm2.dat file ...",/)')
c      GoTo 888
C
C     READ ACTIVE 2-RDM AND TRANSFORM TO NO'S
C
      Write(6,'(/," Reading 2-RDM ...")')
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      Allocate (RDM2(NRDM2))
      RDM2(1:NRDM2)=Zero
C
      If(IEugene.Eq.0) Then
C
      If(LiborNew.Eq.1) Then
C
      Open(10,File='G2.bin',form='unformatted',access='stream',
     $ Status='Old')
      Read(10)I,J,K
C
      Do K1=1,NAc
      Do J1=1,NAc
      Do L1=1,NAc
      Do I1=1,NAc
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
C     If(LiborNew.Eq.1) Then
      Else
C
      Allocate (RDMAB2(NRDM2))
      RDMAB2(1:NRDM2)=Zero
C
      IAAAA=1
      IABBA=0
      Open(10,File='rdmdump.dat',Status='Old')
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
C     If(LiborNew.Eq.1) Then
      EndIf
C
c else of IEugene.Eq.0
      Else
C
      Open(10,File='rdmdump.dat',form='unformatted',access='stream',
     $ Status='Old')
   28 Read(10,End=38) X,I1,I2,I3,I4
C
      If((I3+I4.Ne.0)) Then
      X=X/Two
      I=I1+NInAc
      J=I2+NInAc
      K=I3+NInAc
      L=I4+NInAc
      RDM2(NAddrRDM(I,J,K,L,NBasis))=X
      RDM2(NAddrRDM(J,I,L,K,NBasis))=X
      EndIf
C
      GoTo 28
   38 Close(10)
C
      EndIf
C
      If(IBin.Ge.0) Then
C
      Do I=1,NBasis
      Do J=1,NBasis
      UMOAO((J-1)*NBasis+I)=AUXM(I,J)
      EndDo
      EndDo
C
      If(IUNIT.Eq.0) Call TrRDM2(RDM2,UMOAO,NBasis,NRDM2)
      GoTo 777
C
      EndIf
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
      ETot=ETot+Two*Occ(I)*XKin(II)
      EndDo
C
      Do IP=1,INActive
      Do IQ=1,INActive
      Do IR=1,INActive
      Do IS=1,INActive
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      ETot=ETot+RDM2(IAdd)*TwoEl(NAddr3(IP,IR,IQ,IS))
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
     $ 2.D0*TwoEl(NAddr3(I,I,J,J))-TwoEl(NAddr3(I,J,I,J))
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
      ETot=ETot+RDM2(IAdd)*TwoEl(NAddr3(IP,IR,IQ,IS))
      if(occ(ip).ne.1.d0.and.occ(ir).ne.1.d0.and.occ(iq).ne.1.d0.and.
     $ occ(is).ne.1.d0) then
      etot2=etot2+RDM2(IAdd)*TwoEl(NAddr3(IP,IR,IQ,IS))
      endif
      EndDo
      EndDo
      EndDo
      EndDo
      Write(*,*)'Active Two-Electron Energy',etot2
      write(*,*)'Total DMRG Energy',etot+enuc
c *****************************************

C
  777 Continue
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
c herer!!!
c      If(IJ.Ge.KL) Write(*,'(8I4,F19.12)')
c     $ KKAct,IIAct,LLAct,JJAct,I,J,K,L,
c     $Two*RDM2(NAddrRDM(I,J,K,L,NBasis))
      EndDo
      EndDo
      EndDo
      EndDo
      Close(10)
      Deallocate(RDM2)
      If(IEugene.Eq.0.And.LiborNew.Eq.0) Deallocate(RDMAB2)
      If(IBin.Ge.0) NBasis=NBSave

c herer!!!
  888 Continue
C
C     dump URe, it may be usefull
C
      open(10,file='ure_casno.dat')
      write(10,*)URe
      close(10)
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
      use types
      use sorter
      use tran
c     use Cholesky_old  ! requires AOTWOSORT
      use Cholesky
      use abmat
      use timing
C
      Implicit Real*8 (A-H,O-Z)
      Parameter (Half=0.5D0)
C      Parameter (Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
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
      Integer(8) :: MemSrtSize
      Type(TCholeskyVecs) :: CholeskyVecs
      Real*8, Allocatable :: MatFF(:,:)
      Real*8 Tcpu,Twall
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
C     HAP
      Call create_ind('2RDM',NumOSym,IndInt,MxSym,NBasis)
C
C     LOAD ONE-ELE INTEGS IN AO
      FName(K:K+8)='xone.dat'
C      Call Int1_AO(XKin,NInte1,FName,NumOSym,Nbasis)
C
      Call readoneint_molpro(XKin,'AOONEINT.mol','ONEHAMIL',
     $     .true.,NInte1)
C
C     SET TIMING FOR 2-el integrals
      Call clock('START',Tcpu,Twall)
C
      If(ICholesky==0) Then
C     memory allocation for sorter
      MemSrtSize=MemVal*1024_8**MemType
C     KP: If IFunSR=6 integrals are not needed and are not loaded
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5) Then
      Call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT',MemSrtSize)
      If(ITwoEl.Eq.1) Call LoadSaptTwoEl(3,TwoEl,NBasis,NInte2)
      ElseIf(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
      Call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
      If(ITwoEl.Eq.1) Call LoadSaptTwoEl(4,TwoEl,NBasis,NInte2)
      EndIf
C
      Else If(ICholesky==1) Then
c     If(ICholesky==1) Then
c     Call chol_CoulombMatrix(CholeskyVecs,'AOTWOSORT',ICholeskyAccu)
       Call chol_CoulombMatrix(CholeskyVecs,NBasis,'AOTWOINT.mol',2,
     &                         ICholeskyAccu)
      NCholesky=CholeskyVecs%NCholesky
      EndIf ! ICholesky
      Call clock('2-electron ints',Tcpu,Twall)
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
C
C      Call read_1rdm_molpro(Gamma,InSt(1,1),InSt(2,1),ISpinMs2,
C     $ '2RDM',IWarn,NBasis)
C
C     READ RDMs: NEW
      Wght=One/Float(NStates)
      Do I=1,NStates
      GammaAB(1:NInte1)=Zero
      Call read_1rdm_molpro(GammaAB,InSt(1,I),InSt(2,I),ISpinMs2,
     $ '2RDM',IWarn,NBasis)
      Do K=1,NInte1
      Gamma(K)=Gamma(K)+Wght*GammaAB(K)
      EndDo
      EndDo
C
      Call CpySym(AUXM,Gamma,NBasis)
C
C KP 30.07.2020
      Call read_nact_molpro(nact,'2RDM')
      NAc=nact
C
C     DIAGONALIZE ONLY THE ACTIVE BLOCK OF Gamma TO AVOID THROWING AWAY
C     ACTIVE ORBITAL OF ZERO-OCCUPANCY (which may happen for atoms for some states)
C
      Do I=1,NAc
      Do J=1,NAc
      AUXM1((J-1)*NAc+I)=AUXM(I,J)
      EndDo
      EndDo
      Call Diag8(AUXM1,NAc,NAc,PC,Work)
      Call SortP(PC,AUXM1,NAc)
C
C KP 30.07.2020
      AUXM(1:NBasis,1:NBasis)=Zero
      Do I=1,NAc
      Do J=1,NAc
      AUXM(I,J)=AUXM1((J-1)*NAc+I)
      EndDo
      EndDo
C
      Sum=Zero
      NAc=0
      Do I=1,NBasis
C KP 08.08.2020
C     it may happen that an active orbital has a negative but very small occupation. set it to a positive
      PC(I)=Abs(PC(I))
      Sum=Sum+PC(I)
      If(PC(I).Gt.Zero) NAc=NAc+1
      EndDo
C
C KP 30.07.2020, no need to call read_nact_molpro again
c      Call read_nact_molpro(nact,'2RDM')
      If(NAc.Ne.nact) Then
      Write(6,'(1x,"WARNING! The number of partially occ orbitals
     $ different from nact read from molpro. Some active orbitals
     $ must be unoccupied.",/)')
      NAc=nact
      ISwitch=1
      IWarn=IWarn+1
      EndIf
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
c herer!!!
C     skip canonicalization if CASPIDFTOPT is called
      If(IFunSR.Eq.6) GoTo 543
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
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5.Or.IFunSR.Eq.6) Then
          If(ICholesky==0) Then
          Call FockGen_mithap(FockF,GammaAB,XKin,NInte1,NBasis,
     &                        'AOTWOSORT')
          ElseIf(ICholesky==1) Then
          Call FockGen_CholR(FockF,CholeskyVecs%R(1:NCholesky,1:NInte1),
     &                       GammaAB,XKin,NInte1,NCholesky,NBasis)
          EndIf
      ElseIf (IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
          If(ICholesky==0) Then
             Call FockGen_mithap(FockF,GammaAB,XKin,NInte1,NBasis,
     &                           'AOERFSORT')
          ElseIf(ICholesky==1) Then
             Write(6,'(1x,a)') 'Cholesky not ready for LR-ERF!'
             Stop
          EndIf
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
C KP 15.05.2019
      Do I=1,NInAc
      Do J=1,NBasis
      work1(J,I)=PC(I)
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
c KP 15.05.2019
      open(10,file='fock.dat')
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
      work1(IP,IQ)=PC(IP-NOccup)-work1(IP,IQ)
      write(10,*)ip,iq,work1(ip,iq)
      EndDo
      EndDo
      close(10)
C
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
C
  543 Continue
C
C     If CASPiDFT then skip integral transformation
      If(IFunSR.Eq.6) Then
C
      Do I=1,NBasis
      Do J=1,NBasis
      UAux(IndInt(I),J)=UAOMO(J,I)
      EndDo
      EndDo
      Call MultpM(UAOMO,URe,UAux,NBasis)
      Write(6,'(/," CASPIDFT - skip integral transformation",/)')
C
      Else
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
      Call MultpM(UAOMO,URe,UAux,NBasis)
      Call MatTr(XKin,UAOMO,NBasis)
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
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5) Then
      If (ICholesky==0) Then
C
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
      Call clock('tran4_FOFO',Tcpu,Twall)
C
      ElseIf (ICholesky==1) Then
C
      Allocate(MatFF(NCholesky,NBasis**2))
C     Old 1-step transformation (much slower)
c     Call chol_MOTransf(MatFF,CholeskyVecs,
c    $                   UAux,1,NBasis,
c    $                   UAux,1,NBasis)
C
      Call chol_MOTransf_TwoStep(MatFF,CholeskyVecs,
     $              UAux,1,NBasis,
     $              UAux,1,NBasis,
     $              1500)
C
      Call clock('chol_NOTransf',Tcpu,Twall)
C
      Call chol_ints_fofo(NBasis,Num0+Num1,MatFF,
     $                    NBasis,Num0+Num1,MatFF,
     $                    NCholesky,NBasis,'FOFO')
      Call chol_ints_fofo(NBasis,NBasis,MatFF,
     $                    Num0+Num1,Num0+Num1,MatFF,
     $                    NCholesky,NBasis,'FFOO')
C
      Call clock('chol_FFOOFOFO',Tcpu,Twall)
C
C KP 07.2021: dump MatFF
C
      open(newunit=iunit,file='cholvecs',form='unformatted')
      write(iunit) NCholesky
      write(iunit) MatFF
      close(iunit)
      Deallocate(MatFF)
C
      EndIf
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
c     If(IFunSR.Eq.5) Then
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
C     $ ISpinMs2,'2RDM',IWarn,NAct)
C
C     READ RDMs: NEW
      Wght=One/Float(NStates)
      Do I=1,NStates
      GammaAB(1:NInte1)=Zero
      Call read_2rdm_molpro(HlpRDM2,InSt(1,I),InSt(2,I),ISpinMs2,
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
      If(IFunSR.Ne.6) Then
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
      EndIf
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
C     Reads 1-el integrals in AO and transform to NO
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

      close(iunit77)
      return
      end

      subroutine read1elsym(XKin,UMOAO,NSym,NSymBas,NBasis,NInte1)
C
C     Reads 1-el integrals in AO in symm blocks and transform to NO
      Implicit Real*8 (A-H,O-Z)
C
      Dimension XKin(NInte1),UMOAO(NBasis,NBasis),NSymBas(8)
C
      parameter (mxbuf = 10000)  ! KP
      double precision  dbuf(mxbuf)
      integer ibuf(mxbuf*2),iibuf(mxbuf*2)
      integer iunit77,iunit88, iunit99, ndim, norb, nbas, nfone
      logical iprtvc
      integer  maxrep, naos, lbuf, nibuf, nbits
      common /daltwoel/  maxrep, naos(8), lbuf, nibuf, nbits
      integer nxx
C
      XKin(1:NInte1)=         0.D0
C
      ip=0
      ipq=0
      ipqs=0
      do i=1,NSym
      do ii=1,NSymBas(i)
      ip=ip+1
c
      iq=0
      do j=1,NSym
      do jj=1,NSymBas(j)
      iq=iq+1
c
      if(ip.ge.iq) then
c
      ipq=ipq+1
      if(i.eq.j) then
      ipqs=ipqs+1
      iibuf(ipqs)=ipq
      endif
c
      endif
c
      enddo
      enddo
      enddo
      enddo
C
      iunit77=77
      iunit88=88
      iunit99=99

      OPEN(UNIT=iunit77,FILE='AOONEINT',STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      call initoneel(iunit77)

10    continue
      call readoneel(iunit77, dbuf, ibuf, lbuf, nibuf, nxx)

      do i=1,nxx
      XKin(iibuf(ibuf(i))) = dbuf(i)
      end do

      if (nxx .ge. 0) go to 10
C
C     TRANSFORM THE INTEGRALS
C
      Write(6,'(" Transforming One-electron integrals ...",/)')
        write(55,*) XKin
      Call MatTr(XKin,UMOAO,NBasis)
        write(56,*) XKin

      close(iunit77)

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
      use sorter
C
      Implicit Real*8 (A-H,O-Z)
C
      Type(AOReaderData) :: reader
C
      Integer :: Mon, NInte2, NBasis
      Integer :: IRS, IS, IR, IPQ, IQ, IP
      Integer :: iunit
      Integer(8),external :: NAddr3
      Dimension :: TwoNO(NInte2), Work1(NBasis**2)
      Logical :: empty
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
C      open(newunit=iunit,file=trim(fname),status='OLD',
C     $ access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

      call reader%open(trim(fname))

      IRS=0
      Do IS=1,NBasis
      Do IR=1,IS
      IRS=IRS+1
C
      !read(iunit,rec=IRS) Work1(1:NBasis*(NBasis+1)/2)
      Call reader%getTR(IRS,Work1,empty)
      If(empty) Then
         Work1 = 0d0
      Else
        IPQ=0
        Do IQ=1,NBasis
        Do IP=1,IQ
        IPQ = IPQ + 1
        TwoNO(NAddr3(IP,IQ,IR,IS)) = Work1(IPQ)
C        Write(6,*) IP,IQ,IR,IS, Work1(IPQ)
C        Write(6,*) IP,IQ,IR,IS, TwoNO(NAddr3(IP,IQ,IR,IS))
        EndDo
        EndDo
      End If
C
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

      Call reader%close

      End

*Deck LoadSaptTwoNO
      Subroutine LoadSaptTwoNO(Mon,TwoNO,NBasis,NInte2)
C
      use sorter
C
      Implicit Real*8 (A-H,O-Z)
C
      Type(AOReaderData) :: reader
C
      Integer :: Mon, NInte2, NBasis
      Integer :: IRS, IS, IR, IPQ, IQ, IP
      Integer :: iunit
      Integer(8),external :: NAddr3
      Dimension :: TwoNO(NInte2), Work1(NBasis**2)
      Logical :: empty
      Character*9 :: fname

      If(Mon.Eq.1) Then
      fname='TWOMOAA '
      ElseIf(Mon.Eq.2) Then
      fname='TWOMOBB '
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
C
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
C
      EndDo
      EndDo

      close(iunit)

      End

*Deck CheckSaptTwoEl
      Subroutine CheckSaptTwoEl(Mon,TwoNO,NBasis,NInte2)
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
      TMP=Abs(TwoNO(NAddr3(IP,IQ,IR,IS))-Work1(IPQ))
      If(TMP>1.D-6) Write(*,*) IP,IQ,IR,IS,
     $ TwoNO(NAddr3(IP,IQ,IR,IS)),Work1(IPQ)
C      Write(6,*) IP,IQ,IR,IS, Work1(IPQ)
C      Write(6,*) IP,IQ,IR,IS, TwoNO(NAddr3(IP,IQ,IR,IS))
      EndDo
      EndDo
      EndDo
      EndDo

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
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5) Then
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

      subroutine TwoEneGVBChck(ETwo,Occ,NOccup,NBasis)
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
      Double Precision, Allocatable :: RDM2val(:,:,:,:),
     $                                 work(:),ints(:,:)
      Character(:),Allocatable :: IntJFile
C
C     SET FILES
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5) Then
      IntJFile='FFOO'
      Else
      IntJFile='FFOOERF'
      EndIf
C
      Allocate(work(NBasis**2),ints(NBasis,NBasis))
      Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
C
      Do L=1,NOccup
      Do K=1,NOccup
      Do J=1,NOccup
      Do I=1,NOccup
      RDM2val(I,J,K,L) = FRDM2GVB(I,K,J,L,Occ,NBasis)
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
      end subroutine TwoEneGVBChck

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
      If (IFunSR.Eq.0.Or.IFunSR.Eq.3.Or.IFunSR.Eq.5) Then
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

*Deck SortOrbDal
      Subroutine SortOrbDal(URe1,Occ2,NNIn,NNAct,NSym,IOrbSym,NBasis)
C     sorts the orbitals in URe1 so that the
C     inactive orbitals go first, then active, and secondary orbitals
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 URe1(NBasis*NBasis),URe2(NBasis*NBasis),
     $ Occ1(NBasis),Occ2(NBasis)
      Dimension IOrbSym(8),IActOrb(NBasis),InActOrb(NBasis)
      Dimension LabelAct(NBasis),LabelIAct(NBasis),
     $ ICpy1(NBasis),ICpy2(NBasis)
      Logical FileOcc,FileSIRIFC
C
      Occ2(1:NBasis)=0.0
      Occ1(1:NBasis)=0.0
C
      Inquire(file="occupations.dat",EXIST=FileOcc)
      Inquire(file="SIRIFC",EXIST=FileSIRIFC)
      If (FileOcc) Then
C        read occupations from Dalton output...
         Open(10,File="occupations.dat",Form='Formatted',Status='Old')
         Read(10,*)NNIn,NNAct
         NNIn=NNIn/2
         Read(10,*) (Occ1(I),I=1,NNAct+NNIn)
         Read(10,*) (IActOrb(I),I=1,NSym)
         Read(10,*) (InActOrb(I),I=1,NSym)
         Close(10)
      ElseIf(FileSIRIFC) Then
C        read 1rdm from SIRIFC file...
         Call read_1rdm_dalton(Occ1,IActOrb,InActOrb,NNIn,NNAct,NBasis)
      Else
         Write(6,'(1x,a)') 'Occupation numbers from Dalton not found!'
         Stop
      EndIf ! FileOcc
C
      II=0
      Do I=1,NSym
      Do J=1,IOrbSym(I)
      II=II+1
      LabelAct(II)=0
      If(J.Gt.InActOrb(I).And.J.Le.IActOrb(I)+InActOrb(I))
     $ LabelAct(II)=1
      LabelIAct(II)=0
      If(J.Le.InActOrb(I)) LabelIAct(II)=1
      EndDo
      EndDo
C
      Do I=1,NNIn+NNAct
      ICpy1(I)=0
      ICpy2(I)=0
      EndDo
C
      Do II=1,NNIn+NNAct
C
      Do I=1,NNIn+NNAct
      If(ICpy2(I).Eq.0.And.ICpy1(II).Eq.0.And.Occ1(I).Eq.2.0D0) Then
      ICpy2(I)=1
      ICpy1(II)=1
      Occ2(II)=Occ1(I)
      EndIf
      EndDo
C
      If(ICpy1(II).Eq.0) Then
      Do I=1,NNIn+NNAct
      If(ICpy2(I).Eq.0.And.ICpy1(II).Eq.0) Then
      ICpy2(I)=1
      ICpy1(II)=1
      Occ2(II)=Occ1(I)
      EndIf
      EndDo
      EndIf
C
      EndDo
C
      Do I=1,NBasis
      ICpy1(I)=0
      ICpy2(I)=0
      EndDo
C
      Do II=1,NBasis
C
      Do I=1,NBasis
      If(LabelIAct(I).Eq.1.And.ICpy2(I).Eq.0.And.ICpy1(II).Eq.0) Then
      ICpy2(I)=1
      ICpy1(II)=1
      Do J=1,NBasis
      URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
      EndDo
      EndIf
      EndDo
C
      If(ICpy1(II).Eq.0) Then
      Do I=1,NBasis
      If(LabelAct(I).Eq.1.And.ICpy2(I).Eq.0.And.ICpy1(II).Eq.0) Then
      ICpy2(I)=1
      ICpy1(II)=1
      Do J=1,NBasis
      URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
      EndDo
      EndIf
      EndDo
      EndIf
C
      If(ICpy1(II).Eq.0) Then
      Do I=1,NBasis
      If(ICpy2(I).Eq.0.And.ICpy1(II).Eq.0) Then
      ICpy2(I)=1
      ICpy1(II)=1
      Do J=1,NBasis
      URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
      EndDo
      EndIf
      EndDo
      EndIf
C
C     Do II
      EndDo
C
      Do I=1,NBasis*NBasis
      URe1(I)=URe2(I)
      EndDo
C
      Return
      End

      Subroutine read_1rdm_dalton(Occ,IActOrb,InActOrb,
     $                             NNIn,NNAct,NBasis)
C
C     Purpose: read 1RDM from from Daton together with
C              the number of inactive and active orbitals in each irrep
C
C     Dalton naming convention
C     NASHT  -- no of active orbitals
C     NISHT  -- no of inactive orbitals
C     NNASHX -- triang of active orbitals
C     DV(1:NNASHX) - triangular 1RDM
C
      use tran
      implicit none

      integer,intent(in)  :: NBasis
      integer,intent(out) :: NNIn,NNAct
      integer,intent(out) :: IActOrb(NBasis),InActOrb(NBasis)
      double precision,intent(out) :: Occ(NBasis)

      integer :: isirifc
      integer :: NISHT,NASHT,NNASHX,NSYM
      integer :: MULD2H(8,8),NISH(8),NASH(8)
      integer :: N2ASHX,DUMMY,NDUM(8),HlpDim
      integer :: I
      double precision :: DV(1:NBasis**2)
      double precision,allocatable :: OneAct(:,:),EigAct(:)
      double precision,allocatable :: work(:)


      ! read 1RDM in active orbs from SIRIFC file
      open(newunit=isirifc,file='SIRIFC',status='OLD',
     &     access='SEQUENTIAL',form='UNFORMATTED')
      read(isirifc)
      read(isirifc)
      read(isirifc) NISHT,NASHT,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &              DUMMY,DUMMY,NNASHX,DUMMY,DUMMY,DUMMY,
     &              NSYM,MULD2H,NDUM,NDUM,
     &              NISH,NASH

      read(isirifc)
      read(isirifc)
      read(isirifc) DV(1:NNASHX) ! how to correctly read DS?

C      print*, 'DV',DV(1:NNASHX)
      close(isirifc)

      allocate(EigAct(NASHT))
c      print*, 'Active occ numbers'
      Do I=1,NASHT
         EigAct(I) = DV((I-1)*I/2+I)
C        print*, I,EigAct(i)
      EndDo

C     test if DV diagonal
      If (abs(norm2(EigAct)-norm2(DV)).gt.1D-8) Then

         write(LOUT,'(1x,a)') "Diagonalize DV from Dalton..."
         N2ASHX = NASHT*NASHT
         HlpDim = max(N2ASHX,3*NASHT)
         allocate(OneAct(NASHT,NASHT),work(HlpDim))
         EigAct = 0
         Call triang_to_sq2(DV,OneAct,NASHT)
         Call Diag8(OneAct,NASHT,NASHT,EigAct,work)

         deallocate(OneAct,work)

      End If

C     transfer to GammCor variables
      NNIn  = NISHT
      NNAct = NASHT
      IActOrb(1:NSym)  = NASH(1:NSYM)
      InActOrb(1:NSym) = NISH(1:NSYM)
      Occ = 0d0
      Occ(1:NISHT) = 2d0
      Occ(NISHT+1:NISHT+NASHT) = EigAct(1:NASHT)

      End

      Subroutine ReadCGemDal(CICoef,NELE,INActive,NActive,NBasis)
C
C     Purpose: read geminal coefficients either from Dalton output (coeff.dat)
C              or SIRIFC file
C
      implicit none
C
      integer,intent(in)  :: NELE,NBasis
      integer,intent(out) :: INActive,NActive
      double precision,intent(out) :: CICoef(NBasis)
C
      integer :: i
      logical :: FileCoeff, FileSIRIFC

      Inquire(file="coeff.dat", EXIST=FileCoeff)
      Inquire(file="SIRIFC",EXIST=FileSIRIFC)
C
      If (FileCoeff) Then
C        read coefficients from Dalton output...
         CICoef(1:NBasis)=0.D0
         Open(10,File='coeff.dat',Form='Formatted',Status='Old')
         Read(10,*) NActive
         INActive=NELE-NActive
         Do I=1,INActive
         CICoef(I)=1.D0
         EndDo
         Read(10,*) (CICoef(I+INActive),I=1,2*NActive)
         Close(10)
      ElseIf(FileSIRIFC) Then
C        read coefficients from Dalton SIRIFC file...
         call read_CGEM_dalton(CICoef,INActive,NActive,NBasis)
      Else
         Write(6,'(1x,a)') 'Geminal coefficients from Dalton not found!'
         Stop
      EndIf
C
      End

      Subroutine read_CGEM_dalton(CICoef,INActive,NActive,NBasis)
C
C     Purpose: read geminal coefficients from Dalton SIRIFC file
C              set the number of (in)active geminals
C
      use types
      implicit none

      integer,intent(in)  :: NBasis
      integer,intent(out) :: INActive,NActive
      double precision,intent(out) :: CICoef(NBasis)
C
      integer :: isirifc
      integer :: NGEM,NISHT_G,NASHT_G
      double precision :: CGEM(NBasis)

      CGEM(1:NBasis) = 0.D0
      CICoef(1:NBasis)=0.D0

      open(newunit=isirifc,file='SIRIFC',status='OLD',
     &     access='SEQUENTIAL',form='UNFORMATTED')

      call readlabel(isirifc,'CI+APSG ')
      read(isirifc) NGEM,NISHT_G,NASHT_G
      read(isirifc) CGEM(1:NASHT_G)

C     print*, 'NGEM',   NGEM
C     print*, 'NISHT_G',NISHT_G
C     print*, 'NASHT_G',NASHT_G
      close(isirifc)

      ! set no of (in)active geminals
      INActive = NISHT_G
      NActive  = NASHT_G / 2

      CICoef(1:INActive) = 1.d0
      CICoef(INActive+1:INActive+NASHT_G) = CGEM(1:NASHT_G)
      ! print*,'CICoef-1',CICOef(1:NBasis)

      End

