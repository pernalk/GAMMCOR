*Deck DMSCF 
      Subroutine DMSCF 
     $ (Title,URe,Occ,XKin,XNuc,ENuc,UMOAO,TwoEl,NBasis,NInte1,NInte2,
     $ NGem)
C
C     !!! XKin CONTAINS BOTH KINETIC AND EL-N CONTRIBUTIONS !!!
C     !!! XNuc IS EMPTY 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
C
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(NInte1),TwoEl(NInte2),UMOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension XOne(NInte1),UReSav(Nbasis,NBasis)
C
      Do I=1,NInte1
      XOne(I)=XKin(I)
      EndDo
C
C     PRINT FINAL NO's IN AO's REPRESENTATION 
C
c      Write(6,'(/,X,"NATURAL ORBITALS IN AO BASIS SET")')
cC
c   98 Format(X,10F10.6)
      Call MultpM(UReSav,URe,UMOAO,NBasis)
c      NLine=NBasis/10
c      If(NLine*10-NBasis.Ne.0)NLine=NLine+1
c      Do I=1,NBasis
c      Write(*,'(I3)') I      
cC
c      Do LL=0,NLine-1
c
c      NN=NBasis-10*LL
c      If(NN.Le.10) Then
c      Write(*,98) (UReSav(I,K),K=10*LL+1,NBasis)
c      Else
c      Write(*,98) (UReSav(I,K),K=10*LL+1,10*(LL+1))
c      EndIf
cC
c      EndDo
c      Write(*,*)
c      EndDo
C
      Write(6,'(/,X,'' NATURAL ORBITAL OCCUPANCIES'')')

      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"CICoef",7X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,2E16.6,I6)') I,Occ(I),CICoef(I),IGem(I)
      EndDo
      Write(6,'()')
C
      NDim=NBasis*(NBasis-1)/2
      NDimKer=NBasis*(1+NBasis)*(2+NBasis)*(3+NBasis)/24
C
      If(ICASSCF.Eq.1) Then
C
      Call ACCAS(ETot,ENuc,TwoEl,URe,UReSav,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NGem)
C
      Else
C
C     MH: LATER REPLACE OPTNORB BY ENERGY CHECK!!!
      If(ITwoEl.Eq.1) Then
      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,TwoEl,UMOAO,
     $ NSymMO,Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG)
      EndIf
C
      NGOcc=0
      Call INTERPA(ETot,ENuc,TwoEl,URe,UReSav,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NDim,NGem,NGOcc)
C
      EndIf
C
c      NoEig=3
c      Call ACPINO(ENuc,TwoEl,Occ,XOne,NBasis,NInte1,NInte2,NDim,NGem,
c     $ NoEig)
C
      Return
      End

*Deck Partition
      Subroutine Partition(ULoc,NGOcc,NGem,NBasis)
C     
      Implicit Real*8 (A-H,O-Z)
C     
C     DIVIDE A SET OF OCCUPIED LOCALIZED MO ORBITALS INTO ACTIVE (FRAGMENT A) AND INACTIVE ONES
C     VIRTUALS GO TO THE ACTIVE A, ASSIGN ORBITALS TO GEMINALS
C
      Include 'commons.inc'
C     
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C    
C     UMOAO - A TRANSFORMATION MATRIX FROM AO TO LOCALIZED MO
C 
      Dimension ULoc(NBasis,NBasis)
C    
C     LOCAL ARRAYS
C
C     NAt(I) - AN INDEX OF AN ATOMS ON WHICH THE ITH ORBITALS IS CENTERED
C     NBasisAt(I) - A NUMBER OF ATOMIC ORBITALS CENTERED ON THE ITH ATOM
C     NPair(I,1/2) - INDICES OF THE FIRST AND SECOND ATOMS ON WHICH SPECIES I IS LOCALIZED   
C     IOrbA(I) = 1 - ITH ORBITAL IS ACTIVE
C     IOrbA(I) = 0 - ITH ORBTIAL IS INACTIVE 
C     NoAt - THE NUMBER OF ATOMS IN A MOLECULE
C
      Dimension NAt(NBasis),NBasisAt(5000),NPair(5000,2)
C    
C     BEGINNING OF THE INPUT
C c herer!!! - change NBasisAt if basis set changes
      NoAt=3
      NBasisAt(1)=15
      NBasisAt(2)=5
      NBasisAt(3)=5
C
C     DEFINE AN ACTIVE SYSTEM AS A SET OF ACTIVE BONDS AND ATOMS
C     MAKE SURE THAT NPair(1,I) > NPair(2,I)
C
C     H2O : two bonds are active
C
      NoSpecies=2
C
C     OH1 BOND
      NPair(1,1)=2
      NPair(1,2)=1
C     OH2 bond
      NPair(2,1)=3
      NPair(2,2)=1
C
C     END OF INPUT
C
      Do I=1,NELE
      IOrbA(I)=0
      EndDo
C
C     ALL VIRTUALS WILL BE ACTIVE
C
      Do I=NELE+1,NBasis
      IOrbA(I)=1
      EndDo
C
C     CHECK IF THE NUMBER OF AO's SUM UP TO NBasis
C
      NB=0
      Do I=1,NoAt
      NB=NB+NBasisAt(I)
      EndDo
      If(NB.Ne.NBasis)
     $ Stop 'Fatal Error: No of atomic orbitals different from NBasis!'
C
C     SET A NUMBER OF AN ATOM ON WHICH AN AO IS CENTERED
C
      ICount=1
      Do I=1,NoAt
C
      Do J=1,NBasisAt(I)
      NAt(ICount)=I
      ICount=ICount+1
      EndDo
C
      EndDo
C
C     ASSIGN OCC ORBITALS TO AN ACTIVE OR INACTIVE SPACE
C
      NGemA=0
      Do I=1,NELE
C
      ISpec1=0
      ISpec2=0
C
      Do K=1,NoAt
      AtContr=Zero
      Do J=1,NBasis
      If(NAt(J).Eq.K) AtContr=AtContr+ULoc(I,J)**2
      EndDo
C
      If(AtContr.Gt.0.1D0) Then
      If(ISpec1.Eq.0) Then
      ISpec1=K
      ElseIf(ISpec2.Eq.0) Then
      ISpec2=ISpec1
      ISpec1=K
      Else
      Stop 'Problem with selecting active orbitals!'
      EndIf
      EndIf
C
      EndDo
C
      Do J=1,NoSpecies
      If(NPair(J,1).Eq.ISpec1.And.NPair(J,2).Eq.ISpec2) Then
      IOrbA(I)=1
      NGemA=NGemA+1
      EndIf
      EndDo
C
      If(IOrbA(I).Eq.1) Write(6,*) "ORBITAL NO",I," CENTERED ON ATOMS ",
     $ ISpec1,ISpec2,"IS ACTIVE"
      If(IOrbA(I).Eq.0) Write(6,*) "ORBITAL NO",I," CENTERED ON ATOMS ",
     $ ISpec1,ISpec2,"IS INACTIVE"
C
      EndDo
C
C     ASSIGN INACTIVE ORBITALS TO THE FIRST NGem-NGemA GEMINALS
C
      NGem=NELE
      NGemB=NGem-NGemA
      NGOcc=NGemB
      NActive=NBasis-NGemB
      Write(6,*) 'The number of active orbitals is', NActive
C
      ICount=0
      Do I=1,NBasis
      If(IOrbA(I).Eq.0) Then
      ICount=ICount+1
      IGem(I)=ICount
      EndIf
      EndDo
C
C     ASSIGN ACTIVE ORBITALS TO THE LAST NGemB+1 .... NGem GEMINALS
C
      ICount=NGemB
      Do I=1,NBasis
      If(IOrbA(I).Eq.1) Then
      ICount=ICount+1
      IGem(I)=ICount
      EndIf
      If(ICount.Eq.NGem) ICount=NGemB
      EndDo
C
      Return
      End


*Deck MoldenCAS
      Subroutine MoldenCAS(Occ,UAONO,NBasis)
C
C     writes NO orbitals to cas_ss.molden file, created from the 
C     existing (!) cas.molden file with AOs and UAONO
C     useful when sa-cas is run in molpro (molpro prints to molden 
C     state-averaged cas orbitals 
C
C     WARNING! it is absolutely necessary to use "CARTESIAN" AO's in molpro
C     (in molden cartesian gaussians are assumed)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
      Include 'commons.inc'
C 
      Character*60 FName,FName2
      Character*100 Line,Aux1
      Character*10 Str
      Logical EX
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
      Dimension Occ(NBasis),UAONO(NBasis,NBasis)
C
C     OPEN cas.molden file
C
      Title="cas"
C
      Do I=1,60
      FName(I:I)=' '
      FName2(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      FName2(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+7)='.molden'
      FName2(K:K+11)='_ss.molden'
C
      INQUIRE(file=FName,EXIST=EX)
      If(EX) Then
      Write(6,'(/," ** cas_ss.molden will be created with NOs **",/)')
      Else
      Return
      EndIf
C      
      Open(10,File=FName)
      Open(20,File=Fname2)
C
    2 Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:4).Eq."[MO]") Then
      GoTo 20
      EndIf
      GoTo 2
C
   20 Continue
C      
      Do I=1,NBasis
C
    7 Read(10,'(A100)',End=333)Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:6).Eq." Spin=") Then
      GoTo 27
      EndIf
      GoTo 7
   27 Continue
C
      Read(10,'(A100)')Aux1
      Write(20,'(" Occup=",F12.6)') Two*Occ(I)
C
      Do J=1,NBasis
      Read(10,'(A100)')Aux1
      If(J.Lt.10) Write(20,'(I1,F21.14)')J,UAONO(I,J)
      If(J.Ge.10.And.J.Lt.100) Write(20,'(I2,F21.14)')J,UAONO(I,J) 
      If(J.Ge.100.And.J.Lt.1000) Write(20,'(I3,F21.14)')J,UAONO(I,J) 
      EndDo
C
      EndDo 
C
    9 Read(10,'(A100)',End=333)Aux1
      Write(20,'(A100)') Aux1
      GoTo 9
  333 Continue
C
      Close(10)
      Close(20)
C
      Return
      End

*Deck MoldenPrep
      Subroutine MoldenPrep(Title,Occ,UAONO,NGem,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
      Include 'commons.inc'
C 
      Character*60 FName,FName2
      Character*100 Line,Aux1
      Character*10 Str
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
      Dimension Occ(NBasis),UAONO(NBasis,NBasis)
C
C     OPEN title.molden file
C
      Do I=1,60
      FName(I:I)=' '
      FName2(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      FName2(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+7)='.molden'
      FName2(K:K+11)='_all.molden'
C      
      Open(10,File=FName)
      Open(20,File=Fname2)
C
    2 Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:4).Eq."[MO]") Then
      GoTo 20
      EndIf
      GoTo 2
C
   20 Continue
C      
      Do I=1,NBasis
C
      Do J=1,3
      Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      EndDo
C
      Read(10,'(A100)')Aux1
      Write(20,'(" Occup=",F12.6)') Two*Occ(I)
C
      Do J=1,NBasis
      Read(10,'(A100)')Aux1
      Write(20,'(I3,F18.14)')J,UAONO(I,J)
      EndDo
C
      EndDo 
C
      Close(10)
      Close(20)
C
C     NOW WRITE GEMINALS TO SEPARATE MOLDEN FILES
C
      Do IG=1,NGem
C
      Do I=1,60
      FName2(I:I)=' '
      EndDo  

      Write(Str,'(I10)')IG
      If(IG.Lt.10) FName2=Trim(Title)//'_g'//Str(10:10)//'.molden'
      If(IG.Ge.10) FName2=Trim(Title)//'_g'//Str(9:10)//'.molden'
C
      Open(10,File=FName)
      Open(20,File=Fname2)
C
    3 Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      If(Aux1(1:4).Eq."[MO]") Then
      GoTo 30
      EndIf
      GoTo 3
C
   30 Continue
C
      Do I=1,NBasis
C
      Do J=1,3
      Read(10,'(A100)')Aux1
      Write(20,'(A100)') Aux1
      EndDo
C
      Read(10,'(A100)')Aux1
      If(IGem(I).Eq.IG) Then
      Write(20,'(" Occup=",F12.6)') Two*Occ(I)
      Else
      Write(20,'(" Occup=",F12.6)') Zero
      EndIf
C
      Do J=1,NBasis
      Read(10,'(A100)')Aux1
      Write(20,'(I3,F18.14)')J,UAONO(I,J)
      EndDo
C
      EndDo
C
      Close(20)
      Close(10)
C
      EndDo
C
      Return
      End

*Deck PrintAPSG
      Subroutine PrintAPSG(Title,UMOAO,UNOAO,URe,Occ,
     $           XOne,TwoEl,NGem,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName,Aux1
      Include 'commons.inc'
C 
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
      Dimension UMOAO(NBasis,NBasis),UNOAO(NBasis,NBasis),
     $          URe(NBasis,NBasis),Occ(NBasis),
     $           XOne(NInte1),TwoEl(NInte2)
C
C     LOCAL
C
      Dimension Work(NBasis,NBasis),Work1(NBasis,NBasis),
     $ Gamma(NInte1),Work2(NBasis,NBasis)
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
      FName(K:K+13)='_APSG_all.dat'
C      
      Open(10,File=FName)
C
C     PRINT THE OVERLAP MATRIX S
C
      Open(20,File='s.dat')

      Read(20,'(A10)')Aux1
      Read(20,'(A10)')Aux1
C
      Do I=1,NBasis
      Do J=1,NBasis/5
      Read(20,222) (Work(I,(J-1)*5+K),K=1,5)
  222 Format(5(F15.8,1X))
      EndDo
      MxK=NBasis-(NBasis/5)*5
      IStart=(NBasis/5)*5
      If(MxK.Ne.0) Read(20,222) (Work(I,IStart+K),K=1,MxK)
      EndDo

      Close(20)
C
      Write(10,'(/,X,"*AO OVERLAP MATRIX S")')
      Do I=1,NBasis
      Write(10,*) (Work(I,J),J=1,NBasis)
      EndDo 
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO MO
C
      Write(10,'(/,X,"*MOLECULAR ORBITALS IN AO BASIS SET")')
      Write(10,'(X,"*MOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (UMOAO(I,J),J=1,NBasis)
      EndDo
C
C     TEST0     UMOAO*S*UMOAO^T = 1
c      Do I=1,NBasis
c      Do J=1,NBasis
c      Sum=0.0D0
c      Do K=1,NBasis
c      Do L=1,NBasis
c      Sum=Sum+UMOAO(I,K)*Work(K,L)*UMOAO(J,L)
c      EndDo
c      EndDo
c      Write(*,*)I,J,Sum
c      EndDo
c      EndDo
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO NO
C
      Write(10,'(/,X,"*NATURAL ORBITALS IN AO BASIS SET")')
      Write(10,'(X,"*NOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (UNOAO(I,J),J=1,NBasis)
      EndDo
C
C     PRINT TRANSFROMATION MATRIX FROM AO TO NO
C     
      Write(10,'(/,X,"*NATURAL ORBITALS IN MO BASIS SET")')
      Write(10,'(X,"*NOs IN ROWS")')
      Do I=1,NBasis
      Write(10,*) (URe(I,J),J=1,NBasis)
      EndDo
C 
C     TEST1 UNOAO=URe*UMOAO
C
c      Call MultpM(Work,URe,UMOAO,NBasis)
c      Call DiffM(Work,UNOAO,Work1,NBasis)
c      write(*,*)work1
C
C    PRINT Occ, CICOef's, AND ASSIGMENT TO GEMINALS 
C
      Write(10,'(/,X,"*NATURAL OCCUPATION NUMBERS")')
      Write(10,'(X,
     $ "*IN THE RANGE [0,2]. NORMALIZED TO THE NO OF ELECTRONS")')
      Write(10,*)(2.0D0*Occ(I),I=1,NBasis)
C
      Write(10,'(/,X,"*EXPANSION COEFFICIENTS FOR GEMINALS")')
      Write(10,
     $ '(X,"*(c_i)^2=n_i/2. NORMALIZED TO 1 FOR EACH GEMINAL")')
      Write(10,*)(CICoef(I),I=1,NBasis)
C
      Write(10,'(/,X,"*ASSIGNMENT OF ORBITALS TO GEMINALS")')
      Write(10,
     $ '(X,"*INDICES OF GEMINALS TO WHICH ORBITALS BELONG TO")')
      Write(10,*)(IGem(I),I=1,NBasis)
C
C     PRINT 1RDM IN MO
C
      Call VecTr(Gamma,2.0D0*Occ,URe,NBasis)
      Call CpySym(Work2,Gamma,NBasis)
      Write(10,'(/,X,"*1RDM IN MO")')
      Write(10,'(X,"*1RDM NORMALIZED TO THE NUMBER OF ELECTRONS")')
      Do I=1,NBasis
      Write(10,*) (Work2(I,J),J=1,NBasis)
      EndDo
C
C     PRINT 1RDM IN AO
C
      Do I=1,NBasis  
      Do J=1,I
      Work(I,J)=UMOAO(J,I)
      Work(J,I)=UMOAO(I,J)
      EndDo
      EndDo
      Call MatTr(Gamma,Work,NBasis)
      Call CpySym(Work,Gamma,NBasis)
      Write(10,'(/,X,"*1RDM IN AO")')
      Do I=1,NBasis
      Write(10,*) (Work(I,J),J=1,NBasis)
      EndDo
C
C     TEST2 1RDM_AO=UNOAO^T*2*Occ*UNOAO
C
c      Call VecTr(Gamma,2*Occ,UNOAO,NBasis) 
c      IJ=0
c      Do I=1,NBasis
c      Do J=1,I
c      IJ=IJ+1
c      Write(*,*)I,J,Gamma(IJ)-Work(I,J)
c      EndDo
c      EndDo
C
C     COMPUTE AND PRINT 2-RDM IN NO
C
C     RDM2 IS IN THE NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
      Do I=1,60
      FName(I:I)=' '
      EndDo
      K=0
    6 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 6
      EndIf
      FName(K:K+14)='_APSG_2RDM.dat'
      Open(30,File=FName)
      Write(30,'(X,"2RDM IN THE NO REPRESENTATION")')
      Write(30,'(X,"P Q R S 2RDM_pqrs")')  
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,NBasis
      IPQ=IPQ+1
      IRS=0
      Do IR=1,NBasis
      Do IS=1,NBasis
      IRS=IRS+1
C
      RDM2=0.D0
      If(IP.Eq.IQ.And.IR.Eq.IS.And.IGem(IP).Eq.IGem(IR))
     $ RDM2=RDM2+CICoef(IP)*CICoef(IR)
      If(IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2+2.0D0*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      Write(30,'(4I3,E14.6)')IP,IQ,IR,IS,RDM2
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Close(10)
      Close(30)
C
      Return
      End
