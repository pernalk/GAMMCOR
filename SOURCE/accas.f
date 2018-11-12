*Deck ACCAS
      Subroutine ACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NGem)
C
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
C     DENSITY MATRIX ELEMENTS
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ IndX(NBasis*(NBasis-1)/2),IndN(2,NBasis*(NBasis-1)/2),
     $ IndAux(NBasis),IPair(NBasis,NBasis)
C
C     CONSTRUCT LOOK-UP TABLES
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      ICount=0
C
      Do I=1,NBasis
C
      If(Occ(I).Lt.One.And.Occ(I).Ne.Zero) Then
      IndAux(I)=1
      Write(6,'(X," Active Orbital: ",I4,E14.4)') I, Occ(I)
      ICount=ICount+1
      EndIf
      EndDo
C
      Write(6,'(/,X," In ACCAS: Active Orbitals ",I4,/)')ICount

C     

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
      If( (IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.( Abs(Occ(I)-Occ(J))/Occ(I).Lt.1.D-10 ) 
c     $ .Or.( I.Le.NELE.And.J.Le.NELE.
c     $          And.Occ(I).lt.One.And.Occ(J).Lt.One )
     $)  Then
C
      Write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals  
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
C
      Ind=Ind+1
      IndX(Ind)=Ind
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      IPair(I,J)=1
      IPair(J,I)=1
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
      Write(6,'(2X,"Number of pairs reduced to:",I6)')Ind
      Write(6,'(2X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
      Call RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      Return
      End


* Deck RunACCAS 
      Subroutine RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      use sorter
      use tran    
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
C     (0) srDFT
C
      Real*8, Dimension(:), Allocatable :: TwoElErf
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: XKer
      Dimension NSymNO(NBasis),VSR(NInte1),MultpC(15,15)
      Parameter (Four=4.D0)
C
C     (0) END OF srDFT
C 
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem), EGOne(NGem)
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C     IFlSnd  = 1 - run AC0-CAS (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0-CAS
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1)
     $ Write(6,'(/," *** ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1) 
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1) 
     $ included in ERPA correlation ***",/)')
      EndIf
C
C     CALL AC If IFlAC=1 OR IFlSnd=1
C
      If(IFlAC.Eq.1.Or.IFlSnd.Eq.1) Then
      NGOcc=0
      Call ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      Return
      EndIf
C
C     CALCULATE THE A+B AND A-B MATRICES
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      ACAlpha=One
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
C     (1) BEGINNING OF srDFT
C     ADD CONTRIBUTIONS FROM srKS POTENTIAL 
C
C     IFunSR = 1  : LDA
C              2  : PBE
C
      IFunSR=1
      If(IFunSR.Gt.0) Then
C
C     load NGrid
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,X," The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (TwoElErf(NInte2))
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NBasis*NGrid))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
C
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
C     set/load Alpha - range separation parameter
C
      Alpha=0.5d0
C
C     set/load long-range two-electron integrals and transform to NO
C     TwoElErf MUST BE IN NO !!!!
C
      TwoElErf(1:NInte2)=Zero
      call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT')
      call LoadSaptTwoEl(4,TwoElErf,NBasis,NInte2)
      Write(6,'(/,2x,a)') 'Transforming two-electron LR integrals ...'
      Call TwoNO1(TwoElErf,UNOAO,NBasis,NInte2)
C      print*, norm2(TwoElErf)
C      print*, norm2(TwoNO)
C
C     set/load symmetries of NO's, compute sr potential (vsr=xc+hartree)
C     as a byproduct a sr energy (ensr=sr-xc+sr-hartree) is obtained
C
      NSymNO(1:NBasis)=1
      Call EPotSR(EnSR,VSR,Occ,URe,UNOAO,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NSymNO,TwoNO,TwoElErf,
     $ Alpha,IFunSR,NGrid,NInte1,NInte2,NBasis)
      Do I=1,NInte1
      XOne(I)=XOne(I)+VSR(I)
      EndDo
      Write(6,'(/,"  SR Energy",F15.8)')EnSR
C
C     compute AB matrices with the lr integrals and a modified potential
C
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoElErf,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
      EndIf
C
C     (1) END OF srDFT
C
      If(IFunSR.Eq.0)
     $ Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
C     (2) BEGINNING OF srDFT
C     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
C
      If(IFunSR.Gt.0) Then
C
      Write(6,'(/,"*** Adding a sr-kernel... ***")')
      NDimKer=NBasis*(1+NBasis)*(2+NBasis)*(3+NBasis)/24
      Allocate (XKer(NDimKer))
C     set/load the group multiplication table 
C
      MultpC(1,1)=1
C
C     uncomment Call RhoKernel(RhoVec,SRKer,Alpha,NGrid) in GetKerNO and
C     remove Stop
C     see other changes in GetKerNO
      Call GetKerNO(XKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NDimKer,NBasis,NGrid)
C
      Do I=1,NBasis
      CICoef(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CICoef(I)=-CICoef(I)
      EndDo
C
      Do IRow=1,NDimX
C
      IA=IndN(1,IRow)
      IB=IndN(2,IRow)
      IAB=IndX(IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NDimX
C
      IC=IndN(1,ICol)
      ID=IndN(2,ICol)
      ICD=IndX(ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      XKer1234=XKer(NAddrrK(IA,IB,IC,ID))
      TwoSR=TwoNO(NAddr3(IA,IB,IC,ID))-TwoElErf(NAddr3(IA,IB,IC,ID))
C
      ABMIN((ICol-1)*NDimX+IRow)=ABMIN((ICol-1)*NDimX+IRow)
     $ +Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      EndDo
      EndDo
C
      Write(6,'("*** sr-kernel added. ***")')
C
      EndIf
C
C     (2) END OF srDFT
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,X,"***************************** ")')
      Write(6,'(  X,"*** ERPA-CAS CALCULATIONS *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
C     (3) srDFT
C
      If(IFunSR.Gt.0) Then
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Write(6,'(/," *** LR-CAS-SR-DFT Excitation Energies *** ",/)')
      Do I=1,10
      Write(6,'(I4,4X,E16.6)') I,Eig(I)
      EndDo
      Deallocate(OrbZGrid,OrbYGrid,OrbXGrid,OrbGrid,WGrid)
      Deallocate(TwoElErf)
      Stop
      EndIf
C
C     (3) END OF srDFT
C
C      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)   
      If(IFunSR.Eq.0)
     $ Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
      Write(6,'(/," *** Computing ERPA energy *** ",/)')
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
      ECorr=Ecorr*Half
C
      Write
     $ (6,'(/,1X,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
C

      Return
      End
 
