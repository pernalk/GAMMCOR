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
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
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
      
C      Do I=1,NDimX 
C      Write(6,*) I,IndN(1,I),IndN(2,I)
C      EndDo

C  
      print*, 'ACCAS!'     
      ACAlpha=0.000001 
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)

      TMP=0d0
C      Do I=1,NDimX**2
C      TMP = TMP + ABPLUS(I)**2d0
C      TMP = TMP + ABMIN(I)**2d0
C      EndDo
C      Write(*,*) TMP
C       Do I=1,NDimX
C       Write(*,*) ABPLUS(NDimX*(I-1)+I), ABMIN(NDimX*(I-1)+I)
C       EndDo

C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,X,"***************************** ")')
      Write(6,'(  X,"*** ERPA-CAS CALCULATIONS *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)

!     CHECK FOR SAPT
!      Do I=1,NBasis
!       Write(*,*) CICoef(I)
!      EndDo
C      TMP = 0
C      Do I=1,NDimX**2
C      TMP = TMP + EigVecR(I)**2
C      EndDo
C      Write(6,*) TMP
C
C      Do i=1,NDimX
C       write(*,*) IndN(1,i)
C       write(*,*) IndN(2,i)
C      EndDo
C
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
