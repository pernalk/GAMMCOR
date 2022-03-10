*Deck ERPA
      Subroutine ERPA(TwoMO,TwoMOLR,URe,Occ,XOne,
     $  OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,NDimKer,
     $  NGem,IAPSG,ISERPA,QMAX,Small)
C
C     A ROUTINE FOR EXCITATION ENERGY FROM LINEARIZED ERPA
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     $ Four=4.D0, Big=1.5D0)
     $ Four=4.D0, Big=1.98D0)
C
      Include 'commons.inc'
C
      Character*60 FMultTab
C
      Dimension
     $ URe(NBasis,NBasis),Occ(NBasis),TwoMO(NInte2),TwoMOLR(NInte2),
     $ XOne(NInte1),OrbGrid(NBasis,NGrid),WGrid(NGrid)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ CMAT(NDim*NDim),EMAT(NBasis*NBasis),EMATM(NBasis*NBasis),
     $ DMAT(NDim*NBasis),
     $ DMATK(NDim*NBasis),
     $ IndX(NDim),IndN(2,NDim),
     $ XKer(NDimKer),NSymMO(NBasis),NSymNO(NBasis),
     $ MultpC(15,15)
C
      Do I=1,NDim**2
      CMAT(I)=Zero
      EndDo
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
C     ESTABLISH THE SYMMETRY OF NO's
C
      Write(6,'(2/,X,"SYMMETRY OF NATURAL ORBITALS")')
      Write(6,'(X,"Orbital",3X,"Symmetry",3X,"Gem",3X,"Occupancy")')
      Do I=1,NBasis
      NSymNO(I)=0
C
      Do J=1,NBasis
C
      If(Abs(URe(I,J)).Gt.0.1) Then
C
      ISym=NSymMO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym) 
     $ Write(6,'("Symm of NO cannot be established",I3)')I
      EndIf
C
      EndIf
      EndDo
      Write(6,'(I4,7X,2I4,2X,E16.6)')I,NSymNO(I),IGem(I),Occ(I)
      EndDo
C
C     CALCULATE THE A+B AND A-B MATRICES
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      If(IAPSG.Eq.0) Then
C
      Write(6,'(/,X," ***** LINEARIZED ERPA CALCULATIONS *****")')
C
      Call LinAB_NEST(ABPLUS,ABMIN,URe,Occ,XOne,TwoMOLR,
     $ NBasis,NDim,NInte1,NInte2)
C
C     ABPLUS STORES AMAT, ABMIN STORES BMAT
C
c      IFlag=1
c      Call PrimAB_NEST(ABPLUS,ABMIN,URe,Occ,XOne,TwoMOLR,
c     $ NBasis,NDim,NInte1,NInte2,IFlag)
C
      Else
C
      If (ISERPA.Eq.1) Then
      Write(6,'(/,X," ***** SERPA CALCULATIONS *****")') 
      ElseIf(ISERPA.Eq.0) Then
      Write(6,'(/,X," ***** ERPA CALCULATIONS *****")')
      ElseIf(ISERPA.Eq.2) Then
      Write(6,'(/,X," ***** PINO CALCULATIONS *****")') 
      EndIf
C
      Call APSG_NEST(ABPLUS,ABMIN,CMAT,EMAT,EMATM,DMAT,DMATK,
     $ URe,Occ,XOne,TwoMOLR,
     $ NBasis,NDim,NInte1,NInte2,NGem,ISERPA)
      If(ISERPA.Eq.1.Or.ISERPA.Eq.0) Call CpyV(DMATK,DMAT,NDim*NBasis)
C
      EndIf
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      If(IFunSR.Ne.0) Then
C
      Write(6,'(/,"*** Adding a sr-kernel... ***")')
C
C     ADD A CONTRIBUTION FROM THE SR-XC KERNEL 
C
      Call GetKerNO(XKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NDimKer,NBasis,NGrid)
C
      IAB=0
      Do IA=1,NBasis
      CA=CICoef(IA)
C
      Do IB=1,IA
      CB=CICoef(IB)
C
      If(IB.Lt.IA) IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
      CC=CICoef(IC)
C
      Do ID=1,IC
      CD=CICoef(ID)
C
      XKer1234=XKer(NAddrrK(IA,IB,IC,ID))
      TwoSR=TwoMO(NAddr3(IA,IB,IC,ID))-TwoMOLR(NAddr3(IA,IB,IC,ID))
C
      If(ID.Lt.IC) ICD=ICD+1
C
      IABCD=(ICD-1)*NDim+IAB
      If(IA.Ne.IB.And.IC.Ne.ID)
     $ ABMIN(IABCD)=ABMIN(IABCD)+Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      If(ISERPA.Eq.1) Then
      IABC=(IC-1)*NDim+IAB
      If(IA.Ne.IB.And.IC.Eq.ID)
     $ DMATK(IABC)=DMATK(IABC)+Four*(CA+CB)*CC*(XKer1234+TwoSR)
      EndIf 
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'("*** sr-kernel added. ***")')
C
C     END OF THE CONDITION FOR IFunSR.Ne.0
C
      EndIf 
C
C     REDUCE DIMENSIONS OF THE MATRICES TO TAKE INTO ACCOUNT ONLY
C     NDimX ELEMENTS OF dGamma_ij 
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1  
      IJ=IJ+1
C
      If(J.Le.QMAX) Then
C
      If(Occ(I)+Occ(J).Lt.Big.And.Occ(I)+Occ(J).Gt.Small) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C      
      EndIf 
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
      Write(6,'(2/,X,"Small parameter set to:",E14.4)')Small
      Write(6,'(X,"Big   parameter set to:",E14.4)')Big
C
      Write(6,'(2/,X,"Total number of pairs:",I4)')NDim
      Write(6,'(X,"Reduced to:",I4)')Ind
      Write(6,'(X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
C     REDUCE THE MATRICES
C
      NDimN=NBasis
C
      If(IGVB.Eq.1) Then
      NDimN=2*(NGem-1)
      Write(6,'(2/,X,"NDimN set to:",I4)')NDimN
      EndIf 
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1) 
      CMAT(IJ)=CMAT(IJ1)
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(J-1)*NDim+IndX(I)
      DMAT(IJ)=DMAT(IJ1)
      DMATK(IJ)=DMATK(IJ1)
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(J-1)*NBasis+I
      EMAT(IJ)=EMAT(IJ1)
      EMATM(IJ)=EMATM(IJ1)
      EndDo
      EndDo
C
C     EXCITATION ENERGIES FROM THE ADIABATIC APPROX.
C
      If(ISERPA.Eq.1) Then
C
      Call SERPAExcit(ABPLUS,ABMIN,DMAT,DMATK,EMAT,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN)
C
      ElseIf(ISERPA.Eq.0) Then
C
C     ERPA
C
      Call RPAExcit(ABPLUS,ABMIN,CMAT,Occ,IndN,NSymNO,MultpC,
     $ NBasis,NDimX,IAPSG)
C
      ElseIf(ISERPA.Eq.2) Then
C
C     PINO
C
      Call PINOExcit(ABPLUS,ABMIN,DMAT,DMATK,EMAT,EMATM,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN)
C
      EndIf
C
      Return
      End

*Deck RPAExcit
      Subroutine RPAExcit(ABPLUS,ABMIN,CMAT,Occ,IndN,NSymNO,MultpC,
     $ NBasis,NDimX,IAPSG)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(SMALL=1.D-7,Delta=1.D-2)
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ CMAT(NDimX,NDimX),
     $ HlpAB(NDimX,NDimX),
     $ Q1(NDimX,NDimX),Q2(NDimX,NDimX),
     $ Occ(NBasis),
     $ Eig(NDimX),EigI(NDimX),Work(4*NDimX),
     $ IndN(2,NDimX),XR(NDimX*NDimX),
     $ NSymNO(NBasis),MultpC(15,15),
     $ BAS(NDimX,NDimX),IndS(NDimX)
C
      If(IAPSG.Eq.0) Then
C
C     ADD -2 CMAT TO A+B 
C
      Do I=1,NDimX
      Do J=1,NDimX
c      ABPLUS(I,J)=ABPLUS(I,J)-Two*CMAT(I,J)
      EndDo
      EndDo
C     
C     MULTIPLY A-B AND A+B BY N^-1 FROM THE LEFT
C     
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
C
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I)) Then
C
      Do K=1,NDimX
      Q1(IJ,K)=ABMIN(IJ,K)/(Occ(I)-Occ(J))
      Q2(IJ,K)=ABPLUS(IJ,K)/(Occ(I)-Occ(J))
      EndDo
C
      Else
C
      Q1(IJ,K)=Zero
      Q2(IJ,K)=Zero
      Write(6,'("Set to 0 in N^-1 ",4I3,3X,2E14.6)') 
     $ I,J,NSymNO(I),NSymNO(J),Occ(I),Occ(J)
C
      EndIf
      EndDo
C
      Call MultpM(HlpAB,Q2,Q1,NDimX)
C
C     else of IAPSG.Eq.0
      Else
C
      Do I=1,NDimX
      Do J=1,NDimX
c      ABPLUS(I,J)=ABPLUS(I,J)-Two*CMAT(I,J)
      EndDo
      EndDo


c herer!!!
c      Do IJ=1,NDimX
c      I=IndN(1,IJ)
c      J=IndN(2,IJ)
c      Do KL=1,NDimX
c      K=IndN(1,KL)
c      L=IndN(2,KL)
c      ABPLUS(IJ,KL)=ABPLUS(IJ,KL)*(CICoef(I)+CICoef(J))*
c     $ (CICoef(K)+CICoef(L))
c      ABMIN(IJ,KL)=ABPLUS(IJ,KL)
c      EndDo
c      EndDo


      Call MultpM(HlpAB,ABPLUS,ABMIN,NDimX)
C
C     endif of IAPSG.Eq.0
      EndIf
C
      Write(6,'(/,''Excitation Energies in [au] and [eV]'')')
C
      NISum=0
      Do NSym=1,MxSym
C
      Write(6,'(2/," ********** SYMMETRY CLASS NO ",I2," **********")')
     $ NSym
      Write(6,'(" Excit No        [au]        [eV] ")')
C
C     Find excitation energies for a given symmetry
C
      Index=0
      JS=0
C
      Do J=1,NDimX
      Ind1=IndN(1,J)
      Ind2=IndN(2,J)
      JSym=MultpC(NSymNO(Ind1),NSymNO(Ind2))
C
      If(JSym.Eq.NSym) Then
      JS=JS+1
      IndS(JS)=J
      EndIf
C
      Do I=1,NDimX
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      ISym=MultpC(NSymNO(Ind1),NSymNO(Ind2))
C
      If(ISym.Eq.NSym.And.JSym.Eq.NSym) Then
C
      Index=Index+1
      XR(Index)=HlpAB(I,J)
C
      EndIf
C
      EndDo
      EndDo
C
      NI=Sqrt(Float(Index))
      NISum=NISum+NI
C
      If(NI.Eq.0) GoTo 999
C
      Call DGEEV('N','V',NI,XR,NI,Eig,EigI,
     $           Q2,NI,Q1,NI,Work,4*NI,INFO)
C
      If(INFO.Ne.0) Stop 'INFO FROM DGEEV DIFFERENT FROM 0'
C
      Call SortEig(1,Eig,EigI,Q1,NI)
C
      Do I=1,Min(10,NI)
C
      If(Eig(I).Lt.Zero.And.Abs(Eig(I)).Lt.SMALL) Eig(I)=Zero 
      Im=0
      If(EigI(I).Ne.Zero) Im=1
C
      Write(6,'("SymERPA = ",I1,I5,5X,2F12.5,I2)')NSym,I,Sqrt(Eig(I)),
     $ Sqrt(Eig(I))*27.211,Im
C
      EndDo
      Write(6,*)
C
C     Find the composition of the XR,XI vectors
C
C     COPY EIGENVECTORS TO XR FOR CONVENIENCE 
C
      Call CpyM(XR,Q1,NDimX)
C
      Write(6,'(" Ind1 Ind2 Sym1 Sym2    XR")')
C
      Do I=1,Min(NI,10)
C
      Write(6,'("__________________________________________")')
      Write(6,'(X,"Excit No",I3)')I
C
      Do J=1,NI
C
      Ind1=IndN(1,IndS(J))
      Ind2=IndN(2,IndS(J))
C
      XRR=Abs(XR((I-1)*NI+J))
C
c herer!!!
      If(XRR.Gt.0.01) Then
      Write(*,'(4I3,F12.4)') Ind1,Ind2,NSymNO(Ind1),NSymNO(Ind2),
     $ XRR

      EndIf
C
      EndDo
      EndDo
C
C     End of NSym loop 
C
  999 Continue
      EndDo
C
      If(NISum.Ne.NDimX) Write(*,*) 
     $ 'Sum of Dims of symmetrized matrices different from NDimX'
C
      Return
      End

*Deck SERPAExcit
      Subroutine SERPAExcit(ABPLUS,ABMIN,DMAT,DMATK,EMAT,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ DMAT(NDimX,NDimN),EMAT(NDimN,NDimN),
     $ DMATK(NDimX,NDimN),
     $ Q1(2*(NDimX+NDimN),2*(NDimX+NDimN)),
     $ Q2(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecL(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Occ(NBasis),
     $ Eig(2*(NDimX+NDimN)),EigI(2*(NDimX+NDimN)),
     $ Work(4*2*(NDimX+NDimN)),
     $ IndN(2,NDimX),NSymI(2*(NDimX+NDimN)),
     $ NSymNO(NBasis),MultpC(15,15),IndS(2*NDimX+NDimN)
C
      Do I=1,2*NDimX+NDimN
C
      If(I.Le.NDimX) Then
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind2))
      EndIf
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then 
      Ind1=IndN(1,I-NDimX)
      Ind2=IndN(2,I-NDimX)
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind2))
      EndIf
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
      Ind1=I-2*NDimX
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind1))
      EndIf
C
      If(I.Gt.2*NDimX+NDimN) Then
      Ind1=I-2*NDimX-NDimN
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind1)) 
      EndIf
C
      Do J=1,2*NDimX+NDimN
      Q1(I,J)=Zero
C
C     1ST SUPER ROW
C
      If(I.Le.NDimX) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) Q1(I,J)=ABPLUS(I,J-NDimX)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=DMAT(I,J-2*NDimX)
C
      EndIf
C
C     2ND SUPER ROW
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then
C
      If(J.Le.NDimX) Q1(I,J)=ABMIN(I-NDimX,J)
       If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=DMATK(I-NDimX,J-2*NDimX)
C
      EndIf
C
C     3RD SUPER ROW 
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) 
     $ Q1(I,J)=Two*DMAT(J-NDimX,I-2*NDimX)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=EMAT(I-2*NDimX,J-2*NDimX)
C
      EndIf
C
      EndDo
      EndDo
C
      Write(7,'(/,''Excitation Energies in [au] and [eV]'')')
C
      NISum=0
      Do NSym=1,MxSym
C
      Write(6,'(2/," ********** SYMMETRY CLASS NO ",I2," **********")')
     $ NSym
      Write(6,'(8X," Excit No      [au]        [eV]   Im")')
C
C     Find excitation energies for a given symmetry
C
      Index=0
      JS=0
C
      Do J=1,2*NDimX+NDimN 
      JSym=NSymI(J)
      If(JSym.Eq.NSym) Then
      JS=JS+1
      IndS(JS)=J
      EndIf
C
      Do I=1,2*NDimX+NDimN
C
      ISym=NSymI(I)
C
      If(ISym.Eq.NSym.And.JSym.Eq.NSym) Then
C
      Index=Index+1
      Q2(Index)=Q1(I,J)
C
      EndIf
C
      EndDo
      EndDo
C
      NI=Sqrt(Float(Index))
      NISum=NISum+NI
C
      If(NI.Eq.0) GoTo 999
C
      Call DGEEV('N','V',NI,Q2,NI,Eig,EigI,
     $           EigVecL,NI,EigVecR,NI,Work,4*NI,INFO)
C
      If(INFO.Ne.0) Stop 'INFO FROM DGEEV DIFFERENT FROM 0'
C
      Call SortEig(1,Eig,EigI,EigVecR,NI)
C
      Do I=1,Min(20,NI)
C
      Im=0

      If(EigI(I).Eq.Zero) Then
      Write(6,'("SymSERPA = ",I1,I5,5X,2F12.5,I2)')NSym,I,Eig(I),
     $ Eig(I)*27.211,Im
      Else
      Write(6,'("SymSERPA = ",I1,I5,5X,3F12.5)')NSym,I,Eig(I),
     $ Eig(I)*27.211,EigI(I)*27.211
      EndIf
C
      EndDo
      Write(6,*)
C
C     Find the composition of the EigenVectors
C
      Write(6,'(" Ind1 Ind2 Sym1 Sym2    EigVec")')
C
      ICount=0
      Do I=1,NI
C
C     PRINT THE COMPOSITION OF EXCITATIONS CRRESPONDING TO THE POSITIVE VALUES ONLY
C
      If(Eig(I)*27.211.Gt.0.1D0.And.ICount.Lt.8) Then
C
      ICount=ICount+1
C
      Write(6,'("__________________________________________")')
      Write(6,'(X,"Excit No",I3)')I
C
      Do J=1,NI
C
      If(IndS(J).Le.NDimx) Then
      Ind1=IndN(1,IndS(J))
      Ind2=IndN(2,IndS(J))
      ElseIf(IndS(J).Le.2*NDimx) Then
      Ind1=IndN(1,IndS(J)-NDimX)
      Ind2=IndN(2,IndS(J)-NDimX)
      Else
      Ind1=IndS(J)-2*NDimX
      Ind2=Ind1
      EndIf
C
      XRR=EigVecR((I-1)*NI+J)
C
      If(Abs(XRR).Gt.0.1) Then
      Write(*,'(4I3,F12.4)') Ind1,Ind2,NSymNO(Ind1),NSymNO(Ind2),
     $ XRR
      EndIf
C
      EndDo
C
      EndIf
      EndDo 
C
C     End of NSym loop 
C
  999 Continue
      EndDo
C
      If(NISum.Ne.2*NDimX+NDimN) Write(*,*) 
     $ 'Sum of Dims of symmetrized matrices different from NDimX'
C
      Return
      End


*Deck PINOExcit
      Subroutine PINOExcit(APLUS,AMIN,DPLUS,DMIN,EPLUS,EMIN,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN)
C
C     EXCITATION ENERGIES FROM PINO WITH THE APSG FUNCTIONAL
C     NOTICE THE SIMILARITY WITH SERPA
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ APLUS(NDimX,NDimX),AMIN(NDimX,NDimX),
     $ DPLUS(NDimX,NDimN),EPLUS(NDimN,NDimN),
     $ DMIN(NDimX,NDimN),EMIN(NDimN,NDimN),
     $ Q1(2*(NDimX+NDimN),2*(NDimX+NDimN)),
     $ Q2(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecL(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Occ(NBasis),
     $ Eig(2*(NDimX+NDimN)),EigI(2*(NDimX+NDimN)),
     $ Work(4*2*(NDimX+NDimN)),
     $ IndN(2,NDimX),NSymI(2*(NDimX+NDimN)),
     $ NSymNO(NBasis),MultpC(15,15),IndS(2*(NDimX+NDimN))
C
      Do I=1,2*(NDimX+NDimN)
C
      If(I.Le.NDimX) Then
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind2))
      EndIf
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then
      Ind1=IndN(1,I-NDimX)
      Ind2=IndN(2,I-NDimX)
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind2))
      EndIf
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
      Ind1=I-2*NDimX
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind1))
      EndIf
C
      If(I.Gt.2*NDimX+NDimN) Then
      Ind1=I-2*NDimX-NDimN
      NSymI(I)=MultpC(NSymNO(Ind1),NSymNO(Ind1))
      EndIf
C
      Do J=1,2*(NDimX+NDimN)
      Q1(I,J)=Zero
C
C     1ST SUPER ROW
C
      If(I.Le.NDimX) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) Q1(I,J)=APLUS(I,J-NDimX)
      If(J.Gt.2*NDimX+NDimN) Q1(I,J)=DPLUS(I,J-2*NDimX-NDimN)
C
      EndIf
C
C     2ND SUPER ROW
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then
C
      If(J.Le.NDimX) Q1(I,J)=AMIN(I-NDimX,J)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=DMIN(I-NDimX,J-2*NDimX)
C
      EndIf
C
C     3RD SUPER ROW 
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX)
     $ Q1(I,J)=Two*DPLUS(J-NDimX,I-2*NDimX)
      If(J.Gt.2*NDimX+NDimN) 
     $ Q1(I,J)=EPLUS(I-2*NDimX,J-2*NDimX-NDimN)
C
      EndIf
C
C     4TH SUPER ROW
C
      If(I.Gt.2*NDimX+NDimN) Then
C
      If(J.Le.NDimX) Q1(I,J)=Two*DMIN(J,I-2*NDimX-NDimN)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=EMIN(I-2*NDimX-NDimN,J-2*NDimX)
C
      EndIf
C
      EndDo
      EndDo
C
      Write(7,'(/,''Excitation Energies in [au] and [eV]'')')
C
      NISum=0
      Do NSym=1,MxSym
C
      Write(6,'(2/," ********** SYMMETRY CLASS NO ",I2," **********")')
     $ NSym
      Write(6,'(8X," Excit No      [au]        [eV]   Im")')
C
C     Find excitation energies for a given symmetry
C
      Index=0
      JS=0
C
      Do J=1,2*(NDimX+NDimN)
      JSym=NSymI(J)
      If(JSym.Eq.NSym) Then
      JS=JS+1
      IndS(JS)=J
      EndIf
C
      Do I=1,2*(NDimX+NDimN)
C
      ISym=NSymI(I)
C
      If(ISym.Eq.NSym.And.JSym.Eq.NSym) Then
C
      Index=Index+1
      Q2(Index)=Q1(I,J)
C
      EndIf
C
      EndDo
      EndDo
C
      NI=Sqrt(Float(Index))
      NISum=NISum+NI
C
      If(NI.Eq.0) GoTo 999
C
      Call DGEEV('N','V',NI,Q2,NI,Eig,EigI,
     $           EigVecL,NI,EigVecR,NI,Work,4*NI,INFO)
C
      If(INFO.Ne.0) Stop 'INFO FROM DGEEV DIFFERENT FROM 0'
C
      Call SortEig(1,Eig,EigI,EigVecR,NI)
C
      Do I=1,Min(20,NI)
C
      Im=0

      If(EigI(I).Eq.Zero) Then
      Write(6,'("SymPINO = ",I1,I5,5X,2F12.5,I2)')NSym,I,Eig(I),
     $ Eig(I)*27.211,Im
      Else
      Write(6,'("SymPINO = ",I1,I5,5X,3F12.5)')NSym,I,Eig(I),
     $ Eig(I)*27.211,EigI(I)*27.211
      EndIf
C
      EndDo
      Write(6,*)
C
C     Find the composition of the EigenVectors
C
      Write(6,'(" Ind1 Ind2 Sym1 Sym2    EigVec")')
C
      ICount=0
      Do I=1,NI
C
C     PRINT THE COMPOSITION OF EXCITATIONS CRRESPONDING TO THE POSITIVE VALUES ONLY
C
      If(Eig(I)*27.211.Gt.0.1D0.And.ICount.Lt.8) Then
C
      ICount=ICount+1
C
      Write(6,'("__________________________________________")')
      Write(6,'(X,"Excit No",I3)')I
C
      Do J=1,NI
C
      If(IndS(J).Le.NDimx) Then
      Ind1=IndN(1,IndS(J))
      Ind2=IndN(2,IndS(J))
      ElseIf(IndS(J).Le.2*NDimX) Then
      Ind1=IndN(1,IndS(J)-NDimX)
      Ind2=IndN(2,IndS(J)-NDimX)
      ElseIf (IndS(J).Le.2*NDimX+NDimN) Then
      Ind1=IndS(J)-2*NDimX
      Ind2=Ind1
      Else
      Ind1=IndS(J)-2*NDimX-NDimN
      Ind2=Ind1
      EndIf
C
      XRR=EigVecR((I-1)*NI+J)
C
      If(Abs(XRR).Gt.0.1) Then
      Write(*,'(4I3,F12.4)') Ind1,Ind2,NSymNO(Ind1),NSymNO(Ind2),
     $ XRR
      EndIf
C
      EndDo
C
      EndIf
      EndDo
C
C     End of NSym loop 
C
  999 Continue
      EndDo
C
      If(NISum.Ne.2*(NDimX+NDimN)) Write(*,*)
     $ 'Sum of Dims of symmetrized matrices different from NDimX'
C
      Return
      End

*Deck LinAB_NEST
      Subroutine LinAB_NEST(ABPLUS,ABMIN,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE A+B AND A-B MATRICES IN THE LINEARIZED ERPA
C
C     NESTED COMMUTATOR 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),   
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension EPS(NInte1)
C
C     ONE-ELECTRON MATRIX IN A NO REPRESENTATION
C   
      IJ=0 
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      EPS(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      EPS(IJ)=EPS(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      Do K=1,NBasis
C
      EPS(IJ)=EPS(IJ)+Occ(K)*
     $ (Two*TwoMO(NAddr3(I,J,K,K))-TwoMO(NAddr3(I,K,J,K)))
C
      EndDo
C
      EndDo
      EndDo
C
C     COMPUTE ABPLUS AND ABMIN
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
      TwoEl=TwoMO(NAddr3(IS,IR,IP,IQ))
      BRSPQ=-(Occ(IR)-Occ(IS))*(Occ(IP)-Occ(IQ))*
     $ (Two*TwoEl-TwoMO(NAddr3(IS,IP,IQ,IR)))
      ARSPQ=(Occ(IR)-Occ(IS))*(Occ(IP)-Occ(IQ))*
     $ (Two*TwoEl-TwoMO(NAddr3(IS,IQ,IP,IR)))
C
      If(IS.Eq.IP) BRSPQ=BRSPQ+
     $ (Occ(IP)-Occ(IQ))*EPS((Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR))
      If(IR.Eq.IQ) BRSPQ=BRSPQ-
     $ (Occ(IP)-Occ(IQ))*EPS((Max(IP,IS)*(Max(IP,IS)-1))/2+Min(IP,IS))
C
      If(IR.Eq.IP) ARSPQ=ARSPQ+
     $ (Occ(IP)-Occ(IQ))*EPS((Max(IQ,IS)*(Max(IQ,IS)-1))/2+Min(IQ,IS))
      If(IS.Eq.IQ) ARSPQ=ARSPQ-
     $ (Occ(IP)-Occ(IQ))*EPS((Max(IP,IR)*(Max(IP,IR)-1))/2+Min(IP,IR))
C
C     CHANGE THE SIGN OF ARSPQ TO MAKE IT CONSISTENT WITH THE LINEAR RESPONSE APPROACH
C     (IT IS IMPORTANT TO ADD THE SR-KERNEL TO B+A_TILDE, A_TILDE=-A
C
      ARSPQ=-ARSPQ 
      ABPLUS(IRS,IPQ)=ARSPQ+BRSPQ
      ABMIN(IRS,IPQ)=ARSPQ-BRSPQ
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck ExactABC_NEST
      Subroutine ExactABC_NEST(AMAT,BMAT,CMAT,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C     
C     COMPUTE THE A+B, A-B, AND C MATRICES IN ERPA WITH THE EXACT 2-el 2-RDM
C    
C     NESTED COMMUTATOR 
C     
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C     
      Include 'commons.inc'
C     
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),CMAT(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),Aux(NInte1),ANDMAT(NDim,NBasis),
     $ ADNMAT(NBasis,NDim),ADDMAT(NBasis,NBasis),Aux1(NBasis,NDim),
     $ HNO(NInte1),Aux2(NBasis,NBasis)
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
      Do I=1,NBasis
      If(Occ(I).Lt.Half) Then
      C(I)=-SQRT(Occ(I))
      Else
      C(I)=SQRT(Occ(I))
      EndIf
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      Aux(IJ)=Zero
C
      Do IT=1,NBasis
      Aux(IJ)=Aux(IJ)+C(IT)*TwoMO(NAddr3(IT,I,IT,J))
      EndDo     
C
      EndDo
      EndDo
C
C     COMPUTE THE MATRICES
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
C     ONE-ELECTRON PART
C
      XMAT=Zero
C
      If(IQ.Eq.IR) XMAT=XMAT+(Occ(IR)-Occ(IS))*HNO(IPS)
      If(IS.Eq.IP) XMAT=XMAT-(Occ(IR)-Occ(IS))*HNO(IQR)
C
C     XC PART
C
      XMAT=XMAT+(C(IQ)*C(IR)+C(IP)*C(IS))
     $ *( TwoMO(NAddr3(IP,IS,IQ,IR))+TwoMO(NAddr3(IP,IR,IQ,IS)) )
C
      If(IS.Eq.IP) XMAT=XMAT-C(IR)*Aux(IQR) 
      If(IR.Eq.IQ) XMAT=XMAT-C(IS)*Aux(IPS)
      If(IR.Eq.IP) XMAT=XMAT-C(IP)*Aux(IQS)
      If(IS.Eq.IQ) XMAT=XMAT-C(IQ)*Aux(IPR) 
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
      Do I=1,NDim
      Do J=1,NDim
      A=AMAT(I,J)
      B=BMAT(I,J)
      AMAT(I,J)=A+B
      BMAT(I,J)=A-B
      EndDo
      EndDo
C
C     COMPUTE THE CMAT MATRIX DEFINED AS 
C     C = ANDMAT (ADDMAT)^-1 ADNMAT
C
C     INVERT THE ADDMAT MATRIX
C
      Call CpyM(Aux2,ADDMAT,NBasis)
      Call Diag8(Aux2,NBasis,NBasis,C,Aux)
C
      Do I=1,NBasis
      If(C(I).Gt.1.E-10) Then
      C(I)=One/C(I)
      Else
      C(I)=Zero
      EndIf
      EndDo
C
      Do I=1,NBasis
      Do J=1,I
      ADDMAT(I,J)=Zero
      Do K=1,NBasis
      ADDMAT(I,J)=ADDMAT(I,J)+Aux2(K,I)*C(K)*Aux2(K,J)
      EndDo
      ADDMAT(J,I)=ADDMAT(I,J)
      EndDo
      EndDo
C
      Call MultpMN(Aux1,ADDMAT,ADNMAT,NBasis,NBasis,NBasis,NDim)
      Call MultpMN(CMAT,ANDMAT,Aux1,NDim,NBasis,NBasis,NDim)        
C
      Return
      End

*Deck PrimAB_NEST
      Subroutine PrimAB_NEST(AMAT,BMAT,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2,IFlag)
C     
C     COMPUTE THE A+B AND A-B MATRICES IN ERPA WITH PRIMITIVE FUNCTIONALS
C
C     IFlag = 1 - HF-like XC PART OF 2-RDM:
C                 Gamma^XC_pqrs = - G_pq del_ps del_qr
C           = 0 - Lowdin-Shull-like XC PART of 2-RDM:
C                 Gamma^XC_pqrs = - G_pr del_pq del_rs
C
C     NESTED COMMUTATOR 
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C     
      Include 'commons.inc'
C     
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1)
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
C     COMPUTE THE MATRIX B
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C
      IPQ=0
      Do IP=1,NBasis
C
      ISP=(Max(IS,IP)*(Max(IS,IP)-1))/2+Min(IS,IP)
C
      Do IQ=1,NBasis
      IQR=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)

      If(IP.Gt.IQ) IPQ=IPQ+1
C
      If(IQ.Gt.IP) IQP=(IQ**2-3*IQ+2)/2+IP
C
      If(IP.Ne.IQ) Then
C
C     ONE-ELECTRON AND COULOMB PARTS      
C
      XMAT=Zero
C
      If(IQ.Eq.IR) XMAT=XMAT+(Occ(IR)-Occ(IS))*HNO(ISP)
      If(IS.Eq.IP) XMAT=XMAT-(Occ(IR)-Occ(IS))*HNO(IQR) 
C
      XMAT=XMAT-Half*(Occ(IR)*Occ(IP)+Occ(IS)*Occ(IQ))*
     $  ( Two*TwoMO(NAddr3(IQ,IP,IS,IR))-TwoMO(NAddr3(IQ,IR,IS,IP)) )
C
      If(IQ.Eq.IR) Then
      Sum=Zero
      Do IT=1,NBasis
      Sum=Sum+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IS,IP,IT,IT))-TwoMO(NAddr3(IS,IT,IP,IT)))
      EndDo
      XMAT=XMAT+(Occ(IQ)-Half*Occ(IS))*Sum
      EndIf
C
      If(IS.Eq.IP) Then
      Sum=Zero
      Do IT=1,NBasis
      Sum=Sum+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IQ,IR,IT,IT))-TwoMO(NAddr3(IR,IT,IQ,IT)))
      EndDo
      XMAT=XMAT+(Occ(IP)-Half*Occ(IR))*Sum
      EndIf
C
C     XC PART
C
      If(IFlag.Eq.1) Then
C
      XMAT=XMAT+(-GOCC(Occ(IQ),Occ(IR),0,IQ,IR)
     $           -GOCC(Occ(IP),Occ(IS),0,IP,IS)
     $      +Half*GOCC(Occ(IP),Occ(IR),0,IP,IR)
     $      +Half*GOCC(Occ(IQ),Occ(IS),0,IQ,IS) 
     $  )*( Two*TwoMO(NAddr3(IQ,IP,IS,IR))-TwoMO(NAddr3(IQ,IR,IS,IP)) )
C
      If(IQ.Eq.IR) Then
      Do IT=1,NBasis
      XMAT=XMAT+Half*GOCC(Occ(IS),Occ(IT),0,IS,IT)*
     $ (Two*TwoMO(NAddr3(IS,IP,IT,IT))-TwoMO(NAddr3(IS,IT,IP,IT)))
      EndDo
      EndIf
C
      If(IS.Eq.IP) Then
      Do IT=1,NBasis
      XMAT=XMAT+Half*GOCC(Occ(IR),Occ(IT),0,IR,IT)*
     $ (Two*TwoMO(NAddr3(IQ,IR,IT,IT))-TwoMO(NAddr3(IR,IT,IQ,IT)))
      EndDo
      EndIf
C
      Else
C
      XMAT=XMAT-
     $ ( GOCC(Occ(IQ),Occ(IR),0,IQ,IR)+GOCC(Occ(IP),Occ(IS),0,IP,IS) )*
     $ (TwoMO(NAddr3(IP,IR,IQ,IS))-TwoMO(NAddr3(IP,IS,IQ,IR)) )
C
      EndIf
C
      EndIf
C
      If (IP.Gt.IQ) BMAT(IRS,IPQ)=XMAT
      If (IQ.Gt.IP) AMAT(IRS,IQP)=XMAT
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do I=1,NDim
      Do J=1,NDim
      A=AMAT(I,J)
      B=BMAT(I,J)
      AMAT(I,J)=A+B
      BMAT(I,J)=A-B
      EndDo
      EndDo
C
      Return
      End

*Deck APSG_NEST
      Subroutine APSG_NEST(AMAT,BMAT,CMAT,EMAT,EMATM,
     $ DMAT,DMATM,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2,NGem,ISERPA)
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

*Deck PINOVEC
      Subroutine PINOVEC(EigVecR,Eig,
     $ APLUS,AMIN,DPLUS,DMIN,EPLUS,EMIN,
     $ Occ,NBasis,NDimX,NDimN)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Integer :: DimV1,Max_NDEG
      Integer :: Space1(3,2*(NDimX+NDimN))
C
      Include 'commons.inc'
C
      Dimension
     $ APLUS(NDimX,NDimX),AMIN(NDimX,NDimX),
     $ DPLUS(NDimX,NDimN),EPLUS(NDimN,NDimN),
     $ DMIN(NDimX,NDimN),EMIN(NDimN,NDimN),
     $ Q1(2*(NDimX+NDimN),2*(NDimX+NDimN)),
     $ Q2(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecL(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Occ(NBasis),
     $ Eig(2*(NDimX+NDimN)),EigI(2*(NDimX+NDimN)),
     $ Work(4*2*(NDimX+NDimN))

C
      Do I=1,NDimX
      Do J=I+1,NDimX
      APLUS(I,J)=Half*(APLUS(I,J)+APLUS(J,I))
      APLUS(J,I)=APLUS(I,J)
      AMIN(I,J)=Half*(AMIN(I,J)+AMIN(J,I))
      AMIN(J,I)=AMIN(I,J)
      EndDo
      EndDo
C
      Do I=1,NDimX
      Do J=1,NDimN
      If(Occ(J).Eq.One) Then
      DPLUS(I,J) = Zero
      DMIN(I,J) = Zero
      EndIf
      EndDo
      EndDo
C
      Do I=1,NDimN
      Do J=1,NDimN
      If(Occ(J).Eq.One) Then
      EPLUS(I,J) = Zero
      EMIN(I,J) = Zero
      EndIf
      EndDo
      EndDo
C
CC     ADDITIONAL TEST
C      Write(6,*) 'D(E)MIN = D(E)PLUS'
C      DMIN = DPLUS
C      EMIN = EPLUS
CC
C    
C      Print*, 'DEPLUS,DEMIN' 
C      Print*, norm2(DPLUS),norm2(DMIN)
C      Print*, norm2(EPLUS),norm2(EMIN)
C
      Do I=1,2*(NDimX+NDimN)
      Do J=1,2*(NDimX+NDimN)
      Q1(I,J)=Zero
C
C     1ST SUPER ROW
C
      If(I.Le.NDimX) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) Q1(I,J)=APLUS(I,J-NDimX)
      If(J.Gt.2*NDimX+NDimN) Q1(I,J)=DPLUS(I,J-2*NDimX-NDimN)
C
      EndIf
C
C     2ND SUPER ROW
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then
C
      If(J.Le.NDimX) Q1(I,J)=AMIN(I-NDimX,J)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=DMIN(I-NDimX,J-2*NDimX)
C
      EndIf
C
C     3RD SUPER ROW 
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX)
     $ Q1(I,J)=Two*DPLUS(J-NDimX,I-2*NDimX)
      If(J.Gt.2*NDimX+NDimN)
     $ Q1(I,J)=EPLUS(I-2*NDimX,J-2*NDimX-NDimN)
C
      EndIf
C
C     4TH SUPER ROW
C
      If(I.Gt.2*NDimX+NDimN) Then
C
      If(J.Le.NDimX) Q1(I,J)=Two*DMIN(J,I-2*NDimX-NDimN)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=EMIN(I-2*NDimX-NDimN,J-2*NDimX)
C
      EndIf
C
      EndDo
      EndDo
C
      NI=2*(NDimX+NDimN)
      Index=0
      Do J=1,NI
      Do I=1,NI
      Index=Index+1
      Q2(Index)=Q1(I,J)
      EndDo
      EndDo
C   
      Call DGEEV('N','V',NI,Q2,NI,Eig,EigI,
     $           EigVecL,NI,EigVecR,NI,Work,4*NI,INFO)
C
C     ORTHOGONALISE DEGENERATE VECTORS
      Write(6,*) 'ORTHOGONALIZE DEGENERATE VECTORS IN PINOVEC'
      Call SORT_PINOVECS(EigVecR,Eig,EigI,NI)
      Call CREATE_SPACE(Eig,Space1,NI,DimV1,Max_NDEG)
      Call ORTHO_DEGVEC(EigVecR,Space1,DimV1,NI,Max_NDEG)
C      Call ZERO_DEGVEC(EigVecR,Eig,Space1,DimV1,NDimEx,Max_NDEG)
C    
      NDeg = 0
      Do I=1,NI
      If(Space1(3,I).Gt.1) NDeg = NDeg + 1
      EndDo
      Write(6,*) 'NUMBER OF DEGENERATE VECS:', NDeg 
C      Write(6,*) 'DimV1                    :', DimV1
      Write(6,*) 'MAXIMUM DEGENERACY       :', Max_NDEG
      Write(6,*) ' '
C
C     NORMALIZE [Y,X,W,V] SO THAT 2 Y*X + V*W = 1
C     THE NORMALIZATION IS APPLIED TO EIGVECTORS CORRESPONDING TO POSITIVE OMEGA's
C
      Do NU=1,NI
C
      If(Abs(EigI(NU)).Gt.1.D-12) Then
C
      Write(6,'(X,"Complex PINO Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      EigI(NU)=Zero
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=Zero
      EndDo      
C
      ElseIf(Abs(Eig(NU)).Lt.1.D-2) Then
      Write(6,'(X,"Zero PINO Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=Zero
      EndDo      
C
      EndIf
      EndDo 
C
      Do NU=1,NI
C
      SumNU=Zero
      Do I=1,NDimX+NDimN
      If(I.Le.NDimX) Then
      SumNU=SumNU+Two*
     $ EigVecR((NU-1)*NI+I)*EigVecR((NU-1)*NI+NDimX+I)
      Else
      SumNU=SumNU+EigVecR((NU-1)*NI+NDimX+I)*
     $ EigVecR((NU-1)*NI+NDimX+NDimN+I)
      EndIf
      EndDo
C     
C      OLD
C      If(SumNU.Gt.Zero) Then      
C      SumNU=One/SQRT(SumNU)
C      If(Eig(NU).Lt.Zero) Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
C     $  NU,Eig(NU),SumNU
C      Eig(NU)=Abs(Eig(NU))
C      Else
C      SumNU=Zero
C      Eig(NU)=-Abs(Eig(NU))
C      EndIf
C
C     CHANGE 4
C
      If(SumNU.Gt.Zero) Then
      SumNU=One/SQRT(SumNU)
      If(Eig(NU).Lt.Zero) Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
     $  NU,Eig(NU),SumNU
      Else
      SumNU=Zero
      Eig(NU)=Zero
      EndIf
C
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=EigVecR((NU-1)*NI+I)*SumNU
      EndDo
C
c     enddo NU
      EndDo
C
      Return
      End

*Deck PINOVECRED
      Subroutine PINOVECRED(EigVecR,Eig,INegExcit,
     $ APLUS,AMIN,DPLUS,DMIN,EPLUS,EMIN,
     $ NBasis,NDimX,NDimN)
C
C     A REDUCED NONSYMMETRIC PINO PROBLEM IS SOLVED
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
C     OLD:
     $ Four=4.D0,Small=1.D-10)
C     $ Four=4.D0,Small=1.D-2)
C
      Integer :: DimV1,Max_NDEG
      Integer :: Space1(3,2*(NDimX+NDimN))
C
      Include 'commons.inc'
C
      Dimension
     $ APLUS(NDimX,NDimX),AMIN(NDimX,NDimX),
     $ DPLUS(NDimX,NDimN),EPLUS(NDimN,NDimN),
     $ DMIN(NDimX,NDimN),EMIN(NDimN,NDimN),
     $ Q1(NDimX+NDimN,NDimX+NDimN),
     $ Q2(NDimX+NDimN,NDimX+NDimN),
     $ EigVecR((NDimX+NDimN)*(NDimX+NDimN)),
     $ Eig(NDimX+NDimN),Work(5*(NDimX+NDimN)),
     $ HlpAB(NDimX+NDimN,NDimX+NDimN),
     $ EigI(NDimX+NDimN)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      APLUS(I,J)=Half*(APLUS(I,J)+APLUS(J,I))
      APLUS(J,I)=APLUS(I,J)
      AMIN(I,J)=Half*(AMIN(I,J)+AMIN(J,I))
      AMIN(J,I)=AMIN(I,J)
      EndDo
      EndDo
C
      Do I=1,NDimX+NDimN
      Do J=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
C
      If(J.Le.NDimX) Then
      Q1(I,J)=APLUS(I,J)
      Q2(I,J)=AMIN(I,J)
      Else
      Q1(I,J)=DPLUS(I,J-NDimX)
      Q2(I,J)=DMIN(I,J-NDimX)
      EndIf
C
      Else
C
      If(J.Le.NDimX) Then
      Q1(I,J)=Two*DPLUS(J,I-NDimX)
      Q2(I,J)=Two*DMIN(J,I-NDimX)
      Else
      Q1(I,J)=EPLUS(I-NDimX,J-NDimX)
      Q2(I,J)=EMIN(I-NDimX,J-NDimX)
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NI=NDimX+NDimN
C
      Call MultpM(HlpAB,Q1,Q2,NI)
C
      Call DGEEV('N','V',NI,HlpAB,NI,Eig,EigI,
     $           Q1,NI,EigVecR,NI,Work,5*I,INFO)
C
C     ORTHOGONALISE DEGENERATE VECTORS
      Write(6,*) 'ORTHOGONALIZE DEGENERATE VECTORS IN PINOVEC'
      Call SORT_PINOVECS(EigVecR,Eig,EigI,NI)
      Call CREATE_SPACE(Eig,Space1,NI,DimV1,Max_NDEG)
      Call ORTHO_DEGVEC(EigVecR,Space1,DimV1,NI,Max_NDEG)
C      Call ZERO_DEGVEC(EigVecR,Eig,Space1,DimV1,NDimEx,Max_NDEG)
C    
      NDeg = 0
      Do I=1,NI
      If(Space1(3,I).Gt.1) NDeg = NDeg + 1
      EndDo
      Write(6,*) 'NUMBER OF DEGENERATE VECS:', NDeg 
C      Write(6,*) 'DimV1                    :', DimV1
      Write(6,*) 'MAXIMUM DEGENERACY       :', Max_NDEG
      Write(6,*) ' '
C
      Do NU=1,NI
C
      If(Abs(EigI(NU)).Gt.1.D-12) Then
C
      Write(6,'(X,"Complex PINO Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=Zero
      EndDo
C
      EndIf
      EndDo
C
C     IMPOSE THE NORMALIZATION 2 Y*X + V*W = 1 ON THE EIGENVECTORS
C
      Write(6,'(X,"Threshold for small PINO Eigenvalue ",E12.4)')
     $ Small
C
      Do NU=1,NI
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C    
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NI
C
      If(I.Le.NDimX) Then
C
      Do J=1,NI
      SumNU=SumNU+Two/Eig(NU)*Q2(I,J)*
     $ EigVecR((NU-1)*NI+I)*EigVecR((NU-1)*NI+J)
      EndDo
C
      Else
C
      Do J=1,NI
      SumNU=SumNU+One/Eig(NU)*Q2(I,J)*
     $ EigVecR((NU-1)*NI+I)*EigVecR((NU-1)*NI+J)
      EndDo
C
      EndIf
C
      EndDo
C
C     OLD
C      If(SumNU.Gt.Zero) Then
C      SumNU=One/Sqrt(SumNU)
C      Else
C      WritE(*,*)'Neg Norm in PINOVECSYMM',SumNU
C      SumNU=Zero
C      EndIf
C
CC     CHANGE 5
      If(SumNU.Gt.Zero) Then
      SumNU=One/Sqrt(SumNU)
      Else
      Eig(NU)=-Eig(NU)
      Write(*,*) 'Neg Norm in PINOVECSYMM', NU,Eig(NU),SumNU
      SumNU=One/Sqrt(Abs(SumNU))
      EndIf
C
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=EigVecR((NU-1)*NI+I)*SumNU
      EndDo
C
c     If(Eig(NU).Gt.Small) Then
      Else
C
      Write(6,'(X,"Small or Negative PINO Eigenvalue ",I4,E12.4)')
     $ NU,Eig(NU)
      Eig(NU)=Zero
      Do I=1,NI
      EigVecR((NU-1)*NI+I)=Zero
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
      Return
      End

*Deck PINOVECREDXY
      Subroutine PINOVECREDXY(EigVecY,EigVecX,Eig,INegExcit,
     $ APLUS,AMIN,DPLUS,DMIN,EPLUS,EMIN,IndN,
     $ NBasis,NDimX,NDimN)
C
C     A REDUCED NONSYMMETRIC PINO PROBLEM IS SOLVED
C     UNTILDED! [Y,W]->EigVecY and [X,0]->EigVecX VECTORS ARE RETURNED
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
C     OLD:
     $ Four=4.D0,Small=1.D-10)
C     $ Four=4.D0,Small=1.D-2)
C
      Integer :: DimV1,Max_NDEG
      Integer :: Space1(3,2*(NDimX+NDimN))
C
      Include 'commons.inc'
C
      Dimension
     $ APLUS(NDimX,NDimX),AMIN(NDimX,NDimX),
     $ DPLUS(NDimX,NDimN),EPLUS(NDimN,NDimN),
     $ DMIN(NDimX,NDimN),EMIN(NDimN,NDimN),
     $ Q1(NDimX+NDimN,NDimX+NDimN),
     $ Q2(NDimX+NDimN,NDimX+NDimN),
     $ EigVecY((NDimX+NDimN)*(NDimX+NDimN)),
     $ EigVecX((NDimX+NDimN)*(NDimX+NDimN)),
     $ Eig(NDimX+NDimN),Work(5*(NDimX+NDimN)),
     $ HlpAB(NDimX+NDimN,NDimX+NDimN),
     $ EigI(NDimX+NDimN),
     $ IndN(2,NDimX)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      APLUS(I,J)=Half*(APLUS(I,J)+APLUS(J,I))
      APLUS(J,I)=APLUS(I,J)
      AMIN(I,J)=Half*(AMIN(I,J)+AMIN(J,I))
      AMIN(J,I)=AMIN(I,J)
      EndDo
      EndDo
C
      Do I=1,NDimX+NDimN
      Do J=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
C
      If(J.Le.NDimX) Then
      Q1(I,J)=APLUS(I,J)
      Q2(I,J)=AMIN(I,J)
      Else
      Q1(I,J)=DPLUS(I,J-NDimX)
      Q2(I,J)=DMIN(I,J-NDimX)
      EndIf
C
      Else
C
      If(J.Le.NDimX) Then
      Q1(I,J)=Two*DPLUS(J,I-NDimX)
      Q2(I,J)=Two*DMIN(J,I-NDimX)
      Else
      Q1(I,J)=EPLUS(I-NDimX,J-NDimX)
      Q2(I,J)=EMIN(I-NDimX,J-NDimX)
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NI=NDimX+NDimN
C
      Call MultpM(HlpAB,Q1,Q2,NI)
C
      Call DGEEV('N','V',NI,HlpAB,NI,Eig,EigI,
     $           Q1,NI,EigVecY,NI,Work,5*I,INFO)
C
C     ORTHOGONALISE DEGENERATE VECTORS
      Write(6,*) 'ORTHOGONALIZE DEGENERATE VECTORS IN PINOVECREDXY'
      Call SORT_PINOVECS(EigVecY,Eig,EigI,NI)
      Call CREATE_SPACE(Eig,Space1,NI,DimV1,Max_NDEG)
      Call ORTHO_DEGVEC(EigVecY,Space1,DimV1,NI,Max_NDEG)
C      Call ZERO_DEGVEC(EigVecY,Eig,Space1,DimV1,NDimEx,Max_NDEG)
C
      NDeg = 0
      Do I=1,NI
      If(Space1(3,I).Gt.1) NDeg = NDeg + 1
      EndDo
      Write(6,*) 'NUMBER OF DEGENERATE VECS:', NDeg
C      Write(6,*) 'DimV1                    :', DimV1
      Write(6,*) 'MAXIMUM DEGENERACY       :', Max_NDEG
      Write(6,*) ' '
C
      Do NU=1,NI
C
      If(Abs(EigI(NU)).Gt.1.D-12) Then
C
      Write(6,'(X,"Complex PINO Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      Do I=1,NI
      EigVecY((NU-1)*NI+I)=Zero
      EndDo
C
      EndIf
      EndDo
C
C     IMPOSE THE NORMALIZATION 2 Y*X + V*W = 1 ON THE EIGENVECTORS
C
      Write(6,'(X,"Threshold for small PINO Eigenvalue ",E12.4)')
     $ Small
C
      Do NU=1,NI
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NI
C
      EigVecX((NU-1)*NI+I)=Zero
      Do J=1,NI
      EigVecX((NU-1)*NI+I)=EigVecX((NU-1)*NI+I)+One/Eig(NU)*Q2(I,J)*
     $                                          EigVecY((NU-1)*NI+J)
      EndDo
C
      If(I.Le.NDimX) Then
      SumNU=SumNU+Two*EigVecY((NU-1)*NI+I)*EigVecX((NU-1)*NI+I)
      Else
      SumNU=SumNU+EigVecY((NU-1)*NI+I)*EigVecX((NU-1)*NI+I)
      EndIf
C
      EndDo
C
CC     CHANGE 5
      If(SumNU.Gt.Zero) Then
      SumNU=One/Sqrt(SumNU)
      Else
      Eig(NU)=-Eig(NU)
      Write(*,*) 'Neg Norm in PINOVECSYMM', NU,Eig(NU),SumNU
      SumNU=One/Sqrt(Abs(SumNU))
      EndIf
C
C     VECTORS ARE NORMALIZED TO A POSITIVE NUMEBR FOR POSITIVE AND NEGATIVE Omega
C
      Do I=1,NI
      EigVecY((NU-1)*NI+I)=EigVecY((NU-1)*NI+I)*SumNU
      If(Eig(NU).Gt.Zero) Then
      EigVecX((NU-1)*NI+I)=EigVecX((NU-1)*NI+I)*SumNU
      Else
      EigVecX((NU-1)*NI+I)=-EigVecX((NU-1)*NI+I)*SumNU 
      EndIf
      EndDo
C
c     If(Eig(NU).Gt.Small) Then
      Else
C
      Write(6,'(X,"Small or Negative PINO Eigenvalue ",I4,E12.4)')
     $ NU,Eig(NU)
      Eig(NU)=Zero
      Do I=1,NI
      EigVecY((NU-1)*NI+I)=Zero
      EigVecX((NU-1)*NI+I)=Zero
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
C     DOUBLE CHECKING OF THE NORMALIZATION
C
      Do NU=1,NI
      SumNU=Zero
      If(Eig(NU).Ne.Zero) Then 
      Do I=1,NI
      If(I.Le.NDimX) Then
      SumNU=SumNU+Two*EigVecY((NU-1)*NI+I)*EigVecX((NU-1)*NI+I)
      Else
      SumNU=SumNU+EigVecY((NU-1)*NI+I)*EigVecX((NU-1)*NI+I)
      EndIf
      EndDo
C
      If(SumNU.Le.Zero) Then
      Write(*,*)'NU Omega(NU)',NU,Eig(NU)
      Stop 'Error in PINOVECREDXY: wrong SumN'
      EndIf
C
      EndIf
      EndDo
C
C     FIND UNTILDED Y,X,W, SET V TO ZERO (IT SHOULD HAVE NO CONTRIBUTION TO ANYTHING
C     BECAUSE IT DOES NOT CORRESPOND TO 1-TRDM
C 
C     IN TERMS OF UNTILDED VECTORS WE HAVE:
C  
C     p>q,  1-TRDM_pq = -(n_p - n_q) X_pq
C     p>q,  1-TRDM_qp =  (n_p - n_q) Y_pq
C     p=q   1-TRDM_pp = W_pp 
C
C     EigVecX=[X,0] and EigVecY = [Y,W]
C
      Do K=1,NI
      Do I=1,NI
C
      If(I.Le.NDimX) Then
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
C
      X=EigVecX((K-1)*NI+I)/(CICoef(IP)+CICoef(IQ))
      Y=EigVecY((K-1)*NI+I)/(CICoef(IP)-CICoef(IQ))
C
      EigVecX((K-1)*NI+I)=Half*(X-Y)
      EigVecY((K-1)*NI+I)=Half*(X+Y)
C
      Else
C
      EigVecX((K-1)*NI+I)=Zero
      EigVecY((K-1)*NI+I)=CICoef(I-NDimX)*EigVecY((K-1)*NI+I)
C
      EndIf
C
      EndDo
      EndDo
C
      Return
      End

*Deck Get_APL
      Subroutine Get_APL(IR,IS,IP,IQ,APL,HNO,XMu,TwoMO,AuxH,
     $ NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C     
      Dimension HNO(NInte1),AuxH(NBasis,NInte1),TwoMO(NInte2)
C
C     FOR GIVEN INDICES IR,IS,IP,IQ BELONGING TO THE SAME GEMINAL
C     FIND A_PLUS
C
      ISP=(Max(IS,IP)*(Max(IS,IP)-1))/2+Min(IS,IP)
      IQR=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)
      ISQ=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      APL=TwoMO(NAddr3(IP,IS,IQ,IR))+TwoMO(NAddr3(IP,IR,IQ,IS))
C
      If(IQ.Eq.IR) APL=APL+HNO(ISP)+AuxH(IR,ISP)
      If(IS.Eq.IP) APL=APL+HNO(IQR)+AuxH(IR,IQR)
      If(IS.Eq.IQ) APL=APL+HNO(IPR)+AuxH(IR,IPR)
      If(IP.Eq.IR) APL=APL+HNO(ISQ)+AuxH(IR,ISQ)
      If(IQ.Eq.IR.And.IP.Eq.IS) APL=APL-XMu
      If(IQ.Eq.IS.And.IP.Eq.IR) APL=APL-XMu
C
      Return
      End

*Deck Get_C
      Subroutine Get_C(IP,IQ,IR,IS,CC,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C     
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C      
      Dimension C(NBasis),Occ(NBasis),HNO(NInte1),
     $ AuxH(NBasis,NInte1),AuxXC(NBasis,NInte1),TwoMO(NInte2)
C     
      CC=C(IP)*(TwoMO(NAddr3(IQ,IS,IP,IR))+TwoMO(NAddr3(IQ,IR,IP,IS)))+
     $ Occ(IQ)*(C(IR)-C(IS))*
     $ (TwoMO(NAddr3(IQ,IS,IP,IR))-TwoMO(NAddr3(IQ,IR,IP,IS)))
C
      ISQ=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IQR=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)
      If(IP.Eq.IR) CC=CC+(C(IR)-C(IS))*(HNO(ISQ)+AuxH(IR,ISQ))
     $-AuxXC(IR,ISQ)
      If(IP.Eq.IS) CC=CC-(C(IR)-C(IS))*(HNO(IQR)+AuxH(IR,IQR))
     $-AuxXC(IR,IQR)
C
      Return
      End

*Deck Get_D
      Subroutine Get_D(IP,IQ,IR,IS,DD,C,Occ,HNO,TwoMO,AuxH,AuxXC,
     $ NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C     
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C      
      Dimension C(NBasis),Occ(NBasis),HNO(NInte1),
     $ AuxH(NBasis,NInte1),AuxXC(NBasis,NInte1),TwoMO(NInte2)
C     
      DD=-C(IQ)*(TwoMO(NAddr3(IQ,IS,IP,IR))+TwoMO(NAddr3(IQ,IR,IP,IS)))+
     $ Occ(IP)*(C(IR)+C(IS))*
     $ (Four*TwoMO(NAddr3(IP,IQ,IR,IS))-TwoMO(NAddr3(IQ,IS,IP,IR))
     $ -TwoMO(NAddr3(IQ,IR,IP,IS)))
C
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
      ISP=(Max(IS,IP)*(Max(IS,IP)-1))/2+Min(IS,IP)
      If(IS.Eq.IQ) DD=DD-(C(IR)+C(IS))*(HNO(IPR)+AuxH(IR,IPR))
     $-AuxXC(IR,IPR)
      If(IQ.Eq.IR) DD=DD-(C(IR)+C(IS))*(HNO(ISP)+AuxH(IR,ISP))
     $-AuxXC(IR,ISP)
C
      Return
      End

*Deck SortEig
      Subroutine SortEig(IFlag,Eig,EigI,EigVec,N)
C
C     SORT Eig IN AN ASCENDING ORDER OF THE MODULUS OF EIG 
C     AND CHANGE THE ORDER OF THE IMAGINARY Eig AND EigVec 
C
C     IFlag = 1 - SORT ACC. TO MODULUS OF EIG; EIGENVECTORS STORED IN COLUMNS OF EigVec
C           = 0 - SORT ACC. TO EIG; EIGENVECTORS STORED IN ROWS OF EigVec     
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Eig(N),EigI(N),EigVec(N*N)
C
C     LOCAL ARRAYS
C
      Dimension EigVecOld(N,N),Save(N),Ind(N)
C
      Do I=1,N
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,N
C
      EMin=Abs(Eig(IStart))
      If(IFlag.Eq.0) EMin=Eig(IStart)
      IndMin=IStart
C
      Do J=IStart,N
      If(Abs(Eig(J)).Lt.EMin) Then
      EMin=Abs(Eig(J))
      IndMin=J
      EndIf
      EndDo
C
      If(IFlag.Eq.0) Then
      Do J=IStart,N
      If(Eig(J).Lt.EMin) Then
      EMin=Eig(J)
      IndMin=J
      EndIf
      EndDo
      EndIf
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
      Save(I)=EigI(I)
      Do J=1,N
      IJ=(J-1)*N+I
      EigVecOld(I,J)=EigVec(IJ)
      EndDo
      EndDo
C     
      Do I=1,N
C
      EigI(I)=Save(Ind(I))
C
      If(IFlag.Eq.1) Then
      Do J=1,N
      IJ=(J-1)*N+I
      EigVec(IJ)=EigVecOld(I,Ind(J))
      EndDo
      Else
      Do J=1,N
      IJ=(J-1)*N+I
      EigVec(IJ)=EigVecOld(Ind(I),J)
      EndDo
      EndIf
C
      EndDo
C
      Return
      End
C
C
*Deck EnePINO
      Subroutine EnePINO(ETot,ENuc,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NDimN)
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(SmallE=1.D-2,BigE=1.D20)
c      Parameter(SmallE=1.D-3,BigE=1.D2)
C
C     ONLY EXCITATIONS > SmallE AND < BigE ARE INCLUDED 
C     
      Include 'commons.inc'
C     
      Dimension EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Eig(2*(NDimX+NDimN)),URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),C(NBasis)
C
      NI=2*(NDimX+NDimN)
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
C     ONE_ELECTRON PART + HARTREE + EXCHANGE + np(1-np) <pp|pp>
C     COMPUTE THE APSG ENERGY
C
      ETot=Zero
      EAPSG=Zero
C
      IJ=0
      Do I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*HNO(II)
     $ -Half*Occ(I)*(One-Occ(I))*TwoNO(NAddr3(I,I,I,I))
C
      EAPSG=EAPSG+Two*Occ(I)*HNO(II)
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Occ(I)*Occ(J)*(TwoNO(NAddr3(I,I,J,J))
     $ -Half*TwoNO(NAddr3(I,J,I,J)))
C
      If(IGem(I).Eq.IGem(J)) Then
      EAPSG=EAPSG+FacIJ*C(I)*C(J)*TwoNO(NAddr3(I,J,I,J))
      Else
      EAPSG=EAPSG+FacIJ*Occ(I)*Occ(J)*(Two*TwoNO(NAddr3(I,I,J,J))
     $ -TwoNO(NAddr3(I,J,I,J)))
      EndIf
C
      EndDo
      EndDo
C
      ETot=Two*ETot
C
C     ADD CONTRIBUTIONS FROM THE Y,W EIGENVECTORS
C     COMPUTE INTERGEMINAL CORRELATION TERMS
C
      ECorr1=Zero
      ECorr1CT=Zero
      ECorr2CT=Zero 
      ECorrDispX=Zero
C
C     EIntra - only contributions from terms with all indices belonging to the same geminal
C     EAll   - contributions from all terms 
C
      EIntra=Zero
      EAll=Zero
C
      Do I=1,NDimX+NDimN
C
      If(I.Le.NDimX) Then
      IP=IndN(1,I)
      IR=IndN(2,I)
      Else
      IP=I-NDimX
      IR=IP
      EndIf
C
      Do J=1,NDimX+NDimN
C
      If(J.Le.NDimX) Then
      IQ=IndN(1,J)
      IS=IndN(2,J)
      Else
      IQ=J-NDimX
      IS=IQ
      EndIf
C
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      SumY=Zero
      Do K=1,NI
      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE)
     $ SumY=SumY+EigVecR((K-1)*NI+I)*EigVecR((K-1)*NI+J)
      EndDo
C
      Aux=Two*(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      AuxInterG=Zero
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).And.
     $ IGem(IP).Ne.IGem(IQ)) AuxInterG=Aux 
C
      Aux1CT=Zero
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Ne.IGem(IS).And.
     $ (IGem(IP).Eq.IGem(IQ).Or.IGem(IP).Eq.IGem(IS))) 
     $ Aux1CT=Aux1CT+Aux
      If(IGem(IP).Ne.IGem(IR).And.IGem(IQ).Eq.IGem(IS).And.
     $ (IGem(IP).Eq.IGem(IQ).Or.IGem(IR).Eq.IGem(IQ))) 
     $ Aux1CT=Aux1CT+Aux
C
      Aux2CT=Zero
      If(IGem(IP).Ne.IGem(IR).And.IGem(IQ).Ne.IGem(IS).And.
     $ IGem(IP).Eq.IGem(IQ).And.IGem(IR).Eq.IGem(IS)) Then
      Aux2CT=Aux2CT+Aux
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux2CT=Aux2CT
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
      EndIf
C
      AuxDispX=Zero
      If(IGem(IP).Ne.IGem(IR).And.IGem(IQ).Ne.IGem(IS).And.
     $ IGem(IP).Eq.IGem(IS).And.IGem(IR).Eq.IGem(IQ))
     $ AuxDispX=AuxDispX+Aux
C
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux=Aux
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
C
      EAll=EAll+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ)) EIntra=EIntra
     $ +Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      ETot=ETot+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      ECorr1=ECorr1+AuxInterG*TwoNO(NAddr3(IP,IR,IQ,IS))
      ECorr1CT=ECorr1CT+Aux1CT*TwoNO(NAddr3(IP,IR,IQ,IS))
      ECorr2CT=ECorr2CT+Aux2CT*TwoNO(NAddr3(IP,IR,IQ,IS))
      ECorrDispX=ECorrDispX+AuxDispX*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
      If(IP.Gt.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NI
      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE)
     $ SumY=SumY+EigVecR((K-1)*NI+I)*EigVecR((K-1)*NI+NDimX+J)
      EndDo
C
      EAll=EAll+Four*C(IQ)*(C(IP)+C(IR))*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IP).Eq.IGem(IQ)) EIntra=EIntra
     $ +Four*C(IQ)*(C(IP)+C(IR))*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      ETot=ETot+Four*C(IQ)*(C(IP)+C(IR))*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IP).Ne.IGem(IQ)) ECorr1=
     $ ECorr1+Four*C(IQ)*(C(IP)+C(IR))*SumY*TwoNO(NAddr3(IP,IR,IQ,IQ))
C
      EndIf
C
      If(IP.Eq.IR.And.IQ.Eq.IS) Then
C
      SumY=Zero
      Do K=1,NI
      If(Eig(K).Gt.1.D-3.And.Eig(K).Lt.BigE)
     $ SumY=SumY+EigVecR((K-1)*NI+NDimX+I)*EigVecR((K-1)*NI+NDimX+J)
      EndDo
C
      EAll=EAll+Two*C(IQ)*C(IP)*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))
C
      If(IGem(IP).Eq.IGem(IQ)) EIntra=EIntra
     $ +Two*C(IQ)*C(IP)*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))

      ETot=ETot+Two*C(IQ)*C(IP)*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))
C
      If(IGem(IP).Ne.IGem(IQ)) ECorr1=
     $ ECorr1+Two*C(IQ)*C(IP)*SumY*TwoNO(NAddr3(IP,IP,IQ,IQ))
C
      EndIf
C
      EndDo
      EndDo
C
       Write
     $(6,'(/,1X,''EPINO + ENuc     '', 45X,F15.8)')ETot+ENuc
C
       Write
     $ (6,'(1X,''EAPSG+ENuc, Corr, PINO-APSG '',4X,3F15.8)')EAPSG+ENuc,
     $ Half*(EAll-EIntra),EAPSG+ENuc+Half*(EAll-EIntra)
C
      Return
      End


