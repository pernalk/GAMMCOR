*Deck Exact2
      Subroutine Exact2(TwoMO,TwoMOLR,URe,Occ,XOne,
     $  OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,NDimKer)
C
C     A ROUTINE FOR EXCITATION ENERGY FROM LINEARIZED ERPA
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Eight=8.D0)
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
     $ AMAT(NDim*NDim),CMAT(NDim*NBasis),
     $ EMAT(NBasis*NBasis),
     $ AMATK(NDim*NDim),CMATK(NDim*NBasis),
     $ EMATK(NBasis*NBasis),
     $ IndX(NDim),IndN(2,NDim),
     $ XKer(NDimKer),
     $ NSymMO(NBasis),NSymNO(NBasis),
     $ MultpC(15,15)
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
      Write(6,'(X,"Orbital",3X,"Symmetry")')
      Do I=1,NBasis
      NSymNO(I)=0

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
      Write(6,'(I4,7X,I4)')I,NSymNO(I)
      EndDo
C
C     CALCULATE THE A,C,E MATRICES
C
      Call ACE2el(AMAT,CMAT,EMAT,AMATK,CMATK,EMATK,
     $ URe,Occ,XOne,TwoMOLR,NBasis,NDim,NInte1,NInte2)
C
      If(IFunSR.Ne.0) Then
C
      Write(6,'(/,"*** Adding a sr-kernel... ***")')
C
C     COMPUTE A CONTRIBUTION FROM THE SR-XC KERNEL 
C
      Call GetKerNO(XKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NDimKer,NBasis,NGrid)
C
      IAB=0
      Do IA=1,NBasis
C
      CA=SQRT(Occ(IA))
      If(Occ(IA).Gt.Half) CA=-CA
C
      Do IB=1,IA
C
      CB=SQRT(Occ(IB))
      If(Occ(IB).Gt.Half) CB=-CB
C
      If(IB.Lt.IA) IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
C
      CC=SQRT(Occ(IC))
      If(Occ(IC).Gt.Half) CC=-CC
C
      Do ID=1,IC
C
      CD=SQRT(Occ(ID))
      If(Occ(ID).Gt.Half) CD=-CD
C
      XKer1234=XKer(NAddrrK(IA,IB,IC,ID))
      TwoSR=TwoMO(NAddr3(IA,IB,IC,ID))-TwoMOLR(NAddr3(IA,IB,IC,ID))
C
      If(ID.Lt.IC) ICD=ICD+1
C
      IABCD=(ICD-1)*NDim+IAB
      If(IA.Ne.IB.And.IC.Ne.ID)
     $ AMATK(IABCD)=AMATK(IABCD)+Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      IABC=(IC-1)*NDim+IAB
      If(IA.Ne.IB.And.IC.Eq.ID)
     $ CMATK(IABC)=CMATK(IABC)+Four*(CA+CB)*CC*(XKer1234+TwoSR)
C
      IAC=(IC-1)*NBasis+IA
      If(IA.Eq.IB.And.IC.Eq.ID)
     $ EMATK(IAC)=EMATK(IAC)+Eight*CA*CC*(XKer1234+TwoSR)
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
c herer!!!
C
c      If(I.Gt.MOccTD) Then
c      If(J.Le.MOccTD+MFracTD.And.Occ(J).Gt.1.D-4) Then
c h2
      If(J.Le.nbasis) Then
c      If(J.Le.2.And.(Abs(Occ(I)-Occ(J)).Gt.1.D-4*Occ(I))) Then
c n2
c      If(J.Le.9) Then
c      If(J.Le.10.And.(Abs(Occ(I)-Occ(J)).Gt.1.D-4*Occ(I))) Then
c be
c      If(J.Le.5) Then
c      If(J.Le.2.And.(Abs(Occ(I)-Occ(J)).Gt.1.D-4*Occ(I))) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C      
      EndIf
c      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind       
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
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      AMAT(IJ)=AMAT(IJ1)
      AMATK(IJ)=AMATK(IJ1)
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(J-1)*NDim+IndX(I)
      CMAT(IJ)=CMAT(IJ1)
      CMATK(IJ)=CMATK(IJ1)
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(J-1)*NBasis+I
      EMAT(IJ)=EMAT(IJ1)
      EMATK(IJ)=EMATK(IJ1)
      EndDo
      EndDo
C
C     EXCITATION ENERGIES FROM THE ADIABATIC APPROX.
C
      Call Excit2el(AMAT,CMAT,EMAT,AMATK,CMATK,EMATK,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN) 
C
      Return
      End

*Deck Excit2el
      Subroutine Excit2el(AMAT,CMAT,EMAT,AMATK,CMATK,EMATK,
     $ Occ,IndN,NSymNO,MultpC,NBasis,NDimX,NDimN)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(SMALL=1.D-7,Delta=1.D-4)
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ AMAT(NDimX,NDimX),CMAT(NDimX,NDimN),EMAT(NDimN,NDimN),
     $ AMATK(NDimX,NDimX),CMATK(NDimX,NDimN),EMATK(NDimN,NDimN),
     $ Q1(2*(NDimX+NDimN),2*(NDimX+NDimN)),
     $ Q2(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecR(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ EigVecL(2*(NDimX+NDimN)*2*(NDimX+NDimN)),
     $ Occ(NBasis),
     $ Eig(2*(NDimX+NDimN)),EigI(2*(NDimX+NDimN)),
     $ Work(4*2*(NDimX+NDimN)),
     $ IndN(2,NDimX),NSymI(2*(NDimX+NDimN)),
     $ NSymNO(NBasis),MultpC(15,15)
C
c herer!!!
      Do I=1,2*(NDimX+NDimN)
c      Do I=1,2*NDimX+NDimN
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
c herer!!!
      Do J=1,2*(NDimX+NDimN)
c      Do J=1,2*NDimX+NDimN
      Q1(I,J)=Zero
C
C     1ST SUPER ROW
C
      If(I.Le.NDimX) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) Q1(I,J)=AMAT(I,J-NDimX)
c herer!!!
      If(J.Gt.2*NDimX+NDimN)  Q1(I,J)=CMAT(I,J-2*NDimX-NDimN)
c      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
c     $ Q1(I,J)=CMAT(I,J-2*NDimX)
C
      EndIf
C
C     2ND SUPER ROW
C
      If(I.Gt.NDimX.And.I.Le.2*NDimX) Then
C
      If(J.Le.NDimX) Q1(I,J)=AMATK(I-NDimX,J)
c herer!!!
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=CMATK(I-NDimX,J-2*NDimX)
c       If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
c     $ Q1(I,J)=CMATK(I-NDimX,J-2*NDimX)
C
      EndIf
C
C     3RD SUPER ROW 
C
      If(I.Gt.2*NDimX.And.I.Le.2*NDimX+NDimN) Then
C
      If(J.Gt.NDimX.And.J.Le.2*NDimX) 
     $ Q1(I,J)=Two*CMAT(J-NDimX,I-2*NDimX)
c herer!!!
      If(J.Gt.2*NDimX+NDimN)
     $ Q1(I,J)=EMAT(I-2*NDimX,J-2*NDimX-NDimN)
c      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
c     $ Q1(I,J)=EMATK(I-2*NDimX,J-2*NDimX)
C
      EndIf
C
C     4TH SUPER ROW
C
      If(I.Gt.2*NDimX+NDimN) Then
C
      If(J.Le.NDimX) Q1(I,J)=Two*CMATK(J,I-2*NDimX-NDimN)
      If(J.Gt.2*NDimX.And.J.Le.2*NDimX+NDimN)
     $ Q1(I,J)=EMATK(I-2*NDimX-NDimN,J-2*NDimX)
C
      EndIf
C
      EndDo
      EndDo
c herer!!!
c      NI=2*(NDimX+NDimN)
c      Call DGEEV('N','V',NI,Q1,NI,Eig,EigI,
c     $           EigVecL,NI,EigVecR,NI,Work,4*NI,INFO)
c 
c      do i=1,ni
c      write(*,*)i,eig(i)
c
c      xnorm=zero
c      do j=1,2*NDimX+NDimN
c      xnorm=xnorm+eigvecr((i-1)*NI+j)**2
c      enddo
c      xnorm=one/sqrt(xnorm)
c      
c      do j=1,ni
c      if(j.le.ndimx) Then
c      write(*,*)indn(1,j),indn(2,j),eigvecr((i-1)*NI+j)*xnorm
c      elseif(j.le.2*ndimx) then
c      write(*,*)indn(1,j-ndimx),indn(2,j-ndimx),eigvecr((i-1)*NI+j)
c     $ *xnorm
c      else
c      write(*,*)j-2*ndimx,j-2*ndimx,eigvecr((i-1)*NI+j)*xnorm
c      endif
c      enddo
c      enddo
c
cc      Call SortEig(1,Eig,EigI,EigVecR,NI)
c      stop   
c

C
      Write(6,'(/,''Excitation Energies in [au] and [eV]'')')
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
C
c herer!!!
      Do J=1,2*(NDimX+NDimN)
c      Do J=1,2*NDimX+NDimN
      JSym=NSymI(J)
C
c herer!!!
      Do I=1,2*(NDimX+NDimN)
c      Do I=1,2*NDimX+NDimN

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
      If(EigI(I).Ne.Zero) Im=1
C
      Write(6,'("Sym = ",I1,I5,5X,2F12.5,I2)')NSym,I,Eig(I),
     $ Eig(I)*27.211,Im
C
      EndDo
      Write(6,*)
C
  999 Continue
      EndDo
C
      If(NISum.Ne.2*(NDimX+NDimN)) Write(*,*) 
     $ 'Sum of Dims of symmetrized matrices different from NDimX'
C
      Return
      End

*Deck ACE2el
      Subroutine ACE2el(AMAT,CMAT,EMAT,AMATK,CMATK,EMATK,
     $ URe,Occ,XOne,TwoMO,NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE MATRICES NEEDED FOR 2-EL CASE
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension AMAT(NDim,NDim),CMAT(NDim,NBasis),
     $ EMAT(NBasis,NBasis),
     $ AMATK(NDim,NDim),CMATK(NDim,NBasis),
     $ EMATK(NBasis,NBasis),
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
C     COMPUTE THE ENERGY
C
      ETot=Zero
      Do IP=1,NBasis
C
      CP=SQRT(Occ(IP))
      If(Occ(IP).Gt.Half) CP=-CP
      IPP=IP*(IP-1)/2+IP
      ETot=ETot+Two*Occ(IP)*HNO(IPP)
C
      Do IQ=1,IP
      CQ=SQRT(Occ(IQ))
      If(Occ(IQ).Gt.Half) CQ=-CQ
      Fac=Two
      If(IP.Eq.IQ) Fac=One
C
      TwoEl=TwoMO(NAddr3(IP,IQ,IP,IQ))
      ETot=ETot+Fac*CP*CQ*TwoEl
C
      EMAT(IP,IQ)=TwoEl
      If(IP.Eq.IQ) EMAT(IP,IQ)=EMAT(IP,IQ)+Two*HNO(IPP)
C
      EndDo
C
      EndDo
C
C     COMPUTE THE MATRICES 
C
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IP.Eq.IQ) EMAT(IP,IQ)=EMAT(IP,IQ)-ETot
      EMAT(IQ,IP)=EMAT(IP,IQ)
      EMATK(IP,IQ)=EMAT(IP,IQ)
      EMATK(IQ,IP)=EMAT(IP,IQ)
      EndDo
      EndDo     
C
      IAC=0
      Do IA=1,NBasis
      Do IC=1,IA-1
      IAC=IAC+1
      HAC=HNO((Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC))
C
      Do IB=1,NBasis
      CMAT(IAC,IB)=TwoMO(NAddr3(IA,IB,IC,IB))
      If(IB.Eq.IC.Or.IB.Eq.IA) CMAT(IAC,IB)=CMAT(IAC,IB)+HAC
      CMATK(IAC,IB)=CMAT(IAC,IB)
      EndDo 
C
      IDB=0
      Do ID=1,NBasis
      Do IB=1,ID-1
      IDB=IDB+1
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      IBC=(Max(IB,IC)*(Max(IB,IC)-1))/2+Min(IB,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IDC=(Max(ID,IC)*(Max(ID,IC)-1))/2+Min(ID,IC)
C
      AMAT(IAC,IDB)=TwoMO(NAddr3(IA,ID,IC,IB))
     $ +TwoMO(NAddr3(IA,IB,IC,ID))
C
      If(ID.Eq.IA.And.IB.Eq.IC) AMAT(IAC,IDB)=AMAT(IAC,IDB)-ETot
      If(IB.Eq.IA.And.ID.Eq.IC) AMAT(IAC,IDB)=AMAT(IAC,IDB)-ETot
C
      If(ID.Eq.IC) AMAT(IAC,IDB)=AMAT(IAC,IDB)+HNO(IAB)
      If(ID.Eq.IA) AMAT(IAC,IDB)=AMAT(IAC,IDB)+HNO(IBC)
      If(IB.Eq.IC) AMAT(IAC,IDB)=AMAT(IAC,IDB)+HNO(IAD)
      If(IB.Eq.IA) AMAT(IAC,IDB)=AMAT(IAC,IDB)+HNO(IDC)
C
      AMATK(IAC,IDB)=AMAT(IAC,IDB)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Return
      End











 
      
