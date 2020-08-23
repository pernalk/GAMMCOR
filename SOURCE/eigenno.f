*Deck EigenNO
      Subroutine EigenNO(ETot,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $ NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
C     FOR A GIVEN SET OF NON'S SOLVE THE EIGENEQUATION TO GET NEW NO'S
C
      Parameter(Zero=0.0D0,One=1.0D0,Half=0.5D0)
C
      Parameter(Delta=1.D-4, StLam=0.1D0,FacS=2.5D0, XLMin=1.D-3, 
     $ ETol=1.D-8, MxIt=50)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),XKin(NInte1),
     $          XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension F(NInte1),XU(NBasis,NBasis),UNew(NBasis,NBasis),
     $          EigVal(NBasis),Work(NBasis),FunIJK(NInte1,NBasis),
     $          FunIK(NInte1)
C
C     SORT THE OCCUPATION NUMBERS
C
      Call SortOcc(Occ,URe,NBasis)
C
C     CALCULATE DIFFERENT FUNCTIONS OF OCC
C
      Call GetFunIJ(FunIK,FunIJK,Occ,Delta,NBasis)
C
C     BEGIN THE SCF LOOP
C
      ETotO=Zero
C
      Call CpyM(UNew,URe,NBasis)
C
      XLam=One
      Iter=0
  999 Iter=Iter+1
C
C     CONSTRUCT A DENSITY MATRIX, A FOCK MATRIX IN NO, AND COMPUTE THE TOTAL ENERGY
C
      ETotO=ETot
      Call FockM(ETot,F,UNew,Occ,FunIJK,FunIK,
     $ XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C  
      If(ETot.Lt.ETotO) Then
      Call CpyM(URe,UNew,NBasis)
      If(Iter.Gt.1.And.ETotO-ETot.Le.ETol) GoTo 888
C
      Else
      Call CpyM(UNew,URe,NBasis)
      Call FockM(ETot,F,UNew,Occ,FunIJK,FunIK,
     $ XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
      If(XLam.Gt.StLam+Delta) Then
      XLam=XLam-StLam
      Else
      XLam=Half*XLam
      EndIf
C      
      If(IPrint.Ge.3)
     $ Write(6,'(''Energy increased! XLam reduced to'',F6.2)')XLam
      If(XLam.Lt.XLMin) GoTo 888 
      EndIf
C
C     SHIFT THE LEVELS AND DAMP THE OFF-DIAGONAL ELEMENTS TO F
C
      GradNo=Zero
      IJ=0
      Do 20 I=1,NBasis
      Do 20 J=1,I
      IJ=IJ+1
      If(I.Eq.J) Then
C  
      F(IJ)=I*FacS
      Else
      GradNo=GradNo+((Occ(I)-Occ(J))*F(IJ))**2
      F(IJ)=XLam*F(IJ)
      EndIf
   20 Continue 
      GradNo=Sqrt(GradNo)
C
      If(IPrint.Ge.2)
     $ Write(6,'(/,''EIGENNO: ITERATION:'',I3,3X,''ENERGY'',E16.8,3X,
     $ ''GRAD NORM'',E11.3,3X,''ENERGY ERR'',E10.2)')
     $ Iter,ETot,GradNo,ETot-ETotO
C
      Call CpySym(XU,F,NBasis)    
C
      Call Diag8(XU,NBasis,NBasis,EigVal,Work)
C
C     COMPUTE NEW ORBITALS 
C
      Call MultpM(UNew,XU,URe,NBasis)
C
      If(Iter.Lt.MxIt) GoTo 999
C
  888 Continue
      If(IPrint.Ge.1)
     $ Write(6,'(/,X,''EigNo is quitting with the energy'',
     $ E16.8)')ETot
C
      Return
      End

*Deck FockM 
      Subroutine FockM(ETot,F,URe,Occ,FunIJK,FunIK,
     $ XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C
C     CONSTRUCT A DENSITY MATRIX, A FOCK MATRIX IN NO, AND COMPUTE THE TOTAL ENERGY
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
C
      Dimension F(NInte1),URe(NBasis,NBasis),Occ(NBasis),
     $ FunIJK(NInte1,NBasis),FunIK(NInte1),XKin(NInte1),
     $ XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension Gamma(NInte1),Xikkj(NBasis,NInte1)
C
C     CONSTRUCT NEW DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
C     CALCULATE THE COULOMBIC CONTRIBUTION TO F AND THE INTEGRALS <IK|KJ> IN NO
C
      Call TIKKJ(Xikkj,F,Gamma,URe,TwoEl,NBasis,NInte2)
C
C     ADD THE ONE-ELECTRON CONTRIBUTION TO F
C
      Do 20 I=1,NInte1
   20 F(I)=F(I)+XKin(I)+XNuc(I)
C
C     COMPUTE THE ONE-ELECTRON AND COULOMB COMPONENTS OF ETOT
C
      ETot=Zero
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
   30 ETot=ETot+FacIJ*Gamma(IJ)*(XKin(IJ)+XNuc(IJ)+F(IJ))
C
C     TRANSFORM F FROM AO TO NO
C
      Call MatTr(F,URe,NBasis)
C
C     ADD THE EXCHANGE-CORRELATION PART TO F AND ETot
C
      IJ=0
      Do 40 I=1,NBasis
      IK=I*(I-1)/2
      Do 40 J=1,I
      IJ=IJ+1
C
      Do 40 K=1,NBasis
C
      XKIJ=Xikkj(K,IJ)
C
      F(IJ)=F(IJ)+FunIJK(IJ,K)*XKIJ
C
      If(I.Eq.J.And.K.Le.I) Then
      IK=IK+1
      Fac=Two
      If(I.Eq.K) Fac=One
      ETot=ETot+Fac*FunIK(IK)*XKIJ
      EndIf
C
   40 Continue
C
      Return
      End     

*Deck TIKKJ 
      Subroutine TIKKJ(Xikkj,F,Gamma,URe,TwoEl,NBasis,NInte2)
C
C     CALCULATE THE COULOMBIC CONTRIBUTION TO F AND THE INTEGRALS <IK|KJ> IN NO
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Parameter (Tol2=1.D-12)
C
      Dimension Xikkj(NBasis,NBasis*(NBasis+1)/2),
     $ F(NBasis*(NBasis+1)/2),
     $ Gamma(NBasis*(NBasis+1)/2),URe(NBasis,NBasis),TwoEl(NInte2)
C
C     LOCAL ARRAY     
C
      Dimension Hlp(NBasis,NBasis,NBasis)
C
      Include 'commons.inc'
C
      NInte1=NBasis*(NBasis+1)/2
C
      Do 5 I=1,NBasis
      Do 5 J=1,NBasis
      Do 5 K=1,NBasis
    5 Hlp(K,J,I)=Zero
C
      Do 10 I=1,NInte1
      F(I)=Zero
      Do 10 J=1,NBasis
  10  Xikkj(J,I)=Zero
C
      NAddr=0
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 40 K=1,NBasis
      Do 40 L=1,K
      KL=KL+1
C
      If(IJ.Lt.KL) GoTo 40
      NAddr=NAddr+1
      TwoZet=Two*TwoEl(NAddr)
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      FacIJ=One
      If((I.Eq.J).And.(K.Ne.L)) FacIJ=Two
      FacKL=One
      If((K.Eq.L).And.(I.Ne.J)) FacKL=Two
      FacIL=One
      If((I.Eq.L).And.(K.Ne.J)) FacIL=Two
      FacJK=One
      If((J.Eq.K).And.(I.Ne.L)) FacJK=Two
      FacIK=One
      If((I.Eq.K).And.(L.Ne.J)) FacIK=Two
      FacJL=One
      If((J.Eq.L).And.(I.Ne.K)) FacJL=Two
C
      IL=(Max(I,L)*(Max(I,L)-1))/2+Min(I,L)
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      JL=(Max(J,L)*(Max(J,L)-1))/2+Min(J,L)
C
      If (IFun.Ne.1) Then
      F(IJ)=F(IJ)+FacIJ*Gamma(KL)*TwoZet
      If(IJ.Ne.KL) Then
      F(KL)=F(KL)+FacKL*Gamma(IJ)*TwoZet
      EndIf
      EndIf
C
      Do IA=1,NBasis
C
      Xikkj(IA,IL)=Xikkj(IA,IL)+Half*FacIL*URe(IA,J)*URe(IA,K)*TwoZet
      If(IL.Ne.JK) Xikkj(IA,JK)=Xikkj(IA,JK) 
     $ +Half*FacJK*URe(IA,I)*URe(IA,L)*TwoZet
C
      EndDo
C
      If(K.Eq.L.Or.I.Eq.J) GoTo 40
C
      If(IFun.Ne.1) Then
      F(IJ)=F(IJ)+Gamma(KL)*TwoZet
      If(IJ.Ne.KL) Then
      F(KL)=F(KL)+Gamma(IJ)*TwoZet
      EndIf
      EndIf
C
      Do IA=1,NBasis
C
      Xikkj(IA,IK)=Xikkj(IA,IK)+Half*FacIK*URe(IA,J)*URe(IA,L)*TwoZet
      If(IK.Ne.JL) Xikkj(IA,JL)=Xikkj(IA,JL) 
     $ +Half*FacJL*URe(IA,I)*URe(IA,K)*TwoZet
C
      EndDo
C
      EndIf
C
   40 Continue
   30 Continue
C
      IAB=0
      Do 50 IA=1,NBasis
      Do 50 IB=1,IA
      IAB=IAB+1
      Do 50 I=1,NBasis
      HlpIAB=Xikkj(I,IAB)
C
      Do 60 J=1,NBasis
      Hlp(J,I,IA)=Hlp(J,I,IA)+HlpIAB*URe(J,IB)
      If(IA.Ne.IB) Hlp(J,I,IB)=Hlp(J,I,IB)+HlpIAB*URe(J,IA)
   60 Continue
C
   50 Continue
C
      Do 70 I=1,NInte1
      Do 70 J=1,NBasis
   70 Xikkj(J,I)=Zero
C
      Do 80 I=1,NBasis
      JK=0
      Do 80 J=1,NBasis
      Do 80 K=1,J
      JK=JK+1
      Do 80 IA=1,NBasis
   80 Xikkj(I,JK)=Xikkj(I,JK)+URe(K,IA)*Hlp(J,I,IA)
C
      Return
      End

*Deck SortOcc
      Subroutine SortOcc(Occ,URe,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension UReOld(NBasis,NBasis),Ind(1000),IndOcc(1000)
C
C     SORT THE OCCUPATION NUMBERS IN A DESCENDING ORDER
C
      Do I=1,NBasis
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,NBasis
C
      OccMx=Occ(IStart)
      IndMx=IStart
C
      Do J=IStart,NBasis
      If(Occ(J).Gt.OccMx) Then
      OccMx=Occ(J)
      IndMx=J
      EndIf
      EndDo
C
      Hlp=Occ(IStart)
      IndHlp=Ind(IStart)

      Occ(IStart)=OccMx
      Ind(IStart)=Ind(IndMx)

      Occ(IndMx)=Hlp
      Ind(IndMx)=IndHlp
C
      IStart=IStart+1
C
      EndDo
C
C     SWAP THE ORBITALS
C
      Do I=1,NBasis
      Do J=1,NBasis
      UReOld(J,I)=URe(J,I)
      EndDo
      EndDo
C
      Do J=1,NBasis
      Do I=1,NBasis
      URe(I,J)=UReOld(Ind(I),J)
      EndDo
      EndDo
C
C     CHANGE THE INDICES OF THE BONDING ORBITALS
C
      If(IFun.Eq.6) Then
C
      Do I=1,NBasis
      IndOcc(I)=IType(I)
      EndDo
      Do I=1,NBasis
      IType(I)=IndOcc(Ind(I))
      EndDo
C
      EndIf
C 
C     CHANGE THE ORDER OF GEMINAL INDICES
C
      If(IFun.Eq.13) Then
C
      Do I=1,NBasis
      IndOcc(I)=IGem(I)
      EndDo
      Do I=1,NBasis
      IGem(I)=IndOcc(Ind(I))
      EndDo
C
      Do I=1,NBasis
      UreOld(1,I)=CICoef(I)
      EndDo
      Do I=1,NBasis
      CICoef(I)=UreOld(1,Ind(I))
      EndDo
C
      EndIf
C
      Return
      End

*Deck GetFunIJ
      Subroutine GetFunIJ(FunIK,FunIJK,Occ,Delta,NBasis)
C
C     CALCULATE FUNCTIONS OF OCC
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0)
C
      Dimension FunIK(NBasis*(NBasis+1)/2),
     $ FunIJK(NBasis*(NBasis+1)/2,NBasis),Occ(NBasis)
C
      IJ=0
      Do 40 I=1,NBasis
      Do 40 J=1,I
      IJ=IJ+1
C
      IDer=0
      If(Abs(Occ(I)-Occ(J)).Le.Delta*Occ(I)) IDer=1
C
      IK=I*(I-1)/2
      Do 40 K=1,NBasis
C
      HlpF=GOCC(Occ(I),Occ(K),0,I,K)
      If(K.Le.I) Then
      IK=IK+1
      FunIK(IK)=HlpF
      EndIf
C
      FunIJK(IJ,K)=Zero
      If(IDer.Eq.0) FunIJK(IJ,K)=
     $ (HlpF-GOCC(Occ(J),Occ(K),0,J,K))/(Occ(I)-Occ(J))
C
   40 Continue
C
      Return
      End


