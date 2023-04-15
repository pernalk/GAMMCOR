*Deck GetETot
      Real*8 Function GetETot(Occ,HNO,CoulNO,ExchNO,NBasis)
C
C     COMPUTE THE TOTAL ELECTRONIC ENERGY
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Dimension Occ(NBasis),HNO(NBasis*(NBasis+1)/2),
     $ CoulNO(NBasis*(NBasis+1)/2),ExchNO(NBasis*(NBasis+1)/2)
C
      Include 'commons.inc'
C
      XKu=Zero
      If(IFun.Eq.1) XKu=One
C
      ETot=Zero
      IJ=0
      Do 20 I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*HNO(II)   
C
      Do 20 J=1,I
      IJ=IJ+1
      GemIJ=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(J)) GemIJ=Zero
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
   20 ETot=ETot+FacIJ*(GemIJ*Occ(I)*Occ(J)*CoulNO(IJ)*(One-XKu)
     $ +Half*GOCC(Occ(I),Occ(J),0,I,J)*ExchNO(IJ))
C
      GetETot=Two*ETot
C
      Return
      End

*Deck EneDir
      Real*8 Function EneDir(X,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
C     COMPUTE THE TOTAL ELECTRONIC ENERGY AT A POINT OCC+X*DIR WITH
C     IMPOSING OF THE CORRECT NORMALIZATION (DECREASE X IF NECESSARY)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0, Half=0.5D0,MxIt=25)
      Dimension Dir(NBasis),Occ(NBasis)
C
C     LOCAL ARRAY
C
      Dimension P(NBasis)
C
      Common/PEN/ Pen
C
      It=0 
  999 Continue
      It=It+1
      If(It.Gt.MxIt) Then
      Write(6,'("Fatal error in EneDir!")')
      X=Zero
      Return
      EndIf
C
      Do 10 I=1,NBasis
      XI=ACOS(SQRT(Occ(I))) 
      P(I)=Cos(XI+X*Dir(I))**2
   10 Continue
C
      If(Pen.Eq.Zero) Then
      Call NormN(P,NBasis,IFAIL)
C
      If(IFAIL.Eq.1) Then
      X=Half*X
      GoTo 999
      EndIf
C
      EndIf
C
      EneDir=GetETot(P,HNO,CoulNO,ExchNO,NBasis)
C
C
C     ADD A PENEALTY FUNCTION IF REQUIRED 
C
      If(Pen.Ne.Zero) Then
C
      Sum=Zero
      Do 20 I=1,NBasis
   20 Sum=Sum+P(I)
C
      Sum=(Sum-XELE)**2 
      EneDir=EneDir+Pen*Sum
C
      EndIf
C
      Return
      End

*Deck GradN
      Subroutine GradN(Grad,HNO,Occ,CoulNO,ExchNO,NBasis,NInte1,IFl)
C
C     COMPUTE THE GRADIENT OF ETOT WITH RESPECT TO ni
C  
C     ETot= 2 Sum_ ni h_ii + 4 Sum_ij n_i n_j J_ij + 2 Sum_ij GOCC(n_i,n_j) K_ij 
C
C     IFl = 0  -  ignore (put to zero) gradients corresponding to 
C                 occupancies close to 0 or 1
C         = 1  -  calculate gradients for all occupancies 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, One=1.D0, Two=2.D0, Four=4.D0)
C
C     IF OCC<Small0 OR 1-OCC<Small1 SET GRAD TO ZERO
C
      Parameter(Small0=1.D-8,Small1=1.D-7)
C
      Dimension
     $ Grad(NBasis),HNO(NInte1),Occ(NBasis),CoulNO(NInte1),
     $ ExchNO(NInte1)
C
      Include 'commons.inc'
C
      XKu=Zero
      If(IFun.Eq.1) XKu=One
C
      Do 10, I=1,NBasis
C
      If((Occ(I).Gt.Small0.And.One-Occ(I).Gt.Small1).Or.IFl.Eq.1) Then
C      
      II=I*(I+1)/2    
      Grad(I)=HNO(II)
C
      Do 20 J=1,NBasis
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      GemIJ=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(J)) GemIJ=Zero
      Grad(I)=Grad(I)+Two*GemIJ*Occ(J)*CoulNO(IJ)*(One-XKu)
     $ +GOCC(Occ(I),Occ(J),1,I,J)*ExchNO(IJ)
   20 Continue
C
      Grad(I)=Two*Grad(I)
C 
      Else
      Grad(I)=Zero
      EndIf
C
   10 Continue 
C
      Return
      End

*Deck HessNN
      Subroutine HessNN(Hess,Occ,CoulNO,ExchNO,NBasis,NInte1)
C
C     COMPUTE HESSIAN FOR THE OCCUPATION NUMBERS 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, One=1.D0, Two=2.D0, Four=4.D0)
C
C     IF OCC<Small0 OR 1-OCC<Small1 SET GRAD TO ZERO
C
      Parameter(Small0=1.D-8,Small1=1.D-7)
C
      Dimension
     $ Hess(NBasis,NBasis),Occ(NBasis),CoulNO(NInte1),ExchNO(NInte1)
C
      Include 'commons.inc'
C
      XKu=Zero
      If(IFun.Eq.1) XKu=One
C
      Do 5 I=1,NBasis
      Do 5 J=1,I
    5 Hess(I,J)=Zero
C
      Do 10 I=1,NBasis
      If(Occ(I).Gt.Small0.And.One-Occ(I).Gt.Small1) Then
C
      Do 20 J=1,I
      If(Occ(J).Gt.Small0.And.One-Occ(J).Gt.Small1) Then
C
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      Hess(I,J)=Four*CoulNO(IJ)*(One-XKu)+Two*GOCC(Occ(I),Occ(J),3,I,J)
     $ *ExchNO(IJ)
C
      EndIf
   20 Continue 
C
      EndIf
   10 Continue
C
      Do 30 I=1,NBasis
      If(Occ(I).Gt.Small0.And.One-Occ(I).Gt.Small1) Then
C
      Do 40 J=1,NBasis
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
   40 Hess(I,I)=Hess(I,I)+Two*GOCC(Occ(I),Occ(J),2,I,J)*ExchNO(IJ)
C
      EndIf
   30 Continue
C
      Do 50 I=1,NBasis
      Do 50 J=1,I
   50 Hess(J,I)=Hess(I,J)
C
      Return
      End

*Deck GOCC
      Real*8 Function GOCC(X,Y,IFlag,I,J)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Small=1.D-8, XUpper=0.95D0, XLower=0.05D0)
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Three=3.D0,
     $  Four=4.D0)
C
      Include 'commons.inc'
C
C     IFlag =   0 - calculate F(x,y)
C               1 - calculate dF(x,y)/dx
C               2 - calculate d^2 F(x,y)/d^2 x
C               3 - calculate d^2 F(x,y)/dx/dy
C
      IType(I)=1
      IType(J)=1
      If(IFun.Eq.6) Then
      If(X.Lt.XUpper.And.X.Gt.XLower) IType(I)=0
      If(Y.Lt.XUpper.And.Y.Gt.XLower) IType(J)=0
      EndIf
C
      If(X.Lt.Small.Or.Y.Lt.Small) Then
      GOCC=Zero
      Return
      EndIf
C
C     ** APSG **
C
      If(IFun.Eq.13) Then
C
      If(IGem(I).Eq.IGem(J)) Then
C
      FacX=One
      If(CICoef(I).Lt.Zero) FacX=-One
      FacY=One
      If(CICoef(J).Lt.Zero) FacY=-One
C
C     0th
C
      If(IFlag.Eq.0) Then
C
      GOCC=FacX*FacY*SQRT(X*Y)
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=Half*FacX*FacY*SQRT(Y/X)
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=-FacX*FacY/Four*SQRT(Y/X)/X
C
      Else
      GOCC=FacX*FacY/Four/SQRT(X*Y)
C
      EndIf
C
      EndIf
C
      Return
C
c     else If(IGem(I).Eq.IGem(J))
      Else
C
C     0th
C
      If(IFlag.Eq.0) Then
      GOCC=-X*Y
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Y
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=Zero
C
      Else
      GOCC=-One
C
      EndIf
C
      EndIf
C
      Return
C
c     If(IGem(I).Eq.IGem(J))
      EndIf
C
C     ** End Of APSG **
C
c     If(IFun.Eq.13)
      EndIf
C
C     ** KU **
C
      If(IFun.Eq.1) Then
C
      FacX=One
      If(X.Lt.Half) FacX=-One
      FacY=One
      If(Y.Lt.Half) FacY=-One
C
C     0th
C
      If(IFlag.Eq.0) Then
C
      GOCC=FacX*FacY*SQRT(X*Y)
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=Half*FacX*FacY*SQRT(Y/X) 
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=-FacX*FacY/Four*SQRT(Y/X)/X
C
      Else
      GOCC=FacX*FacY/Four/SQRT(X*Y)
C
      EndIf
C
      EndIf
C
C     ** End Of KU **
C
      Return
      EndIf
C
C     ** GPFBB **
C to be used with ERPA for GPF wavefunction (e.g. CASSCF)
C
      If(IFun.Eq.22) Then
C
C     0th
C
      If(IFlag.Eq.0) Then
C    
      If(IGem(I).Eq.IGem(J)) Then
      GOCC=-SQRT(X*Y)
      Else
      GOCC=-X*Y
      EndIf
C
      Else
C
      Stop 'IFun=22 only available for IFlag=0'
C
      EndIf
C
C     ** End Of GPFBB **
C
      Return
      EndIf
C
C     ** BB **
C
      If(IFun.Eq.2) Then
C
C     0th
C
      If(IFlag.Eq.0) Then
      GOCC=-SQRT(X*Y)
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Half*SQRT(Y/X)
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=One/Four*SQRT(Y/X)/X
C
      Else
      GOCC=-One/Four/SQRT(X*Y)
C
      EndIf
C
      EndIf
C
C     ** End Of BB **
C
      Return
      EndIf
C
C     ** HF **
C
      If(IFun.Eq.12) Then
C
C     0th
C
      If(IFlag.Eq.0) Then
      GOCC=-X*Y
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Y
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=Zero
C
      Else
      GOCC=-One
C
      EndIf
C
      EndIf
C
C     ** End Of HF **
C
      Return
      EndIf
C
C
C     ** BB+HF **
C
      If(IFun.Eq.8) Then
C
C     0th
C
      If(IFlag.Eq.0) Then
      GOCC=-(One-Cmix)*SQRT(X*Y)-Cmix*X*Y
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Half*(One-Cmix)*SQRT(Y/X)-Cmix*Y
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=(One-Cmix)/Four*SQRT(Y/X)/X
C
      Else
      GOCC=-(One-Cmix)/Four/SQRT(X*Y)-Cmix
C
      EndIf
C
      EndIf
C
C     ** End Of BB+HF **
C
      Return
      EndIf
C
C     
C     ** CHF **
C
      If(IFun.Eq.7) Then
C
      SX=X*Abs(One-X)
      SY=Y*Abs(One-Y)
C
C     0th
C     
      If(IFlag.Eq.0) Then
      GOCC=-(X*Y + SQRT(SX*SY))
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Y
      If(SX.Eq.Zero) Return
      GOCC=GOCC-Half*(One-Two*X)*SQRT(SY)/SQRT(SX)
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=Zero
      If(SX.Eq.Zero) Return
      GOCC=One/Four*SQRT(SY)/(SX)**(Three/Two)
C
      Else
      GOCC=Zero
      If(SX.Eq.Zero.Or.SY.Eq.Zero) Return
      GOCC=-(One+Half*Half*(One-Two*X)/SQRT(SX)
     $ *(One-Two*Y)/SQRT(SY))
C     
      EndIf
C     
      EndIf
C     
C     ** End Of CHF **
C
      Return
      EndIf
C
C     ** GU **
C
      If(IFun.Eq.3) Then
C
      If(IFlag.Eq.0) Then
C
C     0th
C
      GOCC=-SQRT(X*Y)
      If(I.Eq.J) Then
      GOCC=GOCC+X-X**2
      EndIf
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Half*SQRT(Y/X)
      If(I.Eq.J) GOCC=-X
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=One/Four*SQRT(Y/X)/X
      If(I.Eq.J) GOCC=-Half
C
      Else
      GOCC=-One/Four/SQRT(X*Y)
      If(I.Eq.J) GOCC=-Half
C
      EndIf
C
      EndIf
C
C     ** End Of GU **
C
      Return
      EndIf
C
C     ** BBC1 **
C
      If(IFun.Eq.4) Then
C
      If(IFlag.Eq.0) Then
C
C     0th
C
      GOCC=-SQRT(X*Y)
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half) GOCC=SQRT(X*Y)
C
      Else
C
C     1st
C
      If(IFlag.Eq.1) Then
      GOCC=-Half*SQRT(Y/X)
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=Half*SQRT(Y/X)
C
C     2nd
C
      ElseIf(IFlag.Eq.2) Then
      GOCC=One/Four*SQRT(Y/X)/X
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=-One/Four*SQRT(Y/X)/X
C
      Else
      GOCC=-One/Four/SQRT(X*Y)
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=One/Four/SQRT(X*Y)
C
      EndIf
C
      EndIf
C
C     ** End Of BBC1
C
      Return
      EndIf
C
C     ** BBC2 OR BBC3 **
C
      If(IFun.Eq.5.Or.IFun.Eq.6) Then
C
      If(IFlag.Eq.0) Then
C
C     *** 0th ***
C
      GOCC=-SQRT(X*Y)
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=SQRT(X*Y)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Ge.Half.And.Y.Ge.Half.And.IType(I)*IType(J).Eq.1)
     $ GOCC=GOCC+Sqrt(X*Y)-X*Y
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(X.Ge.Half.And.IType(I).Eq.1.And.IType(J).Eq.0) 
     $ GOCC=GOCC+Sqrt(X*Y)-X*Y
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(Y.Ge.Half.And.IType(J).Eq.1.And.IType(I).Eq.0) 
     $ GOCC=GOCC+Sqrt(X*Y)-X*Y
C
C     BOTH ORBITALS ARE STRONGLY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) GOCC=GOCC+X-X**2 
C
      Else
C
      If(IFlag.Eq.1) Then
C
C     *** 1st ***
C
      GOCC=-Half*SQRT(Y/X)
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=Half*SQRT(Y/X)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Ge.Half.And.Y.Ge.Half.And.IType(I)*IType(J).Eq.1)
     $ GOCC=-Y
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(X.Ge.Half.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ GOCC=-Y
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(Y.Ge.Half.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ GOCC=-Y
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) GOCC=-X
C
      ElseIf(IFlag.Eq.2) Then
C
C     *** 2nd ***
C
      GOCC=One/Four*SQRT(Y/X)/X
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=-One/Four*SQRT(Y/X)/X
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Ge.Half.And.Y.Ge.Half.And.IType(I)*IType(J).Eq.1)
     $ GOCC=Zero
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(X.Ge.Half.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ GOCC=Zero
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(Y.Ge.Half.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ GOCC=Zero
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) GOCC=-Half
C
C     *** 2nd mixed ***
C
      Else
C
      GOCC=-One/Four/SQRT(X*Y)
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Lt.Half.And.Y.Lt.Half)
     $ GOCC=One/Four/SQRT(X*Y)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.X.Ge.Half.And.Y.Ge.Half.And.IType(I)*IType(J).Eq.1)
     $ GOCC=-One
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(X.Ge.Half.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ GOCC=-One
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(Y.Ge.Half.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ GOCC=-One
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) GOCC=-Half
C
      EndIf
C
      EndIf
C
C     End Of BBC2 BBC3
C
      Return
      EndIf
C  
C
C     ** BBC4 **
C
      If(IFun.Eq.7) Then
C
      IndI=1
      If(X.Lt.Half) IndI=0
      IndJ=1
      If(Y.Lt.Half) IndJ=0
C
      If(IFlag.Eq.0) Then
C
C     0th
C     
      If(I.Ne.J) Then 
      Select Case(IndI+IndJ)
      Case(2) 
              GOCC=-SQRT(X*Y)
      Case(0) 
              GOCC=SQRT(X*Y)
      Case(1) 
              GOCC=FBBC4(0,0,X,Y)
      EndSelect
      Else
      GOCC=FBBC4(0,1,X,X)
      EndIf
C
      ElseIf(IFlag.Eq.1) Then
C
C     1ST DERIVATIVE
C
      If(I.Ne.J) Then
      Select Case(IndI+IndJ)
      Case(2) 
              GOCC=-Half*SQRT(Y/X)
      Case(0) 
              GOCC=Half*SQRT(Y/X)
      Case(1) 
              GOCC=FBBC4(1,0,X,Y)
      EndSelect
      Else
      GOCC=Half*FBBC4(1,1,X,X)
      EndIf
C
C     2ND DERIVATIVES
C
      ElseIf(IFlag.Eq.2) Then
C
      If(I.Ne.J) Then
      Select Case(IndI+IndJ)
      Case(2) 
              GOCC=One/Four*SQRT(Y/X)/X
      Case(0) 
              GOCC=-One/Four*SQRT(Y/X)/X
      Case(1) 
              GOCC=FBBC4(2,0,X,Y)
      EndSelect
      Else
      GOCC=FBBC4(2,1,X,X)/Four
      EndIf
C
      Else
C
      If(I.Ne.J) Then
      Select Case(IndI+IndJ)
      Case(2) 
              GOCC=-One/Four/SQRT(X*Y)
      Case(0) 
              GOCC=One/Four/SQRT(X*Y)
      Case(1) 
              GOCC=FBBC4(3,0,X,Y)
      EndSelect 
      Else
      GOCC=FBBC4(3,1,X,X)/Four
      EndIf
C
      EndIf
C
      EndIf
C
C     ** End Of BBC4
C
      Return
      End 

*Deck FBBC4
      Real*8 Function FBBC4(IFlag1,IFlag2,X,Y)
C
C     CALCULATE VALUES OR DERIVATIVES OF THE BBC4 TERMS FOR THE
C     MIXED (OCCUPIED-VIRTUAL) OR DIAGONAL TERMS
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Three=3.D0,
     $  Four=4.D0)
C
      Include 'commons.inc'
C
C     IFlag1 =  0 - calculate F(x,y)
C               1 - calculate dF(x,y)/dx
C               2 - calculate d^2 F(x,y)/d^2 x
C               3 - calculate d^2 F(x,y)/dx/dy
C
C     IFlag2  = 0 - calculate an off-diagonal term     
C             = 1 - calculate a diagonal term
C
C     OFF-DIAGONAL CASE 
C
      If(IFlag2.Eq.0) Then  
C
      If(X.Gt.Y) Then
      FX=EXP(-AFUN1*(X-BFUN1))
      FXA=AFUN1*FX/(One+FX)
      FXAA=AFUN1*FXA
      FY=EXP(-AFUN2*(Y-BFUN2))
      FYA=AFUN2*FY/(One+FY)
      Else
      FX=EXP(-AFUN2*(X-BFUN2))
      FXA=AFUN2*FX/(One+FX)
      FXAA=AFUN2*FXA
      FY=EXP(-AFUN1*(Y-BFUN1))
      FYA=AFUN1*FY/(One+FY)
      EndIf
C
      HXY=One/(One+FX)/(One+FY)
      HFX=FX/(One+FX)
      HFY=FY/(One+FY)
C
      Select Case(IFlag1)
      Case(0) 
              FBBC4 = -SQRT(X*Y) + HXY*(SQRT(X*Y)-X*Y)

      Case(1)
              FBBC4 = -Half*SQRT(Y/X) 
     $              + HXY*(Half*SQRT(Y/X)-Y+(SQRT(X*Y)-X*Y)*FXA)

      Case(2) 
              FBBC4 = SQRT(Y/X)/Four/X + HXY*FXA
     $              * ( Half*SQRT(Y/X)-Y + (SQRT(X*Y)-X*Y)*FXA )
     $              + HXY*( -SQRT(Y/X)/Four/X + (Half*SQRT(Y/X)-Y)*FXA
     $                     + (SQRT(X*Y)-X*Y)*FXAA*(HFX-One) )

      Case(3) 
              FBBC4 = -One/Four/SQRT(X*Y) + HXY*(
     $                One/Four/SQRT(X*Y)-One
     $              + (Half*SQRT(X/Y)-X)*FXA
     $              + (Half*SQRT(Y/X)-Y)*FYA
     $              + (SQRT(X*Y)-X*Y)*FXA*FYA )
C
      EndSelect
C
      EndIf
C
C     END OF OFF-DIAGONAL CASE 
C 
C     DIAGONAL CASE
C
      If(IFlag2.Eq.1) Then
C
      FX=EXP(-AFUN1*(X-BFUN1))
      HFX=FX/(One+FX)
      FY=EXP(-AFUN2*(Y-BFUN2))
      HFY=FY/(One+FY)
      HXY=One/(One+FX)/(One+FY) 
C
      Select Case(IFlag1)
      Case(0) 
              FBBC4 = -X + (X-X*X)*(One-FX*HXY)

      Case(1) 
              FBBC4 = -One+(One-Two*X)*(One-FX*HXY)
     $              + (X-X*X)*FX*HXY*(AFUN1 - AFUN1*HFX - AFUN2*HFY)
C
      EndSelect
C
      If(IFlag1.Eq.2.Or.IFlag1.Eq.3) Then
      F3=FX*HXY*(AFUN1 - AFUN1*HFX- AFUN2*HFY )
      FBBC4 = -Two*(One-FX*HXY) + Two*F3*(One-Two*X)
     $      + (X-X*X)*(
     $      - F3*(AFUN1-AFUN1*HFX-AFUN2*HFY )
     $      + FX*HXY*( AFUN1*AFUN1*HFX*(One-HFX)
     $               + AFUN2*AFUN2*HFY*(One-HFY) ) )
      EndIf
C
      EndIf      
C
      Return
      End 

*Deck OPTN
      Subroutine OPTN(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,NBasis)
C
C     FIND OPTIMAL OCCUPATION NUMBERS
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
      Common/CPMFT/ MFrac,MOcc,NFrac
C
      Parameter(Half=0.5D0, Zero=0.D0, One=1.0D0, Two=2.D0, Four=4.D0,
     $ Five=5.D0, Ten=10.D0,Ten1=1.D-1, Ten2=1.D-2, Ten3=1.D-3)
C 
      Parameter(TolP=1.D-10)
C
      Logical IZero(NBasis)
C
      Dimension Occ(NBasis),URe(NBasis,NBasis),
     $ XKin(NBasis*(NBasis+1)/2),XNuc(NBasis*(NBasis+1)/2),
     $ CoulNO(NBasis*(NBasis+1)/2),ExchNO(NBasis*(NBasis+1)/2)
C
C     LOCAL ARRAYS
C
      Dimension
     $ Grad(NBasis),HNO(NBasis*(NBasis+1)/2),
     $ Hlp(NBasis*(NBasis+1)/2),P(NBasis),OccP(NBasis),
     $ Ind(NBasis)
C
      NInte1=NBasis*(NBasis+1)/2
      XNBasis=Float(NBasis)
C
      IJ=0
      Do 10 I=1,NBasis
      Do 10 J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do 10 IA=1,NBasis
      Do 10 IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
   10 HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
C
      If(IPrint.Ge.3) Then
      Write(6,'(2X,'' INITIAL ORBITAL OCCUPANCIES '')')
      Do 7 I=1,NBasis
    7 Write(6,'(1X,I3,E16.6)') I,Occ(I)
      EndIf
C
      ETot=GetETot(Occ,HNO,CoulNO,ExchNO,NBasis)
      If(IPrint.Ge.1)
     $Write(6,'(/,2X,'' INITIAL  ENERGY IN OCC OPTIMIZER'',F16.8)')ETot
C
C     CALL POWELL'S MINIMIZATION ALGORITHM
C
      Do 11 I=1,NBasis
   11 OccP(I)=Occ(I)
C
      If(IPrint.Ge.1)Write(6,'(/,X,''Entering Powells optimization'')')
C
      Do I=1,MFrac
      Ind(I)=MOcc+I
      EndDo
C
      Call Powell(TolP,ETotP,IterP,
     $ OccP,HNO,CoulNO,ExchNO,Ind,MFrac,NBasis)
      If(IPrint.Ge.1)
     $ Write(6,'(X,''Powells converged after'',I4,
     $ '' iterations with the energy equals'',F16.8,/)')
     $ IterP,ETotP
C
C     END OF THE OPTIMIZATION 
C
      Do 75 I=1,NBasis
   75 Occ(I)=OccP(I)
      ETot=ETotP
C
C     CALCULATE FINAL GRADIENT AND MIU AFTER A TIGHT NORMALIZATION
C
      Call GradN(Grad,HNO,Occ,CoulNO,ExchNO,NBasis,NInte1,0)
C
C     CHANGE THE GRADIENT TO TAKE INTO ACCOUNT SYMMETRY
C
      Do I=1,NBasis
C
      Hlp(I)=Zero
      ICount=0
C
      Do J=1,NBasis
      Hlp(I)=Hlp(I)+Grad(J)
      ICount=ICount+1
      EndDo
C
      Hlp(I)=Hlp(I)/Float(ICount)
      EndDo
C
      Do I=1,NBasis
      Grad(I)=Hlp(I)
      EndDo
C
      XMiu=Zero
      NDimH=0
      Do 85 I=1,NBasis
      HlpG=Grad(I)
      IZero(I)=.FALSE.
      If(HlpG.Ne.Zero) Then
      IZero(I)=.TRUE.
      XMiu=XMiu+HlpG
      NDimH=NDimH+1
      EndIf
   85 Continue
      XMiu=XMiu/Float(NDimH)
C
      GradNo=Zero
      DevMiu=Zero
C
      Do 90 I=1,NBasis
C
      HlpG=Grad(I)
C  
      If(IZero(I)) Then
      GradNo=GradNo+(HlpG-XMiu)**2*Occ(I)*(One-Occ(I))
      DevMiu=DevMiu+(HlpG-XMiu)**2
      EndIf
C
   90 Continue
      GradNo=Sqrt(GradNo)
      DevMiu=Sqrt(DevMiu/XNBasis)
C
      If(IPrint.Ge.0) Then
      Write(6,'(X,
     $ ''FINAL RESULTS FROM THE OPTIMIZATION OF THE OCCUPANCIES'')')
      Write(6,'(X,''Total Energy'',F16.8,5X,
     $ ''Norm of Grad'',E10.2,5X,''Standard Deviation of dEdni'',E9.2)')
     $ ETot,GradNo,DevMiu
C
      Write(6,'(/,X,'' ORBITAL OCCUPANCIES AND MIU_i'')')
C
      Do 100 I=1,NBasis
  100 Write(6,'(1X,I3,2E16.6)') I,Occ(I),Grad(I)*Half
C
      EndIf
C
      Return
      End



*Deck OptTwo
      Subroutine OptTwo(ETot,URe,Occ,XOne,TwoEl,NSymMO,
     $ NBasis,NInte1,NInte2,NoEig)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SYSTEMS
C     THE WHOLE DENISTY MATRIX IS FOUND IN ONE STEP
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XOne(NInte1),
     $ TwoEl(NInte2),NSymMO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension APA(NBasis*NBasis,NBasis*NBasis),PPA(NBasis*NBasis),
     $ PWorkA(NBasis*NBasis)
       Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $ indx(nbasis,nbasis),hno(ninte1),singl(nbasis*nbasis)

C      goto 999
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
c herer!!!

      goto 999
      IAB=0
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,NBasis
      ICD=ICD+1
C
      APA(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IAC)
      If(IA.Eq.IC) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IBD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call Diag8(APA,NBasis**2,NBasis**2,PPA,PWorkA)

c      write(*,*)PPa
c      kl=0
c      do k=1,nbasis
c      do l=1,nbasis
c      kl=kl+1 
c      write(*,*)"excit",ppa(kl)-ppa(1)
c      ij=0
c      do i=1,nbasis
c      do j=1,nbasis
c      ij=ij+1
c      write(*,*)i,j,apa(kl,ij)
c
c      enddo
c      enddo
c      enddo
c      enddo
c etot
C
      ETot=Zero
      IJ=0
      Do I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*XOne(II)
     $ -Half*Occ(I)*(One-Occ(I))*Twoel(NAddr3(I,I,I,I))
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Occ(I)*Occ(J)*(Twoel(NAddr3(I,I,J,J))
     $ -Half*Twoel(NAddr3(I,J,I,J)))
      enddo
      enddo  
C
      ETot=Two*ETot
      write(*,*)'ETOT 1 ***',etot

       ij=0
       do i=1,nbasis
       do j=1,nbasis
       ij=ij+1
       indx(i,j)=ij
       enddo
       enddo

c select singlets
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       singl(kl)=1
       do ip=1,nbasis
       do ir=1,ip-1
       if(apa(kl,indx(ip,ir))*apa(kl,indx(ir,ip)).lt.-1.D-8) singl(kl)=0
       enddo
       enddo
       enddo
       enddo

       ipr=0
       do ip=1,nbasis
       do ir=1,nbasis
       ipr=ipr+1
       iqs=0
       do iq=1,nbasis
       do is=1,nbasis
       iqs=iqs+1
c
       cr=sqrt(occ(ir))
       if(occ(ir).gt.0.5) cr=-cr
       cq=sqrt(occ(iq))
       if(occ(iq).gt.0.5) cq=-cq
       cp=sqrt(occ(ip))
       if(occ(ip).gt.0.5) cp=-cp
       cs=sqrt(occ(is))
       if(occ(is).gt.0.5) cs=-cs
c
      if(ip.gt.ir.and.iq.gt.is) then
      if(ir.eq.is.and.ip.eq.iq)
     $ etot=etot-(occ(ip)*(one-occ(is))+occ(is)*(one-occ(ip)))*
     $ twoel(NAddr3(Ip,Ir,Iq,Is))
      endif 
c
       if(ip.gt.ir.and.iq.gt.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1 
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))

       enddo  
       enddo
c
       aux=two*(cr+cp)*(cq+cs)*sum

       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.gt.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=4.d0*(cr+cp)*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.eq.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=2.d0*cp*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif

c
       enddo
       enddo
       enddo
       enddo
       write(*,*)'ETOT',etot
c herer!!!
       return
c      stop
  999 continue
 

C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))+TwoEl(NAddr3(IA,ID,IC,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,
     $''DIAGONALIZATION IN THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')') 
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETot=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(X,''Total Energy'',F16.8,/)') ETot
C
      ISym=-1
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(Abs(PWork(IAB)**2).Gt.0.1D0) Then
C
      If(Isym.Eq.-1) Then
      ISym=MultpC(NSymMO(IA),NSymMO(IB))
      Else
      If(ISym.Ne.MultpC(NSymMO(IA),NSymMO(IB))) ISym=0
      EndIf
C
      EndIf
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(URe,PWork,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      Occ(I)=Occ(I)**2
      Sum=Sum+Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=Occ(I)/Sum
      EndDo
C
      Call SortOcc(Occ,URe,NBasis)
C
      Write(6,'(X,'' ORBITAL OCCUPANCIES'')')
      Do 100 I=1,NBasis
  100 Write(6,'(1X,I3,E16.6)') I,Occ(I)
C
      Return
      End

*Deck OptAPSG
      Subroutine OptAPSG(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem)
C
C     OPTIMIZATION ALGORITHM FOR THE c COEFFICIENTS IN THE APSG FUNCTIONAL
C
      Implicit Real*8 (A-H,O-Z)
C
c herer!!!
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0,
c     $ MxIt=10,ETol=1.D-5)
     $ MxIt=10,ETol=1.D-9)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(NInte1),CoulNO(NInte1),ExchNO(NInte1),
     $ NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),AC(NBasis*NBasis),XMu(NBasis),Work(NBasis),
     $ Ind(NBasis)
C
      IJ=0
      Do 10 I=1,NBasis
      Do 10 J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do 10 IA=1,NBasis
      Do 10 IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
   10 HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
C
      EOld=GetETot(Occ,HNO,CoulNO,ExchNO,NBasis)
C
      If(IPrint.Ge.1)
     $Write(6,'(/,X,'' INITIAL  ENERGY IN APSG OPTIMIZER'',F16.8)')EOld
C
      If(IPrint.Ge.3) Then
      Write(6,'(X,'' INITIAL ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem",3X,"Sym")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I),NSymNO(I)
      EndDo
      EndIf
C
C     BEGINNING OF THE MACROLOOP
C
      Do It=1,MxIt
C
      NGem0=NGem
C
C KP 24.11.14
C
C     IF NGem=NELE+1 -> DO NOT OPTIMIZE OCCUPATIONS OF THE FICTITIOUS GEMINAL
C
      If(NGem.Eq.NELE+1) NGem0=NGem-1
C
      Do II=1,NGem0
C
      NDim=0
      Do I=1,NBasis
      If(IGem(I).Eq.II) NDim=NDim+1
      EndDo 
C
C     FIND OPTIMAL COEFFICIENTS FOR THE II'TH GEMINAL
C
      NCol=0
      Do IQ=1,NBasis
C
      If(II.Eq.IGem(IQ)) Then
C
      NCol=NCol+1
      Ind(NCol)=IQ
C
      NRow=0
      Do IP=1,IQ
C
      If(II.Eq.IGem(IP)) Then
C
      NRow=NRow+1
      IndAC=(NCol-1)*NDim+NRow 
C
      IQP=(Max(IP,IQ)*(Max(IP,IQ)-1))/2+Min(IP,IQ)
C
      AC(IndAC)=ExchNO(IQP)
C
      If(IP.Eq.IQ) Then
C
      AC(IndAC)=AC(IndAC)+Two*HNO(IQP)
      Do IR=1,NBasis
      IPR=(Max(IP,IR)*(Max(IP,IR)-1))/2+Min(IP,IR)
      If(II.Ne.IGem(IR)) AC(IndAC)=AC(IndAC)+Two*Occ(IR)*
     $ (Two*CoulNO(IPR)-ExchNO(IPR))
      EndDo
C
      EndIf
C
      AC((NRow-1)*NDim+NCol)=AC(IndAC)
C
C     If(II.Eq.IGem(IP)) 
      EndIf
C
      EndDo
C
C     If(II.Eq.IGem(IQ))      
      EndIf
C
      EndDo
C
      Call Diag8(AC,NDim,NDim,XMu,Work)
C
      Do I=1,NDim
      CICoef(Ind(I))=AC((I-1)*NDim+1)
      Occ(Ind(I))=CICoef(Ind(I))**2
      EndDo
C
c     enddo of II=1,NGem
      EndDo
C
C     CHECK THE CHANGE IN THE ENERGY
C
      ETot=GetETot(Occ,HNO,CoulNO,ExchNO,NBasis)
C
      If(Abs(EOld-ETot).Le.ETol) Then
C
      If(IPrint.Ge.0) Then
      Write(6,'(/,X,'' FINAL ENERGY IN APSG COEFFICIENTS OPTIMIZER'',
     $ F16.8)')ETot
      Write(6,'(X,'' FINAL ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem",3X,"Sym")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I),NSymNO(I)
      EndDo
      EndIf
C
      Return
C
      EndIf
C
      EOld=ETot  
C 
C     END OF MACRO LOOP
      EndDo
C
      Return
      End

*Deck GuessN 
      Subroutine GuessN(Occ,Eps,URe,XKin,XNuc,CoulNO,ExchNO,NBasis)
C
C     GENERATE A GUESS FOR THE OCCUPATION NUMBERS USING HF ORBITAL ENERGIES
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
      Common/PEN/ Pen
C
      Parameter(Half=0.5D0, Zero=0.D0, One=1.D0, Two=2.D0)
C
      Parameter(TolP=1.D-6)
C
      Dimension Occ(NBasis),Eps(NBasis),URe(NBasis,NBasis),
     $ XKin(NBasis*(NBasis+1)/2),XNuc(NBasis*(NBasis+1)/2),
     $ CoulNO(NBasis*(NBasis+1)/2),ExchNO(NBasis*(NBasis+1)/2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NBasis*(NBasis+1)/2),Dir(NBasis)
C
      Pen=Zero
C
      NInte1=NBasis*(NBasis+1)/2
      XNBasis=Float(NBasis)
C
      IJ=0
      Do 40 I=1,NBasis
      Do 40 J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      IAB=Zero
      Do 40 IA=1,NBasis
      Do 40 IB=1,IA
      IAB=IAB+1
      FacAB=Two
      If(IA.Eq.IB) FacAB=One
   40 HNO(IJ)=HNO(IJ)+FacAB*URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
C
C     COMPUTE A PRIMITIVE GUESS
C
      Do 10 I=1,NELE
   10 Occ(I)=One
C
      Do 20 I=NELE+1,NBasis
   20 Occ(I)=Zero
C
      EpsAv=Half*(Eps(NELE)+Eps(NELE+1))
      Diff=Eps(NELE+1)-Eps(NELE)
C
      EpsSum=Zero
      Do 2 I=1,NELE
    2 EpsSum=EpsSum+One/(Eps(I)-EpsAv)
      Coef=-EpsSum
      EpsSum=Zero
      Do 5 I=NELE+1,NBasis
    5 EpsSum=EpsSum+One/(Eps(I)-EpsAv)
      Coef=Coef/EpsSum
C
      Do 79 I=1,NBasis
      If (I.Le.NELE) Then
      Dir(I)=One/(Eps(I)-EpsAv)*Diff
      Else
      Dir(I)=Coef/(Eps(I)-EpsAv)*Diff
      EndIf
   79 Continue
C
      Call BrackOcc(ax,bx,cx,fa,fb,fc,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
      ETot=BrentOcc(ax,bx,cx,TolP,xmin,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      Do 100 I=1,NBasis
      XI=ACOS(SQRT(Occ(I)))
  100 Occ(I)=Cos(XI+xmin*Dir(I))**2
C
      Call NormN(Occ,NBasis,IFAIL) 
C
      If(IPrint.Ge.1) Then
C
      Write(6,'(/,X,''THE ENERGY CORRESPONDING TO A GUESS'',F16.8,/)')
     $ ETot
      Write(6,'(X,'' GUESS ORBITAL OCCUPANCIES '')')
      Do 7 I=1,NBasis
    7 Write(6,'(1X,I3,E16.6)') I,Occ(I)
C
      EndIf
C
      Return
      End

*Deck NormN
      Subroutine NormN(Occ,NBasis,IFAIL)
C
C     CALCULATES NEW OCCUPANCIES IMPOSING PROPER NORMALIZATION
C     ON FRACTIONALLY OCCUPIED ORBITALS (USED IN CPMFT METHOD)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
      Common/CPMFT/ MFrac,MOcc,NFrac
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,
     $ Tol=1.D-10, MaxIt=50)
C
      Dimension Occ(NBasis)
C
C     LOCAL ARRAY
C
      Dimension P(NBasis)
C
      IFAIL=0
C
      XNEL=NFrac
C
      Do 20 I=1,NBasis
      If(Occ(I)-One.Gt.Zero) Occ(I)=One
   20 P(I)=Abs(ACOS(SQRT(Occ(I))))
C
      Iter=0
  100 Iter=Iter+1
C
      If(Iter.Gt.MaxIt) Then
      IFAIL=1
      Return
      EndIf
C
      C2=Zero
      CS=Zero
      Do 30 I=1,MFrac
      C2=C2+Cos(P(MOcc+I))**2
   30 CS=CS+Cos(P(MOcc+I))*Sin(P(MOcc+I))
C
      Help=-C2+Two*C2*C2+CS**2+XNEL-Two*C2*XNEL
C
      If(Help.Ge.Zero) Then
      Del=(-CS+SQRT(-C2+Two*C2*C2+CS**2+XNEL-Two*C2*XNEL))/
     $ (Two*C2-One)
      Else
      Del=(C2-XNEL)/Two/CS
      EndIf
C
      Sum=Zero
      Do 40 I=1,MFrac
      P(MOcc+I)=Abs(P(MOcc+I)+Del)
      Occ(MOcc+I)=Cos(P(MOcc+I))**2
   40 Sum=Sum+Occ(MOcc+I)
C
      If(Abs(XNEL-Sum).Le.Tol) Return
C
      Goto 100
C
      Return
      End

*Deck BrackOcc 
      Subroutine BrackOcc(ax,bx,cx,fa,fb,fc,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Dir(NBasis),Occ(NBasis)
C
      PARAMETER (GOLD=1.618034D0, GLIMIT=100.D0, TINY=1.d-20)
      PARAMETER (ZERO=0.0D0, TEN5=0.05D0, MXIT=20)
C
      ax=ZERO
      bx=TEN5
C
      fa=EneDir(ax,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
      fb=EneDir(bx,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
C
      cx=bx+GOLD*(bx-ax)
      fc=EneDir(cx,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      its=0
    1 its=its+1
C
      if(its.gt.MXIT) Then
      write(*,*) 'brackocc exceeding no of iterations!'
      return
      endif
C
      if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
C
        if((bx-u)*(u-cx).gt.0.)then
C
          fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
C
          u=cx+GOLD*(cx-bx)
          fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)

          endif
C
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
        else
          u=cx+GOLD*(cx-bx)
          fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
        endif
C
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
C
      return
      END

*Deck LinMin1
      Subroutine LinMin1(DirN,NSym,Fret,Occ,HNO,CoulNO,ExchNO,Ind,
     $ NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Tol=1.D-2)
C
      Common/CPMFT/ MFrac,MOcc,NFrac
C
      Dimension DirN(NSym),Occ(NBasis),Ind(NSym)
C
C     LOCAL ARRAY
C
      Dimension Dir(NBasis)
C
      Dir(1:NBasis)=Zero
C
      Do I=1,NSym
      Dir(Ind(I))=DirN(I)
      EndDo
C
      Call BrackOcc(ax,bx,cx,fa,fb,fc,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C      
      Fret=BrentOcc(ax,bx,cx,Tol,xmin,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      Do 10 I=1,NSym
      XI=ACOS(SQRT(Occ(MOcc+I)))
   10 Occ(MOcc+I)=Cos(XI+xmin*Dir(MOcc+I))**2
C
      Call NormN(Occ,NBasis,IFAIL)
C
      Do 20 I=1,NSym
   20 DirN(I)=xmin*DirN(I)
C
      Return
      End

*Deck Powell
      Subroutine Powell(ftol,fret,iter,Occ,HNO,CoulNO,ExchNO,Ind,
     $ NSym,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
      Common/CPMFT/ MFrac,MOcc,NFrac
C
      Dimension Occ(NBasis),Ind(MFrac)
C
C     LOCAL ARRAYS 
C
      Dimension xi(NSym,NSym),p(NSym),pt(NSym),ptt(NSym),
     $ xit(NSym),phlp(NBasis),Dum(NBasis)
C
      Parameter (ITMAX=250)
C
      n=NSym
C
      Do I=1,NBasis
      Dum(I)=0.D0
      phlp(I)=Occ(I)
      EndDo
C
      Do I=1,NSym
      p(I)=ACOS(SQRT(Occ(Ind(I))))
      EndDo
C
      do 5 i=1,n
      do 5 j=1,n
      xi(i,j)=0.d0
      if(i.eq.j) xi(i,j)=1.d0
    5 continue
C
      fret=GetETot(Occ,HNO,CoulNO,ExchNO,NBasis)
C
      do 11 j=1,n
   11 pt(j)=p(j)
C
      iter=0
    1 iter=iter+1
C
      If(IPrint.Ge.3)
     $ Write(6,'(X,''Iteration:'',I2,3X,''Energy'',F16.8)')iter,fret
C
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
   12   xit(j)=xi(j,i)
        fptt=fret
C
        Call LinMin1(xit,n,fret,Occ,HNO,CoulNO,ExchNO,Ind,NBasis)
        Do J=1,NSym
        p(J)=ACOS(SQRT(Occ(Ind(J))))
        EndDo
C
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
C
   13 continue
C
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret))) return
      if(iter.eq.ITMAX) stop 'powell exceeding maximum iterations'
C
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
   14 continue
C
      do 20 k=1,NSym
   20 phlp(Ind(k))=cos(ptt(k))**2
      Call NormN(phlp,NBasis,IFAIL)
      If(IFAIL.Eq.1) GoTo 1
C
      fptt=EneDir(0.d0,Dum,phlp,HNO,CoulNO,ExchNO,NBasis)
C 
      if(fptt.ge.fp)goto 1
      t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
      if(t.ge.0.)goto 1
C
      Call LinMin1(xit,n,fret,Occ,HNO,CoulNO,ExchNO,Ind,NBasis)
      Do J=1,NSym
      p(J)=ACOS(SQRT(Occ(Ind(J))))
      EndDo
C
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
   15 continue
      goto 1
C
      Return
      End


*Deck BrentOcc
      Real*8 Function BrentOcc(ax,bx,cx,tol,xmin,
     $ Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      PARAMETER (ITMAX=100,CGOLD=.3819660D0,ZEPS=1.0D-10)
C
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
C
      fx=EneDir(x,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
C 
        fu=EneDir(u,Dir,Occ,HNO,CoulNO,ExchNO,NBasis)
C
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      write(*,*)'brentocc exceed maximum iterations'
3     xmin=x
      brentocc=fx
      return
      END

      Subroutine OptTwo1(ETot,ENuc,URe,Occ,XOne,TwoEl,NSymMO,
     $ CISapt,NBasis,NInte1,NInte2,NoEig)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SYSTEMS
C     THE WHOLE DENISTY MATRIX IS FOUND IN ONE STEP
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XOne(NInte1),
     $ TwoEl(NInte2),NSymMO(NBasis),CISapt(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension APA(NBasis*NBasis,NBasis*NBasis),PPA(NBasis*NBasis),
     $ PWorkA(NBasis*NBasis)
       Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $ indx(nbasis,nbasis),hno(ninte1),singl(nbasis*nbasis)
c herer!!!

      goto 999
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
      IAB=0
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,NBasis
      ICD=ICD+1
C
      APA(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IAC)
      If(IA.Eq.IC) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IBD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call Diag8(APA,NBasis**2,NBasis**2,PPA,PWorkA)

c      write(*,*)PPa
c      kl=0
c      do k=1,nbasis
c      do l=1,nbasis
c      kl=kl+1 
c      write(*,*)"excit",ppa(kl)-ppa(1)
c      ij=0
c      do i=1,nbasis
c      do j=1,nbasis
c      ij=ij+1
c      write(*,*)i,j,apa(kl,ij)
c
c      enddo
c      enddo
c      enddo
c      enddo
c etot
C
      ETot=Zero
      IJ=0
      Do I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*XOne(II)
     $ -Half*Occ(I)*(One-Occ(I))*Twoel(NAddr3(I,I,I,I))
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Occ(I)*Occ(J)*(Twoel(NAddr3(I,I,J,J))
     $ -Half*Twoel(NAddr3(I,J,I,J)))
      enddo
      enddo  
C
      ETot=Two*ETot
      write(*,*)'ETOT 1 ***',etot

       ij=0
       do i=1,nbasis
       do j=1,nbasis
       ij=ij+1
       indx(i,j)=ij
       enddo
       enddo

c select singlets
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       singl(kl)=1
       do ip=1,nbasis
       do ir=1,ip-1
       if(apa(kl,indx(ip,ir))*apa(kl,indx(ir,ip)).lt.-1.D-8) singl(kl)=0
       enddo
       enddo
       enddo
       enddo

       ipr=0
       do ip=1,nbasis
       do ir=1,nbasis
       ipr=ipr+1
       iqs=0
       do iq=1,nbasis
       do is=1,nbasis
       iqs=iqs+1
c
       cr=sqrt(occ(ir))
       if(occ(ir).gt.0.5) cr=-cr
       cq=sqrt(occ(iq))
       if(occ(iq).gt.0.5) cq=-cq
       cp=sqrt(occ(ip))
       if(occ(ip).gt.0.5) cp=-cp
       cs=sqrt(occ(is))
       if(occ(is).gt.0.5) cs=-cs
c
      if(ip.gt.ir.and.iq.gt.is) then
      if(ir.eq.is.and.ip.eq.iq)
     $ etot=etot-(occ(ip)*(one-occ(is))+occ(is)*(one-occ(ip)))*
     $ twoel(NAddr3(Ip,Ir,Iq,Is))
      endif 
c
       if(ip.gt.ir.and.iq.gt.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1 
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))

       enddo  
       enddo
c
       aux=two*(cr+cp)*(cq+cs)*sum

       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.gt.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=4.d0*(cr+cp)*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.eq.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=2.d0*cp*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif

c
       enddo
       enddo
       enddo
       enddo
       write(*,*)'ETOT',etot
c herer!!!
       return
c      stop

  999 Continue
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))+TwoEl(NAddr3(IA,ID,IC,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,
     $''DIAGONALIZATION IN THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')') 
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETot=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(/,X,''State no'',I3,'' Total Energy'',F16.10,/)') NoEig,
     $ ETot+ENuc
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(URe,PWork,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      CICoef(I)=Occ(I)
      Occ(I)=Occ(I)**2
      Sum=Sum+Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=Occ(I)/Sum
      CICoef(I)=CICoef(I)/Sqrt(Sum)
      EndDo
C
      Call SortOcc(Occ,URe,NBasis)
C
      Write(6,'('' COEFFICIENTS AND ORBITAL OCCUPANCIES'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,2E16.6)') I,CICoef(I),Occ(I)
C
C     WRITE CICoef FOR SAPT
      Do I=1,NBasis
      CISapt(I) = CICoef(I)
      Enddo

C      Do I=1,NBasis
C         Print*, I, IGem(I)
C      EndDo

      Return
      End

*Deck OptTwo
      Subroutine OptTwo2(ETot,ENuc,URe,Occ,XOne,TwoEl,NSymMO,
     $ CISapt,NBasis,NInte1,NInte2,NoEig)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SYSTEMS
C     THE WHOLE DENISTY MATRIX IS FOUND IN ONE STEP
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XOne(NInte1),
     $ CISapt(Nbasis),TwoEl(NInte2),NSymMO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension APA(NBasis*NBasis,NBasis*NBasis),PPA(NBasis*NBasis),
     $ PWorkA(NBasis*NBasis)
       Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $ indx(nbasis,nbasis),hno(ninte1),singl(nbasis*nbasis)
c herer!!!

      goto 999
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
      IAB=0
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,NBasis
      ICD=ICD+1
C
      APA(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IAC)
      If(IA.Eq.IC) APA(IAB,ICD)=APA(IAB,ICD)+XOne(IBD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call Diag8(APA,NBasis**2,NBasis**2,PPA,PWorkA)

c      write(*,*)PPa
c      kl=0
c      do k=1,nbasis
c      do l=1,nbasis
c      kl=kl+1 
c      write(*,*)"excit",ppa(kl)-ppa(1)
c      ij=0
c      do i=1,nbasis
c      do j=1,nbasis
c      ij=ij+1
c      write(*,*)i,j,apa(kl,ij)
c
c      enddo
c      enddo
c      enddo
c      enddo
c etot
C
      ETot=Zero
      IJ=0
      Do I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*XOne(II)
     $ -Half*Occ(I)*(One-Occ(I))*Twoel(NAddr3(I,I,I,I))
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Occ(I)*Occ(J)*(Twoel(NAddr3(I,I,J,J))
     $ -Half*Twoel(NAddr3(I,J,I,J)))
      enddo
      enddo  
C
      ETot=Two*ETot
      write(*,*)'ETOT 1 ***',etot

       ij=0
       do i=1,nbasis
       do j=1,nbasis
       ij=ij+1
       indx(i,j)=ij
       enddo
       enddo

c select singlets
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       singl(kl)=1
       do ip=1,nbasis
       do ir=1,ip-1
       if(apa(kl,indx(ip,ir))*apa(kl,indx(ir,ip)).lt.-1.D-8) singl(kl)=0
       enddo
       enddo
       enddo
       enddo

       ipr=0
       do ip=1,nbasis
       do ir=1,nbasis
       ipr=ipr+1
       iqs=0
       do iq=1,nbasis
       do is=1,nbasis
       iqs=iqs+1
c
       cr=sqrt(occ(ir))
       if(occ(ir).gt.0.5) cr=-cr
       cq=sqrt(occ(iq))
       if(occ(iq).gt.0.5) cq=-cq
       cp=sqrt(occ(ip))
       if(occ(ip).gt.0.5) cp=-cp
       cs=sqrt(occ(is))
       if(occ(is).gt.0.5) cs=-cs
c
      if(ip.gt.ir.and.iq.gt.is) then
      if(ir.eq.is.and.ip.eq.iq)
     $ etot=etot-(occ(ip)*(one-occ(is))+occ(is)*(one-occ(ip)))*
     $ twoel(NAddr3(Ip,Ir,Iq,Is))
      endif 
c
       if(ip.gt.ir.and.iq.gt.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1 
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))

       enddo  
       enddo
c
       aux=two*(cr+cp)*(cq+cs)*sum

       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.gt.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=4.d0*(cr+cp)*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif
c
       if(ip.eq.ir.and.iq.eq.is) then
c
       sum=zero
       kl=0
       do k=1,nbasis
       do l=1,nbasis
       kl=kl+1
       if(singl(kl).eq.1.and.kl.ne.1) sum=sum+
     $ apa(kl,indx(ip,ir))*apa(kl,indx(is,iq))
       enddo
       enddo
c
       aux=2.d0*cp*cq*sum
       etot=etot+aux*TwoEl(NAddr3(Ip,Ir,Iq,Is))
c
       endif

c
       enddo
       enddo
       enddo
       enddo
       write(*,*)'ETOT',etot
c herer!!!
       return
c      stop

  999 Continue
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))+TwoEl(NAddr3(IA,ID,IC,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,
     $''DIAGONALIZATION IN THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')') 
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETot=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(/,X,''State no'',I3,'' Total Energy'',F19.11,/)') NoEig,
     $ ETot+ENuc
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(URe,PWork,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      CICoef(I)=Occ(I)
      Occ(I)=Occ(I)**2
      Sum=Sum+Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=Occ(I)/Sum
      CICoef(I)=CICoef(I)/Sqrt(Sum)
      EndDo
C
      Call SortOcc(Occ,URe,NBasis)
C
      Write(6,'('' COEFFICIENTS AND ORBITAL OCCUPANCIES'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,2E16.6)') I,CICoef(I),Occ(I)
C
C
C     WRITE CICoef FOR SAPT
      Do I=1,NBasis
      CISapt(I) = CICoef(I)
      Enddo


      Return
      End

*Deck OptTwoP
      Subroutine OptTwoP(ETot,ENuc,URe,Occ,
     $ AP,PP,
     $ XOne,TwoEl,NSymMO,
     $ CISapt,NBasis,NInte1,NInte2,NoEig)
C
C     OPTIMIZATION ALGORITHM FOR TWO-ELECTRON SYSTEMS
C     THE WHOLE DENISTY MATRIX IS FOUND IN ONE STEP
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XOne(NInte1),
     $ TwoEl(NInte2),NSymMO(NBasis),CISapt(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1)
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoEl(NAddr3(IA,IC,ID,IB))+TwoEl(NAddr3(IA,ID,IC,IB))
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+XOne(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,
     $''DIAGONALIZATION IN THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETot=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(/,X,''State no'',I3,'' Total Energy'',F16.10,/)') NoEig,
     $ ETot+ENuc
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(URe,PWork,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      CICoef(I)=Occ(I)
      Occ(I)=Occ(I)**2
      Sum=Sum+Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=Occ(I)/Sum
      CICoef(I)=CICoef(I)/Sqrt(Sum)
      EndDo
C
      Call SortOcc(Occ,URe,NBasis)
C
      Write(6,'('' COEFFICIENTS AND ORBITAL OCCUPANCIES'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,2E16.6)') I,CICoef(I),Occ(I)
C
C     WRITE CICoef FOR SAPT
      Do I=1,NBasis
      CISapt(I) = CICoef(I)
      EndDo
C
C     TRANSFORM ALL P VECTORS (STATES) TO NO's OF THE NoEig's STATES
C
      Do INU=1,NInte1 
C
C     SAVE IN PP EXCITATION ENERGIES FROM THE STATE NoEig
C
      PP(INU)=PP(INU)-ETot
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(INU,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      AP(INU,IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      AP(INU,IJ)=AP(INU,IJ)+SQRT(Half)*URe(I,IA)*URe(J,IB)*PWork(IAB)
      EndDo
      EndDo
C
C check
c      If(PP(INU).Lt.Zero) Write(*,*)'check AP',I,J,AP(INU,IJ)
C
      EndDo   
      EndDo
C
C     GET THE NORM AND PRINT (FOR CHECKING)
C
c      PNorm=Zero
c      IAB=0
c      Do IA=1,NBasis
c      Do IB=1,IA
c      IAB=IAB+1
c      Factor=Two
c      If(IA.Eq.IB) Factor=One
c      PNorm=PNorm+Factor*AP(INU,IAB)**2
c      EndDo
c      EndDo
C
c      Write(*,*)'INU, ExcitEn, PNorm',INU,PP(INU),PNorm 
C     INU
      EndDo
C
      Return
      End

*Deck OptTwoPAlph
      Subroutine OptTwoPAlph(ETot,ENuc,URe,Occ,
     $ AP,PP,
     $ XOne,TwoEl,NSymMO,
     $ CISapt,NBasis,NInte1,NInte2,NoEig,ACAlpha)
C
C     DIAGONALIZATION OF THE AC-HAMILTONIAN FOR A GIVEN ALPHA
C
C     RETURNS:
C     AP - P MATRICES IN THE REPRESENTATION OF THE NoEig-TH STATE NO's
C     PP - TRANSITION ENERGIES WITH RESPECT TO THE NoEig-th STATE 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XOne(NInte1),
     $ TwoEl(NInte2),NSymMO(NBasis),CISapt(NBasis),MultpC(15,15)
C
      Dimension AP(NInte1,NInte1),PP(NInte1),PWork(NInte1),
     $H1Alph(NInte1),HNO(NInte1),TwoElAlph(NInte2)
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      Do I=1,NInte1
      H1Alph(I)=XOne(I)
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      H1Alph(IJ)=ACAlpha*H1Alph(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoEl(NAddr3(IT,IT,I,J))-TwoEl(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      H1Alph(IJ)=H1Alph(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN      
C
      Do I=1,NInte2
      TwoElAlph(I)=TwoEl(I)
      EndDo
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ TwoElAlph(NAdd)=ACAlpha*TwoElAlph(NAdd)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      FAB=One
      If(IA.Eq.IB) FAB=SQRT(Half)
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
      ICD=ICD+1
      FCD=One
      If(IC.Eq.ID) FCD=SQRT(Half)
C
      AP(IAB,ICD)=TwoElAlph(NAddr3(IA,IC,ID,IB))
     $                +TwoElAlph(NAddr3(IA,ID,IC,IB))
C
      IAC=(Max(IA,IC)*(Max(IA,IC)-1))/2+Min(IA,IC)
      IAD=(Max(IA,ID)*(Max(IA,ID)-1))/2+Min(IA,ID)
      IBC=(Max(IC,IB)*(Max(IC,IB)-1))/2+Min(IC,IB)
      IBD=(Max(ID,IB)*(Max(ID,IB)-1))/2+Min(ID,IB)
C
      If(IB.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IAC)
      If(IA.Eq.ID) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IBC)
      If(IB.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IAD)
      If(IA.Eq.IC) AP(IAB,ICD)=AP(IAB,ICD)+H1Alph(IBD)
C
      AP(IAB,ICD)=FAB*FCD*AP(IAB,ICD)
c herer!!!
c      If(.Not.(
c     $IGem(IA).Eq.IGem(IB).And.IGem(IB).Eq.IGem(IC).
c     $ And.IGem(IC).Eq.IGem(ID))) AP(IAB,ICD)=Zero
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,
     $''DIAGONALIZATION IN THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Call Diag8(AP,NInte1,NInte1,PP,PWork)
      ETot=PP(NoEig)
C
      Write(6,'(/,X,
     $ ''RESULTS FROM THE TWO-ELECTRON FUNCTIONAL OPTIMIZATION'')')
C
      Write(6,'(/,X,''State no'',I3,'' Total Energy'',F16.10,/)') NoEig,
     $ ETot+ENuc
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(NoEig,IAB)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
      Call CpySym(URe,PWork,NBasis)
      Call Diag8(URe,NBasis,NBasis,Occ,PWork)
C
      Sum=Zero
      Do I=1,NBasis
      CICoef(I)=Occ(I)
      Occ(I)=Occ(I)**2
      Sum=Sum+Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=Occ(I)/Sum
      CICoef(I)=CICoef(I)/Sqrt(Sum)
      EndDo
C
      Call SortOcc(Occ,URe,NBasis)
C
      Write(6,'('' COEFFICIENTS AND ORBITAL OCCUPANCIES'')')
      Do 100 I=1,NBasis
  100 Write(6,'(X,I3,2E16.6)') I,CICoef(I),Occ(I)
C
C     WRITE CICoef FOR SAPT
      Do I=1,NBasis
      CISapt(I) = CICoef(I)
      EndDo
C
C     TRANSFORM ALL P VECTORS (STATES) TO NO's OF THE NoEig's STATES
C
      Do INU=1,NInte1
C
C     SAVE IN PP EXCITATION ENERGIES FROM THE STATE NoEig
C
      PP(INU)=PP(INU)-ETot
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      PWork(IAB)=AP(INU,IAB)
c herer!!!
c      if(abs(PP(INU)+ETot).gt.1.d-10.and.abs(AP(INU,IAB)).gt.1.d-7)
c     $ write(*,*)inu,ia,ib,PP(INU)
C
      If(IA.Eq.IB) PWork(IAB)=SQRT(Two)*PWork(IAB)
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      AP(INU,IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      AP(INU,IJ)=AP(INU,IJ)+SQRT(Half)*URe(I,IA)*URe(J,IB)*PWork(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     INU
      EndDo
C
      Return
      End



