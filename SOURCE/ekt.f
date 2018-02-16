*Deck EKT
      Subroutine EKT(URe,Occ,XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
C     COMPUTES THE EKT IONIZATION POTENTIALS 
C
      Parameter(Zero=0.0D0,One=1.0D0,Half=0.5D0)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),XKin(NInte1),
     $          XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension XLag(NBasis,NBasis),EigVal(NBasis),Work(NBasis)
C
C     SORT THE OCCUPATION NUMBERS
C
      Call SortOcc(Occ,URe,NBasis)
C
      Call Lagrangian(XLag,URe,Occ,XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C
      Do I=1,NBasis
      Do J=1,I
      XLag(J,I)=XLag(I,J)
      EndDo
      EndDo
C
      Call Diag8(XLag,NBasis,NBasis,EigVal,Work) 
C
      write(*,*)eigval  
C
      Return
      End

*Deck Lagrangian 
      Subroutine Lagrangian(XLag,URe,Occ,
     $ XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C
C     CONSTRUCT A LAGRANGIAN MATRIX
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
C
      Dimension XLag(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension F(NInte1),Gamma(NInte1),Xikkj(NBasis,NInte1)
C
C     CONSTRUCT NEW DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
C     CALCULATE THE COULOMBIC CONTRIBUTION TO F AND THE INTEGRALS <IK|KJ> IN NO
C
      Call TIKKJ(Xikkj,F,Gamma,URe,TwoEl,NBasis,NInte1)
C
C     ADD THE ONE-ELECTRON CONTRIBUTION TO F
C
      Do 20 I=1,NInte1
   20 F(I)=F(I)+XKin(I)+XNuc(I)
C
C     TRANSFORM F FROM AO TO NO
C
      Call MatTr(F,URe,NBasis)
C
      Do I=1,NBasis
      Do J=1,NBasis
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      XLag(I,J)=Occ(I)*F(IJ)
C
      Do 40 K=1,NBasis
      XKIJ=Xikkj(K,IJ)
      XLag(I,J)=XLag(I,J)+GOCC(Occ(I),Occ(K),0,I,K)*XKIJ
   40 Continue
C
      XLag(I,J)=XLag(I,J)/Sqrt(Occ(I)*Occ(J))
      if(i.ge.j) write(*,*)i,j,Sqrt(Occ(I)*Occ(J))*xlag(i,j),
     $  Sqrt(Occ(I)*Occ(J))*xlag(j,i)
C
      EndDo
      EndDo
C
      Return
      End     

