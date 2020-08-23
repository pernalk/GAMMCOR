*Deck EKT
      Subroutine EKT(URe,Occ,XOne,TwoNO,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
C     COMPUTES THE EKT IONIZATION POTENTIALS 
C
      Parameter(Zero=0.0D0,One=1.0D0,Half=0.5D0)
      Parameter(toeV=27.21138602D0)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1),
     $          TwoNO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension XLag(NBasis*NBasis),EigVal(NBasis),Work(NBasis)
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
C
      Call LagrCAS(XLag,URe,Occ,XOne,TwoNO,
     $ NBasis,NInte1,NInte2)
C
C     SYMMETRIZE XLag
      Do J=2,NOccup 
      Do I=1,J-1
      Val = Half*(XLag((I-1)*NOccup+J)+XLag((J-1)*NOccup+I))
      XLag((I-1)*NOccup+J) = Val
      XLag((J-1)*NOccup+I) = Val
      EndDo
      EndDo
C
      Call Diag8(XLag,NOccup,NOccup,EigVal,Work) 
C
      Write(6,'(/,
     $ " *** EKT Ionization Energies (a.u., eV) *** ")')
C
      Do I=1,NOccup
      Write(6,'(I4,6X,2E16.6)') I,EigVal(I),toeV*EigVal(I)
      EndDo
C
      Write(6,'(1X,A,F13.8,3X,F13.8)') 'Lowest IP:',
     $ -One*EigVal(NOccup),-One*toeV*EigVal(NOccup)
C
      Return
      End

*Deck LagrCAS
      Subroutine LagrCAS(XLag,URe,Occ,XOne,TwoNO,
     $ NBasis,NInte1,NInte2)
C
C     CONSTRUCT A LAGRANGIAN MATRIX
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Include 'commons.inc'
C
      Dimension XLag(NBasis*NBasis),URe(NBasis,NBasis),Occ(NBasis),
     $ XOne(NInte1),TwoNO(NInte2)
C
C     LOCAL ARRAYS
      Real*8, Allocatable :: RDM2Act(:)
      Dimension HNO(NInte1),Ind2(NBasis)
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
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,'(4I4,F19.12)',End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      GoTo 10
   40 Continue
      Close(10)
C
      XLag = 0
C      
      Do I=1,NOccup
      Do J=1,NOccup
C
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      XLag(I+(J-1)*NOccup)=Occ(I)*HNO(IJ)
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
C
      XLag(I+(J-1)*NOccup)=XLag(I+(J-1)*NOccup)
     $ +FRDM2(IQ,IR,I,IP,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(J,IQ,IP,IR))
C
      EndDo
      EndDo
      EndDo      
C
      XLag(I+(J-1)*NOccup)=XLag(I+(J-1)*NOccup)/Sqrt(Occ(I)*Occ(J))
C 
      EndDo
      EndDo
C
      Deallocate (RDM2Act)
C
      Write(*,*) 'XLag-Ka',norm2(XLag)
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

