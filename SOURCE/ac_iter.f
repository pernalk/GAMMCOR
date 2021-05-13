*Deck ACIter
      Subroutine ACIter(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDim,NGOcc,NGem,
     $ IndN,IndX,NDimX)
C
      use abmat
      use abfofo
C 
C     AC Iteratively
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),
     $ XOne(NInte1),IndAux(NBasis),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ EigVecR(NDim*NDim),
     $ Eig(NDim),EGOne(NGem),
     $ UNOAO(NBasis,NBasis),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension XGrid(100), WGrid(100)
C
C
      Call WIter(ECorr,TwoNO,XOne,URe,Occ,
     $ EGOne,NGOcc,
     $ NBasis,NInte1,NInte2,NDim,NGem,IndAux,
     $ IndN,IndX,NDimX)
      ETot=EGOne(1)
      Write
     $ (6,'(/,2X,''ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
      stop
C
C     GENERATE ABSCISSAS AND WEIGHTS FOR GAUSSIAN-LEGENDRE QUADRATURE
C
      NGrid=5
C
      Call GauLeg(Zero,One,XGrid,WGrid,NGrid)
C 
      ECorr=Zero
      Do I=1,NGrid
C   
      ACAlpha=XGrid(I)
C
      Call WInteg(ECorrA,TwoNO,XOne,URe,Occ,ABPLUS,ABMIN,
     $ EGOne,NGOcc,NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
      Write(*,*)'ACAlpha ',ACAlpha,' W_ALPHA ',ECorrA
C
      ECorr=ECorr+WGrid(I)*ECorrA
C
      EndDo
C
      ETot=EGOne(1)
      Write
     $ (6,'(/,2X,''ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      Call DelInts(ITwoEl)
C
      Return
      End

*Deck WInteg
      Subroutine WInteg(ECorr,TwoNO,XOne,URe,Occ,
     $ ABPLUS,ABMIN,
     $ EGOne,NGOcc,
     $ NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
     $ IndN,IndX,NDimX)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ Occ(NBasis),TwoNO(NInte2),XONe(NInte1),URe(NBasis,NBasis),
     $ IndAux(NBasis),
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ EGOne(NGem),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension IPair(NBasis,NBasis)
      Dimension CMAT(NDim*NDim),AIN(NDim*NDim),COM(NDim*NDim)
      Dimension ipiv(NDim),work(NDim)
      Dimension  XFreq(100),WFreq(100)
C
      PI = 4.0*ATAN(1.0)
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
c
      EGOne(1)=ECASSCF
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1)
      EndDo
      EndDo
C
C     Frequency integration of CMAT 
C
      NGrid=18
      Call FreqGrid(XFreq,WFreq,NGrid)
C   
      COM=0.0 
      Do IGL=1,NGrid
C
      OmI=XFreq(IGL)
      AIN=Zero
      Do I=1,NDimX
      AIN((I-1)*NDimX+I)=One
      EndDo
C     ABPLUS*ABMIN+1 Om^2
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,
     $           ABMIN,NDimX,OmI**2,AIN,NDimX)

c      Call dgetrf(NDimX, NDimX, AIN, NDimX, ipiv, inf1 )
c      Call dgetri(NDimX, AIN, NDimX, ipiv, work, NDimX, inf2 )
c      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
c     $           ABPLUS,NDimX,0.0,CMAT,NDimX)

      CMAT=ABPLUS
      Call dgesv(NDimX,NDimX,AIN,NDimX,ipiv,CMAT,NDimX,inf)
C
      COM=COM+2.D0/PI*CMAT*WFreq(IGL)
C
      EndDo
C
C     W Integrand at Alpha
C
      ECorr=Zero
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IR=IndN(2,I)
C
      Do J=1,NDimX
C
      IQ=IndN(1,J)
      IS=IndN(2,J)
C
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ)).And.IP.Gt.IR.And.IQ.Gt.IS) Then
C
      Aux=(CICoef(IS)+CICoef(IQ))*(CICoef(IP)+CICoef(IR))
     $ *COM((J-1)*NDimX+I)
      If(IR.Eq.IS.And.IP.Eq.IQ)
     $ Aux=Aux
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
      EndDo
      EndDo
C
      Return
      End

*Deck WIter
      Subroutine WIter(ECorr,TwoNO,XOne,URe,Occ,
     $ EGOne,NGOcc,
     $ NBasis,NInte1,NInte2,NDim,NGem,IndAux,
     $ IndN,IndX,NDimX)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Dimension
     $ Occ(NBasis),TwoNO(NInte2),XONe(NInte1),URe(NBasis,NBasis),
     $ IndAux(NBasis),
     $ ABPLUS0(NDim*NDim),ABMIN0(NDim*NDim),
     $ ABPLUS1(NDim*NDim),ABMIN1(NDim*NDim),
     $ EGOne(NGem),
     $ IndX(NDim),IndN(2,NDim)
C
C     LOCAL ARRAYS
C
      Dimension IPair(NBasis,NBasis)
      Dimension AIN(NDimX*NDimX),COM(NDimX*NDimX)
      Dimension ipiv(NDimX),work(NDimX)
      Dimension XFreq(100),WFreq(100)
      Dimension A0(NDimX*NDimX),A0PLUS(NDimX*NDimX),AMBDA(NDimX*NDimX),
     $ A1(NDimX*NDimX),A2(NDimX*NDimX),
     $ C0(NDimX*NDimX),C1(NDimX*NDimX),C2(NDimX*NDimX),C3(NDimX*NDimX)
C
      PI = 4.0*ATAN(1.0)
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
      ACAlpha=Zero
      Call AB_CAS(ABPLUS0,ABMIN0,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
      ACAlpha=One
      Call AB_CAS(ABPLUS1,ABMIN1,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
C
      ABPLUS1=ABPLUS1-ABPLUS0
      ABMIN1=ABMIN1-ABMIN0 
C
      EGOne(1)=ECASSCF
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS0(IJ)=ABPLUS0(IJ1)
      ABMIN0(IJ)=ABMIN0(IJ1)
      ABPLUS1(IJ)=ABPLUS1(IJ1)
      ABMIN1(IJ)=ABMIN1(IJ1)
      EndDo
      EndDo
      write(*,*)'matrices reduced'
C
C     A0=ABPLUS0*ABMIN0      
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,
     $           ABMIN0,NDimX,0.0,A0,NDimX)
C     A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,
     $           ABMIN1,NDimX,0.0,A1,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,
     $           ABMIN0,NDimX,1d0,A1,NDimX)
C     A2=ABPLUS1*ABMIN1
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,
     $           ABMIN1,NDimX,0.0,A2,NDimX)      
C
c      NGrid=18
      NGrid=25
      Call FreqGrid(XFreq,WFreq,NGrid)
C   
      COM=0.0 
      Do IGL=1,NGrid
C
      OmI=XFreq(IGL)

      AIN=Zero
      Do I=1,NDimX
      AIN((I-1)*NDimX+I)=OmI**2
      EndDo
C
      AMBDA=A0+AIN
      Call dgetrf(NDimX, NDimX, AMBDA, NDimX, ipiv, inf1 )
      Call dgetri(NDimX, AMBDA, NDimX, ipiv, work, NDimX, inf2 )
C
C     C0=0.5 LAMBDA*ABPLUS0
      Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,AMBDA,NDimX,
     $           ABPLUS0,NDimX,0.0,C0,NDimX)
C     C1=0.5 LAMBDA*ABLUS1-LAMBDA*A1*C0
c      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,A1,NDimX,
c     $           C0,NDimX,0.0,AIN,NDimX)
c      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AMBDA,NDimX,
c     $           AIN,NDimX,0.0,C1,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AMBDA,NDimX,
     $           A1,NDimX,0.0,AIN,NDimX)  
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
     $           C0,NDimX,0.0,C1,NDimX) 
      Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,AMBDA,NDimX,
     $           ABPLUS1,NDimX,-1.d0,C1,NDimX)
C
C     LAMBDA*A1 in AIN
C     C2=-2 AIN*C1 - 2 LAMBDA*A2*C0
      Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,AMBDA,NDimX,
     $           A2,NDimX,0.0,ABMIN0,NDimX)
C     FROM NOW ON: LAMBDA*A2 in ABMIN0, LAMBDA*A1 in AIN
CCCCCCCCCCCCCCCCCCCCC
CC     C(2)
C      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
C     $           C1,NDimX,0.0,C2,NDimX)
C      Call dgemm('N','N',NDimX,NDimX,NDimX,-2.d0,ABMIN0,NDimX,
C     $           C0,NDimX,-2.d0,C2,NDimX)      
C
CC     C(3)
C      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
C     $           C2,NDimX,0.0,C3,NDimX)
C      Call dgemm('N','N',NDimX,NDimX,NDimX,-6.d0,ABMIN0,NDimX,
C     $           C1,NDimX,-3.d0,C3,NDimX)
C
C      COM=COM+4.D0/PI*(C0+0.5D0*C1+C2/6.d0+C3/24.d0)*WFreq(IGL)
CCCCCCCCCCCCCCCCCCCCC
      WFact=4.D0/PI*WFreq(IGL)
C
      COM=COM+WFact*(C0+0.5D0*C1)
C
      XFactorial=1
      Do N=2,4
C
      XFactorial=XFactorial*N
C     C(n)
      XN1=-N
      XN2=-N*(N-1)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
     $           C1,NDimX,0.0,C2,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,XN2,ABMIN0,NDimX,
     $           C0,NDimX,XN1,C2,NDimX)
C
      FF=WFact/XFactorial/(N+1)
      COM=COM+FF*C2
C
      C0=C1
      C1=C2
C
      EndDo
C
C     IGL
      EndDo
C
C     W Integrand 
C
      ECorr=Zero
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IR=IndN(2,I)
C
      Do J=1,NDimX
C
      IQ=IndN(1,J)
      IS=IndN(2,J)
C
      If(.NOT.(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ)).And.IP.Gt.IR.And.IQ.Gt.IS) Then
C
      Aux=(CICoef(IS)+CICoef(IQ))*(CICoef(IP)+CICoef(IR))
     $ *COM((J-1)*NDimX+I)
      If(IR.Eq.IS.And.IP.Eq.IQ)
     $ Aux=Aux
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
      EndDo
      EndDo
C
      Return
      End

*Deck FreqGrid
      Subroutine FreqGrid(XFreq,WFreq,NFreq)
C
      Implicit Real*8 (A-H,O-Z)
      Dimension XFreq(NFreq),WFreq(NFreq)
C
      Call GauLeg(-1.D0,1.D0,XFreq,WFreq,NFreq) 
C
      X0=0.5D0
C
      Do IGL=1,NFreq
      WFreq(IGL)=2.*X0*WFreq(IGL)/(1.-XFreq(IGL))**2
      XFreq(IGL)=X0*(1.+XFreq(IGL))/(1.-XFreq(IGL))
      EndDo
C
      Return
      End
