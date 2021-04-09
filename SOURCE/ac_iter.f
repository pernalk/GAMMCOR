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
      Call WInteg(ECorrA,TwoNO,XOne,URe,Occ,
     $ ABPLUS,ABMIN,
     $ EGOne,NGOcc,
     $ NBasis,NInte1,NInte2,NDim,NGem,IndAux,ACAlpha,
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
      Dimension XGrid(100), WGrid(100)
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
c      Call ERPASYMM(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
C
C     Frequency integration of CMAT 
C
      NGrid=20
      X0=0.5D0
      Call GauLeg(-1.D0,1.D0,XGrid,WGrid,NGrid)
C   
      COM=0.0 
      Do IGL=1,NGrid
C
      OmI=X0*(One+XGrid(IGL))/(One-XGrid(IGL))
      AIN=Zero
      Do I=1,NDimX
      AIN((I-1)*NDimX+I)=One
      EndDo
C     ABPLUS*ABMIN+1 Om^2
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,
     $           ABMIN,NDimX,OmI**2,AIN,NDimX)

      Call dgetrf(NDimX, NDimX, AIN, NDimX, ipiv, inf1 )
      Call dgetri(NDimX, AIN, NDimX, ipiv, work, NDimX, inf2 )
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,AIN,NDimX,
     $           ABPLUS,NDimX,0.0,CMAT,NDimX)
C
      COM=COM+2.D0/PI*Two*X0*WGrid(IGL)/(One-XGrid(IGL))**2*CMAT
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


