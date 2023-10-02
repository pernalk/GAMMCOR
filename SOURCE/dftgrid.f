*Deck EPotSR
      Subroutine EPotSR(EnSR,EnHSR,VSR,Occ,URe,
     $ UNOAO,Transp,
     $ OrbGrid,OrbXGrid,OrbYGrid,
C    HAP
C     $ OrbZGrid,WGrid,NSymMO,TwoEl,TwoElErf,
     $ OrbZGrid,WGrid,NSymMO,VCoul,
C     $ TwoEl,TwoElErf,
     $ Omega,Flag,NGrid,NInte1,NInte2,NBasis) 
C
C     RETURNS SR ENERGY AND POTENTIAL = SR_XC + SR_H
C
      use abmat
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
C
      Real(8) Omega
      Integer Flag
      Logical Transp 
      Dimension VSR(NInte1),OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $ Occ(NBasis),URe(NBasis,NBasis),
     $ UNOAO(NBasis,NBasis),
C    HAP
     $ TwoEl(NInte2),TwoElErf(NInte2),
     $ VCoul(NInte1),
     $ NSymMO(NBasis)
      Dimension OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),
     $ OrbZGrid(NGrid,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension VHSR(NInte1),Gamma(NInte1)
C
C     HAP:
      Alpha = Omega
      IFunSR = Flag
C
      If(Alpha.Gt.1.D2.And.IFunSR.Eq.1) 
     $ Stop 'Fatal Error in EPotSR: srLDA procedure is not reliable 
     $ for the range-separation parameter alpha > 100'
C
C     COMPUTE THE HARTREE POTENTIAL AND THE ENERGY
C     HAP: VHSR is computed in AO during canonicalization
C
      Call tran_matTr(VCoul,UNOAO,UNOAO,NBasis,Transp)
      VHSR = VCoul
C     Old:
C      Call PotHSR(VHSR,Occ,URe,TwoEl,TwoElErf,NInte1,NInte2,NBasis)
C      Print*,'VHSR-Ka', norm2(VHSR)
C
      Call VecTr(Gamma,Occ,URe,NBasis)
      EnHSR=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Fac=Two
      If(I.Eq.J)Fac=One
      EnHSR=EnHSR+Fac*Gamma(IJ)*VHSR(IJ)
      EndDo
      EndDo
C
C     COMPUTE THE XC PARTS OF THE POTENTIAL AND THE ENERGY 
C
      If(IFunSR.Eq.1) Then
C
      Call PotXCSR_LDA(VSR,Occ,URe,OrbGrid,
     $ OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,NSymMO,NGrid,NInte1,NBasis)
      Call GetExcSR_LDA(EnxcSR,Occ,URe,OrbGrid,WGrid,NGrid,NBasis)
C
      ElseIf(IFunSR.Eq.2) Then
C
      Call GetExcSR_PBE(EnxcSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NSymMO,NGrid,NInte1,NBasis)
C
      ElseIf(IFunSR.Eq.3) Then
C
       Call GetExc_PBE(EnxcSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NSymMO,NGrid,NInte1,NBasis)
C
      EndIf
C
      Do I=1,NInte1
      VSR(I)=VSR(I)+VHSR(I)
      EndDo 
C
      Write(6,'(/," SR_xc Energy",X,F15.8)') EnxcSR
      Write(6,'(" SR_H  Energy",X,F15.8)')   EnHSR
      EnSR=EnxcSR+EnHSR 
C
      Return
      End

*Deck PotXCSR_LDA
      Subroutine PotXCSR_LDA(VxcSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NSymMO,NGrid,NInte1,NBasis)
C
C     RETURNS SR XC POTENTIAL IN A MATRIX REPRESENTATION OF THE ORBITALS
C     USED IN OrbGrid
C     VxcSR_ij = Int OrbGrid_i*OrbGri_j*VxcSR[Rho] d3r
C     WHERE Rho IS CONSTRUCTED FROM THE Occ AND URe*OrbGrid
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension VxcSR(NInte1),OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $ Occ(NBasis),URe(NBasis,NBasis),NSymMO(NBasis)
      Dimension OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),
     $ OrbZGrid(NGrid,NBasis) 
C
      Pi2=ASin(One)
      Pi=Two*Pi2
      Const=Three/Four/Pi
C
      Do I=1,NInte1
      VxcSR(I)=Zero
      EndDo
C
      Do I=1,NGrid
C
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      If(Rho.Eq.Zero) Then
      VxcSRi=Zero
      Else
      Rs=(Const/Rho)**(One/Three)
c      Call LSDSR(Rs,Zero,Alpha,EpsxcSR,VxcSRup,VxcSRdown)
      Call LSDSR(Rs,Zero,Alpha,EpsxcSR,EpsxSR,EpscSR,VxcSRup,VxcSRdown)
      VxcSRi=VxcSRup
      EndIf
C
      JK=0
      Do J=1,NBasis
      Do K=1,J
      JK=JK+1
      If(NSymMO(J).Eq.NSymMO(K))
     $ VxcSR(JK)=VxcSR(JK)+OrbGrid(I,J)*OrbGrid(I,K)*VxcSRi*WGrid(I)
      EndDo
      EndDo 
C
      EndDo
C
      Return
      End

*Deck PotHSR
      Subroutine PotHSR(VHSR,Occ,URe,TwoEl,TwoElErf,
     $ NInte1,NInte2,NBasis)
C
C     RETURNS SR HARTREE (COULOMB) POTENTIAL IN AN MO MATRIX REPRESENTATION 
C     WHERE Rho IS CONSTRUCTED FROM THE Occ AND URe*OrbGrid
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension VHSR(NInte1),Occ(NBasis),URe(NBasis,NBasis),
     $ TwoEl(NInte2),TwoElErf(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension Gamma(NInte1)
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
      VHSR(IAB)=Zero
      Do IC=1,NBasis
      Do ID=1,NBasis
      ICD=(Max(IC,ID)*(Max(IC,ID)-1))/2+Min(IC,ID)
      VHSR(IAB)=VHSR(IAB)+Two*Gamma(ICD)*
     $  ( TwoEl(NAddr3(IA,IB,IC,ID))-TwoElErf(NAddr3(IA,IB,IC,ID)) )
      EndDo
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck GetExcSR_LDA
      Subroutine GetExcSR_LDA(EnxcSR,Occ,URe,OrbGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A VALUE OF THE SHORT_RANGE LDA EXC FUNCTIONAL 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc' 
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),Occ(NBasis),
     $ URe(NBasis,NBasis)
      Dimension RhoGrid(NGrid)
C
      Pi2=ASin(One)
      Pi=Two*Pi2
      Const=Three/Four/Pi
C
C     Short-Range LDA functional
C     EnxcSR = Int rho*EpsxcSR d3r
C
      EnxcSR=Zero
      EnxSR=Zero
      EncSR=Zero
      RhoGrid=Zero
C
      Do I=1,NGrid
C
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      If(Rho.Eq.Zero) Then
      EpsxcSR=Zero
      EpsxSR=Zero
      EpscSR=Zero
C   
      Else
C
C     GET SR XC ENERGY DENSITY FOR HEG
C
      Rs=(Const/Rho)**(One/Three)
c     Call LSDSR(Rs,Zero,Alpha,EpsxcSR,VxcSRup,VxcSRdown)
      Call LSDSR(Rs,Zero,Alpha,EpsxcSR,EpsxSR,EpscSR,VxcSRup,VxcSRdown)
C
      EndIf
C
      EnxcSR=EnxcSR+Rho*EpsxcSR*WGrid(I)
      EnxSR=EnxSR+Rho*EpsxSR*WGrid(I)
      EncSR=EncSR+Rho*EpscSR*WGrid(I)
C
      RhoGrid(I) = Rho
C
      EndDo
C
      Print*, 'EnxSR',EnxSR
      Print*, 'EncSR',EncSR
C
C     Test XCFun
      block
      double precision :: XMu
      double precision :: EpscSRGrid(NGrid)
      XMu=0.4d0
      EpscSRGrid=Zero
      Call SRLDAC(RhoGrid,EpscSRGrid,XMu,NGrid)
      EncSR = Zero
      Do I=1,NGrid
      EncSR = EncSR + EpscSRGrid(I)*WGrid(I)
      EndDo
      Print*, 'EncSR-XCFun',EncSR
      end block
C
      Return
      End

*Deck GetExcSR_PBE
      Subroutine GetExcSR_PBE(EnxcSR,VxcSR,Occ,URe,OrbGrid,OrbXGrid,
     $ OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,NInte1,NBasis)
C
C     RETURNS A VALUE OF THE SHORT_RANGE PBE EXC FUNCTIONAL 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc' 
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $ Occ(NBasis),URe(NBasis,NBasis),VxcSR(NInte1)
      Dimension OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),
     $ OrbZGrid(NGrid,NBasis),NSymMO(NBasis)
C
C     Local  
C
      Dimension RhoGrid(NGrid),Sigma(NGrid),Zk(NGrid)
C
! input
      logical fderiv,open
      double precision rhoo(ngrid)
      double precision sigmaco(ngrid),sigmaoo(ngrid)

! output
      integer igrad
      character*(30) name
      double precision vrhoc(ngrid),vrhoo(ngrid)
      double precision vsigmacc(ngrid),vsigmaco(ngrid),vsigmaoo(ngrid)
C
      EnxcSR=Zero
      FDeriv=.True.
      Open=.False.
C
C     Short-Range PBE functional
C
C     COMPUTE THE DENSITY AND ITS GRAD SQUARED ON THE GRID
C
      Do I=1,NGrid
      vrhoc(I)=Zero
      vsigmacc(i)=Zero
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
      EndDo      
C
C     SR-PBE Exchange
C
      If(IFun.Ne.10) Then      
C
      Call dftfun_exerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,rhoo,
     >                   Sigma,sigmaco,sigmaoo,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
      Do I=1,NGrid
      EnxcSR=EnxcSR+Zk(I)*WGrid(I)
      EndDo
C
      EndIf
C
C     SR-PBE Correlation 
C
      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,rhoo,
     >                   Sigma,sigmaco,sigmaoo,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
      EncSR=Zero
      Do I=1,NGrid
      EncSR=EncSR+Zk(I)*WGrid(I)
      EndDo
C
      EnxcSR=EnxcSR+EncSR
C
C     XC POTENTIAL IN THE MO REPRESENTATION
C
      Do I=1,NInte1
      VxcSR(I)=Zero
      EndDo
C
      Do I=1,NGrid
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      JK=0
      Do J=1,NBasis
      Do K=1,J
      JK=JK+1
C  
      If(NSymMO(J).Eq.NSymMO(K)) Then
C
      VxcSR(JK)=VxcSR(JK)+OrbGrid(I,J)*OrbGrid(I,K)*vrhoc(I)*WGrid(I)
C
C     ADD THE NONLOCAL PART (DERIVATIVE W.R.T. THE GRADIENT OF RHO)
C
      VxcSR(JK)=VxcSR(JK)+Two*WGrid(I)*vsigmacc(I)*
     $ (RhoX*(OrbXGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbXGrid(I,K))
     $ +RhoY*(OrbYGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbYGrid(I,K))  
     $ +RhoZ*(OrbZGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbZGrid(I,K)))
C
      EndIf
C
      EndDo
      EndDo
C
      EndDo
C
      Return
      End

*Deck GetExc_PBE
      Subroutine GetExc_PBE(Enxc,Vxc,Occ,URe,OrbGrid,OrbXGrid,
     $ OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,NInte1,NBasis)
C
C     RETURNS A VALUE OF THE SHORT_RANGE PBE EXC FUNCTIONAL 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc' 
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $ Occ(NBasis),URe(NBasis,NBasis),Vxc(NInte1)
      Dimension OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),
     $ OrbZGrid(NGrid,NBasis),NSymMO(NBasis)
C
C     Local  
C
      Dimension RhoGrid(NGrid),Sigma(NGrid),Zk(NGrid)
      Double Precision vrhoc(NGrid),vsigmacc(NGrid)
      Logical FDeriv,Open

      Enxc=Zero
      FDeriv=.True.
      Open=.False.
C
C     PBE functional
C
C     COMPUTE THE DENSITY AND ITS GRAD SQUARED ON THE GRID
C
      Do I=1,NGrid
      vrhoc(I)=Zero
      vsigmacc(i)=Zero
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)

      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
      EndDo      
C
      If(IFun.Ne.10) Then      
      Call dfun_PBE(RhoGrid,Sigma,Zk,vrhoc,vsigmacc,NGrid)
      EndIf

      Do I=1,NGrid
      Enxc = Enxc + Zk(I)*WGrid(I)
      EndDo
C
C     XC POTENTIAL IN THE MO REPRESENTATION
C
      Vxc = Zero
C
      Do I=1,NGrid
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      JK=0
      Do J=1,NBasis
      Do K=1,J
      JK=JK+1
C  
      If(NSymMO(J).Eq.NSymMO(K)) Then
C
      Vxc(JK)=Vxc(JK)+OrbGrid(I,J)*OrbGrid(I,K)*vrhoc(I)*WGrid(I)
C
C     ADD THE NONLOCAL PART (DERIVATIVE W.R.T. THE GRADIENT OF RHO)
C
      Vxc(JK)=Vxc(JK)+Two*WGrid(I)*vsigmacc(I)*
     $ (RhoX*(OrbXGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbXGrid(I,K))
     $ +RhoY*(OrbYGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbYGrid(I,K))
     $ +RhoZ*(OrbZGrid(I,J)*OrbGrid(I,K)+OrbGrid(I,J)*OrbZGrid(I,K)))
C
      EndIf
C
      EndDo
      EndDo
C
      EndDo
C
      Return
      End

*Deck GetKerNO
      Subroutine GetKerNO(XKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NDimKer,NBasis,NGrid)
C
C     RETURNS a SR-KERNEL IN THE NO's REPRESENTATION
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension XKer(NDimKer),Occ(NBasis),URe(NBasis,NBasis),
     $ OrbGrid(NGrid,NBasis),WGrid(NGrid),NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension SRKer(NGrid),RhoVec(NGrid)
C     ,XNOGrid(NBasis,NGrid)
C
C
      Do I=1,NGrid
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      RhoVec(I)=Rho
      EndDo
C
C     PRODUCE A KERNEL ON A GRID FOR DENSITIES RhoVec
C
c      Stop 'Fatal Error in GetKerNO: xcfun NOT AVAILABLE!'
      Call RhoKernel(RhoVec,SRKer,Alpha,NGrid)
C
C     TRANSFORM MO's ON A GRID TO NO's
C
c      Call TrOrbG(XNOGrid,URe,OrbGrid,NGrid,NBasis)
C
      I1234=0
      Do I1=1,NBasis
      Do I2=1,I1
      Do I3=1,I2
      Do I4=1,I3
      I1234=I1234+1
C
      I1I2S=MultpC(NSymNO(I1),NSymNO(I2))
      I3I4S=MultpC(NSymNO(I3),NSymNO(I4))
      ISym=MultpC(I1I2S,I3I4S)
C
      If(ISym.Eq.1) Then      
c      Call XKerEl(XKer1234,I1,I2,I3,I4,SRKer,XNOGrid,WGrid,NGrid,NBasis)
      Call XKerEl(XKer1234,I1,I2,I3,I4,SRKer,OrbGrid,WGrid,NGrid,NBasis)
      XKer(I1234)=XKer1234
C
      Else
      XKer(I1234)=Zero
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck GetKerNPT
      Subroutine GetKerNPT(SRKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NBasis,NGrid)
C
C     RETURNS a SR-KERNEL IN THE NO's REPRESENTATION
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis,NBasis),SRKer(NGrid),
     $ OrbGrid(NGrid,NBasis),WGrid(NGrid),NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension RhoVec(NGrid)
C     ,XNOGrid(NBasis,NGrid)
C
      Do I=1,NGrid
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      RhoVec(I)=Rho
      EndDo
C
C     PRODUCE A KERNEL ON A GRID FOR DENSITIES RhoVec
C
      Call RhoKernel(RhoVec,SRKer,IFunSR,Alpha,NGrid)
C
      Return
      End


*Deck TrOrbG
      Subroutine TrOrbG(XNOGrid,URe,OrbGrid,NGrid,NBasis)
C     
C     OBTAINS NO'S ON A GRID
C    
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
C      Dimension XNOGrid(NBasis,NGrid),OrbGrid(NBasis,NGrid),
C     $ URe(NBasis,NBasis)
      Dimension XNOGrid(NGrid,NBasis),OrbGrid(NGrid,NBasis),
     $ URe(NBasis,NBasis)
C
      Do J=1,NBasis
      Do I=1,NGrid
      XNOGrid(I,J)=Zero
C
      Do IA=1,NBasis
      XNOGrid(I,J)=XNOGrid(I,J)+OrbGrid(I,IA)*URe(J,IA)
      EndDo
      EndDo
      EndDo
C     
      Return
      End

*Deck XKerEl
      Subroutine XKerEl(XKer1234,I1,I2,I3,I4,SRKer,
     $ XNOGrid,WGrid,NGrid,NBasis)
C     
C     RETURNS AN I1,I2,I3,I4 MATRIX ELEMENT OF THE SR KERNEL IN NO
C    
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension SRKer(NGrid),XNOGrid(NGrid,NBasis),WGrid(NGrid)
C
      XKer1234=Zero
C     
      Do I=1,NGrid
      XKer1234=XKer1234+XNOGrid(I1,I)*XNOGrid(I2,I)*XNOGrid(I3,I)
     $ *XNOGrid(I4,I)*SRKer(I)*WGrid(I)
      EndDo
C     
      Return
      End

*Deck NAddrrK
      Integer Function NAddrrK(I1,I2,I3,I4)
C
      Implicit Real*8 (A-H,O-Z)
C
      J1=Max(Max(I1,I2),Max(I3,I4))
      J4=Min(Min(I1,I2),Min(I3,I4))
C
      If (J1.Eq.Max(I1,I2)) Then
      J2=Max(Min(I1,I2),Max(I3,I4))
      Else
      J2=Max(Max(I1,I2),Min(I3,I4))
      EndIf
C
      If(J4.Eq.Min(I1,I2)) Then
      J3=Min(Max(I1,I2),Min(I3,I4))
      Else
      J3=Min(Min(I1,I2),Max(I3,I4))
      EndIf
C
      NAddrrK=((-1 + J1)*J1*(1 + J1)*(2 + J1))/24. 
     $       +((-1 + J2)*J2*(1 + J2))/6. 
     $       +((-1 + J3)*J3)/2. + J4
C
      Return
      End

*Deck SysAv
      Subroutine SysAv(RsInv,Occ,URe,OrbGrid,WGrid,
     $ NGrid,NInte1,NBasis)
C
C     RETURN SYSTEM AVERAGED 1/rs
C     <1/rs> = Int 1/rs(r) Rho(r)  d3r = (4pi/3)^(1/3) Int Rho(r)^(4/3) d3r
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter (Zero=0.0D0,One=1.D0,Two=2.D0,Three=3.0D0,Four=4.0D0)
      Dimension OrbGrid(NGrid,NBasis),WGrid(NGrid),
     $ Occ(NBasis),URe(NBasis,NBasis)
C
      Pi2=ASin(One)
      Pi=Two*Pi2
      Const=One/(Three/Four/Pi)**(One/Three)
      Const1=Three/Four/Pi
C
      RsInv=Zero
C
      Do I=1,NGrid
C
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      RsInv=RsInv+Const*Rho**(Four/Three)*WGrid(I)
C
      EndDo
C
      RsInv=RsInv/(Two*XELE)
C
      Return
      End

*Deck DenGrid
      Subroutine DenGrid(KX,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
C
C     RETURNS A VALUE OF DENSITY AT THE KTH POINT OF THE GRID
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.0D0,Two=2.D0)
      Dimension OrbGrid(NGrid,NBasis),Occ(NBasis),URe(NBasis,NBasis)
C
      Rho=Zero
      Do I=1,NBasis
C
      OrbNO=Zero
      Do J=1,NBasis
      OrbNO=OrbNO+OrbGrid(KX,J)*URe(J,I)
      EndDo
C
      Rho=Rho+Occ(I)*OrbNO**2
      EndDo
      Rho=Two*Rho
C
      Return
      End   

*Deck DenGrad
      Subroutine DenGrad(KX,DRho,Occ,URe,OrbGrid,DOrbGrid,NGrid,NBasis)
C
C     RETURNS A VALUE OF DENSITY GRADIENT AT THE KTH POINT OF THE GRID
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.0D0,Two=2.D0,Four=4.D0)
      Dimension OrbGrid(NGrid,NBasis),Occ(NBasis),URe(NBasis,NBasis),
     $ DOrbGrid(NGrid,NBasis)

C 
      DRho=Zero
      Do I=1,NBasis
C
      OrbNO=Zero
      DOrbNO=Zero
      Do J=1,NBasis
      OrbNO=OrbNO+OrbGrid(KX,J)*URe(J,I)
      DOrbNO=DOrbNO+DOrbGrid(KX,J)*URe(J,I)
      EndDo
C
      DRho=DRho+Occ(I)*OrbNO*DOrbNO
      EndDo
      DRho=Four*DRho
C
      Return
      End

*Deck GetNGrid
      Subroutine GetNGrid(NGrid,Title)
C
      Implicit Real*8 (A-H,O-Z)
      Character*60 Title,FName
      Character*100 Line,Aux,Aux1,Aux2
C
C     GET THE NUMBER OF GRID POINTS FROM THE MOLPRO OUTPUT
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
      Line="grep 'Grid dimensions' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Read(10,*)Line,Aux,I1,I2,I3
      Close(10)
      NGrid=I3
      Write(*,'(/," GRID SIZE",I8,/)')NGrid 
C
      Return
      End

*Deck CheckNBa 
      Subroutine CheckNBa(NBasis,Title)
C
      Implicit Real*8 (A-H,O-Z)
      Character*60 Title,FName
      Character*100 Line,Aux,Aux1,Aux2
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
C     GET THE NUMBER OF ORBITALS FROM THE MOLPRO OUTPUT
C
      Line="grep 'NUMBER OF CONTRACTIONS' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Read(10,*)Aux,Aux1,Aux2,I1
      Close(10)
      If(I1.Ne.NBasis) Then
      Write(*,'(/,"NBasis = ",I6,"  Molpro : ",I6)'),NBasis,I1
      Stop'Fatal Error: NBasis does not check with Molpro'
      EndIf
C
      Return
      End

*Deck GetAlpha
      Subroutine GetAlpha(Title)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Character*60 Title,FName
      Character*100 Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Line
C
C     GET THE VALUE OF THE OMEGA PARAMETER FROM THE MOLPRO OUTPUT
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
      Line="grep 'long-range erf two-' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Read(10,*) Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Alpha
      Close(10)
C
      Write(*,'(/," RANGE PARAMETER ",F8.3,/)')Alpha
C
      Return
      End

*Deck GetGrid
      Subroutine GetGrid(OrbGrid,WGrid,NSymMO,NGrid,Title,NBasis)
C
C     READS VALUES OF MO ORBITALS ON A GRID (FROM title.cube FILES) 
C     AND INTEGRATION WEIGHTS 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName
      Character*100 Line,Aux,Aux1,Aux2
      Dimension OrbGrid(NBasis,NGrid),WGrid(NGrid),NSymMO(NBasis)
C
      Include 'commons.inc'
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
C     GET THE NUMBER OF ORBITALS FROM THE MOLPRO OUTPUT
C
      Line="grep 'NUMBER OF CONTRACTIONS' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Read(10,*)Aux,Aux1,Aux2,I1
      Close(10)
      If(I1.Ne.NBasis) Stop'Fatal Error 0 in Read_Orbs'
C
C     COLLECT NAMES OF FILES WITH ORBITALS IN A tmp.txt FILE
C
      Line="grep '_orbital_' "//FName(1:K+4)//" > tmp.txt" 
      Call System(Line)
C
      Open(10,File='tmp.txt')
C
      ICounter=0
   20 Read(10,'(A60)',End=30)Line(1:60)
      Do I=1,40
      If(Line(I:I+10).Eq."written to ") I0=I+10
      EndDo
      If(I0.Eq.0) Stop'Fatal Error 2 in Read_Orbs'
      ICounter=ICounter+1
C
C     FIND SYMMETRIES OF ORBITALS
C
      Do II=1,65
      If(Line(II:II).Eq."c") NSymMO(ICounter)=
     $ IChar(Line(II-2:II-2))-48
      EndDo
C
C     OPEN A FILE WITH AN ORBITAL ON A GRID AND LOAD IT
C
      If(IFunSR.Ne.0) Then
C  
      Open(15,File=Line(I0:60))
      Read(15,'(1/)')
      Read(15,*)NoAtoms
      NoAtoms=Abs(NoAtoms)
      Read(15,'(1/)')
      Read(15,*)I
      If(I.Ne.NGrid)Stop'Fatal Error 2.5 in Read_Orbs'
      Do J=1,NoAtoms+1
      Read(15,*)I
      EndDo
      Read(15,*) (OrbGrid(ICounter,I),I=1,NGrid)
      Close(15) 
C
      EndIf
C
      GoTo 20
C
   30 Continue
      Close(10)
C      
      If(ICounter.Ne.NBasis) Stop'Fatal Error 3 in Read_Orbs'        
C
      If(IFunSR.Eq.0) Return
C
C     GET INTEGRATION WEIGHTS FROM THE MOLPRO OUT FILE 
C
      Line="grep ' W ' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Open(10,File='tmp.txt')
      Open(20,File='tmp1.txt')
   50 Read(10,'(A80)',End=40)Line(1:80)
      Write(20,*)Line(3:80)
      GoTo 50
   40 Continue
      Close(10)
      Close(20) 
C
      Open(20,File='tmp1.txt')
      Read(20,*)(WGrid(I),I=1,NGrid)
      Close(20)
C
      Call System('rm tmp.txt tmp1.txt')
C
      Return
      End

* Deck GetGrad
      Subroutine GetGrad(OrbXGrid,OrbYGrid,OrbZGrid,NGrid,Title,NBasis)
C     
C     READS VALUES OF THE X,Y,Z, DERIVATIVES OF THE MO ORBITALS ON A GRID 
C
      Implicit Real*8 (A-H,O-Z)
C     
      Character*60 Title,FName
      Character*100 Line
      Dimension OrbXGrid(NBasis,NGrid),OrbYGrid(NBasis,NGrid),
     $ OrbZGrid(NBasis,NGrid)
C     
      Do I=1,60
      FName(I:I)=' '
      EndDo
C     
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
C     COLLECT NAMES OF FILES WITH DERIVATIVES IN A tmp.txt FILE
C     AND READ DATA FROM THEM
C
      Line="grep '_orbgradx_' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Call ReadGrad(OrbXGrid,NBasis,NGrid)
C
      Line="grep '_orbgrady_' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Call ReadGrad(OrbYGrid,NBasis,NGrid)
C
      Line="grep '_orbgradz_' "//FName(1:K+4)//" > tmp.txt"
      Call System(Line)
      Call ReadGrad(OrbZGrid,NBasis,NGrid) 
C
      Return 
      End

*Deck ReadGrad
      Subroutine ReadGrad(OrbGrad,NBasis,NGrid)
C     
C     READS VALUES OF THE DERIVATIVES FROM FILES LISTED IN tmp.txt 
C
      Implicit Real*8 (A-H,O-Z)
C     
      Character*100 Line
      Dimension OrbGrad(NBasis,NGrid)
C
      Open(10,File='tmp.txt')
C
      ICounter=0
   20 Read(10,'(A65)',End=30)Line(1:65)
      Do I=1,40
      If(Line(I:I+10).Eq."written to ") I0=I+10
      EndDo
      If(I0.Eq.0) Stop'Fatal Error 0 in ReadGrad'
      ICounter=ICounter+1
C
C     OPEN A FILE WITH DERIVATIVES ON A GRID AND LOAD IT
C
      Open(15,File=Line(I0:65))
      Read(15,'(1/)')
      Read(15,*)NoAtoms
      NoAtoms=Abs(NoAtoms)
      Read(15,'(1/)')
      Read(15,*)I
      If(I.Ne.NGrid)Stop'Fatal Error 1 ReadGrad'
      Do J=1,NoAtoms+1
      Read(15,*)I
      EndDo
      Read(15,*) (OrbGrad(ICounter,I),I=1,NGrid)
      Close(15)
C
      GoTo 20
C
   30 Continue
      Close(10)
C
      If(ICounter.Ne.NBasis) Stop'Fatal Error 2 in ReadGrad' 
C
      Return
      End

* Deck GetNO 
      Subroutine GetNO(URe,Occ,UMOAO,Title,NBasis)
C
C     EXTRACT NO's FROM THE MOLPRO OUT
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Two=2.0D0)
C
      Character*60 Title,FName
      Character*100 Line,Aux,Aux1,Aux2
      Dimension Occ(NBasis),URe(NBasis,NBasis),UMOAO(NBasis,NBasis)
      Dimension UMOAOInv(NBasis,NBasis),UNOAO(NBasis,NBasis)
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.out'
C
C     GET THE TRANSFORMATION MATRIX FROM AO TO MO
C
      Open(123,File="orb_hf.dat")
      Do I=1,NBasis
      Read(123,*) (UMOAO(I,K),K=1,NBasis)
      EndDo 
      Close(123)
C
C     GET THE TRANSFORMATION MATRIX FROM AO TO NO
C
      Open(123,File="occ_orb_no.dat")
      Read(123,*)Occ
      Do I=1,NBasis
C sometimes very very small occ produced by molpro are negative
      Occ(I)=Abs(Occ(I))
      Read(123,*) (UNOAO(I,K),K=1,NBasis)
      EndDo
      Close(123) 
      Do I=1,NBasis
      Occ(I)=Occ(I)/Two
      EndDo
C
C     INVERT THE UMOAO MATRIX
C
      Call CpyM(UMOAOInv,UMOAO,NBasis)  
      tol = 1.0d-7
      call minvr(UMOAOInv,tol,det,ier,nbasis)
      if (ier.ne.0) then
        Write(6,'(/,'' ERROR : transformation from ao to '',
     $       ''mo basis is singular'')')
        Stop
      endif
C
C     GET THE TRANSFORMATION MATRIX FROM AO TO NO
C
      Do I=1,NBasis
      Do J=1,NBasis
      URe(I,J)=Zero
      Do K=1,NBasis
      URe(I,J)=URe(I,J)+UNOAO(I,K)*UMOAOInv(K,J)
      EndDo
      EndDo
      EndDo
C
      Return
      End

