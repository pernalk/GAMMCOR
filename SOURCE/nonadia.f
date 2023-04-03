*Deck PrimRSMAT
      Subroutine PrimRSMAT(RMAT,SMAT,URe,Occ,XKin,XNuc,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C
C     CONSTRUCT THE MATRICES RMAT AND SMAT USED IN THE NONADIABATIC APPROXIMATION
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension 
     $ SMAT(NDim,NBasis),RMAT(NBasis,NBasis),
     $ URe(NBasis,NBasis),Occ(NBasis),XKin(NInte1),XNuc(NInte1),
     $ TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension Hlp1(NBasis,NInte1)
C
C     Hlp1(p,qa)=Sum_i F_pi <qa|ii>
C
      Do J=1,NInte1
      Do I=1,NBasis
      Hlp1(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      IQA=0
      Do IQ=1,NBasis
      Do IA=1,IQ
      IQA=IQA+1
C
      TwoInt=TwoMO(NAddr3(IQ,I,IA,I))
C
      Do IP=1,NBasis
      Hlp1(IP,IQA)=Hlp1(IP,IQA)+GOCC(Occ(IP),Occ(I),0,IP,I)*TwoInt
      EndDo
C
      EndDo   
      EndDo
      EndDo
C
C     COMPUTE RMAT AND SMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IQ.Ne.IP) IPQ=IPQ+1
C
      Do IA=1,NBasis
      Do IB=IA,IA
C
      XKernel=
     $ (GOCC(Occ(IP),Occ(IA),0,IP,IA)+GOCC(Occ(IQ),Occ(IB),0,IQ,IB))
     $ *(TwoMO(NAddr3(IP,IA,IQ,IB))+TwoMO(NAddr3(IP,IB,IQ,IA)))
C
      IQA=(Max(IA,IQ)*(Max(IA,IQ)-1))/2+Min(IA,IQ)
      IQB=(Max(IB,IQ)*(Max(IB,IQ)-1))/2+Min(IB,IQ)
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      If(IA.Eq.IP) Then
C
      XKernel=XKernel-Hlp1(IB,IQB)
C
      EndIf
C
      If(IB.Eq.IQ) Then
C
      XKernel=XKernel-Hlp1(IA,IPA)
C
      EndIf
C
      If(IA.Eq.IQ) XKernel=XKernel-Hlp1(IQ,IPB)
      If(IB.Eq.IP) XKernel=XKernel-Hlp1(IP,IQA)
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)-Max(IA,IB)+1
C
      XKernel=-XKernel
C
      If(IP.Eq.IQ) Then
      RMAT(IP,IA)=XKernel
      Else
      SMAT(IPQ,IA)=XKernel
      EndIf
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Return
      End

*Deck NonAdiaPolar
      Subroutine NonAdiaPolar
     $ (Omega,ALPHXX,ALPHYY,ALPHZZ,ABPLUS,ABMIN,CMAT,DMAT,D1MAT,
     $ WMAT,RMAT,SMAT,
     $ Occ,URe,DipX,DipY,DipZ,IndN,
     $ NBasis,NInte1,NDimX,NDimN,NELE1)
C
C     DYNAMIC DIPOLE POLARIZABILITY FOR THE PRIMITIVE FUNCTIONALS
C     THE NONADIABATIC APPROXIMAATION (FIRST THE QUANTITIES XI AND ZI ARE DETERMINED, THEN XR AND ZR)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(SMALL=1.D-10,Delta=1.D-4)
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ CMAT(NDimX,NDimN),SMAT(NDimX,NDimN),
     $ D1MAT(NDimN,NDimX),DMAT(NDimN,NDimX),
     $ WMAT(NDimN,NDimN),RMAT(NDimN,NDimN),
     $ ABPROD(NDimX,NDimX),
     $ Temp(NDimN,NDimN),HlpAB(NDimX,NDimX),
     $ Occ(NBasis),URe(NBasis,NBasis),
     $ DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ UX(NDimX+NDimN),UY(NDimX+NDimN),UZ(NDimX+NDimN),
     $ WX(NDimX+NDimN),WY(NDimX+NDimN),
     $ WZ(NDimX+NDimN),
     $ ZRX(NDimN),ZRY(NDimN),ZRZ(NDimN),
     $ BigMAT(NDimX+NDimN,NDimX+NDimN),BigEig(NDimX+NDimN),
     $ Work(NDimX+NDimN,NDimX+NDimN),
     $ BigInv(NDimX+NDimN,NDimX+NDimN),
     $ HlpA(NDimX,NDimX),
     $ Q2(NDimX,NDimN),Q3(NDimN,NDimX),Q4(NDimN,NDimN),
     $ HlpXN(NDimX,NDimN),HlpNX(NDimN,NDimX),
     $ IndN(2,NDimX)
C
      Om2=Omega**2
C
C     COMPUTE THE DIPOLE MOMENT MATRICES IN NO
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      DipNOX(IJ)=Zero
      DipNOY(IJ)=Zero
      DipNOZ(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      DipNOX(IJ)=DipNOX(IJ)+URe(I,IA)*URe(J,IB)*DipX(IAB)
      DipNOY(IJ)=DipNOY(IJ)+URe(I,IA)*URe(J,IB)*DipY(IAB)
      DipNOZ(IJ)=DipNOZ(IJ)+URe(I,IA)*URe(J,IB)*DipZ(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     MULTIPLY A-B BY N^-1
C
      Do K=1,NDimX
      IJ=0
      Do IJ=1,NDimX
      HlpAB(IJ,K)=ABMIN(IJ,K)
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I))
     $ HlpAB(IJ,K)=HlpAB(IJ,K)/(Occ(I)-Occ(J))
      EndDo
      EndDo
C
C     MULTIPLY A+B BY A-B
C
      Call MultpM(ABPROD,ABPLUS,HlpAB,NDimX)
C
C     ADD OMEGA^2 TO DIAGONALS
C
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      ABPROD(IJ,IJ)=ABPROD(IJ,IJ)-Om2*(Occ(I)-Occ(J))
      EndDo
C
C     GET THE PRODUCT CD AND ADD IT TO ABPROD
C
      Call MultpMN(HlpA,CMAT,DMAT,NDimX,NDimN,NDimN,NDimX)
C
      Do I=1,NDimX
      Do J=1,NDimX
      ABPROD(I,J)=ABPROD(I,J)+Two*HlpA(I,J)
      EndDo
      EndDo
C
C     MULTIPLY  SMAT BY N^-1 FROM THE LEFT
C     
      Do K=1,NDimN
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I))
     $ SMAT(IJ,K)=SMAT(IJ,K)/(Occ(I)-Occ(J))
      EndDo
      EndDo
C     
C     Q1 IS STORED IN -ABPROD
C     
C     Q2
C     C*RMAT->Q2
C     
      Call MultpMN(Q2,CMAT,RMAT,NDimX,NDimN,NDimN,NDimN)
C     
C     (A+B)*SMAT -> HlpXN
C     
      Call MultpMN(HlpXN,ABPLUS,SMAT,NDimX,NDimX,NDimX,NDimN)
C
C     Q2+C -> Q2
C
      Do I=1,NDimX
      Do J=1,NDimN
      Q2(I,J)=Q2(I,J)+HlpXN(I,J)
      EndDo
      EndDo
C
C     Q4
C     WMAT*RMAT->Q4
C
      Call MultpMN(Q4,WMAT,RMAT,NDimN,NDimN,NDimN,NDimN)
      Do I=1,NDimN
      Do J=1,NDimN
      Q4(I,J)=-Q4(I,J)
      EndDo
      EndDo
C
C     2*D1*SMAT->Temp
C
      Call MultpMN(Temp,D1MAT,SMAT,NDimN,NDimX,NDimX,NDimN)
      Do I=1,NDimN
      Do J=1,NDimN
      Temp(I,J)=Two*Temp(I,J)
      EndDo
      EndDo
C
C     Temp-Q4 -> Q4
C
      Do I=1,NDimN
      Do J=1,NDimN
      Q4(I,J)=Temp(I,J)-Q4(I,J)
c herer!!!
c a sign of om2 is changed - why???
      If(I.Eq.J) Q4(I,J)=Q4(I,J)-Om2
      EndDo
      EndDo
C
C     Q3
C     D1*N^-1*(A-B) -> Q3
C
      Call MultpMN(Q3,D1MAT,HlpAB,NDimN,NDimX,NDimX,NDimX)
C
C     -WMAT*D -> D1
C
      Call MultpMN(HlpNX,WMAT,DMAT,NDimN,NDimN,NDimN,NDimX)
C
C     -Q3+2D1 -> Q3
C
      Do I=1,NDimN
      Do J=1,NDimX
      Q3(I,J)=-Two*Q3(I,J)-Two*HlpNX(I,J)
      EndDo
      EndDo
C
C     SET UP BigMAT AND INVERT IT
C
      Do I=1,NDimX
      Do J=1,NDimX
      BigMAT(I,J)=-ABPROD(I,J) 
      EndDo
      EndDo
      Do I=1,NDimX
      Do J=1,NDimN
      BigMAT(I,J+NDimX)=Q2(I,J)
      EndDo
      EndDo
      Do I=1,NDimN
      Do J=1,NDimX
      BigMAT(I+NDimX,J)=Q3(I,J)
      EndDo
      EndDo
      Do I=1,NDimN
      Do J=1,NDimN
      BigMAT(I+NDimX,J+NDimX)=Q4(I,J)
      EndDo
      EndDo
C
      NM=NDimX+NDimN
      Call SVDCMP(BigMAT,NM,NM,NM,NM,BigEig,Work)
C
      EigL=BigEig(1)
      Do I=1,NM
      If(Abs(BigEig(I)).Gt.EigL) EigL=Abs(BigEig(I))
      EndDo
C
      Do I=1,NM
      If(Abs(BigEig(I))/EigL.Gt.SMALL) Then
      BigEig(I)=One/BigEig(I)
      Else
      BigEig(I)=Zero
      EndIf
      EndDo
C
      Do I=1,NM
      Do J=1,NM
C
      BigInv(I,J)=Zero
      Do K=1,NM
      BigInv(I,J)=BigInv(I,J)+Work(I,K)*BigEig(K)*BigMAT(J,K)
      EndDo
C
      EndDo
      EndDo
C
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      WX(IJ)=DipNOX(IJ1)*(Occ(I)-Occ(J))
      WY(IJ)=DipNOY(IJ1)*(Occ(I)-Occ(J))
      WZ(IJ)=DipNOZ(IJ1)*(Occ(I)-Occ(J))
      EndDo
C
      Do I=1,NDimN
      Ind=I+NELE1
      II=(Ind*(Ind+1))/2
      WX(I+NDimX)=DipNOX(II)
      WY(I+NDimX)=DipNOY(II)
      WZ(I+NDimX)=DipNOZ(II)
      EndDo
C
C     CALCULATE XI AND ZI AND STORE IT IN UX, UY, UZ
C
      Do I=1,NM
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
C
      Do J=1,NM
      UX(I)=UX(I)+BigInv(I,J)*WX(J)
      UY(I)=UY(I)+BigInv(I,J)*WY(J)
      UZ(I)=UZ(I)+BigInv(I,J)*WZ(J)
      EndDo
C
c      if(i.gt.ndimx)write(*,*)i,uz(i)*omega*omega
      EndDo
C
C     COPY U TO W
C
      Do I=1,NM
      WX(I)=UX(I)
      WY(I)=UY(I)
      WZ(I)=UZ(I)
      EndDo
C
C     OBTAIN UX,UY,UZ
C
      Do I=1,NDimX
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
C
      Do J=1,NDimX
      UX(I)=UX(I)-HlpAB(I,J)*WX(J)
      UY(I)=UY(I)-HlpAB(I,J)*WY(J)
      UZ(I)=UZ(I)-HlpAB(I,J)*WZ(J)
C
      If(J.Le.NDimN) Then
      UX(I)=UX(I)+SMAT(I,J)*WX(NDimX+J)
      UY(I)=UY(I)+SMAT(I,J)*WY(NDimX+J)
      UZ(I)=UZ(I)+SMAT(I,J)*WZ(NDimX+J)
      EndIf
      EndDo
      EndDo
C
C     OBTAIN ZRX,ZRY,ZRZ
C
      Do I=1,NDimN
      ZRX(I)=Zero
      ZRY(I)=Zero
      ZRZ(I)=Zero
C     
      Do J=1,NDimX
      ZRX(I)=ZRX(I)-Two*DMAT(I,J)*WX(J)
      ZRY(I)=ZRY(I)-Two*DMAT(I,J)*WY(J)
      ZRZ(I)=ZRZ(I)-Two*DMAT(I,J)*WZ(J)
      
      If(J.Le.NDimN) Then
      ZRX(I)=ZRX(I)+RMAT(I,J)*WX(NDimX+J)
      ZRY(I)=ZRY(I)+RMAT(I,J)*WY(NDimX+J)
      ZRZ(I)=ZRZ(I)+RMAT(I,J)*WZ(NDimX+J)
      EndIf
      EndDo
      EndDo
C
C
C     CALCULATE THE DYNAMIC DIPOLE POLARIZABILITY
C
      ALPHXX=Zero
      ALPHYY=Zero
      ALPHZZ=Zero
C
      Do I=1,NDimN
      Ind=I+NELE1
      II=(Ind*(Ind+1))/2
      ALPHXX=ALPHXX+Two*ZRX(I)*DipNOX(II)
      ALPHYY=ALPHYY+Two*ZRY(I)*DipNOY(II)
      ALPHZZ=ALPHZZ+Two*ZRZ(I)*DipNOZ(II)
      EndDo
C
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHXX=ALPHXX+Four*(Occ(I)-Occ(J))*DipNOX(IJ1)*UX(IJ)
      ALPHYY=ALPHYY+Four*(Occ(I)-Occ(J))*DipNOY(IJ1)*UY(IJ)
      ALPHZZ=ALPHZZ+Four*(Occ(I)-Occ(J))*DipNOZ(IJ1)*UZ(IJ)
      EndDo
C
      ALPHXX=-ALPHXX
      ALPHYY=-ALPHYY
      ALPHZZ=-ALPHZZ
C
C     RESTORE SMAT
C
      Do K=1,NDimN
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I))
     $ SMAT(IJ,K)=SMAT(IJ,K)*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
      Return
      End








 
      
