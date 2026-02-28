*Deck EXCIT
      Subroutine EXCITA
     $(URe,TwoMO,Occ,XKin,XNuc,DipX,DipY,DipZ,NBasis,NInte1,NInte2,NDim)
C
C     COMPUTE EXCITATION ENERGIES
C 
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(ITMAX=30, Toler=1.D-5)
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),TwoMO(NInte2),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),
     $ DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ ExcEn(0:2)
C
      Do IExact=1,1
C
      Om1=.28D0
      Call POLARIZ2
     $ (IExact,Om1,ALPHXX1,ALPHYY1,ALPHZZ1,URe,TwoMO,Occ,XKin,XNuc,
     $ DipX,DipY,DipZ,NBasis,NInte1,NInte2,NDim)
      F=One/ALPHYY1/ALPHZZ1
C    
      Om2=Om1
      Do I=1,150
      Om2=Om2+1.D-3
      Call POLARIZ2
     $ (IExact,Om2,ALPHXX2,ALPHYY2,ALPHZZ2,URe,TwoMO,Occ,XKin,XNuc,
     $ DipX,DipY,DipZ,NBasis,NInte1,NInte2,NDim)
C   
      FMid=One/ALPHYY2/ALPHZZ2
      If(FMid*F.Lt.Zero) GoTo 888
C 
      EndDo
      Stop 'Fatal error in EXCITA: Wrong range of omega!' 
C
  888 If(F*FMid.Ge.Zero)
     $ Stop 'Fatal error in EXCITA: Wrong range of omega!'
C
      If(F.Lt.Zero) Then
      RtBis=Om1
      DX=Om2-Om1
      Else
      RtBis=Om2
      DX=Om1-Om2
      EndIf
C
      It=Zero
   10 It=It+1
C
      DX=Half*DX
      XMid=RtBis+DX
C
      Call POLARIZ2
     $ (IExact,XMid,ALPHXX,ALPHYY,ALPHZZ,URe,TwoMO,Occ,XKin,XNuc,
     $ DipX,DipY,DipZ,NBasis,NInte1,NInte2,NDim)
C
      FMid=One/ALPHYY/ALPHZZ 
      If(FMid.Le.Zero) RtBis=XMid
      If(Abs(DX).Lt.Toler.Or.FMid.Eq.Zero) GoTo 999 
C
      If(It.Eq.ITMAX) Stop'Fatal error in EXCITA: no convergence!'
      GoTo 10 
C
  999 Continue 
      ExcEn(IExact)=Xmid  
C
      EndDo
C
      Write(*,*)ExcEn(1),ExcEn(2),ExcEn(0)
C
      Return 
      End

*Deck POLARIZ2
      Subroutine POLARIZ2
     $ (IExact,Omega,ALPHXX,ALPHYY,ALPHZZ,URe,TwoMO,Occ,XKin,XNuc,
     $ DipX,DipY,DipZ,NBasis,NInte1,NInte2,NDim)
C
C     DYNAMIC DIPOLE POLARIZABILITY FOR TWO-ELECTRON SYSTEMS (ONLY!)
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
     $ URe(NBasis,NBasis),TwoMO(NInte2),
     $ ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ AUInv(NDim,NDim),HlpA(NDim,NDim),
     $ ABPROD(NDim,NDim),CMAT(NDim,NBasis),
     $ DMAT(NBasis,NDim),D1MAT(NBasis,NDim),
     $ W1MAT(NBasis,NBasis),W2MAT(NDim,NBasis),
     $ W3MAT(NBasis,NBasis),Temp(NBasis,NBasis),
     $ Q2(NDim,NBasis),Q3(NBasis,NDim),Q4(NBasis,NBasis),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),
     $ DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ EigU(NDim),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ UX(NDim+NBasis),UY(NDim+NBasis),UZ(NDim+NBasis),
     $ WX(NDim+NBasis),WY(NDim+NBasis),WZ(NDim+NBasis),
     $ ZRX(NBasis),ZRY(NBasis),ZRZ(NBasis),
     $ BigMAT(NDim+NBasis,NDim+NBasis),BigEig(NDim+NBasis),
     $ Work(NDim+NBasis,NDim+NBasis),BigInv(NDim+NBasis,NDim+NBasis),
     $ CI(NBasis)
C
      Om2=Omega**2
C
      Do I=1,NBasis
      CI(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CI(I)=-CI(I)
      EndDo
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
C     OBTAIN THE MATRICES A,B,C,D AND STORE THEM IN AUInv,HlpA,CMAT,DMAT
C
      If(IExact.Eq.1) Then
      Call GetABCD(AUInv,HlpA,CMAT,DMAT,D1MAT,W1MAT,W2MAT,W3MAT,
     $ URe,Occ,XKin,XNuc,TwoMO,NBasis,NDim,NInte1,NInte2)
      Else
      Call GetABCD(AUInv,HlpA,CMAT,DMAT,D1MAT,W1MAT,W2MAT,W3MAT,
     $ URe,Occ,XKin,XNuc,TwoMO,NBasis,NDim,NInte1,NInte2)
C
      EndIf
C
C     GET A+B AND A-B
C
      Call AddM(AUInv,HlpA,ABPLUS,NDim)
      Call DiffM(AUInv,HlpA,ABMIN,NDim)
C
C     MULTIPLY A-B BY N^-1
C
      Do K=1,NDim
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I)) 
     $ ABMIN(IJ,K)=ABMIN(IJ,K)/(Occ(I)-Occ(J))
      EndDo
      EndDo
      EndDo
C
C     MULTIPLY A+B BY A-B
C     
      Call MultpM(ABPROD,ABPLUS,ABMIN,NDim)
C     
C     ADD OMEGA^2 TO DIAGONALS
C      
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      ABPROD(IJ,IJ)=ABPROD(IJ,IJ)-Om2*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
      If(IExact.Eq.2) GoTo 777
C
C     GET THE PRODUCT CD AND ADD IT TO ABPROD
C
      Call MultpMN(HlpA,CMAT,DMAT,NDim,NBasis,NBasis,NDim) 
C
      Do I=1,NDim
      Do J=1,NDim
      ABPROD(I,J)=ABPROD(I,J)+Two*HlpA(I,J)
      EndDo
      EndDo
C
      If(IExact.Eq.1) GoTo 999
C
C     INVERT ABPROD
C
      Call SVDCMP(ABPROD,NDim,NDim,NDim,NDim,EigU,HlpA)
C
      EigL=EigU(1)
      Do I=1,NDim
      If(Abs(EigU(I)).Gt.EigL) EigL=Abs(EigU(I))
      EndDo
C
      Do I=1,NDim
      If(Abs(EigU(I))/EigL.Gt.SMALL) Then
      EigU(I)=One/EigU(I)
      Else
      EigU(I)=Zero
      EndIf
      EndDo
C
      Do I=1,NDim
      Do J=1,NDim
C
      AUInv(I,J)=Zero
      Do K=1,NDim
      AUInv(I,J)=AUInv(I,J)+HlpA(I,K)*EigU(K)*ABPROD(J,K)
      EndDo
C
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      WX(IJ)=DipNOX(IJ1)*(Occ(I)-Occ(J))
      WY(IJ)=DipNOY(IJ1)*(Occ(I)-Occ(J))
      WZ(IJ)=DipNOZ(IJ1)*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
C     CALCULATE UX, UY, UZ
C
      Do I=1,NDim
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
C
      Do J=1,NDim
      UX(I)=UX(I)+AUInv(I,J)*WX(J)
      UY(I)=UY(I)+AUInv(I,J)*WY(J)
      UZ(I)=UZ(I)+AUInv(I,J)*WZ(J)
      EndDo
      EndDo
C
C     MULTIPLY U BY A-B AND STORE IN W
C
      Do I=1,NDim
      WX(I)=Zero
      WY(I)=Zero
      WZ(I)=Zero
      Do J=1,NDim
      WX(I)=WX(I)+ABMIN(I,J)*UX(J)
      WY(I)=WY(I)+ABMIN(I,J)*UY(J)
      WZ(I)=WZ(I)+ABMIN(I,J)*UZ(J)
      EndDo
      EndDo
C
C     CALCULATE ZR=2 D U
C
      Do I=1,NBasis
      ZRX(I)=Zero
      ZRY(I)=Zero
      ZRZ(I)=Zero
C
      Do J=1,NDim
      ZRX(I)=ZRX(I)+Two*DMAT(I,J)*UX(J)
      ZRY(I)=ZRY(I)+Two*DMAT(I,J)*UY(J)
      ZRZ(I)=ZRZ(I)+Two*DMAT(I,J)*UZ(J)
      EndDo
      EndDo
C
C     COPY W TO U
C
      Do I=1,NDim
      UX(I)=WX(I)
      UY(I)=WY(I)
      UZ(I)=WZ(I)
      EndDo
C
      GoTo 888 
C
  999 Continue
C
C     MULTIPLY  W2 BY N^-1 FROM THE LEFT
C
      Do K=1,NBasis
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I))
     $ W2MAT(IJ,K)=W2MAT(IJ,K)/(Occ(I)-Occ(J))
      EndDo
      EndDo
      EndDo
C
C     Q1 IS STORED IN -ABPROD
C
C     Q2     
C     C*W1->Q2
C
      Call MultpMN(Q2,CMAT,W1MAT,NDim,NBasis,NBasis,NBasis)
C
C     (A+B)*W2 -> C
C
      Call MultpMN(CMAT,ABPLUS,W2MAT,NDim,NDim,NDim,NBasis)
C
C     Q2+C -> Q2
C
      Do I=1,NDim
      Do J=1,NBasis
      Q2(I,J)=Q2(I,J)+CMAT(I,J)
      EndDo
      EndDo
C
C     Q4
C     W3*W1->Q4
C
      Call MultpMN(Q4,W3MAT,W1MAT,NBasis,NBasis,NBasis,NBasis)
C
C     D1*W2->Temp
C
      Call MultpMN(Temp,D1MAT,W2MAT,NBasis,NDim,NDim,NBasis)
C
C     Temp-Q4 -> Q4
C 
      Do I=1,NBasis
      Do J=1,NBasis
      Q4(I,J)=Temp(I,J)-Q4(I,J)
      If(I.Eq.J) Q4(I,J)=Q4(I,J)+Om2*CI(I)
      EndDo
      EndDo
C
C     Q3
C     D1*(A-B) -> Q3
C
      Call MultpMN(Q3,D1MAT,ABMIN,NBasis,NDim,NDim,NDim)
C
C     W3*D -> D1
C
      Call MultpMN(D1MAT,W3MAT,DMAT,NBasis,NBasis,NBasis,NDim)
C
C     -Q3+2D1 -> Q3
C
      Do I=1,NBasis
      Do J=1,NDim
      Q3(I,J)=-Q3(I,J)+Two*D1MAT(I,J)
      EndDo
      EndDo
C
C     SET UP BigMAT AND INVERT IT
C
      Do I=1,NDim 
      Do J=1,NDim
      BigMAT(I,J)=-ABPROD(I,J)
      EndDo
      EndDo
      Do I=1,NDim
      Do J=1,NBasis
      BigMAT(I,J+NDim)=Q2(I,J)
      EndDo
      EndDo
      Do I=1,NBasis
      Do J=1,NDim
      BigMAT(I+NDim,J)=Q3(I,J)
      EndDo
      EndDo
      Do I=1,NBasis
      Do J=1,NBasis
      BigMAT(I+NDim,J+NDim)=Q4(I,J)
      EndDo
      EndDo
C
      NM=NDim+NBasis
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
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      WX(IJ)=DipNOX(IJ1)*(Occ(I)-Occ(J))
      WY(IJ)=DipNOY(IJ1)*(Occ(I)-Occ(J))
      WZ(IJ)=DipNOZ(IJ1)*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
      Do I=1,NBasis
      II=(I*(I+1))/2
      WX(I+NDim)=Two*Occ(I)*DipNOX(II)
      WY(I+NDim)=Two*Occ(I)*DipNOY(II)
      WZ(I+NDim)=Two*Occ(I)*DipNOZ(II)
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
      Do I=1,NDim
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
C
      Do J=1,NDim
      UX(I)=UX(I)-ABMIN(I,J)*WX(J)
      UY(I)=UY(I)-ABMIN(I,J)*WY(J)
      UZ(I)=UZ(I)-ABMIN(I,J)*WZ(J)
C
      If(J.Le.NBasis) Then
      UX(I)=UX(I)+W2MAT(I,J)*WX(NDim+J)
      UY(I)=UY(I)+W2MAT(I,J)*WY(NDim+J)
      UZ(I)=UZ(I)+W2MAT(I,J)*WZ(NDim+J)
      EndIf
      EndDo
      EndDo
C
C     OBTAIN ZRX,ZRY,ZRZ
C
      Do I=1,NBasis
      ZRX(I)=Zero
      ZRY(I)=Zero
      ZRZ(I)=Zero
C
      Do J=1,NDim
      ZRX(I)=ZRX(I)-Two*DMAT(I,J)*WX(J)
      ZRY(I)=ZRY(I)-Two*DMAT(I,J)*WY(J)
      ZRZ(I)=ZRZ(I)-Two*DMAT(I,J)*WZ(J)

      If(J.Le.NBasis) Then
      ZRX(I)=ZRX(I)+W1MAT(I,J)*WX(NDim+J)
      ZRY(I)=ZRY(I)+W1MAT(I,J)*WY(NDim+J)
      ZRZ(I)=ZRZ(I)+W1MAT(I,J)*WZ(NDim+J)
      EndIf
      EndDo
      EndDo
C
      GoTo 888
C
  777 Continue
C
C     ADIABATIC APPROXIMATION WITH A CORRECT ASYMPTOTIC AT Om->0
C
C     Q1 IS STORED IN ABPROD
C
C     Q2
C     C->Q2
C
      Do I=1,NDim
      Do J=1,NBasis
      Q2(I,J)=CMAT(I,J)
      EndDo
      EndDo
C
C     Q4
C     -W3->Q4
C
      Do I=1,NBasis
      Do J=1,NBasis
      Q4(I,J)=-W3MAT(I,J)
      EndDo
      EndDo
C
C     Q3
C     -D1*(A-B) -> Q3
C
      Call MultpMN(Q3,D1MAT,ABMIN,NBasis,NDim,NDim,NDim)
C
C
C     SET UP BigMAT AND INVERT IT
C
      Do I=1,NDim
      Do J=1,NDim
      BigMAT(I,J)=ABPROD(I,J)
      EndDo
      EndDo
      Do I=1,NDim
      Do J=1,NBasis
      BigMAT(I,J+NDim)=Q2(I,J)
      EndDo
      EndDo
      Do I=1,NBasis
      Do J=1,NDim
      BigMAT(I+NDim,J)=Q3(I,J)
      EndDo
      EndDo
      Do I=1,NBasis
      Do J=1,NBasis
      BigMAT(I+NDim,J+NDim)=Q4(I,J)
      EndDo
      EndDo
C
      NM=NDim+NBasis
      Call SVDCMP(BigMAT,NM,NM,NM,NM,BigEig,Work)
C
      EigL=BigEig(1)
      Do I=1,NM
      If(Abs(BigEig(I)).Gt.EigL) EigL=Abs(BigEig(I))
      EndDo
C
      Do I=1,NM
      If(Abs(Work(NDim+1,I)).Gt.0.8) IS=I
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
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      WX(IJ)=DipNOX(IJ1)*(Occ(I)-Occ(J))
      WY(IJ)=DipNOY(IJ1)*(Occ(I)-Occ(J))
      WZ(IJ)=DipNOZ(IJ1)*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
      Do I=1,NBasis
      II=(I*(I+1))/2
      WX(I+NDim)=Two*Occ(I)*DipNOX(II)
      WY(I+NDim)=Two*Occ(I)*DipNOY(II)
      WZ(I+NDim)=Two*Occ(I)*DipNOZ(II)
      EndDo
C
C
C
      If(IFun.Eq.1) Then 
C
C
C
C     CALCULATE MIU1
C
      XMiu1=Zero
      YMiu1=Zero
      ZMiu1=Zero
      Tr=Zero
      Do I=1,NM
      XMiu1=XMiu1+BigMAT(I,IS)*WX(I)
      YMiu1=YMiu1+BigMAT(I,IS)*WY(I)
      ZMiu1=ZMiu1+BigMAT(I,IS)*WZ(I)
      If(I.Gt.NDim) Tr=Tr+BigMAT(I,IS)*Occ(I-NDim)
      EndDo
      XMiu1=XMiu1/Tr
      YMiu1=YMiu1/Tr
      ZMiu1=ZMiu1/Tr
C
C
C
      Else
C
C
C     CALCULATE MIU1
C
      AX=Zero
      AY=Zero
      AZ=Zero
      AM=Zero
C
      Do I=1,NBasis
      Do J=1,NBasis
      AM=AM+BigInv(NDim+I,NDim+J)*Occ(J)
      EndDo
      EndDo
C
      Do I=1,NBasis
      Do J=1,NM
      AX=AX+BigInv(NDim+I,J)*WX(J)
      AY=AY+BigInv(NDim+I,J)*WY(J)
      AZ=AZ+BigInv(NDim+I,J)*WZ(J)
      EndDo
      EndDo

C
      XMiu1=AX/AM
      YMiu1=AY/AM
      ZMiu1=AZ/AM
C
C
C
      EndIf
C
C
C     CALCULATE (A-B)^-1*N*XR AND ZR AND STORE IT IN UX, UY, UZ
C
      Do J=NDim+1,NM
      WX(J)=WX(J)-XMiu1*Occ(J-NDim)
      WY(J)=WY(J)-YMiu1*Occ(J-NDim)
      WZ(J)=WZ(J)-ZMiu1*Occ(J-NDim)
      EndDo
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
      EndDo
C
C     
C     IF IFun=1 THE MIU TAKES CARE OF THE SINGULARITY OF THE MAIN MATRIX
C     AND THE Z VECTOR HAS TO BE NOW NORMALIZED TO ZERO
C
      If(IFun.Eq.1) Then
C
C     NORMALIZE ZR
C
      Sum=Zero 
      Sum1X=Zero
      Sum1Y=Zero
      Sum1Z=Zero
      Do I=NDim+1,NM
      Sum=Sum+Work(I,IS)
      Sum1X=Sum1X+UX(I)
      Sum1Y=Sum1Y+UY(I)
      Sum1Z=Sum1Z+UZ(I)
      EndDo
C
      AlphX=Sum1X/Sum
      AlphY=Sum1Y/Sum
      AlphZ=Sum1Z/Sum
      Do I=1,NM
      UX(I)=UX(I)-AlphX*Work(I,IS)      
      UY(I)=UY(I)-AlphY*Work(I,IS)
      UZ(I)=UZ(I)-AlphZ*Work(I,IS)
      EndDo
C
C
      EndIf
C
C
C     MULTIPLY U BY A-B AND STORE IN W
C
      Do I=1,NDim
      WX(I)=Zero
      WY(I)=Zero
      WZ(I)=Zero
      Do J=1,NDim
      WX(I)=WX(I)+ABMIN(I,J)*UX(J)
      WY(I)=WY(I)+ABMIN(I,J)*UY(J)
      WZ(I)=WZ(I)+ABMIN(I,J)*UZ(J)
      EndDo
      EndDo
C
C     COPY W TO U AND ZR 
C
      Do I=1,NDim
      UX(I)=WX(I)
      UY(I)=WY(I)
      UZ(I)=WZ(I)
      EndDo
C
      Do I=1,NBasis
      ZRX(I)=UX(NDim+I)
      ZRY(I)=UY(NDim+I)
      ZRZ(I)=UZ(NDim+I)
      EndDo
C
  888 Continue 
C
C     CALCULATE THE DYNAMIC DIPOLE POLARIZABILITY
C     
      ALPHXX=Zero
      ALPHYY=Zero
      ALPHZZ=Zero
C     
      IJ=0
      Do I=1,NBasis
C
C     CONTRIBUTION FROM THE OCCUPATIONS
C
      II=(I*(I+1))/2
      ALPHXX=ALPHXX+Two*ZRX(I)*DipNOX(II)
      ALPHYY=ALPHYY+Two*ZRY(I)*DipNOY(II)
      ALPHZZ=ALPHZZ+Two*ZRZ(I)*DipNOZ(II)
C
      Do J=1,I-1
      IJ=IJ+1
C     
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHXX=ALPHXX+Four*(Occ(I)-Occ(J))*DipNOX(IJ1)*UX(IJ)
      ALPHYY=ALPHYY+Four*(Occ(I)-Occ(J))*DipNOY(IJ1)*UY(IJ)
      ALPHZZ=ALPHZZ+Four*(Occ(I)-Occ(J))*DipNOZ(IJ1)*UZ(IJ)
      EndDo
      EndDo
C
      ALPHXX=-ALPHXX
      ALPHYY=-ALPHYY
      ALPHZZ=-ALPHZZ    
C 
c      Write(6,'(1X,4F10.4)')Omega,ALPHYY,ALPHZZ
C
      Return
      End

*Deck DYNPOLAR 
      Subroutine DYNPOLAR 
     $ (Omega,URe,TwoMO,Occ,XKin,XNuc,DipX,DipY,DipZ,
     $ NBasis,NInte1,NInte2,NDim)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(SMALL=1.D-10)
      Parameter(Zero=0.D0,One=1.D0,Two=2.D0,Three=3.D0,Four=4.D0)
C
      Dimension 
     $ URe(NBasis,NBasis),TwoMO(NInte2),
     $ ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ AUInv(NDim,NDim),HlpA(NDim,NDim),
     $ ABPROD(NDim,NDim),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),
     $ DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ EigU(NDim),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ UX(NDim),UY(NDim),UZ(NDim),
     $ WX(NDim),WY(NDim),WZ(NDim)
C
      Om2=Omega**2
C
C     COMPUTE THE DIPOLE MATRICES IN NO
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
C     OBTAIN THE MATRICES AMAT AND BMAT AND STORE THEM IN AUInv,HlpA
C
      Call GetABMAT(AUInv,HlpA,URe,Occ,XKin,XNuc,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C
C     GET A+B AND A-B
C
      Call AddM(AUInv,HlpA,ABPLUS,NDim) 
      Call DiffM(AUInv,HlpA,ABMIN,NDim)
C
C     MULTIPLY A+B BY A-B
C
      Call MultpM(ABPROD,ABPLUS,ABMIN,NDim)
C
C     ADD OMEGA^2 TO DIAGONALS
C      
      Do I=1,NDim
      ABPROD(I,I)=ABPROD(I,I)-Om2
      EndDo
C
C     INVERT ABPROD
C
      Call SVDCMP(ABPROD,NDim,NDim,NDim,NDim,EigU,HlpA)
C
      EigL=EigU(1)
      Do I=1,NDim
      If(Abs(EigU(I)).Gt.EigL) EigL=Abs(EigU(I))
      EndDo
C
      Do I=1,NDim
      If(Abs(EigU(I))/EigL.Gt.SMALL) Then
      EigU(I)=One/EigU(I)
      Else
      EigU(I)=Zero
      EndIf
      EndDo
C
      Do I=1,NDim
      Do J=1,NDim
C
      AUInv(I,J)=Zero
      Do K=1,NDim
      AUInv(I,J)=AUInv(I,J)+HlpA(I,K)*EigU(K)*ABPROD(J,K)
      EndDo
C
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      UX(IJ)=DipNOX(IJ1)*(Occ(I)-Occ(J))
      UY(IJ)=DipNOY(IJ1)*(Occ(I)-Occ(J))
      UZ(IJ)=DipNOZ(IJ1)*(Occ(I)-Occ(J))
      EndDo
      EndDo
C
C     MULTIPLY W BY A+B
C
      Do I=1,NDim
      WX(I)=Zero 
      WY(I)=Zero
      WZ(I)=Zero
      Do J=1,NDim
      WX(I)=WX(I)+ABPLUS(I,J)*UX(J)
      WY(I)=WY(I)+ABPLUS(I,J)*UY(J)
      WZ(I)=WZ(I)+ABPLUS(I,J)*UZ(J)
      EndDo
      EndDo      
C
C     CALCULATE UX, UY, UZ
C
      Do I=1,NDim
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
C
      Do J=1,NDim
      UX(I)=UX(I)+AUInv(I,J)*WX(J)
      UY(I)=UY(I)+AUInv(I,J)*WY(J)
      UZ(I)=UZ(I)+AUInv(I,J)*WZ(J)
      EndDo
      EndDo
C
      ALPHXX=Zero
      ALPHYY=Zero
      ALPHZZ=Zero
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHXX=ALPHXX+(Occ(I)-Occ(J))*DipNOX(IJ1)*UX(IJ)
      ALPHYY=ALPHYY+(Occ(I)-Occ(J))*DipNOY(IJ1)*UY(IJ)
      ALPHZZ=ALPHZZ+(Occ(I)-Occ(J))*DipNOZ(IJ1)*UZ(IJ)
      EndDo
      EndDo
C
      ALPHXX=-ALPHXX*Four
      ALPHYY=-ALPHYY*Four
      ALPHZZ=-ALPHZZ*Four
C
c      Write(6,'(/,4X,'' Omega     Alph_x    Alph_y    Alph_z '')')
      Write(6,'(1X,4F10.4)')Omega,ALPHYY,ALPHZZ
c     $  Omega,ALPHXX,ALPHYY,ALPHZZ
C
      Return
      End      

*Deck CPDMFT
      Subroutine CPDMFT
     $ (URe,TwoMO,Occ,XKin,XNuc,DipX,DipY,DipZ,
     $ NBasis,NInte1,NInte2,NDim)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Tol=1.D-10, SMALL1=1.D-4, SMALL0=1.D-7,
     $          SMALL=1.D-7,TolA=1.D-3, MxIt=20) 
      Parameter(Zero=0.D0,One=1.D0,Two=2.D0,Three=3.D0,Four=4.D0)
C
      Include 'commons.inc'
c herer!!!
      Common/UZ/ UZZ(10000),ZRZZ(100)
C
      Dimension URe(NBasis,NBasis),TwoMO(NInte2),
     $ AMAT(NDim,NDim),AUInv(NDim,NDim),HlpA(NDim,NDim),
     $ BMATX(NDim),BMATY(NDim),BMATZ(NDim),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),
     $ DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ CMAT(NBasis,NBasis,NBasis),HlpMAT(NBasis,NBasis),
     $ ASV(NBasis,NBasis),AInv(NBasis,NBasis),
     $ EigU(NDim),EigN(NBasis),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ UX(NDim),UY(NDim),UZ(NDim),
     $ Occ1X(NBasis),Occ1Y(NBasis),Occ1Z(NBasis)
C
C     NELE  - half of the number of the electrons
C     NELE1 - a number of fully occupied orbitals 
C     NELE0 - a number of unoccupied orbitals 
C
      Sum=Zero
      NELE0=0
      NELE1=0
      Do I=1,NBasis
      Sum=Sum+Occ(I)
      If(One-Occ(I).Le.SMALL1) NELE1=NELE1+1 
      If(Occ(I).Le.SMALL0) NELE0=NELE0+1
      EndDo
      NELE=Sum+0.1D0
      NDim2=NBasis-NELE1-NELE0
C
C     COMPUTE THE DIPOLE MATRICES IN NO
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
C     SET UP INTIAL VALUES OF THE DERIVS OF OCC
C
      Do I=1,NBasis
      Occ1X(I)=Zero
      Occ1Y(I)=Zero
      Occ1Z(I)=Zero
      EndDo
C
C     BEGIN THE ITERATIONS
C
      ALPHXXO=Zero
      ALPHYYO=Zero
      ALPHZZO=Zero
C
      IFlagA=0
      Do It=1,MxIt 
C
      If(It.Gt.1) IFlagA=1
C
C     SET UP MATRICES BMATX, BMATY, BMATZ, AND AMAT
C
      Call ABMAT(IFlagA,AMAT,BMATX,BMATY,BMATZ,URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,Occ1X,Occ1Y,Occ1Z,
     $ NBasis,NDim,NInte1,NInte2)
C
      If(IFlagA.Eq.0) Then
C
C     DECOMPOSE AMAT USING SVD
C
      Call SVDCMP(AMAT,NDim,NDim,NDim,NDim,EigU,HlpA)
C
      EigL=EigU(1)
      Do I=1,NDim
      If(Abs(EigU(I)).Gt.EigL) EigL=Abs(EigU(I))
      EndDo
C
      Do I=1,NDim
      If(Abs(EigU(I))/EigL.Gt.SMALL) Then
      EigU(I)=One/EigU(I)
      Else
      EigU(I)=Zero
      EndIf
      EndDo
C
C     CALCULATE THE INVERSE OF AMAT
C
      Do I=1,NDim
      Do J=1,NDim
C
      AUInv(I,J)=Zero
      Do K=1,NDim
      AUInv(I,J)=AUInv(I,J)+HlpA(I,K)*EigU(K)*AMAT(J,K)
      EndDo
      EndDo
      EndDo
C
      EndIf
C
      Do I=1,NDim
      UX(I)=Zero
      UY(I)=Zero
      UZ(I)=Zero
      EndDo
C
C     CALCULATE UX, UY, UZ
C
      Do I=1,NDim
      Do J=1,NDim
C
      UX(I)=UX(I)+AUInv(I,J)*BMATX(J)
      UY(I)=UY(I)+AUInv(I,J)*BMATY(J)
      UZ(I)=UZ(I)+AUInv(I,J)*BMATZ(J)
      EndDo
      EndDo
C
C     CALCULATE THE CONTRIBUTION TO ALPH POLARIZABILITY
C     FROM THE PERTUBED ORBITALS 
C
      ALPHXXU=Zero
      ALPHYYU=Zero
      ALPHZZU=Zero
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHXXU=ALPHXXU+(Occ(I)-Occ(J))*DipNOX(IJ1)*UX(IJ)
      ALPHYYU=ALPHYYU+(Occ(I)-Occ(J))*DipNOY(IJ1)*UY(IJ)
      ALPHZZU=ALPHZZU+(Occ(I)-Occ(J))*DipNOZ(IJ1)*UZ(IJ)
C
      EndDo
      EndDo
C
      ALPHXXU=ALPHXXU*Four
      ALPHYYU=ALPHYYU*Four
      ALPHZZU=ALPHZZU*Four
C
      Write(6,'(/,1X,'' ITERATION '',I3)') It      
      Write(6,'(1X,'' ALPHU  '',3F10.4)') ALPHXXU,ALPHYYU,ALPHZZU 
C
C     SET UP NEW MATRICES BMATX, BMATY, BMATZ, AND HlpMAT NEEDED FOR 
C     CALCULATION OF THE DERIVATIVES OF THE OCCUPATION NUMBERS
C
      Call ABMATN(IFlagA,BMATX,BMATY,BMATZ,HlpMAT,
     $ UX,UY,UZ,URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,NBasis,NInte1,NInte2,NDim)
C
      If(IFlagA.Eq.0) Then
C
      Do I=1,NDim2
      Do J=1,NDim2
      ASV(I,J)=HlpMAT(I+NELE1,J+NELE1)
      EndDo
      EndDo
C
      Call SVDCMP(ASV,NDim2,NDim2,NBasis,NBasis,EigN,HlpMAT)
C
      EigS=Abs(EigN(1))
      EigL=EigN(1)
      Do I=1,NDim2
      If(EigN(I).Gt.EigL) EigL=EigN(I)
      If(Abs(EigN(I)).Lt.EigS) Then
      IS=I
      EigS=Abs(EigN(I))
      EndIf
      EndDo
C
      Do I=1,NDim2
      If(EigN(I)/EigL.Gt.SMALL) Then
      EigN(I)=One/EigN(I)
      Else
      EigN(I)=Zero
      EndIf
      EndDo
C
C     CALCULATE THE INVERSE OF ASV
C
      Do I=1,NDim2
      Do J=1,NDim2
      AInv(I,J)=Zero    
      Do K=1,NDim2
      AInv(I,J)=AInv(I,J)+HlpMAT(I,K)*EigN(K)*ASV(J,K)
      EndDo
      EndDo
      EndDo     
C
      EndIf 
C
C     IF THE FUNCTIONAL IFun=1 IS USED THEN THE MATRIX ASV IS SINGULAR AND
C     IT REQUIRES A SPECIAL TREATMENT
C
C
C
      If(IFun.Eq.1) Then
C
C
C
C     CALCULATE MIU1
C
      XMiu1=Zero 
      YMiu1=Zero
      ZMiu1=Zero
      Tr=Zero
      Do I=1,NDim2
      XMiu1=XMiu1-ASV(I,IS)*BMATX(I+NELE1)
      YMiu1=YMiu1-ASV(I,IS)*BMATY(I+NELE1)
      ZMiu1=ZMiu1-ASV(I,IS)*BMATZ(I+NELE1) 
      Tr=Tr+ASV(I,IS)
      EndDo
      XMiu1=XMiu1/Tr
      YMiu1=YMiu1/Tr
      ZMiu1=ZMiu1/Tr
C
C     CALCULATE THE DERIVATIVES OF THE OCC
C
      Do I=1,NBasis
      Occ1X(I)=Zero
      Occ1Y(I)=Zero
      Occ1Z(I)=Zero
      EndDo
C
      Do I=1,NDim2
      Do J=1,NDim2
      Occ1X(I+NELE1)=Occ1X(I+NELE1)+AInv(I,J)*(XMiu1+BMATX(J+NELE1))
      Occ1Y(I+NELE1)=Occ1Y(I+NELE1)+AInv(I,J)*(YMiu1+BMATY(J+NELE1))
      Occ1Z(I+NELE1)=Occ1Z(I+NELE1)+AInv(I,J)*(ZMiu1+BMATZ(J+NELE1))
      EndDo
      EndDo
C
C     ADD A VECTOR CORRESPONDING TO THE ZERO EIGENFUNCTION
C     WITH A PROPER COEFFICIENT TO OBTAIN CORRECT NORM OF Occ1
C 
      SX=Zero
      SY=Zero
      SZ=Zero
      Do I=1,NDim2
      SX=SX+Occ1X(I+NELE1)
      SY=SY+Occ1Y(I+NELE1)
      SZ=SZ+Occ1Z(I+NELE1)
      EndDo
      XLamX=-SX/Tr
      XLamY=-SY/Tr
      XLamZ=-SZ/Tr
C
      Do I=1,NDim2
      Occ1X(I+NELE1)=Occ1X(I+NELE1)+XLamX*ASV(I,IS)
      Occ1Y(I+NELE1)=Occ1Y(I+NELE1)+XLamY*ASV(I,IS)
      Occ1Z(I+NELE1)=Occ1Z(I+NELE1)+XLamZ*ASV(I,IS)
      EndDo
C
C
C
      Else
C
C
C
C
C     CALCULATE MIU1
C
      AX=Zero
      AY=Zero
      AZ=Zero
      AM=Zero
C
      Do I=1,NDim2
      Do J=1,NDim2
      AX=AX-AInv(I,J)*BMATX(J+NELE1)
      AY=AY-AInv(I,J)*BMATY(J+NELE1)
      AZ=AZ-AInv(I,J)*BMATZ(J+NELE1)
      AM=AM+AInv(I,J)
      EndDo
      EndDo
      XMiu1=AX/AM
      YMiu1=AY/AM
      ZMiu1=AZ/AM
C
C     CALCULATE THE DERIVATIVES OF THE OCC
C
      Do I=1,NBasis
      Occ1X(I)=Zero
      Occ1Y(I)=Zero
      Occ1Z(I)=Zero
      EndDo
C
      Do I=1,NDim2
      Do J=1,NDim2
      Occ1X(I+NELE1)=Occ1X(I+NELE1)+AInv(I,J)*(XMiu1+BMATX(J+NELE1))
      Occ1Y(I+NELE1)=Occ1Y(I+NELE1)+AInv(I,J)*(YMiu1+BMATY(J+NELE1))
      Occ1Z(I+NELE1)=Occ1Z(I+NELE1)+AInv(I,J)*(ZMiu1+BMATZ(J+NELE1))
      EndDo
      EndDo
C
C
      EndIf
C
C     CALCULATE THE CONTRIBUTION TO ALPH FROM
C     THE PERTURBED OCCUPANCIES
C
      ALPHXXN=Zero
      ALPHYYN=Zero
      ALPHZZN=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      ALPHXXN=ALPHXXN+Two*OCC1X(I)*DipNOX(II)
      ALPHYYN=ALPHYYN+Two*OCC1Y(I)*DipNOY(II)
      ALPHZZN=ALPHZZN+Two*OCC1Z(I)*DipNOZ(II)
      EndDo
c
      Write(6,'(1X,'' ALPHN  '',3F10.4)') ALPHXXN,ALPHYYN,ALPHZZN
C
      ALPHXX=ALPHXXU+ALPHXXN
      ALPHYY=ALPHYYU+ALPHYYN
      ALPHZZ=ALPHZZU+ALPHZZN
C
      Write(6,'(1X,'' ALPH   '',3F10.4)') ALPHXX,ALPHYY,ALPHZZ 
C
C     CHECK THE CONVERGENCE 
C
      ErrX= Abs((ALPHXX-ALPHXXO)/ALPHXX)
      ErrY= Abs((ALPHYY-ALPHYYO)/ALPHYY)
      ErrZ= Abs((ALPHZZ-ALPHZZO)/ALPHZZ)
C
      Write(6,'(1X,'' ERRORS '',3E10.1)') ErrX,ErrY,ErrZ 
C
      If(ErrX.Lt.TolA.And.ErrY.Lt.TolA.And.ErrZ.Lt.TolA) Goto 999
      ALPHXXO=ALPHXX
      ALPHYYO=ALPHYY
      ALPHZZO=ALPHZZ
C 
      EndDo
C
  999 Write(6,'(/,X,'' *** CONVERGENCE ACHIEVED *** '',/)') 
C
C     CALCULATE THE ALPH_XY POLARIZABILITY
C
      ALPHXY=Zero
C
      IJ=0
      Do I=1,NBasis
      II=(I*(I+1))/2
      ALPHXY=ALPHXY+Two*OCC1Y(I)*DipNOX(II)
      Do J=1,I-1
      IJ=IJ+1
C
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHXY=ALPHXY+Four*(Occ(I)-Occ(J))*DipNOX(IJ1)*UY(IJ)
      EndDo
      EndDo
C
C     CALCULATE THE ALPH_YX POLARIZABILITY
C
      ALPHYX=Zero
C
      IJ=0
      Do I=1,NBasis
      II=(I*(I+1))/2
      ALPHYX=ALPHYX+Two*OCC1X(I)*DipNOY(II)
      Do J=1,I-1
      IJ=IJ+1
C
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
      ALPHYX=ALPHYX+Four*(Occ(I)-Occ(J))*DipNOY(IJ1)*UX(IJ)
      EndDo
      EndDo
C
      Write(6,'(1X,'' ALPHXY ALPHYX  '',2F10.4)') ALPHXY,ALPHYX
C
C     FINAL CHECKING 
C
      Call ABMAT(0,AMAT,BMATX,BMATY,BMATZ,URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,Occ1X,Occ1Y,Occ1Z,
     $ NBasis,NDim,NInte1,NInte2)
C
      ErrUX=Zero
      ErrUY=Zero
      ErrUZ=Zero 
      Do I=1,NDim
      SumX=Zero
      SumY=Zero
      SumZ=Zero
      Do J=1,NDim
      SumX=SumX+AMAT(I,J)*UX(J)
      SumY=SumY+AMAT(I,J)*UY(J)
      SumZ=SumZ+AMAT(I,J)*UZ(J) 
      EndDo
      ErrUX=ErrUX+(SumX-BMATX(I))**2
      ErrUY=ErrUY+(SumY-BMATY(I))**2
      ErrUZ=ErrUZ+(SumZ-BMATZ(I))**2
      EndDo
      ErrUX=SQRT(ErrUX)
      ErrUY=SQRT(ErrUY)
      ErrUZ=SQRT(ErrUZ)
C
      Write(6,'(/,1X,'' ERRORS OF THE U-CPDMFT EQUATIONS'',3E10.1)')
     $ ErrUX,ErrUY,ErrUZ
C
      Call ABMATN(0,BMATX,BMATY,BMATZ,HlpMAT,
     $ UX,UY,UZ,URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,NBasis,NInte1,NInte2,NDim)
C
      ErrNX=Zero
      ErrNY=Zero
      ErrNZ=Zero
      Do I=1+NELE1,NELE1+NDim2
      SumX=Zero
      SumY=Zero
      SumZ=Zero
      Do J=1+NELE1,NELE1+NDim2
      SumX=SumX+HlpMAT(I,J)*Occ1X(J)
      SumY=SumY+HlpMAT(I,J)*Occ1Y(J)
      SumZ=SumZ+HlpMAT(I,J)*Occ1Z(J)
      EndDo
      ErrNX=ErrNX+(SumX-XMiu1-BMATX(I))**2
      ErrNY=ErrNY+(SumY-YMiu1-BMATY(I))**2
      ErrNZ=ErrNZ+(SumZ-ZMiu1-BMATZ(I))**2
      EndDo
      ErrNX=SQRT(ErrNX)
      ErrNY=SQRT(ErrNY)
      ErrNZ=SQRT(ErrNZ)
C
      Write(6,'(/,1X,'' ERRORS OF THE N-CPDMFT EQUATIONS'',3E10.1)')
     $ ErrNX,ErrNY,ErrNZ
C
c herer!!!
      do i=1,ndim
      uzz(i)=uz(i)
      enddo
      do i=1,nbasis
      zrzz(i)=occ1z(i)
      enddo
C
      Return
      End

*Deck GetABCD
      Subroutine GetABCD(AMAT,BMAT,CMAT,DMAT,D1MAT,W1MAT,W2MAT,W3MAT,
     $ URe,Occ,XKin,XNuc,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE MATICES A,B,C,D
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ CMAT(NDim,NBasis),DMAT(NBasis,NDim),
     $ D1MAT(NBasis,NDim),
     $ W1MAT(NBasis,NBasis),W2MAT(NDim,NBasis),
     $ W3MAT(NBasis,NBasis),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension Hlp2(NInte1),CI(NBasis),HNO(NInte1)
C
      Do I=1,NBasis
      CI(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CI(I)=-CI(I)
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
      EndDo
      EndDo
      EndDo
      EndDo 
C
C     Hlp2(a,b)=Sum_d c_d <ab|dd>
C      
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Hlp2(IJ)=Zero
      Do K=1,NBasis
      Hlp2(IJ)=Hlp2(IJ)+CI(K)*TwoMO(NAddr3(I,K,J,K))  
      EndDo
      EndDo
      EndDo
C
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IQ.Ne.IP) IPQ=IPQ+1
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IQA=(Max(IA,IQ)*(Max(IA,IQ)-1))/2+Min(IA,IQ)
      IQB=(Max(IB,IQ)*(Max(IB,IQ)-1))/2+Min(IB,IQ)
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      Hlp1=Zero
      If(IA.Eq.IP) Hlp1=-HNO(IQB)      
      If(IB.Eq.IQ) Hlp1=Hlp1+HNO(IPA)
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)-Max(IA,IB)+1
C
C
      If(IA.Eq.IB.And.IP.Ne.IQ) Then
C
C
      HlpPQA=-Half*(CI(IP)-CI(IQ))/CI(IA)
     $ *TwoMO(NAddr3(IP,IA,IQ,IA))
C
      If(IA.Eq.IQ) Then
      HlpPQA=HlpPQA-Half*
     $ HNO((Max(IP,IQ)*(Max(IP,IQ)-1))/2+Min(IP,IQ))*
     $ (CI(IP)+CI(IQ))/CI(IA)
      ElseIf(IA.Eq.IP) Then
      HlpPQA=HlpPQA+Half*
     $ HNO((Max(IP,IQ)*(Max(IP,IQ)-1))/2+Min(IP,IQ))*
     $ (CI(IP)+CI(IQ))/CI(IA)
      EndIf
C
      CMAT(IPQ,IA)=HlpPQA+Hlp1
C
C
      EndIf
C
C
      Term=(CI(IP)*CI(IA)+CI(IQ)*CI(IB))*
     $ (TwoMO(NAddr3(IP,IB,IQ,IA))+TwoMO(NAddr3(IP,IA,IQ,IB)))
      If(IQ.Eq.IB) Term=Term-CI(IA)*Hlp2(IPA)
      If(IQ.Eq.IA) Term=Term-CI(IA)*Hlp2(IPB)
      If(IP.Eq.IB) Term=Term-CI(IB)*Hlp2(IQA)
      If(IP.Eq.IA) Term=Term-CI(IB)*Hlp2(IQB)
C
      If(IA.Eq.IB) Then
C
C
      If(IP.Gt.IQ) Then
      W2MAT(IPQ,IA)=Term*Half/CI(IA)
      Else
      W1MAT(IP,IA)=Term*Half/CI(IA)
      EndIf 
C
C
      Else 
C
C
      Term=Hlp1*(Occ(IA)-Occ(IB))-Term
C
      If(IP.Eq.IQ) Then
      DMAT(IP,IAB)=Term
      ElseIf(IA.Gt.IB) Then 
      AMAT(IPQ,IAB)=Term
      Else
      BMAT(IPQ,IAB)=-Term
      EndIf 
C
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     OBTAIN D1MAT AND W3MAT 
C
      Do IP=1,NBasis
      IAB=0
      Do IA=1,NBasis
C
      W3MAT(IP,IA)=W1MAT(IP,IA)/CI(IA)/Four
C
      Do IB=1,IA-1
      IAB=IAB+1
C
      D1MAT(IP,IAB)=DMAT(IP,IAB)*(CI(IA)-CI(IB))/(CI(IA)+CI(IB))
C
      EndDo
      EndDo
      EndDo
C
      Return
      End   

*Deck GetABMAT
      Subroutine GetABMAT(AMAT,BMAT,
     $ URe,Occ,XKin,XNuc,TwoMO,
     $ NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE MATRICES B AND A NEEDED FOR CP-DMFT
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),TwoMO(NInte2),
     $ CMAT(NBasis,NBasis,NBasis)
C
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
C     SET UP AN AUXILIARY MATRIX CMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
C
      HIJ=Zero
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HIJ=HIJ+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
      EndDo
      EndDo
C
      Do K=1,NBasis
C
      CMAT(I,J,K)=HIJ*(Occ(I)-Occ(K))
C
      Do L=1,NBasis
C
      Fikl=(FunPair(Occ(I),Occ(L),ALPH,I,L)
     $ -FunPair(Occ(K),Occ(L),ALPH,K,L))
C
      CMAT(I,J,K)=CMAT(I,J,K)
     $ +XKU*Two*(Occ(I)-Occ(K))*Occ(L)*TwoMO(NAddr3(I,J,L,L))
     $ +Fikl*TwoMO(NAddr3(I,L,L,J))
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     SET UP AN AMAT MATRIX
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      KL=0
      Do K=1,NBasis
C
      Fijk=(FunPair(Occ(I),Occ(K),ALPH,I,K)
     $ -FunPair(Occ(J),Occ(K),ALPH,J,K))
C
      Do L=1,K-1
      KL=KL+1
C
      Fijl=(FunPair(Occ(I),Occ(L),ALPH,I,L)
     $ -FunPair(Occ(J),Occ(L),ALPH,J,L))
C
      AMAT(IJ,KL)=-XKU*Two*(Occ(I)-Occ(J))*(Occ(K)-Occ(L))*
     $ TwoMO(NAddr3(I,J,L,K))
     $ -(Fijk-Fijl)*(TwoMO(NAddr3(I,K,L,J)))
C
      BMAT(IJ,KL)=XKU*Two*(Occ(I)-Occ(J))*(Occ(K)-Occ(L))*
     $ TwoMO(NAddr3(I,J,K,L))
     $ +(Fijk-Fijl)*TwoMO(NAddr3(I,L,K,J))
C
      If(L.Eq.J) AMAT(IJ,KL)=AMAT(IJ,KL)+CMAT(I,K,L)
      If(I.Eq.K) AMAT(IJ,KL)=AMAT(IJ,KL)+CMAT(J,L,K)
      If(K.Eq.J) BMAT(IJ,KL)=BMAT(IJ,KL)+CMAT(I,L,K)
      If(I.Eq.L) BMAT(IJ,KL)=BMAT(IJ,KL)+CMAT(J,K,L)
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Return
      End

*Deck GetCMAT
      Subroutine GetCMAT(CMAT,Occ,TwoMO,NBasis,NDim,NInte2)
C     
C     COMPUTE THE MATRIX C
C     
      Implicit Real*8 (A-H,O-Z)
C        
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C        
      Include 'commons.inc'
C     
      Dimension CMAT(NDim,NBasis),Occ(NBasis),TwoMO(NInte2)
C
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
      Do M=1,NBasis
C
      IJ=0
      Do I=1,NBasis
      II=(I*(I+1))/2
C
      Do J=1,I-1
      JJ=(J*(J+1))/2
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
C
      CMAT(IJ,M)=Zero
C
      ExInt=TwoMO(NAddr3(I,M,M,J))
C
      Term1=XKU*Two*(Occ(I)-Occ(J))*TwoMO(NAddr3(I,J,M,M))
     $ +ExInt*(DerFun(Occ(M),Occ(I),1,M,I)-DerFun(Occ(M),Occ(J),1,M,J))
      Term1=Term1*(Occ(I)-Occ(J))
C
      CMAT(IJ,M)=CMAT(IJ,M)+Term1
C
      If(I.Eq.M.Or.J.Eq.M) Then
C
      Do K=1,NBasis
C
      ExInt=TwoMO(NAddr3(I,K,K,J))
C
      If(I.Eq.M) Then
C
      DFik=DerFun(Occ(I),Occ(K),1,I,K)*(Occ(I)-Occ(J))
      CMAT(IJ,M)=CMAT(IJ,M)-ExInt*
     $(-DFik+FunPair(Occ(I),Occ(K),ALPH,I,K)
     $     -FunPair(Occ(J),Occ(K),ALPH,J,K))
C
      Else
      DFjk=DerFun(Occ(J),Occ(K),1,J,K)*(Occ(I)-Occ(J))
      CMAT(IJ,M)=CMAT(IJ,M)-ExInt*
     $ (DFjk-FunPair(Occ(I),Occ(K),ALPH,I,K)
     $     +FunPair(Occ(J),Occ(K),ALPH,J,K))
C
      EndIf
C     
      EndDo
C
      EndIf
C
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck DerFun
      Real*8 Function DerFun(X,Y,IFlag,I,J)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Small=1.D-6) 
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Three=3.D0,
     $  Four=4.D0)
C
      Include 'commons.inc'
C
C     IFlag =	1 - calculate dF(x,y)/dx
C               2 - calculate d^2 F(x,y)/d^2 x	
C               3 - calculate d^2 F(x,y)/dx/dy
C
      If(X.Lt.Small.Or.Y.Lt.Small) Then
      DerFun=Zero 
      Return 
      EndIf
C
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
      If(IFlag.Eq.1) Then
      DerFun=Half*FacX*FacY*SQRT(Y/X)
C
      ElseIf(IFlag.Eq.2) Then
      DerFun=-FacX*FacY/Four*SQRT(Y/X)/X
C
      Else
      DerFun=FacX*FacY/Four/SQRT(X*Y)
C
      EndIf
C
C     ** End of KU **
      Return
      EndIf
C
C
C     ** BB **
C
      If(IFun.Eq.2) Then
C
      If(IFlag.Eq.1) Then       
      DerFun=-Half*SQRT(Y/X)
C
      ElseIf(IFlag.Eq.2) Then      
      DerFun=One/Four*SQRT(Y/X)/X     
C
      Else
      DerFun=-One/Four/SQRT(X*Y)
C
      EndIf
C
C     ** End of BB **
C
      Return
      EndIf
C
C     ** GU **
C
      If(IFun.Eq.3) Then
C
      If(IFlag.Eq.1) Then   
      DerFun=-Half*SQRT(Y/X)
      If(I.Eq.J) DerFun=-X
C     
      ElseIf(IFlag.Eq.2) Then
      DerFun=One/Four*SQRT(Y/X)/X
      If(I.Eq.J) DerFun=-Half
C
      Else
      DerFun=-One/Four/SQRT(X*Y)
      If(I.Eq.J) DerFun=-Half
C
      EndIf
C
C     ** End of GU **
C
      Return
      EndIf
C     
C     ** BBC1 **
C
      If(IFun.Eq.4) Then
C
      If(IFlag.Eq.1) Then
      DerFun=-Half*SQRT(Y/X)
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE)
     $ DFunPair=Half*SQRT(Y/X)
C
      ElseIf(IFlag.Eq.2) Then
      DerFun=One/Four*SQRT(Y/X)/X
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE) 
     $ DerFun=-One/Four*SQRT(Y/X)/X
C
      Else
      DerFun=-One/Four/SQRT(X*Y)
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE)
     $ DerFun=One/Four/SQRT(X*Y)
C
      EndIf
C
C     ** End of BBC1 **
C
      Return
      EndIf
C
C     ** BBC2 OR BBC3 **
C
      If(IFun.Eq.5.Or.IFun.Eq.6) Then
C
C
C     FIRST X DERIVATIVE
C
C
      If(IFlag.Eq.1) Then
C
      DerFun=-Half*SQRT(Y/X)
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE)
     $ DerFun=Half*SQRT(Y/X)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Le.NELE.And.J.Le.NELE.And.IType(I)*IType(J).Eq.1)
     $ DerFun=-Y
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(I.Le.NELE.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ DerFun=-Y
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(J.Le.NELE.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ DerFun=-Y
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) DerFun=-X
C
C
C     SECOND XX DERIVATIVE     
C
C
      ElseIf(IFlag.Eq.2) Then
C
      DerFun=One/Four*SQRT(Y/X)/X
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE)
     $ DerFun=-One/Four*SQRT(Y/X)/X
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Le.NELE.And.J.Le.NELE.And.IType(I)*IType(J).Eq.1)
     $ DerFun=Zero
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(I.Le.NELE.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ DerFun=Zero
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(J.Le.NELE.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ DerFun=Zero
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) DerFun=-Half
C
C     MIXED XY DERIVATIVE
C
      Else
C
      DerFun=-One/Four/SQRT(X*Y)
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE)
     $ DerFun=One/Four/SQRT(X*Y)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Le.NELE.And.J.Le.NELE.And.IType(I)*IType(J).Eq.1)
     $ DerFun=-One
C
      If(IFun.Eq.5) Return
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If(I.Le.NELE.And.IType(I).Eq.1.And.IType(J).Eq.0)
     $ DerFun=-One
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If(J.Le.NELE.And.IType(J).Eq.1.And.IType(I).Eq.0)
     $ DerFun=-One
C
C     ORBITALS I AND J ARE BOTH STRONLGY OR WEAKLY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.IType(I).Eq.1) DerFun=-Half
C
      EndIf
C
      Return
      EndIf
C
      Return
      End

*Deck ABMATN
      Subroutine ABMATN(IFlag,BMATX,BMATY,BMATZ,HlpMAT,
     $ UX,UY,UZ,
     $ URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,NBasis,NInte1,NInte2,NDim)
C
C     COMPUTE THE MATRICES B0 AND A NEEDED FOR CP-DMFT
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter (Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension BMATX(NBasis),BMATY(NBasis),BMATZ(NBasis), 
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),TwoMO(NInte2),
     $ CMAT(NBasis,NBasis,NBasis),HlpMAT(NBasis,NBasis),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ UX(NDim),UY(NDim),UZ(NDim)
C
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
C     SET UP AN AUXILIARY MATRIX CMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
C
      HIJ=Zero
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HIJ=HIJ+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
      EndDo
      EndDo
C
      CMAT(I,J,I)=HIJ
C
      Do L=1,NBasis
C
      Fil=DerFun(Occ(I),Occ(L),1,I,L)
C
      CMAT(I,J,I)=CMAT(I,J,I)+XKU*Two*Occ(L)*TwoMO(NAddr3(I,J,L,L))
     $           +Fil*TwoMO(NAddr3(I,L,L,J))
C
      EndDo
C
      EndDo
      EndDo
C
      If(IFlag.Eq.0) Then
C
C     SET UP A MATRIX A STORED IN HlpMAT 
C
      Do I=1,NBasis
      Do J=1,NBasis
      HlpMAT(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      BMATX(I)=Zero
C      
      Do J=1,NBasis
C
      Hlp=TwoMO(NAddr3(I,J,I,J))
      BMATX(I)=BMATX(I)+Hlp*DerFun(Occ(I),Occ(J),2,I,J) 
      HlpMAT(I,J)=XKU*Two*TwoMO(NAddr3(I,I,J,J))
     $ +Hlp*DerFun(Occ(I),Occ(J),3,I,J)
C      
      EndDo
      EndDo
C
      Do I=1,NBasis
      HlpMAT(I,I)=HlpMAT(I,I)+BMATX(I)       
      EndDo
C
      EndIf
C
C     CONSTRUCT MATRICES BMATX,BMATY,BMATZ 
C
      Do I=1,NBasis
      II=(I*(I+1))/2
      BMATX(I)=-DipNOX(II)
      BMATY(I)=-DipNOY(II)
      BMATZ(I)=-DipNOZ(II) 
C
      JK=0
      Do J=1,NBasis
C
      If(I.Ne.J) Then
C
      IJ=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)-Max(I,J)+1
      FacU=UX(IJ)
      If(J.Gt.I) FacU=-FacU
      BMATX(I)=BMATX(I)-Two*FacU*CMAT(I,J,I) 
      FacU=UY(IJ)
      If(J.Gt.I) FacU=-FacU
      BMATY(I)=BMATY(I)-Two*FacU*CMAT(I,J,I)
      FacU=UZ(IJ)
      If(J.Gt.I) FacU=-FacU
      BMATZ(I)=BMATZ(I)-Two*FacU*CMAT(I,J,I)
C
      EndIf
C
      Do K=1,J-1
C
      JK=JK+1
      Fijk=XKU*Two*(Occ(J)-Occ(K))*TwoMO(NAddr3(I,I,J,K))
     $+(DerFun(Occ(I),Occ(J),1,I,J)-DerFun(Occ(I),Occ(K),1,I,K))
     $ *TwoMO(NAddr3(I,J,K,I))
C
      BMATX(I)=BMATX(I)-Two*Fijk*UX(JK)
      BMATY(I)=BMATY(I)-Two*Fijk*UY(JK)
      BMATZ(I)=BMATZ(I)-Two*Fijk*UZ(JK)
C
      EndDo
      EndDo
      EndDo
C
      Return
      End

*Deck FunPair
      Real*8 Function FunPair(X,Y,ALPH,I,J)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Include 'commons.inc'
C
C     ** KU **
C 
      If(IFun.Eq.1) Then
C
      FacX=One
      If(X.Lt.Half) FacX=-One
      FacY=One
      If(Y.Lt.Half) FacY=-One
      FunPair=FacX*FacY*SQRT(X*Y)
      ALPH=Half
C
      Return
      EndIf
C
C     ** BB **
C
      If(IFun.Eq.2) Then
      FunPair=-SQRT(X*Y)
      ALPH=Half
      Return
      EndIf
C
C     ** GU **
C
      If(IFun.Eq.3) Then
C
      FunPair=-SQRT(X*Y)
      ALPH=Half
      If(I.Eq.J) Then
      FunPair=FunPair+X-X**2     
      ALPH=One
      EndIf 
C
      Return
      EndIf
C
C     ** BBC1 **
C
      If(IFun.Eq.4) Then
C
      FunPair=-SQRT(X*Y)
      ALPH=Half
      If(I.Ne.J.And.I.Gt.NELE.And.J.Gt.NELE) FunPair=SQRT(X*Y)
C
      EndIf
C
C     ** BBC2 OR BBC3 **
C
      If(IFun.Eq.5.Or.IFun.Eq.6) Then
C
      IBOND1=NELE
      IBOND2=IMax
C
      FunPair=-SQRT(X*Y)
      ALPH=Half
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND DIFFERENT
C
      If(I.Ne.J.And.I.Gt.IBOND1.And.J.Gt.IBOND1) FunPair=SQRT(X*Y)
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND DIFFERENT
C
      If(I.Lt.IBOND1.And.J.Lt.IBOND1.And.I.Ne.J) Then
      FunPair=FunPair+Sqrt(X*Y)-X*Y
      ALPH=One
      EndIf
C
C     ORBITAL I IS STRONGLY AND J IS BONDING
C
      If((I.Lt.IBOND1.And.J.Eq.IBOND1).Or.
     $ (I.Lt.IBOND1.And.J.Eq.IBOND2)) Then
      FunPair=FunPair+Sqrt(X*Y)-X*Y
      ALPH=One
      EndIf
C
C     ORBITAL J IS STRONGLY AND I IS BONDING
C
      If((J.Lt.IBOND1.And.I.Eq.IBOND1).Or.
     $ (J.Lt.IBOND1.And.I.Eq.IBOND2)) Then
      FunPair=FunPair+Sqrt(X*Y)-X*Y
      ALPH=One
      EndIf
C
      If(IFun.Eq.5) Return
C
C     BOTH ORBITALS ARE WEAKLY OCCUPIED AND THE SAME 
C
      If(I.Eq.J.And.I.Ne.IBOND1.And.I.Ne.IBOND2) Then
      FunPair=FunPair+X-X**2
      ALPH=One
      EndIf
C
C     ORBITALS I AND J ARE BOTH STRONLGY OCCUPIED AND THE SAME
C
      If(I.Eq.J.And.I.Lt.IBOND1.And.J.Lt.IBOND1) Then
      FunPair=FunPair+X-X**2
      ALPH=One
      EndIf
C
      EndIf
C
      Return
      End

*Deck ABMAT
      Subroutine ABMAT(IFlag,AMAT,BMATX,BMATY,BMATZ,
     $ URe,Occ,XKin,XNuc,TwoMO,
     $ CMAT,DipNOX,DipNOY,DipNOZ,Occ1X,Occ1Y,Occ1Z,
     $ NBasis,NDim,NInte1,NInte2) 
C
C     COMPUTE THE MATRICES B0 AND A NEEDED FOR CP-DMFT
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension AMAT(NDim,NDim),BMATX(NDim),BMATY(NDim),BMATZ(NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XKin(NInte1),XNuc(NInte1),TwoMO(NInte2),
     $ CMAT(NBasis,NBasis,NBasis),
     $ DipNOX(NInte1),DipNOY(NInte1),DipNOZ(NInte1),
     $ Occ1X(NBasis),Occ1Y(NBasis),Occ1Z(NBasis)
C
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
C     SET UP BMAT VECTORS
C
      IJ=0
      Do I=1,NBasis
      II=(I*(I+1))/2
C
      Do J=1,I-1
      JJ=(J*(J+1))/2
      IJ=IJ+1
      IJ1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
C
      BMATX(IJ)=-DipNOX(IJ1)*(Occ(I)-Occ(J))**2
      BMATY(IJ)=-DipNOY(IJ1)*(Occ(I)-Occ(J))**2
      BMATZ(IJ)=-DipNOZ(IJ1)*(Occ(I)-Occ(J))**2
C
      Do K=1,NBasis
C
      ExInt=TwoMO(NAddr3(I,K,K,J))
C
      Term1=XKU*Two*(Occ(I)-Occ(J))*TwoMO(NAddr3(I,J,K,K))
     $ +ExInt*(DerFun(Occ(K),Occ(I),1,K,I)-DerFun(Occ(K),Occ(J),1,K,J))
      Term1=Term1*(Occ(I)-Occ(J))
C
      DFik=DerFun(Occ(I),Occ(K),1,I,K)*(Occ(I)-Occ(J))
      DFjk=DerFun(Occ(J),Occ(K),1,J,K)*(Occ(I)-Occ(J))
C
      Term2X=-Occ1X(I)*DFik+Occ1X(J)*DFjk
      Term2Y=-Occ1Y(I)*DFik+Occ1Y(J)*DFjk
      Term2Z=-Occ1Z(I)*DFik+Occ1Z(J)*DFjk
C
      Fijk=FunPair(Occ(I),Occ(K),ALPH,I,K)
     $     -FunPair(Occ(J),Occ(K),ALPH,J,K)
C
      Term2X=Term2X+Fijk*(Occ1X(I)-Occ1X(J))
      Term2Y=Term2Y+Fijk*(Occ1Y(I)-Occ1Y(J))
      Term2Z=Term2Z+Fijk*(Occ1Z(I)-Occ1Z(J))
C
      BMATX(IJ)=BMATX(IJ)-Occ1X(K)*Term1+Term2X*ExInt
      BMATY(IJ)=BMATY(IJ)-Occ1Y(K)*Term1+Term2Y*ExInt
      BMATZ(IJ)=BMATZ(IJ)-Occ1Z(K)*Term1+Term2Z*ExInt
C
      EndDo
C
c herer!!!
c      BMATX(IJ)=BMATX(IJ)/(Occ(I)-Occ(J))
c      BMATY(IJ)=BMATY(IJ)/(Occ(I)-Occ(J))
c      BMATZ(IJ)=BMATZ(IJ)/(Occ(I)-Occ(J))
      EndDo
      EndDo
C
      If(IFlag.Eq.1) Return
C
C     SET UP AN AUXILIARY MATRIX CMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
C
      HIJ=Zero
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HIJ=HIJ+URe(I,IA)*URe(J,IB)*(XKin(IAB)+XNuc(IAB))
      EndDo
      EndDo
C
      Do K=1,NBasis
C
      CMAT(I,J,K)=HIJ*(Occ(I)-Occ(K))
C
      Do L=1,NBasis
C
      Fikl=(FunPair(Occ(I),Occ(L),ALPH,I,L)
     $ -FunPair(Occ(K),Occ(L),ALPH,K,L))
C 
      CMAT(I,J,K)=CMAT(I,J,K)
     $ +XKU*Two*(Occ(I)-Occ(K))*Occ(L)*TwoMO(NAddr3(I,J,L,L))
     $ +Fikl*TwoMO(NAddr3(I,L,L,J))
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     SET UP AN AMAT MATRIX 
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1      
      IJ=IJ+1
C
      KL=0
      Do K=1,NBasis
C
      Fijk=(FunPair(Occ(I),Occ(K),ALPH,I,K)
     $ -FunPair(Occ(J),Occ(K),ALPH,J,K))
C
      Do L=1,K-1
      KL=KL+1
C
      Fijl=(FunPair(Occ(I),Occ(L),ALPH,I,L)
     $ -FunPair(Occ(J),Occ(L),ALPH,J,L))
C
      AMAT(IJ,KL)=XKU*Two*(Occ(I)-Occ(J))*(Occ(K)-Occ(L))*
     $ (TwoMO(NAddr3(I,J,K,L))+TwoMO(NAddr3(I,J,L,K)))
     $ +(Fijk-Fijl)
     $ *(TwoMO(NAddr3(I,L,K,J))+TwoMO(NAddr3(I,K,L,J)))
C
      If(K.Eq.J) AMAT(IJ,KL)=AMAT(IJ,KL)+CMAT(I,L,K)
      If(L.Eq.J) AMAT(IJ,KL)=AMAT(IJ,KL)-CMAT(I,K,L)
      If(I.Eq.K) AMAT(IJ,KL)=AMAT(IJ,KL)-CMAT(J,L,K)
      If(I.Eq.L) AMAT(IJ,KL)=AMAT(IJ,KL)+CMAT(J,K,L)
C
c herer!!!
      AMAT(IJ,KL)=AMAT(IJ,KL)*(Occ(I)-Occ(J))
C
      EndDo 
      EndDo
C
      EndDo
      EndDo
C
      Return
      End
C
*Deck svdcmp  
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=5000)
CU    USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
C
      if(np.gt.nmax) stop'FATAL ERROR IN svdcmp!'
C   
      g=0.0D0
      scale=0.0D0
      anorm=0.0D0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0D0
        s=0.0D0
        scale=0.0D0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if(scale.ne.0.0D0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0D0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0D0
        s=0.0D0
        scale=0.0D0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if(scale.ne.0.0D0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0D0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0D0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0D0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0D0
            v(j,i)=0.0D0
31        continue
        endif
        v(i,i)=1.0D0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0D0
33      continue
        if(g.ne.0.0D0)then
          g=1.0D0/g
          do 36 j=l,n
            s=0.0D0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0D0
38        continue
        endif
        a(i,i)=a(i,i)+1.0D0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,50
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0D0
          s=1.0D0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0D0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0D0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.500) stop 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0D0*h*y)
          g=pythag(f,1.0D0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0D0
          s=1.0D0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0D0)then
              z=1.0D0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0D0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END

*Deck pythag
      Real*8 FUNCTION pythag(a,b)
      REAL*8 a,b
      REAL*8 absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.D0+(absb/absa)**2)
      else
        if(absb.eq.0.D0)then
          pythag=0.D0
        else
          pythag=absb*sqrt(1.D0+(absa/absb)**2)
        endif
      endif
      return
      END

