*Deck OPTNORB 
      Subroutine OPTNORB(ETot,URe,Occ,XKin,XNuc,
     $ DipX,DipY,DipZ,TwoEl,UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,
     $ NGem,NGOcc,IModG)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C
c herer!!!
      Parameter(GTol=1.D-3,ETol=1.D-4,MxIt=200)
c      Parameter(GTol=1.D-6,ETol=1.D-10,MxIt=100)
C  
      Dimension URe(NBasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(NInte1),DipX(NInte1),DipY(NInte1),DipZ(NInte1),
     $ TwoEl(NInte2),UMOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension F(NGem,NInte1),T(NBasis,NBasis),G(NBasis,NInte1),
     $ A(NBasis,NBasis),FunIJ(NBasis,NBasis),
     $ UReO(NBasis,NBasis),CoulNO(NBasis*(NBasis+1)/2),
     $ ExchNO(NBasis*(NBasis+1)/2),NSymMO(NBasis),NSymNO(NBasis)
C
      Real*8, Allocatable :: TNO(:)
C
      IAlloc=0
C
      MaxX=NBasis*(NBasis-1)/2
C
      If(IFun.Eq.13) Then
C
C     APSG OPTIMIZATION
C
      MaxXV=MaxX
      MxHVec=MaxXV+NBasis
      Call OptAPSGFun(ETot,URe,Occ,XKin,XNuc,
     $ TwoEl,UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,MaxXV,MxHVec,
     $  NGem,NGOcc,IModG)
c herer!!!
c this is needed to generate a guess for a delocalized solution for H4_square !!!
C do not remove
c      IFreeze=0
c      NGOcc=0
c      IModG=1
c      Call OptAPSGFun(ETot,URe,Occ,XKin,XNuc,
c     $ TwoEl,UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,MaxXV,MxHVec,
c     $  NGem,NGOcc,IModG)
c
c      IGVB=1
c      ETot=0.0D0
c      Call OptAPSGFun(ETot,URe,Occ,XKin,XNuc,
c     $ TwoEl,UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,MaxXV,MxHVec,
c     $  NGem,NGOcc,IModG)

      Return
C     
      EndIf
C
      If(IFun.Eq.10) Stop 'Fatal Error: OPTNORB not defined for IFun=10'
C
C     FOR A GIVEN TRANSFORMATION MATRIX URe CONSTRUCT THE MATRICES F,A, AND G, 
C     AND COMPUTE THE TOTAL ENERGY
C
      If(IPrint.Ge.1)
     $ Write(6,'(/,X,''****************************************'',
     $     ''***************************************'')')
C
      Do I=1,NBasis
      NSymNO(I)=Zero
      EndDo
C
c      IDFP=0
      IDFP=1
C
      If(IFun.Eq.13) Then
      Call TrTwoEl(CoulNO,ExchNO,URe,TwoEl,NBasis,NInte2)
      Call OptAPSG(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem)
      EndIf
C
      Do 100 It=1,MxIt
C
      ETotO=ETot
C
      If(It.Gt.1) Call ReWr(1,Occ,URe,Title,NBasis)
C
C     CALL EigenNO ONLY IF THE SYMMETRY IS NOT IMPOSED
C
      If(NoSym.Eq.1.And.IFun.Ne.13) 
     $ Call EigenNO(ETot,URe,Occ,XKin,XNuc,TwoEl,NBasis,NInte1,NInte2)
C
      Call GetFij(FunIJ,Occ,NBasis)
C
      If(IDFP.Eq.0) Then 
C
      If(IAlloc.Eq.0) Then
      Allocate (TNO(NInte2))
      IAlloc=1
      EndIf 
      Call TwoNO(TNO,URe,TwoEl,NBasis,NInte2)
C
      Call FAG(F,A,G,ETot,URe,Occ,FunIJ,XKin,XNuc,TNO,
     $ NBasis,NInte1,NInte2,NGem)
C
      Do 15 I=1,NBasis
      Do 15 J=1,NBasis
   15 T(I,J)=Zero
C
      Call DFPMIN(IDFP,F,A,G,FunIJ,Occ,T,XKin,XNuc,TNO,
     $ MaxX,NBasis,NInte2,NGem,ItDFPMIN)
C
C     NEW URe
C
      Do 20 I=1,NBasis
   20 T(I,I)=T(I,I)+One
      Call CpyM(UReO,URe,NBasis) 
      Call MultpM(URe,T,UReO,NBasis) 
C
C     Else to If(IDFP.Eq.0)
      Else
C
      Call CpyM(UReO,URe,NBasis)
      Call DFPMIN(IDFP,F,A,G,FunIJ,Occ,URe,XKin,XNuc,TwoEl,
     $ MaxX,NBasis,NInte2,NGem,ItDFPMIN)
C      
      EndIf   
C
C     CHECK THE SYMMETRY
C
      If(NoSym.Eq.0) Then
C
      Do I=1,NBasis

      NSymNO(I)=0

      Do J=1,NBasis
C
      If(Abs(URe(I,J)).Ne.Zero) Then
C
      ISym=NSymMO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym)
     $ Write(*,*)'Symm of NO',I,' cannot be established',j
      EndIf
C
      EndIf
      EndDo
      EndDo
C
      EndIf
C
      Call TrTwoEl(CoulNO,ExchNO,URe,TwoEl,NBasis,NInte2)
C
      If(IFun.Ne.13) Then
C
      Call OPTN(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,NBasis)
C
      Else
C
      If(IModG.Eq.0) Then
      Call OptAPSG(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem)
      Else
      Call OptArai(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem,IModG,NGOcc)
      EndIf  
C
      EndIf
      Call SortOcc(Occ,URe,NBasis)
C
      If(ETotO.Lt.Etot.And.It.Gt.1) Then
      If(IPrint.Ge.1)
     $ Write(6,'(''Energy has risen in DFPMIN! IDFP SET TO 1'')')
c      Call CpyM(URe,UReO,NBasis)
      IDFP=1
      If(IAlloc.Eq.1) Deallocate (TNO) 
      IAlloc=0
      EndIf
C
      Call GetFij(FunIJ,Occ,NBasis)
      If(IFun.Ne.13) Then
      Call ETA(ETot,A,URe,Occ,FunIJ,XKin,XNuc,TwoEl,NBasis,NInte2)
      Else
      Call ETAAPSG(ETot,A,URe,Occ,FunIJ,XKin,XNuc,TwoEl,NBasis,NInte2,
     $ NGem)
      EndIf 
C
C     CALCULATE THE NORM OF THE GRADIENT
C
      GNorm=Zero
      IJ=0
      Do 5 I=1,NBasis
      Do 5 J=1,I-1
    5 GNorm=GNorm+(A(I,J)-A(J,I))**2
      GNorm=Sqrt(GNorm)
C
      If(IPrint.Ge.0)
     $ Write(6,'(/,X,''MACRO ITER'',I3,2X,''ENERGY'',F16.8,2X,
     $ ''ENE DIFF '',E10.3,2X,''GRAD NORM '',E9.3)')
     $ It,ETot,ETot-ETotO,GNorm
C
      If(IPrint.Ge.1) 
     $ Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      If(GNorm.Le.GTol.And.ETotO-ETot.Le.ETol) Then
      If(IPrint.Ge.1) Write(6,'(/,10X,''TOTAL CONVERGENCE ATTAINED '')')
      Call ReWr(1,Occ,URe,Title,NBasis)
      Return
c herer!!!
      ElseIf(ETotO-ETot.Le.ETol.And.ItDFPMIN.Eq.1) Then 
c      Write(6,'(/,1X,''CONV CRITERIA NOT MET BUT ItDFPMIN.Eq.1'')')
c      Call ReWr(1,Occ,URe,Title,NBasis)   
c      Return
      EndIf
C
  100 Continue           
C
      Return
      End

*Deck E2
      Subroutine E2(E,F,G,T,FunIJ,B,Occ,TNO,NBasis,NInte2,IFlag,NGem)
C
C     FOR GIVEN MATRICES F,G, AND T - CALCULATE E2 
C
C     IFlag = 1 - CALCULATE A B MATRIX (MAY BE LATER USED IN GRX)
C           = 0 - DO NOT CALCULATE A B MATRIX (TWO TIMES FASTER) 
C      
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Tol2=1.D-12) 
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,
     $ F4=4.D0,F16=16.D0)
C
      Dimension F(NGem,NBasis*(NBasis+1)/2),
     $ G(NBasis,NBasis*(NBasis+1)/2),
     $ T(NBasis,NBasis),FunIJ(NBasis,NBasis),TNO(NInte2),Occ(NBasis),
     $ B(NBasis,NBasis)
C
      Character*60 FMultTab
      Include 'commons.inc'
C 
      If(IFlag.Eq.1) Then
C
      E=Zero
C
      IJ=0
      Do 5 I=1,NBasis
      Do 5 J=1,I
      IJ=IJ+1
      FIIJ=F(IGem(I),IJ)
      FJIJ=F(IGem(J),IJ)
C
      Fac=One
      If(I.Eq.J) Fac=Half
    5 E=E+Fac*( (Occ(I)*FIIJ+G(I,IJ))*T(I,J)
     $         +(Occ(J)*FJIJ+G(J,IJ))*T(J,I) )
C
      E=F4*E
C
      Do 7 I=1,NBasis
      Do 7 J=1,NBasis
      B(J,I)=Zero
    7 Continue
C
      KJ=0
      Do 15 K=1,NBasis
      Do 15 J=1,K
      KJ=KJ+1
C
      Do 15 I=1,NBasis
      FIKJ=F(IGem(I),KJ)
      B(J,I)=B(J,I)+F4*T(I,K)*(Occ(I)*FIKJ+G(I,KJ))
      If(K.Ne.J) B(K,I)=B(K,I)+F4*T(I,J)*(Occ(I)*FIKJ+G(I,KJ))
   15 Continue
C
C     A SECOND PART OF B
C
      Fac16=F16
      If(IFun.Eq.1) Fac16=Zero
      NAddr=0
C
      IJ=0
      Do 35 I=1,NBasis
      Do 35 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 45 K=1,NBasis
C
      GemIK=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(K)) GemIK=Zero
      GemJK=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(K)) GemJK=Zero
C
      Do 45 L=1,K
      KL=KL+1
C
      GemIL=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(L)) GemIL=Zero
      GemJL=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(L)) GemJL=Zero
C
      If(IJ.Lt.KL) GoTo 45
      NAddr=NAddr+1
      TwoZet=TNO(NAddr)
      F16Z=Fac16*TwoZet
      F4Z=F4*TwoZet
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      B(J,I)=B(J,I)+F16Z*GemIK*Occ(I)*Occ(K)*T(K,L)
      B(L,I)=B(L,I)+F4Z*(FunIJ(I,K)*T(K,J)+FunIJ(I,J)*T(J,K))
C
      If(I.Ne.J) Then

      B(I,J)=B(I,J)+F16Z*GemJK*Occ(J)*Occ(K)*T(K,L)
      B(L,J)=B(L,J)+F4Z*(FunIJ(J,K)*T(K,I)+FunIJ(J,I)*T(I,K))
C
      If(K.Ne.L) Then
      B(I,J)=B(I,J)+F16Z*GemJL*Occ(J)*Occ(L)*T(L,K)
      B(K,J)=B(K,J)+F4Z*(FunIJ(J,L)*T(L,I)+FunIJ(J,I)*T(I,L))
      EndIf

      EndIf
C
      If(K.Ne.L) Then
      B(J,I)=B(J,I)+F16Z*GemIL*Occ(I)*Occ(L)*T(L,K)
      B(K,I)=B(K,I)+F4Z*(FunIJ(I,L)*T(L,J)+FunIJ(I,J)*T(J,L))
      EndIf
C
      If(IJ.Ne.KL) Then
C
      B(L,K)=B(L,K)+F16Z*GemIK*Occ(K)*Occ(I)*T(I,J)
      B(J,K)=B(J,K)+F4Z*(FunIJ(K,I)*T(I,L)+FunIJ(K,L)*T(L,I))
C
      If(I.Ne.J) Then

      B(L,K)=B(L,K)+F16Z*GemJK*Occ(K)*Occ(J)*T(J,I)
      B(I,K)=B(I,K)+F4Z*(FunIJ(K,J)*T(J,L)+FunIJ(K,L)*T(L,J))
C
      If(K.Ne.L) Then
      B(K,L)=B(K,L)+F16Z*GemJL*Occ(L)*Occ(J)*T(J,I)
      B(I,L)=B(I,L)+F4Z*(FunIJ(L,J)*T(J,K)+FunIJ(L,K)*T(K,J))
      EndIf

      EndIf
C
      If(K.Ne.L) Then
      B(K,L)=B(K,L)+F16Z*GemIL*Occ(L)*Occ(I)*T(I,J)
      B(J,L)=B(J,L)+F4Z*(FunIJ(L,I)*T(I,K)+FunIJ(L,K)*T(K,I))
      EndIf
C
      EndIf
C
      EndIf
C
   45 Continue
   35 Continue
C
      Do 55 I=1,NBasis
      Do 55 J=1,NBasis
   55 E=E+Half*B(J,I)*T(I,J)
C
      Return
C
      Else
C
      E=Zero
C
      IJ=0
      Do 10 I=1,NBasis
      Do 10 J=1,I
      IJ=IJ+1
C
      FIIJ=F(IGem(I),IJ)
      FJIJ=F(IGem(J),IJ)
C
      Fac=One
      If(I.Eq.J) Fac=Half     
      E=E+Fac*( (Occ(I)*FIIJ+G(I,IJ))*T(I,J)
     $         +(Occ(J)*FJIJ+G(J,IJ))*T(J,I) )
C
      Do 10 K=1,NBasis
      FKIJ=F(IGem(K),IJ)
   10 E=E+T(K,J)*T(K,I)*Fac*(Occ(K)*FKIJ+G(K,IJ))  
C
      E=F4*E
C
C     SECOND PART (QUADRATIC IN T) 
C
      EE=Zero

      Fac4=F4
      If(IFun.Eq.1) Fac4=Zero
      NAddr=0
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 40 K=1,NBasis
C
      GemIK=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(K)) GemIK=Zero
      GemJK=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(K)) GemJK=Zero
C
      Do 40 L=1,K
      KL=KL+1
C
      GemIL=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(L)) GemIL=Zero
      GemJL=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(L)) GemJL=Zero
C
      If(IJ.Lt.KL) GoTo 40
      NAddr=NAddr+1
      TwoZet=TNO(NAddr)
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      Fac=Two
      If(IJ.Eq.KL)Fac=One
C
      T1=FunIJ(I,J)*T(I,L)*T(J,K)
      T2=FunIJ(J,I)*T(J,L)*T(I,K)
C
      EE=EE+(Fac*(Fac4*GemIK*Occ(I)*Occ(K)*T(I,J)*T(K,L)+
     $ FunIJ(I,K)*T(I,L)*T(K,J))+T1)*TwoZet
C
      If(I.Ne.J) Then

      EE=EE+(Fac*(Fac4*GemJK*Occ(J)*Occ(K)*T(J,I)*T(K,L)+
     $ FunIJ(J,K)*T(J,L)*T(K,I))+T2)*TwoZet
C
      If(K.Ne.L) EE=EE+(Fac*(Fac4*GemJL*Occ(J)*Occ(L)*T(J,I)*T(L,K)
     $ +FunIJ(J,L)*T(J,K)*T(L,I))+T1)*TwoZet

      EndIf
C
      If(K.Ne.L) EE=EE+(Fac*(Fac4*GemIL*Occ(I)*Occ(L)*T(I,J)*T(L,K)+
     $ FunIJ(I,L)*T(I,K)*T(L,J))+T2)*TwoZet
C
      If(IJ.Ne.KL) Then
C
      T1=FunIJ(K,L)*T(K,J)*T(L,I)*TwoZet
      T2=FunIJ(K,L)*T(K,I)*T(L,J)*TwoZet
C
      EE=EE+T1
C
      If(I.Ne.J) Then
      EE=EE+T2
      If(K.Ne.L) EE=EE+T1
      EndIf
C
      If(K.Ne.L) EE=EE+T2
C
      EndIf
C
      EndIf
C
   40 Continue
   30 Continue
C
      E=E+Two*EE
C
      EndIf
C
      Return
      End  

*Deck GRX
      Subroutine GRX(GX,F,A,G,T,FunIJ,B,Occ,TNO,NBasis,NInte2,IFlag,
     $ NGem)
C
C     CONSTRUCT A GRADIENT FOR A GIVEN T
C
C     IFlag = 1 - A PART OF THE MATRIX B IS KNOWN (FROM THE E2 CALCULATION)
C           = 0 - CALCULATE A TOTAL B MATRIX 
C           = 2 - GRAD(IJ) = A(IJ) - A(JI)           
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Tol2=1.D-12)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,F4=4.D0,
     $ F16=16.D0)
C
      Dimension GX(NBasis*(NBasis-1)/2),F(NGem,NBasis*(NBasis+1)/2),
     $ A(NBasis,NBasis),G(NBasis,NBasis*(NBasis+1)/2),
     $ Occ(NBasis),B(NBasis,NBasis),T(NBasis,NBasis),
     $ TNO(NInte2),FunIJ(NBasis,NBasis)
C
C     LOCAL ARRAY
C
      Dimension SB(NBasis,NBasis)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     !!! B IS TRANSPOSED!!!
C
      If(IFlag.Eq.1) Then
C
      Do 5 I=1,NBasis
      Do 5 J=1,NBasis
    5 B(J,I)=B(J,I)+A(J,I)
C
      ElseIf(IFlag.Eq.2) Then
C
      Call CpyM(B,A,NBasis)
C
      Do 2 I=1,NBasis
      Do 2 J=1,NBasis
    2 SB(J,I)=Zero     
      GoTo 999 
C
      Else
C       
      Call CpyM(B,A,NBasis)
C
      KJ=0
      Do 20 K=1,NBasis
      Do 20 J=1,K
      KJ=KJ+1
C
      Do 20 I=1,NBasis
      FIKJ=F(IGem(I),KJ)
      B(J,I)=B(J,I)+F4*T(I,K)*(Occ(I)*FIKJ+G(I,KJ))
      If(K.Ne.J) B(K,I)=B(K,I)+F4*T(I,J)*(Occ(I)*FIKJ+G(I,KJ))
   20 Continue
C
C     A SECOND PART OF B 
C
      Fac16=F16
      If(IFun.Eq.1) Fac16=Zero
      NAddr=0
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 40 K=1,NBasis
C
      GemIK=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(K)) GemIK=Zero
      GemJK=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(K)) GemJK=Zero
C
      Do 40 L=1,K
      KL=KL+1
C
      GemIL=One
      If(IFun.Eq.13.And.IGem(I).Eq.IGem(L)) GemIL=Zero
      GemJL=One
      If(IFun.Eq.13.And.IGem(J).Eq.IGem(L)) GemJL=Zero
C
      If(IJ.Lt.KL) GoTo 40
      NAddr=NAddr+1
      TwoZet=TNO(NAddr)
      F16Z=Fac16*TwoZet
      F4Z=F4*TwoZet
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      B(J,I)=B(J,I)+F16Z*GemIK*Occ(I)*Occ(K)*T(K,L)
      B(L,I)=B(L,I)+F4Z*(FunIJ(I,K)*T(K,J)+FunIJ(I,J)*T(J,K))
C
      If(I.Ne.J) Then

      B(I,J)=B(I,J)+F16Z*GemJK*Occ(J)*Occ(K)*T(K,L)
      B(L,J)=B(L,J)+F4Z*(FunIJ(J,K)*T(K,I)+FunIJ(J,I)*T(I,K))
C
      If(K.Ne.L) Then
      B(I,J)=B(I,J)+F16Z*GemJL*Occ(J)*Occ(L)*T(L,K)
      B(K,J)=B(K,J)+F4Z*(FunIJ(J,L)*T(L,I)+FunIJ(J,I)*T(I,L))
      EndIf

      EndIf
C
      If(K.Ne.L) Then
      B(J,I)=B(J,I)+F16Z*GemIL*Occ(I)*Occ(L)*T(L,K)
      B(K,I)=B(K,I)+F4Z*(FunIJ(I,L)*T(L,J)+FunIJ(I,J)*T(J,L))
      EndIf
C
      If(IJ.Ne.KL) Then
C
      B(L,K)=B(L,K)+F16Z*GemIK*Occ(K)*Occ(I)*T(I,J)
      B(J,K)=B(J,K)+F4Z*(FunIJ(K,I)*T(I,L)+FunIJ(K,L)*T(L,I))
C
      If(I.Ne.J) Then

      B(L,K)=B(L,K)+F16Z*GemJK*Occ(K)*Occ(J)*T(J,I)
      B(I,K)=B(I,K)+F4Z*(FunIJ(K,J)*T(J,L)+FunIJ(K,L)*T(L,J))
C
      If(K.Ne.L) Then
      B(K,L)=B(K,L)+F16Z*GemJL*Occ(L)*Occ(J)*T(J,I)
      B(I,L)=B(I,L)+F4Z*(FunIJ(L,J)*T(J,K)+FunIJ(L,K)*T(K,J))
      EndIf

      EndIf
C
      If(K.Ne.L) Then
      B(K,L)=B(K,L)+F16Z*GemIL*Occ(L)*Occ(I)*T(I,J)
      B(J,L)=B(J,L)+F4Z*(FunIJ(L,I)*T(I,K)+FunIJ(L,K)*T(K,I))
      EndIf
C
      EndIf
C 
      EndIf 
C
   40 Continue
   30 Continue
C
      EndIf
C
      Call MultpM(SB,B,T,NBasis)
C
  999 Continue 
C
      IJ=0
      Do 50 I=1,NBasis
      Do 50 J=1,I-1
      IJ=IJ+1
   50 GX(IJ)=-SB(I,J)+SB(J,I)+B(J,I)-B(I,J)
C
      Return
      End

*Deck EXPM
      Subroutine EXPM(T,X,N,IFl)
C
C     CALCULATE Exp(X)-1
C     T = X + (1/2) X**2 + (1/3!) X**3 + ...
C
C     IFl = 0  - GET THE APPROXIMATION TO T UP TO MxPow - ORDER
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.D0,Tol=1.D-12, MxIt=50, MxPow=3)
C
      Dimension T(N*N),X(N*(N-1)/2)
C
C     LOCAL ARRAYS
C
      Dimension XN(N*N),XNN(N*N),X1(N*N)
C
      Do 5 I=1,N*N
      T(I)=Zero
    5 XN(I)=Zero
C
      IJ=0
      Do 10 I=1,N
      II=(I-1)*N+I 
      X1(II)=Zero
      XN(II)=One
C
      Do 10 J=1,I-1
      IJ=IJ+1
      IJ1=(J-1)*N+I
      JI1=(I-1)*N+J
      X1(IJ1)=X(IJ)
   10 X1(JI1)=-X(IJ)
C
      Pow=One
C
      Do 1 It=1,MxIt
C  
      Pow=Pow*Float(It)
      Call MultpM(XNN,XN,X1,N)
C
      Err=Zero
C
      Do 20 I=1,N*N
      T(I)=T(I)+XNN(I)/Pow
      XN(I)=XNN(I) 
   20 Err=Err+XNN(I)**2
      Err=Sqrt(Err) 
C
      If((IFl.Eq.0.And.It.Eq.MxPow).Or.(Err.Lt.Tol)) Return
    1 Continue      
C
      Write(6,'("ERROR IN EXPM TOO HIGH! ",E12.2)')Err
      Stop 'No of iterations exceeded in EXPM!'
C
      Return
      End

*Deck NewU
      Subroutine NewU(U,UOld,Step,N,IFl)
C
C     CALCULATE A NEW U MATRIX
C     U = Exp(Step) * UOld 
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension U(N*N),UOld(N*N),Step(N*(N-1)/2)
C
C     LOCAL ARRAY
C  
      Dimension UHlp(N*N)
C
      Call EXPM(UHlp,Step,N,IFl)
C
      Do 5 I=1,N
      II=N*(I-1)+I
    5 UHlp(II)=UHlp(II)+1.D0 
C
      Call MultpM(U,UHlp,UOld,N)
C
      Return
      End

*Deck NewT
      Subroutine NewT(T,TOld,Step,N,IFl)
C
C     CALCULATE A NEW T MATRIX 
C     T = (TOld+1)*Exp(Step) - 1
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension T(N*N),TOld(N*N),Step(N*(N-1)/2)
C
C     LOCAL ARRAY
C
      Dimension THlp(N*N)
C
      Call EXPM(THlp,Step,N,IFl)
C   
      Call MultpM(T,TOld,THlp,N) 
C
      Do 10 I=1,N*N
   10 T(I)=T(I)+TOld(I)+THlp(I)   
C
      Return
      End

*Deck ETA
      Subroutine ETA(ETot,A,URe,Occ,FunIJ,XKin,XNuc,TwoEl,NBasis,NInte2)
C
C     CALCULATE THE TOTAL ENERGY AND AN A MATRIX FOR ALL FUNCTIONALS EXCEPT FOR APSG
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
      Dimension A(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),
     $ FunIJ(NBasis,NBasis),XKin(NBasis*(NBasis+1)/2),
     $ XNuc(NBasis*(NBasis+1)/2),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension
     $ F(NBasis*(NBasis+1)/2),Gamma(NBasis*(NBasis+1)/2),
     $ Xikkj(NBasis,NBasis*(NBasis+1)/2),
     $ G(NBasis,NBasis*(NBasis+1)/2)
C
      NInte1=NBasis*(NBasis+1)/2
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
      Do 37 I=1,NInte1
      Do 37 J=1,NBasis
   37 G(J,I)=Zero
C
C     THE EXCHANGE-CORRELATION PART
C
      IJ=0
      Do 40 I=1,NBasis
      Do 40 J=1,I
      IJ=IJ+1
C
      Do 40 K=1,NBasis
C
      XKIJ=Xikkj(K,IJ)
C
      Do 50 L=1,NBasis
   50 G(L,IJ)=G(L,IJ)+FunIJ(K,L)*XKIJ
C
      If(I.Eq.J) ETot=ETot+FunIJ(I,K)*XKIJ
C
   40 Continue
C
C     CONSTRUCT A MATRIX A
C
      IJ=0
      Do 60 I=1,NBasis
      Do 60 J=1,I
      IJ=IJ+1
C
      Fij=F(IJ)
      A(J,I)=Four*(Occ(I)*Fij+G(I,IJ))
      If(I.Ne.J) A(I,J)=Four*(Occ(J)*Fij+G(J,IJ))
   60 Continue
C
       Return
       End

*Deck ETAAPSG
      Subroutine ETAAPSG(ETot,A,URe,Occ,FunIJ,XKin,XNuc,TwoEl,NBasis,
     $NInte2,NGem)
C
C     CALCULATE THE TOTAL ENERGY AND AN A MATRIX FOR APSG
C
C     Tol2 - Threshold for two-el integrals
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0,
     $ Tol2=1.D-15)
C
      Dimension A(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),
     $ FunIJ(NBasis,NBasis),XKin(NBasis*(NBasis+1)/2),
     $ XNuc(NBasis*(NBasis+1)/2),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension 
     $ Gamma(NBasis*(NBasis+1)/2),
     $ Gam(NGem,NBasis*(NBasis+1)/2),
     $ GamC(NGem,NBasis*(NBasis+1)/2),
     $ F(NGem*NBasis*(NBasis+1)/2),
     $ GN(NGem*NBasis*(NBasis+1)/2),
     $ GC(NGem*NBasis*(NBasis+1)/2)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      NInte1=NBasis*(NBasis+1)/2 
C
C     CONSTRUCT NEW DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
C     AUXILIARY DENSITY-LIKE MATRICES
C
      Do IG=1,NGem
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA
      IAB=IAB+1
C
      Gam(IG,IAB)=Zero
      GamC(IG,IAB)=Zero
C
      Do K=1,NBasis
      If(IG.Ne.IGem(K)) Gam(IG,IAB)=Gam(IG,IAB)+Occ(K)*
     $ URe(K,IA)*URe(K,IB)
      If(IG.Eq.IGem(K)) GamC(IG,IAB)=GamC(IG,IAB)+CICoef(K)*
     $ URe(K,IA)*URe(K,IB)
      EndDo
C
      EndDo
      EndDo
C
      EndDo
C
      Do IG=1,NGem
      Do I=1,NInte1
      F((IG-1)*NInte1+I)=XKin(I)+XNuc(I)
      GN((IG-1)*NInte1+I)=Zero
      GC((IG-1)*NInte1+I)=Zero
      EndDo
      EndDo
C
C     CONTRACT Gam AND GamC WITH THE TWO-EL INTEGRALS
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
      If(Abs(TwoZet).Lt.Tol2) GoTo 40
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
      Do IG=1,NGem
C
      F((IG-1)*NInte1+IJ)=F((IG-1)*NInte1+IJ)+FacIJ*Gam(IG,KL)*TwoZet
      If(IJ.Ne.KL) F((IG-1)*NInte1+KL)=F((IG-1)*NInte1+KL)
     $ +FacKL*Gam(IG,IJ)*TwoZet
C
      GN((IG-1)*NInte1+IL)=GN((IG-1)*NInte1+IL)
     $ -Half*FacIL*Gam(IG,JK)*TwoZet
      If(IL.Ne.JK) GN((IG-1)*NInte1+JK)=GN((IG-1)*NInte1+JK)
     $ -Half*FacJK*Gam(IG,IL)*TwoZet
      GC((IG-1)*NInte1+IL)=GC((IG-1)*NInte1+IL)
     $ +Half*FacIL*GamC(IG,JK)*TwoZet
      If(IL.Ne.JK) GC((IG-1)*NInte1+JK)=GC((IG-1)*NInte1+JK)
     $ +Half*FacJK*GamC(IG,IL)*TwoZet
C
      EndDo
C
      If(K.Eq.L.Or.I.Eq.J) GoTo 40
C
      Do IG=1,NGem
C 
      F((IG-1)*NInte1+IJ)=F((IG-1)*NInte1+IJ)+Gam(IG,KL)*TwoZet
      If(IJ.Ne.KL) F((IG-1)*NInte1+KL)=F((IG-1)*NInte1+KL)
     $ +Gam(IG,IJ)*TwoZet
C
      GN((IG-1)*NInte1+IK)=GN((IG-1)*NInte1+IK)
     $ -Half*FacIK*Gam(IG,JL)*TwoZet
      If(IK.Ne.JL) GN((IG-1)*NInte1+JL)=GN((IG-1)*NInte1+JL)
     $ -Half*FacJL*Gam(IG,IK)*TwoZet
      GC((IG-1)*NInte1+IK)=GC((IG-1)*NInte1+IK)
     $ +Half*FacIK*GamC(IG,JL)*TwoZet
      If(IK.Ne.JL) GC((IG-1)*NInte1+JL)=GC((IG-1)*NInte1+JL)
     $ +Half*FacJL*GamC(IG,IK)*TwoZet
C
      EndDo
C                                                  
   40 Continue
   30 Continue
C
C     TRANSFORM F, GN, AND GC TO NO
C
      Do IG=1,NGem
C
      Call MatTr(F((IG-1)*NInte1+1),URe,NBasis)
      Call MatTr(GN((IG-1)*NInte1+1),URe,NBasis)
      Call MatTr(GC((IG-1)*NInte1+1),URe,NBasis)
C
      EndDo
C
C     COMPUTE THE ONE-ELECTRON AND COULOMB COMPONENTS OF ETOT
C
      ETot=Zero
C
      IJ=0
      Do I=1,NBasis
      IG=IGem(I)
      ETot=ETot+Occ(I)*F((IG-1)*NInte1+I*(I+1)/2)
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Gamma(IJ)*(XKin(IJ)+XNuc(IJ))
      EndDo
C
      EndDo
C
C     THE EXCHANGE-CORRELATION PART
C
      Do I=1,NBasis
      NAdr=(IGem(I)-1)*NInte1+I*(I+1)/2
      ETot=ETot+Occ(I)*GN(NAdr)+CICoef(I)*GC(NAdr)
      EndDo 
C
C     CONSTRUCT A MATRIX A
C
      IJ=0
      Do 60 I=1,NBasis
      Do 60 J=1,I
      IJ=IJ+1
C
      NAdr=(IGem(I)-1)*NInte1+IJ
      A(J,I)=Four*(Occ(I)*(F(NAdr)+GN(NAdr))+CICoef(I)*GC(NAdr))
      If(I.Ne.J) Then
      NAdr=(IGem(J)-1)*NInte1+IJ
      A(I,J)=Four*(Occ(J)*(F(NAdr)+GN(NAdr))+CICoef(J)*GC(NAdr))
      EndIf
C
   60 Continue
C
      Return
      End

*Deck FAG 
      Subroutine FAG(F,A,G,ETot,URe,Occ,FunIJ,XKin,XNuc,TNO,
     $ NBasis,NInte1,NInte2,NGem)
C
C     CONSTRUCT A FOCK MATRIX IN NO, MATRICES A (TRANSPOSED!) AND G, 
C     AND COMPUTE THE TOTAL ENERGY
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
      Dimension F(NGem,NInte1),URe(NBasis,NBasis),Occ(NBasis),
     $ XKin(NInte1),XNuc(NInte1),TNO(NInte2),
     $ A(NBasis,NBasis),G(NBasis,NInte1),
     $ FunIJ(NBasis,NBasis) 
C
C     LOCAL ARRAY
C
      Dimension Gamma(NInte1),FF(NInte1)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     CONSTRUCT NEW DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
      Do 25 I=1,NInte1
   25 FF(I)=XKin(I)+XNuc(I) 
C
      Call MatTr(FF,URe,NBasis)
C
      Do IG=1,NGem
C
      IJ=0
      Do 35 I=1,NBasis
      Do 35 J=1,I
      IJ=IJ+1
      F(IG,IJ)=FF(IJ)
C
      If(IFun.Ne.1) Then
C
      Do 34 K=1,NBasis
C
      If(IFun.Eq.13) Then
      If(IGem(K).Ne.IG)F(IG,IJ)=F(IG,IJ)+Two*Occ(K)*TNO(NAddr3(I,J,K,K))
      Else
      F(IG,IJ)=F(IG,IJ)+Two*Occ(K)*TNO(NAddr3(I,J,K,K))
      EndIf
C
   34 Continue
C
      EndIf
C
   35 Continue
C
      EndDo
C
C     COMPUTE THE ONE-ELECTRON AND THE COULOMB COMPONENTS OF ETOT
C
      ETot=Zero
C
      IJ=0
      Do 45 I=1,NBasis
      ETot=ETot+Occ(I)*F(IGem(I),I*(I+1)/2)
      Do 45 J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
   45 ETot=ETot+FacIJ*Gamma(IJ)*(XKin(IJ)+XNuc(IJ))
C
      Do 37 I=1,NInte1
      Do 37 J=1,NBasis
   37 G(J,I)=Zero
C
C     THE EXCHANGE-CORRELATION PART 
C
      IJ=0
      Do 40 I=1,NBasis
      Do 40 J=1,I
      IJ=IJ+1
C
      Do 40 K=1,NBasis
C
      XKIJ=TNO(NAddr3(I,K,J,K))
C
      Do 50 L=1,NBasis
   50 G(L,IJ)=G(L,IJ)+FunIJ(K,L)*XKIJ
C
      If(I.Eq.J) ETot=ETot+FunIJ(I,K)*XKIJ
C
   40 Continue
C
C     CONSTRUCT A MATRIX A
C
      IJ=0
      Do 60 I=1,NBasis 
      Do 60 J=1,I
      IJ=IJ+1
C
      FIIJ=F(IGem(I),IJ)
      FJIJ=F(IGem(J),IJ)
      A(J,I)=Four*(Occ(I)*FIIJ+G(I,IJ))
      If(I.Ne.J) A(I,J)=Four*(Occ(J)*FJIJ+G(J,IJ))
   60 Continue
C
      Return
      End     

*Deck GetFij
      Subroutine GetFij(FunIJ,Occ,NBasis)
C
C     CALCULATE F(ni,nj)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension FunIJ(NBasis,NBasis),Occ(NBasis)
C
      Do 10 I=1,NBasis
      Do 10 J=1,I
      FunIJ(J,I)=GOCC(Occ(I),Occ(J),0,I,J)
   10 FunIJ(I,J)=FunIJ(J,I)
C
      Return
      End

*Deck OptArai
      Subroutine OptArai(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem,IModG,NGOcc)
C
C     OPTIMIZATION OF THE ARAI SUBSPACES FOR APSG
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C     
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(NInte1),CoulNO(NInte1),ExchNO(NInte1),
     $ NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Dimension OccSave(NBasis),CISave(NBasis)
C     
C     OPTIMIZE OCCUPATIONS WITHOUT CHANGING ARAI SUBSPACES
C     GEMINALS FROM 1 TO NGOcc CONSIST OF ONE ORBITAL AND THEIR
C     OCCUPATION IS FIXED TO 1
C
      IPrintS=IPrint
      IPrint=-1
C
      Write(6,'(2/,X,'' ENERGY BEFORE MODIFYING ARAI SPACES'',
     $ F16.8)')ETot
      Call OptAPSG(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem) 
C
      Write(6,'(X,'' INITIAL ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem",3X,"Sym")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I),NSymNO(I)
      EndDo
C
      EIni=ETot
C
C     MOVE ORBITALS TO DIFFERENT GEMINALS AND CHECK THE OPTIMAL ENERGY
C
      Do IO=1,NBasis
C
      If(Occ(IO).Lt.Half) Then
C
      Do IG=NGOcc+1,NGem
C
      If(IG.Ne.IGem(IO)) Then
C
C     MOVE THE ORIBTAL TO THE GEMINAL IG
C
      IGSave=IGem(IO)
      Do J=1,NBasis
      OccSave(J)=Occ(J)
      CISave(J)=CICoef(J)
      EndDo
C
      IGem(IO)=IG
      Call OptAPSG(ETotN,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem)
C
      If(ETotN.Lt.ETot) Then
      ETot=ETotN
      Write(6,'('' ORBITAL NO '',I3,'' MOVED FROM GEMINAL'',I3,
     $ '' TO'',I3)')IO,IGSave,IG
      Else
      IGem(IO)=IGSave
      Do J=1,NBasis
      Occ(J)=OccSave(J)
      CICoef(J)=CISave(J)
      EndDo
      EndIf
C
c     endif (IG.Ne.IGSave)
      EndIf
C
C     enddo IG
      EndDo
C
C     end if(IOcc.Gt.Half) 
      EndIf
C
C     end of Do IO
      EndDo
C
      IPrint=IPrintS       
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
      DeltaE=ETot-EIni
      If(DeltaE.Gt.-1.D-6) Then
      IModG=0
      Write(6,'(2/,2X,
     $ ''*** OPTIMIZATION OF ARAI SUBSPACES TURNED OFF ***'')')
      EndIf
C
      Return
C
      End

*Deck GemSpace
      Subroutine GemSpace(Occ,NBasis,NGem)
C       
C     MODIFY ORBITAL SPACES OF APSG GEMINALS
C       
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,One=1.D0,Two=2.D0,Small=1.D-7)       
C
      Dimension Occ(NBasis),OccMin(NGem)
C       
C     FIND A MINIMAL OCCUPANCY FOR EACH GEMINAL 
C      
      Do IG=1,NGem
C
      OccMin(IG)=Two
C
      Do I=1,NBasis
      If(IGem(I).Eq.IG.And.Occ(I).Lt.OccMin(IG)) OccMin(IG)=Occ(I)
      EndDo
C
      EndDo
C
C     FIND A GEMINAL WITH THE HIGHEST MINIMAL OCCUPATION NUMBER 
C     (ITS SPACE WILL BE ENLARGED)
C
      IGMax=0
      OccMax=Zero
C
      Do IG=1,NGem
      If(OccMin(IG).Gt.OccMax) Then
      IGMax=IG
      OccMax=OccMin(IG)
      EndIf 
C
      EndDo
C
C     MOVE ALL ORBITALS THE OCCUPANCIES OF WHICH ARE SMALLER THAN SMALL 
C     TO THE IGMax GEMINAL
C
      Do I=1,NBasis
      If(Occ(I).Lt.Small) IGem(I)=IGmax
      EndDo
C  
      Return
      End

*Deck DFPMIN
      Subroutine DFPMIN(IDFP,FIJ,A,GIJ,FunIJ,Occ,T,XKin,XNuc,TNO,
     $ n,NB,N2,NGem,iter)
C
      Implicit Real*8 (A-H,O-Z)  
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension T(NB*NB)
C
C     LOCAL ARRAYS
C
      Dimension TNew(NB*NB),BIJ(NB*NB),
     $ p(n),dg(n),g(n),hdg(n),pnew(n),xi(n),xin(n)
C
      Real*8, Allocatable :: hessin(:)
C 
      Logical check
C
c herer!!!
      Parameter(ITMAX=100,EPS=3.d-8,TOLX=4.*EPS,StMax=0.1D0,GTol=5.D-4)
c      Parameter(ITMAX=100,EPS=3.d-12,TOLX=4.*EPS,StMax=0.1D0,GTol=1.D-8)
C
      Allocate (hessin(n*(n+1)/2))
C
      nh=n*(n+1)/2
C
      If(IDFP.Eq.0) Then
C
      Call E2(fp,FIJ,GIJ,T,FunIJ,BIJ,Occ,TNO,NB,N2,1,NGem)
      Call GRX(G,FIJ,A,GIJ,T,FunIJ,BIJ,Occ,TNO,NB,N2,1,NGem) 
C
      Else
C
      If(IFun.Ne.13) Then
      Call ETA(fp,A,T,Occ,FunIJ,XKin,XNuc,TNO,NB,N2) 
      Else
      Call ETAAPSG(fp,A,T,Occ,FunIJ,XKin,XNuc,TNO,NB,N2,NGem)
      EndIf
      Call GRX(G,FIJ,A,GIJ,T,FunIJ,BIJ,Occ,TNO,NB,N2,2,NGem)
C
      EndIf
C
      do 10 i=1,nh
   10 hessin(i)=0.D0

      do 12 i=1,n
        p(i)=0.D0
        hessin(i*(i+1)/2)=1.D0
        xi(i)=-g(i)
12    continue
   
      do 27 its=1,ITMAX
        iter=its

        Call LNSRCHX(IDFP,n,p,T,fp,g,xi,pnew,TNew,fret,check,
     $  FIJ,GIJ,FunIJ,BIJ,A,XKin,XNuc,Occ,TNO,NB,N2,NGem)

        if(check.eqv..false.) then
        Deallocate (hessin)
        If(IPrint.Ge.1)
     $  Write(6,'(/,X,''LNSRCHX failed, quitting DFPMIN'')')
        return
        endif

        gnorm=0.d0
        do i=1,n
        gnorm=gnorm+g(i)**2
        enddo
        gnorm=sqrt(gnorm)

        if(IDFP.eq.0) then
        tnorm=0.0
        do i=1,nb*nb
        tnorm=tnorm+t(i)**2
        enddo 
        tnorm=sqrt(tnorm)
        if(tnorm.gt.1.d0) then
        If(IPrint.Ge.1)
     $  Write(6,'(/,X,''Tnorm too large, quitting DFPMIN'')') 
        Deallocate (hessin)
        return
        endif
        endif

       If(IPrint.Ge.2)
     $   Write(6,'(''DPFMIN: Iteration:'',I4,4X,''Function'',E16.8,3X,
     $ ''Grad Norm'',E10.3)')its,fret,gnorm

        fp=fret
        do 13 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
13      continue

        test=0.
        do 14 i=1,n
          temp=abs(xi(i))/max(abs(p(i)),1.D0)
          if(temp.gt.test)test=temp
14      continue

        if(test.lt.TOLX) then
        If(IPrint.Ge.1)
     $  Write(6,'(/,X,''Convergence in DPFMIN attained after '',I3,
     $  '' iterations'')')its
        Deallocate (hessin)
        return
        endif

        do 15 i=1,n
          dg(i)=g(i)
15      continue

        Do 30 I=1,NB*NB
30      T(I)=TNew(I)

        IFG=1
        If(IDFP.Eq.1) IFG=2 
        Call GRX(G,FIJ,A,GIJ,T,FunIJ,BIJ,Occ,TNO,NB,N2,IFG,NGem)
 
        test=0.
        den=max(fret,1.D0)
        do 16 i=1,n
          temp=abs(g(i))*max(abs(p(i)),1.D0)/den
          if(temp.gt.test)test=temp
16      continue

        if(test.lt.gtol)then
        If(IPrint.Ge.1)
     $  Write(6,'(/,X,''Convergence in DPFMIN attained after '',I3,
     $  '' iterations'')')its
        Deallocate (hessin)
        return
        endif

        do 17 i=1,n
          dg(i)=g(i)-dg(i)
17      continue

        do 35 i=1,n
35      hdg(i)=0.D0

        ij=0
        do 40 i=1,n
        do 40 j=1,i
        ij=ij+1
        hij=hessin(ij)
        hdg(i)=hdg(i)+hij*dg(j)
        if(i.ne.j) hdg(j)=hdg(j)+hij*dg(i) 
40      continue

        fac=0.
        fae=0.
        sumdg=0.
        sumxi=0.
        do 21 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
21      continue
c
c calculate a new inversion of a hessian and a step
c
        if(fac**2.gt.EPS*sumdg*sumxi)then
          fac=1./fac
          fad=1./fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue

          do 43 i=1,n
          xin(i)=xi(i)
43        xi(i)=0.D0

          ij=0
          do 45 i=1,n
            do 45 j=1,i
              ij=ij+1
              hessin(ij)=hessin(ij)+fac*xin(i)*xin(j)-fad*hdg(i)*hdg(j)+
     $        fae*dg(i)*dg(j)
              xi(i)=xi(i)-hessin(ij)*g(j)
              if(i.ne.j) xi(j)=xi(j)-hessin(ij)*g(i)
45        continue

        else

        do 63 i=1,n
63        xi(i)=0.D0
          ij=0
          do 65 i=1,n
            do 65 j=1,i
              ij=ij+1
              hij=hessin(ij)
              xi(i)=xi(i)-hij*g(j)
              if(i.ne.j) xi(j)=xi(j)-hij*g(i)
65        continue

        endif
c
c     check if the step is not too large (singular hess) 
c 
      do i=1,n
      if(abs(xi(i)).gt.stmax.and.abs(g(i)).lt.gtol) then
      If(IPrint.Ge.1)
     $ Write(6,'(/,X,''Singular hess in DFPMIN - quitting...'')')
      Deallocate (hessin)
      return 
      endif
      enddo

27    continue

      If(IPrint.Ge.1)Write(6,'(/,X,''Too many iterations in DFPMIN!'')') 
C
      Deallocate (hessin)
C
      Return
      End

*Deck LNSRCHX 
      Subroutine LNSRCHX(IDFP,n,xold,TOld,fold,g,p,x,T,f,check,
     $ FIJ,GIJ,FunIJ,BIJ,A,XKin,XNuc,Occ,TNO,NB,N2,NGem)
C
      Implicit Real*8 (A-H,O-Z)      
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension TOld(NB*NB),T(NB*NB),
     $ g(n),p(n),x(n),xold(n)
C
C     LOCAL ARRAY
C
      Dimension step(n) 
C
      Logical check
C
      Parameter (ALF=1.d-4,TOLX=1.d-7, MXIT=5)
C
      check=.false.
C
C     SCALE THE STEP IF NECESSARY
C
      StepMx=0.1D0

      ScMaxR=0.0
      Do 3 I=1,n
3     ScMaxR=Max(ScMaxR,Abs(p(I)))
C
      If (ScMaxR.Ne.0.0) Then
      ScMaxR=Min(StepMx/ScMaxR,1.d0)
      Else
      ScMaxR=1.d0
      EndIf
C
      If(ScMaxR.Ne.1.d0) Then
      Do 5 I=1,n
5     p(I)=ScMaxR*p(I)
      EndIf

      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.D0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.

      it=0
1     it=it+1 

        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
          step(i)=alam*p(i)
15      continue

        If(IDFP.Eq.0) Then
C
        Call NewT(T,TOld,Step,NB,0) 
        Call E2(f,FIJ,GIJ,T,FunIJ,BIJ,Occ,TNO,NB,N2,1,NGem) 
C
        Else
C
        Call NewU(T,TOld,Step,NB,0)
        If(IFun.Ne.13) Then      
        Call ETA(f,A,T,Occ,FunIJ,XKin,XNuc,TNO,NB,N2)
        Else
        Call ETAAPSG(f,A,T,Occ,FunIJ,XKin,XNuc,TNO,NB,N2,NGem)
        EndIf
C
        EndIf

        if(alam.lt.alamin)then

          do 16 i=1,n
            x(i)=xold(i)
16        continue

          If(IDFP.Eq.0) Then
          Call NewT(T,TOld,Step,NB,1)
          Else
          Call NewU(T,TOld,Step,NB,1)
          EndIf

          Do 20 I=1,NB*NB
            T(I)=TOld(I)
20        Continue
          check=.true.
          return

        else if(f.le.fold+ALF*alam*slope)then
          If(IDFP.Eq.0) Then
          Call NewT(T,TOld,Step,NB,1)
          Else
          Call NewU(T,TOld,Step,NB,1)
          EndIf 
          check=.true.
          return

        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
            if(tmplam.gt.1.) tmplam=0.5
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)

            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              if(disc.le.0.) then
                 tmplam=.5*alam
              elseif(b.le.0.) then
                 tmplam=(-b+sqrt(disc))/(3.*a)
              else
                 tmplam=-slope/(b+sqrt(disc)) 
              endif
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=max(tmplam,.1*alam)
      if(it.lt.MXIT) goto 1

      If(IPrint.Ge.1)
     $ Write(6,'('' No of iterations exceeded in lnsrchx! '')')

      Return
      End



