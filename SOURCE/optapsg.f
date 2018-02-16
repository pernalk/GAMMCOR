*Deck OptAPSGFun
      Subroutine OptAPSGFun(ETot,URe,Occ,XKin,XNuc,
     $ TwoEl,UMOAO,NSymMO,Title,NBasis,NInte1,NInte2,MaxXV,MxHVec,
     $ NGem,NGOcc,IModG)
C
C     OPTMIZATION OF THE APSG FUNCTIONAL (ORBITALS AND EXPANSION COEFFICIENTS
C     SIMULTANEOUSLY)
C
C     GTol - ULTIMATE THRESHOLD FOR THE GRADEINT NORM
C     GTol1 - THRESHOLD FOR THE GRADIENT IN BFGSAPSG USED TILL ARAI OPTIMIZATION IS OFF
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
c herer!!!
       Parameter(GTol=5.D-6,GTol1=2.D-2,GTolGVB=1.D-1,ETol=1.D-8,
     $ MxIt=200,TolBrent=1.D-8)
   
c       Parameter(GTol=5.D-5,GTol1=2.D-2,GTolGVB=1.D-3,ETol=1.D-8,
c     $ MxIt=200,TolBrent=1.D-8)
 
c       Parameter(GTol=1.D-6,GTol1=2.D-2,GTolGVB=5.D-3,ETol=1.D-10,
c     $ MxIt=200,TolBrent=1.D-8)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(NInte1),TwoEl(NInte2),NSymMO(NBasis),UMOAO(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension Grad(MxHVec),Dir(MxHVec),NSymNO(NBasis),
     $ CoulNO(NInte1),ExchNO(NInte1),HNO(NInte1),
     $ IndX(MaxXV),IndN(2,MaxXV),IndAux(NBasis),IndC(NBasis)
C
C     IF IDALTON=1 CHECK THE GRADIENT AND RETURN
C
      If(IDALTON.Eq.1) Then 
C
      Write(6,'(/,X,''  ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,I6)') I,Occ(I),IGem(I)
      EndDo
C
      If(ICASSCF.Eq.1) Return
C
      Call GrAPSG(ETot,Grad,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $ NInte1,NInte2,MaxXV,MxHVec,NGem,1)
      GradNo=Zero
      Do I=1,MxHVec
      GradNo=GradNo+Grad(I)*Grad(I)
      EndDo
      GradNo=Sqrt(GradNo)
C
      Write(6,'(/,X,''  DATA READ FROM DALTON OUTPUT. Electronic ETot = 
     $ '',E16.6)')ETot
      Write(6,'(/,X,''  DATA READ FROM DALTON OUTPUT. GRAD NORM = 
     $ '',E16.6)')GradNo
C
      If(GradNo.Lt.1.D-5) Return
C
      Stop 'DALTON GVB CALCULATIONS HAS NOT BEEN CONVERGED. EXIT.'
C
c     If(IDALTON.Eq.1) 
      EndIf
C
C     CONSTRUCT LOOK-UP TABLES TO REDUCE DIMENSIONS OF GRADIENTS AND HESSIANS
C
      Do I=1,NInActOrb
      IndAux(I)=0
      EndDo
      Do I=1+NInActOrb,NInActOrb+NActOrb
      IndAux(I)=1
      EndDo
      Do I=1+NInActOrb+NActOrb,NBasis
      IndAux(I)=2
      EndDo
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
C     EXCLUDE PAIRS OF TYPE 0-0, 2-2
C
      ISum=IndAux(I)+IndAux(J)
      If(ISum.Ne.0.And.ISum.Ne.4) Then
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
C
C     EXCLUDE ROTATIONS OF FROZEN ORBITALS
C
      If((IFreeze.Eq.1).And.(I.Le.NGOcc.Or.J.Le.NGOcc)) Ind=Ind-1
C
      EndDo
      EndDo
C
      NDimX=Ind
C
      Write(6,'(2/,X,"Total number of X parameters:",I4)')MaxXV
      Write(6,'(X,"Reduced to:",I4)')NDimX
C
C     IF IGVB IS REQUESTED FIND THE NUMBER OF C-COEFFICIENTS TO BE OPTIMIZED
C     IT IS ONLY USED IN GVB CASE in Newton
C
      NDimN=NBasis
      Do I=1,NBasis
      IndC(I)=I
      EndDo
C
      If(IGVB.Eq.1) Then
C
      NDimN=2*NELE
C
      If(IFreeze.Eq.1) NDimN=NDimN-2*NGOcc
C
      Ind=0
      Do I=1,NBasis
C
      If(IFreeze.Eq.1) Then
C
      If(I.Gt.NGOcc.And.I.Le.2*NELE-NGOcc) Then
      Ind=Ind+1
      IndC(Ind)=I
      EndIf
C
      Else
C
      If(I.Le.2*NELE) Then
      Ind=Ind+1
      IndC(Ind)=I
      EndIf
C
      EndIf
C
      EndDo
C
      If(Ind.Ne.NDimN) Stop 'Fatal error: incorrect NDimN'
C
C     endif of IGVB.Eq.1
      EndIf
C
c herer!!!
C
C     FOR TESTING OF GVB ONLY
C     ICrudeGuess = 1  - guess orbitals from crude APSG optimization (with GTolGVB, without Arai opt)
C     ICrudeGuess = 0  - guess orbitals from full APSG optimization (with Arai) with GTol1
      ICrudeGuess=1
C
C     Establish symmetries of NO's
C
      Do I=1,NBasis
      NSymNO(I)=1
      EndDo
C
      If(NoSym.Eq.0) Then
C
      Do I=1,NBasis
C
      NSymNO(I)=0
C
      Do J=1,NBasis
C
      If(Abs(URe(I,J)).Ne.Zero) Then
C
      ISym=NSymMO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym)
     $ Write(*,*)'Symm of NO',I,' cannot be established'
      EndIf
C
      EndIf
      EndDo
      EndDo
C
      EndIf
C
      Iter=0
C
C     BEGINNING OF THE ITERATIONS
C
c herer!!!
      If(ICrudeGuess.Eq.1.And.IGVB.Eq.1) IModG=0

   10 Iter=Iter+1
C
      ETotO=ETot
C
   22 GTolB=GTol
      If(IModG.Ne.0) GTolB=GTol1
c herer!!!
      If(Iter.Eq.1.And.IGVB.Eq.1.And.ICrudeGuess.Eq.1) GTolB=GTolGVB
C
C     BEGIN GVB CALCULATIONS
C
c herer!!!
      If((Iter.Eq.2.And.IGVB.Eq.1.And.ICrudeGuess.Eq.1).Or.
     $ (IGVB.Eq.1.And.ICrudeGuess.Eq.0.And.IModG.Eq.0).Or.
     $ (IGVB.Eq.1.And.IRes.Eq.1.And.Iter.Eq.1) ) Then
C
      ETotO=Zero
      IModG=0
      Call SortAll(Occ,URe,NSymNO,NBasis)
      Call GVBIni(Occ,NBasis,NGem,NGOcc) 
      Call SortAll(Occ,URe,NSymNO,NBasis)
C
      GTolB=GTol
C
      EndIf
C
      NN=NDimX+NBasis
C
      Call BFGSAPSG(ETot,Grad,Occ,URe,XKin,XNuc,TwoEl,NSymNO,NSymMO,
     $ IndX,IndN,NDimX,
     $ NN,NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem,GTolB,MxIt,IConv)
C
C     ARAI SPACE OPTIMIZATION
C
      If(IModG.Ne.0) Then
C
      Call TrTwoEl(CoulNO,ExchNO,URe,TwoEl,NBasis,NInte2)
      Call OptArai(ETot,Occ,URe,XKin,XNuc,CoulNO,ExchNO,
     $ NSymNO,NBasis,NInte1,NGem,IModG,NGOcc)
      Call GrAPSG(ETot,Grad,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $ NInte1,NInte2,MaxXV,MxHVec,NGem,1) 
C
C     IF ARAI IS OFF RETURN TO BFGS WITH TIGHTER CONV CRIT
C
      If(IModG.Eq.0) Then
      Iter=Iter+1
      Goto 22
      EndIf
C
      EndIf
C      
C     CALCULATE THE GRADIENT NORM
C
      GradNo=Zero
C
      Do 20 I=1,MxHVec
   20 GradNo=GradNo+Grad(I)*Grad(I)
      GradNo=Sqrt(GradNo)
C
      Call SortAll(Occ,URe,NSymNO,NBasis)
C
      Call ReWrAPSG(1,Occ,URe,UMOAO,Title,NGem,NBasis)
C
      If(IPrint.Ge.0)
     $ Write(6,'(/,X,''ITER'',I3,2X,''ENERGY'',F16.8,2X,
     $ ''ENE DIFF '',E10.3,2X,''GRAD NORM '',E9.3)')
     $ Iter,ETot,ETot-ETotO,GradNo
      Write(6,'(/,X,''  ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem",3X,"Sym")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I),NSymNO(I)
      EndDo
C
      If(IModG.Eq.1.Or.
     $ ((GradNo.Gt.GTol).And.(Iter.Lt.MxIt).And.ETotO-ETot.Gt.ETol)) 
     $ GoTo 10 
C
      If(Iter.Lt.MxIt) Then
C
      If(GradNo.Lt.GTol) Then
      Write(6,'(/,10X,''TOTAL CONVERGENCE ATTAINED '')')
      Else
      Write(6,'(/,10X,''ENERGY CONVERGED BUT GRADIENT TOO LARGE '')')
      EndIf 
C
      Else
C
      Write(6,'(/,10X,''NO CONVERGENCE ATTAINED '')')
C
      EndIf
C
      Return
      End

*Deck LinMinAPSG
      Subroutine LinMinAPSG(ETot,Dir,URe,Occ,XKin,XNuc,TwoEl,NSymNO,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem,TolBrent)
C
C     LINEAR MINIMIZATION USING BRENT METHOD ALONG Dir
C     ON EXIT URe, Occ, AND CICoef ARE UPDATED
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab 
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
      Dimension Dir(MxHVec),Occ(NBasis),URe(NBasis,NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),NSymNO(NBasis) 
C
C     LOCAL ARRAYS
C
      Dimension Step(MxHVec),UReOld(NBasis,NBasis),CoefOld(NBasis)
C
      Do I=1,NBasis
      CoefOld(I)=CICoef(I)
      Do J=1,NBasis
      UReOld(I,J)=URe(I,J)
      EndDo
      EndDo
C
       Call BrackAPSG(AX,BX,CX,FA,FB,FC, 
     $ Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C   
      ETot=BrentAPSG(AX,BX,CX,TolBrent,XMin,
     $ Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
      Do I=1,MxHVec
      Step(I)=XMin*Dir(I)
      EndDo
C
      Call NewUAPSG(URe,UReOld,Step,NBasis,ICheckU)
      Call NewC(CoefOld,Step(MaxXV+1),NBasis,NGem)
      Do I=1,NBasis
      Occ(I)=CICoef(I)**2
      EndDo
C
      Return 
      End

*Deck GrAPSG
      Subroutine GrAPSG(ETot,Grad,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $NInte1,NInte2,MaxXV,MxHVec,NGem,IFlagG)
C
C     CALCULATE THE GRADIENTS IN NO OF TOTAL APSG ENERGY WITH RESPECT TO 
C     ORBITAL- AND EXPANSION COEFFICIENT- PARAMETERS
C     
C     FLAGS:  IFlagG = 0 ... calculate energy only
C                      1 ... calculate energies and gradients
C
C     Tol2 - Threshold for two-el integrals 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0,
     $ Tol2=1.D-15)
C
      Dimension Grad(MxHVec),URe(NBasis,NBasis),Occ(NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension
     $ Gamma(NInte1),Gam(NGem,NInte1),GamC(NGem,NInte1),
     $ F(NGem*NInte1),GC(NGem*NInte1)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     CONSTRUCT A DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
C     AND AUXILIARY DENSITY-LIKE MATRICES
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
      NAdrIL=(IG-1)*NInte1+IL
      F(NAdrIL)=F(NAdrIL)-Half*FacIL*Gam(IG,JK)*TwoZet
      GC(NAdrIL)=GC(NAdrIL)+Half*FacIL*GamC(IG,JK)*TwoZet
C
      If(IL.Ne.JK) Then
      NAdrJK=(IG-1)*NInte1+JK
      F(NAdrJK)=F(NAdrJK)-Half*FacJK*Gam(IG,IL)*TwoZet
      GC(NAdrJK)=GC(NAdrJK)+Half*FacJK*GamC(IG,IL)*TwoZet
      EndIf
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
      NAdrIK=(IG-1)*NInte1+IK
      F(NAdrIK)=F(NAdrIK)-Half*FacIK*Gam(IG,JL)*TwoZet
      GC(NAdrIK)=GC(NAdrIK)+Half*FacIK*GamC(IG,JL)*TwoZet
C 
      If(IK.Ne.JL) Then
      NAdrJL=(IG-1)*NInte1+JL
      F(NAdrJL)=F(NAdrJL)-Half*FacJL*Gam(IG,IK)*TwoZet
      GC(NAdrJL)=GC(NAdrJL)+Half*FacJL*GamC(IG,IK)*TwoZet
      EndIf
C
      EndDo
C                                                  
   40 Continue
   30 Continue
C
C     TRANSFORM F, AND GC TO NO
C
      Do IG=1,NGem
C
      Call MatTr(F((IG-1)*NInte1+1),URe,NBasis)
      Call MatTr(GC((IG-1)*NInte1+1),URe,NBasis)
C
      EndDo
C
C     COMPUTE THE ENERGY
C
      ETot=Zero
C
      IJ=0
      Do I=1,NBasis
      IG=IGem(I)
      ETot=ETot+Occ(I)*F((IG-1)*NInte1+I*(I+1)/2)
     $ +CICoef(I)*GC((IG-1)*NInte1+I*(I+1)/2)
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
C     ------------------------------------
C     AT THIS POINT THE ENERGIES ARE KNOWN
C     ------------------------------------
C
      If (IFlagG.Eq.0) Return
C
C     CONSTRUCT THE GRADIENTS WITH RESPECT TO X
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ2=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
C
      NAdr=(IGem(I)-1)*NInte1+IJ2
      Grad(IJ)=Four*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      NAdr=(IGem(J)-1)*NInte1+IJ2
      Grad(IJ)=Grad(IJ)-Four*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
C
      EndDo
      EndDo
C
C     CONSTRUCT THE GRADIENTS WITH RESPECT TO x PARAMETERS 
C     RELATED TO EXPANSION COEFFICIENTS c
C
      Do I=1,NBasis
C
      NAdr=(IGem(I)-1)*NInte1+(I*(I+1))/2
      Grad(I+MaxXV)=Four*CICoef(I)*F(NAdr)+Two*GC(NAdr)
C
      Do J=1,NBasis
      If(IGem(I).Eq.IGem(J)) Then
      NAdr=(IGem(J)-1)*NInte1+J*(J+1)/2
      Grad(I+MaxXV)=Grad(I+MaxXV)-Two*CICoef(I)*
     $ (Two*Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      EndIf
      EndDo
C
      EndDo
C
C     -------------------------------------
C     AT THIS POINT THE GRADIENTS ARE KNOWN
C     -------------------------------------
C
C     IF IFreeze=1 SET GRADIENT ELEMENTS CORRESPONDING TO FROZEN ORBS TO ZERO
C
      If(IFreeze.Eq.1) Then
C
      Do IOrb=1,NBasis 
C
      NNGem=0
      Do J=1,NBasis
      If(IGem(J).Eq.IGem(IOrb)) NNGem=NNGem+1
      EndDo
C
C     IF AN ORBITAL IS FULLY OCCUPIED AND A CORRESPONDING GEMINAL 
C     INCLUDES ONLY ONE ORBITAL SET GRAD TO ZERO
C 
      If(Occ(IOrb).Eq.One.And.NNGem.Eq.1) Then
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      If(I.Eq.IOrb.Or.J.Eq.IOrb) Grad(IJ)=Zero
      EndDo
      EndDo
      EndIf
C
      EndDo
C
      EndIf
C
      Return
      End

*Deck HessAPSG
      Subroutine HessAPSG(ETot,Grad,Hess,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $NInte1,NInte2,MaxXV,MxHVec,NGem,IFlagG)
C
C     CALCULATE THE GRADIENTS AND A HESSIAN IN NO OF TOTAL APSG ENERGY WITH RESPECT TO 
C     ORBITAL- AND EXPANSION COEFFICIENT- PARAMETERS
C     
C     FLAGS:  IFlagG = 0 ... calculate energy only
C                      1 ... calculate energy and gradients
C                      2 ... calculate energy, grad, and hess
C
C     Tol2 - Threshold for two-el integrals 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Three=3.D0,
     $ Four=4.D0,Eight=8.0D0,Tol2=1.D-15)
C
      Dimension Grad(MxHVec),URe(NBasis,NBasis),Occ(NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $ Hess(MxHVec*(MxHVec+1)/2)
C
C     LOCAL ARRAYS
C
      Dimension
     $ Gamma(NInte1),Gam(NGem,NInte1),GamC(NGem,NInte1),
     $ F(NGem*NInte1),GC(NGem*NInte1),HNO(NInte1)
C
      Real*8, Allocatable :: TNO(:)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     CONSTRUCT A DENSITY MATRIX
C
      Call VecTr(Gamma,Occ,URe,NBasis)
C
C     IF Hess IS REQUESTED TRANSFORM TO NO AND FIND AUX MATRICES IN NO
C
      If(IFlagG.Eq.2) Then
C
      Allocate (TNO(NInte2))
      Call TwoNO(TNO,URe,TwoEl,NBasis,NInte2)
C
      Do I=1,NInte1
      HNO(I)=XKin(I)+XNuc(I)
      EndDo
      Call MatTr(HNO,URe,NBasis)
C
      Do IG=1,NGem
      Do I=1,NInte1
      F((IG-1)*NInte1+I)=HNO(I)
      GC((IG-1)*NInte1+I)=Zero
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      Do K=1,NBasis
C
      TwoCoul=TNO(NAddr3(I,J,K,K))
      TwoEx=TNO(NAddr3(I,K,J,K))
C
      Do IG=1,NGem
C
      NAd=(IG-1)*NInte1+IJ
C
      If(IGem(K).Eq.IG) Then
      GC(NAd)=GC(NAd)+CICoef(K)*TwoEx
      Else
      F(NAd)=F(NAd)+Occ(K)*(Two*TwoCoul-TwoEx)
      EndIf
C
      EndDo
      EndDo 
C
      EndDo
      EndDo
C
      GoTo 777    
C
C     endif of IFlag.Eq.2
      EndIf
C
C     AND AUXILIARY DENSITY-LIKE MATRICES
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
      NAdrIL=(IG-1)*NInte1+IL
      F(NAdrIL)=F(NAdrIL)-Half*FacIL*Gam(IG,JK)*TwoZet
      GC(NAdrIL)=GC(NAdrIL)+Half*FacIL*GamC(IG,JK)*TwoZet
C
      If(IL.Ne.JK) Then
      NAdrJK=(IG-1)*NInte1+JK
      F(NAdrJK)=F(NAdrJK)-Half*FacJK*Gam(IG,IL)*TwoZet
      GC(NAdrJK)=GC(NAdrJK)+Half*FacJK*GamC(IG,IL)*TwoZet
      EndIf
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
      NAdrIK=(IG-1)*NInte1+IK
      F(NAdrIK)=F(NAdrIK)-Half*FacIK*Gam(IG,JL)*TwoZet
      GC(NAdrIK)=GC(NAdrIK)+Half*FacIK*GamC(IG,JL)*TwoZet
C 
      If(IK.Ne.JL) Then
      NAdrJL=(IG-1)*NInte1+JL
      F(NAdrJL)=F(NAdrJL)-Half*FacJL*Gam(IG,IK)*TwoZet
      GC(NAdrJL)=GC(NAdrJL)+Half*FacJL*GamC(IG,IK)*TwoZet
      EndIf
C
      EndDo
C                                                  
   40 Continue
   30 Continue
C
C     TRANSFORM F, AND GC TO NO
C
      Do IG=1,NGem
C
      Call MatTr(F((IG-1)*NInte1+1),URe,NBasis)
      Call MatTr(GC((IG-1)*NInte1+1),URe,NBasis)
C
      EndDo
C
  777 Continue
C
C     COMPUTE THE ENERGY
C
      ETot=Zero
C
      IJ=0
      Do I=1,NBasis
      IG=IGem(I)
      ETot=ETot+Occ(I)*F((IG-1)*NInte1+I*(I+1)/2)
     $ +CICoef(I)*GC((IG-1)*NInte1+I*(I+1)/2)
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
C     ------------------------------------
C     AT THIS POINT THE ENERGIES ARE KNOWN
C     ------------------------------------
C
      If (IFlagG.Eq.0) Return
C
C     CONSTRUCT THE GRADIENTS WITH RESPECT TO X
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
      IJ2=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)
C
      NAdr=(IGem(I)-1)*NInte1+IJ2
      Grad(IJ)=Four*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      NAdr=(IGem(J)-1)*NInte1+IJ2
      Grad(IJ)=Grad(IJ)-Four*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
C
      EndDo
      EndDo
C
C     CONSTRUCT THE GRADIENTS WITH RESPECT TO x PARAMETERS 
C     RELATED TO EXPANSION COEFFICIENTS c
C
      Do I=1,NBasis
C
      NAdr=(IGem(I)-1)*NInte1+(I*(I+1))/2
      Grad(I+MaxXV)=Four*CICoef(I)*F(NAdr)+Two*GC(NAdr)
C
      Do J=1,NBasis
      If(IGem(I).Eq.IGem(J)) Then
      NAdr=(IGem(J)-1)*NInte1+J*(J+1)/2
      Grad(I+MaxXV)=Grad(I+MaxXV)-Two*CICoef(I)*
     $ (Two*Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      EndIf
      EndDo
C
      EndDo
C
C     -------------------------------------
C     AT THIS POINT THE GRADIENTS ARE KNOWN
C     -------------------------------------
C
      If (IFlagG.Eq.1) Return
C
      Do I=1,MxHVec*(MxHVec+1)/2
      Hess(I)=Zero
      EndDo
C
C     COMPUTE THE Hessian_XX BLOCK OF THE HESSIAN 
C
      IJKL=0
C
      IJ=0
      Do I=1,NBasis      
      Do J=1,I-1
      IJ=IJ+1
C
      KL=0
      Do K=1,NBasis 
      Do L=1,K-1
      KL=KL+1
C
      If(IJ.Ge.KL) Then
C
      IJKL=IJKL+1
      Hess(IJKL)=Zero
C
C     PART 1           
C
      Aux=Zero 
      If(IGem(I).Ne.IGem(K)) Aux=Aux+Occ(I)*Occ(K)
      If(IGem(J).Ne.IGem(K)) Aux=Aux-Occ(J)*Occ(K)
      If(IGem(I).Ne.IGem(L)) Aux=Aux-Occ(I)*Occ(L)
      If(IGem(J).Ne.IGem(L)) Aux=Aux+Occ(J)*Occ(L)
      Aux=Four*Aux*(Four*TNO(NAddr3(I,J,K,L))
     $              -TNO(NAddr3(I,L,K,J))-TNO(NAddr3(I,K,L,J)))
C
      Aux1=Zero
      If(IGem(I).Eq.IGem(K)) Aux1=Aux1+CICoef(I)*CICoef(K)
      If(IGem(J).Eq.IGem(K)) Aux1=Aux1-CICoef(J)*CICoef(K)
      If(IGem(I).Eq.IGem(L)) Aux1=Aux1-CICoef(I)*CICoef(L)
      If(IGem(J).Eq.IGem(L)) Aux1=Aux1+CICoef(J)*CICoef(L)
      Aux1=Four*Aux1*(TNO(NAddr3(I,L,K,J))+TNO(NAddr3(I,K,L,J)))
C
      Hess(IJKL)=Aux+Aux1
C
C     PART 2
C
C     delta_jk 
      If(J.Eq.K) Then
      IL=(Max(I,L)*(Max(I,L)-1))/2+Min(I,L)
      NAdr=(IGem(I)-1)*NInte1+IL
      Hess(IJKL)=Hess(IJKL)+Two*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      NAdr=(IGem(L)-1)*NInte1+IL
      Hess(IJKL)=Hess(IJKL)+Two*(Occ(L)*F(NAdr)+CICoef(L)*GC(NAdr))
      NAdr=(IGem(J)-1)*NInte1+IL
      Hess(IJKL)=Hess(IJKL)-Four*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      EndIf
C
C     delta_jl
      If(J.Eq.L) Then
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      NAdr=(IGem(I)-1)*NInte1+IK
      Hess(IJKL)=Hess(IJKL)-Two*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      NAdr=(IGem(K)-1)*NInte1+IK
      Hess(IJKL)=Hess(IJKL)-Two*(Occ(K)*F(NAdr)+CICoef(K)*GC(NAdr))
      NAdr=(IGem(J)-1)*NInte1+IK
      Hess(IJKL)=Hess(IJKL)+Four*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      EndIf
C
C     delta_il
      If(I.Eq.L) Then
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      NAdr=(IGem(J)-1)*NInte1+JK
      Hess(IJKL)=Hess(IJKL)+Two*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      NAdr=(IGem(K)-1)*NInte1+JK
      Hess(IJKL)=Hess(IJKL)+Two*(Occ(K)*F(NAdr)+CICoef(K)*GC(NAdr))
      NAdr=(IGem(I)-1)*NInte1+JK
      Hess(IJKL)=Hess(IJKL)-Four*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      EndIf
C
C     delta_ik
      If(I.Eq.K) Then
      JL=(Max(J,L)*(Max(J,L)-1))/2+Min(J,L)
      NAdr=(IGem(J)-1)*NInte1+JL
      Hess(IJKL)=Hess(IJKL)-Two*(Occ(J)*F(NAdr)+CICoef(J)*GC(NAdr))
      NAdr=(IGem(L)-1)*NInte1+JL
      Hess(IJKL)=Hess(IJKL)-Two*(Occ(L)*F(NAdr)+CICoef(L)*GC(NAdr))
      NAdr=(IGem(I)-1)*NInte1+JL
      Hess(IJKL)=Hess(IJKL)+Four*(Occ(I)*F(NAdr)+CICoef(I)*GC(NAdr))
      EndIf
C
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C
C     COMPUTE THE Hessian_Xx BLOCK OF THE HESSIAN
C
      Do I=1,NBasis
      Do J=1,NBasis
C
      If(I.Ne.J) Then
C
      Ind1=(Max(I,J)*(Max(I,J)-1))/2+Min(I,J)-Max(I,J)+1
C
      Do K=1,NBasis
C
      Ind2=MaxXV+K
      IJK=(Ind2*(Ind2-1))/2+Ind1
C
      Aijxk=Zero
C
      If(IGem(I).Eq.IGem(K)) Then
C
      NAdr=(IGem(I)-1)*NInte1+
     $ (Max(J,I)*(Max(J,I)-1))/2+Min(J,I)
      Aux=Four*(Two*CICoef(I)*F(NAdr)+GC(NAdr))
      Aux1=-CICoef(I)*CICoef(K)
      If(I.Eq.K) Aux1=Aux1+One
      Aijxk=Aux*Aux1
C
      EndIf
C
      Sum1=Zero
      Sum2=Zero
      Do L=1,NBasis
C 
      Aux=-CICoef(L)*CICoef(K)
      If(L.Eq.K) Aux=Aux+One  
C
      If(IGem(L).Ne.IGem(I).And.IGem(L).Eq.IGem(K)) Then
      Sum1=Sum1+
     $CICoef(L)*(Two*TNO(NAddr3(I,J,L,L))-TNO(NAddr3(I,L,J,L)))*Aux 
      EndIf
C
      If(IGem(L).Eq.IGem(I).And.IGem(L).Eq.IGem(K)) Then
      Sum2=Sum2+TNO(NAddr3(I,L,J,L))*Aux
      EndIf
C
      EndDo      
C
      Aijxk=Aijxk+Eight*Occ(I)*Sum1+Four*CICoef(I)*Sum2
C
      If(I.Gt.J) Then
      Hess(IJK)=Hess(IJK)+Aijxk
      Else
      Hess(IJK)=Hess(IJK)-Aijxk
      EndIf
C
      EndDo
C
C     EndIf I.Ne.J
      EndIf
C
      EndDo
      EndDo
C
C     COMPUTE THE Hessian_xx BLOCK OF THE HESSIAN
C
      Do I=1,NBasis
      Ind1=MaxXV+I  
C
      Do J=1,I
      Ind2=MaxXV+J
      IJ=(Ind1*(Ind1-1))/2+Ind2
      Hess(IJ)=Zero
C
      SumKL=Zero
      SumK=Zero
C
      Do K=1,NBasis
C
      If(IGem(K).Eq.IGem(I).Or.IGem(K).Eq.IGem(J)) Then
C
      Aux=Zero
      If(K.Eq.I.And.IGem(K).Eq.IGem(J)) Aux=-CICoef(J)
      If(K.Eq.J.And.IGem(K).Eq.IGem(I)) Aux=Aux-CICoef(I)  
      If(I.Eq.J.And.IGem(K).Eq.IGem(I)) Aux=Aux-CICoef(K)
      If(IGem(K).Eq.IGem(I).And.IGem(K).Eq.IGem(J)) 
     $ Aux=Aux+Three*CICoef(I)*CICoef(J)*CICoef(K)
C
      NAdr=(IGem(K)-1)*NInte1+K*(K-1)/2+K
      SumK=SumK+Aux*(Four*CICoef(K)*F(NAdr)+Two*GC(NAdr))
C
C     endif of IGem(K).Eq.IGem(I).Or.IGem(K).Eq.IGem(J)
      EndIf 
C
      Do L=1,NBasis
C
      If(IGem(J).Eq.IGem(L).And.IGem(I).Eq.IGem(K)) Then 
C
      Aux=Zero
      If(K.Eq.L) Then
      NAdr=(IGem(K)-1)*NInte1+K*(K-1)/2+K
      Aux=Four*F(NAdr)
      EndIf
C
      If(IGem(K).Eq.IGem(L)) Then
      Aux=Aux+Two*TNO(NAddr3(K,L,K,L))
      Else
      Aux=Aux+Eight*CICoef(K)*CICoef(L)*
     $ (Two*TNO(NAddr3(K,K,L,L))-TNO(NAddr3(K,L,K,L)))
      EndIf
C
      Aux1=-CICoef(J)*CICoef(L)
      If(J.Eq.L) Aux1=Aux1+One
C
      Aux2=-CICoef(K)*CICoef(I)
      If(K.Eq.I) Aux2=Aux2+One
C
      SumKL=SumKL+Aux*Aux1*Aux2
C
      EndIf
C
      EndDo
      EndDo
C
      Hess(IJ)=Hess(IJ)+SumK+SumKL 
C
      EndDo
      EndDo
C
      Deallocate (TNO)
C
      Return
      End

*Deck EDirAPSG
      Real*8 Function EDirAPSG(X,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
C     CALCULATE ETot AT A POINT URe,CICoef moved by Step
C     CICoef AND URe ARE NOT MODIFFIED
C     
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C     
      Dimension Dir(MxHVec),URe(NBasis,NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C     
C     LOCAL ARRAYS
C     
      Dimension Grad(MxHVec),Occ(NBasis),UReNew(NBasis,NBasis),
     $ Step(MxHVec),CoefOld(NBasis)
C
      Do I=1,NBasis
      CoefOld(I)=CICoef(I)
      EndDo
C
      Do I=1,MxHVec
      Step(I)=X*Dir(I)
      EndDo
C
      Call NewUAPSG(UReNew,URe,Step,NBasis,ICheckU)
      Call NewC(CoefOld,Step(MaxXV+1),NBasis,NGem)
      Do I=1,NBasis
      Occ(I)=CICoef(I)**2
      EndDo
C
      Call GrAPSG(ETot,Grad,UReNew,Occ,XKin,XNuc,TwoEl,NBasis,
     $ NInte1,NInte2,MaxXV,MxHVec,NGem,0)
C
C     RESTORE INITIAL CICoef
C
      Do I=1,NBasis
      CICoef(I)=CoefOld(I)
      EndDo
C
      EDirAPSG=ETot
C
      Return 
      End

*Deck BrentAPSG
      Real*8 Function BrentAPSG(ax,bx,cx,tol,xmin,
     $ Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
      Implicit Real*8 (A-H,O-Z)
C
C     BRENT MINIMIZATION ALONG Dir
C     URe AND CICoef ARE NOT MODIFIED
C     THE OPTIMAL POINT RETURNED IN xmin
C
      PARAMETER (ITMAX=100,CGOLD=.3819660D0,ZEPS=1.0D-10)
C 
      Dimension Dir(MxHVec),URe(NBasis,NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
C
      fx=EDirAPSG(x,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
      eold=fx
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
        fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
        if(abs(fu-eold).lt.tol) goto 3
        eold=fu           
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
      write(*,*)'brentapsg exceed maximum iterations'
3     xmin=x
      brentapsg=fx

      write
     $ (6,'(/,"Brent : ",I3," iterations. Final energy: ",E16.8)')
     $ iter,fx
      return
      END

*Deck BrackAPSG
      Subroutine BrackAPSG(ax,bx,cx,fa,fb,fc,
     $ Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
      Implicit Real*8 (A-H,O-Z)
C
      PARAMETER (GOLD=1.618034D0, GLIMIT=100.D0, TINY=1.d-20)
      PARAMETER (ZERO=0.0D0, TEN5=0.05D0, MXIT=20)
C
      Dimension Dir(MxHVec),URe(NBasis,NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C
      ax=ZERO
      bx=TEN5
C
      fa=EDirAPSG(ax,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem) 
      fb=EDirAPSG(bx,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
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
      fc=EDirAPSG(cx,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
      its=0
    1 its=its+1
C
      if(its.gt.MXIT) Then
      write(*,*) 'brackapsg exceeding no of iterations!'
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
          fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
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
          fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
C
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)

          endif
C
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
        else
          u=cx+GOLD*(cx-bx)
          fu=EDirAPSG(u,Dir,URe,XKin,XNuc,TwoEl,
     $ NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem)
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

*Deck NewC
      Subroutine NewC(CoefOld,Step,NBasis,NGem)
C
C     CALCULATE CICoef FROM CoefOld AND STEP
C     
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C     
      Dimension Step(NBasis),CoefOld(NBasis)
C
C     LOCAL ARRAYS
C
      Dimension DenomG(NGem)
C
      Do IG=1,NGem
C
      Sum1=Zero
      Sum2=Zero
C
      Do I=1,NBasis
      If(IG.Eq.IGem(I)) Then
      Sum1=Sum1+CoefOld(I)*Step(I)
      Sum2=Sum2+Step(I)**2
      EndIf 
      EndDo
C
      DenomG(IG)=SQRT(Abs(One+Two*Sum1+Sum2))
C
      EndDo
C
      Do I=1,NBasis
      CICoef(I)=(CoefOld(I)+Step(I))/DenomG(IGem(I))
      EndDo
C
      Return 
      End

*Deck BFGSAPSG
      Subroutine BFGSAPSG(fret,GradNew,Occ,URe,XKin,XNuc,TwoEl,NSymNO,
     $ NSymMO,
     $ IndX,IndN,NDimX,
     $ n,NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem,GTol,MxIt,IConv)
C
      Implicit Real*8 (A-H,O-Z)  
C
C     BFGS MINIMIZATION OF THE APSG FUNCTIONAL
C
C     n = a total naumber of parameters to be optimized
C 
C     GTol - convergence criterium for the gradient       
C
C     POSSIBLE CASES
C 
C     NDimX > 0 and n > NDimX    ... vary X and x
C     NDimX > 0 and n = NDimX    ... vary X only (minimization w.r.t. orb rotations)
C     NDimX = 0 and n = NBasis   ... vary x only (minimization w.r.t. expansion coeffs c)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0)
C    
      Dimension URe(NBasis*NBasis),Occ(NBasis),GradNew(MxHVec),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $ NSymNO(NBasis),NSymMO(NBasis), 
     $ IndX(NDimX),IndN(2,NDimX)
C       
C     LOCAL ARRAYS
C       
      Dimension UReNew(NBasis*NBasis),OccNew(NBasis),CoefOld(NBasis),
     $ p(n),dg(n),g(n),hdg(n),pnew(n),xi(n),xin(n)
c herer!!!
c     $ ,dir(MxHVec),HlpU(NBasis*NBasis),HlpO(NBasis)
C       
      Real*8, Allocatable :: hessin(:)
C 
      Logical check
C       
      Parameter(EPS=3.d-8)
C       
      Allocate (hessin(n*(n+1)/2))
C       
      nh=n*(n+1)/2
      IConv=0
C
C     Find initial energy and gradient
C
      Call GrAPSG(fp,GradNew,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $      NInte1,NInte2,MaxXV,MxHVec,NGem,1)
C
C     CHECK THE NORM
C
        gnorm=Zero
        do i=1,MxHVec
        gnorm=gnorm+GradNew(i)**2
        enddo
        gnorm=sqrt(gnorm)
C
        if(gnorm.lt.GTol) then
        Write(6,'(/,X,''Convergence in BFGS APSG attained after '',I3,
     $ '' iterations'')')0
        Deallocate (hessin)
        IConv=1
        fret=fp
        return
        endif

C
      Do I=1,n
      If(I.Le.NDimX) g(I)=GradNew(IndX(I))
      If(I.Gt.NDimX) g(I)=GradNew(MaxXV+I-NDimX)
      EndDo
C
C     Initialize the inverse hessian to the unit matrix 
C     the initial point to 0, and the initial direction 
C     to the negative grad
C
      do 10 i=1,nh
   10 hessin(i)=Zero

      do 12 i=1,n
        p(i)=Zero
        hessin(i*(i+1)/2)=One
        xi(i)=-g(i)
12    continue
C
C     Main loop over iterations
C
      do 27 its=1,MxIt
C
C       Symmetry-block orbital rotations
C
        If(NoSym.Eq.0) Then 
        Do I=1,NDimX
        IndN(1,I)=IA
        IndN(2,I)=IB
        If (NSymNO(IA).Ne.NSymNO(IB)) xi(I)=Zero
        EndDo
        EndIf
C
        Do I=1,NBasis
        CoefOld(I)=CICoef(I)
        EndDo 
C
        Call LinAPSG(n,p,URe,CoefOld,fp,g,xi,pnew,UReNew,OccNew,fret,
     $  GradNew,check,XKin,XNuc,TwoEl,NBasis,NInte1,NInte2,MaxXV,NGem,
     $  IndX,IndN,NDimX,MxHVec)
C
C       Save the old gradient and get the new one
C
        do 15 i=1,n
          dg(i)=g(i)
          If(I.Le.NDimX) g(I)=GradNew(IndX(I))
          If(I.Gt.NDimX) g(I)=GradNew(MaxXV+I-NDimX)
15      continue
C
        gnorm=Zero
        do i=1,n
        gnorm=gnorm+g(i)**2
        enddo
        gnorm=sqrt(gnorm)

        If(IPrint.Ge.0)
     $   Write(6,'(''BFGS APSG: Iteration:'',I4,4X,''Energy'',E16.8,3X,
     $   ''Grad Norm'',E10.3)')its,fret,gnorm
        If(IPrint.Ge.4) Then
         Write(6,'(/,X,''  ORBITAL OCCUPANCIES '')')
         Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem",3X,"Sym")')
         Do I=1,NBasis
         Write(6,'(X,I3,E16.6,2I6)') I,OccNew(I),IGem(I),NSymNO(I)
         EndDo
        EndIf
C
C       Update the line direction and the current point
C
        fp=fret
        do 13 i=1,n
          xi(i)=pnew(i)-p(i)
          p(i)=pnew(i)
13      continue
        Do 30 I=1,NBasis*NBasis
30      URe(I)=UReNew(I)
        Do 31 I=1,NBasis
        Occ(I)=OccNew(I)
31      Continue
C
        if(gnorm.lt.GTol) then 
        Write(6,'(/,X,''Convergence in BFGS APSG attained after '',I3,
     $ '' iterations'')')its
        Deallocate (hessin)
        IConv=1
        return
        endif
C
C       Calculate a new inversion of the hessian and a new step xi
C
        do 17 i=1,n
          dg(i)=g(i)-dg(i)
17      continue

        do 35 i=1,n
35      hdg(i)=Zero

        ij=0
        do 40 i=1,n
        do 40 j=1,i
        ij=ij+1
        hij=hessin(ij)
        hdg(i)=hdg(i)+hij*dg(j)
        if(i.ne.j) hdg(j)=hdg(j)+hij*dg(i)
40      continue

        fac=Zero
        fae=Zero
        sumdg=Zero
        sumxi=Zero
        do 21 i=1,n
          fac=fac+dg(i)*xi(i)
          fae=fae+dg(i)*hdg(i)
          sumdg=sumdg+dg(i)**2
          sumxi=sumxi+xi(i)**2
21      continue
C
C       Skip update of hessin if fac not sufficiently positive
C
        if(fac**2.gt.EPS*sumdg*sumxi)then

          fac=One/fac
          fad=One/fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue

          do 43 i=1,n
          xin(i)=xi(i)
43        xi(i)=Zero

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
63        xi(i)=Zero
          ij=0
          do 65 i=1,n
            do 65 j=1,i
              ij=ij+1
              hij=hessin(ij)
              xi(i)=xi(i)-hij*g(j)
              if(i.ne.j) xi(j)=xi(j)-hij*g(i)
65        continue

        endif
C
C     End of the main loop over iterations
C
27    continue

      Write(6,'(/,X,''Too many iterations in BFGS APSG'')')
C
      Deallocate (hessin)
C
      Return
      End

*Deck LinAPSG
      Subroutine LinAPSG(n,xold,UReOld,CoefOld,fold,g,p,x,URe,Occ,f,
     $ Grad,check,
     $ XKin,XNuc,TwoEl,NBasis,NInte1,NInte2,MaxXV,NGem,
     $ IndX,IndN,NDimX,MxHVec)
C
C     LINEAR SERCH ALONG p DIRECTION
C     n = a total naumber of parameters to be optimized
C     xold STORES INITIAL X and x parameters
C     UReOld, CoefOld - the URe and c coeff's correspoonding to xold
C     fold - the energy value at xold
C     g - gradient (needed to find slope, g is not modified here)
C     x - final values of X and x
C     URe - final URe matrix
C     Occ - final Occ (final CICoef is in a common block)
C     f - final energy value
C     Grad - new grad 
C     check is false on a normal exit. it is true when x is too close to xold.
C     Parameters: ALF ensures sufficient decrease in function value; 
C     TOLX is the convergence criterion on \Delta x.
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension UReOld(NBasis*NBasis),URe(NBasis*NBasis),
     $ CoefOld(NBasis),Occ(NBasis),
     $ Grad(MxHVec),XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $ g(n),p(n),x(n),xold(n),IndX(NDimX),IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension step(n),StepAux(MaxXV)
C
      Logical check
C
      Parameter (ALF=1.d-4,TOLX=1.d-7, MXIT=5, StepMx=1.57D0/5.)
C
      check=.false.
C
C     SCALE THE INITIAL STEP (BEING p) IF NECESSARY
C     [StepMX = (Pi/2)/5 which should be ok for X and x)
C     taking too large step causes problems in Exp(X) ]
C
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
      If(IPrint.Ge.0)Write(6,'('' Initial step in LinAPSG scaled by '',
     $ F8.6,'' in linear search'')') ScMaxR 
      EndIf
C
C     CHECK IF THE INITIAL STEP IS STILL NOT TOO LONG
C
      Do I=1,MaxXV
      StepAux(I)=Zero
      EndDo
      Do I=1,NDimX
      StepAux(IndX(I))=p(I)
      EndDo
      Call NewUAPSG(URe,UReOld,StepAux,NBasis,ICheckU)
C
      If(ICheckU.Eq.0) Then
      Write(6,'
     $ ('' Initial step in LinAPSG scaled additionally by 0.1D0'')')
      Do I=1,n
      p(I)=0.1D0*p(I)
      EndDo
      EndIf  
C
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
C
C      GET NEW URe
C
       If(NDimX.Gt.0) Then
C
       Do I=1,MaxXV
       StepAux(I)=Zero
       EndDo
       Do I=1,NDimX
       StepAux(IndX(I))=Step(I)
       EndDo
       Call NewUAPSG(URe,UReOld,StepAux,NBasis,ICheckU)
C
       EndIf
C
C      GET NEW CICoef, Occ 
C
       If(n.Gt.NDimX) Call NewC(CoefOld,Step(NDimX+1),NBasis,NGem)
       Do I=1,NBasis
       Occ(I)=CICoef(I)**2
       EndDo
C
C      Calculate the energy and the gradient (to be used in BFGS)
C
       Call GrAPSG(f,Grad,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $      NInte1,NInte2,MaxXV,MxHVec,NGem,1)
C
        if(alam.lt.alamin)then

          do 16 i=1,n
            x(i)=xold(i)
16        continue
          Do 20 I=1,NBasis*NBasis
            URe(I)=UReOld(I)
20        Continue
          Do 21 I=1,NBasis
           CICoef(I)=CoefOld(I)
           Occ(I)=CICoef(I)**2
21        Continue
          check=.true.
          return

        else if(f.le.fold+ALF*alam*slope)then
          If(IPrint.Ge.0) Write(6,'(''Linear Search:'',
     $    I2,'' iter alam='',F7.4 )') it,alam
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
     $ Write(6,'('' No of iterations exceeded in LinAPSG! '')')

      Return
      End

*Deck ConjDir
      Subroutine ConjDir(h,xi,g,n)
C     
C     RETURNS A CONJUGATED GRADIENT DIR IN xi 
C
C     INPUTS:
C     xi - current gradient set to a new direction
C     g  - old gradient
C     h  - old direction (set to new)
C     
      Implicit Real*8 (A-H,O-Z)
C     
      Dimension h(n),xi(n),g(n)
C
      Parameter (Zero=0.D0)
C
      gg=zero
      dgg=zero
C     
      do j=1,n
      gg=gg+g(j)**2
C Fletcher-Reeves
c      dgg=dgg+xi(j)**2
C Polak-Ribiere
      dgg=dgg+(xi(j)-g(j))*xi(j)
      enddo
C    
      if(gg.eq.0) return 
      gam=dgg/gg
      do j=1,n
      g(j)=-xi(j)
      h(j)=g(j)+gam*h(j)
      xi(j)=h(j)
      enddo
C
      Return
      End

*Deck SortAll
      Subroutine SortAll(Occ,URe,NSymNO,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis,NBasis),NSymNO(NBasis)
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
C     CHANGE THE ORDER OF GEMINAL INDICES
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
C     CHANGE SYMMETRY INDICES
C
      Do I=1,NBasis
      IndOcc(I)=NSymNO(I)
      EndDo
      Do I=1,NBasis
      NSymNO(I)=IndOcc(Ind(I))
      EndDo
C
      Return
      End

*Deck ReWrAPSG
      Subroutine ReWrAPSG(IFlag,Occ,URe,UMOAO,Title,NGem,NBasis)
C
C     READ (IFlag=0) OR WRITE (IFlag=1) FROM/TO A RESTART FILE
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName,FMultTab
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis*NBasis),UMOAO(NBasis*NBasis)
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
      FName(K:K+4)='.res'
C
      If(IFlag.Eq.0) Then
C
      Open(10,Form='Unformatted',File=FName,Status='Old')
      Read(10) (CICoef(I),I=1,NBasis)
      Read(10) (Occ(I),I=1,NBasis)
      Read(10) (URe(I),I=1,NBasis*NBasis)
      Read(10) (UMOAO(I),I=1,NBasis*NBasis)
      Read(10) (IGem(I),I=1,NBasis)
      Read(10) NGem
      Close(10)
C
      Else
C
      Open(10,Form='Unformatted',File=FName)
      Write(10) (CICoef(I),I=1,NBasis)
      Write(10) (Occ(I),I=1,NBasis)
      Write(10) (URe(I),I=1,NBasis*NBasis)
      Write(10) (UMOAO(I),I=1,NBasis*NBasis)
      Write(10) (IGem(I),I=1,NBasis)
      Write(10) NGem
      Close(10)
C
      EndIf
C
      Return
      End

*Deck NewUAPSG
      Subroutine NewUAPSG(U,UOld,Step,N,ICheckU)
C
C     CALCULATE A NEW U MATRIX
C     U = Exp(Step) * UOld
C
C     If ICheckU=1 ON EXIT - STEP TOO LONG AND EXPM NOT CONVERGED
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension U(N*N),UOld(N*N),Step(N*(N-1)/2)
C
C     LOCAL ARRAY
C
      Dimension UHlp(N*N)
C
      Call EXPMAPSG(UHlp,Step,N,ICheckU)
C
      Do 5 I=1,N
      II=N*(I-1)+I
    5 UHlp(II)=UHlp(II)+1.D0
C
      Call MultpM(U,UHlp,UOld,N)
C
      Return
      End

*Deck EXPMAPSG
      Subroutine EXPMAPSG(T,X,N,ICheckU)
C
C     CALCULATE Exp(X)-1
C     T = X + (1/2) X**2 + (1/3!) X**3 + ...
C
C     ON EXIT:
C     ICheckU = 1  - CONVERGENCE ACHIEVED WITHIN MxIt 
C             = 0  - EXPANSION NOT CONVERGED
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,One=1.D0,Tol=1.D-12, MxIt=100)
C
      Dimension T(N*N),X(N*(N-1)/2)
C
C     LOCAL ARRAYS
C
      Dimension XN(N*N),XNN(N*N),X1(N*N)
C
      ICheckU=1
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
      If(Err.Lt.Tol) Return
    1 Continue
C
      Write(6,'("ERROR IN EXPM TOO HIGH! ",E12.2)')Err
      ICheckU=0
C
      Return
      End

*Deck GVBIni
      Subroutine GVBIni(Occ,NBasis,NGem,NGOcc) 
C
C     INITIATE GVB CACLUATIONS
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,One=1.D0)
C
      Dimension Occ(NBasis)
C
      Logical ICont
C
C     LOCAL ARRAYS
      Dimension IGemAux(NBasis)
C
      NGem=NELE+1
C
      Do I=1,NELE
      IGemAux(I)=IGem(I)
      EndDo
      Do I=NELE+1,NBasis
      IGemAux(I)=NGem
      EndDo
C
      Do I=1,NELE
C
      ICont=.True.
      Do J=NELE+1,NBasis
      If(ICont.And.IGem(J).Eq.I) Then
      IGemAux(J)=IGem(J)
      ICont=.False.
      EndIf
      EndDo
C     
      EndDo 
C
      Do I=1,NBasis
      IGem(I)=IGemAux(I)
      EndDo
C
C     PUT TO ZERO Occ AND CICoef FOR THE LAST GEMINAL
C
      Do I=1,NBasis
C
      If(IGem(I).Eq.NGem) Then
C
      Occ(I)=Zero
      CICoef(I)=Zero
C
      Else
C
C     MODIFY THE Occ 
C
      If(I.Gt.NELE) Then
C
      Do J=1,NELE
      If(IGem(J).Eq.IGem(I)) Then
      Occ(I)=One-Occ(J)
      If(CICoef(J).Gt.Zero) Then
      CICoef(I)=-SQRT(Occ(I))
      Else
      CICoef(I)=SQRT(Occ(I))
      EndIf
      EndIf
      EndDo
C
      EndIf
C
      EndIf
C
      EndDo
C
C     SET OCCUPANCIES AND CICoef TO 1 IF NGOcc.Ne.0
C
c      ngocc=8

      Do I=1,NGOcc
C
      CICoef(I)=One
      Occ(I)=One
      IG=IGem(I)
c
      Do J=I+1,NBasis
      If(IGem(J).Eq.IG) Then
      Occ(J)=Zero
      CICoef(J)=Zero    
      IGem(J)=NGem
      EndIf
      EndDo
C
      EndDo
C
      Write(6,'(/,X,'' INITIAL ORBITAL OCCUPANCIES IN IGVB'')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I)
      EndDo
C
      Return
      End
