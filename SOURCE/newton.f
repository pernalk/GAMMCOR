*Deck NewtonAPSG
      Subroutine NewtonAPSG(ETot,URe,Occ,XKin,XNuc,TwoEl,
     $ IndX,IndN,NDimX,IndC,NDimN,NGOcc,
     $ NN,NBasis,NInte1,NInte2,MaxXV,MxHVec,NGem,GTol,MxIt)
C
C     NEWTON-RAPHSON FOR APSG  
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Three=3.D0,
     $ Four=4.D0,Five=5.D0,Eight=8.0D0)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $ IndX(NDimX),IndN(2,NDimX),IndC(NDimN)
C
C     LOCAL ARRAYS
C
      Dimension
     $ Grad(MxHVec),Step(MxHVec),UReOld(NBasis,NBasis),OccOld(NBasis),
     $ CICoefOld(NBasis),HessLow(MxHVec*(MxHVec+1)/2),
     $ XProj(MxHVec*(MxHVec+1)/2)
      Dimension
     $ Hess(NN,NN),G(NN),St(NN),EigVal(NN),Work(NN)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     Maximal step for X (should be less than Pi/2, otherwise Exp(X) does not converge)
C
      StepMx=ASin(One)/Five
      NDrop=Zero
C
      Iter=0
C
C     BEGINNING OF THE ITERATIONS
C
      EOld=Zero
C
   10 Iter=Iter+1
C      
C     CALCULATE GRADIENT AND THE ENERGY 
C
      Call HessAPSG(ETot,Grad,HessLow,URe,Occ,XKin,XNuc,TwoEl,NBasis,
     $NInte1,NInte2,MaxXV,MxHVec,NGem,2)
C
C     FIND THE PROJECTION MATRIX AND PROJECT THE HESSIAN
C
      Call Projection(XProj,NBasis,MaxXV,MxHVec,NGem)
      Call ProjectHess(XProj,HessLow,NBasis,MaxXV,MxHVec)
C
      Do I=1,NN
      If(I.Le.NDimX) G(I)=Grad(IndX(I))
      If(I.Gt.NDimX) G(I)=Grad(MaxXV+IndC(I-NDimX))
      EndDo
C
C     CALCULATE THE GRADIENT NORM
C
      GradNo=Zero
C
      Do 20 I=1,NN
   20 GradNo=GradNo+G(I)*G(I)
      GradNo=Sqrt(GradNo)
C
C     COMMENCE LINEAR SEARCH IF ENERGY RISES 
C     OR HESSIAN IS NOT POSITIVE-DEFINITE
C
c      If((ETot.Gt.EOld.Or.INeHEV.Ne.0).
c     $    And.Iter.Ne.1.And.GradNo.Gt.Sqrt(GTol)) Then
cC
cC     end of linear search 
c      EndIf
C
      Do 80 I=1,NBasis
      OccOld(I)=Occ(I)
      CICoefOld(I)=CICoef(I)
      Do 80 J=1,NBasis
   80 UReOld(I,J)=URe(I,J)  
C
      If(IPrint.Ge.0)
     $ Write(6,'(/,X,''QSCF ITER'',I3,2X,''ENERGY'',F16.8,2X,
     $ ''ENE DIFF '',E10.3,2X,''GRAD NORM '',E10.3)')
     $ Iter,ETot,ETot-EOld,GradNo
      If(IPrint.Ge.1) Then
      Write(6,'(/,X,''  ORBITAL OCCUPANCIES '')')
      Write(6,'(2X,"Orb",3X,"Occupancy",7X,"Gem")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,I6)') I,Occ(I),IGem(I)
      EndDo
      EndIf
C
      EOld=ETot
C
C     COPY HessLow TO Hess
C
      Do I=1,NN
      Do J=1,I
C
      If(I.Le.NDimX) Then
C
      II=IndX(I)
C
      If(J.Le.NDimX) Then
      JJ=IndX(J)
      Else
      JJ=MaxXV+IndC(J-NDimX) 
      EndIf
C
C     else to I.Le.NDimx
      Else
C
      II=MaxXV+IndC(I-NDimX)
C
      If(J.Le.NDimX) Then
      JJ=IndX(J)
      Else
      JJ=MaxXV+IndC(J-NDimX)
      EndIf
C
c     endif of I.Le.NDimX
C
      EndIf
C
      IJ=(Max(II,JJ)*(Max(II,JJ)-1))/2+Min(II,JJ)
      Hess(I,J)=HessLow(IJ)
      Hess(J,I)=Hess(I,J)
C
      EndDo
      EndDo
C
C     COMPUTE THE NORMS OF Hxx AND Hvv
C
      HXXNo=Zero
C
      Do 110 I=1,NDimX
      Do 110 J=I,NDimX
      Fac=Two
      If (I.Eq.J) Fac=One
  110 HXXNo=HXXNo+Fac*Hess(I,J)**2
C
      HXXNo=HXXNo**(One/Four)
c      HXXNo=One
C
      HVVNo=Zero
C
      Do 120 I=NDimX+1,NN
      Do 120 J=I,NN
      Fac=Two
      If (I.Eq.J) Fac=One
  120 HVVNo=HVVNo+Fac*Hess(I,J)**2
C
      HVVNo=HVVNo**(One/Four)
c      HVVNo=One
C
C     BLOCK-SCALE THE HESSIAN
C
      Do 130 I=1,NDimX
      Do 130 J=1,NDimX
  130 Hess(I,J)=Hess(I,J)/HXXNo**2
C
      Do 140 I=1,NDimX
      Do 140 J=NDimX+1,NN
      Hess(I,J)=Hess(I,J)/HXXNo/HVVNo
  140 Hess(J,I)=Hess(J,I)/HXXNo/HVVNo
C
      Do 150 I=NDimX+1,NN
      Do 150 J=NDimX+1,NN
  150 Hess(I,J)=Hess(I,J)/HVVNo**2       
C
C     DIAGONALIZE THE HESSIAN 
C
      Call Diag8(Hess,NN,NN,EigVal,Work)
C    
C     COUNT THE NEGATIVE EIGENVALUES OF HESSIAN
C
      INeHEV=0
C
      write(*,*)'negative and smaller then 1.d-5 eigenvalues'
      Do 160 I=1,NN
      If(EigVal(I).Lt.Zero) INeHEV=INeHEV+1
c herer!!!
c      if(ndimx.ne.zero.and.nn.gt.ndimx) then
      if(eigval(i).lt.zero.or.abs(eigval(i)).lt.1.d-5) then
      write(*,*)'eig',i,EigVal(I) 
      endif
c      endif

  160 Continue
C
C     PROJECT OUT ZERO MODES 
C
c herer!!!
c      goto 777
      If(NN.Gt.NDimX) Then
C
      NGemEnd=NGem
      If(NGem.Gt.NELE) NGemEnd=NELE
      IGIni=1
      If(IFreeze.Eq.1) IGIni=NGOcc+1
C
      Do IG=IGIni,NGemEnd
C     
      SMax=Zero 
C
      Do I=1,NN
C
      SS=Zero
      Do J=1,NDimN
      If(IGem(IndC(J)).Eq.IG) 
     $ SS=SS+Abs(CICoef(IndC(J))*Hess(I,NDimX+J))
      EndDo
C
      If(SS.Gt.SMax) Then
      SMax=SS
      IndMax=I
      EndIf
C
      EndDo
C
      If(Eigval(IndMax).Lt.Zero) INeHEV=INeHEV-1
      write(*,*)'projected mode, overlap',indmax,Eigval(IndMax),SMax 
      Eigval(IndMax)=Zero
C
C     enddo of IG
      EndDo
C
C     endif of NN.Gt.NDimX
      EndIf
C
  777 write(*,*)
C
      Do 245 I=1,NN
c
c herer
      If (Abs(Eigval(I)).Lt.1.D-7) Then
c herer!!!
      if(eigval(i).ne.zero)write(*,*)'eigval',i,eigval(i),'put to zero'
      Eigval(I)=Zero
      Else
      Eigval(I)=Abs(One/Eigval(I))
      EndIf
  245 Continue      
C
      If(INeHEV.Ne.0) 
     $ Write(6,'('' The Hessian has '',I4,'' negative eigenvalues'')') 
     $ INeHEV
C
C     BLOCK-SCALE GRADIENT
C
      Do 250 I=1,NN
      If (I.Le.NDimX) Then
      G(I)=G(I)/HXXNo
      Else
      G(I)=G(I)/HVVNo
      EndIf
  250 Continue  
C
C     COMPUTE STEP
C
      Do 260 I=1,NN
  260 St(I)=Zero
C      
      Do 280 K=1,NN
C
      Sum=Zero
      Do 270 J=1,NN
  270 Sum=Sum+Hess(K,J)*G(J)
C
      Sum=Sum*EigVal(K)
C
      Do 280 I=1,NN
  280 St(I)=St(I)-Sum*Hess(K,I)             
C
C     BLOCK-SCALE STEP
C
      Do 290 I=1,NN
      If (I.Le.NDimX) Then
      St(I)=St(I)/HXXNo
      Else
      St(I)=St(I)/HVVNo
      EndIf
  290 Continue
C
C     PROJECT STEP
C
      If(NN.Gt.NDimX) Then
C
      NGemEnd=NGem
      If(NGem.Gt.NELE) NGemEnd=NELE
      IGIni=1
      If(IFreeze.Eq.1) IGIni=NGOcc+1
C
      Do IG=IGIni,NGemEnd
C
C     FIND A PROJECTION OF THE STEP ON THE IG-TH REDUNDANT VECTOR
C
      SS=Zero
      Do J=1,NDimN
      If(IGem(IndC(J)).Eq.IG)
     $ SS=SS+CICoef(IndC(J))*St(NDimX+J)
      EndDo
C
      Do J=1,NDimN
      If(IGem(IndC(J)).Eq.IG)
     $ St(NDimX+J)=St(NDimX+J)-SS*CICoef(IndC(J))
      EndDo
C
C     enddo of IG
      EndDo
C
C     endif of NN.Gt.NDimX
      EndIf
C
C     ASSURE THAT THE ORBITAL ROTATION ANGLE DOES NOT EXCEED Pi/2 
C
      ScMaxR=Zero
      Do 340 I=1,NDimX
  340 ScMaxR=Max(ScMaxR,Abs(St(I)))
C
      If (ScMaxR.Ne.Zero) Then
      ScMaxR=Min(StepMx/ScMaxR,One)
      Else
      ScMaxR=One
      EndIf

c herer!!!
      if(ndimx.ne.zero.and.nn.gt.ndimx) then
      do i=1,NDimX
      II=IndN(1,I)
      JJ=IndN(2,I)
c      write(*,*)'step',ii,jj,St(i),g(i)
      enddo
      do j=1,nbasis
c      write(*,*)'step',j,st(ndimx+j),g(ndimx+j)
      enddo
      endif
C
C     SCALE THE STEP IF NECESSARY
C
      If(ScMaxR.Ne.One) Then
      Do 350 I=1,NN
  350 St(I)=ScMaxR*St(I)
      If(IPrint.Ge.0) 
     $ Write(6,'('' Step scaled by '',F8.6,
     $ '' due to orbital rotation constraints'')') ScMaxR
      EndIf
C
C     COPY St TO Step
C
      Do I=1,MxHVec
      Step(I)=Zero
      EndDo
C
      Do I=1,NDimX
      Step(IndX(I))=St(I)
      EndDo
      If(NN.Gt.NDimX) Then
      Do I=1,NDimN
      Step(MaxXV+IndC(I))=St(NDimX+I)
      EndDo
      EndIf
C
C     CHECK IF THE INITIAL STEP IS STILL NOT TOO LONG
C
      Call NewUAPSG(URe,UReOld,Step,NBasis,ICheckU)
      If(ICheckU.Eq.0) Then
      Write(6,'
     $ ('' Step scaled additionally by 0.1D0'')')
      Do I=1,MxHVec
      Step(I)=0.1D0*Step(I)
      EndDo
      EndIf
C
C     FIND NEW URe, Occ, CICoef    
C
      Call NewUAPSG(URe,UReOld,Step,NBasis,ICheckU)
      Call NewC(CICoefOld,Step(MaxXV+1),NBasis,NGem)
      Do I=1,NBasis
      Occ(I)=CICoef(I)**2
      EndDo
C
C     CHECK THE GRADIENT NORM FOR CONVERGENCE
C
      If((GradNo.Gt.GTol).And.(Iter.Lt.MxIt)) GoTo 10 
C
      If(Iter.Eq.MxIt) Write(*,*)'FATAL ERROR: NO CONVERGENCE IN QSCF!'
C
      Return
      End

*Deck Projection
      Subroutine Projection(XProj,NBasis,MaxXV,MxHVec,NGem)
C
C     A PROJECTION MATRIX IS CONSTRUCTED HERE P=1-O^T*O
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Dimension XProj(MxHVec*(MxHVec+1)/2)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
      IJ=0
      Do I=1,MxHVec
      Do J=1,I
      IJ=IJ+1
      XProj(IJ)=Zero
      If(I.Eq.J) XProj(IJ)=One 
C
      If(I.Gt.MaxXV.And.J.Gt.MaxXV) Then
C
      IP=I-MaxXV    
      IQ=J-MaxXV
      SumQ=Zero
      Do IG=1,NGem
      If(IGem(IP).Eq.IG.And.IGem(IQ).Eq.IG) 
     $ SumQ=SumQ+CICoef(IP)*CICoef(IQ) 
      EndDo
C
      XProj(IJ)=XProj(IJ)-SumQ
C
      EndIf
C
      EndDo
      EndDo
C
      Return
      End      

*Deck ProjectHess
      Subroutine ProjectHess(XProj,HessLow,NBasis,MaxXV,MxHVec)
C     
C     A PROJECTION MATRIX IS CONSTRUCTED HERE P=1-O^T*O
C     
      Implicit Real*8 (A-H,O-Z)
C     
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C     
      Dimension HessLow(MxHVec*(MxHVec+1)/2),XProj(MxHVec*(MxHVec+1)/2)
C
C     LOCAL ARRAY
C
      Dimension Aux1(NBasis,NBasis),Aux2(NBasis,NBasis),
     $ Aux3(NBasis,NBasis)
C
      Character*60 FMultTab
      Include 'commons.inc'
C
C     COPY xx BLOCKS OF HESSIAN AND XProj TO AUXILIARY SQUARE ARRAYS
C
      Do I=1,NBasis
      Do J=1,NBasis
C
      II=MaxXV+I
      JJ=MaxXV+J
      IJ=(Max(II,JJ)*(Max(II,JJ)-1))/2+Min(II,JJ)
C
      Aux1(I,J)=HessLow(IJ)
      Aux1(J,I)=Aux1(I,J)
      Aux2(I,J)=XProj(IJ)
      Aux2(J,I)=Aux2(I,J)
C
      EndDo
      EndDo
C
      Call MultpM(Aux3,Aux1,Aux2,NBasis)
      Call MultpM(Aux1,Aux2,Aux3,NBasis) 
C
C     COPY THE PROJECTED BLOCK BACK TO HESSIAN
C
      Do I=1,NBasis
      Do J=1,I
C
      II=MaxXV+I
      II=MaxXV+J
      IJ=(Max(II,JJ)*(Max(II,JJ)-1))/2+Min(II,JJ)
C
      HessLow(IJ)=Aux1(J,I)
C
      EndDo
      EndDo
C
      Return
      End  


 
