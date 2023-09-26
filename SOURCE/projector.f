*Deck OptGam2
      Subroutine OptGam2(ETot,Occ,URe,XKin,XNuc,TwoEl,TwoElF,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ NInte1,NInte2,NBasis)
C
C     CARRIES OUT ENERGY OPTIMIZATION USING THE PROJECTION GRADIENT METHOD 2
C     Cances&Pernal JCP 128, 134108 (2008), Eq.(19)
C
C     WORKS ONLY FOR THE BB FUNCTIONAL OR BB+HF HYBRID!
C
C     TwoEl  - integrals with erf interaction
C     TwoElF - integrals with full Coulomb interaction 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis,NBasis),
     $          XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),TwoElF(NInte2),
     $          XOne(NInte1)
C
      Dimension WGrid(*),OrbXGrid(*),OrbYGrid(*),OrbZGrid(*),OrbGrid(*),
     $          NSymMO(NBasis),VSR(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension GamSR(NInte1),GR(NInte1),Gam(NInte1),Hcc(NBasis),
     $ GROLD(NInte1),Dir(NInte1)
C 
      Parameter(TolN=1.D-2,TolO=1.D-4,ETol=1.D-5,MxIt=400)
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
C     CHECK THE FUNCTIONAL
C
      If(IFun.Ne.2.And.IFun.Ne.8.And.IFun.Ne.10) 
     $ Stop'Fatal Error in OptGam2: the functional in use is not BB!'
C
C     CONSTRUCT THE INITIAL MATRIX Gamma^(1/2)
C
      Do I=1,NBasis
      Hcc(I)=SQRT(Occ(I))
      EndDo
      Call VecTr(GamSR,Hcc,URe,NBasis)
C
C     BEGINNING OF THE ITERATIONS
C
      IStall=0
      ItS=0
C
      Do It=1,MxIt
C
      ItS=ItS+1
C
C     COMPUTE THE GRADIENT W.R.T. Gamma^(1/2)
C
      ETotO=ETot 
C
      Do I=1,NInte1
      XOne(I)=XKin(I)
      EndDo
C
      If(IFun.Eq.10) Then
      Call EnGrBBsrHF(1,ETot,GR,GamSR,XKin,XNuc,TwoEl,TwoElF,NInte1,
     $ NInte2,NBasis)
      Else
      Call EnGrBB(1,ETot,GR,GamSR,XOne,XNuc,TwoEl,NInte1,NInte2,NBasis)
      EndIf
C
C     COMPUTE THE ERROR AND CHECK THE CONVERGENCE
C
      Call GetError(XMiu,ErrN,Com,GamSR,Occ,URe,GR,NInte1,NBasis)
      Write(6,'(X,''ITER'',I3,2X,''ENERGY'',F16.8,5X,''ENE DIFF'',
     $ 2X,E10.3,2X,''ERRS''2E12.2)') It,ETot,ETot-ETotO,Com,ErrN
C
      If(Com.Le.TolO.And.ErrN.Lt.TolN.And.ETotO-ETot.Le.ETol) Then
C
      Write(6,'(/,X,''CONVERGENCE ATTAINED IN THE PROJECTED GRADIENT 
     $ ALGORITHM'',/)')
C
      Call SortOcc(Occ,URe,NBasis)
      Write(6,'(2X,'' ORBITAL OCCUPANCIES '')')
      Sum=Zero
      Do 7 I=1,NBasis
      Sum=Sum+Occ(I)
    7 Write(6,'(1X,I3,E16.6)') I,Occ(I)
C
      Write(6,'(''Norm='',F13.8,/)')Sum
      Return
C
      EndIf
C
C     OBTAIN A CONJUGATE GRADIENT
C
C note that the gradient is not ok anyway so reset the direction more frequently
C than every nbasis cycles
      If(It.Eq.1.Or.IStall.Eq.1.Or.ItS.Eq.NBasis/2) Then 
      ItS=0
      IStall=0
      Do I=1,NInte1
      Dir(I)=GR(I)
      GROLD(I)=GR(I)
      EndDo
      Else
      Call Conjugate(Dir,GR,GROLD,NInte1)
      EndIf
C
      Call LinMinP(IStall,ETot,GR,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
C
C     PROJECT THE SQUARE ROOT OF GAMMA
C
      Call Project2(GamSR,URe,Occ,NInte1,NBasis)
C
      EndDo 
C
      Return
      End

*Deck EnGrBB 
      Subroutine EnGrBB(IFlag,ETot,GR,GamSR,XKin,XNuc,TwoEl,
     $ NInte1,NInte2,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
C     COMPUTES THE ENERGY AND/OR THE GRADIENT WITH RESPECT TO Gamma^1/2 FOR A BB FUNCTIONAL
C
      Dimension GamSR(NInte1),GR(NInte1),
     $          XKin(NInte1),XNuc(NInte1),TwoEl(NInte2)
C
C     LOCAL ARRAYS (FHF STORES THE HF-EXCHANGE PART OF GRADIENT)
C
      Dimension Gam(NInte1),F(NInte1),Fxc(NInte1),FHF(NInte1)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
C     SET A MIXING OF HF-EXCHANGE (AND REMOVE THE SAME AMOUNT OF BB)
C
      Xmix=Zero
      If(IFun.Eq.8) Xmix=Cmix
C
C     OBTAIN GAMMA
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Gam(IJ)=Zero
C
      Do K=1,NBasis
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      Gam(IJ)=Gam(IJ)+GamSR(IK)*GamSR(JK)
      EndDo
C
      EndDo
      EndDo
C
C     ONE-ELECTRON CONTRIBUTION TO F
C
      IJ=0
      Do 20 I=1,NBasis
      Do 20 J=1,I
      IJ=IJ+1
      Fxc(IJ)=Zero
      FHF(IJ)=Zero
   20 F(IJ)=XKin(IJ)+XNuc(IJ)
C  
C     TWO-ELECTRON CONTRIBUTION TO F AND Fxc 
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
      F(IJ)=F(IJ)+FacIJ*Gam(KL)*TwoZet
      If(IJ.Ne.KL) F(KL)=F(KL)+FacKL*Gam(IJ)*TwoZet
C
      Fxc(IL)=Fxc(IL)-Half*FacIL*GamSR(JK)*TwoZet
      FHF(IL)=FHF(IL)-Half*FacIL*Gam(JK)*TwoZet
C
      If(IL.Ne.JK) Then
C
      Fxc(JK)=Fxc(JK)-Half*FacJK*GamSR(IL)*TwoZet
      FHF(JK)=FHF(JK)-Half*FacJK*Gam(IL)*TwoZet
C
      EndIf
C
      If(K.Eq.L) GoTo 40
C
      If(I.Eq.J) GoTo 40
C
      F(IJ)=F(IJ)+Gam(KL)*TwoZet
      If(IJ.Ne.KL) F(KL)=F(KL)+Gam(IJ)*TwoZet
C
      Fxc(IK)=Fxc(IK)-Half*FacIK*GamSR(JL)*TwoZet
      FHF(IK)=FHF(IK)-Half*FacIK*Gam(JL)*TwoZet
C
      If(IK.Ne.JL) Then
C
      Fxc(JL)=Fxc(JL)-Half*FacJL*GamSR(IK)*TwoZet
      FHF(JL)=FHF(JL)-Half*FacJL*Gam(IK)*TwoZet
C
      EndIf
C
   40 Continue
   30 Continue
C
C     COMPUTE THE TOTAL ENERGY
C
      ETot=Zero
      IJ=0
C
      Do 45 I=1,NBasis
      Do 45 J=1,I
      IJ=IJ+1
C
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
C
   45 ETot=ETot
     $ +FacIJ*Gam(IJ)*(XKin(IJ)+XNuc(IJ)+F(IJ)+Xmix*FHF(IJ))
     $ +FacIJ*(One-Xmix)*GamSR(IJ)*Fxc(IJ)
C
      If(IFlag.Eq.0) Return
C
C     COMPLETE THE GRADIENT       
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I 
      IJ=IJ+1
      Fxc(IJ)=(One-Xmix)*Fxc(IJ)
C
      Do K=1,NBasis
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      Fxc(IJ)=Fxc(IJ)+F(IK)*GamSR(JK)+F(JK)*GamSR(IK)
     $ +Xmix*(FHF(IK)*GamSR(JK)+FHF(JK)*GamSR(IK))
      EndDo
C
      EndDo
      EndDo
C
      Do I=1,NInte1
      GR(I)=Two*Fxc(I)
      EndDo 
C
      Return
      End

*Deck EnGrBBsrHF
      Subroutine EnGrBBsrHF(IFlag,ETot,GR,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NInte1,NInte2,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
C     COMPUTES THE ENERGY AND/OR THE GRADIENT WITH RESPECT TO Gamma^1/2 FOR A BB FUNCTIONAL
C
      Dimension GamSR(NInte1),GR(NInte1),
     $          XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),TwoElF(NInte2)
C
C     LOCAL ARRAYS (FHF STORES THE SR-HF-EXCHANGE PART OF GRADIENT)
C
      Dimension Gam(NInte1),F(NInte1),Fxc(NInte1),FHF(NInte1)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0)
C
C     LR-BB+SR-HF 
C
C     OBTAIN GAMMA
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Gam(IJ)=Zero
C
      Do K=1,NBasis
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      Gam(IJ)=Gam(IJ)+GamSR(IK)*GamSR(JK)
      EndDo
C
      EndDo
      EndDo
C
C     ONE-ELECTRON CONTRIBUTION TO F
C
      IJ=0
      Do 20 I=1,NBasis
      Do 20 J=1,I
      IJ=IJ+1
      Fxc(IJ)=Zero
      FHF(IJ)=Zero
   20 F(IJ)=XKin(IJ)+XNuc(IJ)
C  
C     TWO-ELECTRON CONTRIBUTION TO F AND Fxc 
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
      TwoZetSR=Two*(TwoElF(NAddr)-TwoEl(NAddr))
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
      F(IJ)=F(IJ)+FacIJ*Gam(KL)*TwoZet
      If(IJ.Ne.KL) F(KL)=F(KL)+FacKL*Gam(IJ)*TwoZet
C
      Fxc(IL)=Fxc(IL)-Half*FacIL*GamSR(JK)*TwoZet
      FHF(IL)=FHF(IL)-Half*FacIL*Gam(JK)*TwoZetSR
C
      If(IL.Ne.JK) Then
C
      Fxc(JK)=Fxc(JK)-Half*FacJK*GamSR(IL)*TwoZet
      FHF(JK)=FHF(JK)-Half*FacJK*Gam(IL)*TwoZetSR
C
      EndIf
C
      If(K.Eq.L) GoTo 40
C
      If(I.Eq.J) GoTo 40
C
      F(IJ)=F(IJ)+Gam(KL)*TwoZet
      If(IJ.Ne.KL) F(KL)=F(KL)+Gam(IJ)*TwoZet
C
      Fxc(IK)=Fxc(IK)-Half*FacIK*GamSR(JL)*TwoZet
      FHF(IK)=FHF(IK)-Half*FacIK*Gam(JL)*TwoZetSR
C
      If(IK.Ne.JL) Then
C
      Fxc(JL)=Fxc(JL)-Half*FacJL*GamSR(IK)*TwoZet
      FHF(JL)=FHF(JL)-Half*FacJL*Gam(IK)*TwoZetSR
C
      EndIf
C
   40 Continue
   30 Continue
C
C     COMPUTE THE TOTAL ENERGY
C
      ETot=Zero
      IJ=0
C
      Do 45 I=1,NBasis
      Do 45 J=1,I
      IJ=IJ+1
C
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
C
   45 ETot=ETot
     $ +FacIJ*Gam(IJ)*(XKin(IJ)+XNuc(IJ)+F(IJ)+FHF(IJ))
     $ +FacIJ*GamSR(IJ)*Fxc(IJ)
C
      If(IFlag.Eq.0) Return
C
C     COMPLETE THE GRADIENT       
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I 
      IJ=IJ+1
C
      Do K=1,NBasis
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      Fxc(IJ)=Fxc(IJ)+F(IK)*GamSR(JK)+F(JK)*GamSR(IK)
     $ +(FHF(IK)*GamSR(JK)+FHF(JK)*GamSR(IK))
      EndDo
C
      EndDo
      EndDo
C
      Do I=1,NInte1
      GR(I)=Two*Fxc(I)
      EndDo 
C
      Return
      End

*Deck GetError
      Subroutine GetError(XMiuAv,ErrN,Com,GamSR,Occ,URe,GR,
     $ NInte1,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Common/DERIVA/ DEN(1000),Error
C
      Dimension GamSR(NInte1),GR(NInte1),Hlp(NInte1),
     $ Occ(NBasis),URe(NBasis,NBasis) 
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,
     $          Small=1.D-7)
C
C     CALCULATES THE ERROR FOR THE OCCUPANCIES AND THE ORBITALS (COMMUTATOR)
C
      Do I=1,NInte1
      Hlp(I)=GR(I)
      EndDo
      Call MatTr(Hlp,URe,NBasis)
      XMiuAv=Zero
      ICount=0
      Do I=1,NBasis
      II=I*(I-1)/2+I
      DEN(I)=Hlp(II)
      If(Abs(Occ(I)-One).Gt.Small.And.Occ(I).Gt.Small) Then
      ICount=ICount+1
      XMiuAv=XMiuAv+Hlp(II)/Sqrt(Occ(I))
      EndIf
      EndDo
      XMiuAv=XMiuAv/Float(ICount)
      ErrN=Zero
C
      Do I=1,NBasis
      II=I*(I-1)/2+I
      If(Abs(Occ(I)-One).Gt.Small.And.Occ(I).Gt.Small) 
     $ ErrN=ErrN+(XMiuAv-Hlp(II)/Sqrt(Occ(I)))**2
      EndDo
      ErrN=Sqrt(ErrN)/Float(ICount)
      Error=ErrN
C
C     COMMUTATOR
C
      Com=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      ComIJ=Zero
      Do K=1,NBasis
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      ComIJ=ComIJ+GR(IK)*GamSR(JK)-GamSR(IK)*GR(JK)
      EndDo
      Com=Com+ComIJ**2
      EndDo
      EndDo
      Com=SQRT(Com)
C
      Return
      End

*Deck Conjugate
      Subroutine Conjugate(h,xi,g,n)
C     
C     xi - current gradient and new direction
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
      dgg=dgg+xi(j)**2
C Polak-Ribiere
c      dgg=dgg+(xi(J)+g(j))*xi(j)
      enddo
C     
      gam=dgg/gg
      do j=1,n
      g(j)=xi(j)
      h(j)=g(j)+gam*h(j)
      xi(j)=h(j)
      enddo
C
      Return
      End

*Deck LinMinP
      Subroutine LinMinP(IStall,ETot0,Dir,GamSR,XKin,XNuc,
     $ TwoEl,TwoElF,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0)
C     
      Dimension Dir(NInte1),GamSR(NInte1)
      Dimension XKin(*),XNuc(*),TwoEl(*),TwoElF(*)
C
      Tol=0.5D0
C     
      Call BrackP(istall,ifail,etot0,ax,bx,cx,fa,fb,fc,
     $ Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,NInte2)
C
C     IF ifail=1 TAKE A LONG STEP AND LEAVE 
C     
      If(IFail.Eq.1) Then
      XLam=CX
      GoTo 999
      EndIf
C     
      Fret=BrentP(ax,bx,cx,Tol,XLam,
     $ Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,NInte2)
C
      If(XLam.Eq.Zero.Or.Fret.Gt.ETot0) Then
      IStall=1
      XLam=-1.D-5
      EndIf
C
  999 Continue
      Do I=1,NInte1
      GamSR(I)=GamSR(I)+XLam*Dir(I)
      EndDo
C
      Return
      End

*Deck EneDirP
      Real*8 Function EneDirP(XLam,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
C
C     COMPUTE THE TOTAL ELECTRONIC ENERGY AT A POINT GamSR+XLam*DIR 
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Dimension Dir(NInte1),GamSR(NInte1)
      Dimension XKin(*),XNuc(*),TwoEl(*),TwoElF(*)
C
C     LOCAL ARRAYS
C
      Dimension Hlp(NInte1),URe(NBasis,NBasis),Occ(NBasis),
     $ GR(NInte1)
c
      Do J=1,NInte1
      Hlp(J)=GamSR(J)+XLam*Dir(J)
      EndDo
C
      Call Project2(Hlp,URe,Occ,NInte1,NBasis)
C
      If(IFun.Eq.10) Then
      Call EnGrBBsrHF(0,ETot,GR,Hlp,XKin,XNuc,TwoEl,TwoElF,NInte1,
     $ NInte2,NBasis)
      Else
      Call EnGrBB(0,ETot,GR,Hlp,XKin,XNuc,TwoEl,NInte1,NInte2,NBasis)
      EndIf
C
      EneDirP=ETot
C
      Return
      End

*Deck BrentP
      Real*8 Function BrentP(ax,bx,cx,tol,xmin,
     $ Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension XKin(*),XNuc(*),TwoEl(*),TwoElF(*),Dir(*),GamSR(*)
C
      PARAMETER (ITMAX=100,CGOLD=.3819660D0,ZEPS=1.0D-10)
C
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
C
      fx=EneDirP(x,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,
     $ NInte2)
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
        fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,
     $  NInte2)
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
      stop 'brentocc exceed maximum iterations'
3     xmin=x
      brentp=fx
      Return
      END

*Deck BrackP
      Subroutine BrackP(istall,ifail,etot0,ax,bx,cx,fa,fb,fc,
     $ Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,NInte2)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Dir(*),GamSR(*),XKin(*),XNuc(*),TwoEl(*),TwoElF(*)
C
      PARAMETER (GOLD=1.618034D0, GLIMIT=100.D0, TINY=1.d-20)
      PARAMETER (ZERO=0.0D0, ONE=1.D0, TEN1=0.1D0, TEN5=5.D-2, MXIT=10)
C
      istall=0
      ifail=0
C
      ax=-TEN5 
      fa=EneDirP(ax,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,NBasis,NInte1,
     $ NInte2)
C
      if(fa.gt.etot0) then 
C
      ax=zero
      fa=EneDirP(ax,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
      bx=-0.0005D0
      fb=EneDirP(bx,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)   
      cx=-0.001D0
      fc=EneDirP(cx,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
C
      else
C
      bx=-TEN1
      fb=EneDirP(bx,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
      cx=-0.3D0
      fc=EneDirP(cx,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF, 
     $ NBasis,NInte1,NInte2)
C
      endif   
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
      its=0
    1 its=its+1
C
      if(its.gt.MXIT.or.u.le.-one) then
      istall=1
      ifail=1
      cx=u
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
          fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $   NBasis,NInte1,NInte2)
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
          fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
C
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
C
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
          endif
C
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
        else
          u=cx+GOLD*(cx-bx)
          fu=EneDirP(u,Dir,GamSR,XKin,XNuc,TwoEl,TwoElF,
     $ NBasis,NInte1,NInte2)
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
      Return
      END

*Deck Project2
      Subroutine Project2(GamSR,URe,Occ,NInte1,NBasis)
C     
      Implicit Real*8 (A-H,O-Z)
C     
C     PROJECT GamSR ON THE Q-REPRESENTABLE SPACE
C     ( Tr[GamSR^2] <= N ), Eq.(19) from JCP 128, 134108 (2008)
C
      Common/DERIVA/ DEN(1000),Error
C
      Dimension GamSR(NInte1),URe(NBasis,NBasis),Occ(NBasis)
C     
C     LOCAL ARRAYS
C     
      Dimension Pcc(NBasis),Gam(NInte1),Temp(NBasis)
C     
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0,Four=4.D0,
     $ Small=1.D-9)
C     
      Include 'commons.inc'
C
C     DIAGONALIZE GAMSR TO GET Pcc
C
      Call CpySym(URe,GamSR,NBasis)
      Call Diag8(URe,NBasis,NBasis,Pcc,Temp)
C
C     CHECK IF GamSR IS ALREADY N-REPRESENTABLE
C
      ICheck=0
      Sum=Zero
      Do I=1,NBasis
      Occ(I)=Pcc(I)**2
      Sum=Sum+Pcc(I)**2
      If(Pcc(I).Lt.Zero.Or.Pcc(I).Gt.One) ICheck=1
      EndDo
      If(Sum.Le.XELE.And.ICheck.Eq.0) Return
C
C     PROJECT THE Pcc 
C
C     MODIFIED PROJECTION
C
c      Do I=NBasis,NBasis-NELE+2,-1
c      Pcc(I)=One
c      EndDo
C
      If(Error.Le.1.D0) Then
      Do I=1,NBasis
      If(DEN(I).Lt.Small.And.Pcc(I).Lt.Zero) Then
      Pcc(I)=Abs(Pcc(I))
      EndIf
      EndDo 
      EndIf
C
c      Call PrOcc2(XMiu,One,Pcc,NBasis-NELE+1)
C
c      Call SortOcc(Pcc,URe,NBasis)
c      Do I=NELE+2,NBasis
c      Pcc(I)=Zero
c      EndDo
c      Call PrOcc2(XMiu,XELE,Pcc,NELE+1)

C     ORIGINAL PROJECTION
C
      Call PrOcc2(XMiu,XELE,Pcc,NBasis)
C
C     OBTAIN A NEW MATRIX GamSR
C
      Call VecTr(GamSR,Pcc,URe,NBasis)
C
C     GET NEW OCCUPANCIES
C
      Do I=1,NBasis
      Occ(I)=Pcc(I)**2
      EndDo
C
      Return
      End

*Deck PrOcc1
      Subroutine PrOcc1(XMiu,XNorm,Occ,NB)
C
C     PROJECTION OF Occ VECTOR FOLLOWING
C     Cances&Pernal JCP 128, 134108 (2008), Eq.(17)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Occ(NB)
C
      Parameter (Zero=0.D0,One=1.D0,Two=2.D0,Ten=10.D0,Mx=5)
C
C     CHECK IF XMiu=0 IS OK
C
      Sum=Zero
      Sum2=Zero
      Do I=1,NB
      Sum2=Sum2+Occ(I)
      If(Occ(I).Ge.Zero) Sum=Sum+Occ(I)
      EndDo
C
      If(Sum.Eq.Sum2.And.Sum.Eq.XNorm) Then
      XMiu=Zero
      GoTo 444
      EndIf
C
      XMiu1=Zero
C
      If(XSum1(XNorm,XMiu1,Occ,NB).Gt.Zero) Then
      XMiu2=-0.01
      ICount=0
  333   If(XSum1(XNorm,XMiu2,Occ,NB).Ge.Zero) Then
      XMiu2=XMiu2*Two
      ICount=ICount+1
      If(ICount.Gt.Mx) Stop 'Fatal error2 in PrOcc1'
      GoTo 333
      EndIf
      EndIf
C
      If(XSum1(XNorm,XMiu1,Occ,NB).Lt.Zero) Then
      XMiu2= 0.01
      ICount=0
  336   If(XSum1(XNorm,XMiu2,Occ,NB).Le.Zero) Then
      XMiu2=XMiu2*Two
      ICount=ICount+1
      If(ICount.Gt.Mx) Stop 'Fatal error3 in PrOcc1'
      GoTo 336
      EndIf
      EndIf
C
      XMiu=RtBis(XSum1,XMiu1,XMiu2,XNorm,Occ,NB)
C
  444   Do I=1,NB
C
      If(Occ(I)+Xmiu.Gt.Zero) Then
      Occ(I)=Min(Occ(I)+XMiu,One)
      Else
      Occ(I)=Zero
      EndIf
C
      EndDo
C
      Return
      End

*Deck PrOcc2
      Subroutine PrOcc2(XMiu,XNorm,Pcc,NB)
C
C     FIND XMiu SO THAT THE SUM OF THE PROJECTED Occ IS .LE. XNorm 
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Pcc(NB)
C
      Parameter (Zero=0.D0,One=1.D0,Two=2.D0,Ten=10.D0,Mx=5)
C
C     ORIGINAL PROJECTION
C
C     CHECK IF XMiu=0 IS OK
C
      Sum=Zero
      Do I=1,NB 
      If(Pcc(I).Ge.Zero) Sum=Sum+Min(One,Pcc(I)**2)
      EndDo
C
      If(Sum.Le.XNorm) Then
      XMiu=Zero
      GoTo 444 
      EndIf
C
      XMiu1=Zero
      If(XSum2(XNorm,XMiu1,Pcc,NB).Le.Zero)Stop 'Fatal error1 in PrOcc2'
C
      ICount=0
      XMiu2=Ten
  333 If(XSum2(XNorm,XMiu2,Pcc,NB).Ge.Zero) Then
      XMiu2=XMiu2*Two
      ICount=ICount+1
      If(ICount.Gt.Mx) Stop 'Fatal error2 in PrOcc2'      
      GoTo 333 
      EndIf
C
      XMiu=RtBis(XSum2,XMiu1,XMiu2,XNorm,Pcc,NB)
C
  444 Do I=1,NB
C
      If(Pcc(I)/(One+Xmiu).Gt.Zero) Then
      Pcc(I)=Min(Pcc(I)/(One+XMiu),One)
      Else
      Pcc(I)=Zero
      EndIf
C
      EndDo
C
      Return
      End

*Deck XSum1
      Real*8 Function XSum1(XNorm,XMiu,Occ,NB)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Occ(NB)
C
      Parameter (Zero=0.D0,One=1.D0)
C
      XSum1=Zero
      Do I=1,NB
      If(Occ(I)+Xmiu.Gt.Zero)
     $  XSum1=XSum1+Min(Occ(I)+XMiu,One)
      EndDo
C
      XSum1=XSum1-XNorm
C
      Return
      End

*Deck XSum2
      Real*8 Function XSum2(XNorm,XMiu,Pcc,NB)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Pcc(NB)
C
      Parameter (Zero=0.D0,One=1.D0)
C
      XSum2=Zero
      Do I=1,NB
      If(Pcc(I)/(One+Xmiu).Gt.Zero)
     $  XSum2=XSum2+Min((Pcc(I)/(XMiu+One))**2,One)
      EndDo
C
      XSum2=XSum2-XNorm
C
      Return
      End

*Deck RtBis
      Real*8 Function RtBis(XSum,x1,x2,XNorm,Occ,NB)
C
C     Using bisection, find the root of a function known to lie between
C     x1 and x2. The root, returned as rtbis, will be refined until its 
C     accuracy is xacc.        
C
      Implicit Real*8 (A-H,O-Z)
C
      External XSum
C
      Dimension Occ(NB)
C
      Parameter (Zero=0.D0,Half=0.5D0)
C
C     Maximum allowed number of bisections. 
C     
      Parameter (xacc=1.D-12,JMAX=50)
C     
      fmid=XSum(XNorm,x2,Occ,NB)
      f=XSum(XNorm,x1,Occ,NB)
C     
      if(f*fmid.ge.zero) stop 'root must be bracketed in rtbis'
C     
      if(f.lt.zero) then
C
C     Orient the search so that f>0 lies at x+dx. 
C     
      rtbis=x1
      dx=x2-x1
C     
      else
C     
      rtbis=x2
      dx=x1-x2
C     
      endif
C     
      do j=1,JMAX
C
C     Bisection loop. 
C     
      dx=dx*half 
      xmid=rtbis+dx 
      fmid=XSum(XNorm,xmid,Occ,NB)
C     
      if(fmid.le.zero) rtbis=xmid 
      if(abs(dx).lt.xacc.or.fmid.eq.zero) return
C      
      enddo
C     
      stop 'too many bisections in rtbis!'
C     
      return
      end
