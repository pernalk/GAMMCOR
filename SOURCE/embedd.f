*Deck EmbedDFT
      Subroutine EmbedDFT
     $ (ETot,Title,URe,Occ,XKin,XNuc,ENuc,UMOAO,
     $ DipX,DipY,DipZ,TwoEl,NBasis,NInte1,NInte2,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NSymMO,NGrid,
     $ IModG,NGOcc)
C
C     APSG-in-DFT FUNCTIONAL OPTIMIZATION
C
C     !!! XKin CONTAINS BOTH KINETIC AND EL-N CONTRIBUTIONS 
C     !!! V_XC WILL BE STORED IN XNuc 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
C
      Include 'commons.inc'
C
      Parameter(Zero=0.0D0,Half=0.50D0,One=1.0D0,Two=2.0D0,Four=4.0D0)
C
      Parameter (MXIT=50, ETol=1.D-5, DampG=0.2D0)
C
      Dimension URe(Nbasis,NBasis),Occ(NBasis),XKin(NInte1),
     $ XNuc(*),TwoEl(NInte2)
     $ OrbGrid(NBasis,NGrid),
     $ WGrid(NGrid),UMOAO(Nbasis,NBasis),NSymMO(NBasis)
C
      Dimension OrbXGrid(*),OrbYGrid(*),OrbZGrid(*),
     $ DipX(*),DipY(*),DipZ(*)
C
C     LOCAL ARRAYS
C
      Dimension VXC(NInte1),Gamma(NInte1),
     $ XOne(NInte1),GammaOld(NInte1),Work(NBasis),
     $ NSymNO(NBasis),OccSav(NBasis),UReSav(NBasis,NBasis),
     $ IGemSav(NBasis),IODFTSav(NBasis),OccDFT(NBasis)
C
      OccDFT(1:NBasis)=Zero
      Do I=1,NBasis
      If(IOrbDFT(I).Eq.1) OccDFT(I)=One
      EndDo
C
      Call VecTr(Gamma,OccDFT,URe,NBasis)
C
C     BEGINNING OF THE MACRO CYCLES (V_XC IS COMPUTED WITH DFT DENSITY KEPT FROZEN IN A MACRO CYCLE)
C
      EOld=Zero
C
      Do II=1,MXIT
C
      OccDFT(1:NBasis)=Zero
      Do I=1,NBasis
      If(IOrbDFT(I).Eq.1) OccDFT(I)=One
      EndDo
C
      Call CpyM(UReSav,URe,NBasis)
      Call CpyV(OccSav,Occ,NBasis)
      Do I=1,NBasis
      IGemSav(I)=IGem(I)
      IODFTSav(I)=IOrbDFT(I)
      EndDo
C
C     AVERAGE THE DENSITY MATRIX
C
      If(II.Gt.1) Then
C
      Do L=1,NInte1
      GammaOld(L)=Gamma(L)
      EndDo
C 
      Call VecTr(Gamma,OccDFT,URe,NBasis)
C
      Do L=1,NInte1
      Gamma(L)=DampG*Gamma(L)+(One-DampG)*GammaOld(L)
      EndDo
C
      Call CpySym(URe,Gamma,NBasis)
      Call Diag8(URe,NBasis,NBasis,OccDFT,Work)
C     sometimes very small occ's come out negative so zero them
      Do I=1,NBasis
      OccDFT(I)=Abs(OccDFT(I))
      EndDo
      Call SortOcc(OccDFT,URe,NBasis)
C
C     END OF AVERAGING OF GAMMA
c     endif of If(II.Gt.1) 
      EndIf
C
C     COMPUTE THE XC POTENTIAL IN THE MO REP FOR THE DENSITY U*OccDFT*U+
C
      Call EPotXC(EnXC,VXC,OccDFT,URe,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NSymMO,TwoEl,NGrid,NInte1,NInte2,NBasis)
C
      XVXC=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      Fac=Four
      If(I.Eq.J) Fac=Two
      XVXC=XVXC+Fac*Gamma(IJ)*VXC(IJ)
      EndDo
      EndDo
C
C     BEGIN THE ITERATION WITH THE SAVED Occ AND URe
C
      Call CpyM(URe,UReSav,NBasis)
      Call CpyV(Occ,OccSav,NBasis)
      Do I=1,NBasis
      IGem(I)=IGemSav(I)
      IOrbDFT(I)=IODFTSav(I)      
      EndDo
C
      Do I=1,NInte1
      XNuc(I)=VXC(I)
      EndDo
C
      Call OPTNORB(ETot,URe,Occ,XOne,XNuc,DipX,DipY,DipZ,TwoEl,UMOAO,
     $ NSymMO,Title,NBasis,NInte1,NInte2,NGem,NGOcc,IModG)
C
      Call ReWr(1,Occ,URe,Title,NBasis)
      ENew=ETot+EnXC-XVXC
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
      Write(6,'(X,''APSG-in-DFT ITER'',I3,2X,''ENERGY'',F14.8,2X,
     $ ''ENE DIFF '',E10.3)') II,Enew,ENew-EOld
      Write(6,'(X,''****************************************'',
     $            ''***************************************'',/)')
C
      Err=Abs((EOld-ENew))
C
      If If(Err.Lt.ETol) Return
C
      EOld=ENew
C
C     End Of Loop II=1, MXIT
      EndDo
C
      Write(6,'(/,10X,''NO CONVERGENCE IN DFT/DMFT !!!'')')
C
      Return
      End

*Deck AssignDFT
      Subroutine AssignDFT(URe,Occ,NSymNO,NGOcc,NBasis)
C     
      Implicit Real*8 (A-H,O-Z)
C     
C     DIVIDE A SET OF OCCUPIED LOCALIZED MO ORBITALS INTO THOSE BELONGING
C     TO THE DFT DOMAIN AND THE REST 
C
      Include 'commons.inc'
C     
      Parameter(Zero=0.0D0,One=1.0D0,Two=2.D0)
C    
      Dimension URe(NBasis,NBasis),Occ(NBasis)
C    
C    
      Do I=1,NBasis
      IOrbDFT(I)=0
      EndDo
C
C     H2O in DZ basis
C
      IOrbDFT(1)=1
      IOrbDFT(2)=1
      IOrbDFT(5)=1
      GoTo 555
C
  555 Continue
C
C     SET NGOcc
C
      NGOcc=0
      Do I=1,NBasis
C
      If(IOrbDFT(I).Eq.1) Then
      NGOcc=NGOcc+1
      Occ(I)=Two
      EndIf
C 
      EndDo
C
C     SORT ORBITALS SO THAT THOSE WITH Occ=2 (BELONGING TO THE DFT DOMAIN) ARE PUT ON TOP
C
      Call SortAll(Occ,URe,NSymNO,NBasis)
C
C     LABEL ORBITALS WITH Occ=2 AND SET THE OCCUPATION TO 1
C
      Do I=1,NBasis
      If(Occ(I).Eq.Two) Then
      Occ(I)=One
      IOrbDFT(I)=1
      Else
      IOrbDFT(I)=0
      EndIf
      EndDo
C
      Return
      End
