*Deck INTERPA
      Subroutine INTERPA(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  Title,NBasis,NInte1,NInte2,NDim,
     $  NGem,NGOcc)
C
      use abmat
      use abfofo
      use ab0fofo
C
C      Subroutine INTERPA(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
C     $  Title,OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,
C     $  NGem,IAPSG,ISERPA,QMAX,NGOcc,Small)
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
C     DENSITY MATRIX ELEMENTS
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1)
C     $ XOne(NInte1),OrbGrid(NBasis,NGrid),WGrid(NGrid)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ IndX(NDim),IndN(2,NDim),
     $ IndP(NBasis,NBasis),
     $ NSymMO(NBasis),NSymNO(NBasis),
     $ EigVecR(NDim*NDim),Eig(NDim),
     $ IndAux(NBasis),
     $ ECorrG(NGem), EGOne(NGem),IPair(NBasis,NBasis)
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C
C     IFlSnd  = 1 - run AC0 (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0
C
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation 
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     IFlFrag1 = 1 : run embedding: call FragEcorr
C                0 : do not call FragEcorr
C              USED IN EneERPA 
C
C     IFl12   = 0  : A+,A- MATRICES NOT TRUNCATED
C             = 1  : A+,A- MATRICES ARE TRUNCATED TO INCLUDE ONLY INDICES
C                          OF ORBITALS BELONGING TO GEMINALS NEEDED FOR 
C                          A GIVEN ONE- TWO- THREE- OR FOUR-BODY INTERACTION
C             = 2  : OLD VARIANT OF EERPA (PARTITIONING INTO MONOMERS AND IFlag12 IS SET TO 1)
C             USED IN FragEcorr and OneTwoBody
C             
C             IFlag12 IS ONLY IN EFFECT WHEN IFlFrag1=1
C
C     MH: SWITCHING IT OFF!!!
c      IFl12=1
C
C     RPAX: IFlCore=1, IFlFrag1=0
C     ERPA without core orbitals: IFlCore=0, IFlFrag1=0
C     ERPA and GVB+1,2-body with/without core orbs: IFlCore=1/0, IFl12=1, IFlFrag1=1
C
C
      If(IFlFrag1.Eq.1.And.IFlAC+IFlSnd.Eq.2) 
     $ Stop 'If IFlFrag1.Eq.1 one of the IFlAC, IFlSnd flags 
     $ must be zero!!!'
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1) 
     $ Write(6,'(/," *** ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlSnd.Eq.1.And.ICASSCF.Ne.1)
     $ Write(6,'(/," *** LINEARIZED (AC0-GVB) AC CALCULATIONS ***",/)') 
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1) 
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1) 
     $ included in ERPA correlation ***",/)')
      EndIf
C
      Write(6,'(X," IFlFrag1 ",I4,/)') IFlFrag1
      If(IFlFrag1.Eq.1) Write(6,'(X," IFl12    ",I4,/)') IFl12 
C
      Do I=1,NBasis
      IAuxGem(I)=IGem(I)
      EndDo
C
      Do I=1,NGem
      Do J=1,NGem
c      IConnect(I,J)=1
      IConnect(I,J)=0
      EndDo
      EndDo
C
      NActOrb=1
C
C     CONSTRUCT A LOOK-UP TABLE
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      If(NActOrb.Ne.0) Then
      ICount=0
C
      If(ICASSCF.Eq.0) Then
C
c      ThrAct=0.990
      Write(6,'(/,X," Threshold for Active Orbitals in GVB: ",E14.4)')
     $ ThrGemAct
      Do I=1,NELE
      If(Occ(I).Lt.ThrGemAct) Then
      IndAux(I)=1
      Write(6,'(/,X," Active GVB Orbital: ",I4,E14.4)')
     $ I, Occ(I)
      IndAux(IFindG(I))=1
      Write(6,'(X," Active GVB Orbital: ",I4,E14.4)')
     $ IFindG(I), Occ(IFindG(I))
      ICount=ICount+2
      EndIf
      EndDo
C
      ElseIf(ICASSCF.Eq.1) Then
C
      Do I=1,NBasis
C
      If(Occ(I).Lt.One.And.Occ(I).Ne.Zero) Then
      IndAux(I)=1
      Write(6,'(X," Active Orbital: ",I4,E14.4)') I, Occ(I)
      ICount=ICount+1
      EndIf
      EndDo
C
C     If(ICASSCF.Eq.0)
      EndIf
C     If(NActOrb.Ne.0) Then 
      EndIf
C
      Write(6,'(/,X," In INTERPA: Active Orbitals ",I4,/)')ICount
C
C     CONSTRUCT LOOK-UP TABLES
C
      Write(6,'(2x,a,2e15.5)') 'Threshold for quasi-degeneracy ',
     $ ThrSelAct
C
      IPair(1:NBasis,1:NBasis)=0
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
C
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate active degenerate orbitals if from different geminals
      If((ICASSCF.Eq.0.And.(IGem(I).Ne.IGem(J)).And.(IndAux(I).Eq.1)
     $ .And.(IndAux(J).Eq.1).And.
     $ (Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct))
     $.Or.
     $ (ICASSCF.Eq.1.
     $ .And.(IndAux(I).Eq.1).And.(IndAux(J).Eq.1) 
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct)) ) Then
C
      Write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals  
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      IPair(I,J)=1
      IPair(J,I)=1    
C
      EndIf
C
      EndIf
C
c     If(IndAux(I)+IndAux(J).Ne.0 ...
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
      Write(6,'(/,2X,"Total number of pairs:",I6)')NDim
      Write(6,'(2X,"Reduced to:",I6)')Ind
      Write(6,'(2X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
C     CALL AC If IFlaAC=1 OR IFlSnd=1
C
      If(IFlAC.Eq.1.Or.IFlSnd.Eq.1) Then
      Call ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDim,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      Return
      EndIf
C
C     CALCULATE THE A+B AND A-B MATRICES
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
c      If (ISERPA.Eq.1) Then
c      Write(6,'(/,X," ***** SERPA ENERGY CALCULATIONS *****")') 
c      ElseIf(ISERPA.Eq.0) Then
c      Write(6,'(/,X," ***** ERPA ENERGY CALCULATIONS *****")')
c      ElseIf(ISERPA.Eq.2) Then
c      Write(6,'(/,X," ***** PINO ENERGY CALCULATIONS *****")') 
c      EndIf
C
      ACAlpha=One
      If(ICASSCF.Eq.0) Then
C
      If(ITwoEl.Eq.1) Then

      If (ITrpl .Eq. 0) Then
         ! Singlet GVB response
          Call ACABMAT0(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $                  NBasis,NDim,NInte1,NInte2,NGem,ACAlpha,1)
      ElseIf(ITrpl .Eq. 1) Then
         ! Triplet GVB response
          Call AB_T_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $                IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
      End If
C
C      print*, 'ABPLUS-Ka',norm2(ABPLUS(1:NDim**2))
C      print*, 'ABMIN -Ka',norm2(ABMIN(1:NDim**2))

      ElseIf(ITwoEl.Eq.2) Then

      If(IFlFrag1.Eq.1) Then
C        this should be avoided somehow...
         Call LookUp_mithap(Occ,IndAux,IndP,IndN,IndX,NDimX,NDim,NBasis)
      EndIf
C
      Call ACABMAT0_mithap(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO',0,ACAlpha,1)

C      print*, 'ABPLUS-my',norm2(ABPLUS(1:NDim**2))
C      print*, 'ABMIN -my',norm2(ABMIN(1:NDim**2))

      Call EneGVB_FFFF(ETot,URe,Occ,CICoef,XOne,
     $                 IGem,IndN,NBasis,NInte1,'TWOMO',NDimX,NGem)

      ElseIf(ITwoEl.Eq.3) Then
C      
      If(IFlFrag1.Eq.1) Then
         Call LookUp_mithap(Occ,IndAux,IndP,IndN,IndX,NDimX,NDim,NBasis)
      EndIf
C
      Call ACABMAT0_FOFO(ABPLUS,ABMIN,URe,Occ,XOne,
     $            IndN,IndX,IGem,CICoef,
     $            NActive,NELE,NBasis,NDim,NDimX,NInte1,NGem,
     $            'TWOMO','FFOO','FOFO',0,ACAlpha,1)

      print*, 'ABPLUS-my',norm2(ABPLUS)
      print*, 'ABMIN -my',norm2(ABMIN)

      Call EneGVB_FOFO(NActive,NELE,ETot,URe,Occ,CICoef,XOne,
     $                 IGem,IndN,NBasis,NInte1,'FOFO',NDimX,NGem)

      EndIf
C
      ElseIf(ICASSCF.Eq.1) Then
C
c       Call RDMFT_AB(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
c     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
c      NAct=NELE
c      Do I=1,NELE
c      If(Occ(I).Eq.One) NAct=NAct-1
c      EndDo
c      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
c     $ NAct*2
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS 
c      Call Gamma2_AB(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
c     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
      ETot=ECASSCF
C
      EndIf
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
C     COMPUTE ONE- AND TWO-BODY INTERGEMINAL CORRELATION ENERGIES
C
      EOneTot=Zero
C
      If(IFlFrag1.Eq.1) Then
C
      Write(6,'(/,X,"***** EMBEDDING ERPA (EERPA) CALCULATIONS *****")')
C
c      Write(6,'(/,X," Entering OneTwoBody with IFl12 = ",I4)'),IFl12
c      Call OneTwoBody(ETot,ENuc,ECorrTot,EGOne,EigVecR,Eig,ABPLUS,ABMIN,
c     $ Occ,TwoNO,IndAux,NBasis,NInte1,NInte2,NDim,NGem,NGOcc,IFl12)
c      NFrag=NGem-1
c      Return
C
      Write(6,'(/,X," Entering FragEcorr with IFlCore = ",I4)'),IFlCore
      Write(6,'(X," Entering FragEcorr with IFl12 = ",I4,/)'),IFl12
      Call FragEcorr(ETot,ENuc,ECorrTot,EGOne,EigVecR,Eig,ABPLUS,ABMIN,
     $ UNOAO,Occ,TwoNO,URe,XOne,IndAux,NBasis,NInte1,NInte2,NDim,NGem,
     $ NGOcc,IFl12,NFrag)
C
      If(ITwoEl.ne.1) Call DelInts(ITwoEl)
C      
      Return
C
      IFlFrag1=1
      IStart=1
      EOneTot=Zero
      Do IGG=IStart,NFrag
      EOneTot=EOneTot+EGOne(IGG)
      EndDo
C
      EndIf
C
C     REDUCE DIMENSIONS OF THE MATRICES TO TAKE INTO ACCOUNT ONLY
C     NDimX ELEMENTS OF dGamma_ij 
C
      Write(6,'(/,X,"***************************** ")')
      If(ICASSCF.Eq.1) Then
      Write(6,'(  X,"*** ERPA-CAS CALCULATIONS *** ")')
      Else
      Write(6,'(  X,"*** ERPA-GVB CALCULATIONS *** ")')
      EndIf
      Write(6,'(  X,"***************************** ")') 
C
C     REDUCE THE MATRICES
C
      NDimN=NBasis
C
      If(ITwoEl.Eq.1) Then
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
C     FFFF and FOFO: ABMATs ALREADY TRUNCATED
C
      EndIf
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(ISERPA.Eq.0) Then
C
      Call ERPASYMM(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)

      Write(6,'(/," *** Computing ERPA energy *** ",/)')
C
      If(ITwoEl.Eq.1) Then

      ECorr=EOneTot
      Call EneERPA(ETot,ECorr,ENuc,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)

      ElseIf(ITwoEl.Eq.2) Then

      Call EneERPA_FFFF(ETot,ECorr,ENuc,EigVecR,Eig,Occ,CICoef,
     $                  IGem,IndN,NDimX,NBasis,'TWOMO')

      ElseIf(ITwoEl.Eq.3) Then

      Call EneERPA_FOFO(ECorr,EigVecR,Eig,Occ,CICoef,
     $                  IGem,IndN,NDimX,NELE+NActive,NBasis,'FOFO')

      ECorr = 0.5d0*ECorr
      Write(LOUT,'(1x,a,3f15.8)') 'EGVB+ENuc, Corr, ERPA-GVB',
     $      ETot+ENuc,ECorr,ETot+ENuc+ECorr

      EndIf

C     DELETE MO INTEGRALS 
      If(ITwoel.Eq.2) Then
        Open(newunit=iunit,file='TWOMO',status='OLD')
        Close(iunit,status='DELETE')
      ElseIf(ITwoel.Eq.3) Then
        Open(newunit=iunit,file='FFOO',status='OLD')
        Close(iunit,status='DELETE')
        Open(newunit=iunit,file='FOFO',status='OLD')
        Close(iunit,status='DELETE')
      EndIf

      If (ITrpl .Eq. 0) Then
      Write(6,'(/,
     $ " *** ERPA-GVB Singlet Excitation Energies (a.u., eV) *** ")')
      ElseIf (ITrpl .Eq. 1) Then
      Write(6,'(/,
     $ " *** ERPA-GVB Triplet Excitation Energies (a.u., eV) *** ")')
      End If
C     the purpose of sorting is only to print a few highest (sorted) eigenvectors
      Do I=1,20
      Write(6,'(I4,4X,2F16.8)') I,Eig(I),27.211*Eig(I)
      EndDo
C
c      Write(6,'(/," *** Computing ERPA 2-RDM *** ")')
cC
c      Call RDM2ERPA(EigVecR,Eig,ABMIN,TwoNO,NInte2,IndN,Occ,Title,
c     $ NBasis,NDimX,NGem,NDim)
C
      EndIf
C
      Return
      End

*Deck EneERPA
      Subroutine EneERPA(ETot,ECorr,ENuc,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      Implicit Real*8 (A-H,O-Z)
C    
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
c      Parameter(SmallE=1.D-3,BigE=5.D2)
      Parameter(SmallE=1.D-3,BigE=1.D6)
C
C     ONLY EXCITATIONS SMALLER THAN BigE AND GREATER THAN SmallE ARE INCLUDED 
C     
      Include 'commons.inc'
C     
      Dimension EigVecR(NDimX*NDimX),
     $ Eig(NDimX),URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),C(NBasis),EInterGF2(NGem,NGem,NGem,NGem),
     $ EInterGF3(NGem,NGem,NGem,Ngem),
     $ LabelEig(NDimX),Skipped(NDimX)
C
      ESave=ETot
C
C     it may happen that a certain number of statest with the smallest eigenvalues should be removed
C     (nearly zero eigenvalues of the hessian) set here a number of states to be removed
C     and give the label to such eigenvalues
C
      NStates=0
C
      Do I=1,NDimX
      LabelEig(I)=0
      EndDo
      Do I=1,NStates
      XMin=1.D6
      Do J=1,NDimX
      If(Eig(J).Lt.XMin.And.LabelEig(J).Eq.0) Then
      XMin=Eig(J)
      JMin=J
      EndIf
      EndDo
      LabelEig(JMin)=1
      Write(*,*)'discard', eig(jmin)
      EndDo
C
      Do I=1,NBasis
      C(I)=CICoef(I)
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     ONE_ELECTRON PART + HARTREE + EXCHANGE + np(1-np) <pp|pp>
C     COMPUTE THE APSG ENERGY
C
      ETot=Zero
      EAPSG=Zero
      EOne=Zero
      EIntraGem=Zero
      EInterCoul=Zero
      EInterExch=Zero
C
      Do I=1,NGem
      Do J=1,NGem
      Do K=1,NGem
      Do L=1,NGem 
      EInterGF2(I,J,K,L)=Zero
      EInterGF3(I,J,K,L)=Zero
      EndDo
      EndDo
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
C
      II=(I*(I+1))/2
      ETot=ETot+Occ(I)*HNO(II)
     $ -Half*Occ(I)*(One-Occ(I))*TwoNO(NAddr3(I,I,I,I))
C
      EAPSG=EAPSG+Two*Occ(I)*HNO(II)
      EOne=EOne+Two*Occ(I)*HNO(II)
C
      Do J=1,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
      ETot=ETot+FacIJ*Occ(I)*Occ(J)*(TwoNO(NAddr3(I,I,J,J))
     $ -Half*TwoNO(NAddr3(I,J,I,J)))
C
      If(IGem(I).Eq.IGem(J)) Then
      EAPSG=EAPSG+FacIJ*C(I)*C(J)*TwoNO(NAddr3(I,J,I,J))
      EIntraGem=EIntraGem+FacIJ*C(I)*C(J)*TwoNO(NAddr3(I,J,I,J))
      Else
      EAPSG=EAPSG+FacIJ*Occ(I)*Occ(J)*(Two*TwoNO(NAddr3(I,I,J,J))
     $ -TwoNO(NAddr3(I,J,I,J)))
      EInterCoul=EInterCoul
     $ +FacIJ*Occ(I)*Occ(J)*Two*TwoNO(NAddr3(I,I,J,J))
      EInterExch=EInterExch
     $ -FacIJ*Occ(I)*Occ(J)*TwoNO(NAddr3(I,J,I,J))
      EndIf
C
      EndDo
      EndDo
C
      ETot=Two*ETot
      If(ICASSCF.Ne.1) Then
      Write(6,'(" One-electron energy",25X,F17.8)')EOne
      Write(6,'(" GVB intra-gem electron interaction",10X,F17.8)')
     $ EIntraGem
      Write(6,'(" GVB inter-gem Coulomb interaction",11X,F17.8)')
     $ EInterCoul
      Write(6,'(" GVB inter-gem exchange interaction",10X,F17.8)')
     $ EInterExch
      Write(6,'(" GVB electron interaction energy",13X,F17.8)')
     $ EIntraGem+EInterCoul+EInterExch
      Write(6,'(" Total GVB energy",27X,F18.8)')
     $ EOne+EIntraGem+EInterCoul+EInterExch
      Write(6,'(" Nuclear repulsion",26X,F18.8)')ENuc
      Write(6,'(" E_GVB+ENuc",33X,F18.8)')
     $ EOne+EIntraGem+EInterCoul+EInterExch+ENuc
      EndIf
C
C     ADD CONTRIBUTIONS FROM THE ERPA VECTORS
C     COMPUTE INTERGEMINAL CORRELATION, EInterG
C
      EInterG=Zero
      EDisp=Zero
      ECT=Zero
C
      EIntra=Zero
      EIntraFrag=Zero
      EAll=Zero
C
      EVirt=Zero
      EnonVirt=Zero
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
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      DispPR=Zero
      CTPR=Zero
      If(Occ(IP)+Occ(IR).Ne.Zero) Then
      If(IGem(IP).Eq.IGem(IR).Or.(Occ(IP)*Occ(IR).Eq.Zero))        
     $ DispPR=One
      If(IGem(IP).Ne.IGem(IR).And.(Occ(IP)*Occ(IR).Ne.Zero))
     $ CTPR=One
      EndIf
C
      DispQS=Zero
      CTQS=Zero
      If(Occ(IQ)+Occ(IS).Ne.Zero) Then
      If(IGem(IQ).Eq.IGem(IS).Or.(Occ(IQ)*Occ(IS).Eq.Zero))
     $ DispQS=One
      If(IGem(IQ).Ne.IGem(IS).And.(Occ(IQ)*Occ(IS).Ne.Zero))
     $ CTQS=One
      EndIf
C
      ISkippedEig=0
      SumY=Zero
      Do K=1,NDimX
      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE.And.LabelEig(K).Eq.0) Then
      SumY=SumY+EigVecR((K-1)*NDimX+I)*EigVecR((K-1)*NDimX+J)
      Else
      ISkippedEig=ISkippedEig+1
      Skipped(ISkippedEig)=Eig(K) 
      EndIf
      EndDo
C
      Aux=Two*(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      AuxInterG=Zero
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).And.
     $ IGem(IP).Ne.IGem(IQ)) AuxInterG=Aux
C
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux=Aux
     $ -Occ(IP)*(One-Occ(IS))-Occ(IS)*(One-Occ(IP))
C
      EAll=EAll+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      IntraG=0
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ)) Then
      IntraG=1
      EIntra=EIntra+Aux*TwoNO(NAddr3(IP,IR,IQ,IS)) 
      EndIf
C
      If(IFlFrag1.Eq.1.And.IntraG.Eq.0) Then
      If(IAuxGem(IP).Eq.IAuxGem(IR).And.IAuxGem(IQ).Eq.IAuxGem(IS).
     $ And.IAuxGem(IP).Eq.IAuxGem(IQ)) EIntraFrag=EIntraFrag
     $ +Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
C 
C     DEFINED ONLY FOR GVB      
C
      If(IGVB.Eq.1.And.IGem(IR).Ne.IGem(IS)) Then
      Disp=DispPR*DispQS
      If(Disp.Ne.Zero) then
      EDisp=EDisp+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      endif
      EndIf
C
      If(IGVB.Eq.1.And.IGem(IR).Eq.IGem(IS).And.
     $ Occ(IP)+Occ(IQ).Ne.Zero) 
     $ Then
      CT=CTQS*CTPR+DispPR*CTQS+CTPR*DispQS
      If(CT.Ne.Zero) then
      ECT=ECT+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      endif
      EndIf
C
      ETot=ETot+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EInterG=EInterG+AuxInterG*TwoNO(NAddr3(IP,IR,IQ,IS))
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IF FRAGMENTS ARE "ON" COMPUTE CONTRIBUTIONS TO CORRELATION FROM 
C     PAIRS OF GEMINALS ON TWO DIFFERENT FRAGMENTS
C
      If(IFlFrag1.Eq.1) Then
C
      FVirt=Occ(IP)*Occ(IR)*Occ(IQ)*Occ(IS)
      IGP=IGem(IP)
      IGR=IGem(IR)
      IGQ=IGem(IQ)
      IGS=IGem(IS)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   
      IFP=IAuxGem(IP)
      IFR=IAuxGem(IR)
      IFQ=IAuxGem(IQ)
      IFS=IAuxGem(IS)       
C
      ISavGP=IGP
      ISavGR=IGR
      ISavGQ=IGQ
      ISavGS=IGS
      IPR=NGem*(ISavGR-1)+ISavGP
      IQS=NGem*(ISavGS-1)+ISavGQ
C
      Call IBodySort(IBdyG,IGP,IGR,IGQ,IGS)
      Call IBodySort(IBdyF,IFP,IFR,IFQ,IFS)
C
      If(FVirt.Ne.Zero) Then
C 2-body
      If(IBdyG.Eq.2.And.IBdyF.Eq.2) Then
      If(IPR.Ge.IQS) Then
      EInterGF2(ISavGP,ISavGR,ISavGQ,ISavGS)=
     $ EInterGF2(ISavGP,ISavGR,ISavGQ,ISavGS)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      Else
      EInterGF2(ISavGQ,ISavGS,ISavGP,ISavGR)=
     $ EInterGF2(ISavGQ,ISavGS,ISavGP,ISavGR)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
      EndIf
c 3-body
      If(IBdyG.Eq.3.And.IBdyF.Eq.2) Then
      If(IPR.Ge.IQS) Then
      EInterGF3(ISavGP,ISavGR,ISavGQ,ISavGS)=
     $ EInterGF3(ISavGP,ISavGR,ISavGQ,ISavGS)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      Else
      EInterGF3(ISavGQ,ISavGS,ISavGP,ISavGR)=
     $ EInterGF3(ISavGQ,ISavGS,ISavGP,ISavGR)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
      EndIf
C
      Else
C
      If(IBdyG.Eq.3.And.IBdyF.Eq.3) Then
c IGS is a geminal with virtual orbitals do use IGR,IGQ
      If(IPR.Ge.IQS) Then
      EInterGF2(ISavGP,ISavGR,ISavGQ,ISavGS)=
     $ EInterGF2(ISavGP,ISavGR,ISavGQ,ISavGS)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      Else
      EInterGF2(ISavGQ,ISavGS,ISavGP,ISavGR)=
     $ EInterGF2(ISavGQ,ISavGS,ISavGP,ISavGR)
     $ +Half*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
      EndIf
C
      EndIf
C
C     If(IFlFrag1.Eq.1)
      EndIf
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      If(Occ(IP)*Occ(IR)*Occ(IQ)*Occ(IS).Eq.Zero) Then
C
      IntraG=0
      If(Occ(IP).Eq.Zero.And.IGem(IR).Eq.IGem(IS).
     $ And.IGem(IQ).Eq.IGem(IS)) Then
      IntraG=1
      EVirt=EVirt+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))      
      EndIf
      If(Occ(IQ).Eq.Zero.And.IGem(IR).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IS)) Then
      IntraG=1
      EVirt=EVirt+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
      If(Occ(IQ).Eq.Zero.And.Occ(IP).Eq.Zero.And.IGem(IR).Eq.IGem(IS))
     $ Then
      IntraG=1
      EVirt=EVirt+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf 
C
      If(IFlFrag1.Eq.1.And.IntraG.Eq.0) Then
      If(Occ(IP).Eq.Zero.And.IAuxGem(IR).Eq.IAuxGem(IS).
     $ And.IAuxGem(IQ).Eq.IAuxGem(IS))
     $ EIntraFrag=EIntraFrag+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      If(Occ(IQ).Eq.Zero.And.IAuxGem(IR).Eq.IAuxGem(IS).
     $ And.IAuxGem(IP).Eq.IAuxGem(IS))
     $ EIntraFrag=EIntraFrag+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      If(Occ(IQ).Eq.Zero.And.Occ(IP).Eq.Zero.And.
     $ IAuxGem(IR).Eq.IAuxGem(IS))
     $ EIntraFrag=EIntraFrag+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndIf
C
c     else of Occ(IP)*Occ(IR)*Occ(IQ)*Occ(IS).Eq.Zero
      Else
C
      IPR=0
      If(IGem(IP).Eq.IGem(IR)) IPR=1
      IQS=0
      If(IGem(IQ).Eq.IGem(IS)) IQS=1
      IPQ=0
      If(IGem(IP).Eq.IGem(IQ)) IPQ=1
C
      EnonVirt=EnonVirt+(One-IPR*IQS*IPQ)*Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndIf
C
C     endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
      EndIf
C
      EndDo
      EndDo
C
      Write(6,'(1X,"The number of discarded eigenvalues is",I4)')
     $  ISkippedEig
      Do II=1,ISkippedEig
      Write(6,*)'Skipped',II,Skipped(II)
      EndDo
C
c      Write
c     $ (6,'(/,1X,''EERPA + ENuc '', 49X,F15.8)')ETot+ENuc
C
      If(IGVB.Eq.0) Then
C     
       Write
     $ (6,'(1X,''EAPSG+ENuc, Corr, ERPA-APSG '',4X,3F15.8)')EAPSG+ENuc,
     $ Half*(EAll-EIntra),EAPSG+ENuc+Half*(EAll-EIntra)

       Write
     $ (6,'(1X,''EAPSG+ENuc, IGCorr, APSG+IGCorr '',3F15.8)')EAPSG+ENuc,
     $ Half*EInterG,EAPSG+ENuc+Half*EInterG
C
      Else
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      GoTo 777
      If(IFlFrag1.Eq.1) Then
C
      Write(6,'(/)')
C
      Sum2=Zero
C
      Do I=1,NGem-1
      Do J=1,I-1
      SumIJ=Zero
C
      Do IGG1=1,NGem
      Do IGG2=1,NGem
      Do IGG3=1,NGem
      Do IGG4=1,NGem
C
      II=0
      JJ=0
      If(IGG1.Eq.I.Or.IGG2.Eq.I.Or.IGG3.Eq.I.Or.IGG4.Eq.I) II=1
      If(IGG1.Eq.J.Or.IGG2.Eq.J.Or.IGG3.Eq.J.Or.IGG4.Eq.J) JJ=1
C
      If(II*JJ.Eq.1) Then
      Aux=EInterGF2(IGG1,IGG2,IGG3,IGG4)
      If(Abs(Aux).Gt.1.D-5) Then
      Write(6,'(X,''Gem1, Gem2, Gem3, Gem4, Corr_12 '',5X,4I4,F15.8)')
     $ IGG2,IGG1,IGG4,IGG3,Aux
      Sum2=Sum2+Aux
      SumIJ=SumIJ+Aux
      EndIf
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      If(Abs(SumIJ).Gt.1.D-5) 
     $ Write(6,'(X,''Corr Energy from Gem1, Gem2:'',5X,2I4,F15.8)')
     $ I,J,SumIJ
C
      EndDo 
      EndDo
C
      Write(6,'(/,X,''Sum of inter-fragment two-body corr '',5X,F15.8)')
     $ Sum2
C
      Sum3=Zero
      Do IGG1=1,NGem
      Do IGG2=1,NGem
      Do IGG3=1,NGem
      Do IGG4=1,NGem
      Aux=EInterGF3(IGG1,IGG2,IGG3,IGG4)
      If(Abs(Aux).Gt.1.D-5) Then
      Write(6,'(X,''Gem1, Gem2, Gem3, Gem4, Corr_12 '',5X,4I4,F15.8)')
     $ IGG2,IGG1,IGG4,IGG3,Aux
      Sum3=Sum3+Aux
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
      Write
     $ (6,'(/,X,''Sum of inter-fragment three-body corr '',3X,F15.8)')
     $ Sum3
C
      Write
     $ (6,'(/,X,''Sum of inter-fragment 2- and 3-body corr '',F15.8)')
     $ Sum2+Sum3
      Write(6,'(/)')
C
      EndIf
  777 Continue 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      If(ICASSCF.Eq.1) Then
      Ene=ESave
      Write
     $ (6,'(/,1X,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6X,3F15.8)')
     $ Ene+ENuc,
     $ Half*(EAll-EIntra),Ene+ENuc+Half*(EAll-EIntra)
      Return
C
      EndIf
C
      Ene=EAPSG
      Write
     $ (6,'(1X,''EGVB+ENuc, Corr, ERPA-GVB  '',6X,3F15.8)')Ene+ENuc,
     $ Half*(EAll-EIntra),Ene+ENuc+Half*(EAll-EIntra)
C
      EndIf
C
      Return
      End

C*Deck ERPAVEC
C      Subroutine ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
CC
C      Implicit Real*8 (A-H,O-Z)
CC
C      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
Cc     set small being a square of small in Deck EneERPA
C     $ Four=4.D0, Small=1.D-6)
CC
C      Include 'commons.inc'
CC
C      Dimension
C     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
C     $ HlpAB(NDimX,NDimX),
C     $ EigVecR(NDimX*NDimX),EigVecL(NDimX*NDimX),
C     $ Eig(NDimX),EigI(NDimX),Work(5*NDimX)
CC
C      Call MultpM(HlpAB,ABPLUS,ABMIN,NDimX)
CC
C      Call DGEEV('N','V',NDimX,HlpAB,NDimX,Eig,EigI,
C     $           EigVecL,NDimX,EigVecR,NDimX,Work,5*NDimX,INFO)
CC
CC     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
CC     OMEGA'S 
CC
CC     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
CC
C      Do NU=1,NDimX
C      SumNU=Zero
CC
C      If(Eig(NU).Gt.Small) Then
CC    
C      Eig(NU)=SQRT(Eig(NU))
CC
C      Do I=1,NDimX
C      Do J=1,NDimX
C      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
C     $ EigVecR((NU-1)*NDimX+I)*EigVecR((NU-1)*NDimX+J)
C      EndDo
C      EndDo
CC
C      If(SumNU.Gt.Zero) Then
C      SumNU=One/Sqrt(SumNU)
C      Else
C      SumNU=Zero
C      EndIf
CC
C      Do I=1,NDimX
C      EigVecR((NU-1)*NDimX+I)=EigVecR((NU-1)*NDimX+I)*SumNU
C      EndDo
CC
C      EndIf
Cc     enddo NU
C      EndDo
CC
C      Return
C      End

*Deck ERPAVEC
      Subroutine ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     set small being a square of small in Deck EneERPA
     $ Four=4.D0, Small=1.D-6)
C
      Integer :: DimV1,Max_NDEG
      Integer :: Space1(3,NDimX)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ HlpAB(NDimX,NDimX),
     $ EigVecR(NDimX*NDimX),EigVecL(1),
     $ Eig(NDimX),EigI(NDimX),Work(5*NDimX)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      ABPLUS(I,J)=Half*(ABPLUS(I,J)+ABPLUS(J,I))
      ABPLUS(J,I)=ABPLUS(I,J)
      ABMIN(I,J)=Half*(ABMIN(I,J)+ABMIN(J,I))
      ABMIN(J,I)=ABMIN(I,J)
      EndDo
      EndDo
C
      Call MultpM(HlpAB,ABPLUS,ABMIN,NDimX)
C
      Call DGEEV('N','V',NDimX,HlpAB,NDimX,Eig,EigI,
     $           EigVecL,1,EigVecR,NDimX,Work,5*NDimX,INFO)
C   
C     ORTHOGONALISE DEGENERATE VECTORS
      Call CREATE_SPACE(Eig,Space1,NDimX,DimV1,Max_NDEG)
      Call ORTHO_DEGVEC(EigVecR,Space1,DimV1,NDimX,Max_NDEG)
C     CALL FOR TESTS ONLY:
C      Call ZERO_DEGVEC(EigVecR,Eig,Space1,DimV1,NDimX,Max_NDEG)
C
C     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
C     OMEGA'S 
C
C     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
C
      Do NU=1,NDimX
C
      If(Abs(EigI(NU)).Gt.1.D-12) Then
C
      Write(6,'(X,"Complex ERPA Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=Zero
      EndDo
C
      EndIf
      EndDo
C
      Do NU=1,NDimX
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C   
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NDimX
      Do J=1,NDimX
      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
     $ EigVecR((NU-1)*NDimX+I)*EigVecR((NU-1)*NDimX+J)
      EndDo
      EndDo 
C      If(SumNu.Lt.Zero) Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
C     $  NU,Eig(NU),SumNU
C     CHANGE 1 
      If(SumNu.Lt.Zero) Then
      Eig(NU)=-Eig(NU)
      Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
     $  NU,Eig(NU),SumNU
      EndIf
C
C
C     If SumNU < 0 the eigenvector corresponds to deexcitation, normalize it to -1
C
c      If(SumNU.Gt.Zero) Then
c      SumNU=One/Sqrt(SumNU)
c      Else
c      SumNU=Zero
c      EndIf
      SumNU=One/Sqrt(Abs(SumNU))
C
      ElseIf(Eig(NU).Ne.Zero) Then
      Write(6,'(X,"Negative ERPA Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      SumNU=Zero
C
      EndIf
C
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=EigVecR((NU-1)*NDimX+I)*SumNU
      EndDo
C
c     enddo NU
      EndDo
C
      Return
      End

*Deck SORT_PINOVECS
      Subroutine SORT_PINOVECS(EigVec,Eig,EigI,NDimX)
      Implicit None
C
      Double Precision :: EigVec(NDimX**2)
      Double Precision :: Eig(NDimX),EigI(NDimX) 
      Integer :: NDimX
      Double Precision :: TmpEig,TmpEigI,TmpVec(NDimX)
      Double Precision,Parameter :: Thresh=1.d-9,EigThresh=1.d-15
      Integer :: I,II,J,IndVI,IndVJ,NDeg,NDegSwp
C
      NDeg = 0
      NDegSwp = 0
      Do I=1,NDimX
      Do J=I+1,NDimX
C     CHECK IF EigThresh NECESSARY?
         If(abs(Eig(I)-Eig(J)).Le.Thresh.and.
     $      abs(Eig(I)).Gt.EigThresh) Then
         II = I + 1
         NDeg = NDeg + 1 
         If(J/=II) Then
           NDegSwp = NDegSwp + 1
           IndVI = (II-1)*NDimX 
           IndVJ = (J-1)*NDimX
C          EIGENVECS
           TmpVec = EigVec(IndVI+1:IndVI+NDimX)
           EigVec(IndVI+1:IndVI+NDimX) = 
     $     EigVec(IndVJ+1:IndVJ+NDimX)
           EigVec(IndVJ+1:IndVJ+NDimX) = TmpVec(1:NDimX)
C          EIGENVALS
           TmpEig = Eig(II) 
           Eig(II) = Eig(J)
           Eig(J) = TmpEig
           TmpEigI = EigI(II) 
           EigI(II) = EigI(J)
           EigI(J) = TmpEigI
         EndIf
         EndIf
      EndDo
      EndDo
C
C      Write(6,*) 'NDeg,Swp',NDeg,NDegSwp     
C
      End Subroutine SORT_PINOVECS


*Deck CREATE_SPACE 
      Subroutine CREATE_SPACE(Eig,Space1,NDimX,DimV1,Max_NDEG)
      Implicit None
C    
      Double Precision :: Eig(NDimX) 
      Integer :: NDimX,DimV1,Max_NDEG
      Integer :: Space1(3,NDimX)
      Double Precision :: Tmp
      Double Precision,Parameter :: Thresh=1.d-8
      Integer :: I
C
      Space1(3,:) = 0
C
      I = 1
      DimV1 = 1
      Space1(1,DimV1) = I
      Tmp = Eig(I)
      Do while(I<NDimX)
      I = I + 1
      If(abs(Eig(I)-Tmp)>Thresh) Then
         Space1(2,DimV1) = I-1
         Space1(3,DimV1) = Space1(2,DimV1)-Space1(1,DimV1)+1
         DimV1=DimV1+1
         Space1(1,DimV1) = I
         Tmp=Eig(I)
      EndIf
      EndDo
      Space1(2,DimV1) = NDimX
      Space1(3,DimV1) = Space1(2,DimV1)-Space1(1,DimV1)+1
C  
      Max_NDEG = maxval(Space1(3,:))
C
      End Subroutine CREATE_SPACE

*Deck ORTHO_DEGVEC
      Subroutine ORTHO_DEGVEC(EigVecR,Space1,DimV1,NDimX,Max_NDEG)
      Implicit None
C
      Integer :: DimV1,NDimX,Max_NDEG
      Integer :: Space1(3,DimV1)
      Integer :: NDEG,IVEC,JVEC
      Integer :: I,II,JJ,INFO 
      Double Precision :: Tmp
      Double Precision :: EigVecR(NDimX**2),Eig
      Double Precision :: Smat(Max_NDEG**2),
     $                    Shlp(max(Max_NDEG**2,3*Max_NDEG)),
     $                    Smh(Max_NDEG**2),Sval(Max_NDEG),
     $                    ModVec(NDimX*Max_NDEG)
      Double Precision,external :: ddot
C
      Do I=1,DimV1
C   
         NDEG=Space1(3,I)
         If(NDEG>1) Then
C
         ! LOOP OVER DEGENERATE VECS
         Do JJ=1,NDEG
         Do II=1,JJ
         IVEC=Space1(1,I)+II-1 
         JVEC=Space1(1,I)+JJ-1 
         Smat(II+(JJ-1)*NDEG)=ddot(NDimX,EigVecR((IVEC-1)*NDimX+1),1,
     $                                   EigVecR((JVEC-1)*NDimX+1),1) 
         EndDo
         EndDo
C
         Call DSYEV('V','U',NDEG,Smat,NDEG,Sval,Shlp,3*NDEG,INFO)
C
C        Shlp = (lambda^-1/2).V^T
C        S^-1/2 = V.Shlp^T
         Do II=1,NDEG
         Shlp((II-1)*NDEG+1:II*NDEG)=Smat((II-1)*NDEG+1:II*NDEG)
     $                               /sqrt(Sval(II))
         EndDo
         Call DGEMM('N','T',NDEG,NDEG,NDEG,1d0,Smat,NDEG,Shlp,NDEG,
     $              0d0,Smh,NDEG)
C    
C        ORTHOGONAL VECTORS: V'=S^-1/2.V     
         Call DGEMM('N','N',NDimX,NDEG,NDEG,1d0,
     $             EigVecR((Space1(1,I)-1)*NDimX+1:),NDimX,Smh,NDEG,
     $             0d0,ModVec,NDimX)
C
         EigVecR((Space1(1,I)-1)*NDimX+1:Space1(2,I)*NDimX)=
     $                               ModVec(1:NDimX*NDEG)
C   
         EndIf
      EndDo
C
      End Subroutine ORTHO_DEGVEC

*Deck ZERO_DEGVEC
      Subroutine ZERO_DEGVEC(EigVecR,Eig,Space1,DimV1,NDimX,Max_NDEG)
      Implicit None
C
      Integer :: DimV1,NDimX,Max_NDEG
      Integer :: Space1(3,DimV1)
      Integer :: NDEG,IVEC,JVEC
      Integer :: I,II,JJ,INFO 
      Double Precision :: Tmp
      Double Precision :: EigVecR(NDimX**2),Eig(NDimX)
      Double Precision,external :: ddot
C
      Do I=1,DimV1
C   
         NDEG=Space1(3,I)
         If(NDEG>1) Then
C
         Print*, 'ZERO!'
         EigVecR((Space1(1,I)-1)*NDimX+1:Space1(2,I)*NDimX) = 0
C  
         EndIf
      EndDo
C
C     ADDITIONAL TEST
C      Write(6,*) 'Compare dimensions:',DimV1,NDimX
C      Do I=1,NDimX
C!       if(Eig(I).Gt.0) then
C        Tmp = norm2(EigVecR((I-1)*NDimX+1:(I-1)*NDimX+NDimX))
C        Write(6,*) Eig(I),Tmp
C!       endif
C      EndDo
C
C
      End Subroutine ZERO_DEGVEC

*Deck ERPAVECYX
      Subroutine ERPAVECYX(EigVecR,EigVecL,Eig,ABPLUS,ABMIN,NDimX)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     set small being a square of small in Deck EneERPA
     $ Four=4.D0, Small=1.D-6)
C
      Integer :: DimV1,Max_NDEG
      Integer :: Space1(3,NDimX)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ HlpAB(NDimX,NDimX),
     $ EigVecR(NDimX*NDimX),EigVecL(NDimX*NDimX),
     $ Eig(NDimX),EigI(NDimX),Work(5*NDimX)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      ABPLUS(I,J)=Half*(ABPLUS(I,J)+ABPLUS(J,I))
      ABPLUS(J,I)=ABPLUS(I,J)
      ABMIN(I,J)=Half*(ABMIN(I,J)+ABMIN(J,I))
      ABMIN(J,I)=ABMIN(I,J)
      EndDo
      EndDo
C
      Call MultpM(HlpAB,ABPLUS,ABMIN,NDimX)
C
      Call DGEEV('N','V',NDimX,HlpAB,NDimX,Eig,EigI,
     $           EigVecL,NDimX,EigVecR,NDimX,Work,5*NDimX,INFO)
C
C     ORTHOGONALISE DEGENERATE VECTORS
      Call SORT_PINOVECS(EigVecR,Eig,EigI,NDimX)
      Call CREATE_SPACE(Eig,Space1,NDimX,DimV1,Max_NDEG)
      Call ORTHO_DEGVEC(EigVecR,Space1,DimV1,NDimX,Max_NDEG)
C     For tests only:
c      Call ZERO_DEGVEC(EigVecR,Eig,Space1,DimV1,NDimEx,Max_NDEG)
C
C     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
C     OMEGA'S 
C
C     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
C
      Do NU=1,NDimX
C
      If(Abs(EigI(NU)).Gt.1.D-12) Then
C
      Write(6,'(X,"Complex ERPA Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=Zero
      EndDo
C
      EndIf
      EndDo
C
      Do NU=1,NDimX
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C   
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NDimX
      Do J=1,NDimX
      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
     $ EigVecR((NU-1)*NDimX+I)*EigVecR((NU-1)*NDimX+J)
      EndDo
      EndDo
C
C      If(SumNu.Lt.Zero) Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
C     $  NU,Eig(NU),SumNU 
C     CHANGE 2
      If(SumNu.Lt.Zero) Then
      Eig(NU)=-Eig(NU)
      Write(6,'(X,"Negative Excit Norm",I4,2E12.4)')
     $  NU,Eig(NU),SumNU
      EndIf
C
C      
      Work(NU)=One
      If(SumNU.Lt.Zero) Work(NU)=-One
      SumNU=One/Sqrt(Abs(SumNU))
C
      ElseIf(Eig(NU).Ne.Zero) Then
      Write(6,'(X,"Negative Omega^2 ERPA Eigenvalue",I4,2E12.4)')
     $ NU,Eig(NU),EigI(NU)
      Eig(NU)=Zero
      SumNU=Zero
C
      EndIf
C
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=EigVecR((NU-1)*NDimX+I)*SumNU
      EndDo
C
c     enddo NU
      EndDo
C
C     COMPUTE EigX
C
      Do NU=1,NDimX
C
C 21.07.2020 KP: avoid dividing by 0 if Eig(NU)=0
      If(Eig(NU).Ne.Zero) Then
C    
      Do I=1,NDimX
      EigVecL((NU-1)*NDimX+I)=Zero
      Do J=1,NDimX
C      EigVecL((NU-1)*NDimX+I)=EigVecL((NU-1)*NDimX+I)
C     $ +Work(NU)*One/Eig(NU)*ABMIN(I,J)*EigVecR((NU-1)*NDimX+J)
C     CHANGE 3
      EigVecL((NU-1)*NDimX+I)=EigVecL((NU-1)*NDimX+I)
     $ +One/Eig(NU)*ABMIN(I,J)*EigVecR((NU-1)*NDimX+J)
C
      EndDo
      EndDo
C
C 21.07.2020 KP: put to zero EigVec's corresponding to Eig(NU)=0
      Else
C
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=Zero
      EigVecL((NU-1)*NDimX+I)=Zero
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
C     CHECK XY NORM AND OMEGA
C
      Do NU=1,NDimX
C
      SumNU=0d0
      Do I=1,NDimX
      SumNU=SumNU+
     $ EigVecR((NU-1)*NDimX+I)*EigVecL((NU-1)*NDimX+I)
      EndDo
C
      If(SumNU.Lt.Zero)Write(6,'(X,"Problems with XY Norm!",I4,2E12.4)')
     $   NU,SumNU,Eig(NU)
C
      If(Eig(NU).Lt.Zero)
     $ Write(6,'(X,"Double Check Negative Excit",I4,2E12.4)')
     $ NU,Eig(NU),SumNU
C
      EndDo
C
      Return
      End

*Deck ERPASYMMXY
      Subroutine ERPASYMMXY(EigY,EigX,Eig,ABPLUS,ABMIN,Occ,IndN,
     $ NDimX,NBasis)
C
C     RETURNS X,Y VECTORS (PROPERLY NORMALIZED), WHICH ARE SOLUTIONS
C     OF THE ORIGINAL ERPA PROBLEM:
C     AX + BY = Om N X
C     BX + AY =-Om N Y 
C
C     THIS IS ACHIEVED BY CALLING ERPASYMM0, SOLVING A SYMMETRIZED PROBLEM     
C     A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
C
C     ON EXIT ABPLUS CONTAINS ABPLUS^(1/2) !!!
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ EigY(NDimX*NDimX),EigX(NDimX*NDimX),Eig(NDimX),
     $ Occ(NBasis),IndN(2,NDimX)
C
      Call ERPASYMM0(EigY,EigX,Eig,ABPLUS,ABMIN,NDimX)
C
      Do K=1,NDimX
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
C
      X=EigX((K-1)*NDimX+I)/(CICoef(IP)+CICoef(IQ))
      Y=EigY((K-1)*NDimX+I)/(CICoef(IP)-CICoef(IQ))
C
      EigX((K-1)*NDimX+I)=Half*(X-Y)
      EigY((K-1)*NDimX+I)=Half*(X+Y)
C
      EndDo
C
      EndDo
C
      Return
      End

*Deck ERPAVECTRANS
      Subroutine ERPAVECTRANS(EigY,EigX,Eig,ABPLUS,ABMIN,Occ,IndN,
     $ NDimX,NBasis)
C
C     RETURNS X,Y VECTORS (PROPERLY NORMALIZED), WHICH ARE SOLUTIONS
C     OF THE ORIGINAL ERPA PROBLEM:
C     AX + BY = Om N X
C     BX + AY =-Om N Y 
C
C     BY TRANSFORMING TILDED X,Y OBTAINED IN ERPAVECYX WHERE A
C     NONSYMMETRIC 
C     PROBLEM IS SOLVED
C
C     A+ A- Y_TILDA = om^2 Y_TILDA
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ EigY(NDimX*NDimX),EigX(NDimX*NDimX),Eig(NDimX),
     $ Occ(NBasis),IndN(2,NDimX)
C
      Call ERPAVECYX(EigY,EigX,Eig,ABPLUS,ABMIN,NDimX)
C
      Do K=1,NDimX
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
C
      X=EigX((K-1)*NDimX+I)/(CICoef(IP)+CICoef(IQ))
      Y=EigY((K-1)*NDimX+I)/(CICoef(IP)-CICoef(IQ))
C
      EigX((K-1)*NDimX+I)=Half*(X-Y)
      EigY((K-1)*NDimX+I)=Half*(X+Y)
C
      EndDo
C
      EndDo
C
      Return
      End

*Deck ERPASYMM0
      Subroutine ERPASYMM0(EigY,EigX,Eig,APLSQRT,ABMIN,NDimX)
C
C     ALMOST THE SAME AS ERPASYMM1 BUT BOTH Y AND X ARE RETURNED
C
C     A SYMMETRIZED PROBLEM A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
C
C     ABPLUS IS CHANGED AND TURNS INTO ABPLUS^(1/2) 
C     THIS ALLOWS ONE TO GET RID OF ONE LOCAL BIG ARRAY (COMPARING WITH ERPASYMM)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     set small being a square of small in Deck EneERPA
     $ Four=4.D0, Small=1.D-6)
C
      Include 'commons.inc'
C
C     APLSQRT originally includes ABPLUS 
C
      Dimension
     $ APLSQRT(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ EigY(NDimX*NDimX),EigX(NDimX*NDimX),Eig(NDimX)
C
C     LOCAL ARRAYS
C
      Dimension HlpAB(NDimX,NDimX),Work(5*NDimX)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      APLSQRT(I,J)=Half*(APLSQRT(I,J)+APLSQRT(J,I))
      APLSQRT(J,I)=APLSQRT(I,J)
      ABMIN(I,J)=Half*(ABMIN(I,J)+ABMIN(J,I))
      ABMIN(J,I)=ABMIN(I,J)
      EndDo
      EndDo
C
C     FIND A+^(1/2)
C
      NoNeg=0
      Call CpyM(HlpAB,APLSQRT,NDimX)
      Call Diag8(HlpAB,NDimX,NDimX,Eig,Work)
C
      Do I=1,NDimX
C
      If(Eig(I).Lt.Zero) NoNeg=NoNeg+1
C
      Do J=I,NDimX
C
      APLSQRT(I,J)=Zero
C 
      Do K=1,NDimX
C
      SQRTEig=SQRT(Abs(Eig(K)))
C
      APLSQRT(I,J)=APLSQRT(I,J)+HlpAB(K,I)*SQRTEig*HlpAB(K,J)
C
      EndDo
C
      APLSQRT(J,I)=APLSQRT(I,J)
C
      EndDo
      EndDo
C
      If(NoNeg.Ne.0) Then
      Write(6,*)"The ERPA A+ matrix is not nonnegative definite"
      Write(6,*)"The number of negative eigenvalues of A+ is",NoNeg
      EndIf
C
C      Call MultpM(HlpAB,ABMIN,APLSQRT,NDimX)
C      Call MultpM(EigY,APLSQRT,HlpAB,NDimX)
C
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,
     $           APLSQRT,NDimX,0d0,HlpAB,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,APLSQRT,NDimX,
     $           HlpAB,NDimX,0d0,EigY,NDimX)
C
      Call Diag8(EigY,NDimX,NDimX,Eig,Work)
C
C     COMPUTE Y's
C
      Do I=1,NDimX
      Do J=I,NDimX
      HlpAB(I,J)=EigY(NDimX*(I-1)+J)
      HlpAB(J,I)=EigY(NDimX*(J-1)+I)
      EndDo
      EndDo

C      Call MultpM(EigY,APLSQRT,HlpAB,NDimX)
C
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,APLSQRT,NDimX,
     $           HlpAB,NDimX,0d0,EigY,NDimX)
C
C     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
C     OMEGA'S 
C
C     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
C
      Do NU=1,NDimX
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C    
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NDimX
      Do J=1,NDimX
      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
     $ EigY((NU-1)*NDimX+I)*EigY((NU-1)*NDimX+J)
      EndDo
      EndDo
C
      If(SumNU.Gt.Zero) Then
      SumNU=One/Sqrt(SumNU)
      Else
      SumNU=Zero
      EndIf
C
      Do I=1,NDimX
      EigY((NU-1)*NDimX+I)=EigY((NU-1)*NDimX+I)*SumNU
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
C     COMPUTE EigX
C
      Do NU=1,NDimX
C
      If(Eig(NU).Gt.Small) Then
C    
      Do I=1,NDimX
      EigX((NU-1)*NDimX+I)=Zero
      Do J=1,NDimX
      EigX((NU-1)*NDimX+I)=EigX((NU-1)*NDimX+I)
     $ +One/Eig(NU)*ABMIN(I,J)*EigY((NU-1)*NDimX+J)
      EndDo
      EndDo
C
      Else
C
      Do I=1,NDimX
      EigY((NU-1)*NDimX+I)=Zero
      EigX((NU-1)*NDimX+I)=Zero
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
      Return
      End

*Deck ERPASYMM1
      Subroutine ERPASYMM1(EigVecR,Eig,APLSQRT,ABMIN,NBasis,NDimX)
      use omp_lib
C
C     A SYMMETRIZED PROBLEM A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
C
C     ABPLUS IS CHANGED AND TURNS INTO ABPLUS^(1/2) 
C     THIS ALLOWS ONE TO GET RID OF ONE LOCAL BIG ARRAY (COMPARING WITH ERPASYMM)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     set small being a square of small in Deck EneERPA
     $ Four=4.D0, Small=1.D-6)
C
      Include 'commons.inc'
C
C     APLSQRT originally includes ABPLUS 
C
      Dimension
     $ APLSQRT(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX)
C
C     LOCAL ARRAYS
C
      Dimension HlpAB(NDimX,NDimX),Work(5*NDimX)
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      APLSQRT(I,J)=Half*(APLSQRT(I,J)+APLSQRT(J,I))
      APLSQRT(J,I)=APLSQRT(I,J)
      ABMIN(I,J)=Half*(ABMIN(I,J)+ABMIN(J,I))
      ABMIN(J,I)=ABMIN(I,J)
      EndDo
      EndDo
C
C     FIND A+^(1/2)
C

      NoNeg=0
      Call CpyM(HlpAB,APLSQRT,NDimX)
      Call Diag8(HlpAB,NDimX,NDimX,Eig,Work)    
C
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(APLSQRT,HlpAB,Eig,NDimX)
      Do I=1,NDimX
C
      If(Eig(I).Lt.Zero) NoNeg=NoNeg+1
C
      Do J=I,NDimX
C
      APLSQRT(I,J)=Zero
C 
      Do K=1,NDimX
C
      SQRTEig=SQRT(Abs(Eig(K)))
C
      APLSQRT(I,J)=APLSQRT(I,J)+HlpAB(K,I)*SQRTEig*HlpAB(K,J)
C
      EndDo
C
      APLSQRT(J,I)=APLSQRT(I,J)
C
      EndDo
      EndDo
!$OMP END PARALLEL DO 
C
      If(NoNeg.Ne.0) Then
      Write(6,*)"The ERPA A+ matrix is not nonnegative definite"
      Write(6,*)"The number of negative eigenvalues of A+ is",NoNeg
      Else
      Write(6,*)"ERPA A+ matrix is positive definite"
      EndIf
C 
      Call MultpM(HlpAB,ABMIN,APLSQRT,NDimX)
      Call MultpM(EigVecR,APLSQRT,HlpAB,NDimX)
C
      Write(6,*)"Begin diagonalization of ERPA"
      Call Diag8(EigVecR,NDimX,NDimX,Eig,Work)
      Write(6,*)"Done with diagonalization of ERPA"
C
C     COMPUTE Y's
C
      Do I=1,NDimX
      Do J=I,NDimX
      HlpAB(I,J)=EigVecR(NDimX*(I-1)+J)
      HlpAB(J,I)=EigVecR(NDimX*(J-1)+I)
      EndDo
      EndDo
      Call MultpM(EigVecR,APLSQRT,HlpAB,NDimX)
      Write(6,*)"Y vectors computed"
C
C     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
C     OMEGA'S 
C
C     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
C
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(EigVecR,NDimX,Eig,ABMIN)
      Do NU=1,NDimX
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C    
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NDimX
      Do J=1,NDimX
      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
     $ EigVecR((NU-1)*NDimX+I)*EigVecR((NU-1)*NDimX+J)
      EndDo
      EndDo
C
      If(SumNU.Gt.Zero) Then
      SumNU=One/Sqrt(SumNU)
      Else
      SumNU=Zero
      EndIf
C
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=EigVecR((NU-1)*NDimX+I)*SumNU
      EndDo
C
      EndIf
c     enddo NU
      EndDo
!$OMP END PARALLEL DO 
C
      Write(6,*)"ERPA equation solved"
C
      Return
      End

*Deck ERPASYMM
      Subroutine ERPASYMM(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
C
C     A SYMMETRIZED PROBLEM A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
C
C     MATRICES ABPLUS,ABMIN ARE NOT CHANGED! 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
c     set small being a square of small in Deck EneERPA
     $ Four=4.D0, Small=1.D-6)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX)
C
C     LOCAL ARRAYS
C
      Dimension APLSQRT(NDimX,NDimX),HlpAB(NDimX,NDimX),Work(5*NDimX)
C
c      Write(6,'(/,X,"*** SYMMETRIC ERPA PROBLEM IS SOLVED *** ")')
C
C     SYMMETRIZE A+,A-
C
      Do I=1,NDimX
      Do J=I+1,NDimX
      ABPLUS(I,J)=Half*(ABPLUS(I,J)+ABPLUS(J,I))
      ABPLUS(J,I)=ABPLUS(I,J)
      ABMIN(I,J)=Half*(ABMIN(I,J)+ABMIN(J,I))
      ABMIN(J,I)=ABMIN(I,J)
      EndDo
      EndDo
C 
      Work = 0
C
C     FIND A+^(1/2)
C
      NoNeg=0
      Call CpyM(HlpAB,ABPLUS,NDimX)
      Call Diag8(HlpAB,NDimX,NDimX,Eig,Work) 
C
      Do I=1,NDimX
C
      If(Eig(I).Lt.Zero) NoNeg=NoNeg+1 
C
      Do J=I,NDimX
C
      APLSQRT(I,J)=Zero
C 
      Do K=1,NDimX
C
      If(Eig(K).Ge.Zero) Then
      SQRTEig=SQRT(Eig(K))
      Else
      SQRTEig=SQRT(Abs(Eig(K)))
      EndIf
C
      APLSQRT(I,J)=APLSQRT(I,J)+HlpAB(K,I)*SQRTEig*HlpAB(K,J)
C
      EndDo
C
      APLSQRT(J,I)=APLSQRT(I,J)
C
      EndDo
      EndDo 
C
      If(NoNeg.Ne.0) Then
      Write(6,*)"The ERPA A+ matrix is not nonnegative definite"
      Write(6,*)"The number of negative eigenvalues of A+ is",NoNeg
      EndIf
C
      Call MultpM(HlpAB,ABMIN,APLSQRT,NDimX)
      Call MultpM(EigVecR,APLSQRT,HlpAB,NDimX)
C
C      TEST FOR SAPT
C      Do I=1,NDimX
C      write(6,*) I, Eig(I), EigVecR(NDimX*(I-1)+I)
C      EndDo
C 
      Call Diag8(EigVecR,NDimX,NDimX,Eig,Work)
C
C      Do I=1,NDimX
C      write(6,*) I, Eig(I), EigVecR(NDimX*(I-1)+I)
C      EndDo 
C
C     COMPUTE Y's
C
      Do I=1,NDimX
      Do J=I,NDimX
      HlpAB(I,J)=EigVecR(NDimX*(I-1)+J) 
      HlpAB(J,I)=EigVecR(NDimX*(J-1)+I)
C      Write(6,*) HlpAB(I,J), HlpAB(J,I)
      EndDo
      EndDo
      Call MultpM(EigVecR,APLSQRT,HlpAB,NDimX) 
C
C     IMPOSE THE NORMALIZATION 2 Y*X = 1 ON THE EIGENVECTORS CORRESPONDING TO POSITIVE
C     OMEGA'S 
C
C     SINCE X = Om^-1 ABMIN.Y THEN THE NORMALIZATION READS 2 Om^-1 Y^T AMIN Y = 1
C
      Do NU=1,NDimX
      SumNU=Zero
C
      If(Eig(NU).Gt.Small) Then
C    
      Eig(NU)=SQRT(Eig(NU))
C
      Do I=1,NDimX
      Do J=1,NDimX
      SumNU=SumNU+Two/Eig(NU)*ABMIN(I,J)*
     $ EigVecR((NU-1)*NDimX+I)*EigVecR((NU-1)*NDimX+J)
      EndDo
      EndDo
C
      If(SumNU.Gt.Zero) Then
      SumNU=One/Sqrt(SumNU)
      Else
      SumNU=Zero
      EndIf
C
      Do I=1,NDimX
      EigVecR((NU-1)*NDimX+I)=EigVecR((NU-1)*NDimX+I)*SumNU
      EndDo
C
      EndIf
c     enddo NU
      EndDo
C
C     TEST FOR SAPT
C      Do I=1,NDimX
C      write(6,*) I, EigVecR(NDimX*(I-1)+I)
C      EndDo 

      Return
      End

* Deck OneTwoBody
      Subroutine OneTwoBody(ETot,ENuc,ECorrTot,EGOne,EGOneTwo,
     $ EigVecR,Eig,
     $ ABPLUS,ABMIN,
     $ Occ,TwoNO,IndAux,NBasis,NInte1,NInte2,NDim,NGem,NGOcc,IFlag)
C
C     A ROUTINE FOR COMPUTING ONE- AND TWO-BODY INTERGEMINAL CORRELATION ENERGIES 
C     ONE-BODY IS ONLY DEFINED FOR GVB ! 
C
C     IFlag = 0  : A+,A- MATRICES NOT TRUNCATED
C           = 1  : A+,A- MATRICES ARE TRUNCATED TO INCLUDE ONLY INDICES
C                        OF ORBITALS BELONGING TO GEMINALS NEEDED FOR 
C                        A GIVEN ONE- TWO- THREE- OR FOUR-BODY INTERACTION
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension
     $ EGOne(NGem),EGOneTwo(NGem,NGem),
     $ EigVecR(NDim,NDim),Eig(NDim),
     $ IndAux(NBasis),Occ(NBasis),TwoNO(NInte2),
     $ ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
C
C     LOCAL ARRAYS
C
      Dimension IndN(2,NDim),IndX(NDim),IAux(NBasis),
     $          ABPRed(NDim,NDim),ABMRed(NDim,NDim)
C
C
      ECorrTot=Zero
C
C
      If(IGVB.Eq.1) Then
C
      Write(6,'(2/,X,"*******************************")')
      Write(6,'(   X,"GVB ONE-BODY CORRELATION ENERGY")')
      Write(6,'(   X,"*******************************")')
C
C     ONE-BODY CONTRIBUTIONS
C
      IStart=1
      Do IGG=IStart,NGem-1
C
      Write(6,'(/,X,"COMPUTING CORR FOR THE FRAGMENT: ",I4)')IGG
C      
      Eig(1:NDim) = 0
      EigVecR = 0
C
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
C
C     Check if the fragment is core-type (only core orbitals are included)
C
      ICore=1
      Do I=1,NBasis
      If(IGem(I).Eq.IGG.And.Occ(I).Ne.One) ICore=0
      EndDo
C
      If(ICore.Eq.1) Then
C
      Write(6,'(X,"Only Core Orbitals in the Fragment")')
C     If IFlCore=0 do not compute correlation 
      If(IFlCore.Eq.0) Then
      Write(6,'(X,"IFlCore=0: Do not correlate")')
      NDimX=0
      GoTo 555
      EndIf
C
c     If(ICore.Eq.1)
      Else  
C
      Do I=1,NBasis
C
      If(IndAux(I).Eq.2.And.IGem(I).Ne.IGG.And.Occ(I).Ne.Zero) Then
C     check if fragments are not connected
      If(IConnect(IGG,IGem(I)).Eq.0) Then
      Write(6,'(X," Orbital Included (I, Frag, Occ): ",2I4,E14.4)')
     $ I, IGem(I), Occ(I) 
      IGem(I)=NGem
      EndIf
      EndIf
C
c include active orbitals
      If(IndAux(I).Eq.1.And.Occ(I).Lt.Half.And.IGem(I).Ne.IGG) Then
C     check if fragments are not connected
      If(IConnect(IGG,IGem(I)).Eq.0) Then
       Write(6,'(X," Active Orbital Included (I, Frag, Occ): ",
     $ 2I4,E14.4)') I, IGem(I), Occ(I)
      IGem(I)=NGem
      EndIf 
      EndIf
C
c     Do I=1
      EndDo
C
c     If(ICore.Eq.1)
      EndIf
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate active degenerate orbitals if from different geminals
      If((IAuxGem(I).Ne.IAuxGem(J)).And.(IndAux(I).Eq.1).And.
     $ (IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct) ) Then
C
      Write(*,*)"Discarding nearly degenerate pair",I,J
C
      Else
C
      If(IFlag.Eq.0) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
C     else of IFlag.Eq.0
      Else
C
      If((IGem(I).Eq.NGem.And.IGem(J).Eq.IGG).Or.
     $   (IGem(J).Eq.NGem.And.IGem(I).Eq.IGG).Or.
     $ (IGem(I).Eq.IGG.And.IGem(J).Eq.IGG)) Then
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
C
C     endif IFlag.Eq.0
      EndIf
C
c     do not correlate active degenerate orbitals if from different geminals 
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
      Write(*,*)'NDimX',NDimX
c      Write(6,'(X,"Accepted pairs read:")')
c      Do I=1,Ind
c      Ind1=IndN(1,I)
c      Ind2=IndN(2,I)
c      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
c      EndDo
C
      If(NDimX.Eq.0) Write(6,'(" No allowed pairs for the fragment")')
C 
  555 Continue
C
      If(NDimX.Ne.0) Then
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPRed(I,J) = ABPLUS(IndX(I),IndX(J))
      ABMRed(I,J) =  ABMIN(IndX(I),IndX(J))
      EndDo
      EndDo
c      Print*, 'ABPRed',norm2(ABPRed)
c      Print*, 'ABMRed',norm2(ABMRed)
C
C     COMPUTE THE ONE-BODY CORRELATION
C
C     Do it only once if IFlag.Eq.0
C
      If(IFlSnd.Eq.0) Then
C
      If(IFlag.Eq.1.Or.(IFlag.Eq.0.And.IGG.Eq.IStart)) Then
      Call ERPASYMM(EigVecR(1:NDimX,1:NDimX),Eig,
     $  ABPRed(1:NDimX,1:NDimX),ABMRed(1:NDimX,1:NDimX),NBasis,NDimX)
      EndIf
C
      Else
C
      Do J=1,NDimX
      Do I=1,NDimX
      EigVecR(I,J)=ABPRed(I,J)
      EndDo
      EndDo
C
      EndIf
C
c     print*, 'EigVecR',norm2(EigVecR)
c     print*, 'Eig    ',norm2(Eig)
C
      Call EInterG(ECorr,EigVecR(1:NDimX,1:NDimX),Eig,TwoNO,Occ,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem,2,IGG,IGG,NGem,NGem)
      EGOne(IGG)=ECorr
C
c     If(NDimX.Ne.0) 
      Else
C
      EGOne(IGG)=Zero 
C
c     If(NDimX.Ne.0) 
      EndIf
C
      Write(6,'(X,"1-body correlation: ",F15.8)') EGOne(IGG)
C
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
c     enddo IGG 
      EndDo
C
C     SUM AND PRINT ONE-BODY CONTRIBUTIONS
C
      EOneTot=Zero
      IStart=1
C
      Write(6,'(/,2X,"*** 1-body correlation for each fragment")')
      Do IGG=IStart,NGem-1
      Write(6,'(X,I4,F15.8)') IGG,EGOne(IGG)
      EOneTot=EOneTot+EGOne(IGG)
      EndDo
C
      Write(6,'(/,2X,"*** Total 1-body correlation",2X,F15.8)')EOneTot
      ECorrTot=ECorrTot+EOneTot
C
      Write(6,'(2X,"EGVB + ENuc + 1-body",6X,F15.8)')ETot+ENuc+EOneTot
C
C     TWO-BODY CONTRIBTUIONS
C
      Write(6,'(2/,X,"*******************************")')
      Write(6,'(   X,"GVB TWO-BODY CORRELATION ENERGY")')
      Write(6,'(   X,"*******************************")')
C
      IStart=1
      Do IGG1=IStart,NGem-1
      Do IGG2=IStart,IGG1-1
C
      Write(6,'(/,X,"COMPUTING CORR FOR FRAGMENTS: ",2I4)')IGG1,IGG2
C
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
C
C     Check if one of the fragments is core-type
C
      ICore=0
      Sum=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.IGG1.And.Occ(I).Ne.One) Sum=Sum+Occ(I)
      EndDo
      If(Sum.Eq.Zero) ICore=1
      Sum=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.IGG2.And.Occ(I).Ne.One) Sum=Sum+Occ(I)
      EndDo
      If(Sum.Eq.Zero) ICore=1
C
      If(ICore.Eq.1) Then
C
      Write(6,'(X,"Only Core Orbitals in One of the Fragments")')
C     If IFlCore=0 do not compute correlation 
      If(IFlCore.Eq.0) Then
      Write(6,'(X,"IFlCore=0: Do not correlate")')
      NDimX=0
      GoTo 777
      EndIf 
C
      Else
C
C     only for connected fragments
      If(IConnect(IGG1,IGG2).Eq.1) Then
C
      Write(6,'(X,"Fragments connected")')
C
      Do I=1,NBasis
C
      If(IndAux(I).Eq.2.And.IGem(I).Ne.IGG1.And.IGem(I).Ne.IGG2.
     $ And.Occ(I).Ne.Zero) Then
C
C     check if fragments are not connected
      If(IConnect(IGG1,IGem(I)).Eq.0.And.IConnect(IGG2,IGem(I)).Eq.0) 
     $ Then
      Write(6,'(X," Orbital Included (I, Frag, Occ): ",2I4,E14.4)')
     $ I, IGem(I), Occ(I)
      IGem(I)=NGem
      EndIf
C
      EndIf
C
c include active orbitals
      If(IndAux(I).Eq.1.And.Occ(I).Lt.Half.And.IGem(I).Ne.IGG1
     $ .And.IGem(I).Ne.IGG2) Then
C     check if fragments are not connected
      If(IConnect(IGG1,IGem(I)).Eq.0.And.IConnect(IGG2,IGem(I)).Eq.0) 
     $ Then
       Write(6,'(X," Active Orbital Included (I, Frag, Occ): ",
     $ 2I4,E14.4)') I, IGem(I), Occ(I)
      IGem(I)=NGem
      EndIf
C
      EndIf
C
C     Do I=1
      EndDo
C
c     If(IConnect(IGG1,IGG2).Eq.1)
      Else
      Write(6,'(X,"Fragments disconnected")')
C
C     If(IConnect(IGG1,IGG2).Eq.0) Then
      EndIf
C
c     If(ICore.Eq.1)
      EndIf
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate degenerate active degenerate orbitals if from different geminals
      If((IAuxGem(I).Ne.IAuxGem(J)).And.(IndAux(I).Eq.1).And.
     $ (IndAux(J).Eq.1).And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct) )
     $ Then
C
      Write(*,*)"Discarding nearly degenerate pair",I,J
C
      Else
C
      If(IFlag.Eq.0) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
C     else IFlag.Eq.0
      Else
C
      If((IGem(I).Eq.NGem.And.(IGem(J).Eq.IGG1.Or.IGem(J).Eq.IGG2)).Or.
     $ ((IGem(I).Eq.IGG1.Or.IGem(I).Eq.IGG2).And.
     $ (IGem(J).Eq.IGG1.Or.IGem(J).Eq.IGG2))) Then
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
C
C     endif IFlag.Eq.0
      EndIf
C
c     do not correlate degenerate active degenerate orbitals if from different geminals 
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
c      Write(*,*)'NDimX',NDimX
c      Write(6,'(X,"Accepted pairs read:")')
c      Do I=1,Ind
c      Ind1=IndN(1,I)
c      Ind2=IndN(2,I)
c      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
c      EndDo
C
C
      If(NDimX.Eq.0) Write(6,'(" No allowed pairs for the fragments.")')
C
  777 Continue
C
      If(NDimX.Ne.0) Then
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPRed(I,J) = ABPLUS(IndX(I),IndX(J))
      ABMRed(I,J) =  ABMIN(IndX(I),IndX(J))
      EndDo
      EndDo
C
      Print*, 'ABPRed',norm2(ABPRed)
      Print*, 'ABMRed',norm2(ABMRed)
C
C     COMPUTE THE TWO-BODY CORRELATION
C
C     Call ERPAVEC only if IFlag=1, for IFlag=0 EigVecR is known from calculations of one-body energy
C
      If(IFlSnd.Eq.0) Then
C
      If(IFlag.Eq.1) Then
      Call ERPASYMM(EigVecR(1:NDimX,1:NDimX),Eig,
     $  ABPRed(1:NDimX,1:NDimX),ABMRed(1:NDimX,1:NDimX),NBasis,NDimX)
      EndIf
CC
      Else
C
      Do J=1,NDimX
      Do I=1,NDimX
      EigVecR(I,J)=ABPRed(I,J)
      EndDo
      EndDo
C
      EndIf
C
      print*, 'EigVecR',norm2(EigVecR)
      print*, 'Eig    ',norm2(Eig)
C
      Call EInterG(ECorr,EigVecR(1:NDimX,1:NDimX),Eig,TwoNO,Occ,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem,3,IGG1,IGG2,NGem,NGem)
      EGOneTwo(IGG1,IGG2)=ECorr
C
c     If(NDimX.Ne.0) 
      Else
C
      EGOneTwo(IGG1,IGG2)=Zero
C
c     If(NDimX.Ne.0) 
      EndIf
C
      Write(6,'(X,"2-body correlation: ",F15.8)') EGOneTwo(IGG1,IGG2)
C
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
c     enddo IGG1,IGG2
      EndDo
      EndDo
C
C     SUM AND PRINT TWO-BODY CONTRIBUTIONS
C
      ETwoTot=Zero
      IStart=1
C
      Write
     $ (6,'(/,2X,"*** 2-body correlation for each pair of fragments")')
      Do IGG1=IStart,NGem-1
      Do IGG2=IStart,IGG1-1
      Write(6,'(X,2I4,F15.8)') IGG1,IGG2,EGOneTwo(IGG1,IGG2)
      ETwoTot=ETwoTot+EGOneTwo(IGG1,IGG2)
      EndDo
      EndDo
      Write(6,'(/,2X,"*** Total 2-body correlation",3X,F15.8)')ETwoTot
      ECorrTot=ECorrTot+ETwoTot
cC
c      Write(6,'(/,2X,"1-body + 2-body correlation",F15.8)')
c     $ EOneTot+ETwoTot
cC
      Write(6,'(/,2X,"EGVB + ENuc + 1,2-body",5X,F15.8)')
     $ ETot+ENuc+EOneTot+ETwoTot
C
C    TESTING...
      Return
C
C     THREE-BODY CONTRIBUTIONS
C
      Write(6,'(2/,X,"*********************************")')
      Write(6,'(   X,"GVB THREE-BODY CORRELATION ENERGY")')
      Write(6,'(   X,"*********************************")')
      Write(6,'(/,X,
     $ "THREE-BODY IS COMPUTED ONLY FOR CONNECTED  FRAGMENTS")')
C
      EThreeTot=Zero
C
      IStart=1
      Do IGG1=IStart,NGem-1
      Do IGG2=IStart,IGG1-1
      Do IGG3=IStart,IGG2-1
C
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
C
C     CHECK IF FRAGMENTS ARE DISCONNECTED
C
      If(IConnect(IGG1,IGG2)*IConnect(IGG1,IGG3)
     $ *IConnect(IGG2,IGG3).Ne.1) Then 
      NDimX=0
      GoTo 888
      EndIf
C
      Write(6,'(/,X,"COMPUTING CORR FOR FRAGMENTS: ",3I4)')
     $ IGG1,IGG2,IGG3
      Write(6,'(X,"Fragments connected")') 
C
C     Check if one of the fragments is core-type
C
      ICore=0
      Sum=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.IGG1.And.Occ(I).Ne.One) Sum=Sum+Occ(I)
      EndDo
      If(Sum.Eq.Zero) ICore=1
      Sum=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.IGG2.And.Occ(I).Ne.One) Sum=Sum+Occ(I)
      EndDo
      If(Sum.Eq.Zero) ICore=1
      Sum=Zero
      Do I=1,NBasis
      If(IGem(I).Eq.IGG3.And.Occ(I).Ne.One) Sum=Sum+Occ(I)
      EndDo
      If(Sum.Eq.Zero) ICore=1
C
      If(ICore.Eq.1) Then
C
      Write(6,'(X,"Only Core Orbitals in one of the Fragments")')
C
C     If IFlCore=0 do not compute correlation 
      If(IFlCore.Eq.0) Then
      Write(6,'(X,"IFlCore=0: Do not correlate")')
      NDimX=0
      GoTo 888 
      EndIf
C
      Else
C
      Do I=1,NBasis
C
      If(IndAux(I).Eq.2.And.IGem(I).Ne.IGG1.And.IGem(I).Ne.IGG2.
     $ And.IGem(I).Ne.IGG3.And.Occ(I).Ne.Zero) Then
C     check if fragments are not connected
      If(IConnect(IGG1,IGem(I)).Eq.0.And.IConnect(IGG2,IGem(I)).Eq.0.
     $ And.IConnect(IGG3,IGem(I)).Eq.0)
     $ Then
      Write(6,'(X," Orbital Included (I, Frag, Occ): ",2I4,E14.4)')
     $ I, IGem(I), Occ(I)
      IGem(I)=NGem
      EndIf
      EndIf
C
c include active orbitals
      If(IndAux(I).Eq.1.And.Occ(I).Lt.Half.And.IGem(I).Ne.IGG1
     $ .And.IGem(I).Ne.IGG2.And.IGem(I).Ne.IGG3) Then
C     check if fragments are not connected
      If(IConnect(IGG1,IGem(I)).Eq.0.And.IConnect(IGG2,IGem(I)).Eq.0.
     $ And.IConnect(IGG3,IGem(I)).Eq.0) Then
       Write(6,'(X," Active Orbital Included (I, Frag, Occ): ",
     $ 2I4,E14.4)') I, IGem(I), Occ(I)
      IGem(I)=NGem
      EndIf
      EndIf
C
c     Do I=1 
      EndDo
C
C
c     If(ICore.Eq.1)
      EndIf
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
C     do not correlate idegenerate active degenerate orbitals if from different geminals
      If((IAuxGem(I).Ne.IAuxGem(J)).And.(IndAux(I).Eq.1).And.
     $ (IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct) ) Then
C     
      Write(*,*)"Discarding nearly degenerate pair",I,J
C
      Else
C
      If(IFlag.Eq.0) Then

      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
C     else IFlag.Eq.0
      Else
C
      If(
     $ (IGem(I).Eq.NGem.And.(IGem(J).Eq.IGG1.Or.IGem(J).Eq.IGG2.
     $                       Or.IGem(J).Eq.IGG3)).
     $ Or.
     $ ((IGem(I).Eq.IGG1.Or.IGem(I).Eq.IGG2.Or.IGem(I).Eq.IGG3).And.
     $ ( IGem(J).Eq.IGG1.Or.IGem(J).Eq.IGG2.Or.IGem(J).Eq.IGG3))) Then
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
C
C     endif IFlag.Eq.0
      EndIf
C
      EndIf
C
c     do not correlate degenerate active degenerate orbitals if from different geminals 
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
      If(NDimX.Eq.0) Write(6,'(" No allowed pairs for the fragments.")')
C
  888 Continue
C
      If(NDimX.Ne.0) Then
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPRed(I,J) = ABPLUS(IndX(I),IndX(J))
      ABMRed(I,J) =  ABMIN(IndX(I),IndX(J))
      EndDo
      EndDo
C
C     COMPUTE THE THREE-BODY CORRELATION
C
C     Call ERPAVEC only if IFlag=1, for IFlag=0 EigVecR is known from calculations of one-body energy
C
      If(IFlSnd.Eq.0) Then
C
      If(IFlag.Eq.1) Then
      Call ERPASYMM(EigVecR(1:NDimX,1:NDimX),Eig,
     $ ABPRed(1:NDimX,1:NDimX),ABMRed(1:NDimX,1:NDimX),NBasis,NDimX)
      EndIf
C
      Else
C
      Do J=1,NDimX
      Do I=1,NDimX
      EigVecR(I,J)=ABPRed(I,J)
      EndDo
      EndDo
C
      EndIf
C
      Call EInterG(ECorr,EigVecR(1:NDimX,1:NDimX),Eig,TwoNO,Occ,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem,4,IGG1,IGG2,IGG3,NGem)
C
      EThreeTot=EThreeTot+ECorr
C
      Write(6,'(X,"3-body correlation: ",F15.8)') ECorr
C
c     If(NDimX.Ne.0)
      EndIf
C
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,2X,"*** Total 3-body correlation",3X,F15.8)')EThreeTot
      ECorrTot=ECorrTot+EThreeTot
C
      Write(6,'(/,2X,"EGVB + ENuc,  1,2,3-body corr, Total",3X,3F15.8)')
     $ ETot+ENuc,EOneTot+ETwoTot+EThreeTot,
     $ ETot+ENuc+EOneTot+ETwoTot+EThreeTot
C
      Return 
C
C     else IGVB.Eq.1
C
      Else
C
C     APSG
C
C     TWO-BODY CONTRIBTUIONS
C
      Write(6,'(/,X," APSG TWO-BODY CORRELATION ENERGY")')
C
      IStart=1
c      If(IFreeze.Eq.1) IStart=1+NGOcc
      Do IGG1=IStart,NGem
      Do IGG2=IStart,IGG1-1
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C
      If(IFlag.Eq.0) Then

      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C
C     else IFlag.Eq.0
      Else
C
      If ( (IGem(I).Eq.IGG1.Or.IGem(I).Eq.IGG2).And.
     $     (IGem(J).Eq.IGG1.Or.IGem(J).Eq.IGG2) ) Then
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      EndIf
C
C     endif IFlag.Eq.0
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPRed(I,J) = ABPLUS(IndX(I),IndX(J))
      ABMRed(I,J) =  ABMIN(IndX(I),IndX(J))
      EndDo
      EndDo
C
C     COMPUTE THE TWO-BODY CORRELATION
C
      Call ERPASYMM(EigVecR(1:NDimX,1:NDimX),Eig,
     $  ABPRed(1:NDimX,1:NDimX),ABMRed(1:NDimX,1:NDimX),NBasis,NDimX)
C
      Call EInterG(ECorr,EigVecR(1:NDimX,1:NDimX),Eig,TwoNO,Occ,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem,2,IGG1,IGG1,IGG2,IGG2)
      EGOneTwo(IGG1,IGG2)=ECorr
C
c     enddo IGG
      EndDo
      EndDo
C
C     SUM AND PRINT APSG TWO-BODY CONTRIBUTIONS
C
      ETwoTot=Zero
      IStart=1
c      If(IFreeze.Eq.1) IStart=1+NGOcc
C
      Do IGG1=IStart,NGem
      Do IGG2=IStart,IGG1-1
      Write(6,'(X,2I4,F15.8)') IGG1,IGG2,EGOneTwo(IGG1,IGG2)
      ETwoTot=ETwoTot+EGOneTwo(IGG1,IGG2)
      EndDo
      EndDo
      Write(6,'(2X,"Total 2-body correlation",3X,F15.8)')ETwoTot
      ECorrTot=ECorrTot+ETwoTot
C
      Write(6,'(2X,"EAPSG + ENuc + 2-body",6X,F15.8)')
     $ ETot+ENuc+ETwoTot
C
C     endif IGVB.Eq.1
C
      EndIf
C
      Return
      End

*Deck EInterG
      Subroutine EInterG(ECorr,EigVecR,Eig,TwoNO,Occ,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem,IB,IG1,IG2,IG3,IG4)
C
      use abfofo
      use abmat
C
C     IB = 2 - ONLY TWO-BODY CONTRIBUTIONS (ONE-BODY FOR GVB)
C          3 - ONLY THREE-BODY (TWO-BODY FOR GVB)
C          4   ONLY FOUR-BODY
C     INDICES IP,IQ,IR,IS MUST BELONG TO ONE OF THE G1-G4 GEMINALS
C
      Implicit Real*8 (A-H,O-Z)
c
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Parameter(SmallE=1.D-3,BigE=5.D2)
C
C     ONLY EXCITATIONS SMALLER THAN BigE AND GREATER THAN SmallE ARE INCLUDED
C
      Include 'commons.inc'
C
      Dimension EigVecR(NDimX*NDimX),Eig(NDimX),
     $ Occ(NBasis),TwoNO(NInte2),IndN(2,NDimX)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),Skipped(NDimX)
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C
C     COMPUTE INTERGEMINAL CORRELATION, ECorr
C
      ECorr=Zero
C
      If(ITwoEl.Eq.1) Then
C
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IFlP=0
      If(IGem(IP).Eq.IG1.Or.IGem(IP).Eq.IG2.Or.IGem(IP).Eq.IG3.
     $ Or.IGem(IP).Eq.IG4) IFlP=1
      NoVirtP=1
      If(Occ(IP).Eq.Zero) NoVirtP=0
C
      IR=IndN(2,I)
      IFlR=0
      If(IGem(IR).Eq.IG1.Or.IGem(IR).Eq.IG2.Or.IGem(IR).Eq.IG3.
     $ Or.IGem(IR).Eq.IG4) IFlR=1
      NoVirtR=1
      If(Occ(IR).Eq.Zero) NoVirtR=0
C
      Do J=1,NDimX
C
      IQ=IndN(1,J)
      IFlQ=0
      If(IGem(IQ).Eq.IG1.Or.IGem(IQ).Eq.IG2.Or.IGem(IQ).Eq.IG3.
     $ Or.IGem(IQ).Eq.IG4) IFlQ=1
      NoVirtQ=1
      If(Occ(IQ).Eq.Zero) NoVirtQ=0
C
      IS=IndN(2,J)
      IFlS=0
      If(IGem(IS).Eq.IG1.Or.IGem(IS).Eq.IG2.Or.IGem(IS).Eq.IG3.
     $ Or.IGem(IS).Eq.IG4) IFlS=1
      NoVirtS=1
      If(Occ(IS).Eq.Zero) NoVirtS=0
C
      IFlPRQS=IFlP*IFlR*IFlQ*IFlS
      NoVirt=NoVirtP*NoVirtR*NoVirtQ*NoVirtS
C
      Call IBody(IBdy,IGem(IP),IGem(IR),IGem(IQ),IGem(IS))
C
      ICond=0
C
      If((IP.Gt.IR.And.IQ.Gt.IS.And.IBdy.Eq.IB.And.IFlPRQS.Eq.1))
     $   ICond=1

C     FOR GVB ONLY
C
      If(IGVB.Eq.1.And.IB.Gt.2.And.NoVirt.Eq.1.And.
     $   IP.Gt.IR.And.IQ.Gt.IS.And.IBdy.Eq.IB-1.And.IFlPRQS.Eq.1)
     $ ICond=1
C   
C     FOR GVB ONLY: IF IFrag=1,IB=2 (ONE-BODY), ALLOW ALL CASES 
C     EXCEPT WHEN ALL ORBITALS ARE FROM THE SAME GEMINAL
C
      If(IFrag.Eq.1.And.IGVB.Eq.1.And.IB.Eq.2.And.
     $   IP.Gt.IR.And.IQ.Gt.IS.And.IFlPRQS.Eq.1) Then
  
C     IF IFrag=1 THEN IAuxGem STORES ASSIGNMENTS OF ORBS TO GEMINALS 
C     (IGem STORES ASSIGNMENTS TO FRAGMENTS)
      Call IBody(IBdyG,IAuxGem(IP),IAuxGem(IR),IAuxGem(IQ),IAuxGem(IS))
      If(IBdyG.Ne.1) ICond=1
      EndIf
C
      If(ICond.Eq.1) Then
C
      ISkippedEig=0
C
      If(IFlSnd.EQ.0) Then
C
      SumY=Zero
      Do K=1,NDimX
      If(Eig(K).Gt.SmallE.And.Eig(K).Lt.BigE) Then
      SumY=SumY+EigVecR((K-1)*NDimX+I)*EigVecR((K-1)*NDimX+J)
      Else
      ISkippedEig=ISkippedEig+1
      Skipped(ISkippedEig)=Eig(K)  
      EndIf
      EndDo
C
      Aux=(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      If(IR.Eq.IS.And.IP.Eq.IQ) Aux=Aux
     $ -Half*(Occ(IP)*(One-Occ(IS))+Occ(IS)*(One-Occ(IP)))
C
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      Else
C
      Aux=EigVecR((J-1)*NDimX+I) 
C
      If(IR.Eq.IS.And.IP.Eq.IQ) Then
      If(Occ(IP).Gt.Occ(IR)) Then
      Aux1=Occ(IR)*(Occ(IP)-One)
      Else
      Aux1=Occ(IP)*(Occ(IR)-One)
      EndIf
      Aux=Aux+Two*Aux1
      EndIf
C
      ECorr=ECorr+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))

C     If(IFlSnd.Eq.0) 
      EndIf
C
C     endinf of If(ICond.Eq.1)
      EndIf
C
      EndDo
      EndDo
C
      If(ISkippedEig.Gt.0) Then
      Write(6,'(/,1X,"The number of discarded eigenvalues is",I4)')
     $  ISkippedEig
      Do II=1,ISkippedEig
      Write(6,*)'Skipped',II,Skipped(II)
      EndDo
      EndIf

      ElseIf(ITwoEl.Eq.2) Then
C
      Call EERPA_FFFF(ECorr,EigVecR,Eig,Occ,
     $                C,IGem,IAuxGem,
     $                IG1,IG2,IG3,IG4,IB,
     $                IndN,NDimX,NBasis,
     $                'TWOMO',IFrag)
C
      ElseIf(ITwoEl.Eq.3) Then
C
      Call EERPA_FOFO(ECorr,EigVecR,Eig,Occ,
     $                C,IGem,IAuxGem,
     $                IG1,IG2,IG3,IG4,IB,
     $                IndN,NELE+NActive,
     $                NDimX,NBasis,'FOFO',IFrag)
C
      EndIf
C
      Return
      End 

*Deck IBody
      Subroutine IBody(IBdy,IGP,IGR,IGQ,IGS)
C
      Implicit Real*8 (A-H,O-Z)
C
C     DETERMINES IF IP,IR,IQ,IS BELONG TO:
C     ONE GEMINAL: IBdy=1
C     TWO DIFFERENT GEMINALS: IBdy=2
C     THREE DIFFERENT GEMINALS: IBdy=3
C     FOUR DIFFERENT GEMINALS: IBdy=4
C
      Dimension IG(4)
C
      IG(1)=IGP
      IG(2)=IGR
      IG(3)=IGQ
      IG(4)=IGS
C
      Do I=1,4
C
      Do J=I+1,4
      If(IG(I).Eq.IG(J)) IG(I)=0
      EndDo      
C
C     enddo I
      EndDo
C
      IBdy=0
      Do I=1,4
      If(IG(I).Ne.0) IBdy=IBdy+1 
      EndDo
C
      Return
      End

*Deck IBodySort
      Subroutine IBodySort(IBdy,IGP,IGR,IGQ,IGS)
C
      Implicit Real*8 (A-H,O-Z)
C
C     DETERMINES IF IP,IR,IQ,IS BELONG TO:
C     ONE GEMINAL: IBdy=1
C     TWO DIFFERENT GEMINALS: IBdy=2
C     THREE DIFFERENT GEMINALS: IBdy=3
C     FOUR DIFFERENT GEMINALS: IBdy=4
C
C     SORT IGP,IGR,IGQ
C
      Dimension IG(4)
C
      IG(1)=IGP
      IG(2)=IGR
      IG(3)=IGQ
      IG(4)=IGS
C
      Do I=1,4
C
      Do J=I+1,4
      If(IG(I).Eq.IG(J)) IG(I)=0
      EndDo
C
C     enddo I
      EndDo
C
      IBdy=0
      Do I=1,4
      If(IG(I).Ne.0) IBdy=IBdy+1
      EndDo
C
      Do I=1,4
      MinG=1000
C
      Do J=I,4
      If(IG(J).Lt.MinG) Then
      IMin=J
      MinG=IG(J)
      EndIf
      EndDo
C
      II=IG(I)
      IG(I)=MinG
      IG(IMin)=II
C
      EndDo 
C
      IGP=IG(1)
      IGR=IG(2)
      IGQ=IG(3)
      IGS=IG(4)
C
      Return
      End

*Deck FragEcorr
      Subroutine FragEcorr(ETot,ENuc,ECorrTot,EGOne,EigVecR,Eig,ABPLUS,
     $           ABMIN,UNOAO,Occ,TwoNO,URe,XOne,IndAux,NBasis,NInte1,
     $           NInte2,NDim,NGem,NGOcc,IFlag,NFrag)
C
C     A ROUTINE FOR COMPUTING ONE- AND TWO-BODY INTERGEMINAL CORRELATION ENERGIES 
C     FOR FRAGMENTS! (defined manually)
C     ONE-BODY IS ONLY DEFINED FOR GVB ! 
C
C     IFlag = 0  : A+,A- MATRICES NOT TRUNCATED
C           = 1  : A+,A- MATRICES ARE TRUNCATED TO INCLUDE ONLY INDICES
C                        OF ORBITALS BELONGING TO GEMINALS NEEDED FOR 
C                        A GIVEN ONE- TWO- THREE- OR FOUR-BODY INTERACTION
C           = 2  : OLD VERSION OF EERPA
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0, Big=1.5D0)
C
      Dimension
     $ EGOne(NGem),
     $ EigVecR(2*(NDim+NBasis)*2*(NDim+NBasis)),
     $ Eig(2*(NDim+NBasis)),IndAux(NBasis),ABPLUS(NDim*NDim),
     $ ABMIN(NDim*NDim),Occ(NBasis),TwoNO(NInte2),UNOAO(NBasis,NBasis),
     $ XOne(NInte1),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Character*5 ATCENT(100),CENT
      Dimension HNO(NInte1),Work(NBasis,NBasis),XKinNO(NInte1),
     $ XNucNO(NInte1),IndAA(NBasis),GCHAR(NGEM,100),IFragG(NGEM),
     $ IMonomerG(NBasis),EGOneTwo(NGem*NGem)
C
      If(IGVB.Ne.1) Stop 'Fatal error: FragEcorr works only for GVB'
C
C     SET A VARIABLE FROM commons TO BE USED IN EInterG
C
      IFrag=1
C
C     COPY IGem TO IAuxGem (in commons) TO BE USED IN EInterG
C
      NGemSave=NGem
      Do I=1,NBasis
      IAuxGem(I)=IGem(I)
      EndDo
C
C     DEFINE FRAGMENTS BY ASSIGNING ORBITALS TO THEM
C
c      Open(10,File='frag_assign.dat')
c      Read(10,*) NGem
c      Read(10,*) (IndAA(I),I=1,NELE)
c      Read(10,*) ((IConnect(I,J),I=1,NGem),J=1,NGem)
c      Close(10) 
C
C     NEW: frag_assign.dat is produced in DALTON now. It includes:
C     number_geminal, no_atom, type_atom, charge_on_atom
C
      Open(10,File='frag_assign.dat')
C
C     READ DATA
C
  888 Read(10,'(2I4,A5,F7.2)',END=999)IG,ICENT,CENT,CHARGE
      NoIG=IG
      NoICENT=ICENT
      If(ICENT.Gt.100) Stop 'Fatal Error: too many atomic centers!'
      GCHAR(IG,ICENT)=CHARGE
      ATCENT(ICENT)=CENT
      GoTo 888
  999 Continue
      Close(10)
C
      Write(6,'(2X,"Atomic Center",I3,": ",A5)')
     $ (I,ATCENT(I),I=1,NoICENT)
C
      IFragG(1:NGem)=0
      NoFrag=0
C     first inactive geminal (by assumption the first geminal is inactive)
      InActG=0
      Do I=1,NELE
      If(Occ(I).Eq.One) Then
      NoFrag=1
      IFragG(I)=1
      InActG=InActG+1
      EndIf
      EndDo
C
c XXX
      If(IFlag.Ne.2) Then
C
      Do IG=1,NoIG
C
      If(IFragG(IG+InActG).Eq.0) Then
C
      NoFrag=NoFrag+1
      IFragG(IG+InActG)=NoFrag
C     check if other geminals belong to the same fragment
C     at least two atomic centers must be shared
      Do JG=IG+1,NoIG
C
      IShare=1
      Do IC=1,NoICENT
      If((GCHAR(IG,IC).Gt.1.D-1.And.GCHAR(JG,IC).Lt.1.D-1).
     $Or.(GCHAR(JG,IC).Gt.1.D-1.And.GCHAR(IG,IC).Lt.1.D-1))
     $ IShare=0
      EndDo
      If(IShare.Eq.1) Then
      If(IFragG(JG+InActG).Ne.0) Then
      Stop
     $'Fatal Error: A geminal cannot belong to more than one fragment'
      Else 
      IFragG(JG+InActG)=NoFrag
      EndIf
      EndIf
C
      EndDo
C
      EndIf
C  
      EndDo
C
c XXX
c if IFlag==2 - OLD EERPA => assign geminals to two fragments
      Else
C
      Do IG=1,NoIG
      If(IFragG(IG+InActG).Eq.0) Then
C
      NoFrag=NoFrag+1
      IFragG(IG+InActG)=NoFrag
C
C     check if other geminals belong to the same fragment
      Do II=1,2*NoIG
      Do JG=IG+1,NoIG
      If(IFragG(JG+InActG).Eq.0) Then
      IShare=0
      Do KG=1,NoIG
      If(IFragG(KG+InActG).Eq.NoFrag) Then
      Do IC=1,NoICENT
      If(GCHAR(KG,IC).Gt.0.1.And.GCHAR(JG,IC).Gt.0.1) IShare=1
      EndDo
      EndIf
      EndDo
      If(IShare.Eq.1) IFragG(JG+InActG)=NoFrag
      EndIf
      EndDo
      EndDo
C
      EndIf
      EndDo
C
c if iflag.ne.2
      EndIf
C
C     Assign geminals to monomers (nonoverlapping fragments) 
C
      IMonomerG(1:NBasis)=4
      Do I=1,NELE
      If(Occ(I).Eq.One) Then
      NoFragM=1
      IMonomerG(I)=1
      EndIf
      EndDo
      Do IG=1,NoIG
      IMonomerG(IG+InActG)=0
      EndDo 
C   
      Do IG=1,NoIG
      If(IMonomerG(IG+InActG).Eq.0) Then
C
      NoFragM=NoFragM+1
      IMonomerG(IG+InActG)=NoFragM
C
C     check if other geminals belong to the same fragment
      Do II=1,2*NoIG
      Do JG=IG+1,NoIG
      If(IMonomerG(JG+InActG).Eq.0) Then
      IShare=0
      Do KG=1,NoIG
      If(IMonomerG(KG+InActG).Eq.NoFragM) Then
      Do IC=1,NoICENT
      If(GCHAR(KG,IC).Gt.1.D-1.And.GCHAR(JG,IC).Gt.1.D-1) IShare=1
      EndDo
      EndIf
      EndDo
      If(IShare.Eq.1) IMonomerG(JG+InActG)=NoFragM
      EndIf
      EndDo
      EndDo
C
      EndIf
C
      Do I=1,NELE
      IndAA(I)=IMonomerG(IGem(I)) 
      EndDo
C
      EndDo
      Do I=1,NELE
      IA=IndAA(I)
      IMonomerG(I)=IA
      IMonomerG(IFindG(I))=IA
      EndDo
C
  555 FORMAT(/,2X,'Geminal no',I3,' in fragment:',I3)
C     set connectivity matrix for active fragments
      Do IG=1,NoIG
      Write(6,555)IG+InActG,IFragG(IG+InActG)
      Write(6,'(2X,"Atomic Centers:")')
      Do IC=1,NoICENT
      If(GCHAR(IG,IC).Gt.1.D-1) Write(6,'(16X,A4)')ATCENT(IC)
      EndDo
      EndDo
c XXX
      If(IFlag.Eq.2) Then
      IFlag=1
      IFl12=1
      EndIf
C
      Do IG=1,NoIG
      Do JG=IG,NoIG
C
      IConnect(IFragG(IG+InActG),IFragG(JG+InActG))=0
      Do IC=1,NoICENT
      If(GCHAR(IG,IC).Gt.1.D-1.And.GCHAR(JG,IC).Gt.1.D-1) 
     $ IConnect(IFragG(IG+InActG),IFragG(JG+InActG))=1
      EndDo     
C
      IConnect(IFragG(JG+InActG),IFragG(IG+InActG))=
     $ IConnect(IFragG(IG+InActG),IFragG(JG+InActG))
C
      EndDo
      EndDo
C     inactive fragment is connected with all the others
      If(InActG.Ne.0) Then 
      Do I=1,NoFrag
      IConnect(1,I)=1
      IConnect(I,1)=1
      EndDo
      EndIf
C
      Do I=1,NELE
      IndAA(I)=IFragG(IGem(I))
      EndDo 
C
      NGem=NoFrag
C      Write(*,*)'Orbital assignment to frags read from frag_assign.dat'
      Write(*,*)'No of fragments:',NGem
c      Write(*,*)'Assignemnts to fragments:',(IndAA(I),I=1,NELE)
      Write(*,*)'Fragment connectivity matrix'
      Do I=1,NGem
      Do J=I+1,NGem
      Write(6,*)I,J,IConnect(I,J)
      EndDo
      EndDo
C 
      If(IGVB.Eq.1) Then
      NGem=NGem+1
      Do I=1,NBasis
      IGem(I)=NGem
      EndDo
      EndIf
C
      Do I=1,NELE
      IA=IndAA(I)
      IGem(I)=IA
      IGem(IFindG(I))=IA
      EndDo
C
      Write(6,'(2X,''NUMBER OF FRAGMENTS: '',I3)') NGem-1
c      Write(6,'(2X,"Orb",2X,"Occupancy",X,"Fragment")')
c      Do I=1,NBasis
c      Write(6,'(X,I3,E16.6,I6)') I,Occ(I),IGem(I)
c      EndDo
      Write(6,'(2X,"Orb",2X,"Occupancy",X,"Fragment",X,"Monomer")')
      Do I=1,NBasis
      Write(6,'(X,I3,E16.6,2I6)') I,Occ(I),IGem(I),IMonomerG(I)
      EndDo
C
      Call OneTwoBody(ETot,ENuc,ECorrTot,EGOne,EGOneTwo,EigVecR,Eig,
     $ ABPLUS,ABMIN,Occ,TwoNO,IndAux,
     $ NBasis,NInte1,NInte2,NDim,NGem,NGOcc,IFlag)
      GoTo 444
C
C     FIND THE INTERFRAGMENT EXCHANGE AND COULOMB INTERACTION
C
      Write(6,'(/,X,"INTER-FRAGMENT COULOMB AND EXCHANGE CONTRIBUTIONS"
     $ )')
C
      IStart=1
c      If(NGOcc.Ne.0) Then
c      IStart=NGOcc+1
c      Write(6,'(X,I2," CORE ORBITALS ARE EXCLUDED")')NGOcc
c      EndIf
C
      Do IF1=1,NGem-1
      Do IF2=1,IF1-1
C
      ECoul=Zero
      Exch=Zero
C
      IJ=0
      Do I=IStart,NBasis
C
      Do J=IStart,I
      IJ=IJ+1
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
C
      If((IGem(I).Eq.IF1.And.IGem(J).Eq.IF2).Or.
     $ (IGem(I).Eq.IF2.And.IGem(J).Eq.IF1)) Then
C
      ECoul=ECoul+FacIJ*Occ(I)*Occ(J)*Two*TwoNO(NAddr3(I,I,J,J))
      Exch=Exch-FacIJ*Occ(I)*Occ(J)*TwoNO(NAddr3(I,J,I,J))
C
      EndIf
C
      EndDo
      EndDo
C
      Write(6,'(X,2I4,2F16.8)') IF1,IF2,ECoul,Exch  
C
      EndDo
      EndDo
C
      GoTo 444
C
C     READ THE EKIN AND EPOT AO MATRICES AND GET CONTRIBUTIONS TO EKIN AND EPOT FROM FRAGMENTS
C
      Open(20,File='ekin.dat',Status='Old')
C
      Read(20,'(A10)')Aux1
      Read(20,'(A10)')Aux1
C
      Do I=1,NBasis
      Do J=1,NBasis/5
      Read(20,222) (Work(I,(J-1)*5+K),K=1,5)
  222 Format(5(F15.8,1X))
      EndDo
      MxK=NBasis-(NBasis/5)*5
      IStart=(NBasis/5)*5
      If(MxK.Ne.0) Read(20,222) (Work(I,IStart+K),K=1,MxK)
      EndDo
      Close(20)
C
      IJ=0
      Do I=1,NBasis 
      Do J=1,I
      IJ=IJ+1
      XKinNO(IJ)=Work(I,J)
      EndDo
      EndDo
      Call MatTr(XKinNO,UNOAO,NBasis)
C
C     XNucNO
C       
      Open(20,File='epot.dat',Status='Old')
C
      Read(20,'(A10)')Aux1
      Read(20,'(A10)')Aux1
C
      Do I=1,NBasis
      Do J=1,NBasis/5
      Read(20,222) (Work(I,(J-1)*5+K),K=1,5)
      EndDo
      MxK=NBasis-(NBasis/5)*5
      IStart=(NBasis/5)*5
      If(MxK.Ne.0) Read(20,222) (Work(I,IStart+K),K=1,MxK)
      EndDo
      Close(20)
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      XNucNO(IJ)=Work(I,J)
      EndDo
      EndDo
      Call MatTr(XNucNO,UNOAO,NBasis)
C
C     CHECKING
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      Diff=Abs(HNO(IJ)-XKinNO(IJ)-XNucNO(IJ))
      If(Diff.Gt.1.D-5) Then
      Write(*,*)HNO(IJ),XKinNO(IJ)+XNucNO(IJ)
      Stop 'Fatal Error 2 in FragEcorr'
      EndIf
C
      EndDo
      EndDo
C      
      Write(6,'(/,X,"INTRA-FRAGMENT EKIN AND ENUC CONTRIBUTIONS"
     $ )')
C
      IStart=1
      If(NGOcc.Ne.0) Then
      IStart=NGOcc+1
      Write(6,'(X,I2," CORE ORBITALS ARE EXCLUDED")')NGOcc
      EndIf
C 
      EOneTot=Zero 
      Do IF=1,NGem-1
      EKin=Zero
      EPot=Zero
C
      Do I=IStart,NBasis
      II=(I*(I+1))/2
C
      If(IGem(I).Eq.IF) Then 
      EKin=EKin+Two*Occ(I)*XKinNO(II)
      EPot=EPot+Two*Occ(I)*XNucNO(II)
      EndIf
C
      EndDo
C
      Write(6,'(X,I4,2F16.8)') IF,EKin,EPot
C
      EOneTot=EOneTot+EKin+EPot
      EndDo
      Write(6,'(X,"TOTAL",F16.8)')EOneTot
C
  444 Continue
C
C     If new EERPA is called and partitioning of the correlation energy into monomer-monomer is needed
C     than continue
C
C     WATCH OUT: orbitals are assigned to monomers IMonomerG(I)
C                and to fragments by IFrag(I)
C     ASSIGN FRAGMENTS TO MONOMERS
      Do I=1,NELE
      IFragG(IGem(I))=IMonomerG(I)
      EndDo
C
C     LOOP OVER FRAGMENTS
C
      Write(6,'(/,X, " FRAGMENT  IN   MONOMER")')
      Do IG=1,NGem-1
      Write(6,'(X,2I4)') IG,IFragG(IG)
      EndDo
C
      Write(6,'(/,2X,
     $ "*** 2-body correlation for fragments in different monomers")')
C
      SumTwo=Zero
      Do IGG1=1,NGem-1
      Do IGG2=1,IGG1-1
C
C     if fragments are in different monomers then      
C   
      If(InActG.Ne.0) Then
      If(IFragG(IGG1).Ne.1.And. IFragG(IGG2).Ne.1.
     $ And.IFragG(IGG1).Ne.IFragG(IGG2)) Then
      Write(6,'(X,2I4,F15.8)')IGG1,IGG2,EGOneTwo(IGG1+(IGG2-1)*NGem)
      SumTwo=SumTwo+EGOneTwo(IGG1+(IGG2-1)*NGem)
      EndIf
      EndIf
C
      If(InActG.Eq.0) Then
      If(IFragG(IGG1).Ne.IFragG(IGG2)) Then
      Write(6,'(X,2I4,E17.6)')IGG1,IGG2,EGOneTwo(IGG1+(IGG2-1)*NGem)
      SumTwo=SumTwo+EGOneTwo(IGG1+(IGG2-1)*NGem)
      EndIf
      EndIf
C
      EndDo
      EndDo
C
      Write(6,'(X," Sum : ",1X,F15.8)') SumTwo
C
      NFrag=NGem
      If(IGVB.Eq.1) NFrag=NGem-1
      NGem=NGemSave
      Do I=1,NBasis
      Hlp=IGem(I)
      IGem(I)=IAuxGem(I)
      IAuxGem(I)=Hlp
      EndDo
C
      IFrag=0
C
      Return
      End

*Deck RDM2ERPA
      Subroutine RDM2ERPA(EigVecR,Eig,ABMIN,TwoNO,NInte2,IndN,Occ,Title,
     $ NBasis,NDimX,NGem,NDim) 
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
C     make sure the params below are the same as in EneERPA, otherwise 
C     2RDM does not yield the same energy as calculated in EneERPA
      Parameter(SmallE=1.D-3,BigE=1.D2)
C
      Include 'commons.inc'
C    
      Dimension
     $ ABMIN(NDimX,NDimX),EigVecR(NDimX*NDimX),
     $ Eig(NDimX),Occ(NBasis),IndN(2,NDim),TwoNO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension TNU(NBasis,NBasis,NDimX),AuxPair(NBasis,NBasis)
C
      TNU(1:NBasis,1:NBasis,1:NDimX)=Zero      
      AuxPair(1:NBasis,1:NBasis)=Zero
C
C     COMPUTE TRANSITION DENSITY MATRIX ELEMENTS TNU
C
      Do NU=1,NDimX
C
      If(Eig(NU).Gt.SmallE.And.Eig(NU).Lt.BigE) Then
C   
C     CHECK THE NORM TO BE SURE IT IS OK
C 
      XNormXY=Zero
      Do I=1,NDimX
C
      IP=IndN(1,I)
      IQ=IndN(2,I)
      AuxPair(IP,IQ)=One
      AuxPair(IQ,IP)=One
C
      YPQ=EigVecR((NU-1)*NDimX+I)
C
C     COMPUTE XPQ
C
      XPQ=Zero
      Do J=1,NDimX
      XPQ=XPQ+One/Eig(NU)*ABMIN(I,J)*EigVecR((NU-1)*NDimX+J)
      EndDo
C
      XNormXY=XNormXY+Two*YPQ*XPQ
C
      TNU(IP,IQ,NU)=Half*( (CICoef(IP)+CICoef(IQ))*YPQ
     $ - (CICoef(IP)-CICoef(IQ))*XPQ)
      TNU(IQ,IP,NU)=Half*( (CICoef(IP)+CICoef(IQ))*YPQ
     $ + (CICoef(IP)-CICoef(IQ))*XPQ)
C
      EndDo
C
      If(Abs(XNormXY-One).Gt.1.D-8) 
     $ Write(6,*)'Wrong norm of Y,X!', NU,XNormXY     
C
C     if eig(nu)
C
      EndIf
c     enddo NU
      EndDo
C
C     FULL 2-RDM = Gamma++++ + Gamma+-+- FOR ERPA-APSG/GVB
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
      K=0
    6 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 6
      EndIf
C
      FName(K:K+14)='_ERPA_2RDM.dat'
      Open(30,File=FName)
      Write(30,'(X,"ERPA 2RDM IN THE NO REPRESENTATION")')
      Write(30,'(X,"P Q R S 2RDM_pqrs")')
C
      EeeAPSG=Zero
      Eee=Zero
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
C
      RDM2=Zero
      RDM2APSG=Zero
C
C     INTRA PAIR
C
      If(IP.Eq.IQ.And.IR.Eq.IS.And.IGem(IP).Eq.IGem(IR))
     $ RDM2=RDM2+CICoef(IP)*CICoef(IR)
C
C     INTER COULOMB AND EXCHANGE
C
      If(IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2+2.0D0*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      RDM2APSG=RDM2
C
C     INTER ERPA
C
      If(.NOT.(IGem(IP).Eq.IGem(IQ).And.IGem(IP).Eq.IGem(IR).And.
     $ IGem(IP).Eq.IGem(IS))) Then
C
      If(IP.Eq.IS.And.IQ.Eq.IR) 
     $RDM2=RDM2-Half*Occ(IQ)*(One-Occ(IP))*AuxPair(IP,IR)*AuxPair(IQ,IS)
C
      Do NU=1,NDimX
      RDM2=RDM2+TNU(IP,IR,NU)*TNU(IS,IQ,NU)
      EndDo
C     if(not....
      EndIf
C
      Write(30,'(4I3,E14.6)')IP,IQ,IR,IS,RDM2
C
      EeeAPSG=EeeAPSG+RDM2APSG*TwoNO(NAddr3(IP,IR,IQ,IS)) 
      Eee=Eee+RDM2*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Close(30)
C
      Write(6,
     $ '(/,1X," *** ERPA 2RDM CONSTRUCTED AND SAVED IN A FILE ***")')
C
      Write(6,
     $ '(X,"APSG ELECTRON INTERACTION COMPUTED FROM 2RDM",F17.8)')
     $ EeeAPSG
      Write(6,
     $ '(X,"ERPA CORRELATION COMPUTED FROM 2RDM         ",F17.8)')
     $ Eee-EeeAPSG
C
      Return
      End

*Deck IFindG
      Integer Function IFindG(IO)
C
      Implicit Real*8 (A-H,O-Z)
      Include 'commons.inc'
C
      IFindG=0
      Do I=1,2*NELE
      If(IAuxGem(IO).Eq.IAuxGem(I).And.IO.Ne.I) IFindG=I
      EndDo
C
      If(IFindG.Eq.0) Then 
      IFindG=IO
      EndIf
C
      Return
      End

*Deck AB_CAS
      Subroutine AB_CAS(ABPLUS,ABMIN,ETot,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
C
C     COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
C
C     DOES THE SAME AS Gamma2_AB BUT MORE EFFICIENTLY
C     
C     RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
C     THE FOLLOWING SYMMETRY IS ASSUMED
C     RDM2(ij,kl) = RDM2(kl,ij)
C     ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
C     SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ IPair(NBasis,NBasis),IndN(2,NDim),IndX(NDim)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1)
C
      Do I=1,NDim
      Do J=1,NDim
      ABPLUS(I,J)=Zero
      ABMIN(I,J)=Zero
      EndDo
      EndDo

C
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo

C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I 
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
!$OMP CRITICAL(crit_AB_CAS_1)
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,*,End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X 
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
!$OMP END CRITICAL(crit_AB_CAS_1)
C
C      Write(*,*) 'RDM2Act',norm2(RDM2Act)
C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero 
      Do I=1,NBasis
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*HNO(II)      
      EndDo
C
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETot=ETot+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)')ETot
C
C     CONSTRUCT A+,A- MATRICES FOR A GIVEN 2-RDM
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C      Write(*,*) 'HNO-Ka', norm2(HNO)

C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN      
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C    
C     AUXILIARY MATRIX AuxI  
C   
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      AuxIO(IPQ)=Zero
      Do IT=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) AuxTwo=ACAlpha
      AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndDo
      EndDo
      EndDo
C
C
C      Write(*,*) 'AuxI-Ka',norm2(AuxI)
C      Write(*,*) 'AuxIO-Ka',norm2(AuxIO)
C
C     AUXILIARY MATRIX WMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
      WMAT(I,J)=Zero
      EndDo
      EndDo 
C
      Do IP=1,NBasis
      Do IR=1,NOccup
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) AuxTwo=ACAlpha
C
      WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))
     $ *FRDM2(IW,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)*AuxTwo
     $ +TwoNO(NAddr3(IT,IU,IP,IW))
     $ *FRDM2(IW,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)*AuxTwo
C
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
C      Write(*,*) 'WMAT-Ka',norm2(WMAT)
C
      Do IRow=1,NDimX
C
      IR=IndN(1,IRow)
      IS=IndN(2,IRow)
      IRS=IndX(IRow)  
C
      Do ICol=1,NDimX
C
      IPP=IndN(1,ICol)
      IQQ=IndN(2,ICol)
C
      Do IP=IQQ,IPP,IPP-IQQ
      Do IQ=IQQ,IPP,IPP-IQQ
      If(IP.Ne.IQ) Then
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=One
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=ACAlpha
      AuxPQRS=AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq+Occ(IP)*Occ(IR)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      AuxTwoSQTU=One
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0) AuxTwoSQTU=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*
     $  FRDM2(IP,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoSQTU
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*
     $ FRDM2(IP,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoSQTU
C
      EndDo
      EndDo
C
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxIO(IQS) 
C
      EndIf
      EndIf
C
C     T3+T4 
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq+Occ(IQ)*Occ(IS)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      AuxTwoUTPR=One
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) AuxTwoUTPR=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*
     $ FRDM2(IS,IT,IQ,IU,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoUTPR
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*
     $ FRDM2(IS,IT,IU,IQ,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoUTPR
C
      EndDo
      EndDo
C
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxIO(IPR)
C
      EndIf
C     
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      AuxTwoPTSU=One
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0) AuxTwoPTSU=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*
     $ FRDM2(IT,IU,IQ,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoPTSU
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T6      
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
C
      Else
C
      Do ITT=1,NAct
      Do IUU=1,NAct
      IT=Ind1(ITT)
      IU=Ind1(IUU)
C
      AuxTwoTQUR=One
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0) AuxTwoTQUR=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*
     $ FRDM2(IS,IP,IU,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoTQUR
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      IPQ=IndX(ICol)
C
      If(IR.Gt.IS.And.IP.Gt.IQ) Then
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)+Arspq
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)+Arspq
      EndIf
C    
      If(IR.Gt.IS.And.IQ.Gt.IP) Then
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)+Arspq
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)-Arspq
      EndIf
C
c     If(IP.Ne.IQ)
      EndIf
C     end of IP,IQ LOOPS
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Write(*,*)'ABPLUS,ABMIN CONSTRUCTED'
C
      Do IRow=1,NDimX
C
      IR=IndN(1,IRow)
      IS=IndN(2,IRow)
      IRS=IndX(IRow)
C
      Do ICol=1,NDimX
C
      IP=IndN(1,ICol)
      IQ=IndN(2,ICol)
      IPQ=IndX(ICol)
C
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If((C(IP)-C(IQ))*(C(IR)-C(IS)).Ne.Zero)
     $ ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)/(C(IP)-C(IQ))/(C(IR)-C(IS))
C
      EndDo
      EndDo
C
C      Print*, "AB-Ka",norm2(ABPLUS),norm2(ABMIN)
C
      Deallocate(RDM2Act)
C
      Return
      End

*Deck AB_T_CAS
      Subroutine AB_T_CAS(ABPLUS,ABMIN,ETot,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,ACAlpha)
C
C     COMPUTE TRIPLET A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
C
C     RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
C     THE FOLLOWING SYMMETRY IS ASSUMED
C     RDM2(ij,kl) = RDM2(kl,ij)
C     ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
C     SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ IPair(NBasis,NBasis),IndN(2,NDim),IndX(NDim)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: RDM2Act(:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1)
C
      Do I=1,NDim
      Do J=1,NDim
      ABPLUS(I,J)=Zero
      ABMIN(I,J)=Zero
      EndDo
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind1(I)=INActive+I
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C 
      If(ICASSCF.Eq.1) Then
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,*,End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
C
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
      GoTo 10
   40 Continue
      Close(10)
C
      Else
C
      NOccup=0
      Do I=1,NBasis
      If(Occ(I).Ne.Zero) NOccup=NOccup+1
      EndDo
C
      EndIf
C
C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*HNO(II)
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETot=ETot+FRDM2ST(0,IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)')ETot
C
C
C     CONSTRUCT A+,A- MATRICES FOR A GIVEN 2-RDM
C
      If(ICASSCF.Eq.1) Then
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
      Else
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C
      EndIf  
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C      Write(*,*) 'HNO-Ka', norm2(HNO)

C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRIX AuxI
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      AuxIO(IPQ)=Zero
      Do IT=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) AuxTwo=ACAlpha
      AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndDo
      EndDo
      EndDo
C
C
C      Write(*,*) 'AuxI-Ka',norm2(AuxI)
C      Write(*,*) 'AuxIO-Ka',norm2(AuxIO)
C
C     AUXILIARY MATRIX WMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
      WMAT(I,J)=Zero
      EndDo
      EndDo
C
      Do IP=1,NBasis
      Do IR=1,NOccup
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) AuxTwo=ACAlpha
C
      WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))
     $ *FRDM2ST(0,IW,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)*AuxTwo
     $ +TwoNO(NAddr3(IT,IU,IP,IW))
     $ *FRDM2ST(0,IW,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)*AuxTwo
C
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do IRow=1,NDimX
C
      IR=IndN(1,IRow)
      IS=IndN(2,IRow)
      IRS=IndX(IRow)
C
      Do ICol=1,NDimX
C
      IPP=IndN(1,ICol)
      IQQ=IndN(2,ICol)
C
      Do IP=IQQ,IPP,IPP-IQQ
      Do IQ=IQQ,IPP,IPP-IQQ
      If(IP.Ne.IQ) Then
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=One
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=ACAlpha
      AuxPQRS=AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
c      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
cC
c      Arspq=Arspq+Occ(IP)*Occ(IR)*AuxTwoPQRS*
c     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
c      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
cC
c      Else
C
c      Do ITT=1,NAct
c      Do IUU=1,NAct
c      IT=Ind1(ITT)
c      IU=Ind1(IUU)
      Do ITT=1,NOccup
      Do IUU=1,NOccup
      IT=ITT
      IU=IUU
C
      AuxTwoSQTU=One
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0) AuxTwoSQTU=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*
     $  FRDM2ST(0,IP,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoSQTU
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*
     $ FRDM2ST(1,IP,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoSQTU
C
      EndDo
      EndDo
C
c      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxIO(IQS)
cC
c      EndIf
      EndIf
C
C     T3+T4
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
c      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
cC
c      Arspq=Arspq+Occ(IQ)*Occ(IS)*AuxTwoPQRS*
c     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
c      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
cC
c      Else
C
c      Do ITT=1,NAct
c      Do IUU=1,NAct
c      IT=Ind1(ITT)
c      IU=Ind1(IUU)
      Do ITT=1,NOccup
      Do IUU=1,NOccup
      IT=ITT
      IU=IUU
C
      AuxTwoUTPR=One
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) AuxTwoUTPR=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*
     $ FRDM2ST(0,IS,IT,IQ,IU,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoUTPR
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*
     $ FRDM2ST(1,IS,IT,IU,IQ,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoUTPR
C
      EndDo
      EndDo
C
c      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxIO(IPR)
cC
c      EndIf
C
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
C
c      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
cC
c      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
cC
c      Else
C
c      Do ITT=1,NAct
c      Do IUU=1,NAct
c      IT=Ind1(ITT)
c      IU=Ind1(IUU)
      Do ITT=1,NOccup
      Do IUU=1,NOccup
      IT=ITT
      IU=IUU
C
      AuxTwoPTSU=One
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0) AuxTwoPTSU=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*
     $ FRDM2ST(1,IT,IU,IQ,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoPTSU
C
      EndDo
      EndDo
C
c      EndIf
      EndIf
C
C     T6
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
C
c      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
cC
c      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
cC
c      Else
C
c      Do ITT=1,NAct
c      Do IUU=1,NAct
c      IT=Ind1(ITT)
c      IU=Ind1(IUU)
      Do ITT=1,NOccup
      Do IUU=1,NOccup
      IT=ITT
      IU=IUU
C
      AuxTwoTQUR=One
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0) AuxTwoTQUR=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*
     $ FRDM2ST(1,IS,IP,IU,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *AuxTwoTQUR
C
      EndDo
      EndDo
C
c      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      IPQ=IndX(ICol)
C
      If(IR.Gt.IS.And.IP.Gt.IQ) Then
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)+Arspq
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)+Arspq
      EndIf
C
      If(IR.Gt.IS.And.IQ.Gt.IP) Then
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)+Arspq
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)-Arspq
      EndIf
C
c     If(IP.Ne.IQ)
      EndIf
C     end of IP,IQ LOOPS
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Write(*,*)'ABPLUS,ABMIN CONSTRUCTED'
C
      Do IRow=1,NDimX
C
      IR=IndN(1,IRow)
      IS=IndN(2,IRow)
      IRS=IndX(IRow)
C
      Do ICol=1,NDimX
C
      IP=IndN(1,ICol)
      IQ=IndN(2,ICol)
      IPQ=IndX(ICol)
C
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If((C(IP)-C(IQ))*(C(IR)-C(IS)).Ne.Zero)
     $ ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)/(C(IP)-C(IQ))/(C(IR)-C(IS))
C
      EndDo
      EndDo
C
C      Print*, "AB-Ka",norm2(ABPLUS),norm2(ABMIN)
C
      Deallocate(RDM2Act)
C
      Return
      End

*Deck Gamma2_AB
      Subroutine Gamma2_AB(ABPLUS,ABMIN,ETot,URe,Occ,XOne,TwoNO,IPair,
     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
C     COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
C     
C     RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
C     THE FOLLOWING SYMMETRY IS ASSUMED
C     RDM2(ij,kl) = RDM2(kl,ij)
C     SO ONLY ELEMENTS ij >= kl ARE STORED
C     SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),RDM2(NBasis**2*(NBasis**2+1)/2),
     $ Ind1(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1)
C
      Do I=1,NDim
      Do J=1,NDim
      ABPLUS(I,J)=Zero
      ABMIN(I,J)=Zero
      EndDo
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      RDM2(1:NRDM2)=Zero
C
C     READ 2RDM, COMPUTE THE ENERGY
C
c      INActive=NELE-NAct
c      NOccup=NELE+NAct
c      Do I=1,2*NAct
c      Ind1(I)=INActive+I 
c      EndDo
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Do I=1,NAct
      Ind1(I)=INActive+I 
      EndDo
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,*,End=40)I,J,K,L,X
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
      RDM2(NAddrRDM(J,L,I,K,NBasis))=Half*X
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      Hlp=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  Hlp=Hlp+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ Hlp=Hlp-Occ(IP)*Occ(IQ)
      If(Hlp.Ne.Zero) RDM2(IAdd)=Hlp
      EndDo
      EndDo
      EndDo
      EndDo
C
C     COMPUTE THE ENERGY FOR CHECKING
C
      ETot=Zero 
      Do I=1,NBasis
      II=(I*(I+1))/2
      ETot=ETot+Two*Occ(I)*HNO(II)      
      EndDo
C
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      ETot=ETot+RDM2(IAdd)*TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)')ETot
C
C     CONSTRUCT A+,A- MATRICES FOR A GIVEN 2-RDM
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN      
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C    
C     AUXILIARY MATRIX AuxI  
C   
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      Do IT=1,NBasis
      AuxTwo=One
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) AuxTwo=ACAlpha
      AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRIX WMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
      WMAT(I,J)=Zero
      EndDo
      EndDo 
C
      Do IP=1,NBasis
      Do IR=1,NOccup
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) AuxTwo=ACAlpha
      WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))*RDM2(NAddrRDM(IW,IU,IT,IR,NBasis))
     $ *AuxTwo
     $ +TwoNO(NAddr3(IT,IU,IP,IW))*RDM2(NAddrRDM(IW,IU,IR,IT,NBasis))
     $ *AuxTwo
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C
      If(IPair(IR,IS).Eq.1) Then
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,NBasis
      If(IP.Gt.IQ) IPQ=IPQ+1
C
      If(IPair(IP,IQ).Eq.1) Then
C
      If(IQ.Gt.IP) IQP=(IQ**2-3*IQ+2)/2+IP
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=One
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=ACAlpha
      AuxPQRS=AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq+Occ(IP)*Occ(IR)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoSQTU=One
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0) AuxTwoSQTU=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*RDM2(NAddrRDM(IP,IU,IR,IT,NBasis))
     $  *AuxTwoSQTU
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*RDM2(NAddrRDM(IP,IU,IT,IR,NBasis))
     $  *AuxTwoSQTU
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T3+T4 
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq+Occ(IQ)*Occ(IS)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoUTPR=One
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) AuxTwoUTPR=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*RDM2(NAddrRDM(IS,IT,IQ,IU,NBasis))
     $  *AuxTwoUTPR
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*RDM2(NAddrRDM(IS,IT,IU,IQ,NBasis))
     $  *AuxTwoUTPR
C
      EndDo
      EndDo
C
      EndIf
C     
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoPTSU=One
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0) AuxTwoPTSU=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*RDM2(NAddrRDM(IT,IU,IQ,IR,NBasis))
     $  *AuxTwoPTSU
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T6      
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoTQUR=One
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0) AuxTwoTQUR=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*RDM2(NAddrRDM(IS,IP,IU,IT,NBasis))
     $  *AuxTwoTQUR
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      If(IR.Gt.IS.And.IP.Gt.IQ) Then
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)+Arspq
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)+Arspq
      EndIf
C    
      If(IR.Gt.IS.And.IQ.Gt.IP) Then
      ABPLUS(IRS,IQP)=ABPLUS(IRS,IQP)+Arspq
      ABMIN(IRS,IQP)=ABMIN(IRS,IQP)-Arspq
      EndIf
C
C     If(IPair(IP,IQ).Eq.1)
      EndIf
      EndDo
      EndDo
C     If(IPair(IR,IS).Eq.1)
      EndIf
      EndDo
      EndDo
C
      Write(*,*)'ABPLUS,ABMIN CONSTRUCTED'
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ ABPLUS(IPQ,IRS)=ABPLUS(IPQ,IRS)/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If((C(IP)-C(IQ))*(C(IR)-C(IS)).Ne.Zero)
     $ ABMIN(IPQ,IRS)=ABMIN(IPQ,IRS)/(C(IP)-C(IQ))/(C(IR)-C(IS))
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      Return
      End

*Deck NAddrRDM
      Integer Function NAddrRDM(Ind1,Ind2,Ind3,Ind4,NBasis)
C
C     A POINTER FOR 2-RDM ELEMENT
C
      Implicit Real*8 (A-H,O-Z)
C
      Ind12=(Ind1-1)*NBasis+Ind2
      Ind34=(Ind3-1)*NBasis+Ind4
      NAddrRDM=Max(Ind12,Ind34)*(Max(Ind12,Ind34)-1)/2+
     $       Min(Ind12,Ind34)
C
      Return
      End

*Deck DE_CAS
      Subroutine DE_CAS(DMAT,DMATM,EMAT,EMATM,URe,Occ,XOne,TwoNO,IPair,
     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
C     COMPUTE THE D AND E MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
C     
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension DMAT(NDim,NBasis),EMAT(NBasis,NBasis),
     $ DMATM(NDim,NBasis),EMATM(NBasis,NBasis),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoNO(NInte2),
     $ IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),RDM2(NBasis**2*(NBasis**2+1)/2),
     $ Ind1(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1)
C    
C     START
C
c herer!!! temporary
     $ ,ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ CMAT(NDim*NDim),DMAT2(NDim,NBasis),EMAT2(NBasis,NBasis),
     $ EMAT3(NBasis,NBasis)
C
C     STOP
C
      Do J=1,NBasis
      Do I=1,NDim
      DMAT(I,J)=Zero
      EndDo
      Do I=1,NBasis
      EMAT(I,J)=Zero
      EndDo
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
      NRDM2 = NBasis**2*(NBasis**2+1)/2
      RDM2(1:NRDM2)=Zero
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Do I=1,NAct
      Ind1(I)=INActive+I
      EndDo
C
      Open(10,File="rdm2.dat",Status='Old')
      Write(6,'(/,1X,''Active block of 2-RDM read from rdm2.dat'')')
C
   10 Read(10,*,End=40)I,J,K,L,X
      I=Ind1(I)
      J=Ind1(J)
      K=Ind1(K)
      L=Ind1(L)
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
      RDM2(NAddrRDM(J,L,I,K,NBasis))=Half*X
C
      GoTo 10
   40 Continue
      Close(10)
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NBasis
      Do IS=1,NBasis
      IAdd=NAddrRDM(IP,IQ,IR,IS,NBasis)
      Hlp=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  Hlp=Hlp+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ Hlp=Hlp-Occ(IP)*Occ(IQ)
      If(Hlp.Ne.Zero) RDM2(IAdd)=Hlp
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN      
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C    
C     AUXILIARY MATRIX AuxI  
C   
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      Do IT=1,NBasis
      AuxTwo=One
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.0) AuxTwo=ACAlpha
      AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*AuxTwo*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndDo
      EndDo
      EndDo
C
C     AUXILIARY MATRIX WMAT
C
      Do I=1,NBasis
      Do J=1,NBasis
      WMAT(I,J)=Zero
      EndDo
      EndDo
C
      Do IP=1,NBasis
      Do IR=1,NOccup
      Do IT=1,NOccup
      Do IW=1,NOccup
      Do IU=1,NOccup
      AuxTwo=One
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.0) AuxTwo=ACAlpha
      WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))*RDM2(NAddrRDM(IW,IU,IT,IR,NBasis))
     $ *AuxTwo
     $ +TwoNO(NAddr3(IT,IU,IP,IW))*RDM2(NAddrRDM(IW,IU,IR,IT,NBasis))
     $ *AuxTwo
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR
      If(IS.Lt.IR) IRS=IRS+1
C
      If(IPair(IR,IS).Eq.1.Or.IS.Eq.IR) Then 
C
      Do IP=1,NBasis
      IQ=IP
C
      IQS=(Max(IS,IQ)*(Max(IS,IQ)-1))/2+Min(IS,IQ)
      IPR=(Max(IR,IP)*(Max(IR,IP)-1))/2+Min(IR,IP)
C
      Arspq=Zero
C
      If(IP.Eq.IR) Arspq=Arspq+(Occ(IP)-Occ(IS))*HNO(IQS)
      If(IS.Eq.IQ) Arspq=Arspq+(Occ(IQ)-Occ(IR))*HNO(IPR)
C
      AuxTwoPQRS=One
      If(IGFact(NAddr3(IP,IQ,IR,IS)).Eq.0) AuxTwoPQRS=ACAlpha
      AuxPQRS=AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
C
C     T1+T2
C
      If(Occ(IP)*Occ(IR).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq+Occ(IP)*Occ(IR)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IP.Eq.IR) Arspq=Arspq+Occ(IP)*AuxI(IQS)
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoSQTU=One
      If(IGFact(NAddr3(IS,IQ,IT,IU)).Eq.0) AuxTwoSQTU=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IS,IQ,IT,IU))*RDM2(NAddrRDM(IP,IU,IR,IT,NBasis))
     $  *AuxTwoSQTU
     $ +TwoNO(NAddr3(IS,IU,IT,IQ))*RDM2(NAddrRDM(IP,IU,IT,IR,NBasis))
     $  *AuxTwoSQTU
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T3+T4 
C
      If(Occ(IQ)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq+Occ(IQ)*Occ(IS)*AuxTwoPQRS*
     $ (Two*TwoNO(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IR,IQ,IS)))
      If(IQ.Eq.IS) Arspq=Arspq+Occ(IQ)*AuxI(IPR)
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoUTPR=One
      If(IGFact(NAddr3(IU,IT,IP,IR)).Eq.0) AuxTwoUTPR=ACAlpha
C
      Arspq=Arspq
     $ +TwoNO(NAddr3(IU,IT,IP,IR))*RDM2(NAddrRDM(IS,IT,IQ,IU,NBasis))
     $  *AuxTwoUTPR
     $ +TwoNO(NAddr3(IU,IR,IP,IT))*RDM2(NAddrRDM(IS,IT,IU,IQ,NBasis))
     $  *AuxTwoUTPR
C
      EndDo
      EndDo
C
      EndIf
C     
      EndIf
C
C     T5
C
      If(Occ(IR)*Occ(IQ).Ne.Zero) Then
C
      If(Occ(IQ).Eq.One.Or.Occ(IR).Eq.One) Then
C
      Arspq=Arspq-Occ(IQ)*Occ(IR)*AuxPQRS
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoPTSU=One
      If(IGFact(NAddr3(IP,IT,IS,IU)).Eq.0) AuxTwoPTSU=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IP,IT,IS,IU))*RDM2(NAddrRDM(IT,IU,IQ,IR,NBasis))
     $  *AuxTwoPTSU
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
C     T6      
C
      If(Occ(IP)*Occ(IS).Ne.Zero) Then
C
      If(Occ(IP).Eq.One.Or.Occ(IS).Eq.One) Then
C
      Arspq=Arspq-Occ(IP)*Occ(IS)*AuxPQRS
C
      Else
C
      Do IT=1,NOccup
      Do IU=1,NOccup
C
      AuxTwoTQUR=One
      If(IGFact(NAddr3(IT,IQ,IU,IR)).Eq.0) AuxTwoTQUR=ACAlpha
C
      Arspq=Arspq
     $ -TwoNO(NAddr3(IT,IQ,IU,IR))*RDM2(NAddrRDM(IS,IP,IU,IT,NBasis))
     $  *AuxTwoTQUR
C
      EndDo
      EndDo
C
      EndIf
      EndIf
C
      If(IS.Eq.IQ) Arspq=Arspq-Half*WMAT(IP,IR)
      If(IP.Eq.IR) Arspq=Arspq-Half*WMAT(IQ,IS)
C
      If(IR.Gt.IS) DMAT(IRS,IP)=DMAT(IRS,IP)+Arspq
      If(IR.Eq.IS) EMAT(IR,IP)=EMAT(IR,IP)+Arspq
C    
      EndDo
C     If(IPair(IR,IS).Eq.1)
      EndIf
      EndDo
      EndDo
C
      Write(*,*)'D AND E MATRICES CONSTRUCTED'
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IQ.Lt.IP) IPQ=IPQ+1
C
      Do IR=1,NBasis
      IS=IR 
C
      If(IQ.Lt.IP) Then
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ DMAT(IPQ,IR)=DMAT(IPQ,IR)/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If(Occ(IR).Eq.One) DMAT(IPQ,IR)=Zero 
      DMATM(IPQ,IR)=DMAT(IPQ,IR)
      EndIf
C
      If(IQ.Eq.IP) Then
      If((C(IP)+C(IQ))*(C(IR)+C(IS)).Ne.Zero)
     $ EMAT(IP,IR)=EMAT(IP,IR)/(C(IP)+C(IQ))/(C(IR)+C(IS))
      If(Occ(IP).Eq.One.Or.Occ(IR).Eq.One) EMAT(IP,IR)=Zero
      EMATM(IP,IR)=EMAT(IP,IR)
      EndIf
C
      EndDo
C
      EndDo
      EndDo
C
c     START
c      Return
C
C     LET US TRY TO "MAP" CAS 2-RDM ON GVB 2-RDM BY COUPLING ORBITALS INTO GEMINALS
C
      If(NELE-INActive.Ne.NAct-(NELE-INActive)) Then 
      Write(6,*) 'Fatal Error: mapping of CAS 
     $ on GVB only defined for CAS(m,m)'
      Stop
      EndIf
C
      NGemSave=NGem
      Do I=1,NBasis
      IGemSave(I)=IGem(I)
      EndDo
      NGem=INActive+NAct/2
      Do I=1,NBasis
      If(Occ(I).Ne.Zero) IGem(I)=0
      If(Occ(I).Eq.Zero) IGem(I)=NGem+1
      EndDo
      Do I=1,INActive
      IGem(I)=I
      EndDo
C
C     COUPLE ORBITALS
C
      Write(6,'(/,X,"Mappinng of CAS(n,n) 2-RDM on GVB-like 2-RDM")')
C
      Do IP=INActive+1,NELE
      IGem(IP)=IP
C
      XDevMax=Zero
      Do IQ=NELE+1,NOccup
C
      XPQPQ=Abs(RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis))-
     $ 2.0D0*Occ(IP)*Occ(IQ))
      XPQQP=Abs(RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis))-
     $ (-Occ(IP)*Occ(IQ)))
      XDev=Half*(XPQPQ+XPQQP)/Occ(IQ)
      If(XDev.Gt.XDevMax) Then
      XDevMax=XDev      
      IQMax=IQ
      EndIf
C
      Write(*,*)IP,IQ
      Write(*,*)'PPQQ',RDM2(NAddrRDM(IP,IP,IQ,IQ,NBasis)),
     $ CICoef(IP)*CICoef(IQ)
      Write(*,*)'PQPQ',RDM2(NAddrRDM(IP,IQ,IP,IQ,NBasis)),
     $ 2.0D0*Occ(IP)*Occ(IQ)
      Write(*,*)'PQQP',RDM2(NAddrRDM(IP,IQ,IQ,IP,NBasis)),
     $ -Occ(IP)*Occ(IQ)
C
      EndDo
C
      If(IGem(IQMax).Eq.0) Then 
      IGem(IQMax)=IP
      Write(6,'(X,"**** Orbital",I2," coupled with",I2 )')IP,IQMax
      Else
      Write(6,*)
     $ "Warning: more than 2 orbitals assigned to a geminal no",IP
      EndIf
C
      EndDo
C
      Do I=1,NGem
C
      Write(6,'(/,X,"Geminal no",I2," includes")')I
      Sum=Zero
C
      Do J=1,NBasis
      If(IGem(J).Eq.I) Then
      Sum=Sum+Occ(J)
      Write(6,'(X,"Orbital No: ",I4)')J
      EndIf
      EndDo
      Write(6,'(X,"Norm: ",F12.6)')Sum
C
      EndDo
C
      Do I=INActive+1,NBasis
      Do J=INActive+1,NBasis
      If(IGem(I).Ne.IGem(J)) EMATM(I,J)=EMATM(I,J)+
     $ Four*C(I)*C(J)*(Two*TwoNO(NAddr3(I,I,J,J))
     $ -TwoNO(NAddr3(I,J,J,I)))
      EndDo
      EndDo
C
      Call APSG_NEST(ABPLUS,ABMIN,CMAT,EMAT2,EMAT3,DMAT2,DMATM,
     $ URe,Occ,XOne,TwoNO,
     $ NBasis,NDim,NInte1,NInte2,NGem,2)
C
      NGem=NGemSave
      Do I=1,NBasis
      IGem(I)=IGemSave(I)
      EndDo
C
C     STOP
C
      Return
      End

*Deck RDMFT_AB
      Subroutine RDMFT_AB(ABPLUS,ABMIN,URe,Occ,XOne,TwoMO,
     $ NBasis,NDim,NInte1,NInte2,ACAlpha)
C
C     COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM CORRESPONDING TO A DMFT FUNCTIONAL
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
      Parameter(Delta=1.D-6)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension EPS(NInte1),C(NBasis),AuxXC(NBasis,NInte1),HNO(NInte1),
     $ IGFact(NInte2)
C BB
c       IFun=2
C GPFBB 
      IFun=22
C BBC1
c      IFun=4
c BBC2
c      IFun=5
C
       Write(6,
     $ '(/,2X,"The following kernel is used in ERPA IFun = ",I4,/)')IFun
C
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
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
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=ACAlpha*HNO(IJ)
C
      Else
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoMO(NAddr3(IT,IT,I,J))-TwoMO(NAddr3(IT,I,IT,J)))
      EndDo
C
      Aux=(One-ACAlpha)*Aux
      HNO(IJ)=HNO(IJ)+Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     CONSTRUCT TWO-ELECTRON PART OF AC ALPHA-HAMILTONIAN      
C
      NAdd=Zero
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      KL=0
      Do K=1,NBasis
      Do L=1,K
      KL=KL+1
C
      If(IJ.Ge.KL) Then
      NAdd=NAdd+1
C
      IGFact(NAdd)=1
      If(.Not.(
     $IGem(I).Eq.IGem(J).And.IGem(J).Eq.IGem(K).And.IGem(K).Eq.IGem(L)))
     $ IGFact(NAdd)=0
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      EPS(IJ)=HNO(IJ)
C
      Do K=1,NBasis
      AuxTwo=One
      If(IGFact(NAddr3(I,J,K,K)).Eq.0) AuxTwo=ACAlpha

      EPS(IJ)=EPS(IJ)+Occ(K)*Two*AuxTwo*TwoMO(NAddr3(I,J,K,K))
      EndDo
C
      Do L=1,NBasis
C
      AuxXC(L,IJ)=Zero
C
      Do K=1,NBasis
      AuxTwo=One
      If(IGFact(NAddr3(I,K,J,K)).Eq.0) AuxTwo=ACAlpha
      AuxXC(L,IJ)=AuxXC(L,IJ)
     $ +GOCC(Occ(L),Occ(K),0,L,K)*AuxTwo*TwoMO(NAddr3(I,K,J,K))
      EndDo
C
      EndDo
C
      EndDo
      EndDo
C
C     COMPUTE ABPLUS AND ABMIN
C
      IRS=0
      Do IR=1,NBasis
      Do IS=1,IR-1
      IRS=IRS+1
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
C
      AuxTwo=One
      If(IGFact(NAddr3(IS,IR,IP,IQ)).Eq.0) AuxTwo=ACAlpha
C
      TwoEl=TwoMO(NAddr3(IS,IR,IP,IQ))
      BRSPQ=(Occ(IR)-Occ(IS))*(Occ(IQ)-Occ(IP))*Two*TwoEl*AuxTwo
      BRSPQ=BRSPQ +
     $  (-GOCC(Occ(IR),Occ(IP),0,IR,IP)
     $   +GOCC(Occ(IR),Occ(IQ),0,IR,IQ)
     $   +GOCC(Occ(IS),Occ(IP),0,IS,IP)
     $   -GOCC(Occ(IS),Occ(IQ),0,IS,IQ))
     $ *TwoMO(NAddr3(IS,IP,IQ,IR))*AuxTwo
C
      ARSPQ=(Occ(IR)-Occ(IS))*(Occ(IP)-Occ(IQ))*Two*TwoEl*AuxTwo
      ARSPQ=ARSPQ +
     $  (-GOCC(Occ(IR),Occ(IQ),0,IR,IQ)
     $   +GOCC(Occ(IR),Occ(IP),0,IR,IP)
     $   +GOCC(Occ(IS),Occ(IQ),0,IS,IQ)
     $   -GOCC(Occ(IS),Occ(IP),0,IS,IP))
     $ *TwoMO(NAddr3(IS,IQ,IP,IR))*AuxTwo
C
      IRQ=(Max(IQ,IR)*(Max(IQ,IR)-1))/2+Min(IQ,IR)
      IPS=(Max(IP,IS)*(Max(IP,IS)-1))/2+Min(IP,IS)
      IQS=(Max(IQ,IS)*(Max(IQ,IS)-1))/2+Min(IQ,IS)
      IPR=(Max(IP,IR)*(Max(IP,IR)-1))/2+Min(IP,IR)
C
      If(IS.Eq.IP) BRSPQ=BRSPQ-(Occ(IR)-Occ(IS))*EPS(IRQ)
      If(IR.Eq.IQ) BRSPQ=BRSPQ+(Occ(IR)-Occ(IS))*EPS(IPS)
C
      If(IS.Eq.IP) BRSPQ=BRSPQ-AuxXC(IR,IRQ)+AuxXC(IS,IRQ)
      If(IR.Eq.IQ) BRSPQ=BRSPQ+AuxXC(IR,IPS)-AuxXC(IS,IPS)
C
      If(IR.Eq.IP) ARSPQ=ARSPQ+(Occ(IR)-Occ(IS))*EPS(IQS)
      If(IS.Eq.IQ) ARSPQ=ARSPQ-(Occ(IR)-Occ(IS))*EPS(IPR)
C
      If(IR.Eq.IP) ARSPQ=ARSPQ+AuxXC(IR,IQS)-AuxXC(IS,IQS)
      If(IS.Eq.IQ) ARSPQ=ARSPQ-AuxXC(IR,IPR)+AuxXC(IS,IPR)
C
      ABPLUS(IRS,IPQ)=ARSPQ+BRSPQ
      ABMIN(IRS,IPQ)=ARSPQ-BRSPQ
C
      If(Abs(Abs(C(IP))-Abs(C(IQ))).Le.Delta*Abs(C(IQ)).Or.
     $   Abs(Abs(C(IR))-Abs(C(IS))).Le.Delta*Abs(C(IR))) Then
      ABPLUS(IRS,IPQ)=Zero
      ABMIN(IRS,IPQ)=Zero
      Else
C
      ABPLUS(IRS,IPQ)=ABPLUS(IRS,IPQ)/(C(IR)+C(IS))/(C(IP)+C(IQ))
      ABMIN(IRS,IPQ)=ABMIN(IRS,IPQ)/(C(IR)-C(IS))/(C(IP)-C(IQ))
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

*Deck FRDM2
      Real*8 Function FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C
C     FOR A GIVEN SET OF INDICES AND THE KNOWN PART OF ACTIVE RDM2 
C     RETURNS THE ELEMENT OF RDM2_PQRS FOR CAS
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension RDM2Act(NAct**2*(NAct**2+1)/2),Occ(NBasis),Ind2(NBasis)
C
      RDM2=Zero
      If(IP.Eq.IR.And.IQ.Eq.IS.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  RDM2=RDM2+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
C     ACTIVE PART
C
      If(Ind2(IP)*Ind2(IQ)*Ind2(IR)*Ind2(IS).Ne.Zero) Then 
      RDM2=RDM2Act(NAddrRDM(Ind2(IP),Ind2(IQ),Ind2(IR),Ind2(IS),NAct))
      EndIf
C
      FRDM2=RDM2
C
      Return
      End 

*Deck FRDM2ST
      Real*8 Function FRDM2ST(ISPin,IP,IQ,IR,IS,RDM2Act,Occ,Ind2,
     $ NAct,NBasis)
C
C     FOR A GIVEN SET OF INDICES AND THE KNOWN PART OF ACTIVE RDM2
C     RETURNS THE ELEMENT OF RDM2_PQRS FOR CAS
C
C     ISpin=0  RDM2=1/2 [++++ + ---- + +-+- + +-+-] 
C     ISpin=1  RDM2=1/2 [++++ + ---- - +-+- - +-+-]         
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension RDM2Act(NAct**2*(NAct**2+1)/2),Occ(NBasis),Ind2(NBasis)
C
      RDM2=Zero
C
      If(ICASSCF.Eq.1) Then
C
      If(ISpin.Eq.0.And.IP.Eq.IR.And.IQ.Eq.IS.And. 
c       If(IP.Eq.IR.And.IQ.Eq.IS.And.
     $ (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $  RDM2=RDM2+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And. (Occ(IP).Eq.One.Or.Occ(IQ).Eq.One) )
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
C     ACTIVE PART
C
      stop 
     $ 'Fatal error in FRDM2ST: active part of spin-2RDM not available!'
c      If(Ind2(IP)*Ind2(IQ)*Ind2(IR)*Ind2(IS).Ne.Zero) Then
c      RDM2=RDM2Act(NAddrRDM(Ind2(IP),Ind2(IQ),Ind2(IR),Ind2(IS),NAct))
c      EndIf
C
      Else 
C
C     2RDM for GVB
C
C     INTRA-GEMINAL
C
      If(IP.Eq.IQ.And.IR.Eq.IS.And.IGem(IP).Eq.IGem(IR))
     $ RDM2=RDM2+(One-Two*ISpin)*CICoef(IP)*CICoef(IR)
c     $ RDM2=RDM2+CICoef(IP)*CICoef(IR)
C
C     INTER COULOMB AND EXCHANGE
C
      If(ISpin.Eq.0.And.IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ))
c      If(IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ)) 
     $ RDM2=RDM2+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      EndIf
C
      FRDM2ST=RDM2
C
      Return
      End

*Deck FRDM2GVB
      Real*8 Function FRDM2GVB(IP,IQ,IR,IS,Occ,NBasis)
C
C     FOR A GIVEN SET OF INDICES AND THE KNOWN CICoef AND IGem
C     RETURNS THE ELEMENT OF RDM2_PQRS FOR GVB
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis)
C
      RDM2=Zero
C
C     INTRA PAIR
C
      If(IP.Eq.IQ.And.IR.Eq.IS.And.IGem(IP).Eq.IGem(IR))
     $ RDM2=RDM2+CICoef(IP)*CICoef(IR)
C
C     INTER COULOMB AND EXCHANGE
C
      If(IP.Eq.IR.And.IQ.Eq.IS.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2+Two*Occ(IP)*Occ(IQ)
      If(IP.Eq.IS.And.IQ.Eq.IR.And.IGem(IP).Ne.IGem(IQ))
     $ RDM2=RDM2-Occ(IP)*Occ(IQ)
C
      FRDM2GVB=RDM2
C     
      Return
      End

