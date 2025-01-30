*Deck ACCAS
      Subroutine ACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  Title,BasisSet,NBasis,NInte1,NInte2,NGem)
C
c     use types
      use print_units
      use timing
      use read_external
      use abfofo
C
C     A ROUTINE FOR COMPUTING ELECTRONIC ENERGY USING ERPA TRANSITION
C     DENSITY MATRIX ELEMENTS
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Character(*) :: BasisSet
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1)
C
C     LOCAL ARRAYS
C
      Dimension
     $ IndX(NBasis*(NBasis-1)/2),IndN(2,NBasis*(NBasis-1)/2),
     $ IndAux(NBasis),IPair(NBasis,NBasis),
     $ XMuMat(NBasis,NBasis)
C
      Double Precision  :: Tcpu,Twall
C
C     START TIMING FOR AC PROCEDURES
C
      call gclock('START',Tcpu,Twall)
      If (ICholeskyOTF==1) Write(6,'(1x," BasisSet = ", A)') BasisSet

c      NoEig=1
c      NDimFull=NBasis*(NBasis-1)/2
c      Call ACPINO(ENuc,TwoNO,Occ,XOne,
c     $ NBasis,NInte1,NInte2,NDimFull,NGem,NoEig)
c      Stop
C
C     CONSTRUCT LOOK-UP TABLES
C
      Do I=1,NELE
      IndAux(I)=0
      EndDo
      Do I=1+NELE,NBasis
      IndAux(I)=2
      EndDo
C
      ICount=0
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
      Write(6,'(/,X," In ACCAS: Active Orbitals ",I4,/)')ICount
C
      IPair(1:NBasis,1:NBasis)=0
C
      Write(LOUT,'(2x,a,4x,2e15.5)') 'Threshold for quasi-degeneracy ',
     $ ThrSelAct

      Write(LOUT,'(2x,a,2e15.5)') 'Threshold for quasi-virtual orbital',
     $ ThrQVirt

      Write(LOUT,'(2x,a,2e14.5)')'Threshold for quasi-inactive orbital',
     $ ThrQInact
C
      If(IFlCore.Eq.0) Then
        If(NCoreOrb.Eq.0) Then
           Do I=1,NBasis
              If(Occ(I).Eq.One) NCoreOrb=NCoreOrb+1
           EndDo
        EndIf
        Write(6,'(/,1X,"*** NO. OF CORE ORBITALS :",I3)') 
     $ NCoreOrb
      EndIf
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
      If((IndAux(I).Eq.1).And.(IndAux(J).Eq.1)
     $ .And.(Abs(Occ(I)-Occ(J))/Occ(I).Lt.ThrSelAct)
     $ ) Then
C
      Write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')I,J
C
      Else
C
C     If IFlCore=0 do not include core (inactive) orbitals
      If((IFlCore.Eq.1).Or.
C     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
     $ (IFlCore.Eq.0.And.I.Gt.NCoreOrb.And.J.Gt.NCoreOrb))    
     $ Then
C
      If(Abs(Occ(i)+Occ(j)-Two).Gt.ThrQInact.And.
     $   Abs(Occ(i)+Occ(j)).Gt.ThrQVirt) Then
      Ind=Ind+1
      IndX(Ind)=Ind
      IndN(1,Ind)=I
      IndN(2,Ind)=J
      IPair(I,J)=1
      IPair(J,I)=1
      EndIf
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
      Write(6,'(/2X,"Number of pairs reduced to:",I6)')Ind
      Write(6,'(2X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,3I5,2E14.4)')I,Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
      If(IFlRESPONSE.Eq.1) Then
C
      Write(6,'(/,X,''Polarizability tensor calculation for Om ''
     $ ,F8.4)') Om
C
      If(Max_Cn.Eq.-1) Then
      Call Polariz(FreqOm,UNOAO,XOne,URe,Occ,
     $   IGem,NAcCAS,NInAcCAS,NELE,NBasis,NInte1,NGem,IndAux,
     $   IndN,IndX,NDimX,ICholesky)
      Else
      Write(6,'(/,X,''Expand C(Om) maximally up to order '',I4)') Max_Cn
      Call PolarizAl(FreqOm,UNOAO,XOne,URe,Occ,
     $   IGem,NAcCAS,NInAcCAS,NELE,NBasis,NInte1,NGem,IndAux,
     $   IndN,IndX,NDimX,ICholesky,Max_Cn)
      EndIf
C
      Return
      EndIf
C
C     COMPUTE RESPONSE RDMs FROM AC0-CAS DERIVATIVE-LIKE EXPRESSION
C
      If (IRedVirt.Eq.1) Then
C
      If(ITwoEl.eq.3) Then
C
      Call RDMResp_FOFO(Occ,URe,UNOAO,XOne,IndN,IndX,IndAux,IGem,
     $                  NAcCAS,NInAcCAS,NDimX,NDim,NBasis,NInte1,
     $                  'FFOO','FOFO',ICholesky,IOrbRelax,IOrbIncl)
C
      If (IOrbRelax==1) Call delfile('FFFO')
      Call delfile('FFOO')
      Call delfile('FOFO')
C
      stop
C
      Else
C
      Write(6,'(/,X,
     $ ''Response RDM is available only with FOFO integrals '')')
      Stop
C
      EndIf
C
      EndIf
C
      If(IFunSR.Eq.0) Then
C
      If(IFlAC0D.Eq.1) Then
C
C this is a version of AC0 which takes special care of deexcitations (SA-CAS is required)
C
      Call RunACDEXIT(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      call gclock('ACD model',Tcpu,Twall)

      ElseIf(IDBBSC.Gt.0) Then
C
C this is a version of AC0 with CBS basis set corrections (SR integrals are required)
C
      Call DBBSCH(ETot,ENuc,URe,Occ,XOne,UNOAO,
     $   BasisSet,NBasis,NInte1,IndN,IndX,NDimX)

      Call delfile('cholvecs')
      Call delfile('cholvErf')
      call gclock('DBBSC',Tcpu,Twall)
C
      Else
C
      Call RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)

      call gclock('AC  model',Tcpu,Twall)
C
C     Corr,md correlation energy with local mu, requires Cholesky vecs
C
      If (IFlCorrMD.Eq.1) Then
c     If (ICholesky==0) Stop "Run SRAC0 with Cholesky!"
c      If(ITwoEl.eq.1) then
c      Call LOC_MU_CBS(XMuMAT,URe,UNOAO,Occ,TwoNO,NBasis,NInte2)
c      stop

      Call LOC_MU_CHOL(BasisSet,CorrMD,AvMU,URe,UNOAO,Occ,NBasis)
      Call delfile('cholvecs') ! delete cholesky vecs
      call gclock('LOC_MU_CHOL_v2',Tcpu,Twall)
      EndIf
C
      If (IDBBSC.Eq.1.And.ITwoEl.Eq.1) Then
      Call LOC_MU_CBS(XMuMat,URe,UNOAO,Occ,TwoNO,NBasis,NInte2)
      EndIf
C
      EndIf
C
      ElseIf(IFunSR.Eq.3) Then
C
      Call RunDFOnTop(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      Else
C
      Call RunACCASLR(BasisSet,ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      call gclock('ACLR model',Tcpu,Twall)
C
C     Corr,md correlation energy with local mu, requires Cholesky vecs
C
      If (IFlCorrMD.Eq.1) Then
c     If (ICholesky==0) Stop "Run SRAC0 with Cholesky!"
      Call LOC_MU_CHOL(CorrMD,AvMU,URe,UNOAO,Occ,NBasis,0)
      Call delfile('cholvecs')
      call gclock('LOC_MU_CHOL_v2',Tcpu,Twall)
      EndIf
C
      EndIf
C
      Return
      End

* Deck RunACCAS
      Subroutine RunACCAS(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      use abfofo
      use ab0fofo
      use read_external
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem), EGOne(NGem)
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C     IFlSnd  = 1 - run AC0-CAS (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0-CAS
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1)
     $ Write(6,'(/," *** ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1)
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1)
     $ included in ERPA correlation ***",/)')
      EndIf
C
C     CALL AC If IFlAC=1 OR IFlSnd=1
C
      If(IFlAC.Eq.1.Or.IFlSnd.Eq.1) Then
      NGOcc=0
      If(IFlACFREQ.Eq.0.And.IFlACFREQNTH.Eq.0.
     $ And.IFlAC1FREQNTH.Eq.0) Then
      Call ACECORR(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      Else
      Call ACIter(ETot,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      EndIf
c
c exact AC
c      NoEig=1
c      NDimFull=NBasis*(NBasis-1)/2
c      Call ACPINO(ENuc,TwoNO,Occ,XOne,
c     $ NBasis,NInte1,NInte2,NDimFull,NGem,NoEig)
c
c AC with varying Alpha-dependent RDMs
c      NDimFull=NBasis*(NBasis-1)/2
c      Call ACRDM(ETot,ENuc,TwoNO,Occ,XOne,
c     $ UNOAO,IndN,IndX,IndAux,NDimX,NBasis,NInte1,NInte2,NDimFull,NGem)
C
      If(ITwoEl.Eq.3) Then
      If (ICholesky .Eq. 0) Then
      Call delfile('FFOO')
      Call delfile('FOFO')
      ElseIf (IFlCorrMD.eq.0.and.IDBBSC.Eq.0) Then
      Call delfile('cholvecs')
      EndIf
      EndIf
C
      Return
      EndIf
C
C     CALCULATE THE A+B AND A-B MATRICES
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      ACAlpha=One
c herer!!!
c      ACAlpha=zero
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
      If(ITwoEl.Eq.3) Then
      Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,
     $ NInte1,'FFOO','FOFO',ICholesky,0,ACAlpha,.false.)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
      EndIf
C
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,X,"***************************** ")')
      Write(6,'(  X,"*** ERPA-CAS CALCULATIONS *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(NoSt.Eq.1) Then
c herer!!!
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
c      Call ERPASYMMXY(EigVecR,EigX,Eig,ABPLUS,ABMIN,Occ,IndN,
c     $ NDimX,NBasis)
c       do i=1,ndimx
c       write(*,*)'I',indn(1,i),indn(2,i),eig(i)
c       do j=1,ndimx
c       yy=eigvecr((i-1)*NDimX+j)
c       xx=eigx((i-1)*NDimX+j)
c       if(abs(yy).gt.1.d-7) write(*,*)'y',i,j,yy
c       if(abs(xx).gt.1.d-7) write(*,*)'x',i,j,xx
c       enddo
c       enddo
c
c       stop

      Else
      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      EndIf
      Write(6,'(/,
     $ " *** ERPA-CAS Excitation Energies (a.u., eV) *** ")')
C     the purpose of sorting is only to print a few highest (sorted) eigenvectors
      Call SortEig(1,Eig,ABPLUS,EigVecR,NDimX)
      Do I=1,10
      Write(6,'(I4,4X,2E16.6)') I,Eig(I),27.211*Eig(I)
      EndDo
C
      If(ITwoEl.Eq.1) Then
      Write(6,'(/," *** Computing ERPA 2-RDM *** ")')
      Call RDM2FULL(EigVecR,Eig,ABMIN,TwoNO,NInte2,IndN,
     $ Occ,Title,NBasis,NDimX,NGem,NDim)
      EndIf
C
      Write(6,'(/," *** Computing ERPA energy *** ",/)')
C
      If(ITwoEl.Eq.3) Then
      Call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Occ,
     $ IGem,IndN,IndX,NAcCAS+NInAcCAS,
     $ NDimX,NBasis,'FOFO',ICholesky)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      EndIf
      ECorr=Ecorr*Half
      Write
     $ (6,'(/,1X,''ECASSCF+ENuc, AC1-Corr, ERPA-CASSCF'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
C
      Return
      End

* Deck RunACDEXIT
      Subroutine RunACDEXIT(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
c      use abmat
      use abfofo
      use ab0fofo
      use read_external
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      CHARACTER(100) :: num1char
      Character*32 Str
C
      Include 'commons.inc'
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis),NSymAO(NBasis),
     $ NSymNO(NBasis),MultpC(8,8),NumOSym(8),NumStSym(16),IStSy(16),
     $ IStCAS(2,100),EExcit(100),ICORR(100),ECorrSym(100)
C
C
C     ***************************************
      If(ISymmAC0D.Eq.1) Then
C     ***************************************
C
C     AC0D corrections are computed for states of the same symmetries as
C     the symmetries of SA-CAS states
C
C
C     Symmetry of NO's
C
      Call read_sym_molpro(NSymAO,MxSym,NumOSym,
     $                    'MOLPRO.MOPUN','CASORB  ',NBasis)
C
      Call sym_inf_molpro('2RDM',NumOSym,NSym,NumStSym,IStSy,
     $                    NStSym,NSymAO,NBasis)
C     number of irreps
      write(*,*)'NSym',NSym
C     number of atomic orbitals in each irrep
      write(*,*)'NumOSym',NumOSym(1:NSym)
C     number of irreps in SA-CAS
      write(*,*)'NStSym',NStSym
C     number of states in SA-CAS in each irrep
      write(*,*)'NumStSym',NumStSym(1:NStSym)
C     which irrep (a label) in SA-CAS
      write(*,*)'IStSy',IStSy(1:NStSym)
C
      Write(6,'(2/,X,"SYMMETRY OF NATURAL ORBITALS")')
      Write(6,'(X,"Orbital",3X,"Symmetry",3X,"Occupancy")')
C
      Do I=1,NBasis
      NSymNO(I)=0
C
      Do J=1,NBasis
C
      If(Abs(UNOAO(I,J)).Gt.0.1) Then
C
      ISym=NSymAO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym)
     $ Write(6,'("In RunACDEXIT. Symm of NO cannot be established",
     $ I3)')I
      EndIf
C
      EndIf
      EndDo
      Write(6,'(I4,7X,I4,2X,E16.6)')I,NSymNO(I),Occ(I)
      EndDo
C
C     checking
      Do I=1,MxSym
      II=0
      Do IOrb=1,NBasis
      If(NSymNO(IOrb).Eq.I) II=II+1
      EndDo
      If(II.Ne.NumOSym(I))
     $ Write(*,*) 'In RunACDEXIT. Symmetry of NO cannot be established!'
      EndDo
C
      If(MxSym.Eq.1) Then
      MultpC(1,1)=1
      Else
      Do I=1,MxSym
      Do J=1,I
      MultpC(I,J)=IEOR(I-1,J-1)+1
      MultpC(J,I)=MultpC(I,J)
      EndDo
      EndDo
      EndIf
C
c      write(*,*) 'st, sym',inst(1,1),inst(2,1)
C
C     determine a number of states in the SA-CAS calculation
C
C     grep ' !MCSCF STATE' filename.out > sacas_ene.dat
      Open(10,File="sacas_ene.dat",Status='Old')
      NoStMx=0
      I=0
   20 I=I+1
      Read(10,'(A14,I1,A1,I1,A15,F22.12)',End=40)
     $Str,IStCAS(1,I),Str,IStCAS(2,I),Str,EExcit(I)
      Write(6,'(X,"SA-CAS Energy  ",I1,".",I1,F22.12)')
     $ IStCAS(1,I),IStCAS(2,I),EExcit(I)
      Read(10,*)
      NoStMx=NoStMx+1
      GoTo 20
   40 Continue
      Close(10)
      Write(6,'(/,X,"The number of states in SA-CAS: ",I4,/)')NoStMx
C
      Do I=1,NoStMx
      If(IStCAS(1,I).Eq.inst(1,1).And.IStCAS(2,I).Eq.inst(2,1)) ICAS=I
      EndDo
      Write(6,'(X,"Energy of the NoSt ",I1,".",I1,F22.12)')
     $ IStCAS(1,ICAS),IStCAS(2,ICAS),EExcit(ICAS)
C
      Do I=1,NoStMx
      ICORR(I)=0
      If(EExcit(I).Ge.EExcit(ICAS)) ICORR(I)=1
      EndDo
C
      If(ITwoEl.eq.1) Then
C
      Call AC0DSYMM(ICAS,NoStMx,ICORR,EExcit,IStCAS,NSym,NSymNO,MultpC,
     $ ECorrSym,
     $ ETot,TwoNO,Occ,URe,XOne,
     $ UNOAO,IndN,IndX,NBasis,NDimX,NInte1,NInte2)
C
      ElseIf(ITwoEl.eq.3) Then
C
      Call Y01CASDSYM_FOFO(ICAS,NoStMx,ICORR,EExcit,IStCAS,NSym,NSymNO,
     $ MultpC,
     $ ECorrSym,Occ,URe,XOne,
     $ 'PROP0','PROP1',
     $ 'XY0',UNOAO,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,
     $ NBasis,NDimX,NInte1,NoSt,'EMPTY','FFOO',
     $ 'FOFO',ICholesky,ETot,IFlAC0DP)
C
      EndIf
C
      Write(6,*)
      Write
     $ (6,'(X,''ECASSCF+ENuc, AC0DE-Corr, AC0DE-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorrSym(ICAS),ETot+ENuc+ECorrSym(ICAS)
      Write(6,*)

      Do I=1,NoStMx
      If(ICORR(I).Eq.1.And.I.Ne.ICAS)
     $  Write
     $ (6,'(X,''Dexcitation correction for AC0 for '',
     $ I1,".",I1,"->",I1,".",I1," transition ",F15.8)')
     $ IStCAS(1,I),IStCAS(2,I),IStCAS(1,ICAS),IStCAS(2,ICAS),ECorrSym(I)
      EndDo
C
C     ***************************************
      ElseIf(ISymmAC0D.Eq.0) Then
C     ***************************************
C
C     AC0D corrections will be computed for states which have best overlap
C     with SA-CAS states
C
      NDimD=0
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1
C
      IJ=IJ+1
C
      If(IndAux(I)+IndAux(J).Ne.0.And.IndAux(I)+IndAux(J).Ne.4) Then
C     If IFlCore=0 do not include core (inactive) orbitals
      If((IFlCore.Eq.1).Or.
     $ (IFlCore.Eq.0.And.Occ(I).Ne.One.And.Occ(J).Ne.One)) Then
      If(Abs(Occ(i)+Occ(j)-Two).Gt.1.D-10) Ind=Ind+1
      EndIf
c     If(IndAux(I)+IndAux(J).Ne.0 ...
      EndIf
C
      EndDo
      EndDo
C
      NDimD=Ind
C
C     determine a number of states in the SA-CAS calculation
C
C     grep ' !MCSCF STATE' filename.out > sacas_ene.dat
      Open(10,File="sacas_ene.dat",Status='Old')
      NoStMx=0
      I=0
   22 Read(10,'(A32,F22.12)',End=42) Str,EE
      Write(6,'(X,"SA-CAS Energy  ",I2,F22.12)')
     $ NoStMx+1,EE
      Read(10,*)
      NoStMx=NoStMx+1
      GoTo 22
   42 Continue
      Close(10)
      Write(6,'(/,X,"The number of states in SA-CAS: ",I4,/)')NoStMx
C
      Do IH0St=NoSt,NoStMx
C
      Write(6,'(/,X,86("*"))')
      Write(6,'(X,"**** AC0D calculation for SA-CAS state no",
     $ I3," with H0 constructed for state no",I4,
     $ " ****")')
     $ NoSt,IH0St
      Write(6,'(X,86("*"),/)')
C
      If(ITwoEl.eq.1) Then
C
      Call AC0CDEXCIT(IH0St,ECorr,ETot,TwoNO,Occ,URe,XOne,
     $ UNOAO,IndN,IndX,NBasis,NDimX,NDimD,NInte1,NInte2)
C
      ElseIf(ITwoEl.eq.3) Then
C
      Call Y01CASD_FOFO(IH0St,Occ,URe,XOne,
     $ 'PROP0','PROP1',
     $ 'XY0',UNOAO,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,
     $ NBasis,NDimX,NInte1,NoSt,'EMPTY','FFOO',
     $ 'FOFO',ICholesky,ETot,ECorr)
C
C     ITwoEl
      EndIf
C
      If(IH0St.Eq.NoSt) Then
C
C     ECorr stores AC0 correlation correction for the NoSt state
C     (NoSt is specified in the input.inp)
C
      Write
     $ (6,'(/,X,''ECASSCF+ENuc, AC0DE-Corr, AC0DE-CASSCF '',4X,3F15.8)')
     $ ETot+ENuc,ECorr,ETot+ENuc+ECorr
C
      Else
C
      Write
     $ (6,'(/,X,''Dexcitation correction for AC0 for '',I2," ->",I2,
     $ " transition ",F15.8)')
     $ IH0St,NoSt,ECorr
C
      EndIf
C
C     Do IH0St=NoSt,NoStMx
      EndDo
C
C     ***************************************
      EndIf
C     ***************************************
C
      Return
      End

* Deck RunACCASLR
      Subroutine RunACCASLR(BasisSet,ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
      use sorter
      use abmat
      use abfofo
      use ab0fofo
      use read_external
      use grid_internal
      use timing
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Character(*) :: BasisSet
      Include 'commons.inc'
C
      Character*60 FName
      Real*8, Dimension(:), Allocatable :: TwoEl2
C      Real*8, Dimension(:,:), Allocatable :: OrbGrid
C      Real*8, Dimension(:,:), Allocatable :: OrbXGrid
C      Real*8, Dimension(:,:), Allocatable :: OrbYGrid
C      Real*8, Dimension(:,:), Allocatable :: OrbZGrid
      Real*8,allocatable,target :: PhiGGA(:,:,:)
      Real*8,allocatable,target :: PhiLDA(:,:)
      Real*8,pointer :: OrbGrid(:,:)
      Real*8,pointer :: OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:)

      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: SRKer
      Real*8, Dimension(:), Allocatable :: Work
      Dimension NSymNO(NBasis),VSR(NInte1),MultpC(15,15),NumOSym(15),
     $ IndInt(NBasis)
C
c     Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis),
     $ VCoul(NInte1),Den(NBasis,NBasis),Work2(NInte1)
C
C     LOCAL ARRAYS
C
      Integer(8) :: MemSrtSize
      Integer :: jtsoao(NBasis)
      Integer :: NSymBas(8),NSymOrb(8)
C
      Logical :: doGGA, doGGAdal
C
      Character(*),Parameter :: griddalfile='dftgrid.dat'
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem),EGOne(NGem),
     $ UAux(NBasis,NBasis),VecAux(NBasis)
      Real*8, Dimension(:,:), Allocatable :: Jmat
C
C     IFlAC   = 1 - adiabatic connection formula calculation
C               0 - AC not used
C     IFlSnd  = 1 - run AC0-CAS (linerized in alpha, MP2-like expression for AC is used)
C             = 0 - do not run AC0-CAS
C     IFlCore = 1 - core (inactive) orbitals included in ERPA correlation
C             = 0 - core (inactive) orbitals excluded from ERPA correlation
C
C     PRINT FLAGS
C
      If(IFlAC.Eq.1)
     $ Write(6,'(/," *** LR ADIABATIC CONNECTION CALCULATIONS ***",/)')
C
      If(IFlCore.Eq.0) Then
      Write(6,'(/," *** IFlCore=0: Inactive orbitals (n_p=1)
     $ excluded from ERPA correlation ***",/)')
      Else
      Write(6,'(/," *** IFlCore=1: Inactive orbitals (n_p=1)
     $ included in ERPA correlation ***",/)')
      EndIf
      Write(6,'(/,X,''****************************************'',
     $            ''***************************************'')')
C
      Write(6,'(/,1X,"The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
C
      Allocate  (TwoEl2(NInte2))
C
C     load orbgrid and gradients, and wgrid
C
      If (IFunSR == 4) Then ! POSTCAS
         If (IFunSR2 == 1) doGGA = .false.
         If (IFunSR2 > 1)  doGGA = .true.
      Else ! REGULAR
         If (IFunSR == 1) doGGA = .false.
         If (IFunSR > 1)  doGGA = .true.
      End If
C
      If (InternalGrid==0) then
C
      If (IMOLPRO == 1) then
      Write(LOUT,'(/1x,a)') 'MOLPRO GRID '
C
C     read no of grid pts
      Call molprogrid0(NGrid,NBasis)
C
      Else If (IDALTON == 1) then
      Write(LOUT,'(/1x,a)') 'DALTON GRID '
C     read no of grid pts
      Call daltongrid0(NGrid,griddalfile,doGGAdal,NBasis)
      if (doGGA.neqv.doGGAdal) then
         print*, 'doGGA GammCor =', doGGA
         print*, 'doGGA Dalton  =', doGGAdal
         stop "error in RunACCASLR!"
      endif
C
      EndIf
C
      Write(6,'(1X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (WGrid(NGrid))
      Allocate  (Work(NGrid))
C
c      Call EKT(URe,Occ,XOne,TwoNO,NBasis,NInte1,NInte2)
C
      If (doGGA) then
         allocate(PhiGGA(NGrid,NBasis,4))
         OrbGrid  => PhiGGA(:,:,1)
         OrbXGrid => PhiGGA(:,:,2)
         OrbYGrid => PhiGGA(:,:,3)
         OrbZGrid => PhiGGA(:,:,4)
      ElseIf (.not.doGGA) then
         allocate(PhiLDA(NGrid,NBasis))
         OrbGrid  => PhiLDA
         OrbXGrid => PhiLDA
         OrbYGrid => PhiLDA
         OrbZGrid => PhiLDa
      EndIf
C
      If (IMOLPRO == 1) Then
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $                WGrid,UNOAO,NGrid,NBasis)
      ElseIf (IDALTON == 1) Then
      UAux=transpose(UNOAO)
      If (doGGA) then
         write(lout,'(/1x,a)') 'GRID TYPE: GGA'
         call daltongrid_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                       WGrid,NGrid,griddalfile,NBasis)
         call daltongrid_tran_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                       UAux,0,NGrid,NBasis,.false.)
      Else
         write(lout,'(/1x,a)') 'GRID TYPE: LDA'
         call daltongrid_lda(OrbGrid,WGrid,NGrid,griddalfile,NBasis)
         call daltongrid_tran_lda(OrbGrid,UAux,0,NGrid,NBasis,.false.)
      EndIf ! doGGA
      print*, 'OrbGrid  =',norm2(OrbGrid)
      print*, 'OrbXGrid =',norm2(OrbXGrid)
      print*, 'OrbYGrid =',norm2(OrbYGrid)
      print*, 'OrbZGrid =',norm2(OrbZGrid)
      EndIf ! Interface
C
      ElseIf (InternalGrid==1) then

      Write(6,'(1x,a,i3)') 'IFunSR       =',IFunSR
c     Write(6,'(1x,a,i3)') 'IUnits       =',IUnits
c     Write(6,'(1x,a,i3)') 'InternalGrid =',InternalGrid
      Write(6,'(1x,a,i3)') 'IGridType    =', IGridType
      Write(6,'(1x,a,i3)') 'ORB_ORDERING =', IOrbOrder
C
      If (doGGA) Then
         Write(LOUT,'(/1x,a)') 'INTERNAL GRID GGA'
         Call internal_gga_no_orbgrid(IGridType,BasisSet,
     $                          IOrbOrder,transpose(UNOAO),
     $                          WGrid,PhiGGA,NGrid,NBasis,NBasis,IUnits)
         OrbGrid  => PhiGGA(:,:,1)
         OrbXGrid => PhiGGA(:,:,2)
         OrbYGrid => PhiGGA(:,:,3)
         OrbZGrid => PhiGGA(:,:,4)
      ElseIf(.not.doGGA) Then
         Write(LOUT,'(/1x,a)') 'INTERNAL GRID LDA'
         Call internal_lda_no_orbgrid(IGridType,BasisSet,
     $                          IOrbOrder,transpose(UNOAO),
     $                          WGrid,PhiLDA,NGrid,NBasis,NBasis,IUnits)
         OrbGrid  => PhiLDA
         OrbXGrid => PhiLDA
         OrbYGrid => PhiLDA
         OrbZGrid => PhiLDa
      End If
C
      Allocate  (Work(NGrid))
C
      End If ! InternalGrid
C
C     set/load symmetries of NO's
C
C      Open(10,File='indices_int.dat')
C      Read(10,*) MxSym
C      IStart=0
C      Do I=1,MxSym
C      Read(10,*) X
C      NumOSym(I)=X
C      EndDo
C      Close(10)
C
      If (IMOLPRO == 1) Then
      Call create_ind_molpro('2RDM',NumOSym,IndInt,NSym,NBasis)
C      Print*, 'NumOSym Molpro'
C      Do I=1,15
C         Print*, I,NumOSym(I)
C      EndDo
      ElseIf (IDALTON ==1) Then
      Call read_sym_dalton(NSym,NSymBas,NSymOrb,'SIRIUS.RST','BASINFO ')
c      Print*, 'NSymBas, NSymOrb Dalton'
      NumOSym=0
      NumOSym(1:NSym)=NSymOrb(1:NSym)
c     If (NSym.gt.1) stop "Dalton with Symmetry in RunACCASLR!"
      EndIf
      MxSym=NSym
C
      If (ICholesky==1) Then
         if (IDALTON==1) stop "Finish Dalton+Chol in RunACCASLR!"
         Call read_aosao_map_molpro(jtsoao,
     $              'MOLPRO.MOPUN','CASORBAO',NBasis)
         If (InternalGrid==0.and.MxSym>1.and.IFunSRKer==1) Then
            Stop "In RunACCASLR: Kernel & Sym /= 1 : use INTERNAL GRID!"
         EndIf
      EndIf

      NSymNO(1:NBasis)=0
      IStart=0
      Do I=1,MxSym
      Do J=IStart+1,IStart+NumOSym(I)
C
      If (ICholesky==0) Then
         Do IOrb=1,NBasis
         If(Abs(UNOAO(IOrb,J)).Gt.1.D-1) NSymNO(IOrb)=I
         EndDo
      ElseIf (ICholesky==1) Then
         Do IOrb=1,NBasis
         If(Abs(UNOAO(IOrb,jtsoao(J))).Gt.1.D-1) NSymNO(IOrb)=I
         EndDo
      EndIf
C
      EndDo
      IStart=IStart+NumOSym(I)
      EndDo
C
C     checking
      Do I=1,MxSym
      II=0
      Do IOrb=1,NBasis
      If(NSymNO(IOrb).Eq.I) II=II+1
      EndDo
      If(II.Ne.NumOSym(I))
     $ Write(*,*) 'In RunACCASLR. Symmetry of NO cannot be established!'
      EndDo
C
      If(MxSym.Eq.1) Then
      MultpC(1,1)=1
      Else
      Do I=1,MxSym
      Do J=1,I
      MultpC(I,J)=IEOR(I-1,J-1)+1
      MultpC(J,I)=MultpC(I,J)
      EndDo
      EndDo
      EndIf
C
C     if IFunSR.Gt.0 (but different from 3) TwoEl includes lr-integrals with erf/r
C     load full-range two-electron integrals and transform to NO
C
C      Do I=1,60
C      FName(I:I)=' '
C      EndDo
C      K=0
C    5 K=K+1
C      If (Title(K:K).Ne.' ') Then
C      FName(K:K)=Title(K:K)
C      GoTo 5
C      EndIf
C      FName(K:K+10)='.reg.integ'
C      Call Int2_AO(TwoEl2,NumOSym,MultpC,FName,NInte1,NInte2,NBasis)
      MemSrtSize=MemVal*1024_8**MemType
      If (ICholesky==0) Then
c     this case is run if CBS[H] is computed for AC0 and incore
      If (IDBBSC.Eq.2.And.ITwoEl.Eq.1) Then
      Call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
      Call LoadSaptTwoEl(4,TwoEl2,NBasis,NInte2)
      Else ! regular LR-AC0 calculation
      If (IMOLPRO==1) Then
         Call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT',MemSrtSize)
      ElseIf (IDALTON==1) Then
         Call readtwoint(NBasis,1,'AOTWOINT','AOTWOSORT',MemSrtSize)
      EndIf
      If(ITwoEl.Eq.1) Call LoadSaptTwoEl(3,TwoEl2,NBasis,NInte2)
      EndIf
C
      If(ITwoEl.Eq.1) Then
      If(IDBBSC.Eq.2) Then
      Write(6,'(" Transforming two-electron erf integrals ...")')
      Else
      Write(6,'(" Transforming full two-electron integrals ...")')
      EndIf
      Call TwoNO1(TwoEl2,UNOAO,NBasis,NInte2)
C
      ElseIf(ITwoEl.Eq.3) Then
C     PREPARE POINTERS: NOccup=num0+num1
      Call prepare_nums(Occ,Num0,Num1,NBasis)
C     TRANSFORM J AND K
      UAux=transpose(UNOAO)
      Call tran4_gen(NBasis,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        NBasis,UAux,
     $        'FFOO','AOTWOSORT')
      Call tran4_gen(NBasis,
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        NBasis,UAux,
     $        Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)),
     $        'FOFO','AOTWOSORT')
C
      EndIf ! ITwoEl
      EndIf ! ICholesky
C
      If (IDALTON==1) Then
      NAct=NAcCAS
      INActive=NInAcCAS
      Call PsiHPsi(EPsiHPsi,Occ,XOne,ENuc,INActive,NAct,NInte1,NBasis)
      EndIf
C
C     compute sr potential (vsr=xc+hartree)
C     as a byproduct a sr energy (ensr=sr-xc+sr-hartree) is obtained
C
      IFunSave=IFunSR
      If(IFunSR.Eq.4) IFunSR=IFunSR2
C
      Den=0
      Work2=0
      Do I=1,NBasis
      VecAux=UNOAO(I,:)
      Call dger(NBasis,NBasis,1d0*Occ(I),
     $          VecAux,1,VecAux,1,Den,NBasis)
C     $          UNOAO(I,:),1,UNOAO(I,:),1,Den,NBasis)
      EndDo
C
      IJ = 0
      Do J=1,NBasis
      Do I=1,J
      IJ = IJ + 1
      Work2(IJ) = Den(I,J)
      EndDo
      EndDo
C
      Write(6,'(/,1X,"RANGE PARAMETER ",F8.3)')Alpha
C
C     calculate short-range Coulomb matrix in AO
      If (ICholesky==0) Then
         Call PotCoul_mithap(VCoul,Work2,.true.,'AOERFSORT',NBasis)
C         block
C         double precision :: JMOsr(NBasis,NBasis)
C         call triang_to_sq2(VCoul,JMOsr,NBasis)
C         Print*, 'JMOsr in AO basis =',norm2(JMOsr)
C         do j=1,NBasis
C            write(6,'(*(f13.8))') (JMOsr(i,j),i=1,NBasis)
C         enddo
C         end block
      ElseIf (ICholesky==1) Then
         ! for Cholesky, Jmat(sr) is read from disk
         Allocate(Jmat(NBasis,NBasis))
         Open(newunit=iunit,file='jsrmat',form='unformatted')
         Read(iunit) NBasis2
         If(NBasis2.ne.NBasis) stop "NBasis error in RunACCASLR"
         Read(iunit) Jmat
         Close(iunit,status='delete')
         Call sq_to_triang2(Jmat,VCoul,NBasis)
c         Print*, 'JMOsr in AO basis =',norm2(Jmat)
c         do j=1,NBasis
c            write(6,'(*(f13.8))') (Jmat(i,j),i=1,NBasis)
c         enddo
         Deallocate(Jmat)
      EndIf
C      Print*, 'UNOAO =',norm2(UNOAO)
C      do j=1,NBasis
C         write(6,'(*(f13.8))') (UNOAO(i,j),i=1,NBasis)
C      enddo
C
      Print*, 'VCoul',norm2(VCoul)
      Call EPotSR(EnSR,EnHSR,VSR,Occ,URe,UNOAO,.false.,
     $        OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,
     $        NSymNO,VCoul,Alpha,IFunSR,
C     $        NSymNO,VCoul,TwoEl2,TwoNO,Alpha,IFunSR,
     $        NGrid,NInte1,NInte2,NBasis)
C
      If (IDALTON==1) Then
        Print*, 'Dalton: add VHsr[rho] contribution...'
        XOne = XOne + VSR
      EndIf
C      block
C        double precision :: VHsr(NBasis,NBasis)
C        call triang_to_sq2(VSR,VHsr,NBasis)
C        Print*, 'VHsr (NO) GammCor =',norm2(VHsr)
C        do i=1,NBasis
C           write(6,'(*(f13.8))') (VHsr(i,j),j=1,NBasis)
C        enddo
C      end block
C      Call EPotSR(EnSR,EnHSR,VSR,Occ,URe,OrbGrid,OrbXGrid,OrbYGrid,
C     $ OrbZGrid,WGrid,NSymNO,TwoEl2,TwoNO,NGrid,NInte1,NInte2,NBasis)
      Write(6,'(" SR xcH Energy",F15.8)')EnSR
C
C     CALCULATE THE SR_XC_PBE ENERGY WITH "TRANSLATED" ALPHA AND BETA DENSITIES
C     [as in Gagliardi J. Chem. Phys. 146, 034101 (2017)]
C
      Call SR_PBE_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(/," SR_xc_PBE with translated densities",F15.8,/)')
     $ EXCTOP
C
      Call PBE_ONTOP_MD(PBEMD,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(/," SR_PBE_corr_md ",F15.8,/)') PBEMD
C
      If (IFlCorrMD.Eq.1) Then
      Call PBE_ONTOP_C_MD(PBEmodMD,5.d0,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(1x,"SR_PBE_C_corr_md",F15.8,/)') PBEmodMD
      EndIf
C
c      Call CASPI_SR_PBE(URe,Occ,
c     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
      IFunSR=IFunSave
C
      XVSR=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      XVSR=XVSR+Two*Occ(I)*VSR(II)
      EndDo
C
c ???
      If(IFunSR.Eq.4) Then
c      Write(6,'(X,/,
c     $"*** REMOVING VSR_HXC FROM A ONE-ELECTRON HAMILTONIAN*** ",/)')
      Do I=1,NInte1
c ??? uncomment for srcaspi calculations
c      XOne(I)=XOne(I)-VSR(I)
      EndDo
      eone=zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      eone=eone+Two*Occ(I)*xone(II)
      EndDo
      EndIf
C
      Allocate (SRKer(NGrid))
      SRKer(1:NGrid)=Zero
C
C      IFunSRKer = 1
C
      If(IFunSRKer.Eq.1) Then
C
      Write(6,'(/," *** Generating a sr-kernel on a grid ***")')
C
      Call GetSRKer(SRKer,Occ,URe,OrbGrid,WGrid,NBasis,NGrid)
      Do I=1,NGrid
      Work(I)=WGrid(I)*SRKer(I)
      EndDo
C
      EndIf
C
      If(IFlSnd.Eq.1) Then
C
C     ***** LR-AC0 CALCULATION **********************
C
      Write(6,'(  X,"*** LR-AC0-CAS CALCULATION *** ")')
C
      If(ITwoEl.Eq.3) Then
C
      Call Y01CASLR_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,
     $ MultpC,NSymNO,
     $ SRKer,WGrid,OrbGrid,
     $ 'PROP0','PROP1','XY0',
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,
     $ NGrid,NDimX,NBasis,NDimX,NInte1,NoSt,
     $ 'FOFO','FFOOERF','FOFOERF',ICholesky,0,IFunSRKer,ECASSCF,ECorr)
C
      ElseIf(ITwoEl.Eq.1) Then
C
      If (IDBBSC.Eq.2) Then
C
c if CBS[H] is computed: XOne contains h0 without vrs
c                        TwoNO contains 1/r integrals, TwoEl2: erf/r
      Call AC0CASLR(ECorr,ECASSCF,ENuc,TwoNO,Occ,URe,UNOAO,XOne,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
     $ TwoEl2,OrbGrid,Work,NSymNO,MultpC,NGrid)

      Write
     $ (6,'(1X,''CASSCF+ENuc, AC0-CBS[H], Total'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
C
      Return
C
      EndIf
C
C
      Call AC0CASLR(ECorr,ECASSCF,ENuc,TwoNO,Occ,URe,UNOAO,XOne,
     $ ABPLUS,ABMIN,EigVecR,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
     $ TwoEl2,OrbGrid,Work,NSymNO,MultpC,NGrid)

c      Call AC0CASLR(ECorr,ECASSCF,ENuc,TwoEl2,Occ,URe,XOne,
c     $ ABPLUS,ABMIN,EigVecR,Eig,
c     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
c     $ TwoNO,OrbGrid,Work,NSymNO,MultpC,NGrid)

C
      EndIf
C
      Write(6,*)"lrAC0 Correlation",ECorr
      write(*,*)"One-ele Energy",eone-XVSR
C
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC0-Corr, Total'',6X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC0-Corr,Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr
C
      GoTo 777
C
      EndIf
C
C     ****** LR-AC CALCULATION *******************************************************
C
      If(IFlAC.Eq.1.And.IFlSnd.Eq.0) Then
C
      Write(6,'(  X,"*** LR-AC-CAS CALCULATION *** ")')
C
      If(IFunSRKer.Eq.1) Then
C
C     GENERATE A SR KERNEL AND DUMP IT
C
      If(ITwoEl.Eq.3) Then
      Call ModABMin_FOFO(Occ,SRKer,WGrid,OrbGrid,ABMIN,
     $                   MultpC,NSymNO,
     $                   IndN,IndX,NDimX,NGrid,NBasis,
     $                   NAcCAS,NInAcCAS,'FOFO','FOFOERF',
     $                   .false.,'srdump')
C
      ElseIf(ITwoEl.Eq.1) Then
      Open(20,File="srdump",Form='UNFORMATTED')
C
      Do I=1,NBasis
      CICoef(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CICoef(I)=-CICoef(I)
      EndDo
C
      Do IRow=1,NDimX
C
      IA=IndN(1,IRow)
      IB=IndN(2,IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NDimX
C
      IC=IndN(1,ICol)
      ID=IndN(2,ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      XKer1234=Zero
C
      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
      ISym=MultpC(I1I2S,I3I4S)
      If(ISym.Eq.1) Then
      Do I=1,NGrid
      XKer1234=XKer1234+Work(I)*
CC     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
CC     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
c    $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
c    $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)
     $ OrbGrid(I,IA)*OrbGrid(I,IB)*OrbGrid(I,IC)*OrbGrid(I,ID)
      EndDo
      EndIf
C
      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
C
      Write(20) Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
C
      EndDo
      EndDo
C
      Close(20)
C
C     if TwoEl
      EndIf
C
      EndIf
C
      NGOcc=0
      Call ACECORR(ECASSCF,ENuc,TwoNO,URe,Occ,XOne,UNOAO,
     $ IndAux,ABPLUS,ABMIN,EigVecR,Eig,EGOne,
     $ Title,NBasis,NInte1,NInte2,NDimX,NGOcc,NGem,
     $ IndN,IndX,NDimX)
      ECorr=EGOne(1)
C
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC-Corr, Total'',7X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC-Corr, Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr

C
      GoTo 777
C
      EndIf
C
C     ****** LR-AC1 (aka ERPA1) CALCULATION ******************************************
C
C     compute AB matrices with the lr integrals and a modified potential
C     TwoNO are lr-integrals!  xone already includes vsr (if properly saved in molpro)
C
      ACAlpha=One
C
      If(ITwoEl.Eq.3) Then
      Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,
     $ NInte1,'FFOOERF','FOFOERF',ICholesky,0,ACAlpha,.false.)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call AB_CAS(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne,TwoNO,IPair,
     $ IndN,IndX,NDimX,NBasis,NDimX,NInte1,NInte2,ACAlpha)
C
      EndIf
C
      Write(6,'(/,1X,  ''lrCASSCF+ENuc Energy       '',4X,F15.8)')
     $ ECASSCF-XVSR+ENuc
      Write(6,'(1X,  ''Total lrCASSCF+ENuc+srDF Energy'',F15.8,/)')
     $ ECASSCF-XVSR+EnSR+ENuc
C
      Do I=1,NBasis
      CICoef(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) CICoef(I)=-CICoef(I)
      EndDo
C
      If(IFunSRKer.Eq.1) Then
C
C     ADD A SR-ALDA KERNEL
C
      Call gclock('START',Tcpu,Twall)
C
C      WITHOUT BATCHES:
C
C      Do IRow=1,NDimX
CC
C      IA=IndN(1,IRow)
C      IB=IndN(2,IRow)
C      CA=CICoef(IA)
C      CB=CICoef(IB)
CC
C      Do ICol=1,NDimX
CC
C      IC=IndN(1,ICol)
C      ID=IndN(2,ICol)
C      CC=CICoef(IC)
C      CD=CICoef(ID)
CC
C      XKer1234=Zero
CC
C      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
C      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
C      ISym=MultpC(I1I2S,I3I4S)
CC
C      If(ISym.Eq.1) Then
C      Do I=1,NGrid
C      XKer1234=XKer1234+Work(I)*
CC     $ OrbGrid(IA+(I-1)*NBasis)*OrbGrid(IB+(I-1)*NBasis)*
CC     $ OrbGrid(IC+(I-1)*NBasis)*OrbGrid(ID+(I-1)*NBasis)
C     $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
C     $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)
C      EndDo
C      EndIf
CC
C      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
CC
C      ABMIN((ICol-1)*NDimX+IRow)=ABMIN((ICol-1)*NDimX+IRow)
C     $ +Four*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
CC
C      EndDo
C      EndDo
C
C
      If(ITwoEl.Eq.3) Then
      Call ModABMin_FOFO(Occ,SRKer,WGrid,OrbGrid,ABMIN,
     $                   MultpC,NSymNO,
     $                   IndN,IndX,NDimX,NGrid,NBasis,
     $                   NAcCAS,NInAcCAS,'FOFO','FOFOERF',.false.)
C      Print*,'ABMIN-MY',norm2(ABMIN)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call ModABMinSym(Occ,SRKer,WGrid,OrbGrid,TwoEl2,TwoNO,ABMIN,
     $          MultpC,NSymNO,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
C      Print*,'ABMIN-KA',norm2(ABMIN)
C
      EndIf
C
      Write(6,'(1X," *** sr-kernel added ***")')
      Call gclock('sr-kernel',Tcpu,Twall)
C
      EndIf
C
      Write(6,'(/,1X,"*** LR-ERPA-CAS CALCULATION *** ")')
C
C     FIND EIGENVECTORS (EigVecR) AND COMPUTE THE ENERGY
C
      If(NoSt.Eq.1) Then
      Call ERPASYMM1(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      Else
      Call ERPAVEC(EigVecR,Eig,ABPLUS,ABMIN,NBasis,NDimX)
      EndIf
      Write(6,'(/,
     $ " *** LR-CAS-SR-DFT Excitation Energies (a.u., eV) *** ")')
C     the purpose of sorting is only to print a few highest (sorted) eigenvectors
      Call SortEig(1,Eig,ABPLUS,EigVecR,NDimX)
C
      Do IS=1,Min(5,MxSym)
      NExcit=0
C
      Do I=1,NDimX
c      If(NExcit.Lt.5) Then
      If(NExcit.Lt.20) Then
C
C     find symmetry
      XMax=Zero
      JMax=0
      Do J=1,NDimX
      If(Abs(EigVecR((I-1)*NDimX+J)).Gt.XMax) Then
      XMax=Abs(EigVecR((I-1)*NDimX+J))
      JMax=J
      EndIf
      EndDo
      IP=IndN(1,JMax)
      IQ=IndN(2,JMax)
      ISym=MultpC(NSymNO(IP),NSymNO(IQ))
C
      If(ISym.Eq.IS) Then
      NExcit=NExcit+1
      Write(6,'("Excit",I2,I2,4X,2E16.6)')NExcit,ISym,
     $ Eig(I),27.211*Eig(I)
      EndIf
C
      EndIf
      EndDo
C
      EndDo
C
      Write(6,'(/," *** Computing LR-ERPA energy *** ",/)')

      If(ITwoEl.Eq.3) Then
      Call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Occ,
     $ IGem,IndN,IndX,NAcCAS+NInAcCAS,
     $ NDimX,NBasis,'FOFOERF',ICholesky)
C
      ElseIf(ITwoEl.Eq.1) Then
      Call ACEneERPA(ECorr,EigVecR,Eig,TwoNO,URe,Occ,XOne,
     $ IndN,NBasis,NInte1,NInte2,NDimX,NGem)
C
      EndIf
C
      ECorr=Ecorr*Half
      Write
     $ (6,'(1X,''lrCASSCF+srDF+ENuc, lrAC1-Corr, Total'',6X,3F15.8)')
     $ ECASSCF-XVSR+EnSR+ENuc,ECorr,ECASSCF-XVSR+EnSR+ENuc+ECorr
      Del=EnHSR+EXCTOP
      Write
     $ (6,'(1X,''lrCASSCF+srDF[OnTop]+ENuc, lrAC1-Corr,Total'',3F15.8)')
     $ ECASSCF-XVSR+Del+ENuc,ECorr,ECASSCF-XVSR+Del+ENuc+ECorr
C
  777 Continue
C
      DeAllocate  (TwoEl2)
      DeAllocate  (WGrid)
C      DeAllocate  (OrbGrid)
C      DeAllocate  (OrbXGrid)
C      DeAllocate  (OrbYGrid)
C      DeAllocate  (OrbZGrid)
      DeAllocate  (Work)
      DeAllocate (SRKer)
C
      Call delfile('AOTWOSORT')
      Call delfile('AOERFSORT')
C
      If(ITwoEl.Eq.3) Then
      Call delfile('FFOO')
      Call delfile('FOFO')
      Call delfile('FOFOERF')
      Call delfile('FFOOERF')
      EndIf
C
      Return
      End

*Deck GetSRKer
      Subroutine GetSRKer(SRKer,Occ,URe,OrbGrid,WGrid,NBasis,NGrid)
C
C     RETURNS a SR-ALDA KERNEL ON THE GRID
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension SRKer(NGrid),Occ(NBasis),URe(NBasis,NBasis),
C     $ OrbGrid(NBasis,NGrid),WGrid(NGrid)
     $ OrbGrid(NGrid,NBasis),WGrid(NGrid)
C
C     LOCAL ARRAYS
C
      Dimension RhoVec(NGrid)
C
      Do I=1,NGrid
      Call DenGrid(I,Rho,Occ,URe,OrbGrid,NGrid,NBasis)
      RhoVec(I)=Rho
      EndDo
C
      Call RhoKernel(RhoVec,SRKer,2,Alpha,NGrid)
C
      Return
      End

*Deck AC0CASLR
      Subroutine AC0CASLR(ECorr,ETot,ENuc,TwoNO,Occ,URe,UNOAO,XOne,
     $ ABPLUS,ABMIN,EigY,Eig,
     $ IndN,IndX,NDimX,NBasis,NDim,NInte1,NInte2,
     $ TwoEl2,OrbGrid,SRKerW,NSymNO,MultpC,NGrid)
C
C     A ROUTINE FOR COMPUTING AC INTEGRAND
C
      use timing
C      use abmat
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
c
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Integer,Parameter :: Maxlen = 128
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),
     $ XOne(NInte1),Occ(NBasis),TwoNO(NInte2),
     $ IndAux(NBasis),
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ Eig(NDimX),EigY(NDimX*NDimX),IndX(NDim),IndN(2,NDim),
     $ TwoEl2(NInte2),OrbGrid(NBasis*NGrid),SRKerW(NGrid),
     $ NSymNO(NBasis),MultpC(15,15)
C
C     LOCAL ARRAYS
C
      Integer :: Offset,Batchlen
      Real*8, Allocatable :: RDM2Act(:)
      Double Precision, Allocatable :: Work(:),Batch(:,:)
      Dimension C(NBasis),HNO(NInte1),
     $ IGFact(NInte2),
     $ Ind1(NBasis),Ind2(NBasis),Ind3(NBasis),WMAT(NBasis,NBasis),
     $ AuxI(NInte1),AuxIO(NInte1),IPair(NBasis,NBasis),
     $ EigX(NDimX*NDimX),
     $ IEigAddY(2,NDimX),IEigAddInd(2,NDimX),IndBlock(2,NDimX),
     $ XMAux(NDimX*NDimX),XMuMAT(NBasis,NBasis),work1(NBasis,NBasis)
! delete after tests
     $ ,WorkH(NInte1),TwoElSR(NInte2)
! analysis of AC0
       double precision :: ECorrIJ(6,6)
       integer          :: IGIJ(4,4),IGemNo(6,2)
C
      NGem=MAXVAL(IGem)
      ECorrIJ=0.0d0
      IJ=0
      Do I=1,NGem
         Do J=1,I
            IJ=IJ+1
            IGIJ(I,J)=IJ
            IGIJ(J,I)=IJ
            IGemNo(IJ,1)=I
            IGemNo(IJ,2)=J
        EndDo
      EndDo
C
      IPair(1:NBasis,1:NBasis)=0
      Do II=1,NDimX
      I=IndN(1,II)
      J=IndN(2,II)
      IPair(I,J)=1
      IPair(J,I)=1
      EndDo
C
C     AUXILIARY STUFF LATER NEEDED TO GET A+ AND A- MATRICES FOR ALPHA=0
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
      ETot=ETot+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoNO(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
C
      Write(6,'(/,X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)')ETot
      Do I=1,NBasis
      C(I)=SQRT(Occ(I))
      If(Occ(I).Lt.Half) C(I)=-C(I)
      CICoef(I)=C(I)
      EndDo
C
C     CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA=0 HAMILTONIAN
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Ne.IGem(J)) Then
C
      HNO(IJ)=Zero
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
C     AUXILIARY MATRIX AuxI AND AuxIO
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      IPQ=IPQ+1
      AuxI(IPQ)=Zero
      AuxIO(IPQ)=Zero
      Do IT=1,NOccup
      If(IGFact(NAddr3(IT,IT,IP,IQ)).Eq.1) Then
       AuxI(IPQ)=AuxI(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      If(IT.Le.INActive) AuxIO(IPQ)=AuxIO(IPQ)+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IP,IQ,IT,IT))-TwoNO(NAddr3(IP,IT,IQ,IT)))
      EndIf
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
      If(IGFact(NAddr3(IT,IW,IP,IU)).Eq.1)
     $ WMAT(IP,IR)=WMAT(IP,IR)
     $ +TwoNO(NAddr3(IT,IW,IP,IU))
     $ *FRDM2(IW,IU,IT,IR,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ +TwoNO(NAddr3(IT,IU,IP,IW))
     $ *FRDM2(IW,IU,IR,IT,RDM2Act,Occ,Ind2,NAct,NBasis)
C
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-ACTIVE BLOCK
C
      Write(6,'(" *** ACTIVE-ACTIVE BLOCK ***")')
C
      NFree1=1
      NFree2=1
      NoEig=0
C
      NDimB=0
      Do IQQ=1,NAct
      Do IPP=IQQ+1,NAct
      IP=Ind1(IPP)
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IQQ=1,NAct
      Do IPP=IQQ+1,NAct
      IP=Ind1(IPP)
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      Do ISS=1,NAct
      Do IRR=ISS+1,NAct
      IR=Ind1(IRR)
      IS=Ind1(ISS)
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
C      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
C     ADD A SR KERNEL
C
      If(IFunSRKer.Eq.1) Then
      I1I2S=MultpC(NSymNO(IP),NSymNO(IQ))
      I3I4S=MultpC(NSymNO(IR),NSymNO(IS))
      ISym=MultpC(I1I2S,I3I4S)
      XKer1234=Zero
      If(ISym.Eq.1) Then
      Do I=1,NGrid
      XKer1234=XKer1234+SRKerW(I)*
     $ OrbGrid(I+(IP-1)*NGrid)*OrbGrid(I+(IQ-1)*NGrid)*
     $ OrbGrid(I+(IR-1)*NGrid)*OrbGrid(I+(IS-1)*NGrid)
      EndDo
      EndIf
      TwoSR=TwoEl2(NAddr3(IP,IQ,IR,IS))-TwoNO(NAddr3(IP,IQ,IR,IS))
      ABM=ABM+Four*(C(IP)+C(IQ))*(C(IR)+C(IS))*(XKer1234+TwoSR)
      EndIf
C
C      ABPLUS((ICol-1)*NDimB+IRow)=ABP
C      ABMIN((ICol-1)*NDimB+IRow)=ABM
      ABPLUS((IRow-1)*NDimB+ICol)=ABP
      ABMIN((IRow-1)*NDimB+ICol)=ABM
C
C      EndIf
C
      EndIf
      EndDo
      EndDo
C
      EndIf
      EndDo
      EndDo
C
      If(NDimB.Ne.0) Then
C     SYMMETRIZE BLOCK
      Do I=1,NDimB
      Do J=I+1,NDimB
      ABPLUS((J-1)*NDimB+I)=
     $ Half*(ABPLUS((J-1)*NDimB+I)+ABPLUS((I-1)*NDimB+J))
      ABPLUS((I-1)*NDimB+J)=ABPLUS((J-1)*NDimB+I)
      ABMIN((J-1)*NDimB+I)=
     $ Half*(ABMIN((J-1)*NDimB+I)+ABMIN((I-1)*NDimB+J))
      ABMIN((I-1)*NDimB+J)=ABMIN((J-1)*NDimB+I)
c herer!!!
c      ABMIN((I-1)*NDimB+J)=Zero
c      ABMIN((J-1)*NDimB+I)=Zero
c      ABPLUS((I-1)*NDimB+J)=Zero
c      ABPLUS((J-1)*NDimB+I)=Zero
      EndDo
      EndDo
C
C      Print*, 'ACT-ACT-Ka',NDimB,norm2(ABPLUS(1:NDimB**2)),
C     $ norm2(ABMIN(1:NDimB**2))
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $ NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE ACTIVE-INACTIVE BLOCKS
C
      Write(6,'(" *** ACTIVE-INACTIVE BLOCKS ***")')
C
      Do IQ=1,INActive
C
      NDimB=0
      Do IPP=1,NAct
      IP=Ind1(IPP)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IPP=1,NAct
      IP=Ind1(IPP)
C
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      IS=IQ
      Do IRR=1,NAct
      IR=Ind1(IRR)
C
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
C      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
C      ABPLUS((IRow-1)*NDimB+ICol)=ABP
C      ABMIN((IRow-1)*NDimB+ICol)=ABM
C
C      EndIf
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      If(NDimB.Ne.0) Then
C     SYMMETRIZE BLOCK
      Do I=1,NDimB
      Do J=I+1,NDimB
      ABPLUS((J-1)*NDimB+I)=
     $ Half*(ABPLUS((J-1)*NDimB+I)+ABPLUS((I-1)*NDimB+J))
      ABPLUS((I-1)*NDimB+J)=ABPLUS((J-1)*NDimB+I)
      ABMIN((J-1)*NDimB+I)=
     $ Half*(ABMIN((J-1)*NDimB+I)+ABMIN((I-1)*NDimB+J))
      ABMIN((I-1)*NDimB+J)=ABMIN((J-1)*NDimB+I)
c herer!!!
c      ABMIN((I-1)*NDimB+J)=zero
c      ABMIN((J-1)*NDimB+I)=zero
c      ABPLUS((I-1)*NDimB+J)=Zero
c      ABPLUS((J-1)*NDimB+I)=Zero
      EndDo
      EndDo
C
C      Print*, 'AI-Ka',IP,norm2(ABPLUS(1:NDimB**2)),
C     $ norm2(ABMIN(1:NDimB**2))
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     Do IP
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-ACTIVE BLOCKS
C
      Write(6,'(" *** VIRTUAL-ACTIVE BLOCKS ***")')
C
      Do IP=NOccup+1,NBasis
C
      NDimB=0
      Do IQQ=1,NAct
      IQ=Ind1(IQQ)
      If(IPair(IP,IQ).Eq.1) Then
      NDimB=NDimB+1
      IndBlock(1,NFree1-1+NDimB)=IP
      IndBlock(2,NFree1-1+NDimB)=IQ
      EndIf
      EndDo
C
      Do I=1,NDimB
      IEigAddY(1,NFree1-1+I)=NFree2+(I-1)*NDimB
      IEigAddY(2,NFree1-1+I)=IEigAddY(1,NFree1-1+I)+NDimB-1
      IEigAddInd(1,NFree1-1+I)=NFree1
      IEigAddInd(2,NFree1-1+I)=NFree1+NDimB-1
      EndDo
C
      IRow=0
      Do IQQ=1,NAct
      IQ=Ind1(IQQ)
C
      If(IPair(IP,IQ).Eq.1) Then
C
      IRow=IRow+1
C
      ICol=0
      IR=IP
      Do ISS=1,NAct
      IS=Ind1(ISS)
C
      If(IPair(IR,IS).Eq.1) Then
C
      ICol=ICol+1
C
C      If(IRow.Ge.ICol) Then
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IR,IS,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      ABPLUS((ICol-1)*NDimB+IRow)=ABP
      ABMIN((ICol-1)*NDimB+IRow)=ABM
C      ABPLUS((IRow-1)*NDimB+ICol)=ABP
C      ABMIN((IRow-1)*NDimB+ICol)=ABM
C
C      EndIf
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      If(NDimB.Ne.0) Then
C     SYMMETRIZE BLOCK
      Do I=1,NDimB
      Do J=I+1,NDimB
      ABPLUS((J-1)*NDimB+I)=
     $ Half*(ABPLUS((J-1)*NDimB+I)+ABPLUS((I-1)*NDimB+J))
      ABPLUS((I-1)*NDimB+J)=ABPLUS((J-1)*NDimB+I)
      ABMIN((J-1)*NDimB+I)=
     $ Half*(ABMIN((J-1)*NDimB+I)+ABMIN((I-1)*NDimB+J))
      ABMIN((I-1)*NDimB+J)=ABMIN((J-1)*NDimB+I)
c herer!!!
c      ABMIN((I-1)*NDimB+J)=zero
c      ABMIN((J-1)*NDimB+I)=zero
c      ABPLUS((I-1)*NDimB+J)=Zero
c      ABPLUS((J-1)*NDimB+I)=Zero
      EndDo
      EndDo
C
C      Print*, 'VA-Ka',IP,norm2(ABPLUS(1:NDimB**2)),
C     $ norm2(ABMIN(1:NDimB**2))
      If(NoSt.Eq.1) Then
      Call ERPASYMM0(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $   NDimB)
      Else
      Call ERPAVECYX(EigY(NFree2),EigX(NFree2),Eig(NFree1),ABPLUS,ABMIN,
     $    NDimB)
      EndIf
      EndIf
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
C     Do IP
      EndDo
C
C     FIND THE 0TH-ORDER SOLUTION FOR THE VIRTUAL-INACTIVE BLOCKS
C
      if(idalton.eq.0) then
      open(10,file='fock.dat')
      work1=0
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
      read(10,*)iip,iiq,xx
      work1(iip,iiq)=xx
      EndDo
      EndDo
      close(10)
      endif
C
C     counter for testing if Fock matrix is diagonal...
      icount=0
      Do IP=NOccup+1,NBasis
      Do IQ=1,INActive
C
      NDimB=0
C
      If(IPair(IP,IQ).Eq.1) Then
C
      NDimB=1
      IndBlock(1,NFree1)=IP
      IndBlock(2,NFree1)=IQ
C
      IEigAddY(1,NFree1)=NFree2
      IEigAddY(2,NFree1)=IEigAddY(1,NFree1)
      IEigAddInd(1,NFree1)=NFree1
      IEigAddInd(2,NFree1)=NFree1
C
      Call AB0ELEMENT(ABP,ABM,IP,IQ,IP,IQ,Occ,HNO,IGFact,
     $ TwoNO,AuxI,AuxIO,WMAT,RDM2Act,C,Ind1,Ind2,NAct,NRDM2Act,
     $ NInte1,NInte2,NBasis)
C
      Eig(NFree1)=ABP
C
      If(IDALTON.Eq.0) Then
      If(ABS(ABP-work1(IP,IQ)).Gt.1.d-7) then
      Write(*,*)'ABP inconsistent with eps_a-eps_i for',IP,IQ
      icount=icount+1
      EndIf
      EndIf
C
      EigY(NFree2)=One/Sqrt(Two)
      EigX(NFree2)=One/Sqrt(Two)
C
      NoEig=NoEig+NDimB
      NFree1=NoEig+1
      NFree2=NFree2+NDimB*NDimB
C
      EndIf
C
      EndDo
      EndDo

      if(icount.eq.0) write(6,'("ABP consistent with eps_a-eps_i")')
C
      Write(6,'(/," *** DONE WITH 0TH-ORDER IN AC0-CASSCF ***")')
C
C      Print*, 'Eig,Y,X',norm2(Eig),norm2(EigY),norm2(EigX)
C     DONE 0TH-ORDER CALCULATIONS
C
      Write(6,'(/,
     $" *** COMPUTING ABPLUS(1) AND ABMIN(1) MATRICES ***"
     $ )')
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      If (IDBBSC.Gt.0) Then
C
C     MODIFY integrals with 1/r TO ACCOUNT FOR CBS CORRECTION
C
      Call LOC_MU_CBS(XMuMAT,URe,UNOAO,Occ,TwoNO,NBasis,NInte2)
c      goto 999
C
C     COMPUTE 1st-ORDER ONE-ELECTRON HAMILTONIAN (USING 1/r INTEGRALS)
C
C     SET THE FLAG FOR AB1_CAS
C
      WorkH=XOne
      IHNO1=1
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
C
      If(IGem(I).Eq.IGem(J)) Then
C
      Aux=Zero
C
      Do IT=1,NBasis
      If(IGem(IT).Ne.IGem(I))
     $ Aux=Aux+Occ(IT)*
     $ (Two*TwoNO(NAddr3(IT,IT,I,J))-TwoNO(NAddr3(IT,I,IT,J)))
      EndDo
C
      XOne(IJ)=-Aux
C
      EndIf
C
      EndDo
      EndDo
C
C     COMPUTE SR INTEGRALS
C
      TwoEl2=TwoNO-TwoEl2
      TwoElSR=Zero
C
C     ADD MODIFIED SR INTEGRALS TO FULL-COULOMB INTEGRALS
C
      NAddr=0
C
      IPR=0
      Do IP=1,NBasis
      Do IR=1,IP

      IPR=IPR+1
C
      IQS=0
      Do IQ=1,NBasis
      Do IS=1,IQ

      IQS=IQS+1
C
      If(IPR.Ge.IQS) Then
C
      NAddr=NAddr+1
      AuxSR=Zero
      Do I=1,NBasis
      AuxSR=AuxSR
     $ +TwoEl2(NAddr3(I,IR,IQ,IS))*XMuMAT(IP,I)
     $ +TwoEl2(NAddr3(IP,I,IQ,IS))*XMuMAT(IR,I)
     $ +TwoEl2(NAddr3(IP,IR,I,IS))*XMuMAT(IQ,I)
     $ +TwoEl2(NAddr3(IP,IR,IQ,I))*XMuMAT(IS,I)
      EndDo
C
      AuxSR=AuxSR/Four
C
      TwoElSR(NAddr)=AuxSR
      TwoNO(NAddr)=TwoNO(NAddr)+AuxSR

c      If(.not.(IGem(IP).eq.1.and.
c     $ IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
c     $ And.IGem(IP).Eq.IGem(IQ))) then
c       TwoElSR(NAddr)=AuxSR
c       TwoNO(NAddr)=TwoNO(NAddr)+AuxSR
c       endif
C
C
c comment out the line above and uncomment the lines below to
c include in-active contribution to the CBS correction
c which will be printed as "EIntra"; in practice it is always very small
c      If(IGFact(NAddr3(IP,IR,IQ,IS)).Eq.0) Then
c      TwoNO(NAddr)=TwoNO(NAddr)+AuxSR
c      Else
c      TwoNO(NAddr)=AuxSR
c      EndIf
C
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
c uncomment to include in-active contribution to the CBS correction (see above)
c      IGFact=0
c comment out to include in-active contribution to the CBS correction (see above)
      IHNO1=0
C
      Ind3=0
      Do I=1,NStronglyOccOrb
      Ind3(I)=1
      EndDo
      ETotSR=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
C
      If(.Not.(Ind3(IP).Eq.1.And.
     $ Ind3(IP).Eq.Ind3(IR).And.Ind3(IQ).Eq.Ind3(IS).
     $ And.Ind3(IP).Eq.Ind3(IQ))) Then

      If(IGem(IP).Eq.2.And.
     $ IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ) )  
     $ ETotSR=ETotSR+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *TwoElSR(NAddr3(IP,IR,IQ,IS))

      If(NInAcCAS.Eq.0) Then
      If(IGem(IP).Eq.1.And.
     $ IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ) )
     $ ETotSR=ETotSR+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $  *TwoElSR(NAddr3(IP,IR,IQ,IS))
      EndIf
  
      EndIf

c      If(.Not.(Ind3(IP).Eq.1.And.
c     $ Ind3(IP).Eq.Ind3(IR).And.Ind3(IQ).Eq.Ind3(IS).
c     $ And.Ind3(IP).Eq.Ind3(IQ))) then
c      if(FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis).ne.zero) 
c     $ write(*,*)ip,iq,ir,is
c      endif
 
      EndDo
      EndDo
      EndDo
      EndDo
      Write(6,'(/,X,"1st-order SR energy for",I4," active orbitals")') 
     $ NAcCAS+NInAcCAS-NStronglyOccOrb
      Write(6,'(X," <Ref|H^SR|Ref> = ",F10.6,/)') ETotSR

      goto 444
C    ***** Two-Electron systems investigation
      NDim=NBasis*(NBasis-1)/2
C     AuxI=VSR
      AuxI=Zero 
      NGem=MAXVAL(IGem)
      NoEig=1
      TwoNO=TwoNO-TwoElSR
C
      NAddr=0
      IPR=0
      Do IP=1,NBasis
      Do IR=1,IP
      IPR=IPR+1
      IQS=0
      Do IQ=1,NBasis
      Do IS=1,IQ
      IQS=IQS+1
      If(IPR.Ge.IQS) Then
      NAddr=NAddr+1
C
      If((IGem(IP).Eq.1.And.
     $ IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))) Then 
      If((Ind3(IP).Eq.1.And.
     $ Ind3(IP).Eq.Ind3(IR).And.Ind3(IQ).Eq.Ind3(IS).
     $ And.Ind3(IP).Eq.Ind3(IQ))) TwoElSR(NAddr)=Zero
      EndIf
C
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C Set TwoElSR to Zero to reproduce FCI[H] energy
c      TwoElSR=Zero 
C
C     ETot=<Psi_ref|H0+H'|Psi_ref> so a contribution from the SR interaction must be added
C
      ETotSR=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      ETotSR=ETotSR+FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *TwoElSR(NAddr3(IP,IR,IQ,IS))
      EndDo
      EndDo
      EndDo
      EndDo
      ERef=ETot
      ETot=ETot+ETotSR
      Call ACPINORS(ETot,ENuc,TwoNO,TwoElSR,Occ,WorkH,AuxI,
     $ NBasis,NInte1,NInte2,NDim,NGem,NoEig)
      Write(6,'(/,X," Where, <Ref|H0|Ref> + ENuc = ",F10.6,
     $" <Ref|H^SR|Ref> = ",F10.6)') ERef+ENuc,ETotSR 
      Write(6,'(X," H^SR for ",I2," active orbitals")')
     $ NAcCAS+NInAcCAS-NStronglyOccOrb
c      Call PairD_Vis(URe,UNOAO,Occ,NBasis)
      Stop 
  444 continue
C end of  ***** Two-Electron systems investigation
C
      EndIf
  999 continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      Call AB1_CAS(ABPLUS,ABMIN,URe,Occ,XOne,TwoNO,
     $ RDM2Act,NRDM2Act,IGFact,C,Ind1,Ind2,
     $ IndBlock,NoEig,NDimX,NBasis,NInte1,NInte2)
       Print*, 'AB1-KA',norm2(ABPLUS),norm2(ABMIN)
C
C     ADD A SR KERNEL
C
      If(IFunSRKer.Eq.1) Then
C
      Write(6,'(/," *** ADDING THE SR KERNEL ***" )')
CC
C      Do IRow=1,NoEig
CC
C      IA=IndBlock(1,IRow)
C      IB=IndBlock(2,IRow)
C      CA=CICoef(IA)
C      CB=CICoef(IB)
CC
C      Do ICol=1,NoEig
CC
C      IC=IndBlock(1,ICol)
C      ID=IndBlock(2,ICol)
C      CC=CICoef(IC)
C      CD=CICoef(ID)
CC
C      XKer1234=Zero
CC
C      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
C      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
C      ISym=MultpC(I1I2S,I3I4S)
CC
C      If((ISym.Eq.1).And.(IGFact(NAddr3(IA,IB,IC,ID)).Eq.0)) Then
C      Do I=1,NGrid
C      XKer1234=XKer1234+SRKerW(I)*
C     $ OrbGrid(I+(IA-1)*NGrid)*OrbGrid(I+(IB-1)*NGrid)*
C     $ OrbGrid(I+(IC-1)*NGrid)*OrbGrid(I+(ID-1)*NGrid)
C
C      EndDo
C      EndIf
CC
C      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
CC
C      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
C     $ +Four*(CA+CB)*(CC+CD)*(XKer1234+TwoSR)
CC
C      EndDo
C      EndDo
CC
C     ALTERNATIVELY:
C     1st: ADD INTEGRALS
C
      Do IRow=1,NoEig
C
      IA=IndBlock(1,IRow)
      IB=IndBlock(2,IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NoEig
C
      IC=IndBlock(1,ICol)
      ID=IndBlock(2,ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      TwoSR=TwoEl2(NAddr3(IA,IB,IC,ID))-TwoNO(NAddr3(IA,IB,IC,ID))
C
      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
     $ +Four*(CA+CB)*(CD+CC)*TwoSR
C
      EndDo
      EndDo
C
C     2nd: ADD KERNEL IN BATCHES:
C
      call gclock('START',Tcpu,Twall)
      Allocate(Work(Maxlen),Batch(Maxlen,NBasis))
C
      Do Offset=0,NGrid,Maxlen
      Batchlen=min(NGrid-Offset,Maxlen)
      If(Batchlen==0) exit
C
      Work(1:Batchlen) = SRKerW(Offset+1:Offset+Batchlen)
      Call FILL_BATCH(OrbGrid,Batch,Maxlen,Batchlen,Offset,NGrid,NBasis)
C
      Do IRow=1,NoEig
C
      IA=IndBlock(1,IRow)
      IB=IndBlock(2,IRow)
      CA=CICoef(IA)
      CB=CICoef(IB)
C
      Do ICol=1,NoEig
C
      IC=IndBlock(1,ICol)
      ID=IndBlock(2,ICol)
      CC=CICoef(IC)
      CD=CICoef(ID)
C
      XKer1234=Zero
      I1I2S=MultpC(NSymNO(IA),NSymNO(IB))
      I3I4S=MultpC(NSymNO(IC),NSymNO(ID))
      ISym=MultpC(I1I2S,I3I4S)
C
      If((ISym.Eq.1).And.(IGFact(NAddr3(IA,IB,IC,ID)).Eq.0)) Then
      Do I=1,Batchlen
      XKer1234=XKer1234+Work(I)*
     $ Batch(I,IA)*Batch(I,IB)*Batch(I,IC)*Batch(I,ID)
      EndDo
      EndIf
C
      ABMIN((ICol-1)*NoEig+IRow)=ABMIN((ICol-1)*NoEig+IRow)
     $ +Four*(CA+CB)*(CD+CC)*XKer1234
C
      EndDo
      EndDo
      EndDo
C
      Deallocate(Batch,Work)
CC
      Write(6,'("*** sr-kernel added ***")')
      call gclock('sr-kernel AB(1)',Tcpu,Twall)
CC
      EndIf
C
      Deallocate(RDM2Act)
C
C      Print*, 'ABM-KA',norm2(ABMIN(1:NoEig**2))
C
      Write(6,'(/," *** DONE WITH COMPUTING AB(1) MATRICES ***")')
C
C     1ST-ORDER PART
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      XMAux(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      XMAux(NU+(MU-1)*NoEig)=XMAux(NU+(MU-1)*NoEig)
     $ +ABPLUS(NU+(I-1)*NoEig)*EigX(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      ABPLUS(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
      ABPLUS(NU+(MU-1)*NoEig)=ABPLUS(NU+(MU-1)*NoEig)
     $ +EigX(IStart+II)*XMAux(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      XMAux(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
      XMAux(NU+(MU-1)*NoEig)=XMAux(NU+(MU-1)*NoEig)
     $ +ABMIN(NU+(I-1)*NoEig)*EigY(IStart+II)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      Do NU=1,NoEig
      Do MU=1,NoEig
C
      ABMIN(NU+(MU-1)*NoEig)=Zero
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
      ABMIN(NU+(MU-1)*NoEig)=ABMIN(NU+(MU-1)*NoEig)
     $ +EigY(IStart+II)*XMAux(I+(MU-1)*NoEig)
      II=II+1
      EndDo
C
      EndDo
      EndDo
C
      XMAux(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C KP 03/04/2018
C
      If(Eig(MU).Ne.Zero) Then
C
      Do NU=1,NoEig
C
C KP 03/04/2018
C
      If(Eig(NU).Ne.Zero) Then
C
      IStart=IEigAddY(1,NU)
      II=0
      Do I=IEigAddInd(1,NU),IEigAddInd(2,NU)
C
      XMAux(MU+(I-1)*NoEig)=XMAux(MU+(I-1)*NoEig)+Two*
     $ (ABPLUS(MU+(NU-1)*NoEig)-ABMIN(MU+(NU-1)*NoEig))/
     $ (Eig(MU)+Eig(NU))*EigY(IStart+II)
C
      II=II+1
      EndDo
C
      EndIf
C
      EndDo
C
      EndIf
C
      EndDo
C
      ABPLUS(1:NoEig*NoEig)=Zero
C
      Do MU=1,NoEig
C
      IStart=IEigAddY(1,MU)
      II=0
      Do I=IEigAddInd(1,MU),IEigAddInd(2,MU)
C
      Do J=1,NoEig
      ABPLUS(I+(J-1)*NoEig)=ABPLUS(I+(J-1)*NoEig)+XMAux(MU+(J-1)*NoEig)
     $ *EigY(IStart+II)
      EndDo
C
      II=II+1
      EndDo
      EndDo
C
C     FINALLY THE ENERGY CORRECTION
C
      Print*, 'ABPLUS-energy corr=',norm2(ABPLUS)
C
      EAll=Zero
      EIntra=Zero
C
      Do I=1,NoEig
C
      IP=IndBlock(1,I)
      IR=IndBlock(2,I)
C
      Do J=1,NoEig
C
      IQ=IndBlock(1,J)
      IS=IndBlock(2,J)
C
      If(IP.Gt.IR.And.IQ.Gt.IS) Then
C
      SumY=ABPLUS(I+(J-1)*NoEig)
      Aux=(C(IS)+C(IQ))*(C(IP)+C(IR))*SumY
C
      EAll=EAll+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      If(IGem(IP).Eq.IGem(IR).And.IGem(IQ).Eq.IGem(IS).
     $ And.IGem(IP).Eq.IGem(IQ))
     $ EIntra=EIntra+Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
      ECorrIJ(IGIJ(IGem(IP),IGem(IR)),IGIJ(IGem(IQ),IGem(IS)))= 
     $ ECorrIJ(IGIJ(IGem(IP),IGem(IR)),IGIJ(IGem(IQ),IGem(IS)))
     $              +Aux*TwoNO(NAddr3(IP,IR,IQ,IS))
C
C     endinf of If(IP.Gt.IR.And.IQ.Gt.IS)
      EndIf
C
      EndDo
      EndDo
C
      Print*, "EAll EIntra", EAll,EIntra
C
      ECorr=EAll-EIntra
C
!PRINT CONTRIBUTIONS TO ECorr FROM BLOCKS
      Write(6,'(/,X,
     $ "Contributions to AC0 from (Mu)(Nu) pairs of blocks")')
      IJ=0
      Do I=1,NGem
      Do J=1,I
      IJ=IJ+1
! NGem=3 case
       If(NGem.Eq.3) Then
       IOO=0
       If(IGemNo(IJ,1).Eq.1.And.IGemNo(IJ,2).Eq.1) IOO=1
       IVV=0
       If(IGemNo(IJ,1).Eq.3.And.IGemNo(IJ,2).Eq.3) IVV=1
       IAA1=0
       If(IGemNo(IJ,1).Eq.2.And.IGemNo(IJ,2).Eq.2) IAA1=1
! NGem=2 case (zero inactive orbitals)
       ElseIf(NGem.Eq.2) Then
       IOO=0
       IVV=0
       If(IGemNo(IJ,1).Eq.2.And.IGemNo(IJ,2).Eq.2) IVV=1
       IAA1=0
       If(IGemNo(IJ,1).Eq.1.And.IGemNo(IJ,2).Eq.1) IAA1=1
       EndIf

       If(IVV==0.And.IOO==0) Then
       KL=0
       Do K=1,NGem
       Do L=1,K
         KL=KL+1
         If(NGem.Eq.3) Then
             IOO=0
             If(IGemNo(KL,1).Eq.1.And.IGemNo(KL,2).Eq.1) IOO=1
             IVV=0
             If(IGemNo(KL,1).Eq.3.And.IGemNo(KL,2).Eq.3) IVV=1
             IAA2=0
             If(IGemNo(KL,1).Eq.2.And.IGemNo(KL,2).Eq.2) IAA2=1
         ElseIf(NGem.Eq.2) Then
             IOO=0
             IVV=0
             If(IGemNo(KL,1).Eq.2.And.IGemNo(KL,2).Eq.2) IVV=1
             IAA2=0
             If(IGemNo(KL,1).Eq.1.And.IGemNo(KL,2).Eq.1) IAA2=1
         EndIf
         If(IVV==0.And.IOO==0.And.IAA1+IAA2.Ne.2) Then
         If(IJ.Ge.KL) Then
           EE=ECorrIJ(IJ,KL)
           If(IJ.Ne.KL)EE=EE+ECorrIJ(KL,IJ)
           Write(6,'(X,"(",2I1,")","(",2I1,")",F15.8)')
     $           IGemNo(IJ,1),IGemNo(IJ,2),IGemNo(KL,1),IGemNo(KL,2),EE
         EndIf
         EndIf
      EndDo
      EndDo
      EndIf
      EndDo
      EndDo
C
      Return
      End

*Deck fil_Batch
      Subroutine FILL_BATCH(OrbGrid,Batch,Maxlen,Batchlen,
     $                      Offset,NGrid,NBasis)
C
C     FILLS BATCHES FOR GRID CALCULATIONS
C
      Implicit None
C
      Integer :: Maxlen,Batchlen,Offset,NGrid,NBasis
      Double Precision :: Batch(Maxlen,NBasis),OrbGrid(NGrid,NBasis)
C
      Batch(1:Batchlen,1:NBasis) =
     $ OrbGrid(Offset+1:Offset+Batchlen,1:NBasis)
C
      End Subroutine FILL_BATCH

*Deck GGA_ONTOP
      Subroutine GGA_ONTOP(EXCTOP,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A GGA_XC ENERGY IF THE DENSITY AND SPIN-DENSITY ARE COMPUTED AS:
C     RHO_A/B = 1/2 ( RHO +/- SQRT(RHO^2 - 2 PI )  )
C     WHERE PI IS THE ON-TOP PAIR DENSITY, see Gagliardi JCP 146, 034101 (2017)
C
C     XCFUN IS USED !!!
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
      Dimension Zk(NGrid),RhoA(NGrid),RhoB(NGrid),
     & SigmaAA(NGrid),SigmaAB(NGrid),SigmaBB(NGrid)
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
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid,Occ,URe,OrbGrid,NGrid,NBasis)
C
      If(RhoGrid.Gt.1.D-12) Then
C
      OnTop=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop=OnTop
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      Rho=RhoGrid
      R=Two*OnTop/Rho**2
      XFactor=Zero
      If(R.Lt.One) XFactor=SQRT(One-R)
C
      RhoA(I)=Rho/Two*(One+XFactor)
      RhoB(I)=Rho/Two*(One-XFactor)
C
      RhoXa=RhoX/Two*(One+XFactor)
      RhoXb=RhoX/Two*(One-XFactor)
      RhoYa=RhoY/Two*(One+XFactor)
      RhoYb=RhoY/Two*(One-XFactor)
      RhoZa=RhoZ/Two*(One+XFactor)
      RhoZb=RhoZ/Two*(One-XFactor)
C
      SigmaAA(I)=RhoXa*RhoXa+RhoYa*RhoYa+RhoZa*RhoZa
      SigmaAB(I)=RhoXa*RhoXb+RhoYa*RhoYb+RhoZa*RhoZb
      SigmaBB(I)=RhoXb*RhoXb+RhoYb*RhoYb+RhoZb*RhoZb
C
      Else
      RhoA(I)=Zero
      RhoB(I)=Zero
      SigmaAA(I)=Zero
      SigmaAB(I)=Zero
      SigmaBB(I)=Zero
C
      EndIf
C
      EndDo
C
      Call dfun_GGA_AB(RhoA,RhoB,SigmaAA,SigmaAB,SigmaBB,Zk,NGrid)
C
      EXCTOP=Zero
      Do I=1,NGrid
      EXCTOP=EXCTOP+Zk(I)*WGrid(I)
      EndDo
C
      Return
      End

*Deck SR_PBE_ONTOP
      Subroutine SR_PBE_ONTOP(EXCTOP,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A SR-PBE XC ENERGY IF THE DENSITY AND SPIN-DENSITY ARE COMPUTED AS:
C     RHO_A/B = 1/2 ( RHO +/- SQRT(RHO^2 - 2 PI )  )
C     WHERE PI IS THE ON-TOP PAIR DENSITY, see Gagliardi JCP 146, 034101 (2017)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C     $ WGrid(NGrid),OrbGrid(NBasis,NGrid),OrbXGrid(NBasis,NGrid),
C     $ OrbYGrid(NBasis,NGrid),OrbZGrid(NBasis,NGrid)
C
! input
      Dimension Zk(NGrid),RhoGrid(NGrid),RhoO(NGrid),
     & Sigma(NGrid),SigmaCO(NGrid),SigmaOO(NGrid)
      logical fderiv,open
! output
      integer igrad
      character*(30) name
      double precision vrhoc(ngrid),vrhoo(ngrid)
      double precision vsigmacc(ngrid),vsigmaco(ngrid),vsigmaoo(ngrid)
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
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
C
      If(RhoGrid(I).Gt.1.D-12) Then
C
      OnTop=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop=OnTop
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
C     $ *OrbGrid(IP,I)*OrbGrid(IQ,I)*OrbGrid(IR,I)*OrbGrid(IS,I)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
C
      Rho=RhoGrid(I)
      R=Two*OnTop/Rho**2
      XFactor=Zero
      If(R.Lt.One) XFactor=SQRT(One-R)
C
      Rhoa=Rho/Two*(One+XFactor)
      Rhob=Rho/Two*(One-XFactor)
      RhoGrid(I)=Rhoa+Rhob
      RhoO(I)=Rhoa-Rhob
C
      RhoXa=RhoX/Two*(One+XFactor)
      RhoXb=RhoX/Two*(One-XFactor)
      RhoYa=RhoY/Two*(One+XFactor)
      RhoYb=RhoY/Two*(One-XFactor)
      RhoZa=RhoZ/Two*(One+XFactor)
      RhoZb=RhoZ/Two*(One-XFactor)
C
      RhoXC=RhoXa+RhoXb
      RhoYC=RhoYa+RhoYb
      RhoZC=RhoZa+RhoZb

      RhoXO=RhoXa-RhoXb
      RhoYO=RhoYa-RhoYb
      RhoZO=RhoZa-RhoZb
C
      Sigma(I)=  RhoXC*RhoXC+RhoYC*RhoYC+RhoZC*RhoZC
      SigmaCO(I)=RhoXC*RhoXO+RhoYC*RhoYO+RhoZC*RhoZO
      SigmaOO(I)=RhoXO*RhoXO+RhoYO*RhoYO+RhoZO*RhoZO
C
      Else
      RhoGrid(I)=Zero
      RhoO(I)=Zero
      Sigma(I)=Zero
      SigmaCO(I)=Zero
      SigmaOO(I)=Zero
C
      EndIf
C
      EndDo
C
      EXCTOP=Zero
      FDeriv=.False.
      Open=.True.
C
      Call dftfun_exerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,RhoO,
     >                   Sigma,SigmaCO,SigmaOO,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)

      Do I=1,NGrid
      EXCTOP=EXCTOP+Zk(I)*WGrid(I)
      EndDo
      Exch=EXCTOP
C
      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,RhoO,
     >                   Sigma,SigmaCO,SigmaOO,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
      EnC=Zero
      Do I=1,NGrid
      EnC=EnC+Zk(I)*WGrid(I)
      EndDo
C
      EXCTOP=EXCTOP+EnC
C
      Write(6,'(/," SR_xch_PBE with translated densities",F15.8)')Exch
      Write(6,'(" SR_cor_PBE with translated densities",F15.8)')EnC
C
      Return
      End

*Deck RunDFOnTop
      Subroutine RunDFOnTop(ETot,ENuc,TwoNO,URe,UNOAO,Occ,XOne,
     $  IndAux,IPair,IndN,IndX,NDimX,Title,NBasis,NInte1,NInte2,NGem)
C
C     ETot is calculated from MC-PDFT
C     with PBE xc functional
C     see Eq.(6) in Manni, et al. JCTC 10, 3669-3680 (2014)
C     doi: 10.1021/ct500483t
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 FMultTab,Title
      Include 'commons.inc'
C
      Character*60 FName
      Real*8, Dimension(:), Allocatable :: OrbGrid
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RDM2Act
      Dimension NSymNO(NBasis),VSR(NInte1),MultpC(15,15),NumOSym(15)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Four=4.D0)
C
      Dimension
     $ URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis),
     $ TwoNO(NInte2),XOne(NInte1),IndX(NDimX),IndN(2,NDimX),
     $ IndAux(NBasis),IPair(NBasis,NBasis),Ind1(NBasis),Ind2(NBasis)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),
     $ EigVecR(NDimX*NDimX),Eig(NDimX),
     $ ECorrG(NGem), EGOne(NGem)
C
C     READ 2RDM, COMPUTE THE ENERGY
C
      Write(6,'(/,X," The number of CASSCF Active Orbitals = ",I4)')
     $ NAcCAS
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
C     COMPUTE THE ENERGY FOR CHECKING
C
      EOne=Zero
      Do I=1,NBasis
      II=(I*(I+1))/2
      EOne=EOne+Two*Occ(I)*XOne(II)
      EndDo
      Write(6,'(/,1X,''One-Electron Energy'',6X,F15.8)')EOne
C
C     ITwoEl
      If(ITwoEl.Eq.1) Then
C
      ETot=EOne
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
      ElseIf(ITwoEl.Eq.3) Then
C
      Call TwoEneChck(ETot,RDM2Act,Occ,INActive,NAct,NBasis)
C
C     ITwoEl
      EndIf
C
      Write(6,'(1X,''CASSCF Energy (w/o ENuc)'',X,F15.8)')ETot
      Write(6,'(1X,''Total CASSCF Energy '',5X,F15.8)')ETot+ENuc
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/," The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NBasis*NGrid))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
C
C     load orbgrid and gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
C
      Call GGA_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
     $ OrbZGrid,WGrid,NGrid,NBasis)
      Write(6,'(/," PBE_xc from xcfun with translated densities",
     $ F15.8,/)') EXCTOP
C
c herer!!
c      Alpha=1.D-12
c      Call SR_PBE_ONTOP(EXCTOP,URe,Occ,OrbGrid,OrbXGrid,OrbYGrid,
c     $ OrbZGrid,WGrid,NGrid,NBasis)
c      Write(6,'(/," PBE_xc with translated densities",F15.8,/)')
c     $ EXCTOP
C
C     ITwoEl
      If(ITwoEl.Eq.1) Then
      EnH=Zero
      Do I=1,NOccup
      Do J=1,NOccup
      EnH=EnH+Two*Occ(I)*Occ(J)*TwoNO(NAddr3(I,I,J,J))
      EndDo
      EndDo
C
      ElseIf(ITwoEl.Eq.3) Then
C
      Call TwoEHartree(EnH,RDM2Act,Occ,INActive,NAct,NBasis)
C
      EndIf
C
      ETot=EOne+EnH+EXCTOP
      Write(6,'(1X,''E_Hartree                     '',X,F15.8)') EnH
      Write(6,'(1X,''EOne+EH+E_xc[OnTop] (w/o ENuc)'',X,F15.8)') ETot
      Write(6,'(1X,''EOne+EH+E_xc[OnTop] +ENuc'',5X,F15.8)')ETot+ENuc
C
      DeAllocate (RDM2Act)
      DeAllocate  (WGrid)
      DeAllocate  (OrbGrid)
      DeAllocate  (OrbXGrid)
      DeAllocate  (OrbYGrid)
      DeAllocate  (OrbZGrid)
C
      Return
      End

*Deck CASPI_SR_PBE
      Subroutine CASPI_SR_PBE(URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     COMPUTES A CASPIDFT VERSION OF SR-PBE_CORR ENERGY AS:
C     P(X) eps_corr^SR_PBE
C     X(r) = 2 OnTop(r,r) / Rho^2(r)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis),Ind2(NBasis)
C
! input
      Dimension Zk(NGrid),RhoGrid(NGrid),Sigma(NGrid),OnTop(NGrid)
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
C     READ 2RDM
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
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
      GoTo 10
   40 Continue
      Close(10)
C
      Do I=1,NGrid
C
      vrhoc(I)=Zero
      vsigmacc(i)=Zero
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop(I)=OnTop(I)
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      EndDo
C
      Call dftfun_ecerfpbe(name,FDeriv,Open,igrad,NGrid,RhoGrid,rhoo,
     >                   Sigma,sigmaco,sigmaoo,
     >                   Zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,Alpha)
C
      A=0.2D0
      B=A-One
      C=2.2D0
      G=1.5D0
      D=(C-One)/(One-G)**2
C
      EDYN=Zero
      SR_PBE_c=Zero
      Do I=1,NGrid
C
      If(RhoGrid(I).Ne.Zero) Then
C
      XX=Two*OnTop(I)/RhoGrid(I)**2
      If(XX.Le.One) Then
      PX=A*XX/(One+B*XX)
      Else
      PX=C*XX**0.25-D*(XX-G)**2
      EndIf
C
      EDYN=EDYN+PX*Zk(I)*WGrid(I)
      SR_PBE_c=SR_PBE_c+Zk(I)*WGrid(I)
C
      EndIf
C
      EndDo
C
      Write(6,'(/,1X,''SR-PBE Correlation'',7X,F15.8)')SR_PBE_c
      Write(6,'(1X,''CASPIDFT Correlation'',5X,F15.8,/)')EDYN
C
      Return
      End

*Deck LOC_MU_CHOL_v2
      Subroutine LOC_MU_CHOL(BasisSet,CorrMD,AvMU,
     $                          URe,UNOAO,Occ,NBasis)
C
      use timing
      use types
      use grid_internal
C
C     DGEMM version
C     CAREFUL : memory requirements could still be improved...
C
C     RETURNS SHORT-RANGE CORRELATION ENERGY WITH SR-PBE ONTOP CORRELATION FUNCTIONAL AND LOCAL MU, Giner et al. JCP 152, 174104 (2020)
C
C     IFlFCorr == 0 : calculate CorrMD based only on fCAS function (nact^2 NGrid NCholesky)
C     IFlFCorr == 1 : calculate CorrMD based only on fCAS + fAC0 functions (NBasis^2 NGrid NCholesky)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Half=0.5D0)
c     Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
c    $ Four=4.D0)
C
      Character(*) :: BasisSet
      Include 'commons.inc'
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis)
C
      Real*8, Allocatable :: FPsiB(:),OnTop(:),XMuLoc(:),
     $ Zk(:),RhoGrid(:),Sigma(:),WGrid(:),WorkD(:,:),
     $ Oa(:,:),Oi(:,:),Oia(:,:),Oai(:,:),Q(:,:)
      Real*8,allocatable,target :: PhiGGA(:,:,:)
      Real*8,pointer :: OrbGrid(:,:)
      Real*8,pointer :: OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:)
      Double Precision, Allocatable :: Oaa(:,:),Oii(:,:)
      Double Precision, Allocatable :: OOai(:,:),OOia(:,:)
      Double Precision, Allocatable :: CHVCSAA(:,:),CHVCSII(:,:)
      Double Precision, Allocatable :: CHVCSAI(:,:),CHVCSIA(:,:)
      Double Precision, Allocatable :: RDM2Act(:)
      Double Precision, Allocatable :: RDM2val(:,:,:,:)
      Double Precision, Allocatable :: Work(:,:)
      Dimension IAct(NBasis),Ind2(NBasis)
C     correlation contribution
      Real*8, Allocatable :: FCorr(:)
      Double Precision, Allocatable :: RDM2corr(:,:,:,:)
      Real*8, Allocatable :: QGamII(:,:,:),QGamAA(:,:,:),QGamVV(:,:,:)
      Real*8, Allocatable :: QGamIV(:,:,:),QGamAV(:,:,:),QGamIA(:,:,:)
      Real*8, Allocatable :: OGamI(:,:,:),OGamA(:,:,:),OGamV(:,:,:),
     $                       OnTopCorr(:),XMuLocCorr(:)
      double precision :: Twall,TCpu
C
      Logical :: doGGA
      Logical :: ex ! check rdmcorr file
      Double Precision, Allocatable :: RR(:,:)
      Character(*),Parameter :: griddalfile='dftgrid.dat'
C
      call gclock('START',Tcpu,Twall)
C
C     debug printlevel
      IIPRINT=5
C
C     set constants
C
c     Pi = dacos(-1.d0)
c      SPi = SQRT(Pi)
c     used in Kasia's code:
      SPi=SQRT(3.141592653589793)
      Prefac = SQRT(3.1415)/Two ! SPi/Two
C
      Write(6,'(/,1X,"************** ",
     $ " DGEMM Short-Range Correlation Energy with the Local Mu ")')
      If (IFlFCorr==1) Write(6,'(1X,"Compute fCAS and fAC0!")')
C
C     Load grid (Molpro, Dalton, Internal)
C
      If (InternalGrid==0) then
C
         If (IMOLPRO == 1) Then
            Write(6,'(/,1X,"Read MOLPRO DFT grid file")')
            Call molprogrid0(NGrid,NBasis)
         ElseIf(IDALTON == 1) Then
            Write(6,'(/,1X,"Read DALTON DFT grid file")')
            Call daltongrid0(NGrid,griddalfile,doGGA,NBasis)
            If (.not.doGGA) stop "LOC_MU_CHOL: Run Dalton with SR-PBE!"
         EndIf
C
         Write(6,'(1X,"The number of Grid Points = ",I8)')
     $         NGrid
         Allocate(WGrid(NGrid))
         Allocate(PhiGGA(NGrid,NBasis,4))

      ElseIf (InternalGrid==1) then
C
C     for Orca UNOAO orbitals are read from file
C
!      print*, 'UAONO =',norm2(UNOAO)
!         If (IOrbOrder == 2) then ! orca
!            Allocate(Work(NBasis,NBasis),Oaa(NBasis,NBasis))
!            Oaa=UNOAO
!            Open(newunit=ione,file='uaono.dat',access='sequential',
!     $         form='unformatted',status='old')
!            Read(ione) Work
!            Close(ione)
!            UNOAO = transpose(Work)
!         Endif
         Call internal_gga_no_orbgrid(IGridType,BasisSet,
     $                          IOrbOrder,transpose(UNOAO),
     $                          WGrid,PhiGGA,NGrid,NBasis,NBasis,IUnits)
c
c     for Orca, restore original UNOAO==C(MO,NO)
!     If (IOrbOrder == 2) UNOAO = Oaa
!     If (IOrbOrder == 2) Deallocate(Oaa)
C
      EndIf ! InternalGrid
C
C     use pointers for orbitals and gradients
      OrbGrid  => PhiGGA(:,:,1)
      OrbXGrid => PhiGGA(:,:,2)
      OrbYGrid => PhiGGA(:,:,3)
      OrbZGrid => PhiGGA(:,:,4)
C
      If(InternalGrid==0) Then
C        load orbitals, gradients, and wgrid
         If (IMOLPRO == 1) Then
            Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                      WGrid,UNOAO,NGrid,NBasis)
         ElseIf (IDALTON == 1) Then
            Allocate(Work(NBasis,NBasis))
            Work = transpose(UNOAO)
            Call daltongrid_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                          WGrid,NGrid,griddalfile,NBasis)
            Call daltongrid_tran_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                          Work,0,NGrid,NBasis,.false.)
CC        visualisation
C           Allocate(RR(3,NGrid))
C           Call dalton_grid_coord(NGrid,RR,griddalfile)
           Deallocate(Work)
         EndIf
      EndIf ! InternalGrid
C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=Zero
C
      Open(10,File="rdm2.dat",Status='Old')
C
   10 Read(10,*,End=40)I,J,K,L,X
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=Half*X
      GoTo 10
   40 Continue
      Close(10)
C
      Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
      Do L=1,NOccup
      Do K=1,NOccup
      Do J=1,NOccup
      Do I=1,NOccup
      RDM2val(I,K,J,L)=FRDM2(I,K,J,L,RDM2Act,Occ,Ind2,NAct,NBasis)
      Enddo
      Enddo
      Enddo
      Enddo
      If (IIPRINT.GT.1) Print*, 'RDM2val = ', norm2(RDM2val)
C
      Do I=1,NBasis
      IAct(I)=0
      If(Occ(I).Ne.One.And.Occ(I).Ne.Zero) IAct(I)=1
      EndDo
C
C     LOAD CHOLESKY VECTORS
C
      open(newunit=iunit,file='cholvecs',form='unformatted')
      read(iunit) NCholesky
      Allocate(WorkD(NCholesky,NBasis**2))
      read(iunit) WorkD
      close(iunit)
c     print*,'NCholesky',NCholesky
C
C     truncate CHOLVECS
C
      Allocate (CHVCSAA(NCholesky,NAct**2))
      Do K=1,NCholesky
         IIJJ=0
         Do J=INActive+1,NOccup
            Do I=INActive+1,NOccup
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSAA(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
      Allocate (CHVCSII(NCholesky,INActive**2))
      Do K=1,NCholesky
         IIJJ=0
         Do J=1,INActive
            Do I=1,INActive
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSII(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
      Allocate (CHVCSIA(NCholesky,INActive*NAct))
      Do K=1,NCholesky
         IIJJ=0
         Do J=INActive+1,NOccup
            Do I=1,INActive
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSIA(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
      Allocate (CHVCSAI(NCholesky,NAct*INactive))
      Do K=1,NCholesky
         IIJJ=0
         Do J=1,INActive
            Do I=INActive+1,NOccup
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSAI(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
C
C     COMPUTE f(r)
C
      Allocate (FPsiB(NGrid))
      Allocate (Q(NAct,NAct))
      Allocate (Oa(NCholesky,NAct))
      Allocate (Oai(NCholesky,NAct))
      Allocate (Oi(NCholesky,INActive))
      Allocate (Oia(NCholesky,INActive))
c
      Allocate (Oaa(NAct,NAct))
      Allocate (Oii(INActive,INActive),OOai(NAct,INActive))

      Do IG=1,NGrid
C
      FPsiB(IG)=Zero
C
C     active-active
C
      Call dgemv('N',NCholesky*NAct,NAct,1d0,CHVCSAA,NCholesky*NAct,
     $           OrbGrid(IG,INActive+1:INActive+NAct),1,0d0,Oa,1)
C
C     Q(pq) = \sum_rs \Gamma_pqrs phi_r \phi_s
C
      Q=Zero
      Do IS=INactive+1,NOccup
      Do IR=INactive+1,NOccup
      Do IQ=INactive+1,NOccup
      Do IP=INactive+1,NOccup
      IPP=IP-INActive
      IQQ=IQ-INActive
      Q(IPP,IQQ)=Q(IPP,IQQ)
     &        +RDM2val(IP,IQ,IR,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Call dgemm('T','N',NAct,NAct,NCholesky,1d0,Oa,NCholesky,
     $           Oa,NCholesky,0d0,Oaa,NAct)
      FPsiB(IG)=ddot(NAct*NAct,Oaa,1,Q,1)
C
C     inactive-inactive
C
      If(INActive/=0) Then
C
      Call dgemv('N',NCholesky*INActive,INActive,1d0,CHVCSII,
     $          NCholesky*INActive,OrbGrid(IG,1:INActive),1,0d0,Oi,1)
C
      Call dgemm('T','N',INActive,INActive,NCholesky,1d0,Oi,NCholesky,
     $           Oi,NCholesky,0d0,Oii,INActive)
      Do IP=1,INActive
      Do IQ=1,INActive
      FPsiB(IG)=FPsiB(IG)+Oii(IP,IQ)*OrbGrid(IG,IP)*OrbGrid(IG,IQ)
      EndDo
      EndDo
C
C     active-inactive + inactive-active
C
      Call dgemv('N',NCholesky*NAct,INActive,1d0,CHVCSAI,
     &          NCholesky*NAct,OrbGrid(IG,1:INActive),1,0d0,Oai,1)
C
      Call dgemv('N',NCholesky*INActive,NAct,1d0,CHVCSIA,
     &          NCholesky*INActive,OrbGrid(IG,INactive+1:NOccup),
     &          1,0d0,Oia,1)
c
      Call dgemm('T','N',NAct,INActive,NCholesky,4d0,Oa,NCholesky,
     &           Oi,NCholesky,0d0,OOai,NAct)
c
      Call dgemm('T','N',NAct,INActive,NCholesky,-2d0,Oai,NCholesky,
     &           Oia,NCholesky,1d0,OOai,NAct)
c
      Do IQ=1,INActive
      Do IP=INactive+1,NOccup
      IPP=IP-INActive
      FPsiB(IG)=FPsiB(IG)
     &         +Occ(IP)*OOai(IPP,IQ)*OrbGrid(IG,IP)*OrbGrid(IG,IQ)
      EndDo
      EndDo
C
      EndIf ! INActive.ne.0

      EndDo ! NGrid
      FPsiB = 2d0*FPsiB
      If (IIPRINT.GT.1) Print*, 'Total FPsiB =', norm2(FPsiB)
C
      Deallocate(CHVCSII,CHVCSAA,CHVCSAI,CHVCSIA)
      Deallocate(OOai)
C
      call gclock('fr loop',Tcpu,Twall)
C
      If (IFlFCorr==0) Deallocate(WorkD) ! Chol vecs no longer necessary
C
      If (IFlFCorr==1) Then
c
C     add a contribution from the AC0 correlation Gamma
C
      Allocate(QGamII(NBasis,NBasis,NGrid))
      Allocate(QGamIA(NBasis,NBasis,NGrid))
      Allocate(QGamIV(NBasis,NBasis,NGrid))
      Allocate(QGamAA(NBasis,NBasis,NGrid))
      Allocate(QGamAV(NBasis,NBasis,NGrid))
      Allocate(QGamVV(NBasis,NBasis,NGrid))
!     careful! seems that some savings are
!              possible here :
!     Allocate(QGamVV(NOccup,NOccup,NGrid))
      Allocate(OGamI(NCholesky,NBasis,NGrid))
      Allocate(OGamA(NCholesky,NBasis,NGrid))
      Allocate(OGamV(NCholesky,NBasis,NGrid))
C
      QGamII=Zero
      QGamIA=Zero
      QGamIV=Zero
      QGamAA=Zero
      QGamAV=Zero
      QGamVV=Zero
C
      Allocate(FCorr(NGrid))
      FCorr=Zero
C
      Allocate (RDM2corr(NBasis,NBasis,NOccup,NOccup))
      RDM2corr=Zero
      inquire(file='rdmcorr',EXIST=ex)
      if (.not.ex) stop "no rdmcorr file!"
      Open(10,file='rdmcorr',form='unformatted')
CC
  557 Read(10,End=556) IP,IQ,IR,IS,GCor
c     print*, 'ip > ir, iq > is ',ip,ir,iq,is,GCor
      RDM2corr(IP,IQ,IR,IS) = Half*GCor
c
c     multiply by 1/4 (GCor is used four times) and by 2 because GCor (with proper factors)
c     must satisfy E^AC0 = 1/2 sum_pqrs GCor_pqrs <pqrs>
c     assumed symmetry of GCor pqrs = rqps = psrq = rspq
c
      GoTo 557
  556 Close(10)
C
      If (IIPRINT.GT.1) Print*, 'RDM2corr =',norm2(RDM2corr)
C
      Do IS=1,NBasis
      Do IR=1,NBasis
      IRS=IGem(IR)+(IGem(IS)-1)*3
      select case (IRS)
      case (1)
c        Print*, 'case 1: Inact-Inact'
         Do IG=1,NGrid
         Do IQ=INactive+1,NBasis
         Do IP=INactive+1,NBasis
         QGamII(IP,IQ,IG)=QGamII(IP,IQ,IG)
     &           +RDM2corr(IP,IQ,IR,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
         EndDo
         EndDo
         Enddo
      case (2)
c        Print*, 'case 2: Act - Inact'
         Do IG=1,NGrid
         Do IQ=Inactive+1,NBasis
            Do IP=1,IR-1
            QGamIA(IP,IQ,IG)=QGamIA(IP,IQ,IG)
     &              +RDM2corr(IR,IQ,IP,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
            Do IP=IR+1,NBasis
            QGamIA(IP,IQ,IG)=QGamIA(IP,IQ,IG)
     &              +RDM2corr(IP,IQ,IR,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
         EndDo
         Enddo
      case (3)
c        print*, 'case-3: Virt-Inact'
         Do IG=1,NGrid
         Do IQ=Inactive+1,NBasis
         Do IP=1,NOccup
         QGamIV(IP,IQ,IG)=QGamIV(IP,IQ,IG)
     &           +RDM2corr(IR,IQ,IP,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
         EndDo
         EndDo
         Enddo
c     case (4) IA = (AI)^T
      case (5)
c        Print*, 'case 5: Act-Act'
         Do IG=1,NGrid
         Do IQ=1,IS-1
            Do IP=1,IR-1
            QGamAA(IP,IQ,IG)=QGamAA(IP,IQ,IG)
     &              +RDM2corr(IR,IS,IP,IQ)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
            Do IP=IR+1,NBasis
            QGamAA(IP,IQ,IG)=QGamAA(IP,IQ,IG)
     &              +RDM2corr(IP,IS,IR,IQ)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
         EndDo
         Do IQ=IS+1,NBasis
            Do IP=1,IR-1
            QGamAA(IP,IQ,IG)=QGamAA(IP,IQ,IG)
     &              +RDM2corr(IR,IQ,IP,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
            Do IP=IR+1,NBasis
            QGamAA(IP,IQ,IG)=QGamAA(IP,IQ,IG)
     &              +RDM2corr(IP,IQ,IR,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
         EndDo
         EndDo
      case (6) ! Virt-Act
         Do IG=1,NGrid
         Do IP=1,NOccup
            Do IQ=1,IS-1
               QGamAV(IP,IQ,IG)=QGamAV(IP,IQ,IG)
     &           +RDM2corr(IR,IS,IP,IQ)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
            Do IQ=IS+1,NBasis
               QGamAV(IP,IQ,IG)=QGamAV(IP,IQ,IG)
     &           +RDM2corr(IR,IQ,IP,IS)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
            EndDo
         EndDo
         Enddo
c     case(7) IV = (VI)^T
c     case(8) AV = (VA)^T
      case (9)
c       Print*, 'case 9: Virt-Virt'
         Do IG=1,NGrid
         Do IQ=1,NOccup
         Do IP=1,NOccup
         QGamVV(IP,IQ,IG)=QGamVV(IP,IQ,IG)
     &           +RDM2corr(IR,IS,IP,IQ)*OrbGrid(IG,IR)*OrbGrid(IG,IS)
         EndDo
         EndDo
         Enddo
      end select
      EndDo
      EndDo
C
      If (IIPRINT.GT.1) Then
       Print*, 'QGam-1 II', norm2(QGamII)
       Print*, 'QGam-2 AI', norm2(QGamIA)
       Print*, 'QGam-3 VI', norm2(QGamIV)
       Print*, 'QGam-5 AA', norm2(QGamAA)
       Print*, 'QGam-6 VA', norm2(QGamAV)
       Print*, 'QGam-9 VV', norm2(QGamVV)
      EndIf
c
      Deallocate(RDM2corr)
C
C     OGam matrices
C
      Do IT=1,NBasis
      ITT=(IT-1)*NBasis
      select case (IGem(IT))
      case(1)
         Do IG=1,NGrid
            Do IP=1,NBasis
               Do K=1,NCholesky
                  OGamI(K,IP,IG)=OGamI(K,IP,IG)
     &                          +WorkD(K,ITT+IP)*OrbGrid(IG,IT)
               EndDo
            EndDo
         EndDo
      case(2)
         Do IG=1,NGrid
            Do IP=1,NBasis
               Do K=1,NCholesky
                  OGamA(K,IP,IG)=OGamA(K,IP,IG)
     &                          +WorkD(K,ITT+IP)*OrbGrid(IG,IT)
               EndDo
            EndDo
         EndDo
      case(3)
         Do IG=1,NGrid
            Do IP=1,NBasis
               Do K=1,NCholesky
                  OGamV(K,IP,IG)=OGamV(K,IP,IG)
     &                          +WorkD(K,ITT+IP)*OrbGrid(IG,IT)
               EndDo
            EndDo
         EndDo
      end select
      EndDo
C
      If (IIPRINT.GT.1) Then
       Print*, 'OGamI = ', norm2(OGamI)
       Print*, 'OGamA = ', norm2(OGamA)
       Print*, 'OGamV = ', norm2(OGamV)
      EndIf
C
      call gclock('QGam OGam',Tcpu,Twall)
c
      Deallocate(WorkD)
      Allocate (Work(NBasis,NBasis))
C
C     COMPUTE fcorr(r)
C
      Do IG=1,NGrid
      FCorr(IG)=Zero
c
C     1 inact-inact
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     &  OGamI(:,:,IG),NCholesky,OGamI(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +ddot(NBasis*NBasis,Work,1,QGamII(:,:,IG),1)
c     print*, 'II: IG,Fcorr',IG,FCorr(IG)
c
C     2 act-inact
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     $  OGamA(:,:,IG),NCholesky,OGamI(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +2d0*ddot(NBasis*NBasis,Work,1,QGamIA(:,:,IG),1)
c     print*, 'AI: IG,Fcorr',IG,FCorr(IG)
C
C     3 virt-inact
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     $  OGamV(:,:,IG),NCholesky,OGamI(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +2d0*ddot(NBasis*NBasis,Work,1,QGamIV(:,:,IG),1)
c     print*, 'AI: IG,Fcorr',IG,FCorr(IG)
c
c    5 act-act
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     &  OGamA(:,:,IG),NCholesky,OGamA(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +ddot(NBasis*NBasis,Work,1,QGamAA(:,:,IG),1)
c     print*, 'AA: IG,Fcorr',IG,FCorr(IG)
c
C     6 virt-act
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     &  OGamV(:,:,IG),NCholesky,OGamA(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +2d0*ddot(NBasis*NBasis,Work,1,QGamAV(:,:,IG),1)
c     print*, 'AV: IG,Fcorr',IG,FCorr(IG)
c
c    9 virt-virt
      Call dgemm('T','N',NBasis,NBasis,NCholesky,1d0,
     &  OGamV(:,:,IG),NCholesky,OGamV(:,:,IG),NCholesky,0d0,Work,NBasis)
      FCorr(IG)=FCorr(IG)
     &         +ddot(NBasis*NBasis,Work,1,QGamVV(:,:,IG),1)
c     print*, 'VV: IG,Fcorr',IG,FCorr(IG)
c
      EndDo ! NGrid
      If (IIPRINT.GT.1) Print*, 'Total FCorr =',norm2(FCorr)
c
      Deallocate(OGamV,OGamA,OGamI)
c
      call gclock('FCorr loop',Tcpu,Twall)
C
      EndIf ! IFlFCorr end
c
      If (IFlFCorr==1) Allocate(OnTopCorr(NGrid),XMuLocCorr(NGrid))
C
      Allocate(OnTop(NGrid))
      Allocate(XMuLoc(NGrid))
      Allocate(RhoGrid(NGrid))
      Allocate(Sigma(NGrid))

      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IS=1,NOccup
         ValS=Two*OrbGrid(I,IS)
         Do IR=1,NOccup
            ValRS=OrbGrid(I,IR)*ValS
            Do IQ=1,NOccup
               ValQRS=OrbGrid(I,IQ)*ValRS
               Do IP=1,NOccup
                  OnTop(I)=OnTop(I)
     &                 +RDM2val(IP,IQ,IR,IS)*OrbGrid(I,IP)*ValQRS
               EndDo
            EndDo
         EndDo
      EndDo
C
      XMuLoc(I)=Zero
C
      If(OnTop(I).Gt.1.D-8) Then
c     XMuLoc(I)=Prefac*FPsiB(I)/OnTop(I)
c     used in Kasia's code:
      XMuLoc(I)=SQRT(3.1415)/Two*FPsiB(I)/OnTop(I)
      EndIf
C
      EndDo ! NGrid

      If (IFlFCorr==1) Then
      Do I=1,NGrid
C
      OnTopCorr(I)=Zero
      Do IP=1,NBasis
      Do IQ=1,NBasis
      OnTopCorr(I)=OnTopCorr(I)
     $            +QGamII(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      OnTopCorr(I)=OnTopCorr(I)
     $            +QGamAA(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      OnTopCorr(I)=OnTopCorr(I)
     $            +QGamVV(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      OnTopCorr(I)=OnTopCorr(I)
     $            +2d0*QGamIA(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      OnTopCorr(I)=OnTopCorr(I)
     $            +2d0*QGamIV(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      OnTopCorr(I)=OnTopCorr(I)
     $            +2d0*QGamAV(IP,IQ,I)*OrbGrid(I,IP)*OrbGrid(I,IQ)
      EndDo
      EndDo
C
      XMuLocCorr(I)=Zero
C
      If(OnTop(I).Gt.1.D-8) Then
c     XMuLocCorr(I)=Prefac*(FPsiB(I)+FCorr(I))
c    $ /(OnTop(I)+OnTopCorr(I))
c     used in Kasia's code:
      XMuLocCorr(I)=SQRT(3.1415)/Two*(FPsiB(I)+FCorr(I))
     $ /(OnTop(I)+OnTopCorr(I))
      EndIf
C
      EndDo ! NGrid
      EndIf ! IFlFCorr end
c
      If (IIPRINT.GT.1) Then
       Print*, 'OnTop      =', norm2(OnTop)
       Print*, 'XMuLoc     =', norm2(XMuLoc)
       If(IFlFCorr==1) Print*, 'XMuLocCorr =', norm2(XMuLocCorr)
       If(IFlFCorr==1) Print*, 'OnTopCorr  =', norm2(OnTopCorr)
      EndIf
c
      call gclock('OnTop loop 1',Tcpu,Twall)
c
c     if(allocated(FPsiB)) print*, 'allocated FsiB!'
c     if(allocated(RDM2val)) print*, 'allocated rdm2val!'
      Deallocate(FPsiB)
      Deallocate(RDM2val)
C
      If (IFlFCorr==1) Then
         Deallocate(Work)
         Deallocate(QGamAV,QGamIV,QGamIA)
         Deallocate(QGamVV,QGamAA,QGamII)
         Deallocate(FCorr)
      EndIf
C
      Allocate(Zk(NGrid))
C
      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
C
      SPi=SQRT(3.141592653589793)
      Const=Three/Two/SPi/(One-SQRT(Two))
C
      AvMU=Zero
      CorrMD=Zero
      XEl=Zero
      XPairAv=Zero
      XPairCAv=Zero
      XPairAC0Av=Zero
C
C     Test: averaged mu
c     Alpha = 10.0
c     Print*, 'Use Averaged Alpha = ', Alpha
c     Print*, 'USe C=0.0 ..'
c     XMuLoc = Alpha ! replace local Mu with average Mu
c
      Do I=1,NGrid
C
      XMu=XMuLoc(I)
c     XMu=Alpha ! replace local Mu with average Mu
      AvMU=AvMU+XMu*RhoGrid(I)*WGrid(I)
      XEl=XEl+RhoGrid(I)*WGrid(I)
C
      Bet=Zero
      OnTopC=Zero
      If(XMuLoc(I).Gt.1.D-8) OnTopC=OnTop(I)/(One+Two/SPi/XMu)
      If(OnTopC.Ne.Zero) Then
      Bet=Const*Zk(I)/OnTopC
      C=5.0d0
      CorrMD=CorrMD+Zk(I)/(C*One+Bet*XMu**3)*WGrid(I)
c     print*, 'i,corrmd',I,CorrMD
C
CC     visualisation
C      CorrRZ=Zk(I)/(C*One+Bet*XMu**3)
C      If(Abs(RR(1,I)).Lt.1.D-8.And.Abs(RR(2,I)).Lt.1.D-8)
C     $ Write(6,'(5E15.6)')RR(3,I),RhoGrid(I),XMuLoc(I),
C     $ OnTop(I),CorrRZ
CC
      EndIf
C
      XPairAv=XPairAv+OnTop(I)*WGrid(I)
      XPairCAv=XPairCAv+OnTopC*WGrid(I)
c     XPairAC0Av=XPairAC0Av+(OnTop(I)+OnTopCorr(I))*WGrid(I)
C
      EndDo
      call gclock('CorrMD loop ',Tcpu,Twall)
C
      Write(6,'(/,
     $" Averaged on-top pair densities: CAS, extrapolated OT",
     $ 2F15.4)') XPairAv,XPairCAv
C
      AvMU=AvMU/XEl
C
      If (IFlFCorr==1) Then
C
C     CorrMD with a contribution from Gamma^corr
C
      AvMUCorr=Zero
      CorrMDCorr=Zero
      XEl=Zero
      Do I=1,NGrid
C
      XPairAC0Av=XPairAC0Av+(OnTop(I)+OnTopCorr(I))*WGrid(I)
C
      XMu=XMuLocCorr(I)

      AvMUCorr=AvMUCorr+XMu*RhoGrid(I)*WGrid(I)
      XEl=XEl+RhoGrid(I)*WGrid(I)
C
      Bet=Zero
      OnTopC=Zero
c     If(XMu.Gt.1.D-8) OnTopC=(OnTop(I)+OnTopCorr(I))
      If(XMu.Gt.1.D-8) OnTopC=OnTop(I)
     $ /(One+Two/SPi/XMu)
      If(OnTopC.Ne.Zero) Then
      Bet=Const*Zk(I)/OnTopC
      C=5.0d0
      CorrMDCorr=CorrMDCorr+Zk(I)/(C*One+Bet*XMu**3)*WGrid(I)
C
CC     visualisation
C      CorrRZ=Zk(I)/(C*One+Bet*XMu**3)
C      If(Abs(RR(1,I)).Lt.1.D-8.And.Abs(RR(2,I)).Lt.1.D-8)
C     $ Write(6,'(5E15.6)')RR(3,I),RhoGrid(I),XMuLocCorr(I),
C     $ OnTop(I)+OnTopCorr(I),CorrRZ
C
      EndIf
C
      EndDo
      call gclock('CorrMDCorr loop ',Tcpu,Twall)
C
      Write(6,'(/,
     $" Averaged on-top pair densities: CAS, CAS-AC0, extrapolated OT",
     $ 3F15.4)') XPairAv,XPairAC0Av,XPairCAv
C
      AvMUCorr=AvMUCorr/XEl
      EndIf ! IFlFCorr end
C
      Write(6,'(/," CAS:     Corr_md and average Mu   ",
     $ F12.8,F12.2)') CorrMD,AvMU
      If(IFlFCorr==1) Then
         Write(6,'(  " CAS-AC0: Corr_md and average Mu   ",
     $         F12.8,F12.2/)') CorrMDCorr,AvMUCorr
C
         Deallocate(OnTopCorr,XMuLocCorr)
      EndIf
C
      End Subroutine LOC_MU_CHOL

*Deck PBE_ONTOP_MD
      Subroutine PBE_ONTOP_MD(PBEMD,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A SR-PBE ONTOP CORRELATION DEFINED IN Eq.(25)-(29), Toulouse JCP 150, 084103 (2019)
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
      Dimension Zk(NGrid),RhoGrid(NGrid),Sigma(NGrid),OnTop(NGrid)
C
      Write(6,'(/,1X,"COMPUTING SR-PBE_CORR_MD FOR THE RANGE PARAMETER "
     $ ,F8.3)')Alpha
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
      Open(newunit=iunit,File="rdm2.dat",Status='Old')
C
      Do
      Read(iunit,*,iostat=ios)I,J,K,L,X
      If (ios/=0) Exit
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
      EndDo
      Close(iunit)
C
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop(I)=OnTop(I)
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      EndDo
C
      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
C
      SPi=SQRT(3.141592653589793)
      Const=Three/Two/SPi/(One-SQRT(Two))
C
      PBEMD=Zero
      Do I=1,NGrid
C
      Bet=Zero
      OnTopC=OnTop(I)/(One+Two/SPi/Alpha)
      If(OnTopC.Ne.Zero) Bet=Const*Zk(I)/OnTopC
C
      PBEMD=PBEMD+Zk(I)/(One+Bet*Alpha**3)*WGrid(I)
C
      EndDo
C
      Return
      End

*Deck PBE_ONTOP_C_MD
      Subroutine PBE_ONTOP_C_MD(PBEmodMD,Cfac,URe,Occ,
     $ OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
C
C     RETURNS A SR-PBE modified ONTOP CORRELATION
C     DEFINED IN Eq.(25)-(29), Toulouse JCP 150, 084103 (2019)
C     Difference:
C     Eq. (26): e_{c,md} = e_c^PBE / (C + beta*mu^3)
C     Recommended value: Cfac = 5.0d0
C
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Double Precision, intent(in)  :: Cfac
      Double Precision, intent(out) :: PBEmodMD
C
      Real*8, Allocatable :: RDM2Act(:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),
     $ WGrid(NGrid),OrbGrid(NGrid,NBasis),OrbXGrid(NGrid,NBasis),
     $ OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)
C
      Dimension Zk(NGrid),RhoGrid(NGrid),Sigma(NGrid),OnTop(NGrid)
C
      Write(6,'(1X,"COMPUTING SR-PBE_CORR_MD FOR THE RANGE PARAMETER "
     $ ,F8.3)')Alpha
      Write(6,'(1x, "C FACTOR ", F6.3)')Cfac
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
      Open(10,File="rdm2.dat",Status='Old')
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
      Do I=1,NGrid
C
      Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
      Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
      Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
      Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
      Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
C
      OnTop(I)=Zero
      Do IP=1,NOccup
      Do IQ=1,NOccup
      Do IR=1,NOccup
      Do IS=1,NOccup
      OnTop(I)=OnTop(I)
     $ +Two*FRDM2(IP,IQ,IR,IS,RDM2Act,Occ,Ind2,NAct,NBasis)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      EndDo
C
      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
C
      SPi=SQRT(3.141592653589793)
      Const=Three/Two/SPi/(One-SQRT(Two))
C
      PBEmodMD=Zero
      Do I=1,NGrid
C
      Bet=Zero
      OnTopC=OnTop(I)/(One+Two/SPi/Alpha)
      If(OnTopC.Ne.Zero) Then
         Bet=Const*Zk(I)/OnTopC
c     If(OnTopC.Ne.Zero) Bet=Const*Zk(I)/OnTopC
C
         PBEmodMD=PBEmodMD+Zk(I)/(Cfac*One+Bet*Alpha**3)*WGrid(I)
c        print*, 'i,pbemod',I,PBEmodMD
      EndIf
C
      EndDo
C
      Return
      End

      Subroutine PsiHPSi(Ene,Occ,XOne,ENuc,INActive,NAct,NInte1,NBasis)
C
C     Returns E = <Psi|T+V+W|Psi> , where
C     W is the full 1/r operator
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Integer,Intent(in) :: INActive,NAct,NInte1,NBasis
      Double Precision,Intent(in)  :: ENuc
      Double Precision,Intent(in)  :: Occ(NBasis),XOne(NInte1)
      Double Precision,Intent(out) :: Ene
C
      Integer :: I,J,K,L,II
      Integer :: IFunSR_save
      Integer :: NRDM2Act
      Integer :: Ind1(NBasis)
      Double Precision :: X,EOne,ETwo
      Double Precision, Allocatable :: RDM2Act(:)

c     Write(6,*) 'Entering <Psi|H|Psi>...'
c
      Do I=1,NAct
      Ind1(I)=INActive+I
C     Ind2(INActive+I)=I
      EndDo
C
      NRDM2Act=NAct**2*(NAct**2+1)/2
      Allocate(RDM2Act(NRDM2Act))
      Open(10,File="rdm2.dat",Status='Old')
C
   10 Read(10,*,End=40)I,J,K,L,X
C
C     X IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
C
      RDM2Act(NAddrRDM(J,L,I,K,NAct))=0.5d0*X
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
      EOne=0d0
      Do I=1,NBasis
      II=(I*(I+1))/2
      EOne=EOne+2d0*Occ(I)*XOne(II)
      EndDo
C
      IFunSR_save = IFunSR ! select FOFO integrals
      IFunSR = 0
      Call TwoEneChck(ETwo,RDM2Act,Occ,INActive,NAct,NBasis)
      IFunSR = IFunSR_save
C
      Ene = EOne + ETwo + ENuc
      Print*, '<Psi|T+Vne|Psi>  =',EOne
      Print*, '<Psi|Vee  |Psi>  =',ETwo
      Write(6,'(/1x,a,F14.8)') '<Psi|H|Psi>+ENuc = ',Ene
C
      Deallocate(RDM2Act)

      End

