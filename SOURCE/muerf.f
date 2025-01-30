*Deck DBBSCH
      Subroutine DBBSCH(ETot,ENuc,URe,Occ,XOne,UNOAO,
     $                  BasisSet,NBasis,NInte1,IndN,IndX,NDimX)
C
C     Compute the CBS[H] correction
C     with SR integrals composed from Cholesky vectors
C     CBS[DFT] correction [E. Giner et al. 2018] is a byproduct
C
      use ab0fofo
      use timing
C
      Implicit Real*8 (A-H,O-Z)

      include 'commons.inc'

      character(*)       :: BasisSet
      integer,intent(in) :: NBasis,NInte1,NDimX
      integer,intent(in) :: IndX(NDimX),IndN(2,NDimX)

      double precision,intent(in)  :: ENuc
      double precision,intent(out) :: ETot
      double precision,intent(in) :: Occ(NBasis),XOne(NInte1)
      double precision,intent(in) :: URe(NBasis,NBasis),
     $                               UNOAO(NBasis,NBasis)
C
C     local
C
      integer :: iunit
      integer :: NCholesky,NCholErf
      double precision :: ECASSCF,ECorr
      double precision :: AvMu,ECorrMD
      double precision :: XMuMat(NBasis,NBasis)
      double precision :: ABMIN(NDimX,NDimX),ABPLUS(NDimX,NDimX)
      double precision, allocatable :: CholVecs(:,:),Work(:,:)
C
      double precision :: Tcpu,Twall

      If (ITwoEl.eq.1) then
         write(6,'(1x,a)') "Set PostCAS=.true. and DFunc for DBBSCH!"
         stop "Fix keywords in the Calculation block!"

      ElseIf (IDBBSC.eq.1.and.ITwoEl.eq.3.and.ICholesky.ne.0) Then
C
      Call LOC_MU_CBS_CHOL(XMuMat,ECorrMD,AvMU,URe,UNOAO,Occ,
     $                     BasisSet,NBasis)

      Call AC0CAS_FOFO(ECorr,ECASSCF,Occ,URe,XOne,ABPLUS,ABMIN,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,NInte1,
     $ NoSt,'FFOO','FOFO',ICholesky,IDBBSC,IFlFCorr)
C
      Write
     $ (6,'(/1X,''CASSCF+ENuc, AC0-Corr    , Total'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
      Write
     $ (6,'(1X,''CASSCF+ENuc, AC0-CBS[DFT], Total'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr+ECorrMD,ECASSCF+ENuc+ECorr+ECorrMD
       ETot=ECASSCF+ENuc+ECorr
C
      ElseIf (IDBBSC.eq.2.and.ITwoEl.eq.3.and.ICholesky.ne.0) Then
C
C     Compute local mu(r): XMuMat(p,q) = <p|mu(r)|q>
      Call LOC_MU_CBS_CHOL(XMuMat,CorrMD,AvMU,URe,UNOAO,Occ,
     $                     BasisSet,NBasis)
c     Print*, 'LOC_MU_CBS_CHOL: XMuMat = ',norm2(XMuMat)
C
C     ... test : recover AC0!
c      XMuMat=0d0
c      print*, 'Use FR only !',norm2(XMuMat)

C     transform full-range (FR) cholesky vecs
C     R(k,pq).(q|mu(r)|t) = R(k,pt)
      open(newunit=iunit,file='cholvecs',form='unformatted')
      read(iunit) NCholesky
      Allocate(CholVecs(NCholesky,NBasis**2),Work(NCholesky,NBasis**2))
      read(iunit) CholVecs
      close(iunit)

C     SET TIMING FOR 3-ind transformations of Cholesky vecs
      Call gclock('START',Tcpu,Twall)
      Call dgemm('N','N',NCholesky*NBasis,NBasis,NBasis,1d0,
     $           CholVecs,NCholesky*NBasis,XMuMat,NBasis,
     $           0d0,Work,NCholesky*NBasis)
      Call gclock('3-idx tran FR(k,pq)',Tcpu,Twall)
C     save to disk
      open(newunit=iunit,file='chol1vFR',form='unformatted')
      write(iunit) NCholesky
      write(iunit) Work
      close(iunit)
      deallocate(CholVecs,Work)
C
C     transform long-range (LR) cholesky vecs
      open(newunit=iunit,file='cholvErf',form='unformatted')
      read(iunit) NCholErf
      Allocate(CholVecs(NCholErf,NBasis**2),Work(NCholErf,NBasis**2))
      read(iunit) CholVecs
      close(iunit)
      Call dgemm('N','N',NCholErf*NBasis,NBasis,NBasis,1d0,
     $           CholVecs,NCholErf*NBasis,XMuMat,NBasis,
     $           0d0,Work,NCholErf*NBasis)
      Call gclock('3-idx tran LR(k,pq)',Tcpu,Twall)
C     save to disk
      open(newunit=iunit,file='chol1vLR',form='unformatted')
      write(iunit) NCholErf
      write(iunit) Work
      close(iunit)
      deallocate(CholVecs,Work)

      Call AC0CAS_FOFO(ECorr,ECASSCF,Occ,URe,XOne,ABPLUS,ABMIN,
     $ IndN,IndX,IGem,NAcCAS,NInAcCAS,NDimX,NBasis,NDimX,NInte1,
     $ NoSt,'FFOO','FOFO',ICholesky,IDBBSC,IFlFCorr)

      Call FirstOrderSREne(ETwoSR,NBasis)
      Write(6,'(/,X,"1st-order SR energy for",I4," active orbitals")') 
     $ NAcCAS+NInAcCAS-NStronglyOccOrb
      Write(6,'(X," <Ref|H^SR|Ref> = ",F10.6)') ETwoSR

c     print*, 'ABPLUS-my =',norm2(ABPLUS)
c     print*, 'ABMIN -my =',norm2(ABMIN)

      Write
     $ (6,'(/1X,''CASSCF+ENuc, AC0-CBS[H], Total'',6X,3F15.8)')
     $ ECASSCF+ENuc,ECorr,ECASSCF+ENuc+ECorr
       ETot=ECASSCF+ENuc+ECorr

c     ... delete transformed cholesky vecs
      Open(newunit=iunit,file='chol1vFR',status='OLD')
      Close(iunit,status='DELETE')
      Open(newunit=iunit,file='chol1vLR',status='OLD')
      Close(iunit,status='DELETE')
c
c      Print*, 'ECASSCF = ', ECASSCF
c      Print*, 'ENuc = ', ENuc
c      Print*, 'ECorr = ', ECorr

      EndIf ! ITwoEl, IDBBSC

      End
C *End Subroutine DBBSCH

*Deck LOC_MU_CBS_CHOL
      Subroutine LOC_MU_CBS_CHOL(XMuMat,CorrMD,AvMU,URe,
     $                           UNOAO,Occ,BasisSet,NBasis)
C
      use read_external
      use grid_internal
      use timing
C
C     RETURNS A "BASIS-SET ERROR CORRECTION" WITH SR-PBE ONTOP CORRELATION from Giner et al. JCP 152, 174104 (2020)
C     RETURNS XMuMAT USED TO COMPUTE CBS CORRECTION BY MODIFICATION OF H' IN AC0
C
C     CAREFUL!
C     UNOAO are U(NO,AO) not U(NO,SAO) orbitals!
C
      Implicit Real*8 (A-H,O-Z)
C
C     now defined in GammCor_Integrals
      Parameter(Half=0.5D0)
C     Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
C    $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis)
      Real*8 :: XMuMat(NBasis,NBasis)
      Character(*) :: BasisSet
C
      Real*8 :: UNOSAO(NBasis,NBasis)
      Real*8, Allocatable :: FPsiB(:), OnTop(:), XMuLoc(:)
      Real*8,allocatable,target :: PhiLDA(:,:)
      Real*8,allocatable,target :: PhiGGA(:,:,:)
      Real*8,pointer :: OrbGrid(:,:)
      Real*8,pointer :: OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:)
      Real*8, Allocatable :: Zk(:),RhoGrid(:),Sigma(:),WGrid(:)
      Double Precision, Allocatable :: CHVCSAF(:,:),CHVCSIF(:,:)
      Double Precision, Allocatable :: Q(:,:),OrbGridP(:)
      Double Precision, Allocatable :: Oi(:,:),Oa(:,:),Oaa(:,:)
      Double Precision, Allocatable :: tOi(:),tOa(:)
      Double Precision, Allocatable :: RDM2Act(:)
      Double Precision, Allocatable :: RDM2val(:,:,:,:)
      Double Precision, Allocatable :: Work(:,:),WorkD(:,:)
      Dimension IAct(NBasis),Ind2(NBasis)
      Dimension IndInt(NBasis),NSymNO(NBasis),NumOSym(15),MultpC(15,15)
      Dimension OrbIGGrid(NBasis)
C
      double precision :: Twall,TCpu
C
      Logical :: doGGA
      Double Precision, Allocatable :: RR(:,:)
      Character(*),Parameter :: griddalfile='dftgrid.dat'
C
      call gclock('START',Tcpu,Twall)
C
      If(IFlCore.Eq.0.And.IDBBSC.Eq.1) Then
         ICore = NCoreOrb ! from input.inp
         Do I=1,ICore
         Occ(I)=Zero
         EndDo
      Else
         ICore = 0
      EndIf
C
C     debug printlevel
      IIPRINT=5
C
C     set constants
C
C      Pi = dacos(-1.d0)
C      SPi = SQRT(Pi)
C      Prefac = SQRT(3.1415)/Two ! SPi/Two
C
      Write(lout,'(/1x,5a6)') ('******',i=1,5)
      Write(lout,'(1x,a)') "CBS Correction with local Mu"
      Write(lout,'(1x,5a6)') ('******',i=1,5)
C
      If (InternalGrid==0) then
C
      If (IMOLPRO == 1) Then
         Write(LOUT,'(/1x,a)') 'MOLPRO GRID '
         Call molprogrid0(NGrid,NBasis)
      ElseIf(IDALTON == 1) Then
         Write(LOUT,'(/1x,a)') 'DALTON GRID '
         Call daltongrid0(NGrid,griddalfile,doGGA,NBasis)
         If (.not.doGGA) stop "LOC_MU_CHOL: Run Dalton with SR-PBE!"
      EndIf
      Write(6,'(/,1X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate (WGrid(NGrid))
      Allocate(PhiGGA(NGrid,NBasis,4))
      OrbGrid  => PhiGGA(:,:,1)
      OrbXGrid => PhiGGA(:,:,2)
      OrbYGrid => PhiGGA(:,:,3)
      OrbZGrid => PhiGGA(:,:,4)
C
C     load orbgrid, gradients, and wgrid
C
      If (IMOLPRO == 1) Then
         Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                   WGrid,UNOAO,NGrid,NBasis)
      ElseIf (IDALTON == 1) Then
         Allocate(Work(NBasis,NBasis))
         Work = transpose(UNOAO)
         Call daltongrid_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                       WGrid,NGrid,griddalfile,NBasis)
         Call daltongrid_tran_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     &                       Work,0,NGrid,NBasis,.false.)
CC     visualisation
C        Allocate(RR(3,NGrid))
C        Call dalton_grid_coord(NGrid,RR,griddalfile)
        Deallocate(Work)
      EndIf !Molpro/Dalton grid
c
      ElseIf (InternalGrid==1) then
c        Write(6,'(1x,a,i3)') 'IFunSR       =',IFunSR
c        Write(6,'(1x,a,i3)') 'IGridType    =', IGridType
c        Write(6,'(1x,a,i3)') 'ORB_ORDERING =', IOrbOrder
C
         Allocate(Work(NBasis,NBasis))
         Work = transpose(UNOAO)
         Call internal_gga_no_orbgrid(IGridType,BasisSet,
     $                          IOrbOrder,Work,
     $                          WGrid,PhiGGA,NGrid,NBasis,NBasis,IUnits)
         OrbGrid  => PhiGGA(:,:,1)
         OrbXGrid => PhiGGA(:,:,2)
         OrbYGrid => PhiGGA(:,:,3)
         OrbZGrid => PhiGGA(:,:,4)

         Deallocate(Work)
C
      EndIf ! InternalGrid
C
c      print*, ''
c      print*, 'NGrid    = ', NGrid
c      print*, 'WGrid    = ', norm2(Wgrid)
c      print*, 'OrbGrid  = ', norm2(OrbGrid)
c      print*, 'OrbXGrid = ', norm2(OrbXGrid)
c      print*, 'OrbYGrid = ', norm2(OrbYGrid)
c      print*, 'OrbZGrid = ', norm2(OrbZGrid)
c      print*, ''
C
C     ... symmetry
      If (InternalGrid==1) Then
C     ... internal grid does not use symmetry
c        write(lout,'(/1x,a)') "Internal grid does not use symmetry!"
         NSym=1
         MxSym=1
         NumOSym(1)=NBasis
      ElseIf(InternalGrid==0) Then
         If (IMOLPRO == 1) Then 
            Call create_ind_molpro('2RDM',NumOSym,IndInt,NSym,NBasis)
            MxSym=NSym
         Else
            Stop "Fix NSym in LOC_MU_CBS_CHOL!"
         EndIf
      EndIf
C
C     ...test prints
c     print*, 'NSym    =', NSym
c     print*, 'NumOSym =', NumOSym(1:MxSym)
C
C      print*, 'UNOSAO orbitals:',norm2(UNOAO)
C      do j=1,NBasis
C         write(LOUT,'(*(f13.8))') (UNOAO(i,j),i=1,NBasis)
C      enddo
C
      NSymNO(1:NBasis)=0
      IStart=0
      Do I=1,MxSym
      Do J=IStart+1,IStart+NumOSym(I)
C
      Do IOrb=1,NBasis
      If(Abs(UNOAO(IOrb,J)).Gt.1.D-1) NSymNO(IOrb)=I
      EndDo
C
      EndDo
      IStart=IStart+NumOSym(I)
      EndDo

C      Do IOrb=1,NBasis
C      Print*, 'I,NSymNO = ',IOrb,NSymNO(IOrb)
C      EndDo
C      Print*, ''
C
C     check...
      Do I=1,MxSym
      II=0
      Do IOrb=1,NBasis
      If(NSymNO(IOrb).Eq.I) II=II+1
      EndDo
      If(II.Ne.NumOSym(I))
     $ Write(*,*) 'Symmetry of NO cannot be established!'
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
C     ... end symmetry

C
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=INActive+NAct
      Ind2(1:NBasis)=0
      Do I=1,NAct
      Ind2(INActive+I)=I
      EndDo
C
C     ... more test prints
C      print*, 'NAct    =', NAct
C      print*, 'INActive=', INActive
C      print*, 'NOccup  =', NOccup
C      print*, 'ICore   =', ICore
C
      if (IFlCore.Eq.0) Then
         INActiveC=INactive-ICore ! skip core within inactive
      else
         INActiveC = INActive
      endif
c     PRint*, 'INActiveC =', INActiveC
C
c     Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
c     Call get_RDM2Occ(RDM2val,Occ,NAct,NOccup,NBasis)
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
C     change later from 1212 to 11,22!
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
C
C     truncate CHOLVECS
C
      Allocate (CHVCSAF(NCholesky,NAct*NBasis))
      Do K=1,NCholesky
         IIJJ=0
         Do J=1,NBasis
            Do I=INActive+1,NOccup
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSAF(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
      If (INActiveC.Gt.0) Then
      Allocate (CHVCSIF(NCholesky,INActiveC*NBasis))
      ! use WorkD instead?
      Do K=1,NCholesky
         IIJJ=0
         Do J=1,NBasis
            Do I=ICore+1,INActive
               IIJJ=IIJJ+1
               IJ=I+(J-1)*NBasis
               CHVCSIF(K,IIJJ) = WorkD(K,IJ)
            EndDo
         EndDo
      EndDo
      EndIf !INActive.gt.0

      Deallocate(WorkD)
C
C     COMPUTE f(r)
C
      Allocate (FPsiB(NGrid))
      Allocate (OrbGridP(NAct))
      Allocate (Q(NAct,NAct))
      Allocate (Oa(NCholesky,NAct))
      Allocate (Oaa(NAct,NAct))
      Allocate (tOi(NCholesky),tOa(NCholesky))
      if (INActiveC.gt.0) Allocate (Oi(NCholesky,INActiveC))
c
C     Allocate (Oii(INActive,INActive))
c     Allocate (OOai(NAct,INActive))

      Do IG=1,NGrid
C
      FPsiB(IG)=Zero

      OrbGridP=Zero
      Do IP=1,NAct
      OrbGridP(IP)=Occ(INActive+IP)*OrbGrid(IG,INActive+IP)
      EndDo
C
C     active-active
C
      OrbIGGrid(1:NBasis) = OrbGrid(IG,1:NBasis)
      Call dgemv('N',NCholesky*NAct,NBasis,1d0,CHVCSAF,NCholesky*NAct,
     $           OrbIGGrid,1,0d0,Oa,1)
c    $           OrbGrid(IG,1:NBasis),1,0d0,Oa,1)
c     print*, 'Oa =',norm2(Oa)
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
c     print*, '1st FPsiB =', IG, FPsiB(IG)
C
c     If(IFlCore.Ne.0) Then
      If (INActiveC.Gt.0) Then
C
C     inactive-inactive
C
      OrbIGGrid(1:NBasis) = OrbGrid(IG,1:NBasis)
      Call dgemv('N',NCholesky*INActiveC,NBasis,1d0,CHVCSIF,
     $          NCholesky*INActiveC,OrbIGGrid,1,0d0,Oi,1)
!C
!     ver 1
!      Call dgemm('T','N',INActive,INActive,NCholesky,1d0,Oi,NCholesky,
!     $           Oi,NCholesky,0d0,Oii,INActive)
!c     print*, 'Oi =',norm2(Oi)
!      Do IP=1,INActive
!      Do IQ=1,INActive
!      FPsiB(IG)=FPsiB(IG)+Oii(IP,IQ)*OrbGrid(IG,IP)*OrbGrid(IG,IQ)
!      EndDo
!      EndDo
!c     print*, '2nd FPsiB =', IG, FPsiB(IG)
C
!     ver 2
      OrbIGGrid(1:INActiveC) = OrbGrid(IG,ICore+1:INActive)
      Call dgemv('N',NCholesky,INActiveC,1d0,Oi,
     $          NCholesky,OrbIGGrid,1,0d0,tOi,1)
C    $          NCholesky,OrbGrid(IG,ICore+1:INActive),1,0d0,tOi,1)
      FPsiB(IG)=FPsiB(IG)+ddot(NCholesky,tOi,1,tOi,1)
!
C     active-inactive
C
CC     ver 1 : NChol*NAct*INact scaling
C      Call dgemm('T','N',NAct,INActive,NCholesky,2d0,Oa,NCholesky,
C     &           Oi,NCholesky,0d0,OOai,NAct)
CC
C      Do IQ=1,INActive
C      Do IP=INactive+1,NOccup
C      IPP=IP-INActive
C      FPsiB(IG)=FPsiB(IG)
C     &         +Occ(IP)*OOai(IPP,IQ)*OrbGrid(IG,IP)*OrbGrid(IG,IQ)
C      EndDo
C      EndDo
C
C     ver2 : NChol*NOccup scaling
      Call dgemv('N',NCholesky,NAct,1d0,Oa,
     $          NCholesky,OrbGridP,1,0d0,tOa,1)
      FPsiB(IG)=FPsiB(IG)+2d0*ddot(NCholesky,tOa,1,tOi,1)
C
C     If(IFlCore.Ne.0)
c     EndIf
      EndIf ! INActiveC.Gt.0
C
      EndDo ! IG=1,NGrid
C
      Deallocate(CHVCSAF)
      If (INActiveC.Gt.0) Deallocate(CHVCSIF)
C
      FPsiB = 2d0*FPsiB
      If (IIPRINT.GT.1) Print*, 'Total FPsiB =', norm2(FPsiB)
C
      Allocate(OnTop(NGrid))
      Allocate(XMuLoc(NGrid))
      Allocate(RhoGrid(NGrid))
      Allocate(Sigma(NGrid))
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
      XMuLoc(I)=SQRT(3.1415)/Two*FPsiB(I)/OnTop(I)
      EndIf
C
      EndDo ! NGrid
C
      Print*, 'OnTop', norm2(OnTop)
      print*, 'XMuLoc',norm2(XMuLoc)
C
C     COMPUTE XMuMAT MATRIX in NO
C
      C=1.0d0
      Write(*,*)
      Write(*,*)"*** C= *** ",C
C
      Do IP=1,NBasis
      Do IQ=1,IP
      XMuMAT(IP,IQ)=0d0
C
      ISym=MultpC(NSymNO(IP),NSymNO(IQ))
      If(ISym.Eq.1) Then
C
      Do I=1,NGrid
      XMuMAT(IP,IQ)=XMuMAT(IP,IQ)
     $ +OrbGrid(I,IP)*OrbGrid(I,IQ)*WGrid(I)*Exp(-C*XMuLoc(I))
      EndDo
C
      EndIf
      XMuMAT(IQ,IP)=XMuMAT(IP,IQ)
C
      EndDo
      EndDo
C
C     COMPUTE CBS CORRECTION
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
      Do I=1,NGrid
C
      XMu=XMuLoc(I)
C
      AvMU=AvMU+XMu*RhoGrid(I)*WGrid(I)
      XEl=XEl+RhoGrid(I)*WGrid(I)
C
      Bet=Zero
      OnTopC=Zero
      If(XMuLoc(I).Gt.1.D-8) OnTopC=OnTop(I)/(One+Two/SPi/XMu)
      If(OnTopC.Ne.Zero) Then
      Bet=Const*Zk(I)/OnTopC
      C=1d0
      CorrMD=CorrMD+Zk(I)/(C*One+Bet*XMu**3)*WGrid(I)
      EndIf
C
      EndDo
C
      AvMU=AvMU/XEl
C
      If(IFlCore.Eq.0) Then
      Do I=1,ICore
c     Do I=1,NInAcCAS
      Occ(I)=One
      EndDo
      Write(6,'(/,1X,"*** CBS CORRECTION EXCLUDING CORE ORBITALS ***")')
      Write(6,'(1X,"*** NO. OF CORE ORBITALS :",I3)') ICore
      EndIf
C
      Write(6,'(/,
     $" CBS[DFT] correction, average Mu",
     $ F15.8,F15.3/)') CorrMD,AvMU
C
      call gclock('Esrcmd CBS ',Tcpu,Twall)

      Deallocate(RDM2val)

      End Subroutine LOC_MU_CBS_CHOL

*Deck LOC_MU_CBS
      Subroutine LOC_MU_CBS(XMuMAT,URe,UNOAO,Occ,TwoEl,
     $ NBasis,NInte2)
C
C     INCORE VERSION
C
C     COMPUTES A "BASIS-SET ERROR CORRECTION" WITH SR-PBE ONTOP CORRELATION from Giner et al. JCP 152, 174104 (2020)
C     RETURNS XMuMAT USED TO COMPUTE CBS CORRECTION BY MODIFICATION OF H' IN AC0
C
      use read_external
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension NSymNO(NBasis),MultpC(15,15),NumOSym(15),IndInt(NBasis)
C
      Real*8, Allocatable :: RDM2Act(:),FPsiB(:,:,:,:),
     $ OrbGrid(:,:),OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:),
     $ Zk(:),RhoGrid(:),Sigma(:),WGrid(:),OnTop(:),XMuLoc(:)
C
      Dimension TwoEl(NInte2),URe(NBasis,NBasis),Occ(NBasis),
     $ Ind1(NBasis),Ind2(NBasis),UNOAO(NBasis,NBasis),
     $ XMuMAT(NBasis,NBasis)
C
      Dimension IAct(NBasis)
C
c      If(IFlCore.Eq.0) Then
c      Do I=1,NInAcCAS
c      Occ(I)=Zero
c      EndDo
c      EndIf
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,1X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate (WGrid(NGrid))
      Allocate (OrbGrid(NGrid,NBasis))
      Allocate (OrbXGrid(NGrid,NBasis))
      Allocate (OrbYGrid(NGrid,NBasis))
      Allocate (OrbZGrid(NGrid,NBasis))
      Allocate(OnTop(NGrid),XMuLoc(NGrid))
      Allocate(RhoGrid(NGrid))
      Allocate(Sigma(NGrid))
      Allocate(Zk(NGrid))
C
      Call create_ind_molpro('2RDM',NumOSym,IndInt,NSym,NBasis)
      MxSym=NSym

!     print*, 'NSym =',NSym
!     print*, 'NumOSym =', NumOSym(1:MxSym)

!      print*, 'INCORE: UNOSAO',norm2(UNOAO)
!      do j=1,NBasis
!         write(LOUT,'(*(f13.8))') (UNOAO(i,j),i=1,NBasis)
!      enddo
C
      NSymNO(1:NBasis)=0
      IStart=0
      Do I=1,MxSym
      Do J=IStart+1,IStart+NumOSym(I)
C
      Do IOrb=1,NBasis
      If(Abs(UNOAO(IOrb,J)).Gt.1.D-1) NSymNO(IOrb)=I
      EndDo
C
      EndDo
      IStart=IStart+NumOSym(I)
      EndDo
C
c     Do IOrb=1,NBasis
c     Print*, 'I,NSymNO = ',IOrb,NSymNO(IOrb)
c     EndDo
c     Print*, ''
C
C     checking
      Do I=1,MxSym
      II=0
      Do IOrb=1,NBasis
      If(NSymNO(IOrb).Eq.I) II=II+1
      EndDo
      If(II.Ne.NumOSym(I))
     $ Write(*,*) 'Symmetry of NO cannot be established!'
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
C     load orbgrid, gradients, and wgrid
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
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
      Do I=1,NBasis
      IAct(I)=0
      If(Occ(I).Ne.One.And.Occ(I).Ne.Zero) IAct(I)=1
      EndDo
C
      Allocate (FPsiB(NBasis,NBasis,NOccup,NOccup))
      FPsiB=Zero
C
      Do I=1,NBasis
      Do J=1,NBasis
      Do M=1,NOccup
      Do N=1,NOccup
C
      FPsiB(I,J,M,N)=Zero
      Do K=1,NOccup
      Do L=1,NOccup
      FPsiB(I,J,M,N)=FPsiB(I,J,M,N)+
     $ TwoEl(NAddr3(K,I,L,J))
     $ *Two*FRDM2(K,L,M,N,RDM2Act,Occ,Ind2,NAct,NBasis)
      EndDo
      EndDo
C
      EndDo
      EndDo
      EndDo
      EndDo
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
      XMuLoc(I)=Zero
C
      If(OnTop(I).Gt.1.D-8) Then
C
      Do IP=1,NBasis
      Do IQ=1,NBasis
      Do IR=1,NOccup
      Do IS=1,NOccup
      XMuLoc(I)=XMuLoc(I)+FPsiB(IP,IQ,IR,IS)
     $ *OrbGrid(I,IP)*OrbGrid(I,IQ)*OrbGrid(I,IR)*OrbGrid(I,IS)
      EndDo
      EndDo
      EndDo
      EndDo
C
      XMuLoc(I)=SQRT(3.1415)/Two*XMuLoc(I)/OnTop(I)
C
      EndIf
C
      EndDo
      Print*, 'OnTop', norm2(OnTop)
      print*, 'XMuLoc',norm2(XMuLoc)
C
C     COMPUTE XMuMAT MATRIX in NO
C
      Do IP=1,NBasis
      Do IQ=1,IP
      XMuMAT(IP,IQ)=Zero
C
      ISym=MultpC(NSymNO(IP),NSymNO(IQ))
      If(ISym.Eq.1) Then
C
      Do I=1,NGrid
      XMuMAT(IP,IQ)=XMuMAT(IP,IQ)
     $ +OrbGrid(I,IP)*OrbGrid(I,IQ)*WGrid(I)*Exp(-XMuLoc(I))
      EndDo
C
      EndIf
      XMuMAT(IQ,IP)=XMuMAT(IP,IQ)
C
      EndDo
      EndDo
C
C     COMPUTE CBS CORRECTION
C
      Call PBECor(RhoGrid,Sigma,Zk,NGrid)
C
      SPi=SQRT(3.141592653589793)
      Const=Three/Two/SPi/(One-SQRT(Two))
C
      AvMU=Zero
      XEl=Zero
      CorrMD=Zero
      Do I=1,NGrid
C
      XMu=XMuLoc(I)
C
      AvMU=AvMU+XMu*RhoGrid(I)*WGrid(I)
      XEl=XEl+RhoGrid(I)*WGrid(I)
C
      Bet=Zero
      OnTopC=Zero
      If(XMu.Gt.1.D-8) OnTopC=OnTop(I)/(One+Two/SPi/XMu)
      If(OnTopC.Ne.Zero) Then
      Bet=Const*Zk(I)/OnTopC
      C=1.0
      CorrMD=CorrMD+Zk(I)/(C*One+Bet*XMu**3)*WGrid(I)
      EndIf
C
      EndDo
C
      AvMU=AvMU/XEl
      Write(6,'(/," XEl",F15.2/)') XEl
C
c      If(IFlCore.Eq.0) Then
c      Do I=1,NInAcCAS
c      Occ(I)=One
c      EndDo
c      Write(6,'(/,1X,"*** CBS CORRECTION EXCLUDING CORE ORBITALS ***")')
c      EndIf
C
      Write(6,'(/," CBS[DFT] correction, average Mu",
     $ F15.8,F15.2/)') CorrMD,AvMU
C
      Return
      End
C*End Subroutine LOC_MU_CBS

*Deck get_RDM2Occ
      Subroutine get_RDM2Occ(RDM2val,Occ,INActive,NAct,NOccup,NBasis)
C
C     read active triangle of 2-RDM (stored in rdm2.dat)
C     fill the entire RDM2val(NOccup^4) matrix
C
C     Comment: RDM2Act is 1212, then RDM2val is 1122
C
      implicit none
C
      integer,intent(in) :: INActive,NAct,NOccup,NBasis
      double precision,intent(in)  :: Occ(NBasis)
      double precision,intent(out) :: RDM2val(NOccup,NOccup, 
     $                                        NOccup,NOccup)

      integer :: i,j,k,l
      integer :: iunit,ios
      integer :: NRDM2Act
      integer :: Ind(NBasis)
      double precision :: VAL
      double precision, allocatable :: RDM2Act(:)
      integer,external :: NAddrRDM
      double precision,external :: FRDM2

C     index matrix
      Ind = 0
      do i=1,NAct
         Ind(INActive+i) = i
      enddo
C
      NRDM2Act = NAct**2*(NAct**2+1)/2
      Allocate (RDM2Act(NRDM2Act))
      RDM2Act(1:NRDM2Act)=0d0
C
      Open(newunit=iunit,file='rdm2.dat',status='old')
      Do
         read(iunit,*,iostat=ios) i,j,k,l,VAL
         if(ios/=0) exit
         RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
      Enddo
      Close(iunit)
C
      Do L=1,NOccup
      Do K=1,NOccup
      Do J=1,NOccup
      Do I=1,NOccup
      RDM2val(I,J,K,L)=FRDM2(I,K,J,L,RDM2Act,Occ,Ind,NAct,NBasis)
      Enddo
      Enddo
      Enddo
      Enddo
c     Print*, 'RDM2val = ', norm2(RDM2val)
C
      End
*     End Subroutine get_RDM2Occ

*Deck PairD_Vis
      Subroutine PairD_Vis(URe,UNOAO,Occ,NBasis)
C
C     Visualise electron pair density of HELIUM atom on a sphere of the 
C     radius 0.5 bohr
C     2-RDM is constructed from PMAT read from pmat.dat
C
C     Do not use symmetry!
C
      use read_external
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Three=3.0D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Allocatable ::  OrbGrid(:,:),RhoGrid(:),
     $ OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:),WGrid(:),
     $ RR(:,:)
C
      Dimension URe(NBasis,NBasis),Occ(NBasis),UNOAO(NBasis,NBasis),
     $ PMAT(NBasis,NBasis)
C
      Call molprogrid0(NGrid,NBasis)
      Write(6,'(/,1X,"The number of Grid Points = ",I8)')
     $ NGrid
C
      Allocate (WGrid(NGrid))
      Allocate (OrbGrid(NGrid,NBasis))
      Allocate (OrbXGrid(NGrid,NBasis))
      Allocate (OrbYGrid(NGrid,NBasis))
      Allocate (OrbZGrid(NGrid,NBasis))
      Allocate(RhoGrid(NGrid))
      Allocate(RR(3,NGrid))
C
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $ WGrid,UNOAO,NGrid,NBasis)
      Call molprogrid1(RR,NGrid)
C
      Open(10,file="pmat.dat")
      Read(10,*)X
      Do I=1,NBasis
      Do J=1,NBasis
      Read(10,*) I1,I2,PMAT(I,J)
      EndDo
      EndDo
      Close(10)
C
C     find the radius closest to 0.5 bohr
C
      Radius=Zero
      Err=1.0d5
      Do I1=1,NGrid
C
      R1=RR(1,I1)**2+RR(2,I1)**2+RR(3,I1)**2
      R1=SQRT(R1)
      If(Abs(R1-0.5D0).Lt.Err) Then
      Radius=R1
      Err=Abs(R1-0.5D0)
      EndIf
C
      EndDo
C
      Do I1=1,NGrid
C
      R1=RR(1,I1)**2+RR(2,I1)**2+RR(3,I1)**2
      R1=SQRT(R1)
C
      If(Abs(R1-Radius).Lt.1.D-2) Then
C     
      Do I2=1,NGrid
C
      R2=RR(1,I2)**2+RR(2,I2)**2+RR(3,I2)**2
      R2=SQRT(R2)
C
      If(Abs(R2-Radius).Lt.1.D-2) Then
C
      R12=(RR(1,I1)-RR(1,I2))**2+(RR(2,I1)-RR(2,I2))**2
     $   +(RR(3,I1)-RR(3,I2))**2
      R12=SQRT(R12)
C
      PairD=Zero
      Do IP=1,NBasis
      Do IQ=1,NBasis
      PairD=PairD
     $ +SQRT(Two)*PMAT(IP,IQ)
     $ *OrbGrid(I1,IP)*OrbGrid(I2,IQ)
      EndDo
      EndDo
      PairD=PairD*PairD
C
      Write(*,*)R12,PairD
C
      EndIf
      EndDo
C
      Stop
      EndIf
      EndDo
C
      Write(*,*)'Radius of the sphere =',Radius
C
      Return
      End

*Deck FistOrderSREne
      Subroutine FirstOrderSREne(ETwoSR,NBasis)
C
C     calculate <Ref | H^SR | Ref> = Gamma_prqs * (pq | rs)_sr
C     adopted from TwoEneChol
C
      Implicit Real*8 (A-H,O-Z)
C
      Include 'commons.inc'
C
      Integer NBasis
      Double Precision ETwoSR
C
C     LOCAL ARRAYS
C
      Integer INActive,NAct,NOccup
      Integer iloop,nloop,off
      Integer dimFO,iBatch,BatchSize
      Integer :: Ind(NBasis)
      Double Precision, Allocatable :: FF(:,:), FFtr(:,:)
      Double Precision, Allocatable :: FFErf(:,:), FFErftr(:,:)
      Double Precision, Allocatable :: TmpTr(:,:), TmpErfTr(:,:)
      Double Precision, Allocatable :: RDM2val(:,:,:,:)
      Double Precision, Allocatable :: ints(:,:),work1(:,:)
      Double precision :: AuxICoeff(2,2,2,2),AuxVal,AuxIVal
      Double precision :: AuxActCoeff(3,3,3,3),AuxInactCoeff(3,3,3,3)
c
      Parameter(MaxBatchSize = 120)
C
C     SET DIMENSIONS
      NAct=NAcCAS
      INActive=NInAcCAS
      NOccup=NAct+INActive
c
c     INDEX STRONGLY OCCUPIED (SO) ORBITALS
c     (includes inactive)
      Ind = 1
      Do I=1,NStronglyOccOrb
      Ind(I) = 2
      EndDo
C     Auxiliary matrices
C     act:   IGem(i)=1 Igem(i)=Igem(j)=Igem(k)=Igem(l)
C     inact: IGem(i)=2 Igem(i)=Igem(j)=Igem(k)=Igem(l)
      Do l=1,3
         Do k=1,3
            do j=1,3
               Do i=1,3
                  If((i==2).and.(i==j).and.(j==k).and.(k==l)) then
                     AuxActCoeff(i,j,k,l) = 1
                  Else
                     AuxActCoeff(i,j,k,l) = 0
                  Endif
               Enddo
            enddo
         Enddo
      Enddo
      Do l=1,3
         Do k=1,3
            do j=1,3
               Do i=1,3
                  If((i==1).and.(i==j).and.(j==k).and.(k==l)) then
                     AuxInActCoeff(i,j,k,l) = 1
                  Else
                     AuxInActCoeff(i,j,k,l) = 0
                  Endif
               Enddo
            enddo
         Enddo
      Enddo
c
      Do l=1,2
         Do k=1,2
            do j=1,2
               Do i=1,2
                  If((i==2).and.(i==j).and.(j==k).and.(k==l)) then
                     AuxICoeff(i,j,k,l) = 0
                  Else
                     AuxICoeff(i,j,k,l) = 1
                  Endif
               Enddo
            enddo
         Enddo
      Enddo
C
c     Load Cholesky Vecs
      Open(newunit=iunit,file='cholvecs',form='unformatted')
      Read(iunit) NCholesky
      Allocate(FF(NCholesky,NBasis**2))
      Read(iunit) FF
      Close(iunit)
c
      Open(newunit=iunit,file='chol1vFR',form='unformatted')
      Read(iunit) NCholesky
      Allocate(FFTr(NCholesky,NBasis**2))
      Read(iunit) FFTr
      Close(iunit)
c
      Open(newunit=iunit,file='cholvErf',form='unformatted')
      Read(iunit) NCholErf
      Allocate(FFErf(NCholErf,NBasis**2))
      Read(iunit) FFErf
      Close(iunit)
c
      Open(newunit=iunit,file='chol1vLR',form='unformatted')
      Read(iunit) NCholErf
      Allocate(FFErfTr(NCholErf,NBasis**2))
      Read(iunit) FFErfTr
      Close(iunit)
C
c     print*, 'NCholesk =',NCholesky
c     print*, 'NCholErf =',NCholErf
      Allocate(TmpErfTr(NCholErf,NBasis*NOccup))
      ! p*q : LR
      Do iq=1,NOccup
         Do ip=1,NBasis
            TmpErfTr(:,ip+(iq-1)*NBasis)=FFErfTr(:,iq+(ip-1)*NBasis)
         EndDo
      EndDo
      Allocate(TmpTr(NCholesky,NBasis*NOccup))
      ! p*q : FR
      Do iq=1,NOccup
         Do ip=1,NBasis
            TmpTr(:,ip+(iq-1)*NBasis)=FFTr(:,iq+(ip-1)*NBasis)
         EndDo
      Enddo
C
      Allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))
      CAll get_RDM2Occ(RDM2val,Occ,INActive,NAct,NOccup,NBasis)
C
      dimFO = NBasis*NOccup
      nloop = (dimFO - 1) / MaxBatchSize + 1
C
      Allocate(ints(NBasis,NBasis))
      Allocate(work1(dimFO,MaxBatchSize))
C
      ETwoSR=0
C     EXCHANGE LOOP (FO|FO), use only (OO|OO)
      off = 0
      k   = 0
      l   = 1
      Do iloop=1,nloop

      ! batch size for each iloop; last one is smaller
      BatchSize = min(MaxBatchSize,dimFO-off)
C
      ! assemble (FO|BatchSize)_Sr batch from CholVecs
      ! (pq*|rs) : SR=-LR+FR
      call dgemm('T','N',dimFO,BatchSize,NCholErf,
     $         -0.25d0,FFErfTr,NCholErf,
     $          FFErf(:,off+1:off+BatchSize),NCholErf,0d0,work1,dimFO)
      call dgemm('T','N',dimFO,BatchSize,NCholesky,
     $          0.25d0,FFTr,NCholesky,
     $          FF(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
      ! (pq|rs*): SR=-LR+FR
      call dgemm('T','N',dimFO,BatchSize,NCholErf,
     $         -0.25d0,FFErf,NCholErf,
     $          FFErfTr(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
      call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,FF,NCholesky,
     $          FFTr(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
      ! (p*q|rs)
      call dgemm('T','N',dimFO,BatchSize,NCholErf,
     $         -0.25d0,TmpErfTr,NCholErf,
     $          FFErf(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
      call dgemm('T','N',dimFO,BatchSize,NCholesky,
     $          0.25d0,TmpTr,NCholesky,
     $          FF(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
      ! (pq|r*s)
      call dgemm('T','N',dimFO,BatchSize,NCholErf,
     $         -0.25d0,FFErf,NCholErf,
     $         TmpErfTr(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
      call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,FF,NCholesky,
     $           TmpTr(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
C
      Do iBatch=1,BatchSize

      k = k + 1
      if(k>NBasis) then
         k = 1
         l = l + 1
      endif

      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i,iBatch)
         enddo
      enddo
C
      if(k>NOccup) cycle
C
        Do j=1,NOccup
           Do i=1,NOccup
              AuxVal  = AuxActCoeff(IGem(i),IGem(j),IGem(k),IGem(k))
              AuxIVal = AuxICoeff(Ind(i),Ind(j),Ind(k),Ind(l))
              ETwoSR  = ETwoSR
     $                + AuxVal*AuxIVal*RDM2val(i,j,k,l)*ints(i,j)

              If(NInAcCAS.Eq.0) Then
                AuxVal  = AuxInActCoeff(IGem(i),IGem(j),IGem(k),IGem(k))
                ETwoSR  = ETwoSR
     $                  + AuxVal*AuxIVal*RDM2val(i,j,k,l)*ints(i,j)

              EndIf
          EndDo
        EndDo
C
      EndDo
C
      off = off + MaxBatchSize
C
      EndDo
c
c     print*, 'ETwoSr = ', ETwoSR
C
      Deallocate(ints,work1)
      Deallocate(FF,FFtr,FFErf,FFErfTr)
      Deallocate(TmpTr,TmpErfTr)
      Deallocate(RDM2val)
C
      End
C     End Subroutine FirstOrderSREne

