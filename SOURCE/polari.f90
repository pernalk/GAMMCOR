!#define POLARI_DEBUG 6

module polari

use print_units
use read_external
use tran
use abfofo
use gammcor_integrals

implicit none

contains

subroutine Polariz(FreqOm,ECASSCF,UNOAO,XOne,URe,Occ,&
   IGem,NAct,INActive,NElecBEmb,NELE,NFreq,&
   NBasis,NInte1,IndAux,&
   IndN,IndX,NDimX,ICholesky,IFunSR,IFunSrKer,IntIdx, &
   MemVal,MemType)
!
!  to-do: 
!   -- it wpould be more elegant to it via types
!
! Returns dynamic polarizability tensor for a given frequency
! find C(omega) by inversion
!
! IntIdx = 1 (Molpro), 2 (Dalton), 0 (default)

implicit none
integer,intent(in) :: NBasis,NInte1,NDimX
integer,intent(in) :: NAct,INActive,NElecBEmb,NELE
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),IndAux(NBasis),IGem(NBasis)
integer,intent(in) :: ICholesky,IFunSR,IFunSRKer,IntIdx
integer,intent(in) :: MemVal,MemType
integer,intent(in) :: NFreq
double precision,intent(in) :: FreqOm(NFreq),Occ(NBasis)
double precision,intent(in) :: UNOAO(NBasis,NBasis),URe(NBasis,NBasis)
double precision,intent(inout) :: XONe(NInte1)
double precision,intent(out)   :: ECASSCF

double precision :: CICoef(NBasis),ipiv(NDimX)
double precision :: UAux(NBasis,NBasis) 
double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
double precision :: DipCX(NDimX),DipCY(NDimX),DipCZ(NDimX)
double precision :: ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),AIN(NDimX*NDimX),CMAT(NDimX*NDimX)
double precision :: AXX,AYX,AXY,AZX,AXZ,AYY,AZY,AYZ,AZZ,Om,Alpha
double precision :: AISO
!
double precision :: VsrKS(NBasis,NBasis),Jsr(NBasis,NBasis)
double precision :: JsrTR(NInte1),VsrTR(NInte1)
double precision :: esrDFT(3)

logical :: doGGA
double precision,allocatable :: SRKer(:)
double precision,allocatable :: OrbGrid(:,:),WGrid(:)
double precision,allocatable :: OrbXGrid(:,:),OrbYGrid(:,:),OrbZGrid(:,:)

integer :: NOccup,NGrid
integer :: I,J,II,IJ,im
integer :: IStart,IOrb,inf
integer :: num0,num1
integer :: NSymBas(8),NSymOrb(8)
integer :: NSym,NSymNO(NBasis),NumOSym(15),MultpC(15,15)
integer(8) :: MemSrtSize

character(:),allocatable :: twojfile,twokfile
double precision,external :: ddot

character(*),Parameter :: griddalfile='dftgrid.dat'

!Om=FreqOm
NOccup=NAct+INActive

If (IntIdx == 1) then ! Molpro
   Call ComputeDipoleMom(UNOAO,Occ,'DIP','AOONEINT.mol',NOccup,NBasis)
   Call ReadDip(DipX,DipY,DipZ,UNOAO,'DIP',NBasis)
ElseIf (IntIdx == 2) then ! Dalton
   Call ComputeDipoleMom(UNOAO,Occ,'AOPROPER','AOONEINT',NOccup,NBasis)
   Call ReadDip(DipX,DipY,DipZ,UNOAO,'AOPROPER',NBasis)
ElseIf (IntIdx == 0) then ! Dalton
   Stop "Unknown Interface in Polariz!"
EndIf

If(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
   If (IntIdx == 2) then ! dalton
      call read_vKS_dalton(VsrKS,'dftSRfile.dat',NBasis)
      call read_Jsr_dalton(Jsr,'dftSRfile.dat',NBasis)
      call read_esrDFT_dalton(esrDFT,'dftSRfile.dat')

      write(6,'(/1x,a)') 'SR KS potential read from dftSRfile.dat'
      write(6,'(1x,a)')  'SR Coulomb ints read from dftSRfile.dat'

      Call sq_to_triang2(Jsr,JsrTR,NBasis)
      Call sq_to_triang2(VsrKS,VsrTR,NBasis)
      Call MatTr(JsrTR,UNOAO,NBasis)
      Call MatTr(VsrTR,UNOAO,NBasis)

      XOne = XOne + VsrTR + JsrTR
   ElseIf (IntIdx == 1) then ! molpro
      Stop "Molpro intefrace not ready in Polari!"
   EndIf
EndIf

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

Do IJ=1,NDimX
  I=IndN(1,IJ)
  J=IndN(2,IJ)
  DipCX(IndX(IJ))=(CICoef(I)+CICoef(J))*DipX(I,J)
  DipCY(IndX(IJ))=(CICoef(I)+CICoef(J))*DipY(I,J)
  DipCZ(IndX(IJ))=(CICoef(I)+CICoef(J))*DipZ(I,J)
Enddo

If (IFunSR == 0) Then
   twojfile = 'FFOO'
   twokfile = 'FOFO'
ElseIf (IFunSR > 0) Then
   twojfile = 'FFOOERF'
   twokfile = 'FOFOERF'
EndIf

if(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
if(IFunSRKer.Eq.1) then
  if (IntIdx /= 2) Stop "DFT only with DALTON in Polari!"
  if (ICholesky == 1) Stop "DFT/Cholesky not ready in Polari!"

  call daltongrid0(NGrid,griddalfile,doGGA,NBasis)
  allocate(WGrid(NGrid),OrbGrid(NGrid,NBasis))
  UAux=transpose(UNOAO)
  if (doGGA) then
     write(lout,'(/1x,a)') 'GRID TYPE: GGA'
     allocate(OrbXGrid(NGrid,NBasis))
     allocate(OrbYGrid(NGrid,NBasis))
     allocate(OrbZGrid(NGrid,NBasis))
     call daltongrid_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid, &
                         WGrid,NGrid,griddalfile,NBasis)
     call daltongrid_tran_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid, &
                         UAux,0,NGrid,NBasis,.false.)
  else
     write(lout,'(/1x,a)') 'GRID TYPE: LDA'
     print*, 'Grid ...',NGrid
     call daltongrid_lda(OrbGrid,WGrid,NGrid,griddalfile,NBasis)
     call daltongrid_tran_lda(OrbGrid,UAux,0,NGrid,NBasis,.false.)
  endif

 ! get symmetry of orbitals
  call read_sym_dalton(NSym,NSymBas,NumOSym,'SIRIUS.RST','BASINFO ')

  NSymNO(1:NBasis)=0
  IStart=0
  do I=1,NSym
     do J=IStart+1,IStart+NumOSym(I)
        do IOrb=1,NBasis
           if(Abs(UNOAO(IOrb,J)).Gt.1.D-1) NSymNO(IOrb)=I
        enddo
     enddo
     IStart=IStart+NumOSym(I)
  enddo

!#if POLARI_DEBUG > 5
  print*, 'Polariz: NSym =',NSym
  print*, 'UNOAO orbitals ',norm2(UNOAO)
  do j=1,NBasis
     write(6,'(*(f13.8))') (UNOAO(i,j),i=1,NBasis)
  enddo
  print*, 'NumOSym '
  do i=1,NSym
     print*, i,NumOSym(i)
  enddo

!#endif

! checking
  do I=1,NSym
  II=0
  do IOrb=1,NBasis
  if(NSymNO(IOrb).Eq.I) II=II+1
  enddo
  If(II.Ne.NumOSym(I)) Write(*,*) 'In Polariz: Symmetry of NO cannot be established!'
  EndDo
  if(NSym.Eq.1) Then
     MultpC(1,1)=1
  else
     do I=1,NSym
        do J=1,I
           MultpC(I,J)=IEOR(I-1,J-1)+1
           MultpC(J,I)=MultpC(I,J)
        enddo
     enddo
  endif

 ! get and transform 2-el integrals
 MemSrtSize=MemVal*1024_8**MemType
 call readtwoint(NBasis,1,'AOTWOINT','AOTWOSORT',MemSrtSize)

 write(6,'(" Transforming full two-electron integrals ...")')
 !PREPARE POINTERS: NOccup=num0+num1
 Call prepare_nums(Occ,Num0,Num1,NBasis)
 Call tran4_gen(NBasis,  &
              Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)), &
              Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)), &
              NBasis,UAux, &
              NBasis,UAux, &
              'FFOO','AOTWOSORT')
 Call tran4_gen(NBasis, &
            NBasis,UAux, &
            Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)), &
            NBasis,UAux,&
            Num0+Num1,UAux(1:NBasis,1:(Num0+Num1)), &
            'FOFO','AOTWOSORT')

endif ! IFunSRKer
endif

Alpha=1.0
Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NElecBEmb,&
              NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,ICholesky,0,Alpha,.false.)

!print*, 'ABMIN-before  =', norm2(ABMIN)

if(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) then
   if(IFunSRKer.Eq.1) then
      write(6,'(/," *** Generating a sr-kernel on a grid ***")')
      allocate (SRKer(NGrid))
      call GetKerNPT(SRKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,NBasis,NGrid)
      call ModABMin_FOFO(Occ,SRKer,WGrid,OrbGrid,ABMin, &
                         MultpC,NSymNO, &
                         IndN,IndX,NDimX,NGrid,NBasis, &
                         num0,num1, &
                         'FOFO','FOFOERF',.false.)
   endif
endif

!Print*, 'XOne = ', norm2(XOne)
!Print*, 'ABPLUS =', norm2(ABPLUS)
!Print*, 'ABMIN  =', norm2(ABMIN)
!print*, 'NAct  = ', NAct
!print*, 'INAct = ', INActive
!print*, 'NDimX = ', NDimX

do im=1,NFreq
Om = FreqOm(im)

AIN=0d0
Do I=1,NDimX
    AIN((I-1)*NDimX+I)=1.0
EndDo
!  ABPLUS*ABMIN - 1 Om^2
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,&
           ABMIN,NDimX,-Om**2,AIN,NDimX)
CMAT=0.5d0*ABPLUS
Call dgesv(NDimX,NDimX,AIN,NDimX,ipiv,CMAT,NDimX,inf)

! contract CMAT with dipole moment vectors
Call dgemv('N',NDimX,NDimX,1.d0,CMAT,NDimx,DipCX,1,0.d0,ipiv,1)

AYX=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZX=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AXX=8.d0*ddot(NDimx,DipCX,1,ipiv,1)

Call dgemv('N',NDimX,NDimX,1.d0,CMAT,NDimx,DipCY,1,0.d0,ipiv,1)
AXY=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AZY=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AYY=8.d0*ddot(NDimx,DipCY,1,ipiv,1)

Call dgemv('N',NDimX,NDimX,1.d0,CMAT,NDimx,DipCZ,1,0.d0,ipiv,1)
AXZ=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AYZ=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZZ=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AISO=(AXX+AYY+AZZ)/3d0

Write(6,'(/,X,''Polarizability tensor for frequency '',F8.4)') Om
Write(6,'(/,X,''XX   XY   XZ  '',3F15.8)') AXX, AXY, AXZ
Write(6,'(X,''YX   YY   YZ  '',3F15.8)') AYX, AYY, AYZ
Write(6,'(X,''ZX   ZY   ZZ  '',3F15.8,1/)') AZX, AZY, AZZ
Write(6,'(X,''Isotropic Polarizability'',F15.8)') AISO

enddo ! Om

If(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
  Call EneMCsrDFT(ECASSCF,esrDFT(1),Occ,XOne,JsrTR,VsrTR,NInte1,NBasis)
EndIf

! clean-up
call delfile('FFOO')
call delfile('FOFO')
If(IFunSR.Eq.1.Or.IFunSR.Eq.2.Or.IFunSR.Eq.4) Then
   call delfile('FFOOERF')
   call delfile('FOFOERF')
endif

end subroutine Polariz

subroutine EneMCsrDFT(ECASSCF,EnxcSR,Occ,XOne,Jsr,VsrKS,NInte1,NBasis)
!
! calculate MC-srDFR energy
! E = T+V_ne+VKS^sr+J^sr+... ?

double precision,intent(in) :: ECASSCF,EnxcSR
double precision,intent(in) :: Occ(NBasis)
double precision,intent(in) :: XOne(NInte1),Jsr(NInte1),VsrKS(NInte1)
integer,intent(in) :: NInte1,NBasis

integer :: i,ii,ione
integer :: NSym,NBas(1:NBasis)
double precision ::  EOne,ENuc,EnHSR,EnSR,XVSR
logical :: exione

EOne  = 0d0
EnHSR = 0d0
XVSR  = 0d0
do i=1,NBasis
   ii = (i*(i+1))/2
   EOne  = EOne  + 2d0*Occ(i)*XOne(II)
   XVSR  = XVSR  + 2d0*Occ(i)*(VsrKS(ii)+Jsr(ii))
   EnHSR = EnHSR + Occ(i)*Jsr(ii)
enddo

inquire(File='AOONEINT',EXIST=exione)
if(exione) Then
 open(newunit=ione,File='AOONEINT',access='SEQUENTIAL',Form='UNFORMATTED',Status='OLD')
 read(ione)
 read(ione) NSym,NBas(1:NSym),ENuc
 close(ione)
endif

EnSR = EnxcSR + EnHSR
!Print*, 'ECASSCF', ECASSCF
!print*, 'XVSR', XVSR
write(6,'(1x,"Nuclear repulsion:",T40,F15.8)') ENuc
write(6,'(1x,"sr Coulomb energy:",T40,F15.8)')   EnHSR
write(6,'(1x,"One-electron energy:",T40,F15.8)') EOne-XVSR
write(6,'(1X,"lrCASSCF+ENuc Energy",T40,F15.8)') ECASSCF - XVSR + ENuc
write(6,'(1X,"Total lrCASSCF+ENuc+srDF Energy",T40,F15.8,/)') ECASSCF - XVSR + EnSR + ENuc

end subroutine EneMCsrDFT

subroutine PolarizAl(FreqOm,ECASSCF,UNOAO,XOne,URe,Occ,&
   IGem,NAct,INActive,NElecBEmb,NELE,NFreq,NBasis,NInte1,NGem,IndAux,&
   IndN,IndX,NDimX,BasisSet,ICholesky,Max_Cn,IntIdx)
!
! Returns dynamic polarizability tensor for a given frequency FreqOm
! find C(omega) by expanding around Alpha=0 with a tolerance Eps or
! up to maximal order Max_Cn
!

implicit none
integer,intent(in) :: NBasis,NInte1,NGem,NDimX,Max_cn
integer,intent(in) :: NAct,INActive,NElecBEmb,NELE
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),IndAux(NBasis),IGem(NBasis)
integer,intent(in) :: IntIdx
character(6) :: Source
character(*) :: BasisSet
!double precision,intent(in) :: FreqOm
integer,intent(in) :: NFreq
double precision,intent(in) :: FreqOm(NFreq)
double precision,intent(in) :: Occ(NBasis),XOne(NInte1)
double precision,intent(in) :: UNOAO(NBasis,NBasis),URe(NBasis,NBasis)
double precision :: ECASSCF

integer :: i,j,ij,im,inf,ICholesky,NOccup
double precision :: CICoef(NBasis)
double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis)
double precision :: DipCX(NDimX),DipCY(NDimX),DipCZ(NDimX)
double precision :: ipiv(NDimX)
double precision :: AXX,AYX,AXY,AZX,AXZ,AYY,AZY,AYZ,AZZ,Om,ddot,Alpha
double precision :: AISO
character(:),allocatable ::  XYZPath
character(:),allocatable :: twojfile,twokfile

!Om=FreqOm

! set source
If (IntIdx == 1) then
 Source='MOLPRO'
ElseIf (IntIdx == 2) then
 Source='DALTON'
ElseIf (IntIdx == 0) then
   Stop "Unknown Interface in PolarizAl!"
EndIf

XYZPath="./input.inp"
!Call DipMomOTF_ao(DIpX,DIpY,DipZ,BasisSetPath,XYZPath,IUnits,'MOLPRO')
!Call CompDipMomOTF(AOBasis,System,CAONO,Occ,DipX,DipY,DipZ,NBasis,NBasis)
!stop "stop here..."

NOccup=NAct+INActive
if(trim(Source)=='MOLPRO') then
   Call ComputeDipoleMom(UNOAO,Occ,'DIP','AOONEINT.mol',NOccup,NBasis)
   Call ReadDip(DipX,DipY,DipZ,UNOAO,'DIP',NBasis)
elseif(trim(Source)=='DALTON') then
   Call ComputeDipoleMom(UNOAO,Occ,'AOPROPER','AOONEINT',NOccup,NBasis)
   Call ReadDip(DipX,DipY,DipZ,UNOAO,'AOPROPER',NBasis)
EndIf

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

Do IJ=1,NDimX
  I=IndN(1,IJ)
  J=IndN(2,IJ)
  DipCX(IndX(IJ))=(CICoef(I)+CICoef(J))*DipX(I,J)
  DipCY(IndX(IJ))=(CICoef(I)+CICoef(J))*DipY(I,J)
  DipCZ(IndX(IJ))=(CICoef(I)+CICoef(J))*DipZ(I,J)
Enddo

do im=1,NFreq

Om = FreqOm(im)

Call CFREQPROJ(ipiv,Om,DipCX,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NElecBEmb,&
   NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX)

AYX=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZX=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AXX=8.d0*ddot(NDimx,DipCX,1,ipiv,1)

Call CFREQPROJ(ipiv,Om,DipCY,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NElecBEmb,&
   NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX)

AXY=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AZY=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AYY=8.d0*ddot(NDimx,DipCY,1,ipiv,1)

Call CFREQPROJ(ipiv,Om,DipCZ,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NElecBEmb,&
   NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX)
AXZ=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AYZ=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZZ=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AISO=(AXX+AYY+AZZ)/3d0

Write(6,'(/,X,''Polarizability tensor for frequency '',F8.4)') Om
Write(6,'(/,X,''XX   XY   XZ  '',3F15.8)') AXX, AXY, AXZ
Write(6,'(X,''YX   YY   YZ  '',3F15.8)') AYX, AYY, AYZ
Write(6,'(X,''ZX   ZY   ZZ  '',3F15.8,/)') AZX, AZY, AZZ
Write(6,'(X,''Isotropic Polarizability'',F15.8)') AISO

enddo ! Om

end subroutine PolarizAl

end module polari
