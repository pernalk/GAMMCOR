subroutine Project_DChol(PMat,IndN,NBasis,NDimX)
!
! compute the projector matrix PMat, used optionally in ACFREQ
!
implicit none

integer,intent(in) :: NBasis,NDimX,IndN(2,NDimX)
double precision :: PMat(NDimX,NDimX)

integer :: iunit
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),N,IGL,inf1,inf2,info
double precision,allocatable :: work(:),ints(:,:)
Real*8, Allocatable :: MatFF(:,:),work2(:,:),work3(:),work1(:,:),work4(:),work5(:,:),work6(:,:)
integer :: NCholesky,lwork

open(newunit=iunit,file='cholvecs',form='unformatted')
read(iunit) NCholesky
allocate(MatFF(NCholesky,NBasis**2))
read(iunit) MatFF
close(iunit)

print*,'NCholesky',NCholesky

allocate(work2(NCholesky,NCholesky))
allocate(work1(NCholesky,NDimX))

work1 = 0
do j=1,NDimX
ir=IndN(1,j)
is=IndN(2,j)
irs = is+(ir-1)*NBasis
do i=1,NCholesky
    work1(i,j) = MatFF(i,irs)
enddo
enddo
deallocate(MatFF)

allocate(MatFF(NCholesky,NDimX))
MatFF=work1
deallocate(work1)
! D.D^T
call dgemm('N','T',NCholesky,NCholesky,NDimX,1.0d0,MatFF,NCholesky,MatFF,NCholesky,0d0,work2,NCholesky)
!call dgemm('N','T',NCholesky,NCholesky,NBasis**2,1.0d0,MatFF,NCholesky,MatFF,NCholesky,0d0,work2,NCholesky)

allocate(work3(NCholesky),work4(NBasis),work1(NCholesky,NCholesky))
work1=work2
call Diag8(work2,NCholesky,NCholesky,work3,work4)
do i=1,NCholesky
if(abs(work3(i)).lt.1d-10) then
!      print*, i, work3(i)
   work3(i) = 0
else
   work3(i)=1./work3(i)
endif
enddo
allocate(work5(NCholesky,NCholesky),work6(NCholesky,NCholesky))
work5=0
do i=1,NCholesky
work5(i,i)=work3(i)
enddo
call dgemm('N','N',NCholesky,NCholesky,NCholesky,1.0d0,work5,NCholesky,work2,NCholesky,0d0,work6,NCholesky)
call dgemm('T','N',NCholesky,NCholesky,NCholesky,1.0d0,work2,NCholesky,work6,NCholesky,0d0,work5,NCholesky)

! checking the inverse
call dgemm('N','N',NCholesky,NCholesky,NCholesky,1.0d0,work5,NCholesky,work1,NCholesky,0d0,work6,NCholesky)
print*,'norm of D.D^T.[D.D^T]^-1',norm2(work6(1:NCholesky,1:NCholesky)),'Sqrt(NCholesky)',SQRT(Float(NCholesky))

work2=work5
deallocate(work4,work3,work5,work6,work1)

! D^T.(D.D^T)-1 = work1
!allocate(work1(NBasis**2,NCholesky))
allocate(work1(NDimX,NCholesky))
!call dgemm('T','N',NBasis**2,NCholesky,NCholesky,1.0d0,MatFF,NCholesky,work2,NCholesky,0d0,work1,NBasis**2)
call dgemm('T','N',NDimX,NCholesky,NCholesky,1.0d0,MatFF,NCholesky,work2,NCholesky,0d0,work1,NDimX)
! construct P = work1.D
deallocate(work2)
!allocate(work2(NBasis**2,NBasis**2))
allocate(work2(NDimX,NDimX))
!call dgemm('N','N',NBasis**2,NBasis**2,NCholesky,1.0d0,work1,NBasis**2,MatFF,NCholesky,0d0,work2,NBasis**2)
call dgemm('N','N',NDimX,NDimX,NCholesky,1.0d0,work1,NDimX,MatFF,NCholesky,0d0,work2,NDimX)
deallocate(work1)

!print*, 'PMat',norm2(work2(1:NBasis**2,1:NBasis**2))
print*, 'PMat',norm2(work2(1:NDimX,1:NDimX))
PMAT=work2
deallocate(MatFF)

return

deallocate(work1)
end subroutine Project_DChol

subroutine WIter_D12Chol(ECorr,AC1,Max_Cn,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)
!
!  AC energy cacluation using CHOLESKY VECTORS:
!  (1) expanding AC integrand in alpha around alpha=0, up to Max_Cn order
!  (2) finding C^(n)[omega] and (3) omega integration
!  Max_Cn is passed in input.inp (if missed, the dafault value 3 is used)
!
!  A difference with WIter_DChol: no need to compute COMTildeAct
!
use abfofo
use systemdef
! only to use Y01CAS_FOFO
!use ab0fofo
use sapt_utils

implicit none
integer,intent(in) :: AC1,NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                      IGem(NBasis)
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)

double precision :: ACAlpha
double precision :: ECorr,ECorrAct,EGOne(NGem)
double precision :: XFreq(100),WFreq(100)

integer :: ICholesky
integer :: NCholesky
integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,OmI,XNorm0,XNorm1,ErrMax
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: irdm2

double precision, allocatable :: DChol(:,:),DCholT(:,:),DCholAct(:,:),DCholActT(:,:),WorkD(:,:)
double precision, allocatable :: APLUS0Tilde(:), APLUS1Tilde(:), A1(:), &
                                 COMTilde(:),ABPLUS0(:),ABMIN0(:),ABPLUS1(:),ABMIN1(:), &
                                 C0Tilde(:),C1Tilde(:),C2Tilde(:), &
                                 WORK0(:),WORK1(:)

integer :: nblk
type(EblockData) :: A0blockIV,LambdaIV
type(EblockData),allocatable :: A0block(:),Lambda(:)

interface
subroutine read_D12_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
   double precision, allocatable, intent(out) :: DChol(:,:), DCholAct(:,:)
   integer :: NCholesky
   integer, intent(in) :: NDimX, NBasis, IndN(2,NDimX), IndAux(NBasis)
   double precision, intent(in) :: Occ(NBasis)
end subroutine read_D12_array
end interface

! Get DChol & DCholAct
call read_D12_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
DCholT = transpose(DChol)
DCholActT = transpose(DCholAct)
! ==========================================================================

NGrid=18

If(AC1.Eq.0) Then
Write (6,'(/,X,''AC calculation through W_AC expansion around Alpha=0, Omega Grid = '',I3,&
   '' and max order in C expansion = '',I3,/)') NGrid,Max_Cn
Else
Write (6,'(/,X,''AC1 calculation through W_AC expansion around Alpha=0, Omega Grid = '',I3,&
   '' and max order in C expansion = '',I3,/)') NGrid,Max_Cn
EndIf

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

ICholesky = 1

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

allocate(ABPLUS1(NDimX*NDimX),ABMIN1(NDimX*NDimX))

if(NAct==1) then
  ! active-virtual block
  nblk = NBasis - NAct - INActive
else
  nblk = 1 + NBasis - NAct
endif
!print*, 'outer nblk',nblk

allocate(A0block(nblk))
! AC0BLOCK with ver=0 stores A-(0) and A+(0) matrices
!                            in X and Y, respectively
Call AC0BLOCK(Occ,URe,XOne, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1,'FFOO','FOFO', &
     ICholesky,A0BlockIV,A0Block,nblk,0,'DUMMY',0)

! get AB1PLUS and AB1MIN
ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS1,ABMIN1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,ICholesky,ACAlpha,.false.)

Call sq_symmetrize(ABPLUS1,NDimX)
Call sq_symmetrize(ABMIN1,NDimX)

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1,A0Block,A0BlockIV,-1d0,.false.,nblk,NDimX)
call add_blk_right(ABMIN1, A0Block,A0BlockIV,-1d0,.true., nblk,NDimX)
!print*, 'add_blk_right: ABPLUS1',norm2(ABPLUS1)
!print*, 'add_blk_right: ABMIN1 ',norm2(ABMIN1)

!Calc: A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
allocate(A1(NDimX*NDimX))
call ABPM_HALFTRAN_GEN_L(ABMIN1, A1,0.0d0,A0Block,A0BlockIV,nblk,NDimX,NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1,A1,1.0d0,A0Block,A0BlockIV,nblk,NDimX,NDimX,'X')
!print*, 'A1',norm2(A1)

EGOne(1)=ECASSCF

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(APLUS0Tilde(NDimX*NCholesky))
call ABPM_HALFTRAN_GEN_L(DCholT,APLUS0Tilde,0.0d0,A0Block,A0BlockIV,nblk,NDimX,NCholesky,'Y')

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(APLUS1Tilde(NDimX*NCholesky))
Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholT,NDimX,0.0d0,APLUS1Tilde,NDimX)

deallocate(A0block)
deallocate(A0BlockIV%vec,A0BlockIV%pos)

Call FreqGrid(XFreq,WFreq,NGrid)

! Calc: A0
allocate(A0block(nblk))
! ver=1: store A+(0).A-(0) in blocks 
Call AC0BLOCK(Occ,URe,XOne, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1,'FFOO','FOFO', &
     ICholesky,A0BlockIV,A0Block,nblk,1,'A0BLK',0)

allocate(COMTilde(NDimX*NCholesky))
COMTilde=0.0

allocate(C0Tilde(NDimX*NCholesky),C1Tilde(NDimX*NCholesky),C2Tilde(NDimX*NCholesky),WORK0(NDimX*NCholesky))
allocate(WORK1(NDimX*NCholesky))
allocate(Lambda(nblk))
associate(A => A0BlockIV, L => LambdaIV)
  L%n = A%n
  L%l1 = A%l1
  L%l2 = A%l2
  allocate(L%pos(L%n),L%vec(L%n))
end associate

Do IGL=1,NGrid
   OmI=XFreq(IGL)
   WFact=4.D0/PI*WFreq(IGL)

!  Calc: LAMBDA=(A0+Om^2)^-1
   Call INV_AC0BLK(OmI**2,Lambda,LambdaIV,A0Block,A0BlockIV,nblk,NDimX)

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
   Call ABPM_HALFTRAN_GEN_L(APLUS0Tilde,C0Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')
   C0Tilde = 0.5d0*C0Tilde

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,A1,NDimX,C0Tilde,NDimX,0.0d0,WORK0,NDimX)
   WORK0 = 0.5d0*APLUS1Tilde - WORK0
   Call ABPM_HALFTRAN_GEN_L(WORK0,C1Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')

   COMTilde=COMTilde+WFact*0.5d0*C1Tilde

   XNorm0=1.0d5
   XFactorial=1
   Do N=2,Max_Cn
       XFactorial=XFactorial*N
       XN1=-N
       XN2=-N*(N-1)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,XN2,ABMIN1,NDimX,C0Tilde,NDimX,0.0d0,WORK1,NDimX)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,ABPLUS1,NDimX,WORK1,NDimX,0.0d0,WORK0,NDimX)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,XN1,A1,NDimX,C1Tilde,NDimX,1.0d0,WORK0,NDimX)
       Call ABPM_HALFTRAN_GEN_L(WORK0,C2Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NCholesky,'X')
       FF=WFact/XFactorial/(N+1)
       If(AC1.Eq.1) FF=WFact/XFactorial/2.D0
       XNorm1=Norm2(FF*C2Tilde)
       Write(6,'(X,"Order (n), |Delta_C|",I3,E14.4)')N,XNorm1
       If(XNorm1.Lt.ErrMax) Exit
       If(N.Gt.3.And.XNorm1.Gt.XNorm0) Then
           Write(6,'(X,"Divergence detected. Expansion of C terminated at order ",I3,3F10.4)')N-1
!          Write(6,'(X,"Divergence detected. Continue up to order Max_Cn",I3,3F10.4)')N
           Exit
       EndIf
       XNorm0=XNorm1
       COMTilde=COMTilde+FF*C2Tilde
       C0Tilde=C1Tilde
       C1Tilde=C2Tilde
   EndDo

   Write(6,'(X,"Omega, |C|",I3,2F10.4)')IGL,OmI,Norm2(COMTilde)
   If(IGL.Eq.1) ErrMax=XNorm1
EndDo

deallocate(A1,WORK0,C0Tilde,C1Tilde,C2Tilde,Lambda,APLUS0Tilde,APLUS1Tilde)
deallocate(ABMIN1,ABPLUS1,WORK1)
allocate(WorkD(NDimX,NCholesky))
WorkD=0
WorkD = RESHAPE(COMTilde, (/NDimX, NCholesky/))
ECorr=0
do j=1,NDimX
   do i=1,NCholesky
      ECorr=ECorr+DCholAct(i,j)*WorkD(j,i)
   enddo
enddo

deallocate(WorkD,COMTilde)

Call RELEASE_AC0BLOCK(A0Block,A0blockIV,nblk)

!inquire(file='rdm2.dat',exist=irdm2)
!if(irdm2) then
!   open(newunit=iunit,file='rdm2.dat',status='old')
!   close(iunit,status='delete')
!endif

end subroutine WIter_D12Chol

subroutine WIter_DChol(ECorr,Max_Cn,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)
!
!  AC energy cacluation using CHOLESKY VECTORS:
!  (1) expanding AC integrand in alpha around alpha=0, up to Max_Cn order
!  (2) finding C^(n)[omega] and (3) omega integration
!  Max_Cn is passed in input.inp (if missed, the dafault value 3 is used)
!
use abfofo
use systemdef

implicit none
integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision :: ECorr,ECorrAct,EGOne(NGem)
double precision :: XFreq(100),WFreq(100)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile

double precision, allocatable :: DChol(:,:),DCholT(:,:),DCholAct(:,:),DCholActT(:,:),WorkD(:,:)
double precision, allocatable :: APLUS0Tilde(:), APLUS1Tilde(:),APLUS0TildeAct(:), APLUS1TildeAct(:),  &
                                 A1(:),A2(:), &
                                 COMTilde(:),COMTildeAct(:),ABPLUS0(:),ABMIN0(:),ABPLUS1(:),ABMIN1(:), &
                                 LAMBDA(:), &
                                 C0Tilde(:),C1Tilde(:),C2Tilde(:), &
                                 WORK0(:)
integer :: NCholesky

integer :: nblk
type(EblockData) :: A0blockIV
type(EblockData),allocatable :: A0block(:)

interface
subroutine read_D_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
   double precision, allocatable, intent(out) :: DChol(:,:), DCholAct(:,:)
   integer :: NCholesky
   integer, intent(in) :: NDimX, NBasis, IndN(2,NDimX), IndAux(NBasis)
   double precision, intent(in) :: Occ(NBasis)
end subroutine read_D_array
end interface

! Get DChol & DCholAct
call read_D_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
DCholT = transpose(DChol)
DCholActT = transpose(DCholAct)
! ==========================================================================

NGrid=18
!NGrid=30

Write (6,'(/,X,''AC calculation through W_AC expansion around Alpha=0, Omega Grid = '',I3,&
   '' and max order in C expansion = '',I3,/)') NGrid,Max_Cn

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

allocate(ABPLUS0(NDimX*NDimX),ABMIN0(NDimX*NDimX),ABPLUS1(NDimX*NDimX),ABMIN1(NDimX*NDimX))

ACAlpha=0.D0
call AB_CAS_FOFO(ABPLUS0,ABMIN0,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha,.false.)

ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS1,ABMIN1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha,.false.)

ABPLUS1=ABPLUS1-ABPLUS0
ABMIN1 =ABMIN1 -ABMIN0
EGOne(1)=ECASSCF

!Calc: A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
allocate(A1(NDimX*NDimX))
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,ABMIN1,NDimX,0.0d0,A1,NDimX)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,ABMIN0,NDimX,1d0,A1,NDimX)
deallocate(ABMIN0)

!Calc: A2=ABPLUS1*ABMIN1
allocate(A2(NDimX*NDimX))
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,ABMIN1,NDimX,0.0d0,A2,NDimX)
deallocate(ABMIN1)

!Calc: APLUS0Tilde=ABPLUS0.DChol
allocate(APLUS0Tilde(NDimX*NCholesky))
Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS0,NDimX,DCholT,NDimX,0.0d0,APLUS0Tilde,NDimX)
allocate(APLUS0TildeAct(NDimX*NCholesky))
Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS0,NDimX,DCholActT,NDimX,0.0d0,APLUS0TildeAct,NDimX)
deallocate(ABPLUS0)

!Calc: APLUS1Tilde=ABPLUS1.DChol
allocate(APLUS1Tilde(NDimX*NCholesky))
Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholT,NDimX,0.0d0,APLUS1Tilde,NDimX)
allocate(APLUS1TildeAct(NDimX*NCholesky))
Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholActT,NDimX,0.0d0,APLUS1TildeAct,NDimX)
deallocate(ABPLUS1)

Call FreqGrid(XFreq,WFreq,NGrid)
! Calc: A0
nblk = 1 + NBasis - NAct
allocate(A0block(nblk))
Call AC0BLOCK(Occ,URe,XOne, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1,'FFOO','FOFO', &
     A0BlockIV,A0Block,nblk,'A0BLK',0)

allocate(COMTilde(NDimX*NCholesky),COMTildeAct(NDimX*NCholesky))
COMTilde=0.0
COMTildeAct=0.0

allocate(C0Tilde(NDimX*NCholesky),C1Tilde(NDimX*NCholesky),C2Tilde(NDimX*NCholesky),WORK0(NDimX*NCholesky))
allocate(LAMBDA(NDimX*NDimX))
Do IGL=1,NGrid
   OmI=XFreq(IGL)
   WFact=4.D0/PI*WFreq(IGL)

!  Calc: LAMBDA=(A0+Om^2)^-1
   Call INV_AC0BLK_OLD(OmI**2,LAMBDA,A0Block,A0BlockIV,nblk,NDimX)

   Do K=1,2

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
   If(K==1) Then
     Call dgemm('N','N',NDimX,NCholesky,NDimX,0.5d0,LAMBDA,NDimX,APLUS0Tilde,NDimX,0.0d0,C0Tilde,NDimX)
   Else
     Call dgemm('N','N',NDimX,NCholesky,NDimX,0.5d0,LAMBDA,NDimX,APLUS0TildeAct,NDimX,0.0d0,C0Tilde,NDimX)
   EndIf

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,A1,NDimX,C0Tilde,NDimX,0.0d0,WORK0,NDimX)
   If(K==1) Then
      WORK0=0.5d0*APLUS1Tilde-WORK0
   Else
      WORK0=0.5d0*APLUS1TildeAct-WORK0
   EndIf
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,LAMBDA,NDimX,WORK0,NDimX,0.0d0,C1Tilde,NDimX)

   If(K==1) Then
      COMTilde=COMTilde+WFact*0.5d0*C1Tilde
   Else
      COMTildeAct=COMTildeAct+WFact*0.5d0*C1Tilde
   EndIf

   XFactorial=1
   Do N=2,Max_Cn
       XFactorial=XFactorial*N
       XN1=-N
       XN2=-N*(N-1)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,XN2,A2,NDimX,C0Tilde,NDimX,0.0d0,WORK0,NDimX)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,XN1,A1,NDimX,C1Tilde,NDimX,1.0d0,WORK0,NDimX)
       Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,LAMBDA,NDimX,WORK0,NDimX,0.0d0,C2Tilde,NDimX)
       FF=WFact/XFactorial/(N+1)
       If(K==1) Then
          COMTilde=COMTilde+FF*C2Tilde
       Else
          COMTildeAct=COMTildeAct+FF*C2Tilde
       EndIf
       C0Tilde=C1Tilde
       C1Tilde=C2Tilde
   EndDo

   EndDo

   Write(6,'(X,"Omega, |C|, |C_Act|",I3,3F10.4)')IGL,OmI,Norm2(COMTilde),Norm2(COMTildeAct)

EndDo

allocate(WorkD(NDimX,NCholesky))
WorkD=0
WorkD = RESHAPE(COMTilde, (/NDimX, NCholesky/))
ECorr=0
do j=1,NDimX
   do i=1,NCholesky
      ECorr=ECorr+DChol(i,j)*WorkD(j,i)
   enddo
enddo
! active part of ecorr
WorkD=0
WorkD = RESHAPE(COMTildeAct, (/NDimX, NCholesky/))
ECorrAct=0
do j=1,NDimX
   do i=1,NCholesky
      ECorrAct=ECorrAct+DCholAct(i,j)*WorkD(j,i)
   enddo
enddo
ECorr=ECorr-ECorrAct

deallocate(WorkD,COMTilde,COMTildeAct)

close(iunit)
Call RELEASE_AC0BLOCK(A0Block,A0blockIV,nblk)

end subroutine WIter_DChol

subroutine WIter_FOFO(ECorr,Max_Cn,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)
!
!  AC energy cacluation by: (1) expanding AC integrand in alpha around alpha=0, up to Max_Cn order
!  (2) finding C^(n)[omega] and (3) omega integration
!
use abmat
use abfofo
!use types,only : EblockData
use blocktypes
use print_units

implicit none

integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision,intent(inout) :: ECorr,EGOne(NGem)

double precision :: COM(NDimX*NDimX),XFreq(100),WFreq(100),&
                 ABPLUS0(NDim*NDim),WORK0(NDim*NDim),ABPLUS1(NDim*NDim),WORK1(NDim*NDim),&
                 A1(NDimX*NDimX),A2(NDimX*NDimX),&
                 C0(NDimX*NDimX),C1(NDimX*NDimX),C2(NDimX*NDimX)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs,iipq,iirs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)

integer :: nblk
type(EblockData) :: A0blockIV
type(EblockData),allocatable :: A0block(:)

NGrid=18
Write (6,'(/,X,''AC calculation through W_AC expansion around Alpha=0, Omega Grid = '',I3,&
   '' and max order in C expansion = '',I3,/)') NGrid,Max_Cn

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

ACAlpha=0.D0
call AB_CAS_FOFO(ABPLUS0,WORK0,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha,.false.)

ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha,.false.)

ABPLUS1=ABPLUS1-ABPLUS0
WORK1=WORK1-WORK0
EGOne(1)=ECASSCF

!Calc: A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK1,NDimX,0d0,A1,NDimX)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK0,NDimX,1d0,A1,NDimX)
!Calc: A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0d0,A2,NDimX)

Call FreqGrid(XFreq,WFreq,NGrid)

! Calc: A0
nblk = 1 + NBasis - NAct
allocate(A0block(nblk))
Call AC0BLOCK(Occ,URe,XOne, &
      IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1,'FFOO','FOFO', &
      A0BlockIV,A0Block,nblk,'A0BLK',0)
      !A0BlockIV,A0Block,nblk,1)

COM=0d0
Do IGL=1,NGrid
   OmI=XFreq(IGL)

!  Calc: WORK1=(A0+Om^2)^-1
   Call INV_AC0BLK(OmI**2,WORK1,A0Block,A0BlockIV,nblk,NDimX)
!  Calc: C0=1/2 Lambda.ABPLUS0
   Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
              ABPLUS0,NDimX,0d0,C0,NDimX)
!  Calc: WORK0=Lambda.A1
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK1,NDimX,&
              A1,NDimX,0d0,WORK0,NDimX)
!  Calc: C1=1/2 Lambda.ABPLUS1-WORK0.C0
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK0,NDimX,&
              C0,NDimX,0d0,C1,NDimX)
   Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
              ABPLUS1,NDimX,-1.d0,C1,NDimX)
!  Calc: WORK1=LAMBDA*A2
   Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,WORK1,NDimX,&
              A2,NDimX,0d0,C2,NDimX)
   WORK1=C2

!  FROM NOW ON: Lambda.A2 in WORK1, Lambda.A1 in WORK0
   WFact=4.D0/PI*WFreq(IGL)
   COM=COM+WFact*0.5D0*C1

   XFactorial=1
   Do N=2,Max_Cn
       XFactorial=XFactorial*N
!  Calc: Cn
       XN1=-N
       XN2=-N*(N-1)
       Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK0,NDimX,&
                  C1,NDimX,0d0,C2,NDimX)
       Call dgemm('N','N',NDimX,NDimX,NDimX,XN2,WORK1,NDimX,&
                  C0,NDimX,XN1,C2,NDimX)
       FF=WFact/XFactorial/(N+1)

       COM=COM+FF*C2
       C0=C1
       C1=C2
   EndDo
EndDo

do i=1,NBasis
CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
do k=1,3
   do j=1,3
      do i=1,3
         if((i==j).and.(j==k).and.(k==l)) then
            AuxCoeff(i,j,k,l) = .false.
         endif
      enddo
   enddo
enddo
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis))

ECorr = 0

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
  access='DIRECT',recl=8*NBasis*NOccup)

kl   = 0
SumY = 0
do k=1,NOccup
do l=1,NBasis
   kl = kl + 1
   if(pos(l,k)/=0) then
     irs = pos(l,k)
     ir = l
     is = k
     read(iunit,rec=kl) work(1:NBasis*NOccup)
     do j=1,NOccup
        do i=1,NBasis
           ints(i,j) = work((j-1)*NBasis+i)
        enddo
     enddo
     ints(:,NOccup+1:NBasis) = 0

     do j=1,NBasis
        do i=1,j
           if(pos(j,i)/=0) then
             ipq = pos(j,i)
             ip = j
             iq = i
             Crs = CICoef(l)+CICoef(k)
             Cpq = CICoef(j)+CICoef(i)
             if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then
                Aux = Crs*Cpq*COM((irs-1)*NDimX+ipq)
                ECorr = ECorr + Aux*ints(j,i)
             endif
           endif
        enddo
     enddo

   endif

enddo
enddo

close(iunit)

Call RELEASE_AC0BLOCK(A0Block,A0blockIV,nblk)

deallocate(ints,work)

end subroutine WIter_FOFO

subroutine WInteg_FOFO(ECorr,XOne,URe,Occ,&
   EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,&
   NBasis,NInte1,NDim,NGem,IndAux,ACAlpha,&
   IndN,IndX,NDimX)
!
!  for a given ACAlpha, compute W_AC returned as ECorr by C(omega) integration
!  C(omega) is found by direct inversion of the [A+A- + omega^2] matrix
!
use abfofo

implicit none

integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision,intent(inout) :: ECorr,EGOne(NGem)

double precision :: COM(NDimX*NDimX),ipiv(NDimX),XFreq(100),WFreq(100),&
                 ABPLUS(NDim*NDim),ABMIN(NDim*NDim),AIN(NDimX*NDimX),CMAT(NDimX*NDimX)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf
double precision :: ECASSCF,PI,WFact,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,0,ACAlpha,.false.)
EGOne(1)=ECASSCF

!     Frequency integration of CMAT

NGrid=15
Write (6,'(/,X,''AC Calculation with Omega Grid = '',I3,/)') NGrid
Call FreqGrid(XFreq,WFreq,NGrid)

COM=0d0
Do IGL=1,NGrid
   OmI=XFreq(IGL)
   AIN=0d0
   Do I=1,NDimX
       AIN((I-1)*NDimX+I)=1.0
   EndDo
!     ABPLUS*ABMIN+1 Om^2
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,&
              ABMIN,NDimX,OmI**2,AIN,NDimX)

   CMAT=ABPLUS
   Call dgesv(NDimX,NDimX,AIN,NDimX,ipiv,CMAT,NDimX,inf)

   COM=COM+2.D0/PI*CMAT*WFreq(IGL)
EndDo

do i=1,NBasis
CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
do k=1,3
   do j=1,3
      do i=1,3
         if((i==j).and.(j==k).and.(k==l)) then
            AuxCoeff(i,j,k,l) = .false.
         endif
      enddo
   enddo
enddo
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis))

ECorr = 0

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
  access='DIRECT',recl=8*NBasis*NOccup)

kl   = 0
SumY = 0
do k=1,NOccup
do l=1,NBasis
   kl = kl + 1
   if(pos(l,k)/=0) then
     irs = pos(l,k)
     ir = l
     is = k
     read(iunit,rec=kl) work(1:NBasis*NOccup)
     do j=1,NOccup
        do i=1,NBasis
           ints(i,j) = work((j-1)*NBasis+i)
        enddo
     enddo
     ints(:,NOccup+1:NBasis) = 0

     do j=1,NBasis
        do i=1,j
           if(pos(j,i)/=0) then
             ipq = pos(j,i)
             ip = j
             iq = i
             Crs = CICoef(l)+CICoef(k)
             Cpq = CICoef(j)+CICoef(i)

             if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then

                Aux = Crs*Cpq*COM((irs-1)*NDimX+ipq)
         ! this term is divergent and it has been removed by computing Int C^alpha=0 domega
         !       if(iq.eq.is.and.ip.Eq.ir) then
         !          Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
         !        endif

                ECorr = ECorr + Aux*ints(j,i)
             endif

           endif
        enddo
     enddo

   endif

enddo
enddo

close(iunit)

deallocate(ints,work)

end subroutine WInteg_FOFO

subroutine CIter_FOFO_old(PMat,ECorr,ACAlpha,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)

use abfofo

implicit none

integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha,ACAlpha0,XMix
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision :: ECorr,ECorrAct,EGOne(NGem)

!double precision,intent(in) :: PMat(NDimX,NDimX)
double precision :: PMat(NDimX,NDimX)

double precision :: COM(NDimX*NDimX),ipiv(NDimX),XFreq(100),WFreq(100),&
                 ABPLUS0(NDim*NDim),WORK0(NDim*NDim),ABPLUS1(NDim*NDim),WORK1(NDim*NDim),&
                 A0(NDimX*NDimX),A1(NDimX*NDimX),A2(NDimX*NDimX),&
                 C0(NDimX*NDimX),CMAT(NDimX*NDimX)&
,AIN(NDimX*NDimX),ABPLUS(NDim*NDim),ABMIN(NDim*NDim)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:),DChol(:,:),WorkD(:,:),DCholAct(:,:)
integer :: NCholesky

NGrid=15
Max_Cn=15
XMix=0.6

! safe choice
!Max_Cn=25
!XMix=0.6

Write (6,'(/,X,''AC Iterative Calculation with Omega Grid = '',I3,&
   '' and max order in C expansion = '',I3,/)') NGrid,Max_Cn

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

ACAlpha0=0.D0
call AB_CAS_FOFO(ABPLUS0,WORK0,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha0,.false.)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0d0,A0,NDimX)

call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,1,ACAlpha,.false.)
EGOne(1)=ECASSCF
!A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0d0,A2,NDimX)

A0=0.D0
Do I=1,NDimX
 A0((I-1)*NDimX+I)=A2((I-1)*NDimX+I)
EndDo

A2=A2-A0

! PROJECT A2 WITH PMat : WORK0=PMat.A2
!Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,PMat,NDimX,A2,NDimX,0.0,WORK0,NDimX)
!A0=A0+WORK0
!print*, 'A2 norm beofore projection',norm2(A2(1:NDimX*NDimX))
!A2=A2-WORK0
!print*, 'A2 norm after projection',norm2(A2(1:NDimX*NDimX))


Call FreqGrid(XFreq,WFreq,NGrid)
COM=0d0
Do IGL=1,NGrid
   OmI=XFreq(IGL)

   WORK0=0.D0
   Do I=1,NDimX
       WORK0((I-1)*NDimX+I)=OmI**2
   EndDo

   WORK1=A0+WORK0

! Initial=WORK1=(A0+WORK0)^-1
   Call dgetrf(NDimX, NDimX, WORK1, NDimX, ipiv, inf1 )
   Call dgetri(NDimX, WORK1, NDimX, ipiv, work0, NDimX, inf2 )

!     C0=C(0)
   Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,WORK1,NDimX,&
              ABPLUS1,NDimX,0d0,C0,NDimX)
!     C1=C(1)
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,A2,NDimX,&
              C0,NDimX,0d0,WORK0,NDimX)
   WORK0=ABPLUS1-WORK0
   Call dgemm('N','N',NDimX,NDimX,NDimX,1.0d0,WORK1,NDimX,&
              WORK0,NDimX,0d0,CMAT,NDimX)
   Do N=2,Max_Cn
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,A2,NDimX,&
              CMAT,NDimX,0d0,WORK0,NDimX)
      WORK0=ABPLUS1-WORK0
! damping is needed when active orbitals present: CMAT(n) = (1-XMix)*CMAT(n) + XMix*CMAT(n-1)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1.0d0-XMix,WORK1,NDimX,&
              WORK0,NDimX,XMix,CMAT,NDimX)
   EndDo
   COM=COM+2.D0/PI*CMAT*WFreq(IGL)
EndDo
!
! end of computing the integral of C(Omega)
!
do i=1,NBasis
CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
do k=1,3
   do j=1,3
      do i=1,3
         if((i==j).and.(j==k).and.(k==l)) then
            AuxCoeff(i,j,k,l) = .false.
         endif
      enddo
   enddo
enddo
enddo


open(newunit=iunit,file='cholvecs',form='unformatted')
read(iunit) NCholesky
allocate(WorkD(NCholesky,NBasis**2))
read(iunit) WorkD
close(iunit)

print*,'NCholesky',NCholesky

allocate(DChol(NCholesky,NDimX))
allocate(DCholAct(NCholesky,NDimX))

DChol = 0
DCholAct = 0
do j=1,NDimX
ir=IndN(1,j)
is=IndN(2,j)
irs = is+(ir-1)*NBasis
Crs=CICoef(ir)+CICoef(is)
do i=1,NCholesky
    DChol(i,j) = Crs*WorkD(i,irs)
enddo
if(IndAux(ir)*IndAux(is)==1) then
    do i=1,NCholesky
       DCholAct(i,j) = Crs*WorkD(i,irs)
    enddo
endif
enddo
deallocate(WorkD)

!!!!!!!!!!!!!!!
allocate(WorkD(NDimX,NCholesky))
WorkD=0
! WorkD=C_tilde=C*D^T
Call dgemm('N','T',NDimX,NCholesky,NDimX,1d0,COM,NDimX,&
              DChol,NCholesky,0d0,WorkD,NDimX)
ECorr=0
do j=1,NDimX
do i=1,NCholesky
   ECorr=ECorr+DChol(i,j)*WorkD(j,i)
enddo
enddo
! active part of ecorr
WorkD=0
! WorkD=C_tilde_act=C*D_act^T
Call dgemm('N','T',NDimX,NCholesky,NDimX,1d0,COM,NDimX,&
              DCholAct,NCholesky,0d0,WorkD,NDimX)
ECorrAct=0
do j=1,NDimX
do i=1,NCholesky
   ECorrAct=ECorrAct+DCholAct(i,j)*WorkD(j,i)
enddo
enddo
deallocate(WorkD,DChol,DCholAct)
ECorr=ECorr-ECorrAct

return


!allocate(work(NBasis**2),ints(NBasis,NBasis))
!
!ECorr = 0
!
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBasis*NOccup)
!
!kl   = 0
!do k=1,NOccup
!   do l=1,NBasis
!      kl = kl + 1
!      if(pos(l,k)/=0) then
!        irs = pos(l,k)
!        ir = l
!        is = k
!        read(iunit,rec=kl) work(1:NBasis*NOccup)
!        do j=1,NOccup
!           do i=1,NBasis
!              ints(i,j) = work((j-1)*NBasis+i)
!           enddo
!        enddo
!        ints(:,NOccup+1:NBasis) = 0
!
!        do j=1,NBasis
!           do i=1,j
!              if(pos(j,i)/=0) then
!                ipq = pos(j,i)
!                ip = j
!                iq = i
!                Crs = CICoef(l)+CICoef(k)
!                Cpq = CICoef(j)+CICoef(i)
!
!                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then
!                   Aux = Crs*Cpq*COM((irs-1)*NDimX+ipq)
!         ! this term is divergent and it has been removed by computing Int C^alpha=0 domega
!         !          if(iq.eq.is.and.ip.Eq.ir) then
!         !             Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
!         !          endif
!                  ECorr = ECorr + Aux*ints(j,i)
!                endif
!
!              endif
!           enddo
!        enddo
!
!      endif
!
!   enddo
!enddo
!
!close(iunit)
!
!deallocate(ints,work)

end subroutine CIter_FOFO_old

!subroutine Eccor_iter()
subroutine CIter_FOFO(PMat,ECorr,ACAlpha,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,IndN,IndX,NDimX)
!
!  to do: get rid of multiplication of NDimX.NDimX matrices
!  reduce the number of matrices: do not keep A0, only diagonals of Lambda
!  use an efficient algorithm (lower dimensional) to invert the matrix with PMat
!  turn on projection if the reduction of the norm of A2 after projection is low enough
!
   use abfofo
   use class_IterStats
   use class_CIntegrator
   use class_IterAlgorithm
   use class_IterAlgorithmDIIS
   use class_IterAlgorithmDamping
   use class_LambdaCalculator
   use class_LambdaCalculatorDiag
   use class_LambdaCalculatorBlock
   use class_LambdaCalculatorProjector

   implicit none

   integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
   integer,intent(in) :: NAct,INActive,NELE
   integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),IGem(NBasis)
   double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
   double precision :: ACAlpha
   double precision :: ECorr,ECorrAct,EGOne(NGem)
   double precision,intent(in) :: PMat(NDimX,NDimX)
   double precision :: ABPLUS1(NDim*NDim),A0(NDimX*NDimX),A2(NDimX*NDimX),WORK1(NDim*NDim)
   integer :: NOccup,NGrid
   integer :: i,j
   double precision :: ECASSCF,PI
   character(:),allocatable :: twojfile,twokfile,IntKFile
   double precision, allocatable :: DChol(:,:),DCholT(:,:),DCholAct(:,:),DCholActT(:,:),WorkD(:,:)
   double precision, allocatable :: APlusTilde(:), APlusTildeAct(:), COMTilde(:),COMTildeAct(:)
   integer :: NCholesky
   type(CIntegrator) :: CIntegr
   class(IterAlgorithm), allocatable, target :: iterAlgo
   class(LambdaCalculator), allocatable, target :: lambdaCalc

   integer :: nblk

   interface
   subroutine read_D_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
      double precision, allocatable, intent(out) :: DChol(:,:), DCholAct(:,:)
      integer :: NCholesky
      integer, intent(in) :: NDimX, NBasis, IndN(2,NDimX), IndAux(NBasis)
      double precision, intent(in) :: Occ(NBasis)
   end subroutine read_D_array
   end interface

   ! Get DChol & DCholAct
   ! ==========================================================================
   call read_D_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)
   DCholT = transpose(DChol)
   DCholActT = transpose(DCholAct)
   ! ==========================================================================

   NOccup = NAct + INActive
   PI = 4.0*ATAN(1.0)

   twojfile = 'FFOO'
   twokfile = 'FOFO'
   IntKFile = twokfile

   call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
                  IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                  NInte1,twojfile,twokfile,1,ACAlpha,.false.)
   EGOne(1)=ECASSCF
   ! Calc A2=ABPLUS1*ABMIN1
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0d0,A2,NDimX)

   ! Calc A+Tilde & AAct+Tilde
   ! ==========================================================================================================================
   allocate(APlusTilde(NDimX*NCholesky), APlusTildeAct(NDimX*NCholesky))
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholT,NDimX,0d0,APlusTilde,NDimX)
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholActT,NDimX,0d0,APlusTildeAct,NDimX)
   ! ==========================================================================================================================

   ! Calc CTilde & CTildeAct integrals
  ! ==========================================================================================================================
   allocate(COMTilde(NDimX*NCholesky),COMTildeAct(NDimX*NCholesky))
   COMTilde=0d0
   COMTildeAct=0d0

   ! Create iteration algorithm object
   !iterAlgo = IterAlgorithmDIIS(Threshold=1d-4, DIISN=6, maxIterations=30)
   !!iterAlgo = IterAlgorithmDamping(Threshold=1d-4, XMix=0.2, maxIterations=-1)
   allocate(iterAlgo,SOURCE=IterAlgorithmDIIS(Threshold=1d-4, DIISN=6, maxIterations=30))

   ! Create A0 calculator object
!!    lambdaCalc = LambdaCalculatorDiag(NDimX, A2)
   !LambdaCalc = LambdaCalculatorBlock(NDimX, A2, URe, Occ, XOne, IndN, IndX, IGem, NBasis, NAct, INActive, NInte1)
   allocate(LambdaCalc,SOURCE=LambdaCalculatorBlock(NDimX, A2, URe, Occ, XOne, IndN, IndX, IGem, NBasis, NAct, INActive, NInte1))
!!    LambdaCalc = LambdaCalculatorProjector(NDimX, A2, PMat)

   NGrid = 35
   CIntegr = CIntegrator(iterAlgo=iterAlgo, lambdaCalc=lambdaCalc)
   call CIntegr%setup(NGrid, NDimX, NCholesky, APlusTilde, APlusTildeAct, ACAlpha)
   ! call CIntegr%integrate(COMTilde, COMTildeAct)
   call CIntegr%integrateReverse(COMTilde, COMTildeAct)

   ! Free memory
   call lambdaCalc%clean()
   call CIntegr%clean()
   deallocate(iterAlgo,LambdaCalc)
   ! ==========================================================================================================================

   allocate(WorkD(NDimX,NCholesky))
   WorkD=0
   WorkD = RESHAPE(COMTilde, (/NDimX, NCholesky/))

   ECorr=0
   do j=1,NDimX
      do i=1,NCholesky
         ECorr=ECorr+DChol(i,j)*WorkD(j,i)
      enddo
   enddo
   ! active part of ecorr
   WorkD=0
   WorkD = RESHAPE(COMTildeAct, (/NDimX, NCholesky/))

   ECorrAct=0
   do j=1,NDimX
      do i=1,NCholesky
         ECorrAct=ECorrAct+DCholAct(i,j)*WorkD(j,i)
      enddo
   enddo
   deallocate(WorkD,DChol,DCholAct, APlusTilde, APlusTildeAct, COMTilde, COMTildeAct)
   ECorr=ECorr-ECorrAct

   return

end subroutine CIter_FOFO


subroutine read_D12_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)

   implicit none
   integer, intent(in) :: NDimX, NBasis, IndN(2,NDimX), IndAux(NBasis)
   double precision, intent(in) :: Occ(NBasis)
   double precision, allocatable, intent(out) :: DChol(:,:), DCholAct(:,:)
   integer :: NCholesky, iunit, i, j, ir, is, irs
   double precision, allocatable :: WorkD(:,:)
   double precision :: Crs, CICoef(NBasis)

   open(newunit=iunit,file='cholvecs',form='unformatted')
   read(iunit) NCholesky
   allocate(WorkD(NCholesky,NBasis**2))
   read(iunit) WorkD
   close(iunit)

   print*,'NCholesky',NCholesky

   allocate(DChol(NCholesky,NDimX), DCholAct(NCholesky,NDimX))

   do i=1,NBasis
      CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
   enddo

   DChol = 0
   DCholAct = 0
   do j=1,NDimX
      ir=IndN(1,j)
      is=IndN(2,j)
      irs = is+(ir-1)*NBasis
      Crs=CICoef(ir)+CICoef(is)
      do i=1,NCholesky
!           DChol(i,j) = Crs*WorkD(i,irs)
            DCholAct(i,j) = Crs*WorkD(i,irs)
!           if(IndAux(ir)*IndAux(is)==1) DChol(i,j) = 2.D0*DChol(i,j)
            if(IndAux(ir)*IndAux(is)==1) DCholAct(i,j) = 2.D0*DCholAct(i,j)
      enddo
      if(IndAux(ir)*IndAux(is).ne.1) then
            do i=1,NCholesky
!               DCholAct(i,j) = Crs*WorkD(i,irs)
                DChol(i,j) = Crs*WorkD(i,irs)
            enddo
      endif
   enddo
   deallocate(WorkD)

end subroutine read_D12_array

subroutine read_D_array(NCholesky, DChol, DCholAct, NDimX, NBasis, IndN, Occ, IndAux)

   implicit none
   integer, intent(in) :: NDimX, NBasis, IndN(2,NDimX), IndAux(NBasis)
   double precision, intent(in) :: Occ(NBasis)
   double precision, allocatable, intent(out) :: DChol(:,:), DCholAct(:,:)
   integer :: NCholesky, iunit, i, j, ir, is, irs
   double precision, allocatable :: WorkD(:,:)
   double precision :: Crs, CICoef(NBasis)

   open(newunit=iunit,file='cholvecs',form='unformatted')
   read(iunit) NCholesky
   allocate(WorkD(NCholesky,NBasis**2))
   read(iunit) WorkD
   close(iunit)

   print*,'NCholesky',NCholesky

   allocate(DChol(NCholesky,NDimX), DCholAct(NCholesky,NDimX))

   do i=1,NBasis
      CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
   enddo

   DChol = 0
   DCholAct = 0
   do j=1,NDimX
      ir=IndN(1,j)
      is=IndN(2,j)
      irs = is+(ir-1)*NBasis
      Crs=CICoef(ir)+CICoef(is)
      do i=1,NCholesky
            DChol(i,j) = Crs*WorkD(i,irs)
      enddo
      if(IndAux(ir)*IndAux(is)==1) then
            do i=1,NCholesky
               DCholAct(i,j) = Crs*WorkD(i,irs)
            enddo
      endif
   enddo
   deallocate(WorkD)

end subroutine read_D_array

subroutine AC0BLOCK(Occ,URe,XOne, &
                    IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
                    IntJFile,IntKFile,ICholesky, &
                    A0BlockIV,A0block,nblk,ver,dumpfile,dump)
!
!     A ROUTINE FOR COMPUTING : a) ver=0  ABPLUS^{(0)} and ABMIN^{(0)}
!                                         (stored in matY and matX, respectively)
!                               b) ver=1  A0=ABPLUS^{(0)}.ABMIN^{(0)}
!                 (FOFO VERSION)
!
use abfofo
!use types,only : EblockData
use blocktypes
!
implicit none

integer,intent(in)           :: NAct,INActive,NDimX,NBasis,NDim,NInte1
integer,intent(in)           :: ICholesky
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
integer                      :: nblk
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
character(*)                 :: IntJFile,IntKFile
character(*)                 :: dumpfile
integer,intent(in)           :: ver,dump

integer          :: iunit
integer          :: NOccup
integer          :: i,j,k,l,kl,ii,ip,iq,ir,is,ipq,irs
integer          :: ipos,jpos,iblk,jblk
integer          :: IGem(NBasis),Ind(NBasis),pos(NBasis,NBasis)
integer          :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer          :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                    tmpAV(NAct,INActive+NAct+1:NBasis),&
                    tmpIV(INActive*(NBasis-NAct-INActive))
integer          :: limAA(2),limAI(2,1:INActive),&
                    limAV(2,INActive+NAct+1:NBasis),limIV(2)
double precision :: AuxCoeff(3,3,3,3),Aux,val,ETot
double precision :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision :: C(NBasis)
! hererXXX
!integer :: ICol,IRow
!double precision :: xn1,xn2

type(EblockData) :: A0block(nblk), A0blockIV

! set dimensions
NOccup = NAct + INActive

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

call ABPM0_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
                IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
                IntJFile,IntKFile,ICholesky,ETot)
! hererXXX
!do i=1,NBasis
!   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
!enddo
!xn1=0.25
!call create_pos_ABPL0(pos,IGem,IndN,INActive,NAct,NBasis,NDimX)
!do ICol=1,NDimX
!   ir = IndN(1,ICol)
!   is = IndN(2,ICol)
!   irs = pos(ir,is)
!   do IRow=1,NDimX
!         ip = IndN(1,IRow)
!         iq = IndN(2,IRow)
!         ipq = pos(ip,iq)
!         if(ipq.Eq.irs) then
!            xn2=(Occ(IP)-Occ(IQ))*(Occ(IP)+1.D0-Occ(IQ))*xn1
!            ABPLUS(ipq,irs)=ABPLUS(ipq,irs)+xn2/(C(IP)+C(IQ))**2
!            ABMIN(ipq,irs)=ABMIN(ipq,irs)+xn2/(C(IP)-C(IQ))**2
!         endif
!   enddo
!enddo

call create_blocks_ABPL0(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,pos,&
                         IGem,IndN,INActive,NAct,NBasis,NDimX)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

!allocate(A0block(1+NBasis-NAct))
nblk = 0

!pack AA
if(nAA>0) then
   nblk = nblk + 1
   call pack_A0block(ABPLUS,ABMIN,nAA,limAA(1),limAA(2),tmpAA,A0block(nblk),NDimX,ver)
endif
!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAI(iq),limAI(1,iq),limAI(2,iq),tmpAI(1:nAI(iq),iq),&
                        A0block(nblk),NDimX,ver)
   endif
enddo
!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAV(ip),limAV(1,ip),limAV(2,ip),tmpAV(1:nAV(ip),ip),&
                        A0block(nblk),NDimX,ver)
    endif
enddo
!pack IV
associate(B => A0blockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))

  if(ver==0) then
     do i=1,B%n
        ii = B%l1+i-1
        B%vec(i) = ABPLUS(ii,ii)
     enddo
  elseif(ver==1) then
     do i=1,B%n
        ii = B%l1+i-1
        B%vec(i) = ABPLUS(ii,ii)*ABMIN(ii,ii)
     enddo
  endif

end associate

if(dump==1) then
  ! dump to file
  open(newunit=iunit,file=dumpfile,form='unformatted')
  write(iunit) nblk
  do iblk=1,nblk
     associate(B => A0Block(iblk))
       write(iunit) iblk, B%n, B%l1, B%l2
       write(iunit) B%pos,B%matX
     end associate
  enddo
  associate(B => A0BlockIV)
    write(iunit) B%n,B%l1,B%l2
    write(iunit) B%pos,B%vec
  end associate
  close(iunit)

  !!deallocate blocks
  !do iblk=1,nblk
  !   associate(B => A0block(iblk))
  !     deallocate(B%matX,B%pos)
  !   end associate
  !enddo
  !! (IV part)
  !associate(B => A0blockIV)
  !  deallocate(B%pos,B%vec)
  !end associate
endif

end subroutine AC0BLOCK

subroutine INV_AC0BLK(omega,Lambda,LambdaIV,A0Block,A0BlockIV,nblk,NDimX)
!
! Calculate A0Inv=(A0+Om^2)^-1
! (note: this is the same as calculateLambda_blk in SAPT)
!
use blocktypes

integer,intent(in)           :: nblk,NDimX
double precision,intent(in)  :: omega
type(EblockData)             :: A0block(nblk), A0blockIV

type(EblockData)             :: Lambda(nblk), LambdaIV

integer             :: i,j,ipos,jpos,info
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:)

Lambda = A0block

do iblk=1,nblk
   associate(B => Lambda(iblk))

     if(B%n == 1) then
        ! hartree-fock case
        B%matX(1,1) = 1d0 / (B%matX(1,1) + omega)
     else
        ! more than 1 active orbital
        allocate(ipiv(B%n),work(B%n))

        do i=1,B%n
            B%matX(i,i) = B%matX(i,i) + omega
        enddo

        call dgetrf(B%n,B%n,B%matX,B%n,ipiv,info)
        call dgetri(B%n,B%matX,B%n,ipiv,work,B%n,info)

        deallocate(work,ipiv)
     endif

   end associate
enddo

! IV block (allocations done before)
associate(B => A0blockIV, A => LambdaIV)
  A%pos = B%pos
  do i=1,B%n
     A%vec(i) = 1d0 / (B%vec(i) + omega)
  enddo
end associate

end subroutine INV_AC0BLK

subroutine INV_AC0BLK_OLD(omega,A0Inv,A0Block,A0BlockIV,nblk,NDimX)
!
! Calculate A0Inv=(A0+Om^2)^-1 and return full A0inv matrix!
!
use blocktypes

type(EblockData)   :: A0block(nblk), A0blockIV
type(EblockData)   :: W0block(nblk)
integer,intent(in) :: nblk,NDimX
double precision,intent(in)  :: omega
double precision,intent(out) :: A0Inv(NDimX,NDimX)

integer             :: i,j,ipos,jpos,info
integer,allocatable :: ipiv(:)
double precision,allocatable :: work(:)

W0block = A0block

A0Inv = 0
do iblk=1,nblk
   associate(B => W0block(iblk))

     allocate(ipiv(B%n),work(B%n))

     do i=1,B%n
         B%matX(i,i) = B%matX(i,i) + omega
     enddo

     call dgetrf(B%n,B%n,B%matX,B%n,ipiv,info)
     call dgetri(B%n,B%matX,B%n,ipiv,work,B%n,info)

     ! assign places in work(NDimX,NDimX) matrix
     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           A0Inv(ipos,jpos) = B%matX(i,j)
        enddo
     enddo

     deallocate(work,ipiv)

   end associate
enddo

! IV block
associate(B => A0blockIV)
do i=1,B%n
   ii = B%pos(i)
   A0Inv(ii,ii) = 1d0 / (B%vec(i) + omega)
enddo
end associate

end subroutine INV_AC0BLK_OLD

subroutine RELEASE_AC0BLOCK(A0block,A0blockIV,nblk)

!use types,only : EblockData
use blocktypes

implicit none

integer,intent(in) :: nblk
type(EblockData)   :: A0block(nblk), A0blockIV
integer            :: iblk

! deallocate blocks
do iblk=1,nblk
   associate(B => A0block(iblk))
     deallocate(B%matX,B%pos)
   end associate
enddo
! (IV part)
associate(B => A0blockIV)
  deallocate(B%pos,B%vec)
end associate

end subroutine RELEASE_AC0BLOCK

subroutine READ_AC0BLK(omega,A0Inv,a0blkfile,NDimX)
! an example of reading-in a0blocks from file
! and calculating (A0+omega)^-1
!use types,only : EblockData
use blocktypes

implicit none

double precision    :: omega
integer,intent(in)  :: NDimX
double precision :: A0Inv(NDimX,NDimX)
character(*),intent(in) :: a0blkfile

integer :: i,iblk,iunit
integer :: nblk
type(EBlockData)             :: A0BlockIV
type(EBlockData),allocatable :: A0Block(:)
! for testing:
integer             :: ii,j,ipos,jpos,info
double precision    :: diff
integer,allocatable :: ipiv(:)
double precision,allocatable :: work1(:),work(:,:)

! read from file
open(newunit=iunit,file=a0blkfile,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
allocate(A0Block(nblk))
do iblk=1,nblk
   associate(A => A0Block(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n))
     read(iunit) A%pos,A%matX
   end associate
enddo
associate(A => A0BlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  read(iunit) A%pos,A%vec
end associate
close(iunit)

! invert small blocks
allocate(work(NDimX,NDimX))

work = 0
do iblk=1,nblk
   associate(B => A0block(iblk))

     allocate(ipiv(B%n),work1(B%n))

     do i=1,B%n
         B%matX(i,i) = B%matX(i,i) + omega
     enddo

     call dgetrf(B%n,B%n,B%matX,B%n,ipiv,info)
     call dgetri(B%n,B%matX,B%n,ipiv,work1,B%n,info)
     ! assign places in work(NDimX,NDimX) matrix
     do j=1,B%n
        jpos=B%pos(j)
        do i=1,B%n
           ipos=B%pos(i)
           work(ipos,jpos) = B%matX(i,j)
        enddo
     enddo

     deallocate(work1,ipiv)

   end associate
enddo

! IV block
associate(B => A0blockIV)
do i=1,B%n
   ii = B%pos(i)
   work(ii,ii) = 1d0 / (B%vec(i) + omega)
enddo
end associate

A0Inv=work
!print*, 'test2-2:',norm2(A0Inv)

deallocate(work)

end subroutine READ_AC0BLK

subroutine pack_A0block(ABPLUS,ABMIN,nVal,lim1,lim2,tmpMat,Eblock,NDimX,ver)
!
! packs either :  ver = 0 , ABP, ABM (stored in matY and matX, respectively)
!                 ver = 1 ABP.ABM    (stored in matX)
!
use blocktypes

implicit none

integer,intent(in) :: NDimX
integer,intent(in) :: nVal,lim1,lim2,tmpMat(nVal)
integer,intent(in) :: ver
double precision :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
type(EblockData) :: Eblock(1)
double precision,allocatable :: ABP(:,:),ABM(:,:)

associate(B => Eblock(1))
  B%n  = nVal
  B%l1 = lim1
  B%l2 = lim2
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpMat(1:B%n)

  !allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
  allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
  ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
  if(ver == 0) then
     allocate(B%matX(B%n,B%n),B%matY(B%n,B%n))
     B%matX = ABM
     B%matY = ABP
  elseif(ver == 1) then
     allocate(B%matX(B%n,B%n))
     call dgemm('N','N',B%n,B%n,B%n,1d0,ABP,B%n,ABM,B%n,0d0,B%matX,B%n)
  else
     write(*,*) 'Wrong call for pack_A0block',ver
     stop
  endif
  deallocate(ABM,ABP)

end associate

end subroutine pack_A0block

subroutine Polariz(FreqOm,UNOAO,XOne,URe,Occ,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NGem,IndAux,&
   IndN,IndX,NDimX,ICholesky)
!
! Returns dynamic polarizability tensor for a given frequency
! find C(omega) by inversion
!
use abfofo

implicit none
integer,intent(in) :: NBasis,NInte1,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),IndAux(NBasis),IGem(NBasis)
double precision,intent(in) :: FreqOm,UNOAO(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis),CICoef(NBasis)
double precision :: DipCX(NDimX),DipCY(NDimX),DipCZ(NDimX)
double precision :: ipiv(NDimX),ABPLUS(NDimX*NDimX),ABMIN(NDimX*NDimX),AIN(NDimX*NDimX),CMAT(NDimX*NDimX)
double precision :: ECASSCF,AXX,AYX,AXY,AZX,AXZ,AYY,AZY,AYZ,AZZ,Om,ddot,Alpha
character(:),allocatable :: twojfile,twokfile
integer :: I,J,IJ,inf,ICholesky,NOccup

Om=FreqOm

NOccup=NAct+INActive
Call ComputeDipoleMom(UNOAO,Occ,NOccup,'AOONEINT.mol','DIP',NBasis)

Call ReadDip(DipX,DipY,DipZ,UNOAO,'DIP',NBasis)

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

twojfile = 'FFOO'
twokfile = 'FOFO'

Alpha=1.0
Call AB_CAS_FOFO(ABPLUS,ABMIN,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,ICholesky,Alpha,.false.)
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

Write(6,'(/,X,''Polarizability tensor for frequency '',F8.4)') Om
Write(6,'(/,X,''XX   XY   XZ  '',3F15.8)') AXX, AXY, AXZ
Write(6,'(X,''YX   YY   YZ  '',3F15.8)') AYX, AYY, AYZ
Write(6,'(X,''ZX   ZY   ZZ  '',3F15.8,2/)') AZX, AZY, AZZ

end subroutine Polariz

subroutine PolarizAl(FreqOm,UNOAO,XOne,URe,Occ,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NGem,IndAux,&
   IndN,IndX,NDimX,ICholesky,Max_Cn)
!
! Returns dynamic polarizability tensor for a given frequency FreqOm
! find C(omega) by expanding around Alpha=0 with a tolerance Eps or
! up to maximal order Max_Cn
!
use abfofo

implicit none
integer,intent(in) :: NBasis,NInte1,NGem,NDimX,Max_cn
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),IndAux(NBasis),IGem(NBasis)
double precision,intent(in) :: FreqOm,UNOAO(NBasis,NBasis),URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
double precision :: DipX(NBasis,NBasis),DipY(NBasis,NBasis),DipZ(NBasis,NBasis),CICoef(NBasis)
double precision :: DipCX(NDimX),DipCY(NDimX),DipCZ(NDimX)
double precision :: ipiv(NDimX)
double precision :: ECASSCF,AXX,AYX,AXY,AZX,AXZ,AYY,AZY,AYZ,AZZ,Om,ddot,Alpha
character(:),allocatable :: twojfile,twokfile
integer :: I,J,IJ,inf,ICholesky,NOccup

Om=FreqOm

NOccup=NAct+INActive
Call ComputeDipoleMom(UNOAO,Occ,NOccup,'AOONEINT.mol','DIP',NBasis)

Call ReadDip(DipX,DipY,DipZ,UNOAO,'DIP',NBasis)

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

Call CFREQPROJ(ipiv,Om,DipCX,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX,&
   'FFOO','FOFO')

AYX=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZX=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AXX=8.d0*ddot(NDimx,DipCX,1,ipiv,1)

Call CFREQPROJ(ipiv,Om,DipCY,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX,&
   'FFOO','FOFO')

AXY=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AZY=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)
AYY=8.d0*ddot(NDimx,DipCY,1,ipiv,1)

Call CFREQPROJ(ipiv,Om,DipCZ,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX,&
   'FFOO','FOFO')
AXZ=8.d0*ddot(NDimx,DipCX,1,ipiv,1)
AYZ=8.d0*ddot(NDimx,DipCY,1,ipiv,1)
AZZ=8.d0*ddot(NDimx,DipCZ,1,ipiv,1)

Write(6,'(/,X,''Polarizability tensor for frequency '',F8.4)') Om
Write(6,'(/,X,''XX   XY   XZ  '',3F15.8)') AXX, AXY, AXZ
Write(6,'(X,''YX   YY   YZ  '',3F15.8)') AYX, AYY, AYZ
Write(6,'(X,''ZX   ZY   ZZ  '',3F15.8,2/)') AZX, AZY, AZZ

end subroutine PolarizAl

subroutine CFREQPROJ(COMTilde,OmI,DProj,NProj, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NBasis,NInte1,IndAux,&
   ICholesky,IndN,IndX,NDimX,&
   IntJFile,IntKFile)
!
!  For a given frequency OmI, return a product of the matrices C(Alpha=1,OmI) and DProj
!  where DProj is of the NProj x NDimX size
!  C is found by expanding in alpha around alpha=0
!  with a tolerance Eps or truncating at Max_Cn order
!  
use abfofo
use systemdef
use sapt_utils

implicit none
integer,intent(in) :: NBasis,NInte1,NDimX,NProj
integer,intent(in) :: NAct,INActive,ICholesky
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),IndAux(NBasis),&
                      IGem(NBasis)
double precision :: ACAlpha,Eps
double precision,intent(in) :: DProj(NProj,NDimX),URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
character(*),intent(in) :: IntJfile,IntKFile

double precision,intent(out) :: COMTilde(NDimX*NProj)

integer :: iunit
integer :: ia,ib,ic,id
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: inf1,inf2,Max_Cn
double precision :: XFactorial,XN1,XN2,OmI,XNorm0,XNorm1,ECASSCF,DProjT(NDimX,NProj)

double precision, allocatable :: APLUS0Tilde(:), APLUS1Tilde(:),  &
                                 A1(:),&
                                 ABPLUS0(:),ABMIN0(:),ABPLUS1(:),ABMIN1(:), &
                                 C0Tilde(:),C1Tilde(:),C2Tilde(:), &
                                 WORK0(:),WORK1(:)
integer :: nblk,N
type(EblockData) :: A0blockIV,LambdaIV
type(EblockData),allocatable :: A0block(:),Lambda(:)

! tolerance
Eps=1.d-5

DProjT = transpose(DProj)

allocate(ABPLUS1(NDimX*NDimX),ABMIN1(NDimX*NDimX))

if(NAct==1) then
  ! active-virtual block
  nblk = NBasis - NAct - INActive
else
  nblk = 1 + NBasis - NAct
endif

allocate(A0block(nblk))
! AC0BLOCK with ver=0 stores A-(0) and A+(0) matrices
!                            in X and Y, respectively
Call AC0BLOCK(Occ,URe,XOne, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1, &
     IntJFile,IntKFile, &
     ICholesky,A0BlockIV,A0Block,nblk,0,'DUMMY',0)

! get AB1PLUS and AB1MIN
ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS1,ABMIN1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntJFile,IntKFile,ICholesky,ACAlpha,.false.)

Call sq_symmetrize(ABPLUS1,NDimX)
Call sq_symmetrize(ABMIN1,NDimX)

! AB1 = AB1 - A0
call add_blk_right(ABPLUS1,A0Block,A0BlockIV,-1d0,.false.,nblk,NDimX)
call add_blk_right(ABMIN1, A0Block,A0BlockIV,-1d0,.true., nblk,NDimX)

!Calc: A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
allocate(A1(NDimX*NDimX))
call ABPM_HALFTRAN_GEN_L(ABMIN1, A1,0.0d0,A0Block,A0BlockIV,nblk,NDimX,NDimX,'Y')
call ABPM_HALFTRAN_GEN_R(ABPLUS1,A1,1.0d0,A0Block,A0BlockIV,nblk,NDimX,NDimX,'X')

!Calc: APLUS0Tilde=ABPLUS0.DProj
allocate(APLUS0Tilde(NDimX*NProj))
call ABPM_HALFTRAN_GEN_L(DProjT,APLUS0Tilde,0.0d0,A0Block,A0BlockIV,nblk,NDimX,NProj,'Y')

!Calc: APLUS1Tilde=ABPLUS1.DProj
allocate(APLUS1Tilde(NDimX*NProj))
Call dgemm('N','N',NDimX,NProj,NDimX,1d0,ABPLUS1,NDimX,DProjT,NDimX,0.0d0,APLUS1Tilde,NDimX)

deallocate(A0block)
deallocate(A0BlockIV%vec,A0BlockIV%pos)

! Calc: A0
allocate(A0block(nblk))
! ver=1: store A+(0).A-(0) in blocks
Call AC0BLOCK(Occ,URe,XOne, &
     IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1, &
     IntJFile,IntKFile, &
     ICholesky,A0BlockIV,A0Block,nblk,1,'A0BLK',0)

COMTilde=0.0

allocate(C0Tilde(NDimX*NProj),C1Tilde(NDimX*NProj),C2Tilde(NDimX*NProj),WORK0(NDimX*NProj))
allocate(WORK1(NDimX*NProj))
allocate(Lambda(nblk))
associate(A => A0BlockIV, L => LambdaIV)
  L%n = A%n
  L%l1 = A%l1
  L%l2 = A%l2
  allocate(L%pos(L%n),L%vec(L%n))
end associate

!  Calc: LAMBDA=(A0-Om^2)^-1
Call INV_AC0BLK(-OmI**2,Lambda,LambdaIV,A0Block,A0BlockIV,nblk,NDimX)

!  Calc: C0Tilde=1/2 LAMBDA.APLUS0Tilde
Call ABPM_HALFTRAN_GEN_L(APLUS0Tilde,C0Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NProj,'X')
C0Tilde = 0.5d0*C0Tilde

!  Calc: C1Tilde=LAMBDA.(1/2 APLUS1Tilde - A1.C0Tilde)
Call dgemm('N','N',NDimX,NProj,NDimX,1.d0,A1,NDimX,C0Tilde,NDimX,0.0d0,WORK0,NDimX)
WORK0 = 0.5d0*APLUS1Tilde - WORK0
Call ABPM_HALFTRAN_GEN_L(WORK0,C1Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NProj,'X')

COMTilde=C0Tilde
COMTilde=COMTilde+C1Tilde

XNorm0=1.0d5
XFactorial=1
Do N=2,Max_Cn
    XFactorial=XFactorial*N
    XN1=-N
    XN2=-N*(N-1)
    Call dgemm('N','N',NDimX,NProj,NDimX,XN2,ABMIN1,NDimX,C0Tilde,NDimX,0.0d0,WORK1,NDimX)
    Call dgemm('N','N',NDimX,NProj,NDimX,1.d0,ABPLUS1,NDimX,WORK1,NDimX,0.0d0,WORK0,NDimX)
    Call dgemm('N','N',NDimX,NProj,NDimX,XN1,A1,NDimX,C1Tilde,NDimX,1.0d0,WORK0,NDimX)
    Call ABPM_HALFTRAN_GEN_L(WORK0,C2Tilde,0.0d0,Lambda,LambdaIV,nblk,NDimX,NProj,'X')
    XNorm1=Norm2(C2Tilde/XFactorial)
    Write(6,'(X,"Order (n), |C^(n)/n!|",I3,E14.4)')N,XNorm1
! herer
!    If(XNorm1.Le.Eps) Exit
    If(N.Gt.3.And.XNorm1.Gt.XNorm0) Then
       Write(6,'(X,"Divergence detected. Expansion of C terminated at order ",I3,3F10.4)')N-1
       Exit
    EndIf
    XNorm0=XNorm1
    COMTilde=COMTilde+C2Tilde/XFactorial
    C0Tilde=C1Tilde
    C1Tilde=C2Tilde
EndDo

deallocate(A1,WORK0,C0Tilde,C1Tilde,C2Tilde,Lambda,APLUS0Tilde,APLUS1Tilde)
deallocate(ABMIN1,ABPLUS1,WORK1)

Call RELEASE_AC0BLOCK(A0Block,A0blockIV,nblk)

end subroutine CFREQPROJ 
