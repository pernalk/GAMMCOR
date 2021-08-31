subroutine WChol_FOFO(PMat,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)
!
! compute the projector matrix PMat, used optionally in ACFREQ
!
use abfofo

implicit none

integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision,intent(inout) :: EGOne(NGem)

double precision :: PMat(NDimX,NDimX)

double precision :: COM(NDimX*NDimX),XFreq(100),WFreq(100),&
                 ABPLUS(NDim*NDim),AB(NDim*NDim),ADIAG(NDim)

integer :: iunit
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),N,IGL,inf1,inf2,info,Max_Cn
double precision :: ECASSCF,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)
Real*8, Allocatable :: MatFF(:,:),work2(:,:),work3(:),work1(:,:),work4(:),work5(:,:),work6(:,:)
integer,allocatable :: ipiv(:)
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
   work3(i)=0.0
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
end subroutine WChol_FOFO

subroutine WIter_FOFO(ECorr,XOne,URe,Occ,EGOne,NGOcc,&
   IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
   IndN,IndX,NDimX)
!
!  AC energy cacluation by: (1) expanding AC integrand in alpha around alpha=0, up to Max_Cn order
!  (2) finding C^(n)[omega] and (3) omega integration 
!
use abmat
use abfofo
use types,only : LOUT,EblockData

implicit none

integer,intent(in) :: NGOcc,NBasis,NInte1,NDim,NGem,NDimX
integer,intent(in) :: NAct,INActive,NELE
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),&
                   IGem(NBasis)
double precision :: ACAlpha
double precision,intent(in) :: URe(NBasis,NBasis),Occ(NBasis),XONe(NInte1)
double precision,intent(inout) :: ECorr,EGOne(NGem)

double precision :: COM(NDimX*NDimX),ipiv(NDimX),XFreq(100),WFreq(100),&
                 ABPLUS0(NDim*NDim),WORK0(NDim*NDim),ABPLUS1(NDim*NDim),WORK1(NDim*NDim),&
                 A1(NDimX*NDimX),A2(NDimX*NDimX),&
                 C0(NDimX*NDimX),C1(NDimX*NDimX),C2(NDimX*NDimX)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)

integer :: nblk
type(EblockData) :: A0blockIV
type(EblockData),allocatable :: A0block(:)

NGrid=20
Max_Cn=5
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
              NInte1,twojfile,twokfile,ACAlpha,.false.)
ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,ACAlpha,.false.)
ABPLUS1=ABPLUS1-ABPLUS0
WORK1=WORK1-WORK0
EGOne(1)=ECASSCF

!Calc: A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK1,NDimX,0.0,A1,NDimX)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK0,NDimX,1d0,A1,NDimX)
!Calc: A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0.0,A2,NDimX)

Call FreqGrid(XFreq,WFreq,NGrid)

! Calc: A0 and dump it
nblk = 1 + NBasis - NAct
allocate(A0block(nblk))
Call AC0BLOCK(Occ,URe,XOne, &
      IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,NInte1,'FFOO','FOFO', &
      A0BlockIV,A0Block,nblk,0)
      !A0BlockIV,A0Block,nblk,1)
!
COM=0.0
Do IGL=1,NGrid
   OmI=XFreq(IGL)

!  Calc: WORK1=(A0+Om^2)^-1
   Call INV_AC0BLK(OmI**2,WORK1,A0Block,A0BlockIV,nblk,NDimX)
   !Call READ_AC0BLK(OmI**2,WORK1,'A0BLK',NDimX)

!  Calc: C0=1/2 Lambda.ABPLUS0
   Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
              ABPLUS0,NDimX,0.0,C0,NDimX)
!  Calc: WORK0=Lambda.A1
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK1,NDimX,&
              A1,NDimX,0.0,WORK0,NDimX)
!  Calc: C1=1/2 Lambda.ABPLUS1-WORK0.C0
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK0,NDimX,&
              C0,NDimX,0.0,C1,NDimX)
   Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
              ABPLUS1,NDimX,-1.d0,C1,NDimX)
!  Calc: WORK1=LAMBDA*A2
   Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,WORK1,NDimX,&
              A2,NDimX,0.0,C2,NDimX)
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
                  C1,NDimX,0.0,C2,NDimX)
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
              NInte1,twojfile,twokfile,ACAlpha,.false.)
EGOne(1)=ECASSCF

!     Frequency integration of CMAT 

NGrid=15
Write (6,'(/,X,''AC Calculation with Omega Grid = '',I3,/)') NGrid
Call FreqGrid(XFreq,WFreq,NGrid)

COM=0.0
Do IGL=1,NGrid
   OmI=XFreq(IGL)
   AIN=0.0
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
              NInte1,twojfile,twokfile,ACAlpha0,.false.)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0.0,A0,NDimX)

call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,twojfile,twokfile,ACAlpha,.false.)
EGOne(1)=ECASSCF
!A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0.0,A2,NDimX)

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
COM=0.0
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
              ABPLUS1,NDimX,0.0,C0,NDimX)
!     C1=C(1)
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,A2,NDimX,&
              C0,NDimX,0.0,WORK0,NDimX)
   WORK0=ABPLUS1-WORK0
   Call dgemm('N','N',NDimX,NDimX,NDimX,1.0d0,WORK1,NDimX,&
              WORK0,NDimX,0.0,CMAT,NDimX)
   Do N=2,Max_Cn
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,A2,NDimX,&
              CMAT,NDimX,0.0,WORK0,NDimX)
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
              DChol,NCholesky,0.0,WorkD,NDimX)
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
              DCholAct,NCholesky,0.0,WorkD,NDimX)
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
                  NInte1,twojfile,twokfile,ACAlpha,.false.)
   EGOne(1)=ECASSCF
   ! Calc A2=ABPLUS1*ABMIN1
   Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0.0,A2,NDimX)
   
   ! Calc A+Tilde & AAct+Tilde
   ! ==========================================================================================================================   
   allocate(APlusTilde(NDimX*NCholesky), APlusTildeAct(NDimX*NCholesky))
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholT,NDimX,0.0,APlusTilde,NDimX)
   Call dgemm('N','N',NDimX,NCholesky,NDimX,1d0,ABPLUS1,NDimX,DCholActT,NDimX,0.0,APlusTildeAct,NDimX)
   ! ==========================================================================================================================  

   ! Calc CTilde & CTildeAct integrals
  ! ==========================================================================================================================  
   allocate(COMTilde(NDimX*NCholesky),COMTildeAct(NDimX*NCholesky))
   COMTilde=0.0
   COMTildeAct=0.0

   ! Create iteration algorithm object
   iterAlgo = IterAlgorithmDIIS(Threshold=1d-3, DIISN=16, maxIterations=30)
!    iterAlgo = IterAlgorithmDamping(Threshold=1d-3, XMix=0.2, maxIterations=-1)

   ! Create A0 calculator object
   lambdaCalc = LambdaCalculatorDiag()
!    LambdaCalc = LambdaCalculatorBlock(URe,Occ,XOne,IndN,IndX,IGem,NBasis,NAct,INActive,NInte1,twojfile,twokfile)
!    LambdaCalc = LambdaCalculatorProjector(PMat)

   NGrid = 35
   CIntegr = CIntegrator(iterAlgo=iterAlgo, lambdaCalc=lambdaCalc)
   call CIntegr%setup(NGrid, NDimX, NCholesky, A0, A2, APlusTilde, APlusTildeAct, ACAlpha)
!    call CIntegr%integrate(COMTilde, COMTildeAct)
   call CIntegr%integrateReverse(COMTilde, COMTildeAct)
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
                    IntJFile,IntKFile,A0BlockIV,A0block,nblk,dump)
!
!     A ROUTINE FOR COMPUTING A0=ABPLUS^{(0)}.ABMIN^{(0)}
!                 (FOFO VERSION)
!
use abfofo
use types,only : LOUT,EblockData
!
implicit none

integer,intent(in)           :: NAct,INActive,NDimX,NBasis,NDim,NInte1
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
integer                      :: nblk
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
character(*)                 :: IntJFile,IntKFile
integer,intent(in)           :: dump

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
                IntJFile,IntKFile,ETot)

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
   call pack_A0block(ABPLUS,ABMIN,nAA,limAA(1),limAA(2),tmpAA,A0block(nblk),NDimX)
endif
!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAI(iq),limAI(1,iq),limAI(2,iq),tmpAI(1:nAI(iq),iq),&
                        A0block(nblk),NDimX)
   endif
enddo
!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAV(ip),limAV(1,ip),limAV(2,ip),tmpAV(1:nAV(ip),ip),&
                        A0block(nblk),NDimX)
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
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)*ABMIN(ii,ii)
  enddo

end associate

if(dump==1) then
  ! dump to file
  open(newunit=iunit,file='A0BLK',form='unformatted')
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

subroutine INV_AC0BLK(omega,A0Inv,A0Block,A0BlockIV,nblk,NDimX)
! Calculate A0Inv=(A0+Om^2)^-1
use types,only : EblockData

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
!print*, 'test2-2:',norm2(A0Inv)

end subroutine INV_AC0BLK

subroutine RELEASE_AC0BLOCK(A0block,A0blockIV,nblk)
use types,only : EblockData
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
use types,only : LOUT,EblockData

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

! example: diagonalize small blocks
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

subroutine pack_A0block(ABPLUS,ABMIN,nVal,lim1,lim2,tmpMat,Eblock,NDimX)

use types,only : LOUT,EblockData

implicit none

integer,intent(in) :: NDimX
integer,intent(in) :: nVal,lim1,lim2,tmpMat(nVal)
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
  allocate(B%matX(B%n,B%n))
  allocate(ABP(B%n,B%n),ABM(B%n,B%n))

  ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
  ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)

  call dgemm('N','N',B%n,B%n,B%n,1d0,ABP,B%n,ABM,B%n,0d0,B%matX,B%n)

  deallocate(ABM,ABP)

end associate

end subroutine pack_A0block


