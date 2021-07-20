subroutine WChol_FOFO(PMat,XOne,URe,Occ,EGOne,NGOcc,&
      IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
      IndN,IndX,NDimX)

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

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf1,inf2,info,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)
Real*8, Allocatable :: MatFF(:,:),work2(:,:),work3(:),work1(:,:),work4(:),work5(:,:),work6(:,:)
integer,allocatable :: ipiv(:)
integer :: NCholesky,lwork

NGrid=25

NOccup = NAct + INActive
PI = 4.0*ATAN(1.0)

twojfile = 'FFOO'
twokfile = 'FOFO'
IntKFile = twokfile

open(newunit=iunit,file='cholvecs',form='unformatted')
read(iunit) NCholesky
allocate(MatFF(NCholesky,NBasis**2))
read(iunit) MatFF
close(iunit)

!!!
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
      print*, i, work3(i)
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
!print*,work6(1,1),work6(20,20),work6(1,20)
print*,'norm of D.D^T.[D.D^T]^-1',norm2(work6(1:NCholesky,1:NCholesky)),'Sqrt(NCholesky)',SQRT(Float(NCholesky))

work2=work5
deallocate(work4,work3,work5,work6,work1)

!! (D.D^T)^-1
!allocate(ipiv(NCholesky))
! allocate(work3(NCholesky*NCholesky))
!call dgetrf(NCholesky,NCholesky,work2,NCholesky,ipiv,info)
!!print*, 'info-1',info
!if(info>0) then
!   write(lout,*) 'WChol_FOFO: dgetrf ',info
!   stop
!endif
!call dgetri(NCholesky,work2,NCholesky,ipiv,work3,NCholesky,info)
!!print*, 'info-2',info,work3(1)
!if(info>0) then
!   write(lout,*) 'WChol_FOFO: dgetri ',info
!   stop
!endif
!deallocate(work3,ipiv)

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

allocate(work1(NDimX,NDimX))
work1=work2
! truncate PMat
!work1 = 0
!do j=1,NDimX
!   ir=IndN(1,j)
!   is=IndN(2,j)
!   irs = is+(ir-1)*NBasis
!   do i=1,NDimX
!      ip=IndN(i,1)
!      iq=IndN(i,2)
!      ipq = iq+(ip-1)*NBasis
!      work1(i,j) = work2(ipq,irs)
!   enddo
!enddo
!deallocate(work2)


allocate(work3(NDimX*NDimX))
ACAlpha=1.D0
call AB_CAS_FOFO(ABPLUS,work3,ECASSCF,URe,Occ,XOne, &
                 IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                 NInte1,twojfile,twokfile,ACAlpha,.false.)

Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,work3,NDimX,0.0,AB,NDimX)

ADIAG = 0
do ipq=1,NDimX
   ADIAG(ipq) = AB((ipq-1)*NDimX+ipq)
   AB((ipq-1)*NDimX+ipq) = 0d0
enddo
print*, 'NDIAG and DIAG norms',norm2(AB(1:NDimX*NDimX)),norm2(ADIAG(1:NDimX))

! project
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,AB,NDimX,0.0,work3,NDimX)


do i=1,NDimX
   do j=1,NDimX
!      if(abs(work3(i+(j-1)*ndimx)-ab(i+(j-1)*ndimx)).gt.1.d-6) &
!      if(abs(ab(i+(j-1)*ndimx)).gt.1.d-6) &
!      write(*,*)i,j,work3(i+(j-1)*ndimx)-ab(i+(j-1)*ndimx),ab(i+(j-1)*ndimx)
   enddo
enddo


work3=AB-work3
print*, 'before and after projection',norm2(AB(1:NDimX*NDimX)),norm2(work3(1:NDimX*NDimX))
stop

deallocate(work1)
end subroutine WChol_FOFO

subroutine WIter_FOFO(ECorr,XOne,URe,Occ,EGOne,NGOcc,&
      IGem,NAct,INActive,NELE,NBasis,NInte1,NDim,NGem,IndAux,&
      IndN,IndX,NDimX)

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
                    ABPLUS0(NDim*NDim),WORK0(NDim*NDim),ABPLUS1(NDim*NDim),WORK1(NDim*NDim),&
                    A0(NDimX*NDimX),A1(NDimX*NDimX),A2(NDimX*NDimX),&
                    C0(NDimX*NDimX),C1(NDimX*NDimX),C2(NDimX*NDimX)

integer :: iunit,NOccup
integer :: ia,ib,ic,id,ICol,IRow
integer :: i,j,k,l,kl,ip,iq,ir,is,ipq,irs
integer :: pos(NBasis,NBasis),NGrid,N,IGL,inf1,inf2,Max_Cn
double precision :: ECASSCF,PI,WFact,XFactorial,XN1,XN2,FF,CICoef(NBasis),Cpq,Crs,SumY,Aux,OmI
character(:),allocatable :: twojfile,twokfile,IntKFile
logical :: AuxCoeff(3,3,3,3)
double precision,allocatable :: work(:),ints(:,:)

NGrid=25
Max_Cn=3
Write (6,'(/,X,''AC Iterative Calculation with Omega Grid = '',I3,&
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

!A0=ABPLUS0*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0.0,A0,NDimX)
!A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK1,NDimX,0.0,A1,NDimX)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK0,NDimX,1d0,A1,NDimX)
!A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0.0,A2,NDimX)

Call FreqGrid(XFreq,WFreq,NGrid)

COM=0.0
Do IGL=1,NGrid
      OmI=XFreq(IGL)

      WORK0=0.D0
      Do I=1,NDimX
          WORK0((I-1)*NDimX+I)=OmI**2
      EndDo

      WORK1=A0+WORK0

      Call dgetrf(NDimX, NDimX, WORK1, NDimX, ipiv, inf1 )
      Call dgetri(NDimX, WORK1, NDimX, ipiv, work0, NDimX, inf2 )

!     C0=C(0)
      Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
                 ABPLUS0,NDimX,0.0,C0,NDimX)
!     WORK0=LAMBDA*A1
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK1,NDimX,&
                 A1,NDimX,0.0,WORK0,NDimX)
!     C1=C(1)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK0,NDimX,&
                 C0,NDimX,0.0,C1,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
                 ABPLUS1,NDimX,-1.d0,C1,NDimX)
!     WORK1=LAMBDA*A2
      Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,WORK1,NDimX,&
                 A2,NDimX,0.0,C2,NDimX)
      WORK1=C2

!     FROM NOW ON: LAMBDA*A2 in WORK1, LAMBDA*A1 in WORK0
      WFact=4.D0/PI*WFreq(IGL)

      COM=COM+WFact*(C0+0.5D0*C1)

      XFactorial=1
      Do N=2,Max_Cn
          XFactorial=XFactorial*N
!         C(n)
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
                   if(iq.eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif

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

end subroutine WIter_FOFO

subroutine WInteg_FOFO(ECorr,XOne,URe,Occ,&
      EGOne,NGOcc,&
      IGem,NAct,INActive,NELE,&
      NBasis,NInte1,NDim,NGem,IndAux,ACAlpha,&
      IndN,IndX,NDimX)

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
                   if(iq.eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif

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

subroutine CIter_FOFO(PMat,ECorr,ACAlpha,XOne,URe,Occ,EGOne,NGOcc,&
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
double precision,intent(inout) :: ECorr,EGOne(NGem)

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
double precision,allocatable :: work(:),ints(:,:)

NGrid=15
Max_Cn=10
XMix=0.4

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
call AB_CAS_FOFO(ABPLUS1,WORK1,ECASSCF,URe,Occ,XOne, &
                 IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
                 NInte1,twojfile,twokfile,ACAlpha,.false.)
EGOne(1)=ECASSCF

!A0=ABPLUS0*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0.0,A0,NDimX)
!A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0.0,A2,NDimX)

A0=0.D0
Do I=1,NDimX
    A0((I-1)*NDimX+I)=A2((I-1)*NDimX+I)
EndDo
! herer
Do I=1,0
!Do I=1,NDimX
   IR=IndN(1,I)
   IS=IndN(2,I)
   Do J=1,NDimX
      IP=IndN(1,J)
      IQ=IndN(2,J)
!      if(abs(PMat(I,J)).gt.1.d-1) write(*,*)ir,is,ip,iq,PMat(I,J) 
         If(IndAux(IR).Eq.1.And.IndAux(IS).Eq.1.And.IndAux(IP).Eq.1.And.IndAux(IQ).Eq.1) Then
!!           write(*,*)'act-act',ir,is,ip,iq
           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
         EndIf
         If(IndAux(IR).Eq.1.And.IndAux(IS).Eq.0.And.IndAux(IP).Eq.1.And.IS.Eq.IQ) Then
!           write(*,*)'act-inact',ir,is,ip,iq
           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
         EndIf
!         If(IndAux(IR).Eq.2.And.IndAux(IS).Eq.1.And.IndAux(IQ).Eq.1.And.IR.Eq.IP) Then
!           write(*,*)'virt-nact',ir,is,ip,iq
!           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
 !        EndIf
         If(IndAux(IR).Eq.2.And.IndAux(IS).Eq.1.And.IndAux(IP).Eq.2.And.IS.Eq.IQ) Then
!!           write(*,*)'virt-act',ir,is,ip,iq
           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
         EndIf
         IA=IndAux(IR)*IndAux(IS)*IndAux(IP)*IndAux(IQ)
         IB=0
         SumY=Occ(IP)+Occ(IQ)+Occ(IR)+Occ(IS)
         If(IndAux(IR).Eq.1.And.IndAux(IS).Eq.1.And.IndAux(IP).Eq.1.And.IndAux(IQ).Eq.1) IB=1
         If((IndAux(IR).Eq.1.Or.IndAux(IS).Eq.1.Or.IndAux(IP).Eq.1.Or.IndAux(IQ).Eq.1).And.IA.Eq.2.And.IB.Eq.0) Then
!           write(*,*)ir,is,ip,iq
!           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
         EndIf
         If((IndAux(IR).Eq.1.Or.IndAux(IS).Eq.1.Or.IndAux(IP).Eq.1.Or.IndAux(IQ).Eq.1).And.IA.Eq.4.And.IB.Eq.0) Then
     !      write(*,*)ir,is,ip,iq
!           A0((J-1)*NDimX+I)=A2((J-1)*NDimX+I)
         EndIf

   EndDo
EndDo

A2=A2-A0

!
! WORK0=P.A2
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,PMat,NDimX,A2,NDimX,0.0,WORK0,NDimX)

A0=A0+WORK0
print*, 'A2 norm beofore projection',norm2(A2(1:NDimX*NDimX))
A2=A2-WORK0
print*, 'A2 norm after projection',norm2(A2(1:NDimX*NDimX))

Call FreqGrid(XFreq,WFreq,NGrid)

COM=0.0
Do IGL=1,NGrid
      OmI=XFreq(IGL)

      WORK0=0.D0
      Do I=1,NDimX
          WORK0((I-1)*NDimX+I)=OmI**2
      EndDo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      WORK0=A2+A0+WORK0
!      CMAT=ABPLUS1
!      Call dgesv(NDimX,NDimX,WORK0,NDimX,ipiv,CMAT,NDimX,inf1)
!      COM=COM+2.D0/PI*CMAT*WFreq(IGL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      WORK1=A0+WORK0
      
!      WORK1=A2+A0+WORK0
! LAMBDA=WORK1=(A0+WORK0)^-1
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
!         Call dgemm('N','N',NDimX,NDimX,NDimX,1.0d0,WORK1,NDimX,&
!                 WORK0,NDimX,0.0,CMAT,NDimX)
         Call dgemm('N','N',NDimX,NDimX,NDimX,1.0d0-XMix,WORK1,NDimX,&
                 WORK0,NDimX,XMix,CMAT,NDimX)
      EndDo
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
                   if(iq.eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif

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

end subroutine CIter_FOFO
