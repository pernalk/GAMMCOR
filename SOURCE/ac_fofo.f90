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
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0d0,A0,NDimX)
!A1=ABPLUS0*ABMIN1+ABPLUS1*ABMIN0
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK1,NDimX,0d0,A1,NDimX)
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK0,NDimX,1d0,A1,NDimX)
!A2=ABPLUS1*ABMIN1
Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS1,NDimX,WORK1,NDimX,0d0,A2,NDimX)

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
                 ABPLUS0,NDimX,0d0,C0,NDimX)
!     WORK0=LAMBDA*A1
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK1,NDimX,&
                 A1,NDimX,0d0,WORK0,NDimX)
!     C1=C(1)
      Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,WORK0,NDimX,&
                 C0,NDimX,0d0,C1,NDimX)
      Call dgemm('N','N',NDimX,NDimX,NDimX,0.5d0,WORK1,NDimX,&
                 ABPLUS1,NDimX,-1.d0,C1,NDimX)
!     WORK1=LAMBDA*A2
      Call dgemm('N','N',NDimX,NDimX,NDimX,1.d0,WORK1,NDimX,&
                 A2,NDimX,0d0,C2,NDimX)
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

NGrid=18
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
