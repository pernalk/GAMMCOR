#define RDMCORR_DEBUG -1

module rdmcorr

use print_units
use types
use abfofo
use sapt_utils

contains

subroutine RDMResp_FOFO(Occ,URe,UNOAO,XOne,IndN,IndX,IndAux,IGemIN, &
                        NAct,INActive,NDimX,NDim,NBasis,NInte1,     &
                        AOFile,DipFile,                             &
                        IntJFile,IntKFile,Int3File,IOrbRelax,RDM1)
implicit none
integer,intent(in)           :: IOrbRelax
integer,intent(in)           :: NAct,INActive,NDimX,NDim,NBasis,NInte1
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),IGemIN(NBasis)
character(*)                 :: AOFile,DipFile
character(*)                 :: IntJFile,IntKFile,Int3File
double precision,intent(in)  :: Occ(NBasis),XOne(NInte1)
double precision,intent(in)  :: UNOAO(NBasis,NBasis),URe(NBasis,NBasis)

double precision,intent(out),optional :: RDM1(NBasis,NBasis)

integer                      :: iunit1,iunit2
integer                      :: i,j,k,l,ij,ii,jj,ip,iq,ir,is
integer                      :: ipos,ipos1,ipos2
integer                      :: a,b,c,d,m,ac
integer                      :: NOccup0,NOccup,NVirt,NVirtOld,NOccupResp,NVirtResp
integer                      :: iOccup,iVirt
integer                      :: IGem(NBasis),pos(NBasis,NBasis),IPair(NBasis,NBasis),IPair1(NBasis,NBasis)
integer                      :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer                      :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                                tmpAV(NAct,INActive+NAct+1:NBasis),&
                                tmpIV(INActive*(NBasis-NAct-INActive))
integer                      :: limAA(2),limAI(2,1:INActive),&
                                limAV(2,INActive+NAct+1:NBasis),limIV(2)
integer                      :: IRow,ICol,NDimRed,Max_Cn,NS1
double precision             :: val,xmiu,ETot
double precision             :: UCorr(NBasis,NBasis)
double precision             :: Eps(NBasis,NBasis),CI(NBasis)
double precision             :: AuxMat(NBasis,NBasis),Gamma(NBasis,NBasis),PC(NBasis)
double precision             :: AUXM0(NBasis,NBasis),AUX2(NBasis*NBasis)

integer                      :: XInd(NBasis,NBasis),XInd1(NBasis,NBasis)
integer, allocatable         :: IndBlock(:,:)

double precision,allocatable :: ABPLUS(:,:),ABMIN(:,:)
double precision,allocatable :: AUX1(:),AuxI(:,:)
double precision,allocatable :: ints_J(:,:),ints_K(:,:)
double precision,allocatable :: ints_bi(:,:),ints_bk(:,:),ints_dl(:,:)
double precision,allocatable :: work(:),workSq(:,:)

! set dimensions
NOccup = INActive + NAct

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

XInd  = 0
IPair = 0
do ii=1,NDimX
   i = IndN(1,ii)
   j = IndN(2,ii)
   IPair(i,j) = 1
   IPair(j,i) = 1
   XInd(i,j) = IndX(ii)
   XInd(j,i) = IndX(ii)
enddo

do i=1,NBasis
  CI(i)=SQRT(Occ(i))
  if(Occ(i).Lt.0.5d0) CI(i)=-CI(i)
enddo

! construct ABPLUS(0)
allocate(ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX))

call create_blocks_ABPL0(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,pos,&
                         IGem,IndN,INActive,NAct,NBasis,NDimX)

call ABPM0_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
                IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
                IntJFile,IntKFile,0,ETot)

Eps = 0d0
!get AV
do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = pos(ip,iq)
   Eps(ip,iq) = SQRT(ABPLUS(ipos,ipos)*ABMIN(ipos,ipos))
enddo

! AB(1) PART
call AB_CAS_FOFO(ABPLUS,ABMIN,val,URe,Occ,XOne,&
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntJFile,IntKFile,0,1d0,.true.)

allocate(work(NBasis**2),ints_bi(NBasis,NBasis),ints_bk(NBasis,NBasis))
open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

Gamma = 0d0
ij = 0
do i=1,NOccup
   Gamma(i,i) = Occ(i)
enddo

AUXM0 = 0
NOccup0 = NOccup

NOccupResp = 0
NVirtResp  = 0
Do I=1,NBasis
   AUXM0(I,I)=Occ(I)
!  include active orbitals 
   If(Occ(I).Gt.0.5) NOccupResp=NOccupResp+1  
   If(Occ(I).Gt.0.5) NVirtResp=NVirtResp+1
!  exclude active orbitals
!   If(Occ(I).Gt.0.999) NOccupResp=NOccupResp+1
!   If(Occ(I).Gt.1.D-4) NVirtResp=NVirtResp+1
EndDo
NVirtResp=NVirtResp+1
Write(6,'(/,X,"NOccupResp NVirtResp in unrelaxed part: ",2I5)') NOccupResp,NVirtResp

! NOccup-NOccup block of 1-RDM
do k=1,NOccupResp
   do b=NVirtResp,NBasis

   if(IPair(b,k)==1) then

      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup0)
      do l=1,NOccupResp
         do j=1,NBasis
            ints_bk(j,l) = work((l-1)*NBasis+j)
         enddo
      enddo

      do j=1,NOccupResp
      do i=1,NOccupResp
         do a=NVirtResp,NBasis
              if(IPair(a,i)==1 .and. IPair(a,j)==1) then
                  ipos1 = XInd(a,i)
                  ipos2 = XInd(b,k)
                  Gamma(i,j) = Gamma(i,j) + ints_bk(a,j)* &
                           0.5d0* &
                           ( (CI(I)+CI(A))*(CI(K)+CI(B))*ABPLUS(ipos1,ipos2) &
                            -(CI(I)-CI(A))*(CI(K)-CI(B))* ABMIN(ipos1,ipos2) ) &
                           /(Eps(a,i)+Eps(b,k)) / (Eps(a,j)+Eps(b,k))
                 endif
           enddo
      enddo
      enddo

   endif

   enddo
enddo
!! second part
!do b=NOccup+1,NBasis
!   do k=NOccup+1,b-1
!
!   if(IPair(b,k).eq.1) then
!      ipos2 = XInd(b,k)
!
!      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup)
!      ! ints_bk
!      do l=1,NOccup
!         do j=1,NBasis
!            ints_bk(j,l) = work((l-1)*NBasis+j)
!         enddo
!      enddo
!
!      do i=1,NOccup
!         do a=NOccup+1,NBasis
!
!              if(IPair(a,i).eq.1) then
!                 ipos1 = XInd(a,i)
!                 do j=1,NOccup
!                   if(IPair(a,j).eq.1) then
!                     Gamma(i,j) = Gamma(i,j) + ints_bk(a,j)* &
!                              0.5d0* &
!                              ( (CI(I)+CI(A))*(CI(K)+CI(B))*ABPLUS(ipos1,ipos2) &
!                               -(CI(I)-CI(A))*(CI(K)-CI(B))* ABMIN(ipos1,ipos2) ) &
!                              /(Eps(a,i)+Eps(b,k)) / (Eps(a,j)+Eps(b,k))
!                    endif
!                 enddo
!              endif
!         enddo
!      enddo
!
!   endif
!
!   enddo
!enddo

! Virt-Virt block of 1-RDM

do b=NVirtResp,NBasis
  do k=1,NOccupResp

  if(IPair(b,k).eq.1) then

      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup0)
      do l=1,NOccupResp
         do j=1,NBasis
            ints_bk(j,l) = work((l-1)*NBasis+j)
         enddo
      enddo

      do i=1,NOccupResp
         do c=NVirtResp,NBasis
            if(IPair(c,i).eq.1) then
             do a=NVirtResp,NBasis
                 if(IPair(a,i).eq.1) then
                 ipos1 = XInd(a,i)
                 ipos2 = XInd(b,k)
                 Gamma(a,c) = Gamma(a,c) - ints_bk(c,i)* &
                 0.5d0* &
                 ( (CI(I)+CI(A))*(CI(K)+CI(B))*ABPLUS(ipos1,ipos2) &
                  -(CI(I)-CI(A))*(CI(K)-CI(B))* ABMIN(ipos1,ipos2) ) &
                             /(Eps(a,i)+Eps(b,k)) / (Eps(c,i)+Eps(b,k))
                 endif
             enddo
            endif
         enddo
      enddo

   endif

   enddo
enddo
!
!! second part
!do b=1,NOccup
!  do k=1,b-1
!
!  if(IPair(b,k).eq.1) then
!
!      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup)
!      ! ints_bk
!      do l=1,NOccup
!         do j=1,NBasis
!            ints_bk(j,l) = work((l-1)*NBasis+j)
!         enddo
!      enddo
!
!      do i=1,NOccup
!         do c=NOccup+1,NBasis
!            if(IPair(c,i).eq.1) then
!            do a=NOccup+1,NBasis
!                 if(IPair(a,i).eq.1) then
!                 ipos1 = XInd(a,i)
!                 ipos2 = XInd(b,k)
!                 Gamma(a,c) = Gamma(a,c) - ints_bk(c,i)* &
!                 0.5d0* &
!                 ( (CI(I)+CI(A))*(CI(K)+CI(B))*ABPLUS(ipos1,ipos2) &
!                  -(CI(I)-CI(A))*(CI(K)-CI(B))* ABMIN(ipos1,ipos2) ) &
!                             /(Eps(a,i)+Eps(b,k)) / (Eps(c,i)+Eps(b,k))
!                 endif
!            enddo
!            endif
!         enddo
!      enddo
!
!   endif
!
!   enddo
!enddo
!
!! third part
!do b=NOccup+1,NBasis
!  do k=1,NOccup
!
!  if(IPair(b,k).eq.1) then
!
!      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup0)
!      ! ints_bk
!      do l=1,NOccup0
!!      do j=1,NOccup0
!       do j=1,NBasis
!            ints_bk(j,l) = work((l-1)*NBasis+j)
!         enddo
!      enddo
!
!      do c=NOccup+1,NBasis
!         do a=NOccup+1,NBasis
!            do i=NOccup+1,min(a,c)-1
!               if(IPair(c,i)==1 .and. IPair(a,i)==1) then
!                  ipos1 = XInd(a,i)
!                  ipos2 = XInd(b,k)
!                  Gamma(a,c) = Gamma(a,c) - ints_bk(c,i)* &
!                  0.5d0* &
!                  ( (CI(I)+CI(A))*(CI(K)+CI(B))*ABPLUS(ipos1,ipos2) &
!                  -(CI(I)-CI(A))*(CI(K)-CI(B))* ABMIN(ipos1,ipos2) ) &
!                  /(Eps(a,i)+Eps(b,k)) / (Eps(c,i)+Eps(b,k))
!               endif
!            enddo
!         enddo
!      enddo
!
!   endif
!
!   enddo
!enddo

!
! occ-virt block - relaxation of orbitals
!

if(IOrbRelax==1) then

NOccup=NOccupResp
NS1=NVirtResp-1

AUX2=0.d0
NDimRed=0
XInd1=0
IPair1=0
do d=NS1+1,NBasis
   do l=1,min(NOccup,d-1) 
      if(IPair(d,l)==1) then
          NDimRed=NDimRed+1
          XInd1(d,l)=NDimRed
          IPair1(d,l)=1
       endif
   enddo
enddo
allocate(IndBlock(2,NDimRed))
NDimRed=0
do d=NS1+1,NBasis
   do l=1,min(NOccup,d-1) 
      if(IPair(d,l)==1) then
        NDimRed=NDimRed+1
        IndBlock(1,NDimRed)=d
        IndBlock(2,NDimRed)=l
      endif
   enddo
enddo

allocate(ints_dl(NBasis,NBasis))
! (FF|FO) terms
open(newunit=iunit2,file=trim(Int3File),status='OLD', &
     access='DIRECT',recl=8*NBasis*NBasis)

do d=NS1+1,NBasis
   do l=1,min(d-1,NOccup)

      read(iunit2,rec=(d+(l-1)*NBasis)) work(1:NBasis*NBasis)
      do a=1,NBasis
         do b=1,NBasis
            ints_dl(a,b) = work((a-1)*NBasis+b)
         enddo
      enddo

      do a=NOccup+1,NBasis
      do b=NOccup+1,NBasis
         if(IPair1(d,l)==1) AUX2(XInd1(d,l))=AUX2(XInd1(d,l))-Occ(l)*(gamma(a,b)-AUXM0(a,b))*2.d0*ints_dl(a,b)
         if(IPair1(a,l)==1) AUX2(XInd1(a,l))=AUX2(XInd1(a,l))+Occ(l)*(gamma(d,b)-AUXM0(d,b))*ints_dl(a,b)
      enddo
      enddo

      do j=1,NOccup
      do a=NS1+1,NBasis
      do b=NS1+1,NBasis

      if(a>j.and.b>j) then

      if(IPair(a,j)==1 .and. IPair(d,l)==1) then
      ipos1=XInd(a,j)
      ipos2=XInd(d,l)
      if(IPair1(b,j)==1) AUX2(XInd1(b,j))=AUX2(XInd1(b,j)) &
      -0.5d0* ( (CI(j)+CI(a))*(CI(l)+CI(d))*ABPLUS(ipos1,ipos2) &
              - (CI(j)-CI(a))*(CI(l)-CI(d))*ABMIN(ipos1,ipos2) ) &
              *ints_dl(a,b)/(Eps(a,j)+Eps(d,l))
      endif

      endif

      enddo
      enddo
      enddo

! (FO|FO) terms

      do i=1,NOccup
      do m=1,NOccup
         if(IPair1(d,l)==1) AUX2(XInd1(d,l))=AUX2(XInd1(d,l))-Occ(l)*(gamma(i,m)-AUXM0(i,m))*2.d0*ints_dl(i,m)
         if(IPair1(d,m)==1) AUX2(XInd1(d,m))=AUX2(XInd1(d,m))+Occ(m)*(gamma(i,l)-AUXM0(i,l))*ints_dl(i,m)
      enddo
      enddo

      do b=NS1+1,NBasis
      do i=1,NOccup
      do j=1,NOccup

      if(b>i.and.b>j) then

      if(IPair(b,i)==1 .and. IPair(d,l)==1) then
      ipos1=XInd(b,i)
      ipos2=XInd(d,l)
      if(IPair1(b,j)==1) AUX2(XInd1(b,j))=AUX2(XInd1(b,j)) &
      +0.5d0* ( (CI(i)+CI(b))*(CI(l)+CI(d))*ABPLUS(ipos1,ipos2) &
              - (CI(i)-CI(b))*(CI(l)-CI(d))*ABMIN(ipos1,ipos2) ) &
              *ints_dl(i,j)/(Eps(b,i)+Eps(d,l))
      endif

      endif

      enddo
      enddo
      enddo

    enddo
enddo
close(iunit2)

deallocate(ABPLUS,ABMIN)

allocate(AUX1(NDimRed*NDimRed))
Max_Cn=10
AUX1=0.d0
Write(6,'(/,X,"**** Expansion of the C response matrix up to n = ",I5," ****")') Max_Cn
Call CFREQPROJ(AUX1,0.d0,AUX2,1, &
   Max_Cn,XOne,URe,Occ,&
   IGem,NAct,INActive,NBasis,NInte1,IndAux,&
   0,IndBlock,IndX,NDimRed,&
   IntJFile,IntKFile)

do a=NS1+1,NBasis
   do i=1,min(NOccup,a-1)
     if(IPair(a,i)==1) then
       ! C(omega) is multiplied by 2, see the note
       val=AUX1(XInd1(a,i))*2d0
       Gamma(i,a)=Gamma(i,a)+val
       Gamma(a,i)=Gamma(a,i)+val
     endif
   enddo
enddo

close(iunit1)
deallocate(AUX1,ints_dl)

Write(6,'(/,X,"NOccupResp NVirtResp in relaxed part: ",2I5)') NOccup,NS1+1

! end of orbital relaxation 
endif

Gamma=0.5d0*(Gamma+transpose(Gamma))

deallocate(ints_bk,ints_bi)

AuxMat = Gamma
if (present(RDM1)) RDM1 = AuxMat
call Diag8(AuxMat,NBasis,NBasis,PC,work(1:NBasis))

val = 0d0
if(IOrbRelax==1) then
   write(LOUT,'(/,2x,"Orbitals Relaxed")')
else
   write(LOUT,'(/,2x,"Orbitals UnRelaxed")')
endif

write(LOUT,'(/,2x,"AC0-correlated natural occupation numbers")')
do i=NBasis,1,-1
   write(LOUT,'(X,I3,E16.6,I6)') Nbasis-i+1,PC(i)
   val = val + PC(i)
enddo
write(LOUT,'(/,1x,"Sum of AC0-correlated Occupancies: ",F5.2,/)') val

! compute transformation matrix to correlated NO's and dipole moments
Call MultpM(UCorr,AuxMat,UNOAO,NBasis)
write(LOUT,'(/,x,"Dipole moment with correlated 1-RDM")',advance="no")
if(IOrbRelax==1) then
   write(LOUT,'(1x, " (relaxed)")')
else
   write(LOUT,'(1x, " (unrelaxed)")')
endif
Call ComputeDipoleMom(UCorr,PC,NBasis,AOFile,DipFile,NBasis)

!Call PrOcc1(xmiu,val,PC,NBasis)
!Write(6,'(/,2X,"Projected Occupancies")')
!do i=NBasis,1,-1
!   write(LOUT,'(X,I3,E16.6,I6)') Nbasis-i+1,PC(i)
!enddo
!write(LOUT,'(/,x,"Dipole moment with projected 1-RDM")')
!Call ComputeDipoleMom(UCorr,PC,NBasis,NBasis)

end subroutine RDMResp_FOFO

subroutine sapt_rdm_corr(Mon,Flags,NAO,NBasis)
implicit none

type(FlagsData)     :: Flags
type(SystemBlock)   :: Mon
integer, intent(in) :: NAO, NBasis

integer :: i,j,ij
integer :: iunit
integer :: NInte1,HlpDim

double precision             :: val
double precision             :: EVal(NBasis),Eps(NBasis,NBasis)
double precision             :: CNONO(NBasis,NBasis)
double precision             :: URe(NBasis,NBasis),AuxMat(NBasis,NBasis)
double precision,allocatable :: XOne(:),work(:)

character(8) :: label
character(:),allocatable     :: aofile,dipfile
character(:),allocatable     :: onefile,twojfile,twokfile,two3file

HlpDim = max(NBasis**2,3*NBasis)
NInte1 = NBasis*(NBasis+1)/2

! save old dimensions (NDimX0, NOccup0)
Mon%NDimX0  = Mon%NDimX
Mon%NOccup0 = Mon%num0+Mon%num1

if(Mon%Monomer==1) then
   aofile   = 'AOONEINT_A'
   onefile  = 'ONEEL_A'
   twojfile = 'FFOOAA'
   twokfile = 'FOFOAA'
   two3file = 'FFFOAA'
   if(Flags%IMOLPRO==1) dipfile  = 'DIP_A'
   if(Flags%IDALTON==1) dipfile  = 'AOPROPER_A'
elseif(Mon%Monomer==2) then
   aofile   = 'AOONEINT_B'
   onefile  = 'ONEEL_B'
   twojfile = 'FFOOBB'
   twokfile = 'FOFOBB'
   two3file = 'FFFOBB'
   if(Flags%IMOLPRO==1) dipfile  = 'DIP_B'
   if(Flags%IDALTON==1) dipfile  = 'AOPROPER_B'
endif

! prepare RDM2
if(Flags%ICASSCF==1) then
   if(Mon%Monomer==1) call system('cp rdm2_A.dat rdm2.dat')
   if(Mon%Monomer==2) call system('cp rdm2_B.dat rdm2.dat')
endif

! response requires FFFO integrals
if(Flags%IOrbRelax==1) then
   call tran4_gen(NAO, &
                  NBasis,Mon%CMO, &
                  Mon%NOccup0,Mon%CMO(1:NAO,1:Mon%NOccup0),&
                  NBasis,Mon%CMO, &
                  NBasis,Mon%CMO, &
                  two3file,'AOTWOSORT')
endif

! get H0 matrix
allocate(XOne(NInte1))

open(newunit=iunit,file=onefile,access='sequential',&
     form='unformatted',status='old')

read(iunit)
read(iunit)
read(iunit) label, AuxMat
if(label=='ONEHAMIL') then
   call dgemm('T','N',NBasis,NBasis,NBasis,1d0,Mon%CMO,NBasis,AuxMat,NBasis,0d0,Eps,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,Eps,NBasis,Mon%CMO,NBasis,0d0,AuxMat,NBasis)
   ij = 0
   do j=1,NBasis
      do i=1,j
         ij = ij + 1
         XOne(ij) = AuxMat(i,j)
      enddo
   enddo
else
   write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
   stop
endif

close(iunit)

print*, 'XONE : ',norm2(XOne)

! get correlation contribution to 1-rdm
allocate(work(HlpDim))
allocate(Mon%rdm1c(NBasis,NBasis))

URe = 0d0
do i=1,NBasis
   URe(i,i) = 1d0
enddo
AuxMat = transpose(Mon%CMO)

! return MP2-like rdm1c = CAS + corr
select case(Mon%TwoMoInt)
case(TWOMO_INCORE)
   block
   integer :: NInte2,NVirt
   double precision,allocatable :: TwoMO(:),workTr(:)
   character(:),allocatable :: twofile

   NInte2 = NInte1*(NInte1+1)/2

   if(Mon%monomer==1) then
      twofile  = 'TWOMOAA'
   elseif(Mon%monomer==2) then
      twofile  = 'TWOMOBB'
   endif

   allocate(TwoMO(NInte2),workTr(NInte1))

   call SaptInter(NBasis,Mon,1)

   if(Mon%Monomer==1) call system('cp rdm2_A.dat rdm2.dat')
   if(Mon%Monomer==2) call system('cp rdm2_B.dat rdm2.dat')

   call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBasis,NInte2)
   print*, 'TWOMO',norm2(TWOMO)

   call MP2RDM(TwoMO,Eps,Mon%Occ,URe,AuxMat,XOne,  &
                Mon%IndN,Mon%IndX,Mon%IndAux,Mon%NDimX, &
                NBasis,Mon%NDimX,NInte1,NInte2,NVirt, &
                twofile,Mon%ThrVirt,.false., &
                workTr)

   call triang_to_sq2(workTr,Mon%rdm1c,NBasis)

   end block
case(TWOMO_FOFO)

    if(Mon%RDModel.lt.0) then
       Mon%RDModel = 1
       !write(lout,'(1x,a)') "Warning! Setting RDModel to 0!"
    endif

    if(Mon%RDModel==0) then

      print*, 'MP2RDM_FOFO...'
      call MP2RDM_FOFO(0d0,Eps,Mon%Occ,URe,AuxMat,XOne,&
                    Mon%IndN,Mon%IndX,Mon%IndAux,Mon%IGem,  &
                    Mon%NAct,Mon%INAct,Mon%NDimX,Mon%NDimX,NBasis,NInte1, &
                    twojfile,twokfile,Mon%ThrVirt,Mon%NVZero,Mon%IPrint,  &
                    Mon%rdm1c)

    elseif(Mon%RDModel==1) then

      print*, 'RDMResp_FOFO...'

      call RDMResp_FOFO(Mon%Occ,URe,AuxMat,XOne,&
                        Mon%IndN,Mon%IndX,Mon%IndAux,&
                        Mon%IGem,Mon%NAct,Mon%INAct, &
                        Mon%NDimX,Mon%NDimX,NBasis,NInte1, &
                        aofile,dipfile, &
                        twojfile,twokfile,two3file, &
                        Flags%IOrbRelax,Mon%rdm1c)

   endif

end select

#if RDMCORR_DEBUG > 5

   !print*, 'CAONO',norm2(Mon%CMO)
   !do i=1,NBasis
   !   write(lout,'(*(f12.8))') (AuxMat(i,j),j=1,NBasis)
   !enddo
   print*, '1-RDM MP2 unrelaxed (MO)',norm2(Mon%rdm1c)
   do i=1,NBasis
      write(lout,'(*(f12.8))') (2d0*Mon%rdm1c(i,j),j=1,NBasis)
   enddo
   !block
   !double precision :: DAO(NBasis,NBasis),S(NBasis,NBasis)

   !! try backtransformation
   !call get_one_mat('S',S,MON%Monomer,NBasis)
   !call tran_den(.true.,Mon%CMO,S,DAO,Mon%rdm1c,NBasis,NBasis)

   !print*, 'DAO',norm2(DAO)
   !do i=1,NBasis
   !   write(lout,'(*(f12.8))') (2d0*DAO(i,j),j=1,NBasis)
   !enddo
   !print*,''
   !end block
#endif

! test
!Mon%rdm1c = 0d0
!do i=1,NBasis
!   Mon%rdm1c(i,i) = Mon%Occ(i)
!enddo

CNONO = 0d0
CNONO = Mon%rdm1c

! diagonalize the new 1rdm 
call Diag8(CNONO,NBasis,NBasis,Eval,work)

do i=1,NBasis
   if (EVal(i) < 0d0) then
      print*, 'Negative occupation = ', i, EVal(i)
   endif
enddo

#if RDMCORR_DEBUG > 5
   val = 0d0
   if (Mon%IPrint > 10) write(LOUT,'(2x,"MP2",3x,"Unsorted Occupancy")')
   do i=NBasis,1,-1
      write(LOUT,'(X,I3,E16.6,I6)') i,EVal(i)
      val = val + EVal(i)
   enddo
   write(LOUT,'(/,1x,"Sum of MP2 Occupancies: ",F5.2,/)') val
#endif

call SortOcc(EVal,CNONO,NBasis)

val = 0d0
write(LOUT,'(2x,"MP2",3x,"Sorted Occupancy")')
do i=1,NBasis
   write(LOUT,'(X,I3,E16.6,I6)') i,EVal(i)
   val = val + EVal(i)
enddo
write(LOUT,'(/,1x,"Sum of MP2 Occupancies: ",F5.2,/)') val

! canonicalize
print*, 'WARNING! Add canonicalization for 2nd order!'

! save new orbital coefficients 
! M%CMO = C(AO,NO)
! C(AO,NO).C(NOMP2,NO)^T
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,Mon%CMO,NBasis,CNONO,NBasis,0d0,AuxMat,NBasis)

#if RDMCORR_DEBUG > 5
   print*, 'CAONO:before',norm2(Mon%CMO)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (Mon%CMO(i,j),i=1,NBasis)
   enddo
   print*, 'CAONO:after',norm2(AuxMat)
   do j=1,NBasis
      write(LOUT,'(*(f13.8))') (AuxMat(i,j),i=1,NBasis)
   enddo
#endif

! save old NO to new NO coeffs
allocate(Mon%CNONO(NBasis,NBasis))
Mon%CNONO = transpose(CNONO)

! save new occupation numbers
allocate(Mon%Occ0(NBasis))
!print*, 'test-1: keep CAS occupation numbers'
Mon%Occ0 = Mon%Occ
Mon%Occ  = EVal

!block
!print*, 'test-1: keep CAS occupation numbers for NOccup'
!Mon%Occ(1:Mon%NOccup0) = Mon%Occ0(1:Mon%NOccup0)
!end block

!call ComputeDipoleMom(transpose(Mon%CMO),Mon%Occ0,Mon%NOccup0,aofile,dipfile,NBasis)

! save new C(AO,NO) coefficients
!print*, 'test-2: save old orbitals!'
Mon%CMO = AuxMat

!write(lout,'(1x,a)') "DIPOLE MOMENT WITH MP2-LIKE 1RDM:"
!call ComputeDipoleMom(transpose(Mon%CMO),Mon%Occ,NBasis,aofile,dipfile,NBasis)

deallocate(XOne)

end subroutine sapt_rdm_corr

subroutine prepare_RDM2corr(Mon,NBasis,ver)
!
! prepare RDM2val(NOccup,NOccup,NOccup,NOccup) 
! saved as RDM2val 
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: ver
integer,intent(in) :: NBasis

integer :: i,j,k,l
integer :: ip,iq,ir,is
integer :: NOccup,NOccup0
double precision :: xnorm,refnorm
double precision,allocatable :: work(:,:,:,:)

NOccup0 = Mon%NOccup0
NOccup = Mon%num0+Mon%num1

print*, 'NOccup0 = ',NOccup0
print*, 'NOccup  = ',NOccup
print*, 'ver = ', ver

if(ver==0) then
   if(allocated(Mon%RDM2val)) then
      ! save old RDM2val
      allocate(work(NOccup0,NOccup0,NOccup0,NOccup0))
      work = 0d0
      work = Mon%RDM2val
   else
      print*, 'ERROR in prepare_RDM2corr: RDM2val not available!'
      stop
   endif
endif

if(allocated(Mon%RDM2val)) deallocate(Mon%RDM2val)
allocate(Mon%RDM2val(NOccup,NOccup,NOccup,NOccup))

if(ver==0) then
   ! keep CAS RDM2

   ! transform RDM2val to new NOs
   block
   double precision :: val
   double precision,allocatable :: RDM2(:,:,:,:),Aux(:,:,:,:)
   integer :: NVirt0, NVirt

   !print*, 'test off-NOccup block'
   !NVirt0 = NBasis - NOccup0
   !NVirt  = NBasis - NOccup

   !val = 0d0
   !do j=1,NOccup
   !do i=1,NVirt0
   !! CNONO(NO,NOMP2) / NOccup0,NOccup
   !val = val + Mon%CNONO(NOccup0+i,j)**2
   !enddo
   !enddo
   !print*, 'Blocks V1-I2A2 = ',val
   !val = 0d0
   !do j=1,NVirt
   !do i=1,NOccup0
   !val = val + Mon%CNONO(i,NOccup+j)**2
   !enddo
   !enddo
   !print*, 'Blocks I1A1-V2 = ',val

   ! try full tran
   print*, 'TRY 4-indx transformation of RDM2CAS...'
   allocate(RDM2(NBasis,NBasis,NBasis,NBasis))
   allocate(Aux(NBasis,NBasis,NBasis,NBasis))
   RDM2=0d0
   do l=1,NOccup0
      do k=1,NOccup0
         do j=1,NOccup0
            do i=1,NOccup0
               RDM2(i,j,k,l) = work(i,j,k,l)
            enddo
         enddo
      enddo
   enddo
   call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, RDM2, NBasis, Mon%CNONO, NBasis, 0.d0, Aux, NBasis**3)
   call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, Aux,  NBasis, Mon%CNONO, NBasis, 0.d0, RDM2,NBasis**3)
   call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, RDM2, NBasis, Mon%CNONO, NBasis, 0.d0, Aux, NBasis**3)
   call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, Aux,  NBasis, Mon%CNONO, NBasis, 0.d0, RDM2,NBasis**3)

   Mon%RDM2val = 0d0
   do l=1,NOccup
      do k=1,NOccup
         do j=1,NOccup
            do i=1,NOccup
               Mon%RDM2val(i,j,k,l) = RDM2(i,j,k,l)
            enddo
         enddo
      enddo
   enddo
   deallocate(Aux,RDM2)
   end block

   !do l=1,NOccup0
   !   do k=1,NOccup0
   !      do j=1,NOccup0
   !         do i=1,NOccup0
   !            Mon%RDM2val(i,j,k,l) = work(i,j,k,l) 
   !         enddo
   !      enddo
   !   enddo
   !enddo

   print*, '2RDM from CAS, dimension extended from ', NOccup0, 'to', NOccup
   print*, '2-RDM',norm2(work),norm2(Mon%RDM2val)
   deallocate(work)

elseif(ver==1) then

   ! Hartree-Fock
   ! Gamma(prqs) = 2*np*nq \delta_pr \delta_qs - n_p*nq \delta_ps \delta_qr
   Mon%RDM2val = 0d0
   ! Coulomb
   do iq=1,NOccup
      do ip=1,NOccup
         Mon%RDM2val(ip,ip,iq,iq) = Mon%RDM2val(ip,ip,iq,iq) + 2d0*Mon%Occ(ip)*Mon%Occ(iq)
      enddo
   enddo
   ! exchange
   do iq=1,NOccup
      do ip=1,NOccup
         Mon%RDM2val(ip,iq,iq,ip) = Mon%RDM2val(ip,iq,iq,ip) - Mon%Occ(ip)*Mon%Occ(iq)
      enddo
   enddo
   ! test
   !do is=1,NOccup
   !do iq=1,NOccup
   !do ir=1,NOccup
   !do ip=1,NOccup
   !   if ((ip==ir).and.(iq==is)) then
   !      Mon%RDM2val(ip,ir,iq,is) = Mon%RDM2val(ip,ir,iq,is) + 2d0*Mon%Occ(ip)*Mon%Occ(iq)
   !   endif
   !   if ((ip==is).and.(iq==ir)) then
   !      !Mon%RDM2val(ip,ir,iq,is) = Mon%RDM2val(ip,ir,iq,is) - Mon%Occ(iq)
   !      mon%RDM2val(ip,ir,iq,is) = Mon%RDM2val(ip,ir,iq,is) - Mon%Occ(ip)*Mon%Occ(iq)
   !   endif
   !enddo
   !enddo
   !enddo
   !enddo

   print*, '2RDM non-cumulant: 2 n_p n_q - n_p n_q'
   print*, '2-RDM',norm2(Mon%RDM2val)

elseif(ver==2) then
   ! DMFT
   ! Gamma(prqs) = 2*np*nq \delta_pr \delta_qs - F_pq \delta_ps \delta_qr
   !       1122
   Mon%RDM2val = 0d0
   ! Coulomb (nc part)
   do ip=1,NOccup
      do iq=1,NOccup
         Mon%RDM2val(ip,ip,iq,iq) = Mon%RDM2val(ip,ip,iq,iq) + 2d0*Mon%Occ(ip)*Mon%Occ(iq)
      enddo
   enddo
   ! exchange-corr
   do ip=1,NOccup
      do iq=1,NOccup
         Mon%RDM2val(ip,iq,iq,ip) = Mon%RDM2val(ip,iq,iq,ip) - sqrt(Mon%Occ(ip)*Mon%Occ(iq))
      enddo
   enddo

   print*, '2RDM-corr with BB functional: 2 n_p n_q - sqrt(n_p n_q)'
   print*, '2-RDM',norm2(Mon%RDM2val)

endif

! check norm
refnorm = Mon%XELE*(2d0*Mon%XELE-1)

xnorm = 0d0
do i=1,NOccup
   do j=1,NOccup
      xnorm = xnorm + Mon%RDM2val(i,i,j,j)
   enddo
enddo

if(mon%monomer==1) write(lout,'(/1x,a,i3)') 'Monomer A / NOccup =',NOccup
if(mon%monomer==2) write(lout,'(/1x,a,i3)') 'Monomer B / NOccup =',NOccup
write(lout,'(1x,a,f12.6)',advance="no") '2-RDM2 norm = ', xnorm
write(lout,'(1x,a,f8.3,a)') '(reference =', refnorm, ')'

!if(abs(refnorm-xnorm).gt.1d-5) then
!   Mon%RDM2val = Mon%RDM2val * refnorm / xnorm
!   xnorm = 0d0
!   do i=1,NOccup
!      do j=1,NOccup
!         xnorm = xnorm + Mon%RDM2val(i,i,j,j)
!      enddo
!   enddo
!   write(lout,'(1x,a,f12.6)') 'RDM2 renormalized!',xnorm
!endif

!block
!integer :: ip,iq,ir,is
!double precision :: tmp,ntmp,sumOcc
!double precision :: tst(NBasis,NBasis)
!! test RDM2val
!  tst = 0d0
!  do ip=1,NOccup 
!  do is=1,NOccup
!  do iq=1,NOccup
!     !tst(iq,is) = tst(iq,is) + Mon%RDM2val(iq,is,ip,ip)
!     tst(iq,is) = tst(iq,is) + Mon%RDM2val(ip,ip,iq,is)
!  enddo
!  enddo
!  enddo
!  !tst = tst / (2d0*Mon%XELE-1)
!  tmp = norm2(tst)
!  print*, 'test sum rule',tmp
!  xnorm = 0d0
!  sumOcc = 0d0
!  ntmp=0d0
!  do ip=1,NBasis
!  !   write(lout,'(*(f13.8))') (tst(iq,is),iq=1,NBasis)
!     xnorm = xnorm + tst(ip,ip)**2  
!     sumOcc = sumOcc + tst(ip,ip)
!     ntmp = 2d0*Mon%XELE*Mon%Occ(ip)-Mon%Occ(ip)**2
!     write(lout,'(i3,2f12.6)') ip, tst(ip,ip),ntmp
!     !write(lout,'(i3,2f12.6)') ip, tst(ip,ip),Mon%Occ(ip)
!  enddo
!  xnorm = sqrt(xnorm)
!  print*, 'NORM TEST', abs(tmp-xnorm)
!  print*, 'Summ Occ ', sumOcc
!
!end block

end subroutine prepare_RDM2corr

end module rdmcorr

