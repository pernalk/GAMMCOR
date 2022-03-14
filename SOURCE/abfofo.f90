module abfofo
!use types,only : EblockData
use blocktypes
use tran
use print_units

implicit none

contains

subroutine AB_CAS_FOFO(ABPLUS,ABMIN,ETot,URe,Occ,XOne, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     IntJFile,IntKFile,ACAlpha,AB1)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
!
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! THE FOLLOWING SYMMETRY IS ASSUMED
! RDM2(ij,kl) = RDM2(kl,ij)
! ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
! COULOMB INTEGRALS ARE READ FROM IntJFile IN (FF|OO) FORMAT
! EXCHANGE INTEGRALS ARE READ FROM IntKFile IN (FO|FO) FORMAT
 implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(in)  :: ACAlpha
logical,intent(in) :: AB1

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit1,iunit2,ios
integer :: switch
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,HNOCoef,val
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
ETot   = 0

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

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit1,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit1,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit1)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

if(AB1) then
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
         if(IGem(i)==IGem(j)) HNO(i,j) = 0
      enddo
   enddo
else
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
      enddo
   enddo
endif

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

if(AB1) then
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 0
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
else
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 1
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
endif

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxI  = 0
AuxIO = 0
WMAT  = 0

if(AB1) then
   HNOCoef = -1
else
   HNOCoef = 1 - ACAlpha
endif

if(AB1) then
   switch=1
else
   switch=0
endif

call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
             RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
             INActive,NOccup,NDim,NDimX,NBasis,NInte1,IntJFile,IntKFile,ACAlpha,switch,ETot)

write(LOUT,'(1x,a,5x,f15.8)') "CASSCF Energy (w/o ENuc)", ETot

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = IndX(ICol)
   if(ir/=is) then

      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = IndX(IRow)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
! hererXXX
!         if(ipq.eq.irs) then
!           ABPLUS(ipq,irs) = ABPLUS(ipq,irs)+(Occ(IP)-Occ(IQ))*(Occ(IP)+1.D0-Occ(IQ))*0.25
!           ABMIN(ipq,irs)  = ABMIN(ipq,irs)+(Occ(IP)-Occ(IQ))*(Occ(IP)+1.D0-Occ(IQ))*0.25
!         endif

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo

   endif
enddo

if(AB1) then
   ! symmetrize AB1
   call sq_symmetrize(ABPLUS,NDimX)
   call sq_symmetrize(ABMIN,NDimX)
endif

!print*, "AB-my",norm2(ABPLUS),norm2(ABMIN)

!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!!!$
!!!!$write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!!!!$
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxI=transpose(AuxI)
!call sq_to_triang2(AuxI,work1,NBasis)
!write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!AuxIO=transpose(AuxIO)
!call sq_to_triang2(AuxIO,work1,NBasis)
!write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!
!write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)
deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_FOFO

subroutine MP2RDM_FOFO(PerVirt,Eps,Occ,URe,UNOAO,XOne,IndN,IndX,IndAux,IGemIN, &
                       NAct,INActive,NDimX,NDim,NBasis,NInte1,     &
                       IntJFile,IntKFile,ThrVirt,NVZero,IPrint)
implicit none

integer,intent(in)           :: NAct,INActive,NDimX,NDim,NBasis,NInte1
integer,intent(in)           :: IPrint
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IndAux(NBasis),IGemIN(NBasis)
integer,intent(out)          :: NVZero
double precision,intent(in)  :: PerVirt,ThrVirt
double precision,intent(in)  :: Occ(NBasis),XOne(NInte1),UNOAO(NBasis,NBasis)
double precision,intent(out) :: Eps(NBasis,NBasis)
character(*)                 :: IntJFile,IntKFile

integer                      :: iunit1,iunit2
integer                      :: i,j,k,l,a,b,c,ij,ac,ii,jj,ip,iq,ir,is,ipos
!integer                      :: NOccup,NVirt,NVirtOld,NVZero,iOccup,iVirt
integer                      :: NOccup,NVirt,NVirtOld,iOccup,iVirt
integer                      :: IGem(NBasis),pos(NBasis,NBasis),IPair(NBasis,NBasis)
integer                      :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer                      :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                                tmpAV(NAct,INActive+NAct+1:NBasis),&
                                tmpIV(INActive*(NBasis-NAct-INActive))
integer                      :: limAA(2),limAI(2,1:INActive),&
                                limAV(2,INActive+NAct+1:NBasis),limIV(2)
double precision             :: PerThr
double precision             :: val,ETot
double precision             :: URe(NBasis,NBasis)
double precision             :: AuxMat(NBasis,NBasis),&
                                Gamma(NInte1),PC(NBasis),Fock(NBasis*NBasis)
double precision,allocatable :: ABPLUS(:,:),ABMIN(:,:)
double precision,allocatable :: AuxI(:,:),Aux2(:,:)
double precision,allocatable :: ints_J(:,:),ints_K(:,:),   &
                                ints_bi(:,:),ints_bk(:,:), &
                                work(:),workTr(:),workSq(:,:)
integer          :: EndVirt,NBasisNew
character(100)   :: num1char
!double precision :: xnum1,PerThr

! use PerVirt/100
PerThr = PerVirt

!! TEST!
!if(COMMAND_ARGUMENT_COUNT()==0) then
!   write(LOUT,*)'ERROR, COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
!   stop
!endif
!
!CALL GET_COMMAND_ARGUMENT(1,num1char)
!READ(num1char,*)xnum1
!PerThr = xnum1

! set dimensions
NOccup   = INActive + NAct
NVirtOld = NBasis - NOccup

!print*, 'NOccup',NOccup,INActive,NAct

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

IPair = 0
do ii=1,NDimX
   i = IndN(1,ii)
   j = IndN(2,ii)
   IPair(i,j) = 1
   IPair(j,i) = 1
enddo

! construct ABPLUS(0)
allocate(ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX))

call create_blocks_ABPL0(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,pos,&
                         IGem,IndN,INActive,NAct,NBasis,NDimX)

call ABPM0_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
                IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
                IntJFile,IntKFile,ETot)

Eps = 0
!get AV
do ir=NOccup+1,NBasis
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = pos(ip,iq)
      Eps(ip,iq) = ABPLUS(ipos,ipos)
   enddo
enddo
!get IV
do ii=1,nIV
   i = tmpIV(ii)
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = pos(ip,iq)
   Eps(ip,iq) = ABPLUS(ipos,ipos)
enddo

deallocate(ABMIN,ABPLUS)

allocate(work(NBasis**2),ints_bi(NBasis,NBasis),ints_bk(NBasis,NBasis))
open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

Gamma = 0
ij = 0
do j=1,NOccup
   do i=1,j
      ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
      if(i==j) Gamma(ij) = Occ(i)
   enddo
enddo

do k=1,NOccup
   do b=NOccup+1,NBasis

      read(iunit1,rec=(b+(k-1)*NBasis)) work(1:NBasis*NOccup)
      ! ints_bk
      do l=1,NOccup
         do j=1,NBasis
            ints_bk(j,l) = work((l-1)*NBasis+j)
         enddo
      enddo
      ints_bk(:,NOccup+1:NBasis) = 0

      ij = 0
      do i=1,NOccup
         read(iunit1,rec=(b+(i-1)*NBasis)) work(1:NBasis*NOccup)
         ! ints_bi
         do l=1,NOccup
            do ii=1,NBasis
               ints_bi(ii,l) = work((l-1)*NBasis+ii)
            enddo
         enddo
         ints_bi(:,NOccup+1:NBasis) = 0

         do a=NOccup+1,NBasis
            do j=1,i
               ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
               Gamma(ij) = Gamma(ij) - ints_bk(a,j)*(2d0*ints_bk(a,i)-ints_bi(a,k)) / &
                           (Eps(a,i)+Eps(b,k)) / (Eps(a,j)+Eps(b,k))
            enddo
         enddo

         do c=NOccup+1,NBasis
            do a=NOccup+1,c
               ac = (c*(c-1))/2 + a
               Gamma(ac) = Gamma(ac) + ints_bk(c,i)*(2d0*ints_bk(a,i)-ints_bi(a,k)) / &
                           (Eps(a,i)+Eps(b,k)) / (Eps(c,i)+Eps(b,k))
            enddo
         enddo

      enddo

   enddo
enddo

close(iunit1)
deallocate(ints_bi,ints_bk)

!print*, 'Gamma-my',norm2(Gamma)
!do j=1,NBasis
!   do i=1,j
!      ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
!      print*, i,j,Gamma(ij)
!   enddo
!enddo

call triang_to_sq2(Gamma,AuxMat,NBasis)
call Diag8(AuxMat,NBasis,NBasis,PC,work(1:NBasis))

val = 0d0
if(IPrint>2) write(LOUT,'(2x,"MP2",3x,"Unsorted Occupancy")')
do i=NBasis,1,-1
   if(IPrint>2) write(LOUT,'(X,I3,E16.6,I6)') i,PC(i)
   val = val + PC(i)
enddo
write(LOUT,'(/,1x,"Sum of MP2 Occupancies: ",F5.2,/)') val

! SET Occ of the obitals belonging to NOccup to 1 before sorting
! to make sure that after sorting all NOccup orbitals come first
do i=1,NBasis
   iOccup = 0
   iVirt  = 0
   do j=1,NBasis
      if(AuxMat(i,j)/=0d0.and.j<=NOccup) iOccup = 1
      if(AuxMat(i,j)/=0d0.and.j>NOccup)  iVirt  = 1
   enddo

   if(iOccup==1.and.iVirt==0) PC(i)=1d0
   if(iOccup*iVirt==1.or.iOccup+iVirt==0) then
      write(LOUT,'(1x,a)') 'MP2 1RDM is messed up. Quitting.'
      stop
   endif
enddo

call SortP(PC,AuxMat,NBasis)

! transformation matrix to new virtual orbitals
do i=1,NBasis
   do k=1,NOccup
      AuxMat(k,i) = 0d0
      if(i==k) AuxMat(k,i) = 1d0
   enddo
enddo

NVZero = 0
NVirtOld = NBasis - NOccup

!! occupancy threshold
!do i=NOccup+1,NBasis
!   if(PC(i)<=ThrVirt) then
!      if(IPrint>2) write(LOUT,'(1x,"Virtual MP2 orbital no",i3," of the occup", e14.6," removed")') i,PC(i)
!      AuxMat(i,:) = 0d0
!      NVZero = NVZero + 1
!   endif
!enddo
!write(LOUT,'(1x,"Threshold for virt occup number :",e16.6)') ThrVirt
!write(LOUT,'(/,1x,"Total number of removed orbitals: ",i3)')   NVZero
!!
! percentage threshold
PerThr = PerThr/100d0
val    = floor(dble(NVirtOld)*PerThr)
EndVirt = NBasis - int(val) + 1
do i=NBasis,EndVirt,-1
   if(IPrint>2) write(LOUT,'(1x,"Virtual MP2 orbital no",i3," of the occup", e14.6," removed")') i,PC(i)
   AuxMat(i,:) = 0d0
   NVZero = NVZero + 1
enddo
write(LOUT,'(1x,"Threshold for %virt occup number :",f16.6)') PerThr*100
write(LOUT,'(/,1x,"Total number of removed orbitals: ",i3)')   NVZero

NVirt    = NBasis - NOccup - NVZero
val = (1d0 - dble(NVirt)/dble(NVirtOld) ) * 100d0
write(LOUT,'(1x,"NVirt/NVirtOld",18x,": ",i3," / ",i3)')  NVirt, NVirtOld
write(LOUT,'(1x,"Percent of removed orbitals     : ",f6.2)') Val

NVirt = NBasis - NOccup

! If the new are orbitals are to be used in AC0 they must be
! cannonicalized

allocate(ints_J(NBasis,NBasis),ints_K(NBasis,NBasis))
allocate(AuxI(NVirt,NVirt),Aux2(NVirt,NVirt))
open(newunit=iunit1,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NBasis)
open(newunit=iunit2,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

AuxI = 0
Aux2 = 0

ii = 0
jj = 0
do i=1,NVirt
   do j=1,NVirt
      ii = NOccup + i
      jj = NOccup + j
      ij=(max(ii,jj)*(max(ii,jj)-1))/2+min(ii,jj)
      AuxI(i,j) = AuxI(i,j) + XOne(ij)
      Aux2(i,j) = AuxMat(ii,jj)
   enddo
enddo

do k=1,NOccup

   read(iunit1,rec=k+(k-1)*NOccup) work(1:NBasis*NBasis)
   do j=1,NBasis
      do i=1,NBasis
         ints_J(i,j) = work((j-1)*NBasis+i)
      enddo
   enddo

   work = 0
   do jj=NOccup+1,NBasis
      j = jj - NOccup
      read(iunit2,rec=(jj+(k-1)*NBasis)) work(1:NBasis*NOccup)
      ! ints_bk
      do l=1,NOccup
         do i=1,NBasis
            ints_K(i,l) = work((l-1)*NBasis+i)
         enddo
      enddo
      ints_K(:,NOccup+1:NBasis) = 0

      do ii=NOccup+1,NBasis
         i = ii - NOccup
         AuxI(i,j) = AuxI(i,j) + Occ(k)*(2d0*ints_J(jj,ii)-ints_K(ii,k))
      enddo

   enddo
enddo

close(iunit2)
close(iunit1)

Fock = 0
allocate(workTr(NVirt*(NVirt+1)/2))
! transform Fock to the new basis
workTr = 0
call sq_to_triang2(AuxI,workTr,NVirt)
call MatTr(workTr,Aux2,NVirt)
! diagonalize and sort
call triang_to_sq(workTr,Fock,NVirt)
call Diag8(Fock,NVirt,NVirt,PC,work(1:NVirt))
call SortF(PC,Fock,NVirt)
deallocate(workTr)

!print*, 'Fock,PC',norm2(Fock),norm2(PC(1:NVirt))
!print*, 'Aux2',norm2(AuxI),norm2(Aux2)
!print*, 'AUXM-my',norm2(AuxMat)

do i=1,NVirt
   do j=1,NVirt
      ii = NOccup + i
      jj = NOccup + j
      URe(ii,jj)=Fock((j-1)*NVirt+i)
   endDo
endDo

! Set elements corresponding to the removed orbitals to zero
NVirt=NBasis-NOccup-NVZero

do i=NOccup+NVirt+1,NBasis
   do j=1,NBasis
      URe(i,j) = 0d0
      URe(j,i) = 0d0
   enddo
enddo

! get NO --> NOMP2(canon) tran mat
call MultpM(Eps,URe,AuxMat,NBasis)
!Print*, 'Eps-MY',norm2(Eps)

write(LOUT,'(/,1x,"Integral transformation in progress ... ",/)')
call MatTr(XOne,Eps,NBasis)

! get AO --> NOMP2(canon) tran mat
allocate(workSq(NBasis,NBasis))
workSq=transpose(UNOAO)

call dgemm('N','T',NBasis,NBasis,NBasis,1d0,workSq,NBasis,Eps,NBasis,0d0,AuxMat,NBasis)
NBasisNew = NBasis-NVZero
AuxMat(1:NBasis,NOccup+NVirt+1:NBasis)=0d0
!print*, 'AONOMP2',norm2(AuxMat)
deallocate(workSq)

! transform J and K
call tran4_gen(NBasis,&
               NOccup,AuxMat(1:NBasis,1:NOccup),&
               NOccup,AuxMat(1:NBasis,1:NOccup),&
               NBasis,AuxMat,&
               NBasis,AuxMat,&
               IntJFile,'AOTWOSORT')
call tran4_gen(NBasis,&
               NBasis,AuxMat,&
               NOccup,AuxMat(1:NBasis,1:NOccup),&
               NBasis,AuxMat,&
               NOccup,AuxMat(1:NBasis,1:NOccup),&
               IntKFile,'AOTWOSORT')
!
! restore URe
URe = 0
do i=1,NBasis
   URe(i,i) = 1d0
enddo

! after "truncating" of the virtual space, find a new set of accepted pairs
call AcceptPair(IndN,IndX,NDimX,IndAux,Occ,NOccup,NVirt,NBasis,IPrint)

deallocate(AuxI,Aux2)
deallocate(ints_J,ints_K,work)

end subroutine MP2RDM_FOFO

subroutine ABPM0_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
                      IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
                      IntJFile,IntKFile,ETot)
implicit none

integer,intent(in)           :: NAct,INActive,NDimX,NBasis,NDim,NInte1
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
character(*)                 :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(out),optional :: ETot

integer          :: i,j,k,l
integer          :: ip,iq,ir,is,ipq,irs,ICol,IRow
integer          :: iunit,ios
integer          :: NOccup,NRDM2Act
integer          :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: val,AuxVal,EnOne
double precision :: AuxCoeff(3,3,3,3),C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),&
                    AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
integer,external          :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
EnOne = 0
if(present(ETot)) ETot = 0

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

allocate(work1(NBasis**2),work2(NBasis**2))
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
!write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

deallocate(work1)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
EnOne = EnOne + 2*val

! DEBUG
!print*, 'XOne:', norm2(XOne)
!print*, 'RDM2val:', norm2(RDM2val)
!print*, 'Occ:', norm2(Occ)
!print*, 'HNO:', norm2(HNO)

!write(LOUT,*) 'ONE ELECTRON ENERGY:', EnOne
if(present(ETot)) ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

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

call create_pos_ABPL0(pos,IGem,IndN,INActive,NAct,NBasis,NDimX)

AuxI  = 0
AuxIO = 0
WMAT  = 0

if(present(ETot)) then
   call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
                RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                INActive,NOccup,NDimX,NDimX,NBasis,NInte1,IntJFile,IntKFile,0d0,2,ETot)
else
   call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
                RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                INActive,NOccup,NDimX,NDimX,NBasis,NInte1,IntJFile,IntKFile,0d0,2)
endif

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)

      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

! hererXXX
!         if(ipq.eq.irs) then
!           ABPLUS(ipq,irs) = ABPLUS(ipq,irs)+(Occ(IP)-Occ(IQ))*(Occ(IP)+1.D0-Occ(IQ))*0.25
!           ABMIN(ipq,irs)  = ABMIN(ipq,irs)+(Occ(IP)-Occ(IQ))*(Occ(IP)+1.D0-Occ(IQ))*0.25
!           endif
!         endif

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo
enddo

deallocate(work2)
deallocate(RDM2val)

end subroutine ABPM0_FOFO

subroutine create_pos_ABPL0(pos,IGem,IndN,INActive,NAct,NBasis,NDimX)
implicit none

integer,intent(in)  :: INActive,NAct,NBasis,NDimX
integer,intent(in)  :: IGem(NBasis),IndN(2,NDimX)
integer,intent(out) :: pos(NBasis,NBasis)

integer :: NOccup
integer :: i,ii,ip,iq,ir,is,ipos
integer :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
           tmpAV(NAct,INActive+NAct+1:NBasis),&
           tmpIV(INActive*(NBasis-NAct-INActive))
integer :: limAA(2),limAI(2,1:INActive),&
           limAV(2,INActive+NAct+1:NBasis),limIV(2)

NOccup = NAct + INActive

call create_AC0_blocks_1(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                      IGem,IndN,INActive,NAct,NBasis,NDimX)
call create_AC0_blocks_2(pos,nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,IndN,INActive,NAct,NBasis,NDimX)

end subroutine create_pos_ABPL0

subroutine create_blocks_ABPL0(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                               limAA,limAI,limAV,limIV,pos,&
                               IGem,IndN,INActive,NAct,NBasis,NDimX)
implicit none

integer,intent(in)  :: INActive,NAct,NBasis,NDimX
integer,intent(in)  :: IGem(NBasis),IndN(2,NDimX)
integer,intent(out) :: pos(NBasis,NBasis)

integer :: NOccup
integer :: i,ii,ip,iq,ir,is,ipos
integer :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
           tmpAV(NAct,INActive+NAct+1:NBasis),&
           tmpIV(INActive*(NBasis-NAct-INActive))
integer :: limAA(2),limAI(2,1:INActive),&
           limAV(2,INActive+NAct+1:NBasis),limIV(2)

NOccup = NAct + INActive

call create_AC0_blocks_1(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                      IGem,IndN,INActive,NAct,NBasis,NDimX)
call create_AC0_blocks_2(pos,nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,IndN,INActive,NAct,NBasis,NDimX)

end subroutine create_blocks_ABPL0

subroutine create_AC0_blocks_1(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                               IGem,IndN,INActive,NAct,NBasis,NDimX)
implicit none

integer,intent(in)  :: INActive,NAct,NBasis,NDimX
integer,intent(in)  :: IndN(2,NDimX),IGem(NBasis)
integer,intent(out) :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer,intent(out) :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                       tmpAV(NAct,INActive+NAct+1:NBasis),&
                       tmpIV(INActive*(NBasis-NAct-INActive))
integer :: NOccup
integer :: i,ip,iq

nAA   = 0
nAI   = 0
nAV   = 0
nIV   = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

end subroutine create_AC0_blocks_1

subroutine create_AC0_blocks_2(pos,nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                               limAA,limAI,limAV,limIV,IndN,INActive,NAct,NBasis,NDimX)
implicit none

integer,intent(in)  :: INActive,NAct,NBasis,NDimX
integer,intent(in)  :: IndN(2,NDimX)
integer,intent(in)  :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer,intent(in)  :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                       tmpAV(NAct,INActive+NAct+1:NBasis),&
                       tmpIV(INActive*(NBasis-NAct-INActive))
integer,intent(out) :: pos(NBasis,NBasis)

integer :: NOccup
integer :: i,ii,ip,iq,ir,is,ipos
integer :: limAA(2),limAI(2,1:INActive),&
           limAV(2,INActive+NAct+1:NBasis),limIV(2)

NOccup = NAct + INActive

limAA = 0
limAI = 0
limAV = 0
limIV = 0

pos  = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos
   enddo
   limAI(2,is) = ipos
enddo
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

end subroutine create_AC0_blocks_2

subroutine pack_Eblock(ABPLUS,ABMIN,nVal,lim1,lim2,tmpMat,Eblock,NoSt,NDimX)
!
! 1) solve the ERPA eigen problem for a given block
! 2) save the solution in Eblock type
!
implicit none

integer,intent(in) :: NoSt,NDimX
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

  allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
  allocate(ABP(B%n,B%n),ABM(B%n,B%n))

  ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
  ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
  if(NoSt==1) then
     call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
  elseif(NoSt>1) then
     call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
  endif

  deallocate(ABM,ABP)

end associate

end subroutine pack_Eblock

subroutine dump_Eblock(Eblock,EblockIV,Occ,IndN,nblk,NBasis,NDimX,xy0file)
!
! save blocks to XY0FILE
! XTilde (matX) and YTilde (matY) are saved:
! XTilde = X / 2 / ( c_p + c_q ) - Y / 2 / ( c_p - c_q)
! YTilde = X / 2 / ( c_p + c_q ) + Y / 2 / ( c_p - c_q)
!
implicit none

integer,intent(in)           :: nblk,NBasis,NDimX
integer,intent(in)           :: IndN(2,NDimX)
character(*),intent(in)      :: xy0file
double precision,intent(in)  :: Occ(NBasis)
type(EBlockData),intent(in)  :: EBlock(nblk),EBlockIV
type(EBlockData),allocatable :: SBlock(:)
type(EBlockData)             :: SBlockIV

integer                      :: iunit
integer                      :: i,j,k,ii,ipos,iblk,ip,iq
double precision             :: C(NBasis),fac,valY,valX
double precision,allocatable :: EigY(:,:),EigX(:,:),Eig(:)
double precision,allocatable :: work(:,:)

fac = 1d0/sqrt(2d0)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(SBlock(nblk))

! repack
do iblk=1,nblk
   associate(B => Eblock(iblk),&
             A => Sblock(iblk))

     A%n  = B%n
     A%l1 = B%l1
     A%l2 = B%l2
     allocate(A%pos(A%n),A%vec(A%n))
     A%pos(1:A%n)=B%pos(1:B%n)
     A%vec(1:A%n)=B%vec(1:B%n)

     allocate(A%matX(A%n,A%n),A%matY(A%n,A%n))
     do i=1,B%n
        ipos = B%pos(i)
        ip = IndN(1,ipos)
        iq = IndN(2,ipos)

        valX = 1d0/(C(ip)+C(iq))
        valY = 1d0/(C(ip)-C(iq))
        A%matX(i,1:B%n) = 0.5d0*(B%matX(i,1:B%n)*valX - B%matY(i,1:B%n)*valY)
        A%matY(i,1:B%n) = 0.5d0*(B%matX(i,1:B%n)*valX + B%matY(i,1:B%n)*valY)

     enddo

   end associate
enddo
! repack IV
associate(B => EblockIV, &
          A => SBlockIV )

  A%n  = B%n
  A%l1 = B%l1
  A%l2 = B%l2
  allocate(A%pos(A%n))
  A%pos(1:A%n) = B%pos(1:B%n)

  allocate(A%vec(A%n),A%matX(A%n,1),A%matY(A%n,1))
  do i=1,A%n
     ipos = B%pos(i)
     ip = IndN(1,ipos)
     iq = IndN(2,ipos)

     valX = 1d0/(C(ip)+C(iq))
     valY = 1d0/(C(ip)-C(iq))

     A%matX(i,1) = 0.5d0*fac*(valX-valY)
     A%matY(i,1) = 0.5d0*fac*(valX+valY)

     A%vec(i) = B%vec(i)
  enddo

end associate

! dump to a file
open(newunit=iunit,file=xy0file,form='unformatted')
write(iunit) nblk
do iblk=1,nblk
   associate(B => SBlock(iblk))
     write(iunit) iblk, B%n, B%l1, B%l2
     write(iunit) B%pos,B%matX,B%matY,B%vec
   end associate
enddo
associate(B => SBlockIV)
  write(iunit) B%n,B%l1,B%l2
  write(iunit) B%pos,B%vec
  write(iunit) B%matX,B%matY
end associate
close(iunit)
! deallocate blocks
do iblk=1,nblk
   associate(B => Sblock(iblk))
     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)
   end associate
enddo
! deallocate IV blocks
associate(B => SblockIV)
  deallocate(B%vec)
  deallocate(B%pos)
end associate

deallocate(Sblock)

end subroutine dump_Eblock

subroutine ModABMin_Act_FOFO(Occ,SRKer,Wt,OrbGrid,ABMin,&
                         MultpC,NSymNO,posAA,&
                         IndN,IndX,NDimX,NGrid,NBasis,&
                         NOccup,NAct,INActive,nAA,twokfile,twokerf)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO ABMIN
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NOccup,NAct,INActive,nAA
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),&
                      MultpC(15,15),NSymNO(NBasis),&
                      posAA(NAct*(NAct-1)/2)
character(*),intent(in) :: twokfile,twokerf
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: ABMin(nAA,nAA)

integer :: offset,batchlen,iunit1,iunit2
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
integer :: i1i2s,i3i4s,iSym
integer :: pos(NBasis,NBasis)
double precision :: XKer1234,TwoSR,Cpq,Crs
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

pos=0
do i=1,nAA
   j=posAA(i)
   pos(IndN(1,j),IndN(2,j)) = i
enddo

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(nAA,nAA))

ABKer = 0
do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,nAA
      ipq = posAA(IRow)
      ip = IndN(1,ipq)
      iq = IndN(2,ipq)

      do ICol=1,nAA
         irs = posAA(ICol)
         ir=IndN(1,irs)
         is=IndN(2,irs)

         if(irs.gt.ipq) cycle

         i1i2s = MultpC(NSymNO(ip),NSymNO(iq))
         i3i4s = MultpC(NSymNO(ir),NSymNO(is))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
            batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo

         ABKer(IRow,ICol) = ABKer(IRow,ICol) + XKer1234
         ABKer(ICol,IRow) = ABKer(IRow,ICol)

      enddo
   enddo

enddo

open(newunit=iunit1,file=trim(twokfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)
open(newunit=iunit2,file=trim(twokerf),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl  = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(k>INActive.and.(l>INActive.and.l<=NOccup).and.pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit1,rec=kl) work1(1:NBasis*NOccup)
        read(iunit2,rec=kl) work2(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints1(i,j) = work1((j-1)*NBasis+i)
              ints2(i,j) = work2((j-1)*NBasis+i)
           enddo
        enddo
        ints1(:,NOccup+1:NBasis) = 0
        ints2(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if((j>INActive.and.j<=NOccup).and.i>INActive.and.pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                if(irs.gt.ipq) cycle

                TwoSR = ints1(j,i)-ints2(j,i)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) &
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                ABMIN(irs,ipq) = ABMIN(ipq,irs)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit2)
close(iunit1)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_Act_FOFO

subroutine ModABMin_FOFO(Occ,SRKer,Wt,OrbGrid,ABMin,&
                         MultpC,NSymNO,&
                         IndN,IndX,NDimX,NGrid,NBasis,&
                         NAct,INActive,twokfile,twokerf,&
                         AB1,dfile)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
use timing
implicit none

integer,parameter  :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NAct,INActive
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX),&
                      MultpC(15,15),NSymNO(NBasis)
character(*),intent(in)     :: twokfile,twokerf
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
logical,intent(in)               :: AB1
double precision,intent(inout)   :: ABMin(NDimX,NDimX)
character(*),intent(in),optional :: dfile

integer :: offset,batchlen,iunit1,iunit2
integer :: iunit3,iunit4
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
integer :: i1i2s,i3i4s,iSym
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
integer :: NOccup,IGem(NBasis),AuxCoeff(3,3,3,3)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)
logical :: dump
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

dump = .false.
if(present(dfile)) dump=.true.
print*, dfile

NOccup = NAct + INActive
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

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = 1
if(AB1) then
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 0
               else
                  AuxCoeff(i,j,k,l) = 1
               endif
            enddo
         enddo
      enddo
   enddo
endif

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo


allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

ABKer = 0

!print*, 'ModABMin_FOFO'
do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ip = IndN(1,IRow)
      iq = IndN(2,IRow)
      ipq = IndX(IRow)

      do ICol=1,NDimX
         ir=IndN(1,ICol)
         is=IndN(2,ICol)
         irs=IndX(ICol)
         if(irs.gt.ipq) cycle

         i1i2s = MultpC(NSymNO(ip),NSymNO(iq))
         i3i4s = MultpC(NSymNO(ir),NSymNO(is))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
                       AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))* &
                       batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo

         ABKer(ipq,irs) = ABKer(ipq,irs) + XKer1234
         ABKer(irs,ipq) = ABKer(ipq,irs)

      enddo


   enddo

enddo

call clock('ModABMin:Ker',Tcpu,Twall)
call clock('START',Tcpu,Twall)

open(newunit=iunit1,file=trim(twokfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)
open(newunit=iunit2,file=trim(twokerf),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

if(dump) open(newunit=iunit3,file=trim(dfile),form='unformatted')

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit1,rec=kl) work1(1:NBasis*NOccup)
        read(iunit2,rec=kl) work2(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints1(i,j) = work1((j-1)*NBasis+i)
              ints2(i,j) = work2((j-1)*NBasis+i)
           enddo
        enddo
        ints1(:,NOccup+1:NBasis) = 0
        ints2(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                !if(irs.gt.ipq) cycle

                TwoSR = ints1(j,i)-ints2(j,i)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) &
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                !ABMIN(irs,ipq) = ABMIN(ipq,irs)
                if(dump) write(iunit3) 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))

              endif
           enddo
        enddo

      endif
   enddo
enddo

if(dump) close(iunit3)
close(iunit2)
close(iunit1)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

call clock('ModABMin',Tcpu,Twall)

end subroutine ModABMin_FOFO

subroutine ModABMinSym(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,&
           MultpC,NSymNO,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX), &
                      MultpC(15,15),NSymNO(NBasis)
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

double precision :: CICoef(NBasis)
double precision,allocatable :: work(:),batch(:,:),ABKer(:,:)
integer :: i,j,IRow,ICol,ia,ib,iab,ic,id,icd
integer :: i1i2s,i3i4s,iSym
integer :: offset,batchlen
double precision :: XKer1234,TwoSR,CA,CB,CC,CD
integer,external :: NAddr3,NAddrrK

allocate(work(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

ABKer = 0

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   work(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ia = IndN(1,IRow)
      ib = IndN(2,IRow)
      iab = IndX(IRow)

      do ICol=1,NDimX
         ic=IndN(1,ICol)
         id=IndN(2,ICol)
         icd=IndX(ICol)
         if(icd.gt.iab) cycle

         i1i2s = MultpC(NSymNO(ia),NSymNO(ib))
         i3i4s = MultpC(NSymNO(ic),NSymNO(id))
         iSym = MultpC(i1i2s,i3i4s)
         if(iSym/=1) cycle

         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + work(i)* &
            batch(i,ia)*batch(i,ib)*batch(i,ic)*batch(i,id)
         enddo

         ABKer(iab,icd) = ABKer(iab,icd) + XKer1234
         ABKer(icd,iab) = ABKer(iab,icd)

      enddo
   enddo

enddo

do IRow=1,NDimX
   ia = IndN(1,IRow)
   ib = IndN(2,IRow)
   iab = IndX(IRow)
   CA = CICoef(ia)
   CB = CICoef(ib)

   do ICol=1,NDimX
      ic=IndN(1,ICol)
      id=IndN(2,ICol)
      icd=IndX(ICol)
      CC=CICoef(ic)
      CD=CICoef(id)
      !if(icd.gt.iab) cycle

      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))

      ABMin(iab,icd) = ABMin(iab,icd) &
                       +4.0d0*(CA+CB)*(CD+CC)*(ABKer(iab,icd)+TwoSR)
      !ABMin(icd,iab) = ABMin(iab,icd)

   enddo
enddo

deallocate(ABKer,batch,work)

end subroutine ModABMinSym

subroutine EERPA_FOFO(ECorr,EVec,EVal,Occ,CICoef,IGem,   &
                      IAuxGem,IG1,IG2,IG3,IG4,IB, &
                      IndN,NOccup,NDimX,NBasis,IntKFile,IFrag)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IG1,IG2,IG3,IG4,IB
integer,intent(in) :: IGem(NBasis),IAuxGem(NBasis),IndN(2,NDimX)
integer,intent(in) :: NOccup,IFrag
character(*),intent(in) :: IntKFile
double precision,intent(out) :: ECorr
double precision,intent(in) :: Occ(NBasis),CICoef(NBasis)
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
integer :: OccProd(NBasis,NBasis),AuxIG(NBasis)
integer :: IFlPQRS,IFlP,IFlQ,IFlR,IFlS
integer :: NoVirtP,NoVirtQ,NoVirtR,NoVirtS,NoVirt
integer :: IBdy,IBdyG,ICond
double precision :: Cpq,Crs,SumY,Aux
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 1.d8
integer :: itmp

!do i=1,NBasis
!   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
!enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

OccProd = 0
do j=1,NBasis
   do i=1,NBasis
      if(Occ(i)*Occ(j)/=0d0) OccProd(i,j) = 1
   enddo
enddo

AuxIG = 0
do i=1,NBasis
   if(IGem(i)==IG1.or.IGem(i)==IG2.or.IGem(i)==IG3.or.IGem(i)==IG4) AuxIG(i) = 1
enddo


allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl = 0
itmp = 0
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

                NoVirt  = 0
                IFlPQRS = 0
                if(OccProd(ir,is)*OccProd(ip,iq)/=0) NoVirt = 1
                if(AuxIG(ir)*AuxIG(is)*AuxIG(ip)*AuxIG(iq)/=0) IFlPQRS = 1

                call IBody(IBdy,IGem(ip),IGem(iq),IGem(ir),IGem(is))

                ICond = 0
                if(IFlPQRS==1) then

                   if((IBdy.eq.IB)) ICond = 1

                   ! FOR GVB ONLY
                   if((IB.gt.2.and.NoVirt.Eq.1.and. &
                       IBdy.eq.IB-1)) ICond = 1

                   ! FOR GVB ONLY: IF IFrag=1,IB=2 (ONE-BODY), ALLOW ALL CASES
                   ! EXCEPT WHEN ALL ORBITALS ARE FROM THE SAME GEMINAL
                   if(IFrag.Eq.1.and.IB.Eq.2) then
                      call IBody(IBdyG,IAuxGem(ip),IAuxGem(iq),IAuxGem(ir),IAuxGem(is))
                      If(IBdyG.Ne.1) ICond = 1
                   endif

                endif !IFlPQRS

                if(ICond.Eq.1) then

                   ISkippedEig = 0
                   SumY = 0
                   do kk=1,NDimX
                      if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                         SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                      else
                         ISkippedEig = ISkippedEig + 1
                         Skipped(ISkippedEig) = EVal(kk)
                      endif
                   enddo

                   Aux = Crs*Cpq*SumY

                   if(iq.Eq.is.and.ip.Eq.ir) then
                      Aux = Aux - 0.5d0*(Occ(ip)*(1d0-Occ(is))+Occ(is)*(1d0-Occ(ip)))
                   endif

                   ECorr = ECorr + Aux*ints(j,i)

                endif ! ICond

              endif
           enddo
        enddo

      endif

   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

deallocate(Skipped)
deallocate(ints,work)

end subroutine EERPA_FOFO

subroutine EneGVB_FOFO(NAct,NElHlf,ETot,URe,Occ,C,XOne, &
                       IGem,IndN,NBasis,NInte1,IntKFile,NDimX,NGem)
implicit none

integer,intent(in) :: NAct,NElHlf,NBasis,NInte1,NDimX,NGem
character(*) :: IntKFile
double precision :: ETot

integer :: IndN(2,NDimX),IGem(NBasis)
double precision :: URe(NBasis,NBasis),Occ(NBasis),C(NBasis)
double precision :: XOne(NInte1)

integer :: i,j,ii,ij,ia,ib,iab
integer :: k,l,kl,kk,ll,klround
integer :: iunit
integer :: INActive,NOccup
integer,external :: NAddr3
double precision :: FacIJ,FacIK
double precision :: EOne, EIntraGem, EInterCoul, EInterExch

double precision :: HNO(NBasis,NBasis)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)

write(lout,'(/A)') ' Electronic GVB energy check:'

! set dimensions
INActive = NElHlf - NAct
NOccup = 2*NAct + INActive

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

EOne = 0d0
EIntraGem  = 0d0
EInterCoul = 0d0
EInterExch = 0d0

do i=1,NBasis
   EOne = EOne + 2d0*Occ(i)*HNO(i,i)
enddo

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(l>NOccup) cycle

      read(iunit,rec=kl) work1(1:NBasis*NOccup)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      ints(:,NOccup+1:NBasis) = 0

      if(IGem(k)==IGem(l)) then
         EIntraGem = EIntraGem + C(k)*C(l)*ints(k,l)
      else
         EInterExch = EInterExch - Occ(k)*Occ(l)*ints(k,l)
      endif

      FacIK = 4d0
      if(k==l) then
         do i=k,NOccup
            if(IGem(i)/=IGem(k)) then
            if(i==k) FacIK = 2d0
            EInterCoul = EInterCoul + FacIK*Occ(i)*Occ(k)*ints(i,i)
            endif
         enddo
     endif

   enddo
enddo

close(iunit)

write(LOUT,'(" One-electron energy",24X,F17.8)') EOne
write(LOUT,'(" GVB intra-gem electron interaction",9X,F17.8)') EIntraGem
write(LOUT,'(" GVB inter-gem Coulomb interaction",10X,F17.8)')  EInterCoul
write(LOUT,'(" GVB inter-gem exchange interaction",9X,F17.8)') EInterExch
write(LOUT,'(" Total GVB",33X,F18.8)') EOne+EIntraGem+EInterCoul+EInterExch

ETot = EOne + EIntraGem + EInterCoul + EInterExch

deallocate(ints,work2,work1)

end subroutine EneGVB_FOFO

subroutine check_mp2(EMP2,IntKFIle,INActive,NOccup,NBasis)
! check MP2 energy (only inactive-vritual enter)
! to compare with molpro use {mp2;core,0}
implicit none

integer,intent(in) :: INActive,NOccup,NBasis
character(*),intent(in) :: IntKFile
double precision,intent(inout) :: EMP2

integer :: i,j,k,l,a,b
integer :: iunit
double precision :: val,eps(NBasis,NBasis)
double precision,allocatable :: work(:),ints_bj(:,:),ints_bi(:,:)

allocate(work(NBasis**2),ints_bj(NBasis,NBasis),ints_bi(NBasis,NBasis))

eps = 0
open(10,file='fock.dat')
do i=NOccup+1,NBasis
   do j=1,INActive
      read(10,*) k,l,val
      eps(k,l) = val
   enddo
enddo
close(10)

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

EMP2 = 0
do j=1,INActive
   do b=NOccup+1,NBasis

      read(iunit,rec=(b+(j-1)*NBasis)) work(1:NBasis*NOccup)
      ! ints_bj
      do l=1,NOccup
         do k=1,NBasis
            ints_bj(k,l) = work((l-1)*NBasis+k)
         enddo
      enddo
      ints_bj(:,NOccup+1:NBasis) = 0

      do i=1,INActive
         read(iunit,rec=(b+(i-1)*NBasis)) work(1:NBasis*NOccup)
         ! ints_bi
         do l=1,NOccup
            do k=1,NBasis
               ints_bi(k,l) = work((l-1)*NBasis+k)
            enddo
         enddo
         ints_bi(:,NOccup+1:NBasis) = 0

         do a=NOccup+1,NBasis
            EMP2 = EMP2 - ints_bj(a,i)*(2d0*ints_bj(a,i)-ints_bi(a,j)) / &
                          (eps(a,i)+eps(b,j))
         enddo
      enddo

   enddo
enddo

write(LOUT,'(/,1x,a,13x,f16.8,/)') 'EMP2   Energy', EMP2

deallocate(ints_bi,ints_bj,work)

end subroutine check_mp2

subroutine JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                   INActive,NOccup,NDim,NDimX,NBasis,NInte1,IntJFile,IntKFile,ACAlpha,AB,ETot)

implicit none

integer,intent(in) :: NDim,NDimX,NBasis,NInte1
integer,intent(in) :: INActive,NOccup,AB
integer,intent(in) :: IGem(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision,intent(inout) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
!double precision,intent(inout) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(inout) :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision,intent(in) :: ACAlpha
double precision,intent(in) :: AuxCoeff(3,3,3,3),Occ(NBasis),RDM2val(NOccup,NOccup,NOccup,NOccup)
character(*) :: IntJFile,IntKFile
double precision,intent(inout),optional :: ETot

integer :: i,j,k,l,kk,ll,kl
integer :: ip,iq,ir,is,iu,ipq,irs
integer :: iunit1,iunit2
double precision :: val
double precision :: AuxVal,HNOCoef
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)

!print*, 'start JK loop:',AB

if(AB==1) then
   HNOCoef = -1
elseif(AB==0) then
   HNOCoef = 1 - ACAlpha
elseif(AB==2) then
   HNOCoef = 1
endif

allocate(work1(NBasis**2),ints(NBasis,NBasis))

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      !call triang_to_sq2(work1,ints,NBasis)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = HNOCoef*Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif

         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! exchange
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

      ! AUXILIARY MATRIX WMAT
      if(l<=NOccup) then
         do ir=1,NOccup
            val = 0
            val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                 sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                 sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                 sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                 sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
            WMAT(k,ir) = WMAT(k,ir) + val
         enddo
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then

            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k

         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo

      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      ! COMPUTE THE ENERGY FOR CHECKING
      if(present(ETot)) then
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*HNOCoef*Occ(k)
         if(k>INActive) then
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo

      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1

      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               endif
            enddo

         endif
      endif

      ! T1
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then

         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then

                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l

         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo
enddo

close(iunit2)

deallocate(ints)
deallocate(work1)

end subroutine JK_loop

end module
