module abuksfofo

use print_units
use grid_internal

implicit none

contains

subroutine AB_UKS_FOFO(ABPLUS,ABMin,noa,nva,nob,nvb,NDimX,NBasis, &
                       occa,occb,xfac, &
                       IntJaa,IntJbb,IntKaa,IntKbb,IntKab,&
                       ICholesky)
!
! compute A+B and A-B matrices fro UHF / UKS reference
!
! NDimX = nova + novb
! xfac : HF exchange fraction
!
! COULOMB INTEGRALS ARE READ FROM IntJFile IN (FF|OO) FORMAT
! EXCHANGE INTEGRALS ARE READ FROM IntKFile IN (FO|FO) FORMAT
!
implicit none

integer,intent(in) :: noa,nva,nob,nvb
integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: ICholesky
real(8),intent(in) :: xfac
real(8),intent(in) :: occa(NBasis),occb(NBasis)
character(*)       :: IntJaa,IntJbb
character(*)       :: IntKaa,IntKbb,IntKab
real(8),intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)

integer :: nvoa,nvob
integer :: i,j,a,ipos
real(8),allocatable :: blockP(:,:),blockM(:,:)

nvoa = nva*noa
nvob = nvb*nob

!print*, 'NDimX =', NDimX
!print*, 'nvoa  =', nvoa
!print*, 'nvob  =', nvob

ABPlus = 0d0
ABMin  = 0d0

if (ICholesky==1) then

   stop "Cholesky not ready in AB_UKS_FOFO"
elseif (ICholesky==0) then

   ! alpha-alpha
   allocate(blockP(nvoa,nvoa),blockM(nvoa,nvoa))
   blockP = 0d0
   ipos = 0
   do i=1,noa
      do a=1,nva
         ipos = ipos + 1
         blockP(ipos,ipos) = occa(noa+a)-occa(i)
      enddo
   enddo
   blockM = blockP

   !print*, 'diagonal part:'
   !print*, 'ABPLUS-aa =', norm2(blockP)
   !print*, 'ABMIN -aa =', norm2(blockM)

   !print*, 'ABMIN-diag'
   !do i=1,nvoa
   !  write(6,'(*(f13.8))') (blockM(i,j),j=1,nvoa)
   !enddo
   !write(LOUT,'()')

   ! 2-el part
   call JK_UKS_SameSpin(blockP,blockM,xfac,nva,noa,nvoa,NBasis,IntJaa,IntKaa)

   !print*, 'ABMIN-all'
   !do i=1,nvoa
   !  write(6,'(*(f13.8))') (blockM(i,j),j=1,nvoa)
   !enddo
   !write(LOUT,'()')

   print*, 'alpha-alpha'
   print*, 'ABPLUS-aa =', norm2(blockP)
   print*, 'ABMIN -aa =', norm2(blockM)

   ABPLUS(1:nvoa,1:nvoa) = blockP(1:nvoa,1:nvoa)
   ABMIN(1:nvoa,1:nvoa)  = blockM(1:nvoa,1:nvoa)

   deallocate(blockM,blockP)

   ! beta-beta
   allocate(blockP(nvob,nvob),blockM(nvob,nvob))
   blockP = 0d0
   ipos = 0
   do i=1,nob
      do a=1,nvb
         ipos = ipos + 1
         blockP(ipos,ipos) = occb(nob+a)-occb(i)
      enddo
   enddo
   blockM = blockP

   !print*,  ''
   !print*, 'diagonal part:'
   !print*, 'ABPLUS-bb =', norm2(blockP)
   !print*, 'ABMIN -bb =', norm2(blockM)

   call JK_UKS_SameSpin(blockP,blockM,xfac,nvb,nob,nvob,NBasis,IntJbb,IntKbb)

   !print*, 'beta-beta:'
   print*, 'ABPLUS-bb =', norm2(blockP)
   print*, 'ABMIN -bb =', norm2(blockM)

   ABPLUS(nvoa+1:nvoa+nvob,nvoa+1:nvoa+nvob) = blockP(1:nvob,1:nvob)
   ABMIN(nvoa+1:nvoa+nvob,nvoa+1:nvoa+nvob)  = blockM(1:nvob,1:nvob)

   deallocate(blockM,blockP)

   ! alpha-beta and beta-alfa
   allocate(blockP(nvoa,nvob))

   blockP = 0d0
   call JK_UKS_AlfaBeta(blockP,xfac,nva,noa,nvoa,nvb,nob,nvob,NBasis,IntKab)

   !print*,  ''
   print*, 'alfa-beta:'
   print*, 'ABPLUS-ab =', norm2(blockP)

   ABPLUS(1:nvoa,nvoa+1:nvoa+nvob) = blockP(1:nvoa,1:nvob)
   ABPLUS(nvoa+1:nvoa+nvob,1:nvoa) = transpose(blockP(1:nvoa,1:nvob))

   deallocate(blockP)

   print*, 'ABPLUS total =',norm2(ABPLUS)
   print*, 'ABMIN  total =',norm2(ABMIN)

endif ! ICholesky

end subroutine AB_UKS_FOFO

subroutine JK_UKS_SameSpin(ABPLUS,ABMIN,xfac,nv,no,nvo,NBasis,IntJFile,IntKFile)
!
! (vo | vo)
! two-electron part of the same-spin UKS hessian matrix :
! (A - B)_ai,bj =
!  - xfac * [ (ij | ab) - (ib | aj) ]
!  + (1-xfac) * [ Kxc... ] ! added later
!
! (A + B)_ai,bj = 2 (ia | jb )
!  - xfac * [ (ij | ab) + (ib | aj) ]
!  + (1-xfac) * [ kernel ] ! added later
!
implicit none

integer,intent(in)  :: nv,no,nvo
integer,intent(in)  :: NBasis
real(8),intent(in)  :: xfac
real(8),intent(out) :: ABPLUS(nvo,nvo),ABMIN(nvo,nvo)
character(*) :: IntJFile,IntKFile

integer :: iunit1,iunit2
integer :: i,j,k,l
integer :: ai,bj
integer :: irec

real(8),allocatable :: work1(:),work2(:)
real(8),allocatable :: ints(:,:)

allocate(work1(NBasis**2),ints(NBasis,NBasis))

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
        access='DIRECT',recl=8*NBasis*no)

! exchange loop (FO|FO)
irec = 0
do l=1,no
   irec = irec + no
   do k=1,nv
      irec = irec + 1
      read(iunit1,rec=irec) work1(1:NBasis*no)

      do j=1,no
         do i=1,nv
            ints(i,j) = work1((j-1)*NBasis+no+i)
         enddo
      enddo

      do j=1,no
         do i=1,nv
            ai = (j-1)*nv + i
            bj = (l-1)*nv + k
            ABPLUS(ai,bj) = ABPLUS(ai,bj) + 2d0*ints(i,j)
         enddo
      enddo

      do j=1,no
         do i=1,nv
            ai = (l-1)*nv + i
            bj = (j-1)*nv + k
            ABMIN(ai,bj)  = ABMIN(ai,bj)  + xfac*ints(i,j)
            ABPLUS(ai,bj) = ABPLUS(ai,bj) - xfac*ints(i,j)
         enddo
      enddo

   enddo
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
irec = 0
do l=1,no
   do k=1,no
      irec = irec + 1
      read(iunit2,rec=irec) work1(1:NBasis**2)
      do j=1,nv
         do i=1,nv
            ints(i,j) = work1((no+j-1)*NBasis+no+i)
         enddo
      enddo

      do j=1,nv
         do i=1,nv
            ai = (k-1)*nv + i
            bj = (l-1)*nv + j
            ABMIN(ai,bj)  = ABMIN(ai,bj)  - xfac*ints(i,j)
            ABPLUS(ai,bj) = ABPLUS(ai,bj) - xfac*ints(i,j)
         enddo
      enddo

   enddo
enddo

close(iunit2)

deallocate(ints,work1)

end subroutine JK_UKS_SameSpin

subroutine JK_UKS_AlfaBeta(ABPLUS,xfac,nva,noa,nvoa,nvb,nob,nvob,NBasis,IntKFile)
!
! (vo | vo)
! two-electron part of the same-spin UKS hessian matrix :
! (A + B)_ai,bj = 2 (ia | jb )
!  - xfac * [ (ij | ab) + (ib | aj) ]
!  + (1-xfac) * [ kernel ] ! added later
!
implicit none

integer,intent(in)  :: nva,noa,nvoa,nvb,nob,nvob
integer,intent(in)  :: NBasis
real(8),intent(in)  :: xfac
real(8),intent(out) :: ABPLUS(nvoa,nvob)
character(*) :: IntKFile

integer :: iunit
integer :: i,j,k,l
integer :: ai,bj
integer :: irec

real(8),allocatable :: work(:)
real(8),allocatable :: ints(:,:)

allocate(work(NBasis**2),ints(NBasis,NBasis))

open(newunit=iunit,file=trim(IntKFile),status='OLD', &
        access='DIRECT',recl=8*NBasis*noa)

! exchange loop (FO|FO)
irec = 0
do l=1,nob
   irec = irec + nob
   do k=1,nvb
      irec = irec + 1
      read(iunit,rec=irec) work(1:NBasis*noa)

      do j=1,noa
         do i=1,nva
            ints(i,j) = work((j-1)*NBasis+noa+i)
         enddo
      enddo

      do j=1,noa
         do i=1,nva
            ai = (j-1)*nva + i
            bj = (l-1)*nvb + k
            ABPLUS(ai,bj) = ABPLUS(ai,bj) + 2d0*ints(i,j)
         enddo
      enddo

   enddo
enddo

close(iunit)

deallocate(ints,work)

end subroutine JK_UKS_AlfaBeta

subroutine AB_UKS_KER(ABPLUS,Ca,Cb,Ena,Enb,IndNa,IndNb,xfac,noa,nva,nob,nvb,NDimX,NAO)
implicit none

integer,intent(in) :: noa,nva,nob,nvb
integer,intent(in) :: NDimX,NAO
integer, intent(in) :: IndNa(2,noa*nva),IndNb(2,nob*nvb)
real(8),intent(in) :: xfac
real(8),intent(in) :: Ca(NAO,NAO),Cb(NAO,NAO)
real(8),intent(in) :: Ena(NAO),Enb(NAO)

real(8),intent(inout) :: ABPLUS(NDimX,NDimX)

integer :: nvoa,nvob
integer :: NGrid

integer :: i,j
integer :: mapinv(NAO)
real(8) :: Pa(NAO,NAO),Pb(NAO,NAO)
real(8),allocatable :: blockS(:,:)
real(8),allocatable :: WGrid(:),XKer(:,:),XKerS(:)
real(8),allocatable :: RhoVeca(:),RhoVecb(:)
real(8),allocatable :: OrbGAO(:),OrbGa(:),OrbGb(:)
real(8),allocatable :: OrbXGa(:),OrbYGa(:),OrbZGa(:)
real(8),allocatable :: OrbXGb(:),OrbYGb(:),OrbZGb(:)
real(8) :: URe(NAO,NAO),Occa(NAO),Occb(NAO)
real(8) :: Rhoa,Rhob
real(8),allocatable :: work(:,:)

nvoa = nva*noa
nvob = nvb*nob

! proper way : loop over batches, OrbGrid(AO,batch)-->MO...

! dumb way to do it! load OrbGrid(AO)-->OrbGridMO
write(LOUT,'(/,1x,a)') "MOLPRO GRID"
call molprogrid0(NGrid,NAO)

write(LOUT,'(1x,a,i8)') "The number of Grid Points =",NGrid

allocate(OrbGAO(NGrid*NAO),OrbGa(NGrid*NAO),OrbGb(NGrid*NAO))
!allocate(OrbXGa(NGrid*NAO),OrbYGa(NGrid*NAO),OrbZGa(NGrid*NAO))
allocate(WGrid(NGrid),RhoVeca(NGrid),RhoVecB(NGrid))

call molprogridAO(OrbGAO,mapinv,WGrid,NGrid,NAO)
call dgemm('N','N',NGrid,NAO,NAO,1d0,OrbGAO,NGrid,Ca,NAO,0d0,OrbGa,NGrid)
call dgemm('N','N',NGrid,NAO,NAO,1d0,OrbGAO,NGrid,Cb,NAO,0d0,OrbGb,NGrid)

URe = 0d0
do i=1,NAO
   URe(i,i) = 1d0
enddo
Occa = 0d0; Occb = 0d0
Occa(1:noa)=0.5d0
Occb(1:nob)=0.5d0
RhoVecA=0d0 ; RhoVecB=0d0
do i=1,NGrid
   call DenGrid(i,Rhoa,Occa,URe,OrbGa,NGrid,NAO)
   call DenGrid(i,Rhob,Occb,URe,OrbGb,NGrid,NAO)
   RhoVeca(I)=Rhoa
   RhoVecb(I)=Rhob
enddo

block
real(8) :: rhoint
rhoint=0d0
do i=1,NGrid
   if (RhoVecA(i)>1.0e-10) then 
    rhoint = rhoint + WGrid(i)*RhoVecA(i)
   endif
enddo
print*, 'Rho Alpha =', rhoint
rhoint=0d0
do i=1,NGrid
   if (RhoVecb(i)>1.0e-10) then 
    rhoint = rhoint + WGrid(i)*RhoVecb(i)
   endif
enddo
print*, 'Rho Beta  =', rhoint
end block

allocate(XKer(NGrid,3))
call RhoKernelSpin(XKer,RhoVeca,RhoVecb,NGrid)

allocate(XKerS(NGrid))

!alpha-alpha
allocate(blockS(nvoa,nvoa))

XKerS = 2d0*XKer(:,1)
!blockS(1:nvoa,1:nvoa) = ABPLUS(1:nvoa,1:nvoa)
call AB_UKS_Spin(blockS,XKerS,WGrid,OrbGa,OrbGa,IndNa,IndNa,noa,nva,noa,nva,NGrid,NAO)
print*, 'AB-Ker aa =',norm2(blockS)
do j=1,nvoa
   do i=1,nvoa
      ABPLUS(i,j) = ABPLUS(i,j) + blockS(i,j)
   enddo
enddo
deallocate(blockS)

allocate(blockS(nvob,nvob))
XKerS = 2d0*XKer(:,3)
call AB_UKS_Spin(blockS,XKerS,WGrid,OrbGb,OrbGb,IndNb,IndNb,nob,nvb,nob,nvb,NGrid,NAO)
print*, 'AB-Ker bb =',norm2(blockS)
do j=1,nvob
   do i=1,nvob
      ABPLUS(nvoa+i,nvoa+j) = ABPLUS(nvoa+i,nvoa+j) + blockS(i,j)
   enddo
enddo
deallocate(blockS)

allocate(blockS(nvoa,nvob))
XKerS = 2d0*XKer(:,2)
call AB_UKS_Spin(blockS,XKerS,WGrid,OrbGa,OrbGb,IndNa,IndNb,noa,nva,nob,nvb,NGrid,NAO)
print*, 'AB-Ker ab =',norm2(blockS)
do j=1,nvob
   do i=1,nvoa
      ABPLUS(i,nvoa+j) = ABPLUS(i,nvoa+j) + blockS(i,j)
   enddo
enddo
allocate(work(nvob,nvoa))
work = 0d0
work = transpose(blockS)
do j=1,nvoa
   do i=1,nvob
      ABPLUS(nvoa+i,j) = ABPLUS(nvoa+i,j) + work(i,j)
   enddo
enddo

deallocate(work)
deallocate(blockS)

!! proba w AO ...? 
!!allocate(WGrid(NGrid),OrbGAO(NGrid*NAO))
!!call molprogridAO(OrbGAO,mapinv,WGrid,NGrid,NAO)
!!call get_den(NAO,Ca,Ena,1d0,Pa)
!!call get_den(NAO,Cb,Enb,1d0,Pb)
!!call DenAOGrid(RhoVeca,PA,OrbGAO,NGrid,NAO) ! like dft_rho_sparse
!!call DenAOGrid(RhoVecb,PB,OrbGAO,NGrid,NAO)

!allocate(OrbGAO(NGrid*NAO))
!call molprogrid(OrbGAO,mapinv,WGrid,NGrid,NAO)
!
!
deallocate(XKer)
!
deallocate(WGrid,OrbGAO,OrbGb,OrbGa)

end subroutine AB_UKS_KER

subroutine AB_UKS_Spin(ABPLUS,XKer,Wt,Ca,Cb,IndNa,IndNB,noa,nva,nob,nvb,NGrid,NAO)

implicit none

integer,intent(in) :: NGrid,NAO
integer,intent(in) :: noa,nva,nob,nvb
integer, intent(in) :: IndNa(2,noa*nva),IndNb(2,nob*nvb)
real(8),intent(in) :: XKer(NGrid),Wt(NGrid)
real(8),intent(in) :: Ca(NGrid,NAO),Cb(NGrid,NAO)
real(8),intent(inout) :: ABPLUS(nva*noa,nvb*nob) 

integer :: ai,bj,ig,i,a,j,b
integer :: nvoa,nvob
real(8) :: val

nvoa = noa*nva
nvob = nob*nvb

ABPLUS = 0d0
do ai=1,nvoa
   a = IndNa(1,ai)
   i = IndNa(2,ai)
   !print*, a,i
   do bj=1,nvob
      !if(bj.gt.ai) cycle
      b = IndNb(1,bj)
      j = IndNb(2,bj)
      val = 0d0
      do ig=1,NGrid
         val = val + Wt(ig)*Ca(ig,a)*Ca(ig,i)*Cb(ig,b)*Cb(ig,j)* &
               XKer(ig)
      enddo
      ABPLUS(ai,bj) = ABPLUS(ai,bj) + val
      !ABPLUS(ai,bj) = ABPLUS(bj,ai)
   enddo
enddo

!print*, 'AB Spin =', norm2(ABPLUS)

end subroutine AB_UKS_Spin

subroutine DenAOGrid(RhoVec,D,OrbGrid,NGrid,NAO)
implicit none

integer,intent(in)  :: NGrid,NAO
real(8),intent(in)  :: D(NAO,NAO),OrbGrid(NGrid,NAO)
real(8),intent(out) :: RhoVec(NGrid)

integer :: ig
real(8),allocatable :: Work(:,:)
double precision ,external :: ddot

real(8) :: vec1(NAO),vec2(NAO) 

RhoVec=0d0
allocate(Work(NGrid,NAO))
call dgemm('N','N',NGrid,NAO,NAO,1d0,OrbGrid,NGrid,D,NAO,0d0,Work,NGrid)
do ig=1,NGrid
   vec1=Work(ig,:)
   vec2=OrbGrid(ig,:)
   RhoVec(ig) = 2d0*ddot(NAO,vec2,1,vec1,1)
   !RhoVec(ig) = 2d0*ddot(NAO,Work(ig,:),1,OrbGrid(ig,:),1)
enddo
print*, 'RhoVec=',norm2(RhoVec)

deallocate(Work)

end subroutine DenAOGrid

end module
