subroutine molprogrid0(NGrid,NBasis)
Implicit Real*8 (A-H,O-Z)

integer :: igrid
integer :: npt, ndiff, ntg
integer ::  NBasis, NGrid

 open(newunit=igrid,file='GRID',access='sequential',&
      form='unformatted',status='old')
! read(igrid)
 read(igrid) npt,ndiff,ntg
 close(igrid)
 if(NBasis.Ne.ntg) then
  stop  'Fatal Error in molprogrid: NBasis Ne ntg'
 endif
 NGrid=npt
end subroutine molprogrid0

subroutine molprogrid1(RR,NGrid)
Implicit Real*8 (A-H,O-Z)
double precision RR(3,NGrid)
double precision, allocatable :: wt(:)
integer :: i,j,k,igrid
integer :: npt, ndiff, ntg
integer ::  NGrid
 open(newunit=igrid,file='GRID',access='sequential',&
      form='unformatted',status='old')
! read(igrid)
 read(igrid) npt,ndiff,ntg
 allocate(wt(npt))
 call readgridmolpro(igrid,'GRIDKS  ',npt,RR,wt)
 close(igrid)
 deallocate(wt)
end subroutine molprogrid1

subroutine molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,UMOAO,NGrid,NBasis)
implicit real*8 (A-H,O-Z)

double precision OrbGrid(NGrid,NBasis),WGrid(NGrid),UMOAO(NBasis,NBasis)
double precision OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis),Aux(NGrid,NBasis)
integer ::  NBasis, NGrid
integer, allocatable :: mapinv(:)
double precision, allocatable :: r(:,:),wt(:),orbval(:,:,:)
integer :: i,j,k,igrid
integer :: npt, ndiff, ntg

 open(newunit=igrid,file='GRID',access='sequential',&
      form='unformatted',status='old')
 ! read(igrid)
 read(igrid) npt,ndiff,ntg
 If(NBasis.Ne.ntg) Stop 'Fatal Error in molprogrid: NBasis Ne ntg'
 If(NGrid.Ne.npt) Stop 'Fatal Error in molprogrid: NGrid Ne npt'

 allocate(r(3,npt),wt(npt),mapinv(ntg),orbval(npt,ndiff,ntg))
 !... idftgra: functional contains grad rho terms (1) , del.2 rho (2)
 !... ndiff = (idftgra + 1) * (idftgra + 2) * (idftgra + 3) / 6

 ! get grid points, integration weights and orbitals and its derivatives
 ! wt array contains weights of the grid points
 ! orbval array contains basis functions values and its derivatives at npt
 ! points of the grid
 ! e.g. orbval(k, 1, mapinv(i)) is chi_i(r_k) - i-th basis function in the k-th
 ! grid point value
 ! mapinv - orbital mapping array
 call readgridmolpro(igrid,'GRIDKS  ',npt,r,wt)
 Call CpyV(WGrid,wt,npt)
 call readorbsmolpro(igrid,'ORBVAL  ',mapinv,orbval,npt,ndiff,ntg)

 do J=1,NBasis
    do I=1,NGrid
       OrbGrid(I,J)=orbval(I,1,mapinv(J))
       if(ndiff.gt.1) then
         OrbXGrid(I,J)=orbval(I,2,mapinv(J))
         OrbYGrid(I,J)=orbval(I,3,mapinv(J))
         OrbZGrid(I,J)=orbval(I,4,mapinv(J))
       endif
    enddo
 enddo

! Call CpyV(Aux,OrbGrid,NBasis*NGrid)
! Call TrOrbG(OrbGrid,UMOAO,Aux,NGrid,NBasis)
! Call CpyV(Aux,OrbXGrid,NBasis*NGrid)
! Call TrOrbG(OrbXGrid,UMOAO,Aux,NGrid,NBasis)
! Call CpyV(Aux,OrbYGrid,NBasis*NGrid)
! Call TrOrbG(OrbYGrid,UMOAO,Aux,NGrid,NBasis)
! Call CpyV(Aux,OrbZGrid,NBasis*NGrid)
! Call TrOrbG(OrbZGrid,UMOAO,Aux,NGrid,NBasis)

 Call CpyV(Aux,OrbGrid,NBasis*NGrid)
 call dgemm('N','T',NGrid,NBasis,NBasis,1d0,Aux,NGrid,UMOAO,NBasis,0d0,OrbGrid,NGrid)
 Call CpyV(Aux,OrbXGrid,NBasis*NGrid)
 call dgemm('N','T',NGrid,NBasis,NBasis,1d0,Aux,NGrid,UMOAO,NBasis,0d0,OrbXGrid,NGrid)
 Call CpyV(Aux,OrbYGrid,NBasis*NGrid)
 call dgemm('N','T',NGrid,NBasis,NBasis,1d0,Aux,NGrid,UMOAO,NBasis,0d0,OrbYGrid,NGrid)
 Call CpyV(Aux,OrbZGrid,NBasis*NGrid)
 call dgemm('N','T',NGrid,NBasis,NBasis,1d0,Aux,NGrid,UMOAO,NBasis,0d0,OrbZGrid,NGrid)

 close(igrid)
 deallocate(orbval,wt,r,mapinv)

end subroutine molprogrid


subroutine readgridmolpro(iunit,text,npt,r,wt)
implicit none
integer :: iunit
integer :: ios
integer :: npt
double precision :: r(3,npt),wt(npt)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'
     stop
  endif
  if(label==text) then

     read(iunit) npt
     read(iunit) r(1:3,1:npt)
     read(iunit) wt(1:npt)
     exit
  endif
enddo

end subroutine readgridmolpro

subroutine readalphamolpro(alpha)
implicit none
double precision :: alpha,tmp
integer :: iunit
integer :: ios
character(8) :: text, label
logical :: check

text = 'CHIVAL  '

print*, 'here?'
inquire(file='GRID',exist=check)
print*, 'check',check
if (.not. check) then
   write(6,*) 'WARNING! Molpro GRID file not present!'
   write(6,'(/,1x,a,f12.5)') 'RS PARAM IN INPUT:',  alpha
   return
endif

open(newunit=iunit,file='GRID',access='sequential',&
     form='unformatted',status='old')
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'
     stop
  endif
  if(label==text) then
     read(iunit) tmp 
     exit
  endif
enddo

close(iunit)

if(tmp/=alpha) then
   write(6,'(/,1x,a,f12.5)') 'WARNING! RS PARAM IN INPUT:',  alpha
   write(6,'(1x,a,f12.5)')   '         RS PARAM IN MOLPRO:', tmp
!   write(6,'(1x,a,f12.5)') 'ASSUMING INPUT VALUE!',alpha
   alpha = tmp
endif

end subroutine readalphamolpro

subroutine readorbsmolpro(iunit,text,mapinv,orbval,npt,ndiff,ntg)
implicit none
integer :: iunit
integer :: ios
integer :: npt, ndiff, ntg
integer :: mapinv(ntg)
double precision :: orbval(npt,ndiff,ntg)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'
     stop
  endif
  if(label==text) then

     read(iunit) npt,ndiff,ntg
     read(iunit) orbval(1:npt,1:ndiff,1:ntg)
     read(iunit) mapinv(1:ntg)
     exit
  endif
enddo

end subroutine readorbsmolpro



