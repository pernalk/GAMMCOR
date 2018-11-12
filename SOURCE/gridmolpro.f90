subroutine molprogrid0(NGrid,NBasis)
Implicit Real*8 (A-H,O-Z)

integer :: i,j,k,igrid
integer :: npt, ndiff, ntg
integer ::  NBasis, NGrid
 
 open(newunit=igrid,file='GRID',access='sequential',&
      form='unformatted',status='old')
! read(igrid)
 read(igrid) npt,ndiff,ntg
 close(igrid)
 If(NBasis.Ne.ntg) Stop 'Fatal Error in molprogrid0: NBasis Ne ntg'
 NGrid=npt
end subroutine molprogrid0

!subroutine molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,NBasis)
subroutine molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,UMOAO,NGrid,NBasis)
Implicit Real*8 (A-H,O-Z)

double precision OrbGrid(NBasis,NGrid),WGrid(NGrid),UMOAO(NBasis,NBasis)
double precision OrbXGrid(NBasis,NGrid),OrbYGrid(NBasis,NGrid),OrbZGrid(NBasis,NGrid),Aux(NBasis,NGrid)
integer ::  NBasis, NGrid
double precision, allocatable :: r(:,:),wt(:),mapinv(:),orbval(:,:,:)
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
 Do I=1,NBasis
 Do J=1,NGrid
 OrbGrid(I,J)=orbval(J,1,mapinv(I))
 If(ndiff.Gt.1) Then
 OrbXGrid(I,J)=orbval(J,2,mapinv(I))
 OrbYGrid(I,J)=orbval(J,3,mapinv(I))
 OrbZGrid(I,J)=orbval(J,4,mapinv(I))
 EndIf
 EndDo
 EndDo

 Call CpyV(Aux,OrbGrid,NBasis*NGrid)
 Call TrOrbG(OrbGrid,UMOAO,Aux,NGrid,NBasis)
 Call CpyV(Aux,OrbXGrid,NBasis*NGrid)
 Call TrOrbG(OrbXGrid,UMOAO,Aux,NGrid,NBasis)
 Call CpyV(Aux,OrbYGrid,NBasis*NGrid)
 Call TrOrbG(OrbYGrid,UMOAO,Aux,NGrid,NBasis)
 Call CpyV(Aux,OrbZGrid,NBasis*NGrid)
 Call TrOrbG(OrbZGrid,UMOAO,Aux,NGrid,NBasis)

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

subroutine readorbsmolpro(iunit,text,mapinv,orbval,npt,ndiff,ntg)
implicit none
integer :: iunit
integer :: ios
integer :: npt, ndiff, ntg
double precision :: mapinv(ntg),orbval(npt,ndiff,ntg)
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



