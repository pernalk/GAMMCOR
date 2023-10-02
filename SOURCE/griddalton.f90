subroutine daltongrid0(NGrid,gridfile,DoGGA,NBasis)
implicit none

integer,intent(in)  :: NBasis
character(*),intent(in)  :: gridfile
integer,intent(out) :: NGrid
logical,intent(out) :: DoGGA

integer :: i,j,k,igrid
integer :: npt, ndiff, ntg
character(8) :: label

open(newunit=igrid,file=gridfile,access='sequential',&
     form='unformatted',status='old')
 read(igrid) label
 read(igrid) NGrid
 read(igrid) DoGGA

 close(igrid)
 if(label /= 'DFTGRID ') then
    write(6,'(1x,2a)') 'Wrong label in daltongrid0',label
    stop
 endif

end subroutine daltongrid0

subroutine daltongrid_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,NGrid,gridfile,NBasis)
implicit none

integer,intent(in)           :: NBasis, NGrid
character(*),intent(in)      :: gridfile
double precision,intent(out) :: OrbGrid(NGrid,NBasis),WGrid(NGrid)
double precision,intent(out) :: OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)

integer :: igrid
integer :: i
double precision :: COR(3)
!double precision :: Aux(NGrid,NBasis)

open(newunit=igrid,file=gridfile,access='sequential',&
     form='unformatted',status='old')

rewind(igrid)
read(igrid)
read(igrid)
read(igrid)
do i=1,NGrid
   read(igrid) WGrid(i),COR,OrbGrid(i,1:NBasis), &
               OrbXGrid(i,1:NBasis),OrbYGrid(i,1:NBasis),OrbZGrid(i,1:NBasis)
!   write(6,'(1x,i4,3(a,f12.6,1x))') i,'x = ', COR(1), 'y = ', COR(2), 'z = ',COR(3)
enddo
close(igrid)

end subroutine daltongrid_gga

subroutine daltongrid_lda(OrbGrid,WGrid,NGrid,gridfile,NBasis)
implicit none

integer,intent(in)           :: NBasis, NGrid
character(*),intent(in)      :: gridfile
double precision,intent(out) :: OrbGrid(NGrid,NBasis),WGrid(NGrid)

integer :: i,igrid
double precision :: COR(3)

open(newunit=igrid,file=gridfile,access='sequential',&
     form='unformatted',status='old')
rewind(igrid)
read(igrid)
read(igrid)
read(igrid)
do i=1,NGrid
   read(igrid) WGrid(i),COR,OrbGrid(i,1:NBasis)
!   print*, 'i = ', i, 'out of ', NGrid
!   write(6,'(1x,i4,3(a,f12.6,1x))') i,'x = ', COR(1), 'y = ', COR(2), 'z = ',COR(3)
!   write(6,*) WGrid(i),norm2(OrbGrid(i,1:NBasis)) 
enddo
close(igrid)

print*, 'DaltonGrid-LDA'
print*, 'WGrid', norm2(WGrid)
print*, 'OrbGrid', norm2(OrbGrid)

end subroutine daltongrid_lda

subroutine daltongrid_tran_lda(OrbGrid,CAONO,UCen,NGrid,NBasis,switchAB)
!
! UCen - no of symmetry-independent centers is needed only 
!        for swapping of orbitals in grid
!
implicit none

integer,intent(in)             :: UCen
integer,intent(in)             :: NBasis,NGrid
logical,intent(in)             :: switchAB
double precision               :: CAONO(NBasis,NBasis)
!double precision,intent(in)    :: CAONO(NBasis,NBasis)
double precision,intent(inout) :: OrbGrid(NGrid,NBasis)

double precision :: Aux(NGrid,NBasis)

if (switchAB) then
   call CpyV(Aux,OrbGrid,NBasis*NGrid)
   call swap_orbgrid(Aux,UCen,NGrid,NBasis)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbGrid,NGrid)
else
   call CpyV(Aux,OrbGrid,NBasis*NGrid)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbGrid,NGrid)
endif

!print*, 'OrbGrid-daltongrid',norm2(OrbGrid)
!print*, OrbGrid(1,1:4)
!print*, OrbGrid(10,5:9)

end subroutine daltongrid_tran_lda

subroutine daltongrid_tran_gga(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,CAONO,UCen,NGrid,NBasis,switchAB)
!
! UCen - no of symmetry-independent centers is needed only 
!        for swapping of orbitals in grid
! monomer - swapping is done only for monomer B (monomer=2)
!
implicit none

integer,intent(in)             :: UCen
integer,intent(in)             :: NBasis,NGrid
logical,intent(in)             :: switchAB
double precision,intent(in)    :: CAONO(NBasis,NBasis)
double precision,intent(inout) :: OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: OrbXGrid(NGrid,NBasis),OrbYGrid(NGrid,NBasis),OrbZGrid(NGrid,NBasis)

double precision :: Aux(NGrid,NBasis)

if (switchAB) then
   call CpyV(Aux,OrbGrid,NBasis*NGrid)
   call swap_orbgrid(Aux,UCen,NGrid,NBasis)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbGrid,NGrid)
   call CpyV(Aux,OrbXGrid,NBasis*NGrid)
   call swap_orbgrid(Aux,UCen,NGrid,NBasis)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbXGrid,NGrid)
   call CpyV(Aux,OrbYGrid,NBasis*NGrid)
   call swap_orbgrid(Aux,UCen,NGrid,NBasis)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbYGrid,NGrid)
   call CpyV(Aux,OrbZGrid,NBasis*NGrid)
   call swap_orbgrid(Aux,UCen,NGrid,NBasis)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbZGrid,NGrid)
else
   call CpyV(Aux,OrbGrid,NBasis*NGrid)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbGrid,NGrid)
   call CpyV(Aux,OrbXGrid,NBasis*NGrid)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbXGrid,NGrid)
   call CpyV(Aux,OrbYGrid,NBasis*NGrid)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbYGrid,NGrid)
   call CpyV(Aux,OrbZGrid,NBasis*NGrid)
   call dgemm('N','N',NGrid,NBasis,NBasis,1d0,Aux,NGrid,CAONO,NBasis,0d0,OrbZGrid,NGrid)
endif

!print*, 'OrbGrid -daltongrid',norm2(OrbGrid)
!print*, 'OrbXGrid-daltongrid',norm2(OrbXGrid)
!print*, 'OrbYGrid-daltongrid',norm2(OrbYGrid)
!print*, 'OrbZGrid-daltongrid',norm2(OrbZGrid)

end subroutine daltongrid_tran_gga

subroutine swap_orbgrid(OrbGrid,BUCen,NGrid,NBasis)
!
! swap columns in OrbGrid(NGrid,NBasis)
!
use read_external
!
implicit none

integer,intent(in) :: NGrid,NBasis
integer,intent(in) :: BUCen
double precision,intent(inout) :: OrbGrid(NGrid,NBasis)

integer :: iA,iB,iAB,offset
integer :: iunit
logical :: ex
double precision,allocatable :: work(:,:)

integer :: ANSym,ANOrb,ANSymOrb(8),ANMonBas(8)
integer :: BNSym,BNOrb,BNSymOrb(8),BNMonBas(8)

!print*,'swap_orbgrid...'

! get orbital information for both monomers
ANSymOrb = 0; BNSymOrb = 0
call read_orbinf_dalton('SIRIFC_A',ANSym,ANOrb,ANSymOrb)
call read_orbinf_dalton('SIRIFC_B',BNSym,BNOrb,BNSymOrb)

!print*, 'NSym A, NSym B: ',ANSym, BNSym
!print*, 'NSymOrb A : ',ANSymOrb
!print*, 'NSymOrb B : ',BNSymOrb

ANMonBas = 0; BNMonBas = 0
call read_syminf_dalton(ANSym,BNSym,BUCen,ANSymOrb,BNSymOrb,&
                        ANMonBas,BNMonBas)

!print*, 'NMonBas A : ',ANMonBas 
!print*, 'NMonBas B : ',BNMonBas

! swap cols
call gen_swap_cols(OrbGrid,NGrid,NBasis,BNSym,ANMonBas,BNMonBas)

end subroutine swap_orbgrid

