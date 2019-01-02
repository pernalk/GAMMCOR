subroutine dump_grid(ChiRS)
implicit double precision (a-h, o-z)

include "common/tapes"
include "common/cbas"
include "common/cdft"
!integer :: idiff, ndiff, npt, nmap, ndfx, ndfc, ndf, ndf0
double precision :: ChiRS
double precision :: record, thresh
dimension :: icentres(ntg), map(ntg)
double precision, allocatable :: r(:,:),wt(:),orbval(:,:,:),mapinv(:)
double precision, parameter   ::  zero = 0d0, one = 1d0, tolorb = 1.d-12

write(iout,*) 'DUMPING DFT GRID TO A FILE'

! open recently used grid
record = zero; call grid_open(record)
call grid_get_accuracy(record, thresh); call grid_get_size(0, npt)
write(iout, '(a15, tr5, f8.1, /, a15, tr5, g8.2, /, a15, tr5, i8)') &
    & "grid record:", record, "grid accuracy:", thresh, "grid size:", npt
idiff = idftgra
ndiff = (idiff + 1) * (idiff + 2) * (idiff + 3) / 6

! calculate grid points, integration weights and orbitals and its derivatives
! wt array contains weights of the grid points
! orbval array contains basis functions values and its derivatives at npt points of the grid
! e.g. orbval(k, 1, mapinv(i)) is chi_i(r_k) - i-th basis function in the k-th grid point value
! mapinv, nmap - orbital mapping arrays and summation range
allocate(r(3, npt), wt(npt), orbval(npt, ndiff, ntg), mapinv(ntg))
call grid_obtain(0, npt, r, wt, -1, .true.)
call grid_orbital_initialize(tolorb)
call grid_orbital_value(npt, idiff, ndiff, r, tolorb, map, nmap, orbval, icentres)

forall(i = 1:nmap) mapinv(map(i)) = i
print*, sum(mapinv(1:ntg))

call find_free_unit(ifil)
open(unit=ifil,file='GRID',form='unformatted',status='unknown')
write(ifil) int(npt,kind=4),int(ndiff,kind=4),int(ntg,kind=4)
!grid
write(ifil) 'GRIDKS  '
write(ifil) int(npt,kind=4)
write(ifil) r(1:3,1:npt)
write(ifil) wt(1:npt)
! orbvals
write(ifil) 'ORBVAL  '
write(ifil) int(npt,kind=4), int(ndiff,kind=4), int(ntg,kind=4)
write(ifil) orbval(1:npt,1:ndiff,1:ntg)
write(ifil) mapinv(1:ntg)
write(ifil) 'CHIVAL  '
write(ifil) ChiRS
close(ifil) 

print*, 'wt',sum(wt(1:npt))
print*, 'mapinv',sum(mapinv(1:ntg))
!do i=1,npt
!   write(*,*) i,sum(orbval(i,:,:)) 
!enddo

call grid_orbital_term()
deallocate(r)
deallocate(mapinv,orbval,wt)

end subroutine dump_grid


