module sapt_visual

use types

implicit none

contains

subroutine dump_visual(Flags,A,B,SAPT)
!
! dump visualization data:
! -- NO orbitals
! -- Qmat (E2disp)
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: iunit,i
integer :: NOccupA,NOccupB
character(:),allocatable :: dumpvis

dumpvis = 'SAPTVIS'

! get dimensions
NOccupA = A%num0 + A%num1
NOccupB = B%num0 + B%num1

open(newunit=iunit,file=dumpvis,form="UNFORMATTED")

write(iunit) A%NBasis,B%NBasis
write(iunit) NOccupA,NOccupB
write(iunit)(1.d0,i=1,NOccupA),(0.d0,i=NOccupA+1,A%NBasis),(1.d0,i=1,NOccupB),(0.d0,i=NOccupB+1,B%NBasis)
!write(iunit) A%Occ(1:A%NBasis),B%Occ(1:B%NBasis)
!write(iunit) A%CMO(1:A%NBasis,1:A%NBasis)
!write(iunit) B%CMO(1:B%NBasis,1:B%NBasis)
write(iunit) SAPT%ALOC(1:A%NBasis,1:A%NBasis)
write(iunit) SAPT%BLOC(1:B%NBasis,1:B%NBasis)
write(iunit) SAPT%Qmat(1:NOccupA,1:NOccupB)

close(iunit)
!print*, 'A%CMO: ',norm2(A%CMO)
!print*, 'B%CMO: ',norm2(B%CMO)
print*, 'A%CMO: ',norm2(SAPT%ALOC)
print*, 'B%CMO: ',norm2(SAPT%BLOC)
print*, 'Qmat:  ',norm2(SAPT%Qmat)

end subroutine dump_visual

subroutine dump_elect(Flags,A,B,SAPT)
!
! dump visualization data:
! -- QelA, QelB (E1elst)
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: iunit,i
integer :: NOccupA,NOccupB
character(:),allocatable :: dumpvis

dumpvis = 'SAPT_ELECT'

! get dimensions
NOccupA = A%num0 + A%num1
NOccupB = B%num0 + B%num1

open(newunit=iunit,file=dumpvis,form="UNFORMATTED")

write(iunit) A%NBasis,B%NBasis
write(iunit) NOccupA,NOccupB
write(iunit)(1.d0,i=1,NOccupA),(0.d0,i=NOccupA+1,A%NBasis),(1.d0,i=1,NOccupB),(0.d0,i=NOccupB+1,B%NBasis)
!write(iunit) A%Occ(1:A%NBasis),B%Occ(1:B%NBasis)
!write(iunit) A%CMO(1:A%NBasis,1:A%NBasis)
!write(iunit) B%CMO(1:B%NBasis,1:B%NBasis)
write(iunit) SAPT%ALOC(1:A%NBasis,1:A%NBasis)
write(iunit) SAPT%BLOC(1:B%NBasis,1:B%NBasis)
write(iunit) SAPT%Qmat(1:NOccupA,1:NOccupB)
write(iunit) SAPT%QelA(1:NOccupA)
write(iunit) SAPT%QelB(1:NOccupB)

close(iunit)
!print*, 'A%CMO: ',norm2(A%CMO)
!print*, 'B%CMO: ',norm2(B%CMO)
print*, 'A%CMO: ',norm2(SAPT%ALOC)
print*, 'B%CMO: ',norm2(SAPT%BLOC)
print*, 'Qmat:  ',norm2(SAPT%Qmat)

end subroutine dump_elect

subroutine test_saptvis()
implicit none

integer :: i,j
integer :: iunit
integer :: ANBasis,BNBasis
integer :: ANOccup,BNOccup
character(:),allocatable :: dumpvis
double precision :: e2dv
double precision,allocatable :: ACMO(:,:),BCMO(:,:),AOcc(:),BOcc(:)
double precision,allocatable :: QMat(:,:)

dumpvis = 'SAPTVIS'

open(newunit=iunit,file=dumpvis,form="UNFORMATTED")

read(iunit) ANBasis,BNBasis
read(iunit) ANOccup,BNOccup
allocate(ACMO(ANBasis,BNBasis),BCMO(ANBasis,BNBasis),&
         AOcc(ANBasis),BOcc(BNBasis),&
         QMat(ANOccup,BNOccup))
read(iunit) AOcc(1:ANBasis),BOcc(1:BNBasis)
read(iunit) ACMO(1:ANBasis,1:ANBasis)
read(iunit) BCMO(1:BNBasis,1:BNBasis)
read(iunit) Qmat(1:ANOccup,1:BNOccup)

close(iunit)

print*, 'A%CMO: ',norm2(ACMO)
print*, 'B%CMO: ',norm2(BCMO)
print*, 'Qmat:  ',norm2(Qmat)

e2dv = 0d0
do j=1,BNoccup
do i=1,ANoccup
 e2dv = e2dv + Qmat(i,j)
enddo
enddo

print*, 'sum Qmat = ',e2dv*1000.d0

deallocate(BOcc,AOcc)
deallocate(Qmat,ACMO,BCMO)

end subroutine test_saptvis

end module sapt_visual
