module sapt_interface_chol

use types
use timing
use tran
use sorter
use tran_Chol
use gammcor_integrals, only : TCholeskyVecs, TCholeskyVecsOTF, TAOBasis, TSystem, &
                              chol_CoulombMatrix, chol_MOTransf_TwoStep, &
                              chol_gammcor_Rkab, CholeskyOTF_ao_vecs, &
                              ORBITAL_ORDERING_DALTON, ORBITAL_ORDERING_MOLPRO, &
                              ORBITAL_ORDERING_ORCA, ORBITAL_ORDERING_PYSCF, &
                              BECKE_PARAMS_MEDIUM, auto2e_init, &
                              auto2e_interface_C, &
                              sys_read_xyz, basis_NewAOBasis, &
                              Becke_MolecularGrid, gridfunc_orbitals, &
                              CholeskyOTF_Fock_MO_v2

use abmat
use read_external
use trexio

implicit none

contains

subroutine CholeskyOTF_elpot_AO(Mon,NBasis)
!
! obtain electrostatic potential in AO
! W = V + J
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis

integer      :: ione
logical      :: valid
character(8) :: label
character(:),allocatable :: onefile

double precision :: V(NBasis,NBasis)

if(Mon%Monomer==1) then
  onefile="ONEEL_A"
elseif(Mon%Monomer==2) then
  onefile="ONEEL_B"
endif

V = 0d0
open(newunit=ione,file=onefile,access='sequential',&
     form='unformatted',status='old')
read(ione)
read(ione) label,V
if(label=='POTENTAL') valid=.true.
close(ione)

if(.not.valid) then
   write(LOUT,'(1x,a)') 'V not found in calc_elpot_CholOTF!'
endif

allocate(Mon%WPot(NBasis,NBasis))

 Mon%WPot = V + Mon%Jmat

end subroutine CholeskyOTF_elpot_AO

subroutine chol_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
! transform Cholesky Vecs from AO to NO
! for all 2-index vecs needed in SAPT:
!   FFXX(NCholesky,NBas**2) -- for polarization
!   FFXY(NCholesky,NBas**2) -- for exchange
!   FFYX(NCholesky,NBas**2) -- for exchange
!   OOXX(NCholesky,dimO**2) -- for polarization
!   DCholX(NCholesky,NDimX) -- for polarization
!
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB,dimVA,dimVB,nOVA,nOVB
integer          :: i,j,ip,iq,ipq
double precision :: Cpq
double precision,allocatable :: tmp(:,:)
! test
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

NCholesky = CholeskyVecs%NCholesky
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1
dimVA = A%num1+A%num2
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

print*, 'dimOA',dimOA
print*, 'dimOB',dimOB

 !allocate(A%OV(NCholesky,A%NDimX),B%OV(NCholesky,B%NDimX))

 !allocate(tmp(NCholesky,nOVA))
 !! (OV|AA)
 !call chol_MOTransf(tmp,CholeskyVecs,&
 !                   A%CMO,1,dimOA,&
 !                   A%CMO,A%num0+1,NBasis)
 !A%OV = 0
 !do i=1,A%NDimX
 !   ip  = A%IndN(1,i)
 !   iq  = A%IndN(2,i)
 !   ipq = iq+(ip-A%num0-1)*dimOA
 !   A%OV(:,i)= tmp(:,ipq)
 !enddo
 !print*, 'A-OV',norm2(A%OV)
 !deallocate(tmp)

 !allocate(tmp(NCholesky,nOVB))
 !! (OV|BB)
 !call chol_MOTransf(tmp,CholeskyVecs,&
 !                   B%CMO,1,dimOB,&
 !                   B%CMO,B%num0+1,NBasis)

 !B%OV = 0
 !do i=1,B%NDimX
 !   ip  = B%IndN(1,i)
 !   iq  = B%IndN(2,i)
 !   ipq = iq+(ip-B%num0-1)*dimOB
 !   B%OV(:,i)= tmp(:,ipq)
 !enddo

 !deallocate(tmp)

 allocate(A%OO(NCholesky,dimOA**2),&
          B%OO(NCholesky,dimOB**2) )
 ! (OO|AA)
 !call chol_MOTransf(A%OO,CholeskyVecs,&
 !                  A%CMO,1,dimOA,&
 !                  A%CMO,1,dimOA)
 !                  B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%OO,CholeskyVecs,&
                    A%CMO,1,dimOA,&
                    A%CMO,1,dimOA,&
                    MaxBufferDimMB)
call clock('AOO',Tcpu,Twall)
 ! (OO|BB)
 !call chol_MOTransf(B%OO,CholeskyVecs,&
 !                   B%CMO,1,dimOB,&
 !                   B%CMO,1,dimOB)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%OO,CholeskyVecs,&
                    B%CMO,1,dimOB,&
                    B%CMO,1,dimOB,&
                    MaxBufferDimMB)
call clock('BOO',Tcpu,Twall)

print*, 'A%OO',norm2(A%OO)
print*, 'B%OO',norm2(B%OO)

! if(SAPT%SaptLevel==666) then ! RS2PT2+
    allocate(A%OOAB(NCholesky,dimOA*dimOB), &
             B%OOBA(NCholesky,dimOB*dimOA))

    call chol_MOTransf_TwoStep(A%OOAB,CholeskyVecs,&
                       A%CMO,1,dimOA,&
                       B%CMO,1,dimOB,&
                       MaxBufferDimMB)

    call chol_MOTransf_TwoStep(B%OOBA,CholeskyVecs,&
                       B%CMO,1,dimOB,&
                       A%CMO,1,dimOA,&
                       MaxBufferDimMB)
! endif
print*, 'A%OOAB',norm2(A%OOAB)
print*, 'B%OOBA',norm2(B%OOBA)

 allocate(A%FF(NCholesky,NBasis**2),&
          B%FF(NCholesky,NBasis**2) )
 ! (FF|AA)
 !call chol_MOTransf(A%FF,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   A%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%FF,CholeskyVecs,&
                    A%CMO,1,NBasis,&
                    A%CMO,1,NBasis,&
                    MaxBufferDimMB)
 call clock('AFF',Tcpu,Twall)
 ! (FF|BB)
 !call chol_MOTransf(B%FF,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%FF,CholeskyVecs,&
                    B%CMO,1,NBasis,&
                    B%CMO,1,NBasis,&
                    MaxBufferDimMB)
 call clock('BFF',Tcpu,Twall)

!print*, 'A%FF',norm2(A%FF)
!print*, A%FF(3,:)
!print*, 'B%FF',norm2(B%FF)

 ! DChol(NCholeksy,NDimX)
 allocate(A%DChol(NCholesky,A%NDimX), &
          B%DChol(NCholesky,B%NDimX))

 do j=1,A%NDimX
    ip = A%IndN(1,j)
    iq = A%IndN(2,j)
    ipq = iq + (ip-1)*NBasis
    Cpq = A%CICoef(ip) + A%CICoef(iq)
    A%DChol(:,j) = Cpq*A%FF(:,ipq)
 enddo

 do j=1,B%NDimX
    ip = B%IndN(1,j)
    iq = B%IndN(2,j)
    ipq = iq + (ip-1)*NBasis
    Cpq = B%CICoef(ip) + B%CICoef(iq)
    B%DChol(:,j) = Cpq*B%FF(:,ipq)
 enddo

 if(SAPT%SaptLevel==999) return

 allocate(A%FFAB(NCholesky,NBasis**2),&
          B%FFBA(NCholesky,NBasis**2) )
 ! (FF|AB)
 !call chol_MOTransf(A%FFAB,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   B%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%FFAB,CholeskyVecs,&
                    A%CMO,1,NBasis,&
                    B%CMO,1,NBasis,&
                    MaxBufferDimMB)
 ! (FF|BA)
 !call chol_MOTransf(B%FFBA,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   A%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%FFBA,CholeskyVecs,&
                    B%CMO,1,NBasis,&
                    A%CMO,1,NBasis,&
                    MaxBufferDimMB)

 print*, 'A%FFAB',norm2(A%FFAB)
 print*, 'B%FFBA',norm2(B%FFBA)

 !allocate(A%FO(NCholesky,NBasis*dimOA),&
 !         B%FO(NCholesky,NBasis*dimOA))
 !! (FO|AA)
 !call chol_MOTransf(A%FO,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   A%CMO,1,dimOA)
 !! (FO|BB)
 !call chol_MOTransf(B%FO,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   B%CMO,1,dimOB)

end subroutine chol_sapt_AO2NO_BIN

subroutine chol_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

if(SAPT%InterfaceType==1) then
  write(lout,*) 'Cholesky 3-index AO2NO transformation does not work with DALTON yet!'
  stop
endif

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

allocate(A%OO(NCholesky,dimOA**2),B%OO(NCholesky,dimOB**2))

call chol_gammcor_Rkab(A%OO,A%CAONO,1,dimOA,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('AOO',Tcpu,Twall)

call chol_gammcor_Rkab(B%OO,B%CAONO,1,dimOB,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('BOO',Tcpu,Twall)

print*, 'A%OO',norm2(A%OO)
print*, 'B%OO',norm2(B%OO)

!if(SAPT%SaptLevel==666) then ! RS2PT2+

   allocate(A%OOAB(NCholesky,dimOA*dimOB), &
            B%OOBA(NCholesky,dimOB*dimOA))

   call chol_gammcor_Rkab(A%OOAB,A%CAONO,1,dimOA,B%CAONO,1,dimOB, &
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING_MOLPRO)

   call chol_gammcor_Rkab(A%OOBA,B%CAONO,1,dimOB,A%CAONO,1,dimOA, &
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING_MOLPRO)

print*, 'A%OOAB',norm2(A%OOAB)
print*, 'B%OOBA',norm2(B%OOBA)
!endif

allocate(A%FF(NCholesky,NBasis**2),&
         B%FF(NCholesky,NBasis**2) )

call chol_gammcor_Rkab(A%FF,A%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('AFF',Tcpu,Twall)

call chol_gammcor_Rkab(B%FF,B%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('BFF',Tcpu,Twall)

print*, 'A%FF',norm2(A%FF)
print*, 'B%FF',norm2(B%FF)

! DChol(NCholeksy,NDimX)
allocate(A%DChol(NCholesky,A%NDimX), &
         B%DChol(NCholesky,B%NDimX))

do j=1,A%NDimX
   ip = A%IndN(1,j)
   iq = A%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = A%CICoef(ip) + A%CICoef(iq)
   A%DChol(:,j) = Cpq*A%FF(:,ipq)
enddo

do j=1,B%NDimX
   ip = B%IndN(1,j)
   iq = B%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = B%CICoef(ip) + B%CICoef(iq)
   B%DChol(:,j) = Cpq*B%FF(:,ipq)
enddo

if(SAPT%SaptLevel==999) return

allocate(A%FFAB(NCholesky,NBasis**2),&
         B%FFBA(NCholesky,NBasis**2) )

call chol_gammcor_Rkab(A%FFAB,A%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call chol_gammcor_Rkab(B%FFBA,B%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

print*, 'A%FFAB',norm2(A%FFAB)
print*, 'B%FFBA',norm2(B%FFBA)

end subroutine chol_sapt_AO2NO_OTF

subroutine chol_OO_sapt_AO2NO_BIN(A,B,CholeskyVecs,NBasis,MemVal,MemType)
implicit none

type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB

!   OOXX(NCholesky,dimO**2) -- for polarization

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

NCholesky = CholeskyVecs%NCholesky
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

allocate(A%OO(NCholesky,dimOA**2),&
         B%OO(NCholesky,dimOB**2) )

call chol_MOTransf_TwoStep(A%OO,CholeskyVecs,&
                    A%CMO,1,dimOA,&
                    A%CMO,1,dimOA,&
                    MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%OO,CholeskyVecs,&
                   B%CMO,1,dimOB,&
                   B%CMO,1,dimOB,&
                   MaxBufferDimMB)

allocate(A%OOAB(NCholesky,dimOA*dimOB), &
         B%OOBA(NCholesky,dimOB*dimOA))

call chol_MOTransf_TwoStep(A%OOAB,CholeskyVecs,&
                   A%CMO,1,dimOA,&
                   B%CMO,1,dimOB,&
                   MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%OOBA,CholeskyVecs,&
                       B%CMO,1,dimOB,&
                       A%CMO,1,dimOA,&
                       MaxBufferDimMB)

end subroutine chol_OO_sapt_AO2NO_BIN

subroutine chol_OO_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set orbital ordering
if(SAPT%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(SAPT%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

allocate(A%OO(NCholesky,dimOA**2),B%OO(NCholesky,dimOB**2))

call chol_gammcor_Rkab(A%OO,A%CAONO,1,dimOA,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING)

call clock('AOO',Tcpu,Twall)

call chol_gammcor_Rkab(B%OO,B%CAONO,1,dimOB,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING)

call clock('BOO',Tcpu,Twall)

allocate(A%OOAB(NCholesky,dimOA*dimOB), &
            B%OOBA(NCholesky,dimOB*dimOA))

call chol_gammcor_Rkab(A%OOAB,A%CAONO,1,dimOA,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

call clock('AOOAB',Tcpu,Twall)

call chol_gammcor_Rkab(B%OOBA,B%CAONO,1,dimOB,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

print*, 'A%OOAB',norm2(A%OOAB)
print*, 'A%OOBA',norm2(B%OOBA)

end subroutine chol_OO_sapt_AO2NO_OTF

subroutine chol_FO_sapt_AO2NO_BIN(A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
! prepare (NChol,NBasis*dimO) matrices for FOFO integrals (BIN version)
! (FO|AA), (FO|BB), (FO|AB), (FO|BA)
!
implicit none

type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB

! set dimensions
NCholesky = CholeskyVecs%NCholesky
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

allocate(A%FO(NCholesky,NBasis*dimOA),B%FO(NCholesky,NBasis*dimOB))

call chol_MOTransf_TwoStep(A%FO,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   A%CMO,1,dimOA, &
                   MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%FO,CholeskyVecs,&
                   B%CMO,1,NBasis,&
                   B%CMO,1,dimOB, &
                   MaxBufferDimMB)

allocate(A%FOAB(NCholesky,NBasis*dimOB),B%FOBA(NCholesky,NBasis*dimOA))

call chol_MOTransf_TwoStep(A%FOAB,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   B%CMO,1,dimOB, &
                   MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%FOBA,CholeskyVecs,&
                   B%CMO,1,NBasis,&
                   A%CMO,1,dimOA, &
                   MaxBufferDimMB)

end subroutine chol_FO_sapt_AO2NO_BIN

subroutine chol_FF_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
!   FFXX(NCholesky,NBas**2) -- for polarization
!   FFXY(NCholesky,NBas**2) -- for exchange
!
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB
integer          :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecs%NCholesky

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

allocate(A%FF(NCholesky,NBasis**2),&
        B%FF(NCholesky,NBasis**2) )
call chol_MOTransf_TwoStep(A%FF,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   A%CMO,1,NBasis,&
                   MaxBufferDimMB)
call clock('AFF',Tcpu,Twall)

call chol_MOTransf_TwoStep(B%FF,CholeskyVecs,&
                   B%CMO,1,NBasis,&
                   B%CMO,1,NBasis,&
                   MaxBufferDimMB)
call clock('BFF',Tcpu,Twall)

! DChol(NCholeksy,NDimX)
allocate(A%DChol(NCholesky,A%NDimX), &
         B%DChol(NCholesky,B%NDimX))

do j=1,A%NDimX
   ip = A%IndN(1,j)
   iq = A%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = A%CICoef(ip) + A%CICoef(iq)
   A%DChol(:,j) = Cpq*A%FF(:,ipq)
enddo

do j=1,B%NDimX
   ip = B%IndN(1,j)
   iq = B%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = B%CICoef(ip) + B%CICoef(iq)
   B%DChol(:,j) = Cpq*B%FF(:,ipq)
enddo

if(SAPT%SaptLevel==999) return

allocate(A%FFAB(NCholesky,NBasis**2))
!         B%FFBA(NCholesky,NBasis**2) )

call chol_MOTransf_TwoStep(A%FFAB,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   B%CMO,1,NBasis,&
                   MaxBufferDimMB)
!call chol_MOTransf_TwoStep(B%FFBA,CholeskyVecs,&
!                   B%CMO,1,NBasis,&
!                   A%CMO,1,NBasis,&
!                   MaxBufferDimMB)

print*, 'A%FFAB',norm2(A%FFAB)

end subroutine chol_FF_sapt_AO2NO_BIN

subroutine chol_FO_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
!
! prepare (NChol,NBasis*dimO) matrices for FOFO integrals (OTF version)
! (FO|AA), (FO|BB), (FO|AB), (FO|BA)
!
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set orbital ordering
if(SAPT%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(SAPT%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif

allocate(A%FO(NCholesky,NBasis*dimOA),B%FO(NCholesky,NBasis*dimOB))

call chol_gammcor_Rkab(A%FO,A%CAONO,1,NBasis,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,        &
                   AOBasis,ORBITAL_ORDERING)

call chol_gammcor_Rkab(B%FO,B%CAONO,1,NBasis,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,        &
                   AOBasis,ORBITAL_ORDERING)

call clock('AFO+BFO',Tcpu,Twall)

allocate(A%FOAB(NCholesky,NBasis*dimOB),B%FOBA(NCholesky,NBasis*dimOA))

call chol_gammcor_Rkab(A%FOAB,A%CAONO,1,NBasis,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,          &
                   AOBasis,ORBITAL_ORDERING)

call chol_gammcor_Rkab(B%FOBA,B%CAONO,1,NBasis,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,          &
                   AOBasis,ORBITAL_ORDERING)

call clock('AFOAB+BFOBA',Tcpu,Twall)

end subroutine chol_FO_sapt_AO2NO_OTF

subroutine chol_FOERF_AO2NO_OTF(SAPT,MON,CholErfVecsOTF,AOBasis,Flags,NBasis)
!
! prepare (NCholERF,NBasis*dimO) matrices for ERF FOFO integrals (OTF version)
! (FO | erf(w_A r) | AA)
!
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: MON
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholErfVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholErf
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: dimO
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholErf = CholErfVecsOTF%Chol2Data%NVecs
dimO = MON%num0+MON%num1

! set orbital ordering
if(SAPT%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(SAPT%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   print*, SAPT%InterfaceType
   stop
endif

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif

allocate(MON%FOErf(NCholErf,NBasis*dimO))

call chol_gammcor_Rkab(MON%FOErf,MON%CAONO,1,NBasis,MON%CAONO,1,dimO, &
                   MaxBufferDimMB,CholErfVecsOTF,AOBasis,ORBITAL_ORDERING)

if (MON%Monomer==1) call clock('AFOERF',Tcpu,Twall)
if (MON%Monomer==2) call clock('BFOERF',Tcpu,Twall)

end subroutine chol_FOERF_AO2NO_OTF

subroutine chol_FFERF_AO2NO_OTF(Flags,M,CholErfVecsOTF,AOBasis,NBasis)
!
! performs AO2NO transformation
! to generate FFERF(NCholErf,NBasis**2) vectors
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: M
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholErfVecsOTF
integer,intent(in)     :: NBasis

integer :: NCholErf
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholErf = CholErfVecsOTF%Chol2Data%NVecs

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

! set orbital ordering
if(Flags%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(Flags%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

allocate(M%FFERF(NCholErf,NBasis**2))

call chol_gammcor_Rkab(M%FFErf,M%CAONO,1,NBasis,M%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholErfVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

if (M%Monomer==1) call clock('AFFErf',Tcpu,Twall)
if (M%Monomer==2) call clock('BFFErf',Tcpu,Twall)

end subroutine chol_FFERF_AO2NO_OTF

subroutine chol_FFXX_mon_AO2NO_OTF(Flags,M,CholeskyVecsOTF,AOBasis,NBasis)
!
! performs AO2NO transformation
! to generate 1) FFXX(NCholesky,NBasis**2) vectors
!             2) (c_p+c_q)*FFXX => DChol vectors
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: M
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

! set orbital ordering
if(Flags%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(Flags%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

allocate(M%FF(NCholesky,NBasis**2))

call chol_gammcor_Rkab(M%FF,M%CAONO,1,NBasis,M%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

if (M%Monomer==1) call clock('AFF',Tcpu,Twall)
if (M%Monomer==2) call clock('BFF',Tcpu,Twall)

allocate(M%DChol(NCholesky,M%NDimX))

do j=1,M%NDimX
   ip = M%IndN(1,j)
   iq = M%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = M%CICoef(ip) + M%CICoef(iq)
   M%DChol(:,j) = Cpq*M%FF(:,ipq)
enddo

end subroutine chol_FFXX_mon_AO2NO_OTF

subroutine chol_FFXY_AB_AO2NO_OTF(Flags,A,B,CholeskyVecsOTF,AOBasis,NBasis,abtype)
!
! performs AO2NO transformation
! to generate 1) FFXY(NCholesky,NBasis**2) vectors
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: A,B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
integer,intent(in)     :: NBasis
character(2),intent(in)   :: abtype

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING

integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%Chol2Data%NVecs

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

! set orbital ordering
if(Flags%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(Flags%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

if (abtype == "AB") then

   allocate(A%FFAB(NCholesky,NBasis**2))
   call chol_gammcor_Rkab(A%FFAB,A%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING)

elseif (abtype == "BA") then

   allocate(B%FFBA(NCholesky,NBasis**2))
   call chol_gammcor_Rkab(B%FFBA,B%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING)

endif

end subroutine chol_FFXY_AB_AO2NO_OTF

subroutine chol_JKmat_AO_OTF(Mon,NBasis)
!
! generates J and K matrices in MO
! and backtransforms them to AO
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis

integer :: i,j
integer :: ione
integer :: NOccup,NCholesky
double precision :: val
double precision :: D_no(NBasis,NBasis)
double precision :: SC(NBasis,NBasis),SAO(NBasis,NBasis)
double precision,allocatable :: ints(:),work(:,:)
double precision,allocatable :: Jtmp(:,:),Ktmp(:,:)
character(8) :: label
double precision,external :: ddot

NOccup = Mon%num0+Mon%num1
NCholesky = Mon%NChol

! prepare NO 1-density
D_no = 0d0
do i=1,NOccup
   D_no(i,i) = Mon%Occ(i)
enddo

! prepare Jmat in AO
allocate(Mon%Jmat(NBasis,NBasis))
if(.not.allocated(Mon%Kmat)) allocate(Mon%Kmat(NBasis,NBasis))

allocate(Jtmp(NBasis,NBasis),Ktmp(NBasis,NBasis))
allocate(ints(NBasis**2),work(NBasis,NBasis))

ints = 0d0
Jtmp = 0d0
Ktmp = 0d0
do i=1,NCholesky
   ints(:) = Mon%FF(i,:)
   val = ddot(NBasis**2,ints,1,D_no,1)
   call daxpy(NBasis**2,2d0*val,ints,1,Jtmp,1)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,ints,NBasis, &
              D_no,NBasis,0d0,work,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,-1d0,work,NBasis, &
              ints,NBasis,1d0,Ktmp,NBasis)
enddo

!print*, 'Jtmp', norm2(Jtmp)

! backtransform J to AO
! get S in AO
open(newunit=ione,file='ONEEL_A',access='sequential',&
     form='unformatted',status='old')
read(ione) label, SAO
close(ione)

! J_AO = SC . J_NO . (SC)^T
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SAO,NBasis,Mon%CAONO,NBasis,0d0,SC,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,Jtmp,NBasis,0d0,work,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work,NBasis,SC,NBasis,0d0,mon%Jmat,NBasis)
!print*, 'Jmat-AO', norm2(Mon%Jmat)

! K_AO = SC . K_NO . (SC)^T
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,Ktmp,NBasis,0d0,work,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,-1d0,work,NBasis,SC,NBasis,0d0,mon%Kmat,NBasis)

deallocate(Ktmp,Jtmp)
deallocate(work,ints)

end subroutine chol_JKmat_AO_OTF

end module sapt_interface_chol
