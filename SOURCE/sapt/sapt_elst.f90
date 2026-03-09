module sapt_elst
use types
use sapt_utils

implicit none

contains

subroutine e1elst(A,B,SAPT)
! calculates 1st order electrostatic energy
! in the AO basis
!
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i, j
integer :: NAO,NBas
double precision,allocatable :: PA(:,:),PB(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:),Ja(:,:)
double precision,allocatable :: work(:,:)
double precision :: tmp,ea,eb,eabel,elst
double precision,parameter :: Half=0.5d0
double precision,external  :: trace

! set dimensions
 NAO  = SAPT%NAO
 NBas = A%NBasis

 allocate(PA(NAO,NAO),PB(NAO,NAO),&
          Va(NAO,NAO),Vb(NAO,NAO),Ja(NAO,NAO))
 allocate(work(NAO,NAO))

 call get_den(NAO,NBas,A%CMO,A%Occ,2d0,PA)
 call get_den(NAO,NBas,B%CMO,B%Occ,2d0,PB)

 call get_one_mat('V',Va,A%Monomer,NAO)
 call get_one_mat('V',Vb,B%Monomer,NAO)

 call make_J1(NAO,PA,Ja,'AOTWOSORT')

! Tr[Pa.Va + Pb.Vb + Pb.Ja]
 work=0
 call dgemm('N','N',NAO,NAO,NAO,1d0,PA,NAO,Vb,NAO,0d0,work,NAO)
 ea = trace(work,NAO)
! print*, ea
 call dgemm('N','N',NAO,NAO,NAO,1d0,PB,NAO,Va,NAO,0d0,work,NAO)
 eb = trace(work,NAO)
! print*, eb
 call dgemm('N','N',NAO,NAO,NAO,1d0,PB,NAO,Ja,NAO,0d0,work,NAO)
 eabel = trace(work,NAO)
 elst = ea + eb + eabel + SAPT%Vnn

 call print_en('V_nucB_elA',ea,.false.)
 call print_en('V_nucA_elB',eb,.false.)
 call print_en('V_elA_elB',eabel,.false.)
 call print_en('V_nn',SAPT%Vnn,.false.)
 call print_en('Eelst',elst*1000,.true.)
 SAPT%elst = elst

 deallocate(work)
 deallocate(Ja,Vb,Va,PB,PA)

end subroutine e1elst

subroutine e1elst_o(A,B,SAPT)
!
! unrestricted open-shell electrostatic energy (AO)
! see Eq. (7) in https://doi.org/10.1063/1.4758455
!
implicit none

type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i, j
integer :: NBas

double precision :: ea,eb,elst
double precision,allocatable :: work(:,:)
double precision,allocatable :: Vb(:,:)
double precision,allocatable :: PA(:,:),PB(:,:)
double precision,external  :: trace

! set dimensions
NBas = A%NBasis

allocate(PA(NBas,NBas),PB(NBas,NBas),Vb(NBas,NBas))
allocate(work(NBas,NBas))

call get_den(NBas,NBas,A%UMO(:,:,1),A%UOcc(:,1),1d0,PA)
call get_den(NBas,NBas,A%UMO(:,:,2),A%UOcc(:,2),1d0,work)
PA = PA + work

call get_den(NBas,NBas,B%UMO(:,:,1),B%UOcc(:,1),1d0,PB)
call get_den(NBas,NBas,B%UMO(:,:,2),B%UOcc(:,2),1d0,work)
PB = PB + work

call get_one_mat('V',Vb,B%Monomer,NBas)

! Tr[PA.Vb + PB.WA] + Vnn
call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Vb,NBas,0d0,work,NBas)
ea = trace(work,NBas)

call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,A%WPot,NBas,0d0,work,NBas)
eb = trace(work,NBas)

elst = ea + eb + SAPT%Vnn

call print_en('PA.Vb',ea,.false.)
call print_en('PB.Wa',eb,.false.)
call print_en('V_nn',SAPT%Vnn,.false.)
call print_en('Eelst',elst*1000,.true.)

SAPT%elst = elst

deallocate(work)
deallocate(Vb,PB,PA)

end subroutine e1elst_o

end module sapt_elst
