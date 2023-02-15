module sapt_Chol_exch

use types
use read_external
use sapt_utils

implicit none

contains

subroutine e1exch_Chol(Flags,A,B,SAPT)
!
! E1exch(S2): Eq (9) in SAPT(MC) paper
! doi: 10.1021/acs.jctc.1c00344
! adapted e1exch_NaNb from sapt_pol.f90
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: i, j, k, l, ia, jb
integer :: ij,ipr
integer :: ip,iq,ir,is
integer :: ipq,iu,it
integer :: iunit
integer :: rdm2type
integer :: dimOA,dimOB
integer :: NAO,NBas
double precision :: fac,val,nnS2,tmp
double precision :: tElst,tvk(3),tNa(2),tNb(2),tNaNb
double precision :: exchs2
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:),RDM2Bval(:,:,:,:)
double precision,allocatable :: intA(:,:,:,:),intB(:,:,:,:)
double precision,allocatable :: tmpAB(:,:,:,:)
double precision,allocatable :: work(:,:),ints(:)
double precision,external  :: ddot

print*, 'Testing E1exch Chol...'

! set dimensions
NAO  = SAPT%NAO
NBas = A%NBasis
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

allocate(S(NBas,NBas),Sab(NBas,NBas))
allocate(Va(NBas,NBas),Vb(NBas,NBas),&
         Vabb(NBas,NBas),Vbaa(NBas,NBas),&
         Vaab(NBas,NBas),Vbba(NBas,NBas))

call get_one_mat('V',Va,A%Monomer,NBas)
call get_one_mat('V',Vb,B%Monomer,NBas)

call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)
call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
call tran2MO(Vb,B%CMO,A%CMO,Vbba,NBas)

call get_one_mat('S',S,A%Monomer,NBas)
call tran2MO(S,A%CMO,B%CMO,Sab,NBas)

allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
         RDM2Bval(dimOB,dimOB,dimOB,dimOB))
allocate(intA(dimOA,dimOA,dimOA,dimOB),&
         intB(dimOB,dimOB,dimOA,dimOB))

if(Flags%ICASSCF==1) then
   ! CAS
   RDM2Aval = A%RDM2val
   RDM2Bval = B%RDM2val
elseif(Flags%ICASSCF==0) then
   ! GVB
   RDM2Aval = A%RDM2val
   RDM2Bval = B%RDM2val
endif

call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,RDM2Aval,dimOA**3,Sab,NBas,0d0,intA,dimOA**3)
! careful! intB(B,B,A,B)
do is=1,dimOB
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,RDM2Bval(:,:,:,is),dimOB**2,Sab,NBas,0d0,intB(:,:,:,is),dimOB**2)
enddo

deallocate(RDM2Bval,RDM2Aval)
deallocate(Vb,Va,S)

! n^A n^B Sab Sab
nnS2 = 0
do j=1,dimOB
do i=1,dimOA
   nnS2 = nnS2 + A%Occ(i)*B%Occ(j)*Sab(i,j)**2
enddo
enddo

tElst = 2d0*(SAPT%elst-SAPT%Vnn)*nnS2
print*, 'tELST',tELST*1000

allocate(ints(NBas**2),work(NBas,NBas))

tvk = 0
do iq=1,dimOB
   do ip=1,dimOA
      tvk(1) = tvk(1) + A%Occ(ip)*B%Occ(iq)*Vaab(ip,iq)*Sab(ip,iq)
   enddo
enddo
tvk(1) = -2d0*tvk(1)
print*, 'tvk(1)',tvk(1)*1000

do iq=1,dimOB
   do ip=1,dimOA
      tvk(2) = tvk(2) + A%Occ(ip)*B%Occ(iq)*Vbba(iq,ip)*Sab(ip,iq)
   enddo
enddo
tvk(2) = -2d0*tvk(2)
print*, 'tvk(2)',tvk(2)*1000

! work = PA
call get_den(NAO,NBas,A%CMO,A%Occ,1d0,work)
! tvk = n_p n_q v_pq^qp = PA.K[PB]
do jb=1,NBas
   do ia=1,NBas
      tvk(3) = tvk(3) + work(ia,jb)*B%Kmat(jb,ia)
   enddo
enddo
tvk(3) = -2.0d0*tvk(3)

! tNa
tNa = 0
do it=1,dimOB
   do iq=1,dimOA
      do ir=1,dimOA
         do ip=1,dimOA
            tNa(1) = tNa(1) + B%Occ(it)*intA(ip,ir,iq,it)*Sab(iq,it)*Vbaa(ip,ir)
         enddo
      enddo
   enddo
enddo
tNa(1) = -2d0*tNa(1)

!(OO|OO): (AA|AB)
open(newunit=iunit,file='OOOOAAAB',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOA)

! one loop over integrals
ints = 0
ij = 0
do it=1,dimOB
   do iq=1,dimOA
      ij = ij + 1
      read(iunit,rec=iq+(it-1)*dimOA) ints(1:dimOA*dimOA)

      do ir=1,dimOA
         do ip=1,dimOA
            tNa(2) = tNa(2) + B%Occ(it)*intA(ip,ir,iq,it)*ints(ip+(ir-1)*dimOA)
         enddo
      enddo

   enddo
enddo
tNa(2) = -2d0*tNa(2)
print*, 'tNa-1',tNa(1)*1000
print*, 'tNa-2',tNa(2)*1000

close(iunit)

! tNb
tNb = 0
do iq=1,dimOB
   do it=1,dimOA
      do ir=1,dimOB
         do ip=1,dimOB
            tNb(1) = tNb(1) + A%Occ(it)*intB(ip,ir,it,iq)*Sab(it,iq)*Vabb(ip,ir)
         enddo
      enddo
   enddo
enddo
tNb(1) = -2d0*tNb(1)
!
!(OO|OO): (BB|BA)
open(newunit=iunit,file='OOOOBBBA',status='OLD', &
     access='DIRECT',recl=8*dimOB*dimOB)

! one loop over integrals
ints = 0
do it=1,dimOA
   do iq=1,dimOB
      read(iunit,rec=iq+(it-1)*dimOB) ints(1:dimOB*dimOB)

      do ir=1,dimOB
         do ip=1,dimOB
            tNb(2) = tNb(2) + A%Occ(it)*intB(ip,ir,it,iq)*ints(ip+(ir-1)*dimOB)
         enddo
      enddo

   enddo
enddo
tNb(2) = -2d0*tNb(2)
print*, 'tNb-1',tNb(1)*1000
print*, 'tNb-2',tNb(2)*1000

close(iunit)

open(newunit=iunit,file='OOOOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA*dimOA)

allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))

call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)

val  = 0
work = 0

ints = 0
do ir=1,dimOB
   do ip=1,dimOB
     ipr = ipr + 1
     read(iunit,rec=ip+(ir-1)*dimOB) ints(1:dimOA*dimOA)

     do j=1,dimOA
        do i=1,dimOA
           work(i,j) = ints(i+(j-1)*dimOA)
        enddo
     enddo
     val = val + sum(work(1:dimOA,1:dimOA)*tmpAB(1:dimOA,1:dimOA,ip,ir))

   enddo
enddo
close(iunit)
tNaNb = -2*val
print*, 'tNaNb',tNaNb*1000

exchs2 = tElst + sum(tvk) + sum(tNa) + sum(tNb) + tNaNb
SAPT%exchs2 = exchs2
!write(LOUT,'(/1x,a,f16.8)') 'ExchS2      = ', exchs2*1000d0
call print_en('ExchS2',exchs2*1000,.true.)

deallocate(tmpAB)

deallocate(ints,work)
deallocate(intB,intA)
deallocate(Sab,Vbba,Vaab,Vbaa,Vabb)

call delfile('OOOOAABB')
call delfile('OOOOBBBA')
call delfile('OOOOAAAB')

end subroutine e1exch_Chol

end module sapt_Chol_exch
