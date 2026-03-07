!#define SAPT_OS_DEBUG 10

module sapt_open
use types
use sapt_utils

implicit none

contains

subroutine e1exchs2_os(Flags,A,B,SAPT)
!
! based on e1exch_NaNb
!
type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer          :: NAO,NBas,dimOA,dimOB
integer          :: i,j,ip,iq,ir,it,iu,is
integer          :: iunit
double precision :: fac,val,t1
double precision :: tvk(3),tNa(2),tNb(2),tNaNb
double precision :: exchs2

double precision,allocatable :: gaSga(:,:),gbSgb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Vaab(:,:),Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: Vbab(:,:),Vbba(:,:)
double precision,allocatable :: intAa(:,:,:,:),intAb(:,:,:,:)
double precision,allocatable :: intBa(:,:,:,:),intBb(:,:,:,:)
double precision,allocatable :: int3(:,:,:,:),tmpAB(:,:,:,:)
!double precision,allocatable :: int3a(:,:,:,:),int3b(:,:,:,:)
double precision,allocatable :: intABa(:,:),intABb(:,:)
double precision,allocatable :: intBAa(:,:),intBAb(:,:)
double precision,allocatable :: ints(:)
double precision,allocatable :: work(:,:),tmp(:,:)
! test
integer :: k,l,kl
double precision :: sachrg,saspin

! set dimensions
NAO   = SAPT%NAO
NBas  = A%NBasis
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

!#if SAPT_OS_DEBUG > 1
  print*, 'dimensions in e1exchs2_os:'
  print*, 'NAO, NBas', NAO, NBas
!#endif

allocate(S(NBas,NBas),Va(NBas,NBas),Vb(NBas,NBas))
allocate(Sab(NBas,NBas),&
         Vabb(NBas,NBas),Vbaa(NBas,NBas),&
         Vaab(NBas,NBas),Vbba(NBas,NBas),&
         Vbab(NBas,NBas))

call get_one_mat('V',Va,A%Monomer,NBas)
call get_one_mat('V',Vb,B%Monomer,NBas)

call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)
call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
!call tran2MO(Vb,B%CMO,A%CMO,Vbba,NBas)
call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

call get_one_mat('S',S,A%Monomer,NBas)
call tran2MO(S,A%CMO,B%CMO,Sab,NBas)

!allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
!         RDM2Bval(dimOB,dimOB,dimOB,dimOB))

deallocate(Vb,Va,S)

! *** begin test
sachrg=0
saspin=0
do j=1,NBas
do i=1,NBas
   sachrg = sachrg + (A%g1a(i,j) + A%g1b(i,j))**2
   saspin = saspin + (A%g1a(i,j) - A%g1b(i,j))**2
enddo
enddo
print*, 'Charge A:',sachrg
print*, 'Spin   A', saspin
sachrg=0
saspin=0
do j=1,NBas
do i=1,NBas
   sachrg = sachrg + (B%g1a(i,j) + B%g1b(i,j))**2
   saspin = saspin + (B%g1a(i,j) - B%g1b(i,j))**2
enddo
enddo
print*, 'Charge B:',sachrg
print*, 'Spin   B', saspin
! *** end test

! prepare intermediates
allocate(intAa(dimOA,dimOA,dimOA,dimOB),intAb(dimOA,dimOA,dimOA,dimOB))


! test intA
call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%g2aaba+A%g2bbab,dimOA**3,Sab,NBas,0d0,intAa,dimOA**3)
print*, 'test intA = ',norm2(intAa)

call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%g2aaba,dimOA**3,Sab,NBas,0d0,intAa,dimOA**3)
call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%g2bbab,dimOA**3,Sab,NBas,0d0,intAb,dimOA**3)
print*, 'intAa = ',norm2(intAa)
print*, 'intAb = ',norm2(intAb)

allocate(intBa(dimOB,dimOB,dimOA,dimOB),intBb(dimOB,dimOB,dimOA,dimOB))
! careful! intB(B,B,A,B)
do is=1,dimOB
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%g2aaba(:,:,:,is),dimOB**2,Sab,NBas,0d0,intBa(:,:,:,is),dimOB**2)
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%g2bbab(:,:,:,is),dimOB**2,Sab,NBas,0d0,intBb(:,:,:,is),dimOB**2)
enddo
print*, 'intBa = ',norm2(intBa)
print*, 'intBb = ',norm2(intBb)

! allocate integrals
allocate(ints(NBas*NBas),work(NBas,NBas))

! t1 = (g1a^A.g1a^B + g1b^A.g1b^B) Sab.Sab
allocate(tmp(Nbas,NBas))
allocate(gaSga(NBas,NBas),gbSgb(NBas,NBas))

call dgemm('N','N',NBas,NBas,NBas,1d0,A%g1a,NBas,Sab,NBas,0d0,work,NBas)
call dgemm('N','N',NBas,NBas,NBas,1d0,work,NBas,B%g1a,NBas,0d0,gaSga,NBas)
t1  = 0
do j=1,dimOB
   do i=1,dimOA
      t1 = t1 + gaSga(i,j)*Sab(i,j)
   enddo
enddo

call dgemm('N','N',NBas,NBas,NBas,1d0,A%g1b,NBas,Sab,NBas,0d0,work,NBas)
call dgemm('N','N',NBas,NBas,NBas,1d0,work,NBas,B%g1b,NBas,0d0,gbSgb,NBas)
do j=1,dimOB
   do i=1,dimOA
      t1 = t1 + gbSgb(i,j)*Sab(i,j)
   enddo
enddo

print*, 'SAPT%elst',SAPT%elst
print*, 'SAPT%Vnn ',SAPT%Vnn
print*, 't1',t1
t1 = (SAPT%elst-SAPT%Vnn)*t1
print*, 't1 = ',t1*1000
print*, 'test: ',SAPT%elst*t1
print*, ''

! tvk(1) = (g1a^A.g1a^B + g1b^A.g1b^B) v^A Sab
tvk = 0
do j=1,dimOB
   do i=1,dimOA
      tvk(1) = tvk(1) - gaSga(i,j)*Vaab(i,j) &
                      - gbSgb(i,j)*Vaab(i,j)
   enddo
enddo

print*, 'tvk(1) ', tvk(1)*1000

! tvk(2) = (g1a^A.g1a^B + g1b^A.g1b^B) v^B Sab
call dgemm('N','N',NBas,NBas,NBas,1d0,A%g1a,NBas,Vbab,NBas,0d0,work,NBas)
call dgemm('N','N',NBas,NBas,NBas,1d0,work,NBas,B%g1a,NBas,0d0,gaSga,NBas)

call dgemm('N','N',NBas,NBas,NBas,1d0,A%g1b,NBas,Vbab,NBas,0d0,work,NBas)
call dgemm('N','N',NBas,NBas,NBas,1d0,work,NBas,B%g1b,NBas,0d0,gbSgb,NBas)

do j=1,dimOB
   do i=1,dimOA
      tvk(2) = tvk(2) - gaSga(i,j)*Sab(i,j) &
                      - gbSgb(i,j)*Sab(i,j)
   enddo
enddo
print*, 'tvk(2) ', tvk(2)*1000

! tvk(3) = (g1a^A.g1a^B + g1b^A.g1b^B) v_pq^rs
open(newunit=iunit,file='FOFOABBA',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*NBas*dimOB)

work = 0
ints = 0
kl = 0
do l=1,dimOA
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) ints(1:NBas*dimOB)

      if(k<=dimOB) then
         do j=1,dimOB
            do i=1,dimOA
               work(i,j) = ints(i+(j-1)*NBas)
            enddo
         enddo
         ! k = s
         ! l = q
         do j=1,dimOB
            do i=1,dimOA
               tvk(3) = tvk(3) + A%g1a(i,l)*B%g1a(j,k)*work(i,j) &
                               + A%g1b(i,l)*B%g1b(j,k)*work(i,j)
            enddo
         enddo
      endif

   enddo
enddo
close(iunit)
tvk(3) = -tvk(3)
print*, 'tvk(3)-test : ',tvk(3)*1000

! how to do it in AO basis?
!call dgemm('N','N',NBas,NBas,NBas,1d0,A%CMO,NBas,A%g1a,NBas,0d0,tmp,NBas)
!call dgemm('N','T',NBas,NBas,NBas,1d0,tmp,NBas,A%CMO,NBas,0d0,work,NBas)
!
!call dgemm('N','N',NBas,NBas,NBas,1d0,A%CMO,NBas,A%g1b,NBas,0d0,tmp,NBas)
!call dgemm('N','T',NBas,NBas,NBas,1d0,tmp,NBas,A%CMO,NBas,1d0,work,NBas)
!do j=1,NBas
!   do i=1,NBas
!      tvk(3) = tvk(3) - work(i,j)*B%Kmat(j,i)
!   enddo
!enddo
!print*, 'tvk(3) ', tvk(3)*1000

deallocate(tmp)

! tNa(1) = (g1a^A.Naa^B + g1b^A.Nbb^B) v^A Sab
allocate(intABa(NBas,NBas),intABb(NBas,NBas))
intABa = 0
intABb = 0
do it=1,dimOB
   do iq=1,dimOA
      do ir=1,dimOA
         do ip=1,dimOA
            intABa(iq,it) = intABa(iq,it) + intAa(ip,ir,iq,it)*Vbaa(ip,ir)
            intABb(iq,it) = intABb(iq,it) + intAb(ip,ir,iq,it)*Vbaa(ip,ir)
         enddo
      enddo
   enddo
enddo
! call dgemv('N','N',dimOA,dimOB)

tNa = 0
do iu=1,dimOB
   do it=1,dimOB
      do iq=1,dimOA
         tNa(1) = tNA(1) - B%g1a(it,iu)*intABa(iq,it)*Sab(iq,iu) &
                         - B%g1b(it,iu)*intABb(iq,it)*Sab(iq,iu)
      enddo
   enddo
enddo
print*, 'tNa-1',tNa(1)*1000

!allocate(int3a(dimOA,dimOA,dimOA,dimOB),int3b(dimOA,dimOA,dimOA,dimOB))
allocate(int3(dimOA,dimOA,dimOA,dimOB))

call dgemm('N','N',dimOA**3,dimOB,dimOB,1d0,intAa,dimOA**3,B%g1a,NBas,0d0,int3,dimOA**3)
call dgemm('N','N',dimOA**3,dimOB,dimOB,1d0,intAb,dimOA**3,B%g1b,NBas,1d0,int3,dimOA**3)
!(FO|FO): (AA|AB)
open(newunit=iunit,file='FOFOAAAB',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

ints = 0
do iu=1,dimOB
   do iq=1,dimOA
      read(iunit,rec=iq+(iu-1)*NBas) ints(1:NBas*dimOA)

      do ir=1,dimOA
         do ip=1,dimOA
            tNa(2) = tNa(2) + int3(ip,ir,iq,iu)*ints(ip+(ir-1)*NBas) !&
!                            + int3b(ip,ir,iq,iu)*ints(ip+(ir-1)*NBas)
         enddo
      enddo

   enddo
enddo
tNa(2) = -tNa(2)
print*, 'tNa-2',tNa(2)*1000

close(iunit)

! tNb(1) = (g1a^B.Naa^A + g1b^B.Nbb^A) v^B Sab
intABa = 0
intABb = 0
do iq=1,dimOB
   do it=1,dimOA
      do ir=1,dimOB
         do ip=1,dimOB
            intABa(it,iq) = intABa(it,iq) + intBa(ip,ir,it,iq)*Vabb(ip,ir)
            intABb(it,iq) = intABb(it,iq) + intBb(ip,ir,it,iq)*Vabb(ip,ir)
         enddo
      enddo
   enddo
enddo
tNb = 0
do iu=1,dimOA
   do it=1,dimOA
      do iq=1,dimOB
         tNb(1) = tNb(1) - A%g1a(it,iu)*intABa(it,iq)*Sab(iu,iq) &
                         - A%g1b(it,iu)*intABb(it,iq)*Sab(iu,iq)
      enddo
   enddo
enddo
print*, 'tNb-1',tNb(1)*1000

deallocate(int3)
allocate(int3(dimOB,dimOB,dimOA,dimOB))

do is=1,dimOB
   call dgemm('N','N',dimOB**2,dimOA,dimOA,1d0,intBa(:,:,:,is),dimOB**2,A%g1a,NBas,0d0,int3(:,:,:,is),dimOB**2)
   call dgemm('N','N',dimOB**2,dimOA,dimOA,1d0,intBb(:,:,:,is),dimOB**2,A%g1b,NBas,1d0,int3(:,:,:,is),dimOB**2)
enddo
!(FO|FO): (BB|BA)
open(newunit=iunit,file='FOFOBBBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

ints = 0
do iu=1,dimOA
   do iq=1,dimOB
      read(iunit,rec=iq+(iu-1)*NBas) ints(1:NBas*dimOB)

      do ir=1,dimOB
         do ip=1,dimOB
            tNb(2) = tNb(2) + int3(ip,ir,iu,iq)*ints(ip+(ir-1)*NBas)
         enddo
      enddo

   enddo
enddo
tNb(2) = -tNb(2)
print*, 'tNb-2',tNb(2)*1000

close(iunit)

deallocate(int3)
deallocate(intABb,intABa)

open(newunit=iunit,file='FOFOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA*NBas)
!
allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))
call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intAa,dimOA**2,intBa,dimOB**2,0d0,tmpAB,dimOA**2)
call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intAb,dimOA**2,intBb,dimOB**2,1d0,tmpAB,dimOA**2)

val  = 0
work = 0

ints = 0
do ir=1,dimOB
   do ip=1,dimOB
     read(iunit,rec=ip+(ir-1)*NBas) ints(1:dimOA*NBas)

     do j=1,dimOA
        do i=1,dimOA
           work(i,j) = ints(i+(j-1)*NBas)
        enddo
     enddo
     val = val + sum(work(1:dimOA,1:dimOA)*tmpAB(1:dimOA,1:dimOA,ip,ir))

   enddo
enddo
close(iunit)
tNaNb = -val
print*, 'tNaNb',tNaNb*1000

exchs2 = t1 + sum(tvk) + sum(tNa) + sum(tNb) + tNaNb

close(iunit)
!
call print_en('ExchS2',exchs2*1000,.true.)
SAPT%exchs2 = exchs2

deallocate(intAb,intAa,intBb,intBa)
deallocate(ints,work)
deallocate(tmpAB)
deallocate(gbSgb,gaSga)
deallocate(Sab,Vbba,Vaab,Vbaa,Vabb)
deallocate(Vbab)

end subroutine e1exchs2_os

end module sapt_open
