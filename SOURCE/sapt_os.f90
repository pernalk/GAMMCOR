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

subroutine e2exind_o(Flags,A,B,SAPT)
!
! uncoupled E2exch-ind
! see Eq. (21) in 2012 JCP paper
!
use timing

implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO,NMO
integer :: i

real(8) :: e2exi_ba(4),e2exi_ab(4)
real(8) :: e2exi_ab_sum,e2exi_ba_sum,e2exi

real(8), allocatable :: S(:,:)
real(8), allocatable :: Oa(:,:),Ob(:,:)
real(8), allocatable :: JOa(:,:),JOb(:,:)
real(8), allocatable :: JPBSOa(:,:),JPBSOb(:,:)
real(8), allocatable :: JPASOa(:,:),JPASOb(:,:)
real(8), allocatable :: Paa(:,:),Pab(:,:),Pba(:,:),Pbb(:,:)
real(8), allocatable :: work(:,:)

! set dimensions
NAO = SAPT%NAO
NMO = A%NBasis

allocate(PAa(nao,nao),PAb(nao,nao))
allocate(PBa(nao,nao),PBb(nao,nao))
call get_den(nmo,A%UMO(:,:,1),A%Uocc(:,1),1d0,PAa)
call get_den(nmo,A%UMO(:,:,2),A%Uocc(:,2),1d0,PAb)
call get_den(nmo,B%UMO(:,:,1),B%Uocc(:,1),1d0,PBa)
call get_den(nmo,B%UMO(:,:,2),B%Uocc(:,2),1d0,PBb)

allocate(S(nao,nao),work(nao,nao))
call get_one_mat('S',S,A%Monomer,nao)

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'KAa =', norm2(A%Ka)
   print*, 'KAb =', norm2(A%Kb)
   print*, '-----'
   print*, 'KBa =', norm2(B%Ka)
   print*, 'KBb =', norm2(B%Kb)
   print*, '-----'
endif

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'hAa =', norm2(A%ha)
   print*, 'hAb =', norm2(A%hb)
   print*, '-----'
   print*, 'hBa =', norm2(B%ha)
   print*, 'hBb =', norm2(B%hb)
   print*, '-----'
endif

! o matrices
allocate(Oa(nao,nao),Ob(nao,nao))
call dgemm('N','N',nao,nao,nao,1d0,PAa,nao,S,nao,0d0,work,nao)
call dgemm('N','N',nao,nao,nao,1d0,work,nao,PBa,nao,0d0,Oa,nao)
call dgemm('N','N',nao,nao,nao,1d0,PAb,nao,S,nao,0d0,work,nao)
call dgemm('N','N',nao,nao,nao,1d0,work,nao,PBb,nao,0d0,Ob,nao)

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'Oa =', norm2(Oa)
   print*, 'Ob =', norm2(Ob)
   print*, '-----'
   print*, 'tAa',norm2(A%ta)
   print*, 'tBa',norm2(A%tb)
   print*, '-----'
endif

!!!!!!!!!!!!!!!!!!!!
! E2exch-ind(A<--B)
!!!!!!!!!!!!!!!!!!!!
e2exi_ba=0d0
allocate(JOa(nao,nao),JOb(nao,nao))
allocate(JPBSOa(nao,nao),JPBSOb(nao,nao))
! alpha-alpha
call e2exi_terms_SameSpin(e2exi_ba(1),A%UMO(:,:,1),PAa,PBa,S,A%ta,A%WPot,B%WPot,&
                          A%ha,B%ha,B%Ka,Oa,JOa,JPBSOa,A%NOa,A%NVa,nao,nmo)
! beta-beta
call e2exi_terms_SameSpin(e2exi_ba(2),A%UMO(:,:,2),PAb,PBb,S,A%tb,A%WPot,B%WPot,&
                          A%hb,B%hb,B%Kb,Ob,JOb,JPBSOb,A%NOb,A%NVb,nao,nmo)
! alpha-beta/beta-alpha
call e2exi_terms_OppSpin(e2exi_ba(3),A%UMO(:,:,1),A%UMO(:,:,2),A%ta,A%tb,JOa,JOb,&
                         JPBSOa,JPBSOb,A%NOa,A%NVa,A%NOb,A%NVb,nao,nmo)

e2exi_ba_sum = sum(e2exi_ba)
call print_en('E2exch-ind(A<--B)',e2exi_ba_sum*1.0d3,.true.)
deallocate(JPBSOb,JPBSOa)

!!!!!!!!!!!!!!!!!!!!
! E2exch-ind(A-->B)
!!!!!!!!!!!!!!!!!!!!
e2exi_ab=0d0

! transpose O
work = Oa
Oa = transpose(work)
work = Ob
Ob = transpose(work)

allocate(JPASOa(nao,nao),JPASOb(nao,nao))
! alpha-alpha
call e2exi_terms_SameSpin(e2exi_ab(1),B%UMO(:,:,1),PBa,PAa,S,B%ta,B%WPot,A%WPot,&
                          B%ha,A%ha,A%Ka,Oa,JOa,JPASOa,B%NOa,B%NVa,nao,nmo)
! beta-beta
call e2exi_terms_SameSpin(e2exi_ab(2),B%UMO(:,:,2),PBb,PAb,S,B%tb,B%WPot,A%WPot,&
                          B%hb,A%hb,A%Kb,Ob,JOb,JPASOb,B%NOb,B%NVb,nao,nmo)
! alpha-beta/beta-alpha
call e2exi_terms_OppSpin(e2exi_ab(3),B%UMO(:,:,1),B%UMO(:,:,2),B%ta,B%tb,JOa,JOb,&
                         JPASOa,JPASOb,B%NOa,B%NVa,B%NOb,B%NVb,nao,nmo)

e2exi_ab_sum = sum(e2exi_ab)
call print_en('E2exch-ind(B<--A)',e2exi_ab_sum*1.0d3,.true.)

e2exi = e2exi_ba_sum + e2exi_ab_sum
SAPT%e2exind = e2exi
call print_en('E2exch-ind',e2exi*1.0d3,.false.)

deallocate(JPASOb,JPASOa)
deallocate(Ob,Oa)
deallocate(work,S)
deallocate(PBb,PBa,PAb,PAa)

end subroutine e2exind_o

subroutine e2exi_terms_SameSpin(e2exi_ss,CA,PA,PB,S,tA,WA,WB,hA,hB,KB,O,JO,JPBSO,noa,nva,nao,nmo)
!
! Same spin
! returns E2ind(A<--B)
! together with JO JPBSO matrices
!
! see Eq. (17) in https://doi.org/10.1063/5.0090688
!
integer,intent(in) :: noa,nva
integer,intent(in) :: nao,nmo
real*8,intent(in)  :: tA(noa*nva)
real*8,intent(in)  :: PA(nao,nao),PB(nao,nao)
real*8,intent(in)  :: CA(nao,nmo),S(nao,nao)
real*8,intent(in)  :: WA(nao,nao),WB(nao,nao)
real*8,intent(in)  :: hA(nao,nao),hB(nao,nao)
real*8,intent(in)  :: KB(nao,nao),O(nao,nao)
!
real*8,intent(out) :: e2exi_ss
real*8,intent(out) :: JO(nao,nao),JPBSO(nao,nao)

integer :: i,j,ij
real*8  :: intvo(nva,noa)
real*8  :: PBS(nao,nao),PBSO(nao,nao)
real*8  :: SPA(nao,nao),SPAWB(nao,nao)
real*8  :: WAPBS(nao,nao),WBPAS(nao,nao)
real*8  :: Ko(nao,nao)
real*8  :: intao(nao,nao)
real*8  :: intA(nao,nao),intB(nao,nao)
real*8 ::  tAao(nao,nao),tmp(nao,nva)
real*8,allocatable :: work(:,:)

double precision,external  :: trace

call dgemm('N','N',nao,nao,nao,1d0,PB,nao,S,nao,0d0,PBS,nao)
call dgemm('N','N',nao,nao,nao,1d0,PBS,nao,O,nao,0d0,PBSO,nao)

call dgemm('N','N',nao,nao,nao,1d0,S,nao,PA,nao,0d0,SPA,nao)
call dgemm('N','N',nao,nao,nao,1d0,SPA,nao,WB,nao,0d0,SPAWB,nao)

call dgemm('N','N',nao,nao,nao,1d0,WA,nao,PBS,nao,0d0,WAPBS,nao)
call dgemm('N','T',nao,nao,nao,1d0,WB,nao,SPA,nao,0d0,WBPAS,nao)

call make_J1(nao,O,Jo,'AOTWOSORT')
call make_J1(nao,PBSO,JPBSO,'AOTWOSORT')
call make_K(nao,O,Ko,'AOTWOSORT')

!print*, 'Jo=', norm2(Jo)
!print*, 'Ko=', norm2(Ko)

! test terms 1
intao=KB+Jo-Ko-JPBSO

allocate(work(nao,nao))
! test terms 2
!call dgemm('T','N',nao,nao,nao,1d0,PBS,nao,hA,nao,0d0,intA,nao)
!call dgemm('T','N',nao,nao,nao,1d0,PBS,nao,SPAWB,nao,-1d0,intA,nao)
!call dgemm('T','N',nao,nao,nao,1d0,PBS,nao,WAPBS,nao,-1d0,intA,nao)
!call dgemm('T','T',nao,nao,nao,1d0,PBS,nao,Ko,nao,0d0,intA,nao)

work = transpose(Ko)
work = work + hA - SPAWB- WAPBS
call dgemm('T','N',nao,nao,nao,1d0,PBS,nao,work,nao,0d0,intA,nao)

! test terms 3
!call dgemm('N','N',nao,nao,nao,1d0,hB,nao,PBS,nao,0d0,intB,nao)
!call dgemm('N','N',nao,nao,nao,1d0,WBPAS,nao,PBS,nao,0d0,intB,nao)
!call dgemm('N','N',nao,nao,nao,1d0,Ko,nao,PBS,nao,0d0,intB,nao)
!call dgemm('N','N',nao,nao,nao,1d0,Ko,nao,PBS,nao,0d0,intB,nao)

work = 0d0
work = hB-WBPAS+Ko
call dgemm('N','N',nao,nao,nao,1d0,work,nao,PBS,nao,0d0,intB,nao)

intao=intao+intA+intB

deallocate(work)

!allocate(work(nva,nao))
!! intvo=Cvirt^T(nao,nvirt).inter(nao,nao).Cocc(nao,occ)
!call dgemm('T','N',nva,nao,nao,1d0,CA(:,noa+1:nmo),nao,intao,nao,0d0,work,nva)
!call dgemm('N','N',nva,noa,nao,1d0,work,nva,CA(:,1:noa),nao,0d0,intvo,nva)
!!
!!print*, 'intvo =',norm2(intvo)
!!
!!! e2exi = xA.intvo
!e2exi_ss = 0d0
!ij = 0
!do j=1,noa
!   do i=1,nva
!      ij = ij + 1
!      !print*, 'tA(ij)',ij,tA(ij),intvo(i,j)
!      e2exi_ss = e2exi_ss + tA(ij)*intvo(i,j)
!   enddo
!enddo
!print*, 'e2exch-ind aa =', e2exi_ss
!
!deallocate(work)

call tranMO2AO_vovo(tAao,tA,CA,noa,nva,nao,nmo)

allocate(work(nao,nao))
call dgemm('T','N',nao,nao,nao,1d0,tAao,nao,intao,nao,0d0,work,nao)
e2exi_ss = 0d0
do i=1,nao
      e2exi_ss = e2exi_ss - work(i,i)
enddo
print*, 'e2exch-ind aa =', e2exi_ss

deallocate(work)

end subroutine e2exi_terms_SameSpin

subroutine tranMO2AO_vovo(tout,tin,C,no,nv,nao,nmo)
integer,intent(in) :: no,nv,nao,nmo
real*8,intent(in)  :: tin(nv,no),C(nao,nmo)
real*8,intent(out) :: tout(nao,nao)

integer :: i,j,ij
real*8  :: tmp(nao,nv)

call dgemm('N','T',nao,nv,no,1d0,C(1:nao,1:no),nao,tin,nv,0d0,tmp,nao)
call dgemm('N','T',nao,nao,nv,1d0,tmp,nao,C(1:nao,no+1:nao),nao,0d0,tout,nao)

end subroutine tranMO2AO_vovo

subroutine e2exi_terms_OppSpin(e2exi_os,CAa,CAb,tAa,tAb,JOa,JOb,JPBSOa,JPBSOb,noa_a,nva_a,noa_b,nva_b,nao,nmo)
!
! Opposite ab spin
! returns E2ind(A<--B)
!
integer,intent(in) :: noa_a,nva_a,noa_b,nva_b
integer,intent(in) :: nao,nmo
real*8,intent(in)  :: CAa(nao,nmo),CAb(nao,nao)
real*8,intent(in)  :: tAa(noa_a*nva_a),tAb(noa_b*nva_b)
real*8,intent(in)  :: JOa(nao,nao),JOb(nao,nao)
real*8,intent(in)  :: JPBSOa(nao,nao),JPBSOb(nao,nao)

real*8,intent(out) :: e2exi_os

integer :: i
real*8 :: vala,valb
real*8 :: tAao(nao,nao),intao(nao,nao)
real(8), allocatable :: work(:,:)

! t-alpha
allocate(work(nao,nao))
intao = JOb-JPBSOb
call tranMO2AO_vovo(tAao,tAa,CAa,noa_a,nva_a,nao,nmo)
call dgemm('T','N',nao,nao,nao,1d0,tAao,nao,intao,nao,0d0,work,nao)

vala = 0d0
do i=1,nao
      vala = vala - work(i,i)
enddo
print*, 'e2exch-ind ab =', vala

! t-beta
intao = JOa-JPBSOa
call tranMO2AO_vovo(tAao,tAb,CAb,noa_b,nva_b,nao,nmo)
call dgemm('T','N',nao,nao,nao,1d0,tAao,nao,intao,nao,0d0,work,nao)

valb = 0d0
do i=1,nao
      valb = valb - work(i,i)
enddo
print*, 'e2exch-ind ba =', valb

e2exi_os = vala + valb
print*, 'e2exch-ind os =', e2exi_os

deallocate(work)

end subroutine e2exi_terms_OppSpin

subroutine e2exdisp_o(Flags,A,B,SAPT)
!
! uncoupled E2exch-disp
! see Eq. (23) in https://doi.org/10.1063/1.4758455
!
use timing

implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer :: NBasis

double precision, allocatable :: work(:,:)
double precision, allocatable :: Sa(:,:),Sb(:,:)
double precision, allocatable :: Saboo_aa(:,:),Sabov_aa(:,:),Sabvo_aa(:,:)
double precision, allocatable :: Sbaoo_aa(:,:),Sabvv_aa(:,:),Sbavo_aa(:,:)
double precision, allocatable :: Sat(:,:),Sbt(:,:)
double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)
double precision, allocatable :: Waa_ov(:,:),Wba_ov(:,:)
double precision, allocatable :: ints(:),amps(:)

integer :: i,j
real*8 :: e2xd_vterms(8),e2xd_oterms(4)
real*8 :: e2xd_aa,e2xd_bb,e2xd_ab,e2xd_ba
real*8 :: e2xd_ssv,e2xd_sso
real*8 :: e2exd_unc

NBasis = A%NBasis

! omega potential
allocate(Waa(NBasis,NBasis),Wab(NBasis,NBasis))
allocate(Wba(NBasis,NBasis),Wbb(NBasis,NBasis))

call tran2MO(A%WPot,B%UMO(:,:,1),B%UMO(:,:,1),Waa,NBasis)
call tran2MO(A%WPot,B%UMO(:,:,2),B%UMO(:,:,2),Wab,NBasis)

call tran2MO(B%WPot,A%UMO(:,:,1),A%UMO(:,:,1),Wba,NBasis)
call tran2MO(B%WPot,A%UMO(:,:,2),A%UMO(:,:,2),Wbb,NBasis)

! S matrix
allocate(work(NBasis,NBasis))
allocate(Sa(NBasis,NBasis),Sb(NBasis,NBasis))
allocate(Sat(NBasis,NBasis),Sbt(NBasis,NBasis))

call get_one_mat('S',work,A%Monomer,NBasis)
call tran2MO(work,A%UMO(:,:,1),B%UMO(:,:,1),Sa,NBasis)
call tran2MO(work,A%UMO(:,:,2),B%UMO(:,:,2),Sb,NBasis)

deallocate(work)

e2xd_aa=0d0 ; e2xd_bb=0d0
e2xd_ab=0d0 ; e2xd_ba=0d0
! alpha-alpha
allocate(ints(A%NVa*A%NOa*B%NVa*B%NOa))
allocate(amps(A%NVa*A%NOa*B%NVa*B%NOa))
call load_vovo('OVOVABaa',ints,B%IndNa,A%NVa,A%NOa,B%NVa,B%NOa)
call calc_amps(amps,ints,A%UOrbE(:,1),B%UOrbE(:,1),A%NVa,A%NOa,B%NVa,B%NOa,NBasis)
!print*, 'vovo aaaa = ', norm2(ints)
!print*, 'amps aaaa = ', norm2(amps)
e2xd_ssv=0 ; e2xd_sso=0d0 
call e2exd_vterms_SameSpin(e2xd_ssv,ints,amps,Sa,A%NOa,A%NVa,B%NOa,B%NVa,NBasis)
call e2exd_oterms_SameSpin(e2xd_sso,ints,amps,Waa,Wba,Sa,A%NOa,A%NVa,B%NOa,B%NVa,NBasis)

e2xd_aa = e2xd_ssv + e2xd_sso

deallocate(amps,ints)

! beta-beta
allocate(ints(A%NVb*A%NOb*B%NVb*B%NOb))
allocate(amps(A%NVb*A%NOb*B%NVb*B%NOb))
call load_vovo('OVOVABbb',ints,B%IndNb,A%NVb,A%NOb,B%NVb,B%NOb)
call calc_amps(amps,ints,A%UOrbE(:,2),B%UOrbE(:,2),A%NVb,A%NOb,B%NVb,B%NOb,NBasis)

e2xd_ssv=0 ; e2xd_sso=0d0 
call e2exd_vterms_SameSpin(e2xd_ssv,ints,amps,Sb,A%NOb,A%NVb,B%NOb,B%NVb,NBasis)
call e2exd_oterms_SameSpin(e2xd_sso,ints,amps,Wab,Wbb,Sb,A%NOb,A%NVb,B%NOb,B%NVb,NBasis)

e2xd_bb = e2xd_ssv + e2xd_sso

deallocate(amps,ints)

! alpha-beta
allocate(ints(A%NVa*A%NOa*B%NVb*B%NOb))
allocate(amps(A%NVa*A%NOa*B%NVb*B%NOb))

call load_vovo('OVOVABab',ints,B%IndNb,A%NVa,A%NOa,B%NVb,B%NOb)
call calc_amps(amps,ints,A%UOrbE(:,1),B%UOrbE(:,2),A%NVa,A%NOa,B%NVb,B%NOb,NBasis)
call e2exd_ovterms_OppSpin(e2xd_ab,ints,amps,Wab,Wba,Sa,Sb,A%NOa,A%NVa,A%NOb,A%NVb,&
                           B%NOa,B%NVa,B%NOb,B%NVb,NBasis)

deallocate(amps,ints)

! beta-alpha
allocate(ints(A%NVb*A%NOb*B%NVa*B%NOa))
allocate(amps(A%NVb*A%NOb*B%NVa*B%NOa))

call load_vovo('OVOVABba',ints,B%IndNa,A%NVb,A%NOb,B%NVa,B%NOa)
call calc_amps(amps,ints,A%UOrbE(:,2),B%UOrbE(:,1),A%NVb,A%NOb,B%NVa,B%NOa,NBasis)
call e2exd_ovterms_OppSpin(e2xd_ba,ints,amps,Waa,Wbb,Sb,Sa,A%NOb,A%NVb,A%NOa,A%NVa, &
                           B%NOb,B%NVb,B%NOa,B%NVa,NBasis)

deallocate(amps,ints)

print*, ''
print*, 'E2exd components in Ha:'
print*, 'aa =', e2xd_aa
print*, 'bb =', e2xd_bb
print*, 'ab =', e2xd_ab
print*, 'ba =', e2xd_ba

! sum aa, bb, ab, ba terms
e2exd_unc = e2xd_aa + e2xd_bb + e2xd_ab + e2xd_ba
SAPT%e2exdisp_unc = e2exd_unc
call print_en('E2exch-disp(unc)',e2exd_unc*1.0d3,.true.)

deallocate(Wbb,Waa,Wab,Wba)

end subroutine e2exdisp_o

subroutine e2exd_vterms_SameSpin(e2xd_vsum,ints,amps,Sa,noA,nvA,noB,nvB,NBasis)
!
! compute v-like contributions to E2exch-disp (8),
! see Eq. (23) in https://doi.org/10.1063/1.4758455 
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB,NBasis
real*8, intent(in)  :: Sa(NBasis,NBasis)
real*8, intent(in)  :: ints(nva*noa*nvb*nob)
real*8, intent(in)  :: amps(nva*noa*nvb*nob)

real*8, intent(out) :: e2xd_vsum

integer :: i,j
real*8  :: e2xd_vterms(8)

double precision, allocatable :: Saboo(:,:),Sabov(:,:),Sabvo(:,:)
double precision, allocatable :: Sbaoo(:,:),Sabvv(:,:),Sbavo(:,:)
double precision, allocatable :: Waa_ov(:,:),Wba_ov(:,:)

allocate(Sabov(noa,nvb),Sabvo(nva,nob))
allocate(Saboo(noa,nob),Sbaoo(nob,noa))
allocate(Sabvv(nva,nvb))
allocate(Sbavo(nvb,noa))
! Sov
do j=1,nvb
   do i=1,noa
      Sabov(i,j) = Sa(i,nob+j)
   enddo
enddo
! Svo
do j=1,nob
   do i=1,nva
      Sabvo(i,j) = Sa(noa+i,j)
   enddo
enddo
! Soo
do j=1,nob
   do i=1,noa
      Saboo(i,j) = Sa(i,j)
   enddo
enddo
Sbaoo = transpose(Saboo)
! Svv
do j=1,nvb
   do i=1,nva
      Sabvv(i,j) = Sa(noa+i,nob+j)
   enddo
enddo
!Sbavo
do j=1,noa
   do i=1,nvb
      Sbavo(i,j) = Sa(j,nob+i)
   enddo
enddo

e2xd_vterms = 0d0
call vterm1_e2exd(e2xd_vterms(1),amps,ints,Sabvv,nva,noa,nvb,nob)
print*, 'e2xd T1 =', e2xd_vterms(1)*1d3

call vterm2_e2exd(e2xd_vterms(2),amps,ints,Sabov,nva,noa,nvb,nob)
call vterm4_e2exd(e2xd_vterms(4),amps,ints,Sabvo,nva,noa,nvb,nob)
print*, 'e2xd T2 =', e2xd_vterms(2)*1d3
print*, 'e2xd T4 =', e2xd_vterms(4)*1d3

call vterm3_e2exd(e2xd_vterms(3),amps,ints,nva,noa,nvb,nob,Sbavo,noa)
call vterm5_e2exd(e2xd_vterms(5),amps,ints,nva,noa,nvb,nob,Sabvo,nob)
print*, 'e2xd T3 =', e2xd_vterms(3)*1d3
print*, 'e2xd T5 =', e2xd_vterms(5)*1d3

call vterm6_e2exd(e2xd_vterms(6),amps,ints,Sbaoo,nva,noa,nvb,nob)
print*, 'e2xd T6 =', e2xd_vterms(6)*1d3

call vterm7_e2exd(e2xd_vterms(7),amps,ints,nva,noa,nvb,nob,Saboo,noa)
call vterm8_e2exd(e2xd_vterms(8),amps,ints,nva,noa,nvb,nob,Saboo,nob)
print*, 'e2xd T7 =', e2xd_vterms(7)*1d3
print*, 'e2xd T8 =', e2xd_vterms(8)*1d3

deallocate(Saboo,Sbaoo,Sabvv)
deallocate(Sabvo,Sabov,Sbavo)

e2xd_vsum = sum(e2xd_vterms)

print*, 'E2exd vterms =', e2xd_vsum*1d3

end subroutine e2exd_vterms_SameSpin

subroutine e2exd_oterms_SameSpin(e2xd_osum,ints,amps,Wa,Wb,Sa,noA,nvA,noB,nvB,NBasis)
!
! compute omega-like contributions to E2exch-disp (8),
! see Eq. (23) in https://doi.org/10.1063/1.4758455
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB,NBasis
real*8, intent(in)  :: Sa(NBasis,NBasis)
real*8, intent(in)  :: Wa(NBasis,NBasis),Wb(NBasis,NBasis)
real*8, intent(in)  :: ints(nva*noa*nvb*nob)
real*8, intent(in)  :: amps(nva*noa*nvb*nob)

real*8, intent(out) :: e2xd_osum

integer :: i,j
real*8  :: e2xd_oterms(4)

double precision, allocatable :: Saboo(:,:),Sabov(:,:),Sabvo(:,:)
double precision, allocatable :: Sbaoo(:,:),Sabvv(:,:),Sbavo(:,:)
double precision, allocatable :: Wa_ov(:,:),Wb_ov(:,:)

allocate(Sabov(noa,nvb),Sabvo(nva,nob))
allocate(Saboo(noa,nob),Sbaoo(nob,noa))
allocate(Sabvv(nva,nvb))
allocate(Sbavo(nvb,noa))
! Sov
do j=1,nvb
   do i=1,noa
      Sabov(i,j) = Sa(i,nob+j)
   enddo
enddo
! Svo
do j=1,nob
   do i=1,nva
      Sabvo(i,j) = Sa(noa+i,j)
   enddo
enddo
! Soo
do j=1,nob
   do i=1,noa
      Saboo(i,j) = Sa(i,j)
   enddo
enddo
Sbaoo = transpose(Saboo)
! Svv
do j=1,nvb
   do i=1,nva
      Sabvv(i,j) = Sa(noa+i,nob+j)
   enddo
enddo
!Sbavo
do j=1,noa
   do i=1,nvb
      Sbavo(i,j) = Sa(j,nob+i)
   enddo
enddo

allocate(Wa_ov(nob,nvb),Wb_ov(noa,nva))

do j=1,nvb
   do i=1,nob
      Wa_ov(i,j) = Wa(i,nob+j)
   enddo
enddo
do j=1,nva
   do i=1,noa
      Wb_ov(i,j) = Wb(i,noa+j)
   enddo
enddo

e2xd_oterms = 0d0
call oterms15_e2exd(e2xd_oterms(1),amps,Saboo,Sabvo,Sabvv,Wa_ov,Wb_ov,nva,noa,nvb,nob)
call oterms36_e2exd(e2xd_oterms(3),amps,Saboo,Sabov,Sabvv,Wa_ov,Wb_ov,nva,noa,nvb,nob)
call oterm2_e2exd(e2xd_oterms(2),amps,Wa_ov,nva,noa,nvb,nob,Sbaoo,Sabvo,nob)
call oterm4_e2exd(e2xd_oterms(4),amps,Wb_ov,nva,noa,nvb,nob,Saboo,Sbavo,noa)

e2xd_osum = sum(e2xd_oterms)

print*, 'E2exd oterms =', e2xd_osum*1d3

deallocate(Wb_ov,Wa_ov)
deallocate(Saboo,Sbaoo,Sabvv)
deallocate(Sabvo,Sabov,Sbavo)

end subroutine e2exd_oterms_SameSpin

subroutine e2exd_ovterms_OppSpin(e2xd_os,ints,amps,Wa,Wb,Sa,Sb,&
                    noA_a,nvA_a,noA_b,nvA_b,noB_a,nvB_a,noB_b,nvB_b,NBasis)
!
! compute v- and o-like contributions to E2exch-disp (8),
! see Eq. (24) in https://doi.org/10.1063/1.4758455
!
implicit none

integer, intent(in) :: NBasis
integer, intent(in) :: noA_a,nvA_a,noA_b,nvA_b,noB_a,nvB_a,noB_b,nvB_b
real*8, intent(in)  :: Sa(NBasis,NBasis),Sb(NBasis,NBasis)
real*8, intent(in)  :: Wa(NBasis,NBasis),Wb(NBasis,NBasis)
real*8, intent(in)  :: ints(nva_a*noa_a*nvb_b*nob_b)
real*8, intent(in)  :: amps(nva_a*noa_a*nvb_b*nob_b)

real*8, intent(out) :: e2xd_os

integer :: i,j
real*8  :: e2xd_vterms(4),e2xd_oterms(2)
real*8  :: e2xd_vsum,e2xd_osum

double precision, allocatable :: Saboo_aa(:,:),Sabov_aa(:,:),Sabvo_aa(:,:)
double precision, allocatable :: Saboo_bb(:,:),Sabov_bb(:,:),Sabvo_bb(:,:)
double precision, allocatable :: Sbaoo_aa(:,:),Sabvv_aa(:,:),Sbavo_aa(:,:)
double precision, allocatable :: Sbaoo_bb(:,:),Sabvv_bb(:,:),Sbavo_bb(:,:)
double precision, allocatable :: Wa_ov(:,:),Wb_ov(:,:)

allocate(Sabov_aa(noa_a,nvb_a),Sabvo_aa(nva_a,nob_a))
allocate(Saboo_aa(noa_a,nob_a),Sbaoo_aa(nob_a,noa_a))
allocate(Sbavo_aa(nvb_a,noa_a))
! Sov_aa
do j=1,nvb_a
   do i=1,noa_a
      Sabov_aa(i,j) = Sa(i,nob_a+j)
   enddo
enddo
! Svo_aa
do j=1,nob_a
   do i=1,nva_a
      Sabvo_aa(i,j) = Sa(noa_a+i,j)
   enddo
enddo
! Soo_aa
do j=1,nob_a
   do i=1,noa_a
      Saboo_aa(i,j) = Sa(i,j)
   enddo
enddo
Sbaoo_aa = transpose(Saboo_aa)
!Sbavo_aa
do j=1,noa_a
   do i=1,nvb_a
      Sbavo_aa(i,j) = Sa(j,nob_a+i)
   enddo
enddo

allocate(Sabov_bb(noa_b,nvb_b),Sabvo_bb(nva_b,nob_b))
allocate(Saboo_bb(noa_b,nob_b),Sbaoo_bb(nob_b,noa_b))
allocate(Sbavo_bb(nvb_b,noa_b))
! Sov_bb
do j=1,nvb_b
   do i=1,noa_b
      Sabov_bb(i,j) = Sb(i,nob_b+j)
   enddo
enddo
! Svo_bb
do j=1,nob_b
   do i=1,nva_b
      Sabvo_bb(i,j) = Sb(noa_b+i,j)
   enddo
enddo
! Soo_bb
do j=1,nob_b
   do i=1,noa_b
      Saboo_bb(i,j) = Sb(i,j)
   enddo
enddo
Sbaoo_bb = transpose(Saboo_bb)
!Sbavo_bb
do j=1,noa_b
   do i=1,nvb_b
      Sbavo_bb(i,j) = Sb(j,nob_b+i)
   enddo
enddo


allocate(Wa_ov(nob_b,nvb_b),Wb_ov(noa_a,nva_a))

do j=1,nvb_b
   do i=1,nob_b
      Wa_ov(i,j) = Wa(i,nob_b+j)
   enddo
enddo
do j=1,nva_a
   do i=1,noa_a
      Wb_ov(i,j) = Wb(i,noa_a+j)
   enddo
enddo

! v-terms ab
e2xd_vterms = 0d0
call vterm3_e2exd(e2xd_vterms(1),amps,ints,nva_a,noa_a,nvb_b,nob_b,Sbavo_bb,noa_b)
print*, ''
print*, 'e2xd-v os T3 =', e2xd_vterms(1)*1d3

call vterm5_e2exd(e2xd_vterms(2),amps,ints,nva_a,noa_a,nvb_b,nob_b,Sabvo_aa,nob_a)
print*, 'e2xd-v os T5 =', e2xd_vterms(2)*1d3

call vterm7_e2exd(e2xd_vterms(3),amps,ints,nva_a,noa_a,nvb_b,nob_b,Saboo_bb,noa_b)
print*, 'e2xd-v os T7 =', e2xd_vterms(3)*1d3

call vterm8_e2exd(e2xd_vterms(4),amps,ints,nva_a,noa_a,nvb_b,nob_b,Saboo_aa,nob_a)
print*, 'e2xd-v os T8 =', e2xd_vterms(4)*1d3

! o-terms ab
e2xd_oterms = 0d0
call oterm2_e2exd(e2xd_oterms(1),amps,Wa_ov,nva_a,noa_a,nvb_b,nob_b,Sbaoo_aa,Sabvo_aa,nob_a)
call oterm4_e2exd(e2xd_oterms(2),amps,Wb_ov,nva_a,noa_a,nvb_b,nob_b,Saboo_bb,Sbavo_bb,noa_b)
! o-terms checked...
print*, 'e2xd-o os T2 =', e2xd_oterms(1)*1d3
print*, 'e2xd-o os T4 =', e2xd_oterms(2)*1d3

e2xd_vsum = sum(e2xd_vterms)
e2xd_osum = sum(e2xd_oterms)

print*, 'E2exd vterms ab =', e2xd_vsum*1d3
print*, 'E2exd oterms ab =', e2xd_osum*1d3

e2xd_os = e2xd_vsum + e2xd_osum

deallocate(Wb_ov,Wa_ov)
deallocate(Saboo_bb,Sbaoo_bb)
deallocate(Sabvo_bb,Sabov_bb,Sbavo_bb)
deallocate(Saboo_aa,Sbaoo_aa)
deallocate(Sabvo_aa,Sabov_aa,Sbavo_aa)

end subroutine e2exd_ovterms_OppSpin

subroutine load_vovo(intfile,ints,IndNb,nva,noa,nvb,nob)
!
! load (VO|VO) integrals
!
implicit none

character(*) :: intfile
integer, intent(in) :: nva, noa, nvb, nob
integer, intent(in) :: IndNb(2,nob*nvb)
real*8, intent(out) :: ints(nva,noa,nvb,nob)

integer :: iunit
integer :: nova,novb
integer :: ip,iq,ir,is,irs
double precision :: AuxA(noA*nvA)

novA = noA*nvA
novB = noB*nvB

! (OV|OV) (AA|BB) --> (VO|VO)
open(newunit=iunit,file=intfile,status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*novA)

do irs=1,novB

    ir  = IndNB(1,irs)
    is  = IndNB(2,irs)
    read(iunit,rec=is+(ir-noB-1)*noB) AuxA(1:novA)

    do ip=1,nvA
       do iq=1,noA
          ints(ip,iq,ir-noB,is) = AuxA(iq+(ip-1)*noA)
       enddo
    enddo

enddo

close(iunit)

end subroutine load_vovo

subroutine calc_amps(amps,ints,EnA,EnB,nvA,noA,nvB,noB,n)

integer, intent(in) :: noA,nvA,noB,nvB,n
double precision, intent(in) :: EnA(n),EnB(n)
real*8, intent(in)  :: ints(nva,noa,nvb,nob)
real*8, intent(out) :: amps(nva,noa,nvb,nob)

integer :: ip,iq,ir,is
integer :: ipp,irr
real*8 :: dEnA,dEnB
real*8 :: e2du

amps = 0d0

e2du=0d0
do is=1,noB
   do ir=1,nvB
      irr=ir+noB
      dEnB = EnB(irr)-EnB(is)
      do iq=1,noA
         do ip=1,nvA
            ipp=ip+noA
            dEnA = EnA(ipp)-EnA(iq)
            amps(ip,iq,ir,is) = ints(ip,iq,ir,is)/(dEnA+dEnB)
            e2du = e2du + ints(ip,iq,ir,is)**2/(dEnA+dEnB)
         enddo
      enddo
   enddo
enddo
!print*, 'e2du aa =', e2du*1000
write(*,'(/1x,a)') 'from amps...'
print*, 'e2du ab =', e2du*1000

end subroutine calc_amps

subroutine vterm1_e2exd(ene,amps,ints,Sab,nva,noa,nvb,nob)
!
! Term1: +v(a',i,b',j) . S(a',b) . t(a,i,b,j) . S(b'a)
!
integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Sab(nva,nvb)
real*8, intent(in)  :: ints(nva,noa,nvb,nob)
real*8, intent(in)  :: amps(nva,noa,nvb,nob)
real*8, intent(out) :: ene

integer :: i,j
integer :: novo
real*8,allocatable  :: P(:,:,:,:),Q(:,:,:,:)

novo = noa*nvb*nob

allocate(P(nvb,noa,nvb,nob),Q(nvb,nob,nvb,noa))

call dgemm('T','N',nvb,novo,nva,1d0,Sab,nva,ints,nva,0d0,P,nvb)
call tranP_vterm1(P,Q,nvb*noa,nvb*nob)

deallocate(P)

allocate(P(nva,nob,nvb,noa))
call dgemm('N','N',nva,novo,nvb,1d0,Sab,nva,Q,nvb,0d0,P,nva)

do j=1,nob
   do i=1,noa
      ene = ene + sum(P(:,j,:,i)*amps(:,i,:,j))
   enddo
enddo

deallocate(Q,P)

end subroutine vterm1_e2exd

subroutine tranP_vterm1(P,Q,novab,novbb)
!
! transpose P to Q
!
integer, intent(in) :: novab,novbb
real*8, intent(in)  :: P(novab,novbb)
real*8, intent(out) :: Q(novbb,novab)

Q = transpose(P)

end subroutine tranP_vterm1

subroutine vterm2_e2exd(ene,amps,ints,Sab,nva,noa,nvb,nob)
!
! Term2: -v(a,i',b',j) . S(i',b') . t(a,i,b,j) . S(i,b)
!
integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Sab(noa,nvb)
real*8, intent(in)  :: ints(nva,noa,nvb,nob)
real*8, intent(in)  :: amps(nva,noa,nvb,nob)
real*8, intent(out) :: ene

integer :: j,a
real*8  :: P(nva,nob), Q(nva,nob)

P=0d0 ; Q=0d0
do j=1,nob
   do a=1,nva
      P(a,j) = sum(ints(a,:,:,j)*Sab(:,:))
      Q(a,j) = sum(amps(a,:,:,j)*Sab(:,:))
   enddo
enddo

ene = ene - sum(P(:,:)*Q(:,:))

end subroutine vterm2_e2exd

subroutine vterm3_e2exd(ene,amps,ints,nva,noa,nvb,nob,Sba,nsoa)
!
! SameSpin, OppSpin
! Term2: +v(a,i,b',j) . S(b',i') . t(a,i,b,j) . S(b,i')
!
integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsoa
real*8, intent(in)  :: Sba(nvb,nsoa)
real*8, intent(in)  :: ints(nva*noa,nvb,nob)
real*8, intent(in)  :: amps(nva*noa,nvb,nob)
real*8, intent(out) :: ene

integer :: j
integer :: nova
real*8,allocatable :: P(:,:,:),Q(:,:,:)

nova = noa*nva

allocate(P(nva*noa,nsoa,nob),Q(nva*noa,nsoa,nob))

do j=1,nob
   call dgemm('N','N',nova,nsoa,nvb,1d0,ints(:,:,j),nova,Sba,nvb,0d0,P(:,:,j),nova)
   call dgemm('N','N',nova,nsoa,nvb,1d0,amps(:,:,j),nova,Sba,nvb,0d0,Q(:,:,j),nova)
enddo

ene = ene + sum(P(:,:,:)*Q(:,:,:))

deallocate(Q,P)

end subroutine vterm3_e2exd

subroutine vterm4_e2exd(ene,amps,ints,Sab,nva,noa,nvb,nob)
!
! Term 4: -v(a',i,b,j') . S(a',j') . t(a,i,b,j) . S(a,j)
!
integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Sab(nva,nob)
real*8, intent(in)  :: ints(nva,noa,nvb,nob)
real*8, intent(in)  :: amps(nva,noa,nvb,nob)
real*8, intent(out) :: ene

integer :: i,b
real*8  :: P(noa,nvb), Q(noa,nvb)

P=0d0 ; Q=0d0
do i=1,noa
   do b=1,nvb
      P(i,b) = sum(ints(:,i,b,:)*Sab(:,:))
      Q(i,b) = sum(amps(:,i,b,:)*Sab(:,:))
   enddo
enddo

ene = ene - sum(P(:,:)*Q(:,:))

end subroutine vterm4_e2exd

subroutine vterm5_e2exd(ene,amps,ints,nva,noa,nvb,nob,Sab,nsob)
!
! SameSpin, OppSpin
! Term 5: +v(a',i,b,j) . S(a',j') . t(a,i,b,j) . S(a,j')
!
integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsob
real*8, intent(in)  :: Sab(nva,nsob)
real*8, intent(in)  :: ints(nva,noa*nvb*nob)
real*8, intent(in)  :: amps(nva,noa*nvb*nob)
real*8, intent(out) :: ene

integer :: novo
real*8,allocatable :: P(:,:),Q(:,:)

novo = noa*nvb*nob
allocate(P(nsob,novo),Q(nsob,novo))

call dgemm('T','N',nsob,novo,nva,1d0,Sab,nva,ints,nva,0d0,P,nsob)
call dgemm('T','N',nsob,novo,nva,1d0,Sab,nva,amps,nva,0d0,Q,nsob)

ene = ene + sum(P*Q)

deallocate(Q,P)

end subroutine vterm5_e2exd

subroutine vterm6_e2exd(ene,amps,ints,Sba,nva,noa,nvb,nob)
!
! Term6: +v(a,i',b,j') . S(j',i) . t(a,i,b,j) . S(j,i')
!
integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Sba(nob,noa)
real*8, intent(in)  :: ints(nva*noa*nvb,nob)
real*8, intent(in)  :: amps(nva*noa*nvb,nob)
real*8, intent(out) :: ene

integer :: i,j
integer :: nvov
real*8,allocatable  :: P(:,:,:,:),Q(:,:,:,:)

nvov = nva*noa*nvb

allocate(P(nva,noa,nvb,noa),Q(nva,noa,nvb,noa))

call dgemm('N','N',nvov,noa,nob,1d0,ints,nvov,Sba,nob,0d0,P,nvov)
call dgemm('N','N',nvov,noa,nob,1d0,amps,nvov,Sba,nob,0d0,Q,nvov)

do j=1,noa
   do i=1,noa
      ene = ene + sum(P(:,j,:,i)*Q(:,i,:,j))
   enddo
enddo

deallocate(Q,P)

end subroutine vterm6_e2exd

subroutine vterm7_e2exd(ene,amps,ints,nva,noa,nvb,nob,Sab,nsoa)
!
! SameSpin, OppSpin
! Term 7: -v(a,i,b,j') . t(a,i,b,j) . S(i',j') . S(i',j)
!
integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsoa
real*8, intent(in)  :: Sab(nsoa,nob)
real*8, intent(in)  :: ints(nva*noa*nvb,nob)
real*8, intent(in)  :: amps(nva*noa*nvb,nob)
real*8, intent(out) :: ene

integer :: j,jp
real*8  :: P(nob,nob),Q(nob,nob)

P=0d0
do j=1,nob
   do jp=1,nob
      P(jp,j) = sum(ints(:,jp)*amps(:,j))
   enddo
enddo

call dgemm('T','N',nob,nob,nsoa,1d0,Sab,nsoa,Sab,nsoa,0d0,Q,nob)

ene = ene - sum(P*Q)

end subroutine vterm7_e2exd

subroutine vterm8_e2exd(ene,amps,ints,nva,noa,nvb,nob,Sab,nsob)
!
! SameSpin, OppSpin
! Term 8: -v(a,i',b,j) . t(a,i,b,j) . S(i',j') . S(i,j')
!
integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsob
real*8, intent(in)  :: Sab(noa,nsob)
real*8, intent(in)  :: ints(nva,noa,nvb*nob)
real*8, intent(in)  :: amps(nva,noa,nvb*nob)
real*8, intent(out) :: ene

integer :: i,ip
real*8  :: P(noa,noa),Q(noa,noa)

P=0d0
do i=1,noa
   do ip=1,noa
      P(ip,i) = sum(ints(:,ip,:)*amps(:,i,:))
   enddo
enddo

call dgemm('N','T',noa,noa,nsob,1d0,Sab,noa,Sab,noa,0d0,Q,noa)

ene = ene - sum(P*Q)

end subroutine vterm8_e2exd

subroutine oterms15_e2exd(ene,amps,Soo,Svo,Svv,Wa,Wb,nva,noa,nvb,nob)
!
! omega terms 1 and 5
! term 1 : -t(a,i,b,j) . S(a,j) . S(i,j') . wA(j',b)
! term 5 : -t(a,i,b,j) . S(a,j) . wB(i,a'). S(a',b)
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Soo(noa,nob),Svo(nva,nob),Svv(nva,nvb)
real*8, intent(in)  :: Wa(nob,nvb),Wb(noa,nva)
real*8, intent(in)  :: amps(nva,noa,nvb,nob)
real*8, intent(out) :: ene

integer :: i,b
real*8 ::  P(noa,nvb),Q(noa,nvb)
real*8 :: t1, t5

P=0d0
do b=1,nvb
   do i=1,noa
      P(i,b) = P(i,b) + sum(amps(:,i,b,:)*Svo(:,:))
   enddo
enddo

call dgemm('N','N',noa,nvb,nob,1d0,Soo,noa,Wa,nob,0d0,Q,noa)

t1 = t1 - sum(P*Q)
print*, 'Om term 1 =', t1

call dgemm('N','N',noa,nvb,nva,1d0,Wb,noa,Svv,nva,0d0,Q,noa)

t5 = t5 + sum(P*Q)
print*, 'Om term 5 =', t5

ene = ene + t1 + t5

end subroutine oterms15_e2exd

subroutine oterms36_e2exd(ene,amps,Soo,Sov,Svv,Wa,Wb,nva,noa,nvb,nob)
!
! omega terms 3 and 6
! term 3 : -t(a,i,b,j) . S(i,b) . wB(a,i'). S(i',j)
! term 6 : +t(a,i,b,j) . S(i,b) . S(a,b') . wA(b',j)
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB
real*8, intent(in)  :: Soo(noa,nob),Sov(noa,nvb),Svv(nva,nvb)
real*8, intent(in)  :: Wa(nob,nvb),Wb(noa,nva)
real*8, intent(in)  :: amps(nva,noa,nvb,nob)
real*8, intent(out) :: ene

integer :: j,a
real*8 :: WaT(nvb,nob)
real*8 :: P(nva,nob),Q(nva,nob)
real*8 :: t3, t6

P=0d0
do j=1,nob
   do a=1,nva
      P(a,j) = P(a,j) + sum(amps(a,:,:,j)*Sov(:,:))
   enddo
enddo

call dgemm('T','N',nva,nob,noa,1d0,Wb,noa,Soo,noa,0d0,Q,nva)

t3 = - sum(P*Q)
print*, 'Om term 3 =', t3

WaT=transpose(Wa)
call dgemm('N','N',nva,nob,nvb,1d0,Svv,nva,WaT,nvb,0d0,Q,nva)

t6 = sum(P*Q)
print*, 'Om term 6 =', t6

ene = ene + t3 + t6

end subroutine oterms36_e2exd

subroutine oterm2_e2exd(ene,amps,Wa,nva,noa,nvb,nob,Soo,Svo,nsob)
!
! omega term 2 (SameSpin, OppSpin)
! term 2 : t(a,i,b,j) . wA(b,j). S(a,j').S(j',i)
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsoB
real*8, intent(in)  :: Soo(nsob,noa),Svo(nva,nsob)
real*8, intent(in)  :: Wa(nob,nvb)
real*8, intent(in)  :: amps(nva*noa,nvb*nob)
real*8, intent(out) :: ene

integer :: nova,novb
real*8 :: WaT(nvb,nob)
real*8 :: P(nva,noa),Q(nva,noa)

nova = nva*noa
novb = nvb*nob

WaT=transpose(Wa)
call dgemv('N',nova,novb,1d0,amps,nova,WaT,1,0d0,P,1)
call dgemm('N','N',nva,noa,nsob,1d0,Svo,nva,Soo,nsob,0d0,Q,nva)

ene = sum(P*Q)

print*, 'Om Term 2 =', ene

end subroutine oterm2_e2exd

subroutine oterm4_e2exd(ene,amps,Wb,nva,noa,nvb,nob,Soo,Svo,nsoa)
!
! omega term 4
! term 4 : t(a,i,b,j) . wB(a,i). S(i',b).S(i',j)
!
implicit none

integer, intent(in) :: noA,nvA,noB,nvB
integer, intent(in) :: nsoA
real*8, intent(in)  :: Soo(nsoa,nob),Svo(nvb,nsoa)
real*8, intent(in)  :: Wb(noa,nva)
real*8, intent(in)  :: amps(nva*noa,nvb*nob)
real*8, intent(out) :: ene

integer :: nova,novb
real*8 :: WbT(nva,noa)
real*8 :: P(nvb,nob),Q(nvb,nob)

nova = nva*noa
novb = nvb*nob

WbT=transpose(Wb)
call dgemv('T',nova,novb,1d0,amps,nova,WbT,1,0d0,P,1)
call dgemm('N','N',nvb,nob,nsoa,1d0,Svo,nvb,Soo,nsoa,0d0,Q,nvb)

ene = sum(P*Q)

print*, 'Om Term 4 =', ene

end subroutine oterm4_e2exd

end module sapt_open

