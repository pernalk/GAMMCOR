!#define SAPT_OS_DEBUG 10

module sapt_open
use types
use sapt_utils

implicit none

contains

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

call get_den(NBas,A%UMO(:,:,1),A%UOcc(:,1),1d0,PA)
call get_den(NBas,A%UMO(:,:,2),A%UOcc(:,2),1d0,work)
PA = PA + work

call get_den(NBas,B%UMO(:,:,1),B%UOcc(:,1),1d0,PB)
call get_den(NBas,B%UMO(:,:,2),B%UOcc(:,2),1d0,work)
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

subroutine e2ind_o(Flags,A,B,SAPT)
!
! calculate uncoupled and coupled e2ind
! c.f. Eq (17) in https://doi.org/10.1063/1.4758455
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: Nbasis
double precision :: e2iBAa,e2iBAb,e2iABa,e2iABb
double precision :: e2ab_unc,e2ba_unc,e2ind_unc
double precision :: e2ab,e2ba,e2ind
double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)

NBasis = A%NBasis

allocate(Waa(NBasis,NBasis),Wab(NBasis,NBasis))
allocate(Wba(NBasis,NBasis),Wbb(NBasis,NBasis))

call tran2MO(A%WPot,B%UMO(:,:,1),B%UMO(:,:,1),Waa,NBasis)
call tran2MO(A%WPot,B%UMO(:,:,2),B%UMO(:,:,2),Wab,NBasis)

call tran2MO(B%WPot,A%UMO(:,:,1),A%UMO(:,:,1),Wba,NBasis)
call tran2MO(B%WPot,A%UMO(:,:,2),A%UMO(:,:,2),Wbb,NBasis)

!! uncoupled
!e2iBAa = e2ind_unc_o(Wba,A%UOrbE(:,1),A%NOa,A%NVa,NBasis)
!e2iBAb = e2ind_unc_o(Wbb,A%UOrbE(:,2),A%NOb,A%NVb,NBasis)
!!print*, 'e2ind(A<-B)-a = ', e2iBAa*1000
!!print*, 'e2ind(A<-B)-b = ', e2iBAb*1000
!e2ba_unc = e2iBAa + e2iBAb
!
!e2iABa = e2ind_unc_o(Waa,B%UOrbE(:,1),B%NOa,B%NVa,NBasis)
!e2iABb = e2ind_unc_o(Wab,B%UOrbE(:,2),B%NOb,B%NVb,NBasis)
!!print*, 'e2ind(A->B)-a = ', e2iABa*1000
!!print*, 'e2ind(A->B)-b = ', e2iABb*1000
!e2ab_unc = e2iABa + e2iABb
!
!e2ind_unc = e2ba_unc + e2ab_unc
!
!call print_en('Ind(B<--A,unc)',e2ab_unc*1000d0,.true.)
!call print_en('Ind(A<--B,unc)',e2ba_unc*1000d0,.false.)
!call print_en('E2ind(unc)',e2ind_unc*1000d0,.false.)

! uncoupled + coupled
call solve_cpuhf(A,B%WPot,e2ba_unc,e2ba,Flags,NBasis)
call solve_cpuhf(B,A%WPot,e2ab_unc,e2ab,Flags,NBasis)

e2ind_unc = e2ba_unc + e2ab_unc
e2ind     = e2ba + e2ab

SAPT%e2ind_unc = e2ind_unc
SAPT%e2ind = e2ind

call print_en('Ind(B<--A,unc)',e2ab_unc*1000d0,.true.)
call print_en('Ind(A<--B,unc)',e2ba_unc*1000d0,.false.)
call print_en('E2ind(unc)',e2ind_unc*1000d0,.false.)

call print_en('Ind(B<--A)',e2ab*1000d0,.true.)
call print_en('Ind(A<--B)',e2ba*1000d0,.false.)
call print_en('E2ind',e2ind*1000d0,.false.)

deallocate(Wbb,Wba)
deallocate(Wab,Waa)

contains

function e2ind_unc_o(Wmat,OrbE,no,nv,n) result(res)
implicit none

integer :: no,nv,n
double precision :: Wmat(n,n), OrbE(n)
double precision :: delta_e, res

integer :: ip, iq

res = 0d0
do iq=1,no
   do ip=1,nv
      delta_e = OrbE(iq) - OrbE(no+ip)
      res = res + Wmat(iq,no+ip)*Wmat(no+ip,iq)/delta_e
   enddo
enddo

end function e2ind_unc_o

end subroutine e2ind_o

subroutine e2disp_o(Flags,A,B,SAPT)
!
! calculate unrestricted uncoupled / cpld E20disp
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

if (Flags%SaptLevel==0) then
   call e2disp_unc_o(Flags,A,B,SAPT)
else
   call e2disp_cpld_o(Flags,A,B,SAPT)
endif

end subroutine e2disp_o

subroutine e2disp_unc_o(Flags,A,B,SAPT)
!
! calculated unrestricted uncoupled dispersion energy
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBasis
integer :: ANDimX,BNDimX
double precision :: e2daa,e2dbb,e2dab,e2dba
double precision :: e2du

real*8, allocatable :: OmA(:),OmB(:)
real*8, allocatable :: EVecA(:,:),EVecB(:,:)

NBasis = A%NBasis

!print*, 'ISkipped-A =', A%ISkipped
!print*, 'ISkipped-B =', B%ISkipped

! uncoupled
call e2do_unc(e2daa,A%UOrbE(:,1),B%UOrbE(:,1),B%IndNa,A%NOa,A%NVa,B%NOa,B%NVa,NBasis,'OVOVABaa')
call e2do_unc(e2dbb,A%UOrbE(:,2),B%UOrbE(:,2),B%IndNb,A%NOb,A%NVb,B%NOb,B%NVb,NBasis,'OVOVABbb')
call e2do_unc(e2dab,A%UOrbE(:,1),B%UOrbE(:,2),B%IndNb,A%NOa,A%NVa,B%NOb,B%NVb,NBasis,'OVOVABab')
call e2do_unc(e2dba,A%UOrbE(:,2),B%UOrbE(:,1),B%IndNa,A%NOb,A%NVb,B%NOa,B%NVa,NBasis,'OVOVABba')

if(SAPT%IPrint>=10) write(LOUT,'(/1x,a,f16.8)') 'E2disp(unc,aa) = ', e2daa*1000d0
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,ab) = ', e2dab*1000d0
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,ba) = ', e2dba*1000d0
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,bb) = ', e2dbb*1000d0

e2du = e2daa + e2dbb + e2dab + e2dba
SAPT%e2disp_unc = e2du

! summary unc
call print_en('E2disp(unc)',e2du*1000,.false.)

end subroutine e2disp_unc_o

subroutine e2disp_cpld_o(Flags,A,B,SAPT)
!
! calculated unrestricted coupled dispersion energy
!  (uncoupled as a byproduct)
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NBasis
integer :: ANDimX,BNDimX
integer :: i,j
double precision :: e2daa,e2dbb,e2dab,e2dba
double precision :: e2disp_unc
double precision :: e2disp

real*8, allocatable :: OmA(:),OmB(:)
real*8, allocatable :: EVecA(:,:),EVecB(:,:)
real*8, allocatable :: tmp(:,:),ints(:,:)

integer :: ISkipped
double precision,parameter :: BigE   = 1.D8
double precision,parameter :: SmallE = 1.D-3

NBasis = A%NBasis

ANDimX = A%NOVa+A%NOVb
BNDimX = B%NOVa+B%NOVb

allocate(EVecA(ANDimX,ANdimX),OmA(ANDimX))
allocate(ints(ANDimX,BNdimX))
!allocate(ints(ANDimX,BNdimX),tmp(ANDimX,BNDimX))

call readresp(EVecA,OmA,ANDimX,'EIGPRBLA')

ISkipped = A%ISkipped + B%ISkipped
if (ISkipped /= 0) then
   write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
   write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
   write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
endif

!call assemble_spin_ints(A,B,ints,e2disp_unc,NBasis)
! Assemble integrals
allocate(tmp(A%NOVa,B%NOVa))
call load_vovo('OVOVABaa',tmp,B%IndNa,A%NVa,A%NOa,B%NVa,B%NOa)
ints(1:A%NOVa,1:B%NOVa) = tmp(1:A%NOVa,1:B%NOVa)
call e2do_incore_unc(e2daa,tmp,A%UOrbE(:,1),B%UOrbE(:,1),A%NOa,A%NVa,B%NOa,B%NVa)
deallocate(tmp)

allocate(tmp(A%NOVb,B%NOVb))
call load_vovo('OVOVABbb',tmp,B%IndNb,A%NVb,A%NOb,B%NVb,B%NOb)
ints(A%NOVa+1:ANDimX,B%NOVa+1:BNDimX) = tmp(1:A%NOVb,1:B%NOVb)
call e2do_incore_unc(e2dbb,tmp,A%UOrbE(:,2),B%UOrbE(:,2),A%NOb,A%NVb,B%NOb,B%NVb)
deallocate(tmp)

allocate(tmp(A%NOVa,B%NOVb))
call load_vovo('OVOVABab',tmp,B%IndNb,A%NVa,A%NOa,B%NVb,B%NOb)
ints(1:A%NOVa,B%NOVa+1:BNDimX) = tmp(1:A%NOVa,1:B%NOVb)
call e2do_incore_unc(e2dab,tmp,A%UOrbE(:,1),B%UOrbE(:,2),A%NOa,A%NVa,B%NOb,B%NVb)
deallocate(tmp)

allocate(tmp(A%NOVb,B%NOVa))
call load_vovo('OVOVABba',tmp,B%IndNa,A%NVb,A%NOb,B%NVa,B%NOa)
ints(A%NOVa+1:ANDimX,1:B%NOVa) = tmp(1:A%NOVb,1:B%NOVa)
call e2do_incore_unc(e2dba,tmp,A%UOrbE(:,2),B%UOrbE(:,1),A%NOb,A%NVb,B%NOa,B%NVa)
deallocate(tmp)

if(SAPT%IPrint>=10) then 
  write(LOUT,'(/1x,a,f16.8)') 'E2disp(unc,aa) = ', e2daa*1000d0
   write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,ab) = ', e2dab*1000d0
   write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,ba) = ', e2dba*1000d0
   write(LOUT,'(1x,a,f16.8)') 'E2disp(unc,bb) = ', e2dbb*1000d0
endif

e2disp_unc = e2daa + e2dbb + e2dab + e2dba
SAPT%e2disp_unc = e2disp_unc

call print_en('E2disp(unc)',e2disp_unc*1000,.false.)

! coupled 
allocate(tmp(ANDimX,BNDimX))
call dgemm('T','N',ANDimX,BNDimX,ANDimX,1d0,EVecA,ANDimX,ints,ANDimX,0d0,tmp,ANDimX)

deallocate(EVecA)

allocate(EVecB(BNDimX,BNDimX),OmB(BNDimX))

call readresp(EVecB,OmB,BNDimX,'EIGPRBLB')
print*, 'EVecB', norm2(EVecB)

call dgemm('N','N',ANDimX,BNDimX,BNDimX,1d0,tmp,ANDimX,EVecB,BNDimX,0d0,ints,ANDimX)

e2disp=0d0
if (iskipped==0) then

   do j=1,BNDimX
      do i=1,ANDimX
         e2disp = e2disp - ints(i,j)**2d0/(OmA(i)+OmB(j))
      enddo
   enddo

else ! negative/small eigenvalues present

   do j=1,BNDimX
      if(OmB(j).gt.SmallE.and.OmB(j).lt.BigE) then
         do i=1,ANDimX
            if(OmA(i).gt.SmallE.and.OmA(i).lt.BigE) then
               e2disp = e2disp - ints(i,j)**2d0/(OmA(i)+OmB(j))
            endif
         enddo
      endif
   enddo

endif

SAPT%e2disp = e2disp
 
call print_en('E2disp',e2disp*1000,.false.)

deallocate(tmp,ints)
deallocate(OmB,OmA,EVecB)

end subroutine e2disp_cpld_o

subroutine e2do_unc(e2do,EnA,EnB,IndNB,noA,nvA,noB,nvB,n,intfile)

character(*) :: intfile
integer, intent(in) :: noA,nvA,noB,nvB,n
integer, intent(in) :: IndNB(2,noB*nvB)
double precision, intent(in) :: EnA(n),EnB(n)

double precision, intent(out) :: e2do

integer :: iunit
integer :: nova,novb
integer :: ip,iq,ipq,ir,is,irs

double precision :: val,inv_omega
double precision :: OmA,OmB
double precision :: AuxA(noA*nvA)
double precision :: ints(nvA,noA)

novA = noA*nvA
novB = noB*nvB

! (OV|OV) (AA|BB)
open(newunit=iunit,file=intfile,status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*novA)

ints = 0d0
e2do = 0d0
do irs=1,novB

    ir = IndNB(1,irs) ! virt index
    is = IndNB(2,irs) ! occ index
    !ir = noB + mod(irs - 1, nVB) + 1
    !is = (irs - 1)/nVB + 1
    read(iunit,rec=is+(ir-noB-1)*noB) AuxA(1:novA)

    do ip=1,nvA
       do iq=1,noA
          ints(ip,iq) = AuxA(iq+(ip-1)*noA)
       enddo
    enddo

    OmB = EnB(is) - EnB(ir)

    do iq=1,noA
       do ip=1,nvA
          OmA = EnA(iq) - EnA(noA+ip)
          inv_omega = 1d0 / (OmA + OmB)
          e2do = e2do + ints(ip,iq)**2*inv_omega
       enddo
    enddo

enddo
close(iunit)

!print*, 'e2do = ', e2do*1000

end subroutine e2do_unc

subroutine e2do_incore_unc(e2do,ints,EnA,EnB,noA,nvA,noB,nvB)
!
! uncoupled incore
!
integer, intent(in) :: noA,nvA,noB,nvB
double precision, intent(in) :: ints(nva,noa,nvb,nob)
double precision, intent(in) :: EnA(noA+nvA),EnB(noB+nvB)

double precision, intent(out) :: e2do

integer :: ip,iq,ir,is
integer :: ipq,irs
integer :: novA,novB
integer :: iskipped

double precision :: inv_omega
double precision :: OmA,OmB

real(8),allocatable :: OmA0(:),OmB0(:)
logical,allocatable :: condOmA0(:),condOmB0(:)

double precision,parameter :: BigE   = 1.D8
double precision,parameter :: SmallE = 1.D-3

! check for small/negative eigenvals
novA = noA*nvA
novB = noB*nvB
allocate(OmA0(novA),OmB0(novB))
ipq = 0
do iq=1,noA
   do ip=1,nvA
      ipq = ipq + 1
      OmA0(ipq) = EnA(noA+ip) - EnA(iq)
   enddo
enddo
irs = 0
do is=1,noB
   do ir=1,nvB
      irs = irs + 1
      OmB0(irs) = EnB(noB+ir) - EnB(is)
   enddo
enddo
allocate(condOmA0(novA),condOmB0(novB))
condOmA0 = (OmA0.gt.SmallE.and.OmA0.lt.BigE)
condOmB0 = (OmB0.gt.SmallE.and.OmB0.lt.BigE)

iskipped = count(.not. condOmA0) + count(.not. condOmB0)

if (iskipped .gt. 0) then
  print*, 'Skipped ', iskipped, 'value(s) in E2disp(unc)!'
endif

e2do = 0d0
! no negative/small eigenvals
if (iskipped==0) then
   do is=1,noB
      do ir=1,nvB

       OmB = EnB(is) - EnB(noB+ir)

          do iq=1,noA
             do ip=1,nvA
                OmA = EnA(iq) - EnA(noA+ip)
                inv_omega = 1d0 / (OmA + OmB)
                e2do = e2do + ints(ip,iq,ir,is)**2*inv_omega
             enddo
          enddo

      enddo
   enddo

else ! negative/small/huge eigenvalues present

   irs = 0
   do is=1,noB
      do ir=1,nvB
         irs = irs + 1
         if(condOmB0(irs)) then
           ipq = 0
           do iq=1,noA
              do ip=1,nvA
                 ipq = ipq + 1
                 if(condOmA0(ipq)) then
                    inv_omega = 1d0 / (OmA0(ipq) + OmB0(irs))
                    e2do = e2do - ints(ip,iq,ir,is)**2*inv_omega
                 endif
              enddo
           enddo
        endif
      enddo
   enddo

endif

deallocate(OmB0,OmA0)
deallocate(condOmB0,condOmA0)

end subroutine e2do_incore_unc

subroutine e1exchs2_sq_os(A,B,SAPT)
!
! open-shell unrestricted E1exch(S2)
! in second-quanitzed form (o2v2 cost)
! cf. Eq. (13) in https://doi.org/10.1063/1.4758455
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT


integer :: NBasis
integer :: iunit
integer :: i,j,k
integer :: ip,iq,pq,ir,is,rs
double precision :: val,tmp
double precision :: ex1(2),ex2(2),ex3(2),e1exs2
double precision, allocatable :: Sa(:,:),Sb(:,:)
double precision, allocatable :: Sat(:,:),Sbt(:,:)
double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)
double precision, allocatable :: ints(:,:),Aux(:)
double precision, allocatable :: work(:,:)
double precision,external  :: trace

! set dimensions
NBasis = A%NBasis

!! print dimensions
!print*, 'A occ -alpha =', A%NOa
!print*, 'A virt-alpha =', A%NVa
!print*, 'A occ -beta  =', A%NOb
!print*, 'A virt-beta  =', A%NVa
! 
!print*, 'B occ -alpha =', B%NOa
!print*, 'B virt-alpha =', B%NVa
!print*, 'B occ -beta  =', B%NOb
!print*, 'B virt-beta  =', B%NVb

allocate(Waa(NBasis,NBasis),Wab(NBasis,NBasis))
allocate(Wba(NBasis,NBasis),Wbb(NBasis,NBasis))

call tran2MO(A%WPot,B%UMO(:,:,1),B%UMO(:,:,1),Waa,NBasis)
call tran2MO(A%WPot,B%UMO(:,:,2),B%UMO(:,:,2),Wab,NBasis)

call tran2MO(B%WPot,A%UMO(:,:,1),A%UMO(:,:,1),Wba,NBasis)
call tran2MO(B%WPot,A%UMO(:,:,2),A%UMO(:,:,2),Wbb,NBasis)

!print*, 'Waa =', norm2(Waa)
!print*, 'Wab =', norm2(Wab)
!print*, 'Wba =', norm2(Wba)
!print*, 'Wbb =', norm2(Wbb)

allocate(Sa(NBasis,NBasis),Sb(NBasis,NBasis))
allocate(Sat(NBasis,NBasis),Sbt(NBasis,NBasis))
allocate(work(NBasis,NBasis))

call get_one_mat('S',work,A%Monomer,NBasis)
call tran2MO(work,A%UMO(:,:,1),B%UMO(:,:,1),Sa,NBasis)
call tran2MO(work,A%UMO(:,:,2),B%UMO(:,:,2),Sb,NBasis)

Sat = transpose(Sa)
Sbt = transpose(Sb)

! Waa(ja,ba).S(ia,ba).S(ia,ja)
! alpha-alpha
ex1 = 0d0
do k=1,B%NOa
   do j=1,B%NVa
      val = 0d0
      do i=1,A%NOa
         val = val + Sat(B%NOa+j,i)*Sa(i,k)
      enddo
      ex1(1) = ex1(1) + Waa(k,B%NOa+j)*val
   enddo
enddo
!print*, 'wA.Sa.Sa = ',ex1(1)*1000
if(SAPT%IPrint>=10) write(LOUT,'(/,1x,a,f16.8)') 'ExchS2(T1-a ) = ', ex1(1)*1000d0

! Wab(jb,bb).S(ib,jb).S(ib,bb)
! beta-beta
do k=1,B%NOb
   do j=1,B%NVb
      val = 0d0
      do i=1,A%NOb
         val = val + Sbt(B%NOb+j,i)*Sb(i,k)
      enddo
      ex1(2) = ex1(2) + Wab(k,B%NOb+j)*val
   enddo
enddo
!print*, 'wA.Sb.Sb = ',ex1(2)*1000
!print*, 'sum = ',sum(ex1)*1000
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T1-b ) = ', ex1(2)*1000d0

! better: prepare subrtine with a loop for each of these terms
!         a) dgemm for S(i,b).W^t(b,j) = I(i,j)
!         b) loop for \sum_ij I(i,j).S(i,j)

! Wba(ia,aa).S(ia,ja).S(ja,aa)
! alpha-alpha
ex2 = 0d0
do k=1,A%NOa
   do j=1,A%NVa
      val = 0d0
      do i=1,B%NOa
         val = val + Sa(k,i)*Sat(i,A%NOa+j)
      enddo
      ex2(1) = ex2(1) + Wba(k,A%NOa+j)*val
   enddo
enddo
!print*, 'wB.Sa.Sa = ',ex2(1)*1000
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2-a ) = ', ex2(1)*1000d0
! beta-beta
do k=1,A%NOb
   do j=1,A%NVb
      val = 0d0
      do i=1,B%NOb
         val = val + Sb(k,i)*Sbt(i,A%NOb+j)
      enddo
      ex2(2) = ex2(2) + Wbb(k,A%NOb+j)*val
   enddo
enddo
!print*, 'wB.Sb.Sb = ',ex2(2)*1000
!print*, 'sum = ',sum(ex2)*1000
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2-b ) = ', ex2(2)*1000d0

! alpha-alpha (pq|rs).S(q,r).S(p,s)
allocate(Aux(A%NOVa),ints(A%NVa,A%NOa))
open(newunit=iunit,file='OVOVABaa',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*A%NOVa)

ints = 0d0
ex3 = 0d0
do rs=1,B%NOVa

    ir = B%IndNa(1,rs)
    is = B%IndNa(2,rs)
    !print*, 'r,s,rs',ir,is,rs,B%NOVa
    read(iunit,rec=is+(ir-B%NOa-1)*B%NOa) Aux(1:A%NOVa)

    do ip=1,A%NVa
       do iq=1,A%NOa
          ints(ip,iq) = Aux(iq+(ip-1)*A%NOa)
       enddo
    enddo

    val = 0d0
    do ip=1,A%NVa
       do iq=1,A%NOa
          val = val + ints(ip,iq)*Sa(iq,ir)*Sat(is,A%NOa+ip)
       enddo
    enddo

    ex3(1) = ex3(1) + val

enddo

close(iunit)
deallocate(ints,Aux)
!print*, 'v.Sa.Sa = ',ex3(1)*1000
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T3-a ) = ', ex3(1)*1000d0

allocate(Aux(A%NOVb),ints(A%NVb,A%NOb))
open(newunit=iunit,file='OVOVABbb',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*A%NOVb)

do rs=1,B%NOVb

    ir = B%IndNb(1,rs)
    is = B%IndNb(2,rs)
    read(iunit,rec=is+(ir-B%NOb-1)*B%NOb) Aux(1:A%NOVb)

    do ip=1,A%NVb
       do iq=1,A%NOb
          ints(ip,iq) = Aux(iq+(ip-1)*A%NOb)
       enddo
    enddo

    val = 0d0
    do ip=1,A%NVb
       do iq=1,A%NOb
          val = val + ints(ip,iq)*Sb(iq,ir)*Sbt(is,A%NOb+ip)
       enddo
    enddo

    ex3(2) = ex3(2) + val

enddo
close(iunit)
deallocate(ints)

!print*, 'v.Sb.Sb = ',ex3(2)*1000
!print*, 'sum = ', sum(ex3)*1000
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T3-b ) = ', ex3(2)*1000d0

e1exs2 = sum(ex1) + sum(ex2) + sum(ex3)
e1exs2 = - e1exs2
SAPT%exchs2 = e1exs2

call print_en('E1exch(S2)',e1exs2*1000,.true.)

deallocate(Wbb,Wba,Wab,Waa)
deallocate(Sb,Sa,Sbt,Sat)
deallocate(Aux,work)

end subroutine e1exchs2_sq_os

subroutine e1exch_os(A,B,SAPT)
!
! open-shell unrestricted E1exch (S^infty)
! in the AO basis
! cf. Eq. (8) in https://doi.org/10.1063/1.4758455
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO,NMO
integer :: occa,occb
integer :: i,j,k
integer :: info
integer, allocatable :: ipiv(:)

real(8) :: pex(16),pmix
real(8) :: e1ex

real(8), allocatable :: Sa(:,:),Sb(:,:)
real(8), allocatable :: Va(:,:),Vb(:,:)
real(8), allocatable :: taa(:,:),tbb(:,:)
real(8), allocatable :: tab(:,:),tba(:,:)
real(8), allocatable :: taba(:,:),tabb(:,:),&
                        tbaa(:,:),tbab(:,:)
real(8), allocatable :: pa(:,:),pb(:,:)
real(8), allocatable :: Paa(:,:),Pab(:,:),&
                        Pba(:,:),Pbb(:,:)
real(8), allocatable :: pinva(:,:),pinvb(:,:)
real(8), allocatable :: pinvAAa(:,:),pinvAAb(:,:),&
                        pinvBBa(:,:),pinvBBb(:,:),&
                        pinvABa(:,:),pinvABb(:,:),&
                        pinvBAa(:,:),pinvBAb(:,:)
real(8), allocatable :: haa(:,:),hba(:,:),&
                        hab(:,:),hbb(:,:)
real(8), allocatable :: Kaa(:,:),Kab(:,:),&
                        Kba(:,:),Kbb(:,:)
real(8), allocatable :: Jalpha(:,:),Jbeta(:,:)
real(8), allocatable :: Jtaa(:,:),Ktaa(:,:),&
                        Jtba(:,:),Ktba(:,:),&
                        Jtaba(:,:),Ktaba(:,:)
real(8), allocatable :: Jtab(:,:),Ktab(:,:),&
                        Jtbb(:,:),Ktbb(:,:),&
                        Jtabb(:,:),Ktabb(:,:)
real(8), allocatable :: work(:,:),work2(:,:)

double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)
double precision,external  :: trace

! set dimensions
NAO = SAPT%NAO
NMO = A%NBasis

! occuppied A + occupied B
occa = A%NOa+B%NOa
occb = A%NOb+B%NOb

if (NAO.ne.NMO)  then
   print*, 'NAO =', NAO
   print*, 'NMO =', NMO
   stop "e1exh os!"
endif

allocate(PAa(nao,nao),PAb(nao,nao))
allocate(PBa(nao,nao),PBb(nao,nao))

call get_den(nmo,A%UMO(:,:,1),A%Uocc(:,1),1d0,PAa)
call get_den(nmo,A%UMO(:,:,2),A%Uocc(:,2),1d0,PAb)
call get_den(nmo,B%UMO(:,:,1),B%Uocc(:,1),1d0,PBa)
call get_den(nmo,B%UMO(:,:,2),B%Uocc(:,2),1d0,PBb)

!do i=1,nao
!   write(LOUT,*) i
!   write(LOUT,'(10f13.8)') (Paa(i,j),j=1,nao)
!enddo
!write(LOUT,'()')
!print*,'Pbb =',norm2(Pbb)
!do i=1,nao
!   write(LOUT,*) i
!   write(LOUT,'(10f13.8)') (Pbb(i,j),j=1,nao)
!enddo
!write(LOUT,'()')

! h matrices 
allocate(Va(nao,nao),Vb(nao,nao))
call get_one_mat('V',Va,A%Monomer,nao)
call get_one_mat('V',Vb,B%Monomer,nao)

allocate(haa(nao,nao),hab(nao,nao))
allocate(Kaa(nao,nao),Kab(nao,nao))
allocate(hba(nao,nao),hbb(nao,nao))
allocate(Kba(nao,nao),Kbb(nao,nao))
allocate(Jalpha(nao,nao),Jbeta(nao,nao))

call make_J1(nao,PAa,jalpha,'AOTWOSORT')
call make_J1(nao,PAb,jbeta ,'AOTWOSORT')
call make_K(nao,PAa,Kaa,'AOTWOSORT')
call make_K(nao,PAb,Kab,'AOTWOSORT')
haa = Va + Jalpha + JBeta - Kaa
hab = Va + Jalpha + JBeta - Kab

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'hAa =', norm2(hAa)
   print*, 'Kaa =', norm2(Kaa)
   print*, '-----'
   print*, 'hAb =', norm2(hAb)
   print*, 'Kab =', norm2(Kab)
endif

call make_J1(nao,PBa,jalpha,'AOTWOSORT')
call make_J1(nao,PBb,jbeta ,'AOTWOSORT')
call make_K(nao,PBa,Kba,'AOTWOSORT')
call make_K(nao,PBb,Kbb,'AOTWOSORT')
hba = Vb + Jalpha + JBeta - Kba
hbb = Vb + Jalpha + JBeta - Kbb

! save k matrices for exch-ind
allocate(A%ka(nao,nao),A%kb(nao,nao))
allocate(B%ka(nao,nao),B%kb(nao,nao))
A%Ka = Kaa
A%Kb = Kab
B%Ka = Kba
B%Kb = Kbb

! save h matrices for exch-ind
allocate(A%ha(nao,nao),A%hb(nao,nao))
allocate(B%ha(nao,nao),B%hb(nao,nao))
A%ha = haa
A%hb = hab
B%ha = hba
B%hb = hbb

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'hBa =', norm2(hBa)
   print*, 'KBa =', norm2(Kba)
   print*, '-----'
   print*, 'hBb =', norm2(hBb)
   print*, 'KBb =', norm2(Kbb)
   print*, '-----'
endif

allocate(work(NAO,NAO))
allocate(Sa(NMO,NMO),Sb(NMO,NMO))

call get_one_mat('S',work,A%Monomer,NAO)
call tran2MO(work,A%UMO(:,:,1),B%UMO(:,:,1),Sa,NMO)
call tran2MO(work,A%UMO(:,:,2),B%UMO(:,:,2),Sb,NMO)

! P matrix
!    P = [S + 1]^{-1}-1
! where
!         / 1      sab \ |occA 
!   S   = |            | |
!         \ sba     1  / |occB
!          ____________  
!           occA   occB
! P alpha
allocate(pa(occa,occa))
pa=0d0
do concurrent(i=1:occa)
   pa(i,i)=1d0
enddo
do j=1,B%NOa
   do i=1,A%NOa
      pa(i,A%NOa+j)=Sa(i,j)
      pa(A%NOa+j,i)=Sa(i,j)
   enddo
enddo
! P beta
allocate(pb(occb,occb))
pb=0d0
do concurrent(i=1:occb)
   pb(i,i)=1d0
enddo
do j=1,B%NOb
   do i=1,A%NOb
      pb(i,A%NOb+j)=Sb(i,j)
      pb(A%NOb+j,i)=Sb(i,j)
   enddo
enddo
! invert alpha
allocate(ipiv(occa))
allocate(pinva(occa,occa))
allocate(pinvAAa(A%NOa,A%NOa),pinvBBa(B%NOa,B%NOa), &
         pinvABa(A%NOa,B%NOa),pinvBAa(B%NOa,A%NOa))
pinva=0d0
do concurrent(i=1:occa)
   pinva(i,i)=1d0
enddo
call dgesv(occa,occa,pa,occa,ipiv,pinva,occa,info)
! -1 alpha
do i=1,occa
  pinva(i,i)=pinva(i,i)-1d0
enddo
! monomer blocks P alpha
pinvAAa(1:A%NOa,1:A%NOa)=pinva(1:A%NOa,1:A%NOa)
pinvABa(1:A%NOa,1:B%NOa)=pinva(1:A%NOa,A%NOa+1:occa)
pinvBAa(1:B%NOa,1:A%NOa)=pinva(A%NOa+1:occa,1:A%NOa)
pinvBBa(1:B%NOa,1:B%NOa)=pinva(A%NOa+1:occa,A%NOa+1:occa)

deallocate(ipiv)
deallocate(pinva,pa)

! invert beta
allocate(ipiv(occb))
allocate(pinvb(occb,occb))
allocate(pinvAAb(A%NOb,A%NOb),pinvBBb(B%NOb,B%NOb), &
         pinvABb(A%NOb,B%NOb),pinvBAb(B%NOb,A%NOb))
pinvb=0d0
do concurrent(i=1:occb)
   pinvb(i,i)=1d0
enddo
call dgesv(occb,occb,pb,occb,ipiv,pinvb,occb,info)
! -1 beta
do i=1,occb
  pinvb(i,i)=pinvb(i,i)-1d0
enddo
! monomer blocks P beta
pinvAAb(1:A%NOb,1:A%NOb)=pinvb(1:A%NOb,1:A%NOb)
pinvABb(1:A%NOb,1:B%NOb)=pinvb(1:A%NOb,A%NOb+1:occb)
pinvBAb(1:B%NOb,1:A%NOb)=pinvb(A%NOb+1:occb,1:A%NOb)
pinvBBb(1:B%NOb,1:B%NOb)=pinvb(A%NOb+1:occb,A%NOb+1:occb)

deallocate(ipiv)
deallocate(pinvb,pb)

! Tmatrix alpha
allocate(taa(nao,nao),tba(nao,nao),taba(nao,nao),tbaa(nao,nao))

call tran2mo2ao_gen(pinvAAa,A%NOa,A%NOa,NAO,NMO,A%UMO(:,:,1),A%UMO(:,:,1),taa)
call tran2mo2ao_gen(pinvBBa,B%NOa,B%NOa,NAO,NMO,B%UMO(:,:,1),B%UMO(:,:,1),tba)
call tran2mo2ao_gen(pinvABa,A%NOa,B%NOa,NAO,NMO,A%UMO(:,:,1),B%UMO(:,:,1),taba)
call tran2mo2ao_gen(pinvBAa,B%NOa,A%NOa,NAO,NMO,B%UMO(:,:,1),A%UMO(:,:,1),tbaa)

if(SAPT%IPrint>=50) then
   print*, 'Tmatrix alpha'
   print*, 'tAa  = ', norm2(taa)
   print*, 'tBa  = ', norm2(tba)
   print*, 'tABa = ', norm2(taba)
   print*, 'tBAa = ', norm2(tbaa)
endif

deallocate(pinvBBa,pinvBAa,pinvABa,pinvAAa)

! Tmatrix beta
allocate(tab(nao,nao),tbb(nao,nao),tabb(nao,nao),tbab(nao,nao))

call tran2mo2ao_gen(pinvAAb,A%NOb,A%NOb,NAO,NMO,A%UMO(:,:,2),A%UMO(:,:,2),tab)
call tran2mo2ao_gen(pinvBBb,B%NOb,B%NOb,NAO,NMO,B%UMO(:,:,2),B%UMO(:,:,2),tbb)
call tran2mo2ao_gen(pinvABb,A%NOb,B%NOb,NAO,NMO,A%UMO(:,:,2),B%UMO(:,:,2),tabb)
call tran2mo2ao_gen(pinvBAb,B%NOb,A%NOb,NAO,NMO,B%UMO(:,:,2),A%UMO(:,:,2),tbab)

if(SAPT%IPrint>=50) then
   print*, 'Tmatrix beta'
   print*, 'tAb  = ', norm2(tab)
   print*, 'tBb  = ', norm2(tbb)
   print*, 'tABb = ', norm2(tabb)
   print*, 'tBAb = ', norm2(tbab)
endif

deallocate(pinvBBb,pinvBAb,pinvABb,pinvAAb)

! Coulomb exchange alpha
allocate(Jtaa(nao,nao),Ktaa(nao,nao))
allocate(Jtba(nao,nao),Ktba(nao,nao))
allocate(Jtaba(nao,nao),Ktaba(nao,nao))

! prepare make_JK subroutine
! to avoid multiple passes over aotwosort...

call make_J1(nao,taa, Jtaa, 'AOTWOSORT')
call make_J1(nao,tba, Jtba, 'AOTWOSORT')
call make_J1(nao,taba,Jtaba,'AOTWOSORT')

!print*, 'JTABa'
!do j=1,nao
!do i=1,nao
!   if (abs(JTABa(i,j)).gt.1d-6) then
!      write(6,'(1x,2i3,e16.6)') i,j,JTABa(i,j)
!   endif
!enddo
!enddo

call make_K(nao,taa, Ktaa, 'AOTWOSORT')
call make_K(nao,tba, Ktba, 'AOTWOSORT')
call make_K(nao,taba,Ktaba,'AOTWOSORT')

if(SAPT%IPrint>=50) then
   write(6,'(/1x,a)')'Jmat alpha'
   print*, 'Jtaa   =', norm2(Jtaa)
   print*, 'Jtba   =', norm2(Jtba)
   print*, 'Jtaba  =', norm2(Jtaba)
   print*,'Kmat alpha'
   print*, 'Ktaa   =', norm2(Ktaa)
   print*, 'Ktba   =', norm2(Ktba)
   print*, 'Ktaba  =', norm2(Ktaba)
endif


!print*, 'KTABa'
!do j=1,nao
!do i=1,nao
!   if (abs(KTABa(i,j)).gt.1d-6) then
!      write(6,'(1x,2i3,e16.6)') i,j,KTABa(i,j)
!   endif
!enddo
!enddo

! Coulomb exchange alpha
allocate(Jtab(nao,nao),Ktab(nao,nao))
allocate(Jtbb(nao,nao),Ktbb(nao,nao))
allocate(Jtabb(nao,nao),Ktabb(nao,nao))

call make_J1(nao,tab, Jtab, 'AOTWOSORT')
call make_J1(nao,tbb, Jtbb, 'AOTWOSORT')
call make_J1(nao,tabb,Jtabb,'AOTWOSORT')

call make_K(nao,tab, Ktab, 'AOTWOSORT')
call make_K(nao,tbb, Ktbb, 'AOTWOSORT')
call make_K(nao,tabb,Ktabb,'AOTWOSORT')

write(6,'(/1x,a)') 'Jmat beta '
print*, 'Jtab   =', norm2(Jtab)
print*, 'Jtbb   =', norm2(Jtbb)
print*, 'Jtabb  =', norm2(Jtabb)
print*,'Kmat beta '
print*, 'Ktab   =', norm2(Ktab)
print*, 'Ktbb   =', norm2(Ktbb)
print*, 'Ktabb  =', norm2(Ktabb)

! calculations
pex=0d0
!PART(1): PA_alfa*KB_alfa + PA_beta*KB_beta
call dgemm('N','T',nao,nao,nao,1d0,PAa,nao,Kba,nao,0d0,work,nao)
pex(1) = -trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,PAb,nao,Kbb,nao,0d0,work,nao)
pex(1) = pex(1) - trace(work,nao)
!print*, 'part(1) =', pex(1)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (1) = ', pex(1)

!PART(2): TA_alfa*hB_alfa + TA_beta*hB_beta
call dgemm('N','T',nao,nao,nao,1d0,taa,nao,hba,nao,0d0,work,nao)
pex(2) = trace(work,nao)
!print*, 'part(2)a =', pex(2)
call dgemm('N','T',nao,nao,nao,1d0,tab,nao,hbb,nao,0d0,work,nao)
pmix = trace(work,nao)
!print*, 'part(2)b =', pmix
pex(2) = pex(2) + trace(work,nao)
!print*, 'part(2) =', pex(2)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (2) = ', pex(2)

!PART(3): TB_alfa*hA_alfa + TB_beta*hA_beta
call dgemm('N','T',nao,nao,nao,1d0,tba,nao,haa,nao,0d0,work,nao)
pex(3) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tbb,nao,hab,nao,0d0,work,nao)
pex(3) = pex(3) + trace(work,nao)
!print*, 'part(3) =', pex(3)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (3) = ', pex(3)

!PART(4): T^{AB}_alfa*h^A_alfa + T^{AB}_beta*h^A_beta
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,haa,nao,0d0,work,nao)
pex(4) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,hab,nao,0d0,work,nao)
pex(4) = pex(4) + trace(work,nao)
!print*, 'part(4) =', pex(4)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (4) = ', pex(4)

!PART(5): T^{AB}_alfa*h^B_alfa + T^{AB}_beta*h^B_beta
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,hba,nao,0d0,work,nao)
pex(5) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,hbb,nao,0d0,work,nao)
pex(5) = pex(5) + trace(work,nao)
!print*, 'part(5) =', pex(5)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (5) = ', pex(5)

!PART(6): T^{AB}_alfa(J[T^B_alfa]-K[T^B_alfa])
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtba,nao,0d0,work,nao)
pex(6) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktba,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
! beta
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtbb,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tabb,nao,ktbb,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
!print*, 'part(6) =', pex(6)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (6) = ', pex(6)

!L terms (included in part-6)
! T^AB_alpha (J[T^B_beta]) + T^AB_beta (J[T^B_alpha])
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtbb,nao,0d0,work,nao)
pmix = trace(work,nao) 
pex(6) = pex(6) + pmix
!print*, 'p6mix1 =', pmix
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(6) = pex(6) + pmix
!print*, 'p6mix2 =', pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(6)+mix= ', pex(6)

!PART(7): T^{AB}_alfa*(J[T^A_alfa]-K[T^A_alfa]) + T^{AB}_beta*(J[]-K[])
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtaa,nao,0d0,work,nao)
pex(7) = trace(work,nao)
!print*, 'part(7)-ja =', pex(7)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktaa,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
!pex(16) = trace(work,nao)
!print*, 'part(7)-ka =', pex(16)
!print*, 'part(7)-jka =', pex(7)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tabb,nao,jtab,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
call dgemm('T','N',nao,nao,nao,-1d0,tabb,nao,ktab,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
!print*, 'part(7) =', pex(7)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (7) = ', pex(7)

!L2 terms (included in part-7)
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtab,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(7) = pex(7) + pmix
!print*, 'p7mix1 =', pmix
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtaa,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(7) = pex(7) + pmix
!print*, 'p7mix2 =', pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(7)+mix= ', pex(7)

!PART(8): T^A_alfa(J[T^B_alfa]-K[T^B_alfa]) + T^A_beta(J[T^B_beta]-K[T^B_beta])
call dgemm('N','T',nao,nao,nao, 1d0,taa,nao,jtba,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taa,nao,ktba,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
!print*, 'part(8)a =', pex(8)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tab,nao,jtbb,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tab,nao,ktbb,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
!print*, 'part(8)ab=', pex(8)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (8) = ', pex(8)

! mix 1: T^A_alpha . J^B_beta
call dgemm('N','T',nao,nao,nao, 1d0,taa,nao,jtbb,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(8) = pex(8) + pmix
!print*, 'p8mix1 =', pmix
! mix 2: T^A_beta . J^B_alpha
call dgemm('N','T',nao,nao,nao, 1d0,tab,nao,jtba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(8) = pex(8) + pmix
!print*, 'p8mix2 =', pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(8)+mix= ', pex(8)

!PART(9): TAB_a(J[TAB_a]-K[TAB_a])
call dgemm('N','T',nao,nao,nao, 1d0,taba,nao,jtaba,nao,0d0,work,nao)
pex(9) = trace(work,nao)
!print*, 'part(9)-ja  =', pex(9)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktaba,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
!print*, 'part(9)-jka =', pex(9)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tabb,nao,jtabb,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tabb,nao,ktabb,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
!print*, 'part(9)-jkab=', pex(9)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (9) = ', pex(9)

! mix 2*T^AB_beta.J[T^AB_alpha]
call dgemm('N','T',nao,nao,nao,2d0,tabb,nao,jtaba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(9) = pex(9) + trace(work,nao)
!print*, 'p9-mix1=', pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(9)+mix= ', pex(9)

e1ex = sum(pex)
SAPT%e1exch = e1ex
write(6,*) 'E1ex = ', e1ex*1000

call print_en('E1exch',e1ex*1.0d3,.true.)

deallocate(Pab,Paa)
deallocate(Pba,Pbb)
deallocate(Jtab,Jtbb,Jtabb)
deallocate(Ktab,Ktbb,Ktabb)
deallocate(Jtaa,Jtba,Jtaba)
deallocate(Ktaa,Ktba,Ktaba)

deallocate(Sb,Sa)
deallocate(work)

end subroutine e1exch_os



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
! coupled by scaling : 
!    = fac * E2disp/E2disp(unc)
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
real*8 :: fac,e2exd

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

if (SAPT%SaptLevel==0) return

! coupled E2exch-disp by scaling

fac = SAPT%e2disp/SAPT%e2disp_unc
e2exd = e2exd_unc * fac
SAPT%e2exdisp = e2exd

call print_en('Scaling factor',fac,.false.)
call print_en('E2exch-disp',e2exd*1.0d3,.false.)

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

t1=0d0
t1 = t1 - sum(P*Q)
print*, 'Om term 1 =', t1

call dgemm('N','N',noa,nvb,nva,1d0,Wb,noa,Svv,nva,0d0,Q,noa)

t5=0d0
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

t3 = 0d0
t3 = - sum(P*Q)
print*, 'Om term 3 =', t3

WaT=transpose(Wa)
call dgemm('N','N',nva,nob,nvb,1d0,Svv,nva,WaT,nvb,0d0,Q,nva)

t6 = 0d0
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

