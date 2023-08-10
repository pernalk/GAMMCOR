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
integer :: dimOA,dimOB,NBas
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
call get_den(NBas,A%CMO,A%Occ,1d0,work)
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

subroutine e1exch_Chol_dmft(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j, k, l, ia, jb
integer :: ij,ipr
integer :: ip,iq,ir,is
integer :: ipq,iu,it
integer :: iunit
integer :: rdm2type
integer :: dimOA,dimOB,NBas
double precision :: fac,val,nnS2,tmp
double precision :: tmpELST,tmpDEL
double precision :: e1ex_dmft
double precision :: tvk(3),tNa(3),tNb(3),tNaNb(3)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Saa(:,:),Sbb(:,:),Sab(:,:)
double precision,allocatable :: Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: AlphaA(:),AlphaB(:)
double precision,allocatable :: work(:,:),ints(:)

! set dimensions
 NBas = A%NBasis
 dimOA = A%num0+A%num1
 dimOB = B%num0+B%num1

 allocate(AlphaA(dimOA),AlphaB(dimOB))

 rdm2type = Flags%IRdm2Typ
 print*, 'First-order exchange with RDM2 type =',rdm2Type
 select case(rdm2type)
 case(0)
 ! HF
    AlphaA(1:dimOA) = A%Occ(1:dimOA)
    AlphaB(1:dimOB) = B%Occ(1:dimOB)
 case(1)
    ! BB functional
    do i=1,dimOA
       AlphaA(i) = sqrt(A%Occ(i))
    enddo
    do j=1,dimOB
       AlphaB(j) = sqrt(B%Occ(j))
    enddo
 case(2)
   write(LOUT,*) 'POWER FUNCITONAL NOT READY YET!'
   stop
 end select

 allocate(S(NBas,NBas))
 allocate(Sab(NBas,NBas),Saa(NBas,NBas),Sbb(NBas,NBas))
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
 Saa = 0
 Sbb = 0
 do l=1,NBas
    do k=1,NBas
       do i=1,dimOA
          Sbb(k,l) = Sbb(k,l) + A%Occ(i)*Sab(i,k)*Sab(i,l)
       enddo
       do j=1,dimOB
          Saa(k,l) = Saa(k,l) + B%Occ(j)*Sab(k,j)*Sab(l,j)
       enddo
    enddo
 enddo

 deallocate(Vb,Va,S)

 allocate(work(NBas,NBas),ints(NBas**2))

 ! tvk = n_p n_q (v^A S + v^B S + v_pq^qp)
 open(newunit=iunit,file='OOOOABAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOA*dimOB)

 ipq = 0
 tvk = 0
 print*, 'dimOA',dimOA
 print*, 'dimOB',dimOB
 do iq=1,dimOB
    do ip=1,dimOA
       ipq = ipq + 1
       read(iunit,rec=ipq) ints(1:dimOA*dimOB)

       tvk(3) = tvk(3) + A%Occ(ip)*B%Occ(iq)*ints(ip+(iq-1)*dimOA)

    enddo
 enddo
 tvk(3) = -2d0*tvk(3)
 print*, 'tvk(3)',tvk(3)*1000

 close(iunit)

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

 ! tNa
 tNa = 0
 do iq=1,dimOA
    do ip=1,dimOA
       tNa(1) = tNa(1) + AlphaA(ip)*AlphaA(iq)*Saa(ip,iq)*Vbaa(ip,iq)
    enddo
 enddo
 tNa(1) = 2d0*tNa(1)

!(FO|FO): (AA|AB)
open(newunit=iunit,file='OOOOAAAB',status='OLD', &
     access='DIRECT',recl=8*dimOA*dimOA)

! one loop over integrals
ints = 0
val  = 0
do it=1,dimOB
   do iq=1,dimOA
      read(iunit,rec=iq+(it-1)*dimOA) ints(1:dimOA*dimOA)

      fac = A%Occ(iq)*B%Occ(it)*Sab(iq,it)
      val = 0
      do ip=1,dimOA
         val = val + A%Occ(ip)*ints(ip+(ip-1)*dimOA)
      enddo
      tNa(2) = tNA(2) - 4d0*fac*val

      fac = B%Occ(it)*AlphaA(iq)
      val = 0
      do ip=1,dimOA
         val = val + AlphaA(ip)*Sab(ip,it)*ints(ip+(iq-1)*dimOA)
      enddo
      tNa(3) = tNa(3) + 2d0*fac*val

   enddo
enddo

print*, 'tNa-2',tNa(2)*1000
print*, 'tNa-3',tNa(3)*1000

close(iunit)

 tNb = 0
 do iq=1,dimOB
    do ip=1,dimOB
       tNb(1) = tNb(1) + AlphaB(ip)*AlphaB(iq)*Sbb(ip,iq)*Vabb(ip,iq)
    enddo
 enddo
 tNb(1) = 2d0*tNb(1)
 print*, 'tNb-1',tNb(1)*1000

!(FO|FO): (BB|BA)
open(newunit=iunit,file='OOOOBBBA',status='OLD', &
     access='DIRECT',recl=8*dimOB*dimOB)

! one loop over integrals
ints = 0
val  = 0
do it=1,dimOA
   do iq=1,dimOB
      read(iunit,rec=iq+(it-1)*dimOB) ints(1:dimOB*dimOB)

      fac = A%Occ(it)*B%Occ(iq)*Sab(it,iq)
      val = 0
      do ip=1,dimOB
         val = val + B%Occ(ip)*ints(ip+(ip-1)*dimOB)
      enddo
      tNb(2) = tNb(2) - 4d0*fac*val

      fac = A%Occ(it)*AlphaB(iq)
      val = 0
      do ip=1,dimOB
         val = val + AlphaB(ip)*Sab(it,ip)*ints(ip+(iq-1)*dimOB)
      enddo
      tNb(3) = tNb(3) + 2d0*fac*val

   enddo
enddo

print*, 'tNb-2',tNb(2)*1000
print*, 'tNb-3',tNb(3)*1000

close(iunit)

 !open(newunit=iunit,file='OOOOBBAA',status='OLD',&
 !    access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
 open(newunit=iunit,file='OOOOAABB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

 work  = 0
 tmp   = 0
 tNaNb = 0
 !do iq=1,dimOA
 !   do ip=1,dimOA
 !     read(iunit,rec=ip+(iq-1)*dimOA) work(1:dimOB,1:dimOB)

 !     if(ip==iq) then

 !        val = 0
 !        do it=1,dimOB
 !           val = val + B%Occ(it)*work(it,it)
 !        enddo
 !        val = A%Occ(ip)*val
 !        tmpDEL = tmpDEL - 8d0*nnS2*val
 !        tmp = tmp - 8d0*nnS2*val

 !        val = 0
 !        do iu=1,dimOB
 !           do it=1,dimOB
 !              val = val + AlphaB(it)*AlphaB(iu)*Sbb(it,iu)*work(it,iu)
 !           enddo
 !        enddo
 !        tNaNb(1) = tNaNb(1) - 2d0*val*A%Occ(ip)

 !     endif

 !     val = 0
 !     do it=1,dimOB
 !        val = val + B%Occ(it)*work(it,it)
 !     enddo
 !     tNaNb(2) = tNaNb(2) - 2d0*val*AlphaA(ip)*AlphaA(iq)*Saa(ip,iq)

 !     val = 0
 !     do iu=1,dimOB
 !        do it=1,dimOB
 !           val = val + AlphaB(it)*AlphaB(iu)*Sab(ip,iu)*Sab(iq,it)*work(it,iu)
 !        enddo
 !     enddo
 !     tNaNb(3) = tNaNb(3) + AlphaA(ip)*AlphaA(iq)*val

 !   enddo
 !enddo
 do iu=1,dimOB
    do it=1,dimOB
      read(iunit,rec=it+(iu-1)*dimOB) work(1:dimOA,1:dimOA)

      if(it==iu) then

         val = 0d0
         do iq=1,dimOA
            do ip=1,dimOA
               val = val + AlphaA(ip)*AlphaA(iq)*Saa(ip,iq)*work(ip,iq)
            enddo
         enddo
         ! U^A.T^AA
         tNaNb(1) = tNaNb(1) - 2d0*val*B%Occ(it)

      endif

      val = 0d0
      do ip=1,dimOA
         val = val + A%Occ(ip)*work(ip,ip)
      enddo
      tNaNb(2) = tNaNb(2) - 2d0*val*AlphaB(it)*AlphaB(iu)*Sbb(it,iu)

      val = 0d0
      do iq=1,dimOA
         do ip=1,dimOA
            val = val + AlphaA(ip)*AlphaA(iq)*Sab(ip,iu)*Sab(iq,it)*work(ip,iq)
         enddo
      enddo
      ! v.S11.S11
      tNaNb(3) = tNaNb(3) + AlphaB(it)*AlphaB(iu)*val

    enddo
 enddo

 tNaNB = -2d0*tNaNb
 print*, 'A4',sum(tNaNB)*1000
 print*, 'tmp-A4',(tmp+sum(tNaNB))*1000
 print*, 'tNANB-1',tNaNB(1)*1000
 print*, 'tNANB-2',tNaNB(2)*1000
 print*, 'tNANB-3',tNaNB(3)*1000
 !print*, 'tmpDEL',tmpDEL*1000

 close(iunit)

 e1ex_dmft = sum(tvk)+sum(tNa)+sum(tNb)+sum(tNaNb)
 SAPT%exchs2 = e1ex_dmft

 call print_en('E1exch-DMFT(S2)',e1ex_dmft*1000,.true.)

 deallocate(Vbaa,Vabb,Vbba,Vaab)
 deallocate(AlphaB,AlphaA)
 deallocate(ints,work)
 deallocate(Sbb,Saa,Sab)

end subroutine e1exch_Chol_dmft

end module sapt_Chol_exch
