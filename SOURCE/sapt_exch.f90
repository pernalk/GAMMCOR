module sapt_exch
use types
use tran
use exmisc
use exi
!use timing
use sapt_utils

implicit none

contains

subroutine e1exchs2(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j, k, l, ia, jb
integer :: ij,ipr
integer :: ip, iq, ir, is
integer :: iv, iz, iu, it
integer :: iunit
integer :: dimOA,dimOB
integer :: NAO,NBas
double precision,allocatable :: S(:,:),Sab(:,:)
double precision,allocatable :: USa(:,:),USb(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: JJb(:,:)
double precision,allocatable :: Qab(:,:),Qba(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: tmpA(:,:,:,:),tmpB(:,:,:,:), &
                                tmpAB(:,:,:,:)
double precision,allocatable :: work(:,:)
double precision :: tmp,ea,eb,exchs2
double precision :: t1(2),t2a(4),t2b(2),t2c(2),t2d
double precision :: t1f,t2f
!double precision :: Tcpu,Twall
double precision,parameter   :: Half=0.5d0
double precision,external    :: trace
double precision,allocatable :: work1(:)

! set dimensions
 NAO  = SAPT%NAO
 NBas = A%NBasis

 dimOA = A%num0+A%num1
 dimOB = B%num0+B%num1

! call clock('START',Tcpu,Twall)

 allocate(S(NAO,NAO), Sab(NBas,NBas),&
          PA(NAO,NAO),PB(NAO,NAO), &
          Va(NAO,NAO),Vb(NAO,NAO), &
          PAbb(NBas,NBas),PBaa(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas))
 allocate(USa(NAO,NBas),USb(NAO,NBas), &
          Qab(NAO,NBas),Qba(NAO,NBas))
 allocate(work(NBas,NBas),tmp1(NAO,NAO),tmp2(NAO,NAO))

 call get_den(NAO,NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NAO,NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NAO)
 ! OLD
 !call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran_AO2MO2(S,A%CMO,B%CMO,Sab,NAO,NBas)

 call get_one_mat('V',Va,A%Monomer,NAO)
 call get_one_mat('V',Vb,B%Monomer,NAO)

 call tran_AO2MO2(Va,B%CMO,B%CMO,Vabb,NAO,NBas)
 call tran_AO2MO2(Vb,A%CMO,A%CMO,Vbaa,NAO,NBas)
 ! OLD
 !call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 !call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 !allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
 !         RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 !! 2-RDMs
 !RDM2Aval = A%RDM2val
 !RDM2Bval = B%RDM2val

! USa,USb in AOMO
! old
 !call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,A%CMO,NBas,0d0,USa,NBas)
 !call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,B%CMO,NBas,0d0,USb,NBas)
 call dgemm('N','N',NAO,NBas,NAO,1d0,S,NAO,A%CMO,NAO,0d0,USa,NAO)
 call dgemm('N','N',NAO,NBas,NAO,1d0,S,NAO,B%CMO,NAO,0d0,USb,NAO)

! PA(B), PB(A)
 !call tran2MO(PA,USb,USb,PAbb,NBas)
 !call tran2MO(PB,USa,USa,PBaa,NBas)
 call tran_AO2MO2(PA,USb,USb,PAbb,NAO,NBas)
 call tran_AO2MO2(PB,USa,USa,PBaa,NAO,NBas)

! Qab=0; Qba=0
 ! old
 !call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,USb,NBas,0d0,Qab,NBas)
 !call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,USa,NBas,0d0,Qba,NBas)
 !
 call dgemm('N','N',NAO,NBas,NAO,1d0,PA,NAO,USb,NAO,0d0,Qab,NAO)
 call dgemm('N','N',NAO,NBas,NAO,1d0,PB,NAO,USa,NAO,0d0,Qba,NAO)

! old (too large)
 !print*, 'A: num0, num1',A%num0,A%num1
 !print*, 'B: num0, num1',B%num0,B%num1
 !print*, 'dimOA,dimOB',dimOA,dimOB
! call tran3MO_Q(NBas,dimOA,A%CMO,Qba,'TWOA3B')
! call tran3MO_Q(NBas,dimOB,B%CMO,Qab,'TWOB3A')

 if(Flags%ICholesky==0) then 
 ! old
 !call tran4_gen(NBas,&
 !         dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         dimOA,  Qba(1:NBas,1:(A%num0+A%num1)),&
 !         dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         dimOA,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         'TWOA3B','AOTWOSORT')
 !call tran4_gen(NBas,&
 !         dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         dimOB,  Qab(1:NBas,1:(B%num0+B%num1)),&
 !         dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         dimOB,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'TWOB3A','AOTWOSORT')
 !
 call tran4_gen(NAO,&
          dimOA,A%CMO(1:NAO,1:dimOA),&
          dimOA,  Qba(1:NAO,1:dimOA),&
          dimOA,A%CMO(1:NAO,1:dimOA),&
          dimOA,A%CMO(1:NAO,1:dimOA),&
          'TWOA3B','AOTWOSORT')
 call tran4_gen(NAO,&
          dimOB,B%CMO(1:NAO,1:dimOB),&
          dimOB,  Qab(1:NAO,1:dimOB),&
          dimOB,B%CMO(1:NAO,1:dimOB),&
          dimOB,B%CMO(1:NAO,1:dimOB),&
          'TWOB3A','AOTWOSORT')

 call tran4_gen(NAO, &
        dimOA,A%CMO(1:NAO,1:dimOA), &
        dimOA,A%CMO(1:NAO,1:dimOA), &
        dimOB,B%CMO(1:NAO,1:dimOB), &
        dimOB,B%CMO(1:NAO,1:dimOB), &
       'TMPOOAB','AOTWOSORT')
 endif

 !call make_K(NAO,PB,Kb)

! block
! integer :: ip,iq,ir,is
! integer :: NInte1,NInte2,NOcc
! double precision,allocatable :: TwoMO(:)
! double precision :: ETot
! integer,external :: NAddr3
! double precision,external :: FRDM2
!
! !NOcc=A%NAct+A%INAct
! NOcc=A%NAct+A%INAct
! NInte1 = NBas*(NBas+1)/2
! NInte2 = NInte1*(NInte1+1)/2
!
! allocate(TwoMO(NInte2))
!
! call LoadSaptTwoEl(A%Monomer,TwoMO,NBas,NInte2)
! ETot=0
! do ip=1,NOcc
!    do iq=1,NOcc
!      do ir=1,NOcc
!         do is=1,NOcc
!            ETot=ETot+FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)&
!            *TwoMO(NAddr3(ip,ir,iq,is))
!         enddo
!      enddo
!    enddo
! enddo
! print*, 'Check 2-el part of the energy: ',ETot
!
! deallocate(TwoMO)
!
! end block

! T1a
 t1 = 0
 t1(1) = SAPT%elst
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 print*, '???',SAPT%elst
 tmp1=0
 tmp2=0
 !
 !call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 !call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)

 call dgemm('N','N',NAO,NAO,NAO,1d0,PA,NAO,S,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NAO,NAO,1d0,S,NAO,PB,NAO,0d0,tmp2,NAO)
 do j=1,NAO
    do i=1,NAO
       t1(2) = t1(2) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 t1f = 2d0*t1(1)*t1(2)
 !write(LOUT,*) 'T1 ',t1f
 if(SAPT%IPrint>=10) write(LOUT,'(/,1x,a,f16.8)') 'ExchS2(T1   ) = ', t1f*1000d0

! T2d
 t2d = -2d0*SAPT%Vnn*t1(2)
 !write(LOUT,*) 'T2d',t2d
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2d  ) = ', t2d*1000d0

! T2c
 t2c =0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NAO,NAO,NAO,1d0,Va,NAO,PB,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NAO,NAO,1d0,PA,NAO,S,NAO,0d0,tmp2,NAO)
  do j=1,NAO
    do i=1,NAO
       t2c(1) = t2c(1) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo
 t2c(1) = -2d0*t2c(1)
 !write(LOUT,*) 'T2c(1)',t2c(1)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2c-1) = ', t2c(1)*1000d0

! new
 tmp = 0
 do is=1,dimOB
    do iq=1,dimOB
       !t2c(2) = t2c(2) + sum(RDM2Bval(:,:,iq,is)*Vabb(:,:)*PAbb(is,iq))
       tmp = tmp + sum(B%RDM2val(1:dimOB,1:dimOB,iq,is)*Vabb(1:dimOB,1:dimoB)*PAbb(is,iq))
    enddo
 enddo
 t2c(2) = -2d0*tmp
 !write(LOUT,*) 'T2c(2) ',t2c(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2c-2) = ', t2c(2)*1000d0
!
! T2b
 t2b=0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NAO,NAO,NAO,1d0,Vb,NAO,PA,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NAO,NAO,1d0,PB,NAO,S,NAO,0d0,tmp2,NAO)
  do j=1,NAO
    do i=1,NAO
       t2b(1) = t2b(1) + tmp1(i,j)*tmp2(i,j)
    enddo
  enddo
 t2b(1) = -2d0*t2b(1)
 !print*, 'T2b(1)',t2b(1)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2b-1) = ', t2b(1)*1000d0

!new
 t2b(2)=0
 do is=1,dimOA
    do iq=1,dimOA
       !t2b(2) = t2b(2) + sum(RDM2Aval(:,:,iq,is)*Vbaa(:,:)*PBaa(is,iq))
       !t2b(2) = t2b(2) + sum(RDM2Aval(1:dimOA,1:dimOA,iq,is)*Vbaa(1:dimOA,1:dimOA)*PBaa(is,iq))
       t2b(2) = t2b(2) + sum(A%RDM2val(1:dimOA,1:dimOA,iq,is)*Vbaa(1:dimOA,1:dimOA)*PBaa(is,iq))
    enddo
 enddo
 t2b(2) = -2d0*t2b(2)
 !write(LOUT,*) 'T2b(2) ',t2b(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2b-2) = ', t2b(2)*1000d0

 !block
 !double precision :: t2b2a,t2b2b
 !print*, 'Test kasowania 1'
 !t2b2a = 0d0
 !do iq=1,dimOA
 !do ip=1,dimOA
 !   t2b2a = t2b2a + A%Occ(ip)*Vbaa(ip,ip)*A%Occ(iq)*PBaa(iq,iq)
 !enddo
 !enddo
 !t2b2a = -4d0*t2b2a
 !print*, 'T2b-2 a',t2b2a
 !print*, 'Elst*SS',2d0*SAPT%elALL(1)*t1(2)
 !!
 !t2b2b = 0d0
 !do iq=1,dimOA
 !do ip=1,dimOA
 !  t2b2b = t2b2b + A%Occ(ip)*A%Occ(iq)*PBaa(ip,iq)*Vbaa(ip,iq)
 !enddo
 !enddo
 !t2b2b = 2d0*t2b2b
 !print*, 'T2b-2 b',t2b2b*1000
 !print*, 'T2b-2 a+b',(t2b2a+t2b2b)*1000
 !end block

! T2a
 t2a = 0
 do jb=1,NAO
    do ia=1,NAO
       !t2a(1) = t2a(1) + PA(ia,jb)*Kb(jb,ia)
       t2a(1) = t2a(1) + PA(ia,jb)*B%Kmat(jb,ia)
    enddo
 enddo
 t2a(1) = -2.0d0*t2a(1)
 !write(LOUT,*) 'T2a(1)',t2a(1)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-1) = ', t2a(1)*1000d0

 open(newunit=iunit,file='TWOA3B',status='OLD',&
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
     access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

! new
! Qba
 do ir=1,dimOA
    do ip=1,dimOA
       read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOA,1:dimOA)
           do is=1,dimOA
              do iq=1,dimOA
                     t2a(2) = t2a(2) + work(iq,is)* &
                              !RDM2Aval(ip,ir,iq,is)
                              A%RDM2val(ip,ir,iq,is)
              enddo
           enddo
    enddo
 enddo

!! old
!! Qba
! do ir=1,NBas
!    do ip=1,ir
!      read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOA,1:dimOA)
!
!       if(ip==ir) then
!
!         if(ip<=dimOA) then
!            do is=1,dimOA
!               do iq=1,dimOA
!                    t2a(2) = t2a(2) + work(iq,is)* &
!                             RDM2Aval(ip,ir,iq,is)
!                enddo
!             enddo
!         else
!             do iq=1,dimOA
!                  t2a(2) = t2a(2) + work(iq,iq)* &
!                           2d0*A%Occ(ip)*A%Occ(iq)
!             enddo
!         endif
!
!       else
!
!          if(ir<=dimOA) then
!             do is=1,dimOA
!                do iq=1,dimOA
!                     t2a(2) = t2a(2) + work(iq,is)* &
!                            (RDM2Aval(ip,ir,iq,is)+RDM2Aval(ir,ip,iq,is))
!                enddo
!            enddo
!          endif
!
!       endif
!
!    enddo
! enddo

 t2a(2) = -2*t2a(2)
 ! write(LOUT,*) 'T2a(2) ',t2a(2)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-2) = ', t2a(2)*1000d0
 close(iunit,status='DELETE')
 
 !block
 !double precision :: t2a2a,t2a2b
 !print*, 'Test kasowania 2'
 !t2a2a = 0d0
 !do iq=1,dimOB
 !do ip=1,dimOB
 !   t2a2a = t2a2a + B%Occ(ip)*Vabb(ip,ip)*B%Occ(iq)*PAbb(iq,iq)
 !enddo
 !enddo
 !t2a2a = -4d0*t2a2a
 !print*, 'T2a-2 a',t2a2a
 !print*, 'Elst*SS',2d0*SAPT%elALL(2)*t1(2)
 !!
 !t2a2b = 0d0
 !do iq=1,dimOB
 !do ip=1,dimOB
 !  t2a2b = t2a2b + B%Occ(ip)*B%Occ(iq)*PAbb(ip,iq)*Vabb(ip,iq)
 !enddo
 !enddo
 !t2a2b = 2d0*t2a2b
 !print*, 'T2a-2 b',t2a2b*1000
 !print*, 'T2a-2 a+b',(t2a2a+t2a2b)*1000

 !print*, 'Test kasowania 3'
 !print*, 'T2c a', 2d0*SAPT%elALL(3)*t1(2)*1000
 !end block

 open(newunit=iunit,file='TWOB3A',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
! new
! Qab
 tmp = 0
 do ir=1,dimOB
    do ip=1,dimOB
       read(iunit,rec=ip+(ir-1)*dimOB) work(1:dimOB,1:dimOB)
           do is=1,dimOB
              do iq=1,dimOB
                     tmp = tmp + work(iq,is)* &
                           B%RDM2val(ip,ir,iq,is)
              enddo
           enddo
    enddo
 enddo

!! old
!! Qab
! do ir=1,NBas
!    do ip=1,ir
!      read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOB,1:dimOB)
!
!       if(ip==ir) then
!
!         if(ip<=dimOB) then
!            do is=1,dimOB
!               do iq=1,dimOB
!                    t2a(3) = t2a(3) + work(iq,is)* &
!                             RDM2Bval(ip,ir,iq,is)
!                enddo
!             enddo
!         else
!             do iq=1,dimOB
!                  t2a(3) = t2a(3) + work(iq,iq)* &
!                           2d0*B%Occ(ip)*B%Occ(iq)
!             enddo
!         endif
!
!       else
!
!          if(ir<=dimOB) then
!             do is=1,dimOB
!                do iq=1,dimOB
!                     t2a(3) = t2a(3) + work(iq,is)* &
!                              (RDM2Bval(ip,ir,iq,is)+RDM2Bval(ir,ip,iq,is))
!                enddo
!             enddo
!          endif
!
!       endif
!
!    enddo
! enddo

 t2a(3) = -2*tmp
 !write(LOUT,*) 'T2a(3) ',t2a(3)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-3) = ', t2a(3)*1000d0
 close(iunit,status='DELETE')

! T2a(4)

 ! test
 !dimOA = NBas
 !dimOB = NBas
 ! test

 allocate(tmpA(dimOA,dimOA,dimOA,dimOB),tmpB(dimOB,dimOB,dimOA,dimOB),&
          tmpAB(dimOA,dimOA,dimOB,dimOB))
!
! Full NBas check:
!! N^5
! tmpA = 0
! do iz=1,NBas
!    do ir=1,NBas
!       do iq=1,NBas
!          do ip=1,NBas
!             do is=1,NBas
!                 tmpA(ip,iq,ir,iz) = tmpA(ip,iq,ir,iz) + &
!                                  Sab(is,iz)* &
!                                  FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!! N^5
!tmpB = 0
! do iq=1,NBas
!    do iu=1,NBas
!       do iz=1,NBas
!          do iv=1,NBas
!             do it=1,NBas
!                 tmpB(iv,iz,iu,iq) = tmpB(iv,iz,iu,iq) + &
!                                  Sab(iq,it)* &
!                                  FRDM2(iv,iz,iu,it,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
!! N^6
! tmpAB=0
!
! do iu=1,NBas
!    do iv=1,NBas
!       do ir=1,NBas
!          do ip=1,NBas
!             do iz=1,NBas
!                do iq=1,NBas
!                   tmpAB(ip,ir,iv,iu) = tmpAB(ip,ir,iv,iu) + &
!                                        tmpA(ip,iq,ir,iz)*tmpB(iv,iz,iu,iq)
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
! work=0
! open(newunit=iunit,file='TMPMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
!
! do ir=1,NBas
!    do ip=1,NBas
!       read(iunit,rec=ip+(ir-1)*NBas) work
!
!       do iv=1,NBas
!          do iu=1,NBas
!             t2a(4)=t2a(4)+work(iv,iu)* &
!                    tmpAB(ip,ir,iv,iu)
!          enddo
!       enddo
!
!     enddo
! enddo
! t2a(4) = -2*t2a(4)
! print*, 't2a(4): ',t2a(4)

 ! dimOA, dimOB
 ! old
 !tmpA = 0
 !do iz=1,dimOB
 !   do ir=1,dimOA
 !      do iq=1,dimOA
 !         do ip=1,dimOA
 !            do is=1,dimOA
 !                tmpA(ip,iq,ir,iz) = tmpA(ip,iq,ir,iz) + &
 !                                 Sab(is,iz)* &
 !                                 !FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
 !                                 RDM2Aval(ip,ir,iq,is)
 !            enddo
 !         enddo
 !      enddo
 !   enddo
 !enddo
 ! new
 call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%RDM2val,dimOA**3,Sab,NBas,0d0,tmpA,dimOA**3)
 !print*, 'tmpA',norm2(tmpA)

! old:
! N^5
! tmpB = 0
! do iq=1,dimOA
!    do iu=1,dimOB
!       do iz=1,dimOB
!          do iv=1,dimOB
!             do it=1,dimOB
!                 tmpB(iv,iz,iu,iq) = tmpB(iv,iz,iu,iq) + &
!                                  Sab(iq,it)* &
!                                  !FRDM2(iv,iz,iu,it,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
!                                  RDM2Bval(iv,iu,iz,it)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
! new
 do is=1,dimOB
    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%RDM2val(:,:,:,is),dimOB**2,Sab,NBas,0d0,tmpB(:,:,:,is),dimOB**2)
 enddo
 !print*, 'tmpB',norm2(tmpB)

! old:
! N^6
! tmpAB=0
! do iu=1,dimOB
!    do iv=1,dimOB
!       do ir=1,dimOA
!          do ip=1,dimOA
!             do iz=1,dimOB
!                do iq=1,dimOA
!                   tmpAB(ip,ir,iv,iu) = tmpAB(ip,ir,iv,iu) + &
!                                        tmpA(ip,ir,iq,iz)*tmpB(iv,iu,iz,iq)
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
! new
 call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,tmpA,dimOA**2,tmpB,dimOB**2,0d0,tmpAB,dimOA**2)

 !do is=1,dimOB
 !   do iq=1,dimOB
 !      do ir=1,dimOA
 !         do ip=1,dimOA
 !         write(LOUT,'(1x,a,4i2,f12.6)') 'ip,ir,iq,is',ip,ir,iq,is,tmpAB(ip,ir,iq,is)
 !         enddo
 !      enddo
 !   enddo
 !enddo
 !print*, 'tmpAB-exch',norm2(tmpAB)

 ! test
 !dimOA = A%num0+A%num1
 !dimOB = B%num0+B%num1
 ! test

 work=0
 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)

 tmp = 0
 do ir=1,dimOA
    do ip=1,dimOA
       read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOB,1:dimOB)

      tmp = tmp + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))

    enddo
 enddo
 close(iunit)
 t2a(4) = -2*tmp
 ! write(LOUT,*) 't2a(4): ',t2a(4)
 if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T2a-4) = ', t2a(4)*1000d0

 deallocate(tmpAB,tmpB,tmpA)

 exchs2      = t1f+sum(t2a)+sum(t2b)+sum(t2c)+t2d
 SAPT%exchs2 = exchs2
 !write(LOUT,'(/1x,a,f16.8)') 'ExchS2      = ', exchs2*1000d0
 call print_en('ExchS2',exchs2*1000,.true.)

 deallocate(Vbaa,Vabb,PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
 deallocate(Qba,Qab,USb,USa)
 deallocate(tmp2,tmp1,work)

! call clock('E1exch(S2)',Tcpu,Twall)

end subroutine e1exchs2

subroutine e1exchNN_AO(Flags,A,B,SAPT)
!
! noncumulant part of E1exch(S2)
! Eq (21) in T. Korona, JCP 128, 224104 (2008)
! doi: 10.1063/1.2933312
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO, NBasis
integer :: i,j
double precision :: tvk(6), exchNNs2
double precision,allocatable :: XA(:,:),XB(:,:)
double precision,allocatable :: S(:,:),ASB(:,:),BSA(:,:)
double precision,allocatable :: Aux(:,:),work(:,:)

NAO = SAPT%NAO
NBasis = A%NBasis

write(lout,'(/,1x,a)') 'Noncumulant part of E1exch...'

! get S matrix
allocate(S(NAO,NAO))
call get_one_mat('S',S,A%Monomer,NAO)

! get B matrices: X, Ki, Jr
allocate(XB(NAO,NAO),B%Ki(NAO,NAO),B%Jr(NAO,NAO))
XB = 0d0
do i=1,NBasis
   call dger(NAO,NAO,B%Occ(i),B%CMO(:,i),1,B%CMO(:,i),1,XB,NAO)
enddo
XB = 2d0 * XB

! test 1
!SAPT%Vnn = 0d0
!print*, 'setting Vnn = 0!',SAPT%Vnn

call sapt_Ki_AO(B%Ki,XB,SAPT%Vnn,A%XELE,B%XELE,NAO)
call sapt_Jr_AO(B%Jr,XB,SAPT%Vnn,A%XELE,B%XELE,NAO)

! get A matrices: X, Ko, Jl
allocate(XA(NAO,NAO),A%Ko(NAO,NAO),A%Jl(NAO,NAO))
XA = 0d0
do i=1,NBasis
   call dger(NAO,NAO,A%Occ(i),A%CMO(:,i),1,A%CMO(:,i),1,XA,NAO)
enddo
XA = 2d0 * XA

!block
!double precision :: AS(NAO,NAO),AS2(NAO,NAO)
!double precision :: BS(NAO,NAO),BS2(NAO,NAO)
!double precision :: val
!
!call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,AS,NAO)
!call dgemm('N','N',NAO,NAO,NAO,1d0,AS,NAO,AS,NAO,0d0,AS2,NAO)
!AS2 = 0.5d0*AS2
!print*, 'test AS=(AS)^2'
!val = 0d0
!do j=1,NAO
!do i=1,NAO
!   val = abs(AS2(i,j))-abs(AS(i,j))
!   if(val .gt. 1d-6 ) write(6,'(1x,3f12.6)') AS(i,j),AS2(i,j),val
!enddo
!enddo
!
!call dgemm('N','N',NAO,NAO,NAO,1d0,XB,NAO,S,NAO,0d0,BS,NAO)
!call dgemm('N','N',NAO,NAO,NAO,1d0,BS,NAO,BS,NAO,0d0,BS2,NAO)
!BS2 = 0.5d0*BS2
!print*, 'test BS=(BS)^2'
!val = 0d0
!do j=1,NAO
!do i=1,NAO
!   val = abs(BS2(i,j))-abs(BS(i,j))
!   if(val .gt. 1d-6 ) write(6,'(1x,3f12.6)') BS(i,j),BS2(i,j),val
!enddo
!enddo
!
!end block

call sapt_Ko_AO(A%Ko,XA,SAPT%Vnn,A%XELE,B%XELE,NAO)
call sapt_Jl_AO(A%Jl,XA,SAPT%Vnn,A%XELE,B%XELE,NAO)

allocate(work(NAO,NAO),Aux(NAO,NAO))

! T1 = -1/2 A^T.K(B)
call dgemm('T','N',NAO,NAO,NAO,0.5d0,XA,NAO,B%Ki,NAO,0d0,work,NAO)
tvk = 0d0
do i=1,NAO
   tvk(1) = tvk(1) + work(i,i)
enddo
tvk(1) = -tvk(1)
if(SAPT%IPrint>=10) print*, 'Term-1 = ',tvk(1)*1000

allocate(ASB(NAO,NAO),BSA(NAO,NAO))
call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,XB,NAO,0d0,ASB,NAO)

call dgemm('N','N',NAO,NAO,NAO,1d0,XB,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,XA,NAO,0d0,BSA,NAO)

! T2 = (ASB)^T.[Jr(B)-1/2 Ki(B)]
Aux = B%Jr - 0.5d0*B%Ki
print*, 'B%Ki=',norm2(B%Ki)
print*, 'Aux =',norm2(Aux)
call dgemm('N','N',NAO,NAO,NAO,0.5d0,ASB,NAO,Aux,NAO,0d0,work,NAO)
!print*, 'Jr(B)-1/2Ki(B)',norm2(work)
do i=1,NAO
   tvk(2) = tvk(2) + work(i,i)
enddo
tvk(2) = -tvk(2)
if(SAPT%IPrint>=10) print*, 'Term-2 = ',tvk(2)*1000

! T3 = (BSA).[Jl(A)-1/2 Ko(A)]
Aux = 0d0
Aux = A%Jl - 0.5d0*A%Ko
call dgemm('N','N',NAO,NAO,NAO,0.5d0,BSA,NAO,Aux,NAO,0d0,work,NAO)
!print*, 'Jl(A)-1/2Ko(A)',norm2(work)
do i=1,NAO
   tvk(3) = tvk(3) + work(i,i)
enddo
tvk(3) = -tvk(3)
if(SAPT%IPrint>=10) print*, 'Term-3 = ',tvk(3)*1000

! T4 = (ASBSA).Jr(B)
call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,BSA,NAO,0d0,Aux,NAO)
call dgemm('N','N',NAO,NAO,NAO,0.25d0,Aux,NAO,B%Jr,NAO,0d0,work,NAO)
do i=1,NAO
   tvk(4) = tvk(4) + work(i,i)
enddo
if(SAPT%IPrint>=10) print*, 'Term-4 = ',tvk(4)*1000

! T5 = (BSASB).Jl(A)
call dgemm('N','N',NAO,NAO,NAO,1d0,S,NAO,XB,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,BSA,NAO,work,NAO,0d0,Aux,NAO)
call dgemm('N','N',NAO,NAO,NAO,0.25d0,Aux,NAO,A%Jl,NAO,0d0,work,NAO)
do i=1,NAO
   tvk(5) = tvk(5) + work(i,i)
enddo
if(SAPT%IPrint>=10) print*, 'Term-5 = ',tvk(5)*1000

Aux = 0d0
call sapt_Ko_AO(Aux,ASB,SAPT%Vnn,A%XELE,B%XELE,NAO)
call dgemm('N','N',NAO,NAO,NAO,0.125d0,BSA,NAO,Aux,NAO,0d0,work,NAO)
do i=1,NAO
   tvk(6) = tvk(6) + work(i,i)
enddo
tvk(6) = -tvk(6)
if(SAPT%IPrint>=10) print*, 'Term-6 = ',tvk(6)*1000

exchNNs2 = sum(tvk)
SAPT%exchNNs2 = exchNNs2

call print_en('ExchNNS2(Tania)',exchNNs2*1000,.true.)
!print*, 'E1exch(S2)-AO',exchNNs2*1000

deallocate(ASB,BSA)
deallocate(work,Aux)
deallocate(A%Ko,A%Jl)
deallocate(B%Ki,B%Jr)
deallocate(XB,XA,S)

end subroutine e1exchNN_AO

subroutine e1exchNN_AO_noVnn(Flags,A,B,SAPT)
!
! noncumulant part of E1exch(S2)
! version in which we do not calculate Coulomb parts of 2-RDM
! that should cancel with the first term of Eq (9)
! in doi.org/10.1021/acs.jctc.1c00344, i.e.
! 2*E1elst*\sum_pq np nq Spq^2
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO, NBasis
integer :: i,j
double precision :: term2(3),term3(2),term4(2),term5(3)
double precision :: tvk(6), exchNNs2
double precision,allocatable :: XA(:,:),XB(:,:)
double precision,allocatable :: JA(:,:),KA(:,:),JB(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: S(:,:),ASB(:,:),BSA(:,:)
double precision,allocatable :: Aux(:,:),work(:,:)

NAO = SAPT%NAO
NBasis = A%NBasis

write(lout,'(/,1x,a)') 'Try noncumulant part of E1exch...'

! get S matrix
allocate(S(NAO,NAO),Va(NAO,NAO),Vb(NAO,NAO))
call get_one_mat('S',S,A%Monomer,NAO)
call get_one_mat('V',Va,A%Monomer,NAO)
call get_one_mat('V',Vb,B%Monomer,NAO)

! get B matrices: X
allocate(XB(NAO,NAO),XA(NAO,NAO))
allocate(Aux(NAO,NAO),work(NAO,NAO))
XB = 0d0
do i=1,NBasis
   call dger(NAO,NAO,B%Occ(i),B%CMO(:,i),1,B%CMO(:,i),1,XB,NAO)
enddo
XB = 2d0 * XB

! get A matrices: X, Ko, Jl
XA = 0d0
do i=1,NBasis
   call dger(NAO,NAO,A%Occ(i),A%CMO(:,i),1,A%CMO(:,i),1,XA,NAO)
enddo
XA = 2d0 * XA

allocate(ASB(NAO,NAO),BSA(NAO,NAO))
call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,XB,NAO,0d0,ASB,NAO)

call dgemm('N','N',NAO,NAO,NAO,1d0,XB,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,XA,NAO,0d0,BSA,NAO)

! start with the terms!
tvk = 0d0

! Coulomb part, gets cancelled by Elst_Vb*SS
!call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,Vb,NAO,0d0,work,NAO)
!trXV = 0d0
!do i=1,NAO
!   trXV = trXV + work(i,i)
!enddo
!Aux = trXV * S
!call dgemm('N','N',NAO,NAO,NAO,0.5d0,ASB,NAO,Aux,NAO,0d0,work,NAO)
!do i=1,NAO
!   tvk(2) = tvk(2) + work(i,i)
!enddo
!tvk(2) = -tvk(2)
!print*, 'Term-2 = ',tvk(2)*1000

term2 = 0d0
term3 = 0d0
term4 = 0d0
term5 = 0d0

! T1 = -1/2 ASB.va 
do j=1,NAO
   do i=1,NAO
      term2(1) = term2(1) + ASB(i,j)*Va(i,j)
   enddo
enddo
term2(1) = -0.5d0*term2(1)
print*, 'Term-2a = ',term2(1)*1d3

do j=1,NAO
   do i=1,NAO
      term2(2) = term2(2) + BSA(i,j)*Vb(i,j)
   enddo
enddo
term2(2) = -0.5d0*term2(2)
print*, 'Term-2b = ',term2(2)*1d3

! T1 = -1/2 A^T.K(B)
call dgemm('T','N',NAO,NAO,NAO,1d0,XA,NAO,B%Kmat,NAO,0d0,work,NAO)
do i=1,NAO
   term2(3) = term2(3) + work(i,i)
enddo
term2(3) = -term2(3)
print*, 'Term-2c = ',term2(3)*1d3

! Term3, exchange part : n^B . Gamma^A . v^B . S
call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,Aux,NAO)
call dgemm('N','T',NAO,NAO,NAO,1d0,ASB,NAO,Aux,NAO,0d0,work,NAO)
do j=1,NAO
   do i=1,NAO
      term3(1) = term3(1) + work(i,j)*Vb(i,j)
   enddo
enddo
term3(1) = 0.25d0*term3(1)
print*, 'Term-3a = ',term3(1)*1000

! Term4, exchange part : n^A . Gamma^B . v^A . S
call dgemm('N','N',NAO,NAO,NAO,1d0,XB,NAO,S,NAO,0d0,Aux,NAO)
call dgemm('N','T',NAO,NAO,NAO,1d0,BSA,NAO,Aux,NAO,0d0,work,NAO)
do j=1,NAO
do i=1,NAO
   term4(1) = term4(1) + work(i,j)*Va(i,j)
enddo
enddo
term4(1) = 0.25d0*term4(1)
print*, 'Term-4a = ',term4(1)*1000

allocate(JA(NAO,NAO),JB(NAO,NAO),KA(NAO,NAO))
call make_J1(NAO,XB,JB,'AOTWOSORT')
call make_J1(NAO,XA,JA,'AOTWOSORT')
call make_K(NAO,XA,KA)

! B%Kmat contains 1/2!
Aux = JB - B%Kmat
call dgemm('N','N',NAO,NAO,NAO,1d0,ASB,NAO,Aux,NAO,0d0,work,NAO)
do i=1,NAO
   term3(2) = term3(2) + work(i,i)
enddo
term3(2) = -0.5d0*term3(2)
print*, 'Term-3b = ',term3(2)*1000

Aux = JA - 0.5d0*KA
call dgemm('N','N',NAO,NAO,NAO,1d0,BSA,NAO,Aux,NAO,0d0,work,NAO)
do i=1,NAO
   term4(2) = term4(2) + work(i,i)
enddo
term4(2) = -0.5d0*term4(2)
print*, 'Term-4b = ',term4(2)*1000

! T4 = (ASBSA).J(B)
call dgemm('N','N',NAO,NAO,NAO,1d0,XA,NAO,S,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,work,NAO,BSA,NAO,0d0,Aux,NAO)
call dgemm('N','N',NAO,NAO,NAO,0.25d0,Aux,NAO,JB,NAO,0d0,work,NAO)
do i=1,NAO
   term5(1) = term5(1) + work(i,i)
enddo
print*, 'Term-5a = ',term5(1)*1d3

! T5 = (BSASB).J(A)
call dgemm('N','N',NAO,NAO,NAO,1d0,S,NAO,XB,NAO,0d0,work,NAO)
call dgemm('N','N',NAO,NAO,NAO,1d0,BSA,NAO,work,NAO,0d0,Aux,NAO)
call dgemm('N','N',NAO,NAO,NAO,0.25d0,Aux,NAO,JA,NAO,0d0,work,NAO)
do i=1,NAO
   term5(2) = term5(2) + work(i,i)
enddo
print*, 'Term-5b = ',term5(2)*1d3

call make_K(NAO,ASB,Aux)
call dgemm('N','N',NAO,NAO,NAO,0.125d0,BSA,NAO,Aux,NAO,0d0,work,NAO)
do i=1,NAO
   term5(3) = term5(3) + work(i,i)
enddo
term5(3) = -term5(3)
print*, 'Term-5c = ',term5(3)*1d3

exchNNs2 = sum(term2)+sum(term3)+sum(term4)+sum(term5)

call print_en('ExchNNS2',exchNNs2*1000,.true.)
print*, ''
!print*, 'E1exch-nn =',exchNNs2*1d3

deallocate(JB)
deallocate(work,Aux)
deallocate(BSA,ASB)
deallocate(XB,XA)

end subroutine e1exchNN_AO_noVnn

subroutine e1exch_NaNb(Flags,A,B,SAPT)
!
! E1exch(S2): Eq (9) in SAPT(MC) paper
! doi: 10.1021/acs.jctc.1c00344
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
integer :: NAO,NBas
integer :: dimOA,dimOB
double precision :: fac,val,nnS2,tmp
double precision :: tElst,tvk(3),tNa(2),tNb(2),tNaNb
double precision :: exchs2
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Sab(:,:),Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: intA(:,:,:,:),intB(:,:,:,:)
double precision,allocatable :: tmpAB(:,:,:,:)
double precision,allocatable :: work(:,:),ints(:)
double precision,external  :: ddot

print*, 'Testing E1exch NaNb...'

! set dimensions
NAO  = SAPT%NAO
NBas = A%NBasis

!  for TREXIO!
!dimOA = NBas
!dimOB = NBas
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

allocate(S(NAO,NAO),Sab(NBas,NBas))
allocate(Va(NAO,NAO),Vb(NAO,NAO), &
         Vabb(NBas,NBas),Vbaa(NBas,NBas),&
         Vaab(NBas,NBas),Vbba(NBas,NBas))

call get_one_mat('V',Va,A%Monomer,NAO)
call get_one_mat('V',Vb,B%Monomer,NAO)

!call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
!call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)
!call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
!call tran2MO(Vb,B%CMO,A%CMO,Vbba,NBas)

call tran_AO2MO2(Va,B%CMO,B%CMO,Vabb,NAO,NBas)
call tran_AO2MO2(Vb,A%CMO,A%CMO,Vbaa,NAO,NBas)
call tran_AO2MO2(Va,A%CMO,B%CMO,Vaab,NAO,NBas)
call tran_AO2MO2(Vb,B%CMO,A%CMO,Vbba,NAO,NBas)

call get_one_mat('S',S,A%Monomer,NAO)
!call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
call tran_AO2MO2(S,A%CMO,B%CMO,Sab,NAO,NBas)

allocate(intA(dimOA,dimOA,dimOA,dimOB),&
         intB(dimOB,dimOB,dimOA,dimOB))


call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,A%RDM2val,dimOA**3,Sab,NBas,0d0,intA,dimOA**3)
!call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2Bval,dimOB**3,Sab,NBas,0d0,intB,dimOB**3)

! careful! intB(B,B,A,B)
do is=1,dimOB
   call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,B%RDM2val(:,:,:,is),dimOB**2,Sab,NBas,0d0,intB(:,:,:,is),dimOB**2)
enddo

!deallocate(RDM2Bval,RDM2Aval)
deallocate(Vb,Va,S)

! n^A n^B Sab Sab
nnS2 = 0
do j=1,dimOB
do i=1,dimOA
   nnS2 = nnS2 + A%Occ(i)*B%Occ(j)*Sab(i,j)**2
enddo
enddo

!print*, 'Sab',norm2(Sab)
!do i=1,NBas
!   write(lout,'(i3,2f12.6)') i,A%Occ(i),B%Occ(i)
!enddo

tElst = 2d0*(SAPT%elst-SAPT%Vnn)*nnS2
print*, 'tELST',tELST*1000

allocate(ints(NBas**2),work(NAO,NAO))

! tvk = n_p n_q (v^A S + v^B S + v_pq^qp)
!open(newunit=iunit,file='FFOOABAB',status='OLD',&
!    access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
!ipq = 0
!tvk = 0
!do iq=1,dimOB
!   do ip=1,dimOA
!      ipq = ipq + 1
!      read(iunit,rec=ipq) ints(1:NBas*NBas)
!
!      tvk(3) = tvk(3) + A%Occ(ip)*B%Occ(iq)*ints(ip+(iq-1)*NBas)
!
!   enddo
!enddo
!tvk(3) = -2d0*tvk(3)
!print*, 'tvk(3)',tvk(3)*1000

!open(newunit=iunit,file='FOFOABBA',status='OLD',&
!    access='DIRECT',form='UNFORMATTED',recl=8*NBas*dimOB)
!tvk = 0
!do iq=1,dimOA
!   do ip=1,dimOB
!      ipq = ipq + 1
!      read(iunit,rec=(ip+(iq-1)*NBas)) ints(1:NBas*dimOB)
!
!      tvk(3) = tvk(3) + A%Occ(iq)*B%Occ(ip)*ints(iq+(ip-1)*NBas)
!
!   enddo
!enddo
!tvk(3) = -2d0*tvk(3)
!print*, 'tvk(3)',tvk(3)*1000
!print*, 'HERE: sth wrong with (pq|qp)?'
!print*, 'HERE: use Kb instead??'

!close(iunit)

!! 3rd test...!
!ipq = 0
!tvk = 0
!do iq=1,dimOB
!   do ip=1,dimOA
!      i=(iq-1)*NBas+ip
!     tvk(3) = tvk(3) + A%Occ(ip)*B%Occ(iq)*ddot(A%NChol,A%FFAB(:,i),1,A%FFAB(:,i),1)
!   enddo
!enddo
!tvk(3) = -2d0*tvk(3)
!print*, 'tvk(3)',tvk(3)*1000

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
do jb=1,NAO 
   do ia=1,NAO
      tvk(3) = tvk(3) + work(ia,jb)*B%Kmat(jb,ia)
   enddo
enddo
tvk(3) = -2.0d0*tvk(3)
print*, 'tvk(3)',tvk(3)*1000

!print*, 'B%Kmat',norm2(B%Kmat)

deallocate(work)
allocate(work(NBas,NBas))

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

!(FO|FO): (AA|AB)
open(newunit=iunit,file='FOFOAAAB',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
do it=1,dimOB
   do iq=1,dimOA
      read(iunit,rec=iq+(it-1)*NBas) ints(1:NBas*dimOA)

      do ir=1,dimOA
         do ip=1,dimOA
            tNa(2) = tNa(2) + B%Occ(it)*intA(ip,ir,iq,it)*ints(ip+(ir-1)*NBas)
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

!(FO|FO): (BB|BA)
open(newunit=iunit,file='FOFOBBBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
ints = 0
do it=1,dimOA
   do iq=1,dimOB
      read(iunit,rec=iq+(it-1)*NBas) ints(1:NBas*dimOB)

      do ir=1,dimOB
         do ip=1,dimOB
            tNb(2) = tNb(2) + A%Occ(it)*intB(ip,ir,it,iq)*ints(ip+(ir-1)*NBas)
         enddo
      enddo

   enddo
enddo
tNb(2) = -2d0*tNb(2)
print*, 'tNb-1',tNb(1)*1000
print*, 'tNb-2',tNb(2)*1000

close(iunit)

!open(newunit=iunit,file='TMPOOAB',status='OLD',&
!    access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
open(newunit=iunit,file='FOFOAABB',status='OLD',&
    access='DIRECT',form='UNFORMATTED',recl=8*dimOA*NBas)

allocate(tmpAB(dimOA,dimOA,dimOB,dimOB))

call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,intA,dimOA**2,intB,dimOB**2,0d0,tmpAB,dimOA**2)

val  = 0
work = 0
!do ir=1,dimOA
!   do ip=1,dimOA
    ! read(iunit,rec=ip+(ir-1)*dimOA) work(1:dimOB,1:dimOB)
     !val = val + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))
     !intA(ip,ir,iq,it)*intB(it,iu,is,iq)*ints(it,is)
!   enddo
!enddo

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

end subroutine e1exch_NaNb

!subroutine hl_2el(Flags,A,B,SAPT)
!! calculate Heitler-London energy
!! for 2-electron CAS monomers
!implicit none
!
!type(FlagsData)   :: Flags
!type(SaptData)    :: SAPT
!type(SystemBlock) :: A, B
!
!integer :: i,j,ij,k,l,kl,ia,jb
!integer :: ip,iq,ir,is
!integer :: iunit
!integer :: dimOA,dimOB,NBas
!double precision :: val,fac,tmp
!double precision :: ccs2,nns2,ccaaaa,ccbbbb
!double precision :: t12(3),t34(3),t13(6)
!double precision :: Dfull,Ds2inv,e1,e2,ehl,ehlint
!double precision :: e1s2,e2s2,ehls2,ehls2int
!double precision :: deltaMs2
!double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
!double precision,allocatable :: PA(:,:),PB(:,:), &
!                                PAbb(:,:),PBaa(:,:)
!double precision,allocatable :: Va(:,:),Vb(:,:), &
!                                Ha(:,:),Haa(:,:),&
!                                HB(:,:),Hbb(:,:),Hab(:,:)
!double precision,allocatable :: Kb(:,:),Jmat(:,:)
!double precision,allocatable :: work(:,:),work1(:),ints(:,:)
!double precision,allocatable :: RDM2Aval(:,:,:,:), &
!                                RDM2Bval(:,:,:,:)
!double precision,parameter   :: Half=0.5d0
!double precision,external    :: trace
!! test
!double precision,allocatable :: tmpAB(:,:,:,:),Sab_save(:,:)
!double precision,allocatable :: work2(:),work3(:)
!
!! info
! write(LOUT,'(/,1x,a,/)') 'HEITLER-LONDON ENERGY FOR 2-el MONOMERS'
!
!! set dimensions
! NBas = A%NBasis
! dimOA = A%num0+A%num1
! dimOB = B%num0+B%num1
!
! allocate(S(NBas,NBas),Sab(NBas,NBas),&
!          PA(NBas,NBas),PB(NBas,NBas),&
!          Va(NBas,NBas),Vb(NBas,NBas),&
!          PAbb(NBas,NBas),PBaa(NBas,NBas))
! allocate(work(NBas,NBas),work1(NBas*NBas),&
!          Ha(NBas,NBas),Haa(NBas,NBas),&
!          Hb(NBas,NBas),Hbb(NBas,NBas),&
!          Hab(Nbas,NBas))
! allocate(Sba(NBas,NBas))
!
! call get_den(NBas,A%CMO,A%Occ,1d0,PA)
! call get_den(NBas,B%CMO,B%Occ,1d0,PB)
!
! call get_one_mat('S',S,A%Monomer,NBas)
! call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
! call tran2MO(S,B%CMO,A%CMO,Sba,NBas)
!
! call get_one_mat('V',Va,A%Monomer,NBas)
! call get_one_mat('V',Vb,B%Monomer,NBas)
!
! call get_one_mat('H',Ha,A%Monomer,NBas)
! call get_one_mat('H',Hb,B%Monomer,NBas)
!
! work = 0
! work = Ha + Vb
! !work = Vb
! call tran2MO(work,A%CMO,A%CMO,Haa,NBas)
!
! work = 0
! !work = Va
! work = Hb + Va
! call tran2MO(work,B%CMO,B%CMO,Hbb,NBas)
!
! call tran2MO(work,A%CMO,B%CMO,Hab,NBas)
!
! allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
!          RDM2Bval(dimOB,dimOB,dimOB,dimOB))
!
! ! CAS
! RDM2Aval = A%RDM2val
! RDM2Bval = B%RDM2val
!
! ! denominator
! ccs2 = 0
! nns2 = 0
! do iq=1,NBas
!    do ip=1,NBas
!       nns2 = nns2 + A%Occ(ip)*B%Occ(iq)*Sab(ip,iq)*Sab(ip,iq)
!       ccs2 = ccs2 + A%CICoef(ip)*B%CICoef(iq)*Sab(ip,iq)*Sab(ip,iq)
!    enddo
! enddo
! Dfull  = 1d0 - 2d0 * nns2 + ccs2**2
! Ds2inv = 1d0 + 2d0 * nns2
!
! print*, 'D      ', Dfull
! print*, '1/D(S2)', Ds2inv
! print*, ''
!
! ! <Psi0|HA|Psi0>
! ! 1-el part
! e1   = 0
! e1s2 = 0
! val  = 0
! do ip=1,NBas
!    val = val + A%Occ(ip)*Haa(ip,ip) + B%Occ(ip)*Hbb(ip,ip)
! enddo
! e1 = 2d0*val
!
!! extra terms in 1-el S2
! e1s2 = 4d0*nns2*val
! !!! TEST - elst
! !e1s2 = e1
!
! val = 0
! do iq=1,NBas
!    do ip=1,NBas
!       val = val + A%Occ(ip)*B%Occ(iq)*Hab(ip,iq)*Sab(ip,iq)
!    enddo
! enddo
! e1 = e1 - 4d0*val
!
! val = 0
! do ir=1,NBas
!    do ip=1,NBas
!
!       fac = 2d0*A%CICoef(ip)*A%CICoef(ir)*Haa(ip,ir)
!
!       val = 0
!       do iq=1,NBas
!          val = val + B%Occ(iq)*Sab(ip,iq)*Sab(ir,iq)
!       enddo
!       e1 = e1 - fac*val
!
!    enddo
! enddo
! !print*, 'E1-3',e1
!
! val = 0
! do ir=1,NBas
!    do ip=1,NBas
!
!       fac = 2d0*B%CICoef(ip)*B%CICoef(ir)*Hbb(ip,ir)
!
!       val = 0
!       do iq=1,NBas
!          val = val + A%Occ(iq)*Sab(iq,ip)*Sab(iq,ir)
!       enddo
!       e1 = e1 - fac*val
!
!    enddo
! enddo
!
! ! S2 approx
! e1s2 = e1s2 + e1
!
! val = 0
! do iq=1,NBas
!    do ip=1,NBas
!       val = val + A%CICoef(ip)*B%CICoef(iq)*Hab(ip,iq)*Sab(ip,iq)
!    enddo
! enddo
! e1 = e1 + 4d0*ccs2*val
! print*, '1-el part',e1
!
! ! two-electron part
!
! ! test Coulomb and exchange
! allocate(Kb(NBas,NBas),Jmat(NBas,NBas))
! call make_K(NBas,PB,Kb)
! call make_J1(NBas,PB,Jmat,'AOTWOSORT')
!
! tmp = 0
! do jb=1,NBas
!    do ia=1,NBas
!       tmp = tmp + PA(ia,jb)*Jmat(jb,ia)
!    enddo
! enddo
! tmp = 4.0d0*tmp
! write(LOUT,*) 'Coulomb ',tmp
!
! tmp = 0
! do jb=1,NBas
!    do ia=1,NBas
!       tmp = tmp + PA(ia,jb)*Kb(jb,ia)
!    enddo
! enddo
! tmp = -2.0d0*tmp
! write(LOUT,*) 'Exchange',tmp
! deallocate(Jmat,Kb)
! ! end test
!
!! Coulomb and exchange
! call tran4_gen(NBas,&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          'OOOOAABB','AOTWOSORT')
! call tran4_gen(NBas,&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          'OOOOABAB','AOTWOSORT')
!
! !call tran4_gen(NBas,&
! !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
! !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
! !         NBas,A%CMO,&
! !         NBas,B%CMO,&
! !         'FFOOABAB','AOTWOSORT')
!
!! remaining ints
! call tran4_gen(NBas,&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          'OOOOABBB','AOTWOSORT')
! call tran4_gen(NBas,&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          'OOOOBAAA','AOTWOSORT')
! call tran4_gen(NBas,&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
!          'OOOOAAAA','AOTWOSORT')
! call tran4_gen(NBas,&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
!          'OOOOBBBB','AOTWOSORT')
!
! e2    = 0
! e2s2  = 0
!
! t12 = 0
! t34 = 0
! t13 = 0
!
! allocate(ints(NBas,NBas))
!
!! ! test exchange
!! !(FF|OO):(AB|AB)
!!  open(newunit=iunit,file='FFOOABAB',status='OLD', &
!!     access='DIRECT',recl=8*NBas*NBas)
!!
!! work1 = 0
!! ints  = 0
!! kl    = 0
!! tmp =0
!! do l=1,dimOB
!!    do k=1,dimOA
!!       kl = kl + 1
!!       read(iunit,rec=kl) work1(1:NBas*NBas)
!!
!!          do j=1,NBas
!!             do i=1,NBas
!!                ints(i,j) = work1((j-1)*NBas+i)
!!             enddo
!!          enddo
!!
!!          ip = k
!!          iq = l
!!
!!          tmp = tmp + A%Occ(ip)*B%Occ(iq)*ints(ip,iq)
!!
!!    enddo
!! enddo
!! print*, 'EXCH-TEST',-2d0*tmp
!!
!! close(iunit)
!! !!!!! end test
!
! !(OO|OO):(AA|AA)
! open(newunit=iunit,file='OOOOAAAA',status='OLD', &
!     access='DIRECT',recl=8*dimOA*dimOA)
!
! work1  = 0
! ints   = 0
! ccaaaa = 0
! kl     = 0
! do l=1,dimOA
!    do k=1,dimOA
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOA*dimOA)
!
!       iq = l
!       ip = k
!
!       ij = 0
!       do j=1,dimOA
!          do i=1,dimOA
!             ij = ij + 1
!             !ints(i,j) = work1((j-1)*dimOA+i)
!             ints(i,j) = work1(ij)
!          enddo
!       enddo
!
!       ccaaaa = ccaaaa + A%CICoef(ip)*A%CICoef(iq)*ints(ip,iq)
!
!    enddo
! enddo
! t12(1) = t12(1) + ccaaaa
! close(iunit)
!
!!(OO|OO):(BB|BB)
! open(newunit=iunit,file='OOOOBBBB',status='OLD', &
!     access='DIRECT',recl=8*dimOB*dimOB)
!
! work1  = 0
! ints   = 0
! ccbbbb = 0
! kl     = 0
! do l=1,dimOB
!    do k=1,dimOB
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOB*dimOB)
!
!       iq = l
!       ip = k
!
!       ij = 0
!       do j=1,dimOB
!          do i=1,dimOB
!             ij = ij + 1
!             !ints(i,j) = work1((j-1)*dimOB+i)
!             ints(i,j) = work1(ij)
!          enddo
!       enddo
!
!       ccbbbb = ccbbbb + B%CICoef(ip)*B%CICoef(iq)*ints(ip,iq)
!
!    enddo
! enddo
! t34(1) = t34(1) + ccbbbb
! close(iunit)
!
! !(OO|OO):(AB|BB)
! open(newunit=iunit,file='OOOOABBB',status='OLD', &
!     access='DIRECT',recl=8*dimOA*dimOB)
!
! work1 = 0
! ints  = 0
! kl    = 0
! do l=1,dimOB
!    do k=1,dimOB
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOA*dimOB)
!
!       ir = l
!       iq = k
!
!       fac = B%CICoef(iq)*B%CICoef(ir)
!
!       ij = 0
!       do j=1,dimOB
!          do i=1,dimOA
!             ij = ij + 1
!             !ints(i,j) = work1((j-1)*dimOB+i)
!             ints(i,j) = work1(ij)
!          enddo
!       enddo
!
!       do ip=1,dimOA
!          t13(3) = t13(3) + fac*A%Occ(ip)*Sab(ip,iq)*ints(ip,ir)
!       enddo
!
!    enddo
! enddo
! close(iunit)
! t13(3) = -2d0*t13(3)
! t34(2) = t13(3)
!
! !(OO|OO):(BA|AA)
! open(newunit=iunit,file='OOOOBAAA',status='OLD', &
!     access='DIRECT',recl=8*dimOA*dimOB)
!
! work1 = 0
! ints  = 0
! kl    = 0
! do l=1,dimOA
!    do k=1,dimOA
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOA*dimOB)
!
!       ir = l
!       ip = k
!
!       fac = A%CICoef(ip)*A%CICoef(ir)
!
!       ij = 0
!       do j=1,dimOA
!          do i=1,dimOB
!             ij = ij + 1
!             !ints(i,j) = work1((j-1)*dimOA+i)
!             ints(i,j) = work1(ij)
!          enddo
!       enddo
!
!       do iq=1,dimOB
!          t13(4) = t13(4) + fac*B%Occ(iq)*Sab(ip,iq)*ints(iq,ir)
!       enddo
!
!    enddo
! enddo
! close(iunit)
! t13(4) = -2d0*t13(4)
! t12(2) = t13(4)
!
! ! test!
! open(newunit=iunit,file='TMPOOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
!
! allocate(tmpAB(dimOA,dimOA,dimOB,dimOB),Sab_save(Nbas,Nbas))
!
! tmpAB=0
! do is=1,dimOB
!    do iq=1,dimOB
!       do ir=1,dimOA
!          do ip=1,dimOA
!             tmpAB(ip,ir,iq,is) = tmpAB(ip,ir,iq,is) + &
!                                  A%CICoef(ip)*A%CICoef(ir)*Sab(ip,iq) &
!                                * B%CICoef(iq)*B%CICoef(is)*Sab(ir,is)
!          enddo
!       enddo
!    enddo
! enddo
!
! do is=1,dimOB
!    do iq=1,dimOB
!       do ir=1,dimOA
!          do ip=1,dimOA
!          write(LOUT,'(1x,a,4i2,f12.6)') 'ip,ir,iq,is',ip,ir,iq,is,tmpAB(ip,ir,iq,is)
!          enddo
!       enddo
!    enddo
! enddo
!
! print*, 'tmpAB-hl2e',norm2(tmpAB)
!
! ints = 0
!
! do ir=1,dimOA
!    do ip=1,dimOA
!      read(iunit,rec=ip+(ir-1)*dimOA) ints(1:dimOB,1:dimOB)
!
!      t13(5) = t13(5) + sum(ints(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))
!      !t13(5) = t13(5) + sum(tmpAB(ip,ir,1:dimOB,1:dimOB))
!
!    enddo
! enddo
! t13(5) = -2d0*t13(5)
! !print*, 'T2a-4',t13(5)
! deallocate(tmpAB,Sab_save)
!
! close(iunit)
!
! !(OO|OO):(AA|BB)
! open(newunit=iunit,file='OOOOAABB',status='OLD', &
!     access='DIRECT',recl=8*dimOA*dimOA)
!
! work1 = 0
! ints  = 0
! kl    = 0
! val   = 0
! tmp   = 0
! do l=1,dimOB
!    do k=1,dimOB
!
!!       kl = kl + 1
!!       read(iunit,rec=kl) work1(1:dimOA*dimOA)
!       read(iunit,rec=k+(l-1)*dimOB) ints(1:dimOA,1:dimOA)
!
!       is = l
!       iq = k
!
!       fac = B%CICoef(iq)*B%CICoef(is)
!
!      ! do ir=1,dimOA
!      !    do ip=1,dimOA
!      !    !t13(5) = t13(5) + fac*A%CICoef(ip)*A%CICoef(ir)*Sab(ip,iq)*Sab(ir,is)*ints(ip,ir)
!      !    t13(5) = t13(5) + ints(ip,ir)
!      !    enddo
!      ! enddo
!
!       if(k==l) then
!
!          ip = k
!
!          do iq=1,dimOA
!             t13(1) = t13(1) + A%Occ(iq)*B%Occ(ip)*ints(iq,iq)
!          enddo
!
!       endif
!    enddo
! enddo
! close(iunit)
! t13(1) = 4d0*t13(1)
! !t13(5) = -2d0*t13(5)
!
! !(OO|OO):(AB|AB)
! open(newunit=iunit,file='OOOOABAB',status='OLD', &
!     access='DIRECT',recl=8*dimOA*dimOB)
!
! work1 = 0
! ints  = 0
! kl    = 0
! tmp   = 0
! do l=1,dimOB
!    do k=1,dimOA
!
!       kl = kl + 1
!       read(iunit,rec=kl) work1(1:dimOA*dimOB)
!
!       ip = l
!       iq = k
!
!       ij = 0
!       do j=1,dimOB
!          do i=1,dimOA
!             ij = ij + 1
!             !ints(i,j) = work1((j-1)*dimOB+i)
!             ints(i,j) = work1(ij)
!          enddo
!       enddo
!
!       ! t13b
!       t13(2) = t13(2) + A%Occ(iq)*B%Occ(ip)*ints(iq,ip)
!       ! t12c = t34c
!       t12(3) = t12(3) + A%CICoef(iq)*B%CICoef(ip)*ints(iq,ip)
!
!       is = l
!       ir = k
!
!       fac = A%CICoef(ir)*B%CICoef(is)*Sab(ir,is)
!
!       do iq=1,dimOB
!          do ip=1,dimOA
!             t13(6) = t13(6) + fac*A%CICoef(ip)*B%CICoef(iq)*Sab(ip,iq)*ints(ip,iq)
!          enddo
!       enddo
!
!    enddo
! enddo
! close(iunit)
! t12(3) = ccs2*t12(3)
! t34(3) = t12(3)
!
! t13(2) = -2d0*t13(2)
! t13(6) =  2d0*t13(6)
!
! print*, ''
! print*, 'T12(a)',t12(1)
! print*, 'T12(b)',t12(2)
! print*, 'T12(c)',t12(3)
!
! print*, 'T34(a)',t34(1)
! print*, 'T34(b)',t34(2)
! print*, 'T34(c)',t34(3)
! print*, 'T13(a)[Coul]',t13(1)
! print*, 'T13(b)[Exch]',t13(2)
! print*, 'T13(c)      ',t13(3)
! print*, 'T13(d)      ',t13(4)
! print*, 'T13(e)      ',t13(5)
! print*, 'T13(f)      ',t13(6)
!
! !print*, 'T12-sum',sum(t12)
! !print*, 'T34-sum',sum(t34)
! !print*, 'T13-sum',sum(t13)
!
! e2 = sum(t12) + sum(t34) + sum(t13)
! e2s2 = sum(t12) + sum(t34) - 2d0*t12(3) + sum(t13) - t13(6)
! e2s2 = e2s2 + 2d0*nns2*ccaaaa + 2d0*nns2*ccbbbb + 2d0*nns2*t13(1)
!
! print*, ''
! print*, 'E1    ',e1
! print*, 'E1(S2)',e1s2
!
! print*, 'E2    ',e2
! print*, 'E2(S2)',e2s2
!
! eHL = (e1+e2)/Dfull
! ehls2 = e1s2 + e2s2
!
! print*, 'E(A),E(B)', A%ECASSCF,B%ECASSCF
! print*, 'E(A)+E(B)', A%ECASSCF+B%ECASSCF
!
! print*, ''
! print*, 'eHL    ', ehl
! print*, 'eHL(S2)', ehls2
!
! ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
! !ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn
! ehls2int = ehls2 - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
! print*, ''
! print*, 'eHLint       ', ehlint*1000d0
! print*, 'eHL(S2)int   ', ehls2int*1000d0
! print*, 'ELST+EXCH(s2)',(SAPT%elst+SAPT%exchs2)*1000d0
!
! !print*, 'HL-ELST    ',ehlint-SAPT%elst
! !print*, 'HL(S2)-ELST',ehls2int-SAPT%elst
! !
! !tmp = ehls2int-SAPT%elst-SAPT%exchs2
! !print*, 'test     ',tmp
! tmp = ehlint-SAPT%elst-SAPT%exchs2
! print*, 'delta     ',tmp*1000d0
!
! print*, 'TEST-exch'
! t12(2) = SAPT%exch_part(2)
! t13(2) = SAPT%exch_part(1)
! t34(2) = SAPT%exch_part(3)
! t13(5) = SAPT%exch_part(4)
! t13(6) = SAPT%exch_part(5)
! e2 = sum(t12) + sum(t34) + sum(t13)
! eHL = (e1+e2)/Dfull
! print*, 'eHL-tst', ehl
! ehlint = ehl - A%ECASSCF - B%ECASSCF + SAPT%Vnn + SAPT%monA%PotNuc + SAPT%monB%PotNuc
! print*, 'eHLint       ', ehlint*1000d0
!
! !! TESTY S2
! write(LOUT,'(/,1x,a)') 'TESTY-S2:'
! print*, 'Pb.Ja',t13(1)
! print*, 'nns2 ',nns2
! print*, 'T1   ',2d0*nns2*t13(1)
! print*, 'T2a-1',t13(2)
! print*, 'T2a-2',t12(2)
! print*, 'T2a-3',t34(2)
! print*, 'T2a-4',t13(5)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! test eq 58
! !! MONOMER A
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !write(LOUT,'(/,1x,a)') 'Test Eq. 58'
! !Haa = 0; Hab = 0
! !call tran2MO(Ha,A%CMO,A%CMO,Haa,NBas)
! !call tran2MO(Ha,B%CMO,A%CMO,Hab,NBas)
!
! !! test 1-el
! !tmp=0
! !do i=1,dimOA
! !   tmp = tmp + 2d0*A%Occ(i)*Haa(i,i)
! !enddo
! !print*, 'OneEl(A)',tmp
!
! !work = 0
! !do ir=1,NBas
! !   do iq=1,dimOA
! !   !do iq=1,NBas
!
! !     work(iq,ir) = work(iq,ir) +  A%CICoef(iq)*Hab(ir,iq)
!
! !     val = 0
! !     do ip=1,dimOA
! !        val = val + A%CICoef(ip)*Sab(ip,ir)*Haa(iq,ip)
! !     enddo
! !     work(iq,ir) = work(iq,ir) + val
!
! !   enddo
! !enddo
!
! !!!!(FF|FF):(BA|AA)
! !!call tran4_gen(NBas,&
! !!         NBas,A%CMO,&
! !!         NBas,A%CMO,&
! !!         NBas,B%CMO,&
! !!         NBas,A%CMO,&
! !!         'FFFFBAAA','AOTWOSORT')
!
! !!open(newunit=iunit,file='FFFFBAAA',status='OLD', &
! !!    access='DIRECT',recl=8*NBas*NBas)
!
! !!work1 = 0
! !!ints  = 0
! !!kl    = 0
! !!do l=1,NBas
! !!   do k=1,NBas
! !!      kl = kl + 1
!
! !!      if((l.le.dimOA).and.(k.le.dimOA)) then
!
! !!         read(iunit,rec=kl) work1(1:NBas*NBas)
!
! !!         ip = l
! !!         iq = k
!
! !!         do j=1,NBas
! !!            do i=1,NBas
! !!               ints(i,j) = work1((j-1)*NBas+i)
! !!            enddo
! !!         enddo
!
! !!         do ir=1,NBas
! !!            work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(ir,ip)
! !!         enddo
!
! !!     endif
! !!   enddo
! !!enddo
! !!close(iunit)
! !!call delfile('FFFFBAAA')
! !!!(FO|OO):(BA|AA)
! !call tran4_gen(NBas,&
! !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
! !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
! !         NBas,B%CMO,&
! !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
! !         'FOOOBAAA','AOTWOSORT')
!
! !open(newunit=iunit,file='FOOOBAAA',status='OLD', &
! !    access='DIRECT',recl=8*NBas*dimOA)
!
! !work1 = 0
! !ints  = 0
! !kl    = 0
! !do l=1,dimOA
! !   do k=1,dimOA
! !      kl = kl + 1
!
! !         read(iunit,rec=kl) work1(1:NBas*dimOA)
!
! !         ip = l
! !         iq = k
!
! !         do j=1,NBas
! !            do i=1,NBas
! !               ints(i,j) = work1((j-1)*NBas+i)
! !            enddo
! !         enddo
!
! !         do ir=1,NBas
! !            work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(ir,ip)
! !         enddo
!
! !   enddo
! !enddo
! !close(iunit)
! !call delfile('FOOOBAAA')
!
! !do ir=1,NBas
! !do iq=1,NBas
! !   val = A%ECASSCF*A%CICoef(iq)*Sab(iq,ir)
! !   !print*,'iq,ir',iq,ir,work(iq,ir),val
! !   print*,'iq,ir',iq,ir,work(iq,ir)-val
! !enddo
! !enddo
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! test eq 58
! !! MONOMER B
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hbb = 0; Hab = 0
! call tran2MO(Hb,B%CMO,B%CMO,Hbb,NBas)
! call tran2MO(Hb,A%CMO,B%CMO,Hab,NBas)
!
! ! test 1-el
! tmp=0
! do i=1,dimOB
!    tmp = tmp + 2d0*B%Occ(i)*Hbb(i,i)
! enddo
! print*, 'OneEl(B)',tmp
!
! !work = 0
! !do ir=1,NBas
! !   do iq=1,dimOB
! !   !do iq=1,NBas
!
! !     work(iq,ir) = work(iq,ir) +  B%CICoef(iq)*Hab(iq,ir)
!
! !     val = 0
! !     do ip=1,dimOB
! !        val = val + B%CICoef(ip)*Sab(ip,ir)*Hbb(iq,ip)
! !     enddo
! !     work(iq,ir) = work(iq,ir) + val
!
! !   enddo
! !enddo
!
!!!(FO|OO):(AB|BB)
! !call tran4_gen(NBas,&
! !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
! !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
! !         NBas,A%CMO,&
! !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
! !         'FOOOABBB','AOTWOSORT')
!
! !open(newunit=iunit,file='FOOOABBB',status='OLD', &
! !    access='DIRECT',recl=8*NBas*dimOB)
!
! !work1 = 0
! !ints  = 0
! !kl    = 0
! !do l=1,dimOB
! !   do k=1,dimOB
! !      kl = kl + 1
!
! !         read(iunit,rec=kl) work1(1:NBas*dimOB)
!
! !         ip = l
! !         iq = k
!
! !         do j=1,NBas
! !            do i=1,NBas
! !               ints(i,j) = work1((j-1)*NBas+i)
! !            enddo
! !         enddo
!
! !         do ir=1,NBas
! !            work(iq,ir) = work(iq,ir) + B%CICoef(ip)*ints(ir,ip)
! !         enddo
!
! !   enddo
! !enddo
! !close(iunit)
! !call delfile('FOOOABBB')
!
! !do ir=1,NBas
! !do iq=1,NBas
! !   val = B%ECASSCF*B%CICoef(iq)*Sab(iq,ir)
! !   print*,'iq,ir',iq,ir,work(iq,ir),val
! !enddo
! !enddo
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! test eq 59
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! !Haa = 0; Hab = 0
! !call tran2MO(Ha,A%CMO,A%CMO,Haa,NBas)
! !work = 0
! !!(OO|OO):(AA|AA)
! !open(newunit=iunit,file='OOOOAAAA',status='OLD', &
! !    access='DIRECT',recl=8*dimOA*dimOA)
!
! !work1  = 0
! !ints   = 0
! !kl     = 0
! !do l=1,dimOA
! !   do k=1,dimOA
! !      kl = kl + 1
! !      read(iunit,rec=kl) work1(1:dimOA*dimOA)
!
! !      ir = l
! !      ip = k
!
! !      do j=1,dimOA
! !         do i=1,dimOA
! !            ints(i,j) = work1((j-1)*dimOA+i)
! !         enddo
! !      enddo
!
! !      do iq=1,dimOA
! !         work(iq,ir) = work(iq,ir) + A%CICoef(ip)*ints(iq,ip)
! !      enddo
!
! !   enddo
! !   ir = l
! !   do iq=1,NBas
! !    work(iq,ir) = work(iq,ir) +  A%CICoef(ir)*Haa(iq,ir) + A%CICoef(iq)*Haa(ir,iq)
! !   enddo
!
! !enddo
! !close(iunit)
!
! !do ir=1,Nbas
! !do iq=1,Nbas
! !   print*, 'iq,ir',iq,ir,work(iq,ir)
! !enddo
! !enddo
!
! !do iq=1,dimOA
! !   print*, 'A-iq',iq,work(iq,iq),A%ECASSCF*A%CICoef(iq)
! !enddo
!
! !work = 0
! !!(OO|OO):(BB|BB)
! !open(newunit=iunit,file='OOOOBBBB',status='OLD', &
! !    access='DIRECT',recl=8*dimOB*dimOB)
!
! !work1  = 0
! !ints   = 0
! !kl     = 0
! !do l=1,dimOB
! !   do k=1,dimOB
! !      kl = kl + 1
! !      read(iunit,rec=kl) work1(1:dimOB*dimOB)
!
! !      ir = l
! !      ip = k
!
! !      do j=1,dimOB
! !         do i=1,dimOB
! !            ints(i,j) = work1((j-1)*dimOB+i)
! !         enddo
! !      enddo
!
! !      do iq=1,dimOB
! !         work(iq,ir) = work(iq,ir) + B%CICoef(ip)*ints(iq,ip)
! !      enddo
!
! !   enddo
! !   ir = l
! !   do iq=1,NBas
! !    work(iq,ir) = work(iq,ir) +  B%CICoef(ir)*Hbb(iq,ir) + B%CICoef(iq)*Hbb(ir,iq)
! !   enddo
!
! !enddo
! !close(iunit)
!
! !do iq=1,dimOB
! !   print*, 'B-iq',iq,work(iq,iq),B%ECASSCF*B%CICoef(iq)
! !enddo
! !! end test eq 59
!
! deallocate(ints)
!
!
! deallocate(Hab,Hbb,Hb,Haa,Ha,work1,work)
! deallocate(PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
! deallocate(RDM2Bval,RDM2Aval)
!
! ! delete ints
! call delfile('OOOOABAB')
! call delfile('OOOOAABB')
! call delfile('OOOOABBB')
! call delfile('OOOOBAAA')
! call delfile('OOOOAAAA')
! call delfile('OOOOBBBB')
!
!end subroutine hl_2el

subroutine e2exind(Flags,A,B,SAPT)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer             :: NAO,NBas
integer             :: dimOA,dimOB,dimVA,dimVB,nOVA,nOVB
double precision    :: nelA,nelB
integer,allocatable :: posA(:,:),posB(:,:)
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAaa(:,:),PBbb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: uA(:),uB(:)
double precision,allocatable :: tindA(:),tindB(:),&
                                tindX(:),VindX(:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
! unc
type(EBlockData)             :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable :: SBlockA(:),SBlockB(:)
character(:),allocatable     :: propA,propB

double precision :: termZ_u,termY_u,termX_u
double precision :: e2exi_unc
double precision,allocatable :: uA0(:),uB0(:)
double precision,allocatable :: OmA0(:),OmB0(:)
!double precision,allocatable :: OmA01(:),OmB01(:)
integer                      :: iblk,nblkA,nblkB
! full test
!double precision,allocatable :: AVecX0(:),AVecY0(:), &
!                                BVecX0(:),BVecY0(:)

double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: tmpXA(:),tmpYA(:),&
                                tmpXB(:),tmpYB(:)
integer :: i,j,ipq,ip,iq,irs,ir,is
logical :: both,uncoupled
double precision :: termZ,termY,termX
double precision :: e2exi
double precision :: fact,tmp
! Thresholds
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8

! avoid tran4 in exch-disp
 SAPT%noE2exi = .false.
 both = SAPT%iCpld

 uncoupled = .true.
 if(Flags%ICASSCF==0)   uncoupled = .false.
 if(Flags%ITREXIO==1)   uncoupled = .false.
 if(A%Cubic.or.B%Cubic) uncoupled = .false.

! set dimensions
 NAO  = SAPT%NAO
 NBas = A%NBasis
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 ! SAPT cubic
 ! A: get response properites
 if(A%Cubic) then
   allocate(A%EigX(A%NDimX*A%NDimX),A%EigY(A%NDimX*A%NDimX),&
            A%Eig(A%NDimX))
   if(A%ACAlpha==A%ACAlpha0) then
      propA = 'PROP_A0'
   elseif(A%ACAlpha==A%ACAlpha1) then
      propA = 'PROP_A1'
   elseif(A%ACAlpha==A%ACAlpha2) then
      propA = 'PROP_A2'
   endif
   call readEVecXY(A%EigX,A%EigY,A%NDimX,propA)
   call readEvalXY(A%Eig,A%NDimX,propA)
 endif
 ! B: get response properties
 if(B%Cubic) then
   allocate(B%EigX(B%NDimX*B%NDimX),B%EigY(B%NDimX*B%NDimX),&
            B%Eig(B%NDimX))

   if(B%ACAlpha==B%ACAlpha0) then
      propB = 'PROP_B0'
   elseif(B%ACAlpha==B%ACAlpha1) then
      propB = 'PROP_B1'
   elseif(B%ACAlpha==B%ACAlpha2) then
      propB = 'PROP_B2'
   endif
   call readEVecXY(B%EigX,B%EigY,B%NDimX,propB)
   call readEvalXY(B%Eig,B%NDimX,propB)
 endif

! print thresholds
 if(SAPT%IPrint>5) then
    write(LOUT,'(/,1x,a)') 'Thresholds in E2exch-ind:'
    write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
    write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
 endif

 ! unc
 if(uncoupled) then
    allocate(OmA0(A%NDimX),OmB0(B%NDimX))

    call read_SBlock(SBlockA,SBlockAIV,nblkA,'XY0_A')
    call read_SBlock(SBlockB,SBlockBIV,nblkB,'XY0_B')

    call unpack_Eig(SBlockA,SBlockAIV,nblkA,OmA0,A%NDimX)
    call unpack_Eig(SBlockB,SBlockBIV,nblkB,OmB0,B%NDimX)
 endif

 ! full unc
 !if(Flags%ICASSCF==1) then
 !   allocate(AVecX0(A%NDimX*A%NDimX),OmA0(A%NDimX), &
 !            AVecY0(A%NDimX*A%NDimX), &
 !            BVecX0(B%NDimX*B%NDimX),OmB0(B%NDimX), &
 !            BVecY0(B%NDimX*B%NDimX))

 !   call unpack_XY0_full(AVecX0,AVecY0,OmA0,A%CICoef,A%IndN,A%NDimX,NBas,'XY0_A')
 !   call unpack_XY0_full(BVecX0,BVecY0,OmB0,B%CICoef,B%IndN,B%NDimX,NBas,'XY0_B')
 !endif

! read EigValA_B
! allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX), &
!          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(S(NAO,NAO),&
          PA(NAO,NAO),PB(NAO,NAO),&
          Va(NAO,NAO),Vb(NAO,NAO),&
          Sab(NBas,NBas),Sba(NBas,NBas),&
          PAaa(NBas,NBas),PBbb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas),&
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call get_one_mat('S',S,A%Monomer,NAO)
 !call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 !call tran2MO(S,B%CMO,A%CMO,Sba,NBas)
 call tran_AO2MO2(S,A%CMO,B%CMO,Sab,NAO,NBas)
 call tran_AO2MO2(S,B%CMO,A%CMO,Sba,NAO,NBas)

 call get_den(NAO,NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NAO,NBas,B%CMO,B%Occ,1d0,PB)

 allocate(tmp1(NAO,NBas),tmp2(NAO,NBas))

 call dgemm('N','N',NAO,NBas,NAO,1d0,S,NAO,A%CMO,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NBas,NAO,1d0,S,NAO,B%CMO,NAO,0d0,tmp2,NAO)
 ! PA(B), PB(A)
 !call tran2MO(PA,tmp1,tmp1,PAaa,NBas)
 !call tran2MO(PB,tmp2,tmp2,PBbb,NBas)
 call tran_AO2MO2(PA,tmp1,tmp1,PAaa,NAO,NBas)
 call tran_AO2MO2(PB,tmp2,tmp2,PBbb,NAO,NBas)

 call get_one_mat('V',Va,A%Monomer,NAO)
 call get_one_mat('V',Vb,B%Monomer,NAO)

 !call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 !call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)
 call tran_AO2MO2(Va,B%CMO,B%CMO,Vabb,NAO,NBas)
 call tran_AO2MO2(Vb,A%CMO,A%CMO,Vbaa,NAO,NBas)

 !call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 !call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 !call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)
 call tran_AO2MO2(Va,B%CMO,A%CMO,Vaba,NAO,NBas)
 call tran_AO2MO2(Va,A%CMO,B%CMO,Vaab,NAO,NBas)
 call tran_AO2MO2(Vb,A%CMO,B%CMO,Vbab,NAO,NBas)

 !call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas)
 !call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas)
 call tran_AO2MO2(A%WPot,B%CMO,B%CMO,WaBB,NAO,NBas)
 call tran_AO2MO2(B%WPot,A%CMO,A%CMO,WbAA,NAO,NBas)

 !allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
 !         RDM2Bval(dimOB,dimOB,dimOB,dimOB))


 allocate(posA(NBas,NBas),posB(NBas,NBas))

 if(Flags%ICASSCF==1) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
    enddo
 elseif(Flags%ICASSCF==0) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = i
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = i
    enddo
 endif

 deallocate(tmp2,tmp1)
 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 allocate(tmp1(NAO,NAO),tmp2(NAO,NAO))
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NAO,NAO,NAO,1d0,PA,NAO,S,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NAO,NAO,1d0,S,NAO,PB,NAO,0d0,tmp2,NAO)
 do j=1,NAO
    do i=1,NAO
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 if(uncoupled) termZ_u = 2d0*SAPT%e2ind_unc*termZ

 if(A%Cubic.or.B%Cubic) then
   if(A%ACAlpha==A%ACAlpha0.or.B%ACAlpha==B%ACAlpha0) termZ = 2d0*SAPT%e2ind_a0*termZ
   if(A%ACAlpha==A%ACAlpha1.or.B%ACAlpha==B%ACAlpha1) termZ = 2d0*SAPT%e2ind_a1*termZ
   if(A%ACAlpha==A%ACAlpha2.or.B%ACAlpha==B%ACAlpha2) termZ = 2d0*SAPT%e2ind_a2*termZ
 else
   termZ = 2d0*SAPT%e2ind*termZ
 endif

 deallocate(tmp2,tmp1)

!! Y term

 allocate(tindA(A%NDimX),tindB(B%NDimX))

 ! unc
 if(uncoupled) then

    allocate(uA0(A%NDimX),uB0(B%NDimX))
    uA0 = 0
    tindA = 0
    do i=1,A%NDimX
       ip = A%IndN(1,i)
       iq = A%IndN(2,i)
       ipq = posA(ip,iq)

       tindA(ipq) = tindA(ipq) + (A%Occ(ip)-A%Occ(iq))*WbAA(ip,iq)

    enddo

    call abpm_dgemv_gen(tindA,uA0,SBlockA,SBlockAIV,nblkA,A%NDimX,'YX')
    ! full unc
    !call dgemv('T',A%NDimX,A%NDimX,1d0,AVecY0-AVecX0,A%NDimX,tindA,1,0d0,uA0,1)

    uB0 = 0
    tindB = 0
    do j=1,B%NDimX
       ir = B%IndN(1,j)
       is = B%IndN(2,j)
       irs = posB(ir,is)

       tindB(irs) = tindB(irs) + (B%Occ(ir)-B%Occ(is))*WaBB(ir,is)

    enddo

    call abpm_dgemv_gen(tindB,uB0,SBlockB,SBlockBIV,nblkB,B%NDimX,'YX')
    ! full unc
    !call dgemv('T',B%NDimX,B%NDimX,1d0,BVecY0-BVecX0,B%NDimX,tindB,1,0d0,uB0,1)

    call make_tind_unc(tindA,SBlockA,SblockAIV,Sab,A%Occ,B%Occ,A%IndN,posA,nblkA,A%NDimX,NBas)
    call make_tind_unc(tindB,SblockB,SblockBIV,Sba,B%Occ,A%Occ,B%IndN,posB,nblkB,B%NDimX,NBas)

    ! full unc
    !call make_tind(tindA,AVecX0,AVecY0,Sab,A%Occ,B%Occ,A%IndN,posA,A%NDimX,NBas)
    !call make_tind(tindB,BVecX0,BVecY0,Sba,B%Occ,A%Occ,B%IndN,posB,B%NDimX,NBas)

    termY_u = 0d0
    do i=1,A%NDimX
       if(abs(OmA0(i)).gt.SmallE.and.abs(OmA0(i)).lt.BigE) then
          termY_u = termY_u + (tindA(i)*uA0(i))/OmA0(i)
       endif
    enddo
    do i=1,B%NDimX
       if(abs(OmB0(i)).gt.SmallE.and.abs(OmB0(i)).lt.BigE) then
          termY_u = termY_u + (tindB(i)*uB0(i))/OmB0(i)
       endif
    enddo

    termY_u=-4d0*(SAPT%elst-SAPT%Vnn)*termY_u

 endif

 ! cpld
 if(both) then

    allocate(uA(A%NDimX),uB(B%NDimX))

    uA = 0
    tindA = 0
    do i=1,A%NDimX
       ip = A%IndN(1,i)
       iq = A%IndN(2,i)
       ipq = posA(ip,iq)

       tindA(ipq) = tindA(ipq) + (A%Occ(ip)-A%Occ(iq))*WbAA(ip,iq)

    enddo
    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigY-A%EigX,A%NDimX,tindA,1,0d0,uA,1)
    !print*, 'uA',norm2(uA)

    uB = 0
    tindB = 0
    do j=1,B%NDimX
       ir = B%IndN(1,j)
       is = B%IndN(2,j)
       irs = posB(ir,is)

       tindB(irs) = tindB(irs) + (B%Occ(ir)-B%Occ(is))*WaBB(ir,is)

    enddo
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigY-B%EigX,B%NDimX,tindB,1,0d0,uB,1)
    !print*, 'uB',norm2(uB)

     call make_tind(tindA,A%EigX,A%EigY,Sab,A%Occ,B%Occ,A%IndN,posA,A%NDimX,NBas)
     call make_tind(tindB,B%EigX,B%EigY,Sba,B%Occ,A%Occ,B%IndN,posB,B%NDimX,NBas)

    termY=0d0
    do i=1,A%NDimX
       if(abs(A%Eig(i)).gt.SmallE.and.abs(A%Eig(i)).lt.BigE) then
          termY = termY + (tindA(i)*uA(i))/A%Eig(i)
       endif
    enddo
    do i=1,B%NDimX
       if(abs(B%Eig(i)).gt.SmallE.and.abs(B%Eig(i)).lt.BigE) then
          termY = termY + (tindB(i)*uB(i))/B%Eig(i)
       endif
    enddo

    termY=-4d0*(SAPT%elst-SAPT%Vnn)*termY
    !write(*,*) 'termY',termY

 endif
 deallocate(tindB,tindA)

 ! term A3
 !write(LOUT,'()')
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         'FOFOAABB','AOTWOSORT')
 !! term A1
 !call tran4_gen(NBas,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         NBas,B%CMO,&
 !         'FFOOABAB','AOTWOSORT')
 !! term A2
 !! A2A(B): XX
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,B%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOFOBBBA','AOTWOSORT')
 !call tran4_gen(NBas,&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         'FOFOAAAB','AOTWOSORT')
 !!! A2A(B): YY
 !call tran4_gen(NBas,&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,B%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOFOBBAB','AOTWOSORT')
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,A%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         'FOFOAABA','AOTWOSORT')

 ! A3 prototype
 allocate(tmpXA(A%NDimX),tmpYA(A%NDimX),&
          tmpXB(B%NDimX),tmpYB(B%NDimX))
 tmpXA = 0
 tmpXB = 0
 tmpYA = 0
 tmpYB = 0
 !call exind_A3_XY(A%NDimX,A%NDimX,tmpXA,tmpYA,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
 !                 B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
 !
 !! test B
 !print*, 'test A3-B!'
 !call exind_A3_XY(B%NDimX,B%NDimX,tmpXB,tmpYB,RDM2Bval,RDM2Aval,Sba,nelB,Vbaa,nelA,Vabb,'FOFOBBAA',&
 !                 A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas)

 ! A3 here!!
 call exind_A3_XY_full(A%NDimX,B%NDimX,tmpXA,tmpXB,A%RDM2val,B%RDM2val,Sab, &
                    nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                    B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 tmpYA = -tmpXA
 tmpYB = -tmpXB

! A1
 call exind_A1_AB(A%NDimX,B%NDimX,tmpXA,tmpXB,tmpYA,tmpYB,nelA,nelB,dimOA,dimOB, &
                  A%Occ,B%Occ,PAaa,PBbb,Vaab,Vbab,Sab,A%IndN,B%IndN,posA,posB,NBas)
! A2
 call exind_A2_XX(A%NDimX,B%NDimX,tmpXA,tmpXB,B%RDM2val, &
                  Sab,nelA,Vabb,nelB,Vbab,PAaa,'FOFOBBBA', &
                  B%Occ,A%Occ,B%IndN,A%IndN,posB,posA,  &
                  dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
 call exind_A2_XX(B%NDimX,A%NDimX,tmpXB,tmpXA,A%RDM2val, &
                 Sba,nelB,Vbaa,nelA,Vaba,PBbb,'FOFOAAAB', &
                  A%Occ,B%Occ,A%IndN,B%IndN,posA,posB,  &
                  dimOA,dimOB,A%NDimX,B%NDimX,NBas,.false.)

 !print*, 'XA-tot',norm2(tmpXA)
 !print*, 'XB-tot',norm2(tmpXB)

 call exind_A2_YY(A%NDimX,B%NDimX,tmpYA,tmpYB,B%RDM2val, &
                  Sab,nelA,Vabb,nelB,Vbab,PAaa,'FOFOBBAB', &
                  B%Occ,A%Occ,B%IndN,A%IndN,posB,posA,  &
                  dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)

 !print*, 'YA-tot-1',norm2(tmpYA)
 !print*, 'YB-tot-1',norm2(tmpYB)

 call exind_A2_YY(B%NDimX,A%NDimX,tmpYB,tmpYA,A%RDM2val, &
                  Sba,nelB,Vbaa,nelA,Vaba,PBbb,'FOFOAABA', &
                  A%Occ,B%Occ,A%IndN,B%IndN,posA,posB,  &
                  dimOA,dimOB,A%NDimX,B%NDimX,NBas,.false.)

 !print*, 'YA-tot-2',norm2(tmpYA)
 !print*, 'YB-tot-2',norm2(tmpYB)

 ! unc
 if(uncoupled) then

    allocate(VindX(A%NDimX))
    VindX = 0d0
    call abpm_dgemv_gen(tmpXA,VindX,SBlockA,SBlockAIV,nblkA,A%NDimX,'X')
    call abpm_dgemv_gen(tmpYA,VindX,SBlockA,SBlockAIV,nblkA,A%NDimX,'Y')

    ! full unc
    !call dgemv('T',A%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmpXA,1,0d0,VindX,1)
    !call dgemv('T',A%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmpYA,1,1d0,VindX,1)

    termX_u = 0
    do i=1,A%NDimX
       if(abs(OmA0(i)).gt.SmallE.and.abs(OmA0(i)).lt.BigE) then
          termX_u = termX_u + (VindX(i)*uA0(i))/OmA0(i)
       endif
    enddo

    deallocate(VindX)
    allocate(VindX(B%NDimX))
    VindX = 0d0
    call abpm_dgemv_gen(tmpXB,VindX,SBlockB,SBlockBIV,nblkB,B%NDimX,'X')
    call abpm_dgemv_gen(tmpYB,VindX,SBlockB,SBlockBIV,nblkB,B%NDimX,'Y')

    ! full unc
    !call dgemv('T',B%NDimX,B%NDimX,1d0,BVecX0,B%NDimX,tmpXB,1,0d0,VindX,1)
    !call dgemv('T',B%NDimX,B%NDimX,1d0,BVecY0,B%NDimX,tmpYB,1,1d0,VindX,1)

    do i=1,B%NDimX
       if(abs(OmB0(i)).gt.SmallE.and.abs(OmB0(i)).lt.BigE) then
          termX_u = termX_u + (VindX(i)*uB0(i))/OmB0(i)
       endif
    enddo
    !print*, 'termX_unc-2',termX_u
    termX_u = -2d0*termX_u

    e2exi_unc = termX_u + termY_u + termZ_u

    SAPT%e2exind_unc = e2exi_unc
    !write(LOUT,'(/1x,a,f16.8)') 'E2exch-ind(unc) =', e2exi_unc*1.0d3
    call print_en('E2exch-ind(unc)',e2exi_unc*1.0d3,.true.)

    deallocate(VindX)
    deallocate(uB0,uA0)

 endif

 ! cpld
 if(both) then
    allocate(VindX(A%NDimX))

    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmpXA,1,0d0,VindX,1)
    call dgemv('T',A%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmpYA,1,1d0,VindX,1)

    termX = 0
    do i=1,A%NDimX
       if(abs(A%Eig(i)).gt.SmallE.and.abs(A%Eig(i)).lt.BigE) then
          termX = termX + (VindX(i)*uA(i))/A%Eig(i)
       endif
    enddo
    !print*, 'termX-1',termX

    deallocate(VindX)
    allocate(VindX(B%NDimX))

    VindX=0
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigX,B%NDimX,tmpXB,1,0d0,VindX,1)
    call dgemv('T',B%NDimX,B%NDimX,1d0,B%EigY,B%NDimX,tmpYB,1,1d0,VindX,1)

    do i=1,B%NDimX
       if(abs(B%Eig(i)).gt.SmallE.and.abs(B%Eig(i)).lt.BigE) then
          termX = termX + (VindX(i)*uB(i))/B%Eig(i)
       endif
    enddo
    !print*, 'termX-2',termX
    termX=-2d0*termX

    !print*, 'termX',termX
    !if(SAPT%IPrint>=5) write(LOUT,'(/1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3
    !if(SAPT%IPrint>=5) write(LOUT,'(1x,a,f16.8)')  'term Y      = ',  termY*1.0d3
    !if(SAPT%IPrint>=5) write(LOUT,'(1x,a,f16.8)')  'term X      = ',  termX*1.0d3

    if(SAPT%IPrint>=5) call print_en('term Z',termZ*1.0d3,.true.) 
    if(SAPT%IPrint>=5) call print_en('term Y',termY*1.0d3,.false.)
    if(SAPT%IPrint>=5) call print_en('term X',termX*1.0d3,.false.)

    e2exi = termX + termY + termZ

    if(A%Cubic.or.B%Cubic) then

      ! cubic
      if(A%ACAlpha==A%ACAlpha0.or.B%ACAlpha==B%ACAlpha0) SAPT%e2exi_a0 = e2exi
      if(A%ACAlpha==A%ACAlpha1.or.B%ACAlpha==B%ACAlpha1) SAPT%e2exi_a1 = e2exi
      if(A%ACAlpha==A%ACAlpha2.or.B%ACAlpha==B%ACAlpha2) SAPT%e2exi_a2 = e2exi

      !write(LOUT,'(1x,a,f16.8)') 'A: Alpha     = ',A%ACAlpha
      !write(LOUT,'(1x,a,f16.8)') 'B: Alpha     = ',B%ACAlpha
      write(LOUT,'(1x,a,f16.8)') 'E2exi(Alpha)   = ',e2exi*1000d0

    else

      ! regular
      SAPT%e2exind = e2exi
      !write(LOUT,'(/1x,a,f16.8)') 'E2exch-ind  = ', e2exi*1.0d3
      call print_en('E2exch-ind',e2exi*1.0d3,.true.)

    endif

    deallocate(VindX)
    deallocate(uB,uA)

 endif

 if(uncoupled) then

    ! unc
    ! deallocate SBLOCK
    do iblk=1,nblkA
       associate(A => SBlockA(iblk))
         deallocate(A%matY,A%matX,A%vec)
         deallocate(A%pos)
       end associate
    enddo
    do iblk=1,nblkB
       associate(B => SBlockB(iblk))
         deallocate(B%matY,B%matX,B%vec)
         deallocate(B%pos)
       end associate
    enddo
    ! deallocate IV part
    associate(A => SBlockAIV)
      if(A%n>0) then
         deallocate(A%vec)
         deallocate(A%pos)
      endif
    end associate
    associate(B => SBlockBIV)
      if(B%n>0) then
         deallocate(B%vec)
         deallocate(B%pos)
      endif
    end associate

    deallocate(OmB0,OmA0)

    ! full unc
    ! deallocate(AVecY0,AVecX0,BVecY0,BVecX0)
 endif

 ! SAPT cubic
 if(A%Cubic) deallocate(A%Eig,A%EigY,A%EigX)
 if(B%Cubic) deallocate(B%Eig,B%EigY,B%EigX)

! deallocate(VindB,VindA)
 deallocate(tmpYB,tmpXB)
 deallocate(tmpYA,tmpXA)

 deallocate(posB,posA)
 deallocate(PBbb,PAaa)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vb,Va)
 deallocate(WbAA,WaBB,PB,PA,Sba,Sab,S)

end subroutine e2exind

subroutine e2exdisp(Flags,A,B,SAPT)

use exappr
!use sref_exch

implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer          :: NAO,NBas,NInte1
integer          :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer          :: GdimOA,GdimOB
integer          :: iunit
integer          :: i,j,k,l,kl,ij,ii,jj,pq,rs
integer          :: ip,iq,ir,is,ipq,irs
integer          :: iblk,nblkA,nblkB
logical          :: both,approx,ipropab
logical          :: uncoupled
double precision :: nelA,nelB
double precision :: fact,val,termZ,termY,termX
double precision :: termZ_u,termY_u,termX_u
double precision :: e2exd_u,e2exd
integer,allocatable          :: posA(:,:),posB(:,:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:),&
                                sij(:,:),  &
                                tmp_u(:,:),&
                                work(:),workSq(:,:)
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
type(EBlockData)             :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable :: SBlockA(:),SBlockB(:)
character(:),allocatable     :: propA,propB
! uncoupled full test
double precision,allocatable :: AVecX0(:),AVecY0(:), &
                                BVecX0(:),BVecY0(:), &
                                OmA0(:),OmB0(:)
double precision,external    :: trace,FRDM2,FRDM2GVB
! test for Be
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE   = 1.D8
double precision,parameter :: SmallE = 1.D-3

! set dimensions
 NAO  = SAPT%NAO
 NBas = A%NBasis
 !dimOA = A%INAct+A%NAct
 dimOA = A%num0+A%num1
 GdimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 !dimOB = B%INAct+B%NAct
 dimOB = B%num0+B%num1
 GdimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

! set e2exd_version
 approx    = .false.
 !approx    = .true.
 both      = SAPT%iCpld
 uncoupled = .true.
 if(Flags%ICASSCF==0)   uncoupled = .false.
 if(Flags%ITREXIO==1)   uncoupled = .false.
 if(A%Cubic.or.B%Cubic) uncoupled = .false.

 ! SAPT cubic
 ! A: get response properites
 if(A%Cubic) then
   allocate(A%EigX(A%NDimX*A%NDimX),A%EigY(A%NDimX*A%NDimX),&
            A%Eig(A%NDimX))
   if(A%ACAlpha==A%ACAlpha0) then
      propA = 'PROP_A0'
   elseif(A%ACAlpha==A%ACAlpha1) then
      propA = 'PROP_A1'
   elseif(A%ACAlpha==A%ACAlpha2) then
      propA = 'PROP_A2'
   endif
   call readEVecXY(A%EigX,A%EigY,A%NDimX,propA)
   call readEvalXY(A%Eig,A%NDimX,propA)
 endif
 ! B: get response properties
 if(B%Cubic) then
   allocate(B%EigX(B%NDimX*B%NDimX),B%EigY(B%NDimX*B%NDimX),&
            B%Eig(B%NDimX))

   if(B%ACAlpha==B%ACAlpha0) then
      propB = 'PROP_B0'
   elseif(B%ACAlpha==B%ACAlpha1) then
      propB = 'PROP_B1'
   elseif(B%ACAlpha==B%ACAlpha2) then
      propB = 'PROP_B2'
   endif
   call readEVecXY(B%EigX,B%EigY,B%NDimX,propB)
   call readEvalXY(B%Eig,B%NDimX,propB)
 endif

 ! print thresholds
 if(SAPT%IPrint>5) then
    write(LOUT,'(/,1x,a)') 'Thresholds in E2exch-disp:'
    write(LOUT,'(1x,a,t18,a,e15.4)') 'SmallE','=', SmallE
    write(LOUT,'(1x,a,t18,a,e15.4)') 'BigE',  '=', BigE
 endif

! approximate RDM2
 if(approx) then
    allocate(A%Fmat(NBas,NBas),B%Fmat(NBas,NBas))
    call fill_Fmat(A%Fmat,A%Occ,NBas,1)
    call fill_Fmat(B%Fmat,B%Occ,NBas,1)
 endif

 ! uncoupled - works for CAS only!
 if(uncoupled) then
    allocate(OmA0(A%NDimX),OmB0(B%NDimX))

    call read_SBlock(SBlockA,SBlockAIV,nblkA,'XY0_A')
    call read_SBlock(SBlockB,SBlockBIV,nblkB,'XY0_B')

    call unpack_Eig(SBlockA,SBlockAIV,nblkA,OmA0,A%NDimX)
    call unpack_Eig(SBlockB,SBlockBIV,nblkB,OmB0,B%NDimX)
 endif
 !! uncoupled-ver0
 !allocate(AVecX0(A%NDimX*A%NDimX),OmA0(A%NDimX), &
 !         AVecY0(A%NDimX*A%NDimX), &
 !         BVecX0(B%NDimX*B%NDimX),OmB0(B%NDimX), &
 !         BVecY0(B%NDimX*B%NDimX))
 !
 !call unpack_XY0_full(AVecX0,AVecY0,OmA0,A%CICoef,A%IndN,A%NDimX,NBas,'XY0_A')
 !call unpack_XY0_full(BVecX0,BVecY0,OmB0,B%CICoef,B%IndN,B%NDimX,NBas,'XY0_B')

 allocate(S(NAO,NAO),&
          PA(NAO,NAO),PB(NAO,NAO),&
          Va(NAO,NAO),Vb(NAO,NAO),&
          Sab(NBas,NBas),Sba(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas))
 allocate(tmp1(NAO,NAO),tmp2(NAO,NAO))

 call get_den(NAO,NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NAO,NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NAO)

 call tran_AO2MO2(S,A%CMO,B%CMO,Sab,NAO,NBas)
 call tran_AO2MO2(S,B%CMO,A%CMO,Sba,NAO,NBas)

 call get_one_mat('V',Va,A%Monomer,NAO)
 call get_one_mat('V',Vb,B%Monomer,NAO)

 call tran_AO2MO2(Va,B%CMO,B%CMO,Vabb,NAO,NBas)
 call tran_AO2MO2(Vb,A%CMO,A%CMO,Vbaa,NAO,NBas)

 call tran_AO2MO2(Va,B%CMO,A%CMO,Vaba,NAO,NBas)
 call tran_AO2MO2(Va,A%CMO,B%CMO,Vaab,NAO,NBas)
 call tran_AO2MO2(Vb,A%CMO,B%CMO,Vbab,NAO,NBas)


 allocate(posA(NBas,NBas),posB(NBas,NBas))
 if(Flags%ICASSCF==1) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
    enddo
 elseif(Flags%ICASSCF==0) then
    posA = 0
    do i=1,A%NDimX
       posA(A%IndN(1,i),A%IndN(2,i)) = i
    enddo
    posB = 0
    do i=1,B%NDimX
       posB(B%IndN(1,i),B%IndN(2,i)) = i
    enddo
 endif

 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NAO,NAO,NAO,1d0,PA,NAO,S,NAO,0d0,tmp1,NAO)
 call dgemm('N','N',NAO,NAO,NAO,1d0,S,NAO,PB,NAO,0d0,tmp2,NAO)
 do j=1,NAO
    do i=1,NAO
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 termZ_u = termZ
 termZ_u = 2d0*SAPT%e2disp_unc*termZ_u

 if(A%Cubic.or.B%Cubic) then
    if(A%ACAlpha==A%ACAlpha0.or.B%ACAlpha==B%ACAlpha0) termZ = 2d0*SAPT%e2disp_a0*termZ
    if(A%ACAlpha==A%ACAlpha1.or.B%ACAlpha==B%ACAlpha1) termZ = 2d0*SAPT%e2disp_a1*termZ
    if(A%ACAlpha==A%ACAlpha2.or.B%ACAlpha==B%ACAlpha2) termZ = 2d0*SAPT%e2disp_a2*termZ
 else
    termZ   = 2d0*SAPT%e2disp*termZ
 endif

 !write(LOUT,*) 'termZ ',termZ
 !write(LOUT,'(1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3

 deallocate(tmp2,tmp1)

 allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
          tmp_u(A%NDimX,B%NDimX),sij(A%NDimX,B%NDimX),work(nOVB),workSq(NBas,NBas))
 if(both) allocate(tmp3(A%NDimX,B%NDimX))

 ! term A1
 ! transform J and K
 !if(SAPT%noE2exi) then
 !   call tran4_gen(NBas,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            NBas,A%CMO,&
 !            NBas,B%CMO,&
 !            'FFOOABAB','AOTWOSORT')
 !   endif
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,A%CMO,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         'FOFOABBA','AOTWOSORT')
 ! term A2
 ! A2A(B): XX
 !if(SAPT%noE2exi) then
 !   write(LOUT,'()')
 !   call tran4_gen(NBas,&
 !            NBas,B%CMO,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            NBas,B%CMO,&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            'FOFOBBBA','AOTWOSORT')
 !   call tran4_gen(NBas,&
 !            NBas,A%CMO,&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            NBas,A%CMO,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            'FOFOAAAB','AOTWOSORT')
 !   ! A2A(B): YY
 !   call tran4_gen(NBas,&
 !            NBas,A%CMO,&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            NBas,B%CMO,&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            'FOFOBBAB','AOTWOSORT')
 !   call tran4_gen(NBas,&
 !            NBas,B%CMO,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            NBas,A%CMO,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            'FOFOAABA','AOTWOSORT')
 !endif
 ! term A3
 !if(SAPT%noE2exi) then
 !   call tran4_gen(NBas,&
 !            NBas,B%CMO,&
 !            B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !            NBas,A%CMO,&
 !            A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !            'FOFOAABB','AOTWOSORT')
 !endif
 !! XY and YX, A2
 !write(LOUT,'()')
 !call tran4_gen(NBas,&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
 !         NBas,A%CMO,&
 !         NBas,B%CMO,&
 !         'FFOOABBB','AOTWOSORT')
 !call tran4_gen(NBas,&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
 !         NBas,B%CMO,&
 !         NBas,A%CMO,&
 !         'FFOOBAAA','AOTWOSORT')

 deallocate(work)
 allocate(work(NBas**2),ints(NBas,NBas))

 ! A3: XX
 tmp1=0
 if(approx) then

    call app_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

    call inter_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%RDM2val,B%RDM2val,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                     B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

 ! TERMS XX, YY
 !(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       ir = k
       iq = l

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                do j=1,dimOB
                   do i=1,NBas
                      ints(i,j) = work((j-1)*NBas+i)
                   enddo
                enddo

                fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
                            (B%Occ(ir)-B%Occ(is)) * &
                            ints(ip,is)

                tmp1(ipq,irs) = tmp1(ipq,irs) + fact

             endif
          enddo
       enddo

    enddo
 enddo
 close(iunit)

 sij = tmp1

 ! A1:XX
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XX')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'XX')

 ! A2:XX
 if(approx) then

    call app_A2_XX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat, &
                   Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
                   A%Occ,B%IndN,A%IndN,posB,posA,     &
                   dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_XX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat, &
                   Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
                   B%Occ,A%IndN,B%IndN,posA,posB,     &
                   dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,B%RDM2val,    &
                     Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
                     A%Occ,B%IndN,A%IndN,posB,posA,     &
                     dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,A%RDM2val,    &
                     Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
                     B%Occ,A%IndN,B%IndN,posA,posB,     &
                     dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !X_A.I.X_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,0d0,tmp3,A%NDimX)
 endif

 if(uncoupled) then
    ! UNC
    tmp_u = 0
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'XX')
    !
    ! UNC-test-full
    !!X(0)_A.I.X(0)_B
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecX0,B%NDimX,0d0,tmp_u,A%NDimX)
 endif

  tmp1 = sij

 ! A1:YY
 ! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YY')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'YY')

 ! A2:YY
 if(approx) then

    call app_A2_YY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,B%RDM2val,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
             A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,A%RDM2val,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !Y_A.I.Y_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(uncoupled) then
    ! UNC
    ! Y(0)_A.I.Y(0)_B
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'YY')

    ! UNCOUPLED-test-full
    ! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    ! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecY0,B%NDimX,1d0,tmp_u,A%NDimX)
    !
 endif

! TERMS XY, YX
! A3:XY
 tmp1 = 0
 if(approx) then

    call app_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

   call inter_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%RDM2val,B%RDM2val,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                  B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

! (FF|OO):(AB|AB)
 open(newunit=iunit,file='FFOOABAB',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       iq = k
       is = l

       do j=1,NBas
          do i=1,NBas
             ints(i,j) = work((j-1)*NBas+i)
          enddo
       enddo

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then


                fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
                           (B%Occ(ir)-B%Occ(is)) * &
                           ints(ip,ir)

                tmp1(ipq,irs) = tmp1(ipq,irs) + fact

             endif
            enddo
         enddo

    enddo
 enddo
 close(iunit)

 sij = tmp1

 ! A1:XY
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XY')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab,&
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'XY')

! A2:XY
 if(approx) then

    call app_A2_XY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,B%RDM2val,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,A%RDM2val,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.,B)

 endif

 !X_A.I.Y_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(uncoupled) then
    ! UNC
    !X(0)_A.I.Y(0)_B
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'XY')

    !! UNCOUPLED-test-full
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecX0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecY0,B%NDimX,1d0,tmp_u,A%NDimX)
 endif

 tmp1 = sij

 ! A1:YX
! call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YX')
 call A1_Mod_1el(tmp1,A%XELE,B%XELE,A%Occ,B%Occ,Vaab,Vbab,Sab, &
                 A%IndN,B%IndN,posA,posB,A%NDimX,B%NDimX,NBas,'YX')

 !A2: YX
 if(approx) then

    call app_A2_YX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_XY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,B%RDM2val,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.,B)
    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,A%RDM2val,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !Y_A.I.X_B
 if(both) then
    call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
 endif

 if(uncoupled) then
    ! UNC
    call abpm_tran_gen(tmp1,tmp_u,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                       nblkA,nblkB,A%NDimX,B%NDimX,'YX')

    ! UNC-test-full
    !call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,AVecY0,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
    !call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,BVecX0,B%NDimX,1d0,tmp_u,A%NDimX)
 endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print*, 'sij',norm2(sij)
! print*, 'tmp3-exd',norm2(tmp3)

!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!! SINGLE REFERENCE CODE !!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THERE IS SOME BUG IN A2_XY...
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call exdisp_sref(SAPT%elst,SAPT%Vnn,termZ,sij,Sab,Sba,&
!                 Vaab,Vaba,Vbab,Vabb,Vbaa,posA,posB,NBas,A,B)


 if(uncoupled) then

    ! UNC
    sij = 0
    inquire(file='PROP_AB0',EXIST=ipropab)
    if(ipropab) then
       ! read s_ij
       open(newunit=iunit,file='PROP_AB0',form='UNFORMATTED',&
          access='SEQUENTIAL',status='OLD')
       read(iunit) sij
       close(iunit)

    else

       ! make s_ij uncoupled
       !call make_sij_Y_unc(sij,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
       !                A%num0,B%num0,dimOA,dimOB,nOVB,A%IndN,B%IndN,A%NDimX,B%NDimX,NBas)
       write(LOUT,'(/1x,a)') 'ERROR! make_sij_Yunc not ready yet!'
       write(LOUT,'(1x,a)')  'E2exch-disp(unc) = 0!'
       !stop

    endif

    ! term X(unc)
    termX_u = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then

             termX_u = termX_u + tmp_u(i,j)*sij(i,j)/(OmA0(i)+OmB0(j))

          endif
       enddo
    enddo
    termX_u = -4d0*termX_u

    call make_tij_Y_unc(tmp_u,tmp1,A%Occ,B%Occ,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                        A%IndN,B%IndN,posA,posB,Sab,Sba,nblkA,nblkB,A%NDimX,B%NDimX,NBas)

    ! term Y(unc)
    termY_u = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

          if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
             .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then

             termY_u = termY_u + sij(i,j)*tmp_u(i,j)/(OmA0(i)+OmB0(j))

          endif
       enddo
    enddo

    termY_u = -8d0*(SAPT%elst-SAPT%Vnn)*termY_u

    !if(SAPT%IPrint>2) write(LOUT,'(/1x,a,f16.8)') 'term Z(unc) = ',  termZ_u*1.0d3
    !if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term X(unc) = ',  termX_u*1.0d3
    !if(SAPT%IPrint>2) write(LOUT,'(1x,a,f16.8)')  'term Y(unc) = ',  termY_u*1.0d3

    if(SAPT%IPrint>2) call print_en('term Z(unc)',termZ_u*1.0d3,.true.) 
    if(SAPT%IPrint>2) call print_en('term X(unc)',termX_u*1.0d3,.false.)
    if(SAPT%IPrint>2) call print_en('term Y(unc)',termY_u*1.0d3,.false.)

    e2exd_u = termX_u + termY_u + termZ_u
    SAPT%e2exdisp_unc = e2exd_u
    !write(LOUT,'(/1x,a,f11.8)') 'E2exch-disp(unc) = ', e2exd_u*1.0d3
    call print_en('E2exch-disp(unc)',e2exd_u*1.0d3,.true.)

    ! deallocate SBLOCK
    do iblk=1,nblkA
       associate(A => SBlockA(iblk))
         deallocate(A%matY,A%matX,A%vec)
         deallocate(A%pos)
       end associate
    enddo
    do iblk=1,nblkB
       associate(B => SBlockB(iblk))
         deallocate(B%matY,B%matX,B%vec)
         deallocate(B%pos)
       end associate
    enddo
    ! deallocate IV part
    associate(A => SBlockAIV)
      if(A%n>0) then
         deallocate(A%vec)
         deallocate(A%pos)
      endif
    end associate
    associate(B => SBlockBIV)
      if(B%n>0) then
         deallocate(B%vec)
         deallocate(B%pos)
      endif
    end associate

    deallocate(OmA0,OmB0)
    !deallocate(AVecX0,AVecY0,BVecX0,BVecY0)

 ! end UNC (for CAS only)
 endif

 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(both) then
    ! coupled
    sij = 0
    inquire(file='PROP_AB',EXIST=ipropab)
    if(ipropab) then
       ! read s_ij
       open(newunit=iunit,file='PROP_AB',form='UNFORMATTED',&
          access='SEQUENTIAL',status='OLD')
       read(iunit) sij
       close(iunit)

    else

       ! make s_ij
       if(Flags%ICholesky==0) then
          call make_sij_Y(sij,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
                          A%num0,B%num0,dimOA,dimOB,nOVB,A%IndN,B%IndN, &
                          A%NDimX,B%NDimX,NBas)
       else if(Flags%ICholesky==1) then
          call make_sij_Y_Chol(sij,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
                          A%DChol,B%DChol, &
                          A%num0,B%num0,dimOA,dimOB,nOVB,A%IndN,B%IndN, &
                          A%NDimX,B%NDimX,A%NChol,NBas)
       endif

    endif

    ! term X
    termX = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

      !    if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
      !       .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
             .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then


             termX = termX + tmp3(i,j)*sij(i,j)/(A%Eig(i)+B%Eig(j))

          endif
       enddo
    enddo
    termX = -4d0*termX

    !if(SAPT%IPrint>5) write(LOUT,'(/1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3
    !if(SAPT%IPrint>5) write(LOUT,'(1x,a,f16.8)') 'term X      = ',  termX*1.0d3
   
    if(SAPT%IPrint>5) call print_en('term Z',termZ*1.0d3,.true.)
    if(SAPT%IPrint>5) call print_en('term X',termX*1.0d3,.false.)

    ! termY
    ! term Y: t_ij
    !call make_tij_Y(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
    call make_tij_Y(tmp3,tmp2,tmp1,A%Occ,B%Occ,A%EigY,A%EigX,B%EigY,B%EigX,&
                    A%IndN,B%IndN,posA,posB,Sab,Sba,A%NDimX,B%NDimX,NBas)

    ! term Y
    termY = 0d0
    do j=1,B%NDimX
       do i=1,A%NDimX

         ! if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
         !    .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
             .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then

             !termY = termY + 1d0/(A%Eig(i)+B%Eig(j))
             termY = termY + sij(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))

          endif
       enddo
    enddo

    termY = -8d0*(SAPT%elst-SAPT%Vnn)*termY
    !if(SAPT%IPrint>5) write(LOUT,'(1x,a,f16.8)') 'term Y      = ',  termY*1.0d3
    if(SAPT%IPrint>5) call print_en('term Y',termY*1.0d3,.false.)

    e2exd = termX + termY + termZ

 if(A%Cubic.or.B%Cubic) then

   ! cubic
   if(A%ACAlpha==A%ACAlpha0.or.B%ACAlpha==B%ACAlpha0) SAPT%e2exd_a0 = e2exd
   if(A%ACAlpha==A%ACAlpha1.or.B%ACAlpha==B%ACAlpha1) SAPT%e2exd_a1 = e2exd
   if(A%ACAlpha==A%ACAlpha2.or.B%ACAlpha==B%ACAlpha2) SAPT%e2exd_a2 = e2exd

   !write(LOUT,'(1x,a,f16.8)') 'A: Alpha     = ',A%ACAlpha
   !write(LOUT,'(1x,a,f16.8)') 'B: Alpha     = ',B%ACAlpha
   write(LOUT,'(1x,a,f16.8)') 'E2exd(Alpha)   = ',e2exd*1000d0

 else

    ! regular
    SAPT%e2exdisp     = e2exd
    !write(LOUT,'(/1x,a,f16.8)') 'E2exch-disp = ', e2exd*1.0d3
    call print_en('E2exch-disp',e2exd*1.0d3,.true.)

 endif

 ! end coupled
 endif

 deallocate(sij)

 deallocate(posB,posA,ints)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vabb,Vb,Va,PB,PA,Sab,S)
 deallocate(workSq,work,tmp_u,tmp2,tmp1)
 if(both) deallocate(tmp3)

 ! SAPT cubic
 if(A%Cubic) deallocate(A%Eig,A%EigY,A%EigX)
 if(B%Cubic) deallocate(B%Eig,B%EigY,B%EigX)

end subroutine e2exdisp

subroutine make_tind_unc(tvec,SBlock,SBlockIV,Sab,AOcc,BOcc,AIndN,posA,nblk,NDimX,NBas)
implicit none

integer,intent(in) :: nblk,NBas,NDimX
integer,intent(in) :: AIndN(2,NDimX),posA(NBas,NBas)
double precision,intent(in)    :: Sab(NBas,NBas),&
                                  AOcc(NBas),BOcc(NBas)
type(EBlockData)               :: SBlock(nblk),SBlockIV
double precision,intent(inout) :: tvec(NDimX)

integer :: i,ip,iq,ir,ipq
integer :: ii,j,k
double precision,allocatable :: tmp1(:),tmp2(:)
double precision :: fact,val

tvec = 0

allocate(tmp1(NDimX))
tmp1 = 0
do i=1,NDimX
   ip = AIndN(1,i)
   iq = AIndN(2,i)
   ipq = posA(ip,iq)

   fact = (AOcc(ip)-AOcc(iq))

   val = 0
   do ir=1,NBas
      val = val + BOcc(ir)*Sab(ip,ir)*Sab(iq,ir)
   enddo
   tmp1(ipq) = tmp1(ipq) + fact*val

enddo

call abpm_dgemv_gen(tmp1,tvec,SBlock,SBlockIV,nblk,NDimX,'YX')

deallocate(tmp1)

end subroutine make_tind_unc

subroutine make_tij_Y_unc(tmp2,tmp1,AOcc,BOcc,SBlockA,SBlockAIV,SBlockB,SBlockBIV,&
                          AIndN,BIndN,posA,posB,Sab,Sba,nblkA,nblkB,ANDimX,BNDimX,NBas)
implicit none

type(EBlockData)               :: SBlockA(nblkA),SBlockAIV,&
                                  SBlockB(nblkB),SBlockBIV
integer,intent(in)             :: nblkA,nblkB,ANDimX,BNDimX,NBas
integer,intent(in)             :: posA(NBas,NBas),posB(NBas,NBas),&
                                  AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in)    :: Sab(NBas,NBas),Sba(NBas,NBas),&
                                  AOcc(NBas),BOcc(NBas)
double precision,intent(inout) :: tmp1(ANDimX,BNDimX),&
                                  tmp2(ANDimX,BNDimX)

integer :: i,j,ir,is,irs,ip,iq,ipq
double precision :: fact

tmp1=0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) + Sab(iq,ir)*Sba(is,ip)*fact

   enddo
enddo
tmp2 = 0
!X_A.I.X_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'XX')
!Y_A.I.Y_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'YY')

tmp1 = 0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) - Sab(ip,ir)*Sba(is,iq)*fact

   enddo
enddo
!X_A.I.Y_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'XY')
!Y_A.I.X_B
call abpm_tran_gen(tmp1,tmp2,SBlockA,SBlockAIV,SBlockB,SBlockBIV, &
                   nblkA,nblkB,ANDimX,BNDimX,'YX')

end subroutine make_tij_Y_unc

end module sapt_exch
