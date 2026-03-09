module sapt_exch_e1
use types
use tran
use exmisc
use exi
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
 !


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

subroutine e1exch_dmft(Flags,A,B,SAPT)
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
 case(1,11)
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

 ! n^A n^B Sab Sab
 nnS2 = 0
 do j=1,dimOB
 do i=1,dimOA
    nnS2 = nnS2 + A%Occ(i)*B%Occ(j)*Sab(i,j)**2
 enddo
 enddo

 ! this should = -tmpDEL
 tmpELST = 2d0*(SAPT%elst-SAPT%Vnn)*nnS2
 print*, 'tmpELST',tmpELST*1000

 ! tvk = n_p n_q (v^A S + v^B S + v_pq^qp)
 open(newunit=iunit,file='FFOOABAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)

 ipq = 0
 tvk = 0
 do iq=1,dimOB
    do ip=1,dimOA
       ipq = ipq + 1
       read(iunit,rec=ipq) ints(1:NBas*NBas)

       tvk(3) = tvk(3) + A%Occ(ip)*B%Occ(iq)*ints(ip+(iq-1)*NBas)

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

 tmpDEL = 0
 do i=1,dimOA
    tmpDEL = tmpDEL + A%Occ(i)*Vbaa(i,i)
 enddo
 do i=1,dimOB
    tmpDEL = tmpDEL + B%Occ(i)*Vabb(i,i)
 enddo
 tmpDEL = -4d0*tmpDEL*nnS2

 ! tNa
 tNa = 0
 do iq=1,dimOA
    do ip=1,dimOA
       tNa(1) = tNa(1) + AlphaA(ip)*AlphaA(iq)*Saa(ip,iq)*Vbaa(ip,iq)
    enddo
 enddo
 tNa(1) = 2d0*tNa(1)

!(FO|FO): (AA|AB)
open(newunit=iunit,file='FOFOAAAB',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
val  = 0
do it=1,dimOB
   do iq=1,dimOA
      read(iunit,rec=iq+(it-1)*NBas) ints(1:NBas*dimOA)

      fac = A%Occ(iq)*B%Occ(it)*Sab(iq,it)
      val = 0
      do ip=1,dimOA
         val = val + A%Occ(ip)*ints(ip+(ip-1)*NBas)
      enddo
      tNa(2) = tNA(2) - 4d0*fac*val

      fac = B%Occ(it)*AlphaA(iq)
      val = 0
      do ip=1,dimOA
         val = val + AlphaA(ip)*Sab(ip,it)*ints(ip+(iq-1)*NBas)
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
open(newunit=iunit,file='FOFOBBBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
ints = 0
val  = 0
do it=1,dimOA
   do iq=1,dimOB
      read(iunit,rec=iq+(it-1)*NBas) ints(1:NBas*dimOB)

      fac = A%Occ(it)*B%Occ(iq)*Sab(it,iq)
      val = 0
      do ip=1,dimOB
         val = val + B%Occ(ip)*ints(ip+(ip-1)*NBas)
      enddo
      tNb(2) = tNb(2) - 4d0*fac*val

      fac = A%Occ(it)*AlphaB(iq)
      val = 0
      do ip=1,dimOB
         val = val + AlphaB(ip)*Sab(it,ip)*ints(ip+(iq-1)*NBas)
      enddo
      tNb(3) = tNb(3) + 2d0*fac*val

   enddo
enddo

print*, 'tNb-2',tNb(2)*1000
print*, 'tNb-3',tNb(3)*1000

close(iunit)

 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)

 work  = 0
 tmp   = 0
 tNaNb = 0
 do iq=1,dimOA
    do ip=1,dimOA
      read(iunit,rec=ip+(iq-1)*dimOA) work(1:dimOB,1:dimOB)

      if(ip==iq) then

         val = 0
         do it=1,dimOB
            val = val + B%Occ(it)*work(it,it)
         enddo
         val = A%Occ(ip)*val
         tmpDEL = tmpDEL - 8d0*nnS2*val
         tmp = tmp - 8d0*nnS2*val

         val = 0
         do iu=1,dimOB
            do it=1,dimOB
               val = val + AlphaB(it)*AlphaB(iu)*Sbb(it,iu)*work(it,iu)
            enddo
         enddo
         tNaNb(1) = tNaNb(1) - 2d0*val*A%Occ(ip)

      endif

      val = 0
      do it=1,dimOB
         val = val + B%Occ(it)*work(it,it)
      enddo
      tNaNb(2) = tNaNb(2) - 2d0*val*AlphaA(ip)*AlphaA(iq)*Saa(ip,iq)

      val = 0
      do iu=1,dimOB
         do it=1,dimOB
            val = val + AlphaB(it)*AlphaB(iu)*Sab(ip,iu)*Sab(iq,it)*work(it,iu)
         enddo
      enddo
      tNaNb(3) = tNaNb(3) + AlphaA(ip)*AlphaA(iq)*val

    enddo
 enddo
 tNaNB = -2d0*tNaNb
 print*, 'A4',sum(tNaNB)*1000
 print*, 'tmp-A4',(tmp+sum(tNaNB))*1000
 print*, 'tNANB-1',tNaNB(1)*1000
 print*, 'tNANB-2',tNaNB(2)*1000
 print*, 'tNANB-3',tNaNB(3)*1000

 close(iunit)

 e1ex_dmft = sum(tvk)+sum(tNa)+sum(tNb)+sum(tNaNb)
 SAPT%exchs2 = e1ex_dmft

 if(SAPT%SaptExch==0) then
    call print_en('E1exch-DMFT(S2)',e1ex_dmft*1000,.true.)
 elseif(SAPT%SaptExch==1) then
    if(Flags%IRdm2Typ==0) then
       call print_en('E1exch-DMFT(nn)',e1ex_dmft*1000,.true.)
    elseif(Flags%IRDM2typ==1.or.Flags%IRDM2Typ==11) then
       call print_en('E1exch-DMFT(BB)',e1ex_dmft*1000,.true.)
    endif
 endif

 deallocate(Vbaa,Vabb,Vbba,Vaab)
 deallocate(AlphaB,AlphaA)
 deallocate(ints,work)
 deallocate(Sbb,Saa,Sab)

end subroutine e1exch_dmft

subroutine e1exch_dmft_2(Flags,A,B,SAPT)
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
double precision :: tmpELST,tmpDEL
double precision :: e1ex_dmft
double precision :: tvk(3),tNa(3),tNb(3),tNaNb(3)
double precision,allocatable :: Va(:,:),Vb(:,:),S(:,:)
double precision,allocatable :: Wbb(:,:)
double precision,allocatable :: Saa(:,:),Sbb(:,:),Sab(:,:),Sabh(:,:)
double precision,allocatable :: Vaab(:,:),Vbba(:,:),Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: AlphaA(:),AlphaB(:)
double precision,allocatable :: work(:,:),ints(:)
double precision,external :: ddot

! set dimensions
NBas  = A%NBasis
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
case(1,11)
   ! BB functional
   do i=1,dimOA
      AlphaA(i) = sqrt(A%Occ(i))
   enddo
   do j=1,dimOB
      AlphaB(j) = sqrt(B%Occ(j))
   enddo
case(3)
  write(LOUT,*) 'POWER FUNCITONAL NOT READY YET!'
  stop
end select

allocate(S(NBas,NBas),Sab(NBas,NBas),Sabh(dimOA,dimOB))

call get_one_mat('S',S,A%Monomer,NBas)
call tran2MO(S,A%CMO,B%CMO,Sab,NBas)

Sabh = 0d0
do ip=1,dimOA
   do iu=1,dimOB
      Sabh(ip,iu) = Sabh(ip,iu) + Sab(ip,iu)*AlphaA(ip)*AlphaB(iu)
   enddo
enddo

print*, norm2(Sabh)

allocate(work(dimOB,dimOB))
allocate(Wbb(dimOB,dimOB))

open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)

work  = 0d0
tNaNb = 0d0
do iq=1,dimOA
   do ip=1,dimOA
      read(iunit,rec=ip+(iq-1)*dimOA) work(1:dimOB,1:dimOB)

      Wbb = 0d0
      do iu=1,dimOB
         do it=1,dimOB
            Wbb(it,iu) = Sabh(ip,iu)*Sabh(iq,it)
         enddo
      enddo
      tNaNb(3) = tNaNb(3) + ddot(dimOB**2,work,1,Wbb,1)

   enddo
enddo
tNaNB = -2d0*tNaNb

close(iunit)

print*, '@@@ tNANB-3',tNaNB(3)*1000

deallocate(work,Wbb)
deallocate(AlphaB,AlphaA)
deallocate(Sabh,Sab,S)

end subroutine e1exch_dmft_2

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

allocate(Waa(NBasis,NBasis),Wab(NBasis,NBasis))
allocate(Wba(NBasis,NBasis),Wbb(NBasis,NBasis))

call tran2MO(A%WPot,B%UMO(:,:,1),B%UMO(:,:,1),Waa,NBasis)
call tran2MO(A%WPot,B%UMO(:,:,2),B%UMO(:,:,2),Wab,NBasis)

call tran2MO(B%WPot,A%UMO(:,:,1),A%UMO(:,:,1),Wba,NBasis)
call tran2MO(B%WPot,A%UMO(:,:,2),A%UMO(:,:,2),Wbb,NBasis)

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
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T1-b ) = ', ex1(2)*1000d0

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

if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'ExchS2(T3-b ) = ', ex3(2)*1000d0

e1exs2 = sum(ex1) + sum(ex2) + sum(ex3)
e1exs2 = - e1exs2
SAPT%exchs2 = e1exs2

call print_en('E1exch(S2)',e1exs2*1000,.true.)

deallocate(Wbb,Wba,Wab,Waa)
deallocate(Sb,Sa,Sbt,Sat)
deallocate(Aux,work)

end subroutine e1exchs2_sq_os

end module sapt_exch_e1
