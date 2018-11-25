module sapt_ener
use types
use tran

implicit none

contains

subroutine e1elst(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j
integer :: NBas
double precision,allocatable :: PA(:,:),PB(:,:) 
double precision,allocatable :: Va(:,:),Vb(:,:),Ja(:,:) 
double precision,allocatable :: work(:,:)
double precision :: tmp,ea,eb,elst
double precision,parameter :: Half=0.5d0
double precision,external  :: trace

! set dimensions
 NBas = A%NBasis 

 allocate(PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),Ja(NBas,NBas))
 allocate(work(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,2d0,PA)
 call get_den(NBas,B%CMO,B%Occ,2d0,PB)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call make_J1(NBas,PA,Ja,'AOTWOSORT')

! Tr[Pa.Va + Pb.Vb + Pb.Ja]
 work=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Vb,NBas,0d0,work,NBas)
 ea = trace(work,NBas)
! print*, ea
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Va,NBas,0d0,work,NBas)
 eb = trace(work,NBas) 
! print*, eb
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,Ja,NBas,0d0,work,NBas)
 ea = ea + trace(work,NBas)
! print*, trace(work,NBas) 
 elst = ea + eb + SAPT%Vnn 

 write(LOUT,'(1x,a,f16.8)') 'V_nn        = ', SAPT%Vnn
 write(LOUT,'(1x,a,f16.8)') 'Eelst       = ', elst*1000d0 
 SAPT%elst = elst

!! test wabb
! work = 0
! Ja = 0
! call make_J1(NBas,PB,Ja)
!
! call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,Ja,NBas,0d0,work,NBas)
! print*, 'Pa.Jb', trace(work,NBas)
!
! call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,B%WPot,NBas,0d0,work,NBas)
! print*, 'Pa.Wb', trace(work,NBas)
!

 deallocate(work)
 deallocate(Ja,Vb,Va,PB,PA) 

end subroutine e1elst

subroutine e1exchs2(A,B,SAPT)
implicit none

type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: i, j, k, l, ia, jb
integer :: ip, iq, ir, is
integer :: iv, iz, iu, it
integer :: iunit
integer :: dimOA,dimOB,NBas
integer :: GdimOA,GdimOB
double precision,allocatable :: S(:,:),Sab(:,:)
double precision,allocatable :: USa(:,:),USb(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:) 
double precision,allocatable :: Kb(:,:)
double precision,allocatable :: JJb(:,:)
double precision,allocatable :: Qab(:,:),Qba(:,:) 
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:)
double precision,allocatable :: tmpA(:,:,:,:),tmpB(:,:,:,:), &
                                tmpAB(:,:,:,:)
double precision,allocatable :: work(:,:),RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision :: tmp,ea,eb,exchs2
double precision :: t1(2),t2a(4),t2b(2),t2c(2),t2d
double precision :: t1f,t2f
double precision,parameter :: Half=0.5d0
double precision,external  :: trace,FRDM2

! set dimensions
 NBas = A%NBasis 
 dimOA = A%INAct+A%NAct 
 GdimOA = A%num0+A%num1
 dimOB = B%INAct+B%NAct 
 GdimOB = B%num0+B%num1

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          PAbb(NBas,NBas),PBaa(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas))
 allocate(USa(NBas,NBas),USb(NBas,NBas),&
          Qab(NBas,NBas),Qba(NBas,NBas),&
          Kb(NBas,NBas))
 allocate(work(NBas,NBas),tmp1(NBas,NBas),tmp2(NBas,NBas))
! allocate(JJb(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas) 

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas) 
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas) 

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))
 do l=1,dimOA
    do k=1,dimOA 
       do j=1,dimOA
          do i=1,dimOA
             RDM2Aval(i,j,k,l) = FRDM2(i,k,j,l,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
          enddo
       enddo
    enddo
 enddo
 do l=1,dimOB
    do k=1,dimOB 
       do j=1,dimOB
          do i=1,dimOB
             RDM2Bval(i,j,k,l) = FRDM2(i,k,j,l,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
          enddo
       enddo
    enddo
 enddo

! USa=0; USb=0
! USa,USb in MOAO
! call dgemm('T','N',NBas,NBas,NBas,1d0,A%CMO,NBas,S,NBas,0d0,USa,NBas)
! call dgemm('T','N',NBas,NBas,NBas,1d0,B%CMO,NBas,S,NBas,0d0,USb,NBas)
! USa,USb in AOMO
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,A%CMO,NBas,0d0,USa,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,B%CMO,NBas,0d0,USb,NBas)

! PA(B), PB(A)
 call tran2MO(PA,USb,USb,PAbb,NBas) 
 call tran2MO(PB,USa,USa,PBaa,NBas) 

! Qab=0; Qba=0
 !call dgemm('N','T',NBas,NBas,NBas,1d0,PA,NBas,USb,NBas,0d0,Qab,NBas)
 !call dgemm('N','T',NBas,NBas,NBas,1d0,PB,NBas,USa,NBas,0d0,Qba,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,USb,NBas,0d0,Qab,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,USa,NBas,0d0,Qba,NBas)

 call tran3Q_full(NBas,A%CMO,Qba,'TWOA3B')
 !call tran3Q_full(dimOA,A%CMO(1:dimOA,1:dimOA),Qba(1:dimOA,1:dimOA),'TWOA3B')
 call tran3Q_full(NBas,B%CMO,Qab,'TWOB3A')
! HERE!!!
 ! Design tran3Q_lim

 call make_K(NBas,PB,Kb)

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
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       t1(2) = t1(2) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo
! work=0
! call dgemm('T','N',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'd',trace(work,NBas)
! print*, t1(2)
! t1(2) = trace(work,NBas)

 t1f = 2d0*t1(1)*t1(2) 
 write(LOUT,*) 'T1 ',t1f

! T2d
 t2d = -2d0*SAPT%Vnn*t1(2)
 write(LOUT,*) 'T2d',t2d

! T2c
 t2c=0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,Va,NBas,PB,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp2,NBas)
  do j=1,NBas
    do i=1,NBas
       t2c(1) = t2c(1) + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo
 t2c(1) = -2d0*t2c(1)
 write(LOUT,*) 
 print*, 'T2c(1)',t2c(1) 
! !test
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:NT', -2d0*trace(work,NBas)  
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:TN', -2d0*trace(work,NBas)  

! Full NBas
! do ir=1,NBas
!    do ip=1,NBas
!       do is=1,NBas
!          do iq=1,NBas
!             t2c(2) = t2c(2) + FRDM2(ip,iq,ir,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas) * &
!                    Vabb(ip,ir)*PAbb(is,iq)
!          enddo
!       enddo    
!    enddo
! enddo
! t2c(2) = -2d0*t2c(2)
! write(LOUT,*)
! write(LOUT,*) 'T2c(2) ',t2c(2)

! dimOB
! old
! t2c(2)=0
! do ir=1,dimOB
!    do ip=1,dimOB
!       do is=1,dimOB
!          do iq=1,dimOB
!             t2c(2) = t2c(2) + FRDM2(ip,iq,ir,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas) * &
!                    Vabb(ip,ir)*PAbb(is,iq)
!          enddo
!       enddo    
!    enddo
! enddo
! new
 t2c(2)=0
 do is=1,dimOB
    do iq=1,dimOB
       t2c(2) = t2c(2) + sum(RDM2Bval(:,:,iq,is)*Vabb(:,:)*PAbb(is,iq))
    enddo
 enddo    
 t2c(2) = -2d0*t2c(2)
 write(LOUT,*) 'T2c(2) ',t2c(2)
! 
! T2b
 t2b=0
 tmp1=0
 tmp2=0
 call dgemm('N','N',NBas,NBas,NBas,1d0,Vb,NBas,PA,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,PB,NBas,S,NBas,0d0,tmp2,NBas)
  do j=1,NBas
    do i=1,NBas
       t2b(1) = t2b(1) + tmp1(i,j)*tmp2(i,j)
    enddo
  enddo
 t2b(1) = -2d0*t2b(1)
 write(*,*)
 print*, 'T2b(1)',t2b(1)
! !test
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:NT', -2d0*trace(work,NBas)  
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:TN', -2d0*trace(work,NBas)  

! Full NBas 
! do ir=1,NBas
!    do ip=1,NBas
!       do is=1,NBas
!          do iq=1,NBas
!             t2b(2) = t2b(2) + FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas) * &
!                    Vbaa(ip,ir)*PBaa(is,iq)
!          enddo
!       enddo    
!    enddo
! enddo
! t2b(2) = -2d0*t2b(2)
! write(LOUT,*) 'T2b(2) ',t2b(2)
!
! dimOA
! old:
! t2b(2)=0
! do ir=1,dimOA
!    do ip=1,dimOA
!       do is=1,dimOA
!          do iq=1,dimOA
!             t2b(2) = t2b(2) + FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas) * &
!                    Vbaa(ip,ir)*PBaa(is,iq)
!          enddo
!       enddo    
!    enddo
! enddo
!
!new:
 t2b(2)=0
 do is=1,dimOA
    do iq=1,dimOA
       t2b(2) = t2b(2) + sum(RDM2Aval(:,:,iq,is)*Vbaa(:,:)*PBaa(is,iq))
    enddo
 enddo    
 t2b(2) = -2d0*t2b(2)
 write(LOUT,*) 'T2b(2) ',t2b(2)

!! no FRDM2 T2b, T2c
! block
! double precision :: tmpT2c(3),tmpT2b(3)
!
! tmpT2c = 0
! tmp=0
! do ip=1,B%NAct
!    do iq=1,B%NAct
!       do ir=1,B%NAct
!          do is=1,B%NAct
!             tmpT2c(1) = tmpT2c(1) + B%RDM2Act(ip,iq,ir,is) * &
!                    Vabb(B%INAct+ip,B%INAct+ir)*PAbb(B%INAct+is,B%INAct+iq)
!          enddo
!       enddo    
!    enddo
! enddo
! tmpT2c(1) = -2d0*tmpT2c(1)
!! write(LOUT,*) 'test(2) ',-2d0*tmp
!!
!! inact-inact + inact-act
! do ip=1,B%INAct+B%NAct
!    do iq=1,B%INAct+B%NAct
!       if(ip.gt.B%INAct.and.iq.gt.B%INAct) cycle
!         tmpT2c(2) = tmpT2c(2) + B%Occ(ip)*B%Occ(iq)*Vabb(ip,ip)*PAbb(iq,iq)
!    enddo
! enddo
! tmpT2c(2) = -4d0*tmpT2c(2)
!! print*, 'test(2-1)', tmp
! do ip=1,B%INAct+B%NAct 
!    do iq=1,B%INAct+B%NAct 
!       if(ip.gt.B%INAct.and.iq.gt.B%INAct) cycle
!         tmpT2c(3) = tmpT2c(3) + B%Occ(ip)*B%Occ(iq)*Vabb(ip,iq)*PAbb(ip,iq)
!    enddo
! enddo
! tmpT2c(3) = 2d0*tmpT2c(3)
! print*, 'tmpT2c',sum(tmpT2c)
!
! tmpT2b = 0
! do ip=1,A%NAct
!    do iq=1,A%NAct
!       do ir=1,A%NAct
!          do is=1,A%NAct
!             tmpT2b(1) = tmpT2b(1) + A%RDM2Act(ip,iq,ir,is) * &
!                    Vbaa(A%INAct+ip,A%INAct+ir)*PBaa(A%INAct+is,A%INAct+iq)
!          enddo
!       enddo    
!    enddo
! enddo
! tmpT2b(1) = -2d0*tmpT2b(1)
!! inact-inact + inact-act
! do ip=1,A%INAct+A%NAct
!    do iq=1,A%INAct+A%NAct
!       if(ip.gt.A%INAct.and.iq.gt.A%INAct) cycle
!         tmpT2b(2) = tmpT2b(2) + A%Occ(ip)*A%Occ(iq)*Vbaa(ip,ip)*PBaa(iq,iq)
!    enddo
! enddo
! tmpT2b(2) = -4d0*tmpT2b(2)
! do ip=1,A%INAct+A%NAct 
!    do iq=1,A%INAct+A%NAct 
!       if(ip.gt.A%INAct.and.iq.gt.A%INAct) cycle
!         tmpT2b(3) = tmpT2b(3) + A%Occ(ip)*A%Occ(iq)*Vbaa(ip,iq)*PBaa(ip,iq)
!    enddo
! enddo
! tmpT2b(3) = 2d0*tmpT2b(3)
! print*, 'tmpT2b',sum(tmpT2b)
!
! end block

! T2a 
 t2a=0
 do jb=1,NBas
    do ia=1,NBas
       t2a(1) = t2a(1) + PA(ia,jb)*Kb(jb,ia)
    enddo
 enddo
 t2a(1) = -2.0d0*t2a(1)
 write(LOUT,*),'T2a(1)',t2a(1)

 open(newunit=iunit,file='TWOA3B',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)

! Qba 
 do ir=1,NBas
    do ip=1,ir
       read(iunit,rec=ip+ir*(ir-1)/2) work

       if(ip==ir) then
 
         do is=1,dimOA !NBas
             do iq=1,dimOA !NBas
                  t2a(2) = t2a(2) + work(iq,is)* &
                           FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
             enddo
          enddo

       else
         
          do is=1,dimOA !NBas
             do iq=1,dimOA !NBas
                t2a(2) = t2a(2) + work(iq,is)* &
                        (FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)+ &
                         FRDM2(ir,iq,ip,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas))
             enddo
          enddo
  
       endif

    enddo
 enddo
 t2a(2) = -2*t2a(2)
 print*, 'T2a(2) ',t2a(2)
 close(iunit)

! Qab
 open(newunit=iunit,file='TWOB3A',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)

 do ir=1,NBas
    do ip=1,ir
       read(iunit,rec=ip+ir*(ir-1)/2) work

       if(ip==ir) then
 
         do is=1,dimOB !NBas
             do iq=1,dimOB !NBas
                  t2a(3) = t2a(3) + work(iq,is)* &
                           FRDM2(ip,iq,ir,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
             enddo
          enddo

       else
         
          do is=1,dimOB !NBas
             do iq=1,dimOB !NBas
                t2a(3) = t2a(3) + work(iq,is)* &
                        (FRDM2(ip,iq,ir,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)+ &
                         FRDM2(ir,iq,ip,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas))
             enddo
          enddo
  
       endif

    enddo
 enddo
 t2a(3) = -2*t2a(3)
 print*, 'T2a(3) ',t2a(3)
 close(iunit)

! T2a(4)
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
 call dgemm('N','N',dimOA**3,dimOB,dimOA,1d0,RDM2Aval,dimOA**3,Sab,NBas,0d0,tmpA,dimOA**3)

! N^5
!
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
    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,RDM2Bval(:,:,:,is),dimOB**2,Sab,NBas,0d0,tmpB(:,:,:,is),dimOB**2)
 enddo

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
 call dgemm('N','T',dimOA**2,dimOB**2,dimOA*dimOB,1d0,tmpA,dimOA**2,tmpB,dimOB**2,0d0,tmpAB,dimOA**2)

 work=0
 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*GdimOB**2)

 do ir=1,A%INAct+A%NAct
    do ip=1,A%INAct+A%NAct
     !  print*, ip,ir,ip+(ir-1)*GdimOA 
      read(iunit,rec=ip+(ir-1)*GdimOA) work(1:GdimOB,1:GdimOB)

       do iv=1,B%INAct+B%NAct
          do iu=1,B%INAct+B%NAct
             t2a(4)=t2a(4)+work(iv,iu)* &
                    tmpAB(ip,ir,iv,iu)
          enddo
       enddo

    enddo
 enddo
 close(iunit)
 t2a(4) = -2*t2a(4)
 print*, 't2a(4): ',t2a(4)

! ! TESTYTESTYTESTYTESTY
! ! active-active t2a4?
! block
! integer :: iunit2
! integer :: pq,rs,ip,iq,ir,is
! integer :: dimOA,dimVA,nOVA
! integer :: dimOB,dimVB,nOVB
! double precision,allocatable :: work2(:),work3(:,:)
! double precision,allocatable :: ActDirA(:,:), ActDirB(:,:)
! double precision,allocatable :: ActIndA(:,:,:,:), ActIndB(:,:,:,:)
! double precision :: tst
!
! ! set dimensions
! dimOA = A%num0+A%num1
! dimVA = A%num1+A%num2
! dimOB = B%num0+B%num1
! dimVB = B%num1+B%num2
! nOVA = dimOA*dimVA
! nOVB = dimOB*dimVB
!
!! N^5 
! tmpA = 0
! do iz=1,A%NAct
!    do ir=1,A%NAct
!       do iq=1,A%NAct
!          do ip=1,A%NAct
!             do is=1,A%NAct
!                 tmpA(ip,iq,ir,iz) = tmpA(ip,iq,ir,iz) + &
!                                  Sab(A%INAct+is,B%INAct+iz)* &
!                                  A%RDM2Act(ip,iq,ir,is)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!! N^5
!tmpB = 0
! do iq=1,B%NAct
!    do iu=1,B%NAct
!       do iz=1,B%NAct
!          do iv=1,B%NAct
!             do it=1,B%NAct
!                 tmpB(iv,iz,iu,iq) = tmpB(iv,iz,iu,iq) + &
!                                  Sab(A%INAct+iq,B%INAct+it)* &
!                                  B%RDM2Act(iv,iz,iu,it)
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
!!!! nowe:
!!allocate(ActIndA(dimOA,dimOB,dimOA,dimOB),ActDirA(dimOA,dimOA),&
!!         ActIndB(dimOB,dimOB,dimOB,dimOB),ActDirB(dimOB,dimOB))
!!
!!ActDirA = 0
!!! act-direct
!!do ip=1,A%NAct
!!   do iq=1,A%NAct
!!      do ir=1,A%NAct
!!         do iz=1,B%NAct
!!
!!            ActDirA(ip,ir) = ActDirA(ip,ir) + &
!!                             B%Occ(B%INAct+iz)* &
!!                             Sab(A%INAct+iq,B%INAct+iz)*tmpA(ip,iq,ir,iz)
!!         enddo
!!      enddo
!!   enddo
!!enddo
!!
!!! act-indirect
!!! N^5
!!ActIndA = 0
!!do ip=1,A%NAct
!!   do iq=1,A%NAct
!!      do ir=1,A%NAct
!!         do iz=1,B%NAct
!!            do iv=1,B%NAct  
!!               ActIndA(ip,iq,ir,iv) = ActIndA(ip,iq,ir,iv) + &
!!                         tmpA(ip,iq,ir,iz)*Sab(A%INAct+iq,B%INAct+iv)
!!            enddo
!!         enddo
!!       enddo
!!   enddo
!!enddo
!
!! N^6
! tmpAB=0
! do iu=1,B%NAct
!    do iv=1,B%NAct 
!       do ir=1,A%NAct
!          do ip=1,A%NAct
!             do iz=1,B%NAct 
!                do iq=1,A%NAct
!                   tmpAB(ip,ir,iv,iu) = tmpAB(ip,ir,iv,iu) + &
!                                        tmpA(ip,iq,ir,iz)*tmpB(iv,iz,iu,iq)
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! enddo
!
! allocate(work3(GdimOB,GdimOB))
!
! open(newunit=iunit2,file='TMPOOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*GdimOB**2)
!
! work3=0
! tst=0
! do ir=1,A%INAct+A%NAct
!    do ip=1,A%INAct+A%NAct
!     !  print*, ip,ir,ip+(ir-1)*GdimOA 
!      read(iunit2,rec=ip+(ir-1)*GdimOA) work3(1:GdimOB,1:GdimOB)
!       !zle zrobione!!! 
!       if(ir.gt.A%INAct.and.ip.gt.B%INAct) then
!          do iv=1,B%NAct
!             do iu=1,B%NAct
!                tst = tst + work3(B%INAct+iv,B%INAct+iu)* &
!                      tmpAB(ip-A%INAct,ir-A%INAct,iv,iu)
!             enddo
!          enddo
!       endif
!
!    enddo
! enddo
! print*, 'test: ',-2*tst
!
!! allocate(work2(nOVB))
!! open(newunit=iunit2,file='TWOMOAB',status='OLD',&
!!      access='DIRECT',form='UNFORMATTED',recl=8*nOVB)
!!
!! print*, 'test! Gamma(A)-Gamma(B)'
!! work2=0
!! tst=0
!! do ip=1,A%NAct
!!    do iq=1,A%NAct
!!   ! print*, iq,ip
!!    read(iunit2,rec=(iq+A%INAct)+(ip-1)*dimOA) work2(1:nOVB)
!!
!!    do ir=1,B%NAct
!!       do is=1,B%NAct
!!          tst = tst + tmpAB(ip,iq,ir,is)*&
!!                work2((is+B%INAct)+(ir-1)*dimOB)
!!       enddo
!!    enddo 
!!
!!    enddo
!! enddo
!! print*, 'test: ',-2*tst
!
!! deallocate(work2)
! deallocate(work3)
!! deallocate(ActIndA,ActIndB,ActDirB,ActDirA)
! close(iunit2)
!
! end block


 ! only direct part of RDM2

! ! T2b(2)
! tmp=0
! do i=1,NBas
! do j=1,NBas
!    tmp = tmp + PA(i,j)*Vb(j,i)
! enddo
! enddo 
! write(*,*) 
!! print*, 'Test: ',-4d0*tmp*t1(2)
! print*, 'T2b(2):',t2b(2)+4d0*tmp*t1(2)
! ! T2c(2)
! tmp=0
! do i=1,NBas
! do j=1,NBas
!    tmp = tmp + PB(i,j)*Va(j,i)
! enddo
! enddo 
! print*, 'T2c(2):',t2c(2)+4d0*tmp*t1(2)
!
! Last term: 
! call make_J1(NBas,PB,JJb)
! tmp=0
! do jb=1,NBas
!    do ia=1,NBas
!       tmp = tmp + PA(ia,jb)*JJb(jb,ia)
!    enddo
! enddo
! tmp = -8d0*tmp*t1(2)
! print*, 'Last term:',tmp-t2a(4)
!
!! check direct-(in)direct terms:
! block
! double precision :: PAi(NBas,NBas),PAa(NBas,NBas)
! double precision :: PBi(NBas,NBas),PBa(NBas,NBas)
! double precision :: JAi(NBas,NBas),JAa(NBas,NBas)
! double precision :: JBi(NBas,NBas),JBa(NBas,NBas)
! double precision :: Kasb(NBas,NBas)
! double precision,allocatable :: ASB(:,:)
! double precision :: ddaa,ddii,ddia,ddai
! double precision :: Taa,Tii,Tia,Tai
! double precision :: tst(4)
!
!! inactive dens
! PAi=0
! PBi=0
! do i = 1,A%INAct
!    call dger(NBas, NBas, 1d0*A%Occ(i), A%CMO(:, i), 1, A%CMO(:, i), 1, PAi, NBas)
! enddo
! do i = 1,B%INAct
!    call dger(NBas, NBas, 1d0*B%Occ(i), B%CMO(:, i), 1, B%CMO(:, i), 1, PBi, NBas)
! enddo
!
!! active dens
! PAa=0
! PBa=0
! PAa=PA-PAi
! PBa=PB-PBi
!
!! inactive Coulomb
! JAi=0
! JBi=0
! call make_J1(NBas,PAi,JAi)
! call make_J1(NBas,PBi,JBi)
!
!! active Coulomb
! JBa=0
! call make_J1(NBas,PBa,JBa)
!
!
! tst=0
!! direct-direct
! ddii=0
! ddia=0
! ddai=0
! ddaa=0
! do jb=1,NBas
!    do ia=1,NBas
!       ddii = ddii + PAi(ia,jb)*JBi(jb,ia)
!       ddia = ddia + PAi(ia,jb)*JBa(jb,ia)
!       ddai = ddai + PAa(ia,jb)*JBi(jb,ia)
!       ddaa = ddaa + PAa(ia,jb)*JBa(jb,ia)
!    enddo
! enddo
! print*, ddii,ddia,ddai,ddaa
!
! Tii=0
! do iq=1,A%INAct
!    do it=1,B%INAct
!       Tii = Tii + Sab(iq,it)**2*A%Occ(iq)*B%Occ(it)
!    enddo
! enddo
! !tst(1) = -8*tst(1)*ddii
! Tia=0
! do iq=1,A%INAct
!    do it=1,B%NAct
!       Tia = Tia + Sab(iq,B%INAct+it)**2*A%Occ(iq)*B%Occ(B%INAct+it)
!    enddo
! enddo
! Tai=0
! do iq=1,A%NAct
!    do it=1,B%INAct
!       Tai = Tai + Sab(A%INAct+iq,it)**2*A%Occ(A%INAct+iq)*B%Occ(it)
!    enddo
! enddo
! Taa=0
! do iq=1,A%NAct
!    do it=1,B%NAct
!       Taa = Taa + Sab(A%INAct+iq,B%INAct+it)**2*A%Occ(A%INAct+iq)*B%Occ(B%INAct+it)
!    enddo
! enddo
!! print*, ddii,ddia,ddai,ddaa
! print*, 'direct-direct:',-8*ddii*Tii
! print*, 'direct-direct:',-8*ddii*(Tii+Tia+Tai+Taa)
! tst(1) = ddii*(Tii+Tia+Tai+Taa) + ddia*(Tai+Tii) + ddai*(Tia+Tii) + ddaa*Tii
! print*, -8*tst(1)
! print*, 'aa-ai',-8*ddaa*Tai
! print*, 'aa-aa',-8*ddaa*Taa
! print*, 'ai-ai',-8*ddai*Tai
! print*, 'ai-aa',-8*ddai*Taa
! 
! print*, -8*(ddii+ddaa+ddia+ddai)*(Tia+Tai+Tii+Taa)
! 
!
!! indirect-direct
! allocate(ASB(NBas,NBas))
! tmp1=0
! ASB=0
! call dgemm('N','N',NBas,NBas,NBas,1d0,PAi,NBas,S,NBas,0d0,tmp1,NBas)
! call dgemm('N','N',NBas,NBas,NBas,1d0,tmp1,NBas,PBi,NBas,0d0,ASB,NBas)
! tmp1=0
! tmp2=0
! ! tmp2=mmprod(S,PA)
! call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PAi,NBas,0d0,tmp2,NBas)
! call dgemm('N','N',NBas,NBas,NBas,1d0,ASB,NBas,tmp2,NBas,0d0,tmp1,NBas)
! do j=1,NBas
!    do i=1,NBas
!       tst(2) = tst(2) + tmp1(i,j)*JBi(i,j)
!    enddo
! enddo
! tst(2) = 4*tst(2)
! print*, 'indirect-direct:',tst(2)
!! direct-indirect
! tmp1=0
! tmp2=0
! call dgemm('N','N',NBas,NBas,NBas,1d0,PBi,NBas,S,NBas,0d0,tmp1,NBas)
! call dgemm('N','N',NBas,NBas,NBas,1d0,tmp1,NBas,ASB,NBas,0d0,tmp2,NBas)
! do j=1,NBas
!    do i=1,NBas
!       tst(3) = tst(3) + tmp2(i,j)*JAi(i,j)
!    enddo
! enddo
! tst(3) = 4*tst(3)
! print*, 'direct-indirect:',tst(3)
!
!!indirect-indirect
! Kasb=0
! call make_K(NBas,ASB,Kasb)
! do j=1,NBas
!    do i=1,NBas
!       tst(4) = tst(4) + ASB(i,j)*Kasb(i,j)
!    enddo
! enddo
! tst(4) = -2*tst(4)
! print*, 'indirect-indirect:',tst(4)
!
! print*, 'inac-inact:', sum(tst)
!
! deallocate(ASB)
!
! end block 

 !close(iunit)
 deallocate(tmpAB,tmpB,tmpA)

 exchs2=t1f+sum(t2a)+sum(t2b)+sum(t2c)+t2d
 write(LOUT,*) 'ExchS2: ', exchs2*1000

 deallocate(Vbaa,Vabb,PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
 deallocate(Kb,Qba,Qab,USb,USa) 
 deallocate(tmp2,tmp1,work)
 !deallocate(JJb) 
 deallocate(RDM2Bval,RDM2Aval)

end subroutine e1exchs2

subroutine e2ind(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: e2ba,e2ab,e2iu,e2ic 
double precision :: tmp

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! read EigValA_B
 allocate(EVecA(A%NDimX,A%NDimX),OmA(A%NDimX), &
          EVecB(B%NDimX,B%NDimX),OmB(B%NDimX))
 allocate(AlphaA(A%NDimX,A%NDimX),AlphaB(B%NDimX,B%NDimX), &
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')
 
 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas) 
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas) 

 call calc_resp(EVecA,OmA,AlphaA,0d0,A)
 call calc_resp(EVecB,OmB,AlphaB,0d0,B)

 !tmp=0
 !do i=1,A%NDimX
 !do j=1,A%NDimX
 !   tmp = tmp + AlphaA(i,j)**2
 !enddo
 !enddo
 !print*, 'RMA: ',tmp
 !tmp=0
 !do i=1,B%NDimX
 !do j=1,B%NDimX
 !   tmp = tmp + AlphaB(i,j)**2
 !enddo
 !enddo
 !print*, 'RMB: ',tmp

 e2ba=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    do rs=1,A%NDimX
       ir = A%IndN(1,rs)
       is = A%IndN(2,rs)

       e2ba = e2ba + & 
            WbAA(ip,iq)*AlphaA(pq,rs)*WbAA(ir,is)

    enddo
 enddo
 e2ba = -0.5d0*e2ba
 write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000d0 

 e2ab=0
 do pq=1,B%NDimX
    ip = B%IndN(1,pq)
    iq = B%IndN(2,pq)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       e2ab = e2ab + & 
            WaBB(ip,iq)*AlphaB(pq,rs)*WaBB(ir,is)

    enddo
 enddo
 e2ab = -0.5d0*e2ab
 write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0 

 e2ic = (e2ab + e2ba)
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 
 SAPT%e2ind = e2ic

 deallocate(WaBB,WbAA,AlphaB,AlphaA)
 deallocate(OmB,EVecB,Oma,EVecA)

end subroutine e2ind

subroutine e2ind_apsg(Flags,A,B,SAPT)

implicit none
type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas,ADimEx,BDimEx
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:,:),EVecB(:,:)
double precision,allocatable :: WaBB(:,:),WbAA(:,:)
double precision,allocatable :: AlphaA(:,:),AlphaB(:,:)
integer :: i,j,pq,ip,iq,rs,ir,is
double precision :: termsBA(3), termsAB(3)
double precision :: e2ba,e2ab
double precision :: e2iu,e2ic 
double precision :: e2tmp, tmp

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif
 ADimEx = 2*(A%NDimX + A%NDimN)
 BDimEx = 2*(B%NDimX + B%NDimN)

! read EigValA_B
 allocate(EVecA(ADimEx,BDimEx),OmA(ADimEx), &
          EVecB(BDimEx,BDimEx),OmB(BDimEx))
 allocate(AlphaA(ADimEx,ADimEx),AlphaB(BDimEx,BDimEx), &
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,ADimEx,'PROP_A')
 call readresp(EVecB,OmB,BDimEx,'PROP_B')
 
 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas) 
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas) 

 call calc_resp_apsg(EVecA,OmA,AlphaA,0d0,A)
 call calc_resp_apsg(EVecB,OmB,AlphaB,0d0,B)

 termsBA=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    do rs=1,A%NDimX
       ir = A%IndN(1,rs)
       is = A%IndN(2,rs)

       termsBA(1) = termsBA(1) + & 
            WbAA(ip,iq)*AlphaA(pq,rs)*WbAA(ir,is)

    enddo
 enddo
 termsBA(1) = -0.5d0*termsBA(1)
! write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A) p>q r>s   = ', termsBA(1)*1000 

 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    do ir=1,A%NDimN

      termsBA(2) = termsBA(2) + & 
            WbAA(ip,iq)*AlphaA(pq,A%NDimX+ir)*WbAA(ir,ir)

    enddo
 enddo
 termsBA(2) = -1d0*termsBA(2)

 do ip=1,A%NDimN
    do ir=1,A%NDimN

       termsBA(3) = termsBA(3) + & 
            WbAA(ip,ip)*AlphaA(A%NDimX+ip,A%NDimX+ir)*WbAA(ir,ir)

    enddo
 enddo
 termsBA(3) = -0.5d0*termsBA(3)

 e2ba = sum(termsBA)
 write(LOUT,'(/,1x,a,f16.8)') 'Ind(B--A)   = ', e2ba*1000 

 termsAB=0
 do pq=1,B%NDimX
    ip = B%IndN(1,pq)
    iq = B%IndN(2,pq)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       termsAB(1) = termsAB(1) + & 
            WaBB(ip,iq)*AlphaB(pq,rs)*WaBB(ir,is)


    enddo
 enddo
 termsAB(1) = -0.5d0*termsAB(1)
! write(LOUT,'(1x,a,f16.8)') 'Ind(A--B) p>q r>s  = ', termsAB(1)*1000d0 

 do pq=1,B%NDimX
    ip = B%IndN(1,pq)
    iq = B%IndN(2,pq)
    do ir=1,B%NDimN

       termsAB(2) = termsAB(2) + & 
            WaBB(ip,iq)*AlphaB(pq,B%NDimX+ir)*WaBB(ir,ir)

    enddo
 enddo
 termsAB(2) = -1.0d0*termsAB(2)

 do ip=1,B%NDimN
    do ir=1,B%NDimN

       termsAB(3) = termsAB(3) + & 
            WaBB(ip,ip)*AlphaB(B%NDimX+ip,B%NDimX+ir)*WaBB(ir,ir)

    enddo
 enddo
 termsAB(3) = -0.5d0*termsAB(3)

 e2ab = sum(termsAB)
 write(LOUT,'(1x,a,f16.8)') 'Ind(A--B)   = ', e2ab*1000d0 

 e2tmp = (termsBA(1)+termsAB(1))
 write(LOUT,'(1x,a,f16.8)') 'p>q r>s     = ', e2tmp*1000d0 

 e2ic = (e2ab + e2ba)
 write(LOUT,'(1x,a,f16.8)') 'E2ind       = ', e2ic*1000d0 

end subroutine e2ind_apsg

subroutine e2disp_unc(Flags,A,B,SAPT)
! calculate uncoupled and semi-coupled e2disp
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: kc, ik, ic
double precision,allocatable :: OmA0(:),OmB0(:),&
                                OmA1(:),OmB1(:)
double precision,allocatable :: EVecA0(:),EVecB0(:),& 
                                EVecA1(:),EVecB1(:)
double precision,allocatable :: tmp01(:,:),y0y0(:,:),&  
                                y1y0h(:,:),y1y0(:,:),&
                                y0y1(:,:)
double precision,allocatable :: work(:)
double precision :: fact,dea,deb
double precision :: e2du,e2d
double precision :: e2ds,e2sp,e2ds1,e2ds2
double precision :: e2dw12,e2dw12_sp,e2dsApp
double precision :: inv_omega, inv_om12,tmp
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 
double precision :: Alpha, Beta

! Parameter(SmallE=1.D-3,BigE=1.D8)

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

! read EigValA_B
 allocate(EVecA0(A%NDimX*A%NDimX),OmA0(A%NDimX),&
          EVecB0(B%NDimX*B%NDimX),OmB0(B%NDimX))  
 ! semi-coupled
 if(Flags%IFlag0==0) then
    allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX),&
             EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))  
 endif           

 call readresp(EVecA0,OmA0,A%NDimX,'PROP_A0')
 call readresp(EVecB0,OmB0,B%NDimX,'PROP_B0')
 ! semi-coupled
 if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
    ! CAS
    call readresp(EVecA1,OmA1,A%NDimX,'PROP_A1')
    call readresp(EVecB1,OmB1,B%NDimX,'PROP_B1')
 elseif(Flags%IFlag0==0.and.Flags%ICASSCF==0) then
    ! GVB
    call readEval(OmA1,A%NDimX,'PROP_A1')
    call readEval(OmB1,B%NDimX,'PROP_B1')
 endif

 Alpha = 1.000d0
 Beta  = 1.000d0 

! tran4_gen
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

if(Flags%ISHF==1.and.SAPT%HFCheck) then

write(LOUT,'(/,1x,a)') 'HARTREE-FOCK E2Disp REQUESTED'

allocate(tmp01(A%NDimX,B%NDimX))

 tmp01=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = &
              work(is+(ir-B%num0-1)*dimOB)

       do i=1,A%NDimX

          tmp01(i,rs) = tmp01(i,rs) + &
                       fact * &
                       EVecA1(pq+(i-1)*A%NDimX)
       enddo

    enddo
 enddo

   e2du=0d0
   e2ds2=0d0
   e2ds1=0d0
   tmp=0d0
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      dea = A%OrbE(ip)-A%OrbE(iq)
   !   print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
         deb = B%OrbE(ir)-B%OrbE(is) 
   !      print*, is,ir,is+(ir-B%num0-1)*dimOB,nOVB

       if(OmA0(pq).gt.SmallE.and.OmB0(rs).gt.SmallE&
        .and.OmA0(pq).lt.BigE.and.OmB0(rs).lt.BigE) then

         tmp=0
         do kc=1,B%NDimX
            ik = B%IndN(1,kc)
            ic = B%IndN(2,kc)
            tmp = tmp + work(ic+(ik-B%num0-1)*dimOB)*EVecB1(kc+(rs-1)*B%NDimX)
         enddo

         e2du  = e2du  + work(is+(ir-B%num0-1)*dimOB)**2/(dea+deb)
         e2ds2 = e2ds2 + work(is+(ir-B%num0-1)*dimOB)**2*(Alpha*OmA1(pq)+Beta*OmB1(rs))/(dea+deb)**2
         e2ds1 = e2ds1 + (Beta*tmp+Alpha*tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(dea+deb) 

       endif
 
      enddo
   enddo
    e2sp = 4d0*(e2ds2-e2du)*1000 
    e2ds = (-4d0*e2du-16/sqrt(2d0)*e2ds1+4d0*e2ds2)*1000

    write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0
    write(LOUT,'(1x,a,f16.8)')'E2disp(sp)  = ', e2sp
    write(LOUT,'(1x,a,f16.8)')'E2disp(sc)  = ', e2ds

deallocate(tmp01)
endif

allocate(tmp01(A%NDimX,B%NDimX),y0y0(A%NDimX,B%NDimX))


if(Flags%IFlag0==0) then
 allocate(y1y0h(A%NDimX,B%NDimX),y1y0(A%NDimX,B%NDimX),&
          y0y1(A%NDimX,B%NDimX))
endif

! uncoupled E20disp

tmp01=0
if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
   y1y0h = 0
   y1y0  = 0
   ! unc & semi-cpld
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
  
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
    
         do i=1,A%NDimX
  
            tmp01(i,rs) = tmp01(i,rs) + & 
                         fact * &
                         EVecA0(pq+(i-1)*A%NDimX)
  
            y1y0h(i,rs) = y1y0h(i,rs) + & 
                         fact * &
                         Alpha * &
                         EVecA1(pq+(i-1)*A%NDimX)
                            
         enddo
  
      enddo
   enddo

else 
! unc
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
  
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
    
         do i=1,A%NDimX
  
            tmp01(i,rs) = tmp01(i,rs) + & 
                         fact * &
                         EVecA0(pq+(i-1)*A%NDimX)
                            
         enddo
  
      enddo
   enddo

endif

 
 y0y0=0
 if(Flags%IFlag0==1.or.(Flags%IFlag0==0.and.Flags%ICASSCF==0)) then
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       y0y0(i,j) = y0y0(i,j) + &
                    EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

       enddo
    enddo  
 enddo

 elseif(Flags%IFlag0==0) then

 y0y1 = 0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       y0y0(i,j) = y0y0(i,j) + &
                    EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

       y1y0(i,j) = y1y0(i,j) + &
                   EVecB0(rs+(j-1)*B%NDimX)*y1y0h(i,rs)

       y0y1(i,j) = y0y1(i,j) + &
                   Beta * &
                   EVecB1(rs+(j-1)*B%NDimX)*tmp01(i,rs)

       enddo
    enddo  
 enddo

endif

 e2du=0d0
 e2sp=0d0
 e2ds1=0d0
 e2ds2=0d0
 e2dw12=0d0
 e2dw12_sp=0d0
 e2dsApp=0
 
 
 if(Flags%IFlag0==0.and.Flags%ICASSCF==1) then
 ! CAS
    do j=1,B%NDimX
       do i=1,A%NDimX
          ! remove this if later!!!!
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then
        !  if((OmA0(i)+OmA1(i)).gt.SmallE.and.(OmB0(j)+OmB1(j)).gt.SmallE &
        !   .and.(OmA1(i)+OmA0(i)).lt.BigE.and.(OmB1(j)+OmB0(j)).lt.BigE) then
 
          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          inv_om12 = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du = e2du + y0y0(i,j)**2*inv_omega
          e2ds2 = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(y0y0(i,j)*inv_omega)**2
          !e2ds1 = e2ds1 + tmp02(i,j)*(Alpha*sc10b(i,j)+Beta*sc01b(i,j))*inv_omega
          e2ds1 = e2ds1 + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_omega

          e2dw12 = e2dw12 + y0y0(i,j)**2*inv_om12
          e2dw12_sp = e2dw12_sp + y0y0(i,j)*(y1y0(i,j)+y0y1(i,j))*inv_om12
  
          !print*, i,j,inv_om12,OmA0(i),OmA1(i),OmB0(j),OmB1(j)
          endif
       enddo
    enddo

 elseif(Flags%IFlag0==0.and.Flags%ICASSCF==0) then
 ! GVB
    do j=1,B%NDimX
       do i=1,A%NDimX
          ! remove this if later!!!!
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then
        !  if((OmA0(i)+OmA1(i)).gt.SmallE.and.(OmB0(j)+OmB1(j)).gt.SmallE &
        !   .and.(OmA1(i)+OmA0(i)).lt.BigE.and.(OmB1(j)+OmB0(j)).lt.BigE) then
 
          inv_omega = 1d0/(OmA0(i)+OmB0(j))
          inv_om12 = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du = e2du + y0y0(i,j)**2*inv_omega
          e2ds2 = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(y0y0(i,j)*inv_omega)**2
          e2dw12 = e2dw12 + y0y0(i,j)**2*inv_om12
  
          endif
       enddo
    enddo

 elseif(Flags%IFlag0==1) then

    do j=1,B%NDimX
       do i=1,A%NDimX
          if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
           .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then

          inv_omega = 1d0/(OmA0(i)+OmB0(j))
        !  inv_om12 = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du = e2du + y0y0(i,j)**2*inv_omega
        !  e2dw12 = e2dw12 + y0y0(i,j)**2*inv_om12
   
          endif
       enddo
    enddo

 endif

 SAPT%e2disp_unc = -16d0*e2du
 SAPT%e2disp_sp = -16*e2du+16*e2ds2
 SAPT%e2disp_sc = -16*e2du+16*e2ds2-32*e2ds1

 e2du = -16d0*e2du*1000d0
 e2sp = e2du + 16*e2ds2*1000 
 e2ds = e2du + (16*e2ds2-32*e2ds1)*1000

 e2dw12 = -16*e2dw12*1000
 ! print*, -32*e2dw12_sp
 e2dsApp = e2dw12 -32*e2dw12_sp*1000


 write(LOUT,'(/1x,a,f16.8)')'E2disp(unc) = ', e2du
 write(LOUT,'(1x,a,f16.8)')'E2disp(sp) =  ', e2sp
 write(LOUT,'(1x,a,f16.8)')'E2disp(sc) =  ', e2ds
 write(LOUT,'(/1x,a,f16.8)')'E2disp(w12) = ', e2dw12
 write(LOUT,'(1x,a,f16.8)')'E2disp(sc2) = ', e2dsApp

 close(iunit)
 deallocate(work)



 if(Flags%IFlag0==0) then

    deallocate(y0y1,y1y0,y1y0h)
    deallocate(OmB1,EVecB1,OmA1,EVecA1)
 endif

 deallocate(y0y0,tmp01)
 deallocate(OmB0,EVecB0,OmA0,EVecA0)

end subroutine e2disp_unc

! HERE!!!!
subroutine e2disp(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: kc, ik, ic
double precision,allocatable :: OmA(:), OmB(:), &
                                OmA0(:),OmB0(:),&
                                OmA1(:),OmB1(:)
double precision,allocatable :: EVecA(:), EVecB(:), &
                                EVecA0(:),EVecB0(:),& 
                                EVecA1(:),EVecB1(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),&
                                tmp01(:,:),tmp02(:,:),&  
                                sc10a(:,:),sc10b(:,:),&
                                sc01b(:,:)
double precision,allocatable :: work(:)
double precision,allocatable :: tmp03(:,:) 
double precision :: e2d,fact,tmp
double precision :: e2du,e2sp,dea,deb
double precision :: e2ds,e2ds1,e2ds2
double precision :: inv_omega
! for Be ERPA:
double precision,parameter :: SmallE = 1.D-1
!double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 
double precision :: Alpha, Beta

! Parameter(SmallE=1.D-3,BigE=1.D8)
 
 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 ! print*, A%num0,A%num1,A%num2
 ! print*, nOVA,dimOA,dimVA 
 ! print*, B%num0,B%num1,B%num2

! read EigValA_B
 allocate(EVecA(A%NDimX*A%NDimX),OmA(A%NDimX),&
          EVecB(B%NDimX*B%NDimX),OmB(B%NDimX),&
          EVecA0(A%NDimX*A%NDimX),OmA0(A%NDimX),&
          EVecB0(B%NDimX*B%NDimX),OmB0(B%NDimX))  
 if(Flags%IFlag0==0) then
    allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX),&
             EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))  
 endif           

 call readresp(EVecA,OmA,A%NDimX,'PROP_A')
 call readresp(EVecB,OmB,B%NDimX,'PROP_B')

 if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then
    call readresp(EVecA0,OmA0,A%NDimX,'PROP_A0')
    call readresp(EVecB0,OmB0,B%NDimX,'PROP_B0')

    if(Flags%IFlag0==0) then
       call readresp(EVecA1,OmA1,A%NDimX,'PROP_A1')
       call readresp(EVecB1,OmB1,B%NDimX,'PROP_B1')
    endif
 endif

 !print*, 'Normy:'
 !print*, norm2(OmA),'OmA'
 !print*, norm2(OmB),'OmB'

 !print*, norm2(OmB0)
 !print*, norm2(OmB1)
! print*, norm2(EVecB)
! print*, norm2(EVecB0)
! print*, norm2(EVecA1),'A1'
! print*, norm2(EVecB0),'B0'
! print*, norm2(EVecB1),'B1'
! print*, norm2(EVecB),'B'
! print*, norm2(EVecA),'A'
!
 !print*, 'EXC ENERGIES'
 !do i=1,100!size(OmB)
 !   if(OmB(i).gt.0d0) write(*,*) OmB(i)
 !enddo

! uncoupled
! works with tran4_full
!allocate(work(NInte1))
!open(newunit=iunit,file='TWOMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NInte1)
!
!
!e2du=0d0
!do pq=1,A%NDimX
!   ip = A%IndN(1,pq)
!   iq = A%IndN(2,pq)
!   dea = A%OrbE(ip)-A%OrbE(iq)
!   !print*, iq,ip,iq+ip*(ip-1)/2,NInte1
!   read(iunit,rec=iq+ip*(ip-1)/2) work(1:NInte1)
!   do rs=1,B%NDimX
!      ir = B%IndN(1,rs)
!      is = B%IndN(2,rs)
!      deb = B%OrbE(ir)-B%OrbE(is) 
!      !print*, is,ir,is+ir*(ir-1)/2,NInte1
!      e2du = e2du + work(is+ir*(ir-1)/2)**2/(dea+deb)
!   enddo
!enddo
! write(LOUT,'()') 
! write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0 

! uncoupled
! tran4_full
!allocate(work(NInte1))
!open(newunit=iunit,file='TWOMOAB',status='OLD',&
!     access='DIRECT',form='UNFORMATTED',recl=8*NInte1)

 Alpha = 1.000d0
 Beta  = 1.000d0 

! tran4_gen
allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

if(Flags%ISHF==1.and.SAPT%HFCheck) then

allocate(tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX),&
         tmp03(A%NDimX,B%NDimX))
 tmp01=0
 tmp03=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    !print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = &
              work(is+(ir-B%num0-1)*dimOB)

       do i=1,A%NDimX

          tmp01(i,rs) = tmp01(i,rs) + &
                       fact * &
                       EVecA1(pq+(i-1)*A%NDimX)

!          tmp03(i,rs) = tmp03(i,rs) + &
!                       !fact * sqrt(2d0) * &
!                       fact * sqrt(2d0) * &
!                       EVecA0(pq+(i-1)*A%NDimX)

       enddo

    enddo
 enddo

! tmp02=0
! do j=1,B%NDimX
!    do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
!    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
! 
!       do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!
!       tmp02(pq,j) = tmp02(pq,j) + &
!                    EVecB1(rs+(j-1)*B%NDimX)*&!tmp01(i,rs)
!                    !work(is+(ir-B%num0-1)*dimOB)
!                    tmp03(pq,rs)
!       enddo
!    enddo  
! enddo

   e2du=0d0
   e2ds2=0d0
   e2ds1=0d0
   tmp=0d0
   do pq=1,A%NDimX
      ip = A%IndN(1,pq)
      iq = A%IndN(2,pq)
      dea = A%OrbE(ip)-A%OrbE(iq)
   !   print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
      read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
      do rs=1,B%NDimX
         ir = B%IndN(1,rs)
         is = B%IndN(2,rs)
         deb = B%OrbE(ir)-B%OrbE(is) 
   !      print*, is,ir,is+(ir-B%num0-1)*dimOB,nOVB

       if(OmA0(pq).gt.SmallE.and.OmB0(rs).gt.SmallE&
        .and.OmA0(pq).lt.BigE.and.OmB0(rs).lt.BigE) then

         tmp=0
         do kc=1,B%NDimX
            ik = B%IndN(1,kc)
            ic = B%IndN(2,kc)
            tmp = tmp + work(ic+(ik-B%num0-1)*dimOB)*EVecB1(kc+(rs-1)*B%NDimX)
         enddo

         e2du = e2du + work(is+(ir-B%num0-1)*dimOB)**2/(dea+deb)
         e2ds2 = e2ds2 + work(is+(ir-B%num0-1)*dimOB)**2*(Alpha*OmA1(pq)+Beta*OmB1(rs))/(dea+deb)**2
        ! e2ds1 = e2ds1 + (tmp02(pq,rs)+tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(OmA0(pq)+OmB0(rs))
         e2ds1 = e2ds1 + (Beta*tmp+Alpha*tmp01(pq,rs))*work(is+(ir-B%num0-1)*dimOB)/(dea+deb) !(OmA0(pq)+OmB0(rs))

       endif
 
      enddo
   enddo
    write(LOUT,'()') 
    e2ds = (-4d0*e2du-16/sqrt(2d0)*e2ds1+4d0*e2ds2)*1000
    print*, 'DS1', -16/sqrt(2d0)*e2ds1*1000
    write(LOUT,'(1x,a,f16.8)')'SPOLE:      = ', -4d0*e2du*1000d0 + 4d0*e2ds2*1000d0
    write(LOUT,'(1x,a,f16.8)')'E2disp(sc)  = ', e2ds
    write(LOUT,'(1x,a,f16.8)')'E2disp(unc) = ', -4d0*e2du*1000d0
    write(LOUT,'()') 
    e2ds  = 0
    e2ds1 = 0
    e2ds2 = 0 

deallocate(tmp01,tmp02,tmp03)
! end Hartree-Fock check
endif

allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
        tmp01(A%NDimX,B%NDimX),tmp02(A%NDimX,B%NDimX))

if(Flags%IFlag0==0) then
 allocate(sc10a(A%NDimX,B%NDimX),sc10b(A%NDimX,B%NDimX),&
          sc01b(A%NDimX,B%NDimX))
 sc10a = 0
 sc10b = 0
 sc01b = 0
endif

!print*, 'A-ndimX',A%NDimX

!print*, 'E2disp,CICoeff=1'
!B%CICoef = 0d0
!do i=1,B%NELE !size(B%CICoef)
!   B%CICoef(i)=1d0
!   write(*,*) B%CICoef(i),i  
!enddo


! coupled - 0

if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then

 tmp1=0
 tmp01=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(is+(ir-B%num0-1)*dimOB)
  
       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + & 
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)

          tmp01(i,rs) = tmp01(i,rs) + & 
                       fact * &
                       EVecA0(pq+(i-1)*A%NDimX)

          sc10a(i,rs) = sc10a(i,rs) + & 
                       fact * &
                       Alpha * &
                       EVecA1(pq+(i-1)*A%NDimX)

       enddo

    enddo
 enddo

 tmp2=0
 tmp02=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j) + &
                    EVecB(rs+(j-1)*B%NDimX)*tmp1(i,rs)

       tmp02(i,j) = tmp02(i,j) + &
                    EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)

       sc10b(i,j) = sc10b(i,j) + &
                    EVecB0(rs+(j-1)*B%NDimX)*sc10a(i,rs)

       sc01b(i,j) = sc01b(i,j) + &
                    Beta * &
                    EVecB1(rs+(j-1)*B%NDimX)*tmp01(i,rs)

       enddo
    enddo  
 enddo

elseif(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
    ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOVA
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * &
               work(is+(ir-B%num0-1)*dimOB)
  
       do i=1,A%NDimX
          tmp1(i,rs) = tmp1(i,rs) + & 
                       fact * &
                       EVecA(pq+(i-1)*A%NDimX)

       enddo

    enddo
 enddo

 tmp2=0
 do j=1,B%NDimX
    do i=1,A%NDimX
       do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       tmp2(i,j) = tmp2(i,j) + &
                    EVecB(rs+(j-1)*B%NDimX)*tmp1(i,rs)

       enddo
    enddo  
 enddo

! end GVB select
endif


if(.not.(Flags%ICASSCF==0.and.Flags%ISERPA==0)) then
! uncoupled and semicoupled
 e2du=0d0
! e2ds1=0d0
! e2ds2=0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
       if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
          .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then
       inv_omega = 1d0/(OmA0(i)+OmB0(j))

       e2du = e2du + tmp02(i,j)**2*inv_omega
       e2ds2 = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(tmp02(i,j)*inv_omega)**2
!!       !e2ds1 = e2ds1 + tmp02(i,j)*(Alpha*sc10b(i,j)+Beta*sc01b(i,j))*inv_omega
       e2ds1 = e2ds1 + tmp02(i,j)*(sc10b(i,j)+sc01b(i,j))*inv_omega
!!!       e2ds1 = e2ds1 + tmp02(i,j)*(sc01b(i,j))*inv_omega
!!
       endif
    enddo
 enddo
 SAPT%e2disp_unc = -16d0*e2du
 SAPT%e2disp_sp = -16*e2du+16*e2ds2
 SAPT%e2disp_sc = -16*e2du+16*e2ds2-32*e2ds1

 e2du = -16d0*e2du*1000d0
 e2sp = e2du + 16*e2ds2*1000 
 e2ds= e2du+(16*e2ds2-32*e2ds1)*1000

endif

 e2d = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
       ! remove this if later!!!!
       if(OmA(i).gt.SmallE.and.OmB(j).gt.SmallE&
          .and.OmA(i).lt.BigE.and.OmB(j).lt.BigE) then

       e2d = e2d + tmp2(i,j)**2/(OmA(i)+OmB(j))

       endif
    enddo
 enddo
 SAPT%e2disp = -16d0*e2d

 e2d  = -16d0*e2d*1000d0
 write(LOUT,'(/1x,a,f16.8)')'E2disp(unc) = ', e2du
 write(LOUT,'(1x,a,f16.8)')'E2disp(sp) =  ', e2sp
 write(LOUT,'(1x,a,f16.8)')'E2disp(sc) =  ', e2ds

 write(LOUT,'(/1x,a,f16.8)') 'E2disp      = ',e2d

!! coupled - TEST FULL LOOP
! e2d = 0d0
! do i=1,A%NDimX
!    do j=1,B%NDimX
!    tmp=0
!    do pq=1,A%NDimX
!       ip = A%IndN(1,pq)
!       iq = A%IndN(2,pq)
!       read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
!   
!       do rs=1,B%NDimX
!          ir = B%IndN(1,rs)
!          is = B%IndN(2,rs)
!
!             tmp = tmp + &
!                   (A%CICoef(ip)+A%CICoef(iq)) * &
!                   (B%CICoef(ir)+B%CICoef(is)) * &
!                   EVecA((i-1)*A%NDimX+pq)* & !*work(is+ir*(ir-1)/2)*&
!                   work(is+(ir-B%num0-1)*dimOB)*&
!                   EVecB((j-1)*B%NDimX+rs)
!
!          enddo
!       enddo
!
!          e2d = e2d  + tmp**2/(OmA(i)+OmB(j))
!    enddo
! enddo
! e2d = -16d0*e2d
! 
! print*, 'e2disp: ',e2d*1000

 close(iunit)
 deallocate(work)

 if(Flags%IFlag0==0) then
    deallocate(sc10b,sc01b)
    deallocate(OmB1,EVecB1,OmA1,EVecA1)
 endif

 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,EVecB0,OmA0,EVecA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp



subroutine e2disp_apsg(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ii,jj
integer :: coef,coef2,ADimEx,BDimEx
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact
double precision :: e2d1,e2d2,e2d3,e2d4,e2d,tmp
double precision :: e2du,dea,deb
!double precision,parameter :: SmallE = 1.d-6
double precision,parameter :: SmallE = 1.d-1


 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

! set dimensions
 ADimEx = A%NDimX + A%NDimN
 BDimEx = B%NDimX + B%NDimN
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 dimOA = A%num0+A%num1
 dimFA = NBas
 dimOB = B%num0+B%num1
 dimFB = NBas 
 nOFA = dimOA*dimFA
 nOFB = dimOB*dimFB
 coef  = 1 
 coef2 = 1
 ! Use these coefficients for PINOVEC
 !coef  = 2
 !coef2 = 4

 print*, A%num0,A%num1,A%num2
 print*, B%num0,B%num1,B%num2

! read EigValA_B
 allocate(EVecA(coef2*ADimEx**2),OmA(coef*ADimEx),&
          EVecB(coef2*BDimEx**2),OmB(coef*BDimEx))

 call readresp(EVecA,OmA,coef*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,coef*BDimEx,'PROP_B')

! print*, norm2(EVecA),norm2(EVecB)
! print*, 'OmA'
! do i=1,size(OmA)
!   ! if(OmA(i).gt.0d0)  write(*,*) i,OmA(i)
!    write(*,*) i,",",OmA(i)
! enddo
!
! print*, 'OmB'
! do i=1,size(OmB)
!    !  if(OmB(i).gt.0d0)  write(*,*) i,OmB(i)
!    write(*,*) i,",",OmB(i)
! enddo
 
 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)

 allocate(tmp1(coef*ADimEx,coef*BDimEx),tmp2(coef*ADimEx,coef*BDimEx))

! PART 1: p>q,r>s
 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB

       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
              (B%CICoef(is)+B%CICoef(ir)) * work(is+(ir-1)*dimOB)
 !      print*, work(is+(ir-1)*dimOB)
       !print*, fact
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
          !if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*EVecA((i-1)*coef*ADimEx+pq)

          endif
       enddo
    enddo
 enddo

 tmp2=0
 do j=1,coef*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
    !if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
          !if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
             do rs=1,B%NDimX
                ir = B%IndN(1,rs)
                is = B%IndN(2,rs)
                
                tmp2(i,j) = tmp2(i,j) + &
                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
             
             enddo
          endif
       enddo
    endif
 enddo

!write(*,*) 'tmp1,tmp2:',norm2(tmp1),norm2(tmp2)

 e2d1 = 0d0
 do i=1,coef*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
    !if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
       do j=1,coef*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
          !if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
             e2d1 = e2d1 + tmp2(i,j)**2/(OmA(i)+OmB(j))
             ! print*, OmA(i),OmB(j)
          endif
       enddo
    endif
 enddo

print*, ''
print*, 'PART1: ',-16d0*e2d1*1000d0

! PART 2: p>q,r=s
 tmp1=0
 do pq=1,A%NDimX
    ip = A%IndN(1,pq)
    iq = A%IndN(2,pq)
!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do ir=1,B%NDimN
!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB

       fact = B%CICoef(ir)*(A%CICoef(iq)+A%CICoef(ip)) &
            * work(ir+(ir-1)*dimOB)

       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,ir) = tmp1(i,ir) + &
                  fact*EVecA((i-1)*coef*ADimEx+pq)

          endif
       enddo
    enddo
 enddo
!
 tmp2=0
 do j=1,coef*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do ir=1,B%NDimN
          
                tmp2(i,j) = tmp2(i,j) + &
                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
                     tmp1(i,ir)
                
             enddo
          endif
       enddo
    endif
 enddo

 e2d2 = 0d0
 do i=1,coef*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,coef*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d2 = e2d2 + tmp2(i,j)**2/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 write(LOUT,*) 'PART2: ',-16d0*e2d2*1000d0

! PART 3: p=q,r>s
 tmp1=0
 do ip=1,A%NDimN
    !   ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,B%NDimX
       ir = B%IndN(1,rs)
       is = B%IndN(2,rs)
       !      !print*, is,ir,is+(ir-1)*dimOB,nOFB
       fact = A%CICoef(ip) * &
            (B%CICoef(is)+B%CICoef(ir)) * &
            work(is+(ir-1)*dimOB)
       
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             
             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)
             
          endif
       enddo
    enddo
 enddo
!
 tmp2=0
 do j=1,coef*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do rs=1,B%NDimX
                ir = B%IndN(1,rs)
                is = B%IndN(2,rs)
                
                tmp2(i,j) = tmp2(i,j) + &
                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
                
             enddo
          endif
       enddo
    endif
 enddo
 
 e2d3 = 0d0
 do i=1,coef*ADimEx
    if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
       do j=1,coef*BDimEx
          if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
             e2d3 = e2d3 + tmp2(i,j)**2/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo
 write(LOUT,*) 'PART3: ', -16d0*e2d3*1000d0

! PART 4: p=q,r=s
 tmp1=0
 do ip=1,A%NDimN
    ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
    do ir=1,B%NDimN
       !print*, is,ir,is+(ir-1)*dimOB,nOFB
       fact = A%CICoef(ip)*B%CICoef(ir)* &
              work(ir+(ir-1)*dimOB)
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then

             tmp1(i,ir) = tmp1(i,ir) + &
                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)

          endif
       enddo
    enddo
 enddo
 
! print*, 'tmp1', norm2(tmp1(1:2*ADimEx,1:B%NDimN))
!
 tmp2=0
 do j=1,coef*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             do ir=1,B%NDimN

                tmp2(i,j) = tmp2(i,j) + &
                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
                     tmp1(i,ir)
                
             enddo
          endif
       enddo
    endif
 enddo

! print*, 'tmp2: ',norm2(tmp2)

 e2d4 = 0d0
 do j=1,coef*BDimEx
    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
       do i=1,coef*ADimEx
          if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
             e2d4 = e2d4 + tmp2(i,j)**2/(OmA(i)+OmB(j))
            ! print*, tmp2(i,j)**2,OmA(i),OmB(j)
          endif
       enddo
    endif
 enddo
 print*, 'PART4: ',-16*e2d4*1000

 ! SAPT%e2disp = -16d0*e2d
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000

 close(iunit)
 deallocate(work)
 deallocate(tmp2,tmp1)
 deallocate(OmB,EVecB,OmA,EVecA)

end subroutine e2disp_apsg

subroutine ModABMin_old(Occ,TwoNO,TwoElErf,XKer,ABMin,IndN,IndX,NDimX,NDimKer,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,intent(in) :: NBasis,NDimX,NDimKer,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
double precision,intent(in) :: Occ(NBasis),XKer(NDimKer)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX**2)

double precision :: CICoef(NBasis)
integer :: i,IRow,ICol,ia,ib,iab,ic,id,icd
double precision :: XKer1234,TwoSR,CA,CB,CC,CD
integer,external :: NAddr3,NAddrrK

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do IRow=1,NDimX
   ia = IndN(1,IRow)
   ib = IndN(2,IRow)
   iab = IndX(IRow)
   CA = CICoef(ia)
   CB = CICoef(ib)

   do ICol=1,NDimX
      ic=IndN(1,ICol)
      id=IndN(2,ICol)
      icd=IndX(ICol)
      CC=CICoef(ic)
      CD=CICoef(id)
      
      XKer1234 = XKer(NAddrrK(ia,ib,ic,id))
      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))
      
      ABMin((ICol-1)*NDimX+IRow)=ABMin((ICol-1)*NDimX+IRow) &
                       +4.0d0*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)
   
   enddo
enddo

end subroutine ModABMin_old

subroutine ModABMin(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 64 
integer,intent(in) :: NBasis,NDimX,NGrid,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

double precision :: CICoef(NBasis)
double precision,allocatable :: work(:),batch(:,:),ABKer(:,:)
integer :: i,j,IRow,ICol,ia,ib,iab,ic,id,icd
integer :: offset,batchlen
double precision :: XKer1234,TwoSR,CA,CB,CC,CD
integer,external :: NAddr3,NAddrrK

allocate(work(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

ABKer = 0

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   work(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ia = IndN(1,IRow)
      ib = IndN(2,IRow)
      iab = IndX(IRow)
   
      do ICol=1,NDimX
         ic=IndN(1,ICol)
         id=IndN(2,ICol)
         icd=IndX(ICol)
         if(icd.gt.iab) cycle
    
         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + work(i)* &
            batch(i,ia)*batch(i,ib)*batch(i,ic)*batch(i,id)
         enddo
         
         ABKer(iab,icd) = ABKer(iab,icd) + XKer1234
      
      enddo
   enddo

enddo

do IRow=1,NDimX
   ia = IndN(1,IRow)
   ib = IndN(2,IRow)
   iab = IndX(IRow)
   CA = CICoef(ia)
   CB = CICoef(ib)

   do ICol=1,NDimX
      ic=IndN(1,ICol)
      id=IndN(2,ICol)
      icd=IndX(ICol)
      CC=CICoef(ic)
      CD=CICoef(id)
      if(icd.gt.iab) cycle
 
      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))
      
      ABMin(iab,icd) = ABMin(iab,icd) &
                       +4.0d0*(CA+CB)*(CD+CC)*(ABKer(iab,icd)+TwoSR)
      ABMin(icd,iab) = ABMin(iab,icd)   

   enddo
enddo

deallocate(ABKer,batch,work)

end subroutine ModABMin

subroutine ModABMin_mithap(Occ,SRKer,Wt,OrbGrid,ABMin,IndN,IndX,NDimX,NGrid,NBasis,&
                           twofile,twoerfile)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,intent(in) :: NBasis,NDimX,NGrid
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
character(*),intent(in) :: twofile,twoerfile
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NBasis,NGrid)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

integer :: iunit1,iunit2
integer :: i,j,k,l,kl,irs,ipq,igrd
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(NGrid))

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

do i=1,NGrid
   WtKer(i) = Wt(i)*SRKer(i)
enddo

print*, 'ModABMin_mithap'

open(newunit=iunit1,file=trim(twofile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
open(newunit=iunit2,file=trim(twoerfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1   
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        read(iunit1,rec=kl) work1(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work1,ints1,NBasis)
        read(iunit2,rec=kl) work2(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work2,ints2,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                TwoSR = ints1(i,j)-ints2(i,j)

                XKer1234 = 0
                do igrd=1,NGrid
                   XKer1234 = XKer1234 + WtKer(igrd)* &
                   OrbGrid(l,igrd)*OrbGrid(k,igrd)*OrbGrid(j,igrd)*OrbGrid(i,igrd)
                enddo
          
                ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                               + 4.0d0*Cpq*Crs*(TwoSR+XKer1234)

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

close(iunit1)
close(iunit2)

deallocate(WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_mithap

subroutine readresp(EVec,EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVec(NDim,NDim), EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVec
 read(iunit) EVal

 close(iunit)

end subroutine readresp

subroutine readEval(EVal,NDim,fname)
implicit none

integer :: NDim
double precision :: EVal(NDim)
character(*) :: fname
integer :: iunit

 open(newunit=iunit,file=fname,form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')

 read(iunit) EVal

 close(iunit)

end subroutine readEval

subroutine calc_resp(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec(Mon%NDimX,Mon%NDimX), &
                               EVal(Mon%NDimX)
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha(Mon%NDimX,Mon%NDimX)
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 

 Alpha = 0
 do t=1,Mon%NDimX

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if
    if(EVal(t).gt.SmallE.and.EVal(t).lt.BigE) then

    frac = 8d0*EVal(t)/(EVal(t)**2d0+freq)

    do pq=1,Mon%NDimX
       ip = Mon%IndN(1,pq)
       iq = Mon%IndN(2,pq)
       do rs=1,pq
          ir = Mon%IndN(1,rs)
          is = Mon%IndN(2,rs)

          Alpha(pq,rs) = Alpha(pq,rs) + &
                       (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                       (Mon%CICoef(is)+Mon%CICoef(ir)) * & 
                       EVec(pq,t)*EVec(rs,t)*frac
          Alpha(rs,pq) = Alpha(pq,rs)                

       enddo
    enddo
    endif
 enddo

end subroutine calc_resp

subroutine calc_resp_apsg(EVec,EVal,Alpha,Freq,Mon)
implicit none

type(SystemBlock) :: Mon
double precision,intent(in) :: EVec(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN)), &
                               EVal(2*(Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
double precision,intent(out) :: Alpha(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN))
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
integer :: NDimEx
double precision,parameter :: SmallE = 1.d-6

 NDimEx = 2*(Mon%NDimX + Mon%NDimN)

 Alpha = 0
 ! p>q, r>s
 do t=1,NDimEx

!   if(imag) then
!      v=4d0*omega(p)/(freq**2+omega(p)**2)
!   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
!   end if

    if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
       frac = 8d0*EVal(t)/(EVal(t)**2+freq)

       do pq=1,Mon%NDimX
          ip = Mon%IndN(1,pq)
          iq = Mon%IndN(2,pq)
          do rs=1,pq
             ir = Mon%IndN(1,rs)
             is = Mon%IndN(2,rs)

             Alpha(pq,rs) = Alpha(pq,rs) + &
                          (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                          (Mon%CICoef(is)+Mon%CICoef(ir)) * & 
                          EVec(pq,t)*EVec(rs,t)*frac
             Alpha(rs,pq) = Alpha(pq,rs)                

          enddo
       enddo
    endif
 enddo

! p>q,r=s
! 
  do t=1,NDimEx
 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
!      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if
 
     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)
 
        do pq=1,Mon%NDimX
           ip = Mon%IndN(1,pq)
           iq = Mon%IndN(2,pq)
           do ir=1,Mon%NDimN
 
!              Alpha(pq,ir) = Alpha(pq,ir) + &
!                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
!                          Mon%CICoef(ir) * & 
!                           EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,pq) = Alpha(pq,ir)                

              Alpha(pq,Mon%NDimX+ir) = Alpha(pq,Mon%NDimX+ir) + &
                           (Mon%CICoef(iq)+Mon%CICoef(ip)) * &
                           Mon%CICoef(ir) * & 
                           EVec(pq,t)*EVec(2*Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,pq) = Alpha(pq,Mon%NDimX+ir)                
 
           enddo
        enddo
 
     endif
  enddo
 

 ! p=q,r=s
  do t=1,NDimEx
 
 !   if(imag) then
 !      v=4d0*omega(p)/(freq**2+omega(p)**2)
 !   else
 !      v=4d0*omega(p)/(-freq**2+omega(p)**2)
 !   end if
 
     if(EVal(t).gt.SmallE.and.EVal(t).lt.1d20) then
        frac = 8d0*EVal(t)/(EVal(t)**2+freq)
 
        do ip=1,Mon%NDimN
           do ir=1,ip !Mon%NDimN
 
!              Alpha(ip,ir) = Alpha(ip,ir) + &
!                           Mon%CICoef(ip) * &
!                           Mon%CICoef(ir) * & 
!                           EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
!              Alpha(ir,ip) = Alpha(ip,ir)

              Alpha(Mon%NDimX+ip,Mon%NDimX+ir) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir) + &
                           Mon%CICoef(ip) * &
                           Mon%CICoef(ir) * & 
                           EVec(2*Mon%NDimX+ip,t)*EVec(2*Mon%NDimX+ir,t)*frac
              Alpha(Mon%NDimX+ir,Mon%NDimX+ip) = Alpha(Mon%NDimX+ip,Mon%NDimX+ir)

 
           enddo
        enddo
     endif
  enddo

end subroutine calc_resp_apsg

end module sapt_ener
