module sapt_ener
use types
use tran
use exmisc
use timing

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

subroutine e1exchs2(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
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
double precision :: Tcpu,Twall
double precision,parameter :: Half=0.5d0
double precision,external  :: trace,FRDM2,FRDM2GVB

! set dimensions
 NBas = A%NBasis 
 !dimOA = A%INAct+A%NAct 
 dimOA = A%num0+A%num1
 GdimOA = A%num0+A%num1
 !dimOB = B%INAct+B%NAct 
 dimOB = B%num0+B%num1
 GdimOB = B%num0+B%num1

 call clock('START',Tcpu,Twall)

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

! allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
!          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 allocate(RDM2Aval(dimOA,dimOA,dimOA,dimOA),&
          RDM2Bval(dimOB,dimOB,dimOB,dimOB))

 if(Flags%ICASSCF==1) then
    ! CAS
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val

    !do l=1,dimOA
    !   do k=1,dimOA 
    !      do j=1,dimOA
    !         do i=1,dimOA
    !            RDM2Aval(i,j,k,l) = FRDM2(i,k,j,l,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo
    !do l=1,dimOB
    !   do k=1,dimOB 
    !      do j=1,dimOB
    !         do i=1,dimOB
    !            RDM2Bval(i,j,k,l) = FRDM2(i,k,j,l,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo

 elseif(Flags%ICASSCF==0) then

    ! GVB
    RDM2Aval = A%RDM2val
    RDM2Bval = B%RDM2val

    !do l=1,dimOA
    !   do k=1,dimOA 
    !      do j=1,dimOA
    !         do i=1,dimOA
    !            RDM2Aval(i,j,k,l) = FRDM2GVB(i,k,j,l,A%Occ,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo
    !do l=1,dimOB
    !   do k=1,dimOB 
    !      do j=1,dimOB
    !         do i=1,dimOB
    !            RDM2Bval(i,j,k,l) = FRDM2GVB(i,k,j,l,B%Occ,NBas)
    !         enddo
    !      enddo
    !   enddo
    !enddo

 endif


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

! call tran3Q_full(NBas,dimOA,A%CMO,Qba,'TWOA3B')
! call tran3Q_full(NBas,dimOB,B%CMO,Qab,'TWOB3A')

 call tran3MO_Q(NBas,dimOA,A%CMO,Qba,'TWOA3B')
 call tran3MO_Q(NBas,dimOB,B%CMO,Qab,'TWOB3A')

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
! write(LOUT,*) 'T1 ',t1f

! T2d
 t2d = -2d0*SAPT%Vnn*t1(2)
! write(LOUT,*) 'T2d',t2d

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
! write(LOUT,*) 
! print*, 'T2c(1)',t2c(1) 
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
       !t2c(2) = t2c(2) + sum(RDM2Bval(:,:,iq,is)*Vabb(:,:)*PAbb(is,iq))
       t2c(2) = t2c(2) + sum(RDM2Bval(1:dimOB,1:dimOB,iq,is)*Vabb(1:dimOB,1:dimoB)*PAbb(is,iq))
    enddo
 enddo    
 t2c(2) = -2d0*t2c(2)
! write(LOUT,*) 'T2c(2) ',t2c(2)
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
 !print*, 'T2b(1)',t2b(1)
! !test
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:NT', -2d0*trace(work,NBas)  
! work=0
! call dgemm('N','T',NBas,NBas,NBas,1d0,tmp1,NBas,tmp2,NBas,0d0,work,NBas)
! print*, 'TEST:TN', -2d0*trace(work,NBas)  

!! Full NBas 
! do ir=1,NBas
!    do ip=1,NBas
!       do is=1,NBas
!          do iq=1,NBas
!             !t2b(2) = t2b(2) + FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas) * &
!             t2b(2) = t2b(2) + FRDM2GVB(ip,iq,ir,is,A%Occ,NBas) * &
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
       !t2b(2) = t2b(2) + sum(RDM2Aval(:,:,iq,is)*Vbaa(:,:)*PBaa(is,iq))
       t2b(2) = t2b(2) + sum(RDM2Aval(1:dimOA,1:dimOA,iq,is)*Vbaa(1:dimOA,1:dimOA)*PBaa(is,iq))
    enddo
 enddo    
 t2b(2) = -2d0*t2b(2)
! write(LOUT,*) 'T2b(2) ',t2b(2)

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
! write(LOUT,*) 'T2a(1)',t2a(1)

 open(newunit=iunit,file='TWOA3B',status='OLD',&
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)
     access='DIRECT',form='UNFORMATTED',recl=8*dimOA**2)

!!! Qba
!!! old: 
! If(Flags%ICASSCF==1) then
!    do ir=1,NBas
!       do ip=1,ir
!          !read(iunit,rec=ip+ir*(ir-1)/2) work
!          read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOA,1:dimOA)
!   
!          if(ip==ir) then
!    
!            do is=1,dimOA !NBas
!                do iq=1,dimOA !NBas
!                     t2a(2) = t2a(2) + work(iq,is)* &
!                              FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)
!                enddo
!             enddo
!   
!          else
!            
!             do is=1,dimOA !NBas
!                do iq=1,dimOA !NBas
!                   t2a(2) = t2a(2) + work(iq,is)* &
!                           (FRDM2(ip,iq,ir,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas)+ &
!                            FRDM2(ir,iq,ip,is,A%RDM2,A%Occ,A%Ind2,A%NAct,NBas))
!   
!                enddo
!             enddo
!     
!          endif
!   
!       enddo
!    enddo
!
! elseif(Flags%ICASSCF==0) then
!
!    do ir=1,NBas
!       do ip=1,ir
!          !read(iunit,rec=ip+ir*(ir-1)/2) work
!          read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOA,1:dimOA)
!   
!          if(ip==ir) then
!    
!            do is=1,dimOA !NBas
!                do iq=1,dimOA !NBas
!                     t2a(2) = t2a(2) + work(iq,is)* &
!                              FRDM2GVB(ip,iq,ir,is,A%Occ,NBas)
!                enddo
!             enddo
!   
!          else
!            
!             do is=1,dimOA !NBas
!                do iq=1,dimOA !NBas
!                   t2a(2) = t2a(2) + work(iq,is)* &
!                           (FRDM2GVB(ip,iq,ir,is,A%Occ,NBas)+ &
!                            FRDM2GVB(ir,iq,ip,is,A%Occ,NBas))
!   
!                enddo
!             enddo
!     
!          endif
!   
!       enddo
!    enddo
!
! endif

!! new - but more tests needed!
! Qba 
 do ir=1,NBas
    do ip=1,ir
      !read(iunit,rec=ip+ir*(ir-1)/2) work
      read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOA,1:dimOA)

       if(ip==ir) then

         if(ip<=dimOA) then 
            do is=1,dimOA
               do iq=1,dimOA
                    t2a(2) = t2a(2) + work(iq,is)* &
                             RDM2Aval(ip,ir,iq,is)
                enddo
             enddo
         else
             do iq=1,dimOA
                  t2a(2) = t2a(2) + work(iq,iq)* &
                           2d0*A%Occ(ip)*A%Occ(iq)
             enddo
         endif 
      
       else

          if(ir<=dimOA) then
             do is=1,dimOA 
                do iq=1,dimOA 
                     t2a(2) = t2a(2) + work(iq,is)* &
                            (RDM2Aval(ip,ir,iq,is)+RDM2Aval(ir,ip,iq,is))
                enddo
             enddo
          endif

       endif

    enddo
 enddo

 t2a(2) = -2*t2a(2)
! print*, 'T2a(2) ',t2a(2)
 close(iunit)

! Qab
 open(newunit=iunit,file='TWOB3A',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*dimOB**2)
     !access='DIRECT',form='UNFORMATTED',recl=8*NBas**2)

 If(Flags%ICASSCF==1) then

    do ir=1,NBas
       do ip=1,ir
          !read(iunit,rec=ip+ir*(ir-1)/2) work
          read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOB,1:dimOB)

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

 elseif(Flags%ICASSCF==0) then

    do ir=1,NBas
       do ip=1,ir
          read(iunit,rec=ip+ir*(ir-1)/2) work(1:dimOB,1:dimOB)

          if(ip==ir) then
    
            do is=1,dimOB
                do iq=1,dimOB
                     t2a(3) = t2a(3) + work(iq,is)* &
                              FRDM2GVB(ip,iq,ir,is,B%Occ,NBas)
                enddo
             enddo

          else
            
             do is=1,dimOB
                do iq=1,dimOB
                   t2a(3) = t2a(3) + work(iq,is)* &
                           (FRDM2GVB(ip,iq,ir,is,B%Occ,NBas)+ &
                            FRDM2GVB(ir,iq,ip,is,B%Occ,NBas))
                enddo
             enddo
     
          endif

       enddo
    enddo

 endif
!
!! Qab - new 
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

 t2a(3) = -2*t2a(3)
! print*, 'T2a(3) ',t2a(3)
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
    call dgemm('N','T',dimOB**2,dimOA,dimOB,1d0,RDM2Bval(:,:,:,is),dimOB**2,Sab,NBas,0d0,tmpB(:,:,:,is),dimOB**2)
 enddo

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

 work=0
 open(newunit=iunit,file='TMPOOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*GdimOB**2)

 do ir=1,A%INAct+A%NAct
    do ip=1,A%INAct+A%NAct
     !  print*, ip,ir,ip+(ir-1)*GdimOA 
      read(iunit,rec=ip+(ir-1)*GdimOA) work(1:GdimOB,1:GdimOB)

      t2a(4) = t2a(4) + sum(work(1:dimOB,1:dimOB)*tmpAB(ip,ir,1:dimOB,1:dimOB))

    enddo
 enddo
 close(iunit)
 t2a(4) = -2*t2a(4)
! print*, 't2a(4): ',t2a(4)

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
 SAPT%exchs2 = exchs2
 write(LOUT,'(1x,a,f16.8)') 'ExchS2        = ', exchs2*1000d0


 deallocate(Vbaa,Vabb,PBaa,PAbb,Vb,Va,PB,PA,Sab,S)
 deallocate(Kb,Qba,Qab,USb,USa) 
 deallocate(tmp2,tmp1,work)
 !deallocate(JJb) 
 deallocate(RDM2Bval,RDM2Aval)

 call clock('E1exch(S2)',Tcpu,Twall)

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
integer :: coef,coef2
double precision :: e2ba,e2ab
double precision :: e2iu,e2ic 
double precision :: e2tmp, tmp

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

 coef  = 1
 coef2 = 1
 
 ADimEx = A%NDimX + A%NDimN
 BDimEx = B%NDimX + B%NDimN

! read EigValA_B
 allocate(EVecA(coef*ADimEx,coef*ADimEx),OmA(coef*ADimEx),&
          EVecB(coef*BDimEx,coef*BDimEx),OmB(coef*BDimEx))
 allocate(AlphaA(ADimEx,ADimEx),AlphaB(BDimEx,BDimEx), &
          WaBB(NBas,NBas),WbAA(NBas,NBas))

 call readresp(EVecA,OmA,coef*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,coef*BDimEx,'PROP_B')

 call tran2MO(A%WPot,B%CMO,B%CMO,WaBB,NBas) 
 call tran2MO(B%WPot,A%CMO,A%CMO,WbAA,NBas) 

 print*, 'HERE?'

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

 deallocate(EVecA,OmA,EVecB,OmB)
 deallocate(AlphaA,AlphaB,WaBB,WbAA)

end subroutine e2ind_apsg

subroutine e2disp_unc(Flags,A,B,SAPT)
! calculate uncoupled and semi-coupled e2disp
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(EBlockData)               :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable   :: SBlockA(:),SBlockB(:)
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)

integer :: NBas, NInte1,NInte2
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB,nblkA,nblkB
integer :: iunit,ival
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: iblk,ipos,kc,ik,ic
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
 nblkA  = 1+NBas-A%NAct
 nblkB  = 1+NBas-B%NAct

! read EigValA_B
! allocate(EVecA0(A%NDimX*A%NDimX),OmA0(A%NDimX),&
!          EVecB0(B%NDimX*B%NDimX),OmB0(B%NDimX))  
 allocate(OmA0(A%NDimX),OmB0(B%NDimX))
! semi-coupled
 if(Flags%IFlag0==0) then
    allocate(EVecA1(A%NDimX*A%NDimX),OmA1(A%NDimX),&
             EVecB1(B%NDimX*B%NDimX),OmB1(B%NDimX))  
 endif           

! call readresp(EVecA0,OmA0,A%NDimX,'PROP_A0')
! call readresp(EVecB0,OmB0,B%NDimX,'PROP_B0')

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

 allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
 call convert_XY0_to_Y01(A,Y01BlockA,OmA0,NBas,'XY0_A')
 call convert_XY0_to_Y01(B,Y01BlockB,OmB0,NBas,'XY0_B')

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

tmp01 = 0
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
  
            !tmp01(i,rs) = tmp01(i,rs) + & 
            !             fact * &
            !             EVecA0(pq+(i-1)*A%NDimX)
  
            y1y0h(i,rs) = y1y0h(i,rs) + & 
                         fact * &
                         Alpha * &
                         EVecA1(pq+(i-1)*A%NDimX)
                            
         enddo
 
         associate(Y => Y01BlockA(pq))
            tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
         end associate

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
    
         !do i=1,A%NDimX
  
         !   tmp01(i,rs) = tmp01(i,rs) + & 
         !                fact * &
         !                EVecA0(pq+(i-1)*A%NDimX)
         !                   
         !enddo
 
         associate(Y => Y01BlockA(pq))
            tmp01(Y%l1:Y%l2,rs) = tmp01(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
         end associate

      enddo
   enddo

endif

y0y0 = 0
if(Flags%IFlag0==1.or.(Flags%IFlag0==0.and.Flags%ICASSCF==0)) then

   !do j=1,B%NDimX
   !   do i=1,A%NDimX
   !      do rs=1,B%NDimX
   !      ir = B%IndN(1,rs)
   !      is = B%IndN(2,rs)
   
   !      y0y0(i,j) = y0y0(i,j) + &
   !                   EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)
   
   !      enddo
   !   enddo  
   !enddo
   
   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs))
        call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,y0y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo

elseif(Flags%IFlag0==0) then

y0y1 = 0

   ! unc
   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs)) 
        call dger(A%NDimX,Y%n,1d0,tmp01(:,rs),1,Y%vec0,1,y0y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo

   y1y0 = 0
   do rs=1,B%NDimX
      associate(Y => Y01BlockB(rs)) 
        call dger(A%NDimX,Y%n,1d0,y1y0h(:,rs),1,Y%vec0,1,y1y0(:,Y%l1:Y%l2),A%NDimX)
      end associate
   enddo
   call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp01,A%NDimX,EVecB1,B%NDimX,0d0,y0y1,A%NDimX)

!do j=1,B%NDimX
!   do i=1,A%NDimX
!      do rs=1,B%NDimX
!      ir = B%IndN(1,rs)
!      is = B%IndN(2,rs)
!
!      y0y0(i,j) = y0y0(i,j) + &
!                   EVecB0(rs+(j-1)*B%NDimX)*tmp01(i,rs)
!
!      y1y0(i,j) = y1y0(i,j) + &
!                  EVecB0(rs+(j-1)*B%NDimX)*y1y0h(i,rs)
!
!      y0y1(i,j) = y0y1(i,j) + &
!                  Beta * &
!                  EVecB1(rs+(j-1)*B%NDimX)*tmp01(i,rs)
!
!      enddo
!   enddo  
!enddo

endif

 e2du      = 0d0
 e2sp      = 0d0
 e2ds1     = 0d0
 e2ds2     = 0d0
 e2dw12    = 0d0
 e2dw12_sp = 0d0
 e2dsApp   = 0d0
 
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
          inv_om12  = 1d0/(OmA0(i)+OmA1(i)+OmB0(j)+OmB1(j))
   
          e2du   = e2du + y0y0(i,j)**2*inv_omega
          e2ds2  = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(y0y0(i,j)*inv_omega)**2
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
! e2dsApp = e2dw12 -32*e2dw12_sp*1000

 write(LOUT,'(/1x,a,f16.8)')'E2disp(unc) = ', e2du
 write(LOUT,'(1x,a,f16.8)')'E2disp(sp) =  ' , e2sp
 write(LOUT,'(1x,a,f16.8)')'E2disp(sc) =  ' , e2ds
 write(LOUT,'(/1x,a,f16.8)')'E2disp(w12) = ', e2dw12
! write(LOUT,'(1x,a,f16.8)')'E2disp(sc2) = ' , e2dsApp

 close(iunit)
 deallocate(work)

 ! deallocate Y01Block
 do i=1,A%NDimX
    associate(Y => Y01BlockA(i))
      deallocate(Y%vec0)
    end associate
 enddo
 do i=1,B%NDimX
    associate(Y => Y01BlockB(i))
      deallocate(Y%vec0)
    end associate
 enddo
 deallocate(Y01BlockB,Y01BlockA)

 if(Flags%IFlag0==0) then

    deallocate(y0y1,y1y0,y1y0h)
    deallocate(OmB1,EVecB1,OmA1,EVecA1)
 endif

 deallocate(y0y0,tmp01)
 deallocate(OmB0,OmA0)
 !deallocate(OmB0,EVecB0,OmA0,EVecA0)

end subroutine e2disp_unc

subroutine e2dispCAS(e2d,Flags,A,B,SAPT,ACAlpha,NBasis)
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT
type(EBlockData)               :: SBlockAIV,SBlockBIV
type(EBlockData),allocatable   :: SBlockA(:),SBlockB(:)
type(Y01BlockData),allocatable :: Y01BlockA(:),Y01BlockB(:)
integer,intent(in)           :: NBasis
double precision,intent(in)  :: ACAlpha
double precision,intent(out) :: e2d

integer          :: iunit
integer          :: i,j,pq,rs
integer          :: ip,iq,ir,is
integer          :: dimOA,dimVA, &
                    dimOB,dimVB,nOVA,nOVB
double precision :: fact
integer,allocatable          :: IGemA(:),IGemB(:)
double precision,allocatable :: OmA(:),OmB(:), &
                                EVecA(:),EVecB(:)
double precision,allocatable :: work(:)
double precision,allocatable :: tmp(:,:),tmp1(:,:)
! for Be ERPA:
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: BigE   = 1.D8 
double precision,parameter :: SmallE = 1.D-3

! print thresholds
if(SAPT%IPrint>1) then 
   write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
   write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
   write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
endif

! set dimensions
dimOA = A%num0+A%num1
dimVA = A%num1+A%num2
dimOB = B%num0+B%num1
dimVB = B%num1+B%num2
nOVA  = dimOA*dimVA
nOVB  = dimOB*dimVB

! fix IGem for A
allocate(IGemA(NBasis),IGemB(NBasis))
do i=1,A%INAct
   IGemA(i) = 1
enddo
do i=A%INAct+1,dimOA
   IGemA(i) = 2
enddo
do i=dimOA+1,NBasis
   IGemA(i) = 3
enddo
! fix IGem for B
do i=1,B%INAct
   IGemB(i) = 1
enddo
do i=B%INAct+1,dimOB
   IGemB(i) = 2
enddo
do i=dimOB+1,NBasis
   IGemB(i) = 3
enddo

!do i=1,NBasis
!   print*, i, IGemA(i),IGemA(i)
!enddo

!allocate(EVecA(A%NDimX*A%NDimX),EVecB(B%NDimX*B%NDimX),&
!         OmA(A%NDimX),OmB(B%NDimX))
allocate(OmA(A%NDimX),OmB(B%NDimX))
allocate(tmp(A%NDimX,B%NDimX),tmp1(A%NDimX,B%NDimX))

!call readresp(EVecA,OmA,A%NDimX,'PROP_A')
!call readresp(EVecB,OmB,B%NDimX,'PROP_B')

!new
allocate(Y01BlockA(A%NDimX),Y01BlockB(B%NDimX))
call convert_XY0_to_Y01(A,Y01BlockA,OmA,NBasis,'XY0_A')
call convert_XY0_to_Y01(B,Y01BlockB,OmB,NBasis,'XY0_B')

! check for negative eigenvalues
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative B!',i,OmB(i)
enddo

allocate(work(nOVB))
open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

tmp1 = 0
do pq=1,A%NDimX
   ip = A%IndN(1,pq)
   iq = A%IndN(2,pq)
   read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)
   do rs=1,B%NDimX
      ir = B%IndN(1,rs)
      is = B%IndN(2,rs)

! the version below was our trial of computing dispersion which must be ADDED TO supermoleculr CAS
!      fact1 = (A%CICoef(iq)+A%CICoef(ip)) * &
!              (B%CICoef(is)+B%CICoef(ir)) * &
!              work(is+(ir-B%num0-1)*dimOB)
!
!      fact2 = ACAlpha*fact1
!
      ! remove active part 
!      if(IGemA(ip)==2.and.IGemA(iq)==2.and. &
!         IGemB(ir)==2.and.IGemB(is)==2) then
!         fact2 = fact1
!         fact1 = 0d0
!      endif
!
! dispersion included in supermolecular CAS is computed
      if(IGemA(ip)==2.and.IGemA(iq)==2.and. &
         IGemB(ir)==2.and.IGemB(is)==2) then
         fact = (A%CICoef(iq)+A%CICoef(ip)) * &
                (B%CICoef(is)+B%CICoef(ir)) * &
                work(is+(ir-B%num0-1)*dimOB)
      else
         fact = 0.d0
      endif
 
      !do i=1,A%NDimX
      !   tmp1(i,rs) = tmp1(i,rs) + &
      !                fact1 * &
      !                EVecA(pq+(i-1)*A%NDimX)

      !   tmp2(i,rs) = tmp2(i,rs) + &
      !                fact2 * &
      !                EVecA(pq+(i-1)*A%NDimX)
      !enddo

       associate(Y => Y01BlockA(pq))
         tmp1(Y%l1:Y%l2,rs) = tmp1(Y%l1:Y%l2,rs) + fact * Y%vec0(1:Y%n)
       end associate

   enddo
enddo

close(iunit)
deallocate(work)

!print*, 'tmp1',norm2(tmp1)

! second multiplication
tmp = 0
do rs=1,B%NDimX
   associate(Y => Y01BlockB(rs)) 
     call dger(A%NDimX,Y%n,1d0,tmp1(:,rs),1,Y%vec0,1,tmp(:,Y%l1:Y%l2),A%NDimX)
   end associate
enddo
tmp1 = tmp
! old
!call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp1,A%NDimX,EVecB,B%NDimX,0d0,tmp,A%NDimX)
!tmp1 = tmp

e2d = 0d0
do j=1,B%NDimX
   do i=1,A%NDimX
      if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
         .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

         e2d = e2d + tmp1(i,j)**2/(OmA(i)+OmB(j))

      endif
   enddo
enddo
SAPT%e2disp = -32d0*e2d
e2d  = -32d0*e2d*1000d0

write(LOUT,'(/1x,a,f16.8)') 'E2disp(CAS) = ',e2d

! deallocate Y01Block
do i=1,A%NDimX
   associate(Y => Y01BlockA(i))
     deallocate(Y%vec0)
   end associate
enddo
do i=1,B%NDimX
   associate(Y => Y01BlockB(i))
     deallocate(Y%vec0)
   end associate
enddo
deallocate(Y01BlockB,Y01BlockA)

! old 
!deallocate(A%EigY,A%EigX,A%Eig)
!deallocate(B%EigY,B%EigX,B%Eig)

deallocate(IGemB,IGemA)
deallocate(tmp1,tmp)
deallocate(OmB,OmA)
!deallocate(OmB,OmA,EVecB,EVecA)
 
end subroutine e2dispCAS

! HERE!!!
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
!double precision,parameter :: SmallE = 1.D-1
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

! print thresholds
 if(SAPT%IPrint>1) then 
    write(LOUT,'(/,1x,a)') 'Thresholds in E2disp:'
    write(LOUT,'(1x,a,2x,e15.4)') 'SmallE      =', SmallE
    write(LOUT,'(1x,a,2x,e15.4)') 'BigE        =', BigE
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

  !print*, 'e2disp'
  !print*, A%num0,A%num1,A%num2
  !print*, nOVA,dimOA,dimVA 
  !print*, B%num0,B%num1,B%num2

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


! coupled 
do i=1,A%NDimX
   if(OmA(i)<0d0) write(LOUT,*) 'Negative A!',i,OmA(i)
enddo
do i=1,B%NDimX
   if(OmB(i)<0d0) write(LOUT,*) 'Negative B!',i,OmB(i)
enddo


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
      ! if(OmA0(i).gt.SmallE.and.OmB0(j).gt.SmallE&
      !    .and.OmA0(i).lt.BigE.and.OmB0(j).lt.BigE) then

       if(abs(OmA0(i)).gt.SmallE.and.abs(OmB0(j)).gt.SmallE&
          .and.abs(OmA0(i)).lt.BigE.and.abs(OmB0(j)).lt.BigE) then


       inv_omega = 1d0/(OmA0(i)+OmB0(j))

       e2du = e2du + tmp02(i,j)**2*inv_omega
       e2ds2 = e2ds2 + (Alpha*OmA1(i)+Beta*OmB1(j))*(tmp02(i,j)*inv_omega)**2
       !e2ds1 = e2ds1 + tmp02(i,j)*(Alpha*sc10b(i,j)+Beta*sc01b(i,j))*inv_omega
       e2ds1 = e2ds1 + tmp02(i,j)*(sc10b(i,j)+sc01b(i,j))*inv_omega
       !e2ds1 = e2ds1 + tmp02(i,j)*(sc01b(i,j))*inv_omega

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
       if(abs(OmA(i)).gt.SmallE.and.abs(OmB(j)).gt.SmallE&
          .and.abs(OmA(i)).lt.BigE.and.abs(OmB(j)).lt.BigE) then

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

 ! write amplitude to a file
 call writeampl(tmp2,'PROP_AB')

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

 if(allocated(EVecA0)) print*, 'EVecA0'
 if(allocated(EVecB0)) print*, 'EVecB0'

 deallocate(tmp02,tmp01,tmp2,tmp1)
 deallocate(OmB0,EVecB0,OmA0,EVecA0,OmB,EVecB,OmA,EVecA)

end subroutine e2disp

subroutine term_Y_SR(tpqrs,Sab,elst,vnn,ADimX,BDimX,NBas,A,B)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: ADimX,BDimX,NBas
double precision,intent(in) :: elst,vnn
double precision,intent(in) :: Sab(NBas,NBas),tpqrs(2*ADimX,2*BDimX)

integer :: i,j
integer :: ip,iq,ir,is,ipq,irs
double precision :: termY,tmp(NBas,NBas)

tmp = transpose(Sab)

termY = 0
do j=1,BDimX

   ir = B%IndN(1,j) 
   is = B%IndN(2,j)
   irs = j

   do i=1,ADimX

      ip = A%IndN(1,i)
      iq = A%IndN(2,i)
      ipq = i
  
      termY = termY + tpqrs(ipq,irs)*Sab(ip,is)*tmp(ir,iq) &
                    + tpqrs(ipq,BDimX+irs)*Sab(ip,ir)*tmp(is,iq) &
                    + tpqrs(ADimX+ipq,irs)*Sab(iq,is)*tmp(ir,ip) &
                    + tpqrs(ADimX+ipq,BDimX+irs)*Sab(iq,ir)*tmp(is,ip)

 
   enddo
enddo
!print*, 'termY-000',termY/16d0
termY = 0.5d0*(elst-vnn)*termY

print*, 'termY-t',termY

end subroutine term_Y_SR

subroutine term_X_A1_SR(tpqrs,Sab,Vaab,Vbab,posA,posB,ADimX,BDimX,NBas,A,B)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vaab(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*ADimX,2*BDimX)

integer :: i,j,k,l,kl
integer :: ip,iq,ir,is,ipq,irs
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: termA1,nelA,nelB,tmp(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 allocate(work(NBas*NBas),ints(NBas,NBas))

 !(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 termA1 = 0
 work = 0
 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       ir = k
       iq = l

       call ints_modify(NBas,dimOB,ints,NBas,work,Sab(iq,ir)/nelA,Vaab,Vbab(iq,ir)/nelB,Sab)

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                termA1 = termA1 + tpqrs(ipq,irs)*ints(ip,is)
 
             endif

          enddo
       enddo

       ! YY
       ir = k
       iq = l  
  
       call ints_modify(NBas,dimOB,ints,NBas,work,Vaab(iq,ir)/nelA,Sab,Sab(iq,ir)/nelB,Vbab)

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                termA1 = termA1 + tpqrs(ADimX+ipq,BDimX+irs)*ints(ip,is)

             endif
   
          enddo
       enddo

    enddo
 enddo
 
 close(iunit)

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

       call ints_modify(NBas,NBas,ints,NBas,work,Sab(iq,is)/nelA,Vaab,Vbab(iq,is)/nelB,Sab)

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)
             
             if(ipq/=0.and.irs/=0) then

               termA1 = termA1 + tpqrs(ipq,BDimX+irs)*ints(ip,ir)

             endif

          enddo
       enddo

       ! YX 
       call ints_modify(NBas,NBas,ints,NBas,work,Vaab(iq,is)/nelA,Sab,Sab(iq,is)/nelB,Vbab)

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)
             
             if(ipq/=0.and.irs/=0) then

               termA1 = termA1 + tpqrs(ADimX+ipq,irs)*ints(ip,ir)

             endif

          enddo
       enddo

    enddo
 enddo

 close(iunit)

 termA1 = -0.5d0*termA1 
 print*, 'termA1',termA1

 ! test full
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          NBas,A%CMO,&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFFFABBA','AOTWOSORT')

 !(FF|FF):(AB|BA)
 open(newunit=iunit,file='FFFFABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 termA1 = 0
 work = 0
 ints = 0
 kl = 0
 do l=1,NBas
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       is = k
       ip = l

       call ints_modify(NBas,NBas,ints,NBas,work,Sab(ip,is)/nelA,Vaab,Vbab(ip,is)/nelB,Sab)

       do iq=1,NBas
          do ir=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)
  
             if(ipq/=0.and.irs/=0) then

                termA1 = termA1 + tpqrs(ipq,irs)*ints(iq,ir)

             endif
 
          enddo
       enddo

    enddo
 enddo

 close(iunit,status='delete')

 termA1 = -0.5d0*termA1 
 print*, 'termA1(F)',termA1

 deallocate(ints,work)

end subroutine term_X_A1_SR

subroutine term_X_A3_XX_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbaa,posA,posB,ADimX,BDimX,NBas,A,B,isYY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbaa(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: isYY
character(*),intent(in) :: IntKFile

integer :: i,j,k,l,kl
integer :: iq,ib,ip,ir,ia,ic,iac,ipr
integer :: iunit,offA,offB
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,valTr,val3B,val4B,termXX
double precision :: nelA,nelB
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas)
double precision :: Emat(NBas,NBas),Amat(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

 Amat = 0
 do ib=1,dimOB
    do ir=1,NBas 
       do ip=1,NBas
          Amat(ip,ir) = Amat(ip,ir) + B%Occ(ib)*Sab(ip,ib)*Sab(ir,ib)
       enddo
    enddo
 enddo
 
 valTr = 0
 do iq=1,dimOA
    valTr = valTr + A%Occ(iq)*Amat(iq,iq)
 enddo
 
 Sbb = 0
 do ic=1,NBas
    do ia=1,NBas
       do iq=1,dimOA
          Sbb(ia,ic) = Sbb(ia,ic) + A%Occ(iq)*Sab(iq,ia)*Sab(iq,ic)
       enddo
    enddo
 enddo

 if(isYY) then
    offA = ADimX
    offB = BDimX
 else
    offA = 0
    offB = 0 
 endif

!(FO|FO):(AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

 ! one loop over integrals
 termXX = 0
 Emat = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOA)
 
       call ints_modify(NBas,dimOA,ints,NBas,work,&
                        Vabb(k,l)/nelA,Sxx,Sxx(l,k)/nelB,Vbaa)

       ic = l
       ia = k

       iac = posB(ia,ic)

       if(iac/=0) then

          do i=1,ADimX
             ip = A%IndN(1,i)
             ir = A%IndN(2,i)
             ipr = posA(ip,ir)
             
             ! 1A-1B
             termXX = termXX - 2d0*tpqrs(offA+ipr,offB+iac)*valTr*ints(ip,ir)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,iq)
             enddo

             ! 2A-1B
             termXX = termXX - 2d0*val*tpqrs(offA+ipr,offB+iac)*Amat(ip,ir)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(ip,iq)*Amat(iq,ir)
             enddo

             ! 3A-1B
             termXX = termXX + val*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,ir)*Amat(iq,ip)
             enddo

             ! 4A-1B
             termXX = termXX + val*tpqrs(offA+ipr,offB+iac)

          enddo
       endif

       ! 2B terms
       if(k==l) then

          ib = k
          Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + B%Occ(ib)*ints(1:NBas,1:dimOA)

       endif

       ib = l
       ia = k

       val3B = 0
       do iq=1,dimOA
          val3B = val3B + A%Occ(iq)*ints(iq,iq)
       enddo

       do ic=1,dimOB

          iac = posB(ia,ic)
 
          if(iac/=0) then

             do i=1,ADimX

                ip = A%IndN(1,i)
                ir = A%IndN(2,i)
                ipr = posA(ip,ir)

                ! 1A-3B
                termXX = termXX + ints(ip,ir)*Sbb(ic,ib)*tpqrs(offA+ipr,offB+iac)

                ! 2A-3B
                termXX = termXX + val3B*B%Occ(ib)*Sab(ip,ic)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ic)
                enddo

                ! 3A-3B
                termXX = termXX - val*B%Occ(ib)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)
 
                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ib)
                enddo

                ! 4A-3B
                termXX = termXX - val*B%Occ(ib)*Sab(ip,ic)*tpqrs(offA+ipr,offB+iac)

             enddo

          endif
       enddo

       if(k<=dimOB) then

          ib = k
          ic = l

          val4B = 0
          do iq=1,dimOA
             val4B = val4B + A%Occ(iq)*ints(iq,iq)
          enddo

          do ia=1,NBas

             iac = posB(ia,ic)
 
             if(iac/=0) then

                do i=1,ADimX

                   ip = A%IndN(1,i)
                   ir = A%IndN(2,i)
                   ipr = posA(ip,ir)

                   ! 1A-4B
                   termXX = termXX + B%Occ(ib)*Sbb(ia,ib)*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                   ! 2A-4B
                   termXX = termXX + val4B*B%Occ(ib)*Sab(ip,ib)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ib)
                   enddo 

                   ! 3A-4B
                   termXX = termXX - val*B%Occ(ib)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ia)
                   enddo 

                   ! 4A-4B
                   termXX = termXX - val*B%Occ(ib)*Sab(ip,ib)*tpqrs(offA+ipr,offB+iac)

                enddo
             endif

          enddo
       endif

    enddo
 enddo

 close(iunit)

 ! 2B terms
 do j=1,BDimX

    ia = B%IndN(1,j)
    ic = B%IndN(2,j)
    iac = posB(ia,ic) 

    do i=1,ADimX

       ip = A%IndN(1,i)
       ir = A%IndN(2,i)
       ipr = posA(ip,ir)

       ! 1A-2B 
       termXX = termXX - 2d0*Emat(ip,ir)*Sbb(ia,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,iq) 
       enddo

       ! 2A-2B 
       termXX = termXX - 2d0*val*Sab(ip,ic)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(ip,iq)*Sab(iq,ic)
       enddo

       ! 3A-2B 
       termXX = termXX + val*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,ir)*Sab(iq,ia)
       enddo

       ! 4A-2B 
       termXX = termXX + val*Sab(ip,ic)*tpqrs(offA+ipr,offB+iac)
 
    enddo
 enddo

 print*, 'A3-XX',termXX

 deallocate(ints,work) 

end subroutine term_X_A3_XX_SR

subroutine term_X_A3_XY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbaa,posA,posB,ADimX,BDimX,NBas,A,B,isXY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbaa(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: isXY
character(*),intent(in) :: IntKFile

integer :: i,j,k,l,kl
integer :: iq,ib,ip,ir,ia,ic,iac,ipr
integer :: iunit,offA,offB
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,valTr,val2A,val3B,val4B,termXY
double precision :: nelA,nelB
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas)
double precision :: Emat(NBas,NBas),Amat(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

 Amat = 0
 do ib=1,dimOB
    do ir=1,NBas 
       do ip=1,NBas
          Amat(ip,ir) = Amat(ip,ir) + B%Occ(ib)*Sab(ip,ib)*Sab(ir,ib)
       enddo
    enddo
 enddo
 
 valTr = 0
 do iq=1,dimOA
    valTr = valTr + A%Occ(iq)*Amat(iq,iq)
 enddo
 
 Sbb = 0
 do ic=1,NBas
    do ia=1,NBas
       do iq=1,dimOA
          Sbb(ia,ic) = Sbb(ia,ic) + A%Occ(iq)*Sab(iq,ia)*Sab(iq,ic)
       enddo
    enddo
 enddo

 if(isXY) then
    offA = 0
    offB = BDimX
 else
    offA = ADimX
    offB = 0 
 endif

!(FO|FO):(AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

 ! one loop over integrals
 termXY = 0
 Emat = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOA)

       call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/nelA,Sxx,Sxx(k,l)/nelB,Vbaa)

       ic = l
       ia = k

       iac = posB(ia,ic)

       if(iac/=0) then

          val2A = 0
          do iq=1,dimOA
             val2A = val2A + 2d0*A%Occ(iq)*ints(iq,iq)
          enddo
 
          do i=1,ADimX

             ip = A%IndN(1,i)
             ir = A%IndN(2,i)
             ipr = posA(ip,ir)

             ! 1A-1B
             termXY = termXY - 2d0*valTr*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

             ! 2A-1B
             termXY = termXY - val2A*Amat(ip,ir)*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(ip,iq)*Amat(iq,ir)
             enddo

             ! 3A-1B
             termXY = termXY + val*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,ir)*Amat(ip,iq)
             enddo

             ! 4A-1B
             termXY = termXY + val*tpqrs(offA+ipr,offB+iac)

          enddo
       endif

!       ! 2B terms 
       if(k==l.and.k<=dimOB) then

          ib = k
          Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + B%Occ(ib)*ints(1:NBas,1:dimOA)

       endif

       ia = k
       ib = l

       val4B = 0
       do iq=1,dimOA
          val4B = val4B + A%Occ(iq)*ints(iq,iq)
       enddo
       val4B = val4B*B%Occ(ib)

       do ic=1,dimOB

          iac = posB(ia,ic)

          if(iac/=0) then

             fact = B%Occ(ib)*Sbb(ic,ib)

             do i=1,ADimX

                ip = A%IndN(1,i)
                ir = A%IndN(2,i)
                ipr = posA(ip,ir)

                ! 1A-4B
                termXY = termXY + fact*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                ! 2A-4B
                termXY = termXY + val4B*Sab(ip,ib)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ib)
                enddo

                ! 3A-4B
                termXY = termXY - val*B%Occ(ib)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ic)
                enddo

                ! 4A-4B
                termXY = termXY - val*B%Occ(ib)*Sab(ip,ib)*tpqrs(offA+ipr,offB+iac)

             enddo
          endif

       enddo

       if(k<=dimOB) then

          ib = k
          ic = l

          val3B = 0
          do iq=1,dimOA
             val3B = val3B + A%Occ(iq)*B%Occ(ib)*ints(iq,iq)
          enddo

          do ia=1,NBas

             iac = posB(ia,ic)

             if(iac/=0) then

                fact = B%Occ(ib)*Sbb(ia,ib)

                do i=1,ADimX
 
                   ip = A%IndN(1,i)
                   ir = A%IndN(2,i)
                   ipr = posA(ip,ir)

                   ! 1A-3B
                   termXY = termXY + fact*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                   ! 2A-3B 
                   termXY = termXY + val3B*Sab(ip,ia)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ia)
                   enddo

                   ! 3A-3B 
                   termXY = termXY - val*B%Occ(ib)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ib)
                   enddo

                   ! 4A-3B 
                   termXY = termXY - val*B%Occ(ib)*Sab(ip,ia)*tpqrs(offA+ipr,offB+iac)

                enddo

             endif
          enddo

       endif

    enddo
 enddo

 close(iunit)

 ! 2B terms
 do j=1,BDimX
   
    ia = B%IndN(1,j)
    ic = B%IndN(2,j)
    iac = posB(ia,ic)

    fact = 2d0*Sbb(ia,ic)

    do i=1,ADimX

       ip = A%IndN(1,i)
       ir = A%IndN(2,i)
       ipr = posA(ip,ir)

       ! 1A-2B
       termXY = termXY - fact*Emat(ip,ir)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + 2d0*A%Occ(iq)*Emat(iq,iq)
       enddo 

       ! 2A-2B
       termXY = termXY - val*Sab(ip,ia)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(ip,iq)*Sab(iq,ia)
       enddo

       ! 3A-2B
       termXY = termXY + val*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,ir)*Sab(iq,ic)
       enddo

       ! 4A-2B
       termXY = termXY + val*Sab(ip,ia)*tpqrs(offA+ipr,offB+iac)

    enddo
 enddo

 print*, 'termXY',termXY

 deallocate(ints,work)

end subroutine term_X_A3_XY_SR

subroutine term_X_A2_XX_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,termXX,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|BA) or (AA|AB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ! one loop over integrals
 termXX = 0
 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)
 
       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(l,k)/nelA,Vabb,Vbab(l,k)/nelB,Sxx)
 
       ! 1 and 4
       if(k<=dimOB) then

          iq = k
          iu = l

          do it=1,NBas

             itu = posA(it,iu)
 
             fact = 2d0*Sab(it,iq)*B%Occ(iq)           

             if(itu/=0) then

                do i=1,BDimX 

                   ip = B%IndN(1,i)
                   ir = B%IndN(2,i)
                   ipr = posB(ip,ir)

                   if(trans) then
                      termXX = termXX + fact*tpqrs(ipr,itu)*ints(ip,ir)
                      termXX = termXX - tpqrs(ipr,itu)*ints(ip,iq)*B%Occ(iq)*Sab(it,ir)
                   else
                      termXX = termXX + fact*tpqrs(itu,ipr)*ints(ip,ir)
                      termXX = termXX - tpqrs(itu,ipr)*ints(ip,iq)*B%Occ(iq)*Sab(it,ir)
                   endif

                enddo
             endif 

          enddo
       endif

       ip = k
       iu = l

       do ir=1,dimOB
  
          ipr = posB(ip,ir)
         
          if(ipr/=0) then

             do it=1,NBas

                itu = posA(it,iu) 
                
                if(itu/=0) then
                 
                   ! 2  
                   val = 0
                   do iq=1,dimOB
                      val = val + ints(iq,ir)*B%Occ(iq)*Sab(it,iq)
                   enddo

                   if(trans) then
                      termXX = termXX - val*tpqrs(ipr,itu)
                   else
                      termXX = termXX - val*tpqrs(itu,ipr)
                   endif

                   ! 3
                   val = 0
                   do iq=1,dimOB
                      val = val + 2d0*ints(iq,iq)*B%Occ(iq)
                   enddo

                   if(trans) then
                      termXX = termXX + val*tpqrs(ipr,itu)*Sab(it,ir)
                   else
                      termXX = termXX + val*tpqrs(itu,ipr)*Sab(it,ir)
                   endif
                 
                endif   

             enddo
          endif
       enddo

    enddo
 enddo

 termXX = -0.5d0*termXX

 close(iunit)

 print*, 'termXX',termXX

 deallocate(ints,work)

end subroutine term_X_A2_XX_SR

subroutine term_X_A2_YY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,termYY,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 termYY = 0
 ints = 0
 work = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)
 
       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)
      
       it = k
       iq = l

       do iu=1,dimOA

          itu = posA(it,iu)
          fact = 2d0*B%Occ(iq)*Sab(iu,iq)

          if(itu/=0) then

             do i=1,B%NDimX

                ip = B%IndN(1,i)
                ir = B%IndN(2,i)
                ipr = posB(ip,ir)

                ! 1 and 4
                if(trans) then
                   termYY = termYY + fact*tpqrs(BDimX+ipr,ADimX+itu)*ints(ip,ir)
                   termYY = termYY - tpqrs(BDimX+ipr,ADimX+itu)*ints(ir,iq)*Sab(iu,ip)*B%Occ(iq)
                else
                   termYY = termYY + fact*tpqrs(ADimX+itu,BDimX+ipr)*ints(ip,ir)
                   termYY = termYY - tpqrs(ADimX+itu,BDimX+ipr)*ints(ir,iq)*Sab(iu,ip)*B%Occ(iq)
                endif

             enddo 

          endif 
       enddo

       it = k
       ir = l

       do ip=1,NBas

          ipr = posB(ip,ir)

          if(ipr/=0) then

             do iu=1,dimOA

                itu = posA(it,iu)

                if(itu/=0) then
 
                   val = 0
                   do iq=1,dimOB
                      val = val + ints(ip,iq)*Sab(iu,iq)*B%Occ(iq)
                   enddo

                   ! 2
                   if(trans) then
                      termYY = termYY - val*tpqrs(BDimX+ipr,ADimX+itu)
                   else
                      termYY = termYY - val*tpqrs(ADimX+itu,BDimX+ipr)
                   endif

                   val = 0
                   do iq=1,dimOB
                      val = val + 2d0*ints(iq,iq)*B%Occ(iq)
                   enddo

                   ! 3
                   if(trans) then
                      termYY = termYY + val*tpqrs(BDimX+ipr,ADimX+itu)*Sab(iu,ip)
                   else
                      termYY = termYY + val*tpqrs(ADimX+itu,BDimX+ipr)*Sab(iu,ip)
                   endif

                endif

             enddo
          endif
       enddo

    enddo
 enddo

 termYY = -0.5d0*termYY
 
 close(iunit)
 
 print*, 'termYY',termYY

 deallocate(ints,work)

end subroutine term_X_A2_YY_SR 

subroutine term_X_A2_XY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,termXY,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 ints = 0
 kl = 0
 
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       if(k<=dimOA) then

         call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)

         iu = k 
         iq = l

         do it=1,NBas

            itu = posA(it,iu)

            fact = 2d0*B%Occ(iq)*Sab(it,iq)

            if(itu/=0) then

               do i=1,BDimX

                  ip = B%IndN(1,i)
                  ir = B%IndN(2,i)
                  ipr = posB(ip,ir)

                  ! 1 and 4 
                  if(trans) then
                     termXY = termXY + fact*tpqrs(BDimX+ipr,itu)*ints(ip,ir)
                     termXY = termXY - tpqrs(BDimX+ipr,itu)*ints(ir,iq)*Sab(it,ip)*B%Occ(iq)
                  else
                     termXY = termXY + fact*tpqrs(itu,BDimX+ipr)*ints(ip,ir)
                     termXY = termXY - tpqrs(itu,BDimX+ipr)*ints(ir,iq)*Sab(it,ip)*B%Occ(iq)
                  endif

               enddo 

            endif
         enddo

         iu = k
         ir = l
 
         do ip=1,NBas

            ipr = posB(ip,ir)
 
            if(ipr/=0) then

               do it=1,NBas

                  itu = posA(it,iu)

                  if(itu/=0) then

                     val = 0
                     do iq=1,dimOB
                        val = val + B%Occ(iq)*ints(ip,iq)*Sab(it,iq)
                     enddo

                     ! 2
                     if(trans) then
                        termXY = termXY - val*tpqrs(BDimX+ipr,itu)
                     else
                        termXY = termXY - val*tpqrs(itu,BDimX+ipr)
                     endif

                     val = 0
                     do iq=1,dimOB
                        val = val + 2d0*B%Occ(iq)*ints(iq,iq)
                     enddo

                     ! 3
                     if(trans) then
                        termXY = termXY + val*tpqrs(BDimX+ipr,itu)*Sab(it,ip)
                     else
                        termXY = termXY + val*tpqrs(itu,BDimX+ipr)*Sab(it,ip)
                     endif
                    
                  endif
               enddo
            endif
         enddo

       endif

    enddo
 enddo

 termXY = -0.5d0*termXY

 close(iunit)

 print*, 'termXY',termXY

 deallocate(ints,work)

end subroutine term_X_A2_XY_SR

subroutine term_X_A2_YX_SR(dim1,dim2,tpqrs,IntKFile,IntJFile, &
                           Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile,IntJFile

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,termYX,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 termYX = 0
 work = 0
 ints = 0
 kl = 0
 
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)

       it = k
       iq = l

       do iu=1,dimOA

          itu = posA(it,iu)

          fact = 2d0*B%Occ(iq)*Sab(iu,iq)

          if(itu/=0) then

             do i=1,BDimX

                ip = B%IndN(1,i)
                ir = B%IndN(2,i)
                ipr = posB(ip,ir)

                ! 1 and 4
                if(trans) then
                   termYX = termYX + fact*tpqrs(ipr,ADimX+itu)*ints(ip,ir)
                   termYX = termYX - tpqrs(ipr,ADimX+itu)*ints(ip,iq)*Sab(iu,ir)*B%Occ(iq)
                else
                   termYX = termYX + fact*tpqrs(ADimX+itu,ipr)*ints(ip,ir)
                   termYX = termYX - tpqrs(ADimX+itu,ipr)*ints(ip,iq)*Sab(iu,ir)*B%Occ(iq)
                endif

             enddo

          endif
       enddo

    enddo
 enddo

 close(iunit)

!(FF|OO):(AB|BB) or (BA|AA)
 open(newunit=iunit,file=trim(IntJFile),status='OLD', &
      access='DIRECT',recl=8*NBas*NBas)

 work = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOB
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       call ints_modify(NBas,NBas,ints,NBas,work,&
                        Vabb(k,l)/nelA,Sab,Sxx(k,l)/nelB,Vbab)

       iq = k
       ir = l

       do ip=1,NBas

          ipr = posB(ip,ir)

          if(ipr/=0) then

             do i=1,ADimX

                it = A%IndN(1,i)
                iu = A%IndN(2,i)
                itu = posA(it,iu)

                if(trans) then
                   ! 2
                   termYX = termYX - B%Occ(iq)*tpqrs(ipr,ADimX+itu)*ints(it,ip)*Sab(iu,iq)
                   ! 3
                   if(iq==ir) then 
                      termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ipr,ADimX+itu)*ints(it,ip)*Sab(iu,ir)
                   endif
                else 
                   ! 2
                   termYX = termYX - B%Occ(iq)*tpqrs(ADimX+itu,ipr)*ints(it,ip)*Sab(iu,iq)
                   ! 3
                   if(iq==ir) then 
                      termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ADimX+itu,ipr)*ints(it,ip)*Sab(iu,ir)
                   endif
                endif

             enddo
          endif

       enddo

    enddo
 enddo

 close(iunit)

 termYX = -0.5d0*termYX

 print*, 'termYX',termYX

 deallocate(ints,work)

end subroutine term_X_A2_YX_SR

subroutine e2exdisp(Flags,A,B,SAPT)

implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT

integer :: NBas,NInte1
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
integer :: GdimOA,GdimOB
integer :: iunit
integer :: i,j,k,l,kl,ij,ii,jj,pq,rs
integer :: ip,iq,ir,is,ipq,irs
integer :: kc,ik,ic
logical :: approx
logical :: ipropab 
integer,allocatable :: posA(:,:),posB(:,:)
double precision,allocatable :: OmA(:), OmB(:)
double precision,allocatable :: EVecA(:), EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:),&
                                work(:),workSq(:,:),&
                                sij(:,:),&
                                tpqrs(:,:)
double precision :: nelA,nelB
double precision :: fact,val,termZ,termY,termX
!test
double precision :: fact2
double precision,allocatable :: ints2(:,:)
! 
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision,allocatable :: ints(:,:) 
double precision :: Tcpu,Twall
double precision,external  :: trace,FRDM2
! test for Be
!double precision,parameter :: SmallE = 1.D-1
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 

! set dimensions
 NBas = A%NBasis 
 dimOA = A%INAct+A%NAct
 GdimOA = A%num0+A%num1
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 GdimOB = B%num0+B%num1
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

! set e2exd_version
 !approx=.true.
 approx=.false.

! approximate RDM2
 allocate(A%Fmat(NBas,NBas),B%Fmat(NBas,NBas))
 call fill_Fmat(A%Fmat,A%Occ,NBas,1)
 call fill_Fmat(B%Fmat,B%Occ,NBas,1)

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          Sba(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas))
 allocate(tmp1(NBas,NBas),tmp2(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

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

 allocate(posA(NBas,NBas),posB(NBas,NBas))
 posA = 0
 do i=1,A%NDimX
    posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
 enddo
 posB = 0
 do i=1,B%NDimX
    posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
 enddo

 ! termZ
 ! PA(ab).S(bc)
 ! S(ad).PB(dc)
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 termZ = 2d0*SAPT%e2disp*termZ
 !write(LOUT,*) 'termZ ',termZ
 write(LOUT,'(1x,a,f16.8)') 'term Z      = ',  termZ*1.0d3

 deallocate(tmp2,tmp1)

 allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),&
          tmp3(A%NDimX,B%NDimX),work(nOVB),workSq(NBas,NBas))
 allocate(sij(A%NDimX,B%NDimX))

 ! term A1
 ! transform J and K
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFOOABAB','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOABBA','AOTWOSORT')
 ! term A2
 ! A2A(B): XX
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOBBBA','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAAAB','AOTWOSORT')
 ! A2A(B): YY
 call tran4_gen(NBas,&
          NBas,A%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          'FOFOBBAB','AOTWOSORT')
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAABA','AOTWOSORT')
 ! term A3
 call tran4_gen(NBas,&
          NBas,B%CMO,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          'FOFOAABB','AOTWOSORT')
 ! XY and YX, A2
 call tran4_gen(NBas,&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          B%num0+B%num1,B%CMO(1:NBas,1:(B%num0+B%num1)),&
          NBas,A%CMO,&
          NBas,B%CMO,&
          'FFOOABBB','AOTWOSORT')
 call tran4_gen(NBas,&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          A%num0+A%num1,A%CMO(1:NBas,1:(A%num0+A%num1)),&
          NBas,B%CMO,&
          NBas,A%CMO,&
          'FFOOBAAA','AOTWOSORT')


 deallocate(work)
 allocate(work(NBas**2),ints(NBas,NBas),ints2(NBas,NBas))

 ! A3: XX
 tmp1=0
 if(approx) then

    call app_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

    call inter_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                     B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

 ! TERMS XX, YY
 !(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ints = 0
 ints2 = 0
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
 call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XX')

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

    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Bval,     &
                     Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
                     A%Occ,B%IndN,A%IndN,posB,posA,     &
                     dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Aval,     &
                     Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
                     B%Occ,A%IndN,B%IndN,posA,posB,     &
                     dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !X_A.I.X_B
 call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,0d0,tmp3,A%NDimX)
! print*, 'Xa.I.Xb:',norm2(tmp3)

 tmp1 = sij

 ! A1:YY
 call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YY')

 ! A2:YY
 if(approx) then

    call app_A2_YY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
             A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

 !Y_A.I.Y_B
 call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Ya.I.Yb:',norm2(tmp3)

! TERMS XY, YX
! A3:XY
 tmp1 = 0
 if(approx) then

    call app_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 else

   call inter_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
                  B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)

 endif

! (FF|OO):(AB|AB)
 open(newunit=iunit,file='FFOOABAB',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 ints = 0
 ints2 = 0
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
 call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'XY')

! A2:XY
 if(approx) then

    call app_A2_XY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_YX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.,B)

 endif

 !X_A.I.Y_B
 call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Xa.I.Yb:',norm2(tmp3)

 tmp1 = sij

 ! A1:YX
 call A1_Mod_1el(tmp1,Vaab,Vbab,Sab,posA,posB,A,B,NBas,'YX')


 !A2: YX 
 if(approx) then

    call app_A2_YX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
    call app_A2_XY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 else

    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.,B)
    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)

 endif

! !Y_A.I.X_B
 call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Ya.I.Xb:',norm2(tmp3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WORKING VERSION XX, YY :: 
! ! TERMS XX, YY
! !(FO|FO):(AB|BA)
! open(newunit=iunit,file='FOFOABBA',status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
! tmp1 = 0
! ints = 0
! kl = 0
! do l=1,dimOA
!    do k=1,NBas
!       kl = kl + 1
!       read(iunit,rec=kl) work(1:NBas*dimOB)
!
!       ir = k
!       iq = l
!
!       do ip=1,NBas
!          do is=1,dimOB
!
!             ipq = posA(ip,iq)
!             irs = posB(ir,is)
!
!             if(ipq/=0.and.irs/=0) then
!
!                call ints_modify(NBas,dimOB,ints,NBas,work,Sab(iq,ir)/nelA,Vaab,Vbab(iq,ir)/nelB,Sab)
!
!                fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
!                            (B%Occ(ir)-B%Occ(is)) * &
!                            ints(ip,is)
!
!                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
! 
!             endif
!          enddo
!       enddo
! 
!    enddo
! enddo
! close(iunit)
!
! !X_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,0d0,tmp3,A%NDimX)
! print*, 'Xa.I.Xb:',norm2(tmp3)
!
! !(FO|FO):(AB|BA)
! open(newunit=iunit,file='FOFOABBA',status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
! tmp1 = 0
! ints = 0
! kl = 0
! do l=1,dimOA
!    do k=1,NBas
!       kl = kl + 1
!       read(iunit,rec=kl) work(1:NBas*dimOB)
!
!       ir = k
!       iq = l
!
!       do ip=1,NBas
!          do is=1,dimOB
!
!             ipq = posA(ip,iq)
!             irs = posB(ir,is)
!
!             if(ipq/=0.and.irs/=0) then
!
!                call ints_modify(NBas,dimOB,ints,NBas,work,Vaab(iq,ir)/nelA,Sab,Sab(iq,ir)/nelB,Vbab)
!
!                fact = -2d0*(A%Occ(ip)-A%Occ(iq)) * &
!                            (B%Occ(ir)-B%Occ(is)) * &
!                            ints(ip,is)
!
!                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
! 
!             endif
!          enddo
!       enddo
! 
!    enddo
! enddo
! close(iunit)
!
! !Y_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Ya.I.Yb:',norm2(tmp3)
!
! ! A3: XX
! tmp1=0
! if(approx) then
!
!    call app_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
!                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
!
! else
!
!    call inter_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
!                     B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
!
! endif
!
! !X_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Xa.I.Xb:',norm2(tmp3)
! !Y_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Ya.I.Yb:',norm2(tmp3)
!
! ! A2:XX
! tmp1 = 0
! if(approx) then
!
!    call app_A2_XX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat, &
!                   Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
!                   A%Occ,B%IndN,A%IndN,posB,posA,     &
!                   dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call app_A2_XX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat, & 
!                   Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
!                   B%Occ,A%IndN,B%IndN,posA,posB,     &
!                   dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! else
!
!    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Bval,     &
!                     Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
!                     A%Occ,B%IndN,A%IndN,posB,posA,     &
!                     dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call inter_A2_XX(A%NDimX,B%NDimX,tmp1,RDM2Aval,     &
!                     Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
!                     B%Occ,A%IndN,B%IndN,posA,posB,     &
!                     dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! endif
!
! !X_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Xa.I.Xb:',norm2(tmp3)
!
! ! A2:YY
! tmp1 = 0
! if(approx) then
!
!    call app_A2_YY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call app_A2_YY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! else
!
!    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
!             A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call inter_A2_YY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! endif
!
! !Y_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Ya.I.Yb:',norm2(tmp3)

!! WORKING VERSION: XY, YX
!! TERMS XY, YX
!! (FF|OO):(AB|AB)
! open(newunit=iunit,file='FFOOABAB',status='OLD', &
!     access='DIRECT',recl=8*NBas*NBas)
!
! ints = 0
! tmp1 = 0
! kl = 0
! do l=1,dimOB
!    do k=1,dimOA
!       kl = kl + 1
!       read(iunit,rec=kl) work(1:NBas*NBas)
!
!       iq = k
!       is = l
!
!       call ints_modify(NBas,NBas,ints,NBas,work,Sab(iq,is)/nelA,Vaab,Vbab(iq,is)/nelB,Sab)
!
!       do ir=1,NBas
!          do ip=1,NBas
!
!             ipq = posA(ip,iq)
!             irs = posB(ir,is)
!
!             if(ipq/=0.and.irs/=0) then
!
!
!                fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
!                           (B%Occ(ir)-B%Occ(is)) * &
!                           ints(ip,ir)
!
!                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
! 
!             endif
!            enddo
!         enddo
! 
!    enddo
! enddo
! close(iunit)
!
! !X_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Xa.I.Yb:',norm2(tmp3)
!!
! open(newunit=iunit,file='FFOOABAB',status='OLD', &
!     access='DIRECT',recl=8*NBas*NBas)
!
! ints = 0
! tmp1 = 0
! kl = 0
! do l=1,dimOB
!    do k=1,dimOA
!       kl = kl + 1
!       read(iunit,rec=kl) work(1:NBas*NBas)
!
!       iq = k
!       is = l
!
!       call ints_modify(NBas,NBas,ints,NBas,work,Vaab(iq,is)/nelA,Sab,Sab(iq,is)/nelB,Vbab)
!
!       do ir=1,NBas
!          do ip=1,NBas
!
!             ipq = posA(ip,iq)
!             irs = posB(ir,is)
!
!             if(ipq/=0.and.irs/=0) then
!
!
!                fact = 2d0*(A%Occ(ip)-A%Occ(iq)) * &
!                           (B%Occ(ir)-B%Occ(is)) * &
!                           ints(ip,ir)
!
!                tmp1(ipq,irs) = tmp1(ipq,irs) + fact
! 
!             endif
!            enddo
!         enddo
! 
!    enddo
! enddo
! close(iunit)
! !Y_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Ya.I.Xb:',norm2(tmp3)
!
!! A3:XY
! tmp1 = 0
! if(approx) then
!
!    call app_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
!                   B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
!
! else
!
!   call inter_A3_XY(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,RDM2Aval,RDM2Bval,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
!                  B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
!
! endif
!
! !X_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Xa.I.Yb:',norm2(tmp3)
! !Y_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
! print*, 'Ya.I.Xb:',norm2(tmp3)
!
!! A2:XY
! tmp1 = 0
! if(approx) then
!
!    call app_A2_XY(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call app_A2_YX(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! else
!
!    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA','FFOOBAAA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.,B)
!
! endif
!
! !X_A.I.Y_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigX,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigY,B%NDimX,1d0,tmp3,A%NDimX)
!!! print*, 'Xa.I.Yb:',norm2(tmp3)
!!
!! A2: YX 
! tmp1 = 0
! if(approx) then
!
!    call app_A2_YX(A%NDimX,B%NDimX,tmp1,B%Occ,B%Fmat,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!    call app_A2_XY(A%NDimX,B%NDimX,tmp1,A%Occ,A%Fmat,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! else
!
!    call inter_A2_YX(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBAB','FFOOABBB',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.,B)
!    call inter_A2_XY(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAABA',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
! endif
!
! !Y_A.I.X_B
! call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,A%EigY,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
! call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,B%EigX,B%NDimX,1d0,tmp3,A%NDimX)
!! print*, 'Ya.I.Xb:',norm2(tmp3)

 !! TEST sRef
 !allocate(tpqrs(2*A%NDimX,2*B%NDimX))
 !call prep_tpqrs(sij,tpqrs,posA,posB,A%NDimX,B%NDimX,NBas,A,B) 
 !call term_Y_SR(tpqrs,Sab,SAPT%elst,SAPT%Vnn,A%NDimX,B%NDimX,NBas,A,B)
 !call term_X_A1_SR(tpqrs,Sab,Vaab,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B)

 !call term_X_A2_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBBA',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !call term_X_A2_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAAAB',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.)

 !call term_X_A2_YY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !call term_X_A2_YY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.)

 !call term_X_A2_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !call term_X_A2_YX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA','FFOOBAAA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.)

 !call term_X_A2_YX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB','FFOOABBB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !call term_X_A2_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.)

 !call term_X_A3_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !call term_X_A3_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.true.)

 !call term_X_A3_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.true.)
 !call term_X_A3_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.)
 !
 !deallocate(tpqrs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call make_sij_Y(sij,tmp1,A,B,nOVB,NBas)

 endif

 print*, 'sij',norm2(sij)
 print*, 'tmp3-exd',norm2(tmp3)

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

! write(LOUT,*) 'termX ', termX
 write(LOUT,'(/1x,a,f16.8)') 'term X      = ',  termX*1.0d3

 ! termY
 ! term Y: t_ij 
 call make_tij_Y(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)

 ! term Y
 termY = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
        
      ! if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
      !    .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

       if(abs(A%Eig(i)).gt.SmallE.and.abs(B%Eig(j)).gt.SmallE&
          .and.abs(A%Eig(i)).lt.BigE.and.abs(B%Eig(j)).lt.BigE) then

          termY = termY + sij(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))

       endif
    enddo
 enddo

 termY = -8d0*(SAPT%elst-SAPT%Vnn)*termY
 write(LOUT,'(1x,a,f16.8)') 'term Y      = ',  termY*1.0d3

 write(LOUT,'(/1x,a,f16.8)') 'E2exch-disp = ', (termX+termY+termZ)*1.0d3
 deallocate(sij)

 deallocate(posB,posA,ints)
 deallocate(Vbab,Vaba,Vaab,Vbaa,Vabb,Vb,Va,PB,PA,Sab,S)
 deallocate(workSq,work,tmp2,tmp1)

end subroutine e2exdisp

subroutine e2exd_app(A,B,SAPT)

implicit none

type(SystemBlock),intent(in) :: A, B
type(SaptData),intent(in)    :: SAPT

integer :: iunit
integer :: i,j,k,l
integer :: ip,iq,ir,is,ipq,irs
integer :: NBas,dimOA,dimOB 
double precision :: nelA,nelB,fact,test
double precision :: termZ,termY,termX 
double precision,allocatable :: S(:,:),Sab(:,:),Sba(:,:),&
                                posA(:,:),posB(:,:)
double precision,allocatable :: PA(:,:),PB(:,:), &
                                PAbb(:,:),PBaa(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:), &
                                Vabb(:,:),Vbaa(:,:),&
                                Vaba(:,:),Vaab(:,:),Vbab(:,:)
double precision,allocatable :: RDM2Aval(:,:,:,:), &
                                RDM2Bval(:,:,:,:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),tmp3(:,:),&
                                EVecA(:),EVecB(:), &
                                EigA(:), EigB(:)
double precision,external  :: trace,FRDM2
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 


 ! dimensions
 NBas = A%NBasis 
 dimOA = A%INAct+A%NAct
 dimOB = B%INAct+B%NAct

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 allocate(S(NBas,NBas),Sab(NBas,NBas),&
          Sba(NBas,NBas),&
          PA(NBas,NBas),PB(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Vabb(NBas,NBas),Vbaa(NBas,NBas),&
          Vaba(NBas,NBas),Vaab(NBas,NBas),Vbab(NBas,NBas))
 allocate(tmp1(NBas,NBas),tmp2(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,1d0,PA)
 call get_den(NBas,B%CMO,B%Occ,1d0,PB)

 call get_one_mat('S',S,A%Monomer,NBas)
 call tran2MO(S,A%CMO,B%CMO,Sab,NBas)
 call tran2MO(S,B%CMO,A%CMO,Sba,NBas)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call tran2MO(Va,B%CMO,B%CMO,Vabb,NBas)
 call tran2MO(Vb,A%CMO,A%CMO,Vbaa,NBas)

 call tran2MO(Va,B%CMO,A%CMO,Vaba,NBas)
 call tran2MO(Va,A%CMO,B%CMO,Vaab,NBas)
 call tran2MO(Vb,A%CMO,B%CMO,Vbab,NBas)

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

 allocate(posA(NBas,NBas),posB(NBas,NBas))

 posA = 0
 do i=1,A%NDimX
    posA(A%IndN(1,i),A%IndN(2,i)) = A%IndX(i)
 enddo
 posB = 0
 do i=1,B%NDimX
    posB(B%IndN(1,i),B%IndN(2,i)) = B%IndX(i)
 enddo

 ! termZ
 tmp1 = 0
 tmp2 = 0
 termZ = 0
 call dgemm('N','N',NBas,NBas,NBas,1d0,PA,NBas,S,NBas,0d0,tmp1,NBas)
 call dgemm('N','N',NBas,NBas,NBas,1d0,S,NBas,PB,NBas,0d0,tmp2,NBas)
 do j=1,NBas
    do i=1,NBas
       termZ = termZ + tmp1(i,j)*tmp2(i,j)
    enddo
 enddo

 termZ = 2d0*SAPT%e2disp*termZ
 write(LOUT,*) 'termZ ',termZ

 deallocate(tmp2,tmp1)

 allocate(EVecA(A%NDimX*A%NDimX),EVecB(B%NDimX*B%NDimX))
 allocate(tmp1(A%NDimX,B%NDimX),tmp2(A%NDimX,B%NDimX),tmp3(A%NDimX,B%NDimX))
 allocate(EigA(A%NDimX),EigB(B%NDimX))

 ! termY
 tmp1 = 0
 do j=1,B%NDimX

    ir = B%IndN(1,j)
    is = B%IndN(2,j)
    irs = posB(ir,is)

    do i=1,A%NDimX

       ip = A%IndN(1,i)
       iq = A%IndN(2,i)
       ipq = posA(ip,iq)

       fact = (B%Occ(ir)-B%Occ(is)) * &
              (A%Occ(ip)-A%Occ(iq))

       tmp1(ipq,irs) = tmp1(ipq,irs) + Sab(iq,ir)*Sba(is,ip)*fact

    enddo
 enddo

 call readresp(EVecA,EigA,A%NDimX,'PROP_A')
 call readresp(EVecB,EigB,B%NDimX,'PROP_B')

 !Z_A.I.Z_B
 call dgemm('T','N',A%NDimX,B%NDimX,A%NDimX,1d0,EVecA,A%NDimX,tmp1,A%NDimX,0d0,tmp2,A%NDimX)
 call dgemm('N','N',A%NDimX,B%NDimX,B%NDimX,1d0,tmp2,A%NDimX,EVecB,B%NDimX,0d0,tmp3,A%NDimX)

 tmp2 = 0
 open(newunit=iunit,file='PROP_AB',form='UNFORMATTED',&
    access='SEQUENTIAL',status='OLD')
 read(iunit) tmp2
 close(iunit)

 ! term Y
 termY = 0d0
 do j=1,B%NDimX
    do i=1,A%NDimX
        
       if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
          .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          termY = termY + tmp3(i,j)/(A%Eig(i)+B%Eig(j))
         ! termY = termY + tmp2(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))
         ! termY = termY + sij(i,j)*tmp3(i,j)/(A%Eig(i)+B%Eig(j))

       endif
    enddo
 enddo

 termY = -8d0*(SAPT%elst-SAPT%Vnn)*termY
 write(LOUT,*) 'termY ', termY


 tmp1 = 0
 !! A2
!call inter_XX(A%NDimX,B%NDimX,tmp1,RDM2Bval,Sab,nelA,Vabb,nelB,Vbab,'FOFOBBBA',&
!              A%Occ,B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas,.false.)
!call inter_XX(A%NDimX,B%NDimX,tmp1,RDM2Aval,Sba,nelB,Vbaa,nelA,Vaba,'FOFOAAAB',&
!              B%Occ,A%IndN,B%IndN,posA,posB,dimOA,dimOB,A%NDimX,B%NDimX,NBas,.true.)
!
!tmp1 = -0.5d0*tmp1
!
!! A3
!call app_A3_XX(A%NDimX,B%NDimX,tmp1,A%Occ,B%Occ,A%Fmat,B%Fmat,Sab,nelA,Vabb,nelB,Vbaa,'FOFOAABB',&
!               B%IndN,A%IndN,posB,posA,dimOB,dimOA,B%NDimX,A%NDimX,NBas)
!
!
!test = 0
!
! do j=1,B%NDimX
!    do i=1,A%NDimX
!        
!       if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
!          .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then
!
!          test = test + tmp3(i,j)*sij(i,j)/(A%Eig(i)+B%Eig(j))
!
!       endif
!    enddo
! enddo
!
! print*, 'test',test 
!
!
deallocate(EVecB,EVecA,tmp1)
deallocate(posB,posA)
deallocate(RDM2Bval,RDM2Aval)
deallocate(Vbab,Vaab,Vaba,Vbaa,Vabb,Vb,Va,PB,PA,Sba,Sab,S)

end subroutine e2exd_app

subroutine make_VVmm(OccA,OccB,mat1,mat2,EigA,EigB,OutMat,IndNA,IndNB,NDimXA,NDimXB,NBasis,idx)
implicit none

integer,intent(in) :: NDimXA,NDimXB,NBasis
integer,intent(in) :: idx,IndNA(2,NDimXA),IndNB(2,NDimXB)
double precision,intent(in) :: OccA(NBasis),OccB(NBasis),&
                               mat1(NBasis,NBasis),mat2(NBasis,NBasis),&
                               EigA(NDimXA*NDimXA),EigB(NDimXB*NDimXB)
double precision,intent(out) :: OutMat(NDimXA,NDimXB)

integer :: ip,iq,ir,is,pq,rs,i,j
double precision :: fact
double precision,allocatable :: tmp1(:,:)

allocate(tmp1(NDimXA,NDimXB))

 !first big loop
 tmp1=0
 do pq=1,NDimXA
    ip = IndNA(1,pq)
    iq = IndNA(2,pq)
    do rs=1,NDimXB
       ir = IndNB(1,rs)
       is = IndNB(2,rs)

       if(idx==1) then
          ! qr | sp
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(iq,ir)*mat2(is,ip)
       elseif(idx==2) then
          ! pr | sq
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(ip,ir)*mat2(is,iq)
       elseif(idx==3) then
          ! ps | qr
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(ip,is)*mat2(iq,ir)
       elseif(idx==4) then
          ! qr | ps
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(iq,ir)*mat2(ip,is)
       elseif(idx==5) then
          ! pr | qs
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(ip,ir)*mat2(iq,is)
       elseif(idx==6) then
          ! qs | pr
          fact = (OccA(ip)-OccA(iq)) * &
                 (OccB(ir)-OccB(is)) * &
                  mat1(iq,is)*mat2(ip,ir)
       endif

       do i=1,NDimXA
          tmp1(i,rs) = tmp1(i,rs) + & 
                       fact * EigA(pq+(i-1)*NDimXA)
       enddo

    enddo
 enddo

 OutMat=0
 do j=1,NDimXB
    do i=1,NDimXB
       do rs=1,NDimXB
          ir = IndNB(1,rs)
          is = IndNB(2,rs)
          OutMat(i,j) = OutMat(i,j) + &
                        EigB(rs+(j-1)*NDimXB)*tmp1(i,rs)
       enddo
    enddo  
 enddo

deallocate(tmp1)

end subroutine make_VVmm

subroutine e2disp_pino(Flags,A,B,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: A, B
type(SaptData) :: SAPT
integer :: NBas, NInte1,NInte2
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ipq,irs
integer :: i1,i2,j1,j2,ii,jj,ij
integer :: coef,coef2,ADimEx,BDimEx
integer,allocatable :: AuxT(:,:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact,fpq,frs
double precision :: e2d,tmp
! testy
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
double precision,allocatable :: AVecEx(:),BVecEx(:)

!double precision,parameter :: SmallE = 1.d-6
double precision,parameter :: SmallE = 1.d-12
!double precision,parameter :: SmallE = 1.d-1

 write(LOUT,'(1x,a,e12.4)') 'SmallE in E2Disp PINO:',SmallE

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

 if(allocated(A%PP)) then
    print*, size(A%PP),ADimEX
 else
    print*, 'A%PP not allocated?!'
 endif
 if(allocated(A%PP)) then
    print*, size(B%PP),ADimEx
 endif

 do i=1,coef*ADimEx
    if(A%PP(i)<0d0) then 
       write(LOUT,'(1x,"Monomer A: Negative EVal",i4,f16.8)') i,OmA(i)
    endif
 enddo
 do i=1,coef*BDimEx
    if(B%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,OmB(i)
    endif
 enddo

 allocate(AuxT(2,NInte1))
 AuxT = 0
 ij = 0
 do j=1,NBas
    do i=1,j
       ij = ij + 1
       AuxT(1,ij) = j
       AuxT(2,ij) = i
    enddo
 enddo

 do i=1,coef*ADimEx
    if(A%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer A: Negative EVal",i4,f16.8)') i,A%PP(i)
      !!! test: Zero elements A
      ! do pq=1,ADimEx
      !     A%AP(i,pq) = 0d0 
      ! enddo

    endif
 enddo
 do i=1,coef*BDimEx
    if(B%PP(i)<0d0) then
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,B%PP(i)
      !!! test: Zero elements B
      ! do rs=1,BDimEx
      !     B%AP(i,rs) = 0d0 
      ! enddo
    endif
 enddo

 print*, 'DIMENSIONS ::',ij,NINte1,ADimEx

 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)

 allocate(tmp1(coef*ADimEx,coef*BDimEx),tmp2(coef*ADimEx,coef*BDimEx))

 tmp1=0
 do pq=1,ADimEx
    ip = AuxT(1,pq)
    iq = AuxT(2,pq)

    fpq = 1d0
    if(ip==iq) fpq=0.5d0

    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,BDimEx
       ir = AuxT(1,rs)
       is = AuxT(2,rs)

       frs = 1d0
       if(ir==is) frs=0.5d0

       fact = fpq*frs &
             * (A%CICoef(iq)*B%CICoef(is)+A%CICoef(ip)*B%CICoef(ir) &
             + A%CICoef(iq)*B%CICoef(ir)+A%CICoef(ip)*B%CICoef(is)) & 
             * work(is+(ir-1)*dimOB)

       do i=1,coef*ADimEx

          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*A%AP(i,pq)

          endif
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,coef*BDimEx

    if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
       do i=1,coef*ADimEx
          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
             do rs=1,BDimEx

                tmp2(i,j) = tmp2(i,j) + &
                     B%AP(j,rs)*tmp1(i,rs)

             enddo
          endif
       enddo
    endif
 enddo

 e2d = 0d0
 do i=1,coef*ADimEx

    if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
       do j=1,coef*BDimEx

          if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
             e2d = e2d + tmp2(i,j)**2/(A%PP(i)+B%PP(j))
          endif
       enddo
    endif
 enddo

 print*, ''
 print*, 'TEST: ',-16d0*e2d*1000d0
 SAPT%e2disp = -16d0*e2d
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d)*1000

 call writeampl(tmp2,'PROP_AB')

 deallocate(AuxT)
 deallocate(tmp2,tmp1)
 deallocate(work)

end subroutine e2disp_pino

subroutine e2dispCAS_pino(e2d,Flags,A,B,SAPT,ACAlpha)
implicit none

type(FlagsData)   :: Flags
type(SaptData)    :: SAPT
type(SystemBlock) :: A, B
double precision,intent(in)    :: ACAlpha
double precision,intent(inout) :: e2d

integer :: NBas, NInte1,NInte2
integer :: dimOA,dimFA,dimOB,dimFB,nOFA,nOFB
integer :: iunit
integer :: i,j,pq,rs
integer :: ip,iq,ir,is
integer :: ipq,irs
integer,allocatable          :: AuxT(:,:)
integer,allocatable          :: IGemA(:),IGemB(:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:),& 
                                tmp(:,:)
double precision,allocatable :: work(:)
double precision :: fact1,fact2
! testy
double precision,allocatable :: URe(:,:)
!
double precision,parameter :: SmallE = 1.d-12
!double precision,parameter :: SmallE = 1.d-1

 write(LOUT,'(1x,a,e12.4)') 'SmallE in E2Disp PINO:',SmallE

 if(A%NBasis.ne.B%NBasis) then
    write(LOUT,'(1x,a)') 'ERROR! MCBS not implemented in SAPT!'
    stop
 else
    NBas = A%NBasis 
 endif

 print*, 'Aalpha',norm2(A%AP),norm2(A%PP)
 print*, 'Balpha',norm2(B%AP),norm2(B%PP)

! set dimensions
 NInte1 = NBas*(NBas+1)/2
 dimOA  = NBas
 dimFA  = NBas
 dimOB  = NBas
 dimFB  = NBas 
 nOFA   = dimOA*dimFA
 nOFB   = dimOB*dimFB

 ! fix IGem for A
 allocate(IGemA(NBas),IGemB(NBas))
 do i=1,A%INAct
    IGemA(i) = 1
 enddo
 do i=A%INAct+1,A%INAct+A%NAct
    IGemA(i) = 2
 enddo
 do i=A%INAct+A%NAct+1,NBas
    IGemA(i) = 3
 enddo
 ! fix IGem for B
 do i=1,B%INAct
    IGemB(i) = 1
 enddo
 do i=B%INAct+1,B%INAct+B%NAct
    IGemB(i) = 2
 enddo
 do i=B%INAct+B%NAct+1,NBas
    IGemB(i) = 3
 enddo
 
 !do i=1,NBas
 !   print*, i,IGemA(i),IGemB(i)
 !enddo

 !print*, 'DIMENSIONS ::',NINte1,ADimEx,NBas**2

 ! transform AP
 allocate(tmp1(NInte1,NBas**2),tmp2(NInte1,NBas**2),&
          tmp(NInte1,NInte1))

 ! tran4_gen
 allocate(work(nOFB))
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOFB)
 tmp1 = 0
 tmp2 = 0
 do iq=1,NBas
    do ip=1,NBas
       pq=ip+(iq-1)*NBas
       !pq=iq+(ip-1)*NBas

       read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)

       do is=1,NBas
          do ir=1,NBas
             rs=ir+(is-1)*NBas
             !rs=is+(ir-1)*NBas

! the version below was our trial of computing dispersion which must be ADDED TO supermoleculr CAS
!             fact1 = work(is+(ir-1)*dimOB)
!             fact2 = ACAlpha*fact1
!
!             ! remove excitations form the same group
!              if(IGemA(ip)==IGemA(iq).and. &
!         IGemB(ir)==IGemB(is) .and. IGemA(ip)==IGemB(ir)) then
!                fact2 = fact1
!                fact1 = 0d0
!             endif
!
! dispersion included in supermolecular CAS is computed 
             if(IGemA(ip)==2.and. &
         IGemB(iq)==2.and.IGemA(ir)==2.and.IGemB(is)==2) then
                fact2 = work(is+(ir-1)*dimOB)
                fact1 = work(is+(ir-1)*dimOB)
             else
             fact1=0.d0
             fact2=0.d0
             endif
 

!!!!!!!!!!!!!!!!!!!!

             do i=1,NInte1
                tmp1(i,rs) =  tmp1(i,rs)+fact1*A%AP(i,pq)
             enddo

             do i=1,NInte1
                tmp2(i,rs) =  tmp2(i,rs)+fact2*A%AP(i,pq)
             enddo

          enddo
       enddo

    enddo
 enddo

 close(iunit)
 deallocate(work)

 ! second multiplication
 call dgemm('N','T',NInte1,NInte1,NBas**2,1d0,tmp1,NInte1,B%AP,NInte1,0d0,tmp,NInte1)
 tmp1=0
 tmp1(1:NInte1,1:NInte1) = tmp(1:NInte1,1:NInte1)
 call dgemm('N','T',NInte1,NInte1,NBas**2,1d0,tmp2,NInte1,B%AP,NInte1,0d0,tmp,NInte1)
 tmp2=0
 tmp2(1:NInte1,1:NInte1) = tmp(1:NInte1,1:NInte1)

 e2d = 0
 do j=1,NInte1
    if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
       do i=1,NInte1
          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
             e2d = e2d + tmp1(i,j)*tmp2(i,j)/(A%PP(i)+B%PP(j))
             !e2d = e2d + tmp2(i,j)**2/(A%PP(i)+B%PP(j))
          endif
       enddo 
    endif
 enddo

 SAPT%e2disp = -32d0*e2d
 e2d  = -32d0*e2d*1000d0
 write(LOUT,'(/,1x,a,f16.8)') 'E2disp(CAS) = ',e2d
! write(LOUT,'(/,1x,a,f16.8)') 'E2disp(t) = ',e2d/2d0

 deallocate(tmp,tmp2,tmp1)
 deallocate(IGemB,IGemA)

end subroutine e2dispCAS_pino

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
integer :: ipq,irs
integer :: i1,i2,j1,j2,ii,jj,ij
integer :: coef,coef2,ADimEx,BDimEx
integer,allocatable :: AuxT(:,:)
double precision,allocatable :: OmA(:),OmB(:)
double precision,allocatable :: EVecA(:),EVecB(:)
double precision,allocatable :: tmp1(:,:),tmp2(:,:) 
double precision,allocatable :: work(:)
double precision :: fact
double precision :: e2d1,e2d2,e2d3,e2d4,e2d,tmp
double precision :: e2du,dea,deb
! testy
integer,allocatable :: AIndEx(:,:),BIndEx(:,:)
double precision,allocatable :: AVecEx(:),BVecEx(:)

!double precision,parameter :: SmallE = 1.d-6
double precision,parameter :: SmallE = 1.d-12
!double precision,parameter :: SmallE = 1.d-1

 write(LOUT,'(/,1x,a)') 'Thresholds in E2disp-APSG:'
 write(LOUT,'(1x,a,e12.4)') 'SmallE in E2Disp PINO:',SmallE

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
 ! coef  = 2
 !coef2 = 4

 !print*, A%num0,A%num1,A%num2
 !print*, B%num0,B%num1,B%num2

 !print*, A%IGem
 !print*, ''
 !print*, B%IGem

! read EigValA_B
 allocate(EVecA(coef2*ADimEx**2),OmA(coef*ADimEx),&
          EVecB(coef2*BDimEx**2),OmB(coef*BDimEx))

 call readresp(EVecA,OmA,coef*ADimEx,'PROP_A')
 call readresp(EVecB,OmB,coef*BDimEx,'PROP_B')

! do i=1,coef*BDimEx
!    print*, 'OmB',i,OmB(i)
! enddo

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

! big test of a single loop
 allocate(AIndEx(2,ADimEx),BIndEx(2,BDimEx))
 allocate(AVecEx(coef2*ADimEx**2),BVecEx(coef2*BDimEx**2))

 AIndEx=A%IndNx
 BIndEx=B%IndNx

 ! prepare vectors
 AVecEx = 0
 BVecEx = 0
 do pq=1,ADimEx
    if(pq<=A%NDimX) then

        ip = A%IndN(1,pq)
        iq = A%IndN(2,pq)

        fact = A%CICoef(iq)+A%CICoef(ip)
        !fact = 1d0
        do i=1,coef*ADimEx
           AVecEx((i-1)*coef*ADimEx+pq) = fact * &
                                          EVecA((i-1)*coef*ADimEx+pq)
        enddo
 
    elseif(pq>A%NDimX) then

        ir = pq - A%NDimX
        !fact = 1d0
        fact = A%CICoef(ir)
        do i=1,coef*ADimEx
           AVecEx((i-1)*coef*ADimEx+pq) = fact * &
                                          EVecA((i-1)*coef*ADimEx+pq)
           !! for PINOVEC
           !AVecEx((i-1)*coef*ADimEx+pq) = fact * &
           !                               EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ir)
        enddo

    endif
 enddo

 do rs=1,BDimEx
     if(rs<=B%NDimX) then

        ir = B%IndN(1,rs)
        is = B%IndN(2,rs)

        fact = B%CICoef(ir)+B%CICoef(is)
        !fact = 1d0
        do i=1,coef*BDimEx
           BVecEx((i-1)*coef*BDimEx+rs) = fact * &
                                          EVecB((i-1)*coef*BDimEx+rs)
        enddo
 
    elseif(rs>B%NDimX) then

        ip = rs - B%NDimX
        fact = B%CICoef(ip)
        !fact = 1d0
        do i=1,coef*BDimEx
           BVecEx((i-1)*coef*BDimEx+rs) = fact * &
                                          EVecB((i-1)*coef*BDimEx+rs)
           !! for PINOVEC
           !BVecEx((i-1)*coef*BDimEx+rs) = fact * &
           !          EVecB((i-1)*coef*BDimEx+coef*B%NDimX+ip)
        enddo

    endif
 enddo

 do i=1,coef*ADimEx
    if(OmA(i)<0d0) then 
      ! !do pq=A%NDimX+1,ADimEx
      ! do pq=1,ADimEx
      !     AVecEx((pq-1)*coef*ADimEx+i) = 0d0 
      ! enddo
       write(LOUT,'(1x,"Monomer A: Negative EVal",i4,f16.8)') i,OmA(i)
    endif
 enddo
 do i=1,coef*BDimEx
    if(OmB(i)<0d0) then

       !!do rs=B%NDimX+1,BDimEx
       !do rs=1,BDimEx
       !    BVecEx((i-1)*coef*BDimEx+rs) = 0d0 
       !enddo
       write(LOUT,'(1x,"Monomer B: Negative EVal",i4,f16.8)') i,OmB(i)
    endif
 enddo

! good old stuff
 tmp1=0
 do pq=1,ADimEx
    ip = AIndEx(1,pq)
    iq = AIndEx(2,pq)
    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
    do rs=1,BDimEx
       ir = BIndEx(1,rs)
       is = BIndEx(2,rs)

       fact = work(is+(ir-1)*dimOB)

       do i=1,coef*ADimEx
          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then

             tmp1(i,rs) = tmp1(i,rs) + & 
                  fact*AVecEx((i-1)*coef*ADimEx+pq)

          endif
       enddo

    enddo
 enddo

 tmp2=0
 do j=1,coef*BDimEx
    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
       do i=1,coef*ADimEx
          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
             do rs=1,BDimEx
                ir = BIndEx(1,rs)
                is = BIndEx(2,rs)

                tmp2(i,j) = tmp2(i,j) + &
                     BVecEx((j-1)*coef*BDimEx+rs)*tmp1(i,rs)

             enddo
          endif
       enddo
    endif
 enddo

 e2d1 = 0d0
 do i=1,coef*ADimEx

    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
       do j=1,coef*BDimEx

          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
             e2d1 = e2d1 + tmp2(i,j)**2/(OmA(i)+OmB(j))
          endif
       enddo
    endif
 enddo

allocate(AuxT(2,NInte1))

AuxT = 0
ij = 0
do j=1,NBas
   do i=1,j
      ij = ij + 1
      AuxT(1,ij) = j
      AuxT(2,ij) = i
   enddo
enddo

!print*, 'ij:',ij,NINte1
!
!! try P matrices
! tmp1=0
! do pq=1,ADimEx
!    ip = AuxT(1,pq)
!    iq = AuxT(2,pq)
!
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,BDimEx
!       ir = AuxT(1,rs)
!       is = AuxT(2,rs)
!
!       fact = (A%CICoef(iq)*B%CICoef(is)+A%CICoef(ip)*B%CICoef(ir) &
!               + A%CICoef(iq)*B%CICoef(ir)+A%CICoef(ip)*B%CICoef(is))*work(is+(ir-1)*dimOB)
!       !fact = work(is+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!
!          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + & 
!                  !fact*AVecEx((i-1)*coef*ADimEx+pq)
!                  fact*A%AP(i,pq)
!
!          endif
!       enddo
!
!    enddo
! enddo
!
! tmp2=0
! do j=1,coef*BDimEx
!
!    if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!             do rs=1,BDimEx
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     !BVecEx((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!                     B%AP(j,rs)*tmp1(i,rs)
!
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
! e2d1 = 0d0
! do i=1,coef*ADimEx
!
!    if(abs(A%PP(i)).gt.SmallE.and.abs(A%PP(i)).lt.1d20) then
!       do j=1,coef*BDimEx
!
!          if(abs(B%PP(j)).gt.SmallE.and.abs(B%PP(j)).lt.1d20) then
!             e2d1 = e2d1 + tmp2(i,j)**2/(A%PP(i)+B%PP(j))
!          endif
!       enddo
!    endif
! enddo
!
 print*, ''
 print*, 'TEST: ',-16d0*e2d1*1000d0
 SAPT%e2disp = -16d0*e2d1
 write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000

 call writeampl(tmp2,'PROP_AB')

 e2d1=0

 deallocate(BVecEx,AVecEx)
 deallocate(BIndEx,AIndEx)

 deallocate(AuxT)

!! old code in parts
!! PART 1: p>q,r>s
! tmp1=0
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
!
!       fact = (A%CICoef(iq)+A%CICoef(ip)) * &
!              (B%CICoef(is)+B%CICoef(ir)) * work(is+(ir-1)*dimOB)
! !      print*, work(is+(ir-1)*dimOB)
!       !print*, fact
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + & 
!                  fact*EVecA((i-1)*coef*ADimEx+pq)
!
!          endif
!       enddo
!    enddo
! enddo
!
! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do rs=1,B%NDimX
!                ir = B%IndN(1,rs)
!                is = B%IndN(2,rs)
!                
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!             
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!!write(*,*) 'tmp1,tmp2:',norm2(tmp1),norm2(tmp2)
!
!! e2d1 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d1 = e2d1 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!             ! print*, OmA(i),OmB(j)
!!          endif
!!       enddo
!!    endif
!! enddo
!!
!!print*, ''
!!print*, 'PART1: ',-16d0*e2d1*1000d0
!!
!! !print*, 'test?',dimOB,B%NDimN
!
!! PART 2: p>q,r=s
! tmp1=0
! do pq=1,A%NDimX
!    ip = A%IndN(1,pq)
!    iq = A%IndN(2,pq)
!!   ! print*, iq,ip,iq+(ip-A%num0-1)*dimOA,nOFA
!    read(iunit,rec=iq+(ip-1)*dimOA) work(1:nOFB)
!    do ir=1,B%NDimN
!!      !print*, is,ir,is+(ir-B%num0-1)*dimOB,nOFB
!
!       fact = B%CICoef(ir)*(A%CICoef(iq)+A%CICoef(ip)) &
!            * work(ir+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,ir) = tmp1(i,ir) + &
!                  fact*EVecA((i-1)*coef*ADimEx+pq)
!
!          endif
!       enddo
!    enddo
! enddo
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do ir=1,B%NDimN
!          
!                tmp2(i,j) = tmp2(i,j) + &
!!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
!                     tmp1(i,ir)
!                
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! e2d2 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d2 = e2d2 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!          endif
!!       enddo
!!    endif
!! enddo
!! write(LOUT,*) 'PART2: ',-16d0*e2d2*1000d0
!
!! PART 3: p=q,r>s
! tmp1=0
! do ip=1,A%NDimN
!    !   ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
!    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
!    do rs=1,B%NDimX
!       ir = B%IndN(1,rs)
!       is = B%IndN(2,rs)
!       !      !print*, is,ir,is+(ir-1)*dimOB,nOFB
!       fact = A%CICoef(ip) * &
!            (B%CICoef(is)+B%CICoef(ir)) * &
!            work(is+(ir-1)*dimOB)
!
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,rs) = tmp1(i,rs) + &
!                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)
!
!          endif
!       enddo
!    enddo
! enddo
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do rs=1,B%NDimX
!                ir = B%IndN(1,rs)
!                is = B%IndN(2,rs)
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+rs)*tmp1(i,rs)
!
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! e2d3 = 0d0
!! do i=1,coef*ADimEx
!!    !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!!    if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!!       do j=1,coef*BDimEx
!!          !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!!          if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!!             e2d3 = e2d3 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!!          endif
!!       enddo
!!    endif
!! enddo
!! write(LOUT,*) 'PART3: ', -16d0*e2d3*1000d0
!
!! PART 4: p=q,r=s
! tmp1=0
! do ip=1,A%NDimN
!    ! print*, iq,ip,iq+(ip-1)*dimOA,nOFA
!    read(iunit,rec=ip+(ip-1)*dimOA) work(1:nOFB)
!    do ir=1,B%NDimN
!       !print*, is,ir,is+(ir-1)*dimOB,nOFB
!       fact = A%CICoef(ip)*B%CICoef(ir)* &
!              work(ir+(ir-1)*dimOB)
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!
!             tmp1(i,ir) = tmp1(i,ir) + &
!                  fact*EVecA((i-1)*coef*ADimEx+coef*A%NDimX+ip)
!
!          endif
!       enddo
!    enddo
! enddo
! 
!! print*, 'tmp1', norm2(tmp1(1:2*ADimEx,1:B%NDimN))
!!
!! tmp2=0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             do ir=1,B%NDimN
!
!                tmp2(i,j) = tmp2(i,j) + &
!                     EVecB((j-1)*coef*BDimEx+coef*B%NDimX+ir) * &
!                     tmp1(i,ir)
!                
!             enddo
!          endif
!       enddo
!    endif
! enddo
!
!! print*, 'tmp2: ',norm2(tmp2)
!
! e2d4 = 0d0
! do j=1,coef*BDimEx
!    !if(OmB(j).gt.SmallE.and.OmB(j).lt.1d20) then
!    if(abs(OmB(j)).gt.SmallE.and.abs(OmB(j)).lt.1d20) then
!       do i=1,coef*ADimEx
!          !if(OmA(i).gt.SmallE.and.OmA(i).lt.1d20) then
!          if(abs(OmA(i)).gt.SmallE.and.abs(OmA(i)).lt.1d20) then
!             e2d4 = e2d4 + tmp2(i,j)**2/(OmA(i)+OmB(j))
!            ! print*, tmp2(i,j)**2,OmA(i),OmB(j)
!          endif
!       enddo
!    endif
! enddo
! print*, 'PART4: ',-16*e2d4*1000

! SAPT%e2disp = -16d0*e2d
! write(LOUT,'(1x,a,f16.8)') 'E2disp      = ',-16*(e2d1+e2d2+e2d3+e2d4)*1000

 close(iunit)
 deallocate(work)
 deallocate(tmp2,tmp1)
 deallocate(OmB,EVecB,OmA,EVecA)

end subroutine e2disp_apsg

subroutine ModABMin_old(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES

implicit none

integer,intent(in) :: NBasis,NDimX,NGrid,NInte2
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               !Wt(NGrid),OrbGrid(NBasis,NGrid)
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(in) :: TwoNO(NInte2),TwoElErf(NInte2)
double precision,intent(inout) :: ABMin(NDimX**2)

double precision :: CICoef(NBasis)
double precision,allocatable :: work(:)
integer :: i,IRow,ICol,ia,ib,iab,ic,id,icd
double precision :: XKer1234,TwoSR,CA,CB,CC,CD
integer,external :: NAddr3,NAddrrK

print*, 'MOD-OLD!!!!'

allocate(work(NGrid))

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do i=1,NGrid
   work(i) = Wt(i)*SRKer(i)
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

      XKer1234 = 0
      do i=1,NGrid
         !XKer1234 = XKer1234 + Wt(i)*SRKer(i)* &
         XKer1234 = XKer1234 + work(i)* &
         !OrbGrid(ia,i)*OrbGrid(ib,i)*OrbGrid(ic,i)*OrbGrid(id,i)
         OrbGrid(i,ia)*OrbGrid(i,ib)*OrbGrid(i,ic)*OrbGrid(i,id)
      enddo
      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))

      ABMin((ICol-1)*NDimX+IRow)=ABMin((ICol-1)*NDimX+IRow) &
                       +4.0d0*(CA+CB)*(CD+CC)*(XKer1234+TwoSR)

   enddo
enddo

deallocate(work)

end subroutine ModABMin_old

subroutine ModABMin(Occ,SRKer,Wt,OrbGrid,TwoNO,TwoElErf,ABMin,IndN,IndX,NDimX,NGrid,NInte2,NBasis)
!     ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
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
         ABKer(icd,iab) = ABKer(iab,icd)
      
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
      !if(icd.gt.iab) cycle
 
      TwoSR=TwoNO(NAddr3(ia,ib,ic,id))-TwoElErf(NAddr3(ia,ib,ic,id))
      
      ABMin(iab,icd) = ABMin(iab,icd) &
                       +4.0d0*(CA+CB)*(CD+CC)*(ABKer(iab,icd)+TwoSR)
      !ABMin(icd,iab) = ABMin(iab,icd)   

   enddo
enddo

deallocate(ABKer,batch,work)

end subroutine ModABMin

subroutine ModABMin_mithap(Occ,SRKer,Wt,OrbGrid,ABMin,IndN,IndX,NDimX,NGrid,NBasis,&
                           twofile,twoerfile)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
character(*),intent(in) :: twofile,twoerfile
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

integer :: offset,batchlen,iunit1,iunit2
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

ABKer = 0

!print*, 'ModABMin_mithap'

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ip = IndN(1,IRow)
      iq = IndN(2,IRow)
      ipq = IndX(IRow)
   
      do ICol=1,NDimX
         ir=IndN(1,ICol)
         is=IndN(2,ICol)
         irs=IndX(ICol)
         if(irs.gt.ipq) cycle
    
         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
            batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo
         
         ABKer(ipq,irs) = ABKer(ipq,irs) + XKer1234
         ABKer(irs,ipq) = ABKer(ipq,irs)
      
      enddo
   enddo

enddo

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
                !if(irs.gt.ipq) cycle

                TwoSR = ints1(i,j)-ints2(i,j)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                !ABMIN(irs,ipq) = ABMIN(ipq,irs) 

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

close(iunit1)
close(iunit2)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_mithap

subroutine ACEneERPA_FFFF(ECorr,EVec,EVal,Occ,IGem, &
                          IndN,IndX,NOccup,NDimX,NBasis,IntFile)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IGem(NBasis),IndN(2,NDimX),IndX(NDimX)
integer,intent(in) :: NOccup
character(*),intent(in) :: IntFile
double precision,intent(out) :: ECorr
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)
double precision :: Occ(NBasis)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
logical :: AuxCoeff(3,3,3,3)
double precision :: CICoef(NBasis),Cpq,Crs,SumY,Aux
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 1.d8

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = .false.
            endif
         enddo
      enddo
   enddo
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0

! FULL INTS
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work,ints,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)

                !if(.not.(IGem(ir).eq.IGem(is).and.IGem(ip).eq.IGem(iq)&
                !.and.IGem(ir).eq.IGem(ip)).and.ir.gt.is.and.ip.gt.iq) then
                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then

                   ISkippedEig = 0
                   SumY = 0 
                   do kk=1,NDimX
                      if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                         SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                      else
                         ISkippedEig = ISkippedEig + 1
                         Skipped(ISkippedEig) = EVal(kk)
                      endif
                   enddo
          
                   Aux = 2*Crs*Cpq*SumY
          
                   if(iq.Eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif  

                   ECorr = ECorr + Aux*ints(j,i)
          
                ! endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
                endif

              endif
           enddo
        enddo

      endif

   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

deallocate(Skipped)
deallocate(ints,work)

end subroutine ACEneERPA_FFFF

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
!double precision,intent(in) :: EVec(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN)), &
!                               EVal(2*(Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: EVec((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN)), &
                               EVal((Mon%NDimX+Mon%NdimN))
double precision,intent(in) :: Freq
!double precision,intent(out) :: Alpha(2*(Mon%NDimX+Mon%NDimN),2*(Mon%NDimX+Mon%NDimN))
double precision,intent(out) :: Alpha((Mon%NDimX+Mon%NDimN),(Mon%NDimX+Mon%NDimN))
double precision :: frac
integer :: pq,rs,t,ip,iq,ir,is
integer :: NDimEx
double precision,parameter :: SmallE = 1.d-6

 NDimEx = (Mon%NDimX + Mon%NDimN)

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

subroutine create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig,CICoef,IndN,nblk,NBasis,NDimX)
implicit none

type(EBlockData)   :: Sblock(nblk),SblockIV
type(Y01BlockData) :: Y01Block(NDimX)

integer,intent(in)           :: nblk,NBasis,NDimX
integer,intent(in)           :: IndN(2,NDimX)
double precision,intent(in)  :: CICoef(NBasis)
double precision,intent(out) :: Eig(NDimX)

integer          :: i,ii,ipos,iblk,ip,iq
double precision :: fact

do iblk=1,nblk
   associate(A => SBlock(iblk))

     do i=1,A%n
        ipos = A%pos(i)

        ip   = IndN(1,ipos)
        iq   = IndN(2,ipos)
        fact = CICoef(ip)-CICoef(iq)
        associate(Y => Y01Block(ipos))

          Y%n  = A%n
          Y%l1 = A%l1
          Y%l2 = A%l2
          allocate(Y%vec0(Y%n))
          Y%vec0(1:Y%n) = fact*(A%matY(i,1:Y%n)-A%matX(i,1:Y%n))

        end associate
     enddo
     Eig(A%l1:A%l2) = A%vec(1:A%n)
   end associate
enddo
associate(A => SBlockIV)

  do i=1,A%n
     ii = A%l1+i-1
     ipos = A%pos(i)

     associate(Y => Y01Block(ipos))

       Y%n = 1
       Y%l1 = ii
       Y%l2 = ii
       allocate(Y%vec0(1))
       Y%vec0(1) = 1d0/sqrt(2d0)

     end associate

     Eig(ii) = A%vec(i)
  enddo

end associate

end subroutine create_Y0tilde

subroutine read_SBlock(SBlock,SBlockIV,nblk,xy0file)
implicit none

type(EBlockData)             :: SBlockIV
type(EBlockData),allocatable :: SBlock(:)
integer,intent(out)          :: nblk
character(*),intent(in)      :: xy0file

integer                      :: iunit,i,ival,iblk

! read SBLOCK
open(newunit=iunit,file=xy0file,status='OLD',&
     form='UNFORMATTED')
read(iunit) nblk
!print*, 'ival,nblk',ival,nblk
allocate(SBlock(nblk))
do iblk=1,nblk
   associate(A => SBlock(iblk))
     read(iunit) i, A%n, A%l1, A%l2
     allocate(A%pos(A%n),A%matX(A%n,A%n),A%matY(A%n,A%n),A%vec(A%n))
     read(iunit) A%pos,A%matX,A%matY,A%vec
   end associate
enddo
associate(A => SBlockIV)
  read(iunit) A%n,A%l1,A%l2
  allocate(A%pos(A%n),A%vec(A%n))
  allocate(A%matX(A%n,1),A%matY(A%n,1))
  read(iunit) A%pos,A%vec
  read(iunit) A%matX,A%matY
end associate
close(iunit)

end subroutine read_SBlock

subroutine convert_XY0_to_Y01(Mon,Y01Block,Eig0,NBasis,xy0file)
implicit none

type(SystemBlock)              :: Mon
type(EblockData)               :: SBlockIV
type(EblockData),allocatable   :: SBlock(:)
type(Y01BlockData)             :: Y01Block(Mon%NDimX)

integer,intent(in)           :: NBasis
character(*),intent(in)      :: xy0file
double precision,intent(out) :: Eig0(Mon%NDimX)

integer :: iunit
integer :: i,iblk,nblk,ival,ipos

call read_SBlock(SBlock,SBlockIV,nblk,xy0file)

call create_Y0tilde(Y01Block,SBlock,SBlockIV,Eig0,Mon%CICoef,Mon%IndN,nblk,NBasis,Mon%NDimX)

! deallocate SBLOCK
 do iblk=1,nblk
    associate(A => SBlock(iblk))
      deallocate(A%matY,A%matX,A%vec)
      deallocate(A%pos)
    end associate
 enddo
 ! deallocate IV part 
 associate(A => SblockIV)
   if(A%n>0) then
      deallocate(A%vec)
      deallocate(A%pos)
   endif
 end associate

end subroutine convert_XY0_to_Y01

subroutine writeampl(sij,ampfile)
implicit none

character(*) :: ampfile
double precision :: sij(:,:)
integer :: iunit

 open(newunit=iunit,file=ampfile,form='unformatted')
 write(iunit) sij 
 close(iunit)

end subroutine writeampl

subroutine check_loop(ABPLUS,Occ,IndN,IndBlock,NAct,INActive,NDim,NDimX,NBasis)
implicit none

integer,intent(in) :: NAct,INActive,NDim,NDimX,NBasis
integer,intent(in) :: IndN(2,NDim),IndBlock(2,NDimX)
double precision,intent(in) :: ABPLUS(NDimX,NDimX),Occ(NBasis)

integer :: i,j,k,l,ip,iq,ir,is,ipq,irs
integer :: NOccup 
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis),Cpq,Crs,EAll,EIntra,Aux
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
integer :: iunit
integer :: ii,jj,kl,ipos
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)

! set dimensions
NOccup = NAct + INActive
Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
   pos(IndBlock(1,i),IndBlock(2,i)) = i
enddo

EAll = 0
EIntra = 0
open(newunit=iunit,file='FOFO',status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

kl = 0
do k=1,NOccup
   do l=1,NBasis
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work1(1:NBasis*NOccup)
        do j=1,NOccup
           do i=1,NBasis
              ints(i,j) = work1((j-1)*NBasis+i)
           enddo
        enddo
        ints(:,NOccup+1:NBasis) = 0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = C(ir)+C(is)
                Cpq = C(ip)+C(iq)
!
                Aux = Crs*Cpq*ABPLUS(irs,ipq)
                EAll = EAll + Aux*ints(j,i)

                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

print*, 'FAKE-LOOP',EAll,EIntra
print*, '',EAll-EIntra

deallocate(ints,work2,work1)

end subroutine check_loop

end module sapt_ener
