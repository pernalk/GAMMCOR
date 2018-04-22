module sapt_main
use types
use systemdef
use tran
use sorter
use sapt_ener

implicit none

contains

subroutine sapt_driver(Flags,SAPT)
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT
integer :: i

! TEMPORARY - JOBTYPE_2
! ERPA
 Flags%IFlAC  = 0
 Flags%IFlSnd = 0

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'STARTING SAPT CALCULATIONS'
 write(LOUT,'(8a10)') ('**********',i=1,8)
 
 call sapt_interface(Flags,SAPT)
 write(LOUT,'()')
 call e1elst(SAPT%monA,SAPT%monB,SAPT)
 call e2ind(Flags,SAPT%monA,SAPT%monB,SAPT)

 if(Flags%ISERPA==0) then
    call e2disp(Flags,SAPT%monA,SAPT%monB,SAPT)
 elseif(Flags%ISERPA==2) then
    call e2disp_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)
 endif

 call print_warn(SAPT)
 call free_sapt(SAPT)

end subroutine sapt_driver

subroutine sapt_interface(Flags,SAPT) 
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT

integer :: NBasis,NSq,NInte1,NInte2
integer :: NCMOt, NOrbt, NBasist 
integer :: NSym, NOrb
integer :: NOcc(8),NOrbs(8)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: Ha(:),Hb(:)
double precision,allocatable :: Va(:),Vb(:),S(:)
double precision,allocatable :: Ca(:),Cb(:)
double precision :: potnucA,potnucB
integer :: ione,iorb,isiri,i,j,ij
logical :: exsiri
double precision :: tmp
integer :: p,q
!double precision,allocatable :: work3(:,:)
double precision,allocatable :: work3(:)
integer :: tmp1
double precision ::  potnuc,emy,eactiv,emcscf
integer :: noccA, nvirtA, noccB, nvirtB
integer :: ncen
!integer,parameter :: maxcen = 500
!double precision :: charga(maxcen),xyza(maxcen,3)
!double precision :: chargb(maxcen),xyzb(maxcen,3)

! read basis info
! only DCBS 
 NBasis = 0
 call basinfo(NBasis,'SIRIUS_A.RST')
 if(NBasis==0.and.SAPT%monA%NBasis==0) then
    write(LOUT,'(1x,a)') 'ERROR!!! NBasis NOWHERE TO BE FOUND!'
    stop
 elseif(NBasis==0.and.SAPT%monA%NBasis/=0) then
    ! basis only in input
    NBasis = SAPT%monA%NBasis
 elseif(NBasis/=0.and.SAPT%monA%NBasis==0) then
    ! basis only in SIRIFC
    SAPT%monA%NBasis = NBasis
    SAPT%monB%NBasis = NBasis
 elseif(NBasis/=0.and.SAPT%monA%NBasis/=0) then
    ! choose SIRIFC values
    SAPT%monA%NBasis = NBasis
    SAPT%monB%NBasis = NBasis
 endif

! set dimensions
 NSq = NBasis**2
 NInte1 = NBasis*(NBasis+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 SAPT%monA%NDim = NBasis*(NBasis-1)/2
 SAPT%monB%NDim = NBasis*(NBasis-1)/2

 allocate(work1(NInte1),work2(NBasis**2))
 allocate(Ha(NBasis**2),Hb(NBasis**2))
 allocate(Va(NBasis**2),Vb(NBasis**2),S(NBasis**2))
! read and dump 1-electron integrals 
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
 open(newunit=ione,file='AOONEINT_A',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NOrb,SAPT%monA%PotNuc
 
 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call triang_to_sq(work1,Ha,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call triang_to_sq(work1,work2,NBasis)
 Va(:) = Ha - work2

 call readlabel(ione,'OVERLAP ')
 call readoneint(ione,work1)
 call triang_to_sq(work1,S,NBasis)

 call readlabel(ione,'ISORDK  ')
 read(ione) 
 read(ione) SAPT%monA%charg,ncen,SAPT%monA%xyz 

! print*, ncen
! write(LOUT,*) SAPT%monA%charg(1:ncen)
! write(LOUT,*) SAPT%monA%xyz(1,:)
! write(LOUT,*) SAPT%monA%xyz(2,:)

! write(*,*) 'VA'
! call print_sqmat(Va,NBasis)
! call print_diag(Va,NBasis)

! tmp = 0d0
! do i=1,NDim
!    tmp = tmp + S(i)**2
! enddo
! write(*,*) tmp

 close(ione)
 ! square form
 call writeoneint('ONEEL_A',NSq,S,Va,Ha)
 
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

 open(newunit=ione,file='AOONEINT_B',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NOrb,SAPT%monB%PotNuc

 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call triang_to_sq(work1,Hb,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call triang_to_sq(work1,work2,NBasis)
 Vb(:) = Hb - work2

 call readlabel(ione,'ISORDK  ')
 read(ione) 
 read(ione) SAPT%monB%charg,ncen,SAPT%monB%xyz 

! print*, ncen
! write(LOUT,*) SAPT%monB%charg(1:ncen)
! write(LOUT,*) SAPT%monB%xyz(1,:)
! write(LOUT,*) SAPT%monB%xyz(2,:)

! rearrange in V: (B,A) -> (A,B)
 call arrange_oneint(S,NBasis,SAPT)
 call arrange_oneint(Vb,NBasis,SAPT)
 call arrange_oneint(Hb,NBasis,SAPT)

 !print*, 'VB:'
 !call print_sqmat(Vb,NBasis)

 close(ione)
 ! square form
 call writeoneint('ONEEL_B',NSq,S,Vb,Hb)

! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! read coefficient, occupancies
 
 inquire(file='SIRIFC_A',EXIST=exsiri)
 if(exsiri) then
    open(newunit=isiri,file='SIRIFC_A',status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(isiri,'TRCCINT ')
    read(isiri) NSym,NOrbt,NBasist,NCMOt,NOcc(1:NSym),NOrbs(1:NSym)
!    print*, 'NSym,NOrb,NOcc: '
!    write(*,*) NSym,NOrbt,NBasist,NCMOt
    SAPT%monA%NSym = NSym
    SAPT%monA%NOrb = NOrbt
    SAPT%monA%GFunc(1:NSym) = NOrbs(1:NSym)
 else
    SAPT%monA%NOrb = NOrb
    SAPT%monA%GFunc(1) = NOrb
    NBasist = NBasis
    NCMOt = NOrb*NBasis
 endif 

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
 
 if(Flags%ICASSCF==1.and.Flags%ISHF==0) then
    ! CASSCF
    call readmulti(NBasis,SAPT%monA,.false.,exsiri,isiri,'occupations_A.dat','SIRIUS_A.RST')
 elseif(Flags%ICASSCF==1.and.Flags%ISHF==1) then
    ! Hartree-Fock
    call readmulti(NBasis,SAPT%monA,.true.,exsiri,isiri,'occupations_A.dat','SIRIUS_A.RST')
    call readener(NBasis,SAPT%monA,isiri)
 elseif(Flags%IGVB==1) then 
    ! GVB
    call readgvb(SAPT%monA,NBasis,'coeff_A.dat')
 endif

 if(exsiri) close(isiri) 

! BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
 inquire(file='SIRIFC_B',EXIST=exsiri)
 if(exsiri) then
    open(newunit=isiri,file='SIRIFC_B',status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(isiri,'TRCCINT ')
    read(isiri) NSym,NOrbt,NBasist,NCMOt,NOcc(1:NSym),NOrbs(1:NSym)
 !   write(*,*) 'B: ',NSym,NOrbt,NBasist,NCMOt
 !   write(*,*) NOrb, 'NOrb'
    SAPT%monB%NSym = NSym
    SAPT%monB%NOrb = NOrbt
    SAPT%monB%GFunc(1:NSym) = NOrbs(1:NSym)
 else
    SAPT%monB%NOrb = NOrb
    SAPT%monB%GFunc(1) = NOrb
    
!    NBasist = NBasis
!    NCMOt = NOrb*NBasis
 endif 

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
!    write(*,*)   potnuc,emy,eactiv,emcscf

 if(Flags%ICASSCF==1.and.Flags%ISHF==0) then
    ! CASSCF
    call readmulti(NBasis,SAPT%monB,.false.,exsiri,isiri,'occupations_B.dat','SIRIUS_B.RST')
 elseif(Flags%ICASSCF==1.and.Flags%ISHF==1) then
    ! Hartree-Fock
    call readmulti(NBasis,SAPT%monB,.true.,exsiri,isiri,'occupations_B.dat','SIRIUS_B.RST')
    call readener(NBasis,SAPT%monB,isiri)
 elseif(Flags%IGVB==1) then
    ! GVB 
    call readgvb(SAPT%monB,NBasis,'coeff_B.dat')
 endif

 if(exsiri) close(isiri) 
 call print_occ(NBasis,SAPT,Flags%ICASSCF)

! read orbitals 
! norb.leq.nbas, orbitals mays be deleted due to linear
! dependecies in large basis sets
! ncmot = norb*nbas
 !allocate(Ca(NCMOt),Cb(NCMOt))
 allocate(Ca(NBasis*NBasis),Cb(NBasis*NBasis))
! allocate(work3(NBasis,NBasis))

 call read_mo(Ca,NOrbt,NBasist,'SIRIUS_A.RST','DALTON_A.MOPUN')
 call read_mo(Cb,NOrbt,NBasist,'SIRIUS_B.RST','DALTON_B.MOPUN')

 call arrange_mo(Cb,NBasist,SAPT%monA,SAPT%monB)

! MAYBE: one should print with NOrbt?
! if(SAPT%IPrint.ne.0) call print_mo(Ca,NBasis,'MONOMER A')
! if(SAPT%IPrint.ne.0) call print_mo(Cb,NBasis,'MONOMER B')

! look-up tables
 call select_active(SAPT%monA,NBasis,Flags%ICASSCF,Flags%ISHF,Flags%IFlCore)
 call select_active(SAPT%monB,NBasis,Flags%ICASSCF,Flags%ISHF,Flags%IFlCore)

! ABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB
! read and transform 2-el integrals
 call readtwoint(NBasis,'AOTWOINT_A')
! full 4-idx tran
! call tran4_full(NBasis,Ca,Cb,'TWOMOAB')
 ! integrals stored as (ov|ov)
 call tran4_gen(NBasis,&
                SAPT%monA%num0+SAPT%monA%num1,Ca,&
                SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
                SAPT%monB%num0+SAPT%monB%num1,Cb,&
                SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
                'TWOMOAB')

 if(SAPT%IPrint.gt.100) call print_TwoInt(NBasis)

 call print_active(SAPT,NBasis)

! calculate response
! if(SAPT%SaptLevel.gt.1) then
    call SaptInter(NBasis,SAPT%monA,Flags%ICASSCF)
!   call COMMTST(NBasis) 
    call calc_response(SAPT%monA,Ca,Flags,NBasis,'TWOMOAA',SAPT%EnChck)
! mon B
    call SaptInter(NBasis,SAPT%monB,Flags%ICASSCF)
!   call COMMTST(NBasis) 
    call calc_response(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB',SAPT%EnChck)
! endif

 deallocate(work1,work2)

! maybe better: add writing Ca, Cb to file?!!
 allocate(SAPT%monA%CMO(NBasis,NBasis),SAPT%monB%CMO(NBasis,NBasis)) 
 ij=0
 SAPT%monA%CMO = 0
 SAPT%monB%CMO = 0
 do i=1,NBasis
    do j=1,NBasis
       ij = ij + 1
       SAPT%monA%CMO(j,i) = Ca(ij)
       SAPT%monB%CMO(j,i) = Cb(ij)
    enddo 
 enddo 

 ! calculate electrostatic potential
 call calc_elpot(SAPT%monA,SAPT%monB,NBasis)

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ha,Hb,Va,Vb,S)
 deallocate(Ca,Cb)

end subroutine sapt_interface

subroutine calc_response(Mon,MO,Flags,NBas,fname,EChck)
implicit none

!type(SaptData) :: SAPT
type(SystemBlock) :: Mon
type(FlagsData) :: Flags
double precision :: MO(:)        
integer :: NBas
character(*) :: fname
logical :: EChck
integer :: NSq,NInte1,NInte2
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:) 
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:), &
                                 CMAT(:),EMAT(:),EMATM(:),&
                                 DMAT(:),DMATK(:),&
                                 EigVecR(:), Eig(:)
double precision,parameter :: One = 1d0, Half = 0.5d0
integer :: i,j,ij,ij1,ione,itwo 
double precision :: ACAlpha
double precision :: ECASSCF,ETot,ECorr
character(8) :: label
character(:),allocatable :: onefile,twofile,propfile,rdmfile
double precision :: tmp
double precision,parameter :: SmallE=0d0,BigE=1.D20
     
! set filenames
 if(Mon%Monomer==1) then
    onefile = 'ONEEL_A'
    twofile = 'TWOMOAA'
    propfile = 'PROP_A'
    rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    onefile = 'ONEEL_B'
    twofile = 'TWOMOBB'
    propfile = 'PROP_B'
    rdmfile='rdm2_B.dat'
 endif

! set dimensions
 NSq = NBas**2
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

 allocate(work1(NSq),work2(NSq),XOne(NInte1),URe(NBas,NBas),&
          TwoMO(NInte2))
 
 URe = 0d0
 do i=1,NBas
    URe(i,i) = 1d0
 enddo

! read 1-el
 open(newunit=ione,file=onefile,access='sequential',&
      form='unformatted',status='old')

 read(ione) 
 read(ione)
 read(ione) label, work1
 if(label=='ONEHAMIL') then
    call tran_oneint(work1,MO,MO,Mon%NSym,Mon%GFunc,work2,NBas)
    call sq_to_triang(work1,XOne,NBas) 
 else
    write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
    stop
 endif

! read 2-el
 call tran4_full(NBas,MO,MO,fname)
 call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 ACAlpha=One
 if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

   allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

   call ACABMAT0(ABPlus,ABMin,URe,Mon%Occ,XOne,TwoMO, &
                 NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,ACAlpha,1)
  
   !   do i=1,Mon%NDim
   !       write(LOUT,*) ABPlus(1,i),ABMin(1,i)
   !!      write(LOUT,*) ABPlus(i,i), ABMin(i,i)
   !   enddo

   ! reduce dim
   call reduce_dim('AB',ABPlus,ABMin,Mon)

    tmp=0
    do i=1,Mon%NDimX**2
        tmp = tmp + ABPlus(i)**2
   !    write(*,*) ABPlus(Mon%NDimX*(i-1)+i), ABMin(Mon%NDimX*(i-1)+i)
    enddo
    write(*,*) 'AB+',tmp


   EigVecR = 0
   Eig = 0
   call ERPASYMM(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)

   ! do i=1,Mon%NDimX**2
   !     !tmp = tmp + Eig(i)**2
   !     tmp = tmp + EigVecR(i)**2
   ! enddo
   ! write(*,*) tmp
   ! do i=1,Mon%NDimX
   !    write(LOUT,*) Eig(i),EigVecR(Mon%NDimX*(i-1)+i)
   ! enddo

   if(EChck) then
      write(LOUT,'(/,1x,a)') 'ERPA-GVB ENERGY CHECK REQUESTED:'
      call EneERPA(ETot,ECorr,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,&
           Mon%Occ,XOne,Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
   endif

 ! GVB with TD-APSG response
 elseif(Flags%ICASSCF==0.and.Flags%ISERPA==2) then

   Mon%NDimN=0
   do i=1,NBas
      if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
   enddo

   allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
            CMAT(Mon%NDim**2),EMAT(NBas**2),EMATM(NBas**2),&
            DMAT(Mon%NDim*NBas),DMATK(Mon%NDim*NBas),&
            EigVecR(2*(Mon%NDimX+Mon%NDimN)*2*(Mon%NDimX+Mon%NDimN)),&
            Eig(2*(Mon%NDimX+Mon%NDimN))) 

   CMAT=0
   call APSG_NEST(ABPlus,ABMin,CMAT,EMAT,EMATM,DMAT,DMATK,&
                  URe,Mon%Occ,XOne,TwoMO,&
                  NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,Flags%ISERPA)

   !reduce dimensions
   call reduce_dim('AB',ABPlus,ABMin,Mon)
   call reduce_dim('D',DMAT,DMATK,Mon)
   call reduce_dim('E',EMAT,EMATM,Mon)

  ! tmp=0
  ! do i=1,Mon%NDimX**2
  !     tmp = tmp + ABPlus(i)**2
  ! enddo
  ! write(*,*) 'AB-test:',tmp
  ! tmp=0
  ! do i=1,Mon%NDimN**2
  !     tmp = tmp + EMAT(i)**2
  ! enddo
  ! write(*,*) 'E+',tmp
  ! tmp=0
  ! do i=1,Mon%NDimN**2
  !     tmp = tmp + EMATM(i)**2
  ! enddo
  ! write(*,*) 'E-',tmp
  ! tmp=0d0
  ! do i=1,Mon%NDimN*Mon%NDimX
  !     tmp = tmp + DMAT(i)**2
  ! enddo
  ! write(*,*) 'D+',tmp
  ! tmp=0d0
  ! do i=1,Mon%NDimN*Mon%NDimX
  !     tmp = tmp + DMATK(i)**2
  ! enddo
  ! write(*,*) 'D-',tmp

!    write(*,*) Mon%NGem,Mon%NDimN,Mon%NDimX
!    write(*,*) Mon%Occ,NBas
    EigVecR = 0
    Eig = 0
    call PINOVEC(EigVecR,Eig,ABPlus,ABMin,DMAT,DMATK,EMAT,EMATM, &
                 Mon%Occ,NBas,Mon%NDimX,Mon%NDimN)

    print*, 'TESTY:',Mon%NDimX,Mon%NDimN,Mon%NDimX+Mon%NDimN
    do i=1,size(Eig)
       write(LOUT,*) i,Eig(i)
    enddo

   ! j=0
   ! do i=1,size(Eig)
   !    if(Eig(i).Gt.SmallE.And.Eig(i).Lt.BigE) then
   !    j = j + 1
   !    write(LOUT,*) j,i,Eig(i)
   !    endif
   ! enddo

!   call select_resp(EigVecR,Eig,Mon,NBas,&
!                    Mon%NDimX+Mon%NDimN)

   deallocate(CMAT,EMAT,EMATM,DMAT,DMATK)

 ! CAS-SCF
 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then

   allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

   call execute_command_line('cp '//rdmfile// ' rdm2.dat')
 
!   Gamma_2_AB CAN BE USED 
!   WITH DIFFERENT select_active
!   AND reducing dimensions!
!   call Gamma2_AB(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
!                  NBas,Mon%NDim,NInte1,NInte2,ACAlpha)
   ECASSCF = 0
   call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
               Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)

 !  tmp=0d0
 !  do i=1,Mon%NDimX**2
 !      tmp = tmp + ABMin(i)**2
 !     tmp = tmp + ABPlus(i)**2
 !    write(*,*) ABPlus(Mon%NDimX*(i-1)+i), ABMin(Mon%NDimX*(i-1)+i)
 !  enddo
 !  write(*,*) tmp

!  do i=1,Mon%NDimX
!     write(*,*) ABPlus(Mon%NDimX*(i-1)+i), ABMin(Mon%NDimX*(i-1)+i)
!  enddo

   EigVecR = 0
   Eig = 0
 
   call ERPASYMM1(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)

!   tmp=0
!   do i=1,Mon%NDimX
!      !write(*,*) Eig(i)
!       tmp = tmp + Eig(i)**2
!   enddo
!   write(*,*) tmp
!   tmp=0
!   do i=1,Mon%NDimX**2
!       tmp = tmp + EigVecR(i)**2
!   enddo
!   write(*,*) tmp

   if(EChck) then
      call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
                     Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
      ECorr=Ecorr*0.5d0
  
      write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
            ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr
   endif

 ! CAS-SCF Two-electron systems
 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==2.and.Mon%NELE==1) then
 
   Mon%NDimN=0
   do i=1,NBas
      if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
   enddo

   allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
            CMAT(Mon%NDim**2),EMAT(NBas**2),EMATM(NBas**2),&
            DMAT(Mon%NDim*NBas),DMATK(Mon%NDim*NBas),&
            EigVecR(2*(Mon%NDimX+Mon%NDimN)*2*(Mon%NDimX+Mon%NDimN)),&
            Eig(2*(Mon%NDimX+Mon%NDimN))) 
!            EigVecR(2*(Mon%NDim+NBas)*2*(Mon%NDim+NBas)),&
!            Eig(2*(Mon%NDim+NBas))) 

   CMAT=0
   call APSG_NEST(ABPlus,ABMin,CMAT,EMAT,EMATM,DMAT,DMATK,&
                  URe,Mon%Occ,XOne,TwoMO,&
                  NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,Flags%ISERPA)

   !reduce dimensions
   call reduce_dim('AB',ABPlus,ABMin,Mon)
   call reduce_dim('D',DMAT,DMATK,Mon)
   call reduce_dim('E',EMAT,EMATM,Mon)

  ! tmp=0
  ! do i=1,Mon%NDimX**2
  !      tmp = tmp + ABPlus(i)**2
  !  enddo
  !  write(*,*) 'AB-test:',tmp
  !  tmp=0
  !  do i=1,Mon%NDimN**2
  !      tmp = tmp + EMAT(i)**2
  !  enddo
  !  write(*,*) 'E+',tmp
  !  tmp=0
  !  do i=1,Mon%NDimN**2
  !      tmp = tmp + EMATM(i)**2
  !  enddo
  !  write(*,*) 'E-',tmp
  !  tmp=0d0
  !  do i=1,Mon%NDimN*Mon%NDimX
  !      tmp = tmp + DMAT(i)**2
  !  enddo
  !  write(*,*) 'D+',tmp
  !  tmp=0d0
  !  do i=1,Mon%NDimN*Mon%NDimX
  !      tmp = tmp + DMATK(i)**2
  !  enddo
  !  write(*,*) 'D-',tmp
  ! write(*,*) Mon%NGem,Mon%NDimN,Mon%NDimX
  ! write(*,*) Mon%Occ,NBas

   EigVecR = 0
   Eig = 0
   call PINOVEC(EigVecR,Eig,ABPlus,ABMin,DMAT,DMATK,EMAT,EMATM, &
                Mon%Occ,NBas,Mon%NDimX,Mon%NDimN)

   print*, 'TESTY:',Mon%NDimX,Mon%NDimN,Mon%NDimX+Mon%NDimN
   do i=1,size(Eig)
      write(LOUT,*) i,Eig(i)
   enddo

 !   j=0
 !   do i=1,size(Eig)
 !      if(Eig(i).Gt.SmallE.And.Eig(i).Lt.BigE) then
 !      j = j + 1
 !      write(LOUT,*) j,i,Eig(i)
 !      endif
 !   enddo


 endif

! dump response
 call writeresp(EigVecR,Eig,propfile)

 close(ione)
 deallocate(work1,work2,XOne,TwoMO,URe)
 deallocate(ABPlus,ABMin,EigVecR,Eig)

end subroutine calc_response 

subroutine calc_elpot(A,B,NBas)
implicit none

integer :: NBas
type(SystemBlock) :: A, B
double precision,allocatable :: Pa(:,:),Pb(:,:) 
double precision,allocatable :: Va(:,:),Vb(:,:) 
double precision,allocatable :: Ja(:,:),Jb(:,:) 

 allocate(Pa(NBas,NBas),Pb(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Ja(NBas,NBas),Jb(NBas,NBas))

 call get_den(NBas,A%CMO,2d0*A%Occ,Pa)
 call get_den(NBas,B%CMO,2d0*B%Occ,Pb)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 call make_J2(NBas,Pa,Pb,Ja,Jb)

 allocate(A%WPot(NBas,NBas),B%WPot(NBas,NBas))

 A%WPot = Va + Ja
 B%WPot = Vb + Jb

 deallocate(Jb,Ja,Vb,Va,Pb,Pa)

end subroutine calc_elpot 

subroutine reduce_dim(var,matP,matM,Mon)
implicit none

character(*) :: var
double precision :: matP(:),matM(:)
type(SystemBlock) :: Mon
integer :: i,j,ij,ij1

select case(var) 
case('AB','ab')
   ! matP = ABPlus
   ! matM = ABMin
   ij=0
   ij1=0
   do j = 1,Mon%NDimX
      do i = 1,Mon%NDimX
          ij  = (j-1)*Mon%NDimX + i
          ij1 = (Mon%IndX(j)-1)*Mon%NDim + Mon%IndX(i)
          matP(ij) = matP(ij1)
          matM(ij) = matM(ij1)
      enddo
   enddo

case('D','d')
   ! matP = DMAT
   ! matM = DMATK
   ij=0
   ij1=0
   do j = 1,Mon%NDimN
      do i = 1,Mon%NDimX
         ij  = (j-1)*Mon%NDimX + i
         ij1 = (j-1)*Mon%NDim + Mon%IndX(i)
         matP(ij) = matP(ij1)
         matM(ij) = matM(ij1)
      enddo
   enddo

case('E','e')
   ! matP = EMAT
   ! matM = EMATM
   ij=0
   ij1=0 
   do j=1,Mon%NDimN
      do i=1,Mon%NDimN
         ij  = (j-1)*Mon%NDimN + i
         ij1 = (j-1)*Mon%NBasis + i
         matP(ij) = matP(ij1)
         matM(ij) = matM(ij1)
      enddo
   enddo

end select

end subroutine reduce_dim

subroutine select_resp(EigVec,Eig,Mon,NBas,NDimEx)
implicit none
! is this even ok?
type(SystemBlock) :: Mon
integer :: NBas,NDimEx
double precision :: EigVec(2*NDimEx,2*NDimEx),Eig(2*NDimEx)
double precision,allocatable :: TmpVec(:,:),TmpEig(:)
integer :: i,j
integer :: pq,ip,iq

 !NDimEx = Mon%NDimX+Mon%NDimN
 allocate(TmpEig(NDimEx),TmpVec(NDimEx,NDimEx))

 TmpEig = 0
 TmpVec = 0
 j = 0
 do i=1,2*NDimEx
    if(Eig(i).gt.1d0.and.Eig(i).lt.1d20) then
       j = j + 1
       TmpEig(j) = Eig(i)
       TmpVec(1:NDimEx,j) = EigVec(1:NDimEx,i)
    endif
 enddo  
 
 print*, 'TmpEig'
 do i=1,NDimEx
    write(LOUT,*) i,TmpEig(i)
 enddo

! not too much... 
! this is so complicated...
! maybe a way to simplify it?
  print*, '1st'
  do i=1,NDimEx
    do pq=1,2*NDimEx
       if(pq.le.Mon%NDimX) then
          ip=Mon%IndN(1,pq)
          iq=Mon%IndN(2,pq)
       else
          ip=pq-Mon%NDimX
          iq=ip
       endif
       if(ip.gt.iq) then
          if(Eig(pq).gt.1d0.and.Eig(pq).lt.1d20) then
            write(6,*) ip,iq,EigVec(i,pq)
          endif
       endif
    enddo
  enddo 
 
  print*, '2nd'
! map(pq) here?
  do i=1,NDimEx
    do pq=1,NDimEx
       if(pq.le.Mon%NDimX) then
          ip=Mon%IndN(1,pq)
          iq=Mon%IndN(2,pq)
       else
          ip=pq-Mon%NDimX
          iq=ip
       endif
       if(ip.gt.iq) then
          write(6,*) ip,iq,TmpVec(i,pq)
       endif
    enddo
  enddo 



 deallocate(TmpEig,TmpVec)

end subroutine select_resp

function calc_vnn(A,B) result(Vnn)
implicit none

type(SystemBlock) :: A,B
integer :: ia, ib
double precision :: dx,dy,dz,dist
double precision :: Vnn

 Vnn=0d0
 do ia=1,A%NCen
    do ib=1,B%NCen
       dx = A%xyz(ia,1) - B%xyz(ib,1)
       dy = A%xyz(ia,2) - B%xyz(ib,2)
       dz = A%xyz(ia,3) - B%xyz(ib,3)
       dist = sqrt(dx**2+dy**2+dz**2)
       Vnn = Vnn + A%charg(ia)*B%charg(ib)/dist
    enddo
 enddo

end function calc_vnn

subroutine arrange_oneint(mat,nbas,SAPT)
implicit none

type(SaptData) :: SAPT
integer :: nbas
double precision :: mat(nbas,nbas)
integer :: nA, nB

call read_syminf(SAPT%monA,SAPT%monB,nbas)

call swap_rows(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)
call swap_cols(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)

end subroutine arrange_oneint

subroutine arrange_mo(mat,nbas,A,B)
implicit none

type(SystemBlock) :: A,B
integer :: nbas
double precision :: mat(nbas,nbas)

print*, A%NMonOrb,B%NMonOrb,size(mat)
call swap_rows(A%NMonOrb,B%NMonOrb,mat)

end subroutine arrange_mo

subroutine read_syminf(A,B,nbas)
! reads number of basis functions on each monomer
! from SYMINFO file!
implicit none

type(SystemBlock) :: A, B
integer :: nbas
integer :: iunit,ios
integer :: ibas,icen,last_ibas,last_icen
logical :: ex,dump
integer :: tmp
integer :: ACenTst, ACenBeg, ACenEnd

!print*, A%NCen, B%NCen

inquire(file='SYMINFO_B',EXIST=ex)

if(ex) then
   open(newunit=iunit,file='SYMINFO_B',status='OLD',&
        form='FORMATTED')
   read(iunit,*)
   read(iunit,*)
   dump = .TRUE.
   do
     read(iunit,'(i5,i6)',iostat=ios) ibas,icen
     if((icen.gt.B%NCen).and.dump) then 
        B%NMonOrb = ibas-1
        ACenBeg = icen
        dump = .FALSE.
     elseif(.not.dump.and.ios==0) then
        last_icen = icen
        last_ibas = ibas
     elseif(ios/=0) then
        tmp = last_ibas
        ACenEnd = last_icen
        exit
     endif
     !write(*,*) ibas,icen,dump
   enddo
   
   ACenTst = ACenEnd - ACenBeg + 1
   if(ACenTst/=A%NCen) then
      write(LOUT,'(1x,a)') 'ERROR! MISMATCH IN NUMBER OF ATOMS FOR MONONOMER A!'
      write(LOUT,'(1x,a,i3,4x,a,i3)') 'INPUT: ', A%NCen, 'SYMINFO_B: ', ACenTst
      stop
   endif
   A%NMonOrb = tmp - B%NMonOrb

   close(iunit)
else
   write(LOUT,'(1x,a)') 'ERROR! MISSING SYMINFO_B FILE!'
   stop
endif
! HERE
! maybe add test: tmp vs. NBasis? 

end subroutine read_syminf

subroutine swap_rows(nA,nB,mat)
implicit none

integer :: nA, nB
double precision :: mat(nA+nB,nA+nB)
double precision :: work(nA+nB,nA+nB)

! rows
work = 0d0
work(1:nA,:) = mat(nB+1:nB+nA,:)
work(nA+1:nA+nB,:) = mat(1:nB,:)

mat = work

end subroutine swap_rows

subroutine swap_cols(nA,nB,mat)
implicit none

integer :: nA, nB
double precision :: mat(nA+nB,nA+nB)
double precision :: work(nA+nB,nA+nB)

! columns
work = 0d0
work(:,1:nA) = mat(:,nB+1:nB+nA)
work(:,nA+1:nA+nB) = mat(:,1:nB)

mat = work

end subroutine swap_cols

subroutine read_mo(cmo,norb,nbas,nsiri,nmopun)
! in SAPT orbitals kept in AOMO order!
implicit none

integer :: iunit,norb,nbas
!double precision :: cmo(norb,nbas)
double precision :: cmo(nbas,nbas)
character(*) :: nsiri,nmopun
logical :: isiri
character(60) :: line
integer :: i,j
double precision :: natocc(10)
double precision :: tmp(nbas*norb)

inquire(file=nsiri,EXIST=isiri)

if(isiri) then

   open(newunit=iunit,file=nsiri,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

   call readlabel(iunit,'NEWORB  ')
   !read(iunit) cmo
   read(iunit) tmp

   cmo=0d0 
   do i=1,norb
      do j=1,nbas
         cmo(j,i) = tmp((i-1)*nbas + j)
      end do
   end do

   write(LOUT,'(1x,a)') 'Orbitals read from '//nsiri 

else
   print*, 'Achtung!!!',norb,nbas
   open(newunit=iunit,file=nmopun, &
        form='FORMATTED',status='OLD')
   read(iunit,'(a60)') line
   !do j=1,norb     
   do j=1,nbas     
      read(iunit,'(4f18.14)') (cmo(i,j),i=1,nbas)
   enddo
!   print*, line

   write(LOUT,'(1x,a)') 'Orbitals read from '//nmopun 
endif

close(iunit)

end subroutine read_mo

subroutine writeoneint(mon,ndim,S,V,H)
implicit none

integer :: ione,ndim
character(*) :: mon
double precision,dimension(ndim) :: S, V, H

 open(newunit=ione,file=mon,form='unformatted')
 write(ione) 'OVERLAP ', S 
 write(ione) 'POTENTAL', V
 write(ione) 'ONEHAMIL', H
 close(ione)

 write(LOUT,'(1x,a)') 'One-electron integrals written to file: '//mon
 write(LOUT,'()')

end subroutine writeoneint

subroutine writeresp(EVec,EVal,mon)
implicit none

character(*) :: mon
double precision :: EVec(:), EVal(:)
integer :: iunit

 open(newunit=iunit,file=mon,form='unformatted')
 write(iunit) EVec
 write(iunit) EVal
 close(iunit)

end subroutine writeresp

!subroutine readresp(EVec,EVal,NDim,fname)
!implicit none
!
!integer :: NDim
!double precision :: EVec(NDim,NDim), EVal(NDim)
!character(*) :: fname
!integer :: iunit
!
! open(newunit=iunit,file=fname,form='UNFORMATTED',&
!    access='SEQUENTIAL',status='OLD')
!
! read(iunit) EVec
! read(iunit) EVal
!
! close(iunit)
!
!end subroutine readresp

subroutine readgvb(mon,n,cfile)
implicit none

type(SystemBlock) :: mon
integer :: n
character(*) :: cfile
integer :: iunit
integer :: NAct, NIActive
integer :: i,j
!double precision,allocatable :: CICoef(:), Occ(:)
!integer,allocatable :: IGem(:)

open(newunit=iunit,file=cfile,form='FORMATTED',Status='OLD')
read(iunit,'(i5)') mon%NAct

mon%INAct = mon%NELE - mon%NAct

!write(*,*) mon%NELE, mon%NAct, mon%INAct

allocate(mon%CICoef(n),mon%IGem(n),mon%Occ(n))
mon%CICoef = 0d0

!!!HERE
do i=1,mon%INAct
   mon%CICoef(i) = 1.0d0
   mon%IGem(i) = i
enddo

read(iunit,*) (mon%CICoef(i+mon%INAct),i=1,2*mon%NAct)

do i=mon%INAct+1,mon%NELE
   mon%IGem(i) = i 
   mon%IGem(mon%NELE+i-mon%INAct) = i
enddo
mon%NGem = mon%NELE + 1

do i=1,n   
   if(mon%CICoef(i).eq.0d0) mon%IGem(i) = mon%NGem
   mon%Occ(i) = mon%CICoef(i)**2
enddo

close(iunit)

end subroutine readgvb

subroutine readmulti(nbas,mon,ihf,exsiri,isiri,occfile,occsir)
implicit none 

type(SystemBlock) :: mon
logical :: exsiri, ioccsir
integer :: isiri, nbas
character(*) :: occfile, occsir 
logical :: ihf
logical :: iocc
integer :: iunit,i
integer :: NAct, INAct
double precision :: Occ(nbas), sum1, sum2
double precision :: potnuc, emy, eactiv, emcscf
integer :: istate, ispin, nactel, lsym
integer :: nisht, nasht, nocct, norbt, nbast, nconf, nwopt, nwoph

 allocate(mon%CICoef(nbas),mon%IGem(nbas),mon%Occ(nbas))
! exsiri=.false.
 if(exsiri) then

    rewind(isiri) 
    read (isiri) 
    read (isiri) potnuc,emy,eactiv,emcscf, &
                 istate,ispin,nactel,lsym
    read (isiri) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
   
!    print*,    potnuc,emy,eactiv,emcscf, &
!               istate,ispin,nactel,lsym
!    write (*,*) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
   
    mon%NAct  = nasht
    mon%INAct = nisht

    if(nbast.ne.nbas) then
      write(LOUT,'(1x,a)') 'WARNING! NBasis FROM SIRIFC DOES NOT MATCH!'
      write(LOUT,'(1x,a,i5,1x,a,i5)') 'NBasis: ',nbas, 'SIRIFC: ', nbast 
      write(LOUT,'()')
      mon%IWarn = mon%IWarn + 1 
    endif

   if(.not.ihf) then
      ! CASCF case
      inquire(file=occsir,EXIST=ioccsir)
      if(ioccsir) then
         mon%Occ = 0d0
         open(newunit=iunit,file=occsir,status='OLD', &
              access='SEQUENTIAL',form='UNFORMATTED')
   
         call readlabel(iunit,'NATOCC  ')
         read(iunit) mon%Occ(1:mon%NAct+mon%INAct) 
   
         close(iunit)
      endif
   elseif(ihf) then
      ! Hartree-Fock case
      mon%Occ = 0d0
      mon%Occ(1:mon%NAct+mon%INAct) = 2d0

   endif

   sum1 = 0d0
   do i=1,mon%INAct+mon%NAct
       mon%Occ(i) = mon%Occ(i)/2d0
       sum1 = sum1 + mon%Occ(i)
   enddo
   mon%SumOcc = sum1

 endif
! occupations.dat 
 inquire(file=occfile,EXIST=iocc)
 if(iocc) then
 
    Occ = 0d0
    open(newunit=iunit,file=occfile,form='FORMATTED',status='OLD') 
 
    read(iunit,*) INAct, NAct
    INAct = INAct/2
    read(iunit,*) (Occ(i),i=1,INAct+NAct)
    sum2 = 0d0
    do i=1,INAct+NAct
       Occ(i) = Occ(i)/2d0
       sum2 = sum2 + Occ(i)
    enddo

    if(.not.ioccsir) then
       if(Abs(sum2-mon%XELE).gt.1.0d-8) then
          write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'  
          write(LOUT,'(1x,a,1x,f10.6,5x,a,i3)') 'SUM(OCC): ', sum2, 'MONOMER: ', mon%Monomer
          write(LOUT,'(1x,a)') 'CHECK occupations.dat!'  
          stop
       endif
       mon%INAct = INAct
       mon%NAct  = NAct
       mon%Occ   = Occ
       mon%SumOcc = sum2
       ! SIRIUS.RST not there  
       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occfile 
       write(LOUT,'()')

    else !compare SIRIUS.RST and occupations.dat 
       if(Abs(sum1-mon%XELE).gt.1.0d-8) then

          if(Abs(sum2-mon%XELE).gt.1.0d-8) then
             write(LOUT,'(1x,a,1x,i3)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE! MONOMER: ', mon%Monomer  
             write(LOUT,'(1x,a,1x,f10.6,5x,a,4x,f10.6)') 'OCC(SIRIUS): ', sum1,&
                          'OCC(occupations.dat)', sum2
             stop
          endif

          mon%INAct = INAct
          mon%NAct  = NAct
          mon%Occ   = Occ
          print*, 'sum1 zle' 
       endif

       ! both files correct     
       if(any(abs(Occ-mon%Occ).gt.1.d-9)) then
          write(LOUT,'(1x,a)') 'WARNING! DIFFERENT OCCUPANCIES IN SIRIUS.RST&
                & AND occupations.dat!'
          write(LOUT,'(1x,a)') 'OCCUPANCIES READ FROM occupations.dat!'
          write(LOUT,'()')
          mon%IWarn = mon%IWarn + 1
          mon%INAct = INAct
          mon%NAct  = NAct
          mon%Occ   = Occ
          mon%SumOcc = sum2
       endif

       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occsir 
       write(LOUT,'()')

    endif
 
    close(iunit)

 elseif(ioccsir) then
    if(Abs(sum1-mon%XELE).gt.1.0d-8) then
       write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'  
       write(LOUT,*) 'Occ: ', sum1
       write(LOUT,'(1x,a)') 'CHECK DALTON CALCULATIONS!'  
       stop
    endif 

 elseif(.not.ioccsir) then
     write(LOUT,'(1x,a)') 'ERROR! CANNOT READ OCCUPANCIES!'  
     stop
 endif

 if(mon%INAct==0) then
    mon%NGem = 2
    mon%IGem(1:mon%NAct+mon%INAct) = 1
    mon%IGem(mon%NAct+mon%INAct+1:nbas) = 2
 else
    mon%NGem = 3
    mon%IGem(1:mon%INAct) = 1
    mon%IGem(mon%INAct+1:mon%INAct+mon%NAct) = 2
    mon%IGem(mon%INAct+mon%NAct+1:nbas) = 3
 endif

! construct CICoef
 do i=1,nbas 
    mon%CICoef(i)=sqrt(Occ(I))
    if(mon%Occ(i).lt.0.5d0) mon%CICoef(i)=-mon%CICoef(i)
 enddo

! check
!  write(LOUT,*) mon%Occ

end subroutine readmulti

subroutine select_active(mon,nbas,ICASSCF,ISHF,IFlCore)
implicit none

type(SystemBlock) :: mon
integer :: nbas, ICASSCF, ISHF, IFlCore
integer :: i, j, ij, icnt
integer :: ind, ind_ij
integer :: IAuxGem(nbas) 
integer :: test

 IAuxGem = mon%IGem
 allocate(mon%IndAux(nbas))

 do i=1,mon%NELE
    mon%IndAux(i)=0
 enddo
 do i=1+mon%NELE,nbas
    mon%IndAux(i)=2
 enddo

! active orbitals
 mon%icnt = 0
 if(ICASSCF==0) then
   do i=1,mon%NELE
      if(mon%Occ(i).lt.mon%ThrAct) then
         mon%IndAux(i)=1
         !write(6,'(/,X," Active Orbital: ",I4,E14.4)') &
         !      i, mon%Occ(i)
         mon%IndAux(FindGem(i,mon))=1
         !write(6,'(X," Active Orbital: ",I4,E14.4)') &
         !FindGem(i,mon), mon%Occ(FindGem(i,mon))
         mon%icnt = mon%icnt + 2
      endif
   enddo
 elseif(ICASSCF==1.and.ISHF==0) then
   write(LOUT,'()')
   if(mon%Monomer==1) write(LOUT,'(1x,a)') 'Monomer A' 
   if(mon%Monomer==2) write(LOUT,'(1x,a)') 'Monomer B' 
   do i=1,nbas
      if(mon%Occ(i).lt.1d0.and.mon%Occ(i).ne.0d0) then
         mon%IndAux(i) = 1
         write(6,'(X," Active Orbital: ",I4,E14.4)') i, mon%Occ(i)
         mon%icnt = mon%icnt + 1
     endif
   enddo
 endif

! set generalized "occupied" = num0 + num1
! and "virtual" = num1 + num2 indices
 mon%num0 = 0 
 do i=1,nbas
    if(mon%IndAux(i)/=0) exit
    mon%num0 = mon%num0 + 1
 enddo
 mon%num2 = 0
 do i=nbas,1,-1
    if(mon%IndAux(i)/=2) exit
    mon%num2 = mon%num2 + 1
 enddo
 mon%num1 = nbas - mon%num0 - mon%num2

! print*, '***TESTY***'
! print*, mon%num0,mon%num1,mon%num2
! do i=1,nbas
!    write(LOUT,*) i,mon%IndAux(i)
! enddo

! active pairs
 allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))

 mon%IPair(1:nbas,1:nbas) = 0

 if(ICASSCF==0) then

    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1
   
          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
      !   do not correlate active degenerate orbitals from different geminals 
             if((mon%IGem(i).ne.mon%IGem(j)).and.&
                (mon%IndAux(i)==1).and.(mon%IndAux(j)==1).and.&
                (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-2)) then
      
                write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
              else
                 ! if IFlCore=0 exclude core (inactive) orbitals
                 if(IFlCore==1.or.&
                   (IFlCore==0.and.&
                    mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
      
                    ind = ind + 1
                    mon%IndX(ind) = ij
                    mon%IndN(1,ind) = i
                    mon%IndN(2,ind) = j
                    mon%IPair(i,j) = 1
                    mon%IPair(j,i) = 1
      
                 endif
              endif
      
          endif
   
       enddo
    enddo

 elseif(ICASSCF==1) then

    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1
   
          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
      !   do not correlate active degenerate orbitals from different geminals 
               if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1).and.&
                (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
      
                write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
              else
                 ! if IFlCore=0 exclude core (inactive) orbitals
                 if(IFlCore==1.or.&
                   (IFlCore==0.and.&
                    mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
      
                    ind = ind + 1
                    mon%IndX(ind) = ind
                    mon%IndN(1,ind) = i
                    mon%IndN(2,ind) = j
                    mon%IPair(i,j) = 1
                    mon%IPair(j,i) = 1
      
                 endif
              endif
      
          endif
   
       enddo
    enddo
 
 endif
 mon%NDimX = ind

! old:
! ij=0
! ind = 0
! do i=1,nbas
!    do j=1,i-1
!
!       ij = ij + 1
!       ind_ij = mon%IndAux(i)+mon%IndAux(j)
!       if((ind_ij/=0).and.(ind_ij/=4)) then
!   !   do not correlate active degenerate orbitals from different geminals 
!          if(ICASSCF==0.and.(mon%IGem(i).ne.mon%IGem(j)).and.&
!             (mon%IndAux(i)==1).and.(mon%IndAux(j)==1).and.&
!             (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-2)&
!          .or.&
!             (ICASSCF==1).and.& 
!             (mon%IndAux(i)==1).and.(mon%IndAux(j)==1).and.&
!             (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
!   
!             write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
!           else
!              ! if IFlCore=0 exclude core (inactive) orbitals
!              if(IFlCore==1.or.&
!                (IFlCore==0.and.&
!                 mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
!   
!                 ind = ind + 1
!                 mon%IndX(ind) = ind
!                 !mon%IndX(ind) = ij
!                 mon%IndN(1,ind) = i
!                 mon%IndN(2,ind) = j
!                 mon%IPair(i,j) = 1
!                 mon%IPair(j,i) = 1
!   
!              endif
!           endif
!   
!       endif
!
!    enddo
! enddo
!
! mon%NDimX = ind
! Write(6,'(/,2X,"Total number of pairs:",I6)') mon%NDim !nbas*(nbas-1)/2 
! Write(6,'(2X,"Reduced to:",I6)') mon%NDimX 

contains 

function FindGem(io,mon) result(IFindG)
implicit none

type(SystemBlock) :: mon
integer :: IFindG,i,io

 IFindG = 0
 do i=1,2*mon%NELE
    if((mon%IGem(io).eq.mon%IGem(i)).and.(io.ne.i))  IFindG=i 
 enddo
 
 if(IFindG==0) IFindG = io

end function FindGem

end subroutine select_active


subroutine readener(nbas,mon,isiri)
implicit none

type(SystemBlock) :: mon
integer :: nbas, isiri
integer :: MMORBT
integer :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
           NCDETS, NCMOT,NNASHX,NNASHY,NNORBT 
integer :: i
double precision,allocatable :: fock(:)

! set dimensions
 rewind(isiri)

 read (isiri)
 read (isiri) 
 read (isiri) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
              NCDETS, NCMOT,NNASHX,NNASHY,NNORBT 
 read(isiri)
 read(isiri)
 read(isiri)
 read(isiri)
 read(isiri)

 MMORBT = max(4,NNORBT)
 allocate(fock(MMORBT),mon%OrbE(NORBT))

 read(isiri) fock 
 ! orb energies: diag of Fock
 do i=1,NORBT
    mon%OrbE(i) = fock(i+i*(i-1)/2)
 enddo

 deallocate(fock)

end subroutine readener 

subroutine print_occ(nbas,SAPT,ICASSCF)
implicit none
!!! HERE : Change to A/B monomers!
type(SaptData) :: SAPT
integer :: nbas, ICASSCF
integer :: i

 associate(A => SAPT%monA, B => SAPT%monB)
 if(ICASSCF==0) then
   write(LOUT,'(1x,a)') 'ORBITAL OCCUPANCIES'
   write(LOUT,'(2x,"Orb",3x,"Occupancy-A",6x,"Gem-A",10x,"Occupancy-B",6x,"Gem-B")')
   do i=1,nbas
      !write(6,'(X,i3,2e16.6,i6)') i,A%Occ(i),A%CICoef(i),A%IGem(i)
      write(6,'(x,i3,e16.6,i6,10x,e16.6,i6)') i,A%Occ(i),A%IGem(i),B%Occ(i),B%IGem(i)
      !write(6,'(X,i3,e16.6,i6)') i,B%Occ(i),B%IGem(i)
   enddo
   write(LOUT,'()') 

 else

!   write(LOUT,'(1x,a)') 'NO OF CAS INACTIVE AND ACTIVE ORBITALS'
   write(LOUT, '()')
   write(LOUT,'(1x,a,11x,a,5x,a)') 'CAS ORBITALS','Monomer A',  'Monomer B'
   write(LOUT,'(1x,a,17x,i3,11x,i3)') 'INACTIVE', A%INAct, B%INAct
   write(LOUT,'(1x,a,17x,i3,11x,i3)') 'ACTIVE  ', A%NAct, B%NAct
   write(LOUT, '()')
!   write(LOUT,'(1x,a,1x,i3,i3)') 'NO OF CAS INACTIVE AND ACTIVE ORBITALS:',A%INAct, A%NAct
!   write(LOUT,'(1x,a,1x,i3,i3)') 'NO OF CAS INACTIVE AND ACTIVE ORBITALS:',B%INAct, B%NAct
   write(LOUT,'(1x,a)') 'ORBITAL OCCUPANCIES'
   write(LOUT,'(1x,a,3x,a,4x,a,10x,a,6x,a)') 'CASSCF', 'Occupancy-A', 'Gem-A', 'Occupancy-B','Gem-B'
   do i=1,nbas
      write(LOUT,'(1x,i3,1x,e16.6,1x,i6,7x,e16.6,3x,i6)') i, A%Occ(i),A%IGem(i),B%Occ(i),B%IGem(i) 
   enddo
   write(LOUT,'(2x,a,f8.4,18x,f6.4)') 'SUM OF OCCUPANCIES: ', A%SumOcc, B%SumOcc
   write(LOUT, '()')
 endif
 end associate

end subroutine print_occ

subroutine print_active(SAPT, nbas)
implicit none

type(SaptData) :: SAPT
integer :: nbas, i, ip 
integer :: NDimX

! print orbs  
 write(LOUT,'()')
 write(LOUT,'(27x,a,4x,a)') 'Monomer A', 'Monomer B' 
 do i=1,nbas 
   associate(IndA => SAPT%monA%IndAux(i), &
             OccA => SAPT%monA%Occ(i), &
             IndB => SAPT%monB%IndAux(i), & 
             OccB => SAPT%monB%Occ(i) )
     if(IndA==1.or.IndB==1) then
        write(LOUT,'(1x,a,2x,i2)',advance='no') 'Active orbital: ', i
        if(IndA==1) then
           write(LOUT,'(e14.4)',advance='no') OccA
        else
           write(LOUT,'(14x)',advance='no')
        endif 
        if(IndB==1) then
           write(LOUT,'(e14.4)') OccB
        else
           write(LOUT, '()') 
        endif
     endif
   end associate
 enddo
 write(LOUT,'(1x,6a)') ('--------',i=1,6)
 write(LOUT,'(1x,a,14x,i3,9x,i3)') 'Total Active: ', SAPT%monA%icnt, SAPT%monB%icnt

! print pairs
 if(SAPT%IPrint.gt.0) then
    !NDim = nbas*(nbas-1)/2 
    write(LOUT,'()')
    write(LOUT,'(26x,a,5x,a)') 'Monomer A', 'Monomer B' 
    write(LOUT,'(1x,a,2x,i6,8x,i6)') 'Total number of pairs: ', SAPT%monA%NDim,SAPT%monB%NDim
    write(LOUT,'(1x,a,12x,i6,8x,i6)') 'Reduced to: ', SAPT%monA%NDimX, SAPT%monB%NDimX 
   
    NDimX = max(SAPT%monA%NDimX,SAPT%monB%NDimX) 
    write(LOUT,'()')
    write(LOUT,'(2x,"Accepted pairs:")')
    write(LOUT,'(2x,a,11x,a,28x,a)') 'p  q', 'Monomer A', 'MonomerB'
    write(LOUT,'(2x,8a)',advance='no') ('----',i=1,8)
    write(LOUT,'(5x,8a)') ('----',i=1,8)
    do ip=1,NDimX
       if(ip.gt.SAPT%monA%NDimX) then
          write(LOUT,'(34x)',advance='no')
       else
          associate( idx1 => SAPT%monA%IndN(1,ip), &
                     idx2 => SAPT%monA%IndN(2,ip), &
                     Occ => SAPT%monA%Occ )
            write(LOUT,'(2i3,2e14.4)',advance='no') &
                         idx1,idx2,Occ(idx1),Occ(idx2)
          end associate
       endif
   
       if(ip.gt.SAPT%monB%NDimX) then
          write(LOUT,'(14x)')
       else
          associate( idx1 => SAPT%monB%IndN(1,ip), &
                     idx2 => SAPT%monB%IndN(2,ip), &
                     Occ => SAPT%monB%Occ )
            write(LOUT,'(3x,2i3,2e14.4)') idx1,idx2,Occ(idx1),Occ(idx2)
          end associate
       endif
      
    enddo
 endif

end subroutine print_active

subroutine print_warn(SAPT)
implicit none

type(SaptData) :: SAPT
integer :: cnt,i

cnt = SAPT%monA%IWarn+SAPT%monB%IWarn 
if(cnt.gt.0) then
    write(LOUT,'(1x,a,i2,1x,a)') 'SAPT: CHECK OUTPUT FOR',cnt,'WARNINGS!'
    write(LOUT,'(8a10)') ('**********',i=1,8)
endif

end subroutine print_warn

subroutine free_sapt(SAPT)
implicit none

type(SaptData) :: SAPT

deallocate(SAPT%monA%CICoef,SAPT%monA%IGem,SAPT%monA%Occ, &
           SAPT%monA%IndAux,SAPT%monA%IndX,SAPT%monA%IndN,&
           SAPT%monA%CMO,&
           SAPT%monA%IPair)
deallocate(SAPT%monB%CICoef,SAPT%monB%IGem,SAPT%monB%Occ, &
           SAPT%monB%IndAux,SAPT%monB%IndX,SAPT%monB%IndN,&
           SAPT%monB%CMO,&
           SAPT%monB%IPair)

! HERE - change to SAPTLEVEL?
if(allocated(SAPT%monA%WPot)) then
  deallocate(SAPT%monA%WPot)
endif
if(allocated(SAPT%monB%WPot)) then
  deallocate(SAPT%monB%WPot)
endif

if(allocated(SAPT%monA%OrbE)) then
  deallocate(SAPT%monA%OrbE)
endif
if(allocated(SAPT%monB%OrbE)) then
  deallocate(SAPT%monB%OrbE)
endif

end subroutine free_sapt

!subroutine tranMO(C,nbas)
!implicit none
!
!integer :: nbas
!double precision :: C(nbas,nbas)
!double precision :: tmp(nbas,nbas)
!integer :: i,j
!
!do i=1,nbas
!do j=1,nbas
! tmp(j,i) = C(i,j)
!enddo
!enddo
!
!C = tmp
!
!end subroutine tranMO

end module sapt_main

