module sapt_main
use types
use systemdef
use timing
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
integer :: NBasis
double precision :: Tcpu,Twall

! TEMPORARY - JOBTYPE_2
! ERPA
 Flags%IFlAC  = 0
 Flags%IFlSnd = 0

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'STARTING SAPT CALCULATIONS'
 write(LOUT,'(8a10)') ('**********',i=1,8)


 call clock('START',Tcpu,Twall) 
 call sapt_interface(Flags,SAPT)

 ! SAPT components
 write(LOUT,'()')
 if(Flags%ISERPA==0) then

    call e1elst(SAPT%monA,SAPT%monB,SAPT)
    ! temporary here!!!!
    if(Flags%ICASSCF==1) then
    !   call e1exchs2(SAPT%monA,SAPT%monB,SAPT)
    endif
    if(SAPT%SaptLevel==0) then
       call e2disp_unc(Flags,SAPT%monA,SAPT%monB,SAPT)
   else
       call e2ind(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2disp(Flags,SAPT%monA,SAPT%monB,SAPT)
    endif

 elseif(Flags%ISERPA==2) then

    call e1elst(SAPT%monA,SAPT%monB,SAPT)
    !call e1exchs2(SAPT%monA,SAPT%monB,SAPT)
    call e2ind_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e2disp_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)

 endif

 call print_warn(SAPT)
 call free_sapt(SAPT)

 call clock('SAPT',Tcpu,Twall)

! call interface to PRDMFT

 stop

end subroutine sapt_driver

subroutine sapt_interface(Flags,SAPT) 
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT

integer :: NBasis,NSq,NInte1,NInte2
integer :: NCMOt, NOrbt, NBasist 
integer :: NSym, NBas(8)
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

 allocate(work1(NInte1),work2(NSq))
 allocate(Ha(NSq),Hb(NSq))
 allocate(Va(NSq),Vb(NSq),S(NSq))
! read and dump 1-electron integrals 
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
 open(newunit=ione,file='AOONEINT_A',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NBas(1:NSym),SAPT%monA%PotNuc
 
 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call square_oneint(work1,Ha,NBasis,NSym,NBas)
! call triang_to_sq(work1,Ha,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call square_oneint(work1,work2,NBasis,NSym,NBas)
! call triang_to_sq(work1,work2,NBasis)
 Va(:) = Ha - work2

 call readlabel(ione,'OVERLAP ')
 call readoneint(ione,work1)
 call square_oneint(work1,S,NBasis,NSym,NBas)
! call triang_to_sq(work1,S,NBasis)

 call readlabel(ione,'ISORDK  ')
 read(ione) 
 read(ione) SAPT%monA%charg,ncen,SAPT%monA%xyz 

! print*, 'MONO-A',ncen
! write(LOUT,*) SAPT%monA%charg(1:ncen)
! do i=1,ncen
!    write(LOUT,*) SAPT%monA%xyz(i,:)
! enddo

! write(*,*) 'VA'
! call print_sqmat(Va,NBasis)
! call print_diag(Va,NBasis)

 close(ione)
 
 SAPT%monA%NSym = NSym
 SAPT%monA%NSymBas(1:NSym) = NBas(1:NSym)

 ! square form
 call writeoneint('ONEEL_A',NSq,S,Va,Ha)
 
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

 open(newunit=ione,file='AOONEINT_B',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NBas(1:NSym),SAPT%monB%PotNuc

 !print*, '1B!',NSym,NBas,SAPT%monB%PotNuc

 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call square_oneint(work1,Hb,NBasis,NSym,NBas)
 !call triang_to_sq(work1,Hb,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call square_oneint(work1,work2,NBasis,NSym,NBas)
! call triang_to_sq(work1,work2,NBasis)
 Vb(:) = Hb - work2

 call readlabel(ione,'ISORDK  ')
 read(ione) 
 read(ione) SAPT%monB%charg,ncen,SAPT%monB%xyz 

! print*, 'MONO-B',ncen
! write(LOUT,*) SAPT%monB%charg(1:ncen)
! do i=1,ncen
!    write(LOUT,*) SAPT%monB%xyz(i,:)
! enddo

 SAPT%monB%NSym = NSym
 SAPT%monB%NSymBas(1:NSym) = NBas(1:NSym)

! rearrange in V: (B,A) -> (A,B)

 call read_syminf(SAPT%monA,SAPT%monB,NBasis)

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

    SAPT%monA%NOrb = NOrbt
    SAPT%monA%NSymOrb(1:NSym) = NOrbs(1:NSym)

 else
 !   SAPT%monA%NOrb = NOrb
 !   SAPT%monA%GFunc(1) = NOrb
    NBasist = NBasis
!    NCMOt = NOrb*NBasis
 endif 

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
 
 print*, Flags%ISHF, SAPT%monA%ISHF
 if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monA%NELE/=1.and.(.not.SAPT%monA%ISHF)) then
 !if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monA%NELE/=1) then

    ! CASSCF
    call readmulti(NBasis,SAPT%monA,.false.,exsiri,isiri,'occupations_A.dat','SIRIUS_A.RST')

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monA%NELE==1.and.Flags%ISERPA==0.and.(.not.SAPT%monA%ISHF)) then
    print*, 'here?-why so?'
    ! CASSCF
    ! for 2-el electron case: read from occupations.dat
    call readmulti(NBasis,SAPT%monA,.false.,.false.,isiri,'occupations_A.dat','SIRIUS_A.RST')

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monA%NELE==1.and.Flags%ISERPA==2) then

    call readmulti(NBasis,SAPT%monA,.false.,exsiri,isiri,'occupations_A.dat','SIRIUS_A.RST')

 elseif(Flags%ICASSCF==1.and.(Flags%ISHF==1.or.SAPT%monA%ISHF)) then

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
! if(exsiri) then
    open(newunit=isiri,file='SIRIFC_B',status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(isiri,'TRCCINT ')
    read(isiri) NSym,NOrbt,NBasist,NCMOt,NOcc(1:NSym),NOrbs(1:NSym)
 !   write(*,*) 'B: ',NSym,NOrbt,NBasist,NCMOt
 !   write(*,*) NOrb, 'NOrb'
    SAPT%monB%NOrb = NOrbt
    SAPT%monB%NSymOrb(1:NSym) = NOrbs(1:NSym)
! else
!    SAPT%monB%NOrb = NOrb
!    SAPT%monB%GFunc(1) = NOrb
    
!    NBasist = NBasis
!    NCMOt = NOrb*NBasis
! endif 

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
!    write(*,*)   potnuc,emy,eactiv,emcscf

 if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monB%NELE/=1) then

    ! CASSCF
    call readmulti(NBasis,SAPT%monB,.false.,exsiri,isiri,'occupations_B.dat','SIRIUS_B.RST')

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monB%NELE==1.and.Flags%ISERPA==0) then

    ! CASSCF
    ! for 2-el electron case: read from occupations.dat
    ! info in SIRIFC seems wrong in this case?
     call readmulti(NBasis,SAPT%monB,.false.,.false.,isiri,'occupations_B.dat','SIRIUS_B.RST')

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.SAPT%monB%NELE==1.and.Flags%ISERPA==2) then

    ! CASSCF - TD-APSG
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

 call read_mo(Ca,NBasis,SAPT%monA%NSym,SAPT%monA%NSymBas,SAPT%monA%NSymOrb,&
              'SIRIUS_A.RST','DALTON_A.MOPUN')
 call read_mo(Cb,NBasis,SAPT%monB%NSym,SAPT%monB%NSymBas,SAPT%monB%NSymOrb,&
              'SIRIUS_B.RST','DALTON_B.MOPUN')

 call arrange_mo(Cb,NBasis,SAPT)

! HERE NEW STUFF
 if(SAPT%monA%NSym.gt.1) then
    call sort_sym_mo(Ca,NBasis,SAPT%monA)
 endif
 if(SAPT%monB%NSym.gt.1) then
    call sort_sym_mo(Cb,NBasis,SAPT%monB)
 endif

! MAYBE: one should print with NOrbt?
 !if(SAPT%IPrint.ne.0) call print_mo(Ca,NBasis,'MONOMER A')
 !if(SAPT%IPrint.ne.0) call print_mo(Cb,NBasis,'MONOMER B')

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

! look-up tables
 call select_active(SAPT%monA,NBasis,Flags)
 call select_active(SAPT%monB,NBasis,Flags)

! ABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB
! read and transform 2-el integrals
 call readtwoint(NBasis,'AOTWOINT_A')
! full 4-idx tran
! call tran4_full(NBasis,Ca,Cb,'TWOMOAB')

 if(Flags%ISERPA==0) then

   ! integrals stored as (ov|ov)
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
                  'TWOMOAB')

   ! testing Ecorr ERPASYMM/ERPAVEC
   !call tran4_gen(NBasis,&
   !               SAPT%monA%num0+SAPT%monA%num1,Ca,&
   !               SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
   !               SAPT%monA%num0+SAPT%monA%num1,Ca,&
   !               SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
   !               'TMPMOAA')
   !call tran4_gen(NBasis,&
   !               SAPT%monB%num0+SAPT%monB%num1,Cb,&
   !               SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
   !               SAPT%monB%num0+SAPT%monB%num1,Cb,&
   !               SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
   !               'TMPMOBB')

   ! this is for testing E1exchS2...
   !call tran4_gen(NBasis,&
   !               NBasis,Ca,NBasis,Ca,&
   !               NBasis,Cb,NBasis,Cb,&
   !               'TMPMOAB')

   ! <oo|oo>
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  'TMPOOAB')

 elseif(Flags%ISERPA==2.and.Flags%ISHF==0) then
   ! integrals stored as (oFull,oFull)
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  NBasis,Ca,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  NBasis,Cb,&
                  'TWOMOAB')
   ! this is for testing E1exchS2...
   call tran4_gen(NBasis,&
                  NBasis,Ca,NBasis,Ca,&
                  NBasis,Cb,NBasis,Cb,&
                  'TMPMOAB')


 endif

 if(SAPT%IPrint.gt.100) call print_TwoInt(NBasis)

 call print_active(SAPT,NBasis)

! calculate response
! mon A
 call SaptInter(NBasis,SAPT%monA,Flags%ICASSCF)
! call COMMTST(NBasis)
 if(SAPT%SaptLevel.eq.0) then
     call calc_resp_unc(SAPT%monA,Ca,Flags,NBasis,'TWOMOAA')
 elseif(SAPT%SaptLevel.gt.0) then
   ! call calc_resp_unc(SAPT%monA,Ca,Flags,NBasis,'TWOMOAA')
    call calc_resp_full(SAPT%monA,Ca,Flags,NBasis,'TWOMOAA',SAPT%EnChck) 
 endif

! mon B
    call SaptInter(NBasis,SAPT%monB,Flags%ICASSCF)
!   call COMMTST(NBasis) 
 if(SAPT%SaptLevel.eq.0) then
    call calc_resp_unc(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB')
 elseif(SAPT%SaptLevel.gt.0) then
    !call calc_resp_unc(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB')
    call calc_resp_full(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB',SAPT%EnChck)
 endif

  if(Flags%ISERPA==2.and.Flags%ISHF==1) then

     call tran4_gen(NBasis,&
                      (SAPT%monA%num0+SAPT%monA%num1),SAPT%monA%CMO,&
                      NBasis,SAPT%monA%CMO,&
                      (SAPT%monB%num0+SAPT%monB%num1),SAPT%monB%CMO,&
                      NBasis,SAPT%monB%CMO,&
                      'TWOMOAB')
 
     !call tran4_gen(NBasis,&
     !             NBasis,Ca,NBasis,Ca,&
     !             NBasis,Cb,NBasis,Cb,&
     !             'TMPMOAB')

 endif


 deallocate(work1,work2)

 ! calculate electrostatic potential
 call calc_elpot(SAPT%monA,SAPT%monB,NBasis)

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ha,Hb,Va,Vb,S)
 deallocate(Ca,Cb)

end subroutine sapt_interface

subroutine calc_resp_unc(Mon,MO,Flags,NBas,fname)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags
double precision :: MO(:)        
integer :: NBas
character(*) :: fname
integer :: NSq,NInte1,NInte2
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:) 
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:), &
                                 EigY(:), EigY1(:), Eig(:), Eig1(:)
integer :: i,ione
double precision,parameter :: One = 1d0, Half = 0.5d0
character(8) :: label
character(:),allocatable :: onefile,twofile,propfile0,propfile1,rdmfile
double precision :: tmp

! perform check
 if(Flags%ICASSCF==0) then
   write(LOUT,*) 'ERROR! E2DISP UNC AVAILABLE ONLY FOR CAS/HF!'
   stop
 endif

! set filenames
 if(Mon%Monomer==1) then
    onefile = 'ONEEL_A'
    twofile = 'TWOMOAA'
    propfile0 = 'PROP_A0'
    propfile1 = 'PROP_A1'
    rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    onefile = 'ONEEL_B'
    twofile = 'TWOMOBB'
    propfile0 = 'PROP_B0'
    propfile1 = 'PROP_B1'
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
    call tran_oneint(work1,MO,MO,work2,NBas)
    call sq_to_triang(work1,XOne,NBas) 
 else
    write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
    stop
 endif

 ! transform and read 2-el integrals
 call tran4_full(NBas,MO,MO,fname)
 call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
          EigY(Mon%NDimX**2),EigY1(Mon%NDimX**2), &
          Eig(Mon%NDimX),Eig1(Mon%NDimX))

 ! read 2-RDMs
 call read2rdm(Mon,NBas)
 call system('cp '//rdmfile// ' rdm2.dat')

 EigY  = 0
 EigY1 = 0
 Eig   = 0
 Eig1  = 0

 call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
      EigY,EigY1,Eig,Eig1, &
      Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)

! dump response
! maybe construct smarter resp-files?

 call writeresp(EigY,Eig,propfile0)
 if(Flags%IFlag0==0) then
    call writeresp(EigY1,Eig1,propfile1)
 endif

 close(ione)
 deallocate(TwoMO,URe,XOne,work1,work2)
 deallocate(Eig1,Eig,EigY1,EigY,ABMin,ABPlus)

end subroutine calc_resp_unc

subroutine test_resp_unc(Mon,URe,XOne,TwoMO,NBas,NInte1,NInte2,IFlag0) 
implicit none

type(SystemBlock) :: Mon 
double precision,intent(in) :: TwoMO(NInte2),XOne(NInte1),URe(NBas,NBas)
integer,intent(in) :: NBas,NInte1,NInte2,IFlag0
character(:),allocatable :: propfile0,propfile1
double precision, allocatable :: ABPlus(:),ABMin(:)
double precision, allocatable :: EigY0(:),EigY1(:), &
                                 Eig0(:),Eig1(:)

! set filenames
 if(Mon%Monomer==1) then
    propfile0 = 'PROP_A0'
    propfile1 = 'PROP_A1'
    !rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    propfile0 = 'PROP_B0'
    propfile1 = 'PROP_B1'
    !rdmfile='rdm2_B.dat'
 endif


 allocate(EigY0(Mon%NDimX**2),EigY1(Mon%NDimX**2),&
          Eig0(Mon%NDimX),Eig1(Mon%NDimX), &
          ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2))
 EigY0 = 0
 EigY1 = 0
 Eig0 = 0
 Eig1 = 0

 call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
             EigY0,EigY1,Eig0,Eig1, &
             Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,IFlag0)

 ! dump uncoupled response
 print*, norm2(EigY0)
 call writeresp(EigY0,Eig0,propfile0)
 if(IFlag0==0) then
    call writeresp(EigY1,Eig1,propfile1)
 endif

 deallocate(Eig1,Eig0,EigY1,EigY0)

end subroutine test_resp_unc

subroutine calc_resp_full(Mon,MO,Flags,NBas,fname,EChck)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags
double precision :: MO(:)        
integer :: NBas
character(*) :: fname
logical :: EChck
integer :: NSq,NInte1,NInte2
double precision, allocatable :: work1(:),work2(:),XOne(:),TwoMO(:) 
double precision, allocatable :: ABPlus(:), ABMin(:), URe(:,:), &
                                 CMAT(:),EMAT(:),EMATM(:), &
                                 DMAT(:),DMATK(:), &
                                 EigVecR(:), Eig(:)
double precision, allocatable :: Eig0(:), Eig1(:), EigY0(:), EigY1(:) 
double precision :: Dens(NBas,NBas)
double precision,parameter :: One = 1d0, Half = 0.5d0
integer :: i,j,ij,ij1,ione,itwo 
integer :: ind
double precision :: ACAlpha
double precision :: ECASSCF,ETot,ECorr
character(8) :: label
character(:),allocatable :: onefile,twofile,propfile,rdmfile
character(:),allocatable :: propfile0,propfile1
double precision :: tmp
double precision,parameter :: SmallE=0d0,BigE=1.D20
double precision,external  :: trace
integer :: offset
integer,allocatable :: sort(:)
double precision, allocatable :: EigTmp(:), VecTmp(:)
! testy
 integer :: ii,jj,dimV1
 double precision,parameter :: thresh=1.d-10
 integer :: space1(2,Mon%NDimX)
 integer :: IndN_tmp(2,Mon%NDim)
 double precision :: work(Mon%NDimX)
 double precision,external :: ddot

! set filenames
 if(Mon%Monomer==1) then
    onefile  = 'ONEEL_A'
    twofile  = 'TWOMOAA'
    propfile = 'PROP_A'
    propfile0  = 'PROP_A0' 
    propfile1  = 'PROP_A1' 
    rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    onefile  = 'ONEEL_B'
    twofile  = 'TWOMOBB'
    propfile = 'PROP_B'
    propfile0  = 'PROP_B0' 
    propfile1  = 'PROP_B1' 
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
    call tran_oneint(work1,MO,MO,work2,NBas)
    call sq_to_triang(work1,XOne,NBas) 
 else
    write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
    stop
 endif

 ! transform and read 2-el integrals
 call tran4_full(NBas,MO,MO,fname)
 call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 if(Flags%ISHF==1.and.Flags%ISERPA==2.and.Mon%NELE==1) then

    ! Calculate FCI for 2-el systems
    call calc_fci(NBas,NInte1,NInte2,XOne,URe,TwoMO,Mon)

 endif

 ACAlpha=One
 ! GVB
 if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

   allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

   call ACABMAT0(ABPlus,ABMin,URe,Mon%Occ,XOne,TwoMO, &
                 NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,ACAlpha,1)
  
   ! reduce dim
   call reduce_dim('AB',ABPlus,ABMin,Mon)

   EigVecR = 0
   Eig = 0
   call ERPASYMM(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)

   if(EChck) then
      write(LOUT,'(/,1x,a)') 'ERPA-GVB ENERGY CHECK REQUESTED:'
      call EneERPA(ETot,ECorr,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,&
           Mon%Occ,XOne,Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
   endif

 ! CAS-SCF
 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then

   allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

   ! read 2-RDMs
   call read2rdm(Mon,NBas)

!  call execute_command_line('cp '//rdmfile// ' rdm2.dat')
   call system('cp '//rdmfile// ' rdm2.dat')

!   Gamma_2_AB CAN BE USED 
!   WITH DIFFERENT select_active
!   AND reducing dimensions!
!   call Gamma2_AB(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
!                  NBas,Mon%NDim,NInte1,NInte2,ACAlpha)

   ECASSCF = 0
!
! test semicoupled
!  ACAlpha = 0.000000001
! if(Mon%Monomer==1) then
!   ACAlpha=0.010  
! else
!   ACAlpha=0.0000001
! endif 

   call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
               Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)
 
   EigVecR = 0
   Eig = 0
   call ERPASYMM1(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
 
   !print*, 'Entering ERPAVEC...' 
   !call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
!
!   i = 1
!   dimV1 = 1
!   space1(1,dimV1) = i
!   tmp = Eig(i)
!   do while(i<Mon%NDimX)
!      i = i + 1
!      if(abs(Eig(i)-tmp)>thresh) then
!         space1(2,dimV1) = i-1
!         dimV1=dimV1+1
!         space1(1,dimV1) = i
!         tmp=Eig(i)
!      endif
!   enddo
!   space1(2,dimV1) = Mon%NDimX

  if(EChck) then
      ECorr=0
      call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
                     Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
      ECorr=Ecorr*0.5d0
  
      write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
            ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr
   endif

   ! uncoupled
   ! something deeply wrong with IndN here??????????????? !!!!!!

!   call test_resp_unc(Mon,URe,XOne,TwoMO,NBas,NInte1,NInte2,Flags%IFlag0) 

   allocate(EigY0(Mon%NDimX**2),EigY1(Mon%NDimX**2),&
            Eig0(Mon%NDimX),Eig1(Mon%NDimX))
 
   EigY0 = 0
   EigY1 = 0
   Eig0 = 0
   Eig1 = 0

   allocate(Mon%IndNT(2,Mon%NDim)) 
   Mon%IndNT=0 
   do i=1,Mon%NDim
      Mon%IndNT(1,i) = Mon%IndN(1,i)
      Mon%IndNT(2,i) = Mon%IndN(2,i)
   enddo
   call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
       EigY0,EigY1,Eig0,Eig1, &
       !Mon%IndNT,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)
       Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)

   ! dump uncoupled response
   call writeresp(EigY0,Eig0,propfile0)
   if(Flags%IFlag0==0) then
      call writeresp(EigY1,Eig1,propfile1)
   endif
  
   deallocate(Eig1,Eig0,EigY1,EigY0)

  elseif(Flags%ISERPA==2) then
  ! TD-APSG response: GVB, CASSCF, FCI

      if(Flags%ICASSCF==1.and.MON%NELE/=1) then
         write(LOUT,'(1x,a)') 'ERROR!!! FCI POSSIBLE ONLY FOR 2-EL MONOMERS!'
         stop 
      endif

      if(Flags%ICASSCF==1.and.Flags%ISHF==0) then
        ! read 2-RDMs
         call read2rdm(Mon,NBas)
        ! CAS-SCF
         call execute_command_line('cp '//rdmfile// ' rdm2.dat')
      elseif(Flags%ICASSCF==1.and.Flags%ISHF==1) then
         ! PINO
         call read2rdm(Mon,NBas)
         call init_pino(NBas,Mon,Flags%ICASSCF) 
      endif

      Mon%NDimN=0
      do i=1,NBas
         if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
      enddo
   
      allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
               CMAT(Mon%NDim**2),EMAT(NBas**2),EMATM(NBas**2),&
               DMAT(Mon%NDim*NBas),DMATK(Mon%NDim*NBas),&
               EigVecR(2*(Mon%NDimX+Mon%NDimN)*2*(Mon%NDimX+Mon%NDimN)),&
               Eig(2*(Mon%NDimX+Mon%NDimN))) 
 
      !Mon%NDimN = 0 
      CMAT=0
      call APSG_NEST(ABPlus,ABMin,CMAT,EMAT,EMATM,DMAT,DMATK,&
                     URe,Mon%Occ,XOne,TwoMO,&
                     NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,Flags%ISERPA)
   
      !reduce dimensions
      call reduce_dim('AB',ABPlus,ABMin,Mon)
      call reduce_dim('D',DMAT,DMATK,Mon)
      call reduce_dim('E',EMAT,EMATM,Mon)
   
      !    write(*,*) Mon%NGem,Mon%NDimN,Mon%NDimX
      !    write(*,*) Mon%Occ,NBas
      EigVecR = 0
      Eig = 0
      call PINOVEC(EigVecR,Eig,ABPlus,ABMin,DMAT,DMATK,EMAT,EMATM, &
                   Mon%Occ,NBas,Mon%NDimX,Mon%NDimN)
   
      ETot = 0
      call EnePINO(ETot,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne, &
                   Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NDimN)

      deallocate(CMAT,EMAT,EMATM,DMAT,DMATK)

 endif

! dump response
 call writeresp(EigVecR,Eig,propfile)

 close(ione)
 deallocate(work1,work2,XOne,TwoMO,URe)
 deallocate(ABPlus,ABMin,EigVecR,Eig)

end subroutine calc_resp_full

subroutine calc_fci(NBas,NInte1,NInte2,XOne,URe,TwoMO,Mon)
implicit none

integer,intent(in) :: NBas, NInte1, NInte2
double precision   :: URe(NBas,NBas),XOne(NInte1),TwoMO(NInte2)
type(SystemBlock)  :: Mon
double precision   :: ETot
double precision, allocatable :: Ctmp(:,:),Dtmp(:,:), &
                                 NSymMO(:), TwoMOt(:)
integer :: EigNum

 write(LOUT,'(/,1x,a)') 'ENTERING FCI FOR 2-EL SYSTEMS...'
   
 allocate(Ctmp(NBas,NBas),Dtmp(NBas,NBas),NSymMO(NBas),TwoMOt(NInte2))

 NSymMO = 1
 ETot=0
   
 if(Mon%Monomer==1) then
    EigNum=1
 elseif(Mon%Monomer==2) then
    EigNum=1
 endif

 call OptTwo1(ETot,Mon%PotNuc,URe,Mon%Occ,XOne,TwoMO,NSymMO, &
              Mon%CICoef,NBas,NInte1,NInte2,EigNum)
        
 ! Transform orbitals from MO to NO 
 call dgemm('N','T',NBas,NBas,NBas,1d0,Mon%CMO,NBas,URe,NBas,0d0,Ctmp,NBas)

 Dtmp=0
 ! call get_den(NBas,Ctmp,2d0*Mon%Occ,Dtmp)
 ! print*, 'Trace test',trace(Dtmp,NBas)
 ! call print_mo(Ctmp,NBas,'MONOMER X')

 ! Save NO orbitals
 Mon%CMO = Ctmp

 ! MO->NO transform 2-el ints
 TwoMOt(1:NInte2) = 0
 call TwoNO(TwoMOt,URe,TwoMO,NBas,NInte2)

 ! Save transformed ints
 TwoMO = TwoMOt

 deallocate(TwoMOt,NSymMO,Ctmp,Dtmp)

end subroutine calc_fci

subroutine init_pino(NBas,Mon,ICASSCF)
implicit none

integer,intent(in) :: NBas, ICASSCF
type(SystemBlock) :: Mon
integer :: i,j,ij,ind

 ! set IGem
 Mon%NGem=1
 do j=1,NBas
    Mon%IGem(j)=1
    if(Mon%Occ(1)==0) Mon%IGem(j) = 2 
    if(Mon%Occ(1)==0) Mon%NGem = 2
 enddo

 call SaptInter(NBas,Mon,ICASSCF)

! look-up table 
 do i=1,Mon%NELE
   Mon%IndAux(i)=0
 enddo
 do i=1+Mon%NELE,NBas
   Mon%IndAux(i)=2
 enddo

 Mon%NAct = 1
 if(Mon%NAct.ne.0) then
   Mon%icnt = 0
   do i=1,NBas  
      if(Mon%Occ(i).gt.0d0) then
         Mon%IndAux(i) = 1
         Mon%icnt = Mon%icnt + 1
      endif
   enddo      
 endif

 ! set generalized "occupied" = num0 + num1
 ! and "virtual" = num1 + num2 indices
 Mon%num0 = 0 
 do i=1,nbas
    if(Mon%IndAux(i)/=0) exit
    Mon%num0 = Mon%num0 + 1
 enddo
 Mon%num2 = 0
 do i=nbas,1,-1
    if(Mon%IndAux(i)/=2) exit
       Mon%num2 = Mon%num2 + 1
 enddo
 Mon%num1 = nbas - Mon%num0 - Mon%num2

! do i=1,NBas
!    if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
! enddo

 ij=0
 ind=0
 do i=1,NBas
    do j=1,i-1  

    ij=ij+1

    if(Mon%IndAux(i)+Mon%IndAux(j).Ne.0.and.Mon%IndAux(i)+Mon%IndAux(j).Ne.4) then 

       if((Mon%IGem(i).ne.Mon%IGem(j)).And.(Mon%IndAux(i).Eq.1).And.(Mon%IndAux(j).Eq.1) &
         .and.(Abs(Mon%Occ(i)-Mon%Occ(j))/Mon%Occ(i).Lt.1.D-2) ) Then

          write(*,*)"Discarding nearly degenerate pair",i,j

       else

       ind=ind+1
       Mon%IndX(ind)=ind!ij
       Mon%IndN(1,ind)=i
       Mon%IndN(2,ind)=j
  
       endif
    endif

    enddo
 enddo

 Mon%NDimX = ind
! write(LOUT,*) Mon%NDim,Mon%NDimX,Mon%NDimN
 write(LOUT,*) 'init PINO!!!'


end subroutine init_pino

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

 call get_den(NBas,A%CMO,A%Occ,2d0,Pa)
 call get_den(NBas,B%CMO,B%Occ,2d0,Pb)

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

!subroutine select_resp(EigVec,Eig,Mon,NBas,NDimEx)
!implicit none
!! is this even ok?
!type(SystemBlock) :: Mon
!integer :: NBas,NDimEx
!double precision :: EigVec(2*NDimEx,2*NDimEx),Eig(2*NDimEx)
!double precision,allocatable :: TmpVec(:,:),TmpEig(:)
!integer :: i,j
!integer :: pq,ip,iq
!
! !NDimEx = Mon%NDimX+Mon%NDimN
! allocate(TmpEig(NDimEx),TmpVec(NDimEx,NDimEx))
!
! TmpEig = 0
! TmpVec = 0
! j = 0
! do i=1,2*NDimEx
!    if(Eig(i).gt.1d0.and.Eig(i).lt.1d20) then
!       j = j + 1
!       TmpEig(j) = Eig(i)
!       TmpVec(1:NDimEx,j) = EigVec(1:NDimEx,i)
!    endif
! enddo  
! 
! print*, 'TmpEig'
! do i=1,NDimEx
!    write(LOUT,*) i,TmpEig(i)
! enddo
!
!! not too much... 
!! this is so complicated...
!! maybe a way to simplify it?
!  print*, '1st'
!  do i=1,NDimEx
!    do pq=1,2*NDimEx
!       if(pq.le.Mon%NDimX) then
!          ip=Mon%IndN(1,pq)
!          iq=Mon%IndN(2,pq)
!       else
!          ip=pq-Mon%NDimX
!          iq=ip
!       endif
!       if(ip.gt.iq) then
!          if(Eig(pq).gt.1d0.and.Eig(pq).lt.1d20) then
!            write(6,*) ip,iq,EigVec(i,pq)
!          endif
!       endif
!    enddo
!  enddo 
! 
!  print*, '2nd'
!! map(pq) here?
!  do i=1,NDimEx
!    do pq=1,NDimEx
!       if(pq.le.Mon%NDimX) then
!          ip=Mon%IndN(1,pq)
!          iq=Mon%IndN(2,pq)
!       else
!          ip=pq-Mon%NDimX
!          iq=ip
!       endif
!       if(ip.gt.iq) then
!          write(6,*) ip,iq,TmpVec(i,pq)
!       endif
!    enddo
!  enddo 
!
!
!
! deallocate(TmpEig,TmpVec)
!
!end subroutine select_resp

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

!call read_syminf(SAPT%monA,SAPT%monB,nbas)

call gen_swap_rows(mat,nbas,SAPT%monA%NSym,&
                   SAPT%monA%NMonBas,SAPT%monB%NMonBas)
call gen_swap_cols(mat,nbas,SAPT%monA%NSym,&
                   SAPT%monA%NMonBas,SAPT%monB%NMonBas)

!call swap_rows(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)
!call swap_cols(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)

end subroutine arrange_oneint

subroutine arrange_mo(mat,nbas,SAPT)
implicit none

type(SaptData) :: SAPT 
!integer :: NOrbA,NOrbB
integer :: nbas
double precision :: mat(nbas,nbas)

!print*, NOrbA,NOrbB,size(mat)
call gen_swap_rows(mat,nbas,SAPT%monA%NSym,&
                   SAPT%monA%NMonBas,SAPT%monB%NMonBas)

!call swap_rows(NOrbA,NOrbB,mat)

end subroutine arrange_mo

subroutine read_syminf(A,B,nbas)
! reads number of basis functions on each monomer
! from SYMINFO file!
implicit none

type(SystemBlock) :: A, B
integer :: nbas
integer :: iunit,ios
integer :: ibas,icen,last_ibas,last_icen
integer :: irep,offset
logical :: ex,dump
integer :: tmp
integer :: ACenTst, ACenBeg, ACenEnd

!print*, A%NCen, B%NCen
!print*, A%UCen, B%UCen

inquire(file='SYMINFO_B',EXIST=ex)

if(ex) then
   open(newunit=iunit,file='SYMINFO_B',status='OLD',&
        form='FORMATTED')
   read(iunit,*)
   read(iunit,*)
   offset=0
   irep=1
   read(iunit,'(i5,i6)',iostat=ios) last_ibas,last_icen
   do
     read(iunit,'(i5,i6)',iostat=ios) ibas,icen
     if(ios/=0) then
        A%NMonBas(irep)=last_ibas-offset 
        exit
     elseif(icen/=last_icen) then
          if(last_icen==B%UCen) then
             B%NMonBas(irep) = last_ibas-offset
             offset=last_ibas
          elseif(icen==1) then
             A%NMonBas(irep) = last_ibas-offset
             offset=last_ibas
             irep = irep + 1
          endif
     endif
     last_ibas=ibas
     last_icen=icen
  enddo

!   dump = .TRUE.
!   do
!     read(iunit,'(i5,i6)',iostat=ios) ibas,icen
!     if((icen.gt.B%NCen).and.dump) then 
!        B%NMonOrb = ibas-1
!        ACenBeg = icen
!        dump = .FALSE.
!     elseif(.not.dump.and.ios==0) then
!        last_icen = icen
!        last_ibas = ibas
!     elseif(ios/=0) then
!        tmp = last_ibas
!        ACenEnd = last_icen
!        exit
!     endif
!     !write(*,*) ibas,icen,dump
!   enddo

! check later somehow...?   
!   ACenTst = ACenEnd - ACenBeg + 1
!   if(ACenTst/=A%NCen) then
!      write(LOUT,'(1x,a)') 'ERROR! MISMATCH IN NUMBER OF ATOMS FOR MONONOMER A!'
!      write(LOUT,'(1x,a,i3,4x,a,i3)') 'INPUT: ', A%NCen, 'SYMINFO_B: ', ACenTst
!      stop
!   endif
!   A%NMonOrb = tmp - B%NMonOrb

   close(iunit)
else
   write(LOUT,'(1x,a)') 'ERROR! MISSING SYMINFO_B FILE!'
   stop
endif
! HERE
! maybe add test: tmp vs. NBasis? 

end subroutine read_syminf

subroutine gen_swap_rows(mat,nbas,nsym,nA,nB)
implicit none

integer,intent(in) :: nbas,nsym,nA(8),nB(8)
double precision :: mat(nbas,nbas)
double precision :: work(nbas,nbas)
integer :: irep,iA,iB,iAB,offset

offset = 0

do irep=1,nsym

   iA = nA(irep)
   iB = nB(irep)
   iAB = iA + iB

   work(1:iA,:) = mat(offset+iB+1:offset+iAB,:)
   work(iA+1:iAB,:) = mat(offset+1:offset+iB,:)

   mat(offset+1:offset+iAB,:) = work(1:iAB,:)

   offset = offset + iAB
enddo


end subroutine gen_swap_rows

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

subroutine gen_swap_cols(mat,nbas,nsym,nA,nB)
implicit none

integer,intent(in) :: nbas,nsym,nA(8),nB(8)
double precision :: mat(nbas,nbas)
double precision :: work(nbas,nbas)
integer :: irep,iA,iB,iAB,offset

offset = 0

do irep=1,nsym

   iA = nA(irep)
   iB = nB(irep)
   iAB = iA + iB

   work(:,1:iA) = mat(:,offset+iB+1:offset+iAB)
   work(:,iA+1:iAB) = mat(:,offset+1:offset+iB)

   mat(:,offset+1:offset+iAB) = work(:,1:iAB)
   
   offset = offset + iAB

enddo

end subroutine gen_swap_cols

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

subroutine read_mo(cmo,nbasis,nsym,nbas,norb,nsiri,nmopun)

! in SAPT orbitals kept in AOMO order!
implicit none

integer,intent(in) :: nbasis,nsym,nbas(8),norb(8)
integer :: iunit,irep,idx
integer :: ncmot
!double precision :: cmo(norb,nbas)
double precision :: cmo(nbasis,nbasis)
character(*) :: nsiri,nmopun
logical :: isiri
character(60) :: line
integer :: i,j
integer :: off_i, off_j
double precision :: natocc(10)
double precision :: tmp(nbasis**2)

ncmot = sum(nbas(1:nsym)*norb(1:nsym))

inquire(file=nsiri,EXIST=isiri)

if(isiri) then

   open(newunit=iunit,file=nsiri,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

   call readlabel(iunit,'NEWORB  ')
   !read(iunit) cmo
   read(iunit) tmp(1:ncmot)

!  print*, norb,'NORB!!!!'
   cmo = 0  
   off_i = 0
   off_j = 0
   idx = 0

   do irep=1,nsym
      do j=off_j+1,off_j+norb(irep)
 
         do i=off_i+1,off_i+nbas(irep)   
            idx = idx + 1
            cmo(i,j) = tmp(idx)
         enddo
      
      enddo
      off_i = off_i + nbas(irep)
      off_j = off_j + norb(irep)
   enddo
 
 !   do i=1,norb
 !     do j=1,nbas
 !        cmo(j,i) = tmp((i-1)*nbas + j)
 !     end do
 !  end do

   write(LOUT,'(1x,a)') 'Orbitals read from '//nsiri 

else
    write(LOUT,'(1x,a)') 'FIX READING ORBITALS FROM DALTON.MOPUN'
!   print*, 'Achtung!!!',norb,nbasis
!   open(newunit=iunit,file=nmopun, &
!        form='FORMATTED',status='OLD')
!   read(iunit,'(a60)') line
!   !do j=1,norb     
!   do j=1,nbas     
!      read(iunit,'(4f18.14)') (cmo(i,j),i=1,nbasis)
!   enddo
!!   print*, line
!
!   write(LOUT,'(1x,a)') 'Orbitals read from '//nmopun 
endif

close(iunit)

! call print_sqmat(cmo,nbasis)

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
integer :: iunit,ios,i
integer :: NAct, INAct, NActS(8), INActS(8)
integer :: irep,TotEl,offset
double precision :: Occ(nbas), sum1, sum2
double precision :: potnuc, emy, eactiv, emcscf
integer :: istate, ispin, nactel, lsym
!integer :: nisht, nasht, nocct, norbt, nbast, nconf, nwopt, nwoph
integer :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
           NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
           NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBASM(8) 

 allocate(mon%CICoef(nbas),mon%IGem(nbas),mon%Occ(nbas))
 ioccsir=.false.
 print*, 'SIRI?',exsiri
 !exsiri=.false.
 if(exsiri) then

    rewind(isiri) 
    read (isiri) 
    read (isiri) potnuc,emy,eactiv,emcscf, &
                 istate,ispin,nactel,lsym
!    read (isiri) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
    read (isiri) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
              NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
              NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBASM 
  
!    print*,    potnuc,emy,eactiv,emcscf, &
!               istate,ispin,nactel,lsym
 !   print*, 'READM-TEST'
 !   write (*,*) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
   
    mon%NAct  = nasht
    mon%INAct = nisht
!    print*, nasht, nisht, norbt, nactel 
!    print*, 'nash',nash
!    print*, 'wtf?',sum(nash), sum(nish)

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
    INActS = 0
    NActS  = 0
    open(newunit=iunit,file=occfile,form='FORMATTED',status='OLD') 

    read(iunit,*) INAct, NAct
    INAct = INAct/2
    read(iunit,*) (Occ(i),i=1,INAct+NAct)
    sum2 = 0d0
    do i=1,INAct+NAct
       Occ(i) = Occ(i)/2d0
       sum2 = sum2 + Occ(i)
    enddo

    ! (in)active orbs in each symmetry
    read(iunit,*,iostat=ios) (INActS(i),i=1,mon%NSym)
    if(ios==0) then
       read(iunit,*) (NActS(i),i=1,mon%NSym)
       mon%INActS(1:mon%NSym) = INActS(1:mon%NSym)
       mon%NActS(1:mon%NSym)  = NActS(1:mon%NSym)
    endif

    if(mon%NSym.gt.1) then
      call sort_sym_occ(nbas,mon%NSym,INAct,NAct,INActS,NActS,Occ)
    endif

!    offset = 0
!    sum2 = 0d0
!    do irep=1,mon%NSym
!       read(iunit,*,iostat=ios) INActS(irep), NActS(irep)
!       INActS(irep) = INActS(irep)/2
!       TotEl = INActS(irep) + NActS(irep)
!       ! print*, 'TotEl',TotEl,offset
!       if(TotEl.gt.0) read(iunit,*) (Occ(i),i=offset+1,offset+TotEl)
!       do i=offset+1,offset+TotEl
!          Occ(i) = Occ(i)/2d0
!          sum2 = sum2 + Occ(i)
!       enddo
!       offset = offset + mon%NSymOrb(irep)
!       ! print*, 'offs',offset
!    enddo
!    INAct = sum(INActS)
!    NAct  = sum(NActS)
!    mon%INActS(1:mon%NSym) = INActS(1:mon%NSym)
!    mon%NActS(1:mon%NSym)  = NActS(1:mon%NSym)

    if(.not.ioccsir) then
    !   if(Abs(sum2-mon%XELE).gt.1.0d-8) then
    !      write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'  
    !      write(LOUT,'(1x,a,1x,f10.6,5x,a,i3)') 'SUM(OCC): ', sum2, 'MONOMER: ', mon%Monomer
    !      write(LOUT,'(1x,a)') 'CHECK occupations.dat!'  
    !      stop
    !   endif
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

   ! offset=0
   ! do irep=1,mon%NSym
   !    TotEl = mon%NActS(irep)+mon%INActS(irep)
   !    mon%IGem(offset+1:offset+TotEl) = 1
   !    mon%IGem(offset+TotEl+1:offset+mon%NSymOrb(irep)) = 2
   !     
   !    offset = offset + mon%NSymOrb(irep)
   ! enddo  

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

subroutine select_active(mon,nbas,Flags)
implicit none

type(SystemBlock) :: mon
type(FlagsData) :: Flags
integer :: nbas!, ICASSCF, ISHF, IFlCore
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
 if(Flags%ICASSCF==0) then
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
 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0) then
   write(LOUT,'()')
   if(mon%Monomer==1) write(LOUT,'(1x,a)') 'Monomer A' 
   if(mon%Monomer==2) write(LOUT,'(1x,a)') 'Monomer B' 
   do i=1,nbas
      if(mon%Occ(i).lt.1d0.and.mon%Occ(i).ne.0d0) then
         ! HERE!!! ACTIVE!!!! 
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
 print*, 'TTESSSTT:', mon%num0,mon%num1,mon%num2, nbas

! active pairs
 allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))

 mon%IPair(1:nbas,1:nbas) = 0

 if(Flags%ICASSCF==0) then

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
                 if(Flags%IFlCore==1.or.&
                   (Flags%IFlCore==0.and.&
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

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then


    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1
   
          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
                ! do not correlate active degenerate orbitals from different geminals 
                if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1)  & 
                 .and.&
                 (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
               
                   write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
                else
                  ! if IFlCore=0 exclude core (inactive) orbitals
                  if(Flags%IFlCore==1.or.&
                    (Flags%IFlCore==0.and.&
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

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==2) then

    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1
   
          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
                ! do not correlate active degenerate orbitals from different geminals 
                !if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1)  & 
                ! .and.&
                ! (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
                ! write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
                !else
                 ! if IFlCore=0 exclude core (inactive) orbitals
                 if(Flags%IFlCore==1.or.&
                   (Flags%IFlCore==0.and.&
                    mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
      
                    ind = ind + 1
                    mon%IndX(ind) = ind
                    mon%IndN(1,ind) = i
                    mon%IndN(2,ind) = j
                    mon%IPair(i,j) = 1
                    mon%IPair(j,i) = 1
      
                 endif
              ! endif
      
          endif
   
       enddo
    enddo
 
 
 endif
 mon%NDimX = ind

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

subroutine sort_sym_occ(nbas,nsym,INAct,NAct,INActS,NActS,Occ)
implicit none

integer,intent(in) :: nbas, nsym, INAct, NAct 
integer,intent(in) :: INActS(8), NActS(8)
double precision,intent(inout) :: Occ(nbas)
integer :: TotEl
integer :: i,ii
integer,allocatable :: ICpy1(:),ICpy2(:)
double precision :: OccOrd(nbas)

 TotEl = INAct + NAct
 OccOrd = 0

 allocate(ICpy1(TotEl),ICpy2(TotEl))

 ICpy1 = 0
 ICpy2 = 0

 do ii=1,TotEl 

    ! inactive 
    do i=1,TotEl
       if(ICpy2(i).eq.0.and.ICpy1(ii).eq.0.and.Occ(i).eq.2.0D0) then
          ICpy2(i)  = 1
          ICpy1(ii) = 1
          OccOrd(ii) = Occ(i)
       endif
    enddo

    ! active
    if(ICpy1(ii).eq.0) then
       do i=1,TotEl
          if(ICpy2(i).eq.0.and.ICpy1(ii).eq.0) then
             ICpy2(i)  = 1
             ICpy1(ii) = 1
             OccOrd(ii) = Occ(i)
          endif
       enddo
    endif

 enddo

! check
! do i=1,nbas
!    print*, i,Occ(i),OccOrd(i)
! enddo

 Occ = OccOrd
   
 deallocate(ICpy2,ICpy1)

end subroutine sort_sym_occ

subroutine sort_sym_mo(CMO,nbas,mon)
implicit none

type(SystemBlock) :: mon
integer,intent(in) :: nbas
double precision,intent(inout) :: CMO(nbas,nbas)
integer,allocatable :: ICpy1(:),ICpy2(:)
integer,allocatable :: LabelAct(:),LabelIAct(:)
double precision,allocatable :: COrd(:,:)
integer :: TotEl,TotElIrep,irep,idx
integer :: i,j,ii

 TotEl = mon%INAct + mon%NAct

 allocate(ICpy1(nbas),ICpy2(nbas))
 allocate(LabelAct(nbas),LabelIAct(nbas),COrd(nbas,nbas))

 ICpy1 = 0
 ICpy2 = 0

 ! make labels
 idx = 0
 do irep=1,mon%NSym
    do j=1,mon%NSymOrb(irep)

       idx = idx + 1
       LabelAct(idx) = 0
       TotElIrep = mon%INActS(irep)+mon%NActS(irep)

       if(j.gt.mon%INActS(irep).and.j.le.TotElIrep) then
          LabelAct(idx) = 1
       endif
       LabelIAct(idx)=0

       if(j.le.mon%INActS(irep)) LabelIAct(irep)=1

    enddo
 enddo

 do ii=1,nbas 

   ! inactive 
   do i=1,nbas  
      if(LabelIAct(i).eq.1.and.ICpy2(i).eq.0.and.ICpy1(ii).eq.0) then
         ICpy2(i)  = 1
         ICpy1(ii) = 1

         do j=1,nbas  
!            URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
             COrd(j,ii) = CMO(j,i)
         enddo
      endif
   enddo

   ! active
   if(ICpy1(ii).eq.0) then
      do i=1,nbas  
         if(LabelAct(i).eq.1.and.ICpy2(i).eq.0.and.ICpy1(ii).eq.0) then
            ICpy2(i)  = 1
            ICpy1(ii) = 1
            do j=1,nbas  
               COrd(j,ii) = CMO(j,i)
!               URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
            enddo
         endif
      enddo
   endif

   ! virtual
   if(ICpy1(ii).Eq.0) then
      do i=1,nbas  
         if(ICpy2(i).eq.0.and.ICpy1(ii).eq.0) then
            ICpy2(i)  = 1
            ICpy1(ii) = 1
            do j=1,nbas
               COrd(j,ii) = CMO(j,i)
               ! URe2(II+(J-1)*NBasis)=URe1(I+(J-1)*NBasis)
            enddo
         endif
      enddo
   endif

 enddo

 CMO = COrd

 deallocate(COrd,LabelIAct,LabelAct)
 deallocate(ICpy2,ICpy1)

end subroutine sort_sym_mo

subroutine readener(nbasis,mon,isiri)
implicit none

type(SystemBlock) :: mon
integer :: nbasis, isiri
integer :: MMORBT
integer :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
           NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
           NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBAS(8) 
integer :: i,idx,irep,offset
double precision,allocatable :: fock(:)

! set dimensions
 rewind(isiri)

 read (isiri)
 read (isiri) 
 read (isiri) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
              NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
              NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBAS 
 read(isiri)
 read(isiri)
 read(isiri)
 read(isiri)
 read(isiri)


 MMORBT = max(4,NNORBT)
 allocate(fock(MMORBT),mon%OrbE(NORBT))

 read(isiri) fock 
 ! orb energies: diag of Fock
 offset = 0
 idx = 0
 do irep=1,NSYM

    do i=1,NORB(irep)
       idx = idx + 1
       mon%OrbE(idx) = fock(offset+i*(i+1)/2)
    enddo 

       offset = offset + NORB(irep)*(NORB(irep)+1)/2
 enddo

! print*, mon%OrbE
!
! do i=1,NORBT
!    mon%OrbE(i) = fock(i+i*(i-1)/2)
! enddo

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
    write(LOUT,'()')
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

if(allocated(SAPT%monA%RDM2)) then
  deallocate(SAPT%monA%RDM2,SAPT%monA%RDM2Act)
  deallocate(SAPT%monA%Ind2)
endif
if(allocated(SAPT%monB%RDM2)) then
  deallocate(SAPT%monB%RDM2,SAPT%monB%RDM2Act)
  deallocate(SAPT%monB%Ind2)
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

