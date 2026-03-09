module sapt_interface_io

use types
use timing
use tran
use sorter
use tran_Chol
use gammcor_integrals, only : TCholeskyVecs, TCholeskyVecsOTF, TAOBasis, TSystem, &
                              chol_CoulombMatrix, chol_MOTransf_TwoStep, &
                              chol_gammcor_Rkab, CholeskyOTF_ao_vecs, &
                              ORBITAL_ORDERING_DALTON, ORBITAL_ORDERING_MOLPRO, &
                              ORBITAL_ORDERING_ORCA, ORBITAL_ORDERING_PYSCF, &
                              BECKE_PARAMS_MEDIUM, auto2e_init, &
                              auto2e_interface_C, &
                              sys_read_xyz, basis_NewAOBasis, &
                              Becke_MolecularGrid, gridfunc_orbitals, &
                              CholeskyOTF_Fock_MO_v2

use abmat
use read_external
use trexio

implicit none

contains

subroutine onel_molpro(mon,NBasis,MonBlock,SAPT)
implicit none

type(SaptData)     :: SAPT
type(SystemBlock)  :: MonBlock
integer,intent(in) :: mon,NBasis

integer :: NSq,NInte1
integer                       :: ione,ios,NSym,NBas(8),ncen
double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
double precision, allocatable :: Kmat(:)
double precision, allocatable :: work1(:),work2(:)
character(8)                  :: label
character(:),allocatable      :: infile,outfile

!set dimensions
NSq = NBasis*NBasis
NInte1 = NBasis*(NBasis+1)/2

if(mon==1) then
  infile  = 'AOONEINT_A'
  outfile = 'ONEEL_A'
elseif(mon==2) then
  infile  = 'AOONEINT_B'
  outfile = 'ONEEL_B'
endif

allocate(work1(NInte1),work2(NSq))
allocate(Hmat(NSq),Vmat(NSq),Smat(NSq))
! read and dump 1-electron integrals
open(newunit=ione,file=infile,access='sequential',&
     form='unformatted',status='old')
read(ione)
read(ione) NSym,NBas(1:NSym)
read(ione) MonBlock%PotNuc

do
  read(ione,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL ISORDK   not found!'
     stop
  endif
  if(label=='ISORDK  ') then
     read(ione) ncen
     allocate(MonBlock%charg(ncen),MonBlock%xyz(ncen,3))
     read(ione) MonBlock%charg(1:ncen),MonBlock%xyz(1:ncen,1:3)
     exit
  endif
enddo
 !print*, 'ncen',MonBlock%charg(1:ncen)
 !print*, 'ncen',MonBlock%xyz(1:ncen,1:3)

 close(ione)

 call readoneint_molpro(work1,infile,'ONEHAMIL',.false.,NInte1)
 call square_oneint(work1,Hmat,NBasis,NSym,NBas)
 !call print_sqmat(Hmat,NBasis)

 call readoneint_molpro(work1,infile,'POTENTAL',.false.,NInte1)
 call square_oneint(work1,Vmat,NBasis,NSym,NBas)
 !call print_sqmat(Vmat,NBasis)

 call readoneint_molpro(work1,infile,'OVERLAP ',.false.,NInte1)
 call square_oneint(work1,Smat,NBasis,NSym,NBas)
 !call print_sqmat(Smat,NBasis)

 MonBlock%NSym = NSym
 MonBlock%NSymBas(1:NSym) = NBas(1:NSym)

 !square form
 call writeoneint(outfile,NSq,Smat,Vmat,Hmat)

 deallocate(work2,work1)
 deallocate(Smat,Vmat,Hmat)

end subroutine onel_molpro

subroutine onel_dalton(mon,NBasis,NSq,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData) :: SAPT
 type(SystemBlock) :: MonBlock

 integer,intent(in) :: mon,NBasis,NSq,NInte1

 integer :: ione,NSym,NBas(8),ncen
 integer :: i,ncenA,ncenB
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
 double precision, allocatable :: work1(:),work2(:)
 character(:),allocatable :: infile,outfile

 if(mon==1) then
   infile =  'AOONEINT_A'
   outfile = 'ONEEL_A'
 elseif(mon==2) then
   infile =  'AOONEINT_B'
   outfile = 'ONEEL_B'
 endif

 allocate(work1(NInte1),work2(NSq))
 allocate(Hmat(NSq),Vmat(NSq),Smat(NSq))
! read and dump 1-electron integrals
 open(newunit=ione,file=infile,access='sequential',&
      form='unformatted',status='old')
 read(ione)
 read(ione) NSym,NBas(1:NSym),MonBlock%PotNuc

 ! HERE!!! temp!
 MonBlock%NSymOrb(1:NSym) = NBas(1:NSym)

 call readlabel(ione,'ONEHAMIL')
 call readoneint_dalton(ione,work1)
 call square_oneint(work1,Hmat,NBasis,NSym,NBas)

 call readlabel(ione,'KINETINT')
 call readoneint_dalton(ione,work1)
 call square_oneint(work1,work2,NBasis,NSym,NBas)
 Vmat(:) = Hmat - work2

 call readlabel(ione,'OVERLAP ')
 call readoneint_dalton(ione,work1)
 call square_oneint(work1,Smat,NBasis,NSym,NBas)

 allocate(MonBlock%charg(maxcen),MonBlock%xyz(maxcen,3))
 call readlabel(ione,'ISORDK  ')
 read(ione)
 read(ione) MonBlock%charg,ncen,MonBlock%xyz

 if(mon==2) then

    ! check is B and A monomers were switched
    if(MonBlock%charg(1)/=0d0) MonBlock%switchAB = .true.

    !print*, 'charge before switch'
    !write(LOUT,*) MonBlock%charg(1:ncen)
    ! if Gh(A)-B :
    !  a) adapt xyz coords of B
    !  b) adapt chrge of B
    if(MonBlock%charg(1)==0d0) then
       ncenA = 0
       do i=1,ncen
          if(MonBlock%charg(i)==0d0) ncenA = ncenA + 1
       enddo
       ncenB=ncen-ncenA
       !print*, 'ncenA,ncenB',ncenA,ncenB,ncen
       do i=1,ncenB
          MonBlock%xyz(i,:) = MonBlock%xyz(i+ncenA,:)
          MonBlock%charg(i) = MonBlock%charg(i+ncenA)
       enddo
       MonBlock%charg(ncenB+1:ncen) = 0d0
    endif

 endif

! print*, 'MONO-A',ncen
! write(LOUT,*) SAPT%monA%charg(1:ncen)
! do i=1,ncen
!    write(LOUT,*) SAPT%monA%xyz(i,:)
! enddo

! write(*,*) 'VA'
! call print_sqmat(Va,NBasis)
! call print_diag(Va,NBasis)

 close(ione)

 MonBlock%NSym = NSYm
 MonBlock%NSymBas(1:NSym) = NBas(1:NSym)

 if(mon==2) then
 ! rearrange in V: (B,A) -> (A,B)
    !call read_syminf(SAPT%monA,SAPT%monB,NBasis)
    if(MonBlock%switchAB) then
       call read_syminf_dalton(SAPT%monA%NSym,SAPT%monB%NSym,SAPT%monB%UCen, &
                               SAPT%monA%NSymOrb,SAPT%monB%NSymOrb,&
                               SAPT%monA%NMonBas,SAPT%monB%NMonBas)

       call arrange_oneint(Smat,NBasis,SAPT)
       call arrange_oneint(Vmat,NBasis,SAPT)
       call arrange_oneint(Hmat,NBasis,SAPT)
    endif
 endif

 ! square form
 call writeoneint(outfile,NSq,Smat,Vmat,Hmat)

 deallocate(work2,work1)
 deallocate(Hmat,Vmat,Smat)

end subroutine onel_dalton

subroutine readocc_dalton(NBasis,Mon,Flags)
implicit none

type(SystemBlock)  :: Mon
type(FlagsData)    :: Flags
integer,intent(in) :: NBasis

integer                  :: NSym,NOrbt,NBasist,NCMOt,NOcc(8),NOrbs(8)
integer                  :: i,isiri
double precision         :: potnuc,emy,eactiv,emcscf
logical                  :: exsiri,noSiri,noOccu
character(:),allocatable :: occfile,ifcfile,siriusfile,coefile


 if(Mon%Monomer==1) then
   coefile='coeff_A.dat'
   occfile='occupations_A.dat'
   ifcfile='SIRIFC_A'
   siriusfile='SIRIUS_A.RST'
 elseif(Mon%Monomer==2) then
   coefile='coeff_B.dat'
   occfile='occupations_B.dat'
   ifcfile='SIRIFC_B'
   siriusfile='SIRIUS_B.RST'
 endif

 inquire(file=ifcfile,EXIST=exsiri)
 if(exsiri) then
    call read_orbinf_dalton(ifcfile,NSym,Mon%NOrb,Mon%NSymOrb)
 else
    NBasist = NBasis
 endif
 !print*, 'readocc: NSym, NOrb', NSym,Mon%NOrb

 if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.(.not.Mon%ISHF)) then

    ! CASSCF

    !if(exsiri) close(isiri)

    call readocc_cas_siri(Mon,NBasis,noSiri)
    if(noSiri) call readocc_cas_occu(Mon,NBasis,noOccu)

    if(nosiri.and.nooccu) then
       write(lout,'(1x,a)') &
             'ERROR in readocc_dalton! No SIRIFC or occupations files!'
       stop
    endif

 elseif(Flags%ICASSCF==1.and.(Flags%ISHF==1.or.Mon%ISHF)) then

    if(exsiri) close(isiri)

    ! HARTREE-FOCK
    !call readmulti(NBasis,Mon,.true.,exsiri,isiri,occfile,siriusfile)
    call readocc_hf_siri(Mon,NBasis)

 elseif(Flags%IGVB==1) then

    ! GVB
    call readgvb(Mon,NBasis,coefile)

 endif

 if(Flags%ICASSCF==1) then
    ! construct IGem
    allocate(Mon%IGem(NBasis))
    if(Mon%INAct==0) then
       Mon%NGem = 2
       Mon%IGem(1:Mon%NAct+Mon%INAct)        = 1
       Mon%IGem(Mon%NAct+Mon%INAct+1:NBasis) = 2
    else
       Mon%NGem = 3
       Mon%IGem(1:Mon%INAct) = 1
       Mon%IGem(Mon%INAct+1:Mon%INAct+Mon%NAct) = 2
       Mon%IGem(Mon%INAct+Mon%NAct+1:NBasis)    = 3
    endif

    ! construct CICoef
    allocate(Mon%CICoef(NBasis))
    do i=1,NBasis
       Mon%CICoef(i)=sqrt(Mon%Occ(I))
       if(Mon%Occ(i).lt.0.5d0) Mon%CICoef(i)=-Mon%CICoef(i)
    enddo
 endif

 !if(exsiri) close(isiri)

end subroutine readocc_dalton

subroutine readocc_molpro(NBasis,Mon,OrbAux,OneRdm,Flags)
implicit none
!
! OrbAux  :: on output C(MO,NO)
! OneRdm  :: on output 1-RDM in AO
!
type(SystemBlock)  :: Mon
type(FlagsData)    :: Flags
integer,intent(in) :: NBasis

integer :: NAct,NOccup
integer :: NInte1,HlpDim
integer :: i,info
double precision :: Tmp
double precision :: OrbAux(NBasis,NBasis), &
                    OneRdm(NBasis*(NBasis+1)/2)
double precision,allocatable :: EVal(:)
double precision,allocatable :: work(:)
character(:),allocatable :: mname
character(:),allocatable :: rdmfile

 NInte1 = NBasis*(NBasis+1)/2
 HlpDim = max(NBasis**2,3*NBasis)

 if(Mon%Monomer==1) then
   rdmfile = '2RDMA'
   mname   = 'A'
 elseif(Mon%Monomer==2) then
   rdmfile = '2RDMB'
   mname   = 'B'
 endif

 call read_nact_molpro(NAct,rdmfile)

 allocate(Mon%CICoef(NBasis),Mon%IGem(NBasis),Mon%Occ(NBasis))
 allocate(work(HlpDim),EVal(NBasis))
 OneRdm = 0d0
 EVal   = 0d0
 call read_1rdm_molpro(OneRdm,Mon%InSt(1,1),Mon%InSt(2,1),&
                       Mon%ISpinMs2,rdmfile,Mon%IWarn,NBasis)

 call triang_to_sq2(OneRdm,OrbAux,NBasis)
 call Diag8(OrbAux(1:NAct,1:NAct),NAct,NAct,Eval(1:NAct),work)
 !call Diag8(OrbAux,NBasis,NBasis,Eval,work)

! KP : it may happen that an active orbital has a negative tiny occupation. set it to a positive
 do i=1,NBasis
 Eval(i)=Abs(Eval(i))
 enddo
! call dsyev('V','U',NBasis,OrbAux,NBasis,EVal,work,3*NBasis,info)
 call SortOcc(EVal,OrbAux(1:NAct,1:NAct),NAct)
 !call SortOcc(EVal,OrbAux,NBasis)

! read NAct from 1RDM
 if(Mon%NActFromRDM) Mon%NAct = 0
 Tmp = 0
 do i=1,NBasis
    Tmp = Tmp + EVal(i)
    !if(Mon%NActFromRDM.and.EVal(i)>0.d0) Mon%NAct = Mon%NAct + 1
    if(Mon%NActFromRDM.and.EVal(i)>1d-10) Mon%NAct = Mon%NAct + 1
 enddo

! test NAct from 1RDM
 !call read_nact_molpro(nact,rdmfile)
 if(Mon%NAct/=nact) then
    write(lout,'(1x,2a)') 'Warning! In monomer ', mname
    write(lout,'(1x,"The number of partially occ orbitals '// &
          'different from nact read from molpro. '// &
          'Some active orbitals must be unoccupied.",/)')
    Mon%NAct = nact
    Mon%ISwitchAct = 1  ! change Mon%num0 and Mon%num1 in select_active
    Mon%IWarn = Mon%IWarn + 1
 endif

! Set INAct (also works for open-shells)
 Mon%INAct = Mon%XELE-Tmp+1.d-1
 NOccup = Mon%INAct + Mon%NAct
 Mon%SumOcc = Tmp + Mon%INAct

 Mon%Occ = 0
 do i=1,NOccup
    if(i<=Mon%INAct) then
       Mon%Occ(i) = 1.d0
    else
       Mon%Occ(i) = EVal(i-Mon%INAct)
    endif
 enddo

 if(Mon%INAct==0) then
    Mon%NGem = 2

    Mon%IGem(1:Mon%NAct+Mon%INAct) = 1
    Mon%IGem(Mon%NAct+Mon%INAct+1:NBasis) = 2
 else
    Mon%NGem = 3
    Mon%IGem(1:Mon%INAct) = 1
    Mon%IGem(Mon%INAct+1:Mon%INAct+Mon%NAct) = 2
    Mon%IGem(Mon%INAct+Mon%NAct+1:NBasis) = 3
 endif

! construct CICoef
 do i=1,NBasis
    Mon%CICoef(i)=sqrt(Mon%Occ(i))
    if(Mon%Occ(i).lt.0.5d0) Mon%CICoef(i)=-Mon%CICoef(i)
 enddo

! call print_sqmat(OrbAux,NBasis)

 deallocate(EVal,work)

end subroutine readocc_molpro

subroutine readocc_cas_siri(mon,nbas,noSiriusRst)
!
! From SIRIFC (SIRIUS InterFaCe)
! a) read no. of active inactive orbs for SAPT-DALTON
!    total: NAct and INAct
!    in a given symmetry: INActS(1:NSym), NActS(1:NSym)
! b) read active 1-RDM (DVX)
!
! From SIRIUST.RST
! a') read occupation numbers
!
implicit none

type(SystemBlock)   :: mon
integer,intent(in)  :: nbas
logical,intent(out) :: noSiriusRst

logical           :: ioccsir,exsiri
integer           :: i,iunit,ios
integer           :: isym,off_i,off_a,off_x
integer           :: JACT,JORB,JOFF
integer           :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
                     NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,       &
                     NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBASM(8)

double precision             :: sum1,sum2
double precision,allocatable :: OccX(:)
integer :: MMASHX
double precision,allocatable :: DVX(:)
character(:),allocatable     :: sirfile,sirifcfile

 ! set filnames
 if(Mon%Monomer==1) then
    sirfile  = 'SIRIUS_A.RST'
    sirifcfile = 'SIRIFC_A'
 elseif(Mon%Monomer==2) then
    sirfile  = 'SIRIUS_B.RST'
    sirifcfile = 'SIRIFC_B'
 endif

 inquire(file=sirifcfile,EXIST=exsiri)
 if(exsiri) then
    open(newunit=iunit,file=sirifcfile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(iunit,'TRCCINT ')

    rewind(iunit)
    read (iunit)
    read (iunit)
    read (iunit) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
                 NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
                 NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBASM

    read (iunit)
    read (iunit)

    ! DV = 1-RDM in Dalton
    MMASHX = MAX(4,NNASHX)
    allocate (DVX(MMASHX))
    read (iunit) DVX(1:MMASHX)  ! 1-rdm, triang
    close(iunit)

    mon%INAct = nisht
    mon%NAct  = nasht

    !print*, 'DV (1-rdm Dalton)'
    !do i=1,mmashx
    !   print*, i ,dvx(i)
    !enddo

    if(NSym/=mon%NSym) stop "NSym from SIRIFC and AOONEINT do not match!"

    mon%INActS(1:mon%NSym) = NISH(1:NSym)
    mon%NActS(1:mon%NSym)  = NASH(1:NSym)

    if(nbast.ne.nbas) then
      write(LOUT,'(1x,a)') 'WARNING! NBasis FROM SIRIFC DOES NOT MATCH!'
      write(LOUT,'(1x,a,i5,1x,a,i5)') 'NBasis: ',nbas, 'SIRIFC: ', nbast
      write(LOUT,'()')
      mon%IWarn = mon%IWarn + 1
    endif

 else
    write(lout,'(1x,a)') 'SIRIFC not available!'
    stop
 endif

 !print*, 'NCONF =', NCONF
 !print*, 'INACT =', mon%INACt
 !print*, ' NACT =', mon%NAct
 !print*, 'NISH  =', NISH(1:NSym)
 !print*, 'NASH  =', NASH(1:NSym)

 ! CASSCF
 allocate(OccX(1:norbt))

 if (NCONF.eq.1 .and. NASHT.gt.1) then ! get occupations from 1-rdm
    ! for a single configuration (NCONF=1), e.g., CAS(5,3) for F2 (MS=1/2),
    ! Dalton does not store NATOCC in SIRIUS.RST
    ! In this case, use DV to get occupations
    OccX = 0d0
    JACT = 0
    JOFF = 0
    do ISYM=1,NSym
       JORB = 0
       do I=1,NISH(ISYM)
          JORB = JORB + 1
          OccX(JOFF+JORB) = 2.0D0
       enddo
       do I=1,NASH(ISYM) ! assume natural orbitals
          JORB = JORB + 1
          JACT = JACT + 1
          OccX(JOFF+JORB) = DVX((JACT*JACT+JACT)/2)
       enddo
       JOFF = JOFF + NORB(ISYM)
    enddo
    write(LOUT,'(1x,a,i2,a)') 'Occupancies for monomer',mon%Monomer,' read from 1-RDM (SIRIFC)'

 elseif (NCONF.gt.1) then
    ! use NATOCC label from SIRIUS.RST
    ! to read occupation numbers
    inquire(file=sirfile,EXIST=ioccsir)
    if(ioccsir) then

       noSiriusRst=.false.

       if (mon%Nact.ge.2) then

         open(newunit=iunit,file=sirfile,status='OLD', &
              access='SEQUENTIAL',form='UNFORMATTED')
         call readlabel(iunit,'NATOCC  ')
         read(iunit) OccX(1:NORBT)
         close(iunit)

       elseif(mon%NAct.le.1) then

         write(lout,'(/1x,a,i3)') 'Warning! Number of active orbitals = ',mon%NAct
         write(lout,'(1x,a)') 'Assuming a Hartree-Fock calculation...'
         OccX(1:mon%INAct) = 2d0
         OccX(mon%INAct+1:mon%INAct+mon%NAct) = 1d0

       endif
       write(LOUT,'(1x,a,i2,a)') 'Occupancies for monomer',mon%Monomer,' read from '// sirfile

    else

       noSiriusRst=.true.
       write(lout,'(1x,a)') 'SIRIUS.RST not available!'
       return

    endif

  endif ! NCONF test for OccX

  !Print*, 'Occupations symmetry-packed:'
  !do i=1,norbt
  !  print*, i, occX(i)
  !enddo

 allocate(mon%Occ(nbas))
 ! save occupations in mon%Occ
 ! order from sym ordering to inact-act (ISW/ISX in Dalton)
 mon%Occ = 0d0
 off_i = 0
 off_a = NISHT
 off_x = 0
 do isym=1,NSym
    mon%Occ(off_i+1:off_i+NISH(isym)) = OccX(off_x+1:off_x+NISH(isym))
    mon%Occ(off_a+1:off_a+NASH(isym)) = OccX(off_x+NISH(isym)+1:off_x+NISH(isym)+NASH(isym))
    off_i = off_i + NISH(isym)
    off_a = off_a + NASH(isym)
    off_x = off_x + NORB(isym)
 enddo
 deallocate(OccX)

 sum1 = 0d0
 do i=1,mon%INAct+mon%NAct
     mon%Occ(i) = mon%Occ(i)/2d0
     sum1 = sum1 + mon%Occ(i)
 enddo
 mon%SumOcc = sum1

end subroutine readocc_cas_siri

subroutine readocc_cas_occu(mon,nbas,noOccu)
!
! From occupations.dat
! a) read total number of active (NAct)
!    and inactive (INAct) orbitals for SAPT-DALTON,
!    in a given symmetry: INActS(1:NSym), NActS(1:NSym)
! b) read occupation numbers
!
implicit none

type(SystemBlock)  :: mon
integer,intent(in) :: nbas
logical,intent(out):: noOccu

integer                      :: i
integer                      :: iunit,ios
logical                      :: iocc
double precision             :: sum1
character(:),allocatable     :: occfile

 ! set filenames
 if(Mon%Monomer==1) then
    occfile='occupations_A.dat'
 elseif(Mon%Monomer==2) then
    occfile='occupations_B.dat'
 endif

 allocate(mon%Occ(nbas))
 print*, 'here2?'
 inquire(file=occfile,EXIST=iocc)
 if(iocc) then

    noOccu     = .false.
    mon%Occ    = 0d0
    mon%INActS = 0
    mon%NActS  = 0
    open(newunit=iunit,file=occfile,form='FORMATTED',status='OLD')

    ! read inactive,active,occupations
    read(iunit,*) mon%INAct, mon%NAct
    mon%INAct = mon%INAct/2
    read(iunit,*) (mon%Occ(i),i=1,mon%INAct+mon%NAct)

    sum1 = 0d0
    do i=1,mon%INAct+mon%NAct
       mon%Occ(i) = mon%Occ(i)/2d0
       sum1 = sum1 + mon%Occ(i)
    enddo
    mon%SumOcc = sum1

    ! active and inactive orbs in each symmetry
    read(iunit,*,iostat=ios) (mon%NActS(i),i=1,mon%NSym)
    if(ios==0) then
       read(iunit,*) (mon%INActS(i),i=1,mon%NSym)
    endif

    if(mon%NSym.gt.1) then
      call sort_sym_occ(nbas,mon%NSym,mon%INAct,mon%NAct,mon%Occ)
    endif

    write(LOUT,'(1x,a,i2,a)') 'Occupancies for monomer',mon%Monomer,' read from '// occfile

 else

    noOccu = .true.
    write(lout,'(1x,a)') 'occupations.dat not available!'

 endif

end subroutine readocc_cas_occu

subroutine readocc_hf_siri(mon,nbasis)
implicit none

type(SystemBlock)  :: mon
integer,intent(in) :: nbasis

integer                      :: MMORBT
integer                      :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
                                NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
                                NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBAS(8)
integer                      :: i,iunit,idx,irep,offset
double precision,allocatable :: fock(:)
character(:),allocatable     :: sirifile
logical                      :: exsiri

 ! set filnames
 if(Mon%Monomer==1) then
    sirifile  = 'SIRIFC_A'
 elseif(Mon%Monomer==2) then
    sirifile  = 'SIRIFC_B'
 endif

 inquire(file=sirifile,EXIST=exsiri)
 if(exsiri) then
    open(newunit=iunit,file=sirifile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(iunit,'TRCCINT ')

    rewind(iunit)
    read (iunit)
    read (iunit)
    read (iunit) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
                 NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,      &
                 NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBAS
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)
    read(iunit)

    mon%INAct = nisht
    mon%NAct  = nasht

    !print*, 'nisht,nasht',nisht,nasht
    !print*, 'nisht(1:ns)',nish(1:NSym)
    !print*, 'nasht(1:ns)',nash(1:NSym)

    mon%INActS(1:mon%NSym) = NISH(1:NSym)
    mon%NActS(1:mon%NSym)  = NASH(1:NSym)

    MMORBT = max(4,NNORBT)
    allocate(fock(MMORBT),mon%OrbE(NORBT))

    read(iunit) fock
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
!    print*, mon%OrbE

    deallocate(fock)

 endif

 close(iunit)

 ! Hartree-Fock case
 allocate(mon%Occ(nbasis))

 mon%Occ = 0d0
 mon%Occ(1:mon%NAct+mon%INAct) = 1d0

 !print*, 'Occ:',mon%Occ(1:nisht+nasht)

end subroutine readocc_hf_siri

subroutine read2rdm(Mon,NBas)
!
! Purpose: load rdm2.dat file to memory
!          as Mon%RDM2(NRDM2Act) matrix
!
implicit none

type(SystemBlock) :: Mon
integer, intent(in) :: NBas
character(:),allocatable :: rdmfile
integer :: iunit,ios
integer :: NRDM2Act
integer :: Ind1(NBas),Ind2(NBas)
integer :: i,j,k,l
double precision :: val
double precision,parameter :: Half=0.5d0
integer,external :: NAddrRDM

 if(Mon%Monomer==1) then
    rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    rdmfile='rdm2_B.dat'
 endif

 Ind1=0
 Ind2=0
 do i=1,Mon%NAct
    Ind1(i) = Mon%INAct + i
    Ind2(Mon%INAct+i) = i
 enddo

 NRDM2Act = Mon%NAct**2*(Mon%NAct**2+1)/2

 if(allocated(Mon%RDM2))    deallocate(Mon%RDM2)

 allocate(Mon%RDM2(NRDM2Act))
 Mon%RDM2(1:NRDM2Act)=0

 open(newunit=iunit,file=rdmfile,status='OLD',&
      form='FORMATTED')
 do

   read(iunit,'(4i4,f19.12)',iostat=ios) i,j,k,l,val

!  val IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)

   if(ios==0) then
      Mon%RDM2(NAddrRDM(j,l,i,k,Mon%NAct))=Half*val

    elseif(ios/=0) then
       exit

    endif

 enddo
 close(iunit)

 print*, 'read2rdm: MON%RDM2',norm2(Mon%RDM2)

 if(allocated(Mon%Ind2)) deallocate(Mon%Ind2)
 allocate(Mon%Ind2(NBas))

 Mon%Ind2 = Ind2

end subroutine read2rdm

subroutine read2rdm_spin(Mon,NBas)
!
! Purpose: a) load rdms201.dat file to memory
!          as Mon%RDM201(NRDM2Act) matrix
!          B) if Hartree-Fock, assume = 0
implicit none

type(SystemBlock)   :: Mon
integer, intent(in) :: NBas

character(:),allocatable :: rdmfile
integer :: iunit,ios
integer :: NRDM2Act
integer :: i,j,k,l
double precision :: val
integer,external :: NAddrRDM

if (Mon%Monomer==1) then
   rdmfile='rdms201_A.dat'
elseif (Mon%Monomer==2) then
   rdmfile='rdms201_B.dat'
endif

if(allocated(Mon%RDM201)) deallocate(Mon%RDM201)

NRDM2Act = Mon%NAct**2*(Mon%NAct**2+1)/2

allocate(Mon%RDM201(NRDM2Act))
Mon%RDM201(1:NRDM2Act) = 0

open(newunit=iunit,file=rdmfile,status='OLD',&
     form='FORMATTED')

do
  read(iunit,'(4i4,f19.12)',iostat=ios) i,j,k,l,val

!  val IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)
!  RDM201 = \Gamma^++++ - \Gamma^---- + \Gamma^-+-+ - \Gamma^+-+-


  if(ios==0) then
     Mon%RDM201(NAddrRDM(j,l,i,k,Mon%NAct)) = 0.5d0*val
  elseif(ios/=0) then
     exit
  endif
enddo

close(iunit)

end subroutine read2rdm_spin

subroutine prepare_no_molpro_skip(CMONOAct,CAOMO,INAct,NAct,NAO,NBasis)
implicit none
!
! Purpose: get AO-->NO transformation
!
! CMONO[in] :: on input MOtoNO
! CAOMO[in] :: on input AOtoMO
!     [out] :: on output AOtoNO
!
integer,intent(in) :: NAO,NBasis
integer,intent(in) :: INAct,NAct
double precision   :: CMONOAct(NBasis,NBasis),CAOMO(NAO,NBasis)

integer :: i,j
double precision   :: CMONO(NBasis,NBasis)
double precision   :: work(NAO,NBasis)

! skip canonicalization

 CMONO = 0d0
 forall(i=1:NBasis) CMONO(i,i)=1d0
 do i=1,NAct
    do j=1,NAct
       CMONO(INAct+i,INAct+j) = CMONOAct(i,j)
    enddo
 enddo

 call dgemm('N','N',NAO,NBasis,NBasis,1d0,CAOMO,NAO,CMONO,NBasis,0d0,work,NAO)
 CAOMO = work

end subroutine prepare_no_molpro_skip

subroutine prepare_no_molpro(OrbCAS,OneRdm,CMONOAct,Mon,AOBasis,System, &
                      CholeskyVecs,CholeskyVecsOTF,  &
                      Xgp,Zgk,NGridTHC,NCholeskyTHC, &
                      Flags,NBasis)
implicit none
!
! Prepare C(AO,NO) orbitals by diagonalization of inactive and virtual
! blocks of the Fock matrix
!
! CMONOAct[in]  :: on input  C(MO,NO) in active MOs
!                  (from diagonalization of 1-RDM in MOs)
! OneRDM[inout] :: on input  1-RDM in AO
!                  on output 1-RDM in MO
! OrbCAS[inout] :: on input  C(SAO,MO) from Molpro files
!                  on output C(SAO,NO)
!
! For CholeskyOTF: compute J and K matrices
!
type(FlagsData)        :: Flags
type(SystemBlock)      :: Mon
type(TCholeskyVecs)    :: CholeskyVecs
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(TAOBasis)         :: AOBasis
type(TSystem)          :: System

integer,intent(in) :: NBasis
integer,intent(in) :: NGridTHC,NCholeskyTHC
double precision,intent(in)    :: CMONOAct(NBasis,NBasis)
double precision,intent(in)    :: Xgp(NGridTHC,NBasis),Zgk(NGridTHC,NCholeskyTHC)
double precision,intent(inout) :: OrbCAS(NBasis,NBasis)
double precision,intent(inout) :: OneRdm(NBasis*(NBasis+1)/2)

integer :: NOccup,NVirt,NSym
integer :: NCholesky
integer :: i,j,ia,ib,iab,ioff,idx,NInte1
integer :: itsoao(NBasis),jtsoao(NBasis)

double precision :: CAOMO(NBasis,NBasis),CSAOMO(NBasis,NBasis), &
                    CAONO(NBasis,NBasis),CMONO(NBasis,NBasis)
double precision :: FockSq(NBasis,NBasis),SAO(NBasis,NBasis)
double precision :: work(NBasis,NBasis),SC(NBasis,NBasis)

double precision,allocatable :: H0(:), GammaF(:),Fock(:)
double precision,allocatable :: work1(:),work2(:),work3(:)

character(:),allocatable :: onefile,rdmfile,aoerfile
character(:),allocatable :: orbaofile
! testy
integer :: info

 NInte1 = NBasis*(NBasis+1)/2
 NOccup = Mon%INAct + Mon%NAct
 NVirt = NBasis - Mon%INAct - Mon%NAct

 if(Mon%Monomer==1) then
   onefile   = 'AOONEINT_A'
   rdmfile   = '2RDMA'
   aoerfile  = 'AOERFSORT'
   orbaofile = 'MOLPRO_A.MOPUN'
 elseif(Mon%Monomer==2) then
   onefile   = 'AOONEINT_B'
   rdmfile   = '2RDMB'
   orbaofile = 'MOLPRO_B.MOPUN'
   if(Mon%SameOm) then
      aoerfile = 'AOERFSORT'
   else
      aoerfile = 'AOERFSORTB'
   endif
 endif

 allocate(Mon%NumOSym(15),Mon%IndInt(NBasis))
 allocate(work1(NInte1),work2(NInte1),work3(NBasis))
 allocate(H0(NInte1),GammaF(NInte1),Fock(NBasis**2))

 call create_ind_molpro(rdmfile,Mon%NumOSym,Mon%IndInt,NSym,NBasis)

! COPY C(MO,NO)Act TO CMONO AND OFF SET BY NInAc
 CMONO = 0
 forall(i=1:NBasis) CMONO(i,i)=1d0
 ! with Diag8:
 do i=1,Mon%NAct
    do j=1,Mon%NAct
       CMONO(Mon%INAct+i,Mon%INAct+j) = CMONOAct(i,j)
    enddo
 enddo
 !print*, 'prepare_no: CMONO'
 !do j=1,NBasis
 !   write(6,'(*(f12.6))') (CMONO(i,j),i=1,NBasis)
 !enddo
 ! with dsyev
 !do i=1,Mon%NAct
 !   do j=1,Mon%NAct
 !      URe(Mon%INAct+i,Mon%INAct+j) = OrbAux(NBasis+1-i,NBasis+1-j)
 !   enddo
 !enddo
!print*, norm2(URe)
! call print_sqmat(URe,NBasis)

! FIND CANONICAL INACTIVE AND VIRTUAL ORBITALS

 GammaF = 0
 idx = 0
 do j=1,Mon%INAct
    do i=1,j
       idx = idx + 1
       if(i==j) GammaF(idx) = 1.0d0
    enddo
 enddo
 idx = 0
 do j=1,Mon%NAct
    do i=1,j
       idx = idx + 1
       ioff = (Mon%INAct+j)*(Mon%INAct+j-1)/2 + Mon%INAct
       GammaF(ioff+i) = OneRdm(idx)
    enddo
 enddo

! reorder MOs to no symmetry
! (in Molpro they are arranged by irreps)
do i=1,NBasis
   do j=1,NBasis
      CSAOMO(Mon%IndInt(i),j) = OrbCAS(j,i)
   enddo
enddo

!print*, 'IndInt = '
!do i=1,NBasis
!  print*,i,mon%IndInt(i)
!enddo
!Print*, 'CSAOMO-prepare no =',norm2(CSAOMO)
!do i=1,NBasis
!   write(6,'(*(f13.8))') (CSAOMO(i,j),j=1,NBasis)
!enddo

 iab = 0
 do ia=1,NBasis
    do ib=1,ia
       iab = iab + 1
       OneRdm(iab) = 0.d0
       do i=1,NBasis
          do j=1,NBasis
             idx = max(i,j)*(max(i,j)-1)/2+min(i,j)
             OneRdm(iab) = OneRdm(iab) + CSAOMO(i,ia)*CSAOMO(j,ib)*GammaF(idx)
          enddo
       enddo
    enddo
 enddo

 ! create Fock matrix
 ! H0 = XOne
 call readoneint_molpro(H0,onefile,'ONEHAMIL',.true.,NInte1)
 ! work2 = Fock
 if(Flags%IFunSR==0) then
 ! CASSCF,Hartree-Fock

   if(Flags%ICholeskyBIN==0.and.Flags%ICholeskyOTF==0) then

     call FockGen_mithap(work2,OneRdm,H0,NInte1,NBasis,'AOTWOSORT')

   elseif(Flags%ICholeskyBIN==1) then

     NCholesky = CholeskyVecs%NCholesky
     call FockGen_CholR(work2,CholeskyVecs%R(1:NCholesky,1:NInte1),OneRdm,H0, &
                        NInte1,NCholesky,NBasis)

   elseif(Flags%ICholeskyOTF==1) then

     ! for SAPT with Cholesky OTF compute J and K matrices here
     ! used later in eletrostatic potential (V+J in calc_elpot)
     ! and in E1exch (K for monomer B)
     allocate(Mon%Jmat(NBasis,NBasis),Mon%Kmat(NBasis,NBasis))

     CSAOMO = transpose(CSAOMO)

     call read_caomo_molpro(CAOMO,SAO,itsoao,jtsoao,orbaofile,'CASORBAO',NBasis)

     !call CholeskyOTF_Fock_MO_v1(FockSq,CholeskyVecsOTF,&
     !                      AOBasis,System,mon%Monomer, &
     !                      CAOMO,CSAOMO,H0,GammaF, &
     !                      Flags%MemType,Flags%MemVal,NInte1,NBasis, &
     !                      Mon%Jmat,Mon%Kmat)
     call CholeskyOTF_Fock_MO_v2(FockSq,CholeskyVecsOTF,&
                           AOBasis,System,mon%Monomer,'MOLPRO', &
                           CAOMO,CSAOMO,H0,GammaF, &
                           Xgp,Zgk,NGridTHC,NCholeskyTHC, &
                           Flags%MemType,Flags%MemVal,NInte1,NBasis, &
                           Flags%IH0test, &
                           Mon%Jmat,Mon%Kmat)

     call sq_to_triang2(FockSq,work2,NBasis)

     CSAOMO = transpose(CSAOMO)

   endif

 elseif(Flags%IFunSR>0) then
 ! Kohn-Sham

   ! add and store Coulomb
   ! for RSH short-range Coulomb is stored
   allocate(Mon%VCoul(NInte1))
   call PotCoul_mithap(Mon%VCoul,OneRdm,Mon%doRSH,aoerfile,NBasis)
   ! RSH
   if(Mon%doRSH) then
     ! generate long-range Fock
     call FockGen_mithap(work2,OneRdm,H0,NInte1,NBasis,aoerfile)
     work2 = work2 + Mon%VCoul
   else
   ! non-hybrid DFAs
   !  work2 = H0
     work2 = H0 + Mon%VCoul
   endif

 endif
 if(Flags%ICholeskyOTF==0) call tran_matTr(work2,CSAOMO,CSAOMO,NBasis,.false.)

 Fock = 0
 work3 = 0
 allocate(Mon%OrbE(NBasis))
!INACTIVE
 if(Mon%INAct/=0) then
    do i=1,Mon%INAct
       do j=1,Mon%INAct
          idx = max(i,j)*(max(i,j)-1)/2+min(i,j)
          Fock((j-1)*Mon%INAct+i) = work2(idx)
       enddo
    enddo
    call Diag8(Fock,Mon%INAct,Mon%INAct,work3,work1)
    !call dsyev('V','U',Mon%INAct,Fock,Mon%INAct,work3,work1,3*Mon%INAct,info)
    !print*, 'INACTIVE:',work3(1:Mon%INAct)
    ! test for ICPHF
    Mon%OrbE(1:Mon%INAct) = work3(1:Mon%INAct)

    do i=1,Mon%INAct
      do j=1,Mon%INAct
         CMONO(i,j) = Fock((j-1)*Mon%INAct+i)
      enddo
    enddo
 endif

! VIRTUAL
 if(NVirt/=0) then
    do i=1,NVirt
       do j=1,NVirt
          idx = (max(i+NOccup,j+NOccup)*(max(i+NOccup,j+NOccup)-1))/2 &
              + min(i+NOccup,j+NOccup)
          Fock((j-1)*NVirt+i) = work2(idx)
       enddo
    enddo
    call Diag8(Fock,NVirt,NVirt,work3,work1)
    !call dsyev('V','U',NVirt,Fock,NVirt,work3,work1,3*NVirt,info)
    do i=1,NVirt
       do j=1,NVirt
          CMONO(i+NOccup,j+NOccup) = Fock((j-1)*NVirt+i)
       enddo
    enddo
 endif
 ! test for ICPHF
 !print*, 'work3',work3
 Mon%OrbE(Mon%INAct+1:NBasis)=work3(1:NBasis-Mon%INAct)

! END OF CANONICALIZING

! transform orbitals to (SAO,NO)
! CMONO = C(NO,MO); CSAOMO = C(MO,SAO)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,CMONO,NBasis,CSAOMO,NBasis,0d0,OrbCAS,NBasis)
OrbCAS = transpose(OrbCAS)

if(Flags%ICholeskyOTF==1) then

   allocate(Mon%CAONO(NBasis,NBasis))
   ! CAOMO = C(AO,MO) ; CMONO = C(NO,MO)
   call dgemm('N','T',NBasis,NBasis,NBasis,1d0,CAOMO,NBasis,CMONO,NBasis,0d0,CAONO,NBasis)
   Mon%CAONO = CAONO

   ! transform J/K from MO to AO with SAO
   ! remember: C^T(AO,MO).S(AO).C(AO,MO) = 1
   !   so that C^-1(AO,MO) = C^T.S(AO)
   !           J_MO = C^T . J_AO . C
   !           J_AO = SC . J_MO . (SC)^T
   print*, 'Jmat-MO',norm2(mon%Jmat)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SAO,NBasis,CAOMO,NBasis,0d0,SC,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,mon%Jmat,NBasis,0d0,work,NBasis)
   call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work,NBasis,SC,NBasis,0d0,mon%Jmat,NBasis)

   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,mon%Kmat,NBasis,0d0,work,NBasis)
   call dgemm('N','T',NBasis,NBasis,NBasis,-1d0,work,NBasis,SC,NBasis,0d0,mon%Kmat,NBasis)

   !write(6,*) 'Kmat-AO OTF',mon%Monomer,norm2(mon%Kmat)
   !do j=1,NBasis
   !   write(LOUT,'(*(f13.8))') (Mon%Kmat(i,j),i=1,NBasis)
   !enddo
   !write(LOUT,'()')
   !write(6,*) 'Jmat-AO OTF',mon%Monomer,norm2(mon%Jmat)
   !do j=1,NBasis
   !   write(LOUT,'(*(f13.8))') (Mon%Jmat(i,j),i=1,NBasis)
   !enddo
   !write(LOUT,'()')

endif

deallocate(work3,work2,work1)
deallocate(H0,Fock)
deallocate(Mon%IndInt)

end subroutine prepare_no_molpro

subroutine prepare_rdm2_molpro(Mon,OrbAux,NBasis)
implicit none

type(SystemBlock) :: Mon

integer,intent(in) :: NBasis
double precision,intent(in) :: OrbAux(NBasis,NBasis)
integer :: i,j,k,l,ij,kl,iunit,NRDM2Act
double precision,allocatable :: RDM2Act(:),work1(:)
character(:),allocatable :: rdmfile,outfile
integer,external :: NAddrRDM

 if(Mon%Monomer==1) then
   rdmfile='2RDMA'
   outfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
   rdmfile='2RDMB'
   outfile='rdm2_B.dat'
 endif

 NRDM2Act = Mon%NAct**2*(Mon%NAct**2+1)/2
 allocate(RDM2Act(NRDM2Act),work1(Mon%NAct**2))
 RDM2Act = 0
 call read_2rdm_molpro(RDM2Act,Mon%InSt(1,1),Mon%InSt(2,1),&
                       Mon%ISpinMs2,rdmfile,Mon%IWarn,Mon%NAct)

 if (Mon%NatOrb==1) then
    do i=1,Mon%NAct
       do j=1,Mon%NAct
          work1((j-1)*Mon%NAct+i) = OrbAux(Mon%INAct+j,Mon%INAct+i)
       enddo
    enddo
 else
    do i=1,Mon%NAct
       do j=1,Mon%NAct
          work1((j-1)*Mon%NAct+i) = OrbAux(i,j)
       enddo
    enddo
 endif
 call TrRDM2(RDM2Act,work1,Mon%NAct,NRDM2Act)

 open(newunit=iunit,file=outfile,status='replace',&
      form='formatted')
 do i=1,Mon%NAct
   do j=1,Mon%NAct
      ij = (i-1)*Mon%NAct+j
      do k=1,Mon%NAct
         do l=1,Mon%NAct
            kl = (k-1)*Mon%NAct+l
            if(ij>=kl) then
              write(iunit,'(4i4,f19.12)') &
                k,i,l,j,2d0*RDM2Act(NAddrRDM(i,j,k,l,Mon%NAct))
            endif
         enddo
      enddo
   enddo
 enddo

 close(iunit)

 deallocate(work1,RDM2Act)

end subroutine prepare_rdm2_molpro

subroutine prepare_rdm2_approx(Mon,IRDM2Typ,NBasis)
!
! replace rdm2_A.dat and rdm2_B.dat files
! with approximate density matrices: DMFT or noncumulant
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: IRDM2Typ,NBasis

integer :: i,j,k,l,ij,kl
integer :: NOccup
integer :: iunit,NRDM2Act
double precision :: xnorm
double precision,allocatable :: RDM2val(:,:,:,:)
character(:),allocatable :: rdmfile
integer,external :: NAddrRDM

if(Mon%Monomer==1) then
  rdmfile='rdm2_A.dat'
elseif(Mon%Monomer==2) then
  rdmfile='rdm2_B.dat'
endif

print*, 'REPLACE 2-RDM with APPROXIMATE FORM in ERPA!'

NOccup = Mon%INAct+Mon%NAct
print*, 'NOccup',NOccup
print*, 'FLAG',IRDM2TYP

allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

! Gamma(prqs) = 2*np*nq \delta_pr \delta_qs - F_pq \delta_ps \delta_qr
!       1122
RDM2val = 0d0
! Coulomb (nc part)
do i=1,NOccup
   do j=1,NOccup
      RDM2val(i,i,j,j) = RDM2val(i,i,j,j) + 2d0*Mon%Occ(i)*Mon%Occ(j)
   enddo
enddo
if(IRdm2Typ==0) then
   ! exchange
   do i=1,NOccup
      do j=1,NOccup
         RDM2val(i,j,j,i) = RDM2val(i,j,j,i) - Mon%Occ(i)*Mon%Occ(j)
      enddo
   enddo
elseif(IRdm2Typ==1.or.IRDM2Typ==11) then
   ! exchange-corr
   do i=1,NOccup
      do j=1,NOccup
         RDM2val(i,j,j,i) = RDM2val(i,j,j,i) - sqrt(Mon%Occ(i)*Mon%Occ(j))
      enddo
   enddo
endif
!print*, 'RDM2val =',norm2(RDM2val)

xnorm = 0d0
do i=1,NOccup
   do j=1,NOccup
      xnorm = xnorm + RDM2val(i,i,j,j)
   enddo
enddo

if(mon%monomer==1) write(lout,'(/1x,a)') 'Monomer A'
if(mon%monomer==2) write(lout,'(/1x,a)') 'Monomer B'
write(lout,'(1x,a,f12.6)',advance="no") '2-RDM2 norm = ', xnorm
write(lout,'(1x,a,f8.3,a)') '(reference =', Mon%XELE*(2d0*Mon%XELE-1), ')'

!! re-normalize 2RDM
!print*, 're-normalize 2-RDM...'
!RDM2val = RDM2val * Mon%XELE*(2d0*Mon%XELE-1) / xnorm

open(newunit=iunit,file=rdmfile,status='replace',&
     form='formatted')
do i=1,Mon%NAct
  do j=1,Mon%NAct
     ij = (i-1)*Mon%NAct+j
     do k=1,Mon%NAct
        do l=1,Mon%NAct
           kl = (k-1)*Mon%NAct+l
           if(ij>=kl) then
             write(iunit,'(4i4,f19.12)') &
               !k,i,l,j,RDM2val(Mon%INAct+i,Mon%INAct+j,Mon%INAct+k,Mon%INAct+l)
               k,i,l,j,2d0*RDM2val(Mon%INAct+k,Mon%INAct+i,Mon%INAct+j,Mon%INAct+l)
           endif
        enddo
     enddo
  enddo
enddo

deallocate(RDM2val)

close(iunit)

end subroutine prepare_rdm2_approx

! === GAMMCOR-specific TREXIO subroutines (not in pr-dmft) ===

subroutine onel_trexio(NBasis,NAO,Mon,SAPT)
!
! Purpose:
! reads NAO
! reads Smat, Hmat, Tmat in AOs
! reads geometry: charg,coord
!
 implicit none

 type(SaptData)      :: SAPT
 type(SystemBlock)   :: Mon
 integer,intent(in)  :: NBasis
 integer,intent(out) :: NAO

 integer    :: rc
 integer    :: i,num,offset
 integer(8) :: f
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
 double precision, allocatable :: kinetic(:)
 double precision, allocatable :: charge(:),coord(:,:)
 character(:),allocatable      :: outfile

if(Mon%Monomer==1) then
  outfile = 'ONEEL_A'
elseif(Mon%Monomer==2) then
  outfile = 'ONEEL_B'
endif

f = trexio_open (Mon%TrexFile, 'r', TREXIO_HDF5, rc)

rc = trexio_has_ao_num(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No AO num in file'
end if
rc = trexio_read_ao_num(f,NAO)

allocate(Hmat(NAO**2),Vmat(NAO**2),Smat(NAO**2))
allocate(kinetic(NAO**2))

rc = trexio_has_ao_1e_int_overlap(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No overlap in file'
end if
rc = trexio_read_ao_1e_int_overlap(f, Smat)

rc = trexio_has_ao_1e_int_kinetic(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No kinetic in file'
end if
rc = trexio_read_ao_1e_int_kinetic(f, kinetic)

rc = trexio_has_ao_1e_int_potential_n_e(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No potential_n_e in file'
end if
rc = trexio_read_ao_1e_int_potential_n_e(f, Vmat)

Hmat = 0
Hmat = Vmat + kinetic

call writeoneint(outfile,NAO**2,Smat,Vmat,Hmat)

SAPT%NAO = NAO

rc = trexio_has_nucleus_charge(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No nucleus charge in file'
end if
rc = trexio_has_nucleus_coord(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No nuclei coords in file'
end if

rc = trexio_read_nucleus_num(f,num)
allocate(charge(num),coord(3,num))
rc = trexio_read_nucleus_charge(f,charge)
rc = trexio_read_nucleus_coord(f,coord)

offset = 0
if(Mon%Monomer==2) offset = SAPT%monA%NCen
do i=1,Mon%NCen
   Mon%charg(i)   = charge(offset+i)
   Mon%xyz(i,1:3) = coord(1:3,offset+i)
enddo

rc = trexio_close(f)

deallocate(coord,charge)
deallocate(kinetic)
deallocate(Smat,Vmat,Hmat)

end subroutine onel_trexio

subroutine arrange_oneint(mat,nbas,SAPT)
implicit none

type(SaptData)     :: SAPT
integer,intent(in) :: nbas
double precision,intent(inout) :: mat(nbas,nbas)

!call read_syminf(SAPT%monA,SAPT%monB,nbas)

if(SAPT%monB%switchAB) then
   call gen_swap_rows(mat,nbas,nbas,SAPT%monA%NSym,&
                      SAPT%monA%NMonBas,SAPT%monB%NMonBas)
   call gen_swap_cols(mat,nbas,nbas,SAPT%monA%NSym,&
                      SAPT%monA%NMonBas,SAPT%monB%NMonBas)
endif

!call swap_rows(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)
!call swap_cols(SAPT%monA%NMonOrb,SAPT%monB%NMonOrb,mat)

end subroutine arrange_oneint

subroutine sort_sym_occ(nbas,nsym,INAct,NAct,Occ)
implicit none

integer,intent(in) :: nbas, nsym, INAct, NAct
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
       if(ICpy2(i).eq.0.and.ICpy1(ii).eq.0.and.Occ(i).eq.1.0D0) then
       !if(ICpy2(i).eq.0.and.ICpy1(ii).eq.0.and.Occ(i).eq.2.0D0) then
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

subroutine readgvb(mon,n,cfile)
! set: NAct (number of act. geminals)
!      INAct, CICoef, Occ, IGem
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

subroutine writeoneint(mon,ndim,S,V,H)
implicit none

integer :: ione,ndim
character(*) :: mon
double precision,dimension(ndim) :: S, V, H

 open(newunit=ione,file=mon,form='unformatted')
 write(ione) 'OVERLAP ', S
 write(ione) 'POTENTAL', V
 write(ione) 'ONEHAMIL', H
! write(ione) 'KINETINT', K
 close(ione)

 write(LOUT,'(1x,a)') 'One-electron integrals written to file: '//mon

end subroutine writeoneint

end module sapt_interface_io
