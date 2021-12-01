module sapt_inter

use types
use timing
use tran
use sorter
!use Cholesky_old
use Cholesky
use abmat

implicit none

contains

subroutine sapt_interface(Flags,SAPT,NBasis)
!
! SAPT-DALTON requires SIRIFC and SIRIUS.RST
!                   or SIRIFC and occupations.dat
!
implicit none

type(FlagsData)     :: Flags
type(SaptData)      :: SAPT
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis

integer    :: NSq,NInte1,NInte2
integer    :: NCholesky
integer    :: dimOA,dimOB,dimVA,dimVB,nOVA,nOVB
integer    :: NCMOt, NOrbt, NBasist
integer    :: NSym, NBas(8)
integer    :: NOcc(8),NOrbs(8)
integer    :: ione,iorb,isiri,i,j,ij
integer    :: p,q,pq
integer(8) :: MemSrtSize
integer    :: ICholOLD
double precision :: tmp
double precision :: potnucA,potnucB
double precision :: potnuc,emy,eactiv,emcscf

double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: work(:,:)
double precision,allocatable :: Ha(:),Hb(:)
double precision,allocatable :: Va(:),Vb(:),S(:)
double precision,allocatable :: Ca(:),Cb(:)
double precision,allocatable :: AuxA(:,:),AuxB(:,:)
double precision,allocatable :: OneRdmA(:),OneRdmB(:)

logical :: doRSH
double precision,allocatable :: Sa(:,:),Sb(:,:)
double precision :: Tcpu,Twall

! set monomer print level
 SAPT%monA%IPrint = SAPT%IPrint
 SAPT%monB%IPrint = SAPT%IPrint

 print*, 'Using new sapt_interface...! '

! set dimensions
 NSq = NBasis**2
 NInte1 = NBasis*(NBasis+1)/2
 NInte2 = NInte1*(NInte1+1)/2
 SAPT%monA%NDim = NBasis*(NBasis-1)/2
 SAPT%monB%NDim = NBasis*(NBasis-1)/2

! set RSH
 SAPT%doRSH = .false.
 if(Flags%IFunSR==1.or.Flags%IFunSR==2) SAPT%doRSH = .true.
 doRSH = SAPT%doRSH

! set semicoupled dispersion
 if(Flags%ICASSCF==0) then
    SAPT%SemiCoupled = .false.
 endif

! read and dump 1-electron integrals
 if(SAPT%InterfaceType==1) then
    call onel_dalton(SAPT%monA%Monomer,NBasis,NSq,NInte1,SAPT%monA,SAPT)
    call onel_dalton(SAPT%monB%Monomer,NBasis,NSq,NInte1,SAPT%monB,SAPT)
 elseif(SAPT%InterfaceType==2) then
    call onel_molpro(SAPT%monA%Monomer,NBasis,NSq,NInte1,SAPT%monA,SAPT)
    call onel_molpro(SAPT%monB%Monomer,NBasis,NSq,NInte1,SAPT%monB,SAPT)
 endif

! add empty line
 write(lout,'()')

! read coefficient, occupancies
 if(SAPT%InterfaceType==1) then
    call readocc_dalton(NBasis,SAPT%monA,Flags)
    call readocc_dalton(NBasis,SAPT%monB,Flags)

 !! test readsymfino
 ! allocate(Sa(NBasis,NBasis),Sb(NBasis,NBasis))
 ! call get_one_mat('S',Sa,1,nbasis)
 ! call get_one_mat('S',Sb,2,nbasis)
 ! print*, 'Sb-Sa',norm2(Sb-Sa)
 ! deallocate(Sb,Sa)

 elseif(SAPT%InterfaceType==2) then
    allocate(AuxA(NBasis,NBasis),AuxB(NBasis,NBasis),&
             OneRdmA(NInte1),OneRdmB(NInte1))
    call readocc_molpro(NBasis,SAPT%monA,AuxA,OneRdmA,Flags)
    call readocc_molpro(NBasis,SAPT%monB,AuxB,OneRdmB,Flags)
 endif
 call print_occ(NBasis,SAPT,Flags%ICASSCF)

! read orbitals
! norb.leq.nbas, orbitals mays be deleted due to linear
! dependecies in large basis sets; ncmot = norb*nbas
 allocate(Ca(NBasis*NBasis),Cb(NBasis*NBasis))

 if(SAPT%InterfaceType==1) then

    call read_mo_dalton(Ca,NBasis,SAPT%monA%NSym,SAPT%monA%NSymBas,SAPT%monA%NSymOrb,&
                 'SIRIUS_A.RST','DALTON_A.MOPUN')
    call read_mo_dalton(Cb,NBasis,SAPT%monB%NSym,SAPT%monB%NSymBas,SAPT%monB%NSymOrb,&
                 'SIRIUS_B.RST','DALTON_B.MOPUN')
    call arrange_mo(Cb,NBasis,SAPT)

 elseif(SAPT%InterfaceType==2) then
    call read_mo_molpro(Ca,'MOLPRO_A.MOPUN','CASORB  ',NBasis)
    call read_mo_molpro(Cb,'MOLPRO_B.MOPUN','CASORB  ',NBasis)
 endif

! symmetry sorting
 if(SAPT%InterfaceType==1) then
    if(SAPT%monA%NSym.gt.1) then
       call sort_sym_mo(Ca,NBasis,SAPT%monA)
    endif
    if(SAPT%monB%NSym.gt.1) then
       call sort_sym_mo(Cb,NBasis,SAPT%monB)
    endif
 endif

! read dipole moments
 if(SAPT%InterfaceType==1.and.SAPT%ic6==1) then
    write(LOUT,*) 'DALTON CANNOT BE USED FOR Cn COEFFS!'
 elseif(SAPT%InterfaceType==2.and.SAPT%ic6==1) then
    !call read_dip_molpro(SAPT%monA,'DIP_A',NBasis)
    !call read_dip_molpro(SAPT%monB,'DIP_B',NBasis)
 endif

! read 2-el integrals
 call clock('START',Tcpu,Twall)

! memory allocation for sorter
 MemSrtSize = Flags%MemVal*1024_8**Flags%MemType

! for testing old Cholesky
 ICholOld = 0

 if(Flags%ICholesky==0.or.ICholOld==1) then
    if(SAPT%InterfaceType==1) then
       call readtwoint(NBasis,1,'AOTWOINT_A','AOTWOSORT',MemSrtSize)
    elseif(SAPT%InterfaceType==2) then
       call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT',MemSrtSize)
       if(doRSH) then
          if(SAPT%SameOm) then
             call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
          else
             call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
             call readtwoint(NBasis,2,'AOTWOINT.erfB','AOERFSORTB',MemSrtSize)
          endif
       endif
    endif
 endif

! Cholesky decomposition
 if(Flags%ICholesky==1) then

    !! old Cholesky
    !if(ICholOLD==1) print*, 'old Cholesky transformation...'
    !if(ICholOLD==1) call chol_CoulombMatrix(CholeskyVecs,'AOTWOSORT',Flags%ICholeskyAccu)

    ! new Cholesky
    if(ICholOLD==0) print*, 'new Cholesky transformation...'
    if(SAPT%InterfaceType==1) then
       call chol_CoulombMatrix(CholeskyVecs,NBasis,'AOTWOINT_A',1,Flags%ICholeskyAccu)
    elseif(SAPT%InterfaceType==2) then
       call chol_CoulombMatrix(CholeskyVecs,NBasis,'AOTWOINT.mol',2,Flags%ICholeskyAccu)
    endif

    SAPT%NCholesky  = CholeskyVecs%NCholesky
    SAPT%monA%NChol = SAPT%NCholesky
    SAPT%monB%NChol = SAPT%NCholesky

 endif
 call clock('2ints',Tcpu,Twall)

 if(SAPT%InterfaceType==2) then
    call prepare_no(OneRdmA,AuxA,Ca,SAPT%monA,CholeskyVecs,Flags%IFunSR,Flags%ICholesky,NBasis)
    call prepare_no(OneRdmB,AuxB,Cb,SAPT%monB,CholeskyVecs,Flags%IFunSR,Flags%ICholesky,NBasis)
    call prepare_rdm2_file(SAPT%monA,AuxA,NBasis)
    call prepare_rdm2_file(SAPT%monB,AuxB,NBasis)
 endif

! maybe better: add writing Ca, Cb to file?!
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

 if(Flags%ICholesky==1) then
    ! transform Cholesky Vecs to NO
    call chol_sapt_NOTransf(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call clock('chol_NOTransf',Tcpu,Twall)
 endif

! MAYBE: one should print with NOrbt?
 !if(SAPT%IPrint.ne.0) call print_mo(Ca,NBasis,'MONOMER A')
 !if(SAPT%IPrint.ne.0) call print_mo(Cb,NBasis,'MONOMER B')

! ABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB
 if(SAPT%IPrint.gt.100) call print_TwoInt(NBasis)

 call print_active(SAPT,NBasis)

 if(Flags%ISERPA==2) then
    ! set iPINO
    if(Flags%ICASSCF==1.and.Flags%ISHF==1.and.Flags%SaptLevel/=10) then
       ! 2-electron FCI
       SAPT%iPINO = 0
    elseif(Flags%ICASSCF==1.and.Flags%SaptLevel==10) then
       ! 2-electron e2dispCAS
       SAPT%iPINO = 1
    elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Flags%SaptLevel/=10) then
       ! CAS/LR
       SAPT%iPINO = 2
    elseif(Flags%ICASSCF==0.and.Flags%SaptLevel/=10) then
       ! TEST IS MISSING!
       ! GVB
       SAPT%iPINO = 3
    else
       write(LOUT,'(/,1x,a)') 'UNRECOGNIZED PINO VARIANT!'
    endif
 endif

! calculate exchange K[PB] matrix
 allocate(SAPT%monB%Kmat(NBasis,NBasis))
 allocate(work(NBasis,NBasis))
 ! work = PB
 work = 0
 call get_den(NBasis,SAPT%monB%CMO,SAPT%monB%Occ,1d0,work)
 if(Flags%ICholesky==0) then
    call make_K(NBasis,work,SAPT%monB%Kmat)
 elseif(Flags%ICholesky==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_K_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                      NCholesky,NBasis,work,SAPT%monB%Kmat)
 endif
 deallocate(work)

! calculate electrostatic potential
 call calc_elpot(SAPT%monA,SAPT%monB,CholeskyVecs,Flags%ICholesky,NBasis)

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ca,Cb)
 if(SAPT%InterfaceType==2) then
    deallocate(OneRdmB,OneRdmA,AuxB,AuxA)
 endif

end subroutine sapt_interface

subroutine onel_molpro(mon,NBasis,NSq,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData)     :: SAPT
 type(SystemBlock)  :: MonBlock
 integer,intent(in) :: mon,NBasis,NSq,NInte1

 integer                       :: ione,ios,NSym,NBas(8),ncen
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
 double precision, allocatable :: Kmat(:)
 double precision, allocatable :: work1(:),work2(:)
 character(8)                  :: label
 character(:),allocatable      :: infile,outfile

 if(mon==1) then
   infile  = 'AOONEINT_A'
   outfile = 'ONEEL_A'
 elseif(mon==2) then
   infile  = 'AOONEINT_B'
   outfile = 'ONEEL_B'
 endif

 allocate(work1(NInte1),work2(NSq))
 allocate(Hmat(NSq),Vmat(NSq),Smat(NSq))
 ! test HLONDON
 allocate(Kmat(NSq))
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

 !! test for Heitler-London
 !call readoneint_molpro(work1,infile,'KINETINT',.false.,NInte1)
 !call square_oneint(work1,Kmat,NBasis,NSym,NBas)

 MonBlock%NSym = NSym
 MonBlock%NSymBas(1:NSym) = NBas(1:NSym)

 !square form
 call writeoneint(outfile,NSq,Smat,Vmat,Hmat)

 deallocate(work2,work1)
 deallocate(Smat,Vmat,Hmat)
 deallocate(Kmat)

end subroutine onel_molpro

subroutine onel_dalton(mon,NBasis,NSq,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData) :: SAPT
 type(SystemBlock) :: MonBlock

 integer,intent(in) :: mon,NBasis,NSq,NInte1
 integer :: ione,NSym,NBas(8),ncen
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
 call readoneint(ione,work1)
 call square_oneint(work1,Hmat,NBasis,NSym,NBas)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call square_oneint(work1,work2,NBasis,NSym,NBas)
 Vmat(:) = Hmat - work2

 call readlabel(ione,'OVERLAP ')
 call readoneint(ione,work1)
 call square_oneint(work1,Smat,NBasis,NSym,NBas)

 call readlabel(ione,'ISORDK  ')
 read(ione)
 read(ione) MonBlock%charg,ncen,MonBlock%xyz

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
    call read_syminf(SAPT%monA,SAPT%monB,NBasis)

    call arrange_oneint(Smat,NBasis,SAPT)
    call arrange_oneint(Vmat,NBasis,SAPT)
    call arrange_oneint(Hmat,NBasis,SAPT)

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
character(:),allocatable :: occfile,sirifile,siriusfile,coefile


 if(Mon%Monomer==1) then
   coefile='coeff_A.dat'
   occfile='occupations_A.dat'
   sirifile='SIRIFC_A'
   siriusfile='SIRIUS_A.RST'
 elseif(Mon%Monomer==2) then
   coefile='coeff_B.dat'
   occfile='occupations_B.dat'
   sirifile='SIRIFC_B'
   siriusfile='SIRIUS_B.RST'
 endif

 inquire(file=sirifile,EXIST=exsiri)
 if(exsiri) then
    open(newunit=isiri,file=sirifile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(isiri,'TRCCINT ')
    read(isiri) NSym,NOrbt,NBasist,NCMOt,NOcc(1:NSym),NOrbs(1:NSym)

    Mon%NOrb = NOrbt
    Mon%NSymOrb(1:NSym) = NOrbs(1:NSym)

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
 else
    NBasist = NBasis
 endif

 !if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE/=1.and.(.not.Mon%ISHF)) then
 !  ! CASSCF
 !   call readmulti(NBasis,Mon,.false.,exsiri,isiri,occfile,siriusfile)
 !elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE==1.and.Flags%ISERPA==0.and.(.not.Mon%ISHF)) then
 !  ! CASSCF
 !  ! for 2-el electron case: read from occupations.dat
 !  call readmulti(NBasis,Mon,.false.,.false.,isiri,occfile,siriusfile)

 !elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE==1.and.Flags%ISERPA==2) then
 !
 !   call readmulti(NBasis,Mon,.false.,exsiri,isiri,occfile,siriusfile)

 if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.(.not.Mon%ISHF)) then

    ! CASSCF

    if(exsiri) close(isiri)

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

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

integer :: NBasis
integer :: NInte1,HlpDim,NOccup,nact
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

 allocate(Mon%CICoef(NBasis),Mon%IGem(NBasis),Mon%Occ(NBasis))
 allocate(work(HlpDim),EVal(NBasis))
 OneRdm = 0
 ! HERE! FIRST STATE FOR NOW
 call read_1rdm_molpro(OneRdm,Mon%InSt(1,1),Mon%InSt(2,1),&
                       rdmfile,Mon%IWarn,NBasis)

 call triang_to_sq2(OneRdm,OrbAux,NBasis)
 call Diag8(OrbAux,NBasis,NBasis,Eval,work)
! KP : it may happen that an active orbital has a negative tiny occupation. set it to a positive
 do i=1,Nbasis
 Eval(i)=Abs(Eval(i))
 enddo
! call dsyev('V','U',NBasis,OrbAux,NBasis,EVal,work,3*NBasis,info)
 call SortOcc(EVal,OrbAux,NBasis)

! read NAct from 1RDM
 if(Mon%NActFromRDM) Mon%NAct = 0
 Tmp = 0
 do i=1,NBasis
    Tmp = Tmp + EVal(i)
    !if(Mon%NActFromRDM.and.EVal(i)>0.d0) Mon%NAct = Mon%NAct + 1
    if(Mon%NActFromRDM.and.EVal(i)>1d-10) Mon%NAct = Mon%NAct + 1
 enddo

! test NAct from 1RDM
 call read_nact_molpro(nact,rdmfile)
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

subroutine read_mo_dalton(cmo,nbasis,nsym,nbas,norb,nsiri,nmopun)
!
! reads MOs either from SIRIUS.RST or DALTON.MOPUN
! unpacks symmetry blocks to (NBasis,NOrb) form
! in SAPT orbitals kept in AOMO order!
!
implicit none

integer,intent(in) :: nbasis,nsym,nbas(8),norb(8)
character(*)       :: nsiri,nmopun
double precision   :: cmo(nbasis,nbasis)

integer          :: i,j,idx
integer          :: iunit,irep
integer          :: ncmot
integer          :: off_i,off_j
logical          :: isiri
character(60)    :: line
double precision :: natocc(10)
double precision :: tmp(nbasis**2)

tmp =    0
ncmot    = sum(nbas(1:nsym)*norb(1:nsym))

inquire(file=nsiri,EXIST=isiri)
if(isiri) then

   open(newunit=iunit,file=nsiri,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

   call readlabel(iunit,'NEWORB  ')
   read(iunit) tmp(1:ncmot)

!   cmo = 0
!   off_i = 0
!   off_j = 0
!   idx = 0
!
!   do irep=1,nsym
!      do j=off_j+1,off_j+norb(irep)
!
!         do i=off_i+1,off_i+nbas(irep)
!            idx = idx + 1
!            cmo(i,j) = tmp(idx)
!         enddo
!
!      enddo
!      off_i = off_i + nbas(irep)
!      off_j = off_j + norb(irep)
!   enddo
   write(LOUT,'(1x,a)') 'Orbitals read from '//nsiri

else

   open(newunit=iunit,file=nmopun,form='FORMATTED',status='OLD')
   read(iunit,'(a60)') line
   off_i = 0
   idx   = 0
   do irep=1,nsym
      do j=1,norb(irep)
         read(iunit,'(4f18.14)') (tmp(off_i+i),i=1,nbas(irep))
         off_i = off_i + nbas(irep)
      enddo
   enddo

   write(LOUT,'(1x,a)') 'Orbitals read from '//nmopun
endif

close(iunit)

! unpack orbitals
cmo   = 0
off_i = 0
off_j = 0
idx   = 0

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

! call print_sqmat(cmo,nbasis)

end subroutine read_mo_dalton

subroutine readocc_cas_siri(mon,nbas,noSiri)
!
! From SIRIFC
! a) read number of active inactive orbs for SAPT-DALTON
!    total: NAct and INAct
!    in a given symmetry: INActS(1:NSym), NActS(1:NSym)
! From SIRIUST.RST
! b) read occupation numbers
!
implicit none

type(SystemBlock)   :: mon
integer,intent(in)  :: nbas
logical,intent(out) :: noSiri

logical           :: ioccsir,exsiri
integer           :: i,iunit,ios
integer           :: isym,off_i,off_a,off_x
integer           :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
                     NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,       &
                     NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBASM(8)

double precision             :: sum1,sum2
double precision,allocatable :: OccX(:)
character(:),allocatable     :: sirfile,sirifile

 ! set filnames
 if(Mon%Monomer==1) then
    sirfile  = 'SIRIUS_A.RST'
    sirifile = 'SIRIFC_A'
 elseif(Mon%Monomer==2) then
    sirfile  = 'SIRIUS_B.RST'
    sirifile = 'SIRIFC_B'
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
                 NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
                 NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBASM

    close(iunit)

    mon%INAct = nisht
    mon%NAct  = nasht

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

 ! CASCF
 inquire(file=sirfile,EXIST=ioccsir)
 if(ioccsir) then

    noSiri=.false.
    allocate(mon%Occ(nbas))
    allocate(OccX(1:norbt))

    open(newunit=iunit,file=sirfile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

    call readlabel(iunit,'NATOCC  ')
    read(iunit) OccX(1:NORBT)

    close(iunit)

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

    write(LOUT,'(1x,a,i2,a)') 'Occupancies for monomer',mon%Monomer,' read from '// sirfile

 else

    noSiri=.true.
    write(lout,'(1x,a)') 'SIRIUS.RST not available!'
    return

 endif

 ! Hartree-Fock case
 !mon%Occ = 0d0
 !mon%Occ(1:mon%NAct+mon%INAct) = 2d0

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

!subroutine readmulti(nbas,mon,ihf,exsiri,isiri,occfile,occsir)
!!
!! a) read number of active  inactive orbs for SAPT-DALTON
!! total: NAct and INAct
!! in a given symmetry: INActS(1:NSym), NActS(1:NSym)
!! info is either from SIRIFC or from occupations.dat
!!
!! b) read occupation nubmers either from SIRIUS.RST or from occupations.dat
!!
!! c) set IGem
!!
!implicit none
!
!type(SystemBlock) :: mon
!integer           :: isiri, nbas
!character(*)      :: occfile, occsir
!logical           :: ihf
!logical           :: exsiri, ioccsir
!
!logical           :: iocc
!integer           :: i,iunit,ios
!integer           :: NAct,INAct,NActS(8),INActS(8)
!integer           :: irep,TotEl,offset
!integer           :: istate,ispin,nactel,lsym
!integer           :: isym,off_i,off_a,off_x
!double precision  :: Occ(nbas),sum1,sum2
!double precision  :: potnuc,emy,eactiv,emcscf
!double precision,allocatable :: OccX(:)
!integer           :: NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
!                     NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
!                     NSYM,MULD2H(8,8),NRHF(8),NFRO(8),NISH(8),NASH(8),NORB(8),NBASM(8)
!
! allocate(mon%CICoef(nbas),mon%IGem(nbas),mon%Occ(nbas))
! !ioccsir=.false.
! if(exsiri) then
!
!    rewind(isiri)
!    read (isiri)
!    read (isiri) potnuc,emy,eactiv,emcscf, &
!                 istate,ispin,nactel,lsym
!!    read (isiri) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
!    read (isiri) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,&
!              NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,&
!              NSYM,MULD2H,NRHF,NFRO,NISH,NASH,NORB,NBASM
!
!!    print*,    potnuc,emy,eactiv,emcscf, &
!!               istate,ispin,nactel,lsym
!!   print*, 'READM-TEST'
!!   print*, 'NASHT',nasht
!!   write(lout,*) 'nish',nish(1:nsym)
!!   write(lout,*) 'nash',nash(1:nsym)
!!   stop
!!   write (*,*) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
!
!    mon%INAct = nisht
!    mon%NAct  = nasht
!
!    if(NSym/=mon%NSym) stop "NSym from SIRIUS and AOONEINT do not match!"
!
!    mon%INActS(1:mon%NSym) = NISH(1:NSym)
!    mon%NActS(1:mon%NSym)  = NASH(1:NSym)
!
!    if(nbast.ne.nbas) then
!      write(LOUT,'(1x,a)') 'WARNING! NBasis FROM SIRIFC DOES NOT MATCH!'
!      write(LOUT,'(1x,a,i5,1x,a,i5)') 'NBasis: ',nbas, 'SIRIFC: ', nbast
!      write(LOUT,'()')
!      mon%IWarn = mon%IWarn + 1
!    endif
!
!   if(.not.ihf) then
!
!      ! CASCF case
!      allocate(OccX(1:norbt))
!
!      inquire(file=occsir,EXIST=ioccsir)
!      if(ioccsir) then
!         open(newunit=iunit,file=occsir,status='OLD', &
!              access='SEQUENTIAL',form='UNFORMATTED')
!
!         call readlabel(iunit,'NATOCC  ')
!         read(iunit) OccX(1:NORBT)
!
!         close(iunit)
!
!         ! save occupations in mon%Occ
!         ! order from sym ordering to inact-act (ISW/ISX in Dalton)
!         mon%Occ = 0d0
!         off_i = 0
!         off_a = NISHT
!         off_x = 0
!         do isym=1,NSym
!            mon%Occ(off_i+1:off_i+NISH(isym)) = OccX(off_x+1:off_x+NISH(isym))
!            mon%Occ(off_a+1:off_a+NASH(isym)) = OccX(off_x+NISH(isym)+1:off_x+NISH(isym)+NASH(isym))
!            off_i = off_i + NISH(isym)
!            off_a = off_a + NASH(isym)
!            off_x = off_x + NORB(isym)
!         enddo
!
!      endif
!      deallocate(OccX)
!
!   elseif(ihf) then
!      ! Hartree-Fock case
!      mon%Occ = 0d0
!      mon%Occ(1:mon%NAct+mon%INAct) = 2d0
!
!   endif
!
!   sum1 = 0d0
!   do i=1,mon%INAct+mon%NAct
!       mon%Occ(i) = mon%Occ(i)/2d0
!       sum1 = sum1 + mon%Occ(i)
!   enddo
!   mon%SumOcc = sum1
!
!   write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occsir
!   write(LOUT,'()')
!
! endif
!
!! occupations.dat
! inquire(file=occfile,EXIST=iocc)
! if(iocc) then
!
!    Occ    = 0d0
!    INActS = 0
!    NActS  = 0
!    open(newunit=iunit,file=occfile,form='FORMATTED',status='OLD')
!
!    read(iunit,*) INAct, NAct
!    INAct = INAct/2
!    read(iunit,*) (Occ(i),i=1,INAct+NAct)
!    sum2 = 0d0
!    do i=1,INAct+NAct
!       Occ(i) = Occ(i)/2d0
!       sum2 = sum2 + Occ(i)
!    enddo
!
!    ! active and inactive orbs in each symmetry
!    read(iunit,*,iostat=ios) (INActS(i),i=1,mon%NSym)
!    if(ios==0) then
!       read(iunit,*) (NActS(i),i=1,mon%NSym)
!       mon%NActS(1:mon%NSym)  = NActS(1:mon%NSym)
!       mon%INActS(1:mon%NSym) = INActS(1:mon%NSym)
!    endif
!
!    if(mon%NSym.gt.1) then
!      call sort_sym_occ(nbas,mon%NSym,INAct,NAct,Occ)
!    endif
!
!    if(.not.ioccsir) then
!    !   if(Abs(sum2-mon%XELE).gt.1.0d-8) then
!    !      write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'
!    !      write(LOUT,'(1x,a,1x,f10.6,5x,a,i3)') 'SUM(OCC): ', sum2, 'MONOMER: ', mon%Monomer
!    !      write(LOUT,'(1x,a)') 'CHECK occupations.dat!'
!    !      stop
!    !   endif
!       mon%INAct = INAct
!       mon%NAct  = NAct
!       mon%Occ   = Occ
!       mon%SumOcc = sum2
!       ! SIRIUS.RST not there
!       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occfile
!       write(LOUT,'()')
!
!    else !compare SIRIUS.RST and occupations.dat
!       if(Abs(sum1-mon%XELE).gt.1.0d-8) then
!
!          if(Abs(sum2-mon%XELE).gt.1.0d-8) then
!             write(LOUT,'(1x,a,1x,i3)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE! MONOMER: ', mon%Monomer
!             write(LOUT,'(1x,a,1x,f10.6,5x,a,4x,f10.6)') 'OCC(SIRIUS): ', sum1,&
!                          'OCC(occupations.dat)', sum2
!             stop
!          endif
!
!          mon%INAct = INAct
!          mon%NAct  = NAct
!          mon%Occ   = Occ
!          print*, 'sum1 zle'
!       endif
!
!       ! both files correct
!       if(any(abs(Occ-mon%Occ).gt.1.d-9)) then
!          write(LOUT,'(1x,a)') 'WARNING! DIFFERENT OCCUPANCIES IN SIRIUS.RST&
!                & AND occupations.dat!'
!          write(LOUT,'(1x,a)') 'OCCUPANCIES READ FROM occupations.dat!'
!          write(LOUT,'()')
!          mon%IWarn = mon%IWarn + 1
!          mon%INAct = INAct
!          mon%NAct  = NAct
!          mon%Occ   = Occ
!          mon%SumOcc = sum2
!       endif
!
!       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occsir
!       write(LOUT,'()')
!
!    endif
!
!    close(iunit)
!
! elseif(ioccsir) then
!    if(Abs(sum1-mon%XELE).gt.1.0d-8) then
!       write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'
!       write(LOUT,*) 'Occ: ', sum1
!       write(LOUT,'(1x,a)') 'CHECK DALTON CALCULATIONS!'
!       stop
!    endif
!
! elseif(.not.ioccsir) then
!     write(LOUT,'(1x,a)') 'ERROR! CANNOT READ OCCUPANCIES!'
!     stop
! endif ! occupations.dat ?
!
!
! if(mon%INAct==0) then
!    mon%NGem = 2
!
!    mon%IGem(1:mon%NAct+mon%INAct)      = 1
!    mon%IGem(mon%NAct+mon%INAct+1:nbas) = 2
! else
!    mon%NGem = 3
!    mon%IGem(1:mon%INAct) = 1
!    mon%IGem(mon%INAct+1:mon%INAct+mon%NAct) = 2
!    mon%IGem(mon%INAct+mon%NAct+1:nbas)      = 3
! endif
!
! ! construct CICoef
! do i=1,nbas
!    mon%CICoef(i)=sqrt(mon%Occ(I))
!    if(mon%Occ(i).lt.0.5d0) mon%CICoef(i)=-mon%CICoef(i)
! enddo
!
! !do i=1,NBas
! !   print*, i,mon%Occ(i),Occ(i)
! !enddo
!
!! check
!!  write(LOUT,*) mon%Occ
!
!end subroutine readmulti
!
subroutine arrange_mo(mat,nbas,SAPT)
implicit none

type(SaptData) :: SAPT
!integer :: NOrbA,NOrbB
integer :: nbas
double precision :: mat(nbas,nbas)

call gen_swap_rows(mat,nbas,SAPT%monA%NSym,&
                   SAPT%monA%NMonBas,SAPT%monB%NMonBas)

!call swap_rows(NOrbA,NOrbB,mat)

end subroutine arrange_mo

subroutine read_syminf(A,B,nbas)
! reads number of basis functions on each monomer
! from SYMINFO(B) file!
implicit none

type(SystemBlock) :: A, B
integer :: nbas
integer :: iunit,ios
integer :: ibas,icen,last_ibas,last_icen
integer :: irep,ifun,offset
logical :: ex,dump
integer :: tmp
integer :: ACenTst, ACenBeg, ACenEnd

! sanity checks : in
!print*, A%NCen, B%NCen
!print*, A%UCen, B%UCen
if(A%NSym/=B%NSym) then
  write(lout,*) 'ERROR in read_syminf: NSym different for A and B!'
endif

inquire(file='SYMINFO_B',EXIST=ex)

if(ex) then
   open(newunit=iunit,file='SYMINFO_B',status='OLD',&
        form='FORMATTED')
   read(iunit,*)
   read(iunit,*)

   ! old version: does not work with sym
   ! print*, 'old version'
   ! offset = 0
   ! irep   = 1
   ! read(iunit,'(i5,i6)',iostat=ios) last_ibas,last_icen
   ! do
   !   read(iunit,'(i5,i6)',iostat=ios) ibas,icen
   !   if(ios/=0) then
   !      A%NMonBas(irep)=last_ibas-offset
   !      exit
   !   elseif(icen/=last_icen) then
   !        if(last_icen==B%UCen) then
   !           B%NMonBas(irep) = last_ibas-offset
   !           offset = last_ibas
   !        elseif(icen==1) then
   !           A%NMonBas(irep) = last_ibas-offset
   !           offset = last_ibas
   !           irep   = irep + 1
   !        endif
   !   endif
   !   last_ibas=ibas
   !   last_icen=icen
   !enddo

   ! new version : ok with sym
   do irep=1,B%NSym
      do ifun=1,B%NSymOrb(irep)
         read(iunit,'(i5,i6)',iostat=ios) ibas,icen
         if(icen.le.B%UCen) then
            B%NMonBas(irep) = B%NMonBas(irep) + 1
         else
            A%NMonBas(irep) = A%NMonBas(irep) + 1
         endif
      enddo
   enddo

   close(iunit)
else
   write(LOUT,'(1x,a)') 'ERROR! MISSING SYMINFO_B FILE!'
   stop
endif

! sanity checks : out
do irep=1,B%NSym
   ibas = A%NMonBas(irep)+B%NmonBas(irep)
   if(ibas/=A%NSymOrb(irep)) then
      write(lout,'(1x,a)') 'ERROR in read_syminf!'
      write(lout,'(1x,a,i3,a)') 'For irep =',irep, ':'
      write(lout,*) 'A-NMonBas',A%NMonBas(1:A%NSym)
      write(lout,*) 'B-NMonBas',B%NMonBas(1:B%NSym)
      write(lout,*) 'Sum:     ',A%NMonBas(1:A%NSym)+B%NMonBas(1:B%NSym)
      write(lout,*) 'Should be',A%NSymOrb(1:A%NSym)
      stop
   endif
enddo

end subroutine read_syminf

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

subroutine sort_sym_mo(CMO,nbas,mon)
implicit none

type(SystemBlock)              :: mon
integer,intent(in)             :: nbas
double precision,intent(inout) :: CMO(nbas,nbas)

integer                      :: i,j,ii
integer                      :: TotEl,TotElIrep,irep,idx
integer,allocatable          :: ICpy1(:),ICpy2(:)
integer,allocatable          :: LabelAct(:),LabelIAct(:)
double precision,allocatable :: COrd(:,:)

 TotEl = mon%INAct + mon%NAct

! print*, 'NSym   ',mon%NSym
! print*, 'NSymOrb',mon%NSymOrb
! print*, 'INActS',mon%InActS(1:mon%NSym)
! print*, 'NActS ',mon%NActS(1:mon%NSym)

 allocate(ICpy1(nbas),ICpy2(nbas))
 allocate(LabelAct(nbas),LabelIAct(nbas),COrd(nbas,nbas))

 ICpy1 = 0
 ICpy2 = 0

 ! make labels
 idx = 0
 do irep=1,mon%NSym
    TotElIrep = mon%INActS(irep)+mon%NActS(irep)
    do j=1,mon%NSymOrb(irep)

       idx = idx + 1
       LabelAct(idx) = 0

       if(j.gt.mon%INActS(irep).and.j.le.TotElIrep) then
          LabelAct(idx) = 1
       endif
       LabelIAct(idx) = 0

       if(j.le.mon%INActS(irep)) LabelIAct(idx) = 1

    enddo
 enddo

 do ii=1,nbas

   ! inactive
   do i=1,nbas
      if(LabelIAct(i).eq.1.and.ICpy2(i).eq.0.and.ICpy1(ii).eq.0) then
         ICpy2(i)  = 1
         ICpy1(ii) = 1

         do j=1,nbas
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
            enddo
         endif
      enddo
   endif

 enddo

 CMO = COrd

 deallocate(COrd,LabelIAct,LabelAct)
 deallocate(ICpy2,ICpy1)

end subroutine sort_sym_mo

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

subroutine prepare_no(OneRdm,OrbAux,OrbCAS,Mon,CholeskyVecs,IFunSR,ICholesky,NBasis)
implicit none
!
! OrbCAS[inout] :: on input AOtoCAS
!                  on output AOtoNO
! OrbAux        :: on input CAStoNO
!
type(SystemBlock)   :: Mon
type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: IFunSR,ICholesky,NBasis
double precision   :: OneRdm(NBasis*(NBasis+1)/2)
double precision   :: OrbAux(NBasis,NBasis),OrbCAS(NBasis,NBasis)

integer :: NOccup,NVirt,NSym
integer :: NCholesky
integer :: i,j,ia,ib,iab,ioff,idx,NInte1
double precision,allocatable :: URe(:,:),OrbSym(:,:),Fock(:)
double precision,allocatable :: work1(:),work2(:),work3(:)
character(:),allocatable :: onefile,rdmfile,aoerfile
! testy
integer :: info

 NInte1 = NBasis*(NBasis+1)/2
 NOccup = Mon%INAct + Mon%NAct
 NVirt = NBasis - Mon%INAct - Mon%NAct

 if(Mon%Monomer==1) then
   onefile = 'AOONEINT_A'
   rdmfile = '2RDMA'
   aoerfile = 'AOERFSORT'
 elseif(Mon%Monomer==2) then
   onefile = 'AOONEINT_B'
   rdmfile = '2RDMB'
   if(Mon%SameOm) then
      aoerfile = 'AOERFSORT'
   else
      aoerfile = 'AOERFSORTB'
   endif
 endif

 allocate(Mon%NumOSym(15),Mon%IndInt(NBasis))
 allocate(work1(NInte1),work2(NInte1),work3(NBasis),&
          Fock(NBasis**2),OrbSym(NBasis,NBasis),URe(NBasis,NBasis))

 call create_ind(rdmfile,Mon%NumOSym,Mon%IndInt,NSym,NBasis)

! COPY AUXM TO URe AND OFF SET BY NInAc
 URe = 0
 forall(i=1:NBasis) URe(i,i)=1d0
 ! with Diag8:
 do i=1,Mon%NAct
    do j=1,Mon%NAct
       URe(Mon%INAct+i,Mon%INAct+j) = OrbAux(i,j)
    enddo
 enddo
 ! with dsyev
 !do i=1,Mon%NAct
 !   do j=1,Mon%NAct
 !      URe(Mon%INAct+i,Mon%INAct+j) = OrbAux(NBasis+1-i,NBasis+1-j)
 !   enddo
 !enddo
!print*, norm2(URe)
! call print_sqmat(URe,NBasis)

! FIND CANONICAL INACTIVE AND VIRTUAL ORBITALS
 work1 = 0
 idx = 0
 do j=1,Mon%INAct
    do i=1,j
       idx = idx + 1
       if(i==j) work1(idx) = 1.0d0
    enddo
 enddo
 idx = 0
 do j=1,Mon%NAct
    do i=1,j
       idx = idx + 1
       ioff = (Mon%INAct+j)*(Mon%INAct+j-1)/2 + Mon%INAct
       work1(ioff+i) = OneRdm(idx)
    enddo
 enddo
! do i=1,NInte1
!    print*, i,work1(i)
! enddo

 do i=1,NBasis
    do j=1,NBasis
       OrbSym(Mon%IndInt(i),j) = OrbCAS(j,i)
    enddo
 enddo

 iab = 0
 do ia=1,NBasis
    do ib=1,ia
       iab = iab + 1
       OneRdm(iab) = 0.d0
       do i=1,NBasis
          do j=1,NBasis
             idx = max(i,j)*(max(i,j)-1)/2+min(i,j)
             OneRdm(iab) = OneRdm(iab) &
           + OrbSym(i,ia)*OrbSym(j,ib)*work1(idx)
          enddo
       enddo
    enddo
 enddo

 ! create Fock matrix
 ! work1 = XOne
 call readoneint_molpro(work1,onefile,'ONEHAMIL',.true.,NInte1)
 ! work2 = Fock
 if(IFunSR==0) then
 ! CASSCF,Hartree-Fock

   if(ICholesky==0) then
     call FockGen_mithap(work2,OneRdm,work1,NInte1,NBasis,'AOTWOSORT')
   elseif(ICholesky==1) then
     NCholesky = CholeskyVecs%NCholesky
     call FockGen_CholR(work2,CholeskyVecs%R(1:NCholesky,1:NInte1),OneRdm,work1, &
                        NInte1,NCholesky,NBasis)
   endif

 elseif(IFunSR>0) then
 ! Kohn-Sham

   ! add and store Coulomb
   ! for RSH short-range Coulomb is stored
   allocate(Mon%VCoul(NInte1))
   call PotCoul_mithap(Mon%VCoul,OneRdm,Mon%doRSH,aoerfile,NBasis)
   ! RSH
   if(Mon%doRSH) then
     ! generate long-range Fock
     print*, 'Check: prepare_no:',aoerfile
     call FockGen_mithap(work2,OneRdm,work1,NInte1,NBasis,aoerfile)
     work2 = work2 + Mon%VCoul
   else
   ! non-hybrid DFAs
   !  work2 = work1
     work2 = work1 + Mon%VCoul
   endif

 endif
 call tran_matTr(work2,OrbSym,OrbSym,NBasis,.false.)

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
         URe(i,j) = Fock((j-1)*Mon%INAct+i)
      enddo
    enddo
 endif

! VIRTUAL
 if(NVirt/=0) then
    do i=1,NVirt
       do j=1,NVirt
          idx = max(i+NOccup,j+NOccup)*(max(i+NOccup,j+NOccup)-1)/2&
              + min(i+NOccup,j+NOccup)
          Fock((j-1)*NVirt+i) = work2(idx)
       enddo
    enddo
    call Diag8(Fock,NVirt,NVirt,work3,work1)
    !call dsyev('V','U',NVirt,Fock,NVirt,work3,work1,3*NVirt,info)
    do i=1,NVirt
       do j=1,NVirt
          URe(i+NOccup,j+NOccup) = Fock((j-1)*NVirt+i)
       enddo
    enddo
 endif
 ! test for ICPHF
 !print*, 'work3',work3
 Mon%OrbE(Mon%INAct+1:NBasis)=work3(1:NBasis-Mon%INAct)

! END OF CANONICALIZING

 call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,OrbSym,NBasis,0d0,OrbCAS,NBasis)
 OrbCAS = transpose(OrbCAS)

 deallocate(work3,work2,work1,Fock,OrbSym,URe)
 deallocate(Mon%IndInt)

end subroutine prepare_no

subroutine prepare_rdm2_file(Mon,OrbAux,NBasis)
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
                       rdmfile,Mon%IWarn,Mon%NAct)

 do i=1,Mon%NAct
    do j=1,Mon%NAct
       work1((j-1)*Mon%NAct+i) = OrbAux(i,j)
    enddo
 enddo
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

end subroutine prepare_rdm2_file

subroutine select_active(mon,nbas,Flags)
! set dimensions: NDimX,num0,num1,num2
! set matrices  : IndN,IndX,IPair,IndAux
!
implicit none

type(SystemBlock) :: mon
type(FlagsData) :: Flags
integer :: nbas!, ICASSCF, ISHF, IFlCore
integer :: i, j, ij, icnt
integer :: ind, ind_ij
integer :: IAuxGem(nbas)
integer :: IndHlp(nbas)
integer :: test
character(1) :: mname

 if(mon%Monomer==1) mname='A'
 if(mon%Monomer==2) mname='B'

 IAuxGem = mon%IGem
 allocate(mon%IndAux(nbas))

 do i=1,mon%NELE
    mon%IndAux(i)=0
 enddo
 do i=1+mon%NELE,nbas
    mon%IndAux(i)=2
 enddo

 if(mon%NActOrb/=0) then

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
             ! here!!!
             !if(mon%Occ(i).lt.1d0.and.mon%Occ(i).gt.1d-6) then
             ! HERE!!! ACTIVE!!!!
             mon%IndAux(i) = 1
             write(6,'(X," Active Orbital: ",I4,E14.4)') i, mon%Occ(i)
             mon%icnt = mon%icnt + 1
          endif
       enddo
    endif

 endif

! set generalized "occupied" = num0 + num1
! and "virtual" = num1 + num2 indices
if(Flags%ICASSCF==0) then

   do i=1,mon%NELE
      IndHlp(i)=0
   enddo
   do i=1+mon%NELE,nbas
      IndHlp(i)=2
   enddo

   do i=1,nbas
      if(mon%Occ(i).lt.1d0.and.mon%Occ(i).ne.0d0) then
      IndHlp(i)=1
      EndIf
   enddo

   ! MH: in GVB the num0-num1 choice is not clear to me!
   ! MH: I can't remember why we changed it from IndAux?
   mon%num0 = 0
   do i=1,nbas
      if(IndHlp(i)/=0) exit
      !if(mon%IndAux(i)/=0) exit
      mon%num0 = mon%num0 + 1
   enddo
   mon%num2 = 0
   do i=nbas,1,-1
      if(IndHlp(i)/=2) exit
      !if(mon%IndAux(i)/=2) exit
      mon%num2 = mon%num2 + 1
   enddo
   mon%num1 = nbas - mon%num0 - mon%num2

elseif(Flags%ICASSCF==1) then
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
endif

 ! num0-3 are set based on occupancies
 ! sometimes active orbitals in Molpro have 0.0 occupancy
 ! in which case we have to keep the NAct from Molpro
 if(mon%ISwitchAct==1.and.mon%IPrint>=3) then
    write(lout,'(/1x,3a,i4,a,i4)') 'In monomer ', mname, &
               ' changing num0 from ', mon%num0, ' to', mon%INAct
    write(lout,'(1x,3a,i4,a,i4/)') 'In monomer ', mname, &
               ' changing num1 from ', mon%num1, ' to', mon%NAct
    mon%num0 = mon%InAct
    mon%num1 = mon%NAct
 endif

! active pairs
 allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))

 mon%IPair(1:nbas,1:nbas) = 0

 if(Flags%ICASSCF==0) then
    write(LOUT,'(1x,a,e15.5)') 'Threshold for active orbitals:       ', mon%ThrSelAct
    write(LOUT,'(1x,a,e15.5)') 'Threshold for quasi-virtual orbitals:', mon%ThrQVirt
 ! allocate(mon%IndXh(mon%NDim))

    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1

          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
             ! do not correlate active degenerate orbitals from different geminals
             if((mon%IGem(i).ne.mon%IGem(j)).and.&
                  (mon%IndAux(i)==1).and.(mon%IndAux(j)==1).and.&
                  !(Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-2)) then
                  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.mon%ThrSelAct)) then

                write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j
             else
                ! if IFlCore=0 exclude core (inactive) orbitals
                if(Flags%IFlCore==1.or.&
                     (Flags%IFlCore==0.and.&
                     mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then

                     if(abs(mon%Occ(i)+mon%Occ(j)).lt.mon%ThrQVirt) then
                        write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly virtual-orbitals pair',i,j
                     elseif(abs(mon%Occ(i)+mon%Occ(j)-2d0).gt.1.D-10) then

                       ind = ind + 1
                       mon%IndX(ind) = ij
                       ! mon%IndXh(ind) = ind
                       mon%IndN(1,ind) = i
                       mon%IndN(2,ind) = j
                       mon%IPair(i,j) = 1
                       mon%IPair(j,i) = 1
                     endif

                endif
             endif

          endif

       enddo
    enddo

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then
    write(LOUT,'(1x,a,e15.5)')  'Threshold for active orbitals:       ', mon%ThrSelAct
    write(LOUT,'(1x,a,2e15.5)') 'Threshold for quasi-virtual orbitals:', mon%ThrQVirt

    if(mon%NCen==1.and.mon%ThrSelAct<1.d-3.and.mon%NAct>1) then
       write(LOUT,'(1x,a)') 'Warning! For single atom ThrSelAct should probably have larger value!'
       mon%IWarn = mon%IWarn + 1
    endif

    ij  = 0
    ind = 0
    do i=1,nbas
       do j=1,i-1

          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
             ! do not correlate active degenerate orbitals from different geminals
             if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1)  &
                  .and.&
                  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.mon%ThrSelAct) ) then
                ! here!!!
                !  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1d-3) ) then

                write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j
             else
                ! if IFlCore=0 exclude core (inactive) orbitals
                if(Flags%IFlCore==1.or.&
                     (Flags%IFlCore==0.and.&
                     mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
                     ! exclude pairs of nearly/virtual orbitals
                     if(abs(mon%Occ(i)+mon%Occ(j)).lt.mon%ThrQVirt) then
                     !if(abs(mon%Occ(i)+mon%Occ(j)).lt.1.D-7) then
                        write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly virtual-orbitals pair',i,j
                     elseif(abs(mon%Occ(i)+mon%Occ(j)-2d0).gt.1.D-10) then

                        ind = ind + 1
                        mon%IndX(ind) = ind
                        mon%IndN(1,ind) = i
                        mon%IndN(2,ind) = j
                        mon%IPair(i,j) = 1
                        mon%IPair(j,i) = 1
                     endif

                endif
             endif

          endif

       enddo
    enddo

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==2.and.mon%NELE==1) then

    allocate(mon%IndXh(mon%NDim))
    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1

          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
             !!! do not correlate active degenerate orbitals from different geminals
             !if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1)  &
             ! .and.&
             ! (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
             ! write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j
             !else
             !!! if IFlCore=0 exclude core (inactive) orbitals
             if(Flags%IFlCore==1.or.&
                  (Flags%IFlCore==0.and.&
                  mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then

                ind = ind + 1
                mon%IndX(ind) =  ij !ind
                mon%IndXh(ind) = ij
                mon%IndN(1,ind) = i
                mon%IndN(2,ind) = j
                mon%IPair(i,j) = 1
                mon%IPair(j,i) = 1

             endif
             !endif
          endif
       enddo
    enddo

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==2.and.mon%NELE/=1) then
    write(LOUT,'(1x,a)') 'WARNING!!! Be!'
    write(LOUT,'(1x,a,e15.5)') 'Threshold for active orbitals: ', mon%ThrSelAct
    ! write(*,*) 'Be?, ONLY ACTIVE PAIRS!'
    allocate(mon%IndXh(mon%NDim))

    ij=0
    ind = 0
    do i=1,nbas
       do j=1,i-1

          ij = ij + 1
          ind_ij = mon%IndAux(i)+mon%IndAux(j)
          if((ind_ij/=0).and.(ind_ij/=4)) then
          ! special test
          !if(mon%IndAux(i)==1.and.mon%IndAux(j)==1) then
             ! do not correlate active degenerate orbitals from different geminals
             if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1)  &
              .and.&
              !(Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-10) ) then
              (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.Mon%ThrSelAct) ) then
              write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j
             else
             ! if IFlCore=0 exclude core (inactive) orbitals
             if(Flags%IFlCore==1.or.&
                  (Flags%IFlCore==0.and.&
                  mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then

                ind = ind + 1
                ! active Be
                mon%IndX(ind) = ij !ind
                mon%IndXh(ind) = ij
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

subroutine calc_elpot(A,B,CholeskyVecs,ICholesky,NBas)
implicit none

type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: ICholesky,NBas

integer :: NInte1,NCholesky
double precision,allocatable :: Pa(:,:),Pb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Ja(:,:),Jb(:,:)

 NInte1 = NBas*(NBas+1)/2

 allocate(Pa(NBas,NBas),Pb(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Ja(NBas,NBas),Jb(NBas,NBas))

 call get_den(NBas,A%CMO,A%Occ,2d0,Pa)
 call get_den(NBas,B%CMO,B%Occ,2d0,Pb)

 call get_one_mat('V',Va,A%Monomer,NBas)
 call get_one_mat('V',Vb,B%Monomer,NBas)

 if(ICholesky==0) then
    call make_J2(NBas,Pa,Pb,Ja,Jb)
 elseif(ICholesky==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_J2_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                       Pa,Pb,Ja,Jb,NCholesky,NBas)
 endif

 allocate(A%WPot(NBas,NBas),B%WPot(NBas,NBas))

 A%WPot = Va + Ja
 B%WPot = Vb + Jb

 deallocate(Jb,Ja,Vb,Va,Pb,Pa)

end subroutine calc_elpot

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

subroutine chol_sapt_NOTransf(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
! transform Choleksy Vecs from AO to NO
! for all 2-index vecs needed in SAPT:
!   FFXX(NCholesky,NBas**2) -- for polarization
!   FFXY(NCholesky,NBas**2) -- for exchange
!   FFYX(NCholesky,NBas**2) -- for exchange
!   OOXX(NCholesky,dimO**2) -- for polarization
!   DCholX(NCholesky,NDimX) -- for polarization
!
implicit none

type(SaptData)   :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: dimOA,dimOB,dimVA,dimVB, &
           nOVA,nOVB
integer          :: i,j,ip,iq,ipq
double precision :: Cpq
double precision,allocatable :: tmp(:,:)
! test
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

 ! set buffer size
 if(MemType == 2) then       !MB
    MaxBufferDimMB = MemVal
 elseif(MemType == 3) then   !GB
    MaxBufferDimMB = MemVal * 1024_8
 endif
 write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

 NCholesky = CholeskyVecs%NCholesky
 dimOA = A%num0+A%num1
 dimOB = B%num0+B%num1
 dimVA = A%num1+A%num2
 dimVB = B%num1+B%num2
 nOVA  = dimOA*dimVA
 nOVB  = dimOB*dimVB

 print*, 'dimOA',dimOA
 print*, 'dimOB',dimOB

 !allocate(A%OV(NCholesky,A%NDimX),B%OV(NCholesky,B%NDimX))

 !allocate(tmp(NCholesky,nOVA))
 !! (OV|AA)
 !call chol_MOTransf(tmp,CholeskyVecs,&
 !                   A%CMO,1,dimOA,&
 !                   A%CMO,A%num0+1,NBasis)
 !A%OV = 0
 !do i=1,A%NDimX
 !   ip  = A%IndN(1,i)
 !   iq  = A%IndN(2,i)
 !   ipq = iq+(ip-A%num0-1)*dimOA
 !   A%OV(:,i)= tmp(:,ipq)
 !enddo
 !print*, 'A-OV',norm2(A%OV)
 !deallocate(tmp)

 !allocate(tmp(NCholesky,nOVB))
 !! (OV|BB)
 !call chol_MOTransf(tmp,CholeskyVecs,&
 !                   B%CMO,1,dimOB,&
 !                   B%CMO,B%num0+1,NBasis)

 !B%OV = 0
 !do i=1,B%NDimX
 !   ip  = B%IndN(1,i)
 !   iq  = B%IndN(2,i)
 !   ipq = iq+(ip-B%num0-1)*dimOB
 !   B%OV(:,i)= tmp(:,ipq)
 !enddo

 !deallocate(tmp)

 allocate(A%OO(NCholesky,dimOA**2),&
          B%OO(NCholesky,dimOB**2) )
 ! (OO|AA)
 !call chol_MOTransf(A%OO,CholeskyVecs,&
 !                  A%CMO,1,dimOA,&
 !                  A%CMO,1,dimOA)
 !                  B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%OO,CholeskyVecs,&
                    A%CMO,1,dimOA,&
                    A%CMO,1,dimOA,&
                    MaxBufferDimMB)
call clock('AOO',Tcpu,Twall)
 ! (OO|BB)
 !call chol_MOTransf(B%OO,CholeskyVecs,&
 !                   B%CMO,1,dimOB,&
 !                   B%CMO,1,dimOB)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%OO,CholeskyVecs,&
                    B%CMO,1,dimOB,&
                    B%CMO,1,dimOB,&
                    MaxBufferDimMB)
call clock('BOO',Tcpu,Twall)

 allocate(A%FF(NCholesky,NBasis**2),&
          B%FF(NCholesky,NBasis**2) )
 ! (FF|AA)
 !call chol_MOTransf(A%FF,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   A%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%FF,CholeskyVecs,&
                    A%CMO,1,NBasis,&
                    A%CMO,1,NBasis,&
                    MaxBufferDimMB)
 call clock('AFF',Tcpu,Twall)
 ! (FF|BB)
 !call chol_MOTransf(B%FF,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%FF,CholeskyVecs,&
                    B%CMO,1,NBasis,&
                    B%CMO,1,NBasis,&
                    MaxBufferDimMB)
 call clock('BFF',Tcpu,Twall)

 ! DChol(NCholeksy,NDimX)
 allocate(A%DChol(NCholesky,A%NDimX), &
          B%DChol(NCholesky,B%NDimX))

 do j=1,A%NDimX
    ip = A%IndN(1,j)
    iq = A%IndN(2,j)
    ipq = iq + (ip-1)*NBasis
    Cpq = A%CICoef(ip) + A%CICoef(iq)
    A%DChol(:,j) = Cpq*A%FF(:,ipq)
 enddo

 do j=1,B%NDimX
    ip = B%IndN(1,j)
    iq = B%IndN(2,j)
    ipq = iq + (ip-1)*NBasis
    Cpq = B%CICoef(ip) + B%CICoef(iq)
    B%DChol(:,j) = Cpq*B%FF(:,ipq)
 enddo

 if(SAPT%SaptLevel==999) return

 ! test for Pmat
 !allocate(SAPT%CholVecs(NBasis**2,NCholesky))
 !SAPT%CholVecs = 0
 !do iq=1,NCholesky
 !   do ip=1,NBasis**2
 !      SAPT%CholVecs(ip,iq) = CholeskyVecs%R(iq,ip)
 !   enddo
 !enddo
 !call test_pmat(SAPT%CholVecs,A,B,NBasis,NCholesky)

 allocate(A%FFAB(NCholesky,NBasis**2),&
          B%FFBA(NCholesky,NBasis**2) )
 ! (FF|AB)
 !call chol_MOTransf(A%FFAB,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   B%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(A%FFAB,CholeskyVecs,&
                    A%CMO,1,NBasis,&
                    B%CMO,1,NBasis,&
                    MaxBufferDimMB)
 ! (FF|BA)
 !call chol_MOTransf(B%FFBA,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   A%CMO,1,NBasis)
 !                   B%CMO,1,NBasis)
 !
 call chol_MOTransf_TwoStep(B%FFBA,CholeskyVecs,&
                    B%CMO,1,NBasis,&
                    A%CMO,1,NBasis,&
                    MaxBufferDimMB)

 !allocate(A%FO(NCholesky,NBasis*dimOA),&
 !         B%FO(NCholesky,NBasis*dimOA))
 !! (FO|AA)
 !call chol_MOTransf(A%FO,CholeskyVecs,&
 !                   A%CMO,1,NBasis,&
 !                   A%CMO,1,dimOA)
 !! (FO|BB)
 !call chol_MOTransf(B%FO,CholeskyVecs,&
 !                   B%CMO,1,NBasis,&
 !                   B%CMO,1,dimOB)
 !! test exchange
 !call test_Qmat(A,B,CholeskyVecs,NBasis)

end subroutine chol_sapt_NOTransf

subroutine gen_swap_rows(mat,nbas,nsym,nA,nB)
implicit none

integer,intent(in) :: nbas,nsym,nA(8),nB(8)
double precision   :: mat(nbas,nbas)
double precision   :: work(nbas,nbas)

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
double precision   :: mat(nbas,nbas)
double precision   :: work(nbas,nbas)

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
!
!
!
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

subroutine print_occ(nbas,SAPT,ICASSCF)
implicit none
! HERE : Change to A/B monomers!
type(SaptData)     :: SAPT
integer,intent(in) :: nbas, ICASSCF

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
   write(LOUT,'(2x,a,f8.4,18x,f8.4)') 'SUM OF OCCUPANCIES: ', A%SumOcc, B%SumOcc
   write(LOUT, '()')
 endif
 end associate

end subroutine print_occ

subroutine print_active(SAPT, nbas)
implicit none

type(SaptData) :: SAPT
integer        :: nbas
integer        :: i,ip,NDimX

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
    write(LOUT,'()')

    if(SAPT%IPrint.ge.10) then
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
 endif

end subroutine print_active

end module sapt_inter
