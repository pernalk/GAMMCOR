module sapt_inter

use types
use timing
use tran
use sorter
!use Cholesky_old
use Cholesky
use abmat
use read_external
use trexio

implicit none

contains

subroutine sapt_interface(Flags,SAPT,NBasis)
!
! SAPT-DALTON requires SIRIFC and SIRIUS.RST
!                   (or SIRIFC and occupations.dat)
!                      and rdm2.dat files
! SAPT-MOLPRO requires AOONEINT, AOTWOINT.mol,
!                      DIP, 2RDM files
! SAPT-TREXIO requires .h5 files
!
implicit none

type(FlagsData)     :: Flags
type(SaptData)      :: SAPT
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis

integer    :: NAO,NCholesky
integer    :: NInte1,NInte2
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

! set dimensions
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

! read and dump 1-electron integrals in AO
 if(SAPT%InterfaceType==1) then
    call onel_dalton(SAPT%monA%Monomer,NBasis,NAO,NInte1,SAPT%monA,SAPT)
    call onel_dalton(SAPT%monB%Monomer,NBasis,NAO,NInte1,SAPT%monB,SAPT)
 elseif(SAPT%InterfaceType==2) then
    call onel_molpro(SAPT%monA%Monomer,NBasis,NAO,NInte1,SAPT%monA,SAPT)
    call onel_molpro(SAPT%monB%Monomer,NBasis,NAO,NInte1,SAPT%monB,SAPT)
 elseif(SAPT%InterfaceType==5) then
    call onel_trexio(NBasis,NAO,SAPT%monA,SAPT)
    call onel_trexio(NBasis,NAO,SAPT%monB,SAPT)
 else
    write(lout,'(1x,a)') 'Unrecognized InterfaceType for SAPT!'
    stop
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
    ! AuxA and AuxB are MO-->NO coefficients
    allocate(AuxA(NBasis,NBasis),AuxB(NBasis,NBasis),&
             OneRdmA(NInte1),OneRdmB(NInte1))
    call readocc_molpro(NBasis,SAPT%monA,AuxA,OneRdmA,Flags)
    call readocc_molpro(NBasis,SAPT%monB,AuxB,OneRdmB,Flags)

 elseif(SAPT%InterfaceType==5) then
    allocate(AuxA(NBasis,NBasis),AuxB(NBasis,NBasis),&
             OneRdmA(NInte1),OneRdmB(NInte1))
    call readocc_trexio(NBasis,SAPT%monA,AuxA,OneRdmA,Flags)
    call readocc_trexio(NBasis,SAPT%monB,AuxB,OneRdmB,Flags)

 endif
 call print_occ(NBasis,SAPT,Flags%ICASSCF)

! read orbitals
! Dalton:
! norb.leq.nbas, orbitals mays be deleted due to linear
! dependecies in large basis sets; ncmot = norb*nbas
! in general, NAO /= NBasis (e.g., in TREXIO)

 allocate(Ca(NBasis*NAO),Cb(NBasis*NAO))

 if(SAPT%InterfaceType==1) then

    call read_mo_dalton(Ca,NBasis,SAPT%monA%NSym,SAPT%monA%NSymBas,SAPT%monA%NSymOrb,&
                 'SIRIUS_A.RST','DALTON_A.MOPUN')
    call read_mo_dalton(Cb,NBasis,SAPT%monB%NSym,SAPT%monB%NSymBas,SAPT%monB%NSymOrb,&
                 'SIRIUS_B.RST','DALTON_B.MOPUN')
    call arrange_mo(Cb,NBasis,SAPT)

 elseif(SAPT%InterfaceType==2) then

    call read_mo_molpro(Ca,'MOLPRO_A.MOPUN','CASORB  ',NBasis)
    call read_mo_molpro(Cb,'MOLPRO_B.MOPUN','CASORB  ',NBasis)

 elseif(SAPT%InterfaceType==5) then

    call read_mo_trexio(Ca,SAPT%monA%TrexFile,NAO,NBasis)
    call read_mo_trexio(Cb,SAPT%monB%TrexFile,NAO,NBasis)

 endif

! symmetry sorting for Dalton
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
    elseif(SAPT%InterfaceType==5) then
       call readtwoint(NAO,5,SAPT%monA%TrexFile,'AOTWOSORT',MemSrtSize)
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
    elseif(SAPT%InterfaceType==5) then
       call chol_CoulombMatrix(CholeskyVecs,NAO,SAPT%monA%TrexFile,5,Flags%ICholeskyAccu)
    endif

    SAPT%NCholesky  = CholeskyVecs%NCholesky
    SAPT%monA%NChol = SAPT%NCholesky
    SAPT%monB%NChol = SAPT%NCholesky

 endif
 call clock('2ints',Tcpu,Twall)

 if(SAPT%InterfaceType==2) then
    call prepare_no(OneRdmA,AuxA,Ca,SAPT%monA,CholeskyVecs,Flags%IFunSR,Flags%ICholesky,NBasis)
    call prepare_no(OneRdmB,AuxB,Cb,SAPT%monB,CholeskyVecs,Flags%IFunSR,Flags%ICholesky,NBasis)
    call prepare_rdm2_molpro(SAPT%monA,AuxA,NBasis)
    call prepare_rdm2_molpro(SAPT%monB,AuxB,NBasis)
 elseif(SAPT%InterfaceType==5) then
    call prepare_no_trexio(AuxA,Ca,NAO,NBasis)
    call prepare_no_trexio(AuxB,Cb,NAO,NBasis)
 endif

 allocate(SAPT%monA%CMO(NAO,NBasis),SAPT%monB%CMO(NAO,NBasis))

 call save_CAONO(Ca,SAPT%monA%CMO,NAO,NBasis)
 call save_CAONO(Cb,SAPT%monB%CMO,NAO,NBasis)

! look-up tables
 if(SAPT%InterfaceType==5) then
    call select_active_thresh(SAPT%monA,NBasis,Flags)
    call select_active_thresh(SAPT%monB,NBasis,Flags)
 else
    call select_active(SAPT%monA,NBasis,Flags)
    call select_active(SAPT%monB,NBasis,Flags)
 endif

 ! read & transform 2-rdm / trexio
 if(SAPT%InterfaceType==5) then

   !! test CIPSI energy
   !call energy_trexio_no(SAPT%monA,OneRdmA,AuxA,Ca,NAO,NBasis)
   !call energy_trexio_no(SAPT%monB,OneRdmB,AuxB,Cb,NAO,NBasis)

   ! for testing
   !call truncate_rdm2(SAPT%monA,NAO,NBasis)
   !call truncate_rdm2(SAPT%monB,NAO,NBasis)

   call rw_trexio_rdm2(SAPT%monA,AuxA,NBasis)
   call rw_trexio_rdm2(SAPT%monB,AuxB,NBasis)
 endif

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

! calculate exchange K[PB] matrix in AO
 allocate(SAPT%monB%Kmat(NAO,NAO))
 allocate(work(NAO,NAO))
!get PB density in AO
! work = 0
! do i=1,NBasis
!    call dger(NAO,NAO,SAPT%monB%Occ(i),SAPT%monB%CMO(:,i),1,SAPT%monB%CMO(:,i),1,work,NAO)
! enddo
!
 if(Flags%ICholesky==0) then
!    call make_K(NAO,work,SAPT%monB%Kmat)
   call sapt_Kmat_AO(SAPT%monB%Kmat,SAPT%monB%CMO,SAPT%monB%Occ,NAO,NBasis)
 elseif(Flags%ICholesky==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_K_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                      NCholesky,NAO,work,SAPT%monB%Kmat)
 endif
! print*, 'from sapt_interface',norm2(SAPT%monB%Kmat)
 deallocate(work)

! calculate electrostatic potential
 call calc_elpot(SAPT%monA,SAPT%monB,CholeskyVecs,Flags%ICholesky,NAO,NBasis)

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ca,Cb)
 if(SAPT%InterfaceType==2) then
    deallocate(OneRdmB,OneRdmA,AuxB,AuxA)
 endif

end subroutine sapt_interface

subroutine onel_molpro(mon,NBasis,NAO,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData)      :: SAPT
 type(SystemBlock)   :: MonBlock
 integer,intent(in)  :: mon,NBasis,NInte1
 integer,intent(out) :: NAO

 integer :: NSq
 integer                       :: ione,ios,NSym,NBas(8),ncen
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
 double precision, allocatable :: Kmat(:)
 double precision, allocatable :: work1(:),work2(:)
 character(8)                  :: label
 character(:),allocatable      :: infile,outfile

 NSq = NBasis**2

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

 NAO = NBasis
 SAPT%NAO = NAO

 !square form
 call writeoneint(outfile,NSq,Smat,Vmat,Hmat)

 deallocate(work2,work1)
 deallocate(Smat,Vmat,Hmat)
 deallocate(Kmat)

end subroutine onel_molpro

subroutine onel_dalton(mon,NBasis,NAO,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData) :: SAPT
 type(SystemBlock) :: MonBlock

 integer,intent(in)  :: mon,NBasis,NInte1
 integer,intent(out) :: NAO

 integer :: NSq
 integer :: ione,NSym,NBas(8),ncen
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:)
 double precision, allocatable :: work1(:),work2(:)
 character(:),allocatable :: infile,outfile

 ! set dimensions
 NSq = NBasis**2

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

 MonBlock%NSym = NSym
 MonBlock%NSymOrb(1:NSym) = NBas(1:NSym)

 NAO = NBasis
 SAPT%NAO = NAO

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

! write in square form
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

if(mon%IPrint.gt.10) then
  print*, 'Total number of nuclei :',num
  print*, 'Charges on centers'
  do i=1,num
     write(lout,'(i3,f12.6)') i, Mon%charg(i)
  enddo
  print*, 'Cartesian coordinates (Angstrom)'
  do i=1,num
      write(6,'(i3,3f12.6)') i, Mon%xyz(i,1:3)
  enddo
  print*, ''
endif

rc = trexio_close(f)

deallocate(coord,charge)
deallocate(kinetic)
deallocate(Smat,Vmat,Hmat)

end subroutine onel_trexio

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
!
! output: INActive, NAct
!         NGem, IGem
!
! output: Occ, CICoef
!
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
                       Mon%ISpinMs2,rdmfile,Mon%IWarn,NBasis)

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

subroutine readocc_trexio(NBasis,Mon,OrbAux,OneRdm,Flags)
!
! Purpose:
! read and 1-RDM, get Occ and MO-->NO
! Comment:
! Locally assume all orbitals are active!
!
! [output]: Mon%Occ           : occupation numbers (0:1)
!           Mon%CICoef        : CI coefficients == sqrt(Occ)
!           Mon%IGem          : number of geminals (1=inactive,2=active,3=virtual)
! [output]: OrbAux(N,N)       : MO-->NO coeffs
!           OneRdm(N*(N+1)/2) : triangular 1-RDM in MO 
!
use trexio
implicit none

type(SystemBlock)  :: Mon
type(FlagsData)    :: Flags
integer,intent(in) :: NBasis
double precision   :: OrbAux(NBasis,NBasis), &
                      OneRdm(NBasis*(NBasis+1)/2)

integer    :: i,j,ij,rc,itmp
integer    :: NOccup,HlpDim
integer(8) :: f
double precision :: tmp
double precision :: tmp1,tmp2
double precision,allocatable :: work(:),Eval(:)
double precision,allocatable :: work1d(:)
character(1) :: monlabel

if(mon%monomer==1) monlabel='A'
if(mon%monomer==2) monlabel='B'

HlpDim = max(NBasis**2,3*NBasis)

! allocate occupations, CI coefficients and IGem's
allocate(Mon%CICoef(NBasis),Mon%IGem(NBasis),Mon%Occ(NBasis))

allocate(work(HlpDim),EVal(NBasis))
allocate(work1d(NBasis**2))

f = trexio_open (Mon%TrexFile, 'r', TREXIO_HDF5, rc)

rc = trexio_has_rdm_1e(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No 1-RDM in file'
end if

rc = trexio_read_rdm_1e(f, work1d)

call sq_to_triang(work1d,OneRdm,NBasis)
call Diag8(work1d,NBasis,NBasis,Eval,work)

! sort descending
do i=1,NBasis
   Eval(i) = Abs(Eval(i))
enddo

call SortOcc(EVal,work1d,NBasis)

! in CI all orbitals are active
! they will be truncated later (ThrAct)
Mon%NAct = NBasis

! read NAct from 1RDM
tmp = 0d0
j = 0
do i=1,NBasis
   tmp = tmp + EVal(i) / 2d0
   if((EVal(i)/2d0).lt.Mon%ThrAct) then
      j = j + 1
      if(mon%IPrint.gt.5) write(lout,'(1x,a,i3,e14.4)') 'Warning! Small occupation:',i,Eval(i)
   end if
enddo
if(j.gt.0) write(lout,*) 'Monomoner '//monlabel//': ',j,'orbitals will be treated as virtual'

! Set INAct (also works for open-shells)
Mon%INAct  = Mon%XELE - tmp + 1.d-1
NOccup     = Mon%INAct + Mon%NAct
Mon%SumOcc = tmp + Mon%INAct

Mon%Occ = 0
do i=1,NOccup
   if(i<=Mon%INAct) then
      Mon%Occ(i) = 1d0
   else
      Mon%Occ(i) = EVal(i-Mon%INAct) / 2
   endif
enddo

! construct IGem
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

! switch to MONO order!
ij = 0
do j=1,NBasis
   do i=1,NBasis
      ij = ij + 1
      OrbAux(j,i) = work1d(ij)
   enddo
enddo

rc = trexio_close(f)

deallocate(Eval,work)
deallocate(work1d)

end subroutine readocc_trexio

subroutine readocc_rdm_corr(Mon,NBasis)
!
! set: INActive, NAct
!      NGem, IGem

implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis

integer :: i
integer :: NOccup

Mon%INAct = Mon%num0
Mon%NAct  = Mon%num1
NOccup = Mon%INAct + Mon%NAct

Print*, 'INAct = ',Mon%INAct
Print*, 'NAct = ',Mon%NAct

!do i=1,Mon%NAct
!   print*, 'i,occ,1-occ = ',i,Mon%Occ(i),1d0-Mon%Occ(i)
!enddo

allocate(Mon%IGem0(NBasis))
Mon%IGem0 = Mon%IGem

! this may not work with GVB?
do i=1,Mon%INAct
   Mon%IGem(i) = 1
enddo
do i=Mon%INAct+1,NOccup
   Mon%IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   Mon%IGem(i) = 3
enddo

! construct CICoef
do i=1,NBasis
   Mon%CICoef(i)=sqrt(Mon%Occ(i))
   if(Mon%Occ(i).lt.0.5d0) Mon%CICoef(i)=-Mon%CICoef(i)
enddo

end subroutine readocc_rdm_corr

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

 if(allocated(Mon%Ind2)) deallocate(Mon%Ind2)
 allocate(Mon%Ind2(NBas))

 Mon%Ind2 = Ind2

end subroutine read2rdm

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

subroutine prepare_no_trexio(CMONO,CAOMO,NAO,NBasis)
implicit none
!
! Purpose: get AO-->NO transformation
!
! CMONO[in] :: on input MOtoNO
! CAOMO[in] :: on input AOtoMO
!     [out] :: on output AOtoNO
!
integer,intent(in) :: NAO,NBasis
double precision   :: CMONO(NBasis,NBasis),CAOMO(NAO,NBasis)
double precision   :: work(NAO,NBasis)

! skip canonicalization

 call dgemm('N','N',NAO,NBasis,NBasis,1d0,CAOMO,NAO,CMONO,NBasis,0d0,work,NAO)
 CAOMO = work

end subroutine prepare_no_trexio

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

 !block
 !print*, 'CAOMO: prepare_no',norm2(OrbCAS)
 !do i=1,NBasis
 !   write(lout,'(*(f12.8))') (OrbCAS(i,j),j=1,NBasis)
 !enddo
 !print*,'skip canonicalization...'
 !return
 !end block

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

 call create_ind_molpro(rdmfile,Mon%NumOSym,Mon%IndInt,NSym,NBasis)

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

subroutine rw_trexio_rdm2(Mon,CMONO,NBasis)
!
! Purpose:
! read RDM2 in MO and transform MO2NO (full transformation needed)
! stored in Mon%RDM2val
!
! CAREFUL!!! CURRENTLY NBasis^4 has to fit into memory!
!
!
use trexio
implicit none

type(SystemBlock) :: Mon

integer,intent(in) :: NBasis
double precision,intent(in) :: CMONO(NBasis,NBasis)

integer    :: rc,iunit
integer    :: NOccup
integer    :: NCholesky2RDM
integer(8) :: f, BUFSIZE
integer(8) :: offset,icount
integer    :: i,j,k,l,ij,kl
integer    :: idx_k,idx_l,idx_m,idx_n
double precision :: tmp
integer,allocatable      :: idx_buf(:,:)
double precision,allocatable :: val_buf(:)
double precision,allocatable :: RDM2(:)
double precision,allocatable :: RDM2Chol(:,:,:)
double precision,allocatable :: xnorm
integer :: ichol
double precision             :: tol_chol

f = trexio_open (Mon%TrexFile, 'r', TREXIO_HDF5, rc)

NOccup = NBasis

!print*, 'rw_2rdm_trexio:'
!print*, 'NOccup,NBasis',NOccup,NBasis

allocate(Mon%RDM2val(NOccup,NOccup,NOccup,NOccup))
Mon%RDM2val = 0

if(Mon%Cholesky2rdm) then
   rc = trexio_has_rdm_2e_cholesky(f)
   if(rc /= TREXIO_SUCCESS) then
     stop "No Cholesky-decomposed 2-RDM file!"
   endif
else
   rc = trexio_has_rdm_2e(f)
   if(rc /= TREXIO_SUCCESS) then
     stop "No 2-RDM file!"
   endif
endif

BUFSIZE = NBasis**2

if(Mon%Cholesky2rdm) then
  
   write(lout,'(/,1x,a)') 'WIP1: use 3-ind transformation for 2-RDMs!'
   write(lout,'(1x,a)')   'WIP2: can the matrix be read-in and 3-indx transformed in batches? (rw_trexio_rdm2)'
   allocate(val_buf(BUFSIZE),idx_buf(3,BUFSIZE))
   rc = trexio_read_rdm_2e_cholesky_num(f, NCholesky2rdm)
   write(lout,'(1x,a,i5)') 'Number of Cholesky 2-RDM vectors ', NCholesky2rdm
   print*, '(compared to ',NBasis**2,")"

   allocate(RDM2Chol(NOccup,NOccup,NCholesky2rdm))
   RDM2Chol = 0

   ! read Cholesky-decomposed 2-RDM from TREXIO
   offset  = 0
   icount  = BUFSIZE
   val_buf = 0
   idx_buf = 0
   do while(icount == BUFSIZE)

      rc = trexio_read_rdm_2e_cholesky(f,offset,icount,idx_buf,val_buf)

      do i=1,icount

         idx_k = idx_buf(1,i)
         idx_l = idx_buf(2,i)
         idx_m = idx_buf(3,i)

         RDM2Chol(idx_k,idx_l,idx_m) = val_buf(i)

      enddo

      offset = offset + icount

   enddo

else

   allocate(val_buf(BUFSIZE),idx_buf(4,BUFSIZE))

   ! read full 2-RDM from TREXIO
   offset  = 0
   icount  = BUFSIZE
   val_buf = 0
   idx_buf = 0
   tmp = 0
   xnorm = 0
   do while(icount == BUFSIZE)

      rc = trexio_read_rdm_2e(f,offset,icount,idx_buf,val_buf)

      do i=1,icount

         idx_k = idx_buf(1,i)
         idx_l = idx_buf(2,i)
         idx_m = idx_buf(3,i)
         idx_n = idx_buf(4,i)

         !write(LOUT,'(a,4i3,es15.6)') 'k l m n', idx_k, idx_l, idx_m, idx_n, val_buf(i)
         if (idx_k<=NOccup .and. idx_l<=NOccup .and. idx_m <=NOccup .and. idx_n<=NOccup) then
            Mon%RDM2val(idx_k,idx_m,idx_l,idx_n) = 0.5d0*val_buf(i)
            !Mon%RDM2val(idx_k,idx_l,idx_m,idx_n) = 0.5d0*val_buf(i)
            if(idx_k.eq.idx_m.and.idx_l.eq.idx_n) then
              xnorm = xnorm + 0.5d0*val_buf(i)
            endif
         endif

      enddo

      offset = offset + icount

   enddo

endif
deallocate(val_buf,idx_buf)

! test Cholesky decomposition
if(Mon%Cholesky2rdm) then

  do ichol=1,NCholesky2rdm
     do l=1,NOccup
        do k=1,NOccup
           do j=1,NOccup
              do i=1,NOccup
                 Mon%RDM2val(i,k,j,l) = Mon%RDM2val(i,k,j,l) + RDM2Chol(i,j,ichol)*RDM2Chol(k,l,ichol)
              enddo
           enddo
        enddo
     enddo
  enddo
  Mon%RDM2val = 0.5d0*Mon%RDM2val
  !print*, 'RDM2Chol-norm2 ',norm2(Mon%RDM2val)

  xnorm = 0d0
  do i=1,NOccup
     do j=1,NOccup
        xnorm = xnorm + Mon%RDM2val(i,i,j,j)
     enddo
  enddo

  if(mon%monomer==1) write(lout,'(/1x,a)') 'Monomer A'
  if(mon%monomer==2) write(lout,'(/1x,a)') 'Monomer B'
  write(lout,'(1x,a,f12.6)',advance="no") '2-RDM2 norm = ', xnorm
  write(lout,'(1x,a,f8.3,a)') '(reference =', Mon%XELE*(2d0*Mon%XELE-1), ')'

  ! re-normalize 2RDM
  Mon%RDM2val = Mon%RDM2val * Mon%XELE*(2d0*Mon%XELE-1) / xnorm

endif

call tran_2rdm_trexio(CMONO,Mon%RDM2val,Mon%Occ,   &
                      0,NBasis,NBasis)

rc = trexio_close(f)

! truncate 2-RDM from NBasis^4 to NOccup^4
call truncate_2rdm_trexio(Mon,mon%num0+mon%num1,NBasis)

!print*, 'Gamma-test-NOccup',norm2(Mon%RDM2val)

end subroutine rw_trexio_rdm2

subroutine truncate_2rdm_trexio(Mon,NOccup,NBasis)
!
! truncate 2-RDM from NBasis^4 to NOccup^4
!
implicit none

type(SystemBlock)                :: Mon
integer,intent(in)               :: NOccup,NBasis

integer :: idx_k,idx_m,idx_l,idx_n
double precision,allocatable :: RDM2tru(:,:,:,:)

write(lout,'(1x,a,i1)',advance='no') '2-RDM dim for monomer ', Mon%Monomer
write(lout,'(1x,a,i4,a,i4)') 'truncated from NBasis =', NBasis, ' to NOccup =', NOccup
write(lout,'(1x,a,e13.6,a)') 'according to the', Mon%ThrAct, ' threshold (ThrAct)'

allocate(RDM2tru(NOccup,NOccup,NOccup,NOccup))

RDM2tru = 0d0
do idx_k=1,NOccup
do idx_m=1,NOccup
do idx_l=1,NOccup
do idx_n=1,NOccup
   RDM2tru(idx_k,idx_m,idx_l,idx_n) = Mon%RDM2val(idx_k,idx_m,idx_l,idx_n)
enddo
enddo
enddo
enddo

deallocate(Mon%RDM2val)
allocate(Mon%RDM2val(NOccup,NOccup,NOccup,NOccup))

Mon%RDM2val = RDM2tru

deallocate(RDM2tru)

end subroutine truncate_2rdm_trexio

subroutine tran_2rdm_trexio(CMONO,RDM2val,Occ,INAct,NAct,NBasis)
!
! Purpose: 4-index tran MO2NO of 2-RDM
! Comment: this is now in-core, will be out-of-core!
!
implicit none

integer,intent(in)          :: INAct,NAct,NBasis
double precision,intent(in) :: CMONO(NBasis,NBasis),Occ(NBasis)
double precision,intent(inout)  :: RDM2val(INAct+NAct,INAct+NAct,INAct+NAct,INAct+NAct)

integer :: i,j,k,l
integer :: NOccup,Ind(NBasis)
double precision,allocatable :: work(:,:)

NOccup = INAct + NAct

!print*, 'tran 2rdm:'
!print*, 'INact,NAct',INact,NAct
!print*, 'NOccup',NOccup

Ind = 0
do i=1,NAct
   Ind(INAct+i) =  i
enddo

allocate(work(NBasis,NBasis))

work = 0
!work = transpose(CMONO)
!
!call TrRDM24(RDM2val,work,NOccup,NBasis)
!print*, 'transformed 2-RDM to NO!'
!print*, 'norm-RDM2val',norm2(RDM2val)

work = CMONO
call TrRDM24_dgemm(RDM2val,work,NBasis)
!print*, 'norm-RDM2val',norm2(RDM2val)

deallocate(work)

end subroutine tran_2rdm_trexio

subroutine TrRDM24_dgemm(RDM2,URe,NBasis)
!
!     TRANSFORM RDM2 WITH URe
!
implicit none

integer          :: NBasis
double precision :: URe(NBasis,NBasis),RDM2(NBasis,NBasis,NBasis,NBasis)

double precision :: Aux(NBasis,NBasis,NBasis,NBasis)

write(lout,'(X,"FCI RDM2 TRANSFORMATION TO NO IN PROCESS...")')

call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, RDM2, NBasis, URe, NBasis, 0.d0, Aux, NBasis**3)
call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, Aux,  NBasis, URe, NBasis, 0.d0, RDM2,NBasis**3)
call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, RDM2, NBasis, URe, NBasis, 0.d0, Aux, NBasis**3)
call dgemm('T','N', NBasis**3, NBasis, NBasis, 1.d0, Aux,  NBasis, URe, NBasis, 0.d0, RDM2,NBasis**3)

end subroutine TrRDM24_dgemm

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

end subroutine prepare_rdm2_molpro

subroutine sapt_Kmat_AO(Kmat,CMO,Occ,NAO,NBasis)
implicit none

integer,intent(in) :: NAO, NBasis
double precision,intent(in)    :: Occ(NBasis)
double precision,intent(in)    :: CMO(NBasis,NBasis)
double precision,intent(inout) :: Kmat(NBasis,NBasis)

integer :: i
double precision :: work(NAO,NAO)

work = 0d0
do i=1,NBasis
   call dger(NAO,NAO,Occ(i),CMO(:,i),1,CMO(:,i),1,work,NAO)
enddo

call make_K(NAO,work,Kmat)

print*, 'from sapt_interface',norm2(Kmat)

end subroutine sapt_Kmat_AO

subroutine select_active(mon,nbas,Flags)
!
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
 if(.not.allocated(mon%IndAux)) allocate(mon%IndAux(nbas))

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
          if(mon%Occ(i).lt.mon%ThrGemAct) then
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

! some prints
 print*, 'num0',mon%num0
 print*, 'num1',mon%num1
 print*, 'num2',mon%num2

! active pairs
 if(.not.allocated(mon%IndN)) then
   allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))
 endif

 mon%IPair(1:nbas,1:nbas) = 0

 if(Flags%ICASSCF==0) then
    write(LOUT,'(1x,a,e15.5)') 'Threshold for active orbitals:       ',  mon%ThrSelAct
    write(LOUT,'(1x,a,e15.5)') 'Threshold for quasi-virtual orbitals:',  mon%ThrQVirt
    write(LOUT,'(1x,a,e14.5)') 'Threshold for quasi-inactive orbitals:', mon%ThrQInact
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
                     elseif(abs(mon%Occ(i)+mon%Occ(j)-2d0).gt.mon%ThrQInact) then

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
    write(LOUT,'(1x,a,2e14.5)') 'Threshold for quasi-inactive orbitals:', mon%ThrQInact

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
                     elseif(abs(mon%Occ(i)+mon%Occ(j)-2d0).gt.mon%ThrQInact) then

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

subroutine select_active_thresh(mon,nbas,Flags)
!
! do we really need a separate procedure for TREXIO?
! it would maybe make more sense to set RDMType CI?
!
! Purpose:
! set dimensions: NDimX,num0,num1,num2
! set matrices  : IndN,IndX,IPair,IndAux
!
implicit none

type(SystemBlock)  :: mon
type(FlagsData)    :: Flags
integer,intent(in) :: nbas

integer :: idisPair
integer :: i,j,ij,ind,ind_ij

if(.not.allocated(mon%IndAux)) allocate(mon%IndAux(nbas))

! IndAux = 0 (inactive)
!        = 1 (active)
!        = 2 (virtual)
do i=1,mon%NELE
   mon%IndAux(i) = 0
enddo
do i=1+mon%NELE,nbas
   mon%IndAux(i) = 2
enddo

if(mon%NActOrb/=0) then

   ! select active orbitals based on ThrAct (sets IndAux)
   mon%icnt = 0
   write(LOUT,'()')
   if(mon%Monomer==1) write(LOUT,'(1x,a)') 'Monomer A'
   if(mon%Monomer==2) write(LOUT,'(1x,a)') 'Monomer B'
   do i=1,nbas
      if(abs(1.0d0-mon%Occ(i)).ge.1d-8 .and. mon%Occ(i) .gt. mon%ThrAct) then
      !if(abs(1.0d0-mon%Occ(i)).ge.1d-8 .and. mon%Occ(i) .gt. 1d-10) then
      !if(abs(1.0d0-mon%Occ(i)).ge.1d-8 .and. mon%Occ(i) .gt. 1d-14) then  ! original
      !if(abs(1.0d0-mon%Occ(i)).ge.0d0  .and. mon%Occ(i) .gt. 1d-14) then ! no inactive
         mon%IndAux(i) = 1
         write(6,'(X,"Active Orbital: ",I4,E14.4)') i, mon%Occ(i)
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

if(Mon%IPrint.ge.10) then
  write(lout,'(/1x,a)')     'select_active_thresh:'
  write(lout,'(1x,a,i4)')   'num0 (inactive)',mon%num0
  write(lout,'(1x,a,i4)')   'num1   (active)',mon%num1
  write(lout,'(1x,a,i4,/)') 'num2  (virtual)',mon%num2
endif

! active pairs
if(.not.allocated(mon%IndN)) then
  allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))
endif

mon%IPair(1:nbas,1:nbas) = 0

!print*, 'maybe better call ThrSelAct : Threshold for nearly degenerate pairs?'

write(LOUT,'(/1x,a,e15.5)') 'Threshold for active orbital pairs:       ', mon%ThrSelAct
write(LOUT,'(1x,a,2e15.5)') 'Threshold for quasi-virtual orbital pairs:', mon%ThrQVirt
block
double precision :: ThrInact
ThrInact = 1d0-1d-8
!mon%ThrAct = 1d-10
write(lout,'(1x,a)',advance="no") 'Threshold for selecting active orbitals: '
write(lout,'(1x,e15.5,a,e12.5)') mon%ThrAct, ' < Active < ', ThrInact
end block

if(mon%NCen==1.and.mon%ThrSelAct<1.d-3.and.mon%NAct>1) then
   write(LOUT,'(1x,a)') 'Warning! For single atom ThrSelAct should probably have larger value!'
   mon%IWarn = mon%IWarn + 1
endif

ij  = 0
ind = 0
idisPair = 0
do i=1,nbas
   do j=1,i-1

      ij = ij + 1
      ind_ij = mon%IndAux(i)+mon%IndAux(j)
      if((ind_ij/=0).and.(ind_ij/=4)) then
         ! do not correlate active degenerate orbitals from different geminals
         if((mon%IndAux(i)==1).and.(mon%IndAux(j)==1) &
            .and. (abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.mon%ThrSelAct) ) then
            write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j
         else
            ! if IFlCore=0 exclude core (inactive) orbitals
            if(Flags%IFlCore==1.or.&
                 (Flags%IFlCore==0.and.&
                 mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
                 ! exclude pairs of nearly/virtual orbitals
                 if(abs(mon%Occ(i)+mon%Occ(j)).lt.mon%ThrQVirt) then
                    idisPair = idisPair + 1
                    if(mon%IPrint>10) write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly virtual-orbitals pair',i,j
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
write(LOUT,*) 'Number of discarded nearly degenerate virtual-orbital pairs ', idisPair

mon%NDimX = ind

end subroutine select_active_thresh

subroutine save_CAONO(Cin,Cout,NAO,NBasis)
!
! C(NAO*NBasis) --> CAONO(NAO,NBasis)
!
implicit none

integer,intent(in) :: NAO,NBasis
double precision,intent(in)  :: Cin(:)
double precision,intent(out) :: Cout(NAO,NBasis)

integer :: i,j,ij

Cout = 0
ij = 0
do j=1,NBasis
   do i=1,NAO
      ij = ij + 1
      Cout(i,j) = Cin(ij)
   enddo
enddo

end subroutine save_CAONO

subroutine calc_elpot(A,B,CholeskyVecs,ICholesky,NAO,NBas)
implicit none

type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: ICholesky,NAO,NBas

integer :: ione,i
integer :: NInte1,NCholesky
double precision,allocatable :: Pa(:,:),Pb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Ja(:,:),Jb(:,:)
logical                      :: valid
character(8)                 :: label

 NInte1 = NBas*(NBas+1)/2

 allocate(Pa(NAO,NAO),Pb(NAO,NAO),&
          Va(NAO,NAO),Vb(NAO,NAO),&
          Ja(NAO,NAO),Jb(NAO,NAO))

 !call get_den(NAO,NBas,A%CMO,A%Occ,2d0,Pa)
 !call get_den(NAO,NBas,B%CMO,B%Occ,2d0,Pb)
 Pa = 0d0; Pb = 0d0
 do i=1,NBas
    call dger(NAO,NAO,2d0*A%Occ(i),A%CMO(:,i),1,A%CMO(:,i),1,Pa,NAO)
    call dger(NAO,NAO,2d0*B%Occ(i),B%CMO(:,i),1,B%CMO(:,i),1,Pb,NAO)
 enddo

 !call get_one_mat('V',Va,A%Monomer,NBas)
 !call get_one_mat('V',Vb,B%Monomer,NBas)
 valid=.false.
 Va = 0d0
 open(newunit=ione,file='ONEEL_A',access='sequential',&
      form='unformatted',status='old')
    read(ione)
    read(ione) label,Va
    if(label=='POTENTAL') valid=.true.
 close(ione)
 if(.not.valid) then
    write(LOUT,'(1x,a)') 'Va not found in calc_elpot!'
 endif
 valid=.false.
 Vb = 0d0
 open(newunit=ione,file='ONEEL_B',access='sequential',&
      form='unformatted',status='old')
    read(ione)
    read(ione) label,Vb
    if(label=='POTENTAL') valid=.true.
 close(ione)
 if(.not.valid) then
    write(LOUT,'(1x,a)') 'Vb not found in calc_elpot!'
 endif

 if(ICholesky==0) then
    call make_J2(NAO,Pa,Pb,Ja,Jb)
 elseif(ICholesky==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_J2_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                       Pa,Pb,Ja,Jb,NCholesky,NAO)
 endif

 allocate(A%WPot(NAO,NAO),B%WPot(NAO,NAO))

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

 if(SAPT%SaptLevel==666) then ! RS2PT2+
    allocate(A%OOAB(NCholesky,dimOA*dimOB), &
             B%OOBA(NCholesky,dimOB*dimOA))

    call chol_MOTransf_TwoStep(A%OOAB,CholeskyVecs,&
                       A%CMO,1,dimOA,&
                       B%CMO,1,dimOB,&
                       MaxBufferDimMB)

    call chol_MOTransf_TwoStep(B%OOBA,CholeskyVecs,&
                       B%CMO,1,dimOB,&
                       A%CMO,1,dimOA,&
                       MaxBufferDimMB)
 endif

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

subroutine  square_oneint(tr,sq,nbas,nsym,norb)

implicit none
integer,intent(in) :: nbas,nsym,norb(8)
double precision,intent(in) :: tr(:)
double precision,intent(out) :: sq(nbas,nbas)
integer :: irep,i,j
integer :: offset,idx

sq=0

offset=0
idx=0
do irep=1,nsym
   do j=offset+1,offset+norb(irep)
      do i=offset+1,j

         idx=idx+1
         sq(i,j)=tr(idx)
         sq(j,i)=tr(idx)

      enddo
   enddo
   offset=offset+norb(irep)
enddo

end subroutine square_oneint

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
    write(LOUT,'(1x,a,1x,i6,8x,i6)') 'Total number of pairs: ', SAPT%monA%NDim,SAPT%monB%NDim
    write(LOUT,'(1x,a,12x,i6,8x,i6)') 'Reduced to: ', SAPT%monA%NDimX, SAPT%monB%NDimX
    if (SAPT%monA%NDimX0/=0 .or. SAPT%monB%NDimX0/=0) then
       write(LOUT,'(1x,a,6x,i6,8x,i6,a)') '(without RDM corr:', SAPT%monA%NDimX0, SAPT%monB%NDimX0, ')'
    endif
    write(LOUT,'()')

    if(SAPT%IPrint.gt.10) then
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


subroutine compare_active(Mon, NBas)
!
! compare the current set of occupation numbers (Occ)
! with the initial one (written for RDMCORR = .TRUE. option)
!
implicit none
type(SystemBlock)  :: Mon
integer,intent(in) :: NBas

integer :: i

write(LOUT,'(1x,a)') 'ORBITAL OCCUPANCIES'
write(LOUT,'(1x,a,3x,a,4x,a,10x,a,6x,a)') 'CASSCF', 'Occ-Initial', 'Gem-I', 'Occ-Final  ','Gem-F'
do i=1,nbas
   write(LOUT,'(1x,i3,1x,e16.6,1x,i6,7x,e16.6,3x,i6)') i, Mon%Occ0(i),Mon%IGem0(i),Mon%Occ(i),Mon%IGem(i)
enddo
write(LOUT, '()')

end subroutine compare_active

subroutine print_mo(cmo,n,mon)
implicit none

integer,intent(in) :: n
double precision,intent(in) :: cmo(n,n)
character(*) :: mon
integer :: i,j,ll,nn
integer :: nline

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'NATURAL ORBITALS '//mon
 do i=1,n
    write(LOUT,'(1x,i3)') i
    write(LOUT,'(10f10.6)') cmo(:,i)
    write(LOUT,'()')
 enddo

end subroutine print_mo

subroutine print_TwoInt(NBasis)
! Purpose: for debugging,
!          print trasformed integrals
!
implicit none

integer :: NBasis
integer :: ip,iq,ir,is,irs,ipq
integer :: iunit,i
double precision :: work1(NBasis*NBasis)
double precision :: work2(NBasis*NBasis)


 open(newunit=iunit,file='TWOMOAB',status='OLD', &
      access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'Two-electron integrals in the NO representation:'
 write(LOUT,'(4x,a,12x,a)') 'p   q   r   s', 'Val'
 write(LOUT,'(1x,8a6)') ('------',i=1,8)
 irs=0
 do is=1,NBasis
    do ir=1,is
       irs=irs+1
       read(iunit,rec=irs) work1(1:NBasis*(NBasis+1)/2)
       ipq=0
       do iq=1,NBasis
          do ip=1,iq
             ipq = ipq+1
             write(LOUT,'(1x,4i4,3x,f20.16)') ip,iq,ir,is,work1(ipq)
          enddo
       enddo
    enddo
 enddo

 close(iunit)

end subroutine print_TwoInt

subroutine print_sqmat(mat,ndim)
!
! Print square matrix
!
implicit none

integer,intent(in) :: ndim
double precision,intent(in) :: mat(ndim,ndim)
integer :: i,j

 do i=1,ndim
    write(LOUT,*) i
    write(LOUT,'(10f13.8)') (mat(i,j),j=1,ndim)
 enddo
 write(LOUT,'()')

 return
end subroutine print_sqmat

subroutine print_diag(mat,ndim)
!
! Print diagonal of a square matrix
!
implicit none

integer,intent(in) :: ndim
double precision,intent(in) :: mat(ndim,ndim)
integer :: i

 do i=1,ndim
    write(LOUT,'(10f11.6)') mat(i,i)
 enddo
 write(LOUT,'()')

 return
end subroutine print_diag

end module sapt_inter
