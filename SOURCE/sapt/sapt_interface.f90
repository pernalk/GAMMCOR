module sapt_inter

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
use sapt_interface_chol
use sapt_interface_io

implicit none

contains

subroutine sapt_interface(Flags,SAPT,NBasis,AOBasis,CholeskyVecsOTF)
!
! Possible interface: Dalton, Molpro
!

! Read 1-RDM, get C(AO,NO) and Occ from diagonalization
!     -- canonicalize orbitals, if needed
! Read 2-RDM (active), transform to NO
! Read / Generate 1-, 2-el integrals
!     -- For CBS[H], create local mu(r)
! Transform integrals to NOs
! Calculate elst potential in AO: W=J+K 
! Calculate K[PB] matrix in AO
!
! Comments :
! SAPT-DALTON requires SIRIFC and SIRIUS.RST
! or SIRIFC and occupations.dat
!
implicit none

type(FlagsData)        :: Flags
type(SaptData)         :: SAPT
type(TCholeskyVecs)    :: CholeskyVecs
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(TAOBasis)         :: AOBasis
type(TSystem)          :: System
integer,intent(in)     :: NBasis

integer    :: NSq,NInte1,NInte2
integer    :: NCholesky
integer    :: NCholeskyTHC,NGridTHC
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

character(:),allocatable :: XYZPath
character(:),allocatable :: BasisSetPath

double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: work(:,:)
double precision,allocatable :: Ha(:),Hb(:)
double precision,allocatable :: Va(:),Vb(:),S(:)
double precision,allocatable :: Ca(:),Cb(:)
double precision,allocatable :: AuxA(:,:),AuxB(:,:)
double precision,allocatable :: OneRdmA(:),OneRdmB(:)

! temporary solution, maybe pass Xgp, Zgk through a type
double precision,allocatable :: Xgp(:,:), Zgk(:,:)

logical :: SortAngularMomenta
character(:),allocatable :: BasisSet

logical :: doRSH
logical :: canoni
double precision,allocatable :: Sa(:,:),Sb(:,:)
double precision :: Tcpu,Twall

! set monomer print level
SAPT%monA%IPrint = SAPT%IPrint
SAPT%monB%IPrint = SAPT%IPrint

! set basis set
if (allocated(Flags%BasisSetPath)) then
   write(LOUT,'(/1x,"Flags:BasisSetPath ",a)') Flags%BasisSetPath
   write(LOUT,'(1x, "Flags:BasisSet ",a)')     Flags%BasisSet
   BasisSet = Flags%BasisSetPath // Flags%BasisSet
else
   BasisSet = "Empty"
endif

! set dimensions
NSq = NBasis**2
NInte1 = NBasis*(NBasis+1)/2
NInte2 = NInte1*(NInte1+1)/2
SAPT%monA%NDim = NBasis*(NBasis-1)/2
SAPT%monB%NDim = NBasis*(NBasis-1)/2

! sanity-check orbital ordering
 call check_orbital_ordering(Flags%ICholeskyOTF)

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
    call onel_molpro(SAPT%monA%Monomer,NBasis,SAPT%monA,SAPT)
    call onel_molpro(SAPT%monB%Monomer,NBasis,SAPT%monB,SAPT)
 endif

 SAPT%NAO = NBasis

 if(SAPT%InterfaceType==1) then
 ! (Dalton) read SR Coulomb and V_KS potential (in AO)
    if(doRSH) then
       ! maybe it would be better to calculate Jsr in our code?
       allocate(SAPT%monA%VsrKS(NBasis,NBasis),SAPT%monA%Jsr(NBasis,NBasis))
       call read_vKS_dalton(SAPT%monA%VsrKS,'dftSRfile_A.dat',NBasis)
       call read_Jsr_dalton(SAPT%monA%Jsr  ,'dftSRfile_A.dat',NBasis)
       call read_esrDFT_dalton(SAPT%monA%esrDFT,'dftSRfile_A.dat')

       write(LOUT,'(/1x,a)') 'SR Kohn-Sham potential read from dftSRfile_A.dat'
       write(LOUT,'(1x,a)')  'SR Coulomb integrals   read from dftSRfile_A.dat'

       allocate(SAPT%monB%VsrKS(NBasis,NBasis),SAPT%monB%Jsr(NBasis,NBasis))
       call read_vKS_dalton(SAPT%monB%VsrKS,'dftSRfile_B.dat',NBasis)
       call read_Jsr_dalton(SAPT%monB%Jsr  ,'dftSRfile_B.dat',NBasis)

       call read_esrDFT_dalton(SAPT%monB%esrDFT,'dftSRfile_B.dat')

       call arrange_oneint(SAPT%monB%VsrKS,NBasis,SAPT)
       call arrange_oneint(SAPT%monB%Jsr,NBasis,SAPT)

       write(LOUT,'(/1x,a)') 'SR Kohn-Sham potential read from dftSRfile_B.dat'
       write(LOUT,'(1x,a)')  'SR Coulomb integrals   read from dftSRfile_B.dat'

       !print*, 'JSR from Dalton'
       !do j=1,NBasis
       !   write(6,'(*(f13.8))') (Jsr(i,j),i=1,NBasis)
       !enddo
       !print*, 'VsrKS from Dalton'
       !do j=1,NBasis
       !   write(6,'(*(f13.8))') (VsrKS(i,j),i=1,NBasis)
       !enddo
    endif
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
    !if (allocated(SAPT%monA%dipm)) print*, 'DIP-A allocated in sapt_interface!'
    !if (allocated(SAPT%monB%dipm)) print*, 'DIP-B allocated in sapt_interface!'
      allocate(SAPT%monA%dipm(3,NBasis,NBasis))
      allocate(SAPT%monB%dipm(3,NBasis,NBasis))
    associate (DipA => SAPT%monA%dipm, &
               DipB => SAPT%monB%dipm)
      call read_dip_molpro(DipA(1,:,:),DipA(2,:,:),DipA(3,:,:),'DIP_A',NBasis)
      call read_dip_molpro(DipB(1,:,:),DipB(2,:,:),DipB(3,:,:),'DIP_B',NBasis)
    end associate
 endif

! read 2-el integrals
 call clock('START',Tcpu,Twall)

! memory allocation for sorter
 MemSrtSize = Flags%MemVal*1024_8**Flags%MemType

! for testing old Cholesky
 ICholOld = 0

 if(Flags%ICholesky==0.or.ICholOld==1) then

    if(SAPT%InterfaceType==1) then ! Dalton

       call readtwoint(NBasis,1,'AOTWOINT_A','AOTWOSORT',MemSrtSize)
       if(doRSH) then

         if (SAPT%SameOm) then
            call readtwoint(NBasis,1,'AOERFINT_A','AOERFSORT',MemSrtSize)
         else
            call readtwoint(NBasis,1,'AOERFINT_A','AOERFSORT',MemSrtSize)
            call readtwoint(NBasis,1,'AOERFINT_B','AOERFSORTB',MemSrtSize)
         endif

         if(SAPT%IPrint.gt.1) then
            write(lout,'(/1x,a)') "Sorting AO integrals DALTON interface:"
            write(lout,'(1x,a)') "AOTWOSORT should contain full-range integrals"
            write(lout,'(1x,a)') "AOERFSORT should contain LR integrals"
            !write(lout,'(1x,a)') "AOSR2SORT should contain SR integrals"
         endif
       endif

    elseif(SAPT%InterfaceType==2) then ! Molpro

       call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT',MemSrtSize)
       if(doRSH) then
          if(SAPT%SameOm) then
             call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
          else
             call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT',MemSrtSize)
             call readtwoint(NBasis,2,'AOTWOINT.erfB','AOERFSORTB',MemSrtSize)
          endif

          if(SAPT%IPrint.gt.1) then
             write(lout,'(1x,a)') "Sorting AO integrals Molpro interface:"
             write(lout,'(1x,a)') "AOTWOSORT should contain full-range integrals"
             write(lout,'(1x,a)') "AOERFSORT should contain LR integrals"
          endif
       endif ! doRSH
    endif ! Interface
 endif ! not Cholesky

! Cholesky decomposition
 if(Flags%ICholeskyBIN==1.or.Flags%ICholeskyOTF==1) then

    !! old Cholesky
    !if(ICholOLD==1) print*, 'old Cholesky transformation...'
    !if(ICholOLD==1) call chol_CoulombMatrix(CholeskyVecs,'AOTWOSORT',Flags%ICholeskyAccu)

    ! Cholesky binary
    if(ICholOLD==0) print*, 'new Cholesky transformation...'

    if(Flags%ICholeskyBIN==1) then

       write(lout,'(/1x,3a6)') ('******',i=1,3)
       write(lout,'(1x,a)') 'Cholesky Binary'
       write(lout,'(1x,3a6)') ('******',i=1,3)

       if(SAPT%InterfaceType==1) then
          call chol_CoulombMatrix(CholeskyVecs,NBasis,'AOTWOINT_A',1,Flags%ICholeskyAccu)
       elseif(SAPT%InterfaceType==2) then
          call chol_CoulombMatrix(CholeskyVecs,NBasis,'AOTWOINT.mol',2,Flags%ICholeskyAccu)
       endif

       SAPT%NCholesky  = CholeskyVecs%NCholesky
       SAPT%monA%NChol = SAPT%NCholesky
       SAPT%monB%NChol = SAPT%NCholesky

    ! Cholesky on-the-fly
    elseif(Flags%ICholeskyOTF==1) then

       write(lout,'(/1x,3a6)') ('******',i=1,3)
       write(lout,'(1x,a)') 'Cholesky On-The-Fly'
       write(lout,'(1x,3a6)') ('******',i=1,3)

       call auto2e_init()

       XYZPath = "./input.inp"
       BasisSetPath = BasisSet
       SortAngularMomenta = .true.

       call CholeskyOTF_ao_vecs(CholeskyVecsOTF,AOBasis,System,Flags%IUnits, &
                                 XYZPath,BasisSetPath, &
                                 SortAngularMomenta,Flags%ICholeskyAccu)

       SAPT%NCholesky  = CholeskyVecsOTF%Chol2Data%NVecs
       SAPT%monA%NChol = SAPT%NCholesky
       SAPT%monB%NChol = SAPT%NCholesky
       ! set THC for FockOTF
       NGridTHC = 1
       NCholeskyTHC=1

    elseif(Flags%ICholeskyTHC==1) then

       stop "SAPT not ready with THC!"

    endif ! CholeskyOTF

 endif
 call clock('2ints',Tcpu,Twall)

 if(SAPT%InterfaceType==2) then

    if(SAPT%monA%NatOrb==0) then
       ! create NOs inside GammCor (use canonical CAS orbs)
       call prepare_no_molpro(Ca,OneRdmA,AuxA,SAPT%monA,AOBasis,System, &
                       CholeskyVecs,CholeskyVecsOTF,  &
                       Xgp,Zgk,NGridTHC,NCholeskyTHC, &
                       Flags,NBasis)
    elseif(SAPT%monA%NatOrb==1) then
       print*, 'MONOMER A: use Natural Orbitals from Molpro'
       block
       integer :: ione
       character(8) :: label
       double precision :: CSAOMO(NBasis,NBasis)
       double precision :: SAO(NBasis,NBasis),ttt(NBasis,NBasis)
       CSAOMO = 0d0
       ij = 0
       do j=1,NBasis
          do i=1,NBasis
             ij = ij + 1
             CSAOMO(i,j) = Ca(ij)
          enddo
       enddo
       ! get S in AO
       open(newunit=ione,file='ONEEL_A',access='sequential',&
            form='unformatted',status='old')
       read(ione) label, SAO
       close(ione)


       call read_no_molpro(Ca,SAPT%monA%InSt(1,1),'MOLPRO_A.MOPUN','NATORB  ',NBasis)

       call dgemm('T','N',NBasis,NBasis,NBasis,1d0,CSAOMO,NBasis,SAO,NBasis,0d0,ttt,NBasis)
       call dgemm('N','N',NBasis,NBasis,NBasis,1d0,ttt,NBasis,Ca,NBasis,0d0,AuxA,NBasis)
       end block
    else
       stop "Wrong NatOrb Value!"
    endif

    if(SAPT%monB%NatOrb==0) then
       call prepare_no_molpro(Cb,OneRdmB,AuxB,SAPT%monB,AOBasis,System, &
                       CholeskyVecs,CholeskyVecsOTF,  &
                       Xgp,Zgk,NGridTHC,NCholeskyTHC, &
                       Flags,NBasis)
    elseif(SAPT%monB%NatOrb==1) then

       print*, 'MONOMER B: use Natural Orbitals from Molpro'
       block
       double precision :: CSAOMO(NBasis,NBasis)
       double precision :: SAO(NBasis,NBasis)
       ij = 0
       do j=1,NBasis
          do i=1,NBasis
             ij = ij + 1
             CSAOMO(i,j) = Cb(ij)
          enddo
       enddo

       call read_no_molpro(Cb,SAPT%monB%InSt(1,1),'MOLPRO_B.MOPUN','NATORB  ',NBasis)
       call dgemm('T','N',NBasis,NBasis,NBasis,1d0,CSAOMO,NBasis,SAO,NBasis,0d0,work,NBasis)
       call dgemm('N','N',NBasis,NBasis,NBasis,1d0,work,NBasis,Cb,NBasis,0d0,AuxB,NBasis)
       end block

    else
       stop "Wrong NatOrb Value!"
    endif
!
!       print*, 'Skipping canonicalization...'
!       call prepare_no_molpro_skip(AuxA,Ca,SAPT%monA%INAct,SAPT%monA%NAct,NBasis,NBasis)
!       call prepare_no_molpro_skip(AuxB,Cb,SAPT%monB%INAct,SAPT%monB%NAct,NBasis,NBasis)

       !print*, 'aaaa: AuxA'
       !do j=1,NBasis
       !   write(6,'(*(f12.6))') (AuxA(i,j),i=1,NBasis)
       !enddo

    call prepare_rdm2_molpro(SAPT%monA,AuxA,NBasis)
    call prepare_rdm2_molpro(SAPT%monB,AuxB,NBasis)
 endif

 ! create approximate 2-rdm for SAPT(DMFT)
 if(SAPT%SaptExch==1.and.Flags%IRDM2Typ/=1) call prepare_rdm2_approx(SAPT%monA,Flags%IRDM2Typ,NBasis)
 if(SAPT%SaptExch==1.and.Flags%IRDM2Typ/=1) call prepare_rdm2_approx(SAPT%monB,Flags%IRDM2Typ,NBasis)

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

 !! test
 !write(LOUT,*) 'sapt_interface: CNO'
 !print*, norm2(SAPT%monA%CMO)
 !do j=1,NBasis
 !   print*, j
 !   write(*,'(14f11.6)') (SAPT%monA%CMO(i,j),i=1,nbasis)
 !end do

 if(SAPT%InterfaceType==1.and.Flags%ICholeskyOTF==1) then
    allocate(SAPT%monA%CAONO(NBasis,NBasis),SAPT%monB%CAONO(NBasis,NBasis))
    SAPT%monA%CAONO = SAPT%monA%CMO
    SAPT%monB%CAONO = SAPT%monB%CMO
 endif

! look-up tables
 call select_active(SAPT%monA,NBasis,Flags)
 call select_active(SAPT%monB,NBasis,Flags)

 !print*, 'NACT...'
 !SAPT%monA%NAct=SAPT%monA%num1
 !SAPT%monB%NAct=SAPT%monB%num1

 ! transform Cholesky Vecs to NO
 if(Flags%ICholeskyBIN==1) then
    !call chol_sapt_AO2NO_BIN(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call chol_OO_sapt_AO2NO_BIN(SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call chol_FO_sapt_AO2NO_BIN(SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call chol_FF_sapt_AO2NO_BIN(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call clock('chol_AO2NO_BIN',Tcpu,Twall)
 elseif(Flags%ICholeskyOTF==1) then
    !call chol_sapt_AO2NO_OTF(SAPT,SAPT%monA,SAPT%monB,CholeskyVecsOTF,AOBasis,Flags,NBasis)
    call chol_OO_sapt_AO2NO_OTF(SAPT,SAPT%monA,SAPT%monB,CholeskyVecsOTF,AOBasis,Flags,NBasis)
    call chol_FO_sapt_AO2NO_OTF(SAPT,SAPT%monA,SAPT%monB,CholeskyVecsOTF,AOBasis,Flags,NBasis)
    call clock('chol_AO2NO_OTF',Tcpu,Twall)
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
 if(.not.allocated(SAPT%monB%Kmat)) allocate(SAPT%monB%Kmat(NBasis,NBasis))
 allocate(work(NBasis,NBasis))
 !get PB density in AO
 work = 0
 do i=1,NBasis
 call dger(NBasis,NBasis,SAPT%monB%Occ(i),SAPT%monB%CMO(:,i),1,SAPT%monB%CMO(:,i),1,work,NBasis)
 enddo

 if(Flags%ICholesky==0) then
    call make_K(NBasis,work,SAPT%monB%Kmat,'AOTWOSORT')
 elseif(Flags%ICholeskyBIN==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_K_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                      NCholesky,NBasis,work,SAPT%monB%Kmat)
 endif
 deallocate(work)

! calculate electrostatic potential: W = V + J (in AO)
 if(SAPT%InterfaceType==1.and.Flags%ICholeskyOTF==1) then
    ! Dalton/CholeskyOTF: J,K,W in sapt_mon_ints
 else
    call calc_elpot(SAPT%monA,SAPT%monB,CholeskyVecs,&
                    Flags%ICholesky,Flags%ICholeskyBIN,Flags%ICholeskyOTF,NBasis)
 endif

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ca,Cb)
 if(SAPT%InterfaceType==2) then
    deallocate(OneRdmB,OneRdmA,AuxB,AuxA)
 endif

end subroutine sapt_interface

subroutine saptuks_interface(Flags,SAPT,NBasis,AOBasis,CholeskyVecsOTF)
!
! SAPT(UKS) interface :
!     -- works only with MOLPRO
!     -- reads 1-, 2-el integrals (AO)
!     -- reads C_alpha, C_beta (AO,MO) orbs
!     -- reads orbital energies and occupation numbers
!
implicit none

type(FlagsData)        :: Flags
type(SaptData)         :: SAPT
type(TCholeskyVecs)    :: CholeskyVecs
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(TAOBasis)         :: AOBasis
type(TSystem)          :: System

integer,intent(in)     :: NBasis

integer(8) :: MemSrtSize

character(:),allocatable :: XYZPath
character(:),allocatable :: BasisSetPath
character(:),allocatable :: BasisSet

double precision :: Tcpu,Twall

! check interface
if (SAPT%InterfaceType/=2) stop "Unrestricted SAPT only works with Molpro"

! set basis set
!write(LOUT,'(/1x,"Flags:BasisSetPath ",a)') Flags%BasisSetPath
!write(LOUT,'(1x, "Flags:BasisSet ",a)')     Flags%BasisSet
!BasisSet = Flags%BasisSetPath // Flags%BasisSet

! check where molpro keeps nao
SAPT%NAO = NBasis

! read and dump 1-electron integrals
call onel_molpro(SAPT%monA%Monomer,NBasis,SAPT%monA,SAPT)
call onel_molpro(SAPT%monB%Monomer,NBasis,SAPT%monB,SAPT)

! read unrestricted orbitals
allocate(SAPT%monA%UMO(NBasis,NBasis,2),SAPT%monB%UMO(NBasis,NBasis,2))
call read_umo_molpro(SAPT%monA%UMO,NBasis,'UKSORB  ','MOLPRO_A.MOPUN')
call read_umo_molpro(SAPT%monB%UMO,NBasis,'UKSORB  ','MOLPRO_B.MOPUN')

! read unrestricted occupation numbers
allocate(SAPT%monA%UOcc(NBasis,2),SAPT%monB%UOcc(NBasis,2))
call read_uocc_molpro(SAPT%monA%UOcc,NBasis,'UKSORB  ','MOLPRO_A.MOPUN')
call read_uocc_molpro(SAPT%monB%UOcc,NBasis,'UKSORB  ','MOLPRO_B.MOPUN')

! read unrestricted orbital energies
allocate(SAPT%monA%UOrbE(NBasis,2),SAPT%monB%UOrbE(NBasis,2))
call read_uorbe_molpro(SAPT%monA%UOrbE,NBasis,'UKSORB  ','MOLPRO_A.MOPUN')
call read_uorbe_molpro(SAPT%monB%UOrbE,NBasis,'UKSORB  ','MOLPRO_B.MOPUN')

! unpack to C1 sym
call unpack_uks_sym(SAPT%monA,NBasis)
call unpack_uks_sym(SAPT%monB,NBasis)

! look-up tables
! set unrestricted occ, virt, ov, IndN, ...
call select_uactive(SAPT%monA,SAPT%monB,NBasis)

call print_uocc(NBasis,SAPT)

! read 2-el integrals
call clock('START',Tcpu,Twall)

! memory allocation for sorter
MemSrtSize = Flags%MemVal*1024_8**Flags%MemType

! Cholesky decomposition
if(Flags%ICholeskyBIN==1.or.Flags%ICholeskyOTF==1) then
   stop "SAPT(UKS) with Cholesky not ready yet..."
else
   call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT',MemSrtSize)
endif

!! no need for canonicalization...
!! set iPINO
!SAPT%iPINO = 0

print*, 'saptuks_interface: Skipping K[PB] at this point...'

!! calculate exchange K[PB] matrix in AO
!if(.not.allocated(SAPT%monB%Kmat)) allocate(SAPT%monB%Kmat(NBasis,NBasis))
!allocate(work(NBasis,NBasis))
!!get PB density in AO
!work = 0
!do i=1,NBasis
!   call dger(NBasis,NBasis,SAPT%monB%Occ(i),SAPT%monB%CMO(:,i),1,SAPT%monB%CMO(:,i),1,work,NBasis)
!enddo
!
!if(Flags%ICholesky==0) then
!   call make_K(NBasis,work,SAPT%monB%Kmat)
!else
!  stop "SAPT(UKS) with Cholesky not ready yet..."
!endif

! ABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB

! calculate electrostatic potential: W = V + J (in AO)
call calc_uks_elpot(SAPT%monA,CholeskyVecs,&
                    Flags%ICholesky,Flags%ICholeskyBIN,Flags%ICholeskyOTF,NBasis)
call calc_uks_elpot(SAPT%monB,CholeskyVecs,&
                    Flags%ICholesky,Flags%ICholeskyBIN,Flags%ICholeskyOTF,NBasis)
!
! calc intermolecular repulsion
SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

end subroutine saptuks_interface

subroutine unpack_uks_sym(Mon,NBasis)
!
! destroy symmetry in UHF/UHF occupations, orb. energies and orbitals
! Molpro keeps occ-virt in each irrep
! we reorder to occ1-occ2-...-virt1-...-virtN
!
! Ca = C(SAO,MO) alpha
! Cb = C(SAO,MO) beta

implicit none

type(SystemBlock) :: Mon
integer,intent(in) :: NBasis

integer :: i,j
integer :: NSym
real(8) :: UOcca(NBasis),UOccb(NBasis)
real(8) :: UOrbEa(NBasis),UOrbEb(NBasis)
real(8) :: Ca(NBasis,NBasis),Cb(NBasis,NBasis)

integer :: IndIntA(NBasis),IndIntB(NBasis)
integer :: NumOSymA(15),NumOSymB(15)
character(:),allocatable :: basinfile

if (Mon%Monomer==1) then
   basinfile = 'BASINFA'
elseif (Mon%Monomer==2) then
   basinfile = 'BASINFB'
endif

allocate(Mon%NumOSym(15),Mon%IndInt(NBasis))
! alpha and beta
call create_ind_uks_molpro('A',basinfile,NumOSymA,IndIntA,NSym,NBasis)
call create_ind_uks_molpro('B',basinfile,NumOSymB,IndIntB,NSym,NBasis)

!print*, 'Monomer= ', Mon%Monomer
!print*, 'IndInt = '
!do i=1,NBasis
!  write(lout,'(1x,3i3)') i,IndIntA(i),IndIntB(i)
!enddo

! reorder MO to no symmetry
do i=1,NBasis
   do j=1,NBasis
      Ca(IndIntA(i),j) = mon%UMO(j,i,1)
      Cb(IndIntB(i),j) = mon%UMO(j,i,2)
   enddo
enddo

! reorder Occ to no symmetry
do i=1,NBasis
   UOrbEa(IndInta(i)) = mon%UOrbE(i,1)
   UOrbEb(IndIntb(i)) = mon%UOrbE(i,2)
   UOcca(IndInta(i))  = mon%UOcc(i,1)
   UOccb(IndIntb(i))  = mon%UOcc(i,2)
enddo

!print*,' UOcc  alpha beta'
!do i=1,NBasis
!   write(6,'(i3,2f12.8)') i, UOcca(i), UOccb(i)
!enddo
!
!print*,' OrbEne  alpha beta'
!do i=1,NBasis
!   write(6,'(i3,2f12.8)') i, UOrbEa(i), UOrbEb(i)
!enddo
!
!print*, 'Monomer = ', Mon%Monomer
!Print*, 'CSAOMO-alpha sym unpacked =',norm2(Ca)
!do i=1,NBasis
!   write(6,'(*(f13.8))') (Ca(i,j),j=1,NBasis)
!enddo
!Print*, 'CSAOMO-beta  sym unpacked =',norm2(Cb)
!do i=1,NBasis
!   write(6,'(*(f13.8))') (Cb(i,j),j=1,NBasis)
!enddo

! rewrite
Mon%UOcc(:,1)=UOcca
Mon%UOcc(:,2)=UOccb
!
Mon%UOrbE(:,1)=UOrbEa
Mon%UOrbE(:,2)=UOrbEb
!
mon%UMO = 0
do i=1,NBasis
   do j=1,NBasis
      mon%UMO(i,j,1)=Ca(j,i)
      mon%UMO(i,j,2)=Cb(j,i)
   enddo
enddo

end subroutine unpack_uks_sym

subroutine sapt_interface_spin(Flags,SAPT,NBasis)
!
! Interface for open-shell SAPT (SAPT-OS JobType)
!
! Purpose: construct alpha/beta spin densities in NOs
! from charge/spin densities
!
! 1-RDM in NOs:
! 1/2 * Gamma_{pq} = 1/2 * (GammaChrg^\alpha_{pq} + GammaChrg^\beta_{pq} ) = n_p \delta_pq
!
! charge densities are available as Occ(NBasis) (read in sapt_interface)
!
implicit none

type(FlagsData)     :: Flags
type(SaptData)      :: SAPT
integer,intent(in)  :: NBasis

integer             :: NActA,NActB
integer             :: INActA,INActB
integer             :: i,j
double precision,allocatable :: GChrgA(:,:),GChrgB(:,:)
double precision,allocatable :: GSpinA(:,:),GSpinB(:,:)
double precision,allocatable :: GAAct(:,:),GBAct(:,:)

! dimensions
NActA  = SAPT%monA%NAct
NActB  = SAPT%monB%NAct
INActA = SAPT%monA%INAct
INActB = SAPT%monB%INAct

allocate(GAAct(NActA,NActA),GBAct(NActB,NActB))
allocate(GChrgA(NBasis,NBasis),GChrgB(NBasis,NBasis))
allocate(GSpinA(NBasis,NBasis),GSpinB(NBasis,NBasis))

! charge densities
print*, 'NASHT-A',SAPT%monA%NAct
print*, 'NISHT-A',SAPT%monA%INAct

GChrgA = 0d0
GChrgB = 0d0
do i=1,NBasis
   GChrgA(i,i) = 2.0d0*SAPT%monA%Occ(i)
   GChrgB(i,i) = 2.0d0*SAPT%monB%Occ(i)
enddo

! spin densities
GSpinA = 0d0
GSpinB = 0d0

! active blocks
call read_1rdm_spin_dalton(GAAct,'rdms1_A.dat',NActA,NBasis)
call read_1rdm_spin_dalton(GBAct,'rdms1_B.dat',NActB,NBasis)

! full spin matrices
do j=1,NActA
   do i=1,NActA
      GSpinA(INActA+i,INActA+j) = GAAct(i,j)
   enddo
enddo
do j=1,NActB
   do i=1,NActB
      GSpinB(INActB+i,INActB+j) = GBAct(i,j)
   enddo
enddo

! construct alpha/beta densities
allocate(SAPT%monA%g1a(NBasis,NBasis), &
         SAPT%monB%g1b(NBasis,NBasis))

!SAPT%monA%g1a = 0.5d0 * ( GChrgA + abs(GSpinA) )
!SAPT%monA%g1b = 0.5d0 * ( GChrgA - abs(GSpinA) )
SAPT%monA%g1a = 0.5d0 * ( GChrgA + GSpinA )
SAPT%monA%g1b = 0.5d0 * ( GChrgA - GSpinA )

  !print*, 'G1a = '
  !call print_sqmat(SAPT%monA%g1a,NBasis)
  !print*, 'G1b = '
  !call print_sqmat(SAPT%monA%g1b,NBasis)

SAPT%monB%g1a = 0.5d0 * ( GChrgB + GSpinB )
SAPT%monB%g1b = 0.5d0 * ( GChrgB - GSpinB )

deallocate(GBAct,GAAct)
deallocate(GSpinB,GSpinA)
deallocate(GChrgB,GChrgA)

end subroutine sapt_interface_spin

subroutine sapt_erfint_OTF(Flags,Mon,NBasis,AOBasis,CholErfVecsOTF)
!
! generate Long-Range Cholesky vectors in AO
!
! Comments:
!  - when called in CBS[H], sets Omega=1.0
!  - saves NCholErf to Mon
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: Mon
type(TCholeskyVecsOTF) :: CholErfVecsOTF
type(TAOBasis)         :: AOBasis
integer,intent(in)     :: NBasis

type(TSystem)          :: System

character(:),allocatable :: XYZPath
character(:),allocatable :: BasisSet, BasisSetPath

integer :: i
double precision :: Omega
logical :: doRSH
logical :: SortAngularMomenta

doRSH = .false.
if(Flags%IFunSR==1.or.Flags%IFunSR==2) doRSH = .true.
if(.not.doRSH .and. Flags%IDBBSC/=2) stop "doRSH=F and CBS/=2! in Erf ERIs OTF!"

XYZPath = "./input.inp"
SortAngularMomenta = .true.

! set basis set
BasisSet = Flags%BasisSetPath // Flags%BasisSet

! set RS parameter
if (Flags%IDBBSC==2) then
   Omega = 1.0
else
   Omega = Mon%Omega
endif

write(lout,'(/1x,3a6)') ('******',i=1,3)
write(lout,'(1x,a)') 'Cholesky LR On-The-Fly'
write(lout,'(1x,3a6)') ('******',i=1,3)

call auto2e_init()

!Mon%Omega = 100.0
!print*, 'Mon%OMega' , Mon%OMega
call CholeskyOTF_ao_vecs(CholErfVecsOTF,AOBasis,System,Flags%IUnits, &
                         XYZPath,BasisSet, &
                         SortAngularMomenta,Flags%ICholeskyAccu, &
                         Omega)

Mon%NCholErf = CholErfVecsOTF%Chol2Data%NVecs

end subroutine sapt_erfint_OTF

subroutine internal_orbgrid(Flags,AOBasis,System,Wg,Phi,NPoints,NAO)
implicit none

type(FlagsData) :: Flags
type(TAOBasis)  :: AOBasis
type(TSystem)   :: System
integer,intent(in)  :: NAO
integer,intent(out) :: NPoints

double precision, dimension(:, :), allocatable :: Phi
double precision, dimension(:), allocatable :: Wg

integer :: NAOt
double precision, dimension(:), allocatable :: Xg, Yg, Zg

character(:),allocatable :: XYZPath
character(:),allocatable :: BasisSet, BasisSetPath
logical :: SortAngularMomenta

logical, parameter :: SpherAO = .true.
integer, parameter :: GridType = BECKE_PARAMS_MEDIUM

BasisSet = Flags%BasisSetPath // Flags%BasisSet

! set gridtype : where??
! ...
! set units 
!Units = SYS_UNITS_BOHR

if(Flags%ICholeskyOTF/=1) then
   ! with Cholesky OTF AOBasis and System already avail

   XYZPath = "./input.inp"
   BasisSetPath = BasisSet
   SortAngularMomenta = .true.

   call auto2e_init()
   call sys_Read_XYZ(System, XYZPath)
   !call sys_Read_XYZ(System, XYZPath, Units)
   call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
   if (AOBasis%SpherAO) then
         NAOt = AOBasis%NAOSpher
   else
         NAOt = AOBasis%NAOCart
   end if
   if(NAOt /= NAO) then
     print*, 'NAO =',NAO, 'NAOlib',NAOt
     stop "sth wrong with NAO in internal_orbgrid!"
   endif
endif

! Molecular grid
call becke_MolecularGrid(Xg, Yg, Zg, Wg, NPoints, GridType, System, AOBasis)

! Atomic orbitals on the grid
allocate(Phi(NPoints, NAO))            
call gridfunc_Orbitals(Phi, Xg, Yg, Zg, NPoints, NAO, AOBasis)

end subroutine internal_orbgrid

subroutine internal_tran_orbgrid(OrbGrid,CAONO,Phi,AOBasis,ExternalOrdering,NPoints,NAO,NBasis)
!
! Phi = (NGrid,AO); CAONO(AO,NO)
! output: OrbGrid(NGrid,NO) = Phi.CAONO
!
implicit none

type(TAOBasis)     :: AOBasis
integer,intent(in) :: NPoints, NBasis, NAO
integer,intent(in) :: ExternalOrdering
double precision,intent(in)  :: Phi(NPoints,NAO), CAONO(NAO,NBasis)
double precision,intent(out) :: OrbGrid(NPoints,NBasis)

integer :: NAOt
double precision :: C_ao(NAO,NAO)

if(NBasis/=NAO) stop "NAO.ne.NBasis in internal_tran_orbgrid!"

! AOs from external program -> AOs in the Auto2e format
call auto2e_interface_C(C_ao, CAONO, AOBasis, ExternalOrdering)

print*, 'CAONO',norm2(CAONO)
print*, 'C_ao',norm2(C_ao)

call dgemm('N','N',NPoints,NBasis,NAO,1d0,Phi,NPoints,C_ao,NAO,0d0,OrbGrid,NPoints)
print*, 'OrbGrid',norm2(OrbGrid)

end subroutine internal_tran_orbgrid

subroutine arrange_mo(mat,nbas,SAPT)
implicit none

type(SaptData) :: SAPT
!integer :: NOrbA,NOrbB
integer :: nbas
double precision :: mat(nbas,nbas)

if(SAPT%monB%switchAB) then
   call gen_swap_rows(mat,nbas,nbas,SAPT%monA%NSym,&
                      SAPT%monA%NMonBas,SAPT%monB%NMonBas)
endif

!call swap_rows(NOrbA,NOrbB,mat)

end subroutine arrange_mo

!subroutine read_syminf(A,B,nbas)
!! reads number of basis functions on each monomer
!! from SYMINFO(B) file!
!implicit none
!
!type(SystemBlock) :: A, B
!integer :: nbas
!integer :: iunit,ios
!integer :: ibas,icen,last_ibas,last_icen
!integer :: irep,ifun,offset
!logical :: ex,dump
!integer :: tmp
!integer :: ACenTst, ACenBeg, ACenEnd
!
!! sanity checks : in
!!print*, A%NCen, B%NCen
!!print*, A%UCen, B%UCen
!if(A%NSym/=B%NSym) then
!  write(lout,*) 'ERROR in read_syminf: NSym different for A and B!'
!endif
!
!inquire(file='SYMINFO_B',EXIST=ex)
!
!if(ex) then
!   open(newunit=iunit,file='SYMINFO_B',status='OLD',&
!        form='FORMATTED')
!   read(iunit,*)
!   read(iunit,*)
!
!   ! old version: does not work with sym
!   ! print*, 'old version'
!   ! offset = 0
!   ! irep   = 1
!   ! read(iunit,'(i5,i6)',iostat=ios) last_ibas,last_icen
!   ! do
!   !   read(iunit,'(i5,i6)',iostat=ios) ibas,icen
!   !   if(ios/=0) then
!   !      A%NMonBas(irep)=last_ibas-offset
!   !      exit
!   !   elseif(icen/=last_icen) then
!   !        if(last_icen==B%UCen) then
!   !           B%NMonBas(irep) = last_ibas-offset
!   !           offset = last_ibas
!   !        elseif(icen==1) then
!   !           A%NMonBas(irep) = last_ibas-offset
!   !           offset = last_ibas
!   !           irep   = irep + 1
!   !        endif
!   !   endif
!   !   last_ibas=ibas
!   !   last_icen=icen
!   !enddo
!
!   ! new version : ok with sym
!   do irep=1,B%NSym
!      do ifun=1,B%NSymOrb(irep)
!         read(iunit,'(i5,i6)',iostat=ios) ibas,icen
!         if(icen.le.B%UCen) then
!            B%NMonBas(irep) = B%NMonBas(irep) + 1
!         else
!            A%NMonBas(irep) = A%NMonBas(irep) + 1
!         endif
!      enddo
!   enddo
!
!   close(iunit)
!else
!   write(LOUT,'(1x,a)') 'ERROR! MISSING SYMINFO_B FILE!'
!   stop
!endif
!
!! sanity checks : out
!do irep=1,B%NSym
!   ibas = A%NMonBas(irep)+B%NmonBas(irep)
!   if(ibas/=A%NSymOrb(irep)) then
!      write(lout,'(1x,a)') 'ERROR in read_syminf!'
!      write(lout,'(1x,a,i3,a)') 'For irep =',irep, ':'
!      write(lout,*) 'A-NMonBas',A%NMonBas(1:A%NSym)
!      write(lout,*) 'B-NMonBas',B%NMonBas(1:B%NSym)
!      write(lout,*) 'Sum:     ',A%NMonBas(1:A%NSym)+B%NMonBas(1:B%NSym)
!      write(lout,*) 'Should be',A%NSymOrb(1:A%NSym)
!      stop
!   endif
!enddo
!
!end subroutine read_syminf


subroutine sort_sym_mo(CMO,nbas,mon)
!
! requires: NSym, NSymOrb, INActS, NActS
!
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

! some prints
 print*, 'num0',mon%num0
 print*, 'num1',mon%num1
 print*, 'num2',mon%num2

! active pairs
 allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))

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

subroutine select_uactive(A,B,NBasis)
!
! set : occupied_sigma, virtual_sigma,
!       ov_sigma (=NDimX_sigma)
!       IndN_sigma
!       IGem
!
implicit none
type(SystemBlock)   :: A, B
integer,intent(in)  :: NBasis

integer :: ip,iq,ipq,ir,is,irs

! unrestricted: set active, inactive
A%NAct = 0
B%NAct = 0
A%NActOrb = 0
B%NActOrb = 0
A%INAct = int(sum(A%UOcc))
B%INAct = int(sum(B%UOcc))
A%SumOcc = sum(A%UOcc)
B%SumOcc = sum(B%UOcc)

! unresticed : set occ_sigma, virt_sigma
A%NOa = int(sum(A%UOcc(:,1)))
A%NOb = int(sum(A%UOcc(:,2)))
B%NOa = int(sum(B%UOcc(:,1)))
B%NOb = int(sum(B%UOcc(:,2)))

A%NVa = NBasis - A%NOa
A%NVb = NBasis - A%NOb
B%NVa = NBasis - B%NOa
B%NVb = NBasis - B%NOb

! unrestricted : set NDimX_sigma
A%NOVa = A%NOa*A%NVa
A%NOVb = A%NOb*A%NVb
B%NOVa = B%NOa*B%NVa
B%NOVb = B%NOb*B%NVb

! unrestricted : set IndN_sigma
! alpha
allocate(A%IndNa(2,A%NOVa),B%IndNa(2,B%NOVa))
ipq = 0
do iq=1,A%NOa
   do ip=1,A%NVa
      ipq = ipq + 1
      A%IndNa(1,ipq) = A%NOa + ip
      A%IndNa(2,ipq) = iq
   enddo
enddo
irs = 0
do is=1,B%NOa
   do ir=1,B%NVa
      irs = irs + 1
      B%IndNa(1,irs) = B%NOa + ir
      B%IndNa(2,irs) = is
   enddo
enddo
! beta
allocate(A%IndNb(2,A%NOVb),B%IndNb(2,B%NOVb))
ipq = 0
do iq=1,A%NOb
   do ip=1,A%NVb
      ipq = ipq + 1
      A%IndNb(1,ipq) = A%NOb + ip
      A%IndNb(2,ipq) = iq
   enddo
enddo
irs = 0
do is=1,B%NOb
   do ir=1,B%NVb
      irs = irs + 1
      B%IndNb(1,irs) = B%NOb + ir
      B%IndNb(2,irs) = is
   enddo
enddo

!print*, 'A: OCCUP_alpha = ', A%NOa
!print*, 'A: VIRT_alpha  = ', A%NVa
!print*, 'A: OCCUP_beta  = ', A%NOb
!print*, 'A: OCCUP_beta  = ', A%NVb

!print*, 'B: OCCUP_alpha = ', B%NOa
!print*, 'B: VIRT_alpha  = ', B%NVa
!print*, 'B: OCCUP_beta  = ', B%NOb
!print*, 'B: OCCUP_beta  = ', B%NVb

! unrestricted: set IGem
allocate(A%IGem(NBasis),B%IGem(NBasis))

A%NGem = 2
A%IGem(1:A%NAct+A%INAct) = 1
A%IGem(A%NAct+A%INAct+1:NBasis) = 2

B%NGem = 2
B%IGem(1:B%NAct+B%INAct) = 1
B%IGem(B%NAct+B%INAct+1:NBasis) = 2

end subroutine select_uactive

subroutine calc_elpot(A,B,CholeskyVecs,ICholesky,ICholeskyBIN,ICholeskyOTF,NBas)
implicit none

type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: ICholesky,NBas
integer,intent(in) :: ICholeskyBIN,ICholeskyOTF

integer :: ione,i,j
integer :: NInte1,NCholesky
double precision,allocatable :: Pa(:,:),Pb(:,:)
double precision,allocatable :: Va(:,:),Vb(:,:)
double precision,allocatable :: Ja(:,:),Jb(:,:)
logical                      :: valid
character(8)                 :: label

 NInte1 = NBas*(NBas+1)/2

 allocate(Pa(NBas,NBas),Pb(NBas,NBas),&
          Va(NBas,NBas),Vb(NBas,NBas),&
          Ja(NBas,NBas),Jb(NBas,NBas))

 !call get_den(NBas,A%CMO,A%Occ,2d0,Pa)
 !call get_den(NBas,B%CMO,B%Occ,2d0,Pb)
 Pa = 0d0
 do i=1,NBas
    call dger(NBas,NBas,2d0*A%Occ(i),A%CMO(:,i),1,A%CMO(:,i),1,Pa,NBas)
 enddo

 Pb = 0d0
 do i=1,NBas
    call dger(NBas,NBas,2d0*B%Occ(i),B%CMO(:,i),1,B%CMO(:,i),1,Pb,NBas)
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
    call make_J2(NBas,Pa,Pb,Ja,Jb)
 elseif(ICholeskyBIN==1) then
    NCholesky = CholeskyVecs%NCholesky
    call make_J2_CholR(CholeskyVecs%R(1:NCholesky,1:NInte1), &
                       Pa,Pb,Ja,Jb,NCholesky,NBas)
 elseif(ICholeskyOTF==1) then
   Ja = A%Jmat
   Jb = B%Jmat
   deallocate(A%Jmat,B%Jmat)
 endif

 allocate(A%WPot(NBas,NBas),B%WPot(NBas,NBas))

 A%WPot = Va + Ja
 B%WPot = Vb + Jb

 deallocate(Jb,Ja,Vb,Va,Pb,Pa)

end subroutine calc_elpot

subroutine calc_uks_elpot(M,CholeskyVecs,ICholesky,ICholeskyBIN,ICholeskyOTF,NBas)
!
! calculate open-shell W = V + Ja + Jb (in AO)
!
! CAREFUL, HERE: Pa(b) = P_alpha(beta)
!                Ja(b) = J_alpha(beta)
!
implicit none

type(SystemBlock)   :: M
type(TCholeskyVecs) :: CholeskyVecs

integer,intent(in) :: ICholesky,NBas
integer,intent(in) :: ICholeskyBIN,ICholeskyOTF

integer :: i,iunit

double precision,allocatable :: V(:,:)
double precision,allocatable :: Pa(:,:),Pb(:,:)
double precision,allocatable :: Ja(:,:),Jb(:,:)

logical :: valid
character(8)             :: label
character(:),allocatable :: onefile

if(ICholesky==1) stop "Cholesky not ready in calc_uks_elpot"

if (M%Monomer == 1 ) then
  onefile = "ONEEL_A"
elseif (M%Monomer == 2) then
  onefile = "ONEEL_B"
endif

allocate(Pa(NBas,NBas),Pb(NBas,NBas),&
         V(NBas,NBas),Ja(NBas,NBas),Jb(NBas,NBas))

! alpha dens
Pa = 0d0
do i=1,NBas
   call dger(NBas,NBas,M%UOcc(i,1),M%UMO(:,i,1),1,M%UMO(:,i,1),1,Pa,NBas)
enddo
! beta dens
Pb = 0d0
do i=1,NBas
   call dger(NBas,NBas,M%UOcc(i,2),M%UMO(:,i,2),1,M%UMO(:,i,2),1,Pb,NBas)
enddo

valid=.false.
V = 0d0
open(newunit=iunit,file=onefile,access='sequential',&
     form='unformatted',status='old')
read(iunit)
read(iunit) label,V
if(label=='POTENTAL') valid=.true.
close(iunit)
if(.not.valid) then
   write(LOUT,'(1x,a)') 'V not found in calc_elpot!'
endif

call make_J2(NBas,Pa,Pb,Ja,Jb)

allocate(M%WPot(NBas,NBas))

M%WPot = V + Ja + Jb

!allocate(M%Jos(NBas,NBas,2))

deallocate(Pb,Pa)
deallocate(Jb,Ja,V)

end subroutine calc_uks_elpot


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


!subroutine  square_oneint(tr,sq,nbas,nsym,norb)
!
!implicit none
!integer,intent(in) :: nbas,nsym,norb(8)
!double precision,intent(in) :: tr(:)
!double precision,intent(out) :: sq(nbas,nbas)
!integer :: irep,i,j
!integer :: offset,idx
!
!sq=0
!
!offset=0
!idx=0
!do irep=1,nsym
!   do j=offset+1,offset+norb(irep)
!      do i=offset+1,j
!
!         idx=idx+1
!         sq(i,j)=tr(idx)
!         sq(j,i)=tr(idx)
!
!      enddo
!   enddo
!   offset=offset+norb(irep)
!enddo
!
!end subroutine square_oneint


subroutine check_orbital_ordering(ICholeskyOTF)
!
! 1) Marcin's library uses ORBITAL_ORDERING param
! to distinguish between Orca, Dalton, ... interfaces
! 2) In GammCor ORBITAL ORDERING flas is set in fill_Flags()
! 3) this subroutine checks if the integers assigned to orderings 
!    are the same
!
implicit none

integer,intent(in) :: ICholeskyOTF

integer :: val

val = 0
if (ORBITAL_ORDERING_MOLPRO /= 1) val = 1
if (ORBITAL_ORDERING_ORCA   /= 2) val = 1
if (ORBITAL_ORDERING_DALTON /= 3) val = 1

if (val == 1) then
   write(lout,*) 'ORBITAL_ORDERING inconsistent between GammCor and gammcor-integrals!'
   if (ICholeskyOTF == 1) stop
endif

end subroutine check_orbital_ordering

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

subroutine print_uocc(nbas,SAPT)
!
implicit none

type(SaptData)     :: SAPT
integer,intent(in) :: nbas

integer :: i

 associate(A => SAPT%monA, B => SAPT%monB)
   write(LOUT,'(/1x,a,11x,a,5x,a)')   'UKS ORBITALS','Monomer A',  'Monomer B'
   write(LOUT,'(1x,a,11x,i3,11x,i3)') 'OCCUPIED ALPHA', A%NOa,   B%NOa
   write(LOUT,'(1x,a,11x,i3,11x,i3)') 'OCCUPIED BETA ', A%NOb,   B%NOb
   write(LOUT,'(1x,a,13x,i3,11x,i3)') 'OCCUPIED TOT',   A%INAct, B%INAct
   write(LOUT,'(/1x,a)') 'ORBITAL OCCUPANCIES'
   write(LOUT,'(1x,a,3x,a,4x,a,8x,a,6x,a)') 'UKS', 'Occ-A(alpha)', 'Occ-A(beta)', 'Occ-B(alpha)','Occ-B(beta)'
   do i=1,nbas
      write(6,'(x,i3,f10.6,f10.6,4x,f10.6,f10.6)') i,A%UOcc(i,1),A%UOcc(i,2),B%UOcc(i,1),B%UOcc(i,2)
   enddo
   write(LOUT,'(2x,a,f8.4,18x,f8.4/)') 'SUM OF OCCUPANCIES: ', A%SumOcc, B%SumOcc
 end associate

end subroutine print_uocc

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
