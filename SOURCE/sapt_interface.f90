module sapt_inter

use types
use timing
use tran
use sorter
!use Cholesky_old

use gammcor_integrals
!use Auto2eInterface
!use Cholesky, only : chol_CoulombMatrix, TCholeskyVecs, &
!                     chol_Rkab_ExternalBinary, chol_MOTransf_TwoStep
!use basis_sets
!use sys_definitions
!use CholeskyOTF_interface
!use CholeskyOTF, only : TCholeskyVecsOTF
!use Cholesky_driver, only : chol_Rkab_OTF
!
!use BeckeGrid
!use GridFunctions
!use grid_definitions

use abmat
use read_external

implicit none
!private

!public TCholeskyVecs,TCholeskyVecsOTF

contains

subroutine sapt_interface(Flags,SAPT,NBasis,AOBasis,CholeskyVecsOTF)
!
! SAPT-DALTON requires SIRIFC and SIRIUS.RST
!                   or SIRIFC and occupations.dat
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
write(LOUT,'(/1x,"Flags:BasisSetPath ",a)') Flags%BasisSetPath
write(LOUT,'(1x, "Flags:BasisSet ",a)')     Flags%BasisSet
BasisSet = Flags%BasisSetPath // Flags%BasisSet

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
    call onel_molpro(SAPT%monA%Monomer,NBasis,NSq,NInte1,SAPT%monA,SAPT)
    call onel_molpro(SAPT%monB%Monomer,NBasis,NSq,NInte1,SAPT%monB,SAPT)
 endif

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

    if(SAPT%InterfaceType==1) then ! Dalton

       call readtwoint(NBasis,1,'AOTWOINT_A','AOTWOSORT',MemSrtSize)
       if(doRSH) then

         if (SAPT%SameOm) then
            call readtwoint(NBasis,1,'AOERFINT_A','AOERFSORT',MemSrtSize)
         else
            call readtwoint(NBasis,1,'AOERFINT_A','AOERFSORT',MemSrtSize)
            call readtwoint(NBasis,1,'AOERFINT_B','AOERFSORTB',MemSrtSize)
         endif
         !call readtwoint(NBasis,1,'AOSR2INT','AOSR2SORT',MemSrtSize)  ! skip that for now

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
       endif
    endif
 endif

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

       block
        character(:),allocatable :: XYZPath
        character(:),allocatable :: BasisSetPath

        XYZPath = "./input.inp"
        BasisSetPath = BasisSet
        SortAngularMomenta = .true.

        print*, 'WARNING: move CholeskyOTF_ao_vecs to new module!'
        call CholeskyOTF_ao_vecs(CholeskyVecsOTF,AOBasis,System, &
                                 XYZPath,BasisSetPath, &
                                 SortAngularMomenta,Flags%ICholeskyAccu)

       SAPT%NCholesky  = CholeskyVecsOTF%NVecs
       SAPT%monA%NChol = SAPT%NCholesky
       SAPT%monB%NChol = SAPT%NCholesky

       end block

    endif

 endif
 call clock('2ints',Tcpu,Twall)

 if(SAPT%InterfaceType==2) then

    if(SAPT%monA%NatOrb==0) then
       ! create NOs inside GammCor (use canonical CAS orbs)
       call prepare_no_molpro(Ca,OneRdmA,AuxA,SAPT%monA,AOBasis,System, &
                       CholeskyVecs,CholeskyVecsOTF, &
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
                       CholeskyVecs,CholeskyVecsOTF, &
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
 if(SAPT%SaptExch==1.and.Flags%IRDM2Typ==11) call prepare_rdm2_approx(SAPT%monA,Flags%IRDM2Typ,NBasis)
 if(SAPT%SaptExch==1.and.Flags%IRDM2Typ==11) call prepare_rdm2_approx(SAPT%monB,Flags%IRDM2Typ,NBasis)

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

! orbitals on a grid
block
!integer :: InternalGrid, NGrid
integer :: NGrid
double precision, allocatable :: Phi(:,:),WGrid(:)

SAPT%InternalGrid = 1
if (SAPT%InternalGrid==1) then
   if (Flags%IFunSR/=0) then
      print*,'Internal grid...',Flags%IFunSR
      print*, 'ORBITAL_ORDERING =',Flags%ORBITAL_ORDERING
      call internal_orbgrid(Flags,AOBasis,System,WGrid,Phi,NGrid,NBasis)

      SAPT%NGrid      = NGrid
      SAPT%monA%NGrid = NGrid
      SAPT%monB%NGrid = NGrid
      
      allocate(SAPT%monA%WGrid(NGrid)) ! remove that
      allocate(SAPT%monB%WGrid(NGrid)) ! remove that
      SAPT%monA%WGrid = WGrid
      SAPT%monB%WGrid = WGrid

      allocate(SAPT%monA%CNOGrid(NGrid,NBasis))
      call internal_tran_orbgrid(SAPT%monA%CNOGrid,SAPT%monA%CMO,Phi,AOBasis, &
                                 Flags%ORBITAL_ORDERING,NGrid,NBasis,NBasis)


      allocate(SAPT%monB%CNOGrid(NGrid,NBasis))
      call internal_tran_orbgrid(SAPT%monB%CNOGrid,SAPT%monB%CMO,Phi,AOBasis, &
                                 Flags%ORBITAL_ORDERING,NGrid,NBasis,NBasis)

   endif
endif
end block

! look-up tables
 call select_active(SAPT%monA,NBasis,Flags)
 call select_active(SAPT%monB,NBasis,Flags)

 ! transform Cholesky Vecs to NO
 if(Flags%ICholeskyBIN==1) then
    !call chol_sapt_AO2NO_BIN(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call chol_OO_sapt_AO2NO_BIN(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
    call chol_FO_sapt_AO2NO_BIN(SAPT,SAPT%monA,SAPT%monB,CholeskyVecs,NBasis,Flags%MemVal,Flags%MemType)
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
    call make_K(NBasis,work,SAPT%monB%Kmat)
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
    call read_orbinf_dalton(sirifile,NSym,Mon%NOrb,Mon%NSymOrb)
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

 print*, 'read2rdm: MON%RDM2',norm2(Mon%RDM2)

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
                      CholeskyVecs,CholeskyVecsOTF, &
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
double precision,intent(in)    :: CMONOAct(NBasis,NBasis)
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
print*, 'RDM2val =',norm2(RDM2val)

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

subroutine CholeskyOTF_elpot_AO(Mon,NBasis)
!
! obtain electrostatic potential in AO
! W = V + J
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis

integer      :: ione
logical      :: valid
character(8) :: label
character(:),allocatable :: onefile

double precision :: V(NBasis,NBasis)

if(Mon%Monomer==1) then
  onefile="ONEEL_A"
elseif(Mon%Monomer==2) then
  onefile="ONEEL_B"
endif

V = 0d0
open(newunit=ione,file=onefile,access='sequential',&
     form='unformatted',status='old')
read(ione)
read(ione) label,V
if(label=='POTENTAL') valid=.true.
close(ione)

if(.not.valid) then
   write(LOUT,'(1x,a)') 'V not found in calc_elpot_CholOTF!'
endif

allocate(Mon%WPot(NBasis,NBasis))

 Mon%WPot = V + Mon%Jmat

end subroutine CholeskyOTF_elpot_AO

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

subroutine chol_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
! transform Cholesky Vecs from AO to NO
! for all 2-index vecs needed in SAPT:
!   FFXX(NCholesky,NBas**2) -- for polarization
!   FFXY(NCholesky,NBas**2) -- for exchange
!   FFYX(NCholesky,NBas**2) -- for exchange
!   OOXX(NCholesky,dimO**2) -- for polarization
!   DCholX(NCholesky,NDimX) -- for polarization
!
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB,dimVA,dimVB,nOVA,nOVB
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

print*, 'A%OO',norm2(A%OO)
print*, 'B%OO',norm2(B%OO)

! if(SAPT%SaptLevel==666) then ! RS2PT2+
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
! endif
print*, 'A%OOAB',norm2(A%OOAB)
print*, 'B%OOBA',norm2(B%OOBA)

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

!print*, 'A%FF',norm2(A%FF)
!print*, A%FF(3,:)
!print*, 'B%FF',norm2(B%FF)

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

 print*, 'A%FFAB',norm2(A%FFAB)
 print*, 'B%FFBA',norm2(B%FFBA)

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

end subroutine chol_sapt_AO2NO_BIN

subroutine chol_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

if(SAPT%InterfaceType==1) then
  write(lout,*) 'Cholesky 3-index AO2NO transformation does not work with DALTON yet!'
  stop
endif

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

allocate(A%OO(NCholesky,dimOA**2),B%OO(NCholesky,dimOB**2))

call chol_Rkab_OTF(A%OO,A%CAONO,1,dimOA,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('AOO',Tcpu,Twall)

call chol_Rkab_OTF(B%OO,B%CAONO,1,dimOB,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('BOO',Tcpu,Twall)

print*, 'A%OO',norm2(A%OO)
print*, 'B%OO',norm2(B%OO)

!if(SAPT%SaptLevel==666) then ! RS2PT2+

   allocate(A%OOAB(NCholesky,dimOA*dimOB), &
            B%OOBA(NCholesky,dimOB*dimOA))

   call chol_Rkab_OTF(A%OOAB,A%CAONO,1,dimOA,B%CAONO,1,dimOB, &
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING_MOLPRO)

   call chol_Rkab_OTF(A%OOBA,B%CAONO,1,dimOB,A%CAONO,1,dimOA, &
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING_MOLPRO)

print*, 'A%OOAB',norm2(A%OOAB)
print*, 'B%OOBA',norm2(B%OOBA)
!endif

allocate(A%FF(NCholesky,NBasis**2),&
         B%FF(NCholesky,NBasis**2) )

call chol_Rkab_OTF(A%FF,A%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('AFF',Tcpu,Twall)

call chol_Rkab_OTF(B%FF,B%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call clock('BFF',Tcpu,Twall)

print*, 'A%FF',norm2(A%FF)
print*, 'B%FF',norm2(B%FF)

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

!call chol_Rkab_OTF(MatFF, UAux, a0, a1, UAux, b0, b1,
!                   MaxBufferDimMB, CholeskyVecsOTF,
!                   AOBasis, ORBITAL_ORDERING_MOLPRO)

allocate(A%FFAB(NCholesky,NBasis**2),&
         B%FFBA(NCholesky,NBasis**2) )

call chol_Rkab_OTF(A%FFAB,A%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

call chol_Rkab_OTF(B%FFBA,B%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING_MOLPRO)

print*, 'A%FFAB',norm2(A%FFAB)
print*, 'B%FFBA',norm2(B%FFBA)

end subroutine chol_sapt_AO2NO_OTF

subroutine chol_OO_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB

!   OOXX(NCholesky,dimO**2) -- for polarization

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

allocate(A%OO(NCholesky,dimOA**2),&
         B%OO(NCholesky,dimOB**2) )

call chol_MOTransf_TwoStep(A%OO,CholeskyVecs,&
                    A%CMO,1,dimOA,&
                    A%CMO,1,dimOA,&
                    MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%OO,CholeskyVecs,&
                   B%CMO,1,dimOB,&
                   B%CMO,1,dimOB,&
                   MaxBufferDimMB)

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

end subroutine chol_OO_sapt_AO2NO_BIN

subroutine chol_OO_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set orbital ordering
if(SAPT%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(SAPT%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

allocate(A%OO(NCholesky,dimOA**2),B%OO(NCholesky,dimOB**2))

call chol_Rkab_OTF(A%OO,A%CAONO,1,dimOA,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING)

call clock('AOO',Tcpu,Twall)

call chol_Rkab_OTF(B%OO,B%CAONO,1,dimOB,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,       &
                   AOBasis,ORBITAL_ORDERING)

call clock('BOO',Tcpu,Twall)

allocate(A%OOAB(NCholesky,dimOA*dimOB), &
            B%OOBA(NCholesky,dimOB*dimOA))

call chol_Rkab_OTF(A%OOAB,A%CAONO,1,dimOA,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

call clock('AOOAB',Tcpu,Twall)

call chol_Rkab_OTF(B%OOBA,B%CAONO,1,dimOB,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

print*, 'A%OOAB',norm2(A%OOAB)
print*, 'A%OOBA',norm2(B%OOBA)

end subroutine chol_OO_sapt_AO2NO_OTF

subroutine chol_FO_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
! prepare (NChol,NBasis*dimO) matrices for FOFO integrals (BIN version)
! (FO|AA), (FO|BB), (FO|AB), (FO|BA)
!
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB

! set dimensions
NCholesky = CholeskyVecs%NCholesky
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

allocate(A%FO(NCholesky,NBasis*dimOA),B%FO(NCholesky,NBasis*dimOB))

call chol_MOTransf_TwoStep(A%FO,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   A%CMO,1,dimOA, &
                   MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%FO,CholeskyVecs,&
                   B%CMO,1,NBasis,&
                   B%CMO,1,dimOB, &
                   MaxBufferDimMB)

allocate(A%FOAB(NCholesky,NBasis*dimOB),B%FOBA(NCholesky,NBasis*dimOA))

call chol_MOTransf_TwoStep(A%FOAB,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   B%CMO,1,dimOB, &
                   MaxBufferDimMB)
call chol_MOTransf_TwoStep(B%FOBA,CholeskyVecs,&
                   B%CMO,1,NBasis,&
                   A%CMO,1,dimOA, &
                   MaxBufferDimMB)

end subroutine chol_FO_sapt_AO2NO_BIN

subroutine chol_FF_sapt_AO2NO_BIN(SAPT,A,B,CholeskyVecs,NBasis,MemVal,MemType)
!
!   FFXX(NCholesky,NBas**2) -- for polarization
!   FFXY(NCholesky,NBas**2) -- for exchange
!
implicit none

type(SaptData)      :: SAPT
type(SystemBlock)   :: A, B
type(TCholeskyVecs) :: CholeskyVecs
integer,intent(in)  :: NBasis
integer,intent(in)  :: MemVal,MemType

integer          :: NCholesky
integer          :: MaxBufferDimMB
integer          :: dimOA,dimOB
integer          :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecs%NCholesky

! set buffer size
if(MemType == 2) then       !MB
   MaxBufferDimMB = MemVal
elseif(MemType == 3) then   !GB
   MaxBufferDimMB = MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx Cholesky transformation'

allocate(A%FF(NCholesky,NBasis**2),&
        B%FF(NCholesky,NBasis**2) )
call chol_MOTransf_TwoStep(A%FF,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   A%CMO,1,NBasis,&
                   MaxBufferDimMB)
call clock('AFF',Tcpu,Twall)

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

allocate(A%FFAB(NCholesky,NBasis**2))
!         B%FFBA(NCholesky,NBasis**2) )
 
call chol_MOTransf_TwoStep(A%FFAB,CholeskyVecs,&
                   A%CMO,1,NBasis,&
                   B%CMO,1,NBasis,&
                   MaxBufferDimMB)
!call chol_MOTransf_TwoStep(B%FFBA,CholeskyVecs,&
!                   B%CMO,1,NBasis,&
!                   A%CMO,1,NBasis,&
!                   MaxBufferDimMB)

print*, 'A%FFAB',norm2(A%FFAB)

end subroutine chol_FF_sapt_AO2NO_BIN

subroutine chol_FO_sapt_AO2NO_OTF(SAPT,A,B,CholeskyVecsOTF,AOBasis,Flags,NBasis)
!
! prepare (NChol,NBasis*dimO) matrices for FOFO integrals (OTF version)
! (FO|AA), (FO|BB), (FO|AB), (FO|BA)
!
implicit none

type(SaptData)         :: SAPT
type(SystemBlock)      :: A, B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
type(FlagsData)        :: Flags
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: dimOA,dimOB
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs
dimOA = A%num0+A%num1
dimOB = B%num0+B%num1

! set orbital ordering
if(SAPT%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(SAPT%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif

allocate(A%FO(NCholesky,NBasis*dimOA),B%FO(NCholesky,NBasis*dimOB))

call chol_Rkab_OTF(A%FO,A%CAONO,1,NBasis,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,        &
                   AOBasis,ORBITAL_ORDERING)

call chol_Rkab_OTF(B%FO,B%CAONO,1,NBasis,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,        &
                   AOBasis,ORBITAL_ORDERING)

call clock('AFO+BFO',Tcpu,Twall)

allocate(A%FOAB(NCholesky,NBasis*dimOB),B%FOBA(NCholesky,NBasis*dimOA))

call chol_Rkab_OTF(A%FOAB,A%CAONO,1,NBasis,B%CAONO,1,dimOB, &
                   MaxBufferDimMB,CholeskyVecsOTF,          &
                   AOBasis,ORBITAL_ORDERING)

call chol_Rkab_OTF(B%FOBA,B%CAONO,1,NBasis,A%CAONO,1,dimOA, &
                   MaxBufferDimMB,CholeskyVecsOTF,          &
                   AOBasis,ORBITAL_ORDERING)

call clock('AFOAB+BFOBA',Tcpu,Twall)

end subroutine chol_FO_sapt_AO2NO_OTF

subroutine chol_FFXX_mon_AO2NO_OTF(Flags,M,CholeskyVecsOTF,AOBasis,NBasis)
!
! performs AO2NO transformation
! to generate 1) FFXX(NCholesky,NBasis**2) vectors
!             2) (c_p+c_q)*FFXX => DChol vectors
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: M
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
integer,intent(in)     :: NBasis

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING
integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

! set orbital ordering
if(Flags%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(Flags%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

allocate(M%FF(NCholesky,NBasis**2))

call chol_Rkab_OTF(M%FF,M%CAONO,1,NBasis,M%CAONO,1,NBasis,&
                   MaxBufferDimMB,CholeskyVecsOTF, &
                   AOBasis,ORBITAL_ORDERING)

if (M%Monomer==1) call clock('AFF',Tcpu,Twall)
if (M%Monomer==2) call clock('BFF',Tcpu,Twall)

allocate(M%DChol(NCholesky,M%NDimX))

do j=1,M%NDimX
   ip = M%IndN(1,j)
   iq = M%IndN(2,j)
   ipq = iq + (ip-1)*NBasis
   Cpq = M%CICoef(ip) + M%CICoef(iq)
   M%DChol(:,j) = Cpq*M%FF(:,ipq)
enddo

end subroutine chol_FFXX_mon_AO2NO_OTF

subroutine chol_FFXY_AB_AO2NO_OTF(Flags,A,B,CholeskyVecsOTF,AOBasis,NBasis,abtype)
!
! performs AO2NO transformation
! to generate 1) FFXY(NCholesky,NBasis**2) vectors
!
implicit none

type(FlagsData)        :: Flags
type(SystemBlock)      :: A,B
type(TAOBasis)         :: AOBasis
type(TCholeskyVecsOTF) :: CholeskyVecsOTF
integer,intent(in)     :: NBasis
character(2),intent(in)   :: abtype

integer :: NCholesky
integer :: MaxBufferDimMB
integer :: ORBITAL_ORDERING

integer :: i,j,ip,iq,ipq
double precision :: Cpq

double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
NCholesky = CholeskyVecsOTF%NVecs

! set buffer size
if(Flags%MemType == 2) then       !MB
   MaxBufferDimMB = Flags%MemVal
elseif(Flags%MemType == 3) then   !GB
   MaxBufferDimMB = Flags%MemVal * 1024_8
endif
!write(lout,'(1x,a,i5,a)') 'Using ',MaxBufferDimMB,' MB for 3-indx AO2NO transformation'

! set orbital ordering
if(Flags%InterfaceType==1) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_DALTON
elseif(Flags%InterFaceType==2) then
   ORBITAL_ORDERING = ORBITAL_ORDERING_MOLPRO
else
   print*, 'SAPT with Cholesky OTF does not work with this Interface!'
   stop
endif

if (abtype == "AB") then

   allocate(A%FFAB(NCholesky,NBasis**2))
   call chol_Rkab_OTF(A%FFAB,A%CAONO,1,NBasis,B%CAONO,1,NBasis,&
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING)

elseif (abtype == "BA") then

   allocate(B%FFBA(NCholesky,NBasis**2))
   call chol_Rkab_OTF(B%FFBA,B%CAONO,1,NBasis,A%CAONO,1,NBasis,&
                      MaxBufferDimMB,CholeskyVecsOTF, &
                      AOBasis,ORBITAL_ORDERING)

endif

end subroutine chol_FFXY_AB_AO2NO_OTF

subroutine chol_JKmat_AO_OTF(Mon,NBasis)
!
! generates J and K matrices in MO
! and backtransforms them to AO
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis

integer :: i,j
integer :: ione
integer :: NOccup,NCholesky
double precision :: val
double precision :: D_no(NBasis,NBasis)
double precision :: SC(NBasis,NBasis),SAO(NBasis,NBasis)
double precision,allocatable :: ints(:),work(:,:)
double precision,allocatable :: Jtmp(:,:),Ktmp(:,:)
character(8) :: label
double precision,external :: ddot

NOccup = Mon%num0+Mon%num1
NCholesky = Mon%NChol

! prepare NO 1-density
D_no = 0d0
do i=1,NOccup
   D_no(i,i) = Mon%Occ(i)
enddo

! prepare Jmat in AO
allocate(Mon%Jmat(NBasis,NBasis))
if(.not.allocated(Mon%Kmat)) allocate(Mon%Kmat(NBasis,NBasis))

allocate(Jtmp(NBasis,NBasis),Ktmp(NBasis,NBasis))
allocate(ints(NBasis**2),work(NBasis,NBasis))

ints = 0d0
Jtmp = 0d0
Ktmp = 0d0
do i=1,NCholesky
   ints(:) = Mon%FF(i,:)
   val = ddot(NBasis**2,ints,1,D_no,1)
   call daxpy(NBasis**2,2d0*val,ints,1,Jtmp,1)
   call dgemm('N','N',NBasis,NBasis,NBasis,1d0,ints,NBasis, &
              D_no,NBasis,0d0,work,NBasis)
   call dgemm('N','N',NBasis,NBasis,NBasis,-1d0,work,NBasis, &
              ints,NBasis,1d0,Ktmp,NBasis)
enddo

!print*, 'Jtmp', norm2(Jtmp)

! backtransform J to AO
! get S in AO
open(newunit=ione,file='ONEEL_A',access='sequential',&
     form='unformatted',status='old')
read(ione) label, SAO
close(ione)

! J_AO = SC . J_NO . (SC)^T
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SAO,NBasis,Mon%CAONO,NBasis,0d0,SC,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,Jtmp,NBasis,0d0,work,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work,NBasis,SC,NBasis,0d0,mon%Jmat,NBasis)
!print*, 'Jmat-AO', norm2(Mon%Jmat)

! K_AO = SC . K_NO . (SC)^T
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,SC,NBasis,Ktmp,NBasis,0d0,work,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,-1d0,work,NBasis,SC,NBasis,0d0,mon%Kmat,NBasis)

deallocate(Ktmp,Jtmp)
deallocate(work,ints)

end subroutine chol_JKmat_AO_OTF 

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
