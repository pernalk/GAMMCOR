module sapt_main
use types
use systemdef
use timing
use abmat
use abfofo
use tran
use sorter
use sapt_ener

implicit none

!save
!double precision :: Tcpu,Twall

contains

subroutine sapt_driver(Flags,SAPT)
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT
integer :: i
integer :: NBasis
double precision :: Tcpu,Twall
! testowe

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
       call e1exchs2(SAPT%monA,SAPT%monB,SAPT)
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
!    call e2ind_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e2disp_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)

 endif

 call print_warn(SAPT)
 call free_sapt(SAPT)

 call clock('SAPT',Tcpu,Twall)

! testy molpro grid 
!  call molpro_routines

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
double precision,allocatable :: AuxA(:,:),AuxB(:,:)
double precision,allocatable :: OneRdmA(:),OneRdmB(:)
integer :: tmp1
double precision ::  potnuc,emy,eactiv,emcscf
integer :: noccA, nvirtA, noccB, nvirtB
integer :: ncen
logical :: doRSH

! read basis info
! only DCBS 
 NBasis = 0
 if(SAPT%InterfaceType==1) then
    call basinfo(NBasis,'SIRIUS_A.RST','DALTON')
 elseif(SAPT%InterfaceType==2) then
    call basinfo(NBasis,'AOTWOINT.mol','MOLPRO')
 endif
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

! temporary solution
 doRSH = .false.
 if(Flags%IFunSR==1.or.Flags%IFunSR==2) doRSH = .true.

! read and dump 1-electron integrals 
 if(SAPT%InterfaceType==1) then
    call onel_dalton(SAPT%monA%Monomer,NBasis,NSq,NInte1,SAPT%monA,SAPT)
    call onel_dalton(SAPT%monB%Monomer,NBasis,NSq,NInte1,SAPT%monB,SAPT)
 elseif(SAPT%InterfaceType==2) then
    call onel_molpro(SAPT%monA%Monomer,NBasis,NSq,NInte1,SAPT%monA,SAPT)
    call onel_molpro(SAPT%monB%Monomer,NBasis,NSq,NInte1,SAPT%monB,SAPT)
 endif

! read coefficient, occupancies
 if(SAPT%InterfaceType==1) then
    call readocc_dalton(NBasis,SAPT%monA,Flags) 
    call readocc_dalton(NBasis,SAPT%monB,Flags) 
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
 if(SAPT%monA%NSym.gt.1) then
    call sort_sym_mo(Ca,NBasis,SAPT%monA)
 endif
 if(SAPT%monB%NSym.gt.1) then
    call sort_sym_mo(Cb,NBasis,SAPT%monB)
 endif

! read 2-el integrals
 if(SAPT%InterfaceType==1) then
    call readtwoint(NBasis,1,'AOTWOINT_A','AOTWOSORT')
 elseif(SAPT%InterfaceType==2) then
    call readtwoint(NBasis,2,'AOTWOINT.mol','AOTWOSORT')
    if(doRSH) then
       if(SAPT%SameOm) then
          call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT')
       else 
          call readtwoint(NBasis,2,'AOTWOINT.erf','AOERFSORT')
          call readtwoint(NBasis,2,'AOTWOINT.erfB','AOERFSORTB')
       endif
    endif
 endif

 if(SAPT%InterfaceType==2) then
    call prepare_no(OneRdmA,AuxA,Ca,SAPT%monA,Flags%IFunSR,NBasis)
    call prepare_no(OneRdmB,AuxB,Cb,SAPT%monB,Flags%IFunSR,NBasis)
  
    call prepare_rdm2(SAPT%monA,AuxA,NBasis)
    call prepare_rdm2(SAPT%monB,AuxB,NBasis)
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
! transform 2-el integrals
! full 4-idx tran
! call tran4_full(NBasis,Ca,Cb,'TWOMOAB')

 if(Flags%ISERPA==0) then

   ! integrals stored as (ov|ov)
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
                  'TWOMOAB','AOTWOSORT')

   ! testing Ecorr ERPASYMM/ERPAVEC
   !call tran4_gen(NBasis,&
   !               SAPT%monA%num0+SAPT%monA%num1,Ca,&
   !               SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
   !               SAPT%monA%num0+SAPT%monA%num1,Ca,&
   !               SAPT%monA%num1+SAPT%monA%num2,Ca(NBasis*SAPT%monA%num0+1:NBasis**2),&
   !               'TMPMOAA','AOTWOSORT')
   !call tran4_gen(NBasis,&
   !               SAPT%monB%num0+SAPT%monB%num1,Cb,&
   !               SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
   !               SAPT%monB%num0+SAPT%monB%num1,Cb,&
   !               SAPT%monB%num1+SAPT%monB%num2,Cb(NBasis*SAPT%monB%num0+1:NBasis**2),&
   !               'TMPMOBB','AOTWOSORT')

   ! this is for testing E1exchS2...
   !call tran4_gen(NBasis,&
   !               NBasis,Ca,NBasis,Ca,&
   !               NBasis,Cb,NBasis,Cb,&
   !               'TMPMOAB','AOTWOSORT')

   ! <oo|oo>
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  'TMPOOAB','AOTWOSORT')

 elseif(Flags%ISERPA==2.and.Flags%ISHF==0) then
   ! integrals stored as (oFull,oFull)
   call tran4_gen(NBasis,&
                  SAPT%monA%num0+SAPT%monA%num1,Ca,&
                  NBasis,Ca,&
                  SAPT%monB%num0+SAPT%monB%num1,Cb,&
                  NBasis,Cb,&
                  'TWOMOAB','AOTWOSORT')
   ! this is for testing E1exchS2...
   call tran4_gen(NBasis,&
                  NBasis,Ca,NBasis,Ca,&
                  NBasis,Cb,NBasis,Cb,&
                  'TMPMOAB','AOTWOSORT')


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
    if(Flags%IFunSR/=0) then
       call calc_resp_dft(SAPT%monA,Ca,Flags,NBasis)
    else
       call calc_resp_full(SAPT%monA,Ca,Flags,NBasis,'TWOMOAA',SAPT%EnChck)
    endif
 endif

! mon B
    call SaptInter(NBasis,SAPT%monB,Flags%ICASSCF)
!   call COMMTST(NBasis) 
 if(SAPT%SaptLevel.eq.0) then
    call calc_resp_unc(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB')
 elseif(SAPT%SaptLevel.gt.0) then
    if(Flags%IFunSR/=0) then
       call calc_resp_dft(SAPT%monB,Cb,Flags,NBasis)
    else
       call calc_resp_full(SAPT%monB,Cb,Flags,NBasis,'TWOMOBB',SAPT%EnChck)
   endif
 endif

  if(Flags%ISERPA==2.and.Flags%ISHF==1) then
     call tran4_gen(NBasis,&
                      (SAPT%monA%num0+SAPT%monA%num1),SAPT%monA%CMO,&
                      NBasis,SAPT%monA%CMO,&
                      (SAPT%monB%num0+SAPT%monB%num1),SAPT%monB%CMO,&
                      NBasis,SAPT%monB%CMO,&
                      'TWOMOAB','AOTWOSORT')
 
     !call tran4_gen(NBasis,&
     !             NBasis,Ca,NBasis,Ca,&
     !             NBasis,Cb,NBasis,Cb,&
     !             'TMPMOAB','AOTWOSORT')

 endif


 ! calculate electrostatic potential
 call calc_elpot(SAPT%monA,SAPT%monB,NBasis)

! calc intermolecular repulsion
 SAPT%Vnn = calc_vnn(SAPT%monA,SAPT%monB)

 deallocate(Ca,Cb)
 if(SAPT%InterfaceType==2) then 
    deallocate(OneRdmB,OneRdmA,AuxB,AuxA)
 endif

end subroutine sapt_interface

subroutine onel_molpro(mon,NBasis,NSq,NInte1,MonBlock,SAPT)
 implicit none

 type(SaptData) :: SAPT
 type(SystemBlock) :: MonBlock

 integer,intent(in) :: mon,NBasis,NSq,NInte1
 integer :: ione,ios,NSym,NBas(8),ncen
 double precision, allocatable :: Hmat(:),Vmat(:),Smat(:) 
 double precision, allocatable :: work1(:),work2(:)
 character(8) :: label
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
 
 MonBlock%NSym = NSym
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

end subroutine onel_dalton 

subroutine readocc_dalton(NBasis,Mon,Flags) 
implicit none 

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

integer,intent(in) :: NBasis
integer :: NSym,NOrbt,NBasist,NCMOt,NOcc(8),NOrbs(8)
integer :: isiri
double precision :: potnuc,emy,eactiv,emcscf
character(:),allocatable :: occfile,sirifile,siriusfile,coefile 
logical :: exsiri


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

 else
    NBasist = NBasis
 endif 
    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
 
 if(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE/=1.and.(.not.Mon%ISHF)) then
    ! CASSCF
    call readmulti(NBasis,Mon,.false.,exsiri,isiri,occfile,siriusfile)

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE==1.and.Flags%ISERPA==0.and.(.not.Mon%ISHF)) then
    ! CASSCF
    ! for 2-el electron case: read from occupations.dat
    call readmulti(NBasis,Mon,.false.,.false.,isiri,occfile,siriusfile)

 elseif(Flags%ICASSCF==1.and.Flags%ISHF==0.and.Mon%NELE==1.and.Flags%ISERPA==2) then

    call readmulti(NBasis,Mon,.false.,exsiri,isiri,occfile,siriusfile)

 elseif(Flags%ICASSCF==1.and.(Flags%ISHF==1.or.Mon%ISHF)) then

    ! Hartree-Fock
    call readmulti(NBasis,Mon,.true.,exsiri,isiri,occfile,siriusfile)
    call readener(NBasis,Mon,isiri)

 elseif(Flags%IGVB==1) then 

    ! GVB
    call readgvb(Mon,NBasis,coefile)

 endif

 if(exsiri) close(isiri) 

end subroutine readocc_dalton

subroutine readocc_molpro(NBasis,Mon,OrbAux,OneRdm,Flags)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

integer :: NBasis
integer :: NInte1,HlpDim,NOccup
integer :: i,info
double precision :: Tmp
double precision :: OrbAux(NBasis,NBasis), &
                    OneRdm(NBasis*(NBasis+1)/2)
double precision,allocatable :: EVal(:)
double precision,allocatable :: work(:)
character(:),allocatable :: rdmfile

 NInte1 = NBasis*(NBasis+1)/2
 HlpDim = max(NBasis**2,3*NBasis)

 if(Mon%Monomer==1) then
   rdmfile = '2RDMA'
 elseif(Mon%Monomer==2) then
   rdmfile = '2RDMB'
 endif 

 allocate(Mon%CICoef(NBasis),Mon%IGem(NBasis),Mon%Occ(NBasis))
 allocate(work(HlpDim),EVal(NBasis))
 OneRdm = 0
 ! HERE!!! FIRST STATE FOR NOW
 call read_1rdm_molpro(OneRdm,Mon%InSt(1,1),Mon%InSt(2,1),&
                       rdmfile,Mon%IWarn,NBasis)
 call triang_to_sq2(OneRdm,OrbAux,NBasis)
 call Diag8(OrbAux,NBasis,NBasis,Eval,work)
! call dsyev('V','U',NBasis,OrbAux,NBasis,EVal,work,3*NBasis,info)
! print*,EVal
 call SortOcc(EVal,OrbAux,NBasis)

 Mon%NAct = 0
 Tmp = 0
 do i=1,NBasis
    Tmp = Tmp + EVal(i)
    if(EVal(i)>0.d0) Mon%NAct = Mon%NAct + 1  
 enddo
 !Mon%INAct = Mon%NELE-int(Tmp)
!!! OPEN-SHELL CASE?
 Mon%INAct = Mon%NELE-Tmp+1.d-1
 NOccup = Mon%INAct + Mon%NAct
 Mon%SumOcc = Tmp + Mon%INAct
!
 !print*, Mon%NAct,Mon%INAct
 !print*, Mon%NELE,Tmp

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

subroutine prepare_no(OneRdm,OrbAux,OrbCAS,Mon,IFunSR,NBasis)
implicit none
! OrbCAS[inout] :: on input AOtoCAS
!                  on output AOtoNO
! OrbAux        :: on input CAStoNO

type(SystemBlock) :: Mon

integer :: IFunSR,NBasis
double precision :: OneRdm(NBasis*(NBasis+1)/2)
double precision :: OrbAux(NBasis,NBasis),OrbCAS(NBasis,NBasis)
double precision,allocatable :: URe(:,:),OrbSym(:,:),Fock(:)
double precision,allocatable :: work1(:),work2(:),work3(:)
integer :: NOccup,NVirt,NSym
integer :: i,j,ia,ib,iab,ioff,idx,NInte1
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

   call FockGen_mithap(work2,OneRdm,work1,NInte1,NBasis,'AOTWOSORT')

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
    print*, 'INACTIVE:',work3(1:Mon%INAct)
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
! END OF CANONICALIZING

 call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,OrbSym,NBasis,0d0,OrbCAS,NBasis)
 OrbCAS = transpose(OrbCAS)

 deallocate(work3,work2,work1,Fock,OrbSym,URe)
 deallocate(Mon%IndInt,Mon%NumOSym)

end subroutine prepare_no

subroutine prepare_rdm2(Mon,OrbAux,NBasis)
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

end subroutine prepare_rdm2

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
     EigY(:), EigY1(:), Eig(:), Eig1(:), &
     testY(:),testEig(:)

integer :: i,ione
double precision,parameter :: One = 1d0, Half = 0.5d0
character(8) :: label
character(:),allocatable :: onefile,twofile,propfile0,propfile1,rdmfile

! perform check
! if(Flags%ICASSCF==0) then
!   write(LOUT,*) 'ERROR! E2DISP UNC AVAILABLE ONLY FOR CAS/HF!'
!   stop
! endif

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
 call tran4_full(NBas,MO,MO,fname,'AOTWOSORT')
 call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 ! GVB
 if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

    allocate(EigY(Mon%NDimX**2), &
         Eig(Mon%NDimX),Eig1(Mon%NDimX))
    
    EigY  = 0
    Eig   = 0
    Eig1  = 0

    call Y01GVB(TwoMO,Mon%Occ,URe,XOne, &
         EigY,Eig,Eig1, &
         Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2)
    
    ! dump uncoupled response
    call writeresp(EigY,Eig,propfile0)
    if(Flags%IFlag0==0) then
       call writeEval(Eig1,propfile1)
    endif

    deallocate(Eig1,Eig,EigY)
    
    ! CAS
 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then

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
 
    !call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
    !     EigY,EigY1,Eig,Eig1, &
    !     Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
    !     NBas,Mon%NDim,NInte1,twofile,Flags%IFlag0)
    
    call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
         EigY,EigY1,Eig,Eig1, &
         Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)

    ! dump response
    call writeresp(EigY,Eig,propfile0)
    if(Flags%IFlag0==0) then
       call writeresp(EigY1,Eig1,propfile1)
    endif
    
    deallocate(Eig1,Eig,EigY1,EigY,ABMin,ABPlus)
   
 endif

 close(ione)
 deallocate(TwoMO,URe,XOne,work1,work2)

end subroutine calc_resp_unc

subroutine calc_resp_dft(Mon,MO,Flags,NBas)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

double precision :: MO(:)        
integer :: NBas
integer :: NSq,NInte1,NInte2,NGrid,NDimKer
integer :: NSymNO(NBas),MultpC(15,15)
double precision, allocatable :: work1(:),work2(:)
double precision, allocatable :: XOne(:), &
                                 TwoMO(:), TwoElErf(:), &
                                 WGrid(:),XKer(:),OrbGrid(:),&
                                 OrbXGrid(:),OrbYGrid(:),OrbZGrid(:),&
                                 SRKer(:),SRKerW(:)
double precision, allocatable :: ABPlus(:),ABMin(:),URe(:,:),VSR(:), &
                                 EigY0(:),EigY1(:),Eig0(:),Eig1(:), &
                                 EigVecR(:), Eig(:) 
integer :: i,j,ii,ione
double precision :: ACAlpha,Omega,EnSR,EnHSR,ECorr,ECASSCF,XVSR
double precision :: Tcpu,Twall
character(8) :: label
character(:),allocatable :: onefile,aoerfile,twofile,twoerffile,&
                            twojfile,twokfile,twojerf,twokerf, & 
                            propfile,propfile0,propfile1,rdmfile
double precision,parameter :: One = 1d0, Half = 0.5d0
logical :: doRSH

! temporary RSH solution
 doRSH = .false.
 if(Flags%IFunSR==1.or.Flags%IFunSR==2) doRSH = .true.

! set filenames
 if(Mon%Monomer==1) then
    onefile = 'ONEEL_A'
    twofile = 'TWOMOAA'
    twojfile = 'FFOOAA'
    twokfile = 'FOFOAA'
    twoerffile = 'MO2ERFAA'
    twojerf = 'FFOOERFAA'
    twokerf = 'FOFOERFAA'
    propfile = 'PROP_A'
    propfile0 = 'PROP_A0'
    propfile1 = 'PROP_A1'
    rdmfile='rdm2_A.dat'
    aoerfile = 'AOERFSORT'
 elseif(Mon%Monomer==2) then
    onefile = 'ONEEL_B'
    twofile = 'TWOMOBB'
    twojfile = 'FFOOBB'
    twokfile = 'FOFOBB'
    twoerffile = 'MO2ERFBB'
    twojerf = 'FFOOERFBB'
    twokerf = 'FOFOERFBB'
    propfile = 'PROP_B'
    propfile0 = 'PROP_B0'
    propfile1 = 'PROP_B1'
    rdmfile='rdm2_B.dat'
    if(Mon%SameOm) then
       aoerfile = 'AOERFSORT'
    else
       aoerfile = 'AOERFSORTB'
    endif
 endif

! set dimensions
 NSq = NBas**2
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

 call create_symmats(Mon,MO,NBas)

 allocate(work1(NSq),work2(NSq),XOne(NInte1),URe(NBas,NBas),VSR(NInte1))
 if(Mon%TwoMoInt==1) then
    allocate(TwoMO(NInte2))
    allocate(TwoElErf(NInte2))
    !if(doRSH) allocate(TwoElErf(NInte2))
 endif
 
 URe = 0d0
 do i=1,NBas
    URe(i,i) = 1d0
 enddo

 ! read 1-el
 ! add Coulomb (RSH: sr Coulomb)
 call triang_to_sq(Mon%VCoul,work2,NBas)
 open(newunit=ione,file=onefile,access='sequential',&
      form='unformatted',status='old')

 read(ione) 
 read(ione)
 read(ione) label, work1
 if(label=='ONEHAMIL') then
    work1 = work1 + work2
    call tran_oneint(work1,MO,MO,work2,NBas)
    call sq_to_triang(work1,XOne,NBas) 
 else
    write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
    stop
 endif

 ! load grid
 call molprogrid0(NGrid,NBas)
 write(LOUT,'()')
 write(LOUT,'(1x,a,i8)') "The number of Grid Points =",NGrid

 allocate(WGrid(NGrid),SRKer(NGrid),SRKerW(NGrid), &
          OrbGrid(NGrid*NBas), &
          OrbXGrid(NGrid*NBas),OrbYGrid(NGrid*NBas),OrbZGrid(NGrid*NBas))

 ! load orbgrid and gradients, and wgrid
 call transp_mat1dim(MO,work1,NBas) 
 call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid, &
                                WGrid,work1,NGrid,NBas)

 ! set/load Omega - range separation parameter
 write(LOUT,'(1x,a,f15.8)') "The range-separation parameter =",Mon%Omega
 
 ! transform 2-el integrals
 select case(Mon%TwoMoInt)
 case(TWOMO_INCORE,TWOMO_FFFF) 
    ! full 
    call tran4_full(NBas,MO,MO,twofile,'AOTWOSORT')
 case(TWOMO_FOFO) 
    ! transform J and K
    call tran4_gen(NBas,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             NBas,MO,&
             NBas,MO,&
             twojfile,'AOTWOSORT')
    call tran4_gen(NBas,&
             NBas,MO,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             NBas,MO,&
             Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
             twokfile,'AOTWOSORT')
 end select
 if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 ! trasform LR integrals
 if(doRSH) then
   select case(Mon%TwoMoInt)
   case(TWOMO_INCORE)
      TwoElErf(1:NInte2) = 0
     print*, 'Check: calc_resp_dft:',Mon%SameOm,aoerfile
      call tran4_full(NBas,MO,MO,twoerffile,aoerfile)
   case(TWOMO_FFFF) 
      call tran4_full(NBas,MO,MO,twoerffile,aoerfile)
   case(TWOMO_FOFO) 
      call tran4_gen(NBas,&
              Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
              Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
              NBas,MO,&
              NBas,MO,&
              twojerf,aoerfile)
      call tran4_gen(NBas,&
              NBas,MO,&
              Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
              NBas,MO,&
              Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
              twokerf,aoerfile)
   end select
   if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoEl(Mon%Monomer+4,TwoElErf,NBas,NInte2)
 endif

 !print*, 'TwoMO ',norm2(TwoMO)
 !print*, 'TwoErf',norm2(TwoElerf)

 NSymNO(1:NBas) = 1
 call EPotSR(EnSR,EnHSR,VSR,Mon%Occ,URe,MO,.true.,& 
            OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,WGrid,&
!            NSymNO,TwoMO,TwoElErf,&
            NSymNO,Mon%VCoul,&
            Mon%Omega,Flags%IFunSR,&
            NGrid,NInte1,NInte2,NBas)

 ! MODIFY XOne in PostCAS calculations
 if(Mon%PostCAS) then
    write(LOUT,*) 'TEST POSTCAS!'
    XOne = XOne + VSR 
 endif
 write(LOUT,'(1x,a,f15.8)') "SR Energy: ",EnSR

 XVSR = 0 
 do i=1,NBas
    ii = (i*(i+1))/2
    XVSR = XVSR + 2d0*Mon%Occ(i)*VSR(ii)
 enddo

 ACAlpha=One

 allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2),&
          EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

 ! read 2-RDMs
 call read2rdm(Mon,NBas)
 call system('cp '//rdmfile// ' rdm2.dat')

 ECASSCF = 0
! if(doRSH) then 
 select case(Mon%TwoMoInt)
 case(TWOMO_INCORE)
    call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoElErf,Mon%IPair,&
                Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)
    call EKT(URe,Mon%Occ,XOne,TwoElErf,NBas,NInte1,NInte2)
 case(TWOMO_FFFF)
    call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                NInte1,twoerffile,ACAlpha,.false.)
 case(TWOMO_FOFO)
    call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                NInte1,twojerf,twokerf,ACAlpha,.false.)
!else
! HERE:: ADD SEPARATE PROCEDURE FOR Kohn-Sham! 
!endif
 end select
 write(LOUT,'(/,1x,a,f16.8,a,1x,f16.8)') 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)
 
 ! ADD CONTRIBUTIONS FROM THE (sr)ALDA KERNEL TO AB MATRICES
 MultpC(1,1)=1
 call GetKerNPT(SRKer,Mon%Occ,URe,OrbGrid,WGrid,NSymNO,MultpC, &
                NBas,NGrid)
 call clock('START',Tcpu,Twall)
 select case(Mon%TwoMoInt)
 case(TWOMO_INCORE)
    call ModABMin(Mon%Occ,SRKer,WGrid,OrbGrid,TwoMO,TwoElErf,ABMin,&
                  Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NInte2,NBas)
    !call ModABMin_old(Mon%Occ,SRKer,WGrid,OrbGrid,TwoMO,TwoElErf,ABMin,&
    !                Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NInte2,NBas)
    print*, 'ABMin-Kasia',norm2(ABMin)
 case(TWOMO_FFFF)
    call ModABMin_mithap(Mon%Occ,SRKer,WGrid,OrbGrid,ABMin,&
                         Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NBas,&
                         twofile,twoerffile)
 case(TWOMO_FOFO)
    call ModABMin_FOFO(Mon%Occ,SRKer,WGrid,OrbGrid,ABMin,&
                       Mon%MultpC,Mon%NSymNO,&
                       Mon%IndN,Mon%IndX,Mon%NDimX,NGrid,NBas,&
                       Mon%num0,Mon%num1, & 
                       twokfile,twokerf,.false.)
    print*, 'ABMin-MY',norm2(ABMin)
 end select
 call clock('Mod ABMin',Tcpu,Twall)

 !write(LOUT,'(1x,a)') "*** sr-kernel added. ***"
 ! test true energy
 write(LOUT,'(/,1x,a,f15.8)') "Total lrCASSCF+ENuc+srDF Energy", ECASSCF-XVSR+EnSR+Mon%PotNuc

 EigVecR = 0
 Eig = 0
 if(Mon%NoSt==1) then
    call ERPASYMM1(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
 else
    call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
 endif

 write(LOUT,'(/," *** LR-CAS-SR-DFT Excitation Energies *** ",/)')
 do i=1,10
    write(LOUT,'(i4,4x,e16.6)') i,Eig(i)
 enddo 

 !print*, 'Check! Eig,EigVecR',norm2(Eig),norm2(EigVecR)

 !ECorr=0
 !!call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
 !!     Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
 !  call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Mon%Occ, &
 !                      Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
 !                      Mon%NDimX,NBas,twokfile)
 !ECorr=Ecorr*0.5d0
 !write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
 !     ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr

! dump response
 call writeresp(EigVecR,Eig,propfile)

 ! uncoupled
 allocate(EigY0(Mon%NDimX**2),EigY1(Mon%NDimX**2),&
          Eig0(Mon%NDimX),Eig1(Mon%NDimX))
 
 EigY0 = 0
 EigY1 = 0
 Eig0 = 0
 Eig1 = 0

 !call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
 !       EigY0,EigY1,Eig0,Eig1, &
 !       Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
 !       NBas,Mon%NDim,NInte1,twoerffile,Flags%IFlag0)
 
! call Y01CAS(TwoElErf,Mon%Occ,URe,XOne,ABPlus,ABMin, &
!      EigY0,EigY1,Eig0,Eig1, &
!      Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)
! do i=1,NGrid
!    SRKerW(i) = SRKer(i)*WGrid(i)
! enddo
! call Y01CASLR(TwoElErf,Mon%Occ,URe,XOne,ABPlus,ABMin, &
!        EigY0,EigY1,Eig0,Eig1, &
!        Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0, &
!        TwoMO,OrbGrid,SRKerW,NSymNO,MultpC,NGrid)

 ! dump uncoupled response
 call writeresp(EigY0,Eig0,propfile0)
 if(Flags%IFlag0==0) then
    call writeresp(EigY1,Eig1,propfile1)
 endif

 close(ione)
 deallocate(Eig1,Eig0,EigY1,EigY0)
 
 deallocate(Eig,EigVecR,ABMin,ABPlus)
 deallocate(SRKer,SRKerW,OrbZGrid,OrbYGrid,OrbXGrid,OrbGrid,WGrid)
 if(Mon%TwoMoInt==1) then
    deallocate(TwoMO)
    ! if(doRSH) deallocate(TwoElErf)
    deallocate(TwoElErf)
 endif
 deallocate(VSR,URe,XOne,work1,work2)

end subroutine calc_resp_dft

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
                                 EigVecR(:), Eig(:), &
                                 ABPlusT(:), ABMinT(:)
double precision, allocatable :: Eig0(:), Eig1(:), EigY0(:), EigY1(:) 
double precision :: Dens(NBas,NBas)
double precision,parameter :: One = 1d0, Half = 0.5d0
integer :: i,j,ij,ij1,ione,itwo 
integer :: ind
double precision :: ACAlpha
double precision :: ECASSCF,ETot,ECorr
character(8) :: label
character(:),allocatable :: onefile,twofile,propfile,rdmfile
character(:),allocatable :: twojfile,twokfile
character(:),allocatable :: propfile0,propfile1
double precision :: tmp
double precision,parameter :: SmallE=0d0,BigE=1.D20
double precision,external  :: trace
integer :: INegExcit,offset
integer,allocatable :: sort(:)
double precision, allocatable :: EigTmp(:), VecTmp(:)
! testy
 integer :: ii,jj,dimV1
 double precision,parameter :: thresh=1.d-10
 integer :: space1(2,Mon%NDimX)
 integer :: IndN_tmp(2,Mon%NDim)
 double precision :: work(Mon%NDimX)
 double precision,external :: ddot
 double precision :: tst
 integer :: DimEx

! set filenames
 if(Mon%Monomer==1) then
    onefile  = 'ONEEL_A'
    twofile  = 'TWOMOAA'
    twojfile = 'FFOOAA'
    twokfile = 'FOFOAA'
    propfile = 'PROP_A'
    propfile0  = 'PROP_A0' 
    propfile1  = 'PROP_A1' 
    rdmfile='rdm2_A.dat'
 elseif(Mon%Monomer==2) then
    onefile  = 'ONEEL_B'
    twofile  = 'TWOMOBB'
    twojfile = 'FFOOBB'
    twokfile = 'FOFOBB'
    propfile = 'PROP_B'
    propfile0  = 'PROP_B0' 
    propfile1  = 'PROP_B1' 
    rdmfile='rdm2_B.dat'
 endif

! set dimensions
 NSq = NBas**2
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

 allocate(work1(NSq),work2(NSq),XOne(NInte1),URe(NBas,NBas))
 if(Mon%TwoMoInt==TWOMO_INCORE) then
    allocate(TwoMO(NInte2))
 endif 

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

 ! transform 2-el integrals
 select case(Mon%TwoMoInt)

 case(TWOMO_INCORE,TWOMO_FFFF) 
   ! full - for GVB and CAS
   call tran4_full(NBas,MO,MO,fname,'AOTWOSORT')

 case(TWOMO_FOFO) 
   ! transform J and K
    call tran4_gen(NBas,&
         Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
         Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
         NBas,MO,&
         NBas,MO,&
         twojfile,'AOTWOSORT')
    call tran4_gen(NBas,&
         NBas,MO,&
         Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
         NBas,MO,&
         Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
         twokfile,'AOTWOSORT')
 end select
 if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoEl(Mon%Monomer,TwoMO,NBas,NInte2)

 if(Flags%ISHF==1.and.Flags%ISERPA==2.and.Mon%NELE==1) then

    ! Calculate FCI for 2-el systems
    call calc_fci(NBas,NInte1,NInte2,XOne,URe,TwoMO,Mon)

 endif

 ACAlpha=One
 ! GVB
 if(Flags%ICASSCF==0.and.Flags%ISERPA==0) then

  ! allocate(ABPlus(Mon%NDim**2),ABMin(Mon%NDim**2), &
   allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2), &
            EigVecR(Mon%NDimX**2),Eig(Mon%NDimX))

! HERE: STARTED WORK ON ACABMAT0_mithap!
  ! ACAlpha=sqrt(2d0)/2.d0
   call ACABMAT0_mithap(ABPlus,ABMin,URe,Mon%Occ,XOne,&
                 Mon%IndN,Mon%IndX,Mon%IGem,Mon%CICoef,&
                 Mon%NAct,Mon%INAct,NBas,Mon%NDim,Mon%NDimX,NInte1,Mon%NGem,&
                 twofile,Flags%ISAPT,ACAlpha,1)
!   allocate(ABPlusT(Mon%NDim**2),ABMinT(Mon%NDim**2))
!   call ACABMAT0(ABPlusT,ABMinT,URe,Mon%Occ,XOne,TwoMO, &
!   call ACABMAT0(ABPlus,ABMin,URe,Mon%Occ,XOne,TwoMO, &
!                 NBas,Mon%NDim,NInte1,NInte2,Mon%NGem,ACAlpha,1)

   ! reduce dim
   !call reduce_dim('AB',ABPlus,ABMin,Mon)

   ! TESTY
   !call reduce_dim('AB',ABPlusT,ABMinT,Mon)
   ! write(*,*) 'TEST_PLUS',norm2(ABPlus),norm2(ABPlusT(1:Mon%NDimX**2))
   ! write(*,*) 'TEST_MIN',norm2(ABMin),norm2(ABMinT(1:Mon%NDimX**2))
   ! write(*,*) 'ABPLUS: ',norm2(ABPlus(1:Mon%NDimX**2)-ABPlusT(1:Mon%NDimX**2))
   ! write(*,*) 'ABMIN: ',norm2(ABMin(1:Mon%NDimX**2)-ABMinT(1:Mon%NDimX**2))

   ! do j=1,Mon%NDimX**2
   !   write(*,*) j,ABPlus(j),ABPlusT(j)
   !enddo
   
   !deallocate(ABMinT,ABPlusT)

   EigVecR = 0
   Eig = 0
   call ERPASYMM(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)

   if(EChck) then
      write(LOUT,'(/,1x,a)') 'ERPA-GVB ENERGY CHECK REQUESTED:'
      call EneERPA(ETot,ECorr,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,&
           Mon%Occ,XOne,Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
   endif

   !! uncoupled
   !allocate(EigY0(Mon%NDimX**2),Eig0(Mon%NDimX),Eig1(Mon%NDimX))
 
   !EigY0 = 0
   !Eig0 = 0
   !Eig1 = 0

   !call Y01GVB(TwoMO,Mon%Occ,URe,XOne, &
   !     EigY0,Eig0,Eig1, &
   !     Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2)

   !! dump uncoupled response
   !call writeresp(EigY0,Eig0,propfile0)
   !if(Flags%IFlag0==0) then
   !   call writeEval(Eig1,propfile1)
   !endif

   !deallocate(Eig1,Eig0,EigY0)

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

   !ACAlpha=sqrt(2d0)/2d0
   !print*, 'ACAlpha',ACAlpha
   select case(Mon%TwoMoInt)
   case(TWOMO_FOFO)

      call AB_CAS_FOFO(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                  Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                  NInte1,twojfile,twokfile,ACAlpha,.false.)
   case(TWOMO_FFFF) 

      call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
                  Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
                   NInte1,twofile,ACAlpha,.false.)
   case(TWOMO_INCORE)

      call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
                  Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)

      call EKT(URe,Mon%Occ,XOne,TwoMO,NBas,NInte1,NInte2)

   end select 
   write(LOUT,'(/,1x,a,f16.8,a,1x,f16.8)') 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)

   EigVecR = 0
   Eig = 0
   if(Mon%NoSt==1) then
      call ERPASYMM1(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
   else
      call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
   endif

   !print*, 'Entering ERPAVEC...' 
   !call ERPAVEC(EigVecR,Eig,ABPlus,ABMin,NBas,Mon%NDimX)
   !

  !print*, 'EigVecR',norm2(EigVecR)
  !do i=1,size(Eig)
  !   print*, i,Eig(i)
  !enddo

   if(EChck) then
      ECorr=0
      select case(Mon%TwoMoInt) 
      case(TWOMO_FOFO)
         call ACEneERPA_FOFO(ECorr,EigVecR,Eig,Mon%Occ, &
                              Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
                              Mon%NDimX,NBas,twokfile)
      case(TWOMO_FFFF)
         call ACEneERPA_FFFF(ECorr,EigVecR,Eig,Mon%Occ, &
                              Mon%IGem,Mon%IndN,Mon%IndX,Mon%num0+Mon%num1, &
                              Mon%NDimX,NBas,twofile)
      case(TWOMO_INCORE)
         call ACEneERPA(ECorr,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne,&
                        Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NGem)
      end select
      ECorr=Ecorr*0.5d0

      write(LOUT,'(/,1x,''ECASSCF+ENuc, Corr, ERPA-CASSCF'',6x,3f15.8)') &
           ECASSCF+Mon%PotNuc,ECorr,ECASSCF+Mon%PotNuc+ECorr
   endif

   ! uncoupled
   ! call test_resp_unc(Mon,URe,XOne,TwoMO,NBas,NInte1,NInte2,Flags%IFlag0) 

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

   select case(Mon%TwoMoInt)
   case(TWOMO_FOFO) 
      call Y01CAS_FOFO(Mon%Occ,URe,XOne,ABPlus,ABMin, &
             !EigY0,EigY1,Eig0,Eig1, &
             propfile0,propfile1, &
             Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
             NBas,Mon%NDim,NInte1,Mon%NoSt,twofile,twojfile,twokfile,Flags%IFlag0)
   case(TWOMO_FFFF) 
      call Y01CAS_mithap(Mon%Occ,URe,XOne,ABPlus,ABMin, &
             EigY0,EigY1,Eig0,Eig1, &
             Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX, &
             NBas,Mon%NDim,NInte1,twofile,Flags%IFlag0)
   case(TWOMO_INCORE) 
      call Y01CAS(TwoMO,Mon%Occ,URe,XOne,ABPlus,ABMin, &
           EigY0,EigY1,Eig0,Eig1, &
           !Mon%IndNT,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)
           Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDim,NInte1,NInte2,Flags%IFlag0)
   end select

   ! dump uncoupled response
   if(Mon%TwoMoInt/=TWOMO_FOFO) then 
      call writeresp(EigY0,Eig0,propfile0)
      if(Flags%IFlag0==0) then
         call writeresp(EigY1,Eig1,propfile1)
      endif
   endif
  
   deallocate(Eig1,Eig0,EigY1,EigY0)

elseif(Flags%ISERPA==2) then
   ! TD-APSG response: GVB, CASSCF, FCI

   !if(Flags%ICASSCF==1.and.MON%NELE/=1) then
   !   write(LOUT,'(1x,a)') 'ERROR!!! FCI POSSIBLE ONLY FOR 2-EL MONOMERS!'
   !   stop 
   !endif

   if(Flags%ICASSCF==1.and.Flags%ISHF==0) then
      ! read 2-RDMs
      call read2rdm(Mon,NBas)
      ! CAS-SCF
      !call execute_command_line('cp '//rdmfile// ' rdm2.dat')
      call system('cp '//rdmfile// ' rdm2.dat')
   elseif(Flags%ICASSCF==1.and.Flags%ISHF==1) then
      ! PINO
      call read2rdm(Mon,NBas)
      call init_pino(NBas,Mon,Flags%ICASSCF) 
   endif

   Mon%NDimN=0
   do i=1,NBas
      if(Mon%Occ(i).gt.0d0) Mon%NDimN=Mon%NDimN+1
   enddo
   ! print*, 'NDimN: ',Mon%NDimN
   if(Flags%ICASSCF==1.and.Flags%ISHF==0) then

    !  allocate(ABPlus(Mon%NDimX**2),ABMin(Mon%NDimX**2), &
    !       EMAT(NBas**2),EMATM(NBas**2),&
    !       DMAT(Mon%NDim*NBas),DMATK(Mon%NDim*NBas),&
    !       !EigVecR(2*(Mon%NDimX+Mon%NDimN)*2*(Mon%NDimX+Mon%NDimN)),&
    !       EigVecR((Mon%NDimX+Mon%NDimN)*(Mon%NDimX+Mon%NDimN)),&
    !       !Eig(2*(Mon%NDimX+Mon%NDimN)))
    !       Eig((Mon%NDimX+Mon%NDimN)))

    !  ACAlpha = One
    !  call AB_CAS(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
    !       Mon%IndN,Mon%IndX,Mon%NDimX,NBas,Mon%NDimX,NInte1,NInte2,ACAlpha)
    !  !call AB_CAS_mithap(ABPlus,ABMin,ECASSCF,URe,Mon%Occ,XOne, &
    !  !          Mon%IndN,Mon%IndX,Mon%IGem,Mon%NAct,Mon%INAct,Mon%NDimX,NBas,Mon%NDimX,&
    !  !          NInte1,twofile,ACAlpha,.false.)
    !  write(*,*) 'DE_CAS!'
    !  !Mon%NDimN = 0
    !  call DE_CAS(DMAT,DMATK,EMAT,EMATM,URe,Mon%Occ,XOne,TwoMO,Mon%IPair,&
    !       NBas,Mon%NDim,NInte1,NInte2,ACAlpha)

    !  print*, 'DMAT',norm2(DMAT),norm2(DMATK)
    !  print*, 'EMAT',norm2(EMAT),norm2(EMATM)
    !  ij=0
    !  ij1=0
    !  do j = 1,Mon%NDimN
    !     do i = 1,Mon%NDimX
    !        ij  = ij + 1 !(j-1)*Mon%NDimX + i
    !        ij1 = (j-1)*Mon%NDim + Mon%IndXh(i)
    !        DMAT(ij) = DMAT(ij1)
    !        DMATK(ij) = DMATK(ij1)
    !     enddo
    !  enddo
    !  ij=0
    !  ij1=0 
    !  do j=1,Mon%NDimN
    !     do i=1,Mon%NDimN
    !        ij  = (j-1)*Mon%NDimN + i
    !        ij1 = (j-1)*Mon%NBasis + i
    !        EMAT(ij) = EMAT(ij1)
    !        EMATM(ij) = EMATM(ij1)
    !     enddo
    !  enddo
      !print*, 'DMAT-R',norm2(DMAT(1:Mon%NDimX*Mon%NDimN)),&
      !     norm2(DMATK(1:Mon%NDimX*Mon%NDimN))
      !print*, 'EMAT-R',norm2(EMAT(1:Mon%NDimN**2)),&
      !     norm2(EMATM(1:Mon%NDimN**2))
      !print*, 'ABPlus',norm2(ABPlus),'ABMin',norm2(ABMin)
     ! code for H2,He/Be only!
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
      print*, 'DMAT-R',norm2(DMAT(1:Mon%NDimX*Mon%NDimN)),&
           norm2(DMATK(1:Mon%NDimX*Mon%NDimN))
      print*, 'EMAT-R',norm2(EMAT(1:Mon%NDimN**2)),&
           norm2(EMATM(1:Mon%NDimN**2))
      print*, 'ABPlus',norm2(ABPlus(1:Mon%NDimX**2)),&
           'ABMin',norm2(ABMin(1:Mon%NDimX**2))

   else
  
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

      !print*, 'DMAT-R',norm2(DMAT(1:Mon%NDimX*Mon%NDimN)),&
      !         norm2(DMATK(1:Mon%NDimX*Mon%NDimN))
      !print*, 'EMAT-R',norm2(EMAT(1:Mon%NDimN**2)),&
      !         norm2(EMATM(1:Mon%NDimN**2))
      !print*, 'ABPlus',norm2(ABPlus(1:Mon%NDimX**2)),&
      !         norm2(ABMin(1:Mon%NDimX**2))

      deallocate(CMAT)

   endif
   !    write(*,*) Mon%NGem,Mon%NDimN,Mon%NDimX
   !    write(*,*) Mon%Occ,NBas

   EigVecR = 0
   Eig = 0
  ! call PINOVEC(EigVecR,Eig,ABPlus,ABMin,DMAT,DMATK,EMAT,EMATM, &
  !      Mon%Occ,NBas,Mon%NDimX,Mon%NDimN)
  !

! reduced version
   call PINOVECRED(EigVecR,Eig,INegExcit,ABPlus,ABMin,DMAT,DMATK, &
        EMAT,EMATM,NBas,Mon%NDimX,Mon%NDimN)
!! the lines below are needed if EnePINO is called subsequently
!   EigVecR(1:2*(NDim+NBasis)*2*(NDim+NBasis))=Zero 
!   do i=1,NDimX+NDimN
!     Eig(i+NDimX+NDimN)=Zero
!      do j=1,NDimX+NDimN
!      If(j.Le.NDimX) EigVecR((i-1)*2*(NDimX+NDimN)+j)=
!     $ EigVecRR((i-1)*(NDimX+NDimN)+J)
!      If(j.Gt.NDimX) EigVecR((i-1)*2*(NDimX+NDimN)+J+NDimX)=
!     $ EigVecRR((I-1)*(NDimX+NDimN)+J)
!      enddo
!   enddo 
!
   ETot = 0
!   call EnePINO(ETot,Mon%PotNuc,EigVecR,Eig,TwoMO,URe,Mon%Occ,XOne, &
!        Mon%IndN,NBas,NInte1,NInte2,Mon%NDimX,Mon%NDimN)
   
   deallocate(EMAT,EMATM,DMAT,DMATK)
   
endif

! dump response
 call writeresp(EigVecR,Eig,propfile)

 close(ione)
 
 deallocate(work1,work2,XOne,URe)
 if(Mon%TwoMoInt==1) deallocate(TwoMO)
 deallocate(ABPlus,ABMin,EigVecR,Eig)

 ! delete TWOMO file
 !call delfile(twofile)

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
         ij  = ij + 1 !(j-1)*Mon%NDimX + i
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

subroutine read_mo_dalton(cmo,nbasis,nsym,nbas,norb,nsiri,nmopun)

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

end subroutine read_mo_dalton

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

subroutine writeEval(EVal,mon)
implicit none

character(*) :: mon
double precision :: EVal(:)
integer :: iunit

 open(newunit=iunit,file=mon,form='unformatted')
 write(iunit) EVal
 close(iunit)

end subroutine writeEval

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
 !print*, 'TTESSSTT:', mon%num0,mon%num1,mon%num2, nbas

!! test!
! do i=1,mon%NELE
!    mon%IndAux(i)=0
! enddo
! do i=1+mon%NELE,nbas
!    mon%IndAux(i)=2
! enddo
!
!mon%num0 = 2 
!mon%num1 = 1
!mon%num2 = nbas - mon%num0 - mon%num2

! active pairs
 allocate(mon%IPair(nbas,nbas),mon%IndX(mon%NDim),mon%IndN(2,mon%NDim))

 mon%IPair(1:nbas,1:nbas) = 0

 if(Flags%ICASSCF==0) then
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
                  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1.d-2)) then
                  
                write(LOUT,'(1x,a,2x,2i4)') 'Discarding nearly degenerate pair',i,j 
             else
                ! if IFlCore=0 exclude core (inactive) orbitals
                if(Flags%IFlCore==1.or.&
                     (Flags%IFlCore==0.and.&
                     mon%Occ(i)/=1d0.and.mon%Occ(j)/=1d0) ) then
                   
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
   
       enddo
    enddo

 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==0) then
    write(LOUT,'(1x,a,e15.5)') 'Threshold for active orbitals: ', mon%ThrSelAct
    if(mon%NCen==1.and.mon%ThrSelAct<1.d-3) then 
       write(LOUT,'(1x,a)') 'Warning! For single atom ThrSelAct should probably have larger value!'
       mon%IWarn = mon%IWarn + 1
    endif 
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
                  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.mon%ThrSelAct) ) then
                ! here!!!
                !  (Abs(mon%Occ(i)-mon%Occ(j))/mon%Occ(i).lt.1d-3) ) then
                
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
    
 elseif(Flags%ICASSCF==1.and.Flags%ISERPA==2.and.mon%NELE==1) then

    allocate(mon%IndXh(mon%NDim))
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
                mon%IndX(ind) =  ind !ij 
                mon%IndXh(ind) = ij
                mon%IndN(1,ind) = i
                mon%IndN(2,ind) = j
                mon%IPair(i,j) = 1
                mon%IPair(j,i) = 1
                
             endif
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
   write(LOUT,'(2x,a,f8.4,18x,f8.4)') 'SUM OF OCCUPANCIES: ', A%SumOcc, B%SumOcc
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

!RSH
if(allocated(SAPT%monA%VCoul)) then
  deallocate(SAPT%monA%VCoul)
endif
if(allocated(SAPT%monB%VCoul)) then
  deallocate(SAPT%monB%VCoul)
endif

! delete AOTWOSORT
 call delfile('AOTWOSORT')

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

