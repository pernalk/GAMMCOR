module types
! written by M. Hapka, M. Modrzejewski
use iso_fortran_env

implicit none


integer :: LOUT = output_unit
integer,parameter :: LERR = error_unit

integer, parameter :: INTER_TYPE_DAL = 1
integer, parameter :: INTER_TYPE_MOL = 2
integer, parameter :: INTER_TYPE_OWN = 3

integer, parameter :: TYPE_NO_SYM = 1
integer, parameter :: TYPE_SYM = 0

integer, parameter :: JOB_TYPE_AC    = 1
integer, parameter :: JOB_TYPE_AC0   = 2
integer, parameter :: JOB_TYPE_AC1   = 3
integer, parameter :: JOB_TYPE_ERPA  = 3
integer, parameter :: JOB_TYPE_EERPA = 4
integer, parameter :: JOB_TYPE_SAPT  = 5
integer, parameter :: JOB_TYPE_PDFT  = 6
integer, parameter :: JOB_TYPE_CASPIDFT = 7

integer, parameter :: SAPTLEVEL0 = 0
integer, parameter :: SAPTLEVEL1 = 1
integer, parameter :: SAPTLEVEL2 = 2

integer, parameter :: FLAG_FRAG = 0
integer, parameter :: FLAG_CORE = 1
integer, parameter :: FLAG_NOBASIS = 0
logical, parameter :: FLAG_RESTART = .FALSE.
integer, parameter :: FLAG_PRINT_LEVEL = 0
integer, parameter :: FLAG_DEBUG_FL12 = 1

integer, parameter :: RDM_TYPE_GVB  = 1
integer, parameter :: RDM_TYPE_APSG = 2
integer, parameter :: RDM_TYPE_CAS  = 3
integer, parameter :: RDM_TYPE_DMRG = 4
integer, parameter :: RDM_TYPE_HF   = 5

integer, parameter :: TWOMO_INCORE = 1
integer, parameter :: TWOMO_FFFF   = 2
integer, parameter :: TWOMO_FOFO   = 3

integer, parameter :: DF_NONE       = 0
integer, parameter :: DF_SRLDA      = 1
integer, parameter :: DF_SRPBE      = 2
integer, parameter :: DF_PBE        = 3

integer, parameter :: MONOMER_A = 1
integer, parameter :: MONOMER_B = 2

integer, parameter :: RESP_ERPA = 1
integer, parameter :: RESP_APSG = 2
integer, parameter :: RESP_DFT  = 3

logical, parameter :: FLAG_POSTCAS = .FALSE.

integer,parameter :: maxcen = 500

character(*),parameter :: PossibleInterface(3) = &
[character(8) :: &
'DALTON', 'MOLPRO', 'OWN']

character(*),parameter :: PossibleJobType(7) = &
[character(8) :: &
'AC', 'AC0', 'ERPA', 'AC1', 'SAPT', 'PDFT', 'CASPiDFT']

character(*),parameter :: PossibleRDMType(5) = &
[character(8) :: &
'GVB', 'APSG', 'CASSCF', 'DMRG', 'HF']

character(*),parameter :: PossibleDFAType(3) = &
[character(8) :: &
'srLDA', 'srPBE', 'PBE']

character(*),parameter :: PossibleMonomers(2) = &
[character(8) :: 'A', 'B']

character(:), allocatable :: InputPath
!InputPath = "./input.inp"

type CalculationBlock
      integer :: InterfaceType = INTER_TYPE_DAL
      integer :: NBasis = FLAG_NOBASIS
      integer :: JobType = JOB_TYPE_AC
      integer :: Fragments = FLAG_FRAG
      integer :: RDMType = RDM_TYPE_GVB
      integer :: RDMSource = INTER_TYPE_DAL
      integer :: Response = RESP_ERPA
      integer :: DFApp = DF_NONE
      integer :: Kernel = 1
      integer :: TwoMoInt = TWOMO_INCORE
      integer :: Inactive = FLAG_CORE 
      integer :: SymType = TYPE_NO_SYM
      integer :: SaptLevel = SAPTLEVEL2 
      logical :: Restart = FLAG_RESTART
      logical :: PostCAS = FLAG_POSTCAS
      integer :: IPrint  = 0 !FLAG_PRINT_LEVEL 
      double precision :: RPAThresh = 1.0D-6
      integer :: imon = 1
      character(:), allocatable :: JobTitle
      character(:), allocatable :: IntegralsFilePath
end type CalculationBlock

type SystemBlock
      integer :: NoSt = 1
      integer :: NStates = 1
      integer :: Multiplicity = 1
      integer :: Charge = 0
      integer :: ZNucl   = 0
      integer :: NBasis = 0
      integer :: Monomer = MONOMER_A
      integer :: NELE
      double precision :: XELE
      double precision :: PotNuc 
      double precision :: SumOcc = 0
      double precision :: Omega = 1
      integer :: NSym 
      integer :: NSymBas(8),NSymOrb(8)
      integer :: NOrb, NGem
      integer :: NActOrb = 1
      integer :: NAct, INAct
      integer :: NActS(8), INActS(8)
      integer :: NDim, NDimX
      integer :: NDimN
      integer :: NCen = 0
      integer :: UCen = 0  
      integer :: NMonBas(8) = 0
      integer :: IPrint = 0
      integer :: IWarn = 0
      integer :: icnt
      integer :: num0,num1,num2
      integer :: TwoMoInt = TWOMO_INCORE
      logical :: DeclareTwoMo = .false.
      logical :: DeclareSt  = .false.
      logical :: ISHF = .false. 
      logical :: doRSH = .false., SameOm = .true.
      logical :: PostCAS =.false.
      double precision :: ThrAct = 0.992d0
      double precision :: ThrSelAct = 1.d-8 
      integer,allocatable :: InSt(:,:)
      integer,allocatable :: IGem(:), IndAux(:)
      integer,allocatable :: IndX(:), IndN(:,:), IPair(:,:)
      integer,allocatable :: MultpC(:,:),NSymNO(:)
      integer,allocatable :: IndXh(:)
      integer,allocatable :: NumOSym(:),IndInt(:)
      ! TEST ONLY
      integer,allocatable :: IndNT(:,:)
      integer,allocatable :: Ind2(:)
      double precision,allocatable :: Occ(:), CICoef(:)
      double precision,allocatable :: OrbE(:)
      double precision,allocatable :: TwoMO(:)
      double precision,allocatable :: CMO(:,:)
      double precision,allocatable :: WPot(:,:)
      double precision,allocatable :: VCoul(:)
      double precision,allocatable :: RDM2(:),RDM2Act(:,:,:,:)
      double precision  :: charg(maxcen),xyz(maxcen,3)

end type SystemBlock

type FlagsData
! default setting: ERPA-GVB
     ! mainp.f
     integer :: IDALTON = 1
     integer :: IRes    = 0
     integer :: IAO     = 0
     integer :: INO     = 0
     integer :: NoSym   = 1
     integer :: NoSt    = 1
     integer :: IGVB    = 1
     integer :: ITwoEl  = TWOMO_INCORE 
     integer :: IFun    = 13
     integer :: IFunSR    = 0 
     integer :: IFunSRKer = 0
     double precision :: Alpha = 0
     integer :: IModG   = 1
     integer :: NGOcc   = 0
     integer :: ILoc    = 1
     integer :: IFreeze = 0
     integer :: IAPSG   = 1
     integer :: ISERPA  = 0
     integer :: ISAPT   = 0
     integer :: ISHF    = 0
     character(:), allocatable :: JobTitle
     integer :: JobType = 0
     ! initia.f
     integer :: IA = 1
     integer :: ICASSCF = 0
     integer :: IDMRG   = 0  
     ! interpa.f  
     integer :: IFlAC   = 0
     integer :: IFlSnd  = 0
     integer :: IFlCore = 1
     integer :: IFlFrag = 0
     integer :: IFl12 = 1
     ! sapt_main.f90
     integer :: IFlag0 = 0

end type FlagsData

type InputData

     type(CalculationBlock) :: CalcParams
     type(SystemBlock),allocatable :: SystemInput(:)
     integer :: iflag = 0
     type(FlagsData) :: Flags 

end type InputData

type SaptData

     type(SystemBlock) :: monA, monB
     double precision :: Vnn,elst,e2ind,e2disp
     double precision :: e2disp_sc,e2disp_sp,e2disp_unc
     integer :: InterfaceType = INTER_TYPE_DAL
     integer :: SaptLevel = SAPTLEVEL2
     integer :: IPrint = 1000
     logical :: EnChck = .true., HFCheck=.true. 
     logical :: doRSH = .false., SameOm = .true.

end type SaptData

contains 

subroutine free_System(System)
implicit none
type(SystemBlock) :: System 

deallocate(System%InSt)

end subroutine free_System

subroutine free_Input(Input)
implicit none
type(InputData) :: Input

deallocate(Input%SystemInput(1)%InSt)
deallocate(Input%SystemInput)

end subroutine free_Input

subroutine print_Input(Input)
implicit none
type(InputData) :: Input
integer :: switch
integer :: i,imon

write(LOUT,'()')
!write(LOUT,'(8a10)') ('**********',i=1,8)
write(LOUT,'(1x,a)') 'INPUT '
write(LOUT,'(8a10)') ('**********',i=1,8)
!write(LOUT,'(8a10)') ('----------',i=1,8)

associate( CalcParams => Input%CalcParams)
 ! CALCULATION BLOCK
 if(allocated(CalcParams%JobTitle)) then
    write(LOUT,' (1x,a,4x,a)') "JOB TITLE: ", &
                 CalcParams%JobTitle
 else
    write(LOUT,'(1x,a,4x,a)') "JOB TITLE: ", & 
                 "EMPTY"
 endif
 write(LOUT,' (1x,a,4x,a)') "INTERFACE: ", &
              PossibleInterface(CalcParams%InterfaceType)
 write(LOUT,' (1x,a,5x,a)') "JOB TYPE: ",  &
                      PossibleJobType(CalcParams%JOBtype)
 write(LOUT,' (1x,a,5x,a)') "RDM TYPE: ",  &
              PossibleRDMType(CalcParams%RDMType)
 if(CalcParams%DFApp>0) then
    write(LOUT,' (1x,a,5x,a)') "DFA TYPE: ",  &
                       PossibleDFAType(CalcParams%DFApp)
 endif
 write(LOUT,' (1x,a,3x,a)') "RDM SOURCE: ",  &
              PossibleInterface(CalcParams%RDMSource)
! write(LOUT,' (1x,a,6x,i3)') "NBASIS: ",  &
!              CalcParams%NBasis
 if(CalcParams%Fragments==1) then
   write(LOUT,' (1x,a,4x,a)') "FRAGMENTS: ",  &
              "TRUE"
 endif
! write(LOUT, *) "RPA thresh: ", CalcParams%RPAThresh
 if (allocated(CalcParams%IntegralsFilePath)) then
       write(*, *) "Ints file: ", CalcParams%IntegralsFilePath
 end if

  switch = 0
  if(Input%SystemInput(1)%Monomer==2) switch = 3

 ! SYSTEM BLOCK(S)
 do imon=1,CalcParams%imon
    associate(System => Input%SystemInput(abs(imon-switch)) )
      write(LOUT, '()')
      write(LOUT, '(1x,a,6x,a)') "MONOMER: ", &
                    PossibleMonomers(System%Monomer)
      write(LOUT, '(1x,a,6x,i3)') "ZNUCL: ", System%ZNucl
      write(LOUT, '(1x,a,5x,i3)') "CHARGE: ", System%Charge
      if(System%NCen.gt.0) then
         write(LOUT, '(1x,a,i2)') "NO. OF ATOMS: ", System%NCen 
      endif
      if(System%UCen.gt.0) then
         write(LOUT, '(1x,a,i2)') "NO. OF SYM. EQUIV. ATOMS: ", System%UCen 
      endif
      if(System%ISHF) then
         write(LOUT, '(1x,a,l2)') "HARTREE-FOCK: ", System%ISHF 
      endif

    end associate
 enddo

end associate

! FLAGS BLOCK
 if(Input%iflag==1) then
    write(LOUT, '()')
    write(LOUT, '(1x,a,6x,i3)') "IDALTON ", &
                 (Input%Flags%IDALTON)
     write(LOUT, '(1x,a,6x,i3)') "IRes    ", &
                 (Input%Flags%IRes)
     write(LOUT, '(1x,a,6x,i3)') "IAO     ", &
                 (Input%Flags%IAO)
     write(LOUT, '(1x,a,6x,i3)') "INO     ", &
                 (Input%Flags%INO)
     write(LOUT, '(1x,a,6x,i3)') "NoSym   ", &
                 (Input%Flags%NoSym)
     write(LOUT, '(1x,a,6x,i3)') "IGVB    ", &
                 (Input%Flags%IGVB)
     write(LOUT, '(1x,a,6x,i3)') "IFun    ", &
                 (Input%Flags%IFun)
     write(LOUT, '(1x,a,6x,i3)') "IFunSR  ", &
                 (Input%Flags%IFunSR)
     write(LOUT, '(1x,a,5x,i3)') "IFunSRKer", &
                 (Input%Flags%IFunSRKer)
     write(LOUT, '(1x,a,6x,i3)') "IModG   ", &
                 (Input%Flags%IModG)
     write(LOUT, '(1x,a,6x,i3)') "NGOcc   ", &
                 (Input%Flags%NGOcc)
     write(LOUT, '(1x,a,6x,i3)') "ILoc    ", &
                 (Input%Flags%ILoc)
     write(LOUT, '(1x,a,6x,i3)') "IFreeze ", &
                 (Input%Flags%IFreeze)
     write(LOUT, '(1x,a,6x,i3)') "IAPSG   ", &
                 (Input%Flags%IAPSG)
     write(LOUT, '(1x,a,6x,i3)') "ISERPA  ", &
                 (Input%Flags%ISERPA)
     write(LOUT, '(1x,a,6x,i3)') "IA      ", &
                 (Input%Flags%IA)
     write(LOUT, '(1x,a,6x,i3)') "ICASSCF ", &
                 (Input%Flags%ICASSCF)
     write(LOUT, '(1x,a,6x,i3)') "IDMRG   ", &
                 (Input%Flags%IDMRG)
     write(LOUT, '(1x,a,6x,i3)') "IFlAC   ", &
                 (Input%Flags%IFlAC)
     write(LOUT, '(1x,a,6x,i3)') "IFLSnd  ", &
                 (Input%Flags%IFlSnd)
     write(LOUT, '(1x,a,6x,i3)') "IFlCore ", &
                 (Input%Flags%IFlCore)
     write(LOUT, '(1x,a,6x,i3)') "IFlFrag ", &
                 (Input%Flags%IFlFrag)
     write(LOUT, '(1x,a,6x,i3)') "IFl12   ", &
                 (Input%Flags%IFl12)
     write(LOUT, '(1x,a,6x,i3)') "ISAPT   ", &
                 (Input%Flags%ISAPT)
 endif

!write(LOUT,'(8a10)') ('**********',i=1,8)
!write(LOUT,'(1x,a)') 'END INPUT'
write(LOUT,'()') 
!write(LOUT,'(8a10)') ('**********',i=1,8)

end subroutine print_Input

subroutine print_TwoInt(NBasis)
! print trasformed integrals 
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
implicit none

integer,intent(in) :: ndim
double precision,intent(in) :: mat(ndim,ndim)
integer :: i,j

 do i=1,ndim
    write(LOUT,*) i
    write(LOUT,'(10f11.6)') (mat(i,j),j=1,ndim)
 enddo
 write(LOUT,'()') 
 
 return
end subroutine print_sqmat

subroutine print_diag(mat,ndim)
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

! nLine=n/10
! if(nLine*10-n.Ne.0)nLine=nLine+1
! do i=1,n
!    write(*,'(i3)') i
!
!    do ll=0,nLine-1
!       nn=n-10*ll
!       if(nn.le.10) then
!          write(LOUT,'(10f10.6)') (cmo(i,j),j=10*ll+1,n)
!       else
!          write(LOUT,'(10f10.6)') (cmo(i,j),j=10*ll+1,10*(ll+1))
!       endIf
!    enddo
!    write(LOUT,'()')
! enddo

end subroutine print_mo

!HERE!
subroutine readgridmolpro(iunit,text,npt,r,wt)
implicit none
integer :: iunit
integer :: ios
integer :: npt
double precision :: r(:,:),wt(:)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'  
     stop
  endif
  if(label==text) then

     read(iunit) npt
     read(iunit) r(1:3,1:npt) 
     read(iunit) wt(1:npt)
     exit
  endif 
enddo

end subroutine readgridmolpro

subroutine readorbsmolpro(iunit,text,mapinv,orbval)
implicit none
integer :: iunit
integer :: ios
integer :: npt, ndiff, ntg
double precision :: mapinv(:),orbval(:,:,:)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'  
     stop
  endif
  if(label==text) then

     read(iunit) npt,ndiff,ntg
     read(iunit) orbval(1:npt,1:ndiff,1:ntg)
     read(iunit) mapinv(1:ntg)
     exit
  endif 
enddo

end subroutine readorbsmolpro

subroutine read_2rdm_molpro(twordm,nost,nosym,infile,iwarn,nact)
implicit none
! only active part of 2RDM is kept

integer,intent(in) :: nact,nost,nosym
integer,intent(inout) :: iwarn
character(*),intent(in) :: infile
double precision,intent(out) :: twordm(nact**2*(nact**2+1)/2)

integer :: iunit,ios,ist,ic1d,TrSq
integer :: istsym,isym,nstate,nstsym
integer :: i,j,k,l,ij,kl,ik,jl,ijkl,ikjl,lend
double precision,allocatable :: work(:)
logical :: scanfile
character(8) :: label

 !check if any states declared in input
 scanfile = .false.
 if(nosym>0) scanfile = .true.

 TrSq = nact**2*(nact**2+1)/2

 allocate(work(TrSq))
 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 2RDM     not found!'  
              stop
           endif
           if(label=='2RDM    ') then
              !read(iunit) ic1d,nstate
              read(iunit) ic1d,nstsym
              !print*, ic1d,TrSq,nstsym
              if(scanfile) then 
                 do istsym=1,nstsym
                    read(iunit) isym,nstate
                    do i=1,nstate
                       read(iunit) ist 
                       read(iunit) work(1:TrSq)
                       if(ist==nost.and.isym==nosym) exit fileloop
                    enddo
                 enddo
                 write(LOUT,'(1x,a,i2,a,i1,a)') 'ERROR!!! 2RDM FOR STATE',&
                             & nost,'.',nosym,' NOT PRESENT IN 2RDM FILE!'
                 stop
              else
                 read(iunit) isym,nstate 
                 read(iunit) ist
                 read(iunit) work(1:TrSq)
                 if(isym/=1.or.ist/=1) then
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'WARNING! 2RDM FOR STATE',&
                                & ist,'.',isym,' WILL BE USED IN CALCULATIONS!'
                    iwarn = iwarn + 1
                 endif
                 exit fileloop
              endif
           endif 
         enddo fileloop

 close(iunit)

 ijkl = 0
 do i=1,nact
    do j=1,nact
       do k=1,i
          lend = nact
          if(k==i) lend = j
          do l=1,lend
             ijkl = ijkl + 1

             ik = (i-1)*nact + k
             jl = (j-1)*nact + l
             ikjl = max(ik,jl)*(max(ik,jl)-1)/2+min(ik,jl)
             twordm(ikjl) = work(ijkl)*0.5d0

             ik = (k-1)*nact + i
             jl = (l-1)*nact + j
             ikjl = max(ik,jl)*(max(ik,jl)-1)/2+min(ik,jl)
             twordm(ikjl) = work(ijkl)*0.5d0

          enddo
       enddo
    enddo
 enddo

 deallocate(work)

!!! HERE : in SAPT writing to rdm2.dat will be necessary
!!!! better loop!!!
! ijkl = 0
! open(newunit=iunit,file=outfile,form='unformatted')  
! do i=1,nact
!    do j=1,nact
!       ij = (i-1)*nact+j 
!       do k=1,nact
!          do l=1,nact
!             kl = (k-1)*nact+l
!             if(ij>=kl) then
!               ijkl = ijkl + 1  
!               write(iunit,'(4I4,F19.12)') k,i,l,j,twordm(ijkl)
!             endif        
!          enddo
!       enddo
!    enddo
! enddo
! close(iunit)

end subroutine read_2rdm_molpro

subroutine read_nact_molpro(nact,infile)
implicit none

character(*),intent(in) :: infile
integer :: iunit,ios,ist,nact,nact2,nstate
integer :: i,j,ij,idx
character(8) :: label

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1RDM     not found!'
              stop
           endif
           if(label=='1RDM    ') then
              read(iunit) nact,nact2,nstate
              exit fileloop
!              print*, nact,nact2,nstate
           endif
         enddo fileloop

 close(iunit)

end subroutine read_nact_molpro

subroutine read_1rdm_molpro(onerdm,nost,nosym,infile,iwarn,nbasis)
implicit none

integer,intent(in) :: nbasis,nost,nosym
integer,intent(inout) :: iwarn 
character(*),intent(in) :: infile
double precision,intent(out) :: onerdm(nbasis*(nbasis+1)/2)

integer :: iunit,ios,ist,isym,nact,nact2,nstate,nstsym
integer :: i,j,ij,idx,istsym
double precision,allocatable :: work(:)
logical :: scanfile
character(8) :: label

 !check if any states declared in input
 scanfile = .false.
 if(nosym>0) scanfile = .true.

 allocate(work(NBasis**2))
 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1RDM     not found!'  
              stop
           endif
           if(label=='1RDM    ') then
              read(iunit) nact,nact2,nstsym
              !print*, nact2,nstate
              if(scanfile) then
                 do istsym=1,nstsym
                    read(iunit) isym,nstate 
                    do i=1,nstate
                       read(iunit) ist 
                       read(iunit) work(1:nact2)
                       if(ist==nost.and.isym==nosym) exit fileloop
                    enddo
                 enddo
                 write(LOUT,'(1x,a,i2,a,i1,a)') 'ERROR!!! 1RDM FOR STATE',&
                             & nost,'.',nosym,' NOT PRESENT IN 1RDM FILE!'
                 stop
              else
                 read(iunit) isym,nstate 
                 read(iunit) ist
                 read(iunit) work(1:nact2)
                 if(isym/=1.or.ist/=1) then
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'WARNING! 1RDM FOR STATE',&
                                & ist,'.',isym,' WILL BE USED IN CALCULATIONS!'
                    write(LOUT,'()')
                    iwarn = iwarn + 1
                 endif 
                 exit fileloop
              endif
           endif 
         enddo fileloop

 close(iunit)

 !call sq_to_triang(work,onerdm,nact)
 onerdm = 0
 idx = 0 
 do j=1,nact
    do i=1,j
       idx = idx + 1
       ij = i + (j-1)*nact 
       onerdm(idx) = work(ij)*0.5d0
    enddo
 enddo

 deallocate(work)

end subroutine read_1rdm_molpro

subroutine read_mo_molpro(cmo,infile,text,nbasis)
implicit none

integer,intent(in) :: nbasis
character(*),intent(in) :: infile,text
double precision,intent(out) :: cmo(nbasis,nbasis)
character(8) :: label
integer :: iunit,ios
integer :: nsym,nbas(8),offs(8),ncmot
integer :: i,j,idx,irep,ioff 
double precision ::tmp(nbasis**2)

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 !rewind(iunit)
 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL '//text//' not found!'  
      stop
   endif
   if(label==text) then
      read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
      ncmot = sum(nbas(1:nsym)**2)
      !print*, nsym,nbas(1:nsym),offs(1:nsym)
      !print*, ncmot
      read(iunit) tmp(1:ncmot)
      exit
   endif 
 enddo
 
!! HERE: should there be a norb-offset?
 cmo = 0
 idx = 0
 do irep=1,nsym
    ioff = offs(irep)
    do j=1,nbas(irep)
       do i=1,nbas(irep)
          idx = idx + 1 
          cmo(ioff+i,ioff+j) = tmp(idx)
       enddo
    enddo       
 enddo

  !write(LOUT,*) 'test print'
  !do i=1,NBasis
  !   write(*,'(14f11.6)') (cmo(i,j),j=1,nbasis)
  !end do

 close(iunit)

end subroutine read_mo_molpro 

subroutine readlabel(iunit,text)
! sets file pointer 
! to first data after text
implicit none

integer :: iunit
integer :: ios
character(8) :: text, label(4)

rewind(iunit)
do 

  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! Empty section in AOTWOINT!'
     stop
  endif
  if(label(1)=='********') then
     if(label(4)==text) exit
  endif

enddo

end subroutine readlabel

subroutine readoneint(iunit,ints)
implicit none

integer :: iunit
double precision :: ints(:)
integer,parameter :: lbuf = 600
double precision :: buf(lbuf)
integer :: ibuf(lbuf)
integer :: length,i

ints=0
! information are kept in one record in the 
! order: buf, ibuf, length
! buf: integrals, ibuf: int number
do 
   read(iunit) buf,ibuf,length
   if(length.lt.0) exit
   do i=1,length
      ints(ibuf(i)) = buf(i)
   enddo
enddo

end subroutine readoneint

subroutine readoneint_molpro(mone,infile,text,expand,ntr)
implicit none
! reads one-el integrals from a molpro file
! if expand=.true. destroys symmetry in mone

integer,intent(in) :: ntr
logical,intent(in) :: expand
character(*),intent(in) :: infile,text
character(8) :: label
integer :: iunit,ios
integer :: nsym,nbas(8),offs(8),ntdg
integer :: i,j,idx,irep,ioff,ii 
double precision,intent(out) :: mone(ntr)
double precision :: tmp(ntr)

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 ntdg = 0
 rewind(iunit)
 read(iunit) 
 read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
 read(iunit)
 do irep=1,nsym
    ntdg = ntdg + nbas(irep)*(nbas(irep)+1)/2
 enddo

 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL '//text//' not found!'  
      stop
   endif
   if(label==text) then
      read(iunit) tmp(1:ntdg)
      exit
   endif 
 enddo

 close(iunit)

 if(expand) then
 ! expand to a large triangle
    mone = 0
    ii   = 0
    idx  = 0
    ioff = 0   
    do irep=1,nsym
      ioff = offs(irep)
      do j=1,nbas(irep)
        do i=1,j
           idx = idx + 1
           ii = (ioff+j)*(ioff+j-1)/2 + ioff+i
           mone(ii) = tmp(idx)
        enddo
      enddo 
    enddo
 else
 ! return ntdg
    mone = 0
    mone(1:ntdg) = tmp(1:ntdg)
 endif

end subroutine readoneint_molpro

!subroutine writeoneint(iunit,ndim,S,V,H)
!implicit none 
!
!integer :: iunit, ndim
!double precision, dimension(ndim,ndim) :: S, V, H
!
! write(*,*) 'iunit',iunit
! write(iunit) S 
! !write(iunit) 'POTENTAL', V
! !write(iunit) 'ONEHAMIL', H
!
! write(LOUT,'(1x,i3)') 'One-electron matrices written to record:', iunit 
!
!end subroutine writeoneint

subroutine delfile(filename)
implicit none

character(*) :: filename
integer :: iunit
logical :: iexist

 inquire(file=trim(filename),EXIST=iexist)
 if(iexist) then
   open(newunit=iunit,file=trim(filename),status='OLD')
   close(iunit,status='DELETE')
 else
   write(LOUT,'(1x,a)') 'WARNING! NOT POSSIBLE TO DELETE '// filename //' FILE'
 endif

end subroutine delfile

subroutine read2rdm(Mon,NBas)

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
  !print*, Mon%NAct,NRDM2Act

  allocate(Mon%RDM2(NRDM2Act), &
           Mon%RDM2Act(Mon%NAct,Mon%NAct,Mon%NAct,Mon%NAct))
  Mon%RDM2(1:NRDM2Act)=0
  Mon%RDM2Act=0

  open(newunit=iunit,file=rdmfile,status='OLD',&
       form='FORMATTED')
  do

     read(iunit,'(4i4,f19.12)',iostat=ios) i,j,k,l,val

!    val IS DEFINED AS: < E(IJ)E(KL) > - DELTA(J,K) < E(IL) > = 2 GAM2(JLIK)

     if(ios==0) then
        Mon%RDM2(NAddrRDM(j,l,i,k,Mon%NAct))=Half*val
        !print*, 'old',Mon%RDM2(NAddrRDM(j,l,i,k,Mon%NAct))

!       not all elements would be included ?
!        Mon%RDM2Act(j,l,i,k) = Half*val        
        !print*, 'new',Mon%RDM2Act(j,l,i,k)

        i=Ind1(i)
        j=Ind1(j)
        k=Ind1(k)
        l=Ind1(l)

      elseif(ios/=0) then 
        exit

     endif

  enddo
  close(iunit)

 ! not elegant
  do i=1,Mon%NAct
     do j=1,Mon%NAct
        do k=1,Mon%NAct
           do l=1,Mon%NAct
              Mon%RDM2Act(i,j,k,l) = Mon%RDM2(NAddrRDM(i,j,k,l,Mon%NAct))
           enddo
        enddo
     enddo
  enddo

  allocate(Mon%Ind2(NBas))

  Mon%Ind2 = Ind2

end subroutine read2rdm

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

subroutine get_den(nbas,MO,Occ,Fac,Den)
implicit none

integer,intent(in) :: nbas
double precision, intent(in) :: MO(nbas,nbas)
double precision, intent(in) :: Occ(nbas)
double precision, intent(in) :: Fac
double precision, intent(out) :: Den(nbas,nbas)
integer :: i

Den = 0
do i=1,nbas
   call dger(nbas, nbas, Fac*Occ(i), MO(:, i), 1, MO(:, i), 1, Den, nbas)
enddo

end subroutine get_den

subroutine get_one_mat(var,mat,mono,nbas)
implicit none

character(1),intent(in) :: var
integer,intent(in) :: nbas,mono
double precision,intent(out) :: mat(nbas,nbas)
integer :: ione
logical :: valid
character(8) :: label
character(:),allocatable :: onefile

 if(mono==1) then
    onefile = 'ONEEL_A'
 elseif(mono==2) then
    onefile = 'ONEEL_B'
 else
    write(LOUT,'(1x,a)') 'ERROR!!! ONLY 2 MONOMERS ACCEPTED!'
    stop
 endif

 open(newunit=ione,file=onefile,access='sequential',&
      form='unformatted',status='old')

 valid=.false.
 mat=0
 select case(var)
 case('V','v')

    read(ione)
    read(ione) label,mat 
    if(label=='POTENTAL') valid=.true. 
 
 case('S','s')

    read(ione) label,mat 
    if(label=='OVERLAP ') valid=.true. 

 case('H','h')

    read(ione) 
    read(ione)
    read(ione) label,mat 
    if(label=='ONEHAMIL') valid=.true. 

 case default
    write(LOUT,'()')
    write(LOUT,'(1x,a)') 'ERROR IN get_one_max! TYPE '//var//' NOT AVAILABLE!'
    stop
 end select

 if(.not.valid) then
    write(LOUT,'(1x,a)') 'ERROR!!! LABEL MISMATCH IN get_one_mat!' 
    stop
 endif

 close(ione)

end subroutine get_one_mat

subroutine basinfo(nbasis,basfile,intf)
implicit none

character(*),intent(in) :: basfile,intf
integer,intent(out) :: nbasis
integer :: iunit 
integer :: nsym,nbas(8),norb(8),nrhf(8),ioprhf
logical :: ex

 inquire(file=basfile,EXIST=ex)

 if(ex) then 
    open(newunit=iunit,file=basfile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
   
    if(trim(intf)=='DALTON') then
       ! read basis info
       call readlabel(iunit,'BASINFO ')
   
       read (iunit) nsym,nbas,norb,nrhf,ioprhf
       !write(LOUT,*)  nsym,nbas,norb,nrhf,ioprhf

    elseif(trim(intf)=='MOLPRO') then
       read(iunit) nsym
       read(iunit) nbas(1:nsym)
    endif

    close(iunit)
    nbasis = sum(nbas(1:nsym))

 else
    write(LOUT,'(1x,a)') 'WARNING: '// basfile //' NOT FOUND!'
    write(LOUT,'(1x,a)') 'TRYING TO READ NBasis FROM INPUT!'
 endif


end subroutine basinfo

!function trace(m,n) result(tr)
!implicit none
!
!integer,intent(in) :: n
!double precision,intent(in) :: m(n,n)
!integer :: i
!double precision :: tr
!
! tr = 0
! do i=1,n
!    tr = tr + m(i,i)
! enddo
! 
!end function trace 

function iaddr(IAddr1,IAddr2,IAddr3,IAddr4) result(NAddr3)
! POINTER FOR TWO-ELECTRON INTEGRALS
 implicit none
! parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
 integer    :: IAddr1, IAddr2, IAddr3, IAddr4
 integer(8) :: IAddr12, IAddr34
 integer(8) :: NAddr3

! CHANGE THE ORDER IF NECESSARY

 IAddr12 = Max(IAddr1,IAddr2)*(Max(IAddr1,IAddr2)-1)/2 + &
          Min(IAddr2,IAddr1)
 IAddr34 = Max(IAddr3,IAddr4)*(Max(IAddr3,IAddr4)-1)/2 + &
          Min(IAddr3,IAddr4)

! GET THE POSITION OF THE ELEMEMT (12|34)

 NAddr3 = Max(IAddr12,IAddr34)*(Max(IAddr12,IAddr34)-1)/2 + &
          Min(IAddr12,IAddr34)

end function iaddr

function uppercase(s)
      !
      ! Convert characters to uppercase.
      ! Numbers and special characters are ignored.
      !
      character(:), allocatable :: uppercase
      character(*), intent(in)  :: s
      integer :: idx, k
      character(len=*), parameter :: STR_LETTER_UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      character(len=*), parameter :: STR_LETTER_LOWER = "abcdefghijklmnopqrstuvwxyz"

      uppercase = s
      do k = 1, len_trim(s)
            idx = index(STR_LETTER_LOWER, s(k:k))
            if (idx > 0) then
                  uppercase(k:k) = STR_LETTER_UPPER(idx:idx)
            end if
      end do
end function uppercase

subroutine read_statearray(val,inst,instates,delim)

     character(*), intent(in) :: val
     character(1), intent(in) :: delim
     integer,intent(inout) :: instates

     integer :: ii,k
     logical :: dot
     integer,allocatable :: inst(:,:)
     character(:), allocatable :: w,v,s1,s2

     w = trim(adjustl(val))
     v = trim(adjustl(val))
      
     if (len(w) == 0) then
           write(LOUT,'(1x,a)') 'ERROR!!! NO STATES GIVEN FOR Ensamble!'
           stop
     else
           ! check for dots
           k = index(v,'.')
           if(k /= 0) then 
              dot=.true.
           else
              dot=.false.
           endif
         
           ! get number of states 
           instates = 0
           dimloop: do 
                     k = index(v, delim)
                     instates = instates + 1
                     v = trim(adjustl(v(k+1:)))
                     if (k == 0) exit dimloop
                    enddo dimloop

           ! assign states
           allocate(inst(2,instates))
           instates = 0
           arrloop: do 
                     k = index(w, delim)
                     instates = instates + 1
                     if(k /= 0) then
                         if(dot) then
                            call split(w(1:k-1),s1,s2,'.')
                            read(s1, *) inst(1,instates) 
                            read(s2, *) inst(2,instates) 
                         else
                            s1 = w(1:k-1)
                            read(s1, *) inst(1,instates) 
                            inst(2,instates) = 1  
                         endif
                         w = trim(adjustl(w(k+1:)))
                         !print*, '1 2',inst(1,instates),inst(2,instates)
                     elseif (k == 0) then 
                         !print*, 'last ', w
                         if(dot) then
                            call split(w,s1,s2,'.')
                            read(s1, *) inst(1,instates)
                            read(s2, *) inst(2,instates)
                         else
                            s1 = w
                            read(s1, *) inst(1,instates) 
                            inst(2,instates) = 1  
                         endif
                         !print*, '1 2',inst(1,instates),inst(2,instates)
                         exit arrloop
                     endif 
                  enddo arrloop

     end if

end subroutine read_statearray

subroutine split(s, s1, s2, delimiter)
      !
      ! Split a list of words into two pieces:
      ! "keyword    value1 value2" -> "keyword" + "value1 value2".
      !
      character(*), intent(in)               :: s
      character(:), allocatable, intent(out) :: s1
      character(:), allocatable, intent(out) :: s2
      character(1), intent(in), optional :: delimiter

      integer :: k
      character(:), allocatable :: w
      character(1) :: delim

      if (present(delimiter)) then
            delim = delimiter
      else
            delim = " "
      end if
      
      w = trim(adjustl(s))
      if (len(w) == 0) then
            s1 = ""
            s2 = ""
      else
            k = index(w, delim)
            if (k == 0) then
                  s1 = w
                  s2 = ""
            else
                  s1 = w(1:k-1)
                  s2 = trim(adjustl(w(k+1:)))
            end if
      end if
end subroutine split

subroutine io_text_readline(line, u, eof)
      !
      ! Read a line from a text file. The limit for the line
      ! size is MAXCHUNKS * DEFLEN characters (see the code).
      !
      character(:), allocatable, intent(out) :: line
      integer, intent(in)                    :: u
      logical, optional, intent(out)         :: eof

      character(len=80) :: chunk
      character(len=256) :: errmsg
      integer :: s, ios
      integer :: n
      integer, parameter :: maxchunks = 2**10

      line = ""
      if (present(eof)) eof = .false.
      
      lineloop: do n = 1, maxchunks
            read(u, "(A)", advance="NO", size=s, &
                  iostat=ios, iomsg=errmsg) chunk

            if (s > 0) then
                  line = line // chunk(1:s)
            end if

            if (ios == iostat_end) then
                  if (present(eof)) eof = .true.
                  exit lineloop
            else if (ios == iostat_eor) then
                  exit lineloop
            else if (ios .ne. 0) then
                  write(*, *) "COULD NOT READ LINE"
                  write(*, *) trim(errmsg)
                  stop
            end if
      end do lineloop
end subroutine io_text_readline

function isblank(l)
      logical                      :: isblank
      character(len=*), intent(in) :: l

      if (len_trim(l) .eq. 0) then
            isblank = .true.
      else
            isblank = .false.
      end if
end function isblank

function iscomment(s)
      logical                  :: iscomment
      character(*), intent(in) :: s
      
      character(:), allocatable :: sl

      iscomment = .false.
      if (.not. isblank(s)) then
            sl = adjustl(s)
            if (sl(1:1) == "!") then
                  iscomment = .true.
            end if
      end if            
end function iscomment

!subroutine molpro_routines
!implicit none
!
!integer :: NBasis,NGrid
!
!  call basinfo(NBasis,'SIRIUS_A.RST','DALTON')
!  call molprogrid(NGrid,NBasis)
!  print*, 'here-1?'
!  call molprogrid(NGrid,NBasis)
!  print*, 'here-2?'
! stop
!
!! !... idftgra: functional contains grad rho terms (1) , del.2 rho (2)
!! !... ndiff = (idftgra + 1) * (idftgra + 2) * (idftgra + 3) / 6
!!
!! ! get grid points, integration weights and orbitals and its derivatives
!! ! wt array contains weights of the grid points
!! ! orbval array contains basis functions values and its derivatives at npt
!! ! points of the grid
!! ! e.g. orbval(k, 1, mapinv(i)) is chi_i(r_k) - i-th basis function in the k-th
!! ! grid point value
!!
!!double precision, allocatable :: r(:,:),wt(:),mapinv(:),orbval(:,:,:)
!!integer :: i,j,k,igrid
!!integer(8) :: npt, ndiff, ntg
!!
!! open(newunit=igrid,file='GRID',access='sequential',&
!!      form='unformatted',status='old')
!! read(igrid)
!! read(igrid) npt,ndiff,ntg
!! allocate(r(3,npt),wt(npt),mapinv(ntg),orbval(npt,ndiff,ntg))
!!
!! call readgridmolpro(igrid,'GRIDKS  ',npt,r,wt)
!! call readorbsmolpro(igrid,'ORBVAL  ',mapinv,orbval)
!! 
!! close(igrid)
!!
!! deallocate(orbval,mapinv,wt,r)
!!
!end subroutine molpro_routines

subroutine create_symmats(Mon,MO,NBasis)
implicit none
! create MultpC and NSymNO 
! WARNING!!!! THIS PROCEDURE STILL REQUIRES VERIFICATION
! WITH A SAPT LR-PBE JOB!!!

integer,intent(in) :: NBasis
double precision,intent(in) :: MO(NBasis,NBasis)
type(SystemBlock) :: Mon

integer :: i,j,iorb,istart

allocate(Mon%NSymNO(NBasis),Mon%MultpC(15,15))

if(Mon%NSym==1) then
   Mon%MultpC(1,1)=1
else
   do i=1,Mon%NSym
      do j=1,i
         Mon%MultpC(i,j) = ieor(i-1,j-1)+1
         Mon%MultpC(j,i) = Mon%MultpC(i,j)
      enddo
   enddo
endif

Mon%NSymNO(1:NBasis) = 0
istart = 0
do i=1,Mon%NSym
   do j=istart+1,istart+Mon%NumOSym(i)
      do iorb=1,NBasis
         if(abs(MO(j,iorb).gt.1.d-1)) Mon%NSymNO(iorb) = i
      enddo
   enddo
   istart=istart+Mon%NumOSym(i)
enddo

! check
do i=1,Mon%NSym
   j = 0
   do iorb=1,NBasis
      if(Mon%NSymNO(iorb)==i) j = j + 1
   enddo
   if(j/=Mon%NumOSym(i)) then
      write(LOUT,'(1x,a)') 'ERROR IN create_symmats! SYMMETRY OF NO CANNOT BE ESTABLISHED!'
      stop
   endif
enddo

end subroutine create_symmats

subroutine create_ind(infile,NumOSym,IndInt,NSym,NBasis)
implicit none

character(*) :: infile
integer :: NumOSym(15),IndInt(NBasis),NSym,NBasis
integer :: ifile,ios,i,j,k,isym
integer :: NState,NStSym,IOld,INew
integer :: iclos(8),iact(8),nt(8),ivirt(8),istsy(16),nstats(16)
character(8) :: label

 iclos = 0
 ivirt = 0
 iact = 0
 nstats = 0
 istsy = 0
 NumOSym = 0
 open(newunit=ifile,file=infile,access='sequential',&
      form='unformatted',status='old')

 do 
   read(ifile,iostat=ios) label
    if(ios<0) then
      write(6,*) 'ERROR!!! LABEL ISORDK   not found!'  
      stop
   endif
   if(label=='BASINFO ') then
      read(ifile) NSym
      read(ifile) NStSym
      read(ifile) nstats(1:NStSym)
      read(ifile) istsy(1:NStSym)
      read(ifile) iclos(1:NSym) 
      read(ifile) iact(1:NSym) 
      read(ifile) NumOSym(1:NSym) 
      exit
   endif 
 enddo 

 close(ifile)
 ivirt(1:NSym) = NumOSym(1:NSym)-iclos(1:NSym)-iact(1:NSym)

! print*, 'nact',iact(1:NSym)

 if(NSym>1) then
   ! symmetry 
   IOld = 0
   INew = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
      enddo
   enddo
   IOld = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
      enddo
   enddo
   IOld = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
   enddo

 else
   ! nosym
   do i=1,NBasis  
      IndInt(i) = i
   enddo

 endif

end subroutine create_ind 

!subroutine molprogrid(NGrid,NBasis)
!Implicit Real*8 (A-H,O-Z)
!
!integer ::  NBasis, NGrid
!double precision, allocatable :: r(:,:),wt(:),mapinv(:),orbval(:,:,:)
!integer :: i,j,k,igrid
!integer :: npt, ndiff, ntg
!
! print*, ' '
! print*, 'entergin molprogrid' 
! open(newunit=igrid,file='GRID',access='sequential',&
!      form='unformatted',status='old')
! read(igrid) npt,ndiff,ntg
! !If(NBasis.Ne.ntg) Stop 'Fatal Error in molprogrid: NBasis Ne ntg'
! !If(NGrid.Ne.npt) Stop 'Fatal Error in molprogrid: NGrid Ne npt'
! 
! allocate(r(3,npt),wt(npt),mapinv(ntg),orbval(npt,ndiff,ntg))
! !... idftgra: functional contains grad rho terms (1) , del.2 rho (2)
! !... ndiff = (idftgra + 1) * (idftgra + 2) * (idftgra + 3) / 6
!
! ! get grid points, integration weights and orbitals and its derivatives
! ! wt array contains weights of the grid points
! ! orbval array contains basis functions values and its derivatives at npt
! ! points of the grid
! ! e.g. orbval(k, 1, mapinv(i)) is chi_i(r_k) - i-th basis function in the k-th
! ! grid point value
! ! mapinv - orbital mapping array
!
! call readgridmolpro(igrid,'GRIDKS  ',npt,r,wt)
!! Call CpyV(WGrid,wt,npt)
! call readorbsmolpro(igrid,'ORBVAL  ',mapinv,orbval)
!
! close(igrid)
! deallocate(orbval,wt,r,mapinv)
!
!end subroutine molprogrid

end module types
