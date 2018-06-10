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

integer, parameter :: MONOMER_A = 1
integer, parameter :: MONOMER_B = 2

integer, parameter :: RESP_ERPA = 1
integer, parameter :: RESP_APSG = 2

integer,parameter :: maxcen = 500

character(*),parameter :: PossibleInterface(3) = &
[character(8) :: &
'DALTON', 'MOLPRO', 'OWN']

character(*),parameter :: PossibleJobType(5) = &
[character(8) :: &
'AC', 'AC0', 'ERPA', 'AC1', 'SAPT']

character(*),parameter :: PossibleRDMType(5) = &
[character(8) :: &
'GVB', 'APSG', 'CASSCF', 'DMRG', 'HF']

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
      integer :: Inactive = FLAG_CORE 
      integer :: SymType = TYPE_NO_SYM
      logical :: Restart = FLAG_RESTART
      integer :: IPrint  = 0 !FLAG_PRINT_LEVEL 
      double precision :: RPAThresh = 1.0D-6
      integer :: imon = 1
      character(:), allocatable :: JobTitle
      character(:), allocatable :: IntegralsFilePath
end type CalculationBlock

type SystemBlock
      integer :: Multiplicity = 1
      integer :: Charge = 0
      integer :: ZNucl   = 0
      integer :: NBasis = 0
      integer :: Monomer = MONOMER_A
      integer :: NELE
      double precision :: XELE
      double precision :: PotNuc 
      double precision :: SumOcc = 0
      integer :: NSym 
      integer :: NSymBas(8),NSymOrb(8)
      integer :: NOrb, NGem
      integer :: NAct, INAct
      integer :: NDim, NDimX
      integer :: NDimN
      integer :: NCen = 0
      integer :: UCen = 0  
      integer :: NMonBas(8) = 0
      integer :: IPrint = 0
      integer :: IWarn = 0
      integer :: icnt
      integer :: num0,num1,num2
      logical :: ISHF = .false. 
      double precision :: ThrAct = 0.992d0
      integer,allocatable :: IGem(:), IndAux(:)
      integer,allocatable :: IndX(:), IndN(:,:), IPair(:,:)
      integer,allocatable :: Ind2(:)
      double precision,allocatable :: Occ(:), CICoef(:)
      double precision,allocatable :: OrbE(:)
      double precision,allocatable :: CMO(:,:)
      double precision,allocatable :: WPot(:,:)
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
     integer :: IGVB    = 1
     integer :: IFun    = 13
     integer :: IFunSR    = 0 
     integer :: IFunSRKer = 0
     integer :: IModG   = 1
     integer :: NGOcc   = 0
     integer :: ILoc    = 1
     integer :: IFreeze = 0
     integer :: IAPSG   = 1
     integer :: ISERPA  = 0
     integer :: ISAPT   = 0
     integer :: ISHF    = 0
     character(:), allocatable :: JobTitle
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
     double precision :: Vnn,elst,e2ind,e2disp,e2disp_unc
     integer :: IPrint = 1000
     logical :: EnChck = .true.

end type SaptData

contains 

subroutine free_Input(Input)
implicit none
type(InputData) :: Input

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
 write(LOUT,' (1x,a,3x,a)') "RDM SOURCE: ",  &
              PossibleInterface(CalcParams%RDMSource)
 write(LOUT,' (1x,a,6x,i3)') "NBASIS: ",  &
              CalcParams%NBasis
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
      write(LOUT, '(1x,a,i2)') "MULTIPLICITY: ", System%Multiplicity
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
  print*, Mon%NAct,NRDM2Act

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
do i = 1,nbas
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

subroutine basinfo(nbasis,basfile)
implicit none

character(*),intent(in) :: basfile
integer,intent(out) :: nbasis
integer :: iunit 
integer :: nsym,nbas(8),norb(8),nrhf(8),ioprhf
logical :: ex

 inquire(file=basfile,EXIST=ex)

 if(ex) then 
    open(newunit=iunit,file=basfile,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
   
    ! read basis info
    call readlabel(iunit,'BASINFO ')
   
    read (iunit) nsym,nbas,norb,nrhf,ioprhf
    !write(LOUT,*)  nsym,nbas,norb,nrhf,ioprhf

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

end module types
