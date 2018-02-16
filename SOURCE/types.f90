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
integer, parameter :: JOB_TYPE_ERPA  = 3
integer, parameter :: JOB_TYPE_EERPA = 4
integer, parameter :: JOB_TYPE_SAPT  = 5

integer, parameter :: FLAG_FRAG = 0
integer, parameter :: FLAG_CORE = 1
integer, parameter :: FLAG_NOBASIS = 0
logical, parameter :: FLAG_RESTART = .FALSE.
integer, parameter :: FLAG_PRINT_LEVEL = 0

integer, parameter :: RDM_TYPE_GVB  = 1
integer, parameter :: RDM_TYPE_APSG = 2
integer, parameter :: RDM_TYPE_CAS  = 3
integer, parameter :: RDM_TYPE_DMRG = 4

integer, parameter :: MONOMER_A = 1
integer, parameter :: MONOMER_B = 2

character(*),parameter :: PossibleInterface(*) = &
[character(8) :: &
'DALTON', 'MOLPRO', 'OWN']

character(*),parameter :: PossibleJobType(*) = &
[character(8) :: &
'AC', 'AC0', 'ERPA', 'AC1', 'SAPT']

character(*),parameter :: PossibleRDMType(*) = &
[character(8) :: &
'GVB', 'APSG', 'CASSCF', 'DMRG']

character(*),parameter :: PossibleMonomers(*) = &
[character(8) :: 'A', 'B']

character(:), allocatable :: InputPath
!InputPath = "./input.inp"

type CalculationBlock
      integer :: InterfaceType = INTER_TYPE_DAL
      integer :: NBasis = FLAG_NOBASIS
      integer :: JOBType = JOB_TYPE_AC
      integer :: Fragments = FLAG_FRAG
      integer :: RDMType = RDM_TYPE_GVB
      integer :: RDMSource = INTER_TYPE_DAL
      integer :: Inactive = FLAG_CORE 
      integer :: SymType = TYPE_NO_SYM
      logical :: Restart = FLAG_RESTART
      integer :: IPrint  = FLAG_PRINT_LEVEL 
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

end type FlagsData

type InputData

     type(CalculationBlock) :: CalcParams
     type(SystemBlock),allocatable :: SystemInput(:)
     integer :: iflag = 0
     type(FlagsData) :: Flags 

end type InputData

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
 write(LOUT,' (1x,a,4x,a)') "JOB TITLE: ", &
              CalcParams%JobTitle
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
      write(LOUT, '(1x,a,7x,i3)') "ZNUCL: ", System%ZNucl
      write(LOUT, '(1x,a,5x,i3)') "CHARGE: ", System%Charge
      write(LOUT, '(1x,a,i2)') "MULTIPLICITY: ", System%Multiplicity
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
 endif

!write(LOUT,'(8a10)') ('**********',i=1,8)
!write(LOUT,'(1x,a)') 'END INPUT'
write(LOUT,'()') 
!write(LOUT,'(8a10)') ('**********',i=1,8)

end subroutine print_Input


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
