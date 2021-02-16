module inputfill

use types
implicit none

contains

subroutine read_Input(Input)
implicit none

type(InputData) :: Input
character(:), allocatable :: InputPath

 InputPath = "./input.inp"
 call sapt_scan_inputfile(InputPath, Input%CalcParams)
 call init_Input(InputPath, Input)
 call read_inputfile(InputPath, Input)

 call check_Input(Input)
 call print_Input(Input)

end subroutine read_Input

subroutine init_Input(filename,Input)
implicit none
character(len=*), intent(in) :: filename
type(InputData), intent(inout) :: Input
logical :: EndOfFile
integer :: u
character(:), allocatable :: line
integer :: isys
integer :: current_block
integer, parameter :: block_none    = 0
integer, parameter :: block_system  = 2
integer :: imon

 imon = Input%CalcParams%imon

 !write(LOUT,'(1x,a,i2,1x,a)') 'CHECK: Allocated', imon, 'monomers!'
 allocate(Input%SystemInput(imon))

 open(newunit=u, file=filename, status="old", &
         access="sequential", position="rewind")

 ! check monomers
 select case(Input%CalcParams%JobType)
 ! SAPT case
 case(5)

    current_block = block_none
    isys = 0
    lines: do

             call io_text_readline(line, u, EndOfFile)

             if (EndOfFile) then
                 exit lines
             end if
             if (isblank(line) .or. iscomment(line)) then
                 cycle lines
             end if

             select case (uppercase(line))
             case("SYSTEM")
                current_block = block_system
                cycle lines
             end select

             if(current_block==block_system) then
                call read_sapt_mon(Input,line,isys)
             endif

    enddo lines

    ! check two monomers
    if(isys==1) then
       write(LOUT,'(1x,a)') 'ERROR! DEFINE TWO MONOMERS FOR SAPT!'
       stop
    endif

    ! check different monomers
    associate( monA => Input%SystemInput(1)%Monomer, &
               monB => Input%SystemInput(2)%Monomer )
      if(monA==monB) then
         write(LOUT,'(1x,a)') "ERROR! TWO MONOMERS TYPE "&
                     & // trim(PossibleMonomers(monA)) //" DEFINED&
                     & FOR SAPT!"
         stop
      endif
    end associate

 case default

    current_block = block_none
    isys = 0
    lines2: do

             call io_text_readline(line, u, EndOfFile)

             if (EndOfFile) then
                 exit lines2
             end if
             if (isblank(line) .or. iscomment(line)) then
                 cycle lines2
             end if
             !
             select case (uppercase(line))
             case("SYSTEM")
                current_block = block_system
                isys = isys + 1
                cycle lines2
             end select

    enddo lines2

    if(isys>1) then
       write(LOUT,'(1x,a)') "ERROR! TOO MANY SYSTEM BLOCKS!"
       write(LOUT,'(1x,a)') "ERROR! TOO MANY SYSTEM BLOCKS ???!"
       stop
    endif

 end select

 close(u)

end subroutine init_Input

subroutine read_inputfile(filename, Input)
implicit none
character(len=*), intent(in)   :: filename
type(InputData), intent(inout) :: Input

logical :: EndOfFile
integer :: u
integer :: isys, iflag
integer :: current_block
integer, parameter        :: block_none        = 0
integer, parameter        :: block_calculation = 1
integer, parameter        :: block_system      = 2
integer, parameter        :: block_flags       = 3
character(:), allocatable :: line

 open(newunit=u, file=filename, status="old", &
       access="sequential", position="rewind")
 !
 ! Loop goes over all all lines of the input file
 !
 isys  = 0
 current_block = block_none
 lines: do
       call io_text_readline(line, u, EndOfFile)

       if (EndOfFile) then
             exit lines
       end if
       !
       ! Check if the current line is blank or is a comment
       !
       if (isblank(line) .or. iscomment(line)) then
             cycle lines
       end if
       !
       ! Check if the current line is a block start/block end
       !
       select case (uppercase(line))
       case ("SYSTEM")
             isys = isys + 1
             current_block = block_system
             cycle lines

       case ("CALCULATION")
             current_block = block_calculation
             cycle lines

       case ("FLAGS")
             Input%iflag = Input%iflag + 1
             current_block = block_flags
             cycle lines

       case ("END")
             current_block = block_none
             cycle lines

       end select

       if (current_block == block_system) then
             call read_block_system(Input%SystemInput(isys), line)
       else if (current_block == block_calculation) then
             call read_block_calculation(Input%CalcParams, line)
       else if (current_block == block_flags) then
             call read_block_flags(Input%Flags, line)
       end if
 end do lines

end subroutine read_inputfile

subroutine sapt_scan_inputfile(filename, CalcParams)
! SEARCH INPUT FOR SAPT KEYWORD
implicit none
character(len=*), intent(in)        :: filename
type(CalculationBlock), intent(out) :: CalcParams

integer :: u
integer :: isys
integer :: current_block
integer, parameter :: block_none        = 0
integer, parameter :: block_calculation = 1
integer, parameter :: block_system      = 2
logical            :: EndOfFile
character(:), allocatable :: line

open(newunit=u, file=filename, status="old", &
      access="sequential", position="rewind")
!
! Loop goes over all all lines of the input file
!
 current_block = block_none
 lines: do

          call io_text_readline(line, u, EndOfFile)

          if (EndOfFile) then
              exit lines
          end if
          !
          ! Check if the current line is blank or is a comment
          !
          if (isblank(line) .or. iscomment(line)) then
              cycle lines
          end if
          !
          ! Check for SAPT keyword in input
          !
          select case (uppercase(line))
          case ("CALCULATION")
             current_block = block_calculation
             cycle lines

          end select

          if(current_block==block_calculation) then
             call read_sapt_val(CalcParams, line)
          endif

 enddo lines
close(u)

end subroutine sapt_scan_inputfile

subroutine read_block_calculation(CalcParams, line)
      type(CalculationBlock), intent(inout) :: CalcParams
      character(*), intent(in) :: line

      character(:), allocatable :: key, val

      call split(line, key, val)
      select case (uppercase(key))

      case ("INTERFACE")
           if (uppercase(val) == "DALTON") then
               CalcParams%InterfaceType = INTER_TYPE_DAL
           elseif (uppercase(val) == "MOLPRO") then
               CalcParams%InterfaceType = INTER_TYPE_MOL
               CalcParams%RDMSource = INTER_TYPE_MOL
           elseif (uppercase(val) == "OWN".or.&
                   uppercase(val) == "NONE") then
               CalcParams%InterfaceType = INTER_TYPE_OWN
           elseif (uppercase(val) == "ORCA") then
               CalcParams%InterfaceType = INTER_TYPE_ORCA
               CalcParams%RDMSource = INTER_TYPE_ORCA
               CalcParams%RDMType   = RDM_TYPE_DMRG
           endif

      case ("JOBTYPE")
           if (uppercase(val) == "AC" ) then
               CalcParams%JobType = JOB_TYPE_AC
           elseif (uppercase(val) == "AC0" ) then
               CalcParams%JobType = JOB_TYPE_AC0
           elseif (uppercase(val) == "ERPA" ) then
               CalcParams%JobType = JOB_TYPE_ERPA
           elseif (uppercase(val) == "AC1" ) then
               CalcParams%JobType = JOB_TYPE_ERPA
           elseif (uppercase(val) == "EERPA".or. &
                   uppercase(val) == "ERPA-2") then
               CalcParams%JobType = JOB_TYPE_EERPA
           elseif (uppercase(val) == "EERPA-OLD".or. &
                   uppercase(val) == "ERPA-1") then
               CalcParams%JobType = JOB_TYPE_EERPA_OLD
           elseif (uppercase(val) == "SAPT" ) then
               CalcParams%JobType = JOB_TYPE_SAPT
           elseif (uppercase(val) == "PDFT" ) then
               CalcParams%JobType = JOB_TYPE_PDFT
           elseif (uppercase(val) == "CASPIDFT" ) then
               CalcParams%JobType = JOB_TYPE_CASPIDFT
           elseif (uppercase(val) == "CASPIDFTOPT" ) then
               CalcParams%JobType = JOB_TYPE_CASPIDFTOPT
           elseif (uppercase(val) == "AC0D" ) then
               CalcParams%JobType = JOB_TYPE_AC0D
           elseif (uppercase(val) == "AC0DNOSYMM" ) then
               CalcParams%JobType = JOB_TYPE_AC0DNOSYMM
           elseif (uppercase(val) == "NLOCCORR" ) then
               CalcParams%JobType = JOB_TYPE_NLOCCORR
           endif

     !case ("FRAGMENTS")
     !     if (uppercase(val) == ".TRUE.".or. &
     !         uppercase(val) == "TRUE".or.   &
     !         uppercase(val) == "T") then
     !         CalcParams%Fragments = 1
     !     endif

      case ("CORE")
         read(val, *) CalcParams%Core

      case ("NBASIS")
         read(val, *) CalcParams%NBasis

      case ("RDMTYPE")
           if (uppercase(val) == "GVB") then
              CalcParams%RDMType = RDM_TYPE_GVB
           elseif (uppercase(val) == "APSG" ) then
              CalcParams%RDMType = RDM_TYPE_APSG
           elseif (uppercase(val) == "CASSCF".or.&
                 & uppercase(val) == "CAS") then
              CalcParams%RDMType = RDM_TYPE_CAS
           elseif (uppercase(val) == "HF".or.    &
                 & uppercase(val) == "HFOCK".or. &
                 & uppercase(val) == "HARTREE-FOCK") then
              CalcParams%RDMType = RDM_TYPE_HF
           elseif (uppercase(val) == "DMRG" ) then
              CalcParams%RDMType = RDM_TYPE_DMRG
           endif

      case ("TWOMOINT")
         if(uppercase(val) == "FULL".or. &
            uppercase(val) == "FFFF") then

            CalcParams%TwoMoInt = TWOMO_FFFF

         elseif(uppercase(val) == 'FOFO') then

            CalcParams%TwoMoInt = TWOMO_FOFO

         elseif(uppercase(val) == 'INCORE'.or. &
                uppercase(val) == 'IN-CORE') then

            CalcParams%TwoMoInt = TWOMO_INCORE
         endif

      case ("REDVIRT")
           if (uppercase(val) == ".TRUE.".or. &
               uppercase(val) == "TRUE".or.   &
               uppercase(val) == "T") then
               CalcParams%RedVirt = 1
           endif

      ! here not sure
      case ("RESPONSE")
           if (uppercase(val) == "ERPA-APSG".or.&
               uppercase(val) == "ERPA") then
               CalcParams%Response = RESP_ERPA
           elseif (uppercase(val) == "TD-APSG".or.&
                   uppercase(val) == "FULL") then
               CalcParams%Response = RESP_APSG
           elseif (uppercase(val) == "TD-KS".or.&
                   uppercase(val) == "DFT") then
               CalcParams%Response = RESP_DFT
           endif

      case ("DFUNC")
           if (uppercase(val) == "SRLDA".or.&
               uppercase(val) == "SR-LDA") then
               CalcParams%DFApp = DF_SRLDA
           elseif (uppercase(val) == "SRPBE".or.&
                   uppercase(val) == "SR-PBE") then
               CalcParams%DFApp = DF_SRPBE
           elseif (uppercase(val) == "PBE") then
               CalcParams%DFApp = DF_PBE
           endif

      case ("KERNEL")
           read(val,*) CalcParams%Kernel

      case ("RDMSOURCE")
           if (uppercase(val) == "DALTON") then
              CalcParams%RDMSource = INTER_TYPE_DAL
           elseif (uppercase(val) == "OWN" ) then
              CalcParams%RDMSource = INTER_TYPE_OWN
           endif

      case ("SYMMETRY")
           if (uppercase(val) == "NOSYM") then
              CalcParams%SymType = TYPE_NO_SYM
           elseif (uppercase(val) == "SYM" ) then
              CalcParams%SymType = TYPE_SYM
           endif

      case ("POSTCAS")
           read(val,*) CalcParams%PostCAS

      case ("SAPTLEVEL")
           if (uppercase(val) == "0".or.&
               uppercase(val) == "SAPT0") then
              CalcParams%SaptLevel = SAPTLEVEL0
           elseif (uppercase(val) == "1" ) then
              CalcParams%SaptLevel = SAPTLEVEL1
           elseif (uppercase(val) == "2".or.&
                   uppercase(val) == "SAPT2" ) then
              CalcParams%SaptLevel = SAPTLEVEL2
           elseif (uppercase(val) == "C6") then
              CalcParams%SaptLevel = SAPTLEVEL2
              CalcParams%vdWCoef = 1
           elseif (uppercase(val) == "DISP-CAS") then
              CalcParams%SaptLevel = 10
           endif

      case("RESTART")
           if (uppercase(val) == "TRUE".or.  &
               uppercase(val) == ".TRUE.".or.&
               uppercase(val) == "T") then
               CalcParams%Restart = .TRUE.
           endif

      case ("RPATHRESH")
            read(val, *) CalcParams%RPAThresh

      case ("JOBTITLE")
            CalcParams%JobTitle = val

      case ("INTEGRALSFILEPATH")
            CalcParams%IntegralsFilePath = val

      case ("IPRINT")
            read(val,*) CalcParams%IPrint

      end select
end subroutine read_block_calculation

subroutine read_block_system(SystemParams, line)
implicit none
type(SystemBlock), intent(inout) :: SystemParams

character(*), intent(in)  :: line
character(:), allocatable :: key, val
character(:), allocatable :: first, last

 call split(line, key, val)
 select case (uppercase(key))
 case("STATE")
       SystemParams%DeclareSt = .true.
       call read_statearray(val,SystemParams%InSt,SystemParams%NStates,',')

 case("TRDM")
       SystemParams%DeclareTrSt = .true.
       call read_trstatearray(val,SystemParams%InTrSt,',')
! case ("STATE")
!      read(val, *) SystemParams%NoSt

 case ("NACTIVE")
       read(val, *) SystemParams%NAct
       SystemParams%NActFromRDM = .false.

 case ("CHARGE")
       read(val, *) SystemParams%Charge

 case ("ZNUCL")
       read(val, *) SystemParams%ZNucl

 case ("NATOMS")
       read(val, *) SystemParams%NCen

 case ("ACALPHA")
       read(val, *) SystemParams%ACAlpha

 case ("OMEGA")
       read(val, *) SystemParams%Omega

 case ("EIGFCI")
       read(val, *) SystemParams%EigFCI

 case ("UATOMS")
       read(val, *) SystemParams%UCen
! maybe sth more fancy i.e. swapping monomers
      ! call get_ncen(val,SystemParams)
 case ("THRACT")
       read(val,*) SystemParams%ThrAct

 case ("THRSELACT")
       read(val,*) SystemParams%ThrSelAct

 case ("THRQVIRT")
       read(val,*) SystemParams%ThrQVirt

 case ("THRVIRT")
       read(val,*) SystemParams%ThrVirt

 case ("POSTCAS")
       read(val,*) SystemParams%PostCAS

 case ("TWOMOINT")

      SystemParams%DeclareTwoMo = .true.

      if(uppercase(val) == "FULL".or. &
         uppercase(val) == "FFFF") then

         SystemParams%TwoMoInt = TWOMO_FFFF

      elseif(uppercase(val) == 'FOFO') then

         SystemParams%TwoMoInt = TWOMO_FOFO

      elseif(uppercase(val) == 'INCORE'.or. &
             uppercase(val) == 'IN-CORE') then

         SystemParams%TwoMoInt = TWOMO_INCORE
      endif

 case ("ISHF")
       read(val,*) SystemParams%ISHF

 case ("CUBIC")
       read(val,*) SystemParams%Cubic

 end select
end subroutine read_block_system

subroutine read_block_flags(Flags, line)

type(FlagsData), intent(inout) :: Flags
character(*), intent(in)       :: line

character(:), allocatable :: key, val

      call split(line, key, val)
      select case (uppercase(key))

      case ("IDALTON")
         read(val, *) Flags%IDALTON

      case ("IRES")
         read(val, *) Flags%IRes

      case ("IAO")
         read(val, *) Flags%IAO

      case ("INO")
         read(val, *) Flags%INO

      case ("NOSYM")
         read(val, *) Flags%NoSym

      case ("IGVB")
         read(val, *) Flags%IGVB

      case ("IFUN")
         read(val, *) Flags%IFun

      case ("IFUNSR")
         read(val, *) Flags%IFunSR

      case ("IFUNSRKER")
         read(val, *) Flags%IFunSRKer

      case ("IMODG")
         read(val, *) Flags%IModG

      case ("NGOCC")
         read(val, *) Flags%NGOcc

      case ("ILOC")
         read(val, *) Flags%ILoc

      case ("IFREEZE")
         read(val, *) Flags%IFreeze

      case ("IAPSG")
         read(val, *) Flags%IAPSG

      case ("ISERPA")
         read(val, *) Flags%ISERPA

      case ("IA")
         read(val, *) Flags%IA

      case ("ICASSCF")
         read(val, *) Flags%ICASSCF

      case ("IDMRG")
         read(val, *) Flags%IDMRG

      case ("IFLAC")
         read(val, *) Flags%IFlAC

      case ("IFLSND")
         read(val, *) Flags%IFlSnd

      case ("IFLCORE")
         read(val, *) Flags%IFlCore

      case ("IFLFRAG")
         read(val, *) Flags%IFlFrag1

      case ("IFL12")
         read(val, *) Flags%IFl12

      end select

end subroutine read_block_flags

subroutine read_sapt_val(CalcParams, line)
implicit none

type(CalculationBlock), intent(inout) :: CalcParams
character(*), intent(in)  :: line
character(:), allocatable :: key, val

 call split(line, key, val)
 select case (uppercase(key))
  case ("JOBTYPE")
      if (uppercase(val) == "SAPT" ) then
          CalcParams%JobType = JOB_TYPE_SAPT
          CalcParams%imon = 2
      endif
 end select

end subroutine read_sapt_val 

subroutine read_sapt_mon(Input, line, isys)
implicit none

type(InputData), intent(inout) :: Input
character(*), intent(in)       :: line
integer                        :: isys
character(:), allocatable      :: key, val

 call split(line, key, val)
 select case (uppercase(key))
  case ("MONOMER")
     isys = isys + 1
     if(uppercase(val)=="A") then
         Input%SystemInput(isys)%Monomer = MONOMER_A
     elseif(uppercase(val)=="B") then
         Input%SystemInput(isys)%Monomer = MONOMER_B
     endif
 end select

end subroutine read_sapt_mon

!subroutine get_ncen(val,System)
!
!type(SystemBlock) :: System
!character(*) :: val
!character(:), allocatable :: first, last
!integer :: ncen, AtBeg, AtEnd
!
!call split(val,first,last,":")
!read(first,*) AtBeg
!read(last,*) AtEnd
!ncen = abs(AtEnd - AtBeg) + 1
!System%NCen = ncen
!
!end subroutine get_ncen

subroutine check_Input(Input)
implicit none

type(InputData) :: Input
integer :: imon

 if(Input%iflag.gt.1) then
    write(LOUT,'(1x,a)') 'ERROR! PLACE ALL FLAGS IN ONE BLOCK!'
    stop
 endif

 if(Input%CalcParams%InterfaceType.ne.INTER_TYPE_DAL) then
    if(Input%CalcParams%NBasis==0) then
       !write(LOUT,'(1x,a)') 'FATAL ERROR: NBasis ENTRY MISSING'
       write(LOUT,'(1x,a)') 'NBasis WILL BE READ FROM MOLPRO'
       !stop
    endif
 elseif(Input%CalcParams%NBasis.lt.0) then
     write(LOUT,'(1x,a)') 'FATAL ERROR: INCORRECT ENTRY&
                 & NBasis IN THE INPUT FILE'
     write(LOUT, '(1x,a,3x,i3)') 'NBasis', Input%CalcParams%NBasis
     stop
 endif

 do imon=1,Input%CalcParams%imon
    associate(System => Input%SystemInput(imon))

      if(System%ZNucl==0) then
         write(LOUT,'(1x,a)') 'FATAL ERROR: Znucl ENTRY MISSING&
                    & OR ZERO!'
         stop
      elseif(System%ZNucl.lt.0) then
         write(LOUT,'(1x,a)') 'FATAL ERROR: INCORRECT ENTRY&
                     & ZNucl IN THE INPUT FILE'
         write(LOUT, '(1x,a,3x,i3)') 'ZNucl', System%ZNucl
         stop
      elseif(System%NCen==0.and.Input%CalcParams%imon.gt.1) then
         write(LOUT,'(1x,a)') 'FATAL ERROR: NAtoms ENTRY MISSING!&
                    & HAS TO BE GIVEN FOR SAPT!'
         write(LOUT, '(1x,a,3x,i3)') 'NCen', System%NCen
         stop
      endif
    
      if(System%UCen==0) then
        System%UCen = System%NCen
      elseif(System%UCen.lt.0) then
         write(LOUT,'(1x,a)') 'FATAL ERROR: UAtoms .LT. 0!'
         stop
      endif

      if(.not.System%DeclareSt) then
         allocate(System%InSt(2,1))
      endif
      !write(LOUT, '(1x,a,6x,a)') "Monomer: ", &
      !              PossibleMonomers(System%Monomer)
    end associate
 enddo

end subroutine check_Input

end module inputfill
