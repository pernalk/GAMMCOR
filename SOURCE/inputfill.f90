module inputfill

use types
use print_units

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
integer, parameter        :: block_cholesky    = 4
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

       case ("CHOLESKYBLOCK")
             current_block = block_cholesky
             cycle lines

       case ("FLAGS")
             Input%iflag = Input%iflag + 1
             current_block = block_flags
             cycle lines

       case("XYZ")
             current_block = block_none
             cycle lines

       case ("END")
             current_block = block_none
             cycle lines

       end select

       if (current_block == block_system) then
             call read_block_system(Input%SystemInput(isys), line)
       else if (current_block == block_calculation) then
             call read_block_calculation(Input%CalcParams, line)
       else if (current_block == block_cholesky) then
             call read_block_cholesky(Input%CholeskyParams, line)
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

subroutine read_block_cholesky(CholeskyParams, line)
      type(CholeskyBlock), intent(inout) :: CholeskyParams
      character(*), intent(in) :: line

      character(:), allocatable :: key, val

      call split(line, key, val)
      select case (uppercase(key))

      case ("CHOLESKY")
           if (uppercase(val) == "BINARY".or. &
               uppercase(val) == "BIN"   .or. &
               uppercase(val) == "B") then
              CholeskyParams%Cholesky    = 1
              CholeskyParams%CholeskyBIN = 1
              CholeskyParams%CholeskyOTF = 0
           elseif (uppercase(val) == "ON-THE-FLY".or. &
                   uppercase(val) == "OTF"   .or. &
                   uppercase(val) == "O") then
              CholeskyParams%Cholesky    = 1
              CholeskyParams%CholeskyBIN = 0
              CholeskyParams%CholeskyOTF = 1
           else
              stop "Unknown keyword for Cholesky!"
           endif

      case ("CHOL_ACCU","CHOL_ACCURACY","CHOLESKY_ACCU","CHOLESKY_ACCURACY")
           if (uppercase(val) == "DEFAULT" .or. &
               uppercase(val) == "D" ) then
              CholeskyParams%CholeskyAccu = CHOL_ACCU_DEFAULT
           elseif (uppercase(val) == "TIGHT" .or. &
                   uppercase(val) == "T" ) then
              CholeskyParams%CholeskyAccu = CHOL_ACCU_TIGHT
           elseif (uppercase(val) == "LUDICROUS" .or. &
                   uppercase(val) == "L" ) then
              CholeskyParams%CholeskyAccu = CHOL_ACCU_LUDICROUS
           endif

      case ("H0TEST")
           if (uppercase(val) == ".FALSE.".or. &
               uppercase(val) == "FALSE".or.   &
               uppercase(val) == "F") then
               CholeskyParams%H0test = 0
           endif

      end select

end subroutine read_block_cholesky

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
           elseif (uppercase(val) == "AC0DP" ) then
               CalcParams%JobType = JOB_TYPE_AC0DP
           elseif (uppercase(val) == "ACFREQ" ) then
               CalcParams%JobType = JOB_TYPE_ACFREQ
           elseif (uppercase(val) == "ACFREQNTH" ) then
               CalcParams%JobType = JOB_TYPE_ACFREQNTH
           elseif (uppercase(val) == "AC1FREQNTH" ) then
               CalcParams%JobType = JOB_TYPE_AC1FREQNTH
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

      case ("ORBRELAX")
           if (uppercase(val) == ".FALSE.".or. &
               uppercase(val) == "FALSE".or.   &
               uppercase(val) == "F") then
               CalcParams%OrbRelax = 0
           endif

      case ("ORBINCL", "ORBINCLUDE")
           if (uppercase(val) == ".TRUE.".or. &
               uppercase(val) == "TRUE".or.   &
               uppercase(val) == "T") then
               CalcParams%OrbIncl = 1
           endif

      case ("RDM2APP")
           if (uppercase(val) == "HF".or. &
               uppercase(val) == "HARTREE-FOCK") then
               CalcParams%Rdm2Type = 0
           elseif (uppercase(val) == "BB".or. &
               uppercase(val) == "DMFT") then
               CalcParams%Rdm2Type = 1
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

       case ("TRIPLET")
            if (uppercase(val) == "TRUE".or.  &
                uppercase(val) == ".TRUE.".or.&
                uppercase(val) == "T") then
               CalcParams%Triplet = .TRUE.
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

      case("MAX_CN")
             read(val,*) CalcParams%Max_Cn

      case ("CALPHA")
             read(val,*) CalcParams%CAlpha

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
           elseif (uppercase(val) == "RS".or. &
                   uppercase(val) == "RSPT2") then
              CalcParams%SaptLevel = 999
           elseif (uppercase(val) == "RS+".or. &
                   uppercase(val) == "RSPT2+") then
              CalcParams%SaptLevel = 666
           endif

      case("SAPTEXCH")
           if (uppercase(val) == "DMFT") then
               CalcParams%SaptExch = 1
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

      case ("BASIS")
            CalcParams%BasisSet = val

      case ("BASISPATH")
            CalcParams%BasisSetPath = val

      case("MEMSORT")
            call read_memsrt(val,CalcParams%MemVal,CalcParams%MemType)

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

 case("SPIN")
      SystemParams%DeclareSpin = .true.
      ! spin is kept in Molpro convention: 2*ms
      if(uppercase(val) == "SINGLET" .or. &
         uppercase(val)=="0") then
         SystemParams%ISpinMs2 = 0
      elseif(uppercase(val) == "DOUBLET" .or. &
             uppercase(val) == "1") then
         SystemParams%ISpinMs2 = 1
      elseif(uppercase(val) == "TRIPLET" .or. &
             uppercase(val) == "2") then
         SystemParams%ISpinMs2 = 2
      elseif(uppercase(val) == "QUARTET" .or. &
             uppercase(val) == "3") then
         SystemParams%ISpinMs2 = 3
      elseif(uppercase(val) == "QUINTET" .or. &
             uppercase(val) == "4") then
         SystemParams%ISpinMs2 = 4
      elseif(uppercase(val) == "SEXTET" .or. &
             uppercase(val) == "5") then
         SystemParams%ISpinMs2 = 5
      elseif(uppercase(val) == "SEPTET" .or. &
             uppercase(val) == "6") then
         SystemParams%ISpinMs2 = 6
      endif

 case("TRDM")
       SystemParams%DeclareTrSt = .true.
       call read_trstatearray(val,SystemParams%InTrSt,',')
! case ("STATE")
!      read(val, *) SystemParams%NoSt

 case ("MONOMER")
     if(uppercase(val)=="A") then
         SystemParams%Monomer = MONOMER_A
     elseif(uppercase(val)=="B") then
         SystemParams%Monomer = MONOMER_B
     elseif(uppercase(val)=="AB" .or. &
            uppercase(val)=="DIMER") then
         SystemParams%Monomer = DIMER_AB
     endif

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

 case ("PERVIRT")
       read(val, *) SystemParams%PerVirt

 case ("EIGFCI")
       read(val, *) SystemParams%EigFCI

 case ("UATOMS")
       read(val, *) SystemParams%UCen
! maybe sth more fancy i.e. swapping monomers
      ! call get_ncen(val,SystemParams)
 case ("THRACT")
       read(val,*) SystemParams%ThrAct

 case ("THRSELACT")
       SystemParams%DeclareThrSelAct = .true.
       read(val,*) SystemParams%ThrSelAct

 case ("THRQVIRT")
       SystemParams%DeclareThrQVirt = .true.
       read(val,*) SystemParams%ThrQVirt

 case ("THRQINACT")
       SystemParams%DeclareThrQInact = .true.
       read(val,*) SystemParams%ThrQInact

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

 case ("WEXCIT")
       read(val,*) SystemParams%Wexcit

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

subroutine free_Input(Input)
implicit none
type(InputData) :: Input

deallocate(Input%SystemInput(1)%InSt)
deallocate(Input%SystemInput)

end subroutine free_Input

subroutine print_Input(Input)
!
! A highly imperfect subroutine for Input print
! should be replaced with sth smarter
!
implicit none

type(InputData) :: Input
integer :: i,imon

write(LOUT,'()')
write(LOUT,'(1x,a)') 'INPUT '
write(LOUT,'(8a10)') ('**********',i=1,8)

associate( CalcParams => Input%CalcParams)
 ! CALCULATION BLOCK
 if(allocated(CalcParams%JobTitle)) then
    write(LOUT,' (1x,a,4x,a)') "JOB TITLE: ", &
                 CalcParams%JobTitle
 else
    write(LOUT,'(1x,a,4x,a)') "JOB TITLE: ", &
                 "EMPTY"
 endif

 if(allocated(CalcParams%BasisSet)) &
    write(LOUT,'(1x,a,4x,a)') "BASIS SET: ", &
                 CalcParams%BasisSet
 if(allocated(CalcParams%BasisSet)) &
    write(LOUT,'(1x,a,4x,a)') "BASIS PATH:", &
                 CalcParams%BasisSetPath
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
 if (allocated(CalcParams%IntegralsFilePath)) then
       write(*, *) "Ints file: ", CalcParams%IntegralsFilePath
 end if

associate( CholeskyParams => Input%CholeskyParams)
 ! CHOLESKY BLOCK
 if(CholeskyParams%Cholesky>0) then
    write(LOUT, '()')
    write(LOUT,' (1x,a,5x,a)') "CHOLESKY :",  &
                       ".TRUE."
    if(CholeskyParams%CholeskyBIN>0) then
       write(LOUT,'(1x,a,5x,a)') "ALGORITHM:", &
                       "BINARY"
    endif
    if(CholeskyParams%CholeskyOTF>0) then
       write(LOUT,'(1x,a,5x,a)') "ALGORITHM:", &
                       "ON-THE-FLY"
    endif
    write(LOUT,' (1x,a,a)') "CHOLESKY ACCU: ", &
               PossibleCholAccu(CholeskyParams%CholeskyAccu)
 endif
end associate

 ! SYSTEM BLOCK(S)
 do imon=1,CalcParams%imon
    associate(System => Input%SystemInput(imon))
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
      if(System%DeclareThrSelAct) then
         write(LOUT, '(1x,a,e13.6)') "THRESHOLD SELECT ACTIVE : ", System%ThrSelAct
      endif
      if(System%DeclareThrQVirt) then
         write(LOUT, '(1x,a,e13.6)') "THRESHOLD QUASI-VIRTUAL : ", System%ThrQVirt
      endif
      if(System%DeclareThrQInact) then
         write(LOUT, '(1x,a,e13.6)') "THRESHOLD QUASI-INACTIVE : ", System%ThrQInact
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
     write(LOUT, '(1x,a,6x,i3)') "IFlFrag1 ", &
                 (Input%Flags%IFlFrag1)
     write(LOUT, '(1x,a,6x,i3)') "IFl12   ", &
                 (Input%Flags%IFl12)
     write(LOUT, '(1x,a,6x,i3)') "ISAPT   ", &
                 (Input%Flags%ISAPT)
 endif

write(LOUT,'()')

end subroutine print_Input



subroutine read_memsrt(val,MemVal,MemType)

     character(*), intent(in)  :: val
     integer, intent(out)      :: MemVal
     integer, intent(out)      :: MemType

     character(:), allocatable :: s1,s2

     call split(val,s1,s2)

     read(s1,*) MemVal

     if(trim(uppercase(s2))=='MB') then
        !MemVal = MemVal * 1024_8**2
        MemType = 2
     else if(trim(uppercase(s2))=='GB') then
        !MemVal = MemVal * 1024_8**3
        MemType = 3
     else
        write(lout,'(1x,a)') 'Error in declaration of MemSort!'
        stop
     endif

end subroutine read_memsrt

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

subroutine read_trstatearray(val,intrst,delim)

     character(*), intent(in) :: val
     character(1), intent(in) :: delim

     integer :: ii,k
     logical :: dot
     integer,allocatable :: intrst(:,:)
     integer :: instates
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

           if(instates.gt.1) then
              write(lout,*) 'ONLY SINGLE TRDM POSSIBLE!'
              stop
           endif

           ! assign states
           allocate(intrst(2,instates))
           instates = 0
           arrloop: do
                     k = index(w, delim)
                     instates = instates + 1
                     if(k == 0) then
                         !print*, 'last ', w
                         if(dot) then
                            call split(w,s1,s2,'.')
                            read(s1, *) intrst(1,instates)
                            read(s2, *) intrst(2,instates)
                         else
                            s1 = w
                            read(s1, *) intrst(1,instates)
                            intrst(2,instates) = 1
                         endif
                         !print*, '1 2',inst(1,instates),inst(2,instates)
                         exit arrloop
                     endif
                  enddo arrloop

     end if

end subroutine read_trstatearray

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



end module inputfill
