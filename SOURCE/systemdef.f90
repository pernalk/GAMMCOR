module systemdef
use types

implicit none

contains

subroutine check_Calc(CalcParams)
! perform control and feed commons
!implicit double precision (a-h,o-z)
!implicit integer (i-n)
implicit none

type(CalculationBlock) :: CalcParams
integer :: i

write(LOUT,'()') 
!write(LOUT,'(1x,a)') 'RUN INPUT CHECK'
!write(LOUT,'(8a10)') ('**********',i=1,8)
!write(LOUT,'()') 

!check RDMType
if(CalcParams%RDMType==RDM_TYPE_DMRG) then
   !if(CalcParams%InterfaceType.ne.INTER_TYPE_OWN) then
   !   write(6,'(1x,a)') "ERROR! RDMType DMRG REQUIRES:"
   !   write(6,'(1x,a)') "Interface      OWN"
   !   stop
   !elseif(CalcParams%RDMSource.ne.INTER_TYPE_OWN) then
   !if(CalcParams%RDMSource.ne.INTER_TYPE_OWN) then
   !   write(6,'(1x,a)') "ERROR! RDMType DMRG REQUIRES:"
   !   write(6,'(1x,a)') "RDMSource      OWN"
   !   stop
   !endif

endif

! check JobType
select case(CalcParams%JobType)
!SAPT
case(5)
   
   if(CalcParams%InterfaceType.eq.INTER_TYPE_OWN) then
      write(6,'(1x,a)') "ERROR! SAPT REQUIRES:"
      write(6,'(1x,a)') "Interface    DALTON/MOLPRO"
      stop
   endif

   if(CalcParams%RDMSource.eq.INTER_TYPE_OWN) then
      write(6,'(1x,a)') "ERROR! SAPT REQUIRES:"
      write(6,'(1x,a)') "RDMSource    DALTON/MOLPRO"
      stop
   endif
end select

! check interface
select case(CalcParams%InterfaceType)
! Dalton
case(1)
   if(CalcParams%RDMSource==INTER_TYPE_OWN) then
     ! write(LOUT,'(1x,a)') 'ERROR! Dalton INTERFACE REQUIRES:' 
     ! write(LOUT,'(1x,a)') 'RDMSource     DALTON'
     ! stop
      write(LOUT,'(1x,a)') 'WARNING! Dalton INTERFACE USED WITH:' 
      write(LOUT,'(1x,a)') 'RDMSource     OWN'
   elseif(CalcParams%SymType==TYPE_SYM) then 
      write(LOUT,'(1x,a)') 'ERROR! Dalton INTERFACE REQUIRES:' 
      write(LOUT,'(1x,a)') 'Symmetry      NOSYM'
      stop
   endif
! Molpro
case(2)
   if(CalcParams%RDMSource.ne.INTER_TYPE_MOL) then
      print*, 'RDMSource',CalcParams%RDMSource,INTER_TYPE_MOL
      write(LOUT,'(1x,a)') 'ERROR! Molpro INTERFACE REQUIRES:' 
      write(LOUT,'(1x,a)') 'RDMSource     MOLPRO' 
      stop
    elseif(CalcParams%RDMType==RDM_TYPE_GVB.or.&
           CalcParams%RDMType==RDM_TYPE_DMRG) then
      write(LOUT,'(1x,a)') 'ERROR! INTERFACE Molpro REQUIRES :'
      write(LOUT,'(1x,a)') 'RDMType     CAS' 
      stop
   endif
end select

!! check Fragments (EERPA)
!if(CalcParams%Fragments==1) then
!  write(6,'(1x,a)') 'EERPA CALCULATIONS REQUESTED'
!  write(LOUT,'()') 
!  if(CalcParams%RDMType==RDM_TYPE_CAS.or.&
!     CalcParams%RDMType==RDM_TYPE_DMRG) then
!     write(6,'(1x,a)') 'ERROR! EERPA REQUIRES:'
!     write(6,'(1x,a)') 'RDMType     GVB or APSG'
!     stop 
!  endif
!endif

write(LOUT,'(1x,a)') 'INPUT CHECK PASSED'
write(LOUT,'(8a10)') ('**********',i=1,8)


end subroutine check_Calc

subroutine fill_Flags(Input,Flags)
implicit none

type(InputData) :: Input
type(FlagsData) :: Flags
integer :: i

! if block FLAGS present: CALCULATION block ignored
! default AC-GVB flags are overwritten by the user - good luck!  
if(Input%iflag==1) then

   Flags = Input%Flags
   write(LOUT,'()')
   write(LOUT,'(1x,a)') 'IMPORTANT! FLAGS BLOCK PRESENT: CALCULATION &
               &BLOCK IGNORED!'
   write(LOUT,'(1x,a)') 'INSTEAD, DEFAULT AC-GVB FLAGS OVERWRITTEN BY &
               &VALUES GIVEN IN FLAGS BLOCK!'
   write(LOUT,'()')
   write(LOUT,'(8a10)') ('**********',i=1,8)

else
  ! set JobType
   Flags%JobType = Input%CalcParams%JobType

  ! set TwoEl type
   Flags%ITwoEl = Input%CalcParams%TwoMoInt

  ! Interface 
  select case(Input%CalcParams%InterfaceType)
  case(INTER_TYPE_DAL)
     Flags%IDALTON = 1
     Flags%IAO     = 0
     Flags%INO     = 0
     Flags%NoSym   = 1
     Flags%IA = 1 
     
  case(INTER_TYPE_MOL)
     Flags%IDALTON = 0
     Flags%IAO     = 1
     Flags%INO     = 1
     !Flags%NoSym   = Input%CalcParams%SymType 
     Flags%IA = 1 

  case(INTER_TYPE_OWN)
     Flags%IDALTON = 0
     Flags%IA = 1 
 
  case(INTER_TYPE_ORCA)
     Flags%IDALTON = 1
     Flags%ICASSCF = 1
     Flags%IDMRG   = 1
     Flags%iORCA   = 1
     Flags%IAO     = 0
     Flags%INO     = 0
     Flags%NoSym   = 1
     Flags%IA = 1 
  ! ????
  end select

  if(Input%CalcParams%Restart) Flags%IRes = 1 
  
! RDMType  
  select case(Input%CalcParams%RDMType)
  case(RDM_TYPE_GVB)
     Flags%IGVB = 1
     Flags%ICASSCF = 0

  case(RDM_TYPE_APSG)
     FLags%IGVB = 0
     Flags%ICASSCF = 0

  case(RDM_TYPE_CAS)
     FLags%IGVB = 0
     Flags%ICASSCF = 1

  case(RDM_TYPE_DMRG)
     Flags%ICASSCF = 1
     Flags%IDMRG   = 1

  case(RDM_TYPE_HF)
     FLags%IGVB = 0
     Flags%ICASSCF = 1
     Flags%ISHF = 1
  case default
     write(LOUT,'(1x,a)') 'RDMType not declared! Assuming ICASSCF=1!'
     FLags%IGVB    = 0
     Flags%ICASSCF = 1
  end select

  ! excitations
  Flags%ISERPA = 0
  Flags%IAPSG  = 0

  ! JobType
  select case(Input%CalcParams%JobType)
  case(JOB_TYPE_AC)
     Flags%IFlAC  = 1
     Flags%IFlSnd = 0
     if(Input%CalcParams%DFApp==2) then
        if(Input%CalcParams%PostCAS) then
           Flags%IFunSR = 4
        else
           Flags%IFunSR = 2
        endif
        Flags%IFunSRKer = Input%CalcParams%Kernel 
     endif
!     if(Input%CalcParams%DFApp==2) Flags%IFunSRKer = 1

  case(JOB_TYPE_AC0)
    ! HERE WILL BE CHANGED TO:
    !Flags%IFlAC = 0
     Flags%IFlAC  = 1
     Flags%IFlSnd = 1
     if(Input%CalcParams%DFApp==2) then
        if(Input%CalcParams%PostCAS) then
           Flags%IFunSR = 4
        else
           Flags%IFunSR = 2
        endif
        Flags%IFunSRKer = Input%CalcParams%Kernel 
     endif
!     if(Input%CalcParams%DFApp==2) Flags%IFunSRKer = 1

  ! same as ERPA
  case(JOB_TYPE_AC1)
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0
     if(Input%CalcParams%DFApp==2) then
        if(Input%CalcParams%PostCAS) then
           Flags%IFunSR = 4
        else
           Flags%IFunSR = 2
        endif
        Flags%IFunSRKer = Input%CalcParams%Kernel 
     endif
!     if(Input%CalcParams%DFApp==2) Flags%IFunSRKer = 1

  case(JOB_TYPE_SAPT)
     Flags%ISAPT  = 1
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0

     ! Response for SAPT
     select case(Input%CalcParams%Response)
     case(RESP_ERPA)
        Flags%ISERPA = 0 
     case(RESP_DFT)
        Flags%ISERPA = 0 
        Flags%IFunSR = 1
        if(Input%CalcParams%DFApp>0) then
          Flags%IFunSR = Input%CalcParams%DFApp
        endif 
     case(RESP_APSG)
        Flags%ISERPA = 2
     end select

  case(JOB_TYPE_PDFT)
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0
     Flags%IFunSR = 3 

  case(JOB_TYPE_CASPIDFT)
     Flags%IFunSR = 5 

  case(JOB_TYPE_CASPIDFTOPT)
     Flags%IFunSR = 6

  end select

 ! Inactive
  Flags%IFlCore = Input%CalcParams%Core
 
 !Flags%IFl12 = FLAG_DEBUG_FL12

endif

! SELECT ELECTRONIC STATE
Flags%NoSt = Input%SystemInput(1)%NoSt

if(allocated(Input%CalcParams%JobTitle)) then
  Flags%JobTitle = Input%CalcParams%JobTitle
else
  Flags%JobTitle = 'EMPTY'
endif

! PRINT_ALL_FLAGS / IPrint
if(Input%CalcParams%IPrint.gt.0) then
   call print_Flags(Flags)
endif

end subroutine fill_Flags

subroutine create_System(Input,Flags,System,SAPT)
implicit none

type(InputData) :: Input
type(FlagsData) :: Flags
type(SystemBlock) :: System
type(SaptData) :: SAPT

! fill System
if(Flags%ISAPT.Eq.0) then

   if(Input%CalcParams%imon/=1) then
     write(LOUT,'(a)') 'ERROR! TOO MANY SYSTEM BLOCKS!'
     stop
   endif   
  
   System%NoSt = Input%SystemInput(1)%NoSt
   System%NStates = Input%SystemInput(1)%NStates
   allocate(System%InSt(2,System%NStates)) 
   if(Input%SystemInput(1)%DeclareSt) then
      System%InSt = Input%SystemInput(1)%InSt 
   else
      ! read 1st available state
      System%InSt(:,1) = -255
   endif
   ! only one TrSt alllowed
   allocate(System%InTrSt(2,1))
   if(Input%SystemInput(1)%DeclareTrSt) then
      System%InTrSt = Input%SystemInput(1)%InTrSt
   else
      ! read 1st available state
      System%InTrSt(:,1) = -255
   endif
!   Print*, 'InTrSt',System%InTrSt(1,1),System%InTrSt(2,1)

   System%ZNucl  = Input%SystemInput(1)%ZNucl
   System%Charge = Input%SystemInput(1)%Charge
   System%NBasis = Input%CalcParams%NBasis
   System%Multiplicity = Input%SystemInput(1)%Multiplicity
   System%Omega  = Input%SystemInput(1)%Omega
   System%EigFCI = Input%SystemInput(1)%EigFCI
   System%ThrAct = Input%SystemInput(1)%ThrAct
   System%ThrSelAct = Input%SystemInput(1)%ThrSelAct
   System%TwoMoInt = Input%SystemInput(1)%TwoMoInt
   System%IPrint = Input%CalcParams%IPrint  
  
   System%XELE = (System%ZNucl - System%Charge)/2.0d0
   System%NELE = (System%ZNucl - System%Charge)/2
  
   ! checkif Active Declared
   if(.not.Input%SystemInput(1)%NActFromRDM) then
      System%NAct=Input%SystemInput(1)%NAct 
      print*, 'NAct from input!',System%NAct
   endif   

  ! write(*,*) "ZNucl:",Input%SystemInput(1)%ZNucl
  ! write(*,*) "Sys-val:",System%XELE, System%NELE

! fill SAPT
elseif(Flags%ISAPT.Eq.1) then

 SAPT%InterfaceType = Input%CalcParams%InterfaceType
 SAPT%IPrint = Input%CalcParams%IPrint
 SAPT%SaptLevel = Input%CalcParams%SaptLevel
 SAPT%ic6 = Input%CalcParams%vdWCoef
 if(SAPT%InterfaceType==2) SAPT%HFCheck = .false.
 ! temporary RSH
 if(Flags%IFunSR<3) then
   SAPT%doRSH = .true. 
   SAPT%monA%doRSH = .true.
   SAPT%monB%doRSH = .true.
 endif

 associate( monA => SAPT%monA, &
            monB => SAPT%monB ) 
   select case(Input%SystemInput(1)%Monomer)
   case(1)

      monA%NoSt    = Input%SystemInput(1)%NoSt
      monA%NStates = Input%SystemInput(1)%NStates
      allocate(monA%InSt(2,System%NStates)) 
      if(Input%SystemInput(1)%DeclareSt) then
         monA%InSt = Input%SystemInput(1)%InSt 
      else
         ! assume 1.1 state
         monA%InSt(:,1) = -255
      endif
      monA%ZNucl   = Input%SystemInput(1)%ZNucl
      monA%Charge  = Input%SystemInput(1)%Charge
      monA%ACAlpha = Input%SystemInput(1)%ACAlpha
      monA%Omega   = Input%SystemInput(1)%Omega 
      monA%EigFCI  = Input%SystemInput(1)%EigFCI
      monA%NBasis  = Input%CalcParams%NBasis
      monA%ThrAct       = Input%SystemInput(1)%ThrAct
      monA%ThrSelAct    = Input%SystemInput(1)%ThrSelAct
      monA%DeclareTwoMo = Input%SystemInput(1)%DeclareTwoMo
      monA%TwoMoInt = Input%SystemInput(1)%TwoMoInt
      monA%PostCAS  = Input%SystemInput(1)%PostCAS
      monA%ISHF     = Input%SystemInput(1)%ISHF
      monA%Multiplicity = Input%SystemInput(1)%Multiplicity

      monA%NCen    = Input%SystemInput(1)%NCen
      monA%UCen    = Input%SystemInput(1)%UCen
      monA%Monomer = Input%SystemInput(1)%Monomer
      monA%IPrint  = Input%CalcParams%IPrint  
     
      monA%XELE = (SAPT%monA%ZNucl - SAPT%monA%Charge)/2.0d0 
      monA%NELE = (SAPT%monA%ZNucl - SAPT%monA%Charge)/2 

      if(.not.Input%SystemInput(1)%NActFromRDM) then
         monA%NAct=Input%SystemInput(1)%NAct
         monA%NActFromRDM = Input%SystemInput(1)%NActFromRDM
      endif

      monB%NoSt    = Input%SystemInput(2)%NoSt
      monB%NStates = Input%SystemInput(2)%NStates
      allocate(monB%InSt(2,System%NStates)) 
      if(Input%SystemInput(2)%DeclareSt) then
         monB%InSt = Input%SystemInput(2)%InSt 
      else
         ! assume 1.1 state
         monB%InSt(:,1) = -255
      endif
      monB%ZNucl  = Input%SystemInput(2)%ZNucl
      monB%Charge = Input%SystemInput(2)%Charge
      monB%ACAlpha = Input%SystemInput(2)%ACAlpha
      monB%Omega  = Input%SystemInput(2)%Omega 
      monB%EigFCI = Input%SystemInput(2)%EigFCI
      monB%NBasis = Input%CalcParams%NBasis
      monB%Multiplicity = Input%SystemInput(2)%Multiplicity
      monB%ThrAct = Input%SystemInput(2)%ThrAct
      monB%ThrSelAct = Input%SystemInput(2)%ThrSelAct
      monB%DeclareTwoMo = Input%SystemInput(2)%DeclareTwoMo
      monB%TwoMoInt = Input%SystemInput(2)%TwoMoInt
      monB%PostCAS = Input%SystemInput(2)%PostCAS
      monB%ISHF = Input%SystemInput(2)%ISHF
      monB%NCen = Input%SystemInput(2)%NCen
      monB%UCen = Input%SystemInput(2)%UCen
      monB%Monomer = Input%SystemInput(2)%Monomer
      monB%IPrint = Input%CalcParams%IPrint  
     
      monB%XELE = (monB%ZNucl - monB%Charge)/2.0d0 
      monB%NELE = (monB%ZNucl - monB%Charge)/2 

      if(.not.Input%SystemInput(2)%NActFromRDM) then
         monB%NAct=Input%SystemInput(2)%NAct
         monB%NActFromRDM = Input%SystemInput(2)%NActFromRDM
      endif

   ! write(LOUT,*) monA%ZNucl,'MONO(1)A,case1,dupaaa'
   ! write(LOUT,*) monB%ZNucl,'MONO(1)B,case1'
   
   case(2)
  
      monA%NoSt   = Input%SystemInput(2)%NoSt
      monA%NStates = Input%SystemInput(2)%NStates
      allocate(monA%InSt(2,System%NStates)) 
      if(Input%SystemInput(2)%DeclareSt) then
         monA%InSt = Input%SystemInput(2)%InSt 
      else
         ! assume 1.1 state
         monA%InSt(:,1) = -255
      endif
      monA%ZNucl  = Input%SystemInput(2)%ZNucl
      monA%Charge = Input%SystemInput(2)%Charge
      monA%ACAlpha = Input%SystemInput(2)%ACAlpha
      monA%Omega  = Input%SystemInput(2)%Omega 
      monA%EigFCI = Input%SystemInput(2)%EigFCI
      monA%NBasis = Input%CalcParams%NBasis
      monA%Multiplicity = Input%SystemInput(2)%Multiplicity
      monA%ThrAct = Input%SystemInput(2)%ThrAct
      monA%ThrSelAct = Input%SystemInput(2)%ThrSelAct
      monA%DeclareTwoMo = Input%SystemInput(2)%DeclareTwoMo
      monA%TwoMoInt = Input%SystemInput(2)%TwoMoInt
      monA%PostCAS = Input%SystemInput(2)%PostCAS
      monA%ISHF = Input%SystemInput(2)%ISHF
      monA%NCen = Input%SystemInput(2)%NCen
      monA%UCen = Input%SystemInput(2)%UCen
      monA%Monomer = Input%SystemInput(2)%Monomer
      monA%IPrint = Input%CalcParams%IPrint  
     
      monA%XELE = (monA%ZNucl - monA%Charge)/2.0d0 
      monA%NELE = (monA%ZNucl - monA%Charge)/2 

      if(.not.Input%SystemInput(2)%NActFromRDM) then
         monA%NAct=Input%SystemInput(2)%NAct
         monA%NActFromRDM = Input%SystemInput(2)%NActFromRDM
      endif

      monB%NoSt   = Input%SystemInput(1)%NoSt
      monB%NStates = Input%SystemInput(1)%NStates
      allocate(monB%InSt(2,System%NStates)) 
      if(Input%SystemInput(1)%DeclareSt) then
         monB%InSt = Input%SystemInput(1)%InSt 
      else
         ! assume 1.1 state
         monB%InSt(:,1) = -255
      endif
      monB%ZNucl  = Input%SystemInput(1)%ZNucl
      monB%Charge = Input%SystemInput(1)%Charge
      monB%ACAlpha = Input%SystemInput(1)%ACAlpha
      monB%Omega  = Input%SystemInput(1)%Omega
      monB%EigFCI = Input%SystemInput(1)%EigFCI
      monB%NBasis = Input%CalcParams%NBasis
      monB%Multiplicity = Input%SystemInput(1)%Multiplicity
      monB%ThrAct = Input%SystemInput(1)%ThrAct
      monB%ThrSelAct = Input%SystemInput(1)%ThrSelAct
      monB%DeclareTwoMo = Input%SystemInput(1)%DeclareTwoMo
      monB%TwoMoInt = Input%SystemInput(1)%TwoMoInt
      monB%PostCAS = Input%SystemInput(1)%PostCAS
      monB%ISHF = Input%SystemInput(1)%ISHF
      monB%NCen = Input%SystemInput(1)%NCen
      monB%UCen = Input%SystemInput(1)%UCen
      monB%Monomer = Input%SystemInput(1)%Monomer
      monB%IPrint = Input%CalcParams%IPrint  
     
      monB%XELE = (monB%ZNucl - monB%Charge)/2.0d0 
      monB%NELE = (monB%ZNucl - monB%Charge)/2 
 
      if(.not.Input%SystemInput(1)%NActFromRDM) then
         monB%NActFromRDM = Input%SystemInput(1)%NActFromRDM
         monB%NAct=Input%SystemInput(1)%NAct
      endif

   end select
 end associate

 ! set in/out-of-core
 ! SAPT-CAS: unless defined in input, set FOFO as default!
 if(Flags%ICASSCF==1.and.Flags%ISERPA==0) then
    if(.not.SAPT%monA%DeclareTwoMo) SAPT%monA%TwoMoInt = TWOMO_FOFO 
    if(.not.SAPT%monB%DeclareTwoMo) SAPT%monB%TwoMoInt = TWOMO_FOFO
 endif

 ! set SameOm for SAPT-RSH
 if(SAPT%doRSH) then
    if(SAPT%monA%Omega/=SAPT%monB%Omega) then
      SAPT%SameOm = .false.
      SAPT%monA%SameOm = .false.
      SAPT%monB%SameOm = .false.
    endif
 endif

System = SAPT%monB

endif

call print_System(Flags,System,SAPT)

end subroutine create_System

subroutine print_System(Flags,System,SAPT)
implicit none

type(FlagsData) :: Flags
type(SystemBlock) :: System
type(SaptData) :: SAPT
integer :: i

write(LOUT,'()')
write(LOUT,'(1x,a)') 'SYSTEM '
write(LOUT,'(8a10)') ('**********',i=1,8)

if(Flags%ISAPT.Eq.0) then
   write(LOUT,'(1x,a,1x,i2)') 'NUCLEAR CHARGE: ', System%ZNucl
   write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', System%Charge

elseif(Flags%ISAPT.Eq.1) then
    write(LOUT,'(1x,a)') 'MONOMER A'
    write(LOUT,'(1x,a,1x,i2)') 'NUCLEAR CHARGE: ', SAPT%monA%ZNucl
    write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', SAPT%monA%Charge
    write(LOUT,'(1x,a,3x,i3)') 'NO.OF ATOMS: ', SAPT%monA%NCen
    write(LOUT,'(1x,a,2x,f9.4)') 'THRESH ACTIVE: ', SAPT%monA%ThrAct
    write(LOUT,'()')
    write(LOUT,'(1x,a)') 'MONOMER B'
    write(LOUT,'(1x,a,1x,i2)') 'NUCLEAR CHARGE: ', SAPT%monB%ZNucl
    write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', SAPT%monB%Charge
    write(LOUT,'(1x,a,3x,i3)') 'NO.OF ATOMS: ', SAPT%monB%NCen
    write(LOUT,'(1x,a,2x,f9.4)') 'THRESH ACTIVE: ', SAPT%monB%ThrAct

endif

!   write(LOUT,'()')
!   write(LOUT,'(1x,a,1x,i3)') 'NO. OF CONTRACTIONS: ', System%NBasis

   if(Flags%IRes.Ne.1) then
!      write(LOUT,'()')
!      write(LOUT,'(1x,a)') 'DENSITY MATRIX FUNCTIONAL CALCULATION'
!      write(LOUT,'(1x,a,1x,i3)') 'FUNCTIONAL', Flags%IFun
   elseif(Flags%IRes.Eq.1) then
      write(LOUT,'(1x,a)') 'RESTART REQUESTED' 
   endif   

if(Flags%ISAPT.Eq.0) then
!   write(LOUT,'()')
   write(LOUT,'(1x,a,1x,i5)') 'PRINT LEVEL: ', System%IPrint

elseif(Flags%ISAPT.Eq.1) then
!   write(LOUT,'()')
   write(LOUT,'(1x,a,1x,i5)') 'PRINT LEVEL: ', SAPT%monA%IPrint
 
endif

   write(LOUT,'()')

end subroutine print_System

subroutine print_Flags(Flags)
implicit none
type(FlagsData) :: Flags

write(LOUT, '()')
write(LOUT, '(1x,a,6x,i3)') "IDALTON ", &
             (Flags%IDALTON)
write(LOUT, '(1x,a,6x,i3)') "IRes    ", &
             (Flags%IRes)
write(LOUT, '(1x,a,6x,i3)') "IAO     ", &
             (Flags%IAO)
write(LOUT, '(1x,a,6x,i3)') "INO     ", &
             (Flags%INO)
write(LOUT, '(1x,a,6x,i3)') "NoSym   ", &
             (Flags%NoSym)
write(LOUT, '(1x,a,6x,i3)') "IGVB    ", &
             (Flags%IGVB)
write(LOUT, '(1x,a,6x,i3)') "IFun    ", &
             (Flags%IFun)
write(LOUT, '(1x,a,6x,i3)') "IFunSR  ", &
             (Flags%IFunSR)
write(LOUT, '(1x,a,5x,i3)') "IFunSRKer", &
             (Flags%IFunSRKer)
write(LOUT, '(1x,a,6x,i3)') "IModG   ", &
             (Flags%IModG)
write(LOUT, '(1x,a,6x,i3)') "NGOcc   ", &
             (Flags%NGOcc)
write(LOUT, '(1x,a,6x,i3)') "ILoc    ", &
             (Flags%ILoc)
write(LOUT, '(1x,a,6x,i3)') "IFreeze ", &
             (Flags%IFreeze)
write(LOUT, '(1x,a,6x,i3)') "IAPSG   ", &
             (Flags%IAPSG)
write(LOUT, '(1x,a,6x,i3)') "ISERPA  ", &
             (Flags%ISERPA)
write(LOUT, '(1x,a,6x,i3)') "IA      ", &
             (Flags%IA)
write(LOUT, '(1x,a,6x,i3)') "ICASSCF ", &
             (Flags%ICASSCF)
write(LOUT, '(1x,a,6x,i3)') "IDMRG   ", &
             (Flags%IDMRG)
write(LOUT, '(1x,a,6x,i3)') "IFlAC   ", &
             (Flags%IFlAC)
write(LOUT, '(1x,a,6x,i3)') "IFLSnd  ", &
             (Flags%IFlSnd)
write(LOUT, '(1x,a,6x,i3)') "IFlCore ", &
             (Flags%IFlCore)
write(LOUT, '(1x,a,6x,i3)') "IFlFrag ", &
             (Flags%IFlFrag1)
write(LOUT, '(1x,a,6x,i3)') "IFl12   ", &
              (Flags%IFl12)
write(LOUT, '(1x,a,6x,i3)') "ISAPT   ", &
              (Flags%ISAPT)

end subroutine print_Flags

end module systemdef
