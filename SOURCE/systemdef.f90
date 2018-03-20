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
   if(CalcParams%InterfaceType.ne.INTER_TYPE_OWN) then
      write(6,'(1x,a)') "ERROR! RDMType DMRG REQUIRES:"
      write(6,'(1x,a)') "Interface      OWN"
      stop
   elseif(CalcParams%RDMSource.ne.INTER_TYPE_OWN) then
      write(6,'(1x,a)') "ERROR! RDMType DMRG REQUIRES:"
      write(6,'(1x,a)') "RDMSource      OWN"
      stop
   endif

endif

! check JobType
select case(CalcParams%JobType)
!SAPT
case(5)
   
   if(CalcParams%InterfaceType.ne.INTER_TYPE_DAL) then
      write(6,'(1x,a)') "ERROR! SAPT REQUIRES:"
      write(6,'(1x,a)') "Interface      DALTON"
      stop
   endif

   if(CalcParams%RDMSource.ne.INTER_TYPE_DAL) then
      write(6,'(1x,a)') "ERROR! SAPT REQUIRES:"
      write(6,'(1x,a)') "RDMSource      DALTON"
      stop
   endif
end select

! check interface
select case(CalcParams%InterfaceType)
! Dalton
case(1)
   if(CalcParams%RDMSource==INTER_TYPE_OWN) then
      write(LOUT,'(1x,a)') 'ERROR! Dalton INTERFACE REQUIRES:' 
      write(LOUT,'(1x,a)') 'RDMSource     DALTON'
      stop
   elseif(CalcParams%SymType==TYPE_SYM) then 
      write(LOUT,'(1x,a)') 'ERROR! Dalton INTERFACE REQUIRES:' 
      write(LOUT,'(1x,a)') 'Symmetry      NOSYM'
      stop
   endif
! Molpro
case(2)
   if(CalcParams%RDMSource.ne.INTER_TYPE_OWN) then
      write(LOUT,'(1x,a)') 'ERROR! Molpro INTERFACE REQUIRES:' 
      write(LOUT,'(1x,a)') 'RDMSource     OWN' 
      stop
    elseif(CalcParams%RDMType==RDM_TYPE_CAS.or.&
           CalcParams%RDMType==RDM_TYPE_DMRG) then
      write(LOUT,'(1x,a)') 'ERROR! INTERFACE Molpro REQUIRES :'
      write(LOUT,'(1x,a)') 'RDMType     GVB or APSG' 
      stop
   endif
! OWN ??? 
end select

! check Fragments (EERPA)
if(CalcParams%Fragments==1) then
  write(6,'(1x,a)') 'EERPA CALCULATIONS REQUESTED'
  write(LOUT,'()') 
  if(CalcParams%RDMType==RDM_TYPE_CAS.or.&
     CalcParams%RDMType==RDM_TYPE_DMRG) then
     write(6,'(1x,a)') 'ERROR! EERPA REQUIRES:'
     write(6,'(1x,a)') 'RDMType     GVB or APSG'
     stop 
  endif
endif

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
!    Flags%IAO     = 0
!    Flags%INO     = 0
     Flags%NoSym   = Input%CalcParams%SymType 
     Flags%IA = 1 

  case(INTER_TYPE_OWN)
     Flags%IDALTON = 0
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
  end select

! excitations
  Flags%ISERPA = 0
  Flags%IAPSG  = 0

! JobType
  select case(Input%CalcParams%JobType)
  case(JOB_TYPE_AC)
     Flags%IFlAC  = 1
     Flags%IFlSnd = 0   

  case(JOB_TYPE_AC0)
! HERE WILL BE CHANGED TO:
    !Flags%IFlAC = 0
     Flags%IFlAC  = 1
     Flags%IFlSnd = 1

! same as ERPA  
  case(JOB_TYPE_AC1)
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0
 
  case(JOB_TYPE_SAPT)
     Flags%ISAPT  = 1
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0
  
  end select

 ! Inactive
  Flags%IFlCore = Input%CalcParams%Inactive
 
 ! EERPA
 if(Input%CalcParams%Fragments==1) Flags%IFlFrag = 1

 Flags%IFl12 = FLAG_DEBUG_FL12

endif

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
  
   System%ZNucl  = Input%SystemInput(1)%ZNucl
   System%Charge = Input%SystemInput(1)%Charge
   System%NBasis = Input%CalcParams%NBasis
   System%Multiplicity = Input%SystemInput(1)%Multiplicity
   System%IPrint = Input%CalcParams%IPrint  
  
   System%XELE = (System%ZNucl - System%Charge)/2.0d0 
   System%NELE = (System%ZNucl - System%Charge)/2 
   
  ! write(*,*) "ZNucl:",Input%SystemInput(1)%ZNucl
  ! write(*,*) "Sys-val:",System%XELE, System%NELE

! fill SAPT
elseif(Flags%ISAPT.Eq.1) then

 SAPT%IPrint = Input%CalcParams%IPrint

 associate( monA => SAPT%monA, &
            monB => SAPT%monB ) 
   select case(Input%SystemInput(1)%Monomer)
   case(1)
      monA%ZNucl  = Input%SystemInput(1)%ZNucl
      monA%Charge = Input%SystemInput(1)%Charge
      monA%NBasis = Input%CalcParams%NBasis
      monA%Multiplicity = Input%SystemInput(1)%Multiplicity
      monA%NCen = Input%SystemInput(1)%NCen
      monA%Monomer = Input%SystemInput(1)%Monomer
      monA%IPrint = Input%CalcParams%IPrint  
     
      monA%XELE = (SAPT%monA%ZNucl - SAPT%monA%Charge)/2.0d0 
      monA%NELE = (SAPT%monA%ZNucl - SAPT%monA%Charge)/2 
    
      monB%ZNucl  = Input%SystemInput(2)%ZNucl
      monB%Charge = Input%SystemInput(2)%Charge
      monB%NBasis = Input%CalcParams%NBasis
      monB%Multiplicity = Input%SystemInput(2)%Multiplicity
      monB%NCen = Input%SystemInput(2)%NCen
      monB%Monomer = Input%SystemInput(2)%Monomer
      monB%IPrint = Input%CalcParams%IPrint  
     
      monB%XELE = (monB%ZNucl - monB%Charge)/2.0d0 
      monB%NELE = (monB%ZNucl - monB%Charge)/2 
    
   ! write(LOUT,*) monA%ZNucl,'MONO(1)A,case1,dupaaa'
   ! write(LOUT,*) monB%ZNucl,'MONO(1)B,case1'
   
   case(2)
  
      monA%ZNucl  = Input%SystemInput(2)%ZNucl
      monA%Charge = Input%SystemInput(2)%Charge
      monA%NBasis = Input%CalcParams%NBasis
      monA%Multiplicity = Input%SystemInput(2)%Multiplicity
      monA%NCen = Input%SystemInput(2)%NCen
      monA%Monomer = Input%SystemInput(2)%Monomer
      monA%IPrint = Input%CalcParams%IPrint  
     
      monA%XELE = (monA%ZNucl - monA%Charge)/2.0d0 
      monA%NELE = (monA%ZNucl - monA%Charge)/2 
   
      monB%ZNucl  = Input%SystemInput(1)%ZNucl
      monB%Charge = Input%SystemInput(1)%Charge
      monB%NBasis = Input%CalcParams%NBasis
      monB%Multiplicity = Input%SystemInput(1)%Multiplicity
      monB%NCen = Input%SystemInput(1)%NCen
      monB%Monomer = Input%SystemInput(1)%Monomer
      monB%IPrint = Input%CalcParams%IPrint  
     
      monB%XELE = (monB%ZNucl - monB%Charge)/2.0d0 
      monB%NELE = (monB%ZNucl - monB%Charge)/2 
   
  !  write(LOUT,*) monA%ZNucl,'MONO(1)A,case2'
  !  write(LOUT,*) monB%ZNucl,'MONO(2)B,case2,dupaa'
  
   end select
 end associate

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
   write(LOUT,'(1x,a,1x,i3)') 'NUCLEAR CHARGE: ', System%ZNucl
   write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', System%Charge

elseif(Flags%ISAPT.Eq.1) then
    write(LOUT,'(1x,a)') 'MONOMER A'
    write(LOUT,'(1x,a,1x,i3)') 'NUCLEAR CHARGE: ', SAPT%monA%ZNucl
    write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', SAPT%monA%Charge
    write(LOUT,'(1x,a,3x,i3)') 'NO.OF ATOMS: ', SAPT%monA%NCen
    write(LOUT,'()')
    write(LOUT,'(1x,a)') 'MONOMER B'
    write(LOUT,'(1x,a,1x,i3)') 'NUCLEAR CHARGE: ', SAPT%monB%ZNucl
    write(LOUT,'(1x,a,8x,i3)') 'CHARGE: ', SAPT%monB%Charge
    write(LOUT,'(1x,a,3x,i3)') 'NO.OF ATOMS: ', SAPT%monB%NCen

endif

   write(LOUT,'()')
   write(LOUT,'(1x,a,1x,i3)') 'NO. OF CONTRACTIONS: ', System%NBasis

   if(Flags%IRes.Ne.1) then
!      write(LOUT,'()')
!      write(LOUT,'(1x,a)') 'DENSITY MATRIX FUNCTIONAL CALCULATION'
!      write(LOUT,'(1x,a,1x,i3)') 'FUNCTIONAL', Flags%IFun
   elseif(Flags%IRes.Eq.1) then
      write(LOUT,'(1x,a)') 'RESTART REQUESTED' 
   endif   

if(Flags%ISAPT.Eq.0) then
!   write(LOUT,'()')
   write(LOUT,'(1x,a,1x,i3)') 'PRINT LEVEL: ', System%IPrint

elseif(Flags%ISAPT.Eq.1) then
!   write(LOUT,'()')
   write(LOUT,'(1x,a,1x,i3)') 'PRINT LEVEL: ', SAPT%monA%IPrint
 
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
             (Flags%IFlFrag)
write(LOUT, '(1x,a,6x,i3)') "IFl12   ", &
              (Flags%IFl12)
write(LOUT, '(1x,a,6x,i3)') "ISAPT   ", &
              (Flags%ISAPT)

end subroutine print_Flags

end module systemdef
