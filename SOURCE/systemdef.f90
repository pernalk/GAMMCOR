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
     Flags%IFlAC  = 0
     Flags%IFlSnd = 0
  
  end select

 ! Inactive
  Flags%IFlCore = Input%CalcParams%Inactive
 
 ! EERPA
 if(Input%CalcParams%Fragments==1) Flags%IFlFrag = 1

 Flags%IFl12 = FLAG_DEBUG_FL12

endif

! PRINT_ALL_FLAGS // IPrint
if(Input%CalcParams%IPrint.gt.1) then
   call print_Flags(Flags)
endif

end subroutine fill_Flags

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

end subroutine print_Flags

end module systemdef
