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

  ! load mainp.f Flags
  select case(Input%CalcParams%InterfaceType)
  case(INTER_TYPE_DAL)
     Flags%IDALTON = 1
     Flags%IAO     = 0
     Flags%INO     = 0
     Flags%NoSym   = 1

  case(INTER_TYPE_MOL)
     Flags%IDALTON = 0
!     Flags%IAO     = 0
!     Flags%INO     = 0
     Flags%NoSym   = Input%CalcParams%SymType 

  case(INTER_TYPE_OWN)
     Flags%IDALTON = 0
  ! ????
  end select

   if(Input%CalcParams%Restart) Flags%IRes = 1 
  
! CAREFUL
! HERE !! RDMTypes-related flags 

endif

end subroutine fill_Flags

end module systemdef
