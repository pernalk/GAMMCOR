module sapt_main

use types
use timing
use sapt_inter
use sapt_resp
use sapt_pol
use sapt_Chol_pol
use sapt_Chol_exch
use sapt_exch
use exd_pino
use omp_lib

implicit none

!save
!double precision :: Tcpu,Twall

contains

subroutine sapt_driver(Flags,SAPT)
implicit none

type(FlagsData) :: Flags
type(SaptData)  :: SAPT
integer :: i
integer :: NBasis
double precision :: Tcpu,Twall

 Flags%SaptLevel = SAPT%SaptLevel

! ERPA
 Flags%IFlAC  = 0
 Flags%IFlSnd = 0

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'STARTING SAPT CALCULATIONS'
 write(LOUT,'(8a10)') ('**********',i=1,8)

 ! jump to reduceVirt framework
 if(Flags%IRedVirt==1) call sapt_driver_red(Flags,SAPT)

 call clock('START',Tcpu,Twall)
 call sapt_basinfo(SAPT,NBasis)
 call sapt_interface(Flags,SAPT,NBasis)

 call sapt_mon_ints(SAPT%monA,Flags,NBasis)
 call sapt_mon_ints(SAPT%monB,Flags,NBasis)

 call sapt_response(Flags,SAPT%monA,SAPT%EnChck,NBasis)
 call sapt_response(Flags,SAPT%monB,SAPT%EnChck,NBasis)

 call sapt_ab_ints(Flags,SAPT%monA,SAPT%monB,SAPT%iPINO,NBasis)

 ! SAPT components
 write(LOUT,'()')

 ! switch to Cholesky SAPT
 if(Flags%ICholesky==1) call sapt_Cholesky(Flags,SAPT,Tcpu,TWall,NBasis)

 ! switch to extrapolated SAPT
 if(SAPT%monA%Cubic.or.SAPT%monB%Cubic) call sapt_extrapol(Flags,SAPT,NBasis)

 if(Flags%ISERPA==0.and.SAPT%ic6==1) then

    ! calculate only C6
    call c6_dummy(Flags,SAPT%monA,SAPT%monB,SAPT)

 elseif(Flags%ISERPA==0.and.SAPT%ic6==0) then

    call e1elst(SAPT%monA,SAPT%monB,SAPT)
    call e1exchs2(Flags,SAPT%monA,SAPT%monB,SAPT)
    !call e1exch_NaNb(Flags,SAPT%monA,SAPT%monB,SAPT)
    !if(SAPT%monA%NELE==1.and.SAPT%monB%NELE==1) call hl_2el(Flags,SAPT%monA,SAPT%monB,SAPT)

    if(SAPT%SaptLevel==0) then
       SAPT%iCpld = .false.
       call e2ind_unc(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2disp_unc(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

    elseif(SAPT%SaptLevel==10) then
       SAPT%SemiCoupled = .false.
       call e2dispCAS(Flags,SAPT%monA,SAPT%monB,SAPT,NBasis)

    elseif(SAPT%SaptLevel==2) then
       call e2ind(Flags,SAPT%monA,SAPT%monB,SAPT)
      !call e2ind_resp(Flags,SAPT%monA,SAPT%monB,SAPT) ! version using alpha(0)
      ! call e2ind_hf_icphf(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2disp(Flags,SAPT%monA,SAPT%monB,SAPT)
      ! call e2disp_Cmat(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)
    endif

 elseif(Flags%ISERPA==2) then

    call e1elst(SAPT%monA,SAPT%monB,SAPT)
    !call e1exchs2_pino(Flags,SAPT%monA,SAPT%monB,SAPT)

    if(Flags%ISHF==0) then
       ! APSG procedures use eigenvectors from APSG_NEST
       call e2ind_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2disp_apsg(Flags,SAPT%monA,SAPT%monB,SAPT)
    elseif(Flags%ISHF==1) then
       ! PINO procedures use eigenvectors from OptTwoP
       ! i.e. AP (eigenvec) and PP (eigenval)
       call e2ind_pino(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2disp_pino(Flags,SAPT%monA,SAPT%monB,SAPT)
       call e2exdisp_pino(Flags,SAPT%monA,SAPT%monB,SAPT)
    endif

 endif

 call summary_sapt(SAPT)
 !call summary_sapt_verbose(SAPT)

 call print_warn(SAPT)
 call free_sapt(Flags,SAPT)

 call clock('SAPT',Tcpu,Twall)

 stop

end subroutine sapt_driver

subroutine sapt_driver_red(Flags,SAPT)
! sapt driver with reduced virt space
! Flags%IRedVirt==1
implicit none

type(FlagsData)  :: Flags
type(SaptData)   :: SAPT
integer          :: i
integer          :: NBasis,NBasisRed
integer          :: SAPT_LEVEL_SAVE
double precision :: e2d,e2d_unc,e2dR,e2dR_unc
double precision :: e2exd,e2exd_unc,e2exdR,e2exdR_unc
double precision :: coefDisp,coefExDisp,coefGen
double precision :: Tcpu,Twall
logical          :: onlyDisp

 ! test only?
 onlyDisp = .false.
 !onlyDisp = .true.

 SAPT_LEVEL_SAVE  = SAPT%SaptLevel
 SAPT%SemiCoupled = .false.

 ! 1) uncoupled response
 SAPT%SaptLevel  = 0
 Flags%SaptLevel = 0
 SAPT%iCpld      = .false.
 Flags%IFlag0    = 1

 call clock('START',Tcpu,Twall)
 call sapt_basinfo(SAPT,NBasis)
 call sapt_interface(Flags,SAPT,NBasis)

 call sapt_mon_ints(SAPT%monA,Flags,NBasis)
 call sapt_mon_ints(SAPT%monB,Flags,NBasis)

 call sapt_response(Flags,SAPT%monA,SAPT%EnChck,NBasis)
 call sapt_response(Flags,SAPT%monB,SAPT%EnChck,NBasis)

 call sapt_ab_ints(Flags,SAPT%monA,SAPT%monB,SAPT%iPINO,NBasis)

 write(LOUT,'(/,1x,a)') 'SAPT COMPONENTS'
 write(LOUT,'(8a10)') ('**********',i=1,6)

 call e1elst(SAPT%monA,SAPT%monB,SAPT)
 !call e1exchs2(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2disp_unc(Flags,SAPT%monA,SAPT%monB,SAPT)
 if(onlyDisp.eqv..false.) call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

 e2d_unc   = SAPT%e2disp_unc*1000d0
 if(onlyDisp.eqv..false.) e2exd_unc = SAPT%e2exdisp_unc*1000d0

 call clock('SAPT(FULL SPACE)',Tcpu,Twall)

 call reduce_virt(Flags,SAPT%monA,NBasis)
 call reduce_virt(Flags,SAPT%monB,NBasis)

 ! new
 NBasisRed = NBasis - min(SAPT%monA%NVZero,SAPT%monB%NVZero)

 if(SAPT_LEVEL_SAVE/=0) then
    SAPT%SaptLevel  = 2
    Flags%SaptLevel = 2
    SAPT%iCpld      = .true.
    !Flags%IFlag0    = 1
 endif

 call sapt_response(Flags,SAPT%monA,SAPT%EnChck,NBasis)
 call sapt_response(Flags,SAPT%monB,SAPT%EnChck,NBasis)

 call sapt_ab_ints_red(Flags,SAPT%monA,SAPT%monB,SAPT%iPINO,NBasis,NBasisRed)

 call e2disp_unc(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2disp_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 if(onlyDisp.eqv..false.) call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

 e2dR       = SAPT%e2disp*1000d0
 e2dR_unc   = SAPT%e2disp_unc*1000d0
 if(onlyDisp.eqv..false.) e2exdR     = SAPT%e2exdisp*1000d0
 if(onlyDisp.eqv..false.) e2exdR_unc = SAPT%e2exdisp_unc*1000d0

 write(LOUT,'(/,1x,a)') 'SAPT SCALED SUMMARY'
 write(LOUT,'(8a10)') ('**********',i=1,6)

 write(LOUT,'(1x,a,i3)') 'SAPT level  =', SAPT%SaptLevel
 write(LOUT,'(1x,a,f16.8,f16.8)') 'E2disp(unc) = ',e2d_unc,e2dR_unc
 write(LOUT,'(1x,a,f16.8)')       'E2dispR     = ',e2dR
 if(onlyDisp.eqv..false.) write(LOUT,'(1x,a,f16.8,f16.8)') 'E2exd(unc)  = ',e2exd_unc,e2exdR_unc
 if(onlyDisp.eqv..false.) write(LOUT,'(1x,a,f16.8)')       'E2exdpR     = ',e2exdR

 coefGen    = e2dR/e2dR_unc
 coefDisp   = e2d_unc/e2dR_unc
 if(onlyDisp.eqv..false.) coefExDisp = e2exd_unc/e2exdR_unc
 write(LOUT,'(/,1x,a,f16.8)') 'coefDisp   = ',coefDisp
 if(onlyDisp.eqv..false.) write(LOUT,'(1x,a,f16.8)')   'coefExDisp = ',coefExDisp
 write(LOUT,'(/,1x,a,f16.8)') 'E2disp(SCALED)   = ',e2dR*coefDisp
 if(onlyDisp.eqv..false.) write(LOUT,'(1x,a,f16.8)')   'E2exdisp(SCALED) = ',e2exdR*coefExDisp
 if(onlyDisp.eqv..false.) write(LOUT,'(1x,a,f16.8)')   'E2exdisp(SCALE2) = ',e2exd_unc*coefGen

 call print_warn(SAPT)
 call free_sapt(Flags,SAPT)

 call clock('SAPT',Tcpu,Twall)

 stop

end subroutine sapt_driver_red

subroutine sapt_Cholesky(Flags,SAPT,Tcpu,Twall,NBasis)
!
! sapt driver with Cholesky decompositon
! Flags%isCholesky==1
!
implicit none

type(FlagsData)    :: Flags
type(SaptData)     :: SAPT
integer,intent(in) :: NBasis
double precision,intent(inout) :: Tcpu,Twall
integer :: i

 write(LOUT,'(1x,a)') 'SAPT(MC) with Cholesky decomposition'
 write(LOUT,'(8a10)') ('----------',i=1,8)
 if(SAPT%CAlpha) print*, 'SAPT%CAlpha',SAPT%CAlpha

 if(SAPT%SaptLevel==999 .or. SAPT%SaptLevel==666) then
    if(SAPT%SaptLevel==999)  write(LOUT,'(1x,a,/)') 'RSPT2 calculation requested'
    if(SAPT%SaptLevel==666)  write(LOUT,'(1x,a,/)') 'RSPT2+ calculation requested'

    call e1elst_Chol(SAPT%monA,SAPT%monB,SAPT)
    if(SAPT%SaptLevel==666) call e1exch_Chol(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e2ind_icerpa(Flags,SAPT%monA,SAPT%monB,SAPT)
    if(.not.SAPT%CAlpha) then
      call e2disp_Cmat_Chol_block(Flags,SAPT%monA,SAPT%monB,SAPT)
    else if(SAPT%CAlpha) then
      call e2disp_CAlphaTilde_block(Flags,SAPT%monA,SAPT%monB,SAPT)
    endif

    call summary_rspt(SAPT)

 else

    call e1elst_Chol(SAPT%monA,SAPT%monB,SAPT)
    !call e1exchs2(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e1exch_NaNb(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e2ind(Flags,SAPT%monA,SAPT%monB,SAPT)
    call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
   
    ! test subroutines for E2disp
    !call e2disp_Chol(Flags,SAPT%monA,SAPT%monB,SAPT)
    !call e2disp_Cmat(Flags,SAPT%monA,SAPT%monB,SAPT)

    if(.not.SAPT%CAlpha) then
      ! Adam's test
      !call e2disp_Cmat_Chol(Flags,SAPT%monA,SAPT%monB,SAPT)
      !call e2disp_Cmat_Chol_diag(Flags,SAPT%monA,SAPT%monB,SAPT)
      call e2disp_Cmat_Chol_block(Flags,SAPT%monA,SAPT%monB,SAPT)
      !call e2disp_Cmat_Chol_proj(Flags,SAPT%monA,SAPT%monB,SAPT)
    else if(SAPT%CAlpha) then
      call e2disp_CAlphaTilde_block(Flags,SAPT%monA,SAPT%monB,SAPT)
      !call e2disp_CAlphaTilde_full(Flags,SAPT%monA,SAPT%monB,SAPT)
    endif
    call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)
    call summary_sapt(SAPT)

 endif

 call print_warn(SAPT)
 call free_sapt(Flags,SAPT)

 call clock('SAPT',Tcpu,Twall)

 stop

end subroutine sapt_Cholesky

subroutine sapt_extrapol(Flags,SAPT,NBasis)
implicit none

type(SaptData)     :: SAPT
type(FlagsData)    :: Flags
integer,intent(in) :: NBasis

integer            :: i

 ! test
 if(.not.(SAPT%monA%Cubic).and..not.(SAPT%monB%Cubic)) then
   write(lout,'(1x,a)') 'ERROR! Wrong call of sapt_extrapol!'
   stop
 endif

 SAPT%SemiCoupled = .false.

 ! first order
 call e1elst(SAPT%monA,SAPT%monB,SAPT)
 call e1exchs2(Flags,SAPT%monA,SAPT%monB,SAPT)

 ! E2(cubic) = a*alpha^3 + b*alpha + c

 write(lout,'(/1x,a)') 'SAPT, E(2) cubic'
 write(LOUT,'(8a10)') ('----------',i=1,4)
 ! second order : uncoupled (for c coefficient)
 if(SAPT%monA%Cubic) SAPT%monA%ACAlpha=SAPT%monA%ACAlpha0
 if(SAPT%monB%Cubic) SAPT%monB%ACAlpha=SAPT%monB%ACAlpha0
 write(lout,'(1x,a,f9.6)') 'ACAlpha(A) =',SAPT%monA%ACAlpha
 write(lout,'(1x,a,f9.6)') 'ACAlpha(B) =',SAPT%monB%ACAlpha

 call e2ind_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2disp_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

 !              : semicoupled (for b coefficient)
 if(SAPT%monA%Cubic) SAPT%monA%ACAlpha=SAPT%monA%ACAlpha1
 if(SAPT%monB%Cubic) SAPT%monB%ACAlpha=SAPT%monB%ACAlpha1
 write(lout,'(/1x,a,f9.6)') 'ACAlpha(A) =',SAPT%monA%ACAlpha
 write(lout,'(1x,a,f9.6)') 'ACAlpha(B) =',SAPT%monB%ACAlpha

 call e2ind_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2disp_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

 !              : coupled (for a coefficient)
 if(SAPT%monA%Cubic) SAPT%monA%ACAlpha=SAPT%monA%ACAlpha2
 if(SAPT%monB%Cubic) SAPT%monB%ACAlpha=SAPT%monB%ACAlpha2
 write(lout,'(/1x,a,f9.6)') 'ACAlpha(A) =',SAPT%monA%ACAlpha
 write(lout,'(1x,a,f9.6)') 'ACAlpha(B) =',SAPT%monB%ACAlpha

 call e2ind_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2disp_cpld(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exind(Flags,SAPT%monA,SAPT%monB,SAPT)
 call e2exdisp(Flags,SAPT%monA,SAPT%monB,SAPT)

 call e2_cubic(Flags,SAPT%monA,SAPT%monB,SAPT)

 call summary_sapt(SAPT)
 !call summary_sapt_verbose(SAPT)

 call print_warn(SAPT)
 call free_sapt(Flags,SAPT)

 stop

end subroutine sapt_extrapol

subroutine sapt_basinfo(SAPT,NBasis)
implicit none

type(SaptData)      :: SAPT
integer,intent(out) :: NBasis

 NBasis = 0
 if(SAPT%InterfaceType==1) then
    call basinfo(NBasis,'SIRIUS_A.RST','DALTON')
 elseif(SAPT%InterfaceType==2) then
    call basinfo(NBasis,'AOONEINT_A','MOLPRO')
 endif
 if(NBasis==0.and.SAPT%monA%NBasis==0) then
    write(LOUT,'(1x,a)') 'ERROR! NBasis NOWHERE TO BE FOUND!'
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

end subroutine sapt_basinfo

subroutine sapt_response(Flags,Mon,EnChck,NBasis)
implicit none

type(FlagsData)    :: Flags
type(SystemBlock)  :: Mon
integer,intent(in) :: NBasis
logical,intent(in) :: EnChck

integer          :: i,j,ij
integer          :: SaptLevel
logical          :: regular,extrapolate
double precision :: MO(NBasis*NBasis)

 SaptLevel = Flags%SaptLevel
 if(Mon%Cubic) then
    extrapolate = .true.
    regular     = .false.
 else
    extrapolate = .false.
    regular     = .true.
 endif

 ij = 0
 do j=1,NBasis
    do i=1,NBasis
       ij = ij + 1
       MO(ij) = Mon%CMO(i,j)
    enddo
 enddo

! calculate response
 call SaptInter(NBasis,Mon,Flags%ICASSCF)

 ! prepare RDM2
 if(Flags%ICASSCF==1) then
    call read2rdm(Mon,NBasis)
    if(Mon%Monomer==1) call system('cp rdm2_A.dat rdm2.dat')
    if(Mon%Monomer==2) call system('cp rdm2_B.dat rdm2.dat')
 endif
 call prepare_RDM2val(Mon,Flags%ICASSCF,NBasis)

 if(SaptLevel.eq.1) return

 if(SaptLevel.eq.999 .or. SaptLevel.eq.666) then
    ! for RSPT2 only AB matrices are needed
    call calc_ab_cas(Mon,MO,Flags,NBasis)
    return
 endif

 if(Flags%ISERPA==0) then
    if(SaptLevel.eq.0.or.SaptLevel.eq.10) then

        if(SaptLevel.eq.10) Flags%IFlag0 = 1
        call calc_resp_unc(Mon,MO,Flags,NBasis)

    elseif(SaptLevel.gt.0) then

       if(regular) then

          if(Flags%IFunSR/=0) then
             call calc_resp_dft(Mon,MO,Flags,NBasis)
          else
             print*, 'here?'
             call calc_resp_casgvb(Mon,MO,Flags,NBasis,EnChck)
          endif

       elseif(extrapolate) then

            call calc_resp_extrapolate(Mon,Flags,NBasis)

       endif

    endif
 elseif(Flags%ISERPA==2) then
    call calc_resp_pino(Mon,MO,Flags,NBasis)
   ! call calc_resp_casgvb(Mon,MO,Flags,NBasis,EnChck)
 endif

end subroutine sapt_response

subroutine sapt_ab_ints(Flags,A,B,iPINO,NBasis)
implicit none

type(FlagsData)    :: Flags
type(SystemBlock)  :: A,B
type(AOReaderData) :: reader
integer,intent(in) :: iPINO,NBasis
integer            :: thr_id
integer            :: ntr,iunit_aotwosort

if(Flags%SaptLevel==999) then
  print*, 'In RSPT2 intermoner ints not calculated for now'
  ! deallocate (FF|NCholesky) integrals
  deallocate(A%FF)
  deallocate(B%FF)
  return
endif

if(Flags%SaptLevel==666) then
  print*, 'Get 1-st order exchange ints (RSPT2+)'
  if (Flags%ICholesky==1) then
     call chol_ints_oooo(A%num0+A%num1,A%num0+A%num1,A%OO,&
                         B%num0+B%num1,B%num0+B%num1,B%OO,&
                         A%NChol,'OOOOAABB')
     !call chol_ints_fofo(NBasis,A%num0+A%num1,A%FF,&
     !                    NBasis,B%num0+B%num1,B%FF,&
     !                    A%NChol,NBasis,'FOFOAABB')

     call chol_ints_oooo(B%num0+B%num1,B%num0+B%num1,B%OO,  &
                         B%num0+B%num1,A%num0+A%num1,B%OOBA,&
                         A%NChol,'OOOOBBBA')

     call chol_ints_oooo(A%num0+A%num1,A%num0+A%num1,A%OO,  &
                         A%num0+A%num1,B%num0+B%num1,A%OOAB,&
                         A%NChol,'OOOOAAAB')
  endif
  ! deallocate (FF|NCholesky) integrals
  deallocate(A%FF)
  deallocate(B%FF)
  return
endif

if(Flags%ISERPA==0) then
  if(Flags%ICASSCF==1) then
     if(A%num1/=A%NAct) then
        write(LOUT,'(/1x,a)') 'Warning! Active orbitals in A have either 1.0 or 0.0 occupancy!'
        write(LOUT,'(1x,a/)') 'Correct for the Hartree-Fock reference, proceeding...'
        !stop
        A%num0 = A%INAct
        A%num1 = A%NAct
        A%num2 = NBasis - A%NAct - A%INAct
     endif

     if(B%num1/=B%NAct) then
        write(LOUT,'(/1x,a)') 'Warning! Active orbitals in B have either 1.0 or 0.0 occupancy!'
        write(LOUT,'(1x,a/)') 'Correct for the Hartree-Fock reference, proceeding...'
        !stop
        B%num0 = B%INAct
        B%num1 = B%NAct
        B%num2 = NBasis - B%NAct - B%INAct

     endif
  endif

  if(Flags%SaptLevel/=1.and.Flags%ICholesky==0) then
     ! integrals stored as (ov|ov)
     call tran4_gen(NBasis,&
                    A%num0+A%num1,A%CMO,&
                    A%num1+A%num2,A%CMO(1:NBasis,A%num0+1:NBasis),&
                    B%num0+B%num1,B%CMO,&
                    B%num1+B%num2,B%CMO(1:NBasis,B%num0+1:NBasis),&
                    'TWOMOAB','AOTWOSORT')

  endif

  write(LOUT,'(/1x,a)') 'Transforming E2exch-ind integrals...'

  ! omp tasks
  ! Let's open the AOTWSORT only once, therefore we open it here

  if(Flags%ICholesky==1) then
     !print*, 'dimOA',A%num0+A%num1
     !print*, 'dimOB',B%num0+B%num1
     ! term A3-ind
     call chol_ints_fofo(NBasis,A%num0+A%num1,A%FF,&
                         NBasis,B%num0+B%num1,B%FF,&
                         A%NChol,NBasis,'FOFOAABB')

     ! term A1-ind
     call chol_ints_fofo(NBasis,NBasis,A%FFAB, &
                         A%num0+A%num1,B%num0+B%num1,A%FFAB,&
                         A%NChol,NBasis,'FFOOABAB')

     ! term A2-ind
     call chol_ints_fofo(NBasis,B%num0+B%num1,B%FF,  &
                         NBasis,A%num0+A%num1,B%FFBA,&
                         A%NChol,NBasis,'FOFOBBBA')
     call chol_ints_fofo(NBasis,A%num0+A%num1,A%FF,  &
                         NBasis,B%num0+B%num1,A%FFAB,&
                         A%NChol,NBasis,'FOFOAAAB')
     ! A2A(B): YY
     call chol_ints_fofo(NBasis,B%num0+B%num1,B%FF,  &
                         NBasis,B%num0+B%num1,A%FFAB,&
                         A%NChol,NBasis,'FOFOBBAB')
     call chol_ints_fofo(NBasis,A%num0+A%num1,A%FF,  &
                         NBasis,A%num0+A%num1,B%FFBA,&
                         A%NChol,NBasis,'FOFOAABA')
  else
   ! working stuff...
   !  ! term A3-ind
   !  call tran4_gen(NBasis,&
   !           NBasis,B%CMO,&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           NBasis,A%CMO,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           'FOFOAABB','AOTWOSORT')
   !  ! term A1-ind
   !  call tran4_gen(NBasis,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           NBasis,A%CMO,&
   !           NBasis,B%CMO,&
   !           'FFOOABAB','AOTWOSORT')
   !  ! term A2-ind
   !  ! A2A(B): XX
   !  call tran4_gen(NBasis,&
   !           NBasis,B%CMO,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           NBasis,B%CMO,&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           'FOFOBBBA','AOTWOSORT')
   !  call tran4_gen(NBasis,&
   !           NBasis,A%CMO,&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           NBasis,A%CMO,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           'FOFOAAAB','AOTWOSORT')
   !  ! A2A(B): YY
   !  call tran4_gen(NBasis,&
   !           NBasis,A%CMO,&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           NBasis,B%CMO,&
   !           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
   !           'FOFOBBAB','AOTWOSORT')
   !  call tran4_gen(NBasis,&
   !           NBasis,B%CMO,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           NBasis,A%CMO,&
   !           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
   !           'FOFOAABA','AOTWOSORT')

      call reader%open('AOTWOSORT')

   ! Gianfranco's OMP modification
      write(LOUT,'(/1x,a)') 'Transforming E2exch-ind integrals...'
      write(LOUT,'(1x,a)')  '    (the OMP version is active)     '
      ! - remove `default(shared)`
      !$omp parallel default(shared) private(thr_id)
      !!$print *, "DEBUG: omp num threads: ", omp_get_num_threads()
         !$omp single
         !$omp task
         !$ thr_id = omp_get_thread_num()
      ! term A3-ind
      call new_tran4_gen(NBasis,&
               NBasis,B%CMO,&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               NBasis,A%CMO,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               'FOFOAABB',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp task
         !$ thr_id = omp_get_thread_num()
      ! term A1-ind
      call new_tran4_gen(NBasis,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               NBasis,A%CMO,&
               NBasis,B%CMO,&
               'FFOOABAB',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp task
         !$ thr_id = omp_get_thread_num()
      ! term A2-ind
      ! A2A(B): XX
      call new_tran4_gen(NBasis,&
               NBasis,B%CMO,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               NBasis,B%CMO,&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               'FOFOBBBA',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp task
         !$ thr_id = omp_get_thread_num()
      call new_tran4_gen(NBasis,&
               NBasis,A%CMO,&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               NBasis,A%CMO,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               'FOFOAAAB',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp task
         !$ thr_id = omp_get_thread_num()
      !! A2A(B): YY
      call new_tran4_gen(NBasis,&
               NBasis,A%CMO,&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               NBasis,B%CMO,&
               B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
               'FOFOBBAB',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp task
         !$ thr_id = omp_get_thread_num()
      call new_tran4_gen(NBasis,&
               NBasis,B%CMO,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               NBasis,A%CMO,&
               A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
               'FOFOAABA',reader,iunit_aotwosort,thr_id)
         !$omp end task
         !$omp end single
         !$omp end parallel

      call reader%close
  endif

     write(LOUT,'(/1x,a)') 'Transforming E2exch-disp integrals...'
     if(Flags%ICholesky==1) then
        call chol_ints_fofo(NBasis,B%num0+B%num1,A%FFAB,&
                            NBasis,A%num0+A%num1,B%FFBA,&
                            A%NChol,NBasis,'FOFOABBA')
        ! XY and YX, A2
        call chol_ints_fofo(NBasis,NBasis,A%FFAB, &
                            B%num0+B%num1,B%num0+B%num1,B%FF,&
                            A%NChol,NBasis,'FFOOABBB')
        call chol_ints_fofo(NBasis,NBasis,B%FFBA, &
                            A%num0+A%num1,A%num0+A%num1,A%FF,&
                            A%NChol,NBasis,'FFOOBAAA')
        ! deallocate (FF|NCholesky) integrals
        deallocate(A%FF)
        deallocate(B%FF)
     else
        call tran4_gen(NBasis,&
                 NBasis,B%CMO,&
                 A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
                 NBasis,A%CMO,&
                 B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
                 'FOFOABBA','AOTWOSORT')
        ! XY and YX, A2
        call tran4_gen(NBasis,&
                 B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
                 B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
                 NBasis,A%CMO,&
                 NBasis,B%CMO,&
                 'FFOOABBB','AOTWOSORT')
        call tran4_gen(NBasis,&
                 A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
                 A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
                 NBasis,B%CMO,&
                 NBasis,A%CMO,&
                 'FFOOBAAA','AOTWOSORT')

     endif

  !endif

  ! <oo|oo>
     if(Flags%ICholesky==1) then
        call chol_ints_gen(B%num0+B%num1,B%num0+B%num1,B%OO,&
                           A%num0+A%num1,A%num0+A%num1,A%OO,A%NChol,'TMPOOAB')
     else
        call tran4_gen(NBasis,&
                     A%num0+A%num1,A%CMO,&
                     A%num0+A%num1,A%CMO,&
                     B%num0+B%num1,B%CMO,&
                     B%num0+B%num1,B%CMO,&
                    'TMPOOAB','AOTWOSORT')
     endif

elseif(Flags%ISERPA==2) then

   write(LOUT,'(/,1x,a,i2)') 'Calculate AB integrals for iPINO =', iPINO
   if(iPINO==0.or.iPINO==1) then
      ! FCI
      call tran4_gen(NBasis,&
                     NBasis,A%CMO,&
                     NBasis,A%CMO,&
                     NBasis,B%CMO,&
                     NBasis,B%CMO,&
                     'TWOMOAB','AOTWOSORT')

   else
      ! CAS/LR, GVB(TD-APSG)
      call tran4_gen(NBasis,&
                    (A%num0+A%num1),A%CMO,&
                     NBasis,A%CMO,&
                    (B%num0+B%num1),B%CMO,&
                     NBasis,B%CMO,&
                     'TWOMOAB','AOTWOSORT')

      ! <oo|oo>
      call tran4_gen(NBasis,&
                     A%num0+A%num1,A%CMO,&
                     A%num0+A%num1,A%CMO,&
                     B%num0+B%num1,B%CMO,&
                    B%num0+B%num1,B%CMO,&
                     'TMPOOAB','AOTWOSORT')

   endif

endif

end subroutine sapt_ab_ints

subroutine sapt_ab_ints_red(Flags,A,B,iPINO,NBasis,NBasisRed)
implicit none

type(FlagsData)    :: Flags
type(SystemBlock)  :: A,B
integer,intent(in) :: iPINO,NBasis,NBasisRed

if(Flags%ISERPA==0) then
  if(Flags%ICASSCF==1) then
     if(A%num1/=A%NAct) then
        write(LOUT,'(/1x,a)') 'Warning! Active orbitals in A have either 1.0 or 0.0 occupancy!'
        write(LOUT,'(1x,a/)') 'Correct for the Hartree-Fock reference, proceeding...'
        !stop
        A%num0 = A%INAct
        A%num1 = A%NAct
        A%num2 = NBasis - A%NAct - A%INAct
     endif

     if(B%num1/=B%NAct) then
        write(LOUT,'(/1x,a)') 'Warning! Active orbitals in B have either 1.0 or 0.0 occupancy!'
        write(LOUT,'(1x,a/)') 'Correct for the Hartree-Fock reference, proceeding...'
        !stop
        B%num0 = B%INAct
        B%num1 = B%NAct
        B%num2 = NBasis - B%NAct - B%INAct

     endif
  endif

  ! integrals stored as (ov|ov)
  print*, 'NBASISRED',NBasisRed

  write(LOUT,'(/1x,a)') 'Transforming E2exch-ind integrals...'
  ! term A3-ind
  call tran4_gen(NBasis,&
           NBasis,B%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           NBasis,A%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           'FOFOAABB','AOTWOSORT')
  ! term A1-ind
  call tran4_gen(NBasis,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           NBasis,A%CMO,&
           NBasis,B%CMO,&
           'FFOOABAB','AOTWOSORT')
  ! term A2-ind
  ! A2A(B): XX
  call tran4_gen(NBasis,&
           NBasis,B%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           NBasis,B%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           'FOFOBBBA','AOTWOSORT')
  call tran4_gen(NBasis,&
           NBasis,A%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           NBasis,A%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           'FOFOAAAB','AOTWOSORT')
  !! A2A(B): YY
  call tran4_gen(NBasis,&
           NBasis,A%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           NBasis,B%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           'FOFOBBAB','AOTWOSORT')
  call tran4_gen(NBasis,&
           NBasis,B%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           NBasis,A%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           'FOFOAABA','AOTWOSORT')

  write(LOUT,'(/1x,a)') 'Transforming E2exch-disp integrals...'
  call tran4_gen(NBasis,&
           NBasis,B%CMO,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           NBasis,A%CMO,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           'FOFOABBA','AOTWOSORT')
  ! XY and YX, A2
  call tran4_gen(NBasis,&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           B%num0+B%num1,B%CMO(1:NBasis,1:(B%num0+B%num1)),&
           NBasis,A%CMO,&
           NBasis,B%CMO,&
           'FFOOABBB','AOTWOSORT')
  call tran4_gen(NBasis,&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           A%num0+A%num1,A%CMO(1:NBasis,1:(A%num0+A%num1)),&
           NBasis,B%CMO,&
           NBasis,A%CMO,&
           'FFOOBAAA','AOTWOSORT')

  call tran4_gen(NBasis,&
                 A%num0+A%num1,A%CMO,&
                 A%num1+A%num2,A%CMO(1:NBasis,A%num0+1:NBasisRed),&
                 B%num0+B%num1,B%CMO,&
                 B%num1+B%num2,B%CMO(1:NBasis,B%num0+1:NBasisRed),&
                 'TWOMOAB','AOTWOSORT')

  ! <oo|oo>
  call tran4_gen(NBasis,&
                 A%num0+A%num1,A%CMO,&
                 A%num0+A%num1,A%CMO,&
                 B%num0+B%num1,B%CMO,&
                 B%num0+B%num1,B%CMO,&
                 'TMPOOAB','AOTWOSORT')
!
! have not thought that thru
!elseif(Flags%ISERPA==2) then
!
!   write(LOUT,'(/,1x,a,i2)') 'Calculate AB integrals for iPINO =', iPINO
!   if(iPINO==0.or.iPINO==1) then
!      ! FCI
!      call tran4_gen(NBasis,&
!                     NBasis,A%CMO,&
!                     NBasis,A%CMO,&
!                     NBasis,B%CMO,&
!                     NBasis,B%CMO,&
!                     'TWOMOAB','AOTWOSORT')
!
!   else
!      ! CAS/LR, GVB(TD-APSG)
!      call tran4_gen(NBasis,&
!                    (A%num0+A%num1),A%CMO,&
!                     NBasis,A%CMO,&
!                    (B%num0+B%num1),B%CMO,&
!                     NBasis,B%CMO,&
!                     'TWOMOAB','AOTWOSORT')
!   endif

endif

end subroutine sapt_ab_ints_red

subroutine prepare_RDM2val(Mon,ICASSCF,NBasis)
!
! prepare RDM2val(NOccup,NOccup,NOccup,NOccup) 
! from packed RDM2(NAddr)
!
implicit none

type(SystemBlock)  :: Mon
integer,intent(in) :: ICASSCF
integer,intent(in) :: NBasis

integer :: i,j,k,l
integer :: NOccup
double precision, external :: FRDM2,FRDM2GVB

NOccup = Mon%num0+Mon%num1
if(allocated(Mon%RDM2val)) deallocate(Mon%RDM2val)
allocate(Mon%RDM2val(NOccup,NOccup,NOccup,NOccup))

if(ICASSCF==1) then
! CAS
   do l=1,NOccup
      do k=1,NOccup
         do j=1,NOccup
            do i=1,NOccup
               Mon%RDM2val(i,j,k,l) = &
               FRDM2(i,k,j,l,Mon%RDM2,Mon%Occ,Mon%Ind2,Mon%NAct,NBasis)
            enddo
         enddo
      enddo
   enddo
elseif(ICASSCF==0) then
! GVB
   do l=1,NOccup
     do k=1,NOccup
        do j=1,NOccup
           do i=1,NOccup
              Mon%RDM2val(i,j,k,l) = FRDM2GVB(i,k,j,l,Mon%Occ,NBasis)
           enddo
        enddo
     enddo
   enddo
endif

if(allocated(Mon%RDM2)) deallocate(Mon%RDM2)

end subroutine prepare_RDM2val

subroutine reduce_virt(Flags,Mon,NBas)
 implicit none

 type(SystemBlock)  :: Mon
 type(FlagsData)    :: Flags
 integer,intent(in) :: NBas

 integer            :: i,j,ij
 integer            :: NInte1,NInte2,NVirt
 double precision   :: URe(NBas,NBas)

 double precision,allocatable :: XOne(:),MO(:),TwoMO(:)
 character(:),allocatable     :: onefile,twofile,twojfile,twokfile
 double precision,allocatable :: Eps(:,:),workSq(:,:)

 write(LOUT,'(/,1x,a)') 'TRUNCATE VIRTUAL ORBITAL SPACE'
 write(LOUT,'(8a10)') ('**********',i=1,6)

 ! reduce virt space
 ! set filenames
 if(Mon%Monomer==1) then
    onefile  = 'ONEEL_A'
    twofile  = 'TWOMOAA'
    twojfile = 'FFOOAA'
    twokfile = 'FOFOAA'
 elseif(Mon%Monomer==2) then
    onefile  = 'ONEEL_B'
    twofile  = 'TWOMOBB'
    twojfile = 'FFOOBB'
    twokfile = 'FOFOBB'
 endif

 ! set dimensions
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

 URe = 0d0
 do i=1,NBas
    URe(i,i) = 1d0
 enddo

 allocate(XOne(NInte1),MO(NBas**2))
 if(Mon%TwoMoInt==TWOMO_INCORE) allocate(TwoMO(NInte2))

 ij = 0
 do j=1,NBas
    do i=1,NBas
       ij = ij + 1
       MO(ij) = Mon%CMO(i,j)
    enddo
 enddo

 ! prepare RDM2
 if(Flags%ICASSCF==1) then
    call read2rdm(Mon,NBas)
    if(Mon%Monomer==1) call system('cp rdm2_A.dat rdm2.dat')
    if(Mon%Monomer==2) call system('cp rdm2_B.dat rdm2.dat')
 endif

 ! read 1-el
 call get_1el_h_mo(XOne,MO,NBas,onefile)

 ! INCORE: load 2-el integrals
 if(Mon%TwoMoInt==TWOMO_INCORE) call LoadSaptTwoNO(Mon%Monomer,TwoMO,NBas,NInte2)

 allocate(Eps(NBas,NBas),workSq(NBas,NBas))
 workSq = transpose(Mon%CMO)

 select case(Mon%TwoMoInt)
 case(TWOMO_INCORE)
    call MP2RDM(TwoMO,Eps,Mon%Occ,URe,workSq,XOne,&
                Mon%IndN,Mon%IndX,Mon%IndAux,Mon%NDimX,&
                NBas,Mon%NDim,NInte1,NInte2,NVirt,&
                twofile,Mon%ThrVirt,.false.)

 case(TWOMO_FOFO)

    call MP2RDM_FOFO(Mon%PerVirt,Eps,Mon%Occ,URe,workSq,XOne,&
                     Mon%IndN,Mon%IndX,Mon%IndAux,Mon%IGem,  &
                     Mon%NAct,Mon%INAct,Mon%NDimX,Mon%NDim,NBas,NInte1,&
                     twojfile,twokfile,Mon%ThrVirt,Mon%NVZero,Mon%IPrint)

 case(TWOMO_FFFF)

    write(LOUT,'(1x,a)') 'ERROR! REDVIRT NOT AVAILABLE WITH FFFF (FULL) TWOMOINT'
    stop

 end select

 ! get AO->NO(d) matrix
 call dgemm('N','T',NBas,NBas,NBas,1d0,MO,NBas,Eps,NBas,0d0,workSq,NBas)
 do j=1,NBas
      do i=1,NBas
         Mon%CMO(i,j) = workSq(i,j)
      enddo
 enddo

 deallocate(workSq,Eps)
 deallocate(MO,XOne)
 if(Mon%TwoMoInt==TWOMO_INCORE) deallocate(TwoMO)

end subroutine reduce_virt

subroutine sapt_mon_ints(Mon,Flags,NBas)
implicit none

type(SystemBlock) :: Mon
type(FlagsData) :: Flags

integer :: NBas
integer :: i,j,ij,ione
integer :: NSq,NInte1,NInte2

double precision             :: URe(NBas,NBas),MO(NBas*NBas)
double precision,allocatable :: TwoMO(:)
double precision,allocatable :: work1(:),work2(:),XOne(:)
character(8)                 :: label
character(:),allocatable     :: onefile,twofile
character(:),allocatable     :: twojfile,twokfile
!test
double precision :: Tcpu,Twall

call clock('START',Tcpu,Twall)

! set dimensions
 NSq = NBas**2
 NInte1 = NBas*(NBas+1)/2
 NInte2 = NInte1*(NInte1+1)/2

! set file names
 if(Mon%Monomer==1) then
    onefile  = 'ONEEL_A'
    twofile  = 'TWOMOAA'
    twojfile = 'FFOOAA'
    twokfile = 'FOFOAA'
 elseif(Mon%Monomer==2) then
    onefile  = 'ONEEL_B'
    twofile  = 'TWOMOBB'
    twojfile = 'FFOOBB'
    twokfile = 'FOFOBB'
 endif

 allocate(work1(NSq),work2(NSq),XOne(NInte1))

 URe = 0d0
 do i=1,NBas
    URe(i,i) = 1d0
 enddo

 MO = 0
 ij = 0
 do j=1,NBas
 do i=1,NBas
    ij = ij + 1
    MO(ij) = Mon%CMO(i,j)
 enddo
 enddo

 ! read 1-el
 call get_1el_h_mo(XOne,MO,NBas,onefile)

 if(Flags%SaptLevel==1) return

 ! transform 2-el integrals

 select case(Mon%TwoMoInt)
 case(TWOMO_INCORE,TWOMO_FFFF)
   ! full - for GVB and CAS
   call tran4_full(NBas,MO,MO,twofile,'AOTWOSORT')

 case(TWOMO_FOFO)
   if(Flags%ICholesky==1) then
      call chol_ints_fofo(NBas,NBas,Mon%FF, &
                     Mon%num0+Mon%num1,Mon%num0+Mon%num1,Mon%FF,&
                     Mon%NChol,NBas,twojfile)
      call chol_ints_fofo(NBas,Mon%num0+Mon%num1,Mon%FF,&
                     NBas,Mon%num0+Mon%num1,Mon%FF,&
                     Mon%NChol,NBas,twokfile)
      !call chol_ints_gen(NBas,NBas,Mon%FF, &
      !               Mon%num0+Mon%num1,Mon%num0+Mon%num1,Mon%OO,&
      !               Mon%NChol,twojfile)
      !call chol_ints_gen(NBas,Mon%num0+Mon%num1,Mon%FO,&
      !               NBas,Mon%num0+Mon%num1,Mon%FO,&
      !               Mon%NChol,twokfile)
   else
      ! transform J and K
       call tran4_gen(NBas,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            NBas,MO,&
            NBas,MO,&
            twojfile,'AOTWOSORT')
       call clock('FFOO',Tcpu,Twall)
       call tran4_gen(NBas,&
            NBas,MO,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            NBas,MO,&
            Mon%num0+Mon%num1,MO(1:NBas*(Mon%num0+Mon%num1)),&
            twokfile,'AOTWOSORT')
       call clock('FOFO',Tcpu,Twall)
    endif
 end select

 deallocate(XOne,work2,work1)
 !if(Mon%TwoMoInt==1) deallocate(TwoMO)

end subroutine sapt_mon_ints

subroutine chol_ints_gen(nA,nB,MatAB,nC,nD,MatCD,NCholesky,fname)
! MatAB and MatCD may have any NChol,XX shape
implicit none

integer,intent(in) :: nA,nB,nC,nD,NCholesky
character(*),intent(in)     :: fname
double precision,intent(in) :: MatAB(NCholesky,nA*nB), &
                               MatCD(NCholesky,nC*nD)

integer :: iunit
integer :: nAB,nCD,cd
integer :: i,ic, id, icd
double precision,allocatable :: work(:)

nAB = nA*nB
nCD = nC*nD

allocate(work(nAB))

print*, 'Assemble ',fname,' from Cholesky Vectors'

open(newunit=iunit,file=fname,status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*nAB)

! (FF|OO)
do cd=1,nCD
   call dgemv('T',NCholesky,nAB,1d0,MatAB,NCholesky,MatCD(:,cd),1,0d0,work,1)
   write(iunit,rec=cd) work(1:nAB)
enddo

deallocate(work)
close(iunit)

end subroutine chol_ints_gen

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

subroutine summary_sapt(SAPT)
implicit none

type(SaptData) :: SAPT

integer          :: i,j
double precision :: esapt2

SAPT%esapt2 = SAPT%elst + SAPT%exchs2 + SAPT%e2ind &
              + SAPT%e2exind + SAPT%e2disp + SAPT%e2exdisp

SAPT%esapt0 = SAPT%elst + SAPT%exchs2 + SAPT%e2ind_unc &
              + SAPT%e2exind_unc + SAPT%e2disp_unc + SAPT%e2exdisp_unc

write(LOUT,'(/,8a10)') ('**********',i=1,4)
write(LOUT,'(1x,a)') 'SAPT SUMMARY / milliHartree'
write(LOUT,'(8a10)') ('**********',i=1,4)

write(LOUT,'(1x,a,i3)') 'SAPT level  =', SAPT%SaptLevel

write(LOUT,'(1x,a,t19,a,f16.8)') 'E1elst',    '=', SAPT%elst*1.d03
write(LOUT,'(1x,a,t19,a,f16.8)') 'E1exch(S2)','=', SAPT%exchs2*1.d03

if(SAPT%SaptLevel==2) then
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2ind',      '=', SAPT%e2ind*1.d03
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2exch-ind', '=', SAPT%e2exind*1.0d3
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2disp',     '=', SAPT%e2disp*1.d03
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2exch-disp','=', SAPT%e2exdisp*1.0d3
   write(LOUT,'(1x,a,t19,a,f16.8)') 'Eint(SAPT2)','=', SAPT%esapt2*1.0d3

   if(SAPT%Wexcit) then
     write(lout,'(/1x,a)') 'Dexcitation corrections'
     j = SAPT%monA%InSt(1,1)
     do i=1,size(SAPT%Wind)
        write(lout,'(1x,a,2i1,a,f12.6)') &
             'Wind_',SAPT%monA%InSt(1,1)-1,j,'     =', SAPT%Wind(i)*1000d0
        j = j + 1
     enddo
     j = SAPT%monA%InSt(1,1)
     do i=1,size(SAPT%Wdisp)
        write(lout,'(1x,a,2i1,a,f12.6)') &
             'Wdisp_',SAPT%monA%InSt(1,1)-1,j,'    =', SAPT%Wdisp(i)*1000d0
        j = j + 1
     enddo
   endif

elseif(SAPT%SaptLevel==1) then
   write(LOUT,'(1x,a,t19,a,f16.8)') 'Eint(SAPT1)', '=', SAPT%esapt2*1.0d3

elseif(SAPT%SaptLevel==10) then
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2disp(CAS)', '=', SAPT%e2dispinCAS*1.d03

elseif(SAPT%SaptLevel==0) then
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2ind(unc)',  '=', SAPT%e2ind_unc*1.d03
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2exch-ind',  '=', SAPT%e2exind_unc*1.0d3
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2disp(unc)', '=', SAPT%e2disp_unc*1.d03
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E2exch-disp', '=', SAPT%e2exdisp_unc*1.0d3
   write(LOUT,'(1x,a,t19,a,f16.8)') 'Eint(SAPT0)', '=', SAPT%esapt0*1.0d3
endif

end subroutine summary_sapt

subroutine summary_sapt_verbose(SAPT)
! print all digits
implicit none

type(SaptData) :: SAPT

integer :: i
double precision :: esapt2

SAPT%esapt2 = SAPT%elst + SAPT%exchs2 + SAPT%e2ind &
              + SAPT%e2exind + SAPT%e2disp + SAPT%e2exdisp

write(LOUT,'(/,8a10)') ('**********',i=1,4)
write(LOUT,'(1x,a)') 'SAPT SUMMARY / milliHartree'
write(LOUT,'(8a10)') ('**********',i=1,4)

write(LOUT,'(1x,a,i2)') 'SAPT level  =', SAPT%SaptLevel

write(LOUT,*) 'E1elst      = ', SAPT%elst*1.d03
write(LOUT,*) 'E1exch(S2)  = ', SAPT%exchs2*1.d03

if(SAPT%SaptLevel/=1) then
   write(LOUT,*) 'E2ind       = ', SAPT%e2ind*1.d03
   write(LOUT,*) 'E2exch-ind  = ', SAPT%e2exind*1.0d3
   write(LOUT,*) 'E2disp      = ', SAPT%e2disp*1.d03
   write(LOUT,*) 'E2exch-disp = ', SAPT%e2exdisp*1.0d3
   write(LOUT,*) 'Eint(SAPT2) = ', SAPT%esapt2*1.0d3
elseif(SAPT%SaptLevel==1) then
   write(LOUT,*) 'Eint(SAPT1) = ', SAPT%esapt2*1.0d3
endif

end subroutine summary_sapt_verbose

subroutine summary_rspt(SAPT)
implicit none

type(SaptData) :: SAPT

integer :: i
double precision :: erspt2

erspt2 = SAPT%elst +  SAPT%e2ind + SAPT%e2disp
if(SAPT%SaptLevel==666) then ! RSPT2+
  erspt2 = erspt2 + SAPT%exchs2
endif

write(LOUT,'(/,8a10)') ('**********',i=1,4)
write(LOUT,'(1x,a)') 'SAPT SUMMARY / milliHartree'
write(LOUT,'(8a10)') ('**********',i=1,4)

write(LOUT,'(1x,a)') 'SAPT level  = RSPT2'

write(LOUT,'(1x,a,t19,a,f16.8)') 'E1elst',     '=', SAPT%elst  *1.d03
if(SAPT%SaptLevel==666) then ! RSPT2+
   write(LOUT,'(1x,a,t19,a,f16.8)') 'E1exch',     '=', SAPT%exchs2*1.d03
endif
write(LOUT,'(1x,a,t19,a,f16.8)') 'E2ind',      '=', SAPT%e2ind *1.d03
write(LOUT,'(1x,a,t19,a,f16.8)') 'E2disp',     '=', SAPT%e2disp*1.d03
write(LOUT,'(1x,a,t19,a,f16.8)') 'Eint(RSPT2)','=', erspt2*1.0d3
write(LOUT,'()')

end subroutine summary_rspt

subroutine free_sapt(Flags,SAPT)
implicit none

type(FlagsData) :: Flags
type(SaptData)  :: SAPT

integer         :: ISERPA

deallocate(SAPT%monA%CICoef,SAPT%monA%IGem,SAPT%monA%Occ, &
           SAPT%monA%IndAux,SAPT%monA%IndX,SAPT%monA%IndN,&
           SAPT%monA%CMO,&
           SAPT%monA%IPair)
deallocate(SAPT%monB%CICoef,SAPT%monB%IGem,SAPT%monB%Occ, &
           SAPT%monB%IndAux,SAPT%monB%IndX,SAPT%monB%IndN,&
           SAPT%monB%CMO,&
           SAPT%monB%IPair)

! for PINO
ISERPA = 0

! symmetry matrices
if(allocated(SAPT%monA%NumOSym)) then
   deallocate(SAPT%monA%NumOSym)
endif
if(allocated(SAPT%monB%NumOSym)) then
   deallocate(SAPT%monB%NumOSym)
endif

if(allocated(SAPT%monB%Kmat)) then
  deallocate(SAPT%monB%Kmat)
endif
! test Cholesky
if(allocated(SAPT%monA%DChol)) then
   deallocate(SAPT%monA%Dchol)
endif
if(allocated(SAPT%monB%DChol)) then
   deallocate(SAPT%monB%Dchol)
endif
if(allocated(SAPT%monA%OV)) then
   deallocate(SAPT%monA%OV)
endif
if(allocated(SAPT%monB%OV)) then
   deallocate(SAPT%monB%OV)
endif
if(allocated(SAPT%monA%OO)) then
   deallocate(SAPT%monA%OO)
endif
if(allocated(SAPT%monB%OO)) then
   deallocate(SAPT%monB%OO)
endif
if(allocated(SAPT%monA%FF)) then
   deallocate(SAPT%monA%FF)
endif
if(allocated(SAPT%monB%FF)) then
   deallocate(SAPT%monB%FF)
endif
if(allocated(SAPT%monA%FO)) then
   deallocate(SAPT%monA%FO)
endif
if(allocated(SAPT%monB%FO)) then
   deallocate(SAPT%monB%FO)
endif
if(allocated(SAPT%monA%FFAB)) then
   deallocate(SAPT%monA%FFAB)
endif
if(allocated(SAPT%monB%FFBA)) then
   deallocate(SAPT%monB%FFBA)
endif
! test Pmat
if(allocated(SAPT%CholVecs)) then
   deallocate(SAPT%CholVecs)
endif

! end test Cholesky

if(allocated(SAPT%monA%PP)) deallocate(SAPT%monA%PP)
if(allocated(SAPT%monB%PP)) deallocate(SAPT%monB%PP)

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
   deallocate(SAPT%monA%RDM2)
   deallocate(SAPT%monA%Ind2)
endif
if(allocated(SAPT%monB%RDM2)) then
   deallocate(SAPT%monB%RDM2)
   deallocate(SAPT%monB%Ind2)
endif

if(allocated(SAPT%monA%RDM2val)) then
   deallocate(SAPT%monA%RDM2val)
endif
if(allocated(SAPT%monB%RDM2val)) then
   deallocate(SAPT%monB%RDM2val)
endif

! for PINO only
if(allocated(SAPT%monA%IndNx)) then
   ISERPA = 1
   deallocate(SAPT%monA%IndNx)
endif
if(allocated(SAPT%monB%IndNx)) then
   deallocate(SAPT%monB%IndNx)
endif
if(allocated(SAPT%monA%IndNT)) then
   deallocate(SAPT%monA%IndNT)
endif
if(allocated(SAPT%monB%IndNT)) then
   deallocate(SAPT%monB%IndNT)
endif

! 2nd order exchange
if(allocated(SAPT%monA%Eig)) then
   deallocate(SAPT%monA%Eig,SAPT%monA%EigX,SAPT%monA%EigY)
endif
if(allocated(SAPT%monB%Eig)) then
   deallocate(SAPT%monB%Eig,SAPT%monB%EigX,SAPT%monB%EigY)
endif
! 2-RDM approximations
if(allocated(SAPT%monA%Fmat)) then
   deallocate(SAPT%monA%Fmat)
endif
if(allocated(SAPT%monB%Fmat)) then
   deallocate(SAPT%monB%Fmat)
endif

! c6 coeffs
if(allocated(SAPT%monA%dipm)) then
   deallocate(SAPT%monA%dipm)
endif
if(allocated(SAPT%monB%dipm)) then
   deallocate(SAPT%monB%dipm)
endif

!RSH
if(allocated(SAPT%monA%VCoul)) then
   deallocate(SAPT%monA%VCoul)
endif
if(allocated(SAPT%monB%VCoul)) then
   deallocate(SAPT%monB%VCoul)
endif

if(SAPT%doRSH) then
   call delfile('AOERFSORT')
   call delfile('FFOOERFAA')
   call delfile('FFOOERFBB')
   call delfile('FOFOERFAA')
   call delfile('FOFOERFBB')
   if(.not.SAPT%SameOm) call delfile('AOERFSORTB')
endif

! Wind,Wdisp terms
if(allocated(SAPT%Wind)) then
   deallocate(SAPT%Wind)
endif
if(allocated(SAPT%Wdisp)) then
   deallocate(SAPT%Wdisp)
endif

! cubic dispersion
if(SAPT%monA%Cubic) then
  call delfile('PROP_A0')
  call delfile('PROP_A1')
  call delfile('PROP_A2')
endif
if(SAPT%monB%Cubic) then
  call delfile('PROP_B0')
  call delfile('PROP_B1')
  call delfile('PROP_B2')
endif
! ....

! delete files
call delfile('AOTWOSORT')
if(SAPT%monA%TwoMoInt==TWOMO_INCORE.or.&
   SAPT%monA%TwoMoInt==TWOMO_FFFF) then
   call delfile('TWOMOAA')
endif
if(SAPT%monB%TwoMoInt==TWOMO_INCORE.or.&
   SAPT%monB%TwoMoInt==TWOMO_FFFF) then
   call delfile('TWOMOBB')
endif
call delfile ('ONEEL_A')
call delfile ('ONEEL_B')

call delfile('TMPOOAB')

if(SAPT%SaptLevel/=1) then
   call delfile('TWOMOAB')
endif

! we use them for iterative CERPA
if(.not.SAPT%doRSH) then
   call delfile('ABMAT_A')
   call delfile('ABMAT_B')
endif

if(Flags%ISERPA==0) then
   call delfile('FOFOABBA')
   call delfile('FOFOBBBA')
   call delfile('FOFOAAAB')
   call delfile('FOFOBBAB')
   call delfile('FOFOAABA')
   call delfile('FOFOAABB')
   call delfile('FFOOABAB')
   call delfile('FFOOABBB')
   call delfile('FFOOBAAA')

   if(.not.SAPT%monA%Cubic) call delfile('XY0_A')
   if(.not.SAPT%monB%Cubic) call delfile('XY0_B')
elseif(Flags%ISERPA==0.and.SAPT%SaptLevel==10) then
   call delfile('XY0_A')
   call delfile('XY0_B')
elseif(Flags%ISERPA==2) then
   call delfile('FFOOABAB')
   call delfile('FOFOAABB')
   call delfile('FOFOABBA')
   call delfile('FFFFABBB')
   call delfile('FFFFBAAA')
endif

if(SAPT%monA%TwoMoInt==TWOMO_FOFO) then
   call delfile('FFOOAA')
   call delfile('FOFOAA')
endif
if(SAPT%monB%TwoMoInt==TWOMO_FOFO) then
   call delfile('FFOOBB')
   call delfile('FOFOBB')
endif

if(.not.(SAPT%monA%Cubic.or.SAPT%monB%Cubic)) then
   call delfile('PROP_AB0')
endif

if(SAPT%SaptLevel==2) then
   if(.not.SAPT%monA%Cubic) call delfile('PROP_A')
   if(.not.SAPT%monB%Cubic) call delfile('PROP_B')
   call delfile('PROP_AB')
endif

if(SAPT%SemiCoupled) call delfile('PROP_A1')
if(SAPT%SemiCoupled) call delfile('PROP_B1')

end subroutine free_sapt

end module sapt_main

