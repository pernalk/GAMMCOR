!> \file
!> Interfacing a new module in MOLPRO
!>
!> The mechanism for linking a completely new piece of code into Molpro
!> is as follows. In response to the \c USER data command, the command
!> subroutine leit (src/util/driver.f) simply calls the subroutine
!> \c user without arguments. The default stub in src/util/user.f
!> simply terminates. Thus the procedure is to create  a file user.f
!> in a Molpro develop directory, containing subroutine user which
!> then should run the code which is being introduced.
      subroutine user
      use dumpintao 
      implicit double precision (a-h,o-z)
      include "common/tapes"
      include "common/chirs"
      integer, parameter :: maxf = 5
      character(len = 16) :: field(maxf)
      character(len = 32) :: prdmft_job
      character(len = 8) :: units

      write(iout,*) 'IsRangeHybrid', IsRangeHybrid
      write(iout,*) 'IsChi', IsChi,ChiRS

       call getx(1, field, maxf)
       if(field(2) == "ORB") then
         call dump_grid()
       elseif(field(2) == "INTAO") then
         call dump_aoints 
       elseif(field(2) == "PRDMFT") then
         call dump_grid()
         call dump_aoints 
       elseif(field(2) == "WRK") then
         call dump_test('AOTWOINT.mol')
       ! PTDMFT-SAPT interface
       elseif(field(2) == "PRDMFT-SAPT") then
         if(trim(field(3))=='WRK') then
         call dump_test('AOTWOINT.mol')
         elseif(trim(field(3))=='GRID') then
           call dump_grid() 
         elseif(trim(field(3))=='MONOMERA') then
           call dump_molpro_sapt(1,IsRangeHybrid)
         else
           call dump_molpro_sapt(2,IsRangeHybrid)
         endif
       else
         write(iout,10)
 10   format(' USER command: Please provide keyword for USER')

      call fehler("error")
      return
      endif

      end
