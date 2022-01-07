module timing
! timing procedures adapted from Drake (M. Przybytek)
! with minor changes
implicit none

integer,save :: RefYear = 0

double precision,parameter :: sec      = 1.d0  
double precision,parameter :: millisec = sec * 1.D-3
double precision,parameter :: minute   = sec * 60.d0
double precision,parameter :: hour     = minute * 60.d0 
double precision,parameter :: day      = hour   * 24.d0 

integer,parameter :: string_out_LENGTH = 20

contains

subroutine clock(mode,Tcpu,Twall)
implicit none
! -What watch? -10 watch. -Such much?

character(*),intent(in) :: mode
double precision,intent(inout) :: Tcpu,Twall
double precision :: Tcpu_inter,Twall_inter
integer :: DateTime(8)
character(string_out_LENGTH) :: string_out

call date_and_time(values=DateTime)
if(RefYear==0) RefYear = DateTime(1)

call cpu_time(Tcpu_inter)
call wall_time(Twall_inter,DateTime)

if(adjustl(trim(mode))/='START') then

   string_out = mode(1:min(string_out_LENGTH,len(mode)))

   write(6,'(/,1x,3a)',advance='no') 'CPU  time in ',string_out,' is '
   call print_watch(Tcpu_inter-Tcpu)
   write(6,'(1x,3a)',advance='no')   'Wall time in ',string_out,' is '
   call print_watch(Twall_inter-Twall)

endif

Tcpu  = Tcpu_inter
Twall = Twall_inter

end subroutine clock

subroutine wall_time(time,DateTime)
implicit none

double precision :: time
integer,intent(in) :: DateTime(8)
integer :: imonth,iyear

time = millisec * DateTime(8) &
     + sec      * DateTime(7) &
     + minute   * DateTime(6) &
     + hour     * DateTime(5) &
     + day      *(DateTime(3) - 1)

do imonth=1,DateTime(2)-1
   select case(imonth)
   case(1,3,5,7,8,10,12)
      time = time + day * 31.d0 
   case(4,6,9,11)
      time = time + day * 30.d0 
   case(2)
      if(LeapYear(DateTime(1))) then
         time = time + day * 29.d0
      else
         time = time + day * 28.d0
      endif
   end select
enddo

do iyear=RefYear,DateTime(1)-1
   if(LeapYear(iyear)) then
      time = time + day * 366.d0 
   else
      time = time + day * 365.d0
   endif
enddo

end subroutine wall_time

function LeapYear(year)
implicit none

logical :: leapYear
integer :: year

if(mod(year,4)/=0) then
   leapYear = .false.
elseif(mod(year,100)/=0) then
   leapYear = .true.
elseif(mod(year,400)/=0) then
   leapYear = .false.
else
   leapYear = .true.
endif

end function LeapYear

subroutine print_watch(time)
implicit none

double precision :: time
integer :: chck_hour,chck_minute
double precision :: chck_second

chck_hour   = int(time/hour)
chck_minute = int((time - chck_hour*hour)/minute)
chck_second = time - chck_hour*hour - chck_minute*minute

if(chck_hour/=0) then
   write(6,'(i6,a2, i3,a2, f7.3,a2)') &
        chck_hour,  ' h', &
        chck_minute,' m', &
        chck_second,' s'
elseif(chck_minute/=0) then
   write(6,'(8x, i3,a2, f7.3,a2)') &
        chck_minute,' m', &
        chck_second,' s'
else
   write(6,'(8x, 5x, f7.3,a2)') &
        chck_second,' s'
endif

end subroutine print_watch

end module timing
