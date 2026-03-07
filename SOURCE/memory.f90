module memory
!use precision
implicit none

integer,parameter :: prec    = kind(0.d0)
integer,parameter :: longint = selected_int_kind(18)

integer,parameter :: log_bytesize     = storage_size(.false.  )/8
integer,parameter :: int_bytesize     = storage_size(0        )/8
integer,parameter :: real_bytesize    = storage_size(0._prec  )/8

private
public mem_report
public set_parallel_memory,clear_parallel_memory
public init_thread_memory,collect_thread_memory
public mem_alloc,mem_dealloc

interface mem_alloc
module procedure &
     alloc_log1,alloc_log2,&
     alloc_int1,alloc_int2,&
     alloc_real1,alloc_real2,alloc_real3,alloc_real4,&
     alloc_log1b,alloc_log2b,&
     alloc_int1b,alloc_int2b,&
     alloc_real1b,alloc_real2b,alloc_real3b,alloc_real4b
end interface mem_alloc

interface mem_dealloc
module procedure &
     dealloc_log1,dealloc_log2,&
     dealloc_int1,dealloc_int2,&
     dealloc_real1,dealloc_real2,dealloc_real3,dealloc_real4
end interface mem_dealloc

integer(longint),save :: max_allocated_log  = 0_longint
integer(longint),save :: max_allocated_int  = 0_longint
integer(longint),save :: max_allocated_real = 0_longint
integer(longint),save :: allocated_log  = 0_longint
integer(longint),save :: allocated_int  = 0_longint
integer(longint),save :: allocated_real = 0_longint

logical,save :: parallel_memory = .false.

integer(longint),save :: local_max_allocated_log
integer(longint),save :: local_max_allocated_int
integer(longint),save :: local_max_allocated_real
integer(longint),save :: local_allocated_log
integer(longint),save :: local_allocated_int
integer(longint),save :: local_allocated_real

integer(longint),save :: thread_max_allocated_log
integer(longint),save :: thread_max_allocated_int
integer(longint),save :: thread_max_allocated_real
integer(longint),save :: thread_allocated_log
integer(longint),save :: thread_allocated_int
integer(longint),save :: thread_allocated_real
!$OMP THREADPRIVATE &
!$OMP (thread_max_allocated_log, &
!$OMP  thread_max_allocated_int, &
!$OMP  thread_max_allocated_real,&
!$OMP  thread_allocated_log, &
!$OMP  thread_allocated_int, &
!$OMP  thread_allocated_real)

integer(longint),parameter :: bytes_log  = int(log_bytesize ,kind=longint)
integer(longint),parameter :: bytes_int  = int(int_bytesize ,kind=longint)
integer(longint),parameter :: bytes_real = int(real_bytesize,kind=longint)

integer(longint),parameter :: kilobyte = ishft(ibset(0_longint,0),10)
integer(longint),parameter :: megabyte = ishft(kilobyte,10)
integer(longint),parameter :: gigabyte = ishft(megabyte,10)

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine mem_report
implicit none

write(*,*)

write(*,'(1x,a)',advance='no') 'Max. memory allocated for LOGICAL:  '
call mem_report_number(max_allocated_log)
write(*,'(1x,a)',advance='no') 'Max. memory allocated for INTEGER:  '
call mem_report_number(max_allocated_int)
write(*,'(1x,a)',advance='no') 'Max. memory allocated for REAL:     '
call mem_report_number(max_allocated_real)

if(allocated_log/=0_longint) then
   write(*,'(1x,a)',advance='no') 'Memory still allocated for LOGICAL: '
   call mem_report_number(allocated_log)
endif
if(allocated_int/=0_longint) then
   write(*,'(1x,a)',advance='no') 'Memory still allocated for INTEGER: '
   call mem_report_number(allocated_int)
endif
if(allocated_real/=0_longint) then
   write(*,'(1x,a)',advance='no') 'Memory still allocated for REAL:    '
   call mem_report_number(allocated_real)
endif

end subroutine mem_report

subroutine mem_report_number(number)
implicit none
integer(longint),intent(in) :: number
character(9) :: string

if(number<1000_longint) then
   write(*,'(i9,a)') number,'  b'
elseif(number<1000000_longint) then
   call prepare_string(kilobyte)
   write(*,'(a9,a)') string,' kb'
elseif(number<1000000000_longint) then
   call prepare_string(megabyte)
   write(*,'(a9,a)') string,' Mb'
else
   call prepare_string(gigabyte)
   write(*,'(a9,a)') string,' Gb'
endif

contains

  subroutine prepare_string(base)
  implicit none
  integer(longint),intent(in) :: base

  write(string,'(i5,a1,i3.3)') number/base,'.',&
       (mod(number,base)*1000_longint + base/2_longint)/base

  end subroutine prepare_string

end subroutine mem_report_number

subroutine set_parallel_memory
implicit none

parallel_memory = .true.

local_max_allocated_log  = 0_longint
local_max_allocated_int  = 0_longint
local_max_allocated_real = 0_longint
local_allocated_log  = 0_longint
local_allocated_int  = 0_longint
local_allocated_real = 0_longint

end subroutine set_parallel_memory

subroutine clear_parallel_memory
implicit none

parallel_memory = .false.

max_allocated_log  = &
     max(max_allocated_log ,allocated_log  + local_max_allocated_log)
max_allocated_int  = &
     max(max_allocated_int ,allocated_int  + local_max_allocated_int)
max_allocated_real = &
     max(max_allocated_real,allocated_real + local_max_allocated_real)
allocated_log  = allocated_log  + local_allocated_log
if(allocated_log<0_longint) then
   write(*,*) 'Negative amount of allocated memory?!'
   stop
endif
allocated_int  = allocated_int  + local_allocated_int
if(allocated_int<0_longint) then
   write(*,*) 'Negative amount of allocated memory?!'
   stop
endif
allocated_real = allocated_real + local_allocated_real
if(allocated_real<0_longint) then
   write(*,*) 'Negative amount of allocated memory?!'
   stop
endif

end subroutine clear_parallel_memory

subroutine init_thread_memory
implicit none

thread_max_allocated_log  = 0_longint
thread_max_allocated_int  = 0_longint
thread_max_allocated_real = 0_longint
thread_allocated_log  = 0_longint
thread_allocated_int  = 0_longint
thread_allocated_real = 0_longint

end subroutine init_thread_memory

subroutine collect_thread_memory
implicit none

local_max_allocated_log  = local_max_allocated_log  + thread_max_allocated_log
local_max_allocated_int  = local_max_allocated_int  + thread_max_allocated_int
local_max_allocated_real = local_max_allocated_real + thread_max_allocated_real
local_allocated_log  = local_allocated_log  + thread_allocated_log
local_allocated_int  = local_allocated_int  + thread_allocated_int
local_allocated_real = local_allocated_real + thread_allocated_real

end subroutine collect_thread_memory

subroutine mem_modify_log(xsize,add)
implicit none
integer(longint),intent(in) :: xsize
logical,intent(in) :: add

if(parallel_memory) then

   if(add) then
      thread_allocated_log = thread_allocated_log + xsize
      thread_max_allocated_log = &
           max(thread_max_allocated_log,thread_allocated_log)
   else
      thread_allocated_log = thread_allocated_log - xsize
   endif

else

   if(add) then
      allocated_log = allocated_log + xsize
      max_allocated_log = max(max_allocated_log,allocated_log)
   else
      allocated_log = allocated_log - xsize
      if(allocated_log<0_longint) then
         write(*,*) 'Negative amount of allocated memory?!'
         stop
      endif
   endif

endif

end subroutine mem_modify_log

subroutine mem_modify_int(xsize,add)
implicit none
integer(longint),intent(in) :: xsize
logical,intent(in) :: add

if(parallel_memory) then

   if(add) then
      thread_allocated_int = thread_allocated_int + xsize
      thread_max_allocated_int = &
           max(thread_max_allocated_int,thread_allocated_int)
   else
      thread_allocated_int = thread_allocated_int - xsize
   endif

else

   if(add) then
      allocated_int = allocated_int + xsize
      max_allocated_int = max(max_allocated_int,allocated_int)
   else
      allocated_int = allocated_int - xsize
      if(allocated_int<0_longint) then
         write(*,*) 'Negative amount of allocated memory?!'
         stop
      endif
   endif

endif

end subroutine mem_modify_int

subroutine mem_modify_real(xsize,add)
implicit none
integer(longint),intent(in) :: xsize
logical,intent(in) :: add

if(parallel_memory) then

   if(add) then
      thread_allocated_real = thread_allocated_real + xsize
      thread_max_allocated_real = &
           max(thread_max_allocated_real,thread_allocated_real)
   else
      thread_allocated_real = thread_allocated_real - xsize
   endif

else

   if(add) then
      allocated_real = allocated_real + xsize
      max_allocated_real = max(max_allocated_real,allocated_real)
   else
      allocated_real = allocated_real - xsize
      if(allocated_real<0_longint) then
         write(*,*) 'Negative amount of allocated memory?!'
         stop
      endif
   endif

endif

end subroutine mem_modify_real

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine alloc_log1(x,n1)
implicit none
logical,allocatable :: x(:)
integer,intent(in) :: n1
integer :: ierr
integer(longint) :: xsize

allocate(x(n1),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' log'
   call mem_report
   stop
endif

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.true.)

end subroutine alloc_log1

subroutine alloc_log2(x,n1,n2)
implicit none
logical,allocatable :: x(:,:)
integer,intent(in) :: n1,n2
integer :: ierr
integer(longint) :: xsize

allocate(x(n1,n2),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' x',n2,' log'
   call mem_report
   stop
endif

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.true.)

end subroutine alloc_log2

subroutine alloc_log1b(x,n1)
implicit none
logical,allocatable :: x(:)
integer,intent(in) :: n1(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' log'
   call mem_report
   stop
endif

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.true.)

end subroutine alloc_log1b

subroutine alloc_log2b(x,n1,n2)
implicit none
logical,allocatable :: x(:,:)
integer,intent(in) :: n1(2),n2(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2),n2(1):n2(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' x',n2(1),' :',n2(2),' log'
   call mem_report
   stop
endif

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.true.)

end subroutine alloc_log2b

subroutine dealloc_log1(x)
implicit none
logical,allocatable :: x(:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_log1

subroutine dealloc_log2(x)
implicit none
logical,allocatable :: x(:,:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_log*size(x,kind=longint)
call mem_modify_log(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_log2

!-------------------------------------------------------------------------------

subroutine alloc_int1(x,n1)
implicit none
integer,allocatable :: x(:)
integer,intent(in) :: n1
integer :: ierr
integer(longint) :: xsize

allocate(x(n1),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' int'
   call mem_report
   stop
endif

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.true.)

end subroutine alloc_int1

subroutine alloc_int2(x,n1,n2)
implicit none
integer,allocatable :: x(:,:)
integer,intent(in) :: n1,n2
integer :: ierr
integer(longint) :: xsize

allocate(x(n1,n2),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' x',n2,' int'
   call mem_report
   stop
endif

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.true.)

end subroutine alloc_int2

subroutine alloc_int1b(x,n1)
implicit none
integer,allocatable :: x(:)
integer,intent(in) :: n1(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' int'
   call mem_report
   stop
endif

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.true.)

end subroutine alloc_int1b

subroutine alloc_int2b(x,n1,n2)
implicit none
integer,allocatable :: x(:,:)
integer,intent(in) :: n1(2),n2(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2),n2(1):n2(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' x',n2(1),' :',n2(2),' int'
   call mem_report
   stop
endif

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.true.)

end subroutine alloc_int2b

subroutine dealloc_int1(x)
implicit none
integer,allocatable :: x(:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_int1

subroutine dealloc_int2(x)
implicit none
integer,allocatable :: x(:,:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_int*size(x,kind=longint)
call mem_modify_int(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_int2

!-------------------------------------------------------------------------------

subroutine alloc_real1(x,n1)
implicit none
real(prec),allocatable :: x(:)
integer,intent(in) :: n1
integer :: ierr
integer(longint) :: xsize

allocate(x(n1),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real1

subroutine alloc_real2(x,n1,n2)
implicit none
real(prec),allocatable :: x(:,:)
integer,intent(in) :: n1,n2
integer :: ierr
integer(longint) :: xsize

allocate(x(n1,n2),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' x',n2,' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real2

subroutine alloc_real3(x,n1,n2,n3)
implicit none
real(prec),allocatable :: x(:,:,:)
integer,intent(in) :: n1,n2,n3
integer :: ierr
integer(longint) :: xsize

allocate(x(n1,n2,n3),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' x',n2,' x',n3,' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real3

subroutine alloc_real4(x,n1,n2,n3,n4)
implicit none
real(prec),allocatable :: x(:,:,:,:)
integer,intent(in) :: n1,n2,n3,n4
integer :: ierr
integer(longint) :: xsize

allocate(x(n1,n2,n3,n4),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1,' x',n2,' x',n3,' x',n4,' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real4

subroutine alloc_real1b(x,n1)
implicit none
real(prec),allocatable :: x(:)
integer,intent(in) :: n1(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real1b

subroutine alloc_real2b(x,n1,n2)
implicit none
real(prec),allocatable :: x(:,:)
integer,intent(in) :: n1(2),n2(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2),n2(1):n2(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) n1(1),' :',n1(2),' x',n2(1),' :',n2(2),' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real2b

subroutine alloc_real3b(x,n1,n2,n3)
implicit none
real(prec),allocatable :: x(:,:,:)
integer,intent(in) :: n1(2),n2(2),n3(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2),n2(1):n2(2),n3(1):n3(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) &
        n1(1),' :',n1(2),' x',n2(1),' :',n2(2),' x', &
        n3(1),' :',n3(2),' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real3b

subroutine alloc_real4b(x,n1,n2,n3,n4)
implicit none
real(prec),allocatable :: x(:,:,:,:)
integer,intent(in) :: n1(2),n2(2),n3(2),n4(2)
integer :: ierr
integer(longint) :: xsize

allocate(x(n1(1):n1(2),n2(1):n2(2),n3(1):n3(2),n4(1):n4(2)),stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while allocating memory!'
   write(*,*) &
        n1(1),' :',n1(2),' x',n2(1),' :',n2(2),' x', &
        n3(1),' :',n3(2),' x',n4(1),' :',n4(2),' real'
   call mem_report
   stop
endif

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.true.)

end subroutine alloc_real4b

subroutine dealloc_real1(x)
implicit none
real(prec),allocatable :: x(:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_real1

subroutine dealloc_real2(x)
implicit none
real(prec),allocatable :: x(:,:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_real2

subroutine dealloc_real3(x)
implicit none
real(prec),allocatable :: x(:,:,:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_real3

subroutine dealloc_real4(x)
implicit none
real(prec),allocatable :: x(:,:,:,:)
integer :: ierr
integer(longint) :: xsize

xsize = bytes_real*size(x,kind=longint)
call mem_modify_real(xsize,.false.)

if(.not.allocated(x)) then
   write(*,*) 'Trying to deallocate nonallocated memory!'
   stop
endif
deallocate(x,stat=ierr)
if(ierr/=0) then
   write(*,*) 'Problem while deallocating memory!'
   stop
endif

end subroutine dealloc_real4

end module memory
