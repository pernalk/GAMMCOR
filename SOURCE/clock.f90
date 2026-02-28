module clockMM
      use arithmetic
      use omp_lib

      implicit none

      type tclock
            real(F64) :: cpu_start
            real(F64) :: wall_start
      end type tclock

contains

      subroutine clock_start(t)
            type(tclock), intent(out) :: t
            
            call cpu_time(t%cpu_start)
            t%wall_start = omp_get_wtime()
      end subroutine clock_start

      
      function clock_readcpu(t)
            !
            ! Read CPU time in seconds
            !
            real(F64)                :: clock_readcpu
            type(tclock), intent(in) :: t
            
            real(F64)                :: cpu_end
            call cpu_time(cpu_end)
            clock_readcpu = cpu_end - t%cpu_start
      end function clock_readcpu


      function clock_readwall(t)
            !
            ! Read wallclock time in seconds
            !
            real(F64)                :: clock_readwall
            type(tclock), intent(in) :: t

            real(F64) :: wall_end
            wall_end = omp_get_wtime()
            clock_readwall = wall_end - t%wall_start
      end function clock_readwall
end module clockMM
