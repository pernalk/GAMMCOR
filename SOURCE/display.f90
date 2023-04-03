module display
      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

      implicit none
      !
      ! Standard output
      !
      integer, parameter :: STDOUNIT = OUTPUT_UNIT

contains

      subroutine msg(s)
            !
            ! Output character string to the standard output
            !
            character(len=*), intent(in)  :: s

            write(STDOUNIT, "(1X,A)") s
            flush(STDOUNIT)
      end subroutine msg

      
      subroutine midrule(width)
            integer, intent(in) :: width

            write(STDOUNIT, "(1X,A)") repeat("-", width)
            flush(STDOUNIT)
      end subroutine midrule


      subroutine blankline()
            write(STDOUNIT, "(1X)")
            flush(STDOUNIT)
      end subroutine blankline
end module display
