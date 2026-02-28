module string
      use arithmetic
      implicit none

      interface str
            module procedure :: str_i32
            module procedure :: str_i64
            module procedure :: str_f64
      end interface str

contains

      pure function str_i32(i)
            !
            ! Convert integer to string. The result does not
            ! contain any blanks.
            !
            character(len=:), allocatable :: str_i32
            integer(I32), intent(in) :: i

            character(len=I32_MAXWIDTH) :: t
            integer :: l

            write(t, "(I0)") i
            l = len_trim(t)
            allocate(character(l) :: str_i32)
            str_i32(:) = t(1:l)
      end function str_i32

      
      pure function str_i64(i)
            !
            ! Convert integer to string. The result does not
            ! contain any blanks.
            !
            character(:), allocatable :: str_i64
            integer(I64), intent(in) :: i

            character(len=I64_MAXWIDTH) :: t
            integer :: l

            write(t, "(I0)") i
            l = len_trim(t)
            allocate(character(l) :: str_i64)
            str_i64(:) = t(1:l)
      end function str_i64


      pure function str_f64(f, d)
            !
            ! Convert floating point number to string.
            ! The result does not contain any blanks.
            !
            character(:), allocatable     :: str_f64
            real(F64), intent(in)         :: f
            integer, optional, intent(in) :: d

            character(F64_ES_W) :: t
            character(:), allocatable :: fmt

            if (present(d)) then
                  !
                  ! (ES{F64_ES_W}.dE{F64_ES_E})
                  !
                  fmt = "(ES" // str(F64_ES_W) // "." // &
                        str(min(F64_ES_D, d)) // "E" // str(F64_ES_E) // ")"
            else
                  !
                  ! (ES{F64_ES_W}.{F64_ES_D}E{F64_ES_E})
                  !
                  fmt = "(ES" // str(F64_ES_W) // "." // &
                        str(F64_ES_D) // "E" // str(F64_ES_E) // ")"
            end if
            write(t, fmt) f
            str_f64 = trim(adjustl(t))
      end function str_f64


      pure function lfield(s, w)
            !
            ! Create a left-aligned character string of width W.
            ! Left-aligned cells are useful for printing table headers.
            !
            character(:), allocatable :: lfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k

            allocate(character(w) :: lfield)
            k = min(w, len_trim(s))
            lfield(1:w) = " "
            lfield(1:k) = s(1:k)
      end function lfield


      pure function rfield(s, w)
            !
            ! Create a rightt-aligned character string of width W.
            !
            character(:), allocatable :: rfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k, ls

            ls = len_trim(s)
            if (ls >= w) then
                  rfield = s(1:ls)
            else
                  allocate(character(w) :: rfield)
                  k = min(w, len_trim(s))
                  rfield(1:w) = " "
                  rfield(w-k+1:w) = s(1:k)
            end if
      end function rfield


      pure function cfield(s, w)
            !
            ! Create a centered character string of width W.
            ! Centered cells are useful for printing table headers.
            !
            character(:), allocatable :: cfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k0, k1, l, ls

            ls = len_trim(s)
            if (ls >= w) then
                  cfield = s(1:ls)
            else
                  allocate(character(w) :: cfield)
                  l = min(w, len_trim(s))
                  k0 = (w - l) / 2 + 1
                  k1 = k0 + l - 1
                  cfield(1:w) = " "
                  cfield(k0:k1) = s(1:l)
            end if
      end function cfield
end module string
