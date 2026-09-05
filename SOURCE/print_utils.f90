module print_utils
    implicit none
    private
    public :: print_header, print_section, print_info, print_status, print_kv

    integer, parameter :: line_width  = 70
    integer, parameter :: label_width = 32

    interface print_info
        module procedure print_info_str
        module procedure print_info_int
        module procedure print_info_int_arr
        module procedure print_info_dble
    end interface print_info

contains

    subroutine print_header(title)
        character(len=*), intent(in) :: title
        character(len=line_width) :: border
        border = repeat('=', line_width)
        write(*,'(A)') border
        write(*,'(A)') center_text(title, line_width)
        write(*,'(A)') border
    end subroutine print_header

    subroutine print_section(title)
        character(len=*), intent(in) :: title
        write(*,'(/,A)') '-- ' // trim(title) // ' ' // repeat('-', max(0, line_width - len_trim(title) - 4))
    end subroutine print_section

    subroutine print_info_str(label, msg)
        character(len=*), intent(in) :: label, msg
        write(*,'(2X,A,": ",A)') pad_label(label), trim(msg)
    end subroutine print_info_str

    subroutine print_info_int(label, ival)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ival
        write(*,'(2X,A,": ",I8)') pad_label(label), ival
    end subroutine print_info_int

    subroutine print_info_int_arr(label, ivals)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ivals(:)
        character(len=256) :: buf, item
        integer :: i
        buf = ''
        do i = 1, size(ivals)
            write(item,'(I0)') ivals(i)
            if (i == 1) then
                buf = trim(item)
            else
                buf = trim(buf) // ', ' // trim(item)
            end if
        end do
        write(*,'(2X,A,": ",A)') pad_label(label), trim(buf)
    end subroutine print_info_int_arr

    subroutine print_info_dble(label, dval)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: dval
        write(*,'(2X,A,": ",F18.8)') pad_label(label), dval
    end subroutine print_info_dble

    function pad_label(label) result(lbl)
        character(len=*), intent(in) :: label
        character(len=label_width) :: lbl
        lbl = label
    end function pad_label

    subroutine print_status(msg, ok)
        character(len=*), intent(in) :: msg
        logical, intent(in) :: ok
        character(len=8) :: tag
        tag = merge('[ OK ]  ', '[FAIL]  ', ok)
        write(*,'(2X,A,A)') tag, trim(msg)
    end subroutine print_status

    subroutine print_kv(key, val)
        character(len=*), intent(in) :: key
        real(kind=8), intent(in) :: val
        write(*,'(2X,A20,": ",F16.8)') key, val
    end subroutine print_kv

    function center_text(text, width) result(centered)
        character(len=*), intent(in) :: text
        integer, intent(in) :: width
        character(len=width) :: centered
        integer :: pad
        pad = max(0, (width - len_trim(text)) / 2)
        centered = repeat(' ', pad) // trim(text)
    end function center_text

end module print_utils
