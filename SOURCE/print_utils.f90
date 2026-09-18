module print_utils

    use arithmetic, only : F64
    use display, only : STDOUNIT
    use clock
    use types, only : PRT_MSG
    implicit none
    private
    public :: print_header, print_section, print_info, print_status, print_kv
    public :: print_energy, print_energy_component, print_block_contribution
    public :: print_energy1
    public :: print_memory
    public :: prij, pri3, prix, pr2ix, pr4ix
    public :: mmsg, xmsg, ximsg, tmsg
    public :: print_deb, print_debi
    public :: MSG_DEB, MSG_PRIORITY_THR

    integer, parameter :: line_width  = 70
    integer, parameter :: label_width = 40
    integer, parameter :: block_label_width = 10
    integer, parameter :: deb_label_width = 25

    !
    ! Message priorities: a message is printed only when its priority
    ! is at least MSG_PRIORITY_THR (timings are gated on MSG_TIME alone).
    !
    integer, parameter, private :: MSG_TIME = 0
    integer, parameter, private :: MSG_VERB = 10
    integer, parameter, private :: MSG_NOR  = 50
    integer, parameter, private :: MSG_ERR  = 100

    integer, parameter :: MSG_DEB = 0

    integer, save :: MSG_PRIORITY_THR = 1

    interface print_info
        module procedure print_info_str
        module procedure print_info_int
        module procedure print_info_int_arr
        module procedure print_info_dble
        module procedure print_info_dble2
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
          print*, ''
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

    subroutine print_info_dble2(label, dval, dval2)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: dval, dval2
        write(*,'(2X,A,": ",F18.8,4X,F18.8)') pad_label(label), dval, dval2
    end subroutine print_info_dble2

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

    subroutine print_energy(label, energy_eh, energy_ev)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: energy_eh, energy_ev

        write(*,'(2X,A,":",2X,F18.10,1X,A,4X,"(",F18.10,1X,A,")")') &
            pad_label(label), energy_eh, 'Eh', energy_ev, 'eV'
    end subroutine print_energy

    subroutine print_energy1(label, energy_eh)
          character(len=*), intent(in) :: label
          double precision, intent(in) :: energy_eh
          
          write(*,'(2X,A,":",2X,F18.10,1X,A,4X,"(",F18.10,1X,A,")")') &
                pad_label(label), energy_eh, 'Eh'
    end subroutine print_energy1

    subroutine print_energy_component(label, total, singlet, triplet, triplet_aaaa)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: total, singlet, triplet, triplet_aaaa
        character(len=12) :: lbl

        lbl = adjustl(label)
        write(*,'(2X,A12,4(1X,F14.9))') lbl, total, singlet, triplet, triplet_aaaa
    end subroutine print_energy_component

    subroutine print_block_contribution(block_pair, energy)
        character(len=*), intent(in) :: block_pair
        double precision, intent(in) :: energy
        character(len=block_label_width) :: lbl

        lbl = adjustl(block_pair)
        write(*,'(2X,A,2X,F18.10)') lbl, energy
    end subroutine print_block_contribution

    subroutine print_memory(label, memory_gb)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: memory_gb

        write(*,'(2X,A,": ",ES12.4,1X,A)') pad_label(label), memory_gb, 'GB'
    end subroutine print_memory

    !
    ! ---------------------------------------------------------------
    ! Compact debug printers
    ! ---------------------------------------------------------------
    !
    subroutine prij(st, i, j)
        character(len=*), intent(in) :: st
        integer, intent(in) :: i, j
        write(*,'(2X,A,": ",2I8)') pad_label(st), i, j
    end subroutine prij

    subroutine pri3(st, i, j, k)
        character(len=*), intent(in) :: st
        integer, intent(in) :: i, j, k
        write(*,'(2X,A,": ",3I8)') pad_label(st), i, j, k
    end subroutine pri3

    subroutine prix(st, i, a)
        character(len=*), intent(in) :: st
        integer, intent(in) :: i
        double precision, intent(in) :: a
        write(*,'(2X,A,": ",I8,F25.15)') pad_label(st), i, a
    end subroutine prix

    subroutine pr2ix(st, i, j, a)
        character(len=*), intent(in) :: st
        integer, intent(in) :: i, j
        double precision, intent(in) :: a
        write(*,'(2X,A,": ",2I8,F20.15)') pad_label(st), i, j, a
    end subroutine pr2ix

    subroutine pr4ix(st, i, j, k, l, a)
        character(len=*), intent(in) :: st
        integer, intent(in) :: i, j, k, l
        double precision, intent(in) :: a
        write(*,'(2X,A,": ",4I8,F20.15)') pad_label(st), i, j, k, l, a
    end subroutine pr4ix

    !
    ! ---------------------------------------------------------------
    ! Priority-gated messages
    ! ---------------------------------------------------------------
    !
    subroutine mmsg(s, priority)
        character(len=*), intent(in) :: s
        integer, optional, intent(in) :: priority
        integer :: p

        p = msg_priority(priority)
        if (p >= MSG_PRIORITY_THR) then
            write(STDOUNIT,'(2X,A)') pad_label(s)
            flush(STDOUNIT)
        end if
    end subroutine mmsg

    subroutine xmsg(s, t, priority)
        character(len=*), intent(in) :: s
        real(F64), intent(in) :: t
        integer, optional, intent(in) :: priority
        integer :: p

        p = msg_priority(priority)
        if (p >= MSG_PRIORITY_THR) then
            write(STDOUNIT,'(2X,A,": ",F30.15)') pad_label(s), t
            flush(STDOUNIT)
        end if
    end subroutine xmsg

    subroutine ximsg(s, t, priority)
        character(len=*), intent(in) :: s
        integer, intent(in) :: t
        integer, optional, intent(in) :: priority
        integer :: p

        p = msg_priority(priority)
        if (p >= MSG_PRIORITY_THR) then
            write(STDOUNIT,'(2X,A,": ",I15)') pad_label(s), t
            flush(STDOUNIT)
        end if
    end subroutine ximsg

    subroutine tmsg(s, t, priority)
        character(len=*), intent(in) :: s
        type(tclock), intent(in) :: t
        integer, optional, intent(in) :: priority
        real(F64) :: d
        integer :: p

        p = msg_priority(priority)
        if (p == MSG_TIME) then
            d = clock_readwall(t)
            write(STDOUNIT,'(2X,A,": ",F20.3,1X,A)') pad_label(s), d, '[sec]'
            flush(STDOUNIT)
        end if
    end subroutine tmsg

    function msg_priority(priority) result(p)
        integer, optional, intent(in) :: priority
        integer :: p
        if (present(priority)) then
            p = priority
        else
            p = MSG_NOR
        end if
    end function msg_priority

    !
    ! ---------------------------------------------------------------
    ! Debug prints, gated on PRT_MSG (set in types.f90: 1 with -DDEBUG,
    ! 0 otherwise). With PRT_MSG == 0 both routines return immediately.
    ! The label column is kept even when msg is absent, so that values
    ! printed with and without a label still line up.
    ! ---------------------------------------------------------------
    !
    subroutine print_deb(val, msg)
        double precision, intent(in) :: val
        character(len=*), optional, intent(in) :: msg
        character(len=deb_label_width) :: lbl

        if (PRT_MSG /= 1) return

        lbl = ''
        if (present(msg)) lbl = msg
        write(STDOUNIT,'(2X,A,1X,F22.15)') lbl, val
        flush(STDOUNIT)
    end subroutine print_deb

    subroutine print_debi(val, msg)
        integer, intent(in) :: val
        character(len=*), optional, intent(in) :: msg
        character(len=deb_label_width) :: lbl

        if (PRT_MSG /= 1) return

        lbl = ''
        if (present(msg)) lbl = msg
        write(STDOUNIT,'(2X,A,1X,I22)') lbl, val
        flush(STDOUNIT)
    end subroutine print_debi

end module print_utils
