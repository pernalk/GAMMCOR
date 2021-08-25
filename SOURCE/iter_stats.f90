    module class_IterStats

    implicit none

    type, public :: IterStats
        integer :: count = 0
        integer :: iterNoTotal = 0, iterMax = 0
        double precision :: freqMax = 0, iterAvg = 0
        character(13) :: buffer = ''
        contains
            procedure :: add => add
            procedure :: print => print
            procedure :: dump => dump
    end type IterStats


    contains


    subroutine add(this, freq, iterNo)

        implicit none
        class(IterStats) :: this
        integer :: iterNo
        double precision :: freq

        if(iterNo > this%iterMax) then
            this%iterMax = iterNo
            this%freqMax = freq
        end if

        this%count = this%count + 1
        this%iterNoTotal = this%iterNoTotal + iterNo
        this%iterAvg = float(this%iterNoTotal) / this%count

    end subroutine add


    subroutine print(this, unitOp)

        implicit none
        class(IterStats) :: this
        integer, optional :: unitOp
        integer :: unit

        unit = merge(unitOp, 6, present(unitOp))

        write (unit, *) "------------------------------------------------"
        write (this%buffer, '(i13)') this%iterNoTotal
        write (unit, '(a37, a13)') "total number of iterations:         ", adjustl(this%buffer)
        write (this%buffer, '(f4.1)') this%iterAvg
        write (unit, '(a37, a13)') "average number of iterations:       ", adjustl(this%buffer)
        write (this%buffer, '(i13)') this%iterMax
        write (unit, '(a37, a13)') "max number of iterations:           ", adjustl(this%buffer)
        write (this%buffer, '(e13.6)') this%freqMax
        write (unit, '(a37, a13)') "freq with max number of iterations: ", adjustl(this%buffer)
        write (unit, *) "------------------------------------------------"

    end subroutine print


    subroutine dump(this, fileName, ACAlpha)

        implicit none
        class(IterStats) :: this
        ! integer, intent(in) :: DIISN, NGrid
        double precision, intent(in) :: ACAlpha!, Threshold
        integer :: unit
        character(50) :: fileName

        ! write(fileName,'(a, i0, a, i0, a, e6.1e1)') 'iter_stats_', NGrid, '_', DIISN, '_', Threshold

        open (newunit=unit, file = fileName, position="append")
        write (unit, *) "ACAlpha = ", ACAlpha
        call this%print(unit)
        close(unit)

    end subroutine dump

end module class_IterStats