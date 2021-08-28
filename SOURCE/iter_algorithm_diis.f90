module class_IterAlgorithmDIIS

    use class_IterAlgorithm
    implicit none

    type, public, extends(IterAlgorithm) :: IterAlgorithmDIIS
        integer :: DIISN = 6
        contains
            procedure :: iterate => iterate
            procedure :: getName => getName
            procedure :: getParamsString => getParamsString
            procedure :: toString => toString
    end type IterAlgorithmDIIS


    contains


    subroutine iterate(this, CTilde, NDimX, NCholesky, Lambda, APlusTilde, A2, iStats)

        use diis
        use class_IterStats
        implicit none
        class(IterAlgorithmDIIS) :: this
        type(DIISData) :: DIISBlock

        double precision, intent(out) :: CTilde(NDimX*NCholesky)
        integer,intent(in) :: NDimX, NCholesky
        double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)
        class(IterStats), intent(inout) :: iStats

        integer :: N
        double precision :: A2CTilde(NDimX*NCholesky), CTilde_prev(NDimX*NCholesky)
        double precision :: norm

        !  CTilde = C0
        ! Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,Lambda,NDimX,&
        !             APlusTilde,NDimX,0.0d0,CTilde,NDimX)

        call init_DIIS(DIISBlock,NDimX*NCholesky,NDimX*NCholesky,this%DIISN)
        CTilde_prev = 0.0d0
        N = 1
        do

            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,&
                    CTilde,NDimX,0.0d0,A2CTilde,NDimX)
            A2CTilde=APlusTilde-A2CTilde

            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,Lambda,NDimX,&
                    A2CTilde,NDimX,0.0d0,CTilde,NDimX)

            if(N > 2) then
                call use_DIIS(DIISBlock, CTilde, CTilde - CTilde_prev)
            endif

            norm = norm2(CTilde - CTilde_prev)
            ! print '(a, i3, a, e)', "norm (", N, ") = ", norm

            if ((norm < this%Threshold) .or. (N .ge. this%maxIterations .and. this%maxIterations .ge. 0)) then
                call iStats%addN(N)
                exit
            endif

            CTilde_prev = CTilde
            N = N + 1

        enddo
        call free_DIIS(DIISBlock)
        

    end subroutine iterate


    function getName(this) result(res)

        implicit none
        class(IterAlgorithmDIIS) :: this
        character(:), allocatable :: res
        res = "DIIS"

    end function getName


    function getParamsString(this) result(res)

        implicit none
        class(IterAlgorithmDIIS) :: this
        character(:), allocatable :: res
        character(200) :: buffer
        write(buffer, "(a, i0, a, e6.1e1, a, i0, a)") &
            "{DIISN: ", this%DIISN, ", Threshold: ", this%Threshold, ", maxIterations: ", this%maxIterations, "}"
        res = trim(buffer)

    end function getParamsString


    subroutine toString(this, fileName)

        implicit none
        class(IterAlgorithmDIIS) :: this
        character(50), intent(out) :: fileName
        write(fileName, "(a, i0, a, e6.1e1)") 'diis_n', this%DIISN, '_t', this%Threshold

    end subroutine toString


end module class_IterAlgorithmDIIS