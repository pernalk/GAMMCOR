    module class_IterDIIS

    use class_IterAlgorithm
    implicit none

    type, public, extends(IterAlgorithm) :: IterDIIS
        integer :: NGrid = 35, DIISN = 6, maxIterations = 0
        double precision :: Threshold = 1d-5
        contains
            procedure :: iterate => iterate
            procedure :: createFileName => createFileName
    end type IterDIIS


    contains


    subroutine iterate(this, CTilde, N, NDimX, NCholesky, Lambda, APlusTilde, A2)

        use diis
        implicit none
        class(IterDIIS) :: this
        type(DIISData) :: DIISBlock

        double precision, intent(out) :: CTilde(NDimX*NCholesky)
        integer,intent(in) :: NDimX, NCholesky
        double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)

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
            if (norm < this%Threshold) exit

            CTilde_prev = CTilde
            N = N + 1

        enddo
        call free_DIIS(DIISBlock)

    end subroutine iterate


    subroutine createFileName(this, fileName)

        implicit none
        class(IterDIIS) :: this
        character(50), intent(out) :: fileName
        write(fileName, "(3a, i0, a, i0, a, e6.1e1)") &
            'is_diis_', 'diag', '_g', this%NGrid, '_n', this%DIISN, '_t', this%Threshold

    end subroutine createFileName


end module class_IterDIIS