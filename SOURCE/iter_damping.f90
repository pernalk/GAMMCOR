    module class_Iterdamping

    use class_IterAlgorithm
    implicit none

    type, public, extends(IterAlgorithm) :: Iterdamping
        integer :: NGrid = 35, maxIterations = 0
        double precision :: Threshold = 1d-5, Xmix = 0.2
        contains
            procedure :: iterate => iterate
            procedure :: createFileName => createFileName
    end type Iterdamping


    contains


    subroutine iterate(this, CTilde, N, NDimX, NCholesky, Lambda, APlusTilde, A2)

        implicit none
        class(IterDamping) :: this

        double precision, intent(out) :: CTilde(NDimX*NCholesky)
        integer,intent(in) :: NDimX, NCholesky
        double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)
        
        double precision :: CTilde_prev(NDimX*NCholesky)
        integer :: N
        double precision :: A2CTilde(NDimX*NCholesky)
        double precision :: norm

        ! !     C0=C(0)
        ! Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,Lambda,NDimX,&
        !         APlusTilde,NDimX,0.0d0,C0,NDimX)
        !     C1=C(1)
        Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,&
                CTilde,NDimX,0.0d0,A2CTilde,NDimX)
        A2CTilde=APlusTilde-A2CTilde
        Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,Lambda,NDimX,&
                A2CTilde,NDimX,0.0d0,CTilde,NDimX)
        N = 2
        Do
            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0,A2,NDimX,&
                CTilde,NDimX,0.0d0,A2CTilde,NDimX)
            A2CTilde=APlusTilde-A2CTilde
        ! damping is needed when active orbitals present: CMAT(n) = (1-XMix)*CMAT(n) + XMix*CMAT(n-1)
            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.0d0-this%XMix,Lambda,NDimX,&
                A2CTilde,NDimX,this%XMix,CTilde,NDimX)   

            norm = norm2(CTilde - CTilde_prev)
            ! print '(a, i3, a, e)', "norm (", N, ") = ", norm
            if (norm < this%Threshold) exit

            CTilde_prev = CTilde
            N = N + 1
        EndDo

    end subroutine iterate


    subroutine createFileName(this, fileName)

        implicit none
        class(IterDamping) :: this
        character(50), intent(out) :: fileName
        write(fileName, "(3a, i0, a, f4.2, a, e6.1e1)") &
            'is_damping_', 'diag', '_g', this%NGrid, '_x', this%XMix, '_t', this%Threshold

    end subroutine createFileName


end module class_Iterdamping