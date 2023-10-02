    module class_IterAlgorithmDamping

    use class_IterAlgorithm
    implicit none

    type, public, extends(IterAlgorithm) :: IterAlgorithmDamping
        double precision :: Xmix = 0.2
        contains
            procedure :: iterate => iterate
            procedure :: getName => getName
            procedure :: getParamsString => getParamsString
            procedure :: toString => toString
    end type IterAlgorithmDamping


    contains


        subroutine iterate(this, CTilde, NDimX, NCholesky, Lambda, APlusTilde, A2, iStats)

            use class_IterStats
            implicit none
            class(IterAlgorithmDamping) :: this

            double precision, intent(out) :: CTilde(NDimX*NCholesky)
            integer,intent(in) :: NDimX, NCholesky
            double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)
            class(IterStats), intent(inout) :: iStats
            
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
                if ((norm < this%Threshold) .or. (N .ge. this%maxIterations .and. this%maxIterations .ge. 0)) then
                    call iStats%addN(N)
                    exit
                endif

                CTilde_prev = CTilde
                N = N + 1
            EndDo

            call iStats%addN(N)

        end subroutine iterate


        function getName(this) result(res)

            implicit none
            class(IterAlgorithmDamping) :: this
            character(:), allocatable :: res
            res = "Damping"

        end function getName


        function getParamsString(this) result(res)

            implicit none
            class(IterAlgorithmDamping) :: this
            character(:), allocatable :: res
            character(200) :: buffer
            write(buffer, "(a, f4.2, a, e7.1e1, a, i0, a)") &
                "{XMix: ", this%XMix, ", Threshold: ", this%Threshold, ", maxIterations: ", this%maxIterations, "}"
            res = trim(buffer)

        end function getParamsString


        subroutine toString(this, fileName)

            implicit none
            class(IterAlgorithmDamping) :: this
            character(50), intent(out) :: fileName
            write(fileName, "(a, f4.2, a, e7.1e1)") 'damping_x', this%XMix, '_t', this%Threshold

        end subroutine toString


end module class_IterAlgorithmDamping
