module class_IterAlgorithm

    implicit none
  
    type, public, abstract :: IterAlgorithm
        integer :: maxIterations = -1
        double precision :: Threshold = 1d-5
        contains
            procedure(iterateInterface), deferred :: iterate
            procedure(algGetName), deferred :: getName
            procedure(algGetParamsString), deferred :: getParamsString
            procedure(algToStringInterface), deferred :: toString
    end type IterAlgorithm

    interface
        subroutine iterateInterface(this, CTilde, NDimX, NCholesky, Lambda, APlusTilde, A2, iStats)
            use class_IterStats
            import IterAlgorithm
            implicit none
            class(IterAlgorithm) :: this
            integer,intent(in) :: NDimX, NCholesky
            double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)
            double precision, intent(out) :: CTilde(NDimX*NCholesky)
            class(IterStats), intent(inout) :: iStats
            double precision :: C0(NDimX*NCholesky), A2CTilde(NDimX*NCholesky), CTilde_prev(NDimX*NCholesky)
            double precision :: norm
        end subroutine iterateInterface
    end interface

    interface
        function algGetName(this) result(res)
            import IterAlgorithm
            class(IterAlgorithm) :: this
            character(:), allocatable :: res
        end function algGetName
    end interface

    interface
        function algGetParamsString(this) result(res)
            import IterAlgorithm
            class(IterAlgorithm) :: this
            character(:), allocatable :: res
        end function algGetParamsString
    end interface

    interface
        subroutine algToStringInterface(this, fileName)
            import IterAlgorithm
            class(IterAlgorithm) :: this
            character(50), intent(out) :: fileName
        end subroutine algToStringInterface
    end interface


end module class_IterAlgorithm