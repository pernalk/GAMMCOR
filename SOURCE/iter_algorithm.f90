module class_IterAlgorithm
  
    implicit none
  
    type, public, abstract :: IterAlgorithm
      contains
        procedure(iterateInterface), deferred :: iterate
        procedure(createFileNameInterface), deferred :: createFileName
    end type IterAlgorithm

    interface
        subroutine iterateInterface(this, CTilde, N, NDimX, NCholesky, Lambda, APlusTilde, A2)
            import IterAlgorithm
            class(IterAlgorithm) :: this
            integer,intent(in) :: NDimX, NCholesky
            double precision, intent(in) :: Lambda(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky)
            double precision, intent(out) :: CTilde(NDimX*NCholesky)
            integer :: N
            double precision :: C0(NDimX*NCholesky), A2CTilde(NDimX*NCholesky), CTilde_prev(NDimX*NCholesky)
            double precision :: norm
        end subroutine iterateInterface
    end interface

    interface
        subroutine createFileNameInterface(this, fileName)
            import IterAlgorithm
            class(IterAlgorithm) :: this
            character(50), intent(out) :: fileName
        end subroutine createFileNameInterface
    end interface


end module class_IterAlgorithm