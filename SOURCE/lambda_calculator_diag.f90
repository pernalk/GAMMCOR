module class_LambdaCalculatorDiag
  
    use class_LambdaCalculator
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorDiag
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
    end type LambdaCalculatorDiag


    contains


    subroutine calculateLambda(this, Lambda, OmI, NDimX, A0)

        implicit none
        class(LambdaCalculatorDiag) :: this
        double precision, intent(out) :: Lambda(NDimX*NDimX)
        integer, intent(in) :: NDimX
        double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        integer :: i, j

        Lambda=0.D0
        Do i=1,NDimX
            j = (i-1) * NDimX + i
            Lambda(j) = 1 / ( A0(j) + OmI**2 )
        EndDo

    end subroutine calculateLambda


    subroutine calculateInitialA(this, A0, A2, NDimX)

        implicit none
        class(LambdaCalculatorDiag) :: this
        double precision, intent(out) :: A0(NDimX*NDimX)
        double precision, intent(inout) :: A2(NDimX*NDimX)
        integer, intent(in) :: NDimX
        integer :: i

        A0=0.D0   
        Do i=1,NDimX
           A0((i-1)*NDimX+i)=A2((i-1)*NDimX+i)
        EndDo
        A2=A2-A0

    end subroutine calculateInitialA


    function getName(this) result(res)

        implicit none
        class(LambdaCalculatorDiag) :: this
        character(:), allocatable :: res
        res = "Diag"

    end function getName


    subroutine toString(this, string)

        implicit none
        class(LambdaCalculatorDiag) :: this
        character(50), intent(out) :: string
        string = "diag"

    end subroutine toString


end module class_LambdaCalculatorDiag