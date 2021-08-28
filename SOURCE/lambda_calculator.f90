module class_LambdaCalculator

    implicit none
  
    type, public, abstract :: LambdaCalculator
        contains
            procedure(calculateInitialAInterface), deferred :: calculateInitialA
            procedure(calculateLambdaInterface), deferred :: calculateLambda
            procedure(lamGetName), deferred ::getName
            procedure(lamToStringInterface), deferred :: toString
    end type LambdaCalculator

    interface
        subroutine calculateLambdaInterface(this, Lambda, OmI, NDimX, A0)
            import LambdaCalculator
            implicit none
            class(LambdaCalculator) :: this
            double precision, intent(out) :: Lambda(NDimX*NDimX)
            integer, intent(in) :: NDimX
            double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        end subroutine calculateLambdaInterface
    end interface

    interface
        subroutine calculateInitialAInterface(this, A0, A2, NDimX)
            import LambdaCalculator
            implicit none
            class(LambdaCalculator) :: this
            double precision, intent(out) :: A0(NDimX*NDimX)
            double precision, intent(inout) :: A2(NDimX*NDimX)
            integer, intent(in) :: NDimX
        end subroutine calculateInitialAInterface
    end interface

    interface
        function lamGetName(this) result(res)
            import LambdaCalculator
            class(LambdaCalculator) :: this
            character(:), allocatable :: res
        end function lamGetName
    end interface

    interface
        subroutine lamToStringInterface(this, string)
            import LambdaCalculator
            implicit none
            class(LambdaCalculator) :: this
            character(50), intent(out) :: string
        end subroutine lamToStringInterface
    end interface

end module class_LambdaCalculator