module class_LambdaCalculator

    implicit none
  
    type, public, abstract :: LambdaCalculator
        integer :: NDimX
        double precision, allocatable :: A2(:)
        contains
            procedure(calculateInitialAInterface), deferred :: calculateInitialA
            procedure(calculateLambdaInterface), deferred :: calculateLambda
            procedure(lamGetName), deferred ::getName
            procedure(lamToStringInterface), deferred :: toString
            procedure(cleanInterface), deferred :: clean
    end type LambdaCalculator

    interface
        subroutine calculateLambdaInterface(this, Lambda, OmI)
            import LambdaCalculator
            implicit none
            class(LambdaCalculator) :: this
            double precision, intent(out) :: Lambda(this%NDimX*this%NDimX)
            double precision, intent(in) :: OmI
        end subroutine calculateLambdaInterface
    end interface

    interface
        subroutine calculateInitialAInterface(this)
            import LambdaCalculator
            implicit none
            class(LambdaCalculator) :: this
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

    interface
        subroutine cleanInterface(this)
            import LambdaCalculator
            class(LambdaCalculator) :: this
        end subroutine cleanInterface
    end interface

end module class_LambdaCalculator