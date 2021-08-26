module class_LambdaCalc

    implicit none
  
    type, public, abstract :: LambdaCalc
        contains
            procedure(calculateInitialAInterface), deferred :: calculateInitialA
            procedure(calculateLambdaInterface), deferred :: calculateLambda
    end type LambdaCalc

    interface
        subroutine calculateLambdaInterface(this, InitialC, OmI, NDimX, A0)
            import LambdaCalc
            implicit none
            class(LambdaCalc) :: this
            double precision, intent(out) :: InitialC(NDimX*NDimX)
            integer, intent(in) :: NDimX
            double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        end subroutine calculateLambdaInterface
    end interface

    interface
        subroutine calculateInitialAInterface(this, A0, A2, NDimX)
            import LambdaCalc
            implicit none
            class(LambdaCalc) :: this
            double precision, intent(out) :: A0(NDimX*NDimX)
            double precision, intent(inout) :: A2(NDimX*NDimX)
            integer, intent(in) :: NDimX
        end subroutine calculateInitialAInterface
    end interface

end module class_LambdaCalc