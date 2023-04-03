module class_LambdaCalculatorDiag
  
    use class_LambdaCalculator
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorDiag
        private
        double precision, allocatable :: A0(:)
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
            procedure :: clean => clean
    end type LambdaCalculatorDiag

    interface LambdaCalculatorDiag
        module procedure constructor
    end interface
    private :: constructor


    contains

    
        function constructor(NDimX, A2) result(new)

            implicit none
            type(LambdaCalculatorDiag) :: new
            integer, intent(in) :: NDimX
            double precision :: A2(NDimx*NDimX)

            new%NDimX = NDimX
            new%A2 = A2
            
            allocate(new%A0(new%NDimX*new%NDimX))

        end function


        subroutine calculateLambda(this, Lambda, OmI)

            implicit none
            class(LambdaCalculatorDiag) :: this
            double precision, intent(out) :: Lambda(this%NDimX*this%NDimX)
            double precision, intent(in) :: OmI
            integer :: i, j

            Lambda=0.D0
            Do i=1,this%NDimX
                j = (i-1) * this%NDimX + i
                Lambda(j) = 1 / ( this%A0(j) + OmI**2 )
            EndDo

        end subroutine calculateLambda


        subroutine calculateInitialA(this)

            implicit none
            class(LambdaCalculatorDiag) :: this
            integer :: i

            this%A0 = 0.D0   
            Do i = 1, this%NDimX
                this%A0((i-1)*this%NDimX+i) = this%A2((i-1)*this%NDimX+i)
            EndDo
            this%A2 = this%A2 - this%A0

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


        subroutine clean(this)

            class(LambdaCalculatorDiag) :: this
            deallocate(this%A0, this%A2)

        end subroutine clean


end module class_LambdaCalculatorDiag