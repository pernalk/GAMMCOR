module class_LambdaCalculatorProjector
  
    use class_LambdaCalculator
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorProjector
        double precision, allocatable :: PMat(:,:)
        double precision, allocatable :: A0(:)
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
            procedure :: clean => clean
    end type LambdaCalculatorProjector

    interface LambdaCalculatorProjector
        module procedure constructor
    end interface
    private :: constructor


    contains


        function constructor(NDimX, A2, PMat) result(new)

            implicit none
            type(LambdaCalculatorProjector) :: new
            integer, intent(in) :: NDimX
            double precision, intent(in) :: A2(NDimx*NDimX), PMat(NDimX,NDimX)

            new%NDimX = NDimX
            allocate(new%A2(new%NDimX*new%NDimX))
            allocate(new%Pmat(new%NDimX,new%NDimX))

            new%A2 = A2
            new%PMat = PMat
            
            allocate(new%A0(new%NDimX*new%NDimX))

        end function


        subroutine calculateLambda(this, Lambda, OmI)

            implicit none
            class(LambdaCalculatorProjector) :: this
            double precision, intent(out) :: Lambda(this%NDimX*this%NDimX)
            double precision, intent(in) :: OmI
            integer :: i, inf1, inf2
            double precision :: WORK0(this%NDimX*this%NDimX)
            double precision :: ipiv(this%NDimX)

            WORK0=0.D0
            Do i=1,this%NDimX
                WORK0((i-1)*this%NDimX+i)=OmI**2
            EndDo
            Lambda=this%A0+WORK0
            Call dgetrf(this%NDimX,this%NDimX, Lambda, this%NDimX, ipiv, inf1)
            Call dgetri(this%NDimX, Lambda, this%NDimX, ipiv, work0, this%NDimX, inf2)

        end subroutine calculateLambda


        subroutine calculateInitialA(this)

            implicit none
            class(LambdaCalculatorProjector) :: this
            double precision :: WORK0(this%NDimX*this%NDimX)
            integer :: i

            this%A0=0.D0   
            Do i=1,this%NDimX
                this%A0((i-1)*this%NDimX+i) = this%A2((i-1)*this%NDimX+i)
            EndDo
            this%A2 = this%A2 - this%A0

            ! PROJECT A2 WITH PMat : WORK0=PMat.A2
            Call dgemm('N','N',this%NDimX,this%NDimX,this%NDimX,1d0,this%PMat,this%NDimX,this%A2,this%NDimX,0d0,WORK0,this%NDimX)
            this%A0 = this%A0 + WORK0
            print*, 'A2 norm before projection',norm2(this%A2(1:this%NDimX*this%NDimX))
            this%A2 = this%A2 - WORK0
            print*, 'A2 norm after projection',norm2(this%A2(1:this%NDimX*this%NDimX))

        end subroutine calculateInitialA


        function getName(this) result(res)

            implicit none
            class(LambdaCalculatorProjector) :: this
            character(:), allocatable :: res
            res = "Projector"

        end function getName


        subroutine toString(this, string)

            implicit none
            class(LambdaCalculatorProjector) :: this
            character(50), intent(out) :: string
            string = "projector"

        end subroutine toString


        subroutine clean(this)

            class(LambdaCalculatorProjector) :: this
            deallocate(this%A0, this%A2, this%PMat)

        end subroutine clean


end module class_LambdaCalculatorProjector
