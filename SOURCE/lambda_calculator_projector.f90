module class_LambdaCalculatorProjector
  
    use class_LambdaCalculator
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorProjector
        double precision, allocatable :: PMat(:,:)
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
    end type LambdaCalculatorProjector


    contains


    subroutine calculateLambda(this, Lambda, OmI, NDimX, A0)

        implicit none
        class(LambdaCalculatorProjector) :: this
        double precision, intent(out) :: Lambda(NDimX*NDimX)
        integer, intent(in) :: NDimX
        double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        integer :: i, inf1, inf2
        double precision :: WORK0(NDimX*NDimX)
        double precision :: ipiv(NDimX)

        WORK0=0.D0
        Do i=1,NDimX
            WORK0((i-1)*NDimX+i)=OmI**2
        EndDo
        Lambda=A0+WORK0
        Call dgetrf(NDimX,NDimX, Lambda, NDimX, ipiv, inf1)
        Call dgetri(NDimX, Lambda, NDimX, ipiv, work0, NDimX, inf2)

    end subroutine calculateLambda


    subroutine calculateInitialA(this, A0, A2, NDimX)

        implicit none
        class(LambdaCalculatorProjector) :: this
        double precision, intent(out) :: A0(NDimX*NDimX)
        double precision, intent(inout) :: A2(NDimX*NDimX)
        integer, intent(in) :: NDimX
        double precision :: WORK0(NDimX*NDimX)
        integer :: i

        A0=0.D0   
        Do i=1,NDimX
           A0((i-1)*NDimX+i)=A2((i-1)*NDimX+i)
        EndDo
        A2=A2-A0

        ! PROJECT A2 WITH PMat : WORK0=PMat.A2
        Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,this%PMat,NDimX,A2,NDimX,0.0,WORK0,NDimX)
        A0=A0+WORK0
        print*, 'A2 norm before projection',norm2(A2(1:NDimX*NDimX))
        A2=A2-WORK0
        print*, 'A2 norm after projection',norm2(A2(1:NDimX*NDimX))

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


end module class_LambdaCalculatorProjector