module class_LambdaCalcDiag
  
    use class_LambdaCalc
    implicit none
  
    type, public, extends(LambdaCalc):: LambdaCalcDiag
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
    end type LambdaCalcDiag


    contains


    subroutine calculateLambda(this, InitialC, OmI, NDimX, A0)

        implicit none
        class(LambdaCalcDiag) :: this
        double precision, intent(out) :: InitialC(NDimX*NDimX)
        integer, intent(in) :: NDimX
        double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        integer :: i, j

        InitialC=0.D0
        Do i=1,NDimX
            j = (i-1) * NDimX + i
            InitialC(j) = 1 / ( A0(j) + OmI**2 )
        EndDo

    end subroutine calculateLambda


    subroutine calculateInitialA(this, A0, A2, NDimX)

        implicit none
        class(LambdaCalcDiag) :: this
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


end module class_LambdaCalcDiag