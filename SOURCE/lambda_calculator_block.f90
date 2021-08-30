module class_LambdaCalculatorBlock
  
    use class_LambdaCalculator
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorBlock
        double precision, allocatable :: URe(:,:), Occ(:), XONe(:)
        integer, allocatable :: IndN(:,:), IndX(:), IGem(:)
        integer :: NBasis, NAct, INActive, NInte1
        character(:), allocatable :: twojfile, twokfile
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
    end type LambdaCalculatorBlock


    contains


    subroutine calculateLambda(this, Lambda, OmI, NDimX, A0)

        implicit none
        class(LambdaCalculatorBlock) :: this
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

        use abfofo
        implicit none
        class(LambdaCalculatorBlock) :: this
        double precision, intent(out) :: A0(NDimX*NDimX)
        double precision, intent(inout) :: A2(NDimX*NDimX)
        integer, intent(in) :: NDimX

        double precision :: ABPLUS0(NDimX*NDimX),WORK0(NDimX*NDimX)
        double precision :: ECASSCF, ACAlpha0
    
        ACAlpha0=0.D0
        call AB_CAS_FOFO(ABPLUS0,WORK0,ECASSCF,this%URe,this%Occ,this%XOne, &
                        this%IndN,this%IndX,this%IGem,this%NAct,this%INActive,NDimX,this%NBasis,NDimX,&
                        this%NInte1,this%twojfile,this%twokfile,ACAlpha0,.false.)
        Call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS0,NDimX,WORK0,NDimX,0.0,A0,NDimX)
        A2=A2-A0

    end subroutine calculateInitialA


    function getName(this) result(res)

        implicit none
        class(LambdaCalculatorBlock) :: this
        character(:), allocatable :: res
        res = "Block"

    end function getName


    subroutine toString(this, string)

        implicit none
        class(LambdaCalculatorBlock) :: this
        character(50), intent(out) :: string
        string = "block"

    end subroutine toString


end module class_LambdaCalculatorBlock