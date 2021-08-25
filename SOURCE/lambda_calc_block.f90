module class_LambdaCalcBlock
  
    use class_LambdaCalc
    implicit none
  
    type, public, extends(LambdaCalc):: LambdaCalcBlock
        double precision, allocatable :: URe(:,:), Occ(:), XONe(:)
        integer, allocatable :: IndN(:,:), IndX(:), IGem(:)
        integer :: NBasis, NAct, INActive, NInte1
        character(:), allocatable :: twojfile, twokfile
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
    end type LambdaCalcBlock


    contains


    subroutine calculateLambda(this, InitialC, OmI, NDimX, A0)

        implicit none
        class(LambdaCalcBlock) :: this
        double precision, intent(out) :: InitialC(NDimX*NDimX)
        integer, intent(in) :: NDimX
        double precision, intent(in) :: OmI, A0(NDimX*NDimX)
        integer :: i, inf1, inf2
        double precision :: WORK0(NDimX*NDimX)
        double precision :: ipiv(NDimX)

        WORK0=0.D0
        Do i=1,NDimX
            WORK0((i-1)*NDimX+i)=OmI**2
        EndDo
        InitialC=A0+WORK0
        Call dgetrf(NDimX,NDimX, InitialC, NDimX, ipiv, inf1)
        Call dgetri(NDimX, InitialC, NDimX, ipiv, work0, NDimX, inf2)

    end subroutine calculateLambda


    subroutine calculateInitialA(this, A0, A2, NDimX)

        use abfofo
        implicit none
        class(LambdaCalcBlock) :: this
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


end module class_LambdaCalcBlock