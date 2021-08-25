module class_CIntegrator

    use class_IterAlgorithm
    use class_LambdaCalc
    implicit none

    type, public :: CIntegrator
        integer :: NDimX = 0, NCholesky = 0, NGrid = 0
        double precision :: ACAlpha = 0, XFreq(100) = 0, WFreq(100) = 0
        double precision, allocatable :: A0(:), A2(:), APlusTilde(:), APlusTildeAct(:)
        class(IterAlgorithm), pointer :: iterAlgo
        class(LambdaCalc), pointer :: LambdaCalc
    contains
        procedure :: setup => setup
        procedure :: integrate => integrate
        procedure :: integrateReverse => integrateReverse
    end type CIntegrator


    contains


    subroutine setup(this, NDimX, NCholesky, A0, A2, APlusTilde, APlusTildeAct, XFreq, WFreq, NGrid, ACAlpha)

        implicit none
        class(CIntegrator) :: this
        integer :: NDimX, NCholesky, NGrid
        double precision :: ACAlpha, XFreq(100), WFreq(100)
        double precision :: A0(NDimX*NDimX), A2(NDimX*NDimX), APlusTilde(NDimX*NCholesky), APlusTildeAct(NDimX*NCholesky)
        this%NDimX = NDimX
        this%NCholesky = NCholesky
        this%NGrid = NGrid
        this%ACAlpha = ACAlpha
        this%XFreq = XFreq
        this%WFreq = WFreq
        this%A0 = A0
        this%A2 = A2
        this%APlusTilde = APlusTilde
        this%APlusTildeAct = APlusTildeAct

    end subroutine setup


    subroutine integrate(this, COMTilde, COMTildeAct)

        use class_IterStats
        implicit none
        class(CIntegrator) :: this
        double precision, allocatable :: COMTilde(:), COMTildeAct(:) 

        integer :: IGL, N
        double precision :: PI, OmI, InitialC(this%NDimX*this%NDimX)
        double precision :: CTilde(this%NDimX*this%NCholesky), CTildeAct(this%NDimX*this%NCholesky)
        character(50) :: fileName
        type(IterStats) :: iStats

        PI = 4.0*ATAN(1.0)

        iStats = IterStats()

        associate(NDimX => this%NDimX, NCholesky => this%NCholesky, &
                    NGrid => this%NGrid, ACAlpha => this%ACAlpha, &
                    LambdaCalc => this%LambdaCalc, iterAlgo => this%iterAlgo, &
                    APlusTilde => this%APlusTilde, APlusTildeAct => this%APlusTildeAct, &
                    A0 => this%A0, A2 => this%A2, XFreq => this%XFreq, WFreq => this%WFreq) 

            call LambdaCalc%calculateInitialA(A0, A2, NDimX)

            Do IGL=1,NGrid

                OmI = XFreq(IGL)

                call LambdaCalc%calculateLambda(InitialC, OmI, NDimX, A0)

                Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,InitialC,NDimX,APlusTilde,NDimX,0.0d0,CTilde,NDimX)
                call iterAlgo%iterate(CTilde, N, NDimX, NCholesky, InitialC, APlusTilde, A2)
                call iStats%add(OmI,N)

                Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,InitialC,NDimX,APlusTildeAct,NDimX,0.0d0,CTildeAct,NDimX)
                call iterAlgo%iterate(CTildeAct, N, NDimX, NCholesky, InitialC, APlusTildeAct, A2)
                call iStats%add(OmI,N)

                COMTilde = COMTilde + 2.D0 / PI * CTilde * WFreq(IGL)
                COMTildeAct = COMTildeAct + 2.D0 / PI * CTildeAct * WFreq(IGL)

            EndDo

            call iStats%print()
            call iterAlgo%createFileName(fileName)
            call iStats%dump(fileName, ACAlpha)

        end associate

    end subroutine integrate


    subroutine integrateReverse(this, COMTilde, COMTildeAct)

        use class_IterStats
        implicit none
        class(CIntegrator) :: this
        double precision, allocatable :: COMTilde(:), COMTildeAct(:) 

        integer :: IGL, N
        double precision :: PI, OmI, InitialC(this%NDimX*this%NDimX)
        double precision :: CTilde(this%NDimX*this%NCholesky), CTildeAct(this%NDimX*this%NCholesky)
        character(50) :: fileName
        type(IterStats) :: iStats

        PI = 4.0*ATAN(1.0)

        iStats = IterStats()

        associate(NDimX => this%NDimX, NCholesky => this%NCholesky, &
                    NGrid => this%NGrid, ACAlpha => this%ACAlpha, &
                    LambdaCalc => this%LambdaCalc, iterAlgo => this%iterAlgo, &
                    APlusTilde => this%APlusTilde, APlusTildeAct => this%APlusTildeAct, &
                    A0 => this%A0, A2 => this%A2, XFreq => this%XFreq, WFreq => this%WFreq) 

            call LambdaCalc%calculateInitialA(A0, A2, NDimX)

            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,InitialC,NDimX,APlusTilde,NDimX,0.0d0,CTilde,NDimX)
            Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,InitialC,NDimX,APlusTildeAct,NDimX,0.0d0,CTildeAct,NDimX)

            Do IGL=NGrid,1,-1

                OmI = XFreq(IGL)

                call LambdaCalc%calculateLambda(InitialC, OmI, NDimX, A0)

                call iterAlgo%iterate(CTilde, N, NDimX, NCholesky, InitialC, APlusTilde, A2)
                call iStats%add(OmI,N)

                call iterAlgo%iterate(CTildeAct, N, NDimX, NCholesky, InitialC, APlusTildeAct, A2)
                call iStats%add(OmI,N)

                COMTilde = COMTilde + 2.D0 / PI * CTilde * WFreq(IGL)
                COMTildeAct = COMTildeAct + 2.D0 / PI * CTildeAct * WFreq(IGL)

            EndDo

            call iStats%print()
            call iterAlgo%createFileName(fileName)
            call iStats%dump(fileName, ACAlpha)

        end associate

    end subroutine integrateReverse

end module class_CIntegrator