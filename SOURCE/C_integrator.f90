module class_CIntegrator

    use class_IterAlgorithm
    use class_LambdaCalculator
    use class_IterStats
    implicit none

    type, public :: CIntegrator
        integer :: NDimX = 0, NCholesky = 0, NGrid = 0
        double precision :: ACAlpha = 0, XFreq(100) = 0, WFreq(100) = 0
        double precision, allocatable :: APlusTilde(:), APlusTildeAct(:)
        class(IterAlgorithm), pointer :: iterAlgo
        class(LambdaCalculator), pointer :: lambdaCalc
        type(IterStats) :: iStats = IterStats()
        contains
            procedure :: setup => setup
            procedure :: integrate => integrate
            procedure :: integrateReverse => integrateReverse
            procedure :: showInfo => showInfo
            procedure :: clean => clean
    end type CIntegrator


    contains


        subroutine setup(this, NGrid, NDimX, NCholesky, APlusTilde, APlusTildeAct, ACAlpha)

            implicit none
            class(CIntegrator) :: this
            integer :: NDimX, NCholesky, NGrid
            double precision :: ACAlpha, XFreq(100), WFreq(100)
            double precision :: APlusTilde(NDimX*NCholesky), APlusTildeAct(NDimX*NCholesky)
            this%NDimX = NDimX
            this%NCholesky = NCholesky
            this%NGrid = NGrid
            this%ACAlpha = ACAlpha
            this%APlusTilde = APlusTilde
            this%APlusTildeAct = APlusTildeAct

            call FreqGrid(this%XFreq, this%WFreq, NGrid)

        end subroutine setup


        subroutine integrate(this, COMTilde, COMTildeAct)

            implicit none
            class(CIntegrator) :: this
            double precision, allocatable :: COMTilde(:), COMTildeAct(:) 

            integer :: IGL
            double precision :: PI, OmI, Lambda(this%NDimX*this%NDimX)
            double precision :: CTilde(this%NDimX*this%NCholesky), CTildeAct(this%NDimX*this%NCholesky)

            PI = 4.0*ATAN(1.0)

            associate(NDimX => this%NDimX, NCholesky => this%NCholesky, &
                        NGrid => this%NGrid, ACAlpha => this%ACAlpha, &
                        lambdaCalc => this%lambdaCalc, iterAlgo => this%iterAlgo, iStats => this%iStats, &
                        APlusTilde => this%APlusTilde, APlusTildeAct => this%APlusTildeAct, &
                        XFreq => this%XFreq, WFreq => this%WFreq)

                call showInfo(this, "normal order")

                call iStats%reset()
                iStats%maxIterationsLimit = iterAlgo%maxIterations

                call lambdaCalc%calculateInitialA()

                Do IGL=1,NGrid

                    OmI = XFreq(IGL)

                    call lambdaCalc%calculatelambda(Lambda, OmI)

                    Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,Lambda,NDimX,APlusTilde,NDimX,0.0d0,CTilde,NDimX)
                    call iStats%setFreq(OmI)
                    call iterAlgo%iterate(CTilde, NDimX, NCholesky, Lambda, APlusTilde, lambdaCalc%A2, iStats)

                    Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,Lambda,NDimX,APlusTildeAct,NDimX,0.0d0,CTildeAct,NDimX)
                    call iStats%setFreq(OmI)
                    call iterAlgo%iterate(CTildeAct, NDimX, NCholesky, Lambda, APlusTildeAct, lambdaCalc%A2, iStats)

                    COMTilde = COMTilde + 2.D0 / PI * CTilde * WFreq(IGL)
                    COMTildeAct = COMTildeAct + 2.D0 / PI * CTildeAct * WFreq(IGL)

                EndDo

                call iStats%print()
                ! call saveStatistic(this)

            end associate

        end subroutine integrate


        subroutine integrateReverse(this, COMTilde, COMTildeAct)

            use class_IterStats
            implicit none
            class(CIntegrator) :: this
            double precision, allocatable :: COMTilde(:), COMTildeAct(:) 

            integer :: IGL
            double precision :: PI, OmI, Lambda(this%NDimX*this%NDimX)
            double precision :: CTilde(this%NDimX*this%NCholesky), CTildeAct(this%NDimX*this%NCholesky)

            PI = 4.0*ATAN(1.0)

            associate(NDimX => this%NDimX, NCholesky => this%NCholesky, &
                NGrid => this%NGrid, ACAlpha => this%ACAlpha, &
                lambdaCalc => this%lambdaCalc, iterAlgo => this%iterAlgo, iStats => this%iStats, &
                APlusTilde => this%APlusTilde, APlusTildeAct => this%APlusTildeAct, &
                XFreq => this%XFreq, WFreq => this%WFreq) 

                call showInfo(this, "reverse order")

                call iStats%reset()
                iStats%maxIterationsLimit = iterAlgo%maxIterations

                call lambdaCalc%calculateInitialA()
                call lambdaCalc%calculatelambda(Lambda, XFreq(NGrid))
                Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,Lambda,NDimX,APlusTilde,NDimX,0.0d0,CTilde,NDimX)
                Call dgemm('N','N',NDimX,NCholesky,NDimX,1.d0,Lambda,NDimX,APlusTildeAct,NDimX,0.0d0,CTildeAct,NDimX)

                Do IGL=NGrid,1,-1

                    OmI = XFreq(IGL)

                    if(IGL < NGrid) call lambdaCalc%calculatelambda(Lambda, OmI)

                    call iStats%setFreq(OmI)
                    call iterAlgo%iterate(CTilde, NDimX, NCholesky, Lambda, APlusTilde, lambdaCalc%A2, iStats)
                    ! call iStats%add(OmI,N)

                    call iStats%setFreq(OmI)
                    call iterAlgo%iterate(CTildeAct, NDimX, NCholesky, Lambda, APlusTildeAct, lambdaCalc%A2, iStats)
                    ! call iStats%add(OmI,N)

                    COMTilde = COMTilde + 2.D0 / PI * CTilde * WFreq(IGL)
                    COMTildeAct = COMTildeAct + 2.D0 / PI * CTildeAct * WFreq(IGL)

                EndDo

                call iStats%print()
                ! call saveStatistic(this)

            end associate

        end subroutine integrateReverse


        subroutine showInfo(this, integrationModeName)

            class(CIntegrator) :: this
            character(*) :: integrationModeName
            print ("(a,a)"),        " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            print ("(a,a)"),        " [CIntegrator] C integration info"
            print ("(a,a)"),        " Integration mode:     ", trim(integrationModeName)
            print ("(2a,i0,a)"),    " Integration params:   ", "{NGrid: ", this%NGrid, "}"
            print ("(a,a)"),        " Iteration algorithm:  ", this%iterAlgo%getName()
            print ("(a,a)"),        " Iteration params:     ", this%iterAlgo%getParamsString()
            print ("(a,a)"),        " Lambda calculator:    ", this%lambdaCalc%getName()
            print ("(a,a)"),        " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            
        end subroutine showInfo


        subroutine saveStatistic(this)

            class(CIntegrator) :: this
            character(50) :: iterAlgoString
            character(50) :: lambdaCalcString
            character(100) :: fileName

            call this%iterAlgo%toString(iterAlgoString)
            call this%lambdaCalc%toString(lambdaCalcString)
            write(fileName, "(4a)") "CInteg_", trim(iterAlgoString), "_", trim(lambdaCalcString) 
            call this%iStats%dump(fileName, this%ACAlpha)

        end subroutine saveStatistic


        subroutine clean(this)

            class(CIntegrator) :: this
            deallocate(this%APlusTilde, this%APlusTildeAct)

        endsubroutine clean

end module class_CIntegrator