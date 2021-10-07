module class_LambdaCalculatorBlock
  
    use class_LambdaCalculator
    use abfofo
    implicit none
  
    type, public, extends(LambdaCalculator):: LambdaCalculatorBlock
        double precision, allocatable :: URe(:,:), Occ(:), XONe(:)
        integer, allocatable :: IndN(:,:), IndX(:), IGem(:)
        integer :: NBasis, NAct, INActive, NInte1
        character(:), allocatable :: twojfile, twokfile
        double precision, allocatable :: A0(:)
        integer :: nblk
        type(EblockData) :: A0blockIV
        type(EblockData), allocatable :: A0block(:)
        contains
            procedure :: calculateInitialA => calculateInitialA
            procedure :: calculateLambda => calculateLambda
            procedure :: getName => getName
            procedure :: toString => toString
            procedure :: clean => clean
    end type LambdaCalculatorBlock

    interface LambdaCalculatorBlock
        module procedure constructor
    end interface
    private :: constructor


    contains


        function constructor(NDimX, A2, URe, Occ, XOne, IndN, IndX, IGem, &
            NBasis, NAct, INActive, NInte1) result(new)

            implicit none
            type(LambdaCalculatorBlock) :: new
            integer, intent(in) :: NDimX
            double precision, intent(in) :: A2(NDimx*NDimX)
            double precision, intent(in) :: URe(NBasis,NBasis), Occ(NBasis), XONe(NInte1)
            integer, intent(in) :: IndN(2,NDimX), IndX(NDimX), IGem(NBasis)
            integer, intent(in) :: NBasis, NAct, INActive, NInte1

            new%NDimX = NDimX
            new%A2 = A2
            new%URe = URe
            new%Occ = Occ
            new%XOne = XOne
            new%IndN = IndN
            new%IndX = IndX
            new%IGem = IGem
            new%NBasis = NBasis
            new%NAct = NAct
            new%INActive = INActive
            new%NInte1 = NInte1

            new%twojfile = 'FFOO'
            new%twokfile = 'FOFO'

            new%nblk = 1 + NBasis - NAct
            
            allocate(new%A0(new%NDimX*new%NDimX))
            allocate(new%A0block(new%nblk))

        end function


        subroutine calculateLambda_old(this, Lambda, OmI, NDimX, A0)

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

        end subroutine calculateLambda_old

        subroutine calculateLambda(this, Lambda, OmI)

            implicit none
            class(LambdaCalculatorBlock) :: this
            double precision, intent(out) :: Lambda(this%NDimX*this%NDimX)
            double precision, intent(in) :: OmI
            integer :: i, inf1, inf2
            double precision :: WORK0(this%NDimX*this%NDimX)
            double precision :: ipiv(this%NDimX)

            Call INV_AC0BLK(OmI**2,Lambda,this%A0Block,this%A0BlockIV,this%nblk,this%NDimX)

        end subroutine calculateLambda

        subroutine calculateInitialA(this)

            use abfofo
            implicit none
            class(LambdaCalculatorBlock) :: this

            double precision :: ABPLUS0(this%NDimX*this%NDimX),WORK0(this%NDimX*this%NDimX)
            double precision :: ECASSCF, ACAlpha0
        
            ACAlpha0=0.D0
            call AB_CAS_FOFO(ABPLUS0,WORK0,ECASSCF,this%URe,this%Occ,this%XOne, &
                            this%IndN,this%IndX,this%IGem,this%NAct,this%INActive,this%NDimX,this%NBasis,this%NDimX,&
                            this%NInte1,this%twojfile,this%twokfile,ACAlpha0,.false.)
            Call dgemm('N','N',this%NDimX,this%NDimX,this%NDimX,1d0,ABPLUS0,this%NDimX,WORK0,this%NDimX,0.0,this%A0,this%NDimX)
            this%A2 = this%A2 - this%A0

            Call AC0BLOCK(this%Occ,this%URe,this%XOne, &
                this%IndN,this%IndX,this%IGem,this%NAct,this%INActive,this%NDimX,this%NBasis,this%NDimX,&
                this%NInte1,'FFOO','FOFO', &
                this%A0BlockIV,this%A0Block,this%nblk,'A0BLK',0)

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


        subroutine clean(this)

            class(LambdaCalculatorBlock) :: this

            deallocate(this%A0, this%A2, this%URe, this%Occ, this%XOne, &
                this%IndN, this%IndX, this%IGem, this%twojfile, this%twokfile)

            Call RELEASE_AC0BLOCK(this%A0Block,this%A0blockIV,this%nblk)

        end subroutine clean


end module class_LambdaCalculatorBlock
