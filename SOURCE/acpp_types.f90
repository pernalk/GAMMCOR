module acpp_types

      use iso_fortran_env

      use types

      implicit none


      type TAC0Block
            double precision, dimension(:, :), allocatable :: ASing, ATrip
            double precision, dimension(:), allocatable :: EigS, EigT
            integer, dimension(:), allocatable :: vPlusS, vPlusT
            integer :: order
            integer, dimension(:, :), allocatable :: IndNS, IndNT
            integer, dimension(:, :), allocatable :: IndN1S, IndN1T
            integer, dimension(:, :), allocatable:: IndN2S, IndN2T
            integer, dimension(:, :), allocatable :: MiniMap
            integer :: NDimS, NDimT
            integer :: NDim1S, NDim2S, NDim1T, NDim2T
            character(:), allocatable :: name
            integer :: type
            type(TAC0vaBlock), allocatable :: MiniBlocks(:)
      end type TAC0Block

      type TAC0vaBlock
            double precision, dimension(:,:), allocatable :: MiniAVS, MiniAVT
            double precision, dimension(:), allocatable :: EigS, EigT
            integer, dimension(:), allocatable :: ActListS, ActListT
            integer, dimension(:,:), allocatable :: IndNS, IndNT
            integer :: NDimS, NDimT
      end type TAC0vaBlock


      type TTHCData
            integer :: NTHC, NChol
            double precision, dimension(:,:), allocatable :: Xga, Zgk
            integer :: ExternalOrdering
            double precision, dimension(:), allocatable :: fij, fvw
      end type TTHCData

      type tMxA
            double precision, dimension(:,:), allocatable :: oo
            double precision, dimension(:,:), allocatable :: vv
            double precision, dimension(:,:), allocatable :: aa
            double precision, dimension(:,:), allocatable :: va
            double precision, dimension(:,:), allocatable :: ao
      end type tMxA

      type eigsBlockParams
            double precision, dimension(:), allocatable :: Eig 
            double precision, dimension(:,:), allocatable  :: Eigvec

            integer, dimension(:), allocatable :: v_plus
            !for interaction, vt = 0 is monomer A, vt=1 monomer B, vt=2 mixed                                                                                                       
            integer, dimension(:), allocatable :: vt
            integer, dimension(:, :), allocatable :: IndN 
            integer :: dim
      end type eigsBlockParams

        type :: TLoop
            logical :: run_main = .false.
            logical, dimension(6) :: run_inner = .false.
      end type TLoop



      type TACppData
            integer, dimension(:), allocatable    :: IndX_s, IndX_t
            integer, dimension(:,:), allocatable  :: IndN_s, IndN_t, IndN
            integer, dimension(:), allocatable    :: IndAux
            double precision, dimension(:, :), allocatable :: HNO0, HNO0_THC
            double precision, dimension(:, :, :, :), allocatable :: rdm2_pp
            double precision, dimension(:, :, :, :), allocatable :: rdm2_pm
            double precision, dimension(:), allocatable :: TwoNO
            double precision, dimension(:), allocatable :: XOne
            double precision, allocatable :: Occ(:)
            double precision :: ENuc
            double precision :: ECas
            integer :: general_version, version
            integer :: switch = 0
            integer :: ACType
            integer :: PYSCF = 0 
            integer :: ORCA = 0
            integer :: DALTON = 0
            integer :: NDim, NDim_s, NDim_t
            integer :: NBasis, NA, NI, NIA, NV, NEL
            integer, dimension(:), allocatable    :: IndMod, map
            integer, dimension(:), allocatable    :: IndAuxFirst
            integer, dimension(2) :: nst
            integer :: true_NA, true_NI
            integer :: NInte1, NInte2
            logical :: OnlyEnergy  = .false.
            
      end type TACppData


      type TRdmData
            double precision, dimension(:), allocatable :: R00, R11
            double precision,allocatable :: rdm2_pp(:,:,:,:), rdm2_pm(:,:,:,:)
            ! rdm2_pp1(i,j,k,l)  = rdm2_pp(i, k, j, l)
            double precision,allocatable :: rdm2_pp1(:,:,:,:), rdm2_pm1(:,:,:,:)
            ! rdm2_pp2(i,j,k,l)  = rdm2_pp(i, l, k, j)                                                                                                                          
            double precision,allocatable :: rdm2_pp2(:,:,:,:), rdm2_pm2(:,:,:,:)
            ! like rdm2_pp1 but only active part
            double precision,allocatable :: rdm2_pp1_act(:,:,:,:), rdm2_pm1_act(:,:,:,:)
            ! like rdm2_pp2 but only active part 
            double precision,allocatable :: rdm2_pp2_act(:,:,:,:), rdm2_pm2_act(:,:,:,:)

            double precision,allocatable :: rdm2_pp_act(:,:,:,:), rdm2_pm_act(:,:,:,:)
      end type TRdmData


      type TDA_pphh
            integer, dimension(:,:), allocatable :: IndVO
            integer, dimension(:,:), allocatable :: Ind_O_eq, Ind_O_gt_lt
            integer, dimension(:,:), allocatable :: Ind_V_eq, Ind_V_gt
            integer, dimension(:,:), allocatable :: Ind_O_gt
            integer, dimension(:,:), allocatable :: Indpphh
            integer, dimension(:,:), allocatable :: Indpphhp, Indpphhm
            integer, dimension(:), allocatable :: hh_set
            integer, dimension(:,:), allocatable :: offset
            integer, dimension(:), allocatable :: orb_sym
            integer :: irrep
            integer :: group_order, group_name
            integer :: hh_set_dim
            logical :: allsym
            logical :: nosym = .false. 
            integer :: Npair_virt, Npair_occ, Npair_vvoo, Npair_vo
            integer :: Npair_vvoop, Npair_vvoom
            integer :: vo1, vo2, vo3, vo4
      end type TDA_pphh

      type log_cond
            integer :: Lil, Ljk, Lik, Ljl, Liljk, Likjl
            integer :: Lac, Lbd, Lad, Lbc
            integer :: Lacbd, Ladbc, Lacik, Ladil
            integer :: Lacil, Ladik, Lacjk, Lacjl, Ladjk, Ladjl
            integer :: Lbcbd, Lbdbc, Lbcik, Lbdil
            integer :: Lbcil, Lbdik, Lbcjk, Lbcjl, Lbdjk, Lbdjl

      end type log_cond
contains
      subroutine print_flags2(flags)
            type(FlagsData), intent(in) :: flags

            ! Print integer variables
            write(*,'(A)') '=== FlagsData Configuration ==='
            write(*,'(A,I8)') 'InterfaceType  = ', flags%InterfaceType
            write(*,'(A,I8)') 'IDALTON        = ', flags%IDALTON
            write(*,'(A,I8)') 'IPYSCF         = ', flags%IPYSCF
            write(*,'(A,I8)') 'IMOLPRO        = ', flags%IMOLPRO
            write(*,'(A,I8)') 'IORCA          = ', flags%IORCA
            write(*,'(A,I8)') 'IRes           = ', flags%IRes
            write(*,'(A,I8)') 'IAO            = ', flags%IAO
            write(*,'(A,I8)') 'INO            = ', flags%INO
            write(*,'(A,I8)') 'NoSym          = ', flags%NoSym
            write(*,'(A,I8)') 'NoSt           = ', flags%NoSt
            write(*,'(A,I8)') 'IGVB           = ', flags%IGVB
            write(*,'(A,I8)') 'ITwoEl         = ', flags%ITwoEl
            write(*,'(A,I8)') 'IRedVirt       = ', flags%IRedVirt
            write(*,'(A,I8)') 'IRdm2Typ       = ', flags%IRdm2Typ
            write(*,'(A,I8)') 'IGridType      = ', flags%IGridType
            write(*,'(A,I8)') 'IUnits         = ', flags%IUnits
            write(*,'(A,I8)') 'IFun           = ', flags%IFun
            write(*,'(A,I8)') 'IModG          = ', flags%IModG
            write(*,'(A,I8)') 'NGOcc          = ', flags%NGOcc

            ! Print floating point variables
            write(*,'(A)') '=== Float Values ==='
            write(*,'(A,F12.6)') 'DCholeskyThr   = ', flags%DCholeskyThr
            write(*,'(A,F12.6)') 'DTHCThr        = ', flags%DTHCThr
            write(*,'(A,F12.6)') 'Alpha          = ', flags%Alpha

            ! Print character variables if allocated
            write(*,'(A)') '=== String Values ==='
            if (allocated(flags%JobTitle)) then
                  write(*,'(A,A)') 'JobTitle        = ', trim(flags%JobTitle)
            else
                  write(*,'(A)') 'JobTitle        = <not allocated>'
            endif

            if (allocated(flags%BasisSet)) then
                  write(*,'(A,A)') 'BasisSet        = ', trim(flags%BasisSet)
            else
                  write(*,'(A)') 'BasisSet        = <not allocated>'
            endif

            if (allocated(flags%BasisSetPath)) then
                  write(*,'(A,A)') 'BasisSetPath    = ', trim(flags%BasisSetPath)
            else
                  write(*,'(A)') 'BasisSetPath    = <not allocated>'
            endif
            
      end subroutine print_flags2


          subroutine CalcMem(A, name)
          double precision, dimension(:,:), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2))
          this_mem = this_mem*togb
          write(*, '(A25, A1, F30.20)') name, ',', this_mem
    end subroutine CalcMem

    subroutine CalcMemI(A, name)
          integer, dimension(:,:), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2))
          this_mem = this_mem*togb
          write(*, '(A25, A1, F30.20)') name, ',', this_mem
    end subroutine CalcMemI

    subroutine CalcMem3(A, name)
          double precision, dimension(:,:, :), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2)) * dble(size(A, dim=3))
          this_mem = this_mem*togb
          write(*, '(A25, A1, F30.20)') name, ',', this_mem
    end subroutine CalcMem3

    subroutine CalcMem4(A, name)
          double precision, dimension(:,:, :,:), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2)) * dble(size(A, dim=3)) * dble(size(A, dim=4))
          this_mem = this_mem*togb
          write(*, '(A25, A1, F30.20)') name, ',', this_mem
    end subroutine CalcMem4



end module acpp_types
