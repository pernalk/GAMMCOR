module phac0_thc

      ! =============================================================================
      ! ph-AC0 with THC integrals: A0 blocks (assembly + diagonalization), A1 blocks
      ! (int loop + delta sector + SA-CASSCF symmetrization), then the energy in
      ! thc_energy_ph. With THCData%exact_ints every slab comes from the exact
      ! integrals (thc_exact_ints). The letter o means INACTIVE here.
      ! Notes: .secret_sauce/.skills/ph-ac0/code-thc-a0.md, code-thc-a1.md,
      ! code-thc-energy.md, code-thc-exact-ints.md (read them before changing this module).
      ! =============================================================================

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals
      use types
      use lin
      use clock

      use math_constants
      use pp_utils
      use print_utils
      use ppac_types
      use thc_energy_ph
      use thc_exact_ints, only: exact_V_axby, exact_V_pair
      use tran
      use omp_lib


      implicit none

      ! ph block ids. Prefixed ph* so they cannot be confused with the pp ids
      ! A0oo..A1vaao2 of ppac0_thc (both sets are private to their module).
      ! The A1 ids run contiguously from phA1aavo to NBLOCKS_PH. thc_energy_ph
      ! keeps a copy of all ids - change both together.
      integer, parameter, private :: phA0vo=1, phA0aa=2, phA0va=3, phA0oa=4
      ! A1 block ids phA1<row pair><column pair>, higher space first in each pair
      ! (the block names of theory-a1-blocks/ without the dash). 11-13 are the
      ! transposed partners that symmetrize_A1_ph folds into aova/aoaa/aava and frees.
      ! vovo, vova, vava are never built. 'ao' here = 'oa' in phA0oa (same pair space).
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 1.
      integer, parameter, private :: phA1aavo=5, phA1aova=6, phA1aovo=7
      integer, parameter, private :: phA1aoao=8, phA1aava=9, phA1aoaa=10
      integer, parameter, private :: phA1vaao=11, phA1aaao=12, phA1vaaa=13
      integer, parameter, private :: NBLOCKS_PH = 13

      integer, parameter, private :: tdebug = 0
      integer, parameter, private :: mdebugg = 50

      character(len=8), dimension(NBLOCKS_PH), parameter :: BlockName_ph = (/ &
            'A0vo    ', 'A0aa    ', 'A0va    ', 'A0oa    ', &
            'A1aavo  ', 'A1aova  ', 'A1aovo  ', &
            'A1aoao  ', 'A1aava  ', 'A1aoaa  ', &
            'A1vaao  ', 'A1aaao  ', 'A1vaaa  ' /)

contains


      subroutine ACPH0_THC(THCData, AuxData, CAONO_IN, Flags, E_corr, DeltaH)
            use THC_Gammcor

            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO_IN
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData
            ! ph-AC0 correlation energy of the reference I; MS-AC0: the row
            ! Delta H_JI, J = AuxData%trdm_st(:,1:NJ) (size 0 for plain AC0)
            double precision, intent(out) :: E_corr
            double precision, dimension(:), allocatable, intent(out) :: DeltaH

            type(TAC0Block_ph), dimension(NBLOCKS_PH) :: ACBlocks
            integer, dimension(:,:), allocatable :: pos
            double precision, dimension(:,:), allocatable :: AuxII, AuxAA, BuxII, BuxAA
            double precision, dimension(:,:), allocatable :: Aux3X
            double precision, dimension(:,:), allocatable :: HNO
            double precision, dimension(:,:,:,:), allocatable :: rdm2, rdm2_12, rdm2_13ph
            integer :: i
            type (tclock) :: timer


            if (AuxData%PYSCF==1)then
                  call print_info('External ordering', 'PySCF')
            end if

            call set_frozen_ph(AuxData, Flags)

            associate(NBasis=>AuxData%NBasis, map=>AuxData%map, Occ=>AuxData%Occ, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV)

              allocate(pos(NBasis, NBasis))
              pos = 0

              call init_AC0Blocks_ph(AuxData, ACBlocks, pos)

              allocate(AuxII(NBasis, NBasis))
              allocate(BuxII(NBasis, NBasis))
              allocate(AuxAA(NBasis, NBasis))
              allocate(BuxAA(NBasis, NBasis))
              allocate(Aux3X(NBasis, NBasis))

              AuxII = zero
              BuxII = zero
              AuxAA = zero
              BuxAA = zero
              Aux3X = zero

              ! The spin-summed active 2-RDM, plus the two index-permuted copies
              ! the genuine A1 kernels contract against (see init_rdm_perm_ph).
              allocate(rdm2(NA, NA, NA, NA))
              allocate(rdm2_12(NA, NA, NA, NA))
              allocate(rdm2_13ph(NA, NA, NA, NA))
              call clock_start(timer)
              rdm2 = 0.5d0 * AuxData%rdm2_full
              call init_rdm_perm_ph(rdm2, rdm2_12, rdm2_13ph, NA)
              call tmsg('TIME FOR RDM2 PERMUTED COPIES', timer, tdebug)

              call clock_start(timer)
              call THC_int_loop_ph(THCData, ACBlocks, AuxII, AuxAA, BuxII, BuxAA, &
                    Aux3X, Flags, AuxData, Occ, map, rdm2, rdm2_12, rdm2_13ph, pos, One)
              call tmsg('TIME FOR THC LOOPS', timer, tdebug)

              ! HNO (Dyall operator, update_HNO_AC0) is for the A0 blocks only. The A1 delta
              ! sector needs the BARE h, AuxData%HNO0: HNO is zero on cross-space elements
              ! and its fij/fvw already hold the active dressing (.secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 4).
              allocate(HNO(NBasis, NBasis))
              call clock_start(timer)
              call update_HNO_AC0(AuxData, HNO, THCData%fij, THCData%fvw, AuxII, BuxII)
              call tmsg('TIME FOR UPDATE HNO', timer, tdebug)

              call clock_start(timer)
              call PHERPA0_THC_VO(AuxData, ACBlocks(phA0vo), THCData)
              call tmsg('Hessian and diagonalization '//trim(BlockName_ph(phA0vo)), timer, tdebug)

              call clock_start(timer)
              call PHERPA0_THC_AA(AuxData, ACBlocks(phA0aa), HNO, Aux3X)
              call tmsg('Hessian and diagonalization '//trim(BlockName_ph(phA0aa)), timer, tdebug)

              call clock_start(timer)
              call PHERPA0_THC_VA_Eigs(AuxData, ACBlocks(phA0va), THCData, HNO, Aux3X)
              call tmsg('Hessian and diagonalization '//trim(BlockName_ph(phA0va)), timer, tdebug)

              call clock_start(timer)
              call PHERPA0_THC_OA_Eigs(AuxData, ACBlocks(phA0oa), THCData, HNO, Aux3X, AuxAA, BuxAA)
              call tmsg('Hessian and diagonalization '//trim(BlockName_ph(phA0oa)), timer, tdebug)

              ! A1 blocks. THC_int_loop_ph has already put the delta-free part
              ! (occupation-number t*/d* kernels plus genuine T*/D* kernels) into
              ! every block; PHERPA_A1_THC adds the delta sector built from the
              ! BARE h (AuxData%HNO0, not HNO - see above) plus AuxII / BuxII /
              ! AuxAA / BuxAA / Aux3X. One generic call per block, exactly as pp
              ! does with PPERPA_THC in its do i = 9, 13 loop.
              do i = phA1aavo, NBLOCKS_PH
                    call clock_start(timer)
                    call PHERPA_A1_THC(AuxData, ACBlocks(i), AuxData%HNO0, &
                          AuxII, BuxII, AuxAA, BuxAA, Aux3X)
                    call tmsg('Hessian '//trim(BlockName_ph(i)), timer, tdebug)
              end do

              ! A_RC <- 1/2 (A_RC + A_CR^T) for aova, aoaa, aava; frees the
              ! transposed partners
              call symmetrize_A1_ph(ACBlocks)

              call clock_start(timer)
              call THC_energy_loop_ph(THCData, ACBlocks, Flags, pos, AuxData, rdm2, E_corr, DeltaH)
              call tmsg('Correlation energy', timer, tdebug)

            end associate

      end subroutine ACPH0_THC


      ! ph analogue of init_AC0Blocks (ppac0_thc). Builds the 4 A0 pair lists
      ! and allocates the 6 stored A1 coupling blocks (no vovo, vova, vava -
      ! see the module header) and the 3 transposed partners.
      subroutine init_AC0Blocks_ph(AuxData, AC0Block, pos)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), dimension(NBLOCKS_PH), intent(inout) :: AC0Block
            integer, dimension(:,:), intent(inout) :: pos
            integer, parameter :: SP_OCC = 0, SP_ACT = 1, SP_VIRT = 2
            integer :: i

            do i = 1, 4
                  AC0Block(i)%order = 0
            end do

            do i = 5, NBLOCKS_PH
                  AC0Block(i)%order = 1
            end do

            pos = zero

            call InitA0Blocks_ph(AC0Block(phA0vo), AuxData, pos, SP_VIRT, SP_OCC)
            call InitA0Blocks_ph(AC0Block(phA0aa), AuxData, pos, SP_ACT, SP_ACT)
            call InitA0Blocks_ph(AC0Block(phA0va), AuxData, pos, SP_VIRT, SP_ACT)
            call InitA0Blocks_ph(AC0Block(phA0oa), AuxData, pos, SP_OCC, SP_ACT)

            do i = 1, NBLOCKS_PH
                  AC0Block(i)%name = BlockName_ph(i)
            end do

            call print_section('phAC0 block dimensions')
            write(*,'(2X,A10,2X,A12)') 'Block', 'NDim'
            write(*,'(2X,A)') repeat('-', 26)
            do i = 1, 4
                  write(*,'(2X,A10,2X,I12)') trim(AC0Block(i)%name), AC0Block(i)%NDim
            end do

            ! The two blocks with a vo side, aavo and aovo, are the large ones
            ! (NA**2 * NI*NV and NA*NI * NI*NV doubles, twice).
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aavo, phA0aa, phA0vo)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aova, phA0oa, phA0va)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aovo, phA0oa, phA0vo)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aoao, phA0oa, phA0oa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aava, phA0aa, phA0va)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aoaa, phA0oa, phA0aa)
            ! transposed partners, folded into aova / aoaa / aava by
            ! symmetrize_A1_ph before the energy
            call InitA1Blocks_ph(AC0Block, AuxData, phA1vaao, phA0va, phA0oa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aaao, phA0aa, phA0oa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1vaaa, phA0va, phA0aa)

      end subroutine init_AC0Blocks_ph


      ! ph analogue of InitA1Blocks (ppac0_thc). Rows from parent A0ind1,
      ! columns from parent A0ind2. Coupling matrices, never diagonalised.
      subroutine InitA1Blocks_ph(AC0Block, AuxData, A1ind, A0ind1, A0ind2)
            type(TAC0Block_ph), dimension(:), intent(inout) :: AC0Block
            type(TACppData), intent(in) :: AuxData
            integer, intent(in) :: A1ind, A0ind1, A0ind2
            integer :: NDim1t, NDim2t

            NDim1t = AC0Block(A0ind1)%NDim
            NDim2t = AC0Block(A0ind2)%NDim

            AC0Block(A1ind)%NDim1 = NDim1t
            AC0Block(A1ind)%NDim2 = NDim2t

            allocate(AC0Block(A1ind)%IndN1(2, NDim1t))
            allocate(AC0Block(A1ind)%IndN2(2, NDim2t))
            call CalcMemI(AC0Block(A1ind)%IndN1, 'AC0Block(A1ind)%IndN1')
            call CalcMemI(AC0Block(A1ind)%IndN2, 'AC0Block(A1ind)%IndN2')

            AC0Block(A1ind)%IndN1 = AC0Block(A0ind1)%IndN
            AC0Block(A1ind)%IndN2 = AC0Block(A0ind2)%IndN

            AC0Block(A1ind)%type = A1ind

            allocate(AC0Block(A1ind)%APlus(NDim1t, NDim2t))
            allocate(AC0Block(A1ind)%AMinus(NDim1t, NDim2t))

            call CalcMem(AC0Block(A1ind)%APlus, 'AC0Block(A1ind)%APlus')
            call CalcMem(AC0Block(A1ind)%AMinus, 'AC0Block(A1ind)%AMinus')

            AC0Block(A1ind)%APlus = zero
            AC0Block(A1ind)%AMinus = zero

      end subroutine InitA1Blocks_ph


      ! Frozen core: number of inactive orbitals excluded from the ph pair space.
      !   Core 1                   -> 0
      !   Core 0, NCoreOrb > 0     -> NCoreOrb
      !   Core 0, NCoreOrb = 0     -> NI (all inactive, as phac_dispatch does)
      ! Frozen orbitals stay doubly occupied in the reference: every sum over
      ! occupied orbitals (AuxII, BuxII, fij, HNO, Fel_ph, ECAS) keeps 1..NI.
      ! Only loops over the inactive member of a ph pair start at NFrozen + 1.
      subroutine set_frozen_ph(AuxData, Flags)
            type(TACppData), intent(inout) :: AuxData
            type(FlagsData), intent(in) :: Flags

            if (Flags%IFlCore == 1) then
                  AuxData%NFrozen = 0
            else if (AuxData%NCoreOrb > 0) then
                  AuxData%NFrozen = AuxData%NCoreOrb
            else
                  AuxData%NFrozen = AuxData%NI
            end if

            if (AuxData%NFrozen < 0 .or. AuxData%NFrozen > AuxData%NI) then
                  write(*,'(2X,A,I0,A,I0,A)') 'set_frozen_ph: NCoreOrb = ', AuxData%NFrozen, &
                        ' outside 0..NI (NI = ', AuxData%NI, ')'
                  stop
            end if

            call print_info('Frozen core orbitals (ph pairs)', AuxData%NFrozen)

      end subroutine set_frozen_ph


      ! Pair list, pos(p,q) map and storage of one A0 block (vo, aa, va, oa), first
      ! slot = higher space; pairs with a frozen inactive (< NI0) get pos = 0.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a0.md Sec. 5.
      subroutine InitA0Blocks_ph(A0Block, AuxData, pos, Ip, Iq)
            type(TAC0Block_ph), intent(inout) :: A0Block
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:,:), intent(inout) :: pos
            integer, intent(in) :: Ip, Iq
            integer, parameter :: SP_OCC = 0, SP_ACT = 1, SP_VIRT = 2

            integer :: ind, inda
            integer :: i, j, k, Nbl, Nbl0, Nbl1, NI0

            associate(NA=>AuxData%NA, NV=>AuxData%NV, NI=>AuxData%NI, &
                  NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, Occ=>AuxData%Occ)

              NI0 = AuxData%NFrozen + 1

              if (Ip == SP_VIRT .and. Iq == SP_OCC)then
                    ! ------------- vo block: diagonal -------------
                    A0Block%type = phA0vo
                    A0Block%NDim = NV * (NI - NI0 + 1)
                    allocate(A0Block%IndN(2, A0Block%NDim))

                    ind = 0
                    do i = NIA+1, NBasis
                          do j = NI0, NI
                                ind = ind + 1
                                A0Block%IndN(1, ind) = i
                                A0Block%IndN(2, ind) = j
                                pos(i, j) = ind
                          end do
                    end do

                    ! diagonal block: one column holds the diagonal
                    allocate(A0Block%APlus(A0Block%NDim, 1))
                    allocate(A0Block%AMinus(A0Block%NDim, 1))
                    A0Block%APlus = zero
                    A0Block%AMinus = zero

                    allocate(A0Block%Eig(A0Block%NDim))
                    allocate(A0Block%vPlus(A0Block%NDim))
                    A0Block%Eig = zero
                    A0Block%vPlus = 0

              else if (Ip == SP_ACT .and. Iq == SP_ACT)then
                    ! ------------- aa block: full -------------
                    A0Block%type = phA0aa

                    ind = 0
                    do i = NI+1, NIA
                          do j = NI+1, i-1
                                if (abs(Occ(i)-Occ(j))/Occ(i) .lt. AuxData%ThrSelAct) cycle
                                ind = ind + 1
                          end do
                    end do
                    A0Block%NDim = ind

                    allocate(A0Block%IndN(2, A0Block%NDim))
                    A0Block%IndN = 0

                    ind = 0
                    do i = NI+1, NIA
                          do j = NI+1, i-1
                                if (abs(Occ(i)-Occ(j))/Occ(i) .lt. AuxData%ThrSelAct) cycle
                                ind = ind + 1
                                A0Block%IndN(1, ind) = i
                                A0Block%IndN(2, ind) = j
                                pos(i, j) = ind
                          end do
                    end do

                    allocate(A0Block%APlus(A0Block%NDim, A0Block%NDim))
                    allocate(A0Block%AMinus(A0Block%NDim, A0Block%NDim))
                    call CalcMem(A0Block%APlus, 'A0Block%APlus (aa)')
                    call CalcMem(A0Block%AMinus, 'A0Block%AMinus (aa)')
                    A0Block%APlus = zero
                    A0Block%AMinus = zero

                    allocate(A0Block%Eig(A0Block%NDim))
                    allocate(A0Block%vPlus(A0Block%NDim))
                    A0Block%Eig = zero
                    A0Block%vPlus = 0

              else if (Ip == SP_VIRT .and. Iq == SP_ACT)then
                    ! ------------- va block: mini-blocks per virtual -------------
                    A0Block%type = phA0va
                    Nbl = NV
                    Nbl0 = NIA + 1
                    Nbl1 = NBasis
                    A0Block%NDim = Nbl * NA

                    allocate(A0Block%IndN(2, A0Block%NDim))
                    allocate(A0Block%MiniBlocks(Nbl))

                    ind = 0
                    do i = Nbl0, Nbl1
                          k = i - Nbl0 + 1
                          allocate(A0Block%MiniBlocks(k)%ActList(NA))
                          A0Block%MiniBlocks(k)%ActList = 0
                          A0Block%MiniBlocks(k)%NDim = 0
                          inda = 0
                          do j = NI+1, NIA
                                ind = ind + 1
                                A0Block%IndN(1, ind) = i
                                A0Block%IndN(2, ind) = j
                                pos(i, j) = ind

                                if (abs(Occ(i)+Occ(j)) .gt. AuxData%ThrQVirt)then
                                      inda = inda + 1
                                      A0Block%MiniBlocks(k)%ActList(inda) = j
                                      A0Block%MiniBlocks(k)%NDim = inda
                                end if
                          end do
                    end do

              else if (Ip == SP_OCC .and. Iq == SP_ACT)then
                    ! ------------- oa block: mini-blocks per inactive -------------
                    A0Block%type = phA0oa
                    Nbl0 = NI0
                    Nbl1 = NI
                    Nbl = Nbl1 - Nbl0 + 1
                    A0Block%NDim = Nbl * NA

                    allocate(A0Block%IndN(2, A0Block%NDim))
                    ! indexed by the global inactive orbital, see header
                    allocate(A0Block%MiniBlocks(Nbl0:Nbl1))

                    ind = 0
                    do i = Nbl0, Nbl1
                          k = i
                          allocate(A0Block%MiniBlocks(k)%ActList(NA))
                          A0Block%MiniBlocks(k)%ActList = 0
                          A0Block%MiniBlocks(k)%NDim = 0
                          inda = 0
                          do j = NI+1, NIA
                                ind = ind + 1
                                ! pairs are ordered (act, inact), as in pp A0oa
                                A0Block%IndN(1, ind) = j
                                A0Block%IndN(2, ind) = i
                                pos(j, i) = ind

                                if (abs(Occ(i)+Occ(j)-Two) .gt. AuxData%ThrQInact)then
                                      inda = inda + 1
                                      A0Block%MiniBlocks(k)%ActList(inda) = j
                                      A0Block%MiniBlocks(k)%NDim = inda
                                end if
                          end do
                    end do
              end if

              ! mini-block matrices for the block-diagonal blocks
              if (A0Block%type == phA0va .or. A0Block%type == phA0oa)then
                    do k = lbound(A0Block%MiniBlocks, 1), ubound(A0Block%MiniBlocks, 1)
                          allocate(A0Block%MiniBlocks(k)%MiniAPlus(NA, NA))
                          allocate(A0Block%MiniBlocks(k)%MiniAMinus(NA, NA))
                          A0Block%MiniBlocks(k)%MiniAPlus = zero
                          A0Block%MiniBlocks(k)%MiniAMinus = zero
                          allocate(A0Block%MiniBlocks(k)%Eig(NA))
                          A0Block%MiniBlocks(k)%Eig = zero
                    end do
              end if

            end associate

      end subroutine InitA0Blocks_ph


      ! The 12 integral classes of the ph int loop (at most one virtual index).
      ! (IV|IV), (VA|IV), (VA|VA) and (VV|..) are not run: only the unbuilt vovo,
      ! vova, vava needed them. See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 3.
      subroutine parse_config_ph(loops, TabChol)
            type(TLoop), dimension(6), intent(out) :: loops
            integer, dimension(6), intent(out) :: TabChol
            integer :: i

            do i = 1, 6
                  loops(i)%run_main = .false.
                  loops(i)%run_inner = .false.
            end do

            loops(1)%run_main = .true.        ! (II|..
            loops(1)%run_inner(1) = .true.    ! (II|II)
            loops(2)%run_main = .true.        ! (IA|..
            loops(2)%run_inner(1) = .true.    ! (IA|II)
            loops(2)%run_inner(2) = .true.    ! (IA|IA)
            loops(3)%run_main = .true.        ! (AA|..
            loops(3)%run_inner(1) = .true.    ! (AA|II)
            loops(3)%run_inner(2) = .true.    ! (AA|IA)
            loops(3)%run_inner(3) = .true.    ! (AA|AA)
            loops(4)%run_main = .true.        ! (IV|..
            loops(4)%run_inner(1) = .true.    ! (IV|II)
            loops(4)%run_inner(2) = .true.    ! (IV|IA)
            loops(4)%run_inner(3) = .true.    ! (IV|AA)
            loops(5)%run_main = .true.        ! (VA|..
            loops(5)%run_inner(1) = .true.    ! (VA|II)
            loops(5)%run_inner(2) = .true.    ! (VA|IA)
            loops(5)%run_inner(3) = .true.    ! (VA|AA)

            ! RII, RIA, RAA, RIV, RVA; slot 6 (RVV) is never used in ph.
            TabChol = (/ 1, 1, 1, 1, 1, 0 /)

      end subroutine parse_config_ph


      ! =============================================================================
      ! Integral loop. Per slab V_axby = (a x|b y) the kernels fill AuxII/BuxII/
      ! AuxAA/BuxAA, Aux3X, the genuine A0aa terms and the delta-free part of the
      ! A1 blocks; the "together" loops follow the caseloop. A0aa and A1 are left
      ! unscaled and without delta terms. See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 3.
      ! =============================================================================
      subroutine THC_int_loop_ph(THCData, ACB, AuxII, AuxAA, BuxII, BuxAA, &
            Aux3X, Flags, AuxData, Occ, map, rdm2, rdm2_12, rdm2_13ph, pos, AcAlpha)
            use THC_Gammcor
            use thc_auto_p0_ph
            use thc_auto_p2_ph
            use thc_auto_p3_ph

            type(TTHCData), intent(in) :: THCData
            type(TAC0Block_ph), dimension(:), intent(inout) :: ACB
            double precision, dimension(:,:), intent(inout) :: AuxII, AuxAA
            double precision, dimension(:,:), intent(inout) :: BuxII, BuxAA
            double precision, dimension(:,:), intent(inout) :: Aux3X
            type(FlagsData), intent(in) :: Flags
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:), intent(in) :: Occ
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            ! Index-permuted copies of the same active 2-RDM, the ph analogue of
            ! pp's rdm2_pm_12_act / rdm2_pm_13_act trick: they turn every genuine
            ! contraction into a single ddot over NA**2 (see init_rdm_perm_ph).
            double precision, dimension(:,:,:,:), intent(in) :: rdm2_12, rdm2_13ph
            integer, dimension(:,:), intent(in) :: pos
            double precision, intent(in) :: AcAlpha

            double precision, dimension(:), allocatable :: V_axby
            double precision, dimension(:,:,:), allocatable :: RII, RIA, RAA
            double precision, dimension(:,:,:), allocatable :: RIV, RVA
            integer, dimension(6) :: TabChol
            type(TLoop), dimension(6) :: loops
            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            integer :: NI, NA, NV, NIA, NBasis
            integer :: a, b, ii
            integer :: e1, e2
            type (tclock) :: timer, timer0
            ! Accumulated wall times per integral class, as in pp THC_int_loop:
            ! t_<class>_int = building V_axby (real_aTb_x / real_aTv_x), t_<class>_ker =
            ! all kernels consuming it. tog_* are the genuine "together" loops after
            ! the caseloop. Printed only in a -DDEBUG build.
            double precision :: t_IIII_int, t_IIII_ker, t_IAII_int, t_IAII_ker
            double precision :: t_IAIA_int, t_IAIA_ker, t_AAII_int, t_AAII_ker
            double precision :: t_AAIA_int, t_AAIA_ker, t_AAAA_int, t_AAAA_ker
            double precision :: t_IVII_int, t_IVII_ker, t_IVIA_int, t_IVIA_ker
            double precision :: t_IVAA_int, t_IVAA_ker
            double precision :: t_VAII_int, t_VAII_ker, t_VAIA_int, t_VAIA_ker
            double precision :: t_VAAA_int, t_VAAA_ker
            double precision :: t_tog_AVAA_int, t_tog_AVAA_ker, t_tog_IAAA_int, t_tog_IAAA_ker
            double precision :: t_tog_IIAA_int, t_tog_IIAA_ker, t_tog_IVAA_int, t_tog_IVAA_ker
            double precision :: t_tog_VAAA_int, t_tog_VAAA_ker

            ! no associate for Zgk/Xga: they are unallocated in the exact-integral
            ! run, where no R block is built
            associate(NChol=>THCData%NChol)

              NA = AuxData%NA
              NI = AuxData%NI
              NV = AuxData%NV
              NIA = NI + NA
              NBasis = AuxData%NBasis

              call parse_config_ph(loops, TabChol)
              ! exact-integral run: every slab comes from THCData%TwoEl
              if (THCData%exact_ints) TabChol = 0

              t_IIII_int = zero
              t_IIII_ker = zero
              t_IAII_int = zero
              t_IAII_ker = zero
              t_IAIA_int = zero
              t_IAIA_ker = zero
              t_AAII_int = zero
              t_AAII_ker = zero
              t_AAIA_int = zero
              t_AAIA_ker = zero
              t_AAAA_int = zero
              t_AAAA_ker = zero
              t_IVII_int = zero
              t_IVII_ker = zero
              t_IVIA_int = zero
              t_IVIA_ker = zero
              t_IVAA_int = zero
              t_IVAA_ker = zero
              t_VAII_int = zero
              t_VAII_ker = zero
              t_VAIA_int = zero
              t_VAIA_ker = zero
              t_VAAA_int = zero
              t_VAAA_ker = zero
              t_tog_AVAA_int = zero
              t_tog_AVAA_ker = zero
              t_tog_IAAA_int = zero
              t_tog_IAAA_ker = zero
              t_tog_IIAA_int = zero
              t_tog_IIAA_ker = zero
              t_tog_IVAA_int = zero
              t_tog_IVAA_ker = zero
              t_tog_VAAA_int = zero
              t_tog_VAAA_ker = zero

              allocate(V_axby(max(NA, NI, NV)**2))

              ! No batching: every R block is allocated in full. RVV, the only one
              ! that could need it, is not used at all (see parse_config_ph).
              call print_section('ph THC integral transformation')

              if (TabChol(1) == 1) then
                    allocate(RII(NChol, NI, NI))
                    call CalcMem3(RII, 'RII')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RII, THCData%Xga(:,1:NI), THCData%Xga(:,1:NI), &
                          THCData%Zgk, NI, NI, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RII', timer, tdebug)
              end if
              if (TabChol(2) == 1) then
                    allocate(RIA(NChol, NA, NI))
                    call CalcMem3(RIA, 'RIA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,1:NI), &
                          THCData%Zgk, NA, NI, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RIA', timer, tdebug)
              end if
              if (TabChol(3) == 1) then
                    allocate(RAA(NChol, NA, NA))
                    call CalcMem3(RAA, 'RAA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RAA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,NI+1:NIA), &
                          THCData%Zgk, NA, NA, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RAA', timer, tdebug)
              end if

              if (TabChol(4) == 1) then
                    allocate(RIV(NChol, NV, NI))
                    call CalcMem3(RIV, 'RIV')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RIV, THCData%Xga(:,NIA+1:NBasis), THCData%Xga(:,1:NI), &
                          THCData%Zgk, NV, NI, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RIV', timer, tdebug)
              end if
              if (TabChol(5) == 1) then
                    allocate(RVA(NChol, NA, NV))
                    call CalcMem3(RVA, 'RVA')
                    call clock_start(timer)
                    call thc_gammcor_Rkab_2(RVA, THCData%Xga(:,NI+1:NIA), THCData%Xga(:,NIA+1:NBasis), &
                          THCData%Zgk, NA, NV, NChol, THCData%NTHC)
                    call tmsg('TIME FOR RVA', timer, tdebug)
              end if

              call clock_start(timer0)

            caseloop: do ii = 1, 6
                  if (loops(ii)%run_main)then
                        select case(ii)
                        case(loopII)
                              do a = 1, NI
                                    if (loops(ii)%run_inner(1)) then ! (II|II)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, 1, NI, b, 1, NI, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NI, RII(:,:,a), NChol, RII(:,:,b), NChol, NI, NI, NChol, AcAlpha)
                                                end if
                                                t_IIII_int = t_IIII_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IIII_qstt_qstt(BuxII, a, b, 1, NI, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IIII_qtst_qtst(AuxII, a, b, 1, NI, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                t_IIII_ker = t_IIII_ker + clock_readwall(timer)
                                          end do
                                    end if
                              end do
                        case(loopIA)
                              do a = 1, NI
                                    if (loops(ii)%run_inner(1)) then ! (IA|II)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, 1, NI, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                end if
                                                t_IAII_int = t_IAII_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IAII_qstt_sqtt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IAII_qtst_tsqt(AuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IAII_qtst_tqst(AuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                t_IAII_ker = t_IAII_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(2)) then ! (IA|IA)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                t_IAIA_int = t_IAIA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IAIA_qtst_tqts(AuxII, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call d2a_IA_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t4a_II_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d4a_AI_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d1b_IA_aoao_qrps(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d3b_AI_aoao_qrps(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t5a_IA_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t6a_AI_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_II_aoao_pqrs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_II_aoao_qrps(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call T2_aoao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D5_aoao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                t_IAIA_ker = t_IAIA_ker + clock_readwall(timer)
                                          end do
                                    end if
                              end do
                        case(loopAA)
                              do a = NI+1, NIA
                                    if (loops(ii)%run_inner(1)) then ! (AA|II)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, 1, NI, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                end if
                                                t_AAII_int = t_AAII_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_AAII_qstt_ttqs(BuxAA, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_AAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call t3b_II_aoao_prqs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, 1, NI, pos, Occ, V_axby)
                                                call t5b_IA_aoao_prqs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, 1, NI, pos, Occ, V_axby)
                                                call t6b_AI_aoao_prqs(ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, a, b, NI+1, NIA, 1, NI, pos, Occ, V_axby)
                                                t_AAII_ker = t_AAII_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(2)) then ! (AA|IA)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                t_AAIA_int = t_AAIA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_AAIA_qstt_ttqs(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_AAIA_qstt_ttsq(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_AAIA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_AAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call t4a_AI_aoaa_pqrs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d4a_AI_aoaa_pqrs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t3b_AI_aoaa_prqs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d3b_AI_aoaa_qrps(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t6a_AI_aoaa_pqrs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_AI_aoaa_pqrs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t6b_AI_aoaa_prqs(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_AI_aoaa_qrps(ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call T2_aoaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call T5_aoaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D2_aoaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D5_aoaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                ! Aux3X(inactive, active): the WMAT part of the aa-ao delta sector
                                                call P45_ph_IA(Aux3X, a, b, map, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                ! transposed partner aaao (rows aa, cols ao), symmetrize_A1_ph folds it into aoaa
                                                call d2a_IA_aaao_pqrs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t4a_IA_aaao_pqrs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d1b_IA_aaao_qrps(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t3b_IA_aaao_prqs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t5a_IA_aaao_pqrs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t5b_IA_aaao_prqs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_IA_aaao_pqrs(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_IA_aaao_qrps(ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call T2_aaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call T6_aaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D4_aaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D5_aaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                t_AAIA_ker = t_AAIA_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(3)) then ! (AA|AA)
                                          do b = NI+1, NIA
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                t_AAAA_int = t_AAAA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_AAAA_qstt_qstt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_AAAA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call P45_ph(Aux3X, a, b, map, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                call PhT12rp_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                call PhT12qs_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                call PhT56_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                t_AAAA_ker = t_AAAA_ker + clock_readwall(timer)
                                          end do
                                    end if
                              end do
                        case(loopIV)
                              do a = 1, NI
                                    if (loops(ii)%run_inner(1)) then ! (IV|II)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, 1, NI, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RII(:,:,b), NChol, NV, NI, NChol, AcAlpha)
                                                end if
                                                t_IVII_int = t_IVII_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IVII_qstt_qstt(BuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IVII_qstt_sqtt(BuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IVII_qtst_tsqt(AuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IVII_qtst_tqst(AuxII, a, b, NIA+1, NBasis, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                t_IVII_ker = t_IVII_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(2)) then ! (IV|IA)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RIA(:,:,b), NChol, NV, NA, NChol, AcAlpha)
                                                end if
                                                t_IVIA_int = t_IVIA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IVIA_qtst_tstq(AuxII, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IVIA_qtst_tqts(AuxII, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call d2a_IA_aovo_pqrs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call t4a_II_aovo_pqrs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call d1b_IA_aovo_qrps(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call t5a_IA_aovo_pqrs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_II_aovo_pqrs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_II_aovo_qrps(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                t_IVIA_ker = t_IVIA_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(3)) then ! (IV|AA)
                                          do b = NI+1, NIA
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NIA+1, NBasis, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NV, RIV(:,:,a), NChol, RAA(:,:,b-NI), NChol, NV, NA, NChol, AcAlpha)
                                                end if
                                                t_IVAA_int = t_IVAA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_IVAA_qstt_qstt(BuxAA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_IVAA_qstt_sqtt(BuxAA, a, b, NIA+1, NBasis, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call d2a_IA_aavo_pqrs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call t4a_IA_aavo_pqrs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call t5a_IA_aavo_pqrs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_IA_aavo_pqrs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_AI_aova_qrps(ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                ! transposed partner vaao (rows va, cols ao)
                                                call d6b_IA_vaao_qrps(ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, a, b, NIA+1, NBasis, NI+1, NIA, pos, Occ, V_axby)
                                                t_IVAA_ker = t_IVAA_ker + clock_readwall(timer)
                                          end do
                                    end if
                              end do
                        case(loopVA)
                              do a = NIA+1, NBasis
                                    if (loops(ii)%run_inner(1)) then ! (VA|II)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, 1, NI, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                end if
                                                t_VAII_int = t_VAII_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_VAII_qstt_sqtt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_VAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                call t3b_II_aovo_prqs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NI+1, NIA, 1, NI, pos, Occ, V_axby)
                                                call t5b_IA_aovo_prqs(ACB(phA1aovo)%APlus, ACB(phA1aovo)%AMinus, a, b, NI+1, NIA, 1, NI, pos, Occ, V_axby)
                                                t_VAII_ker = t_VAII_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(2)) then ! (VA|IA)
                                          do b = 1, NI
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                t_VAIA_int = t_VAIA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_VAIA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_VAIA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call t4a_AI_aova_pqrs(ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d1b_IA_aavo_qrps(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t3b_IA_aavo_prqs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t3b_AI_aova_prqs(ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t5b_IA_aavo_prqs(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_AI_aova_pqrs(ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6b_IA_aavo_qrps(ACB(phA1aavo)%APlus, ACB(phA1aavo)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call T5_aova(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D2_aova(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                ! transposed partner vaao (rows va, cols ao)
                                                call t4a_IA_vaao_pqrs(ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call t3b_IA_vaao_prqs(ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call d6a_IA_vaao_pqrs(ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, a, b, NI+1, NIA, NI+1, NIA, pos, Occ, V_axby)
                                                call T6_vaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D4_vaao(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                t_VAIA_ker = t_VAIA_ker + clock_readwall(timer)
                                          end do
                                    end if
                                    if (loops(ii)%run_inner(3)) then ! (VA|AA)
                                          do b = NI+1, NIA
                                                call clock_start(timer)
                                                if (THCData%exact_ints) then
                                                      call exact_V_axby(V_axby, THCData%TwoEl, a, NI+1, NIA, b, NI+1, NIA, AcAlpha)
                                                else
                                                      call real_aTb_x(V_axby, NA, RVA(:,:,a-NIA), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                end if
                                                t_VAAA_int = t_VAAA_int + clock_readwall(timer)
                                                call clock_start(timer)
                                                call p2_VAAA_qstt_sqtt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_VAAA_qstt_qstt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_VAAA_qtst_stqt(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call p2_VAAA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                call P45_ph_VA(Aux3X, a, b, map, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                call T4_aava(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call T5_aava(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D2_aava(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D6_aava(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                ! transposed partner vaaa (rows va, cols aa)
                                                call T4_vaaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call T6_vaaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D4_vaaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2_13ph, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                call D6_vaaa(a, b, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NA)
                                                t_VAAA_ker = t_VAAA_ker + clock_readwall(timer)
                                          end do
                                    end if
                              end do
                        end select
                  end if
            end do caseloop

            ! ---- genuine terms whose two summed indices sit in the SAME integral
            !      pair: (AV|AA). V_axby(t,u) is built for one fixed pair
            !      (e1,e2) at a time, exactly as pp does for the P3c*_gamma kernels.
            do e1 = NI+1, NIA
                  do e2 = NIA+1, NBasis
                        call clock_start(timer)
                        if (THCData%exact_ints) then
                              call exact_V_pair(V_axby, THCData%TwoEl, NI+1, NIA, NI+1, NIA, e1, e2, AcAlpha)
                        else
                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RVA(:, e1-NI, e2-NIA), NChol, NA**2, AcAlpha, Zero)
                        end if
                        t_tog_AVAA_int = t_tog_AVAA_int + clock_readwall(timer)
                        call clock_start(timer)
                        call D1_aava(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        call T3_vaaa(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        call D3_vaaa(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1vaaa)%APlus, ACB(phA1vaaa)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        t_tog_AVAA_ker = t_tog_AVAA_ker + clock_readwall(timer)
                  end do
            end do
            ! ---- genuine terms whose two summed indices sit in the SAME integral
            !      pair: (IA|AA). V_axby(t,u) is built for one fixed pair
            !      (e1,e2) at a time, exactly as pp does for the P3c*_gamma kernels.
            do e1 = 1, NI
                  do e2 = NI+1, NIA
                        call clock_start(timer)
                        if (THCData%exact_ints) then
                              call exact_V_pair(V_axby, THCData%TwoEl, NI+1, NIA, NI+1, NIA, e2, e1, AcAlpha)
                        else
                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RIA(:, e2-NI, e1), NChol, NA**2, AcAlpha, Zero)
                        end if
                        t_tog_IAAA_int = t_tog_IAAA_int + clock_readwall(timer)
                        call clock_start(timer)
                        call T1_aoaa(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        call D1_aoaa(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aoaa)%APlus, ACB(phA1aoaa)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        call D3_aaao(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        ! class AIAA: the same integral, generated with the pair as (e1,e2) = (act, inact)
                        call T1_aaao(e2, e1, NI+1, NIA, NI+1, NIA, ACB(phA1aaao)%APlus, ACB(phA1aaao)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        t_tog_IAAA_ker = t_tog_IAAA_ker + clock_readwall(timer)
                  end do
            end do
            ! ---- genuine terms whose two summed indices sit in the SAME integral
            !      pair: (II|AA). V_axby(t,u) is built for one fixed pair
            !      (e1,e2) at a time, exactly as pp does for the P3c*_gamma kernels.
            do e1 = 1, NI
                  do e2 = 1, NI
                        call clock_start(timer)
                        if (THCData%exact_ints) then
                              call exact_V_pair(V_axby, THCData%TwoEl, NI+1, NIA, NI+1, NIA, e2, e1, AcAlpha)
                        else
                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RII(:, e2, e1), NChol, NA**2, AcAlpha, Zero)
                        end if
                        t_tog_IIAA_int = t_tog_IIAA_int + clock_readwall(timer)
                        call clock_start(timer)
                        call T1_aoao(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aoao)%APlus, ACB(phA1aoao)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        t_tog_IIAA_ker = t_tog_IIAA_ker + clock_readwall(timer)
                  end do
            end do
            ! ---- genuine terms whose two summed indices sit in the SAME integral
            !      pair: (IV|AA). V_axby(t,u) is built for one fixed pair
            !      (e1,e2) at a time, exactly as pp does for the P3c*_gamma kernels.
            do e1 = 1, NI
                  do e2 = NIA+1, NBasis
                        call clock_start(timer)
                        if (THCData%exact_ints) then
                              call exact_V_pair(V_axby, THCData%TwoEl, NI+1, NIA, NI+1, NIA, e2, e1, AcAlpha)
                        else
                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RIV(:, e2-NIA, e1), NChol, NA**2, AcAlpha, Zero)
                        end if
                        t_tog_IVAA_int = t_tog_IVAA_int + clock_readwall(timer)
                        call clock_start(timer)
                        call D1_aova(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aova)%APlus, ACB(phA1aova)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        call D3_vaao(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1vaao)%APlus, ACB(phA1vaao)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        t_tog_IVAA_ker = t_tog_IVAA_ker + clock_readwall(timer)
                  end do
            end do
            ! ---- genuine terms whose two summed indices sit in the SAME integral
            !      pair: (VA|AA). V_axby(t,u) is built for one fixed pair
            !      (e1,e2) at a time, exactly as pp does for the P3c*_gamma kernels.
            do e1 = NIA+1, NBasis
                  do e2 = NI+1, NIA
                        call clock_start(timer)
                        if (THCData%exact_ints) then
                              call exact_V_pair(V_axby, THCData%TwoEl, NI+1, NIA, NI+1, NIA, e2, e1, AcAlpha)
                        else
                              call real_aTv_x(V_axby, RAA(:,:,:), NChol, RVA(:, e2-NI, e1-NIA), NChol, NA**2, AcAlpha, Zero)
                        end if
                        t_tog_VAAA_int = t_tog_VAAA_int + clock_readwall(timer)
                        call clock_start(timer)
                        call T3_aava(e1, e2, NI+1, NIA, NI+1, NIA, ACB(phA1aava)%APlus, ACB(phA1aava)%AMinus, map, pos, rdm2_12, NI+1, NIA, NI+1, NIA, V_axby, NA)
                        t_tog_VAAA_ker = t_tog_VAAA_ker + clock_readwall(timer)
                  end do
            end do

              call tmsg('ph int loop', timer0, tdebug)
!#ifdef DEBUG
call xmsg('t_IIII_int', t_IIII_int, mdebugg)
call xmsg('t_IIII_ker', t_IIII_ker, mdebugg)
call xmsg('t_IAII_int', t_IAII_int, mdebugg)
call xmsg('t_IAII_ker', t_IAII_ker, mdebugg)
call xmsg('t_IAIA_int', t_IAIA_int, mdebugg)
call xmsg('t_IAIA_ker', t_IAIA_ker, mdebugg)
call xmsg('t_AAII_int', t_AAII_int, mdebugg)
call xmsg('t_AAII_ker', t_AAII_ker, mdebugg)
call xmsg('t_AAIA_int', t_AAIA_int, mdebugg)
call xmsg('t_AAIA_ker', t_AAIA_ker, mdebugg)
call xmsg('t_AAAA_int', t_AAAA_int, mdebugg)
call xmsg('t_AAAA_ker', t_AAAA_ker, mdebugg)
call xmsg('t_IVII_int', t_IVII_int, mdebugg)
call xmsg('t_IVII_ker', t_IVII_ker, mdebugg)
call xmsg('t_IVIA_int', t_IVIA_int, mdebugg)
call xmsg('t_IVIA_ker', t_IVIA_ker, mdebugg)
call xmsg('t_IVAA_int', t_IVAA_int, mdebugg)
call xmsg('t_IVAA_ker', t_IVAA_ker, mdebugg)
call xmsg('t_VAII_int', t_VAII_int, mdebugg)
call xmsg('t_VAII_ker', t_VAII_ker, mdebugg)
call xmsg('t_VAIA_int', t_VAIA_int, mdebugg)
call xmsg('t_VAIA_ker', t_VAIA_ker, mdebugg)
call xmsg('t_VAAA_int', t_VAAA_int, mdebugg)
call xmsg('t_VAAA_ker', t_VAAA_ker, mdebugg)
call xmsg('t_tog_AVAA_int', t_tog_AVAA_int, mdebugg)
call xmsg('t_tog_AVAA_ker', t_tog_AVAA_ker, mdebugg)
call xmsg('t_tog_IAAA_int', t_tog_IAAA_int, mdebugg)
call xmsg('t_tog_IAAA_ker', t_tog_IAAA_ker, mdebugg)
call xmsg('t_tog_IIAA_int', t_tog_IIAA_int, mdebugg)
call xmsg('t_tog_IIAA_ker', t_tog_IIAA_ker, mdebugg)
call xmsg('t_tog_IVAA_int', t_tog_IVAA_int, mdebugg)
call xmsg('t_tog_IVAA_ker', t_tog_IVAA_ker, mdebugg)
call xmsg('t_tog_VAAA_int', t_tog_VAAA_int, mdebugg)
call xmsg('t_tog_VAAA_ker', t_tog_VAAA_ker, mdebugg)
!#endif

            end associate

      end subroutine THC_int_loop_ph


      ! Two index-permuted copies of rdm2, so every genuine contraction is one ddot
      ! over NA**2: rdm2_12(i,k,j,l) = rdm2_13ph(i,l,j,k) = rdm2(i,j,k,l).
      ! rdm2_13ph is NOT pp's 13 permutation. See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 5.
      subroutine init_rdm_perm_ph(rdm2, rdm2_12, rdm2_13ph, NA)
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            double precision, dimension(:,:,:,:), intent(out) :: rdm2_12, rdm2_13ph
            integer, intent(in) :: NA
            integer :: i, j, k, l

            do l = 1, NA
                  do k = 1, NA
                        do j = 1, NA
                              do i = 1, NA
                                    rdm2_12(i, k, j, l)   = rdm2(i, j, k, l)
                                    rdm2_13ph(i, l, j, k) = rdm2(i, j, k, l)
                              end do
                        end do
                  end do
            end do

      end subroutine init_rdm_perm_ph


      ! Dressed one-electron element for PHERPA_A1_THC: h + 2 BuxII - AuxII, plus
      ! 2 BuxAA - AuxAA unless the index whose occupation multiplies it (ind) is
      ! active. H0 must be AuxData%HNO0, never HNO. F(v,v') is not available.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 4.
      function Fel_ph(x, y, ind, H0, AuxII, BuxII, AuxAA, BuxAA) result(f)
            integer, intent(in) :: x, y, ind
            double precision, dimension(:,:), intent(in) :: H0
            double precision, dimension(:,:), intent(in) :: AuxII, BuxII, AuxAA, BuxAA
            double precision :: f
            integer, parameter :: SP_ACT = 1

            f = H0(x, y) + two * BuxII(x, y) - AuxII(x, y)
            if (ind .ne. SP_ACT) f = f + two * BuxAA(x, y) - AuxAA(x, y)

      end function Fel_ph


      ! =============================================================================
      ! Delta sector of one A1 block (types c, d, WMAT, T0) on top of the delta-free
      ! part from the int loop. Block-agnostic: everything follows from the spaces of
      ! p,q,r,s; the B groups enter AMinus with the opposite sign. H0 must be
      ! AuxData%HNO0 (bare h), not HNO. The block stays unscaled.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 3-4.
      ! =============================================================================
      subroutine PHERPA_A1_THC(AuxData, ACBlock, H0, AuxII, BuxII, AuxAA, BuxAA, Aux3X)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), intent(inout) :: ACBlock
            double precision, dimension(:,:), intent(in) :: H0
            double precision, dimension(:,:), intent(in) :: AuxII, BuxII
            double precision, dimension(:,:), intent(in) :: AuxAA, BuxAA
            double precision, dimension(:,:), intent(in) :: Aux3X

            integer, parameter :: SP_ACT = 1
            integer :: p, q, r, s, pq, rs
            double precision :: valP, valM, valA, valB

            associate(Occ=>AuxData%Occ, IndAux=>AuxData%IndAux, &
                  NDim1=>ACBlock%NDim1, NDim2=>ACBlock%NDim2, &
                  IndN1=>ACBlock%IndN1, IndN2=>ACBlock%IndN2)

              j_colloops: do pq = 1, NDim2
                    i_rowloops: do rs = 1, NDim1

                          r = IndN1(1, rs)
                          s = IndN1(2, rs)
                          p = IndN2(1, pq)
                          q = IndN2(2, pq)

                          valP = zero
                          valM = zero

                          ! ---- d_pr : T0 + T2c + T1d + WMAT_q,s, from A ----
                          if (p == r) then
                                valA = Occ(p) * Fel_ph(q, s, IndAux(p), H0, AuxII, BuxII, AuxAA, BuxAA) &
                                     - Occ(s) * Fel_ph(q, s, IndAux(s), H0, AuxII, BuxII, AuxAA, BuxAA)
                                if (IndAux(s) == SP_ACT) valA = valA - frac12 * Aux3X(q, s)
                                valP = valP + valA
                                valM = valM + valA
                          end if

                          ! ---- d_sq : T0 + T4c + T3d + WMAT_p,r, from A ----
                          if (s == q) then
                                valA = Occ(q) * Fel_ph(p, r, IndAux(q), H0, AuxII, BuxII, AuxAA, BuxAA) &
                                     - Occ(r) * Fel_ph(p, r, IndAux(r), H0, AuxII, BuxII, AuxAA, BuxAA)
                                if (IndAux(r) == SP_ACT) valA = valA - frac12 * Aux3X(p, r)
                                valP = valP + valA
                                valM = valM + valA
                          end if

                          ! ---- d_qr : D0 + D2c + D1d + WMAT_p,s, from B ----
                          if (q == r) then
                                valB = Occ(q) * Fel_ph(p, s, IndAux(q), H0, AuxII, BuxII, AuxAA, BuxAA) &
                                     - Occ(s) * Fel_ph(p, s, IndAux(s), H0, AuxII, BuxII, AuxAA, BuxAA)
                                if (IndAux(s) == SP_ACT) valB = valB - frac12 * Aux3X(p, s)
                                valP = valP + valB
                                valM = valM - valB
                          end if

                          ! ---- d_sp : D0 + D4c + D3d + WMAT_q,r, from B ----
                          if (s == p) then
                                valB = Occ(p) * Fel_ph(q, r, IndAux(p), H0, AuxII, BuxII, AuxAA, BuxAA) &
                                     - Occ(r) * Fel_ph(q, r, IndAux(r), H0, AuxII, BuxII, AuxAA, BuxAA)
                                if (IndAux(r) == SP_ACT) valB = valB - frac12 * Aux3X(q, r)
                                valP = valP + valB
                                valM = valM - valB
                          end if

                          ACBlock%APlus(rs, pq)  = ACBlock%APlus(rs, pq)  + valP
                          ACBlock%AMinus(rs, pq) = ACBlock%AMinus(rs, pq) + valM

                    end do i_rowloops
              end do j_colloops

            end associate

      end subroutine PHERPA_A1_THC


      ! A_RC <- 1/2 (A_RC + A_CR^T) for aova<-vaao, aoaa<-aaao, aava<-vaaa (APlus and
      ! AMinus), then frees the partners. Needed because the energy doubles one
      ! orientation and A^(1) is symmetric only for a zero state-specific orbital
      ! gradient (not SA-CASSCF). See .secret_sauce/.skills/ph-ac0/code-thc-energy.md Sec. 6.
      subroutine symmetrize_A1_ph(ACB)
            type(TAC0Block_ph), dimension(:), intent(inout) :: ACB

            integer, dimension(2, 3), parameter :: pairs = reshape( &
                  [phA1aova, phA1vaao, phA1aoaa, phA1aaao, phA1aava, phA1vaaa], [2, 3])
            integer :: k, rc, cr, i, j
            double precision :: dp, dm, sp, sm

            call print_section('A1 symmetrization  A_RC <- 1/2 (A_RC + A_CR^T)')
            write(*,'(2X,A16,4A16)') 'block <- partner', '||A+||', '||A+ - A+^T||', &
                  '||A-||', '||A- - A-^T||'
            write(*,'(2X,A)') repeat('-', 82)

            do k = 1, 3
                  rc = pairs(1, k)
                  cr = pairs(2, k)
                  if (ACB(rc)%NDim1 /= ACB(cr)%NDim2 .or. ACB(rc)%NDim2 /= ACB(cr)%NDim1) then
                        call msg('symmetrize_A1_ph: '//trim(BlockName_ph(rc))//' and ' &
                              //trim(BlockName_ph(cr))//' are not transposes')
                        stop
                  end if

                  dp = zero
                  dm = zero
                  sp = zero
                  sm = zero
                  !$omp parallel do default(shared) private(i, j) reduction(+:dp, dm, sp, sm)
                  do j = 1, ACB(rc)%NDim2
                        do i = 1, ACB(rc)%NDim1
                              dp = dp + (ACB(rc)%APlus(i, j) - ACB(cr)%APlus(j, i))**2
                              dm = dm + (ACB(rc)%AMinus(i, j) - ACB(cr)%AMinus(j, i))**2
                              sp = sp + ACB(rc)%APlus(i, j)**2
                              sm = sm + ACB(rc)%AMinus(i, j)**2
                              ACB(rc)%APlus(i, j) = frac12 * (ACB(rc)%APlus(i, j) + ACB(cr)%APlus(j, i))
                              ACB(rc)%AMinus(i, j) = frac12 * (ACB(rc)%AMinus(i, j) + ACB(cr)%AMinus(j, i))
                        end do
                  end do
                  !$omp end parallel do

                  write(*,'(2X,A16,4ES16.4)') trim(BlockName_ph(rc))//' <- '//trim(BlockName_ph(cr)), &
                        sqrt(sp), sqrt(dp), sqrt(sm), sqrt(dm)

                  deallocate(ACB(cr)%APlus, ACB(cr)%AMinus)
            end do

      end subroutine symmetrize_A1_ph


      ! =============================================================================
      ! A0 assembly + diagonalization, one routine per block (vo, aa, va, oa).
      ! Formulas: .secret_sauce/.skills/ph-ac0/theory-a0-blocks.md;
      ! what each routine leaves in the block: .secret_sauce/.skills/ph-ac0/code-thc-a0.md Sec. 3-4.
      ! =============================================================================


      ! vo block: diagonal, A = d_pr d_sq (f_vv - f_ii). The c factors are
      ! exactly 1 here (c_virt = 0, c_inact = 1), so APlus = AMinus = f_v - f_i
      ! and omega = f_v - f_i.
      subroutine PHERPA0_THC_VO(AuxData, ACB_vo, THCData)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), intent(inout) :: ACB_vo
            type(TTHCData), intent(in) :: THCData

            integer :: k, v, i

            do k = 1, ACB_vo%NDim
                  v = ACB_vo%IndN(1, k)
                  i = ACB_vo%IndN(2, k)
                  ACB_vo%APlus(k, 1) = THCData%fvw(v - AuxData%NIA) - THCData%fij(i)
                  ACB_vo%AMinus(k, 1) = ACB_vo%APlus(k, 1)
                  ACB_vo%Eig(k) = ACB_vo%APlus(k, 1)
                  ACB_vo%vPlus(k) = 1
            end do

#ifdef DEBUG
            call print_block_eigs_ph('A0vo', ACB_vo%Eig, ACB_vo%NDim)
#endif

      end subroutine PHERPA0_THC_VO


      ! aa block: add the delta terms (T0 with h^core from HNO, TW1/TW2 with
      ! the genuine WMAT matrix Aux3X) onto the two-electron genuine part
      ! already accumulated in APlus/AMinus by the int loop, with the same
      ! p/q swap bookkeeping as pherpa_incore_spinres; then scale and solve.
      subroutine PHERPA0_THC_AA(AuxData, ACB_aa, HNO, Aux3X)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), intent(inout) :: ACB_aa
            double precision, dimension(:,:), intent(in) :: HNO, Aux3X

            double precision, dimension(:), allocatable :: c
            integer :: i, j, p, q, r, s, pp, qq, nneg
            double precision :: val, den

            associate(Occ=>AuxData%Occ, NBasis=>AuxData%NBasis, nd=>ACB_aa%NDim)

              if (nd == 0) return

              allocate(c(NBasis))
              do i = 1, NBasis
                    c(i) = sqrt(Occ(i))
                    if (Occ(i) .lt. frac12) c(i) = -c(i)
              end do

              do i = 1, nd
                    r = ACB_aa%IndN(1, i)
                    s = ACB_aa%IndN(2, i)
                    do j = 1, nd
                          pp = ACB_aa%IndN(1, j)
                          qq = ACB_aa%IndN(2, j)

                          do p = qq, pp, pp - qq
                                do q = qq, pp, pp - qq
                                      if (p == q) cycle

                                      val = zero
                                      if (p == r) then
                                            val = val + (Occ(p) - Occ(s)) * HNO(q, s)
                                            val = val - frac12 * Aux3X(q, s)
                                      end if
                                      if (s == q) then
                                            val = val + (Occ(q) - Occ(r)) * HNO(p, r)
                                            val = val - frac12 * Aux3X(p, r)
                                      end if

                                      if (p .gt. q) then
                                            ACB_aa%APlus(i, j) = ACB_aa%APlus(i, j) + val
                                            ACB_aa%AMinus(i, j) = ACB_aa%AMinus(i, j) + val
                                      else
                                            ACB_aa%APlus(i, j) = ACB_aa%APlus(i, j) + val
                                            ACB_aa%AMinus(i, j) = ACB_aa%AMinus(i, j) - val
                                      end if
                                end do
                          end do
                    end do
              end do

              ! metric scaling, canonical pairs
              do i = 1, nd
                    r = ACB_aa%IndN(1, i)
                    s = ACB_aa%IndN(2, i)
                    do j = 1, nd
                          p = ACB_aa%IndN(1, j)
                          q = ACB_aa%IndN(2, j)
                          den = (c(p) + c(q)) * (c(r) + c(s))
                          if (den .ne. zero) ACB_aa%APlus(i, j) = ACB_aa%APlus(i, j) / den
                          den = (c(p) - c(q)) * (c(r) - c(s))
                          if (den .ne. zero) ACB_aa%AMinus(i, j) = ACB_aa%AMinus(i, j) / den
                    end do
              end do

              call pherpa_symm_block_ph(ACB_aa%APlus, ACB_aa%AMinus, ACB_aa%Eig, nd, nneg, AuxData%nst(1) > 1)
#ifdef DEBUG
              call print_block_eigs_ph('A0aa', ACB_aa%Eig, nd)
#endif
              if (nneg .gt. 0) call imsg('A0aa negative/zero omega^2 modes', nneg)

            end associate

      end subroutine PHERPA0_THC_AA


      ! va block: one NA x NA eigenproblem per virtual v (p = r = v).
      ! A(s,q) = -n_s HNO(q,s) - 1/2 Aux3X(q,s) + d_sq n_q f_vv, B = 0.
      ! Scaling: (c_v +- c_q)(c_v +- c_s) = c_q c_s in both channels (c_v = 0).
      subroutine PHERPA0_THC_VA_Eigs(AuxData, ACB_va, THCData, HNO, Aux3X)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), intent(inout) :: ACB_va
            type(TTHCData), intent(in) :: THCData
            double precision, dimension(:,:), intent(in) :: HNO, Aux3X

            double precision, dimension(:,:), allocatable :: workp, workm
            double precision, dimension(:), allocatable :: eigw, allgather, c
            integer :: k, v, ii, jj, s, q, nd, nneg, ngather, i
            double precision :: val, den, fvv

            associate(Occ=>AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, &
                  NIA=>AuxData%NIA, NV=>AuxData%NV, NBasis=>AuxData%NBasis)

              allocate(c(NBasis))
              do i = 1, NBasis
                    c(i) = sqrt(Occ(i))
                    if (Occ(i) .lt. frac12) c(i) = -c(i)
              end do

              allocate(allgather(NV * NA))
              ngather = 0

              do k = 1, NV
                    v = NIA + k
                    fvv = THCData%fvw(k)
                    nd = ACB_va%MiniBlocks(k)%NDim
                    if (nd == 0) cycle

                    allocate(workp(nd, nd))
                    allocate(workm(nd, nd))
                    allocate(eigw(nd))

                    do jj = 1, nd
                          q = ACB_va%MiniBlocks(k)%ActList(jj)
                          do ii = 1, nd
                                s = ACB_va%MiniBlocks(k)%ActList(ii)
                                val = -Occ(s) * HNO(q, s) - frac12 * Aux3X(q, s)
                                if (s == q) val = val + Occ(q) * fvv
                                den = c(q) * c(s)
                                workp(ii, jj) = val / den
                                workm(ii, jj) = val / den
                          end do
                    end do

                    call pherpa_symm_block_ph(workp, workm, eigw, nd, nneg, AuxData%nst(1) > 1)

                    ! scatter back, active-local indexing (p = global - NI)
                    do jj = 1, nd
                          q = ACB_va%MiniBlocks(k)%ActList(jj)
                          ACB_va%MiniBlocks(k)%Eig(q - NI) = eigw(jj)
                          do ii = 1, nd
                                s = ACB_va%MiniBlocks(k)%ActList(ii)
                                ACB_va%MiniBlocks(k)%MiniAPlus(s - NI, q - NI) = workp(ii, jj)
                                ACB_va%MiniBlocks(k)%MiniAMinus(s - NI, q - NI) = workm(ii, jj)
                          end do
                    end do

                    allgather(ngather+1:ngather+nd) = eigw(1:nd)
                    ngather = ngather + nd

                    deallocate(workp, workm, eigw)
              end do

#ifdef DEBUG
              call print_block_eigs_ph('A0va', allgather, ngather)
#endif

            end associate

      end subroutine PHERPA0_THC_VA_Eigs


      ! oa block: one NA x NA eigenproblem per correlated inactive i (s = q = i),
      ! i = NI0..NI; frozen inactives have no mini-block (InitA0Blocks_ph).
      ! A(r,p) = (1 - n_r) HNO(p,r) + 2 BuxAA(p,r) - AuxAA(p,r)
      !          - 1/2 Aux3X(p,r) + d_pr (n_p - 1) f_ii,   B = 0.
      ! Scaling: (c_p + 1)(c_r + 1) plus channel, (c_p - 1)(c_r - 1) minus.
      subroutine PHERPA0_THC_OA_Eigs(AuxData, ACB_oa, THCData, HNO, Aux3X, AuxAA, BuxAA)
            type(TACppData), intent(in) :: AuxData
            type(TAC0Block_ph), intent(inout) :: ACB_oa
            type(TTHCData), intent(in) :: THCData
            double precision, dimension(:,:), intent(in) :: HNO, Aux3X
            double precision, dimension(:,:), intent(in) :: AuxAA, BuxAA

            double precision, dimension(:,:), allocatable :: workp, workm
            double precision, dimension(:), allocatable :: eigw, allgather, c
            integer :: k, ii, jj, r, p, nd, nneg, ngather, i
            double precision :: val, den, fii

            associate(Occ=>AuxData%Occ, NI=>AuxData%NI, NA=>AuxData%NA, &
                  NBasis=>AuxData%NBasis)

              allocate(c(NBasis))
              do i = 1, NBasis
                    c(i) = sqrt(Occ(i))
                    if (Occ(i) .lt. frac12) c(i) = -c(i)
              end do

              allocate(allgather(NI * NA))
              ngather = 0

              do k = AuxData%NFrozen + 1, NI
                    fii = THCData%fij(k)
                    nd = ACB_oa%MiniBlocks(k)%NDim
                    if (nd == 0) cycle

                    allocate(workp(nd, nd))
                    allocate(workm(nd, nd))
                    allocate(eigw(nd))

                    do jj = 1, nd
                          p = ACB_oa%MiniBlocks(k)%ActList(jj)
                          do ii = 1, nd
                                r = ACB_oa%MiniBlocks(k)%ActList(ii)
                                val = (One - Occ(r)) * HNO(p, r) &
                                      + Two * BuxAA(p, r) - AuxAA(p, r) &
                                      - frac12 * Aux3X(p, r)
                                if (p == r) val = val + (Occ(p) - One) * fii
                                den = (c(p) + One) * (c(r) + One)
                                workp(ii, jj) = val / den
                                den = (c(p) - One) * (c(r) - One)
                                workm(ii, jj) = val / den
                          end do
                    end do

                    call pherpa_symm_block_ph(workp, workm, eigw, nd, nneg, AuxData%nst(1) > 1)

                    do jj = 1, nd
                          p = ACB_oa%MiniBlocks(k)%ActList(jj)
                          ACB_oa%MiniBlocks(k)%Eig(p - NI) = eigw(jj)
                          do ii = 1, nd
                                r = ACB_oa%MiniBlocks(k)%ActList(ii)
                                ACB_oa%MiniBlocks(k)%MiniAPlus(r - NI, p - NI) = workp(ii, jj)
                                ACB_oa%MiniBlocks(k)%MiniAMinus(r - NI, p - NI) = workm(ii, jj)
                          end do
                    end do

                    allgather(ngather+1:ngather+nd) = eigw(1:nd)
                    ngather = ngather + nd

                    deallocate(workp, workm, eigw)
              end do

#ifdef DEBUG
              call print_block_eigs_ph('A0oa', allgather, ngather)
#endif

            end associate

      end subroutine PHERPA0_THC_OA_Eigs


      ! Ground-state solver for one block (port of pherpa_symm_spinres):
      ! A_+^1/2 A_- A_+^1/2 Z = omega^2 Z, Yt = A_+^1/2 Z, (2/omega) Yt^T A_- Yt = 1.
      ! Out: APlus = Yt, AMinus = symmetrized A_-, Eig (<= 0 flags an unusable mode).
      ! excited -> pherpa_nonsymm_block_ph. See .secret_sauce/.skills/ph-ac0/code-thc-a0.md Sec. 4.
      subroutine pherpa_symm_block_ph(APlus, AMinus, Eig, nd, nneg, excited)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            double precision, dimension(:), intent(out) :: Eig
            integer, intent(in) :: nd
            integer, intent(out) :: nneg
            logical, intent(in) :: excited

            double precision, dimension(:,:), allocatable :: work_a, work_b
            double precision, dimension(:), allocatable :: lambda, vec_tmp
            double precision :: omega_sq, dot_prod, norm_fact
            integer :: j, mode
            double precision, parameter :: small = 1.d-10

            nneg = 0
            if (nd == 0) return

            if (excited) then
                  call pherpa_nonsymm_block_ph(APlus, AMinus, Eig, nd, nneg)
                  return
            end if

            allocate(work_a(nd, nd))
            allocate(work_b(nd, nd))
            allocate(lambda(nd))
            allocate(vec_tmp(nd))

            ! the raw blocks are symmetric only up to the CASSCF gradient
            ! (Brillouin)
            APlus = frac12 * (APlus + transpose(APlus))
            AMinus = frac12 * (AMinus + transpose(AMinus))

            ! A_+^(1/2)
            work_a = APlus
            call symmetric_eigenproblem(lambda, work_a, nd, .true.)
            do j = 1, nd
                  work_b(:, j) = work_a(:, j) * sqrt(abs(lambda(j)))
            end do
            call real_abT(APlus, work_b, work_a)

            ! M = A_+^(1/2) A_- A_+^(1/2)
            work_a = zero
            work_b = zero
            call real_abaT(work_a, APlus, AMinus, work_b)
            call symmetric_eigenproblem(Eig, work_a, nd, .true.)

            ! Ytilde = A_+^(1/2) Z
            call real_ab(work_b, APlus, work_a)
            APlus = work_b

            do mode = 1, nd
                  omega_sq = Eig(mode)
                  if (omega_sq .gt. small) then
                        Eig(mode) = sqrt(omega_sq)
                        call real_Av(vec_tmp, AMinus, APlus(:, mode))
                        call real_vw_x(dot_prod, APlus(:, mode), vec_tmp, nd)
                        norm_fact = (two / Eig(mode)) * dot_prod
                        if (norm_fact .gt. zero) then
                              norm_fact = one / sqrt(norm_fact)
                        else
                              call msg('pherpa_symm_block_ph: zero norm')
                              norm_fact = zero
                        end if
                        call real_scal(APlus(:, mode), norm_fact)
                  else
                        Eig(mode) = -sqrt(abs(omega_sq))
                        nneg = nneg + 1
                  end if
            end do

            deallocate(work_a, work_b, lambda, vec_tmp)

      end subroutine pherpa_symm_block_ph


      ! Excited-state solver for one block (indefinite A_+): A_+ A_- Y = omega^2 Y,
      ! degenerate groups made A_- orthogonal, de-excitations and complex roots
      ! flagged by Eig <= 0. Same in/out contract as pherpa_symm_block_ph.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a0.md Sec. 4.
      subroutine pherpa_nonsymm_block_ph(APlus, AMinus, Eig, nd, nneg)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            double precision, dimension(:), intent(out) :: Eig
            integer, intent(in) :: nd
            integer, intent(out) :: nneg

            double precision, dimension(:,:), allocatable :: work_m, vl, vr
            double precision, dimension(:,:), allocatable :: yg, ayg, sg
            double precision, dimension(:), allocatable :: wr, wi, sval
            integer, dimension(:), allocatable :: order
            logical, dimension(:), allocatable :: is_real
            double precision :: omega
            integer :: i, j, k, g0, g1, ng, tmp
            double precision, parameter :: small = 1.d-10
            double precision, parameter :: small_imag = 1.d-8
            double precision, parameter :: thr_deg = 1.d-8

            nneg = 0
            if (nd == 0) return

            ! the raw blocks are symmetric only up to the CASSCF gradient
            APlus = frac12 * (APlus + transpose(APlus))
            AMinus = frac12 * (AMinus + transpose(AMinus))

            allocate(work_m(nd, nd), vl(nd, nd), vr(nd, nd))
            allocate(wr(nd), wi(nd), order(nd), is_real(nd))

            ! A_+ A_- Y = omega^2 Y
            call real_ab(work_m, APlus, AMinus)
            call nonsymmetric_eigenproblem(wr, wi, vl, vr, work_m)

            do k = 1, nd
                  is_real(k) = abs(wi(k)) <= small_imag * max(one, abs(wr(k)))
            end do

            ! ascending omega^2, so that degenerate roots are adjacent
            order = [(k, k = 1, nd)]
            do i = 2, nd
                  tmp = order(i)
                  j = i - 1
                  do while (j >= 1)
                        if (wr(order(j)) <= wr(tmp)) exit
                        order(j+1) = order(j)
                        j = j - 1
                  end do
                  order(j+1) = tmp
            end do

            APlus = zero
            Eig = zero
            g0 = 1
            do while (g0 <= nd)
                  ! group g0..g1 (sorted positions) of degenerate real roots
                  k = order(g0)
                  g1 = g0
                  if (is_real(k)) then
                        do while (g1 < nd)
                              if (.not. is_real(order(g1+1))) exit
                              if (abs(wr(order(g1+1)) - wr(k)) > thr_deg) exit
                              g1 = g1 + 1
                        end do
                  end if
                  ng = g1 - g0 + 1

                  if (.not. is_real(k)) then
                        ! complex omega^2: dropped
                        Eig(g0) = zero
                        nneg = nneg + 1
                  else if (wr(k) <= small) then
                        ! negative or zero omega^2
                        do i = g0, g1
                              Eig(i) = -sqrt(abs(wr(order(i))))
                              APlus(:, i) = vr(:, order(i))
                        end do
                        nneg = nneg + ng
                  else
                        omega = sqrt(wr(k))
                        allocate(yg(nd, ng), ayg(nd, ng), sg(ng, ng), sval(ng))
                        do i = 1, ng
                              yg(:, i) = vr(:, order(g0+i-1))
                        end do
                        ! S = 2/omega Yg^T A_- Yg, then Yg U
                        call real_ab(ayg, AMinus, yg)
                        call real_aTb(sg, yg, ayg)
                        sg = (two / omega) * frac12 * (sg + transpose(sg))
                        call symmetric_eigenproblem(sval, sg, ng, .true.)
                        call real_ab(ayg, yg, sg)
                        do i = 1, ng
                              j = g0 + i - 1
                              if (sval(i) > small) then
                                    Eig(j) = omega
                                    APlus(:, j) = ayg(:, i) / sqrt(sval(i))
                              else if (sval(i) < -small) then
                                    ! negative norm: de-excitation
                                    Eig(j) = -omega
                                    APlus(:, j) = ayg(:, i) / sqrt(-sval(i))
                                    nneg = nneg + 1
                              else
                                    Eig(j) = zero
                                    nneg = nneg + 1
                              end if
                        end do
                        deallocate(yg, ayg, sg, sval)
                  end if

                  g0 = g1 + 1
            end do

            deallocate(work_m, vl, vr, wr, wi, order, is_real)

      end subroutine pherpa_nonsymm_block_ph


      ! Print the sorted eigenvalues of one A0 block (a.u. and eV), lowest 30.
      subroutine print_block_eigs_ph(name, Eig, ndim)
            character(*), intent(in) :: name
            double precision, dimension(:), intent(in) :: Eig
            integer, intent(in) :: ndim

            double precision, dimension(:), allocatable :: work
            integer, dimension(:), allocatable :: dy
            integer :: i

            if (ndim == 0) then
                  write(*, '(2X,A,A)') name, ': empty block'
                  return
            end if

            allocate(work(ndim))
            allocate(dy(ndim))
            work = Eig(1:ndim)
            do i = 1, ndim
                  dy(i) = i
            end do
            call dsort(work, dy, ndim)

            write(*, '(/,2X,A,A,I8,A)') name, ' eigenvalues (', ndim, ' total; a.u., eV)'
            do i = 1, min(30, ndim)
                  write(*, '(I6, 2F16.6)') i, work(i), 27.2114d0 * work(i)
            end do
            deallocate(work, dy)

      end subroutine print_block_eigs_ph


      ! Genuine WMAT matrix Aux3X(q,r) = 2 G(q,r), q active, class (AA|AA); the
      ! consumers apply -1/2. See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 4.
      subroutine P45_ph(Aux3X, a, b, map, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: Aux3X
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: r, x, y
            double precision :: val

            do r = NI+1, NI+NA
                  do y = y0, y1
                        val = zero
                        do x = x0, x1
                              val = val + rdm2(map(r), map(x), map(b), map(a)) * V_axby(x, y)
                        end do
                        Aux3X(y, r) = Aux3X(y, r) + val
                  end do
                  do x = x0, x1
                        val = zero
                        do y = y0, y1
                              val = val + rdm2(map(r), map(y), map(a), map(b)) * V_axby(x, y)
                        end do
                        Aux3X(x, r) = Aux3X(x, r) + val
                  end do
            end do
      end subroutine P45_ph


      ! Virtual rows of Aux3X, class (VA|AA); needed by ao-vo and aa-va.
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 4.
      subroutine P45_ph_VA(Aux3X, a, b, map, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: Aux3X
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: r, x, y
            double precision :: val

            do r = NI+1, NI+NA
                  val = zero
                  do y = y0, y1
                        do x = x0, x1
                              val = val + rdm2(map(r), map(y), map(x), map(b)) * V_axby(x, y)
                        end do
                  end do
                  Aux3X(a, r) = Aux3X(a, r) + two * val
            end do
      end subroutine P45_ph_VA


      ! Inactive rows of Aux3X, class (AA|IA); needed only by aa-ao (phA1aaao).
      ! See .secret_sauce/.skills/ph-ac0/code-thc-a1.md Sec. 4.
      subroutine P45_ph_IA(Aux3X, a, b, map, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: Aux3X
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: r, x, y
            double precision :: val

            do r = NI+1, NI+NA
                  val = zero
                  do y = y0, y1
                        do x = x0, x1
                              val = val + rdm2(map(r), map(x), map(y), map(a)) * V_axby(x, y)
                        end do
                  end do
                  Aux3X(b, r) = Aux3X(b, r) + two * val
            end do
      end subroutine P45_ph_IA


      ! Route one two-electron contribution arspq into APlus/AMinus of a block.
      ! rs is the (already looked-up) row position; the column pair (p,q) may
      ! come in either order: p>q is the A side (+ into both), p<q is the
      ! B side (B_rspq = A_rsqp: + into APlus, - into AMinus) - the same
      ! bookkeeping as the p/q swap loop of pherpa_incore_spinres. pos entries
      ! are zero for non-canonical or filtered pairs, so those are skipped.
      subroutine accum_A_ph(APlus, AMinus, rs, p, q, pos, val)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            integer, intent(in) :: rs, p, q
            integer, dimension(:,:), intent(in) :: pos
            double precision, intent(in) :: val
            integer :: pq

            if (rs == 0) return
            if (p == q) return
            if (p .gt. q) then
                  pq = pos(p, q)
                  if (pq == 0) return
                  APlus(rs, pq) = APlus(rs, pq) + val
                  AMinus(rs, pq) = AMinus(rs, pq) + val
            else
                  pq = pos(q, p)
                  if (pq == 0) return
                  APlus(rs, pq) = APlus(rs, pq) + val
                  AMinus(rs, pq) = AMinus(rs, pq) - val
            end if
      end subroutine accum_A_ph


      ! Genuine branches of func_T1(r,p,s,q) and func_T2(r,p,s,q) of
      ! pherpa_incore_spinres (tex terms T2 and T1). Row second index s = a,
      ! Gamma second index t = b; per element (r, s=a | p, q):
      !   func_T1: + sum_u rdm2(r,b,u,p) (su|bq)   u = x summed, q = y free
      !   func_T2: + sum_u rdm2(r,b,p,u) (sq|bu)   q = x free, u = y summed
      subroutine PhT12rp_gamma(APlus, AMinus, a, b, map, pos, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: r, p, q, x, y, rs
            double precision :: val

            do r = x0, x1
                  rs = pos(r, a)
                  if (rs == 0) cycle
                  do p = x0, x1
                        do q = x0, x1
                              if (p == q) cycle
                              val = zero
                              do x = x0, x1
                                    val = val + rdm2(map(r), map(b), map(x), map(p)) * V_axby(x, q)
                              end do
                              do y = y0, y1
                                    val = val + rdm2(map(r), map(b), map(p), map(y)) * V_axby(q, y)
                              end do
                              call accum_A_ph(APlus, AMinus, rs, p, q, pos, val)
                        end do
                  end do
            end do
      end subroutine PhT12rp_gamma


      ! Genuine branches of func_T1(s,q,r,p) and func_T2(s,q,r,p) of
      ! pherpa_incore_spinres (tex terms T4 and T3). Row first index r = a,
      ! Gamma second index t = b; per element (r=a, s | p, q):
      !   func_T1(s,q,r,p): + sum_u rdm2(s,b,u,q) (ru|bp)   u = x summed, p = y free
      !   func_T2(s,q,r,p): + sum_u rdm2(s,b,q,u) (rp|bu)   p = x free, u = y summed
      subroutine PhT12qs_gamma(APlus, AMinus, a, b, map, pos, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: s, p, q, x, y, rs
            double precision :: val

            do s = x0, x1
                  rs = pos(a, s)
                  if (rs == 0) cycle
                  do q = x0, x1
                        do p = x0, x1
                              if (p == q) cycle
                              val = zero
                              do x = x0, x1
                                    val = val + rdm2(map(s), map(b), map(x), map(q)) * V_axby(x, p)
                              end do
                              do y = y0, y1
                                    val = val + rdm2(map(s), map(b), map(q), map(y)) * V_axby(p, y)
                              end do
                              call accum_A_ph(APlus, AMinus, rs, p, q, pos, val)
                        end do
                  end do
            end do
      end subroutine PhT12qs_gamma


      ! Genuine func_T34 of pherpa_incore_spinres, entering arspq twice
      ! (tex terms T6 and T5). Both Gamma summation indices are the loop
      ! indices: t = a, u = b; per fixed (a,b) there is no inner sum:
      !   T6 = -func_T34(s,p,r,q): - rdm2(a,b,s,p) (ar|bq)   r = x, q = y
      !   T5 = -func_T34(q,r,p,s): - rdm2(a,b,q,r) (ap|bs)   p = x, s = y
      subroutine PhT56_gamma(APlus, AMinus, a, b, map, pos, rdm2, x0, x1, y0, y1, V_axby, NI, NA)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            integer, intent(in) :: a, b, NI, NA
            integer, dimension(:), intent(in) :: map
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            integer, intent(in) :: x0, x1, y0, y1
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: r, s, p, q, x, y
            double precision :: val

            ! T6: element (r=x, s | p, q=y)
            do s = x0, x1
                  do p = x0, x1
                        val = -rdm2(map(a), map(b), map(s), map(p))
                        if (val == zero) cycle
                        do y = y0, y1
                              do x = x0, x1
                                    call accum_A_ph(APlus, AMinus, pos(x, s), p, y, pos, val * V_axby(x, y))
                              end do
                        end do
                  end do
            end do

            ! T5: element (r, s=y | p=x, q)
            do q = x0, x1
                  do r = x0, x1
                        val = -rdm2(map(a), map(b), map(q), map(r))
                        if (val == zero) cycle
                        do y = y0, y1
                              do x = x0, x1
                                    call accum_A_ph(APlus, AMinus, pos(r, y), x, q, pos, val * V_axby(x, y))
                              end do
                        end do
                  end do
            end do
      end subroutine PhT56_gamma


end module phac0_thc
