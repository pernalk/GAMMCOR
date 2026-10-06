module phac0_thc

      ! =============================================================================
      ! ph-AC0 with THC integrals. Mirrors the structure of ppac0_thc (ACPP_THC),
      ! but only the AC0 part and for the ph (particle-hole) ERPA problem.
      !
      ! THEORY (see .skills/phAC0-theory.tex and .skills/phtheory_a0_blocks.md):
      !   The half-sized problem  A_+ A_- Ytilde = omega^2 Ytilde  is solved per
      !   A0 block, with A_+/- = N^{-1/2}(A+-B)N^{-1/2}. The ph metric is
      !   N_pq = n_p - n_q (NOT the pp metric 1 - n_p - n_q). W^(1) then couples
      !   the zeroth-order vectors through the A1 blocks (eq. with factor 4 in
      !   the tex). A^(1) is symmetric for optimized orbitals (Brillouin), so of
      !   each pair of pair-swapped A1 blocks only one is built.
      !
      ! BLOCK STRUCTURE (see .skills/thc4_loops2.md for the integral classes):
      !   A0 is block diagonal with 4 blocks: vo (VIVI), aa (AAAA), va (VAVA),
      !   oa (AIAI). Here, as everywhere in the Fortran, the letter o means
      !   INACTIVE. There is no all-inactive and no all-virtual ph block.
      !   A1 has 9 blocks: the 3 diagonal-position ones (vovo, vava, oaoa) and
      !   6 off-diagonal couplings (one representative per transposed pair).
      !   There is no aaaa A1 block: the aa-aa part of A sits entirely in A0.
      !
      ! STORAGE (ph analogue of .skills/ppcode_a0_a1_blocks.md):
      !   One matrix pair (APlus/AMinus) per block - no ASing/ATrip/ATripA, no
      !   IndNT/posT/indt: everything that exists twice in pp because of the
      !   triplet has no ph counterpart.
      !     vo      : diagonal, APlus/AMinus(NDim,1) hold one column
      !     aa      : full square APlus/AMinus(NDim,NDim)
      !     va / oa : block-diagonal in the non-active index -> MiniBlocks,
      !               one NA x NA mini matrix per virtual (va) / inactive (oa),
      !               active-local indexing, ThrQVirt/ThrQInact null filters
      !     A1      : rectangular couplings; small ones allocated, the ones
      !               involving the vo pair space (and vava) contracted OTF
      !   The global pair -> position map is the single array pos(p,q)
      !   (ph analogue of posS; there is no posT).
      !
      ! PAIR CONVENTIONS (same first-slot rule as pp: V > A > I, see
      ! InitA0Blocks_ph): vo stores (virt, inact), aa stores (act_i, act_j)
      ! with i > j, va stores (virt, act), oa stores (act, inact).
      ! Null-space filters are the ph ones: aa drops nearly degenerate active
      ! pairs (ThrSelAct, metric n_p - n_q ~ 0); va mini-blocks drop
      ! quasi-virtual actives (ThrQVirt); oa mini-blocks drop quasi-inactive
      ! actives (ThrQInact). These are the same filters as in
      ! acph_spinres_driver (phac_spinres.f90), applied per block.
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
      use thc_energy
      use tran
      use omp_lib


      implicit none

      ! ph block ids. Prefixed ph* so they cannot be confused with the pp ids
      ! A0oo..A1vaao2 of ppac0_thc (both sets are private to their module).
      ! A1 ids 5..9 are the blocks contracted on the fly (OTF), 10..13 the
      ! allocated ones - same split idea as pp, where the vv* A1 blocks are OTF.
      integer, parameter, private :: phA0vo=1, phA0aa=2, phA0va=3, phA0oa=4
      integer, parameter, private :: phA1vovo=5, phA1vavo=6, phA1vooa=7
      integer, parameter, private :: phA1voaa=8, phA1vava=9
      integer, parameter, private :: phA1vaoa=10, phA1vaaa=11, phA1aaoa=12
      integer, parameter, private :: phA1oaoa=13
      integer, parameter, private :: NBLOCKS_PH = 13

      integer, parameter, private :: tdebug = 0

      character(len=8), dimension(NBLOCKS_PH), parameter :: BlockName_ph = (/ &
            'A0vo    ', 'A0aa    ', 'A0va    ', 'A0oa    ', &
            'A1vovo  ', 'A1vavo  ', 'A1vooa  ', 'A1voaa  ', 'A1vava  ', &
            'A1vaoa  ', 'A1vaaa  ', 'A1aaoa  ', 'A1oaoa  ' /)

contains


      subroutine ACPH0_THC(THCData, AuxData, CAONO_IN, Flags)
            use THC_Gammcor
            use ppac0_thc, only: update_HNO_AC0

            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), intent(in) :: CAONO_IN
            type(FlagsData), intent(in) :: Flags
            type(TTHCData), intent(inout) :: THCData

            type(TAC0Block_ph), dimension(NBLOCKS_PH) :: ACBlocks
            integer, dimension(:,:), allocatable :: pos
            double precision, dimension(:,:), allocatable :: AuxII, AuxAA, BuxII, BuxAA
            double precision, dimension(:,:), allocatable :: Aux3X
            double precision, dimension(:,:), allocatable :: HNO


            if (AuxData%PYSCF==1)then
                  call print_info('External ordering', 'PySCF')
            end if

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

              call THC_int_loop_ph(THCData, ACBlocks, AuxII, AuxAA, BuxII, BuxAA, &
                    Aux3X, Flags, AuxData, Occ, map, 0.5d0*AuxData%rdm2_full, pos, One)

              ! HNO: generalized Fock. Active block = HNO0 + 2*BuxII - AuxII
              ! (= h^core), inactive/virtual diagonals from fij/fvw.
              allocate(HNO(NBasis, NBasis))
              call update_HNO_AC0(AuxData, HNO, THCData%fij, THCData%fvw, AuxII, BuxII)

              call PHERPA0_THC_VO(AuxData, ACBlocks(phA0vo), THCData)
              call PHERPA0_THC_AA(AuxData, ACBlocks(phA0aa), HNO, Aux3X)
              call PHERPA0_THC_VA_Eigs(AuxData, ACBlocks(phA0va), THCData, HNO, Aux3X)
              call PHERPA0_THC_OA_Eigs(AuxData, ACBlocks(phA0oa), THCData, HNO, Aux3X, AuxAA, BuxAA)

              ! -------- next steps (not yet implemented) --------
              !A1 blocks (phA1vaoa..phA1oaoa) + OTF ones in the energy loop
              !call THC_energy_loop_ph

            end associate

            call msg('ACPH0_THC: A0 blocks assembled and diagonalized; A1/energy not yet implemented')
            stop

      end subroutine ACPH0_THC


      ! ph analogue of init_AC0Blocks (ppac0_thc). Builds the 4 A0 pair lists
      ! and allocates the small A1 coupling blocks. A1 blocks that involve the
      ! vo pair space, and vava, are too large to store - contracted OTF in the
      ! energy loop, so only their parents' dimensions are recorded here.
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

            ! These blocks are too large and are calculated OTF
            !call InitA1Blocks_ph(AC0Block, AuxData, phA1vovo, phA0vo, phA0vo)
            !call InitA1Blocks_ph(AC0Block, AuxData, phA1vavo, phA0va, phA0vo)
            !call InitA1Blocks_ph(AC0Block, AuxData, phA1vooa, phA0vo, phA0oa)
            !call InitA1Blocks_ph(AC0Block, AuxData, phA1voaa, phA0vo, phA0aa)
            !call InitA1Blocks_ph(AC0Block, AuxData, phA1vava, phA0va, phA0va)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1vaoa, phA0va, phA0oa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1vaaa, phA0va, phA0aa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1aaoa, phA0aa, phA0oa)
            call InitA1Blocks_ph(AC0Block, AuxData, phA1oaoa, phA0oa, phA0oa)

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


      ! ph analogue of InitA0Blocks (ppac0_thc). Builds the pair list of one
      ! A0 block, fills the global pos(p,q) map, and allocates the storage.
      !
      ! Pair lists (first slot = higher space, V > A > I):
      !   vo : all (v, i) pairs. The block is diagonal (canonical orbitals,
      !        diagonal fockian), so APlus/AMinus get one column only.
      !   aa : active pairs (i, j), i > j, skipping nearly degenerate pairs
      !        (ThrSelAct) - those have ph metric n_i - n_j ~ 0 (null space).
      !   va : all (v, a) pairs in the main list; per-virtual MiniBlocks with
      !        ActList of actives passing |n_v + n_a| > ThrQVirt
      !        (n_a ~ 0 gives metric ~ 0).
      !   oa : all (a, i) pairs in the main list; per-inactive MiniBlocks with
      !        ActList of actives passing |n_i + n_a - 2| > ThrQInact
      !        (n_a ~ 1 gives metric ~ 0).
      !
      ! TODO: IFlCore handling (frozen inactive pairs) as in
      ! acph_spinres_driver is not applied here yet.
      subroutine InitA0Blocks_ph(A0Block, AuxData, pos, Ip, Iq)
            type(TAC0Block_ph), intent(inout) :: A0Block
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:,:), intent(inout) :: pos
            integer, intent(in) :: Ip, Iq
            integer, parameter :: SP_OCC = 0, SP_ACT = 1, SP_VIRT = 2

            integer :: ind, inda
            integer :: i, j, k, Nbl, Nbl0, Nbl1

            associate(NA=>AuxData%NA, NV=>AuxData%NV, NI=>AuxData%NI, &
                  NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, Occ=>AuxData%Occ)

              if (Ip == SP_VIRT .and. Iq == SP_OCC)then
                    ! ------------- vo block: diagonal -------------
                    A0Block%type = phA0vo
                    A0Block%NDim = NV * NI
                    allocate(A0Block%IndN(2, A0Block%NDim))

                    ind = 0
                    do i = NIA+1, NBasis
                          do j = 1, NI
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
                    Nbl = NI
                    Nbl0 = 1
                    Nbl1 = NI
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
                    do k = 1, Nbl
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


      ! Which integral classes run, ph-AC0 version of parse_config (ppac0_thc).
      ! For the A0 blocks only three classes are needed (see the module
      ! docstring and .skills/phtheory_a0_blocks.md):
      !   (IA|IA) -> AuxII active block   sum_{t in I} n_t (qt|ts)
      !   (AA|II) -> BuxII active block   sum_{t in I} n_t (qs|tt)
      !   (AA|AA) -> AuxAA/BuxAA, Aux3X (P45_ph) and the aaaa gamma kernels
      ! Everything else (A1 classes) stays off until the A1 terms are coded.
      ! TabChol = [1,1,1,0,0,0]: no R vector with a virtual index is needed yet.
      subroutine parse_config_ph(loops, TabChol)
            type(TLoop), dimension(6), intent(out) :: loops
            integer, dimension(6), intent(out) :: TabChol
            integer :: i

            do i = 1, 6
                  loops(i)%run_main = .false.
                  loops(i)%run_inner = .false.
            end do

            loops(2)%run_main = .true.        ! (IA|..
            loops(2)%run_inner(2) = .true.    ! (IA|IA)

            loops(3)%run_main = .true.        ! (AA|..
            loops(3)%run_inner(1) = .true.    ! (AA|II)
            loops(3)%run_inner(3) = .true.    ! (AA|AA)

            TabChol = (/ 1, 1, 1, 0, 0, 0 /)

      end subroutine parse_config_ph


      ! =============================================================================
      ! ph analogue of THC_int_loop (ppac0_thc). Same skeleton: R?? Cholesky-like
      ! vectors from the THC factors, batch loop, caseloop over the six pair
      ! classes, V_axby = (a x| b y) built per fixed (a,b) by real_aTb_x and
      ! consumed immediately by the kernels - nothing four-index is ever stored.
      !
      ! CURRENT SCOPE: only the A0-necessary terms (ph-AC0, no A1 yet):
      !   (IA|IA)  p2_IAIA_qtst_tqts  -> AuxII active block  sum_{t in I} n_t (qt|ts)
      !   (AA|II)  p2_AAII_qstt_qstt  -> BuxII active block  sum_{t in I} n_t (qs|tt)
      !   (AA|AA)  p2_AAAA_qstt_qstt  -> BuxAA               sum_{t in A} n_t (qs|tt)
      !            p2_AAAA_qtst_qtst  -> AuxAA               sum_{t in A} n_t (qt|ts)
      !            P45_ph             -> Aux3X, the genuine WMAT matrix
      !            PhT12rp_gamma, PhT12qs_gamma, PhT56_gamma
      !                               -> ACB(phA0aa)%APlus/AMinus, the genuine
      !                                  two-electron T terms of the aaaa block
      ! The p2_* kernels are the generated pp ones (thc_auto_p2), called with
      ! the same signatures as in pp THC_int_loop. P45_ph and the Ph*_gamma
      ! kernels are the ph counterparts of pp's P45 / P3c*_gamma.
      !
      ! RDM: ONE spin-summed active 2-RDM array rdm2 (the ph analogue of
      ! pp's rdm2_nsum; today built as rdm2_pp + rdm2_pm in the native storage
      ! layout). No kernel ever sees spin-separated 2-RDMs.
      !
      ! On exit ACB(phA0aa)%APlus/AMinus hold the UNSCALED two-electron genuine
      ! part only (no T0/h^core, no delta-WMAT terms, no metric division) -
      ! exactly the role pp's int loop plays for ASing; the rest is added by
      ! the assembly routines (next step).
      ! =============================================================================
      subroutine THC_int_loop_ph(THCData, ACB, AuxII, AuxAA, BuxII, BuxAA, &
            Aux3X, Flags, AuxData, Occ, map, rdm2, pos, AcAlpha)
            use THC_Gammcor
            use thc_auto_p2

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
            integer, dimension(:,:), intent(in) :: pos
            double precision, intent(in) :: AcAlpha

            double precision, dimension(:), allocatable :: V_axby
            double precision, dimension(:,:,:), allocatable :: RII, RIA, RAA
            integer, dimension(6) :: TabChol
            type(TLoop), dimension(6) :: loops
            integer, parameter :: loopII=1, loopIA=2, loopAA=3, loopIV=4, loopVA=5, loopVV=6
            integer :: NI, NA, NV, NIA, NBasis
            integer :: a, b, ii
            integer :: batch, BatchDim, NBatchI, NBatchA, NBatchV
            type (tclock) :: timer0

            associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga, &
                  NChol=>THCData%NChol, NTHC=>THCData%NTHC)

              NA = AuxData%NA
              NI = AuxData%NI
              NV = AuxData%NV
              NIA = NI + NA
              NBasis = AuxData%NBasis

              call parse_config_ph(loops, TabChol)

              allocate(V_axby(max(NA, NI, NV)**2))

              ! Batching: kept for structural parity with pp. With
              ! TabChol(4:6) = 0 nothing virtual is ever built, so a single
              ! batch does all the work. When the A1 classes arrive, bring in
              ! the RIV/RVA/RVV window refills from pp THC_int_loop - and mind
              ! the BatchDim defects (.skills/ppcode_known_defects.md, A-D).
              BatchDim = 10000
              NBatchI = NI / BatchDim
              NBatchA = NA / BatchDim
              NBatchV = NV / BatchDim
              if (modulo(NI, BatchDim)>0) NBatchI = NBatchI + 1
              if (modulo(NA, BatchDim)>0) NBatchA = NBatchA + 1
              if (modulo(NV, BatchDim)>0) NBatchV = NBatchV + 1

              call print_section('ph THC integral transformation')

              if (TabChol(1) == 1) then
                    allocate(RII(NChol, NI, NI))
                    call CalcMem3(RII, 'RII')
                    call thc_gammcor_Rkab_2(RII, Xga(:,1:NI), Xga(:,1:NI), Zgk, NI, NI, NChol, NTHC)
              end if
              if (TabChol(2) == 1) then
                    allocate(RIA(NChol, NA, NI))
                    call CalcMem3(RIA, 'RIA')
                    call thc_gammcor_Rkab_2(RIA, Xga(:,NI+1:NIA), Xga(:,1:NI), Zgk, NA, NI, NChol, NTHC)
              end if
              if (TabChol(3) == 1) then
                    allocate(RAA(NChol, NA, NA))
                    call CalcMem3(RAA, 'RAA')
                    call thc_gammcor_Rkab_2(RAA, Xga(:,NI+1:NIA), Xga(:,NI+1:NIA), Zgk, NA, NA, NChol, NTHC)
              end if

              call clock_start(timer0)

              batchloop: do batch = 1, max(NBatchI, NBatchA, NBatchV, 1)

                    ! (RIV / RVA / RVV window refills go here - A1 step)

                    caseloop: do ii = 1, 6
                          if (loops(ii)%run_main)then
                                select case(ii)
                                case(loopII)
                                      ! (II|II) - A1 step (phA1oaoa etc. need it)
                                case(loopIA)
                                      do a = 1, NI
                                            ! run_inner(1): (IA|II) - A1 step
                                            if (loops(ii)%run_inner(2)) then ! (IA|IA)
                                                  do b = 1, NI
                                                        call real_aTb_x(V_axby, NA, RIA(:,:,a), NChol, RIA(:,:,b), NChol, NA, NA, NChol, AcAlpha)
                                                        call p2_IAIA_qtst_tqts(AuxII, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                  end do
                                            end if
                                      end do
                                case(loopAA)
                                      do a = NI+1, NIA
                                            if (loops(ii)%run_inner(1)) then ! (AA|II)
                                                  do b = 1, NI
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RII(:,:,b), NChol, NA, NI, NChol, AcAlpha)
                                                        call p2_AAII_qstt_qstt(BuxII, a, b, NI+1, NIA, 1, NI, NI, NIA, NBasis, Occ, V_axby)
                                                  end do
                                            end if
                                            ! run_inner(2): (AA|IA) - A1 step
                                            if (loops(ii)%run_inner(3)) then ! (AA|AA)
                                                  do b = NI+1, NIA
                                                        call real_aTb_x(V_axby, NA, RAA(:,:,a-NI), NChol, RAA(:,:,b-NI), NChol, NA, NA, NChol, AcAlpha)
                                                        call p2_AAAA_qstt_qstt(BuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call p2_AAAA_qtst_qtst(AuxAA, a, b, NI+1, NIA, NI+1, NIA, NI, NIA, NBasis, Occ, V_axby)
                                                        call P45_ph(Aux3X, a, b, map, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        call PhT12rp_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        call PhT12qs_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                        call PhT56_gamma(ACB(phA0aa)%APlus, ACB(phA0aa)%AMinus, a, b, map, pos, rdm2, NI+1, NIA, NI+1, NIA, V_axby, NI, NA)
                                                  end do
                                            end if
                                      end do
                                case(loopIV)
                                      ! (IV|..) classes - A1 step (vo pair space)
                                case(loopVA)
                                      ! (VA|..) classes - A1 step
                                case(loopVV)
                                      ! (VV|..) classes - A1 step
                                end select
                          end if
                    end do caseloop

              end do batchloop

              call tmsg('ph int loop', timer0, tdebug)

            end associate

      end subroutine THC_int_loop_ph


      ! =============================================================================
      ! A0 assembly + diagonalization, ph versions of PPERPA0_THC_OOVV /
      ! PPERPA0_THC_AA / PPERPA0_THC_VA_Eigs / PPERPA0_THC_OA_Eigs.
      !
      ! Per-block A^(0) (see .skills/phtheory_a0_blocks.md; delta terms exactly
      ! as T0/TW1/TW2 of pherpa_incore_spinres at alpha = 0):
      !   vo : d_pr d_sq (f_vv - f_ii), diagonal; f from fij/fvw
      !   aa : [2-el genuine, already in APlus/AMinus from the int loop]
      !        + d_pr (n_p - n_s) HNO(q,s) + d_sq (n_q - n_r) HNO(p,r)
      !        - 1/2 d_sq Aux3X(p,r) - 1/2 d_pr Aux3X(q,s)
      !   va : per virtual v (p = r = v):  A(s,q) =
      !        - n_s HNO(q,s) - 1/2 Aux3X(q,s) + d_sq n_q f_vv
      !   oa : per inactive i (s = q = i): A(r,p) =
      !        (1 - n_r) HNO(p,r) + 2 BuxAA(p,r) - AuxAA(p,r)
      !        - 1/2 Aux3X(p,r) + d_pr (n_p - 1) f_ii
      ! For va and oa the B part is zero at alpha = 0 (all its integrals are
      ! first order), so A_+ and A_- differ only by the metric scaling.
      !
      ! Scaling (same c convention as pherpa_incore_spinres):
      !   APlus  /= (c_p + c_q)(c_r + c_s),   AMinus /= (c_p - c_q)(c_r - c_s)
      ! then the half-sized problem A_+ A_- Ytilde = omega^2 Ytilde is solved by
      ! pherpa_symm_block_ph (port of pherpa_symm_spinres). On exit APlus holds
      ! the normalized Ytilde vectors in columns, AMinus the symmetrized scaled
      ! A_-, Eig the omega (negative entry = -sqrt|omega^2| flags a negative
      ! omega^2 mode). Each routine prints its eigenvalues.
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

            call print_block_eigs_ph('A0vo', ACB_vo%Eig, ACB_vo%NDim)

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

              call pherpa_symm_block_ph(ACB_aa%APlus, ACB_aa%AMinus, ACB_aa%Eig, nd, nneg)
              call print_block_eigs_ph('A0aa', ACB_aa%Eig, nd)
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

                    call pherpa_symm_block_ph(workp, workm, eigw, nd, nneg)

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

              call print_block_eigs_ph('A0va', allgather, ngather)

            end associate

      end subroutine PHERPA0_THC_VA_Eigs


      ! oa block: one NA x NA eigenproblem per inactive i (s = q = i).
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

              do k = 1, NI
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

                    call pherpa_symm_block_ph(workp, workm, eigw, nd, nneg)

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

              call print_block_eigs_ph('A0oa', allgather, ngather)

            end associate

      end subroutine PHERPA0_THC_OA_Eigs


      ! Port of pherpa_symm_spinres (phac_spinres.f90) for one block.
      ! Solves A_+^(1/2) A_- A_+^(1/2) Z = omega^2 Z; Ytilde = A_+^(1/2) Z,
      ! normalized to 2/omega Ytilde^T A_- Ytilde = 1.
      ! In/out: APlus(nd,nd) in: scaled A_+, out: Ytilde vectors in columns;
      ! AMinus(nd,nd) in: scaled A_-, out: its symmetrized version;
      ! Eig(nd) out: omega, with -sqrt|omega^2| flagging negative modes (nneg).
      subroutine pherpa_symm_block_ph(APlus, AMinus, Eig, nd, nneg)
            double precision, dimension(:,:), intent(inout) :: APlus, AMinus
            double precision, dimension(:), intent(out) :: Eig
            integer, intent(in) :: nd
            integer, intent(out) :: nneg

            double precision, dimension(:,:), allocatable :: work_a, work_b
            double precision, dimension(:), allocatable :: lambda, vec_tmp
            double precision :: omega_sq, dot_prod, norm_fact
            integer :: j, mode
            double precision, parameter :: small = 1.d-10

            nneg = 0
            if (nd == 0) return

            allocate(work_a(nd, nd))
            allocate(work_b(nd, nd))
            allocate(lambda(nd))
            allocate(vec_tmp(nd))

            ! the raw blocks are symmetric only up to the CASSCF gradient
            ! (Brillouin, see the Important note in phAC0-theory.tex)
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


      ! ph counterpart of pp's P45: the genuine (all-active) WMAT matrix,
      ! a literal port of the wmat loop of pherpa_incore_spinres restricted to
      ! all-active summation indices:
      !   wmat(p,r) = sum_{tuv in act} [ rdm2(r,v,t,u) (uv|tp)
      !                                + rdm2(r,v,u,t) (up|tv) ]
      ! Class (AA|AA), V(x,y) = (a x| b y):
      !   pattern 1: u=a, v=x, t=b, p=y -> Aux3X(y,r) += sum_x rdm2(r,x,b,a) V(x,y)
      !   pattern 2: u=a, p=x, t=b, v=y -> Aux3X(x,r) += sum_y rdm2(r,y,a,b) V(x,y)
      ! Consumed by the assembly as TW1 = -1/2 d_sq Aux3X(p,r),
      ! TW2 = -1/2 d_pr Aux3X(q,s)  (tex: WMATpr and WMATqs).
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
