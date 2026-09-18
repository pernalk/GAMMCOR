module ppac_types

      use iso_fortran_env
      use basis_definitions, only: TBasisAssignment

      use types
      use print_utils, only: print_memory


      implicit none
      integer, parameter :: I8 = int64
      
      integer, parameter :: VER_PH_TDA_HF = 0
      integer, parameter :: VER_PH_PPHH_TDA_HF = 1
      integer, parameter :: VER_PPHH_TDA_HF = 2
      integer, parameter :: VER_PH_RPA_HF = 3
      integer, parameter :: VER_PP_RPA_MULTI = 4
      integer, parameter :: VER_PH_PPHH_RPA_HF = 5
      integer, parameter :: VER_PPHH_RPA_HF = 6
      integer, parameter :: VER_PH_RPA_TDA_HF = 7
      integer, parameter :: VER_PP_RPA_TDA_HF = 8
      integer, parameter :: VER_PH_PPHH_RPA_TDA_HF = 9
      integer, parameter :: VER_PPHH_RPA_TDA_HF = 10
      integer, parameter :: VER_PP_RPATDA_MULTI = 11
      integer, parameter :: VER_HH_RPA_MULTI = 12

      integer, parameter :: BL_OOOO = 1
      integer, parameter :: BL_OOAO = 2
      integer, parameter :: BL_OOAA = 3
      integer, parameter :: BL_OOVA = 4
      integer, parameter :: BL_OOVV = 5
      integer, parameter :: BL_AOOO = 6
      integer, parameter :: BL_AOAO = 7
      integer, parameter :: BL_AOAA = 8
      integer, parameter :: BL_AOVA = 9
      integer, parameter :: BL_AOVV = 10
      integer, parameter :: BL_AAOO = 11
      integer, parameter :: BL_AAAO = 12
      integer, parameter :: BL_AAAA = 13
      integer, parameter :: BL_AAVA = 14
      integer, parameter :: BL_AAVV = 15
      integer, parameter :: BL_VAOO = 16
      integer, parameter :: BL_VAAO = 17
      integer, parameter :: BL_VAAA = 18
      integer, parameter :: BL_VAVA = 19
      integer, parameter :: BL_VAVV = 20
      integer, parameter :: BL_VVOO = 21
      integer, parameter :: BL_VVAO = 22
      integer, parameter :: BL_VVAA = 23
      integer, parameter :: BL_VVVA = 24
      integer, parameter :: BL_VVVV = 25

      integer, parameter :: GENVER_PP = 0
      integer, parameter :: GENVER_TDA_SING = 1
      integer, parameter :: GENVER_TDA_TRIP = 2
      integer, parameter :: GENVER_RPA_SING = 3
      integer, parameter :: GENVER_RPA_TRIP = 4

      type TAC0Block
            double precision, dimension(:, :), allocatable :: ASing, ATrip
            double precision, dimension(:, :), allocatable :: ATripA
            double precision, dimension(:), allocatable :: EigS, EigT, EigTA
            integer, dimension(:), allocatable :: vPlusS, vPlusT, vPlusTA
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

 type TphAC0Block
            double precision, dimension(:, :), allocatable :: APl, AMn
            double precision, dimension(:), allocatable :: EigPl, EigMn
            integer, dimension(:), allocatable :: vPlusPl, vPlusMn
            integer :: order
            integer, dimension(:, :), allocatable :: IndN
            integer, dimension(:, :), allocatable :: IndN1S, IndN1T
            integer, dimension(:, :), allocatable:: IndN2S, IndN2T
            integer, dimension(:, :), allocatable :: MiniMap
            integer :: NDim
            integer :: NDim1S, NDim2S, NDim1T, NDim2T
            character(:), allocatable :: name
            integer :: type
            type(TAC0MiniBlock), allocatable :: MiniBlocks(:)
      end type TphAC0Block

      type TAC0MiniBlock
            double precision, dimension(:,:), allocatable :: A
            double precision, dimension(:), allocatable :: Eig
            integer, dimension(:), allocatable :: ActList
            integer, dimension(:,:), allocatable :: IndN
            integer :: NDim
      end type TAC0MiniBlock      

       type TAC0vaBlock
            double precision, dimension(:,:), allocatable :: MiniAVS, MiniAVT, MiniAVTA
            double precision, dimension(:), allocatable :: EigS, EigT, EigTA
            integer, dimension(:), allocatable :: ActListS, ActListT, ActListTA
            integer, dimension(:,:), allocatable :: IndNS, IndNT, indNTA
            integer :: NDimS, NDimT, NDimTA
      end type TAC0vaBlock


      type TTHCData
            integer :: NTHC, NChol
            double precision, dimension(:,:), allocatable :: Xga, Zgk
            double precision, dimension(:,:), allocatable :: Xgp, XgpErf
            integer :: ExternalOrdering
            double precision, dimension(:), allocatable :: fij, fvw, ftu
            double precision, dimension(:,:,:), allocatable :: RFO, RFV
            double precision, dimension(:), allocatable :: eorbi, eorba
            integer :: NTHCErf, NCholErf
            double precision, dimension(:,:), allocatable :: XgaErf, ZgkErf
            double precision, dimension(:,:), allocatable :: TXgaErf, TXga
            double precision, dimension(:,:), allocatable :: HNO
            double precision, dimension(:,:), allocatable :: J_SR
            logical :: H0external = .false.
            
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
            integer, dimension(:), allocatable    :: IndX_s, IndX_t, IndX
            integer, dimension(:), allocatable    :: IndX_t_aa, IndX_t_bb
            integer, dimension(:,:), allocatable  :: IndN_s, IndN_t, IndN
            integer, dimension(:,:), allocatable  :: IndN_t_aa, IndN_t_bb
            integer, dimension(:), allocatable    :: IndAux
            integer, dimension(:, :), allocatable    :: IPair
            integer, dimension(:), allocatable :: int_alpha
            double precision :: alpha
            double precision, dimension(:, :), allocatable :: HNO0, HNO0_THC
            double precision, dimension(:, :), allocatable :: CAONO_a, CAONO_b
            double precision, dimension(:, :), allocatable :: CAONO, CMONO
            double precision, dimension(:, :), allocatable :: HNOA
            double precision, dimension(:), allocatable :: Eigs_s, Eigs_t, Eigs
            double precision, dimension(:), allocatable :: Eigs_t_aa, Eigs_t_bb
            integer, dimension(:), allocatable :: vplus_s, vplus_t
            integer, dimension(:), allocatable :: vplus_t_aa, vplus_t_bb
            double precision, dimension(:, :), allocatable :: Eigvec_s, Eigvec_t, Eigvec
            double precision, dimension(:, :), allocatable :: Eigvec_t_aa, Eigvec_t_bb
            double precision, dimension(:, :, :, :), allocatable :: rdm2_pp
            double precision, dimension(:, :, :, :), allocatable :: rdm2_pm
            double precision, dimension(:, :, :, :), allocatable :: rdm2
            double precision, dimension(:, :, :, :), allocatable :: ph_rdm2_aa
            double precision, dimension(:, :, :, :), allocatable :: ph_rdm2_ab
            double precision, dimension(:, :, :, :), allocatable :: rdm2_mm
            double precision, dimension(:, :, :, :), allocatable :: rdm2_mp
            double precision, dimension(:, :, :, :), allocatable :: rdm2_full
            double precision, dimension(:, :, :, :), allocatable :: pp2rdm0
            double precision, dimension(:, :), allocatable :: MxA, MxAt, MxS, MxSt
            double precision, dimension(:, :), allocatable :: ABplus, ABmin
            double precision, dimension(:, :), allocatable :: EigvecX
            double precision, dimension(:, :), allocatable :: EigvX, EigvY
            double precision, allocatable :: ABfull(:,:)
            double precision, allocatable :: Nmetric(:,:)
            type(TBasisAssignment) :: BasisAssign

            double precision, dimension(:, :), allocatable :: rdm1_full
            double precision, dimension(:, :), allocatable :: rdm1_p
            double precision, dimension(:, :), allocatable :: rdm1_m
            integer :: AAnegs, Apnegs, Amnegs


            double precision, dimension(:), allocatable :: TwoNO
            double precision, dimension(:), allocatable :: XOne
            double precision, allocatable :: Occ(:), Occ_rohf(:), cc(:)
            double precision, allocatable :: n_p(:), n_m(:), n(:)
            double precision :: W_aa_act0, W_ab_act0, W_aa0, W_ab0
            double precision :: ENuc
            double precision :: ECAS_read = 0.d0   ! CAS energy read from the external file (PySCF/ORCA)
            double precision :: ECAS_THC = 0.d0    ! CAS energy from THC-factorized integrals
            double precision :: ECAS_calc = 0.d0   ! CAS energy from RDMs + full TwoEl (check_energy_incore[_spinres])
            double precision :: ECAS_oneelectr = 0.d0
            double precision :: EROHF
            logical :: pherpa_print = .false.
            logical :: pperpa_print = .false.
            integer :: HType
            integer :: BatchDim
            integer :: general_version, version
            integer :: switch = 0
            integer :: ACType
            integer :: PYSCF = 0 
            integer :: ORCA = 0
            integer :: DALTON = 0
            integer :: NDim, NDim_s, NDim_t, NDim_t_aa, NDim_t_bb 
            integer :: NBasis, NA, NI, NIA, NV, NEL
            integer, dimension(:), allocatable    :: IndMod, map
            integer, dimension(:), allocatable    :: IndAuxFirst
            integer, dimension(2) :: nst
            integer :: true_NA, true_NI            
            integer :: NInte1, NInte2
            logical :: OnlyEnergy  = .false.
            logical :: pptriplet = .false.
            logical :: triplet = .false.
            integer :: iflmp2
            integer :: NCoreOrb
            double precision :: omega = 1.0
            double precision :: ThrPP
            logical :: spinsep = .true.
            double precision :: ThrSelAct, ThrQVirt, ThrQInact
            double precision :: E_ref_ducc
            integer :: omegaorders = 0
            
      end type TACppData

      type TInts
            integer :: ints1e_dim
            integer(8) :: ints2e_dim
            integer :: NInte2
            double precision, dimension(:, :), allocatable :: ints1e_aa, ints1e_bb
            double precision, dimension(:, :), allocatable :: Aints1e_aa, Aints1e_bb
            double precision, dimension(:), allocatable :: Aints2e_aa, Aints2e_ab
            double precision, dimension(:), allocatable :: ints2e_aa, ints2e_bb, ints2e_ab
            double precision, dimension(:), allocatable :: ints2e
      end type TInts


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

            ! rdm2_pp_12(i,j,k,l)  = rdm2_pp(i, k, j, l)
            double precision,allocatable :: rdm2_pp_12(:,:,:,:), rdm2_pm_12(:,:,:,:)
            ! rdm2_pp_13(i,j,k,l)  = rdm2_pp(i, l, k, j)                                                                                                                          
            double precision,allocatable :: rdm2_pp_13(:,:,:,:), rdm2_pm_13(:,:,:,:)
            ! like rdm2_pp_12 but only active part
            double precision,allocatable :: rdm2_pp_12_act(:,:,:,:), rdm2_pm_12_act(:,:,:,:)
            ! like rdm2_pp_13 but only active part 
            double precision,allocatable :: rdm2_pp_13_act(:,:,:,:), rdm2_pm_13_act(:,:,:,:)


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
#ifdef DEBUG
          call print_memory(name, this_mem)
#endif
    end subroutine CalcMem

    subroutine CalcMemI(A, name)
          integer, dimension(:,:), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2))
          this_mem = this_mem*togb
#ifdef DEBUG
          call print_memory(name, this_mem)
#endif
    end subroutine CalcMemI

    subroutine CalcMem3(A, name)
          double precision, dimension(:,:, :), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2)) * dble(size(A, dim=3))
          this_mem = this_mem*togb
#ifdef DEBUG          
          call print_memory(name, this_mem)
#endif
    end subroutine CalcMem3

    subroutine CalcMem4(A, name)
          double precision, dimension(:,:, :,:), intent(in) :: A
          double precision :: this_mem
          character(*), intent(in) :: name
          double precision :: togb

          togb = 8.d+0/(1024.0)**3

          this_mem = dble(size(A, dim=1)) * dble(size(A, dim=2)) * dble(size(A, dim=3)) * dble(size(A, dim=4))
          this_mem = this_mem*togb
#ifdef DEBUG          
          call print_memory(name, this_mem)
#endif
    end subroutine CalcMem4


      !    integer(8) function gmap(p, q, r, s)                                                               
      !       integer, external :: naddr3                                                                  
      !       integer, intent(in) :: p, q, r, s                                                            

      !       gmap = naddr3(p, q, r, s)                                                                    
      ! end function gmap           

    integer(I8) function gmap(p, q, r, s) result(idx)
        implicit none
        integer, intent(in) :: p, q, r, s
        integer(I8) :: ip, iq, ir, is
        integer(I8) :: ij, kl

        ip = int(p, I8)
        iq = int(q, I8)
        ir = int(r, I8)
        is = int(s, I8)

        ij = max(ip, iq) * (max(ip, iq) - 1_I8) / 2_I8 + min(ip, iq)
        kl = max(ir, is) * (max(ir, is) - 1_I8) / 2_I8 + min(ir, is)

        idx = max(ij, kl) * (max(ij, kl) - 1_I8) / 2_I8 + min(ij, kl)
    end function gmap    


    integer(I8) function gmap_4fold(p,q,r,s,nbasis) result(idx)
          implicit none
          integer, intent(in) :: p,q,r,s,nbasis
          integer(I8) :: ip,iq,ir,is, nb
          integer(I8) :: a1,b1,a2,b2, a,b, c, d

          ip = int(p,I8); iq = int(q,I8); ir = int(r,I8); is = int(s,I8)
          nb = int(nbasis,I8)

          a1 = ip + (iq-1_I8)*nb
          b1 = ir + (is-1_I8)*nb

          a = min(a1,b1)
          b = max(a1,b1)

          a2 = iq + (ip-1_I8)*nb  
          b2 = is + (ir-1_I8)*nb  

          c = min(a2,b2)
          d = max(a2,b2)
          
          if ( (c < a) .or. (c == a .and. d < b) ) then
                a = c
                b = d
          end if

          idx = b*(b-1_I8)/2_I8 + a
    end function gmap_4fold

       function erdm_ppx(rdm2, n, p, q, r, s, IAux, NI)
            double precision :: erdm_ppx
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux
            integer, intent(in) :: NI
            integer, intent(in) :: s, p, q, r

            if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
                  erdm_ppx = rdm2(p - NI, q - NI, r - NI, s - NI)
            else
                  erdm_ppx = 0.d+0
                  if (p==r.and.q==s)then
                        erdm_ppx = erdm_ppx + n(p) * n(q)
                  end if
                  if (p==s.and.q==r)then
                        erdm_ppx = erdm_ppx - n(p) * n(q)
                  end if
            end if

      end function erdm_ppx

      function erdm_pmx(rdm2, n_p, n_m, p, q, r, s, IAux, NI)
            double precision :: erdm_pmx
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            double precision, dimension(:), intent(in) :: n_p, n_m
            integer, dimension(:), intent(in) :: IAux
            integer, intent(in) :: NI
            integer, intent(in) :: s, p, q, r

            if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
                  erdm_pmx = rdm2(p - NI, q - NI, r - NI, s - NI)
            else
                  erdm_pmx = 0.d+0
                  if (p==r.and.q==s)then
                        erdm_pmx = erdm_pmx + n_p(p) * n_m(q)
                  end if
            end if

      end function erdm_pmx



end module ppac_types
