module ppac_subs

      use iso_fortran_env

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals
      use display
      use types
      use lin
      use clock

      use math_constants
      use quadratures
      use ppac_types
      use ppac0
      use pp_utils
      use interface_pp

      implicit none


      integer date_time(8)
      character*10 b(3)

contains

      subroutine ppAC_full(ETot, ECorr_AC, XOne, TWONO, CAONO,  Flags, AuxData)

            double precision, intent(inout) :: ETot, ECorr_AC
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),   intent(in) :: TwoNO
            double precision, dimension(:,:), intent(in) :: CAONO

            type(FlagsData), intent(in) :: Flags
            type(TACppData) :: AuxData

            double precision, dimension(25) :: E_corr_contr_s, E_corr_contr_t, E_corr_contr_ext
            double precision, dimension(25) :: E_corr_contr_big
            double precision :: ACAlpha
            double precision :: Acpp_alpha_s, Acpp_alpha_t
            double precision :: Acpp_alpha_s_wo, Acpp_alpha_t_wo
            double precision :: Acpp_alpha_ext, Acpp_alpha_ext_wo
            type (tclock) :: start, timer, timer0
            double precision :: acalpha0, w_lin, ECorr_AC0, ECorr, ECorr_wo, ECorr_2_n
            integer :: i
            double precision, allocatable :: eorb(:,:)
            type(TDA_pphh) :: Tpphh
            double precision, dimension(25) :: val


            integer :: NGrid
            double precision, dimension(:), allocatable :: XGrid, WGrid
            integer, parameter :: GridD = 100


            if(AuxData%HType == H_DYALL)then
                  print*, 'Dyall'
            else if (AuxData%HType == H_GPF)then
                  print*, 'GPF'
            end if

            ECorr = zero
            ECorr_wo = zero
            ECorr_AC0 = zero
            ECorr_2_n = zero
            ECorr_AC = zero


            associate(Occ=>AuxData%Occ, NI=> AuxData%NI, NA=> AuxData%NA, NV=> AuxData%NV, &
                  NBasis =>AuxData%NBasis, NInte1=> AuxData%NInte1, NInte2=> AuxData%NInte2,&
                  IndAux=>AuxData%IndAux, IndMod => AuxData%IndMod, &
                  ENuc=> AuxData%ENuc, version=>AuxData%version, AC_TYPE=>AuxData%ACType)


              allocate(eorb(NBasis, Nbasis))

              ACAlpha = zero
              E_corr_contr_big = zero

              E_corr_contr_s = zero
              E_corr_contr_t = zero
              E_corr_contr_ext = zero

              print*, 'poczatkowe acalpha', acalpha
              call ppAC_int(AuxData%HType, E_corr_contr_s, E_corr_contr_t, E_corr_contr_ext, &
                    Acpp_alpha_s, Acpp_alpha_t,Acpp_alpha_s_wo, Acpp_alpha_t_wo, &
                    Acpp_alpha_ext, Acpp_alpha_ext_wo, &
                    ETot, ENuc, Occ, XOne, TwoNO, AuxData, IndAux, IndMod, &
                    NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, ACAlpha, version, Flags, AC_TYPE, eorb)

              write(*, '(A10, 4F15.8)') 'wklady', acpp_alpha_s, acpp_alpha_t, Acpp_alpha_s + Acpp_alpha_t ,Acpp_alpha_ext
              write(*, '(A20, F15.8)') 'W(alpha=zero)', Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext
              print*, 'Acpp_alpha_s + Acpp_alpha_t', Acpp_alpha_s , Acpp_alpha_t
              if (abs(Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext).gt.1.d-6)then
                    print*, 'W(alpha = 0) != 0', Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext
                    stop
              end if

              NGrid = 6
              call imsg("NGrid for ACPP", ngrid)
              allocate(XGrid(GridD))
              allocate(WGrid(GridD))

              Call GauLeg(Zero,One,XGrid,WGrid,NGrid)

              call clock_start(timer)
              acalpha0 = 1.d-6
              call ppAC_int(AuxData%HType, E_corr_contr_s, E_corr_contr_t, &
                    E_corr_contr_ext, Acpp_alpha_s, Acpp_alpha_t, Acpp_alpha_s_wo, Acpp_alpha_t_wo, &
                    Acpp_alpha_ext, Acpp_alpha_ext_wo, &
                    ETot, ENuc, Occ, XOne, TwoNO, AuxData, IndAux, IndMod, &
                    NBasis, NA, NI, NV, NInte1, NInte2, &
                    Tpphh, ACAlpha0, version, Flags, AC_TYPE, eorb)
              print*, 0, 'Czas na iteracje: ', clock_readwall(timer)

              w_lin = (Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext)/acalpha0
              print*, Acpp_alpha_s , Acpp_alpha_t , Acpp_alpha_ext, acalpha0
              print*, 'wlin', w_lin
              ECorr_AC0 = frac12 *w_lin


              do i = 1, NGrid
                    call clock_start(timer)
                    ACAlpha = XGrid(i)
                    print*, ''
                    print*, 'iteration',i,  acalpha, wgrid(i)
                    print*, ''
                    E_corr_contr_s = zero
                    E_corr_contr_t = zero
                    E_corr_contr_ext = zero

                    !             version = VER_PP_RPA_TDA_HF
                    call ppAC_int(AuxData%HType, E_corr_contr_s, E_corr_contr_t, &
                          E_corr_contr_ext, Acpp_alpha_s, Acpp_alpha_t, Acpp_alpha_s_wo, Acpp_alpha_t_wo, &
                          Acpp_alpha_ext, Acpp_alpha_ext_wo, &
                          ETot, ENuc,Occ, XOne, TwoNO, AuxData, IndAux, IndMod, &    
                          NBasis, NA, NI, NV, NInte1, NInte2, &
                          Tpphh, ACAlpha, version, Flags, AC_TYPE, eorb)


                    write(*, '(A10, 4F15.8)') 'wklady', acpp_alpha_s, acpp_alpha_t, Acpp_alpha_s + Acpp_alpha_t , Acpp_alpha_ext
                    if (i==1)then
                          write(*, '(A10, F15.8)') 'W(alpha1)', Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext
                          write(*,'(A20, F20.15)')  'zzz', ECorr_AC0
                    else
                          write(*, '(A10, F15.8)') 'W(alpha)', Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext
                    end if

                    write(*,'(A20, 2F20.15)')  'RDSC W(a)', ACAlpha, Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext

                    write(*,'(A20, F20.15)')  'RDSC W1(alpha)', w_lin * ACAlpha


                    ECorr = ECorr + (Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext) * WGrid(i)
                    val = E_corr_contr_s + E_corr_contr_t+E_corr_contr_ext
                    !write(*, '(A10, 5F20.15)')'poln', E_corr_contr_s(BL_VVVV) , E_corr_contr_t(BL_VVVV),E_corr_contr_ext(BL_VVVV),  WGrid(i), (E_corr_contr_s(BL_VVVV) + E_corr_contr_t(BL_VVVV)+E_corr_contr_ext(BL_VVVV)) * WGrid(i)
                    E_corr_contr_big = E_corr_contr_big + (E_corr_contr_s + E_corr_contr_t+E_corr_contr_ext) * WGrid(i)

                    E_corr_contr_t = zero 
                    !call print_contributions(E_corr_contr_s+ E_corr_contr_t+E_corr_contr_ext, E_corr_contr_t)
                    call print_contributions_alpha(val, ACAlpha)
                    print*, E_corr_contr_big(BL_VVVV), 'E_corr_contr_big(vvvv)'

                    ECorr_wo = ECorr_wo + (Acpp_alpha_s_wo + Acpp_alpha_t_wo + Acpp_alpha_ext_wo) * WGrid(i)
                    ECorr_2_n = ECorr_2_n + ((Acpp_alpha_s + Acpp_alpha_t + Acpp_alpha_ext) - w_lin * ACAlpha) * WGrid(i)
                    print*, 'teraz ECorr', ECorr
                    print*, i, 'Czas na iteracje: ', clock_readwall(timer)
                    print*, ''
              end do

              print*, 'ost ECorr', ECorr              
              print*, 'ost ECorr_AC0', ECorr_AC0
              
              print*, ETot+ECorr+ENuc
              print*, 'RDSC CAS', ETot+ENuc
              print*, ''
              print*, 'RDSC ACPP ECORR', ETot+ECorr+ENuc
              print*, ''
              print*, 'RDSC ACPP0 ECORR', ETot+ECorr_AC0+ENuc
              print*, ''
              print*, 0, 'Czas na old version: '//str(clock_readwall(timer0),d=2)
              ECorr_AC = ECorr
              E_corr_contr_s = zero
              call print_contributions(E_corr_contr_big, E_corr_contr_s)


            end associate

      end subroutine ppAC_full



      subroutine pp_erpa(XOne, TWONO, CAONO, Flags, AuxData)

            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),   intent(in) :: TwoNO
            double precision, dimension(:,:), intent(in) :: CAONO

            type(FlagsData), intent(in) :: Flags
            type(TACppData) :: AuxData

            double precision, dimension(:,:), allocatable :: MxA, MxS
            double precision, dimension(:), allocatable :: w, wi
            double precision, dimension(:,:), allocatable :: ltz, rtz
            integer, dimension(:), allocatable :: v_plus
            double precision :: ETot
            integer :: spin_symm
            double precision, allocatable :: eorb(:,:)
            type(TDA_pphh) :: Tpphh
            type (tclock)  :: timer
            integer :: i, k


            print*, AuxData%HType
            if(AuxData%HType==H_DYALL)then
                  print*, 'Dyall'
            else if (AuxData%HType==H_GPF)then
                  print*, 'GPF'
            end if

            if (AuxData%triplet)then
                  spin_symm  = 1
            else
                  spin_symm = 0
            end if



            allocate(MxA(AuxData%NDim, AuxData%NDim))
            allocate(MxS(AuxData%NDim, AuxData%NDim))
            allocate(w(AuxData%NDim))
            allocate(wi(AuxData%NDim))
            allocate(ltz(AuxData%NDim, AuxData%NDim))
            allocate(rtz(AuxData%NDim, AuxData%NDim))
            allocate(v_plus(AuxData%Ndim))


            associate(Occ=>AuxData%Occ, NI=> AuxData%NI, NA=> AuxData%NA, NV=> AuxData%NV, &
                  NBasis =>AuxData%NBasis, NInte1=> AuxData%NInte1, NInte2=> AuxData%NInte2,&
                  IndAux=>AuxData%IndAux, IndMod => AuxData%IndMod, &
                  ENuc=> AuxData%ENuc, version=>AuxData%version, AC_TYPE=>AuxData%ACType)

              allocate(eorb(NBasis, Nbasis))            
              call pp_rowe(AuxData, AuxData%HType, ETot, ENuc, MxA, MxS, Occ, &
                    XOne, TwoNO, AuxData%IndN_s, AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, &
                    NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, zero, spin_symm, version, Flags, eorb)


              call clock_start(timer)

              print*, 'version', version

              print*, "diagonalizacja macierzy niesymetrycznej"
              call nonsymmetric_eigenproblem_pp(w, wi, rtz, MxA, MxS, v_plus,  AuxData%IndN, IndAux, AC_TYPE)


              print*, 'Czas diagonalizacji: ', clock_readwall(timer)
              call clock_start(timer)    
              print*, 'Eigenvalues of ppRPA problem'
              k = 0
              do i = 1, AuxData%NDim_s
                    print*, w(i), wi(i)
                    k = k + 1
              end do
              print*, 'oto bedzie norma', ETot
              call compute_norm(AuxData, w, wi, rtz, v_plus, AuxData%NDim_s, AuxData%IndN_s, AuxData%IndX_s, Occ, 0, Tpphh, spin_symm, ENuc)


              print*, 'Czas normy: ', clock_readwall(timer)
              call date_and_time(b(1), b(2), b(3), date_time)
              print*, 'date: ', date_time(1), '-', date_time(2), '-', date_time(3)
              print*, 'time: ', date_time(5), 'h ', date_time(6), 'min ', date_time(3), 'sec'

              deallocate(MxA)
              deallocate(w)
              deallocate(wi)
              deallocate(ltz)
              deallocate(rtz)

            end associate

      end subroutine pp_erpa


      subroutine ppAC_int_w(HType, ECorr_AC0, &
            ETot, ENuc, n, XOne, TwoNO, AuxData, IndAux, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, version, Flags, AC_TYPE)

            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: HType
            double precision, intent(inout) :: ECorr_AC0, ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type(TDA_pphh) :: Tpphh
            type (tclock) :: timer, timer0
            integer, intent(in) :: version
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: AC_TYPE
            !double precision, dimension(:, :), intent(out) :: eorb
            double precision, allocatable :: eorb(:,:)
            double precision :: Acpp_alpha_ext, Acpp_alpha_ext_wo
            double precision :: E_contr_s, E_contr_t
            double precision, dimension(25) :: E_corr_contr_s, E_corr_contr_t, E_corr_contr_ext

            integer :: spin_symm

            double precision :: dup
            integer :: i
            allocate(eorb(NBasis, NBasis))


            E_contr_s = zero
            E_corr_contr_s = zero
            print*, 'SINGLET PART'
            print*, ''
            call clock_start(timer0)
            call Omega_contr(HType, E_contr_s, E_corr_contr_s, AuxData%NDim_s,  AuxData%IndN_s, AuxData%IndX_s, &
                  ETot, ENuc,  n, XOne, TwoNO, IndAux, IndMod, &
                  NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, version, Flags, AuxData, AC_TYPE, eorb, 0)
            print*, 0, 'Czas na omega contribution singlet: ', clock_readwall(timer0)

            write(*,'(A25, 2F20.15)') 'E_corr_contr(oooo) =, ', E_corr_contr_s(BL_OOOO)
            write(*,'(A25, 2F20.15)') 'E_corr_contr(vava) =, ', E_corr_contr_s(BL_VAVA)

            E_contr_t = zero
            E_corr_contr_t = zero
            print*,'TRIPLET PART'
            print*,	''
            if (AuxData%NDim_t.gt.0)then
                  call clock_start(timer0)
                  call Omega_contr(HType, E_contr_t, E_corr_contr_t, AuxData%NDim_t,  AuxData%IndN_t, AuxData%IndX_t, &
                        ETot, ENuc,  n, XOne, TwoNO, IndAux, IndMod, &
                        NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, version, Flags, AuxData, AC_TYPE, eorb, 1)
                  print*, 0, 'Czas na omega contribution triplet: ', clock_readwall(timer0)
            end if

            ! to jest wyzerowane bo jest kasowanie analityczne z czlonem Y^0
            ! E_corr_contr_ext = zero
            !   call clock_start(timer0)
            !   call ACPPENeERPA_part_ext2(HType, Acpp_alpha_ext, E_corr_contr_ext, AuxData%NDim_s, AuxData%IndN_s, IndAux, TwoNO, n, AC_TYPE)
            !   print*, 0, 'Czas na omega contribution ext: ', clock_readwall(timer0)

            print*, 'E_contr_s', E_contr_s
            print*, 'E_contr_t', E_contr_t
            print*, 'Acpp_alpha_sum', E_contr_s+E_contr_t
            !print*, 'Acpp_alpha_ext', ACpp_alpha_ext
            ECorr_AC0 = E_contr_s + E_contr_t !+ ACpp_alpha_ext       
            write(*,'(A10, F20.15)') 'ECORR_AC0pp', E_contr_s + E_contr_t !+ ACpp_alpha_ext
            print*, ''

            call print_contributions(E_corr_contr_s, E_corr_contr_t)


      end subroutine ppAC_int_w



      subroutine Omega_contr(HType, E_contr, E_corr_contr, NDim,  IndN, IndX, &
            ETot, ENuc, n, XOne, TwoNO, IndA, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, version, Flags, AuxData, AC_TYPE, eorb, spin_symm)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: HType
            double precision, intent(inout) :: E_contr, ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO

            integer, intent(in) :: NDim
            integer, dimension(:,:), intent(in) ::  IndN
            double precision, dimension(:, :), intent(out) :: eorb
            integer, intent(in) :: spin_symm
            integer, dimension(:), intent(in)    :: IndA, IndMod, IndX
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type(TDA_pphh) :: Tpphh
            type (tclock) :: timer, timer0, timer2
            integer, intent(in) :: version
            type(FlagsData), intent(in) :: Flags
            type(TACppData) :: AuxData
            integer, intent(in) :: AC_TYPE

            double precision, dimension(:,:), allocatable :: MxA, MxS, MxAf0, MxAf02
            double precision, dimension(:,:), allocatable :: MxA1, MxAf1, A1, A2
            double precision, dimension(:,:), allocatable :: C_prev, C_prev_prev, C_current
            double precision, dimension(:,:), allocatable :: C_cumul
            double precision, dimension(:,:), allocatable :: A1kw
            double precision :: u
            double precision, dimension(:,:), allocatable :: h
            double precision, dimension(:), allocatable :: Eig_s, Eig_t, Eig_s_i, Eig_t_i
            double precision, dimension(:,:), allocatable :: EigvecR_s, EigvecR_t
            double precision, dimension(:), allocatable :: norm_s, norm_t, tempx
            integer, dimension(:), allocatable :: v_plus_s, v_plus_t

            double precision :: Aux1, e_contr2
            integer :: i, j, k, l, p, q, r, s
            double precision :: norm_temp, dd, temp
            logical :: cond
            !
            ! Variables for omega version
            !

            double precision, dimension(:), allocatable :: wcp, xcp
            integer :: ncp
            logical :: convergedcp


            double precision, dimension(:,:), allocatable:: Z1, Z3
            double precision :: wkl1, wkl2, wkl1b, wkl2b
            double precision :: w1, w2, w3, w4

            double precision, dimension(25), intent(inout) :: E_corr_contr
            double precision :: contr, contr2, dup1, dup2

            double precision, parameter :: PIp = 3.141592653589793238462643383279d+0
            integer, external :: NAddr3
            double precision, dimension(:), allocatable :: Aux2
            integer :: ik, jj, nn, Norders
            integer :: twoint_dim
            double precision :: dupo, dup, oint
            integer :: ai, bj

            twoint_dim = size(TWONO, dim=1)
            allocate(Aux2(NDim*NDim))


            allocate(Z1(NDim, NDim))
            allocate(Z3(NDim, NDim))
            allocate(C_cumul(NDim, NDim))
            allocate(MxS(NDim, NDim))
            allocate(tempx(NDim))





            allocate(MxA(NDim, NDim))
            allocate(MxAf0(NDim, NDim))
            allocate(MxAf02(NDim, NDim))
            allocate(C_prev_prev(NDim, NDim))

            allocate(MxA1(NDim, NDim))
            allocate(MxAf1(NDim, NDim))
            allocate(A1(NDim, NDim))
            allocate(A2(NDim, NDim))
            allocate(A1kw(NDim, NDim))


            allocate(C_prev(NDim, NDim))
            allocate(C_current(NDim, NDim))




            ! CALCULATE (u^2 + (~A^(0))^2)


            ncp = 100

            allocate(wcp(ncp))

            allocate(xcp(ncp))

            call quad_CasimirPolder(xcp, wcp, convergedcp, ncp, 1.0d-14)

            ! calculate A(alpha=0) find diff = (max(eigval_hh) - min(eigval_pp))/2

            ! As a result of this procedur we have A^(0) and S^(-1) matrices
            call clock_start(timer)
            call pp_rowe(AuxData, HType, ETot, ENuc, MxAf0, MxS,  n, XOne, TwoNO, IndN, IndX, IndA, IndMod, NDim, &
                  NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, zero, spin_symm, version, Flags, eorb)
            print*, ''
            print*, 'time for pperpa alpha=0: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)

            !         do ai = 1, size(MxAf0, dim=1)
            !              do bj = 1, size(MxAf0, dim=1)
            !                 if (abs(MxAf0(ai, bj)).gt.1.d+1)then
            ! 	           print*, 'maxf0', ai, bj, MxAf0(ai, bj)
            !                 end if
            !              end do
            !           end do
            ! stop

            ! As a result of this procedur we have A^(1) and S^(-1) matrices
            call pp_rowe(AuxData, HType, ETot, ENuc, MxAf1, MxS,  n, XOne, TwoNO, IndN, IndX, IndA, IndMod, NDim, &
                  NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, one, spin_symm, version, Flags, eorb)
            print*, ''
            print*, 'time for pperpa alpha=1: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)


            MxAf1 = MxAf1 - MxAf0


            print*, ''
            print*, 'bede liczyc macierze', spin_symm, NDim

            ! calculate MxAf02=(A~0 )^2
            call real_ab(MxAf02, MxAf0, MxAf0)
            print*, ''
            print*, 'time for MxAf02^2: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)


            ! (u^2+((A~)^(0))^2) C_prev = A^(1) - ((A~)^(0)(A~)^(1) + (A~)^(1)(A~)^(0))C_prev_prev

            ! calculate A1 = (A~)^(0)(A~)^(1)
            call real_ab(A1, MxAf0, MxAf1)
            print*, ''
            print*, 'time for A0*A1: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)


            ! calculate A2 = (A~)^(1)(A~)^(0)
            call real_ab(A2, MxAf1, MxAf0)
            print*, ''
            print*, 'time for A1*A0: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)

            ! calculate A1kw = (A~)^(1)(A~)^(1)
            call real_ab(A1kw, MxAf1, MxAf1)
            print*, ''
            print*, 'time for A1kw: ', clock_readwall(timer)
            print*, ''
            call clock_start(timer)



            ! calculate (A~)^(0)(A~)^(1) + (A~)^(1)(A~)^(0))C_prev_prev
            ! (A1 + A2) * C_prev_prev <- saved in A2
            A1 = A1 + A2
            A2 = zero

            E_contr = zero
            wkl1b = zero

            wkl2b= zero
            w1 = zero
            w2 = zero
            w3 = zero
            w4 = zero
            oint = zero
            dupo = zero

            print*,''
            print*, 'wklad oooo poczatek', E_corr_contr(BL_OOOO)
            print*, 'wklad vava poczatek', E_corr_contr(BL_VAVA)
            print*, 'wklad vvoo poczatek', E_corr_contr(BL_VVOO)       
            print*, 'wklad oovv poczatek', E_corr_contr(BL_OOVV)

            !       Norders = 2
            Norders = AuxData%omegaorders
            print*, 'RDSC mbpt order', Norders

            call clock_start(timer)

            oint = zero

            !$omp parallel do &
            !$omp default(shared) &
            !$omp private(l, u, jj ,nn) &
            !$omp private(Z1, Z3, C_prev_prev, A2, A1, C_prev, C_current, C_cumul, E_contr) &
            !$omp private(i, j, p, q, r, s, cond, contr, Aux1)
            int_loop: do l = 1, ncp
                  u = xcp(l)

                  Z3 = MxAf02
                  do jj = 1, NDim
                        Z3(jj,jj) = u**2 + Z3(jj,jj)
                  end do


                  ! (MxAf02) .C_prev_prev = MxAf0  solve for C_prev_prev
                  C_prev_prev = MxAf0
                  Z1 = Z3

                  call real_Axb_nonsymmetric_gesv(C_prev_prev, Z1) !C_prev_prev is stored in MxA


                  call real_ab(A2, A1, C_prev_prev)

                  ! (A~)^(1) - ((A~)^(0)(A~)^(1) + (A~)^(1)(A~)^(0))C_prev_prev saved in C_prev
                  C_prev = MxAf1 - A2
                  Z1 = Z3

                  call real_Axb_nonsymmetric_gesv(C_prev, Z1) !C_prev_prev is stored in MxA

                  !          C_current = C_prev / two

                  C_current = C_prev
                  C_cumul = C_prev/two

                  do nn = 2,  Norders
                        A2 = zero
                        C_current = zero

                        call real_ab(C_current, A1, C_prev)
                        call real_ab(A2, A1kw, C_prev_prev)


                        Z1 = C_current
                        C_current = zero
                        C_current = -nn*Z1 -nn*(nn-1) * A2
                        Z1 = Z3

                        call real_Axb_nonsymmetric_gesv(C_current, Z1)
                        C_cumul = C_cumul + C_current / (factorial_tab(nn) * (nn+1))

                        !             C_cumul = C_current / (factorial_tab(nn) * (nn+1))

                        C_prev_prev = C_prev                                                                                                                        
                        C_prev = C_current  
                  end do



                  pq_loop: do i = 1, Ndim
                        rs_loop: do j = 1, NDim
                              p = IndN(1, i)
                              q = IndN(2, i)
                              r = IndN(1, j)
                              s = IndN(2, j)

                              if (HType == H_DYALL)then
                                    cond = ((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)).and.(IndA(p)==1))
                              else if (HType == H_GPF.or.HType==H_MID)then
                                    cond = ((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)))
                              end if


                              if (.not.(cond))then

                                    if (spin_symm == 0)then
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))

                                          if (p==q)then
                                                Aux1 = Aux1 * sqrt(frac12)
                                          end if
                                          if (r==s) then
                                                Aux1 = Aux1 * sqrt(frac12)
                                          end if



                                          !contr = (one/Pi) * wcp(l)* Aux1* MxS(i,i)*(C_prev_prev(i,j)) + frac12 * C_prev(i,j))
                                          contr = (one/Pi) * wcp(l)* Aux1* MxS(i,i)*C_cumul(i,j) !*frac12 !<---to dwa do sprawdzenia
                                          ! if (p==6.and.q==6.and.r==14.and.s==14)then
                                          !    write(*,'(A10, 4F20.12)') 'szesc', Aux1, wcp(l), MxS(i,i), C_cumul(i,j)
                                          ! end if
                                          ! if (r==6.and.s==6.and.p==14.and.q==14)then
                                          !    write(*,'(A10, 4F20.12)') 'czter', Aux1, wcp(l), MxS(i,i), C_cumul(i,j)
                                          ! end if
                                          !                      oint = oint + wcp(l) * C_cumul(i,j)

                                          E_contr = E_contr + contr
                                          !$OMP CRITICAL 
                                          call update_E_corr_contribution(contr, p, q, r, s, IndA, E_corr_contr)
                                          !$OMP end CRITICAL 

                                    else
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
                                          !contr = Three/Pip * wcp(l)* Aux1*MxS(i,i)*(C_prev_prev(i,j)) + frac12 * C_prev(i,j))
                                          contr = Three/Pip * wcp(l)* Aux1*MxS(i,i)*C_cumul(i,j) !*frac12 !<---to dwa do sprawdzenia

                                          !$omp critical
                                          E_contr = E_contr + contr
                                          call update_E_corr_contribution(contr, p, q, r, s, IndA, E_corr_contr)
                                          !$OMP end CRITICAL 
                                    end if
                              end if

                        end do rs_loop
                  end do pq_loop

                  ! print*, ''
                  ! print*, 'time for pqrs loops: ', clock_readwall(timer0)
                  ! print*, ''
                  !         stop
                  ! call clock_start(timer0)
                  ! print*,''
                  ! print*, 'wklad oooo', E_corr_contr(BL_OOOO)
                  ! print*,''
                  ! print*, 'wklad vava', E_corr_contr(BL_VAVA)
                  ! print*,''
                  ! print*, 'wklad vvvv', E_corr_contr(BL_VVVV)
                  ! print*,''
                  ! print*, 'wklad oovv', E_corr_contr(BL_OOVV)
                  ! print*,''
                  ! print*, 'wklad vvoo', E_corr_contr(BL_VVOO)



            end do int_loop
            !$omp end parallel do

            print*, 'wklad oooo koniec', E_corr_contr(BL_OOOO)



            print*, ''
            print*, 'e_contr1', e_contr

            print*, 'Czas omega: ', clock_readwall(timer)


      end subroutine Omega_contr

      subroutine ppAC_int(HType, E_corr_contr_s, E_corr_contr_t, &
            E_corr_contr_ext, Acpp_alpha_s, Acpp_alpha_t, Acpp_alpha_s_wo, Acpp_alpha_t_wo, &
            Acpp_alpha_ext, Acpp_alpha_ext_wo, &
            ETot, ENuc,  n, XOne, TwoNO, AuxData, IndAux, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, ACAlpha, version, Flags, AC_TYPE, eorb)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: HType
            double precision, intent(inout) :: Acpp_alpha_s, Acpp_alpha_t
            double precision, intent(inout) :: Acpp_alpha_s_wo, Acpp_alpha_t_wo
            double precision, intent(inout) :: Acpp_alpha_ext, Acpp_alpha_ext_wo
            double precision, intent(inout) :: ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type(TDA_pphh) :: Tpphh
            type (tclock) :: timer
            integer, intent(in) :: version
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: AC_TYPE
            integer :: spin_symm
            integer :: hhmet

            double precision, dimension(:), intent(inout) :: E_corr_contr_s, E_corr_contr_t, E_corr_contr_ext
            double precision, dimension(:,:), allocatable :: MxA, MxS, MxSs, MxSt
            double precision, dimension(:,:), allocatable :: h
            double precision, dimension(:), allocatable :: Eig_s, Eig_t, Eig_s_i, Eig_t_i
            double precision, dimension(:,:), allocatable :: EigvecR_s, EigvecR_t
            double precision, dimension(:), allocatable :: norm_s, norm_t, tempx
            integer, dimension(:), allocatable :: v_plus_s, v_plus_t

            double precision :: ECorr
            integer :: i, j, k, p, q, r, s
            double precision :: norm_temp, dd
            double precision, dimension(:, :), intent(out) :: eorb
            hhmet = 0


            allocate(MxA(AuxData%NDim_s, AuxData%NDim_s))
            allocate(MxSs(AuxData%NDim_s, AuxData%NDim_s))
            allocate(Eig_s(AuxData%NDim_s))
            allocate(Eig_s_i(AuxData%NDim_s))
            allocate(EigvecR_s(AuxData%NDim_s, AuxData%NDim_s))
            allocate(tempx(AuxData%NDim_s))
            allocate(v_plus_s(AuxData%Ndim_s))

            spin_symm=0
            call clock_start(timer)
            call pp_rowe(AuxData, HType, ETot, ENuc, MxA, MxSs,  n, XOne, TwoNO, AuxData%IndN_s, &
                  AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, &
                  NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, ACAlpha, spin_symm, version, Flags, eorb)
            print*, 'Czas pperpa sing: ', clock_readwall(timer)
!            print*, MxSs


            call clock_start(timer)
            call nonsymmetric_eigenproblem_pp(Eig_s, Eig_s_i, EigvecR_s, MxA, MxSs, v_plus_s,  AuxData%IndN_s, IndAux, AC_TYPE)
            print*, 'Czas diagonalizacji sing: ', clock_readwall(timer)
 !           print*, 'eigssssss'
  !          print*, eig_s

            deallocate(MxA)
            !       deallocate(MxS)
            deallocate(tempx)

            allocate(MxA(AuxData%NDim_t, AuxData%NDim_t))
            allocate(MxSt(AuxData%NDim_t, AuxData%NDim_t))
            allocate(Eig_t(AuxData%NDim_t))
            allocate(Eig_t_i(AuxData%NDim_t))
            allocate(EigvecR_t(AuxData%NDim_t, AuxData%NDim_t))
            allocate(tempx(AuxData%Ndim_t))
            allocate(v_plus_t(AuxData%Ndim_t))
            spin_symm = 1

            Acpp_alpha_t = zero
            if (AuxData%NDim_t.gt.0)then
                  call clock_start(timer)
                  call pp_rowe(AuxData, HType, ETot, ENuc, MxA, MxSt,  n, XOne, TwoNO, AuxData%IndN_t, &
                        AuxData%IndX_t, IndAux, IndMod, AuxData%NDim_t, &
                        NBasis,  NA, NI, NV, NInte1, NInte2, Tpphh, ACAlpha, spin_symm, version, Flags, eorb)
                  print*, 'Czas pperpa_trip: ', clock_readwall(timer)



                  call clock_start(timer)
                  call nonsymmetric_eigenproblem_pp(Eig_t, Eig_t_i, EigvecR_t, MxA, MxSt, v_plus_t,  AuxData%IndN_t, IndAux, AC_TYPE)
                  print*, 'Czas diagonalizacji trip ', clock_readwall(timer)
!                  print*, 'eigt', Eig_t

                  deallocate(MxA)
                  !          deallocate(MxS)


                  call erpa_trip(HType, E_corr_contr_t, Acpp_alpha_t, Eig_t, Eig_t_i, EigvecR_t, EigvecR_t, v_plus_t, &
                        AuxData%NDim_t, AuxData%IndN_t, AuxData%IndAux, MxSt, TwoNO, n, AC_TYPE, 100)
!                  print*, 'E_corr_contr_t', E_corr_contr_t

            end if

            call erpa_sing(HType, E_corr_contr_s, Acpp_alpha_s, Eig_s, Eig_s_i, EigvecR_s, EigvecR_s, v_plus_s, &
                  AuxData%NDim_s, AuxData%IndN_s, AuxData%IndAux, MxSs, TwoNO, n, AC_TYPE, 100)

            print*, 'herehere', ACpp_alpha_s
            call erpa_ext(HType, E_corr_contr_ext, &
                  Acpp_alpha_ext, Acpp_alpha_ext_wo, AuxData%NDim_s, AuxData%IndN_s, IndAux, TwoNO, n, AC_TYPE)

            Acpp_alpha_t_wo = zero
            Acpp_alpha_s_wo = zero

            deallocate(EigvecR_s)
            deallocate(EigvecR_t)
            deallocate(Eig_s)
            deallocate(Eig_t)
            deallocate(Eig_s_i)
            deallocate(Eig_t_i)


      end subroutine ppAC_int

      subroutine erpa_sing(HType, E_corr_contr, ECorr, Eig, Eig_i, EigvecR, EigvecR1, &
            v_plus, NDim, IndN, IndA, MxSs, TwoNO, Occ, AC_TYPE, BlockDim)

            integer, intent(in) :: HType
            double precision, dimension(:), intent(inout) :: E_corr_contr
            double precision, intent(out) :: ECorr
            double precision, dimension(:), intent(in) :: Eig, Eig_i
            integer, dimension(:), intent(in) :: v_plus
            double precision, dimension(:,:), intent(in) :: EigvecR, EigvecR1
            integer, intent(in) :: NDim
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), intent(in) :: IndA
            double precision, dimension(:,:), intent(in) :: MxSs
            double precision, dimension(:),      intent(in) :: TwoNO
            double precision, dimension(:),    intent(in)   :: Occ
            integer, intent(in) :: AC_TYPE
            integer, intent(in) :: BlockDim


            integer :: i, j, k
            integer :: p, q, r, s

            double precision :: SumZ, Aux1, Aux2, Aux3, ecor_aux
            double precision :: pppp, pqpq, pqqp, pqpq_p, pqqp_p
            double precision :: mt4_part, mt4
            type (tclock) :: start, timer

            double precision, parameter :: SmallE=1.d-3, BigE=1.d+3
            integer :: iskipped, pq, rs
            logical :: cond, cond_clos
            double precision :: contr_vv_vv, contr_vv_oo, contr_oo_oo, contr_oo_vv
            double precision :: contr_vv_vv_this, contr_vv_oo_this, contr_oo_oo_this, contr_oo_vv_this
            double precision :: contr_this
            double precision, dimension(:, :), allocatable :: Vl, Vr
            double precision, dimension(:, :), allocatable :: SS
            integer :: cc
            integer :: NBlocks, NVecs
            integer :: i0, i1, j0, j1, bi, bj
            integer :: Ni, Nj
            integer, external :: NAddr3

            integer, dimension(9) :: lst

            logical ::ism, pf, qf


            lst(1) = 3
            lst(2)=5
            lst(3)=12
            lst(4)=13
            lst(5)=18
            lst(6)=21
            lst(7)=24
            lst(8)=26
            lst(9)=29


            call clock_start(timer)
            print *, "------ AC PP ERPA Block algorithm -----"
            NVecs = 0
            do k = 1, NDim
                  if (v_plus(k) == AC_TYPE) NVecs = NVecs + 1
            end do
            allocate(Vl(BlockDim, NVecs))
            allocate(Vr(BlockDim, NVecs))
            allocate(SS(BlockDim, BlockDim))
            NBlocks = NDim / BlockDim
            if (modulo(NDim, BlockDim) > 0) NBlocks = NBlocks + 1
            ECorr = zero
            do bi = 1, NBlocks
                  do bj = 1, NBlocks
                        i0 = 1 + BlockDim * (bi - 1)
                        i1 = min(BlockDim * bi, NDim)
                        j0 = 1 + BlockDim * (bj - 1)
                        j1 = min(BlockDim * bj, NDim)
                        Ni = i1 - i0 + 1
                        Nj = j1 - j0 + 1
                        NVecs = 0
                        if (Ni < NDim) Vl = ZERO
                        if (Nj < NDim) Vr = ZERO
                        do k = 1, NDim
                              if (v_plus(k) == AC_TYPE) then
                                    NVecs = NVecs + 1
                                    Vl(1:Ni, NVecs) = EigvecR(i0:i1, k)
                                    Vr(1:Nj, NVecs) = EigvecR1(j0:j1, k)
                              end if
                        end do
                        !
                        ! SS(1:Ni,1:Nj) = Sum(k) Vl(i0:i1), 1:NVecs) Vr(j0:j1, 1:NVecs)**T
                        !
                        call real_abT(SS, Vl, Vr)                      
                        do i = i0, i1
                              do j = j0, j1
                                    p = IndN(1, i)
                                    q = IndN(2, i)
                                    r = IndN(1, j)
                                    s = IndN(2, j)
                                    !                                  ism = ismixed(lst, p,q,r,s)
                                    if (HType == H_DYALL)then
                                          cond = ((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)).and.(IndA(p)==1))
                                    else if (HType == H_GPF.or.HType==H_MID)then
                                          cond =	((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)))
                                    end if
                                    if (.not.(cond))then
                                          SumZ = SS(i-i0+1, j-j0+1)
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
                                          if (p==q)then
                                                Aux1 = Aux1 * sqrt(frac12)
                                          end if
                                          if (r==s) then
                                                Aux1 = Aux1 * sqrt(frac12)
                                          end if
                                          ! if (ism)then
                                          !    if (abs(Aux1).gt.1.d-6.and.abs(SumZ).gt.1.d-6)then
                                          !       print*, 'ISMIXED', p, q, r, s, Aux1, SumZ
                                          !    end if
                                          ! end if
                                          !contr_this = Aux1* SumZ * MxSs(i,j)*MxSs(i,j)!( Occ(p)+Occ(q)-one )  * ( Occ(r)+Occ(s)-one )
                                          contr_this = Aux1* SumZ * ( Occ(p)+Occ(q)-one )  * ( Occ(r)+Occ(s)-one )
                                          !write(*, '(A20, 2F20.15)')'pluszsing', Aux1, sumz
                                          ECorr = Ecorr +  contr_this
                                          call update_E_corr_contribution(contr_this, p, q, r, s, IndA, E_corr_contr)
                                    end if
                              end do
                        end do
                  end do
            end do
            
            print *, "czas na parts ", clock_readwall(timer)
      end subroutine erpa_sing

      subroutine erpa_trip(HType, E_corr_contr, ECorr, Eig, Eig_i, EigvecR, EigvecR1, &
            v_plus, NDim, IndN, IndA, MxSt, TwoNO, Occ, AC_TYPE, BlockDim)

            integer, intent(in) :: HType
            double precision, intent(out) :: ECorr
            double precision, dimension(:), intent(inout) :: E_corr_contr
            double precision, dimension(:), intent(in) :: Eig, Eig_i
            integer, dimension(:), intent(in) :: v_plus
            double precision, dimension(:,:), intent(in) :: EigvecR, EigvecR1
            integer, intent(in) :: NDim
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), intent(in) :: IndA
            double precision, dimension(:,:), intent(in) :: MxSt
            double precision, dimension(:),      intent(in) :: TwoNO
            double precision, dimension(:),    intent(in)   :: Occ
            integer, intent(in) :: AC_TYPE
            integer, intent(in) :: BlockDim


            integer :: i, j, k
            integer :: p, q, r, s

            double precision :: SumZ, Aux1, Aux2, Aux3, ecor_aux
            double precision :: pppp, pqpq, pqqp, pqpq_p, pqqp_p
            double precision :: mt4_part, mt4
            type (tclock) :: start, timer

            double precision, parameter :: SmallE=1.d-3, BigE=1.d+3
            integer :: iskipped, pq, rs
            logical :: cond, cond_clos
            double precision :: contr_vv_vv, contr_vv_oo, contr_oo_oo, contr_oo_vv
            double precision :: contr_vv_vv_this, contr_vv_oo_this, contr_oo_oo_this, contr_oo_vv_this
            double precision :: contr_this
            double precision, dimension(:, :), allocatable :: Vl, Vr
            double precision, dimension(:, :), allocatable :: SS
            integer :: cc
            integer :: NBlocks, NVecs
            integer :: i0, i1, j0, j1, bi, bj
            integer :: Ni, Nj
            integer, external :: NAddr3

            call clock_start(timer)
            print *, "------ AC PP ERPA Block algorithm -----"
            NVecs = 0
            do k = 1, NDim
                  if (v_plus(k) == AC_TYPE) NVecs = NVecs + 1
            end do
            allocate(Vl(BlockDim, NVecs))
            allocate(Vr(BlockDim, NVecs))
            allocate(SS(BlockDim, BlockDim))
            NBlocks = NDim / BlockDim
            if (modulo(NDim, BlockDim) > 0) NBlocks = NBlocks + 1
            ECorr = zero
            do bi = 1, NBlocks
                  do bj = 1, NBlocks
                        i0 = 1 + BlockDim * (bi - 1)
                        i1 = min(BlockDim * bi, NDim)
                        j0 = 1 + BlockDim * (bj - 1)
                        j1 = min(BlockDim * bj, NDim)
                        Ni = i1 - i0 + 1
                        Nj = j1 - j0 + 1
                        NVecs = 0
                        if (Ni < NDim) Vl = ZERO
                        if (Nj < NDim) Vr = ZERO
                        do k = 1, NDim
                              if (v_plus(k) == AC_TYPE) then
                                    NVecs = NVecs + 1
                                    Vl(1:Ni, NVecs) = EigvecR(i0:i1, k)
                                    Vr(1:Nj, NVecs) = EigvecR1(j0:j1, k)
                              end if
                        end do
                        !
                        ! SS(1:Ni,1:Nj) = Sum(k) Vl(i0:i1), 1:NVecs) Vr(j0:j1, 1:NVecs)**T
                        !
                        call real_abT(SS, Vl, Vr)                      
                        do i = i0, i1
                              do j = j0, j1
                                    p = IndN(1, i)
                                    q = IndN(2, i)
                                    r = IndN(1, j)
                                    s = IndN(2, j)
                                    if (HType == H_DYALL)then
                                          cond = ((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)).and.(IndA(p)==1))
                                    else if (HType == H_GPF.or.HType==H_MID)then
                                          cond =	((IndA(p)==IndA(q)).and.(IndA(p)==IndA(r)).and.(IndA(p)==IndA(s)))
                                    end if
                                    if (.not.(cond))then
                                          SumZ = SS(i-i0+1, j-j0+1)
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q)) -  TwoNO(NAddr3(r,q,s,p)))
                                          !contr_this = three * Aux1* SumZ * MxSt(i,j)*MxSt(i,j)!( Occ(p)+Occ(q)-one )  * ( Occ(r)+Occ(s)-one )
                                          contr_this = three * Aux1* SumZ * ( Occ(p)+Occ(q)-one )  * ( Occ(r)+Occ(s)-one )
                                          ECorr = Ecorr +  contr_this
                                          !write(*, '(A20, 2F20.15)')'plusztrip', Aux1, sumz
                                          call update_E_corr_contribution(contr_this, p, q, r, s, IndA, E_corr_contr)
                                    end if
                              end do
                        end do
                  end do
            end do
            print *, "czas na partt ",clock_readwall(timer)
      end subroutine erpa_trip

      subroutine erpa_ext(HType, E_corr_contr, ecor_aux, ecor_aux_wo, NDim, IndN, IndA, TwoNO, Occ, AC_TYPE)

            integer, intent(in) :: HType
            double precision, intent(out) :: ecor_aux, ecor_aux_wo
            double precision, dimension(:), intent(inout) :: E_corr_contr
            integer, intent(in) :: NDim
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), intent(in) :: IndA
            double precision, dimension(:),      intent(in) :: TwoNO
            double precision, dimension(:),    intent(in)   :: Occ
            integer, intent(in) :: AC_TYPE

            logical :: cond, cond_clos
            double precision :: contr_vv, contr_oo, contr_this = zero
            double precision :: contr_vv_this, contr_oo_this
            integer :: i, j, k
            integer :: p, q, r, s

            double precision :: SumZ, Aux1, Aux2, Aux3
            double precision, parameter :: SmallE=1.d-3, BigE=1.d+8
            integer :: iskipped
            double precision, dimension(:), allocatable :: skipped
            double precision :: pp, pq, pq1, qp1, qp, pq_this, qp_this

            integer, external :: NAddr3


            ecor_aux = zero
            ecor_aux_wo = zero
            contr_vv = zero
            contr_oo = zero
            contr_vv_this = zero
            contr_oo_this = zero

            contr_this = zero

            pp = zero
            pq = zero
            qp = zero

            do i = 1, Ndim
                  p = IndN(1, i)
                  q = IndN(2, i)            

                  if (HType == H_DYALL)then
                        cond = ((IndA(p)==IndA(q)).and.(IndA(p)==1))
                  else if (HType == H_GPF.or.HType==H_MID)then
                        cond = ((IndA(p)==IndA(q)))
                  end if

                  if (.not.(cond))then
                        !               print*, 'yes1'
                        if (AC_TYPE==PPAC)then
                              Aux2 = -two * (1-Occ(q)) * (1-Occ(p)) *( two* TwoNO(NAddr3(p,p,q,q)) -  TwoNO(NAddr3(p,q,p,q)))
                              pq_this = -two * (1-Occ(q)) * (1-Occ(p)) *( two* TwoNO(NAddr3(p,p,q,q)))
                              qp_this = two * (1-Occ(q)) * (1-Occ(p)) *( TwoNO(NAddr3(p,q,p,q)))
                        else if (AC_TYPE == HHAC)then
                              Aux2 = -two * Occ(q) * Occ(p) *( two* TwoNO(NAddr3(p,p,q,q)) -  TwoNO(NAddr3(p,q,p,q)))
                        end if

                        if (p==q)then
                              Aux2 = frac12 * Aux2
                        end if


                        contr_vv_this = zero
                        contr_oo_this = zero

                        !if (abs(Aux2).gt.1.d-8)then

                        ! if (p.gt.nel2.and.q.gt.nel2)then
                        !    contr_vv = contr_vv + Aux2
                        !    contr_vv_this =  Aux2
                        !    !write(*, '(A10, 2I5, F20.10)') 'qvv', p, q, contr_vv_this
                        ! else if (p.le.nel2.and.q.le.nel2)then
                        !    contr_oo = contr_oo + Aux2
                        !    contr_oo_this = Aux2
                        !     !write(*, '(A10, 2I5, F20.10)') 'qoo', p, q,contr_oo_this
                        ! end if

                        !               print*, 'contr',contr_vv_this , contr_oo_this
                        if ((one-Occ(p)-Occ(q)).gt.0.d+0)then
                              contr_vv = contr_vv + Aux2
                              contr_vv_this =  Aux2

                        else if ((one-Occ(p)-Occ(q)).lt.0.d+0)then
                              contr_oo = contr_oo + Aux2
                              contr_oo_this = Aux2

                        end if


                        ecor_aux = ecor_aux + contr_vv_this + contr_oo_this!Aux2
                        call update_E_corr_contribution(contr_vv_this + contr_oo_this, p, q, p, q, IndA, E_corr_contr)
                        ecor_aux_wo = ecor_aux_wo
                  end if
            end do

            print*, 'auxllll', ecor_aux
            !          print*, 'auxllll_wo', ecor_aux_wo

            !       print*, 'contr_vve', contr_vv
            !	print*, 'contr_ooe', contr_oo

      end subroutine erpa_ext


      subroutine update_E_corr_contribution(contr, p, q, r, s, IndA, E_corr_contr)
            double precision, intent(in) :: contr
            integer, intent(in) :: p, q, r, s
            integer, dimension(:), intent(in) :: IndA
            double precision, dimension(25), intent(inout) :: E_corr_contr
            integer :: pq_type, rs_type

            pq_type = IndA(p) * 10 + IndA(q)
            rs_type = IndA(r) * 10 + IndA(s)


            ! write(*,'(A5, 6I4)') 'pqrs', IndA(p), IndA(q), IndA(r), IndA(s), pq_type, rs_type

            select case(pq_type)
            case(0) ! pq is oo
                  select case(rs_type)
                  case(0) ! rs is oo
                        E_corr_contr(BL_OOOO) = E_corr_contr(BL_OOOO) + contr
                        ! if (abs(E_corr_contr(oooo)).gt.1.d+0)then
                        !    print*, 'E_corr_contr(oooo)', E_corr_contr(oooo)
                        ! end if
                  case(10) ! rs is ao
                        E_corr_contr(BL_OOAO) = E_corr_contr(BL_OOAO) + contr
                  case(11) ! rs is aa
                        E_corr_contr(BL_OOAA) = E_corr_contr(BL_OOAA) + contr
                        ! if(abs(contr).gt.1.d-5)then
                        !    write(*, '(A25, 4I5, 2F20.12)') 'E_corr_contr(ooaa)', p,q,r,s,contr, E_corr_contr(BL_OOAA)
                        ! end if
                  case(21) ! rs is va
                        E_corr_contr(BL_OOVA) = E_corr_contr(BL_OOVA) + contr
                  case(22) ! rs is vv
                        E_corr_contr(BL_OOVV) = E_corr_contr(BL_OOVV) + contr
                  case default
                        print*, 'nie wpisuje0'

                  end select

            case(10) ! pq is ao
                  select case(rs_type)
                  case(0) ! rs is oo
                        E_corr_contr(BL_AOOO) = E_corr_contr(BL_AOOO) + contr
                  case(10) ! rs is ao
                        E_corr_contr(BL_AOAO) = E_corr_contr(BL_AOAO) + contr
                  case(11) ! rs is aa
                        E_corr_contr(BL_AOAA) = E_corr_contr(BL_AOAA) + contr
                  case(21) ! rs is va
                        E_corr_contr(BL_AOVA) = E_corr_contr(BL_AOVA) + contr
                  case(22) ! rs is vv
                        E_corr_contr(BL_AOVV) = E_corr_contr(BL_AOVV) + contr
                  case default
                        print*, 'nie wpisuje1'

                  end select

            case(11) ! pq is aa
                  select case(rs_type)
                  case(0) ! rs is oo
                        E_corr_contr(BL_AAOO) = E_corr_contr(BL_AAOO) + contr
                        ! if(abs(contr).gt.1.d-5)then
                        !    write(*, '(A25, 4I5, 2F20.12)') 'E_corr_contr(aaoo)', p,q,r,s,contr, E_corr_contr(BL_AAOO)
                        ! end if
                  case(10) ! rs is ao
                        E_corr_contr(BL_AAAO) = E_corr_contr(BL_AAAO) + contr
                  case(11) ! rs is aa
                        E_corr_contr(BL_AAAA) = E_corr_contr(BL_AAAA) + contr
                  case(21) ! rs is va
                        E_corr_contr(BL_AAVA) = E_corr_contr(BL_AAVA) + contr
                  case(22) ! rs is vv
                        E_corr_contr(BL_AAVV) = E_corr_contr(BL_AAVV) + contr
                  case default
                        print*, 'nie wpisuje2'

                  end select

            case(21) ! pq is va
                  select case(rs_type)
                  case(0) ! rs is oo
                        E_corr_contr(BL_VAOO) = E_corr_contr(BL_VAOO) + contr
                  case(10) ! rs is ao
                        E_corr_contr(BL_VAAO) = E_corr_contr(BL_VAAO) + contr
                  case(11) ! rs is aa
                        E_corr_contr(BL_VAAA) = E_corr_contr(BL_VAAA) + contr
                  case(21) ! rs is va
                        E_corr_contr(BL_VAVA) = E_corr_contr(BL_VAVA) + contr
                        ! if (abs(contr).gt.1.d-7)then
                        !    if (q.ne.s)then
                        !       if (s==3.or.s==5)then
                        !          write(*,'(A20, 4I5)') 'ahha', p, q, r, s
                        !       end if
                        !    end if
                        ! end if
                        !                            if (IndAux(p)==2.and.q==6.and.IndAux(r)==2.and.s==3)then
                  case(22) ! rs is vv
                        E_corr_contr(BL_VAVV) = E_corr_contr(BL_VAVV) + contr
                  case default
                        print*, 'nie wpisuje3'

                  end select

            case(22) ! pq is vv
                  select case(rs_type)
                  case(0) ! rs is oo                      
                        E_corr_contr(BL_VVOO) = E_corr_contr(BL_VVOO) + contr

                  case(10) ! rs is ao
                        E_corr_contr(BL_VVAO) = E_corr_contr(BL_VVAO) + contr
                  case(11) ! rs is aa
                        E_corr_contr(BL_VVAA) = E_corr_contr(BL_VVAA) + contr
                  case(21) ! rs is va
                        E_corr_contr(BL_VVVA) = E_corr_contr(BL_VVVA) + contr
                  case(22) ! rs is vv
                        E_corr_contr(BL_VVVV) = E_corr_contr(BL_VVVV) + contr
                  case default
                        print*, 'nie wpisuje4'
                  end select

            case default
                  ! Handle unexpected case
                  print*, 'nie wpisujea', IndA(p), IndA(q), IndA(r), IndA(s), pq_type, rs_type
            end select


      end subroutine update_E_corr_contribution



      subroutine print_contributions(E_corr_contr_s, E_corr_contr_t)
            double precision, dimension(:), intent(in) :: E_corr_contr_s, E_corr_contr_t
            integer :: i
            double precision :: sum_contr

            sum_contr = zero
            do i = 1, 25
                  sum_contr = sum_contr + E_corr_contr_s(i)+E_corr_contr_t(i)

            end do
            write(*,'(A40, F25.15)') 'E_corr as sum of contributions', sum_contr

            write(*,'(A25, 2F25.15)') 'E_corr_contr(oooo) =, ', E_corr_contr_s(BL_OOOO) + E_corr_contr_t(BL_OOOO) !, E_corr_contr_ext(oooo)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaaa) =, ', E_corr_contr_s(BL_AAAA) + E_corr_contr_t(BL_AAAA) !, E_corr_contr_ext(BL_AAAA)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aoao) =, ', E_corr_contr_s(BL_AOAO) + E_corr_contr_t(BL_AOAO) !, E_corr_contr_ext(BL_AOAO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(vava) =, ', E_corr_contr_s(BL_VAVA) + E_corr_contr_t(BL_VAVA) !, E_corr_contr_ext(BL_VAVA)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvvv) =, ', E_corr_contr_s(BL_VVVV) + E_corr_contr_t(BL_VVVV) !, E_corr_contr_ext(BL_VVVV)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(ooao) =, ', E_corr_contr_s(BL_OOAO) + E_corr_contr_t(BL_OOAO) !, E_corr_contr_ext(ooao)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aooo) =, ', E_corr_contr_s(BL_AOOO) + E_corr_contr_t(BL_AOOO) !, E_corr_contr_ext(BL_AOOO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(ooaa) =, ', E_corr_contr_s(BL_OOAA) + E_corr_contr_t(BL_OOAA) !, E_corr_contr_ext(BL_OOAA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaoo) =, ', E_corr_contr_s(BL_AAOO) + E_corr_contr_t(BL_AAOO) !, E_corr_contr_ext(BL_AAOO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(oova) =, ', E_corr_contr_s(BL_OOVA) + E_corr_contr_t(BL_OOVA) !, E_corr_contr_ext(BL_OOVA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaoo) =, ', E_corr_contr_s(BL_VAOO) + E_corr_contr_t(BL_VAOO) !, E_corr_contr_ext(BL_VAOO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(oovv) =, ', E_corr_contr_s(BL_OOVV) + E_corr_contr_t(BL_OOVV) !, E_corr_contr_ext(BL_OOVV)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvoo) =, ', E_corr_contr_s(BL_VVOO) + E_corr_contr_t(BL_VVOO) !, E_corr_contr_ext(BL_VVOO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aoaa) =, ', E_corr_contr_s(BL_AOAA) + E_corr_contr_t(BL_AOAA) !, E_corr_contr_ext(BL_AOAA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaao) =, ', E_corr_contr_s(BL_AAAO) + E_corr_contr_t(BL_AAAO) !, E_corr_contr_ext(BL_AAAO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aova) =, ', E_corr_contr_s(BL_AOVA) + E_corr_contr_t(BL_AOVA) !, E_corr_contr_ext(BL_AOVA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaao) =, ', E_corr_contr_s(BL_VAAO) + E_corr_contr_t(BL_VAAO) !, E_corr_contr_ext(BL_VAAO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aovv) =, ', E_corr_contr_s(BL_AOVV) + E_corr_contr_t(BL_AOVV) !, E_corr_contr_ext(BL_AOVV)




            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvao) =, ', E_corr_contr_s(BL_VVAO) + E_corr_contr_t(BL_VVAO) !, E_corr_contr_ext(BL_VVAO)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aava) =, ', E_corr_contr_s(BL_AAVA) + E_corr_contr_t(BL_AAVA) !, E_corr_contr_ext(BL_AAVA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaaa) =, ', E_corr_contr_s(BL_VAAA) + E_corr_contr_t(BL_VAAA) !, E_corr_contr_ext(BL_VAAA)


            write(*,'(A25, 2F25.15)') 'E_corr_contr(aavv) =, ', E_corr_contr_s(BL_AAVV) + E_corr_contr_t(BL_AAVV) !, E_corr_contr_ext(BL_AAVV)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvaa) =, ', E_corr_contr_s(BL_VVAA) + E_corr_contr_t(BL_VVAA) !, E_corr_contr_ext(BL_VVAA)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vavv) =, ', E_corr_contr_s(BL_VAVV) + E_corr_contr_t(BL_VAVV) !, E_corr_contr_ext(BL_VAVV)



            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvva) =, ', E_corr_contr_s(BL_VVVA) + E_corr_contr_t(BL_VVVA) !, E_corr_contr_ext(BL_VVVA)

      end subroutine print_contributions


      subroutine print_contributions_alpha(E_corr_contr, alpha)
            double precision, dimension(:), intent(in) :: E_corr_contr
            double precision, intent(in) :: alpha
            integer :: i
            double precision :: sum_contr

            sum_contr = zero
            do i = 1, 25
                  sum_contr = sum_contr + E_corr_contr(i)
            end do
            
            write(*,'(A40, F25.15)') 'E_corr as sum of contributions', sum_contr
            write(*,'(A25, 2F25.15)') 'E_corr_contr(oooo) =, ', E_corr_contr(BL_OOOO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaaa) =, ', E_corr_contr(BL_AAAA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aoao) =, ', E_corr_contr(BL_AOAO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vava) =, ', E_corr_contr(BL_VAVA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvvv) =, ', E_corr_contr(BL_VVVV), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(ooao) =, ', E_corr_contr(BL_OOAO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aooo) =, ', E_corr_contr(BL_AOOO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(ooaa) =, ', E_corr_contr(BL_OOAA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaoo) =, ', E_corr_contr(BL_AAOO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(oova) =, ', E_corr_contr(BL_OOVA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaoo) =, ', E_corr_contr(BL_VAOO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(oovv) =, ', E_corr_contr(BL_OOVV), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvoo) =, ', E_corr_contr(BL_VVOO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aoaa) =, ', E_corr_contr(BL_AOAA), alpha
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aaao) =, ', E_corr_contr(BL_AAAO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aova) =, ', E_corr_contr(BL_AOVA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaao) =, ', E_corr_contr(BL_VAAO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aovv) =, ', E_corr_contr(BL_AOVV), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvao) =, ', E_corr_contr(BL_VVAO), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aava) =, ', E_corr_contr(BL_AAVA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vaaa) =, ', E_corr_contr(BL_VAAA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(aavv) =, ', E_corr_contr(BL_AAVV), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvaa) =, ', E_corr_contr(BL_VVAA), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vavv) =, ', E_corr_contr(BL_VAVV), alpha 
            write(*,'(A25, 2F25.15)') 'E_corr_contr(vvva) =, ', E_corr_contr(BL_VVVA), alpha


            write(*, '(A3, 26F20.15)') 'zra', alpha, E_corr_contr(BL_OOOO), E_corr_contr(BL_AAAA), E_corr_contr(BL_AOAO), &
                  E_corr_contr(BL_VAVA), E_corr_contr(BL_VVVV), E_corr_contr(BL_OOAO), E_corr_contr(BL_AOOO), &
                  E_corr_contr(BL_OOAA), E_corr_contr(BL_AAOO), E_corr_contr(BL_OOVA), E_corr_contr(BL_VAOO), &
                  E_corr_contr(BL_OOVV), E_corr_contr(BL_VVOO), E_corr_contr(BL_AOAA), E_corr_contr(BL_AAAO), &
                  E_corr_contr(BL_AOVA), E_corr_contr(BL_VAAO), E_corr_contr(BL_AOVV), E_corr_contr(BL_VVAO), &
                  E_corr_contr(BL_AAVA), E_corr_contr(BL_VAAA), E_corr_contr(BL_AAVV), E_corr_contr(BL_VVAA), &
                  E_corr_contr(BL_VAVV), E_corr_contr(BL_VVVA)

      end subroutine print_contributions_alpha



      subroutine pp_rowe(AuxData, HType, ETot, ENuc, MxA, MxS,  n, XOne, TwoNO, IndN, IndX, IndAux, IndMod, NDim, &
            NBasis, NA, NI, NV, NInte1, NInte2, Tpphh, ACAlpha, spin_symm, version, Flags, eorb, onlyA, invertS)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            type(TACppData) :: AuxData
            integer, intent(in) :: HType
            double precision, intent(inout) :: ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:,:), intent(inout) :: MxA, MxS
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            integer, dimension(:), intent(in)    ::IndX
            integer, dimension(:,:), intent(in)  :: IndN
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NDim, NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type(TDA_pphh) :: Tpphh
            integer, intent(in) :: spin_symm, version
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            logical, optional, intent(in) :: onlyA, invertS
            logical :: actual_onlyA, actual_invertS


            double precision, dimension(:), allocatable :: HNO, Ha
            double precision, dimension(:), allocatable :: TwoNOA
            double precision, dimension(:), allocatable :: R00, R11
            integer, dimension(:),          allocatable :: Ind1, Ind2
            integer :: NRDM2, NRDM2Act, NOc
            double precision :: Arspq, Arsqp, Brspq, Crspq, dm, Arssum
            double precision, dimension(:,:), allocatable :: MxT
            double precision, dimension(:, :), intent(out) :: eorb
            ! double precision, dimension(:, :), allocatable :: eorb                                                                                                                                 
            type (tclock) :: timer, timeall
            ! External procedures
            !

            integer, external :: NAddrRDM
            integer, external :: NAddr3
            double precision, external :: FRDM2

            integer :: p, q, r, s, t, u, v, pp, qq
            integer :: i ,j, k, l, ij, a, b, ab, kl
            integer :: ii, pq, rs, pr, qs, ps, qr, kk
            integer, parameter :: u1 = 10, u2=20, u3=30
            integer :: error1, error2, error3
            double precision :: num_f, num_h
            double precision :: eh, el, nu, trace, numf, sum
            integer :: ih, il, rrs, ppq
            integer :: c, d, ac, bd, cd, ad, bc, aab, ccd, nn
            integer :: finito 
            integer :: true_nof_lines, true_nof_lines2
            double precision :: foki, twoelp
            integer, dimension(:, :) , allocatable :: indN3pr, IndN3qs
            integer, dimension(:, :) , allocatable :: indXpr, IndXqs
            integer :: tt, uu, vv

            logical :: do_read
            integer :: twoint_dim
            integer, dimension(:), allocatable :: igfact
            double precision :: temp, temp2
            integer :: jj, ai, bj
            double precision, dimension(:,:), allocatable :: miniA
            integer, dimension(9) :: lst

            logical ::ism, pf, qf


            lst(1) = 3
            lst(2)=5
            lst(3)=12
            lst(4)=13
            lst(5)=18
            lst(6)=21
            lst(7)=24
            lst(8)=26
            lst(9)=29

            call clock_start(timeall)

            NRDM2 = NBasis**2*(NBasis**2+1)/2
            !    NRDM2Act  = (true_NA+2)**2*((true_NA+2)**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NOc = NI + NA

            allocate(IndN3pr(3,NDim*NDim))
            allocate(IndN3qs(3,NDim*NDim))
            allocate(IndXpr(2,NDim*NDim))
            allocate(IndXqs(2,NDim*NDim))


            allocate(HNO(NInte1))
            allocate(Ha(NInte1))
            allocate(Ind1(Nbasis))
            allocate(Ind2(Nbasis))
            allocate (R00(NRDM2Act))
            allocate (R11(NRDM2Act))
            ! allocate(eorb(NBasis, NBasis))

            twoint_dim = size(TwoNO, dim=1)
            allocate(TwoNOA(twoint_dim))
            !    print*, 'allocate(R00(NRDM2Act))', NRDM2Act

            call clock_start(timer)
            Ind1 = 0
            Ind2 = 0

            R00 = zero
            R11 = zero

            if (AuxData%general_version == GENVER_PP)then
                  allocate(MxT(NDim, NDim))
                  MxT = zero
            end if
            MxA = zero
            !

            ! Fill one-electron hamiltonian in NO rep.

            !

            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        HNO(ij) = zero
                        a = i
                        b = j
                        ! do a = 1, NBasis
                        !    do b = 1, Nbasis
                        ab = (max(a, b)*(max(a,b)-1))/2 + min(a, b)
                        !              HNO(ij) = HNO(ij) + URe(i, a)*URe(j,b)*XOne(ab)
                        HNO(ij) = XOne(ab)
                        !             end do
                        !         end do

                  end do
            end do
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        temp = HNO(ij)
                        do r = 1, NBasis
                              temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                        end do
                        eorb(i, j) = temp
                  end do
            end do

            print*, 'Number of inacvive orbitals:', NI
            print*, 'Number of active orbitals:  ', NA
            print*, 'Number of virtual orbitals: ', NV


            Ind1 = 0
            Ind2 = 0

            Ind1 = 0
            Ind2 = 0

            k = 1
            do i = 1, NI+NA+NV
                  if(IndAux(i)==1)then
                        Ind2(i) = k
                        k = k+1
                  end if
                  if(IndAux(i)==2)then
                        Ind2(i) = 0
                  end if
            end do
            call read_2rdm("rdm2.dat", R00, NA)
            call read_2rdm("rdms2.dat", R11, NA)



            ETot = zero
            do i = 1, NBasis
                  ii = (i*(i+1))/2
                  ETot = ETot + two* n(i) * HNO(ii)
            end do

            print*, 'etot1', ETot, NOc

            do p = 1, NI+NA
                  do q = 1, NI+NA
                        do r = 1, NI+NA
                              do s = 1, NI+NA
                                    ETot = ETot + FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis) &
                                          * TwoNO(NAddr3(p, r, q, s))
                                    ! ism = ismixed(lst, p,q,r,s)

                                    ! if (ism)then
                                    !    if (abs(TwoNO(NAddr3(p, r, q, s))).gt.1.d-3)then
                                    !       print*, 'TWONOISMIXED', p, q, r, s, TwoNO(NAddr3(p, r, q, s))
                                    !    end if
                                    ! end if
                              end do
                        end do
                  end do
            end do

            print*, 'Time for pocz: ', clock_readwall(timer)
            call clock_start(timer)
            print*, 'RDCS ETot', ETot+Enuc
            print*, 'hno1', hno(1)
            Ha = zero
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        Ha(ij) = ACAlpha * HNO(ij)

                        if (IndAux(i).eq.IndAux(j))then


                              temp = HNO(ij)
                              
                              do r = 1, NBasis                                    

                                    if (HType==H_DYALL)then     	
                                          if (.not.((IndAux(r).eq.IndAux(i)).and.(IndAux(r).eq.1)))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                          end if
                                    else if (HType==H_GPF)then

                                          if (IndAux(r).ne.IndAux(i))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                          end if
                                    end if

                              end do
                              Ha(ij) = Ha(ij) + (one-ACAlpha)*temp
                              
                        end if
                  end do
            end do

            ! ij = 0
            ! do i = 1, Nbasis
            !       do j = 1, i
            !             ij = ij + 1
            !             if (abs(ha(ij)).gt.1.d-5)then
            !                   print*, 'ha', i, j, Ha(ij)
            !             end if
            !       end do
            ! end do



            
            print*, 'Time for przepisywanie hno: ', clock_readwall(timer)
            print*, 'ha', ha(1)
            call clock_start(timer)
            print*, nbasis
            print*, htype
            print*, indaux
            print*, ndim
            print*, spin_symm

            TwoNOA = TwoNO
            ij = 0
            do i = 1, NBasis
                  do j = 1, i
                        ij = ij + 1
                        kl = 0
                        do k = 1, Nbasis
                              do l = 1, k
                                    kl=kl+1

                                    if (HType==H_DYALL)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l).and.IndAux(i)==1)then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                          end if
                                    else if (HType==H_GPF)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l))then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                                if (abs(TwoNOA(NAddr3(i, j, k, l))).gt.1.d-8)then
                                                end if
                                          end if
                                    else if (HType==H_MID)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l).and.IndAux(i).ne.2)then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                          end if

                                    end if

                              end do
                        end do
                  end do
            end do

            print*, 'Time for przepisywanie calek: ', clock_readwall(timer)

            actual_onlyA = .false.
            if (present(onlyA))then
                  if (onlyA==.true.)then
                        actual_onlyA = .true.
                  end if
            end if


            call clock_start(timer)

            !    print*, 'version', version

            if (version==VER_PP_RPA_TDA_HF)then

                  write(*,'(A41, I10)') 'Constructing A matrix for one determinant', NDim

                  print*, 'pp_rpa_yang_hf_tda'
                  call pp_rpa_yang_hf_tda(MxA, MxS, NDim, NBasis, n, IndN, IndX, IndAux, TwoNOA, Ha, NI, NA, spin_symm, ETot)

            else if (version==VER_PP_RPA_MULTI .or. version==VER_PP_RPATDA_MULTI)then
                  !       print*, 'version multi', spin_symm
                  !
                  !----------------------------------MULTIREFERENCE
                  !       
                  if (spin_symm==2)then
                        i_rowloop: do i = 1, NDim
                              r = IndN(1, i)
                              s = IndN(2, i)
                              rs = IndX(i)
                              j_colloop: do j = 1, NDim
                                    p = IndN(1, j)
                                    q = IndN(2, j)
                                    pq = IndX(j)

                                    Arspq = zero
                                    qs = (max(s, q)*(max(s, q)-1))/2 + min(s, q)
                                    pr = (max(r, p)*(max(r, p)-1))/2 + min(r, p)


                                    Arspq = Arspq + TwoNOA(NAddr3(p,r,q,s))* (one - n(p)-n(q)-n(r)-n(s))
                                    if (p==r)then
                                          Arspq = Arspq + Ha(qs) * (one-n(p)-frac12*n(q)-frac12*n(s))
                                    end if
                                    if (q==s)then
                                          Arspq = Arspq + Ha(pr) * (one-n(q)-frac12*n(p)-frac12*n(r))
                                    end if
                                    if (p==r)then
                                          do t = 1, NBasis                                                                                        
                                                Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(q,s,t,t))-TwoNOA(NAddr3(q,t,t,s)))
                                          end do
                                    end if

                                    if (q==s)then
                                          do t = 1, NBasis
                                                Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(p,r,t,t))-TwoNOA(NAddr3(p,t,t,r)))
                                          end do
                                    end if

                                    Arspq = Arspq + P3a(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq + P3a(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq + P3b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq + P3b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq - P3c(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq - P3c(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq - P4a(q, s, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq - P4a(p, r, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq - P4b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq - P4b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq + P5a(q, s, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq + P5a(p, r, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

                                    Arspq = Arspq + P5b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
                                    Arspq = Arspq + P5b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)


                                    MxA(rs, pq) = Arspq / (1-n(r)-n(s))

                              end do j_colloop
                        end do i_rowloop
                  end if



                  if (spin_symm==0.or.spin_symm==1)then
                        
                        !$omp parallel do collapse(2)&
                        !$omp default(shared) &
                        !$omp private(i, j, tt, uu, vv, t, u, v) &
                        !$omp private(p, q, pq, r, s, rs, qs, pr, qr, ps, Arspq, Arsqp, Arssum, num_f, num_h)
                        j_colloops: do j = 1, NDim
                              i_rowloops: do i = 1, NDim

                                    !                num_f = one/(1-n(r)-n(s))


                                    r = IndN(1, i)
                                    s = IndN(2, i)
                                    rs = IndX(i)

                                    p = IndN(1, j)
                                    q = IndN(2, j)
                                    pq = IndX(j)

                                    Arspq = zero
                                    Arsqp = zero
                                    qs = (max(s, q)*(max(s, q)-1))/2 + min(s, q)
                                    pr = (max(r, p)*(max(r, p)-1))/2 + min(r, p)

                                    qr = (max(r, q)*(max(r, q)-1))/2 + min(r, q)
                                    ps = (max(s, p)*(max(s, p)-1))/2 + min(s, p)

                                    !pierwszy
                                    Arspq = Arspq + TwoNOA(NAddr3(p,r,q,s))* (one - n(p)-n(q)-n(r)-n(s))


                                    if (p==r)then
                                          Arspq = Arspq + Ha(qs) * (one-n(p)-frac12*n(q)-frac12*n(s))
                                    end if


                                    if (q==s)then
                                          Arspq = Arspq + Ha(pr) * (one-n(q)-frac12*n(p)-frac12*n(r))
                                    end if

                                    if (p==r)then
                                          do t = 1, NI+NA
                                                Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(q,s,t,t))-TwoNOA(NAddr3(q,t,t,s)))
                                          end do
                                    end if

                                    if (q==s)then
                                          temp = zero
                                          do t = 1, NI+NA
                                                Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(p,r,t,t))-TwoNOA(NAddr3(p,t,t,r)))
                                          end do
                                    end if


                                    ! !                 !P3a-1
                                    if (IndAux(s)==(1).and.IndAux(p)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq + TwoNOA(NAddr3(q,u,r,t))* &
                                                            Gam(t, s, p, u, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do
                                    else
                                          if ((IndAux(s)==(1).and.IndAux(p)==(0)).or.&
                                                (IndAux(s)==(0).and.IndAux(p)==(1)).or.&
                                                (IndAux(s)==(0).and.IndAux(p)==(0)))then
                                                Arspq = Arspq + n(p)*n(s)*TwoNOA(NAddr3(q,s,r,p))
                                          end if
                                    end if

                                    !                 !P3a-2
                                    if (IndAux(r)==(1).and.IndAux(q)==(1))then
                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq + TwoNOA(NAddr3(p,u,s,t))* &
                                                            Gam(t, r, q, u, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do
                                    else
                                          if ((IndAux(r)==(1).and.IndAux(q)==(0)).or.&
                                                (IndAux(r)==(0).and.IndAux(q)==(1)).or.&
                                                (IndAux(r)==(0).and.IndAux(q)==(0)))then
                                                Arspq = Arspq + n(q)*n(r)*TwoNOA(NAddr3(p,r,s,q))
                                          end if
                                    end if


                                    !P3b-1
                                    if (IndAux(s)==(1).and.IndAux(q)==(1))then
                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq + TwoNOA(NAddr3(p,u,r,t))* &
                                                            Gam(t, s, u, q, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do

                                          if (s==q)then  
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + n(t)*n(q)*TwoNOA(NAddr3(p,t,r,t))
                                                end do
                                          end if
                                    else
                                          if (s==q)then
                                                if (IndAux(s).ne.1)then 
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arspq = Arspq + TwoNOA(NAddr3(p,t,r,t))* n(t)*n(q)
                                                      end do
                                                end if
                                          end if
                                    end if


                                    !P3b-2
                                    if (IndAux(r)==(1).and.IndAux(p)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq + TwoNOA(Naddr3(q,u,s,t))* &
                                                            Gam(t, r, u, p, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do
                                          if (r==p)then
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + n(t)*n(p)*TwoNOA(Naddr3(q,t,s,t))
                                                end do
                                          end if

                                    else
                                          if (r==p)then
                                                if (IndAux(r).ne.1)then
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arspq = Arspq + TwoNOA(Naddr3(q,t,s,t))* n(t)*n(p)
                                                      end do
                                                end if
                                          end if
                                    end if

                                    ! !                  ! !P3c-1,2

                                    
                                    if (IndAux(s)==(1).and.IndAux(q)==(1))then
                                          temp = zero
                                          do tt = NI+1, NI+NA  !!!!!!!!!!!!! tego jeszcze nie ma w fofo
                                                do uu = NI+1, NI+NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq -  TwoNOA(NAddr3(p,r,t,u))* &
                                                            (Gam(t, s, u, q, R00, R11, n, Ind2, NA, 0)+Gam(t, s, u, q, R00, R11, n, Ind2, NA, 1))
                                                      temp = temp  -  TwoNOA(NAddr3(p,r,t,u))* &
                                                            (Gam(t, s, u, q, R00, R11, n, Ind2, NA, 0)+Gam(t, s, u, q, R00, R11, n, Ind2, NA, 1))
                                                end do
                                          end do

                                          if (s==q)then
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - two* n(t)*n(s) * TwoNOA(NAddr3(p,r,t,t))
                                                end do
                                          end if
                                    else
                                          if (s==q)then
                                                if (IndAux(s)==0)then
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arspq = Arspq - two* n(t)*n(s) * TwoNOA(NAddr3(p,r,t,t))
                                                      end do
                                                end if
                                          end if
                                    end if


                                    if ((IndAux(s)==0.and.IndAux(q)==1).or.&
                                          (IndAux(s)==1.and.IndAux(q)==0).or.&
                                          (IndAux(s)==0.and.IndAux(q)==0))then
                                          Arspq = Arspq + n(q)*n(s) * TwoNOA(NAddr3(p,r,q,s))
                                    end if

                                    !P3c-3,4
                                    if (IndAux(r)==(1).and.IndAux(p)==(1))then
                                          temp = zero
                                          do tt = NI+1, NI+NA
                                                do uu = NI+1, NI+NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arspq = Arspq -  TwoNOA(Naddr3(q,s,t,u))* &
                                                            (Gam(t, r, u, p, R00, R11, n, Ind2, NA, 0)+Gam(t, r, u, p, R00, R11, n, Ind2, NA, 1))
                                                      temp = temp -  TwoNOA(Naddr3(q,s,t,u))* &
                                                            (Gam(t, r, u, p, R00, R11, n, Ind2, NA, 0)+Gam(t, r, u, p, R00, R11, n, Ind2, NA, 1))
                                                end do
                                          end do

                                          if (r==p)then
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - two* n(t)*n(r) * TwoNOA(Naddr3(q,s,t,t))
                                                end do
                                          end if
                                    else
                                          if (r==p)then
                                                if (IndAux(r)==0)then
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arspq = Arspq - two* n(t)*n(r) * TwoNOA(Naddr3(q,s,t,t))
                                                      end do
                                                end if
                                          end if
                                    end if

                                    if ((IndAux(r)==0.and.IndAux(p)==1).or.&
                                          (IndAux(r)==1.and.IndAux(p)==0).or.&
                                          (IndAux(r)==0.and.IndAux(p)==0))then
                                          Arspq = Arspq + n(p)*n(r) * TwoNOA(Naddr3(q,s,p,r))

                                    end if

                                    !                 !P4a-2
                                    if (p==r)then
                                          temp = 0
                                          if (IndAux(q)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq - frac12  * (TwoNOA(NAddr3(s, t, u, v)) * &
                                                                        Gam(t, u, q, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(q) * n(t) * TwoNOA(NAddr3(s, q, t, t))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(q) * n(t) * TwoNOA(NAddr3(s, q, t, t))

                                                end do
                                          end if
                                    end if


                                    !P4a-1
                                    if (q==s)then

                                          if (IndAux(p)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq - frac12  * (TwoNOA(Naddr3(r, t, u, v)) * &
                                                                        Gam(t, u, p, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(p) * n(t) * TwoNOA(Naddr3(r, p, t, t))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(p) * n(t) * TwoNOA(Naddr3(r, p, t, t))
                                                end do

                                          end if
                                    end if

                                    ! !                 !P4b-1
                                    if (p==r)then

                                          if (IndAux(s)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI+1, NI+NA
                                                            do vv = NI+1, NI+NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq - frac12 * (TwoNOA(NAddr3(q, v, t, u)) * &
                                                                        Gam(t, s, u, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt) 
                                                      Arspq = Arspq - frac12 * n(s) * n(t) * TwoNOA(NAddr3(q, s, t, t))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(s) * n(t) * TwoNOA(NAddr3(q, s, t, t))
                                                end do
                                          end if
                                    end if

                                    !P4b-2
                                    if (q==s)then

                                          if (IndAux(r)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI+1, NI+NA
                                                            do vv = NI+1, NI+NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq - frac12 * (TwoNOA(Naddr3(p, v, t, u)) * &
                                                                        Gam(t, r, u, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(r) * n(t) * TwoNOA(Naddr3(p, r, t, t))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq - frac12 * n(r) * n(t) * TwoNOA(Naddr3(p, r, t, t))
                                                end do
                                                if (pq==3.and.rs==3)then
                                                      print*, '32', Arspq
                                                end if
                                          end if
                                    end if
                                    ! !                 !P5a
                                    if (p==r)then
                                          temp = zero
                                          if (IndAux(q)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq + frac14 * ((TwoNOA(NAddr3(s, u, t, v)) - TwoNOA(NAddr3(s, t, u, v)))&
                                                                        * Gam(t, u, q, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(t) * n(q) * (TwoNOA(NAddr3(s, t, q, t)) -TwoNOA(NAddr3(s, q, t, t)))
                                                end do
                                          else
                                                temp=zero
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(t) * n(q) * (TwoNOA(NAddr3(s, t, q, t)) -TwoNOA(NAddr3(s, q, t, t)))
                                                end do
                                          end if
                                    end if
                                    !P5a
                                    if (q==s)then

                                          if (IndAux(p)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq + frac14 * ((TwoNOA(Naddr3(r, u, t, v)) - TwoNOA(Naddr3(r, t, u, v)))&
                                                                        * Gam(t, u, p, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(t) * n(p) * (TwoNOA(Naddr3(r, t, p, t)) -TwoNOA(Naddr3(r, p, t, t)))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(t) * n(p) * (TwoNOA(Naddr3(r, t, p, t)) -TwoNOA(Naddr3(r, p, t, t)))
                                                end do
                                          end if
                                    end if


                                    ! ! !                 !P5b
                                    if (p==r)then

                                          if (IndAux(s)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq + frac14 * ((TwoNOA(NAddr3(q, u, t, v)) - TwoNOA(NAddr3(q, v, t, u)))&
                                                                        * Gam(t, s, u, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(s) * n(t) * (TwoNOA(NAddr3(q, t, t, s)) -TwoNOA(NAddr3(q, s, t, t)))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(s) * n(t) * (TwoNOA(NAddr3(q, t, t, s)) -TwoNOA(NAddr3(q, s, t, t)))
                                                end do
                                          end if
                                    end if

                                    !P5b
                                    if (q==s)then

                                          if (IndAux(r)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arspq = Arspq + frac14 * ((TwoNOA(Naddr3(p, u, t, v)) - TwoNOA(Naddr3(p, v, t, u)))&
                                                                        * Gam(t, r, u, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(r) * n(t) * (TwoNOA(Naddr3(p, t, t, r)) -TwoNOA(Naddr3(p, r, t, t)))
                                                end do
                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arspq = Arspq + frac12 * n(r) * n(t) * (TwoNOA(Naddr3(p, t, t, r)) -TwoNOA(Naddr3(p, r, t, t)))
                                                end do
                                          end if
                                    end if




                                    ! ! ! ! ! Arsqp
                                    ! !drugi
                                    Arsqp = Arsqp - TwoNOA(NAddr3(q,r,p,s))* (one - n(q)-n(p)-n(r)-n(s))

                                    if (q==r)then
                                          Arsqp = Arsqp - Ha(ps) * (one-n(q)-frac12*n(p)-frac12*n(s))
                                          !                   print*, 'val3', Arsqp
                                    end if

                                    if (p==s)then
                                          Arsqp = Arsqp - Ha(qr) * (one-n(p)-frac12*n(q)-frac12*n(r))
                                          !                  print*, 'val4', Arsqp
                                    end if

                                    if (q==r)then
                                          do t = 1, NBasis
                                                Arsqp = Arsqp - n(t) * (two*TwoNOA(NAddr3(p,s,t,t))-TwoNOA(NAddr3(p,t,t,s)))
                                          end do
                                    end if

                                    if (p==s)then
                                          do t = 1, NBasis
                                                Arsqp = Arsqp - n(t) * (two*TwoNOA(NAddr3(q,r,t,t))-TwoNOA(NAddr3(q,t,t,r)))

                                          end do
                                    end if

                                    !p3a-3

                                    if (IndAux(s)==(1).and.IndAux(q)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp - TwoNOA(NAddr3(p,u,r,t))* &
                                                            Gam(t, s, q, u, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do
                                    else
                                          if ((IndAux(s)==(1).and.IndAux(q)==(0)).or.&
                                                (IndAux(s)==(0).and.IndAux(q)==(1)).or.&
                                                (IndAux(s)==(0).and.IndAux(q)==(0)))then
                                                Arsqp = Arsqp - n(q)*n(s)*TwoNOA(NAddr3(p,s,r,q))
                                          end if
                                    end if

                                    !                 !P3b-3
                                    if (IndAux(s)==(1).and.IndAux(p)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp - TwoNOA(NAddr3(q,u,r,t))* &
                                                            Gam(t, s, u, p, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do

                                          if (s==p)then  
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - n(t)*n(p)*TwoNOA(NAddr3(q,t,r,t))
                                                end do
                                          end if

                                    else
                                          if (s==p)then
                                                if (IndAux(s).ne.1)then 
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arsqp = Arsqp - TwoNOA(NAddr3(q,t,r,t))* n(t)*n(p)
                                                      end do
                                                end if
                                          end if
                                    end if

                                    !p3c-5,6
                                    if (IndAux(s)==(1).and.IndAux(p)==(1))then
                                          do tt = NI+1, NI+NA
                                                do uu = NI+1, NI+NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp +  TwoNOA(NAddr3(q,r,t,u))* &
                                                            (Gam(t, s, u, p, R00, R11, n, Ind2, NA, 0)+Gam(t, s, u, p, R00, R11, n, Ind2, NA, 1))
                                                end do
                                          end do



                                          if (s==p)then
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + two* n(t)*n(s) * TwoNOA(NAddr3(q,r,t,t))
                                                end do
                                          end if
                                    else
                                          if (s==p)then
                                                if (IndAux(s)==0)then
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arsqp = Arsqp + two* n(t)*n(s) * TwoNOA(NAddr3(q,r,t,t))
                                                      end do
                                                end if
                                          end if
                                    end if

                                    if ((IndAux(s)==0.and.IndAux(p)==1).or.&
                                          (IndAux(s)==1.and.IndAux(p)==0).or.&
                                          (IndAux(s)==0.and.IndAux(p)==0))then
                                          Arsqp = Arsqp - n(p)*n(s) * TwoNOA(NAddr3(q,r,p,s))
                                    end if

                                    !p4a
                                    if (q==r)then

                                          if (IndAux(p)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp + frac12  * (TwoNOA(NAddr3(s, t, u, v)) * &
                                                                        Gam(t, u, p, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(p) * n(t) * TwoNOA(NAddr3(s, p, t, t))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(p) * n(t) * TwoNOA(NAddr3(s, p, t, t))
                                                end do

                                          end if
                                    end if

                                    !P4b
                                    if (q==r)then

                                          if (IndAux(s)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI+1, NI+NA
                                                            do vv = NI+1, NI+NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp + frac12 * (TwoNOA(NAddr3(p, v, t, u)) * &
                                                                        Gam(t, s, u, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt) 
                                                      Arsqp = Arsqp + frac12 * n(s) * n(t) * TwoNOA(NAddr3(p, s, t, t))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(s) * n(t) * TwoNOA(NAddr3(p, s, t, t))
                                                end do
                                          end if
                                    end if

                                    !P5a
                                    if (q==r)then

                                          if (IndAux(p)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp - frac14 * ((TwoNOA(NAddr3(s, u, t, v)) - TwoNOA(NAddr3(s, t, u, v)))&
                                                                        * Gam(t, u, p, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(t) * n(p) * (TwoNOA(NAddr3(s, t, p, t)) -TwoNOA(NAddr3(s, p, t, t)))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(t) * n(p) * (TwoNOA(NAddr3(s, t, p, t)) -TwoNOA(NAddr3(s, p, t, t)))
                                                end do
                                          end if
                                    end if

                                    !P5b
                                    if (q==r)then

                                          if (IndAux(s)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp - frac14 * ((TwoNOA(NAddr3(p, u, t, v)) - TwoNOA(NAddr3(p, v, t, u)))&
                                                                        * Gam(t, s, u, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(s) * n(t) * (TwoNOA(NAddr3(p, t, t, s)) -TwoNOA(NAddr3(p, s, t, t)))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(s) * n(t) * (TwoNOA(NAddr3(p, t, t, s)) -TwoNOA(NAddr3(p, s, t, t)))
                                                end do
                                          end if
                                    end if


                                    !                 !P3a-4
                                    if (IndAux(r)==(1).and.IndAux(p)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp - TwoNOA(Naddr3(q,u,s,t))* &
                                                            Gam(t, r, p, u, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do
                                    else
                                          if ((IndAux(r)==(1).and.IndAux(p)==(0)).or.&
                                                (IndAux(r)==(0).and.IndAux(p)==(1)).or.&
                                                (IndAux(r)==(0).and.IndAux(p)==(0)))then
                                                Arsqp = Arsqp - n(p)*n(r)*TwoNOA(Naddr3(q,r,s,p))
                                          end if
                                    end if

                                    !                 !p3b-4
                                    if (IndAux(r)==(1).and.IndAux(q)==(1))then

                                          do tt = NI+1, NI+NA
                                                do uu = NI + 1, NI + NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp - TwoNOA(Naddr3(p,u,s,t))* &
                                                            Gam(t, r, u, q, R00, R11, n, Ind2, NA, 1)
                                                end do
                                          end do

                                          if (r==q)then  
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - n(t)*n(q)*TwoNOA(Naddr3(p,t,s,t))
                                                end do
                                          end if

                                    else
                                          if (r==q)then
                                                if (IndAux(r).ne.1)then 
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arsqp = Arsqp -TwoNOA(Naddr3(p,t,s,t))* n(t)*n(q)
                                                      end do
                                                end if
                                          end if
                                    end if

                                    !p3c-7,8
                                    if (IndAux(r)==(1).and.IndAux(q)==(1))then
                                          do tt = NI+1, NI+NA
                                                do uu = NI+1, NI+NA
                                                      t = IndMod(tt)
                                                      u = IndMod(uu)
                                                      Arsqp = Arsqp +  TwoNOA(Naddr3(p,s,t,u))* &
                                                            (Gam(t, r, u, q, R00, R11, n, Ind2, NA, 0)+Gam(t, r, u, q, R00, R11, n, Ind2, NA, 1))
                                                end do
                                          end do

                                          if (r==q)then
                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + two* n(t)*n(r) * TwoNOA(Naddr3(p,s,t,t))
                                                end do
                                          end if
                                    else
                                          if (r==q)then
                                                if (IndAux(r)==0)then
                                                      do tt = 1, NI+NA
                                                            t = IndMod(tt)
                                                            Arsqp = Arsqp + two* n(t)*n(r) * TwoNOA(Naddr3(p,s,t,t))
                                                      end do
                                                end if
                                          end if
                                    end if

                                    if ((IndAux(r)==0.and.IndAux(q)==1).or.&
                                          (IndAux(r)==1.and.IndAux(q)==0).or.&
                                          (IndAux(r)==0.and.IndAux(q)==0))then
                                          Arsqp = Arsqp - n(q)*n(r) * TwoNOA(Naddr3(p,s,q,r))
                                    end if

                                    !p4a
                                    if (p==s)then

                                          if (IndAux(q)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp + frac12  * (TwoNOA(Naddr3(r, t, u, v)) * &
                                                                        Gam(t, u, q, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(q) * n(t) * TwoNOA(Naddr3(r, q, t, t))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(q) * n(t) * TwoNOA(Naddr3(r, q, t, t))
                                                end do

                                          end if
                                    end if
                                    !p4b
                                    if (p==s)then

                                          if (IndAux(r)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI+1, NI+NA
                                                            do vv = NI+1, NI+NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp + frac12 * (TwoNOA(Naddr3(q, v, t, u)) * &
                                                                        Gam(t, r, u, v, R00, R11, n, Ind2, NA, 1))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt) 
                                                      Arsqp = Arsqp + frac12 * n(r) * n(t) * TwoNOA(Naddr3(q, r, t, t))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp + frac12 * n(r) * n(t) * TwoNOA(Naddr3(q, r, t, t))
                                                end do
                                          end if
                                    end if
                                    !p5a
                                    if (p==s)then

                                          if (IndAux(q)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp - frac14 * ((TwoNOA(Naddr3(r, u, t, v)) - TwoNOA(Naddr3(r, t, u, v)))&
                                                                        * Gam(t, u, q, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(t) * n(q) * (TwoNOA(Naddr3(r, t, q, t)) -TwoNOA(Naddr3(r, q, t, t)))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(t) * n(q) * (TwoNOA(Naddr3(r, t, q, t)) -TwoNOA(Naddr3(r, q, t, t)))
                                                end do
                                          end if
                                    end if
                                    !p5b
                                    if (p==s)then

                                          if (IndAux(r)==(1))then
                                                do tt = NI+1, NI+NA
                                                      do uu = NI + 1, NI + NA
                                                            do vv = NI + 1, NI + NA
                                                                  t = IndMod(tt)
                                                                  u = IndMod(uu)
                                                                  v = IndMod(vv)
                                                                  Arsqp = Arsqp - frac14 * ((TwoNOA(Naddr3(q, u, t, v)) - TwoNOA(Naddr3(q, v, t, u)))&
                                                                        * Gam(t, r, u, v, R00, R11, n, Ind2, NA, 0))
                                                            end do
                                                      end do
                                                end do

                                                do tt = 1, NI
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp -frac12 * n(r) * n(t) * (TwoNOA(Naddr3(q, t, t, r)) -TwoNOA(Naddr3(q, r, t, t)))
                                                end do

                                          else
                                                do tt = 1, NI+NA
                                                      t = IndMod(tt)
                                                      Arsqp = Arsqp - frac12 * n(r) * n(t) * (TwoNOA(Naddr3(q, t, t, r)) -TwoNOA(Naddr3(q, r, t, t)))
                                                end do
                                          end if
                                    end if


                                    
                                    num_f = one
                                    !num_f = one/(1-n(r)-n(s))
                                    
                                    select case(spin_symm)
                                    case(0)
                                          num_h = one
                                          if (p==q.and.r==s)then
                                                num_h = frac12
                                          else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
                                                num_h = sqrt(frac12)
                                          end if

                                          Arssum = Arspq  - Arsqp
                                          !                   MxA(rs, pq) = num_f *num_h * Arssum
                                          MxA(rs, pq) = num_f * num_h * Arssum

                                          ! if	(abs(MxA(rs, pq)).gt.1.d-5)then
                                          !       !                      print*, rs, pq, '  |  ', p, q, r, s, MxA(rs, pq)
                                          !       write(*,'(2I5, A5, 4I5, 2F30.16)') rs, pq,	'  |  ', p, q, r, s, MxA(rs, pq)
                                          ! end if


                                          !                    if (IndAux(p)==2.and.q==6.and.IndAux(r)==2.and.s==3)then
                                          ! !                   if (IndAux(p)==2.and.IndAux(q)==1.and.IndAux(r)==2.and.IndAux(s)==1)then
                                          !                       if (abs(MxA(rs, pq)).gt.1.d-6)then
                                          !                          print*, 'wieksz', p, q, r, s, MxA(rs, pq)
                                          !                       end if
                                          !                    end if

                                          ! if (abs(MxA(rs, pq)).gt.1.d+3)then
                                          !    print*, p, q, r, s, MxA(rs, pq)
                                          !    stop
                                          ! end if

                                    case(1)
                                          if (p.ne.q.and.r.ne.s)then
                                                Arssum = Arspq + Arsqp
                                                MxA(rs, pq) = num_f * Arssum
                                                ! if   (abs(MxA(rs, pq)).gt.1.d-5)then
                                                !    write(*,'(2I5, A5, 4I5, F30.16)') rs, pq, '  |  ', p, q, r, s, MxA(rs, pq)
                                                ! end if

                                          end if
                                    end select

                              end do i_rowloops
                        end do j_colloops
                        !          !$omp end parallel do
                  end if


                  ! print*, 'lllaaa', NDim
                  !        do ai = 1, size(MxA, dim=1)
                  !              do bj = 1, size(MxA, dim=1)
                  !                 if (abs(MxA(ai, bj)).gt.1.d+1)then
                  !                    print*, 'maaxf0', ai, bj, MxA(ai, bj)
                  ! 		end if
                  !              end do
                  !           end do
                  !           print*, 'sraa'




                  print*, ''
                  print*, 'Time for construction of A: ', clock_readwall(timer)
                  print*, ''
                  call clock_start(timer)


                  ! call clock_start(timer)       
                  ! call sss(MxS, NDim, NBasis, n, IndN, IndX, IndAux, TwoNOA, Ha, NI, NA, spin_symm, Tpphh, Etot, version)
                  ! print*, 'Time na sss: ', clock_readwall(timer)

                  !  print*, 'lllaasss', NDim
                  ! do ai = 1, size(MxS, dim=1)
                  !       do bj = 1, size(MxS, dim=1)
                  !          if (abs(MxS(ai, bj)).gt.1.d-5)then
                  !             print*, 'mxs', ai, bj, MxS(ai, bj)
                  !          end if
                  !       end do
                  !    end do
                  !    print*, 'sraa'
                  !    stop

                  ! actual_onlyA = .false.
                  ! if (present(onlyA))then
                  !    if (onlyA==.true.)then
                  !       actual_onlyA = .true.
                  !    end if
                  ! end if

                  ! if (actual_onlyA==.false.)then
                  !    !          print*, 'Czas sss: ', clock_readwall(timer)
                  !    call clock_start(timer)
                  !    MxT = zero
                  !    call clock_start(timer)       
                  !    call dgemm("N", "N", NDim, NDim, NDim, 1.d+0, MxS, NDim, MxA, NDim, zero, MxT, NDim) !< to tak samo

                  !    MxA = MxT
                  !    print*, 'Time na mxa: ', clock_readwall(timer)
                  ! !    print*, 'lllaaaqqq', NDim
                  ! ! do ai = 1, size(MxA, dim=1)
                  ! !       do bj = 1, size(MxA, dim=1)
                  ! !          if (abs(MxA(ai, bj)).gt.1.d+1)then
                  ! !             print*, 'maaxf0', ai, bj, MxA(ai, bj)
                  ! !  	end if
                  ! !       end do
                  ! !    end do
                  ! !    print*, 'sraa'

                  ! else
                  ! end if
                  !       print*, 'Time na reszte: ', clock_readwall(timer)

            end if


            print*, 'Czas na macierz: '//str(clock_readwall(timer),d=2)
            call clock_start(timer)       
            call sss(MxS, NDim, NBasis, n, IndN, IndX, IndAux, TwoNOA, Ha, NI, NA, spin_symm, Tpphh, Etot, version) 
            
            print*, 'Czas sss: '//str(clock_readwall(timer),d=2)
            call clock_start(timer)       
            call dgemm("N", "N", NDim, NDim, NDim, 1.d+0, MxS, NDim, MxA, NDim, zero, MxT, NDim) !< to tak samo
            
            print*, 'Czas na dgemm: '//str(clock_readwall(timer),d=2)
            MxA = MxT


            if (version ==VER_PP_RPATDA_MULTI)then
                  i_rowloopstda: do i = 1, NDim
                        r = IndN(1, i)
                        s = IndN(2, i)
                        rs = IndX(i)
                        j_colloopstda: do j = 1, NDim
                              p = IndN(1, j)
                              q = IndN(2, j)
                              pq = IndX(j)
                              if (not (IndAux(r)==2.and.IndAux(s)==2.and.IndAux(p)==2.and.IndAux(q)==2))then
                                    MxA(rs, pq) = zero
                              end if

                        end do j_colloopstda
                  end do i_rowloopstda
            end if

            do i = 1, NDim
                  MxS(i, i) = one / MxS(i, i)
            end do


            ! do i = 1, NDim
            !       p = IndN(1, i)
            !       q = IndN(2, i)
            !       MxS(i,i) = (1-n(p)-n(q)) !(to macierz S^-1                                                                                                                                                                                   
            ! end do

            print*, 'Time na reszte: ', clock_readwall(timer) 
            ! call clock_start(timer)


            ! actual_invertS = .false.
            ! if (present(invertS))then
            !    if (invertS == .true.) then
            !       actual_invertS = .true.
            !    end if
            ! end if

            ! if (actual_invertS==.false.)then
            !    ! the MxS matrix returned from sss procedure is actually S^(-1), therefore we invert it. The matrix is diagonal.
            !    do i = 1, NDim
            !       MxS(i, i) = one / MxS(i, i)
            !    end do

            ! else

            !    end if
            print*, 'Time na all: ', clock_readwall(timeall)


      end subroutine pp_rowe



      subroutine pp_rpa_yang_hf_tda(MxA, MxS, NDim, NBasis, n, IndN, IndX, IndAux, TwoNO, HNO, NI, NA, spin_symm, ETot)
            use math_constants
            double precision, dimension(:,:), intent(inout) :: MxA, MxS
            double precision, dimension(:),      intent(in) :: TwoNO ,HNO
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in)    ::IndX
            integer, dimension(:,:), intent(in)  :: IndN
            integer, dimension(:), intent(in)    :: IndAux
            integer, intent(in) :: NDim, NBasis
            integer, intent(in) :: NI, NA
            integer, intent(in) :: spin_symm
            double precision, intent(in) :: ETot
            integer :: c, d, cd, a, b, ab, aab, ccd, i, j, k, l
            integer :: ih, il, ac, bd, ad, bc, nn
            double precision :: eh, el, nu, num_f, num2_f
            double precision :: Arspq, Brspq, Crspq
            double precision :: e_trace_A, e_trace_C
            integer, external :: NAddr3

            ! print*, NI
            ! ih = (max(NI, NI)*(max(NI, NI)-1))/2 + min(NI, NI)
            ! print*, 'ihihih', ih

            ! eh = eorb(NI, NI, NI, NI, TWONO, HNO(ih), NI, NA)

            ! il = (max(NI+1, NI+1)*(max(NI+1, NI+1)-1))/2 + min(NI+1, NI+1)
            ! el = eorb(NI+1, NI+1, NI+1, NI+1, TWONO, HNO(il), NI, NA)
            ! nu = (eh+el)/two

            ! print*, 'E(HOMO):', eh
            ! print*, 'E(LUMO):', el
            ! print*, 'nu', nu, sqrt(frac12)
            e_trace_A = zero
            e_trace_C = zero

            MxS = zero
            select case(spin_symm)
            case(1)
                  ! triplet                                                                                                                                                      
                  i_rowloopt: do i = 1, NDim
                        c = IndN(1, i)
                        d = IndN(2, i)
                        cd = IndX(i)
                        j_colloopt: do j = 1, NDim
                              a = IndN(1, j)
                              b = IndN(2, j)
                              ab = IndX(j)

                              if (a > b.and.c > d)then
                                    if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                          IndAux(a)==(2).and.IndAux(b)==(2))then
                                          Arspq = zero
                                          ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                          bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)

                                          Arspq = TwoNO(NAddr3(a,c,b,d))-TwoNO(NAddr3(a,d,b,c))

                                          if (b==d)then
                                                Arspq = Arspq + fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                          end if
                                          if (a==c)then
                                                Arspq = Arspq + fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                          end if


                                          ! if (a==c)then
                                          !    if (b==d)then
                                          !       Arspq = Arspq + eorbs(a, TWONO, HNO(ac), NI, NA)
                                          !       Arspq = Arspq + eorbs(b, TWONO, HNO(bd), NI, NA)
                                          !       !Arspq = Arspq -two*nu
                                          !    end if
                                          ! end if

                                          MxA(cd, ab) = Arspq

                                          if (cd == ab) then
                                                e_trace_A = e_trace_A + Arspq
                                                MxS(cd, ab) = one
                                          end if

                                    end if

                                    if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                          IndAux(a)==(0).and.IndAux(b)==(0))then

                                          Crspq = zero
                                          ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                          bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)

                                          Crspq =  TwoNO(NAddr3(a,c,b,d))-TwoNO(NAddr3(a,d,b,c))

                                          if (b==d)then
                                                Crspq = Crspq - fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                          end if
                                          if (a==c)then
                                                Crspq = Crspq - fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                          end if


                                          ! if (a==c)then
                                          !    if (b==d)then
                                          !       Crspq = Crspq - eorbs(a, TWONO, HNO(ac), NI, NA)
                                          !       Crspq = Crspq - eorbs(b, TWONO, HNO(bd), NI, NA)
                                          !       !Crspq = Crspq + two*nu
                                          !    end if
                                          ! end if

                                          MxA(cd, ab) = -Crspq !!!!DUE TO TDA APPROX

                                          if (cd == ab) then
                                                e_trace_C = e_trace_C + Crspq
                                                MxS(cd, ab) = -one
                                          end if


                                    end if

                                    if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                          IndAux(a)==(2).and.IndAux(b)==(2))then

                                          Brspq = TwoNO(NAddr3(a,c,b,d))-TwoNO(NAddr3(a,d,b,c))

                                          MxA(cd, ab) =  Brspq !!!! DUE TO TDA APPROX
                                    end if
                                    if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                          IndAux(a)==(0).and.IndAux(b)==(0))then

                                          Brspq = TwoNO(NAddr3(a,c,b,d))-TwoNO(NAddr3(a,d,b,c))
                                          MxA(cd, ab) =  -Brspq !!!! DUE TO TDA APPROX
                                    end if

                              end if
                        end do j_colloopt
                  end do i_rowloopt

            case(0)
                  print*, 'tak, robie tutaj'
                  ! singlet
                  e_trace_A = zero
                  e_trace_C = zero

                  i_rowloops: do i = 1, NDim
                        c = IndN(1, i)
                        d = IndN(2, i)
                        cd = IndX(i)
                        j_colloops: do j = 1,NDim
                              a = IndN(1, j)
                              b = IndN(2, j)
                              ab = IndX(j)


                              if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                    IndAux(a)==(2).and.IndAux(b)==(2))then

                                    Arspq = zero
                                    ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                    bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)

                                    ad = (max(a, d)*(max(a, d)-1))/2 + min(a, d)
                                    bc = (max(b, c)*(max(b, c)-1))/2 + min(b, c)

                                    num_f = one
                                    if (a==b)num_f = num_f * sqrt(frac12)
                                    if (c==d)num_f = num_f * sqrt(frac12)

                                    Arspq = TwoNO(NAddr3(a,c,b,d)) +  TwoNO(NAddr3(a,d,b,c))
                                    ! print*, 'qw1', Arspq
                                    ! wersja ogolna z fokianami

                                    if (b==d)then
                                          Arspq = Arspq + fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                    end if
                                    if (a==c)then                   
                                          Arspq = Arspq + fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                    end if


                                    if (b==c)then
                                          Arspq = Arspq + fokian(a, d, TWONO, HNO(ad), n, NI, NA, NBasis)
                                    end if

                                    if (a==d)then
                                          Arspq = Arspq + fokian(b, c, TWONO, HNO(bc), n, NI, NA, NBasis)
                                    end if



                                    !                 if (a==c)then
                                    !                    if (b==d)then
                                    !                       Arspq = Arspq + eorbs(a, TWONO, HNO(ac), NI, NA)
                                    !                       Arspq = Arspq + eorbs(b, TWONO, HNO(bd), NI, NA)
                                    !                       !    Arspq = Arspq -two*nu
                                    !                    end if
                                    !                 end if
                                    !                 if (a==d)then
                                    !                    if (b==c)then
                                    !                       Arspq = Arspq + eorbs(a, TWONO, HNO(ad), NI, NA)
                                    !                       Arspq = Arspq + eorbs(b, TWONO, HNO(bc), NI, NA)
                                    ! !                      Arspq = Arspq -two*nu
                                    !                    end if
                                    !                 end if

                                    MxA(cd, ab) =  num_f*Arspq

                                    if (cd == ab) then
                                          e_trace_A = e_trace_A + num_f*Arspq
                                          MxS(cd, ab) = one
                                    end if


                                    ! if (cd == ab) then
                                    !    print*, 'dodaje', ETot
                                    !    MxA(cd, ab) = MxA(cd, ab) + ETot
                                    ! end if
                                    !                if (abs(MxA(cd, ab)).gt.1.d-4)then
                                    ! print*, ''
                                    ! write(*,'(A8, 2I4, 3F15.8)')	'MxA A', cd, ab, num_f ,  Arspq, MxA(cd, ab)
                                    ! print*,	''
                                    !            end if
                                    !                print*, 'grl', a, b, c, d, cd, ab
                              end if

                              if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                    IndAux(a)==(0).and.IndAux(b)==(0))then

                                    Crspq = zero
                                    ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                    bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)
                                    ad = (max(a, d)*(max(a, d)-1))/2 + min(a, d)
                                    bc = (max(b, c)*(max(b, c)-1))/2 + min(b, c)

                                    num_f = one
                                    if (a==b)num_f = num_f * sqrt(frac12)
                                    if (c==d)num_f = num_f * sqrt(frac12)

                                    Crspq =  TwoNO(NAddr3(a,c,b,d))+TwoNO(NAddr3(a,d,b,c))

                                    ! wersja ogolna z fokianami                                                                                                                          
                                    if (b==d)then
                                          Crspq = Crspq - fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                    end if

                                    if (a==c)then
                                          Crspq = Crspq - fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                    end if


                                    if (b==c)then
                                          Crspq = Crspq - fokian(a, d, TWONO, HNO(ad), n, NI, NA, NBasis)
                                    end if

                                    if (a==d)then
                                          Crspq = Crspq - fokian(b, c, TWONO, HNO(bc), n, NI, NA, NBasis)
                                    end if


                                    ! if (a==c)then
                                    !    if (b==d)then
                                    !       Crspq = Crspq - eorbs(a, TWONO, HNO(ac), NI, NA)
                                    !       Crspq = Crspq - eorbs(b, TWONO, HNO(bd), NI, NA)
                                    !       ! Crspq = Crspq +two*nu
                                    !    end if
                                    ! end if

                                    ! if (a==d)then
                                    !    if (b==c)then
                                    !       Crspq = Crspq - eorbs(a, TWONO, HNO(ad), NI, NA)
                                    !       Crspq = Crspq - eorbs(b, TWONO, HNO(bc), NI, NA)
                                    !       ! Crspq = Crspq +two*nu
                                    !    end if
                                    ! end if

                                    MxA(cd, ab) = -num_f* Crspq!!!! DUE TO TDA APPROX



                                    if (cd == ab) then
                                          e_trace_C = e_trace_C -num_f* Crspq
                                          MxS(cd, ab) = -one
                                    end if

                                    !  if (cd == ab) then
                                    !    print*, 'dodaje', ETot
                                    !    MxA(cd, ab) = MxA(cd, ab) + ETot
                                    ! end if


                                    ! print*,	''
                                    ! write(*,'(A8, 2I4, 3F11.4)')	'MxA C', cd, ab, -num_f ,  Crspq, MxA(cd, ab)
                                    !  print*,	''
                              end if

                              if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                    IndAux(a)==(2).and.IndAux(b)==(2))then

                                    num_f = one
                                    if (a==b) num_f = num_f * sqrt(frac12)
                                    if (c==d)num_f = num_f * sqrt(frac12)
                                    Brspq = TwoNO(NAddr3(a,c,b,d))+TwoNO(NAddr3(a,d,b,c))

                                    MxA(cd, ab) =  num_f * Brspq !!!! DUE TO TDA APPROX

                                    ! print*,	''
                                    ! write(*,'(A8, 2I4, 3F11.4)')	'MxA B_ijab', cd, ab, num_f ,  Brspq, MxA(cd, ab)
                                    !  print*,	''
                              end if

                              if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                    IndAux(a)==(0).and.IndAux(b)==(0))then

                                    num_f = one
                                    if (a==b)num_f = num_f * sqrt(frac12)
                                    if (c==d)num_f = num_f * sqrt(frac12)
                                    Brspq = TwoNO(NAddr3(a,c,b,d))+TwoNO(NAddr3(a,d,b,c))


                                    MxA(cd, ab) =  -num_f * Brspq!!!! DUE TO TDA APPROX
                                    ! print*,	''
                                    ! write(*,'(A8, 2I4, 3F11.4)') 'MxA B_abij', cd, ab, -num_f ,  Brspq, MxA(cd, ab)
                                    !  print*,	''

                              end if
                              ! print*, '-----------------------------------------------------------------------------'
                        end do j_colloops
                  end do i_rowloops
                  print*, 'sniezer', MxA(1,1)
            case(2)
                  print*, 'tak spin sym', spin_symm
                  ! no spin symm
                  e_trace_A = zero
                  e_trace_C = zero
                  i_rowloop: do i = 1, NDim
                        c = IndN(1, i)
                        d = IndN(2, i)
                        cd = IndX(i)
                        ! write(*, '(A5, 4I5)') 'indx', i, cd, c, d
                        j_colloop: do j = 1, NDim
                              a = IndN(1, j)
                              b = IndN(2, j)
                              ab = IndX(j)

                              if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                    IndAux(a)==(2).and.IndAux(b)==(2))then

                                    Arspq = zero
                                    ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                    bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)
                                    Arspq = TwoNO(NAddr3(a,c,b,d)) 

                                    ! wersja ogolna z fokianami      
                                    if (b==d)then
                                          Arspq = Arspq + fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                    end if

                                    if (a==c)then
                                          Arspq = Arspq + fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                    end if


                                    !                 if (a==c)then
                                    !                    if (b==d)then
                                    !                       Arspq = Arspq + eorbs(a, TWONO, HNO(ac), NI, NA)
                                    !                       Arspq = Arspq + eorbs(b, TWONO, HNO(bd), NI, NA)
                                    ! !                      Arspq = Arspq -two*nu
                                    !                    end if
                                    !                 end if

                                    MxA(cd, ab) = Arspq

                                    if (cd == ab) then
                                          e_trace_A = e_trace_A + Arspq
                                          MxS(cd, ab) = one
                                    end if
                                    ! if (cd > 7000.and.ab > 7000)then
                                    !    print*, 'cdab', cd, ab, size(MxA, dim=1), size(MxA, dim=2)
                                    ! end if
                              end if

                              if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                    IndAux(a)==(0).and.IndAux(b)==(0))then

                                    Crspq = zero
                                    ac = (max(a, c)*(max(a, c)-1))/2 + min(a, c)
                                    bd = (max(b, d)*(max(b, d)-1))/2 + min(b, d)

                                    Crspq =  TwoNO(NAddr3(a,c,b,d))!!!! DUE TO TDA APPROX




                                    ! wersja ogolna z fokianami                                                                                                                           
                                    if (b==d)then
                                          Crspq = Crspq - fokian(a, c, TWONO, HNO(ac), n, NI, NA, NBasis)
                                    end if

                                    if (a==c)then
                                          Crspq = Crspq - fokian(b, d, TWONO, HNO(bd), n, NI, NA, NBasis)
                                    end if


                                    ! if (a==c)then
                                    !    if (b==d)then
                                    !       Crspq = Crspq - eorbs(a, TWONO, HNO(ac), NI, NA)
                                    !       Crspq = Crspq - eorbs(b, TWONO, HNO(bd), NI, NA)
                                    !       !Crspq = Crspq + two* nu
                                    !    end if
                                    ! end if

                                    MxA(cd, ab) = -Crspq!!!! DUE TO TDA APPROX

                                    if (cd == ab) then
                                          e_trace_C = e_trace_C + Crspq
                                          MxS(cd, ab) = -one
                                    end if

                              end if

                              if (IndAux(c)==(0).and.IndAux(d)==(0).and.&
                                    IndAux(a)==(2).and.IndAux(b)==(2))then

                                    Brspq = TwoNO(NAddr3(a,c,b,d))

                                    MxA(cd, ab) = Brspq!!!! DUE TO TDA APPROX

                              end if

                              if (IndAux(c)==(2).and.IndAux(d)==(2).and.&
                                    IndAux(a)==(0).and.IndAux(b)==(0))then

                                    Brspq = TwoNO(NAddr3(a,c,b,d))
                                    if (abs(TwoNO(NAddr3(a,c,b,d))).gt.1.d-8)then
                                          write(*, '(A10, 4I5, 2F15.8)') 'integrals', c, d, a, b, TwoNO(NAddr3(a,c,b,d)), TwoNO(NAddr3(a,c,b,d))-TwoNO(NAddr3(a,d,b,c))
                                    end if

                                    MxA(cd, ab) = -Brspq!!!! DUE TO TDA APPROX

                              end if
                        end do j_colloop
                  end do i_rowloop
            end select

            print*, 'MACIERZ', NDim
            do i = 1, NDim
                  do j = 1, NDim
                        !                if (abs(MxA(i, j)).gt.1.d-8)then
                        write(*, '(2I5, F20.10)') i, j, MxA(i, j)
                        !               end if
                        ! write(*, '(36F8.5)') MxA(i, :)
                  end do
            end do

            print*, 'RDSC e_trace_A', e_trace_A
            print*, 'RDSC e_trace_C', e_trace_C

      end subroutine pp_rpa_yang_hf_tda


      subroutine compute_norm(AuxData, w, wi, rtz, v_plus, NDim, IndN, IndX, n, flaghh, Tpphh, spin_symm, Enuc)
            use math_constants
            type(TACppData) :: AuxData
            double precision, dimension(:), intent(inout) :: w, wi
            double precision, dimension(:,:), intent(inout) :: rtz
            integer, dimension(:), intent(in) :: v_plus
            integer, dimension(:), intent(in)    ::IndX
            integer, dimension(:,:), intent(in)  :: IndN
            double precision, dimension(:),    intent(in)   :: n
            integer, intent(in) :: NDim
            integer, intent(in) :: flaghh
            type(TDA_pphh), intent(in) :: Tpphh
            integer, intent(in) :: spin_symm
            double precision, intent(in) :: Enuc


            integer :: i, j
            integer :: r, s, rs, p, q, pq
            double precision, dimension(:,:), allocatable :: MxS, rtztemp
            double precision, dimension(:), allocatable :: tempx, witemp!, norm
            !    double precision :: norm_temp
            double precision, dimension(:), allocatable :: maxv
            integer, dimension(:), allocatable :: maxl
            integer, dimension(:), allocatable :: dy
            integer :: nvec
            integer :: n_start
            integer, parameter :: cf = 20
            integer, parameter :: cfm = 8
            integer, parameter :: cfs = 7
            integer :: offset_vo, offset_vo_vvoop, offset_2vo, offset_2vo_vvoo, offset_2vo_2vvoo
            integer :: offset_2vo_vvoop, offset_2vo_2vvoop, offset_2vo_2vvoop_vvoom
            logical :: cont
            integer :: N1dim
            double precision :: perc_S
            double precision :: e_eigen_plus, e_eigen_minus
            double precision :: sr
            N1dim = 0

            if (AuxData%general_version == GENVER_PP)then
                  allocate(MxS(NDim, NDim))
                  allocate(tempx(NDim))
                  allocate(dy(NDim))
                  allocate(witemp(NDim))
                  allocate(rtztemp(NDim, NDim))
                  nvec  = 7
                  allocate(maxl(nvec))
                  allocate(maxv(nvec))
                  MxS = zero
                  e_eigen_plus = zero
                  e_eigen_minus = zero


                  i_rowloops: do i = 1, NDim
                        dy(i) = i
                        r = IndN(1, i)
                        s = IndN(2, i)
                        rs = IndX(i)
                        j_colloops: do j = 1, NDim
                              p = IndN(1, j)
                              q = IndN(2, j)
                              pq = IndX(j)
                              if (p==r.and.q==s)then
                                    MxS(rs, pq) = 1 - n(r)-n(s)
                              end if
                        end do j_colloops
                  end do i_rowloops

                  select case (flaghh)
                  case(0)
                        call dsort(w, dy, NDim)
                        !          call msg('SORTED in Increasing order')
                  case (1)
                        !          call dsort_rev_pp(w, dy, NDim)
                        !         call msg('Sorted in Decreasing order')
                  end select

                  tempx = zero
                  do i = 1, NDim
                        write(*, '(F20.15)') w(i)
                        witemp(i) = wi(dy(i))
                        rtztemp(:, i) = rtz(:, dy(i))
                  end do
                  print*, 'sorted in eV'


                  print*, 'normedd'
                  do i = 1, NDim
                        if (v_plus(dy(i))==0)then
                              print*, 'odwzbudz sa dla', i, IndN(1, i), IndN(2, i)
                        end if

                  end do

                  write(*, '(2A15, A8, A11, 2A7)') 'RDSC KEY', 'omega', 'norm', 'maxv', 'p', 'q'

                  do i = 1, NDim
                        if (abs(witemp(i)) > 1d-10)then
                              print*, 'complex eingenvalue', witemp(i), ',', v_plus(dy(i))!, norm(i)
                        else
                              select case(flaghh)
                              case(0)
                                    if (v_plus(dy(i)) == 1)then
                                          e_eigen_plus = e_eigen_plus + w(i)
                                          call maxlocval(rtztemp(:,i), NDim, maxv, maxl, nvec, NDim, perc_S)
                                          !                   call msg(lfield('RDSC exc+', cf), w(i), fmt='F20.14'))
                                          do j = 1, nvec
                                                write(*, '(A40, 2F8.5, 2I4)')'', one, maxv(j), IndN(1, maxl(j)), IndN(2, maxl(j))
                                          end do
                                    else
                                          e_eigen_minus = e_eigen_minus + w(i)
                                    end if
                              case(1)
                                    if (v_plus(dy(i))==0)then
                                          call maxlocval(rtztemp(:,i), NDim, maxv, maxl, nvec, N1dim, perc_S)
                                          !                   call msg(lfield('RDSC exc', cf), w(i), fmt='F20.14'))
                                          do j = 1, nvec
                                                !                      call msg(lfield('', 40), 1.d+0, fmt='F8.5')//"   ", maxv(j), &
                                                !                          fmt='F8.5')//"   ", IndN(1, maxl(j)), fmt='I4')//"   ", IndN(2, maxl(j)), fmt='I4'))
                                          end do

                                    end if
                              end select
                        end if
                  end do
                  !       call msg(lfield('RDSC exend', cf))
                  deallocate(MxS)

                  deallocate(tempx)
                  deallocate(dy)
                  deallocate(witemp)
                  print*, 'RDSC e_eigen_plus', e_eigen_plus
                  print*, 'RDSC e_eigen_minus', e_eigen_minus


            else

                  allocate(dy(NDim))
                  allocate(witemp(NDim))
                  allocate(rtztemp(NDim, NDim))
                  nvec  = 10
                  allocate(maxl(nvec))
                  allocate(maxv(nvec))

                  do i = 1, NDim
                        dy(i) = i
                  end do

                  do i = 1, NDim
                        print*, i, w(i)
                  end do


                  call dsort(w, dy, NDim)
                  do i = 1, NDim
                        print*, i, w(i), dy(i)
                  end do

                  ! Print sorted depending on version
                  if (AuxData%general_version == GENVER_TDA_SING .or. AuxData%general_version == GENVER_TDA_TRIP)then

                        do i = 1, NDim
                              write(*, '(F20.15)') w(i)+Enuc
                              witemp(i) = wi(dy(i))
                              rtztemp(:, i) = rtz(:, dy(i))
                        end do

                        print*, 'sorted in eV'
                        do i = 1, NDim
                              print*, w(i) * 27.211399
                        end do
                  else if (AuxData%general_version == GENVER_RPA_SING .or. AuxData%general_version == GENVER_RPA_TRIP) then
                        do i = 1, NDim
                              if (w(i).gt.0.00001d+0)then
                                    write(*, '(F20.15)') w(i)
                              end if
                              witemp(i) = wi(dy(i))
                              rtztemp(:, i) = rtz(:, dy(i))
                        end do

                        print*, 'sorted in eV'
                        do i = 1, NDim
                              if (abs(w(i)).gt.0.000001d+0)then
                                    write(*, '(F20.15)') w(i) * 27.211399
                              end if
                        end do
                  end if


                  if (AuxData%general_version == GENVER_TDA_SING)then
                        offset_vo = 1 + Tpphh%Npair_vo
                  elseif (AuxData%general_version == GENVER_TDA_TRIP)then
                        offset_vo = 1 + Tpphh%Npair_vo
                        offset_vo_vvoop = 1 + Tpphh%Npair_vo+ Tpphh%Npair_vvoop
                  else if (AuxData%general_version  == GENVER_RPA_SING)then
                        offset_vo = Tpphh%Npair_vo
                        offset_2vo = 2 * Tpphh%Npair_vo
                        offset_2vo_vvoo = 2 * Tpphh%Npair_vo + Tpphh%Npair_vvoo
                  else if (AuxData%general_version  == GENVER_RPA_TRIP)then
                        offset_vo = Tpphh%Npair_vo
                        offset_2vo = 2 * Tpphh%Npair_vo
                        offset_2vo_vvoop = 2 * Tpphh%Npair_vo + Tpphh%Npair_vvoop
                        offset_2vo_2vvoop = 2 * Tpphh%Npair_vo + 2* Tpphh%Npair_vvoop
                        offset_2vo_2vvoop_vvoom = 2 * Tpphh%Npair_vo + 2* Tpphh%Npair_vvoop + Tpphh%Npair_vvoom
                  end if
                  write(*, '(2A15, A15, A11, 2A7)') 'RDSC KEY', 'omega   ',  'maxv', 'p', 'q'

                  do i = 1, NDim
                        if (AuxData%general_version == GENVER_RPA_SING.or.AuxData%general_version == GENVER_RPA_TRIP)then
                              cont = .false.
                              if (w(i).gt.0.00001d+0)then
                                    cont = .True.
                              end if
                        else
                              cont = .true.
                        end if

                        if (cont) then
                              if (abs(witemp(i)) > 1d-10)then
                                    print*, 'complex eingenvalue', witemp(i)
                              else

                                    if (AuxData%general_version == GENVER_RPA_SING.or.AuxData%general_version == GENVER_RPA_TRIP)then
                                          N1dim = offset_2vo
                                    else if (AuxData%general_version == GENVER_TDA_SING .or. AuxData%general_version == GENVER_TDA_TRIP)then
                                          N1dim = offset_vo
                                    end if

                                    if(AuxData%version == VER_PH_TDA_HF.or.AuxData%version==VER_PH_RPA_HF)then
                                          N1dim=NDim
                                    end if

                                    call maxlocval(rtztemp(:,i), NDim, maxv, maxl, nvec, N1dim, perc_S)
                                    if (AuxData%general_version == GENVER_RPA_SING.or.AuxData%general_version == GENVER_RPA_TRIP)then
                                          write(*, '(A10, F20.2, F20.14)')'RDSC exc+', perc_S, w(i)

                                    else
                                          write(*, '(A10, F20.2, F20.14)')'RDSC exc+', perc_S,	w(i)+Enuc
                                    end if
                                    do j = 1, nvec

                                          if (AuxData%general_version == GENVER_RPA_SING)then
                                                ! Tpphh%Npair_vo, Tpphh%Npair_vo, Tpphh%Npair_vvoo, Tpphh%Npair_vvoo
                                                if (maxl(j).le.offset_vo)then
                                                      write(*, '(A40, F8.5, 2I5, A10)')'', maxv(j), Tpphh%IndVO(1, maxl(j)), Tpphh%IndVO(2, maxl(j)), ' '
                                                else if (maxl(j).gt.offset_vo.and. maxl(j).le. offset_2vo)then
                                                      write(*, '(A40, F8.5, 2I5, A10)')'', maxv(j), Tpphh%IndVO(1, maxl(j)-offset_vo), Tpphh%IndVO(2, maxl(j)-offset_vo), ' '
                                                else if (maxl(j).gt.offset_2vo.and. maxl(j).le.offset_2vo_vvoo)then

                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j),&
                                                            Tpphh%Indpphh(1, maxl(j)-offset_2vo), &
                                                            Tpphh%Indpphh(2, maxl(j)-offset_2vo),&
                                                            Tpphh%Indpphh(3, maxl(j)-offset_2vo),&
                                                            Tpphh%Indpphh(4, maxl(j)-offset_2vo)

                                                else if (maxl(j).gt.offset_2vo_vvoo)then
                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j), Tpphh%Indpphh(1, maxl(j)-offset_2vo_vvoo), &
                                                            Tpphh%Indpphh(2, maxl(j)-offset_2vo_vvoo),&
                                                            Tpphh%Indpphh(3, maxl(j)-offset_2vo_vvoo),&
                                                            Tpphh%Indpphh(4, maxl(j)-offset_2vo_vvoo)
                                                end if

                                          elseif (AuxData%general_version == GENVER_RPA_TRIP)then

                                                ! Tpphh%Npair_vo, Tpphh%Npair_vo, Tpphh%Npair_vvoop, Tpphh%Npair_vvoop, Tpphh%Npair_vvoom, Tpphh%Npair_vvoom  
                                                if (maxl(j).le.offset_vo)then
                                                      write(*, '(A40, F8.5, 2I5, A10)')'', maxv(j), Tpphh%IndVO(1, maxl(j)), Tpphh%IndVO(2, maxl(j)), ' '
                                                else if (maxl(j).gt.offset_vo.and. maxl(j).le. offset_2vo)then
                                                      write(*, '(A40, F8.5, 2I5, A10)')'', maxv(j), Tpphh%IndVO(1, maxl(j)-offset_vo), Tpphh%IndVO(2, maxl(j)-offset_vo), ' '
                                                else if (maxl(j).gt.offset_2vo.and. maxl(j).le.offset_2vo_vvoop)then
                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j), Tpphh%Indpphhp(1, maxl(j)-offset_2vo), &
                                                            Tpphh%Indpphhp(2, maxl(j)-offset_2vo),&
                                                            Tpphh%Indpphhp(3, maxl(j)-offset_2vo),Tpphh%Indpphhp(4, maxl(j)-offset_2vo)
                                                else if (maxl(j).gt.offset_2vo_vvoop.and. maxl(j).le.offset_2vo_2vvoop)then
                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j), Tpphh%Indpphhp(1, maxl(j)-offset_2vo_vvoop), &
                                                            Tpphh%Indpphhp(2, maxl(j)-offset_2vo_vvoop),&
                                                            Tpphh%Indpphhp(3, maxl(j)-offset_2vo_vvoop),Tpphh%Indpphhp(4, maxl(j)-offset_2vo_vvoop)
                                                else if (maxl(j).gt.offset_2vo_2vvoop.and. maxl(j).le.offset_2vo_2vvoop_vvoom)then
                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j), Tpphh%Indpphhm(1, maxl(j)-offset_2vo_2vvoop), &
                                                            Tpphh%Indpphhm(2, maxl(j)-offset_2vo_2vvoop),&
                                                            Tpphh%Indpphhm(3, maxl(j)-offset_2vo_2vvoop),Tpphh%Indpphhm(4, maxl(j)-offset_2vo_2vvoop)
                                                else if (maxl(j).gt.offset_2vo_2vvoop_vvoom)then
                                                      write(*, '(A40, F8.5, 4I5)')'', maxv(j), Tpphh%Indpphhm(1, maxl(j)-offset_2vo_2vvoop_vvoom), &
                                                            Tpphh%Indpphhm(2, maxl(j)-offset_2vo_2vvoop_vvoom),&
                                                            Tpphh%Indpphhm(3, maxl(j)-offset_2vo_2vvoop_vvoom),&
                                                            Tpphh%Indpphhm(4, maxl(j)-offset_2vo_2vvoop_vvoom)
                                                end if

                                          else if (AuxData%general_version == GENVER_TDA_SING)then
                                                ! 1, Tpphh%Npair_vo, Tpphh%Npair_vvoo
                                                if (maxl(j)==1)then
                                                      write(*, '(A40, F8.5, 2I4)')'', maxv(j), Tpphh%IndVO(1, maxl(j)), Tpphh%IndVO(2, maxl(j))
                                                elseif (maxl(j).le.offset_vo.and.maxl(j)>1)then
                                                      write(*, '(A40, F8.5, 2I4)')'', maxv(j), Tpphh%IndVO(1, maxl(j)-1), Tpphh%IndVO(2, maxl(j)-1)
                                                else if (maxl(j).gt.offset_vo)then
                                                      write(*, '(A40, F8.5, 4I4)')'', maxv(j), Tpphh%Indpphh(1, maxl(j)-offset_vo), &
                                                            Tpphh%Indpphh(2, maxl(j)-offset_vo),&
                                                            Tpphh%Indpphh(3, maxl(j)-offset_vo),Tpphh%Indpphh(4, maxl(j)-offset_vo)
                                                end if

                                          else if (AuxData%general_version == GENVER_TDA_TRIP)then
                                                ! 1, Tpphh%Npair_vo, Tpphh%Npair_vvoop, Tpphh%Npair_vvoom
                                                if (maxl(j)==1)then
                                                      write(*, '(A40, F8.5, 2I4)')'', maxv(j), Tpphh%IndVO(1, maxl(j)), Tpphh%IndVO(2, maxl(j))
                                                elseif (maxl(j).le.offset_vo.and.maxl(j)>1)then
                                                      write(*, '(A40, F8.5, 2I4)')'', maxv(j), Tpphh%IndVO(1, maxl(j)-1), Tpphh%IndVO(2, maxl(j)-1)
                                                else if (maxl(j).gt.offset_vo.and. maxl(j).le.offset_vo_vvoop)then
                                                      write(*, '(A40, F8.5, 4I4)')'', maxv(j), Tpphh%Indpphhp(1, maxl(j)-offset_vo), &
                                                            Tpphh%Indpphhp(2, maxl(j)-offset_vo),&
                                                            Tpphh%Indpphhp(3, maxl(j)-offset_vo),Tpphh%Indpphhp(4, maxl(j)-offset_vo)
                                                else if (maxl(j).gt.offset_vo_vvoop)then
                                                      write(*, '(A40, F8.5, 4I4)')'', maxv(j), Tpphh%Indpphhm(1, maxl(j)-offset_vo_vvoop), &
                                                            Tpphh%Indpphhm(2, maxl(j)-offset_vo_vvoop),&
                                                            Tpphh%Indpphhm(3, maxl(j)-offset_vo_vvoop),Tpphh%Indpphhm(4, maxl(j)-offset_vo_vvoop)
                                                end if
                                          end if

                                    end do
                              end if
                        end if
                  end do
                  ! end if
                  !    call msg(lfield('RDSC exend', cf))
                  deallocate(dy)
                  deallocate(witemp)
                  deallocate(rtztemp)
                  deallocate(maxl)
                  deallocate(maxv)

            end if


      end subroutine compute_norm

      subroutine nonsymmetric_eigenproblem_pp(wr, wi, vr, A, S, v_plus, IndN, IndAux, AC_TYPE)
            double precision, dimension(:), intent(out)      :: wr
            double precision, dimension(:), intent(out)      :: wi
            double precision, dimension(:, :), allocatable   :: vl
            double precision, dimension(:, :), intent(out)   :: vr
            double precision, dimension(:, :), intent(inout) :: A, S
            integer, dimension(:), intent(out) :: v_plus
            integer, dimension(:,:), intent(in) :: indN
            integer, dimension(:),intent(in) :: IndAux
            integer, intent(in) :: AC_TYPE
            !    double precision, optional, intent(out) :: shift
            double precision :: shift, max_hh, min_pp
            double precision, dimension(:),allocatable :: S_diag

            double precision, dimension(1) :: work0
            double precision, dimension(:), allocatable :: work, tempx
            double precision, dimension(:, :), allocatable :: vr_ort
            integer :: lwork, info
            integer :: n, i, j, k
            double precision :: norm, norm_l, norm_rl
            double precision :: dd, ddm, dd_comp
            type (tclock) :: timer
            double precision, parameter :: tol = 1.d-4
            double precision, dimension(:), allocatable :: wr_copy
            integer, dimension(:), allocatable :: StartIdx, EndIdx
            integer, dimension(:), allocatable :: dy
            integer :: count, countj, ss
            logical :: flag_change

            double precision, dimension(:,:), allocatable :: miniA
            double precision, dimension(:,:), allocatable :: vl2, vr2
            double precision, dimension(:), allocatable :: wr2, wi2
            integer :: ii, jj, n2, r, s2, p,q


            external :: dgeev


            n = size(A, dim=1)
            allocate(StartIdx(n))
            allocate(EndIdx(n))
            allocate(dy(n))
            allocate(wr_copy(n))
            allocate(S_diag(n))

            allocate(vl(n,n))
            vl = zero
            vr = zero

            do i = 1, n
                  S_diag(i) = S(i, i)
            end do


            call clock_start(timer)    
            lwork = -1
            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work0, lwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))

            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work, lwork, info)
            !    print*, 'info', info
            if (info /= 0) then
                  print*, "Nonsymmetric matrix eigendecompositino failed with info="
                  ! call msg("Nonsymmetric matrix eigendecompositino failed with info=" // str(info), MSG_ERROR)
                  error stop
            end if
            print*, 'Czas diagonalizacji real: ', clock_readwall(timer)

            allocate(tempx(n))
            call clock_start(timer)
            v_plus = 0

            print*, ''
            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        if (wi(i).gt.0.d+0)then
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)
                              call ddot_norm(vr(:, i+1),S,  vr(:, i+1), n, dd_comp)
                              print*, '>0 dd', i, i+1, dd
                              print*, '>0 dd_comp', dd_comp
                        else
                              call ddot_norm(vr(:, i-1),S,  vr(:, i-1), n, dd)
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd_comp)
                              print*, '<0 dd', i-1, i, dd
                              print*, '<0 dd_comp',	dd_comp

                        end if
                        dd = dd + dd_comp
                        v_plus(i)=2
                  else
                        call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)
                        if (dd.gt.zero)then
                              v_plus(i) = 1
                        end if

                  end if
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'i-norm', i, dd
                  end if
            end do

            if (AC_TYPE == PPAC)then
                  print*, 'AC_TYPE: PP'
            else if (AC_TYPE == HHAC) then
                  print*, 'AC_TYPE: HH'
            end if


            call clock_start(timer)
            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'COMPLEX EIGENVALUES', i, wr(i), wi(i)
                        !         stop
                  end if
            end do

            wr_copy = wr
            do i = 1, n
                  dy(i)=i
            end do
            call dsort(wr_copy, dy, n)

            shift = zero
            do i = 1, n
                  if (v_plus(dy(i))== 0 )then
                        max_hh = wr_copy(i)
                  end if
                  if (v_plus(dy(i))== 1 )then
                        min_pp = wr_copy(i)
                        if (i == 1) then
                              shift = min_pp/two
                              print*, 'nohh'
                        else
                              shift = (abs(min_pp) + abs(max_hh))/two
                              print*, 'max_hh', max_hh
                        end if
                        print*, 'min_pp' , wr_copy(i)
                        print*, 'shift', shift
                        exit
                  else

                  end if
            end do



            wr_copy = zero
            dy = 0
            countj=1
            do i = 1, n
                  if (v_plus(i) == AC_TYPE) then
                        wr_copy(countj) = wr(i)
                        dy(countj) = i
                        countj = countj + 1
                  end if
            end do
            countj = countj-1
            call dsort(wr_copy(1:countj), dy(1:countj), countj)



            allocate(vr_ort(n,n))

            StartIdx = 0
            EndIdx = 0
            count = 1
            j = 0


            do i = 1, countj
                  if (j==0)then
                        StartIdx(count) = i
                        if (i.eq.countj)then
                              EndIdx(count) = i
                        end if
                        j = 1
                  else
                        if (abs(wr_copy(i)-wr_copy(i-1)).lt.tol)then
                              if (i==n)then
                                    EndIdx(count) = i
                              end if
                              if (i==countj)then
                                    EndIdx(count) = i
                              end if
                        else
                              EndIdx(count) = i-1

                              if (i.ne.countj)then
                                    count = count + 1
                                    StartIdx(count)=i
                              else
                                    count = count+1
                                    StartIdx(count) = i
                                    EndIdx(count) = i
                              end if

                        end if
                  end if
            end do


            if (StartIdx(1) ==0)then
                  count = 0
            end if

            call Orthogonalize(vr, count, StartIdx(1:count), EndIdx(1:count), S_diag, dy, AC_TYPE)


      end subroutine nonsymmetric_eigenproblem_pp


      subroutine sss(MxH, NDim, NBasis, n, IndN, IndX, IndAux, TwoNO, HNO, NI, NA, spin_symm, Tpphh, Etot, version)
            use math_constants
            double precision, dimension(:,:), intent(inout) :: MxH
            double precision, dimension(:),      intent(in) :: TwoNO ,HNO
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in)    ::IndX
            integer, dimension(:,:), intent(in)  :: IndN
            integer, dimension(:), intent(in)    :: IndAux
            integer, intent(in) :: NDim, NBasis
            integer, intent(in) :: NI, NA
            integer, intent(in) :: spin_symm
            type(TDA_pphh), intent(in) :: Tpphh
            double precision, intent(in) :: Etot
            integer, intent(in) :: version
            type(log_cond) :: LC
            double precision :: Hpqrs

            integer :: IH1, IH2, IH3, IH4
            logical :: Lik, Lac
            integer :: p, q, r, s, pq, rs, i, j
            integer :: a, b, c, d, k, l
            integer :: NIA
            integer :: offset1, offset2, offset3
            integer :: offset_vo, offset_2vo, offset_2vo_vvoop, offset_2vo_2vvoop, offset_2vo_2vvoop_vvoo
            integer :: offset, offset0
            integer, external :: NAddr3
            integer :: kkk
            double precision :: norm


            ! TA PROCEDURA ZWRACA S^(-1/2) dla wersji tda i S^-1 dla rpa
            NIA = NI + NA

            if (version==VER_PH_PPHH_TDA_HF .or. version==VER_PH_TDA_HF .or. version == VER_PPHH_TDA_HF)then
                  offset = 1+Tpphh%Npair_vo
                  offset0 = 0
                  MxH = zero
                  MxH(1, 1) = one/sqrt(frac12)
                  do IH1 = 2, Tpphh%Npair_vo+1
                        print*, 'wstawiam', IH1, IH1, one
                        MxH(IH1, IH1) = one
                  end do

                  !print*, 'zaczynam block 22'
                  H22_1loop: do IH1 = 1, Tpphh%Npair_vvoo
                        i = Tpphh%Indpphh(1, IH1)
                        j = Tpphh%Indpphh(2, IH1)
                        a = Tpphh%Indpphh(3, IH1)
                        b = Tpphh%Indpphh(4, IH1)

                        H22_2loop: do IH2 = 1, Tpphh%Npair_vvoo
                              k = Tpphh%Indpphh(1, IH2)
                              l = Tpphh%Indpphh(2, IH2)
                              c = Tpphh%Indpphh(3, IH2)
                              d = Tpphh%Indpphh(4, IH2)

                              Hpqrs = zero

                              call conditions(a, b, i, j, c, d, k, l, LC)
                              norm = zero

                              if ((IH1 .ge. 1 .and. IH1 .le. Tpphh%vo1).and.&
                                    (IH2 .ge. 1 .and. IH2 .le. Tpphh%vo1))then

                                    if(LC%Lacbd .and. LC%Likjl)then
                                          norm = norm + two
                                    end if
                              else
                                    if(LC%Lacbd .and. LC%Likjl)then
                                          norm = norm + 1
                                    end if
                                    if(LC%Ladbc .and. LC%Liljk)then
                                          norm = norm + 1
                                    end if
                              end if


                              if (abs(norm).gt.1.d-1)then
                                    norm = sqrt(one/norm)
                              end if

                              MxH(IH1+offset, IH2+offset) = norm
                              ! if (abs(norm).gt.1.d-1)then
                              !    write(*, '(A10, 2I5, F20.12)') 'MxH22', &
                              !         IH1+offset, IH2+offset, norm
                              ! end if
                        end do H22_2loop
                  end do H22_1loop
            else if(version==VER_PH_RPA_HF) then
                  MxH = zero

                  offset = Tpphh%Npair_vo

                  H1loop: do IH1 = 1, Tpphh%Npair_vo
                        MxH(IH1, IH1) = one
                        MxH(IH1+offset, IH1 + offset) = -one
                  end do H1loop
            else if(version==VER_PH_RPA_TDA_HF) then
                  MxH = zero

                  H1loop2: do IH1 = 1, Tpphh%Npair_vo
                        MxH(IH1, IH1) = one
                  end do H1loop2

            else if (version==VER_PP_RPA_MULTI .or.version == VER_HH_RPA_MULTI .or. version == VER_PP_RPA_TDA_HF)then
                  MxH = zero
                  i_rowloop: do i = 1, NDim
                        r = IndN(1, i)
                        s = IndN(2, i)
                        rs = IndX(i)
                        j_colloop: do j = 1, NDim
                              p = IndN(1, j)
                              q = IndN(2, j)
                              pq = IndX(j)
                              if ((q==s).and.(r==p))then
                                    ! write(*, '(A6, 4I5, F20.12)') 'MxH22', &
                                    !      rs, pq, r, s, 1-n(r)-n(s)
                                    MxH(rs, pq) = one/(1-n(r)-n(s))
                                    !                if (abs(MxH(rs, pq)).gt.1.d+0)then
                                    !                   print*, 'zzz', rs, pq, r, s, MxH(rs, pq), n(r), n(s), MxH(119, 119)

                                    !               end if
                              end if
                        end do j_colloop
                  end do i_rowloop

            else if (version == VER_PH_PPHH_RPA_HF .or.version == VER_PH_PPHH_RPA_TDA_HF.or.version == VER_PPHH_RPA_HF.or.version == VER_PPHH_RPA_TDA_HF) then

                  if (spin_symm ==0)then

                        MxH = zero

                        if (version == VER_PH_PPHH_RPA_HF .or.version == VER_PH_PPHH_RPA_TDA_HF)then
                              offset1 = Tpphh%Npair_vo
                              offset2 = 2*Tpphh%Npair_vo
                              offset3 = 2*Tpphh%Npair_vo + Tpphh%Npair_vvoo
                        else if (version == VER_PPHH_RPA_HF .or.version == VER_PPHH_RPA_TDA_HF)then
                              offset1 = 0
                              offset2 = 0
                              offset3 = Tpphh%Npair_vvoo
                        end if
                        if (version == VER_PH_PPHH_RPA_HF .or.version == VER_PH_PPHH_RPA_TDA_HF)then
                              H11_1loop: do IH1 = 1, Tpphh%Npair_vo
                                    a = Tpphh%IndVO(1, IH1)
                                    i = Tpphh%IndVO(2, IH1)
                                    H1_1loop: do IH2 = 1, Tpphh%Npair_vo
                                          c = Tpphh%IndVO(1, IH2)
                                          k = Tpphh%IndVO(2, IH2)

                                          if ((a==c).and.(i==k))then
                                                MxH(IH1, IH2) = one
                                                MxH(IH1+offset1, IH2+offset1) = -one
                                          end if
                                    end do H1_1loop
                              end do H11_1loop
                        end if

                        H22_1loopp: do IH1 = 1, Tpphh%Npair_vvoo
                              i = Tpphh%Indpphh(1, IH1)
                              j = Tpphh%Indpphh(2, IH1)
                              a = Tpphh%Indpphh(3, IH1)
                              b = Tpphh%Indpphh(4, IH1)

                              H22_2loopp: do IH2 = 1, Tpphh%Npair_vvoo
                                    k = Tpphh%Indpphh(1, IH2)
                                    l = Tpphh%Indpphh(2, IH2)
                                    c = Tpphh%Indpphh(3, IH2)
                                    d = Tpphh%Indpphh(4, IH2)

                                    if ((a==c).and. (b==d).and. (i==k).and.(j==l))then
                                          MxH(IH1+offset2, IH2+offset2) = one
                                          MxH(IH1+offset3, IH2+offset3) = -one
                                          !                   write(*,'(A5, 10I5)') 'dup1', IH1, IH2, a, b, c, d, i, j, k, l
                                    end if

                                    if ((a==d).and. (b==c).and. (i==l).and.(j==k))then
                                          MxH(IH1+offset2, IH2+offset2) = one
                                          MxH(IH1+offset3, IH2+offset3) = -one
                                          !                  write(*,'(A5, 10I5)') 'dup2', IH1, IH2, a, b, c, d, i, j, k, l
                                    end if
                              end do H22_2loopp
                        end do H22_1loopp

                  else if (spin_symm == 1)then
                        MxH = zero

                        if (version == VER_PH_PPHH_RPA_HF .or.version == VER_PH_PPHH_RPA_TDA_HF)then
                              offset_vo = Tpphh%Npair_vo
                              offset_2vo = 2*Tpphh%Npair_vo
                              offset_2vo_vvoop = 2*Tpphh%Npair_vo + Tpphh%Npair_vvoop
                              offset_2vo_2vvoop = 2*Tpphh%Npair_vo + 2*Tpphh%Npair_vvoop
                              offset_2vo_2vvoop_vvoo = 2*Tpphh%Npair_vo + 2*Tpphh%Npair_vvoop + Tpphh%Npair_vvoom

                              offset1 = Tpphh%Npair_vo
                              offset2 = 2*Tpphh%Npair_vo
                              offset3 = 2*Tpphh%Npair_vo + Tpphh%Npair_vvoo
                              ! else if (version == VER_PPHH_RPA_HF .or.version == VER_PPHH_RPA_TDA_HF)then
                              !    offset1 = 0
                              !    offset2 = 0
                              !    offset3 = Tpphh%Npair_vvoo
                        end if
                        if (version == VER_PH_PPHH_RPA_HF .or.version == VER_PH_PPHH_RPA_TDA_HF)then

                              H11_1loopa: do IH1 = 1, Tpphh%Npair_vo
                                    a = Tpphh%IndVO(1, IH1)
                                    i = Tpphh%IndVO(2, IH1)
                                    H1_1loopb: do IH2 = 1, Tpphh%Npair_vo
                                          c = Tpphh%IndVO(1, IH2)
                                          k = Tpphh%IndVO(2, IH2)

                                          if ((a==c).and.(i==k))then
                                                MxH(IH1, IH2) = one
                                                MxH(IH1+offset_vo, IH2+offset_vo) = -one
                                          end if
                                    end do H1_1loopb
                              end do H11_1loopa
                        end if

                        H22_1looppa: do IH1 = 1, Tpphh%Npair_vvoop
                              i = Tpphh%Indpphhp(1, IH1)
                              j = Tpphh%Indpphhp(2, IH1)
                              a = Tpphh%Indpphhp(3, IH1)
                              b = Tpphh%Indpphhp(4, IH1)

                              H22_2looppb: do IH2 = 1, Tpphh%Npair_vvoop
                                    k = Tpphh%Indpphhp(1, IH2)
                                    l = Tpphh%Indpphhp(2, IH2)
                                    c = Tpphh%Indpphhp(3, IH2)
                                    d = Tpphh%Indpphhp(4, IH2)

                                    if ((a==c).and. (b==d).and. (i==k).and.(j==l))then
                                          MxH(IH1+offset_2vo, IH2+offset_2vo) = one
                                          MxH(IH1+offset_2vo_vvoop, IH2+offset_2vo_vvoop) = -one
                                    end if

                                    ! if ((a==d).and. (b==c).and. (i==l).and.(j==k))then
                                    !    MxH(IH1+offset_2vo, IH2+offset_2vo) = one
                                    !    MxH(IH1+offset_2vo_vvoop, IH2+offset_2vo_vvoop) = -one
                                    ! end if
                              end do H22_2looppb
                        end do H22_1looppa

                        H22_1looppa1: do IH1 = 1, Tpphh%Npair_vvoom
                              i = Tpphh%Indpphhm(1, IH1)
                              j = Tpphh%Indpphhm(2, IH1)
                              a = Tpphh%Indpphhm(3, IH1)
                              b = Tpphh%Indpphhm(4, IH1)

                              H22_2looppb1: do IH2 = 1, Tpphh%Npair_vvoom
                                    k = Tpphh%Indpphhm(1, IH2)
                                    l = Tpphh%Indpphhm(2, IH2)
                                    c = Tpphh%Indpphhm(3, IH2)
                                    d = Tpphh%Indpphhm(4, IH2)

                                    if ((a==c).and. (b==d).and. (i==k).and.(j==l))then
                                          MxH(IH1+offset_2vo_2vvoop, IH2+offset_2vo_2vvoop) = one
                                          MxH(IH1+offset_2vo_2vvoop_vvoo, IH2+offset_2vo_2vvoop_vvoo) = -one
                                    end if
                                    if ((a==d).and. (b==c).and. (i==l).and.(j==k))then
                                          MxH(IH1+offset_2vo_2vvoop, IH2+offset_2vo_2vvoop) = -one
                                          MxH(IH1+offset_2vo_2vvoop_vvoo, IH2+offset_2vo_2vvoop_vvoo) = one
                                    end if
                              end do H22_2looppb1
                        end do H22_1looppa1
                  end if
            end if


            ! print*, 'mxhhhhhhhhhh', MxH(119, 119)
            ! stop



      end subroutine sss

      function fokian(i, j, TWONO, h, n, NI, NA, NBasis)
            use math_constants
            double precision :: fokian
            integer, intent(in) :: i, j
            double precision, dimension(:), intent(in) :: TwoNO
            double precision, intent(in) :: h
            double precision, dimension(:), intent(in) ::n
            integer, intent(in) :: NI, NA, NBasis
            integer :: t
            integer :: NIA
            integer, external :: NAddr3

            NIA = NI+NA
            fokian = h
            do t = 1, NBasis
                  fokian = fokian + n(t) * (two* TwoNO(NAddr3(i, j, t, t)) - TwoNO(NAddr3(i, t, t, j)))
            end do

      end function fokian

      subroutine conditions(a, b, i, j, c, d, k, l, LC)

            integer, intent(in) :: a, b, i, j, c, d, k, l
            type(log_cond), intent(out) :: LC


            LC%Lil = (i==l)
            LC%Ljk = (j==k)
            LC%Liljk = (LC%Lil .and. LC%Ljk)

            LC%Lik = (i==k)
            LC%Ljl = (j==l)
            LC%Likjl	= (LC%Lik .and. LC%Ljl)

            LC%Lac = (a==c)
            LC%Lbd = (b==d)
            LC%Lacbd	= (LC%Lac .and. LC%Lbd)

            LC%Lad = (a==d)
            LC%Lbc = (b==c)
            LC%Ladbc = (LC%Lad .and. LC%Lbc)

            LC%Lacik = (LC%Lac .and. LC%Lik)
            LC%Ladil = (LC%Lad .and. LC%Lil)
            LC%Lacil = (LC%Lac .and. LC%Lil)
            LC%Ladik = (LC%Lad .and. LC%Lik)

            LC%Lacjk = (LC%Lac .and. LC%Ljk)
            LC%Lacjl = (LC%Lac .and. LC%Ljl)
            LC%Ladjk = (LC%Lad .and. LC%Ljk)
            LC%Ladjl = (LC%Lad .and. LC%Ljl)

            LC%Lbcik = (LC%Lbc .and. LC%Lik)
            LC%Lbdil = (LC%Lbd .and. LC%Lil)
            LC%Lbcil = (LC%Lbc .and. LC%Lil)
            LC%Lbdik = (LC%Lbd .and. LC%Lik)

            LC%Lbcjk = (LC%Lbc .and. LC%Ljk)
            LC%Lbcjl = (LC%Lbc .and. LC%Ljl)
            LC%Lbdjk = (LC%Lbd .and. LC%Ljk)
            LC%Lbdjl = (LC%Lbd .and. LC%Ljl)


      end subroutine conditions


      subroutine maxlocval(v, N, maxv, maxl, k, N1dim, perc_S)
            use math_constants
            integer, dimension(:), intent(out) :: maxl
            double precision, dimension(:), intent(out) :: maxv
            double precision, dimension(:), intent(inout) :: v
            integer, intent(in) :: N1dim
            integer, intent(in) :: N
            integer, intent(in) :: k
            double precision, intent(out) :: perc_S
            integer :: i
            integer :: j

            perc_S = zero

            do j = 1,  N1dim
                  perc_S = perc_S+ v(j)**2
                  !       print*, v(j)
            end do


            do j = 1, k
                  maxl(j) = 1
                  maxv(j) = abs(v(1))


                  do i = 2, N
                        if(abs(v(i)) > maxv(j))then

                              maxl(j) = i
                              maxv(j) = abs(v(i))
                        end if
                  end do

                  v(maxl(j)) = zero
            end do

            do j = 1, k
                  v(maxl(j)) = maxv(j)
            end do

      end subroutine maxlocval


      subroutine fill_IndAuxMod(Occ, IndAux, IndMod, Nbasis, NI, NA, NV)
            use math_constants
            double precision, dimension(:), intent(in) :: Occ
            integer, dimension(:), intent(out) :: IndAux, IndMod
            integer, intent(in) :: Nbasis
            integer, intent(in) :: NI, NA
            integer, intent(out) :: NV

            integer :: i, k
            integer :: gt0
            integer :: u

            IndAux = 0
            IndMod = 0
            do i = 1, NI
                  IndAux(i) = 0
            end do
            do i = NI+1, NI+NA
                  IndAux(i) = 1
            end do
            do i = NI+NA+1, Nbasis
                  IndAux(i)=2
            end do

            k = 1
            do i = 1, nbasis
                  if (IndAux(i)  == 0)then
                        IndMod(k) = i
                        k = k + 1
                  end if
            end do
            do i = 1, nbasis
                  if (IndAux(i) == 1) then
                        IndMod(k) = i
                        k = k + 1
                  end if
            end do
            do i = 1, nbasis
                  if (IndAux(i) == 2) then
                        IndMod(k) = i
                        k = k + 1
                  end if
            end do

            NV = Nbasis-NI-NA

      end subroutine fill_IndAuxMod



end module ppac_subs
