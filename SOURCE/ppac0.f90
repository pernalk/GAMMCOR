module ppac0

      use iso_fortran_env

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals        
      use types
      use lin
      use clock

      use math_constants
      use quadratures
      use ppac_types
      use pp_utils

      implicit none

      double precision, parameter :: tol = 1.d-4
      integer, parameter,private ::  IOO=1, IVV=2, IAA=3, IVA=4, IAO=5
      integer, parameter,private :: IOOT=6, IVVT=7, IAAT=8, IVAT=9, IAOT = 10
      integer, parameter, private :: VVOO=1, VVAO=2, VVAA=3, VAOO=4
      integer, parameter, private :: VAAO=5, VAAA=6, AAOO=7, AAAO=8
      integer, parameter, private :: VAAO2 = 9




contains
      subroutine ACPP0_fast(H_type, ECorr_AC0, &
            ETot, ENuc, n, XOne, TwoNO, AuxData, IndN, IndAux, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, Flags, AC_TYPE)

            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: H_type
            double precision, intent(inout) :: ECorr_AC0, ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), allocatable :: Ind2
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type (tclock) :: timer, timer0
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: AC_TYPE
            double precision, dimension(:), allocatable :: TwoNOA, Ha

            integer :: spin_symm

            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9

            character(len=3), dimension(10) :: bl_name
            character(len=5), dimension(9)  :: bbl_name
            integer, dimension(2,9) :: block_pair
            type(tMxA) :: MxA_s, MxA_t, MxS_s, MxS_t
            integer :: p, q, r, s, i, j
            integer :: twoint_dim
            integer :: dimw, spsym
            integer :: ir, il
            integer, dimension(10) :: dim_st
            double precision, dimension(9) :: E_contr_s, E_contr_t
            type(eigsBlockParams), dimension(10) :: eigsPA

            integer, external :: NAddr3


            ! 1oo, 2vv, 3aa, 4va, 5ao
            ! 6oot, 7vvt, 8aat, 9vat, 10aot

            bbl_name(1) = 'vvoo'
            bbl_name(2) = 'vvao'
            bbl_name(3) = 'vvaa'

            bbl_name(4) = 'vaoo'
            bbl_name(5) = 'vaao'
            bbl_name(6) = 'vaaa'

            bbl_name(7) = 'aaoo'
            bbl_name(8) = 'aaao'
            bbl_name(9) = 'vaao2'

            allocate(Ind2(NBasis))

            allocate(Ha(NInte1))
            twoint_dim = size(TwoNO, dim=1)
            allocate(TwoNOA(twoint_dim))



            call clock_start(timer0)    
            call PPERPA_init(H_type, ETot, ENuc,  n, XOne, TwoNO, Ind2, IndAux, &
                  NBasis, NA, NI, NV, NInte1, NInte2, zero, Flags, TwoNOA, Ha)
            print*, 'TIME na PPERPA_init ', clock_readwall(timer0)

            bl_name(IOO) = 'oo'
            bl_name(IVV) = 'vv'
            bl_name(IAA) = 'aa'
            bl_name(IVA) = 'va'
            bl_name(IAO) = 'ao'

            bl_name(IOOT) = 'oot'
            bl_name(IVVT) = 'vvt'
            bl_name(IAAT) = 'aat'
            bl_name(IVAT) = 'vat'
            bl_name(IAOT) = 'aot'


            allocate(eigsPA(IOO)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVV)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAA)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVA)%IndN(2, NV*NA))
            allocate(eigsPA(IAO)%IndN(2, NA*NI))
            allocate(eigsPA(IOOT)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT)%IndN(2, NA*NI))

            do i = 1, 10
                  eigsPA(i)%IndN = 0
                  eigsPA(i)%dim = 1
            end do

            do i = 1, AuxData%NDim_s
                  p = IndN(1, i)
                  q = IndN(2, i)

                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(p, q, eigsPA(IOO)%IndN, eigsPA(IOOT)%IndN, eigsPA(IOO)%dim, eigsPA(IOOT)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(p, q, eigsPA(IVV)%IndN, eigsPA(IVVT)%IndN, eigsPA(IVV)%dim, eigsPA(IVVT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(p, q, eigsPA(IAA)%IndN, eigsPA(IAAT)%IndN, eigsPA(IAA)%dim, eigsPA(IAAT)%dim)

                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(p, q, eigsPA(IVA)%IndN, eigsPA(IVAT)%IndN, eigsPA(IVA)%dim, eigsPA(IVAT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(p, q, eigsPA(IAO)%IndN, eigsPA(IAOT)%IndN, eigsPA(IAO)%dim, eigsPA(IAOT)%dim)
                  end if
            end do

            do i = 1, 10
                  eigsPA(i)%dim = eigsPA(i)%dim-1
            end do
            ! 1oo, 2vv, 3aa, 4va, 5ao                                                                                                          
            ! 6oot, 7vvt, 8aat, 9vat, 10aot  
            do i = 1, 10
                  dimw = eigsPA(i)%dim
                  allocate(eigsPA(i)%Eig(dimw))
                  if (i == IOO.or.i==IVV.or.i==IOOT.or.i==IVVT)then
                        allocate(eigsPA(i)%Eigvec(1, 1))
                  else
                        allocate(eigsPA(i)%Eigvec(dimw, dimw))
                  end if
                  allocate(eigsPA(i)%v_plus(dimw))
            end do


            call clock_start(timer0)
            do i = 1, 10
                  print*, ''
                  print*, 'teraz block', i
                  if (i.le.5)then
                        spsym = 0
                  else
                        spsym = 1
                  end if
                  if (eigsPA(i)%dim.gt.0)then

                        call eigs0_block(i, eigsPA(i)%Eig, eigsPA(i)%Eigvec, eigsPA(i)%v_plus, n, Ha, TwoNOA, &
                              eigsPA(i)%IndN, Ind2, IndAux, IndMod, eigsPA(i)%dim, NBasis, NA, NI, NV, NInte1, &
                              NInte2, spsym, Flags)

                        if (i>5)then
                              print*, 'wartosci wlasne tego bloku', i
                              do j = 1, eigsPA(i)%dim
                                    if (abs(eigsPA(i)%Eig(j)).gt.1.d-5)then
                                          print*,  eigsPA(i)%Eig(j)
                                    end if
                              end do
                        end if
                  else
                        print*, 'block', bl_name(i), 'dimension is 0'
                  end if
            end do
            print*, 'koniec blokow eigs0'
            print*, 'TIME na all blocks ', clock_readwall(timer0)

            !    stop ! pamietaj odkomentowac parallel

            call clock_start(timer0)
            call PPERPA_init(H_type, ETot, ENuc,  n, XOne, TwoNO, Ind2, IndAux, &
                  NBasis, NA, NI, NV, NInte1, NInte2, one, Flags, TwoNOA, Ha)
            print*, 'TIME na PPERPA_init alphaone ', clock_readwall(timer0)
            print*, 'koniec init'

            do i = 1, 3
                  block_pair(1, i) = 2
            end do
            do i = 4, 6
                  block_pair(1, i) = 4
            end do

            block_pair(1, 7) = 3
            block_pair(1, 8) = 3

            block_pair(2, 1) = 1
            block_pair(2, 2) = 5
            block_pair(2, 3) = 3

            block_pair(2, 4) = 1
            block_pair(2, 5) = 5
            block_pair(2, 6) = 3

            block_pair(2, 7) = 1
            block_pair(2, 8) = 5
            block_pair(1, 9) = 4
            block_pair(2, 9) = 5

            call clock_start(timer)
            ECorr_AC0 = zero
            do i = 1, 9
                  print*, 'teraz robie blok A1 sing', i, bbl_name(i)
                  il  = block_pair(1,i)
                  ir = block_pair(2,i)
                  call clock_start(timer0)
                  print*, bbl_name(i), eigsPA(il)%dim, eigsPA(ir)%dim
                  call calc_block(i, eigsPA(il), eigsPA(ir), n, Ha, TwoNO, TWONOA, Ind2, IndAux, IndMod, &
                        NBasis, NA, NI, NV, NInte1, NInte2, 0, Flags, E_contr_s(i))
                  print*, 'Czas na sing ', bbl_name(i), ' ', clock_readwall(timer0)
                  if (i <9)then
                        write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ')_s=', E_contr_s(i)
                  else
                        write(*,'(A8, A5, A3, F20.15)') 'E_contr(', bbl_name(i), ')_s=', E_contr_s(i)
                  end if
                  ECorr_AC0 = ECorr_AC0 + E_contr_s(i)
                  il = 5 + il
                  ir = 5 + ir
                  print*, 'teraz robie blok A1 trip', i, bbl_name(i)
                  !       stop
                  call clock_start(timer0)

                  call calc_block(i, eigsPA(il), eigsPA(ir), n, Ha, TwoNO, TWONOA, Ind2, IndAux, IndMod, &
                        NBasis, NA, NI, NV, NInte1, NInte2, 1, Flags, E_contr_t(i))
                  print*, 'Czas na trip ', bbl_name(i), ' ', clock_readwall(timer0)

                  if (i <9)then
                        write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ')_t=', E_contr_t(i)
                  else
                        write(*,'(A8, A5, A3, F20.15)') 'E_contr(', bbl_name(i), ')_t=', E_contr_t(i)
                  end if


                  if (i <9)then
                        write(*,'(A8, A4, A2, F20.15)') 'E_contr(', bbl_name(i), ')=', E_contr_s(i) + E_contr_t(i)
                  else
                        write(*,'(A8, A5, A2, F20.15)') 'E_contr(', bbl_name(i), ')=', E_contr_s(i) + E_contr_t(i)
                  end if
                  print*, i, bbl_name(i)
                  ECorr_AC0 = ECorr_AC0 + E_contr_t(i)
            end do
            print*, 'TIME na all blocks ',  clock_readwall(timer)

            print*, ''
            print*, 'RDSC ACPP0 CONTR ECORR ', ECorr_AC0
            print*, ''

      end subroutine ACPP0_fast


      subroutine eigs0_block(x, Eig, Eigvec, v_plus, n, Ha, TwoNOA, &
            IndN, Ind2, IndAux, IndMod, NDim, NBasis, NA, NI, NV, NInte1, &
            NInte2, spin_symm, Flags)

            integer, intent(in) :: x
            double precision, dimension(:), intent(inout) :: Eig
            double precision, dimension(:,:), intent(inout) :: Eigvec
            integer, dimension(:), intent(inout) :: v_plus
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),      intent(in) :: Ha
            double precision, dimension(:),      intent(in) :: TwoNOA
            integer, dimension(:,:), intent(in)  :: IndN
            integer, dimension(:), intent(in) :: Ind2
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NDim
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            integer, intent(in) :: spin_symm
            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:,:), allocatable :: MxA, MxS
            double precision, dimension(:), allocatable :: Eig_i
            type (tclock) :: timer, timer0
            integer :: i
            integer, parameter :: multiply_by_S = 1

            if (NDim .gt.0)then
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))
                  allocate(Eig_i(NDim))

                  !    call clock_start(timer)

                  print*, 'ten wymiar to', NDim
                  call clock_start(timer0)
                  call PPERPA_block(MxA, MxS, n, Ha, TwoNOA, IndN, IndN, Ind2, IndAux, IndMod, NDim, NDim,&
                        NBasis, NA, NI, NV, NInte1, NInte2, zero, spin_symm, Flags, multiply_by_S)
                  print*, 'Czas na block: ', clock_readwall(timer0)
                  print*, 'gwn'
                  ! print*, 'macA'
                  ! if (x==IAA)then
                  !    call geprn(MxA)
                  ! end if
                  ! print*, 'kmacA'


                  !  print*, 'Czas na blok va sing: '//str(clock_readwall(timer),d=2)
                  !   call clock_start(timer)
                  if (x == IOO.or.x == IVV.or. x== IOOT.or.x==IVVT)then
                        ! print*, ''
                        ! print*, 'tak', x, IVV
                        Eigvec = zero
                        do i = 1, NDim
                              Eig(i) = MxA(i, i)
                        end do
                        if (x == IOO.or.x==IOOT)then
                              v_plus = 0
                        else
                              v_plus = 1
                        end if
                  else
                        print*, 'nonsym'
                        call nonsymmetric_eigenproblem_block(Eig, Eig_i, Eigvec, MxA, MxS, v_plus, IndN, IndAux, 1)
                  end if

                  deallocate(MxA)
                  deallocate(MxS)
                  deallocate(Eig_i)
            end if
            ! print*, 'Czas diagonalizacji sing: '//str(clock_readwall(timer),d=2)


      end subroutine eigs0_block
      

      subroutine calc_block(x, ePAl, ePAr, n, Ha, TwoNO,  TWONOA, Ind2, IndAux, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, spin_symm, Flags, E_contr)

            integer, intent(in) :: x
            type(eigsBlockParams), intent(in) :: ePAl, ePAr
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),      intent(in) :: Ha
            double precision, dimension(:),      intent(in) :: TwoNOA, TWONO
            integer, dimension(:), intent(in) :: Ind2
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            integer, intent(in) :: spin_symm
            type(FlagsData), intent(in) :: Flags
            double precision, intent(out) :: E_contr
            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4                                    
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9
            integer :: p, q, r, s, i, j, pq
            integer :: np,km
            double precision :: Npqrs, Aux1, Aux2
            double precision, dimension(:,:), allocatable :: MxA1, MxA, MxS, ttt
            double precision, dimension(:), allocatable :: tempx, tempx2
            double precision :: gr, gr2, gr3, sa, tr1, tr2
            double precision, dimension(:,:), allocatable :: IMxA1, IMxA2
            integer :: z
            integer, external :: NAddrRDM
            integer, external :: NAddr3
            type (tclock) :: timer, timer0, timerall
            integer, parameter :: multiply_by_S = 0
            double precision :: val, temp, pluszek
            external :: dgemv

            E_contr = zero
            gr = zero
            gr2 = zero
            gr3 = zero
            sa = zero
            tr1 = zero
            allocate(MxA1(ePAl%dim, ePAr%dim))
            allocate(MxS(1,1))

            allocate(IMxA2(ePAr%dim, ePAl%dim))
            allocate(ttt(ePAl%dim, ePAr%dim))

            call clock_start(timer)
            call clock_start(timerall)
            ! call PPERPA_block(MxA, MxS, n, Ha, TwoNOA, ePAl%IndN, ePAr%IndN, Ind2, IndAux, IndMod, ePAl%dim, ePAr%dim, &
            !    NBasis, NA, NI, NV, NInte1, NInte2, zero, spin_symm, Flags, 0)
            ! print*, 'koniec bl 0'
            
            call PPERPA_block(MxA1, MxS, n, Ha, TwoNOA, ePAl%IndN, ePAr%IndN, Ind2, IndAux, IndMod, ePAl%dim, ePAr%dim, &
                  NBasis, NA, NI, NV, NInte1, NInte2, one, spin_symm, Flags, multiply_by_S)
            print*, 'koniec bl 1'


            ! MxA1 = MxA1- MxA

            print*, 'Czas na wybrany wycinek macierzy A1: ', clock_readwall(timer)
            deallocate(MxS)
            allocate(IMxA1(ePAl%dim, ePAr%dim))



            do i = 1, ePAl%dim
                  do j = 1, ePAr%dim
                        IMxA2(j, i) = MxA1(i, j)
                  end do
            end do
            !        p = ePAl%IndN(1, i)
            !          q = ePAl%IndN(2, i)

            !          r = ePAr%IndN(1, j)
            !          s = ePAr%IndN(2, j)

            !       write(*, '(6I5, F20.15)') p, q, r, s, i, j, MxA1(i, j)
            !    end do
            ! end do
            allocate(tempx(ePAl%dim))
            allocate(tempx2(ePAr%dim))



            if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then

                  call clock_start(timer)
                  IMxA1 = zero
                  do np = 1, ePAl%dim           
                        do km = 1, ePAr%dim
                              if (ePAl%v_plus(np)==1.and.ePAr%v_plus(km)==0)then

                                    call real_av_x(tempx, MxA1, ePAl%dim,  ePAr%Eigvec(:,km), ePAl%dim, ePAr%dim, one, zero)
                                    call real_vw_x(Aux2, ePAl%Eigvec(:,np), tempx, ePAl%dim)
                                    ! if (abs(Aux2).gt.1.d-5)then

                                    !       write(*, '(A15, 2I5, A30, 4F20.15)') 'dla danego', np, km, 'element Cnk jest rowny', Aux2, ePAl%Eig(np), ePAr%Eig(km), Aux2/(ePAl%Eig(np)- ePAr%Eig(km))
                                    ! end if

                                    Aux2 = Aux2 /(ePAl%Eig(np)-ePAr%Eig(km))

                                    !ttt(np, km) = Aux2
                                    ! if (abs(Aux2).gt.1.d-5)then
                                    !       write(*, '(A15, 2I5, A30, F20.15)') 'dla danego', np, km, 'element Cnk jest rowny', Aux2
                                    ! end if
                                    do j = 1, ePAr%dim
                                          IMxA1(np, j) = IMxA1(np, j) + Aux2 * ePAr%Eigvec(j, km)
                                          ! if (np==16.and.j==2)then
                                          !       write(*,'(3I5, 3F20.15)') np, j, km, Aux2 , ePAr%Eigvec(j, km), IMxA1(np, j)
                                          ! end if
                                    end do
                              end if
                        end do
                  end do
                  print*, 'Czas na pierwsze petle: ', clock_readwall(timer)

            end if


            call clock_start(timer)

            select case(x)

            case(VVAO, VVAA)
                  !    if (x==VVAO.or.x==VVAA)then
                  call clock_start(timer)
                  do i = 1, ePAl%dim
                        do km = 1, ePAr%dim
                              if (ePAr%v_plus(km)==0)then
                                    call real_vw_x(Aux2, IMxA2(:, i), ePAr%Eigvec(:,km), ePAr%dim)
                                    do j = 1, ePAr%dim
                                          p = ePAl%IndN(1, i)
                                          q = ePAl%IndN(2, i)

                                          r = ePAr%IndN(1, j)
                                          s = ePAr%IndN(2, j)

                                          if (spin_symm == 0)then
                                                Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
                                                if (p==q)then
                                                      Aux1 = Aux1 * sqrt(frac12)
                                                end if
                                                if (r==s) then
                                                      Aux1 = Aux1 * sqrt(frac12)
                                                end if
                                          else
                                                Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
                                                Aux1 = Three*Aux1
                                          end if

                                          Npqrs =(one-n(p)-n(q)) * (one-n(r)-n(s))

                                          E_contr = E_contr - Aux2 * Npqrs * Aux1 / (ePAl%Eig(i)-ePAr%Eig(km)) *ePAr%Eigvec(j, km)

                                          !write(*, '(7I5, 6F20.15)') i, km, j, p, q, r, s, Aux2, Npqrs, Aux1, TwoNO(NAddr3(r,p,s,q)), (ePAl%Eig(i)-ePAr%Eig(km)) *ePAr%Eigvec(j, km), E_contr

                                    end do
                              end if
                        end do
                  end do
                  print*, 'Czas na drugie petle: ', clock_readwall(timer)

            case(VAOO, AAOO)

                  call clock_start(timer)

                  do j = 1, ePAr%dim
                        r = ePAr%IndN(1, j)
                        s = ePAr%IndN(2, j)
                        do np = 1, ePAl%dim
                              if (ePAl%v_plus(np)==1)then

                                    ! temp = 0
                                    ! do pq = 1, ePAl%dim 
                                    !       temp = temp + MxA1(pq, j)* ePAl%Eigvec(pq,np)
                                    !       if (abs(MxA1(pq, j)* ePAl%Eigvec(pq,np)).gt.1.d-5)then
                                    !             write(*, '(A10, 5I5, 2F20.15)') 'debil', j, r, s, np, pq, MxA1(pq, j), ePAl%Eigvec(pq,np)
                                    !       end if
                                    ! end do


                                    call real_vw_x(Aux2, MxA1(:, j), ePAl%Eigvec(:,np), ePAl%dim)
                                    !Aux2 = Aux2 / (ePAl%Eig(np)-ePAr%Eig(j))
                                    ! if (abs(Aux2).gt.1.d-5)then
                                    !       write(*,'(A4, 2I5, 2F20.15)') 'jrs', np, j, Aux2, ePAl%Eig(np)
                                    ! end if

                                    ! if (abs(Aux2/ (ePAl%Eig(np)-ePAr%Eig(j))).gt.1.d-5)then
                                    !       write(*, '(A10, 2I5, 4F20.15)') 'niebieski', j, np, Aux2 , ePAl%Eig(np),ePAr%Eig(j)
                                    ! end if

                                    do i = 1, ePAl%dim
                                          p = ePAl%IndN(1, i)
                                          q = ePAl%IndN(2, i)
                                          if (spin_symm == 0)then
                                                if (x == VAAO)then
                                                      Aux1 = TwoNO(NAddr3(r,p,s,q))
                                                else if (x==VAAO2) then
                                                      Aux1 = TwoNO(NAddr3(r,q,s,p))
                                                else
                                                      Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
                                                end if
                                                if (p==q)then
                                                      Aux1 = Aux1 * sqrt(frac12)
                                                end if
                                                if (r==s) then
                                                      Aux1 = Aux1 * sqrt(frac12)
                                                end if
                                          else
                                                if (x == VAAO)then
                                                      Aux1 = TwoNO(NAddr3(r,p,s,q))
                                                else if (x==VAAO2) then
                                                      Aux1 = -TwoNO(NAddr3(r,q,s,p))
                                                else
                                                      Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
                                                end if
                                                Aux1 = Three*Aux1
                                          end if

                                          Npqrs =(one-n(p)-n(q)) * (one-n(r)-n(s))

                                          !write(*, '(A5, 4I5, 6F20.15)')'calk', r, s, p, q, Npqrs, TwoNO(NAddr3(r,p,s,q)) ,  TwoNO(NAddr3(r,q,s,p)), Aux1


                                          E_contr = E_contr - Aux2 * Npqrs * Aux1 / (ePAl%Eig(np)-ePAr%Eig(j)) *ePAl%Eigvec(i,np)
                                          !E_contr = E_contr - Aux2 * Npqrs * Aux1 *ePAl%Eigvec(i,np)

                                          !                                           if (abs(Aux2 * Npqrs * Aux1 *ePAl%Eigvec(i, np)).gt.1.d-5)then
                                          !                             write(*, '(A10, 4I5, F20.15, 3I5, 5F20.15)') 'niebieski', p, q, r, s,  (one-n(p)-n(q)) * (one-n(r)-n(s)), j, i, np, Aux2, ePAl%Eigvec(i,np), E_contr, TwoNO(NAddr3(r,p,s,q)) ,  TwoNO(NAddr3(r,q,s,p))
                                          ! !                                                write(*, '(A7, 7I5, 6F20.15)') 'contre', p, q, r, s, i, j, np, Aux2 , Npqrs * Aux1 ,ePAl%Eigvec(i, np), E_contr
                                          !                                           end if
                                          ! if (abs(Aux1*Npqrs).gt.1.d-5)then
                                          !       if  (np==1)then
                                          !             write(*, '(A10, 4I, 4F20.15)') 'kurwa', p, q, r, s,	Aux1 * Npqrs, Npqrs, TwoNO(NAddr3(r,p,s,q)) , TwoNO(NAddr3(r,q,s,p))
                                          !       end if
                                          ! end if

                                    end do
                              end if
                        end do
                  end do

                  print*, 'Czas na trzecie petle: ', clock_readwall(timer)

            case default 

                  call clock_start(timer)    
                  !    !$omp parallel do collapse(2)&
                  !    !$omp default(shared) &
                  !    !$omp prIVATe(i, j) &
                  !    !$omp prIVATe(p, q, Npqrs, r, s, Aux1, Aux2, E_contr, np, km)
                  i_rowloops: do i = 1, ePAl%dim
                        j_colloops: do j = 1, ePAr%dim
                              p = ePAl%IndN(1, i)
                              q = ePAl%IndN(2, i)

                              r = ePAr%IndN(1, j)
                              s = ePAr%IndN(2, j)

                              if (spin_symm == 0)then
                                    if (x == VAAO)then
                                          Aux1 = TwoNO(NAddr3(r,p,s,q))
                                    else if (x==VAAO2) then
                                          Aux1 = TwoNO(NAddr3(r,q,s,p))
                                    else
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
                                          ! if (abs(TwoNO(NAddr3(r,p,s,q))).gt.1.d-5)then
                                          !       write(*, '(A10, 4I5, 2F20.15)')'calk', p, q, r, s, TwoNO(NAddr3(r,p,s,q)) ,  TwoNO(NAddr3(r,q,s,p))
                                          ! end if
                                    end if
                                    if (p==q)then
                                          Aux1 = Aux1 * sqrt(frac12)
                                    end if
                                    if (r==s) then
                                          Aux1 = Aux1 * sqrt(frac12)
                                    end if
                              else
                                    if (x == VAAO)then
                                          Aux1 = TwoNO(NAddr3(r,p,s,q))
                                    else if	(x==VAAO2) then
                                          Aux1 = -TwoNO(NAddr3(r,q,s,p))
                                    else                
                                          Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
                                    end if
                                    Aux1 = Three*Aux1
                              end if

                              Npqrs =(one-n(p)-n(q)) * (one-n(r)-n(s))

                              if (x==VVOO)then
                                    E_contr = E_contr - MxA1(i, j) * Npqrs * Aux1 / (ePAl%Eig(i)-ePAr%Eig(j))
                              else if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then
                                    do np = 1, ePAl%dim
                                          E_contr = E_contr - Npqrs * Aux1 * ePAl%Eigvec(i,np) * IMxA1(np, j)
                                    end do
                              end if
                              
                                    !                                    write(*, '(A10, 4I5, 5F20.15)') 'witam', p, q, r, s, Npqrs, Aux1, MxA1(i, j), ePAl%Eig(i), ePAr%Eig(j)

                              ! else if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then                                    
                              !       do np = 1, ePAl%dim
                              !             do km = 1, ePAr%dim
                              !                   if (ePAl%v_plus(np)==1.and.ePAr%v_plus(km)==0)then
                              !                         !E_contr = E_contr - Npqrs * Aux1 * ePAl%Eigvec(i,np) * IMxA1(np, j)
                              !                         E_contr = E_contr - Npqrs * Aux1 * ePAl%Eigvec(i,np) * ePAr%Eigvec(j, km) * ttt(np, km)



                              !                         ! if (abs(Npqrs * Aux1 * ePAl%Eigvec(i,np) * ePAr%Eigvec(j, km) * ttt(np, km)).gt.1.d-5)then
                              !                         !       write(*, '(A10, 4I5, F15.7, 4I5, 8F15.7)') 'niebieski', p, q, r, s,  (one-n(p)-n(q)) * (one-n(r)-n(s)), &
                              !                         !             j, i, np, km,  Npqrs * Aux1 , ePAl%Eigvec(i,np) , ttt(np, km), ePAr%Eigvec(j, km)

                              !                         ! end if
                              !                         ! if(p==11.and.q==8.and.r==9.and.s==3)then
                              !                         !       if (abs(ePAl%Eigvec(i,np) ).gt.1.d-5.and.abs(ePAr%Eigvec(j,km)).gt.1.d-5)then
                              !                         !             write(*, '(A10, 4I5, F15.7, 4I5, 8F15.7)') 'niebieskl', p, q, r, s,  (one-n(p)-n(q)) * (one-n(r)-n(s)), &
                              !                         !                   j, i, np, km,  Npqrs * Aux1 , ePAl%Eigvec(i,np) , ttt(np, km), ePAr%Eigvec(j, km), &
                              !                         !                   Npqrs * Aux1 * ePAl%Eigvec(i,np) * ePAr%Eigvec(j, km) * ttt(np, km)
                              !                         !       end if
                              !                         ! end if
                              !                   end if
                              !             end do
                              !       end do
                              ! end if

                        end do j_colloops
                  end do i_rowloops
                  print*, 'Czas na czwarte petle: ', clock_readwall(timer)

            end select
            !       !omp end parallel do

            print*, ''
            print*, 'Czas na timeall  petle: ', clock_readwall(timerall)
            print*, ''

            deallocate(MxA1)
            deallocate(tempx)


      end subroutine calc_block




      subroutine updateIndices(p, q, IndN_s, IndN_t, ind_s, ind_t)
            integer, intent(in) :: p, q
            integer, dimension(:,:), intent(inout) :: IndN_s, IndN_t
            integer, intent(inout) :: ind_s, ind_t

            if (p.ne.q) then
                  IndN_t(1,ind_t) = p
                  IndN_t(2,ind_t) = q
                  ind_t = ind_t + 1
            endif
            IndN_s(1,ind_s) = p
            IndN_s(2,ind_s) = q
            ind_s = ind_s + 1

      end subroutine updateIndices


      subroutine nonsymmetric_eigenproblem_block(wr, wi, vr, A, S, v_plus, IndN, IndAux, AC_TYPE)
            double precision, dimension(:), intent(out)      :: wr
            double precision, dimension(:), intent(out)      :: wi
            double precision, dimension(:, :), allocatable   :: vl
            double precision, dimension(:, :), intent(out)   :: vr
            double precision, dimension(:, :), intent(inout) :: A, S
            integer, dimension(:), intent(out) :: v_plus
            integer, dimension(:,:), intent(in) :: indN
            integer, dimension(:),intent(in) :: IndAux
            integer, intent(in) :: AC_TYPE
            double precision :: shift, max_hh, min_pp
            double precision, dimension(:),allocatable :: S_diag

            double precision, dimension(1) :: work0
            double precision, dimension(:), allocatable :: work, tempx
            integer :: lwork, info
            integer :: n, i, j, k
            double precision :: norm, norm_l, norm_rl
            double precision :: dd, ddm, dd_comp
            type (tclock) :: timer
            double precision, parameter :: tol = 1.d-4
            double precision, dimension(:), allocatable :: wr_plus, wr_minus

            integer, dimension(:), allocatable :: dy_plus, dy_minus
            integer :: count, countj_plus, countj_minus, ss
            logical :: flag_change

            double precision, dimension(:,:), allocatable :: miniA
            double precision, dimension(:,:), allocatable :: vl2, vr2
            double precision, dimension(:), allocatable :: wr2, wi2
            integer :: ii, jj, n2, r, s2, p,q


            external :: dgeev

            n = size(A, dim=1)
            allocate(dy_plus(n))
            allocate(dy_minus(n))
            allocate(wr_plus(n))
            allocate(wr_minus(n))
            allocate(S_diag(n))

            allocate(vl(1,1))
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
            print*, 'info', info
            if (info /= 0) then
                  print*, "Nonsymmetric matrix eigendecompositino failed with info="
                  error stop
            end if
            print*, 'Czas diagonalizacji real: ', clock_readwall(timer)

            allocate(tempx(n))
            call clock_start(timer)
            v_plus = 0


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

            ! if (AC_TYPE == PPAC)then
            !       print*, 'AC_TYPE: PP'
            ! else if (AC_TYPE == HHAC) then
            !       print*, 'AC_TYPE: HH'
            ! end if


            call clock_start(timer)
            print*, 'eigenvalues of singlets'
            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'COMPLEX EIGENVALUES', i, wr(i), wi(i)
                        !         stop
                  end if
            end do

            wr_plus = wr
            wr_minus = wr
            do i = 1, n
                  dy_plus(i)=i
            end do
            call dsort(wr_plus, dy_plus, n)


            shift = zero
            do i = 1, n
                  if (v_plus(dy_plus(i))== 0 )then
                        max_hh = wr_plus(i)
                  end if
                  if (v_plus(dy_plus(i))== 1 )then
                        min_pp = wr_plus(i)
                        if (i == 1) then
                              shift = min_pp/two
                        else
                              shift = (abs(min_pp) + abs(max_hh))/two
                              print*, 'max_hh', max_hh
                        end if
                        print*, 'min_pp' , wr_plus(i)
                        print*, 'shift', shift
                        exit
                  else

                  end if
            end do


            wr_plus = zero
            wr_minus = zero
            dy_plus = 0
            dy_minus = 0
            countj_plus=1
            countj_minus=1

            do i = 1, n
                  if (v_plus(i) == 1) then
                        wr_plus(countj_plus) = wr(i)
                        dy_plus(countj_plus) = i
                        countj_plus = countj_plus + 1
                  else if (v_plus(i) == 0) then
                        wr_minus(countj_minus) = wr(i)
                        dy_minus(countj_minus) = i
                        countj_minus = countj_minus + 1

                  end if
            end do
            countj_plus = countj_plus-1
            countj_minus = countj_minus-1

            call dsort(wr_plus(1:countj_plus), dy_plus(1:countj_plus), countj_plus)
            call dsort(wr_minus(1:countj_minus), dy_minus(1:countj_minus), countj_minus)

            call orthogonalize_degen(n, countj_plus, wr_plus, vr, S_diag, dy_plus, 1)
            call orthogonalize_degen(n, countj_minus, wr_minus, vr, S_diag, dy_minus, 0)

      end subroutine nonsymmetric_eigenproblem_block



      subroutine PPERPA_block(MxA, MxS, n, Ha, TwoNOA, IndN1, IndN2, Ind2, IndAux, IndMod, NDim1,NDim2, &
            NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, spin_symm, Flags, multiply_by_S)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            double precision, dimension(:,:), intent(inout) :: MxA, MxS
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),      intent(in) :: Ha
            double precision, dimension(:),      intent(in) :: TwoNOA
            integer, dimension(:,:), intent(in)  :: IndN1, IndN2
            integer, dimension(:), intent(in) :: Ind2
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NDim1, NDim2, NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            integer, intent(in) :: spin_symm
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), allocatable :: MxT

            double precision, dimension(:), allocatable :: R00, R11
            integer :: NRDM2, NRDM2Act, NOc
            double precision :: Arspq, Arsqp, Brspq, Crspq, dm, Arssum
            type (tclock) :: timer

            !
            ! External procedures
            !    
            integer, external :: NAddrRDM
            integer, external :: NAddr3
            double precision, external :: FRDM2

            integer :: p, q, r, s, t, u, v, pp, qq
            integer :: i ,j, k, l, ij, a, b, ab, kl
            integer :: ii, pq, rs, pr, qs, ps, qr
            double precision :: num_f, num, num_h
            integer :: c, d
            double precision :: temp, temp0, val
            integer :: jj, tt, uu, vv
            integer, dimension(11) :: lst
            logical :: ism
            print*, 'multiply_by_S', multiply_by_S

            lst(1) = 3
            lst(2)=5
            lst(3)=12
            lst(4)=13
            lst(5)=18
            lst(6)=21
            lst(7)=23
            lst(8)=24
            lst(9)=26
            lst(10)=28
            lst(11)=29



            NRDM2 = NBasis**2*(NBasis**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NOc = NI + NA

            allocate (R00(NRDM2Act))
            allocate (R11(NRDM2Act))

            R00 = Zero
            R11 = Zero

            MxA = zero
            MxS = zero

            call read_2rdm("rdm2.dat", R00, NA)
            call read_2rdm("rdms2.dat", R11, NA)

            print*, 'multiply_by_s', multiply_by_s
            call clock_start(timer) 
                       !$omp parallel do collapse(2)&
                       !$omp default(shared) &
                       !$omp private(i, j, tt, uu, vv, t, u, v) &
                       !$omp private(p, q, r, s, qs, pr, qr, ps, Arspq, Arsqp, Arssum, num_f)
            j_colloops: do j = 1, NDim2
                  i_rowloops: do i = 1, NDim1


                        r = IndN1(1, i)
                        s = IndN1(2, i)


                        p = IndN2(1, j)
                        q = IndN2(2, j)


                        Arspq = zero
                        Arsqp = zero
                        qs = (max(s, q)*(max(s, q)-1))/2 + min(s, q)
                        pr = (max(r, p)*(max(r, p)-1))/2 + min(r, p)

                        qr = (max(r, q)*(max(r, q)-1))/2 + min(r, q)
                        ps = (max(s, p)*(max(s, p)-1))/2 + min(s, p)




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
                              do t = 1, NI+NA
                                    Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(p,r,t,t))-TwoNOA(NAddr3(p,t,t,r)))
                              end do
                        end if

                        !                 !P3a-1
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


                        !                  ! !P3c-1,2                                                                           
                        !                  ! !P3c-1,2                                                                                                                                                                                                                                                            
                        if (IndAux(s)==(1).and.IndAux(q)==(1))then
                              do tt = NI+1, NI+NA
                                    do uu = NI+1, NI+NA
                                          t = IndMod(tt)
                                          u = IndMod(uu)
                                          Arspq = Arspq -  TwoNOA(NAddr3(p,r,t,u))* &
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
                              do tt = NI+1, NI+NA
                                    do uu = NI+1, NI+NA
                                          t = IndMod(tt)
                                          u = IndMod(uu)
                                          Arspq = Arspq -  TwoNOA(Naddr3(q,s,t,u))* &
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

                        !P4a-2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                        if (p==r)then
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

                        !P4b-1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
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
                              end if
                        end if
                        !                 !P5a                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                        if (p==r)then
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


                        !                 !P5b                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
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

                        ! ! ! Arsqp                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

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
                        if (multiply_by_s==1)then
                              num_f = one/(1-n(r)-n(s))
                        end if

                        select case(spin_symm)
                        case(0)
                              num_h = one
                              if (p==q.and.r==s)then
                                    num_h = frac12
                              else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
                                    num_h = sqrt(frac12)
                              end if

                              Arssum = Arspq  - Arsqp

                              !                   !$omp critical
                              ! if (abs(Arssum).gt.1.d-5)then
                              !    print*, multiply_by_s, num_f, num_h, Arspq, Arsqp
                              ! end if

                              MxA(i, j) = num_f * num_h *  Arssum
                              !                   !$omp end critical

                              ! if (abs(MxA(i, j)).gt.1.d-5)then
                              !    write(*,'(2I5, A5, 4I5, F20.15)') i, j,       '  |  ', p, q, r, s,MxA(i, j)
                              ! end if

                        case(1)
                              if (p.ne.q.and.r.ne.s)then
                                    Arssum = Arspq + Arsqp
                                    !                      !$omp critical
                                    MxA(i, j) = num_f * Arssum

                                    ! if (abs(MxA(i, j)).gt.1.d-5)then
                                    !       write(*,'(2I5, A5, 4I5, F30.16)') i, j,       '  |  ', p, q, r, s,MxA(i, j)
                                    ! end if

                                    !                      !$omp end critical
                              end if
                        end select

                  end do i_rowloops
            end do j_colloops

                !$omp end parallel do                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

            print*, ''
            print*, 'Time for construction of A: '//str(clock_readwall(timer),d=2)





!!!!!---------------------powyzej kod z acpp

            !     call clock_start(timer)        
            !     ! print*, 'pozdro z PPERPA_block', spin_symm
            !   !$omp parallel do collapse(2)&
            !   !$omp default(shared) &
            !   !$omp prIVATe(i, j) &
            !   !$omp prIVATe(p, q, pq, r, s, rs, qs, pr, qr, ps, Arspq, Arsqp, Arssum, t, num_f)
            !     i_rowloops: do i = 1, NDim1
            !        j_colloops: do j = 1, NDim2

            !           r = IndN1(1, i)
            !           s = IndN1(2, i)

            !           p = IndN2(1, j)
            !           q = IndN2(2, j)

            !           Arspq = zero
            !           Arsqp = zero
            !           qs = (max(s, q)*(max(s, q)-1))/2 + min(s, q)
            !           pr = (max(r, p)*(max(r, p)-1))/2 + min(r, p)

            !           qr = (max(r, q)*(max(r, q)-1))/2 + min(r, q)
            !           ps = (max(s, p)*(max(s, p)-1))/2 + min(s, p)

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'pocz',  Arspq
            !           ! end if

            !           Arspq = Arspq + TwoNOA(NAddr3(p,r,q,s))* (one - n(p)-n(q)-n(r)-n(s))

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*, TwoNOA(NAddr3(p,r,q,s)), n(p), n(r), n(q), n(s)
            !           !    print*,'calk',  Arspq
            !           ! end if

            !           if (p==r)then
            !              Arspq = Arspq + Ha(qs) * (one-n(p)-frac12*n(q)-frac12*n(s))
            !           end if

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'ha-qs',  Arspq
            !           ! end if

            !           if (q==s)then
            !              Arspq = Arspq + Ha(pr) * (one-n(q)-frac12*n(p)-frac12*n(r))
            !           end if

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'ha-pr',  Arspq, Ha(pr)
            !           ! end if

            !           if (p==r)then   
            !              do t = 1, NI+NA
            !                 Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(q,s,t,t))-TwoNOA(NAddr3(q,t,t,s)))
            !              end do
            !           end if
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'nt1',  Arspq
            !           ! end if


            !           if (q==s)then
            !              temp = zero
            !              do t = 1, NI+NA
            !                 Arspq = Arspq + n(t) * (two* TwoNOA(NAddr3(p,r,t,t))-TwoNOA(NAddr3(p,t,t,r)))
            !              end do
            !           end if
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'nt2',  Arspq
            !           ! end if


            !           Arspq = Arspq + P3a(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq + P3a(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P3a',  Arspq
            !           ! end if


            !           Arspq = Arspq + P3b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq + P3b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P3b',  Arspq
            !           ! end if


            !           Arspq = Arspq - P3c(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq - P3c(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P3c',  Arspq
            !           ! end if


            !           Arspq = Arspq - P4a(q, s, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq - P4a(p, r, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P4a',  Arspq
            !           ! end if


            !           Arspq = Arspq - P4b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq - P4b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P4b',  Arspq
            !           ! end if


            !           Arspq = Arspq + P5a(q, s, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq + P5a(p, r, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P5a',  Arspq
            !           ! end if

            !           arspq = Arspq + P5b(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arspq = Arspq + P5b(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           ! if (r==12.and.s==5.and.p==7.and.q==5)then
            !           !    print*,'P5b',  Arspq
            !           ! end if


            !           Arsqp = Arsqp - TwoNOA(NAddr3(q,r,p,s))* (one - n(q)-n(p)-n(r)-n(s))

            !           if (q==r)then
            !              Arsqp = Arsqp - Ha(ps) * (one-n(q)-frac12*n(p)-frac12*n(s))
            !           end if

            !           if (p==s)then
            !              Arsqp = Arsqp - Ha(qr) * (one-n(p)-frac12*n(q)-frac12*n(r))
            !           end if

            !           if (q==r)then
            !              do t = 1, NBasis
            !                 Arsqp = Arsqp - n(t) * (two*TwoNOA(NAddr3(p,s,t,t))-TwoNOA(NAddr3(p,t,t,s)))
            !              end do
            !           end if

            !           if (p==s)then
            !              do t = 1, NBasis
            !                 Arsqp = Arsqp - n(t) * (two*TwoNOA(NAddr3(q,r,t,t))-TwoNOA(NAddr3(q,t,t,r)))
            !              end do
            !           end if


            !           Arsqp = Arsqp - P3a(s, q, p, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp - P3a(r, p, q, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)


            !           Arsqp = Arsqp - P3b(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp - P3b(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           Arsqp = Arsqp + P3c(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp + P3c(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           Arsqp = Arsqp + P4a(p, s, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp + P4a(q, r, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)


            !           Arsqp = Arsqp + P4b(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp + P4b(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)


            !           Arsqp = Arsqp - P5a(p, s, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp - P5a(q, r, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           Arsqp = Arsqp - P5b(s, p, q, r, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)
            !           Arsqp = Arsqp - P5b(r, q, p, s, TwoNOA, R00, R11, n, IndAux, IndMod, Ind2, NI, NA, NBasis, NA)

            !           num_f = one

            !           select case(spin_symm)
            !           case(0)
            ! !               if (p==1.and.q==1.and.r==1.and.s==1)then
            !   !              print*, 'jebane gowno'
            !  !            end if

            !              if (p==q.and.r==s)then
            !                 num_f = frac12
            !              else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
            !                 num_f = sqrt(frac12)
            !              end if

            !              Arssum = Arspq  - Arsqp
            !              MxA(i, j) = num_f * Arssum

            !              ! if (abs(MxA(i, j)).gt.1.d-6)then
            !              !    ism = ismixed(lst, p,q,r,s)
            !              !    if (ism) then
            !              !       !                if (s==3.or.s==5.or.p==3.or.p==5.or.q==5.or.q==3)then
            !              !       ! if (IndAux(p)==0.and.IndAux(q)==0.and.IndAux(r)==0.and.IndAux(s)==0)then
            !              !       write(*, '(A10, 6I5, F20.10)') 'mxamxa', i, j, r, s, p,q, MxA(i, j)
            !              !    end if
            !              ! end if

            !           case(1)
            ! !             if (p==1.and.q==1.and.r==1.and.s==1)then
            !    !             print*, 'tak robie triple'
            !  !            end if
            !              if (p.ne.q.and.r.ne.s)then
            !                 Arssum = Arspq + Arsqp
            !                 MxA(i, j) = Arssum 
            !              end if
            !           end select

            !        end do j_colloops
            !     end do i_rowloops
            !     !$omp end parallel do

            print*, 'Czas tworzenie macierzy A w danym bloku : ', clock_readwall(timer)
            if (multiply_by_s==1)then
                  call clock_start(timer)
                  MxS = zero
                  do i = 1, NDim1          
                        p = IndN1(1, i)
                        q = IndN1(2, i)
                        !          MxS(i,i) = one / (1-n(p)-n(q))
                        MxS(i,i) = (1-n(p)-n(q))
                        !          print*, 'mxs', mxS(i,i), (1-n(p)-n(q))
                  end do
                  !        allocate(MxT(NDim1, NDim1))
                  !        MxT=zero

                  !        call dgemm("N", "N", NDim1, NDim1, NDim1, 1.d+0, MxS, NDim1, MxA, NDim1, zero, MxT, NDim1)
                  !        MxA = MxT

                  ! do i = 1, NDim1
                  !    do j = 1, NDim2
                  !       r = IndN1(1, i)
                  !       s = IndN1(2, i)


                  !       p = IndN2(1, j)
                  !       q = IndN2(2, j)
                  !       if (abs(MxA(i, j)).gt.1.d-5)then                                                                                                                                                                                                                                     
                  !          write(*,'(2I5, A5, 4I5, F20.15)') i, j,       '  |  ', r, s, i, j, MxA(i, j)                                                                                                                                                                                      
                  !       end if
                  !    END do
                  ! end do

                  !           do i = 1, NDim1
                  !              p = IndN1(1, i)
                  !              q = IndN1(2, i)

                  ! !          MxS(i,i) = one / MxS(i,i)
                  !            MxS(i,i) = (1-n(p)-n(q))
                  !        end do
                  print*, ''
                  print*, 'Time for reszta: ', clock_readwall(timer)


                  ! deallocate(MxT)
            end if


      end subroutine PPERPA_block


      function ismixed(lst, p,q,r,s)
            logical :: ismixed
            integer, dimension(9), intent(in) :: lst
            integer, intent(in) :: p, q, r, s
            logical :: pf, qf, rf, sf

            pf = any(lst == p)
            qf = any(lst == q)
            rf = any(lst == r)
            sf = any(lst == s)

            ! if( (pf.and.qf.and.rf.and.sf).or. (.not.pf.and..not.qf.and..not.rf.and..not.sf))then
            !    ismixed = .false.
            ! else
            !    ismixed = .true.
            ! end if

            if ((pf .and. qf .and. rf .and. .not. sf) .or. &
                  (pf .and. qf .and. .not. rf .and. sf) .or. &
                  (pf .and. .not. qf .and. rf .and. sf) .or. &
                  (.not. pf .and. qf .and. rf .and. sf)) then
                  ismixed = .true.
                  ! Sprawdzenie, czy trzy są false, a jeden true
            else if ((.not. pf .and. .not. qf .and. .not. rf .and. sf) .or. &
                  (.not. pf .and. .not. qf .and. rf .and. .not. sf) .or. &
                  (.not. pf .and. qf .and. .not. rf .and. .not. sf) .or. &
                  (pf .and. .not. qf .and. .not. rf .and. .not. sf)) then
                  ismixed = .true.
            else if ((pf .and. qf .and. .not. rf .and. .not. sf) .or. &
                  (pf .and. .not. qf .and. rf .and. .not. sf) .or. &
                  (pf .and. .not. qf .and. .not. rf .and. sf) .or. &
                  (.not. pf .and. qf .and. rf .and. .not. sf) .or. &
                  (.not. pf .and. qf .and. .not. rf .and. sf) .or. &
                  (.not. pf .and. .not. qf .and. rf .and. sf)) then
                  ismixed = .true.
            else
                  ismixed = .false.
            end if

      end function ismixed



      subroutine PPERPA_init(H_type, ETot, ENuc,  n, XOne, TwoNO, Ind2, IndAux, &
            NBasis, NA, NI, NV, NInte1, NInte2,ACAlpha, Flags, TwoNOA, Ha)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: H_type
            double precision, intent(inout) :: ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            integer, dimension(:), intent(in)    :: IndAux
            integer, dimension(:), intent(out)    :: Ind2
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:), intent(out) :: TwoNOA, Ha
            double precision, dimension(:), allocatable :: HNO
            double precision, dimension(:), allocatable :: R00, R11
            integer :: NRDM2, NRDM2Act, NOc
            type (tclock) :: timer

            !
            ! External procedures
            !    
            integer, external :: NAddrRDM
            integer, external :: NAddr3
            double precision, external :: FRDM2

            integer :: p, q, r, s, t, u, v, ii
            integer :: i ,j, k, l, ij, a, b, ab, kl

            integer :: twoint_dim
            double precision :: temp, temp3, ecp, etot1, pppp, pmpm

            NRDM2 = NBasis**2*(NBasis**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NOc = NI + NA

            allocate(HNO(NInte1))
            allocate (R00(NRDM2Act))
            allocate (R11(NRDM2Act))

            ! twoint_dim = size(TwoNO, dim=1)
            ! allocate(TwoNOA(twoint_dim))
            print*, 'allocate(R00(NRDM2Act))', NRDM2Act

            R00 = Zero
            R11 = Zero
            print*, 'lllll', NRDM2Act

            call clock_start(timer)
            !                                                                                                                                      
            ! Fill one-electron hamiltonian in NO rep.                                                                                                
            !
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                        HNO(ij) = XOne(ab)
                        ! HNO(ij) = zero
                        ! do a = 1, NBasis
                        !    do b = 1, Nbasis
                        !       ab = (max(a, b)*(max(a,b)-1))/2 + min(a, b)
                        !       HNO(ij) = HNO(ij) + URe(i, a)*URe(j,b)*XOne(ab)
                        !    end do
                        ! end do
                        !                        if (i.ne.j)then
                        if (abs(hno(ij)).gt.1.d-5)then
                              print*, 'kupa', i, j, hno(ij)
                        end if
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
                        ! if (i.ne.j)then
                        !       if (abs(temp).gt.1.d-5)then
                        !             print*, 'kupa2', i, j, temp
                        !       end if
                        ! end if

                  end do
            end do

            print*, 'Time for fill-one-el: ', clock_readwall(timer)
            call clock_start(timer)


            print*, 'Number of inacvive orbitals:', NI
            print*, 'Number of active orbitals:  ', NA
            print*, 'Number of virtual orbitals: ', NV


            ind2 = 0

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

            R00 = Zero
            R11 = Zero
            ! print*, R00
            !     print*, 'czytam rdm2'
            call read_2rdm("rdm2.dat", R00, NA)
            call read_2rdm("rdms2.dat", R11, NA)
            
            print*, 'rdm++++'
            do p = NI+1, NI+NA
                  do q = NI+1, NI+NA
                        do r = NI+1, NI+NA
                              do s = NI+1, NI+NA
                                    pppp = frac12 * (R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)) &
                                          + R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)))
                                    if (abs(pppp).gt.1.d-5)then
                                          write(*, '(4I5, F20.15)') p, q, r, s, pppp
                                    end if
                              end do
                        end do
                  end do
            end do
            print*, ''
            print*, 'rdm+-+-'
            do p = NI+1, NI+NA
                  do q = NI+1, NI+NA
                        do r = NI+1, NI+NA
                              do s = NI+1, NI+NA
                                    pmpm = frac12 * (R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)) &
                                          - R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)))
                                    if (abs(pmpm).gt.1.d-5)then
                                          write(*, '(4I5, F20.15)') p, q, r, s, pmpm
                                    end if

                              end do
                        end do
                  end do
            end do

            print*, ''
            
            ! print*,'przecz', R00
            ETot1 = zero
            do i = 1, NBasis
                  ii = (i*(i+1))/2
                  ETot1 = ETot1 + two* n(i) * HNO(ii)
                  print*, i, HNO(ii), two* n(i) 
            end do
            
            ecp = etot1
            print*, 'etot1', ETot1, NOc
            print*, 'NOc', NOc

            !    print*, 'ind2'
            !    print*, ind2
            ! print*, 'this is rdm', Ind2
            ! do p = NI+1, NI+NA
            !       do q = NI+1, NI+NA
            !             do r = NI+1, NI+NA
            !                   do s = NI+1, NI+NA
            !                         if (abs(FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis)).gt.1.d-5)then
            !                               if (p.ne.q.and.p.ne.r.and.p.ne.s.and.q.ne.r.and.q.ne.s.and.r.ne.s)then
            !                                     write(*, '(4I5, 2F20.15)')p, q, r, s,  two * FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis)
            !                               end if
            !                         end if
            !                   end do
            !             end do
            !       end do
            ! end do
            ! stop

            etot = zero
            do p = 1, NI+NA
                  do q = 1, NI+NA
                        do r = 1, NI+NA
                              do s = 1, NI+NA
                                    ETot = ETot + FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis) &
                                          * TwoNO(NAddr3(p, r, q, s))
                                    if (abs(TwoNO(NAddr3(p, r, q, s))* FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis)).gt.1.d-5)then
                                          write(*, '(4I5, 3F20.15)')p, q, r, s, &
                                                TwoNO(NAddr3(p, r, q, s)), FRDM2(p, q, r, s, R00, n, Ind2, NA, NBasis), etot
                                    end if

                              end do
                        end do
                  end do
            end do

            etot = etot + etot1
            print*, 'testing...', ETot, ETot+Enuc, Enuc
            print*, 'RDCS ETot', ETot+Enuc
            print*, 'Time for testing etot: ', clock_readwall(timer)
            call clock_start(timer)
!            stop

            Ha = zero
            if (H_type==H_DYALL)then
                  print*, 'H_type=Dyall', H_DYALL
            else if (H_type==H_GPF)then
                  print*, 'H_type=GPF', H_GPF

            end if


            print*, 'ACAlpha', ACAlpha
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        Ha(ij) = ACAlpha * HNO(ij)

                        if (IndAux(i).eq.IndAux(j))then

                              temp = HNO(ij)

                              temp3 = zero
                              do r = 1, NBasis

                                    if (H_type==H_DYALL)then     	
                                          if (.not.((IndAux(r).eq.IndAux(i)).and.(IndAux(r).eq.1)))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                                temp3 = temp3 - n(r) * TwoNO(NAddr3(r,i,r,j))
                                          end if
                                    else if (H_type==H_GPF)then

                                          if (IndAux(r).ne.IndAux(i))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                          end if
                                    end if

                              end do

                              Ha(ij) = Ha(ij) + (One-ACAlpha)*temp

                        end if
                  end do
            end do

            ij = 0
            ! do i = 1, Nbasis
            !       do j = 1, i
            !             ij = ij + 1
            !             if (abs(ha(ij)).gt.1.d-8)then
            !                   print*, 'haakupa', i, j, ha(ij)
            !             end if
            !       end do
            ! end do

            print*, 'Time for Ha: ', clock_readwall(timer)
            call clock_start(timer)


            TwoNOA = TwoNO
            ij = 0
            do i = 1, NBasis
                  do j = 1, i
                        ij = ij + 1
                        kl = 0
                        do k = 1, Nbasis
                              do l = 1, k
                                    kl=kl+1

                                    if (H_type==H_DYALL)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l).and.IndAux(i)==1)then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                          end if
                                          ! if (abs(ACALPHA-one).lt.1.d-5)then
                                          !    write(*, '(A20, 4I5, 4F7.4, F20.15)') 'calka', i, j, k, l, n(i), n(j), n(k), n(l), TwoNOA(NAddr3(i, j, k, l))
                                          ! end if
                                    else if (H_type==H_GPF)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l))then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                                if (abs(TwoNOA(NAddr3(i, j, k, l))).gt.1.d-8)then
                                                end if
                                          end if
                                    end if

                              end do
                        end do
                  end do
            end do
            print*, 'Time for twonoa: ', clock_readwall(timer)
            call clock_start(timer)

            ! if (abs(ACALPHA-one).lt.1.d-5)then
            !    stop
            ! end if

      end subroutine PPERPA_init



end module ppac0
