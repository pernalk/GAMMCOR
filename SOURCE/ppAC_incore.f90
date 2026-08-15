module ppAC_incore

      use iso_fortran_env

      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals        
      use types
      use lin
      use clock

      use math_constants
      use acpp_types


      implicit none

      double precision, parameter :: tol = 1.d-4
      integer, parameter :: pDyall = 0, pGPF = 1
      integer, parameter, private :: VVOO=1, VVAO=2, VVAA=3, VAOO=4
      integer, parameter, private :: VAAO=5, VAAA=6, AAOO=7, AAAO=8
      integer, parameter, private :: VAAO2 = 9




contains

      subroutine ACPP0_incore(AuxData, TwoEl, Flags)
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            type(FlagsData), intent(in) :: Flags
            type (tclock) :: timer, timer0
            double precision, dimension(:), allocatable :: TwoNOA

            integer :: spin_symm

            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9

            character(len=6), dimension(20) :: bl_name
            character(len=5), dimension(9)  :: bbl_name
            integer, dimension(2,9) :: block_pair
            type(tMxA) :: MxA_s, MxA_t, MxS_s, MxS_t
            integer :: p, q, r, s, i, j
            integer :: twoint_dim
            integer :: dimw, spsym
            integer :: ir, il
            integer, dimension(20) :: dim_st
            double precision, dimension(9) :: E_contr_s, E_contr_t
            double precision :: contr_t1, contr_t2, contr_t3
            type(eigsBlockParams), dimension(20) :: eigsPA
            double precision :: ACAlpha,  ECorr_AC0
            integer :: variant
            integer, parameter ::  IOO=1, IVV=2, IAA=3, IVA=4, IAO=5
            integer, parameter :: IOOT=6, IVVT=7, IAAT=8, IVAT=9, IAOT = 10
            integer, parameter :: IOOT_AA=11, IVVT_AA=12, IAAT_AA=13, IVAT_AA=14, IAOT_AA = 15
            integer, parameter :: IOOT_BB=16, IVVT_BB=17, IAAT_BB=18, IVAT_BB=19, IAOT_BB = 20
            logical :: isvv, isoo

            
            

            associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux)

            
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


            call clock_start(timer0)
            ACAlpha = zero
            allocate(TwoNOA(AuxData%NInte2))
            allocate(AuxData%HNOA(NBasis, NBasis))

            call PPERPA_init(AuxData, ACAlpha, Flags, TwoEl, TwoNOA)

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
            bl_name(11) = 'oot_aa'
            bl_name(12) = 'vvt_aa'
            bl_name(13) = 'aat_aa'
            bl_name(14) = 'vat_aa'
            bl_name(15) = 'aot_aa'

            bl_name(16) = 'oot_bb'
            bl_name(17) = 'vvt_bb'
            bl_name(18) = 'aat_bb'
            bl_name(19) = 'vat_bb'
            bl_name(20) = 'aot_bb'


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

            ! Allocate 11-20
            allocate(eigsPA(IOOT_AA)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT_AA)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT_AA)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT_AA)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT_AA)%IndN(2, NA*NI))

            allocate(eigsPA(IOOT_BB)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT_BB)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT_BB)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT_BB)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT_BB)%IndN(2, NA*NI))

            do i = 1, 20
                  eigsPA(i)%IndN = 0
                  eigsPA(i)%dim = 1
            end do

            do i = 1, AuxData%NDim_s                  
                  p = AuxData%IndN(1, i)
                  q = AuxData%IndN(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOO)%IndN, eigsPA(IOO)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVV)%IndN, eigsPA(IVV)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAA)%IndN, eigsPA(IAA)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVA)%IndN, eigsPA(IVA)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAO)%IndN, eigsPA(IAO)%dim)
                  end if
            end do


            do i = 1, AuxData%NDim_t                  
                  p = AuxData%IndN_t(1, i)
                  q = AuxData%IndN_t(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOOT)%IndN, eigsPA(IOOT)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVVT)%IndN, eigsPA(IVVT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAAT)%IndN, eigsPA(IAAT)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVAT)%IndN, eigsPA(IVAT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAOT)%IndN, eigsPA(IAOT)%dim)
                  end if
            end do

            do i = 1, AuxData%NDim_t_aa
                  p = AuxData%IndN_t_aa(1, i)
                  q = AuxData%IndN_t_aa(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOOT_AA)%IndN, eigsPA(IOOT_AA)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVVT_AA)%IndN, eigsPA(IVVT_AA)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAAT_AA)%IndN, eigsPA(IAAT_AA)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVAT_AA)%IndN, eigsPA(IVAT_AA)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAOT_AA)%IndN, eigsPA(IAOT_AA)%dim)
                  end if
            end do

            do i = 1, AuxData%NDim_t_bb
                  p = AuxData%IndN_t_bb(1, i)
                  q = AuxData%IndN_t_bb(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOOT_BB)%IndN, eigsPA(IOOT_BB)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVVT_BB)%IndN, eigsPA(IVVT_BB)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAAT_BB)%IndN, eigsPA(IAAT_BB)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVAT_BB)%IndN, eigsPA(IVAT_BB)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAOT_BB)%IndN, eigsPA(IAOT_BB)%dim)
                  else
                        print*, 'to niewykorzystane', p, q
                  end if
            end do

            

            
            
            do i = 1, 20
                  eigsPA(i)%dim = eigsPA(i)%dim-1
            end do
            ! 1oo, 2vv, 3aa, 4va, 5ao            
            ! 6oot, 7vvt, 8aat, 9vat, 10aot  
            do i = 1, 20
                  dimw = eigsPA(i)%dim
                  allocate(eigsPA(i)%Eig(dimw))
                  if (i == IOO.or.i==IVV.or.i==IOOT.or.i==IVVT &
                        .or. i==IOOT_AA.or. i==IOOT_BB&
                        .or. i==IVVT_AA.or. i==IVVT_BB)then
                        allocate(eigsPA(i)%Eigvec(1, 1))
                  else
                        allocate(eigsPA(i)%Eigvec(dimw, dimw))
                  end if
                  allocate(eigsPA(i)%v_plus(dimw))
            end do

            ! i=IAA
            ! variant=1
            ! spsym = 0
            ! call eigs0_block(.false., .false., eigsPA(i)%Eig, eigsPA(i)%Eigvec, eigsPA(i)%v_plus, AuxData%HNOA, TwoNOA, &            
            !       AuxData, eigsPA(i)%IndN, eigsPA(i)%dim, spsym, Flags, variant)

 !           print*, 'wartosci wlasne tego bloku'!, i, bl_name(i), eigsPA(i)%dim
            ! do j = 1, eigsPA(i)%dim
                  
            !       if (abs(eigsPA(i)%Eig(j)).gt.1.d-5)then
                        !      if (eigsPA(i)%v_plus(j)==1)then
 !                       print*,  eigsPA(i)%Eig(j),eigsPA(i)%v_plus(j)
                        !     end if
            !       end if
            ! end do
!            stop
            
            call clock_start(timer0)
            do i = 1, 20
                  print*, ''
                  print*, 'teraz block', i, bl_name(i)
                  if (i.le.5)then
                        spsym = 0
                        variant = 1
                  else
                        spsym = 1
                        if (i.le.10) then
                              variant = 1
                        else if (i.le.15) then
                              variant = 2
                        else
                              variant = 3
                        end if
                  end if
                  if (eigsPA(i)%dim.gt.0)then
                        isvv = .false.
                        isoo = .false.
                        if (i == IOO.or.i==IOOT .or. i==IOOT_AA.or. i==IOOT_BB)then
                              isoo = .true.
                        else if (i == IVV.or. i==IVVT .or. i==IVVT_AA.or. i==IVVT_BB)then
                              isvv = .true.
                        end if

                        call eigs0_block(isvv, isoo, eigsPA(i)%Eig, eigsPA(i)%Eigvec, eigsPA(i)%v_plus, AuxData%HNOA, TwoNOA, &            
                              AuxData, eigsPA(i)%IndN, eigsPA(i)%dim, spsym, Flags, variant)

                        print*, 'wartosci wlasne tego bloku', i, bl_name(i), eigsPA(i)%dim
                        do j = 1, eigsPA(i)%dim
!                              print*, j, eigsPA(i)%Eig(j), eigsPA(i)%v_plus(j)
                                    if (abs(eigsPA(i)%Eig(j)).gt.1.d-5)then
                                    !      if (eigsPA(i)%v_plus(j)==1)then
                                                print*,  eigsPA(i)%Eig(j),eigsPA(i)%v_plus(j)
                                     !     end if
                                    end if
                              end do
                              !                       end if
                              !if (i==3)stop
                  else
                        print*, 'block', bl_name(i), 'dimension is 0'
                  end if
!                  if (i>5)stop
            end do
            print*, 'koniec blokow eigs0'
!            stop
            print*, 'TIME na all blocks ', clock_readwall(timer0)

            !    stop ! pamietaj odkomentowac parallel

            call clock_start(timer0)
            ACAlpha = one
            call PPERPA_init(AuxData, ACAlpha, Flags, TwoEl, TwoNOA)
            
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


                  call calc_block(AuxData, i, eigsPA(il), eigsPA(ir), AuxData%HNOA, TwoEl, TWONOA,&
                        0, Flags, E_contr_s(i))
                  
                  !print*, 'Czas na sing ', bbl_name(i), ' ', clock_readwall(timer0)
                  if (i <9)then
                        write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ')_s=', E_contr_s(i)
                  else
                        write(*,'(A8, A5, A3, F20.15)') 'E_contr(', bbl_name(i), ')_s=', E_contr_s(i)
                  end if
                  ECorr_AC0 = ECorr_AC0 + E_contr_s(i)
                  il = 5 + il
                  ir = 5 + ir

                  call clock_start(timer0)


                   ! Calculate Triplet Components
                  
                  ! 1. Averaged / Mixed (Variant 1)
                  print*, 'pierwszy block t1 (avg)'
                  call calc_block(AuxData, i, eigsPA(il), eigsPA(ir), AuxData%HNOA, TwoEl, TWONOA, &
                        1, Flags, contr_t1, 1) 

                  ! 2. AAAA (Variant 2)
                  ! Indices +5 from Avg
                  print*, 'drugi block t2 (aaaa)'
                  call calc_block(AuxData, i, eigsPA(il+5), eigsPA(ir+5), AuxData%HNOA, TwoEl, TWONOA, &
                        1, Flags, contr_t2, 2)

                  ! 3. BBBB (Variant 3)
                  ! Indices +10 from Avg
                  print*, 'trzeci block t3 (bbbb)'
                  call calc_block(AuxData, i, eigsPA(il+10), eigsPA(ir+10), AuxData%HNOA, TwoEl, TWONOA, &
                        1, Flags, contr_t3, 3)

                  ! Sum
!                  E_contr_t(i) = three * contr_t1 !+ contr_t2 + contr_t3
                  E_contr_t(i) = contr_t1 + contr_t2 + contr_t3
                  write(*, '(A10, F15.8, A10, F15.8, A10,F15.8)') 't_ABAB=', contr_t1, 't_AAAA=', contr_t2, 't_BBBB=', contr_t3
                  !print*, 'Czas na trip ', bbl_name(i), ' ', clock_readwall(timer0)

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
            write(*,'(A30, F25.15)') 'RDSC ACPP0 CONTR ECORR ', ECorr_AC0
            print*, ''
          end associate
    end subroutine ACPP0_incore

    subroutine ACPP0_incore_mix(AuxData, TwoEl, Flags)
          type(TACppData), intent(inout) :: AuxData
          double precision, dimension(:),      intent(in) :: TwoEl
          type(FlagsData), intent(in) :: Flags
          type (tclock) :: timer, timer0
          double precision, dimension(:), allocatable :: TwoNOA

          integer :: spin_symm

          integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4
          integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
          integer, parameter :: vaao2 = 9

          character(len=6), dimension(15) :: bl_name
          character(len=5), dimension(9)  :: bbl_name
          integer, dimension(2,9) :: block_pair
          type(tMxA) :: MxA_s, MxA_t, MxS_s, MxS_t
          integer :: p, q, r, s, i, j
          integer :: twoint_dim
          integer :: dimw, spsym
          integer :: ir, il
          double precision, dimension(9) :: E_contr_s, E_contr_t
          double precision :: contr_mix1, contr_mix2, contr_mix3, contr_mix4
          double precision :: contr_t_AAAA, contr_t_BBBB

          type(eigsBlockParams), dimension(15) :: eigsPA
          double precision :: ACAlpha,  ECorr_AC0
          integer :: variant
          integer, parameter ::  IOO=1, IVV=2, IAA=3, IVA=4, IAO=5
          integer, parameter :: IOOT_AA=6, IVVT_AA=7, IAAT_AA=8, IVAT_AA=9, IAOT_AA = 10
          integer, parameter :: IOOT_BB=11, IVVT_BB=12, IAAT_BB=13, IVAT_BB=14, IAOT_BB = 15
          integer :: nblock, offs1, offs2
          logical :: isvv, isoo


          associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                NV=>AuxData%NV, NBasis=>AuxData%NBasis, IndAux=>AuxData%IndAux)


            nblock = 15


            bbl_name(1) = 'vvoo'
            bbl_name(2) = 'vvao'
            bbl_name(3) = 'vvaa'

            bbl_name(4) = 'vaoo'
            bbl_name(5) = 'vaao'
            bbl_name(6) = 'vaaa'

            bbl_name(7) = 'aaoo'
            bbl_name(8) = 'aaao'
            bbl_name(9) = 'vaao2'


            call clock_start(timer0)
            ACAlpha = zero
            allocate(TwoNOA(AuxData%NInte2))
            allocate(AuxData%HNOA(NBasis, NBasis))

            call PPERPA_init(AuxData, ACAlpha, Flags, TwoEl, TwoNOA)

            bl_name(IOO) = 'oo'
            bl_name(IVV) = 'vv'
            bl_name(IAA) = 'aa'
            bl_name(IVA) = 'va'
            bl_name(IAO) = 'ao'

            bl_name(IOOT_AA) = 'oot_aa'
            bl_name(IVVT_AA) = 'vvt_aa'
            bl_name(IAAT_AA) = 'aat_aa'
            bl_name(IVAT_AA) = 'vat_aa'
            bl_name(IAOT_AA) = 'aot_aa'

            bl_name(IOOT_BB) = 'oot_bb'
            bl_name(IVVT_BB) = 'vvt_bb'
            bl_name(IAAT_BB) = 'aat_bb'
            bl_name(IVAT_BB) = 'vat_bb'
            bl_name(IAOT_BB) = 'aot_bb'


            allocate(eigsPA(IOO)%IndN(2, NI*NI))
            allocate(eigsPA(IVV)%IndN(2, NV*NV))
            allocate(eigsPA(IAA)%IndN(2, NA*NA))
            allocate(eigsPA(IVA)%IndN(2, 2*NV*NA))
            allocate(eigsPA(IAO)%IndN(2, 2*NA*NI))


            allocate(eigsPA(IOOT_AA)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT_AA)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT_AA)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT_AA)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT_AA)%IndN(2, NA*NI))

            allocate(eigsPA(IOOT_BB)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT_BB)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT_BB)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT_BB)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT_BB)%IndN(2, NA*NI))

            do i = 1, nblock
                  eigsPA(i)%IndN = 0
                  eigsPA(i)%dim = 1
            end do

            do i = 1, AuxData%NDim                  
                  p = AuxData%IndN(1, i)
                  q = AuxData%IndN(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOO)%IndN, eigsPA(IOO)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVV)%IndN, eigsPA(IVV)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAA)%IndN, eigsPA(IAA)%dim)
                  elseif ((IndAux(p)==2 .and. IndAux(q)==1).or.(IndAux(p)==1 .and. IndAux(q)==2)) then
                        call updateIndices(AuxData, p, q, eigsPA(IVA)%IndN, eigsPA(IVA)%dim)
                  elseif ((IndAux(p)==1 .and. IndAux(q)==0).or. (IndAux(p)==0 .and. IndAux(q)==1) )then
                        call updateIndices(AuxData, p, q, eigsPA(IAO)%IndN, eigsPA(IAO)%dim)
                  end if
            end do



            do i = 1, AuxData%NDim_t_aa
                  p = AuxData%IndN_t_aa(1, i)
                  q = AuxData%IndN_t_aa(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOOT_AA)%IndN, eigsPA(IOOT_AA)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVVT_AA)%IndN, eigsPA(IVVT_AA)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAAT_AA)%IndN, eigsPA(IAAT_AA)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVAT_AA)%IndN, eigsPA(IVAT_AA)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAOT_AA)%IndN, eigsPA(IAOT_AA)%dim)
                  end if
            end do

            do i = 1, AuxData%NDim_t_bb
                  p = AuxData%IndN_t_bb(1, i)
                  q = AuxData%IndN_t_bb(2, i)
                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IOOT_BB)%IndN, eigsPA(IOOT_BB)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices(AuxData, p, q, eigsPA(IVVT_BB)%IndN, eigsPA(IVVT_BB)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IAAT_BB)%IndN, eigsPA(IAAT_BB)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices(AuxData, p, q, eigsPA(IVAT_BB)%IndN, eigsPA(IVAT_BB)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices(AuxData, p, q, eigsPA(IAOT_BB)%IndN, eigsPA(IAOT_BB)%dim)
                  else
                        print*, 'to niewykorzystane', p, q
                  end if
            end do


            do i = 1, nblock
                  eigsPA(i)%dim = eigsPA(i)%dim-1
            end do

            do i = 1, nblock
                  dimw = eigsPA(i)%dim

                  allocate(eigsPA(i)%Eig(dimw))
                  if (i == IOO.or.i==IVV.or.i==IOOT_AA.or.i==IVVT_AA &
                        .or. i==IOOT_BB.or. i==IVVT_BB)then
                        allocate(eigsPA(i)%Eigvec(1, 1))
                  else
                        allocate(eigsPA(i)%Eigvec(dimw, dimw))
                  end if
                  allocate(eigsPA(i)%v_plus(dimw))
            end do

            offs1 = 5
            offs2 = 10


            call clock_start(timer0)
            do i = 1, nblock
                  print*, ''
                  print*, 'teraz block', i, bl_name(i)
                  if (i.le.offs1)then
                        spsym = 2
                        variant = 4
                  else
                        spsym = 1
                        if (i.le.offs2) then
                              variant = 2
                        else
                              variant = 3
                        end if
                  end if

                  if (eigsPA(i)%dim.gt.0)then
                        isvv = .false.
                        isoo = .false.
                        if (i == IOO.or. i==IOOT_AA.or. i==IOOT_BB)then
                              isoo = .true.
                        else if (i == IVV .or. i==IVVT_AA .or. i==IVVT_BB)then
                              isvv = .true.
                        end if

                        call eigs0_block(isvv, isoo, eigsPA(i)%Eig, eigsPA(i)%Eigvec, eigsPA(i)%v_plus, AuxData%HNOA, TwoNOA, &            
                              AuxData, eigsPA(i)%IndN, eigsPA(i)%dim, spsym, Flags, variant)

                        print*, 'wartosci wlasne tego bloku', i, bl_name(i)
                        print*, 'o wymiarze', eigsPA(i)%dim

                        do j = 1, eigsPA(i)%dim
                              print*,  j, ',', eigsPA(i)%Eig(j),',', eigsPA(i)%v_plus(j)
                        end do
                  else
                        print*, 'block', bl_name(i), 'dimension is 0'
                  end if

            end do
            print*, 'koniec blokow eigs0'

            !----------------------------------------------------------------------------
            print*, 'koniec blokow eigs0'
            !            stop
            print*, 'TIME na all blocks ', clock_readwall(timer0)

            call clock_start(timer0)
            ACAlpha = one
            call PPERPA_init(AuxData, ACAlpha, Flags, TwoEl, TwoNOA)

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
                  print*, 'teraz robie blok A1 mix', i, bbl_name(i)
                  il = block_pair(1, i)
                  ir = block_pair(2, i)

                  call clock_start(timer0)
                  print*, bbl_name(i), eigsPA(il)%dim, eigsPA(ir)%dim

                  contr_mix1 = zero
                  contr_t_AAAA  = zero
                  contr_t_BBBB  = zero

                  ! mixed ABAB contribution
                  print*, 'mix contrib'
                  call calc_block(AuxData, i, eigsPA(il), eigsPA(ir), AuxData%HNOA, TwoEl, TwoNOA, &
                        2, Flags, contr_mix1, 4)


                  ! AAAA triplet contribution
                  print*, 'AAAA contrib'
                  call calc_block(AuxData, i, eigsPA(il+offs1), eigsPA(ir+offs1), AuxData%HNOA, TwoEl, TwoNOA, &
                        1, Flags, contr_t_AAAA, 2)

                  ! BBBB triplet contribution
                  print*, 'BBBB contrib'
                  call calc_block(AuxData, i, eigsPA(il+offs2), eigsPA(ir+offs2), AuxData%HNOA, TwoEl, TwoNOA, &
                        1, Flags, contr_t_BBBB, 3)

                  E_contr_s(i) = contr_mix1
                  E_contr_t(i) = contr_t_AAAA + contr_t_BBBB

                  write(*, '(A12, F20.15)') 'contr_mix1', contr_mix1
                  write(*, '(A12, F20.15)') 'contr_t_AAAA', contr_t_AAAA
                  write(*, '(A12, F20.15)') 'contr_t_BBBB', contr_t_BBBB

                  if (i < 9) then
                        write(*,'(A8, A4, A7, F20.15)') 'E_contr(', bbl_name(i), ')_mix = ', E_contr_s(i)
                        write(*,'(A8, A4, A5, F20.15)') 'E_contr(', bbl_name(i), ')_t = ', E_contr_t(i)
                        write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ') = ', E_contr_s(i) + E_contr_t(i)
                  else
                        write(*,'(A8, A5, A7, F20.15)') 'E_contr(', bbl_name(i), ')_mix = ', E_contr_s(i)
                        write(*,'(A8, A5, A5, F20.15)') 'E_contr(', bbl_name(i), ')_t = ', E_contr_t(i)
                        write(*,'(A8, A5, A3, F20.15)') 'E_contr(', bbl_name(i), ') = ', E_contr_s(i) + E_contr_t(i)
                  end if

                  print*, i, bbl_name(i)
                  ECorr_AC0 = ECorr_AC0 + E_contr_s(i) + E_contr_t(i)
            end do

            print*, 'TIME na all blocks ', clock_readwall(timer)
            print*, ''
            write(*,'(A30, F25.15)') 'RDSC ACPP0 CONTR ECORR ', ECorr_AC0
            print*, ''
          end associate
    end subroutine ACPP0_incore_mix
    
     subroutine eigs0_block(isvv, isoo, Eig, Eigvec, v_plus, Ha, TwoNOA, &
            AuxData, IndN, NDim, spin_symm, Flags, variant)

            logical, intent(in) :: isvv, isoo
            double precision, dimension(:), intent(inout) :: Eig
            double precision, dimension(:,:), intent(inout) :: Eigvec
            integer, dimension(:), intent(inout) :: v_plus
            double precision, dimension(:,:),      intent(in) :: Ha
            double precision, dimension(:),      intent(in) :: TwoNOA
            type(TACppData), intent(inout) :: AuxData
            integer, dimension(:,:), intent(in)  :: IndN
            integer, intent(in) :: NDim
            integer, intent(in) :: spin_symm
            type(FlagsData), intent(in) :: Flags
            integer, optional, intent(in) :: variant
            double precision, dimension(:,:), allocatable :: MxA, MxS
            double precision, dimension(:), allocatable :: Eig_i
            type (tclock) :: timer, timer0
            integer :: i
            integer, parameter :: multiply_by_S = 1
            integer, parameter ::  IOO=1, IVV=2, IAA=3, IVA=4, IAO=5
            integer, parameter :: IOOT=6, IVVT=7, IAAT=8, IVAT=9, IAOT = 10
            integer, parameter :: IOOT_AA=11, IVVT_AA=12, IAAT_AA=13, IVAT_AA=14, IAOT_AA = 15
            integer, parameter :: IOOT_BB=16, IVVT_BB=17, IAAT_BB=18, IVAT_BB=19, IAOT_BB = 20
            integer :: var
            
            var = 0
            if (present(variant)) var = variant

            if (NDim .gt.0)then
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))
                  allocate(Eig_i(NDim))

                  !    call clock_start(timer)

                  print*, 'ten wymiar to', NDim
                  call clock_start(timer0)
                 
                  if (var == 2) then
                        ! AAAA
                        call pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, Ha, TwoNOA, Ndim, Ndim, &
                              indn, indn, AuxData%indx_t_aa,  multiply_by_S)
                  else if (var == 3) then
                        ! BBBB
                        call pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, Ha, TwoNOA, Ndim, Ndim, &
                              indn, indn, AuxData%indx_t_bb, multiply_by_S)
                  else
                        if (AuxData%spinsep ==.true.)then
                              if (spin_symm == 0)then
                                    call pperpa_incore_opshell(MxA, MxS, AuxData, Ha, TwoNOA, Ndim, Ndim, &
                                          indn, indn, AuxData%indx_s, spin_symm, multiply_by_S)
                              else
                                    call pperpa_incore_opshell(MxA, MxS, AuxData, Ha, TwoNOA, Ndim, Ndim, &
                                          indn, indn, AuxData%indx_t, spin_symm, multiply_by_S)
                              end if
                        else
                              call pperpa_incore_opshell_mix2(MxA, MxS, AuxData, Ha, TwoNOA, Ndim, Ndim, &
                                    indn, indn, AuxData%indx, multiply_by_S)
                        end if
                  end if


                  if (isvv .or. isoo)then
                        Eigvec = zero
                        do i = 1, NDim
                              Eig(i) = MxA(i, i)
                        end do
                        if (isoo)then
                              v_plus = 0
                        else
                              v_plus = 1
                        end if
                  else
                        print*, 'nonsym'
                        call nonsymmetric_eigenproblem_block(Eig, Eig_i, Eigvec, MxA, MxS, v_plus, IndN, AuxData%IndAux, 1)
                  end if

                  deallocate(MxA)
                  deallocate(MxS)
                  deallocate(Eig_i)
            end if
            ! print*, 'Czas diagonalizacji sing: '//str(clock_readwall(timer),d=2)


      end subroutine eigs0_block

      subroutine calc_block(AuxData, x, ePAl, ePAr,  Ha, TwoEl,  TWONOA, &
            spin_symm, Flags, E_contr, variant)
            type(TACppData), intent(inout) :: AuxData
            integer, intent(in) :: x
            type(eigsBlockParams), intent(in) :: ePAl, ePAr
            double precision, dimension(:, :),      intent(in) :: Ha
            double precision, dimension(:),      intent(in) :: TwoEl, TWONOA
            integer, intent(in) :: spin_symm
            type(FlagsData), intent(in) :: Flags
            double precision, intent(out) :: E_contr
            integer, optional, intent(in) :: variant
            integer :: var
            external :: dgemv

            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4                                    
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9
            integer :: p, q, r, s, i, j, pq
            integer :: np,km
            double precision :: Npqrs
            double precision :: Aux1, Aux2, eig_scaling, eig_scaling2
            double precision :: fpq,frs,gpq,grs, Aux3
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
            integer :: xi, yi
            double precision :: xau
            


            associate(ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, Iaux=>AuxData%IndAux, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m)



              E_contr = zero

              
            allocate(MxA1(ePAl%dim, ePAr%dim))
            allocate(MxS(1,1))

            allocate(IMxA2(ePAr%dim, ePAl%dim))
             allocate(ttt(ePAl%dim, ePAr%dim))

            call clock_start(timer)
            call clock_start(timerall)
            
            if (present(variant)) var = variant

            if (var == 2) then
                  ! AAAA
                   call pperpa_incore_opshell_aaaa(MxA1, MxS, AuxData, AuxData%HNOA, TwoNOA, ePAl%dim, ePAr%dim, &
                        ePAl%IndN, ePAr%IndN,  AuxData%IndX_t_aa, multiply_by_S)
            else if (var == 3) then
                  ! BBBB
                   call pperpa_incore_opshell_bbbb(MxA1, MxS, AuxData, AuxData%HNOA, TwoNOA, ePAl%dim, ePAr%dim, &
                        ePAl%IndN, ePAr%IndN,  AuxData%IndX_t_bb, multiply_by_S)
             else
                   
                   if (AuxData%spinsep ==.true.)then
                         
                         if (spin_symm == 0)then
                               call pperpa_incore_opshell(MxA1, MxS, AuxData, AuxData%HNOA, TwoNOA, ePAl%dim, ePAr%dim, &
                                     ePAl%IndN, ePAr%IndN,  AuxData%IndX_s, spin_symm, multiply_by_S)
                         else if (spin_symm == 1)then
                               call pperpa_incore_opshell(MxA1, MxS, AuxData, AuxData%HNOA, TwoNOA, ePAl%dim, ePAr%dim, &
                                     ePAl%IndN, ePAr%IndN,  AuxData%IndX_t, spin_symm, multiply_by_S)
                               
                         end if
                   else
                         call pperpa_incore_opshell_mix2(MxA1, MxS, AuxData, AuxData%HNOA, TwoNOA, ePAl%dim, ePAr%dim, &
                               ePAl%IndN, ePAr%IndN,  AuxData%IndX, multiply_by_S)
                   end if
                   
             end if
                                    
            deallocate(MxS)
            allocate(IMxA1(ePAl%dim, ePAr%dim))
            do i = 1, ePAl%dim
                  do j = 1, ePAr%dim
                        IMxA2(j, i) = MxA1(i, j)
                  end do
            end do

            allocate(tempx(ePAl%dim))
            allocate(tempx2(ePAl%dim))

                  
            if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then

                  call clock_start(timer)
                  IMxA1 = zero
                  do np = 1, ePAl%dim           
                        do km = 1, ePAr%dim
                              if (ePAl%v_plus(np)==1.and.ePAr%v_plus(km)==0)then
                                    call real_av_x(tempx, MxA1, ePAl%dim,  ePAr%Eigvec(:,km), ePAl%dim, ePAr%dim, one, zero)
                                    call real_vw_x(Aux2, ePAl%Eigvec(:,np), tempx, ePAl%dim)
                                    Aux2 = Aux2 /(ePAl%Eig(np)-ePAr%Eig(km))
                                    do j = 1, ePAr%dim

                                          IMxA1(np, j) = IMxA1(np, j) + Aux2 * ePAr%Eigvec(j, km)
                                    end do
                              end if
                        end do
                  end do                  
            end if


            call clock_start(timer)

            select case(x)

            case(VVAO, VVAA)

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


                                          if ((spin_symm == 1 .and. var ==1).or.spin_symm==0)then
                                                Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))

                                                if (spin_symm == 0) then
                                                      fpq = one
                                                      frs = one
                                                      
                                                      if (p == q) then
                                                            fpq = sqrt(frac12)
                                                      end if
                                                      
                                                      if (r == s) then
                                                            frs = sqrt(frac12)
                                                      end if
                                                      
                                                      Aux1 = fpq * frs * Npqrs * (TwoEl(NAddr3(r,p,s,q)) + TwoEl(NAddr3(r,q,s,p)))
                                                      
                                                else if (spin_symm == 1) then
                                                      
                                                      Aux1 = Npqrs  * (TwoEl(NAddr3(r,p,s,q)) -TwoEl(NAddr3(r,q,s,p)))
                                                end if

                                                
                                          else if (spin_symm == 1 .and. var ==2)then

                                                Npqrs = (one-n_p(p)-n_p(q)) * (one - n_p(r)-n_p(s))
                                                Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))

                                          else if (spin_symm == 1 .and. var ==3)then
                                                Npqrs = (one-n_m(p)-n_m(q)) * (one - n_m(r)-n_m(s))
                                                Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                          end if

                                          if (var == 4)then                                                
                                                Npqrs = (one-n_p(p)-n_m(q)) * (one - n_p(r)-n_m(s))
                                                Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q)) 
                                          end if

                                          E_contr = E_contr - Aux2 * Aux1 / (ePAl%Eig(i)-ePAr%Eig(km)) *ePAr%Eigvec(j, km)

                                    end do
                              end if
                        end do
                  end do
!                  print*, 'Czas na drugie petle: ', clock_readwall(timer)

            case(VAOO, AAOO)

                  call clock_start(timer)

                  do j = 1, ePAr%dim
                        r = ePAr%IndN(1, j)
                        s = ePAr%IndN(2, j)
                        do np = 1, ePAl%dim
                              if (ePAl%v_plus(np)==1)then


                                    call real_vw_x(Aux2, MxA1(:, j), ePAl%Eigvec(:,np), ePAl%dim)

                                    do i = 1, ePAl%dim
                                          p = ePAl%IndN(1, i)
                                          q = ePAl%IndN(2, i)
                                          
                                          if ((spin_symm == 1 .and. var ==1).or.spin_symm==0)then
                                                Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))
                                                
                                                if (spin_symm == 0) then
                                                      fpq = one
                                                      frs = one

                                                      if (p == q) then
                                                            fpq = sqrt(frac12)
                                                      end if

                                                      if (r == s) then
                                                            frs = sqrt(frac12)
                                                      end if

                                                      Aux1 = fpq * frs * Npqrs * (TwoEl(NAddr3(r,p,s,q)) + TwoEl(NAddr3(r,q,s,p)))

                                                else if (spin_symm == 1) then
                                                      Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                                end if


                                          else if (spin_symm == 1 .and. var ==2)then

                                                Npqrs = (one-n_p(p)-n_p(q)) * (one - n_p(r)-n_p(s))
                                                Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))

                                          else if (spin_symm == 1 .and. var ==3)then
                                                Npqrs = (one-n_m(p)-n_m(q)) * (one - n_m(r)-n_m(s))
                                                Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                          end if
                                           if (var == 4)then
                                                 Npqrs = (one-n_p(p)-n_m(q)) * (one - n_p(r)-n_m(s))
                                                Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q))
                                          end if


                                          E_contr = E_contr - Aux2 * Aux1 / (ePAl%Eig(np)-ePAr%Eig(j)) *ePAl%Eigvec(i,np)

                                    end do
                              end if
                        end do
                  end do

!                  print*, 'Czas na trzecie petle: ', clock_readwall(timer)

            case default 
                  print*, 'default'
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

                              if (var==4)then

                                    if (x==VAAO)then
                                          if ((Iaux(r)==1.and.Iaux(q)==1).or.(Iaux(p)==1.and.Iaux(s)==1))then
                                                continue
                                          else
                                                cycle
                                          end if
                                    else if (x==VAAO2)then
                                          if ((Iaux(q)==1.and.Iaux(s)==1).or.(Iaux(p)==1.and.Iaux(r)==1))then
                                                continue
                                          else
                                                cycle
                                          end if
                                    end if
                              end if
                                          

                              
                              if ((spin_symm == 1 .and. var ==1).or.spin_symm==0)then
                                    Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))
                                                
                                    if (spin_symm == 0) then
                                          fpq = one
                                          frs = one
                                          
                                          if (p == q) then
                                                fpq = sqrt(frac12)
                                          end if
                                          
                                          if (r == s) then
                                                frs = sqrt(frac12)
                                          end if

                                          if (x == VAAO)then
                                                Aux1 = fpq * frs * Npqrs * TwoEl(NAddr3(r,p,s,q))
                                          else if (x == VAAO2)then
                                                Aux1 = fpq * frs * Npqrs * TwoEl(NAddr3(r,q,s,p))
                                          else
                                                Aux1 = fpq * frs * Npqrs * (TwoEl(NAddr3(r,p,s,q)) + TwoEl(NAddr3(r,q,s,p)))
                                          end if
                                          
                                    else if (spin_symm == 1) then
                                          if (x == VAAO)then
                                                Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q))
                                          else if (x == VAAO2)then
                                                Aux1 = - Npqrs * TwoEl(NAddr3(r,q,s,p))
                                          else
                                                Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                          end if
                                                  
                                    end if
                                    

                              else if (spin_symm == 1 .and. var ==2)then
                                    
                                    Npqrs = (one-n_p(p)-n_p(q)) * (one - n_p(r)-n_p(s))
                                    if (x == VAAO)then
                                          Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q)) 
                                    else if (x == VAAO2) then
                                          Aux1 = -Npqrs *  TwoEl(NAddr3(r,q,s,p))
                                    else
                                          Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                    end if


                              else if (spin_symm == 1 .and. var ==3)then
                                    Npqrs = (one-n_m(p)-n_m(q)) * (one - n_m(r)-n_m(s))
                                    if (x == VAAO)then
                                          Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q)) 
                                    else if (x == VAAO2) then
                                          Aux1 = -Npqrs *  TwoEl(NAddr3(r,q,s,p))
                                    else
                                          Aux1 = Npqrs * (TwoEl(NAddr3(r,p,s,q)) - TwoEl(NAddr3(r,q,s,p)))
                                    end if
                                    
                              end if

                              if (var == 4)then
                                    Npqrs = (one-n_p(p)-n_m(q)) * (one - n_p(r)-n_m(s))
                                    Aux1 = Npqrs * TwoEl(NAddr3(r,p,s,q))
                              end if
                              
                              if (x==VVOO)then
                                    E_contr = E_contr - MxA1(i, j) * Aux1 / (ePAl%Eig(i)-ePAr%Eig(j))
                                    
                              else if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then

                                    do np = 1, ePAl%dim
                                          E_contr = E_contr -  Aux1 * ePAl%Eigvec(i,np) * IMxA1(np, j)
                                    end do
                              end if
                              
                        end do j_colloops
                  end do i_rowloops
!                  print*, 'Czas na czwarte petle: ', clock_readwall(timer)

            end select
            !       !omp end parallel do

            ! print*, 'Czas na timeall  petle: ', clock_readwall(timerall)
            ! print*, ''

            deallocate(MxA1)
            deallocate(tempx)

          end associate
    end subroutine calc_block


      subroutine updateIndices(AuxData, p, q, IndN, ind)
            type(TACppData), intent(inout) :: AuxData
            integer, intent(in) :: p, q
            integer, dimension(:,:), intent(inout) :: IndN
            integer, intent(inout) :: ind

            IndN(1,ind) = p
            IndN(2,ind) = q
            ind = ind + 1

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
            !print*, 'info', info
            if (info /= 0) then
                  print*, "Nonsymmetric matrix eigendecompositino failed with info="
                  error stop
            end if
            !print*, 'Czas diagonalizacji real: ', clock_readwall(timer)

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
            ! print*, 'eigenvalues of singlets'
            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'COMPLEX EIGENVALUES', i, wr(i), wi(i), v_plus(i)
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
                              !print*, 'max_hh', max_hh
                        end if
                        !print*, 'min_pp' , wr_plus(i)
                        !print*, 'shift', shift
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

            ! print*, 'posortowane'
            ! do i = 1, countj_plus
            !       print*, wr_plus(i)
            ! end do

      end subroutine nonsymmetric_eigenproblem_block


      subroutine pperpa_incore_opshell(MxA, MxS, AuxData, h, TwoEl, Ndim1, Ndim2, indn1, indn2, indx, spin_symm, multiply_by_S)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(in) :: h(:,:)
            integer, intent(in) :: NDim1, Ndim2
            integer, dimension(:,:), intent(in) ::indn1, indn2
            integer, dimension(:), intent(in) ::indx
            integer, intent(in) :: spin_symm
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), intent(inout) :: MxA, MxS

            integer :: i, j, t, u, v
            integer p, q, r, s, rs, pq
            double precision :: Arspq, Arsqp
            double precision :: P0, P1, P2, P3, P4
            double precision :: P1a, P1b
            double precision :: P3a1, P3a2, P3b1, P3b2, P3c12, P3c34
            double precision :: P45a1, P45a2, P45b1, P45b2

            double precision :: D0, D1, D2, D3, D4
            double precision :: P3a3, P3a4, P3b3, P3b4, P3c56, P3c78
            double precision :: P45a3, P45a4, P45b3, P45b4
            double precision :: num_f, num_h, Arssum
            double precision :: nnn
            double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
            type (tclock) :: timer0, timer1, timerxx
            integer, external :: NAddr3


            associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m)


              MxA = zero
              MxS = zero
              nnn = zero
!              !$omp parallel do collapse(2)&                                                                                                                                                             
!              !$omp default(shared) &                                                                                                                                                                     
!              !$omp private(i, j, t, u, v) &                                                                                                                                                  
!              !$omp private(p, q, r, s, rs, pq, Arspq, Arsqp, Arssum, num_f, num_h)                           
!              !$omp private(P0, P1, P2, P3, P4, P3a1, P3a2, P3b1, P3b2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2)
!              !$omp private(D0, D1, D2, D3, D4, P3a3, P3a4, P3b3, P3b4, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4)

              tp30 = zero
              tpx0 = zero
              tp31 = zero
              tp32 = zero
              tp41 = zero

              tpplusz = zero
!              print*, 'ndim1', ndim1, ndim2
              call clock_start(timerxx)
              i_rowloops: do i = 1, NDim1
                    j_colloops: do j = 1, NDim2
                          call clock_start(timer0)
                          call clock_start(timer1)                          
                          r = indn1(1, i)
                          s = indn1(2, i)
                          rs = indx(i)

                          p = indn2(1, j)
                          q = indn2(2, j)
                          pq = indx(j)

                          Arspq = zero
                          Arsqp = zero
                          P1 = zero; D1 = zero; P2 = zero; D2 = zero

                          P0 =  two * twoel(gmap(p,r,q,s)) 
                          D0 = -two * twoel(gmap(q,r,p,s)) 
                          
                          P0 = P0 - (n_p(p) + n_m(p) + n_p(q) + n_m(q) + n_p(r) + n_m(r) + n_p(s) + n_m(s))* twoel(gmap(p,r,q,s))
                          D0 = D0 + (n_p(p) + n_m(p) + n_p(q) + n_m(q) + n_p(r) + n_m(r) + n_p(s) + n_m(s))* twoel(gmap(q,r,p,s))

                          P1a = zero
                          P1b = zero
                          if (q == s) then
                                P1a = P1a + h(p, r)*two
                                P1b = P1b + h(p, r) * (- n_m(q) - n_p(q) &
                                      -frac12 * (n_m(p) + n_p(p) + n_m(r) + n_p(r))) 

                                P1 = P1a+P1b
                                ! P1 = P1 + frac12 * h(p, r) * (two - n_m(q) - n_p(q) &
                                !       -frac12 * (n_m(p) + n_p(p) + n_m(r) + n_p(r))) 
                          end if

!                          print*, 'p111111', p1, h(p, r)

                          if (p == r) then
                                P1a = P1a + h(q, s)*two
                                
                                P1b = P1b + h(q, s) * (- n_m(p) - n_p(p) &
                                      -frac12 * (n_m(q) + n_p(q) + n_m(s) + n_p(s)))

                                P1 = P1a+P1b
                                ! P1 = P1 + frac12 * h(q, s) * (two - n_m(p) - n_p(p) &
                                !       -frac12 * (n_m(q) + n_p(q) + n_m(s) + n_p(s)))
                          end if

                          if (p == s) then
                                D1 = D1 - h(q, r) * (two - n_m(p) - n_p(p) &
                                      -frac12 * (n_m(q) + n_p(q) + n_m(r) + n_p(r)) )
                          end if

                          if (q == r) then
                                D1 = D1 - h(p, s) * (two - n_m(q) - n_p(q) &
                                      -frac12 * (n_m(p) + n_p(p) + n_m(s) + n_p(s))) 
                          end if

                          ! if (q == s) then
                          !       P1 = P1 + frac12 * h(p, r) * (one - n_m(q) - n_p(q) &
                          !             -frac12 * (n_m(p) + n_p(p) + n_m(r) + n_p(r))) 
                          ! end if

                          ! if (p == r) then
                          !       P1 = P1 + frac12 * h(q, s) * (one - n_m(p) - n_p(p) &
                          !             -frac12 * (n_m(q) + n_p(q) + n_m(s) + n_p(s)))
                          ! end if

                          ! if (p == s) then
                          !       D1 = D1 - frac12 * h(q, r) * (one - n_m(p) - n_p(p) &
                          !             -frac12 * (n_m(q) + n_p(q) + n_m(r) + n_p(r)) )
                          ! end if

                          ! if (q == r) then
                          !       D1 = D1 - frac12 * h(p, s) * (one - n_m(q) - n_p(q) &
                          !             -frac12 * (n_m(p) + n_p(p) + n_m(s) + n_p(s))) 
                          ! end if

                          P2 = zero; D2 = zero
                          if (q == s) then
                                !p2 = p2  - n(t) * (TwoNOA(NAddr3(p,t,t,r)) - twop* TwoNOA(NAddr3(p,r,t,t)))

                                  do t = 1, NIA
                                        P2 = P2 - (n_p(t) + n_m(t)) * (twoel(gmap(p, t, r, t)) - two *twoel(gmap(p, r, t, t)))

                                        ! if (rs==1.and.pq==1)then
                                        !       if (abs((n_p(t) + n_m(t)) * (twoel(gmap(p, t, r, t)) - two *twoel(gmap(p, r, t, t)))).gt.1.d-5)then
                                        !             write(*, '(A5, 4I5, A1, 4I5, 5F20.10)') 'sliz', p, r, t, t, '|', p, t, t, r,n_p(t)+n_m(t), twoel(gmap(p, r, t, t)), twoel(gmap(p, t, t, r)), &
                                        !                   - (n_p(t) + n_m(t)) * (twoel(gmap(p, t, r, t)) - two *twoel(gmap(p, r, t, t))), p2
                                        !       end if
                                        ! end if                                       
                                end do
                          end if
                          
                          if (p == r) then
                                do t = 1, NIA
                                      P2 = P2 - (n_p(t) + n_m(t)) * (twoel(gmap(q, t, s, t)) - two *twoel(gmap(q, s, t, t)))
                                      ! if (rs==1.and.pq==1)then
                                      !       if (abs((n_p(t) + n_m(t)) * (twoel(gmap(q, t, s, t)) - two *twoel(gmap(q, s, t, t)))).gt.1.d-5)then
                                      !             write(*, '(A5, 3I5, 5F20.10)') 'sliz2', t, q, s,n_p(t)+n_m(t), twoel(gmap(q, s, t, t)), twoel(gmap(q, t, t, s)), &
                                      !                   - (n_p(t) + n_m(t)) * (twoel(gmap(q, t, s, t)) - two *twoel(gmap(q, s, t, t))), p2
                                      !       end if
                                      !   end if                                       

                                end do
                          end if

                          if (p == s) then
                                do t = 1, NIA
                                      D2 = D2 + (n_p(t) + n_m(t)) * (twoel(gmap(q, t, r, t)) - two *twoel(gmap(q, r, t, t))) 
                                end do
                          end if
                          if (q == r) then
                                do t = 1, NIA
                                      D2 = D2 + (n_p(t) + n_m(t)) * (twoel(gmap(p, t, s, t)) - two *twoel(gmap(p, s, t, t))) 
                                end do
                          end if
                          tpx0 = tpx0 + clock_readwall(timer0)

                          P3a1 = zero; P3a2 = zero; P3b1 = zero; P3b2 = zero; P3c12 = zero; P3c34 = zero
                          P3a3 = zero; P3a4 = zero; P3b3 = zero; P3b4 = zero; P3c56 = zero; P3c78 = zero
                          P45a1 = zero; P45a2 = zero; P45b1 = zero; P45b2 = zero
                          P45a3 = zero; P45a4 = zero; P45b3 = zero; P45b4 = zero

                          call clock_start(timer0)

                          P3a1 = P3a1 + func_P3a(AuxData, TwoEl, s, p, q, r)
                          P3a2 = P3a2 + func_P3a(AuxData, TwoEl, r, q, p, s)
                          P3a3 = P3a3 - func_P3a(AuxData, TwoEl, s, q, p, r)
                          P3a4 = P3a4 - func_P3a(AuxData, TwoEl, r, p, q, s)

                          tp30 = tp30 + clock_readwall(timer0)
                          call clock_start(timer0)


                          P3b1 = P3b1 + func_P3b(AuxData, TwoEl, s, q, p, r)
                          P3b2 = P3b2 + func_P3b(AuxData, TwoEl, r, p, q, s)
                          P3b3 = P3b3 - func_P3b(AuxData, TwoEl, s, p, q, r)
                          P3b4 = P3b4 - func_P3b(AuxData, TwoEl, r, q, p, s)
                           

                          tp31 = tp31 + clock_readwall(timer0)
                          call clock_start(timer0)
                          P3c12 = P3c12 - func_P3c(AuxData, TwoEl, s, q, p, r)
                          P3c34 = P3c34 - func_P3c(AuxData, TwoEl, r, p, q, s)
                          P3c56 = P3c56 + func_P3c(AuxData, TwoEl, s, p, q, r)
                          P3c78 = P3c78 + func_P3c(AuxData, TwoEl, r, q, p, s)

                          tp32 = tp32 + clock_readwall(timer0)

                          call clock_start(timer0)
                          if (p == r) then
                                P45a1 = P45a1 - func_P45(AuxData, TwoEl, q, s)
                                P45b1 = P45b1 - func_P45(AuxData, TwoEl, s, q)
                          end if

                          if (q == s) then
                                P45a2 = P45a2 - func_P45(AuxData, TwoEl, p, r)
                                P45b2 = P45b2 - func_P45(AuxData, TwoEl, r, p)
                          end if

                          if (q == r) then
                                P45a3 = P45a3 + func_P45(AuxData, TwoEl, p, s)
                                P45b3 = P45b3 + func_P45(AuxData, TwoEl, s, p)
                          end if

                          if (p == s) then
                                 P45a4 = P45a4 + func_P45(AuxData, TwoEl, q, r)
                                 P45b4 = P45b4 + func_P45(AuxData, TwoEl, r, q)
                           end if

                          tp41 = tp41 + clock_readwall(timer0)
                          !if(r==10.and.s==8.and.p==10.and.q==8)then
                          ! if (rs==5.and.pq==5)then
                          !       write(*,'(A5, F20.10)') 'p0', p0
                          !       write(*,'(A5, F20.10)') 'p1', p1
                          !       write(*,'(A5, F20.10)') 'p2', p2
                          !       write(*,'(A5, F20.10)') 'p3a1', p3a1
                          !       write(*,'(A5, F20.10)') 'p3a2', p3a2
                          !       write(*,'(A5, F20.10)') 'p3b1', p3b1
                          !       write(*,'(A5, F20.10)') 'p3b2', p3b2
                          !       write(*,'(A5, F20.10)') 'p3c12', p3c12
                          !       write(*,'(A5, F20.10)') 'p3c34', p3c34
                          !       write(*,'(A5, F20.10)') 'p45a1', p45a1
                          !       write(*,'(A5, F20.10)') 'p45a2', p45a2
                          !       write(*,'(A5, F20.10)') 'p45b1', p45b1
                          !       write(*,'(A5, F20.10)') 'p45b2', p45b2
                          !       print*, ''
                          !       write(*,'(A5, F20.10)') 'd0', d0
                          !       write(*,'(A5, F20.10)') 'd1', d1
                          !       write(*,'(A5, F20.10)') 'd2', d2
                          !       write(*,'(A5, F20.10)') 'p3a3', p3a3
                          !       write(*,'(A5, F20.10)') 'p3a4', p3a4
                          !       write(*,'(A5, F20.10)') 'p3b3', p3b3
                          !       write(*,'(A5, F20.10)') 'p3b4', p3b4
                          !       write(*,'(A5, F20.10)') 'p3c56', p3c56
                          !       write(*,'(A5, F20.10)') 'p3c78', p3c78
                          !       write(*,'(A5, F20.10)') 'p45a3', p45a3
                          !       write(*,'(A5, F20.10)') 'p45a4', p45a4
                          !       write(*,'(A5, F20.10)') 'p45b3', p45b3
                          !       write(*,'(A5, F20.10)') 'p45b4', p45b4

                          ! end if

                          ! write(*, '(A5, 4A3, 28A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '     P1a', '     P1b', '      P2', '   P3c12',&
                          !       '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
                          !       '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
                          ! write(*, '(A5, 4I3, 28F8.4)') 'plabab', r, s, p, q, &
                          !       P0, P1, P1a, P1b, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                          !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                          !       P3a1, P3a2, P3a3, P3a4

!                           if (any(abs([P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
!              D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                          !              P3a1, P3a2, P3a3, P3a4]) .gt. 1.d-5)) then
 !                          if (rs==1.and.pq==1)then
 !    write(*, '(A5, 4A3, 26A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '      P2', '   P3c12',&
 !        '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
 !        '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
 !    write(*, '(A5, 4I3, 26F8.4)') 'plabab', r, s, p, q, &
 !        P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
 !        D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
 !        P3a1, P3a2, P3a3, P3a4
 ! end if


                          !                           write(*, '(A5, 4A3, 26A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '      P2', '   P3c12',&
                          !       '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
                          !       '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
                          ! write(*, '(A5, 4I3, 26F8.4)') 'plabab', r, s, p, q, &
                          !       P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                          !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                          !       P3a1, P3a2, P3a3, P3a4

                          

                          Arspq = frac12*(P0 + P1 + P2 + P3a1 + P3a2 + P3b1 + P3b2 + P3c12 + P3c34 + P45a1+ P45a2 + P45b1 + P45b2)
                          
                          Arsqp = frac12*(D0 + D1 + D2 + P3a3 + P3a4 + P3b3 + P3b4 + P3c56 + P3c78 + P45a3 + P45a4 + P45b3 + P45b4)
!                          write(8, '(A10, 4I5, 8F12.6)') 'kla', r, s, p, q, p3a1, p3a2, p3a3, p3a4, p3b1, P3b2, P3b3, p3b4
                          
!          write(*, '(4I5, 10F13.5)') r, s, p, q, p3a1, p3a2, p3b1, p3b2, Arspq, Arsqp

                          num_f = one

                          if (multiply_by_S==1)then
                                !num_f = one / (two - (n_p(r) + n_m(s) + n_m(r) + n_p(s)))
                                num_f = one / (frac12* (two - n_p(r) - n_m(s) - n_m(r) - n_p(s)))
!                                write(*, '(A10, 4I5, 4F12.8)')'numf', r, s, r, s, n_p(r),  n_m(s) , n_m(r) , n_p(s)

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

                                MxA(rs, pq) = num_f * num_h * Arssum


                                ! if (abs(MxA(rs, pq)).gt.1.d-5)then
                                !       write(*,'(A10, 2I5, A5, 4 I5,6F20.10)') 'luzluz', rs, pq, '  |  ', p, q, r, s, MxA(rs, pq), num_f, num_h, Arspq, Arsqp, Arssum
                                !       !print*, ''
                                ! end if


                                ! if (abs(MxA(rs, pq)).gt.1.d-5)then
                                !       write(*, '(A10, 2I5, 5F15.8)') 'mxamxa', rs, pq, MxA(rs, pq), Arspq, Arsqp, num_h, num_f
                                ! end if

                                ! if (r==2.and.s==1.and.p==2.and.q==1)then
                                !       write(*,'(6F20.5)') num_h, num_f, arspq, arsqp, arssum, mxa(rs, pq)
                                ! end if

                                !                                if (abs(mxa(rs,pq)).gt.1.d-5)then
                                !write(*,'(2I5, 5F20.5)') rs, pq, mxa(rs, pq), Arspq, Arsqp, num_f, num_h
                                !       write(*, '(A10, 4I5, 4F15.5)') 'singll', r, s, p, q, &
                                !             (one - frac12 * (n_p(r) + n_m(s) + n_m(r) + n_p(s))), num_h, Arspq, Arsqp
                                ! !       write(*,'(I5,A1,I5,A1,I5,A1,I5,A1,I5,A1,I5,A1,F20.5,A1,F20.5,A1,F20.5,A1,F20.5)') rs, ',', pq, ',', r, ',', s,',', p,',', q,',', MxA(rs,pq), ',', num_f, ',', num_h, ',', Arssum
                                !  end if

!                                                             if (abs(MxA(rs, pq)).gt.1.d-5)then
                                !print*, Arspq, Arsqp, num_f, erdm(rdm1_p, r, r, IAux, NI), erdm(rdm1_m, s, s, IAux, NI), Arssum
                                                                   !write(*, '(4I5, 6F15.5)') r, s, p, q, num_f, num_h, Arspq, Arsqp, erdm(rdm1_p, r, r, IAux, NI), erdm(rdm1_m, s, s, IAux, NI)
 !                       write(*, '(4I5, 4F15.5)') r, s, p, q, num_f, num_h, Arspq, Arsqp
  !                              end if
                                

                          case(1)
                                if (p.ne.q.and.r.ne.s)then
                                      Arssum = Arspq + Arsqp
                                      MxA(rs, pq) = num_f * Arssum
                                      ! if (abs(MxA(rs, pq)).gt.1.d-5)then
                                      !       write(*,'(A10, 2I5, A5, 4I5,2F30.16)') 'zulzul', rs, pq, '  |  ', p, q, r, s, MxA(rs, pq)
                                      ! end if
                                      !       write(*, '(A10, 2I5, 5F15.8)') 'mxamxa', rs, pq, MxA(rs, pq), Arspq, Arsqp, num_h, num_f
                                      ! end if

!                                      if (abs(MxA(rs, pq)).gt.1.d-5)then                                      
!                                            write(*, '(A5, 6I5, 3F15.5)') 'gow', r, s, p, q, rs, pq, (one - frac12 * (n_p(r) + n_m(s) + n_m(r) + n_p(s))), MxA(rs, pq), num_f

                                            ! write(*, '(A5, 4A3, 28A12)')   '     ', ' r', ' s', ' p', ' q', &
                                            !       '          P0', '          P1', '         P1a', '         P1b', '          P2', &
                                            !       '       P3c12', '       P3c34', '       P45a1', '       P45a2', '       P45b1', '       P45b2', &
                                            !       '          D0', '          D1', '          D2', &
                                            !       '       P3c56', '       P3c78', '       P45a3', '       P45a4', '       P45b3', '       P45b4', &
                                            !       '        P3b1', '        P3b2', '        P3b3', '        P3b4', &
                                            !       '        P3a1', '        P3a2', '        P3a3', '        P3a4'

                                            ! write(*, '(A5, 4I3, 28F12.4)') 'plabab', r, s, p, q, &
                                            !       P0, P1, P1a, P1b, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                                            !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                                            !       P3a1, P3a2, P3a3, P3a4
                                            
                                            
                                            ! write(*, '(A5, 4A3, 28A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '     P1a', '     P1b', '      P2', '   P3c12',&
                                            !       '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
                                            !       '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
                                            ! write(*, '(A5, 4I3, 28F8.4)') 'plabab', r, s, p, q, &
                                            !       P0, P1, P1a, P1b, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                                            !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                                            !       P3a1, P3a2, P3a3, P3a4

!                                            write(*, '(A5, 4F15.5)') 'gow2', n_p(r), n_p(s),   n_p(r) + n_p(s), one -n_p(r) - n_p(s)
 !                                     end if

                                      ! if (r==18.and.s==5.and.p==5.and.q==4)then                                      
                                      
                                      !end if

!                                      if (abs(Arspq).gt.1.d-5)then
!                                      write(*, '(4I5, 11F10.5)') r, s, p, q, p3a1, p3a2, p3b1, p3b2, P3a3 , &
!                                            P3a4 , P3b3 , P3b4 , Arspq, Arsqp, Arssum
                                      !                                end if
                                      
 !                                     if (abs(MxA(rs, pq)).gt.1.d-5)then
!                                            write(*, '(A10, 4I5, 3F15.5)') 'tripll', r, s, p, q, &
!                                                  (one - frac12 * (n_p(r) + n_m(s) + n_m(r) + n_p(s))), Arspq, Arsqp
!                                                   write(*, '(6I5, F15.5)') r, s, p, q, rs, pq, MxA(rs, pq)
                                                   
                       !        write(*, '(4I5, 8F15.5)') r, s, p, q, n_p(r), n_m(r), n_p(s), n_m(s), num_f, num_h, Arspq, Arsqp
!                                                          end if
                                end if
                          end select
                          tpplusz = tpplusz + clock_readwall(timer1)                          
                    end do j_colloops
              end do i_rowloops
!!$omp end parallel do
              ! print*, 'tp30', tp30
              ! print*, 'tp31', tp31
              ! print*, 'tp32', tp32
              ! print*, 'tp41', tp41
              ! print*, 'tpx0', tpx0
              ! print*, 'tpplusz', tpplusz

              ! print*, 'full time', clock_readwall(timerxx)
              MxS = zero

              if (multiply_by_s==1)then
                    do i = 1, NDim1
                          p = IndN1(1, I)
                          q = IndN1(2, i)
                          MxS(i, i) = frac12* (two - n_p(p) - n_m(q) - n_m(p) - n_p(q))
!                          MxS(i, i) = two - (n_p(p) + n_m(q) + n_m(p) + n_p(q))
!                          frac12*(two - erdm(rdm1_p, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
!                                - erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_p, q, q, IAux, NI))
                    end do
              end if
            end associate

      end subroutine pperpa_incore_opshell


      subroutine pperpa_incore_opshell_aaaa(MxA, MxS, AuxData, h, TwoEl, Ndim1, &
            Ndim2, indn1, indn2, indx, multiply_by_S)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(in) :: h(:,:)
            integer, intent(in) :: NDim1, Ndim2
            integer, dimension(:,:), intent(in) ::indn1, indn2
            integer, dimension(:), intent(in) ::indx
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), intent(inout) :: MxA, MxS

            integer :: i, j, t, u, v
            integer p, q, r, s, rs, pq
            double precision :: Arspq, Arsqp
            double precision :: P0, P1, P2, P3, P4
            double precision :: P1a, P1b
            double precision :: P3a1, P3a2, P3b1, P3b2, P3c12, P3c34
            double precision :: P45a1, P45a2, P45b1, P45b2

            double precision :: D0, D1, D2, D3, D4
            double precision :: P3a3, P3a4, P3b3, P3b4, P3c56, P3c78
            double precision :: P45a3, P45a4, P45b3, P45b4
            double precision :: num_f, num_h, Arssum
            double precision :: nnn
            double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
            type (tclock) :: timer0, timer1, timerxx

            associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, &
                  rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m)


              MxA = zero
              MxS = zero

              call clock_start(timerxx)
              i_rowloops: do i = 1, NDim1
                    j_colloops: do j = 1, NDim2
                          call clock_start(timer0)
                          call clock_start(timer1)                          
                          r = indn1(1, i)
                          s = indn1(2, i)
                          rs = indx(i)

                          p = indn2(1, j)
                          q = indn2(2, j)
                          pq = indx(j)

                          Arspq = zero
                          Arsqp = zero
                          P1 = zero; D1 = zero; P2 = zero; D2 = zero

                          P0 =   twoel(gmap(p,r,q,s))
                          D0 = - twoel(gmap(q,r,p,s)) 
                          
                          P0 = P0 - (n_p(p) + n_p(q) + n_p(r)  + n_p(s))* twoel(gmap(p,r,q,s)) 
                          D0 = D0 + (n_p(p) + n_p(q) + n_p(r)  + n_p(s))* twoel(gmap(q,r,p,s)) 


                          P1a = zero
                          P1b = zero
                          if (q == s) then
                                P1a = P1a + h(p, r) 
                                P1b = P1b + h(p, r) *( - n_p(q) -frac12*( n_p(p) + n_p(r) ) )

                                !P1 = P1 + h(p, r) * (one - n_p(q) -frac12*( n_p(p) + n_p(r) ) )
                                
                          end if

                          if (p == r) then
                                P1a = P1a +  h(q, s) 
                                P1b = P1b +  h(q, s) * ( - n_p(p) -frac12 * (n_p(q)  + n_p(s)))

                                !P1 = P1 +  h(q, s) * (one - n_p(p) -frac12 * (n_p(q)  + n_p(s)))
                          end if
                          P1 = P1a+P1b
                          if (p == s) then
                                D1 = D1 - h(q, r) * (one - n_p(p) -frac12 * ( n_p(q) + n_p(r)) )
                          end if

                          if (q == r) then
                                D1 = D1 - h(p, s) * (one - n_p(q) -frac12 * ( n_p(p)  + n_p(s))) 
                          end if

                          P2 = zero; D2 = zero
                          if (q == s) then                                
                                do t = 1, NIA
                                      P2 = P2 +(n_p(t) + n_m(t)) * twoel(gmap(p, r, t, t))- n_p(t) * twoel(gmap(p, t, r, t))
                                end do
                          end if
                          
                          if (p == r) then
                                do t = 1, NIA
                                      P2 = P2 +(n_p(t) + n_m(t)) * twoel(gmap(q, s, t, t)) - n_p(t) *twoel(gmap(q, t, s, t)) 
                                end do
                          end if

                          if (p == s) then
                                do t = 1, NIA
                                      D2 = D2 - (n_p(t) + n_m(t)) * twoel(gmap(q, r, t, t)) + n_p(t) * twoel(gmap(q, t, r, t))
                                end do
                          end if
                          if (q == r) then
                                do t = 1, NIA
                                      D2 = D2 -  (n_p(t) + n_m(t)) * twoel(gmap(p, s, t, t)) + n_p(t)*twoel(gmap(p, t, s, t))
                                end do
                          end if


                          P3a1 = zero; P3a2 = zero; P3b1 = zero; P3b2 = zero; P3c12 = zero; P3c34 = zero
                          P3a3 = zero; P3a4 = zero; P3b3 = zero; P3b4 = zero; P3c56 = zero; P3c78 = zero
                          P45a1 = zero; P45a2 = zero; P45b1 = zero; P45b2 = zero
                          P45a3 = zero; P45a4 = zero; P45b3 = zero; P45b4 = zero

                          P3b1 = P3b1 + func_P3b_aa(AuxData, TwoEl, s, q, p, r)
                          P3b2 = P3b2 + func_P3b_aa(AuxData, TwoEl, r, p, q, s)
                          P3b3 = P3b3 - func_P3b_aa(AuxData, TwoEl, s, p, q, r)
                          P3b4 = P3b4 - func_P3b_aa(AuxData, TwoEl, r, q, p, s)

                          P3c12 = P3c12 - func_P3c_aa(AuxData, TwoEl, s, q, p, r)
                          P3c34 = P3c34 - func_P3c_aa(AuxData, TwoEl, r, p, q, s)
                          P3c56 = P3c56 + func_P3c_aa(AuxData, TwoEl, s, p, q, r)
                          P3c78 = P3c78 + func_P3c_aa(AuxData, TwoEl, r, q, p, s)

                          if (p == r) then
                                P45a1 = P45a1 - func_P45_aa(AuxData, TwoEl, q, s)
                                P45b1 = P45b1 - func_P45_aa(AuxData, TwoEl, s, q)
                          end if

                          if (q == s) then
                                P45a2 = P45a2 - func_P45_aa(AuxData, TwoEl, p, r)
                                P45b2 = P45b2 - func_P45_aa(AuxData, TwoEl, r, p)
                          end if

                          if (q == r) then
                                P45a3 = P45a3 + func_P45_aa(AuxData, TwoEl, p, s)
                                P45b3 = P45b3 + func_P45_aa(AuxData, TwoEl, s, p)
                          end if

                          if (p == s) then
                                P45a4 = P45a4 + func_P45_aa(AuxData, TwoEl, q, r)
                                P45b4 = P45b4 + func_P45_aa(AuxData, TwoEl, r, q)
                          end if

                          ! P2 = two*P2
                          ! P3b1 = two * P3b1
                          ! P3c12 = two*P3c12
                          ! P45a2 = two * P45a2
                          ! P45b2 = two * P45b2

                          Arspq = P0 + P1 + P2 + P3b1 + P3b2 + P3c12 + P3c34 + P45a1+ P45a2 + P45b1 + P45b2 + &
                                D0 + D1 + D2 +  P3b3 + P3b4 + P3c56 + P3c78 + P45a3 + P45a4 + P45b3 + P45b4

                          num_f = one
                          
                          if (multiply_by_S==1)then
                                ! if (abs(one - n_p(r) - n_p(s)).lt.1.e-5)then
                                !     num_f = zero
                                ! else
                                !num_f = one / (one - frac12*(n_p(r) + n_p(s)))
                                !num_f = one / (one - two*(n_p(r) + n_p(s)))
                                num_f = one / (one - (n_p(r) + n_p(s))) ! this is correct. look at your note pp_triplet in notebook

!                                end if
                          end if
                          
                          if (p.ne.q.and.r.ne.s)then
                                Arssum = Arspq 
                                MxA(rs, pq) = num_f * Arssum 
                                if (abs(MxA(rs, pq)).gt.1.d-5)then
                                      !write(*,'(A10, 2I5, A5, 4I5,2F30.16)') 'zulzul', rs, pq, '  |  ', p, q, r, s, MxA(rs, pq)
                                      ! write(*, '(A5, 6I5, 3F15.5)') 'gow', r, s, p, q, rs, pq, (one - (n_p(r) + n_p(s))), MxA(rs, pq), num_f

                                      ! write(*, '(A5, 4A3, 24A12)')   '     ', ' r', ' s', ' p', ' q', &
                                      !       '          P0', '          P1', '         P1a', '         P1b', '          P2', &
                                      !       '       P3c12', '       P3c34', '       P45a1', '       P45a2', '       P45b1', '       P45b2', &
                                      !       '          D0', '          D1', '          D2', &
                                      !       '       P3c56', '       P3c78', '       P45a3', '       P45a4', '       P45b3', '       P45b4', &
                                      !       '        P3b1', '        P3b2', '        P3b3', '        P3b4'
                                      ! write(*, '(A5, 4I3, 24F12.4)') 'plaaaa', r, s, p, q, P0, P1, P1a, P1b, P2, P3c12, P3c34, P45a1, P45a2, P45b1,&
                                      !       P45b2, D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4

                                      ! write(*, '(A5, 4A3, 24A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '     P1a', '     P1b', '      P2',&
                                      !       '   P3c12', '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56',&
                                      !       '   P3c78', '  P45a3', '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4'
                                      ! write(*, '(A5, 4I3, 24F8.4)') 'plaaaa', r, s, p, q, P0, P1, P1a, P1b, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, D0, D1, D2, P3c56,&
                                      !       P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4


                                 end if
                                
                          end if
                                
                    end do j_colloops
              end do i_rowloops

              MxS = zero

              if (multiply_by_s==1)then
                    do i = 1, NDim1
                          p = IndN1(1, I)
                          q = IndN1(2, i)
                          MxS(i, i) = (one - (n_p(p) + n_p(q)))
                          !if (abs(MxS(i,i)).lt.1.d-5) MxS(i,i) = one
!                          frac12*(two - erdm(rdm1_p, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
!                                - erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_p, q, q, IAux, NI))
                    end do
              end if
            end associate

      end subroutine pperpa_incore_opshell_aaaa

      subroutine pperpa_incore_opshell_bbbb(MxA, MxS, AuxData, h, TwoEl, Ndim1, &
            Ndim2, indn1, indn2, indx, multiply_by_S)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(in) :: h(:,:)
            integer, intent(in) :: NDim1, Ndim2
            integer, dimension(:,:), intent(in) ::indn1, indn2
            integer, dimension(:), intent(in) ::indx
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), intent(inout) :: MxA, MxS

            integer :: i, j, t, u, v
            integer p, q, r, s, rs, pq
            double precision :: Arspq, Arsqp
            double precision :: P0, P1, P2, P3, P4
            double precision :: P3a1, P3a2, P3b1, P3b2, P3c12, P3c34
            double precision :: P45a1, P45a2, P45b1, P45b2

            double precision :: D0, D1, D2, D3, D4
            double precision :: P3a3, P3a4, P3b3, P3b4, P3c56, P3c78
            double precision :: P45a3, P45a4, P45b3, P45b4
            double precision :: num_f, num_h, Arssum
            double precision :: nnn
            double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
            type (tclock) :: timer0, timer1, timerxx

            associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
                  NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, &
                  rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m)


              MxA = zero
              MxS = zero

              call clock_start(timerxx)
              i_rowloops: do i = 1, NDim1
                    j_colloops: do j = 1, NDim2
                          call clock_start(timer0)
                          call clock_start(timer1)                          
                          r = indn1(1, i)
                          s = indn1(2, i)
                          rs = indx(i)

                          p = indn2(1, j)
                          q = indn2(2, j)
                          pq = indx(j)

                          Arspq = zero
                          Arsqp = zero
                          P1 = zero; D1 = zero; P2 = zero; D2 = zero

                          P0 =   twoel(gmap(p,r,q,s))
                          D0 = - twoel(gmap(q,r,p,s)) 
                          
                          P0 = P0 - (n_m(p) + n_m(q) + n_m(r)  + n_m(s))* twoel(gmap(p,r,q,s)) 
                          D0 = D0 + (n_m(p) + n_m(q) + n_m(r)  + n_m(s))* twoel(gmap(q,r,p,s)) 

                          if (q == s) then
                                P1 = P1 + h(p, r) * (one - n_m(q) -frac12*( n_m(p) + n_m(r) ) ) 
                          end if

                          if (p == r) then
                                P1 = P1 +  h(q, s) * (one - n_m(p) -frac12 * (n_m(q)  + n_m(s)))
                          end if

                          if (p == s) then
                                D1 = D1 - h(q, r) * (one - n_m(p) -frac12 * ( n_m(q) + n_m(r)) )
                          end if

                          if (q == r) then
                                D1 = D1 - h(p, s) * (one - n_m(q) -frac12 * ( n_m(p)  + n_m(s))) 
                          end if

                          P2 = zero; D2 = zero
                          if (q == s) then                                
                                do t = 1, NIA
                                      P2 = P2 +(n_p(t) + n_m(t)) * twoel(gmap(p, r, t, t))- n_m(t) * twoel(gmap(p, t, r, t))
                                end do
                          end if
                          
                          if (p == r) then
                                do t = 1, NIA
                                      P2 = P2 +(n_p(t) + n_m(t)) * twoel(gmap(q, s, t, t)) - n_m(t) *twoel(gmap(q, t, s, t)) 
                                end do
                          end if

                          if (p == s) then
                                do t = 1, NIA
                                      D2 = D2 - (n_p(t) + n_m(t)) * twoel(gmap(q, r, t, t)) + n_m(t) * twoel(gmap(q, t, r, t)) 
                                end do
                          end if
                          if (q == r) then
                                do t = 1, NIA
                                      D2 = D2 -  (n_p(t) + n_m(t)) * twoel(gmap(p, s, t, t)) + n_m(t)*twoel(gmap(p, t, s, t))
                                end do
                          end if


                          P3a1 = zero; P3a2 = zero; P3b1 = zero; P3b2 = zero; P3c12 = zero; P3c34 = zero
                          P3a3 = zero; P3a4 = zero; P3b3 = zero; P3b4 = zero; P3c56 = zero; P3c78 = zero
                          P45a1 = zero; P45a2 = zero; P45b1 = zero; P45b2 = zero
                          P45a3 = zero; P45a4 = zero; P45b3 = zero; P45b4 = zero

                          
                          !P3a1 = P3a1 + func_P3a_bb(AuxData, TwoEl, s, q, p, r)
                          !P3a2 = P3a2 + func_P3a_bb(AuxData, TwoEl, r, p, q, s)
                          !P3a3 = P3a3 - func_P3a_bb(AuxData, TwoEl, s, p, q, r)
                          !P3a4 = P3a4 - func_P3a_bb(AuxData, TwoEl, r, q, p, s)

                          P3b1 = P3b1 + func_P3b_bb(AuxData, TwoEl, s, q, p, r)
                          P3b2 = P3b2 + func_P3b_bb(AuxData, TwoEl, r, p, q, s)
                          P3b3 = P3b3 - func_P3b_bb(AuxData, TwoEl, s, p, q, r)
                          P3b4 = P3b4 - func_P3b_bb(AuxData, TwoEl, r, q, p, s)


                          P3c12 = P3c12 - func_P3c_bb(AuxData, TwoEl, s, q, p, r)
                          P3c34 = P3c34 - func_P3c_bb(AuxData, TwoEl, r, p, q, s)
                          P3c56 = P3c56 + func_P3c_bb(AuxData, TwoEl, s, p, q, r)
                          P3c78 = P3c78 + func_P3c_bb(AuxData, TwoEl, r, q, p, s)

                          if (p == r) then
                                P45a1 = P45a1 - func_P45_bb(AuxData, TwoEl, q, s)
                                P45b1 = P45b1 - func_P45_bb(AuxData, TwoEl, s, q)
                          end if

                          if (q == s) then
                                P45a2 = P45a2 - func_P45_bb(AuxData, TwoEl, p, r)
                                P45b2 = P45b2 - func_P45_bb(AuxData, TwoEl, r, p)
                          end if

                          if (q == r) then
                                P45a3 = P45a3 + func_P45_bb(AuxData, TwoEl, p, s)
                                P45b3 = P45b3 + func_P45_bb(AuxData, TwoEl, s, p)
                          end if

                          if (p == s) then

                                ! if (r==18.and.s==5.and.p==5.and.q==4)then
                                !       print*, 'takowoz1', P45b4
                                ! end if
                                P45a4 = P45a4 + func_P45_bb(AuxData, TwoEl, q, r)
                                P45b4 = P45b4 + func_P45_bb(AuxData, TwoEl, r, q)
                                ! if (r==18.and.s==5.and.p==5.and.q==4)then
                                !       print*, 'takowoz2', P45b4
                                ! end if

                          end if



                          Arspq = P0 + P1 + P2 + P3b1+P3b2 + P3c12 + P3c34 + P45a1+ P45a2 + P45b1 + P45b2 + &
                                D0 + D1 + D2 + P3b3 + P3b4+ P3c56 + P3c78 + P45a3 + P45a4 + P45b3 + P45b4

                          ! if (abs(P0) > 1.d-8 .or. abs(P1) > 1.d-8 .or. abs(P2) > 1.d-8 .or. &
                          !       abs(P3c12) > 1.d-8 .or. abs(P3c34) > 1.d-8 .or. &
                          !       abs(P45a1) > 1.d-8 .or. abs(P45a2) > 1.d-8 .or. &
                          !       abs(P45b1) > 1.d-8 .or. abs(P45b2) > 1.d-8) then
                          ! if (abs(D0) > 1.d-8a .or. abs(D1) > 1.d-8 .or. abs(D2) > 1.d-8 .or. &
                          !       abs(P3c56) > 1.d-8 .or. abs(P3c78) > 1.d-8 .or. &
                          !       abs(P45a3) > 1.d-8 .or. abs(P45a4) > 1.d-8 .or. &
                          !       abs(P45b3) > 1.d-8 .or. abs(P45b4) > 1.d-8) then
 !                         if (abs(P3a1) > 1.d-8 .or. abs(P3a2) > 1.d-8 .or. abs(P3a3) > 1.d-8 .or. &
  !                              abs(P3a4) > 1.d-8) then

                          ! write(*, '(A5, 4I5, 9F12.6)')'plbbbb', r, s, p,q, P0 , P1 , P2 , P3c12 , P3c34 , P45a1, P45a2 , P45b1 , P45b2
                          ! write(*, '(A5, 4I5, 9F12.6)')'plbbbb', r, s, p, q, D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4
                          ! write(*, '(A5, 4I5, 9F12.6)')'plbbbb', r, s, p, q, P3a1, P3a2, P3a3, P3a4
!                          end if
                          
                          num_f = one

                          if (multiply_by_S==1)then
                                ! if (abs(one - n_m(r) - n_m(s)).lt.1.e-5)then
                                !       write(*,'(A6, 4I5, 3F12.6)') 'takowoz', r, s, p, q,  n_m(r) , n_m(s), one - n_m(r) - n_m(s)
                                !       num_f = zero
                                ! else
                                !num_f = one / (one - frac12 * (n_m(r) + n_m(s)))
                                num_f = one / (one -  (n_m(r) + n_m(s)))
!                                end if
                          end if
                          
                          if (p.ne.q.and.r.ne.s)then
                                Arssum = Arspq 
                                MxA(rs, pq) = num_f * Arssum 
                                if (abs(MxA(rs, pq)).gt.1.d-5)then
                                      ! write(*, '(A5, 6I5, 3F15.5)') 'gow', r, s, p, q, rs, pq, (one -  (n_m(r) + n_m(s))), MxA(rs, pq), num_f
                                      ! !                                end if
                                      ! !                               if (r==18.and.s==5.and.p==5.and.q==4)then
                                      ! write(*, '(A5, 4A3, 22A12)')   '     ', ' r', ' s', ' p', ' q', &
                                      !       '          P0', '          P1', '          P2', &
                                      !       '       P3c12', '       P3c34', '       P45a1', '       P45a2', '       P45b1', '       P45b2', &
                                      !       '          D0', '          D1', '          D2', &
                                      !       '       P3c56', '       P3c78', '       P45a3', '       P45a4', '       P45b3', '       P45b4', &
                                      !       '        P3b1', '        P3b2', '        P3b3', '        P3b4'
                                      ! write(*, '(A5, 4I3, 22F12.4)') 'plbbbb', r, s, p, q, P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                                      !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4
                                      ! write(*, '(A6,  4F15.5)') 'gow2cz', n_m(r), n_m(s),   n_m(r) + n_m(s), one -n_m(r) - n_m(s)
                                      ! write(*, '(A5, 4A3, 22A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '      P2', '   P3c12', '   P3c34',&
                                      !       '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', '  P45a4',&
                                      !       '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4'
                                      ! write(*, '(A5, 4I3, 22F8.4)') 'plbbbb', r, s, p, q, P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4

                                end if
                                
                          end if
                                
                    end do j_colloops
              end do i_rowloops

              MxS = zero

              if (multiply_by_s==1)then
                    do i = 1, NDim1
                          p = IndN1(1, I)
                          q = IndN1(2, i)
                          !MxS(i, i) = one - frac12 * (n_m(p) + n_m(q))
                          MxS(i, i) = one - (n_m(p) + n_m(q))
                          !                          frac12*(two - erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
                          !                                - erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI))
                    end do
              end if
            end associate
      end subroutine pperpa_incore_opshell_bbbb

      subroutine pperpa_incore_opshell_mix(MxA, MxS, AuxData, h, TwoEl, Ndim1, Ndim2, indn1, indn2, indx, multiply_by_S, mix_case)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(in) :: h(:,:)
            integer, intent(in) :: NDim1, Ndim2
            integer, dimension(:,:), intent(in) ::indn1, indn2
            integer, dimension(:), intent(in) ::indx
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), intent(inout) :: MxA, MxS
            integer, optional, intent(in) :: mix_case
            integer :: mode
            double precision :: a, b, c, d


            integer :: i, j, t, u, v
            integer :: p, q, r, s, rs, pq, rs_shift, pq_shift

            double precision :: P01, P11, P21, P22, P3a1, P3b1, P3c12, P45a1, P45a2, P45b1, P45b2
            double precision :: P02, P12, P23, P24, P3a2, P3b2, P3c34, P45a1_2, P45a2_2, P45b1_2, P45b2_2
            double precision :: D01, D11, D21, D22, P3a3, P3b3, P3c56, P45a3, P45a4, P45b3, P45b4
            double precision :: D02, D12, D23, D24, P3a4, P3b4, P3c78, P45a3_2, P45a4_2, P45b3_2, P45b4_2
            double precision :: P0, D0, P1, D1, P2, D2 

            double precision :: num_f_pm, num_f_mp, num_h, val
            double precision :: Arspq_pm_pm, Arspq_mp_mp, Arspq_pm_mp, Arspq_mp_pm

            associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m)

              MxA = zero

              mode = 0
              if (present(mix_case)) mode = mix_case
              
!              print*, 'ndim1', ndim1, ndim2
              i_rowloops: do i = 1, NDim1
                    j_colloops: do j = 1, NDim2
                          r = indn1(1, i)
                          s = indn1(2, i)
                          rs = indx(i)
                          rs_shift = indx(i) + NDim1

                          p = indn2(1, j)
                          q = indn2(2, j)
                          pq = indx(j)
                          pq_shift = indx(j) + NDim2

                          P01 = (one - n_p(p) - n_m(q) - n_p(r) - n_m(s)) * twoel(gmap(p,r,q,s))
                          P02 = (one - n_m(p) - n_p(q) - n_m(r) - n_p(s)) * twoel(gmap(p,r,q,s))
                          D01 = - (one - n_m(p) - n_p(q) - n_p(r) - n_m(s)) * twoel(gmap(q,r,p,s))
                          D02 = - (one - n_p(p) - n_m(q) - n_m(r) - n_p(s)) * twoel(gmap(q,r,p,s))


                          P11 = zero; P12 = zero; D11 = zero; D12 = zero
                          
                          if (q == s) then
                                P11 = P11 + h(p, r) * (one - n_m(q) - frac12*n_p(p) - frac12*n_p(r))
                                P12 = P12 + h(p, r) * (one - n_p(q) - frac12*n_m(p) - frac12*n_m(r))
                          end if
                          if (p == r) then
                                P11 = P11 + h(q, s) * (one - n_p(p) - frac12*n_m(q) - frac12*n_m(s))
                                P12 = P12 + h(q, s) * (one - n_m(p) - frac12*n_p(q) - frac12*n_p(s))
                          end if

                          if (p == s) then
                                D11 = D11 - h(q, r) * (one - n_m(p) - frac12*n_p(q) - frac12*n_p(r))
                                D12 = D12 - h(q, r) * (one - n_p(p) - frac12*n_m(q) - frac12*n_m(r))
                          end if

                          if (q == r) then
                                D11 = D11 - h(p, s) * (one - n_p(q) - frac12*n_m(p) - frac12*n_m(s))
                                D12 = D12 - h(p, s) * (one - n_m(q) - frac12*n_p(p) - frac12*n_p(s))
                          end if

                          ! P2 & D2
                          P21 = zero; P22 = zero; P23 = zero; P24 = zero
                          D21 = zero; D22 = zero; D23 = zero; D24 = zero
                          if (q == s) then
                                do t = 1, NIA
                                      P21 = P21 +  (n_p(t) + n_m(t)) * twoel(gmap(p, r, t, t)) - n_p(t) * twoel(gmap(p, t, r, t))
                                      P23 = P23 +  (n_p(t) + n_m(t)) * twoel(gmap(p, r, t, t)) - n_m(t) * twoel(gmap(p, t, r, t))
                                end do
                          end if

                          if (p == r) then
                                do t = 1, NIA
                                      P22 = P22 +  (n_p(t) + n_m(t)) * twoel(gmap(q, s, t, t)) - n_m(t) * twoel(gmap(q, t, s, t))
                                      P24 = P24 +  (n_p(t) + n_m(t)) * twoel(gmap(q, s, t, t)) - n_p(t) * twoel(gmap(q, t, s, t))
                                end do
                          end if

                          if (q == r) then
                                do t = 1, NIA
                                      D21 = D21 -  (n_p(t) + n_m(t)) * twoel(gmap(p, s, t, t)) + n_m(t) * twoel(gmap(p, t, s, t))
                                      D23 = D23 -  (n_p(t) + n_m(t)) * twoel(gmap(p, s, t, t)) + n_p(t) * twoel(gmap(p, t, s, t))
                                end do
                          end if

                          if (p == s) then
                                do t = 1, NIA
                                      D22 = D22 -  (n_p(t) + n_m(t)) * twoel(gmap(q, r, t, t)) + n_p(t) * twoel(gmap(q, t, r, t))
                                      D24 = D24 -  (n_p(t) + n_m(t)) * twoel(gmap(q, r, t, t)) + n_m(t) * twoel(gmap(q, t, r, t))
                                end do
                          end if

                          ! P3
                          P3a1 = func_P3a_pm(AuxData, TwoEl, s, p, q, r) + func_P3a_mp(AuxData, TwoEl, r, q, p, s)
                          P3a2 = func_P3a_mp(AuxData, TwoEl, s, p, q, r) + func_P3a_pm(AuxData, TwoEl, r, q, p, s)
                          P3a3 = - (func_P3a_pm(AuxData, TwoEl, s, q, p, r) + func_P3a_mp(AuxData, TwoEl, r, p, q, s))
                          P3a4 = - (func_P3a_mp(AuxData, TwoEl, s, q, p, r) + func_P3a_pm(AuxData, TwoEl, r, p, q, s))

                          P3b1 = func_P3b_pm(AuxData, TwoEl, s, q, p, r) + func_P3b_mp(AuxData, TwoEl, r, p, q, s)
                          P3b2 = func_P3b_mp(AuxData, TwoEl, s, q, p, r) + func_P3b_pm(AuxData, TwoEl, r, p, q, s)
                          P3b3 = - (func_P3b_pm(AuxData, TwoEl, s, p, q, r) + func_P3b_mp(AuxData, TwoEl, r, q, p, s))
                          P3b4 = - (func_P3b_mp(AuxData, TwoEl, s, p, q, r) + func_P3b_pm(AuxData, TwoEl, r, q, p, s))

                          P3c12 = - (func_P3c_mm_pm(AuxData, TwoEl, s, q, p, r) + func_P3c_pp_mp(AuxData, TwoEl, r, p, q, s))
                          P3c56 =    func_P3c_mm_pm(AuxData, TwoEl, s, p, q, r) + func_P3c_pp_mp(AuxData, TwoEl, r, q, p, s)
                          P3c34 = - (func_P3c_pp_mp(AuxData, TwoEl, s, q, p, r) + func_P3c_mm_pm(AuxData, TwoEl, r, p, q, s))
                          P3c78 =    func_P3c_pp_mp(AuxData, TwoEl, s, p, q, r) + func_P3c_mm_pm(AuxData, TwoEl, r, q, p, s)

                          ! P45
                          P45a1 = zero; P45a2 = zero; P45b1 = zero; P45b2 = zero
                          P45a1_2 = zero; P45a2_2 = zero; P45b1_2 = zero; P45b2_2 = zero
                          if (p == r) then
                                P45a1   = - func_P45_mp(AuxData, TwoEl, q, s)
                                P45b1   = - func_P45_mp(AuxData, TwoEl, s, q)
                                P45a1_2 = - func_P45_pm(AuxData, TwoEl, q, s)
                                P45b1_2 = - func_P45_pm(AuxData, TwoEl, s, q)
                          end if
                          if (q == s) then
                                P45a2   = - func_P45_pm(AuxData, TwoEl, p, r)
                                P45b2   = - func_P45_pm(AuxData, TwoEl, r, p)
                                P45a2_2 = - func_P45_mp(AuxData, TwoEl, p, r)
                                P45b2_2 = - func_P45_mp(AuxData, TwoEl, r, p)
                          end if

                          P45a3 = zero; P45a4 = zero; P45b3 = zero; P45b4 = zero
                          P45a3_2 = zero; P45a4_2 = zero; P45b3_2 = zero; P45b4_2 = zero
                          if (q == r) then
                                P45a3   = func_P45_mp(AuxData, TwoEl, p, s)
                                P45b3   = func_P45_mp(AuxData, TwoEl, s, p)
                                P45a3_2 = func_P45_pm(AuxData, TwoEl, p, s)
                                P45b3_2 = func_P45_pm(AuxData, TwoEl, s, p)
                          end if
                          if (p == s) then
                                P45a4   = func_P45_pm(AuxData, TwoEl, q, r)
                                P45b4   = func_P45_pm(AuxData, TwoEl, r, q)
                                P45a4_2 = func_P45_mp(AuxData, TwoEl, q, r)
                                P45b4_2 = func_P45_mp(AuxData, TwoEl, r, q)
                          end if

                          ! Assemble blocks
                          Arspq_pm_pm = P01 + P11 + P21 + P22 + P3a1 + P3b1 + P3c12 + P45a1 + P45a2 + P45b1 + P45b2
                          Arspq_mp_mp = P02 + P12 + P23 + P24 + P3a2 + P3b2 + P3c34 + P45a1_2 + P45a2_2 + P45b1_2 + P45b2_2

                          Arspq_pm_mp = D01 + D11 + D21 + D22 + P3a3 + P3b3 + P3c56 + P45a3 + P45a4 + P45b3 + P45b4
                          Arspq_mp_pm = D02 + D12 + D23 + D24 + P3a4 + P3b4 + P3c78 + P45a3_2 + P45a4_2 + P45b3_2 + P45b4_2


!                          Arspq = frac12*(P0 + P1 + P2 + P3a1 + P3a2 + P3b1 + P3b2 + P3c12 + P3c34 + P45a1+ P45a2 + P45b1 + P45b2)
!                          Arspq = frac12*(D0 + D1 + D2 + P3a3 + P3a4 + P3b3 + P3b4 + P3c56 + P3c78 + P45a3 + P45a4 + P45b3 + P45b4)


                          P0 = P01 + P02
                          D0 = D01 + D02
                          P1 = P11 + P12
                          D1 = D11 + D12

                          P2 = P21 + P22 + P23 + P24
                          D2 = D21 + D22 + D23 + D24
                          P45a1 = P45a1 + P45a1_2
                          P45a2 = P45a2 + P45a2_2
                          P45a3 = P45a3 + P45a3_2
                          P45a4 = P45a4 + P45a4_2

                          P45b1 = P45b1 + P45b1_2
                          P45b2 = P45b2 + P45b2_2
                          P45b3 = P45b3 + P45b3_2
                          P45b4 = P45b4 + P45b4_2

!                                                     if (any(abs([P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
!              D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
!              P3a1, P3a2, P3a3, P3a4]) .gt. 1.d-5)) then
!     write(*, '(A5, 4A3, 26A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '      P2', '   P3c12',&
!         '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
!         '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
!     write(*, '(A5, 4I3, 26F8.4)') 'plabab', r, s, p, q, &
!         P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
!         D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
!         P3a1, P3a2, P3a3, P3a4
! end if

                          
                          ! write(*, '(A5, 4A3, 26A8)')   '     ', ' r', ' s', ' p', ' q', '      P0', '      P1', '      P2', '   P3c12',&
                          !       '   P3c34', '  P45a1', '  P45a2', '  P45b1', '  P45b2', '      D0', '      D1', '      D2', '   P3c56', '   P3c78', '  P45a3', &
                          !       '  P45a4', '  P45b3', '  P45b4', '    P3b1', '    P3b2', '    P3b3', '    P3b4', '    P3a1', '    P3a2', '    P3a3', '    P3a4'
                          ! write(*, '(A5, 4I3, 26F8.4)') 'plabab', r, s, p, q, &
                          !       P0, P1, P2, P3c12, P3c34, P45a1, P45a2, P45b1, P45b2, &
                          !       D0, D1, D2, P3c56, P3c78, P45a3, P45a4, P45b3, P45b4, P3b1, P3b2, P3b3, P3b4, &
                          !       P3a1, P3a2, P3a3, P3a4


                          ! write(*, '(A5, 4A3, 44A8)') '     ', ' r', ' s', ' p', ' q', &
                          !       '     P01', '     P11', '     P21', '     P22', '    P3a1', '    P3b1', '   P3c12', '   P45a1', '   P45a2', '   P45b1', '   P45b2', &
                          !       '     P02', '     P12', '     P23', '     P24', '    P3a2', '    P3b2', '   P3c34', ' P45a1_2', ' P45a2_2', ' P45b1_2', ' P45b2_2', &
                          !       '     D01', '     D11', '     D21', '     D22', '    P3a3', '    P3b3', '   P3c56', '   P45a3', '   P45a4', '   P45b3', '   P45b4', &
                          !       '     D02', '     D12', '     D23', '     D24', '    P3a4', '    P3b4', '   P3c78', ' P45a3_2', ' P45a4_2', ' P45b3_2', ' P45b4_2'
                          ! write(*, '(A5, 4I3, 44F8.4)') 'plabab', r, s, p, q, &
                          !       P01, P11, P21, P22, P3a1, P3b1, P3c12, P45a1, P45a2, P45b1, P45b2, &
                          !       P02, P12, P23, P24, P3a2, P3b2, P3c34, P45a1_2, P45a2_2, P45b1_2, P45b2_2, &
                          !       D01, D11, D21, D22, P3a3, P3b3, P3c56, P45a3, P45a4, P45b3, P45b4, &
                          !       D02, D12, D23, D24, P3a4, P3b4, P3c78, P45a3_2, P45a4_2, P45b3_2, P45b4_2

                          num_h = one
                          if (p==q .and. r==s) then
                                num_h = frac12
                          else if ((p==q .and. r/=s) .or. (p/=q .and. r==s)) then
                                num_h = sqrt(frac12)
                          end if

                          num_f_pm = one
                          num_f_mp = one

                          if (multiply_by_S==1)then
                                num_f_pm = one / (one - n_p(r) - n_m(s))
                                num_f_mp = one / (one - n_m(r) - n_p(s))
                          end if

                          a = num_h * num_f_pm * Arspq_pm_pm
                          d = num_h * num_f_mp * Arspq_mp_mp
                          
                          b = num_h * num_f_pm * Arspq_pm_mp
                          c = num_h * num_f_mp * Arspq_mp_pm

                          
                          select case(mode)
                                
                          case(VVOO)
                                MxA(rs, pq)             = frac12 * (a + d - b - c)
                                MxA(rs, pq_shift)       = frac12 * (a - d + b - c)
                                MxA(rs_shift, pq)       = frac12 * (a - d - b + c)
                                MxA(rs_shift, pq_shift) = frac12 * (a + b + c + d)
                          case(VVAO, VVAA)
                                MxA(rs, pq)             = sqrt(frac12) * (a - c)
                                MxA(rs, pq_shift)       = sqrt(frac12) * (b - d)
                                MxA(rs_shift, pq)       = sqrt(frac12) * (a + c)
                                MxA(rs_shift, pq_shift) = sqrt(frac12) * (b + d)

                          case(VAOO, AAOO)
                                MxA(rs, pq)             = sqrt(frac12) * (a - b)
                                MxA(rs, pq_shift)       = sqrt(frac12) * (a + b)
                                MxA(rs_shift, pq)       = sqrt(frac12) * (c - d)
                                MxA(rs_shift, pq_shift) = sqrt(frac12) * (c + d)
                                
                          ! case(VVAO, VVAA)
                          !       MxA(rs, pq)             = sqrt(frac12) * (a + c)
                          !       MxA(rs, pq_shift)       = sqrt(frac12) * (b + d)
                          !       MxA(rs_shift, pq)       = sqrt(frac12) * (a - c)
                          !       MxA(rs_shift, pq_shift) = sqrt(frac12) * (b - d)

                          ! case(VAOO, AAOO)
                          !       MxA(rs, pq)             = sqrt(frac12) * (a + b)
                          !       MxA(rs, pq_shift)       = sqrt(frac12) * (a - b)
                          !       MxA(rs_shift, pq)       = sqrt(frac12) * (c + d)
                          !       MxA(rs_shift, pq_shift) = sqrt(frac12) * (c - d)
                                
                          case default
                                MxA(rs, pq)             = a
                                MxA(rs, pq_shift)       = b
                                MxA(rs_shift, pq)       = c
                                MxA(rs_shift, pq_shift) = d
                          end select


                          ! MxA(rs, pq)             = num_h * num_f_pm * Arspq_pm_pm
                          ! MxA(rs, pq_shift)       = num_h * num_f_pm * Arspq_pm_mp
                          ! MxA(rs_shift, pq)       = num_h * num_f_mp * Arspq_mp_pm
                          ! MxA(rs_shift,pq_shift)  = num_h * num_f_mp * Arspq_mp_mp

                          ! val = (MxA(rs, pq)-MxA(rs, pq_shift)- MxA(rs_shift, pq)+MxA(rs_shift,pq_shift))/two
!                           if (abs(val).gt.1.d-5)then
!                                 !write(*, '(A10, 2I5, 5F15.8)') 'mxamxa', rs, pq, val,MxA(rs, pq), MxA(rs_shift,pq_shift), MxA(rs, pq_shift), MxA(rs_shift, pq)
!                                 write(*, '(A10, 2I5, 3F15.8)') 'mxamxa1', rs, pq, num_h, num_f_pm, Arspq_pm_pm
!                                 write(*, '(A10, 2I5, 3F15.8)') 'mxamxa2', rs, pq_shift, num_h, num_f_pm, Arspq_pm_mp
!                                 write(*, '(A10, 2I5, 3F15.8)') 'mxamxa3', rs_shift, pq, num_h, num_f_pm, Arspq_mp_pm
!                                 write(*, '(A10, 2I5, 3F15.8)') 'mxamxa4', rs_shift, pq_shift, num_h, num_f_pm, Arspq_mp_mp
                                
! !                                write(*, '(A10, 6I5, 8F15.8)') 'mxamxa', rs, pq, r, s, p, q, val,Arspq_pm_pm, Arspq_mp_mp, Arspq_pm_mp, Arspq_mp_pm, num_h, num_f_pm, num_f_mp
!                           end if
                          
                    end do j_colloops
              end do i_rowloops

              MxS = zero
              if (multiply_by_s==1)then
                    do i = 1, NDim1
                          p = IndN1(1, i)
                          q = IndN1(2, i)
                          MxS(indx(i), indx(i))             = (one - n_p(p) - n_m(q))
                          MxS(indx(i)+NDim1, indx(i)+NDim1) = (one - n_m(p) - n_p(q))
                    end do
              end if

            end associate
      end subroutine pperpa_incore_opshell_mix



      subroutine pperpa_incore_opshell_mix2(MxA, MxS, AuxData, h, TwoEl, Ndim1, Ndim2, indn1, indn2, indx, multiply_by_S)
            type(TACppData), intent(inout) :: AuxData
            double precision, intent(in) :: TwoEl(:)
            double precision, intent(in) :: h(:,:)
            integer, intent(in) :: NDim1, Ndim2
            integer, dimension(:,:), intent(in) ::indn1, indn2
            integer, dimension(:), intent(in) ::indx
            integer, intent(in) :: multiply_by_S
            double precision, dimension(:,:), intent(inout) :: MxA, MxS
            double precision :: a, b, c, d


            integer :: i, j, t, u, v
            integer :: p, q, r, s, rs, pq, rs_shift, pq_shift

            double precision :: P01, P11, P21, P22, P3a1, P3b1, P3c12, P45a1, P45a2, P45b1, P45b2
            double precision :: num_f, val
            double precision :: Arspq

            associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, &
                  IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m)

              MxA = zero

!              print*, 'ndim1', ndim1, ndim2
              i_rowloops: do i = 1, NDim1
                    r = IndN1(1, i)
                    s = IndN1(2, i)
                    rs = indx(i)
                    write(*, '(A10, 2I5, 2F10.6)') 'wybu', r, s, n_p(r), n_m(s)
                    j_colloops: do j = 1, NDim2
                          
                          r = IndN1(1, i)
                          s = IndN1(2, i)
                          rs = indx(i)
                          p = IndN2(1, j)
                          q = IndN2(2, j)
                          pq = indx(j)
                          
                          
                          P01 = (one - n_p(p) - n_m(q) - n_p(r) - n_m(s)) * twoel(gmap(p,r,q,s))


                          P11 = zero
                          
                          if (q == s) then
                                P11 = P11 + h(p, r) * (one - n_m(q) - frac12*n_p(p) - frac12*n_p(r))
                          end if
                          if (p == r) then
                                P11 = P11 + h(q, s) * (one - n_p(p) - frac12*n_m(q) - frac12*n_m(s))
                          end if

                          P21 = zero
                          if (q == s) then
                                do t = 1, NIA
                                      P21 = P21 +  (n_p(t) + n_m(t)) * twoel(gmap(p, r, t, t)) - n_p(t) * twoel(gmap(p, t, r, t))
                                end do
                          end if

                          P22 = zero
                          if (p == r) then
                                do t = 1, NIA
                                      P22 = P22 +  (n_p(t) + n_m(t)) * twoel(gmap(q, s, t, t)) - n_m(t) * twoel(gmap(q, t, s, t))
                                end do
                          end if


                          ! P3
                          P3a1 = func_P3a_pm(AuxData, TwoEl, s, p, q, r) + func_P3a_mp(AuxData, TwoEl, r, q, p, s)

                          P3b1 = func_P3b_pm(AuxData, TwoEl, s, q, p, r) + func_P3b_mp(AuxData, TwoEl, r, p, q, s)

                          P3c12 = - (func_P3c_mm_pm(AuxData, TwoEl, s, q, p, r) + func_P3c_pp_mp(AuxData, TwoEl, r, p, q, s))

                          ! P45
                          P45a1 = zero
                          P45a2 = zero
                          P45b1 = zero
                          P45b2 = zero
                          if (p == r) then
                                P45a1   = - func_P45_mp(AuxData, TwoEl, q, s)
                                P45b1   = - func_P45_mp(AuxData, TwoEl, s, q)
                          end if
                          if (q == s) then
                                P45a2   = - func_P45_pm(AuxData, TwoEl, p, r)
                                P45b2   = - func_P45_pm(AuxData, TwoEl, r, p)
                          end if

                          Arspq = P01 + P11 + P21 + P22 + P3a1 + P3b1 + P3c12 + P45a1 + P45a2 + P45b1 + P45b2



                          num_f = one

                          if (multiply_by_S==1)then
                                num_f = one / (one - n_p(r) - n_m(s))
                          end if

                          MxA(rs, pq) = num_f * Arspq
                          ! if (abs(MxA(rs,pq)).gt.1.d-5)then
                          !       write(*, '(A5, 2I3, 12F15.10)') 'mxa', rs, pq, MxA(rs, pq) , P01 , n_p(p) , n_m(q) , n_p(r) , n_m(s)! P11 , P21 , P22 , P3a1 , P3b1 , P3c12 , P45a1 , P45a2 , P45b1 , P45b2
                          !       write(*, '(A5, 2I3, 12F15.10)') 'mxa', rs, pq, MxA(rs, pq) , P01, P11 , P21 , P22 , P3a1 , P3b1 , P3c12 , P45a1 , P45a2 , P45b1 , P45b2
                          ! end if
                          
                    end do j_colloops
              end do i_rowloops

              MxS = zero
              if (multiply_by_s==1)then
                    do i = 1, NDim1
                          p = IndN1(1, i)
                          q = IndN1(2, i)
                          MxS(indx(i), indx(i)) = (one - n_p(p) - n_m(q))
                    end do
              end if

            end associate
      end subroutine pperpa_incore_opshell_mix2


      
      subroutine pperpa_gen_rdm(AuxData, TwoEl)

            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:), intent(in) :: TwoEl
            double precision :: Epp

            integer :: p, q, r, s
            integer :: i, nu
            integer :: NBasis
            integer :: NDim_s, NDim_t, NDim_taa, NDim_tbb

            double precision :: gamma_pp
            double precision :: trans_s, trans_t, trans_taa, trans_tbb
            double precision :: one_el
            double precision :: int_pqrs
            double precision :: fpq, pair_fac

            double precision, dimension(:,:,:), allocatable :: Gam_s
            double precision, dimension(:,:,:), allocatable :: Gam_t
            double precision, dimension(:,:,:), allocatable :: Gam_taa
            double precision, dimension(:,:,:), allocatable :: Gam_tbb

            double precision, dimension(:,:,:,:), allocatable :: G2_pp
            double precision, dimension(:,:), allocatable :: G1_pp
            double precision :: Eone, Etot
            double precision :: sum_g1
            

            associate( &
                  NBasis   => AuxData%NBasis,        &
                  NDim_s   => AuxData%NDim_s,        &
                  NDim_t   => AuxData%NDim_t,        &
                  NDim_taa => AuxData%NDim_t_aa,     &
                  NDim_tbb => AuxData%NDim_t_bb,     &
                  indn_s   => AuxData%IndN_s,        &
                  indn_t   => AuxData%IndN_t,        &
                  indn_taa => AuxData%IndN_t_aa,     &
                  indn_tbb => AuxData%IndN_t_bb,     &
                  eig_s    => AuxData%Eigvec_s,      &
                  eig_t    => AuxData%Eigvec_t,      &
                  eig_taa  => AuxData%Eigvec_t_aa,   &
                  eig_tbb  => AuxData%Eigvec_t_bb,   &
                  vplus_s  => AuxData%vplus_s,       &
                  vplus_t  => AuxData%vplus_t,       &
                  vplus_taa=> AuxData%vplus_t_aa,    &
                  vplus_tbb=> AuxData%vplus_t_bb,    &
                  n_p      => AuxData%n_p,           &
                  n_m      => AuxData%n_m,           &
                  rdm1_p   => AuxData%rdm1_p,        &
                  IAux     => AuxData%IndAux,        &
                  NI       => AuxData%NI )

              allocate(Gam_s(NBasis, NBasis, NDim_s))
              allocate(Gam_t(NBasis, NBasis, NDim_t))
              allocate(Gam_taa(NBasis, NBasis, NDim_taa))
              allocate(Gam_tbb(NBasis, NBasis, NDim_tbb))
              allocate(G2_pp(NBasis, NBasis, NBasis, NBasis))
              allocate(G1_pp(NBasis, NBasis))
              G1_pp = zero

              Gam_s   = zero
              Gam_t   = zero
              Gam_taa = zero
              Gam_tbb = zero
              G2_pp   = zero

              ! ============================================================
              ! Singlet pp channel: spatially symmetric pair function
              ! ============================================================

              do nu = 1, NDim_s
                    if (vplus_s(nu) == 1) then
                          do i = 1, NDim_s

                                p = indn_s(1,i)
                                q = indn_s(2,i)

                                fpq = one
                                if (p == q) fpq = frac12

                                pair_fac = sqrt( &
                                      frac14 * (two - n_p(p) - n_m(q) - n_m(p) - n_p(q)) &
                                      * (two - n_p(p) - n_m(q) - n_m(p) - n_p(q)) )

                                Gam_s(p,q,nu) = pair_fac / sqrt(fpq) * eig_s(i,nu)
                                Gam_s(q,p,nu) = pair_fac / sqrt(fpq) * eig_s(i,nu)

                          end do
                    end if
              end do

              ! ============================================================
              ! Triplet alpha-beta pp channel: spatially antisymmetric
              ! ============================================================

              do nu = 1, NDim_t
                    if (vplus_t(nu) == 1) then
                          do i = 1, NDim_t

                                p = indn_t(1,i)
                                q = indn_t(2,i)

                                pair_fac = sqrt( &
                                      frac14 * (two - n_m(p) - n_p(p) - n_p(q) - n_m(q)) &
                                      * (two - n_m(p) - n_p(p) - n_p(q) - n_m(q)) )

                                Gam_t(p,q,nu) =  pair_fac * eig_t(i,nu)
                                Gam_t(q,p,nu) = -pair_fac * eig_t(i,nu)

                          end do
                    end if
              end do

              ! ============================================================
              ! Triplet alpha-alpha pp channel
              ! ============================================================

              do nu = 1, NDim_taa
                    if (vplus_taa(nu) == 1) then
                          do i = 1, NDim_taa

                                p = indn_taa(1,i)
                                q = indn_taa(2,i)

                                pair_fac = sqrt((one - n_p(p) - n_p(q)) &
                                      * (one - n_p(p) - n_p(q)))

                                Gam_taa(p,q,nu) =  pair_fac * eig_taa(i,nu)
                                Gam_taa(q,p,nu) = -pair_fac * eig_taa(i,nu)

                          end do
                    end if
              end do

              ! ============================================================
              ! Triplet beta-beta pp channel
              ! ============================================================

              do nu = 1, NDim_tbb
                    if (vplus_tbb(nu) == 1) then
                          do i = 1, NDim_tbb

                                p = indn_tbb(1,i)
                                q = indn_tbb(2,i)

                                pair_fac = sqrt((one - n_m(p) - n_m(q)) &
                                      * (one - n_m(p) - n_m(q)))

                                Gam_tbb(p,q,nu) =  pair_fac * eig_tbb(i,nu)
                                Gam_tbb(q,p,nu) = -pair_fac * eig_tbb(i,nu)

                          end do
                    end if
              end do

              ! ============================================================
              ! Build full pp reconstructed 2-RDM
              ! ============================================================

              Epp = zero

              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis

                                      gamma_pp = zero

                                      ! ---------------------------------------------------------
                                      ! one-electron / delta part
                                      ! ---------------------------------------------------------

                                      one_el = zero

                                      if (s == q) one_el = one_el + erdm(rdm1_p, p, r, IAux, NI)
                                      if (q == r) one_el = one_el - erdm(rdm1_p, p, s, IAux, NI)
                                      if (p == s) one_el = one_el - erdm(rdm1_p, q, r, IAux, NI)
                                      if (p == r) one_el = one_el + erdm(rdm1_p, q, s, IAux, NI)

                                      if (s == q .and. p == r) one_el = one_el - one
                                      if (q == r .and. p == s) one_el = one_el + one

                                      gamma_pp = gamma_pp + one_el

                                      ! ---------------------------------------------------------
                                      ! transition part: all pp channels
                                      ! ---------------------------------------------------------

                                      trans_s   = zero
                                      trans_t   = zero
                                      trans_taa = zero
                                      trans_tbb = zero

                                      do nu = 1, NDim_s
                                            trans_s = trans_s + Gam_s(p,q,nu) * Gam_s(r,s,nu)
                                      end do

                                      do nu = 1, NDim_t
                                            trans_t = trans_t + Gam_t(p,q,nu) * Gam_t(r,s,nu)
                                      end do

                                      do nu = 1, NDim_taa
                                            trans_taa = trans_taa + Gam_taa(p,q,nu) * Gam_taa(r,s,nu)
                                      end do

                                      do nu = 1, NDim_tbb
                                            trans_tbb = trans_tbb + Gam_tbb(p,q,nu) * Gam_tbb(r,s,nu)
                                      end do

                                      gamma_pp = gamma_pp + trans_s + trans_t + trans_taa + trans_tbb

                                      G2_pp(p,q,r,s) = gamma_pp

                                      ! Check index convention here: old dump_ppRDM used gmap(r,p,s,q)
                                      int_pqrs = TwoEl(gmap(r,p,s,q))

                                      Epp = Epp + gamma_pp * int_pqrs

                                end do
                          end do
                    end do
              end do

              do q = 1, NBasis
                    do s = 1, NBasis
                          sum_g1 = zero
                          do p = 1, NBasis
                                sum_g1 = sum_g1 + G2_pp(p, q, p, s)
                          end do
                          G1_pp(q, s) = sum_g1 / (AuxData%Nel - one)
                    end do
              end do

              Eone = zero

              do p = 1, NBasis
                    do q = 1, NBasis
                          Eone = Eone + G1_pp(p,q) * AuxData%HNO0(p,q)
                    end do
              end do
              write(*,'(A30,F20.15)') 'Epp from pp reconstructed 2RDM', Epp
              Etot = AuxData%ENuc + Eone + Epp

              write(*,'(A30,F20.15)') 'E_pperpa one-electron energy', Eone
              write(*,'(A30,F20.15)') 'E_pperpa two-electron energy', Epp
              write(*,'(A30,F20.15)') 'E_pperpa total energy', Etot

              deallocate(Gam_s)
              deallocate(Gam_t)
              deallocate(Gam_taa)
              deallocate(Gam_tbb)
              deallocate(G2_pp)
              deallocate(G1_pp)
            end associate

      end subroutine pperpa_gen_rdm


      subroutine dump_ppRDM(AuxData, TwoEl)
            type(TACppData), intent(in) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer :: p, q, r, s, pq, rs, k
            double precision :: sum_s, sum_t, sum_taa, sum_tbb            
            double precision :: fpq, frs, numf, val, Npqrs
            double precision :: rdm1_trace, rdm2_trace_act, rdm2_trace, sum
            double precision, dimension(:,:,:,:), allocatable :: G2_pp, G2_pp_s, G2_pp_t, G2_pp_taa, G2_pp_tbb
            double precision, dimension(:,:), allocatable :: G1_pp, G1_pp_s, G1_pp_t, work
            double precision, dimension(:,:), allocatable :: G1_pp_taa, G1_pp_tbb
            
            integer :: unit

            double precision :: max_diff_herm, max_diff_perm
            double precision :: diff, val_ref, rdm_ref
            double precision :: aux1, contr_s, contr_t, contr_ext, E_pperpa, E_pperpa_s, E_pperpa_t
            double precision :: E_pperpa_taa, E_pperpa_tbb
            double precision :: E_pperpa_one, E_pperpa_one_s, E_pperpa_one_t, E_pperpa_one_taa, E_pperpa_one_tbb
            logical :: cond
            integer :: bad_p, bad_q, bad_r, bad_s, scale
            logical :: sym_ok            

            associate(NA=>AuxData%NA, NI=>AuxData%NI, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, NDim_s=>AuxData%Ndim_s, &
                  NDim_t=>AuxData%Ndim_t,NDim_taa=>AuxData%Ndim_t_aa, NDim_tbb=>AuxData%Ndim_t_bb, &
                  indn_s=>AuxData%IndN_s,  indn_t=>AuxData%IndN_t, indn_taa=>AuxData%IndN_t_aa, indn_tbb=>AuxData%IndN_t_bb, &
                  IAux=>AuxData%IndAux, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, rdm2_pp=>AuxData%rdm2_pp, rdm2_pm=>AuxData%rdm2_pm, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp)

              allocate(G2_pp(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_pp_s(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_pp_t(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_pp_taa(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_pp_tbb(NBasis, NBasis, NBasis, NBasis))

              allocate(G1_pp(NBasis, NBasis))
              allocate(G1_pp_s(NBasis, NBasis))
              allocate(G1_pp_t(NBasis, NBasis))
              allocate(G1_pp_taa(NBasis, NBasis))
              allocate(G1_pp_tbb(NBasis, NBasis))
              
              G2_pp = zero
              G2_pp_s = zero
              G2_pp_t = zero
              G2_pp_taa = zero
              G2_pp_tbb = zero

              contr_ext = zero
              contr_s = zero
              contr_t = zero
              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis
                                      G2_pp(p, q, r, s) = one_el_contr(AuxData, p, q, r, s)
                                      G2_pp_s(p, q, r, s) = one_el_contr(AuxData, p, q, r, s)
                                      !G2_pp_t(p, q, r, s) = one_el_contr(AuxData, p, q, r, s)
                                      
                                end do
                          end do
                    end do
              end do

              contr_s = zero
              contr_t = zero
              do rs = 1, NDim_s
                    do pq = 1, NDim_s

                          r = indn_s(1, rs)
                          s = indn_s(2, rs)

                          p = indn_s(1, pq)
                          q = indn_s(2, pq)

                          fpq = one
                          frs = one
                          if (p==q) fpq = frac12
                          if (r==s) frs = frac12
                          
                          numf = one / sqrt(fpq * frs)
!                          numf = sqrt(fpq * frs)

                          sum_s = zero

                          Npqrs =  frac14*(two-erdm(rdm1_p, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
                                -erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_p, q, q, IAux, NI)) &
                                * (two-erdm(rdm1_p, r, r, IAux, NI)-erdm(rdm1_m, s, s, IAux, NI) &
                                -erdm(rdm1_m, r, r, IAux, NI)-erdm(rdm1_p, s, s, IAux, NI))
                          
                          do k = 1, NDim_s
                                if (AuxData%vplus_s(k) == 1)then                                      
                                      sum_s = sum_s + AuxData%Eigvec_s(pq, k) * AuxData%Eigvec_s(rs, k)
                                end if
                          end do

                          G2_pp(p, q, r, s) = G2_pp(p, q, r, s) + Npqrs * numf *sum_s
                          if (s /= r) then
                                G2_pp(p, q, s, r) = G2_pp(p, q, s, r) + Npqrs * numf *sum_s
                          end if

                          if (q /= p) then
                                G2_pp(q, p, r, s) = G2_pp(q, p, r, s) + Npqrs * numf *sum_s
                          end if
                          
                          if (q /= p .and. s /= r) then
                                G2_pp(q, p, s, r) = G2_pp(q, p, s, r) + Npqrs * numf *sum_s
                          end if


                          Aux1 = (TwoEl(gmap(r,p,s,q)) +  TwoEl(gmap(r,q,s,p)))
                          contr_s = contr_s + Npqrs * numf *sum_s * Aux1
                    end do
              end do
              G2_pp_s = G2_pp

              G2_pp_t = zero
              do rs = 1, NDim_t
                    do pq = 1, NDim_t

                          r = indn_t(1, rs)
                          s = indn_t(2, rs)

                          p = indn_t(1, pq)
                          q = indn_t(2, pq)

                          sum_t = zero
                          do k = 1, NDim_t
                                if (AuxData%vplus_t(k) == 1)then                                      
                                      sum_t = sum_t + AuxData%Eigvec_t(pq, k) * AuxData%Eigvec_t(rs, k)
                                end if
                          end do


                          Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))
                          
                          ! Npqrs =  frac14 * (two-erdm(rdm1_p, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
                          !       -erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_p, q, q, IAux, NI)) &
                          !       * (two-erdm(rdm1_p, r, r, IAux, NI)-erdm(rdm1_m, s, s, IAux, NI) &
                          !       -erdm(rdm1_m, r, r, IAux, NI)-erdm(rdm1_p, s, s, IAux, NI))

!                          if (abs(Npqrs * sum_t).gt.1.d-5)then
!                                if (r==2.and.s==1.and.p==9.and.q==8)then
!                                      write(*, '(A5, 4I3, 2F12.7)')'trll', r, s, p, q, Npqrs, sum_t
!                                end if
 !                         end if
                          sum_t = Npqrs * sum_t

                          !scale  = three
                          scale = one
                          G2_pp(p, q, r, s) = G2_pp(p, q, r, s) + scale * sum_t 
                          G2_pp(p, q, s, r) = G2_pp(p, q, s, r) - scale * sum_t 

                          G2_pp(q, p, r, s) = G2_pp(q, p, r, s) - scale * sum_t 
                          G2_pp(q, p, s, r) = G2_pp(q, p, s, r) + scale * sum_t

                          Aux1 = (TwoEl(gmap(r,p,s,q)) -  TwoEl(gmap(r,q,s,p)))
                          contr_s = contr_s + Npqrs * scale * sum_t * Aux1


                          G2_pp_t(p, q, r, s) = G2_pp_t(p, q, r, s) + scale * sum_t 
                          G2_pp_t(p, q, s, r) = G2_pp_t(p, q, s, r) - scale * sum_t 

                          G2_pp_t(q, p, r, s) = G2_pp_t(q, p, r, s) - scale * sum_t 
                          G2_pp_t(q, p, s, r) = G2_pp_t(q, p, s, r) + scale * sum_t
                          
                    end do
              end do


              G2_pp_taa = zero
              do rs = 1, NDim_taa
                    do pq = 1, NDim_taa

                          r = indn_taa(1, rs)
                          s = indn_taa(2, rs)

                          p = indn_taa(1, pq)
                          q = indn_taa(2, pq)

                          sum_t = zero
                          do k = 1, NDim_taa
                                if (AuxData%vplus_t_aa(k) == 1)then                                      
                                      sum_t = sum_t + AuxData%Eigvec_t_aa(pq, k) * AuxData%Eigvec_t_aa(rs, k)
                                end if
                          end do
                          !Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))                          
                          Npqrs = (one - n_p(p) - n_p(q)) * (one - n_p(r) - n_p(s))
                          sum_t = Npqrs * sum_t
                          
                          G2_pp(p, q, r, s) = G2_pp(p, q, r, s) + sum_t 
                          G2_pp(p, q, s, r) = G2_pp(p, q, s, r) - sum_t 

                          G2_pp(q, p, r, s) = G2_pp(q, p, r, s) - sum_t 
                          G2_pp(q, p, s, r) = G2_pp(q, p, s, r) + sum_t

                          Aux1 = (TwoEl(gmap(r,p,s,q)) -  TwoEl(gmap(r,q,s,p)))
                          contr_s = contr_s + Npqrs  * sum_t * Aux1


                          G2_pp_taa(p, q, r, s) = G2_pp_taa(p, q, r, s) + sum_t 
                          G2_pp_taa(p, q, s, r) = G2_pp_taa(p, q, s, r) - sum_t 

                          G2_pp_taa(q, p, r, s) = G2_pp_taa(q, p, r, s) - sum_t 
                          G2_pp_taa(q, p, s, r) = G2_pp_taa(q, p, s, r) + sum_t
                          
                    end do
              end do


              G2_pp_tbb = zero
              do rs = 1, NDim_tbb
                    do pq = 1, NDim_tbb

                          r = indn_tbb(1, rs)
                          s = indn_tbb(2, rs)

                          p = indn_tbb(1, pq)
                          q = indn_tbb(2, pq)

                          sum_t = zero
                          do k = 1, NDim_tbb
                                if (AuxData%vplus_t_bb(k) == 1)then                                      
                                      sum_t = sum_t + AuxData%Eigvec_t_bb(pq, k) * AuxData%Eigvec_t_bb(rs, k)
                                end if
                          end do

                          Npqrs = (one - n_m(p) - n_m(q)) * (one - n_m(r) - n_m(s))
                          !Npqrs  = frac14 * (two - n_m(p) - n_p(p) -n_p(q) - n_m(q)) * (two - n_m(r)-n_m(s) - n_p(r)-n_p(s))                          

                          sum_t = Npqrs * sum_t
                          
                          G2_pp(p, q, r, s) = G2_pp(p, q, r, s) + sum_t 
                          G2_pp(p, q, s, r) = G2_pp(p, q, s, r) - sum_t 

                          G2_pp(q, p, r, s) = G2_pp(q, p, r, s) - sum_t 
                          G2_pp(q, p, s, r) = G2_pp(q, p, s, r) + sum_t

                          Aux1 = (TwoEl(gmap(r,p,s,q)) -  TwoEl(gmap(r,q,s,p)))
                          contr_s = contr_s + Npqrs * sum_t * Aux1


                          G2_pp_tbb(p, q, r, s) = G2_pp_tbb(p, q, r, s) + sum_t 
                          G2_pp_tbb(p, q, s, r) = G2_pp_tbb(p, q, s, r) - sum_t 

                          G2_pp_tbb(q, p, r, s) = G2_pp_tbb(q, p, r, s) - sum_t 
                          G2_pp_tbb(q, p, s, r) = G2_pp_tbb(q, p, s, r) + sum_t
                          
                    end do
              end do


              E_pperpa = zero
              E_pperpa_s = zero
              E_pperpa_t = zero
              E_pperpa_taa = zero
              E_pperpa_tbb = zero
              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis
                                      E_pperpa = E_pperpa + frac12 * G2_pp(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_pperpa_s = E_pperpa_s + frac12 * G2_pp_s(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_pperpa_t = E_pperpa_t + frac12 * G2_pp_t(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_pperpa_taa = E_pperpa_taa + frac12 * G2_pp_taa(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_pperpa_tbb = E_pperpa_tbb + frac12 * G2_pp_tbb(p, q, r, s) * TwoEl(gmap(p, r, q, s))
!                         if (abs(frac12 * G2_pp_t(p, q, r, s) * TwoEl(gmap(p, r, q, s))).gt.1.d-5)then
!                               write(*, '(A6, 4I3, 3F12.7)') 'zzz', p, q, r, s, G2_pp_t(p, q, r, s) , TwoEl(gmap(p, r, q, s)), E_pperpa_t
!                         end if
                                            
                                end do
                          end do
                    end do
              end do
              print*, 'E_pperpa two-electron energy = ', E_pperpa
              print*, 'E_pperpa sing-two-electron energy = ', E_pperpa_s
              print*, 'E_pperpa trip-two-electron energyab = ', E_pperpa_t
              print*, 'E_pperpa trip-two-electron energyaa = ', E_pperpa_taa
              print*, 'E_pperpa trip-two-electron energybb = ', E_pperpa_tbb


              rdm2_trace_act = zero
              do p = 1, NIA
                    do q = 1, NIA
                          rdm2_trace_act = rdm2_trace_act + G2_pp(p, q, p, q)
                    end do
              end do



              rdm2_trace = zero
              do p = 1, NIA
                    do q = 1, NIA
                          rdm_ref = + erdm_ppx(rdm2_pp, n_p, p, q, p, q, IAux, NI) &
                                + erdm_ppx(rdm2_mm, n_m, p, q, p, q, IAux, NI)&
                                + erdm_pmx(rdm2_pm, n_p, n_m, p, q, p, q, IAux, NI)&
                                + erdm_pmx(rdm2_mp, n_m, n_p, p, q, p, q, IAux, NI)
                          
                          rdm2_trace = rdm2_trace + rdm_ref
                          diff = abs(rdm_ref- G2_pp(p,q,p,q))
                          ! if (abs(diff).gt.1.d-8)then
                          !       write(*, '(A10, 4I3, 9F12.8)') 'niebieski', p, q, p, q, rdm_ref, G2_pp(p,q,p,q), &
                          !             erdm_ppx(rdm2_pp, n_p, p, q, p, q, IAux, NI), &
                          !             erdm_pmx(rdm2_pm, n_p, n_m, p, q, p, q, IAux, NI)
                          ! end if
                    end do
              end do


              rdm2_trace = zero
              do p = 1, NBasis
                 do q = 1, NBasis
                    rdm2_trace = rdm2_trace + G2_pp(p, q, p, q)
                 end do
              end do

              max_diff_herm = zero
              max_diff_perm = zero
              sym_ok = .true.


              G1_pp = zero
              G1_pp_s = zero
              G1_pp_t = zero
              G1_pp_taa = zero
              G1_pp_tbb = zero
              do q = 1, NBasis
                    do s = 1, NBasis
                          sum = zero
                          sum_s = zero
                          sum_t = zero
                          sum_taa = zero
                          sum_tbb = zero
                          do p = 1, NBasis
                                sum = sum + G2_pp(p, q, p, s)
                                sum_s = sum_s + G2_pp_s(p, q, p, s)
                                sum_t = sum_t + G2_pp_t(p, q, p, s)
                                sum_taa = sum_taa + G2_pp_taa(p, q, p, s)
                                sum_tbb = sum_tbb + G2_pp_tbb(p, q, p, s)
                          end do
                          G1_pp(q, s) = sum / (AuxData%Nel-one)
                          G1_pp_s(q, s) = sum_s / (AuxData%Nel-one)
                          G1_pp_t(q, s) = sum_t / (AuxData%Nel-one)
                          G1_pp_taa(q, s) = sum_taa / (AuxData%Nel-one)
                          G1_pp_tbb(q, s) = sum_tbb / (AuxData%Nel-one)
                    end do
              end do

              rdm1_trace = zero
              do q = 1, NBasis
                    rdm1_trace = rdm1_trace + G1_pp(q,q)
              end do

              E_pperpa_one = zero
              E_pperpa_one_s = zero
              E_pperpa_one_t = zero
              E_pperpa_one_taa = zero
              E_pperpa_one_tbb = zero
              do p = 1, NBasis
                    do q = 1, NBasis
                          E_pperpa_one = E_pperpa_one + G1_pp(p, q) * AuxData%HNO0(p,q)
                          E_pperpa_one_s = E_pperpa_one_s + G1_pp_s(p, q) * AuxData%HNO0(p,q)
!                          if (abs(G1_pp_s(p, q) * AuxData%HNO0(p,q)).gt.1.d-5)then
!                                write(*, '(A5, 2I5, 4F12.6)')'pla', p, q, G1_pp_s(p, q) , AuxData%HNO0(p,q), G1_pp_s(p, q) * AuxData%HNO0(p,q), E_pperpa_one_s
 !                         end if
                          E_pperpa_one_t = E_pperpa_one_t + G1_pp_t(p, q) * AuxData%HNO0(p,q)
                          E_pperpa_one_taa = E_pperpa_one_taa + G1_pp_taa(p, q) * AuxData%HNO0(p,q)
                          E_pperpa_one_tbb = E_pperpa_one_tbb + G1_pp_tbb(p, q) * AuxData%HNO0(p,q)
                    end do
              end do
              print*, 'E_pperpa one-electron energy = ', E_pperpa_one
              print*, 'E_pperpa sing-one-electron energy = ', E_pperpa_one_s
              print*, 'E_pperpa trip-one-electron energyab = ', E_pperpa_one_t
              print*, 'E_pperpa trip-one-electron energyaa = ', E_pperpa_one_taa
              print*, 'E_pperpa trip-one-electron energybb = ', E_pperpa_one_tbb
              print*, 'E_pperpa sing-all-electron energy = ', E_pperpa_one_s + E_pperpa_s
              print*, 'E_pperpa trip-all-electron energy = ', E_pperpa_one_t + E_pperpa_one_taa + E_pperpa_one_tbb+ E_pperpa_t+ E_pperpa_taa+ E_pperpa_tbb
              
              print*, 'total energy                   ', AuxData%ENuc + E_pperpa_one + E_pperpa
              print*, 'nuclear energy               = ', AuxData%ENuc

              
              print*, 'ref trace', rdm2_trace
              print*, 'active trace', rdm2_trace_act



              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis
                                      
                                      val_ref = G2_pp(p, q, r, s)
                                      diff = abs(val_ref - G2_pp(r, s, p, q))
                                      if (diff > max_diff_herm) then
                                            max_diff_herm = diff
                                            if (diff > 1.0d-5) then
                                                  bad_p = p; bad_q = q; bad_r = r; bad_s = s
                                            end if
                                      end if
                                      ! 2. Check Permutation (simultaneous swap): (pq|rs) == (qp|sr)
                                      diff = abs(val_ref - G2_pp(q, p, s, r))
                                      if (diff > max_diff_perm) max_diff_perm = diff

                                end do
                          end do
                    end do
              end do


              if (max_diff_herm > 1.0d-5 .or. max_diff_perm > 1.0d-5) then
                    print *, "WARNING: Significant symmetry violation detected!"
                    print *, "Worst element at:", bad_p, bad_q, bad_r, bad_s
              end if

              call save_2rdm_text('PPrdm2_natural.dat', G2_pp)
              
               if (allocated(AuxData%CMONO)) then
                     print*, 'Transforming G2_pp back to Original/MO basis'
                     allocate(work(NBasis, NBasis))
                     work = transpose(AuxData%CMONO)
                     call rdm2_MO_NO_trans(G2_pp, work, NBasis)
                     deallocate(work)
               end if

               call save_2rdm_text('PPrdm2_original.dat', G2_pp)
               call save_2rdm_text('PPERPA_2RDM.dat', G2_pp)


              print '(A)', "------------------------------------------------------------"
              print '(A, I15)',   " Number of electrons:                ", AuxData%Nel
              print '(A, F15.8)', " Reference 2-RDM trace [N*(N-1)]:    ", real(AuxData%Nel * (AuxData%Nel - 1), F64)
              print '(A, F15.8)', " Calculated 2-RDM trace:             ", rdm2_trace
              print '(A, F15.8)', " Calculated 1-RDM trace:             ", rdm1_trace
              print '(A)',        " --- 2-RDM Symmetry Check ---"
              print '(A, E15.4)', " Max Error Hermiticity (pqrs|rspq):  ", max_diff_herm
              print '(A, E15.4)', " Max Error Permutation (pqrs|qpsr):  ", max_diff_perm
              print '(A)',        " ------------------------------------------------------------"
              print '(A, A)',     "2-RDM matrix saved to:               ", 'PPrdm2.dat'
              print '(A)',        " ------------------------------------------------------------"
              print "(A)",        " --- 2-RDM Storage Convention ---"
              print "(A)",        " Stored as: kjli--> rdm2(l,k,i,j) - original MOLMPS convention"
              print "(A)",        " ------------------------------------------------------------"


              
            end associate
            
      end subroutine dump_ppRDM

      ! ============================================================================
      ! dump_hhRDM: Reconstruct the 2-RDM from the hhERPA (hole-hole) branch
      !
      ! Theory (Eq. 41 from 2RDM-ERPA.md):
      !   Gamma_{pqrs} = sum_{nu in H_{N-2}} [gamma_hh^nu]_{pq} [gamma_hh^nu]_{rs}
      !
      ! Key differences from dump_ppRDM:
      !   1. NO one_el_contr (the hhERPA one_el_contr = 0)
      !   2. Use vplus == 0 (hh modes, normalized to -1) instead of vplus == 1
      ! ============================================================================
      subroutine dump_hhRDM(AuxData, TwoEl, z)
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: z
            integer :: p, q, r, s, pq, rs, k
            double precision :: sum_s, sum_t, sum_taa, sum_tbb, sum_flukt
            double precision :: fpq, frs, numf, val, Npqrs
            double precision :: rdm1_trace, rdm1_flukt_trace, rdm2_trace, sum, rdm2_flukt_trace, rdm2_flukt_trace_act
            double precision, dimension(:,:,:,:), allocatable :: G2_hh, G2_hh_s, G2_hh_t, G2_hh_taa, G2_hh_tbb, G2_hh_flukt
            double precision, dimension(:,:), allocatable :: G1_hh, G1_hh_s, G1_hh_t, G1_hh_taa, G1_hh_tbb, work, G1_hh_flukt

            double precision :: max_diff_herm, max_diff_perm
            double precision :: diff, val_ref
            double precision :: aux1, E_hherpa, E_hherpa_ac, E_hherpa_one, E_hherpa_flukt
            double precision :: E_hherpa_s, E_hherpa_t, E_hherpa_taa, E_hherpa_tbb
            double precision :: E_hherpa_one_flukt, E_hherpa_one_s, E_hherpa_one_t, E_hherpa_one_taa, E_hherpa_one_tbb
            integer :: bad_p, bad_q, bad_r, bad_s, scale
            logical :: sym_ok

            associate(NA=>AuxData%NA, NI=>AuxData%NI, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, NDim_s=>AuxData%Ndim_s, &
                  NDim_t=>AuxData%Ndim_t,NDim_taa=>AuxData%Ndim_t_aa, NDim_tbb=>AuxData%Ndim_t_bb, &
                  indn_s=>AuxData%IndN_s,  indn_t=>AuxData%IndN_t, indn_taa=>AuxData%IndN_t_aa, indn_tbb=>AuxData%IndN_t_bb, &
                  IAux=>AuxData%IndAux, &
                  rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, rdm2_pp=>AuxData%rdm2_pp, rdm2_pm=>AuxData%rdm2_pm, &
                  n_p=>AuxData%n_p, n_m=>AuxData%n_m, rdm2_mm=>AuxData%rdm2_mm, rdm2_mp=>AuxData%rdm2_mp)

              allocate(G2_hh(NBasis, NBasis, NBasis, NBasis))

              allocate(G2_hh_s(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_hh_t(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_hh_taa(NBasis, NBasis, NBasis, NBasis))
              allocate(G2_hh_tbb(NBasis, NBasis, NBasis, NBasis))
              allocate(G1_hh(NBasis, NBasis))
              allocate(G1_hh_flukt(NBasis, NBasis))
              allocate(G1_hh_s(NBasis, NBasis))
              allocate(G1_hh_t(NBasis, NBasis))
              allocate(G1_hh_taa(NBasis, NBasis))
              allocate(G1_hh_tbb(NBasis, NBasis))

              ! hhERPA: NO one_el_contr (it is zero by theory, Eq. 41)
              G2_hh = zero
              G2_hh_s = zero
              G2_hh_t = zero
              G2_hh_taa = zero
              G2_hh_tbb = zero

              ! ================================================================
              ! SINGLET BLOCK: use vplus_s == 0 (hh modes)
              ! ================================================================


              do rs = 1, NDim_s
                    do pq = 1, NDim_s

                          r = indn_s(1, rs)
                          s = indn_s(2, rs)

                          p = indn_s(1, pq)
                          q = indn_s(2, pq)

                          fpq = one
                          frs = one
                          if (p==q) fpq = frac12
                          if (r==s) frs = frac12

                          numf = one / sqrt(fpq * frs)

                          sum_s = zero

                          Npqrs =  frac14*(two-erdm(rdm1_p, p, p, IAux, NI)-erdm(rdm1_m, q, q, IAux, NI)&
                                -erdm(rdm1_m, p, p, IAux, NI)-erdm(rdm1_p, q, q, IAux, NI)) &
                                * (two-erdm(rdm1_p, r, r, IAux, NI)-erdm(rdm1_m, s, s, IAux, NI) &
                                -erdm(rdm1_m, r, r, IAux, NI)-erdm(rdm1_p, s, s, IAux, NI))

                          do k = 1, NDim_s
                                if (AuxData%vplus_s(k) == 0)then
                                      sum_s = sum_s + AuxData%Eigvec_s(pq, k) * AuxData%Eigvec_s(rs, k)
!                                      write(*, '(A10, 3I5, 3F15.8)') 'vect', pq, rs, k, AuxData%Eigvec_s(pq, k) , AuxData%Eigvec_s(rs, k), sum
                                end if
                          end do

                          G2_hh(p, q, r, s) = G2_hh(p, q, r, s) + Npqrs * numf * sum_s
!                          write(*, '(A5, 4I5, 4F12.10)') 'w1', p, q, r, s, G2_hh(p, q, r, s), Npqrs , numf , sum_s
                          
                          G2_hh_s(p, q, r, s) = G2_hh_s(p, q, r, s) + Npqrs * numf * sum_s
                          if (s /= r) then
                                G2_hh(p, q, s, r) = G2_hh(p, q, s, r) + Npqrs * numf * sum_s
                                G2_hh_s(p, q, s, r) = G2_hh_s(p, q, s, r) + Npqrs * numf * sum_s
!                                write(*, '(A5, 4I5, 4F12.10)') 'w2', p, q, s, r, G2_hh(p, q, s, r), Npqrs , numf , sum_s
                          end if

                          if (q /= p) then
                                G2_hh(q, p, r, s) = G2_hh(q, p, r, s) + Npqrs * numf * sum_s
                                G2_hh_s(q, p, r, s) = G2_hh_s(q, p, r, s) + Npqrs * numf * sum_s
!                                write(*, '(A5, 4I5, 4F12.10)') 'w3',q, p, r, s, G2_hh(q, p, r, s), Npqrs , numf , sum_s
                          end if

                          if (q /= p .and. s /= r) then
                                G2_hh(q, p, s, r) = G2_hh(q, p, s, r) + Npqrs * numf * sum_s
                                G2_hh_s(q, p, s, r) = G2_hh_s(q, p, s, r) + Npqrs * numf * sum_s
!                                write(*, '(A5, 4I5, 4F12.10)') 'w3', q, p, s, r, G2_hh(q, p, s, r), Npqrs , numf , sum_s
                          end if

                    end do
              end do

              ! ================================================================
              ! TRIPLET ABAB BLOCK: use vplus_t == 0 (hh modes)
              ! ================================================================
              do rs = 1, NDim_t
                    do pq = 1, NDim_t

                          r = indn_t(1, rs)
                          s = indn_t(2, rs)

                          p = indn_t(1, pq)
                          q = indn_t(2, pq)

                          sum_t = zero
                          do k = 1, NDim_t
                                if (AuxData%vplus_t(k) == 0)then
                                      sum_t = sum_t + AuxData%Eigvec_t(pq, k) * AuxData%Eigvec_t(rs, k)
                                end if
                          end do

                          Npqrs = frac14 * (two - n_m(p) - n_p(p) - n_p(q) - n_m(q)) &
                                * (two - n_m(r) - n_m(s) - n_p(r) - n_p(s))
                          sum_t = Npqrs * sum_t

                          G2_hh(p, q, r, s) = G2_hh(p, q, r, s) + sum_t
                          G2_hh_t(p, q, r, s) = G2_hh_t(p, q, r, s) + sum_t
                          G2_hh(p, q, s, r) = G2_hh(p, q, s, r) - sum_t
                          G2_hh_t(p, q, s, r) = G2_hh_t(p, q, s, r) - sum_t

                          G2_hh(q, p, r, s) = G2_hh(q, p, r, s) - sum_t
                          G2_hh_t(q, p, r, s) = G2_hh_t(q, p, r, s) - sum_t
                          G2_hh(q, p, s, r) = G2_hh(q, p, s, r) + sum_t
                          G2_hh_t(q, p, s, r) = G2_hh_t(q, p, s, r) + sum_t

                    end do
              end do

              ! ================================================================
              ! TRIPLET AAAA BLOCK: use vplus_t_aa == 0 (hh modes)
              ! ================================================================
              do rs = 1, NDim_taa
                    do pq = 1, NDim_taa

                          r = indn_taa(1, rs)
                          s = indn_taa(2, rs)

                          p = indn_taa(1, pq)
                          q = indn_taa(2, pq)

                          sum_t = zero
                          do k = 1, NDim_taa
                                if (AuxData%vplus_t_aa(k) == 0)then
                                      sum_t = sum_t + AuxData%Eigvec_t_aa(pq, k) * AuxData%Eigvec_t_aa(rs, k)
                                end if
                          end do
                          Npqrs = (one - n_p(p) - n_p(q)) * (one - n_p(r) - n_p(s))
                          sum_t = Npqrs * sum_t

                          G2_hh(p, q, r, s) = G2_hh(p, q, r, s) + sum_t
                          G2_hh_taa(p, q, r, s) = G2_hh_taa(p, q, r, s) + sum_t
                          G2_hh(p, q, s, r) = G2_hh(p, q, s, r) - sum_t
                          G2_hh_taa(p, q, s, r) = G2_hh_taa(p, q, s, r) - sum_t

                          G2_hh(q, p, r, s) = G2_hh(q, p, r, s) - sum_t
                          G2_hh_taa(q, p, r, s) = G2_hh_taa(q, p, r, s) - sum_t
                          G2_hh(q, p, s, r) = G2_hh(q, p, s, r) + sum_t
                          G2_hh_taa(q, p, s, r) = G2_hh_taa(q, p, s, r) + sum_t

                    end do
              end do

              ! ================================================================
              ! TRIPLET BBBB BLOCK: use vplus_t_bb == 0 (hh modes)
              ! ================================================================
              do rs = 1, NDim_tbb
                    do pq = 1, NDim_tbb

                          r = indn_tbb(1, rs)
                          s = indn_tbb(2, rs)

                          p = indn_tbb(1, pq)
                          q = indn_tbb(2, pq)

                          sum_t = zero
                          do k = 1, NDim_tbb
                                if (AuxData%vplus_t_bb(k) == 0)then
                                      sum_t = sum_t + AuxData%Eigvec_t_bb(pq, k) * AuxData%Eigvec_t_bb(rs, k)
                                end if
                          end do

                          Npqrs = (one - n_m(p) - n_m(q)) * (one - n_m(r) - n_m(s))
                          sum_t = Npqrs * sum_t

                          G2_hh(p, q, r, s) = G2_hh(p, q, r, s) + sum_t
                          G2_hh_tbb(p, q, r, s) = G2_hh_tbb(p, q, r, s) + sum_t
                          G2_hh(p, q, s, r) = G2_hh(p, q, s, r) - sum_t
                          G2_hh_tbb(p, q, s, r) = G2_hh_tbb(p, q, s, r) - sum_t

                          G2_hh(q, p, r, s) = G2_hh(q, p, r, s) - sum_t
                          G2_hh_tbb(q, p, r, s) = G2_hh_tbb(q, p, r, s) - sum_t
                          G2_hh(q, p, s, r) = G2_hh(q, p, s, r) + sum_t
                          G2_hh_tbb(q, p, s, r) = G2_hh_tbb(q, p, s, r) + sum_t
                          ! if (z==1)then
                          !       if (abs(G2_hh(p, q, r, s)).gt.1.d-5)then
                          !             print*, 'gowno', G2_hh(p, q, r, s), AuxData%pp2rdm0(p, q, r, s)
                          !       end if
                          ! end if
                    end do
              end do
              if (z == 0)then
                    allocate(AuxData%pp2rdm0(NBasis, NBasis, NBasis, NBasis))
                    AuxData%pp2rdm0 = G2_hh
              else
                    allocate(G2_hh_flukt(NBasis, NBasis, NBasis, NBasis))
                    G2_hh_flukt = G2_hh - AuxData%pp2rdm0
              end if
              ! ================================================================
              ! Compute 2-electron energy
              ! ================================================================
              E_hherpa = zero
              E_hherpa_flukt = zero
              E_hherpa_ac = zero
              E_hherpa_s = zero
              E_hherpa_t = zero
              E_hherpa_taa = zero
              E_hherpa_tbb = zero
              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis
                                      if (.not.(Iaux(p)==1.and.Iaux(q)==1.and.Iaux(r)==1.and.Iaux(s)==1))then
                                            E_hherpa_ac = E_hherpa_ac + frac12 * G2_hh(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      end if
                                      E_hherpa = E_hherpa + frac12 * G2_hh(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      if (z ==1)then
                                            E_hherpa_flukt = E_hherpa_flukt + frac12 * G2_hh_flukt(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      end if
                                      E_hherpa_s = E_hherpa_s + frac12 * G2_hh_s(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_hherpa_t = E_hherpa_t + frac12 * G2_hh_t(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_hherpa_taa = E_hherpa_taa + frac12 * G2_hh_taa(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                      E_hherpa_tbb = E_hherpa_tbb + frac12 * G2_hh_tbb(p, q, r, s) * TwoEl(gmap(p, r, q, s))
                                end do
                          end do
                    end do
              end do
              print*, 'E_pperpa two-electron energy_ac = ', E_hherpa_ac
              if (z==1)then
                    print*, 'E_pperpa two-electron energy flukt = ', E_hherpa_flukt
              end if
              print*, 'E_pperpa two-electron energy = ', E_hherpa
              print*, 'E_pperpa sing-two-electron energy = ', E_hherpa_s
              print*, 'E_pperpa trip-two-electron energyab = ', E_hherpa_t
              print*, 'E_pperpa trip-two-electron energyaa = ', E_hherpa_taa
              print*, 'E_pperpa trip-two-electron energybb = ', E_hherpa_tbb

              ! ================================================================
              ! Traces and 1-RDM from partial trace
              ! ================================================================
              rdm2_trace = zero
              rdm2_flukt_trace = zero
              do p = 1, NBasis
                 do q = 1, NBasis
                       rdm2_trace = rdm2_trace + G2_hh(p, q, p, q)
                       if (z ==1)then
                             rdm2_flukt_trace = rdm2_flukt_trace + G2_hh_flukt(p, q, p, q)
                       end if
                 end do
           end do

           if (z ==1)then
                 rdm2_flukt_trace_act = zero
                 do p = NI+1, NIA
                       do q = NI+1, NIA
                             rdm2_flukt_trace_act = rdm2_flukt_trace_act + G2_hh_flukt(p, q, p, q)
                       end do
                 end do
                 print*, 'rdm_flukt_trace_act', rdm2_flukt_trace_act
           end if
           

              G1_hh = zero
              G1_hh_flukt = zero
              G1_hh_s = zero
              G1_hh_t = zero
              G1_hh_taa = zero
              G1_hh_tbb = zero
              do q = 1, NBasis
                    do s = 1, NBasis                          
                          sum = zero
                          sum_flukt = zero
                          sum_s = zero
                          sum_t = zero
                          sum_taa = zero
                          sum_tbb = zero
                          do p = 1, NBasis
                                sum = sum + G2_hh(p, q, p, s)
                                if (z==1)then
                                      sum_flukt = sum_flukt + G2_hh_flukt(p, q, p, s)
                                end if
                                sum_s = sum_s + G2_hh_s(p, q, p, s)
!                                if (abs(G2_hh_s(p, q, p, s)).gt.1.d-5)then
!                                      write(*, '(A5, 3I5, 2F12.6)')'plax', q, s, p, G2_hh_s(p, q, p, s), sum!G1_hh_s(p, q) , AuxData%HNO0(p,q), G1_hh_s(p, q) * AuxData%HNO0(p,q), E_hherpa_one_s
!                                end if
                                sum_t = sum_t + G2_hh_t(p, q, p, s)
                                sum_taa = sum_taa + G2_hh_taa(p, q, p, s)
                                sum_tbb = sum_tbb + G2_hh_tbb(p, q, p, s)
                          end do
                          G1_hh(q, s) = sum / (AuxData%Nel-one)
                          if (z==1)then
                                G1_hh_flukt(q, s) = sum_flukt / (AuxData%Nel-one)
                          end if
!                          write(*, '(A5, 2I3, F15.8)') 'ghh', q, s, G1_hh(q,s)
                          G1_hh_s(q, s) = sum_s / (AuxData%Nel-one)
                          G1_hh_t(q, s) = sum_t / (AuxData%Nel-one)
                          G1_hh_taa(q, s) = sum_taa / (AuxData%Nel-one)
                          G1_hh_tbb(q, s) = sum_tbb / (AuxData%Nel-one)
                    end do
              end do

              rdm1_trace = zero
              rdm1_flukt_trace = zero
              do q = 1, NBasis
                    rdm1_trace = rdm1_trace + G1_hh(q,q)
                    if (z==1)then
                          rdm1_flukt_trace = rdm1_flukt_trace + G1_hh_flukt(q,q)
                    end if
              end do

              E_hherpa_one = zero
              E_hherpa_one_flukt = zero
              E_hherpa_one_s = zero
              E_hherpa_one_t = zero
              E_hherpa_one_taa = zero
              E_hherpa_one_tbb = zero
              do p = 1, NBasis
                    do q = 1, NBasis
                          E_hherpa_one = E_hherpa_one + G1_hh(p, q) * AuxData%HNO0(p,q)
                          E_hherpa_one_s = E_hherpa_one_s + G1_hh_s(p, q) * AuxData%HNO0(p,q)
                          if (z==1)then
                                E_hherpa_one_flukt = E_hherpa_one_flukt + G1_hh_flukt(p, q) * AuxData%HNO0(p,q)
                          end if
                          ! if (abs(G1_hh_s(p, q) * AuxData%HNO0(p,q)).gt.1.d-5)then
                          !       write(*, '(A5, 2I5, 4F12.6)')'pla', p, q, G1_hh_s(p, q) , AuxData%HNO0(p,q), G1_hh_s(p, q) * AuxData%HNO0(p,q), E_hherpa_one_s
                          ! end if

                          E_hherpa_one_t = E_hherpa_one_t + G1_hh_t(p, q) * AuxData%HNO0(p,q)
                          E_hherpa_one_taa = E_hherpa_one_taa + G1_hh_taa(p, q) * AuxData%HNO0(p,q)
                          E_hherpa_one_tbb = E_hherpa_one_tbb + G1_hh_tbb(p, q) * AuxData%HNO0(p,q)
                    end do
              end do
              
              print*, 'E_pperpa one-electron energy = ', E_hherpa_one
              print*, 'E_pperpa sing-one-electron energy = ', E_hherpa_one_s
              if (z==1)then
                    print*, 'E_pperpa sing-one-electron energy flukt= ', E_hherpa_one_flukt
              end if
              print*, 'E_pperpa trip-one-electron energyab = ', E_hherpa_one_t
              print*, 'E_pperpa trip-one-electron energyaa = ', E_hherpa_one_taa
              print*, 'E_pperpa trip-one-electron energybb = ', E_hherpa_one_tbb
              print*, 'E_pperpa sing-all-electron energy = ', E_hherpa_one_s + E_hherpa_s
              print*, 'E_pperpa trip-all-electron energy = ', E_hherpa_one_t + E_hherpa_one_taa + E_hherpa_one_tbb+ E_hherpa_t+ E_hherpa_taa+ E_hherpa_tbb
              
              print*, 'total energy (hh)            ', AuxData%ENuc + E_hherpa_one + E_hherpa
              print*, 'nuclear energy               = ', AuxData%ENuc

              ! ================================================================
              ! Symmetry checks
              ! ================================================================
              max_diff_herm = zero
              max_diff_perm = zero

              do p = 1, NBasis
                    do q = 1, NBasis
                          do r = 1, NBasis
                                do s = 1, NBasis
                                      val_ref = G2_hh(p, q, r, s)
                                      diff = abs(val_ref - G2_hh(r, s, p, q))
                                      if (diff > max_diff_herm) then
                                            max_diff_herm = diff
                                      end if
                                      diff = abs(val_ref - G2_hh(q, p, s, r))
                                      if (diff > max_diff_perm) max_diff_perm = diff
                                end do
                          end do
                    end do
              end do



              ! Spin-summed total: trace N(N-1), partial trace daje gamma (slad N), n_max=2
            ! call check_2rdm_nreppp("hhERPA total (spin-summed)", G2_hh, NBasis, &
            !                      real(AuxData%Nel*(AuxData%Nel-1), F64), &
            !                      real(AuxData%Nel, F64) - 1.0d0, &
            !                      real(AuxData%Nel, F64), 2.0d0, .false.)

            ! Opcjonalnie poszczegolne bloki spinowe:
            !   G2_hh_taa (alpha-alpha): trace Na(Na-1), antysym p<->q tak
            !   G2_hh_s   (singlet ab):  trace Na*Nb,    antysym p<->q nie
            ! block
            !     double precision :: Na
            !     Na = real(AuxData%Nel, F64) * 0.5d0
            !     call check_2rdm_nreppp("hhERPA alpha-alpha", G2_hh_taa, NBasis, &
            !                          Na*(Na-1.0d0), Na - 1.0d0, Na, 1.0d0, .true.)
            !     call check_2rdm_nreppp("hhERPA singlet ab", G2_hh_s, NBasis, &
            !                          Na*Na, Na, Na, 1.0d0, .false.)
            ! end block

              ! ================================================================
              ! Save and transform
              ! ================================================================
              call save_2rdm_text('HHrdm2_natural.dat', G2_hh)

               if (allocated(AuxData%CMONO)) then
                     print*, 'Transforming G2_hh back to Original/MO basis'
                     allocate(work(NBasis, NBasis))
                     work = transpose(AuxData%CMONO)
                     call rdm2_MO_NO_trans(G2_hh, work, NBasis)
                     deallocate(work)
               end if

               call save_2rdm_text('HHrdm2_original.dat', G2_hh)
               call save_2rdm_text('HHERPA_2RDM.dat', G2_hh)

              print '(A)', "------------------------------------------------------------"
              print '(A, I15)',   " Number of electrons:                ", AuxData%Nel
              if (z==1)then
                    print '(A, F15.8)', " Reference 2-RDM trace [N*(N-1)]:    ", real(AuxData%Nel * (AuxData%Nel - 1), F64)
                    print '(A, F15.8)', " Calculated 2-RDM trace (hh):        ", rdm2_trace
                    print '(A, F15.8)', " Calculated 2-RDM trace flukt (hh):        ", rdm2_flukt_trace
              end if
              print '(A, F15.8)', " Calculated 1-RDM trace (hh):        ", rdm1_trace
              print '(A)',        " --- 2-RDM Symmetry Check ---"
              print '(A, E15.4)', " Max Error Hermiticity (pqrs|rspq):  ", max_diff_herm
              print '(A, E15.4)', " Max Error Permutation (pqrs|qpsr):  ", max_diff_perm
              print '(A)',        " ------------------------------------------------------------"
              print '(A, A)',     "2-RDM matrix saved to:               ", 'HHrdm2.dat'
              print '(A)',        " ------------------------------------------------------------"
              print "(A)",        " --- 2-RDM Storage Convention ---"
              print "(A)",        " Stored as: kjli--> rdm2(l,k,i,j) - original MOLMPS convention"
              print "(A)",        " ------------------------------------------------------------"

              deallocate(AuxData%vplus_s)
              deallocate(AuxData%vplus_t)
              deallocate(AuxData%vplus_t_aa)
              deallocate(AuxData%vplus_t_bb)

              deallocate(AuxData%Eigs_s)
              deallocate(AuxData%Eigs_t)
              deallocate(AuxData%Eigs_t_aa)
              deallocate(AuxData%Eigs_t_bb)

              deallocate(AuxData%Eigvec_s)
              deallocate(AuxData%Eigvec_t)
              deallocate(AuxData%Eigvec_t_aa)
              deallocate(AuxData%Eigvec_t_bb)


            end associate

      end subroutine dump_hhRDM




      ! ====================================================================
      ! check_2rdm_nrep — diagnostyka N-reprezentowalności dla 2-RDM
      !
      ! Sprawdza:
      !   1. Sled Sum_pq Gamma(p,q,p,q)        — czy = oczekiwanej
      !   2. Hermitowskosc      G(p,q,r,s) = G(r,s,p,q)
      !   3. Permutacja par     G(p,q,r,s) = G(q,p,s,r)
      !   4. Antysymetria p<->q (opcjonalne, tylko same-spin)
      !   5. Partial trace -> 1-RDM, slad i widmo (czy w [0, n_max])
      !   6. D-condition: wartosci wlasne Gamma jako macierzy par >= 0
      ! ====================================================================
      subroutine check_2rdm_nreppp(label, Gamma, NBasis, &
                                 tr_expected, denom, gamma_tr_expected, &
                                 n_max, check_antisym)
            character(*), intent(in) :: label
            integer, intent(in)      :: NBasis
            double precision, intent(in) :: Gamma(NBasis, NBasis, NBasis, NBasis)
            double precision, intent(in) :: tr_expected, denom, gamma_tr_expected, n_max
            logical, intent(in)      :: check_antisym

            integer :: p, q, r, s, info, dim2, lwork
            double precision :: tr_g, tr_g1
            double precision :: max_herm, max_perm, max_antisym
            double precision :: min_e1, max_e1, min_e2
            double precision, allocatable :: gamma1(:,:), eig(:), work(:), Gflat(:,:)

            print '(A)', "============================================================"
            print '(A,A)', " N-rep diagnostics: ", trim(label)
            print '(A)', "============================================================"

            ! ---- 1) Slad ----
            tr_g = 0.0d0
            do p = 1, NBasis
                  do q = 1, NBasis
                        tr_g = tr_g + Gamma(p,q,p,q)
                  end do
            end do
            print '(A,F16.8,A,F16.8)', " Tr Gamma            = ", tr_g, &
                                       "   expected: ", tr_expected
            print '(A,F16.8)',         " deviation           = ", tr_g - tr_expected

            ! ---- 2) Hermitowskosc ----
            max_herm = 0.0d0
            do p = 1, NBasis
              do q = 1, NBasis
                do r = 1, NBasis
                  do s = 1, NBasis
                    max_herm = max(max_herm, abs(Gamma(p,q,r,s) - Gamma(r,s,p,q)))
                  end do
                end do
              end do
            end do
            print '(A,E16.6)', " Max |G(pqrs)-G(rspq)|  = ", max_herm

            ! ---- 3) Permutacja par (p,q)<->(q,p) wraz z (r,s)<->(s,r) ----
            max_perm = 0.0d0
            do p = 1, NBasis
              do q = 1, NBasis
                do r = 1, NBasis
                  do s = 1, NBasis
                    max_perm = max(max_perm, abs(Gamma(p,q,r,s) - Gamma(q,p,s,r)))
                  end do
                end do
              end do
            end do
            print '(A,E16.6)', " Max |G(pqrs)-G(qpsr)|  = ", max_perm

            ! ---- 4) Antysymetria p<->q (tylko dla same-spin) ----
            if (check_antisym) then
              max_antisym = 0.0d0
              do p = 1, NBasis
                do q = 1, NBasis
                  do r = 1, NBasis
                    do s = 1, NBasis
                      max_antisym = max(max_antisym, abs(Gamma(p,q,r,s) + Gamma(q,p,r,s)))
                    end do
                  end do
                end do
              end do
              print '(A,E16.6)', " Max |G(pqrs)+G(qprs)|  = ", max_antisym
            else
              print '(A)',       " (antysymetria p<->q nie dotyczy mixed-spin)"
            end if

            ! ---- 5) Partial trace -> 1-RDM ----
            allocate(gamma1(NBasis, NBasis))
            gamma1 = 0.0d0
            do q = 1, NBasis
              do s = 1, NBasis
                do p = 1, NBasis
                  gamma1(q,s) = gamma1(q,s) + Gamma(p,q,p,s)
                end do
              end do
            end do
            gamma1 = gamma1 / denom

            tr_g1 = 0.0d0
            do p = 1, NBasis
                  tr_g1 = tr_g1 + gamma1(p,p)
            end do
            print '(A,F16.8,A,F16.8)', " Tr gamma (p-trace)  = ", tr_g1, &
                                       "   expected: ", gamma_tr_expected

            ! Symetryzacja i diagonalizacja
            do p = 1, NBasis
              do q = 1, p-1
                gamma1(p,q) = 0.5d0*(gamma1(p,q) + gamma1(q,p))
                gamma1(q,p) = gamma1(p,q)
              end do
            end do
            lwork = 3*NBasis
            allocate(eig(NBasis), work(lwork))
            call DSYEV('N','U', NBasis, gamma1, NBasis, eig, work, lwork, info)
            min_e1 = minval(eig); max_e1 = maxval(eig)
            print '(A,F16.8,A,F16.8,A,F6.2,A)', &
                  " gamma widmo: [", min_e1, ", ", max_e1, "]  (powinno byc w [0, ", n_max, "])"
            if (min_e1 < -1.0d-6) print '(A)', " *** ujemna obsadnosc!"
            if (max_e1 > n_max + 1.0d-6) print '(A)', " *** przekroczone n_max!"
            deallocate(eig, work)

            ! ---- 6) D-condition: Gamma jako macierz na przestrzeni par ----
            dim2 = NBasis*NBasis
            allocate(Gflat(dim2, dim2), eig(dim2), work(3*dim2))
            do p = 1, NBasis
              do q = 1, NBasis
                do r = 1, NBasis
                  do s = 1, NBasis
                    Gflat((p-1)*NBasis + q, (r-1)*NBasis + s) = Gamma(p,q,r,s)
                  end do
                end do
              end do
            end do
            ! symetryzacja
            do p = 1, dim2
              do q = 1, p-1
                Gflat(p,q) = 0.5d0*(Gflat(p,q) + Gflat(q,p))
                Gflat(q,p) = Gflat(p,q)
              end do
            end do
            call DSYEV('N','U', dim2, Gflat, dim2, eig, work, 3*dim2, info)
            min_e2 = minval(eig)
            print '(A,F16.8,A)', " D-cond min eig     = ", min_e2, "  (powinno byc >= 0)"
            if (min_e2 < -1.0d-6) then
                  print '(A)', " *** D-condition LAMANY ***"
            end if
            deallocate(Gflat, eig, work)
            deallocate(gamma1)

            print '(A)', "============================================================"
      end subroutine check_2rdm_nreppp

      function one_el_contr(AuxData, p, q, r, s)
            double precision :: one_el_contr
            type(TACppData), intent(in) :: AuxData
            integer, intent(in) :: p, q, r, s

            one_el_contr = zero
            if ((s==q).and.(p==r)) one_el_contr = one_el_contr - four
            if ((q==r).and.(p==s)) one_el_contr = one_el_contr + two

            one_el_contr = one_el_contr &
                  + delgam(AuxData%rdm1_p, AuxData, s, q, p, r) &
                  - delgam(AuxData%rdm1_p, AuxData, q, r, p, s) &
                  - delgam(AuxData%rdm1_p, AuxData, p, s, q, r) &                  
                  + delgam(AuxData%rdm1_p, AuxData, p, r, q, s) &
                  
                  + delgam(AuxData%rdm1_m, AuxData, s, q, p, r) &
                  - delgam(AuxData%rdm1_m, AuxData, q, r, p, s) &
                  - delgam(AuxData%rdm1_m, AuxData, p, s, q, r) &
                  + delgam(AuxData%rdm1_m, AuxData, p, r, q, s) &
                  + delgam(AuxData%rdm1_p, AuxData, s, q, p, r) &
                  + delgam(AuxData%rdm1_m, AuxData, p, r, q, s) &
                  + delgam(AuxData%rdm1_m, AuxData, s, q, p, r) &
                  + delgam(AuxData%rdm1_p, AuxData, p, r, q, s) 


      end function one_el_contr

      function delgam2(AuxData, s, q, p, r)
            double precision :: delgam2
            type(TACppData), intent(in) :: AuxData
            integer, intent(in) :: s, q, p, r

            delgam2 = zero

            if (p<=AuxData%NIA.and.r<=AuxData%NIA)then
                  if (s==q) delgam2 = delgam2 + erdm(AuxData%rdm1_p, p, r, AuxData%IndAux, AuxData%NI)
            end if

      end function delgam2

      function delgam(rdm1, AuxData, s, q, p, r)
            double precision :: delgam
            double precision, dimension(:, :), intent(in) ::rdm1
            type(TACppData), intent(in) :: AuxData
            integer, intent(in) :: s, q, p, r

            delgam = zero

            if (p<=AuxData%NIA.and.r<=AuxData%NIA)then
                  if (s==q) delgam = delgam + erdm(rdm1, p, r, AuxData%IndAux, AuxData%NI)
            end if

      end function delgam
      
      function erdm(rdm1, p, q, IAux, NI)
            double precision :: erdm
            double precision, dimension(:,:), intent(in) :: rdm1
            integer, dimension(:), intent(in) :: IAux
            integer, intent(in) :: NI
            integer, intent(in) :: p,q

            if(IAux(p)==1.and.IAux(q)==1)then
                  erdm = rdm1(p - NI, q - NI)
            else if (IAux(p)==0.and. IAux(q) == 0)then
                  if (p==q)then
                        erdm = one
                  else
                        erdm = zero
                  end if
            else
                  erdm = zero
            end if

      end function erdm
      
      function erdm_pp(rdm2, rdm1, p, q, r, s, IAux, NI)
            double precision :: erdm_pp
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            double precision, dimension(:,:), intent(in) :: rdm1
            integer, dimension(:), intent(in) :: IAux
            integer, intent(in) :: NI
            integer, intent(in) :: s, p, q, r

            if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
                  erdm_pp = rdm2(p - NI, q - NI, r - NI, s - NI)
            else
                  erdm_pp = zero
                  if (p==r.and.q==s)then
                        erdm_pp = erdm_pp + erdm(rdm1, p, r, IAux, NI) * erdm(rdm1,q,s, IAux, NI)
                  end if
                  if (p==s.and.q==r)then
                        erdm_pp = erdm_pp - erdm(rdm1, p, s, IAux, NI) * erdm(rdm1, q,r, IAux, NI)
                  end if
            end if

      end function erdm_pp

      function erdm_pm(rdm2, rdm1_p, rdm1_m, p, q, r, s, IAux, NI)
            double precision :: erdm_pm
            double precision, dimension(:,:,:,:), intent(in) :: rdm2
            double precision, dimension(:,:), intent(in) :: rdm1_p, rdm1_m
            integer, dimension(:), intent(in) :: IAux
            integer, intent(in) :: NI
            integer, intent(in) :: s, p, q, r

            if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
                  erdm_pm = rdm2(p - NI, q - NI, r - NI, s - NI)
            else
              	erdm_pm = zero
                  if (p==r.and.q==s)then
                        erdm_pm = erdm_pm + erdm(rdm1_p, p, r, IAux, NI) * erdm(rdm1_m, q,s, IAux, NI)
                  end if
            end if

      end function erdm_pm


      ! function erdm_ppx(rdm2, n, p, q, r, s, IAux, NI)
      !       double precision :: erdm_ppx
      !       double precision, dimension(:,:,:,:), intent(in) :: rdm2
      !       double precision, dimension(:), intent(in) :: n
      !       integer, dimension(:), intent(in) :: IAux
      !       integer, intent(in) :: NI
      !       integer, intent(in) :: s, p, q, r

      !       if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
      !             erdm_ppx = rdm2(p - NI, q - NI, r - NI, s - NI)
      !       else
      !             erdm_ppx = zero
      !             if (p==r.and.q==s)then
      !                   erdm_ppx = erdm_ppx + n(p) * n(q)
      !             end if
      !             if (p==s.and.q==r)then
      !                   erdm_ppx = erdm_ppx - n(p) * n(q) 
      !             end if
      !       end if

      ! end function erdm_ppx

      ! function erdm_pmx(rdm2, n_p, n_m, p, q, r, s, IAux, NI)
      !       double precision :: erdm_pmx
      !       double precision, dimension(:,:,:,:), intent(in) :: rdm2
      !       double precision, dimension(:), intent(in) :: n_p, n_m
      !       integer, dimension(:), intent(in) :: IAux
      !       integer, intent(in) :: NI
      !       integer, intent(in) :: s, p, q, r

      !       if(IAux(p)==1.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
      !             erdm_pmx = rdm2(p - NI, q - NI, r - NI, s - NI)
      !       else
      !         	erdm_pmx = zero
      !             if (p==r.and.q==s)then
      !                   erdm_pmx = erdm_pmx + n_p(p) * n_m(q)
      !             end if
      !       end if

      ! end function erdm_pmx


      subroutine PPERPA_init(AuxData, ACAlpha, Flags, TwoEl, TwoNO)
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            double precision, dimension(:),      intent(out) :: TwoNO
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags

            integer :: twoint_dim, r, l
            integer :: i, j, k, kl, ij, t
            double precision :: temp

            associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
                  NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
                  NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux)

              print*, 'ACAlpha', ACAlpha

              if (.not. allocated(AuxData%HNOA)) allocate(AuxData%HNOA(NBasis, NBasis))

              AuxData%HNOA = ACAlpha * AuxData%HNO0
              
              do i = 1, NBasis
                    do j = 1, NBasis
                          !print*, i, j
                          if (IAux(i)==IAux(j))then
                                AuxData%HNOA(i,j) = AuxData%HNOA(i,j) + (One-ACAlpha) * AuxData%HNO0(i,j)

                                if (IAux(i)==1)then
                                      l = NI
                                else
                                      l = NIA
                                end if
                                temp = zero
                                do t = 1, l
                                      temp = temp + occ(t)* (two*TwoEl(gmap(t,t,i,j))-TwoEl(gmap(t,i,t,j))) 
                                end do

                                AuxData%HNOA(i,j) = AuxData%HNOA(i,j) + (One-ACAlpha) * temp
                                ! if (abs(AuxData%HNOA(i,j)).gt.1.d-8)then
                                !       print*, 'hnohnoa', i, j, AuxData%HNOA(i,j)
                                ! end if
                          end if
                    end do
              end do

!              stop
            TwoNO = TwoEL
            ij = 0
            do i = 1, NBasis
                  do j = 1, i
                        ij = ij + 1
                        kl = 0
                        do k = 1, Nbasis
                              do l = 1, k
                                    kl=kl+1
                                    if ((IAux(i)==IAux(j)).and.(IAux(i)==IAux(k)).and.IAux(i)==IAux(l).and.IAux(i)==1)then                                          
                                          TwoNO(gmap(i, j, k, l)) = TwoNO(gmap(i, j, k, l))
                                    else
                                          TwoNO(gmap(i, j, k, l)) = ACAlpha * TwoNO(gmap(i, j, k, l))
                                    end if

                              end do
                        end do
                  end do
            end do
          end associate
      end subroutine PPERPA_init

           subroutine ddot_norm(vri, S, vrj, n, dd)
            double precision, dimension(:), intent(in) :: vri, vrj
	    double precision, dimension(:, :), intent(in) :: S
            integer, intent(in) :: n
            double precision, intent(out) :: dd
            double precision, dimension(:), allocatable :: tempx

            allocate(tempx(n))
            tempx = zero
            call real_av_x(tempx, S, n,  vrj, n, n, one, zero)

            call real_vw_x(dd, vri, tempx, n)
            deallocate(tempx)

      end subroutine ddot_norm

    subroutine orthogonalize_degen(n, countj, wr, vr, S_diag, dy, x)

            integer, intent(in) :: n
            double precision, dimension(:), intent(in) :: wr
            double precision, dimension(:,:), intent(inout) :: vr
            double precision, dimension(:), intent(in) :: S_diag
            integer, dimension(:), intent(in) :: dy
            integer, intent(in) :: x
            integer, dimension(:), allocatable :: StartIdx, EndIdx
            integer, intent(in) :: countj
            double precision, parameter :: tol = 1.d-4
            integer :: count, j, i


            allocate(StartIdx(n))
            allocate(EndIdx(n))

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
                        if (abs(wr(i)-wr(i-1)).lt.tol)then

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

            if (StartIdx(1) == 0)then
                  count = 0
            end if


            call Orthogonalize(vr, count, StartIdx(1:count), EndIdx(1:count), S_diag, dy, x)
      end subroutine orthogonalize_degen

      function func_P3a(AuxData, TwoEl, s, p, q, r)
            double precision :: func_P3a
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, p, q, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)


              func_P3a = zero

              ! if (s==5.and.p==4.and.q==5.and.r==18)then
              !       write(*, '(A7, F12.8)') 'aixi_1', func_P3a
              ! end if


              if (IAux(s)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_P3a = func_P3a + (erdm_pmx(rdm2_pm, n_p, n_m, t, s, p, u, IAux, NI) &
                                      + erdm_pmx(rdm2_mp, n_m, n_p, t, s, p, u, IAux, NI))* twoel(gmap(q,u,r,t))
                          end do
                    end do

                    ! if (s==11.and.p==12.and.q==11.and.r==12)then
                    !       write(*, '(A7, F12.8)') 'aixiA', func_P3a
                    ! end if

                    ! if (s==12.and.p==11.and.q==12.and.r==11)then
                    !       write(*, '(A7, F12.8)') 'aixiB', func_P3a
                    ! end if

                    ! if (s==11.and.p==11.and.q==12.and.r==12)then
                    !       write(*, '(A7, F12.8)') 'aixiC', func_P3a
                    ! end if

                    ! if (s==12.and.p==12.and.q==11.and.r==11)then
                    !       write(*, '(A7, F12.8)') 'aixiD', func_P3a
                    ! end if

              else
                    func_P3a = (n_p(p) * n_m(s) + n_m(p) * n_p(s)) * TwoEl(gmap(q, s, r, p))

                    
                    ! if (s==11.and.p==12.and.q==11.and.r==12)then
                    !       write(*, '(A7, F12.8)') 'aixiAx', func_P3a
                    ! end if

                    ! if (s==12.and.p==11.and.q==12.and.r==11)then
                    !       write(*, '(A7, F12.8)') 'aixiBx', func_P3a
                    ! end if

                    ! if (s==11.and.p==11.and.q==12.and.r==12)then
                    !       write(*, '(A7, F12.8)') 'aixiCx', func_P3a
                    ! end if

                    ! if (s==12.and.p==12.and.q==11.and.r==11)then
                    !       write(*, '(A7, F12.8)') 'aixiDx', func_P3a
                    ! end if

                    ! if (s==11.and.p==12.and.q==11.and.r==12)then
                    !       !if (s==5.and.p==4.and.q==5.and.r==18)then
                    !       write(*, '(A7, F12.8)') 'aixi_3', func_P3a
                    ! end if

              end if
            end associate
      end function func_P3a


      function func_P3b(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3b
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3b = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3b = func_P3b + (rdm2_pm(t - NI, s - NI, u - NI, q - NI) + &
                                    rdm2_mp(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, u, r, t))
                        end do
                  end do

                    if (s==11.and.p==11.and.q==12.and.r==12)then
                          write(*, '(A7, F12.8)') 'bixiA', func_P3b
                    end if

                    if (s==12.and.p==12.and.q==11.and.r==11)then
                          write(*, '(A7, F12.8)') 'bixiB', func_P3b
                    end if

                    if (s==11.and.p==12.and.q==11.and.r==12)then
                          write(*, '(A7, F12.8)') 'bixiC', func_P3b
                    end if

                    if (s==12.and.p==11.and.q==12.and.r==11)then
                          write(*, '(A7, F12.8)') 'bixiD', func_P3b
                    end if


                  
            end if

            if (s == q) then
                  do t = 1, NIA
                        if (IAux(s) == 0 .or. IAux(t) == 0) then
                              func_P3b = func_P3b + (n_p(t) * n_m(s) + n_m(t) * n_p(s)) * TwoEl(gmap(p, t, r, t))
                        end if
                  end do
            end if
            end associate
      end function func_P3b

      function func_P3a_aa(AuxData, TwoEl, s, p, q, r)
            double precision :: func_P3a_aa
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, p, q, r
            integer :: t, u

            associate(rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)


              func_P3a_aa = zero

              if (IAux(s)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_P3a_aa = func_P3a_aa + rdm2_pp(t-NI, s-NI, p-NI, u-NI)* twoel(gmap(q,u,r,t))
                          end do
                    end do

                    if (s==p)then
                          do t = 1, NI
                                func_P3a_aa = func_P3a_aa  - n_p(t) * n_p(s) * TwoEl(gmap(q, t, r, t))
                          end do
                    end if
              else
                    if (s==p)then
                          do t = 1, NIA
                                func_P3a_aa = func_P3a_aa  - n_p(t) * n_p(s) * TwoEl(gmap(q, t, r, t))
                          end do
                    end if
                    
                    func_P3a_aa = func_P3a_aa + n_p(p) * n_p(s)  * TwoEl(gmap(q, s, r, p))

              end if
            end associate
      end function func_P3a_aa

      
      function func_P3a_bb(AuxData, TwoEl, s, p, q, r)
            double precision :: func_P3a_bb
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, p, q, r
            integer :: t, u

            associate(rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_m=>AuxData%n_m)


              func_P3a_bb = zero

              if (IAux(s)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_P3a_bb = func_P3a_bb + rdm2_mm(t-NI, s-NI, p-NI, u-NI)* twoel(gmap(q,u,r,t))
                          end do
                    end do

                    if (s==p)then
                          do t = 1, NI
                                func_P3a_bb = func_P3a_bb  - n_m(t) * n_m(s) * TwoEl(gmap(q, t, r, t))
                          end do
                    end if
              else
                    if (s==p)then
                          do t = 1, NIA
                                func_P3a_bb = func_P3a_bb  - n_m(t) * n_m(s) * TwoEl(gmap(q, t, r, t))
                          end do
                    end if
                    
                    func_P3a_bb = func_P3a_bb + n_m(p) * n_m(s)  * TwoEl(gmap(q, s, r, p))

              end if
            end associate
      end function func_P3a_bb

      
      function func_P3b_aa(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3b_aa
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3b_aa = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3b_aa = func_P3b_aa + rdm2_pp(t - NI, s - NI, u - NI, q - NI) * TwoEl(gmap(p, u, r, t))
                        end do
                  end do

                  if (s==q)then
                        do t = 1, NI
                              func_P3b_aa = func_P3b_aa + n_p(t) * n_p(s) * TwoEl(gmap(p, t, r, t))
                        end do
                  end if
            else

                  if (s == q) then
                        do t = 1 , NIA
                              func_P3b_aa = func_P3b_aa + n_p(t) * n_p(s) * TwoEl(gmap(p, t, r, t))
                        end do
                  end if
                  func_P3b_aa = func_P3b_aa - n_p(q) * n_p(s)  * TwoEl(gmap(p, s, r, q))
            end if
            end associate
      end function func_P3b_aa

      function func_P3b_bb(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3b_bb
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3b_bb = zero

              ! if (r==18.and.s==5.and.p==5.and.q==4)then
              !       write(*, '(A7, F12.8)') 'xixi_1', func_P3b_bb
              ! end if
              
              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3b_bb = func_P3b_bb + rdm2_mm(t - NI, s - NI, u - NI, q - NI) * TwoEl(gmap(p, u, r, t))
                        end do
                  end do

                  ! if (r==18.and.s==5.and.p==5.and.q==4)then
                  !       write(*, '(A7, F12.8)') 'xixi_2', func_P3b_bb
                  ! end if


                  if (s==q)then
                        do t = 1, NI
                              func_P3b_bb = func_P3b_bb + n_m(t) * n_m(s) * TwoEl(gmap(p, t, r, t))
                        end do
                  end if
                  ! if (r==18.and.s==5.and.p==5.and.q==4)then
                  !       write(*, '(A7, F12.8)') 'xixi_3', func_P3b_bb
                  ! end if

            else

                  if (s == q) then
                        do t = 1 , NIA
                              func_P3b_bb = func_P3b_bb + n_m(t) * n_m(s) * TwoEl(gmap(p, t, r, t))
                        end do
                  end if

                  ! if (r==18.and.s==5.and.p==5.and.q==4)then
                  !       write(*, '(A7, F12.8)') 'xixi_4', func_P3b_bb
                  ! end if

                  func_P3b_bb = func_P3b_bb - n_m(q) * n_m(s)  * TwoEl(gmap(p, s, r, q))

                  ! if (r==18.and.s==5.and.p==5.and.q==4)then
                  !       write(*, '(A7, F12.8)') 'xixi_5', func_P3b_bb
                  ! end if

            end if
            end associate
      end function func_P3b_bb



      function func_P3c_aa(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3c_aa
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, &
                  rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3c_aa = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3c_aa = func_P3c_aa + (rdm2_mp(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_pp(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, r, t, u))
                        end do
                  end do
            end if

            if (s == q) then
                  do t = 1, NIA
                        if (IAux(s) == 0 .or. IAux(t) == 0) then
                              func_P3c_aa = func_P3c_aa + (n_m(t) * n_p(s) + &
                                          n_p(t) * n_p(s)) * TwoEl(gmap(p, r, t, t))
                        end if
                  end do
            end if

            if (IAux(s) == 0 .or. IAux(q) == 0) then
                  func_P3c_aa = func_P3c_aa -  n_p(q) * n_p(s) * TwoEl(gmap(p, r, q, s))
            end if
            end associate
      end function func_P3c_aa

      
      function func_P3c(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3c
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, &
                  rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3c = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3c = func_P3c + (rdm2_pm(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_mp(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_mm(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_pp(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, r, t, u))
                        end do
                  end do
            end if

            if (s == q) then
                  do t = 1, NIA
                        if (IAux(s) == 0 .or. IAux(t) == 0) then
                              func_P3c = func_P3c + (n_p(t) * n_m(s) + n_m(t) * n_p(s) + &
                                         n_m(t) * n_m(s) + n_p(t) * n_p(s)) * TwoEl(gmap(p, r, t, t))
                        end if
                  end do
            end if

            if (IAux(s) == 0 .or. IAux(q) == 0) then
                  func_P3c = func_P3c - (n_m(q) * n_m(s) + n_p(q) * n_p(s)) * TwoEl(gmap(p, r, q, s))
            end if
            end associate
      end function func_P3c

            function func_P3c_bb(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3c_bb
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, &
                  NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_m=>AuxData%n_m, n_p=>AuxData%n_p)

              func_P3c_bb = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3c_bb = func_P3c_bb + (rdm2_pm(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_mm(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, r, t, u))
                        end do
                  end do
            end if

            if (s == q) then
                  do t = 1, NIA
                        if (IAux(s) == 0 .or. IAux(t) == 0) then
                              func_P3c_bb = func_P3c_bb + (n_p(t) * n_m(s) + &
                                          n_m(t) * n_m(s)) * TwoEl(gmap(p, r, t, t))
                        end if
                  end do
            end if

            if (IAux(s) == 0 .or. IAux(q) == 0) then
                  func_P3c_bb = func_P3c_bb -  n_m(q) * n_m(s) * TwoEl(gmap(p, r, q, s))
            end if
            end associate
      end function func_P3c_bb

      function func_P45_bb(AuxData, TwoEl, q, s)
            double precision :: func_P45_bb
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: q, s
            integer :: t, u, v

            associate(rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, &
                  NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_m=>AuxData%n_m, n_p=>AuxData%n_p)

              func_P45_bb = zero
              ! if (q==18.and.s==4)then

              !       write(*, '(A7,F12.8)') 'zizi_1', func_P45_bb
              ! end if

              if (IAux(q) == 1) then
                    do t = NI + 1, NIA
                          do u = NI + 1, NIA
                                do v = NI + 1, NIA
                                      func_P45_bb = func_P45_bb + frac12 * ( &
                                            rdm2_mm(t - NI, u - NI, q - NI, v - NI) + &
                                            rdm2_mp(t - NI, u - NI, q - NI, v - NI) ) * TwoEl(gmap(s, t, v, u))
                                end do
                          end do
                    end do
                    ! if (q==18.and.s==4)then
                    !       write(*, '(A7,F12.8)') 'zizi_2', func_P45_bb
                    ! end if

              end if

              do v = 1, NIA
                    if (IAux(q) == 0 .or. IAux(v) == 0) then

                          func_P45_bb = func_P45_bb + frac12 * (n_m(q) * n_m(v) +  n_m(q) * n_p(v)) * TwoEl(gmap(s, q, v, v))
                          ! if (q==18.and.s==4)then
                          !       print*, n_m(q), n_p(v), n_m(q)* n_p(v)
                          !       write(*, '(A7, 2I3, 6F12.6)')'uouo', q, v, n_m(q), n_p(v), n_m(q) * n_m(v) ,  n_m(q) + n_p(v) , TwoEl(gmap(s, q, v, v)), (n_m(q) * n_m(v) +  n_m(q) + n_p(v)) * TwoEl(gmap(s, q, v, v))
                          ! end if
                    end if
              end do
              ! if (q==18.and.s==4)then
              !       write(*, '(A7,F12.8)') 'zizi_3', func_P45_bb
              ! end if
              do v = 1, NIA
                    if (IAux(q) == 0 .or. IAux(v) == 0) then
                          func_P45_bb = func_P45_bb - frac12 * (n_m(v) * n_m(q)) * TwoEl(gmap(s, v, v, q))
                    end if
              end do
              ! if (q==18.and.s==4)then
              !       write(*, '(A7,F12.8)') 'zizi_4', func_P45_bb
              ! end if


            end associate
      end function func_P45_bb


      function func_P45_aa(AuxData, TwoEl, q, s)
            double precision :: func_P45_aa
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: q, s
            integer :: t, u, v

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, &
                  rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P45_aa = zero


              if (IAux(q) == 1) then
                    do t = NI + 1, NIA
                          do u = NI + 1, NIA
                                do v = NI + 1, NIA
                                      func_P45_aa = func_P45_aa + frac12 * ( &
                                           rdm2_pp(t - NI, u - NI, q - NI, v - NI) + &
                                           rdm2_pm(t - NI, u - NI, q - NI, v - NI) ) * TwoEl(gmap(s, t, v, u))
                                end do
                          end do
                    end do
              end if

              do v = 1, NIA
                    if (IAux(q) == 0 .or. IAux(v) == 0) then

                          func_P45_aa = func_P45_aa + frac12 * (n_p(q) * n_p(v) +  n_p(q) * n_m(v)) * TwoEl(gmap(s, q, v, v))
                          func_P45_aa = func_P45_aa - frac12 * (n_p(v) * n_p(q)) * TwoEl(gmap(s, v, v, q))
                    end if
              end do
            end associate
      end function func_P45_aa


      
      function func_P45(AuxData, TwoEl, q, s)
            double precision :: func_P45
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: q, s
            integer :: t, u, v

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, &
                  rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P45 = zero


            if (IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              do v = NI + 1, NIA
                                    func_P45 = func_P45 + frac12 * ( &
                                          rdm2_mm(t - NI, u - NI, q - NI, v - NI) + &
                                          rdm2_pp(t - NI, u - NI, q - NI, v - NI) + &
                                          rdm2_pm(t - NI, u - NI, q - NI, v - NI) + &
                                          rdm2_mp(t - NI, u - NI, q - NI, v - NI) ) * TwoEl(gmap(s, t, v, u))
                              end do
                        end do
                  end do


            end if

                  do v = 1, NIA
                        if (IAux(q) == 0 .or. IAux(v) == 0) then
                              func_P45 = func_P45 + frac12 * (n_p(q) * n_p(v) + n_p(q) * n_m(v) + n_m(q) * n_p(v)+ n_m(q) * n_m(v)) * TwoEl(gmap(s, q, v, v))
                        end if
                  end do

                  do v = 1, NIA
                        if (IAux(q) == 0 .or. IAux(v) == 0) then
                              func_P45 = func_P45 - frac12 * (n_m(v) * n_m(q) + n_p(v) * n_p(q)) * TwoEl(gmap(s, v, v, q))
                        end if
                  end do
                  ! if (q==18.and.s==4)then
                  !       write(*, '(A7,F12.8)') 'yiyi_4', func_P45
                  ! end if

            end associate
      end function func_P45

      
      subroutine save_2rdm_text(rdm_file, rdm)
            character(len=*), intent(in) :: rdm_file
            real(F64), dimension(:,:,:,:), intent(in) :: rdm
            integer :: i, j, k, l, n

            n = size(rdm, dim=1)
            open(unit=20, file=rdm_file, status='replace', action='write')

!            write(20, '(A)') "2-RDM DATA"
!            write(20, '(A, I5)') "DIMENSION: ", n

            do k = 1, n
                  do j = 1, n
                        do l = 1, n
                              do i = 1, n
                                    if (abs(rdm(l, k, i, j)).gt.1.d-10)then
                                          write(20, '(4I5, E24.16)') l, k, i, j, rdm(l, k, i, j)
                                    end if
                              end do
                        end do
                  end do
            end do

            close(20)
      end subroutine save_2rdm_text

        subroutine rdm2_MO_NO_trans(rdm2, CMONO, NA)

          double precision, dimension(:,:), intent(in) :: CMONO
          integer, intent(in) :: NA
          double precision, dimension(:,:, :,:), intent(inout) :: rdm2
          double precision, dimension(:,:), allocatable :: work1, work2
          integer :: i, j


          allocate(work1(NA,NA))
          allocate(work2(NA, NA))

          do i = 1, NA
                do j = 1, NA
                      work1 =  rdm2(i, j, :,:)
                      call real_ab(work2, work1, CMONO)
                      work1 = zero
                      call real_aTb(work1, CMONO, work2)
                      rdm2(i, j, :,:) = work1
                end do
          end do

          do i = 1, NA
                do j = 1, NA
                      call real_ab(work2, rdm2(:,:,i,j), CMONO)
                      work1 = zero
                      call real_aTb(rdm2(:,:, i, j), CMONO, work2)
                end do
          end do


    end subroutine rdm2_MO_NO_trans

    subroutine rdm1_MO_NO_trans(rdm1, CMONO, NA)
          implicit none
          integer, intent(in) :: NA
          double precision, dimension(:,:), intent(in)    :: CMONO
          double precision, dimension(:,:), intent(inout) :: rdm1
          double precision, allocatable :: work1(:,:), work2(:,:)

          allocate(work1(NA,NA))
          allocate(work2(NA,NA))

          call real_ab(work2, rdm1, CMONO)

          work1 = 0.0d0
          call real_aTb(work1, CMONO, work2)

          rdm1 = work1

          deallocate(work1)
          deallocate(work2)
    end subroutine rdm1_MO_NO_trans


      function func_P3a_pm(AuxData, TwoEl, s, p, q, r)
            double precision :: func_P3a_pm
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, p, q, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3a_pm = zero

              if (IAux(s)==1 .and. IAux(p)==1)then
                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_P3a_pm = func_P3a_pm + rdm2_pm(t - NI, s - NI, p - NI, u - NI) * twoel(gmap(q,u,r,t))
                                !func_P3a_pm = func_P3a_pm + (erdm_pmx(rdm2_pm, n_p, n_m, t, s, p, u, IAux, NI)) * twoel(gmap(q,u,r,t))
                          end do
                    end do
              else
                    func_P3a_pm = (n_p(p) * n_m(s)) * TwoEl(gmap(q, s, r, p))
              end if
            end associate
      end function func_P3a_pm

      function func_P3a_mp(AuxData, TwoEl, s, p, q, r)
            double precision :: func_P3a_mp
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, p, q, r
            integer :: t, u

            associate(rdm2_mp=>AuxData%rdm2_mp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3a_mp = zero

              if (IAux(s)==1 .and. IAux(p)==1)then
                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_P3a_mp = func_P3a_mp + rdm2_mp(t - NI, s - NI, p - NI, u - NI) * twoel(gmap(q,u,r,t))
                                !func_P3a_mp = func_P3a_mp + (erdm_pmx(rdm2_mp, n_m, n_p, t, s, p, u, IAux, NI)) * twoel(gmap(q,u,r,t))
                          end do
                    end do
              else
                    func_P3a_mp = (n_m(p) * n_p(s)) * TwoEl(gmap(q, s, r, p))
              end if
            end associate
      end function func_P3a_mp

      function func_P3b_pm(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3b_pm
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3b_pm = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3b_pm = func_P3b_pm + rdm2_pm(t - NI, s - NI, u - NI, q - NI) * TwoEl(gmap(p, u, r, t))
                        end do
                  end do
              end if

              if (s == q) then
                    do t = 1, NI
                          func_P3b_pm = func_P3b_pm + (n_p(t) * n_m(s)) * TwoEl(gmap(p, t, r, t))
                    end do
                    do t = NI+1, NIA
                          if (IAux(s) == 0) then
                                func_P3b_pm = func_P3b_pm + (n_p(t) * n_m(s)) * TwoEl(gmap(p, t, r, t))
                          end if
                    end do
              end if
            end associate
      end function func_P3b_pm

      function func_P3b_mp(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3b_mp
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_mp=>AuxData%rdm2_mp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3b_mp = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3b_mp = func_P3b_mp + rdm2_mp(t - NI, s - NI, u - NI, q - NI) * TwoEl(gmap(p, u, r, t))
                        end do
                  end do
              end if

              if (s == q) then
                    do t = 1, NI
                          func_P3b_mp = func_P3b_mp + (n_m(t) * n_p(s)) * TwoEl(gmap(p, t, r, t))
                    end do
                    do t = NI+1, NIA
                          if (IAux(s) == 0) then
                                func_P3b_mp = func_P3b_mp + (n_m(t) * n_p(s)) * TwoEl(gmap(p, t, r, t))
                          end if
                    end do
              end if
            end associate
      end function func_P3b_mp


      function func_P3c_pp_mp(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3c_pp_mp
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_mp=>AuxData%rdm2_mp, rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3c_pp_mp = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3c_pp_mp = func_P3c_pp_mp + (rdm2_pp(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_mp(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, r, t, u))
                        end do
                  end do
              end if

              if (s == q) then
                    do t = 1, NI
                          func_P3c_pp_mp = func_P3c_pp_mp + (n_p(t) * n_p(q) + n_m(t) * n_p(q)) * TwoEl(gmap(p, r, t, t))
                    end do
                    do t = NI+1, NIA
                        if (IAux(s) == 0 ) then
                              func_P3c_pp_mp = func_P3c_pp_mp + (n_p(t) * n_p(q) + n_m(t) * n_p(q)) * TwoEl(gmap(p, r, t, t))
                        end if
                  end do
              end if

              if (IAux(s) == 0 .or. IAux(q) == 0) then
                  func_P3c_pp_mp = func_P3c_pp_mp - (n_p(q) * n_p(s)) * TwoEl(gmap(p, r, q, s))
              end if
            end associate
      end function func_P3c_pp_mp

      function func_P3c_mm_pm(AuxData, TwoEl, s, q, p, r)
            double precision :: func_P3c_mm_pm
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: s, q, p, r
            integer :: t, u

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P3c_mm_pm = zero

              if (IAux(s) == 1 .and. IAux(q) == 1) then
                  do t = NI + 1, NIA
                        do u = NI + 1, NIA
                              func_P3c_mm_pm = func_P3c_mm_pm + (rdm2_mm(t - NI, s - NI, u - NI, q - NI) + &
                                         rdm2_pm(t - NI, s - NI, u - NI, q - NI)) * TwoEl(gmap(p, r, t, u))
                        end do
                  end do
              end if

              if (s == q) then
                    do t = 1, NI
                          func_P3c_mm_pm = func_P3c_mm_pm + (n_m(t) * n_m(q) + n_p(t) * n_m(q)) * TwoEl(gmap(p, r, t, t))
                    end do
                    
                    do t = NI+1, NIA
                        if (IAux(s) == 0 .or. IAux(q) == 0) then
                              func_P3c_mm_pm = func_P3c_mm_pm + (n_m(t) * n_m(q) + n_p(t) * n_m(q)) * TwoEl(gmap(p, r, t, t))
                        end if
                  end do
            end if

            if (IAux(s) == 0 .or. IAux(q) == 0) then
                  func_P3c_mm_pm = func_P3c_mm_pm - (n_m(q) * n_m(s)) * TwoEl(gmap(p, r, q, s))
            end if
          end associate
      end function func_P3c_mm_pm


      function func_P45_pm(AuxData, TwoEl, q, s)
            double precision :: func_P45_pm
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: q, s
            integer :: t, u, v

            associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_pp=>AuxData%rdm2_pp, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P45_pm = zero

              if (IAux(q) == 1) then
                    do t = NI + 1, NIA
                          do u = NI + 1, NIA
                                do v = NI + 1, NIA
                                      func_P45_pm = func_P45_pm + frac12 * ( &
                                            rdm2_pp(t - NI, u - NI, q - NI, v - NI) + &
                                            rdm2_pm(t - NI, u - NI, q - NI, v - NI) ) * TwoEl(gmap(s, t, v, u))
                                end do
                          end do
                    end do
              end if

              do v = 1, NI
                      func_P45_pm = func_P45_pm + frac12 * (n_p(q) * n_p(v) + n_p(q) * n_m(v)) * TwoEl(gmap(s, q, v, v))
                      func_P45_pm = func_P45_pm - frac12 * (n_p(v) * n_p(q)) * TwoEl(gmap(s, v, v, q))
                end do
              do v = NI+1, NIA
                    if (IAux(q) == 0 ) then
                          func_P45_pm = func_P45_pm + frac12 * (n_p(q) * n_p(v) + n_p(q) * n_m(v)) * TwoEl(gmap(s, q, v, v))
                          func_P45_pm = func_P45_pm - frac12 * (n_p(v) * n_p(q)) * TwoEl(gmap(s, v, v, q))
                    end if
              end do
            end associate
      end function func_P45_pm

      function func_P45_mp(AuxData, TwoEl, q, s)
            double precision :: func_P45_mp
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:),      intent(in) :: TwoEl
            integer, intent(in) :: q, s
            integer :: t, u, v

            associate(rdm2_mp=>AuxData%rdm2_mp, rdm2_mm=>AuxData%rdm2_mm, NIA=>AuxData%NIA, NI=>AuxData%NI, &
                  IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m)

              func_P45_mp = zero

              if (IAux(q) == 1) then
                    do t = NI + 1, NIA
                          do u = NI + 1, NIA
                                do v = NI + 1, NIA
                                      func_P45_mp = func_P45_mp + frac12 * ( &
                                            rdm2_mm(t - NI, u - NI, q - NI, v - NI) + &
                                            rdm2_mp(t - NI, u - NI, q - NI, v - NI) ) * TwoEl(gmap(s, t, v, u))
                                end do
                          end do
                    end do
              end if

              do v = 1, NI
                    func_P45_mp = func_P45_mp + frac12 * (n_m(q) * n_m(v) + n_m(q) * n_p(v)) * TwoEl(gmap(s, q, v, v))
                    func_P45_mp = func_P45_mp - frac12 * (n_m(v) * n_m(q)) * TwoEl(gmap(s, v, v, q))
              end do
              do v = NI+1, NIA
                    if (IAux(q) == 0 ) then
                          func_P45_mp = func_P45_mp + frac12 * (n_m(q) * n_m(v) + n_m(q) * n_p(v)) * TwoEl(gmap(s, q, v, v))
                          func_P45_mp = func_P45_mp - frac12 * (n_m(v) * n_m(q)) * TwoEl(gmap(s, v, v, q))
                    end if
              end do
            end associate
      end function func_P45_mp

end module ppAC_incore
