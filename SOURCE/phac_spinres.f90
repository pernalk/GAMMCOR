module phac_spinres

      use ppac_types
      use types
      use math_constants
      use real_linalg 
      use THC_Gammcor
      use sort

      ! =============================================================================
      ! MODULE-LEVEL CONVENTION NOTE (read before touching any integral or RDM)
      ! =============================================================================
      !
      ! INTEGRAL STORAGE CONVENTION (FCIDUMP-derived):
      !   ints2e_aa(p,q,r,s)  stores 2*(pq|rs)_physical for αα block
      !   ints2e_bb(p,q,r,s)  stores 2*(pq|rs)_physical for ββ block
      !   ints2e_ab(p,q,r,s)  stores   (pq|rs)_physical for αβ block
      !   where (pq|rs) = ∫ φ_p*(1) φ_q(1) r12^-1 φ_r*(2) φ_s(2) dx1 dx2  [Mulliken]
      !
      !   CONSEQUENCE: any term contracting with ints2e_aa needs frac12 to cancel
      !   the built-in factor of 2. Terms with ints2e_ab need no extra factor.
      !   Example from energy: e2 = 0.25*(e2_aa + e2_bb) + e2_ab
      !   The 0.25 = 0.5 (antisymm) * 0.5 (FCIDUMP factor)
      !
      ! TRANSLATION TO DIRAC NOTATION:
      !   (pq|rs) [Mulliken] = <pr|qs> [Dirac/Pernal paper notation]
      !   So ints2e_aa(p,q,r,s) = 2*<pr|qs>_aa
      !
      ! FULL/FINAL INTEGRALS vs DUCC INTEGRALS:
      !   Full integrals (Ints):
      !     - 8-fold symmetry: (pq|rs)=(qp|sr)=(rs|pq)=(sr|qp)=(pq|sr)=(qp|rs)=(rs|pq)=(sr|qp)
      !     - ints2e_aa = ints2e_bb (spin-free for closed-shell reference)
      !     - stored via gmap_4fold but physically 8-fold symmetric
      !   DUCC integrals (IntsD, Ints_ducc):
      !     - 4-fold symmetry only: (pq|rs)=(qp|sr)=(rs|pq)=(sr|qp)
      !     - (pq|rs) ≠ (pq|sr) in general
      !     - ints2e_aa ≠ ints2e_ab in general (no spin-free simplification)
      !     - defined only in active space (indices 1..NA, global NI+1..NIA)
      !
      ! AC PATH HAMILTONIAN (Eq. 22-24 in Pernal JCP 2018):
      !   H^α = H^(0) + α*H'
      !   Two-electron part for AAAA block: g^α = (1-α)*g^DUCC + α*g^full  [erpa_init_incore_ducc]
      !   Two-electron part for VVVV block: g^α = g^full  (independent of α, already correlated)
      !   Two-electron part for mixed blocks: g^α = α*g^full
      !   One-electron part for AA block:  h^α = (1-α)*h^DUCC + α*h^full
      !   One-electron part for VV block:  h^α = h^full + (1-α)*[contraction with active RDM]
      !
      ! SPIN-FREE ERPA FOR SINGLET (Eq. 13-14 in ERPA_Singletxa notes):
      !   X_{p+q+} = X_{p-q-} ≡ X_pq  (shared spatial amplitude)
      !   Y_{p+q+} = Y_{p-q-} ≡ Y_pq
      !   The ERPA matrix built here is A^{0,0} = A^{++++} + A^{++--} (Eq. 54)
      !   which is spin-free and written in terms of spatial orbitals and
      !   spin-summed 2-RDM: Γ^{0,0}_{pqrs} = Γ^{++++}_{pqrs} + Γ^{+-+-}_{pqrs}
      !
      ! NORMALIZATION OF ERPA EIGENVECTORS (Eq. 221-225 in pherpa.pdf):
      !   Original: Y_ν^T N Y_ν - X_ν^T N X_ν = 1/2,  N_{pq,pq} = n_p - n_q
      !   In tilde basis (pherpa_symm_spinres): 2/ω * Ỹ^T A- Ỹ = 1
      !   Tilde vectors: X̃_pq = (c_p+c_q)(Y+X)_pq,  Ỹ_pq = (c_p-c_q)(Y-X)_pq
      !   where c_p = +√n_p if n_p≥0.5, c_p = -√n_p if n_p<0.5  [ssqrt convention]
      !   Key identity: (c_p+c_q)(c_p-c_q) = c_p²-c_q² = n_p-n_q
      !
      ! 2-RDM STORAGE (rdm2_pp, rdm2_pm etc., active indices 1..NA):
      !   rdm2_pp(p,r,q,s) = Γ^{αααα}_{pqrs} = <0|a†_p a†_q a_s a_r|0>  [AAAA block]
      !   rdm2_pm(p,r,q,s) = Γ^{αβαβ}_{pqrs}                             [AABB block]
      !   Antisymmetry: Γ_{pqrs} = -Γ_{qprs} = -Γ_{pqsr} = Γ_{qpsr}
      !   rdm1_p(p,q) = γ^α_{pq} = n_p δ_{pq}  (diagonal in natural orbital basis)
      ! =============================================================================

      implicit none
contains

  subroutine acph_driver(THCData, AuxData, Flags, Ints, IntsD)
        

        type(TACppData), intent(inout) :: AuxData
        type(FlagsData), intent(in) :: Flags
        type(TTHCData), intent(inout) :: THCData
        type(TInts), intent(inout) :: Ints
        type(TInts), intent(inout),optional :: IntsD
        integer :: ind, indt, i, j, ij, k
        integer :: spin_symm
        integer :: nn
        logical :: rdm_dump

        associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, NInte1=> AuxData%NInte1, &
              NInte2=>AuxData%NInte2, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
              NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux)

          write(LOUT,'(2x,a,4x,2e15.5)') 'Threshold for quasi-degeneracy ', AuxData%ThrSelAct
          write(LOUT,'(2x,a,4x,2e15.5)') 'Threshold for quasi-virtual orbital ', AuxData%ThrQVirt
          write(LOUT,'(2x,a,4x,2e15.5)') 'Threshold for quasi-inactive orbital ', AuxData%ThrQInact

          allocate(AuxData%map(NBasis))
          allocate(AuxData%IndN(2,Nbasis**2))
          allocate(AuxData%IndX(Nbasis**2))
          allocate(AuxData%IPair(Nbasis, NBasis))

          call imsg("NBASIS:", AuxData%NBasis)
          call imsg("NInactive:", AuxData%NI)
          call imsg("NActive:", AuxData%NA)
          call imsg("NVirt:", AuxData%NV)


          ij=0
          ind=0
          do i = 1, NBasis
                do j = 1, i-1
                      ij = ij + 1

                      if(IAux(i)+IAux(j).ne.0.and.IAux(i)+IAux(j).ne.4) then
                            if((IAux(i).eq.1).and.(IAux(j).eq.1).and.&
                                  (abs(Occ(i)-Occ(j))/Occ(i).lt.AuxData%ThrSelAct) ) then

                                  write(6,'(2X,"Discarding nearly degenerate pair ",2I4)')i, j
                            else                                      
                                  if ((Flags%IFlCore.eq.1) .or. (Flags%IFlCore.eq.0.and.i.gt.NI.and.j.gt.NI)) then
                                        if ((abs(Occ(i) + Occ(j) - Two).gt.AuxData%ThrQInact) .and.&
                                              (abs(Occ(i) + Occ(j)).gt.AuxData%ThrQVirt)) then

                                              ind = ind + 1
                                              AuxData%IndX(ind) = ind
                                              AuxData%IndN(1, ind) = i
                                              AuxData%IndN(2, ind) = j

                                              AuxData%IPair(i,j)=1
                                              AuxData%IPair(j,i)=1
                                        end if
                                  end if

                            end if


                      end if

                end do
          end do

          AuxData%NDim = ind
          AuxData%map = zero

          k = 1
          do i = 1, NI+NA+NV
                if(AuxData%IndAux(i)==1)then
                      AuxData%map(i) = k
                      k = k+1
                end if
                if(AuxData%IndAux(i)==2)then
                      AuxData%map(i) = 0
                end if
          end do

          call imsg("Original number of pq pairs", ij)
          call imsg("Reduced number of pq pairs", ind)


          write(*, '(A5, 2A6, 2A14)') '#', 'i', 'j', 'Occ(i)', 'Occ(j)'
          write(*, '(A43)') repeat('-', 43)

          do i = 1, min(100, AuxData%NDim)
                write(*, '(I5, 2I6, 2F15.8)') i, AuxData%IndN(1,i), AuxData%IndN(2,i), &
                      Occ(AuxData%IndN(1,i)), Occ(AuxData%IndN(2,i))
          end do

          select case(Flags%Jobtype)

          case(JOB_TYPE_ERPA)
                select case(Flags%ITwoEl)
                case(1)                          
                      call msg("Starting PhERPA spinres incore")

                      AuxData%HType = H_GPF
                      !                      AuxData%HType =H_DYALL
                      rdm_dump = .false.
                      call pherpa_incore_driver(AuxData, Flags, Ints, IntsD, rdm_dump)

                end select
          case(JOB_TYPE_AC, JOB_TYPE_DUCC)
                select case(Flags%ITwoEl)
                case(1)                          
                      call msg("Starting phAC spinres incore")

                      AuxData%HType = H_GPF
!                      AuxData%HType =H_DYALL
                      call phac_incore_driver(AuxData, Flags, Ints, IntsD)                      
                end select

          end select

          stop
        end associate

  end subroutine acph_driver

  subroutine pherpa_allocate_init(AuxData, Ints)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints

        allocate(Ints%Aints1e_aa(AuxData%NBasis, AuxData%NBasis))
        allocate(AuxData%int_alpha(Ints%Ints2e_dim))
        allocate(AuxData%ABplus(AuxData%NDim, AuxData%NDim))
        allocate(AuxData%ABmin(AuxData%NDim, AuxData%NDim))
        allocate(AuxData%ABfull(2*AuxData%NDim, 2*AuxData%NDim))
        allocate(AuxData%Nmetric(2*AuxData%NDim, 2*AuxData%NDim))
        allocate(AuxData%cc(AuxData%NBasis))
        allocate(AuxData%Eigs(AuxData%NDim))
        allocate(AuxData%Eigvec(AuxData%NDim, AuxData%NDim))
        allocate(AuxData%EigvecX(AuxData%NDim, AuxData%NDim))

        allocate(AuxData%EigvY(AuxData%NDim, AuxData%NDim))
        allocate(AuxData%EigvX(AuxData%NDim, AuxData%NDim))

  end subroutine pherpa_allocate_init

    subroutine pherpa_allocate_init_ducc(AuxData, Ints)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints

        allocate(Ints%Aints2e_aa(Ints%Ints2e_dim))
        allocate(Ints%Aints2e_ab(Ints%Ints2e_dim))

  end subroutine pherpa_allocate_init_ducc


   subroutine pherpa_incore_driver(AuxData, Flags, Ints, Ints_ducc, rdm_dump)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        type(TInts), intent(inout), optional :: Ints_ducc
        type(FlagsData), intent(in) :: Flags
        double precision :: ecorr
        logical, intent(in), optional  :: rdm_dump
        logical :: gen_rdm

        print*, 'herea1'
        if (Flags%JOBTYPE .ne. JOB_TYPE_DUCC)then
              if (Flags%IORCA == 1 .and. Flags%Algorithm == ALG_SPINRES) then
                    ! Integrals already in 4-fold format (ints2e_aa, ints2e_ab)
                    print*, 'tak dobrze'
             else
                   call twoel_8fold_to_ints2e_4fold(Ints, AuxData%NBasis)
             end if
       end if

        call pherpa_allocate_init(AuxData, Ints)

        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then
              call pherpa_allocate_init_ducc(AuxData, Ints)
        end if
        
        AuxData%alpha = one
        AuxData%pherpa_print = .true.
        

        if (present(rdm_dump)) then
              gen_rdm = .true.              
        else
              gen_rdm = .false.
        end if
                
        call pherpa_incore_iter(AuxData, Flags, Ints, ecorr, gen_rdm, Ints_ducc)

        ecorr = ecorr * frac12 ! because  E_corr = frac12 * W^1
        write(*,'(/,1X,''ECASSCF+ENuc, AC1-Corr, ERPA-CASSCF'',6X,3F15.8)'), &
              AuxData%ECAS_calc, ecorr, AuxData%ECAS_calc + ecorr


    end subroutine pherpa_incore_driver

    subroutine phAC_incore_driver(AuxData, Flags, Ints, Ints_ducc)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        type(TInts), intent(inout), optional :: Ints_ducc
        type(FlagsData), intent(in) :: Flags
        double precision :: W, ecorr, Wx, AC0_num
        double precision :: W_0, ecorr_min0
        integer :: ngrid, i
        double precision, dimension(:), allocatable :: XGrid, WGrid


        if (Flags%JOBTYPE .ne. JOB_TYPE_DUCC)then
             if (Flags%IORCA == 1 .and. Flags%Algorithm == ALG_SPINRES) then
                 ! Integrals already in 4-fold format (ints2e_aa, ints2e_ab)
             else
                 call twoel_8fold_to_ints2e_4fold(Ints, AuxData%NBasis)
             end if
        end if
        call pherpa_allocate_init(AuxData, Ints)

        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then
              call pherpa_allocate_init_ducc(AuxData, Ints)
        end if

        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then         

              AuxData%alpha = zero!0.000001!zero
!              AuxData%alpha = frac12 !++!
              call pherpa_incore_iter(AuxData, Flags, Ints, W, .true., Ints_ducc)
              write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X)') &
                    'iter', 'ACalpha', 'wgrid(i)', 'W'

              write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X)') &
                    0, AuxData%alpha,  0.0, W
              AuxData%alpha = 1.d-3 !0.000001!zero
!              AuxData%alpha = frac12 !++!
              call pherpa_incore_iter(AuxData, Flags, Ints, Wx, .true., Ints_ducc)
              AC0_num = frac12*Wx / AuxData%alpha
              print*, 'AC0_num = ', AC0_num
!              stop
              
        else
              AuxData%alpha = zero!0.000001!zero
              AuxData%alpha = frac12 !++!
              call pherpa_incore_iter(AuxData, Flags, Ints, W, .false., Ints_ducc)
              write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X)') &
                    'iter', 'ACalpha', 'wgrid(i)', 'W'
              write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X)') &
                    0, AuxData%alpha,  0.0, W
        end if

        W_0 = W
 !       stop
        ngrid = 50

        
        allocate(xgrid(ngrid))
        allocate(wgrid(ngrid))


        print*, 'ngrid', ngrid

        Call GauLeg(zero, one, xgrid, wgrid, ngrid)

        ecorr = zero
        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then         
              ecorr_min0 = zero
              write(*, '(A5, 13A13)') 'alpha', 'W_aa_act', '-W_aa_act0', 'W_ab_act', '-W_ab_act0', 'W_aa', '-W_aa0', 'W_ab', 'ecorr', 'Waa_act_zwrs', 'Waa_act_dcrs', 'W_aa_act_zw', 'W_aa_act_dc'
              !write(*,'(A5, 9A13)') ' ', 'W_aa_act', 'U_aa_ref_act', 'W_ab_act', 'W_ab_ref_act', 'W_aa', 'W_aa_ref', 'W_ab', 'ecorr', ' '
              !        write(*,'(/,1X,''iter, ACalpha, W, wgrid(i), Ecorr_i,'',6X)')
              !write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X,A15, X,A15, X,A15)') &
              !      'iter', 'ACalpha', 'wgrid(i)', 'W', 'Ecorr_i', 'W-W_0', 'Ecorr_i-W_0'
              !write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X,A15)') &
              !      'iter', 'ACalpha', 'wgrid(i)', 'W', 'Ecorr_i'

              do i = 1, NGrid
                    
                    AuxData%alpha = xgrid(i)
                    call pherpa_incore_iter(AuxData, Flags, Ints, W, .false., Ints_ducc)
                    
                    ecorr = ecorr + W * wgrid(i)
                    !ecorr_min0 = ecorr_min0 + (W-W_0) * wgrid(i)
                    !write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8, 1X,F15.8, 1X,F15.8)') &
                          !      i, AuxData%alpha,  wgrid(i), W, ecorr, W-W_0, ecorr_min0
                    !write(*,'(/,1X,''Ecorr_iter,'',6X, F15.8)') i, AuxData%alpha, W, wgrid(i), ecorr
!                    write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8)') &
                    !i, AuxData%alpha,  wgrid(i), W, ecorr
                    write(*, '(2F15.10, 3I5)')AuxData%alpha,  W, AuxData%AAnegs, AuxData%Apnegs, AuxData%Amnegs

              end do
        else

              write(*,'(/,1X,A5,1X,A15,1X,A15,1X,A15,1X,A15)') &
                    'iter', 'ACalpha', 'wgrid(i)', 'W', 'Ecorr_i'
              do i = 1, NGrid
                    
                    AuxData%alpha = xgrid(i)
                    call pherpa_incore_iter(AuxData, Flags, Ints, W, .false., Ints_ducc)                    
                    ecorr = ecorr + W * wgrid(i)
                    write(*,'(1X,I5,1X,F15.8,1X,F15.8,1X,F15.8,1X,F15.8)') &
                          i, AuxData%alpha,  wgrid(i), W, ecorr
              end do
        end if
              
         if (Flags%JOBTYPE == JOB_TYPE_DUCC)then         
               write(*,'(/,1X,''Eref_DUCC+ENuc, AC-Corr, ERPA-CASSCF'',6X,3F15.8)'), &
                     AuxData%E_ref_ducc, ecorr, AuxData%E_ref_ducc + ecorr
               
               write(*,'(/,1X,''Eref_DUCC+ENuc, AC-Corr--W0, ERPA-CASSCF-W0'',6X,3F15.8)'), &
                     AuxData%E_ref_ducc, ecorr_min0, AuxData%E_ref_ducc + ecorr_min0
         else

               write(*,'(/,1X,''ECASSCF+ENuc, AC-Corr, ERPA-CASSCF'',6X,3F15.8)'), &
                     AuxData%ECAS_calc, ecorr, AuxData%ECAS_calc + ecorr
         end if

  end subroutine phAC_incore_driver


    subroutine pherpa_incore_iter(AuxData, Flags, Ints, ecorr, gen_rdm, IntsD)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        type(TInts), intent(inout), optional :: IntsD
        type(FlagsData), intent(in) :: Flags
        double precision, intent(out) :: ecorr
        logical, intent(in) :: gen_rdm

        integer :: i, nn
        type (tclock) :: timer, timer0
        integer, dimension(:), allocatable :: dy
        double precision, dimension(:), allocatable :: work
        double precision, parameter :: small_e = 1.d-3
        double precision, parameter :: big_e = 1.d+8


        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then
              call erpa_init_incore_ducc(AuxData, Flags, Ints, IntsD)
 !             call write_fcidump_alpha("FCIDUMP-alpha", AuxData, Ints)
!              stop
              call  pherpa_incore_spinres_ducc(AuxData, AuxData%HNO0, &
                    Ints, AuxData%Ndim,  AuxData%Ndim, AuxData%IndN, AuxData%IndN, AuxData%IndX, 1)
              call pherpa_symm_spinres(AuxData)              
              call restore_XY(AuxData)


              ! call  pherpa_incore_spinres_ducc_fullmat(AuxData, AuxData%HNO0, &
              !      Ints, AuxData%Ndim,  AuxData%Ndim, AuxData%IndN, AuxData%IndN, AuxData%IndX, 1)
              ! call pherpa_nonsymm_spinres(AuxData)

              if (abs(AuxData%alpha).lt.1.d-8)then
                    associate (Nbasis => AuxData%NBasis)
                      allocate(AuxData%ph_rdm2_aa(Nbasis, Nbasis, Nbasis, Nbasis))
                      allocate(AuxData%ph_rdm2_ab(Nbasis, Nbasis, Nbasis, Nbasis))
                    end associate
                    !call pherpa_gen_rdm_fluct(AuxData)
              end if
        else
              call erpa_init_incore(AuxData, Flags, Ints)
              !call write_fcidump_alpha_spinres("FCIDUMP-alpha", AuxData, Ints)
              !stop
              call  pherpa_incore_spinres(AuxData, AuxData%HNO0, &
                    Ints, AuxData%Ndim,  AuxData%Ndim, AuxData%IndN, AuxData%IndN, AuxData%IndX, 1)
              call pherpa_symm_spinres(AuxData)
              if (gen_rdm)then
                    associate (Nbasis => AuxData%NBasis)
                      allocate(AuxData%ph_rdm2_aa(Nbasis, Nbasis, Nbasis, Nbasis))
                      allocate(AuxData%ph_rdm2_ab(Nbasis, Nbasis, Nbasis, Nbasis))
                    end associate
                    call restore_XY(AuxData)
                    call pherpa_gen_rdm(AuxData, ecorr, Ints)
                    stop
              end if
        end if


!        if (AuxData%pherpa_print)then
!              Write(6,'(" *** ERPA-CAS SPINRES Excitation Energies (a.u., eV) *** ")')
              allocate(dy(AuxData%NDim))
              allocate(work(AuxData%NDim))
        
              do i = 1, AuxData%NDim
                    dy = i
              end do
              work = AuxData%Eigs
              call dsort(work, dy, AuxData%NDim)
!              print*, 'smllest eigenvalues'
              do i = 1, 10!AuxData%NDim!10
 !                   write(*,'(I4,4X,2F16.6)') i, work(i) , 27.211 * work(i)
                    if (.not.(work(i) > small_e .and. work(i) < big_e)) then
  !                        print*, 'this will be skipped'
                    end if
              end do
              deallocate(dy)
              deallocate(work)

              !        end if

        
        if (Flags%JOBTYPE == JOB_TYPE_DUCC)then
              call pherpa_energy_ducc(AuxData, ecorr, Ints, IntsD)
        else if (Flags%JOBTYPE == JOB_TYPE_AC)then
              call pherpa_energy(AuxData, ecorr, Ints)
        end if
        
  end subroutine pherpa_incore_iter


  subroutine erpa_init_incore(AuxData,  Flags, Ints)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        type(FlagsData), intent(in) :: Flags

        integer :: i, j, k,  t, l
        integer :: i0, i1
        double precision :: temp
        integer(I8) :: nb, npair, pq, rs
        integer :: p,q,r,s
        integer(I8) :: a1,b1,a2,b2, a,b, c,d, idx

!        allocate(Ints%Aints1e_aa(AuxData%NBasis, AuxData%NBasis))
!        allocate(AuxData%int_alpha(Ints%Ints2e_dim))

        associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, &
              NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
              NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)

          ! ===============================
          !  one-electron part 
          ! ===============================
!          print*, 'alpha', alpha
          Aints1e_aa = alpha * AuxData%HNO0

          do i = 1, NBasis
                do j = 1, NBasis
                      if (IAux(i) == IAux(j)) then
                            Aints1e_aa(i, j) = Aints1e_aa(i, j) + (One - alpha) * AuxData%HNO0(i,j)

                            if (AuxData%HType == H_DYALL) then
                            !      print*, 'dyal'
                                  if (IAux(i) == 1) then
                                        i0 = 1
                                        i1 = NI
                                  else  
                                        i0 = 1
                                        i1 = NIA
                                  end if
                            else if (AuxData%HType == H_GPF) then
                                  if (IAux(i) == 0) then
                                        i0 = NI + 1
                                        i1 = NIA
                                  else if (IAux(i) == 1) then
                                        i0 = 1
                                        i1 = NI
                                  else
                                        i0 = 1
                                        i1 = NIA
                                  end if
                            end if

                            temp = zero
                            do t = i0, i1
                                  temp = temp + occ(t) * &
                                        ( ints2e_aa(gmap_4fold(t,t,i,j, NBasis)) &
                                        + ints2e_ab(gmap_4fold(t,t,i,j, NBasis)) &
                                        - ints2e_aa(gmap_4fold(t,i,t,j, NBasis)) )
                            end do

                            Aints1e_aa(i, j) = Aints1e_aa(i, j) + (One - alpha) * temp
!                            if(abs(Aints1e_aa(i,j)).gt.1.d-5)then
!                                  write(*,'(2I3, 2F12.6)') i, j, Aints1e_aa(i,j), Aints1e_aa(j,i)!AuxData%HNO0(i,j), temp!, AuxData%HNO0(i,j)
!                            end if
                      end if
                end do
          end do
!          print*, ''
          ! do i = 1, NBasis
          !       do j = 1, NBasis
          !             if(abs(Aints1e_aa(i,j)).gt.1.d-5)then
          !                   write(*,'(2I3, 2F12.6)') i, j, Aints1e_aa(i,j), Aints1e_aa(j,i)!, AuxData%HNO0(i,j)
          !             end if
          
          !       end do
          ! end do

          ! ===============================
          !  two-electron mask (int_alpha)
          ! ===============================
          
          nb    = int(NBasis, I8)
          npair = nb*nb

          AuxData%int_alpha = 1

          do pq = 1_I8, npair
                q = int((pq-1_I8)/nb + 1_I8)
                p = int(pq - int(q-1,I8)*nb)

                do rs = pq, npair  
                      s = int((rs-1_I8)/nb + 1_I8)
                      r = int(rs - int(s-1,I8)*nb)

                      a1 = pq; b1 = rs
                      a  = min(a1,b1); b = max(a1,b1)

                      a2 = int(q,I8) + (int(p,I8)-1_I8)*nb
                      b2 = int(s,I8) + (int(r,I8)-1_I8)*nb
                      c  = min(a2,b2); d = max(a2,b2)

                      if ( (c < a) .or. (c == a .and. d < b) ) cycle  

                      idx = gmap_4fold(r,s,p,q,NBasis)

                      if (AuxData%HType == H_DYALL) then
                            if (.not. (IAux(p) == 1 .and. IAux(q) == 1 .and. &
                                  IAux(r) == 1 .and. IAux(s) == 1)) then
                                  AuxData%int_alpha(idx) = 0
                            end if

                      else if (AuxData%HType == H_GPF) then
                            if (.not. (IAux(p) == IAux(q) .and. &
                                  IAux(q) == IAux(r) .and. &
                                  IAux(r) == IAux(s))) then
                                  AuxData%int_alpha(idx) = 0
                            end if
                      end if
                end do
          end do

        end associate
  end subroutine erpa_init_incore


  subroutine erpa_init_incore_ducc(AuxData,  Flags, Ints, Ints_ducc)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        type(TInts), intent(inout), optional :: Ints_ducc
        type(FlagsData), intent(in) :: Flags

        integer :: i, j, k,  t, l
        integer :: i0, i1
        double precision :: temp
        integer :: nb, npair, pq, rs
        integer :: p,q,r,s
        integer :: r_g, s_g
        integer :: a1,b1,a2,b2, a,b, c,d, idx
        double precision :: val_full, val_ducc

        associate(Occ=>AuxData%Occ, ENuc=>AuxData%ENuc, &
              NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, &
              NV=>AuxData%NV, NBasis=>AuxData%NBasis, IAux=>AuxData%IndAux, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              Aints2e_aa=>Ints%Aints2e_aa, Aints2e_ab=>Ints%Aints2e_ab, &
              ints2e_ducc_aa=>Ints_ducc%ints2e_aa, ints2e_ducc_ab=>Ints_ducc%ints2e_ab, &
              ints1e_ducc_aa=>Ints_ducc%ints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha, &
              g_aa=>AuxData%rdm1_p, g_bb=>AuxData%rdm1_m)

         
          ! ===============================
          !  Two-electron part mixing
          ! ===============================
          
          ! Assuming Aints2e_aa/ab are allocated
          
          nb    = NBasis
          npair = nb*nb
!          print*, 'calki Aints2e_aa'
          do pq = 1, npair
                q = int((pq-1)/nb + 1)
                p = int(pq - int(q-1)*nb)

                do rs = pq, npair  
                      s = int((rs-1)/nb + 1)
                      r = int(rs - int(s-1)*nb)

                      a1 = pq; b1 = rs
                      a  = min(a1,b1); b = max(a1,b1)

                      a2 = q + (p-1)*nb
                      b2 = s + (r-1)*nb
                      c  = min(a2,b2); d = max(a2,b2)

                      if ( (c < a) .or. (c == a .and. d < b) ) cycle  

                      idx = gmap_4fold(r,s,p,q,NBasis)
                      
                      ! AA Mixing
                      val_full = ints2e_aa(idx)
                      if (IAux(p)==1 .and. IAux(q)==1 .and. IAux(r)==1 .and. IAux(s)==1) then
                           val_ducc = ints2e_ducc_aa(gmap_4fold(r-NI, s-NI, p-NI, q-NI, NA))
                           Aints2e_aa(idx) = (One - alpha) * val_ducc + alpha * val_full
                      else
                           Aints2e_aa(idx) = alpha * val_full
                           if (IAux(p)>=2 .and. IAux(q)>=2 .and. IAux(r)>=2 .and. IAux(s)>=2) then ! VVVV block
                                Aints2e_aa(idx) = val_full
                           end if
                      end if

                      ! AB Mixing
                      val_full = ints2e_ab(idx)
                      if (IAux(p)==1 .and. IAux(q)==1 .and. IAux(r)==1 .and. IAux(s)==1) then
                           val_ducc = ints2e_ducc_ab(gmap_4fold(r-NI, s-NI, p-NI, q-NI, NA))
                           Aints2e_ab(idx) = (One - alpha) * val_ducc + alpha * val_full
                      else
                           Aints2e_ab(idx) = alpha * val_full
                           if (IAux(p)>=2 .and. IAux(q)>=2 .and. IAux(r)>=2 .and. IAux(s)>=2) then ! VVVV block
                                Aints2e_ab(idx) = val_full
                           end if
                      end if

                      
                end do
          end do


          ! ===============================
          !  one-electron part 
          ! ===============================

          ! print*, 'calki ints1e_aa'
          
          ! do i = 1, NBasis
          !       do j = 1, NBasis
          !             if (abs(ints1e_aa(i,j)).gt.1.d-8)then
          !                   write(*,'(ES20.10, 2I5)') ints1e_aa(i,j), i, j
          !             end if
          !       end do
          ! end do

          ! print*, 'calki ints1e_aa_ducc'
          
          ! do i = 1, NA
          !       do j = 1, NA
          !             if (abs(ints1e_ducc_aa(i,j)).gt.1.d-8)then
          !                   write(*,'(F15.10, 2I5)') ints1e_ducc_aa(i,j), i, j
          !             end if
          !       end do
          ! end do
!          print*, 'calki aints', alpha

          Aints1e_aa = alpha * ints1e_aa
          
          do p = 1, NBasis
                do q = 1, NBasis
                     if (IAux(p)==1 .and. IAux(q)==1) then ! AA block
                          Aints1e_aa(p,q) = (One - alpha) * ints1e_ducc_aa(p-NI, q-NI) + alpha * ints1e_aa(p,q)
                     else if (IAux(p)>=2 .and. IAux(q)>=2) then ! VV block
                          Aints1e_aa(p,q) = ints1e_aa(p,q)
                          
                          ! Add contraction terms for VV block
                          ! Term 1: 0.5*(1-alpha) * sum_rs g_aa(r,s) * (pq|rs)_aa
                          ! Term 2: (1-alpha) * sum_rs g_bb(r,s) * (pq|rs)_ab
                          ! Term 3: -0.5*(1-alpha) * sum_rs g_aa(s,r) * (pr|qs)_aa  -> (ps|qr) ? -> (pr|qs) = <pr|qs> = (ps|rq)?
                          ! Python: np.einsum('sr,prqs->pq', g_aa, ints2e_final[0][V,A,V,A])
                          ! prqs -> (pr|qs)_aa. 
                          ! In gmap_4fold(r,s,p,q), indices are p,q,r,s. (pq|rs) = <pq|rs>.
                          ! 2e integral notation: (ij|kl) = <ik|jl>. 
                          ! Python script uses Mulliken (pq|rs).
                          ! Fortran ints2e_aa(gmap_4fold(p,q,r,s)) corresponds to (pq|rs).
                          
                          temp = zero
                          do r = 1, NA ! active indices
                               do s = 1, NA
                                     ! r,s are active indices relative to NActive (1..NA)
                                     ! Need global indices for integrals
                                     r_g = r + NI
                                     s_g = s + NI
                                     
                                     ! Term 1
                                     temp = temp + 0.5d0 * (One - alpha) * g_aa(r,s) * ints2e_aa(gmap_4fold(p,q,r_g,s_g, NBasis))
                                     
                                     ! Term 2
                                     temp = temp + (One - alpha) * g_bb(r,s) * ints2e_ab(gmap_4fold(p,q,r_g,s_g, NBasis))
                                     
                                     ! Term 3
                                     ! (pr|qs)_aa = ints2e_aa(gmap_4fold(i,r_g,j,s_g))  ? No, gmap expects (p,q,r,s).
                                     ! (pr|qs) is Mulliken notation for <pq|rs> (Dirac).
                                     ! Wait, standard chemistry notation: (pq|rs) = \int p*(1) q(1) r*(2) s(2).
                                     ! Python: ints2e_final[0][V,A,V,A] -> indicies p,r,q,s corresponding to (pr|qs)? 
                                     ! pyscf/gammcor integrals are usually (pq|rs).
                                     ! The python code accesses [V,A,V,A]. i=V, r=A, j=V, s=A.
                                     ! so it retrieves (ir|js).
                                     ! And contracts with g_aa(s,r).
                                     ! So term is sum_rs g_aa(s,r) * (ir|js).
                                     ! Using gmap_4fold(i,r_g, j,s_g, NBasis).
                                     
                                     temp = temp - 0.5d0 * (One - alpha) * g_aa(s,r) * ints2e_aa(gmap_4fold(p,r_g,q,s_g, NBasis))
                               end do
                          end do
                          Aints1e_aa(p,q) = Aints1e_aa(p,q) + temp

                    end if
                    ! if (abs(Aints1e_aa(p,q)).gt.1.d-8)then
                    !       write(*,'(ES22.14, 4I5)') Aints1e_aa(p,q), p, q, 0, 0
                    ! end if

                end do
          end do
          
          ! do i = 1, NBasis
          !       do j = 1, NBasis
          !             if (abs(Aints1e_aa(i, j)).gt.1.d-9)then
          !                   write(*, '(A5, 2I3, F20.12)')'Ain', i, j, Aints1e_aa(i, j)
          !             end if
          !       end do
          ! end do



          ! ===============================
          !  two-electron mask (int_alpha)
          ! ===============================
          
          AuxData%int_alpha = 1

          do pq = 1, npair
                q = int((pq-1)/nb + 1)
                p = int(pq - int(q-1)*nb)

                do rs = pq, npair  
                      s = int((rs-1)/nb + 1)
                      r = int(rs - int(s-1)*nb)

                      a1 = pq; b1 = rs
                      a  = min(a1,b1); b = max(a1,b1)

                      a2 = q + (p-1)*nb
                      b2 = s + (r-1)*nb
                      c  = min(a2,b2); d = max(a2,b2)

                      if ( (c < a) .or. (c == a .and. d < b) ) cycle  

                      idx = gmap_4fold(r,s,p,q,NBasis)

                      if (AuxData%HType == H_DYALL) then
                            if (.not. (IAux(p) == 1 .and. IAux(q) == 1 .and. &
                                  IAux(r) == 1 .and. IAux(s) == 1)) then
                                  AuxData%int_alpha(idx) = 0
                            end if

                      else if (AuxData%HType == H_GPF) then
                            if (.not. (IAux(p) == IAux(q) .and. &
                                  IAux(q) == IAux(r) .and. &
                                  IAux(r) == IAux(s))) then
                                  AuxData%int_alpha(idx) = 0
                            end if
                      end if
                end do
          end do

        end associate
  end subroutine erpa_init_incore_ducc

  subroutine pherpa_incore_spinres(AuxData, h, Ints, Ndim1, &
        Ndim2, indn1, indn2, indx, multiply_by_S)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        double precision, intent(in) :: h(:,:)
        integer, intent(in) :: NDim1, Ndim2
        integer, dimension(:,:), intent(in) ::indn1, indn2
        integer, dimension(:), intent(in) ::indx
        integer, intent(in) :: multiply_by_S


        integer :: i, j, t, u, v, w
        integer p, q, r, s, rs, pq
        integer :: pp, qq, rr, ss
        double precision :: Arspq
        double precision, dimension(:,:), allocatable :: auxi_a, auxio_a
        double precision, dimension(:,:), allocatable :: auxi_ab, auxio_ab
        double precision :: scale

        double precision :: T0, T1, T2, T3, T4, T5, T6, TW1, TW2
        double precision :: T1a, T1b, T2a, T2b, T0a, T0b
        double precision, dimension(:,:), allocatable:: wmat
        double precision :: num1, num2, num3
        double precision :: num_f, num_h
        double precision :: nnn, nieb
        double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
        type (tclock) :: timer0, timer1, timerxx
        real(8) :: den, denmin
        integer :: imin, jmin

        associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, &
              NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, &
              IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, &
              rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, &
              rdm2_mp=>AuxData%rdm2_mp, &
              rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)

!          allocate(AuxData%ABplus(NDim1, NDim2))
!          allocate(AuxData%ABmin(NDim1, NDim2))
!          allocate(AuxData%cc(NBasis))


          ! allocate(AuxData%auxi_a (NBasis, NBasis))
          ! allocate(AuxData%auxi_ab(NBasis, NBasis))

          ! allocate(AuxData%auxio_a (NBasis, NBasis))
          ! allocate(AuxData%auxio_ab(NBasis, NBasis))

          ! AuxData%auxi_a   = zero
          ! AuxData%auxio_a  = zero
          ! AuxData%auxi_ab  = zero
          ! AuxData%auxio_ab = zero

          ! ===============================
          !   AuxI pq A =  ∑ N occ​nt *scale * ((pq∣tt)_aa − (pt∣qt)_aa
          ! ===============================

          allocate(wmat(NBasis, NBasis))

          associate(c=>AuxData%cc)
          
          do i = 1, NBasis
                c(i) = sqrt(n(i))
                if (n(i) .lt. frac12) c(i) = -c(i)
          end do


          wmat = zero

          do p = 1, NBasis
                do r = 1, NIA
                do t = 1, NIA
                   do v = 1, NIA
                         do u = 1, NIA
                               num1 = one
                               if (AuxData%int_alpha(gmap_4fold(u,v,t,p, NBasis)) == 0) num1 = alpha
                               wmat(p, r) = wmat(p, r) + num1 * (&
                                     erdm_ppx(rdm2_pp, n_p, r, v, t, u, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, v, t, p, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n_p, n_m, r, v, t, u, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(u, v, t, p, NBasis)) +&
                                     
                                     erdm_ppx(rdm2_pp, n_p, r, v, u, t, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, p, t, v, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n_p, n_m, r, v, u, t, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(v, t, p, u, NBasis)))                             
                       end do
                 end do
           end do
           
            end do
      end do

      !       print *, 'asym wmat = ', sqrt(sum((wmat - transpose(wmat))**2))
      ! print *, 'norm wmat = ', sqrt(sum(wmat*wmat))
      ! print *, 'rel asym wmat = ', sqrt(sum((wmat - transpose(wmat))**2)) / sqrt(sum(wmat*wmat))
      ! print *, 'wmat(2,1), wmat(1,2)=', wmat(2,1), wmat(1,2)


           AuxData%ABplus = zero
           AuxData%ABmin = zero

          i_rowloops: do i = 1, NDim1
            j_colloops: do j = 1, NDim2
                r = indn1(1, i)
                s = indn1(2, i)
                rs = indx(i)
                
                pp = indn2(1, j)
                qq = indn2(2, j)
                pq = indx(j)

                      
                do p = qq, pp, pp - qq
                do q = qq, pp, pp - qq
                      if ( p .ne. q ) then
                                        
                            Arspq = zero

                            T0 = zero
                            T0a = zero
                            T0b = zero
                            if ( p == r ) T0a = T0a + (n_p(p) - n_p(s)) * Aints1e_aa(q, s)
                            if ( s == q ) T0b = T0b + (n_p(q) - n_p(r)) * Aints1e_aa(p, r)
                            T0 = T0a+T0b
                            !                            if (abs(T0).gt.1.d-5)then
                            ! if((p==7.and.q==2.and.r==7.and.s==3).or.(p==7.and.q==3.and.r==7.and.s==2))then
                            !       write(*, '(A10, 4I5, 9F12.7 )') 'T0', r, s, p, q,  T0, n_p(p), n_p(s), n_p(q), n_p(r), &
                            !             Aints1e_aa(q,s), Aints1e_aa(p,r), 
                            ! end if

                            T1a = zero; T1b = zero
                            T1a = T1a + func_T1(AuxData, Ints, r, p, s, q)
                            T1b = T1b + func_T2(AuxData, Ints, r, p, s, q)
                            T1 =  (T1a + T1b)

                            T2a = zero; T2b = zero
                            T2a = T2a + func_T1(AuxData, Ints, s, q, r, p)
                            T2b = T2b + func_T2(AuxData, Ints, s, q, r, p)
                            T2 = (T2a + T2b)
                            
                            T5 = zero; T6 = zero
                            T6 = T6 - func_T34(AuxData, Ints, s, p, r, q)
                            T5 = T5 - func_T34(AuxData, Ints, q, r, p, s)

                            !if (abs(T5).gt.1.d-5)then
                            !      write(*, '(A10, 4I5, F20.15 )') 'niebieski', r, s, p, q,  T5
                            !end if
                            
                            TW1 = zero; TW2 = zero
                            if (s == q) TW1 = TW1 - frac12 * wmat(p, r)
                            if (p == r) TW2 = TW2 - frac12 * wmat(q, s)
                            
                            arspq = T0 + T1 + T2 + T6 + T5 + TW1 + TW2
!                             if ((rs == 3 .and. pq == 2) .or. (rs == 2 .and. pq == 3)) then
!     write(*,'(A,2I4,A,4I4,8ES16.6)') 'DEBUG arspq ', rs, pq, ' | r s p q = ', &
!         r, s, p, q, T0, T1, T2, T6, T5, TW1, TW2, arspq
! end if

!                            if (abs(arspq).gt.1.d-1)then
!                                  write(*, '(A10, 2I5, A1, 4I3, 8F10.5)')'arspq', rs, pq, '|', r,s, p,q, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
!                            end if
                            if ( r .gt. s .and. p .gt. q ) then
                                  !if((rs==3.and.pq==12).or.(rs==12.and.pq==3))then
                                  !if((rs==16.and.pq==20).or.(rs==20.and.pq==16))then
!                                  if (abs(arspq).gt.1.d-5)then
!                                        if (.not.(p==r.and.q==s))then
                                     !         write(*, '(A10, 2I3, A1, 4I3, 8F10.5)')'arspq', rs, pq, '|', r,s, p,q, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
 !                     write(*, '(A10, 4I3, 8F10.5)')'sniez', r,s, p,q, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
  !                                      end if
   !                                     end if

                                  AuxData%ABplus(rs, pq) = AuxData%ABplus(rs, pq) + arspq
                                  AuxData%ABmin (rs, pq) = AuxData%ABmin (rs, pq) + arspq
                            end if

                            if ( r .gt. s .and. q .gt. p ) then
!                                  if((rs==3.and.pq==12).or.(rs==12.and.pq==3))then
                                  !                                        if((rs==16.and.pq==20).or.(rs==20.and.pq==16))then
                            !write(*, '(A10, 2I3, A1, 4I3, 8F10.5)')'ars_qp', rs, pq, '|', r,s, q,p, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
!         end if
 !                                 if (abs(arspq).gt.1.d-5)then
  !                                      if (r>6.and.s>6.and.p>6.and.q>6)then
   !             write(*, '(A10, 4I3, 8F10.5)')'wieksze', r, s, p,q, arspq
    !      end if
   ! end if

                                  AuxData%ABplus(rs, pq) = AuxData%ABplus(rs, pq) + arspq
                                  AuxData%ABmin (rs, pq) = AuxData%ABmin (rs, pq) - arspq
                            end if

                      end if
                end do
          end do
    end do j_colloops
    end do i_rowloops
!call print_asym('BEFORE SCALE ABPLUS', AuxData%abplus, NDim1)
!call print_asym('BEFORE SCALE ABMIN ', AuxData%abmin,  NDim1)

    do i = 1, NDim1
          do j = 1, NDim2

                r = indn1(1, i); s = indn1(2, i); rs = indx(i)
                p = indn2(1, j); q = indn2(2, j); pq = indx(j)

                if ( (c(p) + c(q)) * (c(r) + c(s)) /= zero ) then
                      AuxData%abplus(rs, pq) = AuxData%abplus(rs, pq) / ((c(p) + c(q)) * (c(r) + c(s)))
                end if

                if ( (c(p) - c(q)) * (c(r) - c(s)) /= zero ) then
                      AuxData%abmin(rs, pq) = AuxData%abmin(rs, pq) / ((c(p) - c(q)) * (c(r) - c(s)))
                end if
          end do
    end do
!        call print_asym('AFTER SCALE ABPLUS', AuxData%abplus, NDim1)
!        call print_asym('AFTER SCALE ABMIN ', AuxData%abmin,  NDim1)


!         denmin = huge(1.0d0)

! do i = 1, NDim1
!     do j = 1, NDim2
!         r = indn1(1,i); s = indn1(2,i)
!         p = indn2(1,j); q = indn2(2,j)

!         den = abs((c(p)-c(q))*(c(r)-c(s)))

!         if (den > 0.0d0 .and. den < denmin) then
!             denmin = den
!             imin = i
!             jmin = j
!         end if
!     end do
! end do

! print *, 'min nonzero den_min = ', denmin
! print *, 'at i,j = ', imin, jmin


    !         do i = 1, NDim1
    !       do j = 1, NDim2
    !             r = indn1(1, i); s = indn1(2, i); rs = indx(i)
    !             p = indn2(1, j); q = indn2(2, j); pq = indx(j)

    !             nieb = AuxData%abplus(rs, pq) - AuxData%abplus(pq, rs)
    !             if (abs(nieb).gt.1.d-3)then
    !                   if (r>6.or.s>6.or.p>6.or.q>6)then
    !          write(*, '(A7, 2I3, A1, 4I3, 3F10.5)')'abpls', i, j, '|', r,s,p,q,AuxData%abplus(rs, pq), AuxData%abplus(pq, rs), nieb
    !                   end if
    !             end if
    !       end do
    ! end do


    ! do i = 1, NDim1
    !       do j = 1, NDim2

    !             r = indn1(1, i); s = indn1(2, i); rs = indx(i)
    !             p = indn2(1, j); q = indn2(2, j); pq = indx(j)
               
    !             if (abs(AuxData%abplus(rs, pq)).gt.1.d-5)then
    !                   if (IAux(p)==2.and.IAux(q)==1.and.IAux(r)==1.and.IAux(s)==1)then
    !                         write(*, '(A7, 6I3, F10.5)') 'abpl', i, j, p, q, r, s, AuxData%abplus(rs, pq)
    !                   end if
    !             end if
    !       end do
    ! end do


    deallocate(wmat)
  end associate


  end associate
  
end subroutine pherpa_incore_spinres

  subroutine pherpa_incore_spinres_ducc(AuxData, h, Ints, Ndim1, &
        Ndim2, indn1, indn2, indx, multiply_by_S)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        double precision, intent(in) :: h(:,:)
        integer, intent(in) :: NDim1, Ndim2
        integer, dimension(:,:), intent(in) ::indn1, indn2
        integer, dimension(:), intent(in) ::indx
        integer, intent(in) :: multiply_by_S


        integer :: i, j, t, u, v, w
        integer p, q, r, s, rs, pq
        integer :: pp, qq, rr, ss
        double precision :: Arspq
        double precision, dimension(:,:), allocatable :: auxi_a, auxio_a
        double precision, dimension(:,:), allocatable :: auxi_ab, auxio_ab
        double precision :: scale

        double precision :: T0, T1, T2, T3, T4, T5, T6, TW1, TW2
        double precision :: T1a, T1b, T2a, T2b
        double precision, dimension(:,:), allocatable:: wmat
        double precision :: num1, num2, num3
        double precision :: num_f, num_h
        double precision :: nnn, nieb
        double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
        type (tclock) :: timer0, timer1, timerxx
        real(8) :: den, denmin
        integer :: imin, jmin

        associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, &
              NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, &
              IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, &
              rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, &
              rdm2_mp=>AuxData%rdm2_mp, &
              rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, &
              ints2e_aa=>Ints%Aints2e_aa, ints2e_ab=>Ints%Aints2e_ab, &
              Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)


          allocate(wmat(NBasis, NBasis))

          ! print*, 'inside incore ducc'
          ! do i = 1, Nbasis
          !       write(*, '(A5, I3, 3F12.6)') 'nnn', i, n(i), n_p(i), n_m(i)
          ! end do
          ! print*, ''
          ! do i = 1, NA
          !       do j = 1, NA
          !             write(*, '(A5, 2I3, 2F12.6)') 'rdm', i, j, rdm1_p(i, j), rdm1_m(i, j)
          !       end do
          ! end do

          ! print*, 'auxdata%cc', auxdata%cc
          associate(c=>AuxData%cc)
          
          do i = 1, NBasis
                c(i) = sqrt(n(i))
                if (n(i) .lt. frac12) c(i) = -c(i)
          end do

          num1 = one
          wmat = zero

          do p = 1, NBasis
                do r = 1, NIA
                do t = 1, NIA
                   do v = 1, NIA
                         do u = 1, NIA

                               wmat(p, r) = wmat(p, r) + num1 * (&
                                     frac12 * erdm_ppx(rdm2_pp, n, r, v, t, u, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, v, t, p, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n, n, r, v, t, u, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(u, v, t, p, NBasis)) +&
                                   
                                     frac12 * erdm_ppx(rdm2_pp, n, r, v, u, t, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, p, t, v, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n, n, r, v, u, t, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(v, t, p, u, NBasis)))

                               
                       end do
                 end do
           end do
           
            end do
      end do

      ! print *, 'asym wmat = ', sqrt(sum((wmat - transpose(wmat))**2))
      ! print *, 'norm wmat = ', sqrt(sum(wmat*wmat))
      ! print *, 'rel asym wmat = ', sqrt(sum((wmat - transpose(wmat))**2)) / sqrt(sum(wmat*wmat))
      ! print *, 'wmat(2,1), wmat(1,2)=', wmat(2,1), wmat(1,2)

           AuxData%ABplus = zero
           AuxData%ABmin = zero

          i_rowloops: do i = 1, NDim1
            j_colloops: do j = 1, NDim2
                r = indn1(1, i)
                s = indn1(2, i)
                rs = indx(i)
                
                pp = indn2(1, j)
                qq = indn2(2, j)
                pq = indx(j)

                do p = qq, pp, pp - qq
                      do q = qq, pp, pp - qq

                            if ( p .ne. q ) then
                                        
                            Arspq = zero

                            T0 = zero
                            if ( p == r ) T0 = T0 + (n(p) - n(s)) * Aints1e_aa(q, s)
                            if ( s == q ) T0 = T0 + (n(q) - n(r)) * Aints1e_aa(p, r)

                            ! if ((rs == 3 .and. pq == 2) .or. (rs == 2 .and. pq == 3)) then                            
                            !       write(*,'(A,2I4,A,4I4,10F12.6)') 'DEBUG Tw0pr ', rs, pq, ' | r s p q = ', &
                            !             r, s, p, q, &
                            !             n(p) , n(s) , Aints1e_aa(q, s), (n(p) - n(s)) * Aints1e_aa(q, s)
                            !       write(*,'(A,2I4,A,4I4,10F12.6)') 'DEBUG Tw0sq ', rs, pq, ' | r s p q = ', &
                            !             r, s, p, q, &
                            !             n(q) , n(r) , Aints1e_aa(p, r), (n(q) - n(r)) * Aints1e_aa(p, r)

                            ! end if
                            T1a = zero; T1b = zero
                            T1a = T1a + func_T1_ducc(AuxData, Ints, r, p, s, q)
                            T1b = T1b + func_T2_ducc(AuxData, Ints, r, p, s, q)
                            T1 =  (T1a + T1b)

                            T2a = zero; T2b = zero
                            T2a = T2a + func_T1_ducc(AuxData, Ints, s, q, r, p)
                            T2b = T2b + func_T2_ducc(AuxData, Ints, s, q, r, p)
                            T2 = (T2a + T2b)
                            
                            T5 = zero; T6 = zero
                            T6 = T6 - func_T34_ducc(AuxData, Ints, s, p, r, q)
                            T5 = T5 - func_T34_ducc(AuxData, Ints, q, r, p, s)

                            
                            TW1 = zero; TW2 = zero
                            if (s == q) TW1 = TW1 - frac12 * wmat(p, r)
                            if (p == r) TW2 = TW2 - frac12 * wmat(q, s)

                            arspq = T0 + T1 + T2 + T6 + T5 + TW1 + TW2

!                            if ((rs == 3 .and. pq == 2) .or. (rs == 2 .and. pq == 3)) then
!    write(*,'(A,2I4,A,4I4,8F12.6)') 'DEBUG arspq ', rs, pq, ' | r s p q = ', &
!        r, s, p, q, T0, T1, T2, T6, T5, TW1, TW2, arspq
!end if

! if ((rs == 3 .and. pq == 2) .or. (rs == 2 .and. pq == 3)) then
!     write(*,'(A,2I4,A,4I4,8ES16.6)') 'DEBUG onebody ', rs, pq, ' | r s p q = ', &
!         r, s, p, q, &
!         Aints1e_aa(q,s), Aints1e_aa(p,r), &
!         wmat(p,r), wmat(q,s), &
!         n_p(p), n_p(q), n_p(r), n_p(s)
! end if
                            !if (abs(arspq).gt.1.d-1)then
                            !      write(*, '(A10, 2I5, A1, 4I3, 8F10.5)')'arspq', rs, pq, '|', r,s, p,q, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
                            !end if


                            if ( r .gt. s .and. p .gt. q ) then
                                  ! if((rs==5.and.pq==6).or.(rs==6.and.pq==5))then
                                  !       write(*, '(A10, 2I3, 8F10.5)')'niebieski1', rs, pq, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
                                  ! end if

                                  AuxData%ABplus(rs, pq) = AuxData%ABplus(rs, pq) + arspq
                                  AuxData%ABmin (rs, pq) = AuxData%ABmin (rs, pq) + arspq
                            end if

                            if ( r .gt. s .and. q .gt. p ) then
                                  ! if (abs(arspq).gt.1.d-5)then
                                  !       if(r>6.or.s>6.or.p>6.or.q>6)then
                                  !             write(*, '(A10, 4I3, F10.5)') 'wieksze', r, s, p, q, arspq
                                  !       end if
                                  ! end if
                                  ! if((rs==5.and.pq==6).or.(rs==6.and.pq==5))then
                                  !       write(*, '(A10, 2I3, 8F10.5)')'niebieski2', rs, pq, T0 , T1 , T2 , T6 , T5 , TW1 , TW2, arspq
                                  ! end if

                                  AuxData%ABplus(rs, pq) = AuxData%ABplus(rs, pq) + arspq
                                  AuxData%ABmin (rs, pq) = AuxData%ABmin (rs, pq) - arspq
                            end if


!                                                         if ((rs == 3 .and. pq == 2) .or. (rs == 2 .and. pq == 3)) then
!     write(*,'(A,2I4,A,4I4,10F12.6)') 'DEBUG arspq ', rs, pq, ' | r s p q = ', &
!         r, s, p, q, T0, T1, T2, T6, T5, TW1, TW2, arspq, AuxData%ABplus(rs, pq), AuxData%ABmin(rs, pq)
! end if

                      end if
                end do
          end do
    end do j_colloops
end do i_rowloops


!call print_asym('BEFORE SCALE ABPLUS', AuxData%abplus, NDim1)
!call print_asym('BEFORE SCALE ABMIN ', AuxData%abmin,  NDim1)

!     i = 3
! j = 2

! r = indn1(1,i)
! s = indn1(2,i)
! p = indn2(1,j)
! q = indn2(2,j)

! print *, 'TEST abmin asym i,j=', i,j
! print *, 'r,s,p,q=', r,s,p,q
! print *, 'c(r),c(s),c(p),c(q)=', c(r),c(s),c(p),c(q)
! print *, 'den_min(i,j)=', (c(p)-c(q))*(c(r)-c(s))
! print *, 'abmin(i,j)=', AuxData%abmin(i,j)


    do i = 1, NDim1
          do j = 1, NDim2

                r = indn1(1, i); s = indn1(2, i); rs = indx(i)
                p = indn2(1, j); q = indn2(2, j); pq = indx(j)

                ! if (i==2.and.j==3)then
                !       print*, 'ij23-przed', AuxData%abmin(rs, pq), ((c(p) - c(q)) * (c(r) - c(s)))
                ! end if

                ! if (i==3.and.j==2)then
                !       print*, 'ij32-przed', AuxData%abmin(rs, pq), ((c(p) - c(q)) * (c(r) - c(s)))
                ! end if

                if ( (c(p) + c(q)) * (c(r) + c(s)) /= zero ) then
                      AuxData%abplus(rs, pq) = AuxData%abplus(rs, pq) / ((c(p) + c(q)) * (c(r) + c(s)))
                end if

                if ( (c(p) - c(q)) * (c(r) - c(s)) /= zero ) then
                      AuxData%abmin(rs, pq) = AuxData%abmin(rs, pq) / ((c(p) - c(q)) * (c(r) - c(s)))
                end if
                ! if (i==2.and.j==3)then
                !       print*, 'ij23-po', AuxData%abmin(rs, pq)
                ! end if

                ! if (i==3.and.j==2)then
                !       print*, 'ij32-po', AuxData%abmin(rs, pq)
                ! end if

          end do
    end do

!     denmin = huge(1.0d0)

! do i = 1, NDim1
!     do j = 1, NDim2
!         r = indn1(1,i); s = indn1(2,i)
!         p = indn2(1,j); q = indn2(2,j)

!         den = abs((c(p)-c(q))*(c(r)-c(s)))

!         if (den > 0.0d0 .and. den < denmin) then
!             denmin = den
!             imin = i
!             jmin = j
!         end if
!     end do
! end do

! print *, 'min nonzero den_min = ', denmin
! print *, 'at i,j = ', imin, jmin


! i = 2
! j = 3

! r = indn1(1,i)
! s = indn1(2,i)
! p = indn2(1,j)
! q = indn2(2,j)

! print *, 'TEST abmin asym j,i=', i,j
! print *, 'r,s,p,q=', r,s,p,q
! print *, 'c(r),c(s),c(p),c(q)=', c(r),c(s),c(p),c(q)
! print *, 'den_min(j,i)=', (c(p)-c(q))*(c(r)-c(s))
! print *, 'abmin(j,i)=', AuxData%abmin(i,j)
!     call print_asym('AFTER SCALE ABPLUS', AuxData%abplus, NDim1)
! call print_asym('AFTER SCALE ABMIN ', AuxData%abmin,  NDim1)
!         do i = 1, NDim1
!           do j = 1, NDim2
!                 r = indn1(1, i); s = indn1(2, i); rs = indx(i)
!                 p = indn2(1, j); q = indn2(2, j); pq = indx(j)

!                 nieb = AuxData%abplus(rs, pq) - AuxData%abplus(pq, rs)
! !                 if (abs(nieb).gt.1.d-5)then
! ! !                      if (abs(AuxData%abplus(rs, pq)).gt.1.d-5)then
! !                       write(*, '(A7, 2I3, 3F10.5)')'abpls', i, j, AuxData%abplus(rs, pq), AuxData%abplus(pq, rs), nieb
! !                 end if
!           end do
!     end do

    deallocate(wmat)
  end associate


end associate

! print*, 'pooo auxdata%cc', auxdata%cc
  
end subroutine pherpa_incore_spinres_ducc


subroutine pherpa_incore_spinres_ducc_fullmat(AuxData, h, Ints, Ndim1, &
        Ndim2, indn1, indn2, indx, multiply_by_S)
        !
        ! Builds the FULL matrix problem of Eq. 204:
        !   ( A  B ) (X_n)         ( -N  0 ) (X_n)
        !   ( B  A ) (Y_n) = omega ( 0   N ) (Y_n)
        !
        ! Output: AuxData%ABfull (2*NDim2 x 2*NDim2) and AuxData%Nmetric (2*NDim2 x 2*NDim2)
        ! No c-factor scaling is applied; N stays on the RHS.
        !
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints
        double precision, intent(in) :: h(:,:)
        integer, intent(in) :: NDim1, Ndim2
        integer, dimension(:,:), intent(in) ::indn1, indn2
        integer, dimension(:), intent(in) ::indx
        integer, intent(in) :: multiply_by_S


        integer :: i, j, t, u, v, w
        integer p, q, r, s, rs, pq
        integer :: pp, qq, rr, ss
        double precision :: Arspq
        double precision, dimension(:,:), allocatable :: auxi_a, auxio_a
        double precision, dimension(:,:), allocatable :: auxi_ab, auxio_ab
        double precision :: scale

        double precision :: T0, T1, T2, T3, T4, T5, T6, TW1, TW2
        double precision :: T1a, T1b, T2a, T2b
        double precision, dimension(:,:), allocatable:: wmat
        double precision :: num1, num2, num3
        double precision :: num_f, num_h
        double precision :: nnn, nieb
        double precision :: tp31, tp30, tp41, tp40, tp32, tpx0, tpplusz
        type (tclock) :: timer0, timer1, timerxx
        real(8) :: den, denmin
        integer :: imin, jmin
        integer :: NDimFull

        associate (n=>AuxData%Occ, n_p=>AuxData%n_p, n_m=>AuxData%n_m, NI=>AuxData%NI, &
              NA=>AuxData%NA, NV=>AuxData%NV, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, &
              IAux=>AuxData%IndAux, rdm2_pp=>AuxData%rdm2_pp, &
              rdm2_pm=>AuxData%rdm2_pm, rdm2_mm=>AuxData%rdm2_mm, &
              rdm2_mp=>AuxData%rdm2_mp, &
              rdm1_p=>AuxData%rdm1_p, rdm1_m=>AuxData%rdm1_m, &
              ints2e_aa=>Ints%Aints2e_aa, ints2e_ab=>Ints%Aints2e_ab, &
              Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)


          allocate(wmat(NBasis, NBasis))

          NDimFull = 2 * NDim2

          associate(c=>AuxData%cc)

          do i = 1, NBasis
                c(i) = sqrt(n(i))
                if (n(i) .lt. frac12) c(i) = -c(i)
          end do

          num1 = one
          wmat = zero

          do p = 1, NBasis
                do r = 1, NIA
                do t = 1, NIA
                   do v = 1, NIA
                         do u = 1, NIA

                               wmat(p, r) = wmat(p, r) + num1 * (&
                                     frac12 * erdm_ppx(rdm2_pp, n, r, v, t, u, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, v, t, p, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n, n, r, v, t, u, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(u, v, t, p, NBasis)) +&
                                   
                                     frac12 * erdm_ppx(rdm2_pp, n, r, v, u, t, IAux, NI) *&
                                     ints2e_aa(gmap_4fold(u, p, t, v, NBasis)) +&
                                     erdm_pmx(rdm2_pm, n, n, r, v, u, t, IAux, NI) *&
                                     ints2e_ab(gmap_4fold(v, t, p, u, NBasis)))

                       end do
                 end do
           end do

            end do
      end do


           AuxData%ABfull   = zero
           AuxData%Nmetric  = zero

          i_rowloops: do i = 1, NDim1
            j_colloops: do j = 1, NDim2
                r = indn1(1, i)
                s = indn1(2, i)
                rs = indx(i)

                pp = indn2(1, j)
                qq = indn2(2, j)
                pq = indx(j)

                do p = qq, pp, pp - qq
                      do q = qq, pp, pp - qq

                            if ( p .ne. q ) then

                            Arspq = zero

                            T0 = zero
                            if ( p == r ) T0 = T0 + (n(p) - n(s)) * Aints1e_aa(q, s)
                            if ( s == q ) T0 = T0 + (n(q) - n(r)) * Aints1e_aa(p, r)


                            T1a = zero; T1b = zero
                            T1a = T1a + func_T1_ducc(AuxData, Ints, r, p, s, q)
                            T1b = T1b + func_T2_ducc(AuxData, Ints, r, p, s, q)
                            T1 =  (T1a + T1b)

                            T2a = zero; T2b = zero
                            T2a = T2a + func_T1_ducc(AuxData, Ints, s, q, r, p)
                            T2b = T2b + func_T2_ducc(AuxData, Ints, s, q, r, p)
                            T2 = (T2a + T2b)

                            T5 = zero; T6 = zero
                            T6 = T6 - func_T34_ducc(AuxData, Ints, s, p, r, q)
                            T5 = T5 - func_T34_ducc(AuxData, Ints, q, r, p, s)


                            TW1 = zero; TW2 = zero
                            if (s == q) TW1 = TW1 - frac12 * wmat(p, r)
                            if (p == r) TW2 = TW2 - frac12 * wmat(q, s)

                            arspq = T0 + T1 + T2 + T6 + T5 + TW1 + TW2

                            ! r>s, p>q : same ordering -> contributes to A block
                            if ( r .gt. s .and. p .gt. q ) then
                                  ! A block: top-left and bottom-right
                                  AuxData%ABfull(rs,         pq        ) = &
                                  AuxData%ABfull(rs,         pq        ) + arspq
                                  AuxData%ABfull(rs+NDim2,   pq+NDim2  ) = &
                                  AuxData%ABfull(rs+NDim2,   pq+NDim2  ) + arspq
                            end if

                            ! r>s, q>p : cross ordering -> contributes to B block
                            if ( r .gt. s .and. q .gt. p ) then
                                  ! B block: top-right and bottom-left
                                  AuxData%ABfull(rs,         pq+NDim2  ) = &
                                  AuxData%ABfull(rs,         pq+NDim2  ) + arspq
                                  AuxData%ABfull(rs+NDim2,   pq        ) = &
                                  AuxData%ABfull(rs+NDim2,   pq        ) + arspq
                            end if

                      end if
                end do
          end do
    end do j_colloops
end do i_rowloops


    ! Build the metric N on the RHS: diag(-N, +N)
    do j = 1, NDim2
          p = indn2(1, j); q = indn2(2, j); pq = indx(j)
          AuxData%Nmetric(pq,        pq       ) = -(n(p) - n(q))   ! -N block (X part)
          AuxData%Nmetric(pq+NDim2,  pq+NDim2 ) =  (n(p) - n(q))   ! +N block (Y part)
    end do

    

    deallocate(wmat)
  end associate


end associate

end subroutine pherpa_incore_spinres_ducc_fullmat

subroutine print_asym(label, A, ndim)
  character(len=*), intent(in) :: label
  integer, intent(in) :: ndim
  real(8), intent(in) :: A(ndim,ndim)

  real(8), allocatable :: diff(:,:)
  integer :: loc(2)

  allocate(diff(ndim,ndim))

  diff = A - transpose(A)
  loc = maxloc(abs(diff))

  print *, trim(label)
  print *, 'norm asym = ', sqrt(sum(diff*diff))
  print *, 'norm A    = ', sqrt(sum(A*A))
  print *, 'rel asym  = ', sqrt(sum(diff*diff)) / max(sqrt(sum(A*A)), 1.0d-30)
  print *, 'max asym  = ', maxval(abs(diff))
  print *, 'loc       = ', loc(1), loc(2)
  print *, 'values    = ', A(loc(1),loc(2)), A(loc(2),loc(1))

  deallocate(diff)
end subroutine


  subroutine pherpa_symm_spinres(AuxData)
        !
        ! A SYMMETRIZED PROBLEM A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
        ! ABPLUS IS CHANGED AND TURNS INTO ABPLUS^(1/2)
        !
        type(TACppData), intent(inout) :: AuxData
        double precision, dimension(:,:), allocatable :: work
        real(F64), allocatable, dimension(:, :) :: work_a, work_b, work_c
        real(F64), allocatable, dimension(:)    :: vec_tmp
        real(F64), allocatable, dimension(:)    :: lambda, lambdam
        real(F64)                               :: omega, omega_sq, dot_prod, norm_fact
        integer                                 :: no_neg, no_negm
        integer                                 :: j, k, mode
        double precision, parameter             :: small = 1.d-6
        double precision, allocatable, dimension(:, :) :: diff
        integer, dimension(2) :: loc 
        associate(ndim=>AuxData%NDim)
          
!          allocate(AuxData%Eigs(ndim))
!          allocate(AuxData%Eigvec(ndim, ndim))        
        
          
          allocate(work_a(ndim, ndim))
          allocate(work_b(ndim, ndim))
          allocate(work_c(ndim, ndim))
          allocate(vec_tmp(ndim))
          allocate(lambda(ndim))
          allocate(lambdam(ndim))
          allocate(diff(ndim, ndim))
          
          associate ( abplus => AuxData%ABplus, abmin => AuxData%ABmin, &
                eigvec => AuxData%Eigvec, eigs => AuxData%Eigs )


            ! A = 0.5 * (A + A^T)
            abplus = frac12 * (abplus + transpose(abplus))
            abmin = frac12 * (abmin + transpose(abmin))
            
            !
            ! A+^(1/2) =  U * sqrt(Lambda) * U^T = work_a * sqrt(lambda) * work_a^T
            ! Diagonalize for to get work_a and lambda
            !
            work_a = abplus
            work_c = abmin
            call symmetric_eigenproblem(lambda, work_a, ndim, .true.)
            call symmetric_eigenproblem(lambdam, work_c, ndim, .true.)

            no_neg = 0
            no_negm = 0
            do k = 1, ndim
                  if (lambda(k) <= small) no_neg = no_neg + 1
                  if (lambdam(k) <= small) no_negm = no_negm + 1
            end do
            AuxData%Apnegs = no_neg
            AuxData%Amnegs = no_negm
            
!            if (no_neg /= 0) then
!                  print *, "A+ is not positive definite, count lambda <= small:", no_neg
!                  print *, "min eigenvalue A+ =", minval(lambda)
!            end if

 
          !
          ! Step 2: Calculate H = U * sqrt(lambda) = work_b = work_a * sqrt(Lambda)
          !
          do j = 1, ndim
                omega = sqrt(abs(lambda(j)))
                work_b(:, j) = work_a(:, j) * omega
          end do

          !
          ! Step 3: Calculate A+^(1/2) = H * U^T = work_b * work_a^T -> write to abplus
          !
          call real_abT(abplus, work_b, work_a)

          !
          ! M = A+^(1/2) * A- * A+^(1/2) !! A+^(1/2) is inside abplus now
          !
          work_a = zero
          work_b = zero
          call real_abaT(work_a, abplus, abmin, work_b)


          call symmetric_eigenproblem(eigs, work_a, ndim, .true.)


!           print '(A,F8.4,A)', 'alpha=', AuxData%alpha, 'skipped eigs of A+A-:'
!           print '(6E14.6)', (eigs(j), j=1,ndim)
!           print '(A,F8.4,A,E14.6)', 'alpha=', AuxData%alpha, ' min(lambda A+)=', minval(lambda)
          !           print '(A,F8.4,A)', 'alpha=', AuxData%alpha, 'skipped eigs of A+A-:'

          ! Y = A+^(1/2) * Z
          call real_ab(eigvec, abplus, work_a)
          !print *, "Y vectors computed"

          ! 2 * omega^(-1) * Y^T * A- * Y = 1

          do j = 1, 10
                write(*, '(I3, 3F15.8)')j, eigs(j), lambda(j), lambdam(j)
          end do

          AuxData%AAnegs = 0
          do mode = 1, ndim
                omega = eigs(mode)
                
                if (omega > small) then
                      omega_sq = sqrt(omega)
                      eigs(mode) = omega_sq

                      ! tmp = A- * Y; dot = Y^T * tmp
                      call real_Av(vec_tmp, abmin, eigvec(:, mode))
                      call real_vw_x(dot_prod, eigvec(:, mode), vec_tmp, ndim)

                      ! sum_nu = (2 / omega) * (Y^T A- Y)
                      norm_fact = (two / omega_sq) * dot_prod

                      if (norm_fact > zero) then                            
                            norm_fact = one / sqrt(norm_fact)
                      else
                            print*, 'ALERT!!!! norm is zero', norm_fact, omega, omega_sq
                            norm_fact = zero
                      end if

                      call real_scal(eigvec(:, mode), norm_fact)
                else
                      !write(*, '(F20.10)') omega
                      AuxData%AAnegs = AuxData%AAnegs + 1
                end if
          end do

        end associate

        deallocate(work_a)
        deallocate(work_b)
        deallocate(vec_tmp)
        deallocate(lambda)
 
        
      end associate

end subroutine pherpa_symm_spinres

subroutine pherpa_nonsymm_spinres(AuxData)
      ! Solves the FULL non-symmetric problem of Eq. 204 directly via DGGEV.
      ! Produces EigvX, EigvY, Eigs in the same convention as restore_XY output.
      type(TACppData), intent(inout) :: AuxData
      integer :: NF, info, lwork, k, i, p, q, kept
      integer :: kk
      double precision, allocatable :: M(:,:), Nmet(:,:)
      double precision, allocatable :: alphar(:), alphai(:), beta(:)
      double precision, allocatable :: vl(:,:), vr(:,:), work(:)
      double precision, allocatable :: omega_all(:), Xtmp(:), Ytmp(:)
      double precision :: om, mnorm, scale_fact
      integer, allocatable :: order(:)
      double precision, parameter :: small = 1.d-6
      double precision, parameter :: small_imag = 1.d-8

      associate(ndim=>AuxData%NDim, ABfull=>AuxData%ABfull, Nmetric=>AuxData%Nmetric, &
                EigvX=>AuxData%EigvX, EigvY=>AuxData%EigvY, Eigs=>AuxData%Eigs, &
                IndN=>AuxData%IndN, n=>AuxData%Occ, alpha=>AuxData%alpha)

        NF = 2 * ndim
        allocate(M(NF, NF), Nmet(NF, NF))
        allocate(alphar(NF), alphai(NF), beta(NF))
        allocate(vl(1,1), vr(NF, NF))

        ! DGGEV destroys input matrices; copy first.
        M    = ABfull
        Nmet = Nmetric
        
        ! Explicitly symmetrize the full matrix M to eliminate numerical asymmetric noise.
        ! This matches the symmetrization done on ABplus and ABmin in the symmetric path,
        ! and prevents degenerate roots from splitting into complex conjugate pairs in dggev.
        M = frac12 * (M + transpose(M))

        ! Workspace query
        lwork = -1
        allocate(work(1))
        call dggev('N', 'V', NF, M, NF, Nmet, NF, alphar, alphai, beta, &
                   vl, 1, vr, NF, work, lwork, info)
        lwork = int(work(1))
        deallocate(work); allocate(work(lwork))

        call dggev('N', 'V', NF, M, NF, Nmet, NF, alphar, alphai, beta, &
                   vl, 1, vr, NF, work, lwork, info)
        if (info /= 0) print *, 'DGGEV info =', info

        ! Diagnostic: report any complex roots
        ! do k = 1, NF
        !       if (abs(beta(k)) > small .and. abs(alphai(k)/beta(k)) > small_imag) then
        !             print '(A,F8.4,A,I5,A,2E14.6)', 'alpha=', alpha, &
        !                   ' COMPLEX eigenvalue at k=', k, ' val=', &
        !                   alphar(k)/beta(k), alphai(k)/beta(k)
        !       end if
        ! end do

        ! Select n+ set: real eigenvalues, take positive omega and the matching eigenvector.
        ! For each kept root, normalize so that -X^T N X + Y^T N Y = 1/2
        EigvX = zero
        EigvY = zero
        Eigs  = zero
        kept = 0

        do k = 1, NF
      if (abs(beta(k)) < small) cycle
      if (abs(alphai(k)) > small_imag * abs(beta(k))) cycle

      om = alphar(k) / beta(k)

      ! Compute metric norm using the eigenvector AS-IS
      mnorm = zero
      do i = 1, ndim
            p = IndN(1, i); q = IndN(2, i)
            mnorm = mnorm + (n(p) - n(q)) * &
                  ( vr(i+ndim, k)**2 - vr(i, k)**2 )
      end do

      ! n+ selector: keep ONLY roots with mnorm > 0 (physical normalization)
      ! Drop roots with mnorm < 0 — they are the n- partners of the n+ set.
      ! Drop roots with |mnorm| ~ 0 — singular metric (n_p == n_q).
      if (mnorm < small) cycle    ! covers both negative and ~zero cases

      ! If physical ω_{n+} is negative, take it as-is. The X_n, Y_n
      ! eigenvector is already correctly assigned by DGGEV.
      ! Do NOT take abs(om).

      kept = kept + 1
      if (kept > ndim) then
            print *, 'WARNING: more n+ roots than ndim'
            exit
      end if

      scale_fact = sqrt(frac12 / mnorm)
      Eigs(kept) = om                        ! KEEP the sign of ω
      do i = 1, ndim
            EigvX(i, kept) = scale_fact * vr(i,        k)
            EigvY(i, kept) = scale_fact * vr(i + ndim, k)
      end do
end do

        ! do k = 1, NF
        !       if (abs(beta(k)) < small) cycle                  ! infinite eigenvalue
        !       if (abs(alphai(k)) > small_imag * abs(beta(k))) cycle  ! complex root, skip

        !       om = alphar(k) / beta(k)

        !       ! we want the n+ set; if our convention puts n+ at omega > 0, take those
        !       if (om <= small) cycle

        !       ! Compute metric norm: -X^T N X + Y^T N Y, where N_{pq,pq} = n_p - n_q
        !       mnorm = zero
        !       do i = 1, ndim
        !             p = IndN(1, i); q = IndN(2, i)
        !             mnorm = mnorm + (n(p) - n(q)) * &
        !                   ( vr(i+ndim, k)**2 - vr(i, k)**2 )   ! Y^2 - X^2 weighted by N
        !       end do

        !       if (abs(mnorm) < small) then
        !             print *, 'ALERT: zero metric norm at k=', k, ' omega=', om
        !             cycle
        !       end if

        !       if (mnorm < zero) then
        !             ! sign convention says omega should have been negative; flip
        !             om = -om
        !             mnorm = -mnorm
        !             ! omega is still positive after the flip in magnitude; we keep |om|
        !             om = abs(om)
        !       end if

        !       scale_fact = sqrt(frac12 / mnorm)

        !       kept = kept + 1
        !       if (kept > ndim) then
        !             print *, 'WARNING: more n+ roots than ndim, truncating'
        !             exit
        !       end if

        !       Eigs(kept) = om
        !       do i = 1, ndim
        !             EigvX(i, kept) = scale_fact * vr(i,        k)
        !             EigvY(i, kept) = scale_fact * vr(i + ndim, k)
        !       end do
        ! end do

!        if (kept < ndim) then
!              print '(A,I5,A,I5)', 'pherpa_nonsymm: kept ', kept, ' of expected ', ndim
!        end if

        ! Optionally: sort by ascending omega so order matches the symm version
        ! (skipped here; add if pherpa_energy_ducc cares about ordering)

        deallocate(M, Nmet, alphar, alphai, beta, vl, vr, work)
      end associate
end subroutine pherpa_nonsymm_spinres

  subroutine pherpa_energy(AuxData, ecorr, Ints)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints        
        integer :: p, q, r, s, i, j, ii, k
        integer :: skip
        double precision, intent(out) :: ecorr
        double precision, parameter :: small_e = 1.d-2
        double precision, parameter :: big_e = 1.d+8
        double precision :: temp, aux1, aux2, plusz
        double precision, dimension(:), allocatable :: skipped
        double precision, dimension(:,:), allocatable :: work, temp_mat
        double precision, dimension(:,:), allocatable :: work_xx, work_xy, work_yx, work_yy
        double precision :: aux_xx, aux_yy, aux_xy, aux_yx, ww
        double precision :: W_aa, W_ab, W_ref, W_aa_plusz, W_aa_plusz2, przycz
        double precision, dimension(:,:), allocatable :: work_x, work_y

        associate(NDim=>AuxData%NDim, c=>AuxData%cc, n=>AuxData%occ, indN=>AuxData%IndN, &
              IAux=>AuxData%IndAux, eigs=>AuxData%Eigs, eigvecs=>AuxData%Eigvec, &
              eigvX=>AuxData%EigvX, eigvY=>AuxData%EigvY, rdm2_pp=>AuxData%rdm2_pp, NI=>AuxData%NI)

          skip = 0
          ecorr = zero
          allocate(skipped(NDim))
          allocate(work(NDim, NDim))

          allocate(temp_mat(NDim, NDim))
          work = eigvecs

          do k = 1, ndim
                if (.not.(eigs(k) > small_e .and. eigs(k) < big_e)) then
                      print*, 'this eigval is skipped', eigs(k)
                      skip = skip + 1
                      skipped(skip) = eigs(k)
                      work(:, k) = zero
                end if
          end do

          call real_abT(temp_mat, work, work)

          do i = 1, ndim
                p = indn(1, i)
                q = indn(2, i)

                do j = 1, ndim
                      r = indn(1, j)
                      s = indn(2, j)

                      if ( .not. ( IAux(p) == IAux(q) .and. &
                            IAux(r) == IAux(s) .and. &
                            IAux(p) == IAux(r) ) )then


                            temp = temp_mat(i, j)
                            
                            aux1 = temp_mat(i, j)* (c(s) + c(r)) * (c(p) + c(q)) 
                            
                            ! aux_xx = work_xx(i,j) * (n(p)-n(q)) * (n(r) - n(s))


                            aux2 = zero
                            if (q == s .and. p == r) then
                                  aux2 = aux2 - n(p) * (one - n(s)) - n(s) * (one - n(p))
                            end if
                            
                            ecorr = ecorr + aux1 * (Ints%ints2e_aa(gmap_4fold(p,q,r,s,AuxData%NBasis))) +&
                                  aux1 * Ints%ints2e_ab(gmap_4fold(p,q,r,s,AuxData%NBasis)) + &
                                  aux2 * (Ints%ints2e_aa(gmap_4fold(p,q,r,s,AuxData%NBasis)))
                      end if
                end do
          end do

          deallocate(skipped)
          deallocate(work)
          deallocate(temp_mat)


        end associate
        
  end subroutine pherpa_energy

  subroutine pherpa_energy_ducc(AuxData, ecorr, Ints, IntsD)
        type(TACppData), intent(inout) :: AuxData
        type(TInts), intent(inout) :: Ints, IntsD
        integer :: p, q, r, s, i, j, ii, k
        integer :: skip
        double precision, intent(out) :: ecorr
        double precision, parameter :: small_e = 1.d-3
        double precision, parameter :: big_e = 1.d+8
        double precision :: temp, aux1, aux2, plusz, term_d, term_e
        double precision :: aux_xx, aux_yy, aux_xy, aux_yx
        double precision :: W_aa, W_ab, W_aa_ref, W_ab_ref
        double precision :: W_aa_act_zw, W_aa_act_dc
        double precision :: W_aa_act_zw_rs, W_aa_act_dc_rs
        double precision :: W_aa_act, W_ab_act, W_aa_ref_act, W_ab_ref_act
        double precision :: W_aa_xx, W_aa_xy, W_ab_xy, W_ab_xx, U_aa_ref_act

        double precision :: T_aa_ref_act, T_ab_ref_act, T_aa_ref, T_ab_ref
        
        double precision :: U_aa_ref_act_zw, U_aa_ref_act_dc
        double precision :: W_aa_ref_act_zw, W_aa_ref_act_dc
        double precision, dimension(:), allocatable :: skipped
        integer, dimension(:), allocatable :: skipped_int
        double precision, dimension(:,:), allocatable :: work_xx, work_xy, work_yx, work_yy
        double precision, dimension(:,:), allocatable :: work_x, work_y
        double precision :: aints_aa_rs, aints_aa_sr, aints_ab_rs, aints_ab_sr
        double precision :: ecorr2

        integer :: ip, iq, ir, is, pq, rs

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, NA=>AuxData%NA,&
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              indn=>AuxData%IndN, &
              c=>AuxData%cc, n=>AuxData%Occ,&
              ecorr_calc=>ecorr, &
              NDim=>AuxData%NDim, &
              eigs=>AuxData%Eigs, &
!              EigvecY=>AuxData%Eigvec, &
!              EigvecX=>AuxData%EigvecX)
              EigvecY=>AuxData%EigvY, &
              EigvecX=>AuxData%EigvX)


          skip = 0
          ecorr = zero
          allocate(skipped(NDim))
          allocate(skipped_int(NDim))
          allocate(work_x(NDim, NDim))
          allocate(work_y(NDim, NDim))
          allocate(work_xx(NDim, NDim))
          allocate(work_yy(NDim, NDim))
          allocate(work_xy(NDim, NDim))
          allocate(work_yx(NDim, NDim))
          work_y = EigvecY
          work_x = EigvecX
          skipped = zero
          skipped_int = 0
!          print*, 'eigvals'
!           do k = 1, ndim
!                 if (.not.(eigs(k) > small_e .and. eigs(k) < big_e)) then
! !                      print*, 'this is skipped?', eigs(k)
! !                      print*, 'yes', k
!                       skip = skip + 1
!                       skipped(skip) = eigs(k)
!                       skipped_int(skip) = k
!                       work_x(:, k) = zero
!                       work_y(:, k) = zero
!                 !else
!                 !      print*, k, eigs(k)
!                 end if
!           end do


do k = 1, ndim
      if (.not.(eigs(k) > small_e .and. eigs(k) < big_e)) then
            skip = skip + 1
!            print '(A,F8.4,A,I4,A,E14.6)', &
!                 'alpha=', AuxData%alpha, ' SKIPPED mode ', k, ' eigs=', eigs(k)
            skipped(skip) = eigs(k)
            skipped_int(skip) = k
            work_x(:, k) = zero
            work_y(:, k) = zero
      end if
end do

          
!          do k = 1, skip
!                print*, 'skipped', skipped_int(k), skipped(k)
!          end do
!          print*, ''
          ! do k = 1, ndim
          !       print*, work_x(k, 172)
          ! end do
          ! print*, ''
!           call real_abT(work_xx, work_x, work_x)
!           call real_abT(work_xy, work_x, work_y)
!           call real_abT(work_yx, work_y, work_x)
!           call real_abT(work_yy, work_y, work_y)

         !  do k = 1, ndim
         !        print*, work_xx(k, 172)
         !  end do
         ! print*, ''
         
          work_xx = zero
          work_xy = zero
          work_yx = zero
          work_yy = zero
          do pq = 1, ndim
                do rs = 1, ndim
                      do k = 1, ndim
                            work_xx(pq, rs) =  work_xx(pq, rs) + work_x(pq, k)*work_x(rs,k)
                            ! if (rs==1)then
                            !       if (abs(work_xx(pq, rs)).gt.1.d+2)then
                            !             write(*, '(A5, 3I5, 3F15.7)')'wss', pq, rs, k,   work_x(pq, k), work_x(rs,k), work_xx(pq, rs) 
                            !       end if
                            !end if
                                  work_xy(pq, rs) =  work_xy(pq, rs) + work_x(pq, k)*work_y(rs,k)
                                  work_yx(pq, rs) =  work_yx(pq, rs) + work_y(pq, k)*work_x(rs,k)
                                  work_yy(pq, rs) =  work_yy(pq, rs) + work_y(pq, k)*work_y(rs,k)

                      end do
                end do
          end do
!          print*, 'pluszona'
          W_aa = zero
          W_ab = zero
          W_aa_ref = zero
          W_ab_ref = zero

          W_aa_act = zero
          W_aa_act_zw = zero
          W_aa_act_dc = zero

          W_aa_act_zw_rs = zero
          W_aa_act_dc_rs = zero

          W_ab_act = zero
          W_aa_ref_act = zero



          U_aa_ref_act_zw = zero
          U_aa_ref_act_dc = zero
          
          W_ab_ref_act = zero
          W_aa_xx = zero
          W_aa_xy = zero
          W_ab_xx = zero
          W_ab_xy = zero
          U_aa_ref_act = zero


          T_aa_ref_act = zero
          T_ab_ref_act = zero
          T_aa_ref = zero
          T_ab_ref = zero

          ! print*, Ints%ints2e_aa(gmap_4fold(2,1,2,1,NBasis))
          ! print*, 'W_aact  p  q  r  s   aints_rs    aints_sr  aux_xy    aux_yx   aux_xx    aux_yy     W_aa_act'
          do i = 1, ndim
                p = indn(1, i)
                q = indn(2, i)

                do j = 1, ndim
                      r = indn(1, j)
                      s = indn(2, j)
                      aux_xx = work_xx(i,j) * (n(p)-n(q)) * (n(r) - n(s))
                      aux_xy = work_xy(i,j) * (n(p)-n(q)) * (n(r) - n(s))
                      aux_yx = work_yx(i,j) * (n(p)-n(q)) * (n(r) - n(s))
                      aux_yy = work_yy(i,j) * (n(p)-n(q)) * (n(r) - n(s))


                      
                      if (IAux(p) == IAux(r) .and. IAux(q) == IAux(s) .and. &
                            IAux(p) == IAux(q).and.IAux(p)==1) then


                            aints_aa_rs = Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) - &
                                  IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,r-NI,s-NI,NA))

                            W_aa_act_zw_rs = W_aa_act_zw_rs - frac12 * Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) *  (aux_xy + aux_yx)
                            W_aa_act_dc_rs = W_aa_act_dc_rs + frac12 * IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,r-NI,s-NI,NA)) *  (aux_xy + aux_yx)
                            
                            W_aa_act = W_aa_act - frac12 * aints_aa_rs *  (aux_xy + aux_yx)

                            !AuxData%W_aa0(i,j) = AuxData%W_aa0(i,j) - frac12 * aints_aa_rs *  (aux_xy + aux_yx)
                            
                            U_aa_ref_act = U_aa_ref_act - aints_aa_rs * rdm2_pp(p-NI, r-NI, q-NI, s-NI)

                            !T_aa_ref_act = T_aa_ref_act - frac12 *aints_aa_rs * AuxData%ph_rdm2_aa(p, r, q, s)
                            

                            
                            U_aa_ref_act_zw = U_aa_ref_act_zw - frac12 * Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) *  &
                                  rdm2_pp(p-NI, r-NI, q-NI, s-NI)
                            U_aa_ref_act_dc = U_aa_ref_act_dc + frac12 * IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,r-NI,s-NI,NA)) * &
                                  rdm2_pp(p-NI, r-NI, q-NI, s-NI)

                            
                            aints_aa_sr = Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) - &                                  
                                  IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,s-NI,r-NI,NA))
                            W_aa_act = W_aa_act + frac12 * aints_aa_sr * (aux_xx + aux_yy)
                            
                            !AuxData%W_aa0(i,j) = AuxData%W_aa0(i,j) + frac12 * aints_aa_sr * (aux_xx + aux_yy)

                            W_aa_act_zw = W_aa_act_zw + frac12 * Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) *  (aux_xx + aux_yy)
                            W_aa_act_dc = W_aa_act_dc - frac12 * IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,s-NI,r-NI,NA)) *  (aux_xx + aux_yy)
                            
                            
                            ! if (abs(work_xx(i,j)).gt.1.d+1)then
                            !       write(*, '(A5, 4I5, 7F15.7)')'auxx', r, s, p, q, work_xx(i,j) , (n(p)-n(q)) * (n(r) - n(s)), aux_xx, W_aa_act_zw, W_aa_act_dc, eigs(i), eigs(j)
                            ! end if
                            ! if (abs(work_xy(i,j)).gt.1.d+1)then
                            !       write(*, '(A5, 4I5, 7F15.7)')'auxy', r, s, p, q, work_xy(i,j) , (n(p)-n(q)) * (n(r) - n(s)), aux_xy, W_aa_act_zw_rs, W_aa_act_dc_rs, eigs(i), eigs(j)
                            ! end if
                            ! if (abs(work_yx(i,j)).gt.1.d+1)then
                            !       write(*, '(A5, 4I5, 7F15.7)')'auyx', r, s, p, q, work_yx(i,j) , (n(p)-n(q)) * (n(r) - n(s)), aux_yx, W_aa_act_zw_rs, W_aa_act_dc_rs, eigs(i), eigs(j)
                            ! end if
                            ! if (abs(work_yy(i,j)).gt.1.d+1)then
                            !       write(*, '(A5, 4I5, 7F15.7)')'auyy', r, s, p, q, work_yy(i,j) , (n(p)-n(q)) * (n(r) - n(s)), aux_yy, W_aa_act_zw, W_aa_act_dc, eigs(i), eigs(j)
                            ! end if



                            U_aa_ref_act = U_aa_ref_act -  aints_aa_sr * rdm2_pp(p-NI, s-NI, q-NI, r-NI)

!                            T_aa_ref_act = T_aa_ref_act - frac12 * aints_aa_sr * AuxData%ph_rdm2_aa(p, s, q, r)
                            
                            U_aa_ref_act_zw = U_aa_ref_act_zw - frac12 * Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) *  &
                                  rdm2_pp(p-NI, s-NI, q-NI, r-NI)
                            U_aa_ref_act_dc = U_aa_ref_act_dc + frac12 * IntsD%ints2e_aa(gmap_4fold(p-NI,q-NI,s-NI,r-NI,NA)) * &
                                  rdm2_pp(p-NI, s-NI, q-NI, r-NI)



                            aints_ab_rs = Ints%ints2e_ab(gmap_4fold(p,q,r,s,NBasis)) - &
                                  IntsD%ints2e_ab(gmap_4fold(p-NI,q-NI,r-NI,s-NI,NA))

                            W_ab_act = W_ab_act - aints_ab_rs * (aux_xy + aux_yx)
!                            AuxData%W_ab0(i,j) = AuxData%W_ab0(i,j) - aints_ab_rs * (aux_xy + aux_yx)

!                            T_ab_ref_act = T_ab_ref_act - aints_ab_rs  * AuxData%ph_rdm2_ab(p, r, q, s)

                            aints_ab_sr = Ints%ints2e_ab(gmap_4fold(p,q,s,r,NBasis)) - &
                                  IntsD%ints2e_ab(gmap_4fold(p-NI,q-NI,s-NI,r-NI,NA))

                            W_ab_act = W_ab_act + aints_ab_sr * (aux_xx + aux_yy)
!                            AuxData%W_ab0(i,j) = AuxData%W_ab0(i,j) + aints_ab_sr * (aux_xy + aux_yx)

!                            T_ab_ref_act = T_ab_ref_act + aints_ab_sr  * AuxData%ph_rdm2_ab(p, s, q, r)

                      else
                            if (.not.(IAux(p) == IAux(r) .and. &
                                  IAux(q) == IAux(s) .and. &
                                  IAux(p) == IAux(q)  .and. &
                                  IAux(p)==2)) then                                  

                            aux2 = zero
                            if (q == s .and. p == r) then
                                  aux2 = aux2 + n(p) * (one - n(q)) + n(q) * (one - n(p))
                            end if
                            
                            W_aa = W_aa - &
                                  frac12 * Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) * (aux_xy + aux_yx)
                            W_aa_xy = W_aa_xy - &
                                  frac12 * Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) * (aux_xy + aux_yx)

                            W_aa = W_aa + &
                                  frac12 * Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) * (aux_xx + aux_yy)

                            W_aa_xx = W_aa_xx + &
                                  frac12 * Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) * (aux_xx + aux_yy) 
                                                        

                            W_ab = W_ab - Ints%ints2e_ab(gmap_4fold(p,q,r,s,NBasis)) * (aux_xy + aux_yx)
                            W_ab_xy = W_ab_xy - Ints%ints2e_ab(gmap_4fold(p,q,r,s,NBasis)) * (aux_xy + aux_yx)
                            

                            W_ab = W_ab + Ints%ints2e_ab(gmap_4fold(p,q,s,r,NBasis)) * (aux_xx + aux_yy)
                            W_ab_xx = W_ab_xx + Ints%ints2e_ab(gmap_4fold(p,q,s,r,NBasis)) * (aux_xx + aux_yy)


                            W_aa_ref = W_aa_ref - &                                  
                                  frac12 * Ints%ints2e_aa(gmap_4fold(p,q,s,r,NBasis)) * aux2 
                            
                      end if
                end if
          end do
    end do

    
    do ip = 1, NA
        p = ip + NI
        do iq = 1, NA
            q = iq + NI
            do ir = 1, NA
                r = ir + NI
                do is = 1, NA
                    s = is + NI

                    if (p.ne.q.and.r.ne.s)then
                          W_aa_ref_act = W_aa_ref_act - frac12 * ( &
                                Ints%ints2e_aa(gmap_4fold(p,q,r,s,NBasis)) - &
                                IntsD%ints2e_aa(gmap_4fold(ip,iq,ir,is,NA)) ) * &
                                rdm2_pp(ip, ir, iq, is)
                    end if
                    
                    W_ab_ref_act = W_ab_ref_act - ( &
                        Ints%ints2e_ab(gmap_4fold(p,q,r,s,NBasis)) - &
                        IntsD%ints2e_ab(gmap_4fold(ip,iq,ir,is,NA)) ) * &
                        rdm2_pm(ip, ir, iq, is)

                    
                end do
            end do
        end do
    end do
!     print*, 'W_aa_ref_act', W_aa_ref_act
!     print*, 'W_ab_ref_act', W_ab_ref_act
!     print*, ''
!     print*, 'W_aa_act_zw', W_aa_act_zw
!     print*, 'W_aa_act_dc', W_aa_act_dc
!     print*, 'W_aa_act', W_aa_act
!     print*, ''
!     print*, 'W_ab_act', W_ab_act
!     print*, ''
!     print*, 'U_aa_ref_act', U_aa_ref_act
!     print*, 'T_aa_ref_act', T_aa_ref_act
!     print*, 'T_ab_ref_act', T_ab_ref_act
! !    U_aa_ref_act=-U_aa_ref_act
!     print*, 'U_aa_ref_act_zw', U_aa_ref_act_zw
!     print*, 'U_aa_ref_act_dc', U_aa_ref_act_dc
!     print*, 'U_aa_ref_act_dc + U_aa_ref_act_zw', U_aa_ref_act_dc + U_aa_ref_act_zw
!     print*, ''
!     print*, ''
!     print*, 'W_aa_ref', W_aa_ref
!     print*, 'W_aa_xx', W_aa_xx
!     print*, 'W_aa_xy', W_aa_xy
!     print*, 'W_ab_xy', W_ab_xy
!     print*, 'W_ab_xx', W_ab_xx

    if (abs(AuxData%alpha).lt.1d-5)then

          AuxData%W_aa_act0 = W_aa_act
          AuxData%W_ab_act0 = W_ab_act
          AuxData%W_aa0 = W_aa + W_ab

    end if
    ecorr2 = W_aa_act + U_aa_ref_act +  W_ab_act +  W_ab_ref_act +  W_aa +  W_aa_ref +  W_ab
    ecorr = W_aa_act - AuxData%W_aa_act0 +  W_ab_act -  AuxData%W_ab_act0 +  W_aa +   W_ab - AuxData%W_aa0 - AuxData%W_ab0
    
!    write(*,'(A5, 9A13)') ' ', 'W_aa_act', 'U_aa_ref_act', 'W_ab_act', 'W_ab_ref_act', 'W_aa', 'W_aa_ref', 'W_ab', 'ecorr', ' '
!    write(*, '(A5, 9F13.6)') 'www', W_aa_act,   U_aa_ref_act,        W_ab_act,      W_ab_ref_act,       W_aa,    W_aa_ref,      W_ab,          ecorr2
!    write(*, '(14F13.6)') AuxData%alpha, W_aa_act, - AuxData%W_aa_act0 ,  W_ab_act,  -  AuxData%W_ab_act0 ,  W_aa ,  - AuxData%W_aa0, W_ab,   ecorr, W_aa_act_zw_rs, W_aa_act_dc_rs, W_aa_act_zw, W_aa_act_dc
!    write(*, '(14F13.6)') AuxData%alpha, ecorr
!     print*, ''
!stop          
          ! if (skip .ne. 0) then
          !       write(6, '(/,1x,"the number of total eigenvalues is",i4)') ndim             
          !       write(6, '(/,1x,"the number of discarded eigenvalues is",i4)') skip
          !       do ii = 1, skip
          !             write(6, *) "skipped", ii, skipped(ii)
          !       end do
          ! end if

          deallocate(skipped)
          deallocate(work_xx)
          deallocate(work_xy)
          deallocate(work_yy)
          deallocate(work_yx)



        end associate
        
       end subroutine pherpa_energy_ducc

       subroutine pherpa_gen_rdm(AuxData, Eee, Ints)

             type(TACppData), intent(inout) :: AuxData
             type(TInts), intent(inout) :: Ints
             double precision, intent(out) :: Eee
             
             integer :: p, q, r, s
             integer :: i, nu
             integer :: NBasis, NDim
             integer :: skip
             
             double precision, parameter :: small_e = 1.d-3
             double precision, parameter :: big_e   = 1.d8

             double precision :: TnuTnu, rdm2h,rdm2x, EeeH, EeeX 
             double precision :: gamma_aa, gamma_ab
             double precision :: int_aa, int_ab, rrr, Err

             double precision, dimension(:,:,:), allocatable :: Tnu
             double precision, dimension(:), allocatable :: skipped

             associate( &
                   n_p     => AuxData%n_p,      &
                   n_m     => AuxData%n_m,      &
                   n       => AuxData%Occ,      &
                   indn    => AuxData%IndN,     &
                   eigs    => AuxData%Eigs,     &
                   EigvecX => AuxData%EigvX,    &
                   EigvecY => AuxData%EigvY )
               
               NBasis = AuxData%NBasis
               NDim   = AuxData%NDim


               skip  = 0
               
               allocate(Tnu(NBasis, NBasis, NDim))
               allocate(skipped(NDim))

               Tnu = zero
               skipped = zero

               ! ------------------------------------------------------------
               ! Build full transition-density matrices TNU(p,q,nu)
               !
               ! For pair p > q stored as i:
               !
               ! T_qp^nu = (n(p)-n(q)) Y_pq^nu
               ! T_pq^nu = (n(q)-n(p)) X_pq^nu
               !
               ! T_pp^nu = 0 by construction
               ! ------------------------------------------------------------
               
               do nu = 1, NDim
                     if (eigs(nu) > small_e .and. eigs(nu) < big_e) then
                           do i = 1, NDim
                                 p = indn(1,i)
                                 q = indn(2,i)
                                 
                                 TNU(q,p,nu) = (n(p) - n(q)) * EigvecY(i,nu)
                                 TNU(p,q,nu) = (n(q) - n(p)) * EigvecX(i,nu)

                           end do
                     else
                           skip = skip + 1
                           skipped(skip) = eigs(nu)
                     end if
               end do

               AuxData%ph_rdm2_aa = zero
               AuxData%ph_rdm2_ab = zero
               
               ! ------------------------------------------------------------
               ! Full ERPA-reconstructed 2-RDM
               !
               ! Gamma_pqrs = n_p n_q delta_pr delta_qs + sum_nu T_pr^nu T_sq^nu
               ! - n_q delta_qr delta_ps
               !
               ! Exchange delta term only for same-spin blocks.
               ! ------------------------------------------------------------
               Eee = zero
               Err = zero
               Eeeh  = zero
               Eeex =  zero
               do p = 1, NBasis
                     do q = 1, NBasis
                           do r = 1, NBasis
                                 do s = 1, NBasis
                                       
                                       gamma_aa = zero
                                       gamma_ab = zero

                                       rdm2h = zero
                                       rdm2x = zero
                                       ! Hartree-like part
                                       if (p == r .and. q == s) then
                                             gamma_aa = gamma_aa + n(p) * n(q)
                                             gamma_ab = gamma_ab + n(p) * n(q)
                                             rdm2h = rdm2h +   n(p) * n(q)

                                       end if

                                       if (q == r .and. p == s) then
                                             rdm2x = rdm2x -  n(p) * n(q)                                             
                                       end if

                                       ! Exchange-like part: only same spin
                                       if (q == r .and. p == s) then
                                             gamma_aa = gamma_aa - n(q)
                                       end if

                                       ! Transition-density contribution
                                       TnuTnu = zero
                                       rrr = zero
                                       do nu = 1, NDim
                                             if (eigs(nu) > small_e .and. eigs(nu) < big_e) then
                                                   TnuTnu = TnuTnu + TNU(p,r,nu) * TNU(s,q,nu)
                                                   rrr = rrr +  TNU(p,r,nu) * TNU(s,q,nu)
                                             end if
                                       end do
                                       
                                       gamma_aa = gamma_aa + TnuTnu
                                       gamma_ab = gamma_ab + TnuTnu
                                       
                                       AuxData%ph_rdm2_aa(p,q,r,s) = gamma_aa
                                       AuxData%ph_rdm2_ab(p,q,r,s) = gamma_ab
                                       
                                       int_aa = Ints%ints2e_aa(gmap_4fold(p,r,q,s,NBasis))
                                       int_ab = Ints%ints2e_ab(gmap_4fold(p,r,q,s,NBasis))
                                       
                                       Err = Err + frac12 * rrr * (int_aa  +  int_ab)
                                       Eee = Eee + (gamma_aa * int_aa  + gamma_ab * int_ab)
                                       
                                       Eeeh = Eeeh +  rdm2h *(int_aa  + int_ab)
                                       Eeex = Eeex +  rdm2x * int_aa

                                 end do
                           end do
                     end do
               end do

               write(*, '(A20, F20.15)') "Err from full ERPA-reconstructed 2RDM", Err
               write(*, '(A20, F20.15)') "Eee from full ERPA-reconstructed 2RDM", Eee
               write(*, '(A20, F20.15)') "Eeeh from full ERPA-reconstructed 2RDM", Eeeh
               write(*, '(A20, F20.15)') "Eeex from full ERPA-reconstructed 2RDM", Eeex

               ! N-rep diagnostyka — bloki spinowe osobno
            block
                double precision :: Na
                Na = real(AuxData%Nel, F64) * 0.5d0   ! N_alpha = N/2 dla closed shell

                ! alpha-alpha: trace Na(Na-1), partial trace -> gamma^a, n_max=1
                call check_2rdm_nrep("phERPA alpha-alpha", AuxData%ph_rdm2_aa, NBasis, &
                                     Na*(Na - 1.0d0), Na - 1.0d0, Na, 1.0d0, .true.)

                ! alpha-beta: trace Na*Nb, partial trace daje N_a*gamma^b, n_max=1
                call check_2rdm_nrep("phERPA alpha-beta", AuxData%ph_rdm2_ab, NBasis, &
                                     Na*Na, Na, Na, 1.0d0, .false.)
            end block

               deallocate(TNU)
               deallocate(skipped)

             end associate

       end subroutine pherpa_gen_rdm

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
      subroutine check_2rdm_nrep(label, Gamma, NBasis, &
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
      end subroutine check_2rdm_nrep


       subroutine pherpa_gen_rdm_fluct(AuxData)

             type(TACppData), intent(inout) :: AuxData
             integer :: p, q, r, s
             integer :: i, nu
             integer :: NBasis, NDim
             integer :: skip
             
             double precision, parameter :: small_e = 1.d-3
             double precision, parameter :: big_e   = 1.d8

             double precision :: TnuTnu
             double precision :: gamma_aa, gamma_ab

             double precision, dimension(:,:,:), allocatable :: Tnu
             double precision, dimension(:), allocatable :: skipped

             associate( &
                   n       => AuxData%Occ,      &
                   indn    => AuxData%IndN,     &
                   eigs    => AuxData%Eigs,     &
                   EigvecX => AuxData%EigvX,    &
                   EigvecY => AuxData%EigvY )
               
               NBasis = AuxData%NBasis
               NDim   = AuxData%NDim


               skip  = 0
               
               allocate(Tnu(NBasis, NBasis, NDim))
               allocate(skipped(NDim))

               Tnu = zero
               skipped = zero

               ! ------------------------------------------------------------
               ! Build full transition-density matrices TNU(p,q,nu)
               !
               ! For pair p > q stored as i:
               !
               ! T_qp^nu = (n(p)-n(q)) Y_pq^nu
               ! T_pq^nu = (n(q)-n(p)) X_pq^nu
               !
               ! T_pp^nu = 0 by construction
               ! ------------------------------------------------------------
               
               do nu = 1, NDim
                     if (eigs(nu) > small_e .and. eigs(nu) < big_e) then
                           do i = 1, NDim
                                 p = indn(1,i)
                                 q = indn(2,i)
                                 
                                 TNU(q,p,nu) = (n(p) - n(q)) * EigvecY(i,nu)
                                 TNU(p,q,nu) = (n(q) - n(p)) * EigvecX(i,nu)

                           end do
                     else
                           skip = skip + 1
                           skipped(skip) = eigs(nu)
                     end if
               end do

               AuxData%ph_rdm2_aa = zero
               AuxData%ph_rdm2_ab = zero
               
               ! ------------------------------------------------------------
               ! Full ERPA-reconstructed 2-RDM
               !
               ! Gamma_pqrs_response = sum_nu T_pr^nu T_sq^nu
               !
               ! ------------------------------------------------------------
               
               do p = 1, NBasis
                     do q = 1, NBasis
                           do r = 1, NBasis
                                 do s = 1, NBasis
                                       
                                       gamma_aa = zero
                                       gamma_ab = zero
                                       

                                       ! Transition-density contribution
                                       TnuTnu = zero
                                       do nu = 1, NDim
                                             if (eigs(nu) > small_e .and. eigs(nu) < big_e) then
                                                   TnuTnu = TnuTnu + TNU(p,r,nu) * TNU(s,q,nu)
                                             end if
                                       end do
                                       
                                       gamma_aa = gamma_aa + TnuTnu
                                       gamma_ab = gamma_ab + TnuTnu
                                       
                                       AuxData%ph_rdm2_aa(p,q,r,s) = gamma_aa
                                       AuxData%ph_rdm2_ab(p,q,r,s) = gamma_ab
                                       !if (abs(AuxData%ph_rdm2_aa(p,q,r,s)).gt.1.d-5)then
                                       !      write(*, '(A10, 4I5, F15.8)') 'fluct', p, q, r, s, AuxData%ph_rdm2_aa(p,q,r,s)
                                       !end if
                                       
                                 end do
                           end do
                     end do
               end do

               deallocate(TNU)
               deallocate(skipped)

             end associate

       end subroutine pherpa_gen_rdm_fluct
  
       function func_T1(AuxData, Ints, r, p, s, q)
        double precision :: func_T1
        type(TInts), intent(inout) :: Ints
        type(TACppData), intent(inout) :: AuxData
        double precision :: num
        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)



              func_T1 = zero

              if (IAux(r)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA                                
                                num = one
                                if (AuxData%int_alpha(gmap_4fold(s,u,t,q, NBasis)) == 0) num = alpha

                                func_T1 = func_T1 +  num * (&
                                rdm2_pp(r - NI, t - NI, u - NI, p - NI) * ints2e_aa(gmap_4fold(s,u,t,q,NBasis)) &
                                + rdm2_pm(r - NI, t - NI, u - NI, p - NI) * ints2e_ab(gmap_4fold(s,u,t,q,NBasis)))
                          end do
                    end do
                    if (p == r) then
                          
                          do t = 1, NI
                                num = one
                                if (AuxData%int_alpha(gmap_4fold(s,t,t,q, NBasis)) == 0) num = alpha

                                func_T1 = func_T1 -  num *&
                                      (n_p(r) * n_p(t) )* ints2e_aa(gmap_4fold(s,t,t,q,NBasis))
                          end do

                    end if
              else                    
                    if (p == r) then
                          if (IAux(p) == 0)then
                                do t = 1, NIA
                                      num = one
                                      if (AuxData%int_alpha(gmap_4fold(s,t,t,q, NBasis)) == 0) num = alpha

                                      func_T1 = func_T1 - num * &
                                            (n_p(r) * n_p(t) )* ints2e_aa(gmap_4fold(s,t,t,q,NBasis))
                                end do
                          end if
                    end if
                    
                    if (IAux(r) == 0 .or. IAux(p) == 0) then
                          num = one
                          if (AuxData%int_alpha(gmap_4fold(s,r,p,q, NBasis)) == 0) num = alpha
                                
                          func_T1 = func_T1 + num*(&
                                (n_p(r) * n_p(p)) * ints2e_aa(gmap_4fold(s,r,p,q,NBasis)) +&
                                (n_p(r) * n_m(p)) * ints2e_ab(gmap_4fold(s,r,p,q,NBasis)) )
                    end if
                    
              end if
     
            end associate
      end function func_T1

      function func_T2(AuxData, Ints, r, p, s, q)
            double precision :: func_T2
            type(TInts), intent(inout) :: Ints
            type(TACppData), intent(inout) :: AuxData
            double precision :: num
        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)



              func_T2 = zero

              if (IAux(r)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                num = one
                                if (AuxData%int_alpha(gmap_4fold(s,q,t,u, NBasis)) == 0) num = alpha

                                func_T2 = func_T2 +  num * (&
                                rdm2_pp(r - NI, t - NI, p - NI, u - NI) * ints2e_aa(gmap_4fold(s,q,t,u,NBasis)) &
                                + rdm2_pm(r - NI, t - NI, p - NI, u - NI) * ints2e_ab(gmap_4fold(s,q,t,u,NBasis)))
                          end do
                    end do

                    if (p == r) then                                                                   
                          do t = 1, NI
                                num = one
                                if (AuxData%int_alpha(gmap_4fold(s,q,t,t, NBasis)) == 0) num = alpha

                                func_T2 = func_T2 + num *(&
                                      (n_p(r) * n_p(t) )* ints2e_aa(gmap_4fold(s,q,t,t,NBasis)) + &
                                      (n_p(r) * n_m(t) )* ints2e_ab(gmap_4fold(s,q,t,t,NBasis)) )
                          end do
                    end if
              else                    
                    if (p == r) then
                          if (IAux(p) == 0)then
                                do t = 1, NIA
                                      num = one
                                      if (AuxData%int_alpha(gmap_4fold(s,q,t,t, NBasis)) == 0) num = alpha

                                      func_T2 = func_T2 + num*(&
                                            (n_p(r) * n_p(t) )* ints2e_aa(gmap_4fold(s,q,t,t,NBasis)) + &
                                            (n_p(r) * n_m(t) )* ints2e_ab(gmap_4fold(s,q,t,t,NBasis)) )
                                end do
                          end if
                    end if
                    
                    if (IAux(r) == 0 .or. IAux(p) == 0) then
                          num = one
                          if (AuxData%int_alpha(gmap_4fold(s,q,p,r, NBasis)) == 0) num = alpha

                          func_T2 = func_T2 -  num*&
                                (n_p(r) * n_p(p)) * ints2e_aa(gmap_4fold(s,q,p,r,NBasis))
                    end if
              end if     
            end associate
      end function func_T2

      function func_T34(AuxData, Ints, s, p, r, q)
            double precision :: func_T34
            type(TInts), intent(inout) :: Ints
            type(TACppData), intent(inout) :: AuxData
            double precision :: num

        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              ints2e_aa=>Ints%ints2e_aa, ints2e_ab=>Ints%ints2e_ab, &
              ints1e_aa=>Ints%ints1e_aa, Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)



              func_T34 = zero
              if (IAux(s)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                num = one
                                if (AuxData%int_alpha(gmap_4fold(t,r,u,q, NBasis)) == 0) num = alpha

                                func_T34 = func_T34 + num * (&
                                rdm2_pp(t - NI, u - NI, s - NI, p - NI) * ints2e_aa(gmap_4fold(t,r,u,q,NBasis)) &
                                + rdm2_pm(t - NI, u - NI, s - NI, p - NI) * ints2e_ab(gmap_4fold(t,r,u,q,NBasis)))
                          end do
                    end do
              else
                    num = one
                    if (AuxData%int_alpha(gmap_4fold(s,r,p,q, NBasis)) == 0) num = alpha

                    func_T34 = func_T34 +num * (&
                          (n_p(s) * n_p(p)) * ints2e_aa(gmap_4fold(s,r,p,q,NBasis)) + &
                          (n_p(s) * n_m(p)) * ints2e_ab(gmap_4fold(s,r,p,q,NBasis)) - &
                          (n_p(p) * n_p(s)) * ints2e_aa(gmap_4fold(p,r,s,q,NBasis)))

              end if

            end associate
            
      end function func_T34


      function func_T1_ducc(AuxData, Ints, r, p, s, q)
        double precision :: func_T1_ducc
        type(TInts), intent(inout) :: Ints
        type(TACppData), intent(inout) :: AuxData
        double precision :: num
        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              Aints2e_aa=>Ints%Aints2e_aa, Aints2e_ab=>Ints%Aints2e_ab, &
              Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)

          num = one
          func_T1_ducc = zero
          if (IAux(r)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA                                

                                func_T1_ducc = func_T1_ducc +  num * (&
                             frac12 * rdm2_pp(r - NI, t - NI, u - NI, p - NI) * Aints2e_aa(gmap_4fold(s,u,t,q,NBasis)) &
                             + rdm2_pm(r - NI, t - NI, u - NI, p - NI) * Aints2e_ab(gmap_4fold(s,u,t,q,NBasis)))
                          end do
                    end do

                    if (p == r) then                                                                   
                          do t = 1, NI
                                func_T1_ducc = func_T1_ducc -  num *&
                                      frac12 * (n_p(r) * n_p(t) )* Aints2e_aa(gmap_4fold(s,t,t,q,NBasis))
                          end do
                    end if
              else                    
                    if (p == r) then
                          if (IAux(p) == 0)then
                                do t = 1, NIA

                                      func_T1_ducc = func_T1_ducc - num * &
                                            frac12 * (n_p(r) * n_p(t) )* Aints2e_aa(gmap_4fold(s,t,t,q,NBasis))
                                end do
                          end if

                    end if
                    
                    if (IAux(r) == 0 .or. IAux(p) == 0) then
                                
                          func_T1_ducc = func_T1_ducc + num*(&
                                frac12 * (n_p(r) * n_p(p)) * Aints2e_aa(gmap_4fold(s,r,p,q,NBasis)) +&
                                (n_p(r) * n_m(p)) * Aints2e_ab(gmap_4fold(s,r,p,q,NBasis)) )
                    end if
              end if
     
            end associate
      end function func_T1_ducc

      function func_T2_ducc(AuxData, Ints, r, p, s, q)
            double precision :: func_T2_ducc
            type(TInts), intent(inout) :: Ints
            type(TACppData), intent(inout) :: AuxData
            double precision :: num
        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              Aints2e_aa=>Ints%Aints2e_aa, Aints2e_ab=>Ints%Aints2e_ab, &
              Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)

          num = one

              func_T2_ducc = zero

              if (IAux(r)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA
                                func_T2_ducc = func_T2_ducc +  num * (&
                                frac12 * rdm2_pp(r - NI, t - NI, p - NI, u - NI) * Aints2e_aa(gmap_4fold(s,q,t,u,NBasis)) &
                                + rdm2_pm(r - NI, t - NI, p - NI, u - NI) * Aints2e_ab(gmap_4fold(s,q,t,u,NBasis)))
                          end do
                    end do

                    if (p == r) then                                                                   
                          do t = 1, NI
                                func_T2_ducc = func_T2_ducc + num *(&
                                      frac12 * (n_p(r) * n_p(t) )* Aints2e_aa(gmap_4fold(s,q,t,t,NBasis)) + &
                                      (n_p(r) * n_m(t) )* Aints2e_ab(gmap_4fold(s,q,t,t,NBasis)) )
                          end do
                    end if
              else                    
                    if (p == r) then
                          if (IAux(p) == 0)then
                                do t = 1, NIA
                                      func_T2_ducc = func_T2_ducc + num*(&
                                            frac12 * (n_p(r) * n_p(t) )* Aints2e_aa(gmap_4fold(s,q,t,t,NBasis)) + &
                                            (n_p(r) * n_m(t) )* Aints2e_ab(gmap_4fold(s,q,t,t,NBasis)) )
                                end do
                          end if
                    end if
                    
                    if (IAux(r) == 0 .or. IAux(p) == 0) then

                          func_T2_ducc = func_T2_ducc -  num*&
                                frac12 * (n_p(r) * n_p(p)) * Aints2e_aa(gmap_4fold(s,q,p,r,NBasis))
                    end if
              end if     
            end associate
      end function func_T2_ducc

      function func_T34_ducc(AuxData, Ints, s, p, r, q)
            double precision :: func_T34_ducc
            type(TInts), intent(inout) :: Ints
            type(TACppData), intent(inout) :: AuxData
            double precision :: num

        integer, intent(in) :: r, p, s, q
        integer :: t, u

        associate(rdm2_pm=>AuxData%rdm2_pm, rdm2_mp=>AuxData%rdm2_mp, &
              rdm2_pp=>AuxData%rdm2_pp, rdm2_mm=>AuxData%rdm2_mm,&
              NIA=>AuxData%NIA, NI=>AuxData%NI, &
              IAux=>AuxData%IndAux, n_p=>AuxData%n_p, n_m=>AuxData%n_m, &
              NBasis=>AuxData%NBasis, &
              Aints2e_aa=>Ints%Aints2e_aa, Aints2e_ab=>Ints%Aints2e_ab, &
              Aints1e_aa=>Ints%Aints1e_aa, &
              int_alpha=> AuxData%int_alpha, alpha=>AuxData%alpha)

          
          num = one
              func_T34_ducc = zero
              if (IAux(s)==(1).and.IAux(p)==(1))then

                    do t = NI+1, NIA
                          do u = NI+1, NIA

                                func_T34_ducc = func_T34_ducc + num * (&
                                frac12 * rdm2_pp(t - NI, u - NI, s - NI, p - NI) * Aints2e_aa(gmap_4fold(t,r,u,q,NBasis)) &
                                + rdm2_pm(t - NI, u - NI, s - NI, p - NI) * Aints2e_ab(gmap_4fold(t,r,u,q,NBasis)))
                          end do
                    end do
              else

                    func_T34_ducc = func_T34_ducc +num * (&
                          frac12 * (n_p(s) * n_p(p)) * Aints2e_aa(gmap_4fold(s,r,p,q,NBasis)) + &
                          (n_p(s) * n_m(p)) * Aints2e_ab(gmap_4fold(s,r,p,q,NBasis)) - &
                          frac12 * (n_p(p) * n_p(s)) * Aints2e_aa(gmap_4fold(p,r,s,q,NBasis)))
              end if

            end associate
            
      end function func_T34_ducc


      subroutine twoel_8fold_to_ints2e_4fold(Ints, NBasis)
            type(TInts), intent(inout) :: Ints
            integer, intent(in) :: NBasis

            integer :: p, q, r, s

            Ints%ints2e_dim = int(NBasis,8)**2 * (int(NBasis,8)**2 + 1) / 2

            if (allocated(Ints%ints2e_aa)) deallocate(Ints%ints2e_aa)
            if (allocated(Ints%ints2e_ab)) deallocate(Ints%ints2e_ab)

            allocate(Ints%ints2e_aa(Ints%ints2e_dim))
            allocate(Ints%ints2e_ab(Ints%ints2e_dim))

            Ints%ints2e_aa = zero
            Ints%ints2e_ab = zero

            do p = 1, NBasis
                  do q = 1, NBasis
                        do r = 1, NBasis
                              do s = 1, NBasis
                                    Ints%ints2e_aa(gmap_4fold(p, q, r, s, NBasis)) = Ints%ints2e(gmap(p, q, r, s))
                                    Ints%ints2e_ab(gmap_4fold(p, q, r, s, NBasis)) = Ints%ints2e(gmap(p, q, r, s))
                              end do
                        end do
                  end do
            end do

      end subroutine twoel_8fold_to_ints2e_4fold



      subroutine restore_XY(AuxData)
         type(TACppData), intent(inout) :: AuxData

         integer :: i, k, p, q
         double precision :: X_val, Y_val, omega, norm
         double precision, parameter :: small = 1.d-6

         associate(ndim=>AuxData%NDim, ABmin=>AuxData%ABmin, &
               EigvecY=>AuxData%Eigvec, EigvX=>AuxData%EigvX, &
               EigvY=>AuxData%EigvY, &
               Eigs=>AuxData%Eigs, IndN=>AuxData%IndN, cc=>AuxData%cc, n=>AuxData%Occ)


           EigvX = zero
           EigvY = EigvecY
           do k = 1, ndim
                 omega = eigs(k)
                 !                 print*, k, 'omega', omega
                 if (omega > small) then
                       !                      print*, 'yes'
                       ! X_col = (1/omega) * ABmin * Y_col
                       call real_Av(EigvX(:, k), ABmin, EigvY(:, k))
                       EigvX(:, k) = EigvX(:, k) / omega
                 else
                       EigvX(:, k) = zero
                       EigvY(:, k) = zero 
                 end if
           end do

           ! 3. Transform X_symm, Y_symm -> X, Y
           ! X = X_symm / (c_p + c_q)
           ! Y = Y_symm / (c_p - c_q)
           ! EigvecX = 0.5 * (X - Y)
           ! EigvecY = 0.5 * (X + Y)

           do k = 1, ndim
                 do i = 1, ndim
                       p = IndN(1, i)
                       q = IndN(2, i)

                       X_val = zero
                       Y_val = zero

                       if (abs(cc(p) + cc(q)) > small) then
                             X_val = EigvX(i, k) / (cc(p) + cc(q))
                       end if

                       if (abs(cc(p) - cc(q)) > small) then
                             Y_val = EigvY(i, k) / (cc(p) - cc(q))
                       end if

                       EigvX(i, k) = frac12 * (X_val - Y_val)
                       EigvY(i, k)  = frac12 * (X_val + Y_val)


                 end do
           end do

           do k = 1, ndim
                 norm = 0.0d0
                 do i = 1, ndim
                       p = IndN(1, i)
                       q = IndN(2, i)
                       norm = norm + (n(p) - n(q)) * (EigvY(i, k)**2 - EigvX(i, k)**2)
                       ! if (k==172)then
                       !       write(*,'(5F12.6)') EigvY(i, k)**2, EigvX(i, k)**2,  (EigvY(i, k)**2 - EigvX(i, k)**2), (n(p) - n(q)), norm
                       ! end if
                 end do
                 ! print *, "For vector k=", k, "normalization is", norm
           end do

         end associate

   end subroutine restore_XY

   subroutine write_fcidump_alpha(filename, AuxData, Ints)
        character(len=*), intent(in) :: filename
        type(TACppData), intent(in) :: AuxData
        type(TInts), intent(in) :: Ints

        integer :: p, q, r, s
        integer :: unit, ios
        integer :: s_max
        integer(I8) :: idx
        double precision, parameter :: thresh = 1.0d-9

        unit = 20
        open(unit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            print*, "Error opening file ", filename
            return
        end if

        ! Write Header
        write(unit, '(A, I0, A, I0, A)', advance='no') " &FCI NORB=", AuxData%NBasis, ",NELEC=", AuxData%NEL, ",MS2=0,"
        write(unit, '()') ! newline
        write(unit, '(A)', advance='no') "  ORBSYM="
        do p = 1, AuxData%NBasis
            write(unit, '(I0, A)', advance='no') 1, ","
        end do
        write(unit, '()') ! newline
        write(unit, '(A)') "  ISYM=1,"
        write(unit, '(A)') " &END"

        ! Block 1: 2e_aa
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        if (abs(Ints%Aints2e_aa(idx)) > thresh) then
                            write(unit, '(ES22.14, 4I5)') Ints%Aints2e_aa(idx), p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.14, 4I5)') 0.0d0, 0, 0, 0, 0
        
        ! Block 2: 2e_bb (same as aa for restricted spatial orbitals, since FCIDUMP-full has 5 blocks)
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        if (abs(Ints%Aints2e_aa(idx)) > thresh) then
                            write(unit, '(ES22.14, 4I5)') Ints%Aints2e_aa(idx), p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 3: 2e_ab
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        if (abs(Ints%Aints2e_ab(idx)) > thresh) then
                            write(unit, '(ES22.14, 4I5)') Ints%Aints2e_ab(idx), p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 4: 1e_aa
        do p = 1, AuxData%NBasis
            do q = 1, p
                if (abs(Ints%Aints1e_aa(p, q)) > thresh) then
                    write(unit, '(ES22.14, 4I5)') Ints%Aints1e_aa(p, q), p, q, 0, 0
                end if
            end do
        end do
        write(unit, '(ES22.14, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 5: 1e_bb
        do p = 1, AuxData%NBasis
            do q = 1, p
                if (abs(Ints%Aints1e_aa(p, q)) > thresh) then
                    write(unit, '(ES22.14, 4I5)') Ints%Aints1e_aa(p, q), p, q, 0, 0
                end if
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Core energy
        write(unit, '(ES22.14, 4I5)') AuxData%ENuc, 0, 0, 0, 0

        close(unit)

  end subroutine write_fcidump_alpha
  
   subroutine write_fcidump_alpha_spinres(filename, AuxData, Ints)
        character(len=*), intent(in) :: filename
        type(TACppData), intent(in) :: AuxData
        type(TInts), intent(in) :: Ints

        integer :: p, q, r, s
        integer :: unit, ios
        integer :: s_max
        integer(I8) :: idx
        double precision :: alpha_val, int_val
        double precision, parameter :: thresh = 1.0d-9

        unit = 20
        open(unit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            print*, "Error opening file ", filename
            return
        end if

        ! Write Header
        write(unit, '(A, I0, A, I0, A)', advance='no') " &FCI NORB=", AuxData%NBasis, ",NELEC=", AuxData%NEL, ",MS2=0,"
        write(unit, '()') ! newline
        write(unit, '(A)', advance='no') "  ORBSYM="
        do p = 1, AuxData%NBasis
            write(unit, '(I0, A)', advance='no') 1, ","
        end do
        write(unit, '()') ! newline
        write(unit, '(A)') "  ISYM=1,"
        write(unit, '(A)') " &END"

        ! Block 1: 2e_aa
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        
                        alpha_val = 1.0d0
                        if (AuxData%int_alpha(idx) == 0) alpha_val = AuxData%alpha
                        int_val = alpha_val * Ints%ints2e_aa(idx)
                        
                        if (abs(int_val) > thresh) then
                            write(unit, '(ES22.14, 4I5)') int_val, p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0
        
        ! Block 2: 2e_bb (same as aa for regular spatial/spinres path)
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        
                        alpha_val = 1.0d0
                        if (AuxData%int_alpha(idx) == 0) alpha_val = AuxData%alpha
                        int_val = alpha_val * Ints%ints2e_aa(idx)
                        
                        if (abs(int_val) > thresh) then
                            write(unit, '(ES22.14, 4I5)') int_val, p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 3: 2e_ab
        do p = 1, AuxData%NBasis
            do q = 1, p
                do r = 1, p
                    if (p == r) then
                        s_max = q
                    else
                        s_max = r
                    end if
                    do s = 1, s_max
                        idx = gmap_4fold(p, q, r, s, AuxData%NBasis)
                        
                        alpha_val = 1.0d0
                        if (AuxData%int_alpha(idx) == 0) alpha_val = AuxData%alpha
                        int_val = alpha_val * Ints%ints2e_ab(idx)
                        
                        if (abs(int_val) > thresh) then
                            write(unit, '(ES22.14, 4I5)') int_val, p, q, r, s
                        end if
                    end do
                end do
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 4: 1e_aa
        do p = 1, AuxData%NBasis
            do q = 1, p
                if (abs(Ints%Aints1e_aa(p, q)) > thresh) then
                    write(unit, '(ES22.14, 4I5)') Ints%Aints1e_aa(p, q), p, q, 0, 0
                end if
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Block 5: 1e_bb
        do p = 1, AuxData%NBasis
            do q = 1, p
                if (abs(Ints%Aints1e_aa(p, q)) > thresh) then
                    write(unit, '(ES22.14, 4I5)') Ints%Aints1e_aa(p, q), p, q, 0, 0
                end if
            end do
        end do
        write(unit, '(ES22.1, 4I5)') 0.0d0, 0, 0, 0, 0

        ! Core energy
        write(unit, '(ES22.1, 4I5)') AuxData%ENuc, 0, 0, 0, 0

        close(unit)

   end subroutine write_fcidump_alpha_spinres

   end module phac_spinres
