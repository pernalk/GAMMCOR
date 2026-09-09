module interface_pp

      use ppac_types
      use math_constants
      use basis_sets
      use Cholesky_Gammcor
      use THC_Gammcor
      use OneElectronInts_Gammcor
      use sys_definitions
      use gammcor_integrals
      use geom_input
      use real_linalg
      use sort
      use clock
      use print_utils
      use h5_reader


      implicit none


      character(len=*), parameter :: PYSCF_H5 = 'pyscf_data.h5'
      integer(HID_T) :: pyscf_fid = -1

contains

    subroutine pyscf_h5_open()
          if (pyscf_fid /= -1) return
          call h5_init()
          pyscf_fid = h5_open(PYSCF_H5)
    end subroutine pyscf_h5_open

    subroutine pyscf_h5_close()
          if (pyscf_fid == -1) return
          call h5_close(pyscf_fid)
          call h5_finish()
          pyscf_fid = -1
    end subroutine pyscf_h5_close

    ! "_1.1" (or "" for a single-state job) -> "/POSTHF/STATES/1.1"
    function pyscf_state_group(suffix) result(grp)
          character(len=*), intent(in) :: suffix
          character(len=64) :: grp
          if (len_trim(suffix) == 0) then
                grp = '/POSTHF/STATES/1.1'
          else
                grp = '/POSTHF/STATES/'//suffix(2:len_trim(suffix))
          end if
    end function pyscf_state_group


    ! NBasis for the PySCF interface.  Called from mainp before read_PYSCF, so
    ! it opens and closes the container itself.  Reads /REF/NBASIS when the
    ! HDF5 file is there, otherwise falls back to a legacy auxdata<suffix>.txt.
    subroutine basinfo_pyscf(nbasis)
          integer, intent(out) :: nbasis

          character(len=256) :: filename
          character(len=:), allocatable :: command
          character(len=15) :: temp_file = 'list.tmp'
          integer :: unit, io_status
          integer(HID_T) :: fid
          integer :: iv(1)
          logical :: h5_present

          inquire(file=PYSCF_H5, exist=h5_present)
          if (h5_present) then
                call h5_init()
                fid = h5_open(PYSCF_H5)
                call h5_get(fid, '/REF/NBASIS', iv)
                call h5_close(fid)
                call h5_finish()
                nbasis = iv(1)
                return
          end if

          filename = ''
          command = 'ls -1 auxda*.txt 2>/dev/null | head -n 1 >' // temp_file
          call execute_command_line(command, wait=.true.)

          open(newunit=unit, file=temp_file, status='old', &
               action='read', iostat=io_status)
          if (io_status == 0) then
                read(unit, '(A)', iostat=io_status) filename
                if (io_status /= 0) filename = ''
                close(unit, status='delete')
          end if

          if (trim(filename) == '') then
                print *, "ERROR: found neither ", PYSCF_H5, &
                         " nor a legacy auxdata*.txt file"
                stop
          end if

          open(newunit=unit, file=filename, status='old', action='read', &
               iostat=io_status)
          if (io_status /= 0) then
                print *, "ERROR: cannot open ", trim(filename)
                stop
          end if
          read(unit, *, iostat=io_status) nbasis
          close(unit)
          if (io_status /= 0) then
                print *, "ERROR: no NBasis in ", trim(filename)
                stop
          end if

    end subroutine basinfo_pyscf


    ! NBasis for the ORCA interface, taken from the NORB= field of the
    ! FCIDUMP header.
    subroutine basinfo_orca(nbasis)
          integer, intent(out) :: nbasis
          character(len=256) :: line, filename
          integer :: pos, ios, unit
          logical :: filex1, filex2

          inquire(file='FCIDUMP', exist=filex1)
          inquire(file='FCIDUMP-full', exist=filex2)

          if (filex1) then
                filename = 'FCIDUMP'
          else if (filex2) then
                filename = 'FCIDUMP-full'
          else
                print*, 'NO FCIDUMP file'
                print*, 'interface_pp.f90/basinfo_orca'
                stop
          end if

          open(newunit=unit, file=filename, status="old", &
               action="read", iostat=ios)
          if (ios /= 0) then
                print *, "ERROR: cannot open ", trim(filename)
                stop
          end if
          read(unit, fmt="(a)", iostat=ios) line
          close(unit)
          if (ios /= 0) then
                print *, "FCIDUMP empty"
                stop
          end if

          line = trim(adjustl(line))
          pos = index(line, "NORB=")
          if (pos == 0) then
                print *, "ERROR: no NORB= in the ", trim(filename), " header"
                stop
          end if
          read(line(pos+5:), *, iostat=ios) nbasis
          if (ios /= 0) then
                print *, "ERROR: cannot parse NORB= in ", trim(filename)
                stop
          end if

          print*, 'nbasis', nbasis

    end subroutine basinfo_orca


      !subroutine interface_driver(THCData, AuxData, CAONO, Flags, TwoEl)
      !end subroutine interface_driver


      subroutine read_PYSCF(THCData, AuxData, CAONO, Flags, TwoEl)
            

            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), allocatable, intent(out) :: CAONO
            type(FlagsData), intent(in) :: Flags
            double precision, optional, intent(inout) :: TwoEl(:)

            double precision :: ECAS
            integer :: nbasis, ninactive, nactive, nvirtual
            character(len=128) :: line

            integer :: unit
            integer :: i, j, k, l
            
            double precision, dimension(:), allocatable :: occ_temp
            double precision, dimension(:, :), allocatable :: CAOMO
            double precision, dimension(:, :), allocatable :: work
            double precision, dimension(:, :), allocatable :: mo_occ_int
            logical :: STATE_AV = .false.
            integer :: natural_orb
            integer :: anst(2)
            logical :: existsNosf, existsSf
            integer :: iostat
            character(len=10) :: suffix
            integer :: x, y
            integer :: count_suffix_files
            logical :: file_exists
            character(len=80) :: separator
            !character(len=256) :: rdm1_full_name, rdm2_full_name
            integer, external :: NAddr3

            
            type(TAOBASIS) :: AObasis
            type(TSystem)  :: System
            integer :: ms2

            call pyscf_h5_open()

            if (Flags%ITwoEl.ne.1)then
                  print*, 'Two electron integrals calculated with THC'
                  call GEOM_init(Flags, AuxData, AObasis, System)
            else
                  print*, 'Two electron integrals read from ', PYSCF_H5
                  ! call read_fcidump_header("FCIDUMP", AuxData%NBasis, AuxData%Nel, ms2)
                  ! call read_fcidump_restricted('FCIDUMP', AuxData, TwoEl)

                  call h5_get(pyscf_fid, '/MOINTS/ERI', TwoEl)
            end if

            !-----------------------------------------------------------------------------------------------------
            ! Check if there is a state specified in input
            !-----------------------------------------------------------------------------------------------------
            call print_section('State and job setup')
            call print_info('AuxData%nst', AuxData%nst)
            call print_info('Flags%JobType', Flags%JobType)

            if (Flags%JOBTYPE == JOB_TYPE_MP2 .or. Flags%JOBTYPE == JOB_TYPE_SRMP2) then
                  call load_mp2_aux_data(AuxData)
                  AuxData%PYSCF = 1
            else
                  call resolve_state_suffix(AuxData, suffix)
                  call load_cas_aux_data(suffix, AuxData, natural_orb)
                  AuxData%PYSCF = 1

                  call print_section('CASSCF data loaded')
                  call print_info('NBasis', AuxData%nbasis)
                  call print_info('NInactive', AuxData%NI)
                  call print_info('NActive', AuxData%NA)
                  call print_info('NVirtual', AuxData%NV)
                  call print_info('CASSCF Energy', AuxData%ECAS)
                  call print_info('Nuclear Energy', AuxData%ENuc)
                  call print_info('Number of electrons', AuxData%Nel)
                  call print_info('Natural orbitals used (0/1)', natural_orb)
                  
            end if
            

            associate(NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NBasis=>AuxData%NBasis, NIA=>AuxData%NIA)

              allocate(AuxData%HNO0(NBasis, NBasis))
              allocate(THCData%HNO(NBasis, NBasis))
              THCData%ExternalOrdering = ORBITAL_ORDERING_PYSCF
              THCData%H0external = Flags%H0external
              call print_info('THCData%H0external', merge('TRUE ', 'FALSE', THCData%H0external))
              allocate(THCData%fij(NI))
              allocate(THCData%fvw(NV))

              call load_one_electron_integrals(AuxData, THCData, CAOMO, CAONO, natural_orb, Flags)

              ! do i = 1, nbasis
              !       do j = 1, nbasis
              !             if (abs(CAONO(i,j)).gt.1.d-5)then
              !                   print*, 'caono', i, j, caono(i,j)
              !             end if
              !       end do
              ! end do

              if (Flags%JOBTYPE /= JOB_TYPE_MP2 .and. Flags%JOBTYPE /= JOB_TYPE_SRMP2) then

                    if (Flags%JOBTYPE == JOB_TYPE_AC0.and.Flags%ITwoEl > 1)then
                          call load_density_matrices(suffix, AuxData, THCData, CAOMO, CAONO, Flags, AObasis, System, TwoEl, .true.)
                    else
                          call load_density_matrices(suffix, AuxData, THCData, CAOMO, CAONO, Flags, AObasis, System, TwoEl,  .false.)

                    
                    ! print*, 'rdmyy'
                    ! write(*, '("p", T6, "q", T11, "r", T16, "s", T25, "rdm_pp", T40, "rdm_mm", T55, "rdm_pm", T70, "rdm_mp", T85, "rdm_full")')
                    ! associate(NA=>AuxData%NA)
                    !   do i = 1, NA
                    !         do j = 1, NA
                    !               do k = 1, NA
                    !                     do l = 1, NA
                    !                           if (abs(AuxData%rdm2_pp(i,j,k,l)) > 1.d-5 .or. &
                    !                                 abs(AuxData%rdm2_mm(i,j,k,l)) > 1.d-5 .or. &
                    !                                 abs(AuxData%rdm2_pm(i,j,k,l)) > 1.d-5 .or. &
                    !                                 abs(AuxData%rdm2_mp(i,j,k,l)) > 1.d-5) then

                    !                                 write(*, '(4I5, 5F15.8)') i, j, k, l, &
                    !                                       AuxData%rdm2_pp(i,j,k,l), &
                    !                                       AuxData%rdm2_mm(i,j,k,l), &
                    !                                       AuxData%rdm2_pm(i,j,k,l), &
                    !                                       AuxData%rdm2_mp(i,j,k,l), &
                    !                                       AuxData%rdm2_full(i,j,k,l)
                    !                           end if
                    !                     end do
                    !               end do
                    !         end do
                    !   end do
                    ! end associate

                    end if
                    if (natural_orb == 0)then
                         
                          print*, 'transform to natural orbital basis'
                          if (present(TwoEl)) then
                                ! --- TWOEL and HCOre
                                call Trans2NO_AO(THCData, AuxData, &
                                      CAOMO, CAONO, Flags, TwoEl=TwoEl)
                                ! --- FCIDUMP modification: call MO version instead of AO version ---
                                ! call Trans2NO_MO(THCData, AuxData, &
                                !       CAONO, Flags, TwoEl=TwoEl)
                                

                          else
                                call Trans2NO_AO(THCData, AuxData, &
                                      CAOMO, CAONO, Flags)
                                
!                                AuxData%HNO0 = AuxData%HNO0_THC
                          end if
                          print*, 'writing rdm2 full'
                          call write_rdm2_dat(AuxData%rdm2_full, NA)
                    else
                          print*, 'natural orbitals'
                          call load_occupancy(AuxData)
                          ! No MO->NO transformation is needed, but the later code
                          ! (e.g. AB_CAS_FOFO) still expects rdm2.dat on disk, so
                          ! dump the 2-RDM read from the h5 file here as well.
                          if (allocated(AuxData%rdm2_full)) then
                                print*, 'writing rdm2 full (natural orbitals)'
                                call write_rdm2_dat(AuxData%rdm2_full, NA)
                          end if
                    end if                          
              else
                    call load_rohf_occupation(AuxData)
              end if


              
              if (Flags%ITwoEl > 1)then
                    call THC_init2(Flags, THCData, AuxData,  AObasis, System, CAONO)
                    AuxData%HNO0 = AuxData%HNO0_THC
                    THCData%HNO = AuxData%HNO0
                    
                    ! do i = 1, NBasis
                    !       do j = 1, NBasis
                    !             if (abs(AuxData%HNO0(i,j)).gt.1.d-5)then
                    !                   print*, i, j, AuxData%HNO0(i,j)
                    !             end if
                    !       end do
                    ! end do
                    ! stop

              else
                    allocate(work(NBasis, NBasis))
                    ! do i = 1, NBasis
                    !       do j = 1, NBasis
                    !             if (abs(AuxData%HNO0(i,j)).gt.1.d-5)then
                    !                   write(*, '(A5, 2I3, F12.6)') 'hao', i, j, AuxData%HNO0(i,j)
                    !             end if
                    !       end do
                    ! end do

                    ! --- FCIDUMP modification: skip double transformation ---
                    call real_ab(work, AuxData%HNO0, CAONO)
                    call real_atb(AuxData%HNO0, CAONO, work)
                    ! ------------------------------------------------------
                    THCData%HNO = AuxData%HNO0
                    print*, ''
                    ! do i = 1, NBasis
                    !       do j = 1, NBasis
                    !             if (abs(AuxData%HNO0(i,j)).gt.1.d-5)then
                    !                   write(*, '(A5, 2I3, F12.6)') 'hno', i, j, AuxData%HNO0(i,j)
                    !             end if
                    !       end do
                    ! end do
!                    stop
                    call check_energy_incore(AuxData, TwoEl, spin_sep =.false.)
                    call check_energy_incore(AuxData, TwoEl, spin_sep =.true.)

                    AuxData%spinsep = .true.
                    do i = 1, NBasis
                          if (abs(AuxData%n_p(i)-AuxData%n_m(i)).gt.1.d-8)then
                                AuxData%spinsep = .false.
                          end if
                          write(*, '(A15, I3, 2F12.7)') 'occupancy', i, AuxData%n_p(i),  AuxData%n_m(i)
                    end do
              end if

          end associate

            call pyscf_h5_close()
            
    end subroutine read_PYSCF

    subroutine read_PYSCF_spinres(THCData, AuxData, CAONO, Flags, Ints)
            

            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), allocatable, intent(out) :: CAONO
            type(FlagsData), intent(in) :: Flags
            type(TInts), intent(inout) :: Ints

            double precision :: ECAS
            integer :: nbasis, ninactive, nactive, nvirtual
            character(len=128) :: line

            integer :: unit
            integer :: i, j, k, l
            
            double precision, dimension(:), allocatable :: occ_temp
            double precision, dimension(:, :), allocatable :: CAOMO
            double precision, dimension(:, :), allocatable :: work
            double precision, dimension(:, :), allocatable :: mo_occ_int
            logical :: STATE_AV = .false.
            integer :: natural_orb
            integer :: anst(2)
            logical :: existsNosf, existsSf
            integer :: iostat
            character(len=10) :: suffix
            integer :: x, y
            integer :: count_suffix_files
            logical :: file_exists
            character(len=80) :: separator
            !character(len=256) :: rdm1_full_name, rdm2_full_name

            
            type(TAOBASIS) :: AObasis
            type(TSystem)  :: System
      end subroutine read_PYSCF_spinres

    subroutine PYSCF_wrapper(THCData, AuxData, NInte1, Ninte2, NBasis, &
          ENuc, CAONO, XKin, TwoEl, Occ, NAc, NInAc, anSt, Flags)

          type(TTHCData), intent(out) :: THCData
          type(TACppData) :: AuxData
          integer, intent(in) :: NInte1, Ninte2, NBasis
          integer, intent(out) ::  NAc, NInAc
          double precision, intent(out) :: ENuc
          double precision, dimension(NInte1), intent(out) :: XKin
          double precision, dimension(NInte2), intent(out) :: TwoEl
          double precision, dimension(NBasis), intent(out) :: Occ
          double precision, dimension(NBasis, NBasis), intent(out) :: CAONO
          type(FlagsData), intent(in) :: Flags
          double precision, dimension(:,:), allocatable  :: CAONO_PYSCF
          integer, dimension(:), intent(in) :: anSt
          integer :: unit

          integer :: i, j, ab
          integer, external :: NAddr3 

          AuxData%nst = anSt
          print*, 'inside wrapper'
          print*, Flags%ITwoEl
          AuxData%NInte2 = NInte2
          call read_PYSCF(THCData, AuxData, CAONO_PYSCF, Flags, TwoEl)
          print*, 'zzz', AuxData%ECAS, AuxData%Enuc
          ! if (Flags%ITwoEl==1)then
          !       AuxData%NInte2 = NInte2
          !       unit = 21
          !       open(unit=unit, file='TWOEl.bin', status='old', access='stream', form='unformatted')
          !       read(unit) TwoEl
          !       close(unit)
          !       print*, 'jakto'
          !       call read_PYSCF(THCData, AuxData, CAONO_PYSCF, Flags, TwoEl)
                
          ! else
          !       call read_PYSCF(THCData, AuxData, CAONO_PYSCF, Flags)
          ! end if
          print*, 'po read'
          CAONO = CAONO_PYSCF
          ENuc = AuxData%ENuc
          Occ = AuxData%Occ
          print*, 'AuxData%Occ', AuxData%Occ
          NAc = AuxData%NA
          NInAc = AuxData%NI
!          print*, 'aha1', occ
          ! One-electron integrals are now transformed to NO
          do i = 1, Nbasis
                do j = 1, i
                      ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                      XKin(ab) = AuxData%HNO0(i,j)
                      ! if (i==j)then
                      !       print*, 'jjj', j, AuxData%HNO0(i,j)
                      ! end if
                end do
          end do
          
    end subroutine PYSCF_wrapper


    subroutine ORCA_wrapper(THCData, AuxData, NInte1, Ninte2, NBasis, &
          ENuc, CAONO, XKin, TwoEl, Occ, NAc, NInAc, anSt, Flags)

          type(TTHCData), intent(out) :: THCData
          type(TACppData) :: AuxData
          integer, intent(in) :: NInte1, Ninte2, NBasis
          integer, intent(out) ::  NAc, NInAc
          double precision, intent(out) :: ENuc
          double precision, dimension(NInte1), intent(out) :: XKin
          double precision, dimension(NInte2), intent(out) :: TwoEl
          double precision, dimension(NBasis), intent(out) :: Occ
          double precision, dimension(NBasis, NBasis), intent(out) :: CAONO
          type(FlagsData), intent(in) :: Flags
          double precision, dimension(:,:), allocatable  :: CAONO_PYSCF
          integer, dimension(:), intent(in) :: anSt
          integer :: unit

          integer :: i, j, ab
          integer, external :: NAddr3 

          AuxData%nst = anSt
          print*, 'inside ORCA wrapper'
          print*, Flags%ITwoEl
          AuxData%NInte2 = NInte2
          call read_ORCA(THCData, AuxData, Flags, CAONO_PYSCF, TwoEl)
          print*, 'writing rdm2 full'
          call write_rdm2_dat(AuxData%rdm2_full, AuxData%NA)

          print*, 'zzz', AuxData%ECAS, AuxData%Enuc
          print*, 'po read'
          CAONO = CAONO_PYSCF
          ENuc = AuxData%ENuc

          Occ = AuxData%Occ
          NAc = AuxData%NA
          NInAc = AuxData%NI
!          print*, 'aha1', occ
          ! One-electron integrals are now transformed to NO
          do i = 1, Nbasis
                do j = 1, i
                      ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                      XKin(ab) = AuxData%HNO0(i,j)
                      ! if (i==j)then
                      !       print*, 'jjj', j, AuxData%HNO0(i,j)
                      ! end if
                end do
          end do
          
    end subroutine ORCA_wrapper

    subroutine load_cas_aux_data(suffix, AuxData, natural_orb)
          character(len=*), intent(in) :: suffix
          type(TACppData), intent(inout) :: AuxData
          integer, intent(out) :: natural_orb
          integer :: iv(1)
          double precision :: dv(1)

          ! REF/ holds what auxdata<suffix>.txt used to; the CAS energy is
          ! per state and lives on the state group as the ENERGY attribute.
          call h5_get(pyscf_fid, '/REF/NBASIS', iv); AuxData%Nbasis = iv(1)
          call h5_get(pyscf_fid, '/REF/NI',     iv); AuxData%NI     = iv(1)
          call h5_get(pyscf_fid, '/REF/NA',     iv); AuxData%NA     = iv(1)
          call h5_get(pyscf_fid, '/REF/NV',     iv); AuxData%NV     = iv(1)
          call h5_get(pyscf_fid, '/REF/NEL',    iv); AuxData%NEL    = iv(1)
          call h5_get(pyscf_fid, '/REF/NATORB', iv); natural_orb    = iv(1)
          call h5_get(pyscf_fid, '/REF/ENUC',   dv); AuxData%ENuc   = dv(1)
          call h5_attr(pyscf_fid, trim(pyscf_state_group(suffix)), 'ENERGY', AuxData%ECAS)

          AuxData%NIA = AuxData%NI + AuxData%NA

    end subroutine load_cas_aux_data


    subroutine load_2x2x2x2_array(filename, target_array)
          character(len=*), intent(in) :: filename
          double precision, dimension(:,:,:,:), intent(out) :: target_array
          logical :: exists
          integer :: unit

          unit = 33
          open(unit=unit, file=filename, status="old", access="stream", form="unformatted")
          read(unit) target_array
          close(unit)

    end subroutine load_2x2x2x2_array

    subroutine load_2x2_array_pyscf(filename, target_array)
          character(len=*), intent(in) :: filename
          double precision, dimension(:,:), intent(out) :: target_array
          integer :: unit

          unit = 33
          open(unit=unit, file=filename, status="old", access="stream", form="unformatted")
          read(unit) target_array
          close(unit)

    end subroutine load_2x2_array_pyscf
    

    subroutine load_2x2_array_orca(filename, target_array)
          character(len=*), intent(in) :: filename
          double precision, dimension(:,:), intent(out) :: target_array
          double precision, dimension(:,:), allocatable :: temp
          integer :: unit
          integer :: n, header(3)

          n = size(target_array, dim=1)
          allocate(temp(n,n))

          unit = 33

          open(unit=unit, file=filename, status="old", access="stream", form="unformatted")
          read(unit) header

          read(unit) temp
          target_array = transpose(temp)
          deallocate(temp)
                
          close(unit)

    end subroutine load_2x2_array_orca

    
    
    subroutine load_mp2_aux_data(AuxData)
          type(TACppData), intent(inout) :: AuxData
          integer :: iv(1)
          double precision :: dv(1)

          ! The ROHF export always writes a single state, tagged 1.1.
          call h5_get(pyscf_fid, '/REF/NBASIS', iv); AuxData%Nbasis = iv(1)
          call h5_get(pyscf_fid, '/REF/NI',     iv); AuxData%NI     = iv(1)
          call h5_get(pyscf_fid, '/REF/NA',     iv); AuxData%NA     = iv(1)
          call h5_get(pyscf_fid, '/REF/NV',     iv); AuxData%NV     = iv(1)
          call h5_get(pyscf_fid, '/REF/NEL',    iv); AuxData%NEL    = iv(1)
          call h5_get(pyscf_fid, '/REF/ENUC',   dv); AuxData%ENuc   = dv(1)
          call h5_attr(pyscf_fid, '/POSTHF/STATES/1.1', 'ENERGY', AuxData%EROHF)
          AuxData%NIA = AuxData%NI + AuxData%NA
          allocate(AuxData%IndAux(AuxData%NBasis))
          call h5_get(pyscf_fid, '/SCF/MO_OCC_INT', AuxData%IndAux)
    end subroutine load_mp2_aux_data

    subroutine load_density_matrices(suffix, AuxData, THCData, CAOMO, CAONO, Flags, AObasis, System, TwoEl, only_full)
          ! Loads RDMs from POSTHF/STATES/<state>/ in pyscf_data.h5:
          !   RDM2_AAAA -> rdm2_pp,  RDM2_ABAB -> rdm2_pm
          !   RDM2_BBBB -> rdm2_mm  (fallback: = rdm2_pp for closed-shell)
          !   RDM2_BABA -> rdm2_mp  (fallback: = rdm2_pm)
          !   RDM1_A -> rdm1_p,  RDM1_B -> rdm1_m
          !   RDM1   -> rdm1_full (fallback: rdm1_p = rdm1_m = rdm1_full/2)
          ! rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp  (computed here).
          
          character(len=*), intent(in) :: suffix
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          double precision, dimension(:,:), allocatable, intent(inout) :: CAOMO ! Changed to inout allocatable
          double precision, dimension(:,:), intent(inout) :: CAONO
          type(FlagsData), intent(in) :: Flags
          type(TAOBASIS), intent(in) :: AObasis
          type(TSystem), intent(in) :: System
          double precision, optional, intent(in) :: TwoEl(:)
          logical, optional, intent(in) :: only_full

          character(len=256) :: grp
          character(len=256) :: name_aaaa, name_abab, name_bbbb, name_baba, name_aa, name_bb
          character(len=256) :: name_rdm1, name_rdm2
          double precision, dimension(:), allocatable :: occ_temp
          integer :: unit
          integer :: i, j, k, l
          logical :: exists, ex_bbbb, ex_baba, ex_aa, ex_bb, ex_full
          logical :: use_full

          grp = pyscf_state_group(suffix)
          name_aa   = trim(grp)//'/RDM1_A'
          name_bb   = trim(grp)//'/RDM1_B'
          name_rdm1 = trim(grp)//'/RDM1'
          allocate(AuxData%rdm1_p(AuxData%NA, AuxData%NA))
          allocate(AuxData%rdm1_m(AuxData%NA, AuxData%NA))
          allocate(AuxData%rdm1_full(AuxData%NA, AuxData%NA))
          allocate(AuxData%rdm2_full(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))

          ! only_full asks for the pre-assembled rdm2_full_reordered file. Older dumps
          ! (produced before the pp branch existed) do not have it, but they do have the
          ! spin blocks - in that case fall back to them and assemble rdm2_full below.
          use_full = only_full
          if (only_full)then
             name_rdm2 = trim(grp)//'/RDM2_FULL_REORDERED'
             exists = h5_exists(pyscf_fid, trim(name_rdm2))
             if (exists)then
                   call h5_get(pyscf_fid, trim(name_rdm2), AuxData%rdm2_full)
             else
                   write(*, '(1x,3a)') 'Dataset ', trim(name_rdm2), &
                         ' not found - assembling rdm2_full from the spin blocks'
                   use_full = .false.
             end if
          end if

          if (.not. use_full)then
                name_aaaa = trim(grp)//'/RDM2_AAAA'
                name_abab = trim(grp)//'/RDM2_ABAB'
                name_bbbb = trim(grp)//'/RDM2_BBBB'
                name_baba = trim(grp)//'/RDM2_BABA'
                name_rdm2 = trim(grp)//'/RDM2'

                allocate(AuxData%rdm2_pp(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
                allocate(AuxData%rdm2_pm(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
                allocate(AuxData%rdm2_mm(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
                allocate(AuxData%rdm2_mp(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
                
                
                !call h5_get(pyscf_fid, trim(name_rdm2), AuxData%rdm2_full)
                call h5_get(pyscf_fid, trim(name_aaaa), AuxData%rdm2_pp)
                call h5_get(pyscf_fid, trim(name_abab), AuxData%rdm2_pm)

                ex_bbbb = h5_exists(pyscf_fid, trim(name_bbbb))
                if (ex_bbbb) then
                      call h5_get(pyscf_fid, trim(name_bbbb), AuxData%rdm2_mm)
                else
                      call h5_get(pyscf_fid, trim(name_aaaa), AuxData%rdm2_mm)
                end if

                ex_baba = h5_exists(pyscf_fid, trim(name_baba))
                if (ex_baba) then
                      print*, 'reading baba'
                      call h5_get(pyscf_fid, trim(name_baba), AuxData%rdm2_mp)
                else
                      print*, 'reading abab'
                      call h5_get(pyscf_fid, trim(name_abab), AuxData%rdm2_mp)
                end if
          end if

          ex_full = h5_exists(pyscf_fid, trim(name_rdm1))
          if(ex_full)then
                call h5_get(pyscf_fid, trim(name_rdm1), AuxData%rdm1_full)
          end if

          ex_aa = h5_exists(pyscf_fid, trim(name_aa))

          if (ex_aa) then
                call h5_get(pyscf_fid, trim(name_aa), AuxData%rdm1_p)
          else
                AuxData%rdm1_p = AuxData%rdm1_full / two
          end if
          
          ex_bb = h5_exists(pyscf_fid, trim(name_bb))
          if (ex_bb) then
                print*, 'reading bb'
                call h5_get(pyscf_fid, trim(name_bb), AuxData%rdm1_m)
          else
                AuxData%rdm1_m = AuxData%rdm1_full / two
          end if

          if (.not. use_full)then
                AuxData%rdm2_full = AuxData%rdm2_pp + AuxData%rdm2_mm +&
                      AuxData%rdm2_pm + AuxData%rdm2_mp
          end if

          ! print*, 'rdm1'
          ! do i = 1, AuxData%NA
          !       do j = 1, AuxData%NA
          !             write(*,'(2I5, 3F12.8)') i, j, AuxData%rdm1_p(i,j), AuxData%rdm1_m(i,j), AuxData%rdm1_full(i,j)
          !       end do
          ! end do
          ! AuxData%rdm2_pp = AuxData%rdm2_pp / two
          ! AuxData%rdm2_pm = AuxData%rdm2_pm / two
          ! AuxData%rdm2_mm = AuxData%rdm2_mm / two
          ! AuxData%rdm2_mp = AuxData%rdm2_mp / two
          ! AuxData%rdm2_full = AuxData%rdm2_full / two

          ! AuxData%rdm1_p = AuxData%rdm1_p / two
          ! AuxData%rdm1_m = AuxData%rdm1_m / two
          ! AuxData%rdm1_full = AuxData%rdm1_full / two


          ! Header

           ! associate(NA=>AuxData%NA)
!             do i = 1, NA
!                   do j = 1, NA
!                         write(*, '(A10, 2I5, 4F15.7)')'pluszus', i, j, &
!                               AuxData%rdm1_p(i, j), AuxData%rdm1_m(i, j), &
!                               AuxData%rdm1_p(i, j)+ AuxData%rdm1_m(i, j), &
!                               AuxData%rdm1_full(i, j)
!                   end do
!             end do
                        
          !   do i = 1, NA
          !         do j = 1, NA
          !               do k = 1, NA
          !                     do l = 1, NA
          !                           write(*, '(A10, 4I5, 6F15.7)')'snieznus', i, j, k, l, &
          !                                 AuxData%rdm2_pp(i, j, k, l), AuxData%rdm2_mm(i, j, k, l), &
          !                                 AuxData%rdm2_pm(i, j, k, l), AuxData%rdm2_mp(i, j, k, l), &
          !                                 AuxData%rdm2_pp(i, j, k, l)+ AuxData%rdm2_mm(i, j, k, l) +&
          !                                 AuxData%rdm2_pm(i, j, k, l)+ AuxData%rdm2_mp(i, j, k, l), &
          !                                 AuxData%rdm2_full(i, j, k, l)
          !                     end do
          !               end do
          !         end do
          !   end do
          ! end associate
! stop
!          call ab_occup(AuxData)
          
    end subroutine load_density_matrices

    subroutine load_rohf_occupation(AuxData)
          type(TACppData), intent(inout) :: AuxData
          integer :: unit
          logical :: exists

          allocate(AuxData%Occ_rohf(AuxData%NBasis))
          allocate(AuxData%Occ(AuxData%NBasis))

          call h5_get(pyscf_fid, '/SCF/MO_OCC', AuxData%Occ_rohf)

          AuxData%Occ = AuxData%Occ_rohf / two
    end subroutine load_rohf_occupation

    subroutine load_occupancy(AuxData)          
          type(TACppData), intent(inout) :: AuxData   
          double precision, allocatable :: occ_temp(:)
          integer :: unit
          logical :: exists                           

          ! POSTHF/OCC holds the CAS occupations of all NBasis orbitals; only the
          ! inactive+active head of it is used, exactly as when the same numbers
          ! were streamed out of rdm1.bin / occ.bin.
          allocate(occ_temp(h5_size(pyscf_fid, '/POSTHF/OCC')))
          call h5_get(pyscf_fid, '/POSTHF/OCC', occ_temp)
          ! print*, occ_temp, AuxData%NBasis

          allocate(AuxData%Occ(AuxData%NBasis))
          AuxData%Occ = 0.0d0
          AuxData%Occ(1:AuxData%NI+AuxData%NA) = occ_temp(1:AuxData%NI+AuxData%NA)

          allocate(AuxData%IndAux(AuxData%NBasis))
          AuxData%IndAux = 2
          AuxData%IndAux(1:AuxData%NI) = 0
          AuxData%IndAux(AuxData%NI+1:AuxData%NIA) = 1
!          print*, 'AuxData%Occ', AuxData%Occ

          call load_npm(AuxData)
    end subroutine load_occupancy

    subroutine load_npm(AuxData)
          ! Sets n_p, n_m, n from rdm1_p/m diagonal entries.
          ! Inactive: n_p(i)=n_m(i)=1. Active: from rdm1_p/m(i-NI,i-NI). Virtual: 0.
          ! n(i) = n_p(i) + n_m(i)  [NOT /2; this is what enters ssqrt convention].
          type(TACppData), intent(inout) :: AuxData
          integer :: i
          
          associate(NI=>AuxData%NI, NA=>AuxData%NA, NBasis=>AuxData%NBasis)

            if (.not. allocated(AuxData%n_p)) allocate(AuxData%n_p(NBasis))
            if (.not. allocated(AuxData%n_m)) allocate(AuxData%n_m(NBasis))
            if (.not. allocated(AuxData%n)) allocate(AuxData%n(NBasis))
            AuxData%n_p = zero
            AuxData%n_m = zero
            AuxData%n = zero

            do i =1, NI
                  AuxData%n_p(i) = one!frac12
                  AuxData%n_m(i) = one!frac12
                  AuxData%n(i) = AuxData%n_p(i)  + AuxData%n_m(i)
            end do

            do i = NI+1, NI+NA
                  AuxData%n_p(i) = AuxData%rdm1_p(i-NI, i-NI)
                  AuxData%n_m(i) = AuxData%rdm1_m(i-NI, i-NI)
                  AuxData%n(i) = AuxData%n_p(i)  + AuxData%n_m(i)
            end do
            !print*, 'AuxData%Occ', AuxData%Occ
            !print*, 'AuxData%n_p', AuxData%n_p
            !print*, 'AuxData%n_m', AuxData%n_m
            !print*, 'AuxData%n', AuxData%n
            
          end associate
          
    end subroutine load_npm

    subroutine load_one_electron_integrals(AuxData, THCData, CAOMO, CAONO, natural_orb, Flags)
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          double precision, dimension(:,:), allocatable, intent(out) :: CAOMO, CAONO
          integer, intent(in) :: natural_orb
          type(FlagsData), intent(in) :: Flags
          logical :: file_exists
          integer :: unit

          !----------------------------------------------------------------------------------------                                                                                                                                       
          ! READ 1-el integrals.                                                                                                                                                                                                          
          ! If natural = 1, HCore are in NO basis.                                                                                                                                                                                        
          ! If natural = 0, HCore are in AO basis.                                                                                                                                                                                        
          !----------------------------------------------------------------------------------------
          ! --- FCIDUMP modification: do not load HCore as it's already read from FCIDUMP ---
          call h5_get(pyscf_fid, '/INTS/CORE_HAMILTONIAN', AuxData%HNO0)
          ! -------------------------------------------------------------------------------------

          THCData%HNO = AuxData%HNO0

          ! Separate paths for CAOMO (not Natural) and CAONO (Natural)
          if (natural_orb==1) then

                allocate(CAONO(AuxData%NBasis, AuxData%NBasis))
                call h5_get(pyscf_fid, '/POSTHF/ORB_COEFF', CAONO)
                ! CAOMO is not needed/allocated here
          else
                ! If not natural, we read Canonical MOs into CAOMO for later transformation
                allocate(CAOMO(AuxData%NBasis, AuxData%NBasis))
                call h5_get(pyscf_fid, '/POSTHF/ORB_COEFF', CAOMO)
                ! Allocate CAONO to be filled later by Trans2NO (or copied if MP2)
                allocate(CAONO(AuxData%NBasis, AuxData%NBasis))

                if (Flags%JOBTYPE == JOB_TYPE_MP2 .or. Flags%JOBTYPE == JOB_TYPE_SRMP2) then
                      CAONO = CAOMO
                end if
          endif

    end subroutine load_one_electron_integrals

    subroutine resolve_state_suffix(AuxData, suffix)
          ! Picks the POSTHF/STATES/<n>.<m> group to read from. The requested
          ! state wins; without one (or when it is absent from the file) the
          ! single state present is used, and an ambiguous file is an error.
          type(TACppData), intent(in) :: AuxData
          character(len=*), intent(out) :: suffix
          character(len=256) :: grp
          integer :: x, y, count
          integer :: detected_state(2)

          suffix = ""

          if (AuxData%nst(1) > 0) then
                write(grp, '("/POSTHF/STATES/", I0, ".", I0)') AuxData%nst(1), AuxData%nst(2)
                if (h5_exists(pyscf_fid, trim(grp))) then
                      write(suffix, '("_", I0, ".", I0)') AuxData%nst(1), AuxData%nst(2)
                      return
                end if
                write(*, '(1x,3a)') 'Group ', trim(grp), ' not found in '//PYSCF_H5
          end if

          count = 0
          do x = 1, 10
                do y = 1, 10
                      write(grp, '("/POSTHF/STATES/", I0, ".", I0)') x, y
                      if (h5_exists(pyscf_fid, trim(grp))) then
                            detected_state(1) = x
                            detected_state(2) = y
                            count = count + 1
                      end if
                end do
          end do

          if (count == 0) then
                stop "Error: No POSTHF/STATES group found in pyscf_data.h5."
          end if

          if (count == 1) then
                write(suffix, '("_", I0, ".", I0)') detected_state(1), detected_state(2)
          else
                stop "Error: Several states in pyscf_data.h5 - select one in the input."
          end if
    end subroutine resolve_state_suffix
    

    subroutine read_ORCA(THCData, AuxData, Flags, CAONO, TwoEl)

          implicit none

          type(TTHCData), intent(inout) :: THCData
          type(TACppData), intent(inout) :: AuxData
          type(FlagsData), intent(in) :: Flags
          double precision, dimension(:,:), allocatable, intent(out) :: CAONO
          double precision, optional, intent(inout) :: TwoEl(:)
          double precision, dimension(:, :), allocatable :: CAOMO
          logical :: natural
          integer :: unit
          integer :: Nel
          double precision :: e_nucl
          double precision, dimension(:), allocatable :: occ_temp
          double precision, dimension(:, :), allocatable :: work
          type(TAOBASIS) :: AObasis
          type(TSystem)  :: System
          logical :: file_exists
          integer :: iostat
          double precision :: val, diff
          integer :: ms2
          integer :: i, j, k, l
          integer :: ninte
          logical :: spin_integrals


          spin_integrals = .false.
          
          AuxData%ORCA = 1

          print "(A)", "------------------------------------------------------------"
          print "(A)", " ORCA INTERFACE"
          if (Flags%ITwoEl.ne.1)then
                call GEOM_init(Flags, AuxData, AObasis, System)
          else
                if (spin_integrals)then
                      print "(A)", " Reading two-electron integrals from FCIDUMP spin file"
!                      call read_fcidump_header_spin("FCIDUMP", AuxData%NBasis, AuxData%Nel, ms2)
                else
                      print "(A)", " Reading two-electron integrals from FCIDUMP file"
                      call read_fcidump_header("FCIDUMP", AuxData%NBasis, AuxData%Nel, ms2)
                end if
          end if

!          call load_1rdm_orca_header('G1.bin', AuxData%NA)


          AuxData%NA=4
          print*, AuxData%Nel, 'AuxData%Nel'

          associate(NBasis=>AuxData%NBasis, NA=>AuxData%NA)

            ninte = NBasis*(NBasis+1)/2
            AuxData%Ninte2 = ninte * (ninte+1)/2
              
            allocate(AuxData%rdm2_pp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_pm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_full(NA, NA, NA, NA))
            allocate(AuxData%rdm1_p(NA, NA))
            allocate(AuxData%rdm1_m(NA, NA))
            allocate(AuxData%rdm1_full(AuxData%NA, AuxData%NA))

            
            AuxData%rdm2_pp = zero
            AuxData%rdm2_pm = zero
            AuxData%rdm2_mm = zero
            AuxData%rdm2_mp = zero
            AuxData%rdm1_full = zero
            AuxData%rdm1_p = zero
            AuxData%rdm1_m = zero


            call load_1rdm_orca('G1_uu.bin', AuxData%rdm1_p)
            call load_1rdm_orca('G1_dd.bin', AuxData%rdm1_m)
            inquire(file='G2.bin', exist=file_exists)
            if (file_exists)then
                  call load_1rdm_orca('G1.bin', AuxData%rdm1_full)
            else
                  AuxData%rdm1_full = AuxData%rdm1_p + AuxData%rdm1_m
                  call save_1rdm_orca('G1.bin', AuxData%rdm1_full)
            end if
!            AuxData%rdm1_p = AuxData%rdm1_full/two
 !           AuxData%rdm1_m = AuxData%rdm1_full/two


            call check_natural(AuxData, natural)

            print*, 'calc_occupancy_orca'
            call calc_occupancy_orca(AuxData)
            print*, 'done_calc_occupancy_orca'


            call load_2rdm_orca('G2_uu.bin', AuxData%rdm2_pp, .false.)
            call load_2rdm_orca('G2_dd.bin', AuxData%rdm2_mm, .false.)
            call load_2rdm_orca('G2_ud.bin', AuxData%rdm2_pm, .false.)

            AuxData%rdm2_mp = AuxData%rdm2_pm
            inquire(file='G2.bin', exist=file_exists)
            if (file_exists)then
                  call load_2rdm_orca('G2.bin', AuxData%rdm2_full, .true.)
            else
                  AuxData%rdm2_full = (AuxData%rdm2_pp +AuxData%rdm2_mm +&
                        AuxData%rdm2_pm +AuxData%rdm2_mp )
                  call save_2rdm_orca('G2.bin', AuxData%rdm2_full, .true.)
            end if

            ! AuxData%rdm2_pp = AuxData%rdm2_pp / two
            ! AuxData%rdm2_pm = AuxData%rdm2_pm / two
            ! AuxData%rdm2_mm = AuxData%rdm2_mm / two
            ! AuxData%rdm2_mp = AuxData%rdm2_mp / two
            ! AuxData%rdm2_full = AuxData%rdm2_full / two

            ! AuxData%rdm1_p = AuxData%rdm1_p / two
            ! AuxData%rdm1_m = AuxData%rdm1_m / two
            ! AuxData%rdm1_full = AuxData%rdm1_full / two
            



            !write(*, '("p", T6, "q", T11, "r", T16, "s", T25, "rdm_pp", T40, "rdm_mm", T55, "rdm_pm", T70, "rdm_mp", T85, "rdm_full")')
            print*, 'rdmyyy przed natural'
            associate(NA=>AuxData%NA)
            do i = 1, NA
              do j = 1, NA
                do k = 1, NA
                      do l = 1, NA
                           write(*, '(4I5, 5F15.8)') i, j, k, l, &
                                 AuxData%rdm2_full(i,j,k,l),&
                                  AuxData%rdm2_pp(i,j,k,l), &
                                  AuxData%rdm2_mm(i,j,k,l), &
                                AuxData%rdm2_pm(i,j,k,l), &
                                AuxData%rdm2_mp(i,j,k,l)
                            

                  end do
                end do
              end do
            end do
            end associate
            

            !call print_rdm_orca(AuxData)

            call read_fcidump_restricted('FCIDUMP', AuxData, TwoEl)
            

            !            print*, 'check energy before transforming'
            
            !            call check_energy_incore(AuxData, TwoEl, spin_sep=.true.)          
            if (.not. natural)then
                  allocate(CAONO(AuxData%NBasis, AuxData%NBasis))
                  allocate(CAOMO(AuxData%NBasis, AuxData%NBasis))
                  
                  if (Flags%ITwoEl .ne. 1)then
                        ! Out-core/THC: Read C.bin (AO->MO) and call Trans2NO_AO
                        inquire(file='C.bin', exist=file_exists)
                        if (file_exists) then
                              call load_2x2_array_orca('C.bin', CAOMO)
                        end if
                        print*, 'call transfrom to no'
                        call Trans2NO_AO(THCData, AuxData, &
                              CAOMO, CAONO, Flags, TwoEl=TwoEl)
                        print*, 'calc_occupancy_orca2'
                        call calc_occupancy_orca(AuxData)
                        print*, 'done_calc_occupancy_orca2'
                        
                  else
                        ! In-core: MO Integrals, C.bin not needed, call Trans2NO_MO
                        ! CAOMO is passed but not used/read
                        call Trans2NO_MO(THCData, AuxData, &
                              CAONO, Flags, TwoEl=TwoEl)
                        
                        print*, 'calc_occupancy_orca2'
                        call calc_occupancy_orca(AuxData)
                        print*, 'done_calc_occupancy_orca2'
                        
                  end if
                  allocate(AuxData%CAONO(NBasis, NBasis))
                  AuxData%CAONO = CAONO

                  print*, 'rdmyyy po natural'
                  associate(NA=>AuxData%NA)
                    do i = 1, NA
                          do j = 1, NA
                                do k = 1, NA
                                      do l = 1, NA
                                            write(*, '(4I5, 5F15.8)') i, j, k, l, &
                                                  AuxData%rdm2_full(i,j,k,l),&
                                            AuxData%rdm2_pp(i,j,k,l), &
                                                  AuxData%rdm2_mm(i,j,k,l), &
                                                  AuxData%rdm2_pm(i,j,k,l), &
                                                  AuxData%rdm2_mp(i,j,k,l)
                                            
                                            
                                      end do
                                end do
                          end do
                    end do

                    print*, '1rdm'
                    do i = 1, NA
                          do j = 1, NA
                                write(*, '(2I5, 3F15.8)') i, j, AuxData%rdm1_p(i,j), AuxData%rdm1_m(i,j), AuxData%rdm1_full(i,j)
                          end do
                    end do

                  end associate

            else
                  allocate(CAONO(NBasis, NBasis))
                  CAONO = zero
                  do i = 1, NBASIS
                        CAONO(i,i) = one
                  end do
            end if

                      
            AuxData%spinsep = .true.
            do i = 1, NA
                  if (abs(AuxData%rdm1_p(i,i)-AuxData%rdm1_m(i,i)).gt.1.d-8)then
                        AuxData%spinsep = .false.
                  end if
            end do

            open(10, file='ure_casno.dat')
            write(10, *) CAONO
            close(10)

            print*, 'check energy after'
            call check_energy_incore(AuxData, TwoEl, spin_sep=.true.)          

            if (.not. natural) then
                  print*, 'writing rdm2 full'
                  call write_rdm2_dat(AuxData%rdm2_full, NA)
            end if

          end associate
          
    end subroutine read_ORCA


    subroutine read_ORCA_spinres(THCData, AuxData, Flags, Ints)
          ! ORCA interface, spin-resolved. Reads G1/G2 RDMs, then FCIDUMP.
          ! DUCC path (ducc=.true.):  ints2e_aa has factor 2; ints2e_bb not allocated.
          ! Non-DUCC (ducc=.false.):  separate aa and bb blocks, ints2e_aa has factor 2.
          ! Transforms MO->NO (Trans2NO_MO_spinres) if rdm1_full is not diagonal.

          implicit none

          type(TTHCData), intent(inout) :: THCData
          type(TACppData), intent(inout) :: AuxData
          type(FlagsData), intent(in) :: Flags
          type(TInts), intent(out) :: Ints          
          integer :: unit
          integer :: Nel
          double precision :: e_nucl
          double precision, dimension(:), allocatable :: occ_temp

          type(TAOBASIS) :: AObasis
          type(TSystem)  :: System

          integer :: iostat
          double precision :: val
          logical :: IsDiag, file_exists
          integer :: ms2
          integer :: i, j, k, l
          
          AuxData%ORCA = 1
          print "(A)", "------------------------------------------------------------"
          print "(A)", " ORCA INTERFACE SPIN RESOLVED"

          if (Flags%ITwoEl.ne.1)then
                call GEOM_init(Flags, AuxData, AObasis, System)
          else
                print "(A)", " Reading two-electron integrals from FCIDUMP file"
                call read_fcidump_header("FCIDUMP", AuxData%NBasis, AuxData%Nel, ms2)
          end if

          call load_1rdm_orca_header('G1.bin', AuxData%NA)

          associate(NBasis=>AuxData%NBasis, NA=>AuxData%NA)

              
            allocate(AuxData%rdm2_pp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_pm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_full(NA, NA, NA, NA))
            allocate(AuxData%rdm1_p(NA, NA))
            allocate(AuxData%rdm1_m(NA, NA))
            allocate(AuxData%rdm1_full(NA, NA))

            AuxData%rdm2_pp = zero; AuxData%rdm2_pm = zero; AuxData%rdm2_mm = zero; AuxData%rdm2_mp = zero
            AuxData%rdm1_full = zero; AuxData%rdm2_full = zero
            AuxData%rdm1_p = zero; AuxData%rdm1_m = zero

            call load_1rdm_orca('G1.bin', AuxData%rdm1_full)

            call load_1rdm_orca('G1_uu.bin', AuxData%rdm1_p)

            call load_1rdm_orca('G1_dd.bin', AuxData%rdm1_m)

            call calc_occupancy_orca(AuxData)

            call load_2rdm_orca('G2_uu.bin', AuxData%rdm2_pp, .false.)
            call load_2rdm_orca('G2_dd.bin', AuxData%rdm2_mm, .false.)
            call load_2rdm_orca('G2_ud.bin', AuxData%rdm2_pm, .false.)
            AuxData%rdm2_mp = AuxData%rdm2_pm
            inquire(file='G2.bin', exist=file_exists)
            if (file_exists)then
                  call load_2rdm_orca('G2.bin', AuxData%rdm2_full, .true.)
            else
                  AuxData%rdm2_full = (AuxData%rdm2_pp +AuxData%rdm2_mm +&
                        AuxData%rdm2_pm +AuxData%rdm2_mp )
            end if


            !call print_rdm_orca(AuxData)

            Ints%ints2e_dim = int(NBasis,8)**2 * (int(NBasis,8)**2 + 1) / 2
            if (Flags%JOBTYPE == JOB_TYPE_DUCC)then
                  ! ints_aa=ints_bb
                  call read_fcidump_unrestricted('FCIDUMP', AuxData, Ints, NBasis, .true., .true.)
            else
                  ! ints_aa/=ints_bb
                  print*, 'yes, allokuje'
                  if (.not. allocated(Ints%ints2e)) then
                        print *, 'nienieniekurwie'
                  else
                        print*, 'takokurwie'
                  end if
                        
                  call read_fcidump_unrestricted('FCIDUMP', AuxData, Ints, NBasis, .false., .true.)
            end if

            allocate(AuxData%HNO0(NBasis, NBasis))
            AuxData%HNO0 = Ints%ints1e_aa

            call check_natural(AuxData, IsDiag)
            call check_energy_incore_spinres(AuxData, Ints)
            if (.not. IsDiag) then
                  print*, 'Transforming from MO to NO basis'
                  print*, 'kuper1', AuxData%n_p
                  call Trans2NO_MO_spinres(AuxData, Flags, Ints)
                  print*, 'kuper2', AuxData%n_p
                  print*, 'calc_occupancy_orca2'
                  call calc_occupancy_orca(AuxData)
                  print*, 'done_calc_occupancy_orca2'

            end if

            call check_energy_incore_spinres(AuxData, Ints)
            

          end associate

    end subroutine read_ORCA_spinres

    subroutine read_ORCA_ducc(THCData, AuxData, Flags, Ints_full, Ints_ducc)

          implicit none

          type(TTHCData), intent(inout) :: THCData
          type(TACppData), intent(inout) :: AuxData
          type(FlagsData), intent(in) :: Flags
          type(TInts), intent(out) :: Ints_full
          type(TInts), intent(out) :: Ints_ducc          
          integer :: unit
          integer :: Nel
          double precision :: e_nucl
          double precision, dimension(:), allocatable :: occ_temp

          type(TAOBASIS) :: AObasis
          type(TSystem)  :: System

          integer :: iostat
          double precision :: val
          logical :: natural
          logical :: IsDiag
          integer :: ms2
          integer :: i, j, k, l
          logical :: file_exists
          
          AuxData%ORCA = 1
          print "(A)", "------------------------------------------------------------"
          print "(A)", " ORCA INTERFACE SPIN RESOLVED ducc"

          if (Flags%ITwoEl.ne.1)then
                call GEOM_init(Flags, AuxData, AObasis, System)
          else
                print "(A)", " Reading two-electron integrals from FCIDUMP-full file"
                call read_fcidump_header("FCIDUMP-full", AuxData%NBasis, AuxData%Nel, ms2)
          end if
          print*, 'loading g1'
          call load_1rdm_orca_header('ref/G1.bin', AuxData%NA)

          associate(NBasis=>AuxData%NBasis, NA=>AuxData%NA)
              
            allocate(AuxData%rdm2_pp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_pm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mm(NA, NA, NA, NA))
            allocate(AuxData%rdm2_mp(NA, NA, NA, NA))
            allocate(AuxData%rdm2_full(NA, NA, NA, NA))
            allocate(AuxData%rdm1_p(NA, NA))
            allocate(AuxData%rdm1_m(NA, NA))
            allocate(AuxData%rdm1_full(NA, NA))

            AuxData%rdm2_pp = zero; AuxData%rdm2_pm = zero; AuxData%rdm2_mm = zero; AuxData%rdm2_mp = zero
            AuxData%rdm1_full = zero; AuxData%rdm2_full = zero
            AuxData%rdm1_p = zero; AuxData%rdm1_m = zero
            print*, 'loading g1_uu'

            call load_1rdm_orca('ref/G1.bin', AuxData%rdm1_full)
            call load_1rdm_orca('ref/G1_uu.bin', AuxData%rdm1_p)
            call load_1rdm_orca('ref/G1_dd.bin', AuxData%rdm1_m)
            print*, 'call occupancy'
            call calc_occupancy_orca(AuxData)
            print*, 'load g2'

            print*, 'inside read_orca_ducc ducc'
            do i = 1, Nbasis
                  write(*, '(A5, I3, 3F12.6)') 'nnn', i, AuxData%occ(i), AuxData%n_p(i), AuxData%n_m(i)
            end do
            print*, ''
            do i = 1, NA
                do j = 1, NA
                      write(*, '(A5, 2I3, 2F12.6)') 'rdm', i, j, AuxData%rdm1_p(i, j), AuxData%rdm1_m(i, j)
                end do
          end do

            call load_2rdm_orca('ref/G2.bin', AuxData%rdm2_full, .true.)
            print*, 'load g2uu'
            call load_2rdm_orca('ref/G2_uu.bin', AuxData%rdm2_pp, .false.)
            print*, 'load g2dd'
            call load_2rdm_orca('ref/G2_dd.bin', AuxData%rdm2_mm, .false.)
            print*, 'load g2ud'
            call load_2rdm_orca('ref/G2_ud.bin', AuxData%rdm2_pm, .false.)
            AuxData%rdm2_mp = AuxData%rdm2_pm

            !call print_rdm_orca(AuxData)

            Ints_full%ints2e_dim = int(NBasis,8)**2 * (int(NBasis,8)**2 + 1) / 2
            Ints_ducc%ints2e_dim = int(NA,8)**2 * (int(NA,8)**2 + 1) / 2
            print*, 'Ints_full%ints2e_dim', Ints_full%ints2e_dim
            print*, 'Ints_ducc%ints2e_dim', Ints_ducc%ints2e_dim
            ! ints_aa=ints_bb
            print*, 'load_unrest full'

            inquire(file="FCIDUMP-full", exist=file_exists)
            if (.not. file_exists) then
                  stop "Error: File FCIDUMP-full not found."
            end if
            call read_fcidump_unrestricted('FCIDUMP-full', AuxData, Ints_full, NBasis, .true., .true.)
            print*, 'load_unrest ducc'
            inquire(file="FCIDUMP-ducc", exist=file_exists)
            if (.not. file_exists) then
                  stop "Error: File FCIDUMP-ducc not found."
            end if
            call read_fcidump_unrestricted('FCIDUMP-ducc', AuxData, Ints_ducc, NA, .true., .false.)
            print*, 'check energy'
            call check_energy_ducc_reference(AuxData, Ints_full)
            call check_natural(AuxData, natural)
            print*, 'natural', natural

           if (.not. natural)then
!                 allocate(CAONO(AuxData%NBasis, AuxData%NBasis))
!                 allocate(CAOMO(AuxData%NBasis, AuxData%NBasis))

                 ! call Trans2NO_MO_ducc(AuxData, &
                 !       Flags, Ints_full, Ints_ducc)
                 ! print*, 'calc_occupancy_orca2'
                 ! call calc_occupancy_orca(AuxData)
                 ! print*, 'done_calc_occupancy_orca2'

           end if
            call check_energy_ducc_reference(AuxData, Ints_full)



            !call check_energy_incore_spinres(AuxData, Ints, spin_sep=.true.)


          end associate

    end subroutine read_ORCA_ducc

    subroutine check_energy_ducc_reference(AuxData, Ints)
          ! Verifies reference energy using full integrals and DMRG RDMs.
          ! E = ENuc + E1 + E2 where:
          !   E1 = 2 * sum_{pq in AA} ints1e_aa(p,q) * rdm1_p(p-NI, q-NI)
          !   E2_aa = sum_{pqrs in AAAA} ints2e_aa(pqrs) * rdm2_pp(p-NI,r-NI,q-NI,s-NI)
          !   E2_ab = sum_{pqrs in AAAA} ints2e_ab(pqrs) * rdm2_pm(p-NI,r-NI,q-NI,s-NI)
          !   E2 = 0.5*E2_aa + E2_ab   [frac12 from factor-of-2 in ints2e_aa storage]
          type(TACppData), intent(inout) :: AuxData
          type(TInts), intent(in) :: Ints

          integer :: i, j, k, l
          integer :: p, q, r, s
          double precision :: e1_ref, e2_aa_ref, e2_ab_ref, e2_ref, total_energy
          integer(I8) :: idx
          integer :: NI, NA, NBasis

          associate(NBasis=>AuxData%NBasis, NI=>AuxData%NI, NA=>AuxData%NA)
          
          e1_ref = zero
          do i = 1, NA
              do j = 1, NA
                  p = NI + i
                  q = NI + j
                  e1_ref = e1_ref + two * Ints%ints1e_aa(p,q) * AuxData%rdm1_p(i,j)
                  if (abs(AuxData%rdm1_p(i,j)).gt.1.d-5)then
                        write(*, '(A5, 2I5, F12.8)')'sniez', i, j, AuxData%rdm1_p(i,j)
                  end if
              end do
          end do

          e2_aa_ref = zero
          e2_ab_ref = zero

          do l = 1, NA
              do k = 1, NA
                  do j = 1, NA
                      do i = 1, NA
                          p = NI + i
                          q = NI + j
                          r = NI + k
                          s = NI + l
                          
                          ! AA part
                          idx = gmap_4fold(p,q,r,s, NBasis)
                          e2_aa_ref = e2_aa_ref + Ints%ints2e_aa(idx) * AuxData%rdm2_pp(i,k,j,l)
                          
                          ! AB part
                          e2_ab_ref = e2_ab_ref + Ints%ints2e_ab(idx) * AuxData%rdm2_pm(i,k,j,l)

                      end do
                  end do
              end do
          end do

          e2_ref = 0.5d0 * e2_aa_ref + e2_ab_ref
          total_energy = AuxData%ENuc + e1_ref + e2_ref

          print*, 'Calculate reference energy (DUCC)'
          write(*, '(A, F20.12)') 'e1_ref', e1_ref
          write(*, '(A, F20.12)') 'e2_aa_ref', e2_aa_ref
          write(*, '(A, F20.12)') 'e2_ab_ref', e2_ab_ref
          write(*, '(A, F20.12)') 'Reference energy', total_energy
          AuxData%E_ref_ducc = total_energy

          end associate

    end subroutine check_energy_ducc_reference

    subroutine print_rdm_orca(AuxData)          
          type(TACppData), intent(inout) :: AuxData
          integer :: i, j, k, l

          associate(NA=>AuxData%NA)


            do i =1, NA
                  do j =1, NA
                        if (abs(AuxData%rdm1_full(i,j)).gt.1.d-5)then
                              print*, 'rdm1full', i, j, AuxData%rdm1_full(i,j)
                        end if
                  end do
            end do

            print*, ''
            do i =1, NA
                  do j =1, NA
                        if (abs(AuxData%rdm1_p(i,j)).gt.1.d-5)then
                              print*, 'rdm1p', i, j, AuxData%rdm1_p(i,j)
                        end if
                  end do
            end do

            print*, ''
            do i =1, NA
                  do j =1, NA
                        if (abs(AuxData%rdm1_m(i,j)).gt.1.d-5)then
                              print*, 'rdm1m', i, j, AuxData%rdm1_m(i,j)
                        end if
                  end do
            end do


            print*, "--- RDM2 FULL ---"
            do i = 1, NA
                  do j = 1, NA
                        do k = 1, NA
                              do l = 1, NA
                                    if (abs(AuxData%rdm2_full(i,j,k,l)) > 1.d-5) then
                                          print*, "rdm2full", i, j, k, l, AuxData%rdm2_full(i,j,k,l)
                                    end if
                              end do
                        end do
                  end do
            end do

            print*, "--- RDM2 PP ---"
            do i = 1, NA
                  do j = 1, NA
                        do k = 1, NA
                              do l = 1, NA
                                    if (abs(AuxData%rdm2_pp(i,j,k,l)) > 1.d-5) then
                                          print*, "rdm2pp", i, j, k, l, AuxData%rdm2_pp(i,j,k,l)
                                    end if
                              end do
                        end do
                  end do
            end do

            print*, "--- RDM2 MM ---"
            do i = 1, NA
                  do j = 1, NA
                        do k = 1, NA
                              do l = 1, NA
                                    if (abs(AuxData%rdm2_mm(i,j,k,l)) > 1.d-5) then
                                          print*, "rdm2mm", i, j, k, l, AuxData%rdm2_mm(i,j,k,l)
                                    end if
                              end do
                        end do
                  end do
            end do

            print*, "--- RDM2 PM ---"
            do i = 1, NA
                  do j = 1, NA
                        do k = 1, NA
                              do l = 1, NA
                                    if (abs(AuxData%rdm2_pm(i,j,k,l)) > 1.d-5) then
                                          print*, "rdm2pm", i, j, k, l, AuxData%rdm2_pm(i,j,k,l)
                                    end if
                              end do
                        end do
                  end do
            end do

          end associate
    end subroutine print_rdm_orca
    
    subroutine load_S_packed(fname, S, NBasis)
          implicit none
          character(len=*), intent(in) :: fname
          integer, intent(in) :: NBasis
          double precision, intent(out) :: S(NBasis,NBasis)

          integer :: iu, i, j
          double precision :: x
          integer :: header(3)

          S = 0.0d0

          open(newunit=iu, file=fname, form="unformatted", access="stream", status="old", action="read")
          
          ! Skip 12-byte header
          read(iu) header
          
          do i = 1, NBasis
             do j = 1, i
                   read(iu) x
                   if (abs(x).gt.1.e-5)then
                         print*, i, j, x
                   end if
                S(i,j) = x
                S(j,i) = x
             end do
          end do

          close(iu)
    end subroutine load_S_packed


    subroutine calc_occupancy_orca(AuxData)          
          type(TACppData), intent(inout) :: AuxData
          integer :: p
          double precision :: sum_act

          sum_act = zero

          do p = 1, AuxData%NA
                sum_act = sum_act + AuxData%rdm1_p(p,p) + AuxData%rdm1_m(p,p)                
          end do

          print*, 'sum_act', sum_act
          AuxData%NI = (AuxData%NEL - int(sum_act))/2
          AuxData%NIA = AuxData%NI + AuxData%NA
          AuxData%NV = AuxData%NBasis - AuxData%NIA

          print*,''
          write(*, '(A20, I10)') 'Number of electrons', AuxData%NEL
          write(*, '(A20, I10)') 'Number of inactive', AuxData%NI
          write(*, '(A20, I10)') 'Number of acvitve', AuxData%NA
          write(*, '(A20, I10)') 'Number of occupied', AuxData%NIA
          write(*, '(A20, I10)') 'Number of virtual', AuxData%NV
          write(*, '(A20, I10)') 'Number of basis func', AuxData%NBasis
          print*, ''
          


          if (.not. allocated(AuxData%Occ)) allocate(AuxData%Occ(AuxData%NBasis))
          AuxData%Occ = zero
          
          AuxData%Occ(1:AuxData%NI) = one
          ! print*, 'AuxData%Occ(1:AuxData%NI)', AuxData%Occ(1:AuxData%NI)
          
          do p = 1, AuxData%NA
                AuxData%Occ(AuxData%NI + p) =  (AuxData%rdm1_p(p,p) + AuxData%rdm1_m(p,p))/two
          end do
         

          if (.not. allocated(AuxData%IndAux)) allocate(AuxData%IndAux(AuxData%NBasis))
          AuxData%IndAux = 2
          AuxData%IndAux(1:AuxData%NI) = 0
          AuxData%IndAux(AuxData%NI+1:AuxData%NIA) = 1
          call load_npm(AuxData)
          print*, ''
          print*, 'occupancy loaded'
    end subroutine calc_occupancy_orca
    
    
    subroutine read_fcidump_header(filename, norb, nelec, ms2)
            character(len=*), intent(in) :: filename
            integer, intent(out) :: norb, nelec, ms2
            integer :: unit, ios, tblock
            character(len=256) :: line
            integer ::  pos, i

            tblock = 1
            norb = -1
            nelec = -1
            ms2 = -1
            unit = 20
            open(unit=unit, file=filename, status="old", action="read", iostat=ios)

            header_loop: do i = 1, 5
                  read(unit=unit, fmt="(a)", iostat=ios) line
                  line = trim(adjustl(line))

                  pos = index(line, "NORB=")
                  if (pos > 0) read(line(pos+5:), *) norb

                  pos = index(line, "NELEC=")
                  if (pos > 0) read(line(pos+6:), *) nelec

                  pos = index(line, "MS2=")
                  if (pos > 0) read(line(pos+4:), *) ms2
            end do header_loop

            close(unit)

      end subroutine read_fcidump_header

      subroutine read_fcidump_restricted(filename, AuxData, TwoEl)
            ! Reads a standard spin-free FCIDUMP (one block of 2e integrals).
            ! FORMAT: single block of 2e integrals (no spin separation), then 1e integrals,
            ! 8-fold: (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr)=(rs|pq)=(sr|pq)=(rs|qp)=(sr|qp)
            ! TwoEl(gmap(p,q,r,s)) = (pq|rs)  [physical, no prefactor]
            ! HNO0(p,q) symmetrized. ENuc read from all-zero index line.


            character(len=*), intent(in) :: filename
            type(TACppData), intent(inout) :: AuxData
            real(F64), dimension(:), intent(out) :: TwoEl

            integer :: unit, ios
            character(len=256) :: line
            integer :: p, q, r, s
            integer(I8) :: idx
            real(F64) :: val
            logical :: in_header


            unit = 10
            allocate(AuxData%HNO0(AuxData%NBasis, AuxData%NBasis))
            AuxData%HNO0 = zero
            TwoEl = zero

            open(unit=unit, file=filename, status="old", action="read", iostat=ios)

            in_header = .true.
            do while (in_header)
                  read(unit, '(a)', iostat=ios) line
                  if (index(line, '/') > 0 .or. index(line, '&END') > 0) in_header = .false.
            end do

            read_loop: do
                  read(unit, '(a)', iostat=ios) line
                  line = adjustl(line)
                  if (len_trim(line) == 0) cycle

                  read(line, *) val, p, q, r, s

                  if (p == 0 .and. q == 0 .and. r == 0 .and. s == 0) then
                        AuxData%ENuc = val
                        exit
                  else if (r == 0 .and. s == 0) then
                        AuxData%HNO0(p, q) = val
                        AuxData%HNO0(q, p) = val
                  else
                        idx = gmap(p, q, r, s)
                        TwoEl(idx) = val
                  end if
            end do read_loop

            close(unit)
      end subroutine read_fcidump_restricted

      subroutine read_fcidump_unrestricted(filename, AuxData, Ints, NDim, ducc, rnucl)
            ! Reads UHF-split FCIDUMP (5 blocks: aa, bb, ab, 1e-aa, 1e-bb).
            ! 4-fold: (pq|rs)=(qp|sr)=(rs|pq)=(sr|qp)   [(pq|rs) != (pq|sr)]
            ! NDim = NBasis for full integrals, NA for DUCC active-space-only.
            ! ducc=.true.  -> skip bb blocks; ints2e_bb not allocated; ints2e_bb assumed=ints2e_aa
            ! ducc=.false. -> read all 5 blocks
            ! rnucl=.true. -> read ENuc from terminal non-zero entry
            ! FULL integrals:  ints2e_aa(p,q,r,s) = 2*(pq|rs)_phys  [factor 2 in aa/bb]
            !                  ints2e_ab(p,q,r,s) =   (pq|rs)_phys  [no factor in ab]
            ! DUCC integrals:  ints2e_aa(p,q,r,s) = 2*(pq|rs)_phys  [factor 2 in aa/bb]
            !                  ints2e_ab(p,q,r,s) =   (pq|rs)_phys  [no factor in ab]

            character(len=*), intent(in) :: filename
            type(TACppData), intent(inout) :: AuxData
            type(TInts), intent(inout) :: Ints
            integer, intent(in) :: Ndim
            logical, intent(in) :: ducc, rnucl

            integer :: unit, ios, tblock, i
            integer :: p, q, r, s
            integer(I8) :: idx
            real(F64) :: val
            character(len=256) :: line
            logical :: in_header

            associate(ints2e_dim=>Ints%ints2e_dim)

              if (rnucl) then
                    AuxData%ENuc = zero
              end if
            
            allocate(Ints%ints1e_aa(NDim, NDim))
            allocate(Ints%ints2e_aa(ints2e_dim))
            allocate(Ints%ints2e_ab(ints2e_dim))
            print*, 'zaalokowane teraz'
            if (.not. allocated(Ints%ints2e)) then
                  print *, 'ppp-nie'
            else
                  print *, 'pooo-tak'
            end if
            Ints%ints1e_aa = zero
            Ints%ints2e_aa = zero
            Ints%ints2e_ab = zero
            
            if (.not. ducc) then
                  allocate(Ints%ints1e_bb(NDim, NDim))
                  allocate(Ints%ints2e_bb(ints2e_dim))
                  Ints%ints1e_bb = zero
                  Ints%ints1e_aa = zero
            end if

            unit = 10
            tblock = 1 
            open(unit=unit, file=filename, status="old", action="read", iostat=ios)

            in_header = .true.
            do while (in_header)
                  read(unit, '(a)', iostat=ios) line
                  if (index(line, '/') > 0 .or. index(line, '&END') > 0) in_header = .false.
            end do

            read_loop: do
                  read(unit=unit, fmt="(a)", iostat=ios) line
                  if (ios /= 0) exit read_loop

                  line = trim(adjustl(line))
                  if (len_trim(line) == 0) cycle read_loop

                  read(line, *, iostat=ios) val, p, q, r, s
                  if (ios /= 0) cycle read_loop

                  if (p == 0 .and. q == 0 .and. r == 0 .and. s == 0) then
                        if (abs(val) < 1.0d-12) then
                              tblock = tblock + 1
                              cycle read_loop
                        else
                              if (rnucl) then
                                    AuxData%ENuc = val
                              end if
                              exit read_loop
                        end if
                  end if

                  select case (tblock)
                  case (1) 
                        idx = gmap_4fold(p,q,r,s,NDim)
                        Ints%ints2e_aa(idx) = val
                  case (2)
                        if (.not. ducc)then
                              idx = gmap_4fold(p,q,r,s,NDim)
                              Ints%ints2e_bb(idx) = val
                        end if
                  case (3)
                        idx = gmap_4fold(p,q,r,s,NDim)
                        Ints%ints2e_ab(idx) = val
                  case (4) 
                        if (r == 0 .and. s == 0) then
                              Ints%ints1e_aa(p, q) = val
                              Ints%ints1e_aa(q, p) = val
                        end if
                  case (5) 
                        if (r == 0 .and. s == 0) then
                              if (.not. ducc)then
                                    Ints%ints1e_bb(p, q) = val
                                    Ints%ints1e_bb(q, p) = val
                              end if
                        end if
                  end select

            end do read_loop

            close(unit)

            end associate

      end subroutine read_fcidump_unrestricted


      subroutine load_1rdm_orca(rdm_file, rdm)

            character(len=*), intent(in) :: rdm_file
            real(F64), dimension(:,:), intent(out) :: rdm

            integer :: dim
            integer :: unit, i, j, k
            integer :: header(3)
            real(F64), dimension(:), allocatable :: tril_data
            integer :: num_elements
            integer :: n

            n = size(rdm,dim=1)

            open(newunit=unit, file=rdm_file, form='unformatted', &
                  access='stream', status='old')

            read(unit) header
            print*, 'rdm_file', rdm_file
            print*, header
            dim = header(1)
            if (n.ne.dim)then
                  print*, 'fcidump and rdm1 files incompatible, exiting'
                  print*, 'n in fcidump', n
                  print*, 'n in rdm', dim
                  stop
            end if

            num_elements = dim * (dim + 1) / 2

            allocate(tril_data(num_elements))
            read(unit) tril_data
            close(unit)

            k = 1
            do i = 1, dim
                  do j = 1, i
                        rdm(i, j) = tril_data(k)
                        rdm(j, i) = tril_data(k)
                        k = k + 1
                  end do
            end do
      end subroutine load_1rdm_orca

      subroutine save_1rdm_orca(rdm_file, rdm)

      character(len=*), intent(in) :: rdm_file
      real(F64), dimension(:,:), intent(in) :: rdm

      integer :: dim
      integer :: unit, i, j, k
      integer :: header(3)
      real(F64), dimension(:), allocatable :: tril_data
      integer :: num_elements

      dim = size(rdm, dim=1)

      header(1) = dim
      header(2) = dim
      header(3) = dim*2

      num_elements = dim * (dim + 1) / 2
      allocate(tril_data(num_elements))

      k = 1
      do i = 1, dim
            do j = 1, i
                  tril_data(k) = rdm(i, j)
                  k = k + 1
            end do
      end do

      open(newunit=unit, file=rdm_file, form='unformatted', &
            access='stream', status='replace')

      write(unit) header
      write(unit) tril_data

      close(unit)
      deallocate(tril_data)

end subroutine save_1rdm_orca

      subroutine load_1rdm_orca_header(rdm_file, NA)

            character(len=*), intent(in) :: rdm_file
            integer, intent(out) :: NA

            integer :: unit
            integer :: header(3)
            real(F64), dimension(:), allocatable :: tril_data


            open(newunit=unit, file=rdm_file, form='unformatted', &
                  access='stream', status='old')

            read(unit) header
            NA = header(1)

      end subroutine load_1rdm_orca_header

      subroutine load_2rdm_orca(rdm_file, rdm, with_header)
            character(len=*), intent(in) :: rdm_file
            real(F64), dimension(:,:,:,:), intent(out) :: rdm
            logical, intent(in) :: with_header

            integer :: unit
            integer :: header(3)
            integer :: dim_sq, dim, i, j, k, l
            double precision :: val
            integer :: n

            n = size(rdm,dim=1)
            unit = 10

            open(unit=unit, file=rdm_file, form='unformatted', &
                  access='stream', status='old')

            if (with_header == .true.)then
                  read(unit) header

                  dim_sq = header(1)
                  dim = int(sqrt(real(dim_sq, F64)))
                  if (n.ne.dim)then
                        print*, 'fcidump and rdm files incompatible, exiting'
                        print*, 'header', header
                        print*, 'n in fcidump', n
                        print*, 'n in rdm', dim
                        !stop
                  end if
            end if

            do k = 1, n
                do j = 1, n
                      do l = 1, n
                            do i = 1, n
                                  read(10) val
                                  rdm(l, k, i, j) = val

                            end do
                      end do
                end do
          end do
          close(10)

    end subroutine load_2rdm_orca

    subroutine save_2rdm_orca(rdm_file, rdm, with_header)
      character(len=*), intent(in) :: rdm_file
      real(F64), dimension(:,:,:,:), intent(in) :: rdm
      logical, intent(in) :: with_header

      integer :: unit
      integer :: header(3)
      integer :: n, dim_sq
      integer :: i, j, k, l
      double precision :: val

      n = size(rdm, dim=1)
      dim_sq = n*n
      unit = 10

      open(unit=unit, file=rdm_file, form='unformatted', &
            access='stream', status='replace')

      if (with_header == .true.) then
            header(1) = dim_sq
            header(2) = dim_sq
            header(3) = dim_sq/2
            write(unit) header
      end if

      do k = 1, n
            do j = 1, n
                  do l = 1, n
                        do i = 1, n
                              val = rdm(l, k, i, j)
                              write(unit) val
                        end do
                  end do
            end do
      end do

      close(unit)
end subroutine save_2rdm_orca


    subroutine check_natural(AuxData, IsDiag)
          type(TACppData), intent(inout) :: AuxData
          logical, intent(out) :: IsDiag
          integer :: i, j
          double precision :: err
          double precision :: tol = 1.d-6


          associate(NBasis=>AuxData%NBasis, NA=>AuxData%NA, rdm1=>AuxData%rdm1_full)

            err = zero
            IsDiag = .true.
            do i = 1, NA
                  do j = 1, NA
                        print*, 'rdm1(i,j)', i, j, rdm1(i,j)
                        if (i /= j) then
                              if (abs(rdm1(i, j)) > tol) then
                                    IsDiag = .false.
                                    err = err + abs(rdm1(i, j))
                              end if
                        end if
                  end do
                  if (.not. IsDiag) exit
            end do

            IsDiag = .false.
            if (abs(err) < tol)then
                  print*, 'Natural orbital basis, do not transform', err
                  IsDiag = .true.
            else
                  print*, 'Transform to natural orbital basis', err
            end if
            
          end associate
    end subroutine check_natural
    
    subroutine fix_rdm_basis(AuxData)
          type(TACppData), intent(inout) :: AuxData
          double precision, dimension(:,:), allocatable :: CAONO_small, CAONO_big
          double precision, dimension(:,:), allocatable :: C_s_act, C_b_act
          double precision, dimension(:,:), allocatable :: U_mat, U_t
          double precision, dimension(:,:), allocatable :: work, S_AO, tmp_mat
          double precision, dimension(:,:), allocatable :: U_svd, VT_svd, R_unitary
          double precision, dimension(:,:), allocatable :: A_svd
          double precision, dimension(:), allocatable :: S_vals, work_svd
          integer :: NA, NI, NBasis, info, lwork
          double precision :: work_query(1)
          double precision, dimension(:,:), allocatable :: w1, w2
          integer :: i, j

          external :: dgesvd

          associate(NI=>AuxData%NI, NA=>AuxData%NA, NBasis=>AuxData%NBasis)
            allocate(CAONO_small(NBasis, NBasis))
            allocate(CAONO_big(NBasis, NBasis))
            allocate(S_AO(NBasis, NBasis))

            call load_2x2_array_orca('C_small.bin', CAONO_small)
            call load_2x2_array_orca('C_big.bin', CAONO_big)
            call load_S_packed('S.bin', S_AO, NBasis)



            ! allocate(w1(NBasis, Nbasis))
            ! allocate(w2(NBasis, Nbasis))
            
            ! call real_ab(w1, S_AO, CAONO_big)
            ! call real_aTb(w2, CAONO_big, w1)

            ! do i = 1, NBasis
            !       do j = 1, NBasis
            !             if (abs(w2(i,j)).gt.1.d-5)then
            !                   print*, 'ww', i, j, w2(i, j)
            !             end if
            !       end do
            ! end do

            
            allocate(C_s_act(NBasis, NA))
            allocate(C_b_act(NBasis, NA))
            C_s_act = CAONO_small(:, NI+1:NI+NA)
            C_b_act = CAONO_big(:, NI+1:NI+NA)

            allocate(tmp_mat(NBasis, NA))
            allocate(U_mat(NA, NA))

            call real_ab(tmp_mat, S_AO, C_b_act)
            call real_aTb(U_mat, C_s_act, tmp_mat)

            allocate(U_svd(NA, NA))
            allocate(VT_svd(NA, NA))
            allocate(S_vals(NA))
            allocate(A_svd(NA, NA))

            A_svd = U_mat
            lwork = -1
            call dgesvd('A', 'A', NA, NA, A_svd, NA, S_vals, U_svd, NA, VT_svd, NA, work_query, lwork, info)
            lwork = max(1, int(work_query(1)))

            allocate(work_svd(lwork))
            A_svd = U_mat
            call dgesvd('A', 'A', NA, NA, A_svd, NA, S_vals, U_svd, NA, VT_svd, NA, work_svd, lwork, info)

            if (info .ne. 0) then
                  print*, 'Error in SVD in fix_rdm_basis, info = ', info
                  stop
            end if

            allocate(R_unitary(NA, NA))
            call real_ab(R_unitary, U_svd, VT_svd)

            allocate(U_t(NA, NA))
            U_t = transpose(R_unitary)

            allocate(work(NA, NA))

            call real_ab(work, AuxData%rdm1_p, U_t)
            call real_ab(AuxData%rdm1_p, R_unitary, work)

            call real_ab(work, AuxData%rdm1_m, U_t)
            call real_ab(AuxData%rdm1_m, R_unitary, work)

            call rdm2_MO_NO_trans(AuxData%rdm2_pp, U_t, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_pm, U_t, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mm, U_t, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mp, U_t, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_full, U_t, NA)
          end associate
    end subroutine fix_rdm_basis



    subroutine Trans2NO_AO(THCData, AuxData, CAOMO, CAONO, Flags, TwoEl)

          
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          double precision, dimension(:,:), intent(inout) :: CAOMO
          double precision, dimension(:,:), intent(inout) :: CAONO
          type(FlagsData), intent(in) :: Flags

          double precision, dimension(:), allocatable :: occ_temp
          integer, dimension(:), allocatable :: dy
          double precision, dimension(:,:), allocatable :: CMONO_NA, CMONO, work
          !double precision, optional, intent(inout) :: rdm2_full(:,:,:,:)
          double precision, optional, intent(inout) :: TwoEl(:)

          double precision :: err

          double precision :: tol = 1.d-9
          double precision :: suma, ETot0
          integer :: i, j, NI, NA, NIA
          double precision, dimension(:,:), allocatable :: rdm1_copy
          double precision, dimension(:,:), allocatable :: CAONOt



          NA = AuxData%NA

          associate(NBasis=>AuxData%NBasis, rdm1=>AuxData%rdm1_full)

            !call check_natural(AuxData, IsDiag)
            
            !if (.not. IsDiag)then

            allocate(rdm1_copy(NBasis, NBasis))
            rdm1_copy = rdm1
            allocate(occ_temp(NA))

            call symmetric_eigenproblem(occ_temp, rdm1_copy, NA, .true.)

            allocate(dy(NA))
            do i = 1, NA
                  dy(i) = i
            end do

            call dsort0(occ_temp, dy, NA, -2)
            occ_temp = occ_temp / two
            allocate(CMONO_NA(NA, NA))
          
            CMONO_NA = zero

            do i = 1, NA
                  CMONO_NA(:,i) = rdm1_copy(:, dy(i))
            end do
            !end if

            !----------------------------------------------------------------------------
            ! check for almost zero, or negative occupancies
            !----------------------------------------------------------------------------
            suma = zero
            NA = 0
            do i = 1, AuxData%NA
                  occ_temp(i) = abs(occ_temp(i))
                  suma = suma + occ_temp(i)
                  if (occ_temp(i)>1.d-8)then
                        NA = NA + 1
                  end if
            end do
            AuxData%NA=NA

            AuxData%NI = (AuxData%NEL -(suma*two)+1.d-2)/two
            NI = AuxData%NI

            AuxData%switch = 0
            if (NA.ne.AuxData%NA)then
                  Write(6, '(1x, "WARNING! The number of partially occupied orbitals &
                        different from NActDMRG read from ORCA. Some active orbitals must be unoccupied.", /)')
                  NA = AuxData%NA
                  AuxData%switch = 1
            end if          
          
            !allocate(AuxData%Occ(AuxData%NBasis))
            if (.not. allocated(AuxData%Occ)) allocate(AuxData%Occ(NBasis))
            print*,  'AuxData%NBasis', AuxData%NBasis
            AuxData%Occ = zero
            AuxData%Occ(1: AuxData%NI) = One
            AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp
            print*, 'AuxData%NI', AuxData%NI
            print*, 'AuxData%NA', AuxData%NA

            AuxData%NIA = AuxData%NI + AuxData%NA
            AuxData%NV = AuxData%NBasis - AuxData%NIA
            NIA = AuxData%NIA 
            
            !----------------------------------------------------------------------------
            ! print occupancies
            !----------------------------------------------------------------------------
            write(*, '(A10, A16)')'#ORB',  'Occupancy'
            
            suma = zero
            do i = 1, AuxData%NIA
                  write(*, '(I10, F16.6)') i, AuxData%Occ(i)
                  suma = suma + AuxData%Occ(i)
            end do
            
            write(*,'(A30, F20.15)') "sum of occupancies", suma
            
            
            !----------------------------------------------------------------------------
            ! prepare IndAux
            !----------------------------------------------------------------------------
            
            
            if (.not. allocated(AuxData%IndAux)) allocate(AuxData%IndAux(AuxData%NBasis))
            AuxData%IndAux = 2
            AuxData%IndAux(1:AuxData%NI) = 0
            AuxData%IndAux(AuxData%NI+1:AuxData%NIA) = 1

            if (allocated(AuxData%rdm2_pp)) call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
            if (allocated(AuxData%rdm2_mm)) call rdm2_MO_NO_trans(AuxData%rdm2_mm, CMONO_NA, NA)
            if (allocated(AuxData%rdm2_pm)) call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)
            if (allocated(AuxData%rdm2_mp)) call rdm2_MO_NO_trans(AuxData%rdm2_mp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_full, CMONO_NA, NA)

            call rdm1_MO_NO_trans(AuxData%rdm1_p, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_m, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_full, CMONO_NA, NA)

            ! if (present(rdm2_full)) then
            !       call rdm2_MO_NO_trans(rdm2_full, CMONO_NA, NA)
            ! end if


!            if (IsDiag==.true.)then
!                  print*, 'Data is already in Natural Orbital basis'
!                  CAONO = CAOMO

!            else

                  
                  allocate(CMONO(NBasis, NBasis))

                  CMONO = zero
                  CAONO = zero
                  do i = 1, NBasis
                        CMONO(i,i) = One
                  end do
                  
                  CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 
                  
                  call real_ab(CAONO, CAOMO, CMONO)

                  if (Flags%ITwoel == 1)then
                        print*, 'Transforming 2-el (AO->NO) integrals'
                        ! Use CAONO (AO->NO) because TwoEl is in AO basis (PySCF/Orca-Out)

                        allocate(CAONOt(NBasis, NBasis))
                        CAONOt = transpose(CAONO)
                        !CAONOt = transpose(CMONO)
                        call TwoNO1(TwoEl,CAONOt,NBasis,AuxData%NInte2)
                        allocate(work(NBasis, NBasis))
                        call real_ab(work, AuxData%HNO0, CAONO)
                        call real_atb(AuxData%HNO0, CAONO, work)
                        deallocate(work)
                        !call real_ab(work, AuxData%HNO0, CMONO)
                        !call real_atb(AuxData%HNO0, CMONO, work)
                  end if
 !           end if


            !----------------------------------------------------------------------------
            ! Canonicalization inside THC_init
            ! On Exit CAONO is canonical in inactive/virtual blocks
            !----------------------------------------------------------------------------
            ! if (Flags%ITwoEl > 1)then
            !       call THC_init2(Flags, THCData, AuxData,  AObasis, System, CAONO)
            ! end if
            
            ! allocate(work(NBasis, NBasis))
            ! call real_ab(work, AuxData%HNO0, CAONO)
            ! call real_atb(AuxData%HNO0, CAONO, work)

            
            ! ETot0 = zero
            
            
            ! if (Flags%JOBTYPE .ne. JOB_TYPE_MP2 .and. Flags%JOBTYPE .ne. JOB_TYPE_SRMP2) then                    
            !       do i = 1, NBasis
            !             ETot0 = ETot0 + two* AuxData%Occ(i) * AuxData%HNO0(i,i)                  
            !             !write(*, '(I5, 2F20.15)') i, AuxData%HNO0(i,i), AuxData%Occ(i)
            !       end do
            ! else
            !       do i = 1, NBasis
            !             ETot0 = ETot0 + two* AuxData%Occ_rohf(i) * AuxData%HNO0(i,i)                  
            !       end do
            ! end if
            ! print*, 'etot1 z HNO_ext', ETot0
            
          end associate

    end subroutine Trans2NO_AO


    subroutine Trans2NO_MO(THCData, AuxData, CAONO, Flags, TwoEl)

          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          double precision, dimension(:,:), intent(inout) :: CAONO
          type(FlagsData), intent(in) :: Flags

          double precision, dimension(:), allocatable :: occ_temp
          integer, dimension(:), allocatable :: dy
          double precision, dimension(:,:), allocatable :: CMONO_NA, CMONO, work
          double precision, optional, intent(inout) :: TwoEl(:)

          double precision :: suma
          integer :: i, j, NI, NA, NIA

          double precision, dimension(:,:), allocatable :: rdm1_copy
          double precision, dimension(:,:), allocatable :: CAONOt


          print*, 'AuxData%NA', AuxData%NA
          NA = AuxData%NA
          associate(NBasis=>AuxData%NBasis, rdm1=>AuxData%rdm1_full)
            
            do i = 1, NA
                  do j = 1, NA
                        write(*,'(A7, 2I3, 3F10.7)') 'gooow1', i, j, AuxData%rdm1_full(i, j), AuxData%rdm1_p(i, j), AuxData%rdm1_m(i, j)
                  end do
            end do


            allocate(rdm1_copy(NBasis, NBasis))
            rdm1_copy = rdm1
            allocate(occ_temp(NA))

            call symmetric_eigenproblem(occ_temp, rdm1_copy, NA, .true.)

            allocate(dy(NA))
            do i = 1, NA
                  dy(i) = i
            end do

            call dsort0(occ_temp, dy, NA, -2)
            occ_temp = occ_temp / two
            allocate(CMONO_NA(NA, NA))
          
            CMONO_NA = zero

            do i = 1, NA
                  CMONO_NA(:,i) = rdm1_copy(:, dy(i))
!                  CMONO_NA(:,i) = rdm1_copy(:, i)
            end do

            ! suma = zero
            ! NA = 0
            ! do i = 1, AuxData%NA
            !       occ_temp(i) = abs(occ_temp(i))
            !       suma = suma + occ_temp(i)
            !       if (occ_temp(i)>1.d-8)then
            !             NA = NA + 1
            !       end if
            ! end do
            ! AuxData%NA=NA

            ! AuxData%NI = (AuxData%NEL -(suma*two)+1.d-2)/two
            ! NI = AuxData%NI

            ! AuxData%switch = 0
            ! if (NA.ne.AuxData%NA)then
            !       Write(6, '(1x, "WARNING! The number of partially occupied orbitals &
            !             different from NActDMRG read from ORCA. Some active orbitals must be unoccupied.", /)')
            !       NA = AuxData%NA
            !       AuxData%switch = 1
            ! end if          
          
            print*,  'AuxData%NBasis', AuxData%NBasis
            AuxData%Occ = zero
            AuxData%Occ(1: AuxData%NI) = One
            AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp
            print*, 'AuxData%NI', AuxData%NI
            print*, 'AuxData%NA', AuxData%NA, NA
            
            AuxData%NIA = AuxData%NI + AuxData%NA
            AuxData%NV = AuxData%NBasis - AuxData%NIA
            NIA = AuxData%NIA
            NI = AuxData%NI
            
            write(*, '(A10, A16)')'#ORB',  'Occupancy'
            
            suma = zero
            do i = 1, AuxData%NIA
                  write(*, '(I10, F16.6)') i, AuxData%Occ(i)
                  suma = suma + AuxData%Occ(i)
            end do
            
            write(*,'(A30, F20.15)') "sum of occupancies", suma

            !----------------------------------------------------------------------------
            ! prepare IndAux
            !----------------------------------------------------------------------------
            if (.not. allocated(AuxData%IndAux)) allocate(AuxData%IndAux(AuxData%NBasis))
            AuxData%IndAux = 2
            AuxData%IndAux(1:AuxData%NI) = 0
            AuxData%IndAux(AuxData%NI+1:AuxData%NIA) = 1
            
            call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_full, CMONO_NA, NA)

            call rdm1_MO_NO_trans(AuxData%rdm1_p, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_m, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_full, CMONO_NA, NA)

            allocate(CMONO(NBasis, NBasis))
            CMONO = zero
            do i = 1, NBasis
                  CMONO(i,i) = One
            end do


            do i = 1, NA
                  do j = 1, NA
                        write(*,'(A7, 2I3, 3F10.7)') 'gooow2', i, j, AuxData%rdm1_full(i, j), AuxData%rdm1_p(i, j), AuxData%rdm1_m(i, j)
                  end do
            end do



            CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 

            ! In MO basis case (ORCA In-core), CAONO is set to CMONO for subsequent use
            CAONO = CMONO
            
            ! Save CMONO for back-transformation in rdm_dump
            if (allocated(AuxData%CMONO)) deallocate(AuxData%CMONO)
            allocate(AuxData%CMONO(NBasis, NBasis))
            AuxData%CMONO = CMONO

            if (Flags%ITwoel == 1)then
                  print*, 'Transforming 2-el (MO->NO) integrals'
                  allocate(CAONOt(NBasis, NBasis))
                  ! Use CMONO (MO->NO) because TwoEl is in MO basis (ORCA In-core)
                  CAONOt = transpose(CMONO)
                  call TwoNO1(TwoEl,CAONOt,NBasis,AuxData%NInte2)
                  
                  allocate(work(NBasis, NBasis))
                  call real_ab(work, AuxData%HNO0, CMONO)
                  call real_atb(AuxData%HNO0, CMONO, work)
                  deallocate(work)
            end if

          end associate

    end subroutine Trans2NO_MO

    subroutine twono1_4fold(tno, caono, nbasis, ninte2)
          implicit none
          integer, intent(in) :: nbasis
          integer(8), intent(in) :: ninte2
          double precision, intent(inout) :: tno(ninte2)
          double precision, intent(in) :: caono(nbasis, nbasis)

          integer :: i, j, k, l, a, b, c, d
          integer :: ij, kl, ab, cd
          integer(8) :: naddr
          double precision :: twozet, v, yval
          double precision, parameter :: zero = 0.0d0
          double precision, parameter :: tol2 = 1.d-12
          integer :: ninte1

          double precision, allocatable :: x(:,:,:)
          double precision, allocatable :: y(:,:)

          ninte1 = nbasis * nbasis

          allocate(x(nbasis, nbasis, ninte1))
          allocate(y(ninte1, ninte1))

          ! first index
          x = zero
          
          do ij = 1, ninte1
              ! decode ij -> i, j
              ! ij = i + (j-1)*n
              i = mod(ij-1, nbasis) + 1
              j = (ij-1) / nbasis + 1

              do kl = 1, ij
                   ! decode kl -> k, l
                   k = mod(kl-1, nbasis) + 1
                   l = (kl-1) / nbasis + 1
                   
                   naddr = gmap_4fold(i, j, k, l, nbasis)
                   twozet = tno(naddr)
                   
                   if (abs(twozet) > tol2) then
                        do a = 1, nbasis
                             x(a, j, kl) = x(a, j, kl) + caono(a, i) * twozet
                        end do
                        
                        if (ij /= kl) then
                             do a = 1, nbasis
                                 x(a, l, ij) = x(a, l, ij) + caono(a, k) * twozet
                             end do
                        end if
                   end if
              end do
          end do
          
          ! second index
          y = zero
          
          do kl = 1, ninte1
              do j = 1, nbasis
                   do a = 1, nbasis
                        v = x(a, j, kl)
                        if (abs(v) > tol2) then
                             do b = 1, nbasis
                                  ab = a + (b-1)*nbasis
                                  y(ab, kl) = y(ab, kl) + caono(b, j) * v
                             end do
                        end if
                   end do
              end do
          end do
          
          ! third index
          x = zero
          
          do kl = 1, ninte1
               k = mod(kl-1, nbasis) + 1
               l = (kl-1) / nbasis + 1
               
               do ab = 1, ninte1
                    yval = y(ab, kl)
                    if (abs(yval) > tol2) then
                         do c = 1, nbasis
                              x(c, l, ab) = x(c, l, ab) + caono(c, k) * yval
                         end do
                    end if
               end do
          end do
          
          ! fourth index
          ! Cannot Zero TNO globally because of aliasing (we might overwrite valid computed parts if we zero then calculate).
          ! BUT we are regenerating TNO fully.
          ! Since aliases point to same location, we might write same value multiple times.
          ! This is safe.
          tno = zero
          
          do ab = 1, ninte1
               do cd = 1, ab
                  c = mod(cd-1, nbasis) + 1
                  d = (cd-1) / nbasis + 1
                  a = mod(ab-1, nbasis) + 1
                  b = (ab-1) / nbasis + 1

                  naddr = gmap_4fold(c, d, a, b, nbasis)
                    
                  v = zero
                  do l = 1, nbasis
                      v = v + caono(d, l) * x(c, l, ab)
                  end do
                  
                  tno(naddr) = v
               end do
          end do

          deallocate(x)
          deallocate(y)

    end subroutine twono1_4fold

    subroutine Trans2NO_MO_ducc(AuxData, Flags, Ints_full, Ints_ducc)
              
          type(TACppData), intent(inout) :: AuxData
          type(FlagsData), intent(in) :: Flags
          type(TInts), intent(inout) :: Ints_full
          type(TInts), intent(inout) :: Ints_ducc          


          double precision, dimension(:), allocatable :: occ_temp
          integer, dimension(:), allocatable :: dy
          double precision, dimension(:,:), allocatable :: CMONO_NA, CMONO, work


          double precision :: suma
          integer :: i, j, NI, NA, NIA

          double precision, dimension(:,:), allocatable :: rdm1_copy
          double precision, dimension(:,:), allocatable :: CAONOt, CAONOt_NA


          print*, 'AuxData%NA', AuxData%NA
          NA = AuxData%NA
          associate(NBasis=>AuxData%NBasis, rdm1=>AuxData%rdm1_full)

            allocate(rdm1_copy(NBasis, NBasis))
            rdm1_copy = rdm1
            allocate(occ_temp(NA))

            call symmetric_eigenproblem(occ_temp, rdm1_copy, NA, .true.)

            allocate(dy(NA))
            do i = 1, NA
                  dy(i) = i
            end do

            call dsort0(occ_temp, dy, NA, -2)
            occ_temp = occ_temp / two
            allocate(CMONO_NA(NA, NA))
          
            CMONO_NA = zero

            do i = 1, NA
                  CMONO_NA(:,i) = rdm1_copy(:, dy(i))
            end do


          
            print*,  'AuxData%NBasis', AuxData%NBasis
            AuxData%Occ = zero
            AuxData%Occ(1: AuxData%NI) = One
            AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp
            print*, 'AuxData%NI', AuxData%NI
            print*, 'AuxData%NA', AuxData%NA, NA
            
            AuxData%NIA = AuxData%NI + AuxData%NA
            AuxData%NV = AuxData%NBasis - AuxData%NIA
            NIA = AuxData%NIA
            NI = AuxData%NI
            
            write(*, '(A10, A16)')'#ORB',  'Occupancy'
            
            suma = zero
            do i = 1, AuxData%NIA
                  write(*, '(I10, F16.6)') i, AuxData%Occ(i)
                  suma = suma + AuxData%Occ(i)
            end do
            
            write(*,'(A30, F20.15)') "sum of occupancies", suma
            
            call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_full, CMONO_NA, NA)

            call rdm1_MO_NO_trans(AuxData%rdm1_p, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_m, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_full, CMONO_NA, NA)

            print*, 'rmd1_p'
            do i = 1, NA
                  do j = 1, NA
                        if (abs(AuxData%rdm1_p(i,j)).gt.1.d-5)then
                              write(*, '(A10, 2I5, F12.6)') 'rrr', i, j, AuxData%rdm1_p(i,j)
                        end if
                  end do
            end do

            print*, 'rmd1_m'
            do i = 1, NA
                  do j = 1, NA
                        if (abs(AuxData%rdm1_m(i,j)).gt.1.d-5)then
                              write(*, '(A10, 2I5, F12.6)') 'rrr', i, j, AuxData%rdm1_m(i,j)
                        end if
                  end do
            end do

            allocate(CMONO(NBasis, NBasis))
            CMONO = zero
            do i = 1, NBasis
                  CMONO(i,i) = One
            end do

            CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 

            ! In MO basis case (ORCA In-core), CAONO is set to CMONO for subsequent use
!            CAONO = CMONO
            
            ! Save CMONO for back-transformation in rdm_dump
            if (allocated(AuxData%CMONO)) deallocate(AuxData%CMONO)
            allocate(AuxData%CMONO(NBasis, NBasis))
            AuxData%CMONO = CMONO

 
            print*, 'Transforming 2-el (MO->NO) integrals'
            allocate(CAONOt(NBasis, NBasis))
            allocate(CAONOt_NA(NA,NA))
            ! Use CMONO (MO->NO) because TwoEl is in MO basis (ORCA In-core)
            CAONOt = transpose(CMONO)
            CAONOt_NA = transpose(CMONO_NA)
            call twono1_4fold(Ints_full%ints2e_aa,CAONOt,NBasis,Ints_full%ints2e_dim)
            call twono1_4fold(Ints_full%ints2e_ab,CAONOt,NBasis,Ints_full%ints2e_dim)

            call twono1_4fold(Ints_ducc%ints2e_aa,CAONOt_NA,NA,Ints_ducc%ints2e_dim)
            call twono1_4fold(Ints_ducc%ints2e_ab,CAONOt_NA,NA,Ints_ducc%ints2e_dim)
            
            allocate(work(NBasis, NBasis))
            call real_ab(work, Ints_full%ints1e_aa, CMONO)
            call real_atb(Ints_full%ints1e_aa, CMONO, work)
            deallocate(work)
            allocate(work(NA,NA))
            call real_ab(work, Ints_ducc%ints1e_aa, CMONO_NA)
            call real_atb(Ints_ducc%ints1e_aa, CMONO_NA, work)
            
            deallocate(work)

          end associate

    end subroutine Trans2NO_MO_ducc

    subroutine Trans2NO_MO_spinres(AuxData, Flags, Ints)
          
          type(TACppData), intent(inout) :: AuxData
          type(FlagsData), intent(in) :: Flags
          type(TInts), intent(inout) :: Ints

          double precision, dimension(:), allocatable :: occ_temp
          integer, dimension(:), allocatable :: dy
          double precision, dimension(:,:), allocatable :: CMONO_NA, CMONO, work

          double precision :: suma
          integer :: i, j, NI, NA, NIA

          double precision, dimension(:,:), allocatable :: rdm1_copy
          double precision, dimension(:,:), allocatable :: CAONOt

          print*, 'Trans2NO_MO_spinres: AuxData%NA', AuxData%NA
          NA = AuxData%NA
          associate(NBasis=>AuxData%NBasis, rdm1=>AuxData%rdm1_full)

            allocate(rdm1_copy(NBasis, NBasis))
            rdm1_copy = rdm1
            allocate(occ_temp(NA))

            call symmetric_eigenproblem(occ_temp, rdm1_copy, NA, .true.)

            allocate(dy(NA))
            do i = 1, NA
                  dy(i) = i
            end do

            call dsort0(occ_temp, dy, NA, -2)
            occ_temp = occ_temp / two
            allocate(CMONO_NA(NA, NA))
          
            CMONO_NA = zero

            do i = 1, NA
                  CMONO_NA(:,i) = rdm1_copy(:, dy(i))
            end do
          
            print*,  'AuxData%NBasis', AuxData%NBasis
            AuxData%Occ = zero
            AuxData%Occ(1: AuxData%NI) = One
            AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp
            print*, 'AuxData%NI', AuxData%NI
            print*, 'AuxData%NA', AuxData%NA, NA
            
            AuxData%NIA = AuxData%NI + AuxData%NA
            AuxData%NV = AuxData%NBasis - AuxData%NIA
            NIA = AuxData%NIA
            NI = AuxData%NI
            
            write(*, '(A10, A16)')'#ORB',  'Occupancy'
            
            suma = zero
            do i = 1, AuxData%NIA
                  write(*, '(I10, F16.6)') i, AuxData%Occ(i)
                  suma = suma + AuxData%Occ(i)
            end do
            
            write(*,'(A30, F20.15)') "sum of occupancies", suma
            
            call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_mp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_full, CMONO_NA, NA)

            call rdm1_MO_NO_trans(AuxData%rdm1_p, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_m, CMONO_NA, NA)
            call rdm1_MO_NO_trans(AuxData%rdm1_full, CMONO_NA, NA)

            allocate(CMONO(NBasis, NBasis))
            CMONO = zero
            do i = 1, NBasis
                  CMONO(i,i) = One
            end do
            CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 


            ! Save CMONO for back-transformation in rdm_dump
            if (allocated(AuxData%CMONO)) deallocate(AuxData%CMONO)
            allocate(AuxData%CMONO(NBasis, NBasis))
            AuxData%CMONO = CMONO

 
            print*, 'Transforming 2-el (MO->NO) integrals'
            allocate(CAONOt(NBasis, NBasis))

            ! Use CMONO (MO->NO) because TwoEl is in MO basis (ORCA In-core)
            CAONOt = transpose(CMONO)
            
            call twono1_4fold(Ints%ints2e_aa,CAONOt,NBasis,Ints%ints2e_dim)
            ! Also transform bb and ab
            call twono1_4fold(Ints%ints2e_bb,CAONOt,NBasis,Ints%ints2e_dim)
            call twono1_4fold(Ints%ints2e_ab,CAONOt,NBasis,Ints%ints2e_dim)
            
            allocate(work(NBasis, NBasis))
            call real_ab(work, Ints%ints1e_aa, CMONO)
            call real_atb(Ints%ints1e_aa, CMONO, work)
            ! Also transform bb
            call real_ab(work, Ints%ints1e_bb, CMONO)
            call real_atb(Ints%ints1e_bb, CMONO, work)
            
            deallocate(work)

            ! Update HNO0 to match transformed ints1e_aa
            AuxData%HNO0 = Ints%ints1e_aa

          end associate

    end subroutine Trans2NO_MO_spinres

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


    subroutine DALTON_wrapper(Occ, XKin, ENuc, TwoNO, &
          NBasis, NInte1, NInte2, AuxData, UMOAO, CAONO)

          double precision, dimension(:),   intent(in) :: Occ
          double precision, dimension(:),   intent(in) :: XKin
          double precision, intent(in) :: ENuc
          double precision, dimension(:),intent(in) :: TwoNO
          integer, intent(in) :: Nbasis
          double precision, dimension(NBasis, NBasis), intent(in) :: UMOAO
          double precision, dimension(:,:), allocatable, intent(out) :: CAONO
          integer, intent(in) :: NInte1, NInte2
          type(TACppData), intent(out) :: AuxData
          integer, dimension(:), allocatable :: IActOrb,InActOrb

          AuxData%Nbasis = Nbasis
          AuxData%NInte1 = NInte1
          AuxData%NInte2 = NInte2
          AuxData%ENuc = ENuc
          allocate(AuxData%TwoNO(NInte2))
          allocate(AuxData%XOne(NInte1))
          AuxData%TwoNO = TwoNO
          AuxData%XOne = XKin
          allocate(AuxData%Occ(NBasis))

          allocate(CAONO(NBasis, NBasis))
          CAONO = UMOAO
          
          allocate(IActOrb(NBasis))
          allocate(InActOrb(NBasis))
          call read_1rdm_dalton(AuxData%Occ,IActOrb,InActOrb, AuxData%NI,AuxData%NA,NBasis)

          AuxData%NIA = AuxData%NI + AuxData%NIA
          AuxData%NV = AuxData%NBasis - AuxData%NIA
          AuxData%Occ = Occ

          AuxData%DALTON = 1


          allocate(AuxData%IndAux(NBasis))
          AuxData%IndAux = 2
          AuxData%IndAux(1:AuxData%NI) = 0
          AuxData%IndAux(AuxData%NI+1:AuxData%NIA) = 1



          
    end subroutine DALTON_wrapper

    subroutine GEOM_init(Flags, AuxData, AObasis, System)

          
          type(FlagsData), intent(in) :: Flags
          type(TACppData), intent(inout) :: AuxData

          type(TAOBASIS), intent(out) :: AObasis
          type(TSystem), intent(out) :: System

          character(:), allocatable :: binPath
          integer :: units, nbasis
          logical :: SortAngularMomenta

          units = Flags%IUnits

          SortAngularMomenta = .true.
          call print_section('Geometry and basis setup')
          call print_info('Basis set path', Flags%BasisSetPath)
          call geom_ReadSystemBasis(System, AObasis, Flags, SortAngularMomenta, Flags%IUnits)
          call print_info('Number of atoms', System%NAtoms)
          call print_info('Number of electrons', System%NElectrons)

          call sys_NuclearRepulsion(AuxData%ENuc,System)
          Nbasis = AObasis%NAOSpher
          AuxData%NBasis = AObasis%NAOSpher
          AuxData%NEL = System%NElectrons
          call print_info('AuxData%NEL', AuxData%NEL)
          call print_info('AuxData%NBasis', AuxData%NBasis)
          
    end subroutine GEOM_init

    subroutine check_energy_incore(AuxData, TwoEl, spin_sep)

          type(TACppData), intent(inout) :: AuxData
          double precision, intent(in) :: TwoEl(:)
          logical, intent(in) :: spin_sep
          integer(8), external :: NAddr3
          integer, external :: NAddrRDM
          double precision :: etot, x, val, this
          double precision :: eone, etwo
          double precision, dimension(:), allocatable :: rdm2
          integer :: i, j, k, l, rdm_dim
          integer :: u, ios
          integer :: p, q, r, s
          double precision, dimension(:,:), allocatable :: DM1          
          
          associate(NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, &
                NBasis=>AuxData%NBasis, NIA=>AuxData%NIA, Occ=>AuxData%Occ, &
                IndAux=>AuxData%IndAux)


            allocate(DM1(NBasis, NBasis))
            DM1 = zero
                        
            do i = 1, NI
                  DM1(i, i) = one
            end do

            do i = 1, NA
                  do j = 1, NA
                        DM1(NI+i, NI+j) = (AuxData%rdm1_p(i, j) + AuxData%rdm1_m(i, j))/two
                  end do
            end do

            eone = zero
            do p = 1, NIA
                  do q = 1, NIA
                        if (abs(DM1(p,q)) > 1.d-10) then
                              eone = eone + two * DM1(p, q) * AuxData%HNO0(p, q)
                              ! if (p>NI .and.q>NI)then
                              !       write(*,'(A10, 2I5, 4F15.8)') 'lasa', p, q, DM1(p,q), AuxData%rdm1_m(p-NI,q-NI), AuxData%rdm1_p(p-NI,q-NI), AuxData%Occ(p)
                              ! else
                              !       write(*,'(A10, 2I5, 2F15.8)') 'lasa', p, q, DM1(p,q), AuxData%Occ(p)
                              ! end if
                        end if
                  end do
            end do


            Write(6,'(1X,''CASSCF Energy (one-electron) calculated'',X,F15.8)')eone
            etwo = zero
            if (spin_sep ==.false.)then
                  do p = 1, NI+NA
                        do q = 1, NI+NA
                              do r = 1, NI+NA
                                    do s = 1, NI+NA
                                          this = TwoEl(NAddr3(p, r, q, s))

                                          val = zero
                                          if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                                !val = val+  rdm2(naddrrdm(p-NI, q-NI, r-NI, s-NI, NA))
                                                val = val +  frac12 * AuxData%rdm2_full(p-NI, q-NI, r-NI, s-NI)
                                          else
                                                if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                      val =  val + two * Occ(p) * Occ(q)
                                                end if

                                                if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                      val = val  -Occ(p) * Occ(q)
                                                end if
                                          end if
                                          etwo = etwo + val * this
                                          ! if (abs(val * this).gt.1.d-5)then
                                          !       write(*, '(A20, 4I5, 3F20.10)')'energy-inter', p, q, r, s, val, this, etot
                                          ! end if
                                    end do
                              end do
                        end do
                  end do
            else
                  do p = 1, NI+NA
                        do q = 1, NI+NA
                              do r = 1, NI+NA
                                    do s = 1, NI+NA
                                          this = TwoEl(gmap(p, r, q, s))                                          
                                          val = zero
                                          if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                                val = val +  frac12 * (AuxData%rdm2_pp(p-NI, q-NI, r-NI, s-NI) &
                                                      + AuxData%rdm2_mm(p-NI, q-NI, r-NI, s-NI)&
                                                      + AuxData%rdm2_pm(p-NI, q-NI, r-NI, s-NI)&
                                                      + AuxData%rdm2_mp(p-NI, q-NI, r-NI, s-NI))

                                          else

                                          val = val + two * DM1(p, r) * DM1(q, s)
                                          val = val - DM1(p, s) * DM1(q, r)
                                    end if
                                          etwo = etwo + val * this

                                    end do
                              end do
                        end do
                  end do
            end if
            !print*, 'etot', etot

            Write(6,'(1X,''CASSCF Two-electron Energy calculated'',X,F15.8)')etwo
            Write(6,'(1X,''CASSCF Energy (w/o ENuc-q) calculated'',X,F15.8)')eone+etwo
            Write(6,'(1X,''Total CASSCF Energy calculated '',5X,F15.8)')eone+etwo+ AuxData%enuc


          end associate
                
    end subroutine check_energy_incore

    subroutine THC_init2(Flags, THCData, AuxData,  AObasis, System, CAONO_IN)

          type(FlagsData), intent(in) :: Flags
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          type(TAOBASIS), intent(in)  :: AObasis
          type(TSystem), intent(in) :: System

          double precision, dimension(:,:), intent(inout) :: CAONO_IN
          double precision, dimension(:,:), allocatable ::  CAONO

          double precision :: CholThr, THCThr
          integer :: CholAccu
          integer :: i, j, NAO, ij, ab, u
          double precision, allocatable :: Xgp(:,:)
          double precision, allocatable :: XgpErf(:,:), xmumat(:,:)
          integer :: NA, NI, NV, NIA, nnn
          integer :: unit
          integer :: units, nbasis
          type (tclock) :: timer
          double precision :: val1, val2, val3
          logical :: x2c

          ! Fock matrix diag                                                                                                                                                                                                                                                                                                
          double precision, dimension(:,:), allocatable :: Fockij, Fockvw
          double precision, dimension(:,:), allocatable :: Cpi_extao, Cpa_extao, Cpv_extao
          double precision, dimension(:,:), allocatable :: work, H0_extao

          !---------------------testing energy                                                                                                                                                                                                                                                                              
          integer :: p, q, r, s, k, l, a, b, ii, nn, t
          double precision, allocatable :: Rkab(:,:,:), Rkcd(:,:,:), Rkef(:,:,:)
          double precision :: this, ETot, etot0, val, val_this, thisfr, thistr, thistr1
          double precision ::  this0, that0, then0, thisJ, thisK


          NA = AuxData%NA
          NI = AuxData%NI
          NV = AuxData%NV
          NIA = NI+NA
          nbasis = AuxData%NBasis

          call auto2e_init()
          call clock_start(timer)
          
          CholThr = Flags%DCholeskyThr
          THCThr = Flags%DTHCthr
          CholAccu = Flags%ICholeskyAccu

          call print_section('Cholesky decomposition')

          if (CholThr < zero.or.THCThr <zero)then
                call thc_gammcor_XZ(THCData%Xgp, THCData%Zgk, AOBasis, System, CholAccu)

                if (Flags%IDBBSC == 2 .or. Flags%IFunSR==4)then
                      call thc_gammcor_XZ(THCData%XgpErf, THCData%ZgkErf, AOBasis, System, CholAccu, Omega=AuxData%Omega)
                end if

          else
                print*, 'CholeskyThreshold', CholThr
                print*, 'THCThreshold', THCThr
                call thc_gammcor_XZ(THCData%Xgp, THCData%Zgk, AOBasis, System,CholAccu, CholThr, THCThr)

                if (Flags%IDBBSC == 2 .or.Flags%IFunSR==4)then

                      call thc_gammcor_XZ(THCData%XgpErf, THCData%ZgkErf, AOBasis, System, CholAccu, CholThr, THCThr, Omega=AuxData%Omega)
                end if

                print*, 'CholeskyThreshold', CholThr
                print*, 'THCThreshold', THCThr
          end if
          allocate(CAONO(nbasis, nbasis))
          nao = nbasis
          allocate(work(nao, nbasis))

          if (AuxData%PYSCF==1 .or. AuxData%ORCA==1)then
                CAONO = CAONO_IN
          else
                CAONO = transpose(CAONO_IN)
          end if

          !          x2c = .true.
          !x2c = .false.

          !if (x2c == .true.)then
          if (THCData%H0external == .true.)then
                print*, 'external true'
                !print*, 'przed canonicalize'
                !print*, 'nao', nao, nbasis
                call canonicalize(CAONO, THCData%fij, THCData%fvw, AuxData, THCData, THCData%Xgp, AObasis, System, ext=.true.)

                allocate(H0_extao(nao, nao))
                H0_extao = AuxData%HNO0
                ! do i =1, 5
                !       do j = 1, 5
                !             print*, 'AuxData%HNO0(i,j)=', i, j, AuxData%HNO0(i, j)
                !       end do
                ! end do
                !print*, 'po canonicalize'
                !print*, 'nao', nao, nbasis
                !print*, 'size CAONO', size(CAONO, dim=1), size(CAONO, dim=2)
                !print*, 'size AuxData%HNO0', size(AuxData%HNO0, dim=1), size(AuxData%HNO0, dim=2)
                call real_ab(work, H0_extao, CAONO)
                call real_atb(AuxData%HNO0, CAONO, work)
                allocate(AuxData%HNO0_THC(nbasis, nbasis))
                AuxData%HNO0_THC = AuxData%HNO0 
          else
                !
                ! this is for the regular hamilotnian
                !
                !print*, 'canonicalize'
                !call canonicalize(CAONO, THCData%fij, THCData%fvw, AuxData, THCData, THCData%Xgp, AObasis, System, )
                call canonicalize(CAONO, THCData%fij, THCData%fvw, AuxData, THCData, THCData%Xgp, AObasis, System, ext=.false.)
                !print*, 'po can'
                ! open(unit=10,file='CAONO_can.bin',form='unformatted')
                ! write(10) NBasis
                ! write(10) CAONO
                ! close(10)

                allocate(H0_extao(nao, nao))
                allocate(AuxData%HNO0_THC(nbasis, nbasis))

                call ints1e_gammcor_H0_extao(H0_extao, AObasis, System, THCData%ExternalOrdering)
                !print*, 'h01', H0_extao(1,1)

                call real_ab(work, H0_extao, CAONO)
                call real_atb(AuxData%HNO0_THC, CAONO, work)
                AuxData%HNO0 = AuxData%HNO0_THC
          end if

          THCData%NTHC=size(THCData%Xgp,dim=1)
          THCData%NChol=size(THCData%Zgk,dim=2)
          allocate(THCData%Xga(THCData%NTHC,NBasis))

          Call thc_gammcor_Xga(THCData%Xga, THCData%Xgp, CAONO,&
                AOBasis, THCData%ExternalOrdering)

          if (Flags%IDBBSC == 2 .or. Flags%IFunSR==4)then

                THCData%NTHCErf=size(THCData%XgpErf,dim=1)
                THCData%NCholErf=size(THCData%ZgkErf,dim=2)
                
                allocate(THCData%XgaErf(THCData%NTHCErf,NBasis))
                allocate(THCData%TXgaErf(THCData%NTHCErf,NBasis))
                allocate(THCData%TXga(THCData%NTHC,NBasis))
                THCData%TXgaErf = zero
                THCData%TXga = zero
                Call thc_gammcor_Xga(THCData%XgaErf, THCData%XgpErf, CAONO,&
                      AOBasis, THCData%ExternalOrdering)

                if(Flags%IFunSR==4)then

                      allocate(THCData%J_SR(NBasis, NBasis))
                      call calc_J_SR(CAONO, THCData, AuxData, AObasis, System)
                end if
                
          end if

          
          if (Flags%JOBTYPE .ne. JOB_TYPE_MP2 .and. Flags%JOBTYPE .ne. JOB_TYPE_SRMP2) then
                ETot0 = zero
                do i = 1, NIA
                      ETot0 = ETot0 + two* AuxData%Occ(i) * AuxData%HNO0(i,i)
                      ! if (abs(AuxData%Occ(i) * AuxData%HNO0_THC(i,i)).gt.1.d-1)then
                      !       write(*, '(A10, I5, 3F20.8)')'etot-1', i, AuxData%Occ(i), AuxData%HNO0_THC(i,i), etot0
                      ! end if
                end do
          else
                print*, 'Flags%JOBTYPE', Flags%JOBTYPE
                ETot0 = zero
                do i = 1, NIA
                      ETot0 = ETot0 + AuxData%Occ_rohf(i) * AuxData%HNO0(i,i)
                      ! if (abs(AuxData%Occ(i) * AuxData%HNO0_THC(i,i)).gt.1.d-1)then
                      !       write(*, '(A10, I5, 2F20.15)')'etot-1', i, AuxData%Occ(i), AuxData%HNO0_THC(i,i)
                      ! end if
                end do
          end if
          AuxData%ECAS_oneelectr = ETot0

          AuxData%HNO0 = AuxData%HNO0_THC


          CAONO_IN = CAONO
          AuxData%OnlyEnergy = .false.
          if (AuxData%OnlyEnergy == .true.)then

                ! if (Flags%IDBBSC == 2)then
                !       allocate(Rkab(THCData%NCholErf,nao, nao))
                !       allocate(Rkef(THCData%NCholErf,nao, nao))                                                                                                                                                                                   
                !       Call thc_gammcor_Rkab_2(Rkab, THCData%XgaErf, THCData%XgaErf, THCData%ZgkErf, NBasis, NBasis,&                                
                !             THCData%NCholErf, THCData%NTHCErf)
                ! end if

                allocate(Rkab(THCData%NChol,nao, nao))         
                Call thc_gammcor_Rkab_2(Rkab, THCData%Xga, THCData%Xga, THCData%Zgk, NBasis, NBasis,&                
                      THCData%NChol, THCData%NTHC)
                
                
                associate(Occ=>AuxData%Occ, IndAux=>AuxData%IndAux)
                  ! if (Flags%IDBBSC == 2)then
                  !       allocate(xmumat(nbasis, nbasis))
                  !       open(unit=10,file='xmumat.bin',form='unformatted')
                  !       read(10) nnn
                  !       read(10) XMuMAT
                  !       close(10)
                        
                  !       Call dgemm('N','N',THCData%NCholErf*NBasis,NBasis,NBasis,One,Rkab,THCData%NCholErf*NBasis,XMuMat,NBasis,zero,Rkef,THCData%NCholErf*NBasis)
                  ! end if
                  ! call real_vw_x(this, Rkab(:, 1, 2), Rkab(:, 6, 15), THCData%NCholErf)
                  ! ! call real_vw_x(that0, Rkab(:, 1, 2), Rkef(:, 6, 9), THCData%NCholErf)
                  ! print*, this
                  ! call real_vw_x(this, Rkab(:, 1, 2), Rkab(:, 6, 14), THCData%NCholErf)
                  ! print*, this
                  ! call real_vw_x(this, Rkab(:, 1, 2), Rkab(:, 6, 16), THCData%NCholErf)
                  ! print*, this
                  ! STOP
                  ! ! this = zero
                  ! that0 = zero
                  ! do t = 1, THCData%NCholErf                        
                  !       this = this + Rkab(t, 1, 2) * Rkab(t, 6, 9)
                  !       do p = 1, nbasis
                  !             that0 = that0 + Rkab(t, 1, 2) * Rkab(t, 6, 9)*xmumat(p, 9)
                  !       end do
                  ! end do

                        

                  ! that0 = zero
                  ! do p = 1, nbasis
                  !       sniez = zero
                  !       do t = 1, THCData%NCholErf
                  !             sniez = sniez + Rkab(t, 1, 2) * Rkab(t, 6, p)
                  !       end do
                  !       !that0 = that0 + Rkab(t, 1, 2) * Rkab(t, 6, p)*xmumat(p, 9)
                  !       write(*, '(4I5, 3F20.15)')0, 1, 5, p-1, sniez, xmumat(p, 9), that0
                  !       !print*, p, sniez, xmumat(p, 9), that0
                  !       that0 = that0 + sniez*xmumat(p, 9)
                  ! end do
                  ! print*, 'that', that0
                  ! stop
                  
                  ! this0 = this
                  !  do t = 1, nbasis
                  !        this0 = this0 + xmumat(t,  9)
                  !  end do
                  !  print*, 't1', this0
                  !  this0 = this
                  !  do t = 1, nbasis
                  !        this0 = this0 + xmumat( t, 6)
                  !  end do
                  !  print*, 't2', this0
                  !  this0 = this
                  !  do t = 1, nbasis
                  !        this0 = this0 + xmumat(t, 2)
                  !  end do
                  !  print*, 't3', this0
                  !  this0 = this
                  !  do t = 1, nbasis
                  !        this0 = this0 + xmumat(t, 1)
                  !  end do
                  !  print*, 't4', this0


                  ! print*, 'thisssss', this, that0
                  ! stop
                  ! do p = 1, 1
                  !       do r = 1, 1
                  !             do q = 1, 1
                  !                   do s = 1, nbasis
                  !                         call real_vw_x(this, Rkab(:, p, r), Rkab(:, q, s), THCData%NCholErf)
                  !                         call real_vw_x(thisfr, Rkcd(:, p, r), Rkcd(:, q, s), THCData%NChol)
                  !                         call real_vw_x(thistr, Rkab(:, p, r), Rkef(:, q, s), THCData%NCholErf)
                  !                         call real_vw_x(thistr1, Rkef(:, p, r), Rkab(:, q, s), THCData%NCholErf)
                  !                         ! do t = 1, nbasis
                  !                         !       this = this	+ xmumat(s, t)
                  !                         ! end do

                  !                         if (abs(thisfr).gt.1.d-5)then
                  !                               write(*, '(4I5, 5F20.15)')p-1, r-1, q-1, s-1, thisfr!, thisfr, thisfr-this, thistr, thistr1
                  !                         end if
                  !                   end do
                  !             end do
                  !       end do
                  ! end do
                  ! stop

                  
                  if (Flags%JOBTYPE == JOB_TYPE_MP2 .or. Flags%JOBTYPE == JOB_TYPE_SRMP2) then                    
                        etot = zero
                        do i = 1, NI
                              do j = 1, NI
                                    !  J = (ii|jj)
                                    call real_vw_x(thisJ, Rkab(:, i, i), Rkab(:, j, j), THCData%NChol)
                                    !  K = (ij|ji)
                                    call real_vw_x(thisK, Rkab(:, i, j), Rkab(:, j, i), THCData%NChol)
                                    etot = etot + (two * thisJ - thisK)
                              end do
                        end do

                        do i = 1, NI
                              do t = NI + 1, NI + NA
                                    ! J = (ii|tt)
                                    call real_vw_x(thisJ, Rkab(:, i, i), Rkab(:, t, t), THCData%NChol)

                                    ! K = (it|ti)
                                    call real_vw_x(thisK, Rkab(:, i, t), Rkab(:, t, i), THCData%NChol)

                                    etot = etot + (two * thisJ - thisK)
                              end do
                        end do


                        do t = NI + 1, NI + NA
                              do u = t, NI + NA
                                    !  J = (tt|uu)
                                    call real_vw_x(thisJ, Rkab(:, t, t), Rkab(:, u, u), THCData%NChol)

                                    ! K = (tu|ut)
                                    call real_vw_x(thisK, Rkab(:, t, u), Rkab(:, u, t), THCData%NChol)

                                    if (t == u) then
                                          etot = etot + (thisJ - thisK)
                                    else
                                          etot = etot + two * (thisJ - thisK)
                                    end if
                              end do
                        end do

                        print*, 'etot0', etot0
                        print*, 'etot1', etot
                        print*, 'etot w/o enuc', etot+etot0
                        print*, 'etot', etot+etot0+AuxData%enuc
                        

                        print*, 'erohf from pyscf', AuxData%EROHF                        
                  else
                  
                        
                        etot = zero
!                         do p = 1, NI+NA
!                               do q = 1, NI+NA
!                                     do r = 1, NI+NA
!                                           do s = 1, NI+NA
!                                                 if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
!                       val = (AuxData%rdm2_pp(p-NI, q-NI, r-NI, s-NI) + AuxData%rdm2_pm(p-NI, q-NI, r-NI, s-NI))
!                       if( abs(val).gt.1.d-8)then
!                             print*, val,  AuxData%rdm2_full(p-NI, q-NI, r-NI, s-NI)
!                       end if
!                 end if
!           end do
!     end do
! end do
! end do
! stop
                        do p = 1, NI+NA
                              do q = 1, NI+NA
                                    do r = 1, NI+NA
                                          do s = 1, NI+NA
                                                call real_vw_x(this, Rkab(:, p, r), Rkab(:, q, s), THCData%NChol)
                                                !write(*, '(4I5, F20.15)')p, r, q, s, this

                                                val = zero
                                                if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
           val = val+  (AuxData%rdm2_pp(p-NI, q-NI, r-NI, s-NI) + AuxData%rdm2_pm(p-NI, q-NI, r-NI, s-NI))
!           val = val+  frac12 * AuxData%rdm2_full(p-NI, q-NI, r-NI, s-NI) 
                                  else
                                                      if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                            val =  val + two * Occ(p) * Occ(q)
                                                      end if

                                                      if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                            val = val +  -Occ(p) * Occ(q)
                                                      end if
                                                end if
                                                etot = etot + val * this
               !                                  if (abs(val * this).gt.1.d-1)then
               ! write(*, '(A20, 4I5, 4F15.7)')'energy-inter', p, q, r, s, val, this, val*this, etot
               !                                  end if
                                          end do
                                    end do
                              end do
                        end do

                        print*, 'etoat0', etot0
                        print*, 'etot-twoel', etot
                        print*, 'etot w/o enuc', etot+etot0
                        print*, 'etot', etot+etot0+AuxData%enuc

                        Write(6,'(1X,''CASSCF Energy (one-electron) THC'',X,F15.8)')etot0
                        Write(6,'(1X,''CASSCF Energy (w/o ENuc-q) THC'',X,F15.8)')etot+etot0
                        Write(6,'(1X,''Total CASSCF Energy from THC '',5X,F15.8)')etot + etot0 + AuxData%enuc
                        print*, 'ecas from pyscf', AuxData%Ecas
                  end if
!                  stop
                end associate

          end if

    end subroutine THC_init2

    subroutine calc_J_SR(CAONO, THCData, AuxData, AOBasis, System)
          use Cholesky_Gammcor
          use THC_Gammcor
          use THCFock
          use OneElectronInts_Gammcor
          use basis_sets
          use sys_definitions
          use gammcor_integrals
          type(TAOBASIS) :: AObasis
          type(TSystem) :: System


          double precision, dimension(:,:), intent(inout) :: CAONO
          type(TACppData), intent(in) :: AuxData
          type(TTHCData), intent(inout) :: THCData

          double precision, dimension(:,:, :), allocatable :: Cpo, Cpq
          double precision, dimension(:,:, :), allocatable :: J_LR, J_FR
          double precision, dimension(:,:), allocatable :: Zgh, ZghErf
          double precision :: Nk
          integer, dimension(2) :: NOcc
          integer k


          associate(Zgk=>THCData%Zgk, ZgkErf=>THCData%ZgkErf, Xga=>THCData%Xga, &
                XgaErf=>THCData%XgaErf, NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NBasis=>AuxData%NBasis, &
                Occ=>AuxData%Occ, NTHC=>THCData%NTHC, NTHCErf=>THCData%NTHCErf, ExternalOrdering=>THCData%ExternalOrdering)



            allocate(Cpq(NBasis, NBasis, 1))

            allocate(Cpo(NBasis, NIA, 1))
            Cpo = zero
            call auto2e_interface_C(Cpq(:, :, 1), CAONO, AOBasis, ExternalOrdering)            
            do k = 1, NIA
                  Nk = max(ZERO, Occ(k))
                  Cpo(k, k, 1) = Sqrt(Nk) 
            end do

          allocate(J_LR(NBasis, NBasis, 1))
          allocate(J_FR(NBasis, NBasis, 1))

          NOcc(1) = NIA
          Nocc(2) = 0
          allocate(Zgh(NTHC, NTHC))
          allocate(ZghErf(NTHCErf, NTHCErf))

          call real_abT(Zgh, Zgk, Zgk)
          call real_abT(ZghErf, ZgkErf, ZgkErf)
          
          call thc_Fock_JK(J_LR, Cpo, ZghErf, XgaErf, Nocc, .true., .false., one)
          call thc_Fock_JK(J_FR, Cpo, Zgh, Xga, Nocc, .true., .false., one)
          
          THCData%J_SR = J_FR(:,:,1) -J_LR(:,:,1)
          
        end associate
    end subroutine calc_J_SR

    subroutine canonicalize(CAONO, fij, fvw, AuxData, THCData, Xgp, AObasis, System, ext)
          use Cholesky_Gammcor
          use THC_Gammcor
          use OneElectronInts_Gammcor
          use basis_sets
          use sys_definitions
          use gammcor_integrals

          double precision, dimension(:,:), intent(inout) :: CAONO
          double precision, dimension(:), intent(out) :: fij, fvw            
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(in) :: THCData
          type(TAOBASIS) :: AObasis
          type(TSystem) :: System
          logical, intent(in),optional :: ext


          double precision, dimension(:,:), allocatable :: Cpi_extao, Cpv_extao, Cpa_empty
          double precision, dimension(:,:), allocatable :: Fockij, Fockoo, Fockvw
          double precision, dimension(:,:), allocatable :: H0oo, work
          double precision, dimension(:,:), intent(in) :: Xgp
          double precision, dimension(:,:), allocatable :: Xgt
          double precision, dimension(:,:,:), allocatable :: Rktu
          integer :: i, j
          integer :: t, u, v, w, NChol, NGridTHC
          double precision :: ECASSCF, ECumul, val, this
          logical :: have_spinres, have_rdm

          associate(Zgk=>THCData%Zgk, ExternalOrdering=>THCData%ExternalOrdering)

            ! print*, 'order2', THCData%ExternalOrdering
            ! print*, 'order3', ExternalOrdering
            allocate(Fockij(AuxData%NI, AuxData%NI))
            allocate(Fockoo(AuxData%NIA, AuxData%NIA))
            allocate(Fockvw(AuxData%NV, AuxData%NV))
            allocate(Cpi_extao(AuxData%nbasis, AuxData%NI))
            allocate(Cpv_extao(AuxData%nbasis, AuxData%NV))
            allocate(Cpa_empty(AuxData%nbasis, 0))

            ! call CalcMem(Fockij, 'Fockij')
            ! call CalcMem(Fockvw, 'Fockvw')
            ! call CalcMem(Cpi_extao, 'Cpi_extao')
            ! call CalcMem(Cpv_extao, 'Cpv_extao')

            call print_section('THC Fock matrix')

            if (ext == .false.)then
                  ! print*, 'ext false'
                  call thc_gammcor_F(Fockoo, Fockvw, CAONO(:, 1:AuxData%NIA),&
                        Cpa_empty, &
                        CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
                        AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering)
            else
                  ! print*, 'ext true'
                  ! print*, 'size AuxData%HNO0', size(AuxData%HNO0, dim=1), size(AuxData%HNO0, dim=2)
                  call thc_gammcor_F(Fockoo, Fockvw, CAONO(:, 1:AuxData%NIA),&
                        Cpa_empty, &
                        CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
                        AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering, &
                        AuxData%HNO0)
            end if

            if (AuxData%NI>0)then
                  Fockij = Fockoo(1:AuxData%NI, 1:AuxData%NI)
            end if

            !
            ! Spin-resolved 2-RDM (AC0PP path) or only the full 2-RDM (ph/AC0 path)
            !
            have_spinres = allocated(AuxData%rdm2_pp).and.allocated(AuxData%rdm2_pm)
            have_rdm     = have_spinres.or.allocated(AuxData%rdm2_full)

            if (have_rdm)then
                  allocate(H0oo(AuxData%NIA, AuxData%NIA))
                  if (ext == .false.)then
                        H0oo = AuxData%HNO0(1:AuxData%NIA, 1:AuxData%NIA)
                  else
                        allocate(work(AuxData%NBasis, AuxData%NIA))
                        call real_ab(work, AuxData%HNO0, CAONO(:, 1:AuxData%NIA))
                        call real_atb(H0oo, CAONO(:, 1:AuxData%NIA), work)
                  end if

                  ! E = ENuc + Sum(p) Occ(p) * (H0(p,p) + F(p,p))
                  ECASSCF = AuxData%enuc
                  do i = 1, AuxData%NIA
                        ECASSCF = ECASSCF + AuxData%Occ(i) * (H0oo(i,i) + Fockoo(i,i))
                  end do
            end if

            ! call thc_gammcor_F(Fockij, Fockvw, CAONO(:, 1:AuxData%NI),&
            !       CAONO(:, AuxData%NI+1:AuxData%NIA), &
            !       CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
            !       AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering)

!            print*, 'AuxData%Occ', AuxData%Occ

!            print*, 'fiij', AuxData%NI
!            print*, 'Fockij(1,1)', Fockij(1,1)
            ! do i = 1, AuxData%NI
            !       do i = j, AuxData%NI
            !             write(*,'(2I5, F15.6)') Fockij(i,j)
            !       end do
            ! end do
            ! print*, ''
            if (AuxData%NI>0)then
                  call symmetric_eigenproblem(fij, Fockij, AuxData%NI, .true.)
                  !do i = 1, AuxData%NI
                  !      print*, 'fij', i, fij(i)
                  !end do
            end if
            ! print*, fij
            ! print*, 'fvvww', AuxData%NV
            call symmetric_eigenproblem(fvw, Fockvw, AuxData%NV, .true.)
            !do i = 1, AuxData%NV
            !      print*, 'fvw', i, fvw(i)
            !end do

            ! print*, fvw

            if (AuxData%NI>0)then
                  call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
            end if
            call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)

            if (AuxData%NI>0)then
                  CAONO(:, 1:AuxData%NI)=Cpi_extao
            end if
            CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao

            !----------------------------------------------------------------------------
            ! Active-space cumulant correction to the Fock energy.
            !
            !   E = ENuc + Sum(p) Occ(p)*(H0(p,p) + F(p,p))
            !            + Sum(tuvw) [ Gamma(t,u,v,w) - Gamma_SD(t,u,v,w) ] * (tv|uw)
            !
            !   Gamma_SD(t,u,v,w) = 2*n(t)*n(u)*d(t,v)*d(u,w) - n(t)*n(u)*d(t,w)*d(u,v)
            !
            ! thc_gammcor_F builds F from the occupations only, so ECASSCF above already
            ! contains Gamma_SD everywhere. For CASSCF the cumulant is nonzero ONLY when
            ! all four indices are active, so only the NA**4 block of two-electron
            ! integrals is needed - never the full NBasis**4 Rkab.
            !
            ! Gamma is taken as rdm2_pp + rdm2_pm (AC0PP path), or as rdm2_full/2 when
            ! only the full 2-RDM was loaded (ph/AC0 path, only_full = .true.), since
            ! rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp = 2*(rdm2_pp + rdm2_pm).
            !----------------------------------------------------------------------------
            if (have_rdm)then

                  NGridTHC = size(Xgp, dim=1)
                  NChol    = size(Zgk, dim=2)

                  allocate(Xgt(NGridTHC, AuxData%NA))
                  allocate(Rktu(NChol, AuxData%NA, AuxData%NA))

                  call thc_gammcor_Xga(Xgt, Xgp, CAONO(:, AuxData%NI+1:AuxData%NIA), &
                        AOBasis, ExternalOrdering)
                  call thc_gammcor_Rkab_2(Rktu, Xgt, Xgt, Zgk, AuxData%NA, AuxData%NA, &
                        NChol, NGridTHC)

                  ECumul = zero
                  do t = 1, AuxData%NA
                        do u = 1, AuxData%NA
                              do v = 1, AuxData%NA
                                    do w = 1, AuxData%NA

                                          if (have_spinres)then
                                                val = AuxData%rdm2_pp(t, u, v, w) &
                                                    + AuxData%rdm2_pm(t, u, v, w)
                                          else
                                                val = frac12 * AuxData%rdm2_full(t, u, v, w)
                                          end if

                                          if (t==v.and.u==w) val = val - two * &
                                                AuxData%Occ(AuxData%NI+t) * AuxData%Occ(AuxData%NI+u)

                                          if (t==w.and.u==v) val = val + &
                                                AuxData%Occ(AuxData%NI+t) * AuxData%Occ(AuxData%NI+u)

                                          if (abs(val) > 1.d-12)then
                                                call real_vw_x(this, Rktu(:, t, v), Rktu(:, u, w), NChol)
                                                ECumul = ECumul + val * this
                                          end if

                                    end do
                              end do
                        end do
                  end do

                  deallocate(Xgt, Rktu)

                  !Write(6,'(1X,''CASSCF Energy from Fock (no cumulant)'',X,F15.8)') ECASSCF
                  !Write(6,'(1X,''Active-space cumulant correction'',5X,F15.8)') ECumul
                  AuxData%ECAS_THC = ECASSCF + ECumul
                  call print_energy1('Total CASSCF Energy from Fock', AuxData%ECAS_THC)
            end if
            
          end associate
    end subroutine canonicalize


    subroutine canonicalize_incore(CAONO, rdm1, AuxData, natural, TwoEl)
          
          double precision, dimension(:,:), intent(in) :: rdm1
          double precision, dimension(:,:), intent(inout) :: CAONO
          double precision, dimension(:), allocatable :: fij, fvw            
          type(TACppData), intent(in) :: AuxData
          double precision, dimension(:), intent(in) :: TwoEl
          integer, intent(in) :: natural
          integer :: i, j, k, l, ii, jj
          integer, external :: NAddr3 

          double precision, dimension(:,:), allocatable :: Cpi_extao, Cpv_extao
          double precision, dimension(:,:), allocatable :: Fockij, Fockvw

          associate(NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)
            allocate(Fockij(NI, NI))
            allocate(Fockvw(NV, NV))
            allocate(fij(NI))
            allocate(fvw(NV))
            allocate(Cpi_extao(AuxData%nbasis, NI))
            allocate(Cpv_extao(AuxData%nbasis, NV))

            if (NI > 0)then
                  ! inactive block
                  do i = 1, NI
                        do j = 1, NI
                              Fockij(i, j) = AuxData%HNO0(i, j)
                              do k = 1, NI
                                    Fockij(i, j) = Fockij(i, j) + AuxData%Occ(k)*&
                                          (Two*TwoEl(NAddr3(i, j, k, k))-TwoEl(NAddr3(i, k, j, k)))
                              end do
                              do k = NI+1, NIA
                                    do l = NI+1, NIA
                                          if (natural == 0)then
                                                Fockij(i, j) = Fockij(i, j) * rdm1(k-NI, l-NI)*&
                                                      (Two*TwoEl(NAddr3(i, j, k, l))-TwoEl(NAddr3(i, l, j, k)))
                                          else
                                                Fockij(i, j) = Fockij(i, j) * AuxData%Occ(k)*&
                                                      (Two*TwoEl(NAddr3(i, j, k, l))-TwoEl(NAddr3(i, l, j, k)))                                                
                                          end if
                                    end do
                              end do
                        end do
                  end do
                  ! virtual block
                  do i = 1, NV
                        do j = 1, NV
                              ii = i + NIA
                              jj = j + NIA
                              Fockvw(i, j) = AuxData%HNO0(ii, jj)
                              do k = 1, NI
                                    Fockvw(i, j) = Fockvw(i, j) + AuxData%Occ(k)*&
                                          (Two*TwoEl(NAddr3(ii, jj, k, k))-TwoEl(NAddr3(ii, k, jj, k)))
                              end do
                              do k = NI+1, NIA
                                    do l = NI+1, NIA
                                          if (natural == 0)then
                                                Fockvw(i, j) = Fockvw(i, j) + rdm1(k-NI, l-NI)* &
                                                       (Two*TwoEl(NAddr3(ii, jj, k, l))-TwoEl(NAddr3(ii, l, jj, k)))
                                          else
                                                Fockvw(i, j) = Fockvw(i, j) + AuxData%Occ(k)*&
                                                (Two*TwoEl(NAddr3(ii, jj, k, l))-TwoEl(NAddr3(ii, l, jj, k)))
                                          end if
                                    end do
                              end do
                        end do
                  end do

            end if
            if (AuxData%NI>0)then
                  call symmetric_eigenproblem(fij, Fockij, AuxData%NI, .true.)
                  do i = 1, AuxData%NI
                        print*, 'fij', i, fij(i)
                  end do
            end if


            call symmetric_eigenproblem(fvw, Fockvw, AuxData%NV, .true.)
                        do i = 1, AuxData%NV
                  print*, 'fvw', i, fvw(i)
            end do

            ! print*, fvw

            if (AuxData%NI>0)then
                  call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
            end if


            call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)
            if (AuxData%NI>0)then
                  CAONO(:, 1:AuxData%NI)=Cpi_extao
            end if

            CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao
          end associate
    end subroutine canonicalize_incore


    
    subroutine write_rdm2_dat(rdm2, NA)

        double precision, dimension(NA, NA, NA, NA), intent(in) :: rdm2
        integer, intent(in) :: NA
        integer :: i, j, k, l, unit
        double precision :: val

        open(unit=30, file='rdm2.dat', status='unknown', action='write')
!        print*, 'inside write'
        do l = 1, NA
            do k = 1, NA
                do j = 1, NA
                    do i = 1, NA
                        val = rdm2(i, j, k, l)
                        if (abs(val) > 1.0d-8) then
                              ! write(30, '(I4,1X,I4,1X,I4,1X,I4,1X,F19.12)') &
                              !       j, i, l, k, val
                              write(30, '(I4,1X,I4,1X,I4,1X,I4,1X,F19.12)') &
                                    j, l, i, k, val

                              ! write(*, '(I4,1X,I4,1X,I4,1X,I4,1X,F19.12)') &
                              !       j, i, l, k, val
                        end if
                    end do
                end do
            end do
        end do

        close(30)
    end subroutine write_rdm2_dat          

    
    subroutine check_energy_incore_spinres(AuxData, Ints)

          type(TACppData), intent(inout) :: AuxData
          type(TInts), intent(in) :: Ints
          double precision :: eone, etwo
          integer :: i, j
          integer :: p, q, r, s
          double precision :: this
          integer(I8) :: idx
          double precision :: val
          double precision, dimension(:,:), allocatable :: DM1a, DM1b
          
          associate(NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, &
                NBasis=>AuxData%NBasis, NIA=>AuxData%NIA, Occ=>AuxData%Occ, &
                IndAux=>AuxData%IndAux, &
                ints1e_aa=>Ints%ints1e_aa, ints1e_bb=>Ints%ints1e_bb, &
                ints2e_aa=>Ints%ints2e_aa, ints2e_bb=>Ints%ints2e_bb, ints2e_ab=>Ints%ints2e_ab)

            allocate(DM1a(NBasis, NBasis))
            allocate(DM1b(NBasis, NBasis))
            DM1a = zero
            DM1b = zero

            do i = 1, NI
                  DM1a(i, i) = one
                  DM1b(i, i) = one
            end do

            do i = 1, NA
                  do j = 1, NA
                        DM1a(NI+i, NI+j) = AuxData%rdm1_p(i, j)
                        DM1b(NI+i, NI+j) = AuxData%rdm1_m(i, j)
                  end do
            end do


            ! -----------------------------------
            ! 1. One-Electron Energy
            ! -----------------------------------
            eone = zero
            do p = 1, NIA
                  do q = 1, NIA                        
                        eone = eone + DM1a(p, q) * ints1e_aa(p, q)
                        eone = eone + DM1b(p, q) * ints1e_bb(p, q)
                  end do
            end do

            Write(6,'(1X,''CAS Energy (one-electron) calculated'',X,F15.8)') eone


            ! -----------------------------------
            ! 2. Two-Electron Energy
            ! -----------------------------------
            etwo = zero
            
            do p = 1, NIA
                  do q = 1, NIA
                        do r = 1, NIA
                              do s = 1, NIA
                                    
                                    
                                    idx = gmap_4fold(p, r, q, s, NBasis)
                                    
                                    ! Alpha-Alpha
                                    this = ints2e_aa(idx)
                                    val = zero
                                    if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                          val = AuxData%rdm2_pp(p-NI, q-NI, r-NI, s-NI)
                                    else
                                          val = DM1a(p, r) * DM1a(q, s) - DM1a(p, s) * DM1a(q, r)
                                    end if
                                    etwo = etwo + 0.5d0 * val * this

                                    ! Beta-Beta
                                    this = ints2e_bb(idx)
                                    val = zero
                                    if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                          val = AuxData%rdm2_mm(p-NI, q-NI, r-NI, s-NI)
                                    else
                                          val = DM1b(p, r) * DM1b(q, s) - DM1b(p, s) * DM1b(q, r)
                                    end if
                                    etwo = etwo + 0.5d0 * val * this

                                    ! Alpha-Beta
                                    this = ints2e_ab(idx)
                                    val = zero
                                    if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                          val = AuxData%rdm2_pm(p-NI, q-NI, r-NI, s-NI)
                                    else
                                          val = DM1a(p, r) * DM1b(q, s) 
                                    end if
                                    etwo = etwo + val * this

                              end do
                        end do
                  end do
            end do

            print*, 'AuxData%enuc', AuxData%enuc
            Write(6,'(1X,''CASSCF Two-electron Energy calculated'',X,F15.8)') etwo
            Write(6,'(1X,''CASSCF Total Energy calculated '',5X,F15.8)') eone + etwo + AuxData%enuc
            AuxData%ECas = eone + etwo + AuxData%enuc

          end associate
                
    end subroutine check_energy_incore_spinres

end module interface_pp
