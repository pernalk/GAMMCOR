module interface_pp

      use acpp_types
      use math_constants
      use basis_sets
      use Cholesky_Gammcor
      use THC_Gammcor
      use OneElectronInts_Gammcor
      use sys_definitions
      use gammcor_integrals
      use real_linalg
      use sort
      use clock

      implicit none


contains


      subroutine read_PYSCF(THCData, AuxData, CAONO, Flags, TwoEl)
            

            type(TTHCData), intent(inout) :: THCData
            type(TACppData), intent(inout) :: AuxData
            double precision, dimension(:,:), allocatable, intent(out) :: CAONO
            type(FlagsData), intent(in) :: Flags
            double precision, optional, intent(in) :: TwoEl(:)

            double precision :: ECAS
            integer :: nbasis, ninactive, nactive, nvirtual
            character(len=128) :: line
            character(len=256) :: file_auxdata, file_rdm2_aaaa, file_rdm2_abab, file_rdm1, file_rdm2_full, filename
            integer :: unit
            integer :: i, j, k, l
            
            double precision, dimension(:), allocatable :: occ_temp
            double precision, dimension(:, :), allocatable :: G1, CAOMO
            double precision, dimension(:, :), allocatable :: work
            double precision, dimension(:,:,:,:), allocatable :: rdm2_full
            logical :: STATE_AV = .false.
            integer :: natural
            integer :: anst(2)
            logical :: existsNosf, existsSf
            integer :: iostat
            character(len=10) :: suffix
            integer :: x, y
            integer :: count_suffix_files
            logical :: file_exists
            character(len=80) :: separator
            
            type(TAOBASIS) :: AObasis
            type(TSystem)  :: System

            call GEOM_init(Flags, AuxData, AObasis, System)
            
            !-----------------------------------------------------------------------------------------------------
            ! Check if there is a state specified in input 
            !-----------------------------------------------------------------------------------------------------
            print*, 'AuxData%nst', AuxData%nst

            anst = 0

            if (AuxData%nst(1) < 0)then

                  inquire(file='auxdata.txt', exist=existsNosf, iostat=iostat)
                  if (iostat /= 0) then
                        existsNosf = .false.
                  end if

                  count_suffix_files = 0

                  do x = 1, 10
                        do y = 1, 10

                              write(filename, '("auxdata_", I0, ".", I0, ".txt")') x, y                              
                              inquire(file=trim(filename), exist=existsSf, iostat=iostat)
                              if (iostat /= 0) then
                                    existsSf = .false.
                              end if

                              if (existsSf) then
                                    anst(1) = x
                                    anst(2) = y
                                    count_suffix_files = count_suffix_files + 1
                              end if
                        end do
                  end do


                  if (count_suffix_files == 0) then
                        if (existsNosf == .false.) then
                              print*, ''
                              print*, '-------------------------------------------------------------------------'
                              print *, 'Warning: Neither auxdata.txt nor any auxdata_x.y.txt found, exiting.'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                              stop                              
                        else
                              write(suffix, '("")')
                              anst(1) = 1
                              anst(2) = 1
                              print*, ''
                              print*, '-------------------------------------------------------------------------'                              
                              print *, 'Assuming state ', anst(1), '.', anst(2)
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                        end if
                  else if (count_suffix_files > 1) then
                        if (existsNosf == .false.) then
                              ! no auxdata.txt and multiple other files
                              print*, ''
                              print*, '-------------------------------------------------------------------------'
                              print *, 'Warning: Multiple auxdata_x.y.txt files found. No auxdata.txt file found. State notspecified, exiting'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                              stop
                        else
                              ! auxdata.txt and multiple other files
                              print*, ''
                              print*, '-------------------------------------------------------------------------'
                              print *, 'Warning: Multiple auxdata_x.y.txt files found. Auxdata.txt file found. State notspecified.'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                              write(suffix, '("")')
                              anst(1) = 1
                              anst(2) = 1
                              print*, ''
                              print*, '-------------------------------------------------------------------------'

                              print *, 'Assuming state ', anst(1), '.', anst(2)
                              print *, 'Using auxdata.txt'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                        end if
                  else if (count_suffix_files == 1) then
                        if (existsNosf == .false.) then
                              write(suffix, '("_", I0,".",I0)') anst(1), anst(2)
                              print*, ''
                              print*, '-------------------------------------------------------------------------'                              
                              print *, 'Assuming state ', anst(1), '.', anst(2)
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                        else
                              print*, ''
                              print*, '-------------------------------------------------------------------------'
                              print *, 'Warning: Both auxdata.txt and one auxdata_x.y, exist. State notspecified.'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                              write(suffix, '("")')
                              anst(1) = 1
                              anst(2) = 1
                              print*, ''
                              print*, '-------------------------------------------------------------------------'

                              print *, 'Assuming state ', anst(1), '.', anst(2)
                              print *, 'Using auxdata.txt'
                              print*, '-------------------------------------------------------------------------'
                              print*, ''

                        end if
                  end if


            else if (AuxData%nst(1) > 0) then

                  inquire(file='auxdata.txt', exist=existsNosf, iostat=iostat)
                  if (iostat /= 0) then
                        existsNosf = .false.
                  end if

                  write(suffix, '("_", I0,".",I0)') AuxData%nst(1), AuxData%nst(2)

                  ! Check for "auxdata_x.y.txt"
                  write(filename, '("auxdata", A, ".txt")') trim(suffix)
                  inquire(file=trim(filename), exist=existsSf, iostat=iostat)
                  if (iostat /= 0) then
                        existsSf = .false.
                  end if

                  if (existsNosf .and. existsSf) then
                        print*, ''
                        print*, '-------------------------------------------------------------------------'
                        print *, 'Warning: Both auxdata.txt and', filename, ' specified. Exiting.'
                        print*, '-------------------------------------------------------------------------'
                        print*, ''

                        stop 
                  else if (.not. existsNosf .and. .not. existsSf) then
                        print*, ''
                        print*, '-------------------------------------------------------------------------'
                        print *, 'Warning: Neither auxdata.txt nor auxdata_', trim(suffix), ' found, exiting'
                        print*, '-------------------------------------------------------------------------'
                        print*, ''

                        stop
                  else
                        if (existsNosf)write(suffix, '("")')

                        anst(1) = AuxData%nst(1)
                        anst(2) = AuxData%nst(2)
                        print*, ''
                        print*, '-------------------------------------------------------------------------'
                        print *, 'State', anst(1), '. ', anst(2), 'requested.'
                  end if
            end if

            
            AuxData%PYSCF = 1

            
            print*, 'suffix is', suffix
            write(file_auxdata, '("auxdata", A, ".txt")') trim(suffix)
            print*, 'reading from'
            print*, 'reading from', file_auxdata
            unit = 10    
            open(unit=unit, file=file_auxdata, status='old', action='read')

            read(unit, *) AuxData%Nbasis
            read(unit, *) AuxData%NI
            read(unit, *) AuxData%NA
            read(unit, *) AuxData%NV
            read(unit, *) AuxData%ECAS
            read(unit, *) AuxData%ENuc
            read(unit, *) AuxData%NEL
            read(unit, *) natural
            close(unit)

            print *, 'NBasis:', AuxData%nbasis
            print *, 'NInactive:', AuxData%NI
            print *, 'NActive:', AuxData%NA
            print *, 'NVirtual:', AuxData%NV
            print *, 'CASSCF Energy:', AuxData%ECAS
            print *, 'Nuclear Energy:', AuxData%ENuc
            print *, 'Natural orbitals', natural, '(0 - not used, 1 - used)'


            write(file_rdm2_aaaa, '("rdm2_aaaa", A, ".bin")') trim(suffix)
            write(file_rdm2_abab, '("rdm2_abab", A, ".bin")') trim(suffix)
            if (natural == 0) then
                  write(file_rdm2_full, '("rdm2", A, ".bin")') trim(suffix)
                  write(file_rdm1, '("rdm1", A, ".bin")') trim(suffix)

                  print*, file_auxdata
                  print*, file_rdm2_aaaa
                  print*, file_rdm2_abab
                  print*, file_rdm1
                  print*, file_rdm2_full

            else
                  
                  print*, 'auxddata', file_auxdata
                  print*, 'file_rdm2_aaaa', file_rdm2_aaaa
                  print*, 'rdm2_aaaa', file_rdm2_abab

            end if


            AuxData%NIA = AuxData%NI + AuxData%NA

            associate(NI=>AuxData%NI, NA=>AuxData%NA, NV=>AuxData%NV, NBasis=>AuxData%NBasis, NIA=>AuxData%NIA)
              allocate(AuxData%HNO0(NBasis, NBasis))

              THCData%ExternalOrdering = ORBITAL_ORDERING_PYSCF
              allocate(THCData%fij(NI))
              allocate(THCData%fvw(NV))



            !----------------------------------------------------------------------------------------
            ! READ 1-el integrals.
            ! If natural = 1, HCore are in NO basis.
            ! If natural = 0, HCore are in AO basis.
            !----------------------------------------------------------------------------------------
            unit = 20
            open(unit=unit, file='HCore.bin', status='old', access='stream', form='unformatted')
            read(unit) AuxData%HNO0
            close(unit)
            print*, 'AuxData%HNO0(1,1)', AuxData%HNO0(1,1)
            print*, 'AuxData%HNO0(1,1)', AuxData%HNO0(1,2)


            if (natural == 0) then
                  !----------------------------------------------------------------------------------------
                  ! READ CAOMO. Natural orbitals not used
                  !----------------------------------------------------------------------------------------
                  allocate(CAOMO(nbasis, nbasis))
                  unit = 30
                  inquire(file='C.bin', exist=file_exists)
                  if (file_exists) then
                        open(unit=unit, file='C.bin', status='old', access='stream', form='unformatted')
                  else
                        open(unit=unit, file='CAONO.bin', status='old', access='stream', form='unformatted')
                  endif
                  read(unit) CAOMO                  
                  close(unit)
            else
                  !----------------------------------------------------------------------------------------
                  ! READ CAONO (already in natural orbitals)
                  !----------------------------------------------------------------------------------------

                  allocate(CAONO(nbasis, nbasis))
                  unit = 30
                  inquire(file='C.bin', exist=file_exists)
                  if (file_exists) then
                        open(unit=unit, file='C.bin', status='old', access='stream', form='unformatted')
                  else
                        open(unit=unit, file='CAONO.bin', status='old', access='stream', form='unformatted')
                  endif
                  read(unit) CAONO                  
                  close(unit)
            end if

            
            print *, 'NBasis:', AuxData%nbasis
            print *, 'NInactive:', AuxData%NI
            print *, 'NActive:', AuxData%NA
            print *, 'NVirtual:', AuxData%NV
            print *, 'CASSCF Energy:', AuxData%ECAS
            print *, 'Nuclear Energy:', AuxData%ENuc
            print *, 'Natural orbitals', natural, '(0 - not used, 1 - used)'

            allocate(AuxData%rdm2_pp(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
            allocate(AuxData%rdm2_pm(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
          

            unit = 20
            open(unit=unit, file=file_rdm2_aaaa, status='old', access='stream', form='unformatted')
            read(unit) AuxData%rdm2_pp
            close(unit)

            unit = 21
            open(unit=unit, file=file_rdm2_abab, status='old', access='stream', form='unformatted')
            read(unit) AuxData%rdm2_pm
            close(unit)
            
            if (natural == 0)then
                  allocate(G1(AuxData%NA, AuxData%NA))
                  unit = 20
                  open(unit=unit, file=file_rdm1, status='old', access='stream', form='unformatted')
                  read(unit) G1
                  close(unit)
                  G1 = G1 / Two
                                    
                  !----------------------------------------------------------------------------------------
                  ! If natural orbitals are used, there is no need to read rdm2_full and transform it.
                  ! rdm2.dat is already prepared in NO basis, ready to be read later in AC subroutines.
                  !----------------------------------------------------------------------------------------

                  print*, 'reading rdm2_full', 'from binary file'
                  allocate(rdm2_full(AuxData%NA, AuxData%NA, AuxData%NA, AuxData%NA))
                  unit = 20
                  open(unit=unit, file=file_rdm2_full, status='old', access='stream', form='unformatted')
                  read(unit) rdm2_full
                  close(unit)

                  ! print for debugging
                  ! do i = 1, NA
                  !       do j = 1, NA
                  !             do k = 1, NA
                  !                   do l = 1, NA
                  !                         if (abs(rdm2_full(i, j, k, l)).gt.1.d-5)then
                  !                               write(*,'(4I5, F20.15)') i, j, k, l, rdm2_full(i, j, k, l)
                  !                         end if
                  !                   end do
                  !             end do
                  !       end do
                  ! end do
                  ! print*, ''

                  allocate(CAONO(nbasis, nbasis))

                  if (present(TwoEl)) then
                        call Trans2NO(G1, THCData, AuxData, CAOMO, CAONO, Flags, AObasis, System, rdm2_full=rdm2_full, TwoEl=TwoEl)
                  else
                        call Trans2NO(G1, THCData, AuxData, CAOMO, CAONO, Flags, AObasis, System, rdm2_full=rdm2_full)
                        AuxData%HNO0 = AuxData%HNO0_THC
                  end if


                  ! print for debugging
                  ! print*, 'after trans rdm2'
                  ! do i = 1, NA
                  !       do j = 1, NA
                  !             do k = 1,	NA
                  !                   do l = 1, NA
                  !                         if (abs(rdm2_full(i, j, k, l)).gt.1.d-5)then
                  !                               write(*,'(4I5, F20.15)') i, j, k, l, rdm2_full(i, j, k, l)
                  !                         end if
                  !                   end do
                  !             end do
                  !       end do
                  ! end do
                  ! print*, ''

                  
                  call write_rdm2_dat(rdm2_full, NA)
                  
            else
                  
                  allocate(occ_temp(AuxData%NI+AuxData%NA))

                  unit = 20
                  inquire(file='rdm1.bin', exist=file_exists)
                  if (file_exists) then
                        open(unit=unit, file='rdm1.bin', status='old', access='stream', form='unformatted')
                  else
                        open(unit=unit, file='occ.bin', status='old', access='stream', form='unformatted')
                  endif
                  read(unit) occ_temp
                  close(unit)
                  allocate(AuxData%Occ(AuxData%NBasis))
                  AuxData%Occ = zero
                  AuxData%Occ(1: AuxData%NI+AuxData%NA) = occ_temp
                  print*, AuxData%Occ(1:AuxData%NI+AuxData%NA)
            
                  allocate(AuxData%IndAux(NBasis))
                  AuxData%IndAux = 2
                  AuxData%IndAux(1:NI) = 0
                  AuxData%IndAux(NI+1:NIA) = 1

                  if (Flags%ITwoEl > 1)then
                        call THC_init2(Flags, THCData, AuxData,  AObasis, System, CAONO)
                  else
                        allocate(work(NBasis, NBasis))
                        call real_ab(work, AuxData%HNO0, CAONO)
                        call real_atb(AuxData%HNO0, CAONO, work)
                  end if

            end if            

          end associate
            
    end subroutine read_PYSCF

    subroutine PYSCF_wrapper(THCData, NInte1, Ninte2, NBasis, ENuc, CAONO, XKin, TwoEl, Occ, NAc, NInAc, anSt, Flags)

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

          AuxData%nst = anSt
          print*, 'inside wrapper'
          print*, Flags%ITwoEl
          if (Flags%ITwoEl==1)then
                print*, 'halo'
                unit = 21
                open(unit=unit, file='TWOEl.bin', status='old', access='stream', form='unformatted')
                read(unit) TwoEl
                close(unit)
                call read_PYSCF(THCData, AuxData, CAONO_PYSCF, Flags, TwoEl)
          else
                call read_PYSCF(THCData, AuxData, CAONO_PYSCF, Flags)
          end if
          print*, 'po read'
          CAONO = CAONO_PYSCF
          ENuc = AuxData%ENuc
          Occ = AuxData%Occ
          NAc = AuxData%NA
          NInAc = AuxData%NI

          ! One-electron integrals are now transformed to NO
          do i = 1, Nbasis
                do j = 1, i
                      ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                      XKin(ab) = AuxData%HNO0(i,j)
                end do
          end do
          
    end subroutine PYSCF_wrapper



    subroutine read_ORCA(THCData, AuxData, Flags, CAONO)

          implicit none

          type(TTHCData), intent(inout) :: THCData
          type(TACppData), intent(out) :: AuxData
          type(FlagsData), intent(in) :: Flags
          double precision, dimension(:,:), allocatable, intent(out) :: CAONO
          integer :: NBasis
          integer :: NI, NA, NIA, NEL
          integer :: unit
          integer :: i, j, k, l, n
          double precision, dimension(:), allocatable :: occ_temp
          double precision, dimension(:,:), allocatable :: G1, G1sort, CMONO_NA
          double precision, dimension(:,:), allocatable :: CMONO, CMOAO, CAOMO, work1
          integer, dimension(:), allocatable :: dy

          type(TAOBASIS) :: AObasis
          type(TSystem)  :: System

          integer :: iostat
          double precision :: val, suma
          logical, parameter :: sort_ang_mom = .true.
          integer :: idx
          integer :: n1, n2, kind
          integer :: this_dim, dim1, dim2
          double precision :: Err
          logical :: IsDiag
          double precision :: tol = 1.d-9
          logical :: iex

          double precision, dimension(:), allocatable :: work
          
          AuxData%ORCA = 1

          !----------------------------------------------------------------------------
          ! Read Nbasis and Number of electrons from geometry and basis
          !----------------------------------------------------------------------------

          call GEOM_init(Flags, AuxData, AObasis, System)

          !----------------------------------------------------------------------------
          ! read G1 - rdm1 - spin_traced
          !----------------------------------------------------------------------------
          unit = 10

          open(unit, file='G1.bin', form='unformatted', access='stream', status='old')
          read(unit) NA, n2, kind

          this_dim = NA * (NA +1) / 2
          allocate(work(this_dim))
          read(unit) work
          close(10)
          print*, 'na', NA
          allocate(G1(NA, NA))
          G1 = zero

          idx = 1
          do i = 1, NA
                do j = 1, i
                      G1(i, j) = work(idx)/Two
                      if (i.ne.j)then
                            G1(j, i) = G1(i, j)
                      end if
                      idx = idx + 1
                end do
          end do
          deallocate(work)

          AuxData%NA = NA
          print*, 'AuxData%NA', AuxData%NA
          !----------------------------------------------------------------------------
          ! read G2_uu
          !----------------------------------------------------------------------------

          print* , 'read_g2'
          
          allocate(AuxData%rdm2_pp(NA, NA, NA, NA))
          allocate(AuxData%rdm2_pm(NA, NA, NA, NA))
          AuxData%rdm2_pp = zero
          AuxData%rdm2_pm = zero


          unit = 10
          open(unit, file='G2_uu.bin', form='unformatted', access='stream', status='old')
          
          do k = 1, NA
                do j = 1, NA
                      do l = 1, NA
                            do i = 1, NA
                                  read(10) val
                                  AuxData%rdm2_pp(l, k, i, j) = val
                                  AuxData%rdm2_pp(k, l, j, i) = val
                            end do
                      end do
                end do
          end do
          close(10)

           open(unit=10, file='G2_ud.bin', form='unformatted', access='stream', status='old')

          do k = 1, NA
                do j = 1, NA
                      do l = 1, NA
                            do i = 1, NA
                                  read(10) val
                                  AuxData%rdm2_pm(l, k, i, j) = val
                                  AuxData%rdm2_pm(k, l, j, i) = val
                            end do
                      end do
                end do
          end do
          close(10)

          !----------------------------------------------------------------------------
          ! CAONO
          !----------------------------------------------------------------------------
          print* , 'read_cao'
          allocate(CMOAO(NBasis, NBasis))
          allocate(CAOMO(NBasis, NBasis))
          unit = 10
          open(unit, file='C0.bin', form='unformatted', access='stream', status='old')
          read(unit) dim1, dim2, kind
          print*, dim1, dim2, 'dim1dim2'

          if (dim1.ne.dim2)then
                print*, 'N_row != N_col in DMRG. Linear dependencies in basis. exiting'
                stop
          end if
          read(unit) CMOAO
          close(10)

          CAOMO = transpose(CMOAO)
          deallocate(CMOAO)
          
          allocate(CAONO(NBasis, NBasis))

          print*, 'AuxData%NA2', AuxData%NA

          THCData%ExternalOrdering = ORBITAL_ORDERING_ORCA
          allocate(THCData%fij(AuxData%NI))
          allocate(THCData%fvw(AuxData%NV))


          call Trans2NO(G1, THCData, AuxData, CAOMO, CAONO, Flags, AObasis, System)

          ! allocate(occ_temp(NA))

          ! !----------------------------------------------------------------------------
          ! ! diagonalize rdm1
          ! !----------------------------------------------------------------------------
          ! call symmetric_eigenproblem(occ_temp, G1, NA, .true.)

          ! allocate(dy(NA))
          ! do i = 1, NA
          !       dy(i) = i
          ! end do

          ! !----------------------------------------------------------------------------
          ! ! sort egivals and eigvecs
          ! !----------------------------------------------------------------------------
          ! call dsort0(occ_temp, dy, NA, -2)

          ! allocate(CMONO_NA(NA, NA))
          
          ! CMONO_NA = zero

          ! do i = 1, NA
          !       CMONO_NA(:,i) = G1(:, dy(i))
          ! end do
          ! deallocate(G1)

          ! err = zero
          ! do i = 1, NA
          !       do j = 1, NA
          !             if (i==1)then
          !                   err = err  + abs(One-CMONO_NA(i,i))
          !             else
          !                   err = err + abs(CMONO_NA(i, j))
          !             end if
          !       end do
          ! end do
          ! print*, 'err', err

          

          ! !----------------------------------------------------------------------------
          ! ! check if if G1sort is diagonal - if yes: all data already in Natural Orbitals 
          ! !----------------------------------------------------------------------------
          
          ! IsDiag = .false.
          ! if (abs(err) < tol)then
          !       print*, 'Natural orbital basis, do not transform'
          !       IsDiag = .true.
          ! else
          !       print*, 'Transform to natural orbital basis'
          ! end if
          
          ! !----------------------------------------------------------------------------
          ! ! check for almost zero, or negative occupancies
          ! !----------------------------------------------------------------------------
          ! suma = zero
          ! NA = 0
          ! do i = 1, AuxData%NA
          !       occ_temp(i) = abs(occ_temp(i))
          !       suma = suma + occ_temp(i)
          !       if (occ_temp(i)>1.d-8)then
          !             NA = NA + 1
          !       end if
          ! end do

          ! AuxData%NI = (NEL -(suma*two)+1.d-2)/2
          ! NI = AuxData%NI

          ! AuxData%switch = 0
          ! if (NA.ne.AuxData%NA)then
          !       Write(6, '(1x, "WARNING! The number of partially occupied orbitals different from NActDMRG read from ORCA. Some active orbitals must be unoccupied.", /)')
          !       NA = AuxData%NA
          !       AuxData%switch = 1
          ! end if
          
          
          ! allocate(AuxData%Occ(AuxData%NBasis))
          ! AuxData%Occ = zero
          ! AuxData%Occ(1: AuxData%NI) = One
          ! AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp

          ! AuxData%NIA = AuxData%NI + AuxData%NA
          ! AuxData%NV = AuxData%NBasis - AuxData%NIA
          ! NIA = AuxData%NIA 

          ! !----------------------------------------------------------------------------
          ! ! print occupancies
          ! !----------------------------------------------------------------------------
          ! write(*, '(A10, A16)')'#DMRG',  'Occupancy'
          
          ! suma = zero
          ! do i = 1, AuxData%NIA
          !       write(*, '(I10, F16.6)') i, AuxData%Occ(i)
          !       suma = suma + AuxData%Occ(i)
          ! end do

          ! write(*,'(A30, F20.15)') "sum of occupancies", suma


          ! !----------------------------------------------------------------------------
          ! ! prepare IndAux
          ! !----------------------------------------------------------------------------

          ! print* , 'prepare indaux'
          
          ! allocate(AuxData%IndAux(NBasis))
          ! AuxData%IndAux = 2
          ! AuxData%IndAux(1:NI) = 0
          ! AuxData%IndAux(NI+1:NIA) = 1

          !----------------------------------------------------------------------------
          ! read G2_uu
          !----------------------------------------------------------------------------

          ! print* , 'read_g2'
          
          ! allocate(AuxData%rdm2_pp(NA, NA, NA, NA))
          ! allocate(AuxData%rdm2_pm(NA, NA, NA, NA))
          ! AuxData%rdm2_pp = zero
          ! AuxData%rdm2_pm = zero


          ! unit = 10
          ! open(unit, file='G2_uu.bin', form='unformatted', access='stream', status='old')
          
          ! do k = 1, NA
          !       do j = 1, NA
          !             do l = 1, NA
          !                   do i = 1, NA
          !                         read(10) val
          !                         AuxData%rdm2_pp(l, k, i, j) = val
          !                         AuxData%rdm2_pp(k, l, j, i) = val
          !                   end do
          !             end do
          !       end do
          ! end do
          ! close(10)

          !  open(unit=10, file='G2_ud.bin', form='unformatted', access='stream', status='old')

          ! do k = 1, NA
          !       do j = 1, NA
          !             do l = 1, NA
          !                   do i = 1, NA
          !                         read(10) val
          !                         AuxData%rdm2_pm(l, k, i, j) = val
          !                         AuxData%rdm2_pm(k, l, j, i) = val
          !                   end do
          !             end do
          !       end do
          ! end do
          ! close(10)

          ! call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
          ! call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)

          !----------------------------------------------------------------------------
          ! CAONO
          !----------------------------------------------------------------------------
          ! print* , 'read_cao'
          ! allocate(CMOAO(NBasis, NBasis))
          ! unit = 10
          ! open(unit, file='C0.bin', form='unformatted', access='stream', status='old')
          ! read(unit) dim1, dim2, kind
          ! print*, dim1, dim2, 'dim1dim2'

          ! if (dim1.ne.dim2)then
          !       print*, 'N_row != N_col in DMRG. Linear dependencies in basis. exiting'
          !       stop
          ! end if
          ! read(unit) CMOAO
          ! close(10)


          !----------------------------------------------------------------------------
          ! If IsDiag==true, CMOAO is already in natural orbitals MO=NO. Do not transform.
          ! Else, transform to NO
          !----------------------------------------------------------------------------


          ! if (IsDiag==.true.)then
          !       print*, 'Data from MolMPS is already in Natural Orbital basis'
          !       allocate(CAONO(NBasis, NBasis))
          !       CAONO = transpose(CMOAO)

          ! else
          !       print*, 'Data from MolMPS is not in Natural Orbital basis. Transforming.'


          !       allocate(CMONO(NBasis, NBasis))
          !       CMONO = zero
          !       do i = 1, NBasis
          !             CMONO(i,i) = One
          !       end do

          !       CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 
          !       print* , 'created mono'
                
          !       allocate(CAONO(NBasis, NBasis))

          !       print* , 'calc caono'
          !       call real_atb(CAONO, CMOAO, CMONO)

          ! end if

          !----------------------------------------------------------------------------
          ! One-electron-integrals
          !----------------------------------------------------------------------------
          ! print* , 'open fact'
          ! allocate(AuxData%HNO0(NBasis, NBasis))
          ! unit = 10
          ! open(unit,File='FACT.bin',form='unformatted',access='stream', status='old')

          ! read(unit)n, j, k
          ! print*, n, j, k
          ! this_dim = n * (n+1) /2
          ! allocate(work(this_dim))
          ! read(unit) work
          ! close(unit)
          
          ! if (n.ne.NBasis)then
          !       print*, 'dim of HNO0', n
          !       print*, 'nbasis', NBasis
          !       stop
          ! end if

          ! idx = 1
          ! do i = 1, NBasis
          !       do j = 1, i
          !             AuxData%HNO0(i, j) = work(idx)
          !             if (i.ne.j)then
          !                   AuxData%HNO0(j, i) = AuxData%HNO0(i, j)
          !             end if
          !             idx = idx + 1
	    !   end do
          ! end do
          ! deallocate(work)

          ! print*, 'AuxData%ENuc', AuxData%ENuc

          ! allocate(work1(NBasis, NBasis))
          ! print*, 'mult h01'
          ! call real_ab(work1, AuxData%HNO0, CMONO)
          ! print*, 'mult h02'
          ! call real_atb(AuxData%HNO0, CMONO, work1)
          !call real_atb(AuxData%HNO0, CAONO, work1)

                    
          
    end subroutine read_ORCA


    subroutine Trans2NO(G1, THCData, AuxData, CAOMO, CAONO, Flags, AObasis, System, rdm2_full, TwoEl)

          double precision, dimension(:,:), intent(inout) :: G1
          type(TACppData), intent(inout) :: AuxData
          type(TTHCData), intent(inout) :: THCData
          double precision, dimension(:,:), intent(in) :: CAOMO
          double precision, dimension(:,:), intent(out) :: CAONO
          type(FlagsData), intent(in) :: Flags
          type(TAOBASIS), intent(in) :: AObasis
          type(TSystem), intent(in)  :: System

          double precision, dimension(:), allocatable :: occ_temp
          integer, dimension(:), allocatable :: dy
          double precision, dimension(:,:), allocatable :: CMONO_NA, CMONO, work
          double precision, optional, intent(inout) :: rdm2_full(:,:,:,:)
          double precision, optional, intent(in) :: TwoEl(:)

          double precision :: err
          logical :: IsDiag
          double precision :: tol = 1.d-9
          double precision :: suma, ETot0
          integer :: i, j, NI, NA, NIA
          double precision, dimension(:,:), allocatable :: G1_pure


          print*, 'AuxData%NA-3', AuxData%NA
          NA = AuxData%NA
          associate(NBasis=>AuxData%NBasis)

            allocate(G1_pure(NBasis, NBasis))
            G1_pure = G1
            allocate(occ_temp(NA))
            !----------------------------------------------------------------------------
            ! diagonalize rdm1
            !----------------------------------------------------------------------------
            call symmetric_eigenproblem(occ_temp, G1, NA, .true.)

            allocate(dy(NA))
            do i = 1, NA
                  dy(i) = i
            end do

            !----------------------------------------------------------------------------
            ! sort egivals and eigvecs
            !----------------------------------------------------------------------------
            call dsort0(occ_temp, dy, NA, -2)

            print*, 'occ_temp', occ_temp
            allocate(CMONO_NA(NA, NA))
          
            CMONO_NA = zero

            do i = 1, NA
                  CMONO_NA(:,i) = G1(:, dy(i))
            end do
            
            
            err = zero
            do i = 1, NA
                  do j = 1, NA
                        if (i==1)then
                              err = err  + abs(One-CMONO_NA(i,i))
                        else
                              err = err + abs(CMONO_NA(i, j))
                        end if
                  end do
            end do
            print*, 'err', err



            !----------------------------------------------------------------------------
            ! check if if G1sort is diagonal - if yes: all data already in Natural Orbitals 
            !----------------------------------------------------------------------------
          
            IsDiag = .false.
            if (abs(err) < tol)then
                  print*, 'Natural orbital basis, do not transform'
                  IsDiag = .true.
            else
                  print*, 'Transform to natural orbital basis'
            end if
          
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

            AuxData%NI = (AuxData%NEL -(suma*two)+1.d-2)/2
            NI = AuxData%NI

            AuxData%switch = 0
            if (NA.ne.AuxData%NA)then
                  Write(6, '(1x, "WARNING! The number of partially occupied orbitals different from NActDMRG read from ORCA. Some active orbitals must be unoccupied.", /)')
                  NA = AuxData%NA
                  AuxData%switch = 1
            end if
          
          
            allocate(AuxData%Occ(AuxData%NBasis))
            print*,  'AuxData%NBasis', AuxData%NBasis
            AuxData%Occ = zero
            AuxData%Occ(1: AuxData%NI) = One
            AuxData%Occ(AuxData%NI+1:AuxData%NI+NA) = occ_temp
            print*, 'AuxData%NI', AuxData%NI
            print*, 'AuxData%NA', AuxData%NA, NA
            
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
            
            print* , 'prepare indaux'
            
            allocate(AuxData%IndAux(NBasis))
            AuxData%IndAux = 2
            AuxData%IndAux(1:NI) = 0
            AuxData%IndAux(NI+1:NIA) = 1

          
            call rdm2_MO_NO_trans(AuxData%rdm2_pp, CMONO_NA, NA)
            call rdm2_MO_NO_trans(AuxData%rdm2_pm, CMONO_NA, NA)

            if (present(rdm2_full)) then
                  call rdm2_MO_NO_trans(rdm2_full, CMONO_NA, NA)
            end if


            if (IsDiag==.true.)then
                  print*, 'Data is already in Natural Orbital basis'
                  CAONO = CAOMO

            else
                  print*, 'Data is not in Natural Orbital basis. Transforming.'
                  
                  allocate(CMONO(NBasis, NBasis))
                  CMONO = zero
                  do i = 1, NBasis
                        CMONO(i,i) = One
                  end do
                  
                  CMONO(NI+1:NIA, NI+1:NIA) = CMONO_NA 
                  
                  print* , 'calc caono'
                  call real_ab(CAONO, CAOMO, CMONO)

            end if


            !----------------------------------------------------------------------------
            ! Canonicalization inside THC_init
            ! On Exit CAONO is canonical in inactive/virtual blocks
            !----------------------------------------------------------------------------
            if (Flags%ITwoEl > 1)then
                  call THC_init2(Flags, THCData, AuxData,  AObasis, System, CAONO)
            end if
            
            allocate(work(NBasis, NBasis))
            call real_ab(work, AuxData%HNO0, CAONO)
            call real_atb(AuxData%HNO0, CAONO, work)

            ETot0 = zero
            do i = 1, NBasis
                  ETot0 = ETot0 + two* AuxData%Occ(i) * AuxData%HNO0(i,i)                  
                  !write(*, '(I5, 2F20.15)') i, AuxData%HNO0(i,i), AuxData%Occ(i)
            end do
            print*, 'etot1 z HNO_ext', ETot0
            
          end associate

    end subroutine Trans2NO
    

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

          character(:), allocatable :: basis_path, xyz_path, binPath
          logical, parameter :: sort_ang_mom = .true.
          integer :: units, nbasis
          
          units = Flags%IUnits
          basis_path = Flags%BasisSetPath //Flags%BasisSet

          print*, 'basis_path', basis_path
          xyz_path = "./input.inp"

          call auto2e_init()          
          call sys_Read_XYZ(System, xyz_path, Flags%IUnits)
          print*, 'NAtoms = ', System%NAtoms
          print*, 'NElectrons = ', System%NElectrons
          call basis_newAObasis(AObasis, System, basis_path, .True., sort_ang_mom)

          call sys_NuclearRepulsion(AuxData%ENuc,System)
          Nbasis = AObasis%NAOSpher
          AuxData%NBasis = AObasis%NAOSpher
          AuxData%NEL = System%NElectrons
          print*, 'NEL', AuxData%NEL
          print*, 'NBasis', AuxData%NBasis
          
    end subroutine GEOM_init

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
          integer :: i, j, NAO, ij, ab
          double precision, allocatable :: Xgp(:,:)
          integer :: NA, NI, NV, NIA
          integer :: units, nbasis
          type (tclock) :: timer
          double precision :: val1, val2, val3

          ! Fock matrix diag                                                                                                                                                                                                                                                                                                
          double precision, dimension(:,:), allocatable :: Fockij, Fockvw
          double precision, dimension(:,:), allocatable :: Cpi_extao, Cpa_extao, Cpv_extao
          double precision, dimension(:,:), allocatable :: work, H0_extao

          !---------------------testing energy                                                                                                                                                                                                                                                                              
          integer :: p, q, r, s, k, l, a, b, ii, nn
          double precision, allocatable :: Rkab(:,:,:), Rkcd(:,:,:)
          double precision :: this, ETot, etot0, val, val_this


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

          if (CholThr < zero.or.THCThr <zero)then
                print*, 'CholAccu', CholAccu
                call thc_gammcor_XZ(Xgp, THCData%Zgk, AOBasis, System, CholAccu)
          else
                print*, 'CholeskyThreshold', CholThr
                print*, 'THCThreshold', THCThr
                call thc_gammcor_XZ(Xgp, THCData%Zgk, AOBasis, System,CholAccu, CholThr, THCThr)
                print*, 'CholeskyThreshold', CholThr
                print*, 'THCThreshold', THCThr
          end if
          allocate(CAONO(nbasis, nbasis))
          nao = nbasis
          allocate(work(nao, nbasis))

          if (AuxData%PYSCF==1 .or. AuxData%ORCA==1)then
                print*, 'here', AuxData%NI, AuxData%NIA, AuxData%NBasis
                CAONO = CAONO_IN
                print*, CAONO(1,2), CAONO(2, 1)
          else
                CAONO = transpose(CAONO_IN)
                print*, CAONO(1,2), CAONO(2, 1)
          end if

          call canonicalize(CAONO, THCData%fij, THCData%fvw, AuxData, THCData, Xgp, AObasis, System)


          allocate(H0_extao(nao, nao))
          allocate(AuxData%HNO0_THC(nbasis, nbasis))

          call ints1e_gammcor_H0_extao(H0_extao, AObasis, System, THCData%ExternalOrdering)
          print*, 'h01', H0_extao(1,1)

          call real_ab(work, H0_extao, CAONO)
          call real_atb(AuxData%HNO0_THC, CAONO, work)

          THCData%NTHC=size(Xgp,dim=1)
          THCData%NChol=size(THCData%Zgk,dim=2)
          allocate(THCData%Xga(THCData%NTHC,NBasis))

          Call thc_gammcor_Xga(THCData%Xga, Xgp, CAONO,&
                AOBasis, THCData%ExternalOrdering)

          ETot0 = zero
          do i = 1, NIA
                ETot0 = ETot0 + two* AuxData%Occ(i) * AuxData%HNO0_THC(i,i)
                ! if (abs(AuxData%Occ(i) * AuxData%HNO0_THC(i,i)).gt.1.d-1)then
                !       write(*, '(A10, I5, 2F20.15)')'etot-1', i, AuxData%Occ(i), AuxData%HNO0_THC(i,i)
                ! end if
          end do
          print*, 'etot1 z HNO thc', ETot0

          AuxData%HNO0 = AuxData%HNO0_THC


          CAONO_IN = CAONO

          if (AuxData%OnlyEnergy == .true.)then
                allocate(Rkab(THCData%NChol,nao, nao))                                                                                                                                                                                     
                Call thc_gammcor_Rkab_2(Rkab, THCData%Xga, THCData%Xga, THCData%Zgk, NBasis, NBasis,&                                                                                                                                      
                      THCData%NChol, THCData%NTHC)
                associate(Occ=>AuxData%Occ, IndAux=>AuxData%IndAux)

                  ! do p = 1, 10
                  !       do q = 1, 10
                  !             do r = 1, 10
                  !                   do s = 1, 10
                  !                         call real_vw_x(this, Rkab(:, p, r), Rkab(:, q, s), THCData%NChol)
                  !                         if (abs(this).gt.1.d-5)then
                  !                               write(*, '(4I5, F20.15)')p, r, q, s, this
                  !                         end if
                  !                   end do
                  !             end do
                  !       end do
                  ! end do
                  
                  etot = zero
                  etot = etot0
                  do p = 1, NI+NA
                        do q = 1, NI+NA
                              do r = 1, NI+NA
                                    do s = 1, NI+NA
                                          call real_vw_x(this, Rkab(:, p, r), Rkab(:, q, s), THCData%NChol)
                                          !write(*, '(4I5, F20.15)')p, r, q, s, this

                                          val = zero
                                          if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                                                val = val+  (AuxData%rdm2_pp(p-NI, q-NI, r-NI, s-NI) + AuxData%rdm2_pm(p-NI, q-NI, r-NI, s-NI))
                                          else
                                                if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                      val =  val + two * Occ(p) * Occ(q)
                                                end if

                                                if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                                      val = val +  -Occ(p) * Occ(q)
                                                end if
                                          end if
                                          etot = etot + val * this
                                          if (abs(val * this).gt.1.d-5)then
                                                write(*, '(A20, 4I5, 3F20.15)')'kupsko', p, q, r, s, val, this, etot
                                          end if
                                    end do
                              end do
                        end do
                  end do
                  print*, 'etot0', etot0
                  print*, 'etot1', etot
                  print*, 'etot w/o enuc', etot+etot0
                  print*, 'etot', etot+etot0+AuxData%enuc
                  print*, 'ecas', AuxData%Ecas
                end associate
                stop
          end if

    end subroutine THC_init2

    subroutine canonicalize(CAONO, fij, fvw, AuxData, THCData, Xgp, AObasis, System)
          use Cholesky_Gammcor
          use THC_Gammcor
          use OneElectronInts_Gammcor
          use basis_sets
          use sys_definitions
          use gammcor_integrals

          double precision, dimension(:,:), intent(inout) :: CAONO
          double precision, dimension(:), intent(out) :: fij, fvw            
          type(TACppData), intent(in) :: AuxData
          type(TTHCData), intent(in) :: THCData
          type(TAOBASIS) :: AObasis
          type(TSystem) :: System


          double precision, dimension(:,:), allocatable :: Cpi_extao, Cpv_extao
          double precision, dimension(:,:), allocatable :: Fockij, Fockvw
          double precision, dimension(:,:), intent(in) :: Xgp

          associate(Zgk=>THCData%Zgk, ExternalOrdering=>THCData%ExternalOrdering)

            print*, 'order2', THCData%ExternalOrdering
            print*, 'order3', ExternalOrdering
            allocate(Fockij(AuxData%NI, AuxData%NI))
            allocate(Fockvw(AuxData%NV, AuxData%NV))
            allocate(Cpi_extao(AuxData%nbasis, AuxData%NI))
            allocate(Cpv_extao(AuxData%nbasis, AuxData%NV))

            call CalcMem(Fockij, 'Fockij')
            call CalcMem(Fockvw, 'Fockvw')
            call CalcMem(Cpi_extao, 'Cpi_extao')
            call CalcMem(Cpv_extao, 'Cpv_extao')


            call thc_gammcor_F(Fockij, Fockvw, CAONO(:, 1:AuxData%NI),&
                  CAONO(:, AuxData%NI+1:AuxData%NIA), &
                  CAONO(:, AuxData%NIA+1:AuxData%NBasis), &
                  AuxData%Occ(1:AuxData%NIA), Zgk, Xgp, AOBasis, System, ExternalOrdering)

            call symmetric_eigenproblem(fij, Fockij, AuxData%NI, .true.)
            call symmetric_eigenproblem(fvw, Fockvw, AuxData%NV, .true.)


            call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
            call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)

            CAONO(:, 1:AuxData%NI)=Cpi_extao
            CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao
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
            call symmetric_eigenproblem(fij, Fockij, AuxData%NI, .true.)
            call symmetric_eigenproblem(fvw, Fockvw, AuxData%NV, .true.)

            call real_ab(Cpi_extao, CAONO(:, 1:AuxData%NI), Fockij)
            call real_ab(Cpv_extao, CAONO(:, AuxData%NIA+1:AuxData%NBasis), Fockvw)

            CAONO(:, 1:AuxData%NI)=Cpi_extao
            CAONO(:, AuxData%NIA+1:AuxData%NBasis) = Cpv_extao
          end associate
    end subroutine canonicalize_incore


    
    subroutine write_rdm2_dat(rdm2, NA)

        double precision, dimension(NA, NA, NA, NA), intent(in) :: rdm2
        integer, intent(in) :: NA
        integer :: i, j, k, l, unit
        double precision :: val

        open(unit=30, file='rdm2.dat', status='unknown', action='write')
        print*, 'inside write'
        do l = 1, NA
            do k = 1, NA
                do j = 1, NA
                    do i = 1, NA
                        val = rdm2(i, j, k, l)
                        if (abs(val) > 1.0d-8) then
                            write(30, '(I4,1X,I4,1X,I4,1X,I4,1X,F19.12)') &
                                  j, i, l, k, val
                            write(*, '(I4,1X,I4,1X,I4,1X,I4,1X,F19.12)') &
                                  j, i, l, k, val
                        end if
                    end do
                end do
            end do
        end do

        close(30)
    end subroutine write_rdm2_dat          

    
end module interface_pp
