module ducc_simple

      use acpp_types
      use types
      use math_constants
      use real_linalg !from gammcor integrals                                                                                      
      use THC_Gammcor
      use arithmetic
      use sort  !from gammcor integrals                                                                                            

      implicit none

contains


      subroutine load_1rdm(rdm_file, rdm)

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
      end subroutine load_1rdm

      subroutine load_2rdm(rdm_file, rdm, with_header)
            character(len=*), intent(in) :: rdm_file
            real(F64), dimension(:,:,:,:), intent(out) :: rdm
            logical, intent(in) :: with_header

            integer :: unit
            integer :: header(3)
            integer :: dim_sq, dim, i, j, k, l
            double precision :: val
            integer :: n
            

            n = size(rdm,dim=1)
            print*, 'n', n
            unit = 10
            open(unit=unit, file=rdm_file, form='unformatted', &
                  access='stream', status='old')

            if (with_header == .true.)then
                  read(unit) header
                  dim_sq = header(1)
                  dim = int(sqrt(real(dim_sq, F64)))
                  if (n.ne.dim)then
                        print*, 'fcidump and rdm files incompatible, exiting'
                        print*, 'n in fcidump', n
                        print*, 'n in rdm', dim
                        stop
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

      end subroutine load_2rdm


      subroutine read_rdm12s(pth, suffix, rdm1, rdm1_uu, rdm2, rdm2_uu, rdm2_ud)
            character(len=*), intent(in) :: pth, suffix
            real(F64), dimension(:,:), intent(out) :: rdm1, rdm1_uu
            real(F64), dimension(:,:,:,:), intent(out) :: rdm2, rdm2_uu, rdm2_ud

            character(len=512) :: filename


            filename = trim(pth) // "G1" // trim(suffix) // ".bin"
            call load_1rdm(filename, rdm1)

            filename = trim(pth) // "G1_uu" // trim(suffix) // ".bin"
            call load_1rdm(filename, rdm1_uu)

            filename = trim(pth) // "G2" // trim(suffix) // ".bin"
            call load_2rdm(filename, rdm2, .true.)

            filename = trim(pth) // "G2_uu" // trim(suffix) // ".bin"
            call load_2rdm(filename, rdm2_uu, .false.)

            filename = trim(pth) // "G2_ud" // trim(suffix) // ".bin"
            call load_2rdm(filename, rdm2_ud, .false.)

      end subroutine read_rdm12s


      subroutine read_header(filename, norb, nelec, ms2)
            character(len=*), intent(in) :: filename
            integer, intent(out) :: norb, nelec, ms2
            integer :: unit, ios, tblock
            character(len=256) :: line
            integer ::  pos, i


            tblock = 1
            norb = -1
            nelec = -1
            ms2 = -1

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

      end subroutine read_header


      subroutine read_fcidump_integrals_unrestricted(filename, norb, nelec, ms2, ints1e_aa, &
            ints2e_aa, ints2e_ab, e_nucl)

            implicit none

            character(len=*), intent(in) :: filename
            integer, intent(in) :: norb, nelec, ms2
            real(F64), dimension(:,:), intent(out) :: ints1e_aa
            real(F64), dimension(:,:,:,:), intent(out) :: ints2e_aa, ints2e_ab
            real(F64), intent(out) :: e_nucl


            integer :: unit, ios, tblock
            character(len=256) :: line
            integer :: p, q, r, s, pos, i
            real(F64) :: val


            unit = 10

            tblock = 1
            e_nucl = 0.0d0

            open(unit=unit, file=filename, status="old", action="read", iostat=ios)

            ! skip 5 lines
            header_loop: do i = 1, 5
                  read(unit=unit, fmt="(a)", iostat=ios) line
            end do header_loop

            ints1e_aa = zero
            ints2e_aa = zero
            ints2e_ab = zero

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
                              e_nucl = val
                              exit read_loop
                        end if
                  end if

                  if (tblock == 1 .and. (r /= 0 .or. s /= 0)) then
                        ints2e_aa(p, q, r, s) = val
                  else if (tblock == 3 .and. (r /= 0 .or. s /= 0)) then
                        ints2e_ab(p, q, r, s) = val
                  else if (tblock == 4 .and. r == 0 .and. s == 0) then
                        ints1e_aa(p, q) = val
                  else if (tblock == 2 .or. tblock == 5) then
                        cycle read_loop
                  end if

            end do read_loop

            close(unit)

      end subroutine read_fcidump_integrals_unrestricted


      subroutine ducc_test()
            implicit none

            integer, parameter :: NI = 0
            integer:: NA 
!            integer, parameter :: NIA = NI + NA

            real(F64), dimension(:,:), allocatable :: ints1e_aa_final
            real(F64), dimension(:,:,:,:), allocatable :: ints2e_aa_final, ints2e_ab_final
            real(F64) :: enuc_final
            integer :: NBasis, Nel, ms2_final

            real(F64), dimension(:,:), allocatable :: ints1e_aa_ducc
            real(F64), dimension(:,:,:,:), allocatable :: ints2e_aa_ducc, ints2e_ab_ducc
            real(F64) :: enuc_ducc
            integer :: nelec_ducc, ms2_ducc

            real(F64), dimension(:,:), allocatable :: rdm1_ref, rdm1_aa_ref
            real(F64), dimension(:,:,:,:), allocatable :: rdm2_ref, rdm2_aa_ref, rdm2_ab_ref
            integer :: dim_ref, dim2_ref

            real(F64), dimension(:,:), allocatable :: rdm1_final, rdm1_aa_final
            real(F64), dimension(:,:,:,:), allocatable :: rdm2_final, rdm2_aa_final, rdm2_ab_final
            integer :: dim_final, dim2_final


            real(F64), dimension(:,:), allocatable :: rdm1_bb_temp
            real(F64), dimension(:,:,:,:), allocatable :: rdm2_bb_temp

            real(F64), dimension(:,:,:,:), allocatable :: temp_G

            integer :: p, q, r, s

            real(F64) :: e1_ref, e2_aa_ref, e2_ab_ref, e2_ref, total_dmrg_energy_ref
            real(F64) :: e1_final, e2_aa_final_calc, e2_ab_final_calc, e2_final, total_dmrg_energy_final
            real(F64) :: Expected_correl


            call read_header("FCIDUMP-full-uhf", NBasis, Nel, ms2_final)
            print*, 'reading header complited'

            allocate(ints1e_aa_final(NBasis, NBasis))
            allocate(ints2e_aa_final(NBasis, NBasis, NBasis, NBasis))
            allocate(ints2e_ab_final(NBasis, NBasis, NBasis, NBasis))


            print*, 'reading integrals full'
            call read_fcidump_integrals_unrestricted("FCIDUMP-full-uhf", NBasis, &
                  Nel, ms2_final, &
                  ints1e_aa_final, ints2e_aa_final, ints2e_ab_final, &
                  enuc_final)

            print*, 'reading integrals full completead'
            call read_header("FCIDUMP-ref-uhf", NA, nelec_ducc, ms2_ducc)
            print*, 'reading ducc header completed'
            
            allocate(ints1e_aa_ducc(NA, NA))
            allocate(ints2e_aa_ducc(NA, NA, NA, NA))
            allocate(ints2e_ab_ducc(NA, NA, NA, NA))

            print*, 'reading integrals ducc'
            call read_fcidump_integrals_unrestricted("FCIDUMP-ref-uhf", NA, &
                  nelec_ducc, ms2_ducc, &
                  ints1e_aa_ducc, ints2e_aa_ducc, ints2e_ab_ducc, &
                  enuc_ducc)

            print*, 'reading integrals ducc completed'
            print*, 'NA', NA

            allocate(rdm1_ref(NA, NA))
            allocate(rdm1_aa_ref(NA, NA))
            allocate(rdm2_ref(NA, NA, NA, NA))
            allocate(rdm2_aa_ref(NA, NA, NA, NA))
            allocate(rdm2_ab_ref(NA, NA, NA, NA))

            allocate(rdm1_final(NBasis, NBasis))
            allocate(rdm1_aa_final(NBasis, NBasis))
            allocate(rdm2_final(NBasis, NBasis, NBasis, NBasis))
            allocate(rdm2_aa_final(NBasis, NBasis, NBasis, NBasis))
            allocate(rdm2_ab_final(NBasis, NBasis, NBasis, NBasis))

            print*, 'reading rdm ref'
            
            call read_rdm12s("ref/", "", rdm1_ref, rdm1_aa_ref, rdm2_ref, rdm2_aa_ref, rdm2_ab_ref)


            print*, 'reding ref rdms completed'
            print*, 'reading rdm final'

            call read_rdm12s("final/", "", rdm1_final, rdm1_aa_final, rdm2_final, rdm2_aa_final, rdm2_ab_final)
                  

            print*, 'reading final rdms completed'
            print*, 'nbasis', nbasis, NA


            write(*,*) 'Calculate reference energy'
            e1_ref = two * sum(ints1e_aa_final(1:NA, 1:NA) * rdm1_aa_ref(1:NA, 1:NA))
            
            e2_aa_ref = zero
            do s = 1, NA
                  do r = 1, NA
                        do q = 1, NA
                              do p = 1, NA
                                    e2_aa_ref = e2_aa_ref + ints2e_aa_final(p,q,r,s) * rdm2_aa_ref(p,r,q,s)
                              end do
                        end do
                  end do
            end do
            
            e2_ab_ref = zero
            do s = 1, NA
                  do r = 1, NA
                        do q = 1, NA
                              do p = 1, NA
                                    e2_ab_ref = e2_ab_ref + ints2e_ab_final(p,q,r,s) * rdm2_ab_ref(p,r,q,s)
                              end do
                        end do
                  end do
            end do


            e2_ref = 0.5_F64 * e2_aa_ref + e2_ab_ref
            total_dmrg_energy_ref = enuc_final + e1_ref + e2_ref
            write(*,'(A,F20.12)') 'Reference energy', total_dmrg_energy_ref

            write(*,*) 'Testing Final integrals'
            e1_final = 2.0_F64 * sum(ints1e_aa_final * rdm1_aa_final)

            e2_aa_final_calc = 0.0_F64
            do s = 1, NBasis
                  do r = 1, NBasis
                        do q = 1, NBasis
                              do p = 1, NBasis
                                    e2_aa_final_calc = e2_aa_final_calc + ints2e_aa_final(p,q,r,s) * rdm2_aa_final(p,r,q,s)
                              end do
                        end do
                  end do
            end do

            e2_ab_final_calc = 0.0_F64
            do s = 1, NBasis
                  do r = 1, NBasis
                        do q = 1, NBasis
                              do p = 1, NBasis
                                    e2_ab_final_calc = e2_ab_final_calc + ints2e_ab_final(p,q,r,s) * rdm2_ab_final(p,r,q,s)
                              end do
                        end do
                  end do
            end do

            e2_final = 0.5_F64 * e2_aa_final_calc + e2_ab_final_calc
            total_dmrg_energy_final = enuc_final + e1_final + e2_final
            write(*,'(A,F20.12,A)') "Total DMRG energy for final integrals : ", &
                  total_dmrg_energy_final

            Expected_correl = total_dmrg_energy_final - total_dmrg_energy_ref
            write(*,'(A,F20.12)') 'Correlation should be', Expected_correl

            deallocate(ints1e_aa_final, ints2e_aa_final, ints2e_ab_final)
            deallocate(ints1e_aa_ducc, ints2e_aa_ducc, ints2e_ab_ducc)
            deallocate(rdm1_ref, rdm1_aa_ref, rdm2_ref, rdm2_aa_ref, rdm2_ab_ref)
            deallocate(rdm1_final, rdm1_aa_final, rdm2_final, rdm2_aa_final, rdm2_ab_final)
            print*, 'DUCC finished'
            stop

      end subroutine ducc_test






end module ducc_simple
