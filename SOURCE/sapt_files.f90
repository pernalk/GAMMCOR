module sapt_files
implicit none

character(:),allocatable :: rdm2_file

contains

subroutine set_rdm2_filename(monomer)
!character(len=*), intent(in) :: monomer
integer, intent(in) :: monomer

!    select case (trim(adjustl(monomer)))
    select case (monomer)
    case (1)
    !case ("A", "a")
        rdm2_file = "2RDMA"
    case (2)
    !case ("B", "b")
        rdm2_file = "2RDMB"
    case (3)
    !case ("AB", "ab", "Ab", "aB")
        rdm2_file = "2RDM"
    case default
        print *, "Error: Unrecognized monomer label: ", monomer
        stop 1
    end select

end subroutine set_rdm2_filename

end module sapt_files
