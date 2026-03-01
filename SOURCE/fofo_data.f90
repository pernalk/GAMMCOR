module fofo_data
  implicit none
  double precision, allocatable, target :: IntsFOFO(:,:)  ! (nCD, nAB)
  double precision, allocatable, target :: IntsFFOO(:,:)  ! (nCD, nAB)
  integer :: IFOFO_ram = 0   ! 1=in-memory, 0=disk (default)
end module fofo_data
