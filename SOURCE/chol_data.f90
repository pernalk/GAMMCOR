module chol_data
  implicit none
  double precision, allocatable, target :: CholVecsFF(:,:)
  integer :: NCholesky_stored = 0
end module chol_data
