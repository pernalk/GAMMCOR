module h5_reader
  use hdf5
  use h5lt
  implicit none
  private
  public :: HID_T
  public :: h5_init, h5_finish, h5_open, h5_close, h5_get, h5_attr, h5_exists, h5_size

  interface h5_get
     module procedure get_d1, get_d2, get_d4, get_i1
  end interface

  interface h5_attr
     module procedure attr_d, attr_s
  end interface

contains

  subroutine h5_init()
    integer :: ierr
    call h5open_f(ierr)
  end subroutine

  subroutine h5_finish()
    integer :: ierr
    call h5close_f(ierr)
  end subroutine

  function h5_open(filename) result(fid)
    character(len=*), intent(in) :: filename
    integer(HID_T) :: fid
    integer :: ierr
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, fid, ierr)
    if (ierr /= 0) error stop "h5_open: cannot open "//filename
  end function

  subroutine h5_close(fid)
    integer(HID_T), intent(in) :: fid
    integer :: ierr
    call h5fclose_f(fid, ierr)
  end subroutine

  logical function h5_exists(fid, path)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    integer :: ierr
    call h5ltpath_valid_f(fid, path, .true., h5_exists, ierr)
    if (ierr /= 0) h5_exists = .false.
  end function

  ! Total number of elements of a dataset -- lets the caller allocate before
  ! reading when the size is not known from the input file.
  integer function h5_size(fid, path)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    integer(HSIZE_T) :: fdims(7)
    integer(SIZE_T)  :: tsize
    integer :: tclass, ierr, i
    fdims = 1
    call h5ltget_dataset_info_f(fid, path, fdims, tclass, tsize, ierr)
    if (ierr /= 0) error stop "h5_size: missing dataset "//path
    h5_size = 1
    do i = 1, 7
       if (fdims(i) > 0) h5_size = h5_size * int(fdims(i))
    end do
  end function

  ! --- shape check shared by all ranks -------------------------------
  subroutine check_shape(fid, path, want)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    integer, intent(in) :: want(:)
    integer(HSIZE_T) :: fdims(7)
    integer(SIZE_T)  :: tsize
    integer :: tclass, ierr
    call h5ltget_dataset_info_f(fid, path, fdims, tclass, tsize, ierr)
    if (ierr /= 0) error stop "h5_get: missing dataset "//path
    if (any(fdims(1:size(want)) /= int(want, HSIZE_T))) &
         error stop "h5_get: shape mismatch for "//path
  end subroutine

  ! The dims argument of the h5lt readers is passed through a named local
  ! rather than as int(shape(a), HSIZE_T): an array expression forces the
  ! compiler to build a temporary, which trips -check arg_temp_created
  ! (ifx warning 406) and aborts the coarray debug build.

  subroutine get_d4(fid, path, a)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    double precision, intent(out) :: a(:,:,:,:)
    integer(HSIZE_T) :: dims(4)
    integer :: ierr
    call check_shape(fid, path, shape(a))
    dims = int(shape(a), HSIZE_T)
    call h5ltread_dataset_double_f(fid, path, a, dims, ierr)
    if (ierr /= 0) error stop "h5_get: read failed "//path
  end subroutine

  subroutine get_d2(fid, path, a)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    double precision, intent(out) :: a(:,:)
    integer(HSIZE_T) :: dims(2)
    integer :: ierr
    call check_shape(fid, path, shape(a))
    dims = int(shape(a), HSIZE_T)
    call h5ltread_dataset_double_f(fid, path, a, dims, ierr)
    if (ierr /= 0) error stop "h5_get: read failed "//path
  end subroutine

  subroutine get_d1(fid, path, a)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    double precision, intent(out) :: a(:)
    integer(HSIZE_T) :: dims(1)
    integer :: ierr
    call check_shape(fid, path, shape(a))
    dims = int(shape(a), HSIZE_T)
    call h5ltread_dataset_double_f(fid, path, a, dims, ierr)
    if (ierr /= 0) error stop "h5_get: read failed "//path
  end subroutine

  subroutine get_i1(fid, path, a)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path
    integer, intent(out) :: a(:)
    integer(HSIZE_T) :: dims(1)
    integer :: ierr
    call check_shape(fid, path, shape(a))
    dims = int(shape(a), HSIZE_T)
    call h5ltread_dataset_int_f(fid, path, a, dims, ierr)
    if (ierr /= 0) error stop "h5_get: read failed "//path
  end subroutine

  subroutine attr_d(fid, path, name, val)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path, name
    double precision, intent(out) :: val
    double precision :: buf(1)
    integer :: ierr
    call h5ltget_attribute_double_f(fid, path, name, buf, ierr)
    if (ierr /= 0) error stop "h5_attr: missing "//trim(name)//" on "//path
    val = buf(1)
  end subroutine

  subroutine attr_s(fid, path, name, val)
    integer(HID_T), intent(in) :: fid
    character(len=*), intent(in) :: path, name
    character(len=*), intent(out) :: val
    integer :: ierr
    val = ""
    call h5ltget_attribute_string_f(fid, path, name, val, ierr)
    if (ierr /= 0) error stop "h5_attr: missing "//trim(name)//" on "//path
  end subroutine

end module
