module io
      use arithmetic

      implicit none

      interface io_size_byte
            module procedure :: io_size_byte_rank1_F64
            module procedure :: io_size_byte_rank2_F64
            module procedure :: io_size_byte_rank3_F64
            module procedure :: io_size_byte_rank4_F64
      end interface io_size_byte

contains

      function io_size_byte_rank1_F64(a)
            !
            ! Size of an array A in bytes. The size in bytes
            ! is stored in an integer of the I64 kind to prevent
            ! an overflow for large arrays.
            !
            integer(I64)                        :: io_size_byte_rank1_F64
            real(F64), dimension(:), intent(in) :: a

            io_size_byte_rank1_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank1_F64


      function io_size_byte_rank2_F64(a)
            integer(I64)                           :: io_size_byte_rank2_F64
            real(F64), dimension(:, :), intent(in) :: a

            io_size_byte_rank2_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank2_F64


      function io_size_byte_rank3_F64(a)
            integer(I64)                              :: io_size_byte_rank3_F64
            real(F64), dimension(:, :, :), intent(in) :: a

            io_size_byte_rank3_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank3_F64


      function io_size_byte_rank4_F64(a)
            integer(I64)                                 :: io_size_byte_rank4_F64
            real(F64), dimension(:, :, :, :), intent(in) :: a

            io_size_byte_rank4_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank4_F64
end module io
