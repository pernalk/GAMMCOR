module io
      use arithmetic
      use display

      implicit none

      interface io_size_byte
            module procedure :: io_size_byte_rank1_F64
            module procedure :: io_size_byte_rank2_F64
            module procedure :: io_size_byte_rank3_F64
            module procedure :: io_size_byte_rank4_F64
      end interface io_size_byte
      
      !
      ! The maximum length of an I/O error message
      !
      integer, parameter :: IO_MAX_MSGLEN = 256

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


      function io_text_open(filename, s)
            !
            ! Open a text file.
            !
            integer                            :: io_text_open
            character(*), intent(in)           :: filename
            character(*), intent(in)           :: s

            integer :: open_stat
            character(len=IO_MAX_MSGLEN) :: errmsg

            open(newunit=io_text_open, file=filename, status=s, &
                  access="SEQUENTIAL", iostat=open_stat, iomsg=errmsg)
            
            if (open_stat .ne. 0) then
                  call msg("Could not open file")
                  call msg(trim(adjustl(filename)))
                  call msg(trim(errmsg))
                  error stop
            end if
      end function io_text_open


      subroutine io_text_readline(line, u, eof)
            !
            ! Read a line from a text file. The limit for the line
            ! size is MAXCHUNKS * DEFLEN characters (see the code).
            !
            character(:), allocatable, intent(out) :: line
            integer, intent(in)                    :: u
            logical, optional, intent(out)         :: eof

            character(len=80) :: chunk
            character(len=IO_MAX_MSGLEN) :: errmsg
            integer :: s, ios
            integer :: n
            integer, parameter :: maxchunks = 2**10

            line = ""
            if (present(eof)) eof = .false.
            
            lineloop: do n = 1, maxchunks
                  read(u, "(A)", advance="NO", size=s, &
                        iostat=ios, iomsg=errmsg) chunk

                  if (s > 0) then
                        line = line // chunk(1:s)
                  end if

                  if (ios == iostat_end) then
                        if (present(eof)) eof = .true.
                        exit lineloop
                  else if (ios == iostat_eor) then
                        exit lineloop
                  else if (ios .ne. 0) then
                        call msg("COULD NOT READ LINE")
                        call msg(trim(errmsg))
                        stop
                  end if
            end do lineloop
      end subroutine io_text_readline
end module io
