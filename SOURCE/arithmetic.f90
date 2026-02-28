module arithmetic
      use iso_fortran_env
      
      implicit none
      !
      ! 32-bit IEEE floating point kind
      ! --------------------------------
      !
      integer, parameter :: F32 = selected_real_kind(p=5)
      real(F32), parameter :: F32_ZERO = 0.0_F32
      !
      ! 64-bit IEEE floating point kind
      ! --------------------------------
      !
      integer, parameter :: F64 = selected_real_kind(p=15, r=307)
      real(F64), parameter :: F64_ZERO = 0.0_F64
      !
      ! 128-bit IEEE floating point kind
      !
      integer, parameter :: F128 = selected_real_kind(p=33, r=4931)
      !
      ! 32-bit signed integer kind
      ! ----------------------------
      !
      integer, parameter :: I32 = selected_int_kind(9)
      integer(I32), parameter :: I32_ZERO = 0_I32
      !
      ! 64-bit signed integer kind
      ! ----------------------------
      !
      integer, parameter :: I64 = selected_int_kind(18)
      integer(I64), parameter :: I64_ZERO = 0_I64

      real(F64), parameter :: ZERO = 0.0_F64
      real(F64), parameter :: ONE = 1.0_F64
      real(F64), parameter :: TWO = 2.0_F64
      
      ! -------------------------------------------------------------------
      ! FIXED-WIDTH PRECISION-PRESERVING FORMATS OF FLOATING POINT NUMBERS
      !              ES{FXY_ES_W}.{FXY_ES_D}E{FXY_ES_E}
      ! -------------------------------------------------------------------
      integer, parameter :: F32_ES_D = precision(F32_ZERO) + 1
      integer, parameter :: F32_ES_E = ceiling(log10(real(range(F32_ZERO), F32)))
      !
      ! The width of a character string representing F32 number
      ! in the ES format. An extra digit of precision is printed to
      ! prevent an information loss from rounding.
      !
      integer, parameter :: F32_ES_W = &
            1 + F32_ES_D + 1 &   ! +/- sign + significant digits + one extra digit
            + 1 + 1 + F32_ES_E & ! "E" symbol + sign of the exponent + exponent digits
            + 1                  ! decimal separator

      integer, parameter :: F64_ES_D = precision(F64_ZERO) + 1
      integer, parameter :: F64_ES_E = ceiling(log10(real(range(F64_ZERO), F64)))
      !
      ! The width of a character string representing F64 number
      ! in the ES format. An extra digit of precision is printed to
      ! prevent an information loss from rounding.
      !
      integer, parameter :: F64_ES_W = &
            1 + F64_ES_D + 1 &   ! +/- sign + significant digits + one extra digit
            + 1 + 1 + F64_ES_E & ! "E" symbol + sign of the exponent + exponent digits
            + 1                  ! decimal separator
      !
      ! Maximum width of a character string
      ! representing a variable of I32 kind
      ! in decimal form and without leading
      ! zeros (I0 format). The sign is included.
      !
      integer, parameter :: I32_MAXWIDTH = ceiling(log10(real(huge(I32_ZERO), F64))) + 1
      !
      ! Maximum width of a character string
      ! representing a variable of I64 kind
      ! in decimal form and without leading
      ! zeros (I0 format). The sign is included.
      !
      integer, parameter :: I64_MAXWIDTH = ceiling(log10(real(huge(I64_ZERO), F64))) + 1
end module arithmetic
