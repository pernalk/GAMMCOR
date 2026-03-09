!==============================================================================
! batch_dgemm.F90
!
! Batched DGEMM wrapper for GAMMCOR.
! Three backends selected by preprocessor:
!   USE_MKL_BATCH  -> MKL dgemm_batch_strided (best for Intel MKL)
!   USE_CUDA_BATCH -> cuBLAS cublasDgemmStridedBatched (GPU, future)
!   (neither)      -> OpenMP loop over standard dgemm (any BLAS)
!
! Memory layout: flat 1D arrays with explicit strides.
! Matrix i in batch starts at offset i * stride from array base.
! stride=0 means the same matrix is reused for all batch elements.
!==============================================================================
module batch_dgemm

  implicit none
  private

  public :: batch_dgemm_strided

#ifdef USE_CUDA_BATCH
  public :: gpu_init
  public :: gpu_finalize
  public :: gpu_batch_dgemm_strided
#endif

#ifdef USE_CUDA_BATCH
  use, intrinsic :: iso_c_binding

  integer(c_int), parameter :: CUBLAS_OP_N = 0
  integer(c_int), parameter :: CUBLAS_OP_T = 1
  integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS = 0
  integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
  integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

  type(c_ptr), save :: cublas_handle = c_null_ptr
  logical, save     :: gpu_initialized = .false.

  interface
    integer(c_int) function cublasCreate(handle) &
        bind(C, name='cublasCreate_v2')
      import :: c_ptr, c_int
      type(c_ptr), intent(out) :: handle
    end function

    integer(c_int) function cublasDestroy(handle) &
        bind(C, name='cublasDestroy_v2')
      import :: c_ptr, c_int
      type(c_ptr), intent(in), value :: handle
    end function

    integer(c_int) function cublasDgemmStridedBatched_c( &
        handle, transa, transb, m, n, k,               &
        alpha, A, lda, strideA,                         &
        B, ldb, strideB,                                &
        beta, C, ldc, strideC,                          &
        batchCount) &
        bind(C, name='cublasDgemmStridedBatched')
      import :: c_ptr, c_int, c_double, c_long_long
      type(c_ptr), intent(in), value          :: handle
      integer(c_int), intent(in), value       :: transa, transb
      integer(c_int), intent(in), value       :: m, n, k
      type(c_ptr), intent(in), value          :: alpha
      type(c_ptr), intent(in), value          :: A
      integer(c_int), intent(in), value       :: lda
      integer(c_long_long), intent(in), value :: strideA
      type(c_ptr), intent(in), value          :: B
      integer(c_int), intent(in), value       :: ldb
      integer(c_long_long), intent(in), value :: strideB
      type(c_ptr), intent(in), value          :: beta
      type(c_ptr), intent(in), value          :: C
      integer(c_int), intent(in), value       :: ldc
      integer(c_long_long), intent(in), value :: strideC
      integer(c_int), intent(in), value       :: batchCount
    end function

    integer(c_int) function cudaMalloc(devPtr, nbytes) &
        bind(C, name='cudaMalloc')
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), intent(out) :: devPtr
      integer(c_size_t), intent(in), value :: nbytes
    end function

    integer(c_int) function cudaFree(devPtr) &
        bind(C, name='cudaFree')
      import :: c_ptr, c_int
      type(c_ptr), intent(in), value :: devPtr
    end function

    integer(c_int) function cudaMemcpy(dst, src, nbytes, kind) &
        bind(C, name='cudaMemcpy')
      import :: c_ptr, c_int, c_size_t
      type(c_ptr), intent(in), value :: dst, src
      integer(c_size_t), intent(in), value :: nbytes
      integer(c_int), intent(in), value :: kind
    end function
  end interface
#endif

contains

!==============================================================================
! batch_dgemm_strided -- Main CPU entry point
!
! Performs batch_size operations:
!   C_i = alpha * op(A_i) * op(B_i) + beta * C_i,   i = 0..batch_size-1
!
! A_i starts at A(1 + i*strideA)  (strideA=0 => same A for all)
! B_i starts at B(1 + i*strideB)  (strideB=0 => same B for all)
! C_i starts at C(1 + i*strideC)
!==============================================================================
subroutine batch_dgemm_strided(transa, transb, m, n, k, &
                               alpha, A, lda, strideA,   &
                               B, ldb, strideB,           &
                               beta, C, ldc, strideC,     &
                               batch_size)
  implicit none
  character(1), intent(in)        :: transa, transb
  integer, intent(in)             :: m, n, k, lda, ldb, ldc
  integer, intent(in)             :: strideA, strideB, strideC
  double precision, intent(in)    :: alpha, beta
  double precision, intent(in)    :: A(*)
  double precision, intent(in)    :: B(*)
  double precision, intent(inout) :: C(*)
  integer, intent(in)             :: batch_size

#ifdef USE_MKL_BATCH
  call dgemm_batch_strided(transa, transb, m, n, k, &
       alpha, A, lda, strideA,                       &
       B, ldb, strideB,                               &
       beta, C, ldc, strideC,                         &
       batch_size)
#else
  call batch_dgemm_omp_loop(transa, transb, m, n, k, &
       alpha, A, lda, strideA,                         &
       B, ldb, strideB,                                &
       beta, C, ldc, strideC,                          &
       batch_size)
#endif

end subroutine batch_dgemm_strided


!==============================================================================
! batch_dgemm_omp_loop -- OpenMP loop over standard dgemm
! Portable fallback for non-MKL BLAS (OpenBLAS, BLIS, Netlib).
!==============================================================================
subroutine batch_dgemm_omp_loop(transa, transb, m, n, k, &
                                alpha, A, lda, strideA,   &
                                B, ldb, strideB,           &
                                beta, C, ldc, strideC,     &
                                batch_size)
  implicit none
  character(1), intent(in)        :: transa, transb
  integer, intent(in)             :: m, n, k, lda, ldb, ldc
  integer, intent(in)             :: strideA, strideB, strideC
  double precision, intent(in)    :: alpha, beta
  double precision, intent(in)    :: A(*)
  double precision, intent(in)    :: B(*)
  double precision, intent(inout) :: C(*)
  integer, intent(in)             :: batch_size

  integer :: i, offA, offB, offC

  !$omp parallel do private(i, offA, offB, offC) schedule(static)
  do i = 0, batch_size - 1
     offA = i * strideA + 1
     offB = i * strideB + 1
     offC = i * strideC + 1
     call dgemm(transa, transb, m, n, k, &
          alpha, A(offA), lda, B(offB), ldb, &
          beta, C(offC), ldc)
  end do
  !$omp end parallel do

end subroutine batch_dgemm_omp_loop


#ifdef USE_CUDA_BATCH
!==============================================================================
! gpu_init -- Create cuBLAS handle (call once at program start)
!==============================================================================
subroutine gpu_init()
  implicit none
  integer(c_int) :: stat

  if (gpu_initialized) return

  stat = cublasCreate(cublas_handle)
  if (stat /= CUBLAS_STATUS_SUCCESS) then
     write(6,'(a,i0)') 'ERROR: cublasCreate failed, status = ', stat
     stop 1
  end if
  gpu_initialized = .true.
  write(6,'(a)') 'GPU: cuBLAS handle created'

end subroutine gpu_init


!==============================================================================
! gpu_finalize -- Destroy cuBLAS handle (call at program end)
!==============================================================================
subroutine gpu_finalize()
  implicit none
  integer(c_int) :: stat

  if (.not. gpu_initialized) return

  stat = cublasDestroy(cublas_handle)
  cublas_handle = c_null_ptr
  gpu_initialized = .false.

end subroutine gpu_finalize


!==============================================================================
! gpu_batch_dgemm_strided -- GPU batched DGEMM via cuBLAS
!
! Copies host arrays to device, runs batched dgemm on GPU, copies back.
! TODO: persistent device buffers, CUDA streams, async H2D/D2H overlap
!==============================================================================
subroutine gpu_batch_dgemm_strided(transa, transb, m, n, k, &
                                    alpha, A, lda, strideA,   &
                                    B, ldb, strideB,           &
                                    beta, C, ldc, strideC,     &
                                    batch_size)
  implicit none
  character(1), intent(in)        :: transa, transb
  integer, intent(in)             :: m, n, k, lda, ldb, ldc
  integer, intent(in)             :: strideA, strideB, strideC
  double precision, intent(in)    :: alpha, beta
  double precision, intent(in), target    :: A(*)
  double precision, intent(in), target    :: B(*)
  double precision, intent(inout), target :: C(*)
  integer, intent(in)             :: batch_size

  type(c_ptr) :: d_A, d_B, d_C
  integer(c_size_t) :: sizeA, sizeB, sizeC
  integer(c_int) :: stat, cu_transa, cu_transb
  integer(c_long_long) :: cu_strideA, cu_strideB, cu_strideC
  double precision, target :: h_alpha, h_beta
  integer :: colsA, colsB

  if (.not. gpu_initialized) call gpu_init()

  ! Convert transpose flags
  if (transa == 'T' .or. transa == 't') then
     cu_transa = CUBLAS_OP_T
     colsA = m
  else
     cu_transa = CUBLAS_OP_N
     colsA = k
  end if
  if (transb == 'T' .or. transb == 't') then
     cu_transb = CUBLAS_OP_T
     colsB = n
  else
     cu_transb = CUBLAS_OP_N
     colsB = n
  end if

  ! Compute device allocation sizes
  if (strideA == 0) then
     sizeA = int(lda, c_size_t) * int(colsA, c_size_t) * 8_c_size_t
  else
     sizeA = (int(strideA, c_size_t) * int(batch_size - 1, c_size_t) &
             + int(lda, c_size_t) * int(colsA, c_size_t)) * 8_c_size_t
  end if

  if (strideB == 0) then
     sizeB = int(ldb, c_size_t) * int(colsB, c_size_t) * 8_c_size_t
  else
     sizeB = (int(strideB, c_size_t) * int(batch_size - 1, c_size_t) &
             + int(ldb, c_size_t) * int(colsB, c_size_t)) * 8_c_size_t
  end if

  sizeC = (int(strideC, c_size_t) * int(batch_size - 1, c_size_t) &
          + int(ldc, c_size_t) * int(n, c_size_t)) * 8_c_size_t

  ! Allocate device memory
  stat = cudaMalloc(d_A, sizeA)
  stat = cudaMalloc(d_B, sizeB)
  stat = cudaMalloc(d_C, sizeC)

  ! Host -> Device
  stat = cudaMemcpy(d_A, c_loc(A(1)), sizeA, cudaMemcpyHostToDevice)
  stat = cudaMemcpy(d_B, c_loc(B(1)), sizeB, cudaMemcpyHostToDevice)
  if (beta /= 0.0d0) then
     stat = cudaMemcpy(d_C, c_loc(C(1)), sizeC, cudaMemcpyHostToDevice)
  end if

  h_alpha = alpha
  h_beta  = beta
  cu_strideA = int(strideA, c_long_long)
  cu_strideB = int(strideB, c_long_long)
  cu_strideC = int(strideC, c_long_long)

  ! Execute batched DGEMM on GPU
  stat = cublasDgemmStridedBatched_c( &
       cublas_handle, cu_transa, cu_transb, &
       int(m, c_int), int(n, c_int), int(k, c_int), &
       c_loc(h_alpha), d_A, int(lda, c_int), cu_strideA, &
       d_B, int(ldb, c_int), cu_strideB, &
       c_loc(h_beta), d_C, int(ldc, c_int), cu_strideC, &
       int(batch_size, c_int))

  if (stat /= CUBLAS_STATUS_SUCCESS) then
     write(6,'(a,i0)') 'ERROR: cublasDgemmStridedBatched failed, status=', stat
     stop 1
  end if

  ! Device -> Host (result C only)
  stat = cudaMemcpy(c_loc(C(1)), d_C, sizeC, cudaMemcpyDeviceToHost)

  ! Free device memory
  stat = cudaFree(d_A)
  stat = cudaFree(d_B)
  stat = cudaFree(d_C)

end subroutine gpu_batch_dgemm_strided
#endif

end module batch_dgemm
