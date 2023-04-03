! -------------------------------------------------------------------------------------
!                  Cholesky decomposition of the Coulomb matrix
! -------------------------------------------------------------------------------------
! version 15.04.2021 (Marcin Modrzejewski)
! -------------------------------------------------------------------------------------
! Compute the R factors of the Cholesky decomposition:
!
! V = R**T * R
!
! where
!
! V(pq,rs) = (pq|rs)
!
! and
!
! R(k, pq), k = 1, ... NCholesky, pq=1, ..., NOrbPairs
!
! The Cholesky vectors (rows of R) are generated until one of the following
! conditions is satisfied:
!
! (1) the algorithm has generated a user-defined maximum number of Cholesky vectors;
! (2) a user-defined convergence criterion is satisfied.
!
! The R matrices contain only permutationally unique pq's.
!
! The implementation is based on the algorithm of Ref. 1 with some modifications,
! most importantly, the target decomposition error is Tr(V-R**T*R) as proposed
! by Harbrecht et al. in Ref. 2.
!
! 1. Aquilante et al. Subsection 13.5, Cholesky Decomposition Techniques in Electronic
!    Structure Theory in Linear-Scaling Techniques in Computational Chemistry and
!    Physics: Methods and Applications, 301-343, Springer 2011;
!    doi: 10.1007/978-90-481-2853-2_13
! 2. Harbrecht, H., Peters, M., and Schneider R. Appl. Num. Math. 62, 428 (2012);
!    doi: 10.1016/j.apnum.2011.10.001
!
module Cholesky
      use arithmetic
      use io
      use clockMM
      use string
      use sorter_Cholesky
      use display
      use sort
      
      implicit none

      integer, parameter :: CHOL_ACCURACY_DEFAULT = 1
      integer, parameter :: CHOL_ACCURACY_TIGHT = 2
      integer, parameter :: CHOL_ACCURACY_LUDICROUS = 3

      integer, parameter :: CHOL_TRANSF_BUFFER_DIM = 1000

      integer, parameter :: SIGNIFICANT_PAIRS = 1
      integer, parameter :: QUALIFIED_PAIRS = 2
      integer, parameter :: BASE_PAIRS = 3

      type TCholeskyVecs
            real(F64), dimension(:, :), allocatable :: R
            integer :: NCholesky
            integer :: MaxNCholesky
            integer :: NAO
            integer, dimension(3) :: NOrbPairs
      end type TCholeskyVecs

contains

      subroutine chol_av_x(w, a, lda, v, m, n, alpha, beta)
            !
            ! Perform matrix-vector multiplication w = alpha * Av + beta * w
            !
            real(F64), dimension(*), intent(inout)   :: w
            real(F64), dimension(lda, *), intent(in) :: a
            integer, intent(in)                      :: lda
            real(F64), dimension(*), intent(in)      :: v
            integer, intent(in)                      :: m
            integer, intent(in)                      :: n
            real(F64), intent(in)                    :: alpha
            real(F64), intent(in)                    :: beta

            external :: dgemv
            
            call dgemv("N", m, n, alpha, a, lda, v, 1, beta, w, 1)
      end subroutine chol_av_x
      

      subroutine chol_vwT_x(a, lda, v, w, m, n, alpha)
            !
            ! A <- alpha * v*w**T + A
            !
            real(F64), dimension(lda, *), intent(inout) :: a
            integer, intent(in)                         :: lda
            real(F64), dimension(*), intent(in)         :: v
            real(F64), dimension(*), intent(in)         :: w
            integer, intent(in)                         :: m
            integer, intent(in)                         :: n
            real(F64), intent(in)                       :: alpha

            external :: dger

            call dger(m, n, alpha, v, 1, w, 1, a, lda)
      end subroutine chol_vwT_x


      subroutine chol_abT_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:m, 1:k) B(1:n, 1:k)**T + beta*C(1:m, 1:n)
            !
            real(F64), dimension(ldc, *), intent(out) :: c
            integer, intent(in)                       :: ldc
            real(F64), dimension(lda, *), intent(in)  :: a
            integer, intent(in)                       :: lda
            real(F64), dimension(ldb, *), intent(in)  :: b
            integer, intent(in)                       :: ldb
            integer, intent(in)                       :: m
            integer, intent(in)                       :: n
            integer, intent(in)                       :: k
            real(F64), optional, intent(in)           :: alpha
            real(F64), optional, intent(in)           :: beta
            
            real(F64) :: alpha0, beta0
            external :: dgemm

            if (present(alpha)) then
                  alpha0 = alpha
            else
                  alpha0 = 1.0_F64
            end if

            if (present(beta)) then
                  beta0 = beta
            else
                  beta0 = 0.0_F64
            end if

            call dgemm("N", "C", m, n, k, alpha0, a, lda, b, ldb, beta0, c, ldc)
      end subroutine chol_abT_x


      subroutine chol_aTb_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:k, 1:m)**T B(1:k, 1:n) + beta*C(1:m, 1:n)
            !
            real(F64), dimension(ldc, *), intent(out) :: c
            integer, intent(in)                       :: ldc
            real(F64), dimension(lda, *), intent(in)  :: a
            integer, intent(in)                       :: lda
            real(F64), dimension(ldb, *), intent(in)  :: b
            integer, intent(in)                       :: ldb
            integer, intent(in)                       :: m
            integer, intent(in)                       :: n
            integer, intent(in)                       :: k
            real(F64), optional, intent(in)           :: alpha, beta
            
            real(F64), parameter :: alpha_default = 1.0_F64
            real(F64), parameter :: beta_default = 0.0_F64
            real(F64) :: alpha_param, beta_param
            external :: dgemm

            if (present(alpha)) then
                  alpha_param = alpha
            else
                  alpha_param = alpha_default
            end if
            if (present(beta)) then
                  beta_param = beta
            else
                  beta_param = beta_default
            end if
            call dgemm("C", "N", m, n, k, alpha_param, a, lda, b, ldb, beta_param, c, ldc)
      end subroutine chol_aTb_x


      subroutine chol_aTv_x(w, a, lda, v, m, n, alpha, beta)
            !
            ! Perform matrix-vector multiplication w = alpha * A**T v + beta * w
            !
            real(F64), dimension(*), intent(inout)   :: w
            real(F64), dimension(lda, *), intent(in) :: a
            integer, intent(in)                      :: lda
            real(F64), dimension(*), intent(in)      :: v
            integer, intent(in)                      :: m
            integer, intent(in)                      :: n
            real(F64), intent(in)                    :: alpha
            real(F64), intent(in)                    :: beta

            external :: dgemv
            
            call dgemv("T", m, n, alpha, a, lda, v, 1, beta, w, 1)
      end subroutine chol_aTv_x
      

      subroutine chol_MOTransf(S, Vecs, CA, a0, a1, CB, b0, b1)
            real(F64), dimension(:, :), intent(out) :: S
            type(TCholeskyVecs), intent(in)         :: Vecs
            real(F64), dimension(:, :), intent(in)  :: CA
            integer, intent(in)                     :: a0, a1
            real(F64), dimension(:, :), intent(in)  :: CB
            integer, intent(in)                     :: b0, b1

            integer :: p, q, pq, pq0, pq1
            integer :: NA, NB
            integer :: ldS, ldR
            integer :: BufferDim, NVecs
            real(F64), dimension(:, :), allocatable :: CAT, CBT
            real(F64), dimension(:, :), allocatable :: W

            associate ( &
                  NAO => Vecs%NAO, &
                  NOrbPairs => Vecs%NOrbPairs(BASE_PAIRS), &
                  NCholesky => Vecs%NCholesky, &
                  R => Vecs%R)
                  
                  NA = a1 - a0 + 1
                  NB = b1 - b0 + 1
                  ldS = size(S, dim=1)
                  ldR = size(R, dim=1)
                  BufferDim = CHOL_TRANSF_BUFFER_DIM
                  allocate(CAT(NA, NAO))
                  allocate(CBT(NB, NAO))
                  allocate(W(NA*NB, BufferDim))
                  CAT = transpose(CA(1:NAO, a0:a1))
                  CBT = transpose(CB(1:NAO, b0:b1))
                  S = ZERO
                  NVecs = 0
                  pq = 0
                  do q = 1, NAO
                        do p = 1, q
                              if (NVecs == 0) then
                                    W = ZERO
                              end if
                              NVecs = NVecs + 1
                              pq = pq + 1
                              if (p /= q) then
                                    call chol_vwT_x(W(:, NVecs), NA, CAT(:, p), CBT(:, q), NA, NB, ONE)
                                    call chol_vwT_x(W(:, NVecs), NA, CAT(:, q), CBT(:, p), NA, NB, ONE)
                              else
                                    call chol_vwT_x(W(:, NVecs), NA, CAT(:, p), CBT(:, q), NA, NB, ONE)
                              end if
                              if (pq == NOrbPairs .or. NVecs == BufferDim) then
                                    pq0 = pq - NVecs + 1
                                    pq1 = pq
                                    !
                                    ! S(1:NCholesky, 1:Na*Nb) <- S(1:NCholesky, 1:Na*Nb) + R(1:NCholesky, pq0:pq1)W(1:NA*NB, 1:NVecs)**T
                                    !
                                    call chol_abT_x(S, ldS, R(:, pq0:pq1), ldR, W, NA*NB, NCholesky, NA*NB, NVecs, &
                                                    alpha=ONE, beta=ONE)
                                    NVecs = 0
                              end if
                        end do
                  end do
            end associate
      end subroutine chol_MOTransf


      subroutine chol_MOTransf_TwoStep(S, Vecs, CA, a0, a1, CB, b0, b1, MaxBufferDimMB)
            !
            ! R(k, pq) -> R(k, ab) for a=a0...a1, b=b0...b1, k=1...NCholesky
            !
            ! MaxBufferDimMB Maximum buffer size in megabytes
            !
            real(F64), dimension(:, :), intent(out) :: S
            type(TCholeskyVecs), intent(in)         :: Vecs
            real(F64), dimension(:, :), intent(in)  :: CA
            integer, intent(in)                     :: a0, a1
            real(F64), dimension(:, :), intent(in)  :: CB
            integer, intent(in)                     :: b0, b1
            integer, intent(in)                     :: MaxBufferDimMB

            integer(I64) :: MaxBufferDim
            
            associate ( &
                  NAO => Vecs%NAO, &
                  NCholesky => Vecs%NCholesky, &
                  R => Vecs%R)

                  MaxBufferDim = (int(MaxBufferDimMB,I64)*1024*1024)/(storage_size(S)/8)
                  call chol_MOTransf_TwoStep_x(S, R, CA, CB, a0, a1, b0, b1, &
                        NAO, NCholesky, MaxBufferDim)
            end associate
      end subroutine chol_MOTransf_TwoStep

      
      subroutine chol_MOTransf_TwoStep_x(Rab, Rpq, CA, CB, a0, a1, b0, b1, NAO, NCholesky, MaxBufferDim)
            real(F64), dimension(:, :), intent(out)                    :: Rab
            real(F64), dimension(:, :), intent(in)                     :: Rpq
            real(F64), dimension(:, :), intent(in)                     :: CA
            real(F64), dimension(:, :), intent(in)                     :: CB
            integer, intent(in)                                        :: a0, a1
            integer, intent(in)                                        :: b0, b1
            integer, intent(in)                                        :: NAO
            integer, intent(in)                                        :: NCholesky
            integer(I64), intent(in)                                   :: MaxBufferDim

            integer :: p, q, pq, pq0, pq1
            integer :: NA, NB
            integer :: MaxNk, NPasses, Nk
            integer :: k, a0k, a1k
            integer(I64) :: NW
            real(F64), dimension(:, :), allocatable :: CAT
            real(F64), dimension(:), allocatable :: W
            type(TClock) :: timer, timer_total

            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            allocate(CAT(NA, NAO))
            CAT = transpose(CA(:, a0:a1))
            MaxNk = min(NA, int(MaxBufferDim/(NCholesky*NAO)))
            if (.not. MaxNk > 0) then
                  call msg("Buffer size too small for Rpq->Rab")
                  error stop
            end if
            NPasses = NA / MaxNk
            if (modulo(NA, MaxNk) > 0) then
                  NPasses = NPasses + 1
            end if
            allocate(W(NCholesky*MaxNk*NAO))

            call msg("Two step transformation of Cholesky vectors: R(k,pq)->R(k,ab)")
            call msg("Will perform " // str(NPasses) // " passes over R(k,pq)")
            call clock_start(timer_total)
            Rab = ZERO
            do k = 1, NPasses
                  call clock_start(timer)
                  ! ------------------------------------------------
                  ! First index transformation Rpq -> Raq
                  ! ------------------------------------------------
                  a0k = (k - 1) * MaxNk + a0
                  a1k = min(k * MaxNk, a1)
                  Nk = a1k - a0k + 1
                  call chol_MOTransf_Step1(W, Nk, NAO, CAT, a0k, a1k, a0, a1, Rpq, NCholesky)
                  ! ------------------------------------------------
                  ! Second index transformation Raq -> Rab
                  ! ------------------------------------------------
                  call chol_MOTransf_Step2(Rab, W, CB, a0k, a1k, Nk, a0, a1, b0, b1, NCholesky, NAO)
                  call msg("Completed pass " // str(k) // " in " // str(clock_readwall(timer),d=1) // " seconds")
            end do
            call msg("Completed transformation in " // str(clock_readwall(timer_total),d=1) // " seconds")
      end subroutine chol_MOTransf_TwoStep_x


      subroutine chol_MOTransf_Step1(W, Nk, NAO, CAT, a0k, a1k, a0, a1, Rpq, NCholesky)
            real(F64), dimension(NCholesky, Nk, NAO), intent(out) :: W
            integer, intent(in)                                   :: Nk
            integer, intent(in)                                   :: NAO
            integer, intent(in)                                   :: a0, a1
            real(F64), dimension(a0:a1, 1:NAO), intent(in)        :: CAT
            integer, intent(in)                                   :: a0k, a1k
            real(F64), dimension(:, :), intent(in)                :: Rpq
            integer, intent(in)                                   :: NCholesky

            integer :: pq
            integer :: p, q
            
            W = ZERO
            pq = 0
            do q = 1, NAO
                  do p = 1, q
                        pq = pq + 1
                        !
                        ! W(:, :, p) <- W(:, :, p) + R(1:NVecs, pq) * CA(q, b0k:b1k)**T
                        ! W(:, :, q) <- W(:, :, q) + R(1:NVecs, pq) * CA(p, b0k:b1k)**T
                        !
                        if (p /= q) then
                              call chol_vwT_x(W(:, :, p), NCholesky, Rpq(:, pq), CAT(a0k:a1k, q), NCholesky, Nk, ONE)
                              call chol_vwT_x(W(:, :, q), NCholesky, Rpq(:, pq), CAT(a0k:a1k, p), NCholesky, Nk, ONE)
                        else
                              call chol_vwT_x(W(:, :, p), NCholesky, Rpq(:, pq), CAT(a0k:a1k, q), NCholesky, Nk, ONE)
                        end if
                  end do
            end do
      end subroutine chol_MOTransf_Step1


      subroutine chol_MOTransf_Step2(Rab, W, CB, a0k, a1k, Nk, a0, a1, b0, b1, NCholesky, NAO)
            integer, intent(in)                                           :: a0, a1
            integer, intent(in)                                           :: b0, b1
            real(F64), dimension(NCholesky, a0:a1, b0:b1), intent(inout)  :: Rab
            real(F64), dimension(NCholesky, Nk, NAO), intent(in)          :: W
            real(F64), dimension(:, :), intent(in)                        :: CB
            integer, intent(in)                                           :: a0k, a1k
            integer, intent(in)                                           :: Nk
            integer, intent(in)                                           :: NCholesky
            integer, intent(in)                                           :: NAO

            integer :: m, n, ldW, b
            !
            ! Transform the second AO index by a series of vector-matrix multiplications
            ! --------------------------------------------------------------------------
            !
            ! Explicit loop:
            ! do b = b0, b1
            !       do p = 1, NAO
            !             Rab(:, a0k:a1k, b) = Rab(:, a0k:a1k, b) + W(:, :, p) * CB(p, b)
            !       end do
            ! end do
            !
            ! Call a linear algebra subroutine for each of individual
            ! vector-matrix multiplications
            !
            ! Rab(:, a0k:a1k, b) <- Rab(:, a0k:a1k, b) + Sum(p=1...NAO) W(:, :, p) CB(p, b)
            !
            ! Note that the above cannot be done as a single matrix multiplication due to
            ! the non-contiguous memory layout of Rab.
            !
            m = NCholesky * Nk
            n = NAO
            ldW = NCholesky * Nk
            do b = b0, b1
                  call chol_av_x(Rab(:, a0k:a1k, b), W, ldW, CB(:, b), m, n, ONE, ONE)
            end do
      end subroutine chol_MOTransf_Step2
      

      subroutine chol_Significant_Sorted(Significant, SignificantDim, NSignificant, D, TraceError, &
            NOrbPairs, Vdiag, TargetTraceError, MaxBatchDim)
            
            integer, dimension(:, :), allocatable, intent(out) :: Significant
            integer, dimension(:), allocatable, intent(out)    :: SignificantDim
            integer, intent(out)                               :: NSignificant
            real(F64), dimension(:), intent(out)               :: D
            real(F64), intent(out)                             :: TraceError
            integer, dimension(:), intent(inout)               :: NOrbPairs
            real(F64), dimension(:), intent(in)                :: Vdiag
            real(F64), intent(in)                              :: TargetTraceError
            integer, intent(in)                                :: MaxBatchDim

            real(F64) :: DiscardedTrace, NextContrib
            integer, dimension(:), allocatable :: IdxMap
            integer :: NDiscarded
            integer :: k, l, kl

            call msg("Prescreening diagonal integrals (pq|pq) for the Cholesky decomposition")
            call msg("Sorted orbital pairs will be discarded as long as Tr(V-Vapprox) < " // str(TargetTraceError,d=1))
            allocate(IdxMap(NOrbPairs(BASE_PAIRS)))
            D = Vdiag
            do k = 1, NOrbPairs(BASE_PAIRS)
                  IdxMap(k) = k
            end do
            !
            ! Sort the diagonal matrix elements in decreasing order
            !
            call dsort0(D, IdxMap, NOrbPairs(BASE_PAIRS), -2)
            !
            ! Determine the trace of discarded elements
            !
            DiscardedTrace = ZERO
            NDiscarded = 0
            do k = NOrbPairs(BASE_PAIRS), 1, -1
                  NextContrib = max(ZERO, D(k))
                  if (DiscardedTrace+NextContrib<TargetTraceError) then
                        DiscardedTrace = DiscardedTrace + NextContrib
                        NDiscarded = NDiscarded + 1
                  else
                        exit
                  end if
            end do
            TraceError = DiscardedTrace
            !
            ! Define the set of significant orbital pairs
            !
            NOrbPairs(SIGNIFICANT_PAIRS) = NOrbPairs(BASE_PAIRS) - NDiscarded
            if (NOrbPairs(SIGNIFICANT_PAIRS) > 0) then
                  NSignificant = ceiling(real(NOrbPairs(SIGNIFICANT_PAIRS), F64) / real(MaxBatchDim, F64))
                  allocate(Significant(MaxBatchDim, NSignificant))
                  allocate(SignificantDim(NSignificant))
                  Significant = -1
                  SignificantDim = 0
                  do kl = 1, NOrbPairs(SIGNIFICANT_PAIRS)
                        !
                        ! kl = MaxBatchDim*(k-1)+l
                        !
                        k = (kl - 1) / MaxBatchDim + 1
                        l = kl - MaxBatchDim * (k - 1)
                        Significant(l, k) = IdxMap(kl)
                        SignificantDim(k) = l
                  end do
                  D = ZERO
                  do k = 1, NSignificant
                        do l = 1, SignificantDim(k)
                              D(Significant(l, k)) = Vdiag(Significant(l, k))
                        end do
                  end do
            else
                  call msg("Prescreening failed: Found zero significant Coulomb integrals")
                  error stop
            end if
            call msg("Trace of exact Coulomb matrix: " // str(sum(Vdiag), d=1))
            call msg("Trace of discarded integrals: " // str(TraceError, d=1))
            call msg("Discarded " // str(NDiscarded) // " out of " // str(NOrbPairs(BASE_PAIRS)) // " orbital pairs")
      end subroutine chol_Significant_Sorted


      subroutine chol_Qualified(Qualified, QualifiedDim, NQualified, BufferLoc, NOrbPairs, &
            Significant, SignificantDim, NSignificant, D, Dmin, MaxNQualified, MaxBatchDim)
            
            integer, dimension(:, :), intent(out) :: Qualified
            integer, dimension(:), intent(out)    :: QualifiedDim
            integer, intent(out)                  :: NQualified
            integer, dimension(:, :), intent(out) :: BufferLoc
            integer, dimension(:), intent(inout)  :: NOrbPairs
            integer, dimension(:, :), intent(in)  :: Significant
            integer, dimension(:), intent(in)     :: SignificantDim
            integer, intent(in)                   :: NSignificant
            real(F64), dimension(:), intent(in)   :: D
            real(F64), intent(in)                 :: Dmin
            integer, intent(in)                   :: MaxNQualified
            integer, intent(in)                   :: MaxBatchDim

            integer :: k, l, N
            integer :: BatchIdx
            integer :: NColumns
            real(F64) :: BatchVal
            real(F64), dimension(:), allocatable :: DiagSort, BatchVals
            integer, dimension(:), allocatable :: DiagISort
            
            allocate(DiagSort(NSignificant))
            allocate(DiagISort(NSignificant))
            allocate(BatchVals(MaxBatchDim))
            do k = 1, NSignificant
                  N = SignificantDim(k)
                  if (N > 0) then
                        do l = 1, SignificantDim(k)
                              BatchVals(l) = D(Significant(l, k))
                        end do
                        DiagSort(k) = maxval(BatchVals(1:N))
                  else
                        DiagSort(k) = ZERO
                  end if
                  DiagISort(k) = k
            end do
            !
            ! Sort the diagonal shells in descending order. Take into account only
            ! the largest integral in a shell quartet.
            !
            call dsort0(DiagSort, DiagISort, NSignificant, -2)
            NQualified = 0
            NColumns = 0
            do k = 1, NSignificant
                  BatchVal = DiagSort(k)
                  BatchIdx = DiagISort(k)
                  N = SignificantDim(BatchIdx)
                  if (BatchVal > Dmin .and. NQualified+1 <= MaxNQualified) then
                        Qualified(1:N, k) = Significant(1:N, BatchIdx)
                        QualifiedDim(k) = N
                        NQualified = NQualified + 1
                        do l = NColumns+1, NColumns+N
                              BufferLoc(l-NColumns, k) = l
                        end do
                        NColumns = NColumns + N
                  else
                        exit
                  end if
            end do
            NOrbPairs(QUALIFIED_PAIRS) = sum(QualifiedDim(1:NQualified))
      end subroutine chol_Qualified
      

      subroutine chol_SubtractR(M, R, RQ, NCholesky, NOrbPairs)
            !
            ! Subtract the contribution from all available Cholesky vectors from the Coulomb integrals matrix.
            ! (matrix DeltaPQ, step 5f in Ref. 1)
            !
            ! 1. Aquilante et al. Subsection 13.5, Cholesky Decomposition Techniques in Electronic Structure Theory in
            !    Linear-Scaling Techniques in Computational Chemistry and Physics: Methods and Applications,
            !    301-343, Springer 2011; doi: 10.1007/978-90-481-2853-2_13
            !
            real(F64), dimension(:, :), intent(inout)      :: M
            real(F64), dimension(:, :), intent(in)         :: R
            real(F64), dimension(:, :), intent(in)         :: RQ
            integer, intent(in)                            :: NCholesky
            integer, dimension(:), intent(in)              :: NOrbPairs
            
            integer :: ldR, ldRQ, ldM

            ldR = size(R, dim=1)
            ldRQ = size(RQ, dim=1)
            ldM = size(M, dim=1)
            !
            ! M <- M - R(1:NCholesky, 1:NOrbPairs(BASE_PAIRS))**T R(1:NCholesky, NOrbPairs(QUALIFIED_PAIRS))
            !
            call chol_aTb_x(M, ldM, R, ldR, RQ, ldRQ, &
                  NOrbPairs(BASE_PAIRS), NOrbPairs(QUALIFIED_PAIRS), NCholesky, &
                  alpha=-ONE, beta=ONE)
      end subroutine chol_SubtractR

      
      subroutine chol_Qmax(QmaxLocS, QmaxLocQ, Qmax, D, Qualified, QualifiedDim, NQualified, BufferLoc)
            !
            ! Compute the largest diagonal element and its storage location in the QUALIFIED_PAIRS set.
            !
            integer, intent(out)                 :: QmaxLocS
            integer, intent(out)                 :: QmaxLocQ
            real(F64), intent(out)               :: Qmax
            real(F64), dimension(:), intent(in)  :: D
            integer, dimension(:, :), intent(in) :: Qualified
            integer, dimension(:), intent(in)    :: QualifiedDim
            integer, intent(in)                  :: NQualified
            integer, dimension(:, :), intent(in) :: BufferLoc

            integer :: k, l
            real(F64) ::  Qab

            Qmax = -ONE
            do k = 1, NQualified
                  do l = 1, QualifiedDim(k)
                        Qab = D(Qualified(l, k))
                        if (Qab > Qmax) then
                              Qmax = Qab
                              QmaxLocS = Qualified(l, k)
                              QmaxLocQ = BufferLoc(l, k)
                        end if
                  end do
            end do
      end subroutine chol_Qmax


      subroutine chol_Update_D(D, R, n)
            !
            ! Update the vector of residual diagonal elements
            ! Step 5.h.iv of Ref. 1
            !
            real(F64), dimension(:), intent(inout) :: D
            real(F64), dimension(:), intent(in)    :: R
            integer, intent(in)                    :: n

            integer :: p

            do p = 1, n
                  D(p) = D(p) - R(p)**2
            end do
      end subroutine chol_Update_D


      subroutine chol_KetStorageLoc(KetStorageLoc, BufferLoc, Qualified, QualifiedDim, NQualified)
            integer, dimension(:), intent(out)   :: KetStorageLoc
            integer, dimension(:, :), intent(in) :: BufferLoc
            integer, dimension(:, :), intent(in) :: Qualified
            integer, dimension(:), intent(in)    :: QualifiedDim
            integer, intent(in)                  :: NQualified
            
            integer :: k, l, rs
            
            KetStorageLoc = 0
            do k = 1, NQualified
                  do l = 1, QualifiedDim(k)
                        rs = Qualified(l, k)
                        KetStorageLoc(rs) = BufferLoc(l, k)
                  end do
            end do
      end subroutine chol_KetStorageLoc

      
      subroutine chol_M(M, KetStorageLoc, AOInts)
            real(F64), dimension(:, :), intent(out) :: M
            integer, dimension(:), intent(in)       :: KetStorageLoc
            type(AOReaderChol)                      :: AOInts

            call AOInts%get_TR(KetStorageLoc, M)
      end subroutine chol_M


      subroutine chol_Rk(Rk, Tk, Mk, RQk, R, j0, j1, NOrbPairs, Qmax)
            real(F64), dimension(:), intent(out)   :: Rk
            real(F64), dimension(:), intent(out)   :: Tk            
            real(F64), dimension(:), intent(in)    :: Mk
            real(F64), dimension(:), intent(in)    :: RQk
            real(F64), dimension(:, :), intent(in) :: R
            integer, intent(in)                    :: j0, j1
            integer, dimension(:), intent(in)      :: NOrbPairs
            real(F64), intent(in)                  :: Qmax

            integer :: p, j
            real(F64) :: t

            if (j1 >= j0) then
                  !$omp parallel do private(p, t, j) default(shared)
                  do p = 1, NOrbPairs(BASE_PAIRS)
                        t = ZERO
                        do j = j0, j1
                              t = t + R(j, p) * RQk(j-j0+1)
                        end do
                        Tk(p) = Mk(p) - t
                  end do
                  !$omp end parallel do
            else
                  Tk = Mk
            end if
            Rk = Tk / sqrt(Qmax)
      end subroutine chol_Rk
      

      subroutine chol_MainLoop(R, NCholesky, D, Significant, SignificantDim, NSignificant, PrescreenError, &
            MaxNCholesky, TargetTraceError, TargetTraceErrorPrescreen, TargetMaxError, AOInts, NOrbPairs, &
            MaxNQualified, MaxBatchDim, NAO)
            !
            ! Generate the Cholesky vectors of the Coulomb matrix decomposition V = R**T * R.
            !
            ! 1. Harbrecht, H., Peters, M., and Schneider, R. Appl. Num. Math. 62, 428 (2012);
            !    doi: 10.1016/j.apnum.2011.10.001
            ! 2. Koch, H., de Meras, A.S., and Pedersen, T.B. J. Chem. Phys. 118, 9481 (2003);
            !    doi: 10.1063/1.1578621
            !
            real(F64), dimension(:, :), allocatable, intent(out)      :: R
            integer, intent(out)                                      :: NCholesky
            real(F64), dimension(:), intent(inout)                    :: D
            integer, dimension(:, :), intent(in)                      :: Significant
            integer, dimension(:), intent(in)                         :: SignificantDim
            integer, intent(in)                                       :: NSignificant
            real(F64), intent(in)                                     :: PrescreenError
            integer, intent(in)                                       :: MaxNCholesky
            real(F64), intent(in)                                     :: TargetTraceError
            real(F64), intent(in)                                     :: TargetTraceErrorPrescreen
            real(F64), intent(in)                                     :: TargetMaxError            
            type(AOReaderChol)                                        :: AOInts
            integer, dimension(:), intent(inout)                      :: NOrbPairs
            integer, intent(in)                                       :: MaxNQualified
            integer, intent(in)                                       :: MaxBatchDim
            integer, intent(in)                                       :: NAO

            integer, dimension(:, :), allocatable :: BufferLoc
            real(F64), dimension(:), allocatable :: Rk, Tk, RQk
            real(F64), dimension(:, :), allocatable :: RQ
            real(F64), dimension(:, :), allocatable :: M
            integer, dimension(:, :), allocatable :: Qualified
            integer, dimension(:), allocatable :: QualifiedDim
            integer, dimension(:), allocatable :: KetStorageLoc
            integer :: NQualified
            real(F64) :: Dmax, Dmin
            integer :: i, j0, j1, j
            integer :: k, l, p, pp
            integer :: QmaxLocQ, QmaxLocS
            real(F64) :: Qmax
            real(F64) :: SumD, DiscardError, TraceError
            logical :: Converged
            logical :: DoMicroIters, DoMacroIters
            character(:), allocatable :: line
            real(F64) :: time_rTr, time_2e, time_rk
            type(tclock) :: timer, timer_iter, DeltaT
            real(F64) :: memory_M, memory_R

            call clock_start(timer)
            time_rTr = ZERO
            time_rk = ZERO
            time_2e = ZERO
            call blankline()
            call msg("Starting pivoted Cholesky decomposition")
            call msg("Dimensions of the Coulomb matrix: (" &
                     // str(NOrbPairs(BASE_PAIRS)) // "," // str(NOrbPairs(BASE_PAIRS)) // ")")
            call msg("Min decomposed Dpq: " // str(TargetMaxError, d=1))
            call msg("Target trace error: " // str(TargetTraceError, d=1))
            call msg("Target trace error for prescreening: " // str(TargetTraceErrorPrescreen, d=1))
            call msg("Maximum number of Cholesky vectors: " // str(MaxNCholesky))
            call msg("Dimensions of 2e integral storage buffer: (" &
                     // str(NOrbPairs(BASE_PAIRS)) // "," // str(MaxNQualified*MaxBatchDim) // ")")
            allocate(R(MaxNCholesky, NOrbPairs(BASE_PAIRS)))
            allocate(Qualified(MaxBatchDim, MaxNQualified))
            allocate(QualifiedDim(MaxNQualified))
            allocate(BufferLoc(MaxBatchDim, MaxNQualified))
            allocate(RQ(MaxNCholesky, MaxNQualified*MaxBatchDim))
            allocate(M(NOrbPairs(BASE_PAIRS), MaxNQualified*MaxBatchDim))
            allocate(Rk(NOrbPairs(BASE_PAIRS)))
            allocate(Tk(NOrbPairs(BASE_PAIRS)))
            allocate(KetStorageLoc(NOrbPairs(BASE_PAIRS)))
            allocate(RQk(MaxNQualified*MaxBatchDim))
            memory_R = io_size_byte(R) / real(1024**3, F64)
            memory_M = io_size_byte(M) / real(1024**3, F64)
            call msg("Memory allocation per image")
            call msg("Cholesky vectors (R): " // str(memory_R, d=1) // " gigabytes")
            call msg("two-electron integrals (M):  " // str(memory_M, d=1) // " gigabytes")
            DiscardError = PrescreenError
            R = ZERO
            NCholesky = 0
            i = 0
            SumD = sum(D)
            TraceError = PrescreenError + SumD
            Dmax = maxval(D)
            line = lfield("#", 5) // lfield("Tr(V-R**T*R)", 14) // lfield("Dmin", 12) // lfield("Dmax", 12) &
                  // lfield("NCholesky", 12) // lfield("Time", 12)
            call midrule(width=63)
            call msg(line)
            call midrule(width=63)
            Converged = .false.
            DoMacroIters = (MaxNCholesky > 0)
            MacroIters: do while (DoMacroIters)
                  call clock_start(timer_iter)
                  i = i + 1
                  !
                  ! Estimate the decomposition threshold corresponding to the trace error target.
                  ! This value will get an update in each iteration of the macro loop. When the convergence
                  ! condition is controlled by TargetTraceError, the spread factor varies depending
                  ! on Dmax and is increased near the convergence to minimize the number of Cholesky vectors.
                  ! 
                  ! Small spread factor (0.01): less iterations, more vectors
                  ! Large spread factor (0.1): more iterations, less vectors
                  !
                  if (Dmax < TargetTraceError) then
                        Dmin = max(0.1_F64 * Dmax, TargetMaxError)
                  else
                        Dmin = max(0.01_F64 * Dmax, TargetMaxError)
                  end if
                  call chol_Qualified(Qualified, QualifiedDim, NQualified, BufferLoc, NOrbPairs, Significant, &
                        SignificantDim, NSignificant, D, Dmin, MaxNQualified, MaxBatchDim)
                  call clock_start(DeltaT)                  
                  call chol_KetStorageLoc(KetStorageLoc, BufferLoc, Qualified, QualifiedDim, NQualified)
                  call chol_M(M, KetStorageLoc, AOInts)
                  time_2e = time_2e + clock_readwall(DeltaT)
                  if (NCholesky > 0) then
                        call clock_start(DeltaT)
                        !
                        ! Arrange the columns of R corresponding to the QUALIFIED_PAIRS set of indices
                        ! in a contiguous local memory storage RQ. The matrix RQ will then be
                        ! passed to highly optimized linear algebra subroutines.
                        !
                        do k = 1, NQualified
                              do l = 1, QualifiedDim(k)
                                    pp = BufferLoc(l, k)
                                    p = Qualified(l, k)
                                    RQ(:, pp) = R(:, p)
                              end do
                        end  do
                        call chol_SubtractR(M, R, RQ, NCholesky, NOrbPairs)
                        time_rTr = time_rTr + clock_readwall(DeltaT)
                  end if
                  call chol_Qmax(QmaxLocS, QmaxLocQ, Qmax, D, Qualified, QualifiedDim, NQualified, BufferLoc)
                  j0 = NCholesky + 1
                  j = 0
                  DoMicroIters = (NCholesky < MaxNCholesky .and. NOrbPairs(QUALIFIED_PAIRS) > 0)
                  MicroIters: do while (DoMicroIters)
                        j = j + 1
                        j1 = NCholesky
                        NCholesky = NCholesky + 1
                        call clock_start(DeltaT)
                        if (j1 >= j0) then
                              RQk(1:j1-j0+1) = R(j0:j1, QmaxLocS)
                        end if
                        call chol_Rk(Rk, Tk, M(:, QmaxLocQ), RQk, R, j0, j1, NOrbPairs, Qmax)
                        R(NCholesky, :) = Rk
                        call chol_Update_D(D, Rk, NOrbPairs(BASE_PAIRS))
                        time_rk = time_rk + clock_readwall(DeltaT)
                        call chol_Qmax(QmaxLocS, QmaxLocQ, Qmax, D, Qualified, QualifiedDim, NQualified, BufferLoc)
                        SumD = sum(D)
                        TraceError = PrescreenError + SumD
                        if (TraceError < TargetTraceError) then
                              Converged = .true.
                              DoMicroIters = .false.
                        end if
                        if (Qmax <= Dmin) then
                              DoMicroIters = .false.
                        end if
                        if (j == NOrbPairs(QUALIFIED_PAIRS) .or. NCholesky == MaxNCholesky) then
                              DoMicroIters = .false.
                        end if
                  end do MicroIters
                  line = lfield(str(i), 5) // lfield(str(TraceError,d=2), 14) // lfield(str(Dmin,d=1),12) &
                        // lfield(str(Dmax,d=1),12) // lfield(str(NCholesky), 12) &
                        // lfield(str(clock_readwall(timer_iter),d=1), 12)
                  call msg(line)
                  Dmax = maxval(D)
                  if (TraceError < TargetTraceError .or. Dmax <= TargetMaxError) then
                        Converged = .true.
                        DoMacroIters = .false.
                  end if
                  if (NCholesky == MaxNCholesky .or. NOrbPairs(QUALIFIED_PAIRS) == 0) then
                        DoMacroIters = .false.
                  end if
            end do MacroIters
            call blankline()
            call msg("Cholesky decomposition completed")
            call msg("Computed " // str(NCholesky) // " Cholesky vectors")
            call msg("Trace of discarded matrix elements: " // str(DiscardError, d=1))
            call msg("Largest residual Dpq: " // str(Dmax, d=1))
            TraceError = PrescreenError + sum(D)
            call msg("Total trace error: Tr(Vexact-R**T R) = " // str(TraceError, d=1))
            call msg("NCholesky/NAO = " // str(real(NCholesky,F64)/real(NAO,F64), d=1))
            if (NCholesky == MaxNCholesky .and. .not. Converged) then
                  call msg("Max number of Cholesky vectors reached without convergence")
                  error stop
            end if
            call msg("Cholesky decomposition completed in " // str(clock_readwall(timer),d=1) // " seconds")
            call msg("Detailed timings in seconds")
            call msg("subtract prev contribs I (M-R**T*R): " // str(time_rTr,d=1))
            call msg("subtract prev contribs II (Rk): " // str(time_rk,d=1))
            call msg("two-electron integrals: " // str(time_2e,d=1))
            call blankline()
      end subroutine chol_MainLoop
      

      subroutine chol_CoulombMatrix(Vecs, NAO, IntegralsPath, AOSource, Accuracy)
            !
            ! Compute the matrix of Cholesky vectors, R, defined as
            !
            ! V = R**T * R
            !
            ! or
            !
            ! V(pq, rs) = Sum(k=1,...,NCholesky) R(k, pq) R(k, rs)
            !
            ! where V is the matrix of two electron Coulomb integrals (pq|rs).
            !
            ! Small diagonal Coulomb integrals are removed from the computations
            ! at the prescreening phase. On exit, the leading dimension of R
            ! doesn't equal NCholesky. (R is allocated to the max estimated number
            ! of Cholesky vectors.) Use size(R, dim=1) to get the leading dimension
            ! of R.
            !
            ! The algorithm is based on Ref. 1 with added sparse matrix capability
            ! and the target decomposition error defined as Tr(V-R**T*R).
            ! This convergence condition is used by Harbrecht et al. in Ref. 2.
            !
            ! Accuracy settings, ordered from the lowest accuracy to near-exact:
            ! CHOL_ACCURACY_DEFAULT
            ! CHOL_ACCURACY_TIGHT
            ! CHOL_ACCURACY_LUDICROUS
            ! The constants corresponding to pre-defined accuracy settings are
            ! defined as global parameters in this module.
            !
            ! 1. Aquilante et al. Subsection 13.5, Cholesky Decomposition Techniques in Electronic Structure Theory in
            !    Linear-Scaling Techniques in Computational Chemistry and Physics: Methods and Applications,
            !    301-343, Springer 2011; doi: 10.1007/978-90-481-2853-2_13
            ! 2. Harbrecht, H., Peters, M., and Schneider R. Appl. Num. Math. 62, 428 (2012); doi: 10.1016/j.apnum.2011.10.001
            !
            type(TCholeskyVecs), intent(out)  :: Vecs
            integer, intent(in)               :: NAO
            character(*), intent(in)          :: IntegralsPath
            integer, intent(in)               :: AOSource
            integer, intent(in)               :: Accuracy
            
            real(F64) :: PrescreenError
            real(F64), dimension(:), allocatable :: Vdiag
            real(F64), dimension(:), allocatable :: D
            integer, dimension(:, :), allocatable :: Significant
            integer, dimension(:), allocatable :: SignificantDim
            integer :: NSignificant
            type(AOReaderChol) :: AOInts
            integer, parameter :: MaxBatchDim = 10
            integer, parameter :: MaxNQualified = 100
            real(F64) :: TargetTraceError, TargetTraceErrorPrescreen, TargetMaxError
            integer :: MaxNAOMult

            select case (Accuracy)
            case (CHOL_ACCURACY_DEFAULT)
                  TargetTraceError = 1.0E-2_F64
                  TargetTraceErrorPrescreen = 1.0E-10_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 8
            case (CHOL_ACCURACY_TIGHT)
                  TargetTraceError = 1.0E-3_F64
                  TargetTraceErrorPrescreen = 1.0E-10_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 9
            case (CHOL_ACCURACY_LUDICROUS)
                  TargetTraceError = 1.0E-4_F64
                  TargetTraceErrorPrescreen = 1.0E-10_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 10
            end select
            associate ( &
                  NOrbPairs => Vecs%NOrbPairs, &
                  MaxNCholesky => Vecs%MaxNCholesky, &
                  NCholesky => Vecs%NCholesky)
                  Vecs%NAO = NAO
                  !
                  ! Read diagonal integrals (pq|pq) for all unique pq
                  !
                  allocate(Vdiag(NAO*(NAO+1)/2))
                  call AOInts%open(NAO, AOSource, IntegralsPath, Vdiag)
                  
                  NOrbPairs(BASE_PAIRS) = AOInts%nelm
                  MaxNCholesky = min(MaxNAOMult * NAO, NOrbPairs(BASE_PAIRS))
                  allocate(D(NOrbPairs(BASE_PAIRS)))
                  !
                  ! Prescreen small matrix elements on the diagonal
                  !
                  call chol_Significant_Sorted(Significant, SignificantDim, NSignificant, D, PrescreenError, &
                        NOrbPairs, Vdiag, TargetTraceErrorPrescreen, MaxBatchDim)
                  call chol_MainLoop(Vecs%R, NCholesky, D, Significant, SignificantDim, NSignificant, PrescreenError, &
                        MaxNCholesky, TargetTraceError, TargetTraceErrorPrescreen, TargetMaxError, AOInts, &
                        NOrbPairs, MaxNQualified, MaxBatchDim, NAO)
                  call AOInts%close()
            end associate
      end subroutine chol_CoulombMatrix
end module Cholesky
