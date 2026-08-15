module lin

  use real_linalg
  use types
  implicit none
  save

contains


  subroutine Orthogonalize(V, M, StartIdx, EndIdx, S, IdxMap, ACType)
    double precision, dimension(:, :), intent(inout) :: V
    integer, intent(in)                       :: M
    integer, dimension(:), intent(in)         :: StartIdx
    integer, dimension(:), intent(in)         :: EndIdx
    double precision, dimension(:), intent(in)       :: S
    integer, dimension(:), intent(in)                       :: IdxMap
    integer, intent(in)                       :: ACType

    double precision, dimension(:, :), allocatable :: W, SW
    integer :: N, k, l, i
    integer :: MaxD, D
    logical :: PositiveDefinite
    double precision :: t

    if (ACType == 1) then
       PositiveDefinite = .true.
    else
       PositiveDefinite = .false.
    end if
    N = size(V, dim=1)
    MaxD = 1
    do k = 1, M
       D = EndIdx(k) - StartIdx(k) + 1
       MaxD = max(MaxD, D)            
    end do
    allocate(W(N, MaxD))
    allocate(SW(N, MaxD))
    do k = 1, M
       D = EndIdx(k) - StartIdx(k) + 1
       if (D > 1) then
          do l = 1, D
             i = IdxMap(StartIdx(k) + l - 1)
             W(:, l) = V(:, i)
             SW(:, l) = S(:) * V(:, i)
          end do
          call LowdinOrtho(W, SW, D, N, PositiveDefinite)
          do l = 1, D
             i = IdxMap(StartIdx(k) + l - 1)
             V(:, i) = W(:, l)
          end do
       else
          i = IdxMap(StartIdx(k))
          W(:, 1) = V(:, i)
          SW(:, 1) = S(:) * V(:, i)
          t = dot_product(W(:, 1), SW(:, 1))
          V(:, i) = W(:, 1) / Sqrt(Abs(t))
       end if
    end do
  end subroutine Orthogonalize


  subroutine LowdinOrtho(V, SV, D, N, PositiveDefinite)
    integer, intent(in)                       :: D
    integer, intent(in)                       :: N

    double precision, dimension(N, D), intent(inout) :: V
    double precision, dimension(N, D), intent(inout) :: SV
    logical, intent(in)                       :: PositiveDefinite

    double precision, dimension(D, D) :: VSV, Q, LambdaQ
    double precision, dimension(D) :: Lambda
    integer :: k

    call real_aTb(VSV, V, SV)
    Q = VSV
    call symmetric_eigenproblem(Lambda, Q, D, .true.)
    if (.not. PositiveDefinite) Lambda = -Lambda
    do k = 1, D
       LambdaQ(:, k) = Q(:, k) / Sqrt(Lambda(k))
    end do
    call real_abT(VSV, Q, LambdaQ)
    SV = V
    call real_ab(V, SV, VSV)
  end subroutine LowdinOrtho

end module lin
