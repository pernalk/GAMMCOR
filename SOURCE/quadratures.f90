module quadratures
      use arithmetic
      use math_constants
      use io
      use sort

      implicit none

!      double precision, parameter, private :: PI = 3.141592653589793238462643383279d+0

    contains
      function BesselZero(s)
            !
            ! Compute a rational approximation for the s-th positive zero
            ! of the Bessel function J0(x). Use use the most accurate
            ! formula in Table I of Ref. 1.
            !
            ! 1. Lether, F.G. J. Comp. Appl. Math. 67, 167 (1996);
            !    doi: 10.1016/0377-0427(95)00219-7
            !
            real(F64) :: BesselZero
            integer, intent(in) :: s

            real(F64) :: beta, t, R
            real(F64), parameter :: p1 = 1567450796.0_F64/12539606369.0_F64
            real(F64), parameter :: p2 = 8903660.0_F64/2365861.0_F64
            real(F64), parameter :: p3 = 10747040.0_F64/536751.0_F64
            real(F64), parameter :: p4 = 17590991.0_F64/1696654.0_F64
            real(F64), parameter :: q1 = 1.0_F64
            real(F64), parameter :: q2 = 29354255.0_F64/954518.0_F64
            real(F64), parameter :: q3 = 76900001.0_F64/431847.0_F64
            real(F64), parameter :: q4 = 67237052.0_F64/442411.0_F64

            beta = (s - 1.D0/4.D0) * PI
            t = (1.D0/beta)**2
            R = (p1 + p2*t + p3*t**2 + p4*t**3) / (q1 + q2*t + q3*t**2 + q4*t**3)
            BesselZero = beta + R / beta
          end function BesselZero

            function GLGuess_Tricomi(k, n)
            !
            ! Compute Tricomi's guess for the node arccos(x(k)), x(k)>=0 of the Gauss-Legendre quadrature.
            ! Use Eq. 3.4 in Ref. 1
            !
            ! 1. Hale, N. and Townsend, A. SIAM J. Sci. Comput., 35, A652 (2013);
            ! doi: 10.1137/120889873
            !
            real(F64) :: GLGuess_Tricomi
            integer, intent(in) :: k
            integer, intent(in) :: n

            integer :: kk
            real(F64) :: x, phi

            kk = n + 1 - k
            phi = (kk - 1.D0/4.D0)*PI / (n + 1.D0/2.D0)
            x = (1.D0 - (n - 1.D0) / (8.D0 * n**3) - &
                  (39.0_F64 - 28.0_F64/sin(phi)**2) / (384.0_F64 * n**4)) * cos(phi)
            GLGuess_Tricomi = acos(x)
      end function GLGuess_Tricomi


        function GLGuess_Gatteschi(k, n)
            !
            ! Compute Gatteschi's guess for the node arccos(x(k)), x(k)>=0
            ! of the Gauss-Legendre quadrature. Use Eq. 3.7 in Ref. 1.
            !
            ! 1. Hale, N. and Townsend, A. SIAM J. Sci. Comput., 35, A652 (2013);
            !    doi: 10.1137/120889873
            !
            real(F64) :: GLGuess_Gatteschi
            integer, intent(in) :: k
            integer, intent(in) :: n

            integer :: kk
            real(F64) :: psi

            kk = n + 1 - k
            psi = BesselZero(kk) / (n + 1.D0/2.D0)
            GLGuess_Gatteschi = psi + (psi * cos(psi)/sin(psi) - 1.D0) / (8.D0 * psi * (n + 1.D0/2.D0)**2)
      end function GLGuess_Gatteschi

      subroutine GLGuess(Theta, n)
            !
            ! Compute guess roots of the Gauss-Legendre quadrature.
            ! Use Eq. 3.8 in Ref. 1.
            !
            ! 1. Hale, N. and Townsend, A. SIAM J. Sci. Comput., 35, A652 (2013);
            !    doi: 10.1137/120889873
            !
            real(F64), dimension(:), intent(out) :: Theta
            integer, intent(in) :: n

            integer :: k, l, k0
            real(F64) :: ThetaK

            k0 = n/2 + 1
            k = k0
            if (modulo(n, 2) == 1) then
                  ThetaK = PI/2.D0
            else
                  ThetaK = GLGuess_Tricomi(k, n)
            end if
            Theta(1) = ThetaK
            do while (ThetaK >= PI/3.D0)
                  k = k + 1
                  ThetaK = GLGuess_Tricomi(k, n)
                  Theta(k-k0+1) = ThetaK                  
            end do
            do l = k, n
                  Theta(l-k0+1) = GLGuess_Gatteschi(l, n)
            end do
      end subroutine GLGuess    

      subroutine LegendrePolynomial(Pa, Pb, Pc, x, n)
            !
            ! Compute Legendre polynomials of order n for a vector of abscissas x.
            !
            real(F64), dimension(:), intent(out) :: Pa
            real(F64), dimension(:), intent(out) :: Pb
            real(F64), dimension(:), intent(out) :: Pc
            real(F64), dimension(:), intent(in) :: x
            integer, intent(in) :: n

            integer :: a

            Pb = 0.D0
            Pa = 1.D0
            do a = 1, n
                  Pc = Pb
                  Pb = Pa
                  Pa = (2*a-1.D0)/a * x * Pb - (a-1.D0)/a * Pc
            end do
          end subroutine LegendrePolynomial

          subroutine quad_GaussLegendre(x, w, converged, n, MaxAbsError)
            !
            ! Compute the non-negative roots x and the corresponding
            ! weights w of the n-point Gauss-Legendre quadrature (n >= 2).
            ! The remaining roots are located symmetrically with respect
            ! to the origin.
            !
            ! The arrays x and w should be of length Floor(n/2)+Mod(n,2).
            ! On exit, x contains the non-negative roots of Hn(x), i.e.,
            ! x(Floor(n/2)+1), x(Floor(n/2)+2), ..., x(n).
            ! When n is odd, the first root stored in x is zero.
            !
            ! For 64-bit reals, the recommended threshold is MaxAbsError=1.0E-14;
            ! for that threshold, the Newton root finding successfully
            ! converges for 2 <= n <= 440.
            !
            ! The roots and weights computed with this subroutine
            ! should be used as follows for an even n
            !
            ! int(-1,1) f(x) dx = sum(k=1,n/2) w(k) * (f(x(k)) + f(-x(k)))
            !
            ! For numerical stability, the Newton root finding is carried
            ! out in the theta variable (see Ref. 1). MaxAbsError corresponds
            ! to the errors in the theta variable.
            !
            ! 1. Townsend, A., Trogdon, T., and Olver, S. IMA J. Num. Analysis 36, 337 (2016);
            ! doi: 10.1093/imanum/drv002
            !
            !
            real(F64), dimension(:), intent(out) :: x
            real(F64), dimension(:), intent(out) :: w
            logical, intent(out)                 :: converged
            integer, intent(in)                  :: n
            real(F64), intent(in)                :: MaxAbsError

            real(F64), dimension(:), allocatable :: Pa, Pb, Pc, P_Theta
            real(F64), dimension(:), allocatable :: Theta, ThetaPrev, SinTheta
            integer :: m, l
            integer :: k
            integer, parameter :: maxit = 20

            l = 1 + modulo(n, 2)
            m = n/2 + modulo(n, 2)
            allocate(Pa(m))
            allocate(Pb(m))
            allocate(Pc(m))
            allocate(P_Theta(m))
            allocate(Theta(m))
            allocate(ThetaPrev(m))
            allocate(SinTheta(m))
            call GLGuess(ThetaPrev, n)
            converged = .false.
            if (modulo(n, 2) == 1) then
                  x(1) = 0.D0
                  Theta(1) = PI/2.D0
                  ThetaPrev(1) = PI/2.D0
            end if
            NewtonSteps: do k = 1, maxit
                  x(l:m) = cos(ThetaPrev(l:m))
                  call LegendrePolynomial(Pa, Pb, Pc, x, n)
                  !
                  ! Compute d/dTheta Pn(cos(theta)) = -sin
                  !
                  SinTheta(l:m) = sin(ThetaPrev(l:m))
                  P_Theta(l:m) = n * (x(l:m) * Pa(l:m) - Pb(l:m)) / SinTheta(l:m)
                  !
                  ! Newton step for theta angles
                  !
                  Theta(l:m) = ThetaPrev(l:m) - Pa(l:m) / P_Theta(l:m)
                  if (all(abs(Theta(l:m)-ThetaPrev(l:m)) < MaxAbsError)) then
                        converged = .true.
                        exit NewtonSteps
                  else
                        ThetaPrev = Theta
                  end if
            end do NewtonSteps
            !
            ! Quadrature weights. Formula taken from Numerical Recipes
            ! in Fortran 77
            !
            x(l:m) = cos(Theta(l:m))
            call LegendrePolynomial(Pa, Pb, Pc, x, n)
            SinTheta = sin(Theta)
            P_Theta = n * (x * Pa - Pb) / SinTheta
            !
            ! Weights (Eq. 3.9 in Ref. 1)
            !
            w = 2.D0 / P_Theta**2
      end subroutine quad_GaussLegendre



      subroutine quad_CasimirPolder(x, w, converged, n, MaxAbsError, alpha)
            !
            ! Generate a modified Gauss-Legendre quadrature for Casimir-Polder type
            ! integrals over imaginary frequencies. The conventional variable
            ! of the Gauss-Legendre quadrature is changed to cover the semi-infinite
            ! interval:
            !
            ! xCP(k) = x0 * (1.D0 + xGL(k)) / (1.D0 - xGL(k))
            ! wCP(k) = 2.D0 * wGL(k) * x0 / (1.D0 - xGL(k))**2
            !
            ! This transformation is used in Ref. 1 for the direct RPA correlation
            ! energy and in the SAPT program of Szalewicz et al. [2]. The recommended
            ! value of the scaling parameter is x0=0.5 (see Appendix C of Ref. 1).
            !
            ! The arrays x and w should be of length n. The points and weights
            ! computed with this subroutine should be used as follows (for n>=2):
            !
            ! Int(0, Inf) Pi(u) du = Sum(k=1,n) w(k) * Pi(x(k))
            !
            ! For 64-bit reals, the recommended threshold for the Gauss-Legendre
            ! root finding is MaxAbsError=1.0E-14; the Newton solver successfully
            ! converges for 2 <= n <= 440.
            !
            ! 1. Ren, X., Rinke, P., Blum, V., Wieferink, J., Tkatchenko, A.,
            !    Sanfilippo, A., Reuter, K., and Scheffler, M.,
            !    New J. Phys. 14, 053020 (2012); doi: 10.1088/1367-2630/14/5/053020
            !
            ! 2. SAPT2016: "An Ab Initio Program for Many-Body Symmetry-Adapted
            !    Perturbation Theory Calculations of Intermolecular
            !    Interaction Energies" by R. Bukowski, W. Cencek, P. Jankowski,
            !    B. Jeziorski, M. Jeziorska, T. Korona, S. A. Kucharski, V. F. Lotrich,
            !    M. P. Metz, A. J. Misquitta, R. Moszynski, K. Patkowski, R. Podeszwa,
            !    F. Rob, S. Rybak, K. Szalewicz, H. L. Williams, R. J. Wheatley,
            !    P. E. S. Wormer, and P. S. Żuchowski.
            !
            real(F64), dimension(:), intent(out) :: x
            real(F64), dimension(:), intent(out) :: w
            logical, intent(out)                 :: converged
            integer, intent(in)                  :: n
            real(F64), intent(in)                :: MaxAbsError
            real(F64), optional, intent(in)      :: alpha

            integer :: nGL
            integer :: idx, k
            real(F64) :: t
            real(F64), dimension(:), allocatable :: xGL, wGL
            real(F64), parameter :: x0def = 0.5_F64
            real(F64) :: x0

            nGL = n / 2 + modulo(n, 2)
            allocate(xGL(nGL))
            allocate(wGL(nGL))
            if (present(alpha)) then
                  x0 = alpha
            else
                  x0 = x0def
            end if
            
            call quad_GaussLegendre(xGL, wGL, converged, n, MaxAbsError)

            do k = nGL, 1, -1
                  t = -xGL(k)
                  idx = nGL - k + 1
                  x(idx) = x0 * (1.D0 + t) / (1.D0 - t)
                  w(idx) = 2.D0 * wGL(k) * x0 / (1.D0 - t)**2
            end do

            do k = 1 + modulo(n, 2), nGL
                  t = xGL(k)
                  idx = n / 2 + k
                  x(idx) = x0 * (1.D0 + t) / (1.D0 - t)
                  w(idx) = 2.D0 * wGL(k) * x0 / (1.D0 - t)**2
            end do
      end subroutine quad_CasimirPolder

      subroutine quad_AdiabaticConnection(x, w, n)
            !
            ! The roots and weights computed with this subroutine
            ! should be used as follows for an even n
            !
            ! int(0,1) f(x) dx = sum(k=1,n) w(k) * f(x(k)
            !
            ! For numerical stability, the Newton root finding is carried
            ! out in the theta variable (see Ref. 1). MaxAbsError corresponds
            ! to the errors in the theta variable.
            !
            ! 1. Townsend, A., Trogdon, T., and Olver, S. IMA J. Num. Analysis 36, 337 (2016);
            ! doi: 10.1093/imanum/drv002
            !
            !
            real(F64), dimension(:), intent(out) :: x
            real(F64), dimension(:), intent(out) :: w
            integer, intent(in)                  :: n

            logical :: converged
            integer :: k, k1, k2
            integer :: nGL
            real(F64), dimension(:), allocatable :: x_unsorted, w_unsorted
            integer, dimension(:), allocatable :: permutation
            real(F64), dimension(:), allocatable :: xGL, wGL
            real(F64), parameter :: MaxAbsError = 1.0E-14_F64

            nGL = n / 2 + modulo(n, 2)
            allocate(xGL(nGL))
            allocate(wGL(nGL))
            call quad_GaussLegendre(xGL, wGL, converged, n, MaxAbsError)
            if (.not. converged) then
               print*, "Error while generating the Gauss-Legendre quadrature for adiabatic connection integral"
!                  call msg("Error while generating the Gauss-Legendre quadrature for adiabatic connection integral", MSG_ERROR)
                  error stop
            end if
            allocate(x_unsorted(n))
            allocate(w_unsorted(n))
            if (modulo(n,2) == 1) then
                  w_unsorted(1) = wGL(1) / 2.D0
                  x_unsorted(1) = 1.D0 / 2.D0
                  do k = 2, n/2 + 1
                        k1 = 2*(k - 1)
                        k2 = 2*(k - 1) + 1
                        w_unsorted(k1) = wGL(k) / 2.D0
                        w_unsorted(k2) = wGL(k) / 2.D0
                        x_unsorted(k1) = (1.D0-xGL(k)) / 2.D0
                        x_unsorted(k2) = (1.D0+xGL(k)) / 2.D0
                  end do
            else
                  do k = 1, n/2
                        k1 = 1 + 2*(k - 1)
                        k2 = 2 + 2*(k - 1)
                        w_unsorted(k1) = wGL(k) / 2.D0
                        w_unsorted(k2) = wGL(k) / 2.D0
                        x_unsorted(k1) = (1.D0-xGL(k)) / 2.D0
                        x_unsorted(k2) = (1.D0+xGL(k)) / 2.D0
                  end do
            end if
            allocate(permutation(n))
            do k = 1, n
                  permutation(k) = k
            end do
            call dsort(x_unsorted, permutation, n)
            do k = 1, n
                  x(k) = x_unsorted(k)
                  w(k) = w_unsorted(permutation(k))
            end do
      end subroutine quad_AdiabaticConnection
    end module quadratures
