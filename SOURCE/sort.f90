module sort
      implicit none

contains

      pure subroutine dsort(dx, dy, n)
            !
            ! Sort array of double precision numbers in increasing order.
            ! Carry DY along.
            !
            double precision, dimension(:), intent(inout) :: dx
            integer, dimension(:), intent(inout)          :: dy
            integer, intent(in)                           :: n

            call dsort0(dx, dy, n, 2)
      end subroutine dsort


      PURE SUBROUTINE DSORT0(DX, DY, N, KFLAG)
            !***BEGIN PROLOGUE  DSORT
            !***PURPOSE  Sort an array and optionally make the same interchanges in
            !            an auxiliary array.  The array may be sorted in increasing
            !            or decreasing order.  A slightly modified QUICKSORT
            !            algorithm is used.
            !***LIBRARY   SLATEC
            !***CATEGORY  N6A2B
            !***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
            !***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
            !***AUTHOR  Jones, R. E., (SNLA)
            !           Wisniewski, J. A., (SNLA)
            !***DESCRIPTION
            !
            !   DSORT sorts array DX and optionally makes the same interchanges in
            !   array DY.  The array DX may be sorted in increasing order or
            !   decreasing order.  A slightly modified quicksort algorithm is used.
            !
            !   Description of Parameters
            !      DX - array of values to be sorted   (usually abscissas)
            !      DY - array to be (optionally) carried along
            !      N  - number of values in array DX to be sorted
            !      KFLAG - control parameter
            !            =  2  means sort DX in increasing order and carry DY along.
            !            =  1  means sort DX in increasing order (ignoring DY)
            !            = -1  means sort DX in decreasing order (ignoring DY)
            !            = -2  means sort DX in decreasing order and carry DY along.
            !
            !***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
            !                 for sorting with minimal storage, Communications of
            !                 the ACM, 12, 3 (1969), pp. 185-187.
            !***ROUTINES CALLED  XERMSG
            !***REVISION HISTORY  (YYMMDD)
            !   761101  DATE WRITTEN
            !   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
            !   890531  Changed all specific intrinsics to generic.  (WRB)
            !   890831  Modified array declarations.  (WRB)
            !   891009  Removed unreferenced statement labels.  (WRB)
            !   891024  Changed category.  (WRB)
            !   891024  REVISION DATE from Version 3.2
            !   891214  Prologue converted to Version 4.0 format.  (BAB)
            !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
            !   901012  Declared all variables; changed X,Y to DX,DY; changed
            !           code to parallel SSORT. (M. McClain)
            !   920501  Reformatted the REFERENCES section.  (DWL, WRB)
            !   920519  Clarified error messages.  (DWL)
            !   920801  Declarations section rebuilt and code restructured to use
            !           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
            !***END PROLOGUE  DSORT
            !
            INTEGER, INTENT(IN) :: KFLAG, N
            DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: DX
            INTEGER, DIMENSION(:), INTENT(INOUT) :: DY
            !     .. Local Scalars ..
            DOUBLE PRECISION R, T, TT
            INTEGER :: TY, TTY
            INTEGER I, IJ, J, K, KK, L, M, NN
            !     .. Local Arrays ..
            INTEGER IL(21), IU(21)
            !***FIRST EXECUTABLE STATEMENT  DSORT
            NN = N
            !
            !     Return if nothing is to be done
            !
            IF (N .LE. 1) return

            KK = ABS(KFLAG)
            !
            !     Alter array DX to get decreasing order if needed
            !
            IF (KFLAG .LE. -1) THEN
                  DO I=1,NN
                        DX(I) = -DX(I)
                  END DO
            END IF

            IF (KK .EQ. 2) GO TO 100
            !
            !     Sort DX only
            !
            M = 1
            I = 1
            J = NN
            R = 0.375D0

20          IF (I .EQ. J) GO TO 60
            IF (R .LE. 0.5898437D0) THEN
                  R = R+3.90625D-2
            ELSE
                  R = R-0.21875D0
            END IF

30          K = I
            !
            !     Select a central element of the array and save it in location T
            !
            IJ = I + INT((J-I)*R)
            T = DX(IJ)
            !
            !     If first element of array is greater than T, interchange with T
            !
            IF (DX(I) .GT. T) THEN
                  DX(IJ) = DX(I)
                  DX(I) = T
                  T = DX(IJ)
            END IF
            L = J
            !
            !     If last element of array is less than than T, interchange with T
            !
            IF (DX(J) .LT. T) THEN
                  DX(IJ) = DX(J)
                  DX(J) = T
                  T = DX(IJ)
                  !
                  !        If first element of array is greater than T, interchange with T
                  !
                  IF (DX(I) .GT. T) THEN
                        DX(IJ) = DX(I)
                        DX(I) = T
                        T = DX(IJ)
                  END IF
            END IF
            !
            !     Find an element in the second half of the array which is smaller
            !     than T
            !
40          L = L-1
            IF (DX(L) .GT. T) GO TO 40
            !
            !     Find an element in the first half of the array which is greater
            !     than T
            !
50          K = K+1
            IF (DX(K) .LT. T) GO TO 50
            !
            !     Interchange these elements
            !
            IF (K .LE. L) THEN
                  TT = DX(L)
                  DX(L) = DX(K)
                  DX(K) = TT
                  GO TO 40
            END IF
            !
            !     Save upper and lower subscripts of the array yet to be sorted
            !
            IF (L-I .GT. J-K) THEN
                  IL(M) = I
                  IU(M) = L
                  I = K
                  M = M+1
            ELSE
                  IL(M) = K
                  IU(M) = J
                  J = L
                  M = M+1
            END IF
            GO TO 70
            !
            !     Begin again on another portion of the unsorted array
            !
60          M = M-1
            IF (M .EQ. 0) GO TO 190
            I = IL(M)
            J = IU(M)

70          IF (J-I .GE. 1) GO TO 30
            IF (I .EQ. 1) GO TO 20
            I = I-1

80          I = I+1
            IF (I .EQ. J) GO TO 60
            T = DX(I+1)
            IF (DX(I) .LE. T) GO TO 80
            K = I

90          DX(K+1) = DX(K)
            K = K-1
            IF (T .LT. DX(K)) GO TO 90
            DX(K+1) = T
            GO TO 80
            !
            !     Sort DX and carry DY along
            !
100         M = 1
            I = 1
            J = NN
            R = 0.375D0

110         IF (I .EQ. J) GO TO 150
            IF (R .LE. 0.5898437D0) THEN
                  R = R+3.90625D-2
            ELSE
                  R = R-0.21875D0
            END IF

120         K = I
            !
            !     Select a central element of the array and save it in location T
            !
            IJ = I + INT((J-I)*R)
            T = DX(IJ)
            TY = DY(IJ)
            !
            !     If first element of array is greater than T, interchange with T
            !
            IF (DX(I) .GT. T) THEN
                  DX(IJ) = DX(I)
                  DX(I) = T
                  T = DX(IJ)
                  DY(IJ) = DY(I)
                  DY(I) = TY
                  TY = DY(IJ)
            END IF
            L = J
            !
            !     If last element of array is less than T, interchange with T
            !
            IF (DX(J) .LT. T) THEN
                  DX(IJ) = DX(J)
                  DX(J) = T
                  T = DX(IJ)
                  DY(IJ) = DY(J)
                  DY(J) = TY
                  TY = DY(IJ)
                  !
                  !        If first element of array is greater than T, interchange with T
                  !
                  IF (DX(I) .GT. T) THEN
                        DX(IJ) = DX(I)
                        DX(I) = T
                        T = DX(IJ)
                        DY(IJ) = DY(I)
                        DY(I) = TY
                        TY = DY(IJ)
                  END IF
            END IF
            !
            !     Find an element in the second half of the array which is smaller
            !     than T
            !
130         L = L-1
            IF (DX(L) .GT. T) GO TO 130
            !
            !     Find an element in the first half of the array which is greater
            !     than T
            !
140         K = K+1
            IF (DX(K) .LT. T) GO TO 140
            !
            !     Interchange these elements
            !
            IF (K .LE. L) THEN
                  TT = DX(L)
                  DX(L) = DX(K)
                  DX(K) = TT
                  TTY = DY(L)
                  DY(L) = DY(K)
                  DY(K) = TTY
                  GO TO 130
            END IF
            !
            !     Save upper and lower subscripts of the array yet to be sorted
            !
            IF (L-I .GT. J-K) THEN
                  IL(M) = I
                  IU(M) = L
                  I = K
                  M = M+1
            ELSE
                  IL(M) = K
                  IU(M) = J
                  J = L
                  M = M+1
            END IF
            GO TO 160
            !
            !     Begin again on another portion of the unsorted array
            !
150         M = M-1
            IF (M .EQ. 0) GO TO 190
            I = IL(M)
            J = IU(M)

160         IF (J-I .GE. 1) GO TO 120
            IF (I .EQ. 1) GO TO 110
            I = I-1

170         I = I+1
            IF (I .EQ. J) GO TO 150
            T = DX(I+1)
            TY = DY(I+1)
            IF (DX(I) .LE. T) GO TO 170
            K = I

180         DX(K+1) = DX(K)
            DY(K+1) = DY(K)
            K = K-1
            IF (T .LT. DX(K)) GO TO 180
            DX(K+1) = T
            DY(K+1) = TY
            GO TO 170
            !
            !     Clean up
            !
190         IF (KFLAG .LE. -1) THEN
                  DO I=1,NN
                        DX(I) = -DX(I)
                  END DO
            END IF
      END SUBROUTINE DSORT0


      PURE SUBROUTINE ISORT(IX, IY, N, KFLAG)
            !***  BEGIN PROLOGUE  ISORT
            !***  PURPOSE  Sort an array and optionally make the same interchanges in
            !     an auxiliary array.  The array may be sorted in increasing
            !     or decreasing order.  A slightly modified QUICKSORT
            !     algorithm is used.
            !***  LIBRARY   SLATEC
            !***  CATEGORY  N6A2A
            !***  TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
            !***  KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
            !***  AUTHOR  Jones, R. E., (SNLA)
            !     Kahaner, D. K., (NBS)
            !     Wisniewski, J. A., (SNLA)
            !***  DESCRIPTION
            !     
            !     ISORT sorts array IX and optionally makes the same interchanges in
            !     array IY.  The array IX may be sorted in increasing order or
            !     decreasing order.  A slightly modified quicksort algorithm is used.
            !     
            !     Description of Parameters
            !     IX - integer array of values to be sorted
            !     IY - integer array to be (optionally) carried along
            !     N  - number of values in integer array IX to be sorted
            !     KFLAG - control parameter
            !     =  2  means sort IX in increasing order and carry IY along.
            !     =  1  means sort IX in increasing order (ignoring IY)
            !     = -1  means sort IX in decreasing order (ignoring IY)
            !     = -2  means sort IX in decreasing order and carry IY along.
            !     
            !***  REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
            !     for sorting with minimal storage, Communications of
            !     the ACM, 12, 3 (1969), pp. 185-187.
            !***  ROUTINES CALLED  XERMSG
            !***  REVISION HISTORY  (YYMMDD)
            !     761118  DATE WRITTEN
            !     810801  Modified by David K. Kahaner.
            !     890531  Changed all specific intrinsics to generic.  (WRB)
            !     890831  Modified array declarations.  (WRB)
            !     891009  Removed unreferenced statement labels.  (WRB)
            !     891009  REVISION DATE from Version 3.2
            !     891214  Prologue converted to Version 4.0 format.  (BAB)
            !     900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
            !     901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
            !     920501  Reformatted the REFERENCES section.  (DWL, WRB)
            !     920519  Clarified error messages.  (DWL)
            !     920801  Declarations section rebuilt and code restructured to use
            !     IF-THEN-ELSE-ENDIF.  (RWC, WRB)
            !***  END PROLOGUE  ISORT
            !
            INTEGER, INTENT(IN) :: KFLAG, N
            INTEGER, DIMENSION(:), INTENT(INOUT) :: IX, IY
            !     .. Local Scalars ..
            DOUBLE PRECISION :: R
            INTEGER I, IJ, J, K, KK, L, M, NN, T, TT, TTY, TY
            !     .. Local Arrays ..
            INTEGER IL(21), IU(21)
            !***  FIRST EXECUTABLE STATEMENT  ISORT
            NN = N
            !
            !     Return if nothing is to be done
            !
            IF (N .LE. 1) return

            KK = ABS(KFLAG)
            !     
            !     Alter array IX to get decreasing order if needed
            !     
            IF (KFLAG .LE. -1) THEN
                  DO I=1,NN
                        IX(I) = -IX(I)
                  END DO
            END IF

            IF (KK .EQ. 2) GO TO 100
            !     
            !     Sort IX only
            !     
            M = 1
            I = 1
            J = NN
            R = 0.375D0

20          IF (I .EQ. J) GO TO 60
            IF (R .LE. 0.5898437D0) THEN
                  R = R+3.90625D-2
            ELSE
                  R = R-0.21875D0
            END IF

30          K = I
            !     
            !     Select a central element of the array and save it in location T
            !     
            IJ = I + INT((J-I)*R)
            T = IX(IJ)
            !     
            !     If first element of array is greater than T, interchange with T
            !     
            IF (IX(I) .GT. T) THEN
                  IX(IJ) = IX(I)
                  IX(I) = T
                  T = IX(IJ)
            END IF
            L = J
            !     
            !     If last element of array is less than than T, interchange with T
            !     
            IF (IX(J) .LT. T) THEN
                  IX(IJ) = IX(J)
                  IX(J) = T
                  T = IX(IJ)
                  !     
                  !     If first element of array is greater than T, interchange with T
                  !     
                  IF (IX(I) .GT. T) THEN
                        IX(IJ) = IX(I)
                        IX(I) = T
                        T = IX(IJ)
                  END IF
            END IF
            !     
            !     Find an element in the second half of the array which is smaller
            !     than T
            !     
40          L = L-1
            IF (IX(L) .GT. T) GO TO 40
            !     
            !     Find an element in the first half of the array which is greater
            !     than T
            !     
50          K = K+1
            IF (IX(K) .LT. T) GO TO 50
            !     
            !     Interchange these elements
            !     
            IF (K .LE. L) THEN
                  TT = IX(L)
                  IX(L) = IX(K)
                  IX(K) = TT
                  GO TO 40
            END IF
            !     
            !     Save upper and lower subscripts of the array yet to be sorted
            !     
            IF (L-I .GT. J-K) THEN
                  IL(M) = I
                  IU(M) = L
                  I = K
                  M = M+1
            ELSE
                  IL(M) = K
                  IU(M) = J
                  J = L
                  M = M+1
            END IF
            GO TO 70
            !     
            !     Begin again on another portion of the unsorted array
            !     
60          M = M-1
            IF (M .EQ. 0) GO TO 190
            I = IL(M)
            J = IU(M)

70          IF (J-I .GE. 1) GO TO 30
            IF (I .EQ. 1) GO TO 20
            I = I-1

80          I = I+1
            IF (I .EQ. J) GO TO 60
            T = IX(I+1)
            IF (IX(I) .LE. T) GO TO 80
            K = I

90          IX(K+1) = IX(K)
            K = K-1
            IF (T .LT. IX(K)) GO TO 90
            IX(K+1) = T
            GO TO 80
            !     
            !     Sort IX and carry IY along
            !     
100         M = 1
            I = 1
            J = NN
            R = 0.375E0

110         IF (I .EQ. J) GO TO 150
            IF (R .LE. 0.5898437D0) THEN
                  R = R+3.90625D-2
            ELSE
                  R = R-0.21875D0
            END IF

120         K = I
            !     
            !     Select a central element of the array and save it in location T
            !     
            IJ = I + INT((J-I)*R)
            T = IX(IJ)
            TY = IY(IJ)
            !     
            !     If first element of array is greater than T, interchange with T
            !     
            IF (IX(I) .GT. T) THEN
                  IX(IJ) = IX(I)
                  IX(I) = T
                  T = IX(IJ)
                  IY(IJ) = IY(I)
                  IY(I) = TY
                  TY = IY(IJ)
            END IF
            L = J
            !     
            !     If last element of array is less than T, interchange with T
            !     
            IF (IX(J) .LT. T) THEN
                  IX(IJ) = IX(J)
                  IX(J) = T
                  T = IX(IJ)
                  IY(IJ) = IY(J)
                  IY(J) = TY
                  TY = IY(IJ)
                  !     
                  !     If first element of array is greater than T, interchange with T
                  !     
                  IF (IX(I) .GT. T) THEN
                        IX(IJ) = IX(I)
                        IX(I) = T
                        T = IX(IJ)
                        IY(IJ) = IY(I)
                        IY(I) = TY
                        TY = IY(IJ)
                  END IF
            END IF
            !     
            !     Find an element in the second half of the array which is smaller
            !     than T
            !     
130         L = L-1
            IF (IX(L) .GT. T) GO TO 130
            !     
            !     Find an element in the first half of the array which is greater
            !     than T
            !     
140         K = K+1
            IF (IX(K) .LT. T) GO TO 140
            !     
            !     Interchange these elements
            !     
            IF (K .LE. L) THEN
                  TT = IX(L)
                  IX(L) = IX(K)
                  IX(K) = TT
                  TTY = IY(L)
                  IY(L) = IY(K)
                  IY(K) = TTY
                  GO TO 130
            END IF
            !     
            !     Save upper and lower subscripts of the array yet to be sorted
            !     
            IF (L-I .GT. J-K) THEN
                  IL(M) = I
                  IU(M) = L
                  I = K
                  M = M+1
            ELSE
                  IL(M) = K
                  IU(M) = J
                  J = L
                  M = M+1
            END IF
            GO TO 160
            !     
            !     Begin again on another portion of the unsorted array
            !     
150         M = M-1
            IF (M .EQ. 0) GO TO 190
            I = IL(M)
            J = IU(M)

160         IF (J-I .GE. 1) GO TO 120
            IF (I .EQ. 1) GO TO 110
            I = I-1

170         I = I+1
            IF (I .EQ. J) GO TO 150
            T = IX(I+1)
            TY = IY(I+1)
            IF (IX(I) .LE. T) GO TO 170
            K = I

180         IX(K+1) = IX(K)
            IY(K+1) = IY(K)
            K = K-1
            IF (T .LT. IX(K)) GO TO 180
            IX(K+1) = T
            IY(K+1) = TY
            GO TO 170
            !     
            !     Clean up
            !     
190         IF (KFLAG .LE. -1) THEN
                  DO I=1,NN
                        IX(I) = -IX(I)
                  END DO
            END IF
      END SUBROUTINE ISORT
end module sort
