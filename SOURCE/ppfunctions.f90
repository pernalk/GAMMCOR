module ppfunctions


      use iso_fortran_env
      use real_linalg !from gammcor integrals 
      use types
      use lin
      use clock

      use math_constants
      use quadratures

      implicit none

      integer, parameter, private :: MSG_TIME  = 0
      integer, parameter, private :: MSG_VERB = 10
      integer, parameter, private :: MSG_NOR = 50
      integer, parameter, private :: MSG_ERR = 100

      integer, parameter :: MSG_DEB  = 0

      integer, save :: MSG_PRIORITY_THR = 1


contains




      subroutine prij(st, i, j)
            character(len=*), intent(in) :: st
            integer, intent(in) :: i, j

            write(*, '(A10, 2I5)') st, i, j
      end subroutine prij

      subroutine pri3(st, i, j, k)
            character(len=*), intent(in) :: st
            integer, intent(in) :: i, j, k

            write(*, '(A10, 3I5)') st, i, j, k
      end subroutine pri3

      subroutine pr4ix(st, i, j, k, l, a)
            character(len=*), intent(in) :: st
            integer, intent(in) :: i, j, k, l
            double precision, intent(in) :: a

            write(*, '(A10, 4I5, F20.15)') st, i, j, k, l, a
      end subroutine pr4ix


      subroutine pr2ix(st, i, j, a)
            character(len=*), intent(in) :: st
            integer, intent(in) :: i, j
            double precision, intent(in) :: a

            write(*, '(A10, 2I5, F20.15)') st, i, j,  a
      end subroutine pr2ix

      subroutine prix(st, i, a)
            character(len=*), intent(in) :: st
            integer, intent(in) :: i
            double precision, intent(in) :: a

            write(*, '(A10, I5, F25.15)') st, i, a
      end subroutine prix


      subroutine tmsg(s, t, priority)

            character(len=*), intent(in) :: s
            type(tclock), intent(in)     :: t
            integer, optional            :: priority

            character(len=70) :: line
            character(len=40) :: title
            real(F64) :: d
            integer :: p

            if (present(priority)) then
                  p = priority
            else
                  p = MSG_NOR
            end if


            if (p==MSG_TIME)then
                  ! Compute the elapsed time
                  d = clock_readwall(t)

                  ! Prepare the line
                  line = repeat('-', 70)

                  ! Print the formatted output
                  print *, ''
                  print *, line
                  write(title, '(A)') trim(s)
                  write(STDOUNIT, '(1X,A40,F20.3,A8)') adjustl(title), d, '  [sec]'
                  print *, line
                  print *, ''
                  flush(STDOUNIT)
            end if

      end subroutine tmsg


      subroutine xmsg(s, t, priority)

            character(len=*), intent(in) :: s
            real(F64), intent(in) :: t
            integer, optional            :: priority

            character(len=40) :: title
            integer :: p

            if (present(priority)) then
                  p = priority
            else
                  p = MSG_NOR
            end if


            if (p>=MSG_PRIORITY_THR)then

                  print *, ''                  
                  write(title, '(A)') trim(s)
                  write(STDOUNIT, '(1X,A40,A1, F30.15)') adjustl(title), ',', t
                  flush(STDOUNIT)
            end if

      end subroutine xmsg


      subroutine ximsg(s, t, priority)

            character(len=*), intent(in) :: s
            integer, intent(in) :: t
            integer, optional            :: priority

            character(len=40) :: title
            integer :: p

            if (present(priority)) then
                  p = priority
            else
                  p = MSG_NOR
            end if


            if (p>=MSG_PRIORITY_THR)then

                  print *, ''                  
                  write(title, '(A)') trim(s)
                  write(STDOUNIT, '(1X,A40,I15)') adjustl(title), t
                  flush(STDOUNIT)
            end if

      end subroutine ximsg

      subroutine mmsg(s, priority)

            character(len=*), intent(in) :: s
            integer, optional            :: priority

            character(len=40) :: title
            integer :: p

            if (present(priority)) then
                  p = priority
            else
                  p = MSG_NOR
            end if


            if (p>=MSG_PRIORITY_THR)then

                  print *, ''                  
                  write(title, '(A)') trim(s)
                  write(STDOUNIT, '(1X,A40)') adjustl(title)
                  flush(STDOUNIT)
            end if

      end subroutine mmsg



      subroutine orthogonalize_degen(n, countj, wr, vr, S_diag, dy, x)

            integer, intent(in) :: n
            double precision, dimension(:), intent(in) :: wr
            double precision, dimension(:,:), intent(inout) :: vr
            double precision, dimension(:), intent(in) :: S_diag
            integer, dimension(:), intent(in) :: dy
            integer, intent(in) :: x
            integer, dimension(:), allocatable :: StartIdx, EndIdx
            integer, intent(in) :: countj
            double precision, parameter :: tol = 1.d-4
            integer :: count, j, i


            allocate(StartIdx(n))
            allocate(EndIdx(n))

            StartIdx = 0
            EndIdx = 0
            count = 1
            j = 0

            do i = 1, countj
                  !       print*, 'teraz', i, countj, wr(i), wr(dy(i))                                                                                                                                                                                                                                                                                                                       
                  if (j==0)then
                        !         print*, 'zaczynam liczyć od', i, wr(i)                                                                                                                                                                                                                                                                                                                     
                        StartIdx(count) = i
                        if (i.eq.countj)then
                              EndIdx(count) = i
                        end if
                        j = 1
                  else
                        if (abs(wr(i)-wr(i-1)).lt.tol)then
                              !           print*, 'kontynuuje zliaczanie dla', i, countj, wr(i)                                                                                                                                                                                                                                                                                              
                              if (i==n)then
                                    EndIdx(count) = i
                              end if
                              if (i==countj)then
                                    EndIdx(count) = i
                              end if
                        else
                              !             print*, 'ten juz jest inny', wr(i), 'wiec koncze poprzednim', i-1                                                                                                                                                                                                                                                                                
                              EndIdx(count) = i-1

                              if (i.ne.countj)then
                                    !                print*, 'zaczynam dalej liczyc od', wr(i), i, countj                                                                                                                                                                                                                                                                                    
                                    count = count + 1
                                    StartIdx(count)=i
                                    !                print*, 'i mam count, i', count, i                                                                                                                                                                                                                                                                                                      
                              else
                                    count = count+1
                                    StartIdx(count) = i
                                    EndIdx(count) = i
                              end if

                        end if
                  end if
            end do


            if (StartIdx(1) == 0)then
                  count = 0
            end if
            ! print*, 'startstrop'                                                                                                                                                                                                                                                                                                                                                           
            !     do i =1, count                                                                                                                                                                                                                                                                                                                                                             
            !        print*, StartIdx(i), EndIdx(i)                                                                                                                                                                                                                                                                                                                                          
            !     end do                                                                                                                                                                                                                                                                                                                                                                     

            call Orthogonalize(vr, count, StartIdx(1:count), EndIdx(1:count), S_diag, dy, x)
      end subroutine orthogonalize_degen


      function P3a(s, p, q, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P3a
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, tt, uu
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3

            NOc = NI+NA
            P3a = zero

            if (IAux(s)==(1).and.IAux(p)==(1))then

                  do tt = NI+1, NI+NA
                        do uu = NI + 1, NI + NA
                              t = IMod(tt)
                              u = IMod(uu)
                              P3a = P3a + TwoNO(NAddr3(q,u,r,t))* &
                                    Gam(t, s, p, u, R00, R11, n, I2, true_NA, 1)
                        end do
                  end do
            else
                  if ((IAux(s)==(1).and.IAux(p)==(0)).or.&
                        (IAux(s)==(0).and.IAux(p)==(1)).or.&
                        (IAux(s)==(0).and.IAux(p)==(0)))then
                        P3a = P3a + n(p)*n(s)*TwoNO(NAddr3(q,s,r,p))
                  end if
            end if

      end function P3a


      function P3b(s, q, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P3b
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, tt, uu
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3

            NOc = NI+NA 
            P3b = zero
            if (IAux(s)==(1).and.IAux(q)==(1))then

                  do tt = NI+1, NI+NA
                        do uu = NI + 1, NI + NA
                              t = IMod(tt)
                              u = IMod(uu)
                              P3b = P3b + TwoNO(NAddr3(p,u,r,t))* &
                                    Gam(t, s, u, q, R00, R11, n, I2, true_NA, 1)
                        end do
                  end do

                  if (s==q)then  !s=q, active, t inactive
                        do tt = 1, NI!NI+NA!NBasis
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P3b = P3b + n(t)*n(q)*TwoNO(NAddr3(p,t,r,t))
                              !            end if
                        end do
                  end if

            else
                  if (s==q)then
                        if (IAux(s).ne.1)then !s=q, inactive, t act + inact
                              do tt = 1, NI+NA!NBasis!NI+NA
                                    t = IMod(tt)
                                    !               if	(IAux(t)==0)then
                                    P3b = P3b + TwoNO(NAddr3(p,t,r,t))* n(t)*n(q)
                                    !              end if
                              end do
                        end if
                  end if
            end if

      end function P3b

      function P3c(s, q, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P3c
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, tt, uu
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3

            NOc = NI+NA
            P3c = zero
            if (IAux(s)==(1).and.IAux(q)==(1))then
                  do tt = NI+1, NI+NA
                        do uu = NI+1, NI+NA
                              t = IMod(tt)
                              u = IMod(uu)
                              P3c = P3c +  TwoNO(NAddr3(p,r,t,u))* &
                                    (Gam(t, s, u, q, R00, R11, n, I2, true_NA, 0)+Gam(t, s, u, q, R00, R11, n, I2, true_NA, 1))
                        end do
                  end do

                  if (s==q)then
                        do tt = 1, NI!NI+NA!NBasis!NI
                              t = IMod(tt)
                              !if (IAux(t)==0)then
                              P3c = P3c + two* n(t)*n(s) * TwoNO(NAddr3(p,r,t,t))
                              ! end if
                        end do
                  end if
            else
                  if (s==q)then
                        if (IAux(s)==0)then
                              do tt = 1, NI+NA!NBasis!NI+NA
                                    t = IMod(tt)
                                    !                if	(IAux(t)==0)then
                                    P3c = P3c + two* n(t)*n(s) * TwoNO(NAddr3(p,r,t,t))
                                    !               end if
                              end do
                        end if
                  end if
            end if

            if ((IAux(s)==0.and.IAux(q)==1).or.&
                  (IAux(s)==1.and.IAux(q)==0).or.&
                  (IAux(s)==0.and.IAux(q)==0))then
                  P3c = P3c - n(q)*n(s) * TwoNO(NAddr3(p,r,q,s))
                  !print*, 'p3c2', p3c
            end if
            !    print*, 'p3cwhat', P3c

      end function P3c


      function P4a(q, s, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P4a
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, v, tt, uu, vv
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3


            NOc = NI+NA
            P4a = zero
            if (p==r)then

                  if (IAux(q)==(1))then
                        do tt = NI+1, NI+NA
                              do uu = NI + 1, NI + NA
                                    do vv = NI + 1, NI + NA
                                          t = IMod(tt)
                                          u = IMod(uu)
                                          v = IMod(vv)
                                          P4a = P4a + frac12  * (TwoNO(NAddr3(s, t, u, v)) * &
                                                Gam(t, u, q, v, R00, R11, n, I2, true_NA, 1))
                                    end do
                              end do
                        end do

                        do tt = 1, NI!NI+NA!NBasis!NI
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P4a = P4a + frac12 * n(q) * n(t) * TwoNO(NAddr3(s, q, t, t))
                              !            end if
                        end do

                  else
                        do tt = 1, NI+NA!NBasis!NI+NA
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P4a = P4a + frac12 * n(q) * n(t) * TwoNO(NAddr3(s, q, t, t))
                              !            end if
                        end do

                  end if
            end if
      end function P4a

      function P4b(s, q, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P4b
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, v, tt, uu, vv
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3


            NOc = NI+NA
            P4b = zero
            if (p==r)then

                  if (IAux(s)==(1))then
                        do tt = NI+1, NI+NA
                              do uu = NI+1, NI+NA
                                    do vv = NI+1, NI+NA
                                          t = IMod(tt)
                                          u = IMod(uu)
                                          v = IMod(vv)
                                          P4b = P4b + frac12 * (TwoNO(NAddr3(q, v, t, u)) * &
                                                Gam(t, s, u, v, R00, R11, n, I2, true_NA, 1))
                                    end do
                              end do
                        end do

                        do tt = 1, NI!NI+NA!NBasis!NI
                              t = IMod(tt) 
                              !             if	(IAux(t)==0)then
                              P4b = P4b + frac12 * n(s) * n(t) * TwoNO(NAddr3(q, s, t, t))
                              !            end if
                        end do

                  else
                        do tt = 1, NI+NA!NBasis!NI+NA
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P4b = P4b + frac12 * n(s) * n(t) * TwoNO(NAddr3(q, s, t, t))
                              !            end if
                        end do
                  end if
            end if
      end function P4b


      function P5a(q, s, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P5a
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, v, tt, uu, vv
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3


            NOc = NI+NA
            P5a = zero
            if (p==r)then

                  if (IAux(q)==(1))then
                        do tt = NI+1, NI+NA
                              do uu = NI + 1, NI + NA
                                    do vv = NI + 1, NI + NA
                                          t = IMod(tt)
                                          u = IMod(uu)
                                          v = IMod(vv)
                                          P5a = P5a + frac14 * ((TwoNO(NAddr3(s, u, t, v)) - TwoNO(NAddr3(s, t, u, v)))&
                                                * Gam(t, u, q, v, R00, R11, n, I2, true_NA, 0))
                                    end do
                              end do
                        end do

                        do tt = 1, NI!NI+NA!NBasis!NI
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P5a = P5a + frac12 * n(t) * n(q) * (TwoNO(NAddr3(s, t, q, t)) -TwoNO(NAddr3(s, q, t, t)))
                              !            end if
                        end do

                  else
                        do tt = 1, NI+NA!NBasis!NI+NA
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P5a = P5a + frac12 * n(t) * n(q) * (TwoNO(NAddr3(s, t, q, t)) -TwoNO(NAddr3(s, q, t, t)))
                              !            end if
                        end do
                  end if
            end if
      end function P5a

      function P5b(s, q, p, r, TwoNO, R00, R11, n, IAux, IMod, I2, NI, NA, NBasis, true_NA)
            use math_constants
            double precision :: P5b
            integer, intent(in) :: r, s, p, q
            double precision, dimension(:), intent(in) :: TwoNO, R00, R11
            double precision, dimension(:), intent(in) :: n
            integer, dimension(:), intent(in) :: IAux, IMod, I2
            integer, intent(in) :: NI, NA, NBasis, true_NA
            integer :: t, u, v, tt, uu, vv
            integer :: NOc
            integer, external :: NAddrRDM, NAddr3


            NOc = NI+NA
            P5b = zero
            if (p==r)then

                  if (IAux(s)==(1))then
                        do tt = NI+1, NI+NA
                              do uu = NI + 1, NI + NA
                                    do vv = NI + 1, NI + NA
                                          t = IMod(tt)
                                          u = IMod(uu)
                                          v = IMod(vv)
                                          P5b = P5b + frac14 * ((TwoNO(NAddr3(q, u, t, v)) - TwoNO(NAddr3(q, v, t, u)))&
                                                * Gam(t, s, u, v, R00, R11, n, I2, true_NA, 0))
                                    end do
                              end do
                        end do

                        do tt = 1, NI!NI+NA!NBasis!NI
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P5b = P5b + frac12 * n(s) * n(t) * (TwoNO(NAddr3(q, t, t, s)) -TwoNO(NAddr3(q, s, t, t)))
                              !            end if
                        end do

                  else
                        do tt = 1, NI+NA!NBasis!NI+NA
                              t = IMod(tt)
                              !             if	(IAux(t)==0)then
                              P5b = P5b + frac12 * n(s) * n(t) * (TwoNO(NAddr3(q, t, t, s)) -TwoNO(NAddr3(q, s, t, t)))
                              !            end if
                        end do
                  end if
            end if

      end function P5b


      function Gam(p, q, r, s, R00, R11, n, Ind2, NA, pptype)
            use math_constants
            double precision :: Gam
            integer, intent(in) :: p, q, r, s
            double precision, dimension(:), intent(in) :: R00, R11
            integer, dimension(:), intent(in) :: Ind2
            double precision, dimension(:), intent(in) :: n
            integer, intent(in) :: NA
            integer, intent(in) :: pptype
            integer, external :: NAddrRDM
            double precision :: numf
            double precision :: gamp

            Gam = zero
            gamp = zero

            ! if (pseudo_hf)then
            !    if(p==r.and.q==s)then
            !       gamp = n(p) * n(q)
            !    end if
            !    if (pptype==0)then
            !       if(p==s.and.q==r)then
            !          gamp = Gam - n(p) * n(q)
            !       end if
            !    end if
            !    Gam = gamp
            ! else

            !       type = 0: compute Gam ++++
            !       type = 1: compute Gam +-+-                                                                                                                                   

            if (pptype==0)then
                  numf = one
            else if (pptype==1)then
                  numf = -one
            end if

            Gam = frac12 * (R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)) & 
                  + numf * R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)))

            ! end if

      end function Gam


      function get2RDM(p, q, r, s, R00, R11, Occ, Ind2, IndA, NA, pptype)

            use math_constants
            double precision :: get2RDM
            integer, intent(in) :: p, q, r, s
            double precision, dimension(:), intent(in) :: R00, R11
            integer, dimension(:), intent(in) :: Ind2, IndA
            double precision, dimension(:), intent(in) :: Occ
            integer, intent(in) :: NA
            integer, intent(in) :: pptype
            integer, external :: NAddrRDM
            double precision :: numf
            double precision :: temp

            get2RDM = zero
            temp = zero

            if (IndA(p)==IndA(q).and.IndA(q)==IndA(r).and.IndA(r)==IndA(s).and.IndA(p)==1)then

                  !type = 0: compute get2RDM ++++
                  !type = 1: compute get2RDM +-+-

                  if (pptype==0)then
                        numf = one
                  else if (pptype==1)then
                        numf = -one
                  end if

                  get2RDM = frac12 * (R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)) &
                        + numf * R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)))


                  ! if (abs(R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA))).gt.1.d-5)then
                  !    if (p.ne.q.and.q.ne.r.and.p.ne.r.and.p.ne.s.and.r.ne.s)then
                  !       if (pptype==1)then
                  !       write(*,'(5I5, 5F20.15)') pptype, p, q, r,s,  R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)), R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(r),Ind2(s), NA)), &
                  !            R00(NAddrRDM(Ind2(p),Ind2(q),Ind2(s),Ind2(r), NA)), R11(NAddrRDM(Ind2(p),Ind2(q),Ind2(s),Ind2(r), NA)), get2RDM
                  !    end if
                  !    end if
                  ! end if


            else

                  if(p==r.and.q==s)then
                        get2RDM = Occ(p) * Occ(q)
                  end if

                  if (pptype==0)then
                        if(p==s.and.q==r)then
                              get2RDM = get2RDM - Occ(p) * Occ(q)
                        end if
                  end if


            end if

      end function Get2RDM




      subroutine ddot_norm(vri, S, vrj, n, dd)
            double precision, dimension(:), intent(in) :: vri, vrj
            double precision, dimension(:, :), intent(in) :: S
            integer, intent(in) :: n
            double precision, intent(out) :: dd
            double precision, dimension(:), allocatable :: tempx

            allocate(tempx(n))
            tempx = zero
            call real_av_x(tempx, S, n,  vrj, n, n, one, zero)

            call real_vw_x(dd, vri, tempx, n)
            deallocate(tempx)

      end subroutine ddot_norm

      subroutine read_2rdm(filename, RR, NAct)
            character(len=*), intent(in) :: filename
            double precision, dimension(:), intent(out) :: RR
            integer, intent(in) :: NAct
            integer :: error1, error2, error3

            integer, parameter :: u2 = 20
            integer :: finito
            integer :: true_nof_lines
            integer :: i, j, k, l, ii
            double precision :: dm
            double precision, parameter :: Frac12=0.5+0
            integer, external :: NAddr3
            integer, external :: NAddrRDM

            open(u2, file=filename, status='old', access='sequential', form = 'formatted', action='read')

            finito = 0
            true_nof_lines = 0
            do while (finito==0)
                  read(u2,*, iostat=error2) i, j, k, l, dm
                  select case (error2)
                  case(0)
                        true_nof_lines = true_nof_lines + 1
                  case(iostat_end)
                        exit
                        finito = 1
                        exit
                  case Default

                        finito=1
                  end select
            end do
            close(u2)
            !            print*, true_nof_lines
            open(u2, file=filename, status='old', access='sequential', form = 'formatted', action='read')
            !           print*, 'true_nof_lines', true_nof_lines                                                                                                      
            do ii = 1, true_nof_lines
                  read(u2,*, iostat=error2) i, j, k, l, dm
                  !               write(*, '(A10, 4I5, F20.10)') 'sra', i, j, k, l, dm, nact                                                                                 
                  !              print*,'grr',  NAddrRDM(j,l,i,k, NAct)                                                                                                     
                  RR(NAddrRDM(j,l,i,k, NAct))=Frac12*dm
                  !               write(*, '(2I5, A6, 4I5, F20.10)') ii, NAddrRDM(j,l,i,k, NAct), '---', j, l, i, k, Frac12 * dm                                               
            end do

            close(u2)

      end subroutine read_2rdm

      function factorial_tab(i)
            double precision :: factorial_tab
            integer, intent(in) :: i
            double precision, dimension(13) :: ft

            ft(1) = 1.0d+0
            ft(2) = 2.0d+0
            ft(3) = 6.0d+0
            ft(4) = 24.0d+0
            ft(5) = 120.0d+0
            ft(6) = 720.0d+0
            ft(7) = 5040.0d+0
            ft(8) = 40320.0d+0
            ft(9) = 362880.0d+0
            ft(10) = 3628800.0d+0
            ft(11) =39916800.0d+0
            ft(12) =479001600.0d+0
            ft(13) =6227020800.0d+0

            if (i.le.13)then
                  factorial_tab = ft(i)
            else
                  print*, 'too large n for factorial'
                  stop
            end if

      end function factorial_tab






end module ppfunctions
