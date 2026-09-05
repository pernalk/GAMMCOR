module thc_auto_p0
      use math_constants
      use omp_lib
      implicit none

contains


      subroutine bare_int_prqs_vvvv_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, ind_virt, ndim_virt, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1, ndim_virt
            integer, dimension(:,:), intent(in) :: posS, posT, ind_virt
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, i

            double precision :: val1, val2

            !$omp parallel default(shared) private(y, x, i, pq, rs, val1, val2)            
            !$omp do

            do i = 1, ndim_virt
                  x = ind_virt(1, i)
                  y = ind_virt(2, i)

                  ! if (a >=b)then
                  !       do y = y0, y1
                  !             do x = x0, x1   
                  pq = posS(a, b)
                  rs = posS(x, y)
                  val1 = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y)) * V_axby(x,y)
                  val2 = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y)) * V_axby(y,x)
                  if (pq>0.and.rs>0)then

                        ASing(rs, pq) = ASing(rs, pq) + val1+val2
                        pq = posT(a, b)
                        rs = posT(x, y) 
                        if (pq>0.and.rs>0)then
                              ATrip(rs, pq) = ATrip(rs, pq) +val1 -val2
                              ATripA(rs, pq) = ATripA(rs, pq) +val1 -val2
                        end if
                  end if
            end do
            !$omp end do
            !$omp end parallel
            !end do
            !end if


      end subroutine bare_int_prqs_vvvv_prqs


      subroutine bare_int_prqs_oooo_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_oooo_prqs

      subroutine bare_int_prqs_aooo_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aooo_prqs

      subroutine bare_int_prqs_ooao_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, min(y1, a)
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(a, y) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(a, y) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_ooao_prqs

      subroutine bare_int_prqs_aaoo_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aaoo_prqs

      subroutine bare_int_prqs_ooaa_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(a, b) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(a, b) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_ooaa_prqs

      subroutine bare_int_prqs_vaoo_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_vaoo_prqs

      subroutine bare_int_prqs_oova_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(x, y)
                              rs = posS(a, b) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(a, b) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_oova_prqs

      subroutine bare_int_prqs_vvoo_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_vvoo_prqs

      subroutine bare_int_prqs_oovv_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(a, b) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(a, b) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_oovv_prqs

      subroutine bare_int_prqs_aoao_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aoao_prqs

      subroutine bare_int_prqs_aaao_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aaao_prqs

      subroutine bare_int_prqs_aoaa_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, min(y1, a)
                  do x = x0, x1   
                        pq = posS(a, y)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_aoaa_prqs

      subroutine bare_int_prqs_aaaa_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aaaa_prqs

      subroutine bare_int_prqs_vaao_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, x1
                        pq = posS(x, b)
                        rs = posS(a, y) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(a, y) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_vaao_prqs

      subroutine bare_int_prqs_aova_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, x1   
                        pq = posS(a, y)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_aova_prqs

      subroutine bare_int_prqs_vvao_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, min(y1, a)
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(a, y) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(a, y) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_vvao_prqs

      subroutine bare_int_prqs_aovv_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, min(y1, a)
                  do x = x0, x1   
                        pq = posS(a, y)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))

                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_aovv_prqs

      subroutine bare_int_prqs_vaaa_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = max(x0, b), x1
                        pq = posS(x, b)
                        rs = posS(a, y) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))

                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(a, y) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_vaaa_prqs

      subroutine bare_int_prqs_aava_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aava_prqs

      subroutine bare_int_prqs_vvaa_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(a, b) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(a, b) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_vvaa_prqs

      subroutine bare_int_prqs_aavv_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = max(x0, y), x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_aavv_prqs

      subroutine bare_int_prqs_vava_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_vava_prqs

      subroutine bare_int_prqs_vvva_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = max(x0, b), x1   
                        pq = posS(a, y)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_prqs_vvva_prqs

      subroutine bare_int_prqs_vavv_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(a, b)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) +val *  V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) +val *  V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_prqs_vavv_prqs


      subroutine bare_int_psqr_oooo_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_oooo_psqr

      subroutine bare_int_psqr_aooo_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(b, a)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(b, a)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aooo_psqr

      subroutine bare_int_psqr_ooao_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = max(y0, a), y1
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(y, a) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(y, a) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_ooao_psqr

      subroutine bare_int_psqr_aaoo_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aaoo_psqr

      subroutine bare_int_psqr_aoao_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(y, a) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(y, a) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_aoao_psqr

      subroutine bare_int_psqr_ooaa_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(b, a) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(b, a) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_ooaa_psqr

      subroutine bare_int_psqr_vaoo_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(b, a)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(b, a)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_vaoo_psqr

      subroutine bare_int_psqr_oova_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(x, y)
                              rs = posS(b, a) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(b, a) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_oova_psqr

      subroutine bare_int_psqr_vvoo_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_vvoo_psqr

      subroutine bare_int_psqr_oovv_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(b, a) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(b, a) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_oovv_psqr

      subroutine bare_int_psqr_aoaa_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = max(y0, a), y1
                  do x = x0, x1   
                        pq = posS(y, a)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(y, a)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_aoaa_psqr

      subroutine bare_int_psqr_aaao_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aaao_psqr

      subroutine bare_int_psqr_vaao_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(b, a)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(b, a)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_vaao_psqr

      subroutine bare_int_psqr_aova_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(y, a) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(y, a) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_aova_psqr

      subroutine bare_int_psqr_aaaa_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aaaa_psqr

      subroutine bare_int_psqr_vaaa_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, min(x1, b)   
                        pq = posS(b, x)
                        rs = posS(a, y) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(b, x)
                              rs = posT(a, y) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_vaaa_psqr

      subroutine bare_int_psqr_aava_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aava_psqr

      subroutine bare_int_psqr_vvao_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = max(y0, a), y1
                  do x = x0, x1   
                        pq = posS(x, b)
                        rs = posS(y, a) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(x, b)
                              rs = posT(y, a) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_vvao_psqr

      subroutine bare_int_psqr_aovv_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = max(y0, a), y1
                  do x = x0, x1   
                        pq = posS(y, a)
                        rs = posS(x, b) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(y, a)
                              rs = posT(x, b) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_aovv_psqr

      subroutine bare_int_psqr_vvaa_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = max(x0,y), x1   
                              pq = posS(x, y)
                              rs = posS(b, a) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(x, y)
                                    rs = posT(b, a) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_vvaa_psqr

      subroutine bare_int_psqr_vava_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, x1   
                        pq = posS(a, y)
                        rs = posS(b, x) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(b, x) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_vava_psqr

      subroutine bare_int_psqr_aavv_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a >=b)then
                  do y = y0, y1
                        do x = x0, min(x1, y)   
                              pq = posS(a, b)
                              rs = posS(y, x) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(a, b)
                                    rs = posT(y, x) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_aavv_psqr

      subroutine bare_int_psqr_vavv_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            if (a <=b)then
                  do y = y0, y1
                        do x = x0, x1   
                              pq = posS(b, a)
                              rs = posS(x, y) 
                              if (pq>0.and.rs>0)then
                                    val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                                    ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                                    pq = posT(b, a)
                                    rs = posT(x, y) 
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                          ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                                    end if
                              end if
                        end do
                  end do
            end if

      end subroutine bare_int_psqr_vavv_psqr

      subroutine bare_int_psqr_vvva_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val
            do y = y0, y1
                  do x = x0, min(x1, b)   
                        pq = posS(a, y)
                        rs = posS(b, x) 
                        if (pq>0.and.rs>0)then
                              val = (One -Occ(a)-Occ(b)-Occ(x)-Occ(y))
                              ASing(rs, pq) = ASing(rs, pq) + val * V_axby(x,y)
                              pq = posT(a, y)
                              rs = posT(b, x) 
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) - val * V_axby(x,y)
                                    ATripA(rs, pq) = ATripA(rs, pq) - val * V_axby(x,y)
                              end if
                        end if
                  end do
            end do

      end subroutine bare_int_psqr_vvva_psqr




end module thc_auto_p0
