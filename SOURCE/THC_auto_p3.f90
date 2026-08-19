module THC_auto_p3

      implicit none


contains


      subroutine p3a1_II_oooo_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = x0, min(x1, b)
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(b, x)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(y, a)
                              rs = posT(b, x)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a1_II_oooo_ps_qsrp

      subroutine p3a1_II_aooo_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(a) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a1_II_aooo_ps_qsrp

      subroutine p3a1_AI_ooao_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(x) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a1_AI_ooao_ps_qsrp

      subroutine p3a1_AI_aoao_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(x) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a1_AI_aoao_ps_qsrp

      subroutine p3a1_AI_ooaa_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(y, x)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a1_AI_ooaa_ps_qsrp

      subroutine p3a1_AI_aoaa_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(x) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a1_AI_aoaa_ps_qsrp

      subroutine p3a1_IA_aaoo_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(b, a)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(b, a)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a1_IA_aaoo_ps_qsrp

      subroutine p3a1_IA_vaoo_ps_qsrp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(a) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a1_IA_vaoo_ps_qsrp

      subroutine p3a2_II_oooo_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a2_II_oooo_qr_prsq

      subroutine p3a2_II_ooao_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a2_II_ooao_qr_prsq

      subroutine p3a2_AI_ooaa_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a2_AI_ooaa_qr_prsq

      subroutine p3a2_AI_oova_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a2_AI_oova_qr_prsq

      subroutine p3a2_IA_aooo_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a2_IA_aooo_qr_prsq

      subroutine p3a2_IA_aaoo_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a2_IA_aaoo_qr_prsq

      subroutine p3a2_IA_aoao_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a2_IA_aoao_qr_prsq

      subroutine p3a2_IA_aaao_qr_prsq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a2_IA_aaao_qr_prsq

      subroutine p3a3_II_oooo_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = x0, min(x1, b)
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(b, x)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(a, y)
                              rs = posT(b, x)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a3_II_oooo_qs_psrq

      subroutine p3a3_II_aooo_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_II_aooo_qs_psrq

      subroutine p3a3_II_ooao_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_II_ooao_qs_psrq

      subroutine p3a3_II_aoao_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(b) * Occ(a)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a3_II_aoao_qs_psrq

      subroutine p3a3_AI_ooaa_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_AI_ooaa_qs_psrq

      subroutine p3a3_AI_aoaa_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(x) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(y, x)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_AI_aoaa_qs_psrq

      subroutine p3a3_AI_oova_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_AI_oova_qs_psrq

      subroutine p3a3_AI_aova_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_AI_aova_qs_psrq

      subroutine p3a3_IA_aaoo_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_IA_aaoo_qs_psrq

      subroutine p3a3_IA_vaoo_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_IA_vaoo_qs_psrq

      subroutine p3a3_IA_aaao_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_IA_aaao_qs_psrq

      subroutine p3a3_IA_vaao_qs_psrq(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a3_IA_vaao_qs_psrq

      subroutine p3a4_II_oooo_pr_qrsp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a4_II_oooo_pr_qrsp

      subroutine p3a4_AI_ooao_pr_qrsp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(x) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a4_AI_ooao_pr_qrsp

      subroutine p3a4_AI_ooaa_pr_qrsp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(y, x)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a4_AI_ooaa_pr_qrsp

      subroutine p3a4_IA_aooo_pr_qrsp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3a4_IA_aooo_pr_qrsp

      subroutine p3a4_IA_aaoo_pr_qrsp(ASing, ATrip, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3a4_IA_aaoo_pr_qrsp

      subroutine p3c1_II_oooo_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_II_oooo_qs_prqs

      subroutine p3c1_II_aooo_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_II_aooo_qs_prqs

      subroutine p3c1_II_ooao_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(b) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c1_II_ooao_qs_prqs

      subroutine p3c1_II_aoao_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_II_aoao_qs_prqs

      subroutine p3c1_AI_ooaa_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_AI_ooaa_qs_prqs

      subroutine p3c1_AI_aoaa_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(b)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c1_AI_aoaa_qs_prqs

      subroutine p3c1_AI_oova_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_AI_oova_qs_prqs

      subroutine p3c1_AI_aova_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(b)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c1_AI_aova_qs_prqs

      subroutine p3c1_IA_aaoo_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_IA_aaoo_qs_prqs

      subroutine p3c1_IA_vaoo_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_IA_vaoo_qs_prqs

      subroutine p3c1_IA_aaao_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c1_IA_aaao_qs_prqs

      subroutine p3c1_IA_vaao_qs_prqs(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(b) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c1_IA_vaao_qs_prqs

      subroutine p3c3_II_oooo_pr_qspr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(b, a)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(b, a)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c3_II_oooo_pr_qspr

      subroutine p3c3_AI_ooao_pr_qspr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(x) * Occ(a)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) + val

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c3_AI_ooao_pr_qspr

      subroutine p3c3_AI_ooaa_pr_qspr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(y, x)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c3_AI_ooaa_pr_qspr

      subroutine p3c3_IA_aooo_pr_qspr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(a) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c3_IA_aooo_pr_qspr

      subroutine p3c3_IA_aaoo_pr_qspr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(b, a)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) + val

                                    pq = posT(b, a)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c3_IA_aaoo_pr_qspr

      subroutine p3c4_II_oooo_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c4_II_oooo_ps_qrps

      subroutine p3c4_II_aooo_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c4_II_aooo_ps_qrps

      subroutine p3c4_AI_ooao_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(x) * Occ(a)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c4_AI_ooao_ps_qrps

      subroutine p3c4_AI_aoao_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(b)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c4_AI_aoao_ps_qrps

      subroutine p3c4_AI_ooaa_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(y, x)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c4_AI_ooaa_ps_qrps

      subroutine p3c4_AI_aoaa_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(b)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c4_AI_aoaa_ps_qrps

      subroutine p3c4_IA_aaoo_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c4_IA_aaoo_ps_qrps

      subroutine p3c4_IA_vaoo_ps_qrps(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c4_IA_vaoo_ps_qrps

      subroutine p3c5_II_oooo_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_II_oooo_qr_psqr

      subroutine p3c5_II_ooao_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(b) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c5_II_ooao_qr_psqr

      subroutine p3c5_AI_ooaa_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_AI_ooaa_qr_psqr

      subroutine p3c5_AI_oova_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(b)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_AI_oova_qr_psqr

      subroutine p3c5_IA_aooo_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(x)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_IA_aooo_qr_psqr

      subroutine p3c5_IA_aaoo_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_IA_aaoo_qr_psqr

      subroutine p3c5_IA_aoao_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(b) * Occ(y)* V_axby(x,y)        
                              ASing(rs, pq) = ASing(rs, pq) - val

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATrip(rs, pq) = ATrip(rs, pq) + val
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3c5_IA_aoao_qr_psqr

      subroutine p3c5_IA_aaao_qr_psqr(ASing, ATrip, ATripA, a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(y)* V_axby(x,y)        
                                    ASing(rs, pq) = ASing(rs, pq) - val

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATrip(rs, pq) = ATrip(rs, pq) + val
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3c5_IA_aaao_qr_psqr

      !------------------------------------------------
            subroutine p3b1_II_oooo_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = x0, min(x1, b)
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(b, x)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(y, a)
                              rs = posT(b, x)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b1_II_oooo_ps_qsrp

      subroutine p3b1_II_aooo_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(a) * Occ(y)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b1_II_aooo_ps_qsrp

      subroutine p3b1_AI_ooao_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(x) * Occ(y)* V_axby(x,y)        

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b1_AI_ooao_ps_qsrp

      subroutine p3b1_AI_aoao_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(x, b)
                        rs = posS(a, y)
                        if (pq>0.and.rs>0)then
                              val = Occ(x) * Occ(y)* V_axby(x,y)        

                              pq = posT(x, b)
                              rs = posT(a, y)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b1_AI_aoao_ps_qsrp

      subroutine p3b1_AI_ooaa_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(y, x)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b1_AI_ooaa_ps_qsrp

      subroutine p3b1_AI_aoaa_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(x) * Occ(b)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b1_AI_aoaa_ps_qsrp

      subroutine p3b1_IA_aaoo_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(b, a)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(b, a)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b1_IA_aaoo_ps_qsrp

      subroutine p3b1_IA_vaoo_ps_qsrp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(a) * Occ(y)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b1_IA_vaoo_ps_qsrp

      subroutine p3b2_II_oooo_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b2_II_oooo_qr_prsq

      subroutine p3b2_II_ooao_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b2_II_ooao_qr_prsq

      subroutine p3b2_AI_ooaa_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b2_AI_ooaa_qr_prsq

      subroutine p3b2_AI_oova_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b2_AI_oova_qr_prsq

      subroutine p3b2_IA_aooo_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b2_IA_aooo_qr_prsq

      subroutine p3b2_IA_aaoo_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b2_IA_aaoo_qr_prsq

      subroutine p3b2_IA_aoao_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(a, y)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b2_IA_aoao_qr_prsq

      subroutine p3b2_IA_aaao_qr_prsq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(a, b)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b2_IA_aaao_qr_prsq

      subroutine p3b3_II_oooo_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = x0, min(x1, b)
                  do y = y0, min(y1, a)

                        pq = posS(a, y)
                        rs = posS(b, x)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(a, y)
                              rs = posT(b, x)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b3_II_oooo_qs_psrq

      subroutine p3b3_II_aooo_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_II_aooo_qs_psrq

      subroutine p3b3_II_ooao_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_II_ooao_qs_psrq

      subroutine p3b3_II_aoao_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(b) * Occ(a)* V_axby(x,y)        

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b3_II_aoao_qs_psrq

      subroutine p3b3_AI_ooaa_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_AI_ooaa_qs_psrq

      subroutine p3b3_AI_aoaa_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(x) * Occ(b)* V_axby(x,y)        

                                    pq = posT(y, x)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_AI_aoaa_qs_psrq

      subroutine p3b3_AI_oova_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_AI_oova_qs_psrq

      subroutine p3b3_AI_aova_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(x, y)
                              rs = posS(b, a)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(x, y)
                                    rs = posT(b, a)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_AI_aova_qs_psrq

      subroutine p3b3_IA_aaoo_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_IA_aaoo_qs_psrq

      subroutine p3b3_IA_vaoo_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_IA_vaoo_qs_psrq

      subroutine p3b3_IA_aaao_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(a, b)
                              rs = posS(y, x)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(a, b)
                                    rs = posT(y, x)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_IA_aaao_qs_psrq

      subroutine p3b3_IA_vaao_qs_psrq(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = y0, y1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(a) * Occ(y)* V_axby(x,y)        

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b3_IA_vaao_qs_psrq

      subroutine p3b4_II_oooo_pr_qrsp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b4_II_oooo_pr_qrsp

      subroutine p3b4_AI_ooao_pr_qrsp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(x, b)
                        rs = posS(y, a)
                        if (pq>0.and.rs>0)then
                              val = -Occ(x) * Occ(y)* V_axby(x,y)        

                              pq = posT(x, b)
                              rs = posT(y, a)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b4_AI_ooao_pr_qrsp

      subroutine p3b4_AI_ooaa_pr_qrsp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (b<=a)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do x = x0, x1
                        do y = max(y0, x), y1


                              pq = posS(y, x)
                              rs = posS(a, b)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(y) * Occ(a)* V_axby(x,y)        

                                    pq = posT(y, x)
                                    rs = posT(a, b)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b4_AI_ooaa_pr_qrsp

      subroutine p3b4_IA_aooo_pr_qrsp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val


!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
            do x = max(x0, b), x1
                  do y = max(y0, a), y1

                        pq = posS(y, a)
                        rs = posS(x, b)
                        if (pq>0.and.rs>0)then
                              val = -Occ(y) * Occ(x)* V_axby(x,y)        

                              pq = posT(y, a)
                              rs = posT(x, b)
                              if (pq>0.and.rs>0)then
                                    ATripA(rs, pq) = ATripA(rs, pq) + val
                              end if
                        end if
                  end do
            end do
!!$omp end do
!!$omp end parallel

      end subroutine p3b4_IA_aooo_pr_qrsp

      subroutine p3b4_IA_aaoo_pr_qrsp(ATripA,  a, b, x0, x1, y0, y1, posS, posT, Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: ATripA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, dimension(:,:), intent(in) :: posS, posT
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs

            double precision :: val

            if (a<=b)then
!!$omp parallel default(shared)
!!$omp private(x, y, pq, rs, val)
!!$omp do collapse(2)
                  do y = y0, y1
                        do x = max(x0, y), x1


                              pq = posS(b, a)
                              rs = posS(x, y)
                              if (pq>0.and.rs>0)then
                                    val = -Occ(b) * Occ(x)* V_axby(x,y)        

                                    pq = posT(b, a)
                                    rs = posT(x, y)
                                    if (pq>0.and.rs>0)then
                                          ATripA(rs, pq) = ATripA(rs, pq) + val
                                    end if
                              end if
                        end do
                  end do
!!$omp end do
!!$omp end parallel
            end if
      end subroutine p3b4_IA_aaoo_pr_qrsp



end module THC_auto_p3
