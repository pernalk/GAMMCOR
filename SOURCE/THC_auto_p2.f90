module THC_auto_p2
      use math_constants
      implicit none

contains



      subroutine p2_IIII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = 1, NI
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IIII_qstt_qstt

      subroutine p2_AAII_qstt_ttqs(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do y = 1, NI
                  Bux_AA(b, y) = Bux_AA(b, y) + Occ(t) * V_axby(t,y)
            end do
      end subroutine p2_AAII_qstt_ttqs

      subroutine p2_AAII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_AAII_qstt_qstt

      subroutine p2_IAII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IAII_qstt_qstt

      subroutine p2_IAII_qstt_sqtt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_II(x, a) = Bux_II(x, a) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IAII_qstt_sqtt

      subroutine p2_AAIA_qstt_ttqs(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do y = NI+1, NIA
                  Bux_AA(b, y) = Bux_AA(b, y) + Occ(t) * V_axby(t,y)
            end do
      end subroutine p2_AAIA_qstt_ttqs

      subroutine p2_AAIA_qstt_ttsq(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do y = NI+1, NIA
                  Bux_AA(y, b) = Bux_AA(y, b) + Occ(t) * V_axby(t,y)
            end do
      end subroutine p2_AAIA_qstt_ttsq

      subroutine p2_IVII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IVII_qstt_qstt

      subroutine p2_IVII_qstt_sqtt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_II(x, a) = Bux_II(x, a) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IVII_qstt_sqtt

      subroutine p2_IVAA_qstt_qstt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_AA(a, x) = Bux_AA(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IVAA_qstt_qstt

      subroutine p2_IVAA_qstt_sqtt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_AA(x, a) = Bux_AA(x, a) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_IVAA_qstt_sqtt

      subroutine p2_AAAA_qstt_qstt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_AA(a, x) = Bux_AA(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_AAAA_qstt_qstt

      subroutine p2_VAII_qstt_sqtt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_II(x, a) = Bux_II(x, a) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VAII_qstt_sqtt

      subroutine p2_VAII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VAII_qstt_qstt

      subroutine p2_VAAA_qstt_sqtt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_AA(x, a) = Bux_AA(x, a) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VAAA_qstt_sqtt

      subroutine p2_VAAA_qstt_qstt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NI+1, NIA
                  Bux_AA(a, x) = Bux_AA(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VAAA_qstt_qstt

      subroutine p2_VVII_qstt_qstt(Bux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_II(a, x) = Bux_II(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VVII_qstt_qstt

      subroutine p2_VVAA_qstt_qstt(Bux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Bux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t=b
            do x = NIA+1, NBasis
                  Bux_AA(a, x) = Bux_AA(a, x) + Occ(t) * V_axby(x,t)
            end do

      end subroutine p2_VVAA_qstt_qstt

      subroutine p2_IIII_qtst_qtst(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = 1, NI
                  Aux_II(a, b) = Aux_II(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_IIII_qtst_qtst

      subroutine p2_IAII_qtst_tsqt(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do x = NI+1, NIA
                  Aux_II(b, x) = Aux_II(b, x) + Occ(t) *  V_axby(x,t)
            end do

      end subroutine p2_IAII_qtst_tsqt

      subroutine p2_IAII_qtst_tqst(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do x = NI+1, NIA
                  Aux_II(x, b) = Aux_II(x, b) + Occ(t) *  V_axby(x,t)
            end do

      end subroutine p2_IAII_qtst_tqst

      subroutine p2_IVII_qtst_tsqt(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do x = NIA+1, NBasis
                  Aux_II(b, x) = Aux_II(b, x) + Occ(t) *  V_axby(x,t)
            end do

      end subroutine p2_IVII_qtst_tsqt

      subroutine p2_IVII_qtst_tqst(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            t = a
            do x = NIA+1, NBasis
                  Aux_II(x, b) = Aux_II(x, b) + Occ(t) *  V_axby(x,t)
            end do

      end subroutine p2_IVII_qtst_tqst

      subroutine p2_IAIA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_IAIA_qtst_qtst

      subroutine p2_IAIA_qtst_tqts(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            if (a==b)then
                  t = a
                  do y = NI+1, NIA
                        do x = NI+1, NIA
                              Aux_II(x, y) = Aux_II(x, y) + Occ(t) *  V_axby(x,y)
                        end do
                  end do
            end if

      end subroutine p2_IAIA_qtst_tqts

      subroutine p2_AAIA_qtst_stqt(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(b, a) = Aux_AA(b, a) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_AAIA_qtst_stqt

      subroutine p2_AAIA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_AAIA_qtst_qtst

      subroutine p2_VAIA_qtst_stqt(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(b, a) = Aux_AA(b, a) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_VAIA_qtst_stqt

      subroutine p2_VAIA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_VAIA_qtst_qtst

      subroutine p2_IVIA_qtst_tstq(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            if (a==b)then
                  t = a
                  do y = NI+1, NIA
                        do x = NIA+1, NBasis
                              Aux_II(y, x) = Aux_II(y, x) + Occ(t) *  V_axby(x,y)
                        end do
                  end do
            end if

      end subroutine p2_IVIA_qtst_tstq

      subroutine p2_IVIA_qtst_tqts(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            if (a==b)then
                  t = a
                  do y = NI+1, NIA
                        do x = NIA+1, NBasis
                              Aux_II(x, y) = Aux_II(x, y) + Occ(t) *  V_axby(x,y)
                        end do
                  end do
            end if

      end subroutine p2_IVIA_qtst_tqts

      subroutine p2_AAAA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_AAAA_qtst_qtst

      subroutine p2_VAAA_qtst_stqt(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(b, a) = Aux_AA(b, a) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_VAAA_qtst_stqt

      subroutine p2_VAAA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_VAAA_qtst_qtst

      subroutine p2_IVIV_qtst_tqts(Aux_II, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_II
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            if (a==b)then
                  t = a
                  do y = NIA+1, NBasis
                        do x = NIA+1, NBasis
                              Aux_II(x, y) = Aux_II(x, y) + Occ(t) *  V_axby(x,y)
                        end do
                  end do
            end if

      end subroutine p2_IVIV_qtst_tqts

      subroutine p2_VAVA_qtst_qtst(Aux_AA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
            double precision, dimension(:,:), intent(inout) :: Aux_AA
            integer, intent(in) :: a, b
            integer, intent(in) :: x0, x1, y0, y1
            integer, intent(in) ::NI, NIA, NBasis
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
            integer :: x, y, pq, rs, t


            do t = NI+1, NIA
                  Aux_AA(a, b) = Aux_AA(a, b) + Occ(t) *  V_axby(t,t)
            end do

      end subroutine p2_VAVA_qtst_qtst


end module THC_auto_p2
