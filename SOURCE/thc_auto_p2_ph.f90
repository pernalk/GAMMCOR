module thc_auto_p2_ph
      use math_constants
      implicit none

contains

subroutine p2_IIII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = 1, NI
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IIII_qstt_qstt

subroutine p2_AAII_qstt_ttqs(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do y = 1, NI
        BuxAA(b, y) = BuxAA(b, y) + Occ(t) * V_axby(t,y)
        end do
end subroutine p2_AAII_qstt_ttqs

subroutine p2_AAII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_AAII_qstt_qstt

subroutine p2_IAII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IAII_qstt_qstt

subroutine p2_IAII_qstt_sqtt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxII(x, a) = BuxII(x, a) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IAII_qstt_sqtt

subroutine p2_AAIA_qstt_ttqs(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do y = NI+1, NIA
        BuxAA(b, y) = BuxAA(b, y) + Occ(t) * V_axby(t,y)
        end do
end subroutine p2_AAIA_qstt_ttqs

subroutine p2_AAIA_qstt_ttsq(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do y = NI+1, NIA
        BuxAA(y, b) = BuxAA(y, b) + Occ(t) * V_axby(t,y)
        end do
end subroutine p2_AAIA_qstt_ttsq

subroutine p2_IVII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IVII_qstt_qstt

subroutine p2_IVII_qstt_sqtt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxII(x, a) = BuxII(x, a) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IVII_qstt_sqtt

subroutine p2_IVAA_qstt_qstt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxAA(a, x) = BuxAA(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IVAA_qstt_qstt

subroutine p2_IVAA_qstt_sqtt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxAA(x, a) = BuxAA(x, a) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_IVAA_qstt_sqtt

subroutine p2_AAAA_qstt_qstt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxAA(a, x) = BuxAA(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_AAAA_qstt_qstt

subroutine p2_VAII_qstt_sqtt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxII(x, a) = BuxII(x, a) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VAII_qstt_sqtt

subroutine p2_VAII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VAII_qstt_qstt

subroutine p2_VAAA_qstt_sqtt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxAA(x, a) = BuxAA(x, a) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VAAA_qstt_sqtt

subroutine p2_VAAA_qstt_qstt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NI+1, NIA
        BuxAA(a, x) = BuxAA(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VAAA_qstt_qstt

subroutine p2_VVII_qstt_qstt(BuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxII(a, x) = BuxII(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VVII_qstt_qstt

subroutine p2_VVAA_qstt_qstt(BuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: BuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t=b
        do x = NIA+1, NBasis
        BuxAA(a, x) = BuxAA(a, x) + Occ(t) * V_axby(x,t)
        end do
        
end subroutine p2_VVAA_qstt_qstt


subroutine p2_IIII_qtst_qtst(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = 1, NI
        AuxII(a, b) = AuxII(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_IIII_qtst_qtst

subroutine p2_IAII_qtst_tsqt(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do x = NI+1, NIA
        AuxII(b, x) = AuxII(b, x) + Occ(t) *  V_axby(x,t)
        end do
        
end subroutine p2_IAII_qtst_tsqt

subroutine p2_IAII_qtst_tqst(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do x = NI+1, NIA
        AuxII(x, b) = AuxII(x, b) + Occ(t) *  V_axby(x,t)
        end do
        
end subroutine p2_IAII_qtst_tqst

subroutine p2_IVII_qtst_tsqt(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do x = NIA+1, NBasis
        AuxII(b, x) = AuxII(b, x) + Occ(t) *  V_axby(x,t)
        end do
        
end subroutine p2_IVII_qtst_tsqt

subroutine p2_IVII_qtst_tqst(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        t = a
        do x = NIA+1, NBasis
        AuxII(x, b) = AuxII(x, b) + Occ(t) *  V_axby(x,t)
        end do
        
end subroutine p2_IVII_qtst_tqst

subroutine p2_IAIA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_IAIA_qtst_qtst

subroutine p2_IAIA_qtst_tqts(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
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
        AuxII(x, y) = AuxII(x, y) + Occ(t) *  V_axby(x,y)
        end do
        end do
        end if
        
end subroutine p2_IAIA_qtst_tqts

subroutine p2_AAIA_qtst_stqt(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(b, a) = AuxAA(b, a) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_AAIA_qtst_stqt

subroutine p2_AAIA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_AAIA_qtst_qtst

subroutine p2_VAIA_qtst_stqt(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(b, a) = AuxAA(b, a) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_VAIA_qtst_stqt

subroutine p2_VAIA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_VAIA_qtst_qtst

subroutine p2_IVIA_qtst_tstq(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
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
        AuxII(y, x) = AuxII(y, x) + Occ(t) *  V_axby(x,y)
        end do
        end do
        end if
        
end subroutine p2_IVIA_qtst_tstq

subroutine p2_IVIA_qtst_tqts(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
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
        AuxII(x, y) = AuxII(x, y) + Occ(t) *  V_axby(x,y)
        end do
        end do
        end if
        
end subroutine p2_IVIA_qtst_tqts

subroutine p2_AAAA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_AAAA_qtst_qtst

subroutine p2_VAAA_qtst_stqt(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(b, a) = AuxAA(b, a) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_VAAA_qtst_stqt

subroutine p2_VAAA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_VAAA_qtst_qtst

subroutine p2_IVIV_qtst_tqts(AuxII, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxII
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
        AuxII(x, y) = AuxII(x, y) + Occ(t) *  V_axby(x,y)
        end do
        end do
        end if
        
end subroutine p2_IVIV_qtst_tqts

subroutine p2_VAVA_qtst_qtst(AuxAA, a, b, x0, x1, y0, y1, NI, NIA, NBasis,  Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: AuxAA
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, intent(in) ::NI, NIA, NBasis
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs, t

    
        do t = NI+1, NIA
        AuxAA(a, b) = AuxAA(a, b) + Occ(t) *  V_axby(t,t)
        end do
        
end subroutine p2_VAVA_qtst_qtst
end module thc_auto_p2_ph
