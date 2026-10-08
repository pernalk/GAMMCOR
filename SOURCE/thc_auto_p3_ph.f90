module thc_auto_p3_ph
      use math_constants
      use real_linalg
      implicit none

contains

subroutine T5_aova(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! T5 blok aova: - sum_tu Gamma_tuqr (su|pt), klasa VAIA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, c)
                  rs = pos(d, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T5_aova

subroutine D2_aova(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! D2 blok aova: + sum_tu Gamma_trqu (su|pt), klasa VAIA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, d)
                  rs = pos(c, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D2_aova

subroutine D1_aova(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! D1 blok aova: + sum_tu Gamma_truq (sp|tu), klasa IVAA, slots (a,x,b,y) = sptu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(e2, d)
                  rs = pos(c, e1)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D1_aova

subroutine T4_vava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! T4 blok vava: + sum_tu Gamma_tsqu (ru|pt), klasa VAVA, slots (a,x,b,y) = rupt, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(b, c)
                  rs = pos(a, d)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T4_vava

subroutine T3_vava(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! T3 blok vava: + sum_tu Gamma_tsuq (pr|tu), klasa VVAA, slots (a,x,b,y) = prtu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(e1, d)
                  rs = pos(e2, c)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T3_vava

subroutine D6_vava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! D6 blok vava: - sum_tu Gamma_tuqs (ur|pt), klasa VAVA, slots (a,x,b,y) = rupt, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(b, d)
                  rs = pos(a, c)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D6_vava

subroutine T2_aoao(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! T2 blok aoao: + sum_tu Gamma_trpu (su|qt), klasa IAIA, slots (a,x,b,y) = suqt, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(c, b)
                  rs = pos(d, a)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T2_aoao

subroutine T1_aoao(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! T1 blok aoao: + sum_tu Gamma_trup (sq|tu), klasa IIAA, slots (a,x,b,y) = sqtu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(d, e2)
                  rs = pos(c, e1)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T1_aoao

subroutine D5_aoao(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! D5 blok aoao: - sum_tu Gamma_tupr (su|qt), klasa IAIA, slots (a,x,b,y) = suqt, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(d, b)
                  rs = pos(c, a)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D5_aoao

subroutine T4_aava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! T4 blok aava: + sum_tu Gamma_tsqu (ru|pt), klasa VAAA, slots (a,x,b,y) = ptru, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, d)
                  rs = pos(b, c)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T4_aava

subroutine T3_aava(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! T3 blok aava: + sum_tu Gamma_tsuq (pr|tu), klasa VAAA, slots (a,x,b,y) = prtu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(e1, d)
                  rs = pos(e2, c)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T3_aava

subroutine T5_aava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! T5 blok aava: - sum_tu Gamma_tuqr (su|pt), klasa VAAA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, c)
                  rs = pos(d, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T5_aava

subroutine D2_aava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! D2 blok aava: + sum_tu Gamma_trqu (su|pt), klasa VAAA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, d)
                  rs = pos(c, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D2_aava

subroutine D1_aava(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! D1 blok aava: + sum_tu Gamma_truq (sp|tu), klasa AVAA, slots (a,x,b,y) = sptu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(e2, d)
                  rs = pos(c, e1)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D1_aava

subroutine D6_aava(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! D6 blok aava: - sum_tu Gamma_tuqs (ur|pt), klasa VAAA, slots (a,x,b,y) = ptru, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, c)
                  rs = pos(b, d)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D6_aava

subroutine T2_aoaa(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! T2 blok aoaa: + sum_tu Gamma_trpu (su|qt), klasa AAIA, slots (a,x,b,y) = qtsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(d, a)
                  rs = pos(c, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T2_aoaa

subroutine T1_aoaa(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! T1 blok aoaa: + sum_tu Gamma_trup (sq|tu), klasa IAAA, slots (a,x,b,y) = sqtu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(d, e2)
                  rs = pos(c, e1)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T1_aoaa

subroutine T5_aoaa(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! T5 blok aoaa: - sum_tu Gamma_tuqr (su|pt), klasa AAIA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, c)
                  rs = pos(d, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine T5_aoaa

subroutine D2_aoaa(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_13ph, x0, x1, y0, y1, V_axby, NA)
! D2 blok aoaa: + sum_tu Gamma_trqu (su|pt), klasa AAIA, slots (a,x,b,y) = ptsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_13ph
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(a, d)
                  rs = pos(c, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_13ph(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D2_aoaa

subroutine D1_aoaa(e1, e2, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2_12, x0, x1, y0, y1, V_axby, NA)
! D1 blok aoaa: + sum_tu Gamma_truq (sp|tu), klasa IAAA, slots (a,x,b,y) = sptu, together
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: e1, e2, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2_12
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(e2, d)
                  rs = pos(c, e1)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2_12(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  + val
                        AMinus(rs, pq) = AMinus(rs, pq) - val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D1_aoaa

subroutine D5_aoaa(a, b, c0, c1, d0, d1, APlus, AMinus, map, pos, rdm2, x0, x1, y0, y1, V_axby, NA)
! D5 blok aoaa: - sum_tu Gamma_tupr (su|qt), klasa AAIA, slots (a,x,b,y) = qtsu, split
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b, NA
integer, intent(in) :: c0, c1, d0, d1
integer, dimension(:), intent(in) :: map
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:,:,:,:), intent(in) :: rdm2
integer, intent(in) :: x0, x1, y0, y1
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: c, d, pq, rs
double precision :: val
    
      !$omp parallel do default(shared) private(c, d, pq, rs, val) collapse(2)
      do c = c0, c1
            do d = d0, d1
                  pq = pos(c, a)
                  rs = pos(d, b)
                  if (pq>0.and.rs>0)then
                        call real_vw_x(val, rdm2(:,:,map(c), map(d)), V_axby, NA**2)
                        APlus(rs, pq)  = APlus(rs, pq)  - val
                        AMinus(rs, pq) = AMinus(rs, pq) + val
                  end if
            end do
      end do
      !$omp end parallel do

end subroutine D5_aoaa
end module thc_auto_p3_ph
