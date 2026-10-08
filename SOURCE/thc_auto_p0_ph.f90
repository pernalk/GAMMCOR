module thc_auto_p0_ph
      use math_constants
      implicit none

contains

subroutine d2a_IA_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d2a_IA_aoao_pqrs

subroutine d2a_IA_aovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d2a_IA_aovo_pqrs

subroutine d2a_IA_aavo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = y0, min(y1, b-1)

        pq = pos(x, a)
        rs = pos(b, y)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d2a_IA_aavo_pqrs

subroutine t4a_II_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_II_aoao_pqrs

subroutine t4a_AI_aoaa_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_AI_aoaa_pqrs

subroutine t4a_II_aovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_II_aovo_pqrs

subroutine t4a_IA_aavo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = y0, min(y1, b-1)

        pq = pos(x, a)
        rs = pos(b, y)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_IA_aavo_pqrs

subroutine t4a_II_vovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_II_vovo_pqrs

subroutine t4a_AI_aova_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_AI_aova_pqrs

subroutine t4a_AI_vova_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_AI_vova_pqrs

subroutine d4a_AI_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d4a_AI_aoao_pqrs

subroutine d4a_AI_aoaa_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d4a_AI_aoaa_pqrs

subroutine d1b_IA_aoao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(a) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d1b_IA_aoao_qrps

subroutine d1b_IA_aovo_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(x, b)
        rs = pos(y, a)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d1b_IA_aovo_qrps

subroutine d1b_IA_aavo_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = max(y0, x+1), y1

                    
        pq = pos(a, b)
        rs = pos(y, x)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d1b_IA_aavo_qrps

subroutine t3b_II_aoao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t3b_II_aoao_prqs

subroutine t3b_AI_aoaa_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = y0, min(y1, a-1)

        pq = pos(a, y)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(y) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t3b_AI_aoaa_prqs

subroutine t3b_II_aovo_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t3b_II_aovo_prqs

subroutine t3b_IA_aavo_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do y = y0, y1
                  do x = max(x0, y+1), x1

                    
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t3b_IA_aavo_prqs

subroutine t3b_AI_aova_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = y0, min(y1, a-1)

        pq = pos(a, y)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(y) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t3b_AI_aova_prqs

subroutine t3b_II_vovo_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t3b_II_vovo_prqs

subroutine t3b_AI_vova_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = y0, min(y1, a-1)

        pq = pos(a, y)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(y) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t3b_AI_vova_prqs

subroutine d3b_AI_aoao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(y) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d3b_AI_aoao_qrps

subroutine d3b_AI_aoaa_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = -One * Occ(y) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d3b_AI_aoaa_qrps

subroutine t5a_IA_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t5a_IA_aoao_pqrs

subroutine t5a_IA_aovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t5a_IA_aovo_pqrs

subroutine t5a_IA_aavo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = y0, min(y1, b-1)

        pq = pos(x, a)
        rs = pos(b, y)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t5a_IA_aavo_pqrs

subroutine t5b_IA_aoao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t5b_IA_aoao_prqs

subroutine t5b_IA_aovo_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t5b_IA_aovo_prqs

subroutine t5b_IA_aavo_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do y = y0, y1
                  do x = max(x0, y+1), x1

                    
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t5b_IA_aavo_prqs

subroutine t6a_AI_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t6a_AI_aoao_pqrs

subroutine t6a_AI_aoaa_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t6a_AI_aoaa_pqrs

subroutine d6a_II_aoao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_II_aoao_pqrs

subroutine d6a_AI_aoaa_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_AI_aoaa_pqrs

subroutine d6a_II_aovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_II_aovo_pqrs

subroutine d6a_IA_aavo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = y0, min(y1, b-1)

        pq = pos(x, a)
        rs = pos(b, y)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_IA_aavo_pqrs

subroutine d6a_II_vovo_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, a+1), x1
            do y = max(y0, b+1), y1

        pq = pos(x, a)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_II_vovo_pqrs

subroutine d6a_AI_aova_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_AI_aova_pqrs

subroutine d6a_AI_vova_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(a, x)
        rs = pos(y, b)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_AI_vova_pqrs

subroutine t6b_AI_aoao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t6b_AI_aoao_prqs

subroutine t6b_AI_aoaa_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = y0, min(y1, a-1)

        pq = pos(a, y)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t6b_AI_aoaa_prqs

subroutine d6b_II_aoao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6b_II_aoao_qrps

subroutine d6b_AI_aoaa_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6b_AI_aoaa_qrps

subroutine d6b_II_aovo_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(x, b)
        rs = pos(y, a)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(a)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6b_II_aovo_qrps

subroutine d6b_IA_aavo_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = max(y0, x+1), y1

                    
        pq = pos(a, b)
        rs = pos(y, x)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d6b_IA_aavo_qrps

subroutine d6b_II_vovo_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(y, a)
        rs = pos(x, b)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6b_II_vovo_qrps

subroutine d6b_AI_aova_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = max(y0, a+1), y1

        pq = pos(x, b)
        rs = pos(y, a)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(a)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6b_AI_aova_qrps

subroutine d6b_AI_vova_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(y, x)
        rs = pos(a, b)
        if (pq>0.and.rs>0)then
        val = One * Occ(x) * Occ(b)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d6b_AI_vova_qrps

subroutine d1b_IA_aaao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = max(y0, x+1), y1

                    
        pq = pos(a, b)
        rs = pos(y, x)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d1b_IA_aaao_qrps


subroutine d2a_IA_aaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = Two * Occ(b) * Occ(a)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d2a_IA_aaao_pqrs


subroutine d6a_IA_aaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_IA_aaao_pqrs


subroutine d6b_IA_aaao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = max(y0, x+1), y1

                    
        pq = pos(a, b)
        rs = pos(y, x)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d6b_IA_aaao_qrps


subroutine t3b_IA_aaao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do y = y0, y1
                  do x = max(x0, y+1), x1

                    
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t3b_IA_aaao_prqs


subroutine t4a_IA_aaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = Two * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_IA_aaao_pqrs


subroutine t5a_IA_aaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(b) * Occ(a)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t5a_IA_aaao_pqrs


subroutine t5b_IA_aaao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (b<a)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do y = y0, y1
                  do x = max(x0, y+1), x1

                    
        pq = pos(a, b)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine t5b_IA_aaao_prqs


subroutine d6a_IA_vaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = -Two * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine d6a_IA_vaao_pqrs


subroutine d6b_IA_vaao_qrps(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        if (a<b)then
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, x1
                  do y = y0, y1

                
        pq = pos(b, a)
        rs = pos(x, y)
        if (pq>0.and.rs>0)then
        val = One * Occ(a) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) - val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        end if
end subroutine d6b_IA_vaao_qrps


subroutine t3b_IA_vaao_prqs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = max(x0, b+1), x1
            do y = y0, min(y1, a-1)

        pq = pos(x, b)
        rs = pos(a, y)
        if (pq>0.and.rs>0)then
        val = -One * Occ(b) * Occ(y)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t3b_IA_vaao_prqs


subroutine t4a_IA_vaao_pqrs(APlus, AMinus, a, b, x0, x1, y0, y1, pos, Occ, V_axby)
double precision, dimension(:,:), intent(inout) :: APlus, AMinus
integer, intent(in) :: a, b
integer, intent(in) :: x0, x1, y0, y1
integer, dimension(:,:), intent(in) :: pos
double precision, dimension(:), intent(in) :: Occ
double precision, dimension(x0:x1, y0:y1), intent(in) :: V_axby
integer :: x, y, pq, rs

double precision :: val
    
        
        !!$omp parallel default(shared)
        !!$omp private(x, y, pq, rs, val)
        !!$omp do collapse(2)
                    do x = x0, min(x1, a-1)
            do y = max(y0, b+1), y1

        pq = pos(y, b)
        rs = pos(a, x)
        if (pq>0.and.rs>0)then
        val = Two * Occ(b) * Occ(x)* V_axby(x,y)        
        APlus(rs, pq) = APlus(rs, pq) + val

        AMinus(rs, pq) = AMinus(rs, pq) + val
        end if
        end do
        end do
        !!$omp end do
        !!$omp end parallel
        
end subroutine t4a_IA_vaao_pqrs
end module thc_auto_p0_ph
