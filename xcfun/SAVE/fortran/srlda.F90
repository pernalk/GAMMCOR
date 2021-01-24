program xc_example

!  this example contains all f90 interface routines
!  that are needed to "talk to" the xcfun library
!  and demonstrates how to use them

   use xcfun

   implicit none

   character(1000)      :: text
   integer              :: id, order, ilen, olen
   integer              :: i, k, ipoint, nr_points, block_length, max_block_length
   real(8)              :: derivative_nn_ab, derivative_ss_ab
   real(8)              :: derivative_nn_rs, derivative_ss_rs
   real(8), allocatable :: density_variables(:, :), derivatives(:, :)

!  print some info and copyright about the library
!  please always include this info in your code
!   call xcfun_splash(text)
!   print *, text(1:len_trim(text))

!  create a new functional
!  we need this for interacting with the library
   id = xc_new_functional()
!  use sr-lda
   call xc_set_param(id, XC_LDAERFX, 1.0d0)
   call xc_set_param(id, XC_LDAERFC, 1.0d0)
   call xc_set_param(id, XC_RANGESEP_MU, 0.3d0) 
   !  we switch to total/spin density mode (set the second variable to zero for closed-shell)
   call xc_set_mode(id, XC_VARS_NS)
    
   order=2

   ilen = xc_input_length(id)
   olen = xc_output_length(id, order)

   allocate(density_variables(ilen, 1))
   allocate(derivatives(olen, 1))

   density_variables(1, 1) = 9.0136283818039951E-006
   density_variables(2, 1) = 0.0d0

   call xc_eval(id, order, 1, density_variables, derivatives)
   derivative_nn_rs = derivatives(XC_D20, 1)
   derivative_ss_rs = derivatives(XC_D02, 1)
   write(*,*)'potential',derivatives(XC_D10, 1)


   deallocate(density_variables)
   deallocate(derivatives)

end program
