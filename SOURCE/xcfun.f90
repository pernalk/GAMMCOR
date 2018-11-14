subroutine RhoKernel(Rho,SRKer,IFunSR,XMu,NGrid)
!
! computes short-range kernel for densities stored in a vector Rho
!
   use xcfun

   implicit none

   character(1000)      :: text
   integer              :: id, order, ilen, olen
   real(8), allocatable :: density_variables(:, :), derivatives(:, :)

   integer              :: NGrid, IFunSR, i
   real(8)              :: XMu

   real(8)              :: Rho,SRKer
   dimension Rho(NGrid),SRKer(NGrid)

!  print some info and copyright about the library
!  please always include this info in your code
!   call xcfun_splash(text)
!   print *, text(1:len_trim(text))

!  create a new functional
   id = xc_new_functional()
!  use sr-lda
   if(IFunSR>2) then
     call xc_set_param(id, XC_SLATERX, 1.0d0)
     call xc_set_param(id, XC_VWN5C, 1.0d0)
   else
     call xc_set_param(id, XC_LDAERFX, 1.0d0)
     call xc_set_param(id, XC_LDAERFC, 1.0d0)
     call xc_set_param(id, XC_RANGESEP_MU, XMu)
   endif
! use XC_VARS_N for closed-shell calculations
   call xc_set_mode(id, XC_VARS_N)
!   call xc_set_mode(id, XC_VARS_NS)

   order=2

   ilen = xc_input_length(id)
   olen = xc_output_length(id, order)
    
   allocate(density_variables(ilen, NGrid))
   allocate(derivatives(olen, NGrid))

   do i=1,NGrid
   density_variables(1, i) = Rho(i)
!   density_variables(2, i) = 0.d0
   enddo 

   call xc_eval(id, order, NGrid, density_variables, derivatives)

   do i=1,NGrid
!   SRKer(i)=derivatives(XC_D10, i)
   SRKer(i)=derivatives(XC_D2, i)
   enddo  

   deallocate(density_variables)
   deallocate(derivatives)

   call xc_free_functional(id) 

return
end

subroutine LYP(Rho,Sigma,Ene,NGrid)
!
! computes short-range kernel for densities stored in a vector Rho
!
   use xcfun

   implicit none

   character(1000)      :: text
   integer              :: id, order, ilen, olen
   real(8), allocatable :: density_variables(:, :), derivatives(:, :)

   integer              :: NGrid, i

   real(8)              :: Rho,Sigma,Ene
! ,VKS
   dimension Rho(NGrid),Sigma(NGrid),Ene(NGrid)
! ,VKS(NGrid)

!  create a new functional
   id = xc_new_functional()
!  use the following functionals
   call xc_set_param(id, XC_LYPC, 1.0d0)
! use XC_VARS_N for closed-shell calculations
   call xc_set_mode(id, XC_VARS_N)
!   call xc_set_mode(id, XC_VARS_NS)
! order=0 -> only xc energy on a grid
   order=0

   ilen = xc_input_length(id)
   olen = xc_output_length(id, order)

   allocate(density_variables(ilen, NGrid))
   allocate(derivatives(olen, NGrid))

   do i=1,NGrid
   density_variables(1, i) = Rho(i)
   density_variables(2, i) = Sigma(i)
   enddo

! evaluate what you need
   call xc_eval(id, order, NGrid, density_variables, derivatives)

   do i=1,NGrid
   Ene(i)=derivatives(1, i)
!   VKS(i)=derivatives(2, i) 
   enddo

   deallocate(density_variables)
   deallocate(derivatives)

   call xc_free_functional(id)

return
end subroutine LYP 

subroutine dfun_PBE(Rho,Sigma,Ene,vrhoc,vsigmacc,NGrid)
!
! computes stuff 
!
   use xcfun

   implicit none

   integer,intent(in)  :: NGrid
   double precision,intent(in)  :: Rho(NGrid),Sigma(NGrid)
   double precision,intent(out) :: vrhoc(NGrid),vsigmacc(NGrid),Ene(NGrid)

   character(1000)      :: text
   integer              :: i,id,order,ilen,olen
   double precision, allocatable :: density_variables(:, :), derivatives(:, :)

!  create a new functional
   id = xc_new_functional()
!  use the following functionals
   call xc_set_param(id, XC_PBEX, 1.0d0)
   call xc_set_param(id, XC_PBEC, 1.0d0)
! use XC_VARS_N for closed-shell calculations
   call xc_set_mode(id, XC_VARS_N)
!   call xc_set_mode(id, XC_VARS_NS)
! order=0 -> only xc energy on a grid
   order=1

   ilen = xc_input_length(id)
   olen = xc_output_length(id, order)
   print *, 'Length of output (how many derivatives):', olen

   allocate(density_variables(ilen, NGrid))
   allocate(derivatives(olen, NGrid))

   do i=1,NGrid
      density_variables(1, i) = Rho(i)
      density_variables(2, i) = Sigma(i)
   enddo

! evaluate what you need
   call xc_eval(id, order, NGrid, density_variables, derivatives)

   do i=1,NGrid
      Ene(i)=derivatives(XC_D00, i)
      vrhoc(i)=derivatives(XC_D10, i) 
      vsigmacc(i)=derivatives(XC_D01, i) 
   enddo

   deallocate(density_variables)
   deallocate(derivatives)

   call xc_free_functional(id)

end subroutine dfun_PBE

