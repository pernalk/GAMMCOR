module mp

      use acpp_types
      use types
      use math_constants
      use real_linalg !from gammcor integrals
      use THC_Gammcor
      use sort  !from gammcor integrals         

      implicit none
      
      contains

            subroutine mp2_driver(BasisSet, THCData, AuxData, CAONO, Flags)
                  Character(*) :: BasisSet
                  type(TACppData), intent(in) :: AuxData
                  double precision, dimension(:,:), intent(in) :: CAONO
                  type(FlagsData), intent(in) :: Flags
                  type(TTHCData), intent(inout) :: THCData


                  double precision, dimension(:,:),  allocatable  :: XMuMat
                  double precision, dimension(:,:),  allocatable  :: URe, CAONOT
                  double precision :: AvMu,ECorrMD, emp2
                  double precision, allocatable :: Rkab(:,:,:), Rkcd(:,:,:), Rkef(:,:,:)
                  double precision :: ETot, etot0, val, val_this, thisfr, thistr, thistr1
                  double precision, dimension(:,:,:), allocatable :: RAV, RVI
                  double precision, dimension(:,:,:), allocatable :: RVI_FR, RVI_FR1, RVI_FR2
                  double precision, dimension(:,:,:), allocatable :: RVI_LR, RVI_LR1, RVI_LR2

                  double precision, dimension(:,:,:), allocatable :: RII_FR, RII_FR1, RII_FR2
                  double precision, dimension(:,:,:), allocatable :: RII_LR, RII_LR1, RII_LR2  

                  double precision, dimension(:,:), allocatable :: this_VV
                  double precision, dimension(:,:), allocatable :: this, this_FR, this_FR1, this_FR2, this_FR3, this_FR4
                  double precision, dimension(:,:), allocatable :: this_SR, this_LR1, this_LR2, this_LR3, this_LR4, this_LR, work

                  double precision :: numerator, denominator, aux
                  integer :: x, y, a, b, i, j, a_idx, b_idx, t, u, nnn
                  double precision, dimension(:), allocatable :: occ_frozen


                  integer :: NCholErf1, iunit

                  integer i0, licznik




                  
                  associate(Zgk=>THCData%Zgk, Xga=>THCData%Xga,  TXga=>THCData%TXga, ExternalOrdering=>THCData%ExternalOrdering, NChol=>THCData%NChol, NTHC=>THCData%NTHC, &
                        XgaErf=>THCData%XgaErf, TXgaErf=>THCData%TXgaErf, ZgkErf=>THCData%ZgkErf, NCholErf=>THCData%NCholErf, NTHCErf=>THCData%NTHCErf, NBasis=>AuxData%NBasis, &
                        NI=>AuxData%NI, NA=>AuxData%NA, NIA=>AuxData%NIA, NV=>AuxData%NV)




                    allocate(occ_frozen(NBasis))
                    occ_frozen = AuxData%occ
                    

!                    allocate(RVI(NChol, NV, NI))
!                    allocate(RII(Nchol, NI, NI))

                    allocate(RVI_FR(NChol, NV, NI))
                    allocate(RII_FR(Nchol, NI, NI))
                    allocate(RVI_LR(NCholErf, NV, NI))
                    allocate(RII_LR(NcholErf, NI, NI))

                    
                    allocate(RVI_FR1(NChol, NV, NI))
                    allocate(RVI_FR2(NChol, NV, NI))


                    allocate(RII_FR1(Nchol, NI, NI))
                    allocate(RII_FR2(Nchol, NI, NI))


                    allocate(RVI_LR1(NCholErf, NV, NI))
                    allocate(RVI_LR2(NCholErf, NV, NI))


                    allocate(RII_LR1(NcholErf, NI, NI))
                    allocate(RII_LR2(NcholErf, NI, NI))

                    
                    !allocate(RAV(NChol, NV, NA))

                    ! For frozen core                                                                                                                                                                                                         
                    if (AuxData%NCoreOrb > 0)then
                          i0 = 1+ AuxData%NCoreOrb
                          do i = 1, i0-1
                                occ_frozen(i) = zero
                          end do
                    else
                          i0 = 1
                    end if



!                    print*, 'THCData%fij', THCData%fij
                    
!                    print*, 'THCData%fvw', THCData%fvw

                     
                    if (Flags%JobType == JOB_TYPE_SRMP2)then
                          allocate(URe(NBasis, NBasis))
                          allocate(CAONOT(NBasis, NBasis))
                           
                          URe = zero
                          do i = 1, NBasis
                                URe(i, i) = one
                          end do
                          CAONOT = transpose(CAONO)
                          allocate(XMuMat(nbasis, nbasis))
                           
                          !Call LOC_MU_CBS_CHOL(XMuMat,ECorrMD,AvMU,URe, CAONOT,Occ_frozen,BasisSet,NBasis)

                          call real_abt(THCData%TXgaErf, THCData%XgaErf, XMuMAT)
                          call real_abt(THCData%TXga, THCData%Xga, XMuMAT)                                 

                     end if

                    
                    
                     call thc_gammcor_Rkab_2(RVI_FR, Xga(:,NIA+1:NV), Xga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
                     call thc_gammcor_Rkab_2(RII_FR, Xga(:,1:NI), Xga(:,1:NI), Zgk, NI, NI, NChol, NTHC)

                     if (Flags%JobType == JOB_TYPE_SRMP2)then

!                           allocate(RRVI_FR(NChol, NV, NI))
!                           allocate(RVI_FRt(NChol, NV, NI))

                           
                           call thc_gammcor_Rkab_2(RVI_FR1, TXga(:,NIA+1:NV), Xga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
                           call thc_gammcor_Rkab_2(RVI_FR2, Xga(:,NIA+1:NV), TXga(:,1:NI), Zgk, NV, NI, NChol, NTHC)
                           
                           call thc_gammcor_Rkab_2(RII_FR1, TXga(:,1:NI), Xga(:,1:NI), Zgk, NI, NI, NChol, NTHC)
                           call thc_gammcor_Rkab_2(RII_FR2, Xga(:,1:NI), TXga(:,1:NI), Zgk, NI, NI, NChol, NTHC)
                           
                           call thc_gammcor_Rkab_2(RVI_LR, XgaErf(:,NIA+1:NV), XgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)
                           
                           call thc_gammcor_Rkab_2(RVI_LR1, TXgaErf(:,NIA+1:NV), XgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)
                           call thc_gammcor_Rkab_2(RVI_LR2, XgaErf(:,NIA+1:NV), TXgaErf(:,1:NI), ZgkErf, NV, NI, NCholErf, NTHCErf)

                           call thc_gammcor_Rkab_2(RII_LR1, TXgaErf(:,1:NI), XgaErf(:,1:NI), ZgkErf, NI, NI, NCholErf, NTHCErf)
                           call thc_gammcor_Rkab_2(RII_LR2, XgaErf(:,1:NI), TXgaErf(:,1:NI), ZgkErf, NI, NI, NCholErf, NTHCErf)


                     end if

                     allocate(this_VV(NV, NV))
                     allocate(this_SR(NV, NV))
                     allocate(this_FR(NV, NV))

                    emp2 = zero
                    
                    if (Flags%JobType == JOB_TYPE_SRMP2)then
                          !if (Flags%IDBBSC == 2)then
                          licznik = 0
                                do i = i0, NI
                                      do j = i0, NI
                                            
                                            call real_aTb_x(this_FR, NV, RVI_FR(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, One)
                                            
                                            call real_aTb_x(this_SR, NV, RVI_FR1(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, FRAC14)
                                            call real_aTb_x(this_SR, NV, RVI_FR2(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, FRAC14, One)
                                            call real_aTb_x(this_SR, NV, RVI_FR(:,:,i), NChol, RVI_FR1(:,:,j), NChol, NV, NV, NChol, FRAC14, One)
                                            call real_aTb_x(this_SR, NV, RVI_FR(:,:,i), NChol, RVI_FR2(:,:,j), NChol, NV, NV, NChol, FRAC14, One)
                                            
                                            call real_aTb_x(this_SR, NV, RVI_LR1(:,:,i), NCholErf, RVI_LR(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                                            call real_aTb_x(this_SR, NV, RVI_LR2(:,:,i), NCholErf, RVI_LR(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                                            call real_aTb_x(this_SR, NV, RVI_LR(:,:,i), NCholErf, RVI_LR1(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                                            call real_aTb_x(this_SR, NV, RVI_LR(:,:,i), NCholErf, RVI_LR2(:,:,j), NCholErf, NV, NV, NCholErf, -FRAC14, One)
                                            

                                            this_VV = this_FR + this_SR

                                            do a = 1, NV
                                                  do b = 1, NV
                                                        numerator  =  this_VV(a,b) *(two * this_VV(a,b) - this_VV(b,a))
                                                        denominator = THCData%fij(i)+THCData%fij(j)-THCData%fvw(a)-THCData%fvw(b)
                                                        emp2  = emp2 + numerator/denominator
                                                  end do
                                            end do
                                      end do
                                end do
                                print*, 'SRMp2', emp2
                          else

                          do i = i0, NI
                                do j = i0, NI
                                      call real_aTb_x(this_VV, NV, RVI_FR(:,:,i), NChol, RVI_FR(:,:,j), NChol, NV, NV, NChol, One)
                                      do a = 1, NV
                                            do b = 1, NV
                                                  numerator  =  ( two * this_VV(a,b) - this_VV(b,a))
                                                  denominator = THCData%fij(i)+THCData%fij(j)-THCData%fvw(a)-THCData%fvw(b)
                                                  aux = numerator/denominator
                                                  emp2  = emp2 + this_VV(a,b) * numerator/denominator
                                            end do
                                      end do
                                end do
                          end do
                          print*, 'emp2-kod1', emp2
                    end if

            stop
          end associate
      end subroutine mp2_driver

      
end module mp
