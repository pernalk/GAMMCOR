module ppAC0_fofo

      use iso_fortran_env
      use real_linalg !from gammcor integrals
      use sort  !from gammcor integrals      
      use types
      use lin
      use clock

      use tran
      use math_constants
      use quadratures
      use acpp_types
      use ppfunctions
      use ppAC0
      use ppAC0_THC

      implicit none

      integer, parameter,private ::  IOO=1, IVV=2, IAA=3, IVA=4, IAO=5
      integer, parameter,private :: IOOT=6, IVVT=7, IAAT=8, IVAT=9, IAOT = 10
      integer, parameter, private :: VVOO=1, VVAO=2, VVAA=3, VAOO=4
      integer, parameter, private :: VAAO=5, VAAA=6, AAOO=7, AAAO=8
      integer, parameter, private :: VAAO2 = 9

      integer, dimension(:), allocatable, private  :: IGem(:)




contains

      subroutine test_pluszowy(UMOAO, Flags, Sys, Occ, XOne, ninte1,  enuc, NI, NA,IndAux)

            use Cholesky_Gammcor
            use THC_Gammcor
            use OneElectronInts_Gammcor
            use basis_sets
            use sys_definitions
            use gammcor_integrals
            double precision, dimension(:,:), intent(in) :: UMOAO
            double precision, dimension(:), intent(in) :: Occ, XOne
            character(:), allocatable :: basis_path, xyz_path, binPath
            double precision, intent(in) :: enuc
            integer, intent(in) ::IndAux(:)
            integer, intent(in) ::NA, NI
            type(TAOBASIS) :: AObasis
            type(TSystem) :: System
            type(FlagsData), intent(in) :: Flags
            type(SystemBlock), intent(in) :: Sys
            integer, intent(in) :: ninte1
            integer, parameter :: units  = SYS_UNITS_ANGSTROM
            integer, parameter :: ExternalOrdering = ORBITAL_ORDERING_DALTON
            logical, parameter :: sort_ang_mom = .true.
            double precision, dimension(:,:), allocatable :: H0_extao
            integer :: i, j, nao
            integer :: NCholesky, NGridTHC
            double precision, allocatable :: Xgp(:,:), Zgk(:,:)
            double precision, allocatable :: Xga(:,:)
            double precision, allocatable ::  UAux(:,:)
            double precision, allocatable :: Rkab(:,:,:), Rkcd(:,:,:)
            double precision, allocatable :: hno(:)
            double precision :: this, ETot
            integer, dimension(:), allocatable :: map
            double precision, dimension(:), allocatable :: R00
            integer :: p, q, r, s, k, l, a, b, ab, ij, ii, nn
            integer, external :: NAddrRDM
            double precision, external :: FRDM2


            integer :: nbasis, io_status

            Character*60 FName,Aux1

            ! nbasis = 14

            ! allocate(XOne(nbasis, nbasis))

            ! binPath = "/home/aleksandra.tucholska/test/test1_ola/xone.bin"

            ! open(unit=10, file='xone.bin', form='unformatted', access='stream', status='old', iostat=io_status)

            ! if (io_status /= 0) then
            !    print *, "Error opening file. Iostat = ", io_status
            !    stop
            ! end if

            ! read(10) XOne

            ! close(10)

            ! print *, "XOne matrix:"
            ! do i = 1, nbasis
            !    do j = 1, nbasis
            !       write(*, '(F12.6)', advance='no') XOne(i,j)
            !    end do
            ! end do

            ! stop


            ! MOLPRO, ORCA

            ! do i =1, ninte1
            !       if(abs(XOne(i)).gt.1.d-5)then
            !             write(*, '(I8, F20.15)') i, XOne(i)
            !       end if
            ! end do
            ! stop
            ! print*, 'asdf'
            basis_path = Flags%BasisSetPath //Flags%BasisSet
            nbasis = Sys%NBasis
            print*, basis_path
            print*, nbasis
!            xyz_path = "/home/aleksandra.tucholska/test/test1_ola/ne/input.inp"

            xyz_path = "./input.inp"
            
            call auto2e_init()

            call sys_Read_XYZ(System, xyz_path, units)
            print*, System%NAtoms
            print*, System%NElectrons
            !stop
            call basis_newAObasis(AObasis, System, basis_path, .True., sort_ang_mom)

            nao = AObasis%NAOSpher

            !----------THC---------------
            print*, 'a1', nao


            !CHOL_ACCU_DEFAULT   = 1

            !CHOL_ACCU_TIGHT     = 2

            !CHOL_ACCU_LUDICROUS = 3

            
!            call thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, CHOL_ACCU_DEFAULT)
!            call thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, CHOL_ACCU_TIGHT)
            call thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, CHOL_ACCU_LUDICROUS)
            print*, 'a2'
            NGridTHC=size(Xgp,dim=1)
            NCholesky=size(Zgk,dim=2)
            print*, 'ngrid', 'nchol', ngridthc, ncholesky
            print*, 'a3'
            allocate(Xga(NGridTHC,NBasis))
            print*, 'xga-g', size(Xga, dim=1)
            print*, 'xga-2', size(Xga, dim=2)
            print*, 'xgp-g', size(Xgp, dim=1)
            print*,	'xgp-p', size(Xgp, dim=2)
            print*,	size(UMOAO, dim=1)
            print*,	size(UMOAO, dim=2)

            allocate(UAux(nbasis, nbasis))
            UAux  = transpose(UMOAO)

            Call thc_gammcor_Xga(Xga, Xgp, UAux,&
                  AOBasis, ExternalOrdering)
            print*, 'a4'
            allocate(Rkab(NCholesky,nao, nao))
            allocate(Rkcd(NCholesky, nao, nao))
            Call thc_gammcor_Rkab_2(Rkab, Xga, Xga, Zgk, NBasis, NBasis,&
                  NCholesky, NGridTHC)
            print*, 'a5'
            Call thc_gammcor_Rkab_2(Rkcd, Xga, Xga, Zgk, NBasis, NBasis,&
                  NCholesky, NGridTHC)


            ! calc calki 1-el



            allocate(hno(ninte1))
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                        HNO(ij) = XOne(ab)
                  end do
            end do

            ETot = zero
            do i = 1, NBasis
                  ii = (i*(i+1))/2
                  ETot = ETot + two* Occ(i) * HNO(ii)
                  print*, Occ(i) , HNO(ii)
            end do

            print*, 's1'
            Nn = NA**2*(NA**2+1)/2
            allocate(R00(nn))
            print*, 'nn', nn
            print*, 'etot1', ETot

            R00 = Zero
            call read_2rdm("rdm2.dat", R00, NA)

            allocate(map(nbasis))

            k = 1
            map = 0
            do i = 1, Nbasis
                  if(IndAux(i)==1)then
                        Map(i) = k
                        k = k+1
                  end if
                  if(IndAux(i)==2)then
                        Map(i) = 0
                  end if
            end do
            print*, map

            do p = 1, NI+NA
                  do q = 1, NI+NA
                        do r = 1, NI+NA
                              do s = 1, NI+NA
                                    call real_vw_x(this, Rkab(:, p, r), Rkcd(:, q, s), NCholesky)
                                    ETot = ETot + FRDM2(p, q, r, s, R00, Occ, Map, NA, NBasis) &
                                          * this
                              end do
                        end do
                  end do
            end do
            print*, 'etot', etot+enuc, enuc, etot
            stop





            do p = 1, NBasis
                  do q = 1, NBasis
                        do r = 1, NBasis
                              do s = 1, NBasis
                                    call real_vw_x(this, Rkab(:, p, q), Rkcd(:, r, s), NCholesky)
                                    if (abs(this).gt.1.d-5)then
                                          write(*, '(4I5, F20.15)') p, q, r, s, this
                                    end if
                              end do
                        end do
                  end do
            end do

            stop


            print*, nbasis
            allocate(H0_extao(nbasis, nbasis))

            ! !call ints1e_OverlapMatrix(overlap, AObasis)                                                                                                                                                                                                                                
            call ints1e_gammcor_H0_extao(H0_extao, AObasis, System, ExternalOrdering)

            do i = 1, nbasis
                  do j  = 1, nbasis
                        if (abs(H0_extao(i,j)).gt.1.d-5)then
                              write(*, '(I5, A5, I5, A5, F20.15)')i, ',', j, ',', H0_extao(i,j)
                        end if
                  end do
            end do
            stop
            ! ! suma op e. kin i op culomb oddz el-jadr                                                                                                                                                                                                                                   
            ! deallocate(overlap)


      end subroutine test_pluszowy




      Subroutine PP0_JK_loop(ASing, ATrip, H_type, AuxAA, BuxAA, Aux3A, NBasis, NA, NI, NV, IGem, Occ, &
            HNO, AuxCoeff, AuxInd, pos, &
            RDM2_pppp_perm, RDM2_pmpm_perm, SpinSymm, IntJFile, IntKFile, IntVV, ACAlpha)

            double precision, dimension(:,:), intent(inout) :: AuxAA, BuxAA
            double precision, dimension(:,:), intent(inout) :: Aux3A

            integer, intent(in) ::H_type
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:, :), intent(inout) :: HNO
            double precision, dimension(:), intent(in) :: Occ
            integer, intent(in) :: SpinSymm
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, dimension(:),intent(in) :: IGem
            double precision, dimension(:,:, :, :), intent(in) :: AuxCoeff
            integer, dimension(:,:), intent(in) :: AuxInd(3,3), pos
            double precision, dimension(:,:,:,:), intent(in) :: RDM2_pppp_perm, RDM2_pmpm_perm
            character(*)                 :: IntJFile,IntKFile,IntVV
            double precision, intent(in) :: ACAlpha
            integer :: iunit1, iunit2, iunit3
            double precision, dimension(:), allocatable :: TwoNOA, Ha
            double precision :: AuxVal
            double precision,allocatable :: ints(:,:)
            integer :: k, l, kl, i, j
            integer :: kk, ll, ii, jj
            integer :: t, u, x
            integer :: ipq, irs
            integer :: ip, iq, ir, is
            integer :: NIA
            double precision :: HNOCoef, val, temp, temp1, val1, val2, numf
            double precision :: num_spin

            allocate(ints(NBasis, Nbasis))

            AuxAA = zero
            BuxAA = zero
            Aux3A = zero

            NIA = NI + NA
            print*, NI, NA, NV, NBasis, 'sniez'
            HNOCoef = One- ACAlpha
            print*, 'NBasis', NBasis


            !-------------------------------------------------------------------------
            open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis*NIA)

            do l=1, NIA
                  do k=1,NBasis

                        kl = k + NBasis * (l-1)

                        read(iunit1, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NIA)

                        ints(:,NIA+1:NBasis) = 0

                        ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

                        !exchange
                        val = HNOCoef*Occ(l)
                        if(IGem(l)==1) then
                              if (H_type==pDyall)then
                                    if(IGem(k)==1) then
                                          do i=1, NI
                                                HNO(i,k) = HNO(i,k) - val*ints(i,l)
                                          end do
                                    end if
                                    if(IGem(k)==2) then
                                          do i=NI+1,NIA
                                                HNO(i,k) = HNO(i,k) - val*ints(i,l)
                                          end do
                                    else if(IGem(k)==3) then
                                          do i=NIA+1,NBasis
                                                HNO(i,k) = HNO(i,k) - val*ints(i,l)
                                          end do
                                    end if
                              end if
                        end if

                        if(IGem(l)==2) then
                              if(IGem(k)==1) then
                                    do i=1,NI
                                          HNO(i,k) = HNO(i,k) - val*ints(i,l)
                                    end do
                              else if(IGem(k)==3) then
                                    do i=NIA+1,NBasis
                                          HNO(i,k) = HNO(i,k) - val*ints(i,l)
                                    end do
                              end if
                        end if
                  end do
            end do

            close(iunit1)

            open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis*NIA)

            intk_l_loop: do l=NI+1, NIA
                  intk_k_loop: do k=NI+1,NIA

                        kl = k + NBasis * (l-1)

                        read(iunit1, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NIA)

                        ints(:,NIA+1:NBasis) = 0

                        ip = k
                        t  = l

                        if (IGem(t)==2)then
                              do ir = NI+1, NIA
                                    AuxAA(ip, ir) = AuxAA(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(ir,t)
                              end do
                        end if

                        ip = k
                        ir = l
                        do	t =NI+1, NIA
                              BuxAA(ip, ir) = BuxAA(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(t,t)
                        end do

                        ! ------
                        ! n_q n_s (ps|rq)   P3a-3

                        ir = k
                        iq = l

                        do ip = iq, NIA
                              !                if (ip==iq) numf=numf*sqrt(frac12)
                              do is = NI+1, ir     ! r >=s
                                    !                  if (ir==is) numf=numf*sqrt(frac12)
                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)

                                    if (ipq>0.and.irs>0)then
                                          AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                                          val = zero
                                          if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq) * Occ(is)

                                          ASing(irs, ipq) = ASing(irs, ipq) - AuxVal * val * ints(ip, is)

                                          ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ip, is)
                                    end if
                              end do
                        end do


                        ir = k
                        iq = l

                        do ip =	iq, NIA
                              do is = NI+1, NIA

                                    ! P3c-5,6
                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)
                                    if (ipq>0.and.irs>0)then
                                          val = AuxCoeff(IGem(ir),IGem(iq),2,2)* &
                                                (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,ip)*ints(NI+1:NIA,NI+1:NIA)) + &
                                                sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,ip)*ints(NI+1:NIA,NI+1:NIA)))

                                          ASing(irs, ipq) = ASing(irs, ipq)  - val
                                          ATrip(irs, ipq) = ATrip(irs, ipq)  + val
                                    end if
                              end do
                        end do

                        ip = k
                        is = l

                        do ir = is, NIA
                              do iq = NI+1, NIA
                                    !P3c-7,8
                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)
                                    if (ipq>0.and.irs>0)then
                                          val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                                                (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,ir,iq)*ints(NI+1:NIA,NI+1:NIA)) + &
                                                sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,ir,iq)*ints(NI+1:NIA,NI+1:NIA)))

                                          ASing(irs, ipq) = ASing(irs, ipq)  - val
                                          ATrip(irs, ipq) = ATrip(irs, ipq)  + val
                                    end if
                              end do
                        end do


                        ir = k
                        t = l

                        do iq =	NI+1, NIA
                              do is = NI+1, NIA
                                    do ip = iq, NIA
                                          !P3a-3 P3b-1
                                          ipq = pos(ip, iq)
                                          irs = pos(ir, is)
                                          if (irs>0.and.ipq>0)then
                                                val1 = zero
                                                val2 = zero
                                                val1 = - AuxCoeff(IGem(ip),IGem(t),IGem(ir),2)* sum(RDM2_pmpm_perm(NI+1:NIA,is,iq,t)*ints(ip, NI+1:NIA)) !P3a-3
                                                val2 =   AuxCoeff(IGem(ip),IGem(t),IGem(ir),2) * sum(RDM2_pmpm_perm(NI+1:NIA,t,iq,is)*ints(ip, NI+1:NIA)) !P3b-1
                                                ASing(irs, ipq) = ASing(irs, ipq)  - val1 + val2
                                                ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                                          end if
                                    end do
                              end do
                        end do


                        ! W Aux3A sa wyrazy z P4 i P5
                        if(l<=NIA) then
                              do ip=1,NIA                 
                                    val = zero
                                    val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                                          sum(ints(NI+1:NIA,NI+1:NIA)*RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,l,ip))

                                    val1 = zero
                                    val1 = val1 + AuxCoeff(IGem(k),IGem(l),2,2)* &
                                          sum(ints(NI+1:NIA,NI+1:NIA)*RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,l, ip))

                                    Aux3A(k,ip) = Aux3A(k,ip) + val + val1
                              enddo
                        endif


                  end do intk_k_loop
            end do intk_l_loop



            close(iunit1)

            open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis**2)

            do l=1,NIA
                  do k=1,NIA
                        kl = k + NIA * (l-1)

                        read(iunit2, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NBasis)


                        if(k==l)then
                              val = 2*HNOCoef*Occ(k)
                              if (k<=NI)then
                                    do ip = NI+1, NIA
                                          do iq = NI+1, NIA
                                                HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                                          end do
                                    end do
                              end if
                              do ip = 1, NI
                                    do iq = 1, NI
                                          HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                                    end do
                              end do
                              do ip = NIA+1, NBasis
                                    do iq = NIA+1, NBasis
                                          HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                                    end do
                              end do
                        end if
                  end do
            end do

            close(iunit2)

            open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis**2)

            intj_l_loop: do l=NI+1,NIA
                  intj_k_loop: do k=NI+1,NIA
                        kl = k + NIA * (l-1)

                        read(iunit2, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NBasis)


                        ! ------
                        ! n_p n_r (ps|rq)  P3a-4
                        ! n_p n_s (ps|rq)  P3c-5

                        ip = k
                        is = l

                        do ir = is, NIA
                              do iq = NI+1, ip 

                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)
                                    if (ipq>0.and.irs>0)then

                                          AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                                          val = zero
                                          if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip) * Occ(is)

                                          if((AuxInd(IGem(ip),IGem(ir))==1).and.(ir <= NIA)) val = val - Occ(ip)*Occ(ir)


                                          ASing(irs, ipq) = ASing(irs, ipq) - AuxVal * val * ints(ir, iq)
                                          ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ir, iq)

                                    end if
                              end do
                        end do


                        ip = k
                        ir = l

                        if(IGem(ip)==2.and.Igem(ir)==2)then

                              do iq = NI+1, ip
                                    do is = NI+1, ir
                                          !P3c-1,2
                                          ipq = pos(ip,iq)
                                          irs = pos(ir, is)
                                          if (ipq>0.and.irs>0)then

                                                val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                                      (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-1,2
                                                      sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)))
                                                ASing(irs, ipq) = ASing(irs, ipq)  - val
                                                ATrip(irs, ipq) = ATrip(irs, ipq)  - val
                                          end if

                                    end do
                              end do

                        end if

                        if (k==l)then

                              iq = k
                              is = l

                              do ip = max(NI+1, iq), NIA
                                    do ir = max(NI+1, is), NIA
                                          !P3c-3,4
                                          ipq = pos(ip,iq)
                                          irs = pos(ir, is)
                                          if (ipq>0.and.irs>0)then

                                                val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                                                      (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,ir,ip)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-3,4
                                                      sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,ir,ip)*ints(NI+1:NIA,NI+1:NIA)))
                                                ASing(irs, ipq) = ASing(irs, ipq)  - val
                                                ATrip(irs, ipq) = ATrip(irs, ipq)  - val
                                          end if

                                    end do
                              end do


                              is = k
                              u = l

                              do iq = NI+1, NIA
                                    do ip = iq, NBasis
                                          do ir = max(NI+1, is), NIA

                                                ipq = pos(ip,iq)
                                                irs = pos(ir, is)
                                                if (ipq>0.and.irs>0)then

                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = AuxCoeff(IGem(ip),IGem(is),IGem(u),2)* &
                                                            (sum(RDM2_pmpm_perm(NI+1:NIA,ir,iq,u)*ints(NI+1:NIA,ip))) !P3a-2
                                                      val2 = - AuxCoeff(IGem(ip),IGem(is),IGem(u),2)* &
                                                            sum(RDM2_pmpm_perm(NI+1:NIA,u,iq,ir)*ints(NI+1:NIA,ip)) !P3b-4

                                                      ASing(irs, ipq) = ASing(irs, ipq)  + val1 - val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                                                end if
                                          end do
                                    end do
                              end do


                              ! do iq = 1, NIA
                              !       do ip = max(iq, NI+1), NIA
                              !             do ir = max(NI+1, is), NIA
                              !                   ipq = pos(ip,iq)
                              !                   irs = pos(ir, is)
                              !                   if (ipq>0.and.irs>0)then

                              !                         val1 = zero
                              !                         val2 = zero
                              !                         val1 = -AuxCoeff(IGem(iq),IGem(is),IGem(u),2)* sum(RDM2_pmpm_perm(NI+1:NIA,ir,ip,u)*ints(NI+1:NIA,iq)) !P3a-4
                              !                         val2 = AuxCoeff(IGem(iq),IGem(is),IGem(u),2)* sum(RDM2_pmpm_perm(NI+1:NIA,u,ip,ir)*ints(NI+1:NIA,iq)) !P3b-2

                              !                         ASing(irs, ipq) = ASing(irs, ipq)  - val1 + val2
                              !                         ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                              !                   end if
                              !             end do
                              !       end do
                              ! end do

                              iq = k
                              t = l

                              do is = NI+1, NIA
                                    do ip = max(NI+1, iq), NIA
                                          do ir = is, NIA

                                                ipq = pos(ip,iq)
                                                irs = pos(ir, is)
                                                if (ipq>0.and.irs>0)then

                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = AuxCoeff(IGem(iq),IGem(ir),IGem(t),2)*sum(RDM2_pmpm_perm(NI+1:NIA,ip,is,t)*ints(NI+1:NIA,ir)) !P3a-1
                                                      val2 = -AuxCoeff(IGem(iq),IGem(ir),IGem(t),2)*sum(RDM2_pmpm_perm(NI+1:NIA,t,is,ip)*ints(NI+1:NIA,ir)) !P3b-3

                                                      ASing(irs, ipq) = ASing(irs, ipq)   + val1 - val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)   + val1 + val2

                                                end if

                                          end do
                                    end do
                              end do

                        end if


                        ll = l
                        jj = k

                        if (jj<=ll)then
                              call bare_int_loops(ll, jj, pos, NI+1, NIA, ASing, ATrip, ints, AuxCoeff, Occ)
                        end if


                  end do intj_k_loop
            end do intj_l_loop
            !       !$OMP END PARALLEL DO                                                                                                                                                                                      

            close(iunit2)


      end Subroutine PP0_JK_loop



      Subroutine PP_JK_loop(ASing, ATrip, H_type, AuxII, AuxAA, BuxII, BuxAA, Aux3A, NBasis, NA, NI, NV, IGem, Occ, &
            HNO, AuxCoeff, AuxInd, pos, &
            RDM2_pppp_perm, RDM2_pmpm_perm, SpinSymm, IntJFile, IntKFile, IntVV, ACAlpha)

            double precision, dimension(:,:), intent(inout) :: AuxII, AuxAA, BuxII, BuxAA
            double precision, dimension(:,:), intent(inout) :: Aux3A

            integer, intent(in) ::H_type
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:, :), intent(inout) :: HNO
            double precision, dimension(:), intent(in) :: Occ
            integer, intent(in) :: SpinSymm
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, dimension(:),intent(in) :: IGem
            double precision, dimension(:,:, :, :), intent(in) :: AuxCoeff
            integer, dimension(:,:), intent(in) :: AuxInd(3,3), pos
            double precision, dimension(:,:,:,:), intent(in) :: RDM2_pppp_perm, RDM2_pmpm_perm
            character(*)                 :: IntJFile,IntKFile,IntVV
            double precision, intent(in) :: ACAlpha
            integer :: iunit1, iunit2, iunit3
            double precision, dimension(:), allocatable :: TwoNOA, Ha
            double precision :: AuxVal
            double precision,allocatable :: ints(:,:)
            real :: start_time, end_time, total_io_time

            type (tclock) :: timer, timer0
            integer :: k, l, kl, i, j
            integer :: kk, ll, ii, jj
            integer :: t, u, x
            integer :: ipq, irs
            integer :: ip, iq, ir, is
            integer :: NIA
            double precision :: HNOCoef, val, temp, temp1, val1, val2, numf

            allocate(ints(NBasis, Nbasis))

            AuxAA = zero
            BuxAA = zero
            AuxII = zero
            BuxII = zero
            Aux3A = zero

            NIA = NI + NA
            print*, NI, NA, NV, NBasis, 'sniez'
            HNOCoef = One- ACAlpha
            print*, 'NBasis', NBasis


            !-------------------------------------------------------------------------
            open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis*NIA)

            call clock_start(timer)
            kl = 0
            intk_l_loop: do l=1,NIA
                  intk_k_loop: do k=1,NBasis
                        kl = kl + 1
                        read(iunit1, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NIA)

                        ints(:,NIA+1:NBasis) = 0


                        ll = k
                        jj = l

                        if (ll>NIA)then
                              if (jj<=ll)then
                                    call bare_int_loops_fofo(ll, jj, pos, NI, NA, NBasis, ASing, ATrip, ints, AuxCoeff, Occ)
                              end if
                        end if


            !             ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

            !             !exchange
            !             val = HNOCoef*Occ(l)
            !             if(IGem(l)==1) then
            !                   if (H_type==pDyall)then
            !                         if(IGem(k)==1) then
            !                               do i=1, NI
            !                                     HNO(i,k) = HNO(i,k) - val*ints(i,l)
            !                               end do
            !                         end if
            !                         if(IGem(k)==2) then
            !                               do i=NI+1,NIA
            !                                     HNO(i,k) = HNO(i,k) - val*ints(i,l)
            !                               end do
            !                         else if(IGem(k)==3) then
            !                               do i=NIA+1,NBasis
            !                                     HNO(i,k) = HNO(i,k) - val*ints(i,l)
            !                               end do
            !                         end if
            !                   end if
            !             end if

            !             if(IGem(l)==2) then
            !                   if(IGem(k)==1) then
            !                         do i=1,NI
            !                               HNO(i,k) = HNO(i,k) - val*ints(i,l)
            !                         end do
            !                   else if(IGem(k)==3) then
            !                         do i=NIA+1,NBasis
            !                               HNO(i,k) = HNO(i,k) - val*ints(i,l)
            !                         end do
            !                   end if
            !             end if


            !             ip = k
            !             t  = l

            !             if (IGem(t)==1)then
            !                   do ir = 1, NBasis
            !                         AuxII(ip, ir) = AuxII(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(ir,t)
            !                   end do

            !             else if (IGem(t)==2)then
            !                   do ir = 1, NBasis
            !                         AuxAA(ip, ir) = AuxAA(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(ir,t)
            !                   end do
            !             end if


            !             ip = k
            !             ir = l
            !             do t =1, NI
            !                   BuxII(ip, ir) = BuxII(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(t,t)
            !             end do
            !             do	t =NI+1, NIA
            !                   BuxAA(ip, ir) = BuxAA(ip, ir) + AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(t,t)
            !             end do


            !             ! ------

                        ! n_q n_s (ps|rq)   P3a-3

                        ir = k
                        iq = l
                        numf = one

                        do ip = iq, NBasis   ! p>=q
                              if (ip==iq) numf=numf*sqrt(frac12)
                              do is = 1, ir     ! r >=s
                                    if (ir==is) numf=numf*sqrt(frac12)
                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)

                                    if (ipq>0.and.irs>0)then
                                          AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                                          val = zero
                                          if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq) * Occ(is)

                                          ASing(irs, ipq) = ASing(irs, ipq) - AuxVal * val * ints(ip, is)
                                          ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ip, is)
                                    end if
                              end do
                        end do
                        

            !             if(k>NI)then

            !                   ir = k
            !                   iq = l

            !                   do ip =	iq, NIA
            !                         do is = NI+1, NIA

            !                               ! P3c-5,6
            !                               ipq = pos(ip,iq)
            !                               irs = pos(ir, is)
            !                               if (ipq>0.and.irs>0)then
            !                                     val = AuxCoeff(IGem(ir),IGem(iq),2,2)* &
            !                                           (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,ip)*ints(NI+1:NIA,NI+1:NIA)) + &
            !                                           sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,ip)*ints(NI+1:NIA,NI+1:NIA)))
            !                                     ASing(irs, ipq) = ASing(irs, ipq)  - val
            !                                     ATrip(irs, ipq) = ATrip(irs, ipq)  + val
            !                               end if
            !                         end do
            !                   end do

            !                   ip = k
            !                   is = l

            !                   do ir = is, NIA
            !                         do iq = NI+1, NIA
            !                               !P3c-7,8
            !                               ipq = pos(ip,iq)
            !                               irs = pos(ir, is)
            !                               if (ipq>0.and.irs>0)then
            !                                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
            !                                           (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,ir,iq)*ints(NI+1:NIA,NI+1:NIA)) + &
            !                                           sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,ir,iq)*ints(NI+1:NIA,NI+1:NIA)))
            !                                     ASing(irs, ipq) = ASing(irs, ipq)  - val
            !                                     ATrip(irs, ipq) = ATrip(irs, ipq)  + val
            !                               end if
            !                         end do
            !                   end do

            !             end if


                        if (k>NI.and.l>NI)then

                              ir = k
                              t = l

                              do iq =	NI+1, NIA
                                    do is = NI+1, NIA
                                          do ip = iq, NBasis
                                                !P3a-3 P3b-1
                                                ipq = pos(ip, iq)
                                                irs = pos(ir, is)
                                                if (irs>0.and.ipq>0)then
                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = - AuxCoeff(IGem(ip),IGem(t),IGem(ir),2)* sum(RDM2_pmpm_perm(NI+1:NIA,is,iq,t)*ints(ip, NI+1:NIA)) !P3a-3
                                                      val2 =   AuxCoeff(IGem(ip),IGem(t),IGem(ir),2) * sum(RDM2_pmpm_perm(NI+1:NIA,t,iq,is)*ints(ip, NI+1:NIA)) !P3b-1
                                                      ASing(irs, ipq) = ASing(irs, ipq)  - val1 + val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                                                end if
                                          end do
                                    end do
                              end do
                        end if


            !             ! W Aux3A sa wyrazy z P4 i P5
            !             if(l<=NIA) then
            !                   do ip=1,NIA                 
            !                         val = zero
            !                         val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
            !                               sum(ints(1:NI,1:NI)*RDM2_pmpm_perm(1:NI,1:NI,l, ip))
            !                         val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
            !                               sum(ints(NI+1:NIA,1:NI)*RDM2_pmpm_perm(NI+1:NIA,1:NI,l, ip))
            !                         val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
            !                               sum(ints(1:NI,NI+1:NIA)*RDM2_pmpm_perm(1:NI,NI+1:NIA,l, ip))
            !                         val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
            !                               sum(ints(NI+1:NIA,NI+1:NIA)*RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,l,ip))

            !                         val1 = zero

            !                         val1 = val1 + AuxCoeff(IGem(k),IGem(l),1,1)* &
            !                               sum(ints(1:NI,1:NI)*RDM2_pppp_perm(1:NI,1:NI,l, ip))

            !                         val1 = val1 + AuxCoeff(IGem(k),IGem(l),2,1)* &
            !                               sum(ints(NI+1:NIA,1:NI)*RDM2_pppp_perm(NI+1:NIA,1:NI,l, ip))
            !                         val1 = val1 + AuxCoeff(IGem(k),IGem(l),1,2)* &
            !                               sum(ints(1:NI,NI+1:NIA)*RDM2_pppp_perm(1:NI,NI+1:NIA,l, ip))
            !                         val1 = val1 + AuxCoeff(IGem(k),IGem(l),2,2)* &
            !                               sum(ints(NI+1:NIA,NI+1:NIA)*RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,l, ip))

            !                         Aux3A(k,ip) = Aux3A(k,ip) + val + val1
            !                   enddo
            !             endif


                  end do intk_k_loop
            end do intk_l_loop

            print*, 'TIME na unit1 ', clock_readwall(timer)

            close(iunit1)
            call clock_start(timer)
            open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
                  access='DIRECT',recl=8*NBasis**2)
            kl = 0
            intj_l_loop: do l=1,NIA
                  intj_k_loop: do k=1,NIA
                        kl = kl + 1
                        read(iunit2, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NBasis)


                        ! if(k==l.and.k<=NIA) then
                        !       val = 2*HNOCoef*Occ(k)
                        !       if (k<=NI)then
                        !             do ip = NI+1, NIA
                        !                   do iq = NI+1, NIA
                        !                         HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                        !                   end do
                        !             end do
                        !       end if
                        !       !                else if (k>NI)then
                        !       do ip = 1, NI
                        !             do iq = 1, NI
                        !                   HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                        !             end do
                        !       end do
                        !       do ip = NIA+1, NBasis
                        !             do iq = NIA+1, NBasis
                        !                   HNO(ip,iq) = HNO(ip, iq) + val * ints (ip,iq)
                        !             end do
                        !       end do
                        !       !               end if
                        ! end if


                        ! ! ------                                                                                                                                                                                                                                                                  
                        ! n_p n_r (ps|rq)  P3a-4
                        ! n_p n_s (ps|rq)  P3c-5

                        ip = k
                        is = l

                        do ir = is, NBasis
                              do iq = 1, ip 

                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)
                                    if (ipq>0.and.irs>0)then

                                          AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                                          val = zero
                                          if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip) * Occ(is)

                                          if((AuxInd(IGem(ip),IGem(ir))==1).and.(ir <= NIA)) val = val - Occ(ip)*Occ(ir)


                                          ASing(irs, ipq) = ASing(irs, ipq) - AuxVal * val * ints(ir, iq)
                                          ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ir, iq)

                                    end if
                              end do
                        end do



                        ! ! ------
                        ! ! n_q n_r (ps|rq) P3c-7

                        ir = k
                        iq = l

                        if(AuxInd(IGem(iq),IGem(ir))==1)then
                              val =  - Occ(iq) * Occ(ir)

                              do ip = iq, NBasis
                                    do is = 1, ir

                                          ipq = pos(ip,iq)
                                          irs = pos(ir, is)
                                          if (ipq>0.and.irs>0)then

                                                AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))


                                                ASing(irs, ipq) = ASing(irs, ipq) - AuxVal * val * ints(ip, is)
                                                ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ip, is)
                                          end if
                                    end do
                              end do
                        end if



                        ! ! ------                                                                                                                                                                                                                                                    
                        ! ! n_q n_s (pr|qs) P3c1
                        ! n_p n_s (pr|qs) P3a-1
                        ! ! n_q n_r (pr|qs) P3a-2

                        iq = k
                        is = l

                        do ip = iq, NBasis
                              do ir = is, NBasis

                                    ipq = pos(ip,iq)
                                    irs = pos(ir, is)

                                    if (ipq>0.and.irs>0)then
                                          AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                                          val = zero
                                          if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq) * Occ(is)
                                          if((AuxInd(IGem(ip),IGem(is))==1).and.(ip <= NIA)) val = val + Occ(ip)*Occ(is)
                                          if((AuxInd(IGem(iq),IGem(ir))==1).and.(ir <= NIA)) val = val + Occ(iq)*Occ(ir)
                                          ASing(irs, ipq) = ASing(irs, ipq) + AuxVal * val * ints(ip, ir)
                                          ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(ip, ir)
                                    end if
                              end do
                        end do


                        ! ------
                        ! n_p n_r (pr|qs)      P3c -3                                                                                                                                                                                                                                               

                        ip = k
                        ir = l

                        if(AuxInd(IGem(ip),IGem(ir))==1) then

                              val  =  Occ(ip) * Occ(ir)

                              do iq = 1, ip
                                    do is = 1, ir

                                          ipq = pos(ip,iq)
                                          irs = pos(ir, is)

                                          if (ipq>0.and.irs>0)then
                                                AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))


                                                ASing(irs, ipq) = ASing(irs, ipq) + AuxVal * val * ints(iq, is)
                                                ATrip(irs, ipq) = ATrip(irs, ipq) + AuxVal * val * ints(iq, is)
                                          end if
                                    end do
                              end do
                        end if

                        ! ip = k
                        ! ir = l

                        ! if(IGem(ip)==2.and.Igem(ir)==2)then

                        !       do iq = NI+1, ip
                        !             do is = NI+1, ir
                        !                   !P3c-1,2

                        !                   ipq = pos(ip,iq)
                        !                   irs = pos(ir, is)
                        !                   if (ipq>0.and.irs>0)then

                        !                         val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                        !                               (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-1,2
                        !                               sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)))
                        !                         ASing(irs, ipq) = ASing(irs, ipq)  - val
                        !                         ATrip(irs, ipq) = ATrip(irs, ipq)  - val
                        !                   end if

                        !             end do
                        !       end do


                        ! end if

                        ! if (k==l)then

                        !       t = l
                        !       do ip =1, NBasis
                        !             do ir = NIA+1, NBasis
                        !                   val = AuxCoeff(IGem(ip),IGem(t),IGem(ir),IGem(t))* Occ(t) * ints(ip,ir)
                        !                   if(t<=NI)BuxII(ip, ir) = BuxII(ip, ir) + val
                        !                   if(t>NI)BuxAA(ip, ir) = BuxAA(ip, ir) + val
                        !             end do
                        !       end do
                        ! end if

                        ! !             if(k>NI.and.l>NI)then

                        ! iq = k
                        ! is = l

                        ! do ip = max(NI+1, iq), NIA
                        !       do ir = max(NI+1, is), NIA
                        !             !P3c-3,4
                        !             ipq = pos(ip,iq)
                        !             irs = pos(ir, is)
                        !             if (ipq>0.and.irs>0)then

                        !                   val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                        !                         (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,ir,ip)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-3,4
                        !                         sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,ir,ip)*ints(NI+1:NIA,NI+1:NIA)))
                        !                   ASing(irs, ipq) = ASing(irs, ipq)  - val
                        !                   ATrip(irs, ipq) = ATrip(irs, ipq)  - val
                        !             end if

                        !       end do
                        ! end do

                        !HERE
                        if (l>NI)then
                              is = k
                              u = l

                              do iq = NI+1, NIA
                                    do ip = iq, NBasis
                                          do ir = max(NI+1, is), NIA

                                                ipq = pos(ip,iq)
                                                irs = pos(ir, is)
                                                if (ipq>0.and.irs>0)then

                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = AuxCoeff(IGem(ip),IGem(is),IGem(u),2)* &
                                                            (sum(RDM2_pmpm_perm(NI+1:NIA,ir,iq,u)*ints(NI+1:NIA,ip))) !P3a-2
                                                      val2 = - AuxCoeff(IGem(ip),IGem(is),IGem(u),2)* &
                                                            sum(RDM2_pmpm_perm(NI+1:NIA,u,iq,ir)*ints(NI+1:NIA,ip)) !P3b-4

                                                      ASing(irs, ipq) = ASing(irs, ipq)  + val1 - val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                                                end if
                                          end do
                                    end do
                              end do


                              do iq = 1, NIA
                                    do ip = max(iq, NI+1), NIA
                                          do ir = max(NI+1, is), NIA
                                                ipq = pos(ip,iq)
                                                irs = pos(ir, is)
                                                if (ipq>0.and.irs>0)then

                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = -AuxCoeff(IGem(iq),IGem(is),IGem(u),2)* sum(RDM2_pmpm_perm(NI+1:NIA,ir,ip,u)*ints(NI+1:NIA,iq)) !P3a-4
                                                      val2 = AuxCoeff(IGem(iq),IGem(is),IGem(u),2)* sum(RDM2_pmpm_perm(NI+1:NIA,u,ip,ir)*ints(NI+1:NIA,iq)) !P3b-2

                                                      ASing(irs, ipq) = ASing(irs, ipq)  - val1 + val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)  + val1 + val2
                                                end if
                                          end do
                                    end do
                              end do


                              !ZHERE

                              iq = k
                              t = l

                              do is = NI+1, NIA
                                    do ip = max(NI+1, iq), NIA
                                          do ir = is, NBasis

                                                ipq = pos(ip,iq)
                                                irs = pos(ir, is)
                                                if (ipq>0.and.irs>0)then

                                                      val1 = zero
                                                      val2 = zero
                                                      val1 = AuxCoeff(IGem(iq),IGem(ir),IGem(t),2)*sum(RDM2_pmpm_perm(NI+1:NIA,ip,is,t)*ints(NI+1:NIA,ir)) !P3a-1
                                                      val2 = -AuxCoeff(IGem(iq),IGem(ir),IGem(t),2)*sum(RDM2_pmpm_perm(NI+1:NIA,t,is,ip)*ints(NI+1:NIA,ir)) !P3b-3

                                                      ASing(irs, ipq) = ASing(irs, ipq)   + val1- val2
                                                      ATrip(irs, ipq) = ATrip(irs, ipq)   + val1 + val2
                                                      ! if(ip==2.and.iq==1.and.ir==2.and.is==2)then
                                                      !       write(*, '(6I5, )')
                                                end if
                                                
                                          end do
                                    end do
                              end do

                        end if


                        ll = l
                        jj = k

                        if (jj<=ll)then
                              call bare_int_loops(ll, jj, pos, 1, NBasis, ASing, ATrip, ints, AuxCoeff, Occ)
                        end if


                  end do intj_k_loop
            end do intj_l_loop
            !       !$OMP END PARALLEL DO                                                                                                                                                                                      

            close(iunit2)

            print*, 'TIME na unit2 ', clock_readwall(timer)

            call clock_start(timer)

            open(newunit=iunit3,file=trim(IntVV),status='OLD', &                                                                                                                                                      
                  access='DIRECT',recl=8*NBasis*NBasis)                                                                                                                                                                

            kl = 0
            total_io_time = zero
            intk_l_loop3: do l=NI+1, NBasis !NIA + 1, NBasis                                                                                                                                                                     
                  intk_k_loop3: do k=NIA+1,NBasis                                                                                                                                                                     
                        kl =  k + NBasis * (l-1)
                        call cpu_time(start_time)  
                        read(iunit3, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NBasis)
                        call cpu_time(end_time)
                        total_io_time = total_io_time + (end_time - start_time)

                        if (l>NIA)then

                              ll = l
                              jj = k

                              if (jj<=ll)then
                                    write(*,'(4I5)' )k, l
                                    call bare_int_loops(ll, jj, pos, 1, NBasis, ASing, ATrip, ints, AuxCoeff, Occ)
                              end if


            !                   ip = k
            !                   ir = l

            !                   do iq =	NI+1, NIA
            !                         do is = NI+1, NIA
            !                               !P3c-1,2                                                                                                                                                                                    
            !                               ipq = pos(ip,iq)
            !                               irs = pos(ir, is)
            !                               if (ipq>0.and.irs>0)then                         
            !                                     val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
            !                                           (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-1,2
            !                                           sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)))
            !                                     ASing(irs, ipq) = ASing(irs, ipq)  - val
            !                                     ATrip(irs, ipq) = ATrip(irs, ipq)  - val
            !                               end if
            !                         end do
            !                   end do
            !             else

            !                   ip = k
            !                   ir = l

            !                   do iq = NI+1, NIA
            !                         do is = NI+1, min(NIA, l)
            !                               !P3c-1,2
            !                               ipq = pos(ip,iq)
            !                               irs = pos(ir, is)
            !                               if (ipq>0.and.irs>0)then
            !                                     val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
            !                                           (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-1,2                                                                                       
            !                                           sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)))
            !                                     ASing(irs, ipq) = ASing(irs, ipq)  - val
            !                                     ATrip(irs, ipq) = ATrip(irs, ipq)  - val
            !                               end if
            !                         end do
            !                   end do

            !                   ip = l
            !                   ir = k

            !                   do iq = NI+1, min(NIA, l)
            !                         do is = NI+1, NIA
            !                               !P3c-1,2                                                                                                                                                                                     
            !                               ipq = pos(ip,iq)
            !                               irs = pos(ir, is)
            !                               if (ipq>0.and.irs>0)then
            !                                     val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
            !                                           (sum(RDM2_pmpm_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)) + & !P3c-1,2
            !                                           sum(RDM2_pppp_perm(NI+1:NIA,NI+1:NIA,is,iq)*ints(NI+1:NIA,NI+1:NIA)))
            !                                     ASing(irs, ipq) = ASing(irs, ipq)  - val
            !                                     ATrip(irs, ipq) = ATrip(irs, ipq)  - val
            !                               end if
            !                         end do
            !                   end do

                        end if

                  end do intk_k_loop3
            end do intk_l_loop3

            close(iunit3)

            ! print*, 'TIME na unit3 ', clock_readwall(timer)
            ! print *, "Całkowity czas spędzony na operacjach I/O: ", total_io_time, " sekund"  

      end Subroutine PP_JK_loop


      subroutine bare_int_loops(ll, jj, pos, N0, N1, ASing, ATrip, ints, AuxCoeff, Occ)

            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: ll, jj, N0, N1
            integer :: kk, ii


            if (jj==ll)then

                  ! w1
                  kk = ll
                  ii = ll
                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)

                  do ii = N0, ll-1
                        !w2
                        kk = ll
                        call update_ASing5(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  end do

                  do kk = ll+1, N1
                        !w4
                        ! print*, 'w4', kk, ii, ll, jj
                        ! print*, ''
                        ii = ll
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing3(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        ! print*, ''
                        !w6

                        ii = kk
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)


                        do ii = N0, ll-1
                              !w14
                              call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              !                     call update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        end do

                        do ii = ll + 1, kk - 1
                              !w15
                              call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                              call update_ASing3(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        end do

                  end do

            else
                  !w5
                  kk = ll
                  ii = jj
                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                  call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)

                  do ii = jj+1, ll-1
                        !w8
                        kk = ll
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing5(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  end do

                  do kk = ll + 1, N1
                        !w9
                        ii = kk
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        !w16
                        ii = ll
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing3(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)

                        !w13
                        ii = jj
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)

                        do ii = N0, jj-1
                              !w10
                              call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              call update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        end do

                        do ii = jj+1, ll-1
                              !w11
                              call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                              call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                              call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                              call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        end do

                        do ii = ll+1, kk-1
                              !w12
                              call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                              call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                              call update_ASing3(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                              call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        end do
                  end do
            end if

      end subroutine bare_int_loops

      subroutine bare_int_loops_fofo(ll, jj, pos, NI, NA, NBasis, ASing, ATrip, ints, AuxCoeff, Occ)

            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:), intent(in) :: Occ
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: ll, jj, NI, NA, NBasis
            integer :: kk, ii, NIA
            NIA = NI + NA


            !w5
            ! print*, 'w5'
            kk = ll
            ii = jj
            call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
            call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

            call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
            call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
            call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)

            do ii = jj+1, min(NIA, ll-1)
                  !w8
                  ! print*, 'w8'
                  kk = ll
                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                  call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                  call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing5(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  !               print*, 'w8', ASing(56,56)

            end do

            do kk = ll + 1, NBasis
                  !w13
                  ! print*, 'w13'
                  ii = jj
                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                  call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                  call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  call update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                  !                print*, 'w13', ASing(56,56)

                  do ii = 1, jj-1
                        !w10
                        ! print*, 'w10'
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        !                    print*, 'w10', ASing(56,56)

                  end do

                  do ii = jj+1, min(NIA, ll-1)
                        !w11
                        ! print*, 'w11'
                        call update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)
                        call update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 0, Occ)

                        call update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        call update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, 1, Occ)
                        !                     print*, 'w11', ASing(56,56)

                  end do

            end do


      end subroutine bare_int_loops_fofo


      subroutine ACPP0_fast_fofo(H_type, ECorr_AC0, &
            ETot, ENuc, Occ, XOne, TwoNO, AuxData, IndN, IndAux, IndMod, &
            NBasis, NA, NI, NV, NInte1, NInte2, Flags, AC_TYPE, UMOAO, IntJFile, IntKFile, IntVV)

            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: H_type
            double precision, intent(inout) :: ECorr_AC0, ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: Occ
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            type(TACppData), intent(in) :: AuxData
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), allocatable :: Map
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            type (tclock) :: timer, timer0
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: AC_TYPE
            character(*)                 :: IntJFile,IntKFile, IntVV
            double precision, dimension(:,:), intent(in) :: UMOAO
            double precision, dimension(:, :), allocatable :: HNO, HNO0

            integer :: SpinSymm
            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9

            character(len=3), dimension(10) :: bl_name
            character(len=5), dimension(9)  :: bbl_name
            integer, dimension(2,9) :: block_pair
            integer :: p, q, r, s, i, j, k, l
            integer :: twoint_dim
            integer :: dimw, spsym
            integer :: ir, il, pq, rs
            integer, dimension(10) :: dim_st
            double precision, dimension(9) :: E_contr_s, E_contr_t
            type(eigsBlockParams), dimension(10) :: eigsPA

            ! fofo
            integer          :: Ind(NBasis) !, IGem(NBasis)
            integer :: NIA
            integer :: NRDM2, NRDM2Act

            double precision, dimension(:,:), allocatable :: MxS, ASing, pluszon, ATrip, Ure
            type(TRdmData) :: RdmData
            double precision :: AuxCoeff(3,3,3,3)
            integer :: AuxInd(3,3), pos(NBasis,NBasis), pos_block_s(NBasis, NBasis), pos_block_t(NBasis, NBasis)

            double precision, dimension(:,:), allocatable :: Aux3A, Aux3B
            double precision, dimension(:,:), allocatable :: AuxII, AuxAA, BuxII, BuxAA, work
            double precision, dimension(:), allocatable :: work1, work2
            integer, dimension(2,4) :: limits
            double precision :: ACAlpha, val


            NRDM2 = NBasis**2*(NBasis**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NIA = NI + NA

            allocate(HNO(NBasis, NBasis))
            allocate(HNO0(NBasis, NBasis))
            allocate(Ure(NBasis, NBasis))
            allocate(work1(NBasis**2),work2(NBasis**2))

            allocate (RdmData%R00(NRDM2Act))
            allocate (RdmData%R11(NRDM2Act))

            allocate(RdmData%rdm2_pp(NIA,NIA,NIA,NIA))
            allocate(RdmData%rdm2_pm(NIA,NIA,NIA,NIA))

            allocate(RdmData%rdm2_pp_12(NIA,NIA,NIA,NIA))
            allocate(RdmData%rdm2_pm_12(NIA,NIA,NIA,NIA))

            allocate(RdmData%rdm2_pp_13(NIA,NIA,NIA,NIA))
            allocate(RdmData%rdm2_pm_13(NIA,NIA,NIA,NIA))

            allocate(RdmData%rdm2_pp_act(NA,NA,NA,NA))
            allocate(RdmData%rdm2_pm_act(NA,NA,NA,NA))


            allocate(RdmData%rdm2_pp_12_act(NA,NA,NA,NA))
            allocate(RdmData%rdm2_pm_12_act(NA,NA,NA,NA))

            allocate(RdmData%rdm2_pp_13_act(NA,NA,NA,NA))
            allocate(RdmData%rdm2_pm_13_act(NA,NA,NA,NA))

            allocate(Aux3A(NBasis, NBasis))
            allocate(Aux3B(NBasis, NBasis))
              
            allocate(AuxII(NBasis, NBasis))
            allocate(AuxAA(NBasis, NBasis))
            allocate(BuxII(NBasis, NBasis))
            allocate(BuxAA(NBasis, NBasis))

            Allocate(ASing(AuxData%NDim_s, AuxData%NDim_s))
            Allocate(ATrip(AuxData%NDim_s, AuxData%NDim_s))
            allocate(work(AuxData%NDim_s, AuxData%NDim_s))
            allocate(IGem(NBasis))

            Allocate(pluszon(AuxData%NDim_s, AuxData%NDim_s))
            Allocate(MxS(AuxData%NDim_s, AuxData%NDim_s))


            RdmData%R00 = Zero
            RdmData%R11 = Zero

            AuxII = zero
            BuxII = zero
            AuxAA = zero
            BuxAA = zero
            Aux3A = zero

            call read_2rdm("rdm2.dat", RdmData%R00, NA)
            call read_2rdm("rdms2.dat", RdmData%R11, NA)

            allocate(map(NBasis))

            k = 1
            do i = 1, NI+NA+NV
                  if(IndAux(i)==1)then
                        map(i) = k
                        k = k+1
                  end if
                  if(IndAux(i)==2)then
                        map(i) = 0
                  end if
            end do


            associate(R00 => RdmData%R00, R11 =>RdmData%R11)

              do l=1,NIA
                    do k=1,NIA
                          do j=1,NIA
                                do i=1,NIA
                                      RdmData%rdm2_pp(i, j, k, l) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 0)
                                      RdmData%rdm2_pm(i, j, k, l) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 1)
                                      RdmData%rdm2_pp_12(i, k, j, l) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 0)
                                      RdmData%rdm2_pm_12(i, k, j, l) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 1)
                                      RdmData%rdm2_pp_13(i, l, k, j) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 0)
                                      RdmData%rdm2_pm_13(i, l, k, j) = get2rdm(i, j, k, l, R00, R11, Occ, map, IndAux, NA, 1)
                                enddo
                          enddo
                    enddo
              enddo
            end associate

            do l=NI+1, NIA
                  do k=NI+1, NIA
                        do j=NI+1, NIA
                              do i=NI+1, NIA
                                    RdmData%rdm2_pp_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pp(i, j, k, l)
                                    RdmData%rdm2_pm_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pm(i, j, k, l)

                                    RdmData%rdm2_pp_12_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pp_12(i, j, k, l)
                                    RdmData%rdm2_pm_12_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pm_12(i, j, k, l)
                                    RdmData%rdm2_pp_13_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pp_13(i, j, k, l)
                                    RdmData%rdm2_pm_13_act(map(i), map(j), map(k), map(l)) = RdmData%rdm2_pm_13(i, j, k, l)
                              enddo
                        enddo
                  enddo
            enddo


            URe = Zero
            do i = 1, Nbasis
                  URe(i,i) = One
            end do

            call triang_to_sq(XOne,work1,NBasis)
            call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
            call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO0,NBasis)
            call sq_symmetrize(HNO0,NBasis)
            deallocate(ure)

            ! print*, 'hno'
            ! do i = 1, NBasis
            !       do j = 1, NBasis
            !             if (abs(HNO0(i,j)).gt.1.d-5)then
            !                   print*, i, j, HNO0(i, j)
            !             end if
            !       end do
            ! end do


            !     write(*, '(4I5, 6F20.15)') 9,6, 3,2, RDM2_pm(9,6,3,2), RDM2_pm(3,2,9,6), RDM2_pm(6,9,2,3)
            !     write(*, '(4I5, 6F20.15)') 9,6, 3,2, RDM2_pm(9,3,6,2), RDM2_pm(6,2,9,3), RDM2_pm(3,9,2,6)
            !     stop
            !     print*, 'rdmki'
            !     do l=NI+1,NIA
            !        do k=NI+1,NIA
            !           do j=NI+1,NIA
            !              do i=NI+1,NIA
            !                 if (abs(rdmData%RDM2_pp(i, j, k, l)).gt.1.d-5) then
            !                    if (i.ne.j.and.j.ne.k.and.i.ne.k.and.i.ne.l.and.k.ne.l)then
            !                          write(*, '(4I5, 6F10.6)') i, j, k, l,RdmData%rdm2_pp(i, j, k, l) , RdmData%rdm2_pp_12(i, k, j, l), RdmData%rdm2_pp_13(i, l, k, j)

            !                    end if
            !              end if
            !              if (abs(rdmData%RDM2_pm(i, j, k, l)).gt.1.d-5) then
            !                    if (i.ne.j.and.j.ne.k.and.i.ne.k.and.i.ne.l.and.k.ne.l)then
            !                          write(*, '(4I5, 6F10.6)') i, j, k, l,&
            !                                RdmData%rdm2_pm(i, j, k, l) ,      RdmData%rdm2_pm_12(i, k, j, l), RdmData%rdm2_pm_13(i, l, k, j)
            !                    end if
            !                 end if

            !              enddo
            !           enddo
            !        enddo
            !     enddo


            ! stop    


            NIA = NI+NA

            do i=1,NI
                  IGem(i) = 1
            end do

            do i=NI+1, NIA
                  IGem(i) = 2
            end do

            do i=NIA+1,NBasis
                  IGem(i) = 3
            end do

            ACAlpha = zero

            AuxInd = 0
            AuxInd(1:2,1:2) = 1
            AuxInd(2,2) = 2

            pos = 0
            do i = 1, AuxData%NDim_s
                  pos( AuxData%IndN_s(1,i), AuxData%IndN_s(2,i)) =  AuxData%IndX_s(i)
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             HNO = HNO0
             ASing = zero
             ATrip = zero
             ACAlpha = 0.0356d+0!one

             call init_pperpa(ACAlpha, HNO, AuxCoeff, IGem, NBasis)


            ASing = zero
            ATrip = zero
            AuxAA = zero
            AuxII = zero
            BuxAA = zero
            BuxII = zero
            Aux3A = zero
            Aux3B = zero
            ! call clock_start(timer)

            !------------------testing THC-----------------
            ! print*, 'starePP'
            ! call PP_JK_loop(ASing, ATrip, H_type, AuxII, AuxAA, BuxII, BuxAA, Aux3A, NBasis, NA, NI, NV, IGem, Occ, &
            !       HNO, AuxCoeff, AuxInd, pos, &
            !       RdmData%rdm2_pp_12, RdmData%rdm2_pm_12, SpinSymm, IntJFile, IntKFile, IntVV, ACAlpha)
            !     print*, 'Czas na PP_JK: '//str(clock_readwall(timer),d=2)
            ! print*, 'block'
            ! call PPERPA_block_fofo(ASing, ATrip, MxS, Occ, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B, &
            !       HNO, AuxData%IndN_s, AuxData%IndN_s, map, AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, AuxData%NDim_s, IGem, AuxInd, &
            !       NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, 0, Flags, 1)



            ! !--------------------------THC code--------------
!           ASing = zero
            ATrip = zero
            pluszon = zero
!            print*, 'start THC'

!            call THC_int_loop(ATrip, UMOAO, Flags, AuxData, Occ, IGem, map, AuxCoeff, AuxInd, RdmData, pos)

            call clock_start(timer)


!            call THC_int_loop(pluszon, ATrip, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B, UMOAO, Flags, AuxData, Occ, IGem, map,  AuxInd, RdmData, pos, ACalpha)
            ! print*, ''
            ! print*, 'Time for thc loop: ', clock_readwall(timer)
            ! print*, ''
            ! call clock_start(timer)

            ! print*, 'wynik thc'
            ! do i = 1, AuxData%NDim_s
            !       do j = 1, AuxData%NDim_s
                        
            !             r = AuxData%IndN_s(1, i)
            !             s = AuxData%IndN_s(2, i)
            !             rs = AuxData%IndX_s(i)
                        
            !             p = AuxData%IndN_s(1, j)
            !             q = AuxData%IndN_s(2, j)
            !             pq = AuxData%IndX_s(j)
                        
            !             if (abs(ASing(i, j)-pluszon(i, j)).gt.1.d-5)then
            !                   write(*,'(2I5, A5, 4I5,F30.16)') i, j, '  |  ', p, q, r, s, ASing(i, j)-pluszon(i,j)

            !             end if

            !       end do
            ! end do



            
            ! call PPERPA_block_fofo(pluszon, ATrip, MxS, Occ, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B,&
            !       HNO, AuxData%IndN_s, AuxData%IndN_s, map, AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, AuxData%NDim_s, IGem, AuxInd, &
            !       NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, 0, Flags, 1)

            ! print*, ''
            ! print*, 'Time for construction of A: ', clock_readwall(timer)
            ! print*, ''

            stop



            ! stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            HNO = HNO0
            ASing = zero
            ATrip = zero

            call init_pperpa(ACAlpha, HNO, AuxCoeff, IGem, NBasis)

            call clock_start(timer0)    
            call PP0_JK_loop(ASing, ATrip, H_type, AuxAA, BuxAA, Aux3A, NBasis, NA, NI, NV, IGem, Occ, &
                  HNO, AuxCoeff, AuxInd, pos, &
                  RdmData%rdm2_pp_12, RdmData%rdm2_pm_12, SpinSymm, IntJFile, IntKFile, IntVV, ACAlpha)

            print*, 'Czas na PP0: ', clock_readwall(timer0)

            call clock_start(timer)


            call PPERPA_block_fofo(ASing, ATrip, MxS, Occ, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B, &
                  HNO, AuxData%IndN_s, AuxData%IndN_s, map, AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, AuxData%NDim_s, IGem, AuxInd, &
                  NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, 0, Flags, 1)
            print*, 'Czas na fofo pp-block: ', clock_readwall(timer)




            ! Na tym etapie utworzona jest duza macierz A singletowa ASing, oraz duza trypletowa ATrip.


            !     stop


            ! 1oo, 2vv, 3aa, 4va, 5ao
            ! 6oot, 7vvt, 8aat, 9vat, 10aot

            bbl_name(1) = 'vvoo'
            bbl_name(2) = 'vvao'
            bbl_name(3) = 'vvaa'

            bbl_name(4) = 'vaoo'
            bbl_name(5) = 'vaao'
            bbl_name(6) = 'vaaa'

            bbl_name(7) = 'aaoo'
            bbl_name(8) = 'aaao'
            bbl_name(9) = 'vaao2'


            bl_name(IOO) = 'oo'
            bl_name(IVV) = 'vv'
            bl_name(IAA) = 'aa'
            bl_name(IVA) = 'va'
            bl_name(IAO) = 'ao'

            bl_name(IOOT) = 'oot'
            bl_name(IVVT) = 'vvt'
            bl_name(IAAT) = 'aat'
            bl_name(IVAT) = 'vat'
            bl_name(IAOT) = 'aot'


            allocate(eigsPA(IOO)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVV)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAA)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVA)%IndN(2, NV*NA))
            allocate(eigsPA(IAO)%IndN(2, NA*NI))
            allocate(eigsPA(IOOT)%IndN(2, (NI+1)*NI/2))
            allocate(eigsPA(IVVT)%IndN(2, (NV+1)*NV/2))
            allocate(eigsPA(IAAT)%IndN(2, (NA+1)*NA/2))
            allocate(eigsPA(IVAT)%IndN(2, NV*NA))
            allocate(eigsPA(IAOT)%IndN(2, NA*NI))


            do i = 1, 10
                  eigsPA(i)%IndN = 0
                  eigsPA(i)%dim = 1
            end do

            pos_block_s = 0
            pos_block_t = 0
            do i = 1, AuxData%NDim_s
                  p = IndN(1, i)
                  q = IndN(2, i)

                  if (IndAux(p)==0 .and. IndAux(q)==0) then
                        call updateIndices_fofo(p, q, pos_block_s, pos_block_t, eigsPA(IOO)%IndN, eigsPA(IOOT)%IndN, eigsPA(IOO)%dim, eigsPA(IOOT)%dim)
                  elseif (IndAux(p)==2 .and. IndAux(q)==2) then
                        call updateIndices_fofo(p, q, pos_block_s, pos_block_t, eigsPA(IVV)%IndN, eigsPA(IVVT)%IndN, eigsPA(IVV)%dim, eigsPA(IVVT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==1) then
                        call updateIndices_fofo(p, q, pos_block_s, pos_block_t, eigsPA(IAA)%IndN, eigsPA(IAAT)%IndN, eigsPA(IAA)%dim, eigsPA(IAAT)%dim)

                  elseif (IndAux(p)==2 .and. IndAux(q)==1) then
                        call updateIndices_fofo(p, q, pos_block_s, pos_block_t, eigsPA(IVA)%IndN, eigsPA(IVAT)%IndN, eigsPA(IVA)%dim, eigsPA(IVAT)%dim)
                  elseif (IndAux(p)==1 .and. IndAux(q)==0) then
                        call updateIndices_fofo(p, q, pos_block_s, pos_block_t, eigsPA(IAO)%IndN, eigsPA(IAOT)%IndN, eigsPA(IAO)%dim, eigsPA(IAOT)%dim)
                  end if
            end do


            ! do p = 1, NI
            !    do q = 1, p          
            !       print*, p,q, 'sa na pozycji', pos_block_t(p,q), 'oo'
            !    end do
            ! end do
            ! print*, ''
            ! do p = NI+1, NIA
            !    do q = 1, NI
            !       print*, p,q, 'sa na pozycji', pos_block_t(p,q), 'ao'
            !    end do
            ! end do
            ! print*, ''
            ! do p = NI+1, NIA
            !    do q = NI+1, p
            !       print*, p,q, 'sa na pozycji', pos_block_t(p,q), 'aa'
            !    end do
            ! end do
            ! print*, ''
            ! do p = NIA+1, NBasis
            !    do q = NIA+1, p
            !       print*, p,q, 'sa na pozycji', pos_block_t(p,q), 'vv'
            !    end do
            ! end do
            ! print*, ''
            ! do p = NIA+1, NBasis
            !    do q = NI+1, NIA
            !       print*, p,q, 'sa na pozycji', pos_block_t(p,q), 'va'
            !    end do
            ! end do
            ! print*, ''





            ! stop

            do i = 1, 10
                  eigsPA(i)%dim = eigsPA(i)%dim-1
            end do
            ! 1oo, 2vv, 3aa, 4va, 5ao                                                                                                          
            ! 6oot, 7vvt, 8aat, 9vat, 10aot  
            do i = 1, 10
                  dimw = eigsPA(i)%dim
                  allocate(eigsPA(i)%Eig(dimw))
                  if (i == IOO.or.i==IVV.or.i==IOOT.or.i==IVVT)then
                        allocate(eigsPA(i)%Eigvec(1, 1))
                  else
                        allocate(eigsPA(i)%Eigvec(dimw, dimw))
                  end if
                  allocate(eigsPA(i)%v_plus(dimw))
            end do


            call clock_start(timer)
            do i = 1, 10
                  print*, ''
                  print*, 'teraz blok', bl_name(i)
                  if (i.le.5)then
                        spsym = 0
                  else
                        spsym = 1
                  end if
                  if (eigsPA(i)%dim.gt.0)then
                        call eigs0_block_fofo(i, ASing, ATrip, pos, pos_block_s, pos_block_t, &
                              eigsPA(i)%Eig, eigsPA(i)%Eigvec, eigsPA(i)%v_plus, Occ, &
                              eigsPA(i)%IndN, IndAux, eigsPA(i)%dim, &
                              spsym, Flags)

                        print*, 'wartosci wlasne tego bloku'
                        do j = 1, eigsPA(i)%dim
                              if (abs(eigsPA(i)%Eig(j)).gt.1.d-5)then
                                    print*, j, eigsPA(i)%Eig(j)
                              end if
                        end do

                  else
                        print*, 'block', bl_name(i), 'dimension is 0'
                  end if
            end do
            print*, 'koniec blokow eigs0'
            stop
            print*, 'Czas na cut i diag: ', clock_readwall(timer)
            print*, 'time for all 0 blocks: ', clock_readwall(timer0)

            ! now the A1 matrix

            HNO = HNO0
            ASing = zero
            ATrip = zero
            ACAlpha = one

            call init_pperpa(ACAlpha, HNO, AuxCoeff, IGem, NBasis)

            print*, 'oblicz cale A'
            call clock_start(timer)
            call PP_JK_loop(ASing, ATrip, H_type, AuxII, AuxAA, BuxII, BuxAA, Aux3A, NBasis, NA, NI, NV, IGem, Occ, &
                  HNO, AuxCoeff, AuxInd, pos, &
                  RdmData%rdm2_pp_12, RdmData%rdm2_pm_12, SpinSymm, IntJFile, IntKFile, IntVV, ACAlpha)
            print*, 'Czas na PP_JK: ', clock_readwall(timer)
            call clock_start(timer)

            call PPERPA_block_fofo(ASing, ATrip, MxS, Occ, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B,&
                  HNO, AuxData%IndN_s, AuxData%IndN_s, map, AuxData%IndX_s, IndAux, IndMod, AuxData%NDim_s, AuxData%NDim_s, IGem, AuxInd, &
                  NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, 0, Flags, 1)
            print*, 'Czas na A1: ', clock_readwall(timer)
            do i = 1, 3
                  block_pair(1, i) = 2
            end do
            do i = 4, 6
                  block_pair(1, i) = 4
            end do

            block_pair(1, 7) = 3
            block_pair(1, 8) = 3

            block_pair(2, 1) = 1
            block_pair(2, 2) = 5
            block_pair(2, 3) = 3

            block_pair(2, 4) = 1
            block_pair(2, 5) = 5
            block_pair(2, 6) = 3

            block_pair(2, 7) = 1
            block_pair(2, 8) = 5
            block_pair(1, 9) = 4
            block_pair(2, 9) = 5


            ECorr_AC0 = zero
            do i = 1, 9
                  print*, 'teraz robie blok-lacza sing', i
                  il  = block_pair(1,i)
                  ir = block_pair(2,i)
                  call clock_start(timer)
                  print*, ''       
                  print*, bbl_name(i), eigsPA(il)%dim, eigsPA(ir)%dim
                  print*, ''
                  print*, 'wycinek macierzy A1 singlet'
                  call calc_block_fofo(i, eigsPA(il), eigsPA(ir), ASing, pos, pos_block_s, Occ, AuxData%IndN_s, map, IndAux, IndMod, &
                        NBasis, NA, NI, NV, 0, Flags, E_contr_s(i), IntJFile, IntKFile, IntVV)
                  print*, 'Czas na sing ', bbl_name(i), ' ', clock_readwall(timer)
                  write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ')_s=', E_contr_s(i)
                  ECorr_AC0 = ECorr_AC0 + E_contr_s(i)
                  il = 5 + il
                  ir = 5 + ir

                  call clock_start(timer)

                  print*,	'wycinek macierzy A1 triplet'
                  call calc_block_fofo(i, eigsPA(il), eigsPA(ir), ATrip, pos, pos_block_t, Occ, AuxData%IndN_s, map, IndAux, IndMod, &
                        NBasis, NA, NI, NV, 1, Flags, E_contr_t(i), IntJFile, IntKFile, IntVV)
                  print*, 'Czas na trip ', bbl_name(i), ' ',clock_readwall(timer)

                  write(*,'(A8, A4, A3, F20.15)') 'E_contr(', bbl_name(i), ')_t=', E_contr_t(i)

                  if (i <9)then
                        write(*,'(A8, A4, A2, F20.15)') 'E_contr(', bbl_name(i), ')=', E_contr_s(i) + E_contr_t(i)
                  else
                        write(*,'(A8, A5, A2, F20.15)') 'E_contr(', bbl_name(i), ')=', E_contr_s(i) + E_contr_t(i)
                  end if
                  print*, i, bbl_name(i)
                  ECorr_AC0 = ECorr_AC0 + E_contr_t(i)
            end do

            print*, ''
            print*, 'RDSC ACPP0 CONTR ECORR ', ECorr_AC0
            print*, ''


      end subroutine ACPP0_fast_fofo


      subroutine eigs0_block_fofo(x, ASing, ATrip, pos, pos_block_s, pos_block_t, Eig, Eigvec, v_plus, Occ, &
            IndN, IndAux, NDim, &
            SpinSymm, Flags)

            integer, intent(in) :: x
            double precision, dimension(:,:), intent(in) :: ASing, ATrip
            integer, dimension(:,:), intent(in) :: pos, pos_block_s, pos_block_t
            double precision, dimension(:), intent(inout) :: Eig
            double precision, dimension(:,:), intent(inout) :: Eigvec
            integer, dimension(:), intent(inout) :: v_plus
            double precision, dimension(:),    intent(in)   :: Occ
            integer, dimension(:,:), intent(in)  :: IndN
            integer, dimension(:), intent(in)    :: IndAux
            integer, intent(in) :: NDim
            integer, intent(in) :: SpinSymm
            type (tclock) :: timer, timer0

            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:,:), allocatable :: MxA, MxS
            double precision, dimension(:), allocatable :: Eig_i
            integer :: i

            if (NDim .gt.0)then
                  allocate(MxA(NDim, NDim))
                  allocate(MxS(NDim, NDim))
                  allocate(Eig_i(NDim))


                  call clock_start(timer)
                  if (SpinSymm==0)then
                        call cut_A_piece(ASing, pos, MxA, MxS, NDim, NDim, IndN, IndN, Occ, .true.)
                  else if (SpinSymm==1)then
                        call cut_A_piece(ATrip, pos, MxA, MxS, NDim, NDim, IndN, IndN, Occ, .true.)
                  end if
                  print*, 'Czas na cut ', clock_readwall(timer)
                  if (x == IOO.or.x == IVV.or. x== IOOT.or.x==IVVT)then
                        Eigvec = zero
                        do i = 1, NDim
                              Eig(i) = MxA(i, i)
                        end do
                        if (x == IOO.or.x==IOOT)then
                              v_plus = 0
                        else
                              v_plus = 1
                        end if
                  else
                        call clock_start(timer)

                        call nonsymmetric_eigenproblem_block(Eig, Eig_i, Eigvec, MxA, MxS, v_plus, IndN, IndAux, 1)
                        print*, 'Czas na diag ', clock_readwall(timer)

                  end if

                  deallocate(MxA)
                  deallocate(MxS)
                  deallocate(Eig_i)
            end if
            ! print*, 'Czas diagonalizacji sing: ', clock_readwall(timer)


      end subroutine eigs0_block_fofo

      subroutine cut_A_piece(big_A, pos, small_A, small_S, NDim1, NDim2, IndN1, IndN2, Occ, norm)

            double precision, dimension(:,:), intent(in) :: big_A
            double precision, dimension(:,:), intent(out) :: small_A, small_S
            integer, dimension(:, :), intent(in) :: pos
            integer, intent(in) :: NDim1, NDim2
            integer, dimension(:,:), intent(in) :: IndN1, IndN2
            double precision, dimension(:), intent(in) :: Occ
            logical, intent(in) :: norm

            integer :: p, q, r, s
            integer :: i, j, big_i, big_j
            double precision :: denom

            small_A = zero

            do i = 1, NDim1
                  r = IndN1(1, i)
                  s = IndN1(2, i)
                  denom = one
                  if (norm) then
                        denom = (one-Occ(r)-Occ(s))
                        small_S(i,i) = denom
                  end if
                  do j = 1, NDim2
                        p = IndN2(1, j)
                        q = IndN2(2, j)

                        big_i = pos(r, s)
                        big_j = pos(p, q)

                        small_A(i, j) = big_A(big_i, big_j)/denom
                  end do
            end do


      end subroutine cut_A_piece



      subroutine calc_block_fofo(x, ePAl, ePAr, MxA1_big, pos, pos_block, Occ, IndN, map, IndAux, IndMod, &
            NBasis, NA, NI, NV, SpinSymm, Flags, E_contr, IntJFile, IntKFile, IntVV)

            integer, intent(in) :: x
            type(eigsBlockParams), intent(in) :: ePAl, ePAr
            double precision, dimension(:),    intent(in)   :: Occ
            double precision, dimension(:, :), intent(in) :: MxA1_big
            integer, dimension(:,:), intent(in) :: pos, pos_block
            integer, dimension(:,:), intent(in) :: IndN
            integer, dimension(:), intent(in) :: map
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: SpinSymm
            type(FlagsData), intent(in) :: Flags
            double precision, intent(out) :: E_contr
            character(*)                 :: IntJFile,IntKFile,IntVV
            integer, parameter :: vvoo=1, vvao=2, vvaa=3, vaoo=4                                    
            integer, parameter :: vaao=5, vaaa=6, aaoo=7, aaao=8
            integer, parameter :: vaao2 = 9
            integer :: p, q, r, s, i, j
            integer :: np,km
            double precision :: Npqrs, Aux1, Aux2
            double precision, dimension(:,:), allocatable :: MxA1, MxS
            double precision, dimension(:), allocatable :: tempx, tempx2, tempx3
            double precision :: gr, gr2, gr3, sa, tr1, tr2
            double precision, dimension(:,:), allocatable :: IMxA1, IMxA2

            integer :: z
            integer, external :: NAddrRDM
            integer, external :: NAddr3
            type (tclock) :: timer, timer0, timerall
            real :: start_time, end_time, total_io_time
            integer :: kk
            integer :: NIA
            integer :: iunit1
            integer :: ipq, irs, kl, k, l, jrs
            double precision, dimension(:,:), allocatable :: ints
            double precision :: val
            external :: dgemv

            E_contr = zero
            gr = zero
            gr2 = zero
            gr3 = zero
            sa = zero
            tr1 = zero

            NIA = NI+NA

            allocate(MxA1(ePAl%dim, ePAr%dim))
            allocate(MxS(1,1))
            allocate(ints(NBasis, Nbasis))


            allocate(IMxA2(ePAr%dim, ePAl%dim))

            ! Liczymy wybrany wycinek macierzy A1.
            call clock_start(timer)

            kk = size(MxA1_big, dim=1)

            call cut_A_piece(MxA1_big, pos, MxA1, MxS, ePAl%dim, ePAr%dim, ePAl%IndN, ePAr%IndN, Occ, .false.)

            print*, 'wycinek macierzy MxA1', ePAl%dim, ePAr%dim, SpinSymm
            ! do i = 1, ePAl%dim
            !    do j = 1, ePAr%dim
            !       if (abs(MxA1(i, j)).gt.1.d-5)then
            !          write(*, '(2I5, F20.15)') i, j, MxA1(i,j)
            !       end if
            !    end do
            ! end do

            ! MxA1 = MxA1- MxA

            print*, 'Czas na wybrany wycinek macierzy A1: ', clock_readwall(timer)
            deallocate(MxS)
            allocate(IMxA1(ePAl%dim, ePAr%dim))



            do i = 1, ePAl%dim
                  do j = 1, ePAr%dim
                        IMxA2(j, i) = MxA1(i, j)
                  end do
            end do

            allocate(tempx(ePAl%dim))
            allocate(tempx2(ePAr%dim))


            if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then

                  call clock_start(timer)
                  IMxA1 = zero
                  do np = 1, ePAl%dim
                        do km = 1, ePAr%dim
                              if (ePAl%v_plus(np)==1.and.ePAr%v_plus(km)==0)then
                                    call real_av_x(tempx, MxA1, ePAl%dim,  ePAr%Eigvec(:,km), ePAl%dim, ePAr%dim, one, zero)
                                    call real_vw_x(Aux2, ePAl%Eigvec(:,np), tempx, ePAl%dim)
                                    Aux2 = Aux2 /(ePAl%Eig(np)-ePAr%Eig(km))

                                    do j = 1, ePAr%dim
                                          IMxA1(np, j) = IMxA1(np, j) + Aux2 * ePAr%Eigvec(j, km)
                                    end do
                              end if
                        end do
                  end do
                  print*, 'Czas na pierwsze petle: ', clock_readwall(timer)
            end if


            call clock_start(timer)

            select case(x)

            case(VVAO, VVAA)

                  allocate(tempx3(NV))

                  open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
                        access='DIRECT',recl=8*NBasis*NIA)


                  do k = NIA + 1, NBasis

                        do km = 1, ePAr%dim
                              if (ePAr%v_plus(km)==0)then

                                    p = k
                                    do q = NIA + 1, p

                                          ipq = pos_block(p, q)
                                          if (ipq>0)then
                                                call real_vw_x(Aux2, IMxA2(:, ipq), ePAr%Eigvec(:,km), ePAr%dim)
                                                tempx3(q-NIA) = Aux2
                                          end if
                                    end do


                                    if (x==VVAO)then
                                          do l = 1, NI
                                                kl = k + NBasis * (l-1)
                                                !                      call cpu_time(start_time)
                                                read(iunit1, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NIA)
                                                !                      call cpu_time(end_time)

!                                                total_io_time = total_io_time + (end_time - start_time)
                                                ints(:,NIA+1:NBasis) = 0

                                                s = l
                                                do q = NIA + 1, p
                                                      do r = NI+1, NIA

                                                            val = ints(q, r)

                                                            Aux1 = one
                                                            if (p==q)then
                                                                  Aux1 = Aux1 * sqrt(frac12)
                                                            end if
                                                            if (SpinSymm == 1)then
                                                                  Aux1 = -Three * Aux1
                                                            end if


                                                            ipq = pos_block(p, q)
                                                            jrs = pos_block(r, s)

                                                            if (ipq>0.and.jrs>0)then

                                                                  Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))

                                                                  E_contr = E_contr - tempx3(q-NIA) * Npqrs * Aux1 * val / (ePAl%Eig(ipq)-ePAr%Eig(km)) *ePAr%Eigvec(jrs, km)
                                                            end if
                                                      end do
                                                end do
                                          end do
                                    end if

                                    do l = NI+1, NIA

                                          kl = k + NBasis * (l-1)
                                          read(iunit1, rec=kl) ((ints(i, j), i = 1, NBasis), j = 1, NIA)                   	
                                          ints(:,NIA+1:NBasis) = 0

                                          r = l
                                          do q = NIA + 1, p
                                                if (x==VVAO)then
                                                      do s = 1, NI
                                                            val = ints(q, s)

                                                            Aux1 = one
                                                            if (p==q)then
                                                                  Aux1 = Aux1 * sqrt(frac12)
                                                            end if
                                                            if (SpinSymm == 1)then
                                                                  Aux1 = Three * Aux1
                                                            end if
                                                            ipq = pos_block(p, q)
                                                            jrs = pos_block(r, s)

                                                            if (ipq>0.and.jrs>0)then
                                                                  Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))

                                                                  E_contr = E_contr - tempx3(q-NIA) * Npqrs * Aux1 * val / (ePAl%Eig(ipq)-ePAr%Eig(km)) *ePAr%Eigvec(jrs, km)
                                                                  !write(*, '(7I5, 6F20.15)') ipq, km, jrs, p, q, r, s, tempx3(q-NIA), Npqrs, Aux1*val, val, (ePAl%Eig(ipq)-ePAr%Eig(km)) *ePAr%Eigvec(jrs, km), E_contr
                                                            end if
                                                      end do
                                                else if (x==VVAA)then
                                                      do s = NI+1, NIA
                                                            val = ints(q, s)

                                                            Aux1 = one
                                                            if (p==q)then
                                                                  Aux1 = Aux1 * sqrt(frac12)
                                                            end if
                                                            if (r==s)then
                                                                  Aux1 = Aux1 * sqrt(frac12)
                                                            end if

                                                            if (SpinSymm == 1)then
                                                                  Aux1 = Three * Aux1
                                                            end if
                                                            ipq = pos_block(p, q)
                                                            jrs = pos_block(r, s)

                                                            if (ipq>0.and.jrs>0)then
                                                                  Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))

                                                                  E_contr = E_contr - tempx3(q-NIA) * Npqrs * Aux1 * val / (ePAl%Eig(ipq)-ePAr%Eig(km)) *ePAr%Eigvec(jrs, km)
                                                                  !                                  write(*, '(7I5, 6F20.15)') ipq, km, jrs, p, q, r, s, tempx3(q-NIA), Npqrs, Aux1*val, val, (ePAl%Eig(ipq)-ePAr%Eig(km)) *ePAr%Eigvec(jrs, km), E_contr                        
                                                            end if

                                                      end do
                                                end if
                                          end do
                                    end do
                              end if
                        end do
                  end do

                  !print *, "Całkowity czas spędzony na operacjach I/O: ", total_io_time, " sekund"
                  close(iunit1)
            end select

            ! call clock_start(timer)
            ! do i = 1, ePAl%dim
            !    do km = 1, ePAr%dim
            !       if (ePAr%v_plus(km)==0)then
            !          call real_vw_x(Aux2, IMxA2(:, i), ePAr%Eigvec(:,km), ePAr%dim)
            !          do j = 1, ePAr%dim
            !             p = ePAl%IndN(1, i)
            !             q = ePAl%IndN(2, i)

            !             r = ePAr%IndN(1, j)
            !             s = ePAr%IndN(2, j)

            !             if (SpinSymm == 0)then
            !                if (x == VAAO)then
            !                   Aux1 = TwoNO(NAddr3(r,p,s,q))
            !                else if (x==VAAO2) then
            !                   Aux1 = TwoNO(NAddr3(r,q,s,p))
            !                else
            !                   Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
            !                end if
            !                if (p==q)then
            !                   Aux1 = Aux1 * sqrt(frac12)
            !                end if
            !                if (r==s) then
            !                   Aux1 = Aux1 * sqrt(frac12)
            !                end if
            !             else
            !                if (x == VAAO)then
            !                   Aux1 = TwoNO(NAddr3(r,p,s,q))
            !                else if    (x==VAAO2) then
            !                   Aux1 = -TwoNO(NAddr3(r,q,s,p))
            !                else
            !                   Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
            !                end if
            !                Aux1 = Three*Aux1
            !             end if

            !             Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))

            !             E_contr = E_contr - Aux2 * Npqrs * Aux1 / (ePAl%Eig(i)-ePAr%Eig(km)) *ePAr%Eigvec(j, km)
            !          end do
            !       end if
            !    end do
            ! end do
            ! print*, 'Czas na drugie petle: ', clock_readwall(timer)

            !  case(VAOO, AAOO)
            !    call clock_start(timer)

            !     do j = 1, ePAr%dim
            !        r = ePAr%IndN(1, j)
            !        s = ePAr%IndN(2, j)
            !        do np = 1, ePAl%dim
            !           if (ePAl%v_plus(np)==1)then
            !              call real_vw_x(Aux2, MxA1(:, j), ePAl%Eigvec(:,np), ePAl%dim)

            !           do i = 1, ePAl%dim
            !              p = ePAl%IndN(1, i)
            !              q = ePAl%IndN(2, i)
            !              if (SpinSymm == 0)then
            !                 if (x == VAAO)then
            !                    Aux1 = TwoNO(NAddr3(r,p,s,q))
            !                 else if (x==VAAO2) then
            !                    Aux1 = TwoNO(NAddr3(r,q,s,p))
            !                 else
            !                    Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
            !                 end if
            !                 if (p==q)then
            !                    Aux1 = Aux1 * sqrt(frac12)
            !                 end if
            !                 if (r==s) then
            !                    Aux1 = Aux1 * sqrt(frac12)
            !                 end if
            !              else
            !                 if (x == VAAO)then
            !                    Aux1 = TwoNO(NAddr3(r,p,s,q))
            !                 else if (x==VAAO2) then
            !                    Aux1 = -TwoNO(NAddr3(r,q,s,p))
            !                 else
            !                    Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
            !                 end if
            !                 Aux1 = Three*Aux1
            !              end if

            !              Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))
            !              E_contr = E_contr - Aux2 * Npqrs * Aux1 / (ePAl%Eig(np)-ePAr%Eig(j)) *ePAl%Eigvec(i,np)
            !           end do
            !        end if
            !        end do
            !     end do

            !     print*, 'Czas na trzecie petle: ', clock_readwall(timer)

            !  case default
            !     call clock_start(timer)

            !     i_rowloops: do i = 1, ePAl%dim
            !        j_colloops: do j = 1, ePAr%dim
            !           p = ePAl%IndN(1, i)
            !           q = ePAl%IndN(2, i)

            !           r = ePAr%IndN(1, j)
            !           s = ePAr%IndN(2, j)

            !           if (SpinSymm == 0)then
            !              if (x == VAAO)then
            !                 Aux1 = TwoNO(NAddr3(r,p,s,q))
            !              else if (x==VAAO2) then
            !                 Aux1 = TwoNO(NAddr3(r,q,s,p))
            !              else
            !                 Aux1 = (TwoNO(NAddr3(r,p,s,q)) +  TwoNO(NAddr3(r,q,s,p)))
            !              end if
            !      	if (p==q)then
            !                 Aux1 = Aux1 * sqrt(frac12)
            !              end if
            !              if (r==s) then
            !                 Aux1 = Aux1 * sqrt(frac12)
            !              end if
            !           else
            !              if (x == VAAO)then
            !                 Aux1 = TwoNO(NAddr3(r,p,s,q))
            !              else if (x==VAAO2) then
            !                 Aux1 = -TwoNO(NAddr3(r,q,s,p))
            !      	else
            !                 Aux1 = (TwoNO(NAddr3(r,p,s,q))-  TwoNO(NAddr3(r,q,s,p)))
            !              end if
            !      	Aux1 = Three*Aux1
            !           end if

            !           Npqrs =(one-Occ(p)-Occ(q)) * (one-Occ(r)-Occ(s))

            !           if (x==VVOO)then
            !              E_contr = E_contr - MxA1(i, j) * Npqrs * Aux1 / (ePAl%Eig(i)-ePAr%Eig(j))
            !           else if (x==VAAO.or.x==VAAA.or.x==AAAO.or.x==VAAO2)then
            !                 do np = 1, ePAl%dim
            !                    E_contr = E_contr - Npqrs * Aux1 * ePAl%Eigvec(i,np) * IMxA1(np, j)
            !                 end do
            !              end if

            !        end do j_colloops
            !     end do i_rowloops
            !      print*, 'Czas na czwarte petle: ', clock_readwall(timer)

            !  end select
            !  !       !omp end parallel do                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

            !  print*, ''
            !  print*, 'Czas na timeall  petle: ', clock_readwall(timerall)
            !  print*, ''

            !  deallocate(MxA1)
            !  deallocate(tempx)



      end subroutine calc_block_fofo




      subroutine updateIndices_fofo(p, q, pos_block_s, pos_block_t, IndN_s, IndN_t, ind_s, ind_t)
            integer, intent(in) :: p, q
            integer, dimension(:,:), intent(out) :: pos_block_s, pos_block_t
            integer, dimension(:,:), intent(inout) :: IndN_s, IndN_t
            integer, intent(inout) :: ind_s, ind_t


            if (p.ne.q) then
                  pos_block_t(p,q) = ind_t
                  IndN_t(1,ind_t) = p
                  IndN_t(2,ind_t) = q
                  ind_t = ind_t + 1
            endif
            pos_block_s(p,q) = ind_s
            IndN_s(1,ind_s) = p
            IndN_s(2,ind_s) = q
            ind_s = ind_s + 1

      end subroutine updateIndices_fofo


      subroutine nonsymmetric_eigenproblem_block_fofo(wr, wi, vr, A, S, v_plus, IndN, IndAux, AC_TYPE)
            double precision, dimension(:), intent(out)      :: wr
            double precision, dimension(:), intent(out)      :: wi
            double precision, dimension(:, :), allocatable   :: vl
            double precision, dimension(:, :), intent(out)   :: vr
            double precision, dimension(:, :), intent(inout) :: A, S
            integer, dimension(:), intent(out) :: v_plus
            integer, dimension(:,:), intent(in) :: indN
            integer, dimension(:),intent(in) :: IndAux
            integer, intent(in) :: AC_TYPE
            double precision :: shift, max_hh, min_pp
            double precision, dimension(:),allocatable :: S_diag

            double precision, dimension(1) :: work0
            double precision, dimension(:), allocatable :: work, tempx
            integer :: lwork, info
            integer :: n, i, j, k
            double precision :: norm, norm_l, norm_rl
            double precision :: dd, ddm, dd_comp
            type (tclock) :: timer
            double precision, parameter :: tol = 1.d-4
            double precision, dimension(:), allocatable :: wr_plus, wr_minus

            integer, dimension(:), allocatable :: dy_plus, dy_minus
            integer :: count, countj_plus, countj_minus, ss
            logical :: flag_change

            double precision, dimension(:,:), allocatable :: miniA
            double precision, dimension(:,:), allocatable :: vl2, vr2
            double precision, dimension(:), allocatable :: wr2, wi2
            integer :: ii, jj, n2, r, s2, p,q


            external :: dgeev

            n = size(A, dim=1)
            allocate(dy_plus(n))
            allocate(dy_minus(n))
            allocate(wr_plus(n))
            allocate(wr_minus(n))
            allocate(S_diag(n))

            allocate(vl(1,1))
            vl = zero
            vr = zero

            do i = 1, n
                  S_diag(i) = S(i, i)
            end do

            ! diagonalizacja
            call clock_start(timer)    
            lwork = -1
            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work0, lwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))
            call dgeev("n", "V", n, A, n, wr, wi, vl, n, vr, n, work, lwork, info)
            print*, 'info', info
            if (info /= 0) then
                  print*, "Nonsymmetric matrix eigendecompositino failed with info="
                  ! call msg("Nonsymmetric matrix eigendecompositino failed with info=" // str(info), MSG_ERROR)
                  error stop
            end if
            print*, 'Czas diagonalizacji real: ', clock_readwall(timer)

            allocate(tempx(n))
            call clock_start(timer)
            v_plus = 0
            ! print*, 'eigenval na poczatku'
            ! do i = 1, n
            !    print*, wr(i)
            ! end do
            ! print*, ''
            ! print*, 'eigvec na poczatku'
            ! do i = 1, n
            !    print*, 'wetor przed', i
            !    do j = 1, n
            !       print*, vr(j, i)
            !    end do

            ! end do

            ! do i = 1, n
            !    do j = 1, n
            !       call ddot_norm(vr(:, i),S,  vr(:, j), n, dd)
            !       print*, 'normerr', i, j, dd
            !    end do
            ! end do


            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        if (wi(i).gt.0.d+0)then
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)
                              call ddot_norm(vr(:, i+1),S,  vr(:, i+1), n, dd_comp)
                              print*, '>0 dd', i, i+1, dd
                              print*, '>0 dd_comp', dd_comp
                        else
                              call ddot_norm(vr(:, i-1),S,  vr(:, i-1), n, dd)
                              call ddot_norm(vr(:, i),S,  vr(:, i), n, dd_comp)
                              print*, '<0 dd', i-1, i, dd
                              print*, '<0 dd_comp',	dd_comp

                        end if
                        dd = dd + dd_comp
                        v_plus(i)=2
                  else
                        call ddot_norm(vr(:, i),S,  vr(:, i), n, dd)
                        if (dd.gt.zero)then
                              v_plus(i) = 1
                        end if

                  end if
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'i-norm', i, dd
                  end if
            end do

            if (AC_TYPE == PPAC)then
                  print*, 'AC_TYPE: PP'
            else if (AC_TYPE == HHAC) then
                  print*, 'AC_TYPE: HH'
            end if
            !    print*, 'v_plus', v_plus

            print*, 'Czas r1: ', clock_readwall(timer)
            call clock_start(timer)
            print*, 'eigenvalues of singlets'
            do i = 1, n
                  if(abs(wi(i)).gt.1.d-8)then
                        print*, 'COMPLEX EIGENVALUES', i, wr(i), wi(i)
                        !         stop
                  end if
            end do

            wr_plus = wr
            wr_minus = wr
            do i = 1, n
                  dy_plus(i)=i
            end do
            call dsort(wr_plus, dy_plus, n)


            shift = zero
            do i = 1, n
                  if (v_plus(dy_plus(i))== 0 )then
                        max_hh = wr_plus(i)
                  end if
                  if (v_plus(dy_plus(i))== 1 )then
                        min_pp = wr_plus(i)
                        if (i == 1) then
                              shift = min_pp/two
                        else
                              shift = (abs(min_pp) + abs(max_hh))/two
                              print*, 'max_hh', max_hh
                        end if
                        print*, 'min_pp' , wr_plus(i)
                        print*, 'shift', shift
                        exit
                  else

                  end if
            end do


            wr_plus = zero
            wr_minus = zero
            dy_plus = 0
            dy_minus = 0
            countj_plus=1
            countj_minus=1

            do i = 1, n
                  if (v_plus(i) == 1) then
                        wr_plus(countj_plus) = wr(i)
                        dy_plus(countj_plus) = i
                        countj_plus = countj_plus + 1
                  else if (v_plus(i) == 0) then
                        wr_minus(countj_minus) = wr(i)
                        dy_minus(countj_minus) = i
                        countj_minus = countj_minus + 1

                  end if
                  !     dy(i) = i
            end do
            countj_plus = countj_plus-1
            countj_minus = countj_minus-1


            ! print*, 'wr_plus przed', wr_plus, dy_plus
            call dsort(wr_plus(1:countj_plus), dy_plus(1:countj_plus), countj_plus)
            ! print*, 'wr_plus po', wr_plus, dy_plus

            ! print*, 'wr_minus przed', wr_minus, dy_minus
            call dsort(wr_minus(1:countj_minus), dy_minus(1:countj_minus), countj_minus)
            ! print*, 'wr_minus po', wr_minus, dy_minus

            print*, 'countj_plus jest',countj_plus
            print*, 'countj_minus jest',countj_minus

            ! do i = 1, n
            !    print*, 'wetor przedtuz', i
            !    do j = 1, n
            !       print*, vr(j, i)
            !    end do

            ! end do

            call orthogonalize_degen(n, countj_plus, wr_plus, vr, S_diag, dy_plus, 1)
            ! do i = 1, n
            !    print*, 'wetor po1', i
            !    do j = 1, n
            !       print*, vr(j, i)
            !    end do

            ! end do

            call orthogonalize_degen(n, countj_minus, wr_minus, vr, S_diag, dy_minus, 0)
            ! do i = 1, n
            !    print*, 'wetor po2', i
            !    do j = 1, n
            !       print*, vr(j, i)
            !    end do

            ! end do


            !     do i = 1, n
            !    do j = 1, n
            !       call ddot_norm(vr(:, i),S,  vr(:, j), n, dd)
            !       print*, 'normerr2', i, j, dd
            !    end do
            ! end do



      end subroutine nonsymmetric_eigenproblem_block_fofo


      subroutine orthogonalize_degen_fofo(n, countj, wr, vr, S_diag, dy, x)

            integer, intent(in) :: n
            double precision, dimension(:), intent(in) :: wr
            double precision, dimension(:,:), intent(inout) :: vr
            double precision, dimension(:), intent(in) :: S_diag
            integer, dimension(:), intent(in) :: dy
            integer, intent(in) :: x
            integer, dimension(:), allocatable :: StartIdx, EndIdx
            integer, intent(in) :: countj
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
      end subroutine orthogonalize_degen_fofo



      subroutine PPERPA_block_fofo(ASing, ATrip, MxS, Occ, AuxII, AuxAA, BuxII, BuxAA, Aux3A, Aux3B,&
            HNO, IndN1, IndN2, map, IndX, IndAux, IndMod, NDim1,NDim2, IGem, AuxInd, &
            NBasis, NA, NI, NV, NInte1, NInte2,  ACAlpha, SpinSymm, Flags, dos)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip, MxS
            double precision, dimension(:),    intent(in)   :: Occ
            double precision, dimension(:,:), intent(in) :: AuxII, AuxAA, BuxII, BuxAA,Aux3A, Aux3B
            double precision, dimension(:,:), intent(inout) :: HNO
            integer, dimension(:,:), intent(in)  :: IndN1, IndN2
            integer, dimension(:), intent(in) :: map, IndX
            integer, dimension(:), intent(in)    :: IndAux, IndMod
            integer, intent(in) :: NDim1, NDim2, NBasis, NA, NI, NV
            integer, dimension(:), intent(in) :: IGem
            integer, dimension(:,:), intent(in) :: AuxInd
            integer, intent(in) :: NInte1, NInte2
            integer, intent(in) :: SpinSymm
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            integer, intent(in) :: dos
            double precision, dimension(:,:), allocatable :: MxT

            integer :: NRDM2, NRDM2Act, NOc
            double precision :: Arspq, Arsqp, dm, Arssum
            type (tclock) :: timer, timer1, timer2, timer3, timer4

            !
            ! External procedures
            !    

            integer :: p, q, r, s, t, u, v, pp, qq
            integer :: i ,j, k, l, ij, a, b, ab, kl
            integer :: ii, pq, rs, pr, qs, ps, qr
            double precision :: NumF, NumH
            integer :: c, d
            double precision :: temp
            integer :: jj
            logical :: ism
            double precision :: val, vals, valt, valh



            NRDM2 = NBasis**2*(NBasis**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NOc = NI + NA

                MxS = zero

                call clock_start(timer)

                do q = 1, NBasis
                      do s = 1, Nbasis
                            if (IGem(q)==Igem(s))then
                                  if(IGem(q)==1.or.IGem(q)==3)then
                                        HNO(q,s) = HNO(q,s) + (One/ACAlpha -One) * ( Two * BuxAA(q,s)- AuxAA(q,s))
                                  end if                            
                                  HNO(q,s) = HNO(q,s) + (One/ACAlpha-one) * (Two * BuxII(q, s)  - AuxII(q,s))
                            end if
                      end do
                end do 
                 print*, 'mandar'
                 print*, 'hno',  clock_readwall(timer)
                 call clock_start(timer)
!               !$omp parallel do collapse(2)&
!               !$omp default(shared) &
!               !$omp prIVATe(i, j) &
!               !$omp prIVATe(p, q, pq, r, s, rs, val, vals, valt, NumF)


                  j_colloops: do j = 1, NDim2
                        i_rowloops: do i = 1, NDim1

                            r = IndN1(1, i)
                            s = IndN1(2, i)
                            rs = IndX(i)

                            p = IndN2(1, j)
                            q = IndN2(2, j)
                            pq = IndX(j)


                            val = zero
                            vals = zero
                            

                            if (p==r) then


                                  val = val + two * (BuxII(q, s) + BuxAA(q, s)) !P2a-1
                                  val = val - AuxII(q, s) - AuxAA(q, s)   !P2a-1
                                  val = val + (one- Occ(p)-frac12 * Occ(q) - frac12*Occ(s)) * HNO(q,s) !P1
                                  if (IGem(q)==2) val = val - frac12 * Aux3A(q, s)
                                  if (IGem(s)==2) val = val - frac12 * Aux3B(s, q)
                                  ! val = val - frac12 * (WMAT(q,s) + WMAT(s, q))  !P4a2 + P5a2
                                  
                                  if (Igem(q)==1) val =	val  - (BuxII(q,s) + BuxAA(q, s)) + frac12 * (AuxII(q,s) + AuxAA(q, s))
                                  if (Igem(q)==2) val = val  - Occ(q) * (BuxII(q, s) - frac12 * AuxII(q,s))

                                  if (Igem(s)==1) val = val  - (BuxII(q,s) + BuxAA(q, s))	+ frac12 * (AuxII(q,s) + AuxAA(q, s))
                                  if (Igem(s)==2) val = val  - Occ(s) * (BuxII(q, s) - frac12 * AuxII(q,s))


                                  select case(AuxInd(IGem(p),IGem(r)))
                                  case(1)
                                        val = val + Occ(p)*(AuxII(q,s) + AuxAA(q,s) - two * BuxII(q, s) -two * BuxAA(q, s))  ! p3B-1, P3c-3,4
                                  case(2)
                                        val = val + Occ(p)*(AuxII(q,s)- two * BuxII(q, s))
                                  end select
                            end if

                            if (q==s)then

                                  
                                  val = val + two * (BuxII(p,r) + BuxAA(p,r)) !p2A-2
                                  val = val - AuxII(p,r) - AuxAA(p,r) !p2A-2
                                  val = val+ (one- Occ(q)-frac12 * Occ(p) - frac12*Occ(r)) * HNO(p,r) !p1-2
                                   if (IGem(p)==2) val = val - frac12 * Aux3A(p, r)
                                   if (IGem(r)==2) val = val - frac12 * Aux3B(r, p)

                                  ! val = val - frac12 * (WMAT(p,r) + WMAT(r, p)) !p4A-1 p5A1

                                  if (Igem(p)==1) val = val  - (BuxII(p,r) + BuxAA(p,r)) + frac12 * (AuxII(p,r) + AuxAA(p,r))
                                  if (Igem(p)==2) val = val  - Occ(p) * (BuxII(p,r) - frac12 * AuxII(p,r))

                                  if (Igem(r)==1) val = val  - (BuxII(p,r) + BuxAA(p,r))   + frac12 * (AuxII(p,r) + AuxAA(p,r))
                                  if (Igem(r)==2) val = val  - Occ(r) * (BuxII(p,r) - frac12 * AuxII(p,r))


                                  select case(AuxInd(IGem(q),IGem(s)))
                                  case(1)
                                        val = val + Occ(q)*(AuxII(p,r) + AuxAA(p,r) - two* BuxII(p, r) - two * BuxAA(p, r)) ! p3B-2 p3C-1,2
                                  case(2)
                                        val = val + Occ(q)*(AuxII(p,r)- two *BuxII(p, r))                
                                  end select

                            end if
                            vals = val
                            valt = val
                            val = zero
                            if (p==s)then
                                  
                                  val = val - two	* (BuxII(q,r) + BuxAA(q,r )) !p2a-3
                                  val = val + (AuxII(q,r) + AuxAA(q,r))!P2a3
                                  val = val - (one- Occ(p)-frac12 * Occ(q) - frac12*Occ(r)) * HNO(q,r) !P1-3
                                  if (IGem(r)==2) val = val + frac12 * Aux3B(r, q)
                                  if (IGem(q)==2) val = val + frac12 * Aux3A(q, r)

                                  ! val = val + frac12 * (WMAT(q,r) + WMAT(r, q)) !P4a-3 P5a-3

                                  if (Igem(q)==1) val = val  + (BuxII(q,r) + BuxAA(q,r)) - frac12 * (AuxII(q,r) + AuxAA(q,r))
                                  if (Igem(q)==2) val = val  + Occ(p) * (BuxII(q,r) - frac12 * AuxII(q,r))
                                  
                                  if (Igem(r)==1) val = val  + (BuxII(q,r) + BuxAA(q,r))   - frac12 * (AuxII(q,r) + AuxAA(q,r))
                                  if (Igem(r)==2) val = val  + Occ(r) * (BuxII(q,r) - frac12 * AuxII(q,r))

                                  select case(AuxInd(IGem(p),IGem(s)))
                                  case(1)
                                        val = val - Occ(p)*(AuxII(q,r) + AuxAA(q,r) - two * BuxII(q, r) - two * BuxAA(q, r)) !p3B-3 p3C-5,6
                                  case(2)
                                        val = val - Occ(p)*(AuxII(q,r) - two * BuxII(q, r))
                                  end select
                            end if

                            if (q==r)then


                                  val = val - two	* (BuxII(p,s) + BuxAA(p,s)) !P2a-4
                                  val = val + (AuxII(p,s) + AuxAA(p,s))!P2a-4
                                  val = val - (one- Occ(q)-frac12 * Occ(p) - frac12*Occ(s)) * HNO(p,s) !P1-4
                                  if (IGem(s)==2) val = val + frac12 * Aux3B(s, p)
                                  if (IGem(p)==2) val = val + frac12 * Aux3A(p, s)

                                  !  val = val + frac12 * (WMAT(p,s) + WMAT(s, p)) !P4a3 P5a3

                                  if (Igem(p)==1) val = val  + (BuxII(p,s) + BuxAA(p,s)) - frac12 * (AuxII(p,s) + AuxAA(p,s))
                                  if (Igem(p)==2) val = val  + Occ(p) * (BuxII(p,s) - frac12 * AuxII(p,s))

                                  if (Igem(s)==1) val = val  + (BuxII(p,s) + BuxAA(p,s))   - frac12 * (AuxII(p,s) + AuxAA(p,s))
                                  if (Igem(s)==2) val = val  + Occ(s) * (BuxII(p,s) - frac12 * AuxII(p,s))


                                  select case(AuxInd(IGem(q),IGem(r)))
                                  case(1)
                                        val = val - Occ(q)*(AuxII(p,s) + AuxAA(p,s) - two * BuxII(p, s) - two * BuxAA(p, s)) !p3B-4 p3C-7,8
                                  case(2)
                                        val = val - Occ(q)*(AuxII(p,s) - two * BuxII(p, s))
                                  end select

                            end if

                            vals = vals -val
                            valt = valt + val

                            NumF = one
                            !          NumF = one/(1-Occ(r)-Occ(s))

                            select case(SpinSymm)
                            case(0)
                                  NumH = one
                                  if (p==q.and.r==s)then
                                        NumH = frac12
                                  else if ((p==q.and.r.ne.s).or.(p.ne.q.and.r==s))then
                                        NumH = sqrt(frac12)
                                  end if

                                  ! vals = zero
                                  ! valt = zero
                                  ASing(rs, pq) = NumF * NumH * ASing(rs, pq) + NumF * NumH*vals
                                  ATrip(rs, pq) = NumF * ATrip(rs, pq) + NumF  * valt


                                               if (abs(ASing(rs, pq)).gt.1.d-5)then
                                                     write(*,'(2I5, A5, 4I5,2F30.16)') rs, pq, '  |  ', p, q, r, s, ASing(rs, pq)
                                               end if
                                  !              if (abs(ATrip(rs, pq)).gt.1.d-5)then
                                  !                 if (p.ne.q.and.r.ne.s)then
                                  ! !                   write(*,'(2I5, A5, 4I5, F30.16)') rs, pq, '  |  ', p, q, r, s, ATrip(rs, pq)
                                  !                 end if
                                  !              end if

                                  ! case(1)

                                  !    if (p.ne.q.and.r.ne.s)then
                                  !       Arssum = Arspq + Arsqp
                                  !       MxA(i, j) = Arssum 
                                  !    end if
                            end select

                      end do i_rowloops
                end do j_colloops
!                !$omp end parallel do

                print*, 'petle',  clock_readwall(timer)
      end subroutine PPERPA_block_fofo


      subroutine PPERPA_init_fofo(H_type, ETot, ENuc,  n, XOne, TwoNO, map, IndAux, &
            NBasis, NA, NI, NV, NInte1, NInte2,ACAlpha, Flags, TwoNOA, Ha)
            use math_constants
            Use, intrinsic :: iso_fortran_env, Only : iostat_end
            integer, intent(in) :: H_type
            double precision, intent(inout) :: ETot
            double precision, intent(in) :: ENuc
            double precision, dimension(:),    intent(in)   :: n
            double precision, dimension(:),   intent(in) :: XOne
            double precision, dimension(:),      intent(in) :: TwoNO
            integer, dimension(:), intent(in)    :: IndAux
            integer, dimension(:), intent(out)    :: map
            integer, intent(in) :: NBasis, NA, NI, NV
            integer, intent(in) :: NInte1, NInte2
            double precision, intent(in) :: ACAlpha
            type(FlagsData), intent(in) :: Flags
            double precision, dimension(:), intent(out) :: TwoNOA, Ha
            double precision, dimension(:), allocatable :: HNO
            double precision, dimension(:), allocatable :: R00, R11
            integer :: NRDM2, NRDM2Act, NOc

            !
            ! External procedures
            !    
            integer, external :: NAddrRDM
            integer, external :: NAddr3
            double precision, external :: FRDM2

            integer :: p, q, r, s, t, u, v, ii
            integer :: i ,j, k, l, ij, a, b, ab, kl

            integer :: twoint_dim
            double precision :: temp

            NRDM2 = NBasis**2*(NBasis**2+1)/2
            NRDM2Act = NA**2*(NA**2+1)/2
            NOc = NI + NA

            allocate(HNO(NInte1))
            allocate (R00(NRDM2Act))
            allocate (R11(NRDM2Act))

            ! twoint_dim = size(TwoNO, dim=1)
            ! allocate(TwoNOA(twoint_dim))
            print*, 'allocate(R00(NRDM2Act))', NRDM2Act

            R00 = Zero
            R11 = Zero
            print*, 'lllll', NRDM2Act

            !                                                                                                                                      
            ! Fill one-electron hamiltonian in NO rep.                                                                                                
            !
            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        ab = (max(i, j)*(max(i,j)-1))/2 + min(i, j)
                        HNO(ij) = XOne(ab)
                        ! HNO(ij) = zero
                        ! do a = 1, NBasis
                        !    do b = 1, Nbasis
                        !       ab = (max(a, b)*(max(a,b)-1))/2 + min(a, b)
                        !       HNO(ij) = HNO(ij) + URe(i, a)*URe(j,b)*XOne(ab)
                        !    end do
                        ! end do
                  end do
            end do

            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        temp = HNO(ij)
                        do r = 1, NBasis
                              temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                        end do
                  end do
            end do

            print*, 'Number of inacvive orbitals:', NI
            print*, 'Number of active orbitals:  ', NA
            print*, 'Number of virtual orbitals: ', NV


            map = 0

            k = 1
            do i = 1, NI+NA+NV
                  if(IndAux(i)==1)then
                        map(i) = k
                        k = k+1
                  end if
                  if(IndAux(i)==2)then
                        map(i) = 0
                  end if
            end do

            R00 = Zero
            R11 = Zero
            ! print*, R00
            !     print*, 'czytam rdm2'
            call read_2rdm("rdm2.dat", R00, NA)
            call read_2rdm("rdms2.dat", R11, NA)

            ! print*,'przecz', R00
            ETot = zero
            do i = 1, NBasis
                  ii = (i*(i+1))/2
                  ETot = ETot + two* n(i) * HNO(ii)
            end do

            print*, 'etot1', ETot, NOc
            print*, 'NOc', NOc
            stop
            !    print*, 'map'
            !    print*, map

            do p = 1, NI+NA
                  do q = 1, NI+NA
                        do r = 1, NI+NA
                              do s = 1, NI+NA
                                    !               print*, p, q, r, s
                                    ETot = ETot + FRDM2(p, q, r, s, R00, n, map, NA, NBasis) &
                                          * TwoNO(NAddr3(p, r, q, s))
                              end do
                        end do
                  end do
            end do

            print*, 'testing...', ETot, ETot+Enuc, Enuc
            print*, 'RDCS ETot', ETot+Enuc


            Ha = zero
            if (H_type==pDyall)then
                  print*, 'H_type=Dyall', pDyall
            else if (H_type==pGPF)then
                  print*, 'H_type=GPF', pGPF
            end if

            ij = 0
            do i = 1, Nbasis
                  do j = 1, i
                        ij = ij + 1
                        Ha(ij) = ACAlpha * HNO(ij)

                        if (IndAux(i).eq.IndAux(j))then

                              temp = HNO(ij)
                              do r = 1, NBasis

                                    if (H_type==pDyall)then     	
                                          if (.not.((IndAux(r).eq.IndAux(i)).and.(IndAux(r).eq.1)))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                          end if
                                    else if (H_type==pGPF)then

                                          if (IndAux(r).ne.IndAux(i))then
                                                temp = temp + n(r) * (two*TwoNO(NAddr3(r,r,i,j))-TwoNO(NAddr3(r,i,r,j)))
                                          end if
                                    end if

                              end do
                              Ha(ij) = Ha(ij) + (one-ACAlpha)*temp

                        end if
                  end do
            end do


            TwoNOA = TwoNO
            ij = 0
            do i = 1, NBasis
                  do j = 1, i
                        ij = ij + 1
                        kl = 0
                        do k = 1, Nbasis
                              do l = 1, k
                                    kl=kl+1

                                    if (H_type==pDyall)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l).and.IndAux(i)==1)then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                          end if
                                          ! if (abs(ACALPHA-one).lt.1.d-5)then
                                          !    write(*, '(A20, 4I5, 4F7.4, F20.15)') 'calka', i, j, k, l, n(i), n(j), n(k), n(l), TwoNOA(NAddr3(i, j, k, l))
                                          ! end if
                                    else if (H_type==pGPF)then
                                          if ((IndAux(i)==IndAux(j)).and.(IndAux(i)==IndAux(k)).and.IndAux(i)==IndAux(l))then
                                                TwoNOA(NAddr3(i, j, k, l)) = TwoNO(NAddr3(i, j, k, l))
                                          else
                                                TwoNOA(NAddr3(i, j, k, l)) = ACAlpha * TwoNO(NAddr3(i, j, k, l))
                                                if (abs(TwoNOA(NAddr3(i, j, k, l))).gt.1.d-8)then
                                                end if
                                          end if
                                    end if

                              end do
                        end do
                  end do
            end do

            ! if (abs(ACALPHA-one).lt.1.d-5)then
            !    stop
            ! end if

      end subroutine PPERPA_init_fofo


      subroutine update_ASing1(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints
            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ
            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = kk
            i = ii
            l = ll
            j = jj

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !                   if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                  if (rs == 37.and.pq == 78)then
                        write(*, '(A10, 6I5, 2F20.15, I5)') '1temptemp', rs, pq, k, i, l, j,  val, ASing(rs,pq), lr
                        print*, 'val', val
                        print*, 'occ', one- Occ(k) -Occ(i)-Occ(l)-Occ(j)
                        print*,  AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) 
                        print*, kk, ii, ints(kk,ii)
                  end if

            end if


      end subroutine update_ASing1


      subroutine update_ASing2(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = kk
            i = ii
            l = jj
            j = ll

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                   if (rs == 37.and.pq == 78)then
                         write(*, '(A10, 6I5, 2F20.15, I5)') '2temptemp', rs, pq, k, i, l, j,  val, ASing(rs,pq), lr
                         print*, 'val', val
                         print*, 'occ', one- Occ(k) -Occ(i)-Occ(l)-Occ(j)
                         print*,  AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j))
                         print*, kk, ii, ints(kk,ii)
                   end if

            end if

      end subroutine update_ASing2

      subroutine update_ASing3(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = ii
            i = kk
            l = ll
            j = jj

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                      write(*, '(A10, 6I5, 3F20.15, I5)') '3temptemp', rs, pq, k, i, l, j,  val,  ASing(rs,pq), ints(kk,ii), lr
                  end if

            end if


      end subroutine update_ASing3

      subroutine update_ASing4(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = ii
            i = kk
            l = jj
            j = ll

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                     write(*, '(A10, 6I5, F20.15, I5)') '4temptemp', rs, pq, k, i, l, j,  val, lr
                  end if

            end if


      end subroutine update_ASing4

      subroutine update_ASing5(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = ll
            i = jj
            l = kk
            j = ii

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                      write(*, '(A10, 6I5, F20.15, I5)') '5temptemp', rs, pq, k, i, l, j,  val, lr
                  end if

            end if

      end subroutine update_ASing5

      subroutine update_ASing6(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = ll
            i = jj
            l = ii
            j = kk

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                      write(*, '(A10, 6I5, 3F20.15, I5)') '6temptemp', rs, pq, k, i, l, j,  val,  ASing(rs,pq), ints(kk,ii), lr
                  end if

            end if

      end subroutine update_ASing6

      subroutine update_ASing7(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = jj
            i = ll
            l = kk
            j = ii

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  !  if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                     write(*, '(A10, 6I5, F20.15, I5)') '7temptemp', rs, pq, k, i, l, j,  val, lr
                  end if

            end if

      end subroutine update_ASing7

      subroutine update_ASing8(kk, ii, ll, jj, pos, ASing, ATrip, ints, AuxCoeff, lr, Occ)

            integer, intent(in) :: kk, ii, ll, jj
            integer, dimension(:,:), intent(in) :: pos
            double precision, dimension(:,:,:,:), intent(in) :: AuxCoeff
            double precision, dimension(:,:), intent(inout) :: ASing, ATrip
            double precision, dimension(:,:), intent(in) :: ints

            integer, intent(in) :: lr
            double precision, dimension(:), intent(in) :: Occ

            ! lr is integer that is 0 if we update A_rspq, and 1 if A_rsqp
            double precision :: val, s
            integer :: pq, rs, i, j, k, l

            k = jj
            i = ll
            l = ii
            j = kk

            val = AuxCoeff(IGem(k),IGem(i),IGem(l),IGem(j)) * ints(kk,ii)
            if (lr==0)then
                  rs = pos(i, j)
                  s = one
            else if (lr==1)then
                  s = -one
                  rs = pos(j, i)
            end if

            pq = pos(k,l)

            if (pq > 0 .and. rs > 0) then
                  ASing(rs,pq) = ASing(rs,pq) + val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ATrip(rs,pq) = ATrip(rs,pq) + s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))
                  ! if (abs(s * val * (one- Occ(k) -Occ(i)-Occ(l)-Occ(j))).gt.1.d-5) then
                                    if (rs == 37.and.pq == 78)then
                     write(*, '(A10, 6I5, 2F20.15, I5)') '8temptemp', rs, pq, k, i, l, j,  val, ASing(rs,pq), lr
                  end if
            end if

      end subroutine update_ASing8


      subroutine init_pperpa(ACAlpha, HNO, AuxCoeff, IGem, NBasis)

            double precision, intent(in) :: ACAlpha
            double precision, dimension(:,:), intent(inout) :: HNO
            double precision, dimension(:,:,:,:), intent(inout) :: AuxCoeff
            integer, dimension(:), intent(in) :: IGem
            integer, intent(in) :: NBasis

            integer :: i, j, k, l

            do j=1,NBasis
                  do i=1,NBasis
                        if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
                  enddo
            enddo

            do l=1,3
                  do k=1,3
                        do j=1,3
                              do i=1,3
                                    if((i==j).and.(j==k).and.(k==l).and.(i==2)) then
                                          AuxCoeff(i,j,k,l) = 1
                                    else
                                          AuxCoeff(i,j,k,l) = ACAlpha
                                    endif
                              enddo
                        enddo
                  enddo
            enddo

      end subroutine init_pperpa




end module ppAC0_fofo
