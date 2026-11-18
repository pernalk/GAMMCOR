module interface_sapt
      ! SAPT with the PySCF interface: reads pyscf_data_A.h5 and pyscf_data_B.h5.
      use types
      use h5_reader
      !use sapt_inter, only : writeoneint
      use sapt_inter

      implicit none

      character(len=*), parameter :: PYSCF_H5_A = 'pyscf_data_A.h5'
      character(len=*), parameter :: PYSCF_H5_B = 'pyscf_data_B.h5'

contains

      subroutine sapt_pyscf_driver(Flags,SAPT,NAO,NMO)
            ! Entry point from mainp for JobType SAPT + Interface PYSCF.
            implicit none
            type(FlagsData) :: Flags
            type(SaptData)  :: SAPT
            integer,intent(out) :: NAO,NMO

            call h5_init()
            call read_monomer_pyscf(SAPT%monA,Flags,PYSCF_H5_A,'ONEEL_A',NAO,NMO)
            call read_monomer_pyscf(SAPT%monB,Flags,PYSCF_H5_B,'ONEEL_B',NAO,NMO)

            call print_occ(NMO,SAPT,Flags%ICASSCF)
            call sapt_eris(SAPT,PYSCF_H5_A,PYSCF_H5_B,NAO)
            call sapt_test_energy(SAPT,SAPT%monA,NAO,NMO)

            call h5_finish()

            ! TODO: Cholesky vectors, then hand over to the sapt_driver steps
            !stop 'sapt_pyscf_driver: reading done, rest not implemented'

      end subroutine sapt_pyscf_driver

      subroutine read_monomer_pyscf(Mon,Flags,h5file,onefile,NAO,NMO)
            ! Dimensions, 1-el integrals, nuclei, orbitals and occupations of one monomer.
            implicit none
            type(SystemBlock)            :: Mon
            type(FlagsData)              :: Flags
            character(len=*), intent(in) :: h5file,onefile
            integer, intent(out)         :: NAO, NMO

            integer(HID_T) :: fid
            integer :: i,j
            integer :: NBasis
            integer :: NCenAB,NCen,NEl
            integer :: natural_orb
            integer, allocatable :: dims(:)
            integer, allocatable :: icharg(:)
            double precision, allocatable :: S(:,:),V(:,:),H(:,:),xyz(:,:)
            character(:), allocatable :: rdm2file
            character(len=16) :: eri_basis

            ! set filenames
            if(Mon%Monomer==1) rdm2file = 'rdm2_A.dat'
            if(Mon%Monomer==2) rdm2file = 'rdm2_B.dat'

            fid = h5_open(h5file)

            ! dimensions and nuclear repulsion
            call h5_get(fid,'/REF/NBASIS',NBasis)
            call h5_get(fid,'/REF/NI',Mon%NI)
            call h5_get(fid,'/REF/NA',Mon%NA)
            call h5_get(fid,'/REF/NV',Mon%NV)
            call h5_get(fid,'/REF/NEL',NEl)
            call h5_get(fid,'/REF/ENUC',Mon%PotNuc)
            call h5_get(fid,'/REF/NATORB',natural_orb)
            Mon%NBasis = NBasis
            Mon%NIA    = Mon%NI + Mon%NA
            Mon%NELE   = NEl/2
            Mon%XELE   = NEl/2d0
            ! variable names used in SAPT
            Mon%INAct  = Mon%NI
            Mon%NAct   = Mon%NA
            Mon%NOccup = Mon%NIA
            print*, 'natural_orb = ', natural_orb

            allocate(Mon%IGem(NBasis))
            if(Mon%INAct==0) then
               Mon%NGem = 2
               Mon%IGem(1:Mon%NAct+Mon%INAct)        = 1
               Mon%IGem(Mon%NAct+Mon%INAct+1:NBasis) = 2
            else
               Mon%NGem = 3
               Mon%IGem(1:Mon%INAct) = 1
               Mon%IGem(Mon%INAct+1:Mon%INAct+Mon%NAct) = 2
               Mon%IGem(Mon%INAct+Mon%NAct+1:NBasis)    = 3
            endif

      !      ! S, V (own nuclei only), H = T + V; stored like ONEEL_A/B
            allocate(S(NBasis,NBasis),V(NBasis,NBasis),H(NBasis,NBasis))
            call h5_get(fid,'/INTS/OVERLAP',S)
            call h5_get(fid,'/INTS/POTENTIAL',V)
            call h5_get(fid,'/INTS/CORE_HAMILTONIAN',H)
            call writeoneint(onefile,NBasis**2,S,V,H)

            allocate(Mon%HNO0(NBasis,NBasis),Mon%V(NBasis,NBasis))
            Mon%HNO0 = H
            Mon%V    = V

            ! nuclei: ghost atoms have charge 0; GEOMETRY (natm,3) arrives as (3,natm)
            NCenAB = h5_size(fid,'/SYSTEM/ATOM_CHARGES')
            allocate(icharg(NCenAB),xyz(3,NCenAB))
            call h5_get(fid,'/SYSTEM/ATOM_CHARGES',icharg)
            call h5_get(fid,'/SYSTEM/GEOMETRY',xyz)

            NCen = 0
            do i=1,NCenAB
               if (icharg(i) .gt. 0) NCen = NCen + 1 
            enddo
            print*, 'NCen =', NCen
            print*, 'icharg', icharg
            !print*, 'xyz =', xyz
            allocate(Mon%charg(NCen),Mon%xyz(NCen,3))
            j = 0
            do i=1,NCenAB
               if (icharg(i) .gt. 0) then 
                  j = j + 1
                  Mon%charg(j) = dble(icharg(i))
                  Mon%xyz(j,1) = xyz(1,i) 
                  Mon%xyz(j,2) = xyz(2,i) 
                  Mon%xyz(j,3) = xyz(3,i) 
               endif
            enddo

            ! orbitals C(ao,mo) and occupations 0..1 (natural orbitals)

            !allocate(Mon%CMO(NBasis,NBasis),Mon%Occ(NBasis))
            !call h5_get(fid,'/POSTHF/ORB_COEFF',Mon%CMO)
            !call h5_get(fid,'/POSTHF/OCC',Mon%Occ)

            call h5_shape(fid, '/POSTHF/ORB_COEFF', dims)

            if (size(dims) /= 2) then
               error stop 'POSTHF/ORB_COEFF must be a rank-2 dataset'
            end if
            
            NAO = dims(1)
            NMO = dims(2)

            if (NAO /= NMO) print*, 'Warning! NAO .ne. NMO'
            
            allocate(Mon%CMO(NAO,NMO))
            call h5_get(fid, '/POSTHF/ORB_COEFF', Mon%CMO)
            
            deallocate(dims)

            allocate(Mon%Occ(NBasis))
            call h5_get(fid,'/POSTHF/OCC',Mon%Occ)

            ! construct CICoef
            allocate(Mon%CICoef(NBasis))
            do i=1,NBasis
               Mon%CICoef(i)=sqrt(Mon%Occ(I))
               if(Mon%Occ(i).lt.0.5d0) Mon%CICoef(i)=-Mon%CICoef(i)
            enddo

            ! check D 
            block
            integer :: j
            real(8) :: val
            real(8) :: D(nbasis,nbasis),C_occ(nbasis,nbasis)

            do i=1,nmo
               C_occ(1:nao,i) = 2d0*Mon%Occ(i) * Mon%CMO(1:nao,i)
            end do
            call dgemm('N','T',nao,nao,nmo,1.d0,C_occ,nao,Mon%CMO,nao,0.d0,D,nao)

            val = 0d0
            do j=1,nao
            do i=1,nao
               val = val + D(i,j)*S(i,j)
            enddo
            enddo
            print*, 'Tr(S.rho) = ', val


            end block

            ! check : CT.S.C = 1
            !print*, 'Monomer =', Mon%Monomer, norm2(Mon%CMO)
            !block
            ! integer :: N
            !N = NBasis
            !real(8) :: C(N,N),CTSC(N,N)
            !real(8) :: work(N,N)
            !C = Mon%CMO
            !! Work = S * C
            !call dgemm('N','N',NBasis,NBasis,NBasis, &
            !           1.0d0, S, NBasis, C, NBasis, 0.0d0, Work, NBasis)
            !! CTSC = C^T * Work
            !call dgemm('T', 'N', NBasis, NBasis, NBasis, &
            !           1.0d0, C, NBasis, Work, NBasis, 0.0d0, CTSC, NBasis)
            !print*, 'CTSC'
            !call print_sqmat(CTSC,NBasis)
            !end block

            Mon%SumOcc = sum(Mon%Occ)
            call print_sqmat(Mon%CMO,NBasis)
            
            ! symmetry sorting
            
            ! load density matrices
            call read_1rdm_pyscf(Mon,fid,.false.)
            call read_2rdm_pyscf(Mon,fid,.true.)
            call dump_2rdmfile(Mon,rdm2file)

            ! read dipole moment integrals in AO
            allocate(Mon%dipm(NBasis,NBasis,3))
            call h5_get(fid,'/INTS/DIPOLE_X',Mon%dipm(:,:,1))
            call h5_get(fid,'/INTS/DIPOLE_Y',Mon%dipm(:,:,2))
            call h5_get(fid,'/INTS/DIPOLE_Z',Mon%dipm(:,:,3))

            call compute_dipmom_pyscf(Mon,NBasis)

            ! look-up tables
            call select_active(Mon,NBasis,Flags)
            print*, 'NDimX = ', Mon%NDimX

            call h5_close(fid)

            write(LOUT,'(1x,a,i5,a,i3,a,i3)') 'Read '//h5file//': NBasis',NBasis, &
                  '  INAct',Mon%INAct,'  NAct',Mon%NAct

            deallocate(V,S,H)

      end subroutine read_monomer_pyscf

      subroutine read_1rdm_pyscf(Mon,pyscf_fid,only_full)
            ! read 1-RDM from PySCF
            implicit none

          type(SystemBlock) :: Mon
          integer(HID_T) :: pyscf_fid
          logical, optional, intent(in) :: only_full

          !character(len=256) :: grp
          character(len=256) :: name_aa, name_bb
          character(len=256) :: name_rdm1
          logical :: exists, ex_bbbb, ex_baba, ex_aa, ex_bb, ex_full
          logical :: use_full

          use_full = only_full

          !grp = pyscf_state_group(suffix)
          name_aa   = 'POSTHF/STATES/1.1/RDM1_A'
          name_bb   = 'POSTHF/STATES/1.1/RDM1_B'
          name_rdm1 = 'POSTHF/STATES/1.1/RDM1'
          allocate(Mon%rdm1_p(Mon%NAct, Mon%NAct))
          allocate(Mon%rdm1_m(Mon%NAct, Mon%NAct))
          allocate(Mon%rdm1_full(Mon%NAct, Mon%NAct))

          ex_full = h5_exists(pyscf_fid, trim(name_rdm1))
          if(ex_full)then
                call h5_get(pyscf_fid, trim(name_rdm1), Mon%rdm1_full)
          end if
          print*, '1-rdm', norm2(Mon%rdm1_full)
          print*, Mon%rdm1_full

          ex_aa = h5_exists(pyscf_fid, trim(name_aa))

          if (ex_aa) then
                call h5_get(pyscf_fid, trim(name_aa), Mon%rdm1_p)
          else
                Mon%rdm1_p = Mon%rdm1_full / 2d0
          end if

          ex_bb = h5_exists(pyscf_fid, trim(name_bb))
          if (ex_bb) then
                call h5_get(pyscf_fid, trim(name_bb), Mon%rdm1_m)
          else
                Mon%rdm1_m = Mon%rdm1_full / 2d0
          end if

      end subroutine read_1rdm_pyscf

      subroutine read_2rdm_pyscf(Mon,pyscf_fid,only_full)
      !
      ! for now: read spin-summed 2-RDM
      !
      ! to-do  : use only_full to switch 
      !          for spin-resolved 2-rdms
      !
          type(SystemBlock) :: Mon
          integer(HID_T) :: pyscf_fid
          logical, optional, intent(in) :: only_full

          logical :: exists, ex_bbbb, ex_baba, ex_full
          logical :: use_full

          character(len=256) :: name_rdm2
   
          use_full = only_full

          allocate(Mon%rdm2_full(Mon%NAct,Mon%NAct,Mon%NAct,Mon%NAct))

          if (only_full) then
             name_rdm2 = '/POSTHF/STATES/1.1/RDM2_FULL_REORDERED'
             exists = h5_exists(pyscf_fid, trim(name_rdm2))
             if (exists) then
                   call h5_get(pyscf_fid, trim(name_rdm2), Mon%rdm2_full)
                   print*, 'rdm2-full =', norm2(Mon%rdm2_full)
             else
                   write(*, '(1x,3a)') 'Dataset ', trim(name_rdm2), &
                         ' not found - assembling rdm2_full from the spin blocks'
                   use_full = .false.
             end if
          end if

      end subroutine read_2rdm_pyscf

      subroutine dump_2rdmfile(Mon,rdm2file)
      ! 
      ! save spin-summed 2RDM to rdm2.dat files
      !   (used in Hessian subroutines)
      ! 
          implicit none

          type(SystemBlock) :: Mon
          character(*) :: rdm2file

          integer :: i,j,k,l,ij,kl
          integer :: iunit

          ! save to rdm2.dat files
          open(newunit=iunit,file=rdm2file,status='replace',&
               form='formatted')
          do i=1,Mon%NAct
            do j=1,Mon%NAct
               ij = (i-1)*Mon%NAct+j
               do k=1,Mon%NAct
                  do l=1,Mon%NAct
                     kl = (k-1)*Mon%NAct+l
                     if(ij>=kl) then
                       write(iunit,'(4i4,f19.12)') &
                       !write(6,'(4i4,f19.12)') &
                         k,i,l,j,Mon%rdm2_full(i,j,k,l)
                     endif
                  enddo
               enddo
            enddo
          enddo
          close(iunit)

      end subroutine dump_2rdmfile

      subroutine compute_dipmom_pyscf(Mon,NBasis)
      !
      ! compute dipole moment
      !   a) AO2NO transform Dip
      !
      implicit none

          type(SystemBlock)  :: Mon
          integer,intent(in) :: NBasis

          integer :: i, ncen, NOccup
          double precision :: nuc_dmx, nuc_dmy, nuc_dmz 
          double precision :: dm_x, dm_y, dm_z, dxyz
          double precision :: work(NBasis,NBasis)

          !stop "compute_dip_pyscf not ready!"
          NOccup = Mon%INAct + Mon%NAct

          associate(dipz=>Mon%dipm(:,:,3),dipy=>Mon%dipm(:,:,2),dipx=>Mon%dipm(:,:,1),&
                    UAONO=>Mon%CMO,Occ=>Mon%Occ,xyz=>Mon%xyz,charg=>Mon%charg)

          ! transform to NO
          call dgemm('T','N',NBasis,NBasis,NBasis,1d0,UAONO,NBasis,dipz,NBasis,0d0,work,NBasis)
          call dgemm('N','N',NBasis,NBasis,NBasis,1d0,work,NBasis,UAONO,NBasis,0d0,dipz,NBasis)
          call dgemm('T','N',NBasis,NBasis,NBasis,1d0,UAONO,NBasis,dipy,NBasis,0d0,work,NBasis)
          call dgemm('N','N',NBasis,NBasis,NBasis,1d0,work,NBasis,UAONO,NBasis,0d0,dipy,NBasis)
          call dgemm('T','N',NBasis,NBasis,NBasis,1d0,UAONO,NBasis,dipx,NBasis,0d0,work,NBasis)
          call dgemm('N','N',NBasis,NBasis,NBasis,1d0,work,NBasis,UAONO,NBasis,0d0,dipx,NBasis)

          ncen = size(charg)
          !print*, 'ncen =', ncen
          !print*, 'noccup =', noccup
          nuc_dmx=0d0; nuc_dmy=0d0; nuc_dmz=0d0
          do i=1,ncen
             nuc_dmx = nuc_dmx + charg(i)*xyz(i,1)
             nuc_dmy = nuc_dmy + charg(i)*xyz(i,2)
             nuc_dmz = nuc_dmz + charg(i)*xyz(i,3)
          enddo

          write(6,'(/,1X,A,3f14.8)') 'Nuclear Dipole Moment   ', &
                                nuc_dmx,nuc_dmy,nuc_dmz
!
          dm_x=0d0; dm_y=0d0; dm_z=0d0
          do i=1,noccup
             dm_x = dm_x + 2d0*occ(i)*dipx(i,i)
             dm_y = dm_y + 2d0*occ(i)*dipy(i,i)
             dm_z = dm_z + 2d0*occ(i)*dipz(i,i)
          enddo
          write(6,'(1x,a,3f14.8)')   'Electronic dipole moment', &
                                dm_x,dm_y,dm_z

          write(6,'(1X,A,3f14.8,/)') 'Total Dipole Moment     ', &
                                nuc_dmx+dm_x,nuc_dmy+dm_y,nuc_dmz+dm_z

          dxyz=sqrt((nuc_dmx+dm_x)**2+(nuc_dmy+dm_y)**2+(nuc_dmz+dm_z)**2)
          write(6,'(1X,A,2f14.8,/)') '|dipole moment| a.u./D', dxyz,dxyz/0.393456

          end associate


      end subroutine compute_dipmom_pyscf

      subroutine sapt_eris(SAPT,h5fileA,h5fileB,NAO) 
      !
      !     read 2-el integrals
      !       -- read ERIs in AO basis in triangle form
      !
            implicit none

            type(SaptData)               :: SAPT
            integer,intent(in)           :: NAO
            character(len=*), intent(in) :: h5fileA,h5fileB

            integer(HID_T) :: fidA, fidB, fidERI
            integer        :: i,j,p,q,r,s
            integer        :: NInte1,NInte2
            character(len=16) :: eri_basis
            character(len=:), allocatable :: eri_file

            fidA = h5_open(h5fileA)
            fidB = h5_open(h5fileB)

            ! read 2-el integrals
            ! prefer A; try B when A has no ERIs

            if (h5_exists(fidA, '/MOINTS/ERI')) then
               fidERI = fidA
               eri_file = h5fileA
            else if (h5_exists(fidB, '/MOINTS/ERI')) then
               fidERI = fidB
               eri_file = h5fileB
            else
               print *, 'ERROR: no MOINTS/ERI dataset in either file:'
               print *, trim(h5fileA)
               print *, trim(h5fileB)
               error stop
            end if

            call h5_attr(fidERI, '/MOINTS/ERI', 'BASIS', eri_basis)

            if (trim(eri_basis) /= 'AO') then
               print *, 'MOINTS/ERI: expected AO integrals, found ', trim(eri_basis)
               error stop
            end if

            NInte1 = NAO*(NAO+1)/2
            NInte2 = NInte1*(NInte1+1)/2

            allocate(SAPT%TwoEl(NInte2))
            call h5_get(fidERI, '/MOINTS/ERI', SAPT%TwoEl)
            print*, 'Two electron integrals read from ', eri_file

            call h5_close(fidA)
            call h5_close(fidB)

            call unpack_twoel(SAPT,NAO)
            ! keep only unpacked TwoEl4
            deallocate(SAPT%TwoEl)

            ! example fofo
            !! FOFO AB|BA
            !call tran4gen_incore(SAPT%TwoEl,SAPT%TwoEl4,NBasis,&
            !                     SAPT%monA%CMO,NBasis, &
            !                     SAPT%monB%CMO,(SAPT%monB%num0+SAPT%monB%num1), &
            !                     SAPT%monB%CMO,NBasis, &
            !                     SAPT%monA%CMO,(SAPT%monA%num0+SAPT%monA%num1))
            !call write2disk(SAPT%TwoEl4, &
            !                NBasis, &
            !                SAPT%monB%num0+SAPT%monB%num1, &
            !                NBasis, &
            !                SAPT%monA%num0+SAPT%monA%num1, &
            !                'FOFOABBA')

      end subroutine sapt_eris

      subroutine unpack_twoel(SAPT,NAO)
      !
      ! unpack NAddr3 to NAO**4
      !
      implicit none

            type(SaptData)     :: SAPT
            integer,intent(in) :: nao

            integer :: p,q,r,s
            integer,external :: NAddr3

            allocate(SAPT%TwoEl4(nao,nao,nao,nao))
            do s=1,NAO
               do r=1,NAO
                  do q=1,NAO
                     do p=1,NAO
                        SAPT%TwoEl4(p,q,r,s) = SAPT%TwoEl(NAddr3(p,q,r,s))
                     enddo
                  enddo
               enddo
            enddo

      end subroutine unpack_twoel

      subroutine tran4gen_incore(A,C,nao,Ca,nmo_a,Cb,nmo_b,Cc,nmo_c,Cd,nmo_d)

       implicit none 

       integer,intent(in)   :: nao,nmo_a,nmo_b,nmo_c,nmo_d
       real(8),intent(in)   :: A(nao**4)
       real(8),intent(in)   :: Ca(nao,nmo_a), Cb(nao,nmo_b), Cc(nao,nmo_c), Cd(nao,nmo_d)
       real(8),intent(out)  :: C(nmo_a*nmo_b*nmo_c*nmo_d)

       integer(8) :: sz_w1, sz_w2
       real(8), allocatable :: W1(:), W2(:)
       
       sz_w1 = max(int(nao, 8)**3 * nmo_a, int(nao, 8) * nmo_a * nmo_b * nmo_c)
       sz_w2 = int(nao, 8)**2 * nmo_a * nmo_b
       
       allocate(W1(sz_w1), W2(sz_w2))
       
       ! Step 1: A(nao, nao, nao, nao) -> W1(nao, nao, nao, nmo_a)
       call dgemm('T', 'N', nao*nao*nao,       nmo_a, nao, 1.d0, A,  nao, Ca, nao, 0.d0, W1, nao*nao*nao)
       
       ! Step 2: W1(nao, nao, nao, nmo_a) -> W2(nao, nao, nmo_a, nmo_b)
       call dgemm('T', 'N', nao*nao*nmo_a,     nmo_b, nao, 1.d0, W1, nao, Cb, nao, 0.d0, W2, nao*nao*nmo_a)
       
       ! Step 3: W2(nao, nao, nmo_a, nmo_b) -> W1(nao, nmo_a, nmo_b, nmo_c)
       call dgemm('T', 'N', nao*nmo_a*nmo_b,   nmo_c, nao, 1.d0, W2, nao, Cc, nao, 0.d0, W1, nao*nmo_a*nmo_b)
       
       ! Step 4: W1(nao, nmo_a, nmo_b, nmo_c) -> C(nmo_a, nmo_b, nmo_c, nmo_d)
       call dgemm('T', 'N', nmo_a*nmo_b*nmo_c, nmo_d, nao, 1.d0, W1, nao, Cd, nao, 0.d0, C,  nmo_a*nmo_b*nmo_c)
       
       deallocate(W1, W2)

     end subroutine tran4gen_incore
     
     subroutine write2disk(C,nmo_a,nmo_b,nmo_c,nmo_d,fname)

        implicit none
        integer,intent(in)   :: nmo_a,nmo_b,nmo_c,nmo_d
        real(8),intent(in)   :: C(nmo_a*nmo_b,nmo_c*nmo_d)
        character(*),intent(in) :: fname

        integer :: iunit
        integer :: i,j
        integer :: nAB,nCD

        nAB = nmo_a*nmo_b
        nCD = nmo_c*nmo_d

        open(newunit=iunit,file=fname,status='REPLACE',&
            access='DIRECT',form='UNFORMATTED',recl=8*nAB)

        do j=1,nCD
           write(iunit,rec=j) (C(i,j),i=1,nAB)
        enddo

        close(iunit)

     end subroutine write2disk

     subroutine sapt_test_energy(SAPT,Mon,NAO,NMO)

            implicit none

            type(SaptData)      :: SAPT
            type(SystemBlock)   :: Mon
            integer, intent(in) :: NAO, NMO

            integer :: i, j, p, q ,r, s
            integer :: iunit, irec
            real(8) :: eone, etwo, val
            real(8) :: HNO0_mo(NAO,NAO)
            real(8) :: AUXM(NMO,NAO), ints(NMO,NMO)
            real(8) :: work(NMO**2)
            real(8),allocatable :: W1(:,:,:,:)

          associate(NI=>mon%INAct, NA=>mon%NAct, NIA=>mon%NOccup, &
                    HNO0=>mon%HNO0, Occ=>mon%Occ,UAONO=>mon%CMO, &
                    IndAux=>mon%IndAux)

            ! FFOO
            print*, 'NOccupA =', NIA
            print*, 'check =', mon%num0+mon%num1
            allocate(W1(NMO,NMO,NIA,NIA))
            call tran4gen_incore(SAPT%TwoEl4,W1,NAO,&
                                 mon%CMO,NMO, &
                                 mon%CMO,NMO, &
                                 mon%CMO,(mon%num0+mon%num1),&
                                 mon%CMO,(mon%num0+mon%num1))
                     print*, 'W1 =', norm2(W1)
            call write2disk(W1, &
                            NMO, &
                            NMO, &
                            mon%num0+mon%num1, &
                            mon%num0+mon%num1, &
                            'FFOOAA')

            call dgemm('T','N',NMO,NAO,NAO,1d0,UAONO,NAO,HNO0,NAO,0d0,AUXM,NMO)
            call dgemm('N','N',NMO,NMO,NAO,1d0,AUXM,NMO,UAONO,NAO,0d0,HNO0_mo,NMO)
            print*, 'XOne =', norm2(HNO0_mo)

            eone = 0d0
            do p=1,NIA
               eone = eone + 2d0*Occ(p)*HNO0_mo(p,p)
            enddo
            write(6,'(1X,''CASSCF Energy (one-electron) calculated'',X,F15.8)')eone

            etwo = 0d0
            do s=1,NIA
               do r=1,NIA
                  do q=1,NIA
                     do p=1,NIA
                        val = 0d0
                        if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                              val = val +  0.5d0 * SAPT%monA%rdm2_full(p-NI, q-NI, r-NI, s-NI)
                        else
                              if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                  val =  val + 2d0 * Occ(p) * Occ(q)
                              end if

                              if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                                  val = val  -Occ(p) * Occ(q)
                              end if
                        end if
                        etwo = etwo + val * W1(p,r,q,s)
                     enddo
                  enddo
               enddo
            enddo
            write(6,'(1X,''CASSCF Two-electron Energy calculated'',X,F15.8)')etwo
            write(6,'(1X,''CASSCF Energy (w/o ENuc-q) calculated'',X,F15.8)')eone+etwo

            deallocate(W1)

          print*, 'test FFOO file...'

          open(newunit=iunit,file='FFOOAA',status='OLD', &
               access='DIRECT',recl=8*NMO**2)

          etwo=0d0 
          ! COULOMB LOOP (FF|OO)
          do r=1,NIA
          do p=1,NIA
             irec = p + (r-1)*NIA
             read(iunit,rec=irec) work(1:NMO**2)
             Do j=1,NMO
             Do i=1,NMO
             ints(i,j) = work((j-1)*NMO+i)
             EndDo
             EndDo

             do s=1,NIA
             do q=1,NIA
                val = 0d0
                if (IndAux(p)==1.and.IndAux(q)==1.and.IndAux(r)==1.and.IndAux(s)==1)then
                      val = val +  0.5d0 * SAPT%monA%rdm2_full(p-NI, q-NI, r-NI, s-NI)
                else
                      if (p==r.and.q==s.and.(IndAux(p)==0.or. IndAux(q)==0))then
                          val =  val + 2d0 * Occ(p) * Occ(q)
                      end if

                      if (p==s.and.q==r.and.(IndAux(p)==0.or. IndAux(q)==0))then
                          val = val  -Occ(p) * Occ(q)
                      end if
                end if
                etwo = etwo + val * ints(q,s)
             enddo
             enddo

          enddo
          enddo
          write(6,'(1X,''CASSCF Two-electron Energy 2         '',X,F15.8)')etwo

          close(iunit,status='DELETE')

          end associate

     end subroutine sapt_test_energy

end module interface_sapt
