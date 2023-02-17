subroutine write_molpro_2_trexio(NBasis)
use trexio
use types
use tran
use read_external
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NBasis

character*(32), parameter:: filename = 'test.trexio.text'
integer    :: rc  ! return code of trexio functions
integer(8) :: f   ! TREXIO file

integer :: i,j,k,l,m,kl
integer :: NInte1

!buffers of integrals
!integer(8), parameter :: BUFSIZE = NBasis**2
!integer :: buffer_index(4,BUFSIZE)
!double precision :: buffer_values(BUFSIZE)
integer(8) :: BUFSIZE
integer,allocatable          :: buffer_index(:,:)
double precision,allocatable :: buffer_values(:)
integer(8) :: offset, icount
integer :: num_basis
double precision :: val
logical :: empty

double precision, allocatable :: work(:)
double precision, allocatable :: Smat(:,:),overlap(:,:)
double precision, allocatable :: Kmat(:,:),kinetic(:,:)

print*, 'TreX-io test'

NInte1 = NBasis*(NBasis+1)/2
BUFSIZE = NInte1
allocate(buffer_index(4,BUFSIZE),buffer_values(BUFSIZE))

f = trexio_open (filename, 'w', TREXIO_TEXT, rc)

! write
rc = trexio_write_mo_num(f, NBasis)
if (rc /= TREXIO_SUCCESS) then
  stop 'Error writing MO num'
end if

! test overlap
allocate(work(NInte1),Smat(NBasis,NBasis),Kmat(NBasis,NBasis))

call readoneint_molpro(work,'AOONEINT.mol','OVERLAP',.true.,NInte1)
call triang_to_sq2(work,Smat,NBasis)
print*, 'Smat   ',norm2(SMat)

rc = trexio_write_mo_1e_int_overlap(f, Smat)
if (rc /= TREXIO_SUCCESS) then
  stop 'Sth wrong in overlap'
end if

! test kinetic
call readoneint_molpro(work,'AOONEINT.mol','KINETINT',.true.,NInte1)
call triang_to_sq2(work,Kmat,NBasis)
print*, 'Kmat   ',norm2(KMat)

rc = trexio_write_mo_1e_int_kinetic(f, Kmat)
if (rc /= TREXIO_SUCCESS) then
  stop 'Sth wrong in kinetic'
end if

! test ERIs
! l > k and j > i
call reader%open(trim("AOTWOSORT"))

icount = BUFSIZE
offset = 0
val = 0
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1

      call reader%getTR(kl,work,empty)
      if(empty) cycle

      m = 0
      do j=1,NBasis
         do i=1,j
            m = m + 1
            !'<i k|j l> = buffer_values(m)
             buffer_index(1,m) = i
             buffer_index(2,m) = k
             buffer_index(3,m) = j
             buffer_index(4,m) = l
             buffer_values(m)  = work(m)
             val = val + buffer_values(m)**2
         enddO
      enddo
      offset = offset + m
      rc = trexio_write_mo_2e_int_eri(f, offset, icount, buffer_index, buffer_values)

   enddo
enddo
print*, 'offset',offset
print*, 'val-1',val

call reader%close

rc = trexio_close(f)

deallocate(Kmat,Smat,work)
deallocate(buffer_values,buffer_index)

end subroutine write_molpro_2_trexio

subroutine read_molpro_2_trexio()
use trexio
use types
implicit none

character*(32), parameter:: filename = 'test.trexio.text'
integer    :: rc
integer(8) :: f

integer :: i,j,k,l,m

!buffers of integrals
integer(8) :: BUFSIZE
integer,allocatable          :: buffer_index(:,:)
double precision,allocatable :: buffer_values(:)
integer(8) :: offset, icount

integer          :: num_basis
double precision :: val

double precision, allocatable :: ints(:,:)
double precision, allocatable :: overlap(:,:)
double precision, allocatable :: kinetic(:,:)

f = trexio_open (filename, 'r', TREXIO_TEXT, rc)
rc = trexio_read_mo_num(f, num_basis)

print *,  'Number of MOs: ', num_basis
allocate(overlap(num_basis,num_basis), kinetic(num_basis,num_basis))

rc = trexio_has_mo_1e_int_overlap(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No overlap in file'
end if

rc = trexio_read_mo_1e_int_overlap(f, overlap)
print*, 'overlap',norm2(overlap)

rc = trexio_has_mo_1e_int_kinetic(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No kinetic in file'
end if

rc = trexio_read_mo_1e_int_kinetic(f, kinetic)
print*, 'kinetic',norm2(kinetic)

! eri
allocate(ints(num_basis,num_basis))
!BUFSIZE = 2*num_basis**2
BUFSIZE = 600
allocate(buffer_index(4,BUFSIZE),buffer_values(BUFSIZE))

icount = BUFSIZE
offset = 0
val = 0
do while(icount == BUFSIZE)

   rc = trexio_read_mo_2e_int_eri(f, offset, icount, buffer_index, buffer_values)

   !if(rc /= TREXIO_SUCCESS) then
   !   stop 'Error reading ERI'
   !end if

   offset = offset + icount

   do i=1,icount
      val = val + buffer_values(i)**2
   enddo

enddo

print*, 'offset',offset
print*, 'val-2',val

rc = trexio_close(f)
deallocate(kinetic,overlap)

deallocate(buffer_values,buffer_index)

end subroutine read_molpro_2_trexio

subroutine testonel_trexio(NBasis)
 use trexio
 implicit none

 integer,intent(in) :: NBasis

 integer    :: rc
 integer    :: NAO, NMO, IsCartesian
 integer(8) :: f
 character(:),allocatable     :: TreXFile
 double precision,allocatable :: Vao(:),Vmo(:)
 double precision,allocatable :: Coeff(:)
 double precision,allocatable :: work(:),Vao2mo(:)

TreXfile = "h2o-xxx.h5"

f = trexio_open (TrexFile, 'r', TREXIO_HDF5, rc)

rc = trexio_has_ao_num(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No AO num in file'
end if
rc = trexio_read_ao_num(f,NAO)

rc = trexio_has_mo_num(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No MO num in file'
end if
rc = trexio_read_mo_num(f,NMO)

print*, 'AO num:',NAO
print*, 'MO num:',NMO

rc = trexio_read_ao_cartesian(f,IsCartesian)
print*, 'IsCartesian: ', IsCartesian

allocate(Vao(NAO*NAO),Vmo(NMO*NMO))
allocate(Coeff(NAO*NMO))

rc = trexio_has_ao_1e_int_potential_n_e(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No AO potential_n_e in file'
end if
rc = trexio_read_ao_1e_int_potential_n_e(f, Vao)

print*, 'Vmat-AO'
call print_sqmat_cols(Vao,NAO)

rc = trexio_has_mo_coefficient(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No AOMO coefficients in file'
end if
rc = trexio_read_mo_coefficient(f, Coeff)

call print_sqmat_cols(Coeff,NAO,NMO)

allocate(work(NAO*NMO),Vao2mo(NBasis**2))

! CT.A.C
call dgemm('T','N',NMO,NAO,NAO,1d0,Coeff,NAO,Vao,NAO,0d0,work,NMO)
call dgemm('N','N',NMO,NMO,NAO,1d0,work,NMO,Coeff,NAO,0d0,Vao2mo,NMO)

print*, 'Vmat-AO2MO : CT.A.C'
call print_sqmat_cols(Vao2mo,NMO)

! C.A.CT
call dgemm('N','N',nbasis,nbasis,nbasis,1d0,Coeff,nbasis,Vao,nbasis,0d0,work,nbasis)
call dgemm('N','T',nbasis,nbasis,nbasis,1d0,work,nbasis,Coeff,nbasis,0d0,Vao2mo,nbasis)

print*, 'Vmat-AO2MO : C.A.CT'
call print_sqmat_cols(Vao2mo,nbasis)

deallocate(Vao2mo,work)

rc = trexio_has_mo_1e_int_potential_n_e(f)
if (rc /= TREXIO_SUCCESS) then
  stop 'No MO potential_n_e in file'
end if
rc = trexio_read_mo_1e_int_potential_n_e(f, Vmo)

print*, 'Vmat-MO'
call print_sqmat_cols(Vmo,Nbasis)

rc = trexio_close(f)

deallocate(Vao,Vmo,Coeff)

end subroutine testonel_trexio

subroutine print_sqmat_cols(mat,ndim)
implicit none

integer,intent(in) :: ndim
double precision,intent(in) :: mat(ndim,ndim)

integer :: i,j

 ! print columns
 do j=1,ndim
    write(6,*) j
    write(6,'(10f13.8)') (mat(i,j),i=1,ndim)
 enddo
 write(6,'()')

 return
end subroutine print_sqmat_cols

subroutine tran4_rdm_2e(NBas,nA,CA,nB,CB,nC,CC,nD,CD,fname,srtfile)
implicit none

integer,intent(in) :: NBas
integer,intent(in) :: nA,nB,nC,nD
! CA(NBas*nA)
double precision,intent(in) :: CA(*), CB(*), CC(*), CD(*)
character(*) :: fname,srtfile
double precision, allocatable :: work1(:), work2(:), work3(:,:)
integer :: iunit,iunit2,iunit3
integer :: nsq,ntr,nAB,nCD,nloop
integer,parameter :: cbuf=576
integer :: i,rs,ab
logical :: empty
!  test
integer :: l,k,kl

 write(6,'(1x,a)') 'Transforming integrals for '//fname

 ntr = NBas*(NBas+1)/2
 nsq = NBas**2
 nAB = nA*nB
 nCD = nC*nD

! rbuf = sqrt(cbuf)
! sbuf = sqrt(cbuf)

 ! set no. of squares in buffer
 !nloop = (ntr-1)/cbuf+1
 nloop = (nsq-1)/cbuf+1

 allocate(work1(NBas*NBas),work2(NBas*NBas))
 allocate(work3(cbuf,nAB))

 !call reader%open(trim(srtfile))

 ! half-transformed file
 open(newunit=iunit2,file='TMPMO',status='REPLACE',&
     access='DIRECT',form='UNFORMATTED',recl=8*cbuf)

! (ab|
 do i=1,nloop

!    call read_pass_trexio_2rdm(i,rbuf,sbuf,NBas,RDM2val)

    !! loop over cbuf
    !do rs=(i-1)*cbuf+1,min(i*cbuf,ntr)

    !   work1 = 0
    !   !read(iunit,rec=rs) work1(1:ntr)
    !   call reader%getTR(rs,work1,empty)
    !   if(empty) then

    !      work3(rs-(i-1)*cbuf,1:nAB) = 0

    !   else

    !      call triang_to_sq(work1,work2,NBas)
    !      ! work1=CA^T.work2
    !      ! work2=work1.CB
    !      call dgemm('T','N',nA,NBas,NBas,1d0,CA,NBas,work2,NBas,0d0,work1,nA)
    !      call dgemm('N','N',nA,nB,NBas,1d0,work1,nA,CB,NBas,0d0,work2,nA)
    !      ! transpose
    !      work3(rs-(i-1)*cbuf,1:nAB) = work2(1:nAB)

    !   endif

    !enddo

    !do ab=1,nAB
    !   write(iunit2,rec=(i-1)*nAB+ab) work3(1:cbuf,ab)
    !enddo

 enddo

! close(iunit)
! call reader%close

!! |cd)
! open(newunit=iunit3,file=fname,status='REPLACE',&
!     access='DIRECT',form='UNFORMATTED',recl=8*nCD)
!
! do ab=1,nAB
!
!    do i=1,nloop
!       ! get all (ab| for given |rs)
!       read(iunit2,rec=(i-1)*nAB+ab) work1((i-1)*cbuf+1:min(i*cbuf,ntr))
!    enddo
!
!    call triang_to_sq(work1,work2,NBas)
!    ! work1=CC^T.work2
!    ! work2=work1.CD
!    call dgemm('T','N',nC,NBas,NBas,1d0,CC,NBas,work2,NBas,0d0,work1,nC)
!    call dgemm('N','N',nC,nD,NBas,1d0,work1,nC,CD,NBas,0d0,work2,nC)
!    write(iunit3,rec=ab) work2(1:nCD)
!    !work1 = 0
!    !kl = 0 
!    !do l=1,nD
!    !do k=1,nC
!    !   kl = kl + 1
!    !   work1(kl) = work2((k-1)*NBas+l)
!    !enddo
!    !enddo
!    !write(iunit3,rec=ab) work1(1:nCD)
!
! enddo
!
! deallocate(work1,work2,work3)
! close(iunit3)
! close(iunit2,status='DELETE')

end subroutine tran4_rdm_2e

subroutine get_h0_no(CAONO,Hmat,onefile,NAO,NBasis)
! read Hmat in AO, transform to NO
! [output] :: Hmat(triang)
!
use types
use tran
implicit none

integer,intent(in) :: NAO,NBasis
double precision,intent(in)  :: CAONO(NAO*NBasis)
double precision,intent(out) :: Hmat(NBasis*(NBasis+1)/2)
character(*) :: onefile

integer          :: iunit,i,j
character(8)     :: label
double precision :: workAO(NAO**2),workNO(NBasis**2)

Hmat = 0

open(newunit=iunit,file=onefile,access='sequential', &
     form='unformatted',status='old')

read(iunit)
read(iunit)
read(iunit) label, workAO
if(label=='ONEHAMIL') then
   call tran_AO2MO(workAO,CAONO,CAONO,workNO,NAO,NBasis)
   call sq_to_triang(workNO,Hmat,NBasis)

else
   write(LOUT,'(a)') 'ERROR! ONEHAMIL NOT FOUND IN '//onefile
   stop
endif

close(iunit)

end subroutine get_h0_no

subroutine AB_CIPSI_FOFO(ABPLUS,ABMIN,RDM2val,ETot,URe,Occ,XOne, &
     IndN,IndX,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     IntJFile,IntKFile,ACAlpha,AB1)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM ...
!
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! NO SYMMETRY IS ASSUMED
! COULOMB INTEGRALS ARE READ FROM IntJFile IN (FF|OO) FORMAT
! EXCHANGE INTEGRALS ARE READ FROM IntKFile IN (FO|FO) FORMAT
! 
 use abfofo
 implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim)
double precision,intent(in)  :: ACAlpha
!double precision,intent(in)  :: RDM2val(NBasis,NBasis,NBasis,NBasis)
double precision,intent(in)  :: RDM2val(INActive+NAct,INActive+NAct,INActive+NAct,INActive+NAct)
logical,intent(in) :: AB1

integer :: i,j,k,l,ij,kl,kk,ll
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit1,iunit2,ios
integer :: switch
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,HNOCoef,val
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,parameter :: ThreshCIPSI = 0d0
!double precision,parameter :: ThreshCIPSI = 1d-14

ABPLUS = 0
ABMIN  = 0
ETot   = 0

! set dimensions
NOccup = NAct + INActive
!print*, 'NOccup',NOccup
Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
!call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val
write(lout,'(1x,a,9x,f15.8)') 'CIPSI One-electron Energy', ETot

if(AB1) then
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
         if(IGem(i)==IGem(j)) HNO(i,j) = 0
      enddo
   enddo
else
   do j=1,NBasis
      do i=1,NBasis
         if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
      enddo
   enddo
endif

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

if(AB1) then
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 0
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
else
   do l=1,3
      do k=1,3
         do j=1,3
            do i=1,3
               if((i==j).and.(j==k).and.(k==l)) then
                  AuxCoeff(i,j,k,l) = 1
               else
                  AuxCoeff(i,j,k,l) = ACAlpha
               endif
            enddo
         enddo
      enddo
   enddo
endif

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxI  = 0
AuxIO = 0
WMAT  = 0

if(AB1) then
   HNOCoef = -1
else
   HNOCoef = 1 - ACAlpha
endif

if(AB1) then
   switch=1
else
   switch=0
endif

call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
             RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
             INActive,NOccup,NDim,NDimX,NBasis,NInte1,IntJFile,IntKFile,ACAlpha,switch,ETot)

!print*, 'after JK Loop + ', norm2(ABPLUS)
!print*, 'after JK Loop - ', norm2(ABMIN)

write(LOUT,'(1x,a,5x,f15.8)') "CIPSI Total Energy (w/o ENuc)", ETot
!print*, "CIPSI Energy (w/o ENuc)", ETot

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = IndX(ICol)
   if(ir/=is) then

      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = IndX(IRow)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         ! CIPSI check : removing this condition
         !if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         if(abs(val) .gt. ThreshCIPSI) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         !if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         if(abs(val) .gt. ThreshCIPSI) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo

   endif
enddo

if(AB1) then
   ! symmetrize AB1
   call sq_symmetrize(ABPLUS,NDimX)
   call sq_symmetrize(ABMIN,NDimX)
endif

deallocate(ints,work2,work1)

end subroutine AB_CIPSI_FOFO

subroutine ABPM0_CIPSI_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,RDM2val, &
                      IndN,IndX,IGemIN,NAct,INActive,NOccup, &
                      NDimX,NBasis,NDim,NInte1, &
                      IntJFile,IntKFile,ETot)
use abfofo
implicit none

integer,intent(in)           :: NAct,INActive,NOccup,NDimX,NBasis,NDim,NInte1
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)

double precision,intent(in)  :: RDM2val(NOccup,NOccup,NOccup,NOccup)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
character(*)                 :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(out),optional :: ETot

integer          :: i,j,k,l
integer          :: ip,iq,ir,is,ipq,irs,ICol,IRow
integer          :: iunit,ios
integer          :: NRDM2Act
integer          :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: val,AuxVal,EnOne
double precision :: AuxCoeff(3,3,3,3),C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),&
                    AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision,allocatable :: work1(:),work2(:)
integer,external          :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
EnOne  = 0
if(present(ETot)) ETot = 0

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

allocate(work1(NBasis**2),work2(NBasis**2))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

deallocate(work1)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
EnOne = EnOne + 2*val

!write(LOUT,*) 'ONE ELECTRON ENERGY:', EnOne
if(present(ETot)) ETot = ETot + 2*val

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = 0d0
   enddo
enddo

AuxInd = 0
AuxInd(1:2,1:2) = 1
AuxInd(2,2) = 2

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

call create_pos_ABPL0(pos,IGem,IndN,INActive,NAct,NBasis,NDimX)

AuxI  = 0
AuxIO = 0
WMAT  = 0

if(present(ETot)) then
   call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
                RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                INActive,NOccup,NDimX,NDimX,NBasis,NInte1,IntJFile,IntKFile,0d0,2,ETot)
else
   call JK_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,&
                RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                INActive,NOccup,NDimX,NDimX,NBasis,NInte1,IntJFile,IntKFile,0d0,2)
endif

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

do ICol=1,NDimX
   ir = IndN(1,ICol)
   is = IndN(2,ICol)
   irs = pos(ir,is)

      do IRow=1,NDimX
         ip = IndN(1,IRow)
         iq = IndN(2,IRow)
         ipq = pos(ip,iq)

         val = 0

         if(ip==ir) then
            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
            select case(AuxInd(IGem(ip),IGem(ir)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,is)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,is)
            end select
         endif

         if(iq==is) then
            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
            select case(AuxInd(IGem(iq),IGem(is)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,ir)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         val = 0

         if(iq==ir) then
            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
            select case(AuxInd(IGem(iq),IGem(ir)))
            case(1)
               val = val + Occ(iq)*AuxI(ip,is)
            case(2)
               val = val + Occ(iq)*AuxIO(ip,is)
            end select
         endif

         if(ip==is) then
            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
            select case(AuxInd(IGem(ip),IGem(is)))
            case(1)
               val = val + Occ(ip)*AuxI(iq,ir)
            case(2)
               val = val + Occ(ip)*AuxIO(iq,ir)
            end select
         endif

         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo
enddo

deallocate(work2)

end subroutine ABPM0_CIPSI_FOFO

subroutine AC0BLOCK_CIPSI(RDM2val,Occ,URe,XOne, &
                    IndN,IndX,IGemIN,NAct,INActive,NOccup, &
                    NDimX,NBasis,NDim,NInte1, &
                    IntJFile,IntKFile,A0BlockIV,A0block,nblk,dumpfile,dump)
!
!     A ROUTINE FOR COMPUTING A0=ABPLUS^{(0)}.ABMIN^{(0)}
!                 (FOFO VERSION)
!
use abfofo
use print_units
use blocktypes
!
implicit none

integer,intent(in)           :: NAct,INActive,NOccup,NDimX,NBasis,NDim,NInte1
integer,intent(in)           :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
integer                      :: nblk
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
double precision,intent(in)  :: RDM2val(NOccup,NOccup,NOccup,NOccup)
character(*)                 :: IntJFile,IntKFile
character(*)                 :: dumpfile
integer,intent(in)           :: dump

integer          :: iunit
integer          :: i,j,k,l,kl,ii,ip,iq,ir,is,ipq,irs
integer          :: ipos,jpos,iblk,jblk
integer          :: nblkIN
integer          :: IGem(NBasis),Ind(NBasis),pos(NBasis,NBasis)
integer          :: nAA,nAI(INActive),nAV(INActive+NAct+1:NBasis),nIV
integer          :: tmpAA(NAct*(NAct-1)/2),tmpAI(NAct,1:INActive),&
                    tmpAV(NAct,INActive+NAct+1:NBasis),&
                    tmpIV(INActive*(NBasis-NAct-INActive))
integer          :: limAA(2),limAI(2,1:INActive),&
                    limAV(2,INActive+NAct+1:NBasis),limIV(2)
double precision :: AuxCoeff(3,3,3,3),Aux,val,ETot
double precision :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision :: C(NBasis)

type(EblockData) :: A0block(nblk), A0blockIV

nblkIN = nblk

Ind = 0
do i=1,NAct
   Ind(INActive+i) = i
enddo

! fix IGem
do i=1,INActive
   IGem(i) = 1
enddo
do i=INActive+1,NOccup
   IGem(i) = 2
enddo
do i=NOccup+1,NBasis
   IGem(i) = 3
enddo

do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = 0
            endif
         enddo
      enddo
   enddo
enddo

call ABPM0_CIPSI_FOFO(Occ,URe,XOne,ABPLUS,ABMIN,RDM2val, &
                IndN,IndX,IGemIN,NAct,INActive,NOccup,&
                NDimX,NBasis,NDim,NInte1, &
                IntJFile,IntKFile,ETot)

call create_blocks_ABPL0(nAA,nAI,nAV,nIV,tmpAA,tmpAI,tmpAV,tmpIV,&
                         limAA,limAI,limAV,limIV,pos,&
                         IGem,IndN,INActive,NAct,NBasis,NDimX)

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

!allocate(A0block(1+NBasis-NAct))
nblk = 0

!pack AA
if(nAA>0) then
   nblk = nblk + 1
   call pack_A0block(ABPLUS,ABMIN,nAA,limAA(1),limAA(2),tmpAA,A0block(nblk),NDimX)
endif
!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAI(iq),limAI(1,iq),limAI(2,iq),tmpAI(1:nAI(iq),iq),&
                        A0block(nblk),NDimX)
   endif
enddo
!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then
      nblk = nblk + 1
      call pack_A0block(ABPLUS,ABMIN,nAV(ip),limAV(1,ip),limAV(2,ip),tmpAV(1:nAV(ip),ip),&
                        A0block(nblk),NDimX)
    endif
enddo

! test nblk
print*, 'nblk  ',nblk
print*, 'nblkIN',nblkIN
if((nblk.ne.nblkIN).and.(NAct.ne.0)) then
 print*, 'ERROR! nblocks wrong in AC0BLOCK_CIPSI!'
 print*, 'nblk  ',nblk
 print*, 'nblkIN',nblkIN
 stop
endif

!pack IV
associate(B => A0blockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)*ABMIN(ii,ii)
  enddo

end associate

if(dump==1) then
  ! dump to file
  open(newunit=iunit,file=dumpfile,form='unformatted')
  write(iunit) nblk
  do iblk=1,nblk
     associate(B => A0Block(iblk))
       write(iunit) iblk, B%n, B%l1, B%l2
       write(iunit) B%pos,B%matX
     end associate
  enddo
  associate(B => A0BlockIV)
    write(iunit) B%n,B%l1,B%l2
    write(iunit) B%pos,B%vec
  end associate
  close(iunit)

endif 

end subroutine AC0BLOCK_CIPSI

subroutine ERPASYMMXY_CIPSI(EigY,EigX,Eig,ABPLUS,ABMIN,CICoef,IndN,NDimX,NBasis)
!
! RETURNS X,Y VECTORS (PROPERLY NORMALIZED), WHICH ARE SOLUTIONS
! OF THE ORIGINAL ERPA PROBLEM:
! AX + BY = Om N X
! BX + AY =-Om N Y
!
! THIS IS ACHIEVED BY CALLING ERPASYMM0, SOLVING A SYMMETRIZED PROBLEM
! A+^(1/2) A- A+^(1/2) [A+^(-1/2)] Y = om^2 [A+^(-1/2)] Y IS SOLVED
!
! ON EXIT ABPLUS CONTAINS ABPLUS^(1/2) !!!
!
implicit none
integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IndN(2,NDimX)
double precision,intent(in) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(in) :: CICoef(NBasis)

double precision,intent(out) :: EigY(NDimX*NDimX),EigX(NDimX*NDimX),Eig(NDimX)

integer :: i, k, ip, iq
double precision :: val, X, Y
!double precision,parameter :: ThreshCIPSI = 1.0d-12
double precision,parameter :: ThreshCIPSI = 0.0d0

call ERPASYMM0(EigY,EigX,Eig,ABPLUS,ABMIN,NDimX)

do k=1,NDimX
   do i=1,NDimX

      ip = IndN(1,i)
      iq = IndN(2,i)

      val = CiCoef(ip) + CiCoef(iq)
      if(abs(val).gt.ThreshCIPSI) X = EigX((k-1)*NDimX+i) / val

      val = CiCoef(ip) - CiCoef(iq)
      if(abs(val).gt.ThreshCIPSI) Y = EigY((k-1)*NDimX+i) / val

      EigX((k-1)*NDimX+i) = 0.5d0 * (X - Y)
      EigY((k-1)*NDimX+i) = 0.5d0 * (X + Y)

   enddo
enddo

end subroutine ERPASYMMXY_CIPSI

subroutine trace_AB(A,B,ndim,res)
implicit none

integer,intent(in) :: ndim
double precision,intent(in) :: A(ndim,ndim),B(ndim,ndim)
double precision,intent(out) :: res

integer :: i,j

res = 0d0
do j=1,ndim
   do i=1,ndim
      res = res + A(i,j)*B(j,i)
   enddo
enddo

end subroutine trace_AB

