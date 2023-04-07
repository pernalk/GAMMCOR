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

!call print_sqmat_cols(Coeff,NAO)

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

