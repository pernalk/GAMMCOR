module read_external
!
! common subroutines for reading/writing external
! files (RDMs, MOs, ...) from Dalton, Molpro, ...

use print_units

contains

subroutine read_2rdm_molpro(twordm,nost,nosym,noms2,infile,iwarn,nact)
!
! Purpose: read active part of spin-summed 2RDM from Molpro
!          in MO representation
!
implicit none

integer,intent(in) :: nact,nost,nosym,noms2
integer,intent(inout) :: iwarn
character(*),intent(in) :: infile
double precision,intent(out) :: twordm(nact**2*(nact**2+1)/2)

integer :: iunit,ios,ist,ic1d,TrSq
integer :: istsym,isym,nstate,ims2,nstsym
integer :: i,j,k,l,ij,kl,ik,jl,ijkl,ikjl,lend
double precision,allocatable :: work(:)
logical :: scanfile,scanms2
character(8) :: label

 !check if any states declared in input
 scanfile = .false.
 scanms2  = .false.
 if(nosym>0)  scanfile = .true.
 if(noms2>=0) scanms2  = .true.

 TrSq = nact**2*(nact**2+1)/2

 allocate(work(TrSq))
 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 2RDM     not found!'
              stop
           endif
           if(label=='2RDM    ') then
              !read(iunit) ic1d,nstate
              read(iunit) ic1d,nstsym
              if(scanfile) then
                 do istsym=1,nstsym
                    read(iunit) isym,nstate,ims2
                    !print*, 'isym,nstate,ims2',isym,nstate,ims2
                    do i=1,nstate
                       read(iunit) ist
                       read(iunit) work(1:TrSq)
                       if(scanms2) then
                          if(ist==nost.and.isym==nosym.and.ims2==noms2) exit fileloop
                       else
                          if(ist==nost.and.isym==nosym) exit fileloop
                       endif
                    enddo
                 enddo
                 if(scanms2) then
                    write(LOUT,'(1x,a,i2,a,i1,a,i1,a)') 'ERROR!!! 2RDM FOR STATE',&
                                & nost,'.',nosym,' AND MS2 = ',noms2,' NOT PRESENT IN 2RDM FILE!'
                 else
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'ERROR!!! 2RDM FOR STATE',&
                                & nost,'.',nosym,' NOT PRESENT IN 2RDM FILE!'
                 endif
                 stop
              else
                 read(iunit) isym,nstate
                 read(iunit) ist
                 read(iunit) work(1:TrSq)
                 if(isym/=1.or.ist/=1) then
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'WARNING! 2RDM FOR STATE',&
                                & ist,'.',isym,' WILL BE USED IN CALCULATIONS!'
                    iwarn = iwarn + 1
                 endif
                 exit fileloop
              endif
           endif
         enddo fileloop

 close(iunit)

 ijkl = 0
 do i=1,nact
    do j=1,nact
       do k=1,i
          lend = nact
          if(k==i) lend = j
          do l=1,lend
             ijkl = ijkl + 1

             ik = (i-1)*nact + k
             jl = (j-1)*nact + l
             ikjl = max(ik,jl)*(max(ik,jl)-1)/2+min(ik,jl)
             twordm(ikjl) = work(ijkl)*0.5d0

             ik = (k-1)*nact + i
             jl = (l-1)*nact + j
             ikjl = max(ik,jl)*(max(ik,jl)-1)/2+min(ik,jl)
             twordm(ikjl) = work(ijkl)*0.5d0

          enddo
       enddo
    enddo
 enddo

 deallocate(work)

!!! HERE : in SAPT writing to rdm2.dat will be necessary
!!!! better loop!!!
! ijkl = 0
! open(newunit=iunit,file=outfile,form='unformatted')
! do i=1,nact
!    do j=1,nact
!       ij = (i-1)*nact+j
!       do k=1,nact
!          do l=1,nact
!             kl = (k-1)*nact+l
!             if(ij>=kl) then
!               ijkl = ijkl + 1
!               write(iunit,'(4I4,F19.12)') k,i,l,j,twordm(ijkl)
!             endif
!          enddo
!       enddo
!    enddo
! enddo
! close(iunit)

end subroutine read_2rdm_molpro

subroutine read_1rdm_molpro(onerdm,nost,nosym,noms2,infile,iwarn,nbasis)
!
! Purpose: read 1RDM (ntriang) file from Molpro
!
implicit none

integer,intent(in) :: nbasis,nost,nosym,noms2
integer,intent(inout) :: iwarn
character(*),intent(in) :: infile
double precision,intent(out) :: onerdm(nbasis*(nbasis+1)/2)

integer :: iunit,ios,ist,isym,nact,nact2,nstate,ims2,nstsym
integer :: i,j,ij,idx,istsym
double precision,allocatable :: work(:)
logical :: scanfile,scanms2
character(8) :: label

 !check if any states declared in input
 scanfile = .false.
 scanms2  = .false.
 if(nosym>0)  scanfile = .true.
 if(noms2>=0) scanms2  = .true.

 allocate(work(NBasis**2))
 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1RDM     not found!'
              stop
           endif
           if(label=='1RDM    ') then
              read(iunit) nact,nact2,nstsym
              !print*, nact2,nstate
              if(scanfile) then
                 do istsym=1,nstsym
                    read(iunit) isym,nstate,ims2
                    do i=1,nstate
                       read(iunit) ist
                       read(iunit) work(1:nact2)
                       if(scanms2) then
                          !print*, 'ist,isym,ims2',ist,isym,ims2
                          if(ist==nost.and.isym==nosym.and.ims2==noms2) exit fileloop
                       else
                          if(ist==nost.and.isym==nosym) exit fileloop
                       endif
                    enddo
                 enddo
                 if(scanms2) then
                    write(LOUT,'(1x,a,i2,a,i1,a,i1,a)') 'ERROR!!! 1RDM FOR STATE',&
                                & nost,'.',nosym,' AND MS2 = ',noms2,' NOT PRESENT IN 1RDM FILE!'
                 else
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'ERROR!!! 1RDM FOR STATE',&
                                & nost,'.',nosym,' NOT PRESENT IN 1RDM FILE!'
                 endif
                 stop
              else
                 read(iunit) isym,nstate
                 read(iunit) ist
                 read(iunit) work(1:nact2)
                 if(isym/=1.or.ist/=1) then
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'WARNING! 1RDM FOR STATE',&
                                & ist,'.',isym,' WILL BE USED IN CALCULATIONS!'
                    write(LOUT,'()')
                    iwarn = iwarn + 1
                 endif
                 exit fileloop
              endif
           endif
         enddo fileloop

 close(iunit)

 !call sq_to_triang(work,onerdm,nact)
 onerdm = 0
 idx = 0
 do j=1,nact
    do i=1,j
       idx = idx + 1
       ij = i + (j-1)*nact
       onerdm(idx) = work(ij)*0.5d0
    enddo
 enddo

 deallocate(work)

end subroutine read_1rdm_molpro

subroutine read_1trdm_molpro(onerdm,stbrIn,stketIn,infile,nbasis)
implicit none

integer,intent(in) :: nbasis,stbrIn,stketIn
character(*),intent(in) :: infile
double precision,intent(out) :: onerdm(nbasis,nbasis)
!double precision,intent(out) :: onerdm(nbasis*(nbasis+1)/2)

integer :: iunit,ios,ist,isym,nact,nact2,nstate,nstsym
integer :: stbr,stket,istbr,istket
integer :: i,j,ij,idx,istsym
double precision,allocatable :: work(:)
logical :: scanfile
character(8) :: label

 !order states
 if(stbrIn.gt.stketIn) then
    stbr  = stketIn
    stket = stbrIn
 else
    stbr =  stbrIn
    stket = stketIn
 endif

 write(LOUT,'(1x,a,i1,a,i1,a)') &
       'Reading <',stbr,'|',stket,'> 1-TRDM...'

 !check if any states declared in input
 scanfile = .false.
 if(stket>0) scanfile = .true.

 allocate(work(NBasis**2))
 work=0
 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1TRMD    not found!'
              stop
           endif
           if(label=='1TRDM    ') then
              read(iunit) nact,nact2,nstsym
              !print*, nact2,nstate
              if(scanfile) then
                 do istsym=1,nstsym
                    read(iunit) isym,nstate
                    do j=1,nstate
                       do i=1,j-1
                          read(iunit) istbr,istket
                          read(iunit) work(1:nact2)
                          if(istbr==stbr.and.istket==stket) exit fileloop
                       enddo
                    enddo
                 enddo
                 write(LOUT,'(1x,a,i2,a,i1,a)') 'ERROR!!! 1TRDM FOR STATE',&
                             & stbr,'.',stket,' NOT PRESENT IN 1RDM FILE!'
                 stop
              else
                 read(iunit) isym,nstate
                 read(iunit) istbr,istket
                 read(iunit) work(1:nact2)
                 if(isym/=1.or.ist/=1) then
                    write(LOUT,'(1x,a,i2,a,i1,a)') 'WARNING! 1TRDM FOR STATE',&
                                & istbr,'.',istket,' WILL BE USED IN CALCULATIONS!'
                    write(LOUT,'()')
                 endif
                 exit fileloop
              endif
           endif
         enddo fileloop

 close(iunit)

 !call sq_to_triang(work,onerdm,nact)
 onerdm = 0
 do j=1,nact
    do i=1,nact
       ij = i + (j-1)*nact
       onerdm(i,j) = work(ij)*0.5d0
    enddo
 enddo

 deallocate(work)

end subroutine read_1trdm_molpro

subroutine readoneint_molpro(mone,infile,text,expand,ntr)
!
! reads one-el integrals from a molpro file
! if expand=.true. destroys symmetry in mone
!
implicit none

integer,intent(in) :: ntr
logical,intent(in) :: expand
character(*),intent(in) :: infile,text
character(8) :: label
integer :: iunit,ios
integer :: nsym,nbas(8),offs(8),ntdg
integer :: i,j,idx,irep,ioff,ii
double precision,intent(out) :: mone(ntr)
double precision :: tmp(ntr)

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')
 ntdg = 0
 rewind(iunit)
 read(iunit)
 read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
 read(iunit)
 do irep=1,nsym
    ntdg = ntdg + nbas(irep)*(nbas(irep)+1)/2
 enddo

 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL '//text//' not found!'
      stop
   endif
   if(label==text) then
      read(iunit) tmp(1:ntdg)
      exit
   endif
 enddo

 close(iunit)

 if(expand) then
 ! expand to a large triangle
    mone = 0
    ii   = 0
    idx  = 0
    ioff = 0
    do irep=1,nsym
      ioff = offs(irep)
      do j=1,nbas(irep)
        do i=1,j
           idx = idx + 1
           ii = (ioff+j)*(ioff+j-1)/2 + ioff+i
           mone(ii) = tmp(idx)
        enddo
      enddo
    enddo
 else
 ! return ntdg
    mone = 0
    mone(1:ntdg) = tmp(1:ntdg)
 endif

end subroutine readoneint_molpro

subroutine read_mo_molpro(cmo,infile,text,nbasis)
implicit none

integer,intent(in) :: nbasis
character(*),intent(in) :: infile,text
double precision,intent(out) :: cmo(nbasis,nbasis)
character(8) :: label
integer :: iunit,ios
integer :: nsym,nbas(8),offs(8),ncmot
integer :: i,j,idx,irep,ioff
double precision ::tmp(nbasis**2)

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 !rewind(iunit)
 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL '//text//' not found!'
      stop
   endif
   if(label==text) then
      read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
      ncmot = sum(nbas(1:nsym)**2)
      !print*, nsym,nbas(1:nsym),offs(1:nsym)
      !print*, ncmot
      read(iunit) tmp(1:ncmot)
      exit
   endif
 enddo

!! HERE: should there be a norb-offset?
 cmo = 0
 idx = 0
 do irep=1,nsym
    ioff = offs(irep)
    do j=1,nbas(irep)
       do i=1,nbas(irep)
          idx = idx + 1
          cmo(ioff+i,ioff+j) = tmp(idx)
       enddo
    enddo
 enddo

  !write(LOUT,*) 'test print'
  !do i=1,NBasis
  !   write(*,'(14f11.6)') (cmo(i,j),j=1,nbasis)
  !end do

 close(iunit)

end subroutine read_mo_molpro

subroutine readgridmolpro(iunit,text,npt,r,wt)
!
! Purpose: read grid points and weights from Molpro
!
implicit none

integer :: iunit
integer :: ios
integer :: npt
double precision :: r(:,:),wt(:)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'
     stop
  endif
  if(label==text) then

     read(iunit) npt
     read(iunit) r(1:3,1:npt)
     read(iunit) wt(1:npt)
     exit
  endif
enddo

end subroutine readgridmolpro

subroutine readorbsmolpro(iunit,text,mapinv,orbval)
!
! Purpose: read orbitals vals on a grid
!
implicit none

integer :: iunit
integer :: ios
integer :: npt, ndiff, ntg
double precision :: mapinv(:),orbval(:,:,:)
character(8) :: text, label

rewind(iunit)
do
  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! LABEL '//text//' not found!'
     stop
  endif
  if(label==text) then

     read(iunit) npt,ndiff,ntg
     read(iunit) orbval(1:npt,1:ndiff,1:ntg)
     read(iunit) mapinv(1:ntg)
     exit
  endif
enddo

end subroutine readorbsmolpro



subroutine read_nact_molpro(nact,infile)
!
! Purpose: reads the number of active electrons
!          which is kept in 2RDM file
!
implicit none

character(*),intent(in) :: infile
integer :: iunit,ios,ist,nact,nact2,nstate
integer :: i,j,ij,idx
character(8) :: label

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1RDM     not found!'
              stop
           endif
           if(label=='1RDM    ') then
              read(iunit) nact,nact2,nstate
              exit fileloop
!              print*, nact,nact2,nstate
           endif
         enddo fileloop

 close(iunit)

end subroutine read_nact_molpro

subroutine read_NoSt_molpro(NoSt,infile)
!
! when State variable not given in input, read the state number (NoSt)
! for which calculations will be performed
!
implicit none

integer,intent(inout) :: NoSt
character(*),intent(in) :: infile

integer :: iunit,ios,isym,nstate
character(8) :: label

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 fileloop: do
           read(iunit,iostat=ios) label
           if(ios<0) then
              write(LOUT,*) 'ERROR!!! LABEL 1RDM     not found!'
              stop
           endif
           if(label=='1RDM    ') then
              read(iunit)
              read(iunit) isym,nstate
              read(iunit) NoSt
              exit fileloop
           endif
         enddo fileloop

 close(iunit)

end subroutine read_NoSt_molpro

subroutine read_dip_molpro(matdx,matdy,matdz,infile,nbasis)
!
! Purpose: read dipole integrals and unpack withour symmetry
!
implicit none

!type(SystemBlock) :: mon

integer,intent(in) :: nbasis
character(*),intent(in) :: infile
double precision :: matdx(nbasis,nbasis),matdy(nbasis,nbasis),matdz(nbasis,nbasis)
double precision,allocatable :: dipx(:),dipy(:),dipz(:)

integer :: iunit,ios,ictrl
integer :: nsym,nbas(8),offs(8),ntqg
integer :: isx,isy,isz,isym_p,isym_q
integer :: i,j,ij,irep
logical :: scanfile
character(8) :: label

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 ntqg = 0
 rewind(iunit)
 read(iunit)
 read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
 do irep=1,nsym
    ntqg = ntqg + nbas(irep)**2
 enddo

 allocate(dipx(ntqg),dipy(ntqg),dipz(ntqg))

 ictrl = 0
 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL DIPMOMX not found!'
      stop
   endif
   if(label=='DIPMOMX ') then
      ictrl = ictrl + 1
      read(iunit) isx
      read(iunit) dipx(1:ntqg)
   elseif(label=='DIPMOMY ') then
      ictrl = ictrl + 1
      read(iunit) isy
      read(iunit) dipy(1:ntqg)
   elseif(label=='DIPMOMZ ') then
      ictrl = ictrl + 1
      read(iunit) isz
      read(iunit) dipz(1:ntqg)
   endif
   if(ictrl==3) exit
 enddo

 close(iunit)

 !print*, 'dmx',dipx(1:ntqg)
 !print*,''
 !print*, 'dmy',dipy(1:ntqg)
 !print*,''
 !print*, 'dmz',dipz(1:ntqg)

 !allocate(mon%dipm(3,nbasis,nbasis))
 !mon%dipm=0d0
 matdx=0
 matdy=0
 matdz=0

 ij = 0
 do isym_p=1,nsym
    isym_q = ieor(isx-1,isym_p-1)+1
    if(nbas(isym_p)>0.and.nbas(isym_q)>0) then
       !print*, isym_p,isym_q
       !print*, 'offs_p(q)',offs(isym_p),offs(isym_q)
       !print*, 'nbas_p(q)',nbas(isym_p),nbas(isym_q)
       do j=1,nbas(isym_q)
          do i=1,nbas(isym_p)
             ij = ij + 1
             !print*, i,j,dipx(ij)
             !mon%dipm(1,offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipx(ij)
             matdx(offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipx(ij)
          enddo
       enddo
    endif
 enddo

 ij = 0
 do isym_p=1,nsym
    isym_q = ieor(isy-1,isym_p-1)+1
    if(nbas(isym_p)>0.and.nbas(isym_q)>0) then
       !print*, isym_p,isym_q
       do j=1,nbas(isym_q)
          do i=1,nbas(isym_p)
             ij = ij + 1
             !mon%dipm(2,offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipy(ij)
             matdy(offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipy(ij)
          enddo
       enddo
    endif
 enddo

 ij = 0
 do isym_p=1,nsym
    isym_q = ieor(isz-1,isym_p-1)+1
    if(nbas(isym_p)>0.and.nbas(isym_q)>0) then
       do j=1,nbas(isym_q)
          do i=1,nbas(isym_p)
             ij = ij + 1
             !mon%dipm(3,offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipz(ij)
             matdz(offs(isym_p)+i,offs(isym_q)+j) = -1d0*dipz(ij)
          enddo
       enddo
    endif
 enddo

 deallocate(dipx,dipy,dipz)

end subroutine read_dip_molpro

subroutine read_sym_molpro(NSymMO,nsym,nbas,infile,text,nbasis)
!
! Purpose: read irreps of orbitals (NSymMO)
!          and how many orbitals in each irrep (nbas)
!
implicit none

integer,intent(in) :: nbasis
character(*),intent(in) :: infile,text
integer :: NSymMO(nbasis)
character(8) :: label
integer :: iunit,ios
integer :: nsym,nbas(8),offs(8),ncmot
integer :: i,j,idx,irep,ioff
double precision ::tmp(nbasis**2)

 open(newunit=iunit,file=infile,status='OLD', &
      access='SEQUENTIAL',form='UNFORMATTED')

 do
   read(iunit,iostat=ios) label
   if(ios<0) then
      write(6,*) 'ERROR!!! LABEL '//text//' not found!'
      stop
   endif
   if(label==text) then
      read(iunit) nsym,nbas(1:nsym),offs(1:nsym)
      ncmot = sum(nbas(1:nsym)**2)
      read(iunit) tmp(1:ncmot)
      exit
   endif
 enddo

 close(iunit)

 do irep=1,nsym
    ioff = offs(irep)
    do j=1,nbas(irep)
      NSymMO(ioff+j)=irep
    enddo
 enddo

end subroutine read_sym_molpro

subroutine sym_inf_molpro(infile,NumOSym,NSym,nstats,istsy,NStSym,NSymAO,NBasis)
!
! Purpose: reads number of irreps (NSym)
!          number of atomic orbitals in each irrep (NumOSym)
!          number of irreps in SA-CAS (NStSym)
!          number of states in SA-CAS in each irrep (NumStSym)
!          which irrep (a label) in SA-CAS (IStSy)
!
implicit none

character(*) :: infile
integer :: NumOSym(8),NSym,NBasis
integer :: ifile,ios,i,j,k,isym
integer :: NState,NStSym,IOld,INew,ioff
integer :: iclos(8),iact(8),nt(8),ivirt(8),istsy(16),nstats(16),NSymAO(NBasis)
character(8) :: label

 iclos = 0
 ivirt = 0
 iact = 0
 nstats = 0
 istsy = 0
 NumOSym = 0
 open(newunit=ifile,file=infile,access='sequential',&
      form='unformatted',status='old')

 do
   read(ifile,iostat=ios) label
    if(ios<0) then
      write(6,*) 'ERROR!!! LABEL ISORDK   not found!'
      stop
   endif
   if(label=='BASINFO ') then
      read(ifile) NSym
      read(ifile) NStSym
      read(ifile) nstats(1:NStSym)
      read(ifile) istsy(1:NStSym)
      read(ifile) iclos(1:NSym)
      read(ifile) iact(1:NSym)
      read(ifile) NumOSym(1:NSym)
      exit
   endif
 enddo

ioff=0
do i=1,Nsym
    do j=1,NumOSym(i)
      ioff=ioff+1
      NSymAO(ioff)=i
    enddo
 enddo

end subroutine sym_inf_molpro

subroutine create_ind_molpro(infile,NumOSym,IndInt,NSym,NBasis)
!
! Purpose: reads number of atomic orbitals in each irrep (NumOSym)
!          creates index array IndInt(NBasis)
!          used to reorder orbitals (all closed first, then active, then virt)
!
implicit none

character(*) :: infile
integer :: NumOSym(15),IndInt(NBasis),NSym,NBasis
integer :: ifile,ios,i,j,k,isym
integer :: NState,NStSym,IOld,INew
integer :: iclos(8),iact(8),nt(8),ivirt(8),istsy(16),nstats(16)
character(8) :: label

 iclos = 0
 ivirt = 0
 iact = 0
 nstats = 0
 istsy = 0
 NumOSym = 0
 open(newunit=ifile,file=infile,access='sequential',&
      form='unformatted',status='old')

 do
   read(ifile,iostat=ios) label
    if(ios<0) then
      write(6,*) 'ERROR!!! LABEL ISORDK   not found!'
      stop
   endif
   if(label=='BASINFO ') then
      read(ifile) NSym
      read(ifile) NStSym
      read(ifile) nstats(1:NStSym)
      read(ifile) istsy(1:NStSym)
      read(ifile) iclos(1:NSym)
      read(ifile) iact(1:NSym)
      read(ifile) NumOSym(1:NSym)
      exit
   endif
 enddo

 close(ifile)
 ivirt(1:NSym) = NumOSym(1:NSym)-iclos(1:NSym)-iact(1:NSym)

! print*, 'nact',iact(1:NSym)

 if(NSym>1) then
   ! symmetry
   IOld = 0
   INew = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
      enddo
   enddo
   IOld = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
      enddo
   enddo
   IOld = 0
   do isym=1,NSym
      do k=1,iclos(isym)
         IOld = IOld + 1
      enddo
      do j=1,iact(isym)
         IOld = IOld + 1
      enddo
      do i=1,ivirt(isym)
         IOld = IOld + 1
         INew = INew + 1
         IndInt(IOld) = INew
         !print*, IOld,INew
      enddo
   enddo

 else
   ! nosym
   do i=1,NBasis
      IndInt(i) = i
   enddo

 endif

end subroutine create_ind_molpro


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dalton subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readlabel(iunit,text)
!
! Purpose: sets file pointer
!          to first data after text
!          (works with Dalton files)
!
implicit none

integer :: iunit
integer :: ios
character(8) :: text, label(4)

rewind(iunit)
do

  read(iunit,iostat=ios) label
  if(ios<0) then
     write(6,*) 'ERROR!!! Empty section in AOTWOINT!'
     stop
  endif
  if(label(1)=='********') then
     if(label(4)==text) exit
  endif

enddo

end subroutine readlabel

subroutine readoneint_dalton(iunit,ints)
!
! Purpose: read 1-el ints from Dalton
! information are kept in one record in the
! order: buf, ibuf, length
! buf: integrals, ibuf: int number

implicit none

integer :: iunit
double precision :: ints(:)
integer,parameter :: lbuf = 600
double precision :: buf(lbuf)
integer :: ibuf(lbuf)
integer :: length,i

ints=0
do
   read(iunit) buf,ibuf,length
   if(length.lt.0) exit
   do i=1,length
      ints(ibuf(i)) = buf(i)
   enddo
enddo

end subroutine readoneint_dalton

subroutine read_mo_dalton(cmo,nbasis,nsym,nbas,norb,nsiri,nmopun)
!
! reads MOs either from SIRIUS.RST or DALTON.MOPUN
! unpacks symmetry blocks to (NBasis,NOrb) form
! in SAPT orbitals kept in AOMO order!
!
implicit none

integer,intent(in) :: nbasis,nsym,nbas(8),norb(8)
character(*)       :: nsiri,nmopun
double precision   :: cmo(nbasis,nbasis)

integer          :: i,j,idx
integer          :: iunit,irep
integer          :: ncmot
integer          :: off_i,off_j
logical          :: isiri
character(60)    :: line
double precision :: natocc(10)
double precision :: tmp(nbasis**2)

tmp =    0
ncmot    = sum(nbas(1:nsym)*norb(1:nsym))

inquire(file=nsiri,EXIST=isiri)
if(isiri) then

   open(newunit=iunit,file=nsiri,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

   call readlabel(iunit,'NEWORB  ')
   read(iunit) tmp(1:ncmot)

!   cmo = 0
!   off_i = 0
!   off_j = 0
!   idx = 0
!
!   do irep=1,nsym
!      do j=off_j+1,off_j+norb(irep)
!
!         do i=off_i+1,off_i+nbas(irep)
!            idx = idx + 1
!            cmo(i,j) = tmp(idx)
!         enddo
!
!      enddo
!      off_i = off_i + nbas(irep)
!      off_j = off_j + norb(irep)
!   enddo
   write(LOUT,'(1x,a)') 'Orbitals read from '//nsiri

else

   open(newunit=iunit,file=nmopun,form='FORMATTED',status='OLD')
   read(iunit,'(a60)') line
   off_i = 0
   idx   = 0
   do irep=1,nsym
      do j=1,norb(irep)
         read(iunit,'(4f18.14)') (tmp(off_i+i),i=1,nbas(irep))
         off_i = off_i + nbas(irep)
      enddo
   enddo

   write(LOUT,'(1x,a)') 'Orbitals read from '//nmopun
endif

close(iunit)

! unpack orbitals
cmo   = 0
off_i = 0
off_j = 0
idx   = 0

do irep=1,nsym
   do j=off_j+1,off_j+norb(irep)

      do i=off_i+1,off_i+nbas(irep)
         idx = idx + 1
         cmo(i,j) = tmp(idx)
      enddo

   enddo
   off_i = off_i + nbas(irep)
   off_j = off_j + norb(irep)
enddo

! call print_sqmat(cmo,nbasis)

end subroutine read_mo_dalton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Eugene subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readoneint_eugene(mone,enuc,infile,ntr,ipos)
!
! Purpose: read one-el integrals from a eugene's file
!
implicit none

integer,intent(in) :: ntr
integer(8),intent(in) :: ipos
character(*),intent(in) :: infile
double precision,intent(out) :: mone(ntr),enuc

integer :: iunit
integer :: i,j,k,l,ind
double precision :: val

 open(newunit=iunit,file=trim(infile),form='unformatted',access='stream',status='old')
 read(iunit,pos=ipos)
 do
    read(iunit) val,i,j,k,l
    if(i+j==0) then
       enuc = val
       exit
    else
       ind=(max(i,j)*(max(i,j)-1))/2+min(i,j)
       mone(ind) = val
    endif
 enddo

 close(iunit)

end subroutine readoneint_eugene

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General file handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine delfile(filename)
implicit none

character(*) :: filename
integer :: iunit
logical :: iexist

 inquire(file=trim(filename),EXIST=iexist)
 if(iexist) then
   open(newunit=iunit,file=trim(filename),status='OLD')
   close(iunit,status='DELETE')
 else
   write(LOUT,'(1x,a)') 'WARNING! NOT POSSIBLE TO DELETE '// filename //' FILE'
 endif

end subroutine delfile

end module read_external
