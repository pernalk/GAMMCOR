module sapt_main
use types
use systemdef
use tran
!use sorter

implicit none

contains

subroutine sapt_driver(Flags,SAPT)
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT
integer :: i

! TEMPORARY - JOBTYPE_2
! AC0
 Flags%IFlAC  = 1
 Flags%IFlSnd = 1

 write(LOUT,'()')
 write(LOUT,'(1x,a)') 'STARTING SAPT CALCULATIONS'
 write(LOUT,'(8a10)') ('**********',i=1,8)

 call sapt_interface(Flags,SAPT)

 call print_warn(SAPT)
 call free_sapt(SAPT)

end subroutine sapt_driver

subroutine sapt_interface(Flags,SAPT) 
implicit none

type(FlagsData) :: Flags
type(SaptData) :: SAPT

integer :: NBasis,NDim,NInte1,NInte2
integer :: NCMOt, NOrbt, NBasist 
integer :: NSym, NOrb
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: Ha(:),Hb(:)
double precision,allocatable :: Va(:),Vb(:),S(:)
double precision,allocatable :: Ca(:),Cb(:)
double precision :: potnucA,potnucB
integer :: ione,iorb,isiri,i,j
logical :: exsiri
double precision :: tmp
!integer :: K,LL,NN,NLine
integer :: p,q
double precision,allocatable :: work3(:,:)
character(8) :: dupa
integer :: tmp1
double precision ::  potnuc,emy,eactiv,emcscf



! read basis info
! only dimer basis
 call basinfo(NBasis,'SIRIUS_A.RST')
 if(NBasis==0) then
   NBasis = SAPT%monA%NBasis
 endif
! set dimensions
  NDim = NBasis**2
  NInte1 = NBasis*(NBasis+1)/2

 allocate(work1(NInte1),work2(NBasis**2))
 allocate(Ha(NBasis**2),Hb(NBasis**2))
 allocate(Va(NBasis**2),Vb(NBasis**2),S(NBasis**2))
! read 1-electron integrals
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
 open(newunit=ione,file='AOONEINT_A',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NOrb,potnucA
 
 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call triang_to_sq(work1,Ha,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call triang_to_sq(work1,work2,NBasis)
 Va(:) = Ha - work2

 call readlabel(ione,'OVERLAP ')
 call readoneint(ione,work1)
 call triang_to_sq(work1,S,NBasis)

 call print_sqmat(S,NBasis)
 call print_diag(Va,NBasis)

 tmp = 0d0
 do i=1,NDim
    tmp = tmp + S(i)**2
 enddo
 write(*,*) tmp

 close(ione)
 ! square form
 call writeoneint('ONEEL_A',NDim,S,Va,Ha)
 
!BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

 open(newunit=ione,file='AOONEINT_B',access='sequential',&
      form='unformatted',status='old')
 read(ione) 
 read(ione) NSym,NOrb,potnucB

 call readlabel(ione,'ONEHAMIL')
 call readoneint(ione,work1)
 call triang_to_sq(work1,Hb,NBasis)

 call readlabel(ione,'KINETINT')
 call readoneint(ione,work1)
 call triang_to_sq(work1,work2,NBasis)
 Vb(:) = Hb - work2

 work1 = 0
 call readlabel(ione,'OVERLAP ')
 call readoneint(ione,work1)
 call triang_to_sq(work1,S,NBasis)

 tmp = 0d0
 do i=1,NBasis**2
    tmp = tmp + S(i)**2
 enddo
 write(*,*) tmp


 write(*,*) potnucB,NOrb,NBasis
! call print_sqmat(S,NBasis)

 close(ione)
 ! square form
 call writeoneint('ONEEL_B',NDim,S,Vb,Hb)

! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
 ! read orbitals, coefficient, occupancies
 
 inquire(file='SIRIFC_A',EXIST=exsiri)
 if(exsiri) then
    open(newunit=isiri,file='SIRIFC_A',status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')
    call readlabel(isiri,'TRCCINT ')
    read(isiri) NSym,NOrbt,NBasist,NCMOt
 !   write(*,*) NSym,NOrbt,NBasist,NCMOt
 !   write(*,*) NOrb, 'NOrb'
    SAPT%monA%NOrb = NOrbt
 else
    SAPT%monA%NOrb = NOrb
    NBasist = NBasis
    NCMOt = NOrb*NBasis
 endif 

    rewind(isiri)
    read (isiri)
    read (isiri) potnuc,emy,eactiv,emcscf
    write(*,*)   potnuc,emy,eactiv,emcscf


 if(Flags%ICASSCF==1) then
    call readmulti(NBasis,SAPT%monA,exsiri,isiri,'occupations.dat','SIRIUS_A.RST')
 elseif(Flags%IGVB==1) then 
    call readgvb(SAPT%monA,NBasis,'coeff_A.dat')
 endif

 if(exsiri) close(isiri) 

 call print_occ(NBasis,SAPT,Flags%ICASSCF)

 ! norb.leq.nbas, orbitals mays be deleted due to linear
 ! dependecies in large basis sets
 ! ncmot = norb*nbas
 allocate(Ca(NCMOt))
 allocate(work3(NBasis,NBasis))

 call read_mo(Ca,NOrbt,NBasist,'SIRIUS_A.RST','DALTON_A.MOPUN')
 if(SAPT%IPrint.ne.0) call print_mo(Ca,NOrbt)

! HERE!!!!


! allocate(work3(NBasis,NBasis))

! k = 0
! work3 = 0
! do p=1,NBasis
!    do q=1,NBasis
!       k = k + 1
!       work3(p,q) = Ca(k)
!    enddo
! enddo
 
 
 deallocate(work1,work2)

!!!! check stuff!
!  S = 0
!  Va = 0
!  Ha = 0
!  open(newunit=ione,file='ONEEL_A',access='SEQUENTIAL', &
!        form='UNFORMATTED',status='OLD')
! 
!  read(ione) dupa, S 
!  read(ione) dupa, Va
!  read(ione) dupa, Ha
!
!  write(*,*) 'Ha',Ha(1:NBasis) 
!  write(*,*) 'S',S(1:NBasis) 
!  print*, dupa
! 
!  close(ione) 

! don't forget to write Ca, Cb to file!!!

 deallocate(Ha,Hb,Va,Vb,S)
 deallocate(Ca)

end subroutine sapt_interface

subroutine read_mo(cmo,norb,nbas,nsiri,nmopun)
! in SAPT orbitals kept in AOMO order!
implicit none

integer :: iunit,norb,nbas
double precision :: cmo(norb,nbas)
character(*) :: nsiri,nmopun
logical :: isiri
character(60) :: line
integer :: i,j
double precision :: natocc(10)

inquire(file=nsiri,EXIST=isiri)

if(isiri) then

   open(newunit=iunit,file=nsiri,status='OLD', &
         access='SEQUENTIAL',form='UNFORMATTED')

   call readlabel(iunit,'NEWORB  ')
   read(iunit) cmo

   write(LOUT,'(1x,a)') 'ORBITALS READ FROM SIRIUS.RST' 
else

   open(newunit=iunit,file=nmopun, &
        form='FORMATTED',status='OLD')
   read(iunit,'(a60)') line
   do j=1,norb     
      read(iunit,'(4f18.14)') (cmo(i,j),i=1,nbas)
   enddo
!   print*, line

   write(LOUT,'(1x,a)') 'ORBITALS READ FROM DALTON.MOPUN' 
endif

close(iunit)

end subroutine read_mo

subroutine print_mo(cmo,n)
implicit none

integer,intent(in) :: n
double precision,intent(in) :: cmo(n,n) 
integer :: i,j,ll,nn
integer :: nline

 do i=1,n
    write(LOUT,'(1x,i3)') i
    write(LOUT,'(10f10.6)') cmo(:,i)
    write(LOUT,'()')
 enddo

! nLine=n/10
! if(nLine*10-n.Ne.0)nLine=nLine+1
! do i=1,n
!    write(*,'(i3)') i
!
!    do ll=0,nLine-1
!       nn=n-10*ll
!       if(nn.le.10) then
!          write(LOUT,'(10f10.6)') (cmo(i,j),j=10*ll+1,n)
!       else
!          write(LOUT,'(10f10.6)') (cmo(i,j),j=10*ll+1,10*(ll+1))
!       endIf
!    enddo
!    write(LOUT,'()')
! enddo
 

end subroutine print_mo

subroutine writeoneint(mon,ndim,S,V,H)
implicit none

integer :: ione,ndim
character(*) :: mon
double precision,dimension(ndim) :: S, V, H

 open(newunit=ione,file=mon,form='unformatted')
 write(ione) 'OVERLAP ', S 
 write(ione) 'POTENTAL', V
 write(ione) 'ONEHAMIL', H
 close(ione)

 write(LOUT,'(1x,a)') 'One-electron integrals written to file: '//mon
 write(LOUT,'()')

end subroutine writeoneint


subroutine readgvb(mon,n,cfile)
implicit none

type(SystemBlock) :: mon
integer :: n
character(*) :: cfile
integer :: iunit
integer :: NAct, NIActive
integer :: i,j
!double precision,allocatable :: CICoef(:), Occ(:)
!integer,allocatable :: IGem(:)

open(newunit=iunit,file=cfile,form='FORMATTED',Status='OLD')
read(iunit,'(i5)') mon%NAct

mon%INAct = mon%NELE - mon%NAct

write(*,*) mon%NELE, mon%NAct, mon%INAct

allocate(mon%CICoef(n),mon%IGem(n),mon%Occ(n))
mon%CICoef = 0d0

!!!HERE
do i=1,mon%INAct
   mon%CICoef(i) = 1.0d0
   mon%IGem(i) = i
enddo

read(iunit,*) (mon%CICoef(i+mon%INAct),i=1,2*mon%NAct)

do i=mon%INAct+1,mon%NELE
   mon%IGem(i) = i 
   mon%IGem(mon%NELE+i-mon%INAct) = i
enddo
mon%NGem = mon%NELE + 1

do i=1,n
   if(mon%CICoef(i).eq.0d0) mon%IGem(i) = mon%NGem
   mon%Occ(i) = mon%CICoef(i)**2
enddo

close(iunit)

end subroutine readgvb

subroutine readmulti(nbas,mon,exsiri,isiri,occfile,occsir)
implicit none 

type(SystemBlock) :: mon
logical :: exsiri, ioccsir
integer :: isiri, nbas
character(*) :: occfile, occsir 
logical :: iocc
integer :: iunit,i
integer :: NAct, INAct
double precision :: Occ(nbas), sum1, sum2
double precision :: potnuc, emy, eactiv, emcscf
integer :: istate, ispin, nactel, lsym
integer :: nisht, nasht, nocct, norbt, nbast, nconf, nwopt, nwoph

 allocate(mon%CICoef(nbas),mon%IGem(nbas),mon%Occ(nbas))
 if(exsiri) then

    rewind(isiri) 
    read (isiri) 
    read (isiri) potnuc,emy,eactiv,emcscf, &
                 istate,ispin,nactel,lsym
    read (isiri) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
   
!    print*,    potnuc,emy,eactiv,emcscf, &
!               istate,ispin,nactel,lsym
!    write (*,*) nisht,nasht,nocct,norbt,nbast !,nconf,nwopt,nwoph
   
    mon%NAct  = nasht
    mon%INAct = nisht

    if(nbast.ne.nbas) then
      write(LOUT,'(1x,a)') 'WARNING! NBasis FROM SIRIFC DOES NOT MATCH!'
      write(LOUT,'(1x,a,i5,1x,a,i5)') 'NBasis: ',nbas, 'SIRIFC: ', nbast 
      write(LOUT,'()')
      mon%IWarn = mon%IWarn + 1 
    endif

   inquire(file=occsir,EXIST=ioccsir)
   if(ioccsir) then
      mon%Occ = 0d0
      open(newunit=iunit,file=occsir,status='OLD', &
           access='SEQUENTIAL',form='UNFORMATTED')

      call readlabel(iunit,'NATOCC  ')
      read(iunit) mon%Occ(1:mon%NAct+mon%INAct) 

      sum1 = 0d0
      do i=1,mon%INAct+mon%NAct
         mon%Occ(i) = mon%Occ(i)/2d0
         sum1 = sum1 + mon%Occ(i)
      enddo
      mon%SumOcc = sum1
      close(iunit)
   endif

 endif
! occupations.dat 
 inquire(file=occfile,EXIST=iocc)
 if(iocc) then
 
    Occ = 0d0
    open(newunit=iunit,file=occfile,form='FORMATTED',status='OLD') 
 
    read(iunit,*) INAct, NAct
    INAct = INAct/2
    read(iunit,*) (Occ(i),i=1,INAct+NAct)
    sum2 = 0d0
    do i=1,INAct+NAct
       Occ(i) = Occ(i)/2d0
       sum2 = sum2 + Occ(i)
    enddo

    if(.not.ioccsir) then
       if(Abs(sum2-mon%XELE).gt.1.0d-8) then
          write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'  
          write(LOUT,'(1x,a,1x,f10.6,5x,a,i3)') 'SUM(OCC): ', sum2, 'MONOMER: ', mon%Monomer
          write(LOUT,'(1x,a)') 'CHECK occupations.dat!'  
          stop
       endif
       mon%INAct = INAct
       mon%NAct  = NAct
       mon%Occ   = Occ
       ! SIRIUS.RST not there  
       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occfile 
       write(LOUT,'()')

    else !compare SIRIUS.RST and occupations.dat 
       if(Abs(sum1-mon%XELE).gt.1.0d-8) then

          if(Abs(sum2-mon%XELE).gt.1.0d-8) then
             write(LOUT,'(1x,a,1x,i3)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE! MONOMER: ', mon%Monomer  
             write(LOUT,'(1x,a,1x,f10.6,5x,a,4x,f10.6)') 'OCC(SIRIUS): ', sum1,&
                          'OCC(occupations.dat)', sum2
             stop
          endif

          mon%INAct = INAct
          mon%NAct  = NAct
          mon%Occ   = Occ
          print*, 'sum1 zle' 
       endif

       ! both files correct     
       if(any(abs(Occ-mon%Occ).gt.1.d-9)) then
          write(LOUT,'(1x,a)') 'WARNING! DIFFERENT OCCUPANCIES IN SIRIUS.RST&
                & AND occupations.dat!'
          write(LOUT,'(1x,a)') 'OCCUPANCIES READ FROM occupations.dat!'
          write(LOUT,'()')
          mon%IWarn = mon%IWarn + 1
          mon%INAct = INAct
          mon%NAct  = NAct
          mon%Occ   = Occ
          mon%SumOcc = sum2
       endif

       write(LOUT,'(1x,a,i2,a)') 'OCCUPANCIES FOR MONOMER',mon%Monomer,' READ FROM '// occsir 
       write(LOUT,'()')

    endif
 
    close(iunit)

 elseif(ioccsir) then
    if(Abs(sum1-mon%XELE).gt.1.0d-8) then
       write(LOUT,'(1x,a)') 'ERROR! OCCUPANCIES DO NOT SUM TO NELE!'  
       write(LOUT,*) 'Occ: ', sum1
       write(LOUT,'(1x,a)') 'CHECK DALTON CALCULATIONS!'  
       stop
    endif 

 elseif(.not.ioccsir) then
     write(LOUT,'(1x,a)') 'ERROR! CANNOT READ OCCUPANCIES!'  
     stop
 endif

 if(mon%INAct==0) then
    mon%NGem = 2
    mon%IGem(1:mon%NAct+mon%INAct) = 1
    mon%IGem(mon%NAct+mon%INAct+1:nbas) = 2
 else
    mon%NGem = 3
    mon%IGem(1:mon%INAct) = 1
    mon%IGem(mon%INAct+1:mon%INAct+mon%NAct) = 2
    mon%IGem(mon%INAct+mon%NAct+1:nbas) = 3
 endif

end subroutine readmulti

subroutine print_occ(nbas,SAPT,ICASSCF)
implicit none
!!! HERE : Change to A/B monomers!
type(SaptData) :: SAPT
integer :: nbas, ICASSCF
integer :: i

 associate(A => SAPT%monA, B => SAPT%monB)
 if(ICASSCF==0) then
   write(LOUT,'(2X,"Orb",3X,"Occupancy",7x,"CICoef",7X,"Gem")')
   do i=1,nbas
      write(6,'(X,i3,2e16.6,i6)') i,A%Occ(i),A%CICoef(i),A%IGem(i)
   enddo
   write(LOUT,'()') 

 else

   write(LOUT,'(1x,a,1x,i3,i3)') 'NO OF CAS INACTIVE AND ACTIVE ORBITALS:',A%INAct, A%NAct
   write(LOUT,'(1x,a,3x,a,4x,a)') 'CASSCF', 'Occupancy', 'Gem'
   do i=1,nbas
      write(LOUT,'(1x,i3,e16.6,i6)') i, A%Occ(i), A%IGem(i)
   enddo
   write(LOUT,'(2x,a,e16.6)') 'SUM OF OCCUPANCIES: ', A%SumOcc
 endif
 end associate

end subroutine print_occ

subroutine print_warn(SAPT)
implicit none

type(SaptData) :: SAPT
integer :: cnt,i

cnt = SAPT%monA%IWarn+SAPT%monB%IWarn 
if(cnt.gt.0) then
    write(LOUT,'(1x,a,i2,1x,a)') 'SAPT: CHECK OUTPUT FOR',cnt,'WARNINGS!'
    write(LOUT,'(8a10)') ('**********',i=1,8)
endif

end subroutine print_warn

subroutine free_sapt(SAPT)
implicit none

type(SaptData) :: SAPT

deallocate(SAPT%monA%CICoef, SAPT%monA%IGem, SAPT%monA%Occ)

end subroutine free_sapt

end module sapt_main

