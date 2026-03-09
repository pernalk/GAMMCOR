module sapt_exch_os
use types
use tran
use exmisc
use exi
use sapt_utils

implicit none

contains

subroutine e1exch_os(A,B,SAPT)
!
! open-shell unrestricted E1exch (S^infty)
! in the AO basis
! cf. Eq. (8) in https://doi.org/10.1063/1.4758455
!
implicit none

type(FlagsData)   :: Flags
type(SystemBlock) :: A, B
type(SaptData)    :: SAPT

integer :: NAO,NMO
integer :: occa,occb
integer :: i,j,k
integer :: info
integer, allocatable :: ipiv(:)

real(8) :: pex(16),pmix
real(8) :: e1ex

real(8), allocatable :: Sa(:,:),Sb(:,:)
real(8), allocatable :: Va(:,:),Vb(:,:)
real(8), allocatable :: taa(:,:),tbb(:,:)
real(8), allocatable :: tab(:,:),tba(:,:)
real(8), allocatable :: taba(:,:),tabb(:,:),&
                        tbaa(:,:),tbab(:,:)
real(8), allocatable :: pa(:,:),pb(:,:)
real(8), allocatable :: Paa(:,:),Pab(:,:),&
                        Pba(:,:),Pbb(:,:)
real(8), allocatable :: pinva(:,:),pinvb(:,:)
real(8), allocatable :: pinvAAa(:,:),pinvAAb(:,:),&
                        pinvBBa(:,:),pinvBBb(:,:),&
                        pinvABa(:,:),pinvABb(:,:),&
                        pinvBAa(:,:),pinvBAb(:,:)
real(8), allocatable :: haa(:,:),hba(:,:),&
                        hab(:,:),hbb(:,:)
real(8), allocatable :: Kaa(:,:),Kab(:,:),&
                        Kba(:,:),Kbb(:,:)
real(8), allocatable :: Jalpha(:,:),Jbeta(:,:)
real(8), allocatable :: Jtaa(:,:),Ktaa(:,:),&
                        Jtba(:,:),Ktba(:,:),&
                        Jtaba(:,:),Ktaba(:,:)
real(8), allocatable :: Jtab(:,:),Ktab(:,:),&
                        Jtbb(:,:),Ktbb(:,:),&
                        Jtabb(:,:),Ktabb(:,:)
real(8), allocatable :: work(:,:),work2(:,:)

double precision, allocatable :: Waa(:,:),Wab(:,:)
double precision, allocatable :: Wba(:,:),Wbb(:,:)
double precision,external  :: trace

! set dimensions
NAO = SAPT%NAO
NMO = A%NBasis

! occuppied A + occupied B
occa = A%NOa+B%NOa
occb = A%NOb+B%NOb

if (NAO.ne.NMO)  then
   print*, 'NAO =', NAO
   print*, 'NMO =', NMO
   stop "e1exh os!"
endif

allocate(PAa(nao,nao),PAb(nao,nao))
allocate(PBa(nao,nao),PBb(nao,nao))

call get_den(nmo,nmo,A%UMO(:,:,1),A%Uocc(:,1),1d0,PAa)
call get_den(nmo,nmo,A%UMO(:,:,2),A%Uocc(:,2),1d0,PAb)
call get_den(nmo,nmo,B%UMO(:,:,1),B%Uocc(:,1),1d0,PBa)
call get_den(nmo,nmo,B%UMO(:,:,2),B%Uocc(:,2),1d0,PBb)

! h matrices
allocate(Va(nao,nao),Vb(nao,nao))
call get_one_mat('V',Va,A%Monomer,nao)
call get_one_mat('V',Vb,B%Monomer,nao)

allocate(haa(nao,nao),hab(nao,nao))
allocate(Kaa(nao,nao),Kab(nao,nao))
allocate(hba(nao,nao),hbb(nao,nao))
allocate(Kba(nao,nao),Kbb(nao,nao))
allocate(Jalpha(nao,nao),Jbeta(nao,nao))

call make_J1(nao,PAa,jalpha,'AOTWOSORT')
call make_J1(nao,PAb,jbeta ,'AOTWOSORT')
call make_K(nao,PAa,Kaa)
call make_K(nao,PAb,Kab)
haa = Va + Jalpha + JBeta - Kaa
hab = Va + Jalpha + JBeta - Kab

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'hAa =', norm2(hAa)
   print*, 'Kaa =', norm2(Kaa)
   print*, '-----'
   print*, 'hAb =', norm2(hAb)
   print*, 'Kab =', norm2(Kab)
endif

call make_J1(nao,PBa,jalpha,'AOTWOSORT')
call make_J1(nao,PBb,jbeta ,'AOTWOSORT')
call make_K(nao,PBa,Kba)
call make_K(nao,PBb,Kbb)
hba = Vb + Jalpha + JBeta - Kba
hbb = Vb + Jalpha + JBeta - Kbb

if(SAPT%IPrint>=50) then
   print*, '-----'
   print*, 'hBa =', norm2(hBa)
   print*, 'KBa =', norm2(Kba)
   print*, '-----'
   print*, 'hBb =', norm2(hBb)
   print*, 'KBb =', norm2(Kbb)
   print*, '-----'
endif

allocate(work(NAO,NAO))
allocate(Sa(NMO,NMO),Sb(NMO,NMO))

call get_one_mat('S',work,A%Monomer,NAO)
call tran2MO(work,A%UMO(:,:,1),B%UMO(:,:,1),Sa,NMO)
call tran2MO(work,A%UMO(:,:,2),B%UMO(:,:,2),Sb,NMO)

! P matrix
!    P = [S + 1]^{-1}-1
allocate(pa(occa,occa))
pa=0d0
do concurrent(i=1:occa)
   pa(i,i)=1d0
enddo
do j=1,B%NOa
   do i=1,A%NOa
      pa(i,A%NOa+j)=Sa(i,j)
      pa(A%NOa+j,i)=Sa(i,j)
   enddo
enddo
! P beta
allocate(pb(occb,occb))
pb=0d0
do concurrent(i=1:occb)
   pb(i,i)=1d0
enddo
do j=1,B%NOb
   do i=1,A%NOb
      pb(i,A%NOb+j)=Sb(i,j)
      pb(A%NOb+j,i)=Sb(i,j)
   enddo
enddo
! invert alpha
allocate(ipiv(occa))
allocate(pinva(occa,occa))
allocate(pinvAAa(A%NOa,A%NOa),pinvBBa(B%NOa,B%NOa), &
         pinvABa(A%NOa,B%NOa),pinvBAa(B%NOa,A%NOa))
pinva=0d0
do concurrent(i=1:occa)
   pinva(i,i)=1d0
enddo
call dgesv(occa,occa,pa,occa,ipiv,pinva,occa,info)
! -1 alpha
do i=1,occa
  pinva(i,i)=pinva(i,i)-1d0
enddo
! monomer blocks P alpha
pinvAAa(1:A%NOa,1:A%NOa)=pinva(1:A%NOa,1:A%NOa)
pinvABa(1:A%NOa,1:B%NOa)=pinva(1:A%NOa,A%NOa+1:occa)
pinvBAa(1:B%NOa,1:A%NOa)=pinva(A%NOa+1:occa,1:A%NOa)
pinvBBa(1:B%NOa,1:B%NOa)=pinva(A%NOa+1:occa,A%NOa+1:occa)

deallocate(ipiv)
deallocate(pinva,pa)

! invert beta
allocate(ipiv(occb))
allocate(pinvb(occb,occb))
allocate(pinvAAb(A%NOb,A%NOb),pinvBBb(B%NOb,B%NOb), &
         pinvABb(A%NOb,B%NOb),pinvBAb(B%NOb,A%NOb))
pinvb=0d0
do concurrent(i=1:occb)
   pinvb(i,i)=1d0
enddo
call dgesv(occb,occb,pb,occb,ipiv,pinvb,occb,info)
! -1 beta
do i=1,occb
  pinvb(i,i)=pinvb(i,i)-1d0
enddo
! monomer blocks P beta
pinvAAb(1:A%NOb,1:A%NOb)=pinvb(1:A%NOb,1:A%NOb)
pinvABb(1:A%NOb,1:B%NOb)=pinvb(1:A%NOb,A%NOb+1:occb)
pinvBAb(1:B%NOb,1:A%NOb)=pinvb(A%NOb+1:occb,1:A%NOb)
pinvBBb(1:B%NOb,1:B%NOb)=pinvb(A%NOb+1:occb,A%NOb+1:occb)

deallocate(ipiv)
deallocate(pinvb,pb)

! Tmatrix alpha
allocate(taa(nao,nao),tba(nao,nao),taba(nao,nao),tbaa(nao,nao))

call tran2mo2ao_gen(pinvAAa,A%NOa,A%NOa,NAO,NMO,A%UMO(:,:,1),A%UMO(:,:,1),taa)
call tran2mo2ao_gen(pinvBBa,B%NOa,B%NOa,NAO,NMO,B%UMO(:,:,1),B%UMO(:,:,1),tba)
call tran2mo2ao_gen(pinvABa,A%NOa,B%NOa,NAO,NMO,A%UMO(:,:,1),B%UMO(:,:,1),taba)
call tran2mo2ao_gen(pinvBAa,B%NOa,A%NOa,NAO,NMO,B%UMO(:,:,1),A%UMO(:,:,1),tbaa)

if(SAPT%IPrint>=50) then
   print*, 'Tmatrix alpha'
   print*, 'tAa  = ', norm2(taa)
   print*, 'tBa  = ', norm2(tba)
   print*, 'tABa = ', norm2(taba)
   print*, 'tBAa = ', norm2(tbaa)
endif

deallocate(pinvBBa,pinvBAa,pinvABa,pinvAAa)

! Tmatrix beta
allocate(tab(nao,nao),tbb(nao,nao),tabb(nao,nao),tbab(nao,nao))

call tran2mo2ao_gen(pinvAAb,A%NOb,A%NOb,NAO,NMO,A%UMO(:,:,2),A%UMO(:,:,2),tab)
call tran2mo2ao_gen(pinvBBb,B%NOb,B%NOb,NAO,NMO,B%UMO(:,:,2),B%UMO(:,:,2),tbb)
call tran2mo2ao_gen(pinvABb,A%NOb,B%NOb,NAO,NMO,A%UMO(:,:,2),B%UMO(:,:,2),tabb)
call tran2mo2ao_gen(pinvBAb,B%NOb,A%NOb,NAO,NMO,B%UMO(:,:,2),A%UMO(:,:,2),tbab)

if(SAPT%IPrint>=50) then
   print*, 'Tmatrix beta'
   print*, 'tAb  = ', norm2(tab)
   print*, 'tBb  = ', norm2(tbb)
   print*, 'tABb = ', norm2(tabb)
   print*, 'tBAb = ', norm2(tbab)
endif

deallocate(pinvBBb,pinvBAb,pinvABb,pinvAAb)

! Coulomb exchange alpha
allocate(Jtaa(nao,nao),Ktaa(nao,nao))
allocate(Jtba(nao,nao),Ktba(nao,nao))
allocate(Jtaba(nao,nao),Ktaba(nao,nao))

call make_J1(nao,taa, Jtaa, 'AOTWOSORT')
call make_J1(nao,tba, Jtba, 'AOTWOSORT')
call make_J1(nao,taba,Jtaba,'AOTWOSORT')

call make_K(nao,taa, Ktaa)
call make_K(nao,tba, Ktba)
call make_K(nao,taba,Ktaba)

if(SAPT%IPrint>=50) then
   write(6,'(/1x,a)')'Jmat alpha'
   print*, 'Jtaa   =', norm2(Jtaa)
   print*, 'Jtba   =', norm2(Jtba)
   print*, 'Jtaba  =', norm2(Jtaba)
   print*,'Kmat alpha'
   print*, 'Ktaa   =', norm2(Ktaa)
   print*, 'Ktba   =', norm2(Ktba)
   print*, 'Ktaba  =', norm2(Ktaba)
endif

! Coulomb exchange beta
allocate(Jtab(nao,nao),Ktab(nao,nao))
allocate(Jtbb(nao,nao),Ktbb(nao,nao))
allocate(Jtabb(nao,nao),Ktabb(nao,nao))

call make_J1(nao,tab, Jtab, 'AOTWOSORT')
call make_J1(nao,tbb, Jtbb, 'AOTWOSORT')
call make_J1(nao,tabb,Jtabb,'AOTWOSORT')

call make_K(nao,tab, Ktab)
call make_K(nao,tbb, Ktbb)
call make_K(nao,tabb,Ktabb)

write(6,'(/1x,a)') 'Jmat beta '
print*, 'Jtab   =', norm2(Jtab)
print*, 'Jtbb   =', norm2(Jtbb)
print*, 'Jtabb  =', norm2(Jtabb)
print*,'Kmat beta '
print*, 'Ktab   =', norm2(Ktab)
print*, 'Ktbb   =', norm2(Ktbb)
print*, 'Ktabb  =', norm2(Ktabb)

! calculations
pex=0d0
!PART(1): PA_alfa*KB_alfa + PA_beta*KB_beta
call dgemm('N','T',nao,nao,nao,1d0,PAa,nao,Kba,nao,0d0,work,nao)
pex(1) = -trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,PAb,nao,Kbb,nao,0d0,work,nao)
pex(1) = pex(1) - trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (1) = ', pex(1)

!PART(2): TA_alfa*hB_alfa + TA_beta*hB_beta
call dgemm('N','T',nao,nao,nao,1d0,taa,nao,hba,nao,0d0,work,nao)
pex(2) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tab,nao,hbb,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(2) = pex(2) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (2) = ', pex(2)

!PART(3): TB_alfa*hA_alfa + TB_beta*hA_beta
call dgemm('N','T',nao,nao,nao,1d0,tba,nao,haa,nao,0d0,work,nao)
pex(3) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tbb,nao,hab,nao,0d0,work,nao)
pex(3) = pex(3) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (3) = ', pex(3)

!PART(4): T^{AB}_alfa*h^A_alfa + T^{AB}_beta*h^A_beta
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,haa,nao,0d0,work,nao)
pex(4) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,hab,nao,0d0,work,nao)
pex(4) = pex(4) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (4) = ', pex(4)

!PART(5): T^{AB}_alfa*h^B_alfa + T^{AB}_beta*h^B_beta
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,hba,nao,0d0,work,nao)
pex(5) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,hbb,nao,0d0,work,nao)
pex(5) = pex(5) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (5) = ', pex(5)

!PART(6): T^{AB}_alfa(J[T^B_alfa]-K[T^B_alfa])
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtba,nao,0d0,work,nao)
pex(6) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktba,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
! beta
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtbb,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tabb,nao,ktbb,nao,0d0,work,nao)
pex(6) = pex(6) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (6) = ', pex(6)

!L terms (included in part-6)
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtbb,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(6) = pex(6) + pmix
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(6) = pex(6) + pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(6)+mix= ', pex(6)

!PART(7): T^{AB}_alfa*(J[T^A_alfa]-K[T^A_alfa]) + T^{AB}_beta*(J[]-K[])
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtaa,nao,0d0,work,nao)
pex(7) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktaa,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tabb,nao,jtab,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tabb,nao,ktab,nao,0d0,work,nao)
pex(7) = pex(7) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (7) = ', pex(7)

!L2 terms (included in part-7)
call dgemm('N','T',nao,nao,nao,1d0,taba,nao,jtab,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(7) = pex(7) + pmix
call dgemm('N','T',nao,nao,nao,1d0,tabb,nao,jtaa,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(7) = pex(7) + pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(7)+mix= ', pex(7)

!PART(8): T^A_alfa(J[T^B_alfa]-K[T^B_alfa]) + T^A_beta(J[T^B_beta]-K[T^B_beta])
call dgemm('N','T',nao,nao,nao, 1d0,taa,nao,jtba,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taa,nao,ktba,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tab,nao,jtbb,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tab,nao,ktbb,nao,0d0,work,nao)
pex(8) = pex(8) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (8) = ', pex(8)

! mix 1: T^A_alpha . J^B_beta
call dgemm('N','T',nao,nao,nao, 1d0,taa,nao,jtbb,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(8) = pex(8) + pmix
! mix 2: T^A_beta . J^B_alpha
call dgemm('N','T',nao,nao,nao, 1d0,tab,nao,jtba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(8) = pex(8) + pmix
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(8)+mix= ', pex(8)

!PART(9): TAB_a(J[TAB_a]-K[TAB_a])
call dgemm('N','T',nao,nao,nao, 1d0,taba,nao,jtaba,nao,0d0,work,nao)
pex(9) = trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,taba,nao,ktaba,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
! beta
call dgemm('N','T',nao,nao,nao, 1d0,tabb,nao,jtabb,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
call dgemm('N','T',nao,nao,nao,-1d0,tabb,nao,ktabb,nao,0d0,work,nao)
pex(9) = pex(9) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'part (9) = ', pex(9)

! mix 2*T^AB_beta.J[T^AB_alpha]
call dgemm('N','T',nao,nao,nao,2d0,tabb,nao,jtaba,nao,0d0,work,nao)
pmix = trace(work,nao)
pex(9) = pex(9) + trace(work,nao)
if(SAPT%IPrint>=10) write(LOUT,'(1x,a,f16.8)') 'pt(9)+mix= ', pex(9)

e1ex = sum(pex)
SAPT%e1exch = e1ex
write(6,*) 'E1ex = ', e1ex*1000

call print_en('E1exch',e1ex*1.0d3,.true.)

deallocate(Pab,Paa)
deallocate(Pba,Pbb)
deallocate(Jtab,Jtbb,Jtabb)
deallocate(Ktab,Ktbb,Ktabb)
deallocate(Jtaa,Jtba,Jtaba)
deallocate(Ktaa,Ktba,Ktaba)

deallocate(Sb,Sa)
deallocate(work)

end subroutine e1exch_os

end module sapt_exch_os
