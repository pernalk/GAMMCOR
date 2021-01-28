module sref_exch
use types
use exmisc

contains

subroutine exdisp_sref(elst,Vnn,termZ,sij,Sab,Sba,Vaab,Vaba,Vbab,Vabb,Vbaa,posA,posB,NBas,A,B)

implicit none

type(SystemBlock) :: A, B
integer,intent(in) :: NBas
double precision,intent(in) :: elst,Vnn,termZ 
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas),Sba(NBas,NBas)
double precision,intent(in) :: sij(A%NDimX,B%NDimX)
double precision,intent(in) :: Vaab(NBas,NBas),Vaba(NBas,NBas),&
                               Vbab(NBas,NBas),Vabb(NBas,NBas),&
                               Vbaa(NBas,NBas)

double precision :: val,termY,e2exd
double precision,allocatable :: tpqrs(:,:)

 allocate(tpqrs(2*A%NDimX,2*B%NDimX))

 e2exd = 0
 e2exd = e2exd + termZ
 call prep_tpqrs(sij,tpqrs,posA,posB,A%NDimX,B%NDimX,NBas,A,B) 
 call term_Y_SR(tpqrs,Sab,elst,Vnn,A%NDimX,B%NDimX,NBas,A,B,termY)
 e2exd = e2exd + termY

 call term_X_A1_SR(tpqrs,Sab,Vaab,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,val)
 e2exd = e2exd + val
!
 call term_X_A2_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBBA',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 call term_X_A2_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAAAB',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.,val)
 e2exd = e2exd + val

 call term_X_A2_YY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 call term_X_A2_YY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.,val)
 e2exd = e2exd + val

 call term_X_A2_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 call term_X_A2_YX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA','FFOOBAAA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.,val)
 e2exd = e2exd + val

 call term_X_A2_YX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOBBAB','FFOOABBB',Sab,Vabb,Vbab,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 call term_X_A2_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABA',Sba,Vbaa,Vaba,posB,posA,B%NDimX,A%NDimX,NBas,B,A,.true.,val)
 e2exd = e2exd + val

 call term_X_A3_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 call term_X_A3_XX_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.true.,val)
 e2exd = e2exd + val

 call term_X_A3_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.true.,val)
 e2exd = e2exd + val
 call term_X_A3_XY_SR(A%NDimX,B%NDimX,tpqrs,'FOFOAABB',Sab,Vabb,Vbaa,posA,posB,A%NDimX,B%NDimX,NBas,A,B,.false.,val)
 e2exd = e2exd + val
 
 write(6,'(1x,a,f16.8)') 'Single-reference code =',e2exd*1000d0

 deallocate(tpqrs)

end subroutine exdisp_sref

subroutine prep_tpqrs(sij,tpqrs,posA,posB,ADimX,BDimX,NBas,A,B)

implicit none 

type(SystemBlock) :: A, B

integer,intent(in) :: ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas) 
double precision,intent(in)  :: sij(ADimX,BDimX)
double precision,intent(out) :: tpqrs(2*ADimX,2*BDimX)

integer :: i,j,k,l,iunit
integer :: ip,iq,ir,is,ipq,irs
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,test
double precision,allocatable :: work(:),tmp(:,:)
double precision,parameter :: SmallE = 1.D-3
double precision,parameter :: BigE = 1.D8 

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 allocate(tmp(ADimX,BDimX),work(noVB))

 tmp=0
 do j=1,B%NDimX
    do i=1,A%NDimX
        
     !  if(A%Eig(i).gt.SmallE.and.B%Eig(j).gt.SmallE&
     !     .and.A%Eig(i).lt.BigE.and.B%Eig(j).lt.BigE) then

          tmp(i,j) = 16d0*sij(i,j)/(A%Eig(i)+B%Eig(j))

     !  endif
    enddo
 enddo

 tpqrs = 0
 do l=1,BDimX

    ir = B%IndN(1,l)
    is = B%IndN(2,l)
    irs = posB(ir,is)

    do k=1,ADimX

       ip = A%IndN(1,k)
       iq = A%IndN(2,k)
       ipq = posA(ip,iq)

       fact = (A%Occ(ip)-A%Occ(iq)) * &
              (B%Occ(ir)-B%Occ(is))

       ! XX
       val = 0
       do j=1,BDimX
          do i=1,ADimX

             val = val - fact * tmp(i,j) * &
                              A%EigX(ipq+(i-1)*A%NDimX) * &
                              B%EigX(irs+(j-1)*B%NDimX)

          enddo
       enddo

       tpqrs(ipq,irs) = tpqrs(ipq,irs) + val

       ! YY

       fact = (A%Occ(iq)-A%Occ(ip)) * &
              (B%Occ(is)-B%Occ(ir))

       val = 0
       do j=1,BDimX
          do i=1,ADimX

             val = val - fact * tmp(i,j) * &
                              A%EigY(ipq+(i-1)*A%NDimX) * &
                              B%EigY(irs+(j-1)*B%NDimX)

          enddo
       enddo
        
       tpqrs(A%NDimX+ipq,B%NDimX+irs) = tpqrs(A%NDimX+ipq,B%NDimX+irs) + val

       ! XA.YB

       fact = (A%Occ(ip)-A%Occ(iq)) * &
              (B%Occ(is)-B%Occ(ir))

       val = 0
       do j=1,BDimX
          do i=1,ADimX

             val = val - fact * tmp(i,j) * &
                              A%EigX(ipq+(i-1)*A%NDimX) * &
                              B%EigY(irs+(j-1)*B%NDimX)

          enddo
       enddo

       tpqrs(ipq,B%NDimX+irs) = tpqrs(ipq,B%NDimX+irs) + val

       ! YA.XB

       fact = (A%Occ(iq)-A%Occ(ip)) * &
              (B%Occ(ir)-B%Occ(is))

       val = 0
       do j=1,BDimX
          do i=1,ADimX

             val = val - fact * tmp(i,j) * &
                              A%EigY(ipq+(i-1)*A%NDimX) * &
                              B%EigX(irs+(j-1)*B%NDimX)

          enddo
       enddo

       tpqrs(A%NDimX+ipq,irs) = tpqrs(A%NDimX+ipq,irs) + val

    enddo
 enddo

 ! E2disp-test 
 open(newunit=iunit,file='TWOMOAB',status='OLD',&
      access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

 test = 0
 do ipq=1,ADimX

    ip = A%IndN(1,ipq)
    iq = A%IndN(2,ipq)
    read(iunit,rec=iq+(ip-A%num0-1)*dimOA) work(1:nOVB)

    do irs=1,BDimX

       ir = B%IndN(1,irs)
       is = B%IndN(2,irs)


       test = test + tpqrs(ipq,irs)*work(is+(ir-B%num0-1)*dimOB)  &
                   + tpqrs(ipq,B%NDimX+irs)*work(is+(ir-B%num0-1)*dimOB) &
                   + tpqrs(A%NDimX+ipq,irs)*work(is+(ir-B%num0-1)*dimOB) &
                   + tpqrs(ADimX+ipq,BDimX+irs)*work(is+(ir-B%num0-1)*dimOB)

    enddo
 enddo

 close(iunit)
 print*, 'E2disp/test',test*1000

 print*, 'tpqrs0',norm2(tpqrs)

deallocate(work,tmp)

end subroutine prep_tpqrs 

subroutine term_Y_SR(tpqrs,Sab,elst,vnn,ADimX,BDimX,NBas,A,B,termY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: ADimX,BDimX,NBas
double precision,intent(in) :: elst,vnn
double precision,intent(in) :: Sab(NBas,NBas),tpqrs(2*ADimX,2*BDimX)
double precision, intent(out) :: termY

integer :: i,j
integer :: ip,iq,ir,is,ipq,irs
double precision :: tmp(NBas,NBas)

tmp = transpose(Sab)

termY = 0
do j=1,BDimX

   ir = B%IndN(1,j) 
   is = B%IndN(2,j)
   irs = j

   do i=1,ADimX

      ip = A%IndN(1,i)
      iq = A%IndN(2,i)
      ipq = i
  
      termY = termY + tpqrs(ipq,irs)*Sab(ip,is)*tmp(ir,iq) &
                    + tpqrs(ipq,BDimX+irs)*Sab(ip,ir)*tmp(is,iq) &
                    + tpqrs(ADimX+ipq,irs)*Sab(iq,is)*tmp(ir,ip) &
                    + tpqrs(ADimX+ipq,BDimX+irs)*Sab(iq,ir)*tmp(is,ip)

 
   enddo
enddo
!print*, 'termY-000',termY/16d0
termY = 0.5d0*(elst-vnn)*termY

end subroutine term_Y_SR

subroutine term_X_A1_SR(tpqrs,Sab,Vaab,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,termA1)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vaab(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*ADimX,2*BDimX)
double precision,intent(out) :: termA1

integer :: i,j,k,l,kl
integer :: ip,iq,ir,is,ipq,irs
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: nelA,nelB,tmp(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 allocate(work(NBas*NBas),ints(NBas,NBas))

 !(FO|FO):(AB|BA)
 open(newunit=iunit,file='FOFOABBA',status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 termA1 = 0
 work = 0
 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       ir = k
       iq = l

       call ints_modify(NBas,dimOB,ints,NBas,work,Sab(iq,ir)/nelA,Vaab,Vbab(iq,ir)/nelB,Sab)

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                termA1 = termA1 + tpqrs(ipq,irs)*ints(ip,is)
 
             endif

          enddo
       enddo

       ! YY
       ir = k
       iq = l  
  
       call ints_modify(NBas,dimOB,ints,NBas,work,Vaab(iq,ir)/nelA,Sab,Sab(iq,ir)/nelB,Vbab)

       do ip=1,NBas
          do is=1,dimOB

             ipq = posA(ip,iq)
             irs = posB(ir,is)

             if(ipq/=0.and.irs/=0) then

                termA1 = termA1 + tpqrs(ADimX+ipq,BDimX+irs)*ints(ip,is)

             endif
   
          enddo
       enddo

    enddo
 enddo
 
 close(iunit)

! (FF|OO):(AB|AB)
 open(newunit=iunit,file='FFOOABAB',status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOA
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

       iq = k
       is = l

       call ints_modify(NBas,NBas,ints,NBas,work,Sab(iq,is)/nelA,Vaab,Vbab(iq,is)/nelB,Sab)

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)
             
             if(ipq/=0.and.irs/=0) then

               termA1 = termA1 + tpqrs(ipq,BDimX+irs)*ints(ip,ir)

             endif

          enddo
       enddo

       ! YX 
       call ints_modify(NBas,NBas,ints,NBas,work,Vaab(iq,is)/nelA,Sab,Sab(iq,is)/nelB,Vbab)

       do ir=1,NBas
          do ip=1,NBas

             ipq = posA(ip,iq)
             irs = posB(ir,is)
             
             if(ipq/=0.and.irs/=0) then

               termA1 = termA1 + tpqrs(ADimX+ipq,irs)*ints(ip,ir)

             endif

          enddo
       enddo

    enddo
 enddo

 close(iunit)

 termA1 = -0.5d0*termA1 
! print*, 'termA1',termA1

 !! test full
 !call tran4_gen(NBas,&
 !         NBas,B%CMO,&
 !         NBas,A%CMO,&
 !         NBas,A%CMO,&
 !         NBas,B%CMO,&
 !         'FFFFABBA','AOTWOSORT')

 !!(FF|FF):(AB|BA)
 !open(newunit=iunit,file='FFFFABBA',status='OLD', &
 !    access='DIRECT',recl=8*NBas*NBas)

 !termA1 = 0
 !work = 0
 !ints = 0
 !kl = 0
 !do l=1,NBas
 !   do k=1,NBas
 !      kl = kl + 1
 !      read(iunit,rec=kl) work(1:NBas*NBas)

 !      is = k
 !      ip = l

 !      call ints_modify(NBas,NBas,ints,NBas,work,Sab(ip,is)/nelA,Vaab,Vbab(ip,is)/nelB,Sab)

 !      do iq=1,NBas
 !         do ir=1,NBas

 !            ipq = posA(ip,iq)
 !            irs = posB(ir,is)
 ! 
 !            if(ipq/=0.and.irs/=0) then

 !               termA1 = termA1 + tpqrs(ipq,irs)*ints(iq,ir)

 !            endif
 !
 !         enddo
 !      enddo

 !   enddo
 !enddo

 !close(iunit,status='delete')

 !termA1 = -0.5d0*termA1 
 !print*, 'termA1(F)',termA1

 deallocate(ints,work)

end subroutine term_X_A1_SR

subroutine term_X_A3_XX_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbaa,posA,posB,ADimX,BDimX,NBas,A,B,isYY,termXX)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbaa(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: isYY
character(*),intent(in) :: IntKFile
double precision,intent(out) :: termXX

integer :: i,j,k,l,kl
integer :: iq,ib,ip,ir,ia,ic,iac,ipr
integer :: iunit,offA,offB
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,valTr,val3B,val4B
double precision :: nelA,nelB
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas)
double precision :: Emat(NBas,NBas),Amat(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)
! test
double precision :: valTmp

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

 Amat = 0
 do ib=1,dimOB
    do ir=1,NBas 
       do ip=1,NBas
          Amat(ip,ir) = Amat(ip,ir) + B%Occ(ib)*Sab(ip,ib)*Sab(ir,ib)
       enddo
    enddo
 enddo
 
 valTr = 0
 do iq=1,dimOA
    valTr = valTr + A%Occ(iq)*Amat(iq,iq)
 enddo
 
 Sbb = 0
 do ic=1,NBas
    do ia=1,NBas
       do iq=1,dimOA
          Sbb(ia,ic) = Sbb(ia,ic) + A%Occ(iq)*Sab(iq,ia)*Sab(iq,ic)
       enddo
    enddo
 enddo

 if(isYY) then
    offA = ADimX
    offB = BDimX
 else
    offA = 0
    offB = 0 
 endif

!(FO|FO):(AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

 ! one loop over integrals
 termXX = 0
 Emat = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOA)
 
       call ints_modify(NBas,dimOA,ints,NBas,work,&
                        Vabb(k,l)/nelA,Sxx,Sxx(l,k)/nelB,Vbaa)

       ic = l
       ia = k

       iac = posB(ia,ic)

       if(iac/=0) then

          do i=1,ADimX
             ip = A%IndN(1,i)
             ir = A%IndN(2,i)
             ipr = posA(ip,ir)
             
             ! 1A-1B
             termXX = termXX - 2d0*tpqrs(offA+ipr,offB+iac)*valTr*ints(ip,ir)
             valTmp = valTmp - 2d0*tpqrs(offA+ipr,offB+iac)*valTr*ints(ip,ir)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,iq)
             enddo

             ! 2A-1B
             termXX = termXX - 2d0*val*tpqrs(offA+ipr,offB+iac)*Amat(ip,ir)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(ip,iq)*Amat(iq,ir)
             enddo

             ! 3A-1B
             termXX = termXX + val*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,ir)*Amat(iq,ip)
             enddo

             ! 4A-1B
             termXX = termXX + val*tpqrs(offA+ipr,offB+iac)

          enddo
       endif

       ! 2B terms
       if(k==l) then

          ib = k
          Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + B%Occ(ib)*ints(1:NBas,1:dimOA)

       endif

       ib = l
       ia = k

       val3B = 0
       do iq=1,dimOA
          val3B = val3B + A%Occ(iq)*ints(iq,iq)
       enddo

       do ic=1,dimOB

          iac = posB(ia,ic)
 
          if(iac/=0) then

             do i=1,ADimX

                ip = A%IndN(1,i)
                ir = A%IndN(2,i)
                ipr = posA(ip,ir)

                ! 1A-3B
                termXX = termXX + ints(ip,ir)*Sbb(ic,ib)*tpqrs(offA+ipr,offB+iac)

                ! 2A-3B
                termXX = termXX + val3B*B%Occ(ib)*Sab(ip,ic)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ic)
                enddo

                ! 3A-3B
                termXX = termXX - val*B%Occ(ib)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)
 
                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ib)
                enddo

                ! 4A-3B
                termXX = termXX - val*B%Occ(ib)*Sab(ip,ic)*tpqrs(offA+ipr,offB+iac)

             enddo

          endif
       enddo

       if(k<=dimOB) then

          ib = k
          ic = l

          val4B = 0
          do iq=1,dimOA
             val4B = val4B + A%Occ(iq)*ints(iq,iq)
          enddo

          do ia=1,NBas

             iac = posB(ia,ic)
 
             if(iac/=0) then

                do i=1,ADimX

                   ip = A%IndN(1,i)
                   ir = A%IndN(2,i)
                   ipr = posA(ip,ir)

                   ! 1A-4B
                   termXX = termXX + B%Occ(ib)*Sbb(ia,ib)*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                   ! 2A-4B
                   termXX = termXX + val4B*B%Occ(ib)*Sab(ip,ib)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ib)
                   enddo 

                   ! 3A-4B
                   termXX = termXX - val*B%Occ(ib)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ia)
                   enddo 

                   ! 4A-4B
                   termXX = termXX - val*B%Occ(ib)*Sab(ip,ib)*tpqrs(offA+ipr,offB+iac)

                enddo
             endif

          enddo
       endif

    enddo
 enddo

 close(iunit)

 ! 2B terms
 do j=1,BDimX

    ia = B%IndN(1,j)
    ic = B%IndN(2,j)
    iac = posB(ia,ic) 

    do i=1,ADimX

       ip = A%IndN(1,i)
       ir = A%IndN(2,i)
       ipr = posA(ip,ir)

       ! 1A-2B 
       termXX = termXX - 2d0*Emat(ip,ir)*Sbb(ia,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,iq) 
       enddo

       ! 2A-2B 
       termXX = termXX - 2d0*val*Sab(ip,ic)*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(ip,iq)*Sab(iq,ic)
       enddo

       ! 3A-2B 
       termXX = termXX + val*Sab(ir,ia)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,ir)*Sab(iq,ia)
       enddo

       ! 4A-2B 
       termXX = termXX + val*Sab(ip,ic)*tpqrs(offA+ipr,offB+iac)
 
    enddo
 enddo

! print*, 'A3-XX',termXX
! print*, 'valTmp(XX):',valTmp

 deallocate(ints,work) 

end subroutine term_X_A3_XX_SR

subroutine term_X_A3_XY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbaa,posA,posB,ADimX,BDimX,NBas,A,B,isXY,termXY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbaa(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: isXY
character(*),intent(in) :: IntKFile
double precision,intent(out) :: termXY

integer :: i,j,k,l,kl
integer :: iq,ib,ip,ir,ia,ic,iac,ipr
integer :: iunit,offA,offB
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,valTr,val2A,val3B,val4B
double precision :: nelA,nelB
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas)
double precision :: Emat(NBas,NBas),Amat(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)
! test
double precision :: valTmp

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

 Amat = 0
 do ib=1,dimOB
    do ir=1,NBas 
       do ip=1,NBas
          Amat(ip,ir) = Amat(ip,ir) + B%Occ(ib)*Sab(ip,ib)*Sab(ir,ib)
       enddo
    enddo
 enddo
 
 valTr = 0
 do iq=1,dimOA
    valTr = valTr + A%Occ(iq)*Amat(iq,iq)
 enddo
 
 Sbb = 0
 do ic=1,NBas
    do ia=1,NBas
       do iq=1,dimOA
          Sbb(ia,ic) = Sbb(ia,ic) + A%Occ(iq)*Sab(iq,ia)*Sab(iq,ic)
       enddo
    enddo
 enddo

 if(isXY) then
    offA = 0
    offB = BDimX
 else
    offA = ADimX
    offB = 0 
 endif

!(FO|FO):(AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

 ! one loop over integrals
 termXY = 0
 Emat = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOA)

       call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/nelA,Sxx,Sxx(k,l)/nelB,Vbaa)

       ic = l
       ia = k

       iac = posB(ia,ic)

       if(iac/=0) then

          val2A = 0
          do iq=1,dimOA
             val2A = val2A + 2d0*A%Occ(iq)*ints(iq,iq)
          enddo
 
          do i=1,ADimX

             ip = A%IndN(1,i)
             ir = A%IndN(2,i)
             ipr = posA(ip,ir)

             ! 1A-1B
             termXY = termXY - 2d0*valTr*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)
             valTmp = valTmp - 2d0*valTr*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

             ! 2A-1B
             termXY = termXY - val2A*Amat(ip,ir)*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(ip,iq)*Amat(iq,ir)
             enddo

             ! 3A-1B
             termXY = termXY + val*tpqrs(offA+ipr,offB+iac)

             val = 0
             do iq=1,dimOA
                val = val + A%Occ(iq)*ints(iq,ir)*Amat(ip,iq)
             enddo

             ! 4A-1B
             termXY = termXY + val*tpqrs(offA+ipr,offB+iac)

          enddo
       endif

!       ! 2B terms 
       if(k==l.and.k<=dimOB) then

          ib = k
          Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + B%Occ(ib)*ints(1:NBas,1:dimOA)

       endif

       ia = k
       ib = l

       val4B = 0
       do iq=1,dimOA
          val4B = val4B + A%Occ(iq)*ints(iq,iq)
       enddo
       val4B = val4B*B%Occ(ib)

       do ic=1,dimOB

          iac = posB(ia,ic)

          if(iac/=0) then

             fact = B%Occ(ib)*Sbb(ic,ib)

             do i=1,ADimX

                ip = A%IndN(1,i)
                ir = A%IndN(2,i)
                ipr = posA(ip,ir)

                ! 1A-4B
                termXY = termXY + fact*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                ! 2A-4B
                termXY = termXY + val4B*Sab(ip,ib)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ib)
                enddo

                ! 3A-4B
                termXY = termXY - val*B%Occ(ib)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

                val = 0
                do iq=1,dimOA
                   val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ic)
                enddo

                ! 4A-4B
                termXY = termXY - val*B%Occ(ib)*Sab(ip,ib)*tpqrs(offA+ipr,offB+iac)

             enddo
          endif

       enddo

       if(k<=dimOB) then

          ib = k
          ic = l

          val3B = 0
          do iq=1,dimOA
             val3B = val3B + A%Occ(iq)*B%Occ(ib)*ints(iq,iq)
          enddo

          do ia=1,NBas

             iac = posB(ia,ic)

             if(iac/=0) then

                fact = B%Occ(ib)*Sbb(ia,ib)

                do i=1,ADimX
 
                   ip = A%IndN(1,i)
                   ir = A%IndN(2,i)
                   ipr = posA(ip,ir)

                   ! 1A-3B
                   termXY = termXY + fact*ints(ip,ir)*tpqrs(offA+ipr,offB+iac)

                   ! 2A-3B 
                   termXY = termXY + val3B*Sab(ip,ia)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(ip,iq)*Sab(iq,ia)
                   enddo

                   ! 3A-3B 
                   termXY = termXY - val*B%Occ(ib)*Sab(ir,ib)*tpqrs(offA+ipr,offB+iac)

                   val = 0
                   do iq=1,dimOA
                      val = val + 0.5d0*A%Occ(iq)*ints(iq,ir)*Sab(iq,ib)
                   enddo

                   ! 4A-3B 
                   termXY = termXY - val*B%Occ(ib)*Sab(ip,ia)*tpqrs(offA+ipr,offB+iac)

                enddo

             endif
          enddo

       endif

    enddo
 enddo

 close(iunit)

 ! 2B terms
 do j=1,BDimX
   
    ia = B%IndN(1,j)
    ic = B%IndN(2,j)
    iac = posB(ia,ic)

    fact = 2d0*Sbb(ia,ic)

    do i=1,ADimX

       ip = A%IndN(1,i)
       ir = A%IndN(2,i)
       ipr = posA(ip,ir)

       ! 1A-2B
       termXY = termXY - fact*Emat(ip,ir)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + 2d0*A%Occ(iq)*Emat(iq,iq)
       enddo 

       ! 2A-2B
       termXY = termXY - val*Sab(ip,ia)*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(ip,iq)*Sab(iq,ia)
       enddo

       ! 3A-2B
       termXY = termXY + val*Sab(ir,ic)*tpqrs(offA+ipr,offB+iac)

       val = 0
       do iq=1,dimOA
          val = val + A%Occ(iq)*Emat(iq,ir)*Sab(iq,ic)
       enddo

       ! 4A-2B
       termXY = termXY + val*Sab(ip,ia)*tpqrs(offA+ipr,offB+iac)

    enddo
 enddo

! print*, 'term A3-XY',termXY
! print*, 'valTmp(XY):',valTmp

 deallocate(ints,work)

end subroutine term_X_A3_XY_SR

subroutine term_X_A2_XX_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans,termXX)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile
double precision,intent(out) :: termXX

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|BA) or (AA|AB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

 ! one loop over integrals
 termXX = 0
 ints = 0
 kl = 0
 do l=1,dimOA
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)
 
       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(l,k)/nelA,Vabb,Vbab(l,k)/nelB,Sxx)
 
       ! 1 and 4
       if(k<=dimOB) then

          iq = k
          iu = l

          do it=1,NBas

             itu = posA(it,iu)
 
             fact = 2d0*Sab(it,iq)*B%Occ(iq)           

             if(itu/=0) then

                do i=1,BDimX 

                   ip = B%IndN(1,i)
                   ir = B%IndN(2,i)
                   ipr = posB(ip,ir)

                   if(trans) then
                      termXX = termXX + fact*tpqrs(ipr,itu)*ints(ip,ir)
                      termXX = termXX - tpqrs(ipr,itu)*ints(ip,iq)*B%Occ(iq)*Sab(it,ir)
                   else
                      termXX = termXX + fact*tpqrs(itu,ipr)*ints(ip,ir)
                      termXX = termXX - tpqrs(itu,ipr)*ints(ip,iq)*B%Occ(iq)*Sab(it,ir)
                   endif

                enddo
             endif 

          enddo
       endif

       ip = k
       iu = l

       do ir=1,dimOB
  
          ipr = posB(ip,ir)
         
          if(ipr/=0) then

             do it=1,NBas

                itu = posA(it,iu) 
                
                if(itu/=0) then
                 
                   ! 2  
                   val = 0
                   do iq=1,dimOB
                      val = val + ints(iq,ir)*B%Occ(iq)*Sab(it,iq)
                   enddo

                   if(trans) then
                      termXX = termXX - val*tpqrs(ipr,itu)
                   else
                      termXX = termXX - val*tpqrs(itu,ipr)
                   endif

                   ! 3
                   val = 0
                   do iq=1,dimOB
                      val = val + 2d0*ints(iq,iq)*B%Occ(iq)
                   enddo

                   if(trans) then
                      termXX = termXX + val*tpqrs(ipr,itu)*Sab(it,ir)
                   else
                      termXX = termXX + val*tpqrs(itu,ipr)*Sab(it,ir)
                   endif
                 
                endif   

             enddo
          endif
       enddo

    enddo
 enddo

 termXX = -0.5d0*termXX

 close(iunit)

! print*, 'termXX',termXX

 deallocate(ints,work)

end subroutine term_X_A2_XX_SR

subroutine term_X_A2_YY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans,termYY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile
double precision,intent(out) :: termYY

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 termYY = 0
 ints = 0
 work = 0
 kl = 0
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)
 
       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)
      
       it = k
       iq = l

       do iu=1,dimOA

          itu = posA(it,iu)
          fact = 2d0*B%Occ(iq)*Sab(iu,iq)

          if(itu/=0) then

             do i=1,B%NDimX

                ip = B%IndN(1,i)
                ir = B%IndN(2,i)
                ipr = posB(ip,ir)

                ! 1 and 4
                if(trans) then
                   termYY = termYY + fact*tpqrs(BDimX+ipr,ADimX+itu)*ints(ip,ir)
                   termYY = termYY - tpqrs(BDimX+ipr,ADimX+itu)*ints(ir,iq)*Sab(iu,ip)*B%Occ(iq)
                else
                   termYY = termYY + fact*tpqrs(ADimX+itu,BDimX+ipr)*ints(ip,ir)
                   termYY = termYY - tpqrs(ADimX+itu,BDimX+ipr)*ints(ir,iq)*Sab(iu,ip)*B%Occ(iq)
                endif

             enddo 

          endif 
       enddo

       it = k
       ir = l

       do ip=1,NBas

          ipr = posB(ip,ir)

          if(ipr/=0) then

             do iu=1,dimOA

                itu = posA(it,iu)

                if(itu/=0) then
 
                   val = 0
                   do iq=1,dimOB
                      val = val + ints(ip,iq)*Sab(iu,iq)*B%Occ(iq)
                   enddo

                   ! 2
                   if(trans) then
                      termYY = termYY - val*tpqrs(BDimX+ipr,ADimX+itu)
                   else
                      termYY = termYY - val*tpqrs(ADimX+itu,BDimX+ipr)
                   endif

                   val = 0
                   do iq=1,dimOB
                      val = val + 2d0*ints(iq,iq)*B%Occ(iq)
                   enddo

                   ! 3
                   if(trans) then
                      termYY = termYY + val*tpqrs(BDimX+ipr,ADimX+itu)*Sab(iu,ip)
                   else
                      termYY = termYY + val*tpqrs(ADimX+itu,BDimX+ipr)*Sab(iu,ip)
                   endif

                endif

             enddo
          endif
       enddo

    enddo
 enddo

 termYY = -0.5d0*termYY
 
 close(iunit)
 
! print*, 'termYY',termYY

 deallocate(ints,work)

end subroutine term_X_A2_YY_SR 

subroutine term_X_A2_XY_SR(dim1,dim2,tpqrs,IntKFile,Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans,termXY)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile
double precision,intent(out) :: termXY

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 termXY = 0
 ints = 0
 kl = 0
 
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       if(k<=dimOA) then

         call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)

         iu = k 
         iq = l

         do it=1,NBas

            itu = posA(it,iu)

            fact = 2d0*B%Occ(iq)*Sab(it,iq)

            if(itu/=0) then

               do i=1,BDimX

                  ip = B%IndN(1,i)
                  ir = B%IndN(2,i)
                  ipr = posB(ip,ir)

                  ! 1 and 4 
                  if(trans) then
                     termXY = termXY + fact*tpqrs(BDimX+ipr,itu)*ints(ip,ir)
                     termXY = termXY - tpqrs(BDimX+ipr,itu)*ints(ir,iq)*Sab(it,ip)*B%Occ(iq)
                  else
                     termXY = termXY + fact*tpqrs(itu,BDimX+ipr)*ints(ip,ir)
                     termXY = termXY - tpqrs(itu,BDimX+ipr)*ints(ir,iq)*Sab(it,ip)*B%Occ(iq)
                  endif

               enddo 

            endif
         enddo

         iu = k
         ir = l
 
         do ip=1,NBas

            ipr = posB(ip,ir)
 
            if(ipr/=0) then

               do it=1,NBas

                  itu = posA(it,iu)

                  if(itu/=0) then

                     val = 0
                     do iq=1,dimOB
                        val = val + B%Occ(iq)*ints(ip,iq)*Sab(it,iq)
                     enddo

                     ! 2
                     if(trans) then
                        termXY = termXY - val*tpqrs(BDimX+ipr,itu)
                     else
                        termXY = termXY - val*tpqrs(itu,BDimX+ipr)
                     endif

                     val = 0
                     do iq=1,dimOB
                        val = val + 2d0*B%Occ(iq)*ints(iq,iq)
                     enddo

                     ! 3
                     if(trans) then
                        termXY = termXY + val*tpqrs(BDimX+ipr,itu)*Sab(it,ip)
                     else
                        termXY = termXY + val*tpqrs(itu,BDimX+ipr)*Sab(it,ip)
                     endif
                    
                  endif
               enddo
            endif
         enddo

       endif

    enddo
 enddo

 termXY = -0.5d0*termXY

 close(iunit)

! print*, 'termXY',termXY

 deallocate(ints,work)

end subroutine term_X_A2_XY_SR

subroutine term_X_A2_YX_SR(dim1,dim2,tpqrs,IntKFile,IntJFile, &
                           Sab,Vabb,Vbab,posA,posB,ADimX,BDimX,NBas,A,B,trans,termYX)

implicit none

type(SystemBlock) :: A, B

integer,intent(in) :: dim1,dim2,ADimX,BDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas)
double precision,intent(in) :: Sab(NBas,NBas), &
                               Vabb(NBas,NBas),Vbab(NBas,NBas), &
                               tpqrs(2*dim1,2*dim2)
logical,intent(in) :: trans
character(*),intent(in) :: IntKFile,IntJFile
double precision,intent(out) :: termYX

integer :: i,j,k,l,kl
integer :: iq,ip,ir,it,iu,itu,ipr
integer :: iunit
integer :: dimOA,dimVA,dimOB,dimVB,nOVA,nOVB
double precision :: fact,val,nelA,nelB
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)
double precision,allocatable :: AuxMat(:,:) 

 ! dimensions
 dimOA = A%INAct+A%NAct
 dimVA = A%num1+A%num2
 dimOB = B%INAct+B%NAct
 dimVB = B%num1+B%num2
 nOVA = dimOA*dimVA
 nOVB = dimOB*dimVB

 nelA = 2d0*A%XELE
 nelB = 2d0*B%XELE

 Sxx = 0
 do i=1,NBas
    Sxx(i,i) = 1d0
 enddo

 allocate(work(NBas*NBas),ints(NBas,NBas))
 allocate(AuxMat(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
 open(newunit=iunit,file=trim(IntKFile),status='OLD', &
      access='DIRECT',recl=8*NBas*dimOB)

 termYX = 0
 work = 0
 ints = 0
 kl = 0
 
 do l=1,dimOB
    do k=1,NBas
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*dimOB)

       call ints_modify(NBas,dimOB,ints,NBas,work,&
                        Sab(k,l)/nelA,Vabb,Vbab(k,l)/nelB,Sxx)

       it = k
       iq = l

       do iu=1,dimOA

          itu = posA(it,iu)

          fact = 2d0*B%Occ(iq)*Sab(iu,iq)

          if(itu/=0) then

             do i=1,BDimX

                ip = B%IndN(1,i)
                ir = B%IndN(2,i)
                ipr = posB(ip,ir)

                ! 1 and 4
                if(trans) then
                   termYX = termYX + fact*tpqrs(ipr,ADimX+itu)*ints(ip,ir)
                   termYX = termYX - tpqrs(ipr,ADimX+itu)*ints(ip,iq)*Sab(iu,ir)*B%Occ(iq)
                else
                   termYX = termYX + fact*tpqrs(ADimX+itu,ipr)*ints(ip,ir)
                   termYX = termYX - tpqrs(ADimX+itu,ipr)*ints(ip,iq)*Sab(iu,ir)*B%Occ(iq)
                endif

             enddo

          endif
       enddo

    enddo
 enddo

 close(iunit)

!(FF|OO):(AB|BB) or (BA|AA)
 open(newunit=iunit,file=trim(IntJFile),status='OLD', &
      access='DIRECT',recl=8*NBas*NBas)

 AuxMat = 0
 work = 0
 ints = 0
 kl = 0
 do l=1,dimOB
    do k=1,dimOB
       kl = kl + 1
       read(iunit,rec=kl) work(1:NBas*NBas)

      ! call ints_modify(NBas,NBas,ints,NBas,work,&
      !                  Vabb(k,l)/nelA,Sab,Sxx(k,l)/nelB,Vbab)
      ! iq = k
      ! ir = l

      ! any difference?
      call ints_modify(NBas,NBas,ints,NBas,work,&
                       Sxx(k,l)/nelB,Vbab,Vabb(k,l)/nelA,Sab)
      iq = l
      ir = k

       do ip=1,NBas

          ipr = posB(ip,ir)

          if(ipr/=0) then

             do i=1,ADimX

                it = A%IndN(1,i)
                iu = A%IndN(2,i)
                itu = posA(it,iu)

                ! sth wrong here...?
                if(trans) then
                   ! 2
                   termYX = termYX - B%Occ(iq)*tpqrs(ipr,ADimX+itu)*ints(it,ip)*Sab(iu,iq)
                   ! 3
                   if(iq==ir) AuxMat = ints
                   termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ipr,ADimX+itu)*AuxMat(it,ip)*Sab(iu,ir)
                      !termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ipr,ADimX+itu)*ints(it,ip)*Sab(iu,ir)
                else 
                   ! 2
                   termYX = termYX - B%Occ(iq)*tpqrs(ADimX+itu,ipr)*ints(it,ip)*Sab(iu,iq)
                   ! 3
                   ! new?
                   if(iq==ir) AuxMat = ints 
                   termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ADimX+itu,ipr)*AuxMat(it,ip)*Sab(iu,ir)
                   ! old:
                   !if(iq==ir) then 
                   !   termYX = termYX + 2d0*B%Occ(iq)*tpqrs(ADimX+itu,ipr)*ints(it,ip)*Sab(iu,ir)
                   !endif
                endif

             enddo
          endif

       enddo

    enddo
 enddo

 close(iunit)

 termYX = -0.5d0*termYX

!print*, 'termYX',termYX

 deallocate(AuxMat)
 deallocate(ints,work)

end subroutine term_X_A2_YX_SR

end module sref_exch
