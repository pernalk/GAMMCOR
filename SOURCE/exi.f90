module exi

implicit none

private ints_modify

contains

subroutine ints_modify(n1,n2,ints,n,raw_ints,alphaX,X,alphaY,Y)
implicit none

integer,intent(in) :: n1,n2,n
double precision,intent(out) :: ints(n1,n2)
double precision,intent(in) :: alphaX,alphaY,raw_ints(n,n),X(n,n),Y(n,n)

ints(1:n1,1:n2) = raw_ints(1:n1,1:n2) + alphaX*X(1:n1,1:n2) + alphaY*Y(1:n1,1:n2)
!ints(1:n1,1:n2) = raw_ints(1:n1,1:n2)

end subroutine ints_modify

subroutine make_tind(tvec,Xmat,Ymat,Sab,AOcc,BOcc,AIndN,posA,NDimX,NBas)
implicit none

integer,intent(in) :: NBas,NDimX
integer,intent(in) :: AIndN(2,NDimX),posA(NBas,NBas)
double precision,intent(in) :: Xmat(NDimX,NDimX),Ymat(NDimX,NDimX)
double precision,intent(in) :: Sab(NBas,NBas),&
                               AOcc(NBas),BOcc(NBas)
double precision,intent(inout) :: tvec(NDimX)

integer :: i,ip,iq,ir,ipq
integer :: ii,j,k
double precision,allocatable :: tmp1(:),tmp2(:)
double precision :: fact,val

tvec = 0

!Sba=transpose(Sab)

allocate(tmp1(NDimX))
tmp1 = 0
do i=1,NDimX
   ip = AIndN(1,i)
   iq = AIndN(2,i)
   ipq = posA(ip,iq)

   fact = (AOcc(ip)-AOcc(iq))

   val = 0
   do ir=1,NBas
      !tmp1(ipq) = tmp1(ipq) + fact*BOcc(ir)*Sab(iq,ir)*Sba(ir,ip)
      val = val + BOcc(ir)*Sab(ip,ir)*Sab(iq,ir)
   enddo
   tmp1(ipq) = tmp1(ipq) + fact*val

enddo
!print*, 'tmp1',norm2(tmp1) 

call dgemv('T',NDimX,NDimX,1d0,Ymat-Xmat,NDimX,tmp1,1,0d0,tvec,1)

deallocate(tmp1)

end subroutine make_tind

subroutine make_tind2(tvec,Xmat,Ymat,Sab,Pmat,Occ,IndN,pos,dimOB,NDimX,NBas)
implicit none

integer,intent(in) :: dimOB,NBas,NDimX
integer,intent(in) :: IndN(2,NDimX),pos(NBas,NBas)
double precision,intent(in) :: Xmat(NDimX,NDimX),Ymat(NDimX,NDimX)
double precision,intent(in) :: Sab(NBas,NBas),Pmat(NBas,NBas),&
                               Occ(NBas)
double precision,intent(inout) :: tvec(NDimX)

integer :: i,ip,iq,ir,ipq
integer :: ii,j,k
double precision :: Sba(NBas,NBas)
double precision,allocatable :: tmp1(:,:),tmp2(:)
double precision :: fact,val

 tvec = 0
 Sba = transpose(Sab)

 allocate(tmp1(NBas,NBas),tmp2(NDimX))

 ! P(r,s).Sba(s,p)
 call dgemm('N','N',NBas,NBas,NBas,1d0,Pmat,NBas,Sba,NBas,0d0,tmp1,NBas)

 tmp2 = 0
 do i=1,NDimX
    ip = IndN(1,i)
    iq = IndN(2,i)
    ipq = pos(ip,iq)

    fact=(Occ(ip)-Occ(iq))
 
    val = 0
    !do ir=1,NBas
    do ir=1,dimOB
       !tmp2(ipq) = tmp2(ipq) + fact*Sab(iq,ir)*tmp1(ir,ip)
       val  = val + Sab(iq,ir)*tmp1(ir,ip)
    enddo
    tmp2(ipq) = tmp2(ipq) + fact*val

 enddo

 !! test
 !do i=1,NDimX
 !   ip = IndN(1,i)
 !   iq = IndN(2,i)
 !   ipq = pos(ip,iq)

 !   do ii=1,NDimX
 !     tvec(i) = tvec(i) + (Ymat(i,ii)-Xmat(i,ii))*tmp2(ii)
 !   enddo

 !enddo

 !print*, 'tmp2',norm2(tmp2)
 call dgemv('T',NDimX,NDimX,1d0,Ymat-Xmat,NDimX,tmp2,1,0d0,tvec,1)

 deallocate(tmp2,tmp1)

end subroutine make_tind2

subroutine exind_A3_XY_full(dim1,dim2,intA,intB,RDM2Aval,RDM2Bval,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX)
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intA(dim1),intB(dim2)
double precision,intent(in) :: AXELE,BXELE
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)
double precision,intent(in) :: RDM2Aval(dimOA,dimOA,dimOA,dimOA),RDM2Bval(dimOB,dimOB,dimOb,dimOB)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iae,ipt
integer :: it,is,id,ie
double precision :: factA,factB,fact,valX,valB
double precision :: Sxx(NBas,NBas),Sba(NBas,NBas),&
                    ABmat(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),&
                                ints(:,:),intsr(:,:),tints(:,:)
double precision,allocatable :: intAB(:,:,:,:),intV(:,:,:,:), &
                                calkA(:,:,:,:),calkB(:,:,:,:),&
                                tmpA(:,:,:,:),tmpB(:,:,:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

Sba = transpose(Sab)

allocate(intAB(dimOA,dimOA,dimOB,dimOB),intV(dimOA,dimOA,dimOB,dimOB))

allocate(calkA(dimOA,dimOA,dimOB,NBas))
!,calkB(dimOB,dimOB,dimOA,NBas))

allocate(tmp(dimOA**2,dimOB**2))
allocate(work(NBas*NBas),ints(NBas,NBas),intsr(dimOA,dimOA))

! create (qs,bd) intermediate
intAB = 0
do id=1,dimOB
   do ib=1,dimOB
      call dger(dimOA,dimOA,1d0,Sab(:,id),1,Sab(:,ib),1,intAB(:,:,ib,id),dimOA)
   enddo
enddo

! intAB: RDM2A(tr,qs).I(qs,bd).RDM2B(ec,bd)^T
call dgemm('N','N',dimOA**2,dimOB**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intAB,dimOA**2,0d0,tmp,dimOA**2)
call dgemm('N','N',dimOA**2,dimOB**2,dimOB**2,1d0,tmp,dimOA**2,RDM2Bval,dimOB**2,0d0,intAB,dimOA**2)

deallocate(tmp)
allocate(tmp(dimOA,dimOA))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
calkA = 0
intV = 0
ints = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ic = l
      ia = k

      do ie=1,dimOB

         iae = posB(ia,ie)
   
         if(iae/=0) then

            !1-1'(B)
            !print*, '1-1p B'
            do ir=1,dimOA
               do ip=1,dimOA
                  intB(iae) = intB(iae) + intAB(ip,ir,ie,ic)*ints(ip,ir)
               enddo
            enddo

         endif

      enddo

      ! calkA
      intsr(1:dimOA,1:dimOA) = ints(1:dimOA,1:dimOA)
      call dgemv('T',dimOA**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intsr,1,0d0,calkA(:,:,ic,ia),1)

      if(k<=dimOB) then

         ! intV
         !if(ia<=dimOB) then
           call dger(dimOA**2,dimOB**2,1d0,calkA(:,:,ic,ia),1,RDM2Bval(:,:,ia,ic),1,intV,dimOA**2)
         !endif 

         ! A
         do i=1,ADimX

            ip = AIndN(1,i)
            it = AIndN(2,i)
            ipt = posA(ip,it)

            valX = 0
            if(ip<=dimOA) then

               !X: 3-1' /  Y: 1-1'
               !print*, '3-1p A'
               do ir=1,dimOA
                 valX = valX - intAB(ir,ip,ia,ic)*ints(ir,it)
               enddo

            endif

            !X: 1-1' / Y: 3-1'
            !print*, '1-1p A'
            do ir=1,dimOA
               valX = valX + intAB(it,ir,ia,ic)*ints(ip,ir)
            enddo

            intA(ipt) = intA(ipt) + valX

         enddo
 
         ! B
         ie = l
         ic = k

         do ia=1,dimOB

            iae = posB(ia,ie)

            if(iae/=0) then

               !1-3'(B)
               !print*, '1-3p B'
               do ir=1,dimOA
                  do ip=1,dimOA
                     intB(iae) = intB(iae) - intAB(ip,ir,ic,ia)*ints(ip,ir)
                  enddo
               enddo

            endif

         enddo

      endif

   enddo
enddo

close(iunit)

!print*, 'intV',norm2(intV)

! A
do i=1,ADimX

   ip = AIndN(1,i)
   it = AIndN(2,i)
   ipt = posA(ip,it)

   valX = 0
   if(ip<=dimOA) then

     do id=1,dimOB
        do ib=1,dimOB

            ! 4-1'(X)
            !print*, '4-1p A'
            do iq=1,dimOA 
               valX = valX - intV(iq,ip,ib,id)*Sab(it,ib)*Sab(iq,id)
            enddo

        enddo
     enddo

   endif

   do id=1,dimOB
      do ib=1,dimOB

         ! 2-1'(X)
         !print*, '2-1p A'
         do is=1,dimOA
            valX = valX + intV(it,is,ib,id)*Sab(is,ib)*Sab(ip,id)
         enddo

      enddo 
   enddo 

   intA(ipt) = intA(ipt) + valX

enddo

! B
do j=1,BDimX

   ia = BIndN(1,j)
   ie = BIndN(2,j)
   iae = posB(ia,ie)

   valB = 0
   if(ia<=dimOB) then

      ! 1-4' 
      !print*, '1-4p B' 
      do ib=1,dimOB
         do is=1,dimOA
            do iq=1,dimOA
               valB = valB - intV(iq,is,ib,ia)*Sab(iq,ie)*Sab(is,ib)
            enddo
         enddo
      enddo

   endif

   ! 1-2'(B)
   !print*, '1-2p B' 
   do id=1,dimOB
       do is=1,dimOA
          do iq=1,dimOA
            valB = valB + intV(iq,is,ie,id)*Sab(iq,id)*Sab(is,ia)  
          enddo
       enddo
   enddo

   intB(iae) = intB(iae) + valB

enddo

intA = -2d0*intA
intB = -2d0*intB

print*, ''
print*, 'exind_A3_XY_full:'
print*, 'intA',norm2(intA)
print*, 'intB',norm2(intB)

deallocate(intsr,ints,work,tmp)
deallocate(calkA)
deallocate(intV,intAB)

end subroutine exind_A3_XY_full

subroutine app_exi_A3_XY_full(dim1,dim2,intA,intB,occA,occB,IRDM2Type,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: IRDM2Type
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX)
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intA(dim1),intB(dim2)
double precision,intent(in) :: AXELE,BXELE
double precision,intent(in) :: occA(NBas),occB(NBas)
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
integer :: it,is,id,ie
double precision :: AlphaA(NBas),AlphaB(NBas) 
double precision :: vAA,vBA,vAB,vBB
double precision :: fact,val
double precision :: Sxx(NBas,NBas),S10(NBas,NBas)
double precision :: TA(dimOA,dimOA),TB(dimOB,dimOB),U(dimOA,dimOB)
double precision :: intOAl(NBas,NBas),intOAr(NBas,NBas)
double precision :: intOBl(NBas,NBas),intOBr(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:),tmp(:,:)
double precision,allocatable :: intOA(:,:),intOB(:,:),  &
                                intVaa(:,:),intVbb(:,:),&
                                intVi(:,:),intVo(:,:)
double precision,external :: ddot

Sxx = 0d0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

AlphaA = 0d0
AlphaB = 0d0
if(IRDM2Type==0) then
   do i=1,dimOA
      AlphaA(i) = occA(i)
   enddo
   do i=1,dimOB
      AlphaB(i) = occB(i)
   enddo
elseif(IRDM2Type==1) then
   do i=1,dimOA
      AlphaA(i) = sqrt(occA(i))
   enddo
   do i=1,dimOB
      AlphaB(i) = sqrt(occB(i))
   enddo
endif

allocate(intOA(NBas,NBas),intOB(NBas,NBas))
allocate(intVaa(NBas,NBas),intVbb(NBas,NBas))
allocate(intVi(NBas,NBas),intVo(NBas,NBas))

S10 = 0d0
do ia=1,NBas
   do ip=1,dimOA
      S10(ip,ia) = AlphaA(ip)*Sab(ip,ia)
   enddo
enddo

! intermediate from overlap integrals
! intOA(a,a') = sum_b S^01(a,b).S^01(a',b)
intOA = 0d0
do ir=1,NBas
   do ip=1,NBas
      do ib=1,dimOB
         intOA(ip,ir) = intOA(ip,ir) + occB(ib)*Sab(ip,ib)*Sab(ir,ib)
      enddo
   enddo
enddo
! intOB(b,b') = sum_a S^10(a,b).S^10(a,b')
intOB = 0d0
do ib=1,NBas
   do ia=1,NBas
      do iq=1,dimOA
         intOB(ia,ib) = intOB(ia,ib) + occA(iq)*Sab(iq,ia)*Sab(iq,ib)
      enddo
   enddo
enddo

TB = 0d0
do ib=1,dimOB
   do ia=1,dimOB
      fact = AlphaB(ia)*AlphaB(ib)
      do iq=1,dimOA
         TB(ia,ib) = TB(ia,ib) + occA(iq)*fact*Sab(iq,ia)*Sab(iq,ib)
      enddo
   enddo
enddo
TA = 0d0
do iq=1,dimOA
   do ip=1,dimOA
      fact = AlphaA(ip)*AlphaA(iq)
      do ib=1,dimOB
         TA(ip,iq) = TA(ip,iq) + occB(ib)*fact*Sab(ip,ib)*Sab(iq,ib)
      enddo
   enddo
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))
allocate(tmp(NBas,NBas))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints   = 0d0
intVi  = 0d0
intVo  = 0d0
intVaa = 0d0
intVbb = 0d0
U      = 0d0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ib = l
      ia = k

      fact = 0.5d0*AlphaB(ia)*AlphaB(ib)

      ! A part
      if (k <= dimOB) then

         if (ia == ib) then

            do ir=1,dimOA
               do ip=1,NBas
                  intVaa(ip,ir) = intVaa(ip,ir) + occB(ib)*ints(ip,ir)
               enddo
            enddo

         endif

         ! 1-1b
         do i=1,ADimX

            ip = AIndN(1,i)
            ir = AIndN(2,i)
            ipr = posA(ip,ir)

            val = 0.5d0*(occA(ir)-occA(ip))

            intA(ipr) = intA(ipr) - val*ints(ip,ir)*TB(ia,ib)

         enddo

         do iq=1,dimOA

            ! for 2-1b
            U(ia,ib) = U(ia,ib) + fact*occA(iq)*ints(iq,iq)

            ! outer: Sum_{a'b'} (aa'|b'b) sqrt(n_a') sqrt(n_b') S(a',b')
            ! used in A part
            do ip=1,NBas
               intVo(ip,ib) = intVo(ip,ib) + AlphaA(iq)*AlphaB(ia)*ints(ip,iq)*Sab(iq,ia)
            enddo

         enddo

      endif ! k <=dimOB

      ! B part
      do ip=1,dimOA
         intVbb(ia,ib) = intVbb(ia,ib) + occA(ip)*ints(ip,ip)
      enddo

      ! inner: Sum_{a'b'} (a'a|bb') sqrt(n_a') sqrt(n_b') S(a',b')
      ! used in B part
      do iq=1,dimOA
         do ip=1,dimOA
            intVi(iq,ia) = intVi(iq,ia) + AlphaA(ip)*AlphaB(ib)*ints(ip,iq)*Sab(ip,ib)
         enddo
      enddo

      ic = l
      ia = k

      iac = posB(ia,ic)

      ! (B) 1-1b
      if(iac/=0) then

         val = 0.5d0*(occB(ic)-occB(ia))

         do iq=1,dimOA
            do ip=1,dimOA
               intB(iac) = intB(iac) - val*ints(ip,iq)*TA(ip,iq)
            enddo
         enddo

      endif

   enddo
enddo

!print*, 'A:11b',norm2(intA)
!print*, 'B:11b',norm2(intB)

close(iunit)

! A
vAA = 0d0
vBA = 0d0
do iq=1,dimOA
   vAA = vAA + occA(iq)*intOA(iq,iq)
   vBA = vBA + occA(iq)*intVaa(iq,iq)
enddo
! B
vAB = 0d0
vBB = 0d0
do ib=1,dimOB
   vAB = vAB + occB(ib)*intOB(ib,ib)
   vBB = vBB + occB(ib)*intVbb(ib,ib)
enddo
! A
intOAl = 0d0
intOAr = 0d0
do ir=1,NBas
   do iq=1,dimOA
      intOAl(iq,ir) = AlphaA(iq)*intOA(iq,ir)
      intOAr(ir,iq) = AlphaA(iq)*intOA(ir,iq)
   enddo
enddo
! B
intOBl = 0d0
intOBr = 0d0
do ic=1,NBas
   do ib=1,dimOB
      intOBl(ib,ic) = AlphaB(ib)*intOB(ib,ic)
      intOBr(ic,ib) = AlphaB(ib)*intOB(ic,ib)
   enddo
enddo

! A part
do i=1,ADimX

   ip = AIndN(1,i)
   ir = AIndN(2,i)
   ipr = posA(ip,ir)

   fact = occA(ir) - occA(ip)

   ! 1-1a
   intA(ipr) = intA(ipr) + fact*intVaa(ip,ir)*vAA

   ! 2-1a
   intA(ipr) = intA(ipr) + fact*intOA(ip,ir)*vBA

   ! 1-1b is done in ints loop

   ! 2-1b
   val = 0d0
   do ib=1,dimOB
      do ia=1,dimOB
         val = val - U(ia,ib)*Sab(ip,ia)*Sab(ir,ib)
      enddo
   enddo
   intA(ipr) = intA(ipr) + fact*val

enddo
!print*, 'A:11a:21a:21b',norm2(intA)

! B part
call dgemm('N','N',NBas,NBas,NBas,1d0,intVaa,NBas,S10,NBas,0d0,tmp,NBas)
do i=1,BDimX

   ia = BIndN(1,i)
   ic = BIndN(2,i)
   iac = posB(ia,ic)

   fact = occB(ic) - occB(ia)

   ! (B) 1-1a
   intB(iac) = intB(iac) + fact*intVbb(ia,ic)*vAB

   ! (B) 2-1a
   intB(iac) = intB(iac) + fact*intOB(ia,ic)*vBB

   ! (B) 1-1b is done in ints loop

   ! (B) 2-1b
   !call dgemm('N','N',NBas,NBas,NBas,fact,intVaa,NBas,S10,NBas,0d0,tmp,NBas)
   !call dgemm('T','N',NBas,NBas,NBas,-0.5d0,S10,NBas,tmp,NBas,1d0,intB,NBas)
   val = 0d0
   do ip=1,dimOA
      val = val + S10(ip,ia)*tmp(ip,ic)
   enddo
   intB(iac) = intB(iac) -0.5d0*fact*val

enddo
!print*, 'B:11a:21a:21b',norm2(intB)

! A part : sqrt(n_r) - sqrt(n_p)
do i=1,ADimX

   ip = AIndN(1,i)
   ir = AIndN(2,i)
   ipr = posA(ip,ir)

   fact = AlphaA(ip) - AlphaA(ir)

   ! 3-1a
   intA(ipr) = intA(ipr) + 0.5d0*fact*dot_product(intVaa(ip,:),intOAl(:,ir))
   !!intA(ipr) = intA(ipr) + 0.5d0*fact*ddot(NBas,intVaa(ip,1),NBas,intOAl(:,ir),1)

   ! 4-1a
   intA(ipr) = intA(ipr) + 0.5d0*fact*dot_product(intOAr(ip,:),intVaa(:,ir))
   !!intA(ipr) = intA(ipr) + 0.5d0*fact*ddot(NBas,intOAr(ip,1),NBas,intVaa(:,ir),1)

   ! 3-1b
   val = 0d0
   do ib=1,dimOB
       val = val + AlphaB(ib)*intVo(ip,ib)*Sab(ir,ib)
   enddo
   intA(ipr) = intA(ipr) - 0.25d0*fact*val

   ! 4-1b
   val = 0d0
   do ib=1,dimOB
       val = val + AlphaB(ib)*intVo(ir,ib)*Sab(ip,ib)
   enddo
   intA(ipr) = intA(ipr) - 0.25d0*fact*val

enddo
!print*, 'A:31a:41a:31b:41b',norm2(intA)
!print*, 'A:31b:41b',norm2(intA)
!print*, 'A:31a:41a',norm2(intA)

! B part : sqrt(n_a) - sqrt(n_c)
do i=1,BDimX

   ia = BIndN(1,i)
   ic = BIndN(2,i)
   iac = posB(ia,ic)

   fact = AlphaB(ia) - AlphaB(ic)

   ! (B) 3-1a
   intB(iac) = intB(iac) + 0.5d0*fact*dot_product(intVbb(ia,:),intOBl(:,ic))
   !intB(iac) = intB(iac) + 0.5d0*fact*ddot(NBas,intVbb(ia,1),NBas,intOBl(:,ic),1)

   ! (B) 4-1a
   intB(iac) = intB(iac) + 0.5d0*fact*dot_product(intObr(ia,:),intVbb(:,ic))
   !intB(iac) = intB(iac) + 0.5d0*fact*ddot(NBas,intObr(ia,1),NBas,intVbb(:,ic),1)

   ! (B) 3-1b
   val = 0d0
   do iq=1,dimOA
      val = val + AlphaA(iq)*intVi(iq,ia)*Sab(iq,ic)
   enddo
   intB(iac) = intB(iac) - 0.25d0*fact*val

   ! (B) 4-1b
   val = 0d0
   do iq=1,dimOA
      val = val + AlphaA(iq)*intVi(iq,ic)*Sab(iq,ia)
   enddo
   intB(iac) = intB(iac) - 0.25d0*fact*val

enddo
!print*, 'B:31a:41a:31b:41b',norm2(intB)
!print*, 'B:31b:41b',norm2(intB)
!print*, 'B:31a:41a',norm2(intB)

intA = -8d0*intA
intB = -8d0*intB

!print*, ''
!print*, 'app_exi_A3_XY_full:'
!print*, 'intA',norm2(intA)
!print*, 'intB',norm2(intB)

deallocate(intOB,intOA)
deallocate(intVo,intVi,intVbb,intVaa)
deallocate(tmp,ints,work)

end subroutine app_exi_A3_XY_full

subroutine exind_A3_XY(dim1,dim2,intX,intY,RDM2Aval,RDM2Bval,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX)
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intX(dim1),intY(dim2)
double precision,intent(in) :: AXELE,BXELE
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)
double precision,intent(in) :: RDM2Aval(dimOA,dimOA,dimOA,dimOA),RDM2Bval(dimOB,dimOB,dimOb,dimOB)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iae,ipt
integer :: it,is,id,ie
double precision :: factA,factB,fact,valX,valY
double precision :: Sxx(NBas,NBas),Sba(NBas,NBas),&
                    ABmat(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),&
                                ints(:,:),intsr(:,:),tints(:,:)
double precision,allocatable :: intA(:,:,:,:),intB(:,:,:,:),  &
                                intAB(:,:,:,:),intV(:,:,:,:), &
                                calkA(:,:,:,:),calkB(:,:,:,:),&
                                tmpA(:,:,:,:),tmpB(:,:,:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

Sba = transpose(Sab)

allocate(intAB(dimOA,dimOA,dimOB,dimOB),intV(dimOA,dimOA,dimOB,dimOB))
!,intB(dimOB,dimOB,dimOA,NBas),intA(dimOA,dimOA,dimOB,NBas),

allocate(calkA(dimOA,dimOA,dimOB,NBas))
!,calkB(dimOB,dimOB,dimOA,NBas))

allocate(tmp(dimOA**2,dimOB**2))
allocate(work(NBas*NBas),ints(NBas,NBas),tints(dimOA,NBas),intsr(dimOA,dimOA))

! create (qs,bd) intermediate
intAB = 0
do id=1,dimOB
   do ib=1,dimOB
      call dger(dimOA,dimOA,1d0,Sab(:,id),1,Sab(:,ib),1,intAB(:,:,ib,id),dimOA)
   enddo
enddo

! intAB: RDM2A(tr,qs).I(qs,bd).RDM2B(ec,bd)^T
call dgemm('N','N',dimOA**2,dimOB**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intAB,dimOA**2,0d0,tmp,dimOA**2)
call dgemm('N','N',dimOA**2,dimOB**2,dimOB**2,1d0,tmp,dimOA**2,RDM2Bval,dimOB**2,0d0,intAB,dimOA**2)


deallocate(tmp)
allocate(tmp(dimOA,dimOA))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
calkA = 0
!calkB = 0
intV = 0
tints = 0
ints = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      if(k<=dimOB) then

         ic = l 
         ia = k

         ! calkA
         intsr(1:dimOA,1:dimOA) = ints(1:dimOA,1:dimOA)
         call dgemv('T',dimOA**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intsr,1,0d0,calkA(:,:,ic,ia),1)

         ! intV
         !if(ia<=dimOB) then
           call dger(dimOA**2,dimOB**2,1d0,calkA(:,:,ic,ia),1,RDM2Bval(:,:,ia,ic),1,intV,dimOA**2)
         !endif 

         do i=1,ADimX

            ip = AIndN(1,i)
            it = AIndN(2,i)
            ipt = posA(ip,it)

            valX = 0
            valY = 0
            if(ip<=dimOA) then

               !X: 3-1' /  Y: 1-1'
               do ir=1,dimOA 
                 valX = valX - intAB(ir,ip,ia,ic)*ints(ir,it)
                 valY = valY + intAB(ip,ir,ia,ic)*ints(it,ir)
               enddo

            endif

            !X: 1-1' / Y: 3-1'
            do ir=1,dimOA
               valX = valX + intAB(it,ir,ia,ic)*ints(ip,ir)
               valY = valY - intAB(ir,it,ia,ic)*ints(ip,ir)
            enddo

            intX(ipt) = intX(ipt) + valX
            intY(ipt) = intY(ipt) + valY

         enddo

      endif

   enddo
enddo

close(iunit)

!print*, 'intV',norm2(intV)

do i=1,ADimX

   ip = AIndN(1,i)
   it = AIndN(2,i)
   ipt = posA(ip,it)

   valX = 0
   valY = 0
   if(ip<=dimOA) then

     do id=1,dimOB
        do ib=1,dimOB

            ! 4-1'(X)
            do iq=1,dimOA 
               valX = valX - intV(iq,ip,ib,id)*Sab(it,ib)*Sab(iq,id)
            enddo

            !2-1'(Y)
            do is=1,dimOA 
               valY = valY + intV(ip,is,ib,id)*Sab(it,id)*Sab(is,ib)
            enddo

        enddo
     enddo

   endif

   do id=1,dimOB
      do ib=1,dimOB

         ! 2-1'(X)
         do is=1,dimOA
            valX = valX + intV(it,is,ib,id)*Sab(is,ib)*Sab(ip,id)
         enddo

         ! 4-1'(Y)
         do iq=1,dimOA
            valY = valY - intV(iq,it,ib,id)*Sab(iq,id)*Sab(ip,ib)
         enddo

      enddo 
   enddo 

   intX(ipt) = intX(ipt) + valX
   intY(ipt) = intY(ipt) + valY

enddo

intX = -2d0*intX
intY = -2d0*intY

print*, 'intX',norm2(intX)
print*, 'intY',norm2(intY)

deallocate(intsr,tints,ints,work,tmp)
deallocate(calkA)
deallocate(intV,intAB)

end subroutine exind_A3_XY

subroutine exind_A2_XX(dim1,dim2,intXA,intXB,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,PAaa,IntKFile,&
                    BOcc,AOcc,BIndN,AIndN,posB,posA,&
                    dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intXA(dim1),intXB(dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas),&
                               PAaa(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val
double precision :: Sxx(NBas,NBas),Amat(dimOB,dimOA),Bmat(NBas,NBas)
double precision,allocatable :: tmp1(:,:,:,:),tmp2(:,:,:,:),&
                                work(:),workSq(:,:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp1(dimOB,dimOB,dimOB,NBas),tmp2(dimOB,dimOB,dimOB,NBas),&
         work(NBas*NBas),workSq(NBas,NBas),ints(NBas,NBas))

! reduce dims of tmp1(dimOB,dimOB,dimOB,dimOA)?
call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp1,dimOB**3)

! P(p,w)^T.Sab(p,s)
call dgemm('T','N',NBas,NBas,NBas,1d0,PAaa,NBas,Smat,NBas,0d0,workSq,NBas)
!! test for o3v matrix
!tmp2=0
!do is=1,dimOB
!   do it=1,NBas
!      do iq=1,dimOB
!         do iw=1,dimOB
!            do ir=1,dimOB
!               tmp2(ir,iw,iq,it) = tmp2(ir,iw,iq,it) &
!                                + RDM2val(ir,iw,iq,is)*ints(it,is) !Smat(it,is)
!!!                                !+
!!!                                FRDM2(ir,iq,iw,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)*Smat(it,is)
!            enddo
!         enddo
!      enddo
!   enddo
!enddo
!print*, 'tmp2-1', norm2(tmp2)
!call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,workSq,NBas,0d0,tmp2,dimOB**3)
call dgemm('N','T',dimOB**3,NBas,dimOB,2d0,RDM2val,dimOB**3,workSq,NBas,0d0,tmp2,dimOB**3)

!do is=NBasis
!   do iw=1,dimOA
!      vec(is) = vec(is) + AOcc(iw)*Smat(iw,is)
!   enddo
!enddo
!tmp2=0
!do is=dimOB
!   tmp2() = tmpRDM2val(iu,ir,iq,is)*vec(is)
!enddo
!call dgemv('N',)...

!(FO|FO):(BB|BA) or (AA|AB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
Amat = 0
ints = 0
kl = 0
do l=1,dimOA
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(l,k)/AXELE,Vabb,Vbab(l,k)/BXELE,Sxx)

      if(k<=dimOB) then

         iu = l
         iq = k

         do it=1,NBas

            itu = posA(it,iu)

            fact = 2d0*(AOcc(iu) - AOcc(it))

            if(itu/=0) then

              ! XA-1
              do ir=1,dimOB 
                 do ip=1,dimOB
                    intXA(itu) = intXA(itu) - fact*tmp1(ip,ir,iq,it)*ints(ip,ir)
                 enddo
              enddo

            endif

         enddo

         iw = l
         iq = k

         ! XB-1,3
         do i=1,BDimX

            it = BIndN(1,i)
            iu = BIndN(2,i)
            itu = posB(it,iu)

            !fact = 2d0

            if(it<=dimOB) then

              ! 3
              do ir=1,dimOB
                 !intXB(itu) = intXB(itu) + fact*tmp2(ir,it,iq,iw)*ints(ir,iu)
                 intXB(itu) = intXB(itu) + tmp2(ir,it,iq,iw)*ints(ir,iu)
              enddo

            endif

            ! 1 
            do ir=1,dimOB
               !intXB(itu) = intXB(itu) - fact*tmp2(iu,ir,iq,iw)*ints(it,ir)
               intXB(itu) = intXB(itu) - tmp2(iu,ir,iq,iw)*ints(it,ir)
            enddo

         enddo

         ! 4
         !print*, 'XA-4'
         do it=1,dimOB
            do ir=1,dimOB
               do is=1,dimOB
                  Amat(it,iw) = Amat(it,iw) + RDM2val(is,ir,iq,it)*ints(is,ir)
               enddo
            enddo
         enddo

      endif

    iw = l
    it = k

    ! XB-2
    do iu=1,dimOB
 
       itu=posB(it,iu)

       !fact = 2d0*(BOcc(iu) - BOcc(it))
       fact = 2d0

       if(itu/=0) then

          do ir=1,dimOB
             do iq=1,dimOB
                !intXB(itu) = intXB(itu) - fact*tmp2(iq,ir,iu,iw)*ints(iq,ir)
                intXB(itu) = intXB(itu) - tmp2(iq,ir,iu,iw)*ints(iq,ir)
             enddo
          enddo

       endif

    enddo 

   enddo
enddo
close(iunit)

!4-th term
do j=1,BDimX

   it = BIndN(1,j)
   iu = BIndN(2,j)
   itu = posB(it,iu)
   
   !fact = 2d0*(BOcc(iu) - BOcc(it)) 
   fact = 2d0

   if(it<=dimOB) then
   
      do iw=1,dimOA
         intXB(itu) = intXB(itu) + fact*Amat(it,iw)*workSq(iw,iu)
      enddo

   endif
enddo

!!print*, 'workSq,Amat',norm2(workSq),norm2(Amat)
print*, ''
print*, 'A2 intXA',norm2(intXA)
print*, 'A2 intXB',norm2(intXB)

deallocate(ints,work,workSq,tmp2,tmp1)

end subroutine exind_A2_XX

subroutine app_exi_A2_XX(dim1,dim2,intXA,intXB,IRDM2Type,Smat, &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
                    occB,occA,BIndN,AIndN,posB,posA,&
                    dimOB,dimOA,BDimX,ADimX,NBas)
!
! compute A2A and A2B contributions to the X part
!
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
integer,intent(in) :: IRDM2Type
character(*) :: IntKFile
double precision,intent(inout) :: intXA(dim1),intXB(dim2)
double precision,intent(in) :: AXELE,BXELE,occA(NBas),occB(NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,ir,ipr,ia,ib,ic,iac
double precision :: AlphaA(NBas),AlphaB(NBas)
double precision :: fact,val
double precision :: Sxx(NBas,NBas)
double precision :: S20(NBas,NBas),S21(NBas,NBas),S22(NBas,NBas)
double precision :: intR(NBas,NBas),intO(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

AlphaA = 0d0
AlphaB = 0d0
if(IRDM2Type==0) then
   do i=1,dimOA
      AlphaA(i) = occA(i)
   enddo
   do i=1,dimOB
      AlphaB(i) = occB(i)
   enddo
elseif(IRDM2Type==1) then
   do i=1,dimOA
      AlphaA(i) = sqrt(occA(i))
   enddo
   do i=1,dimOB
      AlphaB(i) = sqrt(occB(i))
   enddo
endif

S20 = 0d0
S21 = 0d0
S22 = 0d0
do ib=1,dimOB
   do ia=1,dimOA
      S20(ia,ib) = occA(ia)*Smat(ia,ib)
      S21(ia,ib) = occA(ia)*AlphaB(ib)*Smat(ia,ib)
      S22(ia,ib) = occA(ia)*occB(ib)*Smat(ia,ib)
   enddo
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|BA) or (AA|AB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
intR = 0d0
intO = 0d0
ints = 0d0
kl = 0
do l=1,dimOA
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(l,k)/AXELE,Vabb,Vbab(l,k)/BXELE,Sxx)

      ia = l ! A
      ib = k ! B

      !if(ib<=dimOB) then
      do ic=1,dimOB
         intR(ib,ia) = intR(ib,ia) + occB(ic)*ints(ic,ic)
      enddo

      if(ib<=dimOB) then
      do ic=1,NBas
         intO(ic,ia) = intO(ic,ia) + AlphaB(ib)*ints(ic,ib)
      enddo
      endif

      ip = l ! A
      ib = k ! B

      ! 3b
      do i=1,BDimX

         ia = BIndN(1,i)
         ic = BIndN(2,i)
         iac = posB(ia,ic)

         fact = occB(ic) - occB(ia)

         intXB(iac) = intXB(iac) - 4d0*fact*ints(ia,ic)*S22(ip,ib)

      enddo

      ip = l ! A
      ia = k ! B

      ! 6b
      do ic=1,dimOB

         iac = posB(ia,ic)

         if(iac/=0) then

            fact = AlphaB(ic) - AlphaB(ia)
            
            val = 0d0
            do ib=1,dimOB
               val = val + ints(ib,ic)*S21(ip,ib)
            enddo

            intXB(iac) = intXB(iac) + 2.0d0*fact*val

         endif

      enddo


   enddo
enddo
close(iunit)

! intXA part
do i=1,ADimX

   ip  = AIndN(1,i)
   ir  = AIndN(2,i)
   ipr = posA(ip,ir)

   fact = occA(ir)-occA(ip)

   ! 1a
   val = 0d0
   do ib=1,dimOB
      val = val + intR(ib,ir)*occB(ib)*Smat(ip,ib) 
   enddo
   intXA(ipr) = intXA(ipr) - 4d0*fact*val

   ! 2a
   val = 0d0
   do ia=1,dimOB
      val = val + intO(ia,ir)*AlphaB(ia)*Smat(ip,ia)
   enddo
   intXA(ipr) = intXA(ipr) + 2.0d0*fact*val

enddo

! intXB part
do i=1,BDimX

   ia = BIndN(1,i)
   ic = BIndN(2,i)
   iac = posB(ia,ic)

   fact = occB(ic)-occB(ia)

   ! 4b
   val = 0d0
   do ip=1,dimOA 
      val = val + intR(ia,ip)*S20(ip,ic)
   enddo
   intXB(iac) = intXB(iac) - 4d0*fact*val

   fact = AlphaB(ic)-AlphaB(ia)

   ! 5b
   val = 0d0
   do ip=1,dimOA 
      val = val + intO(ia,ip)*S20(ip,ic)
   enddo
   intXB(iac) = intXB(iac) + 2.0d0*fact*val

enddo

!print*, ''
!print*, 'app_exi_A2_XX '
!print*, 'intXA',norm2(intXA)
!print*, 'intXB',norm2(intXB)

!intXA = 0d0
!intXB = 0d0

deallocate(work,ints)

end subroutine app_exi_A2_XX

subroutine exind_A2_YY(dim1,dim2,intYA,intYB,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,PAaa,IntKFile,&
                    BOcc,AOcc,BIndN,AIndN,posB,posA,&
                    dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intYA(dim1),intYB(dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas),&
                               PAaa(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val
double precision :: Sxx(NBas,NBas),Amat(dimOB,dimOA),Bmat(NBas,NBas)
double precision,allocatable :: tmp1(:,:,:,:),tmp2(:,:,:,:),&
                                work(:),workSq(:,:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp1(dimOB,dimOB,dimOB,NBas),tmp2(dimOB,dimOB,dimOB,NBas),&
         work(NBas*NBas),workSq(NBas,NBas),ints(NBas,NBas))

! reduce dims of tmp1(dimOB,dimOB,dimOB,dimOA)?
call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp1,dimOB**3)

! P(p,w)^T.Sab(p,s)
call dgemm('T','N',NBas,NBas,NBas,1d0,PAaa,NBas,Smat,NBas,0d0,workSq,NBas)
!call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,workSq,NBas,0d0,tmp2,dimOB**3)
call dgemm('N','T',dimOB**3,NBas,dimOB,2d0,RDM2val,dimOB**3,workSq,NBas,0d0,tmp2,dimOB**3)
!print*,'tmp2-2',norm2(tmp2)

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
Amat = 0
ints = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)

      iq = l
      it = k

      do iu=1,dimOA

         itu = posA(it,iu)

         fact = 2d0*(AOcc(it) - AOcc(iu))

         if(itu/=0) then

            do ir=1,dimOB
               do ip=1,dimOB
                  intYA(itu) = intYA(itu) - fact*tmp1(ip,ir,iq,iu)*ints(ip,ir)
               enddo
            enddo

         endif

      enddo

      if(k<=dimOA) then

         iq = l
         iw = k

         do j=1,BDimX

            it = BIndN(1,j)
            iu = BIndN(2,j)
            itu = posB(it,iu)

            !fact = 2d0

            if(it<=dimOB) then

               ! YB-1
               do ir=1,dimOB
                  !intYB(itu) = intYB(itu) - fact*tmp2(it,ir,iq,iw)*ints(iu,ir)
                  intYB(itu) = intYB(itu) - tmp2(it,ir,iq,iw)*ints(iu,ir)
               enddo

            endif

            ! YB-3
            do ir=1,dimOB
               !intYB(itu) = intYB(itu) + fact*tmp2(ir,iu,iq,iw)*ints(it,ir)
               intYB(itu) = intYB(itu) + tmp2(ir,iu,iq,iw)*ints(it,ir)
            enddo

         enddo

         ! YB-4
         do iu=1,dimOB
            do ir=1,dimOB
               do is=1,dimOB
                  Amat(iu,iw) = Amat(iu,iw) + RDM2val(is,ir,iq,iu)*ints(is,ir)
               enddo
            enddo
         enddo

         iu = l
         iw = k

         ! YB-2
         do it=1,dimOB

            itu = posB(it,iu)

            if(itu/=0) then

               do ir=1,dimOB
                  do iq=1,dimOB
                     !intYB(itu) = intYB(itu) - fact*tmp2(iq,ir,it,iw)*ints(iq,ir)
                     intYB(itu) = intYB(itu) - tmp2(iq,ir,it,iw)*ints(iq,ir)
                  enddo
               enddo

            endif

         enddo

      endif

   enddo
enddo

close(iunit)

!YB-4
do j=1,BDimX

   it = BIndN(1,j)
   iu = BIndN(2,j)
   itu = posB(it,iu)

   fact = 2d0

   do iw=1,dimOA
     intYB(itu) = intYB(itu) + fact*Amat(iu,iw)*workSq(iw,it)
   enddo

enddo

print*,  'intYA',norm2(intYA)
print*,  'intYB',norm2(intYB)

deallocate(ints,workSq,work,tmp2,tmp1)

end subroutine exind_A2_YY

subroutine app_exi_A2_YY(dim1,dim2,intYA,intYB,IRDM2Type,Smat, &
                         AXELE,Vabb,BXELE,Vbab,IntKFile, &
                         occB,occA,BIndN,AIndN,posB,posA,&
                         dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
integer,intent(in) :: IRDM2Type
character(*) :: IntKFile
double precision,intent(inout) :: intYA(dim1),intYB(dim2)
double precision,intent(in) :: AXELE,BXELE,occA(NBas),occB(NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,ir,ipr,ia,ib,ic,iac
double precision :: AlphaA(NBas),AlphaB(NBas)
double precision :: fact,val
double precision :: Sxx(NBas,NBas)
double precision :: S20(NBas,NBas),S21(NBas,NBas),S22(NBas,NBas)
double precision :: intR(NBas,NBas),intI(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

AlphaA = 0d0
AlphaB = 0d0
if(IRDM2Type==0) then
   do i=1,dimOA
      AlphaA(i) = occA(i)
   enddo
   do i=1,dimOB
      AlphaB(i) = occB(i)
   enddo
elseif(IRDM2Type==1) then
   do i=1,dimOA
      AlphaA(i) = sqrt(occA(i))
   enddo
   do i=1,dimOB
      AlphaB(i) = sqrt(occB(i))
   enddo
endif

S20 = 0d0
S21 = 0d0
S22 = 0d0
do ib=1,NBas
   do ia=1,dimOA
      S20(ia,ib) = occA(ia)*Smat(ia,ib)
   enddo
enddo
do ib=1,dimOB
   do ia=1,dimOA
      S21(ia,ib) = occA(ia)*AlphaB(ib)*Smat(ia,ib)
      S22(ia,ib) = occA(ia)*occB(ib)*Smat(ia,ib)
   enddo
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
intI = 0d0
intR = 0d0
ints = 0d0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)

      ib = l
      ia = k

      do ic=1,dimOB 
         intR(ia,ib) = intR(ia,ib) + occB(ic)*ints(ic,ic)
      enddo

      do ic=1,dimOB
         intI(ic,ia) = intI(ic,ia) + AlphaB(ib)*ints(ib,ic)
      enddo

      ib = l ! B
      ip = k ! A

      ! 3b
      do i=1,BDimX

         ia = BIndN(1,i)
         ic = BIndN(2,i)
         iac = posB(ia,ic)

         fact = occB(ia) - occB(ic)

         intYB(iac) = intYB(iac) - 4d0*fact*ints(ia,ic)*S22(ip,ib)

      enddo

      ic = l ! B
      ip = k ! A

      ! 6b
      do ia=1,NBas

         iac = posB(ia,ic)

         if(iac/=0) then

            fact = AlphaB(ia) - AlphaB(ic)

            val = 0d0 
            do ib=1,dimOB
               val = val + ints(ia,ib)*S21(ip,ib)
            enddo

            intYB(iac) = intYB(iac) + 2.0d0*fact*val

         endif

      enddo

   enddo
enddo

close(iunit)

! intYA part
do i=1,ADimX

   ip  = AIndN(1,i)
   ir  = AIndN(2,i)
   ipr = posA(ip,ir)

   fact = occA(ip)-occA(ir)

   ! 1a 
   val = 0d0
   do ib=1,dimOB
      val = val + intR(ip,ib)*occB(ib)*Smat(ir,ib)
   enddo
   intYA(ipr) = intYA(ipr) - 4d0*fact*val

   ! 2a
   val = 0d0
   do ia=1,dimOB
      val = val + intI(ia,ip)*AlphaB(ia)*Smat(ir,ia)
   enddo
   intYA(ipr) = intYA(ipr) + 2.0d0*fact*val

enddo

! intYB part
do i=1,BDimX

   ia = BIndN(1,i)
   ic = BIndN(2,i)
   iac = posB(ia,ic)

   fact = occB(ia)-occB(ic)

   ! 4b
   val = 0d0
   do ip=1,dimOA
      val = val + intR(ip,ic)*S20(ip,ia)
   enddo
   intYB(iac) = intYB(iac) - 4d0*fact*val

   fact = AlphaB(ia)-AlphaB(ic)

   ! 5b
   val = 0d0
   do ip=1,dimOA
      val = val + intI(ic,ip)*S20(ip,ia)
   enddo
   intYB(iac) = intYB(iac) + 2.0d0*fact*val

enddo

!intYA = -4d0*intYA
!intYB = -4d0*intYB

!print*, ''
!print*, 'app_exi_A2_YY'
!print*, 'intYA',norm2(intYA)
!print*, 'intYB',norm2(intYB)

deallocate(ints,work)

end subroutine app_exi_A2_YY

subroutine exind_A1_AB(ANDimX,BNDimX,intXA,intXB,intYA,intYB,AXELE,BXELE,dimOA,dimOB, &
                       AOcc,BOcc,PAaa,PBbb,Vaab,Vbab,Sab,AIndN,BIndN,posA,posB,NBas)
implicit none
 
integer,intent(in)          :: ANDimX,BNDimX,dimOA,dimOB,NBas
double precision,intent(in) :: AXELE,BXELE
integer,intent(in)          :: posA(NBas,NBas),posB(NBas,NBas),&
                               AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in)    :: Sab(NBas,NBas),&
                                  PAaa(NBas,NBas),PBbb(NBas,NBas),&
                                  Vaab(NBas,NBas),Vbab(NBas,NBas),&
                                  AOcc(NBas),BOcc(NBas)
double precision,intent(inout) :: intXA(ANDimX),intXB(BNDimX),&
                                  intYA(ANDimX),intYB(BNDimX)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ipq,ir,is,irs
double precision :: fact,ielA,ielB
double precision,allocatable :: work(:),&
                                ints(:,:),intsr(:,:),tints(:,:)

ielA = 1d0/(AXELE)
ielB = 1d0/(BXELE)
!ielA = 1d0/(2d0*AXELE)
!ielB = 1d0/(2d0*BXELE)

allocate(work(NBas*NBas),ints(NBas,NBas))

!(FF|OO):(AB|AB)
open(newunit=iunit,file='FFOOABAB',status='OLD', &
    access='DIRECT',recl=8*NBas*NBas)

ints = 0
kl = 0
do l=1,dimOB
   do k=1,dimOA
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*NBas)

      call ints_modify(NBas,NBas,ints,NBas,work,&
                      Sab(k,l)*ielA,Vaab,Vbab(k,l)*ielB,Sab)

      ir = l
      iq = k

      ! XA
      do ip=1,NBas

         ipq = posA(ip,iq)

         fact = 2d0*(AOcc(ip)-AOcc(iq))

         if(ipq/=0) then

            !do is=1,dimOB
            !   intXA(ipq) = intXA(ipq) + fact*ints(ip,is)*PBbb(ir,is)
            !enddo
            !print*, 'XA'
            intXA(ipq) = intXA(ipq) + fact*BOcc(ir)*ints(ip,ir)

         endif 
        
      enddo

      ! YB
      is = l
      iq = k

      do ir=1,NBas

         irs = posB(ir,is)

         fact = 2d0*(BOcc(ir)-BOcc(is))

         if(irs/=0) then

            !do ip=1,dimOA
            !   intYB(irs) = intYB(irs) - fact*ints(ip,ir)*PAaa(ip,iq)
            !enddo
 
            !print*, 'YB'
            intYB(irs) = intYB(irs) - fact*AOcc(iq)*ints(iq,ir)

         endif 

     enddo

     ints=0
     call ints_modify(NBas,NBas,ints,NBas,work, &
                       Vaab(k,l)*ielA,Sab,Sab(k,l)*ielB,Vbab)

     is = l
     ip = k

     ! XB 
     do ir=1,NBas

        irs = posB(ir,is)

        fact = 2d0*(BOcc(ir)-BOcc(is))

        if(irs/=0) then

           !do iq=1,dimOA
           !   intXB(irs) = intXB(irs) + fact*ints(iq,ir)*PAaa(ip,iq)
           !enddo

           !print*, 'XB'
           intXB(irs) = intXB(irs) + fact*AOcc(ip)*ints(ip,ir) 

        endif 

    enddo

    ! YA
    is = l
    iq = k

    do ip=1,NBas

        ipq = posA(ip,iq)

        fact = 2d0*(AOcc(ip)-AOcc(iq))

        if(ipq/=0) then

           !do ir=1,dimOB
           !   intYA(ipq) = intYA(ipq) - fact*ints(ip,ir)*PBbb(ir,is)
           !enddo
 
           !print*, 'YA'
           intYA(ipq) = intYA(ipq) - fact*BOcc(is)*ints(ip,is)

        endif 

    enddo

  enddo
enddo

!! test
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,dimOA
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*NBas)
!
!      call ints_modify(NBas,NBas,ints,NBas,work,&
!                       Vaab(k,l)/AXELE,Sab,Sab(k,l)/BXELE,Vbab)
!
!      !call ints_modify(NBas,NBas,ints,NBas,work,&
!      !                 Sab(k,l)/AXELE,Vaab,Vbab(k,l)/BXELE,Sab)
!
!      is = l
!      iq = k
!
!      do ip=1,NBas
!
!         ipq = posA(ip,iq)
!
!         fact = AOcc(ip)-AOcc(iq)
!
!         if(ipq/=0) then
!
!            do ir=1,dimOB
!               intXA(ipq) = intXA(ipq) + fact*ints(ip,ir)*PBbb(ir,is)
!            enddo
!
!         endif 
!        
!      enddo
!
!   enddo
!enddo
!! end test

close(iunit)

!print*, 'A1 intXA',norm2(intXA)
!print*, 'A1 intXB',norm2(intXB)
!print*, 'A1 intYA',norm2(intYA)
!print*, 'A1 intYB',norm2(intYB)

!! test-2
!open(newunit=iunit,file='FOFOABAB',status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
!intA = 0
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       Vaab(k,l)/AXELE,Sab,Sab(k,l)/BXELE,Vbab)
!
!      is = l
!      ip = k
!
!      ! XA
!      do iq=1,dimOA
!
!         ipq = posA(ip,iq)
!
!         fact = AOcc(ip)-AOcc(iq)
!
!         if(ipq/=0) then
!
!            do ir=1,dimOB
!               intXA(ipq) = intXA(ipq) + fact*ints(iq,ir)*PBbb(ir,is)
!            enddo
!
!         endif 
!        
!      enddo
!
!   enddo
!enddo
!close(iunit)
!print*, 'A1*intXA',norm2(intXA)

!!test-3
!intYb = 0
!open(newunit=iunit,file='FOFOABBA',status='OLD', &
!    access='DIRECT',recl=8*NBas*dimOB)
!
!ints = 0
!kl = 0
!do l=1,dimOA
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       !Vaab(k,l)/AXELE,Sab,Sab(k,l)/BXELE,Vbab)
!                       Vaab(l,k)/AXELE,Sab,Sab(l,k)/BXELE,Vbab)
!
!
!      ir = k
!      ip = l
!
!      do is=1,dimOB
!
!         fact = BOcc(ir) - BOcc(is)
!
!         irs = posB(ir,is)
!         if(irs/=0) then
!
!            do iq=1,dimOA
!               intYB(irs) = intYB(irs) - fact*ints(iq,is)*PAaa(ip,iq)
!            enddo
!
!         endif
!
!      enddo
!
!   enddo
!enddo
!close(iunit)
!
!intYB = 2d0*intYB
!print*, 'A1*intYB',norm2(intYB)

deallocate(ints,work)

end subroutine exind_A1_AB

end module exi

