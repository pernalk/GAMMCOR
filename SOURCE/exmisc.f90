module exmisc 

implicit none

contains

subroutine ints_modify(n1,n2,ints,n,raw_ints,alphaX,X,alphaY,Y)
implicit none

integer,intent(in) :: n1,n2,n
double precision,intent(out) :: ints(n1,n2)
double precision,intent(in) :: alphaX,alphaY,raw_ints(n,n),X(n,n),Y(n,n)

ints(1:n1,1:n2) = raw_ints(1:n1,1:n2) + alphaX*X(1:n1,1:n2) + alphaY*Y(1:n1,1:n2)
!ints(1:n1,1:n2) = raw_ints(1:n1,1:n2)

end subroutine ints_modify


!subroutine A1_Mod_1el(tmp,Vaab,Vbab,Sab,posA,posB,A,B,NBas,inter)
subroutine A1_Mod_1el(tmp,AXELE,BXELE,AOcc,BOcc,Vaab,Vbab,Sab,&
                      AIndN,BIndN,posA,posB,ANDimX,BNDimX,NBas,inter)

implicit none

!type(SystemBlock) :: A, B
integer,intent(in) :: ANDimX,BNDimX,NBas
double precision,intent(in) :: AXELE,BXELE
character(*),intent(in) :: inter
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas),&
                      AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in) :: Sab(NBas,NBas),&
                               Vaab(NBas,NBas),Vbab(NBas,NBas),&
                               AOcc(NBas),BOcc(NBas)
double precision,intent(inout) :: tmp(ANDimX,BNDimX)

integer :: i,j,ip,iq,ipq,ir,is,irs
double precision :: fact,ielA,ielB

 ielA = 1d0/(2d0*AXELE)
 ielB = 1d0/(2d0*BXELE)

select case(inter)
case('XX','xx')

   do j=1,BNDimX
   
      ir = BIndN(1,j)
      is = BIndN(2,j)
      irs = posB(ir,is) 
   
      do i=1,ANDimX
   
         ip = AIndN(1,i)
         iq = AIndN(2,i)
         ipq = posA(ip,iq)
   
         fact = -2d0*(AOcc(ip)-AOcc(iq)) * &
                     (BOcc(ir)-BOcc(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(ip,is)*Sab(iq,ir) + fact*ielB*Vbab(iq,ir)*Sab(ip,is)
   
      enddo
   enddo

case('YY','yy')

   do j=1,BNDimX
   
      ir = BIndN(1,j)
      is = BIndN(2,j)
      irs = posB(ir,is) 
   
      do i=1,ANDimX
   
         ip = AIndN(1,i)
         iq = AIndN(2,i)
         ipq = posA(ip,iq)
   
         fact = -2d0*(AOcc(ip)-AOcc(iq)) * &
                     (BOcc(ir)-BOcc(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(iq,ir)*Sab(ip,is) + fact*ielB*Sab(iq,ir)*Vbab(ip,is)
   
      enddo
   enddo

case('XY','xy')

   do j=1,BNDimX

      ir = BIndN(1,j)
      is = BIndN(2,j)
      irs = posB(ir,is) 
 
      do i=1,ANDimX
 
         ip = AIndN(1,i)
         iq = AIndN(2,i)
         ipq = posA(ip,iq)

         fact = 2d0*(AOcc(ip)-AOcc(iq)) * &
                    (BOcc(ir)-BOcc(is))
   
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(ip,ir)*Sab(iq,is) + fact*ielB*Sab(ip,ir)*Vbab(iq,is)

      enddo
   enddo

case('YX','yx')

   do j=1,BNDimX

      ir = BIndN(1,j)
      is = BIndN(2,j)
      irs = posB(ir,is) 
 
      do i=1,ANDimX
 
         ip = AIndN(1,i)
         iq = AIndN(2,i)
         ipq = posA(ip,iq)

         fact = 2d0*(AOcc(ip)-AOcc(iq)) * &
                    (BOcc(ir)-BOcc(is))
 
         tmp(ipq,irs) = tmp(ipq,irs) + fact*ielA*Vaab(iq,is)*Sab(ip,ir) + fact*ielB*Sab(iq,is)*Vbab(ip,ir)
 
      enddo
   enddo

end select

end subroutine A1_Mod_1el

!subroutine make_sij_Y(sij,tmp1,A,B,nOVB,NBas)
subroutine make_sij_Y(sij,tmp1,AOcc,BOcc,AEigY,AEigX,BEigY,BEigX,&
                      Anum0,Bnum0,dimOA,dimOB,nOVB,AIndN,BIndN,ANDimX,BNDimX,NBas)
implicit none

!type(SystemBlock) :: A, B

integer,intent(in) :: Anum0,Bnum0,dimOA,dimOB,nOVB,ANDimX,BNDimX,NBas
integer,intent(in) :: AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in) :: AOcc(NBas),BOcc(NBas),&
                               AEigY(ANDimX*ANDimX),AEigX(ANDimX*ANDimX),&
                               BEigY(BNDimX*BNDimX),BEigX(BNDimX*BNDimX)
double precision,intent(inout) :: tmp1(ANDimX,BNDimX),&
                                  sij(ANDimX,BNDimX)

integer :: iunit
!integer :: dimOA,dimOB
integer :: i,j,ir,is,irs,ip,iq,ipq
double precision :: fact
double precision,allocatable :: work(:) 

!dimOA = A%INAct+A%NAct
!dimOB = B%INAct+B%NAct

allocate(work(nOVB))

open(newunit=iunit,file='TWOMOAB',status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*nOVB)

tmp1=0
do ipq=1,ANDimX
   ip = AIndN(1,ipq)
   iq = AIndN(2,ipq)
   read(iunit,rec=iq+(ip-Anum0-1)*dimOA) work(1:nOVB)
   do irs=1,BNDimX
      ir = BIndN(1,irs)
      is = BIndN(2,irs)

      fact = (AOcc(ip)-AOcc(iq)) * &
             (BOcc(ir)-BOcc(is)) * &
              work(is+(ir-Bnum0-1)*dimOB)

      do i=1,ANDimX
         tmp1(i,irs) = tmp1(i,irs) + &
                      fact * &
                      (AEigY(ipq+(i-1)*ANDimX)-AEigX(ipq+(i-1)*ANDimX))
      enddo

   enddo
enddo
close(iunit)

sij=0
do j=1,BNDimX
   do i=1,ANDimX
      do irs=1,BNDimX
         ir = BIndN(1,irs)
         is = BIndN(2,irs)
         sij(i,j) = sij(i,j) + &
                    (BEigY(irs+(j-1)*BNDimX)-BEigX(irs+(j-1)*BNDimX))*tmp1(i,irs)
      enddo
   enddo  
 enddo

deallocate(work)

end subroutine make_sij_Y

!subroutine make_tij_Y(tmp3,tmp2,tmp1,posA,posB,Sab,Sba,A,B,NBas)
subroutine make_tij_Y(tmp3,tmp2,tmp1,AOcc,BOcc,AEigY,AEigX,BEigY,BEigX,&
                      AIndN,BIndN,posA,posB,Sab,Sba,ANDimX,BNDimX,NBas)
implicit none

integer,intent(in) :: ANDimX,BNDimX,NBas
integer,intent(in) :: posA(NBas,NBas),posB(NBas,NBas),&
                      AIndN(2,ANDimX),BIndN(2,BNDimX)
double precision,intent(in) :: Sab(NBas,NBas),Sba(NBas,NBas)
double precision,intent(in) :: AOcc(NBas),BOcc(NBas),&
                               AEigY(ANDimX*ANDimX),AEigX(ANDimX*ANDimX),&
                               BEigY(BNDimX*BNDimX),BEigX(BNDimX*BNDimX)
double precision,intent(inout) :: tmp1(ANDimX,BNDimX),&
                                  tmp2(ANDimX,BNDimX),&
                                  tmp3(ANDimX,BNDimX)

integer :: i,j,ir,is,irs,ip,iq,ipq
double precision :: fact

tmp1=0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) + Sab(iq,ir)*Sba(is,ip)*fact

   enddo
enddo
! CPLD
!X_A.I.X_B
call dgemm('T','N',ANDimX,BNDimX,ANDimX,1d0,AEigX,ANDimX,tmp1,ANDimX,0d0,tmp2,ANDimX)
call dgemm('N','N',ANDimX,BNDimX,BNDimX,1d0,tmp2,ANDimX,BEigX,BNDimX,0d0,tmp3,ANDimX)
!print*, 'Xa.I.Xb:',norm2(tmp3)
!Y_A.I.Y_B
call dgemm('T','N',ANDimX,BNDimX,ANDimX,1d0,AEigY,ANDimX,tmp1,ANDimX,0d0,tmp2,ANDimX)
call dgemm('N','N',ANDimX,BNDimX,BNDimX,1d0,tmp2,ANDimX,BEigY,BNDimX,1d0,tmp3,ANDimX)
!print*, 'Ya.I.Yb:,test1',norm2(tmp3)

tmp1=0
do j=1,BNDimX

   ir = BIndN(1,j)
   is = BIndN(2,j)
   irs = posB(ir,is)

   do i=1,ANDimX

      ip = AIndN(1,i)
      iq = AIndN(2,i)
      ipq = posA(ip,iq)

      fact = (BOcc(ir)-BOcc(is)) * &
             (AOcc(ip)-AOcc(iq))

      tmp1(ipq,irs) = tmp1(ipq,irs) - Sab(ip,ir)*Sba(is,iq)*fact

   enddo
enddo
!X_A.I.Y_B
call dgemm('T','N',ANDimX,BNDimX,ANDimX,1d0,AEigX,ANDimX,tmp1,ANDimX,0d0,tmp2,ANDimX)
call dgemm('N','N',ANDimX,BNDimX,BNDimX,1d0,tmp2,ANDimX,BEigY,BNDimX,1d0,tmp3,ANDimX)
!Y_A.I.X_B
call dgemm('T','N',ANDimX,BNDimX,ANDimX,1d0,AEigY,ANDimX,tmp1,ANDimX,0d0,tmp2,ANDimX)
call dgemm('N','N',ANDimX,BNDimX,BNDimX,1d0,tmp2,ANDimX,BEigX,BNDimX,1d0,tmp3,ANDimX)
!print*, 'test1:',norm2(tmp3)

end subroutine make_tij_Y

!subroutine make_tind(tvec,Xmat,Ymat,Sab,Pmat,Occ,IndN,pos,NDimX,NBas)
!implicit none
!
!integer,intent(in) :: NBas,NDimX
!integer,intent(in) :: IndN(2,NDimX),pos(NBas,NBas)
!double precision,intent(in) :: Xmat(NDimX,NDimX),Ymat(NDimX,NDimX)
!double precision,intent(in) :: Sab(NBas,NBas),Pmat(NBas,NBas),&
!                               Occ(NBas)
!double precision,intent(inout) :: tvec(NDimX)
!
!integer :: i,ip,iq,ir,ipq
!integer :: ii,j,k
!double precision :: Sba(NBas,NBas)
!double precision,allocatable :: tmp1(:,:),tmp2(:)
!double precision :: fact
!
! Sba = transpose(Sab)
!
! allocate(tmp1(NBas,NBas),tmp2(NDimX))
!
! ! P(r,s).Sba(s,p)
! call dgemm('N','N',NBas,NBas,NBas,1d0,Pmat,NBas,Sba,NBas,0d0,tmp1,NBas)
!
! tmp2=0
! do i=1,NDimX
!    ip = IndN(1,i)
!    iq = IndN(2,i)
!    ipq = pos(ip,iq)
!
!    fact=(Occ(ip)-Occ(iq))
! 
!    do ir=1,NBas
!       tmp2(i) = tmp2(i) + fact*Sab(iq,ir)*tmp1(ir,ip)
!    enddo
!
! enddo
!
! call dgemv('N',NDimX,NDimX,1d0,Ymat-Xmat,NDimX,tmp2,1,0d0,tvec,1)
! !print*, 'tvec-2',norm2(tvec)
!
! deallocate(tmp2,tmp1)
!
!end subroutine make_tind

subroutine inter_A3_XY(dim1,dim2,intXY,AOcc,BOcc,RDM2Aval,RDM2Bval,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intXY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas)
double precision,intent(in) :: RDM2Aval(dimOA,dimOA,dimOA,dimOA),RDM2Bval(dimOB,dimOB,dimOb,dimOB)
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iae,ipt
integer :: it,is,id,ie
double precision :: factA,factB,fact,val,val2,valTr
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

allocate(intAB(dimOA,dimOA,dimOB,dimOB),intB(dimOB,dimOB,dimOA,NBas),&
         intA(dimOA,dimOA,dimOB,NBas),intV(dimOA,dimOA,dimOB,dimOB))
allocate(calkA(dimOA,dimOA,dimOB,NBas),calkB(dimOB,dimOB,dimOA,NBas))
allocate(tmp(dimOA**2,dimOB**2))
allocate(work(NBas*NBas),ints(NBas,NBas),tints(dimOA,NBas),intsr(dimOA,dimOA))

! I(qs,bd) 
intAB = 0
do id=1,dimOB
   do ib=1,dimOB
      call dger(dimOA,dimOA,1d0,Sab(:,id),1,Sab(:,ib),1,intAB(:,:,ib,id),dimOA)
   enddo
enddo

! intAB: RDM2A(tr,qs).I(qs,bd).RDMB(ec,bd)^T
call dgemm('N','N',dimOA**2,dimOB**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intAB,dimOA**2,0d0,tmp,dimOA**2)
call dgemm('N','N',dimOA**2,dimOB**2,dimOB**2,1d0,tmp,dimOA**2,RDM2Bval,dimOB**2,0d0,intAB,dimOA**2)
!print*, 'intAB', norm2(intAB)

! intB
allocate(tmpB(dimOB,dimOB,dimOA,NBas))
tmpB = 0
do ip=1,NBas
   do is=1,dimOA
      call dger(dimOB,dimOB,1d0,Sba(:,is),1,Sba(:,ip),1,tmpB(:,:,is,ip),dimOB)
   enddo
enddo
call dgemm('N','N',dimOB**2,dimOA*NBas,dimOB**2,1d0,RDM2Bval,dimOB**2,tmpB,dimOB**2,0d0,intB,dimOB**2)
!print*, 'intB:', norm2(intB)
deallocate(tmpB)

! intA
allocate(tmpA(dimOA,dimOA,dimOB,NBas))
tmpA = 0
do ia=1,NBas
   do id=1,dimOB
      call dger(dimOA,dimOA,1d0,Sab(:,id),1,Sab(:,ia),1,tmpA(:,:,id,ia),dimOA)
   enddo
enddo
call dgemm('N','N',dimOA**2,dimOB*NBas,dimOA**2,1d0,RDM2Aval,dimOA**2,tmpA,dimOA**2,0d0,intA,dimOA**2)
!print*, 'intA:', norm2(intA)

deallocate(tmpA)

deallocate(tmp)
allocate(tmp(dimOA,dimOA))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
calkA = 0
calkB = 0
intV  = 0
tints = 0
ints  = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ic = l 
      ia = k

      ! calkA
      intsr(1:dimOA,1:dimOA) = ints(1:dimOA,1:dimOA)
      call dgemv('T',dimOA**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intsr,1,0d0,calkA(:,:,ic,ia),1)

      if(ia<=dimOB) then
         
         do i=1,dimOB 
            do j=1,dimOB 
               do is=1,dimOA 
                  do ip=1,dimOA 
                     intV(ip,is,i,j) = intV(ip,is,i,j) + calkA(ip,is,l,k)*RDM2Bval(k,l,i,j)
                  enddo
               enddo
            enddo
         enddo

      endif

      do ie=1,dimOB

         iae = posB(ia,ie)

         if(iae/=0) then

            do i=1,ADimX

               ip = AIndN(1,i)
               it = AIndN(2,i)
               ipt = posA(ip,it)

               val = 0
               if(ip<=dimOA) then
                  ! 3-3'(0)
                  do ir=1,dimOA
                     val = val + intAB(ir,ip,ic,ie)*ints(ir,it)
                  enddo
               endif     

               ! 1-3'(0)
               do ir=1,dimOA
                  val = val - intAB(it,ir,ic,ie)*ints(ip,ir)
               enddo

               intXY(ipt,iae) = intXY(ipt,iae) + val

           enddo

         endif
      enddo

      if(k<=dimOB) then

         ic = l
         ib = k

         do ir=1,dimOA
            do ip=1,NBas
               tints(ir,ip) = ints(ip,ir)
            enddo
         enddo

         ! calkB
         call dger(dimOB**2,dimOA*NBas,1d0,RDM2Bval(:,:,ib,ic),1,tints,1,calkB,dimOB**2)

         ic = l
         ie = k

         do ia=1,dimOB

            iae = posB(ia,ie)

            if(iae/=0) then

               do i=1,ADimX

                  ip = AIndN(1,i)
                  it = AIndN(2,i)
                  ipt = posA(ip,it)

                  val = 0
                  ! 3-1'(0)
                  if(ip<=dimOA) then
                     do ir=1,dimOA
                        val = val - intAB(ir,ip,ia,ic)*ints(ir,it)
                     enddo
                  endif
 
                  ! 1-1'(0)
                  do ir=1,dimOA
                     val = val + intAB(it,ir,ia,ic)*ints(ip,ir)
                  enddo
                  
                  intXY(ipt,iae) = intXY(ipt,iae) + val

               enddo
            endif

         enddo
      endif

   enddo
enddo
close(iunit)

!print*, 'calkA',norm2(calkA)
!print*, 'calkB',norm2(calkB)
!print*, 'intV:', norm2(intV)

! 1 and 2
ABmat = 0
do j=1,BDimX

   ia = BIndN(1,j)
   ie = BIndN(2,j)
   iae = posB(ia,ie) 

   do i=1,ADimX

      ip = AIndN(1,i)
      it = AIndN(2,i)
      ipt = posA(ip,it) 

      val = 0
      ! 2-3'(1)
      do is=1,dimOA
         do ic=1,dimOB
            val = val - calkA(it,is,ic,ia)*intB(ic,ie,is,ip)
         enddo
      enddo
  
      ! 1-4'(1) - ?
       do ir=1,dimOA
         do ib=1,dimOB
            val = val - intA(ir,it,ib,ia)*calkB(ib,ie,ir,ip)
         enddo
      enddo

      ! 2-4'(2)
      do ib=1,dimOB
         do is=1,dimOA
            val = val - intV(it,is,ib,ie)*Sab(ip,ia)*Sab(is,ib)
         enddo
      enddo

      if(ia<=dimOB) then

         ! 2-1'(1)
         do ic=1,dimOB
            do is=1,dimOA
               val = val + calkA(it,is,ic,ie)*intB(ia,ic,is,ip)
            enddo
         enddo

         ! 1-2'(1)
         do id=1,dimOB
            do ir=1,dimOA
               val = val + intA(it,ir,id,ie)*calkB(ia,id,ir,ip)
            enddo
         enddo

         ! 2-2'(2) 
         do id=1,dimOB
            do is=1,dimOA
               val = val + intV(it,is,ia,id)*Sab(ip,id)*Sab(is,ie)
            enddo
         enddo

         if(ip<=dimOA) then

            ! 4-1'(1)
            do ic=1,dimOB
               do iq=1,dimOA
                 val = val - calkA(iq,ip,ic,ie)*intB(ia,ic,it,iq)
               enddo 
            enddo 

            !3-2'(1)
            do id=1,dimOB
               do ir=1,dimOA
                  val = val - intA(ir,ip,id,ie)*calkB(ia,id,it,ir)
               enddo 
            enddo 

            ! 4-2'(2)
            do id=1,dimOB
               do iq=1,dimOA
                  !ABmat(ip,ia) = ABmat(ip,ia) + intV(iq,ip,ia,id)*Sab(iq,id)
                  val = val - intV(iq,ip,ia,id)*Sab(iq,id)*Sab(it,ie)
               enddo
            enddo
!
!            !intXY(ipt,iae) = intXY(ipt,iae) - ABmat(ip,ia)*Sab(it,ie)
!            !intXY(ipt,iae) = intXY(ipt,iae) + val 
!
         endif
      endif

      if(ip<=dimOA) then

         ! 4-3'(1)
         do ic=1,dimOB
            do iq=1,dimOA
               val = val + calkA(iq,ip,ic,ia)*intB(ic,ie,it,iq)
            enddo 
         enddo 

         ! 3-4'(1)
         do ir=1,dimOA
            do ib=1,dimOB
               val = val + intA(ip,ir,ib,ia)*calkB(ib,ie,it,ir)
            enddo
         enddo

         ! 4-4'(2)
         do ib=1,dimOB
            do iq=1,dimOA
               val = val + intV(iq,ip,ib,ie)*Sab(iq,ia)*Sab(it,ib)
            enddo
         enddo

      endif

      intXY(ipt,iae) = intXY(ipt,iae) + val 
   
   enddo
enddo

intXY = -2d0*intXY
!print*, 'intXY-A3',norm2(intXY)
!print*, ''

deallocate(intsr,tints,ints,work)
deallocate(calkB,calkA)
deallocate(intV,intA,intB,intAB)

end subroutine inter_A3_XY

!subroutine app_A3_XY(dim1,dim2,intXY,AOcc,BOcc,AFmat,BFmat,Sab, &
!                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
!                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
!character(*) :: IntKFile
!double precision,intent(inout) :: intXY(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
!                               AFmat(NBas,NBas),BFmat(NBas,NBas),&
!                               Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)
!
!integer :: iunit
!integer :: i,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
!double precision :: factA,factB,fact,val,val2,valTr
!double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas),&
!                    Amat(NBas,NBas),Emat(NBas,dimOA)
!double precision,allocatable :: work(:),ints(:,:)
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(work(NBas*NBas),ints(NBas,NBas))
!
!Amat = 0
!do ib=1,dimOB
!   do ir=1,NBas
!      do ip=1,NBas
!         Amat(ip,ir) = Amat(ip,ir) + BOcc(ib)*Sab(ip,ib)*Sab(ir,ib)
!      enddo
!   enddo
!enddo
!
!valTr = 0
!do iq=1,dimOA
!   valTr = valTr + AOcc(iq)*Amat(iq,iq)
!enddo
!
!Sbb = 0
!do ic=1,NBas
!   do ia=1,NBas
!      do iq=1,dimOA
!         Sbb(ia,ic) = Sbb(ia,ic) + AOcc(iq)*Sab(iq,ia)*Sab(iq,ic)
!      enddo
!   enddo
!enddo
!
!!(FO|FO): (AA|BB)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOA)
!
!! one loop over integrals
!ints = 0
!Emat = 0
!kl = 0
!do l=1,dimOB
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOA)
!
!      call ints_modify(NBas,dimOA,ints,NBas,work,&
!                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)
!
!      ic = l
!      ia = k
!
!      iac = posB(ia,ic)
!
!      if(iac/=0) then
!
!         do i=1,ADimX
!            ip = AIndN(1,i)
!            ir = AIndN(2,i)
!            ipr = posA(ip,ir)
!            
!            fact  = 8d0*(AOcc(ir)-AOcc(ip))*(BOcc(ia)-BOcc(ic))
!            factB = 4d0*(BOcc(ia)-BOcc(ic))
!
!            ! 1A-1B
!            intXY(ipr,iac) = intXY(ipr,iac) - fact*valTr*ints(ip,ir)
!
!            val = 0
!            do iq=1,dimOA
!               val = val + AOcc(iq)*ints(iq,iq)  
!            enddo
!
!            !2A-1B 
!            intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Amat(ip,ir)
!
!            val = 0
!            do iq=1,dimOA
!               val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Amat(iq,ir)
!            enddo
!
!            !3A-1B
!            intXY(ipr,iac) = intXY(ipr,iac) - factB*val
!
!            val = 0
!            do iq=1,dimOA
!               val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Amat(ip,iq)
!            enddo
!
!            !4A-1B
!            intXY(ipr,iac) = intXY(ipr,iac) - factB*val
!
!         enddo
!
!      endif
!
!      ! all 2B terms 
!      if(k==l.and.k<=dimOB) then
!
!         ib = k
!         Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + BOcc(ib)*ints(1:NBas,1:dimOA)
!
!      endif
!
!      ib = l
!      ia = k 
!   
!      val = 0 
!      do iq=1,dimOA
!         val = val + AOcc(iq)*ints(iq,iq)
!      enddo
! 
!      do ic=1,dimOB
!
!         iac = posB(ia,ic)
!
!         if(iac/=0) then
!
!            factB = 2d0*(BFmat(ib,ic)-BFmat(ib,ia))
!
!            do i=1,ADimX
!
!               ip = AIndN(1,i)
!               ir = AIndN(2,i)
!               ipr = posA(ip,ir)
!
!               fact = 2d0*factB*(AOcc(ir)-AOcc(ip))
!
!               !1A-3B 
!               intXY(ipr,iac) = intXY(ipr,iac) - fact*Sbb(ic,ib)*ints(ip,ir)
!               
!               !2A-3B 
!               intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ib)*Sab(ir,ic)
!
!               val2 = 0
!               do iq=1,dimOA
!                  val2 = val2 + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ib)
!               enddo
!
!               !3A-3B 
!               intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ir,ic)
!
!               val2 = 0
!               do iq=1,dimOA
!                  val2 = val2 + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ic)
!               enddo
!
!               !4A-3B 
!               intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ip,ib)
!
!            enddo 
!         endif
!      enddo
!
!      if(k<=dimOB) then
!
!         ic = l
!         ib = k 
!
!         val = 0 
!         do iq=1,dimOA
!            val = val + AOcc(iq)*ints(iq,iq)
!         enddo
!
!         do ia=1,NBas
!
!            iac = posB(ia,ic)
!
!            if(iac/=0) then
!
!               factB = 2d0*(BFmat(ic,ib)-BFmat(ia,ib))
!
!               do i=1,ADimX
!
!                  ip = AIndN(1,i)
!                  ir = AIndN(2,i)
!                  ipr = posA(ip,ir)
!
!                  fact = 2d0*factB*(AOcc(ir)-AOcc(ip)) 
!
!                  !1A-4B  
!                  intXY(ipr,iac) = intXY(ipr,iac) - fact*ints(ip,ir)*Sbb(ia,ib)
!
!                  !2A-4B 
!                  intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ia)*Sab(ir,ib)
!
!                  val2 = 0
!                  do iq=1,dimOA
!                     val2 = val2 + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ia)
!                  enddo
!
!                  !3A-4B
!                  intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ir,ib)
!
!                  val2 = 0 
!                  do iq=1,dimOA
!                     val2 = val2 + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ib)
!                  enddo
!
!                  !4A-4B
!                  intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ip,ia)
!
!               enddo
!            endif
!
!         enddo
!      endif
!
!   enddo
!enddo
!
!close(iunit)
!
!! 2B terms 
!do l=1,BDimX
!
!   ia = BIndN(1,l)
!   ic = BIndN(2,l)
!   iac = posB(ia,ic)
!
!   factB = 4d0*(BOcc(ia)-BOcc(ic))
!
!   do k=1,ADimX
!      ip = AIndN(1,k)
!      ir = AIndN(2,k)
!      ipr = posA(ip,ir)
!
!      fact = 2d0*factB*(AOcc(ir)-AOcc(ip))
!
!      !1A-2B
!      intXY(ipr,iac) = intXY(ipr,iac) - fact*Emat(ip,ir)*Sbb(ia,ic)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + AOcc(iq)*Emat(iq,iq)
!      enddo
!
!      !2A-2B
!      intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ia)*Sab(ir,ic)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + (AFmat(ip,iq)-AFmat(ir,iq))*Emat(ip,iq)*Sab(iq,ia)
!      enddo
!
!      !3A-2B
!      intXY(ipr,iac) = intXY(ipr,iac) - factB*val*Sab(ir,ic)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + (AFmat(iq,ip)-AFmat(iq,ir))*Emat(iq,ir)*Sab(iq,ic)
!      enddo
!
!      !4A-2B
!      intXY(ipr,iac) = intXY(ipr,iac) - factB*val*Sab(ip,ia)
!
!   enddo 
!enddo
!
!deallocate(ints,work)
!
!print*, 'appA3_XY', norm2(intXY)
!
!end subroutine app_A3_XY

subroutine inter_A3_XX(dim1,dim2,intXX,AOcc,BOcc,RDM2Aval,RDM2Bval,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intXX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas)
double precision,intent(in) :: Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)
double precision,intent(in) :: RDM2Aval(dimOA,dimOA,dimOA,dimOA),RDM2Bval(dimOB,dimOB,dimOb,dimOB)

integer :: iunit
integer :: i,j,k,l,kl,ip,iq,ir,ia,ib,ic,iae,ipt
integer :: it,is,id,ie
double precision :: factA,factB,fact,val
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

allocate(intAB(dimOA,dimOA,dimOB,dimOB),intB(dimOB,dimOB,dimOA,NBas),&
         intA(dimOA,dimOA,dimOB,NBas),intV(dimOA,dimOA,dimOB,dimOB))
allocate(calkA(dimOA,dimOA,dimOB,NBas),calkB(dimOB,dimOB,dimOA,NBas))
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
!call dgemm('N','T',dimOA**2,dimOB**2,dimOB**2,1d0,tmp,dimOA**2,RDM2Bval,dimOB**2,0d0,intAB,dimOA**2)
call dgemm('N','N',dimOA**2,dimOB**2,dimOB**2,1d0,tmp,dimOA**2,RDM2Bval,dimOB**2,0d0,intAB,dimOA**2)
!print*, 'intAB', norm2(intAB)

! intB
allocate(tmpB(dimOB,dimOB,dimOA,NBas))
tmpB = 0
do ip=1,NBas
   do is=1,dimOA
      call dger(dimOB,dimOB,1d0,Sba(:,is),1,Sba(:,ip),1,tmpB(:,:,is,ip),dimOB)
   enddo
enddo
call dgemm('N','N',dimOB**2,dimOA*NBas,dimOB**2,1d0,RDM2Bval,dimOB**2,tmpB,dimOB**2,0d0,intB,dimOB**2)
!print*, 'intB:', norm2(intB)
deallocate(tmpB)

! intA
allocate(tmpA(dimOA,dimOA,dimOB,NBas))
tmpA = 0
do ia=1,NBas
   do id=1,dimOB
      call dger(dimOA,dimOA,1d0,Sab(:,id),1,Sab(:,ia),1,tmpA(:,:,id,ia),dimOA)
   enddo
enddo
call dgemm('N','N',dimOA**2,dimOB*NBas,dimOA**2,1d0,RDM2Aval,dimOA**2,tmpA,dimOA**2,0d0,intA,dimOA**2)
!print*, 'intA:', norm2(intA)
deallocate(tmpA)

deallocate(tmp)
allocate(tmp(dimOA,dimOA))

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
calkA = 0
calkB = 0
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

      ic = l 
      ia = k

      ! calkA
      intsr(1:dimOA,1:dimOA) = ints(1:dimOA,1:dimOA)
      call dgemv('T',dimOA**2,dimOA**2,1d0,RDM2Aval,dimOA**2,intsr,1,0d0,calkA(:,:,ic,ia),1)

      if(ia<=dimOB) then

        call dger(dimOA**2,dimOB**2,1d0,calkA(:,:,ic,ia),1,RDM2Bval(:,:,ia,ic),1,intV,dimOA**2)
         
         !do j=1,dimOB
         !   do i=1,dimOB
         !      do is=1,dimOA
         !         do ip=1,dimOA
         !           intV(ip,is,i,j) = intV(ip,is,i,j) + calkA(ip,is,l,k)*RDM2Bval(k,l,i,j)
         !         enddo
         !      enddo
         !   enddo
         !enddo

      endif

      do ie=1,dimOB

         iae = posB(ia,ie)
    
         if(iae/=0) then

            do i=1,ADimX

               ip = AIndN(1,i)
               it = AIndN(2,i)
               ipt = posA(ip,it)

               val = 0
               ! 3-1'(0)
               if(ip<=dimOA) then
                  do ir=1,dimOA
                     val = val - intAB(ir,ip,ie,ic)*ints(ir,it)
                  enddo
               endif

               ! 1-1'(0)
               do ir=1,dimOA
                  val = val + intAB(it,ir,ie,ic)*ints(ip,ir)
               enddo

               intXX(ipt,iae) = intXX(ipt,iae) + val

            enddo
         endif
     enddo

     if(k<=dimOB) then

        ic = l
        ib = k

        do ir=1,dimOA
           do ip=1,NBas
              tints(ir,ip) = ints(ip,ir)
           enddo
        enddo

        ! calkB
        call dger(dimOB**2,dimOA*NBas,1d0,RDM2Bval(:,:,ib,ic),1,tints,1,calkB,dimOB**2)

        ie = l
        ic = k

        do ia=1,NBas

           iae = posB(ia,ie)

           if(iae/=0) then

              do i=1,ADimX

                 ip = AIndN(1,i)
                 it = AIndN(2,i)
                 ipt = posA(ip,it)

                 val = 0
                 if(ia<=dimOB) then

                    ! 1-3'(0)
                    do ir=1,dimOA
                       val = val - intAB(it,ir,ic,ia)*ints(ip,ir)
                    enddo
 
                    if(ip<=dimOA) then

                       ! 3-3'(0)
                       do ir=1,dimOA
                          val = val + intAB(ir,ip,ic,ia)*ints(ir,it)
                       enddo

                    endif
                 endif 

                 intXX(ipt,iae) = intXX(ipt,iae) + val

              enddo

           endif
        enddo

     endif

  enddo
enddo

close(iunit)

!print*, 'calkA',norm2(calkA)
!print*, 'calkB',norm2(calkB)
!print*, 'intV:', norm2(intV)

! 2XX
ABmat = 0
do j=1,BDimX

   ia = BIndN(1,j)
   ie = BIndN(2,j)
   iae = posB(ia,ie)

   do i=1,ADimX

      ip = AIndN(1,i)
      it = AIndN(2,i)
      ipt = posA(ip,it)

      val = 0
      ! 2-1'(1) 
      do is=1,dimOA
         do ic=1,dimOB
            val = val + calkA(it,is,ic,ia)*intB(ie,ic,is,ip)
         enddo
      enddo

      ! 1-2'(1)
      do ir=1,dimOA
         do id=1,dimOB
            val = val + intA(it,ir,id,ia)*calkB(ie,id,ir,ip)
         enddo
      enddo

      if(ip<=dimOA) then

         ! 4-1'(1)
         do iq=1,dimOA
            do ic=1,dimOB
               val = val - calkA(iq,ip,ic,ia)*intB(ie,ic,it,iq)
            enddo
         enddo

         ! 3-2'(1)
         do ir=1,dimOA
            do id=1,dimOB
               val = val - intA(ir,ip,id,ia)*calkB(ie,id,it,ir)
            enddo
         enddo

         ! 4-2'(2)
         do id=1,dimOB
            do iq=1,dimOA
               val = val - intV(iq,ip,ie,id)*Sab(iq,id)*Sab(it,ia)
            enddo
         enddo

      endif

      ! 2-2'(2)
      do id=1,dimOB
         do is=1,dimOA
            val = val + intV(it,is,ie,id)*Sab(is,ia)*Sab(ip,id)
         enddo
      enddo

      if(ia<=dimOB) then

         ! 2-3'(1)
         do is=1,dimOA
            do ic=1,dimOB
               val = val - calkA(it,is,ie,ic)*intB(ic,ia,is,ip)
            enddo
         enddo

         ! 1-4'(1)
         do ir=1,dimOA
            do ib=1,dimOB
               val = val - intA(it,ir,ie,ib)*calkB(ib,ia,ir,ip)
            enddo  
         enddo

!         ! HERE
!         intXX(ipt,iae) = intXX(ipt,iae) + val
!
         !val = 0 
         ! 2-4'(2)
         do ib=1,dimOB
            do is=1,dimOA
               val = val - intV(it,is,ib,ia)*Sab(is,ib)*Sab(ip,ie)
            enddo
         enddo

         !intXX(ipt,iae) = intXX(ipt,iae) + val*Sab(ip,ie)

         !val = 0
         if(ip<=dimOA) then

            ! 4-3'(1)
            do iq=1,dimOA
               do ic=1,dimOB
                  val = val + calkA(iq,ip,ie,ic)*intB(ic,ia,it,iq)
               enddo
            enddo

            ! 3-4'(1)
            do ir=1,dimOA
               do ib=1,dimOB
                  val = val + intA(ir,ip,ie,ib)*calkB(ib,ia,it,ir)
               enddo  
            enddo

            ! 4-4'(2)
            do ib=1,dimOB
               do iq=1,dimOA
                  val = val + intV(iq,ip,ib,ia)*Sab(iq,ie)*Sab(it,ib)
               enddo
            enddo

         endif
      endif

      intXX(ipt,iae) = intXX(ipt,iae) + val 

   enddo
enddo

deallocate(intsr,tints,ints,work,tmp)
deallocate(calkB,calkA)
deallocate(intV,intA,intB,intAB)

intXX = -2d0*intXX

!print*, 'intXX-A3',norm2(intXX)
!print*, ''

end subroutine inter_A3_XX

!subroutine app_A3_XX(dim1,dim2,intXX,AOcc,BOcc,AFmat,BFmat,Sab, &
!                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
!                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
!character(*) :: IntKFile
!double precision,intent(inout) :: intXX(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
!                               AFmat(NBas,NBas),BFmat(NBas,NBas),&
!                               Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)
!
!integer :: iunit
!integer :: i,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
!double precision :: factA,factB,fact,val,val3B,val4B,valTr
!double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas),&
!                    Amat(NBas,NBas),Emat(NBas,dimOA)
!double precision,allocatable :: work(:),ints(:,:)
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(work(NBas*NBas),ints(NBas,NBas))
!
!Amat = 0
!do ib=1,dimOB
!   do ir=1,NBas !dimOA
!      do ip=1,NBas
!         Amat(ip,ir) = Amat(ip,ir) + BOcc(ib)*Sab(ip,ib)*Sab(ir,ib)
!      enddo
!   enddo
!enddo
!
!valTr = 0
!do iq=1,dimOA
!   valTr = valTr + AOcc(iq)*Amat(iq,iq)
!enddo
!
!Sbb = 0
!do ic=1,NBas
!   do ia=1,NBas
!      do iq=1,dimOA
!         Sbb(ia,ic) = Sbb(ia,ic) + AOcc(iq)*Sab(iq,ia)*Sab(iq,ic)
!      enddo
!   enddo
!enddo
!
!!(FO|FO): (AA|BB)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOA)
!
!! one loop over integrals
!ints = 0
!Emat = 0
!kl = 0
!do l=1,dimOB
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOA)
!
!      call ints_modify(NBas,dimOA,ints,NBas,work,&
!                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)
!
!      ic = l
!      ia = k
!
!      iac = posB(ia,ic)
!
!      if(iac/=0) then
!
!         factB = 4d0*(BOcc(ic)-BOcc(ia))
!
!         do i=1,ADimX
!            ip = AIndN(1,i)
!            ir = AIndN(2,i)
!            ipr = posA(ip,ir)
!            
!            fact  = 2d0*(AOcc(ir)-AOcc(ip))*factB
!
!            ! 1A-1B
!            intXX(ipr,iac) = intXX(ipr,iac) - fact*valTr*ints(ip,ir)
!
!            val = 0
!            do iq=1,dimOA
!               val = val + AOcc(iq)*ints(iq,iq)
!            enddo
!
!            !2A-1B
!            intXX(ipr,iac) = intXX(ipr,iac) - fact*val*Amat(ip,ir)
!
!            val = 0
!            do iq=1,dimOA
!               val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Amat(iq,ir)
!            enddo
!
!            !3A-1B
!            intXX(ipr,iac) = intXX(ipr,iac) - factB*val
!
!            val = 0
!            do iq=1,dimOA
!               val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Amat(ip,iq)
!            enddo
!
!            !4A-1B
!            intXX(ipr,iac) = intXX(ipr,iac) - factB*val
!
!         enddo
!
!      endif
!
!      ! all 2B terms 
!      if(k==l) then
!
!         ib = k
!         Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + BOcc(ib)*ints(1:NBas,1:dimOA)
!
!      endif
!
!      ib = l
!      ia = k
!   
!      val3B = 0 
!      do iq=1,dimOA
!         val3B = val3B + AOcc(iq)*ints(iq,iq)
!      enddo
! 
!      do ic=1,dimOB
!
!         iac = posB(ia,ic)
!
!         if(iac/=0) then
!
!            factB = 2d0*(BFmat(ia,ib)-BFmat(ic,ib))
!
!            do i=1,ADimX
!
!               ip = AIndN(1,i)
!               ir = AIndN(2,i)
!               ipr = posA(ip,ir)
!
!               fact = 2d0*factB*(AOcc(ir)-AOcc(ip))
!
!               !1A-3B
!               intXX(ipr,iac) = intXX(ipr,iac) - fact*Sbb(ic,ib)*ints(ip,ir)
!               
!               !2A-3B
!               intXX(ipr,iac) = intXX(ipr,iac) - fact*val3B*Sab(ip,ic)*Sab(ir,ib)
!
!               val = 0
!               do iq=1,dimOA
!                  val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ic)
!               enddo
!
!               !3A-3B
!               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ib)
!
!               val = 0
!               do iq=1,dimOA
!                  val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ib)
!               enddo
!
!               !4A-3B
!               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ic)
!
!            enddo
!         endif
!      enddo
!
!      if(k<=dimOB) then
!
!         ic = l
!         ib = k
!
!         val4B = 0
!         do iq=1,dimOA
!            val4B = val4B + AOcc(iq)*ints(iq,iq)
!         enddo
!
!         do ia=1,NBas
!
!            iac = posB(ia,ic)
!
!            if(iac/=0) then
!
!               factB = 2d0*(BFmat(ib,ia)-BFmat(ib,ic))
!
!               do i=1,ADimX
!
!                  ip = AIndN(1,i)
!                  ir = AIndN(2,i)
!                  ipr = posA(ip,ir)
!
!                  fact = 2d0*factB*(AOcc(ir)-AOcc(ip))
!
!                  !1A-4B
!                  intXX(ipr,iac) = intXX(ipr,iac) - fact*ints(ip,ir)*Sbb(ib,ia)
!
!                  !2A-4B
!                  intXX(ipr,iac) = intXX(ipr,iac) - fact*val4B*Sab(ip,ib)*Sab(ir,ia)
!
!                  val = 0
!                  do iq=1,dimOA
!                     val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ib)
!                  enddo
!
!                  !3A-4B
!                  intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ia)
!
!                  val = 0
!                  do iq=1,dimOA
!                     val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ia)
!                  enddo
!
!                  !4A-4B
!                  intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ib)
!
!               enddo
!            endif
!
!         enddo
!      endif
!
!   enddo
!enddo
!
!close(iunit)
!
!! 2B terms
!do l=1,BDimX
!
!   ia = BIndN(1,l)
!   ic = BIndN(2,l)
!   iac = posB(ia,ic)
!
!   factB = 4d0*(BOcc(ic)-BOcc(ia))
!
!   do k=1,ADimX
!
!      ip = AIndN(1,k)
!      ir = AIndN(2,k)
!      ipr = posA(ip,ir)
!
!      fact = 2d0*(AOcc(ir)-AOcc(ip))*factB
!
!      !1A-2B
!      intXX(ipr,iac) = intXX(ipr,iac) - fact*Emat(ip,ir)*Sbb(ic,ia)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + AOcc(iq)*Emat(iq,iq)
!      enddo
!
!      !2A-2B
!      intXX(ipr,iac) = intXX(ipr,iac) - fact*val*Sab(ip,ic)*Sab(ir,ia)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + (AFmat(ip,iq)-AFmat(ir,iq))*Emat(ip,iq)*Sab(iq,ic)
!      enddo
!
!      !3A-2B
!      intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ia)
!
!      val = 0
!      do iq=1,dimOA
!         val = val + (AFmat(iq,ip)-AFmat(iq,ir))*Emat(iq,ir)*Sab(iq,ia)
!      enddo
!
!      !4A-2B
!      intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ic)
!
!   enddo 
!enddo
!
!deallocate(ints,work)
!
!print*, 'appA3_XX',norm2(intXX)
!
!end subroutine app_A3_XX
!
!
!subroutine app_A2_YX(dim1,dim2,intYX,BOcc,Fmat,Smat,  &
!                  AXELE,Vabb,BXELE,Vbab,IntKFile,IntJFile,&
!                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),&
!                      posA(NBas,NBas),posB(NBas,NBas)
!logical,intent(in) :: trans
!character(*) :: IntKFile,IntJFile
!double precision,intent(inout) :: intYX(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
!                               BOcc(NBas),Fmat(NBas,NBas), & 
!                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)
!
!integer :: iunit,i,j,k,l,kl
!integer :: ip,iq,ir,is,it,iu,itu,ipr
!integer :: iw,ipw
!double precision :: fact,factA,factF,val
!double precision :: Sxx(NBas,NBas)
!double precision,allocatable :: work(:),ints(:,:),tmp(:,:)
!double precision,external  :: FRDM2
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(tmp(dimOB,dimOB),work(NBas*NBas),ints(NBas,NBas))
!
!!(FO|FO):(BB|AB) or (AA|BA)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
!! one loop over integrals
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)
!
!      iq = l
!      it = k
!
!      ! terms 1 and 3
!      do i=1,BDimX
! 
!         ip = BIndN(1,i)
!         ir = BIndN(2,i)
!         ipr = posB(ip,ir)
!
!         fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)*ints(ip,ir)
!         factF = (Fmat(ip,iq)-Fmat(ir,iq))*ints(ip,iq)
!
!         do iu=1,dimOA
! 
!            itu = posA(it,iu)
!          
!            if(itu/=0) then
!
!               factA = 2d0*(AOcc(it)-AOcc(iu))
!
!               if(trans) then 
!                  intYX(ipr,itu) = intYX(ipr,itu) - factA*fact*Smat(iu,iq) &
!                                                  - factA*factF*Smat(iu,ir)
!               else
!                  intYX(itu,ipr) = intYX(itu,ipr) - factA*fact*Smat(iu,iq) &
!                                                  - factA*factF*Smat(iu,ir)
!               endif
!
!            endif
!         enddo
!      enddo
!
!   enddo
!enddo
!
!close(iunit)
!
!! 2nd term
!!(FF|OO): (AB|BB) or (BA|AA)
!open(newunit=iunit,file=trim(IntJFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*NBas)
!
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,dimOB
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*NBas)
!
!      call ints_modify(NBas,NBas,ints,NBas,work,&
!                       Sxx(k,l)/BXELE,Vbab,Vabb(k,l)/AXELE,Smat)
!
!      iq = l
!      ir = k
!
!      ! term 4
!      if(trans) then
!         do i=1,ADimX
!
!            it = AIndN(1,i)
!            iu = AIndN(2,i)
!            itu = posA(it,iu)
!
!            factA = 2d0*(AOcc(it)-AOcc(iu))
!
!            do ip=1,NBas
!
!               ipr = posB(ip,ir)
!
!               if(ipr/=0) then
!
!                  factF = Fmat(iq,ip)-Fmat(iq,ir)
!                  intYX(ipr,itu) = intYX(ipr,itu) - factA*factF*Smat(iu,iq)*ints(it,ip)
!
!               endif
!
!            enddo
!         enddo
!
!      else
!
!         do ip=1,NBas
!
!            ipr = posB(ip,ir)
!
!            if(ipr/=0) then
!
!               factF = Fmat(iq,ip)-Fmat(iq,ir)
!
!               do i=1,ADimX
!
!                  it = AIndN(1,i)
!                  iu = AIndN(2,i)
!                  itu = posA(it,iu)
!
!                  factA = 2d0*(AOcc(it)-AOcc(iu))
!                  intYX(itu,ipr) = intYX(itu,ipr) - factA*factF*Smat(iu,iq)*ints(it,ip)
!
!               enddo
!            endif
!         enddo
!      endif 
!
!     ! term 2
!     if(k==l) then
!
!        iq = l
!
!        if(trans) then
!
!           do j=1,ADimX
!              it = AIndN(1,j)
!              iu = AIndN(2,j)
!              itu = posA(it,iu)
!
!              factA = 2d0*(AOcc(it)-AOcc(iu))
!
!              do i=1,BDimX
!                 ip = BIndN(1,i)
!                 ir = BIndN(2,i)
!                 ipr = posB(ip,ir)
!
!                 fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)
!
!                 intYX(ipr,itu) = intYX(ipr,itu) - fact*factA*Smat(iu,ir)*ints(it,ip)
!
!              enddo
!           enddo   
!
!        else
!           do i=1,BDimX
!              ip = BIndN(1,i)
!              ir = BIndN(2,i)
!              ipr = posB(ip,ir)
!
!              fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)
!
!              do j=1,ADimX
!                 it = AIndN(1,j)
!                 iu = AIndN(2,j)
!                 itu = posA(it,iu)
!
!                 factA = 2d0*(AOcc(it)-AOcc(iu))
!
!                 intYX(itu,ipr) = intYX(itu,ipr) - fact*factA*Smat(iu,ir)*ints(it,ip)
!
!              enddo
!           enddo
!        endif
!
!     endif 
!
!   enddo
!enddo
!
!close(iunit)
!
!deallocate(ints,work,tmp)
!
!print*,'appYX', norm2(intYX)
!
!end subroutine app_A2_YX 
!
!subroutine app_A2_XY(dim1,dim2,intXY,BOcc,Fmat,Smat,  &
!                  AXELE,Vabb,BXELE,Vbab,IntKFile,&
!                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
!logical,intent(in) :: trans
!character(*) :: IntKFile
!double precision,intent(inout) :: intXY(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
!                               BOcc(NBas),Fmat(NBas,NBas),&
!                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)
!
!integer :: iunit
!integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
!double precision :: factA,factF,fact,val
!double precision :: Sxx(NBas,NBas)
!double precision,allocatable :: tmp(:,:),work(:),ints(:,:)
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(tmp(dimOA,dimOB),work(NBas*NBas),ints(NBas,NBas))
!
!!(FO|FO):(BB|AB) or (AA|BA)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
!! one loop over integrals
!tmp  = 0
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,dimOA
!      kl = kl + 1
!      read(iunit,rec=k+(l-1)*NBas) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)
!
!      iq = l
!      iu = k
!
!      ! terms 1 and 3  
!      do i=1,BDimX
!          
!         ip = BIndN(1,i)
!         ir = BIndN(2,i)
!         ipr = posB(ip,ir)
!         
!         fact = 2d0*(BOcc(ip)-BOcc(ir))*BOcc(iq)*ints(ip,ir)
!         factF = (Fmat(ir,iq)-Fmat(ip,iq))*ints(ir,iq)
!
!         do it=1,NBas
!            
!            itu = posA(it,iu)
!
!            if(itu/=0) then
!             
!               factA = 2d0*(AOcc(iu)-AOcc(it))
!
!               if(trans) then
!                  intXY(ipr,itu) = intXY(ipr,itu) - factA*fact*Smat(it,iq) &
!                                                  - factA*factF*Smat(it,ip)
!               else                               
!                  intXY(itu,ipr) = intXY(itu,ipr) - factA*fact*Smat(it,iq) &
!                                                  - factA*factF*Smat(it,ip)
!               endif
!
!            endif
!         enddo
!      enddo
!
!      ! term 2
!      ir = l
!      iu = k
!     
!      do iq=1,dimOB
!         tmp(iu,ir) = tmp(iu,ir) + BOcc(iq)*ints(iq,iq)
!      enddo
!
!      ! term4
!      do it=1,NBas
!
!         itu = posA(it,iu) 
!
!         if(itu/=0) then
!
!            factA = 2d0*(AOcc(iu)-AOcc(it))
! 
!            do ip=1,NBas
!
!               ipr = posB(ip,ir)
! 
!               if(ipr/=0) then
!
!                  if(trans) then            
!                     do iq=1,dimOB
!                        intXY(ipr,itu) = intXY(ipr,itu) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
!                                       * Smat(it,iq)*ints(ip,iq)
!                     enddo
!                  else
!                     do iq=1,dimOB
!                        intXY(itu,ipr) = intXY(itu,ipr) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
!                                       * Smat(it,iq)*ints(ip,iq)
!                     enddo
!                  endif
!
!               endif
!            enddo 
!         endif
!      enddo
!
!   enddo
!enddo
!
!close(iunit)
!
!if(trans) then
!
!   do l=1,ADimX
!      
!      it = AIndN(1,l)
!      iu = AIndN(2,l)
!      itu = posA(it,iu)
!
!      factA = 2d0*(AOcc(iu)-AOcc(it))
! 
!      do k=1,BDimX
!       
!         ip = BIndN(1,k)
!         ir = BIndN(2,k)
!         ipr = posB(ip,ir)
!
!         fact = 2d0*(BOcc(ip)-BOcc(ir))
!
!         intXY(ipr,itu) = intXY(ipr,itu) - fact*factA*tmp(iu,ir)*Smat(it,ip)
!
!      enddo
!   enddo
!
!else
!
!   do l=1,BDimX
!      
!      ip = BIndN(1,l)
!      ir = BIndN(2,l)
!      ipr = posB(ip,ir)
!
!      fact = 2d0*(BOcc(ip)-BOcc(ir)) 
!
!      do k=1,ADimX
!   
!         it = AIndN(1,k)
!         iu = AIndN(2,k)
!         itu = posA(it,iu)
!
!         factA = 2d0*(AOcc(iu)-AOcc(it))
!         intXY(itu,ipr) = intXY(itu,ipr) - fact*factA*tmp(iu,ir)*Smat(it,ip)
!   
!      enddo
!   enddo
!
!endif
!
!deallocate(ints,work,tmp)
!
!print*, 'appXY',norm2(intXY)
!
!end subroutine app_A2_XY
!
!subroutine app_A2_YY(dim1,dim2,intYY,BOcc,Fmat,Smat,  &
!                  AXELE,Vabb,BXELE,Vbab,IntKFile,&
!                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
!logical,intent(in) :: trans
!character(*) :: IntKFile
!double precision,intent(inout) :: intYY(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
!                               BOcc(NBas),Fmat(NBas,NBas),&
!                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)
!
!integer :: iunit
!integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
!double precision :: factA,factF,fact,val
!double precision :: Sxx(NBas,NBas)
!double precision,allocatable :: tmp(:,:),work(:),ints(:,:)
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(tmp(NBas,dimOB),work(NBas*NBas),ints(NBas,NBas))
!
!!(FO|FO):(BB|AB) or (AA|BA)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
!! one loop over integrals
!tmp  = 0
!ints = 0
!kl = 0
!do l=1,dimOB
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)
!
!      iq = l
!      it = k
!
!      ! term 1 and 3
!      do i=1,BDimX
!
!         ip = BIndN(1,i) 
!         ir = BIndN(2,i) 
!         ipr = posB(ip,ir)
!
!         fact = 2d0*BOcc(iq)*(BOcc(ip)-BOcc(ir))*ints(ip,ir)
!         factF = (Fmat(ir,iq)-Fmat(ip,iq))*ints(ir,iq)
!
!         do iu=1,dimOA
!
!            itu = posA(it,iu)
!
!            if(itu/=0) then
!
!              factA = 2d0*(AOcc(it)-AOcc(iu))
!              if(trans) then 
!                 intYY(ipr,itu) = intYY(ipr,itu) - fact*factA*Smat(iu,iq) &
!                                                 - factF*factA*Smat(iu,ip)
!              else
!                 intYY(itu,ipr) = intYY(itu,ipr) - fact*factA*Smat(iu,iq) &
!                                                 - factF*factA*Smat(iu,ip)
!              endif
!              
!            endif
!         enddo
!      enddo
!
!      ir = l 
!      it = k
!
!      ! term 2
!      do iq=1,dimOB  
!         tmp(it,ir) = tmp(it,ir) + BOcc(iq)*ints(iq,iq)
!      enddo
!
!      ! term 4
!      do iu=1,dimOA
!         
!         itu = posA(it,iu)
!
!         if(itu/=0) then
!            
!            factA = 2d0*(AOcc(it)-AOcc(iu))
!
!            do ip=1,NBas
!
!               ipr=posB(ip,ir)
!
!               if(ipr/=0) then
!
!                  if(trans) then
!
!                     do iq=1,dimOB
!                        intYY(ipr,itu) = intYY(ipr,itu) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
!                                       * Smat(iu,iq)*ints(ip,iq) 
!                     enddo
! 
!                  else
!
!                     do iq=1,dimOB
!                        intYY(itu,ipr) = intYY(itu,ipr) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
!                                       * Smat(iu,iq)*ints(ip,iq) 
!                     enddo
!
!                  endif
!               endif
!            enddo
!
!         endif
!      enddo
!
!   enddo
!enddo
!
!close(iunit)
!
!if(trans) then
!
!   do l=1,ADimX
!
!      it = AIndN(1,l)
!      iu = AIndN(2,l)
!      itu = posA(it,iu)
!
!      factA = 2d0*(AOcc(it)-AOcc(iu))
!
!      do k=1,BDimX
!         ip = BIndN(1,k)
!         ir = BIndN(2,k)
!         ipr = posB(ip,ir)
!
!         fact = 2d0*(BOcc(ip)-BOcc(ir))
!         intYY(ipr,itu) = intYY(ipr,itu) - fact*factA*tmp(it,ir)*Smat(iu,ip)
!
!      enddo
!   enddo
!
!else
!
!   do l=1,BDimX
!
!      ip = BIndN(1,l)
!      ir = BIndN(2,l)
!      ipr = posB(ip,ir)
!   
!      fact = 2d0*(BOcc(ip)-BOcc(ir))
!   
!      do k=1,ADimX
!
!         it = AIndN(1,k) 
!         iu = AIndN(2,k) 
!         itu = posA(it,iu)
!         factA = 2d0*(AOcc(it)-AOcc(iu))
!   
!         intYY(itu,ipr) = intYY(itu,ipr) - fact*factA*tmp(it,ir)*Smat(iu,ip)
!   
!      enddo
!   enddo
!
!endif
!
!deallocate(ints,work,tmp)
!
!print*, 'appYY',norm2(intYY)
!
!end subroutine app_A2_YY
!
!subroutine app_A2_XX(dim1,dim2,intXX,BOcc,Fmat,Smat,  &
!                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
!                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
!implicit none
!
!integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
!integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
!logical,intent(in) :: trans
!character(*) :: IntKFile
!double precision,intent(inout) :: intXX(dim1,dim2)
!double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
!                               BOcc(NBas),Fmat(NBas,NBas),&
!                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)
!
!integer :: iunit
!integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
!double precision :: factA,factF,fact,val
!double precision :: Sxx(NBas,NBas)
!double precision,allocatable :: tmp(:,:),work(:),ints(:,:)
!
!Sxx = 0
!do i=1,NBas
!   Sxx(i,i) = 1d0
!enddo
!
!allocate(tmp(NBas,dimOA),work(NBas*NBas),ints(NBas,NBas))
!
!!do i=1,NBas
!!    print*, i,posA(i,i),posB(i,i)
!!enddo
!
!!(FO|FO):(BB|BA) or (AA|AB)
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBas*dimOB)
!
!! one loop over integrals
!ints = 0
!tmp  = 0
!kl = 0
!do l=1,dimOA
!   do k=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                       Smat(l,k)/AXELE,Vabb,Vbab(l,k)/BXELE,Sxx)
!
!      if(k<=dimOB) then
!         iu = l
!         iq = k
!
!         ! terms: 1, 3
!         do i=1,BDimX
!            ip = BIndN(1,i)
!            ir = BIndN(2,i)
!            ipr = posB(ip,ir)
!
!            fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)*ints(ip,ir)
!            factF = (Fmat(ip,iq)-Fmat(ir,iq))*ints(ip,iq)
!
! 
!            do it=1,NBas
!
!               itu = posA(it,iu)
!
!               if(itu/=0) then
!
!                  factA = 2d0*(AOcc(iu)-AOcc(it))
!
!                  if(trans) then
!                     intXX(ipr,itu) = intXX(ipr,itu) - factA*fact*Smat(it,iq) &
!                                                     - factA*factF*Smat(it,ir)
!                  else
!                     intXX(itu,ipr) = intXX(itu,ipr) - factA*fact*Smat(it,iq) &
!                                                     - factA*factF*Smat(it,ir)
!                  endif
!
!               endif
!            enddo
!
!         enddo
!      endif
!   
!      iu = l
!      ip = k
!
!      ! term 2
!      do iq=1,dimOB
!         tmp(ip,iu) = tmp(ip,iu) + BOcc(iq)*ints(iq,iq) 
!      enddo
!
!      ! term 4
!      do ir=1,dimOB
!
!         ipr = posB(ip,ir)
!         
!         if(ipr/=0) then
!
!            do it=1,NBas
!               
!               itu = posA(it,iu)
!
!               if(itu/=0) then
!
!                  factA = 2d0*(AOcc(iu)-AOcc(it))
!
!                  val = 0 
!                  do iq=1,dimOB
!                     val = val - (Fmat(iq,ip)-Fmat(iq,ir)) * &
!                           ints(iq,ir)*Smat(it,iq)
!                  enddo
!
!                  if(trans) then
!                     intXX(ipr,itu) = intXX(ipr,itu) + factA*val
!                  else
!                     intXX(itu,ipr) = intXX(itu,ipr) + factA*val
!                  endif 
!  
!               endif               
!            enddo
!         endif
!      enddo       
!
!   enddo
!enddo
!
!close(iunit)
!
!if(trans) then
!
!   do l=1,ADimX
!
!      it = AIndN(1,l)
!      iu = AIndN(2,l)
!      itu = posA(it,iu)
!
!      factA = 4d0*(AOcc(iu)-AOcc(it))
!
!      do k=1,BDimX
!         ip = BIndN(1,k)
!         ir = BIndN(2,k)
!         ipr = posB(ip,ir)
!
!         intXX(ipr,itu) = intXX(ipr,itu) - factA*(BOcc(ir)-BOcc(ip))*tmp(ip,iu)*Smat(it,ir)
!
!      enddo
!   enddo
!
!else
!   do l=1,BDimX
!   
!      ip = BIndN(1,l)
!      ir = BIndN(2,l)
!      ipr = posB(ip,ir)
!   
!      fact = 4d0*(BOcc(ir)-BOcc(ip))
!   
!      do k=1,ADimX
!   
!         it = AIndN(1,k)
!         iu = AIndN(2,k)
!         itu = posA(it,iu)
!   
!         factA = fact*(AOcc(iu)-AOcc(it))
!   
!         intXX(itu,ipr) = intXX(itu,ipr) - factA*tmp(ip,iu)*Smat(it,ir)
!   
!      enddo
!   enddo
!endif
!
!print*,'appXX:', norm2(intXX)
!
!deallocate(ints,work,tmp)
!
!end subroutine app_A2_XX

subroutine inter_A2_YX(dim1,dim2,intYX,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,IntJFile,&
                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
!                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans,B)
implicit none

!type(SystemBlock) :: B
integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),&
                      posA(NBas,NBas),posB(NBas,NBas)
logical,intent(in) :: trans
character(*) :: IntKFile,IntJFile
double precision,intent(inout) :: intYX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit,i,j,k,l,kl
integer :: ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val,val1,val2 
double precision :: Sxx(NBas,NBas),Amat(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:),tmp(:,:,:,:)
double precision,external  :: FRDM2

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOB,dimOB,dimOB,NBas),work(NBas*NBas),ints(NBas,NBas))

call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp,dimOB**3)

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

      ! terms 1 and 3
      do iu=1,dimOA

         itu = posA(it,iu)
 
         if(itu/=0) then

            fact = 2d0*(AOcc(it)-AOcc(iu))

            do i=1,BDimX
               ip = BIndN(1,i)
               iw = BindN(2,i)
               ipw = posB(ip,iw)

               val = 0
               do ir=1,dimOB
                  val = val - tmp(iw,ir,iq,iu)*ints(ip,ir)
               enddo

               if(ip<=dimOB) then
                  do ir=1,dimOB
                     val = val + tmp(ir,ip,iq,iu)*ints(ir,iw)
                  enddo
               endif 

               if(trans) then
                  intYX(ipw,itu) = intYX(ipw,itu) + fact*val
               else
                  intYX(itu,ipw) = intYX(itu,ipw) + fact*val
               endif

            enddo
         endif       
      enddo

      ! term 4
      do ip=1,dimOB
         do is=1,dimOB
            do ir=1,dimOB
               Amat(ip,it) = Amat(ip,it) + RDM2val(ir,is,ip,iq)*ints(ir,is)
            enddo 
         enddo 
      enddo 

   enddo
enddo

! 4th term

do l=1,BDimX

   ip = BIndN(1,l)
   iw = BIndN(2,l)
   ipw = posB(ip,iw)

   do k=1,ADimX

      it = AIndN(1,k)
      iu = AIndN(2,k)
      itu = posA(it,iu)

      fact = 2d0*(AOcc(it)-AOcc(iu))

      if(trans) then
         intYX(ipw,itu) = intYX(ipw,itu) + fact*Amat(ip,it)*Smat(iu,iw)
      else
         intYX(itu,ipw) = intYX(itu,ipw) + fact*Amat(ip,it)*Smat(iu,iw)
      endif

   enddo
enddo

close(iunit)

deallocate(tmp)
allocate(tmp(NBas,dimOB,dimOB,dimOB))

!tmp=0
!! test for o3v matrix
!   do iq=1,dimOB
!      do ir=1,dimOB
!         do iw=1,dimOB
!            do iu=1,NBas 
!               do is=1,dimOB
!               tmp(iu,iw,ir,iq) = tmp(iu,iw,ir,iq) &
!                                + RDM2val(is,iw,ir,iq)*Smat(iu,is)
!                               ! +  FRDM2(is,iw,ir,iq,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)*Smat(iu,is)
!            enddo
!         enddo
!      enddo
!   enddo
!enddo
!print*, ''
!print*, 'tmp1', norm2(tmp)
!
!tmp=0
! sth wrong here!?
call dgemm('N','N',NBas,dimOB**3,dimOB,1d0,Smat,NBas,RDM2val,dimOB,0d0,tmp,NBas)
!print*, 'tmp2', norm2(tmp)

! 2nd term
!(FF|OO): (AB|BB) or (BA|AA)
open(newunit=iunit,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBas*NBas)

ints = 0
kl = 0
do l=1,dimOB
   do k=1,dimOB
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*NBas)

      call ints_modify(NBas,NBas,ints,NBas,work,&
                       Sxx(k,l)/BXELE,Vbab,Vabb(k,l)/AXELE,Smat)

      iq = l
      ir = k

      if(trans) then

         do j=1,ADimX

            it = AIndN(1,j)
            iu = AIndN(2,j)
            itu = posA(it,iu)

            fact = 2d0*(AOcc(it)-AOcc(iu))

            do i=1,BDimX
   
               ip = BIndN(1,i)
               iw = BIndN(2,i)
               ipw = posB(ip,iw)

               intYX(ipw,itu) = intYX(ipw,itu) - fact*tmp(iu,iw,ir,iq)*ints(it,ip)
      
            enddo
         enddo

      else

         do i=1,BDimX

            ip = BIndN(1,i)
            iw = BIndN(2,i)
            ipw = posB(ip,iw)

            do j=1,ADimX

               it = AIndN(1,j)
               iu = AIndN(2,j)
               itu = posA(it,iu)

               fact = 2d0*(AOcc(it)-AOcc(iu))

               intYX(itu,ipw) = intYX(itu,ipw) - fact*tmp(iu,iw,ir,iq)*ints(it,ip)

            enddo
         enddo

      endif

   enddo
enddo

close(iunit)

!print*, 'intYX',norm2(intYX)

deallocate(ints,work,tmp)

end subroutine inter_A2_YX

subroutine inter_A2_XY(dim1,dim2,intXY,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),&
                      posA(NBas,NBas),posB(NBas,NBas)
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intXY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit,i,k,l,kl
integer :: ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val
double precision :: Sxx(NBas,NBas),Amat(dimOB,dimOA)
double precision,allocatable :: work(:),ints(:,:),tmp(:,:,:,:)
double precision,external  :: FRDM2

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOB,dimOB,dimOB,NBas),work(NBas*NBas),ints(NBas,NBas))

call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp,dimOB**3)

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
Amat = 0
ints = 0
do l=1,dimOB
   do k=1,dimOA
      read(iunit,rec=k+(l-1)*NBas) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)

      iq = l
      iu = k

      ! terms 1 and 3
      do it=1,NBas

         itu = posA(it,iu)
 
         if(itu/=0) then

            fact = 2d0*(AOcc(iu)-AOcc(it))
 
            do i=1,BDimX

               ip = BIndN(1,i)
               iw = BIndN(2,i)
               ipw = posB(ip,iw)
 
               val = 0
               do ir=1,dimOB
                  val = val + tmp(ir,iw,iq,it)*ints(ip,ir)
               enddo

               if(ip<=dimOB) then
                  do ir=1,dimOB
                     val = val - tmp(ip,ir,iq,it)*ints(ir,iw)
                  enddo 
               endif                

               if(trans) then
                  intXY(ipw,itu) = intXY(ipw,itu) + fact*val
               else
                  intXY(itu,ipw) = intXY(itu,ipw) + fact*val
               endif

            enddo
         endif
      enddo 

      ! term 4
      do iw=1,dimOB
         do is=1,dimOB
            do ir=1,dimOB
               Amat(iw,iu) = Amat(iw,iu) + RDM2val(ir,is,iw,iq)*ints(ir,is)
            enddo
         enddo
      enddo

     ! term 2
     iw = l
     iu = k

     do it=1,NBas

        itu = posA(it,iu)
 
        if(itu/=0) then

           fact = 2d0*(AOcc(iu)-AOcc(it))

           do ip=1,dimOB

              ipw = posB(ip,iw)

              if(ipw/=0) then
                 
                 val = 0
                 do iq=1,dimOB
                    do ir=1,dimOB
                       val = val - tmp(iq,ir,ip,it)*ints(ir,iq)
                    enddo  
                 enddo

                 if(trans) then
                    intXY(ipw,itu) = intXY(ipw,itu) + fact*val
                 else
                    intXY(itu,ipw) = intXY(itu,ipw) + fact*val
                 endif

              endif

           enddo
        endif
     enddo

   enddo
enddo

! 4th term
do l=1,BDimX
   ip = BIndN(1,l)
   iw = BIndN(2,l)
   ipw = posB(ip,iw)

   do k=1,ADimX
      it = AIndN(1,k)
      iu = AIndN(2,k)
      itu = posA(it,iu)

      fact = 2d0*(AOcc(iu)-AOcc(it))

      if(trans) then
         intXY(ipw,itu) = intXY(ipw,itu) + fact*Amat(iw,iu)*Smat(it,ip)
      else
         intXY(itu,ipw) = intXY(itu,ipw) + fact*Amat(iw,iu)*Smat(it,ip)
      endif

   enddo

enddo

!print*, 'intXY:',norm2(intXY)

close(iunit)

deallocate(ints,work,tmp)

! terms 2 and 4 tested

end subroutine inter_A2_XY

subroutine inter_A2_YY(dim1,dim2,intYY,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),&
                      posA(NBas,NBas),posB(NBas,NBas)
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intYY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit,i,k,l,kl
integer :: ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val
double precision :: Sxx(NBas,NBas),Amat(dimOB,NBas)
double precision,allocatable :: work(:),ints(:,:),tmp(:,:,:,:)
double precision,external  :: FRDM2

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOB,dimOB,dimOB,NBas),work(NBas*NBas),ints(NBas,NBas))

tmp = 0
call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp,dimOB**3)

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

      ! terms 1 and 3
      do iu=1,dimOA

         itu = posA(it,iu)

         if(itu/=0) then 
          
            fact = 2d0*(AOcc(it)-AOcc(iu))

            do i=1,BDimX
               ip = BIndN(1,i)
               iw = BIndN(2,i)
               ipw = posB(ip,iw)

               val = 0
               do ir=1,dimOB
                 val = val + tmp(ir,iw,iq,iu)*ints(ip,ir)
               enddo
              
               if(ip<=dimOB) then
                  do ir=1,dimOB
                     val = val - tmp(ip,ir,iq,iu)*ints(ir,iw)
                  enddo
               endif

               if(trans) then
                  intYY(ipw,itu) = intYY(ipw,itu) + fact*val
               else
                  intYY(itu,ipw) = intYY(itu,ipw) + fact*val
               endif

            enddo
         endif
      enddo

      !! 4th term
      do iw=1,dimOB
         do is=1,dimOB
            do ir=1,dimOB
               Amat(iw,it) = Amat(iw,it) + RDM2val(ir,is,iw,iq)*ints(ir,is)
            enddo 
         enddo 
      enddo

      ! 2nd term
      iw = l
      it = k

      do iu=1,dimOB !NBas 

         itu=posA(it,iu)

         if(itu/=0) then

            fact = 2d0*(AOcc(it)-AOcc(iu))

            do ip=1,dimOB

               ipw = posB(ip,iw)
               
               if(ipw/=0) then

                 val = 0 
                 do iq=1,dimOB
                    do ir=1,dimOB
                       val = val + tmp(iq,ir,ip,iu)*ints(ir,iq)
                    enddo
                 enddo
        
                if(trans) then
                   intYY(ipw,itu) = intYY(ipw,itu) - fact*val
                else
                   intYY(itu,ipw) = intYY(itu,ipw) - fact*val
                endif
               endif
           enddo   

         endif
      enddo 

   enddo
enddo

!! 4th term!
!! really big one? 
do l=1,BDimX

   ip = BIndN(1,l)
   iw = BIndN(2,l)
   ipw = posB(ip,iw)

      do k=1,ADimX
 
         it = AIndN(1,k)
         iu = AIndN(2,k)
         itu = posA(it,iu)
         fact = 2d0*(AOcc(it)-AOcc(iu))

         if(trans) then
            intYY(ipw,itu) = intYY(ipw,itu) + fact*Amat(iw,it)*Smat(iu,ip)
         else
            intYY(itu,ipw) = intYY(itu,ipw) + fact*Amat(iw,it)*Smat(iu,ip)
         endif

      enddo
enddo

!print*, 'intYY',norm2(intYY)

!!!! test 4th term -- very bizzarre!!!
!intYY = 0
!ints = 0
!kl = 0
!do iq=1,dimOB
!   do it=1,NBas
!      kl = kl + 1
!      read(iunit,rec=kl) work(1:NBas*dimOB)
!
!      call ints_modify(NBas,dimOB,ints,NBas,work,&
!                      Smat(it,iq)/AXELE,Vabb,Vbab(it,iq)/BXELE,Sxx)
!
!       do iu=1,NBas
!
!          itu=posA(it,iu)
!        
!          if(itu/=0) then
!             
!             fact = 2d0*(AOcc(it)-AOcc(iu))
!             do i=1,BDimX
!                ip = BIndN(1,i)
!                iw = BIndN(2,i)
!                ipw=posB(ip,iw) 
!
!                val = 0
!                do is=1,dimOB
!                   do ir=1,dimOB
!                      val = val + RDM2val(ir,is,iq,iw)*ints(ir,is)*Smat(iu,ip)
!                    enddo
!                enddo 
!
!                intYY(itu,ipw) = intYY(itu,ipw) + fact*val
!
!             enddo
!          endif
!       enddo
!   enddo
!enddo
!print*, 'intYY',norm2(intYY)

close(iunit)
! term2 tested

deallocate(ints,work,tmp)

end subroutine inter_A2_YY

subroutine inter_A2_XX(dim1,dim2,intXX,RDM2val,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intXX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               RDM2val(dimOB,dimOB,dimOB,dimOB),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,is,it,iu,iw,itu,ipw
double precision :: fact,val
double precision :: Sxx(NBas,NBas),Amat(dimOB,dimOA),Bmat(NBas,NBas)
double precision,allocatable :: tmp(:,:,:,:),work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOB,dimOB,dimOB,NBas),work(NBas*NBas),ints(NBas,NBas))

! test sym
!do l=1,dimOB
!   do k=1,dimOB
!      print*, k,l,norm2(RDM2val(:,:,k,l)-transpose(RDM2val(:,:,l,k)))
!   enddo
!enddo

!tmp = 0
! test for o3v matrix
!do is=1,dimOB
!   do it=1,NBas
!      do iq=1,dimOB
!         do iw=1,dimOB
!            do ir=1,dimOB
!               tmp(ir,iw,iq,it) = tmp(ir,iw,iq,it) &
!                                + RDM2val(ir,iw,iq,is)*Smat(it,is)
!                                !+  FRDM2(ir,iq,iw,is,B%RDM2,B%Occ,B%Ind2,B%NAct,NBas)*Smat(it,is)
!            enddo
!         enddo
!      enddo
!   enddo
!enddo
!print*, 'tmp2', norm2(tmp)
call dgemm('N','T',dimOB**3,NBas,dimOB,1d0,RDM2val,dimOB**3,Smat,NBas,0d0,tmp,dimOB**3)

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

         ! terms: 1, 3
         do it=1,NBas

            itu = posA(it,iu)

            if(itu/=0) then

               fact = 2d0*(AOcc(iu)-AOcc(it))

               do i=1,BDimX

                  ip = BIndN(1,i)
                  iw = BIndN(2,i)
                  ipw = posB(ip,iw)

                  val = 0
                  do ir=1,dimOB
                     val = val - tmp(iw,ir,iq,it)*ints(ip,ir)
                  enddo

                  if(ip<=dimOB) then
                    !print*, 'aaaa?'
                    do ir=1,dimOB
                       val = val + tmp(ir,ip,iq,it)*ints(ir,iw)
                    enddo
                  endif

                  if(trans) then
                     intXX(ipw,itu) = intXX(ipw,itu) + fact*val
                  else
                     intXX(itu,ipw) = intXX(itu,ipw) + fact*val
                  endif

               enddo
            endif
         enddo

         ! term: 4 
         do ip=1,dimOB
            do is=1,dimOB
               do ir=1,dimOB
                  Amat(ip,iu) = Amat(ip,iu) + RDM2val(ir,is,ip,iq)*ints(ir,is)
               enddo
            enddo 
         enddo

      endif  

      ! 2nd term
      iu = l
      ip = k

      do iw=1,dimOB

         ipw = posB(ip,iw)
      
         if(ipw/=0) then

            do it=1,NBas 

               itu = posA(it,iu)

               if(itu/=0) then

                  fact = 2d0*(AOcc(iu)-AOcc(it))
                  val = 0 
                  do iq=1,dimOB
                     do ir=1,dimOB
                        val = val + tmp(iq,ir,iw,it)*ints(ir,iq)
                     enddo
                  enddo

                  if(trans) then
                     intXX(ipw,itu) = intXX(ipw,itu) - fact*val
                  else
                     intXX(itu,ipw) = intXX(itu,ipw) - fact*val
                  endif

               endif

            enddo
         endif
      enddo

   enddo
enddo
close(iunit)
!close(iunit,status='DELETE')

!! 4th term
do l=1,BDimX

   ip = BIndN(1,l)
   iw = BIndN(2,l)
   ipw = posB(ip,iw)

   if(ip<=dimOB) then
      do k=1,ADimX
 
         it = AIndN(1,k)
         iu = AIndN(2,k)
         itu = posA(it,iu)
         fact = 2d0*(AOcc(iu)-AOcc(it))

         if(trans) then
            intXX(ipw,itu) = intXX(ipw,itu) + fact*Amat(ip,iu)*Smat(it,iw)
         else
            intXX(itu,ipw) = intXX(itu,ipw) + fact*Amat(ip,iu)*Smat(it,iw)
         endif

      enddo
   endif
enddo

!print*,'intXX:', norm2(intXX)

deallocate(ints,work,tmp)

end subroutine inter_A2_XX

end module exmisc
