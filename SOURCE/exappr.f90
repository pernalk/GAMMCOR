module exappr 
use exmisc

implicit none

contains

subroutine app_A3_XY(dim1,dim2,intXY,AOcc,BOcc,AFmat,BFmat,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intXY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
                               AFmat(NBas,NBas),BFmat(NBas,NBas),&
                               Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
double precision :: factA,factB,fact,val,val2,valTr
double precision :: prefac,factC
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas),&
                    Amat(NBas,NBas),Emat(NBas,dimOA)
double precision,allocatable :: work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

Amat = 0
do ib=1,dimOB
   do ir=1,NBas
      do ip=1,NBas
         Amat(ip,ir) = Amat(ip,ir) + BOcc(ib)*Sab(ip,ib)*Sab(ir,ib)
      enddo
   enddo
enddo

valTr = 0
do iq=1,dimOA
   valTr = valTr + AOcc(iq)*Amat(iq,iq)
enddo

Sbb = 0
do ic=1,NBas
   do ia=1,NBas
      do iq=1,dimOA
         Sbb(ia,ic) = Sbb(ia,ic) + AOcc(iq)*Sab(iq,ia)*Sab(iq,ic)
      enddo
   enddo
enddo

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
Emat = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ic = l
      ia = k

      iac = posB(ia,ic)

      if(iac/=0) then

         do i=1,ADimX
            ip = AIndN(1,i)
            ir = AIndN(2,i)
            ipr = posA(ip,ir)
            
            fact  = 8d0*(AOcc(ir)-AOcc(ip))*(BOcc(ia)-BOcc(ic))
            factB = 4d0*(BOcc(ia)-BOcc(ic))

            ! 1A-1B
            intXY(ipr,iac) = intXY(ipr,iac) - fact*valTr*ints(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + AOcc(iq)*ints(iq,iq)  
            enddo

            !2A-1B 
            intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Amat(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Amat(iq,ir)
            enddo

            !3A-1B
            intXY(ipr,iac) = intXY(ipr,iac) - factB*val

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Amat(ip,iq)
            enddo

            !4A-1B
            intXY(ipr,iac) = intXY(ipr,iac) - factB*val

         enddo

      endif

      ! all 2B terms 
      if(k==l.and.k<=dimOB) then

         ib = k
         Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + BOcc(ib)*ints(1:NBas,1:dimOA)

      endif

      ib = l
      ia = k 
   
      val = 0 
      do iq=1,dimOA
         val = val + AOcc(iq)*ints(iq,iq)
      enddo
 
      do ic=1,dimOB

         iac = posB(ia,ic)

         if(iac/=0) then

            factB = 2d0*(BFmat(ib,ic)-BFmat(ib,ia))

            do i=1,ADimX

               ip = AIndN(1,i)
               ir = AIndN(2,i)
               ipr = posA(ip,ir)

               fact = 2d0*factB*(AOcc(ir)-AOcc(ip))

               !1A-3B 
               intXY(ipr,iac) = intXY(ipr,iac) - fact*Sbb(ic,ib)*ints(ip,ir)
               
               !2A-3B 
               intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ib)*Sab(ir,ic)

               val2 = 0
               do iq=1,dimOA
                  val2 = val2 + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ib)
               enddo

               !3A-3B 
               intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ir,ic)

               val2 = 0
               do iq=1,dimOA
                  val2 = val2 + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ic)
               enddo

               !4A-3B 
               intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ip,ib)

            enddo 
         endif
      enddo

      if(k<=dimOB) then

         ic = l
         ib = k 

         val = 0 
         do iq=1,dimOA
            val = val + AOcc(iq)*ints(iq,iq)
         enddo

         do ia=1,NBas

            iac = posB(ia,ic)

            if(iac/=0) then

               factB = 2d0*(BFmat(ic,ib)-BFmat(ia,ib))

               do i=1,ADimX

                  ip = AIndN(1,i)
                  ir = AIndN(2,i)
                  ipr = posA(ip,ir)

                  fact = 2d0*factB*(AOcc(ir)-AOcc(ip)) 

                  !1A-4B  
                  intXY(ipr,iac) = intXY(ipr,iac) - fact*ints(ip,ir)*Sbb(ia,ib)

                  !2A-4B 
                  intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ia)*Sab(ir,ib)

                  val2 = 0
                  do iq=1,dimOA
                     val2 = val2 + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ia)
                  enddo

                  !3A-4B
                  intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ir,ib)

                  val2 = 0 
                  do iq=1,dimOA
                     val2 = val2 + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ib)
                  enddo

                  !4A-4B
                  intXY(ipr,iac) = intXY(ipr,iac) - factB*val2*Sab(ip,ia)

               enddo
            endif

         enddo
      endif

   enddo
enddo

close(iunit)

prefac = 0
do iq=1,dimOA
   prefac = prefac + AOcc(iq)*Emat(iq,iq)
enddo

! 2B terms 
do l=1,BDimX

   ia = BIndN(1,l)
   ic = BIndN(2,l)
   iac = posB(ia,ic)

   factB = 4d0*(BOcc(ia)-BOcc(ic))

   do k=1,ADimX
      ip = AIndN(1,k)
      ir = AIndN(2,k)
      ipr = posA(ip,ir)

      fact  = 2d0*factB*(AOcc(ir)-AOcc(ip))
      factC = prefac*fact

      !1A-2B
      intXY(ipr,iac) = intXY(ipr,iac) - fact*Emat(ip,ir)*Sbb(ia,ic)

      !val = 0
      !do iq=1,dimOA
      !   val = val + AOcc(iq)*Emat(iq,iq)
      !enddo

      !2A-2B
      !intXY(ipr,iac) = intXY(ipr,iac) - fact*val*Sab(ip,ia)*Sab(ir,ic)
      intXY(ipr,iac) = intXY(ipr,iac) - factC*Sab(ip,ia)*Sab(ir,ic)

      val = 0
      do iq=1,dimOA
         val = val + (AFmat(ip,iq)-AFmat(ir,iq))*Emat(ip,iq)*Sab(iq,ia)
      enddo

      !3A-2B
      intXY(ipr,iac) = intXY(ipr,iac) - factB*val*Sab(ir,ic)

      val = 0
      do iq=1,dimOA
         val = val + (AFmat(iq,ip)-AFmat(iq,ir))*Emat(iq,ir)*Sab(iq,ic)
      enddo

      !4A-2B
      intXY(ipr,iac) = intXY(ipr,iac) - factB*val*Sab(ip,ia)

   enddo 
enddo

deallocate(ints,work)

print*, 'appA3_XY', norm2(intXY)

end subroutine app_A3_XY

subroutine app_nn_A3_XX(dim1,dim2,intXX,AOcc,BOcc,AFmat,BFmat,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intXX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
                               AFmat(NBas,NBas),BFmat(NBas,NBas),&
                               Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
double precision :: factA,factB,fact,val,val3B,val4B,valTr
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas),&
                    Amat(NBas,NBas),Emat(NBas,dimOA)
double precision,allocatable :: work(:),ints(:,:)
! test
double precision :: preQ,preQa,factAn,factBn
double precision,allocatable :: AOccN(:),BOccN(:),  &
                                SNab(:,:),intV(:,:),&
                                intE(:,:),intU(:,:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(AOccN(NBas),BOccN(NBas))

! BB functional
do i=1,NBas
   AOccN(i) = sqrt(AOcc(i))
   BOccN(i) = sqrt(BOcc(i))
enddo

allocate(SNab(NBas,NBas),intV(NBas,NBas),intE(NBas,NBas))
allocate(intU(ADimX,dimOA,dimOB))

SNab = 0
do ib=1,dimOB
   do iq=1,dimOA
      SNab(iq,ib) = SNab(iq,ic) + AOccN(iq)*Sab(iq,ib)*BOccN(ib)
   enddo
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

Amat = 0
do ib=1,dimOB
   do ir=1,NBas !dimOA
      do ip=1,NBas
         Amat(ip,ir) = Amat(ip,ir) + BOcc(ib)*Sab(ip,ib)*Sab(ir,ib)
      enddo
   enddo
enddo

valTr = 0
do iq=1,dimOA
   valTr = valTr + AOcc(iq)*Amat(iq,iq)
enddo

Sbb = 0
do ic=1,NBas
   do ia=1,NBas
      do iq=1,dimOA
         Sbb(ia,ic) = Sbb(ia,ic) + AOcc(iq)*Sab(iq,ia)*Sab(iq,ic)
      enddo
   enddo
enddo

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
Emat = 0
intV = 0
intU = 0
kl   = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      !ints(dimOA+1:NBas,dimOA+1:NBas) = 0

      ic = l
      ia = k

      iac = posB(ia,ic)

      if(iac/=0) then

         factB = 4d0*(BOcc(ic)-BOcc(ia))

         do i=1,ADimX
            ip = AIndN(1,i)
            ir = AIndN(2,i)
            ipr = posA(ip,ir)
            
            fact  = 2d0*(AOcc(ir)-AOcc(ip))*factB

            ! 1A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - fact*valTr*ints(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + AOcc(iq)*ints(iq,iq)
            enddo

            !2A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - fact*val*Amat(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Amat(iq,ir)
            enddo

            !3A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - factB*val

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Amat(ip,iq)
            enddo

            !4A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - factB*val

         enddo

      endif

      ! all 2B terms 
      if(k==l) then

         ib = k
         Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + BOcc(ib)*ints(1:NBas,1:dimOA)

      endif

      ib = l
      ia = k
   
      val3B = 0 
      do iq=1,dimOA
         val3B = val3B + AOcc(iq)*ints(iq,iq)
      enddo
 
      do ic=1,dimOB

         iac = posB(ia,ic)

         if(iac/=0) then

            factB = 2d0*(BFmat(ia,ib)-BFmat(ic,ib))

            do i=1,ADimX

               ip = AIndN(1,i)
               ir = AIndN(2,i)
               ipr = posA(ip,ir)

               fact = 2d0*factB*(AOcc(ir)-AOcc(ip))

               !1A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - fact*Sbb(ic,ib)*ints(ip,ir)
               
               !2A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - fact*val3B*Sab(ip,ic)*Sab(ir,ib)

               val = 0
               do iq=1,dimOA
                  val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ic)
               enddo

               !3A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ib)

               val = 0
               do iq=1,dimOA
                  val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ib)
               enddo

               !4A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ic)

            enddo
         endif
      enddo

      if(k<=dimOB) then

         ic = l
         ib = k

         val4B = 0
         do iq=1,dimOA
            val4B = val4B + AOcc(iq)*ints(iq,iq)
         enddo

         ! new: 3A-4B
         call dgemv('N',NBas,dimOA,1d0,ints,NBas,SNab(:,ib),1,1d0,intV(:,ic),1)
         !do ip=1,NBas
            !val = 0
            !do iq=1,dimOA
               !val = val + ints(ip,iq)*SNab(iq,ib)
            !enddo
            !intV(ip,ic) = intV(ip,ic) + val
         !enddo

         ! new: 4A-4B
         do i=1,ADimX

            ip = AIndN(1,i)
            ir = AIndN(2,i)
            ipr = posA(ip,ir)

            do iq=1,dimOA 
               intU(ipr,iq,ic) = intU(ipr,iq,ic) + BOccN(ib)*ints(iq,ir)*Sab(ip,ib)
            enddo

         enddo


         do ia=1,NBas

            iac = posB(ia,ic)

            if(iac/=0) then

               factB = 2d0*(BFmat(ib,ia)-BFmat(ib,ic))

               do i=1,ADimX

                  ip = AIndN(1,i)
                  ir = AIndN(2,i)
                  ipr = posA(ip,ir)

                  fact = 2d0*factB*(AOcc(ir)-AOcc(ip))

                  !1A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - fact*ints(ip,ir)*Sbb(ib,ia)

                  !2A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - fact*val4B*Sab(ip,ib)*Sab(ir,ia)

               enddo
            endif

         enddo
      endif

   enddo
enddo

close(iunit)

! 2B terms

preQ = 0
do iq=1,dimOA
   preQ = preQ + AOcc(iq)*Emat(iq,iq)
enddo

intE = 0
do ic=1,NBas
   do ip=1,NBas
      do iq=1,dimOA
         intE(ip,ic) = intE(ip,ic) - AOccN(iq)*Sab(iq,ic)*Emat(ip,iq)
      enddo
   enddo
enddo

do l=1,BDimX

   ia = BIndN(1,l)
   ic = BIndN(2,l)
   iac = posB(ia,ic)

   factB  = 4d0*(BOcc(ic)-BOcc(ia))
   factBn = 2d0*(BOccN(ic)-BOccN(ia))

   do k=1,ADimX

      ip = AIndN(1,k)
      ir = AIndN(2,k)
      ipr = posA(ip,ir)

      fact   = 2d0*(AOcc(ir)-AOcc(ip))*factB
      factAn = (AOccN(ir)-AOccN(ip))

      !1A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - fact*Emat(ip,ir)*Sbb(ic,ia)

      !2A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - preQ*fact*Sab(ip,ic)*Sab(ir,ia)

      !3A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - factAn*factB*intE(ip,ic)*Sab(ir,ia)

      !4A-2B
      ! not sure why this works - ?
      intXX(ipr,iac) = intXX(ipr,iac) - factAn*factB*intE(ir,ia)*Sab(ip,ic)

      ! new: 3A-4B
      intXX(ipr,iac) = intXX(ipr,iac) - factAn*factBn*intV(ip,ic)*Sab(ir,ia)

      val = 0
      do iq=1,dimOA
         val = val + AOccN(iq)*intU(ipr,iq,ic)*Sab(iq,ia)
      enddo

      ! new: 4A-4B
      intXX(ipr,iac) = intXX(ipr,iac) - factAn*factBn*val

   enddo 
enddo

deallocate(ints,work)
deallocate(intU)
deallocate(intV,intE)
deallocate(BOccN,AOccN)

print*, 'appA3_nn_XX',norm2(intXX(1:dim1,1:dim2))

end subroutine app_nn_A3_XX

!check: a) is (OO|OO) here enough?
!       b) Fpq = n_p^a * n_q^a
subroutine app_A3_XX(dim1,dim2,intXX,AOcc,BOcc,AFmat,BFmat,Sab, &
                    AXELE,Vabb,BXELE,Vbaa,IntKFile,&
                    BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
character(*) :: IntKFile
double precision,intent(inout) :: intXX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),BOcc(NBas),&
                               AFmat(NBas,NBas),BFmat(NBas,NBas),&
                               Sab(NBas,NBas),Vabb(NBas,NBas),Vbaa(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,ia,ib,ic,iac,ipr
double precision :: factA,factB,fact,val,val3B,val4B,valTr
double precision :: prefac,factC
double precision :: Sxx(NBas,NBas),Sbb(NBas,NBas),&
                    Amat(NBas,NBas),Emat(NBas,dimOA)
double precision,allocatable :: work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(work(NBas*NBas),ints(NBas,NBas))

Amat = 0
do ib=1,dimOB
   do ir=1,NBas !dimOA
      do ip=1,NBas
         Amat(ip,ir) = Amat(ip,ir) + BOcc(ib)*Sab(ip,ib)*Sab(ir,ib)
      enddo
   enddo
enddo

valTr = 0
do iq=1,dimOA
   valTr = valTr + AOcc(iq)*Amat(iq,iq)
enddo

Sbb = 0
do ic=1,NBas
   do ia=1,NBas
      do iq=1,dimOA
         Sbb(ia,ic) = Sbb(ia,ic) + AOcc(iq)*Sab(iq,ia)*Sab(iq,ic)
      enddo
   enddo
enddo

!(FO|FO): (AA|BB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOA)

! one loop over integrals
ints = 0
Emat = 0
kl = 0
do l=1,dimOB
   do k=1,NBas
      kl = kl + 1
      read(iunit,rec=kl) work(1:NBas*dimOA)

      call ints_modify(NBas,dimOA,ints,NBas,work,&
                       Vabb(k,l)/AXELE,Sxx,Sxx(k,l)/BXELE,Vbaa)

      ic = l
      ia = k

      iac = posB(ia,ic)

      if(iac/=0) then

         factB = 4d0*(BOcc(ic)-BOcc(ia))

         do i=1,ADimX
            ip = AIndN(1,i)
            ir = AIndN(2,i)
            ipr = posA(ip,ir)
            
            fact  = 2d0*(AOcc(ir)-AOcc(ip))*factB

            ! 1A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - fact*valTr*ints(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + AOcc(iq)*ints(iq,iq)
            enddo

            !2A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - fact*val*Amat(ip,ir)

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Amat(iq,ir)
            enddo

            !3A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - factB*val

            val = 0
            do iq=1,dimOA
               val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Amat(ip,iq)
            enddo

            !4A-1B
            intXX(ipr,iac) = intXX(ipr,iac) - factB*val

         enddo

      endif

      ! all 2B terms 
      if(k==l) then

         ib = k
         Emat(1:NBas,1:dimOA) = Emat(1:NBas,1:dimOA) + BOcc(ib)*ints(1:NBas,1:dimOA)

      endif

      ib = l
      ia = k
   
      val3B = 0 
      do iq=1,dimOA
         val3B = val3B + AOcc(iq)*ints(iq,iq)
      enddo
 
      do ic=1,dimOB

         iac = posB(ia,ic)

         if(iac/=0) then

            factB = 2d0*(BFmat(ia,ib)-BFmat(ic,ib))

            do i=1,ADimX

               ip = AIndN(1,i)
               ir = AIndN(2,i)
               ipr = posA(ip,ir)

               fact = 2d0*factB*(AOcc(ir)-AOcc(ip))

               !1A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - fact*Sbb(ic,ib)*ints(ip,ir)
               
               !2A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - fact*val3B*Sab(ip,ic)*Sab(ir,ib)

               val = 0
               do iq=1,dimOA
                  val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ic)
               enddo

               !3A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ib)

               val = 0
               do iq=1,dimOA
                  val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ib)
               enddo

               !4A-3B
               intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ic)

            enddo
         endif
      enddo

      if(k<=dimOB) then

         ic = l
         ib = k

         val4B = 0
         do iq=1,dimOA
            val4B = val4B + AOcc(iq)*ints(iq,iq)
         enddo

         do ia=1,NBas

            iac = posB(ia,ic)

            if(iac/=0) then

               factB = 2d0*(BFmat(ib,ia)-BFmat(ib,ic))

               do i=1,ADimX

                  ip = AIndN(1,i)
                  ir = AIndN(2,i)
                  ipr = posA(ip,ir)

                  fact = 2d0*factB*(AOcc(ir)-AOcc(ip))

                  !1A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - fact*ints(ip,ir)*Sbb(ib,ia)

                  !2A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - fact*val4B*Sab(ip,ib)*Sab(ir,ia)

                  val = 0
                  do iq=1,dimOA
                     val = val + (AFmat(ip,iq)-AFmat(ir,iq))*ints(ip,iq)*Sab(iq,ib)
                  enddo

                  !3A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ia)

                  val = 0
                  do iq=1,dimOA
                     val = val + (AFmat(iq,ip)-AFmat(iq,ir))*ints(iq,ir)*Sab(iq,ia)
                  enddo

                  !4A-4B
                  intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ib)

               enddo
            endif

         enddo
      endif

   enddo
enddo

close(iunit)

prefac = 0
do iq=1,dimOA
   prefac = prefac + AOcc(iq)*Emat(iq,iq)
enddo

! 2B terms
do l=1,BDimX

   ia = BIndN(1,l)
   ic = BIndN(2,l)
   iac = posB(ia,ic)

   factB = 4d0*(BOcc(ic)-BOcc(ia))

   do k=1,ADimX

      ip = AIndN(1,k)
      ir = AIndN(2,k)
      ipr = posA(ip,ir)

      fact  = 2d0*(AOcc(ir)-AOcc(ip))*factB
      factC = prefac*fact

      !1A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - fact*Emat(ip,ir)*Sbb(ic,ia)

      !val = 0
      !do iq=1,dimOA
      !   val = val + AOcc(iq)*Emat(iq,iq)
      !enddo

      !2A-2B
      !intXX(ipr,iac) = intXX(ipr,iac) - fact*val*Sab(ip,ic)*Sab(ir,ia)
      intXX(ipr,iac) = intXX(ipr,iac) - factC*Sab(ip,ic)*Sab(ir,ia)

      val = 0
      do iq=1,dimOA
         val = val + (AFmat(ip,iq)-AFmat(ir,iq))*Emat(ip,iq)*Sab(iq,ic)
      enddo

      !3A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ir,ia)

      val = 0
      do iq=1,dimOA
         val = val + (AFmat(iq,ip)-AFmat(iq,ir))*Emat(iq,ir)*Sab(iq,ia)
      enddo

      !4A-2B
      intXX(ipr,iac) = intXX(ipr,iac) - factB*val*Sab(ip,ic)

   enddo 
enddo

deallocate(ints,work)

print*, 'appA3_XX',norm2(intXX(1:dim1,1:dim2))

end subroutine app_A3_XX


subroutine app_A2_YX(dim1,dim2,intYX,BOcc,Fmat,Smat,  &
                  AXELE,Vabb,BXELE,Vbab,IntKFile,IntJFile,&
                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),&
                      posA(NBas,NBas),posB(NBas,NBas)
logical,intent(in) :: trans
character(*) :: IntKFile,IntJFile
double precision,intent(inout) :: intYX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               BOcc(NBas),Fmat(NBas,NBas), & 
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit,i,j,k,l,kl
integer :: ip,iq,ir,is,it,iu,itu,ipr
integer :: iw,ipw
double precision :: fact,factA,factF,val
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: work(:),ints(:,:),tmp(:,:)
double precision,external  :: FRDM2

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOB,dimOB),work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
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
      do i=1,BDimX
 
         ip = BIndN(1,i)
         ir = BIndN(2,i)
         ipr = posB(ip,ir)

         fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)*ints(ip,ir)
         factF = (Fmat(ip,iq)-Fmat(ir,iq))*ints(ip,iq)

         do iu=1,dimOA
 
            itu = posA(it,iu)
          
            if(itu/=0) then

               factA = 2d0*(AOcc(it)-AOcc(iu))

               if(trans) then 
                  intYX(ipr,itu) = intYX(ipr,itu) - factA*fact*Smat(iu,iq) &
                                                  - factA*factF*Smat(iu,ir)
               else
                  intYX(itu,ipr) = intYX(itu,ipr) - factA*fact*Smat(iu,iq) &
                                                  - factA*factF*Smat(iu,ir)
               endif

            endif
         enddo
      enddo

   enddo
enddo

close(iunit)

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

      ! term 4
      if(trans) then
         do i=1,ADimX

            it = AIndN(1,i)
            iu = AIndN(2,i)
            itu = posA(it,iu)

            factA = 2d0*(AOcc(it)-AOcc(iu))

            do ip=1,NBas

               ipr = posB(ip,ir)

               if(ipr/=0) then

                  factF = Fmat(iq,ip)-Fmat(iq,ir)
                  intYX(ipr,itu) = intYX(ipr,itu) - factA*factF*Smat(iu,iq)*ints(it,ip)

               endif

            enddo
         enddo

      else

         do ip=1,NBas

            ipr = posB(ip,ir)

            if(ipr/=0) then

               factF = Fmat(iq,ip)-Fmat(iq,ir)

               do i=1,ADimX

                  it = AIndN(1,i)
                  iu = AIndN(2,i)
                  itu = posA(it,iu)

                  factA = 2d0*(AOcc(it)-AOcc(iu))
                  intYX(itu,ipr) = intYX(itu,ipr) - factA*factF*Smat(iu,iq)*ints(it,ip)

               enddo
            endif
         enddo
      endif 

     ! term 2
     if(k==l) then

        iq = l

        if(trans) then

           do j=1,ADimX
              it = AIndN(1,j)
              iu = AIndN(2,j)
              itu = posA(it,iu)

              factA = 2d0*(AOcc(it)-AOcc(iu))

              do i=1,BDimX
                 ip = BIndN(1,i)
                 ir = BIndN(2,i)
                 ipr = posB(ip,ir)

                 fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)

                 intYX(ipr,itu) = intYX(ipr,itu) - fact*factA*Smat(iu,ir)*ints(it,ip)

              enddo
           enddo   

        else
           do i=1,BDimX
              ip = BIndN(1,i)
              ir = BIndN(2,i)
              ipr = posB(ip,ir)

              fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)

              do j=1,ADimX
                 it = AIndN(1,j)
                 iu = AIndN(2,j)
                 itu = posA(it,iu)

                 factA = 2d0*(AOcc(it)-AOcc(iu))

                 intYX(itu,ipr) = intYX(itu,ipr) - fact*factA*Smat(iu,ir)*ints(it,ip)

              enddo
           enddo
        endif

     endif 

   enddo
enddo

close(iunit)

deallocate(ints,work,tmp)

print*,'appYX', norm2(intYX)

end subroutine app_A2_YX 

subroutine app_A2_XY(dim1,dim2,intXY,BOcc,Fmat,Smat,  &
                  AXELE,Vabb,BXELE,Vbab,IntKFile,&
                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intXY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               BOcc(NBas),Fmat(NBas,NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
double precision :: factA,factF,fact,val
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(dimOA,dimOB),work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
tmp  = 0
ints = 0
kl = 0
do l=1,dimOB
   do k=1,dimOA
      kl = kl + 1
      read(iunit,rec=k+(l-1)*NBas) work(1:NBas*dimOB)

      call ints_modify(NBas,dimOB,ints,NBas,work,&
                       Smat(k,l)/AXELE,Vabb,Vbab(k,l)/BXELE,Sxx)

      iq = l
      iu = k

      ! terms 1 and 3  
      do i=1,BDimX
          
         ip = BIndN(1,i)
         ir = BIndN(2,i)
         ipr = posB(ip,ir)
         
         fact = 2d0*(BOcc(ip)-BOcc(ir))*BOcc(iq)*ints(ip,ir)
         factF = (Fmat(ir,iq)-Fmat(ip,iq))*ints(ir,iq)

         do it=1,NBas
            
            itu = posA(it,iu)

            if(itu/=0) then
             
               factA = 2d0*(AOcc(iu)-AOcc(it))

               if(trans) then
                  intXY(ipr,itu) = intXY(ipr,itu) - factA*fact*Smat(it,iq) &
                                                  - factA*factF*Smat(it,ip)
               else                               
                  intXY(itu,ipr) = intXY(itu,ipr) - factA*fact*Smat(it,iq) &
                                                  - factA*factF*Smat(it,ip)
               endif

            endif
         enddo
      enddo

      ! term 2
      ir = l
      iu = k
     
      do iq=1,dimOB
         tmp(iu,ir) = tmp(iu,ir) + BOcc(iq)*ints(iq,iq)
      enddo

      ! term4
      do it=1,NBas

         itu = posA(it,iu) 

         if(itu/=0) then

            factA = 2d0*(AOcc(iu)-AOcc(it))
 
            do ip=1,NBas

               ipr = posB(ip,ir)
 
               if(ipr/=0) then

                  if(trans) then            
                     do iq=1,dimOB
                        intXY(ipr,itu) = intXY(ipr,itu) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
                                       * Smat(it,iq)*ints(ip,iq)
                     enddo
                  else
                     do iq=1,dimOB
                        intXY(itu,ipr) = intXY(itu,ipr) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
                                       * Smat(it,iq)*ints(ip,iq)
                     enddo
                  endif

               endif
            enddo 
         endif
      enddo

   enddo
enddo

close(iunit)

if(trans) then

   do l=1,ADimX
      
      it = AIndN(1,l)
      iu = AIndN(2,l)
      itu = posA(it,iu)

      factA = 2d0*(AOcc(iu)-AOcc(it))
 
      do k=1,BDimX
       
         ip = BIndN(1,k)
         ir = BIndN(2,k)
         ipr = posB(ip,ir)

         fact = 2d0*(BOcc(ip)-BOcc(ir))

         intXY(ipr,itu) = intXY(ipr,itu) - fact*factA*tmp(iu,ir)*Smat(it,ip)

      enddo
   enddo

else

   do l=1,BDimX
      
      ip = BIndN(1,l)
      ir = BIndN(2,l)
      ipr = posB(ip,ir)

      fact = 2d0*(BOcc(ip)-BOcc(ir)) 

      do k=1,ADimX
   
         it = AIndN(1,k)
         iu = AIndN(2,k)
         itu = posA(it,iu)

         factA = 2d0*(AOcc(iu)-AOcc(it))
         intXY(itu,ipr) = intXY(itu,ipr) - fact*factA*tmp(iu,ir)*Smat(it,ip)
   
      enddo
   enddo

endif

deallocate(ints,work,tmp)

print*, 'appXY',norm2(intXY)

end subroutine app_A2_XY

subroutine app_A2_YY(dim1,dim2,intYY,BOcc,Fmat,Smat,  &
                  AXELE,Vabb,BXELE,Vbab,IntKFile,&
                  AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intYY(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               BOcc(NBas),Fmat(NBas,NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
double precision :: factA,factF,fact,val
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(NBas,dimOB),work(NBas*NBas),ints(NBas,NBas))

!(FO|FO):(BB|AB) or (AA|BA)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
tmp  = 0
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

      ! term 1 and 3
      do i=1,BDimX

         ip = BIndN(1,i) 
         ir = BIndN(2,i) 
         ipr = posB(ip,ir)

         fact = 2d0*BOcc(iq)*(BOcc(ip)-BOcc(ir))*ints(ip,ir)
         factF = (Fmat(ir,iq)-Fmat(ip,iq))*ints(ir,iq)

         do iu=1,dimOA

            itu = posA(it,iu)

            if(itu/=0) then

              factA = 2d0*(AOcc(it)-AOcc(iu))
              if(trans) then 
                 intYY(ipr,itu) = intYY(ipr,itu) - fact*factA*Smat(iu,iq) &
                                                 - factF*factA*Smat(iu,ip)
              else
                 intYY(itu,ipr) = intYY(itu,ipr) - fact*factA*Smat(iu,iq) &
                                                 - factF*factA*Smat(iu,ip)
              endif
              
            endif
         enddo
      enddo

      ir = l 
      it = k

      ! term 2
      do iq=1,dimOB  
         tmp(it,ir) = tmp(it,ir) + BOcc(iq)*ints(iq,iq)
      enddo

      ! term 4
      do iu=1,dimOA
         
         itu = posA(it,iu)

         if(itu/=0) then
            
            factA = 2d0*(AOcc(it)-AOcc(iu))

            do ip=1,NBas

               ipr=posB(ip,ir)

               if(ipr/=0) then

                  if(trans) then

                     do iq=1,dimOB
                        intYY(ipr,itu) = intYY(ipr,itu) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
                                       * Smat(iu,iq)*ints(ip,iq) 
                     enddo
 
                  else

                     do iq=1,dimOB
                        intYY(itu,ipr) = intYY(itu,ipr) - factA*(Fmat(iq,ir)-Fmat(iq,ip)) &
                                       * Smat(iu,iq)*ints(ip,iq) 
                     enddo

                  endif
               endif
            enddo

         endif
      enddo

   enddo
enddo

close(iunit)

if(trans) then

   do l=1,ADimX

      it = AIndN(1,l)
      iu = AIndN(2,l)
      itu = posA(it,iu)

      factA = 2d0*(AOcc(it)-AOcc(iu))

      do k=1,BDimX
         ip = BIndN(1,k)
         ir = BIndN(2,k)
         ipr = posB(ip,ir)

         fact = 2d0*(BOcc(ip)-BOcc(ir))
         intYY(ipr,itu) = intYY(ipr,itu) - fact*factA*tmp(it,ir)*Smat(iu,ip)

      enddo
   enddo

else

   do l=1,BDimX

      ip = BIndN(1,l)
      ir = BIndN(2,l)
      ipr = posB(ip,ir)
   
      fact = 2d0*(BOcc(ip)-BOcc(ir))
   
      do k=1,ADimX

         it = AIndN(1,k) 
         iu = AIndN(2,k) 
         itu = posA(it,iu)
         factA = 2d0*(AOcc(it)-AOcc(iu))
   
         intYY(itu,ipr) = intYY(itu,ipr) - fact*factA*tmp(it,ir)*Smat(iu,ip)
   
      enddo
   enddo

endif

deallocate(ints,work,tmp)

print*, 'appYY',norm2(intYY)

end subroutine app_A2_YY

subroutine app_A2_XX(dim1,dim2,intXX,BOcc,Fmat,Smat,  &
                    AXELE,Vabb,BXELE,Vbab,IntKFile,&
                    AOcc,BIndN,AIndN,posB,posA,dimOB,dimOA,BDimX,ADimX,NBas,trans)
implicit none

integer,intent(in) :: dim1,dim2,dimOA,dimOB,ADimX,BDimX,NBas
integer,intent(in) :: AIndN(2,ADimX),BIndN(2,BDimX),posA(NBas,NBas),posB(NBas,NBas) 
logical,intent(in) :: trans
character(*) :: IntKFile
double precision,intent(inout) :: intXX(dim1,dim2)
double precision,intent(in) :: AXELE,BXELE,AOcc(NBas),&
                               BOcc(NBas),Fmat(NBas,NBas),&
                               Smat(NBas,NBas),Vabb(NBas,NBas),Vbab(NBas,NBas)

integer :: iunit
integer :: i,k,l,kl,ip,iq,ir,is,it,iu,itu,ipr
double precision :: factA,factF,fact,val
double precision :: Sxx(NBas,NBas)
double precision,allocatable :: tmp(:,:),work(:),ints(:,:)

Sxx = 0
do i=1,NBas
   Sxx(i,i) = 1d0
enddo

allocate(tmp(NBas,dimOA),work(NBas*NBas),ints(NBas,NBas))

!do i=1,NBas
!    print*, i,posA(i,i),posB(i,i)
!enddo

!(FO|FO):(BB|BA) or (AA|AB)
open(newunit=iunit,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBas*dimOB)

! one loop over integrals
ints = 0
tmp  = 0
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
         do i=1,BDimX
            ip = BIndN(1,i)
            ir = BIndN(2,i)
            ipr = posB(ip,ir)

            fact = 2d0*(BOcc(ir)-BOcc(ip))*BOcc(iq)*ints(ip,ir)
            factF = (Fmat(ip,iq)-Fmat(ir,iq))*ints(ip,iq)

 
            do it=1,NBas

               itu = posA(it,iu)

               if(itu/=0) then

                  factA = 2d0*(AOcc(iu)-AOcc(it))

                  if(trans) then
                     intXX(ipr,itu) = intXX(ipr,itu) - factA*fact*Smat(it,iq) &
                                                     - factA*factF*Smat(it,ir)
                  else
                     intXX(itu,ipr) = intXX(itu,ipr) - factA*fact*Smat(it,iq) &
                                                     - factA*factF*Smat(it,ir)
                  endif

               endif
            enddo

         enddo
      endif
   
      iu = l
      ip = k

      ! term 2
      do iq=1,dimOB
         tmp(ip,iu) = tmp(ip,iu) + BOcc(iq)*ints(iq,iq) 
      enddo

      ! term 4
      do ir=1,dimOB

         ipr = posB(ip,ir)
         
         if(ipr/=0) then

            do it=1,NBas
               
               itu = posA(it,iu)

               if(itu/=0) then

                  factA = 2d0*(AOcc(iu)-AOcc(it))

                  val = 0 
                  do iq=1,dimOB
                     val = val - (Fmat(iq,ip)-Fmat(iq,ir)) * &
                           ints(iq,ir)*Smat(it,iq)
                  enddo

                  if(trans) then
                     intXX(ipr,itu) = intXX(ipr,itu) + factA*val
                  else
                     intXX(itu,ipr) = intXX(itu,ipr) + factA*val
                  endif 
  
               endif               
            enddo
         endif
      enddo       

   enddo
enddo

close(iunit)

if(trans) then

   do l=1,ADimX

      it = AIndN(1,l)
      iu = AIndN(2,l)
      itu = posA(it,iu)

      factA = 4d0*(AOcc(iu)-AOcc(it))

      do k=1,BDimX
         ip = BIndN(1,k)
         ir = BIndN(2,k)
         ipr = posB(ip,ir)

         intXX(ipr,itu) = intXX(ipr,itu) - factA*(BOcc(ir)-BOcc(ip))*tmp(ip,iu)*Smat(it,ir)

      enddo
   enddo

else
   do l=1,BDimX
   
      ip = BIndN(1,l)
      ir = BIndN(2,l)
      ipr = posB(ip,ir)
   
      fact = 4d0*(BOcc(ir)-BOcc(ip))
   
      do k=1,ADimX
   
         it = AIndN(1,k)
         iu = AIndN(2,k)
         itu = posA(it,iu)
   
         factA = fact*(AOcc(iu)-AOcc(it))
   
         intXX(itu,ipr) = intXX(itu,ipr) - factA*tmp(ip,iu)*Smat(it,ir)
   
      enddo
   enddo
endif

print*,'appXX:', norm2(intXX)

deallocate(ints,work,tmp)

end subroutine app_A2_XX

end module exappr
