module abchol

implicit none

contains

subroutine JK_Chol_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                   INActive,NOccup,NDim,NDimX,NBasis,NInte1,IntJFile,IntKFile,ACAlpha,AB,ETot)

implicit none

integer,intent(in) :: NDim,NDimX,NBasis,NInte1
integer,intent(in) :: INActive,NOccup,AB
integer,intent(in) :: IGem(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision,intent(inout) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(inout) :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision,intent(in) :: ACAlpha
double precision,intent(in) :: AuxCoeff(3,3,3,3),Occ(NBasis),RDM2val(NOccup,NOccup,NOccup,NOccup)
character(*) :: IntJFile,IntKFile
double precision,intent(inout),optional :: ETot

integer :: i,j,k,l,kk,ll,kl
integer :: ip,iq,ir,is,iu,ipq,irs
integer :: iunit
integer :: iloop,nloop,off,dimFO,dimOO,dimFF
integer :: mloop
integer :: iBatch
integer :: BatchSize,MaxBatchSize = 120
integer :: NCholesky
double precision :: val
double precision :: AuxVal,HNOCoef
double precision,allocatable :: work1(:,:),work2(:,:)
double precision,allocatable :: ints(:,:),MatFF(:,:)

!print*, 'start JK Chol:'

if(AB==1) then
   HNOCoef = -1
elseif(AB==0) then
   HNOCoef = 1 - ACAlpha
elseif(AB==2) then
   HNOCoef = 1
endif

!print*, 'NDimX = ', NDimX
!print*, 'Occ =',norm2(Occ)
!print*, 'HNO =',norm2(HNO)
!print*, 'WMAT=',norm2(WMAT)
!print*, 'AuxI=',norm2(AuxI)
!print*, 'AuxIO',norm2(AuxIO)
!print*, 'RDM2val',norm2(RDM2val)

! read cholesky (FF|K) vectors
open(newunit=iunit,file='cholvecs',form='unformatted')
read(iunit) NCholesky
allocate(MatFF(NCholesky,NBasis**2))
read(iunit) MatFF
close(iunit)

! set number of loops over integrals
dimFO = NOccup*NBasis
nloop = (dimFO - 1) / MaxBatchSize + 1

!print*, 'nloop  ',nloop
!print*, 'dimFO',dimFO

allocate(work1(dimFO,MaxBatchSize),ints(NBasis,NBasis))

off = 0
k   = 0
l   = 1
! exchange loop (FO|FO)
do iloop=1,nloop

   ! batch size for each iloop; last one is smaller
   BatchSize = min(MaxBatchSize,dimFO-off)

   ! assemble (FO|BatchSize) batch from CholVecs
   call dgemm('T','N',dimFO,BatchSize,NCholesky,1d0,MatFF,NCholesky, &
              MatFF(:,off+1:off+BatchSize),NCholesky,0d0,work1,dimFO)

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NBasis) then
         k = 1
         l = l + 1
      endif

      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0
      !kl = (l-1)*NBasis + k
      !print*, 'k,l,kl = ',k,l,kl
      !print*, 'ints',norm2(ints)

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = HNOCoef*Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif

         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! exchange
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

      ! AUXILIARY MATRIX WMAT
      if(l<=NOccup) then
         do ir=1,NOccup
            val = 0
            val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                 sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                 sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                 sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                 sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
            WMAT(k,ir) = WMAT(k,ir) + val
         enddo
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then

            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k

         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo

   off = off + MaxBatchSize

enddo

deallocate(work1)

!Print*, 'from JK_Chol_loop'
!Print*, 'ABPLUS-after K',norm2(ABPLUS)
!Print*, 'ABMIN -after K',norm2(ABMIN)
!ABPLUS = 0
!ABMIN  = 0
!HNO = HNO_save
!WMAT = 0
!AuxI = 0
!AuxIO = 0

dimOO = NOccup*NOccup
nloop = (dimOO - 1) / MaxBatchSize + 1

allocate(work1(NBasis**2,MaxBatchSize))
allocate(work2(NCholesky,MaxBatchSize))

off = 0
k   = 0
l   = 1
! Coulomb loop (FF|OO)
do iloop=1,nloop

   ! batch size
   BatchSize = min(MaxBatchSize,dimOO-off)

   kk = k
   ll = l
   do iBatch=1,BatchSize
      kk = kk + 1
      if(kk>NOccup) then
         kk = 1
         ll = ll + 1
      endif
      kl = (ll - 1)*NBasis + kk
      work2(:,iBatch) = MatFF(:,kl)
   enddo

   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,1d0,MatFF,NCholesky, &
              work2,NCholesky,0d0,work1,NBasis**2)

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NOccup) then
         k = 1
         l = l + 1
      endif

      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(k>NOccup.or.l>NOccup) cycle
      kl = (l - 1)*NOccup + k
      !print*, 'k,l,kl',k,l,kl
      !print*, norm2(ints)

      ! COMPUTE THE ENERGY FOR CHECKING
      if(present(ETot)) then
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*HNOCoef*Occ(k)
         if(k>INActive) then
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo

      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1

      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               endif
            enddo

         endif
      endif

      ! T1
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then

         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then

                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l

         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo

   off = off + MaxBatchSize

enddo

!Print*, 'ABPLUS-after J',norm2(ABPLUS)
!Print*, 'ABMIN -after J',norm2(ABMIN)

!ABPLUS = 0
!ABMIN  = 0
!HNO = HNO_save
!WMAT = 0
!AuxI = 0
!AuxIO = 0

!print*, 'HNO   on exit  =',norm2(HNO)
!print*, 'WMAT  on exit  =',norm2(WMAT)
!print*, '2*WMA on exit  =',2d0*norm2(WMAT)
!print*, 'AuxI  on exit  =',norm2(AuxI)
!print*, 'AuxIO on exit  =',norm2(AuxIO)

deallocate(work2,work1,ints)

end subroutine JK_Chol_loop

subroutine JK_SR_Chol_loop(ABPLUS,ABMIN,HNO,AuxI,AuxIO,WMAT,RDM2val,Occ,AuxCoeff,IGem,AuxInd,pos,&
                        INActive,NOccup,NDim,NDimX,NBasis,NInte1,IntJFile,IntKFile,ACAlpha,AB,ETot)

implicit none

integer,intent(in) :: NDim,NDimX,NBasis,NInte1
integer,intent(in) :: INActive,NOccup,AB
integer,intent(in) :: IGem(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision,intent(inout) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(inout) :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision,intent(in) :: ACAlpha
double precision,intent(in) :: AuxCoeff(3,3,3,3),Occ(NBasis),RDM2val(NOccup,NOccup,NOccup,NOccup)
character(*) :: IntJFile,IntKFile
double precision,intent(inout),optional :: ETot

integer :: i,j,k,l,kk,ll,kl
integer :: ip,iq,ir,is,iu,ipq,irs
integer :: iunit
integer :: iloop,nloop,off,dimFO,dimOO,dimFF
integer :: mloop
integer :: iBatch
integer :: BatchSize,MaxBatchSize = 120
integer :: NCholesky,NCholErf
double precision :: val
double precision :: AuxVal,HNOCoef
double precision,allocatable :: FF(:,:),FFErf(:,:)
double precision,allocatable :: FFTr(:,:),FFErfTr(:,:)
double precision,allocatable :: TmpTr(:,:),TmpErfTr(:,:)
double precision,allocatable :: ints(:,:),intsFR(:,:)
double precision,allocatable :: work1(:,:),work2(:,:)
double precision,allocatable :: work3(:,:),work4(:,:)
double precision,allocatable :: work6(:,:),work8(:,:)
double precision,allocatable :: work10(:,:),work12(:,:)
!
! compute Hessian matrices with modified integrals : 
! a) assemble all required integrals
! b) symmetrize : g^SR = 1/4 * ( <p*q|rs>_SR + <pq*|rs>_SR + <pq|r*s>_SR + <pq|rs*>_SR )
! c) use : <pq|rs>_mod =  <pq|rs>^full-range + g^SR
!

if(AB==1) then
   HNOCoef = -1
elseif(AB==0) then
   HNOCoef = 1 - ACAlpha
elseif(AB==2) then
   HNOCoef = 1
endif

! read regular cholesky (k|r|FF) vecs
open(newunit=iunit,file='cholvecs',form='unformatted')
read(iunit) NCholesky
allocate(FF(NCholesky,NBasis**2))
read(iunit) FF
close(iunit)

! read transformed cholesky (k|r|FF*) vecs
open(newunit=iunit,file='chol1vFR',form='unformatted')
read(iunit) NCholesky
allocate(FFTr(NCholesky,NBasis**2))
read(iunit) FFTr
close(iunit)

! read LR cholesky (k|erf|FF) vecs
open(newunit=iunit,file='cholvErf',form='unformatted')
read(iunit) NCholErf
allocate(FFErf(NCholErf,NBasis**2))
read(iunit) FFErf
close(iunit)

! read transformed LR cholesky (k|erf|FF*) vecs
open(newunit=iunit,file='chol1vLR',form='unformatted')
read(iunit) NCholErf
allocate(FFErfTr(NCholErf,NBasis**2))
read(iunit) FFErfTr
close(iunit)

allocate(TmpErfTr(NCholErf,NBasis*NOccup))
! p*q : LR
do iq=1,NOccup
   do ip=1,NBasis
      TmpErfTr(:,ip+(iq-1)*NBasis)=FFErfTr(:,iq+(ip-1)*NBasis)
   enddo
enddo
allocate(TmpTr(NCholesky,NBasis*NOccup))
! p*q : FR
do iq=1,NOccup
   do ip=1,NBasis
      TmpTr(:,ip+(iq-1)*NBasis)=FFTr(:,iq+(ip-1)*NBasis)
   enddo
enddo

! set number of loops over integrals
dimFO = NBasis*NOccup
nloop = (dimFO - 1) / MaxBatchSize + 1

allocate(work1(dimFO,MaxBatchSize),work2(dimFO,MaxBatchSize))
allocate(ints(NBasis,NBasis))
allocate(intsFR(NBasis,NBasis))

off = 0
k   = 0
l   = 1

  kk = 0
  ll = 1
! exchange loop (FO|FO)
!print*,'exchange loop (FO|FO)'
do iloop=1,nloop

   ! batch size for each iloop; last one is smaller
   BatchSize = min(MaxBatchSize,dimFO-off)

   ! (pq*|rs) : SR=-LR+FR
   call dgemm('T','N',dimFO,BatchSize,NCholErf,-0.25d0,FFErfTr,NCholErf, &
              FFErf(:,off+1:off+BatchSize),NCholErf,0d0,work1,dimFO)
   call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,FFTr,NCholesky, &
              FF(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
   ! (pq|rs*): SR=-LR+FR
   call dgemm('T','N',dimFO,BatchSize,NCholErf,-0.25d0,FFErf,NCholErf, &
              FFErfTr(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
   call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,FF,NCholesky, &
              FFTr(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
   ! (p*q|rs)
   call dgemm('T','N',dimFO,BatchSize,NCholErf,-0.25d0,TmpErfTr,NCholErf, &
              FFErf(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
   call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,TmpTr,NCholesky, &
              FF(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)
   ! (pq|r*s)
   call dgemm('T','N',dimFO,BatchSize,NCholErf,-0.25d0,FFErf,NCholErf, &
              TmpErfTr(:,off+1:off+BatchSize),NCholErf,1d0,work1,dimFO)
   call dgemm('T','N',dimFO,BatchSize,NCholesky,0.25d0,FF,NCholesky, &
              TmpTr(:,off+1:off+BatchSize),NCholesky,1d0,work1,dimFO)

   !! test: SR bez transformacji
   !! full-range
   !call dgemm('T','N',dimFO,BatchSize,NCholesky,1d0,FF,NCholesky, &
   !           FF(:,off+1:BatchSize),NCholesky,0d0,work1,dimFO)
   !! long-range
   !call dgemm('T','N',dimFO,BatchSize,NCholErf,-1d0,FFErf,NCholErf, &
   !           FFErf(:,off+1:BatchSize),NCholErf,1d0,work1,dimFO)
   !!! work1=short-range

   ! regular
   call dgemm('T','N',dimFO,BatchSize,NCholesky,1d0,FF,NCholesky, &
              FF(:,off+1:off+BatchSize),NCholesky,0d0,work2,dimFO)
   work1 = work1 + work2
   ! work1 = regular + short-range* 

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NBasis) then
         k = 1
         l = l + 1
      endif

      do j=1,NOccup
         do i=1,NBasis
            ints(i,j)   = work1((j-1)*NBasis+i,iBatch)
            intsFR(i,j) = work2((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0
      intsFR(:,NOccup+1:NBasis) = 0
      !kl = (l-1)*NBasis + k
      !print*, 'k,l,kl = ',k,l,kl
      !print*, 'ints',norm2(ints)

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = HNOCoef*Occ(l)
         if(IGem(l)==1) then
            if(IGem(k)==2) then
               do i=INActive+1,NOccup
                  HNO(i,k) = HNO(i,k) - val*intsFR(i,l)
                  !HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*intsFR(i,l)
                  !HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif

         if(IGem(l)==2) then
            if(IGem(k)==1) then
               do i=1,INActive
                  HNO(i,k) = HNO(i,k) - val*intsFR(i,l)
                  !HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            elseif(IGem(k)==3) then
               do i=NOccup+1,NBasis
                  HNO(i,k) = HNO(i,k) - val*intsFR(i,l)
                  !HNO(i,k) = HNO(i,k) - val*ints(i,l)
               enddo
            endif
         endif
      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! exchange
      if(l<=INActive) then
         do i=1,NBasis
            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif
      if(l<=NOccup) then
         do i=1,NBasis
            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
         enddo
      endif

      ! AUXILIARY MATRIX WMAT
      if(l<=NOccup) then
         do ir=1,NOccup
            val = 0
            val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
                 sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
                 sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
                 sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
            val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
                 sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
            WMAT(k,ir) = WMAT(k,ir) + val
         enddo
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
      ! Coulomb
      if(k>INActive.and.l<=NOccup) then
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then

            do iq=1,NOccup
               do ip=INActive+1,NBasis
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                     val = 0
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     val = 0
                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                     val = 2*AuxVal*val*ints(ip,iq)

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            enddo

         endif
      endif

      ! exchange
      if(l<=NOccup) then
         if(k>INActive) then
            ip = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
                        val = - AuxVal*val*ints(ir,iq)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               endif
            enddo
         endif
      endif

      ! T1 p->q
      if((k>INActive).and.(l<=NOccup)) then
         ip = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if((l<=NOccup).and.(k>INActive)) then
         iq = l
         ir = k

         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then
                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                  endif
               enddo
            endif
         enddo

      endif

      ! T4+T6: AuxInd = 2
      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
         ir = k
         iu = l

         ! T4
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T6
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do is=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo

   off = off + MaxBatchSize

enddo

deallocate(work1,work2)
deallocate(TmpErfTr,TmpTr)

!Print*, 'from JK_SR_Chol_loop'
!Print*, 'ABPLUS-after K',norm2(ABPLUS)
!Print*, 'ABMIN -after K',norm2(ABMIN)

dimOO = NOccup*NOccup
nloop = (dimOO - 1) / MaxBatchSize + 1

allocate(TmpErfTr(NCholErf,NBasis*NBasis))
! p*q : LR
do iq=1,NBasis
   do ip=1,NBasis
      TmpErfTr(:,ip+(iq-1)*NBasis)=FFErfTr(:,iq+(ip-1)*NBasis)
   enddo
enddo
allocate(TmpTr(NCholesky,NBasis*NBasis))
! p*q : FR
do iq=1,NBasis
   do ip=1,NBasis
      TmpTr(:,ip+(iq-1)*NBasis)=FFTr(:,iq+(ip-1)*NBasis)
   enddo
enddo

allocate(work1(NBasis**2,MaxBatchSize))
allocate(work3(NBasis**2,MaxBatchSize))
!
allocate(work2(NCholesky,MaxBatchSize))
allocate(work4(NCholErf,MaxBatchSize))
!
allocate(work6(NCholesky,MaxBatchSize))
allocate(work8(NCholErf,MaxBatchSize))
!
allocate(work10(NCholesky,MaxBatchSize))
allocate(work12(NCholErf,MaxBatchSize))

off = 0
k   = 0
l   = 1
! Coulomb loop (FF|OO)
!print*, 'Coulomb loop (FF|OO)...'
do iloop=1,nloop

   ! batch size
   BatchSize = min(MaxBatchSize,dimOO-off)

   kk = k
   ll = l
   do iBatch=1,BatchSize
      kk = kk + 1
      if(kk>NOccup) then
         kk = 1
         ll = ll + 1
      endif
      kl = (ll - 1)*NBasis + kk
      work2(:,iBatch) = FF(:,kl)
      work4(:,iBatch) = FFErf(:,kl)
!
      work6(:,iBatch) = FFTr(:,kl)
      work8(:,iBatch) = FFErfTr(:,kl)
!
      work10(:,iBatch) = TmpTr(:,kl)
      work12(:,iBatch) = TmpErfTr(:,kl)
!
   enddo

   ! (pq*|rs) : SR=-LR+FR
   call dgemm('T','N',NBasis**2,BatchSize,NCholErf,-0.25d0,FFErfTr,NCholErf, &
              work4,NCholErf,0d0,work1,NBasis**2)
   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,0.25d0,FFtr,NCholesky, &
              work2,NCholesky,1d0,work1,NBasis**2)
   ! (pq|rs*) : SR=-LR+FR
   call dgemm('T','N',NBasis**2,BatchSize,NCholErf,-0.25d0,FFErf,NCholErf, &
              work8,NCholErf,1d0,work1,NBasis**2)
   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,0.25d0,FF,NCholesky, &
              work6,NCholesky,1d0,work1,NBasis**2)
   ! (p*q|rs) : SR=-LR+FR
   call dgemm('T','N',NBasis**2,BatchSize,NCholErf,-0.25d0,TmpErfTr,NCholErf, &
              work4,NCholErf,1d0,work1,NBasis**2)
   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,0.25d0,TmpTr,NCholesky, &
              work2,NCholesky,1d0,work1,NBasis**2)
   ! (pq|r*s) : SR=-LR+FR
   call dgemm('T','N',NBasis**2,BatchSize,NCholErf,-0.25d0,FFErf,NCholErf, &
              work12,NCholErf,1d0,work1,NBasis**2)
   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,0.25d0,FF,NCholesky, &
              work10,NCholesky,1d0,work1,NBasis**2)

   !!! test: SR bez transformacji
   !! full-range
   !call dgemm('T','N',NBasis**2,BatchSize,NCholesky,1d0,FF,NCholesky, &
   !           work2,NCholesky,0d0,work1,NBasis**2)
   !! long-range
   !call dgemm('T','N',NBasis**2,BatchSize,NCholErf,-1d0,FFErf,NCholErf, &
   !           work4,NCholErf,1d0,work1,NBasis**2)
   !!! work1 = short-range

   ! regular
   call dgemm('T','N',NBasis**2,BatchSize,NCholesky,1d0,FF,NCholesky, &
              work2,NCholesky,0d0,work3,NBasis**2)
   work1 = work1 + work3
   !! work1 = regular + short-range* 

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NOccup) then
         k = 1
         l = l + 1
      endif

      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i,iBatch)
            intsFR(i,j) = work3((j-1)*NBasis+i,iBatch)
         enddo
      enddo

      if(k>NOccup.or.l>NOccup) cycle
      kl = (l - 1)*NOccup + k
      !print*, 'k,l,kl',k,l,kl
      !print*, norm2(ints)

      ! COMPUTE THE ENERGY FOR CHECKING
      if(present(ETot)) then
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
      endif

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*HNOCoef*Occ(k)
         if(k>INActive) then
            do j=1,INActive
               do i=1,INActive
                  HNO(i,j) = HNO(i,j) + val*intsFR(i,j)
                  !HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         if(k<=INActive) then
            do j=INActive+1,NOccup
               do i=INActive+1,NOccup
                  HNO(i,j) = HNO(i,j) + val*intsFR(i,j)
                  !HNO(i,j) = HNO(i,j) + val*ints(i,j)
               enddo
            enddo
         endif
         do j=NOccup+1,NBasis
            do i=NOccup+1,NBasis
               HNO(i,j) = HNO(i,j) + val*intsFR(i,j)
               !HNO(i,j) = HNO(i,j) + val*ints(i,j)
            enddo
         enddo

      endif

      ! AUXILIARY MATRIX AuxI AND AuxIO

      ! Coulomb
      if(k==l) then
         val = 2*Occ(k)
         if(k<=INActive) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
               enddo
            enddo
         endif
         if(k>INActive.and.k<=NOccup) then
            do j=1,NBasis
               do i=1,NBasis
                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
               enddo
            enddo
         endif
      endif

      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
      ! T1+T2+T3+T4+T5+T6: AuxInd = 1

      ! exchange
      if(l<=NOccup) then
         if(k<=NOccup) then
            iq = k
            is = l

            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then

                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))

                        val = 0
                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
                        val = - AuxVal*val*ints(ip,ir)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               endif
            enddo

         endif
      endif

      ! T1
      if((k<=NOccup).and.(l<=NOccup)) then
         iq = k
         is = l

         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do ip=INActive+1,NOccup
                  ipq = pos(ip,iq)
                  if(ipq>0) then

                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))

                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                  endif
               enddo
            endif
         enddo

      endif

      ! T3
      if(IGem(k)==2.and.IGem(l)==2) then

         do is=INActive+1,NOccup
            do ir=INActive+1,NBasis
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     do ip=INActive+1,NBasis
                        ipq = pos(ip,iq)
                        if(ipq>0) then

                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                                RDM2val(iq,is,k,l)*ints(ip,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      endif

      ! T2+T5: AuxInd = 2
      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
         is = k
         iu = l

         ! T2
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

         ! T5
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=INActive+1,NOccup
                  do ip=INActive+1,NBasis
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

                     endif
                  enddo
               enddo
            endif
         enddo
         ! p->q
         do ir=INActive+1,NOccup
            irs = pos(ir,is)
            if(irs>0) then
               do iq=1,NOccup
                  do ip=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo
            endif
         enddo

      endif

   enddo

   off = off + MaxBatchSize

enddo

!deallocate(work4,work3,work2,work1)
deallocate(work4,Work2,work1)
deallocate(work6,work8)
deallocate(ints)
deallocate(intsFR)

!Print*, 'from JK_SR_Chol_loop'
!Print*, 'ABPLUS-after JK',norm2(ABPLUS)
!Print*, 'ABMIN -after JK',norm2(ABMIN)

deallocate(FF,FFErf)
deallocate(FFTr,FFErfTr)

!print*, 'HNO   on exit  =',norm2(HNO)
!print*, 'WMAT  on exit  =',norm2(WMAT)
!print*, '2*WMA on exit  =',2d0*norm2(WMAT)
!print*, 'AuxI  on exit  =',norm2(AuxI)
!print*, 'AuxIO on exit  =',norm2(AuxIO)

end subroutine JK_SR_Chol_loop

subroutine Chol_AC0ECorr(ECorr,ABPLUS,IndN,IndX,Occ,INActive,NOccup,NDimX,NBasis,CholVecFile)
implicit none

integer,intent(in)          :: INActive,NOccup,NDimX,NBasis
integer,intent(in)          :: IndN(2,NDimX),IndX(NDimX)
double precision,intent(in) :: Occ(NBasis)
double precision,intent(in) :: ABPLUS(NDimX,NDimX)
character(*),intent(in)     :: CholVecFile
double precision,intent(out) :: ECorr

integer :: iunit
integer :: i,j,k,l
integer :: ip,iq,ipq,ir,is,irs
integer :: NCholesky
integer :: dimFO,iBatch,BatchSize
integer :: iloop,nloop,off
integer :: MaxBatchSize = 100
integer :: IGem(NBasis),AuxCoeff(3,3,3,3)
integer :: pos(NBasis,NBasis)
double precision :: EAll, EIntra
double precision :: C(NBasis)
double precision :: Cpq, Crs, Aux
double precision,allocatable :: MatFF(:,:)
double precision,allocatable :: ints(:,:),work(:,:)

!print*, 'Entering Chol_AC0ECor...'

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

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

do i=1,NBasis
   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

! read cholesky (FF|K) vectors
open(newunit=iunit,file=CholVecFile,form='unformatted')
read(iunit) NCholesky
allocate(MatFF(NCholesky,NBasis**2))
read(iunit) MatFF
close(iunit)

! set number of loops over integrals
dimFO = NOccup*NBasis
nloop = (dimFO - 1) / MaxBatchSize + 1

allocate(work(dimFO,MaxBatchSize),ints(NBasis,NBasis))

ECorr  = 0d0
EAll   = 0d0
EIntra = 0d0

off = 0
k   = 0
l   = 1
! exchange loop (FO|FO)
do iloop=1,nloop

   ! batch size for each iloop; last one is smaller
   BatchSize = min(MaxBatchSize,dimFO-off)

   ! assemble (FO|BatchSize) batch from CholVecs
   call dgemm('T','N',dimFO,BatchSize,NCholesky,1d0,MatFF,NCholesky, &
              MatFF(:,off+1:BatchSize),NCholesky,0d0,work,dimFO)

   ! loop over integrals
   do iBatch=1,BatchSize

      k = k + 1
      if(k>NBasis) then
         k = 1
         l = l + 1
      endif

      if(pos(k,l)==0) cycle
       ir = k
       is = l
       irs = pos(k,l)

       do j=1,NOccup
          do i=1,NBasis
             ints(i,j) = work((j-1)*NBasis+i,iBatch)
          enddo
       enddo

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

       do j=1,NBasis
          do i=1,j
             if(pos(j,i)/=0) then
               ipq = pos(j,i)
               ip = j
               iq = i
               Crs = C(ir)+C(is)
               Cpq = C(ip)+C(iq)

               Aux = Crs*Cpq*ABPLUS(ipq,irs)
               EAll = EAll + Aux*ints(j,i)

               if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)

             endif
          enddo
       enddo

   enddo ! iBatch

   off = off + MaxBatchSize

enddo

print*, 'EAll-Chol  ', EAll
print*, 'EIntra-Chol', EIntra

ECorr = EAll-EIntra

deallocate(MatFF)
deallocate(work,ints)

end subroutine Chol_AC0ECorr

end module
