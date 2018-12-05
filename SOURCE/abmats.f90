module abmat
use types,only : LOUT
use tran
implicit none

contains 

subroutine AB_CAS_mithap(ABPLUS,ABMIN,ETot,URe,Occ,XOne, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1,IntFileName,ACAlpha,&
     AB1)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
! 
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! THE FOLLOWING SYMMETRY IS ASSUMED
! RDM2(ij,kl) = RDM2(kl,ij)
! ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
! SIZE: NRDM2 = NBasis**2*(NBasis**2+1)/2
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntFileName
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(in)  :: ACAlpha
logical,intent(in) :: AB1

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,HNOCoef,val
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
ETot   = 0

! set dimensions
NOccup = NAct + INActive
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
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)


call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

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

open(newunit=iunit,file=trim(IntFileName),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do ll=1,NBasis
   do kk=1,ll
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
      do klround=1,merge(1,2,kk==ll)
         select case(klround)
         case(1)
            k = kk
            l = ll
         case(2)
            k = ll
            l = kk
         end select

! COMPUTE THE ENERGY FOR CHECKING
         if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))

! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
         ! Coulomb
         if(IGem(k)==IGem(l)) then

            val = 0
            do i=1,NBasis
               if(IGem(i)/=IGem(k)) val = val + Occ(i)*ints(i,i)
            enddo
            val = 2*HNOCoef*val
            HNO(k,l) = HNO(k,l) + val
            
         ! exchange
         else

            val = HNOCoef*Occ(k)
            do i=1,NBasis
               if(IGem(i)==IGem(l)) HNO(i,l) = HNO(i,l) - val*ints(i,k)
            enddo

         endif

! AUXILIARY MATRIX AuxI AND AuxIO
         ! Coulomb
         val = 0
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),1,1)
         do i=1,INActive
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxIO(k,l) = AuxIO(k,l) + 2*val
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),2,2)
         do i=INActive+1,NOccup
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxI(k,l) = AuxI(k,l) + 2*val

         ! exchange
         if(k<=INActive) then
            do i=1,NBasis
               AuxIO(i,l) = AuxIO(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
            enddo
         endif
         if(k<=NOccup) then
            do i=1,NBasis
               AuxI(i,l) = AuxI(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
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
                        val = 2*AuxVal*val*ints(iq,ip)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo

            endif
         endif

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
                           val = - AuxVal*val*ints(iq,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
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
         if((k>INActive).and.(l>INActive)) then
            ip = k
            ir = l

            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
            
         endif

         if((k<=NOccup).and.(l>INActive)) then
            iq = k
            ir = l

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
                                sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(INActive+1:NOccup,ip))

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
                                sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(INActive+1:NOccup,ip))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo

         endif
         
      enddo
   enddo
enddo

close(iunit)

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
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo

   endif
enddo

!print*, "AB-my",norm2(ABPLUS),norm2(ABMIN)

!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$HNO=transpose(HNO)
!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!!$
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxI=transpose(AuxI)
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxIO=transpose(AuxIO)
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_mithap

subroutine AB_CAS_FOFO(ABPLUS,ABMIN,ETot,URe,Occ,XOne, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     IntJFile,IntKFile,ACAlpha,AB1)
!
! COMPUTE THE A+B AND A-B MATRICES FOR 2-RDM READ FROM A rdm2.dat FILE
! 
! RDM2 IS IN NO REPRESENTATION. IT IS PARTIALLY SPIN-SUMMED
! THE FOLLOWING SYMMETRY IS ASSUMED
! RDM2(ij,kl) = RDM2(kl,ij)
! ONLY ELEMENTS ij >= kl ARE STORED BUT ONLY FOR THE ACTIVE ELEMENTS
! COULOMB INTEGRALS ARE READ FROM IntJFile IN (FF|OO) FORMAT
! EXCHANGE INTEGRALS ARE READ FROM IntKFile IN (FO|FO) FORMAT
 implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDim,NDim),ABMIN(NDim,NDim)
double precision,intent(out) :: ETot
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis)
double precision,intent(in)  :: ACAlpha
logical,intent(in) :: AB1

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit1,iunit2,ios
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,HNOCoef,val
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2

ABPLUS = 0
ABMIN  = 0
ETot   = 0

! set dimensions
NOccup = NAct + INActive
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
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit1,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit1,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit1)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)


call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

val = 0
do i=1,NOccup
   val = val + Occ(i)*HNO(i,i)
enddo
ETot = ETot + 2*val

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

! does not work in Y01CAS?
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           HNOCoef,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      !call triang_to_sq2(work1,ints,NBasis)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

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
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      ! COMPUTE THE ENERGY FOR CHECKING
      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))

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
enddo

close(iunit2)

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
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val

      enddo

   endif
enddo

!print*, "AB-my",norm2(ABPLUS),norm2(ABMIN)

!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$HNO=transpose(HNO)
!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot
!!$
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxI=transpose(AuxI)
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxIO=transpose(AuxIO)
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_FOFO  

!subroutine AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           HNOCoef,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)
!!  ASSEMBLES AB MATRICES FOR CAS USING FOFO/FFOO
!!  FILES. SHARED BY AB_CAS_FOFO AND Y01CAS_FOFO
!implicit none
!
!integer,intent(in) :: NDimX,NDim,INActive,NOccup,NBasis
!character(*) :: IntJFile,IntKFile
!double precision,intent(inout) :: ABPLUS(:,:),ABMIN(:,:),HNO(NBasis,NBasis)
!double precision,intent(in)  :: Occ(NBasis),AuxCoeff(3,3,3,3)
!double precision,intent(in)  :: HNOCoef
!double precision,intent(in)  :: RDM2val(NOccup,NOccup,NOccup,NOccup)
!integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGem(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
!
!double precision :: AuxVal,val,ETot
!integer :: iunit1,iunit2 
!integer :: i,j,k,l,ij,kl,kk,ll,klround
!integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
!double precision,allocatable :: AuxI(:,:),AuxIO(:,:),WMAT(:,:),C(:)
!double precision,allocatable :: work1(:),work2(:)
!double precision,allocatable :: ints(:,:)
!
!allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))
!allocate(WMAT(NBasis,NBasis),AuxIO(NBasis,NBasis),AuxI(NBasis,NBasis),C(NBasis))
!
!AuxI = 0
!AuxIO = 0
!WMAT = 0
!
!open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBasis*NOccup)
!
!! exchange loop (FO|FO)
!kl = 0
!do ll=1,NOccup
!   do kk=1,NBasis
!      kl = kl + 1
!      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
!      do j=1,NOccup
!         do i=1,NBasis
!            ints(i,j) = work1((j-1)*NBasis+i)
!         enddo
!      enddo
!      k = kk
!      l = ll
!
!      if(l>NOccup) cycle
!      ints(:,NOccup+1:NBasis) = 0
!
!      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
!
!      ! exchange
!      if(l<=NOccup) then
!         val = HNOCoef*Occ(l)
!         if(IGem(l)==1) then
!            if(IGem(k)==2) then
!               do i=INActive+1,NOccup
!                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
!               enddo
!            elseif(IGem(k)==3) then
!               do i=NOccup+1,NBasis
!                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
!               enddo
!            endif
!         endif
!         
!         if(IGem(l)==2) then
!            if(IGem(k)==1) then
!               do i=1,INActive
!                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
!               enddo
!            elseif(IGem(k)==3) then
!               do i=NOccup+1,NBasis
!                  HNO(i,k) = HNO(i,k) - val*ints(i,l)
!               enddo
!            endif
!         endif
!      endif
!         
!      ! AUXILIARY MATRIX AuxI AND AuxIO
!      
!      ! exchange 
!      if(l<=INActive) then
!         do i=1,NBasis
!            AuxIO(i,k) = AuxIO(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
!         enddo
!      endif
!      if(l<=NOccup) then
!         do i=1,NBasis
!            AuxI(i,k) = AuxI(i,k) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(l)*ints(i,l)
!         enddo
!      endif
!
!      ! AUXILIARY MATRIX WMAT
!      if(l<=NOccup) then
!         do ir=1,NOccup
!            val = 0
!            val = val + AuxCoeff(IGem(k),IGem(l),1,1)* &
!                 sum(ints(1:INActive,1:INActive)*RDM2val(1:INActive,1:INActive,ir,l))
!            val = val + AuxCoeff(IGem(k),IGem(l),2,1)* &
!                 sum(ints(INActive+1:NOccup,1:INActive)*RDM2val(INActive+1:NOccup,1:INActive,ir,l))
!            val = val + AuxCoeff(IGem(k),IGem(l),1,2)* &
!                 sum(ints(1:INActive,INActive+1:NOccup)*RDM2val(1:INActive,INActive+1:NOccup,ir,l))
!            val = val + AuxCoeff(IGem(k),IGem(l),2,2)* &
!                 sum(ints(INActive+1:NOccup,INActive+1:NOccup)*RDM2val(INActive+1:NOccup,INActive+1:NOccup,ir,l))
!            WMAT(k,ir) = WMAT(k,ir) + val
!         enddo
!      endif
!
!      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
!      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
!      ! Coulomb
!      if(k>INActive.and.l<=NOccup) then
!         ir = k
!         is = l
!         irs = pos(ir,is)
!         if(irs>0) then
!            
!            do iq=1,NOccup
!               do ip=INActive+1,NBasis
!                  ipq = pos(ip,iq)
!                  if(ipq>0) then
!
!                     AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
!                     
!                     val = 0
!                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
!                     if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
!                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
!                     if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
!                     val = 2*AuxVal*val*ints(ip,iq)
!                  
!                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!
!                     val = 0
!                     if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
!                     if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
!                     if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
!                     if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
!                     val = 2*AuxVal*val*ints(ip,iq)
!                     
!                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                     
!                  endif
!               enddo
!            enddo
!
!         endif
!      endif
!
!      ! exchange
!      if(l<=NOccup) then
!         if(k>INActive) then
!            ip = k
!            is = l
!
!            do ir=INActive+1,NBasis
!               irs = pos(ir,is)
!               if(irs>0) then                  
!                  do iq=1,NOccup
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        
!                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
!                        
!                        val = 0
!                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val + Occ(iq)*Occ(ir)
!                        if(AuxInd(IGem(ip),IGem(is))==1) val = val + Occ(ip)*Occ(is)
!                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val - Occ(ip)*Occ(ir)
!                        if(AuxInd(IGem(iq),IGem(is))==1) val = val - Occ(iq)*Occ(is)
!                        val = - AuxVal*val*ints(ir,iq)
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                        
!                     endif
!                  enddo
!               endif
!            enddo
!         endif
!      endif
!
!      ! T1 p->q
!      if((k>INActive).and.(l<=NOccup)) then
!         ip = k
!         is = l
!         
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=INActive+1,NOccup
!                  ipq = pos(ip,iq)
!                  if(ipq>0) then
!                     
!                     val = AuxCoeff(IGem(ip),IGem(is),2,2)* &
!                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
!                     
!                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                     
!                  endif
!               enddo
!            endif
!         enddo
!         
!      endif
!
!      ! T3
!      if((l<=NOccup).and.(k>INActive)) then
!         iq = l
!         ir = k
!         
!         do is=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do ip=INActive+1,NOccup
!                  ipq = pos(ip,iq)
!                  if(ipq>0) then
!                     val = AuxCoeff(IGem(iq),IGem(ir),2,2)* &
!                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,is)*ints(INActive+1:NOccup,INActive+1:NOccup))
!                     
!                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                     ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                     
!                  endif
!               enddo
!            endif
!         enddo
!         
!      endif
!
!      ! T4+T6: AuxInd = 2
!      if((k>INActive).and.(l>INActive.and.l<=NOccup)) then
!         ir = k
!         iu = l
!
!         ! T4
!         do is=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=INActive+1,NOccup
!                  do ip=INActive+1,NBasis
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(ip,INActive+1:NOccup))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                           
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         ! p->q
!         do is=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=1,NOccup
!                  do ip=INActive+1,NOccup
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,ip,is,iu)*ints(INActive+1:NOccup,iq))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         
!         ! T6
!         do is=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=1,NOccup
!                  do ip=INActive+1,NOccup
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then                        
!                        val = - AuxCoeff(IGem(iq),IGem(ir),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,ip,iu,is)*ints(INActive+1:NOccup,iq))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         ! p->q
!         do is=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=INActive+1,NOccup
!                  do ip=INActive+1,NBasis
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then  
!                        val = - AuxCoeff(IGem(ip),IGem(ir),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(ip,INActive+1:NOccup))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         
!      endif
!      
!   enddo
!enddo
!
!close(iunit1)
!
!open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
!     access='DIRECT',recl=8*NBasis**2)
!
!! Coulomb loop (FF|OO)
!kl = 0
!do ll=1,NOccup
!   do kk=1,NOccup
!      kl = kl + 1
!      read(iunit2,rec=kl) work1(1:NBasis**2)
!      do j=1,NBasis
!         do i=1,NBasis
!            ints(i,j) = work1((j-1)*NBasis+i)
!         enddo
!      enddo
!      
!      k = kk
!      l = ll
!
!      if(k>NOccup.or.l>NOccup) cycle
!
!      ! COMPUTE THE ENERGY FOR CHECKING
!      if((k<=NOccup).and.(l<=NOccup)) ETot = ETot + sum(RDM2val(:,:,k,l)*ints(1:NOccup,1:NOccup))
!
!      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
!
!      ! Coulomb
!      if(k==l.and.k<=NOccup) then
!         val = 2*HNOCoef*Occ(k)
!         if(k>INActive) then 
!            do j=1,INActive
!               do i=1,INActive
!                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
!               enddo
!            enddo
!         endif
!         if(k<=INActive) then
!            do j=INActive+1,NOccup
!               do i=INActive+1,NOccup
!                  HNO(i,j) = HNO(i,j) + val*ints(i,j)
!               enddo
!            enddo
!         endif
!         do j=NOccup+1,NBasis
!            do i=NOccup+1,NBasis
!               HNO(i,j) = HNO(i,j) + val*ints(i,j)
!            enddo
!         enddo
!         
!      endif
!          
!      ! AUXILIARY MATRIX AuxI AND AuxIO
!
!      ! Coulomb
!      if(k==l) then
!         val = 2*Occ(k)
!         if(k<=INActive) then
!            do j=1,NBasis
!               do i=1,NBasis
!                  AuxIO(i,j) = AuxIO(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
!                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),1,1)
!               enddo
!            enddo
!         endif
!         if(k>INActive.and.k<=NOccup) then
!            do j=1,NBasis
!               do i=1,NBasis
!                  AuxI(i,j) = AuxI(i,j) + val*ints(i,j)*AuxCoeff(IGem(i),IGem(j),2,2)
!               enddo
!            enddo
!         endif
!      endif
!
!      ! CONSTRUCT TWO-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
!      ! T1+T2+T3+T4+T5+T6: AuxInd = 1
!      
!      ! exchange
!      if(l<=NOccup) then
!         if(k<=NOccup) then
!            iq = k
!            is = l
!            
!            do ir=INActive+1,NBasis
!               irs = pos(ir,is)
!               if(irs>0) then
!                  do ip=INActive+1,NBasis
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        
!                        AuxVal = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))
!                        
!                        val = 0
!                        if(AuxInd(IGem(ip),IGem(ir))==1) val = val + Occ(ip)*Occ(ir)
!                        if(AuxInd(IGem(iq),IGem(is))==1) val = val + Occ(iq)*Occ(is)
!                        if(AuxInd(IGem(iq),IGem(ir))==1) val = val - Occ(iq)*Occ(ir)
!                        if(AuxInd(IGem(ip),IGem(is))==1) val = val - Occ(ip)*Occ(is)
!                        val = - AuxVal*val*ints(ip,ir)
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                        
!                     endif
!                  enddo
!               endif
!            enddo
!               
!         endif
!      endif
!
!      ! T1 
!      if((k<=NOccup).and.(l<=NOccup)) then
!         iq = k
!         is = l
!         
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do ip=INActive+1,NOccup
!                  ipq = pos(ip,iq)
!                  if(ipq>0) then
!                     
!                     val = AuxCoeff(IGem(iq),IGem(is),2,2)* &
!                          sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,ip,ir)*ints(INActive+1:NOccup,INActive+1:NOccup))
!                     
!                     ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                     ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                     
!                  endif
!               enddo
!            endif
!         enddo
!            
!      endif
!         
!      ! T3
!      if(IGem(k)==2.and.IGem(l)==2) then
!         
!         do is=INActive+1,NOccup
!            do ir=INActive+1,NBasis
!               irs = pos(ir,is)
!               if(irs>0) then
!                  do iq=INActive+1,NOccup
!                     do ip=INActive+1,NBasis
!                        ipq = pos(ip,iq)
!                        if(ipq>0) then
!                           
!                           val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
!                                RDM2val(iq,is,k,l)*ints(ip,ir)
!                           
!                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                           ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                           
!                        endif
!                     enddo
!                  enddo
!               endif
!            enddo
!         enddo
!      endif
!
!      ! T2+T5: AuxInd = 2
!      if((k<=NOccup).and.(l>INActive.and.l<=NOccup)) then
!         is = k
!         iu = l
!         
!         ! T2
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=1,NOccup
!                  do ip=INActive+1,NOccup
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,ip,ir,iu)*ints(INActive+1:NOccup,iq))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         ! p->q
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=INActive+1,NOccup
!                  do ip=INActive+1,NBasis
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,iq,ir,iu)*ints(INActive+1:NOccup,ip))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!
!         ! T5
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=INActive+1,NOccup
!                  do ip=INActive+1,NBasis
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = - AuxCoeff(IGem(ip),IGem(is),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,iq,iu,ir)*ints(INActive+1:NOccup,ip))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         ! p->q
!         do ir=INActive+1,NOccup
!            irs = pos(ir,is)
!            if(irs>0) then
!               do iq=1,NOccup
!                  do ip=INActive+1,NOccup
!                     ipq = pos(ip,iq)
!                     if(ipq>0) then
!                        val = - AuxCoeff(IGem(iq),IGem(is),IGem(iu),2)* &
!                             sum(RDM2val(INActive+1:NOccup,ip,iu,ir)*ints(INActive+1:NOccup,iq))
!                        
!                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!                        
!                     endif
!                  enddo
!               enddo
!            endif
!         enddo
!         
!      endif
!      
!   enddo
!enddo
!
!close(iunit2)
!
!do i=1,NBasis
!   C(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
!enddo
!
!do ICol=1,NDimX
!   ir = IndN(1,ICol)
!   is = IndN(2,ICol)
!   irs = IndX(ICol)
!   if(ir/=is) then
!      
!      do IRow=1,NDimX
!         ip = IndN(1,IRow)
!         iq = IndN(2,IRow)
!         ipq = IndX(IRow)
!
!         val = 0
!
!         if(ip==ir) then
!            val = val + (Occ(ip)-Occ(is))*HNO(iq,is) - WMAT(iq,is)
!            select case(AuxInd(IGem(ip),IGem(ir)))
!            case(1)
!               val = val + Occ(ip)*AuxI(iq,is)
!            case(2)
!               val = val + Occ(ip)*AuxIO(iq,is)
!            end select
!         endif
!
!         if(iq==is) then
!            val = val + (Occ(iq)-Occ(ir))*HNO(ip,ir) - WMAT(ip,ir)
!            select case(AuxInd(IGem(iq),IGem(is)))
!            case(1)
!               val = val + Occ(iq)*AuxI(ip,ir)
!            case(2)
!               val = val + Occ(iq)*AuxIO(ip,ir)
!            end select
!         endif
!
!         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
!
!         val = 0
!
!         if(iq==ir) then
!            val = val + (Occ(iq)-Occ(is))*HNO(ip,is) - WMAT(ip,is)
!            select case(AuxInd(IGem(iq),IGem(ir)))
!            case(1)
!               val = val + Occ(iq)*AuxI(ip,is)
!            case(2)
!               val = val + Occ(iq)*AuxIO(ip,is)
!            end select
!         endif
!         
!         if(ip==is) then
!            val = val + (Occ(ip)-Occ(ir))*HNO(iq,ir) - WMAT(iq,ir)
!            select case(AuxInd(IGem(ip),IGem(is)))
!            case(1)
!               val = val + Occ(ip)*AuxI(iq,ir)
!            case(2)
!               val = val + Occ(ip)*AuxIO(iq,ir)
!            end select
!         endif
!
!         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
!         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
!
!         val = (C(ip) + C(iq))*(C(ir) + C(is))
!         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
!         val = (C(ip) - C(iq))*(C(ir) - C(is))
!         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
!
!      enddo
!
!   endif
!enddo
!
!deallocate(C,AuxI,AuxIO,WMAT)
!deallocate(ints,work2,work1)
!
!end subroutine AB_CAS_assemb

subroutine Y01CAS_mithap(Occ,URe,XOne,ABPLUS,ABMIN, &
     EigY,EigY1,Eig,Eig1, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1,IntFileName,IFlag0)
!
!     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
!     IN THE 1ST-ORDER APPROXIMATION
!
!     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
!              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntFileName
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(out) :: EigY(NDimX,NDimX),EigY1(NDimX,NDimX)
double precision,intent(out) :: Eig(NDimX),Eig1(NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis),IFlag0

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
double precision,parameter :: Thresh = 1.D-12

ABPLUS = 0
ABMIN  = 0

! set dimensions
NOccup = NAct + INActive
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
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

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

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

AuxI  = 0
AuxIO = 0
WMAT  = 0

open(newunit=iunit,file=trim(IntFileName),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do ll=1,NBasis
   do kk=1,ll
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
      do klround=1,merge(1,2,kk==ll)
         select case(klround)
         case(1)
            k = kk
            l = ll
         case(2)
            k = ll
            l = kk
         end select

! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
         !! Coulomb
         if(IGem(k)==IGem(l)) then

            val = 0
            do i=1,NBasis
               if(IGem(i)/=IGem(k)) val = val + Occ(i)*ints(i,i)
            enddo
            HNO(k,l) = HNO(k,l) + 2*val

         !! exchange
         else

            val = Occ(k)
            do i=1,NBasis
               if(IGem(i)==IGem(l)) HNO(i,l) = HNO(i,l) - val*ints(i,k)
            enddo

         endif


! AUXILIARY MATRIX AuxI AND AuxIO
         ! Coulomb
         val = 0
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),1,1)
         do i=1,INActive
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxIO(k,l) = AuxIO(k,l) + 2*val
         
         AuxVal = AuxCoeff(IGem(k),IGem(l),2,2)
         do i=INActive+1,NOccup
            val = val + AuxVal*Occ(i)*ints(i,i)
         enddo
         AuxI(k,l) = AuxI(k,l) + 2*val

         ! exchange
         if(k<=INActive) then
            do i=1,NBasis
               AuxIO(i,l) = AuxIO(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
            enddo
         endif
         if(k<=NOccup) then
            do i=1,NBasis
               AuxI(i,l) = AuxI(i,l) - AuxCoeff(IGem(i),IGem(i),IGem(k),IGem(l))*Occ(k)*ints(i,k)
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
                        val = 2*AuxVal*val*ints(iq,ip)

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

                     endif
                  enddo
               enddo

            endif
         endif

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
                           val = - AuxVal*val*ints(iq,ir)

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
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
         if((k>INActive).and.(l>INActive)) then
            ip = k
            ir = l

            do is=INActive+1,NOccup
               irs = pos(ir,is)
               if(irs>0) then
                  do iq=INActive+1,NOccup
                     ipq = pos(ip,iq)
                     if(ipq>0) then
                        val = AuxCoeff(IGem(ip),IGem(ir),2,2)* &
                             sum(RDM2val(INActive+1:NOccup,INActive+1:NOccup,iq,is)*ints(INActive+1:NOccup,INActive+1:NOccup))

                        ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                        ABMIN(ipq,irs) = ABMIN(ipq,irs) + val
                        
                     endif
                  enddo
               endif
            enddo
            
         endif

         if((k<=NOccup).and.(l>INActive)) then
            iq = k
            ir = l

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
                                sum(RDM2val(INActive+1:NOccup,iq,is,iu)*ints(INActive+1:NOccup,ip))

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
                                sum(RDM2val(INActive+1:NOccup,iq,iu,is)*ints(INActive+1:NOccup,ip))

                           ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
                           ABMIN(ipq,irs) = ABMIN(ipq,irs) - val
                           
                        endif
                     enddo
                  enddo
               endif
            enddo

         endif         
         
      enddo
   enddo
enddo

close(iunit)

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

         !if(ipq>=irs) then
         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) + val

         ! ABPLUS(irs,ipq) = ABPLUS(ipq,irs) 
         ! ABMIN(irs,ipq) = ABMIN(ipq,irs) 
         !endif

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

         ! if(ipq>=irs) then
         ABPLUS(ipq,irs) = ABPLUS(ipq,irs) + val
         ABMIN(ipq,irs) = ABMIN(ipq,irs) - val

         ! ABPLUS(irs,ipq) = ABPLUS(ipq,irs) 
         ! ABMIN(irs,ipq) = ABMIN(ipq,irs)
         ! endif

         ! if(ipq>=irs) then
         val = (C(ip) + C(iq))*(C(ir) + C(is))
         if(val/=0d0) ABPLUS(ipq,irs) = ABPLUS(ipq,irs)/val
         val = (C(ip) - C(iq))*(C(ir) - C(is))
         if(val/=0d0) ABMIN(ipq,irs)  = ABMIN(ipq,irs)/val
         
         ! ABPLUS(irs,ipq) = ABPLUS(ipq,irs) 
         !  ABMIN(irs,ipq) = ABMIN(ipq,irs)
         ! endif
         
      enddo

enddo

EigY1 = 0
EigY = 0

i = limAA(1)
j = limAA(2)
!! write(*,*) 'AB0-my',norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!$write(*,*) 'AB+',norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!$write(*,*) 'AB-',norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))

if(nAA>0) then
   allocate(ABP(nAA,nAA),ABM(nAA,nAA),EigYt(nAA,nAA),EigXt(nAA,nAA),Eigt(nAA))

   ABP = ABPLUS(i:j,i:j)
   ABM = ABMIN(i:j,i:j)
   call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAA)
   
   do ii=1,nAA
      ipos = tmpAA(ii)
      EigY(ipos,i:j) = EigYt(ii,1:nAA)
      if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAA)
   enddo
   Eig(i:j) = Eigt(1:nAA)

   deallocate(Eigt,EigXt,EigYt,ABM,ABP)

   do iq=1,INActive
      i = limAI(1,iq)
      j = limAI(2,iq)

!!$   write(*,*) 'ABai-my',iq,norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!$   write(*,*) 'AB+',iq,norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!$   write(*,*) 'AB-',iq,norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))

      allocate(ABP(nAI(iq),nAI(iq)),ABM(nAI(iq),nAI(iq)),&
           EigYt(nAI(iq),nAI(iq)),EigXt(nAI(iq),nAI(iq)),Eigt(nAI(iq)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAI(iq))
   
      do ii=1,nAI(iq)
         ipos = tmpAI(ii,iq)
         EigY(ipos,i:j) = EigYt(ii,1:nAI(iq))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAI(iq))
      enddo
      Eig(i:j) = Eigt(1:nAI(iq))
      
      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
      
   enddo

   do ip=NOccup+1,NBasis
      i = limAV(1,ip)
      j = limAV(2,ip)
!!$   write(*,*) 'ABav-my',ip,norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!$   write(*,*) 'AB+',ip,norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!$   write(*,*) 'AB-',ip,norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))

      allocate(ABP(nAV(ip),nAV(ip)),ABM(nAV(ip),nAV(ip)),&
           EigYt(nAV(ip),nAV(ip)),EigXt(nAV(ip),nAV(ip)),Eigt(nAV(ip)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAV(ip))

      do ii=1,nAV(ip)
         ipos = tmpAV(ii,ip)
         EigY(ipos,i:j) = EigYt(ii,1:nAV(ip))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAV(ip))
      enddo
      Eig(i:j) = Eigt(1:nAV(ip))
   
      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
      
   enddo
endif

do i=limIV(1),limIV(2)
   ii = tmpIV(i-limIV(1)+1)
   ip = IndN(1,ii)
   iq = IndN(2,ii)
   ! write(*,*) 'ABiv-my',ip,iq,ABPLUS(i,i)
   Eig(i) = ABPLUS(i,i)
enddo
do ii=1,nIV
   i = tmpiV(ii)
   EigY(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
   if(IFlag0==0) EigY1(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
enddo

deallocate(work1)

if(IFlag0==1) return


! AB(1) PART
call AB_CAS_mithap(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntFileName,1d0,.true.)

allocate(work1(NDimX**2))
! work1=ABPLUS.EigX
! ABPLUS=work1.EigX
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,EigY1,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)

! work1=ABMIN.EigY
! ABMIN=work1.EigY
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,EigY,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY,NDimX,0d0,ABMIN,NDimX)

do i=1,NDimX
   Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
enddo

EigY1 = 0
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) then
            val = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            if(Abs(Eig(i)-Eig(j))>Thresh) then
               val = val + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
            endif
            do ii=1,NDimX
               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
            enddo
         endif
      enddo
   endif
enddo

! TESTS
!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-m1', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$HNO=transpose(HNO)
!!$call sq_to_triang2(HNO,work1,NBasis)
!!$write(LOUT,*) 'HNO-m2', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$!do i=1,NBasis*(NBasis+1)/2
!!$!   print*, i, work1(i)
!!$!enddo
!!$
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxI=transpose(AuxI)
!!$call sq_to_triang2(AuxI,work1,NBasis)
!!$write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$AuxIO=transpose(AuxIO)
!!$call sq_to_triang2(AuxIO,work1,NBasis)
!!$write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!!$
!!$write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)
!!$write(LOUT,*) 'WMAT-norm',norm2(WMAT(1:NOccup,1:NOccup)-transpose(WMAT(1:NOccup,1:NOccup)))


deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine Y01CAS_mithap

subroutine Y01CAS_FOFO(Occ,URe,XOne,ABPLUS,ABMIN, &
     EigY,EigY1,Eig,Eig1, &
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     IntFileName,IntJFile,IntKFile,IFlag0)
!
!     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
!     IN THE 1ST-ORDER APPROXIMATION
!
!     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
!              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1
character(*) :: IntFileName,IntJFile,IntKFile
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
double precision,intent(out) :: EigY(NDimX,NDimX),EigY1(NDimX,NDimX)
double precision,intent(out) :: Eig(NDimX),Eig1(NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1)
integer,intent(in) :: IndN(2,NDim),IndX(NDim),IGemIN(NBasis),IFlag0

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: iunit1,iunit2
integer :: NOccup,NRDM2Act
integer :: IGem(NBasis),Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: C(NBasis)
double precision :: HNO(NBasis,NBasis),AuxI(NBasis,NBasis),AuxIO(NBasis,NBasis),WMAT(NBasis,NBasis)
double precision :: AuxCoeff(3,3,3,3),AuxVal,val
double precision :: EnDummy
double precision,allocatable :: RDM2val(:,:,:,:),RDM2Act(:)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
double precision,parameter :: Thresh = 1.D-12

ABPLUS = 0
ABMIN  = 0

! set dimensions
NOccup = NAct + INActive
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
allocate(RDM2val(NOccup,NOccup,NOccup,NOccup))

NRDM2Act = NAct**2*(NAct**2+1)/2
allocate(RDM2Act(NRDM2Act))

RDM2Act = 0
open(newunit=iunit,file='rdm2.dat',status='old')
write(LOUT,'(/,1x,''Active block of 2-RDM read from rdm2.dat'')')
do
   read(iunit,*,iostat=ios) i,j,k,l,val
   if(ios/=0) exit
   RDM2Act(NAddrRDM(i,k,j,l,NAct)) = 0.5d0*val
enddo
close(iunit)

do l=1,NOccup
   do k=1,NOccup
      do j=1,NOccup
         do i=1,NOccup
            RDM2val(i,j,k,l) = FRDM2(i,k,j,l,RDM2Act,Occ,Ind,NAct,NBasis)
         enddo
      enddo
   enddo
enddo

deallocate(RDM2Act)

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

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

nAA = 0
nAI = 0
nAV = 0
nIV = 0
tmpAA = 0
tmpAI = 0
tmpAV = 0
tmpIV = 0
limAA = 0
limAI = 0
limAV = 0
limIV = 0

do i=1,NDimX
   ip = IndN(1,i)
   iq = IndN(2,i)
   if(IGem(ip)==2.and.IGem(iq)==2) then
      nAA = nAA + 1
      tmpAA(nAA) = i
    elseif(IGem(ip)==2.and.IGem(iq)==1) then
       nAI(iq) = nAI(iq) + 1
       tmpAI(nAI(iq),iq) = i
    elseif(IGem(ip)==3.and.IGem(iq)==2) then
       nAV(ip) = nAV(ip) + 1
       tmpAV(nAV(ip),ip) = i
    elseif(IGem(ip)==3.and.IGem(iq)==1) then
       nIV = nIV + 1
       tmpIV(nIV) = i
   endif
enddo

pos = 0
ipos = 0
limAA(1) = ipos + 1
do ii=1,nAA
   i = tmpAA(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limAA(2) = ipos
do is=1,INActive
   limAI(1,is) = ipos + 1
   do ii=1,nAI(is)
      i = tmpAI(ii,is)
      ! print*, ii,is,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAI(2,is) = ipos
enddo
!print*, ' '
do ir=NOccup+1,NBasis
   limAV(1,ir) = ipos + 1
   do ii=1,nAV(ir)
      i = tmpAV(ii,ir)
      ! print*, ii,ir,i
      ip = IndN(1,i)
      iq = IndN(2,i)
      ipos = ipos + 1
      pos(ip,iq) = ipos      
   enddo
   limAV(2,ir) = ipos
enddo
limIV(1) = ipos + 1
do ii=1,nIV
   i = tmpIV(ii)
   ! print*, ii,i
   ip = IndN(1,i)
   iq = IndN(2,i)
   ipos = ipos + 1
   pos(ip,iq) = ipos
enddo
limIV(2) = ipos

!not working... ???
!call AB_CAS_assemb(ABPLUS,ABMIN,Occ,RDM2val,HNO,pos, &
!           1d0,AuxInd,AuxCoeff, &
!           IGem,IndN,IndX,NDimX,NDim,INActive,NOccup,NBasis, &
!           IntJFile,IntKFile)

AuxI  = 0
AuxIO = 0
WMAT  = 0


open(newunit=iunit1,file=trim(IntKFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*NOccup)

! exchange loop (FO|FO)
kl = 0
do ll=1,NOccup
   do kk=1,NBasis
      kl = kl + 1
      read(iunit1,rec=kl) work1(1:NBasis*NOccup)
      !call triang_to_sq2(work1,ints,NBasis)
      do j=1,NOccup
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      k = kk
      l = ll

      if(l>NOccup) cycle
      ints(:,NOccup+1:NBasis) = 0

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! exchange
      if(l<=NOccup) then
         val = Occ(l)
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
enddo

close(iunit1)

open(newunit=iunit2,file=trim(IntJFile),status='OLD', &
     access='DIRECT',recl=8*NBasis**2)

! Coulomb loop (FF|OO)
kl = 0
do ll=1,NOccup
   do kk=1,NOccup
      kl = kl + 1
      read(iunit2,rec=kl) work1(1:NBasis**2)
      do j=1,NBasis
         do i=1,NBasis
            ints(i,j) = work1((j-1)*NBasis+i)
         enddo
      enddo
      
      k = kk
      l = ll

      if(k>NOccup.or.l>NOccup) cycle

      ! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN

      ! Coulomb
      if(k==l.and.k<=NOccup) then
         val = 2*Occ(k)
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
enddo

close(iunit2)

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


EigY1 = 0
EigY = 0

i = limAA(1)
j = limAA(2)

if(nAA>0) then
   allocate(ABP(nAA,nAA),ABM(nAA,nAA),EigYt(nAA,nAA),EigXt(nAA,nAA),Eigt(nAA))

   ABP = ABPLUS(i:j,i:j)
   ABM = ABMIN(i:j,i:j)
   call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAA)
   
   do ii=1,nAA
      ipos = tmpAA(ii)
      EigY(ipos,i:j) = EigYt(ii,1:nAA)
      if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAA)
   enddo
   Eig(i:j) = Eigt(1:nAA)

   deallocate(Eigt,EigXt,EigYt,ABM,ABP)

   do iq=1,INActive
      i = limAI(1,iq)
      j = limAI(2,iq)

      allocate(ABP(nAI(iq),nAI(iq)),ABM(nAI(iq),nAI(iq)),&
           EigYt(nAI(iq),nAI(iq)),EigXt(nAI(iq),nAI(iq)),Eigt(nAI(iq)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAI(iq))
   
      do ii=1,nAI(iq)
         ipos = tmpAI(ii,iq)
         EigY(ipos,i:j) = EigYt(ii,1:nAI(iq))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAI(iq))
      enddo
      Eig(i:j) = Eigt(1:nAI(iq))
      
      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
      
   enddo

   do ip=NOccup+1,NBasis
      i = limAV(1,ip)
      j = limAV(2,ip)

      allocate(ABP(nAV(ip),nAV(ip)),ABM(nAV(ip),nAV(ip)),&
           EigYt(nAV(ip),nAV(ip)),EigXt(nAV(ip),nAV(ip)),Eigt(nAV(ip)))

      ABP = ABPLUS(i:j,i:j)
      ABM = ABMIN(i:j,i:j)
      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAV(ip))

      do ii=1,nAV(ip)
         ipos = tmpAV(ii,ip)
         EigY(ipos,i:j) = EigYt(ii,1:nAV(ip))
         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAV(ip))
      enddo
      Eig(i:j) = Eigt(1:nAV(ip))
   
      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
      
   enddo
endif

do i=limIV(1),limIV(2)
   ii = tmpIV(i-limIV(1)+1)
   ip = IndN(1,ii)
   iq = IndN(2,ii)
   ! write(*,*) 'ABiv-my',ip,iq,ABPLUS(i,i)
   Eig(i) = ABPLUS(i,i)
enddo
do ii=1,nIV
   i = tmpiV(ii)
   EigY(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
   if(IFlag0==0) EigY1(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
enddo

deallocate(work1)


if(IFlag0==1) return
!return

! AB(1) PART
call AB_CAS_FOFO(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntJFile,IntKFile,1d0,.true.)

! here!!!! can this be made cheaper?
allocate(work1(NDimX**2))
! work1=ABPLUS.EigX
! ABPLUS=work1.EigX
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,EigY1,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)

! work1=ABMIN.EigY
! ABMIN=work1.EigY
call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,EigY,NDimX,0d0,work1,NDimX) 
call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY,NDimX,0d0,ABMIN,NDimX)

do i=1,NDimX
   Eig1(i)=ABPLUS(i,i)+ABMIN(i,i)
enddo

EigY1 = 0
do j=1,NDimX
   if(Eig(j)/=0d0) then
      do i=1,NDimX
         if(Eig(i)/=0d0) then
            val = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
            if(Abs(Eig(i)-Eig(j))>Thresh) then
               val = val + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
            endif
            do ii=1,NDimX
               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
            enddo
         endif
      enddo
   endif
enddo

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine Y01CAS_FOFO

subroutine ACABMAT0_mithap(AMAT,BMAT,URe,Occ,XOne,&
               IndN,IndX,IGem,C, &
               NAct,INActive,NBasis,NDim,NDimX,NInte1,NGem,IntFileName,ISAPT,ACAlpha,IFlag)
!     IFlag = 1 - AMAT AND BMAT WILL CONTAIN (A+B)/C+/C+ AND (A-B)/C-/C-, RESPECTIVELY
!             0 - AMAT AND BMAT WILL CONTAIN A ANB B MATRICES, RESPECTIVELY
!
!     ACAlpha - Alpha-connection parameter in AC
!     HNO AND TwoMO are modified to correspond to an alpha-Hamiltonian
!
!     STRAIGHTFORWARD IMPLEMENTATION OF THE AMAT AND BMAT DEFINITIONS (NO SPECIAL CASES CONSIDERED) 
!     
!     COMPUTE THE A+B AND A-B MATRICES IN ERPA WITH APSG APPROXIMATION
!
!     NESTED COMMUTATOR 
!
implicit none

integer,intent(in) :: NAct,INActive,NBasis,NDim,NDimX,NInte1,ISAPT,NGem
character(*) :: IntFileName
double precision,intent(out) :: AMAT(NDimX,NDimX),BMAT(NDimX,NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1),C(NBasis)
double precision,intent(in)  :: ACAlpha
integer,intent(in) :: IndN(2,NDim),IndX(NDim)
integer,intent(in) :: IGem(NBasis),IFlag

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: NOccup
integer :: IGemType
integer :: Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: HNO(NBasis,NBasis),HNOCoef
double precision :: AuxCoeff(NGem,NGem,NGem,NGem),AuxVal,val
double precision :: OccProd(NBasis,NBasis),CProd(NBasis,NBasis)
double precision :: AuxH(NBasis,NBasis,NGem),AuxXC(NBasis,NBasis,NGem)
double precision :: SaveA,SaveB
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,parameter :: Delta = 1.d-6

print*, 'Wszystko od nowa!'
if(ISAPT==1) then
   write(6,'(1x,a)') 'Computing response-my'
else
   write(6,'(/,X,"***** COMPUTING AMAT, BMAT IN ACABMAT0 *****",/)')
endif

AMAT = 0
BMAT = 0

! set dimensions
NOccup = NAct + INActive

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)/=IGem(j)) HNO(i,j) = ACAlpha*HNO(i,j)
   enddo
enddo

do l=1,NGem
   do k=1,NGem
      do j=1,NGem
         do i=1,NGem
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = 1
            else
               AuxCoeff(i,j,k,l) = ACAlpha
            endif
         enddo
      enddo
   enddo
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

CProd = 0
OccProd = 0
do j=1,NBasis
   do i=1,NBasis
      if(IGem(i)==IGem(j)) then
         CProd(i,j) = C(i)*C(j)         
      else
         OccProd(i,j) = Occ(i)*Occ(j)
      endif
   enddo
enddo

AuxH = 0
AuxXC = 0

HNOCoef = 1 - ACAlpha

open(newunit=iunit,file=trim(IntFileName),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do ll=1,NBasis
   do kk=1,ll
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
      do klround=1,merge(1,2,kk==ll)
         select case(klround)
         case(1)
            k = kk
            l = ll
         case(2)
            k = ll
            l = kk
         end select

! CONSTRUCT ONE-ELECTRON PART OF THE AC ALPHA-HAMILTONIAN
         ! Coulomb
         if(IGem(k)==IGem(l)) then

            val = 0
            do i=1,NBasis
               if(IGem(i)/=IGem(k)) val = val + Occ(i)*ints(i,i)
            enddo
            val = 2*HNOCoef*val
            HNO(k,l) = HNO(k,l) + val
            
            ! exchange
         else

            val = HNOCoef*Occ(k)
            do i=1,NBasis
               if(IGem(i)==IGem(l)) HNO(i,l) = HNO(i,l) - val*ints(i,k)
            enddo

         endif

         do IGemType=1,NGem
            
            ! Coulomb
            val = 0
            do i=1,NBasis
               if(IGem(i)/=IGemType) val = val + &
                    AuxCoeff(IGem(k),IGem(l),IGem(i),IGem(i))*Occ(i)*ints(i,i)
            enddo
            AuxH(k,l,IGemType) = AuxH(k,l,IGemType) + 2*val

            ! exchange
            if(IGem(k)==IGemType) then
               do i=1,NBasis
                  AuxXC(i,l,IGemType) = AuxXC(i,l,IGemType) + AuxCoeff(IGem(i),IGem(k),IGem(k),IGem(l))*C(k)*ints(i,k)
               enddo
            else
               do i=1,NBasis
                  AuxH(i,l,IGemType) = AuxH(i,l,IGemType) - AuxCoeff(IGem(i),IGem(k),IGem(k),IGem(l))*Occ(k)*ints(i,k)
               enddo
            endif
            
         enddo

         ! INTERGEMINAL PART
         ir = k
         is = l
         irs = pos(ir,is)
         if(irs>0) then
            do ipq=1,NDimX
               ip = IndN(1,ipq)
               iq = IndN(2,ipq)

               val = 2*AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                    *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ip,iq)   
               BMAT(ipq,irs) = BMAT(ipq,irs) + val
               AMAT(ipq,irs) = AMAT(ipq,irs) - val
                   
            enddo
         endif

         ip = k
         is = l
         do iq=1,ip-1
            ipq = pos(ip,iq)
            if(ipq>0) then
               do ir=is+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then
                     
                     val = -AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ir,iq)
                     BMAT(ipq,irs) = BMAT(ipq,irs) + val
                          
                  endif
               enddo
            endif
         enddo

         ! p->q
         iq = k
         is = l
         do ip=iq+1,NBasis
            ipq = pos(ip,iq)
            if(ipq>0) then
               do ir=is+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then
                     
                     val = -AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(-OccProd(ip,ir)+OccProd(ip,is)+OccProd(iq,ir)-OccProd(iq,is))*ints(ir,ip)
                     AMAT(ipq,irs) = AMAT(ipq,irs) - val
                          
                  endif
               enddo
            endif
         enddo

         ! INTRAGEMINAL PART
         iq = k
         ir = l
         do ip=iq+1,NBasis
            ipq = pos(ip,iq)
            if(ipq>0) then
               do is=1,ir-1 !NBasis
                  irs = pos(ir,is)
                  if(irs>0) then

                     val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(CProd(iq,ir)+CProd(ip,is))*ints(ip,is)
                     BMAT(ipq,irs) = BMAT(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
         
         ! p->q
         ip = k
         ir = l
         do iq=1,ip-1
            ipq = pos(ip,iq)
            if(ipq>0) then
               do is=1,ir-1
                  irs = pos(ir,is)
                  if(irs>0) then

                     val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(CProd(ip,ir)+CProd(iq,is))*ints(iq,is)
                     AMAT(ipq,irs) = AMAT(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo

         iq = k
         is = l
         do ip=iq+1,NBasis
            ipq = pos(ip,iq)
            if(ipq>0) then
               do ir=is+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then

                     val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(CProd(iq,ir)+CProd(ip,is))*ints(ip,ir)
                     BMAT(ipq,irs) = BMAT(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
         
         ! p->q       
         ip = k
         is = l
         do iq=1,ip-1
            ipq = pos(ip,iq)
            if(ipq>0) then
               do ir=is+1,NBasis
                  irs = pos(ir,is)
                  if(irs>0) then

                     val = AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is)) &
                          *(CProd(ip,ir)+CProd(iq,is))*ints(iq,ir)
                     AMAT(ipq,irs) = AMAT(ipq,irs) + val
                     
                  endif
               enddo
            endif
         enddo
         
      enddo
   enddo
enddo

close(iunit)

do ICol=1,NDimX
   ip = IndN(1,ICol)
   iq = IndN(2,ICol)
   ipq = ICol
 
   if(ip/=iq) then
      do IRow=1,NDimX
         ir = IndN(1,IRow)
         is = IndN(2,IRow)
         irs = IRow

            val = 0
            if(iq==ir) val = val &
                 + (Occ(ir)-Occ(is))*HNO(ip,is) &
                 + Occ(iq)*AuxH(ip,is,IGem(iq)) &
                 - Occ(is)*AuxH(ip,is,IGem(is)) &
                 - C(is)*AuxXC(ip,is,IGem(is))
            if(ip==is) val = val &
                 - (Occ(ir)-Occ(is))*HNO(iq,ir) &
                 + Occ(ip)*AuxH(iq,ir,IGem(ip)) &
                 - Occ(ir)*AuxH(iq,ir,IGem(ir)) &
                 - C(ir)*AuxXC(iq,ir,IGem(ir))
            if(ip==ir) val = val &
                 - C(ip)*AuxXC(iq,is,IGem(ip))
            if(iq==is) val = val &
                 - C(iq)*AuxXC(ip,ir,IGem(iq))
            BMAT(irs,ipq) = BMAT(irs,ipq) + val
         
            val = 0
            if(ip==ir) val = val &
                 + (Occ(ir)-Occ(is))*HNO(iq,is) &
                 + Occ(ip)*AuxH(iq,is,IGem(ip)) &
                 - Occ(is)*AuxH(iq,is,IGem(is)) &
                 - C(is)*AuxXC(iq,is,IGem(is))
            if(iq==is) val = val &
                 - (Occ(ir)-Occ(is))*HNO(ip,ir) &
                 + Occ(iq)*AuxH(ip,ir,IGem(iq)) &
                 - Occ(ir)*AuxH(ip,ir,IGem(ir)) & 
                 - C(ir)*AuxXC(ip,ir,IGem(ir))
            if(iq==ir) val = val &
                 - C(iq)*AuxXC(ip,is,IGem(iq))
            if(ip==is) val = val &
                 - C(ip)*AuxXC(iq,ir,IGem(ip))
            AMAT(irs,ipq) = AMAT(irs,ipq) + val

            if(IFlag/=0) then
!!$               if(abs(abs(C(ip))-abs(C(iq)))<=Delta*abs(C(iq)).or.& 
!!$                    abs(abs(C(ir))-abs(C(is)))<=Delta*abs(C(ir))) then
!!$               
!!$                  AMAT(irs,ipq) = 0
!!$                  BMAT(irs,ipq) = 0
!!$               else
                  
                  SaveA = AMAT(irs,ipq)
                  SaveB = BMAT(irs,ipq)
                  
                  val = (C(ip) + C(iq))*(C(ir) + C(is))
                  if(val==0) then
                     AMAT(irs,ipq) = 0
                  else
                     AMAT(irs,ipq) = (SaveA+SaveB)/val
                  endif
                  val = (C(ip) - C(iq))*(C(ir) - C(is))
                  if(val==0) then
                     BMAT(irs,ipq) = 0
                  else
                     BMAT(irs,ipq) = (SaveA-SaveB)/val
                  endif
               endif
!!$            endif
                        
         enddo
      endif
   enddo

print*, "AB-my",norm2(AMAT),norm2(BMAT)

call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
HNO=transpose(HNO)
call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

!!$do IGemType=1,NGem
!!$   call sq_to_triang2(AuxH(:,:,IGemType),work1,NBasis)
!!$   write(LOUT,*) 'AuxH-my', IGemType,norm2(work1(1:NBasis*(NBasis+1)/2))
!!$   AuxH(:,:,IGemType)=transpose(AuxH(:,:,IGemType))
!!$   call sq_to_triang2(AuxH(:,:,IGemType),work1,NBasis)
!!$   write(LOUT,*) 'AuxH-tr', IGemType,norm2(work1(1:NBasis*(NBasis+1)/2))
!!$enddo
!!$
!!$print*, ''
!!$do IGemType=1,NGem
!!$   call sq_to_triang2(AuxXC(:,:,IGemType),work1,NBasis)
!!$   write(LOUT,*) 'AuxXC-my', IGemType,norm2(work1(1:NBasis*(NBasis+1)/2))
!!$   AuxXC(:,:,IGemType)=transpose(AuxXC(:,:,IGemType))
!!$   call sq_to_triang2(AuxXC(:,:,IGemType),work1,NBasis)
!!$   write(LOUT,*) 'AuxXC-tr', IGemType,norm2(work1(1:NBasis*(NBasis+1)/2))
!!$enddo

deallocate(ints,work2,work1)

end subroutine ACABMAT0_mithap

subroutine FockGen_mithap(Fock,OneRdm,XOne,NInte1,NBasis,IntFileName)
!
!     GENERALIZED FOCK MATRIX
!
implicit none

integer,intent(in) :: NInte1,NBasis
character(*) :: IntFileName
double precision,intent(in) :: OneRdm(NInte1),XOne(NInte1)
double precision,intent(out) :: Fock(Ninte1)
integer :: iunit,kk,ll,kl,klround,k,l
double precision,allocatable :: OneRdmSq(:,:),FockSq(:,:),ints(:,:),work1(:)


allocate(OneRdmSq(NBasis,NBasis),FockSq(NBasis,NBasis),ints(NBasis,NBasis),work1(NInte1))

call triang_to_sq2(OneRdm,OneRdmSq,NBasis)
call triang_to_sq2(XOne,FockSq,NBasis)

open(newunit=iunit,file=trim(IntFileName),status='OLD',&
     access='DIRECT',form='UNFORMATTED',recl=8*NInte1)

kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
       
      if(k==l) then
        call daxpy(NBasis**2,2.d0*OneRdmSq(k,l),ints,1,FockSq,1)
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,k),1,1.d0,FockSq(:,l),1) 
      else
        call daxpy(NBasis**2,2.d0*(OneRdmSq(k,l)+OneRdmSq(l,k)),ints,1,FockSq,1)
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,k),1,1.d0,FockSq(:,l),1) 
        call dgemv('N',NBasis,NBasis,-1.d0,ints,NBasis,OneRdmSq(:,l),1,1.d0,FockSq(:,k),1) 
      endif

   enddo
enddo

close(iunit)

call sq_to_triang2(FockSq,Fock,NBasis) 

deallocate(work1,ints,FockSq,OneRdmSq)

end subroutine FockGen_mithap

subroutine PotHSR_mithap(VHSR,Occ,MO,NBasis)
!
!     RETURNS SR HARTREE (COULOMB) POTENTIAL IN AN MO MATRIX REPRESENTATION 
!
implicit none

integer,intent(in) :: NBasis
double precision,intent(in) :: Occ(NBasis),MO(NBasis,NBasis)
double precision,intent(out) :: VHSR(NBasis*(NBasis+1)/2)
integer :: i
double precision,allocatable :: J(:,:),Jerf(:,:),work(:,:)

allocate(Jerf(NBasis,NBasis),J(NBasis,NBasis),work(NBasis,NBasis))

work = 0
do i=1,NBasis
   call dger(NBasis,NBasis,2d0*Occ(i),MO(i,:),1,MO(i,:),1,work,NBasis)
enddo

call make_J1(NBasis,work,J,'AOTWOSORT')
call make_J1(NBasis,work,Jerf,'AOERFSORT')
work = J - Jerf
call sq_to_triang2(work,VHSR,NBasis) 
call tran_matTr(VHSR,MO,MO,NBasis,.false.)

deallocate(work,J,Jerf)

end subroutine PotHSR_mithap

subroutine PotCoul_mithap(VHSR,OneRdm,doRSH,NBasis)
!
!     RETURNS SR HARTREE (COULOMB) POTENTIAL IN AN MO MATRIX REPRESENTATION 
!
implicit none

integer,intent(in) :: NBasis
logical,intent(in) :: doRSH
double precision,intent(in) :: OneRdm(NBasis*(NBasis+1)/2)
double precision :: VHSR(NBasis*(NBasis+1)/2)
double precision,allocatable :: J(:,:),Jerf(:,:),work(:,:)

allocate(Jerf(NBasis,NBasis),J(NBasis,NBasis),work(NBasis,NBasis))

work = 0
call triang_to_sq2(OneRdm,work,NBasis)

if(doRSH) then
  call make_J1(NBasis,work,J,'AOTWOSORT')
  call make_J1(NBasis,work,Jerf,'AOERFSORT')
  work = 2.0d0*J - 2.0d0*Jerf
else
  call make_J1(NBasis,work,J,'AOTWOSORT')
  work = 2.0d0*J
endif
call sq_to_triang2(work,VHSR,NBasis) 

deallocate(work,J,Jerf)

end subroutine PotCoul_mithap

end module abmat
