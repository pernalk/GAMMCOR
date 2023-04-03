module abmat

use print_units
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

call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
HNO=transpose(HNO)
call sq_to_triang2(HNO,work1,NBasis)
write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

write(LOUT,'(/,1X,''CASSCF Energy (w/o ENuc)'',5X,F15.8)') ETot

call sq_to_triang2(AuxI,work1,NBasis)
write(LOUT,*) 'AuxI-my', norm2(work1(1:NBasis*(NBasis+1)/2))
AuxI=transpose(AuxI)
call sq_to_triang2(AuxI,work1,NBasis)
write(LOUT,*) 'AuxI-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

call sq_to_triang2(AuxIO,work1,NBasis)
write(LOUT,*) 'AuxIO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
AuxIO=transpose(AuxIO)
call sq_to_triang2(AuxIO,work1,NBasis)
write(LOUT,*) 'AuxIO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

write(LOUT,*) 'WMAT-my', 2d0*norm2(WMAT)

deallocate(RDM2val)
deallocate(ints,work2,work1)

end subroutine AB_CAS_mithap

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
     propfile0,propfile1, & 
     y01file,&
     IndN,IndX,IGemIN,NAct,INActive,NDimX,NBasis,NDim,NInte1, &
     NoSt,IntFileName,IFlag0)
!
!     A ROUTINE FOR COMPUTING Y VECTORS AND EIGENVALUES OF ERPA 
!     IN THE 1ST-ORDER APPROXIMATION
!
!     IFlag0 = 1 - compute only 0th-order Y [EigY] and 0th-order omega [Eig] 
!              0 - compute both 0th-order and 1st-order Y [EigY1] and omega [Eig1]
implicit none

integer,intent(in) :: NAct,INActive,NDimX,NBasis,NDim,NInte1,NoSt
character(*) :: IntFileName
character(*) :: propfile0,propfile1,y01file
double precision,intent(out) :: ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX)
!double precision,intent(out) :: EigY(NDimX,NDimX),EigY1(NDimX,NDimX)
!double precision,intent(out) :: Eig(NDimX),Eig1(NDimX)
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
double precision,allocatable :: work1(:),work2(:),workA(:,:)
double precision,allocatable :: ints(:,:)
double precision,allocatable :: EigY(:,:),EigY1(:,:)
double precision,allocatable :: Eig(:),Eig1(:)
double precision,allocatable :: EigYt(:,:),EigXt(:,:),Eigt(:)
double precision,allocatable :: ABP(:,:),ABM(:,:)
integer,external :: NAddrRDM
double precision,external :: FRDM2
integer :: ii,jj,ipos,jpos,iblk,jblk,nblk
integer :: nAA,tmpAA(NAct*(NAct-1)/2),limAA(2)
integer :: nAI(INActive),tmpAI(NAct,1:INActive),limAI(2,1:INActive)
integer :: nAV(INActive+NAct+1:NBasis),tmpAV(NAct,INActive+NAct+1:NBasis),limAV(2,INActive+NAct+1:NBasis)
integer :: nIV,tmpIV(INActive*(NBasis-NAct-INActive)),limIV(2)
!
type(EblockData) :: Eblock(1+NBasis-NAct)
type(EblockData) :: EblockIV
!
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

!! MH: 18.12.2019 : old code
!EigY1 = 0
!EigY = 0
!
!i = limAA(1)
!j = limAA(2)
!!! write(*,*) 'AB0-my',norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!!$write(*,*) 'AB+',norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!!$write(*,*) 'AB-',norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))
!
!if(nAA>0) then
!   allocate(ABP(nAA,nAA),ABM(nAA,nAA),EigYt(nAA,nAA),EigXt(nAA,nAA),Eigt(nAA))
!
!   ABP = ABPLUS(i:j,i:j)
!   ABM = ABMIN(i:j,i:j)
!   call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAA)
!   
!   do ii=1,nAA
!      ipos = tmpAA(ii)
!      EigY(ipos,i:j) = EigYt(ii,1:nAA)
!      if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAA)
!   enddo
!   Eig(i:j) = Eigt(1:nAA)
!
!   deallocate(Eigt,EigXt,EigYt,ABM,ABP)
!endif
!   
!   do iq=1,INActive
!      i = limAI(1,iq)
!      j = limAI(2,iq)
!      if(nAI(iq)>0) then
!!!$   write(*,*) 'ABai-my',iq,norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!!$   write(*,*) 'AB+',iq,norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!!$   write(*,*) 'AB-',iq,norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))
!
!      allocate(ABP(nAI(iq),nAI(iq)),ABM(nAI(iq),nAI(iq)),&
!           EigYt(nAI(iq),nAI(iq)),EigXt(nAI(iq),nAI(iq)),Eigt(nAI(iq)))
!
!      ABP = ABPLUS(i:j,i:j)
!      ABM = ABMIN(i:j,i:j)
!      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAI(iq))
!   
!      do ii=1,nAI(iq)
!         ipos = tmpAI(ii,iq)
!         EigY(ipos,i:j) = EigYt(ii,1:nAI(iq))
!         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAI(iq))
!      enddo
!      Eig(i:j) = Eigt(1:nAI(iq))
!      
!      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
!      endif 
!   enddo
!
!   do ip=NOccup+1,NBasis
!      i = limAV(1,ip)
!      j = limAV(2,ip)
!!!$   write(*,*) 'ABav-my',ip,norm2(ABPLUS(i:j,i:j)),norm2(ABMIN(i:j,i:j))
!!!$   write(*,*) 'AB+',ip,norm2(ABPLUS(i:j,i:j)-transpose(ABPLUS(i:j,i:j)))
!!!$   write(*,*) 'AB-',ip,norm2(ABMIN(i:j,i:j)-transpose(ABMIN(i:j,i:j)))
!      if(nAV(ip)>0) then
!      allocate(ABP(nAV(ip),nAV(ip)),ABM(nAV(ip),nAV(ip)),&
!           EigYt(nAV(ip),nAV(ip)),EigXt(nAV(ip),nAV(ip)),Eigt(nAV(ip)))
!
!      ABP = ABPLUS(i:j,i:j)
!      ABM = ABMIN(i:j,i:j)
!      call ERPASYMM0(EigYt,EigXt,Eigt,ABP,ABM,nAV(ip))
!
!      do ii=1,nAV(ip)
!         ipos = tmpAV(ii,ip)
!         EigY(ipos,i:j) = EigYt(ii,1:nAV(ip))
!         if(IFlag0==0) EigY1(ipos,i:j) = EigXt(ii,1:nAV(ip))
!      enddo
!      Eig(i:j) = Eigt(1:nAV(ip))
!   
!      deallocate(Eigt,EigXt,EigYt,ABM,ABP)
!      endif
!   enddo
!!endif
!
!do i=limIV(1),limIV(2)
!   ii = tmpIV(i-limIV(1)+1)
!   ip = IndN(1,ii)
!   iq = IndN(2,ii)
!   ! write(*,*) 'ABiv-my',ip,iq,ABPLUS(i,i)
!   Eig(i) = ABPLUS(i,i)
!enddo
!do ii=1,nIV
!   i = tmpiV(ii)
!   EigY(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
!   if(IFlag0==0) EigY1(i,limIV(1)+ii-1) = 1d0/sqrt(2d0)
!enddo
!
!deallocate(work1)
!
!if(IFlag0==1) return
!
!! AB(1) PART
!call AB_CAS_mithap(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
!              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
!              NInte1,IntFileName,1d0,.true.)
!
!allocate(work1(NDimX**2))
!! work1=ABPLUS.EigX
!! ABPLUS=work1.EigX
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABPLUS,NDimX,EigY1,NDimX,0d0,work1,NDimX) 
!call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)
!
!! work1=ABMIN.EigY
!! ABMIN=work1.EigY
!call dgemm('N','N',NDimX,NDimX,NDimX,1d0,ABMIN,NDimX,EigY,NDimX,0d0,work1,NDimX) 
!call dgemm('T','N',NDimX,NDimX,NDimX,1d0,work1,NDimX,EigY,NDimX,0d0,ABMIN,NDimX)

! MH: 18.12.2019, new code
nblk = 0

! put pack into separate subroutine
!pack AA
if(nAA>0) then

   nblk = nblk + 1
   associate(B => Eblock(nblk))
     B%n  = nAA
     B%l1 = limAA(1)
     B%l2 = limAA(2)
     allocate(B%pos(B%n))
     B%pos(1:B%n) = tmpAA(1:B%n)

     allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
     allocate(ABP(B%n,B%n),ABM(B%n,B%n))

     ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
     ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
     if(NoSt==1) then
        call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     elseif(NoSt>1) then
        call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
     endif
   
     deallocate(ABM,ABP)
   
   end associate

endif

!pack AI
do iq=1,INActive
   if(nAI(iq)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAI(iq)
        B%l1 = limAI(1,iq)
        B%l2 = limAI(2,iq)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAI(1:B%n,iq)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif
      
        deallocate(ABM,ABP)

      end associate

   endif
enddo

!pack AV
do ip=NOccup+1,NBasis
   if(nAV(ip)>0) then

      nblk = nblk + 1
      associate(B => Eblock(nblk))
        B%n  = nAV(ip)
        B%l1 = limAV(1,ip)
        B%l2 = limAV(2,ip)
        allocate(B%pos(B%n))
        B%pos(1:B%n) = tmpAV(1:B%n,ip)

        allocate(B%vec(B%n),B%matX(B%n,B%n),B%matY(B%n,B%n))
        allocate(ABP(B%n,B%n),ABM(B%n,B%n))
  
        ABP = ABPLUS(B%l1:B%l2,B%l1:B%l2)
        ABM = ABMIN(B%l1:B%l2,B%l1:B%l2)
        if(NoSt==1) then
           call ERPASYMM0(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        elseif(NoSt>1) then
           call ERPAVECYX(B%matY,B%matX,B%vec,ABP,ABM,B%n)
        endif
      
        deallocate(ABM,ABP)

      end associate

    endif
enddo

!pack IV
associate(B => EblockIV)

  B%l1 = limIV(1)
  B%l2 = limIV(2)
  B%n  = B%l2-B%l1+1
  allocate(B%pos(B%n))
  B%pos(1:B%n) = tmpIV(1:B%n)

  allocate(B%vec(B%n))
  do i=1,B%n
     ii = B%l1+i-1
     B%vec(i) = ABPLUS(ii,ii)
  enddo

end associate

if(IFlag0==1) return

! AB(1) PART
call AB_CAS_mithap(ABPLUS,ABMIN,EnDummy,URe,Occ,XOne,&
              IndN,IndX,IGem,NAct,INActive,NDimX,NBasis,NDimX,&
              NInte1,IntFileName,1d0,.true.)

! B = A.X
! C = X^T.B
allocate(workA(NDimX,NDimX))

call ABPM_TRAN(ABPLUS,workA,EBlock,EBlockIV,nblk,NDimX,.true.)
ABPLUS=workA
call ABPM_TRAN(ABMIN,workA,EBlock,EBlockIV,nblk,NDimX,.false.)
ABMIN=workA

deallocate(workA)

!! EBlock --> Y01Block
! block
! type(Y01BlockData) :: Y01Block(NDimX)
!
! do iblk=1,nblk
!   associate(B => Eblock(iblk))
!
!     do i=1,B%n
!        ipos = B%pos(i)
!        associate(Y => Y01Block(ipos))
!
!          Y%n = B%n
!          Y%l1 = B%l1
!          Y%l2 = B%l2
!          allocate(Y%vec0(Y%n))
!          Y%vec0(1:Y%n) = B%matY(i,1:Y%n)
!
!        end associate
!     enddo
!
!   end associate
! enddo
!
! associate(B => EblockIV)
! 
!   do i=1,B%n
!      ii = B%l1+i-1
!      ipos = B%pos(i)
!      associate(Y => Y01Block(ipos))
!        
!        Y%n = 1
!        Y%l1 = ii
!        Y%l2 = ii
!        allocate(Y%vec0(1))
!        Y%vec0(1) = 1d0/sqrt(2d0)
!      end associate
!   enddo
! 
! end associate
!
! ! dump to a file
! open(newunit=iunit,file=y01file,form='unformatted')
! do i=1,NDimX
!    associate(Y => Y01Block(i))
!      write(iunit) i, Y%n, Y%l1, Y%l2
!      write(iunit) Y%vec0
!    end associate
! enddo
! close(iunit)
!
! do i=1,NDimX
!    associate(Y => Y01Block(i))
!      deallocate(Y%vec0)
!    end associate
! enddo
!
! end block

 allocate(EigY(NDimX,NDimX),Eig(NDimX))
 if(IFlag0==0) allocate(EigY1(NDimX,NDimX),Eig1(NDimX),workA(NDimX,NDimX))
 
 EigY = 0
 Eig  = 0

! unpack Eig (1)
do iblk=1,nblk
   associate(B => Eblock(iblk))

     do i=1,B%n
        ipos = B%pos(i)
        EigY(ipos,B%l1:B%l2) = B%matY(i,1:B%n)
        if(IFlag0==0) EigY1(ipos,B%l1:B%l2) = B%matX(i,1:B%n)
     enddo
     Eig(B%l1:B%l2) = B%vec(1:B%n)

     deallocate(B%matY,B%matX,B%vec)
     deallocate(B%pos)

   end associate
enddo
!unpack Eig (2, IV part)
associate(B => EblockIV)

  do i=1,B%n
     ii = B%l1+i-1
     ipos = B%pos(i)
     EigY(ipos,ii) = 1d0/sqrt(2d0)
     if(IFlag0==0) EigY1(ipos,ii) = 1d0/sqrt(2d0)
     Eig(ii) = B%vec(i)
  enddo

  deallocate(B%vec)
  deallocate(B%pos)

end associate

if(IFlag0==0) then

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
   !!EigY1 = 0
   !workA=0
   !do j=1,NDimX
   !   if(Eig(j)/=0d0) then
   !      do i=1,NDimX
   !         if(Eig(i)/=0d0) then
   !            workA(i,j) = (ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
   !            if(Abs(Eig(i)-Eig(j))>Thresh) then
   !               workA(i,j) = workA(i,j) + (ABPLUS(i,j)+ABMIN(i,j))/(Eig(j)-Eig(i))
   !            endif
   !         endif
   !      enddo
   !   endif
   !enddo
   !
   !call dgemm('N','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,workA,NDimX,0d0,EigY1,NDimX)

   deallocate(workA)
endif

! dump response to a file!
open(newunit=iunit,file=propfile0,form='unformatted')
write(iunit) EigY
write(iunit) Eig
close(iunit)
if(IFlag0==0) then

   open(newunit=iunit,file=propfile1,form='unformatted')
   write(iunit) EigY1
   write(iunit) Eig1
   close(iunit)
endif

deallocate(EigY,Eig)
if(IFlag0==0) deallocate(EigY1,Eig1)

!! deallocate blocks
!do iblk=1,nblk
!   associate(B => Eblock(iblk))
!
!     deallocate(B%matY,B%matX,B%vec)
!     deallocate(B%pos)
!
!   end associate
!enddo
!! (IV part)
!associate(B => EblockIV)
!
!  deallocate(B%vec)
!  deallocate(B%pos)
!
!end associate

!! for AC0Corr
!EigY1 = 0
!do j=1,NDimX
!   if(Eig(j)/=0d0) then
!      do i=1,NDimX
!         if(Eig(i)/=0d0) then
!            val = 2d0*(ABPLUS(i,j)-ABMIN(i,j))/(Eig(i)+Eig(j))
!            do ii=1,NDimX
!               EigY1(ii,j) = EigY1(ii,j) + val*EigY(ii,i)
!            enddo
!         endif
!      enddo
!   endif
!enddo
!
!! second loop ...
!call dgemm('T','N',NDimX,NDimX,NDimX,1d0,EigY,NDimX,EigY1,NDimX,0d0,ABPLUS,NDimX)
!
!! energy loop
!open(newunit=iunit,file=trim(IntKFile),status='OLD', &
!     access='DIRECT',recl=8*NBasis*NOccup)
!!
!kl = 0
!do k=1,NOccup
!   do l=1,NBasis
!      kl = kl + 1
!      if(pos(l,k)/=0) then
!        irs = pos(l,k)
!        ir = l
!        is = k
!        read(iunit,rec=kl) work(1:NBasis*NOccup)
!        do j=1,NOccup
!           do i=1,NBasis
!              ints(i,j) = work((j-1)*NBasis+i)
!           enddo
!        enddo
!        ints(:,NOccup+1:NBasis) = 0
!
!        do j=1,NBasis
!           do i=1,j
!              if(pos(j,i)/=0) then
!                ipq = pos(j,i)
!                ip = j
!                iq = i
!                Crs = CICoef(l)+CICoef(k)
!                Cpq = CICoef(j)+CICoef(i)
!!
!                Aux = Crs*Cpq*ABPLUS(ipq,irs)
!                EAll = EAll + Aux*ints(j,i)
!
!                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))==1) EIntra = EIntra + Aux*ints(j,i)
!
!              endif
!           enddo
!        enddo
!
!      endif
!   enddo
!enddo
!
!close(iunit)

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



subroutine ACABMAT0_mithap(AMAT,BMAT,URe,Occ,XOne,&
               IndN,IndX,IGem,C, &
               NBasis,NDim,NDimX,NInte1,NGem,IntFileName,ISAPT,ACAlpha,IFlag)
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

integer,intent(in) :: NBasis,NDim,NDimX,NInte1,ISAPT,NGem
character(*) :: IntFileName
double precision,intent(out) :: AMAT(NDimX,NDimX),BMAT(NDimX,NDimX)
double precision,intent(in)  :: URe(NBasis,NBasis),Occ(NBasis),XOne(NInte1),C(NBasis)
double precision,intent(in)  :: ACAlpha
integer,intent(in) :: IndN(2,NDim),IndX(NDim)
integer,intent(in) :: IGem(NBasis),IFlag

integer :: i,j,k,l,ij,kl,kk,ll,klround
integer :: ip,iq,ir,is,it,iu,iw,ipq,irs,ICol,IRow
integer :: iunit,ios
integer :: IGemType
integer :: Ind(NBasis),AuxInd(3,3),pos(NBasis,NBasis)
double precision :: HNO(NBasis,NBasis),HNOCoef
double precision :: AuxCoeff(NGem,NGem,NGem,NGem),AuxVal,val
double precision :: valPQ,valRS
double precision :: OccProd(NBasis,NBasis),CProd(NBasis,NBasis)
double precision :: AuxH(NBasis,NBasis,NGem),AuxXC(NBasis,NBasis,NGem)
double precision :: SaveA,SaveB
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)
double precision,parameter :: Delta = 1.d-6

if(ISAPT==1) then
   write(6,'(1x,a)') 'Computing response-my'
else
   write(6,'(/,X,"***** COMPUTING AMAT, BMAT IN ACABMAT0: FFFF *****",/)')
endif

AMAT = 0
BMAT = 0

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
               !! Kasia's way
               !if(abs(abs(C(ip))-abs(C(iq)))<=Delta*abs(C(iq)).or.& 
               !     abs(abs(C(ir))-abs(C(is)))<=Delta*abs(C(ir))) then
               !
               !   AMAT(irs,ipq) = 0
               !   BMAT(irs,ipq) = 0

               !else
               !   
               !   SaveA = AMAT(irs,ipq)
               !   SaveB = BMAT(irs,ipq)
               !   
               !   val = (C(ip) + C(iq))*(C(ir) + C(is))
               !   AMAT(irs,ipq) = (SaveA+SaveB)/val
               !   val = (C(ip) - C(iq))*(C(ir) - C(is))
               !   BMAT(irs,ipq) = (SaveA-SaveB)/val

               ! our way
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

               !endif ! KP
            endif
                        
         enddo
      endif
   enddo

print*, "AB-my",norm2(AMAT),norm2(BMAT)

!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))
!
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

subroutine ModABMin_mithap(Occ,SRKer,Wt,OrbGrid,ABMin,IndN,IndX,NDimX,NGrid,NBasis,&
                           twofile,twoerfile)
! ADD CONTRIBUTIONS FROM THE srALDA KERNEL TO AB MATRICES
implicit none

integer,parameter :: maxlen = 128
integer,intent(in) :: NBasis,NDimX,NGrid
integer,intent(in) :: IndN(2,NDimX),IndX(NDimX)
character(*),intent(in) :: twofile,twoerfile
double precision,intent(in) :: Occ(NBasis),SRKer(NGrid), &
                               Wt(NGrid),OrbGrid(NGrid,NBasis)
double precision,intent(inout) :: ABMin(NDimX,NDimX)

integer :: offset,batchlen,iunit1,iunit2
integer :: i,j,k,l,kl,ip,iq,ir,is,irs,ipq,igrd
integer :: IRow,ICol
double precision :: XKer1234,TwoSR,Cpq,Crs
integer :: pos(NBasis,NBasis)
double precision :: CICoef(NBasis)
double precision,allocatable :: work1(:),work2(:),WtKer(:)
double precision,allocatable :: batch(:,:),ABKer(:,:)
double precision,allocatable :: ints1(:,:),ints2(:,:)

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

allocate(work1(NBasis**2),work2(NBasis**2),&
         ints1(NBasis,NBasis),ints2(NBasis,NBasis),&
         WtKer(maxlen),batch(maxlen,NBasis),ABKer(NDimX,NDimX))

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

ABKer = 0

!print*, 'ModABMin_mithap'

do offset=0,NGrid,maxlen
   batchlen = min(NGrid-offset,maxlen)
   if(batchlen==0) exit

   WtKer(1:batchlen) = Wt(offset+1:offset+batchlen)*SRKer(offset+1:offset+batchlen)
   batch(1:batchlen,1:NBasis) = OrbGrid(offset+1:offset+batchlen,1:NBasis)

   do IRow=1,NDimX
      ip = IndN(1,IRow)
      iq = IndN(2,IRow)
      ipq = IndX(IRow)
   
      do ICol=1,NDimX
         ir=IndN(1,ICol)
         is=IndN(2,ICol)
         irs=IndX(ICol)
         if(irs.gt.ipq) cycle
    
         XKer1234 = 0
         do i=1,batchlen
            XKer1234 = XKer1234 + WtKer(i)* &
            batch(i,ip)*batch(i,iq)*batch(i,ir)*batch(i,is)
         enddo
         
         ABKer(ipq,irs) = ABKer(ipq,irs) + XKer1234
         ABKer(irs,ipq) = ABKer(ipq,irs)
      
      enddo
   enddo

enddo

open(newunit=iunit1,file=trim(twofile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
open(newunit=iunit2,file=trim(twoerfile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1   
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        read(iunit1,rec=kl) work1(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work1,ints1,NBasis)
        read(iunit2,rec=kl) work2(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work2,ints2,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)
                !if(irs.gt.ipq) cycle

                TwoSR = ints1(i,j)-ints2(i,j)

                ABMIN(ipq,irs) = ABMIN(ipq,irs) & 
                               + 4.0d0*Cpq*Crs*(TwoSR+ABKer(ipq,irs))
                !ABMIN(irs,ipq) = ABMIN(ipq,irs) 

              endif  
           enddo
        enddo

      endif 
   enddo
enddo 

close(iunit1)
close(iunit2)

deallocate(ABKer,batch,WtKer,ints2,ints1,work2,work1)

end subroutine ModABMin_mithap

subroutine FockGen_mithap(Fock,OneRdm,XOne,NInte1,NBasis,IntFileName)
!
!     GENERALIZED FOCK MATRIX
!
implicit none

type(AOReaderData) :: reader

integer,intent(in) :: NInte1,NBasis
character(*) :: IntFileName
double precision,intent(in)  :: OneRdm(NInte1),XOne(NInte1)
double precision,intent(out) :: Fock(Ninte1)
integer :: iunit,kl,k,l
logical :: empty
double precision,allocatable :: OneRdmSq(:,:),FockSq(:,:),ints(:,:),work1(:)


allocate(OneRdmSq(NBasis,NBasis),FockSq(NBasis,NBasis),ints(NBasis,NBasis),work1(NInte1))

call triang_to_sq2(OneRdm,OneRdmSq,NBasis)
call triang_to_sq2(XOne,FockSq,NBasis)

call reader%open(trim(IntFileName))

kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      !read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call reader%getTR(kl,work1,empty)
      if(empty) cycle
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

!close(iunit)
call reader%close

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

subroutine PotCoul_mithap(VHSR,OneRdm,doRSH,aoerfile,NBasis)
!
!     RETURNS SR HARTREE (COULOMB) POTENTIAL IN AN MO MATRIX REPRESENTATION 
!
implicit none

integer,intent(in) :: NBasis
logical,intent(in) :: doRSH
character(*),intent(in) :: aoerfile
double precision,intent(in) :: OneRdm(NBasis*(NBasis+1)/2)
double precision :: VHSR(NBasis*(NBasis+1)/2)
double precision,allocatable :: J(:,:),Jerf(:,:),work(:,:)

allocate(Jerf(NBasis,NBasis),J(NBasis,NBasis),work(NBasis,NBasis))

work = 0
call triang_to_sq2(OneRdm,work,NBasis)

if(doRSH) then
  call make_J1(NBasis,work,J,'AOTWOSORT')
  call make_J1(NBasis,work,Jerf,aoerfile)
  work = 2.0d0*J - 2.0d0*Jerf
else
  call make_J1(NBasis,work,J,'AOTWOSORT')
  work = 2.0d0*J
endif
call sq_to_triang2(work,VHSR,NBasis) 

deallocate(work,J,Jerf)

end subroutine PotCoul_mithap

subroutine EneGVB_FFFF(ETot,URe,Occ,C,XOne, &
                        IGem,IndN,NBasis,NInte1,IntFileName,NDimX,NGem)
implicit none

integer,intent(in) :: NBasis,NInte1,NDimX,NGem
character(*) :: IntFileName
double precision :: ETot

integer :: IndN(2,NDimX),IGem(NBasis)
double precision :: URe(NBasis,NBasis),Occ(NBasis),C(NBasis)
double precision :: XOne(NInte1)

integer :: i,j,ii,ij,ia,ib,iab
integer :: k,l,kl,kk,ll,klround
integer :: iunit
integer,external :: NAddr3
double precision :: FacIJ,FacIK
double precision :: EOne, EIntraGem, EInterCoul, EInterExch

integer :: pos(NBasis,NBasis)
double precision :: HNO(NBasis,NBasis)
double precision,allocatable :: work1(:),work2(:)
double precision,allocatable :: ints(:,:)

write(lout,'(/A)') ' Electronic GVB energy check:'

allocate(work1(NBasis**2),work2(NBasis**2),ints(NBasis,NBasis))

call triang_to_sq(XOne,work1,NBasis)
call dgemm('N','N',NBasis,NBasis,NBasis,1d0,URe,NBasis,work1,NBasis,0d0,work2,NBasis)
call dgemm('N','T',NBasis,NBasis,NBasis,1d0,work2,NBasis,URe,NBasis,0d0,HNO,NBasis)
call sq_symmetrize(HNO,NBasis)

!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-my', norm2(work1(1:NBasis*(NBasis+1)/2))
!HNO=transpose(HNO)
!call sq_to_triang2(HNO,work1,NBasis)
!write(LOUT,*) 'HNO-tr', norm2(work1(1:NBasis*(NBasis+1)/2))

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

EOne = 0d0
EIntraGem  = 0d0
EInterCoul = 0d0
EInterExch = 0d0

do i=1,NBasis
   EOne = EOne + 2d0*Occ(i)*HNO(i,i)
enddo

open(newunit=iunit,file=trim(IntFileName),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)

kl = 0
do ll=1,NBasis
   do kk=1,ll
      kl = kl + 1
      read(iunit,rec=kl) work1(1:NBasis*(NBasis+1)/2)
      call triang_to_sq2(work1,ints,NBasis)
      k = kk
      l = ll

      FacIJ = 2d0
      if(k==l) FacIJ = 1d0
      if(IGem(k)==IGem(l)) then
         EIntraGem = EIntraGem + FacIJ*C(k)*C(l)*ints(k,l)
      else
        EInterExch = EInterExch - FacIJ*Occ(k)*Occ(l)*ints(k,l)
      endif

      FacIK = 2d0
      if(k==l) then
         do i=k,NBasis
            if(IGem(i)/=IGem(k)) then
            if(i==k) FacIK = 1d0
            EInterCoul = EInterCoul + 2d0*FacIK*Occ(i)*Occ(k)*ints(i,i)
            endif
         enddo
      endif

   enddo
enddo

close(iunit)

write(LOUT,'(" One-electron energy",25X,F17.8)') EOne
write(LOUT,'(" GVB intra-gem electron interaction",10X,F17.8)') EIntraGem
write(LOUT,'(" GVB inter-gem Coulomb interaction",11X,F17.8)')  EInterCoul
write(LOUT,'(" GVB inter-gem exchange interaction",10X,F17.8)') EInterExch
write(LOUT,'(" Total GVB",34X,F18.8)') EOne+EIntraGem+EInterCoul+EInterExch

ETot = EOne + EIntraGem + EInterCoul + EInterExch

deallocate(ints,work2,work1)

end subroutine EneGVB_FFFF

subroutine EERPA_FFFF(ECorr,EVec,EVal,Occ,CICoef,IGem,   &
                      IAuxGem,IG1,IG2,IG3,IG4,IB, &
                      IndN,NDimX,NBasis,IntFile,IFrag)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IG1,IG2,IG3,IG4,IB 
integer,intent(in) :: IGem(NBasis),IAuxGem(NBasis),IndN(2,NDimX)
integer,intent(in) :: IFrag
character(*),intent(in) :: IntFile
double precision,intent(out) :: ECorr
double precision,intent(in) :: Occ(NBasis)
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
integer :: OccProd(NBasis,NBasis),AuxIG(NBasis)
integer :: IFlPQRS,IFlP,IFlQ,IFlR,IFlS
integer :: NoVirtP,NoVirtQ,NoVirtR,NoVirtS,NoVirt
integer :: IBdy,IBdyG,ICond
double precision :: CICoef(NBasis),Cpq,Crs,SumY,Aux
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 5.d2
integer :: itmp1,itmp2

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

OccProd = 0
do j=1,NBasis
   do i=1,NBasis
      if(Occ(i)*Occ(j)/=0d0) OccProd(i,j) = 1
   enddo
enddo

AuxIG = 0
do i=1,NBasis
   if(IGem(i)==IG1.or.IGem(i)==IG2.or.IGem(i)==IG3.or.IGem(i)==IG4) AuxIG(i) = 1
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0

itmp1 = 0
itmp2 = 0
! FULL INTS
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work,ints,NBasis)

        !IFlS=0
        !If(IGem(IS).Eq.IG1.Or.IGem(IS).Eq.IG2.Or.IGem(IS).Eq.IG3.Or.IGem(IS).Eq.IG4) IFlS=1
        !NoVirtS=1
        !If(Occ(IS).Eq.0d0) NoVirtS=0

        !IFlR=0
        !If(IGem(IR).Eq.IG1.Or.IGem(IR).Eq.IG2.Or.IGem(IR).Eq.IG3.Or.IGem(IR).Eq.IG4) IFlR=1
        !NoVirtR=1
        !If(Occ(IR).Eq.0d0) NoVirtR=0

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i

                !IFlP=0
                !If(IGem(ip).Eq.IG1.Or.IGem(ip).Eq.IG2.Or.IGem(ip).Eq.IG3.Or.IGem(ip).Eq.IG4) IFlP=1
                !NoVirtP=1
                !If(Occ(IP).Eq.0d0) NoVirtP=0

                !IFlQ=0
                !If(IGem(iq).Eq.IG1.Or.IGem(iq).Eq.IG2.Or.IGem(iq).Eq.IG3.Or.IGem(iq).Eq.IG4) IFlQ=1
                !NoVirtQ=1
                !If(Occ(IQ).Eq.0d0) NoVirtQ=0

                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)

                NoVirt  = 0
                IFlPQRS = 0
                if(OccProd(ir,is)*OccProd(ip,iq)/=0) NoVirt = 1
                if(AuxIG(ir)*AuxIG(is)*AuxIG(ip)*AuxIG(iq)/=0) IFlPQRS = 1
               
                call IBody(IBdy,IGem(ip),IGem(iq),IGem(ir),IGem(is))

                ICond = 0
                if(IFlPQRS==1.and.IBdy.eq.IB) ICond = 1

                !! FOR GVB ONLY
                if(IB.gt.2.and.NoVirt.Eq.1.and. &
                    IBdy.eq.IB-1.and.IFlPQRS==1) ICond = 1

                ! FOR GVB ONLY: IF IFrag=1,IB=2 (ONE-BODY), ALLOW ALL CASES 
                ! EXCEPT WHEN ALL ORBITALS ARE FROM THE SAME GEMINAL
                if(IFrag.Eq.1.and.IB.Eq.2.and.IFlPQRS==1) then
                   call IBody(IBdyG,IAuxGem(ip),IAuxGem(iq),IAuxGem(ir),IAuxGem(is))
                   If(IBdyG.Ne.1) ICond = 1
                endif

                if(ICond.Eq.1) then

                   ISkippedEig = 0
                   SumY = 0
                   do kk=1,NDimX
                      if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                         SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                      else
                         ISkippedEig = ISkippedEig + 1
                         Skipped(ISkippedEig) = EVal(kk)
                      endif
                   enddo

                   Aux = Crs*Cpq*SumY
 
                   if(iq.Eq.is.and.ip.Eq.ir) then
                      Aux = Aux - 0.5d0*(Occ(ip)*(1d0-Occ(is))+Occ(is)*(1d0-Occ(ip)))
                   endif  

                   ECorr = ECorr + Aux*ints(j,i)

                endif ! ICond

              endif
           enddo
        enddo

      endif

   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

deallocate(Skipped)
deallocate(ints,work)

end subroutine EERPA_FFFF

subroutine EneERPA_FFFF(ETot,ECorr,ENuc,EVec,EVal,Occ,CICoef,IGem,   &
                        IndN,NDimX,NBasis,IntFile)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IGem(NBasis),IndN(2,NDimX)
character(*),intent(in) :: IntFile
double precision,intent(in)    :: ENuc
double precision,intent(inout) :: ETot,ECorr
double precision,intent(in) :: CICoef(NBasis),Occ(NBasis)
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
integer :: OccProd(NBasis,NBasis)
double precision :: Cpq,Crs,SumY,Aux,EIntra
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 5.d2

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = i
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0
EIntra= 0

! FULL INTS
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work,ints,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i

                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)

                ISkippedEig = 0
                SumY = 0
                do kk=1,NDimX
                   if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                      SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                   else
                      ISkippedEig = ISkippedEig + 1
                      Skipped(ISkippedEig) = EVal(kk)
                   endif
                enddo

                Aux = 2d0*Crs*Cpq*SumY

                if(iq.Eq.is.and.ip.Eq.ir) then
                   Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                endif

                ECorr = ECorr + Aux*ints(j,i)

                if(IGem(ip).eq.IGem(iq).and.IGem(ir).eq.IGem(is).and.IGem(ip).eq.IGem(ir)) then
                   EIntra = EIntra + Aux*ints(j,i)
                endif

              endif
           enddo
        enddo

      endif
   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

ECorr = 0.5d0*(ECorr-EIntra)
write(lout,'(1x,a,3f15.8)') 'EGVB+ENuc, Corr, ERPA-GVB',ETot+ENuc,ECorr,ETot+ENuc+ECorr

deallocate(Skipped)
deallocate(ints,work)

end subroutine EneERPA_FFFF

subroutine ACEneERPA_FFFF(ECorr,EVec,EVal,Occ,IGem, &
                          IndN,IndX,NOccup,NDimX,NBasis,IntFile)
implicit none

integer,intent(in) :: NDimX,NBasis
integer,intent(in) :: IGem(NBasis),IndN(2,NDimX),IndX(NDimX)
integer,intent(in) :: NOccup
character(*),intent(in) :: IntFile
double precision,intent(out) :: ECorr
double precision,intent(in) :: EVec(NDimX,NDimX),EVal(NDimX)
double precision :: Occ(NBasis)

integer :: i,j,k,l,kl,kk,ip,iq,ir,is,ipq,irs
integer :: iunit,ISkippedEig
integer :: pos(NBasis,NBasis)
logical :: AuxCoeff(3,3,3,3)
double precision :: CICoef(NBasis),Cpq,Crs,SumY,Aux
double precision,allocatable :: work(:),ints(:,:),Skipped(:)
double precision,parameter :: SmallE = 1.d-3,BigE = 1.d8

do i=1,NBasis
   CICoef(i) = sign(sqrt(Occ(i)),Occ(i)-0.5d0)
enddo

pos = 0
do i=1,NDimX
   pos(IndN(1,i),IndN(2,i)) = IndX(i)
enddo

AuxCoeff = .true.
do l=1,3
   do k=1,3
      do j=1,3
         do i=1,3
            if((i==j).and.(j==k).and.(k==l)) then
               AuxCoeff(i,j,k,l) = .false.
            endif
         enddo
      enddo
   enddo
enddo

allocate(work(NBasis**2),ints(NBasis,NBasis),Skipped(NDimX))

ISkippedEig = 0
ECorr = 0

! FULL INTS
open(newunit=iunit,file=trim(IntFile),status='OLD', &
     access='DIRECT',recl=8*NBasis*(NBasis+1)/2)
kl = 0
do l=1,NBasis
   do k=1,l
      kl = kl + 1
      if(pos(l,k)/=0) then
        irs = pos(l,k)
        ir = l
        is = k
        read(iunit,rec=kl) work(1:NBasis*(NBasis+1)/2)
        call triang_to_sq2(work,ints,NBasis)

        do j=1,NBasis
           do i=1,j
              if(pos(j,i)/=0) then
                ipq = pos(j,i)
                ip = j
                iq = i
                Crs = CICoef(l)+CICoef(k)
                Cpq = CICoef(j)+CICoef(i)

                !if(.not.(IGem(ir).eq.IGem(is).and.IGem(ip).eq.IGem(iq)&
                !.and.IGem(ir).eq.IGem(ip)).and.ir.gt.is.and.ip.gt.iq) then
                if(AuxCoeff(IGem(ip),IGem(iq),IGem(ir),IGem(is))) then

                   ISkippedEig = 0
                   SumY = 0 
                   do kk=1,NDimX
                      if(EVal(kk).gt.SmallE.and.EVal(kk).lt.BigE) then
                         SumY = SumY + EVec(ipq,kk)*EVec(irs,kk)
                      else
                         ISkippedEig = ISkippedEig + 1
                         Skipped(ISkippedEig) = EVal(kk)
                      endif
                   enddo
          
                   Aux = 2*Crs*Cpq*SumY
          
                   if(iq.Eq.is.and.ip.Eq.ir) then
                      Aux = Aux - Occ(ip)*(1d0-Occ(is))-Occ(is)*(1d0-Occ(ip))
                   endif  

                   ECorr = ECorr + Aux*ints(j,i)
          
                ! endinf of If(IP.Gt.IR.And.IQ.Gt.IS) 
                endif

              endif
           enddo
        enddo

      endif

   enddo
enddo

close(iunit)

if(ISkippedEig/=0) then
  write(LOUT,'(/,1x,"The number of discarded eigenvalues is",i4)') &
       ISkippedEig
  do i=1,ISkippedEig
     write(LOUT,'(1x,a,i4,f15.8)') 'Skipped',i,Skipped(i)
  enddo
endif

deallocate(Skipped)
deallocate(ints,work)

end subroutine ACEneERPA_FFFF



end module abmat
