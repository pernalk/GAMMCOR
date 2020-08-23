*fordeck dftfun_ecerfpbe $Revision: 2006.3 $
      subroutine dftfun_ecerfpbe(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,mu)
C************************************************************************
C     Short-range PBE correlation energy functional for erf interaction
C
C
C************************************************************************
      implicit none
c      include "common/chirs"
      include "debug"
      include "constants.h"


! input
      character*(*) name
      logical fderiv,open
      integer igrad,npt
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)

! output
      double precision zk(*)
      double precision vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)

! local
      double precision tol
      parameter(tol=1d-12)

      character*(30) namedummy

      double precision zkcerflda(npt)
      double precision vrhoccerflda(npt)
      double precision vrhoocerflda(npt)

      double precision zkclda(npt)
      double precision vrhocclda(npt)
      double precision vrhooclda(npt)

      integer i
      double precision mu
      double precision rho,drho2,rhoa,rhob
      double precision ecerflda,decerfldadrho
      double precision eclda,decldadrho
      double precision ecerfpbe,decerfpbedrho,decerfpbedrhoo
      double precision decerfpbeddrho2
      double precision arglog,arglogs,arglogss,alpha,beta,betas,gamma
      double precision Aa,Ab,Ac,Aas,tq,tqs,tqss,decerfpur,decpur
      double precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision t11,t12,t13,t14,t15,t16,t17,t18,t19
      double precision zeta,phi,phi2,phi3,phi4,phis,arglogsc

c herer!!!
      do i=1,npt
      zk(i)=z0   
      end do
C debug beg
c      here='dftfun_ecerfpbe'
c      ldebug = index(debuglist,trim(here)) .ne. 0
c     .    .or. index(debuglist,'all') .ne. 0
c      if(ldebug) print*,trim(here),': entering'
C debug end
      ldebug=.False. 

C Description of the functional for output of molpro
      name='SR PBE correlation erf'

C Compute exchange and correlation LDA part
      do 10 i=1,npt
         zkcerflda(i)=z0
         vrhoccerflda(i)=z0
         vrhoocerflda(i)=z0
         zkclda(i)=z0
         vrhocclda(i)=z0
         vrhooclda(i)=z0
 10   continue

      call dftfun_ecerf(namedummy,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   zkcerflda,vrhoccerflda,vrhoocerflda,mu)

      call dftacg_pw92c(namedummy,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   zkclda,vrhocclda,vrhooclda,
     >                   vsigmacc,vsigmaco,vsigmaoo)

C First-type gradient functional
      igrad=1

      alpha=2.78d0
      gamma=3.1091d-2

      if(.not.open) then

C     Loop over the grid
      do 1 i=1,npt

!        test on density
         if (dabs(rhoc(i)).lt.tol) goto 1

!        Density
         rho = rhoc(i)

!        Square of density gradient
         drho2 = sigmacc(i)

         if(ldebug) write(*,*)"rho=",rho
         if(ldebug) write(*,*)"drho2=",drho2

!        LDA energy density
         ecerflda = zkcerflda(i)
         eclda = zkclda(i)

         if ((ecerflda/eclda).le.z0) then
            beta=z0
         else
            beta=6.6725d-2*(ecerflda/eclda)**alpha
         endif
         tq=drho2*6.346820607d-2*rho**(-z7/z3)
         Ab=dexp(-ecerflda/(rho*gamma))-z1
         if (dabs(Ab).le.dabs(beta*tol)) then
            ecerfpbe=ecerflda
         else
            Aa=beta/(gamma*Ab)
            if (Aa.lt.tol) Aa=tol
            Ac=z1+Aa*tq+Aa**2*tq**2
            arglog=z1+beta*(z1-z1/Ac)/(gamma*Aa)
            ecerfpbe=ecerflda+rho*gamma*dlog(arglog)
         end if

         zk(i) = ecerfpbe

         if(ldebug) write(*,*)"ecerfpbe=",ecerfpbe

C Derivative

      if (fderiv) then

!        LDA energy density derivative
         decerfldadrho = vrhoccerflda(i)
         decldadrho = vrhocclda(i)

         decerfpur=(decerfldadrho-ecerflda/rho)/rho
         decpur=(decldadrho-eclda/rho)/rho
         betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
         if (dabs(Ab).le.dabs(beta*tol)) then
            decerfpbedrho=decerfldadrho
         else
            Aas=betas/(gamma*Ab)+Aa*(z1+z1/Ab)*decerfpur/gamma
            tqs=-z7*tq/(z3*rho)
            arglogs=betas*tq*(z1+Aa*tq)/(Ac*gamma)
     &       +beta*tqs*(z1+Aa*tq)/(Ac*gamma)
     &       -beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(z2+Aa*tq)/(Ac**2*gamma)
            decerfpbedrho=decerfldadrho+gamma*dlog(arglog)
     &       +gamma*rho*arglogs/arglog
         end if

         if (dabs(Ab).le.dabs(beta*tol)) then
            decerfpbeddrho2=0.0d0
         else
            arglogsc=Ab*(Aa+z2*Aa*Aa*tq)/(Ac*Ac)
            tqss=6.346820607d-2*rho**(-z7/z3)
            arglogss=tqss*arglogsc
            decerfpbeddrho2=rho*gamma*arglogss/arglog
         end if

C        derivatives
         vrhoc(i) = vrhoc(i) + decerfpbedrho
         vsigmacc(i) = vsigmacc(i) + decerfpbeddrho2

         if(ldebug) write(*,*)" decerfpbedrho=",decerfpbedrho
         if(ldebug) write(*,*)" decerfpbeddrho2=",decerfpbeddrho2

C     End of derivative
      end if

C     End of loop over the grid
1     continue

      else
      do 2 i=1,npt

!        test on density
         if (dabs(rhoc(i)).lt.tol) goto 2

!        Density
         rho = rhoc(i)

!        Square of density gradient
         drho2 = sigmacc(i)

!        Spin polarisation
         rhoa=max((rhoc(i)+rhoo(i))*.5d0,1.0d-15)
         rhob=max((rhoc(i)-rhoo(i))*.5d0,1.0d-15)
         zeta = (rhoa-rhob)/(rhoa+rhob)

         if(ldebug) write(*,*)"rhoc=",rho
         if(ldebug) write(*,*)"drho2=",drho2
         if(ldebug) write(*,*)"zeta=",zeta

!        LDA energy density
         ecerflda = zkcerflda(i)
         eclda = zkclda(i)

         if ((ecerflda/eclda).le.z0) then
            beta=z0
         else
            beta=6.6725d-2*(ecerflda/eclda)**alpha
         endif
         phi=((z1+zeta)**(z2/z3)+(z1-zeta)**(z2/z3))/z2
         phi2=phi*phi
         phi3=phi2*phi
         phi4=phi3*phi
         tq=drho2*6.346820607d-2*rho**(-z7/z3)/phi2
         Ab=dexp(-ecerflda/(rho*gamma*phi3))-z1
         if (dabs(Ab).le.dabs(beta*tol)) then
            ecerfpbe=ecerflda
         else
            Aa=beta/(gamma*Ab)
            Ac=z1+Aa*tq+Aa**2*tq**2
            if (Aa.lt.tol) Aa=tol
            arglog=z1+beta*(z1-z1/Ac)/(gamma*Aa)
            ecerfpbe=ecerflda+rho*phi3*gamma*dlog(arglog)
         end if

         zk(i) = ecerfpbe

         if(ldebug) write(*,*)"ecerfpbe=",ecerfpbe

C Derivative

      if (fderiv) then

!        LDA energy density derivative
         decerfldadrho = vrhoccerflda(i)
         decldadrho = vrhocclda(i)

         decerfpur=(decerfldadrho-ecerflda/rho)/rho
         decpur=(decldadrho-eclda/rho)/rho
         betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
         phis=((rhoa - rhob)
     &     *((rhoa/(rhoa + rhob))**f13 - (rhob/(rhoa + rhob))**f13))
     &     /(3.*2**f13*(rhoa/(rhoa + rhob))**f13
     &     *(rhob/(rhoa + rhob))**f13
     &     *(rhoa + rhob)**2)
         if (dabs(Ab).le.dabs(beta*tol)) then
            decerfpbedrho=decerfldadrho
         else
            Aas=betas/(gamma*Ab)+Aa*(z1+z1/Ab)*
     &       (decerfpur/phi3-z3*phis*ecerflda/(rho*phi4))/gamma
            tqs=-z7*tq/(z3*rho)-z2*tq*phis/phi
            arglogs=betas*tq*(z1+Aa*tq)/(Ac*gamma)
     &       +beta*tqs*(z1+Aa*tq)/(Ac*gamma)
     &       -beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(z2+Aa*tq)/(Ac**2*gamma)
            decerfpbedrho=decerfldadrho+gamma*(phi3*dlog(arglog)
     &       +z3*rho*phis*phi2*dlog(arglog)+rho*phi3*arglogs/arglog)
         end if

c     Alternative calculation of arglog (use beta=Aa*Ab*gamma)
c        arglogsa=Ab*(tq+z2*Aa*tq*tq)/(Ac*Ac)
c        arglogsb=-(z1-z1/Ac)*(Ab+z1)/(gamma)
c        arglogsc=Ab*(Aa+z2*Aa*Aa*tq)/(Ac*Ac)
c        arglogs=arglogsa*Aas
c    &    +arglogsb*(decerfpur/phi3-z3*phis*ecerflda/(rho*phi4))
c    &    +arglogsc*tqs

         if (dabs(Ab).le.dabs(beta*tol)) then
            decerfpbeddrho2=0.0d0
         else
            arglogsc=Ab*(Aa+z2*Aa*Aa*tq)/(Ac*Ac)
            tqss=6.346820607d-2*rho**(-z7/z3)/phi2
            arglogss=tqss*arglogsc
            decerfpbeddrho2=rho*gamma*phi3*arglogss/arglog
         end if

!        LDA energy density derivative
         decerfldadrho = vrhoocerflda(i)
         decldadrho = vrhooclda(i)

         decerfpur=decerfldadrho/rho
         decpur=decldadrho/rho
         betas=alpha*beta*(decerfpur*rho/ecerflda-decpur*rho/eclda)
         phis=(rhob*(rhoa/(rhoa + rhob))**(2*f13)
     &     -rhoa*(rhob/(rhoa + rhob))**(2*f13))
     &     /(3.*2**f13*rhoa*rhob)

         if (dabs(Ab).le.dabs(beta*tol)) then
            decerfpbedrhoo=decerfldadrho
         else
            Aas=betas/(gamma*Ab)+Aa*(z1+z1/Ab)*
     &       (decerfpur/phi3-z3*phis*ecerflda/(rho*phi4))/gamma
            tqs=-z2*tq*phis/phi
            arglogs=betas*tq*(z1+Aa*tq)/(Ac*gamma)
     &       +beta*tqs*(z1+Aa*tq)/(Ac*gamma)
     &       -beta*tq*Aa*tq*(Aas*tq+Aa*tqs)*(z2+Aa*tq)/(Ac**2*gamma)
            decerfpbedrhoo=decerfldadrho+gamma*
     &       (z3*rho*phis*phi2*dlog(arglog)+rho*phi3*arglogs/arglog)
         end if

C        derivatives
         vrhoc(i) = vrhoc(i) + decerfpbedrho
         vrhoo(i) = vrhoo(i) + decerfpbedrhoo
         vsigmacc(i) = vsigmacc(i) + decerfpbeddrho2

         if(ldebug) write(*,*)" decerfpbedrho=",decerfpbedrho
         if(ldebug) write(*,*)" decerfpbedrhoo=",decerfpbedrhoo
         if(ldebug) write(*,*)" decerfpbeddrho2=",decerfpbeddrho2

C     End of derivative
      end if

C     End of loop over the grid
2     continue

      end if

C debug beg
      if(ldebug) then
       print*,trim(here),': exiting'
      end if
C debug end

      return
      end
