*fordeck dftfun_exerfpbe $Revision: 2006.0 $
C************************************************************************
      subroutine dftfun_exerfpbe(name,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   sigmacc,sigmaco,sigmaoo,
     >                   zk,vrhoc,vrhoo,
     >                   vsigmacc,vsigmaco,vsigmaoo,mu)
C************************************************************************
C     Short-range PBE exchange energy functional for erf interaction
C
C
C************************************************************************
      implicit none
c      include "common/chirs"
      include "debug"
      include "constants.h"
      include "big"

! input
      character*(*) name
      logical fderiv,open
      integer igrad,npt
      double precision rhoc(*),rhoo(*)
      double precision sigmacc(*),sigmaco(*),sigmaoo(*)
      double precision mu

! output
      double precision zk(*)
      double precision vrhoc(*),vrhoo(*)
      double precision vsigmacc(*),vsigmaco(*),vsigmaoo(*)

! function
      double precision berf
      double precision dberfda

! local
      double precision tol
      parameter(tol=1d-12)

      character*(30) namedummy

      double precision zkxerflda(npt)
      double precision vrhocxerflda(npt)
      double precision vrhooxerflda(npt)

      integer i
      integer izk, icorr
      double precision rho,drho2
      double precision exerflda,dexerfldadrho
      double precision exerfpbe,dexerfpbedrho,dexerfpbeddrho2
      double precision t1,t2,t3,t4
      double precision kappa,sq,sqs,sqss,fx,fxs,ksig
      double precision q

c herer!!!
      izk=0
c izk=icorr(npt)
      do i=1,npt
      zk(i)=z0
      end do
      
      ldebug=.False.
C debug beg
      here='dftfun_exerfpbe'
      ldebug = index(debuglist,trim(here)) .ne. 0
     .    .or. index(debuglist,'all') .ne. 0
      if(ldebug) print*,trim(here),': entering'
C debug end

C Description of the functional for output of molpro
      name='SR PBE exchange erf'

      if(.not.open) then

C Compute LDA part
      do i=1,npt
         zkxerflda(i)=z0
         vrhocxerflda(i)=z0
         vrhooxerflda(i)=z0
      end do

      call dftfun_exerf(namedummy,fderiv,open,igrad,npt,rhoc,rhoo,
     >                   zkxerflda,vrhocxerflda,vrhooxerflda,mu)

C     Loop over the grid
      do 1 i=1,npt

!        density
         rho = rhoc(i)
!        square of density gradient
         drho2 = sigmacc(i)

!        test on density
         if (rho.lt.tol) goto 1

         if(ldebug) write(*,*)"rho=",rho
         if(ldebug) write(*,*)"drho2=",drho2

!        LDA energy density
         exerflda = zkxerflda(i)

         kappa=0.804d0
         sq=drho2*2.6121172985233599567768d-2*rho**(-z8/z3)
         fx=z1+kappa
     &    -kappa/(z1+berf(1.616204596739954813d
     &    -1*mu*rho**(-f13))*sq/kappa)
         exerfpbe=exerflda*fx
c herer!!!
c         if(i.eq.1538) then
c         write(*,*)i,zkxerflda(i),fx
c         endif

         zk(i) = exerfpbe

         if(ldebug) write(*,*)"exerfpbe=",exerfpbe

C Derivative

         if (fderiv) then

!        LDA energy density derivative
         dexerfldadrho = vrhocxerflda(i)

         sqs=-z8*sq/(z3*rho)
         fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-z4*f13)/z3
     &    *dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq
     &    +berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)
     &    /(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
         dexerfpbedrho=dexerfldadrho*fx+exerflda*fxs
         sqss=2.6121172985233599567768d-2*rho**(-z8/z3)
         dexerfpbeddrho2=exerflda*berf(1.616204596739954813d-1*mu
     &    *rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa
     &    +berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

C      derivatives
       vrhoc(i) = vrhoc(i) + dexerfpbedrho
       vsigmacc(i) = vsigmacc(i) + dexerfpbeddrho2

       if(ldebug) write(*,*)" dexerfpbedrho=",dexerfpbedrho
       if(ldebug) write(*,*)" dexerfpbeddrho2=",dexerfpbeddrho2

C     End of derivative
      end if

C     End of loop over the grid
1     continue

      else ! else of (.not.open)

C Compute LDA part
      do i=1,npt
         zkxerflda(i)=z0
         vrhocxerflda(i)=z0
         vrhooxerflda(i)=z0
      end do

      do i=1,npt
         q(izk-1+i)=max((rhoc(i)+rhoo(i))*.5d0,0d0)*2.D0
      end do

      call dftfun_exerf(namedummy,fderiv,.false.,igrad,npt,q(izk),rhoo,
     >                   zkxerflda,vrhocxerflda,vrhooxerflda,mu)

C     Loop over the grid
      do 2 i=1,npt

!        Alpha density
         rho=max((rhoc(i)+rhoo(i))*.5d0,0d0)*2.D0
!        Square of alpha density gradient
         drho2=.25d0*sigmacc(i)+.25d0*sigmaoo(i)+.5d0*sigmaco(i)
         drho2=max(drho2,0d0)*4.0d0

!        test on density
         if (rho.lt.tol) goto 2

!        LDA energy density
         exerflda = zkxerflda(i)

         kappa=0.804d0
         sq=drho2*2.6121172985233599567768d-2*rho**(-z8/z3)
         fx=z1+kappa
     &    -kappa/(z1+berf(1.616204596739954813d
     &    -1*mu*rho**(-f13))*sq/kappa)
         exerfpbe=exerflda*fx

         zk(i) = exerfpbe*0.5d0

C Derivative

         if (fderiv) then

!        LDA energy density derivative
         dexerfldadrho = vrhocxerflda(i)

         sqs=-z8*sq/(z3*rho)
         fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-z4*f13)/z3
     &    *dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq
     &    +berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)
     &    /(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
         dexerfpbedrho=dexerfldadrho*fx+exerflda*fxs
         sqss=2.6121172985233599567768d-2*rho**(-z8/z3)
         dexerfpbeddrho2=exerflda*berf(1.616204596739954813d-1*mu
     &    *rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa
     &    +berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

C        derivatives
         vrhoc(i) = vrhoc(i) + dexerfpbedrho*0.5d0
         vrhoo(i) = vrhoo(i) + dexerfpbedrho*0.5d0
         vsigmacc(i) = vsigmacc(i) + dexerfpbeddrho2*0.5d0
         vsigmaoo(i) = vsigmaoo(i) + dexerfpbeddrho2*0.5d0
         vsigmaco(i) = vsigmaco(i) + dexerfpbeddrho2

C        End of derivative
         end if

C     End of loop over the grid
2     continue

C Compute LDA part
      do i=1,npt
         zkxerflda(i)=z0
         vrhocxerflda(i)=z0
         vrhooxerflda(i)=z0
      end do

      do i=1,npt
         q(izk-1+i)=max((rhoc(i)-rhoo(i))*.5d0,0d0)*2.D0
      end do

      call dftfun_exerf(namedummy,fderiv,.false.,igrad,npt,q(izk),rhoo,
     >                   zkxerflda,vrhocxerflda,vrhooxerflda,mu)

C     Loop over the grid
      do 3 i=1,npt

!        Beta density
         rho=max((rhoc(i)-rhoo(i))*.5d0,0d0)*2.D0
!        Square of beta density gradient
         drho2=.25d0*sigmacc(i)+.25d0*sigmaoo(i)-.5d0*sigmaco(i)
         drho2=max(drho2,0d0)*4.0d0

!        test on density
         if (rho.lt.tol) goto 3

!        LDA energy density
         exerflda = zkxerflda(i)

         kappa=0.804d0
         sq=drho2*2.6121172985233599567768d-2*rho**(-z8/z3)
         fx=z1+kappa
     &    -kappa/(z1+berf(1.616204596739954813d
     &    -1*mu*rho**(-f13))*sq/kappa)
         exerfpbe=exerflda*fx

         zk(i) = zk(i) + exerfpbe*0.5d0

C Derivative

         if (fderiv) then

!        LDA energy density derivative
         dexerfldadrho = vrhocxerflda(i)

         sqs=-z8*sq/(z3*rho)
         fxs=kappa**2*(-1.616204596739954813d-1*mu*rho**(-z4*f13)/z3
     &    *dberfda(1.616204596739954813d-1*mu*rho**(-f13))*sq
     &    +berf(1.616204596739954813d-1*mu*rho**(-f13))*sqs)
     &    /(kappa+berf(1.616204596739954813d-1*mu*rho**(-f13))*sq)**2
         dexerfpbedrho=dexerfldadrho*fx+exerflda*fxs
         sqss=2.6121172985233599567768d-2*rho**(-z8/z3)
         dexerfpbeddrho2=exerflda*berf(1.616204596739954813d-1*mu
     &    *rho**(-1.d0/3.d0))*sqss*kappa**2/(kappa
     &    +berf(1.616204596739954813d-1*mu*rho**(-1.d0/3.d0))*sq)**2

C        derivatives
         vrhoc(i) = vrhoc(i) + dexerfpbedrho*0.5d0
         vrhoo(i) = vrhoo(i) - dexerfpbedrho*0.5d0
         vsigmacc(i) = vsigmacc(i) + dexerfpbeddrho2*0.5d0
         vsigmaoo(i) = vsigmaoo(i) + dexerfpbeddrho2*0.5d0
         vsigmaco(i) = vsigmaco(i) - dexerfpbeddrho2

C        End of derivative
         end if

C     End of loop over the grid
3     continue

      end if ! end  of (.not.open)

C debug beg
      if(ldebug) then
       print*,trim(here),': exiting'
      end if
C debug end

C First-type gradient functional
      igrad=1

      return
      end

!-------------------------------------------
      function berf(a)
!-------------------------------------------
!  Second-order exchange gradient expansion coefficient for erf
!  interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none

      double precision a
      double precision eta,fak,berf
      include "constants.h"

! function
      double precision erf

      eta=19.0d0
      fak=2.540118935556d0*dexp(-eta*a*a)

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-z7+72.d0*a*a)
     >        /(27.d0*(-z3-24.d0*a*a+32.d0*a**4+z8*dsqrt(pi)*a))

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)
     >        - 23.d0/(358400.d0*a**6)

      else


!      Code generated by Mathematica
      include "berf.inc"

      end if

      berf=berf*fak

      return
      end

!-------------------------------------------
      function dberfda(a)
!-------------------------------------------
!  Derivative of second-order exchange gradient
!  expansion coefficient for erf interaction
!  a = mu/(2*kF)
!
!  Author : J. Toulouse
!  Date   : 10-03-04
!-------------------------------------------
      implicit none

      double precision a
      double precision eta,fak,dfakda,berf,dberfda
      double precision t1,t2,t3,t4,t5
      include "constants.h"

! function
      double precision erf

      eta=19.0d0
      fak=2.540118935556d0*dexp(-eta*a*a)
      dfakda=-2.0d0*eta*a*fak

      if(a .lt. 0.075d0) then
!      expansion for small mu to avoid numerical problems
!      denominator becomes zero for a approximately 0.4845801308
!      (and for one negative and two complex values of a)
       berf = (-z7+72.d0*a*a)
     >        /(27.d0*(-z3-24.d0*a*a+32.d0*a**4+z8*dsqrt(pi)*a))
       dberfda = (z8*(-96.d0*a + 112.d0*a**3 - 576.d0*a**5
     >  + z7*dsqrt(pi) + 72.d0*a**2*dsqrt(pi)))/
     >  (27.d0*(z3 + 24.d0*a**2 - 32.d0*a**4 - z8*a*dsqrt(pi))**2)

      else if(a .gt. 50.d0) then
       berf = 1.d0/(72.d0*a*a)-1.d0/(17280.d0*a**4)
     >        - 23.d0/(358400.d0*a**6)
       dberfda = - 1.d0/(36.d0*a**3) +  1.d0/(4320.d0*a**5)
     >        + 69.d0/(179200.d0*a**7)


      else

!      Code generated by Mathematica
      include "berf.inc"
      include "dberfda.inc"

      end if

      dberfda=dberfda*fak+berf*dfakda

      return
      end

