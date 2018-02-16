*fordeck dftfun_exerf $Revision: 2006.0 $
C*************************************************************************
      subroutine dftfun_exerf(name,fderiv,open,igrad,npt,rhoc,
     >            rhoo,zk,vrhoc,vrhoo,amu)
C*************************************************************************
C
C     Subroutine dftfun_exerf for short-range lda exchange functional
C     complementary to erf iteraction
C
C     Originally written by T. Leininger
C     Rewritten by Julien Toulouse on 09-05-03
C
C************************************************************************

      implicit double precision (a-h,o-z)

c      include "common/chirs"
      include "debug"

      logical fderiv,open
      character*(*) name
      dimension zk(*)
      dimension vrhoc(*),rhoc(*)
      dimension vrhoo(*),rhoo(*)

C debug beg
c      here='dftfun_exerf'
c      ldebug = index(debuglist,here(1:lenstr(here))) .ne. 0
c     .    .or. index(debuglist,'all') .ne. 0
c      if(ldebug) print*,here(1:lenstr(here)),': entering'
C debug end
      ldebug=.false.

C Description of the functional for output of molpro
      name='Short-range LDA exchange for erf'

C Double precision numbers
      z0  = 0.D0
      z1  = 1.D0
      z2  = 2.D0
      z3  = 3.D0
      z4  = 4.D0
      z6  = 6.D0
      z8  = 8.D0
      z16 = 16.D0
      z24 = 24.D0
      z96 = 96.D0
      f12 = z1/z2
      f13 = z1/z3
      f14 = z1/z4
      f32 = z3/z2

      igrad=0

C Constants
      pi = dacos(-1.d0)
      pisqrt = dsqrt(pi)
      ckf = (z3*pi*pi)**f13

      if (.not.open) then

C       Loop over the grid
        do i=1,npt

C         Density and kF
          rho = rhoc(i)
          akf = ckf*(rho**f13)
          if(ldebug) write(*,*)"rho=",rho
          a = amu/(z2*akf)
          a2 = a*a
          a3 = a2*a

C         Test on the value of a

C         Limit for small a (expansion not so important as for large a)
          if (a.lt.1.d-9) then
            zk(i) = -z3/z8*rho*(z24*rho/pi)**f13
            if (fderiv) then
              vrhoc(i) = vrhoc(i) - ((z3/pi)*rho)**f13
            end if

C         Intermediate values of a
          elseif (a.le.100) then
            zk(i) = - (rho*(z24*rho/pi)**f13)
     >         *(z3/z8-a*(pisqrt*erf(f12/a)+
     >         (z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
            if (fderiv) then
              vrhoc(i) = vrhoc(i) - (z3*rho/pi)**f13
     >          + z2*a*amu/pi*(dexp(-f14/a2)-z1)+amu/pisqrt
     >          * erf(f12/a)
            end if

C         Expansion for large a
          elseif (a.lt.1.d+9) then
            zk(i) = - (rho*(z24*rho/pi)**f13)
     >                *z1/(z96*a2)
            if (fderiv) then
              vrhoc(i) = vrhoc(i)-pi*rho/(z2*amu*amu)
            end if

C         Limit for large a
          else
            zk(i) = z0
          end if

          if(ldebug) write(*,*)"zk=",zk(i)
          if(ldebug) write(*,*)i,"vrhoc=",vrhoc(i)

C       End of loop over the grid
        end do

      else ! else of (.not.open)

C       Loop over the grid
        do i=1,npt

C         Density and kF
          rhoa=max((rhoc(i)+rhoo(i))*.5d0,0d0)*2.D0
          akf = ckf*(rhoa**f13)
          a = amu/(z2*akf)
          a2 = a*a
          a3 = a2*a

C         Test on the value of a

C         Limit for small a (expansion not so important as for large a)
          if (a.lt.1.d-9) then
            zka = -z3/z8*rhoa*(z24*rhoa/pi)**f13
            if (fderiv) then
              vrhoc(i) = vrhoc(i) - ((z3/pi)*rhoa)**f13*f12
              vrhoo(i) = vrhoo(i) - ((z3/pi)*rhoa)**f13*f12
            end if

C         Intermediate values of a
          elseif (a.le.100) then
            zka = - (rhoa*(z24*rhoa/pi)**f13)
     >         *(z3/z8-a*(pisqrt*erf(f12/a)+
     >         (z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
            if (fderiv) then
              vrhoc(i) = vrhoc(i) + (-(z3*rhoa/pi)**f13
     >          + z2*a*amu/pi*(dexp(-f14/a2)-z1)+amu/pisqrt
     >          * erf(f12/a))*f12
              vrhoo(i) = vrhoo(i) + (-(z3*rhoa/pi)**f13
     >          + z2*a*amu/pi*(dexp(-f14/a2)-z1)+amu/pisqrt
     >          * erf(f12/a))*f12
            end if

C         Expansion for large a
          elseif (a.lt.1.d+9) then
            zka = - (rhoa*(z24*rhoa/pi)**f13)
     >                *z1/(z96*a2)
            if (fderiv) then
              vrhoc(i) = vrhoc(i)-pi*rhoa/(z2*amu*amu)*f12
              vrhoo(i) = vrhoo(i)-pi*rhoa/(z2*amu*amu)*f12
            end if

C         Limit for large a
          else
            zka = z0
          end if

C         Density and kF
          rhob=max((rhoc(i)-rhoo(i))*.5d0,0d0)*2.D0
          akf = ckf*(rhob**f13)
          a = amu/(z2*akf)
          a2 = a*a
          a3 = a2*a

C         Test on the value of a

C         Limit for small a (expansion not so important as for large a)
          if (a.lt.1.d-9) then
            zkb = -z3/z8*rhob*(z24*rhob/pi)**f13
            if (fderiv) then
              vrhoc(i) = vrhoc(i) - ((z3/pi)*rhob)**f13*f12
              vrhoo(i) = vrhoo(i) + ((z3/pi)*rhob)**f13*f12
            end if

C         Intermediate values of a
          elseif (a.le.100) then
            zkb = - (rhob*(z24*rhob/pi)**f13)
     >         *(z3/z8-a*(pisqrt*erf(f12/a)+
     >         (z2*a-z4*a3)*dexp(-f14/a2)-z3*a+z4*a3))
            if (fderiv) then
              vrhoc(i) = vrhoc(i) + (-(z3*rhob/pi)**f13
     >          + z2*a*amu/pi*(dexp(-f14/a2)-z1)+amu/pisqrt
     >          * erf(f12/a))*f12
              vrhoo(i) = vrhoo(i) - (-(z3*rhob/pi)**f13
     >          + z2*a*amu/pi*(dexp(-f14/a2)-z1)+amu/pisqrt
     >          * erf(f12/a))*f12
            end if

C         Expansion for large a
          elseif (a.lt.1.d+9) then
            zkb = - (rhob*(z24*rhob/pi)**f13)
     >                *z1/(z96*a2)
            if (fderiv) then
              vrhoc(i) = vrhoc(i)-pi*rhob/(z2*amu*amu)*f12
              vrhoo(i) = vrhoo(i)+pi*rhob/(z2*amu*amu)*f12
            end if

C         Limit for large a
          else
            zkb = z0
          end if

          zk(i) = (zka+zkb) / 2.D0

C       End of loop over the grid
        end do

      end if ! end if of (.not.open)

C debug beg
c      if(ldebug) then
c       print*,here(1:lenstr(here)),': exiting'
c      end if
C debug end

      return
      end
