*fordeck dftfun_ecerf $Revision: 2006.0 $
      subroutine dftfun_ecerf(name,fderiv,open,igrad,npt,
     >    rhoc,rhoo,zk,vrhoc,vrhoo,amu)
C************************************************************************
C
C     Short Range local correlation functional complememtary to erf
C     with the exact mu -> infinity limit imposed
C     Author: Paola Gori-Giorgi
C     Date: 04-01-06
C
C************************************************************************

      implicit double precision (a-h,o-z)

c      include "common/chirs"
      include "debug"

      parameter(tol=1d-12)
      logical fderiv,open
      character*(*) name
      character*(30) xxxx

      dimension zk(*)
      dimension vrhoc(*),rhoc(*)
      dimension vrhoo(*),rhoo(*)

C debug beg
c      here='dftfun_ecerf'
c      ldebug = index(debuglist,here(1:lenstr(here))) .ne. 0
c     .    .or. index(debuglist,'all') .ne. 0
c      if(ldebug) print*,here(1:lenstr(here)),': entering'
C debug end
      ldebug=.false.

C Description of the functional for output of molpro
      name='SR LDA correlation for erf with exact large mu limit'

C Double precision numbers
      f13 = 1.0d0/3.0d0
      pi=dacos(-1.d0)
      rsfac = (3.0d0/(4.0d0*pi))**f13

      igrad=0

      if(.not.open) then

C       Loop over the grid
        do 1 i=1,npt

C         Test on density
          if (dabs(rhoc(i)).lt.tol) goto 1

          if(ldebug) write(*,*)"rho=",rhoc(i)

          rs = rsfac/(rhoc(i)**f13)

c     rs=(3.d0/(4.d0*pi*rhoc(i)))**f13

          call ecPW(rs,0.0d0,ec,ecd,ecz)
          call ecorrlr(rs,0.0d0,amu,eclr)
          zk(i)=(ec-eclr)*rhoc(i)

          if(ldebug) write(*,*)"zk=",zk(i)

C         Derivative
          if (fderiv) then


            vcup=ec-rs/3.d0*ecd+ecz
            vcdown=ec-rs/3.d0*ecd-ecz
            call vcorrlr(rs,0.0d0,amu,vclrup,vclrdown)
            vrhoc(i)=vrhoc(i)+0.5d0*(vcup-vclrup+vcdown-vclrdown)

            if(ldebug) write(*,*)" vrhoc=",
     >      0.5d0*(vcup-vclrup+vcdown-vclrdown)

C         End of derivative
          end if

C       End of loop over the grid
1       continue

      else !else of (.not.open)

C       Loop over the grid
        do 2 i=1,npt

C         Test on density
          if (dabs(rhoc(i)).lt.tol) goto 2

          rs=rsfac/(rhoc(i)**f13)
          rhoa=max((rhoc(i)+rhoo(i))*.5d0,1.0d-15)
          rhob=max((rhoc(i)-rhoo(i))*.5d0,1.0d-15)
          z=(rhoa-rhob)/(rhoa+rhob)

          call ecPW(rs,z,ec,ecd,ecz)
          call ecorrlr(rs,z,amu,eclr)
          zk(i)=(ec-eclr)*rhoc(i)

C         Derivative
          if (fderiv) then

            vcup=ec-rs/3.d0*ecd-(z-1.d0)*ecz
            vcdown=ec-rs/3.d0*ecd-(z+1.d0)*ecz
            call vcorrlr(rs,z,amu,vclrup,vclrdown)
            vrhoc(i)=vrhoc(i)+0.5d0*(vcup-vclrup+vcdown-vclrdown)
            vrhoo(i)=vrhoo(i)+0.5d0*(vcup-vclrup-vcdown+vclrdown)

C         End of derivative
          end if

C       End of loop over the grid
2       continue

      end if !end of (.not.open)

C debug beg
c      if(ldebug) then
c       print*,here(1:lenstr(here)),': exiting'
c      end if
C debug end

      return
      end

      double precision function g0f(x)
ccc on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
ccc x -> rs
      implicit none
      double precision C0f,D0f,E0f,F0f,x
      C0f             = 0.0819306d0
      D0f             = 0.752411d0
      E0f             = -0.0127713d0
      F0f             = 0.00185898d0
      g0f=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+
     $     F0f*x**4)*exp(-abs(D0f)*x)/2.d0
      return
      end


