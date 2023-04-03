C
C   WRITTEN BY PAOLA GORI-GIORGI (DOWNLOADED FROM HER WEB SITE) 
C
C   COMMENT OF KP: ESR=Int Rho*excsr dr
C                  VxcSR=d ESR/d Rho 
C
ccc example: use the subroutine lsdsr to compute the complementary
ccc short-range exchange-correlation energy 'excsr' and
ccc the corresponding up and down potentials 'vxcsrup','vxcsrdown'
ccc at polarization z=0.7, cutoff mu=0.5, and for 0.2 < rs < 20,
ccc and write them on a file
c      program testex
c      implicit none
c      double precision z,rs,mu
c      double precision excsr,vxcsrup,vxcsrdown
c      integer i
c      open(9,file='testex',status='unknown')
c      z=0.7d0
c      mu=0.5d0
c      do i=1,100
c         rs=0.2*i
c         call lsdsr(rs,z,mu,excsr,vxcsrup,vxcsrdown)
c         write(9,*) rs,excsr,vxcsrup,vxcsrdown
c      enddo
c      stop
c      end


      subroutine lsdsr(rs,z,mu,excsr,vxcsrup,vxcsrdown)
ccc Hartree atomic units used
ccc for given density parameter 'rs', spin polarization 'z'
ccc and cutoff parameter 'mu' 
ccc gives the complementary  short-range exchange-correlation
ccc energy  (i.e., xc energy of jellium minus xc energy of long-range
ccc interacting electron gas) => 'excsr'
ccc and the corresponding exchange-correlation potentials for
ccc spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
ccc from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601353
      implicit none
      double precision rs,z,mu,excsr,vxcsrup,vxcsrdown
      double precision eclr,exlr,ec,ecd,ecz,ex
      double precision vclrup,vclrdown,vxlrup,vxlrdown
      double precision vxup,vxdown,vcup,vcdown
      double precision pi,alpha,cf
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      ex=-3.d0*cf/rs/8.d0/pi*((1.d0+z)**(4.d0/3.d0)+
     $     (1.d0-z)**(4.d0/3.d0))

      vxup=-(1.d0+z)**(1.d0/3.d0)*(3.d0/2.d0/pi)**(2.d0/3.d0)/rs
      vxdown=-(1.d0-z)**(1.d0/3.d0)*(3.d0/2.d0/pi)**(2.d0/3.d0)/rs

      call ecPW(rs,z,ec,ecd,ecz)
      vcup=ec-rs/3.d0*ecd-(z-1.d0)*ecz
      vcdown=ec-rs/3.d0*ecd-(z+1.d0)*ecz

      call exchangelr(rs,z,mu,exlr)
      call vexchangelr(rs,z,mu,vxlrup,vxlrdown)

      call ecorrlr(rs,z,mu,eclr)
      call vcorrlr(rs,z,mu,vclrup,vclrdown)

      excsr=ex+ec-(exlr+eclr)
      vxcsrup=vxup+vcup-(vxlrup+vclrup)
      vxcsrdown=vxdown+vcdown-(vxlrdown+vclrdown)

      return
      end



      subroutine ecorrlr(rs,z,mu,eclr)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff parameter mu 
ccc gives the correlation energy of the LR gas
ccc  => eclr
      implicit none
      double precision rs,z,mu,eclr,ec,ecd,ecz
      double precision pi,alpha,cf,phi
      double precision g0,dpol,d2anti,d3anti,Qrpa
      double precision coe2,coe3,coe4,coe5
      double precision a1,a2,a3,a4,b0
      double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
cc parameters from the fit
      adib   = 0.784949d0   
      q1a    = -0.388d0   
      q2a    = 0.676d0   
      q3a    = 0.547d0   
      t1a    = -4.95d0   
      t2a    = 1.d0    
      t3a    = 0.31d0   

      b0=adib*rs

      d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
      d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2) 
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
      endif

      call ecPW(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.d0*ec)
      
      eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

      return
      end

      subroutine vcorrlr(rs,z,mu,vclrup,vclrdown)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the correlation LSD potential for LR interaction
ccc  => vclrup (spin-up electrons), vclrdown (spin-down electrons)
      implicit none
      double precision rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
      double precision ec,ecd,ecz
      double precision pi,alpha,cf,phi
      double precision g0,dpol,d2anti,d3anti,Qrpa
      double precision g0d,dpold,d2antid,d3antid,Qrpad,x
      double precision coe2,coe3,coe4,coe5
      double precision coe2rs,coe3rs,coe4rs,coe5rs
      double precision coe2z,coe3z,coe4z,coe5z
      double precision a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs,
     $     b0rs,a1z,a2z,a3z,a4z,a5z,b0z
      double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
      
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
cc parameters from the fit
      adib   = 0.784949d0   
      q1a    = -0.388d0   
      q2a    = 0.676d0   
      q3a    = 0.547d0   
      t1a    = -4.95d0   
      t2a    = 1.d0    
      t3a    = 0.31d0   

      b0=adib*rs

      d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
      d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2

      d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/
     -    rs**2)*exp(-q3a*rs)
      d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/
     -    rs**3)*exp(-rs*t3a)

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0(rs)-0.5d0)
      coe2rs=-3.d0/8.d0/rs**3*(1.d0-z**2)*g0d(rs)+
     $     9.d0/8.d0/rs**4*(1.d0-z**2)*(g0(rs)-0.5d0)
      coe2z=-3.d0/8.d0/rs**3*(-2.d0*z)*(g0(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**3)
      coe3rs=-(1.d0-z**2)*g0d(rs)/(sqrt(2.d0*pi)*rs**3)+
     $    3.d0*(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**4) 
      coe3z=2.d0*z*g0(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(dpold(rs)
     $        +2.d0*cf**2*2**(5.d0/3.d0)/5.d0/rs**3)
        coe4z=-9.d0/64.d0/rs**3*(dpol(rs)-rs/6.d0*dpold(rs)-2.d0*d2anti
     $       -4.d0/15.d0/rs**2*cf**2*2.d0**(5.d0/3.d0))*z
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)
        coe5rs=-3.d0/rs*coe5-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpold(rs)
        coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(dpol(rs)-rs/6.d0*
     $       dpold(rs)-2.d0*d3anti)*z

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)
         coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.d0-z)**(8.d0/3.d0))/rs**3)
         coe4z=-9.d0/64.d0/rs**3*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d2anti-
     $        4.d0/15.d0/rs**2*cf**2*((1.d0+z)**(5.d0/3.d0)-
     $        (1.d0-z)**(5.d0/3.d0)))

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
         coe5rs=-3.d0/rs*coe5-9.d0/(40.d0*sqrt(2.d0*pi)*rs**3)*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3antid)
         coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d3anti)

      endif

      call ecPW(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a1rs=24.d0*adib*b0**5*coe3+4.d0*b0**6*coe3rs+8.d0*adib*b0**7*
     $     coe5+b0**8*coe5rs
      a1z=4.d0*b0**6*coe3z+b0**8*coe5z

      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a2rs=24.d0*adib*b0**5*coe2+4.d0*b0**6*coe2rs+8.d0*adib*b0**7*
     $     coe4+b0**8*coe4rs+24.d0*adib*b0**3*ec+6.d0*b0**4*ecd
      a2z=4.d0*b0**6*coe2z+b0**8*coe4z+6.d0*b0**4*ecz

      a3=b0**8*coe3
      a3rs=8.d0*adib*b0**7*coe3+b0**8*coe3rs
      a3z=b0**8*coe3z

      a4=b0**6*(b0**2*coe2+4.d0*ec)
      a4rs=8.d0*adib*b0**7*coe2+b0**8*coe2rs+24.d0*adib*b0**5*ec+
     $     4.d0*b0**6*ecd
      a4z=b0**6*(b0**2*coe2z+4.d0*ecz)

      a5=b0**8*ec
      a5rs=8.d0*adib*b0**7*ec+b0**8*ecd
      a5z=b0**8*ecz

      x=mu*sqrt(rs)/phi

      eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+a5*mu**8)/((1.d0+b0**2*mu**2)**4)
      
      eclrrs=-4.d0/(1.d0+b0**2*mu**2)*2.d0*adib*b0*mu**2*eclr+
     $     1.d0/((1.d0+b0**2*mu**2)**4)*(phi**2*mu/(2.d0*sqrt(rs))
     $     *Qrpad(x)+
     $     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)


      if(z.eq.1.d0) then
         vclrup=eclr-rs/3.d0*eclrrs
         vclrdown=0.d0
      elseif(z.eq.-1.d0) then
         vclrup=0.d0
         vclrdown=eclr-rs/3.d0*eclrrs
      else

         eclrz=(phi**2*((1.d0+z)**(-1.d0/3.d0)-(1.d0-z)**(-1.d0/3.d0))
     $        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1.d0+z)**(-1.d0/3.d0)
     $        -(1.d0-z)**(-1.d0/3.d0))/3.d0+
     $        a1z*mu**3+a2z*mu**4+a3z*mu**5+
     $        a4z*mu**6+a5z*mu**8)/((1.d0+b0**2*mu**2)**4)

         vclrup=eclr-rs/3.d0*eclrrs-(z-1.d0)*eclrz
         vclrdown=eclr-rs/3.d0*eclrrs-(z+1.d0)*eclrz
      endif
      return
      end


      double precision function g0(x)
ccc on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
ccc x -> rs
      implicit none
      double precision C0f,D0f,E0f,F0f,x
      C0f             = 0.0819306d0    
      D0f             = 0.752411d0     
      E0f             = -0.0127713d0   
      F0f             = 0.00185898d0   
      g0=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+
     $     F0f*x**4)*exp(-abs(D0f)*x)/2.d0
      return
      end

      double precision function g0d(rs)
ccc derivative of on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
      implicit none
      double precision Bg0,Cg0,Dg0,Eg0,Fg0,rs
      Cg0             = 0.0819306d0    
      Fg0             = 0.752411d0     
      Dg0             = -0.0127713d0   
      Eg0             = 0.00185898d0
      Bg0             =0.7317d0-Fg0
      g0d=(-Bg0+2*Cg0*rs+3*Dg0*rs**2+4*Eg0*rs**3)/2.d0*exp(-Fg0*rs)
     -   - (Fg0*(1 - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/
     -   2.d0*exp(-Fg0*rs)
      return
      end


      double precision function dpol(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0 
      p3p    = 0.4319d0   
      dpol=2.d0**(5.d0/3.d0)/5.d0*cf**2/rs**2*(1.d0+(p3p-0.454555d0)*rs)
     $     /(1.d0+p3p*rs+p2p*rs**2)
      return
      end

      double precision function dpold(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0 
      p3p    = 0.4319d0   
      dpold=2.d0**(5.d0/3.d0)/5.d0*cf**2*
     - (-2. + (0.454555 - 4.*p3p)*rs + 
     -    (-4.*p2p + 
     -       (0.90911 - 2.*p3p)*p3p)*rs**2
     -      + p2p*(1.363665 - 3.*p3p)*
     -     rs**3)/
     -  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      return
      end

      double precision function Qrpa(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0 
      c2              = 3.91744d0 
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpa=Acoul*log((1.d0+a2*x+b2*x**2+c2*x**3)/(1.d0+a2*x+d2*x**2))
      return
      end

      double precision function Qrpad(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0 
      c2              = 3.91744d0 
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpad=Acoul*((x*(b2*(2.d0 + a2*x) + 
     -      c2*x*(3.d0 + 2.d0*a2*x) + 
     -      d2*(-2.d0 - a2*x + c2*x**3)))/
     -  ((1.d0 + a2*x + d2*x**2)*
     -    (1.d0 + a2*x + b2*x**2 + c2*x**3)))
      return
      end

      subroutine exchangelr(rs,z,mu,exlr)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the exchange energy of the LR gas
ccc  => exlr
      implicit none
      double precision rs,z,mu,exlr
      double precision pi,alpha,fx,y,erf
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
c KP 21/09/2010
c for very large y exlr=0. impose it to avoid problems (adding too huge numbers
c may not result in 0 due to the lack of precision) 
c      y=mu*alpha*rs/2.d0/2.d0**(1.d0/3.d0)
c      if(y.gt.1.D5) then
c      exlr=0.d0
c      return
c      endif

      if(abs(z).eq.1.d0) then
         y=mu*alpha*rs/2.d0/2.d0**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=mu*fx
      else
         y=mu*alpha*rs/2.d0/(1.+z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=(1.d0+z)*mu*fx/2.d0
         y=mu*alpha*rs/2.d0/(1.-z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=exlr+(1.d0-z)*mu*fx/2.d0
      endif
      return
      end

      subroutine vexchangelr(rs,z,mu,vxlrup,vxlrdown)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the exchange LSD potential for LR interaction
ccc  => vxlrup (spin-up electrons), vxlrdown (spin-down electrons)
      implicit none
      double precision rs,z,mu,vxlrup,vxlrdown
      double precision pi,alpha,fx,fx1,y,exlr,derrs,derz,erf
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      if(z.eq.1.d0) then
         vxlrup=(rs*alpha*mu**2)/
     -   (2**(1.d0/3.d0)*pi) - (rs*alpha*mu**2)/(2**(1.d0/3.d0)*pi)*
     -     exp(-2**(2.d0/3.d0)/(rs**2*alpha**2*mu**2)) - 
     -  (mu*erf(2**(1.d0/3.d0)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrdown=0.d0
      elseif(z.eq.-1.d0) then
         vxlrdown=(rs*alpha*mu**2)/
     -   (2**(1.d0/3.d0)*pi) - (rs*alpha*mu**2)/(2**(1.d0/3.d0)*pi)*
     -     exp(-2**(2.d0/3.d0)/(rs**2*alpha**2*mu**2)) - 
     -  (mu*erf(2**(1.d0/3.d0)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrup=0.d0
      else       
         y=mu*alpha*rs/2.d0/(1.+z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         fx1=(3.d0*(1 + (-4.d0 + 4.d0*exp(-1.d0/(4.d0*y**2)))*y**2))/pi
         derrs=1.d0/4.d0*(1.d0+z)**(2.d0/3.d0)*mu**2*alpha*fx1
         derz=1.d0/2.d0*mu*fx-1.d0/6.d0*fx1*mu*y
         vxlrup=rs/3.d0*derrs+(z-1.d0)*derz
         vxlrdown=rs/3.d0*derrs+(z+1.d0)*derz
         
         y=mu*alpha*rs/2.d0/(1.-z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         fx1=(3.d0*(1 + (-4.d0 + 4.d0*exp(-1.d0/(4.d0*y**2)))*y**2))/pi
         derrs=1.d0/4.d0*(1.d0-z)**(2.d0/3.d0)*mu**2*alpha*fx1
         derz=-1.d0/2.d0*mu*fx+1.d0/6.d0*fx1*mu*y
         vxlrup=vxlrup+rs/3.d0*derrs+(z-1.d0)*derz
         vxlrdown=vxlrdown+rs/3.d0*derrs+(z+1.d0)*derz
      
         call exchangelr(rs,z,mu,exlr)
         vxlrup=exlr-vxlrup
         vxlrdown=exlr-vxlrdown
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c correlation energy and its derivative w.r.t. rs and z at mu=infinity
c Perdew & Wang PRB 45, 13244 (1992)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ecPW(x,y,ec,ecd,ecz)
c in Hartree; ec=ec(rs,zeta)
c x -> rs; y -> zeta
ccc ecd is d/drs ec
ccc ecz is d/dz ec
      implicit none
      double precision pi,f02,ff,x,y,ec,ecd,ec0,ec0d,ec1,ec1d,
     $     aaa,G,Gd,alfac,alfacd,ecz
      pi=dacos(-1.d0)
      
      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)-
     $     2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-log(2.d0))/pi**2
      call  GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,
     $     1.6382d0,0.49294d0,G,Gd)
      ec0=G
      ec0d=Gd

      aaa=aaa/2.d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,
     $     3.3662d0,0.62517d0,G,Gd)
      ec1=G
      ec1d=Gd
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,
     $     0.88026d0,0.49671d0,G,Gd)
      alfac=-G
      alfacd=-Gd

      ec=ec0+alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1.d0-y**4)+(ec1d-ec0d)*
     $     ff*y**4
      ecz=alfac*(-4.d0*y**3)*ff/f02+alfac*(1.d0-y**4)/f02*
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)+(ec1-ec0)*(4.d0*y**3*ff+
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)*y**4)

      return
      end

      subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
ccc Gd is d/drs G
      implicit none
      double precision G,Gd,Ac,alfa1,beta1,beta2,beta3,beta4,x
      G=-2.d0*Ac*(1.d0+alfa1*x)*dlog(1.d0+1.d0/(2.d0*
     $     Ac*(beta1*x**0.5d0+
     $     beta2*x+beta3*x**1.5d0+beta4*x**2)))
      Gd=(1.d0+alfa1*x)*(beta2+beta1/(2.d0*sqrt(x))+3.d0*beta3*
     $     sqrt(x)/2.d0+2.d0*beta4*x)/((beta1*sqrt(x)+beta2*x+
     $     beta3*x**(3.d0/2.d0)+beta4*x**2)**2*(1.d0+1.d0/
     $     (2.d0*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2))))-2.d0*Ac*alfa1*dlog(1.d0+1.d0/(2.d0*Ac*
     $     (beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2)))
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*Deck Erf
      double precision function erf(x)
C
C error function in double precision
C
      implicit real*8 (a - h, o - z)
      dimension a(0 : 64), b(0 : 64)
      data (a(i), i = 0, 12) / 
     &    0.00000000005958930743d0, -0.00000000113739022964d0, 
     &    0.00000001466005199839d0, -0.00000016350354461960d0, 
     &    0.00000164610044809620d0, -0.00001492559551950604d0, 
     &    0.00012055331122299265d0, -0.00085483269811296660d0, 
     &    0.00522397762482322257d0, -0.02686617064507733420d0, 
     &    0.11283791670954881569d0, -0.37612638903183748117d0, 
     &    1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) / 
     &    0.00000000002372510631d0, -0.00000000045493253732d0, 
     &    0.00000000590362766598d0, -0.00000006642090827576d0, 
     &    0.00000067595634268133d0, -0.00000621188515924000d0, 
     &    0.00005103883009709690d0, -0.00037015410692956173d0, 
     &    0.00233307631218880978d0, -0.01254988477182192210d0, 
     &    0.05657061146827041994d0, -0.21379664776456006580d0, 
     &    0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) / 
     &    0.00000000000949905026d0, -0.00000000018310229805d0, 
     &    0.00000000239463074000d0, -0.00000002721444369609d0, 
     &    0.00000028045522331686d0, -0.00000261830022482897d0, 
     &    0.00002195455056768781d0, -0.00016358986921372656d0, 
     &    0.00107052153564110318d0, -0.00608284718113590151d0, 
     &    0.02986978465246258244d0, -0.13055593046562267625d0, 
     &    0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) / 
     &    0.00000000000382722073d0, -0.00000000007421598602d0, 
     &    0.00000000097930574080d0, -0.00000001126008898854d0, 
     &    0.00000011775134830784d0, -0.00000111992758382650d0, 
     &    0.00000962023443095201d0, -0.00007404402135070773d0, 
     &    0.00050689993654144881d0, -0.00307553051439272889d0, 
     &    0.01668977892553165586d0, -0.08548534594781312114d0, 
     &    0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) / 
     &    0.00000000000155296588d0, -0.00000000003032205868d0, 
     &    0.00000000040424830707d0, -0.00000000471135111493d0, 
     &    0.00000005011915876293d0, -0.00000048722516178974d0, 
     &    0.00000430683284629395d0, -0.00003445026145385764d0, 
     &    0.00024879276133931664d0, -0.00162940941748079288d0, 
     &    0.00988786373932350462d0, -0.05962426839442303805d0, 
     &    0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) / 
     &    -0.00000000029734388465d0, 0.00000000269776334046d0, 
     &    -0.00000000640788827665d0, -0.00000001667820132100d0, 
     &    -0.00000021854388148686d0, 0.00000266246030457984d0, 
     &    0.00001612722157047886d0, -0.00025616361025506629d0, 
     &    0.00015380842432375365d0, 0.00815533022524927908d0, 
     &    -0.01402283663896319337d0, -0.19746892495383021487d0, 
     &    0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) / 
     &    -0.00000000001951073787d0, -0.00000000032302692214d0, 
     &    0.00000000522461866919d0, 0.00000000342940918551d0, 
     &    -0.00000035772874310272d0, 0.00000019999935792654d0, 
     &    0.00002687044575042908d0, -0.00011843240273775776d0, 
     &    -0.00080991728956032271d0, 0.00661062970502241174d0, 
     &    0.00909530922354827295d0, -0.20160072778491013140d0, 
     &    0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) / 
     &    0.00000000003147682272d0, -0.00000000048465972408d0, 
     &    0.00000000063675740242d0, 0.00000003377623323271d0, 
     &    -0.00000015451139637086d0, -0.00000203340624738438d0, 
     &    0.00001947204525295057d0, 0.00002854147231653228d0, 
     &    -0.00101565063152200272d0, 0.00271187003520095655d0, 
     &    0.02328095035422810727d0, -0.16725021123116877197d0, 
     &    0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) / 
     &    0.00000000002319363370d0, -0.00000000006303206648d0, 
     &    -0.00000000264888267434d0, 0.00000002050708040581d0, 
     &    0.00000011371857327578d0, -0.00000211211337219663d0, 
     &    0.00000368797328322935d0, 0.00009823686253424796d0, 
     &    -0.00065860243990455368d0, -0.00075285814895230877d0, 
     &    0.02585434424202960464d0, -0.11637092784486193258d0, 
     &    0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) / 
     &    -0.00000000000367789363d0, 0.00000000020876046746d0, 
     &    -0.00000000193319027226d0, -0.00000000435953392472d0, 
     &    0.00000018006992266137d0, -0.00000078441223763969d0, 
     &    -0.00000675407647949153d0, 0.00008428418334440096d0, 
     &    -0.00017604388937031815d0, -0.00239729611435071610d0, 
     &    0.02064129023876022970d0, -0.06905562880005864105d0, 
     &    0.09084526782065478489d0 / 
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + 
     &        a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + 
     &        a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + 
     &        a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + 
     &        a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + 
     &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + 
     &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + 
     &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + 
     &        b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      erf = y
      end function

