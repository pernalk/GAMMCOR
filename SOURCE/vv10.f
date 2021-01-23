      Subroutine VV10(URe,UNOAO,Occ,NBasis)
C
C     WRITTEN BY REZA JAN/2021
C
C     COMPUTATION OF THE VV10 NONLOCAL CORRELATION
C
      use types
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Real*8, Dimension(:), Allocatable :: OrbXGrid
      Real*8, Dimension(:), Allocatable :: OrbYGrid
      Real*8, Dimension(:), Allocatable :: OrbZGrid
      Real*8, Dimension(:), Allocatable :: WGrid
      Real*8, Dimension(:), Allocatable :: RhoGrid
      Real*8, Dimension(:), Allocatable :: Sigma
      Real*8, Allocatable :: OrbGrid(:,:)
      Real*8, Allocatable :: RR(:,:)
C
      Real*8, Dimension(:), Allocatable :: omega_p_2
      Real*8, Dimension(:), Allocatable :: omega_g_2
      Real*8, Dimension(:), Allocatable :: kappa
      Real*8, Dimension(:), Allocatable :: gfunc
      Real*8, Dimension(:), Allocatable :: omega_0
      Real*8, Dimension(:), Allocatable :: integ_energy
      Real*8, Dimension(:), Allocatable :: ffunc
C
      Dimension URe(NBasis,NBasis),UNOAO(NBasis,NBasis),Occ(NBasis)
      integer igrad, I, K, L, M
      CHARACTER(100) :: num1char,num2char
C
      Call molprogrid0(NGrid,NBasis)
C
      Allocate  (WGrid(NGrid))
      Allocate  (OrbGrid(NGrid,NBasis))
      Allocate  (OrbXGrid(NBasis*NGrid))
      Allocate  (OrbYGrid(NBasis*NGrid))
      Allocate  (OrbZGrid(NBasis*NGrid))
      Allocate  (RhoGrid(NGrid))
      Allocate  (Sigma(NGrid))
      Allocate  (RR(3,NGrid))
C     add new vectors
      Allocate  (omega_p_2(NGrid))
      Allocate  (omega_g_2(NGrid))
      Allocate  (kappa(NGrid))
      Allocate  (gfunc(NGrid))
      Allocate  (omega_0(NGrid))
      Allocate  (integ_energy(NGrid))
      Allocate  (ffunc(NGrid))
C
C     initializing vectors
      omega_p_2 = 0.0
      omega_g_2 = 0.0
      kappa = 0.0
      gfunc = 0.0
      omega_0 = 0.0
      integ_energy = 0.0
      Sigma = 0.0
      ffunc = 0.0
C
      IF(COMMAND_ARGUMENT_COUNT().Eq.0) THEN
      Bo = 5.9
      Co = 0.0093
      ELSE
      CALL GET_COMMAND_ARGUMENT(1,num1char)   !B first
      CALL GET_COMMAND_ARGUMENT(2,num2char)   !C second
      READ(num1char,*)xnum1          
      READ(num2char,*)xnum2  
      Bo = xnum1
      Co = xnum2
      Write(6,'(/,X,"Parameters B and C from a command line: ",
     $ 2F10.6)')Bo,Co
      ENDIF

      PI = 4.0*ATAN(1.0)
C
      thresh = 1.e-8
      thresh_g = 1.e-8
C
      Call molprogrid1(RR,NGrid)
      Call molprogrid(OrbGrid,OrbXGrid,OrbYGrid,OrbZGrid,
     $                                         WGrid,UNOAO,NGrid,NBasis)
C
      Do I=1,NGrid
            Call DenGrid(I,RhoGrid(I),Occ,URe,OrbGrid,NGrid,NBasis)
            Call DenGrad(I,RhoX,Occ,URe,OrbGrid,OrbXGrid,NGrid,NBasis)
            Call DenGrad(I,RhoY,Occ,URe,OrbGrid,OrbYGrid,NGrid,NBasis)
            Call DenGrad(I,RhoZ,Occ,URe,OrbGrid,OrbZGrid,NGrid,NBasis)
            Sigma(I)=RhoX**2+RhoY**2+RhoZ**2
      EndDo
C     beta is density-independent constant from the uniform electron gas
C     equation 11
C
      beta = (0.03125)*((3.0/(Bo*Bo))**(0.75))
C
C     The local plasma frequency omega_p_2 in equation 5
C
      Do I=1,NGrid
            omega_p_2(I) = 4.0 * PI * RhoGrid(I)
      EndDo
C
C     The local band gap omega_g_2  equation 6
C
      Do I=1, NGrid
            If(RhoGrid(I).Gt.thresh)
     $              omega_g_2(I) = ( Co * (Sigma(I)**2))/(RhoGrid(I)**4)
C
C     omega_0  equation 5
            omega_0(I) = SQRT(omega_g_2(I) + (omega_p_2(I)/3.0))
C
C     Kappa formula equation 14
            kappa(I) = (Bo * PI * 1.5)*((RhoGrid(I)/(9.0*PI))**
     $                                                        (1.0/6.0))
      EndDo
C
C     The nonlocal correlation energy equation 1
      Do K=1, NGrid
            Do L=1, NGrid
                  dist_2 = 0.0
                  Do M=1,3
                        dist_2 = dist_2 + ((RR(M,K)-RR(M,L))**2.0)
                  EndDo
C
                  func_g = 0.0
                  func_g_prim = 0.0
                  func_phi = 0.0
C
                  func_g = omega_0(K) * dist_2 + kappa(K)
                  func_g_prim = omega_0(L) * dist_2 + kappa(L)
                  If((func_g_prim.gt.thresh_g).and.(func_g.gt.thresh_g))
     $            func_phi = (-1.5) / ((func_g_prim + func_g) *
     $                                           (func_g_prim * func_g))
C
C Q to Reza: should it be "or" or "and" below?
                  If((RhoGrid(K).Gt.thresh).or.(RhoGrid(L).Gt.thresh))
     $            ffunc(K) = ffunc(K) +
     $            (0.5 * func_phi * RhoGrid(L) * WGrid(L) * RhoGrid(K))
            EndDo
            ffunc(K) = ffunc(K) + (beta * RhoGrid(K))
      EndDo
C
      VV10_Energy = 0.0
      Do I=1, NGrid
            VV10_Energy = VV10_Energy + (ffunc(I) * WGrid(I))
      EndDo
C
      Write(6,'(/,X,"***************** VV10 NONLOCAL FUNCTIONAL
     $ ******************* ")')
      write (*,*)
      write (*,*) 'NonLocal Correlation (VV10) Energy  = ' , VV10_Energy

      Return
      End
