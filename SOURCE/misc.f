*Deck ReWr
      Subroutine ReWr(IFlag,Occ,URe,Title,NBasis)
C
C     READ (IFlag=0) OR WRITE (IFlag=1) FROM/TO A RESTART FILE
C
      Implicit Real*8 (A-H,O-Z)
C
      Character*60 Title,FName
C
      Include 'commons.inc'
C
      Dimension Occ(NBasis),URe(NBasis*NBasis)
C
      Do I=1,60
      FName(I:I)=' '
      EndDo
C
      K=0
    5 K=K+1
      If (Title(K:K).Ne.' ') Then
      FName(K:K)=Title(K:K)
      GoTo 5
      EndIf
      FName(K:K+4)='.res'
C
      If(IFlag.Eq.0) Then
C
      Open(10,Form='Unformatted',File=FName,Status='Old')
      Read(10) NBasis
      Read(10) IFun
      Read(10) (Occ(I),I=1,NBasis)
      Read(10) (URe(I),I=1,NBasis*NBasis)
      If(IFun.Eq.6) Read(10) (IType(I),I=1,NBasis)
C
      Close(10)
C
      Else
C
      Open(10,Form='Unformatted',File=FName)
      Write(10) NBasis
      Write(10) IFun
      Write(10) (Occ(I),I=1,NBasis)
      Write(10) (URe(I),I=1,NBasis*NBasis)
      If(IFun.Eq.6) Write(10) (IType(I),I=1,NBasis)
      Close(10)
C
      EndIf
C
      Return
      End

*Deck TrTwoEl
      Subroutine TrTwoEl(CoulNO,ExchNO,URe,TwoEl,NBasis,NInte2)
C
C     CALCULATE COULOMB AND EXCHANGE INTEGRALS IN NO REPRESENTATION
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0)
C
      Parameter(Tol=1.D-12)
C
      Dimension
     $ CoulNO(NBasis*(NBasis+1)/2),ExchNO(NBasis*(NBasis+1)/2),
     $ URe(NBasis,NBasis),TwoEl(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension Hlp1(NBasis,NBasis*(NBasis+1)/2),
     $ Hlp2(NBasis,NBasis*(NBasis+1)/2)
C
      NInte1=NBasis*(NBasis+1)/2
C
      Do 10 I=1,NInte1
      CoulNO(I)=Zero
      ExchNO(I)=Zero
      Do 10 J=1,NBasis
      Hlp1(J,I)=Zero
  10  Hlp2(J,I)=Zero
C
      NAddr=0
C
      IJ=0
      Do 20 I=1,NBasis
      Do 20 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 30 K=1,NBasis
      Do 30 L=1,K
      KL=KL+1
C
      If(IJ.Lt.KL) GoTo 30
      NAddr=NAddr+1
      TwoZet=Two*TwoEl(NAddr)
C
      If(Abs(TwoZet).Gt.Tol) Then
C
      FacIJ=One
      If((I.Eq.J).And.(K.Ne.L)) FacIJ=Two
      FacKL=One
      If((K.Eq.L).And.(I.Ne.J)) FacKL=Two
      FacIL=One
      If((I.Eq.L).And.(K.Ne.J)) FacIL=Two
      FacJK=One
      If((J.Eq.K).And.(I.Ne.L)) FacJK=Two
      FacIK=One
      If((I.Eq.K).And.(L.Ne.J)) FacIK=Two
      FacJL=One
      If((J.Eq.L).And.(I.Ne.K)) FacJL=Two
C
      IL=(Max(I,L)*(Max(I,L)-1))/2+Min(I,L)
      JK=(Max(J,K)*(Max(J,K)-1))/2+Min(J,K)
      IK=(Max(I,K)*(Max(I,K)-1))/2+Min(I,K)
      JL=(Max(J,L)*(Max(J,L)-1))/2+Min(J,L)
C
      Do 40 IA=1,NBasis
C
      Hlp1(IA,IL)=Hlp1(IA,IL)+Half*FacIL*URe(IA,J)*URe(IA,K)*TwoZet
      If(IL.Ne.JK) Hlp1(IA,JK)=Hlp1(IA,JK)
     $ +Half*FacJK*URe(IA,I)*URe(IA,L)*TwoZet
C
      Hlp2(IA,IJ)=Hlp2(IA,IJ)+Half*FacIJ*URe(IA,K)*URe(IA,L)*TwoZet
      If(IJ.Ne.KL) Hlp2(IA,KL)=Hlp2(IA,KL)
     $ +Half*FacKL*URe(IA,I)*URe(IA,J)*TwoZet
C
   40 Continue
C
      If(K.Eq.L.Or.I.Eq.J) GoTo 30
C
      Do 50 IA=1,NBasis
C
      Hlp1(IA,IK)=Hlp1(IA,IK)+Half*FacIK*URe(IA,J)*URe(IA,L)*TwoZet
      If(IK.Ne.JL) Hlp1(IA,JL)=Hlp1(IA,JL)
     $ +Half*FacJL*URe(IA,I)*URe(IA,K)*TwoZet
C
      Hlp2(IA,IJ)=Hlp2(IA,IJ)+Half*URe(IA,K)*URe(IA,L)*TwoZet
      If(IJ.Ne.KL) Hlp2(IA,KL)=Hlp2(IA,KL)
     $ +Half*URe(IA,I)*URe(IA,J)*TwoZet
C
   50 Continue
C
      EndIf
C
   30 Continue
   20 Continue
C
      IJ=0
      Do 60 I=1,NBasis
      Do 60 J=1,I
      IJ=IJ+1
C
      IAB=0
      Do 60 IA=1,NBasis
      Do 60 IB=1,IA
      IAB=IAB+1
      FacAB=Two
      If (IA.Eq.IB) FacAB=One
      CoulNO(IJ)=CoulNO(IJ)+FacAB*URe(J,IA)*URe(J,IB)*Hlp2(I,IAB)
   60 ExchNO(IJ)=ExchNO(IJ)+FacAB*URe(J,IA)*URe(J,IB)*Hlp1(I,IAB)
C
      Return
      End

*Deck TwoNO1
      Subroutine TwoNO1(TNO,URe,NBasis,NInte2)
C
C     FULL TRANSFORMATION OF THE TWO-ELECTRON INTEGRALS TO NO
C     INTEGRALS STORED IN TNO ARE TRANSFORMED AND RETURNED IN TNO
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Parameter (Tol2=1.D-12)
C
      Dimension TNO(NInte2),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: X(:,:,:)
      Real*8, Allocatable :: Y(:,:)
C
      Allocate (X(NBasis,NBasis,NBasis*(NBasis+1)/2))
      Allocate (Y(NBasis*(NBasis+1)/2,NBasis*(NBasis+1)/2))
C
      NInte1=NBasis*(NBasis+1)/2
C
C     FIRST INDEX
C
      Do 10 I=1,NInte1
      Do 10 J=1,NBasis
      Do 10 K=1,NBasis
   10 X(K,J,I)=Zero
C
      NAddr=0
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 40 K=1,NBasis
      Do 40 L=1,K
      KL=KL+1
C
      If(IJ.Lt.KL) GoTo 40
      NAddr=NAddr+1
      TwoZet=TNO(NAddr)
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      Do 50 IA=1,NBasis
      X(IA,J,KL)=X(IA,J,KL)+URe(IA,I)*TwoZet
      If(I.Ne.J) X(IA,I,KL)=X(IA,I,KL)+URe(IA,J)*TwoZet
   50 Continue
C
      If(IJ.Eq.KL) GoTo 40
C
      Do 60 IA=1,NBasis
      X(IA,L,IJ)=X(IA,L,IJ)+URe(IA,K)*TwoZet
      If(K.Ne.L) X(IA,K,IJ)=X(IA,K,IJ)+URe(IA,L)*TwoZet
   60 Continue
C
      EndIf
C
   40 Continue
   30 Continue
C
C     SECOND INDEX
C
      Do 70 I=1,NInte1
      Do 70 J=1,NInte1
   70 Y(J,I)=Zero
C
      IJ=0
      Do 80 I=1,NBasis
      Do 80 J=1,I
      IJ=IJ+1
      Do 80 K=1,NBasis
C
      IAB=0
      Do 80 IA=1,NBasis
      XAKIJ=X(IA,K,IJ)
C
      If(Abs(XAKIJ).Gt.Tol2) Then
      IAA=IA*(IA-1)/2
      Do 95 IB=1,IA
      IAB=IAA+IB
   95 Y(IAB,IJ)=Y(IAB,IJ)+URe(IB,K)*XAKIJ
      EndIf
C
   80 Continue
C
C     THIRD INDEX
C
      Do 100 I=1,NInte1
      Do 100 J=1,NBasis
      Do 100 K=1,NBasis
  100 X(K,J,I)=Zero
C
      Do 105 I=1,NInte1
      Do 105 J=1,I-1
      Hlp=Y(J,I)
      Y(J,I)=Y(I,J)
  105 Y(I,J)=Hlp

      IAB=0
      Do 110 IA=1,NBasis
      Do 110 IB=1,IA
      IAB=IAB+1
C
      IJ=0
      Do 110 I=1,NBasis
      Do 110 J=1,I
      IJ=IJ+1
C
      YABIJ=Y(IJ,IAB)

      If(Abs(YABIJ).Gt.Tol2) Then
      Do 120 IC=1,NBasis
      X(IC,J,IAB)=X(IC,J,IAB)+URe(IC,I)*YABIJ
      If(I.Ne.J) X(IC,I,IAB)=X(IC,I,IAB)+URe(IC,J)*YABIJ
  120 Continue
      EndIf
C
  110 Continue
C
C     FOURTH INDEX
C
      Do 125 I=1,NInte2
  125 TNO(I)=Zero
C
      NAddr=0
      IAB=0
      Do 130 IA=1,NBasis
      Do 130 IB=1,IA
      IAB=IAB+1
      ICD=0
      Do 130 IC=1,NBasis
      Do 130 ID=1,IC
      ICD=ICD+1
C
      If(IAB.Ge.ICD) Then
      NAddr=NAddr+1
      Do 140 J=1,NBasis
  140 TNO(NAddr)=TNO(NAddr)+URe(ID,J)*X(IC,J,IAB)
      EndIf
C
  130 Continue
C
      Deallocate(X)
      Deallocate(Y)
C
      Return
      End

*Deck TwoNO
      Subroutine TwoNO(TNO,URe,TwoEl,NBasis,NInte2)
C
C     FULL TRANSFORMATION OF THE TWO-ELECTRON INTEGRALS TO NO
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      Parameter (Tol2=1.D-12)
C
      Dimension TNO(NInte2),TwoEl(NInte2),URe(NBasis,NBasis)
C
C     LOCAL ARRAYS
C
      Real*8, Allocatable :: X(:,:,:)
      Real*8, Allocatable :: Y(:,:)
C
      Allocate (X(NBasis,NBasis,NBasis*(NBasis+1)/2))
      Allocate (Y(NBasis*(NBasis+1)/2,NBasis*(NBasis+1)/2))
C
      NInte1=NBasis*(NBasis+1)/2
C
C     FIRST INDEX
C
      Do 10 I=1,NInte1
      Do 10 J=1,NBasis
      Do 10 K=1,NBasis
   10 X(K,J,I)=Zero
C
      NAddr=0
C
      IJ=0
      Do 30 I=1,NBasis
      Do 30 J=1,I
      IJ=IJ+1
C
      KL=0
      Do 40 K=1,NBasis
      Do 40 L=1,K
      KL=KL+1
C
      If(IJ.Lt.KL) GoTo 40
      NAddr=NAddr+1
      TwoZet=TwoEl(NAddr)
C
      If(Abs(TwoZet).Gt.Tol2) Then
C
      Do 50 IA=1,NBasis
      X(IA,J,KL)=X(IA,J,KL)+URe(IA,I)*TwoZet
      If(I.Ne.J) X(IA,I,KL)=X(IA,I,KL)+URe(IA,J)*TwoZet
   50 Continue
C
      If(IJ.Eq.KL) GoTo 40
C
      Do 60 IA=1,NBasis
      X(IA,L,IJ)=X(IA,L,IJ)+URe(IA,K)*TwoZet
      If(K.Ne.L) X(IA,K,IJ)=X(IA,K,IJ)+URe(IA,L)*TwoZet
   60 Continue
C
      EndIf
C
   40 Continue
   30 Continue
C
C     SECOND INDEX
C
      Do 70 I=1,NInte1
      Do 70 J=1,NInte1
   70 Y(J,I)=Zero
C
      IJ=0
      Do 80 I=1,NBasis
      Do 80 J=1,I
      IJ=IJ+1
      Do 80 K=1,NBasis
C
      IAB=0
      Do 80 IA=1,NBasis
      XAKIJ=X(IA,K,IJ)
C
      If(Abs(XAKIJ).Gt.Tol2) Then
      IAA=IA*(IA-1)/2
      Do 95 IB=1,IA
      IAB=IAA+IB
   95 Y(IAB,IJ)=Y(IAB,IJ)+URe(IB,K)*XAKIJ
      EndIf
C
   80 Continue
C
C     THIRD INDEX
C
      Do 100 I=1,NInte1
      Do 100 J=1,NBasis
      Do 100 K=1,NBasis
  100 X(K,J,I)=Zero
C
      Do 105 I=1,NInte1
      Do 105 J=1,I-1
      Hlp=Y(J,I)
      Y(J,I)=Y(I,J)
  105 Y(I,J)=Hlp

      IAB=0
      Do 110 IA=1,NBasis
      Do 110 IB=1,IA
      IAB=IAB+1
C
      IJ=0
      Do 110 I=1,NBasis
      Do 110 J=1,I
      IJ=IJ+1
C
      YABIJ=Y(IJ,IAB)

      If(Abs(YABIJ).Gt.Tol2) Then
      Do 120 IC=1,NBasis
      X(IC,J,IAB)=X(IC,J,IAB)+URe(IC,I)*YABIJ
      If(I.Ne.J) X(IC,I,IAB)=X(IC,I,IAB)+URe(IC,J)*YABIJ
  120 Continue
      EndIf
C
  110 Continue
C
C     FOURTH INDEX
C
      Do 125 I=1,NInte2
  125 TNO(I)=Zero
C
      NAddr=0
      IAB=0
      Do 130 IA=1,NBasis
      Do 130 IB=1,IA
      IAB=IAB+1
      ICD=0
      Do 130 IC=1,NBasis
      Do 130 ID=1,IC
      ICD=ICD+1
C
      If(IAB.Ge.ICD) Then
      NAddr=NAddr+1
      Do 140 J=1,NBasis
  140 TNO(NAddr)=TNO(NAddr)+URe(ID,J)*X(IC,J,IAB)
      EndIf
C
  130 Continue
C
      Deallocate(X)
      Deallocate(Y)
C
      Return
      End

C*Deck NAddr3
C      Integer(8) Function NAddr3(IAddr1,IAddr2,IAddr3,IAddr4)
CC
CC     POINTER FOR TWO-ELECTRON INTEGRALS
CC     NAddr3 WILL POINT TO THE INTEGRAL (IAddr1,IAddr2,IAddr3,IAddr4)
CC
C      Implicit Real*8 (A-H,O-Z)
CC
C      Parameter(Zero=0.0D0,One=1.0D0,Two=2.0D0)
CC
C      Addr1=IAddr1
C      Addr2=IAddr2
C      Addr3=IAddr3
C      Addr4=IAddr4
CC
C      NAddr3=Zero
CC
CC     CHANGE THE ORDER IF NECESSARY
CC
C      Addr12=Max(Addr1,Addr2)*(Max(Addr1,Addr2)-1)/2+
C     $       Min(Addr2,Addr1)
C      Addr34=Max(Addr3,Addr4)*(Max(Addr3,Addr4)-1)/2+
C     $       Min(Addr3,Addr4)
CC
CC     GET THE POSITION OF THE ELEMEMT (12|34)
CC
C      NAddr3=Max(Addr12,Addr34)*(Max(Addr12,Addr34)-1)/2+
C     $       Min(Addr12,Addr34)
CC
C      Return
C      End
C
*Deck NAddr3
      Integer(8) Function NAddr3(IAddr1,IAddr2,IAddr3,IAddr4)
C
C     POINTER FOR TWO-ELECTRON INTEGRALS
C     NAddr3 WILL POINT TO THE INTEGRAL (IAddr1,IAddr2,IAddr3,IAddr4)
C
      Integer(8) :: IAddr12,IAddr34
C
C     CHANGE THE ORDER IF NECESSARY
C
      IAddr12=Max(IAddr1,IAddr2)*(Max(IAddr1,IAddr2)-1)/2+
     $       Min(IAddr2,IAddr1)
      IAddr34=Max(IAddr3,IAddr4)*(Max(IAddr3,IAddr4)-1)/2+
     $       Min(IAddr3,IAddr4)
C
C     GET THE POSITION OF THE ELEMEMT (12|34)
C
      NAddr3=Max(IAddr12,IAddr34)*(Max(IAddr12,IAddr34)-1)/2+
     $       Min(IAddr12,IAddr34)
C
      Return
      End

*Deck Diag8
      Subroutine Diag8_old(C,NDim,NVar,AII,AJJ)

      Real*8 C,AII,AJJ
!
!     GENERAL PURPOSE DIAGONALIZATION ROUTINE
!
      Call Tred8(NDim,NVar,AII,AJJ,C)
      Call Tql8(NDim,NVar,AII,AJJ,C,IErr)
      If(IErr.Eq.0) Return

      Stop
      End

      subroutine Diag8(A, nmax, n, eigenvalues, Work_unused)
      ! replaces the old diag8 above
      implicit none
      integer, intent(in)            :: nmax,n
      double precision, intent(inout):: A(nmax,n)
      double precision, intent(out)  :: eigenvalues(n)
      double precision, intent(in)   :: work_unused(*)

      double precision, allocatable  :: work(:,:)
      integer         , allocatable  :: iwork(:)
      integer                        :: lwork, liwork, info
      integer                        :: i,j

      lwork = 1
      liwork = 1
      allocate (work(lwork,1),iwork(liwork))

      ! Determine optimal size for work arrays
      lwork = -1
      liwork = -1
      call DSYEVD( 'V', 'L', n, A, nmax, eigenvalues, work, lwork,
     $             iwork, liwork, info )
      if (info < 0) then
         deallocate (work,iwork)
         print *, 'Diag8: ',
     $     ' DSYEVD: the ',-info,'-th argument had an illegal value'
         stop 2
      endif

      lwork  = max(int(work(1,1)), 2*n*n + 6*n+ 1)
      liwork = max(iwork(1), 5*n + 3)

!     /!\ liwork becomes negative when > 2147483648 (integer*4 overflow)
      if ((liwork < 0) .or. (lwork < 0)) then
         deallocate (work,iwork)
         print *, 'Diag8: ',
     $     ' Required work space too large'
         stop 3
      endif

      ! Allocate temporary arrays
      deallocate (work,iwork)
      allocate (work(lwork,1),iwork(liwork))

      ! Diagonalize
      call DSYEVD( 'V', 'L', n, A, nmax, eigenvalues, work, lwork,
     $      iwork, liwork, info)

      deallocate(work,iwork)

      if (info < 0) then
         deallocate (work,iwork)
         print *, 'Diag8:',
     $      ': DSYEVD: the ',-info,'-th argument had an illegal value'
         stop 2
      else if( info > 0 ) then
         deallocate (work,iwork)
         write(*,*)'DSYEVD Failed'
         stop 1
      end if

      ! Transpose result
      allocate(work(size(A,1),size(A,2)))
      work(:,:) = A(:,:)
      do j=1,n
        do i=1,n
          A(i,j) = work(j,i)
        enddo
      enddo
      deallocate(work)

! Suggestion :
!   Remove the transposition, and use the eigenvectors as column vectors
!   outside of this routine

      end


*Deck Tred8
      Subroutine Tred8(NM,N,D,E,Z)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension D(N),E(N),Z(NM,N)
C
      Parameter(Zero=0.0D0,One=1.0D0)
C
      If(N.Eq.1) GoTo 320
C
      Do 300 II=2,N
      I=N+2-II
      L=I-1
      H=Zero
      Scale=Zero
C
      If(L.Lt.2) GoTo 130
C
      Do 120 K=1,L
  120 Scale=Scale+Abs(Z(K,I))
C
      If(Scale.Ne.Zero) GoTo 140
  130 E(I)=Z(L,I)
      GoTo 290
C
  140 Do 150 K=1,L
      Z(K,I)=Z(K,I)/Scale
  150 H=H+Z(K,I)*Z(K,I)
C
      F=Z(L,I)
      G=-Sign(Sqrt(H),F)
      E(I)=Scale*G
      H=H-F*G
      Z(L,I)=F-G
      F=Zero
      Do 240 J=1,L
      Z(I,J)=Z(J,I)/(Scale*H)
      G=Zero
C
      Do 180 K=1,J
  180 G=G+Z(K,J)*Z(K,I)
C
      JP1=J+1
      If(L.Lt.JP1) GoTo 220
C
      Do 200 K=JP1,L
  200 G=G+Z(J,K)*Z(K,I)
C
  220 E(J)=G/H
      F=F+E(J)*Z(J,I)
  240 Continue
C
      HH=F/(H+H)
C
      Do 260 J=1,L
      F=Z(J,I)
      G=E(J)-HH*F
      E(J)=G
C
      Do 260 K=1,J
  260 Z(K,J)=Z(K,J)-F*E(K)-G*Z(K,I)
C
      Do 280 K=1,L
  280 Z(K,I)=Scale*Z(K,I)
C
  290 D(I)=H
  300 Continue
  320 D(1)=Zero
C
      E(1)=Zero
C
C
      Do 500 I=1,N
      L=I-1
      If(D(I).Eq.Zero) GoTo 380
C
      Do 360 J=1,L
      G=Zero
C
      Do 340 K=1,L
  340 G=G+Z(K,I)*Z(J,K)
C
      Do 360 K=1,L
  360 Z(J,K)=Z(J,K)-G*Z(I,K)
C
  380 D(I)=Z(I,I)
C
      Z(I,I)=One
      If(L.Lt.1) GoTo 500
C
      Do 400 J=1,L
      Z(I,J)=Zero
  400 Z(J,I)=Zero
C
  500 Continue
C
      Return
      End

      Subroutine Tql8(NM,N,D,E,Z,IErr)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension D(N),E(N),Z(NM,N)
      Real*8 MachEp
C
      Parameter(Zero=0.0D0,One=1.0D0)
      Parameter(MaxIt=30,Toler=1.0D-32,MachEp=2.0D0**(-24))
C
      IErr=0
      If(N.Eq.1) GoTo 1001
C
      Do 100 I=2,N
  100 E(I-1)=E(I)
C
      F=Zero
      B=Zero
      E(N)=Zero
C
      Do 240 L=1,N
      J=0
      H=MachEp*(Abs(D(L))+Abs(E(L)))
      If(B.Lt.H) B=H
C
      Do 110 M=L,N
      If(Abs(E(M)).Le.B) GoTo 120
  110 Continue
C
  120 If(M.Eq.L) GoTo 220
C
  130 If(J.Eq.MaxIt) GoTo 1000
      J=J+1
      P=(D(L+1)-D(L))/(E(L)+E(L))
      R=Sqrt(P*P+One)
      H=D(L)-E(L)/(P+Sign(R,P))
C
      Do 140 I=L,N
  140 D(I)=D(I)-H
C
      F=F+H
      P=D(M)
      C=One
      S=Zero
      MML=M-L
C
      Do 200 II=1,MML
      I=M-II
      G=C*E(I)
      H=C*P
      If(Abs(P).Lt.Abs(E(I))) GoTo 150
      C=E(I)/P
      R=Sqrt(C*C+One)
      E(I+1)=S*P*R
      S=C/R
      C=One/R
      GoTo 160
  150 C=P/E(I)
C
      R=Sqrt(C*C+One)
      E(I+1)=S*E(I)*R
      S=One/R
      C=C*S
  160 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
C
      Do 180 K=1,N
      H=Z(I+1,K)
      Z(I+1,K)=S*Z(I,K)+C*H
  180 Z(I,K)=C*Z(I,K)-S*H
C
  200 Continue
C
      E(L)=S*P
      D(L)=C*P
      If(Abs(E(L)).Gt.B) GoTo 130
  220 D(L)=D(L)+F
C
  240 Continue
C
      Do 300 II=2,N
      I=II-1
      K=I
      P=D(I)
C
      Do 260 J=II,N
      If(D(J).Ge.P) GoTo 260
      K=J
      P=D(J)
  260 Continue
C
      D(K)=D(I)
      D(I)=P
      T=One
      If(Z(K,1).Lt.Toler) T=-One
C
      Do 280 J=1,N
      P=Z(I,J)
      Z(I,J)=Z(K,J)*T
  280 Z(K,J)=P
  300 Continue
C
      GoTo 1001
 1000 IErr=L
C
 1001 Return
      End

*Deck GauLeg
      Subroutine GauLeg(x1,x2,x,w,n)
C
C     From "Numerical Recipes, Ch.4.5:
C     Given the lower and upper limits of integration x1,x2 and given n,
C     this routine scales the range of integration from (x1,x2) to
C     (-1,1) and provides abscissas x(1:n) and weights w(1:n)
C     of the Gauss-Legendre n-point quadrature formula
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension x(n),w(n)
      Parameter (EPS=3.d-14)
C
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
C
      Return
      End


















