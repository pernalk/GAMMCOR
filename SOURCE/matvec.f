*Deck AddM
      Subroutine AddM(X,Y,Z,N)
C
C     ADD X TO Y: X+Y=Z
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension X(N*N),Y(N*N),Z(N*N)
C
      Do 10 I=1,N*N
   10 Z(I)=X(I)+Y(I) 
C
      Return
      End

*Deck DiffM 
      Subroutine DiffM(X,Y,Z,N)
C
C     TAKE A DIFFERENCE OF X AND Y: X-Y=Z
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension X(N*N),Y(N*N),Z(N*N)
C
      Do 10 I=1,N*N
   10 Z(I)=X(I)-Y(I)
C
      Return
      End

*Deck CpyV
      Subroutine CpyV(X,Y,N)
C
C     COPY A VECTOR Y TO X
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension X(N),Y(N)
C
      Do 10 I=1,N
   10 X(I)=Y(I)
C
      Return
      End
     
*Deck CpyM
      Subroutine CpyM(X,Y,N)
C
C     COPY A MATRIX Y TO X
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension X(N*N),Y(N*N)
C
      Do 10 I=1,N*N
   10 X(I)=Y(I)
C
      Return
      End

*Deck CpySym
      Subroutine CpySym(X,Y,N)
C
C     COPY A SYMMETRIC MATRIX Y TO X
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension X(N,N),Y(N*(N+1)/2)
C
      IJ=0
      Do 10 I=1,N
      Do 10 J=1,I
      IJ=IJ+1
      X(I,J)=Y(IJ)
   10 X(J,I)=Y(IJ)
C
      Return
      End

*Deck VecTr
      Subroutine VecTr(X,V,U,N)
C
C     MULTIPLY A VECTOR V BY U AND U+
C     X = U+ V U  ,
C     WHERE X IS TRIANGULAR
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0)
C
      Dimension X(N*(N+1)/2),V(N),U(N,N)
C
      IJ=0
      Do 10 I=1,N
      Do 10 J=1,I
      IJ=IJ+1
      X(IJ)=Zero
      Do 10 K=1,N
   10 X(IJ)=X(IJ)+V(K)*U(K,I)*U(K,J)
C
      Return
      End

*Deck MatTr
      Subroutine MatTr(X,U,N)
C
C     TRANSFORM A SYMMETRIC MATRIX X WITH U
C     X -> U X U+
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.0D0)
C
      Dimension X(N*(N+1)/2),U(N,N),Hlp(N,N)
C
      Do 10 I=1,N
      Do 10 J=1,N
   10 Hlp(J,I)=Zero
C
      IJ=0
      Do 20 I=1,N
      Do 20 J=1,I
      IJ=IJ+1
      XIJ=X(IJ)
C
      Do 20 K=1,N
      Hlp(K,J)=Hlp(K,J)+U(K,I)*XIJ
      If(I.Ne.J) Hlp(K,I)=Hlp(K,I)+U(K,J)*XIJ
   20 Continue
C
      IJ=0
      Do 30 I=1,N
      Do 30 J=1,I
      IJ=IJ+1
      X(IJ)=Zero
      Do 30 K=1,N
   30 X(IJ)=X(IJ)+Hlp(I,K)*U(J,K)
C
      Return
      End

*Deck MultpM
      Subroutine MultpM(X,Y,U,N)
C
C     MATRIX MULTIPLICATION
C     X = Y*U
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0)
C
      Dimension X(N*N),Y(N*N),U(N*N)
C
      Do 10 I=1,N*N
   10 X(I)=Zero
C
      Do 20 I=1,N
      II=(I-1)*N
C
      Do 20 J=1,N
      JI=II+J
C
      Do 20 K=1,N
      KI=(I-1)*N+K
      JK=(K-1)*N+J
   20 X(JI)=X(JI)+Y(JK)*U(KI)
C
      Return
      End

*Deck MultpMN
      Subroutine MultpMN(X,Y,U,N1,M1,M2,N2)
C
C     MATRIX MULTIPLICATION
C     X = Y*U
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0)
C
      Dimension X(N1*N2),Y(N1*M1),U(M2*N2)
C
      If(M1.Ne.M2) Stop'Fatal error in MultpMN!'
C
      Do 10 I=1,N1*N2
   10 X(I)=Zero
C
      Do 20 I=1,N2
      II=(I-1)*N1
C
      Do 20 J=1,N1
      JI=II+J
C
      Do 20 K=1,M1
      KI=(I-1)*M1+K
      JK=(K-1)*N1+J
   20 X(JI)=X(JI)+Y(JK)*U(KI)
C
      Return
      End
 
