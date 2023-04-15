*Deck SRLRDynamic
      Subroutine SRLRDynamic(TwoMO,TwoMOLR,URe,Occ,XOne,
     $  OrbGrid,WGrid,NSymMO,NBasis,NInte1,NInte2,NGrid,NDim,NDimKer,
     $  QMAX)
C
C     A ROUTINE FOR EXCITATION ENERGY FROM THE SRDF+LRDMF 
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0,Big=1.99D0, Small=1.D-6)
C
      Include 'commons.inc'
C
      Dimension
     $ URe(NBasis,NBasis),Occ(NBasis),TwoMO(NInte2),TwoMOLR(NInte2),
     $ XOne(NInte1),OrbGrid(NBasis,NGrid),WGrid(NGrid)
C
C     LOCAL ARRAYS
C
      Dimension
     $ ABPLUS(NDim*NDim),ABMIN(NDim*NDim),
     $ CMAT(NDim*NBasis),DMAT(NBasis*NDim),
     $ D1MAT(NBasis*NDim),WMAT(NBasis*NBasis),
     $ IndX(NDim),IndN(2,NDim),
     $ XKer(NDimKer),NSymMO(NBasis),NSymNO(NBasis),
     $ MultpC(15,15),NSymExc(NDim)
C
C     READ THE DIRECT MULTIPLICATION TABLE
C
      Open(10,File=FMultTab)
      Do I=1,MxSym
      Read(10,*)(MultpC(I,J),J=1,I)
      Do J=1,I
      MultpC(J,I)=MultpC(I,J)
      EndDo
      EndDo
      Close(10)
C
C     ESTABLISH THE SYMMETRY OF NO's
C
      Write(6,'(2/,X,"SYMMETRY OF NATURAL ORBITALS")')
      Write(6,'(X,"Orbital",3X,"Symmetry")')
      Do I=1,NBasis
      NSymNO(I)=0

      Do J=1,NBasis
C
      If(Abs(URe(I,J)).Gt.0.1) Then
C
      ISym=NSymMO(J)
      If(NSymNO(I).Eq.0) Then
      NSymNO(I)=ISym
      Else
      If(NSymNO(I).Ne.ISym) 
     $ Write(6,'("Symm of NO cannot be established",I3)')I
      EndIf
C
      EndIf
      EndDo
      Write(6,'(I4,7X,I4)')I,NSymNO(I)
      EndDo
C
C     OBTAIN THE KERNEL, I.E. THE MATRICES A,B,C,D,D1, AND W
C
c herer!!!
c      Call PrimKer(ABPLUS,ABMIN,CMAT,DMAT,D1MAT,WMAT,
c     $ URe,Occ,XOne,TwoMOLR,NBasis,NDim,NInte1,NInte2)
      Call PrimKer2(ABPLUS,ABMIN,CMAT,DMAT,D1MAT,WMAT,
     $ URe,Occ,XOne,TwoMOLR,NBasis,NDim,NInte1,NInte2)
c      Call KerAPSG(ABPLUS,ABMIN,   
c     $ URe,Occ,XOne,TwoMOLR,NBasis,NDim,NInte1,NInte2)
C
      If(IFunSR.Ne.0) Then
C
      Write(6,'(/,"*** Adding a sr-kernel... ***")')
C
C     ADD A CONTRIBUTION FROM THE SR-XC KERNEL 
C
      Call GetKerNO(XKer,Occ,URe,OrbGrid,WGrid,NSymNO,MultpC,
     $ NDimKer,NBasis,NGrid)
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA-1
      IAB=IAB+1
C
      ICD=0
      Do IC=1,NBasis
      Do ID=1,IC
C
      XKer1234=XKer(NAddrrK(IA,IB,IC,ID))
      TwoSR=TwoMO(NAddr3(IA,IB,IC,ID))-TwoMOLR(NAddr3(IA,IB,IC,ID))
C
      If(IC.Eq.ID) Then
C
      IABC=(IC-1)*NDim+IAB
      ICAB=(IAB-1)*NBasis+IC
      CMAT(IABC)=CMAT(IABC)-Two*(Occ(IA)-Occ(IB))*(XKer1234+TwoSR)
      D1MAT(ICAB)=CMAT(IABC)
C
      Else     
C
      ICD=ICD+1
      IABCD=(ICD-1)*NDim+IAB
      ABPLUS(IABCD)=ABPLUS(IABCD)
     $ -Four*(Occ(IA)-Occ(IB))*(Occ(IC)-Occ(ID))*(XKer1234+TwoSR)
C
      EndIf
C
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do I=1,NBasis
      Do J=1,I
C
      XKer1234=XKer(NAddrrK(I,I,J,J))
      IJ=(J-1)*NBasis+I
      JI=(I-1)*NBasis+J
C
      TwoSR=TwoMO(NAddr3(I,I,J,J))-TwoMOLR(NAddr3(I,I,J,J)) 
      WMAT(IJ)=WMAT(IJ)-Two*(XKer1234+TwoSR)
      WMAT(JI)=WMAT(IJ)
C
      EndDo
      EndDo
C
      Write(6,'("*** sr-kernel added. ***")')
C
C     END OF THE CONDITION FOR IFunSR.Ne.0
C
      EndIf 
C
C     REDUCE DIMENSIONS OF THE MATRICES TO TAKE INTO ACCOUNT ONLY
C     NDimX ELEMENTS OF dGamma_ij AND NDimN ELEMENTS OF dn_i
C
      NELE0=0
      NELE1=0
C
C     NELE1 SHOULD BE KEPT AS ZERO, OTHERWISE THE OMEGA->0
C     ASYMPTOTIC OF THE NONADIABATIC APPROXIMATION WILL BE WRONG
C
      NDimN=NBasis-NELE1-NELE0
C
C     CONSTRUCT LOOK-UP TABLES
C
      IJ=0
      Ind=0
      Do I=1,NBasis
      Do J=1,I-1  
      IJ=IJ+1
C
      If(J.Le.QMAX) Then

      If(Occ(I)+Occ(J).Lt.Big.And.Occ(I)+Occ(J).Gt.Small) Then
C
      Ind=Ind+1
      IndX(Ind)=IJ
      IndN(1,Ind)=I
      IndN(2,Ind)=J
C      
      EndIf
      EndIf
C
      EndDo
      EndDo
C
      NDimX=Ind       
C
      Write(6,'(2/,X,"Total number of pairs:",I4)')NDim
      Write(6,'(X,"Reduced to:",I4)')Ind
      Write(6,'(X,"Accepted pairs read:")')
      Do I=1,Ind
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      Write(6,'(2X,2I3,2E14.4)')Ind1,Ind2,Occ(Ind1),Occ(Ind2)
      EndDo
C
C     REDUCE THE MATRICES
C
      Do J=1,NDimX
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(IndX(J)-1)*NDim+IndX(I)
      ABPLUS(IJ)=ABPLUS(IJ1)
      ABMIN(IJ)=ABMIN(IJ1) 
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimX
      IJ=(J-1)*NDimX+I
      IJ1=(J+NELE1-1)*NDim+IndX(I)
      CMAT(IJ)=CMAT(IJ1)    
      EndDo
      EndDo
C
      Do J=1,NDimX
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(IndX(J)-1)*NBasis+I+NELE1
      D1MAT(IJ)=D1MAT(IJ1)  
      DMAT(IJ)=DMAT(IJ1)
      EndDo
      EndDo
C
      Do J=1,NDimN
      Do I=1,NDimN
      IJ=(J-1)*NDimN+I
      IJ1=(J+NELE1-1)*NBasis+I+NELE1
      WMAT(IJ)=WMAT(IJ1)
      EndDo
      EndDo
C
C     EXCITATION ENERGIES FROM THE ADIABATIC APPROX.
C
      Call PrimExcit(0,ABPLUS,ABMIN,CMAT,D1MAT,WMAT,Occ,
     $ IndN,NSymNO,MultpC,NSymExc,NBasis,NDimX,NDimN,NELE1)
C
      Return
      End

*Deck PrimExcit
      Subroutine PrimExcit(IFlag,ABPLUS,ABMIN,CMAT,D1MAT,
     $ WMAT,Occ,IndN,NSymNO,MultpC,NSymExc,NBasis,NDimX,NDimN,NELE1)
C
C     IFlag = 0 Direct Adiabatic Approximation (DAA)
C           = 1 Corrected Adiabatic Approximation (CAA)
C
C     CALCULATE THE EXCITATION ENERGIES FROM THE SYMMETRIZED EQUATIONS
C     IFlag=1 : THE ADIABATIC APPROXIMAATION WITH THE CORRECT Om->0 ASYMPTOTIC IS EMPLOYED
C               Eq.(47) in PCCP 9, 5956 (2007)
C     IFlag=0 : Eq.(47) without C*T^-1*C^T    
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(SMALL=1.D-7,Delta=1.D-4)
      Parameter(Zero=0.D0,Half=0.5D0,One=1.D0,Two=2.D0,Three=3.D0,
     $ Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension
     $ ABPLUS(NDimX,NDimX),ABMIN(NDimX,NDimX),
     $ CMAT(NDimX,NDimN),D1MAT(NDimN,NDimX),
     $ WMAT(NDimN,NDimN),HlpAB(NDimX,NDimX),
     $ Q3((NDimX+NDimN)*(NDimX+NDimN)),WInv(NDimN,NDimN),
     $ ABMINS(NDimX,NDimX),
     $ Occ(NBasis),Temp(NDimN,NDimN),
     $ Eig(NDimX),Work(NDimX),
     $ IndN(2,NDimX),XR(NDimX*NDimX),ZR(NBasis,NDimX),
     $ NSymNO(NBasis),MultpC(15,15),NSymExc(NDimX),
     $ BAS(NDimX,NDimX),IndS(NDimX)
C
C     SYMMETRIZE AND CHANGE THE SIGN OF THE A-B MATRIX
C
      Do I=1,NDimX
      Do J=1,I
      HlpAB(I,J)=-(ABMIN(I,J)+ABMIN(J,I))/Two
      HlpAB(J,I)=HlpAB(I,J)
      EndDo
      EndDo
C
      INeg=0
      Call Diag8(HlpAB,NDimX,NDimX,Eig,Work)
C
      Do I=1,NDimX
      If(Eig(I).Lt.Zero) Then
      Eig(I)=SQRT(-Eig(I))
      INeg=INeg+1
      Else
      Eig(I)=SQRT(Eig(I))
      EndIf
      EndDo
C
      If(INeg.Ne.0) Write(6,'(/,''Negative Eigenvalues of A-B:'',I3,/)')
     $ INeg
C     
      Do I=1,NDimX
      Do J=1,I
      ABMINS(I,J)=Zero
      BAS(I,J)=Zero
      Do K=1,NDimX
      ABMINS(I,J)=ABMINS(I,J)+HlpAB(K,I)*Eig(K)*HlpAB(K,J)
      BAS(I,J)=BAS(I,J)+HlpAB(K,I)/Eig(K)*HlpAB(K,J)
      EndDo
      ABMINS(J,I)=ABMINS(I,J)
      BAS(J,I)=BAS(I,J)
      EndDo
      EndDo
C
C     MULTIPLY Sqrt(B-A) BY N^-1 FROM THE RIGHT
C     
      Do IJ=1,NDimX
      I=IndN(1,IJ)
      J=IndN(2,IJ)
C
      If(Abs(Occ(I)-Occ(J)).Gt.Delta*Occ(I)) Then
C
      Do K=1,NDimX
      ABMINS(K,IJ)=ABMINS(K,IJ)/(Occ(I)-Occ(J))
      EndDo
C
      Else
C
      Do K=1,NDimX
      ABMINS(K,IJ)=Zero
      EndDo
C
      Write(6,'("Set to 0 in N^-1 ",4I3,3X,2E14.6)') 
     $ I,J,NSymNO(I),NSymNO(J),Occ(I),Occ(J)
C
      EndIf
      EndDo
C
      If(IFlag.Eq.1) Then
C
C     INVERT WMAT
C
      Do I=1,NDimN
      Do J=1,I
      Temp(I,J)=(WMAT(I,J)+WMAT(J,I))/Two
      Temp(J,I)=Temp(I,J)
      EndDo
      EndDo
C
      Call Diag8(Temp,NDimN,NDimN,Eig,Work)
C
      EigL=Eig(1)
      Do I=1,NDimN
      If(Abs(Eig(I)).Gt.EigL) EigL=Abs(Eig(I))
      EndDo
C
      Do I=1,NDimN
      If(Abs(Eig(I))/EigL.Gt.SMALL) Then
      Eig(I)=One/Eig(I)
      Else
      Eig(I)=Zero
      Write(*,*)'Zero eigenvalue of the W matrix!'
      EndIf
      EndDo
C
      Do I=1,NDimN
      Do J=1,I
      WInv(I,J)=Zero
      Do K=1,NDimN
      WInv(I,J)=WInv(I,J)+Temp(K,I)*Eig(K)*Temp(K,J)
      EndDo
      WInv(J,I)=WInv(I,J)
      EndDo
      EndDo
C
C     WInv*D1MAT -> Q3
C
      Call MultpMN(Q3,WInv,D1MAT,NDimN,NDimN,NDimN,NDimX)
C
C     CMAT*Q3 -> HlpAB
C
      Call MultpMN(HlpAB,CMAT,Q3,NDimX,NDimN,NDimN,NDimX)
C
C     End of IFlag.Eq.1
C
      EndIf
C
      Do I=1,NDimX
      Do J=1,I
C
      HlpAB(I,J)=-ABPLUS(I,J)+Two*HlpAB(I,J)*Float(IFlag)
      HlpAB(J,I)=HlpAB(I,J)
      EndDo
      EndDo
C
      Call MultpMN(Q3,ABMINS,HlpAB,NDimX,NDimX,NDimX,NDimX)
C
C     COMPUTE TRANSPOSED ABMINS, i.e. N^-1*(B-A)^1/2
C
      Do I=1,NDimX
      Do J=1,I
      Hlp=ABMINS(I,J)
      ABMINS(I,J)=ABMINS(J,I)
      ABMINS(J,I)=Hlp
      EndDo
      EndDo
C
      Call MultpMN(HlpAB,Q3,ABMINS,NDimX,NDimX,NDimX,NDimX)
C
      Do I=1,NDimX
      Do J=1,I
      HlpAB(I,J)=(HlpAB(I,J)+HlpAB(J,I))/Two
      HlpAB(J,I)=HlpAB(I,J)
      EndDo
      EndDo
C
      Write(6,'(/,''Excitation Energies in [au] and [eV]'')')
      NISum=0
      Do NSym=1,MxSym
C
      Write(6,'(2/," ********** SYMMETRY CLASS NO ",I2," **********")')
     $ NSym
      Write(6,'(" Excit No        [au]        [eV] ")')
C
C     Find excitation energies for a given symmetry
C
      Index=0
C
      Do I=1,NDimX
      Ind1=IndN(1,I)
      Ind2=IndN(2,I)
      ISym=MultpC(NSymNO(Ind1),NSymNO(Ind2))
C
      Do J=1,NDimX
      Ind1=IndN(1,J)
      Ind2=IndN(2,J)
      JSym=MultpC(NSymNO(Ind1),NSymNO(Ind2))
C
      If(ISym.Eq.NSym.And.JSym.Eq.NSym) Then
C
      Index=Index+1
      XR(Index)=HlpAB(I,J)
C
      If(Index.Eq.1) IHlp=I
      If(I.Eq.IHlp) IndS(Index)=J
C
      EndIf
C
      EndDo
      EndDo
C
      NI=Sqrt(Float(Index))
      NISum=NISum+NI
      Call Diag8(XR,NI,NI,Eig,Work)
C
      Do I=1,Min(10,NI)
      If(Eig(I).Lt.Zero.And.Abs(Eig(I)).Lt.SMALL) Eig(I)=Zero 
      Write(6,'(I5,5X,2F12.5)')I,Sqrt(Eig(I)),
     $ Sqrt(Eig(I))*27.211
      EndDo
      Write(6,*)
C
C     Find the composition of the XR,XI vectors
C
      Write(6,'(" Ind1 Ind2 Sym1 Sym2    GamR    GamI")')
C
      Do I=1,Min(NI,4)
C
      Write(6,'("__________________________________________")')
      Write(6,'(X,"Excit No",I3)')I
C
      Do J=1,NI
      SumR=Zero
      SumI=Zero
      Do K=1,NI
      IK=NI*(K-1)+I
      SumR=SumR+ABMINS(IndS(J),IndS(K))*XR(IK)
      SumI=SumI+Sqrt(Eig(I))*BAS(IndS(J),IndS(K))*XR(IK)
      EndDo
C
      Ind1=IndN(1,IndS(J))
      Ind2=IndN(2,IndS(J))
C
      SumR=SumR*(Occ(Ind1)-Occ(Ind2))
      SumI=SumI*(Occ(Ind1)-Occ(Ind2))
      

      SumRI=Abs(SumR)+Abs(SumI)
C
      If(SumRI.Gt.0.05) Then
      Write(*,'(4I3,2F12.4)') Ind1,Ind2,NSymNO(Ind1),NSymNO(Ind2),
     $ SumR,SumI
      EndIf
C
      EndDo
      EndDo
C
C     End of NSym loop 
C
      EndDo
C
      If(NISum.Ne.NDimX) Write(*,*) 
     $ 'Sum of Dims of symmetrized matrices different from NDimX'
C
      If(IFlag.Eq.1) Then
C
      Stop 'Enable calculations of the Z vector'
C
C     Get components of Z=-2*T^-1*C^T*XR 
C
      Write(*,*)'Components of the Z vector'
C
      Call MultpMN(Q3,WInv,D1MAT,NDimN,NDimN,NDimN,NDimX)
      Call MultpMN(ZR,Q3,XR,NBasis,NDimX,NDimX,NDimX)

      Sum=Zero
      Do I=1,NDimX
      Write(*,*)'Excit=',Sqrt(Eig(I)),Sqrt(Eig(I))*27.211
      Do J=1,NBasis
      If(Abs(ZR(J,I)).Gt.1.D-3) Write(*,*)J,-Two*ZR(J,I)
      Sum=Sum+ZR(J,I)
      EndDo
      EndDo
C
      Write(*,*)'Sum of Z = ',Sum
C
      EndIf
C
      Return
      End

*Deck PrimKer 
      Subroutine PrimKer(ABPLUS,ABMIN,CMAT,DMAT,D1MAT,WMAT,
     $ URe,Occ,XOne,TwoMO,NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE ADIABATIC KERNEL FOR THE PRIMITIVE FUNCTIONALS 
C     AND STORE ITS PARTS IN ABPLUS, ABMIN, CMAT, DMAT
C     OBTAIN ALSO D1MAT AND WMAT 
C
C     ASSUME THE <II|JJ> INTEGRALS IN THE XC PART OF THE FUNCTIONALS!!!
C
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ CMAT(NDim,NBasis),DMAT(NBasis,NDim),
     $ D1MAT(NBasis,NDim),WMAT(NBasis,NBasis),
     $ URe(NBasis,NBasis),
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),Hlp1(NBasis,NInte1),Hlp2(NBasis,NInte1),
     $ Hlp3(NInte1),AMAT(NDim,NDim),BMAT(NDim,NDim),
     $ temp(Nbasis,NBasis),eig(nbasis),work(nbasis)
C
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
C     ONE-ELECTRON MATRIX IN A NO REPRESENTATION
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     Hlp1(p,qa)=Sum_i F_pi <qa|ii>
C
      Do J=1,NInte1
      Hlp3(J)=Zero
      Do I=1,NBasis
      Hlp1(I,J)=Zero
      Hlp2(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      IQA=0
      Do IQ=1,NBasis
      Do IA=1,IQ
      IQA=IQA+1
C
      TwoInt=TwoMO(NAddr3(IQ,I,IA,I))
C
      Hlp3(IQA)=Hlp3(IQA)+Occ(I)*TwoMO(NAddr3(IQ,IA,I,I))
C
      Do IP=1,NBasis
      Hlp1(IP,IQA)=Hlp1(IP,IQA)+GOCC(Occ(IP),Occ(I),0,IP,I)*TwoInt
      Hlp2(IP,IQA)=Hlp2(IP,IQA)+GOCC(Occ(IP),Occ(I),1,IP,I)*TwoInt
      EndDo
C
      EndDo   
      EndDo
      EndDo
C
C     OBTAIN CMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
      IPQ1=(Max(IP,IQ)*(Max(IP,IQ)-1))/2+Min(IP,IQ)
      Do IA=1,NBasis
C
      CMAT(IPQ,IA)=Zero
C
      If(IA.Eq.IP) 
     $ CMAT(IPQ,IA)=CMAT(IPQ,IA)-HNO(IPQ1)
     $ -Two*XKU*Hlp3(IPQ1)-Hlp2(IA,IPQ1)
C
      If(IA.Eq.IQ) 
     $ CMAT(IPQ,IA)=CMAT(IPQ,IA)+HNO(IPQ1)
     $ +Two*XKU*Hlp3(IPQ1)+Hlp2(IA,IPQ1)
C
      CMAT(IPQ,IA)=CMAT(IPQ,IA)
     $ -Two*XKU*(Occ(IP)-Occ(IQ))*TwoMO(NAddr3(IP,IQ,IA,IA))
     $ -(GOCC(Occ(IA),Occ(IP),1,IA,IP)-GOCC(Occ(IA),Occ(IQ),1,IA,IQ))
     $ *TwoMO(NAddr3(IP,IA,IQ,IA)) 
C
      EndDo
      EndDo
      EndDo
C
C     COMPUTE AMAT,BMAT, AND DMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IQ.Ne.IP) IPQ=IPQ+1
C
      Do IA=1,NBasis
      Do IB=1,NBasis
C
      If(IA.Ne.IB) Then
C
      XKernel=Two*XKU*(Occ(IP)-Occ(IQ))*(Occ(IA)-Occ(IB))
     $ *TwoMO(NAddr3(IP,IQ,IA,IB))
     $ +(GOCC(Occ(IP),Occ(IA),0,IP,IA)+GOCC(Occ(IQ),Occ(IB),0,IQ,IB))
     $ *(TwoMO(NAddr3(IP,IA,IQ,IB))+TwoMO(NAddr3(IP,IB,IQ,IA)))
C
      IQA=(Max(IA,IQ)*(Max(IA,IQ)-1))/2+Min(IA,IQ)
      IQB=(Max(IB,IQ)*(Max(IB,IQ)-1))/2+Min(IB,IQ)
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      Part1=Zero
      If(IA.Eq.IP) Then
C
      Part1=-HNO(IQB)
C
      XKernel=XKernel+Two*XKU*(Occ(IA)-Occ(IB))*Hlp3(IQB)
     $ -Hlp1(IB,IQB)
C
      EndIf
C
      If(IB.Eq.IQ) Then
C
      Part1=Part1+HNO(IPA)
C
      XKernel=XKernel-Two*XKU*(Occ(IA)-Occ(IB))*Hlp3(IPA)
     $ -Hlp1(IA,IPA)
C
      EndIf
C
      If(IA.Eq.IQ) XKernel=XKernel-Hlp1(IQ,IPB)
      If(IB.Eq.IP) XKernel=XKernel-Hlp1(IP,IQA)
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)-Max(IA,IB)+1
C
      XKernel=Part1*(Occ(IA)-Occ(IB))-XKernel
C
      If(IP.Eq.IQ) Then
      DMAT(IP,IAB)=XKernel
      ElseIf(IA.Gt.IB) Then
      AMAT(IPQ,IAB)=XKernel
      Else
      BMAT(IPQ,IAB)=-XKernel
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     GET WMAT
C
      Do IA=1,NBasis
      Do IP=1,NBasis
C
      WMAT(IP,IA)=-Two*XKU*TwoMO(NAddr3(IP,IP,IA,IA))
     $ -GOCC(Occ(IA),Occ(IP),3,IA,IP)*TwoMO(NAddr3(IP,IA,IP,IA))
C
      EndDo
C
      Do IQ=1,NBasis
C
      WMAT(IA,IA)=WMAT(IA,IA)
     $ -GOCC(Occ(IA),Occ(IQ),2,IA,IQ)*TwoMO(NAddr3(IA,IQ,IA,IQ))
C
      EndDo
C
      EndDo
C
C     CONSTRUCT A D1MAT MATRIX
C
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA-1 
      IAB=IAB+1
      Do IP=1,NBasis
C
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      D1MAT(IP,IAB)=Two*XKU*(Occ(IB)-Occ(IA))*TwoMO(NAddr3(IP,IP,IA,IB))
     $+ (GOCC(Occ(IP),Occ(IB),1,IP,IB)-GOCC(Occ(IP),Occ(IA),1,IP,IA))
     $ *TwoMO(NAddr3(IP,IA,IP,IB))
C
      If(IP.Eq.IA) 
     $ D1MAT(IP,IAB)=D1MAT(IP,IAB)-HNO(IPB)
     $ -Two*XKU*Hlp3(IPB)-Hlp2(IP,IPB)
C
      If(IP.Eq.IB) 
     $ D1MAT(IP,IAB)=D1MAT(IP,IAB)+HNO(IPA)
     $ +Two*XKU*Hlp3(IPA)+Hlp2(IP,IPA)
C
      EndDo
      EndDo
      EndDo
C
C     GET A+B AND A-B
C
      Call AddM(AMAT,BMAT,ABPLUS,NDim)
      Call DiffM(AMAT,BMAT,ABMIN,NDim)
C
      Return
      End

*Deck PrimKer2
      Subroutine PrimKer2(ABPLUS,ABMIN,CMAT,DMAT,D1MAT,WMAT,
     $ URe,Occ,XOne,TwoMO,NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE ADIABATIC KERNEL FOR THE PRIMITIVE FUNCTIONALS
C     AND STORE ITS PARTS IN ABPLUS, ABMIN, CMAT, DMAT
C     OBTAIN ALSO D1MAT AND WMAT
C
C     THE EXCHANGE INTEGRALS <IJ|JI> ARE ASSUMED IN THE EX PART 
C     OF THE FUNCTIONALS!
C
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ CMAT(NDim,NBasis),DMAT(NBasis,NDim),
     $ D1MAT(NBasis,NDim),WMAT(NBasis,NBasis),
     $ URe(NBasis,NBasis),   
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),Hlp1(NBasis,NInte1),Hlp2(NBasis,NInte1),
     $ Hlp3(NInte1),Hlp4(NBasis,NInte1),AMAT(NDim,NDim),BMAT(NDim,NDim)
C   
      XKU=One
      If(IFun.Eq.1) XKU=Zero
C
C     ONE-ELECTRON MATRIX IN A NO REPRESENTATION
C   
      IJ=0 
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
      EndDo
      EndDo
C
C     Hlp1(p,qa)=Sum_i F_pi <qa|ii>
C
      Do J=1,NInte1
      Hlp3(J)=Zero
      Do I=1,NBasis
      Hlp1(I,J)=Zero
      Hlp2(I,J)=Zero
      Hlp4(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      IQA=0
      Do IQ=1,NBasis
      Do IA=1,IQ
      IQA=IQA+1
C
      TwoInt=TwoMO(NAddr3(IQ,I,IA,I))
C
      Hlp3(IQA)=Hlp3(IQA)+Occ(I)*TwoMO(NAddr3(IQ,IA,I,I))
      If(IFun.Eq.13) Then
      Do J=1,NBasis
      If(IGem(J).Ne.IGem(I)) Hlp4(J,IQA)=Hlp4(J,IQA)
     $ +Occ(I)*TwoMO(NAddr3(IQ,IA,I,I))
      EndDo
      EndIf
C
      Do IP=1,NBasis
      Hlp1(IP,IQA)=Hlp1(IP,IQA)+GOCC(Occ(IP),Occ(I),0,IP,I)*TwoInt
      Hlp2(IP,IQA)=Hlp2(IP,IQA)+GOCC(Occ(IP),Occ(I),1,IP,I)*TwoInt
      EndDo
C
      EndDo
      EndDo
      EndDo
C
C     OBTAIN CMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
      IPQ1=(Max(IP,IQ)*(Max(IP,IQ)-1))/2+Min(IP,IQ)
      Do IA=1,NBasis
C
      CMAT(IPQ,IA)=Zero
C
      If(IA.Eq.IP)
     $ CMAT(IPQ,IA)=CMAT(IPQ,IA)-HNO(IPQ1)
     $ -Two*XKU*Hlp3(IPQ1)-Hlp2(IA,IPQ1)
C
      If(IA.Eq.IQ)
     $ CMAT(IPQ,IA)=CMAT(IPQ,IA)+HNO(IPQ1)
     $ +Two*XKU*Hlp3(IPQ1)+Hlp2(IA,IPQ1)
C
      CMAT(IPQ,IA)=CMAT(IPQ,IA)
     $ -Two*XKU*(Occ(IP)-Occ(IQ))*TwoMO(NAddr3(IP,IQ,IA,IA))
     $ -(GOCC(Occ(IA),Occ(IP),1,IA,IP)-GOCC(Occ(IA),Occ(IQ),1,IA,IQ))
     $ *TwoMO(NAddr3(IP,IA,IQ,IA))
C
      EndDo
      EndDo
      EndDo
C
C     COMPUTE AMAT,BMAT, AND DMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP
      If(IQ.Ne.IP) IPQ=IPQ+1
C
      Do IA=1,NBasis
      Do IB=1,NBasis
C
      If(IA.Ne.IB) Then
C
      XKernel=(GOCC(Occ(IP),Occ(IA),0,IP,IA)
     $ +GOCC(Occ(IQ),Occ(IB),0,IQ,IB)
     $ - GOCC(Occ(IP),Occ(IB),0,IP,IB)-GOCC(Occ(IQ),Occ(IA),0,IQ,IA))
     $ *TwoMO(NAddr3(IP,IA,IQ,IB))
C
      If(IFun.Ne.13) Then 
      XKernel=XKernel+Two*XKU*(Occ(IP)-Occ(IQ))*(Occ(IA)-Occ(IB))
     $ *TwoMO(NAddr3(IP,IQ,IA,IB))
      Else
      GPA=Zero
      If(IGem(IP).Ne.IGem(IA)) GPA=One
      GPB=Zero
      If(IGem(IP).Ne.IGem(IB)) GPB=One
      GQA=Zero
      If(IGem(IQ).Ne.IGem(IA)) GQA=One
      GQB=Zero
      If(IGem(IQ).Ne.IGem(IB)) GQB=One
      XKernel=XKernel+Two*
     $ (GPA*Occ(IP)*Occ(IA)-GPB*Occ(IP)*Occ(IB)
     $ -GQA*Occ(IQ)*Occ(IA)+GQB*Occ(IQ)*Occ(IB))
     $ *TwoMO(NAddr3(IP,IQ,IA,IB))
      EndIf
C
      IQA=(Max(IA,IQ)*(Max(IA,IQ)-1))/2+Min(IA,IQ)
      IQB=(Max(IB,IQ)*(Max(IB,IQ)-1))/2+Min(IB,IQ)
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      Part1=Zero
      If(IA.Eq.IP) Then
C
      Part1=-HNO(IQB)
C
      XKernel=XKernel-Hlp1(IB,IQB)+Hlp1(IA,IQB)
C
      If(IFun.Ne.13) Then
      XKernel=XKernel+Two*XKU*(Occ(IA)-Occ(IB))*Hlp3(IQB)
      Else
      XKernel=XKernel+Two*(Occ(IA)*Hlp4(IA,IQB)-Occ(IB)*Hlp4(IB,IQB))
      EndIf
C
      EndIf
C
      If(IB.Eq.IQ) Then
C
      Part1=Part1+HNO(IPA)
C
      XKernel=XKernel-Hlp1(IA,IPA)+Hlp1(IB,IPA)
C
      If(IFun.Ne.13) Then
      XKernel=XKernel-Two*XKU*(Occ(IA)-Occ(IB))*Hlp3(IPA)
      Else
      XKernel=XKernel-Two*(Occ(IA)*Hlp4(IA,IPA)-Occ(IB)*Hlp4(IB,IPA))
      EndIf
C
      EndIf
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)-Max(IA,IB)+1
C
      XKernel=Part1*(Occ(IA)-Occ(IB))-XKernel
C
      If(IP.Eq.IQ) Then
      DMAT(IP,IAB)=XKernel
      ElseIf(IA.Gt.IB) Then
      AMAT(IPQ,IAB)=XKernel
      Else
      BMAT(IPQ,IAB)=-XKernel
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      EndDo
      EndDo
C
C     GET WMAT
C
      Do IA=1,NBasis
      Do IP=1,NBasis
C
      WMAT(IP,IA)=-Two*XKU*TwoMO(NAddr3(IP,IP,IA,IA))
     $ -GOCC(Occ(IA),Occ(IP),3,IA,IP)*TwoMO(NAddr3(IP,IA,IP,IA))
C
      EndDo
C
      Do IQ=1,NBasis
C
      WMAT(IA,IA)=WMAT(IA,IA)
     $ -GOCC(Occ(IA),Occ(IQ),2,IA,IQ)*TwoMO(NAddr3(IA,IQ,IA,IQ))
C
      EndDo
C
      EndDo
C
C
C     CONSTRUCT A D1MAT MATRIX
C
      Do IP=1,NBasis
      IAB=0
      Do IA=1,NBasis
      Do IB=1,IA-1
      IAB=IAB+1
C
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      D1MAT(IP,IAB)=Two*XKU*(Occ(IB)-Occ(IA))*TwoMO(NAddr3(IP,IP,IA,IB))
     $+ (GOCC(Occ(IP),Occ(IB),1,IP,IB)-GOCC(Occ(IP),Occ(IA),1,IP,IA))
     $ *TwoMO(NAddr3(IP,IA,IP,IB))
C
      If(IP.Eq.IA)
     $ D1MAT(IP,IAB)=D1MAT(IP,IAB)-HNO(IPB)
     $ -Two*XKU*Hlp3(IPB)-Hlp2(IP,IPB)
C
      If(IP.Eq.IB)
     $ D1MAT(IP,IAB)=D1MAT(IP,IAB)+HNO(IPA)
     $ +Two*XKU*Hlp3(IPA)+Hlp2(IP,IPA)
C
      EndDo
      EndDo
      EndDo
C
C     GET A+B AND A-B
C
      Call AddM(AMAT,BMAT,ABPLUS,NDim)
      Call DiffM(AMAT,BMAT,ABMIN,NDim)
C
      Return
      End

*Deck KerAPSG
      Subroutine KerAPSG(ABPLUS,ABMIN,
     $ URe,Occ,XOne,TwoMO,NBasis,NDim,NInte1,NInte2)
C
C     COMPUTE THE ADIABATIC KERNEL FOR THE APSG FUNCTIONAL
C     AND STORE ITS PARTS IN ABPLUS, ABMIN
C
C     THE EXCHANGE INTEGRALS <IJ|JI> ARE ASSUMED IN THE EX 
C     (DIFFERENT GEMINALS) PART AND <II|JJ> IN THE CORRELATION 
C     (SAME GEMINALS) PART 
C     SHOULD BE EQUIAVELENT TO ERPA WITH APSG
C
      Implicit Real*8 (A-H,O-Z)
C
      Parameter(Zero=0.D0, Half=0.5D0, One=1.D0, Two=2.D0, Four=4.D0)
C
      Include 'commons.inc'
C
      Dimension ABPLUS(NDim,NDim),ABMIN(NDim,NDim),
     $ URe(NBasis,NBasis),   
     $ Occ(NBasis),XOne(NInte1),TwoMO(NInte2)
C
C     LOCAL ARRAYS
C
      Dimension HNO(NInte1),Hlp1(NBasis,NInte1),
     $ Hlp4(NBasis,NInte1),Hlp2(NBasis,NInte1),
     $ AMAT(NDim,NDim),BMAT(NDim,NDim),C(NBasis)
C
      Do I=1,NBasis
      C(I)=CICoef(I)
      EndDo
C   
C     ONE-ELECTRON MATRIX IN A NO REPRESENTATION
C   
      IJ=0 
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      HNO(IJ)=Zero
C
      Do IA=1,NBasis
      Do IB=1,NBasis
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)
      HNO(IJ)=HNO(IJ)+URe(I,IA)*URe(J,IB)*XOne(IAB)
      EndDo
      EndDo
      EndDo
      EndDo
C
      Do J=1,NInte1
      Do I=1,NBasis
      Hlp1(I,J)=Zero
      Hlp2(I,J)=Zero
      Hlp4(I,J)=Zero
      EndDo
      EndDo
C
      Do I=1,NBasis
      IQA=0
      Do IQ=1,NBasis
      Do IA=1,IQ
      IQA=IQA+1
C
      TwoInt=TwoMO(NAddr3(IQ,I,IA,I))
C
      Do J=1,NBasis
C
      If(IGem(J).Ne.IGem(I)) Then
      Hlp4(J,IQA)=Hlp4(J,IQA)
     $ +Occ(I)*TwoMO(NAddr3(IQ,IA,I,I))
      Hlp1(J,IQA)=Hlp1(J,IQA)-Occ(I)*TwoInt
      Else
      Hlp2(J,IQA)=Hlp2(J,IQA)+C(I)*TwoInt
      EndIf
C
      EndDo
C
      EndDo
      EndDo
      EndDo
C
C     COMPUTE AMAT AND BMAT
C
      IPQ=0
      Do IP=1,NBasis
      Do IQ=1,IP-1
      IPQ=IPQ+1
C
      Do IA=1,NBasis
      Do IB=1,NBasis
C
      If(IA.Ne.IB) Then
C
      GPA=Zero
      If(IGem(IP).Eq.IGem(IA)) GPA=One
      GQB=Zero
      If(IGem(IQ).Eq.IGem(IB)) GQB=One
      XKernel=(GPA*C(IP)*C(IA)+GQB*C(IQ)*C(IB))
     $ *(TwoMO(NAddr3(IP,IA,IQ,IB))+TwoMO(NAddr3(IP,IB,IQ,IA)))
C
      GPA=Zero
      If(IGem(IP).Ne.IGem(IA)) GPA=One
      GPB=Zero
      If(IGem(IP).Ne.IGem(IB)) GPB=One
      GQA=Zero
      If(IGem(IQ).Ne.IGem(IA)) GQA=One
      GQB=Zero
      If(IGem(IQ).Ne.IGem(IB)) GQB=One
      XKernel=XKernel+
     $ (GPA*Occ(IP)*Occ(IA)-GPB*Occ(IP)*Occ(IB)
     $ -GQA*Occ(IQ)*Occ(IA)+GQB*Occ(IQ)*Occ(IB))
     $ *(Two*TwoMO(NAddr3(IP,IQ,IA,IB))-TwoMO(NAddr3(IP,IA,IQ,IB)))
C
      IQA=(Max(IA,IQ)*(Max(IA,IQ)-1))/2+Min(IA,IQ)
      IQB=(Max(IB,IQ)*(Max(IB,IQ)-1))/2+Min(IB,IQ)
      IPA=(Max(IA,IP)*(Max(IA,IP)-1))/2+Min(IA,IP)
      IPB=(Max(IB,IP)*(Max(IB,IP)-1))/2+Min(IB,IP)
C
      If(IA.Eq.IP) XKernel=XKernel-C(IB)*Hlp2(IB,IQB)
      If(IB.Eq.IQ) XKernel=XKernel-C(IA)*Hlp2(IA,IPA)
      If(IA.Eq.IQ) XKernel=XKernel-C(IQ)*Hlp2(IQ,IPB)
      If(IB.Eq.IP) XKernel=XKernel-C(IP)*Hlp2(IP,IQA)
C
      Part1=Zero
      If(IA.Eq.IP) Then
C
      Part1=-HNO(IQB)
C
      XKernel=XKernel-Occ(IB)*Hlp1(IB,IQB)+Occ(IA)*Hlp1(IA,IQB)
     $ +Two*(Occ(IA)*Hlp4(IA,IQB)-Occ(IB)*Hlp4(IB,IQB))
C
      EndIf
C
      If(IB.Eq.IQ) Then
C
      Part1=Part1+HNO(IPA)
C
      XKernel=XKernel-Occ(IA)*Hlp1(IA,IPA)+Occ(IB)*Hlp1(IB,IPA)
     $ -Two*(Occ(IA)*Hlp4(IA,IPA)-Occ(IB)*Hlp4(IB,IPA))
C
      EndIf
C
      IAB=(Max(IA,IB)*(Max(IA,IB)-1))/2+Min(IA,IB)-Max(IA,IB)+1
C
      XKernel=Part1*(Occ(IA)-Occ(IB))-XKernel
C
      If(IA.Gt.IB) Then
      AMAT(IPQ,IAB)=XKernel
      Else
      BMAT(IPQ,IAB)=-XKernel
      EndIf
C
      EndIf
C
      EndDo
      EndDo
C
      EndDo
      EndDo
c herer!!!
c      open(10,file="ab.dat")
c      do i=1,ndim
c      do j=1,ndim
c      read(10,*)amat(i,j),bmat(i,j) 
c      if(AMAT(i,j).ne.zero.or.BMAT(i,j).ne.zero)
c     $ write(*,*)i,j,AMAT(i,j),BMAT(i,j)
c      enddo
c      enddo
c      close(10)
C
C
C     GET A+B AND A-B
C
      Call AddM(AMAT,BMAT,ABPLUS,NDim)
      Call DiffM(AMAT,BMAT,ABMIN,NDim)
C
      Return
      End









 
      
