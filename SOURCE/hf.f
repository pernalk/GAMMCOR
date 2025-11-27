*Deck HF
      Subroutine HF(URe,EOrb,ENuc,Occ,XKin,XNuc,TwoEl,TwoElSR,NSymMO,
     $ NInte1,NInte2,NBasis,NOcc,IPrint,IEps)
C
      Implicit Real*8 (A-H,O-Z)
C
C     CARRIES OUT ORDINARY FIRST-ORDER HF CALCULATIONS
C
C     IEps = 1 - OBTAIN ORBITAL ENERGIES AND RETURN
C            0 - HF CALCULATIONS
C
      Parameter(TolOrE=5.0D-8,DumpF=0.25D0,MxItSC=30)
C
      Dimension URe(NBasis,NBasis),EOrb(NBasis),
     $          XKin(NInte1),XNuc(NInte1),TwoEl(NInte2),
     $          TwoElSR(NInte2),
     $          NSymMO(NBasis),Occ(NBasis)
C
C     LOCAL ARRAYS 
C
      Dimension P(NInte1)
C
      Real*8 F(NBasis,NBasis),EigVal(NBasis),Work(NBasis)
C
      Parameter(Zero=0.0D0,Half=0.5D0,One=1.0D0,Two=2.0D0)
C
      NAddr=0
      IPR=0
      Do IP=1,NBasis
      Do IR=1,IP
      IPR=IPR+1
      IQS=0
      Do IQ=1,NBasis
      Do IS=1,IQ
      IQS=IQS+1
      If(IPR.Ge.IQS) Then
      NAddr=NAddr+1
C
c      If(Occ(IP)+Occ(IQ)+Occ(IR)+Occ(IS).lt.3) Then
c      TwoEl(NAddr)=TwoEl(NAddr)+TwoElSR(NAddr)
c      EndIf
C set integrals with core indices to zero
c      If(IP.Eq.1.Or.IQ.Eq.1.Or.IR.Eq.1.Or.IS.Eq.1) TwoElSR(NAddr)=Zero
C
      EndIf
      EndDo
      EndDo
      EndDo
      EndDo
C
c      E2=Zero
c      do i=1,nbasis
c      do j=1,nbasis
c      E2=E2+occ(i)*occ(j)*(
c     $ 2.*TwoEl(NAddr3(i,i,j,j))-TwoEl(NAddr3(i,j,i,j)))
c      enddo
c      enddo
c      write(*,*)'Two-el energy 01',E2
C
C     INITIALIZE DENSITY MATRIX
C
      Do 1 I=1,NInte1
    1 P(I)=Zero
c      Do 2 I=1,NOcc
c      II=I*(I+1)/2
c    2 P(II)=One 
      Do 2 I=1,NBasis
      II=I*(I+1)/2
    2 P(II)=Occ(I)  
C
C     START OF THE SCF LOOP
C
      Iter=0
  100 Iter=Iter+1
      IFlag=0
C
C     CONSTRUCT P AND F 
C
      If(Iter.Gt.1) Then
C
      IJ=0
      Do 50 I=1,NBasis
      Do 50 J=1,I
      IJ=IJ+1
      Help=Zero
      Do 55 K=1,NOcc
   55 Help=Help+F(K,I)*F(K,J)
   50 P(IJ)=(One-DumpF)*Help+DumpF*P(IJ)  
C
      EndIf
C
C     ONE-ELECTRON CONTRIBUTION TO F
C
      IJ=0
      Do 20 I=1,NBasis
      Do 20 J=1,I
      IJ=IJ+1
   20 F(I,J)=XKin(IJ)+XNuc(IJ)    
C
C     TWO-ELECTRON CONTRIBUTION TO F
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
      TwoZet=Two*TwoEl(NAddr)
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
      F(I,J)=F(I,J)+FacIJ*P(KL)*TwoZet
      If(IJ.Ne.KL) F(K,L)=F(K,L)+FacKL*P(IJ)*TwoZet      
      F(I,L)=F(I,L)-Half*FacIL*P(JK)*TwoZet
C
      If(IL.Ne.JK) Then
C
      If(J.Ge.K) Then 
      F(J,K)=F(J,K)-Half*FacJK*P(IL)*TwoZet
      Else
      F(K,J)=F(K,J)-Half*FacJK*P(IL)*TwoZet
      EndIf
C
      EndIf
C
      If(K.Eq.L) GoTo 40
C
      If(I.Eq.J) GoTo 40
C
      F(I,J)=F(I,J)+P(KL)*TwoZet
      If(IJ.Ne.KL) F(K,L)=F(K,L)+P(IJ)*TwoZet
      F(I,K)=F(I,K)-Half*FacIK*P(JL)*TwoZet
C
      If(IK.Ne.JL) Then
C
      If(J.Ge.L) Then
      F(J,L)=F(J,L)-Half*FacJL*P(IK)*TwoZet
      Else
      F(L,J)=F(L,J)-Half*FacJL*P(IK)*TwoZet
      EndIf
C
      EndIf      
C                                                  
   40 Continue
   30 Continue
C
      If(IEps.Eq.1) Then
C
      IJ=0
      Do I=1,NBasis
      Do J=1,I
      IJ=IJ+1
      If(I.Eq.J) EOrb(I)=F(I,I)
      EndDo
      EndDo
C
      Call SortHF(EOrb,Occ,URe,NSymMO,NBasis) 
C
      Return
C
      EndIf
C
C     SYMMETRIZE F AND COMPUTE THE TOTAL ENERGY
C
      ETot=Zero
      IJ=0
C
      Do 45 I=1,NBasis
      Do 45 J=1,I
C
      IJ=IJ+1
C
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
C
      F(J,I)=F(I,J)
   45 ETot=ETot+FacIJ*P(IJ)*(XKin(IJ)+XNuc(IJ)+F(I,J))
C
      Write(6,'(/,'' SCF: ITERATION '',I3,'' TOTAL ENERGY '',E21.14)')
     $ Iter,ETot+ENuc
C
C     DIAGONALIZE F
C
      Call Diag8(F,NBasis,NBasis,EigVal,Work)
C
      If(IPrint.Ge.4) Then
      Write(6,'(/,'' ORBITAL ENERGIES '')')
C
      Do 776 I=1,NOcc
  776 Write(6,'(1X,I3,4X,E21.14)') I,EigVal(I)
      Write(6,*)
      Do 777 I=NOcc+1,NBasis
  777 Write(6,'(1X,I3,4X,E21.14)') I,EigVal(I)
      EndIf
C
      Do 90 I=1,NBasis
      If(Abs(EigVal(I)-EOrb(I)).Ge.TolOrE*Abs(EigVal(I))) IFlag=1
   90 EOrb(I)=EigVal(I)
C
      If((Iter.Lt.MxItSC).And.(IFlag.Eq.1)) GoTo 100
C
      If(Iter.Eq.MxItSC) Stop 'FATAL ERROR: NO CONVERGENCE IN SCF!'
C
      If((IPrint.Eq.1).Or.(IPrint.Eq.2))
     $ Write(6,'(/,'' SCF: Convergence attained after '',I3,
     $ '' iterations'')') Iter
C
  999 Continue
C
C     COMPUTE THE FINAL P AND THE ENERGIES
C
      EKin=Zero
      ETot=Zero
C
      IJ=0
      Do 110 I=1,NBasis
      Do 110 J=1,I
C
      FacIJ=Two
      If(I.Eq.J) FacIJ=One
C
      IJ=IJ+1   
      Help=Zero
C
      Do 120 K=1,NOcc
  120 Help=Help+F(K,I)*F(K,J)
      P(IJ)=Help
      ETot=ETot+FacIJ*P(IJ)*(XKin(IJ)+XNuc(IJ))
  110 EKin=EKin+FacIJ*P(IJ)*XKin(IJ)
C
      EKin=Two*EKin
C
      Do 130 K=1,NOcc
  130 ETot=ETot+EOrb(K)       
C
      If(IPrint.Ge.1) 
     $ Write(6,'(/,'' TOTAL ENERGY   '',E21.14)') ETot+ENuc
C
      If(IPrint.Ge.2) Then
      Write(6,'(/,'' ORBITAL ENERGIES '')')
C
      Do 778 I=1,NOcc
  778 Write(6,'(1X,I3,4X,E21.14)') I,EigVal(I)
      Write(6,*)
      Do 779 I=NOcc+1,NBasis
  779 Write(6,'(1X,I3,4X,E21.14)') I,EigVal(I)
      EndIf
C 
      ECBSPT=Zero
      ECBSAC1=Zero
C
      Do I=1,NOcc
      Do IA=NOcc+1,NBasis
C
      SumJ=Zero
      SumJ1=Zero
      Do J=1,NOcc
C
      SumJ=SumJ+Two*TwoElSR(NAddr3(IA,I,J,J))
     $ -TwoElSR(NAddr3(IA,J,I,J)) 
C
      SumJ1=SumJ1+TwoEl(NAddr3(IA,J,I,J))
C
      EndDo
C
      SumJ1=SumJ1*SumJ
      SumJ=SumJ**2 
      ECBSPT=ECBSPT+SumJ/(EOrb(I)-EOrb(IA))
      ECBSAC1=ECBSAC1+SumJ1/(EOrb(I)-EOrb(IA))
C
      EndDo
      EndDo
      ECBSPT=Two*ECBSPT
      Write(6,'(/,'' HF, ECBS_PT(S), HF+ECBS_PT(S)     '',3F16.8)') 
     $ ETot+ENuc,ECBSPT,ECBSPT+ETot+ENuc
       ECBSAC1=Two*ECBSAC1+ECBSPT
       Write(6,'(/,'' HF, ECBS_AC(S1), HF+ECBS_AC(S1)   '',3F16.8)')
     $ ETot+ENuc,ECBSAC1,ECBSAC1+ETot+ENuc
C
c      Do 900 I=1,NBasis
c      Do 900 J=1,NBasis
c  900 URe(J,I)=F(J,I)
C

C
      Return
      End

*Deck SortHF
      Subroutine SortHF(Eps,Occ,URe,NSymMO,NBasis)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension Eps(NBasis),Occ(NBasis),URe(NBasis,NBasis),
     $ NSymMO(NBasis)
C
C     LOCAL ARRAYS
C
      Dimension EOcc(NBasis),UReOld(NBasis,NBasis),Ind(1000),
     $ IndOcc(1000)
C
C     COPY NEGATIVE EPSILONS TO EOcc
C
      Do I=1,NBasis
      EOcc(I)=-Eps(I)
      EndDo
C
C     SORT THE OCCUPATION NUMBERS IN A DESCENDING ORDER
C
      Do I=1,NBasis
      Ind(I)=I
      EndDo
C
      IStart=1
C
      Do I=1,NBasis
C
      OccMx=EOcc(IStart)
      IndMx=IStart
C
      Do J=IStart,NBasis
      If(EOcc(J).Gt.OccMx) Then
      OccMx=EOcc(J)
      IndMx=J
      EndIf
      EndDo
C
      Hlp=EOcc(IStart)
      IndHlp=Ind(IStart)

      EOcc(IStart)=OccMx
      Ind(IStart)=Ind(IndMx)

      EOcc(IndMx)=Hlp
      Ind(IndMx)=IndHlp
C
      IStart=IStart+1
C
      EndDo
C
C     SWAP THE ORBITALS
C
      Do I=1,NBasis
      Do J=1,NBasis
      UReOld(J,I)=URe(J,I)
      EndDo
      EndDo
C
      Do J=1,NBasis
      Do I=1,NBasis
      URe(I,J)=UReOld(Ind(I),J)
      EndDo
      EndDo
C 
C     CHANGE THE ORDER OF EPSILONS, OCC, AND SYMMETRY INDICES
C
      Do I=1,NBasis
      UreOld(1,I)=Eps(I)
      EndDo
      Do I=1,NBasis
      Eps(I)=UreOld(1,Ind(I))
      EndDo
C
      Do I=1,NBasis
      UreOld(1,I)=Occ(I)
      EndDo
      Do I=1,NBasis
      Occ(I)=UreOld(1,Ind(I))
      EndDo
C
      Do I=1,NBasis
      UreOld(1,I)=NSymMO(I)
      EndDo
      Do I=1,NBasis
      NSymMO(I)=UreOld(1,Ind(I))
      EndDo
C
      Return
      End

