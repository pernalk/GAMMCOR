C
      Character*60 FMultTab
C
      Common/MISCFL/ XELE,IFun,IFunSR,IFunSRKer,NELE,IRes,INO,IPrint,
     $ IInteg,IAO,IGVB,IDALTON,IFreeze,ILoc,NInActOrb,NActOrb,NVirtOrb,
     $ ISAPT,IWarn,ITwoEl,ITrpl,
     $ ICholesky,ICholeskyAccu,MemVal,MemType
      Common/CAS/ ICASSCF,NAcCAS,NInAcCAS,NStates,
     $ InSt(2,100),InTrSt(2,1),ISpinMs2
      Common/BBC3/ IType(5000)
      Common/DFTSR/ Alpha,Cmix
      Common/Sym/ FMultTab,NoSym,MxSym
      Common/APSG/ IGem(5000),CICoef(5000)
      Common/EMBEDD/ IFrag,NGemSave,IAuxGem(5000),
     $ IGemSave(5000), IOrbA(5000),NGemB, NActive,
     $ IConnect(5000,5000)
      Common/AC/ IFlAC,IFlSnd,IFlAC0D,IFlAC0DP,IFlCore,IFlFrag1,IFl12,
     $ NoSt,ISymmAC0D,IFlACFREQ,IFlACFREQNTH,IFlAC1FREQNTH,Max_Cn,
     $ IFlRESPONSE,FreqOm,IRedVirt
      Common/THRGEM/ ThrGemAct
      Common/THR/ ThrSelAct,ThrAct,ThrQVirt,ThrQInact
C
