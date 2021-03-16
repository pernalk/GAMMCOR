#! /usr/bin/env python
import os,sys,string,math
import numpy as np
import subprocess
from numpy.testing import assert_allclose
from numpy.testing import assert_approx_equal

## user defined 
#RUN_GAMMCOR = "/home/kasia/GAMMCOR/gammcor"
RUN_GAMMCOR = "/home/kasia/GAMMCOR/prdmft.exe

### functions reading components of SAPT
## SAPT ENERGY TAGS
ELST_TAG = "E1elst"
E1EX_TAG = "E1exch(S2)"
IND_TAG  = "E2ind "
DSP_TAG  = "E2disp "
DCS_TAG  = "E2disp(CAS)"
DRV_TAG  = "E2disp(SCALED)  "
WI1_TAG  = "Wind_12"
WI2_TAG  = "Wind_13"
WD1_TAG  = "Wdisp_12"
WD2_TAG  = "Wdisp_13"
EXI_TAG  = "E2exch-ind "
EXD_TAG  = "E2exch-disp"
XDRV_TAG = "E2exdisp(SCALED)"
EINT_TAG = "Eint(SAPT2)"
SLVL_TAG = "SAPT level"

IND0_TAG = "E2ind(unc)"
DSP0_TAG = "E2disp(unc)"
SP_TAG   = "E2disp(sp)"
SC_TAG   = "E2disp(sc)"
EXD0_TAG = "E2exch-disp(unc)"
EXI0_TAG = "E2exch-ind(unc)"
INT0_TAG = "Eint(SAPT0)"
#
def get_jobname(log_path,log_name):
   output = os.path.join(log_path,log_name)
   f = open(output)
   for line in f:
      line_adjl = line.lstrip()
      if line_adjl.startswith("JobTitle"):
         job_name = line_adjl.replace("JobTitle","").lstrip()

   return job_name
#
def get_sapt2_energies(log_path,log_name):
    output = os.path.join(log_path,log_name)
    e_elst = 0.0
    e_exch = 0.0
    e_ind  = 0.0
    e_exi  = 0.0
    e_disp = 0.0
    e_exd  = 0.0
    e_int  = 0.0
    wd_12  = 0.0
    wd_13  = 0.0
    f = open(output)
    SAPT2 = []
    for line in f:

        if SLVL_TAG in line:
           for line in f:
               line_adjl = line.lstrip("! ")
               if line_adjl.startswith(ELST_TAG):
                   e_elst = float(line_adjl.split("=")[1])
                   SAPT2.append([ELST_TAG,e_elst])
               if line_adjl.startswith(E1EX_TAG):
                   e_exch = float(line_adjl.split("=")[1])
                   SAPT2.append([E1EX_TAG,e_exch])
               if line_adjl.startswith(IND_TAG):
                   e_ind  = float(line_adjl.split("=")[1])
                   SAPT2.append([IND_TAG,e_ind])
               if line_adjl.startswith(EXI_TAG):
                   e_exi  = float(line_adjl.split("=")[1])
                   SAPT2.append([EXI_TAG,e_exi])
               if line_adjl.startswith(DSP_TAG):
                   e_disp = float(line_adjl.split("=")[1])
                   SAPT2.append([DSP_TAG,e_disp])
               if line_adjl.startswith(DCS_TAG):
                   e_disp = float(line_adjl.split("=")[1])
                   SAPT2.append([DCS_TAG,e_disp])
               if line_adjl.startswith(DRV_TAG):
                   e_disp = float(line_adjl.split("=")[1])
                   SAPT2.append([DRV_TAG,e_disp])
               if line_adjl.startswith(EXD_TAG):
                   e_exd  = float(line_adjl.split("=")[1])
                   SAPT2.append([EXD_TAG,e_exd])
               if line_adjl.startswith(XDRV_TAG):
                   e_exd  = float(line_adjl.split("=")[1])
                   SAPT2.append([XDRV_TAG,e_exd])
               if line_adjl.startswith(EINT_TAG):
                   e_int  = float(line_adjl.split("=")[1])
                   SAPT2.append([EINT_TAG,e_int])
               if line_adjl.startswith(WD1_TAG):
                   wd_12  = float(line_adjl.split("=")[1])
                   SAPT2.append([WD1_TAG,wd_12])
               if line_adjl.startswith(WD2_TAG):
                   wd_13  = float(line_adjl.split("=")[1])
                   SAPT2.append([WD2_TAG,wd_13])

    return SAPT2
#
def get_sapt0_energies(log_path,log_name):
    output = os.path.join(log_path,log_name)
    e_elst = 0.0
    e_exch = 0.0
    e_ind  = 0.0
    e_exi  = 0.0
    e_disp = 0.0
    e_exd  = 0.0
    f = open(output)
    SAPT0 = []
    for line in f:

        if SLVL_TAG in line:
           for line in f:
            line_adjl = line.lstrip("! ")
            if line_adjl.startswith(ELST_TAG):
                e_elst = float(line_adjl.split("=")[1])
                SAPT0.append([ELST_TAG,e_elst])
            if line_adjl.startswith(E1EX_TAG):
                e_exch = float(line_adjl.split("=")[1])
                SAPT0.append([E1EX_TAG,e_exch])
            if line_adjl.startswith(IND0_TAG):
                e_ind  = float(line_adjl.split("=")[1])
                SAPT0.append([IND0_TAG,e_ind])
            if line_adjl.startswith(EXI0_TAG):
                e_exi  = float(line_adjl.split("=")[1])
                SAPT0.append([EXI0_TAG,e_exi])
            if line_adjl.startswith(DSP0_TAG):
                e_disp  = float(line_adjl.split("=")[1])
                SAPT0.append([DSP0_TAG,e_disp])
            if line_adjl.startswith(EXD0_TAG):
                e_exd  = float(line_adjl.split("=")[1])
                SAPT0.append([EXD0_TAG,e_exd])

    f.close()
    e_int = e_elst + e_exch + e_ind + e_exi + e_disp + e_exd
    SAPT0.append([INT0_TAG,e_int])

    return SAPT0
#
def ac0_en(file):
   f=open(file,'r')
   s=f.read()
   f.close()
   a=s.find(" ECASSCF+ENuc, AC0-Corr, AC0-CASSCF")
   S=s[a+70:a+83]
   return float(S)
#
def ac_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find("ECASSCF+ENuc, AC-Corr, AC-ERPA-CASSCF")
        S=s[a+75:a+89]
        return float(S)
#
def ac1_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find(" ECASSCF+ENuc, AC1-Corr, ERPA-CASSCF")
        S=s[a+74:a+88]
        return float(S)
#
def ac0d_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find(" Dexcitation correction for AC0 ")
        S=s[a+60:a+71]
        return float(S)
#
def ac0gvb_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find("  EGVB+ENuc, 0th+1st-order ECorr, AC0-GVB")
        S=s[a+78:a+91]
        return float(S)
#
def acgvb_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find("  EGVB+ENuc, Corr, AC-ERPA-GVB")
        S=s[a+66:a+80]
        return float(S)
#
def ac1gvb_en(file):
   f=open(file,'r')
   for line in f:
      if line[0:26] == " EGVB+ENuc, Corr, ERPA-GVB" :
         l=line.split(" ")
         le=len(l)-1
         l=l[le]
         S=l[0:len(l)-1]
   return float(S)
#
def eerpagvb_en(file):
        f=open(file,'r')
        for line in f:
                if line[0:24] == "  EGVB + ENuc + 1,2-body" :
                        l=line.split(" ")
                        le=len(l)-1
                        l=l[le]
                        S=l[0:len(l)-1]
        return float(S)
#
def caspidft_en(file):
        f=open(file,'r')
        s=f.read()
        f.close()
        a=s.find(" PiDFT Correlation")
        S=s[a+30:a+43]
        return float(S)
#
#
def run_single(dir,fun):
# run calc in dir and extract data with fun
   cwd=os.getcwd()
   wrkdir=cwd+"/"+dir
   os.chdir(wrkdir)
  # os.system(cwd+"/gammcor > "+cwd+"/test.out")
   cmd = "{execpath} > {output}".format(execpath=RUN_GAMMCOR, output=cwd+"/test.out")
   subprocess.call(cmd, shell=True)
   os.chdir(cwd)
   t=fun("test.out")
   os.system("rm test.out")
   return t
#
def reference(dirl,fun):
        ref_val=[]
        for dir in dirl:
                ref_val.append(fun(dir+"/gammcor.out"))
        return ref_val
#
def run_test(dirl,fun):
	ref_val=reference(dirl,fun)
	print('ref :',ref_val)
#
	test_val=[]
	for dir in dirl:
		test_val.append(run_single(dir,fun))
#	
	print('test:',test_val)
	try:
		assert_allclose(test_val,ref_val,rtol=0,atol=tolene)
		print('passed with tol=',tolene,'[Ha] \n')
	except:
		print(' failed with errors:')
		print(np.subtract(test_val,ref_val))
	return 
#
def run_test_sapt(dirl,slevel):
#
   tolsapt=1.e-5
   ref_val =[]
   test_val=[]
   cwd=os.getcwd()
   for dir in dirl:
      wrkdir=cwd+"/"+dir
      os.chdir(wrkdir)
      cmd = "{execpath} > {output}".format(execpath=RUN_GAMMCOR, output=wrkdir+"/test.out")
      subprocess.call(cmd, shell=True)
      os.chdir(cwd)
      if slevel=='sapt2' :
         ref_tab = get_sapt2_energies(dir,"gammcor.out")
         tst_tab = get_sapt2_energies(dir,"test.out")
      elif slevel=='sapt0' :
         ref_tab = get_sapt0_energies(dir,"gammcor.out")
         tst_tab = get_sapt0_energies(dir,"test.out")


      ref_lbl = [ref_tab[i][0] for i in range(len(ref_tab))]
      ref_val = [ref_tab[i][1] for i in range(len(ref_tab))]
      tst_lbl = [tst_tab[i][0] for i in range(len(tst_tab))]
      tst_val = [tst_tab[i][1] for i in range(len(tst_tab))]

      jname = get_jobname(dir,"input.inp")
      print("Job name: ",jname)
      print("{ECONT:>16} {CALC:>12} {REF:>12} {DIFF:>12} {IPASS}".format(ECONT="Energy",
            CALC="ref_val",REF="test_val",DIFF="ref-tst",IPASS="Check "))
      resid = np.zeros(len(ref_val))
      for i in range(len(ref_val)):
         resid[i] = abs(tst_val[i] - ref_val[i])
         print("{ECONT}  {CALC:>12}  {REF:>12}  {DIFF}  {IPASS}".format(ECONT="%16s"%ref_lbl[i], 
               CALC="%11.8f"%tst_val[i], REF="%11.8f"%ref_val[i], DIFF="%10.8f"%resid[i],IPASS=resid[i]<=tolsapt))
      try:
         assert_allclose(tst_val,ref_val,rtol=0,atol=tolsapt)
         print('passed with tol=',tolsapt,'[mHa] \n')
      except:
         print(' failed with errors:')
         print(np.subtract(tst_val,ref_val))

   return
#
###### BEGIN TESTS ######

subprocess.call("export MKL_NUM_THREADS=1", shell=True)
#
## SAPT ##
dir_sapt2=["TESTS/SAPT/TEST1","TESTS/SAPT/TEST2","TESTS/SAPT/TEST3","TESTS/SAPT/TEST4","TESTS/SAPT/TEST5","TESTS/SAPT/TEST6","TESTS/SAPT/TEST7","TESTS/SAPT/GVB"]
dir_sapt0=["TESTS/SAPT/SAPT0"]
#
tolene=1.e-7
## AC0 ##
dir_ac0=["TESTS/AC0/TEST1","TESTS/AC0/TEST1/INCORE_INTEG" \
        ]
#         ,"TESTS/AC0/TEST2", \
#         "TESTS/AC0/TEST2/INCORE_INTEG"]
### AC ##
dir_ac=["TESTS/AC/TEST1","TESTS/AC/TEST1/INCORE_INTEG"]
### AC1 ##
dir_ac1=["TESTS/AC1/TEST1","TESTS/AC1/TEST1/INCORE_INTEG"]
### AC0D ##
dir_ac0d=["TESTS/AC0D/TEST2"]
### CASPIDFT ##
dir_caspidft=["TESTS/CASPIDFT/TEST1"]
# "TESTS/CASPIDFT/TEST2"] - this calculation is numerically unstable
### AC0-GVB ##
dir_ac0gvb=["TESTS/AC0GVB/TEST1","TESTS/AC0GVB/TEST1/INCORE_INTEG"]
### AC-GVB ##
dir_acgvb=["TESTS/ACGVB/TEST1","TESTS/ACGVB/TEST1/INCORE_INTEG"]
### AC1-GVB ##
dir_ac1gvb=["TESTS/AC1GVB/TEST1","TESTS/AC1GVB/TEST1/INCORE_INTEG"]
## EERPA-GVB ##
dir_eerpagvb=["TESTS/EERPAGVB/TEST1","TESTS/EERPAGVB/TEST2","TESTS/EERPAGVB/TEST3","TESTS/EERPAGVB/TEST4"]
### AC0_DALTON WITH CASSCF ##
dir_ac0dalton=["TESTS/AC0_DALTON/TEST1","TESTS/AC0_DALTON/TEST1/INCORE","TESTS/AC0_DALTON/TEST2","TESTS/AC0_DALTON/TEST2/INCORE"]
### AC0_DALTON WITH HF WAVEFUNCTION ##
dir_ac0_hf_dalton=["TESTS/AC0_HF_DALTON"]
#
# run tests
#
print("* testing SAPT energy components \n")
run_test_sapt(dir_sapt2,'sapt2')
#run_test_sapt(dir_sapt0,'sapt0')
#
print("* testing EERPA-GVB energy")
run_test(dir_eerpagvb,eerpagvb_en)
#
print ("* testing AC0-GVB energy")
run_test(dir_ac0gvb,ac0gvb_en)
#
print ("* testing AC-GVB energy")
run_test(dir_acgvb,acgvb_en)
##
print ("* testing AC1-GVB energy")
run_test(dir_ac1gvb,ac1gvb_en)
##
print("* testing AC0-CAS[molpro] energy")
run_test(dir_ac0,ac0_en)
##
print ("* testing AC0-CAS[dalton] energy")
run_test(dir_ac0dalton,ac0_en)
##
print ("* testing AC0-HF[dalton] energy")
run_test(dir_ac0_hf_dalton,ac0_en)
##
print ("* testing AC-CAS[molpro] energy")
run_test(dir_ac,ac_en)
##
print ("* testing AC1-CAS[molpro] energy")
run_test(dir_ac1,ac1_en)
##
print ("* testing PIDFT energy")
run_test(dir_caspidft,caspidft_en)
##
print ("* testing AC0D-CAS[molpro]")
run_test(dir_ac0d,ac0d_en)


