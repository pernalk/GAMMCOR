#! /usr/bin/env python
import os,sys,string,math,numpy
from numpy.testing import assert_allclose
from numpy.testing import assert_approx_equal

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
	wrkdir=dir
        os.chdir(wrkdir)
	os.system("prdmft.exe > "+cwd+"/test.out")
	os.chdir(cwd)
	t=fun("test.out")
	os.system("rm test.out")
	return t
#
def reference(dirl,fun):
        ref_val=[]
        for dir in dirl:
		dirn="/home/kasia/GAMMCOR/"+dir
                ref_val.append(fun(dirn+"/gammcor.out"))
        return ref_val
#
def run_test(dirl,fun):
	ref_val=reference(dirl,fun)
	print 'ref :',ref_val
#
	test_val=[]
	for dir in dirl:
		dirn="/home/kasia/GAMMCOR/"+dir
		test_val.append(run_single(dirn,fun))
#	
	print 'test:',test_val
	try:
		assert_allclose(test_val,ref_val,rtol=0,atol=tolene)
		print 'passed with tol=',tolene,'\n'
	except:
		print test_name,' failed with errors:'
		print numpy.subtract(test_val,ref_val)	
	return 

###### BEGIN TESTS ######
tolene=1.e-7
## AC0 ##
dir_ac0=["TESTS/AC0/TEST1","TESTS/AC0/TEST1/INCORE_INTEG" \
        ]
#         ,"TESTS/AC0/TEST2", \
#         "TESTS/AC0/TEST2/INCORE_INTEG"]
## AC ##
dir_ac=["TESTS/AC/TEST1","TESTS/AC/TEST1/INCORE_INTEG"]
## AC1 ##
dir_ac1=["TESTS/AC1/TEST1","TESTS/AC1/TEST1/INCORE_INTEG"]
## AC0D ##
dir_ac0d=["TESTS/AC0D/TEST2"]
## CASPIDFT ##
dir_caspidft=["TESTS/CASPIDFT/TEST1","TESTS/CASPIDFT/TEST2"]
## AC0-GVB ##
dir_ac0gvb=["TESTS/AC0GVB/TEST1","TESTS/AC0GVB/TEST1/INCORE_INTEG"]
## AC-GVB ##
dir_acgvb=["TESTS/ACGVB/TEST1","TESTS/ACGVB/TEST1/INCORE_INTEG"]
## AC1-GVB ##
dir_ac1gvb=["TESTS/AC1GVB/TEST1","TESTS/AC1GVB/TEST1/INCORE_INTEG"]
## EERPA-GVB
dir_eerpagvb=["TESTS/EERPAGVB/TEST1","TESTS/EERPAGVB/TEST2","TESTS/EERPAGVB/TEST3","TESTS/EERPAGVB/TEST4"]
## AC0_DALTON
dir_ac0dalton=["TESTS/AC0_DALTON/TEST1","TESTS/AC0_DALTON/TEST1/INCORE","TESTS/AC0_DALTON/TEST2","TESTS/AC0_DALTON/TEST2/INCORE"]
#
# run tests
#
print "* testing EERPA-GVB energy"
run_test(dir_eerpagvb,eerpagvb_en)
#
print "* testing AC0-GVB energy"
run_test(dir_ac0gvb,ac0gvb_en)
#
print "* testing AC-GVB energy"
run_test(dir_acgvb,acgvb_en)
#
print "* testing AC1-GVB energy"
run_test(dir_ac1gvb,ac1gvb_en)
#
print "* testing AC0-CAS[molpro] energy"
run_test(dir_ac0,ac0_en)
#
print "* testing AC0-CAS[dalton] energy"
run_test(dir_ac0dalton,ac0_en)
#
print "* testing AC-CAS[molpro] energy"
run_test(dir_ac,ac_en)
#
print "* testing AC1-CAS[molpro] energy"
run_test(dir_ac1,ac1_en)
#
print "* testing PIDFT energy"
run_test(dir_caspidft,caspidft_en)
#
print "* testing AC0D-CAS[molpro]"
run_test(dir_ac0d,ac0d_en)


