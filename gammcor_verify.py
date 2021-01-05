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
	os.system(cwd+"/gammcor > "+cwd+"/test.out")
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
	print 'ref :',ref_val
#
	test_val=[]
	for dir in dirl:
		test_val.append(run_single(dir,fun))
#	
	print 'test:',test_val
	try:
		assert_allclose(test_val,ref_val,rtol=0,atol=tolene)
		print 'passed with tol',tolene
	except:
		print 'test failed with errors:'
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
#
print "testing AC0-CAS energy"
run_test(dir_ac0,ac0_en)
#
print "testing AC-CAS energy"
run_test(dir_ac,ac_en)
#
print "testing AC1-CAS energy"
run_test(dir_ac1,ac1_en)
#
print "testing PIDFT energy"
run_test(dir_caspidft,caspidft_en)
#
print "testing AC0D"
run_test(dir_ac0d,ac0d_en)

