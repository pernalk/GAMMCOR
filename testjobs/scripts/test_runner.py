import os
import numpy as np
import subprocess
from numpy.testing import assert_allclose
from numpy.testing import assert_approx_equal
from scripts.colors import *

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
def ac1n_en(file):
    f=open(file,'r')
    s=f.read()
    f.close()
    a=s.find("  ECASSCF+ENuc, AC1-Corr, AC1-CASSCF")
    S=s[a+58:a+72]
    return float(S)
#
def acn_en(file):
    f=open(file,'r')
    s=f.read()
    f.close()
    a=s.find("  ECASSCF+ENuc, ACn-Corr, ACn-CASSCF")
    S=s[a+58:a+72]
    return float(S)
#

def run_single(dir,fun,RUN_GAMMCOR):
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
def run_test(units,fun,RUN_GAMMCOR):
    no_passed = no_failed = 0
    failed_tests = []
    
    dirl = []
    atol = []
    rtol = []
    for unit in units:
        dirl.append(unit.get('name'))
        atol.append(unit.get('atol'))

    ref_val=reference(dirl,fun)
    print('ref :',ref_val)
    #
    test_val=[]
    for dir in dirl:
        test_val.append(run_single(dir,fun,RUN_GAMMCOR))
    #	
    print('test:',test_val)
    for i in range(len(dirl)):
        try:
            assert_allclose(test_val[i],ref_val[i],rtol=0,atol=atol[i])
            print(f"{OKGREEN}{i+1}. passed with tol= {atol[i]} [Ha]{ENDC}")
            no_passed += 1
        except:
            print(f"{FAIL}{i+1}. failed with errors:")
            print(f"{np.subtract(test_val[i],ref_val[i])}{ENDC}")
            no_failed += 1
            failed_tests.append(f"{fun.__name__}: {dirl[i]}")
    
    print()
    return {'passed': no_passed, 'failed': no_failed, 'failed_tests': failed_tests}
