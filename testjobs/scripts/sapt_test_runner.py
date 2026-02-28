import os
import numpy as np
import subprocess
from numpy.testing import assert_allclose
from numpy.testing import assert_approx_equal
from scripts.colors import *

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

    return job_name.rstrip()
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
def run_test_sapt(units,slevel,RUN_GAMMCOR):
    no_passed = no_failed = 0
    failed_tests = []

    ref_val =[]
    test_val=[]
    cwd=os.getcwd()
    for unit in units:
        dir = unit.get('name')
        tolsapt = unit.get('atol')
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
        print(f"{BOLD}Job name: {jname}{ENDC}")
        print("{ECONT:>16} {CALC:>12} {REF:>12} {DIFF:>12} {IPASS}".format(ECONT="Energy",
            CALC="ref_val",REF="test_val",DIFF="ref-tst",IPASS="Check "))
        resid = np.zeros(len(ref_val))
        for i in range(len(ref_val)):
            resid[i] = abs(tst_val[i] - ref_val[i])
            print("{ECONT}  {CALC:>12}  {REF:>12}  {DIFF}  {IPASS}".format(ECONT="%16s"%ref_lbl[i], 
                CALC="%11.8f"%tst_val[i], REF="%11.8f"%ref_val[i], DIFF="%10.8f"%resid[i],IPASS=resid[i]<=tolsapt))
        
        try:
            assert_allclose(tst_val,ref_val,rtol=0,atol=tolsapt)
            print(f"{OKGREEN}passed with tol=',tolsapt,'[mHa] \n{ENDC}")
            no_passed += 1
        except:
            print(f"{FAIL} failed with errors:")
            print(f"{np.subtract(tst_val,ref_val)}{ENDC}")
            no_failed += 1
            failed_tests.append(f"SAPT: {jname}")
    
    return {'passed': no_passed, 'failed': no_failed, 'failed_tests': failed_tests}