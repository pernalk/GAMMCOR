from scripts.test_runner import *
from scripts.sapt_test_runner import *

def get():
    tests_list = []

    ## SAPT ###
    #tests_list.append({
    #    'name': '* testing SAPT energy components',
    #    'fun': 'sapt2',
    #    'runner': run_test_sapt,
    #    'units': {
    #        "SAPT/Dalton/TEST6":   {'atol': 1.e-5, 'level': 'short'},
    #        "SAPT/Dalton/TEST7":   {'atol': 1.e-5, 'level': 'short'},
    #        "SAPT/QPackage/TEST1": {'atol': 1.e-5, 'level': 'short'},
    #    }
    #})

    ### AC0 ###
    tests_list.append({
        'name': '* testing AC0-CAS energy',
        'fun': ac0_en,
        'runner': run_test,
        'units': {
            "AC0/ethene":  {'atol': 1.e-7, 'level': 'short'},
        }
    })

    ### AC0, N2 CAS(6,6) state-specific (the same CASSCF in PySCF and Molpro) ###
    # compare: pairs that must agree with each other (fresh results, no extra runs).
    # FullA1 tests are in test-SA: for a state-specific reference A1 is symmetric.
    AC0 = "pyscf-tests/test-AC0/"
    tests_list.append({
        'name': '* testing AC0-CAS energy, state-specific N2 (THC, FOFO, INCORE; PySCF and Molpro)',
        'fun': ac0_en,
        'runner': run_test,
        'units': {
            AC0+"thc-alg":                         {'atol': 1.e-7, 'level': 'short', 'fun': thc_ac0_en},
            AC0+"thc-exact-alg":                   {'atol': 1.e-7, 'level': 'short', 'fun': thc_ac0_en},
            AC0+"fofo-alg-pyscf":                  {'atol': 1.e-7, 'level': 'short'},
            AC0+"fofo-alg-molpro":                 {'atol': 1.e-7, 'level': 'short'},
            AC0+"incore-alg-pyscf":                {'atol': 1.e-7, 'level': 'short'},
            AC0+"incore-alg-molpro":               {'atol': 1.e-7, 'level': 'short'},
            AC0+"test-ac0-incore-natorb":          {'atol': 1.e-7, 'level': 'short'},
            AC0+"test-ac0-incore-canonical":       {'atol': 1.e-7, 'level': 'short'},
            AC0+"test-ac0-incore-canonical-nosym": {'atol': 1.e-7, 'level': 'short'},
            AC0+"test-ac0-incore-rotated-active":  {'atol': 1.e-7, 'level': 'short'},
        },
        'compare': [
            (AC0+"thc-alg",                 AC0+"fofo-alg-pyscf",           1.e-8),  # 2e-10: same data
            (AC0+"thc-alg",                 AC0+"fofo-alg-molpro",          1.e-6),  # 9e-7: PySCF THC E_CAS from THC integrals
            (AC0+"thc-alg",                 AC0+"thc-exact-alg",            1.e-5),  # 8e-7: THC error
            (AC0+"incore-alg-pyscf",        AC0+"incore-alg-molpro",        1.e-7),  # 1e-8
        ],
    })

    ### AC0, LiF CAS(6,6) SA2, state 2.1 (the same SA-CASSCF in PySCF and Molpro) ###
    # A1 is not symmetric for SA: INCORE without FullA1 keeps one triangle, THC,
    # FOFO and INCORE FullA1 average it. So INCORE is compared only with INCORE,
    # and THC-exact with INCORE FullA1.
    SA = "pyscf-tests/test-SA/"
    tests_list.append({
        'name': '* testing AC0-CAS energy, SA LiF (THC, FOFO, INCORE, INCORE FullA1; PySCF and Molpro)',
        'fun': ac0_en,
        'runner': run_test,
        'units': {
            SA+"thc-alg":                   {'atol': 1.e-7, 'level': 'short', 'fun': thc_ac0_en},
            SA+"thc-exact-alg":             {'atol': 1.e-7, 'level': 'short', 'fun': thc_ac0_en},
            SA+"fofo-alg-pyscf":            {'atol': 1.e-7, 'level': 'short'},
            SA+"fofo-alg-molpro":           {'atol': 1.e-7, 'level': 'short'},
            SA+"incore-alg-pyscf":          {'atol': 1.e-7, 'level': 'short'},
            SA+"incore-alg-molpro":         {'atol': 1.e-7, 'level': 'short'},
            SA+"incore-alg-A1full-pyscf":   {'atol': 1.e-7, 'level': 'short'},
            SA+"incore-alg-A1full-molpro":  {'atol': 1.e-7, 'level': 'short'},
        },
        'compare': [
            (SA+"thc-alg",                 SA+"fofo-alg-pyscf",           1.e-8),  # 1e-10: same data
            (SA+"thc-alg",                 SA+"fofo-alg-molpro",          1.e-6),  # 7e-7: PySCF THC E_CAS from THC integrals
            (SA+"thc-alg",                 SA+"thc-exact-alg",            1.e-5),  # 8e-7: THC error
            (SA+"thc-exact-alg",           SA+"incore-alg-A1full-pyscf",  1.e-7),  # 4e-8: same data
            (SA+"thc-exact-alg",           SA+"incore-alg-A1full-molpro", 1.e-7),  # 8e-9
            (SA+"incore-alg-pyscf",        SA+"incore-alg-molpro",        1.e-6),  # 9e-8: PySCF vs Molpro SA-CASSCF (E_CAS 2e-7)
            (SA+"incore-alg-A1full-pyscf", SA+"incore-alg-A1full-molpro", 1.e-7),  # 3e-8
        ],
    })

    tests_list.append({
        'name': '* testing PySCF SS AC0-CAS energy (FOFO)',
        'fun': ac0_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/test-SS-fofo":  {'atol': 1.e-7, 'level': 'short'},
        }
    })

    ### pyscf AC0 SS, THC ###
    tests_list.append({
        'name': '* testing PySCF SS AC0-CAS energy (THC)',
        'fun': thc_ac0_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/test-SS-thc":  {'atol': 1.e-7, 'level': 'short'},
        }
    })


    ### MS-AC0, LiF CAS(6,6) SA2, state 1.1 (the same SA-CASSCF in PySCF and Molpro) ###
    # values: [E_total, H(1,1), H(2,1)]. INCORE/FOFO MSAC0 read the Molpro 2RDM
    # file only, so those units use the Molpro export.
    MS = "pyscf-tests/test-MSAC0/"
    tests_list.append({
        'name': '* testing MS-AC0 energy, H(1,1), H(2,1) (THC, FOFO, INCORE, INCORE FullA1)',
        'fun': msac0_molpro_en,
        'runner': run_test,
        'units': {
            # reference from 2026-10-13; the current code gives 1.8e-6 lower E_total
            MS+"thc-alg":            {'atol': 1.e-5, 'level': 'short', 'fun': msac0_en},
            MS+"thc-exact-alg":      {'atol': 1.e-7, 'level': 'short', 'fun': msac0_en},
            MS+"fofo-alg-molpro":    {'atol': 1.e-7, 'level': 'short'},
            MS+"incore-alg":         {'atol': 1.e-7, 'level': 'short'},
            MS+"incore-alg-A1full":  {'atol': 1.e-7, 'level': 'short'},
        },
        'compare': [
            (MS+"thc-alg",        MS+"fofo-alg-molpro",   1.e-5),  # 1.3e-6: THC error, PySCF vs Molpro
            (MS+"thc-alg",        MS+"thc-exact-alg",     1.e-5),  # 1.1e-6: THC error
            (MS+"thc-exact-alg",  MS+"incore-alg-A1full", 1.e-6),  # 2.2e-7: PySCF vs Molpro SA-CASSCF
        ],
    })

    ### pyscf lrAC0 ###
    tests_list.append({
        'name': '* testing PySCF CAS-srPBPE+lrAC0 energy',
        'fun': lrac0_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/test-sr-pbe":  {'atol': 1.e-5, 'level': 'short'},
        }
    })

    tests_list.append({
        'name': '* testing PySCF MP2 energy',
        'fun': mp2_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/CBSH/test-mp2":              {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-frozen":       {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cbsh":         {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-frozen-cbsh":  {'atol': 1.e-5, 'level': 'short'},
        }
    })

    tests_list.append({
        'name': '* testing PySCF ROMP2 energy',
        'fun': mp2_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/CBSH/test-romp2":             {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-romp2-frozen":      {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-romp2-cbsh":        {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-romp2-frozen-cbsh": {'atol': 1.e-5, 'level': 'short'},
        }
    })

    tests_list.append({
        'name': '* testing PySCF CAS MP2 energy (FOFO)',
        'fun': ac0_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/CBSH/test-mp2-cas-fofo":               {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-frozen":        {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-opshell":       {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-frozen-opshell": {'atol': 1.e-5, 'level': 'short'},
        }
    })

    tests_list.append({
        'name': '* testing PySCF CAS MP2 energy (THC)',
        'fun': thc_ac0_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/CBSH/test-mp2-cas-thc":               {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-thc-frozen":        {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-thc-opshell":       {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-thc-frozen-opshell": {'atol': 1.e-5, 'level': 'short'},
        }
    })

    tests_list.append({
        'name': '* testing PySCF CAS MP2 CBS[H] energy (FOFO)',
        'fun': ac0_cbsh_en,
        'runner': run_test,
        'units': {
            "pyscf-tests/CBSH/test-mp2-cas-fofo-cbsh":               {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-frozen-cbsh":        {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-cbsh-opshell":       {'atol': 1.e-5, 'level': 'short'},
            "pyscf-tests/CBSH/test-mp2-cas-fofo-frozen-cbsh-opshell": {'atol': 1.e-5, 'level': 'short'},
        }
    })



    ### ACn ###
    tests_list.append({
        'name': '* testing ACn-CAS[molpro] energy',
        'fun': acn_en,
        'runner': run_test,
        'units': {
            "ACn/AR_CAS":   {'atol': 1.e-5, 'level': 'short'},
            "ACn/AR_HF":    {'atol': 1.e-5, 'level': 'short'},
            "ACn/H2O_MOLMPS": {'atol': 1.e-5, 'level': 'short'}
        }
    })

    #### rdmcorr ###
    #tests_list.append({
    #    'name': '* testing AC0-1RDM[orca] energy',
    #    'fun': ac0_rdm,
    #    'runner': run_test,
    #    'units': {
    #        "RDMCORR/H2O_MOLMPS": {'atol': 1.e-5, 'level': 'short'}
    #    }
    #})


    return tests_list
