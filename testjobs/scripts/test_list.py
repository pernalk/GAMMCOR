from scripts.test_runner import *
from scripts.sapt_test_runner import *

def get():
    tests_list = []

    ## SAPT ###
    tests_list.append({
        'name': '* testing SAPT energy components',
        'fun': 'sapt2',
        'runner': run_test_sapt,
        'units': {
            "SAPT/Dalton/TEST7":   {'atol': 1.e-5, 'level': 'short'},
            "SAPT/QPackage/TEST1":   {'atol': 1.e-5, 'level': 'short'},
        }
    })

    ### ACn ###
    tests_list.append({
        'name': '* testing ACn-CAS[molpro] energy',
        'fun': acn_en,
        'runner': run_test,
        'units': {
            "ACn/AR_CAS":   {'atol': 1.e-5, 'level': 'short'},
            "ACn/AR_HF":    {'atol': 1.e-5, 'level': 'short'}
        }
    })

    return tests_list
