PySCF-GAMMCOR Interface Guide
This guide explains how to run calculations using the PySCF-GAMMCOR interface.

Quick Start: Required Basis Set Files
For the interface to work, you must provide two basis set files, both located in the path/to/gammcor/bazy-do-pyscf directory:

A .txt file for GAMMCOR (e.g., mybasis.txt).

A .nw file for PySCF (e.g., mybasis.nw).

Python Script Configuration (PySCF)
In your Python script, you need to configure the path and call the data generation function.

Add the GAMMCOR scripts to your Python path:

import sys
# Add the path to the directory containing Pyscf2Gammcor.py
sys.path.append('/path/to/gammcor/')
from Pyscf2Gammcor import get_data_for_gammcor

Set the basis_file variable:
Point this to the .nw basis set file.

basis_file = '/path/to/gammcor/bazy-do-pyscf/mybasis.nw'

# ... your PySCF mol, myhf, mycas setup ...

Generate GAMMCOR data:
Call this function at the very end of your script to save the calculation data.

# ... your PySCF calculations ...

get_data_for_gammcor(mol, myhf, mycas)

Run the script:
Execute your PySCF script as usual.

python3.12 your_script_name.py
or python3 whatever you have 


GAMMCOR Input Configuration
In your GAMMCOR input file, specify the interface and the path to your basis sets. The Basis keyword should point to the .txt file.

Interface PYSCF

BasisPath /path/to/gammcor/bazy-do-pyscf/
Basis mybasis.txt

# Recommended Cholesky settings
Choleskyblock
Cholesky OTF
Cholesky_Accuracy L
end

Optional: How to Add a New Basis Set if it’s not in bazy-do-pyscf directory
If you need to use a basis set that is not already in the directory, follow these steps.

Go to the Basis Set Exchange.

Find and select your desired basis set.

Download the basis set in two formats using the same base filename:

Select NWChem format and save it (e.g., mybasis.nw).

Select GAMESS (US) format and save it with a .txt extension (e.g., mybasis.txt).

Place both files into your basis set directory (e.g., /path/to/gammcor/bazy-do-pyscf/).