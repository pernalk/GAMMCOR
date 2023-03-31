import os
import sys, getopt
import shutil
import re
import time
from copy import deepcopy
import numpy


# ------------------------------- PARAMS FOR SIROPT ------------------------------------------------------------------
GAMMCOR_TAG = "!$GAMMCOR --- GammCor interface, Pernal et al. 2023 ---\n"

"!write interface to GammCor program by Kasia Pernal\n"

SIROPT_STRING1 = """
C
C
      IF (GAMMCOR) THEN
      IF (NCONF .GT. 1) THEN
C     DUMP RHO2(i,j,k,l) for GammCor //Michal Hapka 2022
         KFREEG = 1
         LFREEG = LWRKG
         CALL DUMP_RHO2(WRK(KCREF),WRK(KCINDX),WRK(KWRKG),KFREEG,LFREEG)
      END IF
      END IF
"""

# ------------------------------- PARAMS FOR SRDFT ------------------------------------------------------------------

SRDFT_STRING1 = """
!     Should the grid information be dumped for GammCor program?
!     (DUMP_GAMMCOR is used to only do it once)

      IF (MYNUM .NE. 0) DUMP_GAMMCOR = .FALSE. ! only master should dump grid information
      IF (DUMP_GAMMCOR) THEN
         DUMP_GAMMCOR = GAMMCOR
      END IF
"""
SRDFT_STRING2 = """
      IF (DUMP_GAMMCOR) THEN ! dump grid information for GammCor
         CALL DUMP_DFTGRID_POINT(WGHT,COR,GSO)
      END IF
"""
SRDFT_STRING3 = """
      IF (DUMP_GAMMCOR) THEN
         CALL DUMP_DFTGRID_END()
         DUMP_GAMMCOR = .FALSE. ! only dump grid information once
      END IF
"""

# ------------------------------- END OF PARAMS  ------------------------------------------------------------------


################################################################################################################
#                                         START OF MAIN SCRIPT
################################################################################################################


### Set executable path
if len(sys.argv) < 2:
    raise Exception("Please provide Dalton path!")

### read user provided path
user_dalton_path = sys.argv[1]

### check if path correct
if not os.path.isdir(os.path.join(user_dalton_path, "DALTON")):
    print('')
    print("Wrong path: The DALTON directory does not exist in the path.")
    print('')
    sys.exit(1)


### set dalton path and script path
local_dalton_path = user_dalton_path.rstrip('/')
script_dir = os.path.dirname(os.path.abspath(__file__))

### set paths to the files that will be modified

gnrinf = os.path.join(local_dalton_path, "DALTON", "include", "gnrinf.h")
dalgnr = os.path.join(local_dalton_path, "DALTON", "main", "dalgnr.F")
gamm   = os.path.join(local_dalton_path, "DALTON", "sirius", "gammcor.F")
sirfck = os.path.join(local_dalton_path, "DALTON", "sirius", "sirfck.F")
siropt = os.path.join(local_dalton_path, "DALTON", "sirius", "siropt.F")
srdftjt= os.path.join(local_dalton_path, "DALTON", "srdft", "srdftjt.F")
dalton_source = os.path.join(local_dalton_path, "cmake", "SourcesDALTON.cmake")

### define border for a nice welcome
def print_with_border(output):
    border = "=" * 80
    double_border = "#" + border + "#"
    title = " GAMMCOR - DALTON - COMPATIBILITY "
    title_length = len(title) + 2
    tag = "!$GAMMCOR --- Gammcorr interface, Pernal et al. 2023 ---"
    tag_length = len(tag) + 2
    title_frame = "#" + "=" * ((84 - title_length) // 2) + title + "=" * ((80 - title_length) // 2 + (80 - title_length) % 2) + "#"
    tag_frame = "#" + " " * ((80 - tag_length) // 2) + tag + " " * ((80 - tag_length) // 2 + (80 - tag_length) % 2) + "#"
    print('')
    print(double_border)
    print(title_frame)
    print(double_border)
    print("|" + " " * 80 + "|")
    for line in output:
        print("| " + line.ljust(78) + " |")
    print("|" + " " * 80 + "|")
    print(tag_frame)
    print(double_border)
    
    

###  Define a list of files to be modified and added
files_to_modify = [gnrinf, dalgnr, sirfck, siropt, srdftjt, dalton_source]
files_to_add = [gamm]

### Define a description of the script
description1 = "This script modifies Dalton files to make them compatible with GammCorr."
description2 = "Remember to compile Dalton after running this script."

### Print the output with a border
output = [description1, description2, "", "FILES MODIFIED:"]
output.extend(files_to_modify)
output.append("")
output.append("FILES ADDED:")
output.extend(files_to_add)
print_with_border(output)


# ############################### COPY GAMMCOR FILE AND MODIFY DALTON_SOURCE ##########################################

gammcor_src = os.path.join(script_dir, "gammcor.F")
shutil.copy(gammcor_src, gamm)

with open(dalton_source, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if "DALTON/sirius/dypc.F" in line:
        if "!$GAMMCOR" not in lines[i+1]:
            lines.insert(i+1, f"    DALTON/sirius/gammcor.F #{GAMMCOR_TAG}")
            break
        else:
            break

with open(dalton_source, "w") as f:
    f.writelines(lines)


f.close()


# ############################### CHANGE OF THE GNRINF COMMON FILE ###################################################

with open(gnrinf, "r") as f:
    lines = f.readlines()
    
for i, line in enumerate(lines):
    if 'DKHINT' in line:
        if "!$GAMMCOR" not in lines[i+1]:
            lines.insert(i+1, f"//{GAMMCOR_TAG}")
            lines.insert(i+2, f"     &        GAMMCOR,                                                  &\n")
            lines.insert(i+3, f"//{GAMMCOR_TAG}")
            break
        else:
            break

with open(gnrinf, "w") as f:
    f.writelines(lines)

f.close()

# ############################### CHANGE OF THE DALGNR FILE ##########################################################

with open(dalgnr, "r") as f:
    lines = f.readlines()

last_label = 0
break_out_flag = False
for i, line in enumerate(lines):

    match = re.search(r'\bNTABLE\s*=\s*(\d+)', line)
    if match:
        ntable_value = int(match.group(1))
        new_ntable_value = str(ntable_value + 1)
        if "!$GAMMCOR" not in lines[i]:
            lines[i] = line.replace('NTABLE = ' + str(ntable_value), 'NTABLE = ' + new_ntable_value).rstrip('\n')
            lines[i] += f"{GAMMCOR_TAG}"

    if 'DATA TABLE' in line:
        for j in range(i+1, len(lines)):
            if '.FDE' in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f"     &            '.GAMMCO', {GAMMCOR_TAG}")
                    break_out_flag = True
                    break
                else:
                    break_out_flag = True
                    break
        if break_out_flag:
            break

break_out_flag = False
for i, line in enumerate(lines):
    if 'Initialize /GNRINF/' in line:
        for j in range(i+1, len(lines)):
            if 'WESTA' in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1,f"      GAMMCOR= .FALSE. {GAMMCOR_TAG}")
                    break_out_flag = True
                    break
                else:
                    break_out_flag = True
                    break
        if break_out_flag:
            break

break_out_flag = False
break_out_flag2 = False
for i, line in enumerate(lines):
    if " 100 READ (LUCMD" in line:
        for j in range(i+1, len(lines)):
            if "GO TO (" in lines[j]:
                for k in range(j+1, len(lines)):                    
                    if ")" in lines[k]:
                        if "!$GAMMCOR" not in lines[k]:
                            num_list = lines[k].split(',')
                            lsplit = lines[k].split(')')
                            for elem in num_list:
                                if ")" in elem:
                                    last_label = int(elem.rstrip(")"))
                                    new_label = last_label + 1
                                    new_line = lsplit[0] + f",{new_label})" + lsplit[1]
                                    lines[k] = lsplit[0] + ","
                                    lines.insert(k+1, f"\n     &                {new_label}) {GAMMCOR_TAG}")
                                    lines.insert(k+2, f"     &                  "+ lsplit[1])
                                    brake_out_flag = True
                                    break
                        else:
                            num_list = lines[k].split(',')
                            lsplit = lines[k].split(')')
                            for elem in num_list:
                                if ")" in elem:                                    
                                    match = re.search(r'\b(\d+)\)\s*!\$GAMMCOR', lines[k])
                                    if match:
                                        new_label = match.group(1)
                                        last_label = str(int(new_label) -1)
                            brake_out_flag = True
                            break
                if brake_out_flag:
                    brake_out_flag2 = True
                    break
        if brake_out_flag2:
            break

brake_out_flag = False
for i, line in enumerate(lines):
    if f"{last_label}     CONTINUE" in line:
        for j in range(i+1, len(lines)):
            if "GOTO 100" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f" {GAMMCOR_TAG}")
                    lines.insert(j+2,  f" {new_label}     CONTINUE  \n")
                    lines.insert(j+3, f"            GAMMCOR = .TRUE. \n")
                    lines.insert(j+4, f"         GOTO 100 ")
                    lines.insert(j+5, f"\n {GAMMCOR_TAG}")
                    brake_out_flag = True
                    break
                else:
                    brake_out_flag = True
                    break
        if brake_out_flag:
            break

for i, line in enumerate(lines):
    if "Information for WESTA" in line:
        if "!$GAMMCOR" not in lines[i+1]:
            lines.insert(i+1, f"{GAMMCOR_TAG}")
            lines.insert(i+2, f"      IF (GAMMCOR) WRITE (LUPRI,\'(4X,A)\') 'Information'//\n")
            lines.insert(i+3, f"     &  'for GammCor will be calculated and written to files.'")
            lines.insert(i+4, f"\n{GAMMCOR_TAG}")
            break
        else:
            break


with open(dalgnr, "w") as f:
    f.writelines(lines)

f.close()


# ############################### CHANGE OF THE SIRFCK FILE ##########################################################

with open(sirfck, "r") as f:
    lines = f.readlines()


for i, line in enumerate(lines):
    if "GNRINF: SRINTS, CHIVAL" in line:
        if "!$GAMMCOR" not in line:
            lines[i] = f"C   GNRINF: SRINTS, CHIVAL, GAMMCOR {GAMMCOR_TAG}"

    if "EJSR = EJSR - EJKVSR" in line:
        if "!$GAMMCOR" not in lines[i+2]:
            lines.insert(i+2, f"{GAMMCOR_TAG}")
            lines.insert(i+3, "\n            IF (GAMMCOR) THEN ! dump Jsr + HFXFAC*Ksr to GammCor program\n")
            lines.insert(i+4, f"               CALL DUMP_JSR(WORK(KFSRAO))\n")
            lines.insert(i+5, f"            END IF\n")
            lines.insert(i+6, f"{GAMMCOR_TAG}")

    break_out_flag = False
    if "Add Short-range integrals to FC" in line:
        for j in range(i+1, len(lines)):
            if "EJKVSR = D0" in lines[j]:
                for k in range(j+1, len(lines)):
                    if "!" in lines[k]:
                        if "!$GAMMCOR" not in lines[k]:
                            lines.insert(k, f"{GAMMCOR_TAG}")
                            lines.insert(k+1, f"            IF (DOERG .AND. GAMMCOR) THEN ! dump Jsr + HFXFAC*Ksr to GammCor program \n")
                            lines.insert(k+2, f"               CALL DUMP_JSR(WORK(KFSRAO)) \n")
                            lines.insert(k+3, f"            END IF \n")
                            lines.insert(k+4, f"{GAMMCOR_TAG}")
                            break_out_flag = True
                            break
                        else:
                            break_out_flag = True
                            break                    
                if break_out_flag:
                    break
    
with open(sirfck, "w") as f:
    f.writelines(lines)

f.close()

# ############################### CHANGE OF THE SIROPT FILE ##########################################################

with open(siropt, "r") as f:
    lines = f.readlines()
        
for i, line in enumerate(lines):

    break_out_flag = False
    if "only EKT output" in line:
        for j in range(i+1, len(lines)):
            if "END IF" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f'{GAMMCOR_TAG}')
                    lines.insert(j+2, SIROPT_STRING1)
                    lines.insert(j+3, f'{GAMMCOR_TAG}')
                    break_out_flag = True
                    break
                else:
                    break_out_flag = True
                    break  
        if break_out_flag:
            break


with open(siropt, "w") as f:
    f.writelines(lines)

f.close()

# ############################### CHANGE OF THE SRDFTJT FILE ##########################################################


with open(srdftjt, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):

    if "logical, intent(in) :: DOLND," in line:
        for j in range(i+1, len(lines)):
            if "logical :: FROMVX" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f"{GAMMCOR_TAG}")
                    lines.insert(j+2, f"      logical, save :: DUMP_GAMMCOR = .TRUE. \n")
                    lines.insert(j+3, f"{GAMMCOR_TAG}")
                    break
                else:
                    break


    if "SRDFT1JT: triplet only for response" in line:
        for j in range(i+1, len(lines)):
            if "END IF" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f"{GAMMCOR_TAG}")
                    lines.insert(j+2, SRDFT_STRING1)
                    lines.insert(j+3, f"{GAMMCOR_TAG}")
                    break
                else:
                    break

    if "CALL GPOPEN" in line:
        for j in range(i+1, len(lines)):
            if ")" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f" {GAMMCOR_TAG}")
                    lines.insert(j+2, f"      IF (DUMP_GAMMCOR) CALL DUMP_DFTGRID_INIT(NTOT_DFTGRID,DOGGA) \n")
                    lines.insert(j+3, f"{GAMMCOR_TAG}")
                    break
                else:
                    break

    if "CALL GETSOS" in line:
        for j in range(i+1, len(lines)):
            if ")" in lines[j]:
                if "!$GAMMCOR" not in lines[j+1]:
                    lines.insert(j+1, f"{GAMMCOR_TAG}")
                    lines.insert(j+2, SRDFT_STRING2)
                    lines.insert(j+3, f"{GAMMCOR_TAG}")
                    break
                else:
                    break

    if "CALL GPCLOSE" in line:
        for j in range(i+1, len(lines)):
            if "RETURN" in lines[j]:
                if "!$GAMMCOR" not in lines[j-1]:
                    lines.insert(j, f"\n{GAMMCOR_TAG}")
                    lines.insert(j, SRDFT_STRING3)
                    lines.insert(j, f"{GAMMCOR_TAG}")
                    break
                else:
                    break


with open(srdftjt, "w") as f:
    f.writelines(lines)

f.close()

print('')
print('#################################### SUCCESS #####################################')
print('')

