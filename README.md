# GAMMCOR

### Compilation
Unpack the gammcor.tar.gz package and inside the directory run
`make`

To compile with gfortran change flags in Makefie.

### Electron integrals
The code:
- is interfaced either  with Dalton (CASSCF, GVB) or with Molpro (CASSCF),
- requires: 1- and 2-electron integrals in the AO representation
  in Dalton we use the AOONEINT and AOTWOINT files available, e.g, through the "-get" option,
  in Molpro damping of the integrals is done through the "user" interface
           (based on the internal "printints1-3" Molpro procedures).

All AO integrals are sorted in our code and saved on disk (currently we keep all integrals, i.e.
numerical zero values are not discarded).
Two-electron AO integrals are stored in a (triang | triang) convention, where
triang =  nbas*(nbas+1)/2 and (11 | 22) "chemists'  notation is used,
i.e.,  the AO file contains a "triang" number of "rs"  records and
each "rs" records contains a "pq" triangle of AO integrals.

  Both AC0 and SAPT calculations require integrals in the NO representation.
We transform and store two types of NO 2-electron integrals on disk:  (FF | OO) and (FO | FO),
where the (11 | 22) chemists' notation is used ,
"F" denotes all NO-indices and
"O" refers to the Occupied (=inactive + active) NO indices.
Note that SAPT calculations require the use of a general 4-index transformation, in which each index
can be transformed by a different NO-coefficient matrix (belonging to different monomers).

## EERPA
* JobType
  * EERPA     (alternatively: EERPA-2)
  * EERPA-OLD (alternatively: EERPA-1)

EERPA or EERPA-2 sets     IFl12 = 1. EERPA-OLD or EERPA-1 sets IFl12 = 2. Both keywords set IFlAC=0 so ERPA-GVB is carried our for fragments.

## AC,AC0,AC1
Calculation block
* Core
  * 0  - inactive orbitals not correlated
  * 1  - inactive orbitals contribute to correlation
* JobType
  * AC
  * AC0
  * AC1
 
By deafault integrals are kept in memory ("in-core" version). To use disc variant use:
   ``TwoMoInt   FOFO``

WARNING: AC0 correlation energies from in-core and fofo variants may vary if SA-CAS RDM's are used!
         This is due to the fact that 1st-order ABPLUS/MIN matrices are not symmetric. The in-core 
         version uses only the upper triangles of ABPLUS/MIN which are copied to lower triangle. 
         In the fofo version full AB matrices are computed and symmetrized.  

## AC0D : AC0 deexcitation corrections
* JobType
  * AC0D: corrections for states selected based on symmetry ISymmAC0D set to 1
  * AC0DNOSYMM: corrections for states selected base on overlap with SA-CAS states ISymmAC0D set to 0   

Output: 
* AC0D corrections
* for state 1.1 transition dipole moments in the (0)- and (1)-order ERPA approximations
* for AC0DNOSYMM for state 1.1 transition dipole moments in CAS approximation (based on transition matrices from molpro)

## CASPiDFT
* JobType: CASPIDFTOPT

PiDFT optimal parameters from _Hapka et al., J. Phys. Chem. Lett. 2020, 11, 5883_ are used
or they can be read from the command line in the order: A B G C D. 

Input: 1RDM, 2RDM, GRID, AO Orb on grid, CAS orbs.
Input files (generated from molpro): 2RDM, GRID, MOLPRO.MOPUN, 
AOONEINT.mol, AOTWOINT.mol (must be present). xcfun library is used.

Output: PiDFT and LYP Correlation Energy and the Ionicity Index.

## SAPT 

### GENERAL REMARKS

* Current SAPT implementation works only in the dimer-centered basis (DCBS).
* Extension to open-shell monomers does not work yet.
* In default settings AO integrals are sorted and transformed out-of-core. 
   Storage of integrals requires a considerable amount of disc space.
   Largest calculations to date were performed for 384 basis functions 
   (PCCP dimer in aug-cc-pVTZ) and took ca. 100 GB of disc space. 
   No SCRATCH directory option has been added so far, so either make sure 
   that you have enough space in the current working directory, or run jobs 
   straight from your scratch directory.
* The default settings do not make use of symmetry.
* Currently, only 1st and 2nd order Rayleigh-Schrodinger energy contributions
   are available (1st order ExchS2 is also implemented, but not made avaiable).

### SAPT: INPUT FILE
Several remarks on the input file:

* for now, the SAPT input file should have a fixed name "input.inp" and be 
   present in the current working directory, 
* definitions of three blocks are required: one Calculation block 
   and two System blocks,
* neither the order of blocks nor order of keywords within a block matters,
* each block terminates with the "end" keyword,
* keywords are not case sensitive,
* comments are possible and follow the "!" character (see the example below), 
* although the geometry is defined separately either in Molpro or Dalton files, 
   SAPT calculation requires defining both the number of atoms and nuclear charge 
   for each monomer (NAtoms and ZNucl keywords in the System block, respectively),
* default settings in the Calculation block are set to SAPT-CAS calculation
   via Dalton interface.

Available keywords:
* Calculation block
  * JobTitle:
                put anything, it will be printed in the output 
  * Interface:
                 Dalton
                 Molpro
  * JobType:    
               AC
               AC0
               AC1
               ERPA 
               SAPT
  * RDMType: 
               CAS (default)
               HF
               GVB 

      *note that GVB works only with Dalton interface,

  * Response: 
               ERPA (default)
               TD-APSG
               DFT
               
      *for SAPT-CAS and SAPT-GVB calculations ERPA is the default,
      **SAPT-DFT calculations are still under development and require
        both the Molpro interface and xcfun library.

  * DFunc:
            srLDA
            srPBE
            PBE 

      *functional choice is possible only with the "Response DFT" option,
      **srLDA: doi:10.1103/PhysRevB.73.155111
      **srPBE: no idea, adapted from Molpro
      ***PBE:  doi:10.1103/PhysRevLett.77.3865

  * SaptLevel:
                SAPT0
                SAPT2 (default)
      *SAPT0 calculates only uncoupled/single-pole and semicoupled E2disp          
               
* System block
  * Monomer:
              A or B
  * NAtoms (required)
  * ZNucl  (required)
  * Charge (default=0)
  * Multiplicity (default=0)
  * TwoMoInt:
               InCore (default for SAPT-GVB)
               Full
               FOFO (default for SAPT-CAS/HF)
      *chooses between in-core and out-of-core 2-el MO transformation
       FOFO creates (FF|OO) and (FO|FO) integrals and requires less disc
       space than (FF|FF)                 

  * ThrSelAct (default=1.d-8 for CAS
                 default=1.d-2 for GVB)
      *threshold for discarding pairs of degenerate or nearly
       degenerate orbitals which may cause instabilities in hessian
       matrices (will appear in the output as "Skipped" eigenvalues

### INPUT EXAMPLES

##### Input for SAPT-CAS calculation:

* input for SAPT-CAS with Dalton interface
``
  Calculation
  JobTitle   H2-H2
  Interface  DALTON
  JobType    SAPT
  RDMType    CAS
  IPrint     1  
  end

  ! define monA
  System
  Monomer A
  NAtoms  2
  ZNucl   2
  Charge  0
  Multiplicity 1
  end

  ! define monB
  System
  Monomer B
  NAtoms  2
  ZNucl   2
  Charge  0
  Multiplicity 1
  end
``

* input for SAPT-GVB calculation:
  ! H2O-H2O SAPT-GVB with Dalton
  Calculation
  JobTitle   H2O-H2O
  Interface  DALTON
  JobType    SAPT
  RDMType    GVB
  IPrint     1
  end

  ! define monA
  System
  Monomer A
  NAtoms  3
  ZNucl   10
  Charge  0
  Multiplicity 1
  end

  ! define monB
  System
  Monomer B
  NAtoms  3
  ZNucl   10
  Charge  0
  Multiplicity 1
  end

### RUNNING A SAPT-DALTON CALCULATION

To run a SAPT-CAS, SAPT-GVB or SAPT-HF calculation based on Dalton files:
1) Prepare *.mol and .dal files for each monomer (see examples). 
   a) unfortunately, in Dalton2015 ghost (dummy) atoms cannot be placed 
      before regular atoms, so the active monomer should always 
      be defined first,
   b) basis functions in *.mol files have to be defined explicitly (Dalton2015 
      does not offer an easy way to define ghost atoms); basis sets in Dalton
      format can be found in the /basis subdirectory of your Dalton directory.
2) Prepare "input.inp" file for SAPT calculation.
3) Run Dalton CAS or GVB calculations for each of the monomers and store 
   the following files on disc: 
   - for GVB calculation: AOTWOINT AOONEINT AMFI_SYMINFO.TXT SIRIFC SIRIUS.RST DALTON.MOPUN ,
   - for CAS calculations additionally store: rdm2.dat .

   An example Dalton run of job.mol and job.dal files:
   yourdaltondir/build/dalton -get "AOTWOINT AOONEINT AMFI_SYMINFO.TXT SIRIFC SIRIUS.RST DALTON.MOPUN " job  > skas
4) Run SAPT:
   yourprdmftdir/prdmft.exe > OUTPUT

The /EXAMPLES directory contains a bash script (run_Dalton_PRDMFT.sh) 
that handles both the Dalton and SAPT jobs.

To summarize, the Dalton+SAPT run requires 5 input files: 4 files for Dalton
calculations (two *.dal and two *.mol files for every monomer) 
and 1 input.inp file for SAPT-CAS/GVB calculation.

### RUNNING A SAPT-MOLPRO CALCULATION

To run a SAPT-CAS, SAPT-HF or SAPT-DFT calculation based on Molpro files:
1) Prepare input for Molpro. To dump the necessary AOs and RDMs use 
   the "user,prdmft-sapt,..." keywords (see examples).
2) Prepare "input.inp" file for SAPT calculation.  
3) Run Molpro calculation with "-d ." option to store scratch files:
   yourmolprodir/src/prdmft/molpro -d . run.inp 
4) Run SAPT:
   yourprdmftdir/prdmft.exe > OUTPUT

To summarize, the Molpro+SAPT run requires 2 input files: 
1 file for Molpro calculation and 1 "input.inp" file for SAPT-CAS/DFT calculaton.


