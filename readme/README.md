# GAMMCOR

GAMMCOR computes electron correlation energies and intermolecular interaction
energies for multireference wave functions. Starting from a reference such as
CASSCF, GVB or DMRG (computed by an external program), it adds the missing
correlation with methods based on the adiabatic connection (AC, AC0) and
extended random-phase approximation (ERPA), and computes interaction energies
with symmetry-adapted perturbation theory (SAPT).

GAMMCOR does not compute the reference wave function itself. It reads it from
one of the supported quantum-chemistry programs:

* **PySCF**
* **Molpro**
* **Dalton**
* **ORCA**

<!-- TODO: Add a short description of the typical GAMMCOR workflow. -->

# Documentation

* [Theory Primer](00_theory_primer.md)
* [Setup](01_setup.md)
* [Input](02_input_structure.md)
* [PySCF interface](03_pyscf_interface.md)
* [Examples](04_examples.md)
* [Citation](05_citation.md)

Full keyword reference: [GAMMCOR user manual](https://qchem.gitlab.io/gammcor-manual/).

# Citation

If you use GAMMCOR, please cite both the software and the following paper:

1. K. Pernal, M. Hapka, M. Modrzejewski, A. Sokół, and A. Tucholska,
   [GAMMCOR](https://github.com/pernalk/GAMMCOR), computer software.
2. M. Hapka, A. Krzemińska, M. Modrzejewski, M. Przybytek, and K. Pernal,
   “Efficient Calculation of the Dispersion Energy for Multireference Systems
   with Cholesky Decomposition: Application to Excited-State Interactions,”
   *J. Phys. Chem. Lett.* **14**, 6895–6903 (2023),
   [doi:10.1021/acs.jpclett.3c01568](https://doi.org/10.1021/acs.jpclett.3c01568).

In addition, cite the publications relevant to the methods used in your
calculation. See [Citation](05_citation.md) for the method-specific references
and ready-to-use BibTeX entries.

# Authors

* Katarzyna Pernal
* Michał Hapka
* Aleksandra Tucholska
* Marcin Modrzejewski (`gammcor-integrals` submodule)

<!-- TODO: names, ORCID links and what each person contributed, as in:
* [Name Surname](https://orcid.org/...) (what they did) -->

# License

This program is distributed under the GNU General Public License, version 3.
For full licensing details, see the [LICENSE](../LICENSE) file.
