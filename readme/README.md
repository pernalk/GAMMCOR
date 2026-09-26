# GAMMCOR

GAMMCOR computes electron correlation energies and intermolecular interaction
energies for multireference wave functions. Starting from a reference such as
CASSCF, GVB or DMRG (computed by an external program), it adds the missing
correlation with methods based on the adiabatic connection (AC, AC0) and
extended random-phase approximation (ERPA), and computes interaction energies
with symmetry-adapted perturbation theory (SAPT).

GAMMCOR does not compute the reference wave function itself. It reads it from
one of the supported quantum-chemistry programs:

* **PySCF** – through the `gammcor_tools` Python package, which writes a `pyscf_data.h5` file
* **Molpro**
* **Dalton**
* **ORCA**

<!-- TODO: Add a short description of the typical GAMMCOR workflow and cite
the relevant papers once the references have been confirmed. -->

# Documentation

* [Theory Primer](00_theory_primer.md)
* [Setup](01_setup.md)
* [Input](02_input_structure.md)
* [Examples](03_examples.md)
* [Literature](04_literature.md)

Full keyword reference: [GAMMCOR user manual](https://qchem.gitlab.io/gammcor-manual/).

# Authors

* Katarzyna Pernal
* Michał Hapka
* Aleksandra Tucholska
* Marcin Modrzejewski (`gammcor-integrals` submodule)

<!-- TODO: names, ORCID links and what each person contributed, as in:
* [Name Surname](https://orcid.org/...) (what they did) -->

<!-- TODO: Add the Citation section once the primary GAMMCOR references have
been confirmed. -->

# License

This program is distributed under the GNU General Public License, version 3.
For full licensing details, see the [LICENSE](../LICENSE) file.
