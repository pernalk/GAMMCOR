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

* [Setup](docs/01_setup.md)
* [Input](docs/02_input_structure.md)
* [PySCF interface](docs/03_pyscf_interface.md)
* [Examples](docs/04_examples.md)
* [Citation](docs/05_citation.md)

Full keyword reference: [GAMMCOR user manual](https://qchem.gitlab.io/gammcor-manual/).

# Citation

If you use GAMMCOR, please cite both the software and the following paper:

<table>
<thead>
<tr><th>Reference</th><th>BibTeX</th></tr>
</thead>
<tbody>
<tr>
<td><strong>GAMMCORSoftware</strong><br>K. Pernal, M. Hapka, M. Modrzejewski, A. Sokół, and A. Tucholska, <em>GAMMCOR</em>, computer software.<br><a href="https://github.com/pernalk/GAMMCOR">github.com/pernalk/GAMMCOR</a></td>
<td><details><summary>BibTeX</summary><pre><code>@misc{GAMMCORSoftware,
  author       = {Pernal, Katarzyna and Hapka, Micha{\l} and Modrzejewski, Marcin and Sok{\'o}{\l}, Adam and Tucholska, Aleksandra},
  title        = {{GAMMCOR}},
  howpublished = {Computer software},
  url          = {https://github.com/pernalk/GAMMCOR}
}</code></pre></details></td>
</tr>
<tr>
<td><strong>Hapka2023GAMMCOR</strong><br>M. Hapka, A. Krzemińska, M. Modrzejewski, M. Przybytek, and K. Pernal, “Efficient Calculation of the Dispersion Energy for Multireference Systems with Cholesky Decomposition: Application to Excited-State Interactions,” <em>J. Phys. Chem. Lett.</em> <strong>14</strong>, 6895–6903 (2023).<br><a href="https://doi.org/10.1021/acs.jpclett.3c01568">doi:10.1021/acs.jpclett.3c01568</a></td>
<td><details><summary>BibTeX</summary><pre><code>@article{Hapka2023GAMMCOR,
  author  = {Hapka, Micha{\l} and Krzemi{\'n}ska, Agnieszka and Modrzejewski, Marcin and Przybytek, Micha{\l} and Pernal, Katarzyna},
  title   = {Efficient Calculation of the Dispersion Energy for Multireference Systems with Cholesky Decomposition: Application to Excited-State Interactions},
  journal = {The Journal of Physical Chemistry Letters},
  year    = {2023},
  volume  = {14},
  number  = {30},
  pages   = {6895--6903},
  doi     = {10.1021/acs.jpclett.3c01568},
  url     = {https://doi.org/10.1021/acs.jpclett.3c01568}
}</code></pre></details></td>
</tr>
</tbody>
</table>

In addition, cite the publications relevant to the methods used in your
calculation. See [Citation](docs/05_citation.md) for the method-specific references
and ready-to-use BibTeX entries.

# Authors

* Katarzyna Pernal
* Michał Hapka
* Aleksandra Tucholska
* Marcin Modrzejewski

<!-- TODO: names, ORCID links and what each person contributed, as in:
* [Name Surname](https://orcid.org/...) (what they did) -->

# License

This program is distributed under the GNU General Public License, version 3.
For full licensing details, see the [LICENSE](LICENSE) file.
