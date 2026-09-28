# RDMs in `Pyscf2Gammcor_mini.py` and `interface_pp.f90`

## 0. Convention names used in this document

| Name | Definition |
|---|---|
| **pyscf** | `0123` |
| **gammcor** | `0213`|
| **other** | -|

---

## 1. Arrays produced by `Pyscf2Gammcor_mini.py`

### 1a. `get_data_for_gammcor` — state-specific branch (`mycas.weights is None`)

| array (legacy file) | HDF5 dataset | convention | expression |
|---|---|---|---|
| `rdm2_aaaa.bin` | `.../RDM2_AAAA` | gammcor | `(G_aaaa + G_bbbb)/2` |
| `rdm2_bbbb.bin` | `.../RDM2_BBBB` | gammcor | `(G_aaaa + G_bbbb)/2` |
| `rdm2_abab.bin` | `.../RDM2_ABAB` | gammcor | `(G_abab + G_baba)/2` |
| `rdm2_baba.bin` | `.../RDM2_BABA` | gammcor | `(G_abab + G_baba)/2` |
| `rdm2_full.bin` | `.../RDM2` | pyscf | `Γ_aaaa + 2·Γ_abab + Γ_bbbb` |
| `rdm2_full_reordered.bin` | `.../RDM2_FULL_REORDERED` | gammcor | `2·(G_spin_avg + G_ab_avg)` |
| `rdm2.dat` | — | pyscf (ASCII) | `Γ_aaaa + Γ_abab + Γ_baba + Γ_bbbb` |
| `rdm1a.bin` | `.../RDM1_A` | — | `(γ_α + γ_β)/2` |
| `rdm1b.bin` | `.../RDM1_B` | — | `(γ_α + γ_β)/2` |
| — (no file) | `.../RDM1` | — | `γ_α + γ_β` |
| `rdm1p.bin` | `.../OCC_A` | other | `[1.0]×NI ⊕ diag(γ_α) ⊕ [0.0]×NV` |
| `rdm1m.bin` | `.../OCC_B` | other | `[1.0]×NI ⊕ diag(γ_β) ⊕ [0.0]×NV` |
| `rdm1.bin` | `POSTHF/OCC` | other | `mycas.mo_occ / 2` |

with `G_x = reorder_rdm(Γ_x)`, `G_baba = G_abab.transpose(1,0,3,2)`.

### 1b. `get_data_for_gammcor` — state-averaged branch (suffix `_<state>.<sym>`)

| array (legacy file) | HDF5 dataset | convention | expression |
|---|---|---|---|
| `rdm2_aaaa_S.bin` | `.../RDM2_AAAA` | gammcor | `(G_aaaa + G_bbbb)/2` |
| `rdm2_abab_S.bin` | `.../RDM2_ABAB` | gammcor | `(G_abab + G_baba)/2` |
| *(no `bbbb`/`baba` files)* | — | — | — |
| `rdm2_full_S.bin` | `.../RDM2` | pyscf | `Γ_aaaa + 2·Γ_abab + Γ_bbbb` |
| `rdm2_full_reordered_S.bin` | `.../RDM2_FULL_REORDERED` | gammcor | `2·(G_spin_avg + G_ab_avg)` |
| `rdm1a_S.bin` | `.../RDM1_A` | — | `γ_α` |
| `rdm1b_S.bin` | `.../RDM1_B` | — | `γ_β` |
| `rdm1_full_S.bin` | `.../RDM1` | — | `γ_α + γ_β` |
| *(no `rdm2.dat`, no `rdm1.bin`, no `rdm1p/m.bin`)* | | | |

### 1c. `get_rohf_for_gammcor`

| array (legacy file) | HDF5 dataset | convention | expression |
|---|---|---|---|
| `rdm2_rohf.dat` | `.../RDM2` | other (ASCII) | `X[I,J,K,L] = 2·Γ_JL,IK` from `n_α = 0.5`, `n_β = 0.0` |
| `mo_occ_int.bin` | `SCF/MO_OCC_INT` | other | `int32` `mo_occ` (0/1/2) |
| `occ_rohf.bin` | `SCF/MO_OCC` | other | `myhf.mo_occ` (0/1/2) |

---

## 2. `read_PYSCF` in `SOURCE/interface_pp.f90`

`only_full` is **`.true.` iff `JobType == AC0` and `ITwoEl > 1`**; `.false.` otherwise.
`only_full = .true.` first looks for `rdm2_full_reordered<sfx>.bin`; if that file is absent (older dumps made before the pp branch) it silently falls back to the spin blocks and assembles `rdm2_full` from them, exactly as in the `only_full = .false.` case.
`natural` is `natural_orb`, the 8th line of `auxdata*.txt` (= `int(mycas.natorb)`).

### 2.1 `natural = 0`, `only_full = .false.`  (the usual CASSCF path)

| variable | status | file / expression |
|---|---|---|
| `AuxData%rdm2_pp` | read | `rdm2_aaaa<sfx>.bin` |
| `AuxData%rdm2_pm` | read | `rdm2_abab<sfx>.bin` |
| `AuxData%rdm2_mm` | read | `rdm2_bbbb<sfx>.bin`  (fallback: `rdm2_aaaa<sfx>.bin`) |
| `AuxData%rdm2_mp` | read | `rdm2_baba<sfx>.bin`  (fallback: `rdm2_abab<sfx>.bin`) |
| `AuxData%rdm1_full` | read | `rdm1_full<sfx>.bin` |
| `AuxData%rdm1_p` | read | `rdm1a<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm1_m` | read | `rdm1b<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm2_full` | calculated | `rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp` |
| `AuxData%Occ(NI+1:NIA)` | calculated | `eig(rdm1_full)/2`, sorted descending |
| all `rdm2_*`, all `rdm1_*` | transformed | `Uᵀ Γ U` (MO→NO), `U = CMONO_NA` = eigenvectors of `rdm1_full` |
| `rdm2.dat` | saved | `write_rdm2_dat(rdm2_full, NA)` |

**Important:** `rdm2_full<sfx>.bin` is *not* read (the call is commented out). `NA` is shrunk to the orbitals with occupancy > 1e-8 and `NI` is recomputed from `NEL` before the transformation.

### 2.2 `natural = 0`, `only_full = .true.`  (AC0 + THC)

| variable | status | file / expression |
|---|---|---|
| `AuxData%rdm2_full` | read | `rdm2_full_reordered<sfx>.bin` |
| `AuxData%rdm1_full` | read | `rdm1_full<sfx>.bin` |
| `AuxData%rdm1_p` | read | `rdm1a<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm1_m` | read | `rdm1b<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%Occ(NI+1:NIA)` | calculated | `eig(rdm1_full)/2`, sorted descending |
| `rdm2_full`, `rdm1_p/m/full` | transformed | `Uᵀ Γ U` (MO→NO), `U = CMONO_NA` |
| `rdm2.dat` | saved | `write_rdm2_dat(rdm2_full, NA)` |

If `rdm2_full_reordered<sfx>.bin` is missing, the routine falls back to §2.1: it reads the four spin blocks and sets `rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp`.

**Important:** with the reordered file present, `rdm2_pp/mm/pm/mp` are never allocated and `canonicalize` uses `Γ = rdm2_full/2` for the active-space cumulant correction instead of `rdm2_pp + rdm2_pm`. On the fallback path they *are* allocated, so `canonicalize` takes the `rdm2_pp + rdm2_pm` branch — the two are equal because `rdm2_mm = rdm2_pp` and `rdm2_mp = rdm2_pm` in these dumps.

### 2.3 `natural = 1`, `only_full = .false.`

| variable | status | file / expression |
|---|---|---|
| `AuxData%rdm2_pp` | read | `rdm2_aaaa<sfx>.bin` |
| `AuxData%rdm2_pm` | read | `rdm2_abab<sfx>.bin` |
| `AuxData%rdm2_mm` | read | `rdm2_bbbb<sfx>.bin`  (fallback: `rdm2_aaaa<sfx>.bin`) |
| `AuxData%rdm2_mp` | read | `rdm2_baba<sfx>.bin`  (fallback: `rdm2_abab<sfx>.bin`) |
| `AuxData%rdm1_full` | read | `rdm1_full<sfx>.bin` |
| `AuxData%rdm1_p` | read | `rdm1a<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm1_m` | read | `rdm1b<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm2_full` | calculated | `rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp` |
| `AuxData%Occ(1:NIA)` | read | `rdm1.bin`  (fallback: `occ.bin`) |
| `AuxData%n_p`, `n_m` | calculated | inactive: `1`; active: `rdm1_p(i-NI,i-NI)`, `rdm1_m(i-NI,i-NI)`; virtual: `0` |
| `AuxData%n` | calculated | `n_p + n_m` |

**Important:** nothing is transformed — the RDMs stay exactly as read, and `rdm2.dat` is **not** written. `rdm1_full` is loaded but unused.

### 2.4 `natural = 1`, `only_full = .true.`

| variable | status | file / expression |
|---|---|---|
| `AuxData%rdm2_full` | read | `rdm2_full_reordered<sfx>.bin` |
| `AuxData%rdm1_full` | read | `rdm1_full<sfx>.bin` |
| `AuxData%rdm1_p` | read | `rdm1a<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%rdm1_m` | read | `rdm1b<sfx>.bin`  (fallback: `rdm1_full/2`) |
| `AuxData%Occ(1:NIA)` | read | `rdm1.bin`  (fallback: `occ.bin`) |
| `AuxData%n_p`, `n_m` | calculated | inactive: `1`; active: `rdm1_p/m(i-NI,i-NI)`; virtual: `0` |
| `AuxData%n` | calculated | `n_p + n_m` |

If `rdm2_full_reordered<sfx>.bin` is missing, the routine falls back to §2.3: it reads the four spin blocks and sets `rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp`.

**Important:** nothing is transformed, nothing is saved. With the reordered file present `rdm2_pp/mm/pm/mp` are never allocated and the cumulant correction uses `rdm2_full/2`.

---

## 3. `rdm2.dat` conventions

The GAMMCOR readers (`AB_CAS_*` in `interpa.f`, `read_2rdm` in `pp_utils.f90`) all do, for a line `I J K L X`:

```
RDM2Act(NAddrRDM(J,L,I,K,NAct)) = X/2      with  X = <E_IJ E_KL> − δ_JK <E_IL> = 2·Γ2(J,L,I,K)
```

| writer | line layout | resulting `RDM2Act` |
|---|---|---|
| **Python** `calc_full_dm2` | `i j k l  Γ_tot[i,j,k,l]` (pyscf order, all `NA⁴` lines) | `RDM2Act(p,q,r,s) = Γ[r,p,s,q]/2` |
| **Fortran** `write_rdm2_dat` | `j l i k  rdm2_full(i,j,k,l)` (gammcor-ordered array, lines with `|val| ≤ 1e-8` skipped) | `RDM2Act(p,q,r,s) = rdm2_full(s,r,q,p)/2` |

The two agree: with `rdm2_full(a,b,c,d) = Γ[a,c,b,d]`, both give `Γ[r,p,s,q]/2` after using `Γ[i,j,k,l] = Γ[k,l,i,j]`.

---

## 4. Notes

**N1 — `n_p`/`n_m` unset when `natural = 0`.** `load_npm` is called only from `load_occupancy` (the `natural = 1` branch). On the `natural = 0` path `AuxData%n_p`/`n_m` stay unallocated, yet `read_PYSCF` reads them in the `ITwoEl == 1` branch (`interface_pp.f90:232-238`).

**N2 — `HNO0` transformed twice when `natural = 0` and `ITwoEl == 1`.** `Trans2NO_AO` already does `HNO0 ← CAONOᵀ HNO0 CAONO` (`interface_pp.f90:2246-2249`), and `read_PYSCF` repeats it (`interface_pp.f90:210-212`).

**N3 — `rdm1_full.bin` is not written in the state-specific legacy branch.** The `dump.array(f'{sgrp}/RDM1', ...)` call has no `legacy=` argument, so `AuxData%rdm1_full` is left uninitialized — and it is what `Trans2NO_AO` diagonalizes when `natural = 0`.

**N4 — Written but never read by GAMMCOR:** `rdm2_full[_S].bin`, `rdm1p.bin`, `rdm1m.bin`, `rdm2_rohf.dat`.

**N5 — `HCore.bin` is always AO**, contrary to the comment in `load_one_electron_integrals`.

**N6 — `TWOEl.bin` basis differs between branches.** State-specific dumps MO/NO (`doMOtrans=True`), state-averaged dumps AO (`doMOtrans=False`). With `natural = 1` GAMMCOR performs no `TwoNO1` transformation, so SA + `natorb=1` leaves the ERIs in AO.

**N7 — MS ≠ 0.** The state-specific branch writes spin-averaged `aaaa`/`bbbb` and `abab`/`baba` pairs, so genuine spin polarization is lost. The SA branch keeps true `γ_α`, `γ_β` in `rdm1a/b` but still spin-averages the 2-RDM blocks.
