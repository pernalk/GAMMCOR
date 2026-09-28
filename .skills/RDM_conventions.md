# RDMs in `Pyscf2Gammcor_mini.py` and `interface_pp.f90`

## 0. Convention names used in this document

| Name | Definition |
|---|---|
| **pyscf** | `0123` — exactly as `fcisolver.make_rdm12(s)` returns it: `Γ[p,q,r,s] = <E_pq E_rs> − δ_qr <E_ps>` (chemists' pair order `(pq|rs)`). |
| **gammcor** | `0213` — `dm2.transpose(0,2,1,3)`, i.e. `G[p,q,r,s] = Γ[p,r,q,s]`. This is `reorder_rdm()` in the script, applied to every `rdm2_*.bin` except `rdm2_full*.bin`. |
| **other** | anything that is not a 4-index 2-RDM (occupation vectors, ASCII lists, ROHF model RDM). |

**Storage note (matters for the 4-index arrays).** All `.bin` files are raw `numpy.tofile()` streams: no header, little-endian, C-contiguous, `float64` (only `mo_occ_int.bin` is `int32`). Fortran `access='stream'` reads them column-major, so the index order is *reversed* on read: `F(a,b,c,d) = numpy[d,c,b,a]`.
For the 2-index arrays the script pre-transposes on dump (`CAONO.T`, `HCore.T`), so Fortran gets them un-transposed; the RDM1 blocks are symmetric so it does not matter.
For the 4-index arrays there is no pre-transpose, but because `Γ[i,j,k,l] = Γ[k,l,i,j] = Γ[j,i,l,k]` (real wavefunction), the reversal composed with the `0213` python transpose gives back exactly `0213`: the Fortran array equals `Γ[a,c,b,d]`. Net effect is the intended **gammcor** convention.

---

## 1. Arrays produced by `Pyscf2Gammcor_mini.py`

### 1a. `get_data_for_gammcor` — state-specific branch (`mycas.weights is None`)

| array (legacy file) | HDF5 dataset | convention | notes |
|---|---|---|---|
| `rdm2_aaaa.bin` | `.../RDM2_AAAA` | gammcor | `(G_aaaa + G_bbbb)/2` — **spin-averaged** |
| `rdm2_bbbb.bin` | `.../RDM2_BBBB` | gammcor | *same array as* `rdm2_aaaa.bin` (byte-identical) |
| `rdm2_abab.bin` | `.../RDM2_ABAB` | gammcor | `(G_abab + G_baba)/2`, with `G_baba = G_abab.transpose(1,0,3,2)` — **spin-averaged** |
| `rdm2_baba.bin` | `.../RDM2_BABA` | gammcor | *same array as* `rdm2_abab.bin` |
| `rdm2_full.bin` | `.../RDM2` | **pyscf** | `Γ_aaaa + 2·Γ_abab + Γ_bbbb` (spin-traced). **Never read by GAMMCOR** — the read is commented out in `load_density_matrices` |
| `rdm2_full_reordered.bin` | `.../RDM2_FULL_REORDERED` | gammcor | `2·(spin_avg + ab_avg)` = same spin-traced sum, reordered |
| `rdm2.dat` | — | pyscf (ASCII) | `calc_full_dm2`: 1-based `i j k l val`, `val = Γ_aaaa+Γ_abab+Γ_baba+Γ_bbbb`, **all** `NA⁴` lines (no threshold) |
| `rdm1a.bin` | `.../RDM1_A` | — | `(γ_α + γ_β)/2` — **spin-averaged** |
| `rdm1b.bin` | `.../RDM1_B` | — | *same array as* `rdm1a.bin` |
| — (no file) | `.../RDM1` | — | `γ_α + γ_β`. **HDF5 only — `rdm1_full.bin` is not written in legacy mode** (see note N3) |
| `rdm1p.bin` | `.../OCC_A` | other | length-`NBasis` vector: `1.0` inactive, `diag(γ_α)` active, `0` virtual. **Never read by GAMMCOR** |
| `rdm1m.bin` | `.../OCC_B` | other | same with `diag(γ_β)`. **Never read by GAMMCOR** |
| `rdm1.bin` | `POSTHF/OCC` | other | `mycas.mo_occ / 2` (per-spin occupations), length `NBasis` |

### 1b. `get_data_for_gammcor` — state-averaged branch (suffix `_<state>.<sym>`)

| array (legacy file) | HDF5 dataset | convention | notes |
|---|---|---|---|
| `rdm2_aaaa_S.bin` | `.../RDM2_AAAA` | gammcor | `(G_aaaa + G_bbbb)/2` — spin-averaged |
| `rdm2_abab_S.bin` | `.../RDM2_ABAB` | gammcor | `(G_abab + G_baba)/2` — spin-averaged |
| *(no `bbbb`/`baba` files)* | — | — | GAMMCOR falls back: `rdm2_mm ← rdm2_aaaa`, `rdm2_mp ← rdm2_abab` |
| `rdm2_full_S.bin` | `.../RDM2` | **pyscf** | `Γ_aaaa + 2·Γ_abab + Γ_bbbb`. Never read |
| `rdm2_full_reordered_S.bin` | `.../RDM2_FULL_REORDERED` | gammcor | `2·(spin_avg + ab_avg)` |
| `rdm1a_S.bin` | `.../RDM1_A` | — | **true `γ_α` of the state — NOT spin-averaged** (differs from 1a) |
| `rdm1b_S.bin` | `.../RDM1_B` | — | **true `γ_β`** |
| `rdm1_full_S.bin` | `.../RDM1` | — | `γ_α + γ_β` |
| *(no `rdm2.dat`, no `rdm1.bin`, no `rdm1p/m.bin`)* | | | |

### 1c. `get_rohf_for_gammcor`

| array (legacy file) | HDF5 dataset | convention | notes |
|---|---|---|---|
| `rdm2_rohf.dat` | `.../RDM2` | other (ASCII) | `X[I,J,K,L] = 2·Γ_JL,IK` built from effective occupations `n_α=0.5`, `n_β=0.0` on the `NA` singly-occupied orbitals. **Never read by GAMMCOR** |
| `mo_occ_int.bin` | `SCF/MO_OCC_INT` | other | **`int32`**, length `NBasis`, values `0/1/2`. Read straight into `AuxData%IndAux` |
| `occ_rohf.bin` | `SCF/MO_OCC` | other | `myhf.mo_occ` (0/1/2), `float64`, length `NBasis` |

### 1d. Non-RDM payload (for completeness)

| file | contents | format |
|---|---|---|
| `C.bin` | `CAONO.T` (`mo_coeff` transposed on dump ⇒ Fortran reads AO×MO) | `NBasis²` f64 |
| `HCore.bin` | `HCore.T` — **always AO basis**, both for `natorb=0` and `natorb=1` | `NBasis²` f64 |
| `TWOEl.bin` | ERI packed with `NAddr3`; **AO basis** in the SA branch (`doMOtrans=False`), **MO/NO basis** in the state-specific branch (`doMOtrans=True`) | `ninte2` f64 |
| `auxdata[_S].txt` | `NBasis, NI, NA, NV, ECAS, ENuc, NEL, natorb[, frozen]` — one per line | ASCII |
| `auxdata_rohf.txt` | `NBasis, NI, NA, NV, EROHF, ENuc, NEL` | ASCII |

---

## 2. `read_PYSCF` in `SOURCE/interface_pp.f90`

`only_full` is **`.true.` iff `JobType == AC0` and `ITwoEl > 1`**; `.false.` otherwise.
`only_full = .true.` first looks for `rdm2_full_reordered<sfx>.bin`; if that file is absent (older dumps made before the pp branch) it silently falls back to the spin blocks and assembles `rdm2_full` from them, exactly as in the `only_full = .false.` case.
`natural` is `natural_orb`, read as the 8th line of `auxdata*.txt` (= `int(mycas.natorb)`).

### 2.0 Common to all four cases

| step | file → variable | convention on arrival |
|---|---|---|
| `load_cas_aux_data` | `auxdata<sfx>.txt` → `NBasis, NI, NA, NV, ECAS, ENuc, NEL, natural_orb`; sets `NIA = NI+NA` | — |
| `read_PYSCF` (if `ITwoEl==1`) | `TWOEl.bin` → `TwoEl(:)` | `NAddr3` packed, AO or MO per §1d |
| `load_one_electron_integrals` | `HCore.bin` → `AuxData%HNO0`; `THCData%HNO = HNO0` | AO basis |
| `load_density_matrices` | `rdm1_full<sfx>.bin` → `rdm1_p = rdm1_m = rdm1_full/2` | fallback only, if the `a`/`b` files are missing |
| `load_density_matrices` | `rdm1a<sfx>.bin` → `rdm1_p`; `rdm1b<sfx>.bin` → `rdm1_m` | as dumped |

### 2.1 `natural = 0`, `only_full = .false.`  (the usual CASSCF path)

| what | detail |
|---|---|
| **read (orbitals)** | `C.bin` → `CAOMO`. `CAONO` allocated, filled later. |
| **read (2-RDM)** | `rdm2_aaaa<sfx>.bin`→`rdm2_pp`, `rdm2_abab<sfx>.bin`→`rdm2_pm`, `rdm2_bbbb<sfx>.bin`→`rdm2_mm` (fallback `aaaa`), `rdm2_baba<sfx>.bin`→`rdm2_mp` (fallback `abab`) — all **gammcor** convention |
| **read (1-RDM)** | `rdm1_full<sfx>.bin`→`rdm1_full`, `rdm1a`→`rdm1_p`, `rdm1b`→`rdm1_m` |
| **calculated** | `rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp` |
| **transformed** (`Trans2NO_AO`) | diagonalize `rdm1_full` (NA block) → eigenvalues sorted desc., `/2` → `AuxData%Occ(NI+1:NIA)`; eigenvectors → `CMONO_NA`. `NA` shrunk to orbitals with occ > 1e-8; `NI` recomputed from `NEL`; `NV`, `NIA`, `IndAux` rebuilt. Then **MO→NO** rotation with `CMONO_NA` of: `rdm2_pp, rdm2_mm, rdm2_pm, rdm2_mp, rdm2_full, rdm1_p, rdm1_m, rdm1_full`. `CAONO = CAOMO · CMONO`. If `ITwoEl==1`: `TwoEl` AO→NO via `TwoNO1(CAONOᵀ)`, and `HNO0 ← CAONOᵀ HNO0 CAONO`. |
| **saved** | `write_rdm2_dat(rdm2_full, NA)` → **`rdm2.dat`** (see §3) |
| **not set** | `n_p, n_m, n` are **never allocated** on this path (`load_npm` is not called) — see note N1 |
| **then** | `ITwoEl>1`: `THC_init2`, `HNO0 ← HNO0_THC`. `ITwoEl==1`: `HNO0 ← CAONOᵀ HNO0 CAONO` again (note N2), `check_energy_incore` with `spin_sep=.false.` (uses `rdm2_full`) and `.true.` (uses `rdm2_pp/mm/pm/mp`) |

### 2.2 `natural = 0`, `only_full = .true.`  (AC0 + THC)

| what | detail |
|---|---|
| **read (orbitals)** | `C.bin` → `CAOMO` |
| **read (2-RDM)** | `rdm2_full_reordered<sfx>.bin` → `rdm2_full` (**gammcor**); `rdm2_pp/mm/pm/mp` are *not allocated*. If that file is missing, the four spin blocks are read instead (as §2.1) and `rdm2_full` is assembled from them |
| **read (1-RDM)** | as §2.0 (`rdm1_full`, `rdm1a`→`rdm1_p`, `rdm1b`→`rdm1_m`) |
| **calculated** | nothing (`rdm2_full` comes from file) |
| **transformed** | `Trans2NO_AO` as in §2.1, but the `rdm2_pp/mm/pm/mp` rotations are skipped (guarded by `allocated(...)`); `rdm2_full`, `rdm1_p/m/full` are rotated MO→NO |
| **saved** | `rdm2.dat` (same call) |
| **then** | always `ITwoEl>1` here → `THC_init2`, `HNO0 ← HNO0_THC`. `canonicalize` uses `Γ = rdm2_full/2` for the cumulant correction (instead of `rdm2_pp+rdm2_pm`). `check_energy_incore` is **not** called |

### 2.3 `natural = 1`, `only_full = .false.`

| what | detail |
|---|---|
| **read (orbitals)** | `C.bin` → **`CAONO`** directly (`CAOMO` never allocated) |
| **read (2-RDM)** | same four files as §2.1 → `rdm2_pp, rdm2_pm, rdm2_mm, rdm2_mp` |
| **read (1-RDM)** | `rdm1_full`, `rdm1_p`, `rdm1_m` |
| **calculated** | `rdm2_full = rdm2_pp + rdm2_mm + rdm2_pm + rdm2_mp` |
| **transformed** | **nothing** — no `Trans2NO_*`, RDMs stay exactly as read |
| **read (occupations)** | `load_occupancy`: `rdm1.bin` (fallback `occ.bin`) → `Occ(1:NIA)`; `IndAux` = 0/1/2 by block; `load_npm` → `n_p(i)=n_m(i)=1` inactive, `= rdm1_p/m(i-NI,i-NI)` active, `n = n_p + n_m` |
| **saved** | **nothing** — `rdm2.dat` is *not* written |
| **not used** | `rdm1_full` is loaded but unused on this path |
| **then** | `ITwoEl>1`: `THC_init2`. `ITwoEl==1`: `HNO0 ← CAONOᵀ HNO0 CAONO` (correct here — `HCore.bin` is AO), `check_energy_incore` both variants, `spinsep` flag from `|n_p−n_m|` |

### 2.4 `natural = 1`, `only_full = .true.`

| what | detail |
|---|---|
| **read (orbitals)** | `C.bin` → `CAONO` |
| **read (2-RDM)** | `rdm2_full_reordered<sfx>.bin` → `rdm2_full`; if missing, the four spin blocks are read instead (as §2.3) and `rdm2_full` is assembled from them |
| **read (1-RDM)** | `rdm1_full`, `rdm1_p`, `rdm1_m` |
| **calculated / transformed** | nothing |
| **read (occupations)** | `load_occupancy` + `load_npm`, exactly as §2.3 |
| **saved** | nothing |
| **then** | `ITwoEl>1` → `THC_init2`; cumulant correction from `rdm2_full/2`; no energy check |

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

## 4. Notes / things to be aware of

**N1 — `n_p`/`n_m` unset when `natural = 0`.** `load_npm` is only called from `load_occupancy` (the `natural = 1` branch) and from `calc_occupancy_orca`. On the `natural = 0` path `AuxData%n_p`/`n_m` stay unallocated, yet `read_PYSCF` reads them in the `ITwoEl == 1` branch (`interface_pp.f90:232-238`, the `spinsep` loop and the `occupancy` printout).

**N2 — `HNO0` transformed twice when `natural = 0` and `ITwoEl == 1`.** `Trans2NO_AO` already does `HNO0 ← CAONOᵀ HNO0 CAONO` inside its `if (Flags%ITwoel == 1)` block (`interface_pp.f90:2246-2249`), and `read_PYSCF` repeats the same two lines afterwards (`interface_pp.f90:210-212`, under the comment *"FCIDUMP modification: skip double transformation"* that is not actually skipping it).

**N3 — `rdm1_full.bin` is not written in the state-specific legacy branch.** `dump.array(f'{sgrp}/RDM1', dm1s[0]+dm1s[1])` is called without a `legacy=` argument, so with `dump_hdf5=False` no `rdm1_full.bin` appears. `load_density_matrices` then leaves `AuxData%rdm1_full` **uninitialized**, and `Trans2NO_AO` diagonalizes it on the `natural = 0` path. The SA branch does pass the legacy name and is fine.

**N4 — Files written but never read by GAMMCOR:** `rdm2_full[_S].bin` (the read is commented out), `rdm1p.bin`, `rdm1m.bin`, `rdm2_rohf.dat`. Also `INTS/OVERLAP`, `INTS/KINETIC`, `INTS/POTENTIAL` in HDF5.

**N5 — `HCore.bin` is always AO**, contrary to the comment in `load_one_electron_integrals` (*"If natural = 1, HCore are in NO basis"*). Both Python branches dump `myhf.get_hcore().T`.

**N6 — `TWOEl.bin` basis differs between branches.** State-specific uses `doMOtrans=True` (integrals already in the `CAONO` basis); state-averaged uses `doMOtrans=False` (AO). With `natural = 1` no `TwoNO1` transformation is performed in GAMMCOR, so the SA + `natorb=1` combination would leave the ERIs in AO while everything else is in the NO basis.

**N7 — MS ≠ 0.** Both branches print a warning when `‖Γ_aaaa − Γ_bbbb‖² > 1e-5`, but the state-specific branch still writes the spin-averaged `aaaa`/`bbbb` and `abab`/`baba` pairs, so genuine spin polarization is lost in the `.bin` files there. The SA branch keeps true `γ_α`, `γ_β` in `rdm1a/b` but still spin-averages the 2-RDM blocks.
