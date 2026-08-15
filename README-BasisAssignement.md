# `BasisAssignement` — mixed basis sets in GammCor

`BasisAssignement` lets you assign **different basis sets to different atoms** of the same
molecule. The classic single-basis input keeps working exactly as before.

---

## 1. The classic way (still supported)

Give the basis-set library directory and the basis file in the `Calculation` block:

```
Calculation
 BasisPath /home/../basis-folder
 Basis     cc-pvdz.txt
end
```

All atoms get `cc-pvdz.txt`. Nothing changed — old inputs run unmodified.

---

## 2. The new way — `BasisAssignement`

> **`BasisPath` is still required** and must stay in the `Calculation` block.
> If `BasisAssignement` is present, the `Basis` keyword is **ignored**.

```
Calculation
 BasisPath /home/../basis-folder
 Basis     cc-pvdz.txt      <-- ignored when BasisAssignement is present
end

BasisAssignement
 * cc-pvdz.txt
end
```

Each line inside `BasisAssignement` has the form:

```
<selector>  <basis-file>
```

where `<basis-file>` is a **file name inside `BasisPath`** (write it with its
extension, e.g. `cc-pvqz.txt`).

### Selectors

| Selector | Meaning | Example |
|----------|---------|---------|
| `*` | default basis for **every** atom | `* cc-pvdz.txt` |
| element symbol | all atoms of that element | `Fe cc-pvqz.txt` |
| atom number | one specific atom (line number in `xyz`, counting from 1) | `10 cc-pvtz.txt` |
| atom range | consecutive atoms `n-m` | `6-10 cc-pvtz.txt` |

### Priority

When several rules match an atom, the most specific one wins:

```
   atom number / range   >   element symbol   >   *  (default)
        (highest)                                  (lowest)
```

Order of the lines inside the block does not matter.

---

## 3. Examples

### 3.1 One common basis for all atoms

```
BasisAssignement
 * cc-pvdz.txt
end
```

### 3.2 Default basis + special basis for chosen elements

```
BasisAssignement
 * cc-pvdz.txt
 Fe cc-pvqz.txt
end
```

or, with more elements:

```
BasisAssignement
 * cc-pvdz.txt
 Fe cc-pvqz.txt
 S  cc-pvtz.txt
end
```

### 3.3 Default + per-element + one individual atom

```
BasisAssignement
 * cc-pvdz.txt
 S  cc-pvtz.txt
 Fe cc-pvqz.txt
 10 cc-pvtz.txt
end
```

`10` is the **line number of the atom in the `xyz` block**, counting from 1:

```
xyz
10
H  -0.3060000 -1.9589999 -1.8649999
S  -1.3759939 -1.3483813 -1.3382063
H   0.1980000  1.9369999 -1.8409999
S   1.1355853  1.1515255 -1.2925355
H  -0.1120000 -1.9629999  2.0249998
S   0.9297809 -1.4714457  1.3418066
H   0.5120000  1.8939999  2.0439998
S  -0.7443586  1.6232089  1.6653868
Fe -0.3157074  0.0326084  0.1116605
H  -1.8552819  0.2595304  0.2109527
end
```

Resulting assignment:

| atom # | element | basis | rule that applied |
|--------|---------|-------|-------------------|
| 1, 3, 5, 7 | H | `cc-pvdz.txt` | `*` (default) |
| 2, 4, 6, 8 | S | `cc-pvtz.txt` | element rule `S` |
| 9 | Fe | `cc-pvqz.txt` | element rule `Fe` |
| **10** | H | **`cc-pvtz.txt`** | atom rule `10` — beats the default |

### 3.4 Ranges of atoms

If the atoms you want to single out sit next to each other in the `xyz` block,
use a range `n-m`:

```
xyz
10
Fe -0.3157074  0.0326084  0.1116605
H  -0.3060000 -1.9589999 -1.8649999
H   0.1980000  1.9369999 -1.8409999
H  -0.1120000 -1.9629999  2.0249998
H   0.5120000  1.8939999  2.0439998
S  -1.3759939 -1.3483813 -1.3382063
S   1.1355853  1.1515255 -1.2925355
S  -0.7443586  1.6232089  1.6653868
S   0.9297809 -1.4714457  1.3418066
H  -1.8552819  0.2595304  0.2109527
end
```

```
BasisAssignement
 * cc-pvdz.txt
 1 cc-pvqz.txt
 6-10 cc-pvtz.txt
end
```

| atom # | element | basis |
|--------|---------|-------|
| 1 | Fe | `cc-pvqz.txt` |
| 2–5 | H | `cc-pvdz.txt` (default) |
| 6–9 | S | `cc-pvtz.txt` |
| 10 | H | `cc-pvtz.txt` |

---

