# Geometry input in GammCor

The molecular geometry can be written **directly in the GammCor input file**, as
before, or kept in a **separate `.xyz` file** that the input only points to.
Both forms are read by the same `xyz` directive.

---

## 1. Geometry inside the input file (classic)

Open an `xyz` block, give the number of atoms, then one line per atom, and close
the block with `end`:

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

Coordinates are in **Ångström** by default.

---

## 2. Geometry in a separate xyz file

Put the file name **on the same line as the `xyz` keyword**:

```
xyz geom.xyz
end
```

A path may be relative or absolute:

```
xyz /home/.../structures/geom.xyz
end
```

`geom.xyz` is then a **plain, standard xyz file** — atom count, a comment line
(may be empty), and the atom lines:

```
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
```

This is exactly what Avogadro, VMD, Molden, OpenBabel, ORCA, … write out, so
structures can be used without editing.


