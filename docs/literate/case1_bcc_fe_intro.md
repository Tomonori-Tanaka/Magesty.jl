# Case 1: bcc Fe

*Download the runnable workflow below as a
[Jupyter notebook](case1_bcc_fe.ipynb) or a [plain script](case1_bcc_fe.jl).*

This worked example builds a spin-cluster expansion (SCE) model for
body-centered cubic (bcc) iron and walks through the full workflow, from the
noncollinear spin-DFT reference data to the fitted exchange interactions.

!!! note
    This page is a work in progress. The Magesty workflow sections below are
    still being written.

## Overview

The reference system is bcc Fe with conventional lattice constant
``a = 2.8298`` Å. We use a ``4 \times 4 \times 4`` supercell of the conventional
(two-atom) cell, i.e. **128 atoms**. The training data is generated with
noncollinear, spin-DFT calculations in which the **direction** of each local
moment is constrained, so that for an arbitrary spin configuration we obtain the
total energy and the per-atom magnetic torque. These (energy, torque) pairs are
exactly what the SCE design matrix is fitted against.

## Preparing the inputs

The reference data comes from VASP calculations run with the constrained-direction
(penalty) method: the moment direction on every atom is held fixed while the
electronic structure relaxes, and the resulting energy and constraining torques
are recorded. Sampling many spin configurations (here via the mean-field
sampling at ``\tau = 0.05``) and collecting their energies and torques produces
the training set.

You prepare four VASP input files. The structure (`POSCAR`) and the run settings
(`INCAR`, `KPOINTS`) are shown here in abridged form and are downloadable in
full; the pseudopotential is described in words because of licensing.

Download the complete input files:
[`INCAR`](case1_inputs/INCAR), [`KPOINTS`](case1_inputs/KPOINTS),
[`POSCAR`](case1_inputs/POSCAR).

### POSCAR — structure

The ``4 \times 4 \times 4`` bcc Fe supercell.

```text:POSCAR
Fe444
2.8298
4.0  0.0  0.0
0.0  4.0  0.0
0.0  0.0  4.0
   Fe
  128
Direct
     0.000000000         0.000000000         0.000000000
     0.125000000         0.125000000         0.125000000
     0.000000000         0.000000000         0.250000000
     ...
     (128 atoms in total)
```

### POTCAR — pseudopotential

The pseudopotential file is not redistributed here for licensing reasons. These
calculations used the PBE PAW **Fe potential with 8 valence electrons**.

### KPOINTS — k-point mesh

A ``\Gamma``-centered ``3 \times 3 \times 3`` mesh for the supercell.

```text:KPOINTS
Fetest
0
Gamma
3 3 3
0 0 0
```

### INCAR — run settings

A single static (`NSW = 0`, `IBRION = -1`) noncollinear calculation without
spin-orbit coupling (`LNONCOLLINEAR = .TRUE.`, `LSORBIT = .FALSE.`). The moment
directions are constrained with the penalty method
(`I_CONSTRAINED_M`, penalty strength `LAMBDA`, integration radius
`RWIGS`); `M_CONSTR` sets the target direction on each atom. Symmetry is switched
off (`ISYM = 0`), as required for arbitrary noncollinear configurations.

The `M_CONSTR` and `MAGMOM` lines below show the collinear ``+z`` reference (all
moments along ``+z``); during sampling these per-atom directions are replaced for
each spin configuration.

```text:INCAR
NCORE = 4

!ICHARG = 1

ENCUT = 350.0
PREC = accurate
NBANDS = 1700
GGA = PE
GGA_COMPAT = False
LASPH = True

EDIFF = 1.28E-6
NELM = 60

IMIX = 4
AMIX      = 0.05
BMIX      = 0.001
AMIX_MAG  = 0.05
BMIX_MAG  = 0.001
LMAXMIX = 4

NSW = 0
IBRION = -1

LREAL = .FALSE.
LWAVE = .False.
LCHARG = .TRUE.

LORBIT = 0

LNONCOLLINEAR = .TRUE.
LSORBIT = .FALSE.
I_CONSTRAINED_M = 4
RWIGS = 1.22533
LAMBDA = 1
M_CONSTR = 0 0 3.0   0 0 3.0   ...   (one triplet per atom, 128 in total)

ISYM = 0
MAGMOM = 0 0 3.0   0 0 3.0   ...   (one triplet per atom, 128 in total)
```

!!! warning "Increase `LAMBDA` on restart"
    The INCAR above sets a small penalty, `LAMBDA = 1`, only so that the first
    SCF converges stably; starting directly with a large penalty often runs into
    convergence trouble. Once this run has converged, **do not forget to restart**
    the calculation with a much larger penalty (e.g. `LAMBDA = 50`), reading the
    converged charge density by setting `ICHARG = 1` (which reads `CHGCAR`). The
    constraining torques are only
    reliable once the direction constraint is tightly enforced. See
    [Penalty Term Dependence](@ref) for how a too-small `LAMBDA` degrades the
    fitted model.

### From VASP to Magesty

The structure and the reference data feed the SCE workflow through two files:

- the **input TOML**, which encodes the structure and the interaction settings
  (built from `POSCAR`), and
- the **`EMBSET`** file, which holds the sampled spin configurations together
  with their DFT energies and torques.

See the [Tools](@ref) page for the `magesty vasp toml` and `magesty vasp embset`
converters that produce these from the VASP outputs.

The remaining sections run the Magesty workflow itself.
