# Case 1: bcc Fe

*Download the runnable workflow below as a
[Jupyter notebook](case1_bcc_fe.ipynb) or a [plain script](case1_bcc_fe.jl).*

This worked example builds a spin-cluster expansion (SCE) model for
body-centered cubic (bcc) iron and walks through the full workflow, from the
noncollinear spin-DFT reference data to the fitted exchange interactions.

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

## Sampling spin configurations

!!! info
    This step uses the `magesty` command-line tool; if it is not set up yet, see
    [Installation](@ref "Command-Line Interface").

With the template `INCAR` in hand, draw the spin configurations to compute. We
sample from the mean-field-approximation (MFA) thermal distribution at the
reduced temperature ``\tau = T / T_{\mathrm{C}}^{\mathrm{MFA}} = 0.05``. This is
deep in the low-temperature regime, far below the mean-field Curie temperature,
so the drawn configurations are small fluctuations about the ferromagnetic ground
state. We draw **50 configurations**:

```bash
magesty vasp mfa INCAR tau --start 0.05 --stop 0.05 --num-points 1 --num-samples 50
```

This writes one INCAR per drawn configuration (`sample-NN.INCAR`), each with
`MAGMOM` and `M_CONSTR` set to the sampled directions and all other keys copied
from the template. The arguments control the sweep:

- `--start`, `--stop`, and `--num-points` set the values of the control variable
  (here ``\tau``); the sweep is `range(start, stop; length = num_points)`. With
  `--start 0.05 --stop 0.05 --num-points 1` this collapses to the single
  temperature ``\tau = 0.05``.
- `--num-samples` is the number of configurations drawn **per** sweep value, so
  `--num-samples 50` yields 50 configurations at ``\tau = 0.05`` (the total is
  `num_points` × `num_samples`).

See
[Mean-Field Sampling](@ref "Mean-Field Sampling: Generating Thermal Spin Configurations")
for the underlying distribution and the role of ``\tau``.

## Running the DFT calculations

Run VASP on each `sample-NN.INCAR` on your own compute resources (a supercomputer
or cluster) to obtain the energy and the constraining torques for every
configuration. As noted above, **do not forget to use a sufficiently large
`LAMBDA`**: converge each run first with the small `LAMBDA = 1`, then restart with
a large value (e.g. `LAMBDA = 50`, reading `CHGCAR` via `ICHARG = 1`) so that the
torques are reliable.

## From VASP to Magesty

The reference data and the structure feed the SCE workflow as two files: the
`EMBSET` training data and the input TOML.

### Building the EMBSET

Each calculation produces an `OSZICAR` holding the energy and the constraining
field. Convert the runs into a single `EMBSET` training-data file with
`magesty vasp embset`, in either of two ways.

- **Convert each run, then merge.** Convert one `OSZICAR` at a time,

  ```bash
  magesty vasp embset OSZICAR --output EMBSET
  ```

  and combine the contents of all the resulting `EMBSET` files into one.

- **Convert all at once.** Gather the `OSZICAR` files in one place, e.g. renamed
  `01.oszicar`, `02.oszicar`, …, `50.oszicar`, and convert them in a single call.
  Each file becomes one configuration block, numbered in the given order:

  ```bash
  magesty vasp embset *.oszicar --output EMBSET
  ```

!!! tip "Skip the DFT for this example"
    The 50 constrained calculations are expensive. To follow along without
    running them, download the precomputed [`EMBSET`](case1_inputs/EMBSET) for
    this example and use it directly in the steps below.

### Generating the input TOML

Build the Magesty input TOML from `POSCAR` with `magesty vasp toml`:

```bash
magesty vasp toml POSCAR --output input.toml
```

This fills the `[general]`, `[symmetry]`, `[interaction]`, and `[structure]`
tables from the structure. The interaction settings are written as placeholders,
so **edit them** to the basis you want before fitting. For this example the
edited file is:

```toml:input.toml
[general]
name = "Fe444"
nat  = 128
kd   = ["Fe"]
periodicity = [true, true, true]

[symmetry]
tolerance = 1.0e-5
isotropy  = true

[interaction]
nbody = 2

[interaction.body1]
lmax.Fe = 0

[interaction.body2]
cutoff."Fe-Fe" = -1
lsum = 2

[structure]
lattice = [
  [11.3192, 0.0,     0.0],
  [0.0,     11.3192, 0.0],
  [0.0,     0.0,     11.3192],
]
kd_list  = [1, 1, 1, 1, ...]       # 128 entries, element index per atom
position = [
  [0.000, 0.000, 0.000],
  [0.125, 0.125, 0.125],
  # ... (128 atoms, direct coordinates)
]
```

- `[general]` — the system name, atom count (`nat = 128`), the element list
  `kd` (a single species here), and periodicity in all three directions.
- `[symmetry]` — `tolerance` for spacegroup detection and `isotropy = true`,
  which restricts the basis to rotationally invariant (``L_f = 0``) terms, i.e.
  isotropic exchange.
- `[interaction]` — `nbody = 2` keeps up to pair terms. `body1`'s `lmax.Fe = 0`
  is the on-site angular-momentum order; `body2`'s `lsum = 2` is the cutoff on
  the summed angular momentum of a pair basis function, and
  `cutoff."Fe-Fe"` is the Fe–Fe pair cutoff radius in Å (`-1` includes all pairs,
  i.e. no distance cutoff). Together these define an isotropic two-body (pair)
  model.
- `[structure]` — the supercell `lattice` (the 11.3192 Å cube), `kd_list` (the
  element index of each atom), and the 128 atomic `position`s in direct
  coordinates, all taken straight from `POSCAR`.

See [Input Keys](../input_keys.md) for the full key reference and the
[Tools](../tools.md) page for the converter options. The edited
[`input.toml`](case1_inputs/input.toml) is downloadable.

With the input TOML and the combined `EMBSET` in hand, the remaining sections run
the Magesty workflow itself.
