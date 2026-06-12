# Case 1: bcc Fe

*Download the runnable workflow below as a
[Jupyter notebook](case1_bcc_fe.ipynb) or a [plain script](case1_bcc_fe.jl).*

This worked example builds a spin-cluster expansion (SCE) model for
body-centered cubic (bcc) iron and walks through the full workflow, from the
reference data to the fitted exchange interactions.

!!! note
    This page is a skeleton. The narrative, code, and numbers are added in a
    later revision.

## Overview

Describe the bcc Fe system: lattice, magnetic ground state, and why it is a good
first case.

(to be written)

## Preparing the inputs

Before running Magesty you prepare the noncollinear spin-DFT reference data and
the configuration files. This section explains what to set up and why, including
the DFT input files (e.g. `INCAR`, `POSCAR`) used to generate the reference
energies and torques, and how they map onto the Magesty input TOML and the
`EMBSET` training-data file.

(to be written)

The remaining sections run the Magesty workflow itself.
