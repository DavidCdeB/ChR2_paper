# ChR2 paper

<div align="center">
  <img src="./TOC-ChR2-multiphoton.png" width="400px" alt="TOC" />
  <p><strong>Multiphoton Absorption Spectra of Channelrhodopsin-2 via Multiscale Simulation Methods</strong></p>
</div>

This repository hosts the scripts and files used to reproduce the data presented here:

**Carrasco-Busturia D***, Linares M, Norman P, Haugaard Olsen JM. Multiphoton Absorption Spectra of Channelrhodopsin-2 via Multiscale Simulation Methods. ChemRxiv. 2025; https://doi.org/10.26434/chemrxiv-2025-vzz1r-v2 

 
## Table of Contents

- [Merging trajectories and centering to QM region](#merging-trajectories-and-centering-to-qm-region)
- [MiMiC files](#mimic-files)
- [Split trajectory into pqr files](#split-trajectory-into-pqr-files)
- [Creation of potentials](#creation-of-potentials)
    - [mfcc approach](#mfcc-approach)
    - [cp3 approach](#cp3-approach)
        - [QM-sampled moiety; LYR-472](#qm-sampled-moiety-lyr-472)
        - [MM-sampled moiety; LYR-225](#mm-sampled-moiety-lyr-225)
- [Autocorrelation function](#autocorrelation-function)
- [Placing ECPs](#placing-ecps)
- [OPA, TPA, 3PA dal files](#opa-tpa-3pa-dal-files)


## Merging trajectories and centering to QM region

The first step in the workflow requires to:

1. Merge the trajectories comming from different restarts to a single trajectory, and 

2. Center the resulting trajectory with respect to the QM-sampled moiety (LYR-472) so that we can run the PE framwork (Fragment-based PE and PE-TD-DFT MPA spectra, steps 4 and 5 in the workflow, respectively). 

Follow the instructions [here](./Merging_and_centering_trajectory.ipynb).

## MiMiC files

Inside the MiMiC directory, you may find the input files to run the Berendsen equilibration and Nosé-Hoover production:

```bash
MiMiC_files/
└── Berendsen_equilibration
    ├── equilibration.inp
    └── Nose-Hoover_production
        └── equilibration.inp
```



## Split trajectory into pqr files:

From the QM/MM MD trajectory, we will generate a .pqr file for each snapshot. To ensure compatibility with PyFraME, we use a forked version of MDAnalysis that provides the required .pqr format. 

```bash
mkdir chr2-paper && chr2-paper
git clone git@github.com:DavidCdeB/mdanalysis.git
```

Create a virtual environment and activate it:

```bash
cd ../
python -m venv mda-venv
. mda-venv/bin/activate
```

Install:

```bash
python -m pip install ./chr2-paper/mdanalysis
```

Within the activated the virtual environment, now run [trr_to_pqr_snapshots.py](./trr_to_pqr_snapshots.py)

```bash
python trr_to_pqr_snapshots.py
```

This will generate the 200 .pqr files necessary to reproduce the results presented in this paper.

## Autocorrelation function

The script needed to create the autocorrelation function of the TPA cross sections: [tpa_acf.py](./tpa_acf.py)


## Creation of potentials

The creation of the `.pot` files has been generated through two approaches (see publication for full details):

- The __Molecular Fractionation with Conjugate Caps__ approach, i.e. "*mfcc*"
- __Standard potenitals__,  i.e. . "*cp3*" 


### mfcc approach

Execute one after the other, the following scripts:

1. [mfcc_ChR2_and_close-water_create_inputs.py](./mfcc_ChR2_and_close-water_create_inputs.py) creates the input fragment files

2. [mfcc_ChR2_and_close-water_write_potential.py](./mfcc_ChR2_and_close-water_write_potential.py) writes the embedding potential in the form of a .pot file


### cp3 approach

This approach was used for production, where, as described in the paper, we will produce the spectra for the QM- and MM-sampled moieties:

#### QM-sampled moiety; LYR-472

Execute one after the other, the following scripts:

1. [cp3_core_QM_7480_LYR-472_create_inputs.py](./cp3_core_QM_7480_LYR-472_create_inputs.py)

2. [cp3_core_QM_7480_LYR-472_run_embedding.py](./cp3_core_QM_7480_LYR-472_run_embedding.py)

#### MM-sampled moiety; LYR-225

Execute one after the other, the following scripts:


1. [cp3_core_MM_3523_LYR-225_create_inputs.py](./cp3_core_MM_3523_LYR-225_create_inputs.py)

2. [cp3_core_MM_3523_LYR-225_run_embedding.py](./cp3_core_MM_3523_LYR-225_run_embedding.py)


## Placing ECPs:

This script modifies the mol file and places ECPs in the surrounding MM atoms (6 Ang from core): [modify_mol.py](./modify_mol.py)


## OPA, TPA, 3PA dal files:

[opa.dal](./opa.dal)

[tpa.dal](./tpa.dal)

[3pa.dal](./3pa.dal)

