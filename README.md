# ChR2 paper

<div align="center">
  <img src="./TOC-ChR2-multiphoton.png" width="100px" alt="TOC" />
  <p><strong>Multiphoton Absorption Spectra of Channelrhodopsin-2 via Multiscale Simulation Methods</strong></p>
</div>

This repository hosts the scripts and files used to reproduce the data presented in https://doi.org/10.26434/chemrxiv-2025-vzz1r

## Table of Contents

- [Split trajectory into pqr files](#split-trajectory-into-pqr-files)
- [Creation of potentials](#creation-of-potentials)
    - [mfcc approach](#mfcc-approach)
    - [cp3 approach](#cp3-approach)
        - [QM-sampled moiety; LYR-472](#qm-sampled-moiety-lyr-472)
        - [MM-sampled moiety; LYR-225](#mm-sampled-moiety-lyr-225)
- [Placing ECPs](#placing-ecps)
- [OPA, TPA, 3PA dal files](#opa-tpa-3pa-dal-files)


## Split trajectory into pqr files:

```
trr_to_pqr_snapshots.py
```

## Creation of potentials

The creation of the `.pot` files has been generated through two approaches (see publication for full details):

- The __Molecular Fractionation with Conjugate Caps__ approach, i.e. "*mfcc*"
- __Standard potenitals__,  i.e. . "*cp3*" 


### mfcc approach

Execute one after the other, the following scripts:

1. `mfcc_ChR2_and_close-water_create_inputs.py` creates the input fragment files

2. `mfcc_ChR2_and_close-water_write_potential.py` writes the embedding potential in the form of a .pot file


### cp3 approach

This approach was used for production, where, as described in the paper, we will produce the spectra for the QM- and MM-sampled moieties:

#### QM-sampled moiety; LYR-472

Execute one after the other, the following scripts:

1. `cp3_core_QM_7480_LYR-472_create_inputs.py`

2. `cp3_core_QM_7480_LYR-472_run_embedding.py`

#### MM-sampled moiety; LYR-225

Execute one after the other, the following scripts:


1. `cp3_core_MM_3523_LYR-225_create_inputs.py`

2. `cp3_core_MM_3523_LYR-225_run_embedding.py`


## Placing ECPs:

This script modifies the mol file and places ECPs in the surrounding MM atoms (6 Ang from core):

```
modify_mol.py
```

## OPA, TPA, 3PA dal files:

```
opa.dal
tpa.dal
3pa.dal
```
