# ChR2 paper

This repository hosts the scripts and files used to reproduce the data presented in https://doi.org/10.26434/chemrxiv-2025-vzz1r

## MDAnalysis script to split trajectory into .pqr files:

```python
import MDAnalysis as mda

u = mda.Universe('./Data/mimic.tpr', 'Data/mimic_R1-R8_40ps_settime_centered.xtc')
print (u)

# spacing = 2 means write trajectory every 2 steps, i.e.
# similar to: u.trajectory[0::2]
frame_number = 0
spacing = 40
for frame in u.trajectory[0::spacing]:
    print (frame, frame_number)
    u.atoms.write(f"R1-R8_40ps_settime_centered_frame_{frame_number}.pqr",\
                  frames = u.trajectory[frame_number:frame_number+1])
    frame_number = frame_number + spacing
```

## Creation of potentials

The creation of the `.pot` files has been generated through two approaches (see publication for full details):

- The __Molecular Fractionation with Conjugate Caps__ approach, i.e. "*mfcc*"
- __Standard potenitals__,  i.e. . "*cp3*" 


### mfcc approach

```
mfcc_ChR2_and_close-water_create_inputs.py
mfcc_ChR2_and_close-water_write_potential.py
```

### cp3 approach

This approach was used for production, 

#### QM core = LYR-472

```
cp3_core_QM_7480_LYR-472_create_inputs.py
cp3_core_QM_7480_LYR-472_run_embedding.py
```

#### QM core = LYR-225

```
cp3_core_MM_3523_LYR-225_create_inputs.py
cp3_core_MM_3523_LYR-225_run_embedding.py
```