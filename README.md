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

### $\emph{mfcc}$