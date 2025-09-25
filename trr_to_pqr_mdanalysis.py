#

import MDAnalysis as mda
import sys

u = mda.Universe('../../mimic.tpr', 'mimic_R1-R8_40ps_settime_centered.xtc')
print (u)

# spacing = 2 means write trajectory every 2 steps, i.e.
# similar to: u.trajectory[0::2]
frame_number = 0
# spacing = 400
spacing = 40
for frame in u.trajectory[0::spacing]:
    print (frame, frame_number)
    u.atoms.write(f"R1-R8_40ps_settime_centered_frame_{frame_number}.pqr",\
                  frames = u.trajectory[frame_number:frame_number+1])
    frame_number = frame_number + spacing


