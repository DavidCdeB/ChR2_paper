#

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import os


def get_sigma_tpa(fname):
    with open(fname) as f:
        for line in f:
            word = "Sym  No  Energy  Polarization"
            S = []
            if word in line:
                next(f)
                target=next(f)
                S1 = float(target.split()[7])
                S.append(S1)

                next(f)
                target=next(f)
                S2 = float(target.split()[7])
                S.append(S2)

                next(f)
                target=next(f)
                S3 = float(target.split()[7])
                S.append(S3)

                next(f)
                target=next(f)
                S4 = float(target.split()[7])
                S.append(S4)

                next(f)
                target=next(f)
                S5 = float(target.split()[7])
                S.append(S5)

                aux = np.array(S)
                return aux

# Autocorrelation function:
def autocorrelation_3(scalar_properties):
    """
    Calculate the autocorrelation function of scalar properties.

    :param scalar_properties: A list of scalar values at different time steps.
    :return: A list of autocorrelation values.
    """
    n = len(scalar_properties)
    mean = np.mean(scalar_properties)
    # Variance:
    var = np.var(scalar_properties)
    acf = np.correlate(scalar_properties - mean, scalar_properties - mean, mode='full') / (var * n)
    return acf[n-1:]

# Exponential fit:
def exp_decay(t, A, tau):
    """
    Exponential decay function.

    :param t: The time variable.
    :param A: Amplitude at t=0.
    :param tau: Characteristic decay time.
    :return: Value of the function at time t.
    """
    return A * np.exp(-t / tau)


# Last restart (R8) contains 10,000 frames in mimic output
# 10,000 frames * 0.5fs = 5000 fs = 5 ps

# We'll always use the gromacs mimic_centered.trr, that saves every 10 steps,
# or every 5 fs
# (Each frame separated by 10 frames*0.5 fs = 5 fs)
# so: 10,000/10 = 1,000 frames * 5 fs = 5,000 fs

# frame 0; time = 0 fs
# frame 1; time = 5 fs
# frame 2; time = 10 fs
# frame 3; time = 15 fs
# ...
# frame 1000; time = 5,000 fs

# Frames-domain:
first_frame = 0
last_frame = 1000
spacing = 1

# Times-domain (in fs):
first_time = 0
last_time = 5000
spacing_time = 5

folders = [i for i in range(first_frame, last_frame+1, spacing)]
times_in_fs = [i for i in range(first_time, last_time+1, spacing_time)]

# Gather TPA cross-sections:
all_S = []
tpa_out = "tpa.out"

for i in folders:
    fname = f"R8_mimic_centered_frame_{i}/OPA/{tpa_out}"
    S = get_sigma_tpa(fname)[0]
    all_S.append(S)

# Plot cross-sections:
plt.plot(np.array(times_in_fs),  np.array(all_S),\
             'o-', color='blue')

plt.ticklabel_format(useOffset=False)
plt.ylabel(r"Sigma TPA, $\sigma^{\mathrm{2PA}}$ (GM)", fontsize=18)
plt.xlabel(r'Time (fs)', fontsize=18)

plt.grid()
plt.savefig(f'sigma_tpa.pdf', bbox_inches='tight')
plt.savefig(f'sigma_tpa.png', bbox_inches='tight')


# Autocorrelation function (possibility 3 outlined in "correlation_sigmas.ipynb"):
c = autocorrelation_3(all_S)

# Fit the autocorrelation function to an exponential decay
# Note: times_in_fs should be the array of times corresponding to the autocorrelation values in c
params, covariance = curve_fit(exp_decay, np.array(times_in_fs), np.array(c))

# Extract fitted parameters
A_fit, tau_fit = params

# Generate fitted curve
t_fit = np.linspace(min(times_in_fs), max(times_in_fs), 1000)
c_fit = exp_decay(t_fit, A_fit, tau_fit)

# Plotting
fig, ax = plt.subplots()
ax.plot(times_in_fs, c, 'o-', color='black', markersize=2, linewidth=0.3, label='Data')
ax.plot(t_fit, c_fit, '-', color='red', label=f'Fit: A={A_fit:.2f}, $\\tau$={tau_fit:.2f}')

ax.ticklabel_format(useOffset=False)
ax.set_ylabel(r"$C(t)$", fontsize=18)
ax.set_xlabel(r'Time [fs]', fontsize=18)
# ax.grid(linewidth=0.3)
ax.grid(axis='x', color='0.80')
ax.grid(axis='y', color='0.80')
ax.legend(loc='upper right', fontsize=12)

ax.tick_params(axis='y', labelsize=14)
ax.tick_params(axis='x', labelsize=14)

# Inset
ax_inset = inset_axes(ax, width="40%", height="30%", loc='center right', bbox_to_anchor=(-0.2, 0.03, 1.2, 1.2), bbox_transform=ax.transAxes)
ax_inset.plot(times_in_fs, c, 'o-', markersize=2, linewidth=0.3, color='black')
ax_inset.plot(t_fit, c_fit, '-', color='red')
ax_inset.set_xlim([-20, 1200])
ax_inset.set_ylim([-0.186, 1.048])
ax_inset.set_xticks([0, 200, 400, 600, 800, 1000])  # Set specific xticks
ax_inset.grid(axis='x', color='0.80')
ax_inset.grid(axis='y', color='0.80')

# Save the figures
plt.savefig('C_ti_book_sigma_tpa_fit_with_inset.pdf', bbox_inches='tight')
plt.savefig('C_ti_book_sigma_tpa_fit_with_inset.png', bbox_inches='tight')

plt.show()


