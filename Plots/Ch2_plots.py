# ======================================================================== #
# Imports

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np 
from scipy.special import erf
import functions_matricial as ft

# ======================================================================== #
# CONFIG

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

Et = {
      'R' : 0.925,
      'n' : 2.2926446,
      'd' : 251.63e-6,
      'theta' : 0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 56.6,
      'angle' : 0.23
      }
# ======================================================================== #
# FIGURE -> Example of Etalon Transmission Profile + Prefilter

l0 = 6173
Wavelengths = np.linspace(6167, 6179, 1500)

def etalon_funct(wvls, l0, Et, dn):
    m  = round(2 * Et['n'] * Et['d'] / l0)    
    dh = (m * l0 - 2 * Et['n'] * Et['d']) / (2 * Et['n'])     
    thick = Et['d'] + dh

    a =  thick * (Et['n'] + dn) * Et['cos_th']

    return 1 / (1 + ft.Finesse(Et['R']) * np.sin(2 * np.pi * a / wvls) ** 2)

def smoothed_top_hat(x, center, width, smoothing):
    """
    Creates a smoothed top-hat function using error functions.

    Parameters:
    x (numpy array): The input array over which to evaluate the function.
    center (float): The center of the top-hat function.
    width (float): The width of the top-hat function.
    smoothing (float): The smoothing parameter, determining how smooth the transitions are.

    Returns:
    numpy array: The values of the smoothed top-hat function.
    """
    left = center - width / 2
    right = center + width / 2
    
    # Use erf to smooth the edges
    smoothed_left = 0.5 * (1 + erf((x - left) / smoothing))
    smoothed_right = 0.5 * (1 - erf((x - right) / smoothing))
    
    return smoothed_left * smoothed_right


etalon = etalon_funct(Wavelengths * 1E-10, l0 * 1E-10, Et, 0)
etalon_shifted = etalon_funct(Wavelengths * 1E-10, l0 * 1E-10, Et, 0.085)
prefilter = smoothed_top_hat(Wavelengths, l0, 2.5, 0.75)

fig, axs = plt.subplots(figsize = (10.5, 5))
axs.plot(Wavelengths - l0, etalon, c = 'crimson', lw = 1.5, ls = '--')
axs.plot(Wavelengths - l0, etalon * prefilter, c = 'crimson', lw = 2, ls = '-', label = f"n : {round(Et['n'], 2)}")

axs.plot(Wavelengths - l0, etalon_shifted, c = 'orange', lw = 1.5, ls = '--')
axs.plot(Wavelengths - l0, etalon_shifted * prefilter, c = 'orange', lw = 2, ls = '-',label = f"n : {round(Et['n'] + 0.085, 2)}")

axs.plot(Wavelengths - l0, prefilter, c = 'm', lw = 0.5)
axs.grid(True, c = 'k', alpha = 0.3)

axs.set_xlabel(r"$\lambda - \lambda _0\ \ [\AA]$")
axs.set_ylabel("Transmisivity")

axs.set_xlim(-6, 6)


axs.fill_between(Wavelengths - l0, np.zeros(len(Wavelengths)), prefilter, color = 'magenta', alpha = 0.1, label = 'Pre-filter')


axs2 = axs.twinx()

axs2.plot([10, 10], [0.5, 0.6], ls = '-', c = 'k', lw = 2, label = "Transmission profile")
axs2.plot([10, 10], [0.5, 0.6], ls = '--', c = 'k', lw = 1.5, label = "Filtered transmission")


axs.legend(loc = 'upper right', edgecolor = 'k')
axs2.legend(loc = 'upper left', edgecolor = 'k')

axs2.set_yticks([])

plt.tight_layout()

plt.savefig(f"Plots/plots/Etalon_and_prefilter_example.pdf", bbox_inches = 'tight')

plt.savefig(f"figures/Introduction_to_spectropolarimeters/Etalon_and_prefilter_example.pdf", bbox_inches = 'tight')
