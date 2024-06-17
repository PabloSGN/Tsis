# ======================================================================== #
# Imports

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np 
from scipy.special import erf

import functions_matricial as ft
import functions_tele as ft_tele
import functions as ft_imperf
# ======================================================================== #
# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

Et = {
      'R' : 0.925,
      'n' : 2.2926446,
      'd' : 251.63e-6,
      'theta' : 0.0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 56.6,
      'angle' : 0.23
      }
# ======================================================================== #
# Etalon Profiles Plot
fig, axs = plt.subplots(figsize = (10.5, 5))

def etalon_funct(wvls, l0, Et, dn):
    m  = round(2 * Et['n'] * Et['d'] / l0)    
    dh = (m * l0 - 2 * Et['n'] * Et['d']) / (2 * Et['n'])     
    thick = Et['d'] + dh

    a =  thick * (Et['n'] + dn) * Et['cos_th']

    return 1 / (1 + ft.Finesse(Et['R']) * np.sin(2 * np.pi * a / wvls) ** 2)

wavelengths = np.linspace(6173, 6174, 1000)
l0 = 6173.5

collimated = etalon_funct(wavelengths * 1E-10, l0 * 1E-10, Et, 0)
tele_perfe = ft_tele.Psi_single(wavelengths * 1E-10, 1, l0*1E-10, Et)

da_ref = ft_imperf.da_tuning(Et, l0 * 1E-10)
tele_imperfe = ft_imperf.Psi(wavelengths * 1E-10, 6173.67 * 1E-10, 1, da_ref, Et, 0.3)

axs.plot(wavelengths - l0, collimated, c = 'crimson', lw = 3)
axs.plot(wavelengths - l0, tele_perfe, c = 'orange', lw = 3)
axs.plot(wavelengths - l0, tele_imperfe, c = 'indigo', lw = 3)

axs.set_xlim(-0.5, 0.5)
axs.set_ylabel("Amplitude")
axs.set_xlabel(r"$\lambda$ - $\lambda _ 0$ [$\AA$]")
axs.grid(True, c = 'k', alpha = 0.3)
plt.tight_layout()
plt.show()

