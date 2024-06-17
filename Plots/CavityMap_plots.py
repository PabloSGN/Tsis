# ======================================================================== #
# Imports

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np 
from scipy.special import erf

from scipy.integrate import simps, nquad

import etalon_funtions as etf
import functions_matricial as ft
import functions_tele as ft_tele
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

def PSF(wvls, l0, angle1, angle2, xi, eta, Da, Et):

    """
    Function that calculates the Etalon transmission profile.

    Parameters
    ----------
    angle1 : float
        Projection of the angle in X axis
    angle2 : float
        Projection of the angle in Y axis
    xi : float
        Coordinate 1 in detector plane in meters
    eta : float
        Coordinate 2 in detector plane in meters
    Da : float
        Relative thickness variation of etalon
    Et : dict
        Etalon Properties.

    Returns
    -------
    Inum : array
        Etalon Profile.

    """

    # Configuration
    
    Nl = int(len(wvls))
    Inum = np.zeros(Nl)
    tau  = etf.transm(Et['R'])
    
    xi0  = np.sin(np.deg2rad(angle1) * Et['fnum'])
    eta0 = np.sin(np.deg2rad(angle2) * Et['fnum'])
    
    xi  += xi0
    eta += eta0
    
    theta3 = 0
    f = Et['fnum'] * 2 
    lims = [[0, 1], [0, 2 * np.pi]]  # Integral limits
    accur = {'limit': 50, 'epsabs':0,'epsrel': 1.49e-8}  # Optional arguments
    j = -1
    F = etf.F(Et['R'])

    # Loop over wavelengths. 
    for wvli in wvls:
        
        j += 1
        k = 2 * np.pi / wvli
        
        params = (tau, Et['R'], F, xi, eta, xi0, eta0, k, f, wvli, theta3, 
                  Et['n'], Et['fnum'], Et['d'], Da, l0)
        H11tr = nquad(etf.H11pintr_, lims, args=params, opts=accur)
        H11ti = nquad(etf.H11pinti_, lims, args=params, opts=accur)
        H11t = (H11tr[0] + 1j * H11ti[0]) / (np.pi * 1 ** 2)
        Inum[j] = np.real(H11t * np.conj(H11t))

    return Inum


wavelengths = np.linspace(6173, 6174, 1000)
l0 = 6173.5

collimated = etalon_funct(wavelengths * 1E-10, l0 * 1E-10, Et, 0)
tele_perfe = ft_tele.Psi_single(wavelengths * 1E-10, 1, l0*1E-10, Et)

tele_imperfe = PSF(wavelengths * 1E-10, l0 * 1E-10, 0, 0.3, 0, 0, 1, Et)

axs.plot(wavelengths - l0, collimated, c = 'crimson', lw = 3, label = "Collimated")
axs.plot(wavelengths - l0, tele_perfe, c = 'orange', lw = 3, label = "Perfect telecentric")
axs.plot(wavelengths - l0, tele_imperfe, c = 'indigo', lw = 3, label = "Imperfect telecentric")

axs.legend(edgecolor = "k")
axs.set_xlim(-0.25, 0.25)
axs.set_ylabel("Amplitude")
axs.set_xlabel(r"$\lambda$ - $\lambda _ 0$ [$\AA$]")
axs.grid(True, c = 'k', alpha = 0.3)
plt.tight_layout()
plt.savefig(f"Plots/plots/etalon_setups_profiles.pdf", bbox_inches = "tight")

