# ======================================================================== #
# Imports

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np 
from scipy.special import erf

from scipy.integrate import simps, nquad

from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


import etalon_funtions as etf
import functions_matricial as ft
import functions_tele as ft_tele
#import etalon as et
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

mosaic = """AAB"""

fig, axs = plt.subplot_mosaic(mosaic, figsize = (15, 5))

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


"""etalon_model = et.Etalon(checkpoint='Plots/2024-01-10-09_25_54.best.pth', gpu=0, verbose=True)

l0 = 6173.5 
wavelengths = np.arange(-500, 505, 1.5) * 1E-3 + l0

collimated = etalon_funct(wavelengths * 1E-10, l0 * 1E-10, Et, 0)
tele_perfe = ft_tele.Psi_single(wavelengths * 1E-10, 1, l0*1E-10, Et)

npix = 50
n_models = npix * npix

s = 300e-6
angle2 = np.ones(n_models) * 0
xi, eta = np.meshgrid(np.linspace(-s, s, npix), np.linspace(-s, s, npix))
Da  =    np.ones(n_models) * 1

et_asym_mod = etalon_model.evaluate(np.ones(n_models) * 0.3, np.ones(n_models) * 0.0, xi.flatten(), eta.flatten(), Da, wavelengths)

tele_imperfe = et_asym_mod.reshape(npix, npix , len(wavelengths))# PSF(wavelengths * 1E-10, l0 * 1E-10, 0, 0.3, 0, 0, 1, Et)
print('k?')
axs["A"].plot(wavelengths - l0, collimated, c = 'crimson', lw = 3, label = "Collimated")
axs["A"].plot(wavelengths - l0, tele_perfe, c = 'orange', lw = 3, label = "Perfect telecentric")
axs["A"].plot(wavelengths - l0,  savgol_filter(tele_imperfe[24, 24], 15, 3) , c = 'indigo', lw = 3, label = "Imperfect telecentric")

axs["A"].legend(edgecolor = "k")
axs["A"].set_xlim(-0.25, 0.25)
axs["A"].set_ylabel("Amplitude")
axs["A"].set_xlabel(r"$\lambda$ - $\lambda _ 0$ [$\AA$]")
axs["A"].grid(True, c = 'k', alpha = 0.3)

im = axs["B"].imshow(tele_imperfe[:, :, len(wavelengths) // 2], norm = 'log', cmap = 'inferno')
plt.colorbar(im, fraction = 0.046, pad = 0.04, label = "Amplitude")

axs["B"].set_xticks([])
axs["B"].set_yticks([])

axs["A"].set_title("Spectral PSFs")
axs["B"].set_title(r"Imperfect telecentric spatial PSF at $\lambda _ 0$")

plt.tight_layout()"""
#plt.savefig(f"Plots/plots/etalon_setups_profiles.pdf", bbox_inches = "tight")


# ============================================================ #

#plt.close("all")

fig, axs = plt.subplots(figsize = (10.5, 5))

Et = {
      'R' : 0.925,
      'n' : 2.2926446,
      'd' : 251.63e-6,
      'theta' : 0.0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 56.6,
      'angle' : 0
      }

Spectrum, Wavelengths = ft.fts_spectra(6172.5, 6175)
Spectrum = Spectrum / np.max(Spectrum)
Wavelengths_increased = np.arange(6172.5, 6174.3, 0.0005)
spc_interp = interp1d(Wavelengths, Spectrum, kind = 'cubic')
spc_increased = spc_interp(Wavelengths_increased)
Wavelengths = Wavelengths_increased
Spectrum = spc_increased

l0 = 6173.344

axs.plot(Wavelengths - l0, Spectrum, c = 'k', lw = '3', label = 'Solar spectrum')

wls = list(np.linspace(-120, 120, 5))
wls.append(300)
wls = np.array(wls) * 1E-3 + l0

colors = ["crimson", "darkorange", "deeppink", "indigo", "dodgerblue", "forestgreen"]

profile = []

for ind, wl in enumerate(wls):
    etalon = etalon_funct(Wavelengths * 1E-10, wl * 1E-10, Et, 0)
    axs.plot(Wavelengths - l0, etalon * Spectrum, lw = 2, c = colors[ind])
    
    profile.append(np.sum(etalon * Spectrum))

profile = np.array(profile) / profile[-1]

for ind, wl in enumerate(wls):
    axs.scatter(wl - l0, profile[ind], s = 120, marker = 'P', c = colors[ind], zorder = 10)


axs.scatter(10, 0.5, s = 75, marker = 'P', c = 'k', label = 'Measured Intensity')
axs.plot([10, 11], [0.5, 0.5], c = 'k', lw = 2)


axs.plot(wls - l0, profile, marker ='x', ls = "--", color = 'indigo', zorder = 0)
axs.set_xlim(-0.2, 0.4)
axs.set_ylabel("Intensity [Norm.]")
axs.set_xlabel(r"$\lambda - \lambda _ 0$ [$\AA$]")
axs.grid(True, c = 'k', alpha = 0.3)

# Multicolored legend
from matplotlib.collections import PatchCollection
class MulticolorPatch(object):
    def __init__(self, colors):
        self.colors = colors
        
# define a handler for the MulticolorPatch object
class MulticolorPatchHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        width, height = handlebox.width, handlebox.height
        height *=0.35
        patches = []
        for i, c in enumerate(orig_handle.colors):
            patches.append(plt.Rectangle([width/len(orig_handle.colors) * i - handlebox.xdescent, 
                                          -handlebox.ydescent],
                           width / len(orig_handle.colors),
                           height, 
                           facecolor=c, 
                           edgecolor='none'))

        patch = PatchCollection(patches,match_original=True)

        handlebox.add_artist(patch)
        return patch


h, l = axs.get_legend_handles_labels()
h.append(MulticolorPatch(colors))
l.append(r"$\psi ^{\lambda _ s} (\lambda)\times\mathcal{O}(\lambda)$")


axs.legend(h, l, loc='lower right', edgecolor = 'k' ,
         handler_map={MulticolorPatch: MulticolorPatchHandler()})
plt.tight_layout()


plt.savefig(f"Plots/plots/ProfileMeasurement.pdf", bbox_inches = "tight")
plt.savefig(f"figures/EtalonPaper/ProfileMeasurement.pdf", bbox_inches = "tight")