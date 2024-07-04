# ======================================================================== #
# Imports

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np 
from scipy.special import erf
import functions_matricial as ft

from astropy.io import fits
import matplotlib.patches as patches

import etalon_mancha as et

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

etalon = et.Etalon(checkpoint='Plots/Inputs/2024-01-10-09_25_54.best.pth', gpu=0, verbose=True)
l0 = 6173.5 
wvls = np.arange(-500, 505, 5) * 1E-3 + l0


# ======================================================================== #

mosaic = """ABC
            DDE"""


fig, axs = plt.subplot_mosaic(mosaic, figsize = (12, 7))

# Symetric vs Asymetric

npix = 50
n_models = npix * npix

s = 300e-6
angle2 = np.ones(n_models) * 0
xi, eta = np.meshgrid(np.linspace(-s, s, npix), np.linspace(-s, s, npix))
Da  =    np.ones(n_models) * 1

et_sym_mod = etalon.evaluate(np.ones(n_models) * 0, angle2, xi.flatten(), eta.flatten(), Da, wvls)
et_asym_mod = etalon.evaluate(np.ones(n_models) * 0.5, np.ones(n_models) * 0, xi.flatten(), eta.flatten(), Da, wvls)
et_sym = et_sym_mod.reshape(npix, npix , len(wvls))
et_asym = et_asym_mod.reshape(npix, npix , len(wvls))

axs["E"].set_title("Transmission profiles")
axs["E"].plot(wvls- l0, et_sym[24, 24], c = 'indigo', lw = 3, label = 'Symmetric')
axs["E"].plot(wvls - l0, et_asym[24, 24], c = 'crimson', lw = 3, label = 'Asymmetric')
axs["E"].legend(edgecolor = 'k')
axs["E"].set_ylabel("Amplitude")
axs["E"].set_xlabel(r"$\lambda$ - $\lambda _ 0$ [$\AA$]")
axs["E"].grid(True, c = 'k', alpha = 0.3)
axs["E"].set_xlim(-0.2, 0.2)
axs["E"].set_ylabel("Amplitude [Norm.]")
# Simulation Inputs

# ONLY LOCAL PC IAA
data_file =         "/home/pablo/Desktop/EtalonCorrection/ObsSimulation/data.fits"
flat = fits.getdata("/home/pablo/Desktop/EtalonCorrection/ObsSimulation/flats.fits")
cav = fits.getdata( "/home/pablo/Desktop/EtalonCorrection/ObsSimulation/CavityMap_Amstng.fits")

l0 = 6173.5 

gain = flat[0, 500:500+561, 500:500+561]
cav  = cav[300:300+561, 700:700+561]

cav -= np.mean(cav)
cav *= 2

header = fits.getheader(data_file)
data = fits.getdata(data_file)


im = axs["A"].imshow(data[:, :, 1, 100], cmap = 'hot')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs["B"].imshow(cav, cmap = 'bwr')
plt.colorbar(im, fraction=0.046, pad=0.04, label = r"$\Delta \lambda$ [$\AA$]")

im = axs["C"].imshow(gain, cmap = 'gist_heat')
plt.colorbar(im, fraction=0.046, pad=0.04)

axs["A"].set_title("Simulation")
axs["B"].set_title("Cavity Map")
axs["C"].set_title("Gain")

for i in ["A", "B", "C"]:
    axs[i].set_xticks([])
    axs[i].set_yticks([])

sqx = [40, 100]
sqy = [40, 100]
axs["A"].plot([sqx[0], sqx[1]], [sqy[0], sqy[0]], ls = '-', c = 'k', lw = 2)
axs["A"].plot([sqx[0], sqx[1]], [sqy[1], sqy[1]], ls = '-', c = 'k', lw = 2)
axs["A"].plot([sqx[0], sqx[0]], [sqy[0], sqy[1]], ls = '-', c = 'k', lw = 2)
axs["A"].plot([sqx[1], sqx[1]], [sqy[0], sqy[1]], ls = '-', c = 'k', lw = 2)
rect = patches.Rectangle((40, 40), 60, 60, linewidth=1, edgecolor='k', facecolor='cyan', alpha = 0.5, label = "Averaged Area")
axs["A"].add_patch(rect)
axs["A"].legend(edgecolor = 'k', loc = 'lower right')





def psf_shift(l0, psf, wvls):

    fixed_l0 = 6173.5 

    new_psf = np.zeros(np.shape(psf))

    up   = np.max(wvls)
    down = np.min(wvls)

    shift = l0 - fixed_l0

    new_axis = wvls + shift

    if shift < 0:
        lower_cut = np.argmin(abs(new_axis - down))

        lenprof = round(len(wvls) - lower_cut)
        
        new_psf[:, :, :lenprof] = psf[: , :, lower_cut:(lenprof + lower_cut)]
        
        for i in range(lenprof, len(wvls)):
            new_psf[:, :, i] = psf[:, :, -1]

    elif shift > 0:

        upper_cut = np.argmin(abs(new_axis - up))

        lenprof = round(len(wvls) - upper_cut)

        for i in range(lenprof):
            new_psf[:, :, i] = psf[:, :, 0]

        new_psf[:, :, lenprof - 1:] = psf[: , :, :upper_cut + 1]
    
    else:
        new_psf = psf
        
    return new_psf


axs["D"].plot(np.arange(-500, 505, 10) * 1E-3 + l0,  np.mean(data[40:100, 40:100, 0, :], axis = (0, 1)) / np.max(np.mean(data[40:100, 40:100, 0, :], axis = (0, 1))), zorder = 2, lw = 3, c = 'indigo', label = 'Average profile') 


colors = ['darkorange', 'crimson', 'crimson', 'cyan', 'orange', 'mediumspringgreen', 'm', 'royalblue', 'forestgreen']
wvls_scan = np.arange(6173.5 - 0.2, 6173.5 + 0.25, 0.01)
ic = 0
labels = ["Lower limit", "Upper limit"]
for i, wvl in enumerate(wvls_scan):
    if i == 0:
        axs["D"].plot([wvl, wvl], [0, 1], lw = 1, ls = '--', c = 'k', alpha = 0.5, label = "Scanned positions")
    
    axs["D"].plot([wvl, wvl], [0, 1], ls = '--', c = 'k', alpha = 0.5)

    if i == 0 or i == len(wvls_scan) -1:
        psf_shifted = psf_shift(wvl, et_sym, wvls)
        axs["D"].plot(wvls, psf_shifted[25, 25, :] / np.max(psf_shifted[25, 25, :]), lw = 3, alpha = 1, c = colors[ic], label = labels[ic])
        ic += 1
        if ic == len(colors):
            ic = 0

axs["D"].grid(True, c = 'k', alpha = 0.3)
axs["D"].legend(edgecolor = 'k', loc = "lower left")
axs["D"].set_xlabel(r"Wavelengths [$\AA$]")


plt.tight_layout()

plt.savefig("Plots/plots/Inputs_mancha.pdf", bbox_inches = "tight")
plt.savefig("Plots/plots/Inputs_mancha.png", bbox_inches = "tight")