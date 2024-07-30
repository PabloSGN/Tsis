#%%
# ============================ IMPORTS ====================================== #
 
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib as mp
from scipy.signal import savgol_filter

from astropy.io import fits

from scipy.integrate import simps

from scipy.interpolate import interp1d
from scipy import signal

from scipy.fft import fft2, ifft2, fftshift, fft, ifft, ifftshift


import functions_old as ft
 
# ============================= CONFIG ====================================== #
 
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

GainFile = 'Plots/InputGain.fits'

N = 100

wavelength = 6173.34388

Et = {
      'R' : 0.892,
      'n' : 2.3268,
      'd' : 281e-6,
      'theta' : 0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 60,
      'angle1' : 0.25,
      'angle2' : 0.45
      }
#%%
# ======================= GAIN INPUT AND Da ================================= #

MaxAngle = 0.00001
AngleMap = np.zeros((N, N))
MaxNorm = np.sqrt(2) * (N/2)
for i in range(N):
    for j in range(N):
        Norm = np.sqrt(abs(i - N/2) ** 2 + abs(j - N/2) ** 2)
        AngleMap[i, j] = (MaxAngle / MaxNorm) * Norm 

AngleMap -= 0.000005

Da = AngleMap * wavelength

fig, axs = plt.subplots(1, 2, figsize = (10.5, 5))

axs[1].set_title(r'$\Delta a$ [$\AA$]')
im = axs[1].imshow(Da, cmap = 'Spectral_r')

divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="5%", pad=0.05)

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

plt.colorbar(im, cax =cax)#, format=ticker.FuncFormatter(fmt))

Gain = fits.getdata(GainFile)[300:400, 300:400]
axs[0].set_title('Gain')
im = axs[0].imshow(Gain, cmap = 'gray')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax =cax)

plt.tight_layout()
plt.savefig("Plots/plots/Gain_Da_Inputs.pdf", bbox_inches='tight')

#%%
# =========================================================================== #


"""


lw = 6172
hw = 6174

Nwls = 5

Spectrum , Wavelengths = ft.fts_spectra(lw, hw)

wls = list(np.linspace(-105, 105, Nwls))
wls.append(300)
wls = np.array(wls) * 1E-3 + wavelength

ind = np.argmin(abs(Wavelengths - wls[-1]))
Spectrum = Spectrum / Spectrum[ind]





Colors = ['crimson', 'royalblue', 'deeppink', 'forestgreen', 'coral', 'steelblue']

plt.figure(figsize = (5, 4.85))
plt.plot(Wavelengths - wavelength, Spectrum, 'k', lw = 3, label = 'Solar Spectrum')
plt.xlim([-0.2, 0.4])

plt.grid(True, c = 'k', alpha = 0.2)

Prof = ft.Prof_col(Spectrum, wls, Wavelengths, 1, 1, Et, -1)

zorders = [1, 2, 5, 2, 1, 0]

for ind, wl in enumerate(wls):
    
    Etalon = ft.Psi_col(Wavelengths * 1E-10, 1, wl * 1E-10, Et)
    
    
    if ind == 0:
        plt.plot(Wavelengths - wavelength, savgol_filter(Etalon * Spectrum, 15, 3), zorder = zorders[ind],
                 color = Colors[ind], lw = 2, label = r'$\Psi ^ {\lambda _ s}(\lambda) \times O(\lambda)$')
        plt.scatter(wls[ind] - wavelength, Prof[ind], marker = 'P', 
                    color = Colors[ind], s = 150, zorder = 10, label = 'Measured Intensity')
    else:
        plt.plot(Wavelengths - wavelength, savgol_filter(Etalon * Spectrum, 15, 3),  zorder = zorders[ind],
                 color = Colors[ind], lw = 2)
        plt.scatter(wls[ind] - wavelength, Prof[ind], marker = 'P', 
                    color = Colors[ind], s = 150, zorder = 10)
        
plt.plot(wls - wavelength, Prof, ls = '--', lw = 1, color = 'navy', zorder = 0, label = 'Measured profile')
plt.xlabel(r'$\lambda - \lambda _ 0$ [$\AA$]')
plt.ylabel('Intensity Normalized')
plt.legend(edgecolor = 'k')
plt.tight_layout()"""
#plt.savefig("Plots_standard/ProfileMeasurement.pdf", bbox_inches ='tight')


#%% 

lw = 6172
hw = 6175
wavelength = 6173.34388
Nwls = 11

Spectrum , Wavelengths = ft.fts_spectra(lw, hw)
Spectrum = Spectrum / np.max(Spectrum)

wls = list(np.linspace(-200, 200, 11))
wls = np.array(wls) * 1E-3 + wavelength


Prof = ft.Prof_col(Spectrum, wls, Wavelengths, 1, 1, Et, -1)

# Extend me
# Extend measured profile to avoid boundary effects
wls_convolution  = np.zeros((len(wls) + 4))
Prof_convolution = np.zeros((len(wls) + 4))

wls_convolution[0]  = wls[0] - 0.8
wls_convolution[1]  = wls[0] - 0.4
wls_convolution[-2]  = wls[-1] + 0.4
wls_convolution[-1]  = wls[-1] + 0.8
Prof_convolution[0]  = Prof[0]
Prof_convolution[1]  = Prof[0]
Prof_convolution[-2]  = Prof[-1]
Prof_convolution[-1] = Prof[-1]

for ind, val in enumerate(Prof):
    wls_convolution[ind + 2] = wls[ind]
    Prof_convolution[ind + 2] = Prof[ind]
    




plt.figure(figsize = (10.5, 5))
plt.plot(Wavelengths - wavelength, Spectrum, 'indigo', lw = 3, label = 'Solar Spectrum')
# plt.scatter(wls - wavelength, Prof, marker = 'X', color = 'navy', zorder = 0, label = 'Measured profile')
plt.xlim([-0.95, 0.95])

plt.scatter(wls_convolution - wavelength, Prof_convolution, marker = 'x', 
            color = 'dodgerblue', zorder = 0, s = 150, label = 'Mean profile')

f2 = interp1d(wls_convolution, Prof_convolution, kind='cubic')

new_wls = np.linspace(wls_convolution[0], wls_convolution[-1], 1000)

plt.plot(new_wls - wavelength, f2(new_wls), ls = '--', color = 'dodgerblue', label = 'Mean profile interpolation')

Etalon = ft.Psi_col(new_wls * 1E-10, 
                    1, wavelength * 1E-10, Et)

plt.plot(new_wls- wavelength, Etalon, lw = 2, ls = '-', color = 'crimson', label = 'Etalon transmission profile')


def Deconvolve(Signal, filter, eps, N):

    """ WIENER-HELSTORM filter """
    
    H = fftshift(fft(fftshift(filter)))
    S = fftshift(fft(fftshift(Signal)))

    H_star = np.conj(H)
    H_ps = abs(H * H_star)

    W = H_star / (H_ps + eps)
    G = W * S

    Dec = np.real(ifftshift(ifft(ifftshift(G)))[N : 2* N])

    Dec /= np.max(Dec)

    return Dec


H = fftshift(fft(fftshift(Etalon)))
S = fftshift(fft(fftshift(f2(new_wls))))

Eps = 10
H_star = np.conj(H)
H_ps = abs(H * H_star)

W = H_star / (H_ps + Eps)
G = W * S
Dec = np.real(ifftshift(ifft(ifftshift(G))))

# plt.plot(Wavelengths, Spectrum, color = 'forestgreen', label = 'Spectrum')
# plt.plot(new_wls, f2(new_wls), c = 'crimson', label = 'Convolution')
plt.plot(new_wls - wavelength, Dec / np.max(Dec), color = 'orange', label = 'Deconvolution', lw = 3)
plt.grid(True, color = 'k', alpha = 0.2)
plt.xlim([-0.3, 0.3])
plt.xlabel(r'$\lambda - \lambda _ 0$ [$\AA$]')
plt.ylabel('Intensity Normalized')
plt.legend(edgecolor = 'k', loc = 'lower right')

Proffff = ft.Prof_col(Dec / np.max(Dec), wls, new_wls, 1, 1, Et, -1)

plt.tight_layout()
plt.savefig("Plots/plots/Deconvolution.pdf", bbox_inches='tight')


#%%


"""Wavelengths = np.linspace(6172.5, 6173.5, 10000)

Collimated = ft.Psi_col(Wavelengths * 1E-10, 1, 6173*1E-10, Et)
Telecentric = ft.Psi_tel(Wavelengths * 1E-10, 1, 6173*1E-10, Et)
Numerical = ft.Numerical_TelecentricEtalon(Wavelengths * 1E-10, 6173*1E-10, 1, Et)



plt.figure(figsize = (5, 4.85))
plt.plot(Wavelengths -  6173, Collimated, c = 'crimson', label = 'Collimated', lw = 3)
plt.plot(Wavelengths -  6173, Telecentric, c = 'limegreen', label = 'Telecentric Perfect', lw = 3, alpha = 0.8)
plt.plot(Wavelengths -  6173, Numerical, c = 'royalblue', label = 'Telecentric Imperfect', lw = 3)


plt.xlim(-0.2, 0.2)
plt.xlabel(r'$\lambda - \lambda _ 0$ [$\AA$]')
plt.ylabel('Amplitude')
plt.grid(True, alpha = 0.3, c = 'k')
plt.legend(edgecolor = 'k')

plt.tight_layout()"""
#plt.savefig("Plots_standard/Profiles.pdf")


#%%

# DERIVATIVES

wavelength = 6173.34388
# Wavelengths = np.linspace(6172.5, 6173.5, 10000)
wls = list(np.linspace(-500, 500, 500))
wls = np.array(wls) * 1E-3 + wavelength

Spectrum, Wavelengths = ft.fts_spectra(6172, 6175)

Et['R'] = 0.892


perfil = ft.Prof_col(Spectrum, wls, Wavelengths, 0.95, 1, Et, -1)


def xi(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont, R):
    
    Et['R'] = R
    
    dentro = ft.fun_col(phi_obs, scan_wls, etalon_wls, g, da, phi_real, Et, Ncont)

    return dentro ** 2


R_vec = np.linspace(0.85, 0.99, 500)
h = 0.003

   
xi_plus  = xi(perfil, wls, Wavelengths, 0.81, 1, Spectrum, Et, -1, 0.92 + h) 
xi_minus = xi(perfil, wls, Wavelengths, 0.81, 1, Spectrum, Et, -1,  0.92 - h) 
    
DxiR = (xi_plus - xi_minus) / (2*h)
    


h = 0.003
xi_plus = xi(perfil, wls, Wavelengths, 0.7 + h, 1, Spectrum, Et, -1, 0.892) 
xi_minus = xi(perfil, wls, Wavelengths, 0.7 - h, 1, Spectrum, Et, -1,  0.892) 

Dxig=(xi_plus - xi_minus) / (2*h)


h = 0.000005

    
xi_plus = xi(perfil, wls, Wavelengths, 0.95, 0.999995 + h, Spectrum, Et, -1, 0.892) 
xi_minus = xi(perfil, wls, Wavelengths, 0.95, 0.999995 - h, Spectrum, Et, -1,  0.892) 

DxiDa = (xi_plus - xi_minus) / (2*h)
    


plt.figure(figsize = (6, 5.85))
# plt.plot(DxiDa /np.max(DxiDa))
plt.plot(wls - wavelength, Dxig, c = 'indigo', lw = 2,label = r'$\partial \varepsilon / \partial g$')
plt.plot(wls - wavelength, DxiR, c = 'orange', lw = 2,label = r'$\partial \varepsilon / \partial R$')
plt.plot(wls - wavelength, savgol_filter(DxiDa / abs(np.min(DxiDa)), 11, 3), c = 'crimson', lw = 2, label = r'$\partial \varepsilon / \partial \Delta a$  [Norm]')

plt.xlabel(r'$\lambda - \lambda _ 0$ [$\AA$]')
plt.grid(True, c = 'k', alpha = 0.2)
plt.legend(edgecolor = 'k', loc = 'upper right')

plt.ylabel('Amplitude [Arbitary Units]')


plt.xlim(-0.4, 0.4)

plt.tight_layout()
plt.savefig("Plots/plots/Derivadas.pdf", bbox_inches='tight')


