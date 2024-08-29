import matplotlib.pyplot as plt
import glob
import matplotlib as mp
import numpy as np

from images_reader import image_reader

import functions_old as ft

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

sun_obs = "/media/pablo/T7 Touch/SUNRISE_3/OldTests/TuMAG/INTA_2021/12122021/Sun_scan_voltage_51727/Scan 3/"
pre_folder = "/media/pablo/T7 Touch/SUNRISE_3/OldTests/TuMAG/INTA_2021/02122021/Spectroscopic E2E/Iodine cell/without cell/517_30/"

all_ims = sorted(glob.glob(f"{sun_obs}*0.img"))
pref_ims = sorted(glob.glob(f"{pre_folder}*0.img"))

intensity = []
prefilter = []

for i, im in enumerate(all_ims):
    print(i)
    I = image_reader(im ,PlotFlag= False)
    intensity.append(np.mean(I[500:1500, 500:1500]))

for i, im in enumerate(pref_ims):
    print(i)
    I = image_reader(im ,PlotFlag= False)
    prefilter.append(np.mean(I[500:1500, 500:1500]))

Volts_vec = np.arange(-4000, 4000, 50)
prefilter = np.array(prefilter)

intensity = np.array(intensity)

prefilter /= np.max(prefilter)
intensity /= np.max(intensity)


Spectrum, Wavelengths = ft.fts_spectra(5171, 5174)
Spectrum = Spectrum / np.max(Spectrum)

fig, axs = plt.subplots(1, 2, figsize = (10.5, 5))


axs[0].plot(Volts_vec, prefilter, color = 'indigo', lw = 3, label = "Pre-filter measurement")
axs[0].plot(Volts_vec, intensity / prefilter, c = 'darkorange', lw = 3, label = "517 nm line observation.")

axs[0].legend(edgecolor = 'k', loc = 'upper right', ncol = 1)

axs[1].plot([-100, -100], [0, 0], color = 'crimson', lw = 3, label = r"Spectrum fit-$1^{st}$ order.")
axs[1].plot([-100, -100], [0, 0], color = 'seagreen', lw = 3, label = r"Spectrum fit-$2^{nd}$ order.")

axs[1].legend(edgecolor = 'k', loc = 'upper right', ncol = 1)



axs[1].plot(Volts_vec, prefilter, color = 'indigo', lw = 3)
axs[1].plot(Volts_vec, intensity / prefilter, c = 'darkorange', lw = 3)

axs[0].set_xlim(-4000, 0)
axs[1].set_xlim(0, 4000)



ax2 = axs[0].inset_axes([0, 0, 1, 1], transform=axs[0].transAxes, facecolor='none') # overlay ax2 on top of ax1
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="crimson")
ax2.tick_params(axis='y', colors="crimson")
ax2.set_yticks([])

ax2.set_xlabel(r"Wavelengths [$\AA$]", color = "crimson")

ax2.plot(Wavelengths, Spectrum, c = 'crimson', lw = 3 )
ax2.set_xlim(5172.17, 5173.35)
ax2.set_ylim(0, 0.9)

ax2 = axs[1].inset_axes([0, 0, 1, 1], transform=axs[1].transAxes, facecolor='none') # overlay ax2 on top of ax1
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="seagreen")
ax2.tick_params(axis='y', colors="seagreen")

ax2.set_xlabel(r"Wavelengths [$\AA$]", color = "seagreen")
ax2.plot(Wavelengths, Spectrum, c = 'seagreen', lw = 3 )
ax2.set_xlim(5171.55, 5172.785)
ax2.set_ylim(0.05, 0.9)

ax2.set_ylabel("Spectrum intensity [Norm. to contiuum]", color = "seagreen")
axs[0].set_xlabel("Etalon voltage [V]")
axs[1].set_xlabel("Etalon voltage [V]")

axs[0].set_ylabel("Measured intensity [Norm. to pre-filter]")

#axs[0].grid(True, c = 'k', alpha = 0.3)
#axs[1].grid(True, c = 'k', alpha = 0.3)


axs[0].set_ylim(0.14, 1.42)
axs[1].set_ylim(0.14, 1.42)

axs[1].set_yticks([])
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.savefig("plots/secondorder.pdf", bbox_inches = "tight")





