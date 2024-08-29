import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from utils import read_Tumag

from scipy.io import loadmat 
from scipy.signal import convolve, find_peaks, resample
from scipy.optimize import minimize


import functions_old as ft

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

Et = {
      'R' : 0.895,
      'n' : 2.2926446,
      'd' : 251.63e-6,
      'theta' : 0,
      'cos_th' : np.cos(np.deg2rad(0)),
      'fnum' : 56.6,
      'angle' : 0.23
      }

# IODINE CELL 

# PostFlight
I2_517_pos = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V_CELL/F517"
I2_502_pos = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V_CELL/F52502"
I2_506_pos = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V_CELL/F52506"

# PREFILTERS
PF_517_mps = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V/F517"
PF_502_mps = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V/F52502"
PF_506_mps = "/media/pablo/T7 Touch/PREFILTERS/MPS_2023/SCAN_V/F52506"

PF_517_kiru = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2024/F517"
PF_502_kiru = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2024/F52502"
PF_506_kiru = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2024/F52506"

PF_517_kiru_2022 = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2022/Scan1/517_27"
PF_502_kiru_2022 = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2022/Scan1/525_02"
PF_506_kiru_2022 = "/media/pablo/T7 Touch/PREFILTERS/KIRU_2022/Scan1/525_06"

mps = [PF_517_mps, PF_502_mps, PF_506_mps]
kiru_2024 = [PF_517_kiru, PF_502_kiru, PF_506_kiru]
kiru_2022 = [PF_517_kiru_2022, PF_502_kiru_2022, PF_506_kiru_2022]


int_mps = np.zeros((3, 143))
vol_mps = np.zeros((3, 143))

int_k24 = np.zeros((3, 143))
vol_k24 = np.zeros((3, 143))

int_k22 = np.zeros((3, 143))
vol_k22 = np.zeros((3, 143))

for i in range(3):
    all_mps = sorted(glob.glob(f"{mps[i]}/*_0_*.img"))
    all_k24 = sorted(glob.glob(f"{kiru_2024[i]}/*_0_*.img"))
    all_k22 = sorted(glob.glob(f"{kiru_2022[i]}/*_0_*.img"))


    for ind, img in enumerate(all_mps):
        print(f"mps:{ind}")
        I, H = read_Tumag(img)

        int_mps[i, ind] = np.mean(I[800:1200, 800:1200])
        vol_mps[i, ind] = (-1) ** (int(H['EtalonSign']) + 1) \
                            * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1))
        
    for ind, img in enumerate(all_k24):
        print(f"k24:{ind}")
        I, H = read_Tumag(img)

        int_k24[i, ind] = np.mean(I[800:1200, 800:1200])
        vol_k24[i, ind] = (-1) ** (int(H['EtalonSign']) + 1) \
                            * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1))
        
    for ind, img in enumerate(all_k22):
        print(f"k22:{ind}")
        I, H = read_Tumag(img)

        int_k22[i, ind] = np.mean(I[800:1200, 800:1200])
        vol_k22[i, ind] = (-1) ** (int(H['EtalonSign']) + 1) \
                            * (int(H['EtalonDN']) * 4999 / (2 ** 12 - 1))

rest = [-2241, 2070, -2673]


mosaic = """
            ADD
            BEE
            CFF
"""

indexes = ["A", "B", "C"]
fig, axs = plt.subplot_mosaic(mosaic=mosaic, figsize = (12.5, 11.5))
titles = ["517 nm pre-filter.", "525.02 nm pre-filter.", "525.06 nm pre-filter."]
for i in range(3):
    axs[indexes[i]].plot(vol_mps[i][:-2], int_mps[i][:-2] / np.max(int_mps[i][:-2]), c = 'indigo', lw = 3, label = "PFI AIV 2024")
    axs[indexes[i]].plot(vol_k24[i][:-2], int_k24[i][:-2] / np.max(int_k24[i][:-2]), c = 'orange', lw = 3, label = "System AIV 2024")
    axs[indexes[i]].plot(vol_k22[i][:-2], int_k22[i][:-2] / np.max(int_k22[i][:-2]), c = 'crimson', lw = 3, label = "System AIV 2022")

    axs[indexes[i]].set_title(titles[i])
    axs[indexes[i]].grid(True, c = 'k', alpha = 0.3)
    axs[indexes[i]].set_xlabel("Voltage [V]")
    axs[indexes[i]].set_ylabel("Intensity [Norm.]")

    

    axs[indexes[i]].plot([rest[i], rest[i]], [0, 1.5], c = 'k', ls = '--', lw = 3, label = 'Line core')

axs["A"].set_xlim(-3500, 0)
axs["B"].set_xlim(0, 3500)
axs["C"].set_xlim(-3500, 0)    

axs["A"].set_ylim(0.3, 1.05)
axs["B"].set_ylim(0.7, 1.05)
axs["C"].set_ylim(0.8, 1.05)

axs["A"].legend(loc = "lower left", edgecolor = 'k')
#plt.savefig("plots/prefilters_tumag.pdf", bbox_inches = "tight")




# INTA POSTFLIGHT
PF517_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V/INTA_SCAN_F517_CAM1.mat"
PF52502_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V/INTA_SCAN_F52502_CAM1.mat"
PF52506_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V/INTA_SCAN_F52506_CAM1.mat"

cell517_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V_CELL/INTA_SCAN_cell_F517_CAM1.mat"
cell52502_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V_CELL/INTA_SCAN_cell_F52502_CAM1.mat"
cell52506_file = "/home/pablo/Desktop/TuMag/Data/TuMag integration at MPS/2023-11-21 Alberto/SCAN_V_CELL/INTA_SCAN_cell_F52506_CAM1.mat"


AtlasFile = '/media/pablo/T7 Touch/SUNRISE_3/OldTests/TuMAG/INTA_2021/E2E_port/Spectral_Calibration/INPUTS/Large_I_Spectrum.txt' # File containing Spectrum reference
# Loading I2 Spectrum
Reference_data = np.loadtxt(AtlasFile)
Wavelength_vector =  1 / Reference_data[:, 0] * 1e8
# Displacement = 5251.407174 - 5249.953321
Displacement = 0
Wavelength_vector -= Displacement
# Wavelength range and central wavelengths for the Spectrum reference [A]
low_lim_25 = 5250.6 - 1.5
up_lim_25 = 5250.6 + 1.5

low_lim_17 = 5171.3 + 0.25
up_lim_17 = 5175 + 0.25

idx_25 = (np.where(up_lim_25 > Wavelength_vector[np.where(low_lim_25 < Wavelength_vector)]) \
       + np.where(low_lim_25 < Wavelength_vector)[0][0])[0]
idx_17 = (np.where(up_lim_17 > Wavelength_vector[np.where(low_lim_17 < Wavelength_vector)]) \
       + np.where(low_lim_17 < Wavelength_vector)[0][0])[0]
    
wls_25 = Wavelength_vector[idx_25]
Spectrum_25 = Reference_data[idx_25, 1]

Spectrum_25 = Spectrum_25[::-1]
Ref_wls_25 = wls_25[::-1]

wls_17 = Wavelength_vector[idx_17]
Spectrum_17 = Reference_data[idx_17, 1]

Spectrum_17 = Spectrum_17[::-1]
Ref_wls_17 = wls_17[::-1]

Et['R'] = 0.892 
central_wls = (Ref_wls_25[-1] - Ref_wls_25[0]) / 2 + Ref_wls_25[0] + 0.014
Etalon = ft.Psi_col(Ref_wls_25 * 1E-10, 1, central_wls * 1E-10, Et)

# f = FilterProperties['CWL2']['Interp'](wls_25)

Convolution_25 = convolve(Spectrum_25, (Etalon ** 2) ,  mode = 'same')
Convolution_25 /= np.max(Convolution_25)

Et['R'] = 0.892
central_wls = (Ref_wls_17[-1] - Ref_wls_17[0]) / 2 + Ref_wls_17[0] + 0.014
Etalon = ft.Psi_col(Ref_wls_17 * 1E-10, 1, central_wls * 1E-10, Et)

# f = FilterProperties['CWL3']['Interp'](wls_17)

Convolution_17 = convolve(Spectrum_17, (Etalon ** 2) ,  mode = 'same')
Convolution_17 /= np.max(Convolution_17)




pf517 = loadmat(PF517_file)
pf52502 = loadmat(PF52502_file)
pf52506 = loadmat(PF52506_file)

cell517 = loadmat(cell517_file)
cell52502 = loadmat(cell52502_file)
cell52506 = loadmat(cell52506_file)


Volts_vec = np.arange(-4000, 4000, 50)






ax2 = axs["D"].inset_axes([0, 0, 1, 1], transform=axs["D"].transAxes, facecolor='none') # overlay ax2 on top of ax1
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="crimson")
ax2.tick_params(axis='y', colors="crimson")


axs["D"].plot(cell517['HV'], cell517['Mwo'][:,3] / pf517['Mwo'][:,3], color = 'darkorange', lw = 3, label = "PFI AIV measure.")
ax2.plot(Ref_wls_17 - 1.4, Convolution_17, color = "indigo", label = 'Simulation', linewidth = 3)
axs["D"].set_xlim(-3500, 1000)
ax2.set_xlim(5172.38, 5173.73)
axs["D"].set_ylim(0.58, 0.87)
ax2.set_ylim(0.5, 1.1)


ax2.set_xlabel(r"Wavelengths [$\AA$]", color = "crimson")
axs["D"].set_xlabel("Volts [V]")
axs["D"].set_ylabel("Intensity [Norm.]")
axs["D"].plot([2000, 2000], [0, 1], c = 'indigo', lw = 3, label = "Simulation")
axs["D"].legend(edgecolor = 'k')
ax2.grid(True, c = 'k', alpha = 0.3)


ax2 = ax2 = axs["E"].inset_axes([0, 0, 1, 1], transform=axs["E"].transAxes, facecolor='none') # overlay ax2 on top of ax1
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="crimson")
ax2.tick_params(axis='y', colors="crimson")

ax2.plot(Ref_wls_25, Convolution_25, color = "indigo", label = 'Simulation', linewidth = 3)
axs["E"].plot(cell52502['HV'], cell52502['Mwo'][:,3] / pf52502['Mwo'][:,3], color = 'darkorange', lw = 3)

axs["E"].set_xlim(0, 3200)
ax2.set_xlim(5251, 5251.97)

axs["E"].set_ylim(0.58, 0.87)
ax2.set_ylim(0.5, 1.1)

ax2.set_xlabel(r"Wavelengths [$\AA$]", color = "crimson")
axs["E"].set_xlabel("Volts [V]")
axs["E"].set_ylabel("Intensity [Norm.]")
ax2.grid(True, c = 'k', alpha = 0.3)



ax2 = ax2 = axs["F"].inset_axes([0, 0, 1, 1], transform=axs["F"].transAxes, facecolor='none') # overlay ax2 on top of ax1
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="crimson")
ax2.tick_params(axis='y', colors="crimson")


ax2.plot(Ref_wls_25, Convolution_25, color = "indigo", label = 'Simulation', linewidth = 3)
axs["F"].plot(cell52506['HV'], cell52506['Mwo'][:,3] / pf52506['Mwo'][:,3], color = 'darkorange', lw = 3)


axs["F"].set_xlim(-2000, 1000)
ax2.set_xlim(5250.4, 5251.3)

axs["F"].set_ylim(0.6, 0.81)
ax2.set_ylim(0.6, 1.05)


ax2.set_xlabel(r"Wavelengths [$\AA$]", color = "crimson")
axs["F"].set_xlabel("Volts [V]")
axs["F"].set_ylabel("Intensity [Norm.]")
ax2.grid(True, c = 'k', alpha = 0.3)


plt.tight_layout()
plt.savefig("plots/Spectroscopic_calibration.pdf", bbox_inches = "tight")

