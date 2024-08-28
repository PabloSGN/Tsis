import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from utils import read_Tumag

# CONFIG

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
fig, axs = plt.subplots(1, 3, figsize = (15.5, 5))
titles = ["517 nm pre-filter.", "525.02 nm pre-filter.", "525.06 nm pre-filter."]
for i in range(3):
    axs[i].plot(vol_mps[i][:-2], int_mps[i][:-2] / np.max(int_mps[i][:-2]), c = 'indigo', lw = 3, label = "PFI AIV 2024")
    axs[i].plot(vol_k24[i][:-2], int_k24[i][:-2] / np.max(int_k24[i][:-2]), c = 'orange', lw = 3, label = "System AIV 2024")
    axs[i].plot(vol_k22[i][:-2], int_k22[i][:-2] / np.max(int_k22[i][:-2]), c = 'crimson', lw = 3, label = "System AIV 2022")

    axs[i].set_title(titles[i])
    axs[i].grid(True, c = 'k', alpha = 0.3)
    axs[i].set_xlabel("Voltage [V]")
    axs[i].set_ylabel("Intensity [Norm.]")

    

    axs[i].plot([rest[i], rest[i]], [0, 1.5], c = 'k', ls = '--', lw = 2, label = 'Line core')



axs[0].set_xlim(-3500, 0)
axs[1].set_xlim(0, 3500)
axs[2].set_xlim(-3500, 0)    



axs[0].set_ylim(0.3, 1.05)
axs[1].set_ylim(0.7, 1.05)
axs[2].set_ylim(0.7, 1.05)

plt.tight_layout()
axs[0].legend(loc = "lower left", edgecolor = 'k')
plt.show()






