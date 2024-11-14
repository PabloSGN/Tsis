from astropy.io import fits
import matplotlib.pyplot as plt
import os
import glob
import matplotlib as mp

import sys
# Own Libs
sys.path.append(r"C:\Users\pablo\OneDrive\Escritorio\TuMag\TuMags_Reduction_Pipeline")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate, demodulate_v2
import alignment as al
from pd_functions_v22 import restore_ima

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "9"
import numpy as np

#def colorbar(ax):



all_files = glob.glob(os.path.join(r"C:\Users\pablo\OneDrive\Escritorio\Trabajo\Tumag\fits", "*"))


zkes = [ 0,0,0,0.47935868 , 0.08074037, -0.84290989, 0.0295268  , 0.15844359,  0.34845456,  0.45557947,  0.26175956,  0.50482008, 0.14049373 , 0.02425028 ,-0.0321224  -0.09221082 , 0.02892796, -0.11530846, -0.01152968, -0.03183491,  0.01621198]

timestamps = ["2024_07_12_12_38_33_090", "2024_07_15_22_07_58_037", "2024_07_13_17_34_46_418", "2024_07_10_19_14_12_463", "2024_07_13_11_44_01_385", "2024_07_13_23_12_18_331", "2024_07_11_06_52_13_264",  "2024_07_11_02_22_14_193",  "2024_07_15_08_43_41_810"]
Timelines = ["AR_1", "AR_5", "CLV1", "EMF_1", "FS_1", "PL_1", "QSDC_1", "SP_4", "SP_4"]

fig, axs = plt.subplots(3, 3, figsize = (11, 10))
# Loop through the subplots using a single index

for i in range(9):
    # Calculate row and column indices
    row = i // 3
    col = i % 3

    text = f"Timeline : {Timelines[i]}\nT.stamp: {timestamps[i]}"

    data = fits.getdata(all_files[i])[0]

    reconstructed, _ = restore_ima(data[200:-250, 200:-250], np.array(zkes))

    print(np.shape(reconstructed))

    
    axs[row,col].imshow(reconstructed, cmap = "gray")
    #plt.colorbar(im, fraction=0.046, pad=0.04)
    axs[row, col].set_xticks([])
    axs[row, col].set_yticks([])
    #axs[row, col].set_title(text)

    # axs[row, col].text(0.05, 0.95, text, transform=axs[row, col].transAxes,
    #                    verticalalignment='top', horizontalalignment='left', fontsize=10, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

# Adjust layout
plt.tight_layout()
plt.show()
#plt.savefig(r"C:\Users\pablo\OneDrive\Escritorio\TuMag\timelines_mosaic.pdf", bbox_inches = "tight")
