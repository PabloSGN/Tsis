import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os

import utils as ut

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "9"


folder = r"C:\Users\pablo\OneDrive\Escritorio\Tésis\TuMag data"



dark, _ = ut.read_Tumag(os.path.join(folder, "Dark.img"))
lpol, _ = ut.read_Tumag(os.path.join(folder, "Lpol.img"))
mpol, _ = ut.read_Tumag(os.path.join(folder, "Micro.img"))
phol, _ = ut.read_Tumag(os.path.join(folder, "Pinhole.img"))
stry , _= ut.read_Tumag(os.path.join(folder, "Straylight.img"))
PD, _ = ut.read_Tumag(os.path.join(folder, "PD", "2024_07_11_22_27_33_491_1_2501.img"))
F4, _ = ut.read_Tumag(os.path.join(folder, "PD", "2024_07_11_22_27_41_738_1_2541.img"))

ff = np.load(os.path.join(folder, "Flat_example.npy"))

np.shape(dark)


fig, axs = plt.subplots(2, 4, figsize = (10.5, 4.4))

im = axs[0, 0].imshow(dark, cmap = "gray", vmin = np.mean(dark) - 2 * np.std(dark),  vmax = np.mean(dark) + 2 * np.std(dark))
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[1, 0].imshow(ff, cmap = 'inferno')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[0, 1].imshow(phol, cmap = 'inferno', vmax = 3000)
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[1, 1].imshow(stry, cmap = 'inferno')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[0, 2].imshow(mpol, cmap = 'inferno')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[1, 2].imshow(lpol, cmap = 'inferno')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[0, 3].imshow(F4, cmap = 'gray')
plt.colorbar(im, fraction=0.046, pad=0.04)

im = axs[1, 3].imshow(PD, cmap = 'gray')
plt.colorbar(im, fraction=0.046, pad=0.04)

axs[0, 0].set_title("Dark current")
axs[1, 0].set_title("Flat-field")

axs[0, 1].set_title("Pinhole")
axs[1, 1].set_title("Straylight")

axs[0, 2].set_title("Micropolarizers")
axs[1, 2].set_title("Linear polarizer")

axs[0, 3].set_title("PD pair - focused")
axs[1, 3].set_title("PD pair - blurred")

for i in range(2):
    for j in range(4):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

plt.tight_layout()
plt.savefig(os.path.join("plots", "cal_modes_examples.pdf"), bbox_inches = "tight")




