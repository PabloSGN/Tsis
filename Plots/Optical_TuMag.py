import matplotlib.pyplot as plt
import glob
import matplotlib as mp
import numpy as np

from images_reader import image_reader
from utils import read_Tumag

import functions_old as ft

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"


folder = "/media/pablo/T7 Touch/SUNRISE_3/OldTests/TuMAG/INTA_2021/02122021/ImageE2E/USAF throughfocus_best position/Without PD/*img"

images = glob.glob(folder)



fig, axs = plt.subplots(2, 3, figsize = (10.5, 5.5))


I = image_reader(images[0], PlotFlag=False)
im = axs[0, 0].imshow(I[850:1025, 800:975], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[0, 0],  fraction=0.046, pad=0.04)


I = image_reader(images[1], PlotFlag=False)
im = axs[1, 0].imshow(I[850:1025, 1070:1245], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[1, 0],  fraction=0.046, pad=0.04)


I = image_reader(images[4], PlotFlag=False)
im = axs[0, 1].imshow(I[850:1025, 800:975], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[0, 1],  fraction=0.046, pad=0.04)


I = image_reader(images[5], PlotFlag=False)
im = axs[1, 1].imshow(I[850:1025, 1070:1245], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[1, 1],  fraction=0.046, pad=0.04)

I = image_reader(images[8], PlotFlag=False)
im = axs[0, 2].imshow(I[850:1025, 800:975], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[0, 2],  fraction=0.046, pad=0.04)


I = image_reader(images[9], PlotFlag=False)
im = axs[1, 2].imshow(I[850:1025, 1070:1245], cmap = 'inferno', vmax =1200)
fig.colorbar(im, ax=axs[1, 2],  fraction=0.046, pad=0.04)

for i in range(2):
    for j in range(3):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])


axs[0, 0].set_title("517 nm Pre-filter")
axs[0, 1].set_title("525.06 nm Pre-filter")
axs[0, 2].set_title("525.02 nm Pre-filter")

axs[0, 0].set_ylabel("Camera 1")
axs[1, 0].set_ylabel("Camera 2")

plt.tight_layout()
plt.savefig("plots/USAF_E2E.pdf", bbox_inches = 'tight')