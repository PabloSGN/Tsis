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


folder = "/media/pablo/T7 Touch/SUNRISE_3/OldTests/PostFlight/Opticos/PD/Random Dot/*"

images = glob.glob(folder)

fig, axs = plt.subplots(2, 2, figsize = (9, 8))

I, H = read_Tumag(images[1])
im = axs[0, 0].imshow(I, cmap = "inferno", vmax = 1000)
fig.colorbar(im, ax=axs[0, 0],  fraction=0.046, pad=0.04)

im = axs[0, 1].imshow(I[700:700+640, 700:700+640], cmap = "inferno", vmax = 1000)
fig.colorbar(im, ax=axs[0, 1],  fraction=0.046, pad=0.04)

I, H = read_Tumag(images[5])
im = axs[1, 0].imshow(I, cmap = "inferno", vmax = 500)
fig.colorbar(im, ax=axs[1, 0],  fraction=0.046, pad=0.04)

im = axs[1, 1].imshow(I[700:700+640, 700:700+640], cmap = "inferno", vmax = 500)
fig.colorbar(im, ax=axs[1, 1],  fraction=0.046, pad=0.04)

for i in range(2):
    for j in range(2):
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])


def plot_square(ax, low, up):
    ax.plot([low, low], [low, up], c = 'w', ls = "--", lw = 1, label = 'Zoomed region')
    ax.plot([up, up], [low, up], c = 'w', ls = "--", lw = 1)
    ax.plot([low, up], [low, low], c = 'w', ls = "--", lw = 1)
    ax.plot([low, up], [up, up], c = 'w', ls = "--", lw = 1)

plot_square(axs[0, 0], 700, 700 + 640)
plot_square(axs[1, 0], 700, 700 + 640)

axs[0,0].legend(edgecolor = 'w')

axs[0, 0].set_title("Full-sized image")
axs[0, 1].set_title("Zoomed-in region")

axs[0, 0].set_ylabel("Without PD")
axs[1, 0].set_ylabel("With PD")

plt.tight_layout()
plt.savefig("plots/random_dot.pdf", bbox_inches = "tight")