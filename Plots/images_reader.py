# ============================ IMPORTS ====================================== #

# Built-in modules
import os
import sys
import numpy as np
from matplotlib import pyplot as plt

# Own Libs
sys.path.append('./functions')
sys.path.append('./imread')
from SPGCam_lib import *

# ============================= CONFIG ====================================== #

plt.style.use('dark_background')
# ImagePath = 'images/2021_10_07_07_28_36_230_0_0.img'

# all_images = os.listdir('images')
# PATHS = ['images/' + x for x in all_images]

# =========================================================================== #

def image_reader(ImagePath, PlotFlag):
    
    #Read image
    # image = getHeader(ImagePath)

    try:
        image = np.fromfile(ImagePath, dtype='<i4')
        image = image.reshape([2048,2048])
        
    except:
        image = np.fromfile(ImagePath, dtype='<i2')
        image = image.reshape([2048,2048])

    if PlotFlag:
        plt.imshow(image, cmap = 'magma')
        plt.colorbar()
        plt.show()
        
    return image

def multiple_images_reader(PATHS, Ndim):
    
    Nimages = len(PATHS)
    
    Images = np.zeros((Nimages, Ndim, Ndim))
    
    
    
    for index, path in enumerate(PATHS):
        im = image_reader(path, False)
        Images[index] = im
        
    return Images


"""
# Reading one image and plotting 
Im = image_reader(ImagePath, True)

# # Reading all images and storing them
All_im = multiple_images_reader(PATHS)
"""

