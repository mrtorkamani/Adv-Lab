#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 14:35:09 2023

@author: mr
"""

import numpy as np
import pandas as pd
import scipy.optimize as so
import matplotlib.pyplot as plt
import scipy.signal as sig
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit
import glob
from astropy.io import fits

plt.style.use('custom_style.mplstyle')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

colors = [(0, 0, 0), (1, 1, 1)] # first color is black, last is white
cm = LinearSegmentedColormap.from_list(
        'custom_colormap', colors, N=4000)


FIT_files_DARK = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/DARK/*.FIT')
FIT_files_FLAT_B = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/FLAT_B/*.FIT')
FIT_files_FLAT_R = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/FLAT_R/*.FIT')
FIT_files_FLAT_V = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/FLAT_V/*.FIT')
FIT_files_SCIENCE_B = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_B/*.FIT')
FIT_files_SCIENCE_R = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_R/*.FIT')
FIT_files_SCIENCE_V = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_V/*.FIT')

data_BIAS, header_BIAS = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_BIAS.FIT', header=True, ext=0)
storage = np.zeros((3,2048,3072),dtype=np.float64)

storage[1] = data_BIAS

for DATA in FIT_files_DARK:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")
    
for DATA in FIT_files_FLAT_B:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")
    
for DATA in FIT_files_FLAT_R:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")
    
for DATA in FIT_files_FLAT_V:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")

for DATA in FIT_files_SCIENCE_B:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")


for DATA in FIT_files_SCIENCE_R:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")


for DATA in FIT_files_SCIENCE_V:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[0] = data
    storage[2] = storage[0] - storage[1]
    fits.writeto(DATA, storage[2] , header = header, overwrite = True ,
                 output_verify = "ignore")

    