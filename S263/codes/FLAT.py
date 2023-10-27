#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 14:57:45 2023

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


FIT_files_FLAT = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/FLAT_V/*.FIT')

storage = np.zeros((len(FIT_files_FLAT),2048,3072),dtype=np.float64)
DARK = np.zeros((2048,3072),dtype=np.float64)

data, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_DARK.FIT', header=True, ext=0)
DARK = data

i = 0
for DATA in FIT_files_FLAT:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[i] = data
    storage[i] = storage[i]-DARK*.5/30
    i += 1

sumation = 0
for i in range(len(FIT_files_FLAT)):
    sumation = sumation + storage[i]
avarage = sumation/len(FIT_files_FLAT)

print(np.average(avarage))

fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Master_FLAT_V.FIT', avarage , header = header, overwrite = True ,
output_verify = "ignore")
fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Normal_Master_FLAT_V.FIT', avarage/np.average(avarage) , header = header, overwrite = True ,
output_verify = "ignore")


fig, ax = plt.subplots(figsize=(16, 12))
ax.set_title(r'\textbf{Master FLAT V}')
# Add colorbar with custom ticks
cax =ax.imshow(avarage,cmap=cm,origin='lower')
cbar = plt.colorbar(cax)
cbar.set_label(r'\textbf{Color Scale}')
fig.savefig(r'Master_FLAT_V.png')

fig1, ax1 = plt.subplots(figsize=(16, 12))
ax1.set_title(r'\textbf{Normalised Master FLAT V}')
# Add colorbar with custom ticks
cax1 =ax1.imshow(avarage/np.average(avarage),cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Normal_Master_FLAT_V.png')

