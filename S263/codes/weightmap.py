#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 15:27:45 2023

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

data = np.zeros((3,2048,3072),dtype=np.float64)

data[0], header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_FLAT_B.FIT', header=True, ext=0)
data[1], header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_FLAT_R.FIT', header=True, ext=0)
data[2], header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_FLAT_V.FIT', header=True, ext=0)

fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/weightmap_FLAT_B.FIT', data[0]/np.max(data[0]) , header = header, overwrite = True ,
output_verify = "ignore")
fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/weightmap_FLAT_R.FIT', data[1]/np.max(data[1]) , header = header, overwrite = True ,
output_verify = "ignore")
fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/weightmap_FLAT_V.FIT', data[2]/np.max(data[2]) , header = header, overwrite = True ,
output_verify = "ignore")


fig, ax = plt.subplots(figsize=(16, 12))
ax.set_title(r'\textbf{Weight map FLAT B}')
# Add colorbar with custom ticks
cax =ax.imshow(data[0]/np.max(data[0]),cmap=cm,origin='lower')
cbar = plt.colorbar(cax)
cbar.set_label(r'\textbf{Color Scale}')
fig.savefig(r'weightmap_FLAT_B.png')

fig1, ax1 = plt.subplots(figsize=(16, 12))
ax1.set_title(r'\textbf{Weight map FLAT R}')
# Add colorbar with custom ticks
cax1 =ax1.imshow(data[1]/np.max(data[1]),cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'weightmap_FLAT_R.png')

fig2, ax2 = plt.subplots(figsize=(16, 12))
ax2.set_title(r'\textbf{Weight map FLAT V}')
# Add colorbar with custom ticks
cax2 =ax2.imshow(data[2]/np.max(data[2]),cmap=cm,origin='lower')
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'weightmap_FLAT_V.png')

