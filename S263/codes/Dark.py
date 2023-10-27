#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 21:37:20 2023

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


FIT_files = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/DARK/*.FIT')

storage = np.zeros((len(FIT_files),2048,3072),dtype=np.float64)

i = 0
for DATA in FIT_files:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage[i] = data
    i += 1

sumation = 0
for i in range(len(FIT_files)):
    sumation = sumation + storage[i]
avarage = sumation/len(FIT_files)

print(np.min(avarage))

fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Master_DARK.FIT', avarage , header = header, overwrite = True ,
output_verify = "ignore")


fig, ax = plt.subplots(figsize=(16, 12))
ax.set_title(r'\textbf{Master DARK}')
# Add colorbar with custom ticks
cax =ax.imshow(avarage,cmap=cm,origin='lower')
cbar = plt.colorbar(cax)
cbar.set_label(r'\textbf{Color Scale}')

fig.savefig(r'Master_DARK.png')
