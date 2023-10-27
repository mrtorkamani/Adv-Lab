#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 21:32:24 2023

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

saved_address_B = []
saved_address_R = []
saved_address_V = []
FIT_files_SCIENCE_B = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_B/*.FIT')
FIT_files_SCIENCE_R = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_R/*.FIT')
FIT_files_SCIENCE_V = glob.glob('/home/mr/Desktop/Lab/S263/Photometry/M34/SCIENCE_V/*.FIT')
for i in FIT_files_SCIENCE_B:
    saved_address_B.append('/home/mr/Desktop/Lab/S263/codes/produced data/SCIENCE/SCIENCE_B/'+i[51:])    
for i in FIT_files_SCIENCE_R:
    saved_address_R.append('/home/mr/Desktop/Lab/S263/codes/produced data/SCIENCE/SCIENCE_R/'+i[51:])    
for i in FIT_files_SCIENCE_V:
    saved_address_V.append('/home/mr/Desktop/Lab/S263/codes/produced data/SCIENCE/SCIENCE_V/'+i[51:])    


storage_B = np.zeros((len(FIT_files_SCIENCE_B),2048,3072),dtype=np.float64)
storage_R = np.zeros((len(FIT_files_SCIENCE_R),2048,3072),dtype=np.float64)
storage_V = np.zeros((len(FIT_files_SCIENCE_V),2048,3072),dtype=np.float64)

DARK = np.zeros((2048,3072),dtype=np.float64)

FLAT_B = np.zeros((2048,3072),dtype=np.float64)
FLAT_R = np.zeros((2048,3072),dtype=np.float64)
FLAT_V = np.zeros((2048,3072),dtype=np.float64)

data_DARK, header_DARK = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Master_DARK.FIT', header=True, ext=0)
DARK = data_DARK

data_FLAT_B, header_FLAT_B = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Normal_Master_FLAT_B.FIT', header=True, ext=0)
FLAT_B = data_FLAT_B

data_FLAT_R, header_FLAT_R = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Normal_Master_FLAT_R.FIT', header=True, ext=0)
FLAT_R = data_FLAT_R

data_FLAT_V, header_FLAT_V = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Normal_Master_FLAT_V.FIT', header=True, ext=0)
FLAT_V = data_FLAT_V


i = 0
for DATA in FIT_files_SCIENCE_B:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_B[i] = data
    storage_B[i] = (storage_B[i]-DARK)/FLAT_B
    i += 1


i = 0
for DATA in FIT_files_SCIENCE_R:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_R[i] = data
    storage_R[i] = (storage_R[i]-DARK)/FLAT_R
    i += 1
    

i = 0
for DATA in FIT_files_SCIENCE_V:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_V[i] = data
    storage_V[i] = (storage_V[i]-DARK)/FLAT_V
    i += 1

for i in range(len(storage_B)):
    fits.writeto(saved_address_B[i], storage_B[i] , header = header, overwrite = True ,
                 output_verify = "ignore")

for i in range(len(storage_R)):
    fits.writeto(saved_address_R[i], storage_R[i] , header = header, overwrite = True ,
                 output_verify = "ignore")

for i in range(len(storage_V)):
    fits.writeto(saved_address_V[i], storage_V[i] , header = header, overwrite = True ,
                 output_verify = "ignore")



fig1, ax1 = plt.subplots(figsize=(16, 12))
ax1.set_title(r'\textbf{Reduced Science Frame B}')
# Add colorbar with custom ticks
cax1 =ax1.imshow(storage_B[3],cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Reduced_Science_Frame_B.png')
