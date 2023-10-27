#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 22:35:13 2023

@author: mr
"""

import scipy.ndimage as sciim
from image_registration import chi2_shift
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


directory = '/home/mr/Desktop/Lab/S263/codes/produced data/SCIENCE'

FIT_files_SCIENCE_B = glob.glob(directory+'/SCIENCE_B/*.FIT')
FIT_files_SCIENCE_R = glob.glob(directory+'/SCIENCE_R/*.FIT')
FIT_files_SCIENCE_V = glob.glob(directory+'/SCIENCE_V/*.FIT')





storage_B = np.zeros((len(FIT_files_SCIENCE_B),2048,3072),dtype=np.float64)
storage_R = np.zeros((len(FIT_files_SCIENCE_R),2048,3072),dtype=np.float64)
storage_V = np.zeros((len(FIT_files_SCIENCE_V),2048,3072),dtype=np.float64)



i = 0
for DATA in FIT_files_SCIENCE_B:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_B[i] = data
    i += 1

i = 0
for DATA in FIT_files_SCIENCE_R:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_R[i] = data
    i += 1

i = 0
for DATA in FIT_files_SCIENCE_V:
    data, header = fits.getdata(DATA, header=True, ext=0)
    storage_V[i] = data
    i += 1


refrence_data = storage_B[0]

i = 0
for data in storage_B:  
    xoff ,yoff ,exoff ,eyoff = chi2_shift(data , refrence_data , return_error = True)
    storage_B[i] = sciim.interpolation.shift(data ,(yoff , xoff))
    i += 1 

sumation = 0
for i in range(len(FIT_files_SCIENCE_B)):
    sumation = sumation + storage_B[i]
avarage_B = sumation/len(FIT_files_SCIENCE_B)
################
i = 0
for data in storage_R:  
    xoff ,yoff ,exoff ,eyoff = chi2_shift(data , refrence_data , return_error = True)
    storage_R[i] = sciim.interpolation.shift(data ,(yoff , xoff))
    i += 1 

sumation = 0
for i in range(len(FIT_files_SCIENCE_R)):
    sumation = sumation + storage_R[i]
avarage_R = sumation/len(FIT_files_SCIENCE_R)
##################
i = 0
for data in storage_V:  
    xoff ,yoff ,exoff ,eyoff = chi2_shift(data , refrence_data , return_error = True)
    storage_V[i] = sciim.interpolation.shift(data ,(yoff , xoff))
    i += 1 

sumation = 0
for i in range(len(FIT_files_SCIENCE_V)):
    sumation = sumation + storage_V[i]
avarage_V = sumation/len(FIT_files_SCIENCE_V)

fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Co_added_B', avarage_B , header = header, overwrite = True ,
                 output_verify = "ignore")
fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Co_added_R', avarage_R , header = header, overwrite = True ,
                 output_verify = "ignore")
fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/Co_added_V', avarage_V , header = header, overwrite = True ,
                 output_verify = "ignore")

fig1, ax1 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax1 =ax1.imshow(avarage_B,cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Co-added_reduced_master_B.png')

fig2, ax2 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax2 =ax2.imshow(avarage_R,cmap=cm,origin='lower')
cbar2 = plt.colorbar(cax1)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'Co-added_reduced_master_R.png')

fig3, ax3 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax3 =ax3.imshow(avarage_V,cmap=cm,origin='lower')
cbar3 = plt.colorbar(cax1)
cbar3.set_label(r'\textbf{Color Scale}')
fig3.savefig(r'Co-added_reduced_master_V.png')