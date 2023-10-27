#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 13:40:40 2023

@author: mr
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import astropy
import csv
from scipy.optimize import curve_fit
from astropy.io import fits
from pybaselines.utils import gaussian
from pybaselines.polynomial import poly

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


plt.style.use('custom_style.mplstyle')

galaxy_file = '/home/mr/Desktop/Lab/S264/S624_Data_A7/sto25_S7_spec_685.csv'

data = pd.read_csv(galaxy_file,sep=' ') 
channel = pd.to_numeric(data[data.columns[0]])
T = pd.to_numeric(data[data.columns[6]])


non_peaks = (
    (channel < 1450) | (channel > 1600)
)
    
x_masked = channel[non_peaks]
y_masked = T[non_peaks]

# fit only the masked x and y
_, params = poly(y_masked, x_masked, poly_order=3, return_coef=True)
# recreate the polynomial using numpy and the full x-data
baseline = np.polynomial.Polynomial(params['coef'])(channel)
#print(params)
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(channel, T, color = 'black')
ax.plot(channel, baseline,'--', color= 'red')
#plt.xlim(1100,1900)
ax.set_xlabel('Channel')
ax.set_ylabel(r'$T_{A}$ [K]')
ax.grid()
fig.savefig('calibration.png')

T_sub = T-baseline
fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(channel, T_sub , color ='black')
ax1.set_xlabel('Channel')
ax1.set_ylabel(r'$T_{A}$ [K]')
ax1.grid()
fig1.savefig('calibration_baseline_removal.png')

errbase= np.polynomial.Polynomial(params['coef'])(x_masked)


sam = (y_masked - errbase)**2
S=sum(sam)
dT = np.sqrt(S/len(errbase))
print('Delta T_A =',dT)
print('T_A =',max(T_sub))

alpha = 96/max(T_sub)
dalpha = 96/max(T_sub)**2 * dT

print(fr'alpha = {alpha:0.2f} +- {dalpha:0.2f}')
