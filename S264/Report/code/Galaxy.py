#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 09:26:16 2023

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
from scipy import integrate

plt.style.use('custom_style.mplstyle')


F0 = 1420.41e6
c = 299792.458
def FtoV(x):
    return (1-x/F0)*c

def VtoF(x):
    return (1 - x/c)*F0

galaxy_file = '/home/mr/Desktop/Lab/S264/S624_Data_A7/sto25_NGC2403_spec_686.csv'
offset_file = '/home/mr/Desktop/Lab/S264/S624_Data_A7/sto25_NGC2403_spec_688.csv'
data = pd.read_csv(galaxy_file,sep=' ') 
channel = pd.to_numeric(data[data.columns[0]])
T = pd.to_numeric(data[data.columns[6]])
F = pd.to_numeric(data[data.columns[1]])
channel_galaxy = channel[:700]
T_galaxy = T[:700]
F_galaxy = F[:700]

dataoff = pd.read_csv(offset_file,sep=' ') 
channeloff = pd.to_numeric(dataoff[dataoff.columns[0]])
Toff = pd.to_numeric(dataoff[dataoff.columns[6]])
Foff = pd.to_numeric(dataoff[dataoff.columns[1]])
channel_off = channeloff[:700]
T_off = Toff[:700]
F_off = Foff[:700]

non_peaks_off = (
    (channel_off < 1500) | (channel_off > 1620)
)
    
x_masked_off = channel_off[non_peaks_off]
y_masked_off = T_off[non_peaks_off]

# fit only the masked x and y
_, params = poly(y_masked_off, x_masked_off, poly_order=5, return_coef=True)
# recreate the polynomial using numpy and the full x-data
baseline_off = np.polynomial.Polynomial(params['coef'])(channel_galaxy)

figoff, axoff = plt.subplots(figsize=(10, 6))

axoff.plot(F_off, T_off, color = 'black')
axoff.plot(F_off, baseline_off,'--', color= 'red')
#plt.xlim(1100,1900)
axoff.set_xlabel('Frequency')
axoff.set_ylabel(r'$T_{A}$ [K]')
axoff.grid()
figoff.savefig('offset.png')

T_sub_offset = T_off-baseline_off
figoff1, axoff1 = plt.subplots(figsize=(10, 6))
axoff1.plot(F_off, T_sub_offset , color ='black')
secax = axoff1.secondary_xaxis('top', functions=(FtoV, VtoF))
secax.set_xlabel(r'Velocity $\mathrm{[kms^{-1}]}$')
axoff1.set_xlabel('Frequency')
axoff1.set_ylabel(r'$T_{A}$ [K]')
axoff1.grid()
figoff1.savefig('offset_baseline_removal.png')



#x = np.linspace(1200, 1900, 1400)

#real_baseline = 5 + 15 * np.exp(-x / 400)


# bitwise "or" (|) and "and" (&) operators for indexing numpy array
non_peaks_galaxy = (
    (channel_galaxy < 1500) | (channel_galaxy > 1790)
)
    
x_masked_galaxy = channel_galaxy[non_peaks_galaxy]
y_masked_galaxy = T_galaxy[non_peaks_galaxy]

# fit only the masked x and y
_, params = poly(y_masked_galaxy, x_masked_galaxy, poly_order=5, return_coef=True)
# recreate the polynomial using numpy and the full x-data
baseline_galaxy = np.polynomial.Polynomial(params['coef'])(channel_galaxy)


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(F_galaxy, T_galaxy, linewidth = 1, color='black')
ax.plot(F_galaxy, baseline_galaxy ,'--', color='red')
ax.set_xlabel('Frequency')
ax.set_ylabel(r'$T_{A}$ [K]')
ax.grid()
fig.savefig('galaxy.png')


T_sub_galaxy = T_galaxy-baseline_galaxy
fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(F_galaxy, T_sub_galaxy , color ='black')
secax = ax1.secondary_xaxis('top', functions=(FtoV, VtoF))
secax.set_xlabel(r'Velocity $\mathrm{[kms^{-1}]}$')
ax1.set_xlabel('Frequency')
ax1.set_ylabel(r'$T_{A}$ [K]')
ax1.grid()
fig1.savefig('galaxy_baseline_removal.png')


sub = T_sub_galaxy-T_sub_offset

fig2, ax2 = plt.subplots(figsize=(10, 6))
ax2.plot(F_off, sub, color='black')

ax2.axvline(x=1419.17e6,ymin=-1.4,ls='--', color='blue', label="Galaxy's border")
ax2.axvline(x=VtoF(-15),ymin=-1.4,ls='--', color='blue', label="Galaxy's border")
ax2.axvline(x=VtoF(125),ymin=-1.4,ls='--', color='red', label="Galaxy's radial velocity")
secax = ax2.secondary_xaxis('top', functions=(FtoV, VtoF))
secax.set_xlabel(r'Velocity $\mathrm{[kms^{-1}]}$')
ax2.set_xlabel('Frequency')
ax2.set_ylabel(r'$T_{A}$ [K]')
ax2.legend()
ax2.grid()
#ax2.set_ylim(0,.1)
fig2.savefig('sub.png')


mask = ((F_off >  VtoF(27)) & (F_off < VtoF(-15))) 
sub[mask] = 0.048


fig3, ax3 = plt.subplots(figsize=(10, 6))
ax3.plot(F_off, sub*7.71/0.09, color='black')

ax3.axvline(x=1419.15e6,ymin=-1.4,ls='--', color='blue', label="Galaxy's border")
ax3.axvline(x=VtoF(-15),ymin=-1.4,ls='--', color='blue', label="Galaxy's border")
ax3.axvline(x=VtoF(125),ymin=-1.4,ls='--', color='red', label="Galaxy's radial velocity")
secax = ax3.secondary_xaxis('top', functions=(FtoV, VtoF))
secax.set_xlabel(r'Velocity $\mathrm{[kms^{-1}]}$')
ax3.set_xlabel('Frequency')
ax3.set_ylabel(r'$S$ [Jy]')
ax3.annotate('Interpolated data',
             xy = (VtoF(10), 4.25),
             xytext = (VtoF(105), 7),
             arrowprops = dict(facecolor = 'black', width = 0.2, headwidth = 8),
             horizontalalignment = 'center')
ax3.legend()
ax3.grid()
ax3.set_ylim(0,8)
fig3.savefig('subsub.png')

mask1 = ((F_off > 1419.15e6) & (F_off < VtoF(-15)))
area = sub[mask1]*7.71/0.09
Area = sum(area)
print(Area)
print(len(area))
A = integrate.trapezoid(area)
print(A)
#plt.savefig('galaxy.png')