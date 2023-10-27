#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 20:25:28 2023

@author: mr
"""

from datetime import datetime
from time import time
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pickle
import sys
datain = []


try:
	datatab = ascii.read("/home/mr/Desktop/Lab/S262/Data/A7/A7/stilldata_out.txt") #Put your Still Scan data file here
except Exception as e:
	print(e)
    
    
#########Getting the Sky Coordinates and inputting the Fluxes from the two telescopes###################

coords = [SkyCoord(d[0]+" "+d[1],unit=(u.hourangle,u.deg))  for d in zip(datatab['col3'],datatab['col4'])]

fluxa = list(datatab['col5'])
fluxb = list(datatab['col8'])
I = list(datatab['col6'])
Q = list(datatab['col7'])
nr = list(datatab['col8'])
coordpair = zip(coords,nr)
loadpickle = 0

pfluxa =[list(f) for f in zip(fluxa,range(len(fluxa)))]
pfluxb =[ list(f) for f in zip(fluxb,range(len(fluxb)))]
pQ =[list(f) for f in  zip(Q,range(len(Q)))]
pI =[list(f) for f in  zip(I,range(len(I)))]

##########################Loop to remove noise from the data###########################
i=8
while i:
	ind = np.where(abs(fluxa - np.roll(fluxa,1**i)) > 60)[0]
	pfluxa = np.delete(pfluxa,ind,axis=0)
	fluxa = np.delete(fluxa,ind)
	print(len(ind))
	ind = np.where(abs(fluxb - np.roll(fluxb,1**i)) > 60)[0]
	pfluxb = np.delete(pfluxb,ind,axis=0)
	fluxb = np.delete(fluxb,ind)

	ind = np.where(abs(Q - np.roll(Q,1**i)) > 30)[0]
	pQ = np.delete(pQ,ind,axis=0)
	Q = np.delete(Q,ind)
	
	ind = np.where(abs(I - np.roll(I,1**i)) > 30)[0]
	pI = np.delete(pI,ind,axis=0)
	I = np.delete(I,ind)


	i-=1

plt.style.use('custom_style.mplstyle')


plt.plot(fluxa)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Amplitude}')
plt.title(r'\textbf{Denoised Data antena A}')
plt.savefig('still_Denoised_Data_antenna_A.png')
plt.show()

plt.plot(fluxb)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Amplitude}')
plt.title(r'\textbf{Denoised Data antenna B}')
plt.savefig('still_Denoised_Data_antena_B.png')
plt.show()