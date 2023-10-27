#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:01:17 2023

@author: mr
"""

from time import time
import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.interpolate import griddata
import pickle

import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

datain = []

try:
	datatab = ascii.read('/home/mr/Desktop/Lab/S262/Data/A7/A7/movedata_out.txt') #Put your Declination Scan datafile in here
except Exception as e:
	print(e)
    
##################################Define a Coordinate System###############################################

coords = [SkyCoord(d[0]+" "+d[1],unit=(u.hourangle,u.deg))  for d in zip(datatab['col3'],datatab['col4'])]

###########################################################################################################

######Input the Flux Values##########
fluxa = list(datatab['col5'])
fluxb = list(datatab['col8'])
I = list(datatab['col6'])
Q = list(datatab['col7'])
nr = list(datatab['col8'])

coordpair = zip(coords,nr)
loadpickle = 0
######################################

#create lists which hold the recorded values associated with an
#index number
pfluxa =[list(f) for f in zip(fluxa,range(len(fluxa)))]
pfluxb =[ list(f) for f in zip(fluxb,range(len(fluxb)))]
pQ =[list(f) for f in  zip(Q,range(len(Q)))]
pI =[list(f) for f in  zip(I,range(len(I)))]
#######################################################################

########extracting data from coords in deg############
ra = [coords[i].ra.deg for i in range(len(coords))]
dec = [coords[i].dec.deg for i in range(len(coords))]

#ploting#
plt.style.use('custom_style.mplstyle')

plt.plot(ra,dec)
plt.xlabel(r'\textbf{RA ($\mathrm{deg}$)}')
plt.ylabel(r'\textbf{Dec ($\mathrm{deg}$)}')
#plt.title(r'\textbf{}')
plt.savefig('decscan.png')
plt.show()

plt.plot(fluxa)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Intensity}')
plt.title(r'\textbf{Raw Data antenna A}')
plt.savefig('Raw_Data_antena_A.png')
plt.show()

plt.plot(fluxb)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Intensity}')
plt.title(r'\textbf{Raw Data antenna B}')
plt.savefig('Raw_Data_antena_B.png')
plt.show()


###############Loop to remove noise from the data#########################
i=20
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
############################################################################
    
#plotting#


plt.plot(fluxa)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Amplitude}')
plt.title(r'\textbf{Denoised Data antenna A}')
plt.savefig('Denoised_Data_antena_A.png')
plt.show()

plt.plot(fluxb)
plt.xlabel(r'\textbf{Time}')
plt.ylabel(r'\textbf{Amplitude}')
plt.title(r'\textbf{Denoised Data antenna B}')
plt.savefig('Denoised_Data_antena_B.png')
plt.show()

#Selecting the data over the range of observation with the Sun being in the field of view
try:
	print("=== Choose data range after looking at the plot ===")
	data_start = int(datatab['col9'][5000]) #Put in the datarange here
	data_end =  int(datatab['col9'][50000]) #Put in the datarange here
	fluxa = fluxa[data_start:data_end] 
	fluxb = fluxb[data_start:data_end]
	Q     = Q[data_start:data_end]
	I     = I[data_start:data_end]
	pfluxa = pfluxa[data_start:data_end]
	pfluxb = pfluxb[data_start:data_end]
	pQ     = pQ[data_start:data_end]
	pI     = pI[data_start:data_end]
	coords = coords[data_start:data_end]
except Exception as e:
	print(e)
    
meanfluxa = np.mean(fluxa)
meanfluxb = np.mean(fluxb)
meanQ = np.mean(Q)
meanI = np.mean(I)
#using the value index pairs interpolate over the missing data 

ifluxa = np.interp(range(data_start,data_end),[p[1] for p in pfluxa],[p[0]-meanfluxa for p in pfluxa])
ifluxb = np.interp(range(data_start,data_end),[p[1] for p in pfluxb],[p[0]-meanfluxb for p in pfluxb])
iQ = np.interp(range(data_start,data_end),[p[1] for p in pQ],[p[0]-meanQ for p in pQ])
iI = np.interp(range(data_start,data_end),[p[1] for p in pI],[p[0]-meanI for p in pI])

datdic = {}
for c,fa,fb,ii,qq in zip(coords,ifluxa,ifluxb,iI,iQ):
	try:
		datdic[(c.ra.value,c.dec.value)].append([fa,fb,ii,qq])
		#print "a"
	except:
		datdic[(c.ra.value,c.dec.value)]=[[fa,fb,ii,qq]]

mediandat = {}
dcoords = []
dvalues = []
for d in datdic:
	medi = np.median(datdic[d],axis=0)
	mediandat[d] = medi
	dcoords.append(d)
	dvalues.append(medi)

#load file
if loadpickle:
	with open("cleansig.txt",'r') as fin:
		print("Loading pickled file")
		dvalues = pickle.load(fin)
		dcoords = [coords[c[1]] for c in dvalues[0]]
        
        
#####################Create grid to plot the image###########################
ramax = max([c.ra.value for c  in coords])
ramin = min([c.ra.value for c  in coords])
decmax = max([c.dec.value for c  in coords])
decmin = min([c.dec.value for c  in coords])
#create the grid
grid_ra,grid_dec = np.mgrid[ramin:ramax:50j,decmin:decmax:50j]
#use a time stamp for newly created files so we don't overwrite:
##############################Produce Image from Antenna A#############################################
plt.figure(figsize=(4,6))
grid1 = griddata(dcoords,[d[0] for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-2,2:-2],extent = [ramin,ramax,decmin,decmax],cmap= 'inferno')
plt.title("Antenna A Sun Scan")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.grid(True)
plt.savefig("7.AntennaAscan.png", dpi=300, bbox_inches='tight')
plt.show()
###############################Produce Image from Antenna B#############################################
plt.figure(figsize=(4,6))
grid1 = griddata(dcoords,[d[1] for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-2,2:-2],extent = [ramin,ramax,decmin,decmax],cmap= 'inferno')
plt.title("Antenna B Sun Scan")

plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.grid(True)
plt.savefig("8.AntennaBscan.png", dpi=300, bbox_inches='tight')
plt.show()

#####################################Interferometric Image###################################################
plt.figure(figsize=(4,6))
grid_ra,grid_dec = np.mgrid[ramin:ramax:150j,decmin:decmax:100j]
grid1 = griddata(dcoords,[abs(d[2]+d[3]*1j)**2 for d in dvalues],(grid_ra,grid_dec),method = 'cubic')
plt.imshow(grid1.T[2:-1,2:-1],extent = [ramin,ramax,decmin,decmax],cmap= 'inferno')
plt.tricontour([d[0] for d in dcoords],[d[1] for d in dcoords],[abs(d[2]+d[3]*1j)**2 for d in dvalues])
plt.title("Interferometer Scan")
plt.xlabel("RA [deg]")
plt.ylabel("Dec [deg]")
plt.axvline(x=147,color='black')
plt.axvline(x=146.35 ,color='black')
plt.grid(True)
plt.savefig("9.Interfscan.png", dpi=300, bbox_inches='tight')
plt.show()
