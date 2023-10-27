#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 20:04:12 2023

@author: mr
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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

'''
data1, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Master_BIAS.FIT', header=True, ext=0)

fig1, ax1 = plt.subplots(figsize=(16, 12))
ax1.set_title(r'\textbf{Master BIAS}')
cax1 =ax1.imshow(data1,cmap=cm,origin='lower',vmin= 1931 ,vmax=1974)
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig('Master_BIAS.png',dpi = 100)

data2, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Master_DARK.FIT', header=True, ext=0)

fig2, ax2 = plt.subplots(figsize=(16, 12))
ax2.set_title(r'\textbf{Master DARK}')
cax2 =ax2.imshow(data2,cmap=cm,origin='lower',vmin= -16 ,vmax=16)
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig('Master_DARK.png',dpi = 100)

data1, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Co_added_B', header=True, ext=0)

fig1, ax1 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax1 =ax1.imshow(data1,cmap=cm,origin='lower',vmin= 53 ,vmax=102)
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Co-added_reduced_master_B.png',dpi = 100)

data2, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Co_added_R', header=True, ext=0)

fig2, ax2 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax2 =ax2.imshow(data2,cmap=cm,origin='lower',vmin= 133 ,vmax=203)
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'Co-added_reduced_master_R.png',dpi = 100)

data3, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Co_added_V', header=True, ext=0)

fig3, ax3 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
cax3 =ax3.imshow(data3,cmap=cm,origin='lower',vmin= 103 ,vmax=203)
cbar3 = plt.colorbar(cax3)
cbar3.set_label(r'\textbf{Color Scale}')
fig3.savefig(r'Co-added_reduced_master_V.png',dpi = 100)
'''
'''
data1, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Master_FLAT_B.FIT', header=True, ext=0)

fig1, ax1 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax1.set_title(r'\textbf{Master FLAT B')
cax1 =ax1.imshow(data1,cmap=cm,origin='lower',vmin= 27772 ,vmax=35321)
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Master_FLAT_B.png',dpi = 100)

data2, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Master_FLAT_R.FIT', header=True, ext=0)

fig2, ax2 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax2.set_title(r'\textbf{Master FLAT R}')
cax2 =ax2.imshow(data2,cmap=cm,origin='lower',vmin= 32772 ,vmax=40321)
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'Master_FLAT_R.png',dpi = 100)

data3, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Master_FLAT_V.FIT', header=True, ext=0)

fig3, ax3 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax3.set_title(r'\textbf{Master FLAT V}')
cax3 =ax3.imshow(data3,cmap=cm,origin='lower',vmin= 36772 ,vmax=44321)
cbar3 = plt.colorbar(cax3)
cbar3.set_label(r'\textbf{Color Scale}')
fig3.savefig(r'Master_FLAT_V.png',dpi = 100)
'''

'''
data1, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Normal_Master_FLAT_B.FIT', header=True, ext=0)

fig1, ax1 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax1.set_title(r'\textbf{Normalised Master FLAT B')
cax1 =ax1.imshow(data1,cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'Normal_Master_FLAT_B.png',dpi = 100)

data2, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Normal_Master_FLAT_R.FIT', header=True, ext=0)

fig2, ax2 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax2.set_title(r'\textbf{Normalised Master FLAT R}')
cax2 =ax2.imshow(data2,cmap=cm,origin='lower')
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'Normal_Master_FLAT_R.png',dpi = 100)

data3, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Normal_Master_FLAT_V.FIT', header=True, ext=0)

fig3, ax3 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax3.set_title(r'\textbf{Normalised Master FLAT V}')
cax3 =ax3.imshow(data3,cmap=cm,origin='lower')
cbar3 = plt.colorbar(cax3)
cbar3.set_label(r'\textbf{Color Scale}')
fig3.savefig(r'Normal_Master_FLAT_V.png',dpi = 100)
'''

data1, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/weightmap_FLAT_B.FIT', header=True, ext=0)

fig1, ax1 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax1.set_title(r'\textbf{Weight map FLAT B')
cax1 =ax1.imshow(data1,cmap=cm,origin='lower')
cbar1 = plt.colorbar(cax1)
cbar1.set_label(r'\textbf{Color Scale}')
fig1.savefig(r'weightmap_FLAT_B.png',dpi = 100)

data2, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/weightmap_FLAT_R.FIT', header=True, ext=0)

fig2, ax2 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax2.set_title(r'\textbf{Weight map FLAT R}')
cax2 =ax2.imshow(data2,cmap=cm,origin='lower')
cbar2 = plt.colorbar(cax2)
cbar2.set_label(r'\textbf{Color Scale}')
fig2.savefig(r'weightmap_FLAT_R.png',dpi = 100)

data3, header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/weightmap_FLAT_V.FIT', header=True, ext=0)

fig3, ax3 = plt.subplots(figsize=(16, 12))
# Add colorbar with custom ticks
ax3.set_title(r'\textbf{Weight map FLAT V}')
cax3 =ax3.imshow(data3,cmap=cm,origin='lower')
cbar3 = plt.colorbar(cax3)
cbar3.set_label(r'\textbf{Color Scale}')
fig3.savefig(r'weightmap_FLAT_V.png',dpi = 100)
