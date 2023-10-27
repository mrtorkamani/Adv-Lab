#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:24:36 2023

@author: mr
"""

from astropy.io import fits
from astropy.stats import sigma_clipped_stats , SigmaClip
from photutils.segmentation import detect_threshold , detect_sources
from photutils.utils import circular_footprint
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as colorss
import numpy as np



plt.style.use('custom_style.mplstyle')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

colors = [(0, 0, 0), (1, 1, 1)] # first color is black, last is white
cm = LinearSegmentedColormap.from_list(
        'custom_colormap', colors, N=4000)


data , header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced data/Co_added_V' , header = True , ext = 0)



sigma_clip = SigmaClip(sigma = 4.0 , maxiters = 11)
threshold = detect_threshold(data , nsigma = 2.0 , sigma_clip = sigma_clip)
segment_img = detect_sources(data , threshold , npixels = 5)
footprint = circular_footprint(radius = 11)
mask = segment_img.make_source_mask(footprint = footprint)
mean , median , std = sigma_clipped_stats(data , sigma = 4.0 , mask = mask)
backsub_data = data - mean

fits.writeto('/home/mr/Desktop/Lab/S263/codes/produced data/backsub_V.FIT', backsub_data , header = header, overwrite = True ,
                 output_verify = "ignore")
