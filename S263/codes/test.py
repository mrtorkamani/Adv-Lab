#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 14:38:22 2023

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


data , header = fits.getdata('/home/mr/Desktop/Lab/S263/codes/produced_data/Astrometric_Calibration_B.fits' , header = True , ext = 0)

print(header)
