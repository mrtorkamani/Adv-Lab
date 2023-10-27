#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 11:06:24 2023

@author: mr
"""

import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import scipy.signal as sig
from matplotlib.pyplot import cm

plt.style.use('custom_style.mplstyle')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

### calibration from previous task
c = 479.9449423815623  ### in MHz s^-1
c_err =  0.24163352321285145 ### in MHz s^-1
k = 1.380649 * 10**-23
c_light = 299792458
u = 1.66053906660 * 10**-27

data5 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 3/Exercise3 CH2.csv', skiprows = 21, delimiter = ',')
data7 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Beer Law/7 cm Vapor Cell/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data2 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Beer Law/2 cm Vapor cell/tek0001CH2.csv', skiprows = 21, delimiter = ',')


def background(x, a, b, c, d, e,f):
    return a*x**5 + b*x**4 +c*x**3 +d*x**2 +e *x + f                    ## fitparameters and definition of functions, used 

def gaussian(v, v_d, v_0, amp):
    return amp*np.exp(-(v - v_0)**2/(2*v_d**2))

def fitfunc(x, a, b,c,d,e,f , a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, e1, e2, e3, f1, f2, f3):
    return background(x, a, b, c,d,e,f)+ gaussian(x, a1, a2, a3) + gaussian(x, b1, b2, b3) + gaussian(x, c1, c2, c3) + gaussian(x, d1, d2, d3) + gaussian(x, e1, e2, e3) + gaussian(x, f1, f2, f3)
    


time = data2[:,0]

mask = ((time > -0.0145) & (time < 0.009))
#time=time[mask]
#T =T[mask]
#frequency = c*np.array(time)-c*time[0]

color=iter(cm.Dark2(np.linspace(0,1,8)))
color2=iter(cm.Accent(np.linspace(0,1,4)))

coll = ['lightgray', 'lightpink', 'lightblue']

initvals = [0,0,0,0,0,0
            ,-0.26, 1.900,-.26
            ,-0.25, 2.300,-.25
            ,-0.325,3.200,-325
            ,-0.3,6.000,-.3
            ,-0.16,7.700,-.16
            ,-0.2,8.200,-.2]
fig, ax = plt.subplots(figsize = (10,6))
length = np.array([5, 7, 2])
v = np.linspace(0,11,10000)

amplitudes = np.zeros((3,6))
amps_errs = np.zeros((3,6))
z=0
for dataset in data5, data7, data2:
    col=next(color)
    time = dataset[mask,0]
    frequency = c*time -c*time[0]
    T = dataset[mask,1]
    p, perr = so.curve_fit(fitfunc, frequency, T, p0=initvals, maxfev = 15000)
    err = np.sqrt(np.diag(perr))
    ax.plot(v, fitfunc(v, *p),  zorder=2, color =col, label = '%i cm cell' %(length[z]))
    col2=next(color2)
    ax.errorbar(frequency, T, fmt = '.', zorder=1, markersize = 1, color = coll[z])

    print('For the %i cm vapor cell' %length[z])
    for i in range(6):
        print('line %i & $%.5f \pm %.5f$ & $%.4f \pm %.4f$ & $%.5f \pm %.5f$ //' %((i+1), p[3*i+8], err[3*i+8], p[3*i+7], err[3*i+7], -p[3*i+8], err[3*i+8]))
        amplitudes[z,i] = p[3*i+8]
        amps_errs[z,i] = err[3*i+8]
    z=z+1

ax.legend()
ax.set_xlabel(r'$\nu- \nu_0$ $\mathrm{[GHz]}$ ')
ax.set_ylabel('Transmission  [arb. units]')
fig.savefig('All_spectra.png')

def lin_func(x, a, b):
    return a*x+b
color=iter(cm.Dark2(np.linspace(0,1,6)))

fig, ax = plt.subplots(figsize = (10,6))
for i in range(6):
    p, perr = so.curve_fit(lin_func, length, amplitudes[:,i], sigma = amps_errs[:,i])
    err = np.sqrt(np.diag(perr))
    c=next(color)
    ax.errorbar(length, amplitudes[:,i], yerr = amps_errs[:,i], label = 'line %i' %(i+1), fmt = 'o', color = c)
    ax.plot(length, lin_func(length, *p), color = c)
    ax.legend()
ax.set_xlabel(r'length $\mathrm{[cm]}$')
ax.set_ylabel('Absorption depth [arb. units]')
fig.savefig('Lambeert_Beer.png')

