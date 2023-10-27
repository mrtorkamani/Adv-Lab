#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 21:24:32 2023

@author: mr
"""

import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import scipy.signal as sig

plt.style.use('custom_style.mplstyle')

data = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 3/Excericse 3 CH3.csv', skiprows = 21, delimiter = ',')  ### read in data
#data = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 3/Exercise3 CH2.csv', skiprows = 21, delimiter = ',')  ### read in data
time = data[:,0]
T = data[:,1]

mask = ((time > -0.0145) & (time < -0.0126))
time=time[mask]
T =T[mask]

#plt.plot(time,T)
#plt.xlim(-.014,-0.009)

#print(max(T))

def fitfunc(x, a, b, c, d, e):                   ## fitfunction as transmission of FPI
    return a / (1+b*np.sin(c*x+d)**2) + e
initvals = [1,1,9000,1,11]
p, pcov = so.curve_fit(fitfunc, time, T, p0= initvals)       ### fitted function
parerrors = np.sqrt(np.diag(pcov))

peaks = sig.find_peaks(T, height = 0.07, distance = 500)     ### find peaks of spectrum
a = np.zeros(5)
for i in range(len(peaks[0])-1):
    a[i] = time[peaks[0][i+1]]-time[peaks[0][i]]

mean_time = np.mean(a)              ### mean time in between peaks
mean_time_var = np.var(a)

FSR = 2 * mean_time
FSR_err = 2*mean_time_var
print('The mean separation of peaks is', mean_time*1000,'+-', mean_time_var*1000,'ms')
print('The t_FSR', FSR*1000,'+-', FSR_err*1000,'ms')


dt = 2 /np.sqrt(p[1]*p[2]**2)         ### FWHM calculation
dt_err = np.sqrt((dt/(2*p[1])*parerrors[1])**2 + (dt/p[2] * parerrors[2])**2)



print('FWHM is', dt*1000,'+-',  dt_err*1000, 'ms')

F = FSR/dt               ### Finesse calculation
F_err = np.sqrt((F/FSR *FSR_err)**2 + (F/dt*dt_err)**2)
print('The calculated Finesse is: F =', F,'+-', F_err)


fig, ax = plt.subplots(figsize = (10, 6))  ### plot of spectrum and fit
ax.errorbar(time*1000, T, fmt ='+', markersize = 1, linestyle = '', color = 'black', label = 'data')
ax.plot(time*1000, fitfunc(time*1000, p[0], p[1], p[2]/1000, p[3], p[4]), color = 'red', label = 'fitted function')
ax.set_xlabel(r'time  $[\mathrm{ms}]$')
ax.set_ylabel(r'voltage [arb. units]')
ax.set_ylim(-.15,.2)
ax.legend()
ax.grid()
fig.savefig('FPI.png')
print('fitted parameters')
print('a =',p[0],'±',parerrors[0])
print('b =',p[1],'±',parerrors[1])
print('c =',p[2],'±',parerrors[2])
print('d =',p[3],'±',parerrors[3])
print('e =',p[4],'±',parerrors[4])

### calibration


Dnu = 149.9348
Dnu_err = 0.0002
mean_time = mean_time
mean_time_err = mean_time_var


c = Dnu / mean_time
c_err = np.sqrt((c/Dnu*Dnu_err)**2 + (c/mean_time*mean_time_err)**2)

print('Calibration constant is c =',c/1000,'+-', c_err/1000,'MHz ms^-1')

