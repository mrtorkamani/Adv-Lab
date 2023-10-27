#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 12:27:55 2023

@author: mr
"""

import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm

plt.style.use('custom_style.mplstyle')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})
### calibration from previous task
c = 402.6176154672397  ### in MHz s^-1
c_err =  0.24163352321285145 ### in MHz s^-1
k = 1.380649 * 10**-23
c_light = 299792458
u = 1.66053906660 * 10**-27

data = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0/tek0000CH2.csv', skiprows = 21, delimiter = ',')

time = data[:,0]
T = data[:,1]

mask = ((time > -0.0145) & (time < +0.009))
time=time[mask]
ch3 =T[mask]
frequency = c*time -c*time[0]



#plt.plot(freqb,ch3b)

#print(frequency[5080:5240])
def background(x, a, b, c, d, e,f):
    return a*x**5 + b*x**4 +c*x**3 +d*x**2 +e *x + f  
#Doppler Lorentzian fit
def doppler(v, FWHM, v0, A, B, C):
    return A * (FWHM/2)**2 / ((v-v0)**2+(FWHM/2)**2) + B*v + C 



initials1 = [0.07, 1.45, 0.0, 0,0]
p1, pcov1 = curve_fit(doppler, frequency[3420:3840], ch3[3420:3840], p0 =initials1, maxfev = 15000)
perr1 = np.sqrt(np.diag(pcov1))

initials2 = [0.01, 2.1, -0.00, 0,0]
p2, pcov2 = curve_fit(doppler, frequency[5000:5380], ch3[5000:5380], p0 =initials2, maxfev = 15000)
perr2 = np.sqrt(np.diag(pcov2))

initials3 = [0.01, 2.63, -0.00, 0,0]
p3, pcov3 = curve_fit(doppler, frequency[6330:6740], ch3[6330:6740], p0 =initials3, maxfev = 15000)
perr3 = np.sqrt(np.diag(pcov3))

initials4 = [0.01, 2.93, 0, 10,0]
p4, pcov4 = curve_fit(doppler, frequency[7020:7750], ch3[7020:7750], p0 =initials4, maxfev = 15000)
perr4 = np.sqrt(np.diag(pcov4))

initials5 = [0.01, 4.86, 0, 0,0]
p5, pcov5 = curve_fit(doppler, frequency[11900:12200], ch3[11900:12200], p0 =initials5, maxfev = 15000)
perr5 = np.sqrt(np.diag(pcov5))

initials6 = [0.06, 5.15, 1.0, 0,0]
p6, pcov6 = curve_fit(doppler, frequency[12500:13100], ch3[12500:13100], p0 =initials6, maxfev = 15000)
perr6 = np.sqrt(np.diag(pcov6))

#initials7 = [2.56, 6.60, 0.0, 0,0]
#p7, pcov7 = curve_fit(doppler, frequency[16110:16600], ch3[16110:16600], p0 =initials7, maxfev = 15000)
#perr7 = np.sqrt(np.diag(pcov7))

initials8 = [.06, 7.05, 0.0, 0,0]
p8, pcov8 = curve_fit(doppler, frequency[17200:17800], ch3[17200:17800], p0 =initials8, maxfev = 15000)
perr8 = np.sqrt(np.diag(pcov8))
fig, ax = plt.subplots(figsize=(9, 9))
ax.errorbar(frequency, ch3, fmt ='+', markersize = 1, linestyle = '', color = 'gray', label = 'data')
'''
v1 = np.linspace(1.40,1.50,1000)
ax.plot(v1,doppler(v1,*p1),color='blue')

v2 = np.linspace(2.05,2.12,1000)
ax.plot(v2,doppler(v2,*p2),color='blue',label = 'Fitted function')

v3 = np.linspace(2.55,2.7,1000)
ax.plot(v3,doppler(v3,*p3),color='blue')

v4 = np.linspace(2.84,3.1,1000)
ax.plot(v4,doppler(v4,*p4),color='blue')

v5 = np.linspace(4.8,4.90,1000)
ax.plot(v5,doppler(v5,*p5),color='blue')

v6 = np.linspace(5.05,5.25,1000)
ax.plot(v6,doppler(v6,*p6),color='blue')

#v7 = np.linspace(6.5,6.75,1000)
#ax.plot(v7,doppler(v7,*p7),color='blue')
v8 = np.linspace(6.99,7.150,1000)
plt.plot(v8,doppler(v8,*p8),color='blue')
'''
ax.legend()
ax.grid()
ax.set_xlabel(r'$Frequency$ $\mathrm{[GHz]}$ ')
ax.set_ylabel('Transmission  [arb. units]')

fig.savefig('dip_analysis.png')

#FWHM = [p1[0],p2[0],p3[0],p4[0],p5[0],p6[0],p7[0],p8[0]]
#FWHM_err =[perr1[0],perr2[0],perr2[0],perr2[0],perr2[0],perr2[0],perr2[0],perr2[0]]

#positon = [p1[1],p2[1],p3[1],p4[1],p5[1],p6[1],p7[1],p8[1]]
#positon_err = [perr1[1],perr2[1],perr3[1],perr4[1],perr5[1],perr6[1],perr7[1],perr8[1]]


### HFS splitting
Delta_nu_S85 = p6[1]-p3[1]
Delta_nu_P85 = p5[1]-p4[1]
Delta_nu_S87 = 6.6-p1[1]
Delta_nu_P87 = p8[1]-p2[1]

DS85_err = np.sqrt((perr6[1]**2+perr4[1]**2))
DP85_err = np.sqrt((perr5[1]**2+perr3[1]**2))
DS87_err = np.sqrt((.2**2+perr1[1]**2))
DP87_err = np.sqrt((perr8[1]**2+perr2[0]**2))
                  

A_S85 = Delta_nu_S85/3
A_S85_err = DS85_err/3

A_P85 = Delta_nu_P85/3
A_P85_err = DP85_err/3

A_S87 = Delta_nu_S87/2
A_S87_err = DS87_err/2

A_P87 = Delta_nu_P87/2
A_P87_err = DP87_err/2

print('The hyperfine coupling constant for the ground state of 85 Rb is h*', A_S85*4,'+-' , A_S85_err*4,'GHz')
print('The hyperfine coupling constant for the excited state of 85 Rb is h*', A_P85/2,'+-' , A_P85_err/2,'GHz')
print('The hyperfine coupling constant for the ground state of 87 Rb is h*', A_S87*4,'+-' , A_S87_err*4,'GHz')
print('The hyperfine coupling constant for the excited state of 87 Rb is h*', A_P87/2,'+-' , A_P87_err/2,'GHz')

c0 = 299792458 #m/s
l0 = 794e-9
v_85 = c0/l0 # MHz

#Spectral resolution:
#v_85 = 377.107385690 * 10 ** 12

FWHM = np.array([p1[0], p2[0], p3[0], p4[0], p5[0], p6[0],p8[0]])*1000
dFWHM = np.array([perr1[0], perr2[0], perr3[0], perr4[0], perr5[0], perr6[0], perr8[0]])*1000

spectral_resolution = []
dspectral_resolution = []
for i in FWHM:
    spectral_resolution.append(i*10**6/v_85)
for i in dFWHM:
    dspectral_resolution.append(i*10**6/v_85)
    
Rval = np.mean(spectral_resolution)
R_err = np.std(spectral_resolution)

print(fr'leaser frequency = {v_85}')
print('The spectral resolution is R = %.9f+-%.9f'%(Rval, R_err))