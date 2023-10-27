#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 23:00:45 2023

@author: mr
"""

import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt

plt.style.use('custom_style.mplstyle')


### calibration from previous task
c = 479.9449423815623  ### in MHz s^-1
c_err =  0.24163352321285145 ### in MHz s^-1
k = 1.380649 * 10**-23
c_light = 299792458
u = 1.66053906660 * 10**-27


v_85 = 377.107385690 * 10 ** 12
v_87 = 377.1074635 * 10 ** 12

dv = v_87 - v_85

A_85_s =  1.0119108130 *10**9 
A_85_p =  120.527 *10**6

A_87_s = 3.41734130545215 * 10 ** 9
A_87_p = 408.328 * 10 ** 6

### theoretical spectra
m_85 = 84.9118 * u
m_87 = 86.9092 * u
T = 300 
def v_d(T, m, v_0):
    return v_0 * np.sqrt(k*T/(m*c_light**2))


v_d_85 = v_d(T, m_85, v_85)
v_d_87 = v_d(T, m_87, v_87)

def gaussian(v, v_d, v_0, amp):
    return amp*np.exp(-(v - v_0)**2/(2*v_d**2))
v = np.linspace(-5*10**9, 6*10**9, 1000000)

### 85
amp_33 = 5/7/3
amp_32 = 1/3
amp_23 = 1/3
amp_22 = 2/7/3

v_33 = -A_85_s*5/4 + 5/4*A_85_p
v_32 = +A_85_s*7/4 + 5/4*A_85_p
v_23 = -A_85_s*5/4 - 7/4*A_85_p
v_22 = +A_85_s*7/4 - 7/4*A_85_p



### 87

amp_22_87 = 1/3 * 27.83/72.17
amp_21 = 1/3 * 27.83/72.17
amp_12 = 1/3 * 27.83/72.17
amp_11 = 1/5/3 * 27.83/72.17

v_22_87 = -A_87_s*3/4 + 3/4*A_87_p + dv
v_21 = +A_87_s*5/4 + 3/4*A_87_p + dv
v_12 = -A_87_s*3/4 - 5/4*A_87_p + dv
v_11 = +A_87_s*5/4 - 5/4*A_87_p + dv
fig, ax = plt.subplots(figsize = (10,6))

ax.plot(v, 1 - gaussian(v, v_d_85, v_33, amp_33), label = '3 -> 3, 85')
ax.plot(v, 1 - gaussian(v, v_d_85, v_32, amp_32), label = '2 -> 3, 85')
ax.plot(v, 1 - gaussian(v, v_d_85, v_23, amp_23), label = '3 -> 2, 85')
ax.plot(v, 1 - gaussian(v, v_d_85, v_22, amp_22), label = '2 -> 2, 85')
ax.plot(v, 1 - gaussian(v, v_d_87, v_22_87, amp_22_87), label = '2 -> 2, 87')
ax.plot(v, 1 - gaussian(v, v_d_87, v_21, amp_21), label = '1 -> 2, 87')
ax.plot(v, 1 - gaussian(v, v_d_87, v_12, amp_12), label = '2 -> 1, 87')
ax.plot(v, 1 - gaussian(v, v_d_87, v_11, amp_11), label = '1 -> 1, 87')
ax.plot(v, 1- gaussian(v, v_d_85, v_33, amp_33) - gaussian(v, v_d_85, v_32, amp_32)- gaussian(v, v_d_85, v_23, amp_23) - gaussian(v, v_d_85, v_22, amp_22) - gaussian(v, v_d_87, v_22_87, amp_22_87)-
        gaussian(v, v_d_87, v_21, amp_21) - gaussian(v, v_d_87, v_12, amp_12) - gaussian(v, v_d_87, v_11, amp_11), label = 'total spectrum', color = 'black')
ax.legend()
ax.set_xlabel(r'$\nu-\nu_0$ $\mathrm{[GHz]}$')
ax.set_ylabel('Transmission')
fig.savefig('Theoretical_spectrum.png')

#data = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 3/Excericse 3 CH3.csv', skiprows = 21, delimiter = ',')  ### read in data
data = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 3/Exercise3 CH2.csv', skiprows = 21, delimiter = ',')  ### read in data
time = data[:,0]
T = data[:,1]

mask = ((time > -0.0145) & (time < 0.009))
time=time[mask]
T =T[mask]



def background(x, a, b, c, d, e,f):
    return a*x**5 + b*x**4 +c*x**3 +d*x**2 +e *x + f                    ## fitparameters and definition of functions, used 

def fitfunc(x, a, b,c,d,e,f , a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, e1, e2, e3, f1, f2, f3):
    return background(x, a, b, c,d,e,f) + gaussian(x, a1, a2, a3) + gaussian(x, b1, b2, b3) + gaussian(x, c1, c2, c3) + gaussian(x, d1, d2, d3) + gaussian(x, e1, e2, e3) + gaussian(x, f1, f2, f3)

frequency = c*np.array(time)-c*time[0]

initvals = [0,0,0,0,0,0
            ,-0.26, 1.800,-.26
            ,-0.25, 2.300,-.25
            ,-0.325,3.200,-325
            ,-0.3,6.000,-.3
            ,-0.16,7.700,-.16
            ,-0.2,8.200,-.2]
p, perr = so.curve_fit(fitfunc, frequency, T, p0=initvals, maxfev = 15000)

#plt.errorbar(frequency, T, fmt ='+', markersize = 1, linestyle = '', color = 'black', label = 'data')

v = np.linspace(0,11,10000)
fig1, ax1 = plt.subplots(figsize = (10,6))
ax1.plot(v, fitfunc(v, *p), color = 'blue', zorder=2,label = 'Fitted curve')
ax1.errorbar(frequency, T, color = 'gray', fmt = '+', zorder=1, markersize = 1,label = 'Data')
ax1.plot(v,background(v,p[0],p[1],p[2],p[3],p[4],p[5]),'--',color='red',label = 'Background')
#ax.plot(v, fitfunc(v, *initvals), color = 'black')
ax1.set_xlabel(r'$\nu- \nu_0$ $\mathrm{[GHz]}$')
ax1.set_ylabel('Transmission [arb. units]')
ax1.set_xlim(0,11)
ax1.set_ylim(-.45,-.26)
ax1.legend()
ax1.grid()
fig1.savefig('Experimental_spectrum.png')

err = np.sqrt(np.diag(perr))


for i in range(6):
    print('amp%i =' %(i+1), p[3*i+8],'+-',err[3*i+8])
    print('pos%i =' %(i+1), p[3*i+7],'+-',err[3*i+7])
    print('width%i =' %(i+1), p[3*i+6],'+-',err[3*i+6])

########background removal##########
fig2, ax2 = plt.subplots(figsize = (10,6))
TT = T - background(frequency,p[0],p[1],p[2],p[3],p[4],p[5])
ax2.plot(v, fitfunc(v, *p)-background(v,p[0],p[1],p[2],p[3],p[4],p[5]), color = 'blue', zorder=2,label = 'Fitted curve')
ax2.errorbar(frequency, TT, color = 'gray', fmt = '+', zorder=1, markersize = 1,label = 'Data')
#ax1.plot(v,background(v,p[0],p[1],p[2],p[3],p[4],p[5]),'--',color='red',label = 'Background')
#ax.plot(v, fitfunc(v, *initvals), color = 'black')
ax2.set_xlabel(r'$\nu- \nu_0$ $\mathrm{[GHz]}$')
ax2.set_ylabel('Transmission [arb. units]')
ax2.set_xlim(0,11)
ax2.set_ylim(-.15,.1)
ax2.legend()
ax2.grid()
fig2.savefig('Experimental_spectrum_backgrundremove.png')


### HFS splitting
Delta_nu_85 = p[16]-p[13]
Delta_nu_87 = np.array([p[22]-p[10],p[19]-p[7]])
D85_err = np.sqrt((err[16]**2 + err[13]**2))
D87_err = np.array([np.sqrt(err[22]**2+err[10]**2), np.sqrt(err[19]**2+ err[10]**2)])
                  
print(Delta_nu_85, Delta_nu_87)

A_85 = Delta_nu_85/3
A_85_err = D85_err/3
A_87 = np.mean(Delta_nu_87)/2
A_87_err = np.sqrt(D87_err[0]**2+D87_err[1]**2)
print('The hyperfine coupling constant for the ground state of 85 Rb is h*', A_85,'+-' , A_85_err,'GHz')
print('The hyperfine coupling constant for the ground state of 87 Rb is h*', A_87,'+-' , A_87_err,'GHz')



### FWHM

def FWHM(sigma):
    return np.abs(2.3548 * sigma)

def R(sigma):
    return FWHM(sigma) * 10**6/v_85
for i in range(6):
    print('The FWHM of line %i is %.1f +- %.1f MHz' %(i+1, FWHM(p[3*i+6])*1000,FWHM(err[3*i+6])*1000))

for i in range(6):
    print('The spectral resolution of line %i is %.9f +- %.9f' %(i+1, R(p[3*i+6])*1000,R(err[3*i+6])*1000))
    
Rval = np.mean(R(p[6::3])*1000)
R_err = np.std(R(p[6::3])*1000)

print(R(p[6::3])*1000)
print('The spectral resolution is R = %.9f+-%.9f'%(Rval, R_err))

