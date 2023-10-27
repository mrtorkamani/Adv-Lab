#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:17:30 2023

@author: mr
"""

import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.optimize import curve_fit
from matplotlib.pyplot import cm

#plt.style.use('custom_style.mplstyle')

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

c = 402.6176154672397  ### in MHz s^-1
c_err =  0.09999531103017982 ### in MHz s^-1
k = 1.380649 * 10**-23
c_light = 299792458
u = 1.66053906660 * 10**-27

data  = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data1 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.037/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data2 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.321/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data3 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.512/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data4 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.955/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data5 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.037 + 0.321/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data6 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.321 + 0.512/tek0000CH2.csv', skiprows = 21, delimiter = ',')
data7 = np.loadtxt('/home/mr/Desktop/Lab/A246/A264/Exercise 4/Saturation/0.321 + 0.955/tek0000CH2.csv', skiprows = 21, delimiter = ',')





t_1, c3_1 = data[:,0],data[:,1]
mask = ((t_1 > -0.0145) & (t_1 < +0.009))
t_1=t_1[mask]
c3_1 =c3_1[mask]
f_1 = c*np.array(t_1)-c*t_1[0]

t_2, c3_2 = data1[:,0],data1[:,1]
mask = ((t_2 > -0.0145) & (t_2 < +0.009))
t_2=t_2[mask]
c3_2 =c3_2[mask]
f_2 = c*np.array(t_2)-c*t_2[0]

t_3, c3_3 = data2[:,0],data2[:,1]
mask = ((t_3 > -0.0145) & (t_3 < +0.009))
t_3=t_3[mask]
c3_3 =c3_3[mask]
f_3 = c*np.array(t_3)-c*t_3[0]

t_4, c3_4 = data3[:,0],data3[:,1]
mask = ((t_4 > -0.0145) & (t_4 < +0.009))
t_4=t_4[mask]
c3_4 =c3_4[mask]
f_4 = c*np.array(t_4)-c*t_4[0]

t_5, c3_5 = data4[:,0],data4[:,1]
mask = ((t_5 > -0.0145) & (t_5 < +0.009))
t_5=t_5[mask]
c3_5 =c3_5[mask]
f_5 = c*np.array(t_5)-c*t_5[0]

t_6, c3_6 = data5[:,0],data5[:,1]
mask = ((t_6 > -0.0145) & (t_6 < +0.009))
t_6=t_6[mask]
c3_6 =c3_6[mask]
f_6 = c*np.array(t_6)-c*t_6[0]

t_7, c3_7 = data6[:,0],data6[:,1]
mask = ((t_7 > -0.0145) & (t_7 < +0.009))
t_7=t_7[mask]
c3_7 =c3_7[mask]
f_7 = c*np.array(t_7)-c*t_7[0]

t_8, c3_8 = data7[:,0],data7[:,1]
mask = ((t_8 > -0.0145) & (t_8 < +0.009))
t_8=t_8[mask]
c3_8 =c3_8[mask]
f_8 = c*np.array(t_8)-c*t_8[0]

#Doppler Lorentzian fit
def doppler(v, FWHM, v0, A, B, C):
    return A * (FWHM/2)**2 / ((v-v0)**2+(FWHM/2)**2) + B*v + C



initial =[.01,1.44,0,0,0]
p_1, pcov_1 = curve_fit(doppler, f_1[3420:3840], c3_1[3420:3840],p0=initial, maxfev = 15000)
yfit_1 = []
for i in f_1[3420:3840]:
    yfit_1.append(p_1[2] * (p_1[0]/2)**2 / ((i-p_1[1])**2+(p_1[0]/2)**2) + p_1[3]*i + p_1[4])
perr_1 = np.sqrt(np.diag(pcov_1))

p_2, pcov_2 = curve_fit(doppler, f_2[3420:3760], c3_2[3420:3760],p0=initial, maxfev = 15000)
yfit_2 = []
for i in f_2[3420:3760]:
    yfit_2.append(p_2[2] * (p_2[0]/2)**2 / ((i-p_2[1])**2+(p_2[0]/2)**2) + p_2[3]*i + p_2[4])
perr_2 = np.sqrt(np.diag(pcov_2))

p_3, pcov_3 = curve_fit(doppler, f_3[3420:3800], c3_3[3420:3800],p0=initial, maxfev = 15000)
yfit_3 = []
for i in f_3[3420:3800]:
    yfit_3.append(p_3[2] * (p_3[0]/2)**2 / ((i-p_3[1])**2+(p_3[0]/2)**2) + p_3[3]*i + p_3[4])
perr_3 = np.sqrt(np.diag(pcov_3))

p_4, pcov_4 = curve_fit(doppler, f_4[3420:3840], c3_4[3420:3840],p0=initial, maxfev = 15000)
yfit_4 = []
for i in f_4[3420:3840]:
    yfit_4.append(p_4[2] * (p_4[0]/2)**2 / ((i-p_4[1])**2+(p_4[0]/2)**2) + p_4[3]*i + p_4[4])
perr_4 = np.sqrt(np.diag(pcov_4))

p_5, pcov_5 = curve_fit(doppler, f_5[3420:3800], c3_5[3420:3800],p0=initial, maxfev = 15000)
yfit_5 = []
for i in f_5[3420:3800]:
    yfit_5.append(p_5[2] * (p_5[0]/2)**2 / ((i-p_5[1])**2+(p_5[0]/2)**2) + p_5[3]*i + p_5[4])
perr_5 = np.sqrt(np.diag(pcov_5))

p_6, pcov_6 = curve_fit(doppler, f_6[3420:3700], c3_6[3420:3700],p0=initial, maxfev = 15000)
yfit_6 = []
for i in f_6[3420:3700]:
    yfit_6.append(p_6[2] * (p_6[0]/2)**2 / ((i-p_6[1])**2+(p_6[0]/2)**2) + p_6[3]*i + p_6[4])
perr_6 = np.sqrt(np.diag(pcov_6))

p_7, pcov_7 = curve_fit(doppler, f_7[3420:3800], c3_7[3420:3800],p0=initial, maxfev = 15000)
yfit_7 = []
for i in f_7[3420:3800]:
    yfit_7.append(p_7[2] * (p_7[0]/2)**2 / ((i-p_7[1])**2+(p_7[0]/2)**2) + p_7[3]*i + p_7[4])
perr_7 = np.sqrt(np.diag(pcov_7))

p_8, pcov_8 = curve_fit(doppler, f_8[3420:3800], c3_8[3420:3800],p0=initial, maxfev = 15000)
yfit_8 = []
for i in f_7[3420:3800]:
    yfit_8.append(p_8[2] * (p_8[0]/2)**2 / ((i-p_8[1])**2+(p_8[0]/2)**2) + p_8[3]*i + p_8[4])
perr_8 = np.sqrt(np.diag(pcov_8))




Attenuator = [0.000, 0.037, 0.321, 0.512, 0.955, 0.321*0.037, 0.512*0.321, 0.321*0.955]#Attenuators
#Plotting:
fig3, axes3 = plt.subplots(2,4,figsize=(15, 6))



axes3[0][0].errorbar(f_1[3420:3840], c3_1[3420:3840], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[0][0].plot(f_1[3420:3840], yfit_1, zorder=2, color = 'red', label='Fitted curve')
axes3[0][0].legend()
axes3[0][0].set_xlabel('Frequnecy [GHz]')
axes3[0][0].set_ylabel('Transmission [arb. units]')
axes3[0][0].set_title(f'''Att. {Attenuator[0]:0.3f}''')

axes3[0][1].errorbar(f_2[3420:3760], c3_2[3420:3760], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[0][1].plot(f_2[3420:3760], yfit_2, zorder=2, color = 'red', label='Fitted curve')
axes3[0][1].legend()
axes3[0][1].set_xlabel('Frequnecy [GHz]')
axes3[0][1].set_ylabel('Transmission [arb. units]')
axes3[0][1].set_title(f'''Att. {Attenuator[1]:0.3f}''')

axes3[0][2].errorbar(f_3[3420:3800], c3_3[3420:3800], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[0][2].plot(f_3[3420:3800], yfit_3, zorder=2, color = 'red', label='Fitted curve')
axes3[0][2].legend()
axes3[0][2].set_xlabel('Frequnecy [GHz]')
axes3[0][2].set_ylabel('Transmission [arb. units]')
axes3[0][2].set_title(f'''Att. {Attenuator[2]:0.3f}''')

axes3[0][3].errorbar(f_4[3420:3840], c3_4[3420:3840], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[0][3].plot(f_4[3420:3840], yfit_4, zorder=2, color = 'red', label='Fitted curve')
axes3[0][3].legend()
axes3[0][3].set_xlabel('Frequnecy [GHz]')
axes3[0][3].set_ylabel('Transmission [arb. units]')
axes3[0][3].set_title(f'''Att. {Attenuator[3]:0.3f}''')

axes3[1][0].errorbar(f_5[3420:3800], c3_5[3420:3800], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[1][0].plot(f_5[3420:3800], yfit_5, zorder=2, color = 'red', label='Fitted curve')
axes3[1][0].legend()
axes3[1][0].set_xlabel('Frequnecy [GHz]')
axes3[1][0].set_ylabel('Transmission [arb. units]')
axes3[1][0].set_title(f'''Att. {Attenuator[4]:0.3f}''')

axes3[1][1].errorbar(f_6[3420:3700], c3_6[3420:3700], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[1][1].plot(f_6[3420:3700], yfit_6, zorder=2, color = 'red', label='Fitted curve')
axes3[1][1].legend()
axes3[1][1].set_xlabel('Frequnecy [GHz]')
axes3[1][1].set_ylabel('Transmission [arb. units]')
axes3[1][1].set_title(f'''Att. {Attenuator[5]:0.3f}''')

axes3[1][2].errorbar(f_7[3420:3800], c3_7[3420:3800], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[1][2].plot(f_7[3420:3800], yfit_7, zorder=2, color = 'red',label='Fitted curve')
axes3[1][2].legend()
axes3[1][2].set_xlabel('Frequnecy [GHz]')
axes3[1][2].set_ylabel('Transmission [arb. units]')
axes3[1][2].set_title(f'''Att. {Attenuator[6]:0.3f}''')

axes3[1][3].errorbar(f_8[3420:3800], c3_8[3420:3800], color = 'black', fmt = '.', zorder=1, markersize = 1)
axes3[1][3].plot(f_8[3420:3800], yfit_8, zorder=2, color = 'red',label='Fitted curve')
axes3[1][3].legend()
axes3[1][3].set_xlabel('Frequnecy [GHz]')
axes3[1][3].set_ylabel('Transmission [arb. units]')
axes3[1][3].set_title(f'''Att. {Attenuator[7]:0.3f}''')

fig3.tight_layout()
fig3.savefig("saturation_fit.png")

def tryfit(I, I0, a):
    return(a*(I-I0))

I = 115 #mA 
P = [tryfit(I, 54, 10.41)]
dP = [np.sqrt(((I-54)*0.1)**2 + 2*(10.41*0.1)**2)] #Gaussian error propagation
for i in Attenuator[1:]:
    P.append(i*P[0])
    dP.append(i*dP[0])

print('Attenuator & FWHM &  nu_0 & A & B & C  ')
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[0], p_1[0], perr_1[0], p_1[1], perr_1[1], p_1[2], perr_1[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[1], p_2[0], perr_2[0], p_2[1], perr_2[1], p_2[2], perr_2[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[2], p_3[0], perr_3[0], p_3[1], perr_3[1], p_3[2], perr_3[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[3], p_4[0], perr_4[0], p_4[1], perr_4[1], p_4[2], perr_4[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[4], p_5[0], perr_5[0], p_5[1], perr_5[1], p_5[2], perr_5[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[5], p_6[0], perr_6[0], p_6[1], perr_6[1], p_6[2], perr_6[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[6], p_7[0], perr_7[0], p_7[1], perr_7[1], p_7[2], perr_7[2]))
print('$%.4f$ & $%.4f \pm %.4f$& $%.4f \pm %.4f$& $%.4f \pm %.4f$ \ \ ' %(Attenuator[7], p_8[0], perr_8[0], p_8[1], perr_8[1], p_8[2], perr_8[2]))


A = np.array([p_1[2], p_2[2], p_3[2], p_4[2], p_5[2], p_6[2], p_7[2],p_8[2]])
dA = np.array([perr_1[2], perr_2[2], perr_3[2], perr_4[2], perr_5[2], perr_6[2], perr_7[2],perr_8[2]])

v2 = np.array([p_1[0], p_2[0], p_3[0], p_4[0], p_5[0], p_6[0], p_7[0], p_8[0]])*4
dv2 = np.array([perr_1[0], perr_2[0], perr_3[0], perr_4[0], perr_5[0], perr_6[0], perr_7[0], perr_8[0]])*4

P = np.sort(P)
dP = np.sort(dP)
A = np.sort(A)
dA = np.sort(dA)
v2 = np.sort(v2)
dv2 = np.sort(dv2)

def amp(P, Psat, C1, C2):
    return(C1 * P/Psat/(np.sqrt(1+P/Psat)) + C2)

def dnusquared(P, Psat, dnu):
    return(dnu**2 * (1+P/Psat))

LLL = np.linspace(0,700,1200)
p_10, pcov_10 = curve_fit(amp, np.array(P), A, maxfev = 15000)
yfit_10 = amp(LLL, *p_10)
#yfit_7 = []
#for i in f_7[820:890]:
#    yfit_7.append(p_7[2] * (p_7[0]/2)**2 / ((i-p_7[1])**2+(p_7[0]/2)**2) + p_7[3]*i + p_7[4])
perr_10 = np.sqrt(np.diag(pcov_10))

p_11, pcov_11 = curve_fit(dnusquared, np.array(P), v2**2, maxfev = 15000)
yfit_11 = dnusquared(LLL, *p_11)
#yfit_7 = []
#for i in f_7[820:890]:
#    yfit_7.append(p_7[2] * (p_7[0]/2)**2 / ((i-p_7[1])**2+(p_7[0]/2)**2) + p_7[3]*i + p_7[4])
perr_11 = np.sqrt(np.diag(pcov_11))

#Plotting:

fig4, axes4 = plt.subplots(1,2,figsize=(15, 6))
axes4[0].errorbar(P, A, xerr=dP, yerr=dA, label='Data', color = 'black', fmt = '.', zorder=1, markersize = 1)
axes4[0].errorbar(LLL, yfit_10, zorder=2, color = 'blue', label='Fitted curve')
axes4[0].set_xlabel('laser power [$\mu$W]')
axes4[0].set_ylabel('amplitude [arb. units]')
axes4[0].legend(loc='upper left')
#axes4[0].set_ylim(-0.5, 0.2)

axes4[1].errorbar(P, v2**2, xerr=dP, yerr=dv2**2, label='Data', color = 'black', fmt = '.', zorder=1, markersize = 1)
axes4[1].errorbar(LLL, yfit_11, zorder=2, color = 'blue', label='Fitted curve')
axes4[1].set_ylabel(f'''line width $\Delta$v$^2$ [GHz]''')
axes4[1].set_xlabel('laser power [$\mu$W]')
axes4[1].legend()
#axes4[1].set_ylim(0.0, 0.004)
fig4.savefig("saturation-power.png")

print(fr'p_sat_amp = {p_10[0]}+-{perr_10[0]/100}')
print(fr'p_sat_linewidth = {np.sqrt(p_11[0])}+-{np.sqrt(perr_11[0]/100)}')






