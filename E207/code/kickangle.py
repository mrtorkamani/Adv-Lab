#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 09:01:42 2023

@author: mr
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as so


plt.style.use('custom_style.mplstyle')


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

def linear(x, m, b):
    return m*x+b

def chisqure(a,b):
    # no of hours a student studies
    # in a week vs expected no of hours
    observed_data = a
    expected_data = b
      
      
    # determining chi square goodness of fit using formula
    chi_square_test_statistic1 = 0
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            np.square((observed_data[i]-expected_data[i]))/expected_data[i]
    return chi_square_test_statistic1



C0 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Kick Angle Calibration/C0.csv')
C1 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Kick Angle Calibration/C1.csv')
errC0=[.01,.01,.01,.01,.01,.01,.01,.01,.01,.01]
errC1x = [.05,.05,.05,.05,.05,.05,.05,.05,.05]
errC1y = [.1,.1,.1,.1,.1]

'''
a = 305
poptx, pcovx = so.curve_fit(linear, a*C0['Ix'], C0['x'], sigma=errC0)
perrx = np.sqrt(np.diag(pcovx))

b = 305
popty, pcovy = so.curve_fit(linear, b*C0['Iy'], C0['y'], sigma=errC0)
perry = np.sqrt(np.diag(pcovy))

c = 305
poptx1, pcovx1 = so.curve_fit(linear, c*C1['Ix'], C1['x'], sigma=errC1x)
perrx1 = np.sqrt(np.diag(pcovx1))

d = 305
popty1, pcovy1 = so.curve_fit(linear, d*C1['Iy'][:5], C1['y'][:5], sigma=errC1y)
perry1 = np.sqrt(np.diag(pcovy1))




fig, axs = plt.subplots(2, 2,figsize=(9, 9))
axs[0, 0].errorbar(a*C0['Ix'],C0['x'],yerr= errC0,fmt=".")
axs[0, 0].plot(a*C0['Ix'],linear(a*C0['Ix'], *poptx),label = fr"$\chi^2$ ={np.abs(chisqure(C0['x'], linear(a*C0['Ix'], *poptx))):0.2}")
axs[0, 0].legend()
axs[0, 0].set_title(r'Corrector 0 X-direction')

axs[0, 1].errorbar(b*C0['Iy'],C0['y'],yerr= errC0,fmt=".")
axs[0, 1].plot(b*C0['Iy'],linear(b*C0['Iy'], *popty),label = fr"$\chi^2$ ={np.abs(chisqure(C0['y'], linear(b*C0['Iy'], *popty))):0.2}")
axs[0, 1].legend()
axs[0, 1].set_title(r'Corrector 0 Z-direction')

axs[1, 0].errorbar(c*C1['Ix'],C1['x'],yerr= errC1x,fmt=".")
axs[1, 0].plot(c*C1['Ix'],linear(c*C1['Ix'], *poptx1),label = fr"$\chi^2$ ={np.abs(chisqure(C1['x'], linear(c*C1['Ix'], *poptx1))):0.2}")
axs[1, 0].legend()
axs[1, 0].set_title(r'Corrector 1 X-direction')

axs[1, 1].errorbar(d*C1['Iy'][:5],C1['y'][:5],yerr= errC1y,fmt=".")
axs[1, 1].plot(d*C1['Iy'],linear(d*C1['Iy'], *popty1),label = fr"$\chi^2$ ={np.abs(chisqure(C1['y'][:5], linear(d*C1['Iy'][:5], *popty1))):0.2}")
axs[1, 1].legend()
axs[1, 1].set_title(r'Corrector 1 Z-direction')

for i in fig.get_axes():
    i.set_xlabel(r'$\alpha \mathrm{[mrad]}$')
    i.set_ylabel(r'Offset displacement $\mathrm{[mm]}$')

fig.tight_layout()

#fig.savefig('displacement_by_angle_kickangle.png')

print('Corrector & Slope X-direction $\mathrm{[m]}$ & Slope Z-direction $\mathrm{[m]}$ & L $\mathrm{[m]}$ \\\ ')
print('C0 & $%.3f \pm %.3f$ & $%.3f \pm %.3f$ & $0.24425$ \\\ ' %(poptx[0] ,perrx[0] , popty[0],perry[0]))
print('C1 & $%.3f \pm %.3f$ & $%.2f \pm %.2f$ & $0.37525$ \\\ ' %(poptx1[0] ,perrx1[0] , popty1[0],perry1[0]))

L0 = 0.24425
L1 = 0.37525

fig, axs = plt.subplots(2, 2,figsize=(9, 9))
axs[0, 0].errorbar(a*C0['Ix'],C0['x']/L0,yerr= errC0,fmt=".")
axs[0, 0].plot(a*C0['Ix'],linear(a*C0['Ix'], *poptx)/L0,label = fr"$\chi^2$ ={np.abs(chisqure(C0['x']/L0, linear(a*C0['Ix'], *poptx)/L0)):0.2}")
axs[0, 0].legend()
axs[0, 0].set_title(r'Corrector 0 X-direction')

axs[0, 1].errorbar(b*C0['Iy'],C0['y']/L0,yerr= errC0,fmt=".")
axs[0, 1].plot(b*C0['Iy'],linear(b*C0['Iy'], *popty)/L0,label = fr"$\chi^2$ ={np.abs(chisqure(C0['y']/L0, linear(b*C0['Iy'], *popty)/L0)):0.2}")
axs[0, 1].legend()
axs[0, 1].set_title(r'Corrector 0 Z-direction')

axs[1, 0].errorbar(c*C1['Ix'],C1['x']/L1,yerr= errC1x,fmt=".")
axs[1, 0].plot(c*C1['Ix'],linear(c*C1['Ix'], *poptx1)/L1,label = fr"$\chi^2$ ={np.abs(chisqure(C1['x']/L1, linear(c*C1['Ix'], *poptx1)/L1)):0.2}")
axs[1, 0].legend()
axs[1, 0].set_title(r'Corrector 1 X-direction')

axs[1, 1].errorbar(d*C1['Iy'][:5],C1['y'][:5]/L1,yerr= errC1y,fmt=".")
axs[1, 1].plot(d*C1['Iy'],linear(d*C1['Iy'], *popty1)/L1,label = fr"$\chi^2$ ={np.abs(chisqure(C1['y'][:5]/L1, linear(d*C1['Iy'][:5], *popty1)/L1)):0.2}")
axs[1, 1].legend()
axs[1, 1].set_title(r'Corrector 1 Z-direction')

for i in fig.get_axes():
    i.set_xlabel(r'$\alpha_{obs} \mathrm{[mrad]}$')
    i.set_ylabel(r'$\alpha_{exp} \mathrm{[mrad]}$')

fig.tight_layout()

#fig.savefig('angle_by_angle_kickangle.png')

#print('Corrector & Slope X-direction $\mathrm{[m]}$ & Slope Z-direction $\mathrm{[m]}$ \\\ ')
#print('C0 & $%.3f \pm %.3f$ & $%.3f \pm %.3f$  \\\ ' %(poptx[0]/L0 ,perrx[0]/L0 , popty[0]/L0,perry[0]/L0))
#print('C1 & $%.2f \pm %.2f$ & $%.2f \pm %.2f$ \\\ ' %(poptx1[0]/L1 ,perrx1[0]/L1 , popty1[0]/L1,perry1[0]/L1))
'''
for i in range(len(C1['x'])):
    print('$',C1['Ix'][i],'$ & $',round(305*C1['Ix'][i],1),'$ & $',C1['x'][i],'\pm', '0.01$','& $',C1['Iy'][i],'$ & $',round(305*C1['Iy'][i],1),'$ & $',C1['y'][i],'\pm', '0.01$','\\\ ')

