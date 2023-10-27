#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 18:09:19 2023

@author: mr
"""
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from IPython.display import Markdown as md
import scipy.stats as stats
import scipy.optimize as so

import os

plt.style.use('custom_style.mplstyle')


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

def linear(x, m, b):
    return m*x+b

def poly2(x, a, b, c):
    return a*x**2 + b*x + c

def xsq(x, s, a, b):
    return a*(x-s)**2 + b
def chisqure(a,b):
    # no of hours a student studies
    # in a week vs expected no of hours
    observed_data = a
    expected_data = b
      
      
    # determining chi square goodness of fit using formula
    chi_square_test_statistic1 = 0
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            (np.square(observed_data[i]-expected_data[i]))/expected_data[i]
    return chi_square_test_statistic1

##momentum calculation
epot = 511
ekin = 25
gamma = (epot + ekin)/epot
c = 299792458
velec = np.sqrt(1-(1/gamma)**2) *c 
pelec = 9.10938e-31 *gamma*velec  / 5.34428599e-28

##loading the data
#qs2 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/part2/Quadrupole Scan/s2.csv')
qs2 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/part2/Quadrupole Scan/s3.csv')

#print(qs2['k'])

I = qs2['k']
sigma_x = qs2['sigmax']
sigma_y = qs2['sigmay'] 
err = [.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2]

k = 128027 * I/pelec *10**3
errk = np.abs(150* I/pelec *10**3)

sigma_x2 = sigma_x**2
sigma_y2 = sigma_y**2

'''
poptx, pcovx = so.curve_fit(poly2, k, sigma_x2, sigma=err)
perrx = np.sqrt(np.diag(pcovx))

poptz, pcovz = so.curve_fit(poly2, k, sigma_y2, sigma=err)
perrz = np.sqrt(np.diag(pcovz))


fig, axs = plt.subplots(1, 2,figsize=(9, 4.5))

Domainx = np.arange(-50, 100, 0.01)
axs[0].set_title(r'X-direction')
axs[0].plot(Domainx,poly2(Domainx, *poptx),label= r'$\chi^2 =$ %.2f' %(chisqure(sigma_x2, poly2(k, *poptx))))
axs[0].errorbar(k, sigma_x2,yerr=err,xerr = errk, fmt=".")
axs[0].set_xlabel(r'k $\mathrm{[1/m^2]}$')
axs[0].set_ylabel(r'$\sigma_x^2$$\mathrm{[mm^2]}$')
axs[0].legend()

Domainz = np.arange(-100, 40, 0.01)
axs[1].set_title(r'Z-direction')
axs[1].plot(Domainz,poly2(Domainz, *poptz),label= r'$\chi^2 =$ %.2f' %(chisqure(sigma_y2, poly2(k, *poptz))))
axs[1].errorbar(k, sigma_y2,yerr=err,xerr = errk, fmt=".")
axs[1].set_xlabel(r'k $\mathrm{[1/m^2]}$')
axs[1].set_ylabel(r'$\sigma_z^2$$\mathrm{[mm^2]}$')
axs[1].legend()
'''


#plt.savefig('Quadrupole_3.png')
  
#print(r'Screen & Direction &a $\marhrm{mm^2m^{-4}}$&b $\marhrm{mm^2m^{-2}}$&c $\marhrm{mm^2}$ \\ ')
#print('S3& x & $%.3f \pm %.3f$ & $%.2f \pm %.2f$ & $%.1f \pm %.1f$ \\\ ' %(poptx[0],perrx[0],poptx[1],perrx[1],poptx[2],perrx[2]))
#print('S3& z & $%.4f \pm %.4f$ & $%.3f \pm %.3f$ & $%.1f \pm %.1f$ \\\ ' %(poptz[0],perrz[0],poptz[1],perrz[1],poptz[2],perrz[2]))


#print(poptx)
#print(poptz)
#print(perrx)
#print(perrz)

for i in range(len(k)):
    print('$',round(k[i],1),'$ & $',round(sigma_x[i],2),'\pm 0.1 $ & $',round(sigma_y[i],2),'\pm 0.1 $ \\\ ')


