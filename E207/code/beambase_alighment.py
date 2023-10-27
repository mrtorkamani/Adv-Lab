#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 12:59:26 2023

@author: mr
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import csv

plt.style.use('custom_style.mplstyle')


def chisqure(a,b):
    # no of hours a student studies
    # in a week vs expected no of hours
    observed_data = a
    expected_data = b
      
      
    # determining chi square goodness of fit using formula
    chi_square_test_statistic1 = 0
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            np.square((observed_data[i]-expected_data[i])/expected_data[i])
    return chi_square_test_statistic1

def linear(x,m,b):
    return m*x +b

data1 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Beam Alignment/Q1.csv')
data2 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Beam Alignment/Q2.csv')
data3 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Beam Alignment/Q3.csv')
data4 = pd.read_csv('/home/mr/Desktop/Lab/E207/E207/Part 1/Beam Alignment/Q4.csv')



fig, axs = plt.subplots(2, 2,figsize=(9, 9))

currentx = data1.Ix
dposx = data1.x
currenty = data1.Iy
dposy = data1.y
D = np.arange(-4,5,1)
error = [.01,.01,.01,.01,.01]

popt01, pcov01 = curve_fit(linear, -dposx, currentx,sigma= error)
perr01 = np.sqrt(np.diag(pcov01))

yfit1 = linear(D, *popt01)

popt02, pcov02 = curve_fit(linear, dposy, currenty,sigma= error)
perr02 = np.sqrt(np.diag(pcov02))

yfit2 = linear(D, *popt02)


chi01=chisqure(currentx, linear(-dposx, *popt01))
chi02=chisqure(currenty, linear(dposy, *popt02))


axs[0,0].set_title(r'C0')
axs[0,0].errorbar(-dposx ,currentx, xerr=error,fmt=".")
axs[0,0].plot(D,yfit1,label='X direction')
axs[0,0].errorbar(dposy ,currenty, xerr=error,fmt=".")
axs[0,0].plot(D,yfit2,label='Z direction')
axs[0,0].legend()

currentx = data2.Ix
dposx = data2.x
currenty = data2.Iy
dposy = data2.y
D = np.arange(-7.5,5.5,1)
error = [.2,.2,.2,.2,.2]

popt11, pcov11 = curve_fit(linear, dposx, currentx,  sigma= error)
perr11 = np.sqrt(np.diag(pcov11))
yfit1 = linear(D, *popt11)

popt12, pcov12 = curve_fit(linear, dposy, currenty,  sigma= error)
perr12 = np.sqrt(np.diag(pcov12))
yfit2 = linear(D, *popt12)

chi11=chisqure(currentx, linear(dposx, *popt11))
chi12=chisqure(currenty, linear(dposy, *popt12))

axs[0,1].set_title(r'C1')
axs[0,1].errorbar(dposx ,currentx, xerr=error,fmt=".")
axs[0,1].plot(D,yfit1,label='X direction')
axs[0,1].errorbar(dposy ,currenty, xerr=error,fmt=".")
axs[0,1].plot(D,yfit2,label='Z direction')
axs[0,1].legend()

currentx = data3.Ix
dposx = data3.x
currenty = data3.Iy
dposy = data3.y
D = np.arange(-4,5,1)
error = [.2,.2,.2,.2,.2]

popt21, pcov21 = curve_fit(linear, dposx, currentx,  sigma= error)
perr21 = np.sqrt(np.diag(pcov21))
yfit1 = linear(D, *popt21)

popt22, pcov22 = curve_fit(linear, -dposy, currenty,  sigma= error)
perr22 = np.sqrt(np.diag(pcov22))
yfit2 = linear(D, *popt22)

chi21=chisqure(currentx, linear(dposx, *popt21))
chi22=chisqure(currenty, linear(-dposy, *popt22))

axs[1,0].set_title(r'C2')
axs[1,0].errorbar(dposx ,currentx, xerr=error,fmt=".")
axs[1,0].plot(D,yfit1,label='X direction')
axs[1,0].errorbar(-dposy ,currenty, xerr=error,fmt=".")
axs[1,0].plot(D,yfit2,label='Z direction')
axs[1,0].legend()

currentx = data4.Ix
dposx = data4.x
currenty = data4.Iy
dposy = data4.y
D = np.arange(-2.5,1.5,1)
error = [.2,.2,.2,.2,.2]

popt31, pcov31 = curve_fit(linear, dposx, currentx,  sigma= error)
perr31 = np.sqrt(np.diag(pcov31))
yfit1 = linear(D, *popt31)

popt32, pcov32 = curve_fit(linear, -dposy, currenty,  sigma= error)
perr32 = np.sqrt(np.diag(pcov32))
yfit2 = linear(D, *popt32)

chi31=chisqure(currentx, linear(dposx, *popt31))
chi32=chisqure(currenty, linear(-dposy, *popt32))

axs[1,1].set_title(r'C3')
axs[1,1].errorbar(dposx ,currentx, xerr=error,fmt=".")
axs[1,1].plot(D,yfit1,label='X direction')
axs[1,1].errorbar(-dposy ,currenty, xerr=error,fmt=".")
axs[1,1].plot(D,yfit2,label='Z direction')
axs[1,1].legend()

for i in fig.get_axes():
    i.set_xlabel(r'Displacement $\mathrm{[mm]}$')
    i.set_ylabel(r'Current $\mathrm{[mA]}$')

fig.tight_layout()

plt.savefig('beambasealinghment.png')

print('Corrector & $I_x \mathrm{[mA]}$& $\chi^2_x$ &$I_y \mathrm{[mA]}$& $\chi^2_y$ \\\ ')
print('C0 & $%f \pm %f$ & $%f$ & $%f \pm %f$ & $%f$ \\\ ' %(popt01[1],perr01[1],chi01,popt02[1],perr02[1],chi02))
print('C1 & $%f \pm %f$ & $%f$ & $%f \pm %f$ & $%f$ \\\ ' %(popt11[1],perr11[1],chi11,popt12[1],perr12[1],chi12))
print('C2 & $%f \pm %f$ & $%f$ & $%f \pm %f$ & $%f$ \\\ ' %(popt21[1],perr21[1],chi21,popt22[1],perr22[1],chi22))
print('C3 & $%f \pm %f$ & $%f$ & $%f \pm %f$ & $%f$ \\\ ' %(popt31[1],perr31[1],chi31,popt32[1],perr32[1],chi32))

#for i in range(len(data1.Ix)):
#    print('$',round(data1.Ix[i],3),'$ & $', round(data1.x[i],2), '\pm 0.2 $ & $',round(data1.Iy[i],3),'$ & $', round(data1.y[i],2), '\pm 0.2 $ \\\ ')

