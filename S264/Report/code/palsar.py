#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 22:53:06 2023

@author: mr
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy import linalg as la
import astropy.units as u

#Fitfunctions:

#Chi squared goodness of the fit:
def Goodnessofthefit(y, ndof):
    return(chisquare(y)[0]/ndof)

#Gaussian:
def Gaussfit(x, A, sigma, m, d):
    return (A/(np.sqrt(2*np.pi)*sigma))*np.exp(-(x-m)**2/(2*sigma**2)) + d

#linear fit:
def Linfit(d,m,n):
    return(m*d+n)
    
    
#Loading data Pulsar:
file = open('/home/mr/Desktop/Lab/S264/S624_Data_A7/0329_1246_8.ascii', 'r')
channel = []
sub1 = []
sub2 = []
sub3 = []
sub4 = []
sub5 = []
sub6 = []
sub7 = []
sub8 = []

for line in file:
    data = line.strip('\n').split(' ')
    if len(data)==9:
        channel.append(float(line.strip('\n').split(' ')[0]))
        sub1.append(float(line.strip('\n').split(' ')[1]))
        sub2.append(float(line.strip('\n').split(' ')[2]))
        sub3.append(float(line.strip('\n').split(' ')[3]))
        sub4.append(float(line.strip('\n').split(' ')[4]))
        sub5.append(float(line.strip('\n').split(' ')[5]))
        sub6.append(float(line.strip('\n').split(' ')[6]))
        sub7.append(float(line.strip('\n').split(' ')[7]))
        sub8.append(float(line.strip('\n').split(' ')[8]))
        
t_channel = 714.577/256 #ms   total measuring time divided by number of channels
bandwidth = 98.4/8.0 #MHz   spectrum subdivided into 8 frequnecy bands
central_frequency = []
for i in range(1,9):
    central_frequency.append(1430.391 - bandwidth * (i - 0.5)) #MHz
t = t_channel * np.asarray(channel) #ms

print('t_channel:',t_channel)
print('band width sub-bands:',bandwidth)


plt.style.use('custom_style.mplstyle')

#Plotting:
fig, axes = plt.subplots(figsize=(15, 6))
#ax = fig.add_axes([0.55, 0.6, 0.15, 0.2])
#ax.errorbar(t[20:45], sub1[20:45], label=f'''Band 1 $v_1$={central_frequency[0]:0.3f}''', color='b')
#ax.errorbar(t[20:45], sub2[20:45], label=f'''Band 2 $v_1$={central_frequency[0]:0.3f}''', color='g')
#ax.errorbar(t[20:45], sub3[20:45], label=f'''Band 3 $v_1$={central_frequency[0]:0.3f}''', color='r')
#ax.errorbar(t[20:45], sub4[20:45], label=f'''Band 4 $v_1$={central_frequency[0]:0.3f}''', color='c')
#ax.errorbar(t[20:45], sub5[20:45], label=f'''Band 5 $v_1$={central_frequency[0]:0.3f}''', color='m')
#ax.errorbar(t[20:45], sub6[20:45], label=f'''Band 6 $v_1$={central_frequency[0]:0.3f}''', color='y')
#ax.errorbar(t[20:45], sub7[20:45], label=f'''Band 7 $v_1$={central_frequency[0]:0.3f}''', color='k')
#ax.errorbar(t[20:45], sub8[20:45], label=f'''Band 8 $v_1$={central_frequency[0]:0.3f}''', color='#fd978b')
#ax.grid(which='both')
#ax.set_title('Pulsar signal')
axes.errorbar(t, sub1, label=f'''Band 1 $v_1$={central_frequency[0]:0.3f}''', color='b')
axes.errorbar(t, sub2, label=f'''Band 2 $v_2$={central_frequency[1]:0.3f}''', color='g')
axes.errorbar(t, sub3, label=f'''Band 3 $v_3$={central_frequency[2]:0.3f}''', color='r')
axes.errorbar(t, sub4, label=f'''Band 4 $v_4$={central_frequency[3]:0.3f}''', color='c')
axes.errorbar(t, sub5, label=f'''Band 5 $v_5$={central_frequency[4]:0.3f}''', color='m')
axes.errorbar(t, sub6, label=f'''Band 6 $v_6$={central_frequency[5]:0.3f}''', color='y')
axes.errorbar(t, sub7, label=f'''Band 7 $v_7$={central_frequency[6]:0.3f}''', color='k')
axes.errorbar(t, sub8, label=f'''Band 8 $v_8$={central_frequency[7]:0.3f}''', color='#fd978b')
axes.set_xlabel('Time [ms]')
axes.set_ylabel('Intensity [a.U.]')
axes.set_xlim(0,500)
axes.legend()
axes.grid(which='both')
fig.savefig("Data_frequency_bands.png")


TT = np.linspace(30,180,1000)
#gaussian fits:
popt1, pcov1 = curve_fit(Gaussfit, t[10:65], sub1[10:65],[10**7, 10, 90, 0])
yfit1 = Gaussfit(TT, *popt1)
perr1 = np.sqrt(np.diag(pcov1))

popt2, pcov2 = curve_fit(Gaussfit, t[10:65], sub2[10:65],[6*10**7, 10, 90, 0])
yfit2 = Gaussfit(TT, *popt2)
perr2 = np.sqrt(np.diag(pcov2))

popt3, pcov3 = curve_fit(Gaussfit, t[10:65], sub3[10:65],[6*10**7, 10, 90, 0])
yfit3 = Gaussfit(TT, *popt3)
perr3 = np.sqrt(np.diag(pcov3))

popt4, pcov4 = curve_fit(Gaussfit, t[10:65], sub4[10:65],[6*10**7, 10, 90, 0])
yfit4 = Gaussfit(TT, *popt4)
perr4 = np.sqrt(np.diag(pcov4))

popt5, pcov5 = curve_fit(Gaussfit, t[10:65], sub5[10:65],[6*10**7, 10, 90, 0])
yfit5 = Gaussfit(TT, *popt5)
perr5 = np.sqrt(np.diag(pcov5))

popt6, pcov6 = curve_fit(Gaussfit, t[10:65], sub6[10:65],[6*10**7, 10, 90, 0])
yfit6 = Gaussfit(TT, *popt6)
perr6 = np.sqrt(np.diag(pcov6))

popt7, pcov7 = curve_fit(Gaussfit, t[10:65], sub7[10:65],[6*10**7, 10, 90, 0])
yfit7 = Gaussfit(TT, *popt7)
perr7 = np.sqrt(np.diag(pcov7))

popt8, pcov8 = curve_fit(Gaussfit, t[10:65], sub8[10:65],[6*10**7, 10, 90, 0])
yfit8 = Gaussfit(TT, *popt8)
perr8 = np.sqrt(np.diag(pcov8))


'''
fig, axes = plt.subplots(3,3,figsize=(15, 15))

axes[0][0].errorbar(t[10:65], sub1[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[0][0].plot(TT, yfit1, color='red', label=f''Gaussian fit'')
axes[0][0].set_xlabel('Time [ms]')
axes[0][0].set_ylabel('Intensity [a.U.]')
axes[0][0].legend()
axes[0][0].grid(which='both')
axes[0][0].set_title(f''Frequency band 1'')

axes[0][1].errorbar(t[10:65], sub2[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[0][1].plot(TT, yfit2, color='red', label=f''Gaussian fit'')
axes[0][1].set_xlabel('Time [ms]')
axes[0][1].set_ylabel('Intensity [a.U.]')
axes[0][1].legend()
axes[0][1].grid(which='both')
axes[0][1].set_title(f''Frequency band 2'')

axes[0][2].errorbar(t[10:65], sub3[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[0][2].plot(TT, yfit3, color='red', label=f''Gaussian fit'')
axes[0][2].set_xlabel('Time [ms]')
axes[0][2].set_ylabel('Intensity [a.U.]')
axes[0][2].legend()
axes[0][2].grid(which='both')
axes[0][2].set_title(f''Frequency band 3'')

axes[1][0].errorbar(t[10:65], sub4[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[1][0].plot(TT, yfit4, color='red', label=f''Gaussian fit'')
axes[1][0].set_xlabel('Time [ms]')
axes[1][0].set_ylabel('Intensity [a.U.]')
axes[1][0].legend()
axes[1][0].grid(which='both')
axes[1][0].set_title(f''Frequency band 4'')

axes[1][1].errorbar(t[10:65], sub8[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[1][1].plot(TT, yfit8, color='red', label=f''Gaussian fit'')
axes[1][1].set_xlabel('Time [ms]')
axes[1][1].set_ylabel('Intensity [a.U.]')
axes[1][1].legend()
axes[1][1].grid(which='both')
axes[1][1].set_title(f''Frequency band 8'')

fig.delaxes(axes[1, 2])

axes[2][0].errorbar(t[10:65], sub5[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[2][0].plot(TT, yfit5, color='red', label=f''Gaussian fit'')
axes[2][0].set_xlabel('Time [ms]')
axes[2][0].set_ylabel('Intensity [a.U.]')
axes[2][0].legend()
axes[2][0].grid(which='both')
axes[2][0].set_title(f''Frequency band 5'')

axes[2][1].errorbar(t[10:65], sub6[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[2][1].plot(TT, yfit6, color='red', label=f''Gaussian fit'')
axes[2][1].set_xlabel('Time [ms]')
axes[2][1].set_ylabel('Intensity [a.U.]')
axes[2][1].legend()
axes[2][1].grid(which='both')
axes[2][1].set_title(f''Frequency band 6'')

axes[2][2].errorbar(t[10:65], sub7[10:65], fmt = '.', zorder=1, markersize = 8, color='black')
axes[2][2].plot(TT, yfit7, color='red', label=f''Gaussian fit'')
axes[2][2].set_xlabel('Time [ms]')
axes[2][2].set_ylabel('Intensity [a.U.]')
axes[2][2].legend()
axes[2][2].grid(which='both')
axes[2][2].set_title(f''Frequency band 7'')


fig.tight_layout()
fig.savefig("Gaussian_good_fits_pulsar.png")
'''

#fig, axes1 = plt.subplots(1,3,figsize=(15, 6))


#fig.tight_layout()
#fig.savefig("Gaussian_bad_fits_pulsar.png")


peak_position = [popt1[2], popt2[2], popt3[2], popt4[2], 90, 90, 90, popt8[2]]
dpeak_position = [perr1[2], perr2[2], perr3[2], perr4[2], perr5[2], perr6[2], perr7[2], perr8[2]]

#popt, pcov = curve_fit(Linfit, np.asarray(central_frequency), peak_position, sigma=dpeak_position)
#yfit = Linfit(np.asarray(central_frequency), *popt)
#perr = np.sqrt(np.diag(pcov))

#fig, axes = plt.subplots(figsize=(15, 6))
#axes.errorbar(central_frequency, peak_position, yerr=dpeak_position, label='Data', fmt = '.', zorder=1, markersize = 8, color='black')
#axes.plot(central_frequency, yfit, color='red', label=f'''Goodness of the fit: {Goodnessofthefit(yfit, 6):0.2f}''')
#axes.set_ylabel('Peak position [ms]')
#axes.set_xlabel(f'''$v_i$ [MHz]''')
#axes.legend()
#axes.grid(which='both')
#fig.savefig("Pulsar_linfit.pdf")

#print('Fit Peak\,Position = (\\num{',popt[0],'+-',perr[0],'}) \cdot v_i\,\mathrm{[MHz]} + (\\num{',popt[1],'+-',perr[1],'})')

t1 = (popt1[2])*10**(-3)
dt1 = (np.sqrt((perr1[2])**2))*10**(-3)
t8 = (popt8[2])*10**(-3)
dt8 = (np.sqrt((perr8[2])**2))*10**(-3)
v1 = central_frequency[0] 
v8 = central_frequency[7]
dv1 = 0
dv8 = 0

#Dispersion measure DM:
DM = 1/(4.15*10**3) * (t8-t1)/  (v8**(-2)-v1**(-2))
dDM = DM * (dt8+dt1)/(t8-t1)

#Distance to Pulsar:
ne = 0.03 #1/cm**3 electron density
D = DM/ne
dD = dDM/ne

print(peak_position)
print(fr't1 = {popt1[2]:.2f} +-{perr1[2]:.2f}')
print(fr't2 = {popt2[2]:.2f} +-{perr2[2]:.2f}')
print(fr't3 = {popt3[2]:.2f} +-{perr3[2]:.2f}')
print(fr't4 = {popt4[2]:.2f} +-{perr4[2]:.2f}')
print(fr't5 = {popt5[2]:.2f} +-{perr5[2]:.2f}')
print(fr't6 = {popt6[2]:.2f} +-{perr6[2]:.2f}')
print(fr't7 = {popt7[2]:.2f} +-{perr7[2]:.2f}')
print(fr't8 = {popt8[2]:.2f} +-{perr8[2]:.2f}')

print('t1 =',t1*10**3,'+-',dt1*10**3,'ms')
print('t8 =',t8*10**3,'+-',dt8*10**3,'ms')
print('DM =',DM,'+-',dDM,'pc cm^-3')
print('D = ',D,'+-',dD,'pc')






