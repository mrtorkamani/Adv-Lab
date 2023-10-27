#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:11:49 2023

@author: mr
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import csv
from scipy.optimize import curve_fit
plt.style.use('custom_style.mplstyle')

def linear_function(x, m, b):
    return m * x + b

address= '/home/mr/Desktop/Lab/A246/A264/Exercise 5/Diode_Laser_Charectristics.csv'

data = pd.read_csv(address)

#data.plot(kind = 'scatter', x = 'Current (mA)', y = 'Power (µW)',color = 'black')
#plt.savefig('laserchar.png')

current = (data['Current (mA)'].tolist())
power = (data['Power (µW)'].tolist())

params, covariance = curve_fit(linear_function,current[14:],power[14:])

m_err = np.sqrt(covariance[0, 0])
b_err = np.sqrt(covariance[1, 1])

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(current,power, color = 'black')
x = np.linspace(54,100,100)
ax.plot(x,linear_function(x,params[0],params[1]))
ax.set_xlabel(r'Injected current $[\mathrm{mA}]$')
ax.set_ylabel(r'Output power $[\mathrm{\mu W}]$')
ax.grid()
fig.savefig('Laserchar.png')

print('fiting parameters')
print(fr'slope = {params[0]:.2f} ± {m_err:.2f}')
print(fr'intercept = {params[1]:.1f} ± {b_err:.1f}')
I_thr = (power[14] - params[1])/params[0]
yerr = 1
I_thr_err = I_thr*((b_err+yerr) /(power[14] - params[1]) + m_err/params[0])
print(fr'I_thr = {I_thr:.0f}±{I_thr_err:.0f}')


#########quantum efficency######
h= sp.constants.h
e = sp.constants.e
c = sp.constants.c
lambdaa = 749e-9
E = h*c/lambdaa
qe = []
for i in range(14,len(current)):
    qe.append((power[i]*e) / (current[i]*E) *10**-3)    

fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.scatter(current[14:],qe, color = 'black')
ax1.set_xlabel(r'Injected current $[\mathrm{mA}]$')
ax1.set_ylabel(r'quantum efficency')
ax1.grid()
fig1.savefig('Quantum_efficency.png')
