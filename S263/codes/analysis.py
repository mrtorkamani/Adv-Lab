#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:20:16 2023

@author: mr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('custom_style.mplstyle')


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

def residualcalculator():
    res = [] #list for storing residual values
    shift= [] #list for storing shift values
    #loop over different shift a in lenghts of 0.01
    for a in np.arange(7.5,9.5,0.01):
        r = 0
        #loop over iso data  
        for i in range(len(data1[8][:253])): 
            VV = V - a
            #loop over data point and finding the right interval for calculating
            #the residual
            for j in range(len(B_V)):
                if data1[8][i]>B_V[j] and data1[8][i+1]<B_V[j]:
                    r += np.abs(VV[j]-.5*(data1[6][i]+data1[6][i+1]))
                else:
                    continue
        res.append(r)
        shift.append(a)
        print(r,a)
    print(shift[res.index(min(res))])#print the shift regarding the minimum residual 
    print(min(res)) #printing the residual
    
data = pd.read_csv('/home/mr/Desktop/Lab/S263/codes/sex/script/result/Merge.csv')

B_V =  (data['MAG_AUTO_B']+24.4999) - (data['MAG_AUTO_V']+24.6943)
V = data['MAG_AUTO_V'] + 24.6943 

# Open the file for reading
with open('/home/mr/Desktop/Lab/S263/Photometry/GC_and_isochrones/isochrones/020/iso_c020_0850.UBVRIJHKLM', 'r') as file:
    # Initialize an empty list for each column
    data1 = [[] for _ in range(21)]  # Assuming there are 21 columns
    
    # Read lines from the file
    lines = file.readlines()
    
    # Iterate through lines
    for line in lines:
        # Skip comment lines starting with '#'
        if not line.startswith('#'):
            # Split the line by whitespace and convert values to float
            values = [float(val) for val in line.split()]
            
            # Append each value to the respective column list
            for i in range(len(values)):
                data1[i].append(values[i])

residualcalculator()
'''
V = V - 8.38

plt.plot(data1[8][:235],data1[6][:235],color = 'black')
plt.gca().invert_yaxis()
plt.scatter(B_V, V)
plt.ylabel(r'\textbf{abs mag(V)}')
plt.xlabel(r'\textbf{B-V}')
plt.title(r'Z= 0.02 and Age $10^{8.50}$ yrs')
plt.xlim(-.5,1.2)
#plt.ylim(8,-2)

plt.savefig('iso_z020_850.png')
'''

