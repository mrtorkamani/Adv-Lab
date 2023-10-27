#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 18:36:07 2023

@author: mr
"""
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('custom_style.mplstyle')


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

data = pd.read_csv('/home/mr/Desktop/Lab/S263/codes/sex/script/result/Merge2.csv')

print(data['NUMBER'])

B_V =  (data['MAG_AUTO_B']+24.4999) - (data['MAG_AUTO_V']+24.6943)
V = data['MAG_AUTO_V'] + 24.6943

plt.figure(figsize=(9,6))
plt.scatter(B_V, V)
plt.ylabel(r'\textbf{app mag(V)}')
plt.xlabel(r'\textbf{B-V}')
plt.gca().invert_yaxis()
plt.xlim(-3,3)
plt.savefig('(5,5)CMD.png')
