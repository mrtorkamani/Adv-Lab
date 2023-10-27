#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 21:00:40 2023

@author: mr
"""
import numpy as np
import pandas as pd
import sympy as sp

import matplotlib
import matplotlib.pyplot as plt

from IPython.display import Markdown as md

import scipy.optimize as so

import os


data =pd.read_csv('/home/mr/Desktop/Lab/E207/E207/part2/multiple screen/a.csv')

sigmax2 = [data['sigmax'][i]**2 for i in range(len(data['sigmax']))]
sigmaz2 = [data['sigmay'][i]**2 for i in range(len(data['sigmay']))]

sigmax2 = sigmax2[:]
sigmaz2 = sigmaz2[:]
def linear(x, m, b):
    return m*x+b

def poly2(x, a, b, c):
    return a*x**2 + b*x + c

def xsq(x, s, a, b):
    return a*(x-s)**2 + b

def quadrupole(k, L, dim="z"):
    if dim == "x":
        k = -k
    out = np.array( [ [1, L], [k * L, 1] ] )
    return out

def drift(L):
    return np.array( [ [1, L], [0, 1] ] )

def kickangle(mrad):
    return np.array([[0], [mrad]]).T

def solveDet(dictionary):
    a = float(dictionary[ea])
    b = float(dictionary[eb])
    g = float(dictionary[eg])
    emitt_sq = b * g - a**2
    return np.sqrt(emitt_sq)

##momentum calculation
epot = 511
ekin = 25
gamma = (epot + ekin)/epot
c = 299792458
velec = np.sqrt(1-(1/gamma)**2) *c 
pelec = 9.10938e-31 *gamma*velec  / 5.34428599e-28
I = np.array([0.024,-0.054,0.03,0.043])
k = 128027 * I/pelec *10**3
errk = np.abs(150* I/pelec *10**3)


ds0 = .0785/2 +.185 +0.02 # to S0
dq1 = ds0 + .04+.06 # to Q1

dqs = .03 +.0785+.316+0.02 # from Qx to Sx
dsq = .02+.06 # from Sx to Qx+1

Ms0 = drift(ds0)
Ms0q1 = drift(dq1) # drift til q1
    # Mq1
Mq1s1 = drift(dqs)   # up to here for s1
Ms1q2 = drift(dsq)
    # Mq2
Mq2s2 = drift(dqs)   # up to here for s2
Ms2q3 = drift(dsq)
    # Mq3
Mq3s3 = drift(dqs)   # up to here for s3
Ms3q4 = drift(dsq)
    # Mq4
Mq4s4 = drift(dqs)   # up to here for s4

    # dict for final matrices per dimension x/z
lgs = {"x" : [], "z" : []}
for d in "x", "z":
    # quadrupole matrices with used k's
    Mq = [quadrupole(k, .074, d) for k in k]
    # whole LAB in order of appearance
    M = [Ms0, Ms0q1, Mq[0], Mq1s1, Ms1q2, Mq[1], Mq2s2, Ms2q3, Mq[2], Mq3s3, Ms3q4, Mq[3], Mq4s4]
    # for all screens calculate the final matrix
    screen_posx = [0, 3+1, 6+1, 9+1, 12+1] #
    screen_posy = [0, 3+1, 6+1, 9+1, 12+1] #, 
    if d == "x":
        screen_pos = screen_posx[1:]
    else:
        screen_pos = screen_posy[1:]
    for N in screen_pos:
        # matrix up to first screen
        tmp = M[0]
        # multiply up to screen
        for i in range(1, N):
            tmp = np.matmul(M[i], tmp)
        #print("Dimension", d, f"screen {(N-1)/3}", tmp, "\n")
        # calculate vector for lgs
        m_element_vec = np.array([  tmp[0, 0]**2, -2*tmp[0,0]*tmp[0,1], tmp[0,1]  ])
        # store result in dict
        lgs[d].append(m_element_vec)
        
#print(lgs['x'])

ea, eb, eg = sp.symbols('\epsilon_a, \epsilon_b, \epsilon_g')
x = np.array([[eb], [ea], [eg]]).T
#print(data['sigmax'])

sigma = {"x" : np.array(sigmax2), "z" : np.array(sigmaz2)}


final_lgs = {}


for d in "x", "z":
    LGS = np.array(lgs[d])

    lhs = LGS.T * sigma[d].reshape(sigma[d].size, 1).T
    rhs = np.matmul(LGS.T, LGS) * x

    final_lgs[d] = [sp.Eq(eq, 0) for eq in rhs.sum(axis=1) - lhs.sum(axis=1)]
    solution = sp.solve(final_lgs[d], dict=False)
        
    print(d, sp.solve(final_lgs[d], dict=False), "\n")

#     print(0*"### with BBA Quadrupole strength" +"\n" + \
#        f"$$\epsilon_x = {solveDet(sp.solve(final_lgs['x'], dict=False))}$$ $$ \epsilon_z = {solveDet(sp.solve(final_lgs['z'], dict=False))}$$")
print(solveDet(sp.solve(final_lgs['x'])), solveDet(sp.solve(final_lgs['z'])))

print(k)
