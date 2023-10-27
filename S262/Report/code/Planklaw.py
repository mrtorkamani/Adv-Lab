#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 20:16:18 2023

@author: mr
"""

import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})


def planck_law_frequency(frequency, temperature):
    h = 6.626e-34  # Planck's constant (Joule second)
    kB = 1.38e-23  # Boltzmann constant (Joules per Kelvin)
    
    numerator = 2 * h * frequency**3 / (c**2)
    denominator = np.exp((h * frequency) / (kB * temperature)) - 1
    
    return numerator / denominator

frequency_range = np.linspace(1e12, 3e15, 1000)  # Frequency range from 1 THz to 3 PHz
temperature = 6000  # Temperature in Kelvin
c = 3.0e8  # Speed of light (meters per second)

wavelengths = c / frequency_range  # Convert frequency to wavelength


plt.figure(figsize=(10, 6))
plt.style.use('custom_style.mplstyle')
plt.plot(frequency_range * 1e-12, planck_law_frequency(frequency_range, 1000), label =r'1000 $\mathrm{K}$')
plt.plot(frequency_range * 1e-12, planck_law_frequency(frequency_range, 4000), label =r'4000 $\mathrm{K}$')
plt.plot(frequency_range * 1e-12, planck_law_frequency(frequency_range, 8000), label =r'8000 $\mathrm{K}$')
plt.plot(frequency_range * 1e-12, planck_law_frequency(frequency_range, 5777), label =r'5777 $\mathrm{K}(T_\odot)$')
plt.xlabel(r'\textbf{Frequency ($\mathrm{THz}$)}')
plt.ylabel(r'\textbf{Spectral energy density ($\mathrm{W/m^2/Hz/sr)}$}')
plt.title(r'\textbf{Planck\'s Law (Frequency Domain)}')
plt.legend()
plt.grid(True)
plt.savefig('Plank_law.png')