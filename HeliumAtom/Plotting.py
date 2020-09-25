#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:46:13 2020

@author: lordemomo
"""

import numpy as np
from matplotlib import pyplot as plt


#%% Sweep
plot_all = 1
plot_energy_levels = 1

eigenvalues, u0s = np.genfromtxt('u0s.dat', unpack = 'True')

plt.figure(figsize=(16,9))
plt.title('Boundary Term vs Eigenvalue - '
          +str(np.size(eigenvalues))+' points', fontsize=20)
plt.plot(eigenvalues, u0s)
plt.xlabel('Eigenvalue', fontsize=18)
plt.ylabel('u(0)', fontsize=18)
plt.grid()

spec='global'

if not plot_all:
    #Plot energy levels
    if plot_energy_levels:
        plt.xlim(-0.6,0)
        plt.ylim(-0.25,0.25)
        spec = 'bound_levels'
        
        for i in range(1,5):
            plt.vlines(-0.5/i**2, -0.2, 0.2,
                       color='r', label='n='+str(i))
            plt.text(-0.5/i**2, 0.2, 'n='+str(i),
                     horizontalalignment='center', fontsize=14,
                     bbox={'facecolor': 'red', 'alpha': 1, 'pad': 5})
    else:
        plt.xlim(0,5)
        plt.ylim(-0.15,0.15)
        spec='unknown'

plt.tight_layout()
plt.savefig(dpi = 200,fname = 'BoundaryTerm_vs_Eigenvalue--'+spec+'.png')