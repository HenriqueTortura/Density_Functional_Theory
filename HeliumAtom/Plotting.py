#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:46:13 2020

@author: lordemomo
"""

import numpy as np
from matplotlib import pyplot as plt


#%%
if 1:
    r, u = np.genfromtxt('Hydrogen_u.dat', unpack = 'True', skip_footer=400000)
    prob = u**2/(2*np.pi)
    
    theta = np.linspace(0,2*np.pi,50)
    A, B = np.meshgrid(theta, r)
    c = np.zeros(A.shape)
    
    T, probabilidade_2D = np.meshgrid(theta, prob)

#%%
    figure, ax = plt.subplots(1, subplot_kw=dict(projection='polar'))
    figure.set_size_inches((16,16))
    plot = ax.contourf(theta, r, probabilidade_2D, 100, cmap=plt.cm.inferno)
    plt.colorbar(plot, ax=ax)
    
    plt.title('Probability 2D', fontsize=20)
    ax.set_ylim(0,3)
    #ax.grid(False)
    plt.savefig('hydrogen2D.png', dpi=200, bbox_inches='tight')


    plt.figure(figsize=(16,9))
    plt.title('Probability 1D', fontsize=20)
    plt.plot(r, prob*2*np.pi)
    plt.xlabel('r', fontsize=18)
    plt.ylabel('$u(r)^2$', fontsize=18)
    plt.xlim(0,10)
    plt.xticks(np.arange(0, 10, step=1))
    plt.grid()
    plt.savefig('hydrogen1D.png', dpi=200, bbox_inches='tight')

#%% Sweep
if 1:
    plot_all = 1
    plot_energy_levels = 1
    
    eigenvalues, u0s = np.genfromtxt('Hydrogen_u0s.dat', unpack = 'True')
    
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
    plt.savefig('BoundaryTerm_vs_Eigenvalue--'+spec+'.png',
                dpi=200, bbox_inches='tight')