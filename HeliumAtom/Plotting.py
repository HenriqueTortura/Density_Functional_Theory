#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:46:13 2020

@author: lordemomo
"""

import numpy as np
from matplotlib import pyplot as plt

eigenvalues, u0s = np.genfromtxt('u0s.dat', unpack = 'True')

#%%
plt.figure(figsize=(16,9))
plt.title('Boundary Term vs Eigenvalue')
plt.plot(eigenvalues, u0s)
plt.xlabel('Eigenvalue')
plt.ylabel('u(0)')
plt.grid()
plt.tight_layout()
plt.savefig(dpi = 200,fname = 'BoundaryTerm_vs_Eigenvalue--'
            +str(np.size(eigenvalues))+'points--from'
            +str(np.min(eigenvalues))+
            'to'+str(np.max(eigenvalues))
            +'.png')