import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import sys, os

# Change path according to your machine
sys.path.append('/home/lordemomo/Documentos/Density_Functional_Theory/Thijssen_Chapter5/')

# Use on terminal to create pydft module: $f2py -c -m pydft pydft.f95
import pydft
print(pydft.__doc__)

#%%
def func(x, a, b, c):
    return a/(x**c) + b

# Parameters
r_range = [0, 20]
Eigenvalue_range = [-1, 0.]

KS_int_max = 1

Eigenvalue_tol = 0.0001
u0_tol= 0.0001

Uniform_Numerov = [False, True]

h = 0.0001

write_data = False

rp = 50

j_max = np.arange(1000, 10000, 100, dtype=int)
delta = rp/j_max

u0 = np.zeros(np.size(j_max))

for i in range(0,int(np.size(u0))):
    print(i)
    print(delta[i]*j_max[i])
    u0[i] = np.abs(pydft.dft.hydrogenatom(r_range, Eigenvalue_range, KS_int_max, Eigenvalue_tol, u0_tol,
                  Uniform_Numerov, h, j_max[i], delta[i], write_data))

#%% Fitting
popt, pcov = curve_fit(func, j_max, u0)

#%%
plt.figure(figsize=(16,9))
plt.title('Termo de fronteira $u_0$ vs $N$', fontsize=20)
plt.plot(j_max, u0, label='Calculado', linewidth=4)
plt.plot(j_max, func(j_max, *popt), '-.', color='r',
          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt), linewidth = 4)
plt.xlabel('$N=j_{max}$', fontsize=18)
plt.ylabel('$|u_0|$', fontsize=18)
plt.xlim(j_max[0],j_max[-1])

# plt.yscale('log')
plt.xscale('log')
plt.legend(loc='best', fontsize=14)
plt.grid()
plt.savefig('/home/lordemomo/Documentos/Density_Functional_Theory/Thijssen_Chapter5/img/u0vsN.png', dpi=200, bbox_inches='tight')

