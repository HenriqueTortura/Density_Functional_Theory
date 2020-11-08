import numpy as np
from matplotlib import pyplot as plt
import sys, os

sys.path.append('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/')
import pydft
print(pydft.hydrogenatom.__doc__)

#%% 5.1e
if 1:
    r_range = [0, 50]
    Eigenvalue_range = [-5, 0]

    KS_int_max = 100

    Eigenvalue_tol = 0.00001
    u0_tol= 0.001

    Uniform_Numerov = [False, True]

    h = 0.0001

    j_max = 1000
    delta = 0.0001

    write_data = False
    eigenvalue = 0
    
    j_max = np.arange(100, 500, 10, dtype=int)
    erro = np.zeros(np.size(j_max))
    
    for i in range(0,int(np.size(erro))):
        print(i)
        erro[i] = pydft.hydrogenatom(r_range, Eigenvalue_range, KS_int_max, Eigenvalue_tol, u0_tol,
                      Uniform_Numerov, h, j_max[i], delta, write_data)
#%%
    an = erro/(1/j_max**4)
    plt.figure(figsize=(16,9))
    plt.title('Função de onda radial $u(r) = r\cdot\psi(r)$', fontsize=20)
    plt.plot(j_max, erro, label='Exata: $u(r)=2re^{-r}$', linewidth=4)
    plt.plot(j_max, an, label='a', linewidth=4)
    plt.xlabel('$N$', fontsize=18)
    plt.ylabel('$|\epsilon-\epsilon_{calc}|$', fontsize=18)
    plt.legend(loc='best', fontsize=14)
    plt.grid()

#%% Wave function and Potential
if 0:
    r, u = np.genfromtxt('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/Hydrogen_u.dat', unpack = 'True')
    r2, potentialU = np.genfromtxt('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/Hydrogen_Potential_U.dat', unpack = 'True')
    
    deltau = u - 2*r*np.exp(-r)
    deltaPotential = potentialU + (r+1)*np.exp(-2*r)-1
    
    #%%
    plt.figure(figsize=(16,9))
    plt.title('Função de onda radial $u(r) = r\cdot\psi(r)$', fontsize=20)
    plt.plot(r, 2*r*np.exp(-r), label='Exata: $u(r)=2re^{-r}$', linewidth=4)
    plt.plot(r, u, '-.', label='Calculada', linewidth=4, color='red')
    plt.xlabel('$r$', fontsize=18)
    plt.ylabel('$u(r)$', fontsize=18)
    plt.xlim(0,10)
    plt.xticks(np.arange(0, 11, step=1))
    plt.legend(loc='best', fontsize=14)
    plt.grid()
    plt.savefig('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/img/RadialWaveFunction.png', dpi=200, bbox_inches='tight')
    
    plt.figure(figsize=(16,9))
    plt.title('Potencial $U(r) = r\cdot V_H(r)$', fontsize=20)
    plt.plot(r, -(r+1)*np.exp(-2*r)+1, label='Exato: $U(r) = -(r+1)e^{-2r}+1$', linewidth=4)
    plt.plot(r, potentialU, '-.', label='Calculado', linewidth=4, color='red')
    plt.xlabel('$r$', fontsize=18)
    plt.ylabel('$U(r)$', fontsize=18)
    plt.xlim(0,10)
    plt.xticks(np.arange(0, 11, step=1))
    plt.legend(loc='best', fontsize=14)
    plt.grid()
    plt.savefig('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/img/RadialPotential.png', dpi=200, bbox_inches='tight')

#%% Probability
if 0:
    r, u = np.genfromtxt('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/Hydrogen_u.dat', unpack = 'True')
    prob = u**2/(2*np.pi)
    
    theta = np.linspace(0,2*np.pi,50)
    A, B = np.meshgrid(theta, r)
    c = np.zeros(A.shape)
    
    T, probabilidade_2D = np.meshgrid(theta, prob)
    
    figure, ax = plt.subplots(1, subplot_kw=dict(projection='polar'))
    figure.set_size_inches((16,16))
    plot = ax.contourf(theta, r, probabilidade_2D, 100, cmap=plt.cm.inferno)
    plt.colorbar(plot, ax=ax)
    
    plt.title('Probability 2D', fontsize=20)
    ax.set_ylim(0,3)
    #ax.grid(False)
    plt.savefig('img/hydrogen2D.png', dpi=200, bbox_inches='tight')


    plt.figure(figsize=(16,9))
    plt.title('Probability 1D', fontsize=20)
    plt.plot(r, prob*2*np.pi)
    plt.xlabel('r', fontsize=18)
    plt.ylabel('$u(r)^2$', fontsize=18)
    plt.xlim(0,10)
    plt.xticks(np.arange(0, 10, step=1))
    plt.grid()
    plt.savefig('img/hydrogen1D.png', dpi=200, bbox_inches='tight')

#%% Sweep Hydrogen
if 0:
    plot_all = 1
    plot_energy_levels = 1
    
    eigenvalues, u0s = np.genfromtxt('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/Hydrogen_u0s.dat', unpack = 'True')
    
    plt.figure(figsize=(16,9))
    plt.title('Termo de fronteira $u(0)$ vs Autovalor - '
              +str(np.size(eigenvalues))+' pontos', fontsize=20)
    plt.plot(eigenvalues, u0s)
    plt.xlabel('Autovalor', fontsize=18)
    plt.ylabel('$u(0)$', fontsize=18)
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
    plt.savefig('/home/lordemomo/Documentos/Density_Functional_Theory/HeliumAtom/img/BoundaryTerm_vs_Eigenvalue--'+spec+'.png',
                dpi=200, bbox_inches='tight')