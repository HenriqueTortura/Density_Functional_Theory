import numpy as np
from matplotlib import pyplot as plt

# Change path according to your machine
path = '/home/lordemomo/Documentos/Density_Functional_Theory/Thijssen_Chapter5/'

plot_Hydrogen_Wave_and_Potential = 1
plot_Sweep_Hydrogen = 1

#%% Wave function and Potential
if plot_Hydrogen_Wave_and_Potential:
    r, u = np.genfromtxt(path+'Hydrogen_u.dat', unpack = 'True')
    r2, potentialU = np.genfromtxt(path+'Hydrogen_Potential_U.dat', unpack = 'True')
    
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
    plt.savefig(path+'img/RadialWaveFunction.png', dpi=200, bbox_inches='tight')
    
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
    plt.savefig(path+'img/RadialPotential.png', dpi=200, bbox_inches='tight')


#%% Sweep Hydrogen
if plot_Sweep_Hydrogen:
    plot_all = 1
    plot_energy_levels = 1
    
    eigenvalues, u0s = np.genfromtxt(path+'Hydrogen_u0s.dat', unpack = 'True')
    
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
        #Plot positive eigenvalues
        else:
            plt.xlim(0,5)
            plt.ylim(-0.15,0.15)
            spec='unknown'
    
    plt.tight_layout()
    plt.savefig(path+'img/BoundaryTerm_vs_Eigenvalue--'+spec+'.png',
                dpi=200, bbox_inches='tight')