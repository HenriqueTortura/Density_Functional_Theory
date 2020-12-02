import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import sys

# Change path according to your machine
path = '/home/lordemomo/Documentos/Density_Functional_Theory/Thijssen_Chapter5/'
sys.path.append(path)

# Use on terminal to create pydft module: $f2py -c -m pydft dft.f95
import pydft
print(pydft.__doc__)
# print(pydft.dft.__doc__)

def yes_or_no(question): #Straight from https://stackoverflow.com/questions/47735267/while-loop-with-yes-no-input-python
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return 1
    elif reply[0] == 'n':
        return 0
    else:
        return yes_or_no("Please Enter")

#%% Section 5.5
run_Section55 = yes_or_no('Run section 5.5 (be patient with the plotting) ?')
if run_Section55:
    print('#####################')
    print('#### Section 5.5')
    print('#### A density functional program for the helium atom')
    
    #%% 5.5.1 Solving the radial equation
    print('\n')
    print('####### 5.5.1 Solving the radial equation')
    
    # Parameters
    r_range = [0, 50]
    Eigenvalue_range = [-5, 0.]
    
    u0_tol= 0.001
    Eigenvalue_tol = 0.00001
    h = 0.00001
    
    KS_int_max = 100
    
    Uniform_Numerov = [True, True] # Starting with uniform grid
    
    verbose = False
    write_data = True # F2PY doesn't allow to allocate arrays inside a
                      # subroutine, so I chose to write u(r) and U(r)
                      # in a file
    
    # Calculating
    hydrogen_eigenvalue, u0, HE, EE, CE = np.abs(pydft.dft.hydrogenatom(r_range,
                                          Eigenvalue_range,
                                          KS_int_max, Eigenvalue_tol,
                                          u0_tol, Uniform_Numerov, h,
                                          200000, 0.0001,
                                          verbose, write_data, path+'data/'))
    
    print('Eigenvalue: '+str(hydrogen_eigenvalue))
    print('Expected: 0.5')
    print('Absolute error: '+str(abs(0.5-hydrogen_eigenvalue)))
    
    # Plotting
    print('\n')
    print('Checking u(r) for exact hydrogen solution:')
    r, u, potentialU = np.genfromtxt(path+'data/Hydrogen_data.dat', unpack = 'True')
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
    plt.show()
    plt.savefig(path+'img/RadialWaveFunction.png', dpi=200, bbox_inches='tight')
    
    #%% 5.5.2 Including the Hartree potential
    print('\n')
    print('####### 5.5.2 Including the Hartree potential')
    
    print('Checking U(r) for exact hydrogen solution:')
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
    plt.show()
    plt.savefig(path+'img/RadialPotential.png', dpi=200, bbox_inches='tight')
    #%%
    print('\n')
    print('Checking helium energy with self-interaction correction (without exchange):')
    
    # Parameters
    SelfCons_int_max = 100
    SelfCons_tol = 0.0001
    
    # Helium
    Z = 2
    N_electrons = 2
    
    Exchange_Method = 0 # No exchange
    SelfInteractionCorrection = True # Self-interaction correction
                                     # proposed by Thijssen
    h = 0.0001 # Guarantees convergence
    #h = 0.00001 # For a better result in energy, but
                 #without convergence in 100 self-consistency interactions
    
    Energy, eigenvalue, HC, EC, CC = pydft.dft.heliumatom(r_range,
                                     Eigenvalue_range, SelfCons_int_max,
                                     KS_int_max, Eigenvalue_tol, u0_tol,
                                     SelfCons_tol, Uniform_Numerov, h,
                                     200000, 0.0001, verbose,
                                     SelfInteractionCorrection,
                                     Exchange_Method, 0, Z, N_electrons)
    
    print('## Helium eigenvalue')
    print('Calculated: '+str(eigenvalue))
    print('Expected by Thijssen: -0.923')
    
    print('## Hartree Correction')
    print('Calculated: '+str(HC))
    print('Expected by Thijssen: -1.0155')
    
    print('## Helium energy')
    print('Calculated: '+str(Energy))
    print('Expected by Thijssen: -2.861')
    print('Experimental value: -2.903')
    
    #%% 5.5.3 The local density exchange potential
    print('\n')
    print('####### 5.5.3 The local density exchange potential')
    
    # Parameters
    Exchange_Method = 1 # LDA Exchange
    SelfInteractionCorrection = False # No self-interaction correction
    
    Energy, eigenvalue, HC, EC, CC = pydft.dft.heliumatom(r_range,
                         Eigenvalue_range, SelfCons_int_max, KS_int_max,
                         Eigenvalue_tol, u0_tol, SelfCons_tol,
                         Uniform_Numerov, h, 200000, 0.0001, verbose,
                         SelfInteractionCorrection,
                         Exchange_Method, 0, Z, N_electrons)
    
    print('## Helium eigenvalue')
    print('Calculated: '+str(eigenvalue))
    print('Expected by Thijssen: -0.52')
    
    print('## Helium energy')
    print('Calculated: '+str(Energy))
    print('Expected by Thijssen: -2.72')

#%% Extra: Extra: Sweep hydrogen eigenvalues
plot_Sweep_Hydrogen = yes_or_no('Extra: Sweep hydrogen eigenvalues?')
if plot_Sweep_Hydrogen:
    plot_all = 0
    plot_energy_levels = 1
    
    pydft.dft.kohnshamsweep([0.0,50.0], [-1.0,5.0],
                            0.0001, 1001, path+'data/')
    eigenvalues, u0s = np.genfromtxt(path+'data/Hydrogen_u0s.dat', unpack = 'True')
    
    plt.figure(figsize=(16,9))
    plt.title('Termo de fronteira $u(0)$ vs Autovalor - '
              +str(np.size(eigenvalues))+' pontos', fontsize=20)
    plt.plot(eigenvalues, u0s)
    plt.xlabel('Autovalor', fontsize=18)
    plt.ylabel('$u(0)$', fontsize=18)
    plt.grid()
    spec='global'
    if not plot_all:
        if plot_energy_levels: #Plot energy levels
            plt.xlim(-0.6,0)
            plt.ylim(-0.25,0.25)
            spec = 'bound_levels'
            for i in range(1,5):
                plt.vlines(-0.5/i**2, -0.2, 0.2,
                           color='r', label='n='+str(i))
                plt.text(-0.5/i**2, 0.2, 'n='+str(i),
                         horizontalalignment='center', fontsize=14,
                         bbox={'facecolor': 'red', 'alpha': 1, 'pad': 5})
        else: #Plot positive eigenvalues
            plt.xlim(0,5)
            plt.ylim(-0.15,0.15)
            spec='unknown'
    plt.tight_layout()
    plt.show()
    plt.savefig(path+'img/BoundaryTerm_vs_Eigenvalue--'+spec+'.png',
                dpi=200, bbox_inches='tight')

#%% Exercise 5.1
run_Exercise51 = yes_or_no('Run exercise 5.1?')
if run_Exercise51:
    print('\n')
    print('#####################')
    print('#### Exercise 5.1')
    #%%
    print('\n####### Hydrogen')
    # Parameters
    r_range = [0, 50]
    Eigenvalue_range = [-5, 0.]
    
    u0_tol= 0.001
    Eigenvalue_tol = 0.00001
    h = 0.0001
    
    j_max = 500000 
    delta = 0.0001
    
    KS_int_max = 100
    
    verbose = False
    write_data = False
    
    Uniform_Numerov = [True, True] # Uniform grid
    
    hydrogen_eigenvalue, u0, HE, EE, CE = np.abs(pydft.dft.hydrogenatom(r_range,
                                          Eigenvalue_range,
                                          KS_int_max, Eigenvalue_tol,
                                          u0_tol, Uniform_Numerov, h,
                                          j_max, delta,
                                          verbose, write_data, path+'data/'))
    print('## Uniform grid')
    print('u(0): '+str(u0))
    print('Eigenvalue: '+str(hydrogen_eigenvalue))
    print('Absolute error: '+str(abs(0.5-hydrogen_eigenvalue)))
    
    Uniform_Numerov = [False, False] # Non-uniform grid (Runge-Kutta)
    
    hydrogen_eigenvalue, u0, HE, EE, CE = np.abs(pydft.dft.hydrogenatom(r_range,
                                          Eigenvalue_range,
                                          KS_int_max, Eigenvalue_tol,
                                          u0_tol, Uniform_Numerov, h,
                                          j_max, delta,
                                          verbose, write_data, path+'data/'))
    print('\n## Non-uniform grid (Runge-Kutta)')
    print('u(0): '+str(u0))
    print('Eigenvalue: '+str(hydrogen_eigenvalue))
    print('Absolute error: '+str(abs(0.5-hydrogen_eigenvalue)))
    
    Uniform_Numerov = [False, True] # Non-uniform grid (Numerov)
    
    hydrogen_eigenvalue, u0, HE, EE, CE = np.abs(pydft.dft.hydrogenatom(r_range,
                                          Eigenvalue_range,
                                          KS_int_max, Eigenvalue_tol,
                                          u0_tol, Uniform_Numerov, h,
                                          j_max, delta,
                                          verbose, write_data, path+'data/'))
    print('\n## Non-uniform grid (Numerov)')
    print('u(0): '+str(u0))
    print('Eigenvalue: '+str(hydrogen_eigenvalue))
    print('Absolute error: '+str(abs(0.5-hydrogen_eigenvalue)))
    
    print('\n####### Helium')
    r_range = [0, 20]
    
    u0_tol= 0.0001
    Eigenvalue_tol = 0.0001
    SelfCons_tol = 0.0001
    
    j_max = 200000
    h = 0.0001
    
    SelfCons_int_max = 100
    
    Uniform_Numerov = [True, True] # Uniform grid
    
    Energy, eigenvalue, HC, EC, CC = pydft.dft.heliumatom(r_range,
                         Eigenvalue_range, SelfCons_int_max, KS_int_max,
                         Eigenvalue_tol, u0_tol, SelfCons_tol,
                         Uniform_Numerov, h, j_max, delta, verbose,
                         False, 1, 0, 2, 2)
    print('## Uniform grid')
    print('Eigenvalue: '+str(eigenvalue))
    print('Deviation from -0.52: '+str(abs(eigenvalue+0.52)))
    print('Energy: '+str(Energy))
    print('Deviation from -2.72: '+str(abs(Energy+2.72)))
    
    Uniform_Numerov = [False, False] # Non-uniform grid (Runge-Kutta)
    
    Energy, eigenvalue, HC, EC, CC = pydft.dft.heliumatom(r_range,
                         Eigenvalue_range, SelfCons_int_max, KS_int_max,
                         Eigenvalue_tol, u0_tol, SelfCons_tol,
                         Uniform_Numerov, h, j_max, delta, verbose,
                         False, 1, 0, 2, 2)
    print('\n## Non-uniform grid (Runge-Kutta)')
    print('Eigenvalue: '+str(eigenvalue))
    print('Deviation from -0.52: '+str(abs(eigenvalue+0.52)))
    print('Energy: '+str(Energy))
    print('Deviation from -2.72: '+str(abs(Energy+2.72)))
    
    Uniform_Numerov = [False, True] # Non-uniform grid (Numerov)
    
    Energy, eigenvalue, HC, EC, CC = pydft.dft.heliumatom(r_range,
                         Eigenvalue_range, SelfCons_int_max, KS_int_max,
                         Eigenvalue_tol, u0_tol, SelfCons_tol,
                         Uniform_Numerov, h, j_max, delta, verbose,
                         False, 1, 0, 2, 2)
    print('\n## Non-uniform grid (Numerov)')
    print('Eigenvalue: '+str(eigenvalue))
    print('Deviation from -0.52: '+str(abs(eigenvalue+0.52)))
    print('Energy: '+str(Energy))
    print('Deviation from -2.72: '+str(abs(Energy+2.72)))
    
    #%%
    print('\n####### Error Analysis')
    def func(x, a, b, c):
        return a/(x**c) + b
    
    # Parameters
    r_range = [0, 20]
    Eigenvalue_range = [-1, 0.]
    KS_int_max = 1
    Eigenvalue_tol = 0.0001
    u0_tol= 0.0001
    h = 0.0001
    
    Uniform_Numerov = [False, True]
    write_data = False
    
    rp = 50
    j_max = np.arange(1000, 10000, 100, dtype=int)
    delta = rp/j_max
    u0 = np.zeros(np.size(j_max))
    
    for i in range(0,int(np.size(u0))):
        # print(i)
        # print(delta[i]*j_max[i])
        hydrogen_eigenvalue, u0[i], HE, EE, CE = np.abs(pydft.dft.hydrogenatom(r_range,
                                          Eigenvalue_range,
                                          KS_int_max, Eigenvalue_tol,
                                          u0_tol, Uniform_Numerov, h,
                                          j_max[i], delta[i],
                                          verbose, write_data, path+'data/'))
    
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
    plt.show()
    plt.savefig(path+'img/u0vsN.png', dpi=200, bbox_inches='tight')