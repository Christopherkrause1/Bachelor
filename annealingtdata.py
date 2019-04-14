import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)*60   #frequnecy factor
N_Y_inf = 4.84 * 10**(-2) * 1.4*10**(13) #min/cm
g_c = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a = 1.59 * 10**(-2) #cm**(-1)   introduction rate
g_y = 5.16*10**(-2)   #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13)*60 #1/min   frequnecy factor


t, T = np.genfromtxt('tdata.txt', unpack=True)


#Änderung der effektiven Dotierungskonzentration für WE-25k$\Omega$cm


def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                          #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*T)))

def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*T)))

def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    return phi * g_a * np.exp(-t/tau_A(T))

def N_Y(t, phi, T):                                    #longterm annealing
    return N_Y_inf(phi) * (1- 1/(1 + t/tau_Y(T)))


def N_eff(t, phi, T):                                  #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)




plt.gcf().subplots_adjust(bottom=0.18)
plt.semilogx(t/60, N_eff(t/60, 5*10**(15), T+273.15), 'r.', label='Änderung N_eff R1', Markersize=6)
plt.semilogx(t/60, N_eff(t/60, 5*10**(15), 80+273.15), 'b.', label='Änderung N_eff 80°C', Markersize=6)
plt.title('Annealingeffekt für R1')
plt.legend()
plt.grid()
#plt.ylim(1.52*10**(14), 1.56*10**(14))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
plt.savefig('build/annealingtdata.pdf')
plt.clf()
