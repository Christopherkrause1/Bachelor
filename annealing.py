import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
c = 1
tau_1 = 24.1 #min
tau = 880 #min
N_Y_inf = 4.84 * 10**(-2) * 1.4*10**(13) #min/cm
g_C = 1.49 * 10**(-2)  #cm**(-1)
g_a_1 = 1.59 * 10**(-2) #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j
k_0a = 2.4 *10**(13)*60 #1/min

t, phi, T = np.genfromtxt('daten.txt', unpack=True)

stable = N_C0 *(1 - np.exp(-c*phi)) + g_C * phi
shortterm = stable + phi * g_a_1 * np.exp(-t * k_0a * (np.exp(-E_aa/(k_B * T))))
longterm = stable + N_Y_inf * (1- 1/(1 + t/tau))

def N_eff():
    return (N_C0 *(1 - np.exp(-c*phi)) + g_C * phi +              #N_C
    phi * g_a_1 * np.exp(-t * k_0a * (np.exp(-E_aa/(k_B * T)))) + #N_A
    N_Y_inf * (1- 1/(1 + t/tau)))                                 #N_Y

#print(N_eff(1.3*10**(13), 100, 60+273.15))
#x_plot = np.linspace(10, 10000, 5)
#params, covariance_matrix = curve_fit(N_eff, t, phi, T)
#plt.plot(x_plot, N_eff(x_plot, *params), 'k-', label='Anpassungsfunktion', linewidth=0.5)
#print(params)
#print(np.sqrt(np.diag(covariance_matrix)))

plt.gcf().subplots_adjust(bottom=0.18)

plt.semilogx(t , N_eff(), 'r.', label='Änderung N_eff', Markersize=6)
plt.semilogx(t, stable, 'k-', label='stable Damage', Markersize=4)
plt.semilogx(t, shortterm, 'g--', label='short term', Markersize=4)
plt.semilogx(t, longterm, 'b--', label='long term', Markersize=4)
#plt.semilogx(t, shortterm + longterm - stable, 'k--', label='wahrer Verlaus', Markersize=4)
plt.title('annealing effekt')
plt.legend()
plt.grid()
plt.ylim(0, 10*10**(11))
plt.xlabel(r'annealing Zeit / $\mathrm{min}$')
plt.ylabel(r'$\Delta$N_eff /$\mathrm{cm^{-3}} $')
plt.savefig('build/annealing.pdf')
plt.clf()

#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13)*60 #1/min
E_I = 1.11 * 1.6 * 10**(-19) #j
a_0 = -8.9*10**(-17) + 4.6*10**(-14) * 1/T    #A/cm
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #s

def damage():
    return (a_I * np.exp(-t* k_0I* np.exp(-E_I/(k_B*T)))  +     #shortterm
    a_0 - b * np.log(t/t_0))                                    #longterm

plt.gcf().subplots_adjust(bottom=0.18)

plt.semilogx(t , damage(), 'r.', label='damage rate', Markersize=6)
plt.title('current related damage rate')
plt.legend()
plt.grid()
#plt.ylim(0, 10*10**(11))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\alpha  / \mathrm{A cm^{-1}} $')
plt.savefig('build/damage.pdf')
