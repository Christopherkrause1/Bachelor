import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
c = 1
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)*60   #frequnecy factor
#tau_y = k_0y *np.exp(E_y/(k_B*T)) #time constant
N_Y_inf = 4.84 * 10**(-2) * 1.4*10**(13) #min/cm
g_C = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a_1 = 1.59 * 10**(-2) #cm**(-1)   introduction rate
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13)*60 #1/min   frequnecy factor

t, phi, T, T_2, T_3, T_4 = np.genfromtxt('daten.txt', unpack=True)

#Änderung der effektiven Dotierungskonzentration für WE-25k$\Omega$cm

stable = N_C0 *(1 - np.exp(-c*phi)) + g_C * phi
shortterm = stable + phi * g_a_1 * np.exp(-t * k_0a * (np.exp(-E_aa/(k_B * T))))
longterm = stable + N_Y_inf * (1- 1/(1 + t*k_0y *np.exp(-E_y/(k_B*T))))

def N_eff(T):
    return (N_C0 *(1 - np.exp(-c*phi)) + g_C * phi +              #N_C
    phi * g_a_1 * np.exp(-t * k_0a * (np.exp(-E_aa/(k_B * T)))) + #N_A
    4.84*10**(-2)* phi * (1- 1/(1 + t*k_0y *np.exp(-E_y/(k_B*T)))))          #N_Y


plt.gcf().subplots_adjust(bottom=0.18)

plt.semilogx(t , N_eff(T), 'r.', label='Änderung N_eff 60°C', Markersize=6)
plt.semilogx(t , N_eff(T=T_2), 'b.', label='Änderung N_eff 21°C', Markersize=6)
plt.semilogx(t , N_eff(T=T_3), 'g.', label='Änderung N_eff 80°C', Markersize=6)
plt.semilogx(t, stable, 'k-', label='stable Damage', Markersize=4)
plt.semilogx(t, shortterm, 'g--', label='short term', Markersize=4)
plt.semilogx(t, longterm, 'b--', label='long term', Markersize=4)
#plt.semilogx(t, shortterm + longterm - stable, 'k--', label='wahrer Verlauf', Markersize=4)
plt.title('Annealingeffekt für WE-25k$\Omega$cm')
plt.legend()
plt.grid()
plt.ylim(0, 11*10**(11))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
#plt.savefig('build/annealing.pdf')
plt.clf()




#dotierungs Konzentration bei WI-4kOhm cm

N_C0 = 1.3*10**11 #1/cm**3
c = 1
tau_1 = 24.1 #min
tau = 1010 #min
N_Y_inf = 4.92 * 10**(-2) * 1.4*10**(13) #min/cm
g_C = 1.49 * 10**(-2)  #cm**(-1)
g_a_1 = 2.01 * 10**(-2) #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j
k_0a = 2.4 *10**(13)*60 #1/min


def N_eff_2(T):
    return (N_C0 *(1 - np.exp(-c*phi)) + g_C * phi +              #N_C
    phi * g_a_1 * np.exp(-t * k_0a * (np.exp(-E_aa/(k_B * T)))) + #N_A
    4.92*10**(-2)* phi * (1- 1/(1 + t*k_0y *np.exp(-E_y/(k_B*T)))))                                 #N_Y


plt.gcf().subplots_adjust(bottom=0.18)

plt.semilogx(t , N_eff_2(T), 'r.', label='Änderung N_eff 60°C', Markersize=6)
plt.semilogx(t , N_eff_2(T=T_2), 'b.', label='Änderung N_eff 21°C', Markersize=6)
plt.semilogx(t , N_eff_2(T=T_3), 'g.', label='Änderung N_eff 80°C', Markersize=6)
plt.semilogx(t, stable, 'k-', Markersize=4)
plt.semilogx(t, shortterm, 'g--', Markersize=4)
plt.semilogx(t, longterm, 'b--', Markersize=4)
plt.title('Annealingeffekt für WI-4k$\Omega$cm')
plt.legend()
plt.grid()
plt.ylim(0, 11*10**(11))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
plt.savefig('build/annealing2.pdf')
plt.clf()
