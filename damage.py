import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13)*60 #1/min
E_I = 1.11 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #min
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante

t, phi, T, T_2, T_3, T_4 = np.genfromtxt('daten.txt', unpack=True)




def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def a_0(T):                                       #part of the longterm annealing
    return -8.9*10**(-17) + 4.6*10**(-14) * 1/T

def damage(t, T):                                #damage rate
    return ((a_I * np.exp(-t/tau_I(T))) + a_0(T) - b * np.log(t/t_0))

plt.gcf().subplots_adjust(bottom=0.18)

plt.semilogx(t , damage(t, T), 'r.', label='damage rate 60°C', Markersize=6)
plt.semilogx(t , damage(t, T=T_2), 'b.', label='damage rate 21°C', Markersize=6)
plt.semilogx(t , damage(t, T=T_3), 'k.', label='damage rate 80°C', Markersize=6)
plt.semilogx(t , damage(t, T=T_4), 'g.', label='damage rate 106°C', Markersize=6)
plt.title('current related damage rate')
plt.legend()
plt.grid()
#plt.ylim(0, 10*10**(11))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\alpha  / \mathrm{A cm^{-1}} $')
plt.savefig('build/damage.pdf')
# damage(t=80, T=60°C) = 4.0236277*10**(-17) A/cm, gemessener Wert 4.01 +/- 0.04
plt.clf()
