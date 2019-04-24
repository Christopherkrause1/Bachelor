import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13) #1/s
E_I = 1.11 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #min
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante

t, phi, T, T_2, T_3, T_4 = np.genfromtxt('daten.txt', unpack=True)




def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def gett_I(t, tau_I0, T):
    timediff_I = np.zeros(len(t))
    timediff_I = np.ediff1d(t, to_begin=0)
    tau_I0 = np.roll(tau_I0, shift=1) # shifting array by one to the right
    timediff_I /= tau_I0
    t_I = np.zeros(len(t))
    for i in range(0, len(t)):
        t_I[i] = np.sum(timediff_I[0:i+1])
    return t_I

def a_0(T):                                       #part of the longterm annealing
    return -8.9*10**(-17) + 4.6*10**(-14) * 1/T

def damage(t, T):                                      #damage rate
    tau_I0 = tau_I(T)                                  #tau_I0 = Array [egal, tau_I(T[0]), tau_I(T[1]),...]
    t_I = gett_I(t, tau_I0, T)
    return ((a_I * np.exp(-t_I) + a_0(T)) - b * np.log(t/t_0))

plt.gcf().subplots_adjust(bottom=0.18)

#fig = plt.figure()
#ax1 = fig.add_subplot(111)
##ax1.plot(t, damage(t, T))
#ax1.set_ylabel(r'$ \alpha  / \mathrm{A cm^{-1}}$')
#
#ax2 = ax1.twinx()
##ax2.plot(t, T, 'r-')
#ax2.set_ylabel(r'$T/K$')
#for tl in ax2.get_yticklabels():
#    tl.set_color('k')

plt.semilogx(t/60 , damage(t, T), 'r.', label='damage rate 60°C', Markersize=6)
plt.semilogx(t/60 , damage(t, T_2), 'b.', label='damage rate 21°C', Markersize=6)
plt.semilogx(t/60 , damage(t, T_3), 'k.', label='damage rate 80°C', Markersize=6)
plt.semilogx(t/60 , damage(t, T_4), 'g.', label='damage rate 106°C', Markersize=6)
plt.title('current related damage rate')
plt.legend()
plt.grid()
#plt.ylim(0, 10*10**(11))
plt.xlabel(r't / $\mathrm{min}$')
plt.ylabel(r'$\alpha  / \mathrm{A cm^{-1}} $')
plt.savefig('build/damagekorrektur.pdf')
plt.clf()
