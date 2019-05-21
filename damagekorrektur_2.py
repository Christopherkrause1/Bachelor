import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math

#current related damage rate
a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13) #1/s
E_I = 1.11 * 1.6 * 10**(-19) #j
E_I2 = 1.3 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
b_0 = 4.6*10**(-14) #AK/cm

t_0 = 60 #s
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
T_ref = 322.15
t_d, phi, T, T_2, T_3, T_4 = np.genfromtxt('Daten/daten.txt', unpack=True)
t_1, T_1 = np.genfromtxt('Daten/tdata.txt', unpack=True)



def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def gett_I(t, tau_I0, T):
    timediff_I = np.zeros(len(t))
    timediff_I = np.ediff1d(t, to_begin=0)
    tau_I0 = np.roll(tau_I0, shift=1) # shifting array by one to the right
    tau_I1 = tau_I(T)
    timediff_I /= (tau_I0 + tau_I1)/2
    t_I = np.zeros(len(t))
    for i in range(0, len(t)):
        t_I[i] = np.sum(timediff_I[0:i+1])
    return t_I

def a_0(T):                                       #part of the longterm annealing
    return -8.9*10**(-17) + 4.6*10**(-14) * 1/T

def a_02():                                  #Temperatur unabhängige parametrisiserung
    return -8.9*10**(-17) + (b_0/ T_ref)

def theta(T):
    return np.exp(-E_I2/k_B *(1/T - 1/T_ref))


def gett_theta(t, theta_0, T):
    timediff_theta = np.zeros(len(t))
    timediff_theta = np.ediff1d(t, to_begin=0)
    theta_0 = np.roll(theta_0, shift=1) # shifting array by one to the right
    theta_1 = theta(T)
    timediff_theta *= (theta_0 + theta_1)/2
    timediff_theta[0]=10**(-90)
    t_theta = np.zeros(len(t))
    for i in range(0, len(t)):
        t_theta[i] = np.sum(timediff_theta[0:i+1])
    return t_theta

def damage(t, T):                                      #damage rate
    tau_I0 = tau_I(T)                                  #tau_I0 = Array [egal, tau_I(T[0]), tau_I(T[1]),...]
    t_I = gett_I(t, tau_I0, T)
    theta_0 = theta(T)                                  #tau_I0 = Array [egal, tau_I(T[0]), tau_I(T[1]),...]
    t_theta = gett_theta(t, theta_0, T)
    return a_I * np.exp(-t_I) + a_02() - b * np.log(t_theta /t_0)

def interpolation_t(t, T):
    t_int = np.array(t[0])
    T_min = min(T)
    for i in range(1, len(T)):
        n = math.ceil((0.05*abs(T[i-1]- T_min) + 0.2))
        for j in range(1, n+1):
            t_int = np.append(t_int, t[i-1] + abs(t[i-1]-t[i])/n *j)
    return t_int;


def interpolation_T(t, T):
    T_int = np.array(T[0])
    T_min = min(T)
    for i in range(1, len(T)):
        n = math.ceil((0.05*abs(T[i-1]- T_min) + 0.2))
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))
    return T_int


fig, ax1 = plt.subplots()
plt.semilogx(interpolation_t(t_1, T_1)/60 , interpolation_T(t_1, T_1), 'r.', label='interpolierte Temperatur', Markersize=6)
plt.semilogx(t_1/60 , T_1, 'g.', label='Temperatur', Markersize=6)
#plt.ylim(290, 300)
#ax1.bar()
#ax1.scatter()
ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Zeit / min")
ax1.legend(loc='best')


ax2 = ax1.twinx()
#plt.semilogx(t_1/60 , damage(t_1, T_1+273.15), 'b.', label='Schadensrate', Markersize=6)
#plt.semilogx(t_1/60 , damage(t_1, 49+273.15), 'b.', label='Schadensrate 49°C', Markersize=6)
plt.semilogx(interpolation_t(t_1, T_1)/60 , damage(interpolation_t(t_1, T_1), interpolation_T(t_1, T_1)+273.15), 'b.', label='interpolierte Schadensrate', Markersize=6)
plt.semilogx(t_1/60 , damage(t_1, T_1+273.15), 'k.', label='Schadensrate', Markersize=6)
#plt.semilogx(t/60 , damage(t, T_3), 'k.', label='damage rate 80°C', Markersize=6)
#plt.semilogx(t/60 , damage(t, T_4), 'g.', label='damage rate 106°C', Markersize=6)
ax2.set_ylabel(r"$\alpha  / \mathrm{A cm^{-1}} $",color='blue')
ax2.tick_params('y',colors='blue')
#ax2.set_yscale('log')
#ax2.scatter()
plt.ylim(0.1*10**(-16), 0.9*10**(-16))
ax1.grid()
ax2.legend(loc='lower center')



#plt.gcf().subplots_adjust(bottom=0.18)
#plt.semilogx(t/60 , damage(t, T), 'r.', label='damage rate 60°C', Markersize=6)
#plt.semilogx(t/60 , damage(t, T_2), 'b.', label='damage rate 21°C', Markersize=6)
#plt.semilogx(t/60 , damage(t, T_3), 'k.', label='damage rate 80°C', Markersize=6)
#plt.semilogx(t/60 , damage(t, T_4), 'g.', label='damage rate 106°C', Markersize=6)
#plt.title('current related damage rate')
#plt.legend()
#plt.grid()
#plt.xlabel(r't / $\mathrm{min}$')
#plt.ylabel(r'$\alpha  / \mathrm{A cm^{-1}} $')
plt.savefig('build/damagekorrektur_2.pdf')
plt.clf()
