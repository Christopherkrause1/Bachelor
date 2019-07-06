import matplotlib.pyplot as plt
import numpy as np
from configuration import *





k_B = 1.38064852 * 10**(-23)             #Boltzmann constant in J/K
t_0 = 60                                 #s



def tau_I(T):                            #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def gett_I(t, T_s, T):                      #creates an approximation for t/tau_I
    timediff_I = np.zeros(len(t))           #creates an array of zeros to work with
    timediff_I = np.ediff1d(t, to_begin=0)  #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_s = np.roll(T_s, shift=1)             #shifting tau_I array by one to the right
    T_n = T
    timediff_I /= tau_I((T_s+T_n)/2)        #dividing each element by the mean of 2 adjacent tau_I elements
    t_I = np.zeros(len(t))                  #create an array of zeros to work with
    for i in range(0, len(t)):
        t_I[i] = np.sum(timediff_I[0:i+1])  #writes in each element the sum of the cooresponding timediff_I values
    return t_I                               #now looks like [0, 0 + t[1]-t[0]/(tau_I[0] + tau_I[1])/2, ...]


def a_02():                                  #temperature independent
    return a_0 + (b_0/ T_ref)

def theta(T):                                #scaling factor for the time
    return np.exp(-E_I2/k_B *(1/T - 1/T_ref))


def gett_theta(t, T_t, T):                            #creates an approximation for theta*t
    timediff_theta = np.zeros(len(t))                 #creates an array of zeros to work with
    timediff_theta = np.ediff1d(t, to_begin=0)        #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_t = np.roll(T_t, shift=1)                       #shifting theta array by one to the right
    T_n = T
    timediff_theta *= theta((T_n+T_t)/2)              #dividing each element by the mean of 2 adjacent theta elements
    timediff_theta[0]=10**(-90)                       #to avoid dividing by zero in the logarithm
    t_theta = np.zeros(len(t))                        #create an array of zeros to work with
    for i in range(0, len(t)):
        t_theta[i] = np.sum(timediff_theta[0:i+1])    #writes in each element the sum of the cooresponding timediff values
    return t_theta                                   #now looks like [0, 0 + t[1]-t[0]*(theta[0] + theta[1])/2, ...]

def damage(t, T):                                     #damage rate
    T_s = T                                           #assigning array for the t/tau_I approximation
    t_I = gett_I(t, T_s, T)                           #assigning name to the new approximated time
    T_t = T                                           #assigning array theta*t for the approximation
    t_theta = gett_theta(t, T_t, T)                   #assigning name for the approximation
    return a_I * np.exp(-t_I) + a_02() - beta * np.log(t_theta /t_0)

fig, ax1 = plt.subplots()
plt.semilogx(t_1/60 , T_1, 'r.', label='Temperatur', Markersize=6)
ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Zeit / min")
ax1.legend(loc='upper left')


ax2 = ax1.twinx()
plt.semilogx(t_1/60, damage(t_1, T_1+273.15), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ von R1', Markersize=6)
ax2.set_ylabel(r"$\alpha $ /$\mathrm{A cm^{-1}} $",color='blue')
ax2.tick_params('y',colors='blue')
plt.ylim(0, 1*10**(-16))
ax1.grid()
ax2.legend(loc='best')
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylim(0, 1*10**(-16))
plt.savefig('build/damagekorrektur.pdf')
#plt.show()
plt.clf()
