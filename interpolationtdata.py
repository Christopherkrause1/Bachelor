import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)   #frequency factor
g_c = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a = 1.59 * 10**(-2) #cm**(-1)   introduction rate
g_y = 5.16*10**(-2)   #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13) #1/s   frequency factor




def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, tau_Y0, T):
    timediff_Y = np.zeros(len(t))
    timediff_Y = np.ediff1d(t, to_begin=0)
    tau_Y0 = np.roll(tau_Y0, shift=1) # shifting array by one to the right
    timediff_Y /= tau_Y0
    t_Y = np.zeros(len(t))
    for i in range(0, len(t)):
        t_Y[i] = np.sum(timediff_Y[0:i+1])
    return t_Y


def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, tau_A0, T):                              #sum of time differences divided by tau(T)
    timediff_A = np.zeros(len(t))
    timediff_A = np.ediff1d(t, to_begin=0)
    tau_A0 = np.roll(tau_A0, shift=1) # shifting array by one to the right, t[1]-t[0] soll ja durch tau[0] geteilt werden
    timediff_A /= tau_A0
    t_A = np.zeros(len(t))
    for i in range(0, len(t)):
        t_A[i] = np.sum(timediff_A[0:i+1])
    return t_A


def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    tau_A0 = tau_A(T)                                  #tau_A0 = Array [egal, tau_A(T[0]), tau_A(T[1]),...]
    t_A = gett_A(t, tau_A0, T)                         #Vektor t_1 - t_0/tau_A(0)
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                    #longterm annealing
    tau_Y0 = tau_Y(T)                                  #tau_Y0 = Array [egal, tau_Y(T[0]), tau_Y(T[1]),...]
    t_Y = gett_Y(t, tau_Y0, T)                         #Vektor t_1 - t_0/tau_Y(0)
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)


#Änderung der effektiven Dotierungskonzentration für R1
t, T_1 = np.genfromtxt('tdata.txt', unpack=True)   #R1 daten
plt.semilogx(t/60 , T_1, 'r.', label='Temperatur', Markersize=6)
plt.title('Interpolation der Temperaturen')
plt.legend()
plt.grid()
plt.xlabel(r'$t$ /$\mathrm{min} $')
plt.ylabel(r'T / $^\mathrm{\circ}$C')
plt.savefig('build/interpolationtdata.pdf')
plt.clf()

#T_max = max(T_1)
#for i in range(1, len(T_1)):
#    if T_1[i]-T_1[i-1] <= 0.05*(T_1[i-1]- T_max) + 0.02*273.15:
#        n = 0.05*(T_1[i-1]- T_max) + 0.02*273.15
#        t_I = (t[i]-t[i-1])/n
#        print(i)
#        print(t_I)
#print('---------------------------')
#T_max = max(T_1)
#for i in range(1, len(T_1)):
#    print(T_1[i]-T_1[i-1])


# Generate some random data
import scipy.interpolate
import scipy as sp
y = T_1
x = t

# Interpolate the data using a linear spline
new_x = np.zeros(4*len(x)-4)
for i in range(0,len(x)-1):
    for j in range(0,4): # mit 3 nur die dazwischen
        new_x[4*i+j] = x[i] + (x[i+1] - x[i])*(j+1)/4
new_y = sp.interpolate.interp1d(x, y, kind='linear')(new_x)
new_x = np.insert(new_x, 0, 0)
new_y = np.insert(new_y, 0, T_1[0])
print(len(new_x))
print(new_y)
plt.plot(x/60, y+5, 'bo-')
plt.plot(new_x/60, new_y, 'ro-')
plt.title('Using 1D linear Spline Interpolation')
plt.xlabel(r'$t$ /$\mathrm{min} $')
plt.ylabel(r'T / $^\mathrm{\circ}$C')
plt.savefig('build/interpolationtdata.pdf')
plt.clf()
#plt.subplot(2,1,2)
