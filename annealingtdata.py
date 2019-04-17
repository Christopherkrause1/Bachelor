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
k_0a = 2.4 *10**(13) #1/min   frequency factor

t, T = np.genfromtxt('tdata.txt', unpack=True)   #R1 daten
t_unix, T3 = np.genfromtxt('2018-09-22_11_21_40_Annealingtest_1950.txt', usecols=(0, 1), unpack=True)  #unix daten
t_s = t_unix-t_unix[0]          #t_s = vergangene Zeit in Sekunden


#Änderung der effektiven Dotierungskonzentration für R1


def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, tau_Y):
    tmp_Y = np.zeros(len(t))
    tmp_Y = np.ediff1d(t, to_begin=0)
    tau_Y = np.roll(tau_Y, shift=1) # shifting array by one to the right
    tmp_Y /= tau_Y
    t_Y = np.zeros(len(t))
    for i in range(0, len(t)):
        t_Y[i] = np.sum(tmp_Y[0:i+1])
    return t_Y
tau_Y = tau_Y(T)
t_Y = gett_Y(t, tau_Y)

def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, tau_A):
    tmp_A = np.zeros(len(t))
    tmp_A = np.ediff1d(t, to_begin=0)
    tau_A = np.roll(tau_A, shift=1) # shifting array by one to the right, zweiter Term soll ja durch tau[1] geteilt werden
    tmp_A /= tau_A
    t_A = np.zeros(len(t))
    for i in range(0, len(t)):
        t_A[i] = np.sum(tmp_A[0:i+1])
    return t_A
tau_A = tau_A(T)
t_A = gett_A(t, tau_A)

def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t_A, phi):                                    #shortterm annealing
    return phi * g_a * np.exp(-t_A)


def N_Y(t_Y, phi):                                    #longterm annealing
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def N_eff(t, phi):                                  #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t_A, phi) + N_Y(t_Y, phi)



#print(N_A(t_A, 5*10**(15)))
#print(N_Y(t_Y, 5*10**(15)))
#print((t[1]-t[0])/tau_Y[0] + (t[2]-t[1])/tau_Y[1])
plt.gcf().subplots_adjust(bottom=0.18)
plt.semilogx(t/60, N_eff(t, 5*10**(15)), 'r.', label='Änderung N_eff R1', Markersize=6)
plt.semilogx(t/60, N_C(5*10**(15))+N_A(t_A, 5*10**(15)) + N_Y(t_Y, 5*10**(15)), 'k-', label='Änderung N_eff R1', Markersize=6)
plt.semilogx(t/60, N_C(5*10**(15))+N_A(t_A, 5*10**(15)), 'b-', label='Änderung N_A', Markersize=6)
plt.semilogx(t/60, N_C(5*10**(15))+N_Y(t_Y, 5*10**(15)), 'g-', label='Änderung N_A', Markersize=6)
#plt.semilogx(t/60, N_eff(t/60, 5*10**(15)), 'b.', label='Änderung N_eff 80°C', Markersize=6)
plt.title('Annealingeffekt für R1')
plt.legend()
plt.grid()
#plt.ylim(1.52*10**(14), 1.56*10**(14))
plt.xlabel(r'Zeit / $\mathrm{min}$')
plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
plt.savefig('build/annealingtdata.pdf')
plt.clf()



#Änderung der effektiven Dotierungskonzentration für Diode mit Unix Zeiten

t = t_s
T = T3
#plt.gcf().subplots_adjust(bottom=0.18)
#plt.semilogx(t/60, N_eff(t/60, 1*10**(15)), 'r.', label='Änderung N_eff', Markersize=6)
##plt.semilogx(t_s/60, N_eff(t_s/60, 1*10**(15)), 'b.', label='Änderung N_eff 60°C', Markersize=6)
#plt.title('Annealingeffekt für Unix Timestamps')
#plt.legend()
#plt.grid()
##plt.ylim(1.52*10**(14), 1.56*10**(14))
#plt.xlabel(r'Zeit / $\mathrm{min}$')
#plt.ylabel(r'$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $')
#plt.savefig('build/annealingunix.pdf')
#plt.clf()
