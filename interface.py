from tkinter import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate
import scipy as sp
import math
import string



N_C0 = 1.3*10**11 #1/cm**3
E_y = 1.33*1.6*10**(-19)    #resulting activation Energy
k_0y = 1.5 * 10**(15)   #frequency factor
g_c = 1.49 * 10**(-2)  #cm**(-1)    Acceptor introduction Rate
g_a = 1.59 * 10**(-2) #cm**(-1)   introduction rate
g_y = 5.16*10**(-2)   #cm**(-1)
k_B = 1.38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 * 1.6* 10**(-19) #j   activation Energy
k_0a = 2.4 *10**(13) #1/s   frequency factor

a_I = 1.23 * 10**(-17) #A/cm
k_0I = 1.2 * 10**(13)*60 #1/min
E_I = 1.11 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
t_0 = 1 #min


def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))


def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))


def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    return phi * g_a * np.exp(-t/tau_A(T))


def N_Y(t, phi, T):                                    #longterm annealing
    return N_Y_inf(phi) * (1- 1/(1 + t/tau_Y(T)))

def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)



def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def a_0(T):                                       #part of the longterm annealing
    return -8.9*10**(-17) + 4.6*10**(-14) * 1/T

def damage(t, T):                          #damage rate
    return a_I * np.exp(-t/tau_I(T)) + a_0(T) - b * np.log(t/t_0)



master = Tk()
master.title("Annealing effects at constant temperatures")
Label(master, text=r'Insert a time for the annealing to end, a temperature to anneal with and a fluence.').grid(row=0)
Label(master, text=r'Click on "plot" to create a plot of the effective doping concentration.').grid(row=1)
Label(master, text=r'Time [min]').grid(row=2)
Label(master, text="Temperature [°C]").grid(row=3)
Label(master, text=r'Fluence [10^15 /cm^2]').grid(row=4)

t_q = Entry(master)
t_q.grid(row=2, column=1)


T_q = Entry(master)
T_q.grid(row=3, column=1)


phi_q = Entry(master)
phi_q.grid(row=4, column=1)



def plot():
    t_1 = float(t_q.get())
    T_1 = float(T_q.get())
    phi = float(phi_q.get())
    new_t = np.array(0)
    for i in range(0, int(t_1+1)):
        for j in range(0, 60):
            new_t = np.append(new_t, 60*i+j)
    np.delete(new_t, 0)
    new_T = T_1
    new_phi = phi * 10**(15)
    plt.semilogx(new_t, N_eff(new_t, new_phi, new_T), 'r.')
    plt.grid()
    plt.ylabel(r'$N_{\mathrm{eff}} /\mathrm{cm}^2$')
    plt.xlabel(r'$Time / $min')
    plt.show()


z = Button(master, text="plot", width=10, command=plot)
z.grid(row=5, column=1)


Label(master, text=r'').grid(row=6)
Label(master, text=r'Annealing of the leakage current').grid(row=7)
Label(master, text=r'Insert a time for the annealing to end and a temperature to anneal with.').grid(row=8)
Label(master, text=r'Time [min]').grid(row=9)
Label(master, text="Temperature [°C]").grid(row=10)

t_a = Entry(master)
t_a.grid(row=9, column=1)


T_a = Entry(master)
T_a.grid(row=10, column=1)

def plot_2():
    t_1 = float(t_a.get())
    T_1 = float(T_a.get())
    new_t2 = np.array(0)
    for i in range(0, int(t_1)):
        for j in range(0, 60):
            new_t2 = np.append(new_t2, 60*i+j)
    np.delete(new_t2, 0)
    print(new_t2)
    new_t2[0] = 0.5
    new_T2 = T_1 +273.15
    plt.semilogx(new_t2, damage(new_t2, new_T2), 'r.')
    plt.grid()
    plt.ylabel(r'$\alpha /\mathrm{Acm}^2$')
    plt.xlabel(r'$Time / $min')
    plt.show()



p = Button(master, text="plot", width=10, command=plot_2)
p.grid(row=11, column=1)

mainloop()
