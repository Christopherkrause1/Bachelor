from tkinter import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate
import scipy as sp
import math
import string


#self.b2.place(x=200, y=150)
#self.lbl3.place(x=100, y=200)
#self.t3.place(x=200, y=200)
    #def add(self):
    #    t_q=int(self.t1.get())
    #    T_q=int(self.t2.get())
    #    self.t3.insert(END, str(result))
    #def sub(self, event):
    #    num1=int(self.t1.get())
    #    num2=int(self.t2.get())
    #    result=num1-num2
    #    self.t3.insert(END, str(result))



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






master = Tk()
master.title("Annealing effect of the doping concentration")
Label(master, text=r'Insert a time for the annealing to end, a temperature to anneal with and a fluence.').grid(row=0)
Label(master, text=r'Time [min]').grid(row=1)
Label(master, text="Temperature [°C]").grid(row=2)
Label(master, text=r'Fluence [10^15 /cm^2]').grid(row=3)

t_q = Entry(master)
t_q.grid(row=1, column=1)
#t_q.pack()

T_q = Entry(master)
T_q.grid(row=2, column=1)
#T_q.pack()
#T_q.insert(0, '0')

phi_q = Entry(master)
phi_q.grid(row=3, column=1)
#phi_q.pack()

#t_q.focus_set()
#T_q.focus_set()
#phi_q.focus_set()


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


b = Button(master, text="plot", width=10, command=plot)
b.grid(row=4, column=1)
#b.pack()

mainloop()
