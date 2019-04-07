import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

N_C0 = 1.3*10**11 #1/cm**3
c = 1
tau_1 = 24.1 #min
tau = 1060 #min
N_Y_inf = 5.16 * 10**(-2) * 1.4*10**(13) #min/cm
g_C = 1.49 * 10**(-2)  #cm**(-1)
g_a_1 = 1.81 * 10**(-2) #cm**(-1)
k_B = 1,38064852 * 10**(-23) #Boltzmann Konstante
E_aa = 1.09 ##!/usr/bin/env python
k_0a = 2.4 *10**(-13) #1/s

def N_eff(phi, t, T):
    #tau_a = 1/(k_0a * (np.exp(-E_aa/k_B*T))
    return (N_C0 *(1 - np.exp(-c*phi)) + g_C * phi +              #N_C
    phi * g_a_1 * np.exp(-t/(k_0a * (np.exp(-E_aa/(k_B * T))))) + #N_A
    N_Y_inf * (1- 1/(1 + t/tau)))                                 #N_Y

print(N_eff(1.3*10**(13), 100, 5))
