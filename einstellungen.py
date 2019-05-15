import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate
import scipy as sp
import math

#constant parameters of the N_eff function
N_C0 = 1.3*10**11                   #stable Damage amplitude in 1/cm**3
E_y = 1.33*1.6*10**(-19)            #resulting activation Energy in j
k_0y = 1.5 * 10**(15)               #frequency factor in 1/s
g_c = 1.49 * 10**(-2)               #Acceptor introduction Rate in cm**(-1)
g_a = 1.59 * 10**(-2)               #introduction rate in 1/cm
g_y = 5.16*10**(-2)                 #1/cm
E_aa = 1.09 * 1.6* 10**(-19)        #activation Energy in j
k_0a = 2.4 *10**(13)                #frequency factor in 1/s
phi = 5*10**(15)                    #constant fluence in 1/cm**2



#txt file with time and temperature values
t, T_1 = np.genfromtxt('Daten/tdata.txt', unpack=True) # t: Zeit, T_1: Temperatur
t_2, T_2 =np.genfromtxt('Daten/mareike_annealing.txt', usecols=(0, 16), unpack=True)
t_3, T_3 =np.genfromtxt('Daten/mareike_annealing2.txt', usecols=(0, 16), unpack=True)
t_4, T_4 =np.genfromtxt('Daten/20190306_mareike_annealing_R3.txt', usecols=(0, 16), unpack=True)
t_2 -= t_2[0]
t_3 -= t_3[0]
t_4 -= t_4[0]
