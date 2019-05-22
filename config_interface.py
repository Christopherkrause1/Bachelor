from tkinter import *
import matplotlib.pyplot as plt
import numpy as np
import string
import locale
k_B = 1.38064852 * 10**(-23)    #Boltzmann constant


#constant parameters of the N_eff function with a "WE-25k Ohm cm" Diode
N_C0 = 1.1*10**11               #stable Damage amplitude in 1/cm**3
E_y = 1.33*1.6*10**(-19)        #resulting activation Energy in j
k_0y = 1.5 * 10**(15)           #frequency factor in 1/s
g_c = 1.58 * 10**(-2)           #Acceptor introduction Rate in cm**(-1)
g_a = 1.59 * 10**(-2)           #introduction rate in cm**(-1)
g_y = 4.84*10**(-2)             #cm**(-1)
E_aa = 1.09 * 1.6* 10**(-19)    #activation Energy in j
k_0a = 2.4 *10**(13)            #frequency factor in 1/s
c = 75 * 10**(-14)              #fit parameter in 1/cm**2


#constant parameters of the damage rate
a_I = 1.23 * 10**(-17)         #A/cm
a0 = -8.9*10**(-17)            #fit parameter in A/cm
k_0I = 1.2 * 10**(13)*60       #1/min
E_I = 1.11 * 1.6 * 10**(-19)   #j
beta = 3.07*10**(-18)          #A/cm
b_0 = 4.6*10**(-14)            #fit parameter in A*K/cm
t_0 = 1                        #min
