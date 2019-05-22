import matplotlib.pyplot as plt
import numpy as np
import math
import shutil

# merging different txt.files into one:
##################################################
#with open('output_file.txt','wb') as wfd:
#    for f in ['data1.txt','data2.txt','data3.txt', ...]:
#        with open(f,'rb') as fd:
#            shutil.copyfileobj(fd, wfd)
##################################################


#txt file with time and temperature values (no rows with nans)
t_1, T_1 = np.genfromtxt('Daten/tdata.txt', usecols=(0, 1), unpack=True) #adjust columns if necessary
#t_1: Time in seconds/unix timestamps, T_1: Temperature in degree celsius




#constant fluence in 1/cm**2
phi = 5*10**(15)

#constant parameters of the N_eff function with a "WE-25k Ohm cm" Diode
N_C0 = 1.1*10**11                   #stable Damage amplitude in 1/cm**3
g_c = 1.58 * 10**(-2)               #acceptor introduction Rate in cm**(-1)
E_y = 1.33*1.6*10**(-19)            #resulting activation Energy in j
k_0y = 1.5 * 10**(15)               #frequency factor in 1/s
g_y = 4.84*10**(-2)                 #average introduction rate in 1/cm
g_a = 1.59 * 10**(-2)               #introduction rate in 1/cm
E_aa = 1.09 * 1.6* 10**(-19)        #activation Energy in j
k_0a = 2.4 *10**(13)                #frequency factor in 1/s
c = 75 * 10**(-14)                  #fit parameter in 1/cm**2

#constant parameters of the damage rate
a_I = 1.23 * 10**(-17)       #amplitude in A/cm
a_0 = -8.9*10**(-17)         #fit parameter in A/cm
k_0I = 1.2 * 10**(13)        #1/s
E_I = 1.11 * 1.6 * 10**(-19) #j
E_I2 = 1.3 * 1.6 * 10**(-19) #j
beta = 3.07*10**(-18)        #A/cm
b_0 = 4.6*10**(-14)          #fit parameter in A*K/cm
T_ref = 322.15               #reference temperature in kelvin

#parameters for the numbers of linear interpolation intervalls n = x * (T-T_min) + y
x_int = 0.05
y_int = 0.2






#constant parameters of the N_eff function with a "WI-4k Ohm cm" diode
#N_C0 = 0.9*10**11                   #stable Damage amplitude in 1/cm**3
#g_c = 1.80 * 10**(-2)               #acceptor introduction Rate in cm**(-1)
#E_y = 1.33 *1.6 *10**(-19)          #resulting activation Energy in j
#g_y = 4.92*10**(-2)                 #average introduction rate 1/cm
#k_0y = 1.5 * 10**(15)*60            #frequency factor in 1/s
#g_a = 2.01 * 10**(-2)               #introduction rate in 1/cm
#E_aa = 1.09 * 1.6* 10**(-19)        #activation Energy in j
#k_0a = 2.4 *10**(13)*60             #frequency factor in 1/s
#c = 31 * 10**(14)                   #fit parameter 1/cm**2
