import matplotlib.pyplot as plt
import numpy as np
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

#constant parameters of the  damage rate
a_I = 1.23 * 10**(-17) #A/cm
a_0 = -8.9*10**(-17) #A/cm
k_0I = 1.2 * 10**(13) #1/s
E_I = 1.11 * 1.6 * 10**(-19) #j
E_I2 = 1.3 * 1.6 * 10**(-19) #j
b = 3.07*10**(-18)    #A/cm
b_0 = 4.6*10**(-14) #AK/cm
t_0 = 60 #s
T_ref = 322.15 #reference temperature in kelvin


#txt file with time and temperature values
t_1, T_1 = np.genfromtxt('Daten/tdata.txt', unpack=True) # t_1: Zeit, T_1: Temperatur
t_2, T_2 =np.genfromtxt('Daten/mareike_annealing.txt', usecols=(0, 16), unpack=True)
t_3, T_3 =np.genfromtxt('Daten/mareike_annealing2.txt', usecols=(0, 16), unpack=True)
t_4, T_4 =np.genfromtxt('Daten/20190306_mareike_annealing_R3.txt', usecols=(0, 16), unpack=True)
t_2 -= t_2[0]
t_3 -= t_3[0]
t_4 -= t_4[0]


#merging different txt.files into one:
##################################################
#import shutil
#with open('output_file.txt','wb') as wfd:
#    for f in ['seg1.txt','seg2.txt','seg3.txt']:
#        with open(f,'rb') as fd:
#            shutil.copyfileobj(fd, wfd)
##################################################
