import interpolationmareike as inter
import damagekorrektur_2 as dmg
from einstellungen import *

#########################################################
# DONT TOUCH!!!!!
t_1 -= t_1[0]  #converts unix time stamps to seconds

inter.plot_N_eff(t_1, phi, T_1)    #function of the doping concentration
dmg.plot_damage_rate(t_1, T_1)     #function of the damage rate
#########################################################








#t_2, T_2 =np.genfromtxt('Daten/mareike_annealing.txt', usecols=(0, 16), unpack=True)
#t_3, T_3 =np.genfromtxt('Daten/mareike_annealing2.txt', usecols=(0, 16), unpack=True)
#t_4, T_4 =np.genfromtxt('Daten/20190306_mareike_annealing_R3.txt', usecols=(0, 16), unpack=True)
