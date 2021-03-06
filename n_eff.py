from configuration import *

t_1 -= t_1[0]                                          #converts unix time stamps to seconds
k_B = 1.38064852 * 10**(-23)                           #Boltzmann constant in J/K

print(t_extra)

def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                          #Time constant of the longterm annealing
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, T_s_1, T):                     #creating an approximation for t/tau_Y  (T_s ist the Temperature shifted to the right)
    timediff_Y = np.zeros(len(t))            #creates an array of zeros to work with
    timediff_Y = np.ediff1d(t, to_begin=0)   #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_s_1 = np.roll(T_s_1, shift=1)          #shifting tau_Y array by one to the right
    T_n = T                                  #Temperature not shifted to the right
    timediff_Y /= tau_Y((T_s_1+T_n)/2)       #dividing each element by the mean of 2 adjacent Temperature elements (inside tau)
    t_Y = np.zeros(len(t))                   #create an array of zeros to work with
    for i in range(0, len(t)):
        t_Y[i] = np.sum(timediff_Y[0:i+1])   #writes in each element the sum of the cooresponding timediff_Y values
    return t_Y                               #now looks like [0, 0 + t[1]-t[0]/(tau_Y[0] + tau_Y[1])/2, ...]


def tau_A(T):                                #Time constant of the short term annealing
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, T_s_2, T):                    #creating an approximation for t/tau_A
    timediff_A = np.zeros(len(t))           #creates an array of zeros to work with
    timediff_A = np.ediff1d(t, to_begin=0)  #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_s_2 = np.roll(T_s_2, shift=1)         #shifting tau_A array by one to the right
    T_n = T                                 #Temperature not shifted to the right
    timediff_A /= tau_A((T_s_2+T_n)/2)      #dividing each element by the mean of 2 adjacent tau_A elements
    t_A = np.zeros(len(t))                  #create an array of zeros to work with
    for i in range(0, len(t)):
        t_A[i] = np.sum(timediff_A[0:i+1])  #writes in each element the sum of the cooresponding timediff_Y values
    return t_A                              #now looks like [0, 0 + t[1]-t[0]/(tau_A[0] + tau_A[1])/2, ...]


def N_C(phi):                                        #stable damage
    return N_C0 *(1 - np.exp(-c * phi)) + g_c * phi

def N_A(t, phi, T):                                  #shortterm annealing
    T_s_2 = T                                        #assigning array for the approximation
    t_A = gett_A(t, T_s_2, T)                        #assigning name to the new approximated time
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                  #longterm annealing
    T_s_1 = T                                        #assigning array for the approximation
    t_Y = gett_Y(t, T_s_1, T)                        #assigning name to the new approximated time
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def interpolation_t(t, T):                                            #linear interpolation of the time for more data
    t_int = np.array(t[0])                                            #create new time starting with arrays of zeros
    T_max = max(T_1)                                                   #number of intervalls depend on maximum temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T_1[i-1]- T_1[i])/(T_max-(T_1[i-1]+T_1[i])/2 +y_int) + z_int))             #function for the number of intervalls
        for j in range(1, n+1):
            t_int = np.append(t_int, t[i-1] + abs(t[i-1]-t[i])/n *j)  #new interpolated times (includes initial times)
    return t_int;


def interpolation_T(t, T):                                            #linear interpolation of the temperature for more data
    T_int = np.array(T[0])                                            #create new temperature starting with arrays of zeros
    T_max = max(T_1)                                                    #number of intervalls depend on minimal temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T_1[i-1]- T_1[i])/(T_max-(T_1[i-1]+T_1[i])/2 +y_int) + z_int))           #function for the number of intervalls
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))#new interpolated temperatures (includes initial ones)
    return T_int




def N_eff(t, phi, T):  #change of the doping concentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)




#function that plots n_eff
def plot_N_eff(t, phi, T):
    if 'show_temperature_curve' in globals():
        fig, ax1 = plt.subplots()   #temperature axis
        plt.semilogx(interpolation_t(t, T)/60 , interpolation_T(t, T), 'g.', label='interpolated temperature', Markersize=6)
        plt.semilogx(t/60 , T, 'r.', label='temperature', Markersize=6)
        ax1.set_ylabel(r"Temperature / $^{\circ}$C", color = 'red', size=13)
        ax1.tick_params('y',colors='red')
        ax1.set_xlabel("Time / min", size=13)
        ax1.legend(loc='best')
        ax2 = ax1.twinx()          #n_eff axis
        plt.semilogx(interpolation_t(t, T)/60, N_eff(interpolation_t(t, T), phi, interpolation_T(t, T)), 'b.', label=r'interpolated $\Delta N_{\mathrm{eff}}$', Markersize=6)
        plt.semilogx(t/60, N_eff(t, phi, T), 'k.', label=r'$\Delta N_{\mathrm{eff}}$', Markersize=6)
        if 'T_const' in globals():
            plt.semilogx(t/60, N_eff(t, phi, T_const), 'c.', label=r'$\Delta N_{\mathrm{eff}} @$'+ str(T_const) + '°C', Markersize=6)
        ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue', size=13)
        ax2.tick_params('y',colors='blue')
        ax1.grid()
        ax2.legend(loc='lower left')
        plt.axvline(x=4.917)
        plt.axvline(x=12.4)
        ax2.axhline(y=2.2112*10**(13))
        plt.show()
        #plt.savefig('build/mareike_R9_2.pdf')
        plt.clf()

    else:
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.semilogx(interpolation_t(t, T)/60, N_eff(interpolation_t(t, T), phi, interpolation_T(t, T)), 'b.', label=r'interpolated $\Delta N_{\mathrm{eff}}$', Markersize=6)
        plt.semilogx(t/60, N_eff(t, phi, T), 'k.', label=r'$\Delta N_{\mathrm{eff}}$', Markersize=6)
        if 'T_const' in globals():
            plt.semilogx(t/60, N_eff(t, phi, T_const), 'c.', label=r'$\Delta N_{\mathrm{eff}} @$'+ str(T_const) + '°C', Markersize=6)
        plt.xlabel("Time / min", size=13)
        plt.ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue', size=13)
        plt.tick_params('y',colors='blue')
        plt.grid()
        plt.legend(loc='best')
        plt.show()
        plt.clf()

plot_N_eff(t_1, phi, T_1)    #function of the doping concentration
