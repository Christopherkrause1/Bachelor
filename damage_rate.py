from configuration import *
t_1 -= t_1[0]                            #converts unix time stamps to seconds
k_B = 1.38064852 * 10**(-23)             #Boltzmann constant in J/K
t_0 = 60                                 #in s



def tau_I(T):                               #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def gett_I(t, T_s, T):                      #creates an approximation for t/tau_I
    timediff_I = np.zeros(len(t))           #creates an array of zeros to work with
    timediff_I = np.ediff1d(t, to_begin=0)  #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_s = np.roll(T_s, shift=1)             #shifting tau_I array by one to the right
    T_n = T
    timediff_I /= tau_I((T_s+T_n)/2)        #dividing each element by the mean of 2 adjacent tau_I elements
    t_I = np.zeros(len(t))                  #create an array of zeros to work with
    for i in range(0, len(t)):
        t_I[i] = np.sum(timediff_I[0:i+1])  #writes in each element the sum of the corresponding timediff_I values
    return t_I                              #now looks like [0, 0 + t[1]-t[0]/(tau_I[0] + tau_I[1])/2, ...]


def a_02():                                  #temperature independent
    return a_0 + (b_0/ T_ref)

def theta(T):                                #scaling factor for the time
    return np.exp(-E_I2/k_B *(1/T - 1/T_ref))


def gett_theta(t, T_t, T):                            #creates an approximation for theta*t
    timediff_theta = np.zeros(len(t))                 #creates an array of zeros to work with
    timediff_theta = np.ediff1d(t, to_begin=0)        #creates array = [0, t[1]-t[0], t[2]-t[1], ...]
    T_t = np.roll(T_t, shift=1)                       #shifting theta array by one to the right
    T_n = T
    timediff_theta *= theta((T_n+T_t)/2)              #dividing each element by the mean of 2 adjacent theta elements
    timediff_theta[0]=10**(-90)                       #to avoid dividing by zero in the logarithm
    t_theta = np.zeros(len(t))                        #create an array of zeros to work with
    for i in range(0, len(t)):
        t_theta[i] = np.sum(timediff_theta[0:i+1])    #writes in each element the sum of the cooresponding timediff values
    return t_theta                                    #now looks like [0, 0 + t[1]-t[0]*(theta[0] + theta[1])/2, ...]

def damage(t, T):                                     #damage rate
    T_s = T                                           #assigning array for the t/tau_I approximation
    t_I = gett_I(t, T_s, T)                           #assigning name to the new approximated time
    T_t = T                                           #assigning array theta*t for the approximation
    t_theta = gett_theta(t, T_t, T)                   #assigning name for the approximation
    return a_I * np.exp(-t_I) + a_02() - beta * np.log(t_theta /t_0)

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
    T_max = max(T_1)                                                    #number of intervalls depend on maximum temperature
    for i in range(1, len(T)):
        n = math.ceil((x_int*abs(T_1[i-1]- T_1[i])/(T_max-(T_1[i-1]+T_1[i])/2 +y_int) + z_int))           #function for the number of intervalls
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))#new interpolated temperatures (includes initial ones)
    return T_int                                         #new interpolated temperatures (includes initial temperatures)


#function that plots the damage rate
def plot_damage_rate(t, T):
    if 'show_temperature_curve' in globals():
        fig, ax1 = plt.subplots()
        plt.semilogx(interpolation_t(t, T)/60 , interpolation_T(t, T), 'g.', label='Interpolierte Temperatur', Markersize=6)
        plt.semilogx(t/60 , T, 'r.', label='Temperatur', Markersize=6)
        ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red', size=13)
        ax1.tick_params('y',colors='red')
        ax1.set_xlabel("Zeit / min", size=13)
        ax1.legend(loc='best')
        ax2 = ax1.twinx() #second axis (damage rate)
        plt.semilogx(interpolation_t(t, T)/60 , damage(interpolation_t(t, T), interpolation_T(t, T)+273.15), 'b.', label='Interpolierte Schadensrate', Markersize=6)
        plt.semilogx(t/60 , damage(t, T+273.15), 'k.', label='Schadensrate', Markersize=6)

        if 'T_const' in globals():
            plt.semilogx(t/60 , damage(t/60, T_const+273.15), 'c.', label=r'$\alpha @$'+ str(T_const) + '°C', Markersize=6)
        ax2.set_ylabel(r"$\alpha  / \mathrm{A cm^{-1}} $",color='blue', size=13)
        ax2.tick_params('y',colors='blue')
        #ax2.set_ylim(None, 1.2*damage(interpolation_t(t, T), interpolation_T(t, T)+273.15)[1])
        ax2.set_ylim(0, 1*10**(-16))
        ax1.grid()
        ax2.legend(loc='lower center')
        plt.savefig('build/damage_interpolation.pdf')
        #plt.show()
        plt.clf()

    else:
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.semilogx(interpolation_t(t, T)/60, damage(interpolation_t(t, T), interpolation_T(t, T)+273.15), 'b.', Markersize=6)
        plt.semilogx(t/60, damage(t, T+273.15), 'k.', label=r'Schadensrate', Markersize=6)
        plt.semilogx(t_extra[0]/60 ,  4.44689013*10**(-17), 'co', label='Gemessene Schadensrate', Markersize=6)
        plt.semilogx(t_extra[1]/60 ,  4.31653104*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[2]/60 ,  4.23412643*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[3]/60 ,  4.14873766*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[5]/60 ,  4.04255188*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[8]/60 ,  3.93018074*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[13]/60 , 3.81504309*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[23]/60 , 3.62997850*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[38]/60 , 3.55766500*10**(-17), 'co', Markersize=6)
        plt.semilogx(t_extra[59]/60 , 3.44813677*10**(-17), 'co', Markersize=6)
        if 'T_const' in globals():
            plt.semilogx(t/60, damage(t/60, T_const+273.15), 'c.', label=r'$\alpha @$'+ str(T_const) + '°C', Markersize=6)
        plt.xlabel("Zeit / min", size=13)
        plt.xlim(10, None)
        plt.legend(loc='best')
        plt.ylabel(r"$\alpha$ /$\mathrm{A cm^{-1}} $", size=13)
        plt.ylim(None, 1.2*damage(interpolation_t(t, T), interpolation_T(t, T)+273.15)[1])
        #plt.tick_params('y',colors='blue')
        plt.grid()
        #plt.show()
        plt.savefig('build/damage_P_3.pdf')
        plt.clf()


#shows the plot of the given data
plot_damage_rate(t_1, T_1)
