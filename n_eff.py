from einstellungen import *
k_B = 1.38064852 * 10**(-23)        #Boltzmann constant in J/K



def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))

def gett_Y(t, tau_Y0, T):
    timediff_Y = np.zeros(len(t))
    timediff_Y = np.ediff1d(t, to_begin=0)
    tau_Y0 = np.roll(tau_Y0, shift=1) # shifting array by one to the right
    tau_Y1 = tau_Y(T)
    timediff_Y /= (tau_Y0+ tau_Y1)/2
    t_Y = np.zeros(len(t))
    for i in range(0, len(t)):
        t_Y[i] = np.sum(timediff_Y[0:i+1])
    return t_Y


def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))

def gett_A(t, tau_A0, T):                              #sum of time differences divided by tau(T)
    timediff_A = np.zeros(len(t))
    timediff_A = np.ediff1d(t, to_begin=0)
    tau_A0 = np.roll(tau_A0, shift=1) # shifting array by one to the right, t[1]-t[0] soll ja durch tau[0] geteilt werden
    tau_A1 = tau_A(T)
    timediff_A /= (tau_A0 + tau_A1)/2
    t_A = np.zeros(len(t))
    for i in range(0, len(t)):
        t_A[i] = np.sum(timediff_A[0:i+1])
    return t_A


def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    tau_A0 = tau_A(T)                         #tau_A0 = Array [egal, tau_A(T[0]), tau_A(T[1]),...]
    t_A = gett_A(t, tau_A0, T)                   #Vektor t_1 - t_0/tau_A(0)
    return phi * g_a * np.exp(-t_A)


def N_Y(t, phi, T):                                    #longterm annealing
    tau_Y0 = tau_Y(T)                             #tau_Y0 = Array [egal, tau_Y(T[0]), tau_Y(T[1]),...]
    t_Y = gett_Y(t, tau_Y0, T)                         #Vektor t_1 - t_0/tau_Y(0)
    return N_Y_inf(phi) * (1- 1/(1 + t_Y))


def interpolation_t(t, T):
    t_int = np.array(t[0])
    T_min = min(T)
    for i in range(1, len(T)):
        n = math.ceil((0.05*abs(T[i-1]- T_min) + 0.2))
        for j in range(1, n+1):
            t_int = np.append(t_int, t[i-1] + abs(t[i-1]-t[i])/n *j)
    return t_int;


def interpolation_T(t, T):
    T_int = np.array(T[0])
    T_min = min(T)
    for i in range(1, len(T)):
        n = math.ceil((0.05*abs(T[i-1]- T_min) + 0.2))
        for j in range(1, n+1):
            T_int = np.append(T_int, T[i-1] + (T[i]-T[i-1])/(n) * (j))
    return T_int




def N_eff(t, phi, T):  #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)








#erster Datensatz
def plot_N_eff(t, phi, T):
    fig, ax1 = plt.subplots()
    plt.semilogx(interpolation_t(t, T)/60 , interpolation_T(t, T), 'r.', label='interpolierte Temperatur', Markersize=6)
    plt.semilogx(t/60 , T, 'g.', label='Temperatur', Markersize=6)
    ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
    ax1.tick_params('y',colors='red')
    ax1.set_xlabel("Zeit / min")
    ax1.legend(loc=6)


    ax2 = ax1.twinx()
    plt.semilogx(interpolation_t(t, T)/60, N_eff(interpolation_t(t, T), phi, interpolation_T(t, T)), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
    plt.semilogx(t/60, N_eff(t, phi, T), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
    ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
    ax2.tick_params('y',colors='blue')
    ax1.grid()
    ax2.legend(loc='best')
    plt.savefig('build/interpolationmareike.pdf')
    plt.clf()



#zweiter Datensatz
#fig, ax1 = plt.subplots()
#plt.semilogx(interpolation_t(t_3, T_3)/60 , interpolation_T(t_3, T_3), 'r.', label='interpolierte Temperatur', Markersize=6)
#plt.semilogx(t_3/60 , T_3, 'g.', label='Temperatur', Markersize=6)
#ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Zeit / min")
#ax1.legend(loc=6)
#
#
#ax2 = ax1.twinx()
#plt.semilogx(interpolation_t(t_3, T_3)/60, N_eff(interpolation_t(t_3, T_3), phi, interpolation_T(t_3, T_3)), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
#plt.semilogx(t_3/60, N_eff(t_3, 5*10**(15), T_3), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#ax1.grid()
#ax2.legend(loc='best')
#plt.savefig('build/interpolationmareike2.pdf')
#plt.clf()
#
#
##erster Datensatz R3
#fig, ax1 = plt.subplots()
#plt.semilogx(interpolation_t(t_4, T_4)/60 , interpolation_T(t_4, T_4), 'r.', label='interpolierte Temperatur', Markersize=6)
#plt.semilogx(t_4/60 , T_4, 'g.', label='Temperatur', Markersize=6)
#ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
#ax1.tick_params('y',colors='red')
#ax1.set_xlabel("Zeit / min")
#ax1.legend(loc=6)
#
#
#ax2 = ax1.twinx()
#plt.semilogx(interpolation_t(t_4, T_4)/60, N_eff(interpolation_t(t_4, T_4), phi, interpolation_T(t_4, T_4)), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
#plt.semilogx(t_4/60, N_eff(t_4, 5*10**(15), T_4), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
#ax2.tick_params('y',colors='blue')
#ax1.grid()
#ax2.legend(loc='best')
#plt.savefig('build/interpolationmareikeR3_1.pdf')
#plt.clf()
