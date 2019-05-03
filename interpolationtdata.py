from einstellungen import *





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


def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)




new_t = np.array(t[0])
new_T = np.array(T_1[0])

T_min = min(T_1)

for i in range(1, len(T_1)):
    n = math.ceil((0.05*abs(T_1[i-1]- T_min) + 0.2))
    print(n)
    for j in range(1, n+1):
        new_T = np.append(new_T, T_1[i-1] + (T_1[i]-T_1[i-1])/(n) * (j))
        new_t = np.append(new_t, t[i-1] + abs(t[i-1]-t[i])/n *j)



fig, ax1 = plt.subplots()
plt.semilogx(new_t/60 , new_T, 'r.', label='interpolierte Temperatur', Markersize=6)
plt.semilogx(t/60 , T_1, 'g.', label='Temperatur', Markersize=6)
ax1.set_ylabel(r"Temperatur / $^{\circ}$C", color = 'red')
ax1.tick_params('y',colors='red')
ax1.set_xlabel("Zeit / min")
ax1.legend(loc=6)


ax2 = ax1.twinx()
plt.semilogx(new_t/60, N_eff(new_t, 5*10**(15), new_T), 'b.', label=r'$\Delta N_{\mathrm{eff}}$ mit interpolation', Markersize=6)
plt.semilogx(t/60, N_eff(t, 5*10**(15), T_1), 'k.', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#plt.semilogx(new_t/60, N_C(5*10**(15))+N_Y(new_t, 5*10**(15), new_T), 'y-', label=r'$\Delta N_{\mathrm{eff}}$ ohne interpolation', Markersize=6)
#plt.semilogx(new_x/60, N_eff(new_x, 5*10**(15), 80), 'k--', label=r'$\Delta N_{\mathrm{eff}}$ für 80°C', Markersize=6)
ax2.set_ylabel(r"$\Delta N_{eff}$ /$\mathrm{cm^{-3}} $",color='blue')
ax2.tick_params('y',colors='blue')
ax1.grid()
ax2.legend(loc='best')
#plt.xlim(0, 2*10**2)
plt.savefig('build/interpolationtdata.pdf')




# Interpolate the data using a linear spline
#y = T_1
#x = t
#new_x = np.zeros(5*len(x)-5)
#for i in range(0,len(x)-1):
#    for j in range(0,5): # mit 3 nur die dazwischen
#        new_x[5*i+j] = x[i] + (x[i+1] - x[i])*(j+1)/5
#new_y = sp.interpolate.interp1d(x, y, kind='linear')(new_x)
#new_x = np.insert(new_x, 0, t[0])
#new_y = np.insert(new_y, 0, T_1[0])
#plt.plot(x/60, y+5, 'bo-')
#plt.plot(new_x/60, new_y, 'ro-')
#plt.title('Using 1D linear Spline Interpolation')
#plt.xlabel(r'$t$ /$\mathrm{min} $')
#plt.ylabel(r'T / $^\mathrm{\circ}$C')
#plt.savefig('build/interpolationtdata.pdf')
#plt.clf()
