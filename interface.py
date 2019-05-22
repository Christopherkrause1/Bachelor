from config_interface import *


def N_Y_inf(phi):                                      #longterm annealing amplitude
    return g_y * phi

def tau_Y(T):                                        #Time constant
    return 1/(k_0y *np.exp(-E_y/(k_B*(T+273.15))))


def tau_A(T):                                          #Time constant
    return 1/(k_0a *np.exp(-E_aa/(k_B*(T+273.15))))


def N_C(phi):                                          #stable damage
    return N_C0 *(1 - np.exp(-phi)) + g_c * phi

def N_A(t, phi, T):                                    #shortterm annealing
    return phi * g_a * np.exp(-t/tau_A(T))


def N_Y(t, phi, T):                                    #longterm annealing
    return N_Y_inf(phi) * (1- 1/(1 + t/tau_Y(T)))

def N_eff(t, phi, T):                                #Änderung der Dotierungskonzentration
    return N_C(phi) + N_A(t, phi, T) + N_Y(t, phi, T)



def tau_I(T):                                     #time constant
    return 1/(k_0I* np.exp(-E_I/(k_B*T)))

def a_0(T):                                       #part of the longterm annealing
    return a0 + b_0 * 1/T

def damage(t, T):                          #damage rate
    return a_I * np.exp(-t/tau_I(T)) + a_0(T) - beta * np.log(t/t_0)



master = Tk()
master.title("Annealing effects at constant temperatures for a 'WE-4k' diode")
Label(master, text=r'Insert a time for the annealing to end, a temperature to anneal with and a fluence.').grid(row=0)
Label(master, text=r'Click on "plot" to create a plot of the effective doping concentration.').grid(row=1)
Label(master, text=r'Time [min]').grid(row=2)
Label(master, text="Temperature [°C]").grid(row=3)
Label(master, text=r'Fluence [10^15 /cm^2]').grid(row=4)

t_q = Entry(master)
t_q.grid(row=2, column=1)


T_q = Entry(master)
T_q.grid(row=3, column=1)


phi_q = Entry(master)
phi_q.grid(row=4, column=1)




def plot():
    t_1 = float(t_q.get())
    T_1 = float(T_q.get())
    phi = float(phi_q.get())
    if t_1 < 2:
        t_1=2
        new_t = np.logspace(np.log10(0.1), np.log10(int(t_1)), np.floor(np.log10(int(t_1))*20))
        new_T = T_1
        new_phi = phi * 10**(15)
        plt.semilogx(new_t, N_eff(new_t*60, new_phi, new_T), 'r.')
        plt.grid()
        plt.ylabel(r'$N_{\mathrm{eff}} /\mathrm{cm}^2$', size=25)
        plt.xlabel(r'$Time / $min', size=25)
        plt.show()
    else:
        new_t = np.logspace(np.log10(0.1), np.log10(int(t_1)), np.floor(np.log10(int(t_1))*20))
        new_T = T_1
        new_phi = phi * 10**(15)
        plt.semilogx(new_t, N_eff(new_t*60, new_phi, new_T), 'r.')
        plt.grid()
        plt.ylabel(r'$N_{\mathrm{eff}} /\mathrm{cm}^2$', size=25)
        plt.xlabel(r'$Time / $min', size=25)
        plt.show()


z = Button(master, text="plot", width=10, command=plot)
z.grid(row=5, column=1)


Label(master, text=r'').grid(row=6)
Label(master, text=r'Annealing of the leakage current').grid(row=7)
Label(master, text=r'Insert a time for the annealing to end and a temperature to anneal with.').grid(row=8)
Label(master, text=r'Time [min]').grid(row=9)
Label(master, text="Temperature [°C]").grid(row=10)

t_a = Entry(master)
t_a.grid(row=9, column=1)


T_a = Entry(master)
T_a.grid(row=10, column=1)

def plot_2():
    t_2 = float(t_a.get())
    T_2 = float(T_a.get())
    if t_2 < 2:
        t_2 = 2
        new_t2 = np.logspace(np.log10(0.1), np.log10(int(t_2)), np.floor(np.log10(int(t_2))*20))
        new_T2 = T_2 +273.15
        plt.semilogx(new_t2, damage(new_t2, new_T2), 'r.')
        plt.grid()
        plt.ylabel(r'$\alpha /\mathrm{Acm}^2$', size=25)
        plt.xlabel(r'$Time / $min', size=25)
        plt.show()
    else:
        new_t2 = np.logspace(np.log10(0.1), np.log10(int(t_2)), np.floor(np.log10(int(t_2))*20))
        new_T2 = T_2 +273.15
        plt.semilogx(new_t2, damage(new_t2, new_T2), 'r.')
        plt.grid()
        plt.ylabel(r'$\alpha /\mathrm{Acm}^2$', size=25)
        plt.xlabel(r'$Time / $min', size=25)
        plt.show()



p = Button(master, text="plot", width=10, command=plot_2)
p.grid(row=11, column=1)

mainloop()
