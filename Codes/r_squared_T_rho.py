import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

r_squared_list = []
rho_list = []

for i in range(0, 6):
    N_H_pre = np.loadtxt('Rates'+str(i)+'/'+'SINK_9pc'+'_N_tri_H1_for_rate.csv', delimiter=',')
    N_Mg_pre = np.loadtxt('Rates'+str(i)+'/'+'SINK_9pc'+'_N_tri_Mg1_for_rate.csv', delimiter=',')
    T, rho = np.loadtxt('Rates'+str(i)+'/'+'Toy_params.txt')
    print(T, rho)
    N_H = [N_H_pre[i] for i in range(len(N_H_pre)) if N_H_pre[i]!=0 and N_Mg_pre[i]!=0]
    N_Mg = [N_Mg_pre[i] for i in range(len(N_Mg_pre)) if N_H_pre[i]!=0 and N_Mg_pre[i]!=0]

    N_H_log = np.log10(N_H)
    N_Mg_log = np.log10(N_Mg)

    r_squared = (stats.linregress(N_H_log, N_Mg_log).rvalue)**2
    print(r_squared)

    r_squared_list.append(r_squared)
    rho_list.append(rho)

    plt.plot(N_H, N_Mg, '.')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('H I Column density [g/cm$^{-2}]$)')
    plt.ylabel('Mg II Column density [g/cm$^{-2}]$)')
    plt.show()

plt.plot(rho_list, r_squared_list, '.')
plt.xlabel('Density of toy model [g/cm$^{-3}]$)')
plt.xscale('log')
plt.ylabel('r_squared')
plt.show()
