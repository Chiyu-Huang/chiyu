import numpy as np
import matplotlib.pyplot as plt


Ns = np.loadtxt('Reals/Rates_real_12/'+'SINK_9pc'+'_N_tri_O1_for_rate.csv', delimiter=',')
Ew = np.loadtxt('Reals/Rates_real_12/'+'SINK_9pc'+'_EW_O1_for_rate.csv', delimiter=',')
b = np.loadtxt('Reals/Rates_real_12/'+'SINK_9pc'+'_b.csv', delimiter=',')
rates = np.loadtxt('Reals/Rates_real_12/'+'SINK_9pc'+'_N_rates_sph.csv', delimiter=',')
b = np.loadtxt('ray_vec copy/b_12.csv', delimiter = ',')

plt.hist(b)
plt.show()


# plt.plot([1e17, 1e19], [1e1, 1e3], '-')
plt.plot(Ns.flatten(), Ew.flatten(), '.')
plt.xlabel('Column density by integrating over sightline [cm$^{-2}$]')
# plt.ylabel('Column density by fitting [cm$^{-2}$]')
plt.ylabel('Equivalent width by fitting [A]')
plt.xscale('log')
plt.yscale('log')
# plt.title('Mg II 1240.4')
# plt.xlim(0, 1e13)
# plt.ylim(0, 1e13)
plt.show()


plt.plot(b.flatten(), Ns.flatten(), '.')
plt.xlabel('impact param [kpc]')
# plt.ylabel('Column density by fitting [cm$^{-2}$]')
plt.ylabel('Column density')
plt.xscale('log')
plt.yscale('log')
# plt.title('Mg II 1240.4')
# plt.xlim(0, 1e13)
# plt.ylim(0, 1e13)
plt.show()

plt.plot(b.flatten(), rates.flatten(), '.')
plt.xlabel('impact param [kpc]')
# plt.ylabel('Column density by fitting [cm$^{-2}$]')
plt.ylabel('Outflow rate')
plt.xscale('log')
plt.yscale('log')
# plt.title('Mg II 1240.4')
# plt.xlim(0, 1e13)
# plt.ylim(0, 1e13)
plt.show()

plt.plot(Ns.flatten(), rates.flatten(),  '.')
plt.xlabel('Column density')
# plt.ylabel('Column density by fitting [cm$^{-2}$]')
plt.ylabel('rates')
plt.xscale('log')
plt.yscale('log')
# plt.title('Mg II 1240.4')
# plt.xlim(1e7, 1e14)
# plt.ylim(1e-5, 1e-2)
plt.show()