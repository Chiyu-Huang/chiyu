
#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import matplotlib.animation as ani
import Codes.outflow_rate as outflow_rate
import Codes.fitting as fitting
import Codes.Random_rays as Random_rays
from Codes.Ray_maker import find_rate, to_cl, ew, to_kpc
from datetime import datetime
yt.set_log_level(50)


N_H = []
EW_H = []
N_O = []
EW_O = []
N_Mg = []
EW_Mg = []
N_C = []
EW_C = []
spe = ['H I 1216', 'O VI 1038', 'C IV 1548', 'Mg II 1240']
spe = ['H I 1216', 'H I 973']
spe = ['O VI 1032']
spe = ['H I 1216']


_rho = 3e-25
_Temp = 1e4
_vel_r = 200*1000

ds, info = outflow_rate.load_data(203, file = 'SINK_9pc')

print('------------------')
print(ds.all_data()[('gas', 'theta')])


print('aaaaaaaaaaaa')

#%%

_rho = 8e-25
_Temp = 1e4
_vel_r = 7000*1000


_rho = 3e-25
_Temp = 1e4
_vel_r = 200*1000

def _temperature_mu_2(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc')
    rv = []
    mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 41))
    rv = np.where(mask, _Temp, 1e6)
    return yt.YTArray(rv, 'K')
def _density_2(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc')
    rv = data['gas', 'density'].in_units('g/cm**3').value
    def den(r):
        return _rho
    mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 41))
    rv = np.where(mask, den(r), 1e-29)
    return yt.YTArray(rv, 'g/cm**3')
def _vel(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc')
    rv = data['gas', 'radial_velocity'].in_units('cm/s').value
    T = data[('gas', 'temperature_mu')]
    mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 41))
    rv = np.where(mask, _vel_r, 0)
    return yt.YTArray(rv, 'cm/s')

# ds.add_field(
#     ("gas", "temperature_mu"),
#     sampling_type="cell",
#     function=_temperature_mu_2,
#     units='K',
#     force_override = True
# )



# ds.add_field(
#     ("gas", "density_2"),
#     sampling_type="cell",
#     function=_density_2,
#     units='g/cm**3',
#     force_override = True
# )


# ds.add_field(
#     ("gas", "vel_rad"),
#     sampling_type="cell",
#     function=_vel,
#     units='cm/s',
#     force_override = True
# )

#%%

spe = ['O VI 1032']

start_l = np.loadtxt('ray_vec copy/start_0.csv', delimiter=',')
end_l = np.loadtxt('ray_vec copy/end_0.csv', delimiter=',')
b_l = np.loadtxt('ray_vec copy/b_0.csv', delimiter=',')

for i in range (len(start_l)):

    print(i+1)

    start, end, b = start_l[i], end_l[i], b_l[i]

    # start, end = to_cl(np.array([-15.09493946+75, -14.68964895+75 ,-75. +75       ])) ,to_cl(np.array([-15.09493946 +75,-14.68964895 +75, 75.   +75     ]))
    # start, end = to_cl(np.array([  9.06661739+75 ,-15.40089089+75, -75. +75       ])) ,to_cl(np.array([  9.06661739 +75,-15.40089089 +75, 75.   +75     ]))
    # start, end = np.array([0.6,0.5,0]),np.array([0.6,0.5,1])
    print(to_kpc(start)-75, to_kpc(end)-75)
    print('impact param = ', to_kpc(b))

    

    ew_fit, ew_raw, Ns, Ns_tri, ray = ew(ds, start, end,
            spe,\
            plot = True, crit = 'ew', b = 1.1,\
            MilkyWay = False, QSO = False, Noise = False, plot_fit = False)
    rate = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b),\
            mode = 'ray_temp', ray = ray, start_end = (start, end), spe = 'O VI')
    D = end-start
    cos_theta = D[2]/np.linalg.norm(D)
    ad = ray.all_data()
    vel = np.average(ad[('gas', 'velocity_los')])/cos_theta

    # print(vel)
    print('rate = ', rate)
    print('N = ', Ns_tri[spe[0]])
    EW_H.append(rate)
    # EW_H.append(ew_raw[spe[0]])
    N_H.append(Ns_tri[spe[0]])
    # EW_H.append(ew_raw[spe[0]])
    plt.plot(N_H, EW_H, '.')
    plt.plot(N_H[-1], EW_H[-1], 'x')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    # N_O.append(Ns_tri[spe[1]])
    # EW_O.append(ew_raw[spe[1]])

    # N_C.append(Ns_tri[spe[2]])
    # EW_C.append(ew_raw[spe[2]])

    # N_Mg.append(Ns_tri[spe[3]])
    # EW_Mg.append(ew_raw[spe[3]])


# plt.plot(N_H, EW_H, '.')
# plt.plot(N_H[-1], EW_H[-1], 'x')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# plt.plot(N_O, EW_O, '.')
# plt.plot(N_O[-1], EW_O[-1], 'x')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# plt.plot(N_C, EW_C, '.')
# plt.plot(N_C[-1], EW_C[-1], 'x')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

# plt.plot(N_Mg, EW_Mg, '.')
# plt.plot(N_Mg[-1], EW_Mg[-1], 'x')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()
# %%
