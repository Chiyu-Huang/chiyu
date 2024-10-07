#%%
import yt
import os
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
import scipy.stats as stats


def k_inclination(i):

    # i = 90 => edge on
    # i = 0 => face on
    if i == 90:
        z = 0
        x = np.random.random()
        y = np.random.random()
        return [x, y, z]
    elif i < 90:
        z = 1
        r = np.tan(i*np.pi/180)
        phi = np.random.random()*2*np.pi
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        return [x, y, z]
    elif i > 90:
        z = -1
        r = np.tan(i*np.pi/180)
        phi = np.random.random()*2*np.pi
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        return [x, y, z]


def gen_start(inclination, nums):
    _start_l = []
    _end_l = []
    _b_l = []
    for i in range(nums):
        # k = [0, 0, 1]
        # k = np.random.random(3)-0.5
        if str(inclination) == 'rand':
            # print('random')
            k = np.random.random(3)-0.5
        else:
            k = k_inclination(inclination)
            # print(k)
        start, end, b = Random_rays.start_end(k, maxb = 40, minb = 5)
        print(to_kpc(b))
        _start_l.append(start)
        _end_l.append(end)
        _b_l.append(b)
        np.savetxt('ray_vec/start_'+str(inclination)+'.csv', _start_l, delimiter=',')
        np.savetxt('ray_vec/end_'+str(inclination)+'.csv', _end_l, delimiter=',')
        np.savetxt('ray_vec/b_'+str(inclination)+'.csv', _b_l, delimiter=',')
    return None


def find_EW_rates (file, bmax, rho = None, temp = None, vel = None, QSO = False,\
                    ite = None, _density_2 = None, _temperature_mu_2 = None,\
                    _vel = None, mode_alt = None, inclination = 'rand'):
    print(ite)
    rates = []
    bs = []
    nums = np.arange(122, 204, step = 3, dtype = int)
    nums = [203]
    
    print(nums)
    for sims in file:
        print(sims)
        ew_O1 = []
        ew_C1 = []
        ew_Mg1 = []
        ew_H1 = []
        ew_O2 = []
        ew_C2 = []
        ew_Mg2 = []
        ew_H2 = []

        N_O1 = []
        N_C1 = []
        N_Mg1 = []
        N_H1 = []
        N_O2 = []
        N_C2 = []
        N_Mg2 = []
        N_H2 = []

        N_tri_O1 = []
        N_tri_C1 = []
        N_tri_Mg1 = []
        N_tri_H1 = []
        N_tri_O2 = []
        N_tri_C2 = []
        N_tri_Mg2 = []
        N_tri_H2 = []
        non_det_O1 = []
        non_det_O2 = []
        non_det_C1 = []
        non_det_C2 = []
        non_det_H1 = []
        non_det_H2 = []
        non_det_Mg1 = []
        non_det_Mg2 = []

        impact_params = []

        rates_sph = []
        rates_pla = []
        rates_ray_H = []
        rates_ray_O = []
        rates_ray_C = []
        rates_ray_Mg = []

        starts = []
        ends = []
        vel_outs = []

        for ts in nums:
            print(ts)
            ds, info = outflow_rate.load_data(ts, file = sims)
            
            if _density_2 is not None:
                print('aj')
                ds.add_field(
                    ("gas", "density_2"),
                    sampling_type="cell",
                    function=_density_2,
                    units='g/cm**3',
                    force_override = True
                )
            if _temperature_mu_2 is not None:
                print('aj')
                ds.add_field(
                    ("gas", "temperature_mu"),
                    sampling_type="cell",
                    function=_temperature_mu_2,
                    units='K',
                    force_override = True
                )
            if _vel is not None:
                print('aj')
                ds.add_field(
                    ("gas", "vel_rad"),
                    sampling_type="cell",
                    function=_vel,
                    units='km/s',
                    force_override = True
                )

            start_l = np.loadtxt('ray_vec copy/start_'+str(inclination)+'.csv', delimiter=',')
            end_l = np.loadtxt('ray_vec copy/end_'+str(inclination)+'.csv', delimiter=',')
            b_l = np.loadtxt('ray_vec copy/b_'+str(inclination)+'.csv', delimiter=',')
            for i in range (len(b_l)):
                start_t = datetime.now()
                print(i+1, '/', len(b_l))
                start, end, b = start_l[i], end_l[i], b_l[i]

                print('impact param = ', to_kpc(b))
                ew_fit, ew_raw, Ns, Ns_tri, ray = ew(ds, start, end,
                        ['O VI 1032', 'Mg II 1240',\
                        'C IV 1548', 'H I 1216',\
                        'O VI 1038', 'Mg II 1240.4',\
                        'C IV 1551', 'H I 973'],\
                        plot = False, crit = 'ew', b = 1.1,\
                        MilkyWay = False, QSO = QSO, Noise = False, plot_fit = False)
                D = end-start
                ad = ray.all_data()
                vel_out = np.average(np.abs(ad[('gas', 'velocity_los')].value))

                rate_sp = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'sphere')
                rate_pl = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'plate')
                # rate_sp = 0
                # print('sphe', rate_sp)
                # rate_pl = 0
                rate_ra_H = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'ray_temp', ray = ray, start_end = (start, end), spe = 'H I')
                rate_ra_O = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'ray_temp', ray = ray, start_end = (start, end), spe = 'O VI')
                rate_ra_C = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'ray_temp', ray = ray, start_end = (start, end), spe = 'C IV')
                rate_ra_Mg = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'ray_temp', ray = ray, start_end = (start, end), spe = 'Mg II')
                


                vel_outs.append(vel_out)
                ew_O1.append(ew_raw['O VI 1032'])
                non_det_O1.append(ew_raw['O VI 1032_is_upper_lim'])
                ew_Mg1.append(ew_raw['Mg II 1240'])
                non_det_Mg1.append(ew_raw['Mg II 1240_is_upper_lim'])
                ew_C1.append(ew_raw['C IV 1548'])
                non_det_C1.append(ew_raw['C IV 1548_is_upper_lim'])
                ew_H1.append(ew_raw['H I 1216'])
                non_det_H1.append(ew_raw['H I 1216_is_upper_lim'])
                ew_O2.append(ew_raw['O VI 1038'])
                non_det_O2.append(ew_raw['O VI 1038_is_upper_lim'])
                ew_Mg2.append(ew_raw['Mg II 1240.4'])
                non_det_Mg2.append(ew_raw['Mg II 1240.4_is_upper_lim'])
                ew_C2.append(ew_raw['C IV 1551'])
                non_det_C2.append(ew_raw['C IV 1551_is_upper_lim'])
                ew_H2.append(ew_raw['H I 973'])
                non_det_H2.append(ew_raw['H I 973_is_upper_lim'])

                N_O1.append(Ns['O VI 1032'])
                N_Mg1.append(Ns['Mg II 1240'])
                N_C1.append(Ns['C IV 1548'])
                N_H1.append(Ns['H I 1216'])
                N_O2.append(Ns['O VI 1038'])
                N_Mg2.append(Ns['Mg II 1240.4'])
                N_C2.append(Ns['C IV 1551'])
                N_H2.append(Ns['H I 973'])

                N_tri_O1.append(Ns_tri['O VI 1032'])
                N_tri_Mg1.append(Ns_tri['Mg II 1240'])
                N_tri_C1.append(Ns_tri['C IV 1548'])
                N_tri_H1.append(Ns_tri['H I 1216'])
                N_tri_O2.append(Ns_tri['O VI 1038'])
                N_tri_Mg2.append(Ns_tri['Mg II 1240.4'])
                N_tri_C2.append(Ns_tri['C IV 1551'])
                N_tri_H2.append(Ns_tri['H I 973'])

                rates_sph.append(rate_sp)
                rates_pla.append(rate_pl)
                rates_ray_H.append(rate_ra_H)
                rates_ray_O.append(rate_ra_O)
                rates_ray_C.append(rate_ra_C)
                rates_ray_Mg.append(rate_ra_Mg)

                impact_params.append(to_kpc(b))

                starts.append(start)
                ends.append(end)
                print('plate time = ', datetime.now()-start_t)
                
                

        name = str(sims)+'_Reals/Rates' + str(ite) + '_' + str(inclination) +'/'
        print(name)
        os.makedirs(name, exist_ok=True)
        # print(rates_ray_H)
        np.savetxt(name+sims+'_vel_out.csv', vel_outs, delimiter=',')
        np.savetxt(name+sims+'_EW_C1_for_rate.csv', ew_C1, delimiter=',')
        np.savetxt(name+sims+'_EW_O1_for_rate.csv', ew_O1, delimiter=',')
        np.savetxt(name+sims+'_EW_H1_for_rate.csv', ew_H1, delimiter=',')
        np.savetxt(name+sims+'_EW_Mg1_for_rate.csv', ew_Mg1, delimiter=',')
        np.savetxt(name+sims+'_EW_C2_for_rate.csv', ew_C2, delimiter=',')
        np.savetxt(name+sims+'_EW_O2_for_rate.csv', ew_O2, delimiter=',')
        np.savetxt(name+sims+'_EW_H2_for_rate.csv', ew_H2, delimiter=',')
        np.savetxt(name+sims+'_EW_Mg2_for_rate.csv', ew_Mg2, delimiter=',')
        
        np.savetxt(name+sims+'_non_det_O1_for_rate.csv', non_det_O1, delimiter=',')
        np.savetxt(name+sims+'_non_det_O2_for_rate.csv', non_det_O2, delimiter=',')
        np.savetxt(name+sims+'_non_det_C1_for_rate.csv', non_det_C1, delimiter=',')
        np.savetxt(name+sims+'_non_det_C2_for_rate.csv', non_det_C2, delimiter=',')
        np.savetxt(name+sims+'_non_det_H1_for_rate.csv', non_det_H1, delimiter=',')
        np.savetxt(name+sims+'_non_det_H2_for_rate.csv', non_det_H2, delimiter=',')
        np.savetxt(name+sims+'_non_det_Mg1_for_rate.csv', non_det_Mg1, delimiter=',')
        np.savetxt(name+sims+'_non_det_Mg2_for_rate.csv', non_det_Mg2, delimiter=',')

        np.savetxt(name+sims+'_N_C1_for_rate.csv', N_C1, delimiter=',')
        np.savetxt(name+sims+'_N_O1_for_rate.csv', N_O1, delimiter=',')
        np.savetxt(name+sims+'_N_H1_for_rate.csv', N_H1, delimiter=',')
        np.savetxt(name+sims+'_N_Mg1_for_rate.csv', N_Mg1, delimiter=',')
        np.savetxt(name+sims+'_N_C2_for_rate.csv', N_C2, delimiter=',')
        np.savetxt(name+sims+'_N_O2_for_rate.csv', N_O2, delimiter=',')
        np.savetxt(name+sims+'_N_H2_for_rate.csv', N_H2, delimiter=',')
        np.savetxt(name+sims+'_N_Mg2_for_rate.csv', N_Mg2, delimiter=',')

        np.savetxt(name+sims+'_N_tri_C1_for_rate.csv', N_tri_C1, delimiter=',')
        np.savetxt(name+sims+'_N_tri_O1_for_rate.csv', N_tri_O1, delimiter=',')
        np.savetxt(name+sims+'_N_tri_H1_for_rate.csv', N_tri_H1, delimiter=',')
        np.savetxt(name+sims+'_N_tri_Mg1_for_rate.csv', N_tri_Mg1, delimiter=',')
        np.savetxt(name+sims+'_N_tri_C2_for_rate.csv', N_tri_C2, delimiter=',')
        np.savetxt(name+sims+'_N_tri_O2_for_rate.csv', N_tri_O2, delimiter=',')
        np.savetxt(name+sims+'_N_tri_H2_for_rate.csv', N_tri_H2, delimiter=',')
        np.savetxt(name+sims+'_N_tri_Mg2_for_rate.csv', N_tri_Mg2, delimiter=',')

        np.savetxt(name+sims+'_N_rates_sph.csv', rates_sph, delimiter=',')
        np.savetxt(name+sims+'_N_rates_pla.csv', rates_pla, delimiter=',')
        np.savetxt(name+sims+'_N_rates_ray_H.csv', rates_ray_H, delimiter=',')
        np.savetxt(name+sims+'_N_rates_ray_O.csv', rates_ray_O, delimiter=',')
        np.savetxt(name+sims+'_N_rates_ray_C.csv', rates_ray_C, delimiter=',')
        np.savetxt(name+sims+'_N_rates_ray_Mg.csv', rates_ray_Mg, delimiter=',')

        np.savetxt(name+sims+'_b.csv', impact_params, delimiter=',')

        np.savetxt(name+sims+'_start.csv', start_l, delimiter=',')
        np.savetxt(name+sims+'_end.csv', end_l, delimiter=',')
        
        np.savetxt(name+'Toy_params.txt', [str(temp), str(rho), str(vel), str(mode_alt)], delimiter=',', fmt='%s')

        print('FIN.')
    return 0


def plate (file, bmax, rho = None, temp = None, vel = None, QSO = False,\
                    ite = None, _density_2 = None, _temperature_mu_2 = None,\
                    _vel = None, mode_alt = None, inclination = 'rand'):
    print(ite)
    rates = []
    bs = []
    nums = np.arange(122, 204, step = 3, dtype = int)
    nums = [203]
    
    print(nums)
    for sims in file:
        print(sims)
        rates_pla = []

        for ts in nums:
            print(ts)
            ds, info = outflow_rate.load_data(ts, file = sims)
            
            if _density_2 is not None:
                print('aj')
                ds.add_field(
                    ("gas", "density_2"),
                    sampling_type="cell",
                    function=_density_2,
                    units='g/cm**3',
                    force_override = True
                )
            if _temperature_mu_2 is not None:
                print('aj')
                ds.add_field(
                    ("gas", "temperature_mu"),
                    sampling_type="cell",
                    function=_temperature_mu_2,
                    units='K',
                    force_override = True
                )
            if _vel is not None:
                print('aj')
                ds.add_field(
                    ("gas", "vel_rad"),
                    sampling_type="cell",
                    function=_vel,
                    units='km/s',
                    force_override = True
                )

            start_l = np.loadtxt('ray_vec copy/start_'+str(inclination)+'.csv', delimiter=',')
            end_l = np.loadtxt('ray_vec copy/end_'+str(inclination)+'.csv', delimiter=',')
            b_l = np.loadtxt('ray_vec copy/b_'+str(inclination)+'.csv', delimiter=',')
            for i in range (len(b_l)):
                start_t = datetime.now()
                print(i+1, '/', len(b_l))
                start, end, b = start_l[i], end_l[i], b_l[i]

                print('impact param = ', to_kpc(b))

                rate_pl = outflow_rate.outflow_rate(ds, to_kpc(b), -to_kpc(b), mode = 'plate')
                rates_pla.append(rate_pl)
                print('plate time = ', datetime.now()-start_t)
                
                

        name = str(sims)+'_Reals/Rates' + str(ite) + '_' + str(inclination) +'/'
        print(name)
        os.makedirs(name, exist_ok=True)
       
        np.savetxt(name+sims+'_N_rates_pla.csv', rates_pla, delimiter=',')

        print('FIN.')
    return 0


def Ew_plot_rates(file, mode, QSO = False, file_read = None, file_rates = None, inclin = None, species = ['H1', 'O1', 'C1', 'Mg1']):


    if inclin is not None:
        file_read = file + file_read+str(inclin)+'/'
    else:
        file_read = file + file_read

    if mode == 'sphere':
        # print('sph')
        rates_file_H = '_N_rates_sph.csv'
        rates_file_O = '_N_rates_sph.csv'
        rates_file_C = '_N_rates_sph.csv'
        rates_file_Mg = '_N_rates_sph.csv'
    elif mode == 'plate':
        # print('pla')
        rates_file_H = '_N_rates_pla.csv'
        rates_file_O = '_N_rates_pla.csv'
        rates_file_C = '_N_rates_pla.csv'
        rates_file_Mg = '_N_rates_pla.csv'

    elif mode == 'ray':
        # print('ray')
        rates_file_H = '_N_rates_ray_H.csv'
        rates_file_O = '_N_rates_ray_O.csv'
        rates_file_C = '_N_rates_ray_C.csv'
        rates_file_Mg = '_N_rates_ray_Mg.csv'



    for line in species:
        if line == 'H1':
            if mode == 'ray':
                lower_thresh = 2e11
                rates_file_H = '_N_rates_ray_H.csv'
            elif mode == 'sphere':
                lower_thresh = 5e11
                rates_file_H = '_N_rates_sph.csv'
            elif mode == 'plate':
                lower_thresh = 1e11
                rates_file_H = '_N_rates_pla.csv'
            rates = np.loadtxt(file_read+file+rates_file_H, delimiter=',')

        elif line == 'O1':
            if mode == 'ray':
                lower_thresh = 2e7
                rates_file_O = '_N_rates_ray_O.csv'
            elif mode == 'sphere':
                lower_thresh =5e12
                rates_file_O = '_N_rates_sph.csv'
            elif mode == 'plate':
                lower_thresh = 1e8
                rates_file_O = '_N_rates_pla.csv'
            rates = np.loadtxt(file_read+file+rates_file_O, delimiter=',')

        elif line == 'C1':
            if mode == 'ray':
                lower_thresh = 4e7
                rates_file_C = '_N_rates_ray_C.csv'
            elif mode == 'sphere':
                lower_thresh =1e12
                rates_file_C = '_N_rates_sph.csv'
            elif mode == 'plate':
                lower_thresh = 1e8
                rates_file_C = '_N_rates_pla.csv'
                
            rates = np.loadtxt(file_read+file+rates_file_C, delimiter=',')

        elif line == 'Mg1':
            if mode == 'ray':
                lower_thresh = 1e7
                rates_file_Mg = '_N_rates_ray_Mg.csv'
            elif mode == 'sphere':
                lower_thresh =1e7
                rates_file_Mg = '_N_rates_sph.csv'
            elif mode == 'plate':
                lower_thresh = 1e7
                rates_file_Mg = '_N_rates_pla.csv'

            rates = np.loadtxt(file_read+file+rates_file_Mg, delimiter=',')

        b = (np.loadtxt(file_read + file + '_b.csv', delimiter=','))
        O1dat = (np.loadtxt(file_read+file+'_N_tri_'+line+'_for_rate.csv', delimiter = ','))
        vel_out = np.loadtxt(file_read+file+'_vel_out.csv', delimiter = ',')
        O1dat2 = []
        rates2 = []
        b2 = []
        print(len(O1dat))
        for i in range(len(b)):
            # if b[i] > 0:
            if O1dat[i] >= 1:
                O1dat2.append(O1dat[i])
                rates2.append(rates[i])
                b2.append(b[i])

        O1dat_log = []
        rates_log = []
        counts_N = 0
        counts_Rates = 0
        for i in range(len(rates2)):
            if O1dat2[i] >= lower_thresh and rates2[i]>0:
                rates_log.append(np.log10(rates2[i]))
                O1dat_log.append(np.log10(O1dat2[i]))
                column_den = np.array([lower_thresh, np.max(O1dat2)])
        # print('number of points 0 column density = ', counts_N)
        # print('number of points 0 rate but non 0 column density = ', counts_Rates)


        #Fitting linear:
        lin_regression = stats.linregress((O1dat_log), (rates_log))
        print('r^2 = ',(lin_regression.rvalue)**2)
        print('slope = ', (lin_regression.slope))
        intercept = lin_regression.intercept
        slope = lin_regression.slope
        
        M_dot = (10**intercept) * column_den**slope
        plt.scatter(np.array(O1dat2), rates2, c = b2, cmap='viridis', label = file)
        plt.plot(column_den, M_dot, 'r-')
        plt.vlines(column_den[0], ymin = np.min(rates2), ymax = np.max(rates2), linestyles = '--', colors = 'Blue')

        plt.xlabel('Column density [cm$^{-2}$]')
        plt.ylabel('Outflow rate [Msun/yr]')
        cbar = plt.colorbar()
        cbar.set_label('impact parameter [kpc]')
        if line == 'H1':
            plt.title('H I, ' + mode + ' outflow'  + ', ' + str(inclin)+'$^{o}$')
            # plt.xlim(1e11, 1e23)
        elif line == 'O1':
            plt.title('O VI, ' + mode + ' outflow'  + ', ' + str(inclin)+'$^{o}$')
            # plt.xlim(column_den[0], column_den[-1])
        elif line == 'C1':
            plt.title('C IV, ' + mode + ' outflow'  + ', ' + str(inclin)+'$^{o}$')
            # plt.xlim(1e11, 5e14)
            # plt.xlim(column_den[0], column_den[-1])
        elif line == 'Mg1':
            plt.title('Mg II, ' + mode + ' outflow'  + ', ' + str(inclin)+'$^{o}$')
            # plt.xlim(1e11, 5e14)
        # plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('Plots/'+mode+line+' column_density.png', dpi = 900)
        plt.show()
        plt.clf()




    species = ['H1', 'H2', 'O1', 'O2', 'C1', 'C2', 'Mg1', 'Mg2']
    for line in species:
        if line == 'H1' or line == 'H2':
            rates = np.loadtxt(file_read+file+rates_file_H, delimiter=',')
        elif line == 'O1' or line == 'O2':
            rates = np.loadtxt(file_read+file+rates_file_O, delimiter=',')
        elif line == 'C1' or line == 'C2':
            rates = np.loadtxt(file_read+file+rates_file_C, delimiter=',')
        elif line == 'Mg1' or line == 'Mg2':
            rates = np.loadtxt(file_read+file+rates_file_Mg, delimiter=',')
        b = to_kpc(np.loadtxt(file_read + file + '_b.csv', delimiter=','))
        O1dat = (np.loadtxt(file_read+file+'_EW_'+line+'_for_rate.csv', delimiter = ','))
        vel_out = np.loadtxt(file_read+file+'_vel_out.csv', delimiter = ',')

        if 'GTT' in file:
            col = 'C0'
        else:
            col = 'C1'
        if '18' in file:
            lw = 1
        elif '9' in file:
            lw = 1.5
        elif '4' in file:
            lw = 2
        # plt.plot(O1dat_s, rates_s, '-', linewidth = lw, label = file)
        # plt.fill_betweenx(rates, O1_up, O1_lo, color = 'C0', alpha = 0.2, interpolate=True, label='_nolegend_')
        plt.scatter(O1dat, rates, c = b, cmap='viridis', label = file)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Equivalend width [A]')
        plt.ylabel('Outflow rate [Msun/yr]')
        if line == 'O1':
            plt.title(file + ' O VI 1032 ' + mode)
        elif line == 'O2':
            plt.title(file + ' O VI 1038 ' + mode)
        elif line == 'C1':
            plt.title(file + ' C IV 1548 ' + mode)
        elif line == 'C2':
            plt.title(file + ' C IV 1551 ' + mode)
        elif line == 'H1':
            plt.title(file + ' H I 1216 ' + mode)
        elif line == 'H2':
            plt.title(file + ' H I 973 ' + mode)
        elif line == 'Mg1':
            plt.title(file + ' Mg II 1240 ' + mode)
        elif line == 'Mg2':
            plt.title(file + ' Mg II 1240.4 ' + mode)
        cbar = plt.colorbar()
        cbar.set_label('impact parameter [kpc]')
        # plt.legend()
        # plt.savefig('Plots/'+mode+line+' Equivalent_width.png', dpi = 900)
        # plt.show()
        plt.clf()


def r_squared(line, mode, file_read, sims):

    if mode == 'sphere':
        rates_file_H = '_N_rates_sph.csv'
        rates_file_O = '_N_rates_sph.csv'
        rates_file_C = '_N_rates_sph.csv'
        rates_file_Mg = '_N_rates_sph.csv'
    elif mode == 'plate':
        rates_file_H = '_N_rates_pla.csv'
        rates_file_O = '_N_rates_pla.csv'
        rates_file_C = '_N_rates_pla.csv'
        rates_file_Mg = '_N_rates_pla.csv'
    elif mode == 'ray':
        rates_file_H = '_N_rates_ray_H.csv'
        rates_file_O = '_N_rates_ray_O.csv'
        rates_file_C = '_N_rates_ray_C.csv'
        rates_file_Mg = '_N_rates_ray_Mg.csv'
    
    if line == 'H1':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_H, delimiter=',')
        lower_lim = 2e11
    elif line == 'O1':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_O, delimiter=',')
        lower_lim = 2e7
    elif line == 'C1':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_C, delimiter=',')
        lower_lim = 4e7
    elif line == 'Mg1':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_Mg, delimiter=',')
        lower_lim = 1e7

    if line == 'H1' and mode == 'sphere':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_H, delimiter=',')
        lower_lim = 5e11
    elif line == 'O1' and mode == 'sphere':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_O, delimiter=',')
        lower_lim = 5e12
    elif line == 'C1' and mode == 'sphere':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_C, delimiter=',')
        lower_lim = 1e12
    elif line == 'Mg1' and mode == 'sphere':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_Mg, delimiter=',')
        lower_lim = 1e7

    if line == 'H1' and mode == 'plate':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_H, delimiter=',')
        lower_lim = 5e11
    elif line == 'O1' and mode == 'plate':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_O, delimiter=',')
        lower_lim = 5e12
    elif line == 'C1' and mode == 'plate':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_C, delimiter=',')
        lower_lim = 1e12
    elif line == 'Mg1' and mode == 'plate':
        rates = np.loadtxt(file_read+'SINK_9pc'+rates_file_Mg, delimiter=',')
        lower_lim = 1e7

    
    O1dat = (np.loadtxt(file_read+'SINK_9pc'+'_N_tri_'+line+'_for_rate.csv', delimiter = ','))
    b = (np.loadtxt(file_read + sims + '_b.csv', delimiter=','))
    O1dat2 = []
    rates2 = []
    b2 = []
    for i in range(len(b)):
        # if b[i] < 50:
        if O1dat[i]>=lower_lim:
            O1dat2.append(O1dat[i])
            rates2.append(rates[i])
            b2.append(b[i])
    O1dat_log = []
    rates_log = []
    for i in range(len(rates2)):
        if O1dat2[i] > 0 and rates2[i]>0:
            rates_log.append(np.log(rates2[i]))
            O1dat_log.append(np.log(O1dat2[i]))
        else:
            pass
    return (stats.linregress((O1dat_log), (rates_log)).rvalue)**2,\
            stats.linregress((O1dat_log), (rates_log)).slope,\
            stats.linregress((O1dat_log), (rates_log)).stderr
files = ['SINK_9pc']

# gen_start(90, 1)

start, end, b = Random_rays.start_end([1, 0.5, 0.5], maxb = 40, minb = 5)
print((b))
print(np.sqrt((start[2]-0.5)**2 + (start[0]-0.5)**2 + (start[1]-0.5)**2))

#%% Real



find_EW_rates(file = files, bmax = 25, temp = None, vel = None, rho = None, QSO = False, ite = '_'+str('real'),\
            _temperature_mu_2 = None,\
                _density_2 = None, _vel = None, mode_alt = 'real', inclination = 'rand')

# incls = np.append(np.arange(0, 92, 4), 90)
# print(incls)
# files = ['SINK_4pc', 'GTT_4pc']
# for i in [0]:
#     print(i)
#     find_EW_rates(file = files, bmax = 25, temp = None, vel = None, rho = None, QSO = False, ite = '_'+str('real'),\
#                 _temperature_mu_2 = None,\
#                     _density_2 = None, _vel = None, mode_alt = 'real', inclination = i)


# files = ['SINK_4pc', 'GTT_4pc', 'GTT_9pc', 'SINK_18pc', 'GTT_18pc']
# find_EW_rates(file = files, bmax = 25, temp = None, vel = None, rho = None, QSO = False, ite = '_'+str('real'),\
#             _temperature_mu_2 = None,\
#                 _density_2 = None, _vel = None, mode_alt = 'real', inclination = 'rand')


#%%
incls = np.append(np.arange(0, 92, 4), 90)
files = ['SINK_9pc']
for i in incls:
    print(i)
    plate(file = files, bmax = 25, temp = None, vel = None, rho = None, QSO = False, ite = '_'+str('real'),\
                _temperature_mu_2 = None,\
                    _density_2 = None, _vel = None, mode_alt = 'real', inclination = i)
    
files = ['GTT_9pc', 'SINK_18pc', 'GTT_18pc']
for i in [0]:
    print(i)
    plate(file = files, bmax = 25, temp = None, vel = None, rho = None, QSO = False, ite = '_'+str('real'),\
                _temperature_mu_2 = None,\
                    _density_2 = None, _vel = None, mode_alt = 'real', inclination = i)










#%% In between
print('uniform Cone') 
_rho = 3e-25
_Temp = 1e4
_vel_r = 500*1000

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

find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = 'between(e25)_with_ball_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination=0)
incs = inclinations
for i in incs:
    print(i)
    find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = 'between(e25)_with_ball_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination=i)


#%% Sphere
print('Inverse square') 
_rho = 8e-25
_Temp = 1e4
_vel_r = 200*1000

def _temperature_mu_2(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc')
    rv = []
    mask = ((r < 42))
    rv = np.where(mask, _Temp, 1e6)
    return yt.YTArray(rv, 'K')
def _density_2(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc').value
    rv = data['gas', 'density'].in_units('g/cm**3').value
    def den(r):
        return _rho/r**2
    mask = ((r < 42) & (r > 0))
    rv = np.where(mask, den(r), 1e-29)
    return yt.YTArray(rv, 'g/cm**3')
def _vel(field, data):
    theta = data['gas', 'theta']
    r = data['gas', 'radius'].in_units('kpc')
    rv = data['gas', 'radial_velocity'].in_units('cm/s').value
    T = data[('gas', 'temperature_mu')]
    mask = ((r < 42))
    rv = np.where(mask, _vel_r, 0)
    return yt.YTArray(rv, 'cm/s')

find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = '_sph_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination='rand')
incs = inclinations
for i in incs:
    print(i)
    find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = '_sph_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination=i)



#%% Inverse square
print('Inverse square') 
_rho = 8e-25
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
    r = data['gas', 'radius'].in_units('kpc').value
    rv = data['gas', 'density'].in_units('g/cm**3').value
    def den(r):
        return _rho/r**2
    mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 41) & (r > 0))
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

find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = '_inverse_sq_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination='rand')
incs = inclinations
for i in incs:
    print(i)
    find_EW_rates(file = files, bmax = 25, temp = _Temp, rho = _rho, vel = _vel_r,\
              QSO = False, ite = '_inverse_sq_'+str(int(_vel_r/1000)) + '_' + str(int(_Temp)),\
              _temperature_mu_2 = _temperature_mu_2,\
                _density_2 = _density_2, _vel = _vel, mode_alt = 'uniform', inclination=i)



#%%

incls = np.append(np.arange(0, 92, 4), 90)

files = ['SINK_9pc', 'GTT_9pc', 'SINK_18pc', 'GTT_18pc', 'SINK_4pc', 'GTT_4pc']

# for file in files:
#     Ew_plot_rates(file, 'sphere', QSO = False, file_read = file+'_Reals/Rates_real_0/', file_rates = file+'_Reals/Rates_real_0/', species = ['C1'])

# Ew_plot_rates('SINK_9pc', 'sphere', QSO = False, file_read = 'Reals/Rates_real_0/', file_rates = 'Reals/Rates_real_0/')

for i in files:
    print(i)
    Ew_plot_rates(i, 'ray', QSO = False, file_read = '_Reals/Rates_real_', inclin = 'rand', species = ['O1'])



# TODO: verify that large b points exist but has 0 N
# TODO: do slope of truncated part against the inclination
# try slope, maybe egde one would be higher

#%%
incls = np.append(np.arange(0, 92, 4), 90)
print(len(incls))
for ele in ['Mg1', 'H1', 'C1', 'O1']:
    sqs = []
    slopes = []
    stderr = []

    for incl in incls:
        # Ew_plot_rates(files, 'ray', QSO = False, file_read = 'Rates_real_'+str(incl)+'/')
        # Ew_plot_rates(files, 'sphere', QSO = False, file_read = 'Rates_real/')
        r_sq, slo, err = r_squared(ele, 'ray', 'SINK_9pc_Reals/Rates_real_'+str(incl)+'/', sims = 'SINK_9pc')
        sqs.append(r_sq)
        slopes.append(slo)
        stderr.append(err)

    plt.errorbar(incls, slopes, yerr = stderr, fmt='.', linewidth = 1, capsize = 1.4, ecolor='black', markersize = 0.1, zorder = 1)
    plt.scatter(incls, slopes, c = sqs, cmap='viridis', zorder = 3, vmin = 0, vmax = 1)
    # plt.colorbar()
    cbar = plt.colorbar()
    cbar.set_label('r^2')
    plt.xlabel('inclination (degrees)')
    plt.ylabel('slope')
    if ele == 'Mg1':
        plt.title('ray Mg II')
    if ele == 'H1':
        plt.title('ray H I')
    if ele == 'O1':
        plt.title('ray O VI')
    if ele == 'C1':
        plt.title('ray C IV')
    plt.show()

for ele in ['Mg1', 'H1', 'C1', 'O1']:
    sqs = []
    slopes = []
    stderr = []

    for incl in incls:
        # Ew_plot_rates(files, 'ray', QSO = False, file_read = 'Rates_real_'+str(incl)+'/')
        # Ew_plot_rates(files, 'sphere', QSO = False, file_read = 'Rates_real/')
        r_sq, slo, err = r_squared(ele, 'sphere', 'SINK_9pc_Reals/Rates_real_'+str(incl)+'/', sims = 'SINK_9pc')
        sqs.append(r_sq)
        slopes.append(slo)
        stderr.append(err)

    plt.errorbar(incls, slopes, yerr = stderr, fmt='.', linewidth = 1, capsize = 1.4, ecolor='black', markersize = 0.1, zorder = 1)
    plt.scatter(incls, slopes, c = sqs, cmap='viridis', zorder = 3, vmin = 0, vmax = 1)
    cbar = plt.colorbar()
    cbar.set_label('r^2')
    plt.xlabel('inclination (degrees)')
    plt.ylabel('slope')
    if ele == 'Mg1':
        plt.title('sphere Mg II')
    if ele == 'H1':
        plt.title('sphere H I')
    if ele == 'O1':
        plt.title('sphere O VI')
    if ele == 'C1':
        plt.title('sphere C IV')
    plt.show()



# for ele in ['Mg1', 'H1', 'C1', 'O1']:
#     sqs = []
#     slopes = []
#     stderr = []

#     for incl in incls:
#         r_sq, slo, err = r_squared(ele, 'plate', 'SINK_9pc_Reals/Rates_real_'+str(incl)+'/', sims = 'SINK_9pc')
#         sqs.append(r_sq)
#         slopes.append(slo)
#         stderr.append(err)

#     plt.errorbar(incls, slopes, yerr = stderr, fmt='.', linewidth = 1, capsize = 1.4, ecolor='black', markersize = 0.1, zorder = 1)
#     plt.scatter(incls, slopes, c = sqs, cmap='viridis', zorder = 3)
#     plt.colorbar()
#     plt.xlabel('inclination (degrees)')
#     plt.ylabel('slope')
#     plt.title('plate ' + ele)
#     plt.show()


# %%
print('Ratesbetween(e25)_with_ball_200_10000_0/')
Ew_plot_rates(files, 'ray', QSO = False, file_read = 'Ratesbetween(e25)_with_ball_200_10000_rand/')
# Ew_plot_rates(files, 'sphere', QSO = False, file_read = 'Ratesbetween(e25)_with_ball_200_10000_45/')

#%%
print('Rates_inverse_sq_200_10000_rand/')
Ew_plot_rates(files, 'ray', QSO = False, file_read = 'Rates_inverse_sq_200_10000_/')
# Ew_plot_rates(files, 'sphere', QSO = False, file_read = 'Rates_inverse_sq_200_10000/')

# %%
inclinations = np.arange(0, 91, 5)
print((inclinations))
Ew_plot_rates(files, 'ray', QSO = False, file_read = 'Rates_real_35/')



#%%

Ns = np.loadtxt('SINK_9pc_Reals/Rates_real_0/SINK_9pc_N_C1_for_rate.csv', delimiter = ',')
start = np.loadtxt('ray_vec copy/start_0.csv', delimiter = ',')
end = np.loadtxt('ray_vec copy/end_0.csv', delimiter = ',')
b = np.loadtxt('ray_vec copy/b_0.csv', delimiter = ',')


ds, info = outflow_rate.load_data(203, file = 'SINK_9pc')
line_list = ['O', 'H', 'C', 'Mg']
trident.add_ion_fields(ds, ions=['Mg II'])
trident.add_ion_fields(ds, ions=['O VI'])
trident.add_ion_fields(ds, ions=['C IV'])
trident.add_ion_fields(ds, ions=['H I'])

start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<2e11]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=2e11]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'H_p0_density'))
p3.annotate_contour(("gas", "H_p0_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('H I, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/H I face on regions.png')
p3.show()

start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<2e7]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=2e7]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'O_p5_density'))
p3.annotate_contour(("gas", "O_p5_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('O VI, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/O VI face on regions.png')
p3.show()




start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<4e7]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=4e7]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'C_p3_density'))
p3.annotate_contour(("gas", "C_p3_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('C IV, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/C IV face on regions.png')
p3.show()


# %%

Ns = np.loadtxt('SINK_9pc_Reals/Rates_real_90/SINK_9pc_N_C1_for_rate.csv', delimiter = ',')
start = np.loadtxt('ray_vec copy/start_90.csv', delimiter = ',')
end = np.loadtxt('ray_vec copy/end_90.csv', delimiter = ',')
b = np.loadtxt('ray_vec copy/b_90.csv', delimiter = ',')
for i in range(1000):
    print(start[i][2]-0.5)
    print(b)

x_rot = [np.sqrt(b[i]*b[i] - (start[i][2]-0.5)*(start[i][2]-0.5)) for i in range(len(b))]
print(x_rot)

#%%


ds, info = outflow_rate.load_data(203, file = 'SINK_9pc')
line_list = ['O', 'H', 'C', 'Mg']
trident.add_ion_fields(ds, ions=['Mg II'])
trident.add_ion_fields(ds, ions=['O VI'])
trident.add_ion_fields(ds, ions=['C IV'])
trident.add_ion_fields(ds, ions=['H I'])

start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<2e11]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=2e11]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'H_p0_density'))
p3.annotate_contour(("gas", "H_p0_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('H I, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/H I face on regions.png')
p3.show()

start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<2e7]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=2e7]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'O_p5_density'))
p3.annotate_contour(("gas", "O_p5_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('O VI, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/O VI face on regions.png')
p3.show()




start_reg_1 = [start[i] for i in range(len(Ns)) if Ns[i]<4e7]
start_reg_2 = [start[i] for i in range(len(Ns)) if Ns[i]>=4e7]
end_reg_1 = [end[i] for i in range(len(Ns)) if Ns[i]>1e11]
print(len(start_reg_1))

p3 = yt.ProjectionPlot(ds, 'z', ('gas', 'C_p3_density'))
p3.annotate_contour(("gas", "C_p3_density"), clim = (0.8e-11, 0.801e-11))
for i in range(len(start_reg_1)):
    # print(i)
    # ray = trident.make_simple_ray(ds,
    #                             start_position=start_reg_1[i],
    #                             end_position=end_reg_1[i],
    #                             data_filename="ray.h5",
    #                             lines = line_list,
    #                             )
    p3.annotate_title('C IV, Face on')
    p3.annotate_marker((start_reg_1[i][0], start_reg_1[i][1], start_reg_1[i][2]), coord_system="data", color = 'turquoise')
for i in range(len(start_reg_2)):
    p3.annotate_marker((start_reg_2[i][0], start_reg_2[i][1], start_reg_2[i][2]), coord_system="data", color = 'white')
p3.zoom(1.7)
p3.save('regions_illus/C IV face on regions.png')
p3.show()