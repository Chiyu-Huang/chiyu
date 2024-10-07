#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import outflow_rate as outflow_rate
import fitting as fitting
import matplotlib.animation as ani
import Random_rays as Random_rays
#%%
yt.set_log_level(50)

def ew(ds, start, end, lines, crit, plot = False, Noise = True, QSO = True, MilkyWay = True):
    '''
    Function that finds the outflow rate and the equivalent width given a ray

    Args:
        time_step (int) : Time step of the simulation
        start (1D array) : Coordinate of the start in code units (0.5 being the centre)
        end (1D array) : Coordinate of the end similarly
        lines (tuple) : The specific lines you want to fit
        plot (Bool, optional) : Whether or not to plot the ray trajectory (default False)
        Noise (Bool, optional) : Whether or not to add Gaussian Noise (SNR = 20, but can change) (default True)
        QSO (Bool, optional) : Whether or not to add a QSO continuum (default True)
        MilkyWay (Bool, optional) : Whether or not to add Milky Way foreground spectrum (default True)

    Returns:
        float : Outflow rate in Msun/yr
        dictionary : Fitted equivalent width of all populations specified
        dictionary : Raw data equivalent widths
    '''

    #   Creating the light ray
    ray_start = start
    ray_end = end
    ray = trident.make_simple_ray(ds,
                                start_position=ray_start,
                                end_position=ray_end,
                                data_filename="ray.h5",
                                lines = ['H', 'C', 'N', 'O', 'Mg', 'Si'])
    if plot:
        p = yt.ProjectionPlot(ds, 'y', ('gas', 'O_p5_number_density'))
        p.annotate_ray(ray, arrow=True)
        p.show()
        p = yt.ProjectionPlot(ds, 'y', ('gas', 'C_p3_number_density'))
        p.annotate_ray(ray, arrow=True)
        p.show()
        p = yt.ProjectionPlot(ds, 'y', ('gas', 'Mg_p1_number_density'))
        p.annotate_ray(ray, arrow=True)
        p.show()
    
    wavelength = []
    flux = []
    # plt.plot(wavelength, flux)
    # plt.show()

    #   specifying background then fit
    background = np.ones(len(flux)) #need to change later
    #ew_fit, ew_raw, N = fitting.fitting(wavelength, flux, lines, bg=background, ray=ray, crit=crit, plot_fit = False)
    N = fitting.fitting(wavelength, flux, lines, bg=background, ray=ray, crit=crit, plot_fit = False, column = True)
    return (N)

def find_rate(time_step, file):
    #   Determine outflow rate
    ds, info = outflow_rate.load_data(timestep = time_step, file = file)
    rate = outflow_rate.outflow_rate(ds, 10, -10)


    '''Add ion fields, creating the light ray, and finding the equivalent width'''
    #   Adding ion fields
    trident.add_ion_fields(ds, ions=['Mg II'])
    trident.add_ion_fields(ds, ions=['O VI'])
    trident.add_ion_fields(ds, ions=['C IV'])
    trident.add_ion_fields(ds, ions=['H I'])
    return ds, rate

def to_cl(kpc):
    return kpc*3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22) #have to change accordingly with the specific file

def find_N_rates (file):
    rates = []
    bs = []
    nums = np.arange(122, 204, step = 1, dtype = int)
    for sims in file:
        print(sims)
        N_O1 = []
        N_C1 = []
        N_Mg1 = []
        N_H1 = []
        rates = []
        for ts in nums:
            print(ts)
            ds, rate = find_rate(ts, file = sims)
            temp_O1 = []
            temp_C1 = []
            temp_Mg1 = []
            temp_H1 = []
            for i in range(101):
                print(i, '/100')
                k = np.random.rand(3)
                # print(k)
                start, end, b = Random_rays.start_end(k)
                #print(start, end)
                N = ew(ds, start, end,
                        ['O VI 1032', 'Mg II 1240',\
                        'C IV 1548', 'H I 1216'],
                        plot = False, crit = 'chi', MilkyWay = False, QSO = False, Noise = False)
                # ew_fit, ew_raw = ew(ds, start, end,
                #                             ['H I 1216', 'H I 973'],\
                #                                 plot = False, crit = 'chi', MilkyWay = False, QSO = False, Noise = False)
                temp_O1.append(N['O VI 1032'])
                temp_Mg1.append(N['Mg II 1240'])
                temp_C1.append(N['C IV 1548'])
                temp_H1.append(N['H I 1216'])
            N_O1.append(temp_O1)
            N_Mg1.append(temp_Mg1)
            N_C1.append(temp_C1)
            N_H1.append(temp_H1)
            rates.append(rate)
        np.savetxt('CSVs_column_rates/'+sims+'C1_for_rate.csv', N_C1, delimiter=',')
        np.savetxt('CSVs_column_rates/'+sims+'O1_for_rate.csv', N_O1, delimiter=',')
        np.savetxt('CSVs_column_rates/'+sims+'H1_for_rate.csv', N_H1, delimiter=',')
        np.savetxt('CSVs_column_rates/'+sims+'Mg1_for_rate.csv', N_Mg1, delimiter=',')
        np.savetxt('CSVs_column_rates/'+sims+'rates.csv', rates, delimiter=',')
        print('FIN.')
    return 0

#%%
# find_N_rates(file = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'SINK_9pc', 'SINK_18pc'])

# %%
def N_plot_rates(file, name):

    def sort_by_rates(rates, dat):
        Big_list = np.array([rates, dat])
        
        return(Big_list[:, Big_list[0].argsort()][1])
    
    def med_upper_lower(a):
            print(a)
            print(len(a))
            print(len(a[1]))
            median = []
            upper = []
            lower = []
            cal = []
            for i in range (len(a)):
                med = np.median(a[i])
                low_per_ind = int(np.round(len(a[i])*0.159))
                high_per_ind = int(np.round(len(a[i])*0.841))
                median.append(med)
                upper.append(np.sort(a[i])[high_per_ind])
                lower.append(np.sort(a[i])[low_per_ind])
                cal.append(a[i][6])

            print(median)
            return (median, upper, lower, cal)
    colm = []
    # O VI
    for sims in file:
        rates = np.loadtxt('CSVs_column_rates/'+sims+'rates.csv', delimiter=',')
        O1dat = (np.loadtxt('CSVs_column_rates/'+sims+'O1_for_rate.csv', delimiter = ','))
        O1_med, O1_up, O1_lo, cal = med_upper_lower(O1dat)
        colm.append(O1_med)
        plt.plot(cal, rates, '.')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
        O1_med = sort_by_rates(rates, O1_med)
        O1_up = sort_by_rates(rates, O1_up)
        O1_lo = sort_by_rates(rates, O1_lo)
        rates = np.sort(rates)
        if 'GTT' in sims:
            col = 'C0'
        else:
            col = 'C1'
        if '18' in sims:
            lw = 1
        elif '9' in sims:
            lw = 1.5
        elif '4' in sims:
            lw = 2
        plt.plot(O1_med, rates, col, '.', linewidth = lw, label = sims)
        # plt.plot(cal, rates, col, '.', linewidth = lw, label = sims)
        plt.fill_betweenx(rates, O1_up, O1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Columns density ' + "[cm$^{-2}$]")
    plt.ylabel('Outflow rate [Msun/yr]')
    plt.title('O VI')
    # plt.legend()
    plt.savefig(name+' O VI Column_rates.png', dpi = 900)
    plt.show()
    plt.clf()

    # C IV
    
    for sims in file:
        rates = np.loadtxt('CSVs_column_rates/'+sims+'rates.csv', delimiter=',')
        C1dat = (np.loadtxt('CSVs_column_rates/'+sims+'C1_for_rate.csv', delimiter = ','))
        C1_med, C1_up, C1_lo, cal = med_upper_lower(C1dat)
        colm.append(C1_med)
        plt.plot(cal, rates, '.')
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
        C1_med = sort_by_rates(rates, C1_med)
        C1_up = sort_by_rates(rates, C1_up)
        C1_lo = sort_by_rates(rates, C1_lo)
        cal = sort_by_rates(rates, cal)
        rates = np.sort(rates)
        if 'GTT' in sims:
            col = 'C0'
        else:
            col = 'C1'
        if '18' in sims:
            lw = 1
        elif '9' in sims:
            lw = 1.5
        elif '4' in sims:
            lw = 2
        plt.plot(C1_med, rates, col, '-', linewidth = lw, label = sims)
        # plt.plot(cal, rates, col, marker = '.', linewidth = lw, label = sims)
        plt.fill_betweenx(rates, C1_up, C1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Columns density ' + "[cm$^{-2}$]")
    plt.ylabel('Outflow rate [Msun/yr]')
    plt.title('C IV')
    # plt.legend()
    plt.savefig(name+' C IV Column_rates.png', dpi = 900)
    plt.show()
    plt.clf()

    # H I
    for sims in file:
        rates = np.loadtxt('CSVs_column_rates/'+sims+'rates.csv', delimiter=',')
        H1dat = (np.loadtxt('CSVs_column_rates/'+sims+'H1_for_rate.csv', delimiter = ','))
        H1_med, H1_up, H1_lo, cal = med_upper_lower(H1dat)
        colm.append(H1_med)
        H1_med = sort_by_rates(rates, H1_med)
        H1_up = sort_by_rates(rates, H1_up)
        H1_lo = sort_by_rates(rates, H1_lo)
        rates = np.sort(rates)
        if 'GTT' in sims:
            col = 'C0'
        else:
            col = 'C1'
        if '18' in sims:
            lw = 1
        elif '9' in sims:
            lw = 1.5
        elif '4' in sims:
            lw = 2
        plt.plot(H1_med, rates, col, '-', linewidth = lw, label = sims)
        plt.fill_betweenx(rates, H1_up, H1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Columns density ' + "[cm$^{-2}$]")
    plt.ylabel('Outflow rate [Msun/yr]')
    plt.title('H I')
    # plt.legend()
    plt.savefig(name + ' H I Column_rates.png', dpi = 900)
    plt.show()
    plt.clf()

    # Mg II
    for sims in file:
        rates = np.loadtxt('CSVs_column_rates/'+sims+'rates.csv', delimiter=',')
        Mg1dat = (np.loadtxt('CSVs_column_rates/'+sims+'Mg1_for_rate.csv', delimiter = ','))
        Mg1_med, Mg1_up, Mg1_lo, cal = med_upper_lower(Mg1dat)
        colm.append(Mg1_med)
        Mg1_med = sort_by_rates(rates, Mg1_med)
        Mg1_up = sort_by_rates(rates, Mg1_up)
        Mg1_lo = sort_by_rates(rates, Mg1_lo)
        rates = np.sort(rates)
        if 'GTT' in sims:
            col = 'C0'
        else:
            col = 'C1'
        if '18' in sims:
            lw = 1
        elif '9' in sims:
            lw = 1.5
        elif '4' in sims:
            lw = 2
        plt.plot(Mg1_med, rates, col, '-', linewidth = lw, label = sims)
        print('aaaaa')
        plt.fill_betweenx(rates, Mg1_up, Mg1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Columns density ' + "[cm$^{-2}$]")
    plt.ylabel('Outflow rate [Msun/yr]')
    plt.title('Mg II')
    # plt.legend()
    plt.savefig(name + ' Mg II Column_rates.png', dpi = 900)
    plt.show()

    rates = np.loadtxt('CSVs_column_rates/'+'GTT_9pc'+'rates.csv', delimiter=',')
    nums = np.arange(122, 204, step = 1, dtype = int)
    plt.plot(nums, rates, '.')
    plt.show()
    plt.plot(nums, colm[0], '.')
    plt.show()


    return 0

files = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'SINK_9pc', 'SINK_18pc']
files = ['GTT_9pc']
find_N_rates(file = files)

# %%
# files = ['GTT_9pc']
# # find_N_rates(files)
# N_plot_rates(files, 'GTT')
# # %%
