#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import outflow_rate as outflow_rate
import fitting as fitting
import matplotlib.animation as ani
import Random_rays as Random_rays
from datetime import datetime
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
    ds, info = outflow_rate.load_data(time_step, file)
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

def N_find(file):
    rates = []
    bs = []
    nums = np.arange(122, 204, step = 1, dtype = int)
    nums = [203]

    temp_O1 = []
    temp_C1 = []
    temp_Mg1 = []
    temp_H1 = []
    temp_rates = []
    temp_b = []
    params_log = (np.random.rand(30))*1.75
    params = np.sort(10**params_log)
    ds, rate = find_rate(203, file)
    for i in range(40):
        print(i)
        k = np.random.rand(3)
        # print(k)
        temp_temp_O1 = []
        temp_temp_C1 = []
        temp_temp_Mg1 = []
        temp_temp_H1 = []
        for j in params:
            start, end, b = Random_rays.start_end(k, to_cl(j))
            b = j
            #print(start, end)
            # print('a')
            N = ew(ds, start, end,
                    ['O VI 1032', 'Mg II 1240',\
                    'C IV 1548', 'H I 1216'],
                    plot = False, crit = 'chi', MilkyWay = False, QSO = False, Noise = False)
            # print('b')
            temp_temp_O1.append(N['O VI 1032'])
            temp_temp_Mg1.append(N['Mg II 1240'])
            temp_temp_C1.append(N['C IV 1548'])
            temp_temp_H1.append(N['H I 1216'])
            # print('c')
        temp_O1.append(temp_temp_O1)
        temp_Mg1.append(temp_temp_Mg1)
        temp_C1.append(temp_temp_C1)
        temp_H1.append(temp_temp_H1)
        temp_b = params

        bs.append(np.average(temp_b))
    np.savetxt('C1_'+file+'.csv', temp_C1, delimiter=',')
    np.savetxt('O1_'+file+'.csv', temp_O1, delimiter=',')
    np.savetxt('H1_'+file+'.csv', temp_H1, delimiter=',')
    np.savetxt('Mg1_'+file+'.csv', temp_Mg1, delimiter=',')
    np.savetxt('b_'+file+'.csv', temp_b, delimiter=',')
    print('FIN.')
    return 0
#%%

def N_plot(file, name):
    def med_upper_lower(a):
            median = []
            upper = []
            lower = []
            for i in range (len(a)):
                med = np.median(a[i])
                low_per_ind = int(np.round(len(a[i])*0.159))
                high_per_ind = int(np.round(len(a[i])*0.841))
                median.append(med)
                upper.append(np.sort(a[i])[high_per_ind])
                lower.append(np.sort(a[i])[low_per_ind])
                #print(np.min(a[i]))
            return (median, upper, lower)
    # first you have to make the figure
    fig, axs = plt.subplots(2, 2, sharex='col')

    # O
    for sims in file:
        b = np.loadtxt('b_'+sims+'.csv', delimiter=',')
        O1dat = np.transpose(np.loadtxt('O1_'+sims+'.csv', delimiter = ','))
        O1_med, O1_up, O1_lo = med_upper_lower(O1dat)
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

        axs[1,1].plot(b, O1_med, col, '-', linewidth = lw, label = sims)
        axs[1,1].fill_between(b, O1_up, O1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[1,1].set_yscale('log')
    axs[1,1].yaxis.set_label_position("right")
    axs[1,1].set_ylim(1e8, 1e15)
    axs[1,1].set(xlabel='impact parameter [kpc]')
    #axs[1,1].set(ylabel='Columns density ' + "[cm$^{-2}$]")
    axs[1,1].yaxis.tick_right()

    #     plt.plot(b, O1_med, col, '-', linewidth = lw, label = sims)
    #     plt.fill_between(b, O1_up, O1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    # # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylabel('Columns density ' + "[cm$^{-2}$]")
    # plt.xlabel('impact parameter [kpc]')
    # plt.title('O VI')
    # #plt.legend()
    # plt.savefig(name+' O VI Column_b.pdf', dpi = 900)
    # plt.show()
    # plt.clf()

    # C
    for sims in file:
        b = np.loadtxt('b_'+sims+'.csv', delimiter=',')
        C1dat = np.transpose(np.loadtxt('C1_'+sims+'.csv', delimiter = ','))
        C1_med, C1_up, C1_lo = med_upper_lower(C1dat)
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

        axs[1,0].plot(b, C1_med, col, '-', linewidth = lw, label = sims)
        axs[1,0].fill_between(b, C1_up, C1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylim(1e8, 1e15)
    axs[1,0].set(xlabel='impact parameter [kpc]', ylabel='Columns density ' + "[cm$^{-2}$]")

    #     plt.plot(b, C1_med, col, '-', linewidth = lw, label = sims)
    #     plt.fill_between(b, C1_up, C1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    # # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylabel('Columns density ' + "[cm$^{-2}$]")
    # plt.xlabel('impact parameter [kpc]')
    # plt.title('C IV')
    # plt.legend()
    # plt.savefig(name+' C IV Column_b.pdf', dpi = 900)
    # plt.show()
    # plt.clf()

    # H
    for sims in file:
        b = np.loadtxt('b_'+sims+'.csv', delimiter=',')
        H1dat = np.transpose(np.loadtxt('H1_'+sims+'.csv', delimiter = ','))
        H1_med, H1_up, H1_lo = med_upper_lower(H1dat)
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
        axs[0,0].plot(b, H1_med, col, '-', linewidth = lw, label = sims)
        axs[0,0].fill_between(b, H1_up, H1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[0,0].set_yscale('log')
    axs[0,0].set_ylim(1e11, 1e19)
    axs[0,0].set(ylabel='Columns density ' + "[cm$^{-2}$]")
    # axs[1,0].set_title('C VI')



    #     plt.plot(b, H1_med, col, '-', linewidth = lw, label = sims)
    #     plt.fill_between(b, H1_up, H1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    # # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylabel('Columns density ' + "[cm$^{-2}$]")
    # plt.xlabel('impact parameter [kpc]')
    # plt.title('H I')
    # plt.legend()
    # plt.savefig(name + ' H I Column_b.pdf', dpi = 900)
    # plt.show()
    # plt.clf()

    # Mg
    for sims in file:
        b = np.loadtxt('b_'+sims+'.csv', delimiter=',')
        Mg1dat = np.transpose(np.loadtxt('Mg1_'+sims+'.csv', delimiter = ','))
        Mg1_med, Mg1_up, Mg1_lo = med_upper_lower(Mg1dat)
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
        axs[0,1].plot(b, Mg1_med, col, '-', linewidth = lw, label = sims)
        axs[0,1].fill_between(b, Mg1_up, Mg1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[0,1].set_yscale('log')
    axs[0,1].set_ylim(2e6, 2e14)
    axs[0,1].yaxis.tick_right()
    axs[0,1].yaxis.set_label_position("right")
    #axs[0,1].set(ylabel='Columns density ' + "[cm$^{-2}$]")


    axs[1,1].set_xlim(4.11, 41.1)
    axs[1,0].set_xlim(4.11, 41.1)
    axs[0,1].set_xlim(4.11, 41.1)
    axs[0,0].set_xlim(4.11, 41.1)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('subplots.pdf', dpi = 900)
    plt.show()



    #     plt.plot(b, Mg1_med, col, '-', linewidth = lw, label = sims)
    #     print('aaaaa')
    #     plt.fill_between(b, Mg1_up, Mg1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    # # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylabel('Columns density ' + "[cm$^{-2}$]")
    # plt.xlabel('impact parameter [kpc]')
    # plt.title('Mg II')
    # plt.legend()
    # plt.savefig(name + ' Mg II Column_b.pdf', dpi = 900)
    # plt.show()
    return 0
# %%

files = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'SINK_9pc', 'SINK_18pc']
# for i in files:
#     print(i)
#     N_find(i)


# %%
# files = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc']
# N_plot(files, 'GTT')
# %%
# files = ['SINK_4pc', 'SINK_9pc', 'SINK_18pc']
# N_plot(files, 'SINK')
# %%
files = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'SINK_9pc', 'SINK_18pc']
N_plot(files, 'All')

