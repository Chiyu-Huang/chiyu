#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import Codes.outflow_rate as outflow_rate
import Codes.fitting as fitting
import matplotlib.animation as ani
import Codes.Random_rays as Random_rays
from datetime import datetime
from Codes.Ray_maker import find_rate, to_cl, ew, to_kpc
#%%
yt.set_log_level(50)
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
    N_temp_O1 = []
    N_temp_C1 = []
    N_temp_Mg1 = []
    N_temp_H1 = []
    temp_rates = []
    temp_b = []
    params_log = (np.random.rand(15))*1.75
    params = np.sort(10**params_log)
    ds, rate = find_rate(203, file)
    for i in range(20):
        print(i)
        k = np.random.rand(3)-0.5
        # print(k)
        temp_temp_O1 = []
        temp_temp_C1 = []
        temp_temp_Mg1 = []
        temp_temp_H1 = []
        N_temp_temp_O1 = []
        N_temp_temp_C1 = []
        N_temp_temp_Mg1 = []
        N_temp_temp_H1 = []
        count = 0
        for j in params:
            print(count)
            count+=1
            start, end, b = Random_rays.start_end(k, to_cl(j))
            b = j
            # print(start, end)
            start_t=datetime.now()
            print('b = ', b)
            # print('a')
            # print('theoretically EW?')
            ew_fit, ew_raw, Ns, Ns_tri = ew(ds, start, end,
                        ['O VI 1032', 'Mg II 1240',\
                        'C IV 1548', 'H I 1216',\
                        'O VI 1038', 'Mg II 1240.4',\
                        'C IV 1551', 'H I 973'],\
                        plot = False, crit = 'chi',\
                        MilkyWay = False, QSO = True, Noise = False, column = False)
            # print(ew_raw['O VI 1032'] + ew_raw['O VI 1038'])
            # print('theoretically N??')
            # ew_fit, ew_raw, Ns, Ns_tri = ew(ds, [0.41204907, 0.29913888, 0.        ], [0.61337747, 0.69514484, 1.        ],
            #             ['O VI 1032',
            #             'O VI 1038'],\
            #             plot = False, crit = 'chi',\
            #             MilkyWay = False, QSO = True, Noise = False, column = False)
            # print('b')
            # print(Ns)
            temp_temp_O1.append(Ns['O VI 1032'] + Ns['O VI 1038'])
            temp_temp_Mg1.append(Ns['Mg II 1240'] + Ns['Mg II 1240.4'])
            temp_temp_C1.append(Ns['C IV 1548'] + Ns['C IV 1551'])
            temp_temp_H1.append(Ns['H I 1216'] + Ns['H I 973'])
            N_temp_temp_O1.append(Ns_tri['O VI 1032'])
            N_temp_temp_Mg1.append(Ns_tri['Mg II 1240'])
            N_temp_temp_C1.append(Ns_tri['C IV 1548'])
            N_temp_temp_H1.append(Ns_tri['H I 1216'])
            print ('runtime = ', datetime.now()-start_t)
        temp_O1.append(temp_temp_O1)
        temp_Mg1.append(temp_temp_Mg1)
        temp_C1.append(temp_temp_C1)
        temp_H1.append(temp_temp_H1)
        N_temp_O1.append(N_temp_temp_O1)
        N_temp_Mg1.append(N_temp_temp_Mg1)
        N_temp_C1.append(N_temp_temp_C1)
        N_temp_H1.append(N_temp_temp_H1)
        temp_b = params

        bs.append(np.average(temp_b))
    np.savetxt('C1_'+file+'.csv', temp_C1, delimiter=',')
    np.savetxt('O1_'+file+'.csv', temp_O1, delimiter=',')
    np.savetxt('Mg1_'+file+'.csv', temp_Mg1, delimiter=',')
    np.savetxt('H1_'+file+'.csv', temp_H1, delimiter=',')
    np.savetxt('N_C1_'+file+'.csv', N_temp_C1, delimiter=',')
    np.savetxt('N_O1_'+file+'.csv', N_temp_O1, delimiter=',')
    np.savetxt('N_Mg1_'+file+'.csv', N_temp_Mg1, delimiter=',')
    np.savetxt('N_H1_'+file+'.csv', N_temp_H1, delimiter=',')
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
        N_O1dat = np.transpose(np.loadtxt('N_O1_'+sims+'.csv', delimiter = ','))
        N_O1_med, N_O1_up, N_O1_lo = med_upper_lower(N_O1dat)
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

        # axs[1,1].plot(b, N_O1dat, 'C0', '.-', label = sims)
        # axs[1,1].plot(b, O1dat, 'C1', 'x-', label = sims)
        axs[1,1].plot(b, O1_med, 'C0', '-', linewidth = lw, label = sims)
        axs[1,1].fill_between(b, O1_up, O1_lo, color = 'C0', alpha = 0.2, interpolate=True, label='_nolegend_')
        axs[1,1].plot(b, N_O1_med, 'C1', '-', linewidth = lw, label = sims)
        axs[1,1].fill_between(b, N_O1_up, N_O1_lo, color = 'C1', alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[1,1].set_yscale('log')
    axs[1,1].yaxis.set_label_position("right")
    # axs[1,1].set_ylim(1e8, 1e15)
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
        N_C1dat = np.transpose(np.loadtxt('N_C1_'+sims+'.csv', delimiter = ','))
        C1_med, C1_up, C1_lo = med_upper_lower(C1dat)
        N_C1_med, N_C1_up, N_C1_lo = med_upper_lower(N_C1dat)
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

        axs[1,0].plot(b, C1_med, 'C0', '-', linewidth = lw, label = sims)
        axs[1,0].fill_between(b, C1_up, C1_lo, color = 'C0', alpha = 0.2, interpolate=True, label='_nolegend_')
        axs[1,0].plot(b, N_C1_med, 'C1', '-', linewidth = lw, label = sims)
        axs[1,0].fill_between(b, N_C1_up, N_C1_lo, color = 'C1', alpha = 0.2, interpolate=True, label='_nolegend_')
        # axs[1,0].plot(b, C1_med, col, '-', linewidth = lw, label = sims)
        # axs[1,0].fill_between(b, C1_up, C1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[1,0].set_yscale('log')
    # axs[1,0].set_ylim(1e8, 1e15)
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
        N_H1dat = np.transpose(np.loadtxt('N_H1_'+sims+'.csv', delimiter = ','))
        H1dat = np.transpose(np.loadtxt('H1_'+sims+'.csv', delimiter = ','))
        H1_med, H1_up, H1_lo = med_upper_lower(H1dat)
        N_H1_med, N_H1_up, N_H1_lo = med_upper_lower(N_H1dat)
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
        axs[0,0].plot(b, H1_med, 'C0', '-', linewidth = lw, label = sims)
        axs[0,0].fill_between(b, H1_up, H1_lo, color = 'C0', alpha = 0.2, interpolate=True, label='_nolegend_')
        axs[0,0].plot(b, N_H1_med, 'C1', '-', linewidth = lw, label = sims)
        axs[0,0].fill_between(b, N_H1_up, N_H1_lo, color = 'C1', alpha = 0.2, interpolate=True, label='_nolegend_')
        # axs[0,0].plot(b, H1_med, col, '-', linewidth = lw, label = sims)
        # axs[0,0].fill_between(b, H1_up, H1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[0,0].set_yscale('log')
    # axs[0,0].set_ylim(1e11, 1e19)
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
        N_Mg1dat = np.transpose(np.loadtxt('N_Mg1_'+sims+'.csv', delimiter = ','))
        Mg1dat = np.transpose(np.loadtxt('Mg1_'+sims+'.csv', delimiter = ','))
        Mg1_med, Mg1_up, Mg1_lo = med_upper_lower(Mg1dat)
        N_Mg1_med, N_Mg1_up, N_Mg1_lo = med_upper_lower(N_Mg1dat)
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
        axs[0,1].plot(b, Mg1_med, 'C0', '-', linewidth = lw, label = sims)
        axs[0,1].fill_between(b, Mg1_up, Mg1_lo, color = 'C0', alpha = 0.2, interpolate=True, label='_nolegend_')
        axs[0,1].plot(b, N_Mg1_med, 'C1', '-', linewidth = lw, label = sims)
        axs[0,1].fill_between(b, N_Mg1_up, N_Mg1_lo, color = 'C1', alpha = 0.2, interpolate=True, label='_nolegend_')
        # axs[0,1].plot(b, Mg1_med, col, '-', linewidth = lw, label = sims)
        # axs[0,1].fill_between(b, Mg1_up, Mg1_lo, color = col, alpha = 0.2, interpolate=True, label='_nolegend_')
    axs[0,1].set_yscale('log')
    # axs[0,1].set_ylim(2e6, 2e14)
    axs[0,1].yaxis.tick_right()
    axs[0,1].yaxis.set_label_position("right")
    #axs[0,1].set(ylabel='Columns density ' + "[cm$^{-2}$]")


    # axs[1,1].set_xlim(4.11, 41.1)
    # axs[1,0].set_xlim(4.11, 41.1)
    # axs[0,1].set_xlim(4.11, 41.1)
    # axs[0,0].set_xlim(4.11, 41.1)
    axs[1,0].set_xscale('log')
    axs[1,1].set_xscale('log')
    axs[0,1].set_xscale('log')
    axs[0,0].set_xscale('log')
    # plt.subplots_adjust(wspace=0, hspace=0)
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

# files = ['SINK_9pc', 'SINK_18pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'GTT_4pc']
files = ['SINK_9pc']
for i in files:
    print(i)
    N_find(i)


# %%
files = ['SINK_9pc']
N_plot(files, 'SINK')