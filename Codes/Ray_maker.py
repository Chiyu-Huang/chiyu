#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import matplotlib.animation as ani
import Codes.outflow_rate as outflow_rate
import Codes.fitting as fitting
import Codes.Random_rays as Random_rays
yt.set_log_level(50)
from datetime import datetime
#%% Defining functions for EW and rate
def ew(ds, start, end, lines, crit, b, plot = False, Noise = True, QSO = True,
       MilkyWay = True, plot_fit = False, _density_2 = None, _temperature_mu_2 = None,
       _vel = None):
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

    line_list = ['O', 'H', 'C', 'Mg']
    trident.add_ion_fields(ds, ions=['Mg II'])
    trident.add_ion_fields(ds, ions=['O VI'])
    trident.add_ion_fields(ds, ions=['C IV'])
    trident.add_ion_fields(ds, ions=['H I'])
    # line_list = ['H']
    # line_list = 'all'
    
    # print('start spec')
    ray_start = start
    ray_end = end
    start_t = datetime.now()
    ray = trident.make_simple_ray(ds,
                                start_position=ray_start,
                                end_position=ray_end,
                                data_filename="ray.h5",
                                lines = line_list)
    
    if plot:
        p = yt.SlicePlot(ds, 'y', ('gas', 'O_p5_number_density'))
        p.annotate_ray(ray, arrow=True)
        p.zoom(1)
        p.show()
        p = yt.SlicePlot(ds, 'x', ('gas', 'O_p5_number_density'))
        p.annotate_ray(ray, arrow=True)
        p.show()
        p = yt.ProjectionPlot(ds, 'z', ('gas', 'H_p0_number_density'))
        p.zoom(5)
        p.annotate_ray(ray, arrow=True)
        p.show()
    #   Generating spectrum
    
    sg = trident.SpectrumGenerator(lambda_min = 800, lambda_max = 1600, dlambda=50/10000) #Change from_central as change dlambda
    sg.make_spectrum(ray, lines = line_list)
    if QSO == True:
        sg.add_qso_spectrum()
    if Noise == True:
        sg.add_gaussian_noise(20)
    if MilkyWay == True:
        sg.add_milky_way_foreground()
    
    wavelength = np.array(sg.lambda_field.value)
    flux = np.array(sg.flux_field)
    #   specifying background then fit
    # start_t = datetime.now()

    background = np.ones(len(flux)) #TODO: need to change later
    ew_fit, ew_raw, Ns, Ns_tri = fitting.fitting(wavelength, flux, lines,\
                                        bg=background, ray=ray, crit=crit,\
                                        b = b, plot_fit = plot_fit, column = False)
    # print('runtime = ', datetime.now() - start_t)
    return (ew_fit, ew_raw, Ns, Ns_tri, ray)

def find_rate(time_step, file, h):
    #   Determine outflow rate
    ds, info = outflow_rate.load_data(timestep = time_step, file = file)
    rate = outflow_rate.outflow_rate(ds, h, -h)


    '''Add ion fields, creating the light ray, and finding the equivalent width'''
    #   Adding ion field
    return ds, rate

def to_cl(kpc):
    return kpc*3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22) #have to change accordingly with the specific file

def to_kpc(cl):
    return cl/(3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22)) #have to change accordingly with the specific file


#%% Iterating and loading widths and rates
def plot_ani():
    rates = []
    ews_fits_O1 = []
    ews_fits_O2 = []
    ews_raws_O1 = []
    ews_raws_O2 = []
    ews_fits_Mg1 = []
    ews_fits_Mg2 = []
    ews_raws_Mg1 = []
    ews_raws_Mg2 = []
    ews_fits_C1 = []
    ews_fits_C2 = []
    ews_raws_C1 = []
    ews_raws_C2 = []
    nums = np.arange(122, 204, step = 1, dtype = int)
    # nums = [141]
    nums = [1,5,9,13,25,29,41,53,73,89,100,105, 108,112, 120,135,143,150,165,173,185,197,203]
    for i in nums:
        print('timestep = ', i)
        ds, rate = find_rate(i, file = 'GTT_9pc')
        # ew_fit, ew_raw = ew(ds, [1, 0.5, 1], [0, 0.5, 0.1],
                                    # ['O VI 1032', 'O VI 1038', 'Mg II 1240',\
                                    #     'Mg II 1240.4', 'C IV 1548', 'C IV 1551'],
                                    #     plot = False, crit = 'chi', MilkyWay = False, QSO = False, Noise = False)
        ew_fit, ew_raw = ew(ds, [1, 0.5, 1], [0, 0.5, 0.1],
                                       [\
                                        'O VI 1032', 'O VI 1038',\
                                        'Mg II 1240', 'Mg II 1240.4',\
                                        'C IV 1548', 'C IV 1551'],
                                        plot = False, crit = 'chi', MilkyWay = False, QSO = False, Noise = False)
        rates.append(rate)
        ews_fits_O1.append(ew_fit['O VI 1032'])
        ews_fits_O2.append(ew_fit['O VI 1038'])
        ews_raws_O1.append(ew_raw['O VI 1032'])
        ews_raws_O2.append(ew_raw['O VI 1038'])
        ews_fits_Mg1.append(ew_fit['Mg II 1240'])
        ews_fits_Mg2.append(ew_fit['Mg II 1240.4'])
        ews_raws_Mg1.append(ew_raw['Mg II 1240'])
        ews_raws_Mg2.append(ew_raw['Mg II 1240.4'])
        ews_fits_C1.append(ew_fit['C IV 1548'])
        ews_fits_C2.append(ew_fit['C IV 1551'])
        ews_raws_C1.append(ew_raw['C IV 1548'])
        ews_raws_C2.append(ew_raw['C IV 1551'])
    print('FIN.')



    ########################
    ########################
    ########################
    ####        To change max allowed error for a fit to be used


    diag_x = np.arange(0, 0.175, 0.005)
    plt.plot(diag_x, diag_x, '-')
    plt.plot(ews_fits_C1, ews_raws_C1, '.')
    plt.plot(ews_fits_C2, ews_raws_C2, '.')
    plt.plot(ews_fits_Mg1, ews_raws_Mg1, '.')
    plt.plot(ews_fits_Mg2, ews_raws_Mg2, '.')
    plt.plot(ews_fits_O1, ews_raws_O1, '.')
    plt.plot(ews_fits_O2, ews_raws_O2, '.')
    plt.xlabel('W from fitting (A)')
    plt.ylabel('W from data (A)')
    plt.legend(('line of equal width', 'C IV 1548', 'C IV 1551', 'Mg II 1240', 'Mg II 1240.6', 'O VI 1032', 'O VI 1038'))
    plt.savefig('fit_vs_raw_funky_EW_15_pure.PNG')
    plt.show()
    ##%% Plotting the rates with width
    width = []
    for i in range(len(ews_fits_C1)):
        width.append(ews_fits_C1[i]+ews_fits_C2[i]+\
                    ews_fits_Mg1[i]+ews_fits_Mg2[i]+\
                    ews_fits_O1[i]+ews_fits_O2[i])
    plt.plot(width, rates, '.')
    plt.xlabel('equivalent width (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.savefig('rate_ew_funky_EW_15_pure.PNG')
    plt.show()


    Wo1 = []
    Wo2 = []
    Wc1 = []
    Wc2 = []
    Wm1 = []
    Wm2 = []
    for i in range(len(ews_fits_C1)):
        Wc1.append(ews_raws_C1[i])
        Wo1.append(ews_raws_O1[i])
        Wm1.append(ews_raws_Mg1[i])
        Wc2.append(ews_raws_C2[i])
        Wo2.append(ews_raws_O2[i])
        Wm2.append(ews_raws_Mg2[i])

    t = []
    for i in nums:
        ds, info = outflow_rate.load_data(timestep = i)
        t.append(info['t_myr'])


    fig = plt.figure()
    plt.xlabel('equivalent width of C 1548 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    p = plt.scatter(np.zeros(len(Wc1)), np.zeros(len(Wc1)), marker='.', c=t, cmap='viridis')  # Create the line
    plt.ylim(0.0325, 0.052)
    plt.xlim(0, 0.123)
    plt.title('EW from data directly')
    cbar = plt.colorbar()
    cbar.set_label('Time (Myr)')

    def buildmechart(i=int):
        x = np.concatenate((Wo1[0:i], np.zeros(len(Wo1)-i)), axis = None)
        y = np.concatenate((rates[0:i], np.zeros(len(Wo1)-i)), axis = None)
        #p = plt.scatter(Wo1[0:i], rates[0:i], marker = '.', c=t, cmap='viridis')
        p = plt.scatter(x, y, marker = '.', c=t, cmap='viridis')
        return p

    animator=ani.FuncAnimation(fig,buildmechart,frames=range(len(Wo1)),interval=100,repeat=False)
    animator.save(r'animation.gif', dpi = 700)
    plt.show()

    plt.scatter(Wo1, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of O 1032 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__O1032_pure.PNG')
    plt.show()

    plt.scatter(Wo2, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of O 1038 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__O1038_pure.PNG')
    plt.show()

    plt.scatter(Wc1, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of C 1548 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__C1548_pure.PNG')
    plt.show()

    plt.scatter(Wc2, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of C 1551 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__C1551_pure.PNG')
    plt.show()

    plt.scatter(Wm1, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of Mg 1240 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__Mg1240_pure.PNG')
    plt.show()

    plt.scatter(Wm2, rates, marker='.', c=t, cmap='viridis')
    plt.colorbar()
    plt.xlabel('equivalent width of Mg 1240.4 (A)')
    plt.ylabel('outflow rate (Msun/yr)')
    plt.title('EW from data directly')
    plt.savefig('rate_ew_funky_N_QSO_EW__Mg12404_pure.PNG')
    plt.show()

    return 0

# %%
