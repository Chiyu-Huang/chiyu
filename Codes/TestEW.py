#%%
import yt
import matplotlib.pyplot as plt
import numpy as np
import trident
import Codes.outflow_rate as outflow_rate
import Codes.fitting as fitting

def rate_ew_TEST(time_step, start, end, lines, crit, plot = False, Noise = True, QSO = True, MilkyWay = True):
    '''
    Function that finds the outflow rate and the equivalent width given a ray

    Args:
        time_step (int) : Time step of the simulation
        start (1D array) : Coordinate of the start in code units (0.5 being the centre)
        end (1D array) : Coordinate of the end similarly
        lines (tuple) : The specific lines you want to fit
        plot (Bool) : Whether or not to plot the ray trajectory (default False)

    Returns:
        float : Outflow rate in Msun/yr
        dictionary : Fitted equivalent width of all populations specified
        dictionary : Raw data equivalent widths
    '''
    #   Determine outflow rate
    ds = outflow_rate.load_data(timestep = time_step)


    '''Add ion fields, creating the light ray, and finding the equivalent width'''
    ### No noise case
    trident.add_ion_fields(ds, ions=['Mg II'])
    trident.add_ion_fields(ds, ions=['O VI'])
    trident.add_ion_fields(ds, ions=['C IV'])

    #   Creating the light ray
    ray_start = start
    ray_end = end
    ray = trident.make_simple_ray(ds,
                                start_position=ray_start,
                                end_position=ray_end,
                                data_filename="ray.h5",
                                lines = ['O', 'C', 'H', 'S', 'N', 'Mg'])
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
    
    #   Generating spectrum
    sg = trident.SpectrumGenerator(lambda_min = 800, lambda_max = 1600, dlambda=10/1000)
    sg.make_spectrum(ray, lines = 'all')
    # sg.add_gaussian_noise(20)
    # sg.add_qso_spectrum()
    #sg.add_milky_way_foreground()
    wavelength = sg.lambda_field.value
    flux = sg.flux_field
    # plt.plot(wavelength, flux)
    # plt.title('pure')
    # plt.xlim(1545, 1552)
    # plt.show()
    #   specifying background then fit
    background = np.ones(len(flux)) #need to change later
    ew_fit, ew_raw = fitting.fitting_op(wavelength, flux, lines, bg=background, ray=ray, ts = time_step, crit=crit)
    print('all, ', ew_fit)


    ### Now for the noisy cases
    sgN = sg
    sgN.add_gaussian_noise(20)
    sgN.add_qso_spectrum()
    # sgN.add_milky_way_foreground()
    wavelengthN = sgN.lambda_field.value
    fluxN = sgN.flux_field
    #   specifying background then fit
    background = np.ones(len(fluxN)) #need to change later
    ew_fit_N, ew_raw = fitting.fitting_op(wavelengthN, fluxN, lines, bg=background, ray=ray, ts = str(time_step)+'_noisy', crit=crit)
    # plt.plot(wl_cut_noisy, flux_raw_noisy)
    # plt.plot(wl_cut, flux_raw_pure)
    # plt.title('N')
    # plt.show()
    # print('not all, ', ew_fit_N)
    # print('O = ', ew_fit['O VI 1032']+ew_fit['O VI 1038'] - ew_fit_N['O VI 1032'] + ew_fit_N['O VI 1038'])
    # print('Mg = ',ew_fit['Mg II 1240']+ew_fit['Mg II 1240.4'] - ew_fit_N['Mg II 1240']+ew_fit_N['Mg II 1240.4'])
    # print('C = ',ew_fit['C IV 1548']+ew_fit['C IV 1551'] - ew_fit_N['C IV 1548']+ew_fit_N['C IV 1551'])

    return (ew_fit_N, ew_fit)

Ew__O = []
Ew_N__O = []
Ew__C = []
Ew_N__C = []
Ew__Mg = []
Ew_N__Mg = []
nums = [1,5,9,13,25,29,41,53,73,89,100,105, 108,112, 120,135,143,150,165,173,185,197,203]
#nums = [120]
diag_x = np.arange(0, 0.13, 0.005)
for i in nums:
    print(i)
    ew_fit, ew_fit_N = rate_ew_TEST(i, [1, 0.5, 1], [0, 0.5, 0.1],
                                    ['O VI 1032', 'O VI 1038', 'Mg II 1240',\
                                        'Mg II 1240.4', 'C IV 1548', 'C IV 1551'],
                                        plot = False, crit = 'sq')
    Ew__O.append(ew_fit['O VI 1032']+ew_fit['O VI 1038'])
    Ew__Mg.append(ew_fit['Mg II 1240']+ew_fit['Mg II 1240.4'])
    Ew__C.append(ew_fit['C IV 1548']+ew_fit['C IV 1551'])
    Ew_N__O.append(ew_fit_N['O VI 1032'] + ew_fit_N['O VI 1038'])
    Ew_N__Mg.append(ew_fit_N['Mg II 1240']+ew_fit_N['Mg II 1240.4'])
    Ew_N__C.append(ew_fit_N['C IV 1548']+ew_fit_N['C IV 1551'])
    # plt.plot(diag_x, diag_x, '-')
    # plt.plot(Ew__C, Ew_N__C, 'g.')
    # plt.show()
plt.savefig('waste.png')

#%%
diag_x = np.arange(0, 0.18, 0.005)
plt.plot(diag_x, diag_x, '-')
plt.plot(Ew__O, Ew_N__O, 'r.')
plt.plot(Ew__Mg, Ew_N__Mg, 'b.')
plt.plot(Ew__C, Ew_N__C, 'g.')
plt.title('Noise+QSO+MW vs pure')
plt.xlabel('Noise_QSO_MW (A)')
plt.ylabel('pure (A)')
plt.legend(('line of equal width', 'O', 'Mg', 'C'))
plt.xlim(0, 0.12)
plt.ylim(0, 0.12)
plt.savefig('Noise_sq_one_997.png')
plt.show()

# %%
