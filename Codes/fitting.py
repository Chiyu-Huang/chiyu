import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib import colormaps
import scipy as sp
from astropy.cosmology import FlatLambdaCDM
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
from datetime import datetime
from collections import defaultdict



def fitting(wavelength, flux, lines, ray, crit, b, ts=None, bg = [], plot_fit = False, column = False):
    '''
    Function to fit absorption peaks given the specific lines (down to wavelength), then gives the equivalent widths or column densities

    Args:
        wavelength (1D array) : List of wavelength of the spectrum
        flux (1D array) : List of flux corresponding to the wavelength
        lines (tuple) : all the lines wanted to fit
        ray (trident simple ray object) : ray for which we find spectrum
        crit (str) : criteria for evaluating goondess of fit
        ts (int, optional) : time step of the file, if not none then save file (default = None)
        bg (1D array, optional) : list of background profile with same length as flux (default = None)
        plot_fit (bool, optional) : whether or not to plot the fit vs raw (default = False)
        column (bool, optional) : whether or not you are interested in column density or equivalent width (default = False)
    
    Returns:
        - dictionary - equivalent width of fit
        - dictionary - equivalent width of raw data
        - or is column == True
        - dictionary - column density of the raw data
    '''
    if len(flux) != len(bg):
        raise Exception('Background and flux do not have the same length')

    ad = ray.all_data()

    # n_den_dict = {
    #     'Mg II 1240' : ad[('gas', 'Mg_p1_number_density')],
    #     'Mg II 1240.4' : ad[('gas', 'Mg_p1_number_density')],
    #     'C IV 1548' : ad[('gas', 'C_p3_number_density')],
    #     'C IV 1551' : ad[('gas', 'C_p3_number_density')],
    #     'O VI 1032' : ad[('gas', 'O_p5_number_density')],
    #     'O VI 1038' : ad[('gas', 'O_p5_number_density')],
    #     'H I 1216' : ad[('gas', 'H_p0_number_density')],
    #     'H I 973' : ad[('gas', 'H_p0_number_density')],
    # }


    e_w_fit = {}
    e_w_raw = {}
    Ns = {}
    Ns_tri = {}
    
    '''
    This next bit iterates through all the lines specified. We isolate the lines
    by using the find_bounds function which returns the lower and upper index of
    the line and the index where the flux is a minimum. Then, to give some
    benefits of doubts the boundary is extended by certain amount of units.
    This is what wavelength/flux _cut represent: the data for an isolated line.

    It also take from find_bounds() the FWHM to estimate the Doppler b,
    the average background level, and the flux with background removed.
    The flux with background removed is the flux used to find the equivalent width.

    '''

    # plt.plot(wavelength, flux)
    # plt.xlim(980, 1060)
    # plt.show()
    for feature in lines:
        if feature == 'Mg II 1240' or feature == 'Mg II 1240.4':
            n_density = ad[('gas', 'Mg_p1_number_density')]
        elif feature == 'C IV 1548' or feature == 'C IV 1551':
            n_density = ad[('gas', 'C_p3_number_density')]
        elif feature == 'O VI 1032' or feature == 'O VI 1038':
            n_density = ad[('gas', 'O_p5_number_density')]
        elif feature == 'H I 1216' or feature == 'H I 973':
            n_density = ad[('gas', 'H_p0_number_density')]

        dl = ad[('gas', 'dl')]
        init_N = np.sum(n_density * dl)
        Ns_tri.update({feature : init_N.value})
        speciesDicts[feature]['init_N'] = init_N
        central_wl = speciesDicts[feature]['wavelength'][0]
        wl_range = speciesDicts[feature]['wl_range']

        # Redshift
        if all([x == 0 for x in n_density]):
            z_max_n = ad[('gas', 'velocity_los')].in_units('m/s').value[np.argmax(n_density * dl)]/(3e8)
        else:
            z_max_n = np.average(ad[('gas', 'redshift_dopp')].value, weights = n_density.value * dl.value)
        extra_lambda = central_wl*z_max_n
        central_wl = central_wl + extra_lambda
        ###########  Isolating the line  ###########
        if b > 4:
            speciesDicts['H I 1216']['wl_range'] = 10
        lower, upper, ind_min_flux, FWHM, bg_level, flux_no_bg = find_bounds(
            wavelength, flux,
            feature, bg, plot_reg = plot_fit)
        
        if flux_no_bg is None:
            e_w_raw.update({feature + '_is_upper_lim' : True})
            e_w_fit.update({feature + '_is_upper_lim' : True})
            upper = np.searchsorted(wavelength, central_wl+wl_range)
            lower = np.searchsorted(wavelength, central_wl-wl_range)
            upper_lim = equivalent_width(wavelength[lower:upper], flux[lower:upper])
            e_w_fit.update({feature : upper_lim})
            e_w_raw.update({feature : upper_lim})
            Ns.update({feature : 0})
            # print(feature, ' : Trough not pronounced enough, upper limit found = ', upper_lim)
            continue

        wavelength_cut = wavelength[lower : upper]
        flux_cut_normed = flux_no_bg
        '''This might need to be switched out as for now the second bg_cut_normed assume no QSO, but that may not be an issue'''
        # print('*********')
        # plt.plot(wavelength_cut, flux_cut_normed)
        # plt.plot(wavelength_cut, bg_cut_normed)
        # plt.show()
        # print('*********')
        ###########################################

        #############   Fitting while finding/using the optimal initial guess   #############
        # generating initial guesses
        maxb = 120*FWHM + 100    #empirical
        speciesDicts[feature]['maxb'] = maxb
        fitLim = (flux_no_bg[ind_min_flux])+ 0.005     #to ensure we only save the main feature
        complexLim = 1.1
        # trying 5 different 'init_b' parameters in speciesDicts
        # and savingthe fit with the closest fit (given by lowest
        # difference(fitted_flux, flux_cut))
        diff = []
        fits = []
        fitted = []
        # print('aaaa')
        # for j in range (3):
        #     maxb = 120*FWHM + (j+1)*33    #empirical
        #     speciesDicts[feature]['maxb'] = maxb
        #     for i in range (3):
        #         speciesDicts[feature]['init_b'] = (maxb/3)*(i+1)
        #         for k in range (3):
        #             speciesDicts[feature]['init_N'] = (init_N/1.5)*(k+1)
        #             fitted_lines, fitted_flux, cBounds = generate_total_fit(
        #                 wavelength_cut, flux_cut_normed,
        #                 [feature], speciesDicts, fitLim = fitLim,
        #                 first = 1, last = 1,
        #                 complexLim = complexLim
        #                 )
        #             if crit == 'sq':
        #                 diff.append(difference(fitted_flux, flux_cut_normed))
        #             elif crit == 'ew':
        #                 diff.append(np.abs(equivalent_width(wavelength_cut, fitted_flux)-\
        #                             equivalent_width(wavelength_cut, flux_cut_normed)))
        #             elif crit == 'chi':
        #                 diff.append(chisq(fitted_flux, flux_cut_normed))
        #             fits.append(fitted_flux)
        #             fitted.append(fitted_lines)
        #
        # fitted_lines = fitted[np.argmin(diff)]
        # print('bbbb')

        maxb = 120*FWHM + 100
        speciesDicts[feature]['maxb'] = maxb
        speciesDicts[feature]['init_b'] = (maxb/5)
        speciesDicts[feature]['init_N'] = (init_N)
        fitted_lines, fitted_flux = generate_total_fit(
            wavelength_cut, flux_cut_normed,
            [feature], speciesDicts, fitLim = fitLim,
            first = 1, last = 1,
            complexLim = complexLim
            )
        
        #####################################################################################

        if plot_fit:
            plt.plot(wavelength_cut, flux_cut_normed)
            plt.plot(wavelength_cut, fitted_flux)
            plt.title(feature)
            plt.show()
            print('min flux = ', np.min(flux_cut_normed))

        #finding equivalent width
        e_w_raw.update({feature + '_is_upper_lim' : False})
        e_w_fit.update({feature + '_is_upper_lim' : False})
        e_w_fit.update({feature : equivalent_width(wavelength_cut, fitted_flux)})
        e_w_raw.update({feature : equivalent_width(wavelength_cut, flux_cut_normed)})
        Ns.update({feature : np.sum(fitted_lines[feature]['N'])})

    return (e_w_fit, e_w_raw, Ns, Ns_tri)

def is_close(a, b, th):
    return np.abs(a-b)<=th


def find_bounds(wavelength, flux, feature, bg, b = 1.1, plot_reg = False):
    '''
    It estimates the lower and upper bound considered for the fit, the FWHM,
    the background level, and the flux with background removed through a linear
    approximation.

    Args:
        wavelength (1D array) : list of wavelength of the spectrum
        flux (1D array) : list of corresponding flux
        feature (tuple) : tuple of strings specifying the lines
        bg (1D array) : continuum level
        b (float) : for determining if there are multiple seperated peaks or not, helps estimate the FWHM
    
    Returns:
        - int - index of the lower bound of the raw data (input wavelength)
        - int - index of the upper bound of the raw data (input wavelength)
        - int - index of the peak in the raw data
        - flaot - Estimated FWHM of the line
        - float - maximum value of the continuum (averaged)
        - 1D array - raw flux with continuum removed (cut already)
    '''
    def linear(x, m, c):
        # m is the slope and c is y-intercept
        return m*x + c
    
    if len(flux) != len(bg):
        raise Exception('Background and flux do not have the same length')


    #   0th Quick cut to rolling average the fluxes within wl_range
    central_wavelength = speciesDicts[feature]['wavelength'][0]
    wl_range = speciesDicts[feature]['wl_range']
    ind_min_raw = np.searchsorted(wavelength, (central_wavelength-wl_range))
    ind_max_raw = np.searchsorted(wavelength, (central_wavelength+wl_range))
    #   Take an average to smooth out noise
    flux_avg = []
    wl_avg = []
    bg_avg = []
    for i in range(ind_min_raw, ind_max_raw):
        bin_size = int(np.round((ind_max_raw - ind_min_raw)/70))
        flux_avg.append(np.average(flux[i-bin_size:i+bin_size]))
        wl_avg.append(wavelength[i])
        bg_avg.append(np.average(bg[i-bin_size:i+bin_size]))

    # print('----------')
    # plt.plot(wl_avg, flux_avg, '.')
    # plt.plot(wavelength[ind_min_raw:ind_max_raw], flux[ind_min_raw:ind_max_raw])
    # plt.legend(('average', 'raw'))
    # plt.show()
    # print('----------')
    
    
    # Do the 1st cut again as above with the average data that we cubic spline.
    # Thus the rough_regions is the 1st cut, used for second cut later
    N = 100000
    cs = sp.interpolate.CubicSpline(wl_avg, flux_avg)
    xs = np.linspace(wl_avg[0], wl_avg[-1], N)
    flux_c = cs(xs)
    xs = xs
    ind_min_spli = np.searchsorted(xs, central_wavelength-wl_range)
    ind_max_spli = np.searchsorted(xs, central_wavelength+wl_range)
    flux_rough_region = flux_c[ind_min_spli:ind_max_spli]
    wl_rough_region = xs[ind_min_spli:ind_max_spli]

    # Linear approximation to remove background
    wl_start_avg = np.average(wl_rough_region[0:3])
    wl_end_avg = np.average(wl_rough_region[-3:-1])
    flux_start_avg = np.average(flux_rough_region[0:3])
    flux_end_avg = np.average(flux_rough_region[-3:-1])

    m_1st_cut = (flux_end_avg - flux_start_avg) / (wl_end_avg - wl_start_avg)
    c_1st_cut = flux_start_avg - m_1st_cut * wl_start_avg
    bg_linear = linear(wl_rough_region, m_1st_cut, c_1st_cut)
    flux_c_no_bg = flux_rough_region/bg_linear
    flux_no_bg = flux[ind_min_raw:ind_max_raw]/linear(wavelength[ind_min_raw:ind_max_raw], m_1st_cut, c_1st_cut)
    
    if plot_reg:
        print('OKOKOKOKOK')
        plt.plot(wl_rough_region, flux_rough_region)
        plt.plot(wl_rough_region, flux_c_no_bg)
        plt.plot(wavelength[ind_min_raw:ind_max_raw], flux_no_bg)
        plt.xlim(np.min(wl_rough_region), np.max(wl_rough_region))
        plt.title('before v after bg correction')
        plt.legend(('avg', 'corrected', 'raw'))
        #plt.ylim(0, 1.2)
        plt.show()
        print('OKOKOKOKOK')
        print('min flux = ', np.min(flux_no_bg))


    '''
    ind_upper = ind_min_flux.copy()
    ind_lower = ind_min_flux.copy()
    while not is_close(flux_c_no_bg[ind_upper], 1, 0.0000005) and flux_c_no_bg[ind_upper] < 1.00001:
        if ind_upper >= len(flux_c_no_bg)-1:
            break
        else:
            ind_upper += 1
    while not is_close(flux_c_no_bg[ind_lower], 1, 0.0000005) and flux_c_no_bg[ind_lower] < 1.00001:
        if ind_lower <= 0:
            break
        else:
            ind_lower -= 1
    '''
    if feature == 'H I 1216':
        ind_min_flux = np.searchsorted(wl_rough_region, central_wavelength)
    else:
        ind_min_flux = np.argmin(flux_c_no_bg)

    if ind_min_flux<=1 or ind_min_flux>=len(flux_c_no_bg)-1:
        '''None detection'''
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None
    # print('CUBIC')
    # plt.plot(wavelength[ind_min_raw:ind_max_raw], flux[ind_min_raw:ind_max_raw])
    # #plt.plot(xs, flux_c, '.')
    # # plt.plot(
    # #     xs[ind_cent-from_central:ind_cent+from_central],
    # #     flux_c[ind_cent-from_central:ind_cent+from_central], 'r-'
    # #     )
    # plt.plot(wl_rough_region, flux_c_no_bg)
    # plt.plot(
    #     wl_rough_region,
    #     flux_c_no_bg, 'r-'
    #     )
    # plt.vlines(wl_rough_region[np.argmin(flux_c_no_bg)], np.min(flux_c_no_bg), flux[ind_min_raw:ind_max_raw][0], 'C0')
    # plt.vlines(wl_rough_region[ind_min_flux], np.min(flux_c_no_bg), flux[ind_min_raw:ind_max_raw][0], 'C1')
    # # plt.xlim(wl_rough_region[0], wl_rough_region[-1])
    # plt.legend(('raw', 'cubic_avged', 'cent'))
    # plt.show()
    # print('CUBIC')
    
    
    if feature == 'H I 1216':
        ind_lower = np.argmax(flux_c_no_bg[0:ind_min_flux])
        ind_upper = np.argmax(flux_c_no_bg[ind_min_flux:-1])+ind_min_flux
    else:
        ind_lower = ind_min_flux.copy()
        ind_upper = ind_min_flux.copy()
        
        while not is_close(flux_c_no_bg[ind_upper], 1, 1e-7) and flux_c_no_bg[ind_upper] < 1.0000001:
            if ind_upper >= len(flux_c_no_bg)-2:
                break
            ind_upper += 2
        while not is_close(flux_c_no_bg[ind_lower], 1, 1e-7) and flux_c_no_bg[ind_lower] < 1.0000001:
            if ind_lower <= 2:
                break
            ind_lower -= 2

    if ind_lower >= ind_upper:
        # print('None 2')
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None


    m_processed = (np.average(flux_rough_region[ind_upper-3: ind_upper])-np.average(flux_rough_region[ind_lower: ind_lower+3]))\
        /(np.average(wl_rough_region[ind_upper-3: ind_upper])-np.average(wl_rough_region[ind_lower: ind_lower+3]))
    c_processed = np.average(flux_rough_region[ind_lower: ind_lower+3])-m_processed*np.average(wl_rough_region[ind_lower: ind_lower+3])

    # print('#########')
    # print(ind_lower, ind_upper, ind_min_flux)
    # plt.plot(wl_rough_region[ind_lower:ind_upper], flux_c_no_bg[ind_lower:ind_upper])
    # # plt.plot(wl_rough_region[ind_upper-500: ind_upper], flux_c_no_bg[ind_upper-500: ind_upper])
    # # plt.plot(wl_rough_region[ind_lower: ind_lower+500], flux_c_no_bg[ind_lower: ind_lower+500])
    # plt.title('Bit of spectrum fitted'+feature)
    # #plt.ylim(0, 1.2)
    # plt.show()
    # print('#########')

    flux_c_no_bg = flux_c_no_bg[ind_lower:ind_upper]
    wl_c_no_bg = wl_rough_region[ind_lower:ind_upper]
    ind_min_flux = np.argmin(flux_c_no_bg)
    # print(ind_min_flux, wl_c_no_bg[ind_min_flux])
    # plt.plot(wl_c_no_bg, flux_c_no_bg)
    # plt.show()

    if is_close(np.min(flux_c_no_bg), 1, 1e-6) or\
        is_close(ind_min_flux, 0, 1) or\
        is_close(ind_min_flux, len(flux_c_no_bg), 1):
        #Trough too shallow
        # print('None 5')
        # print(is_close(np.min(flux_c_no_bg), 1, 1e-6))
        # print(is_close(ind_min_flux, 0, 1))
        # print(is_close(ind_min_flux, len(flux_c_no_bg), 1))
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None

    # Finding the FWHM more precisely using rolling averages
    if len(flux_c_no_bg) <= 1:
        # print('None 3')
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None
    wl_min_flux = wl_c_no_bg[ind_min_flux]
    ind_fwhm_upper = ind_min_flux.copy()
    ind_fwhm_lower = ind_fwhm_upper.copy()
    while not is_close(flux_c_no_bg[ind_fwhm_upper],
                       0.5*(1+np.min(flux_c_no_bg)), 0.0000005):
        if ind_fwhm_upper >= len(flux_c_no_bg)-4*4:
            break
        ind_fwhm_upper += 4
    while not is_close(flux_c_no_bg[ind_fwhm_lower],
                       0.5*(1+np.min(flux_c_no_bg)), 0.0000005):
        if ind_fwhm_lower <= 4*4:
            break
        ind_fwhm_lower -= 4
    ind_fwhm_upper = np.searchsorted(wl_avg, wl_c_no_bg[ind_fwhm_upper])
    ind_fwhm_lower = np.searchsorted(wl_avg, wl_c_no_bg[ind_fwhm_lower])
    ind_extra = np.searchsorted(wavelength, wl_avg[0])
    ind_fwhm_upper += ind_extra
    ind_fwhm_lower += ind_extra
    if b < 1:
        FWHM = (wavelength[ind_fwhm_upper]-wavelength[ind_fwhm_lower])/5
    else:
        FWHM = (wavelength[ind_fwhm_upper]-wavelength[ind_fwhm_lower])/1.5
    # Translate the relevant values (ind_upper, lower, ind_min_flux, no_bg flux, etc)
    # into the raw data from the cubic splined averaged data
    
    
    ind_upper = np.searchsorted(wavelength, wl_rough_region[ind_upper])
    ind_lower = np.searchsorted(wavelength, wl_rough_region[ind_lower])
    ind_upper += 0  #  +/- 0 for safe measures
    ind_lower -= 0
    bg_linear = linear(wavelength[ind_lower:ind_upper], m_processed, c_processed)
    # flux_raw_no_bg = flux[ind_lower:ind_upper]/bg_linear
    flux_raw_no_bg = flux[ind_lower:ind_upper]
    if np.abs(ind_lower-ind_upper) < 10:
        # print('None 4')
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None
    ind_min_flux = np.searchsorted(wavelength[ind_lower:ind_upper], wl_min_flux)

    # plt.plot(wavelength[ind_lower:ind_upper], flux_raw_no_bg)
    # plt.title('no bg raw')
    # plt.show()

    if ind_min_flux >= len(flux_raw_no_bg):
        return 0, 0, np.min(flux_c_no_bg), 0, 0, None
    
    return (ind_lower, ind_upper, ind_min_flux, FWHM, np.max(bg_avg), flux_raw_no_bg)



def difference(a, b):
    '''returns the absolute value of the sum of the square differences'''
    diff = 0
    for i in range(len(a)):
        diff += (a[i]-b[i])**2
    return np.abs(diff)

def chisq(a, b):
    '''a is fit, b is data'''
    chi_squared = 0
    for i in range(len(a)):
        chi_squared += (a[i]-b[i])**2/(b[i])
    return chi_squared

def equivalent_width(wavelength, flux):
    Fs = flux
    dL = wavelength[1]-wavelength[0]
    integrand = []
    for i in range(len(Fs)):
        integrand.append((1 - Fs[i])*dL)
    ew = np.sum(integrand)
    return ew



'''
Here defined the fitting parameters of each lines, these are mostly pre-existing
ones from Trident, aside from 'wl_range'. 'wl_range' is the defines the range
of wavelength (from the central) we would first consider. Within this range we
would later find the bounds used for fitting using rolling averages.
'''
H_I_alpha_parameters = {'name':'H I 1216',
        'f': [4.1641e-01],
        'Gamma':[4.6986e+08],
        'wavelength':[1216],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 500, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':700,
        'init_N':1E14,
        'wl_range':25}
H_I_973_parameters = {'name':'H I 973',
        'f': [2.9006e-02],
        'Gamma':[1.2785e+07],
        'wavelength':[972.517],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14,
        'wl_range':1.8}
C_IV_1548_parameters = {'name':'C IV 1548',
        'f': [1.9e-1],
        'Gamma':[2.65e+08],
        'wavelength':[1548.187],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 100, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E15,
        'wl_range':1.5}
C_IV_1551_parameters = {'name':'C IV 1551',
        'f': [9.52e-02],
        'Gamma':[2.64e+08],
        'wavelength':[1550.772],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14,
        'wl_range':1.5}
O_VI_1032_parameters = {'name':'O VI 1032',
        'f': [1.33e-01],
        'Gamma':[4.16e+08],
        'wavelength':[1031.912],
        'numLines':1,
        'maxN': 1E15, 'minN':1E8,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E12,
        'wl_range':1.2}
O_VI_1038_parameters = {'name':'O VI 1038',
        'f': [6.60e-02],
        'Gamma':[4.09e+08],
        'wavelength':[1037.7],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14,
        'wl_range':0.8}
Mg_II_1240_parameters_1 = {'name':'Mg II 1240',
        'f': [6.21e-04],
        'Gamma':[1.35e+06],
        'wavelength':[1239.92],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14,
        'wl_range':0.26}
Mg_II_1240_parameters_2 = {'name':'Mg II 1240.4',
        'f': [6.21e-04],
        'Gamma':[1.35e+06],
        'wavelength':[1240.4],
        'numLines':1,
        'maxN': 1E22, 'minN':1E6,
        'maxb': 6000, 'minb':1,
        'maxz': 6, 'minz':0,
        'init_b':30,
        'init_N':1E14,
        'wl_range': 0.26}
speciesDicts = {'C IV 1548': C_IV_1548_parameters,
                'C IV 1551': C_IV_1551_parameters,
                'H I 1216': H_I_alpha_parameters,
                'H I 973': H_I_973_parameters,
                'O VI 1032': O_VI_1032_parameters,
                'O VI 1038': O_VI_1038_parameters,
                'Mg II 1240': Mg_II_1240_parameters_1,
                'Mg II 1240.4': Mg_II_1240_parameters_2
                }
