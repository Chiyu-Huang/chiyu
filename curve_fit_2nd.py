#%%
import numpy as np
import matplotlib.pyplot as plt
from fitting import fitting_from_spec
from fitting import fitting
import json
#%%
def equivalent_width(wavelength, flux):
    Fs = flux
    dL = wavelength[1]-wavelength[0]
    integrand = []
    for i in range(len(Fs)):
        integrand.append((1 - Fs[i])*dL)
    ew = np.sum(integrand)
    return ew

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

wls = {}
wlps = {}
fls = {}
flps = {}
ws = {}
Ns = {}
nums = [53,73,89,100,105, 108,112, 120,135,143,150,165,173,185,197,203]
for i in nums:
    wl, fl, widths, columns = np.loadtxt('CSVs/C IV 1551_'+str(i)+'_noisy.csv', delimiter=',', unpack = True)
    wlp, flp, widthsp, columnsp = np.loadtxt('CSVs/C IV 1551_'+str(i)+'.csv', delimiter=',', unpack = True)
    name = str(i)+'_1551'
    wls[name] = wl
    wlps[name] = wlp
    fls[name] = fl
    flps[name] = flp
    ws[name] = widths[0]
    Ns[name] = columns[0]

for i in nums:
    wl, fl, widths, columns = np.loadtxt('CSVs/C IV 1548_'+str(i)+'_noisy.csv', delimiter=',', unpack = True)
    wlp, flp, widthsp, columnsp = np.loadtxt('CSVs/C IV 1548_'+str(i)+'.csv', delimiter=',', unpack = True)
    name = str(i)+ '_1548'
    wls[name] = wl
    wlps[name] = wlp
    fls[name] = fl
    flps[name] = flp
    ws[name] = widths[0]
    Ns[name] = columns[0]



wls_o = {}
wlps_o = {}
fls_o = {}
flps_o = {}
ws_o = {}
Ns_o = {}
nums = [41,73,89,100,105, 108,112, 120,135,143,150,165,173,185,197,203]
for i in nums:
    wl, fl, widths, columns = np.loadtxt('CSVs/O VI 1032_'+str(i)+'_noisy.csv', delimiter=',', unpack = True)
    wlp, flp, widthsp, columnsp = np.loadtxt('CSVs/O VI 1032_'+str(i)+'.csv', delimiter=',', unpack = True)
    name = str(i)+ '_1032'
    wls_o[name] = wl
    wlps_o[name] = wlp
    fls_o[name] = fl
    flps_o[name] = flp
    ws_o[name] = widths[0]
    Ns_o[name] = columns[0]
for i in nums:
    wl, fl, widths, columns = np.loadtxt('CSVs/O VI 1038_'+str(i)+'_noisy.csv', delimiter=',', unpack = True)
    wlp, flp, widthsp, columnsp = np.loadtxt('CSVs/O VI 1038_'+str(i)+'.csv', delimiter=',', unpack = True)
    name = str(i)+ '_1038'
    wls_o[name] = wl
    wlps_o[name] = wlp
    fls_o[name] = fl
    flps_o[name] = flp
    ws_o[name] = widths[0]
    Ns_o[name] = columns[0]

plt.hist(ws.values(), bins = 10)
plt.show()


fwhms = {}
ns = {}
maxns = {}
maxbs = {}
initbs = {}

par_dict = {
    'FWHM' : fwhms,
    'N' : ns,
    'maxN' : maxns,
    'maxb' : maxbs,
    'init_b' : initbs
}

o_fwhms = {}
o_ns = {}
o_maxns = {}
o_maxbs = {}
o_initbs = {}

par_dict_o = {
    'FWHM' : o_fwhms,
    'N' : o_ns,
    'maxN' : o_maxns,
    'maxb' : o_maxbs,
    'init_b' : o_initbs
}

def min_diff(ts, FWHM, N, wl, fl, wlp, flp, par_dict, crit = 'ew', n_tries = 7):
    diff = np.full((n_tries,n_tries, n_tries), float(300.0))
    print(FWHM, N)
    for i in range(n_tries):
        print(i)
        maxb = 120*FWHM + (100/n_tries)*(i+1)
        for j in range(n_tries):
            init_b = (maxb/n_tries)*(j+1)
            for k in range(n_tries):
                maxN = (N/n_tries)*(k+1)
                fitted_flux, raw_flux = fitting_from_spec(maxb, init_b, maxN, 'C IV 1551', wl, fl, wlp, flp, crit = 'ew')
                #Find percentage difference in ew of pure and fit
                dew = 100*np.abs(equivalent_width(wl, fitted_flux)-\
                                equivalent_width(wlp, flp))/(equivalent_width(wlp, flp))
                plt.plot(wl, raw_flux)
                plt.plot(wl, fitted_flux)
                plt.plot(wlp, flp)
                plt.legend(('raw', 'fit', 'pure'))
                plt.show()
                diff[i, j, k] = dew
                if k >= 4:
                    if diff[i, j, k] > diff[i, j, k-1] or diff[i, j, k] > diff[i, j, k-2] or\
                    diff[i, j, k] > diff[i, j, k-3] or diff[i, j, k] > diff[i, j, k-4]:
                        break


    min_ind = np.unravel_index(diff.argmin(), diff.shape)

    b = 120*FWHM + (100/n_tries)*(min_ind[0]+1)
    i_b = (b/n_tries)*(min_ind[1]+1)
    maxN = (2*N/n_tries)*(min_ind[2]+1)
    par_dict_o['FWHM'][ts] = FWHM
    par_dict_o['N'][ts] = N
    par_dict_o['maxN'][ts] = maxN
    par_dict_o['maxb'][ts] = b
    par_dict_o['init_b'][ts] = i_b

    return diff

print(ws.keys())



# %%
nums_str = [str(i) for i in ws.keys()]
for ts in nums_str:
    diff = min_diff(ts, ws[ts], Ns[ts], wls[ts], fls[ts], wlps[ts], flps[ts], par_dict, crit = 'ew', n_tries = 3)
print(diff)

json.dump( par_dict, open( "par_dict.json", 'w' ) )

nums_str = [str(i) for i in ws_o.keys()]
for ts in nums_str:
    diff = min_diff(ts, ws[ts], Ns[ts], wls[ts], fls[ts], wlps[ts], flps[ts], par_dict, crit = 'ew', n_tries = 3)
print(diff)

json.dump( par_dict_o, open( "par_dict_o.json", 'w' ) )

# %%
par_dict = json.load( open( "par_dict.json" ) )
a = list(par_dict['FWHM'].values())
key_closest_fwhm = list(par_dict['FWHM'].keys())[min(range(len(a)), key=lambda i: abs(a[i]-0.02))]
maxb = par_dict['maxb'][key_closest_fwhm]
init_b = par_dict['init_b'][key_closest_fwhm]

Ns = list(par_dict['N'].values())
key_closest_N = list(par_dict['N'].keys())[min(range(len(Ns)), key=lambda i: abs(Ns[i]-0.02))]
init_N = par_dict['maxN'][key_closest_N]

