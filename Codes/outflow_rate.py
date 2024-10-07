#%% import libraries
import yt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib import colormaps
import scipy as sp
import pandas as pd
import os
import glob
from astropy.cosmology import FlatLambdaCDM
import pickle
import trident
import h5py
from trident.absorption_spectrum.absorption_spectrum_fit import generate_total_fit
from datetime import datetime
from unyt import unyt_quantity
yt.set_log_level(50)
#from ratatouille.readNsave import read_info, get_ramses_index

#%% Required functions

g2Msun = 1./1.9891e33
kpc2cm = 3.086e21
Myr2s = (303*31536000+97*31622400)/400*1e6  # Gregorian calendar: 400 year cycle of 303 common years + 97 leap years.

kB     = 1.3806200e-16          # Boltzmann constant [erg/K]

mH     = 1.6735575e-24          # Hydrogen atom mass [g]

def read_info(RamsesDir, timestep, alt_path=None):
    """ Returns info from the simulation (dict). """
    #savepath = get_alt_path(RamsesDir, alt_path=alt_path)
    #check_ratapath(savepath)
    if alt_path:
        savepath = alt_path+RamsesDir.split('rey', 1)[1]
    else:
        savepath = RamsesDir
    check_ratapath(savepath)
    path_info = f'{savepath}/ratadat/info_dict.pkl'
    if not os.path.exists(path_info):                # If first time, save general data
        info_dict = read_info_glob(RamsesDir, timestep=timestep)
        with open(path_info, 'wb') as file: 
            pickle.dump({'00000': info_dict}, file)
    with open(path_info, 'rb') as file:              # Otherwise, load common data
        all_info = pickle.load(file)
        info_dict = all_info.get('00000', {}).copy() # Important to copy, otherwise '00000' will be modified
        if timestep != 0:
            t_group = f'{timestep:05d}'
            if t_group in all_info:                  # And add timestep specific data.
                info_dict = all_info[t_group]
            else:
                info_dict = read_info_tstep(RamsesDir, timestep, info_dict)
                all_info[t_group] = info_dict
                with open(path_info, 'wb') as file:
                    pickle.dump(all_info, file)
    return info_dict

def check_ratapath(savepath):
    """ Check if ratadat folder exists, creates it otherwise. """
    ratapath = f'{savepath}/ratadat/'
    if not os.path.exists(ratapath):
        try: 
            os.makedirs(ratapath, exist_ok=True)  # Can create more than one level of directories.
        except PermissionError:
            raise PermissionError(f"Cannot mkdir {savepath}/ratadat. Add the argument alt_path='/general/path/to/allsimus'.")

def read_info_glob(RamsesDir, timestep=1):
    """ Read info parameters common to all timesteps of a given simulation. """
    tstep = f'{timestep:05d}'
    path_out = f'{RamsesDir}/output_{tstep}/'
    
    # Initialise some variables which might not be defined.
    info_dict = {}
    info_dict['new_format'] = get_format(f'{path_out}header_{tstep}.txt')                     # Get whether post-2017 or not
    for key in ['isH2', 'delayed_cooling', 'momentum_feedback']:
        info_dict[key] = False

    # Namelist
    info_dict  = param2dict(info_dict, get_nml(RamsesDir, timestep))                                    # Namelist
    info_dict  = param2dict(info_dict, f'{path_out}info_rt_{tstep}.txt')                                # info_rt_xxxxx.txt
    info_dict  = hydro2dict(info_dict, f'{path_out}hydro_file_descriptor.txt', info_dict['new_format']) # hydro_file_descriptor
    info_dict  = get_nvar  (info_dict, path_out, tstep)                                                 # (after hydro2dict)

    # Corrections
    info_dict['X_fraction'] = round(info_dict['X_fraction'], 6) # Correct precision error in info_rt_xxxxx.txt
    info_dict['Y_fraction'] = round(info_dict['Y_fraction'], 6)

    # Determine simulation type
    # Note: original way is the following BUT have to read info_xxxxx.txt (resp header_xxx.txt), and some variables there depends on timestep
    #       info_dict['is_cosmo'] = info_dict['omega_m']!=1.0 and info_dict['aexp']!=1 and info_dict['H0']!=1 
    #       info_dict['is_zoom']  = info_dict['is_cosmo'] and info_dict['ndm'] != (2*info_dict['levelmin'])*info_dict['ndim']
    info_dict['is_cosmo'] = info_dict['cosmo']
    info_dict['is_zoom']  = info_dict['is_cosmo'] and 'initfile(2)' in info_dict
    info_dict['is_cool_refine'] = 'cooling_refine' in info_dict and any(value != -1 for value in info_dict['cooling_refine'])
    info_dict['nvarnoadvect'] = info_dict['cooling_time_ivar'] if info_dict['is_cool_refine'] else 0

    return info_dict


def read_info_tstep(RamsesDir, timestep, info_dict):
    """ Read info parameters for a given timestep anc add some units. """
    tstep = f'{timestep:05d}'
    path_out = f'{RamsesDir}/output_{tstep}/'
    info_dict  = param2dict(info_dict, f'{path_out}info_{tstep}.txt')               # info_xxxxx.txt
    info_dict  = headr2dict(info_dict, f'{path_out}header_{tstep}.txt', info_dict['new_format']) # header_xxxxx.txt

    # Units
    info_dict['unit_m']      = info_dict['unit_d']*info_dict['unit_l']**3
    info_dict['unit_v']      = info_dict['unit_l']/info_dict['unit_t']
    info_dict['unit_P']      = info_dict['unit_d']*info_dict['unit_v']**2
    info_dict['unit_T2']     = mH/kB*info_dict['unit_v']**2
    info_dict['unit_nH']     = info_dict['X_fraction']*info_dict['unit_d']/mH
    info_dict['unit_nHe']    = info_dict['unit_nH'] * info_dict['Y_fraction']/info_dict['X_fraction'] * 0.25
    info_dict['cu2cm']       = info_dict['unit_l'] * info_dict['boxlen']
    info_dict['boxlen_cMpc'] = info_dict['boxlen']*info_dict['unit_l']/info_dict['aexp']/kpc2cm/1e3 if info_dict['is_cosmo'] else None
    info_dict['redshift']    = 1./info_dict['aexp']-1
    info_dict['t_myr']       = info_dict['time']*info_dict['unit_t']/Myr2s if not info_dict['is_cosmo']\
        else FlatLambdaCDM(H0=info_dict['H0'],Om0=info_dict['omega_m']).age(1/info_dict['aexp']-1).value*1e3
    return info_dict


def get_nml(RamsesDir, timestep):
    """Return the path to the nml."""
    # Check the output directory
    path_nml = f'{RamsesDir}/output_{timestep:05d}/namelist.txt'    
    if os.access(path_nml, os.R_OK):            # Check it's readable
        with open(path_nml, 'r') as file:
            for line in file:
                if '&RUN_PARAMS' in line:       
                    return path_nml
    # If the file is doesn't exist/is corrupted, check the RamsesDir
    nml_files = glob.glob(f"{RamsesDir}/*.nml")
    path_nml = next((f for f in nml_files if os.access(f, os.R_OK)), None) # First nml with read access
    if path_nml is None:
        raise Exception(f"No accessible .nml file found in {RamsesDir}.")
    return path_nml


def headr2dict(info_dict, file_path, new_format):
    """ Add header variables (i.e. number of particles and particle fields) to info_dict. """
    header = {
        "Total number of particles": "npart",       # Old format
        "Total number of dark matter particles": "ndm",
        "Total number of star particles": "nstar",
        "Total number of sink particles": "nsink",
        "DM": "ndm",                                # New format
        "star": "nstar",
    }
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if "Particle fields" in line:
                info_dict["particle_fields"] = next(file).strip().split()
                break  # No more info after particle fields
            parts = line.split()
            if new_format and len(parts) == 2 and parts[0] in header:
                info_dict[header[parts[0]]] = int(parts[1])
            elif not new_format and 'Total' in line:
                info_dict[header[line]] = int(next(file).strip())
    return info_dict


def hydro2dict(info_dict, path_hydrofd, new_format):
    """ Read hydro_file_descriptor and add the indexes of the hydro parameters to info_dict. """
    with open(path_hydrofd, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('#'): # skip the first in new format
                parts = line.split(',') if new_format else line.replace('#', ':').split(':')
                index1, index2 = (0, 1) if new_format else (-2, -1)
                if len(parts) >= 2:
                    info_dict[parts[index2].strip()] = int(parts[index1].strip())
    rename_map = {
        'density': 'idens',
        'velocity_x': 'ivel',
        'thermal_pressure': 'iPT',
        'pressure': 'iPT'
    }
    for old_name, new_name in rename_map.items():
        if old_name in info_dict:
            info_dict[new_name] = info_dict.pop(old_name)
    return info_dict


def get_nvar(info_dict, path_out, tstep):
    """ Try to get nvar from hydro_xxxxx.out00001, and in the hydro_file_descriptor if it fails. 
        First option is favoured as the hydro_file_descriptor can be badly coded. """
    hydro_file = f'{path_out}hydro_{tstep}.out00001'
    if os.path.exists(hydro_file):          # Try to read hydro_xxxxx.out00001
        with open(f'{path_out}hydro_{tstep}.out00001', 'rb') as f:  
            info_dict['nvar'] = np.fromfile(f, dtype=np.int32, count=5)[4]
    elif not 'dm' in path_out.lower():      # Else, skip if DM-only simulation and take the highest value from the hydro_file_descriptor.
        info_dict['nvar'] = info_dict[max((key for key in info_dict if 'scalar_' in key), key=lambda x: int(x.rsplit('_', 1)[1]), default=None)]
    return info_dict
def get_ramses_index(info, include_rt=False):                                                         # TODO: CR, MHD
    """ 
    Determine the index of the variables in the ramses output.
    """
    # Variables obtained in info
    dict_index = {'Density':       info['idens'],
                  'x-velocity':       info['ivel'],
                  'y-velocity':       info['ivel']+1                    if info['ndim']>=2 else None,
                  'z-velocity':       info['ivel']+2                    if info['ndim']>=3 else None,
                  'Pressure':        info['iPT']+info['rt_isIRtrap'],
                  'Metallicity':        info['iPT']+info['rt_isIRtrap']+1 if info['metal']   else None,
                  'xHI':      info['iIons']                     if info['isH2']    else None,
                  'xHII':     info['iIons']+info['isH2']        if info['rt']      else None,
                  'xHeII':    info['iIons']+info['isH2']+1      if info['rt']      else None,
                  'xHeIII':   info['iIons']+info['isH2']+2      if info['rt']      else None,
                }

    # Additional variables
    counter = info['iPT'] + info['rt_isIRtrap'] + info['metal'] + info['nIons']  # info['isH2'] already counted in nIons
    if info['delayed_cooling']:   dict_index['DC_var'],         counter = counter + 1, counter + 1  # delayed cooling
    if info['momentum_feedback']: dict_index['KR_turb'],        counter = counter + 1, counter + 1  # patch mom2 by Kretschmer
    if info['is_zoom']:           dict_index['zoom_var'],       counter = counter + 1, counter + 1  # zoom-in simulations
    if info['is_cool_refine']:    dict_index['cooling_length'], counter = counter + 1, counter + 1  # cooling length refinement

    # RT variables: N_photons and flux_nx/y/z for each group n.
    if info['rt'] and include_rt:
        for i in range(1, info['nGroups'] + 1):
            dict_index[f'N_photons{i}'] = info['nvar'] + 1 + (i - 1) * 4 if info['nGroups'] >= i else None
            for axis in ['x', 'y', 'z']:
                dict_index[f'flux_{i}{axis}'] = info['nvar'] + 2 + "xyz".index(axis) + (i - 1) * 4 if info['nGroups'] >= i else None

    return dict_index

def get_values(var, dictionary):    # TODO: restrict to one value, not T or temperature. Then, remove this function -> everything will be more readable.
    """ Helper function to be able to format dictionary keys as tuples. """
    for key in dictionary:
        if var in key:
            return dictionary[key]
    return None


def get_format(path_hdr):
    """ Returns whether the code is formatted following RAMSES post-2017 or not. """
    with open(path_hdr, 'r') as file:
        first_line = file.readline().strip()
        if 'Total number of particles' in first_line: return False
        elif '#      Family     Count' in first_line: return True
        else: raise ValueError("Unrecognized header format.")


def fort2py(value):
    """ Convert Fortran values to Python. """
    if value.lower() in ('.true.', '.false.'):  return value.lower() == '.true.'    # Booleans
    elif (value.startswith("'") and value.endswith("'")): return value[1:-1]        # Strings
    else:
        try: return float(value) if 'e' in value or '.' in value else int(value)    # Numerals
        except ValueError:  return value                                            # Other strings


def param2dict(info_dict, path_file):
    """ General function to get parameters from a file and add them to a dictionary. """
    with open(path_file, 'r') as file:
        for line in file:
            line = line.split('!')[0].strip()       # Remove comments and whitespace
            if '=' in line:
                key, value = map(str.strip, line.split('='))
                if key!='movie_vars_txt':           # Replace 'd' with 'e' for float compatibility
                    value = value.replace('d', 'e')
                if ',' in value or '*' in value:    # Convert lists
                    value_list = []
                    for element in value.split(','):
                        if '*' in element:
                            count, val = element.split('*')
                            value_list.extend([int(float(val.strip()))] * int(count))
                        else:
                            value_list.append(fort2py(element.strip()))
                    value = value_list
                else:
                    value = fort2py(value)          # Convert single values
                info_dict[key] = value
            elif "DOMAIN" in line or "Photon group properties" in line:
                break
    return info_dict


#%% Loading data

def load_data(timestep, file = 'GTT_9pc'):
    RamsesDir = '/minkewhale/kimm/rey/G8_normal/'+str(file)
    alt_path='/minkewhale/chiyu'
    timestep = timestep
    info = read_info(RamsesDir, timestep, alt_path=alt_path)
    #print(info['t_myr'])
    dict_index = get_ramses_index(info)

    extra_particle_fields = [("particle_birth_time", "float64"), ("particle_metallicity", "float64"), ("particle_imass", "float64")]
        #extra fields included in the data that yt does not automatically detect (I assume to save time?)
    center_dat=[0.5, 0.5, 0.5]
    rad_dat = 0.2739726724728318
    bbox = [[c-rad_dat for c in center_dat], [c+rad_dat for c in center_dat]]
    ds = yt.load(f'{RamsesDir}/output_{timestep:05d}', extra_particle_fields=extra_particle_fields, fields = list(dict_index.keys()), bbox = bbox)
    
    return (ds, info)

#%% Defining outflow rate

def to_kpc(cl):
    return cl/(3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22)) #have to change accordingly with the specific file

def to_cl(kpc):
    return kpc*3.08567758128E+21/(0.150000000000000E+03*0.308567758128200E+22) #have to change accordingly with the specific file

def is_close(a, b, th):
    return np.abs(a-b)<=th


def outflow_rate(ds, upper, lower, mode = 'plate', ray = None, start_end = None, spe = None, temp_range = None):
    #TODO: try different outflow rate definitions
    # 1. outflow rate at the absorption
    # 2. same plate but at height = impact params of light ray
    # 3. sphere at radius = impact params
    '''
    Function to find the outflow rate given the upper and lower lid

    Args:
        ds (yt data object) : yt data (can be loaded from a file)
        upper (float) : The distance to upper lid to calculate outflow rate from in kpc
        upper (float) : The distance to upper lid to calculate outflow rate from in kpc

    Returns:
        float : Outflow rate in Msun/yr
    '''
    if mode == 'plate':
        #   upper lid
        r_u = ds.r[:, :, (upper+to_kpc(0.5), "kpc")]
        # start_t = datetime.now()
        rho_u = r_u[('gas', 'density')].in_units('Msun/pc**3')
        v_z_u = r_u[('gas', 'velocity_z')].in_units('pc/yr')
        dx = r_u[('gas', 'dx')].in_units('pc')
        dy = r_u[('gas', 'dy')].in_units('pc')
        dA_u = dx*dy
        # print('plate time = ', datetime.now()-start_t)
        #Ensure only accounting for outflows and no inflows by setting rho to 0 (no contribution to sum)
        mask = (v_z_u < 0)
        rho_u = np.where(mask, rho_u, unyt_quantity(0, units = 'Msun/pc**3'))
        #abs to ensure it is positive
        rate_u = np.abs(np.sum([rho_u[i]*v_z_u[i]*dA_u[i] for i in range (len(dA_u))]))

        #   lower lid
        r_l = ds.r[:, :, (lower+to_kpc(0.5), "kpc")]
        rho_l = r_l[('gas', 'density')].in_units('Msun/pc**3')
        v_z_l = r_l[('gas', 'velocity_z')].in_units('pc/yr')
        dy = r_l[('gas', 'dy')].in_units('pc')
        dx = r_l[('gas', 'dx')].in_units('pc')
        dA_l = dx*dy
        mask = (v_z_l >  0)
        rho_l = np.where(mask, rho_l, unyt_quantity(0, units = 'Msun/pc**3'))
        rate_l = np.abs(np.sum([rho_l[i]*v_z_l[i]*dA_l[i] for i in range (len(dA_l))]))
        return rate_l + rate_u  #in units (Msun/yr)

    if mode == 'sphere':

        ad = ds.sphere([0.5, 0.5, 0.5], (upper+1, "kpc"))
        # ad = ds.cut_region(sp, ["obj[('index', 'spherical_radius')].in_units('kpc') > " + str((upper-0.01))])
        dx = ad[('gas', 'dx')].in_units('pc')
        dy = ad[('gas', 'dy')].in_units('pc')
        dz = ad[('gas', 'dz')].in_units('pc')
        d_dist = np.sqrt(dx*dx + dy*dy + dz*dz)
        Radius = ad[('gas', 'radius')].in_units('pc')
        vel_r = ad[('gas', 'vel_rad')].in_units('kpc/yr')
        rho = ad[('gas', 'density')].in_units('Msun/pc**3')
        vols = ad[('gas', 'volume')].in_units('kpc**3')
        mask = (is_close(Radius.value, np.abs(upper)*1000, d_dist)) & (vel_r > 0)
        rates = np.where(mask, np.abs(rho*vel_r*vols/(2*d_dist)).value, 0)


        return np.sum(rates)
        
            
    if mode == 'ray':
        Msun = 1.989e33 #in grams
        pc = 3.08567758e18
        M_proton = 1.67262192e-27
        ad = ds.ray(start_end[0], start_end[1])
        rho = ad[('gas', 'density_2')]
        if spe == 'H I':
            n = ray.all_data()[('gas', 'H_p0_number_density')].in_units('1/pc**3').value
            rho = n*M_proton*1/Msun
        elif spe == 'O VI':
            n = ray.all_data()[('gas', 'O_p5_number_density')].in_units('1/pc**3').value
            rho = n*M_proton*16/Msun
        elif spe == 'C IV':
            n = ray.all_data()[('gas', 'C_p3_number_density')].in_units('1/pc**3').value
            rho = n*M_proton*12/Msun
        elif spe == 'Mg II':
            n = ray.all_data()[('gas', 'Mg_p1_number_density')].in_units('1/pc**3').value
            rho = n*M_proton*24.3/Msun

        v_z = ray.all_data()[('gas', 'relative_velocity_z')].in_units('pc/yr')
        z = ad[('gas', 'z')].in_units('pc')
        dA = ad[('gas', 'dx')].in_units('pc')*ad[('gas', 'dy')].in_units('pc')
        #Ensure only accounting for outflows and no inflows by setting rho to 0 (no contribution to sum)
        mask = (((z>75) & (v_z<0)) | ((z<75) & (v_z>0))) & (n > np.percentile(n, 10))
        rho = np.where(mask, 0, rho)

        #abs to ensure it is positive
        rate = np.sum(  np.abs([rho[i]*v_z[i]*dA[i] for i in range (len(dA))])  )
        return rate
        
    if mode == 'ray_temp':
        #TODO: finish writing this
        Msun = 1.989e33
        pc = 3.08567758e18
        M_proton = 1.67262192e-27
        ad = ds.ray(start_end[0], start_end[1])
        # print(start_end)

        rho = ad[('gas', 'density_2')].in_units('g/pc**3')
        T = ad[('gas', 'temperature_mu')].in_units('K')
        theta = ad[('gas', 'theta')]
        vel_rad = ad[('gas', 'vel_rad')].in_units('pc/yr')
        dx = ad[('gas', 'dx')].in_units('pc')
        dy = ad[('gas', 'dy')].in_units('pc')
        z = ad[('gas', 'z')].in_units('kpc')

        np.set_printoptions(suppress=True)
        # print(np.array(theta))

        rho = sort_by_sight_line(ad['t'], rho)
        T = sort_by_sight_line(ad['t'], T)
        theta = (sort_by_sight_line(ad['t'], theta))
        vel_rad = sort_by_sight_line(ad['t'], vel_rad)
        dx = sort_by_sight_line(ad['t'], dx)
        dy = sort_by_sight_line(ad['t'], dy)
        z = sort_by_sight_line(ad['t'], z)

        dA = dx*dy
        v_z = vel_rad*np.cos(theta)
        

        if spe == 'H I':
            n = ray.all_data()[('gas', 'H_p0_number_density')].in_units('1/pc**3').value
        elif spe == 'O VI':
            n = ray.all_data()[('gas', 'O_p5_number_density')].in_units('1/pc**3').value
        elif spe == 'C IV':
            n = ray.all_data()[('gas', 'C_p3_number_density')].in_units('1/pc**3').value
        elif spe == 'Mg II':
            n = ray.all_data()[('gas', 'Mg_p1_number_density')].in_units('1/pc**3').value
        
        #Ensure only accounting for outflows and no inflows by setting rho to 0 (no contribution to sum)
        if temp_range == 'hot':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T > 1e6)
        elif temp_range == 'warm':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T <= 1e6) & (T >= 1e5)
        elif temp_range == 'cold':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T < 1e4)
        else:
            # mask = (((z>75) & (v_z>0)) | ((z<75) & (v_z<0))) & (n > np.percentile(n, 10))
            mask = (((z>75) & (v_z>0)) | ((z<75) & (v_z<0))) & (n > 0)
        rho = np.where(mask, rho, 0)
        
        #abs to ensure it is positive
        rate = np.sum([np.abs(rho[i]*v_z[i]*dA[i]) for i in range (len(dA))])/Msun
        return rate
    

    if mode == 'ahoj':
        #TODO: finish writing this
        Msun = 1.989e33
        pc = 3.08567758e18
        M_proton = 1.67262192e-27
        ad = ds.ray(start_end[0], start_end[1])
        # print(start_end)

        rho = ad[('gas', 'density_2')].in_units('g/pc**3')
        T = ad[('gas', 'temperature_mu')].in_units('K')
        theta = ad[('gas', 'theta')]
        vel_rad = ad[('gas', 'vel_rad')].in_units('pc/yr')
        dx = ad[('gas', 'dx')].in_units('pc')
        dy = ad[('gas', 'dy')].in_units('pc')
        z = ad[('gas', 'z')].in_units('kpc')

        np.set_printoptions(suppress=True)
        # print(np.array(theta))

        rho = sort_by_sight_line(ad['t'], rho)
        T = sort_by_sight_line(ad['t'], T)
        theta = (sort_by_sight_line(ad['t'], theta))
        vel_rad = sort_by_sight_line(ad['t'], vel_rad)
        dx = sort_by_sight_line(ad['t'], dx)
        dy = sort_by_sight_line(ad['t'], dy)
        z = sort_by_sight_line(ad['t'], z)

        dA = dx*dy
        v_z = vel_rad*np.cos(theta)
        

        if spe == 'H I':
            n = ray.all_data()[('gas', 'H_p0_number_density')].in_units('1/pc**3').value
        elif spe == 'O VI':
            n = ray.all_data()[('gas', 'O_p5_number_density')].in_units('1/pc**3').value
        elif spe == 'C IV':
            n = ray.all_data()[('gas', 'C_p3_number_density')].in_units('1/pc**3').value
        elif spe == 'Mg II':
            n = ray.all_data()[('gas', 'Mg_p1_number_density')].in_units('1/pc**3').value
        
        #Ensure only accounting for outflows and no inflows by setting rho to 0 (no contribution to sum)
        if temp_range == 'hot':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T > 1e6)
        elif temp_range == 'warm':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T <= 1e6) & (T >= 1e5)
        elif temp_range == 'cold':
            mask = (((z>0) & (v_z>0)) | ((z<0) & (v_z<0))) & (n > np.percentile(n, 10)) & (T < 1e4)
        else:
            # mask = (((z>75) & (v_z>0)) | ((z<75) & (v_z<0))) & (n > np.percentile(n, 10))
            mask = (((z>75) & (v_z>0)) | ((z<75) & (v_z<0))) & (n > 0)
        rho = np.where(mask, rho, 0)
        # print('z, ', z)
        # print('vr,', vel_rad)
        # print('vz, ', v_z)
        # plt.plot(z, n)
        # plt.show()
        
        #abs to ensure it is positive
        rate = np.sum( [np.abs(rho[i]*v_z[i]*dA[i]) for i in range (len(dA))])/Msun
        return rate



def sort_by_sight_line(t, dat):
    Big_list = np.array([t, dat])
    return(Big_list[:, Big_list[0].argsort()][1])

# %%
