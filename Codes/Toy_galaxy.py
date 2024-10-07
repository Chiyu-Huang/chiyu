#%%
from Codes.outflow_rate import load_data
import numpy as np
from copy import deepcopy
import yt
from yt.units import YTArray
import h5py
import trident

ds, info = load_data(203, 'GTT_9pc')
# ad = ds.all_data()
# ad[('gas', 'density')] = yt.YTArray(np.ones(len(ad[('gas', 'density')])), 'g/cm**3')
# ds_new = yt.load('ah.h5')
# ad_new = ds_new.ray([0, 0, 0], [1, 1, 1])
ad = ds.ray([0.5,0.5,1], [0.5,0.5,0])

def _density_2(field, data):
    theta = data['index', 'spherical_theta']
    r = data['index', 'spherical_radius'].in_units('kpc')

    rv = data['gas', 'density'].in_units('g/cm**3').value

    # n_density = 3.16227766e-5
    # m_spe = 1.67262192e-24*16
    # density = n_density*m_spe

    density = 4e-20

    if len(np.shape(theta)) == 1:
        # Create a mask for the condition
        mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 20))
        
        # Apply the condition using NumPy's `where` function for vectorized assignment
        rv = np.where(mask, density, rv)
    return yt.YTArray(rv, 'g/cm**3')


def _vel(field, data):
    theta = data['index', 'spherical_theta']
    r = data['index', 'spherical_radius'].in_units('kpc')

    rv = data['gas', 'radial_velocity'].in_units('km/s').value

    # n_density = 3.16227766e-5
    # m_spe = 1.67262192e-24*16
    # density = n_density*m_spe

    v_rad = 250

    if len(np.shape(theta)) == 1:
        # Create a mask for the condition
        mask = ((theta > 2.44) | (theta < 0.698)) & ((r < 20))
        
        # Apply the condition using NumPy's `where` function for vectorized assignment
        rv = np.where(mask, v_rad, rv)
    return yt.YTArray(rv, 'km/s')

ds.add_field(
    ("gas", "vel_rad"),
    sampling_type="cell",
    function=_vel,
    units='km/s',
    force_override = True
)

#%%

print(ds.derived_field_list)
print(ad[('gas', 'vel_rad')])

# %%
