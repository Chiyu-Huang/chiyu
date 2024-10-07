import yt
import unyt
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

ds, info = outflow_rate.load_data(203, file = 'SINK_9pc')

def _temperature_mu_2(field, ds):
    theta = ds.all_data()['index', 'spherical_theta']
    r = ds.all_data()['index', 'spherical_radius'].in_units('kpc')

    rv = ds.all_data()['gas', 'temperature'].in_units('K').value

    if len(np.shape(theta)) == 1:
        # Create a mask for the condition
        mask = (theta > 2.44) | (theta < 0.698) | (r < 10)
        
        # Apply the condition using NumPy's `where` function for vectorized assignment
        rv = np.where(mask, 1e4, rv)

    # rv = yt.YTArray(np.ones(len(ds.all_data()['gas', 'density']))*1e4, 'K')
    return yt.YTArray(rv, 'K')

temps = _temperature_mu_2(('gas', 'density'), ds)
print(len(temps))
print(np.max(temps))

# %%
