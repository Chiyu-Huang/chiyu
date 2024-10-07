#%% import libraries
from Codes.fitting import fitting
from Codes.outflow_rate import outflow_rate, load_data
import trident
import matplotlib.pyplot as plt

ds, a = load_data(timestep = 203)
trident.add_ion_fields(ds, ions=['Mg II'])
trident.add_ion_fields(ds, ions=['O VI'])
trident.add_ion_fields(ds, ions=['C IV'])

#   Creating the light ray
ray_start, ray_end = [1, 0.5, 1], [0, 0.5, 0.1]
ray = trident.make_simple_ray(ds,
                            start_position=ray_start,
                            end_position=ray_end,
                            data_filename="ray.h5",
                            lines = ['H'])


#   Generating spectrum
sg = trident.SpectrumGenerator(lambda_min = 1100, lambda_max = 1600, dlambda=10/1000)
sg.make_spectrum(ray, lines = ['H'])
# sg.add_milky_way_foreground()
sg.add_qso_spectrum()
wavelength = sg.lambda_field.value
flux = sg.flux_field
plt.plot(wavelength, flux)
plt.xlim(1214, 1220)
plt.show()

# %%
