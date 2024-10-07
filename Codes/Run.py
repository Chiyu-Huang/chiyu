#%%
from Column_Rate import find_N_rates, N_plot_rates
from EW_rates import find_EW_rates, Ew_plot_rates

#%%
files = ['GTT_4pc', 'GTT_9pc', 'GTT_18pc', 'SINK_4pc', 'SINK_9pc']
files = ['GTT_9pc']
# 
find_EW_rates(files)
# find_N_rates(files)

# #%%
# files = ['SINK_9pc']
# # N_plot_rates(files, 'GTT')
# Ew_plot_rates(files, 'GTT')
# %%
