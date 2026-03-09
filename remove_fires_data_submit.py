import xarray as xr
import glob
import os
import pdb
import numpy as np
import sys

#sectors = ['Agriculture', 'Energy_Use', 'Off_Road_Transport', 'Shipping', 'TOTALS_No_ANTHRO_FIRES', 'Waste', 'Aviation', 'Industry', 'Road_Transport', 'TOTALS_No_FIRES']
#sys.argv[1]
sectors = ['TOTALS_No_FIRES']
path = '/ptmp/acasallasg/Sectors/'
list_files = sorted(glob.glob(f'{path}TOTALS/wrfchemi*')) #wrfchemi for no anthro

for sector in sectors:
    arch = f'/ptmp/acasallasg/Sectors/{sector}/'
    for filename in list_files:
        #nombre = filename[-36:] #fire
        nombre = filename[-16:] #chem
        print(f'----------> {nombre}') 
        ds = xr.open_dataset(filename)
        #os.system(f"rm {filename}")
        for var in ds.data_vars:
            if np.issubdtype(ds[var].dtype, np.number):  # Is the variable numeric?
                 ds[var] = ds[var] * 0 ## the 0 removes the fires
        ds.to_netcdf(f'{arch}{nombre}')
        #pdb.set_trace()

