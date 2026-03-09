import xarray as xr
import numpy as np
import pdb

### Modify the sector
path = '/home/netapp-clima/scratch/acasallas/MPI/TOTALS/'
pathout = '/home/tompkins-archive2/acasallas/Reg_Col_WRF_chem/Results/TOTALS/'

filenum = np.arange(1,100)

for filen in filenum:
    try:
        if filen < 10:  
            ds = xr.open_dataset(pathout+'PM2_5_file0'+str(filen)+'.nc')
        else:
            ds = xr.open_dataset(pathout+'PM2_5_file'+str(filen)+'.nc')
    except:
        print('File number: '+str(filen))
        if filen < 10:
            ds = xr.open_dataset(path+'wrfout_file0'+str(filen)+'.nc')
        else:
            ds = xr.open_dataset(path+'wrfout_file'+str(filen)+'.nc')
        pm = ds.PM2_5_DRY[:,0,:,:]
        ozone = ds.o3[:,0,:,:]
        del(ds)
        if filen < 10:
            pm.to_netcdf(pathout+'PM2_5_file0'+str(filen)+'.nc')
            ozone.to_netcdf(pathout+'Ozone_file0'+str(filen)+'.nc')
        else:
            pm.to_netcdf(pathout+'PM2_5_file'+str(filen)+'.nc')
            ozone.to_netcdf(pathout+'Ozone_file'+str(filen)+'.nc')
