import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs
import time as time_lib
import cartopy.crs as ccrs

from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr

domainID = 'NVK'
var      = 'u10'

t_start  = (2020,1,23)
t_end    = (2021,1,4)

clim_t_start  = (2000,1,1)
clim_t_end    = (2021,1,4)

#process_hindcast     = True
#process_observations = False
#verify = True

#high_res = False

#print('Process hindcast')
#if process_hindcast:

#    print('\tLoad hindcast')
#    hindcast = ECMWF_S2SH(high_res=high_res)\
#                    .load(var,t_start,t_end,domainID)[var]-272.15

#    print('\tInterpolate hindcast to point locations')
#    hindcast = hindcast.sortby(['time','step'])\
#                    .interpolate_na(
#                                    dim='lon',
#                                    method='nearest',
#                                    fill_value="extrapolate"
#                                )\
#                        .interp(
#                                lon=observations.lon,
#                                lat=observations.lat,
#                                method='linear'
#                            )
#    print('\tApply 7D running mean')
#    hindcast = hindcast\
#                    .rolling(step=7,center=True).mean()\
#                        .dropna('step')

#    print('\tAssign validation time')
#    hindcast = xh.assign_validation_time(hindcast)

#    hindcast.to_netcdf(config['VALID_DB']+'/h_temp.nc')

#hindcast = xr.open_dataset(config['VALID_DB']+'/h_temp.nc')



air = xr.tutorial.open_dataset("air_temperature").air

p = air.isel(time=0).plot(subplot_kws=dict(projection=ccrs.NorthPolarStereo(central_longitude=0)), transform=ccrs.PlateCarree(),)


p.axes.set_global()

p.axes.coastlines()

