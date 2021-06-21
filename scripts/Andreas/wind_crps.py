import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import datetime
import xskillscore as xs


from S2S.local_configuration import config
from S2S.data_handler        import BarentsWatch, ERA5, ECMWF_S2SH, Archive

import S2S.xarray_helpers    as xh
import S2S.models            as models
import S2S.graphics.graphics as gr

from S2S.graphics import mae,crps,brier,spread_error
from S2S import location_cluster

latitude=60
longitude=4.5
lead_time=21 #Options: 7,14,21,28,35,42
start_slice='2001-01-01'
end_slice='2019-01-01'
#start_slice='2018-01-01'
#end_slice='2019-01-01'
season_cutoff=False
Season='DJF'
standard_bias_correction=True
nird=True

if nird==False:
	xd1 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
	xd1_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
	xd1_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
	xd2 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')
	xd3 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absoulte_wind_era.nc')
elif nird==True:
	xd1 = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_anomalies_hindcast.nc')
	xd1_mean = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_era_mean.nc')
	xd1_std = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_era_std.nc')
	xd2 = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_anomalies_era.nc')
	xd3 = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absoulte_wind_era.nc')

hindcast_anomaly = xd1.abs_wind
observation_mean = xd1_mean.abs_wind
observation = xd3.abs_wind
observation_std = xd1_std.abs_wind
observation_anomaly = xd2.abs_wind
if standard_bias_correction==True:
	hindcast_sbc = hindcast_anomaly*observation_std+observation_mean # sbc = standard_bias_correction
	observation_sbc = observation_anomaly*observation_std+observation_mean

 crps.skill_agg(
 observation,
 hindcast_sbc,
 observation_mean,
 observation_std,
 title='ERA Simple Bias Adjustment',
 filename='crps_ERA_sb',
 dim='validation_time.month'
 ) 

