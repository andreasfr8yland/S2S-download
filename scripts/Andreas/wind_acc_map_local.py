import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import datetime
import cartopy
import cartopy.crs as ccrs
import xskillscore as xs

import S2S.graphics.graphics as gr
import S2S.xarray_helpers as xh
from S2S.local_configuration import config
#from scripts.Henrik.acc import ACC
from S2S.models import bias_adjustment_torralba
import S2S.scoring as sc
from scripts.Andreas.functions import ACC_anom, ACC_anom_error

latitude=57
longitude=3
lead_time=7 #Options: 7,14,21,28,35,42
start_slice='1999-07-01'
end_slice='2019-06-25'
SEASONAL=True
Season='DJF'
BIAS_TORABLA=True
# RollingMean=True

hindcast_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
era_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
era_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
era_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')

hindcast_anomaly = hindcast_anomalies.abs_wind
observation_mean = era_mean.abs_wind
observation_std = era_std.abs_wind
observation_anomaly = era_anomalies.abs_wind
if BIAS_TORABLA==True:
	hindcast_anomaly_toralba = bias_adjustment_torralba(hindcast_anomaly.sel(lat=latitude, lon=longitude), observation_anomaly.sel(lat=latitude, lon=longitude))

#hindcast_anomaly = hindcast_anomaly*observation_std+observation_mean
#observation_anomaly = observation_anomaly*observation_std+observation_mean
if SEASONAL==True:
	hindcast_anomaly = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True)
	observation_mean = observation_mean.where(hindcast_anomaly.time.dt.season==Season,drop=True)
	observation_std = observation_std.where(hindcast_anomaly.time.dt.season==Season,drop=True)
	observation_anomaly = observation_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True)

#if RollingMean==True:
#	hindcast_anomaly = hindcast_anomaly.rolling('lat',center=True).mean()
#	#.rolling(lon,center=True,window=4).mean()
#	observation_anomaly = observation_anomaly.rolling('lat',center=True).mean()
#	#.rolling(lon,center=True,window=4).mean()

ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').values
oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').time.values
lat_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').lat.values
lon_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').lon.values
time_oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).time.values


# weights = np.cos(np.deg2rad(hindcast_anomaly.lat))
# weights.name = "weights"
#
#
# hindcast_anomaly_weighted = hindcast_anomaly.weighted(weights)
# observation_anomaly_weighted = observation_anomaly.weighted(weights)
#
# hindcast_anomaly_weighted_mean = hindcast_anomaly_weighted.mean(("lon", "lat"))
# observation_anomaly_weighted_mean = observation_anomaly_weighted.mean(("lon", "lat"))
#
# wmh = hindcast_anomaly_weighted_mean.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').values
# wmo =observation_anomaly_weighted_mean.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
#
# Anom_Corr_weighted = ACC_anom(wmh,wmo)
#
# print(Anom_Corr_weighted)

#Anom_Corr = np.empty_like(ha[0,:,:])
#for y in range(len(lat_ha)):
#	for x in range(len(lon_ha)):
#		Anom_Corr[x,y] = ACC_anom(ha[:,x,y],oa[:,x,y])

Anom_Corr = ACC_anom(ha,oa)

fig, ax = plt.subplots(figsize=(8,8))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=10))
lvs = np.arange(-0.1, 0.8, 0.1)
lvs_lines = np.arange(-1, 1, 0.1)
CS = plt.contourf(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs, extend='both', transform=ccrs.PlateCarree())
CS2 = plt.contour(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs_lines, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.5, linestyles='--')
plt.colorbar(CS)
#plt.colorbar(CS, fraction=0.035, pad=0.04)
plt.clabel(CS2,fmt='%1.2f')
ax.coastlines(resolution='50m')
ax.gridlines()
ax.add_feature(cartopy.feature.BORDERS)
if SEASONAL==True:
	ax.set_title('Anomaly Correlation - Leadtime {:01d} days for {}'.format(lead_time, Season), fontsize=14)
elif SEASONAL==False:
	ax.set_title('Anomaly Correlation - Leadtime {:01d}'.format(lead_time), fontsize=14)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.clabel(CS2, inline=1, fontsize=10, fmt='%1.2f')
if SEASONAL==True:
	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/ACC_map_leadtime_{:01d}_{}.png'.format(lead_time, Season), transparent=True)
elif SEASONAL==False:
	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/ACC_map_leadtime_{:01d}.png'.format(lead_time), transparent=True)
plt.show()

Spatial_Anomaly_Correlation = sc.ACc(hindcast_anomaly,observation_anomaly,centered=False)

print(Spatial_Anomaly_Correlation.where(hindcast_anomaly.time.dt.season==Season,drop=True).mean('time'))
print(Spatial_Anomaly_Correlation.step)

plt.figure()
plt.plot(
	Spatial_Anomaly_Correlation.mean('time').step.values
	,Spatial_Anomaly_Correlation.mean('time').values
	, marker='.', color='MidnightBlue', label='Spatial ACC')
# plt.plot(
# 	[7,14,21,28,35,42]
# 	,Spatial_Anomaly_Correlation.mean('time').values
# 	,fmt='-', marker='.', color='MidnightBlue', label='Spatial ACC')


plt.show()

#print(observation_std)
#print(xd1)

#hindcast_anomaly_vd = xh.assign_validation_time(hindcast_anomaly)
hindcast_anomaly_vd = hindcast_anomaly.expand_dims('validation_time').assign_coords(validation_time=hindcast_anomaly.time+hindcast_anomaly.step).drop('time')
observation_anomaly_vd = observation_anomaly.expand_dims('validation_time').assign_coords(validation_time=observation_anomaly.time+observation_anomaly.step).drop('time')

print(hindcast_anomaly)
print(observation_anomaly)


#Anom_Corr_Henrik = ACC(hindcast_anomaly_vd.transpose("member","step","validation_time","lon","lat").values,observation_anomaly_vd.transpose("step","validation_time","lon","lat").values,np.ones_like(hindcast_anomaly_vd.values))

#gr.skill_map(
#	hindcast_anomaly_vd,
#	observation_anomaly_vd,
#	dim='validation_time.month',
#	title='Skill map 10 m wind',
#	lead_time=[7,14,21,28,35],
#	filename='skill_map_wind'
#	)
