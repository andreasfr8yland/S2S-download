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
from scripts.Andreas.functions import ACC_anom, ACC_anom_error, ACC_anom_nan
from S2S.xarray_helpers import assign_validation_time

latitude=60
longitude=3
lead_time=21 #Options: 7,14,21,28,35,42
# For Monthly mean
lead_time_start=7
lead_time_end=21
# For Running_Spatial_Mean
mean_grid=4 # Number of grid points to take mean of in both lon and lat direction
start_slice='1999-07-01'
end_slice='2019-06-25'
# start_slice='1900-07-01'
# end_slice='2025-06-25'
SEASONAL=False
Season='DJF'
BIAS_TORABLA=False
Monthly_mean=False
Running_Spatial_Mean=False


hindcast_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_anomalies_concat.nc')
era_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_anomalies_concat.nc')


# hindcast_anomalies = hindcast_anomalies.sortby('time')
# era_anomalies = era_anomalies.sortby('time')

# hindcast_anomalies.sortby(hindcast_anomalies.time)
# era_anomalies.sortby(hindcast_anomalies.time)

# hindcast_anomalies = assign_validation_time(hindcast_anomalies)
# era_anomalies = assign_validation_time(hindcast_anomalies)

hindcast_anomaly = hindcast_anomalies.U10
era_anomaly = era_anomalies.U10

# hindcast_anomaly = hindcast_anomaly.rolling(time=4, center=True).mean()
# era_anomaly = era_anomaly.rolling(time=4, center=True).mean()


if SEASONAL==True:
	hindcast_anomaly = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True)
	era_anomaly = era_anomaly.where(era_anomaly.time.dt.season==Season,drop=True)
# hindcast_anomaly.sel(lon=(hindcast_anomaly.lon < 40))
# era_anomaly.sel(lon=(era_anomaly.lon < 40))

# hindcast_anomaly = hindcast_anomaly.assign_coords(lon=(((hindcast_anomaly.lon + 180) % 360) - 180))
# era_anomaly= era_anomaly.assign_coords(lon=(((era_anomaly.lon + 180) % 360) - 180))
#
# hindcast_anomaly = hindcast_anomaly.sortby(hindcast_anomaly.lon)
# era_anomaly = era_anomaly.sortby(hindcast_anomaly.lon)

# hindcast_anomaly.sortby(hindcast_anomaly.time)
if Monthly_mean==False:
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).values
	oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
	# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).time.values
	# validation_time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).validation_time.values
	lat_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lat.values
	lon_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lon.values
	time_oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).time.values
elif Monthly_mean==True:
	ha = hindcast_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).values
	oa = era_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).values
	# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
	time_ha = hindcast_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).time.values
	# validation_time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).validation_time.values
	lat_ha = hindcast_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lat.values
	lon_ha = hindcast_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lon.values
	time_oa = era_anomaly.sel(step=slice(pd.Timedelta(lead_time_start,'D'), pd.Timedelta(lead_time_end,'D'))).mean(dim='step').sel(time=slice(start_slice,end_slice)).time.values

if Running_Spatial_Mean==True:
	weights_ha = np.cos(np.deg2rad(hindcast_anomaly.lat))
	weights_ha.name = "weights"
	weights_oa = np.cos(np.deg2rad(era_anomaly.lat))
	weights_oa.name = "weights"
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).rolling(lat=mean_grid, center=True).mean().rolling(lon=mean_grid, center=True).mean().values
	oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).rolling(lat=4, center=True).mean().rolling(lon=mean_grid, center=True).mean().values
	# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).rolling(lat=mean_grid, center=True).mean().rolling(lon=mean_grid, center=True).mean().time.values
	# validation_time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).validation_time.values
	lat_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).rolling(lat=mean_grid, center=True).mean().rolling(lon=mean_grid, center=True).mean().lat.values
	lon_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).rolling(lat=mean_grid, center=True).mean().rolling(lon=mean_grid, center=True).mean().lon.values
Anom_Corr = ACC_anom_nan(ha,oa)


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
	if Monthly_mean==False:
		if Running_Spatial_Mean==False:
			ax.set_title('10m wind ACC - Leadtime {:01d} days for {}'.format(lead_time, Season), fontsize=14)
		elif Running_Spatial_Mean==True:
			ax.set_title('10m wind ACC - {}x{} grid mean Leadtime {:01d} days for {}'.format(mean_grid, mean_grid, lead_time, Season), fontsize=14)
	if Monthly_mean==True:
		if Running_Spatial_Mean==False:
			ax.set_title('10m wind ACC - Mean leadtime {:01d} to {:01d} days for {}'.format(lead_time_start,lead_time_end,Season), fontsize=14)
		elif Running_Spatial_Mean==True:
			ax.set_title('10m wind ACC - {}x{} grid mean Mean leadtime {:01d} to {:01d} days for {}'.format(mean_grid, mean_grid, lead_time_start,lead_time_end,Season), fontsize=14)
elif SEASONAL==False:
	if Monthly_mean==False:
		if Running_Spatial_Mean==False:
			ax.set_title('10m wind ACC - Leadtime {:01d}'.format(lead_time), fontsize=14)
		elif Running_Spatial_Mean==True:
			ax.set_title('10m wind ACC - {}x{} grid mean Leadtime {:01d}'.format(mean_grid, mean_grid, lead_time), fontsize=14)
	if Monthly_mean==True:
		if Running_Spatial_Mean==False:
			ax.set_title('10m wind ACC - Mean leadtime {:01d} to {:01d}'.format(lead_time_start,lead_time_end), fontsize=14)
		elif Running_Spatial_Mean==True:
			ax.set_title('10m wind ACC - {}x{} grid mean Mean leadtime {:01d} to {:01d}'.format(mean_grid, mean_grid, lead_time_start,lead_time_end), fontsize=14)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.clabel(CS2, inline=1, fontsize=10, fmt='%1.2f')
if SEASONAL==True:
	if Monthly_mean==False:
		if Running_Spatial_Mean==False:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_leadtime_{:01d}_{}.png'.format(lead_time, Season), transparent=True)
		elif Running_Spatial_Mean==True:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_{}x{}_grid_mean_leadtime_{:01d}_{}.png'.format(mean_grid, mean_grid, lead_time, Season), transparent=True)
	elif Monthly_mean==True:
		if Running_Spatial_Mean==False:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_mean_leadtime_{:01d}_to_{:01d}_{}.png'.format(lead_time_start,lead_time_end, Season), transparent=True)
		elif Running_Spatial_Mean==True:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_{}x{}_grid_mean_mean_leadtime_{:01d}_to_{:01d}_{}.png'.format(mean_grid, mean_grid, lead_time_start,lead_time_end, Season), transparent=True)
elif SEASONAL==False:
	if Monthly_mean==False:
		if Running_Spatial_Mean==False:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_leadtime_{:01d}.png'.format(lead_time), transparent=True)
		elif Running_Spatial_Mean==True:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_{}x{}_grid_mean_leadtime_{:01d}.png'.format(mean_grid, mean_grid, lead_time), transparent=True)
	elif Monthly_mean==True:
		if Running_Spatial_Mean==False:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_mean_leadtime_{:01d}_to_{:01d}.png'.format(lead_time_start, lead_time_end), transparent=True)
		elif Running_Spatial_Mean==True:
			fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_{}x{}_grid_mean_mean_leadtime_{:01d}_to_{:01d}.png'.format(lead_time_start, lead_time_end), transparent=True)
plt.show()

# exit()

ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').values
oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').time.values
time_oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).time.values

print(len(time_ha))
print(len(time_oa))

Anom_Corr = ACC_anom_nan(ha,oa)

print(time_ha.shape,ha.shape,time_oa.shape,oa.shape)

plt.figure(figsize=(20,5))
plt.plot_date(time_oa,oa,fmt='-', marker='.', color='DarkOrange', label='Era5 anomaly')
plt.plot_date(time_ha,ha,fmt='-', marker='.', color='MidnightBlue', label='Hindcast anomaly')
plt.plot([], [], ' ', label='ACC = %.2f' %(Anom_Corr))
#plt.plot([], [], ' ', label='Bias = %.2f' %(Bias))
#plt.plot([], [], ' ', label='RMSE = %.2f' %(RMSE))
plt.grid()
#plt.show()
plt.ylabel('10 m wind speed [m s$^{-1}$]')
plt.xlabel('Time')
plt.title('S2S extended range forecast [Lat:%.1f, Lon:%.1f, Lead time: %i days]' %(latitude, longitude, lead_time))
#plt.title('%s climatological anomaly at Lat=%.2f, Lon=%.2f' %(vname, Slat[iy], Slon[ix]))
plt.legend()
#fig.savefig('%s at Lon=%.2f, Lat=%.2f.png' %(vname, Slon[ix], Slat[iy]))
# if Lat_lon_int==False:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.1f_Lon_%.1f_leadtime_%i_days.png' %(latitude,longitude,lead_time))
# if Lat_lon_int==True:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.2f_Lon_%.0f_leadtime_%i_days_5x5_mean.png' %(np.mean(latitude_int),np.mean(longitude_int),lead_time))
plt.show()

# Apply bootstraps
n = len(time_ha)
reps = 1000
ha_bootstrap = np.random.choice(ha, (n, reps))
# oa_bootstrap = np.random.choice(oa, (n, reps))

ACC_bootstrap = np.zeros(reps, dtype=float)
print(len(ha_bootstrap[:,0]))
print(len(ha_bootstrap[0,:]))
for i in range(reps):
	ACC_bootstrap[i] = ACC_anom(ha_bootstrap[:,i], oa)

ACC_5Perc = np.percentile(ACC_bootstrap, [2.5])
ACC_95Perc = np.percentile(ACC_bootstrap, [97.5])

print(ACC_5Perc)
print(ACC_95Perc)
