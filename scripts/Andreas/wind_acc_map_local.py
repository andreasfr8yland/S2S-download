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
from scripts.Henrik.acc import ACC


latitude=57
longitude=3
lead_time=28 #Options: 7,14,21,28,35,42
start_slice='1999-07-01'
end_slice='2019-06-25'
Season='DJF'

xd1 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
xd1_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
xd1_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
xd2 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')


hindcast_anomaly = xd1.abs_wind
observation_mean = xd1_mean.abs_wind
observation_std = xd1_std.abs_wind
observation_anomaly = xd2.abs_wind
#hindcast_anomaly = hindcast_anomaly*observation_std+observation_mean
#observation_anomaly = observation_anomaly*observation_std+observation_mean

ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').values
oa = observation_anomaly.where(observation_anomaly.time.dt.season==Season,drop=True) .sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
climatology= observation_mean.where(observation_mean.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
time_ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').time.values
lat_ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').lat.values
lon_ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member').lon.values
time_oa = observation_anomaly.where(observation_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).time.values



#def ACC_anom(FC_anom,OBS_anom):
#	top=np.sum(FC_anom*OBS_anom)
#	bottom=np.sqrt(np.sum(np.square(FC_anom))*np.sum(np.square(OBS_anom)))
#	ACC_anom=np.divide(top,bottom)
#	return ACC_anom

#def ACC_anom(FC_anom,OBS_anom):
#	top = np.mean((FC_anom)*(OBS_anom))
#	bottom = np.sqrt(np.mean((FC_anom)**2)*np.mean((OBS_anom)**2))
#	ACC = top/bottom
#	return ACC_anom

def ACC_anom(FC_anom,OBS_anom):
	top=np.sum(FC_anom*OBS_anom)
	bottom=np.sqrt(np.sum(np.square(FC_anom))*np.sum(np.square(OBS_anom)))
	ACC_anom=np.divide(top,bottom)
	return ACC_anom

#def ACC(FC,OBS,CL):
#	top = np.mean((FC-CL)*(OBS-CL))
#	bottom = np.sqrt(np.mean((FC-CL)**2)*np.mean((OBS-CL)**2))
#	ACC = top/bottom
#	return ACC

#def ACC(FC,OBS,CL):
#	top = np.mean((FC-CL)*(OBS-CL))
#	bottom = np.sqrt(np.mean((FC-CL)**2)*np.mean((OBS-CL)**2))
#	ACC = top/bottom
#	return ACC

Anom_Corr = np.empty_like(ha[0,:,:])
for y in range(len(lat_ha)):
	for x in range(len(lon_ha)):
		Anom_Corr[x,y] = ACC_anom(ha[:,x,y],oa[:,x,y])

#Anomaly_Correlation = ACC(ha,oa,climatology)
##Anomaly_Correlation = ACC_anom(ha,oa)
##Anomaly_Correlation = np.correlate(ha, oa)
#Error = np.subtract(oa, ha)
#Bias = np.mean(Error)
#RMSE = np.sqrt(Bias)

fig, ax = plt.subplots(figsize=(12, 8))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
#ax = plt.axes(projection=ccrs.Mercator())
#ax = plt.axes(projection=ccrs.Mercator(central_longitude=5.0,min_latitude=40,max_latitude=80))
#plt.contourf(Slon, Slat, Anom_Corr, 60, transform=ccrs.PlateCarree())
#plt.contourf(Slat, Slon, Anom_Corr, 60, transform=ccrs.Mercator(central_longitude=5.0,min_latitude=40,max_latitude=80))
lvs = np.arange(-0.1, 0.8, 0.1)
lvs_lines = np.arange(-1, 1, 0.1)
#CS = plt.contourf(lat_ha, lon_ha, Anom_Corr, levels=lvs, extend='both', transform=ccrs.PlateCarree())
#CS2 = plt.contour(lat_ha, lon_ha, Anom_Corr, levels=lvs_lines, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.5, linestyles='--')
CS = plt.contourf(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs, extend='both', transform=ccrs.PlateCarree())
CS2 = plt.contour(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs_lines, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.5, linestyles='--')
#DT = plt.plot([Slon, Slat][Anom_Corr_p_value<=0.05], Anom_Corr_p_value[Anom_Corr_p_value<=0.05], transform=ccrs.PlateCarree())
plt.colorbar(CS, fraction=0.035, pad=0.04)
plt.clabel(CS2,fmt='%1.2f')
ax.coastlines(resolution='50m')
ax.gridlines()
ax.add_feature(cartopy.feature.BORDERS)
#if Monthly == 'yes':
#	#ax.set_title('Anomaly Correlation - Leadtime={:01d} - Month={:02d}'.format(leadtime, WhichMonth), fontsize=14)
#	ax.set_title('Anomaly Correlation - Leadtime={:01d} - Month={}'.format(leadtime, CurrentMonth), fontsize=14)
#elif Monthly == 'no':
#	ax.set_title('Anomaly Correlation - Leadtime={:01d}'.format(leadtime), fontsize=14)
#elif Monthly == 'DJF' or Monthly == 'MAM' or Monthly == 'JJA' or Monthly == 'SON':
#	ax.set_title('Anomaly Correlation - Leadtime={:01d} - Month={}'.format(leadtime, Monthly))
ax.set_title('Anomaly Correlation - Leadtime={:01d} days '.format(lead_time))
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.clabel(CS2, inline=1, fontsize=10, fmt='%1.2f')
#if Monthly == 'yes':
#	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/Correlation/ACC_map_Europe_leadtime_{:01d}_month_{:02d}.png'.format(leadtime, WhichMonth), transparent=True)
#elif Monthly == 'no':
#	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/Correlation/ACC_map_Europe_leadtime_{:01d}.png'.format(leadtime), transparent=True)
#elif Monthly == 'DJF' or Monthly == 'MAM' or Monthly == 'JJA' or Monthly == 'SON':
#	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/Correlation/ACC_map_Europe_leadtime_{:01d}_month_{}.png'.format(leadtime, Monthly),transparent=True)
plt.show()

print(observation_std)
print(xd1)

hindcast_anomaly_vd = xh.assign_validation_time(hindcast_anomaly)
observation_anomaly_vd = xh.assign_validation_time(observation_anomaly)
Anom_Corr_Henrik = ACC(hindcast_anomaly_vd.transpose("member","step","validation_time","lon","lat").values,observation_anomaly_vd.transpose("step","validation_time","lon","lat").values,np.ones_like(hindcast_anomaly.values))

gr.skill_map(
	hindcast_anomaly,
	observation_anomaly,
	dim='validation_time.month',
	title='SS',
	lead_time=[7,14,21,28,35],
	filename='skill_map_wind'
	)
