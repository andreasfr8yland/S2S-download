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
lead_time=7 #Options: 7,14,21,28,35,42
start_slice='1999-07-01'
end_slice='2019-06-25'
start_slice_lat=59
end_slice_lat=62
start_slice_lon=1
end_slice_lon=4
# start_slice='1900-07-01'
# end_slice='2025-06-25'
SEASONAL=False
Season='DJF'
BIAS_TORABLA=False
Weighted=True

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
ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).values
oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).values
time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).time.values
# validation_time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).validation_time.values
lat_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lat.values
lon_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).mean(dim='member', skipna=True).lon.values
time_oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice)).time.values

# Anom_Corr = ACC_anom_nan(ha,oa)
#
#
# fig, ax = plt.subplots(figsize=(8,8))
# ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=10))
# lvs = np.arange(-0.1, 0.8, 0.1)
# lvs_lines = np.arange(-1, 1, 0.1)
# CS = plt.contourf(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs, extend='both', transform=ccrs.PlateCarree())
# CS2 = plt.contour(lon_ha, lat_ha, np.transpose(Anom_Corr), levels=lvs_lines, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.5, linestyles='--')
# plt.colorbar(CS)
# #plt.colorbar(CS, fraction=0.035, pad=0.04)
# plt.clabel(CS2,fmt='%1.2f')
# ax.coastlines(resolution='50m')
# ax.gridlines()
# ax.add_feature(cartopy.feature.BORDERS)
# if SEASONAL==True:
# 	ax.set_title('Anomaly Correlation - Leadtime {:01d} days for {}'.format(lead_time, Season), fontsize=14)
# elif SEASONAL==False:
# 	ax.set_title('Anomaly Correlation - Leadtime {:01d}'.format(lead_time), fontsize=14)
# ax.set_xlabel('Longitude')
# ax.set_ylabel('Latitude')
# ax.clabel(CS2, inline=1, fontsize=10, fmt='%1.2f')
# if SEASONAL==True:
# 	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_leadtime_{:01d}_{}.png'.format(lead_time, Season), transparent=True)
# elif SEASONAL==False:
# 	fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/NorthAtlantic/ACC_map_NA_leadtime_{:01d}.png'.format(lead_time), transparent=True)
# plt.show()
#
# exit()

if Weighted==True:
	weights = np.cos(np.deg2rad(hindcast_anomaly.lat))
	weights.name = "weights"
	# hindcast_anomaly = hindcast_anomaly.weighted(weights)
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').weighted(weights).mean(dim='lat').mean(dim='lon')
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').weighted(weights).mean(dim='lat').mean(dim='lon').time.values
	weights = np.cos(np.deg2rad(era_anomaly.lat))
	weights.name = "weights"
	# era_anomaly = era_anomaly.weighted(weights)
	oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).weighted(weights).mean(dim='lat').mean(dim='lon')
	time_oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).weighted(weights).mean(dim='lat').mean(dim='lon').time.values

if Weighted==False:
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon')
	oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon')
	# climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').time.values
	time_oa = era_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').time.values

# print(len(time_ha))
# print(len(time_oa))

Anom_Corr = ACC_anom(ha.values,oa.values)
Anom_Corr_DJF = ACC_anom(ha.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).values, oa.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).values)
Anom_Corr_MAM = ACC_anom(ha.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).values, oa.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).values)
Anom_Corr_JJA = ACC_anom(ha.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).values, oa.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).values)
Anom_Corr_SON = ACC_anom(ha.where(hindcast_anomaly.time.dt.season=='SON',drop=True).values, oa.where(hindcast_anomaly.time.dt.season=='SON',drop=True).values)

print(time_ha.shape,ha.shape,time_oa.shape,oa.shape)

plot_latitude = np.mean(hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lon').lat.values)
plot_longitude = np.mean(hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').lon.values)

plt.figure(figsize=(20,5))
plt.plot_date(time_oa,oa.values,fmt='-', marker='.', color='DarkOrange', label='Era5 anomaly')
plt.plot_date(time_ha,ha.values,fmt='-', marker='.', color='MidnightBlue', label='Hindcast anomaly')
plt.plot([], [], ' ', label='ACC = %.2f' %(Anom_Corr))
#plt.plot([], [], ' ', label='Bias = %.2f' %(Bias))
#plt.plot([], [], ' ', label='RMSE = %.2f' %(RMSE))
plt.grid()
#plt.show()
plt.ylabel('10 m wind speed [m s$^{-1}$]')
plt.xlabel('Time')
plt.title('S2S extended range forecast spatial mean {Lon: [%.1f, %.1f] Lat: [%.1f, %.1f] Lead time: %i days}' %(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat, lead_time))
#plt.title('%s climatological anomaly at Lat=%.2f, Lon=%.2f' %(vname, Slat[iy], Slon[ix]))
plt.legend()
#fig.savefig('%s at Lon=%.2f, Lat=%.2f.png' %(vname, Slon[ix], Slat[iy]))
# if Lat_lon_int==False:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.1f_Lon_%.1f_leadtime_%i_days.png' %(latitude,longitude,lead_time))
# if Lat_lon_int==True:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.2f_Lon_%.0f_leadtime_%i_days_5x5_mean.png' %(np.mean(latitude_int),np.mean(longitude_int),lead_time))
plt.show()

leadtime_series = [7, 14, 21, 28, 35]

ha_series={}
oa_series={}
acc_series=np.zeros_like(leadtime_series, dtype=float)

for idx, x in enumerate(leadtime_series):
	# print(x)
	# print(idx)
	# print(leadtime_series[idx])
	ha_series[idx] = hindcast_anomaly.sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').values
	oa_series[idx] = era_anomaly.sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').values
	acc_series[idx] = ACC_anom(ha_series[idx], oa_series[idx])

plt.figure()
plt.plot(leadtime_series, acc_series, marker='o')
plt.title('ACC for spatial mean wind within Lon [%.1f, %.1f] Lat [%.1f, %.1f]' %(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat))
plt.xlabel('Lead time [days]')
plt.ylabel('Anomaly correlation')
plt.show()

# Seasonal version
leadtime_series = [7, 14, 21, 28, 35]

ha_series_DJF={}
ha_series_MAM={}
ha_series_JJA={}
ha_series_SON={}

oa_series_DJF={}
oa_series_MAM={}
oa_series_JJA={}
oa_series_SON={}

acc_series_DJF=np.zeros_like(leadtime_series, dtype=float)
acc_series_MAM=np.zeros_like(leadtime_series, dtype=float)
acc_series_JJA=np.zeros_like(leadtime_series, dtype=float)
acc_series_SON=np.zeros_like(leadtime_series, dtype=float)

for idx, x in enumerate(leadtime_series):
	# print(x)
	# print(idx)
	# print(leadtime_series[idx])
	ha_series_DJF[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').values
	ha_series_MAM[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').values
	ha_series_JJA[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').values
	ha_series_SON[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='SON',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').mean(dim='lat').mean(dim='lon').values
	oa_series_DJF[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').values
	oa_series_MAM[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').values
	oa_series_JJA[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').values
	oa_series_SON[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='SON',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).mean(dim='lat').mean(dim='lon').values
	acc_series_DJF[idx] = ACC_anom(ha_series_DJF[idx], oa_series_DJF[idx])
	acc_series_MAM[idx] = ACC_anom(ha_series_MAM[idx], oa_series_MAM[idx])
	acc_series_JJA[idx] = ACC_anom(ha_series_JJA[idx], oa_series_JJA[idx])
	acc_series_SON[idx] = ACC_anom(ha_series_SON[idx], oa_series_SON[idx])



plt.figure()
plt.plot(leadtime_series, acc_series, marker='o', label='Full year')
plt.plot(leadtime_series, acc_series_DJF, marker='o', label='DJF')
plt.plot(leadtime_series, acc_series_MAM, marker='o', label='MAM')
plt.plot(leadtime_series, acc_series_JJA, marker='o', label='JJA')
plt.plot(leadtime_series, acc_series_SON, marker='o', label='SON')
plt.title('ACC for spatial mean wind within Lon [%.1f, %.1f] Lat [%.1f, %.1f]' %(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat))
plt.xlabel('Lead time [days]')
plt.ylabel('Anomaly correlation')
plt.legend()
plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_within_Lon_{:01d}_to_{:01d}_Lat_{:01d}_to_{:01d}.png'.format(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat), transparent=True, dpi=350)
plt.show()

ha_series_max={}
oa_series_max={}
acc_series_max=np.zeros_like(leadtime_series, dtype=float)

for idx, x in enumerate(leadtime_series):
	# print(x)
	# print(idx)
	# print(leadtime_series[idx])
	ha_series_max[idx] = hindcast_anomaly.sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').max(dim='lat').max(dim='lon').values
	oa_series_max[idx] = era_anomaly.sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).max(dim='lat').max(dim='lon').values
	acc_series_max[idx] = ACC_anom(ha_series_max[idx], oa_series_max[idx])

plt.figure()
plt.plot(leadtime_series, acc_series_max, marker='o')
plt.title('ACC for spatial max wind within Lon [%.1f, %.1f] Lat [%.1f, %.1f]' %(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat))
plt.xlabel('Lead time [days]')
plt.ylabel('Anomaly correlation')
plt.show()
# exit()


# Seasonal version
leadtime_series = [7, 14, 21, 28, 35]

ha_series_DJF_max={}
ha_series_MAM_max={}
ha_series_JJA_max={}
ha_series_SON_max={}

oa_series_DJF_max={}
oa_series_MAM_max={}
oa_series_JJA_max={}
oa_series_SON_max={}

acc_series_DJF_max=np.zeros_like(leadtime_series, dtype=float)
acc_series_MAM_max=np.zeros_like(leadtime_series, dtype=float)
acc_series_JJA_max=np.zeros_like(leadtime_series, dtype=float)
acc_series_SON_max=np.zeros_like(leadtime_series, dtype=float)

for idx, x in enumerate(leadtime_series):
	# print(x)
	# print(idx)
	# print(leadtime_series[idx])
	ha_series_DJF_max[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').max(dim='lat').max(dim='lon').values
	ha_series_MAM_max[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').max(dim='lat').max(dim='lon').values
	ha_series_JJA_max[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').max(dim='lat').max(dim='lon').values
	ha_series_SON_max[idx] = hindcast_anomaly.where(hindcast_anomaly.time.dt.season=='SON',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice),lat=slice(start_slice_lat,end_slice_lat),lon=slice(start_slice_lon, end_slice_lon)).mean(dim='member').max(dim='lat').max(dim='lon').values
	oa_series_DJF_max[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='DJF',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).max(dim='lat').max(dim='lon').values
	oa_series_MAM_max[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='MAM',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).max(dim='lat').max(dim='lon').values
	oa_series_JJA_max[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='JJA',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).max(dim='lat').max(dim='lon').values
	oa_series_SON_max[idx] = era_anomaly.where(hindcast_anomaly.time.dt.season=='SON',drop=True).sel(step=pd.Timedelta(x,'D'),time=slice(start_slice,end_slice), lat=slice(start_slice_lat,end_slice_lat), lon=slice(start_slice_lon, end_slice_lon)).max(dim='lat').max(dim='lon').values
	acc_series_DJF_max[idx] = ACC_anom(ha_series_DJF_max[idx], oa_series_DJF_max[idx])
	acc_series_MAM_max[idx] = ACC_anom(ha_series_MAM_max[idx], oa_series_MAM_max[idx])
	acc_series_JJA_max[idx] = ACC_anom(ha_series_JJA_max[idx], oa_series_JJA_max[idx])
	acc_series_SON_max[idx] = ACC_anom(ha_series_SON_max[idx], oa_series_SON_max[idx])



plt.figure()
plt.plot(leadtime_series, acc_series_max, marker='o', label='Full year')
plt.plot(leadtime_series, acc_series_DJF_max, marker='o', label='DJF')
plt.plot(leadtime_series, acc_series_MAM_max, marker='o', label='MAM')
plt.plot(leadtime_series, acc_series_JJA_max, marker='o', label='JJA')
plt.plot(leadtime_series, acc_series_SON_max, marker='o', label='SON')
plt.title('ACC for spatial maximum wind within Lon [%.1f, %.1f] Lat [%.1f, %.1f]' %(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat))
plt.xlabel('Lead time [days]')
plt.ylabel('Anomaly correlation')
plt.legend()
plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_max_wind_within_Lon_{:01d}_to_{:01d}_Lat_{:01d}_to_{:01d}.png'.format(start_slice_lon, end_slice_lon, start_slice_lat, end_slice_lat), transparent=True, dpi=350)
plt.show()
