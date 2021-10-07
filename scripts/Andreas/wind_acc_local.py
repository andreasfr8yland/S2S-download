import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import datetime

latitude=60
longitude=3
#latitude_int=[57,58.5,60,61.5,63]
latitude_int=np.linspace(latitude-3,latitude+3,num=5)
longitude_int=np.linspace(longitude-3,longitude+3,num=5)
#longitude_int=[1.5,3,4.5,6,7.5]
lead_time=21 #Options: 7,14,21,28,35,42
start_slice='2001-01-01'
end_slice='2019-01-01'
#start_slice='2018-01-01'
#end_slice='2019-01-01'
season_cutoff=False
# Season cutoff is not included in titles and savefig at the moment!!
Season='DJF'
standard_bias_correction=False
nird=False
Lat_lon_int=False
BIAS_TORABLA=False
single_member=True
member_nr=0

if nird==False:
	xd1 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
	xd1_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
	xd1_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
	xd2 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')
elif nird==True:
	xd1 = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_anomalies_hindcast.nc')
	xd1_mean = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_era_mean.nc')
	xd1_std = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_era_std.nc')
	xd2 = xr.open_dataset('/nird/projects/NS9853K/DATA/processed/wind_S2S/wind10m/absolute_wind_anomalies_era.nc')


hindcast_anomaly = xd1.abs_wind
observation_mean = xd1_mean.abs_wind
observation_std = xd1_std.abs_wind
observation_anomaly = xd2.abs_wind
if standard_bias_correction==True:
	hindcast_anomaly = hindcast_anomaly*observation_std+observation_mean
	observation_anomaly = observation_anomaly*observation_std+observation_mean

if BIAS_TORABLA==True:
	hindcast_anomaly_toralba = bias_adjustment_torralba(hindcast_anomaly.sel(lat=latitude, lon=longitude), observation_anomaly.sel(lat=latitude, lon=longitude))


if season_cutoff==True:
	ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').values
	oa = observation_anomaly.where(observation_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	climatology= observation_mean.where(observation_mean.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	time_ha = hindcast_anomaly.where(hindcast_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').time.values
	time_oa = observation_anomaly.where(observation_anomaly.time.dt.season==Season,drop=True).sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).time.values
elif season_cutoff==False:
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').values
	oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).mean(dim='member').time.values
	time_oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).time.values
	if Lat_lon_int==True:
		ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int).mean(dim='member')
		weights = np.cos(np.deg2rad(ha.lat))
		weights.name = "weights"
		ha_weighted = ha.weighted(weights)
		ha = ha_weighted.mean(("lon", "lat")).values
		oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
		weights = np.cos(np.deg2rad(oa.lat))
		weights.name = "weights"
		oa_weighted = oa.weighted(weights)
		oa = oa_weighted.mean(("lon", "lat")).values
		climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
		climatology = climatology.mean(dim='lat').mean(dim='lon').values
		time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int).mean(dim='member')
		time_ha = time_ha.mean(dim='lat').mean(dim='lon').time.values
		time_oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
		time_oa = time_oa.mean(dim='lat').mean(dim='lon').time.values

if single_member==True:
	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude, member=member_nr).values
	oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).values
	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude, member=member_nr).time.values
	time_oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude, lon=longitude).time.values
	# if Lat_lon_int==True:
	# 	ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int).mean(dim='member')
	# 	weights = np.cos(np.deg2rad(ha.lat))
	# 	weights.name = "weights"
	# 	ha_weighted = ha.weighted(weights)
	# 	ha = ha_weighted.mean(("lon", "lat")).values
	# 	oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
	# 	weights = np.cos(np.deg2rad(oa.lat))
	# 	weights.name = "weights"
	# 	oa_weighted = oa.weighted(weights)
	# 	oa = oa_weighted.mean(("lon", "lat")).values
	# 	climatology= observation_mean.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
	# 	climatology = climatology.mean(dim='lat').mean(dim='lon').values
	# 	time_ha = hindcast_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int).mean(dim='member')
	# 	time_ha = time_ha.mean(dim='lat').mean(dim='lon').time.values
	# 	time_oa = observation_anomaly.sel(step=pd.Timedelta(lead_time,'D'),time=slice(start_slice,end_slice), lat=latitude_int, lon=longitude_int)
	# 	time_oa = time_oa.mean(dim='lat').mean(dim='lon').time.values



#def ACC_anom(FC_anom,OBS_anom):
#	top=np.sum(FC_anom*OBS_anom)
#	bottom=np.sqrt(np.sum(np.square(FC_anom))*np.sum(np.square(OBS_anom)))
#	ACC_anom=np.divide(top,bottom)
#	return ACC_anom

# ACC with subtracted error
def ACC_anom(FC_anom,OBS_anom):
	top=np.mean((FC_anom-np.mean(FC_anom))*(OBS_anom-np.mean(OBS_anom)))
	bottom=np.sqrt(np.mean(np.square(FC_anom-np.mean(FC_anom)))*np.mean(np.square(OBS_anom-np.mean(OBS_anom))))
	ACC_anom=np.divide(top,bottom)
	return ACC_anom

#def ACC_anom(FC_anom,OBS_anom):
#	top = np.mean(FC_anom*OBS_anom)
#	bottom = np.sqrt(np.mean(FC_anom**2)*np.mean(OBS_anom**2))
#	ACC = top/bottom
#	return ACC_anom
def ACC(FC,OBS,CL):
	top = np.nanmean((FC-CL)*(OBS-CL))
	bottom = np.sqrt(np.nanmean((FC-CL)**2)*np.nanmean((OBS-CL)**2))
	ACC = top/bottom
	return ACC

if standard_bias_correction==True:
	Anomaly_Correlation = ACC(ha,oa,climatology)
elif standard_bias_correction==False:
	Anomaly_Correlation = ACC_anom(ha,oa)
#Anomaly_Correlation = ACC_anom(ha,oa)
#Anomaly_Correlation = np.correlate(ha, oa)
Error = np.subtract(oa, ha)
Bias = np.mean(Error)
RMSE = np.sqrt(Bias**2)

plt.figure(figsize=(20,5))
plt.plot_date(time_ha,ha,fmt='-', marker='.', color='MidnightBlue', label='Hindcast anomaly')
plt.plot_date(time_oa,oa,fmt='-', marker='.', color='DarkOrange', label='Era5 anomaly')
plt.plot([], [], ' ', label='ACC = %.2f' %(Anomaly_Correlation))
#plt.plot([], [], ' ', label='Bias = %.2f' %(Bias))
#plt.plot([], [], ' ', label='RMSE = %.2f' %(RMSE))
plt.grid()
#plt.show()
plt.ylabel('10 m wind speed [m s$^{-1}$]')
plt.xlabel('Time')
if Lat_lon_int==False:
	plt.title('S2S extended range forecast [Lat:%.1f, Lon:%.1f, Lead time: %i days]' %(latitude, longitude, lead_time))
elif Lat_lon_int==True:
	plt.title('S2S extended range forecast [Lat:%.1f, Lon:%.1f, Lead time: %i days, 5x5 grid mean]' %(np.mean(latitude_int),np.mean(longitude_int),lead_time))
#plt.title('%s climatological anomaly at Lat=%.2f, Lon=%.2f' %(vname, Slat[iy], Slon[ix]))
plt.legend()
#fig.savefig('%s at Lon=%.2f, Lat=%.2f.png' %(vname, Slon[ix], Slat[iy]))
# if Lat_lon_int==False:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.1f_Lon_%.1f_leadtime_%i_days.png' %(latitude,longitude,lead_time))
# if Lat_lon_int==True:
# 	plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.2f_Lon_%.0f_leadtime_%i_days_5x5_mean.png' %(np.mean(latitude_int),np.mean(longitude_int),lead_time))
plt.show()
#plt.savefig('Tidsserie.pdf'y
