import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import datetime

latitude=60
longitude=4.5
lead_time=21 #Options: 7,14,21,28,35,42
start_slice='2001-01-01'
end_slice='2019-01-01'
#start_slice='2018-01-01'
#end_slice='2019-01-01'
season_cutoff=False
Season='DJF'
standard_bias_correction=False

xd1 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
xd1_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
xd1_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
xd2 = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')


hindcast_anomaly = xd1.abs_wind
observation_mean = xd1_mean.abs_wind
observation_std = xd1_std.abs_wind
observation_anomaly = xd2.abs_wind
if standard_bias_correction==True:
	hindcast_anomaly = hindcast_anomaly*observation_std+observation_mean
	observation_anomaly = observation_anomaly*observation_std+observation_mean

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



def ACC_anom(FC_anom,OBS_anom):
	top=np.sum(FC_anom*OBS_anom)
	bottom=np.sqrt(np.sum(np.square(FC_anom))*np.sum(np.square(OBS_anom)))
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
RMSE = np.sqrt(Bias)

plt.figure(figsize=(20,5))
plt.plot_date(time_ha,ha,fmt='-', marker='.', color='MidnightBlue', label='Hindcast anomaly')
plt.plot_date(time_oa,oa,fmt='-', marker='.', color='DarkOrange', label='Era5 anomaly')
plt.plot([], [], ' ', label='ACC = %.2f' %(Anomaly_Correlation))
plt.plot([], [], ' ', label='Bias = %.2f' %(Bias))
plt.plot([], [], ' ', label='RMSE = %.2f' %(RMSE))
plt.grid()
#plt.show()
plt.ylabel('10 m wind speed [m s$^{-1}$]')
plt.xlabel('Time')
plt.title('S2S extended range forecast at Lat=%.2f Lon=%.2f Lead time: %.2f days' %(latitude, longitude, lead_time))
#plt.title('%s climatological anomaly at Lat=%.2f, Lon=%.2f' %(vname, Slat[iy], Slon[ix]))
plt.legend()
#fig.savefig('%s at Lon=%.2f, Lat=%.2f.png' %(vname, Slon[ix], Slat[iy]))
plt.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/Tidsserie_S2S_anomaly_ACC_Lat_%.2f_Lon_%.0f_leadtime_%.2f_days.png' %(latitude,longitude,lead_time))
plt.show()
#plt.savefig('Tidsserie.pdf'y
