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

latitude=60
longitude=3.0
lead_time=14
Select_Loc=True

# # Norwegian Coast
# hindcast_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_hindcast.nc')
# hindcast_absolute = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_hindcast.nc')
# era_absolute = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era.nc')
# era_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_mean.nc')
# era_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_era_std.nc')
# era_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/absolute_wind_anomalies_era.nc')

# North Atlantic
hindcast_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_anomalies_concat.nc')
hindcast_absolute = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_absolute_concat.nc')
era_absolute = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_absolute_concat.nc')
era_mean = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_mean_concat.nc')
era_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_std_concat.nc')
era_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_anomalies_concat.nc')

# hindcast_anomalies = xh.assign_validation_time(hindcast_anomalies)
# hindcast_absolute = xh.assign_validation_time(hindcast_absolute)
# era_absolute = xh.assign_validation_time(era_absolute)
# era_anomalies = xh.assign_validation_time(era_anomalies)
# era_mean = xh.assign_validation_time(era_mean)
# era_std = xh.assign_validation_time(era_std)


# # Norwegian Coast
# hindcast_anomaly = hindcast_anomalies.abs_wind
# hindcast = hindcast_absolute.abs_wind
# era = era_absolute.abs_wind
# observation_mean = era_mean.abs_wind
# observation_std = era_std.abs_wind
# observation_anomaly = era_anomalies.abs_wind

# North Atlantic
hindcast_anomaly = hindcast_anomalies.U10
hindcast = hindcast_absolute.U10
era = era_absolute.U10
era_mean = era_mean.U10
era_std = era_std.U10
era_anomaly = era_anomalies.U10

# Select location
if Select_Loc==True:
    era_anomaly_loc = era_anomaly.sel(lat=latitude, lon=longitude)
    hindcast_anomaly_loc = hindcast_anomaly.sel(lat=latitude, lon=longitude)
    hindcast_loc = hindcast.sel(lat=latitude, lon=longitude)
    era_loc = era.sel(lat=latitude, lon=longitude)
    era_mean_loc = era_mean.sel(lat=latitude, lon=longitude)
    era_std_loc = era_std.sel(lat=latitude, lon=longitude)

CRPS_loc = xs.crps_ensemble(era_loc, hindcast_loc, dim='time')
CRPS_era_loc = xs.crps_gaussian(era_loc, era_mean_loc, era_std_loc, dim='time')
CRPS_2_loc = sc.crps_ensemble(era_loc, hindcast_loc, fair=True)
print(CRPS_loc)
print(CRPS_2_loc.mean(dim='time'))
print(CRPS_era_loc)

Ones = xr.ones_like(CRPS_loc, dtype=float)
Skill_score = Ones - (CRPS_loc / CRPS_era_loc)
# Skill_score_M2 = (CRPS_era_loc)
print(Skill_score)
exit()

# CRPS = sc.crps_ensemble(observation_anomaly, hindcast_anomaly, dim=hindcast_anomaly.time)
CRPS = xs.crps_ensemble(era, hindcast, dim='time')
CRPS_era = xs.crps_gaussian(era, era_mean, era_std, dim='time')
CRPS_2 = sc.crps_ensemble(era, hindcast, fair=True)
# def skill_plot(in_mod,in_clim,dim='validation_time.month',
#                                 filename='',
#                                 title='',
#                                 ylab=''):
#
#         if dim.split('.')[0]=='validation_time':
#             mod = xh.assign_validation_time(mod)
#             clim = xh.assign_validation_time(clim)
#
#         fname    = 'skill/'+filename+'_'+\
#                                 'location'
#         suptitle = title+' '+'location'
#
#         dim_name = dim.split('.')[0]
#         ###########################
#         #### Initialize figure ####
#         ###########################
#         latex.set_style(style='white')
#         subplots = fg(mod,dim)
#         fig,axes = plt.subplots(subplots[0],subplots[1],
#                         figsize=latex.set_size(width='thesis',
#                             subplots=(subplots[0],subplots[1]))
#                         )
#         axes = axes.flatten()
#         ###########################
#
#         x_group = list(mod.groupby(dim))
#         y_group = list(clim.groupby(dim))
#
#         for n,(xlabel,xdata) in enumerate(x_group):
#
#             ylabel,ydata = y_group[n]
#             # this approach is not bulletproof
#             # check manually that right groups go together
#             print('\tgraphics.skill_plot: matched groups ',xlabel,ylabel)
#
#             xdata = xdata.unstack().sortby(['time','step'])
#             ydata = ydata.unstack().sortby(['time','step'])
#
#             xdata,ydata = xr.align(xdata,ydata)
#
#             ss = 1-(xdata/ydata).mean('time',skipna=True)
#
#             zx = xdata.mean('time',skipna=True)
#             zy = ydata.mean('time',skipna=True)
#
#             x = np.array([td.days for td in ss.step.to_pandas()])
#             z = ss.values
#
#             ax = axes[n]
#
#             ax.set_title(month(xlabel))
#             ax.set_xlabel('lead time [D]')
#             ax.set_ylabel(ylab)
#
#             ax.plot(x,z,'o-',alpha=0.4,ms=1,label='ss')
#             ax.plot(x,zx,'o-',alpha=0.4,ms=1,label='fc')
#             ax.plot(x,zy,'o-',alpha=0.4,ms=1,label='clim')
#             ax.legend()
#             # ax.set_ylim((-1,1))
#             ax.plot(
#                     [0, 1],[0.5, 0.5],'k',
#                     transform=ax.transAxes,alpha=0.7,
#                     linewidth=0.6
#                     )
#
#         fig.suptitle(suptitle)
#         save_fig(fig,fname)
# skill_plot(CRPS_2, CRPS_era)

# exit()
# print(CRPS_2)
CRPS_2_mean = CRPS_2.mean(dim='time').mean(dim='lon').mean(dim='lat')
# print(CRPS_2_mean)
# CRPS_loc = CRPS.sel(lat=latitude, lon=longitude, step=pd.Timedelta(lead_time,'D'))
# print(CRPS_loc)
CRPS_mean = CRPS.mean(dim='lon').mean(dim='lat')
# print(CRPS_mean)

CRPS_era_mean = CRPS_era.mean(dim='lon').mean(dim='lat')
# print(CRPS_era_mean)

from S2S.graphics.crps import ss as crps_ss

CRPSS = crps_ss(hindcast, era)

print(CRPSS.sel(lat=latitude,lon=longitude))

# CRPSS = xr.ones_like(CRPS_2_mean) - (CRPS_2_mean/CRPS_era_mean)
# # CRPSS_numpy = np.subtract(1, np.divide(CRPS_2_mean.values, CRPS_era_mean.values))
# print(CRPSS)

# plt.figure()
# plt.plot(CRPS_mean.step.values, CRPSS)
# plt.title('Atempt CRPSS')
# plt.show()

# CRPSS.plot.line(x='step', color='purple', marker='o')
#
CRPS_2_mean.plot.line(x='step', color='purple', marker='o')
CRPS_era_mean.plot.line(x='step', color='orange', marker='o')
plt.title('Mean CRPS for all lats and lons')
plt.show()
#
# CRPS_era_mean.plot.line(x='step', color='red', marker='o')
# plt.show()
#
# print(CRPS_mean)
# print(CRPS_mean.step.values)

# plt.figure()
# # plt.plot(pd.to_datetime(CRPS_mean.step).values, CRPS_mean.values)
# plt.plot(CRPS_mean.step.values, CRPS_mean.values)
#
# plt.show()

print(CRPS_mean.sel(step=pd.Timedelta(lead_time,'D')))

# exit()
# print(CRPS_loc.mean(dim='time'))

MAE = xs.mae(era_anomaly, hindcast_anomaly, dim=['lat', 'lon'], skipna=True)
MSE = xs.mse(era_anomaly, hindcast_anomaly, dim=['lat', 'lon'], skipna=True)
RMSE = xs.rmse(era_anomaly, hindcast_anomaly, dim=['lat', 'lon'], skipna=True)

threshold_brier_score = xs.threshold_brier_score(era, hindcast, [3,5,10,15,20,25], dim=None)
threshold_brier_score_era = xs.threshold_brier_score(era, era_mean, [3,5,10,15,20,25], dim=None)
print('Threshold brier score:')
print(threshold_brier_score)
print('MAE')
print(MAE)
print('MSE')
print(MSE)
print('RMSE')
print(RMSE)
# brier_score = xs.brier_score(obs3>.5, (fct3>.5).mean('member'))

# plt.figure()
# plt.plot(hindcast_anomaly.step.values, threshold_brier_score)
# plt.show()


CRPS2 = xs.crps_ensemble(era, hindcast,dim='time')
# print(CRPS2.sel(lat=latitude, lon=longitude,step=pd.Timedelta(lead_time,'D')))
