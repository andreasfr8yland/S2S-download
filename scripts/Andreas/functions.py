import xarray as xr
import numpy as np

# Functions for calculating ACC on the temporal scale, for each grid point:

# ACC with subtracted error - input must be numpy arrays
def ACC_anom_error(FC_anom,OBS_anom):
	top=np.mean((FC_anom-np.mean(FC_anom, axis=0))*(OBS_anom-np.mean(OBS_anom, axis=0)), axis=0)
	bottom=np.sqrt(np.mean(np.square(FC_anom-np.mean(FC_anom, axis=0)), axis=0)*np.mean(np.square(OBS_anom-np.mean(OBS_anom, axis=0)), axis=0))
	ACC_anom=np.divide(top,bottom)
	return ACC_anom

# ACC without subtracted error - input must be numpy arrays
def ACC_anom(FC_anom,OBS_anom):
	top=np.mean((FC_anom*OBS_anom), axis=0)
	bottom=np.sqrt(np.mean(np.square(FC_anom), axis=0)*np.mean(np.square(OBS_anom), axis=0))
	ACC_anom=np.divide(top,bottom)
	return ACC_anom

# ACC without subtracted error - input must be numpy arrays - using nanmean
def ACC_anom_nan(FC_anom,OBS_anom):
	top=np.nanmean((FC_anom*OBS_anom), axis=0)
	bottom=np.sqrt(np.nanmean(np.square(FC_anom), axis=0)*np.nanmean(np.square(OBS_anom), axis=0))
	ACC_anom=np.divide(top,bottom)
	return ACC_anom

# # Not completely sure about this one
# def ACC_anom_C_version(FC_anom,FC_anom_WM,OBS_anom, OBS_anom_WM):
# 	top=np.sum(np.subtract(FC_anom,FC_anom_WM)*np.subtract(OBS_anom,OBS_anom_WM))
# 	bottom=np.sqrt(np.sum(np.square(np.subtract(FC_anom,FC_anom_WM)))*np.sum(np.square(OBS_anom-OBS_anom_WM)))
# 	ACC_anom=np.divide(top,bottom)
# 	return ACC_anom_C_version

# # ACC plotting routine
# def ACC_map(lon,lat,acc,filename,central_longitude=10)
#     fig, ax = plt.subplots(figsize=(8,8))
#     ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=10))
#     lvs = np.arange(-0.1, 0.8, 0.1)
#     lvs_lines = np.arange(-1, 1, 0.1)
#     CS = plt.contourf(lon, lat, np.transpose(acc), levels=lvs, extend='both', transform=ccrs.PlateCarree())
#     CS2 = plt.contour(lon, lat, np.transpose(acc), levels=lvs_lines, extend='both', transform=ccrs.PlateCarree(), colors='k', linewidths=0.5, linestyles='--')
#     plt.colorbar(CS)
#     plt.clabel(CS2,fmt='%1.2f')
#     ax.coastlines(resolution='50m')
#     ax.gridlines()
#     ax.add_feature(cartopy.feature.BORDERS)
#     if SEASONAL==True:
#         ax.set_title('Anomaly Correlation - Leadtime {:01d} days for {}'.format(lead_time, Season), fontsize=14)
#     elif SEASONAL==False:
#         ax.set_title('Anomaly Correlation - Leadtime {:01d}'.format(lead_time), fontsize=14)
#     ax.set_xlabel('Longitude')
#     ax.set_ylabel('Latitude')
#     ax.clabel(CS2, inline=1, fontsize=10, fmt='%1.2f')
#     if SEASONAL==True:
#         fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/ACC_map_leadtime_{:01d}_{}.png'.format(lead_time, Season), transparent=True)
#     elif SEASONAL==False:
#         fig.savefig('/lustre/storeB/project/fou/om/Sesongvind/Figurer/S2S/ACC_map/ACC_map_leadtime_{:01d}.png'.format(lead_time), transparent=True)
#     plt.show()
#     return
