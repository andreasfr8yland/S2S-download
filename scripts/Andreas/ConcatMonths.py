import xarray as xr

months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
hindcast_anomalies_month={}
for x in months:
    hindcast_anomalies_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_anomalies_{}.nc'.format(x))
    hindcast_anomalies_month[x] = hindcast_anomalies_month[x].where(hindcast_anomalies_month[x].time.dt.month==x, drop=True)

hindcast_anomalies = xr.concat([hindcast_anomalies_month[1], hindcast_anomalies_month[2], hindcast_anomalies_month[3], hindcast_anomalies_month[4], hindcast_anomalies_month[5], hindcast_anomalies_month[6], hindcast_anomalies_month[7], hindcast_anomalies_month[8], hindcast_anomalies_month[9], hindcast_anomalies_month[10], hindcast_anomalies_month[11], hindcast_anomalies_month[12]], dim='time')

hindcast_anomalies = hindcast_anomalies.assign_coords(lon=(((hindcast_anomalies.lon + 180) % 360) - 180))
hindcast_anomalies = hindcast_anomalies.sortby(hindcast_anomalies.lon)
hindcast_anomalies = hindcast_anomalies.sortby(hindcast_anomalies.time)


hindcast_anomalies.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_anomalies_concat.nc', mode="w")

hindcast_absolute_month={}
for x in months:
    hindcast_absolute_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_absolute_{}.nc'.format(x))
    hindcast_absolute_month[x] = hindcast_absolute_month[x].where(hindcast_absolute_month[x].time.dt.month==x, drop=True)

hindcast_absolute = xr.concat([hindcast_absolute_month[1], hindcast_absolute_month[2], hindcast_absolute_month[3], hindcast_absolute_month[4], hindcast_absolute_month[5], hindcast_absolute_month[6], hindcast_absolute_month[7], hindcast_absolute_month[8], hindcast_absolute_month[9], hindcast_absolute_month[10], hindcast_absolute_month[11], hindcast_absolute_month[12]], dim='time')

hindcast_absolute = hindcast_absolute.assign_coords(lon=(((hindcast_absolute.lon + 180) % 360) - 180))
hindcast_absolute = hindcast_absolute.sortby(hindcast_absolute.lon)
hindcast_absolute = hindcast_absolute.sortby(hindcast_absolute.time)


hindcast_absolute.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_absolute_concat.nc', mode="w")

era_anomalies_month={}
for x in months:
    era_anomalies_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_anomalies_{}.nc'.format(x))
    era_anomalies_month[x] = era_anomalies_month[x].where(era_anomalies_month[x].time.dt.month==x, drop=True)


era_anomalies = xr.concat([era_anomalies_month[1], era_anomalies_month[2], era_anomalies_month[3], era_anomalies_month[4], era_anomalies_month[5], era_anomalies_month[6], era_anomalies_month[7], era_anomalies_month[8], era_anomalies_month[9], era_anomalies_month[10], era_anomalies_month[11], era_anomalies_month[12]], dim='time')


era_anomalies = era_anomalies.assign_coords(lon=(((era_anomalies.lon + 180) % 360) - 180))
era_anomalies = era_anomalies.sortby(era_anomalies.lon)
era_anomalies = era_anomalies.sortby(era_anomalies.time)

era_anomalies.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_anomalies_concat.nc', mode="w")

era_std_month={}
for x in months:
    era_std_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_std_{}.nc'.format(x))
    era_std_month[x] = era_std_month[x].where(era_std_month[x].time.dt.month==x, drop=True)


era_std = xr.concat([era_std_month[1], era_std_month[2], era_std_month[3], era_std_month[4], era_std_month[5], era_std_month[6], era_std_month[7], era_std_month[8], era_std_month[9], era_std_month[10], era_std_month[11], era_std_month[12]], dim='time')


era_std = era_std.assign_coords(lon=(((era_std.lon + 180) % 360) - 180))
era_std = era_std.sortby(era_std.lon)
era_std = era_std.sortby(era_std.time)

era_std.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_std_concat.nc', mode="w")

era_absolute_month={}
for x in months:
    era_absolute_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_absolute_{}.nc'.format(x))
    era_absolute_month[x] = era_absolute_month[x].where(era_absolute_month[x].time.dt.month==x, drop=True)


era_absolute = xr.concat([era_absolute_month[1], era_absolute_month[2], era_absolute_month[3], era_absolute_month[4], era_absolute_month[5], era_absolute_month[6], era_absolute_month[7], era_absolute_month[8], era_absolute_month[9], era_absolute_month[10], era_absolute_month[11], era_absolute_month[12]], dim='time')


era_absolute = era_absolute.assign_coords(lon=(((era_absolute.lon + 180) % 360) - 180))
era_absolute = era_absolute.sortby(era_absolute.lon)
era_absolute = era_absolute.sortby(era_absolute.time)

era_absolute.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_absolute_concat.nc', mode="w")

era_mean_month={}
for x in months:
    era_mean_month[x] = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_mean_{}.nc'.format(x))
    era_mean_month[x] = era_mean_month[x].where(era_mean_month[x].time.dt.month==x, drop=True)


era_mean = xr.concat([era_mean_month[1], era_mean_month[2], era_mean_month[3], era_mean_month[4], era_mean_month[5], era_mean_month[6], era_mean_month[7], era_mean_month[8], era_mean_month[9], era_mean_month[10], era_mean_month[11], era_mean_month[12]], dim='time')


era_mean = era_mean.assign_coords(lon=(((era_mean.lon + 180) % 360) - 180))
era_mean = era_mean.sortby(era_mean.lon)
era_mean = era_mean.sortby(era_mean.time)

era_mean.to_netcdf('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_mean_concat.nc', mode="w")
