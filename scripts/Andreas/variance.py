import xarray as xr

hindcast_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_anomalies_concat.nc')
hindcast = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/hindcast_absolute_concat.nc')
era_anomalies = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_anomalies_concat.nc')
era = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_absolute_concat.nc')
era_std = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_std_concat.nc')
era_climatology = xr.open_dataset('/lustre/storeB/project/fou/om/Sesongvind/Data/S2S/NorthAtlantic/era_mean_concat.nc')
