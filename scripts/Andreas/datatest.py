#import properscoring as ps
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import xskillscore as xs

from scripts.Henrik.data_handler import ERA5, ECMWF_S2SH
#import scripts.Henrik.xarray_helpers as xh
from S2S.local_configuration import config
#import scripts.Henrik.models as models

#domainID = 'NVK'
#var      = 'u10'

#t_start  = (2000,6,25)
#t_end    = (2019,8,10)

#clim_t_start  = (2000,1,1)
#clim_t_end    = (2021,1,4)

#observations	= ERA5().load(var,clim_t_start,clim_t_end,domainID)[var]
#hindcast	= ECMWF_S2SH().load(var,t_start,t_end,domainID)[var]

open_data = xr.open_dataset('/nird/projects/NS9853K/DATA/S2S/hindcast/ECMWF/sfc/u10/u10_CY46R1_2020-06-25_cf.grb',engine='cfgrib')
